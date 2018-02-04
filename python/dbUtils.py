#!/usr/bin/env python

'''
    Utils to to get informations from GEOMETRY db
    After installing mysql server
    sudo mysqld_safe
    mysql -u root
    mysql> CREATE DATABASE GEOMETRY;
    mysql> USE GEOMETRY;
    mysql> CREATE USER 'acqilc'@'localhost'
    mysql> GRANT ALL PRIVILEGES ON GEOMETRY.* TO 'acqilc'@'localhost';
    # ON lyosdhcal10
    mysqldump -u acqilc -pRPC_2008  GEOMETRY > geometry_m3_2017.sql
    #locally
    mysqldump -u acqilc -pRPC_2008  GEOMETRY < geometry_m3_2017.sql
    if not working:
    mysql -u acqilc -pRPC_2008
    use GEOMETRY;
    source geometry_m3_2017;
'''


from __future__ import print_function  # import print function from py3 if running py2.x
import sys
import pymysql as pmsql


# -----------------------------------------------------------------------------
def selectListOfTestBeams(dbcursor):
    '''
    '''
    dbcursor.execute("SELECT * FROM VERSIONS")
    return dbcursor.fetchall()


# -----------------------------------------------------------------------------
def dumpDifList(dbcursor, testBeamId):
    ''' Dump list of dif in TestBeam db
    '''
    rows = selectDifList(dbcursor, testBeamId).fetchall()
    print ("\n\n--------------- Dumping Dif List ---------------")
    print ("------------------------------------------------")
    print (" - - - - - - - - - - - - - - - - - - - - - - - -")
    print ("|   Dif   \t| Layer\t|   X\t|   Y   \t|")
    print ("| - - - - - - - - - - - - - - - - - - - - - - - |")
    for row in rows:
        print ("|  ", row[0], "   \t| ", row[1], " \t| ", row[2], "\t| ", row[3], " \t|")
    print ("| - - - - - - - - - - - - - - - - - - - - - - - |")


# -----------------------------------------------------------------------------
def dumpLayerList(dbcursor, testBeamId):
    ''' Dump list of chamber in TestBeam db
    '''
    rows = selectLayerList(dbcursor, testBeamId).fetchall()
    print ("--- Dumping Chambers ---")
    print (" - - - - - - - - - - - - - - - -")
    print ("| Layer\t|   X\t|   Y\t|   Z\t|")
    print ("| - - - - - - - - - - - - - - - |")
    for row in rows:
        print ("|  ", row[0], "\t| ", row[1], "\t| ", row[2], "\t| ", row[3], "\t|")
    print (" - - - - - - - - - - - - - - - -")


# -----------------------------------------------------------------------------
def dumpLayerPositions(dbcursor, testBeamId):
    '''
    '''
    rows = selectLayerPositionList(dbcursor, testBeamId).fetchall()
    print ("--- Dumping Chambers Position ---")
    print (" - - - - - - - - - - - - - - - -")
    for row in rows:
        print ("|  ", row[0], "\t| ", row[3], "\t|")
    print (" - - - - - - - - - - - - - - - -")


# -----------------------------------------------------------------------------
def selectTestBeam(dbcursor, testBeamName):
    ''' Get Idx for selected TestBeam
    '''
    rows = selectListOfTestBeams(dbcursor)
    found = False
    for row in rows:
        if testBeamName in row:
            found = True
            dbcursor.execute("SELECT IDX FROM VERSIONS WHERE TESTNAME = %s", testBeamName)
            return dbcursor.fetchone()[0]  # testBeamIdx
    if found is False:
        print ("[dbUtils.py] - Selected Test Beam '%s' not found in database..." % testBeamName)
        rows = selectListOfTestBeams(dbcursor)
        print ("[dbUtils.py] - List of available Test Beam :")
        for line in rows:
            print ("\t\t- '%s' " % str(line[1]))
        sys.exit(0)


# -----------------------------------------------------------------------------
def selectLayerList(dbcursor, testBeamId):
    ''' Get list of layers ( Not very useful now...)
    '''
    dbcursor.execute("SELECT NUM,X0,Y0,Z0 FROM PLANS WHERE VERSIONID = (%s)", testBeamId)
    return dbcursor


# -----------------------------------------------------------------------------
def selectLayerPositionList(dbcursor, testBeamId):
    '''Get list of layers with positions
    '''
    dbcursor.execute("select NUM,X0,Y0,Z0,X1,Y1,Z1,TYPE,(SELECT NUM FROM PLANS WHERE PLANS.IDX=CHAMBERS.PLANID) FROM CHAMBERS WHERE (SELECT NUM FROM PLANS WHERE PLANS.IDX=CHAMBERS.PLANID) IS NOT NULL AND VERSIONID = (%s)", testBeamId)
    return dbcursor


# -----------------------------------------------------------------------------
def selectDifList(dbcursor, testBeamId):
    '''Get List of difs/Chambers/difposition associated with testBeam
    '''
    dbcursor.execute("select NUM,(SELECT NUM FROM CHAMBERS WHERE CHAMBERS.IDX=DIFS.CHAMBERID),DI,DJ,POLI,POLJ FROM DIFS WHERE (SELECT NUM FROM CHAMBERS WHERE CHAMBERS.IDX=DIFS.CHAMBERID) IS NOT NULL AND VERSIONID = (%s)", testBeamId)
    return dbcursor


# -----------------------------------------------------------------------------
def selectEnergyFromRun(dbcursor, runNumber):
    '''Get energy associated with runNumber
    '''
    dbcursor.execute("SELECT ENERGY FROM LOGBOOK WHERE RUN = (%s)" % runNumber)
    return dbcursor


# -----------------------------------------------------------------------------
def selectFilePathFromRun(dbcursor, runNumber, isCompressed=False):
    '''Get filePath associated with runNumber
    '''
    dbcursor.execute("select LOCATION from FILES WHERE RUN=(SELECT RUN FROM RUNS WHERE RUN = (%s))  AND COMPRESS = (%s)" % (runNumber, isCompressed))
    # print (dbcursor.fetchall())
    return dbcursor


# -----------------------------------------------------------------------------
def selectFileNameFromRun(dbcursor, runNumber, isCompressed=False):
    '''Get fileName associated with runNumber
    '''
    filePath = selectFilePathFromRun(dbcursor, runNumber, isCompressed)
    fileList = []
    for row in filePath.fetchall()[0]:
        fileName = str(row).split('/')
        fileList.append(fileName[-1])
    print (fileList)
    return dbcursor


# -----------------------------------------------------------------------------
def xmlValueBuilder(rootHandle, xmlHandleName, value, parameterType=None, option=None, optionValue=None, xmlParList=None):
    '''
    '''
    xmlHandle = etree.SubElement(rootHandle, "parameter", name=xmlHandleName)
    if parameterType is not None:
        xmlHandle.set("type", parameterType)
    if option is not None:
        xmlHandle.set(option, optionValue)
    xmlHandle.text = value
    if xmlParList is not None:
        xmlParList[xmlHandleName] = [value]


# -----------------------------------------------------------------------------
def buildList(dbList, rowElements):
    '''
    '''
    rowList = ['\n']
    for difRow in dbList:
        row = []
        for i in range(0, rowElements):
            row.append(str(difRow[i]))
        # print ("row:", row)
        rowList.append((',').join(row))

    rowList.append('\n')
    return '\n'.join(rowList)


'''
    difList = selectDifList(dbcursror, testBeamIdx).fetchall()
    chamberList = selectLayerList(dbcursror, testBeamIdx).fetchall()
'''


# -----------------------------------------------------------------------------
def createGeomXml(xmlFileName, difList, layerList):
    '''
    '''
    rootHandle = etree.Element('setup_geom')

    difToPrint = buildList(difList, 6)
    xmlValueBuilder(rootHandle, "DifGeom", difToPrint)

    # layerToPrint = buildList(layerList, 4)
    # xmlValueBuilder(rootHandle, "ChamberGeom", layerToPrint)

    # pretty string
    s = etree.tostring(rootHandle, pretty_print=True)
    print (s)
    with open(xmlFileName, "w") as outFile:
        outFile.write(s)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
def main():
    '''
    '''
    try:
        db = pmsql.connect(host='localhost', user='acqilc', passwd='RPC_2008', db='GEOMETRY')
        testName = "SPS_06_2015"
        cur = db.cursor()

        print ("[dbUtils.py] - Selected TestBeam: '%s'" % testName)

        testBeamIdx = selectTestBeam(cur, testName)
        print ("[dbUtils.py] - TestBeam index :", testBeamIdx)

        # dumpDifList(cur, testBeamIdx)

        difList = selectDifList(cur, testBeamIdx).fetchall()
        layerList = selectLayerList(cur, testBeamIdx).fetchall()

        xmlName = "geometry_%s.xml" % testName
        createGeomXml(xmlName, difList, layerList)

        # selectPathFileFromRun(cur, '726254')
        selectEnergyFromRun(cur, '726254')
        print ("[dbUtils.py] - Energy for run '%s' = %s GeV" % ('726254', cur.fetchone()[0]))
        selectFilePathFromRun(cur, '726254')
        print (cur.fetchall())
        selectFileNameFromRun(cur, '726254')
        print (cur.fetchall())

        # dumpListOfTestBeams(cur)
        # dumpDifList(cur, testBeamIdx)

        # Get files location for given Energy and TestBeamId,
        # cur.execute("select LOCATION from FILES WHERE RUN  in (SELECT RUN FROM RUNS WHERE ENERGY = (%s) AND VERSIONID = (%s) ) AND COMPRESS= %s", (energy, versionId, isCompressed)

        # Retrieve all rows at once
        # print row[0], row[1]

        # print ('[Streamout.py] - Running Marlin...OK, - ', end='')

        # Retrieve data one row at a time
        # for i in range(cur.rowcount):
        #     row = cur.fetchone()
        #     print row

        # create Dict Cursor -> Return data as python dict instead of python list
        # -> can use print row["Id"], row["Name"]

        # cur =db.cursor(pmsql.cursors.DictCursor)
        # cur.execute("SELECT * FROM Writers LIMIT 4")

        # rows = cur.fetchall()

        # for row in rows:
        #     print row["Id"], row["Name"]

        # Close all cursors
        cur.close()
        # Close all databases
        db.close()

    except pmsql.Error as e:
        print ("*** *** Error %d: %s" % (e.args[0], e.args[1]))
        sys.exit(1)


if __name__ == '__main__':
    main()
