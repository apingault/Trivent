"""
    Configuration module for Trivent Marlin processor
    Generate the xml file with data imported from external config file
"""

from __future__ import print_function # import print function from py3 if running py2.x

import os
import sys
import time
import subprocess
import shlex # split properly command line for Popen
from lxml import etree

import MySQLdb as mdb
import dbUtils as dbu # custom interfaces to TestBeam databases

# Import default config file
# Not needed here just for dumb editor not to complain about config not existing
import config_trivent as config


'''
    Create xml geometry file from testBeam database
'''
# -----------------------------------------------------------------------------
def generateXMLGeometryFile(testBeamPeriod):
    db = mdb.connect(host='localhost', user='acqilc', passwd='RPC_2008', db='GEOMETRY')
    cur = db.cursor()

    print ("[Trivent.py] - Selected TestBeam: '%s'" % testBeamPeriod)

    testBeamIdx = dbu.selectTestBeam(cur, testBeamPeriod)
    print ("[Trivent.py] - TestBeam index :", testBeamIdx)

    difList = dbu.selectDifList(cur, testBeamIdx).fetchall()
    layerList = dbu.selectLayerList(cur, testBeamIdx).fetchall()

    xmlName = config.geomFile
    dbu.createGeomXml(xmlName, difList, layerList)

# -----------------------------------------------------------------------------
'''
'''
def xmlValueBuilder(rootHandle, xmlHandleName, value, parameterType=None, option=None, optionValue=None, xmlParList=None):
    xmlHandle = etree.SubElement(rootHandle, "parameter", name=xmlHandleName)
    if parameterType is not None:
        xmlHandle.set("type", parameterType)
    if option is not None:
        xmlHandle.set(option, optionValue)
    xmlHandle.text = value
    xmlParList[xmlHandleName] = [value]
    
    
# -----------------------------------------------------------------------------
'''
''' 
def generateXML(inputFiles, outputFile, parList):
    # generate XML
    # List of parameter:value to easily print them afterwards
    ## TAG: Execute
    marlin = etree.Element('marlin')
    execute = etree.SubElement(marlin, "execute")
    for proc in config.processorList:
        processor = etree.SubElement(execute, "processor", name=proc[0])

    ## TAG: Global
    glob = etree.SubElement(marlin, "global")

    # -- Processor Parameters
    xmlValueBuilder(glob, "LCIOInputFiles", inputFiles, xmlParList=parList)
    # --- Max number of evts to process
    xmlValueBuilder(glob, "MaxRecordNumber", str(config.maxEvt), xmlParList=parList)
    # --- Number of evts to skip
    xmlValueBuilder(glob, "SkipNEvents", str(config.nSkipEvt), xmlParList=parList)
    # --- Verbosity
    xmlValueBuilder(glob, "Verbosity", config.verbosity, option="options", optionValue="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT", xmlParList=parList)

    ## TAG: Processors
    for proc in config.processorList:
        processor = etree.SubElement(marlin, "processor", name=proc[0])
        processor.set("type", proc[1])

        # --- Verbosity
        xmlValueBuilder(processor, "Verbosity", config.verbosity, parameterType="string", option="options", optionValue="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT", xmlParList=parList)
        
        xmlValueBuilder(processor, "InputCollectionNames", config.inputCollectionName,       "StringVec", xmlParList=parList)
        xmlValueBuilder(processor, "OutputCollectionName", config.outputCollectionName,      "string",    xmlParList=parList)
        xmlValueBuilder(processor, "LCIOOutputFile",       outputFile + ".slcio",            "string",    xmlParList=parList)
        xmlValueBuilder(processor, "ROOTOutputFile",       outputFile + ".root",             "string",    xmlParList=parList)
        xmlValueBuilder(processor, "PlotFolder",           config.plotPath,                  "string",    xmlParList=parList)
        xmlValueBuilder(processor, "SetupGeometry",        config.geomFile,                  "string",    xmlParList=parList)
        xmlValueBuilder(processor, "GainCorrectionMode",   str(config.gainCorrection),       "bool",      xmlParList=parList)
        xmlValueBuilder(processor, "ElectronicNoiseCut",   str(config.electronicNoiseCut),   "int",       xmlParList=parList)
        xmlValueBuilder(processor, "LayerGap",             str(config.layerGap),             "double",    xmlParList=parList)
        xmlValueBuilder(processor, "LayerCut",             str(config.layerCut),             "int",       xmlParList=parList)
        xmlValueBuilder(processor, "NoiseCut",             str(config.noiseCut),             "int",       xmlParList=parList)
        xmlValueBuilder(processor, "TimeWin",              str(config.timeWin),              "int",       xmlParList=parList)
        xmlValueBuilder(processor, "HasCerenkovDif",       str(config.hasCerenkovDif),       "bool",      xmlParList=parList)
        xmlValueBuilder(processor, "CerenkovDifId",        str(config.cerenkovDifId),        "int",       xmlParList=parList)
        xmlValueBuilder(processor, "CerenkovTimeWindow",   str(config.cerenkovTimeWin),      "int",       xmlParList=parList)

    # pretty string
    s = etree.tostring(marlin, pretty_print=True)
    with open(config.xmlFile, "w") as outFile:
        outFile.write(s)
    
'''
'''
# -----------------------------------------------------------------------------
def findFile(folder, findString):
    found = False
    command_line = "ls %s" % folder
    args = shlex.split(command_line)
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout_list = proc.communicate()[0].split('\n')
    for line in stdout_list:
        if findString in line:
            found = True
    if found is False:
        print ("\n[Trivent.py] - File Not found, available files: ")
        for line in stdout_list:
            print ("\t%s" % line)
    return found

'''
'''
# -----------------------------------------------------------------------------
def elapsedTime(startTime):
    t_sec = round(time.time() - startTime)
    (t_min, t_sec) = divmod(t_sec, 60)
    (t_hour, t_min) = divmod(t_min, 60)
    print('Time passed: {:.0f} hour {:.0f} min {:.0f} sec'.format(t_hour, t_min, t_sec))


def checkPeriod(runNumber, runPeriod, configFile):
    if runNumber < '726177' and (runPeriod != 'SPS_08_2012' or runPeriod != 'SPS_11_2012'):
        print ("[Trivent.py] - RunNumber '%s' is from TestBeam 'SPS_08_2012' or 'SPS_11_2012', you selected '%s' in configFile '%s'" % (runNumber, config.runPeriod, configFile))
        sys.exit(0)

    if runNumber >= '726177' and runNumber <= '726414' and runPeriod != 'SPS_12_2014':
        print ("[Trivent.py] - RunNumber '%s' is from TestBeam 'SPS_12_2014', you selected '%s' in configFile '%s'" % (runNumber, config.runPeriod, configFile))
        sys.exit(0)

    if runNumber >= '727760' and runNumber <= '728456' and runPeriod != 'SPS_04_2015':
        print ("[Trivent.py] - RunNumber '%s' is from TestBeam 'SPS_04_2015', you selected '%s' in configFile '%s'" % (runNumber, config.runPeriod, configFile))
        sys.exit(0)

    if runNumber >= '728501' and runNumber <= '728682' and runPeriod != 'PS_06_2015':
        print ("[Trivent.py] - RunNumber '%s' is from TestBeam 'PS_06_2015', you selected '%s' in configFile '%s'" % (runNumber, config.runPeriod, configFile))
        sys.exit(0)

    if runNumber >= '730436' and runNumber <= '730926' and runPeriod != 'SPS_10_2015':
        print ("[Trivent.py] - RunNumber '%s' is from TestBeam 'SPS_10_2015', you selected '%s' in configFile '%s'" % (runNumber, config.runPeriod, configFile))
        sys.exit(0)

    if runNumber >= '730927':
        print ("[Trivent.py] - RunNumber '%s' is from TestBeam 'UNKNOWN', you selected '%s' in configFile '%s'" % (runNumber, config.runPeriod, configFile))
        sys.exit(0)




'''
'''
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
def main():
    runNumber = 0
    configFile = "config_trivent"

    runListArg = None
    # --- Parse CLI arguments
    if len(sys.argv) > 1:
        # --- Load configuration File
        configFile = sys.argv[1]
        try:
            exec("import %s as config" % configFile)
        except (ImportError, SyntaxError):
            print ("[Trivent.py] - Cannot import config file '%s'" % configFile)
            return
        # --- /Load configuration File
        if len(sys.argv) > 2:
            # --- Load runList
            runListArg = sys.argv[2]
    else:
        print ("[Trivent.py] - Please give : configFile - runNumber(s)(optional if set up in configFile)")
        return
    # --- /Parse CLI arguments



    # --- Load List of runs
    if runListArg is None: # if no runs on CLI, load list from configFile
        try:
            runList = config.runList
        except AttributeError:
            print ("[Trivent.py] - No runs specified at command line or in configFile...exiting")
            return
    else:
        runList = runListArg.split(',')
    # --- /Load List of runs



    for run in runList:
        checkPeriod(run, config.runPeriod, configFile)
        runNumber = int(run)
        print ('[Trivent.py] - Looking for files to process for run \'%d\' in \'%s\'... ' % (runNumber, config.inputPath), end="")

        stringToFind = 'DHCAL_%d_SO' % runNumber
        found = findFile(config.inputPath, stringToFind)
        # if found is False:
            # return
        print ('OK')
                
        inputFiles = config.inputFile % (config.inputPath, runNumber)
        outputFile = config.outputFile % (config.outputPath, runNumber) # extension slcio/root added in the xml generator
        print ("[Trivent.py] - output file : %s.slcio" % outputFile)
        print ("[Trivent.py] - MARLIN_DLL: %s" % (config.marlinLib))



        # Looking for or generating xml Geometry File
        if os.path.exists(config.geomFile) is False:
            print ("[Trivent.py] - No geometry file found, creating one from database for period '%s'..." % config.runPeriod)
            generateXMLGeometryFile(config.runPeriod)
            print ("[Trivent.py] - No geometry file found, creating one from database for period '%s'...OK -> '%s' " % (config.runPeriod, config.geomFile))
        else:
            print ("[Trivent.py] - Found geometry file '%s'" % config.geomFile)



        xmlParameterList = {}
        generateXML(inputFiles, outputFile, xmlParameterList)

        print("\n[Trivent.py] ========================")
        print("[Trivent.py] --- Dumping xml parameters: ")
        for par in xmlParameterList:
            print ("[Trivent.py] \t\t%s:\t%s" % (par, xmlParameterList[par]))
        print("[Trivent.py] ========================\n")
        
        log = open(config.logFile % (config.logPath, runNumber), "w", 1) # line-buffered
        
        print("\n[Trivent.py] ========================")
        print ('[Trivent.py] --- Running Marlin...')
        print ('[Trivent.py] --- Ouput is logged to \'%s\'' % log)
        beginTime = time.time()
        # subprocess.call(["Marlin", config.xmlFile], env=dict(os.environ, MARLIN_DLL=config.marlinLib, MARLINDEBUG="1"), stdout=log)
        print ('[Trivent.py] - Running Marlin...OK, - ', end='')
        elapsedTime(beginTime)
        print("[Trivent.py] ========================\n")
        
        print ('[Trivent.py] - Removing xmlFile...', end='')
        # subprocess.Popen(["rm", config.xmlFile])
        print ("OK")

if  __name__ == '__main__':
    main()
