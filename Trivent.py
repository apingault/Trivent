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
# from lxml import etree
import yaml

# from ganga import *

import MySQLdb as mdb
import dbUtils as dbu # custom interfaces to TestBeam databases
from marlin import Marlin
# import marlinUtils as mu
# Import default config file
# Not needed here just for dumb editor not to complain about config not existing
import config_trivent as conf


'''
    Create xml geometry file from testBeam database
'''
# -----------------------------------------------------------------------------
def generateXMLGeometryFile(testBeamPeriod):
    db = mdb.connect(host='localhost', user='acqilc', passwd='RPC_2008', db='GEOMETRY')
    cur = db.cursor()

    print ("[{0}] - Selected TestBeam: '{1}'".format(os.path.basename(__file__), testBeamPeriod) )

    testBeamIdx = dbu.selectTestBeam(cur, testBeamPeriod)
    print ("[{0}] - TestBeam index : '{1}'".format(os.path.basename(__file__), testBeamIdx) )

    difList = dbu.selectDifList(cur, testBeamIdx).fetchall()
    layerList = dbu.selectLayerList(cur, testBeamIdx).fetchall()

    xmlName = conf.geomFile
    dbu.createGeomXml(xmlName, difList, layerList)


# # -----------------------------------------------------------------------------
# '''
# '''
# def xmlValueBuilder(rootHandle, xmlHandleName, value, parameterType=None, option=None, optionValue=None, xmlParList=None):
#     xmlHandle = etree.SubElement(rootHandle, "parameter", name=xmlHandleName)
#     if parameterType is not None:
#         xmlHandle.set("type", parameterType)
#     if option is not None:
#         xmlHandle.set(option, optionValue)
#     xmlHandle.text = value
#     xmlParList[xmlHandleName] = [value]


# # -----------------------------------------------------------------------------
# '''
# '''
# def generateXML(inputFiles, outputFile, parList):
#     # generate XML
#     # List of parameter:value to easily print them afterwards
#     ## TAG: Execute
#     marlin = etree.Element('marlin')
#     execute = etree.SubElement(marlin, "execute")
#     for proc in conf.processorList:
#         processor = etree.SubElement(execute, "processor", name=proc[0])

#     ## TAG: Global
#     glob = etree.SubElement(marlin, "global")

#     # -- Processor Parameters
#     xmlValueBuilder(glob, "LCIOInputFiles", inputFiles, xmlParList=parList)
#     # --- Max number of evts to process
#     xmlValueBuilder(glob, "MaxRecordNumber", str(conf.maxEvt), xmlParList=parList)
#     # --- Number of evts to skip
#     xmlValueBuilder(glob, "SkipNEvents", str(conf.nSkipEvt), xmlParList=parList)
#     # --- Verbosity
#     xmlValueBuilder(glob, "Verbosity", conf.verbosity, option="options", optionValue="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT", xmlParList=parList)

#     ## TAG: Processors
#     for proc in conf.processorList:
#         processor = etree.SubElement(marlin, "processor", name=proc[0])
#         processor.set("type", proc[1])

#         # --- Verbosity
#         xmlValueBuilder(processor, "Verbosity", conf.verbosity, parameterType="string", option="options", optionValue="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT", xmlParList=parList)

#         xmlValueBuilder(processor, "InputCollectionNames", conf.inputCollectionName,       "StringVec", xmlParList=parList)
#         xmlValueBuilder(processor, "OutputCollectionName", conf.outputCollectionName,      "string",    xmlParList=parList)
#         xmlValueBuilder(processor, "LCIOOutputFile",       outputFile + ".slcio",            "string",    xmlParList=parList)
#         xmlValueBuilder(processor, "ROOTOutputFile",       outputFile + ".root",             "string",    xmlParList=parList)
#         xmlValueBuilder(processor, "PlotFolder",           conf.plotPath,                  "string",    xmlParList=parList)
#         xmlValueBuilder(processor, "SetupGeometry",        conf.geomFile,                  "string",    xmlParList=parList)
#         xmlValueBuilder(processor, "GainCorrectionMode",   str(conf.gainCorrection),       "bool",      xmlParList=parList)
#         xmlValueBuilder(processor, "ElectronicNoiseCut",   str(conf.electronicNoiseCut),   "int",       xmlParList=parList)
#         xmlValueBuilder(processor, "LayerGap",             str(conf.layerGap),             "double",    xmlParList=parList)
#         xmlValueBuilder(processor, "LayerCut",             str(conf.layerCut),             "int",       xmlParList=parList)
#         xmlValueBuilder(processor, "NoiseCut",             str(conf.noiseCut),             "int",       xmlParList=parList)
#         xmlValueBuilder(processor, "TimeWin",              str(conf.timeWin),              "int",       xmlParList=parList)
#         xmlValueBuilder(processor, "HasCerenkovDif",       str(conf.hasCerenkovDif),       "bool",      xmlParList=parList)
#         xmlValueBuilder(processor, "CerenkovDifId",        str(conf.cerenkovDifId),        "int",       xmlParList=parList)
#         xmlValueBuilder(processor, "CerenkovTimeWindow",   str(conf.cerenkovTimeWin),      "int",       xmlParList=parList)

#     # pretty string
#     s = etree.tostring(marlin, pretty_print=True)
#     with open(conf.xmlFile, "w") as outFile:
#         outFile.write(s)

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
        print ("\n[{0}] - File Not found, available files: ".format(os.path.basename(__file__)))
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


'''
'''
# -----------------------------------------------------------------------------
def checkPeriod(runNumber, runPeriod, configFile):
    ''' Check runNumber is associated to correct runPeriod
        Abort execution if not
    '''
    def periodError(goodPeriod):
        return "[{0}] - RunNumber '{1}' is from TestBeam '{2}', you selected '{3}' in configFile '{4}'".format(os.path.basename(__file__), runNumber, goodPeriod, conf.runPeriod, configFile)

    if runNumber < '726177' and (runPeriod != 'SPS_08_2012' or runPeriod != 'SPS_11_2012'):
        sys.exit(periodError('SPS_08_2012 or SPS_11_2012'))

    if runNumber >= '726177' and runNumber <= '726414' and runPeriod != 'SPS_12_2014':
        sys.exit(periodError('SPS_12_2014'))

    if runNumber >= '727760' and runNumber <= '728456' and runPeriod != 'SPS_04_2015':
        sys.exit(periodError('SPS_04_2015'))

    if runNumber >= '728501' and runNumber <= '728682' and runPeriod != 'PS_06_2015':
        sys.exit(periodError('PS_06_2015'))

    if runNumber >= '730436' and runNumber <= '730926' and runPeriod != 'SPS_10_2015':
        sys.exit(periodError('SPS_10_2015'))

    if runNumber >= '730927' and runNumber <= '732909' and runPeriod != 'SPS_06_2016':
        sys.exit(periodError('SPS_06_2016'))

    if runNumber >= '734927' and runPeriod != 'SPS_10_2016':
        sys.exit(periodError('SPS_10_2016'))

# -----------------------------------------------------------------------------
def createJob(executable, args = [], name='', comment='', backend='Local', backendCE='', voms=''):
    ''' Create Ganga job. Default backend is Local
    '''
    j = Job()
    j.application = Executable(exe=File(executable), args=args)
    j.name = name
    j.comment = comment

    j.backend = backend
    if backend == 'CREAM':
        j.backend.CE = backendCE
        try:
            gridProxy.voms = voms
        except NameError: # ganga > 6.3 no longer has the gridProxy credentials system
            j.backend.credential_requirements = VomsProxy(vo=voms)
    return j



# -----------------------------------------------------------------------------
def setCliOptions(marlin, xmlSection):
    ''' properly set cliOptions from configfile for marlin.py
    '''
    sectionName = xmlSection.name + "."
    for param, value in vars(xmlSection).items():
        if param != 'name':
            marlin.setCliOption(sectionName + param, value)




# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
def main():
    '''
    '''
    scriptName = os.path.basename(__file__) # For logging clarity

    runNumber = 0
    runListArg = None
    configFile = "config_trivent"

    # --- Parse CLI arguments
    if len(sys.argv) > 1:
        # --- Load configuration File
        configFile = sys.argv[1]
        # try:
        exec("import {0} as conf".format(configFile))
        # except (ImportError, SyntaxError):
            # sys.exit("[{0}] - Cannot import config file '{1}'".format(scriptName, configFile))
        # --- /Load configuration File
        if len(sys.argv) > 2:
            # --- Load runList
            runListArg = sys.argv[2]
    else:
        sys.exit("[{0}] - Please give : configFile - runNumber(s)(optional if set up in configFile)".format(scriptName))
    # --- /Parse CLI arguments



    # --- Load List of runs
    if runListArg is None: # if no runs on CLI, load list from configFile
        try:
            runList = conf.runList
        except AttributeError:
            sys.exit("[{0}] - No runs specified at command line or in configFile...exiting".format(scriptName))
    else:
        runList = runListArg.split(',')
    # --- /Load List of runs



    for run in runList:
        checkPeriod(str(run), conf.runPeriod, configFile)
        runNumber = int(run)
        print ("[{0}] - Looking for files to process for run '{1}' in {2}... ".format(scriptName, runNumber, conf.inputPath), end="")
        #   TODO Remove Hardcoding here
        # stringToFind = 'DHCAL_%d_SO_Antoine' % runNumber
        stringToFind = 'DHCAL_%d_SO' % runNumber
        if os.path.exists(conf.inputPath) is False:
            sys.exit("\n[{0}] - Folder '{1}' does not exist...exiting".format(scriptName, conf.inputPath))

        found = findFile(conf.inputPath, stringToFind)
        if found is False:
            sys.exit('Exiting')
        print ('OK')

        conf.glob.LCIOInputFiles = conf.inputFile.format(conf.inputPath, runNumber)
        outputFile = conf.outputFile.format(conf.outputPath, runNumber)
        conf.triventProc.LCIOOutputFile = outputFile + ".slcio"
        conf.triventProc.ROOTOutputFile = outputFile + ".root"
        print ("[{0}] - output file : {1}.slcio".format(scriptName, outputFile))
        print ("[{0}] - MARLIN_DLL: {1}".format(scriptName, conf.marlinLib))

        # Checking if trivent alredy run
        if os.path.exists("{0}.slcio".format(outputFile)) is True:
            sys.exit("[{0}] - OutputFile already present...exiting".format(scriptName))


        # Looking for or generating xml Geometry File
        if os.path.exists(conf.geomFile) is False:
            print ("[{0}] - No geometry file found, creating one from database for period '{1}'...".format(scriptName, conf.runPeriod))
            generateXMLGeometryFile(conf.runPeriod)
            print ("[{0}] - No geometry file found, creating one from database for period '{1}'...OK -> '{2}' ".format(scriptName, conf.runPeriod, conf.geomFile))
        else:
            print ("[{0}] - Found geometry file '{1}'".format(scriptName, conf.geomFile))



        # xmlParameterList = {}
        # generateXML(inputFiles, outputFile, xmlParameterList)


        # Printing modification to xml file
        print("\n[{0}] ========================".format(scriptName))
        print("[{0}] --- Dumping modified xml parameters: ".format(scriptName))
        for par, value in vars(conf.glob).items():
            if par != 'name':
                print ("[{0}] \t\t{1}:\t{2}".format(scriptName, par, value))
        for par, value in vars(conf.triventProc).items():
            if par != 'name':
                print ("[{0}] \t\t{1}:\t{2}".format(scriptName, par, value))
        print("[{0}] ========================\n".format(scriptName))


        # Marlin configuration
        marlinCfgFile = conf.marlinCfgFile.format(runNumber)
        marlin = Marlin()
        marlin.setXMLConfig(conf.xmlFile)
        marlin.setLibraries(conf.marlinLib)
        marlin.setILCSoftScript(conf.initILCSoftScript)
        setCliOptions(marlin, conf.glob) # TODO: Move to Marlin.py
        setCliOptions(marlin, conf.triventProc)
        marlin.writeConfigFile(marlinCfgFile)




        # Running locally
        if conf.runOnGrid is False:
            log = open(conf.logFile.format(conf.logPath, runNumber), "w", 1) # line-buffered

            print("\n[{0}] ========================".format(scriptName))
            print ('[{0}] --- Running Marlin...'.format(scriptName))
            print ("[{0}] --- Ouput is logged to '{1}'".format(scriptName, log))
            beginTime = time.time()
            # subprocess.call(["Marlin", conf.xmlFile], env=dict(os.environ, MARLIN_DLL=conf.marlinLib, MARLINDEBUG="1"), stdout=log)
            subprocess.call(['python', 'run_marlin.py', marlinCfgFile], stdout=log, stderr=log)
            # subprocess.call(['python', 'run_marlin.py', marlinCfgFile])
            print ('[{0}] - Running Marlin...OK - '.format(scriptName), end='')
            try:
                elapsedTime(beginTime)
            except ValueError:
                print ("Can't print time...")
            print("[{0}] ========================\n".format(scriptName))

            # print ('[{0}] - Removing xmlFile...', end='')
            # subprocess.Popen(["rm", conf.xmlFile])
            # print ("OK")

        else:
            with open(marlinCfgFile, 'r') as ymlfile:
                cfg = yaml.load(ymlfile)
            gridSection = {}
            cfg['Grid'] = gridSection
            gridSection['downloader'] = conf.gridDownloader
            gridSection['uploader'] = conf.gridUploader
            gridSection['LCG_CATALOG_TYPE'] = conf.LCG_CATALOG_TYPE
            gridSection['LFC_HOST'] = conf.LFC_HOST
            gridSection['inputFiles'] = conf.gridInputFiles
            with open(marlinCfgFile, 'w') as ymlfile:
                ymlfile.write(yaml.dump(cfg, default_flow_style=False))

            try:
                print ("[{0}] --- Submiting Job ... ".format(scriptName))
                # Navigate jobtree, if folder doesn't exist create it
                treePath = conf.runPeriod + '/' + scriptName[:-3] # Remove .py at end of scriptName
                # if jobtree.exists(treePath) is False: # Always returns true...
                jobtree.cd('/') # make sure we are in root folder
                try :
                    jobtree.cd(treePath)
                except TreeError:
                    try:
                        jobtree.mkdir(treePath)
                    except: # mkdir should write all missing folder if any....apparently not true
                        print("WhatThe?")
                        jobtree.mkdir(conf.runPeriod)
                        jobtree.mkdir(treePath)
                    jobtree.cd(treePath)

                eos_installation ='/afs/cern.ch/project/eos/installation/user/'
                eos_home='/eos/user/a/apingaul/CALICE/'

                # Update ganga configuration for eos access
                config.Output.MassStorageFile['defaultProtocol'] = 'root://eosuser.cern.ch'
                config.Output.MassStorageFile['uploadOptions']['cp_cmd'] = eos_installation + 'bin/eos.select cp'
                config.Output.MassStorageFile['uploadOptions']['ls_cmd'] = eos_installation + 'bin/eos.select ls'
                config.Output.MassStorageFile['uploadOptions']['mkdir_cmd'] = eos_installation + 'bin/eos.select mkdir'
                config.Output.MassStorageFile['uploadOptions']['path'] = eos_home
                # Print it
                print (config.Output.MassStorageFile['defaultProtocol'])
                print (config.Output.MassStorageFile['uploadOptions']['cp_cmd'])
                print (config.Output.MassStorageFile['uploadOptions']['ls_cmd'])
                print (config.Output.MassStorageFile['uploadOptions']['mkdir_cmd'])
                print (config.Output.MassStorageFile['uploadOptions']['path'])




                inputFiles = []
                for f in conf.gridInputFiles:
                    inputFiles.append(LocalFile(f))
                inputFiles.append(LocalFile(marlinCfgFile))
                print ('inputFiles:\n', inputFiles)

                inputData=[]
                inputDataFileList = []
                inputDataFileList.append(conf.glob.LCIOInputFiles)
                for item in [inputDataFileList]:
                    l = []
                    print ("\nitem=",item,"\n")
                    for f in item:
                        l.append(MassStorageFile(f))
                    print ("\nl=",l,"\n")
                    for f in l:
                        print ("\nf=",f,"\n")

                inputData = GangaDataset(treat_as_inputfiles=False, files=[f for f in l])

                print ('\n\ninputDataType:\n', type(inputData))
                print ('\n\ninputData:\n', inputData)
                for item in inputData:
                    print (type(item))
                    print (item)

                j = createJob(executable='run_marlin.py', args=[marlinCfgFile], name=str(runNumber), backend=conf.backend, backendCE=conf.CE, voms=conf.voms)
                j.comment = "Trivent " + conf.runPeriod
                j.outputfiles = [MassStorageFile(namePattern="*.*", outputfilenameformat='GridOutput/Trivent/{jid}/{fname}')]
                j.inputfiles = inputFiles
                j.inputdata = inputData

                # Save job as txtfile locally
                # export(jobs(j.id), 'my_job.txt')

                jobtree.add(j)
                j.submit()
                print ("\n[{0}] ... submitting job done.\n".format(scriptName))

                #queues.add(j.submit)
            except:
                print ("[{0}] --- Failed to submit job ".format(scriptName))
                raise

if  __name__ == '__main__':
    main()
