#!/usr/bin/env python
"""
    Configuration module for Marlin processor
    export PYTHONPATH=/cvmfs/ganga.cern.ch/Ganga/install/LATEST/python:${PYTHONPATH}
"""

from __future__ import print_function  # import print function from py3 if running py2.x
from __future__ import absolute_import

import os
import sys
import time
import subprocess

from marlin import Marlin
import pymysql as pmsql

import dbUtils as dbu  # custom interfaces to TestBeam databases

import yaml
try:
    import Ganga
except ImportError:
    print("Ganga not found on system, won't run on grid")
    gangaNotFound = True
else:
    gangaNotFound = False
    from Ganga.Core.exceptions import IncompleteJobSubmissionError
    from Ganga.GPI import *
    # Job, jobtree, TreeError, config, LocalFile, MassStorageFile, GangaDataset

sys.path.append(".")


# -----------------------------------------------------------------------------
def checkForConfigFile(configFile):
    '''check that all needed parameters are defined in the configFile,
        Exit if a parameter is not found
    '''
    try:
        check = configFile.runPeriod
        # check = configFile.runList
        check = configFile.inputPath
        check = configFile.outputPath
        check = configFile.outputFile
        check = configFile.processorPath
        check = configFile.marlinLib
        check = configFile.marlinCfgFile
        check = configFile.marlinCfgPath
        check = configFile.initILCSoftScript

        if 'Streamout' not in configFile.processorType:
            check = configFile.geomFile

        check = configFile.marlinProc
        check = configFile.glob
        check = configFile.glob.LCIOInputFiles

        # Parameters to scp files if not present
        check = configFile.serverName
        check = configFile.serverDataPath

        if configFile.runOnGrid is True:  # Don't check for grid param if running locally
            check = configFile.voms
            check = configFile.lcg_catalog_type
            check = configFile.lfc_host
            check = configFile.gridInputFiles
            check = configFile.lfc_host
            check = configFile.SE
            check = configFile.backend
            check = configFile.CE
        else:
            check = configFile.logPath
            check = configFile.logFile

    except AttributeError as e:
        par = e.message.split(' ')[-1]
        sys.exit('[{}] - Parameter {} not set in configFile {} '.format(os.path.basename(__file__), par, configFile))


# -----------------------------------------------------------------------------
def elapsedTime(startTime):
    ''' Print time to finish job
    '''
    t_sec = round(time.time() - startTime)
    (t_min, t_sec) = divmod(t_sec, 60)
    (t_hour, t_min) = divmod(t_min, 60)
    print('Time passed: {:.0f} hour {:.0f} min {:.0f} sec'.format(t_hour, t_min, t_sec))


# -----------------------------------------------------------------------------
def generateXMLGeometryFile(testBeamPeriod, fName):
    db = pmsql.connect(host='localhost', user='acqilc', passwd='RPC_2008', db='GEOMETRY')
    cur = db.cursor()

    print("[{}] - Selected TestBeam: '{}'".format(os.path.basename(__file__), testBeamPeriod))

    testBeamIdx = dbu.selectTestBeam(cur, testBeamPeriod)
    print("[{}] - TestBeam index : '{}'".format(os.path.basename(__file__), testBeamIdx))

    difList = dbu.selectDifList(cur, testBeamIdx).fetchall()
    layerList = dbu.selectLayerList(cur, testBeamIdx).fetchall()

    dbu.createGeomXml(fName, difList, layerList)


# -----------------------------------------------------------------------------
def checkPeriod(runNumber, runPeriod, configFile):
    ''' Check runNumber is associated to correct runPeriod
        Abort execution if not
    '''
    def periodError(goodPeriod):
        return "[{}] - RunNumber '{}' is from TestBeam '{}', you selected '{}' in configFile '{}'".format(
            os.path.basename(__file__), runNumber, goodPeriod, conf.runPeriod, configFile)

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

    if runNumber >= '733626' and runNumber <= '733759' and runPeriod != 'SPS_10_2016':
        sys.exit(periodError('SPS_10_2016'))

    if runNumber >= '736500' and runNumber <= '736575' and runPeriod != 'SPS_09_2017':
        sys.exit(periodError('SPS_09_2017'))


# -----------------------------------------------------------------------------
def scp(runNumber, fName, serverName, serverPath, localPath):
    ''' Download file from serverName:serverPath to localPath
    '''
    print(
        "[{}] - Downloading run '{}' from {}:{}".format(
            os.path.basename(__file__), runNumber, serverName, serverPath
        )
    )
    fName.format(runNumber)
    scpPath = serverName + ":" + serverPath + fName.format(runNumber)
    print(scpPath)
    try:
        subprocess.check_call(['scp', scpPath, localPath])
    except subprocess.CalledProcessError:
        sys.exit("[{}] - Something wrong happened while downloading with scp: return code".format(
            os.path.basename(__file__)))


# -----------------------------------------------------------------------------
def createJob(executable, args=[], name='', comment='', backend='Local', backendCE='', voms=''):
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
        except NameError:  # ganga > 6.3 no longer has the gridProxy credentials system
            print("using new cred system")
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
def makeArxiv(folderToArxiv, arxivName, excludeList=[]):
    '''Make a tgz arxiv from folderToArxiv skipping excludeList
    '''
    cmd = "tar"
    excludeList.append(arxivName)  # Always exclude the arxivName id the arxiv
    for f in excludeList and excludeList:
        cmd += " --exclude " + f
    cmd += " -zcvf " + arxivName + " " + folderToArxiv
    # print(cmd + [a for a in args])
    realCmd = cmd.split(' ')
    print('Making arxiv of folder \'{}\' with cmd \'{}\''.format(folderToArxiv, realCmd))
    try:
        subprocess.check_output(realCmd)
    except subprocess.CalledProcessError as e:
        sys.exit(e.output)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
def main():
    '''
    '''
    scriptName = os.path.basename(__file__)  # For logging clarity

    runNumber = None
    runListArg = None
    configFile = "configFile"  # Default config file if none is given on cli

    # --- Parse CLI arguments
    if len(sys.argv) > 1:
        # --- Load configuration File
        configFile = sys.argv[1]
        try:
            exec("import {0} as conf".format(configFile))
        except (ImportError, SyntaxError):
            sys.exit("[{}] - Cannot import config file '{}'".format(scriptName, configFile))
        # --- /Load configuration File
        if len(sys.argv) > 2:
            # --- Load runList
            runListArg = sys.argv[2]

    elif conf.runList:  # runNumber configured in configFile
        runList = conf.runList
        print("[{}] : Running with configuration file '{}' on runs '{}'".format(scriptName, configFile, runList))
    else:
        sys.exit("Please give : configFile - runNumber(s)(optional if set up in configFile)")
    # --- /Parse CLI arguments

    checkForConfigFile(conf)
    # --- Load List of runs
    if runListArg is None:  # if no runs on CLI, load list from configFile
        try:
            runList = conf.runList
        except AttributeError:
            sys.exit("[{}] - No runs specified at command line or in configFile...exiting".format(scriptName))
    else:
        runList = runListArg.split(',')
    # --- /Load List of runs

    # TODO modify the logic of the loop to make subjobs from the runList on the grid
    # For now each run will make a job! Looooot of clutter

    if conf.runOnGrid is True:
        makeArxiv('./', 'processor.tgz', ['./.git', './.vscode', './build', './lib', './*.pyc', './*.yml'])

    for run in runList:
        # Check if runNumber match given period in configFile
        checkPeriod(str(run), conf.runPeriod, configFile)
        runNumber = int(run)
        #   Create list of files to run on
        if conf.runOnGrid is False:
            print("[{}] - Looking for files to process for run '{}' in '{}'... ".format(scriptName, runNumber, conf.inputPath))
            if os.path.exists(conf.inputPath) is False:
                sys.exit("\n[{}] - Folder '{}' does not exist...exiting".format(scriptName, conf.inputPath))

            inputDataFileList = [f for f in os.listdir(conf.inputPath) if str(runNumber) in f]
            if not inputDataFileList:
                doScp = raw_input("[{}] - No file found...download it from '{}' ? (y/n)".format(scriptName, conf.serverName))
                if doScp == 'y':
                    scp(runNumber, conf.inputFile, conf.serverName, conf.serverDataPath, conf.inputPath)
                else:
                    sys.exit('Exiting')
        else:
            gridInputPath = '/grid/' + conf.voms + '/' + conf.inputPath
            try:
                inputDataFileList = [
                    f for f in subprocess.check_output(['lfc-ls', gridInputPath]).splitlines() if str(runNumber) in f
                ]
            except (subprocess.CalledProcessError, OSError):
                sys.exit("\n[{}] - Folder '{}' not found...exiting".format(scriptName, gridInputPath))

        if conf.runOnGrid is False:
            outputFile = conf.outputPath + conf.outputFile.format(runNumber)
            marlinLib = conf.processorPath + 'lib/' + conf.marlinLib
            inputDataFileList = [conf.inputPath + f for f in inputDataFileList]
        else:  # file is uploaded to the local folder on the WorkerNode
            marlinLib = 'lib/' + conf.marlinLib
            outputFile = conf.outputFile.format(runNumber)
        conf.glob.LCIOInputFiles = ' '.join(inputDataFileList)

        conf.marlinProc.LCIOOutputFile = outputFile + ".slcio"
        conf.marlinProc.ROOTOutputFile = outputFile + ".root"

        if conf.runOnGrid is False:
            print("[{}] - output file : {}.slcio".format(scriptName, outputFile))
            print("[{}] - MARLIN_DLL: {}".format(scriptName, marlinLib))

            if os.path.exists("{0}.slcio".format(outputFile)) is True:
                sys.exit("[{}] - OutputFile already present : Delete or move it before running again...exiting".format(
                    scriptName))

        # Looking for or generating xml Geometry File - Not necessary for streamout
        try:
            if os.path.exists(conf.geomPath + conf.geomFile) is False:
                print("[{}] - No geometry file found, creating one from database for period '{}'...".format(
                    scriptName, conf.runPeriod))
                generateXMLGeometryFile(conf.runPeriod, conf.geomPath + conf.geomFile)
                print("[{}] - No geometry file found, creating one from database for period '{}'...OK -> '{}' ".format(
                    scriptName, conf.runPeriod, conf.geomPath + conf.geomFile))
            else:
                print("[{}] - Found geometry file '{}'".format(scriptName, conf.geomPath + conf.geomFile))
        except AttributeError:
            try:
                if 'Streamout' in conf.processorType:
                    pass
            except AttributeError:
                sys.exit("[{}] - No geometryFile attribute found in configFile '{}'".format(scriptName, configFile))

        # Printing modification to xml file
        if conf.runOnGrid is False:
            print("\n[{}] ========================".format(scriptName))
            print("[{}] --- Dumping modified xml parameters: ".format(scriptName))
            for par, value in vars(conf.glob).items():
                if par != 'name':
                    print("[{}] \t\t{}:\t{}".format(scriptName, par, value))
            for par, value in vars(conf.marlinProc).items():
                if par != 'name':
                    print("[{}] \t\t{}:\t{}".format(scriptName, par, value))
            print("[{}] ========================\n".format(scriptName))

        # Marlin configuration
        # marlinCfgFile = conf.marlinCfgFile.format(conf.processorPath, runNumber)
        marlinCfgFile = conf.marlinCfgFile.format(runNumber)
        marlin = Marlin()
        marlin.setXMLConfig(conf.xmlFile)
        marlin.setLibraries(marlinLib)
        marlin.setILCSoftScript(conf.initILCSoftScript)
        setCliOptions(marlin, conf.glob)  # TODO: Move to Marlin.py
        setCliOptions(marlin, conf.marlinProc)
        marlin.writeConfigFile(marlinCfgFile)

        # Running locally
        if conf.runOnGrid is False:
            log = open(conf.logFile.format(conf.logPath, runNumber), "w", 1)  # line-buffered
            print("\n[{}] ========================".format(scriptName))
            print('[{}] --- Running Marlin...'.format(scriptName))
            print("[{}] --- Output is logged to '{}'".format(scriptName, log))
            beginTime = time.time()
            if conf.logToFile is True:
                subprocess.call(['python', 'run_marlin.py', conf.marlinCfgPath + marlinCfgFile], stdout=log, stderr=log)
            else:
                subprocess.call(['python', 'run_marlin.py', conf.marlinCfgPath + marlinCfgFile])

            print('[{}] - Running Marlin...OK - '.format(scriptName), end='')
            try:
                elapsedTime(beginTime)
            except ValueError:
                print("Can't print time...")
            print("[{}] ========================\n".format(scriptName))

            # print ('[{}] - Removing xmlFile...'.format(scriptName), end='')
            # subprocess.Popen(["rm", xmlFile])

        # Add gridInfo to Marlin configuration file + make Dic
        else:
            if gangaNotFound is True:
                sys.exit('Ganga python binding not found, can not run on grid')
            with open(conf.marlinCfgPath + marlinCfgFile, 'r') as ymlfile:
                cfg = yaml.load(ymlfile)
            gridSection = {}
            cfg['Grid'] = gridSection
            # gridSection['downloader'] = conf.gridDownloader
            # gridSection['uploader'] = conf.gridUploader
            gridSection['LCG_CATALOG_TYPE'] = conf.lcg_catalog_type
            gridSection['LFC_HOST'] = conf.lfc_host
            gridSection['inputFiles'] = conf.gridInputFiles
            with open(conf.marlinCfgPath + marlinCfgFile, 'w') as ymlfile:
                ymlfile.write(yaml.dump(cfg, default_flow_style=False))

            try:
                print("[{}] --- Submiting Job ... ".format(scriptName))
                # Navigate jobtree, if folder doesn't exist create it
                treePath = conf.runPeriod + '/' + scriptName[:-3]  # Remove .py at end of scriptName
                # print("[{}] --- treePath: {}".format(scriptName, treePath))
                # if jobtree.exists(treePath) is False: # Always returns true...
                jobtree.cd('/')  # make sure we are in root folder
                try:
                    jobtree.cd(treePath)
                except TreeError:
                    try:
                        jobtree.mkdir(treePath)
                    except TreeError:  # mkdir should write all missing folder if any....apparently not true
                        print("WhatThe?")
                        jobtree.mkdir(conf.runPeriod)
                        jobtree.mkdir(treePath)
                    jobtree.cd(treePath)

                # Update ganga configuration for eos access
                # config.Output.MassStorageFile['defaultProtocol'] = 'root://eosuser.cern.ch'
                # config.Output.MassStorageFile['fileExtensions'] = ["*.slcio", "*.root"]
                # config.Output.MassStorageFile['uploadOptions']['path'] = conf.eos_home

                # Print it
                # print(config.Output.MassStorageFile['defaultProtocol'])
                # print(config.Output.MassStorageFile['uploadOptions']['cp_cmd'])
                # print(config.Output.MassStorageFile['uploadOptions']['ls_cmd'])
                # print(config.Output.MassStorageFile['uploadOptions']['mkdir_cmd'])
                # print(config.Output.MassStorageFile['uploadOptions']['path'])

                inputFiles = []
                for f in conf.gridInputFiles:
                    inputFiles.append(LocalFile(f))
                inputFiles.append(LocalFile(conf.marlinCfgPath + marlinCfgFile))
                # print('inputFiles:\n', inputFiles)

                # TODO: Check that files already exists on grid before submitting
                inputData = []
                for item in [inputDataFileList]:
                    inputList = []
                    for f in item:
                        inputList.append(
                            LCGSEFile(
                                namePattern=f,
                                se_rpath=conf.inputPath,
                                lfc_host=conf.lfc_host,
                                se=conf.SE,
                                credential_requirements=VomsProxy(vo=conf.voms)))
                inputData = GangaDataset(treat_as_inputfiles=True, files=[f for f in inputList])

                j = createJob(
                    executable='generateEnv.sh',
                    args=['run_marlin.py', marlinCfgFile],
                    comment=conf.processorType + ' ' + conf.runPeriod,
                    name=str(runNumber),
                    backend=conf.backend,
                    backendCE=conf.CE,
                    voms=conf.voms)
                # j.outputfiles = [
                #     MassStorageFile(namePattern="*.*", outputfilenameformat='GridOutput/Streamout/{jid}/{fname}')
                # ]
                # NamePattern for LCGSEFile needs to be exact and not use wildcard otherwise the guid of the file is
                # not set at the end of the job https://github.com/ganga-devs/ganga/issues/1186
                # Only the first file will have its locations attribute updated (others fail because of the \n at the begining of the line...)
                # TODO: make script to recover name from the files that were not uploaded and make the proper alias
                j.outputfiles = [
                    LCGSEFile(
                        namePattern=outputFile + ".slcio",
                        se_rpath=conf.outputPath,
                        lfc_host=conf.lfc_host,
                        se=conf.SE,
                        credential_requirements=VomsProxy(vo=conf.voms)),
                    LCGSEFile(
                        namePattern=outputFile + ".root",
                        se_rpath=conf.outputPath,
                        lfc_host=conf.lfc_host,
                        se=conf.SE,
                        credential_requirements=VomsProxy(vo=conf.voms))
                ]
                j.inputfiles = inputFiles
                j.inputdata = inputData

                # Add a postprocessor that will move the outputfiles on the grid to their proper location
                # This needs to be present in the same folder
                j.postprocessors.append(CustomChecker(module='GangaCheckerMoveFileOnGrid.py'))
                jobtree.add(j)

                queues.add(j.submit)
                # j.submit()
                print("\n[{}] ... submitting job done.\n".format(scriptName))

            except IncompleteJobSubmissionError:
                sys.exit("[{}] --- Failed to submit job ".format(scriptName))


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
if __name__ == '__main__':
    main()
