#!/usr/bin/env python2
#
import os
import sys
from subprocess import call
import yaml


class Marlin(object):
    def __init__(self):

        self.xmlConfig = None
        self.runOnGrid = False
        self.libraries = None
        self.ilcSoftInitScript = None
        self.uploadScript = None
        self.downloadScript = None
        self.outputPath = None
        self.outputFiles = []
        self.cliOptions = {}

    def setXMLConfig(self, xmlFile):
        ''' Set xml config file
        '''
        self.xmlConfig = xmlFile

    def setRunOnGrid(self, runOnGrid):
        ''' Set xml config file
        '''
        self.runOnGrid = runOnGrid

    def setLibraries(self, libraries):
        ''' Set dll libraries
        '''
        self.libraries = libraries

    def setCliOptions(self, options):
        ''' Set all options to be superseded in xml
        '''
        for k, v in options:
            self.setCliOption(k, v)

    def setCliOption(self, option, value):
        ''' Set One option to be superseded in xml
        '''
        self.cliOptions[option] = value

    def setILCSoftScript(self, script):
        ''' Set init_ilcsoft.sh script to source for proper software linking
        '''
        self.ilcSoftInitScript = script

    def setUploadScript(self, script):
        ''' Set upload script to grid
        '''
        self.uploadScript = script

    def setDownloadScript(self, script):
        ''' Set Download script to grid
        '''
        self.downloadScript = script

    def setOutputPath(self, path):
        ''' Set Download script to grid
        '''
        self.outputPath = path

    def setOutputFiles(self, files):
        ''' Set Download script to grid
        '''
        self.outputFiles = files

    def writeConfigFile(self, cfgFile):
        ''' Write marlin specific config into configFile
            update file if found, create it otherwise
        '''
        try:
            with open(cfgFile, 'r') as ymlFile:
                cfg = yaml.load(ymlFile)
        except IOError:
            # print "Creating new configuration file '{0}'".format(cfgFile)
            cfg = {}
        else:  # found existing config file
            if cfg is None:  # empty config file
                cfg = {}
        with open(cfgFile, 'w') as ymlFile:
            marlinSection = {}
            cfg['Marlin'] = marlinSection
            marlinSection['xmlConfig'] = self.xmlConfig
            marlinSection['libraries'] = self.libraries
            marlinSection['ilcSoftInitScript'] = self.ilcSoftInitScript
            marlinSection['uploadScript'] = self.uploadScript
            marlinSection['downloadScript'] = self.downloadScript
            marlinSection['outputPath'] = self.outputPath
            marlinSection['outputFiles'] = self.outputFiles

            if self.cliOptions is not None:
                cliOptionSubSection = {}
                marlinSection['cliOptions'] = cliOptionSubSection

                for key, value in self.cliOptions.items():
                    cliOptionSubSection[key] = value

            ymlFile.write(yaml.dump(cfg, default_flow_style=False))

    def readConfigFile(self, cfgFile):
        ''' Read Marlin specific config in cfgFile
            Exit if cfgFile not found
        '''
        try:
            with open(cfgFile, "r") as ymlFile:
                cfg = yaml.load(ymlFile)
        except IOError:
            sys.exit("[run_marlin.py] --- ERROR: Config file '{0}' not found....exiting".format(cfgFile))

        try:
            marlinSection = cfg['Marlin']
            self.xmlConfig = marlinSection['xmlConfig']
            self.libraries = marlinSection['libraries']
            self.ilcSoftInitScript = marlinSection['ilcSoftInitScript']
            self.uploadScript = marlinSection['uploadScript']
            self.downloadScript = marlinSection['downloadScript']
            self.outputPath = marlinSection['outputPath']
            self.outputFiles = marlinSection['outputFiles']
            self.setCliOptions(marlinSection['cliOptions'].items())
        except KeyError as exc:
            print "Key {0} not found in cfgFile".format(exc)

    def checkConfig(self, cfgFile):
        ''' Check that core config are set
            Exit otherwise
        '''
        missingConfig = []
        for k, v in vars(self).items():
            if v is None:
                missingConfig.append(k)
            if missingConfig:
                sys.exit("[run_marlin.py] --- ERROR: Some configuration are missing from '{0}'' : {1} ... exiting".format(cfgFile, missingConfig))

    def run(self, cfgFile):
        self.checkConfig(cfgFile)
        # source ilcsoft for proper environment
        cmd = "source {0}; ".format(self.ilcSoftInitScript)
        cmd += "Marlin "

        # add eventual replacement options to configFile
        # print self.cliOptions
        if self.cliOptions:
            for key, value in self.cliOptions.items():
                cmd += '--{0}="{1}" '.format(key, value)

        # Add configuration file
        cmd += self.xmlConfig
        print "[run_marlin.py] --- running Marlin with cmd '{0}'".format(cmd)
        return call(cmd, env=dict(os.environ, MARLIN_DLL=self.libraries), shell=True)

    def uploadFiles(self):
        if self.runOnGrid:
            for f in self.outputFiles:
                cmd = "source {} {} ./ {}".format(self.uploadScript, f, self.outputPath)
                print "[marlin.py] --- uploading {} to {} with cmd '{}'".format(f, self.outputPath, cmd)
                call(cmd, env=dict(os.environ), shell=True)
