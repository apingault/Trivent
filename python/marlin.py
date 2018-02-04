#!/usr/bin/python
#
import os, sys
from subprocess import call
import yaml

class Marlin(object):
    def __init__(self):
        
        self.xmlConfig = None
        self.libraries = None
        self.ilcSoftInitScript = None
        self.cliOptions = {}
    
    def setXMLConfig(self, xmlFile):
        ''' Set xml config file
        '''
        self.xmlConfig = xmlFile
        
    def setLibraries(self, librarires):
        ''' Set dll librarires
        '''
        self.libraries = librarires
        
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
        
    def writeConfigFile(self, cfgFile):
        ''' Write marlin specific config into configFile
            update file if found, create it otherwise
        '''
        try:
            with open(cfgFile, 'r') as ymlfile:
                cfg = yaml.load(ymlfile)
        except IOError:
            # print "Creating new configuration file '{0}'".format(cfgFile)
            cfg = {}
        else: # found existing config file 
            if cfg is None: # empty config file
                cfg = {}
        with open(cfgFile, 'w') as ymlfile:
            marlinSection = {}
            cfg['Marlin'] = marlinSection
            marlinSection['xmlConfig'] = self.xmlConfig
            marlinSection['libraries'] = self.libraries
            marlinSection['ilcSoftInitScript'] = self.ilcSoftInitScript
            
            if self.cliOptions is not None:
                cliOptionSubSection = {}
                marlinSection['cliOptions'] = cliOptionSubSection
                
                for key, value in self.cliOptions.items():
                    cliOptionSubSection[key] = value
                    
            ymlfile.write(yaml.dump(cfg, default_flow_style=False))


    def readConfigFile(self, cfgFile):
        ''' Read Marlin specific config in cfgFile
            Exit if cfgFile not found
        '''
        try:
            with open(cfgFile, "r") as ymlfile:
                cfg = yaml.load(ymlfile)
        except IOError:
            sys.exit("[run_marlin.py] --- ERROR: Config file '{0}' not found....exiting".format(cfgFile))            
        
        try:
            marlinSection = cfg['Marlin']
            self.xmlConfig = marlinSection['xmlConfig'] 
            self.libraries = marlinSection['libraries'] 
            self.ilcSoftInitScript = marlinSection['ilcSoftInitScript']      
            self.setCliOptions(marlinSection['cliOptions'].items())
        except KeyError, exc:
            print "Key {0} not found in cfgFile".format(exc)

    def checkConfig(self,cfgFile):
        ''' Check that core config are set
            Exit otherwise
        '''
        missingConfig = []
        for k, v in vars(self).items():
            if v is None:
                missingConfig.append(k)
            if missingConfig:
               sys.exit("[run_marlin.py] --- ERROR: Some configuration are missing from '{0}'' : {1} ... exiting".format(cfgFile,missingConfig))


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
