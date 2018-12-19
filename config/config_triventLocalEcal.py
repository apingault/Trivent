#!/usr/bin/env python2
'''

'''

import sys
import platform

runList = [
    744211, 744214, 744218, 744221, 744226, 744230, 744233, 744237, 744242, 744246, 744249, 744254, 744258, 744263,
    744269, 744273, 744278, 744283, 744317
]

maxNumberOfTrigger = 0  # Max Number of trigger to analyse, 0 = all
numberOfTriggerToSkip = 1  #

###
#  Paths
runPeriod = 'SPS_09_2018'
eos_sdhcal = '/eos/project/s/sdhcal/data/'

processorPath = '/eos/user/a/apingaul/CALICE/Software/Trivent/'
dataPath = '/eos/user/a/apingaul/CALICE/Data/' + runPeriod + '/'
inputPath = eos_sdhcal + runPeriod + '/Raw/'  # input folder with slcio
outputPath = dataPath + 'Ecal/'  # output folder for slcio/root

inputFile = 'DHCAL_{}_I0_0.slcio'  # .format(runNumber)
outputFile = 'TDHCAL_Ecal_{}'  # .format(runNumber)

plotPath = dataPath + 'Plots/'
logPath = dataPath + 'Logs/'

ilcSoftVersion = "v01-19-05"
ilcSoftPath = '/opt/ilcsoft/'

runOnLxplus = False  # Automatically adapt ilcsoft paths when on lxplus
if platform.node().find("lxplus") != -1:
    runOnLxplus = True

if runOnLxplus:
    ilcSoftVersion = 'current/CI_gcc'
    ilcSoftPath = '/cvmfs/clicdp.cern.ch/iLCSoft/builds/'
initILCSoftScript = ilcSoftPath + ilcSoftVersion + "/init_ilcsoft.sh"
'''
    # Should not need to edit part below
'''

xmlFile = 'config/TriventProcessor.xml'
geomFile = 'sdhcal_SPS_09_2018.json'
geomPath = processorPath + 'DifGeom/'  # folder with json geometry files


# Marlin options
class xmlOptionSection(object):
    def __init__(self, optionName):
        self.name = optionName


# Global param
###
processorType = 'Trivent'
glob = xmlOptionSection('global')
glob.MaxRecordNumber = maxNumberOfTrigger  # Max Number of event to process
glob.SkipNEvents = numberOfTriggerToSkip  # Number of event to skip
glob.LCIOInputFiles = []

# Processor param
###
marlinProc = xmlOptionSection('MyTriventProcessor')
marlinProc.SetupGeometry = geomPath + geomFile
marlinProc.TreeName = 'sdhcal'

# scp
serverName = 'lxplus'
serverDataPath = eos_sdhcal + runPeriod + '/Raw/'

runOnGrid = False
logToFile = True
logFile = '{}/triventLog_{}'  # .format(logPath, runNumber)

libExt = 'so'
if sys.platform == 'darwin':
    libExt = 'dylib'
marlinLib = "libTriventProc.{}".format(libExt)  # Marlin library for processor
marlinCfgFileName = "marlinCfg_{}.yml"  # .format(runNumber) cfgFile name written by script to properly run marlin with all variables set
marlinCfgFile = processorPath + 'config/marlinCfg/' + marlinCfgFileName  # Where to put the configFiles
