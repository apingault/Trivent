import sys
import os
from subprocess import *
import xml.etree.cElementTree as ET

print "==================================================="

call("date") 
call (["uname", "-a"])

print "\n================ Exporting MarlinDLL================"

os.environ['MARLIN_DLL'] = "/Users/antoine/CALICE/DataAnalysis/triventArnaud/lib/libTrivent.dylib"

print "\n================ Creating Steering File ================"
if len(sys.argv)<2:
	print "\n *** You did not provide a run number to analyse ***"
	print "\nAvailable Files"
	call(["ls ../data/ | grep SO.slcio", "-a"], shell=True)
	sys.exit(0)

# Global
_runNumber = sys.argv[1]
_maxRecordNumber=5
_skipEventNumber = 1
_verbosity = "MESSAGE"
_supressCheck = "false"

# Trivent
_setupGeometry = "/Users/antoine/CALICE/DataAnalysis/triventArnaud/DifGeom/m3_bonneteau_avril2015.xml"
_gainCorrection = "false"
_DIFMapping = "./DifGeom/Slot1_39_Geom.txt"
_electronicNoiseCut=100000
_layerGap=2.8
_layerCut=7
_noiseCut=7
_timeWin=2
# _time2prevEventCut
# _zCut
_cerenkovWindow = 15
_cerenkovLength = 1


tree = ET.parse('steer/LCIO_runProcess.xml')
root = tree.getroot()

# Global Variables
glob = root.find('global')

if not glob:
	print "\n *** Global section not found ***"

else:
	LCIOInputFiles = glob.find(".//*[@name='LCIOInputFiles']")
	LCIOInputFiles.text = "../data/DHCAL_" + str(_runNumber) + "_SO.slcio"
	if not os.path.isfile(LCIOInputFiles.text):
		print "\n *** File " + LCIOInputFiles.text + " not found ***"
		print "\nAvailable Files:"
		call(["ls ../data/ | grep SO.slcio", "-a"], shell=True)
		sys.exit(0)

	MaxRecordNumber = glob.find(".//*[@name='MaxRecordNumber']")
	MaxRecordNumber.text = str(_maxRecordNumber)

	SkipEventNumber = glob.find(".//*[@name='SkipNEvents']")
	SkipEventNumber.text = str(_skipEventNumber)

	Verbosity = glob.find(".//*[@name='Verbosity']")
	Verbosity.text = _verbosity

	SupressCheck = glob.find(".//*[@name='SupressCheck']")
	SupressCheck.text = _supressCheck

# Trivent Variables
triventProcessor = root.find('processor', ".//*[@name='MyTriventProc']")
if not triventProcessor: 
    print "Process  not found, or element has no subelements"

else:
	DIFMapping = triventProcessor.find(".//*[@name='DIFMapping']")
	DIFMapping.text = _DIFMapping

	SetupGeometry = triventProcessor.find(".//*[@name='SetupGeometry']")
	SetupGeometry.text = _setupGeometry

	GAIN_CORRECTION_MODE = triventProcessor.find(".//*[@name='GAIN_CORRECTION_MODE']")
	GAIN_CORRECTION_MODE.text = str(_gainCorrection)

	ElectronicNoiseCut = triventProcessor.find(".//*[@name='ElectronicNoiseCut']")
	ElectronicNoiseCut.text = str(_electronicNoiseCut)

	LayerGap = triventProcessor.find(".//*[@name='LayerGap']")
	LayerGap.text = str(_layerGap)

	LayerCut = triventProcessor.find(".//*[@name='LayerCut']")
	LayerCut.text = str(_layerCut)

	NoiseCut = triventProcessor.find(".//*[@name='NoiseCut']")
	NoiseCut.text = str(_noiseCut)

	TimeWin = triventProcessor.find(".//*[@name='TimeWin']")
	TimeWin.text = str(_timeWin)

	LCIOOutputFile = triventProcessor.find(".//*[@name='LCIOOutputFile']")
	LCIOOutputFile.text = "./output/DHCAL_" + str(_runNumber) + "_Trivent.slcio"



tree.write('LCIO.xml')

print "\n================  Running Trivent ======================="

call("Marlin LCIO.xml", shell=True)
call("rm LCIO.xml", shell=True)



