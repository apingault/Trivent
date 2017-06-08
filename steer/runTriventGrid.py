from Ganga import *


# runPeriod = "SPS_12_2014"
# runPeriod = "PS_06_2015"
# runPeriod = "SPS_04_2015"
runPeriod = "SPS_10_2015"
# runPeriod = "SPS_06_2016"
# runPeriod = "SPS_10_2016"
scriptDir = "/eos/user/a/apingaul/CALICE/Software/Trivent"
scriptName = scriptDir + "/steer/runTrivent.sh"
marlinDLL = scriptDir + '/lib/libTrivent.so'
# geometryFile = scriptDir + '/DifGeom/m3_bonneteau.xml'
geometryFile = scriptDir + '/DifGeom/m3_bonneteau_avril2015.xml'

carefulLoaderDir = "/eos/user/a/apingaul/CALICE/script/"
lfcHost = "grid-lfc.desy.de"
# lfcSE = "srm-public.cern.ch"
lfcSE = "lyogrid06.in2p3.fr"


dataDir = "/eos/user/a/apingaul/CALICE/CALICE/Data/" + runPeriod + "/"
gridDir = "/grid/calice/SDHCAL/pingault/TB/CERN/" + runPeriod + "/"
outputDir = dataDir + "Trivent/"
outputGridDir = gridDir + "Trivent/"


inputFiles = [carefulLoaderDir + files for files in ['carefulUpload.sh', 'carefulDownload.sh']]
inputFiles.append(marlinDLL)
inputFiles.append(LocalFile(geometryFile))
inputFiles.append("~/init_ilcsoft.sh")

j = Job()
j.application = Executable(exe=File(scriptName), args=['arg1', 'arg2'])

# j.backend='Local'
j.comment = "Trivent " + runPeriod
j.inputfiles = inputFiles


j.backend='CREAM'
# j.backend.CE='lyogrid07.in2p3.fr:8443/cream-pbs-calice'
j.backend.CE='grid-cr0.desy.de:8443/cream-pbs-desy'

# j.outputfiles = [LocalFile(namePattern="*.root", localDir=outputDir)]
j.outputfiles = [LCGSEFile(namePattern="*.root", se_rpath=outputGridDir, lfc_host=lfcHost, se=lfcSE)] # localDir=outputGridDir


#SPS_12_2014
# runList=[726373]
# ,726407,726408,726409,726411,726412,726413,726414,726289,726301,726307,726310,726328]

#SPS_04_2015
# runList=[728338]

#PS_06_2015
# runList=[728645]

#SPS_10_2015
# runList=[730677]
runList=[730657, 730648, 730655, 730705, 730709, 730677, 730668, 730713]
# runList=[733656]
# runList = [733598, 733615] 
# 733626,
#  733628,
#  733631, 733632, # 733633 FAILED,
#  733636,
#  733637, 733641, 733642, 733643, 733644, 733645, 733646, 733647, 733650, 733653, 733654, 733655,
#  733656,
#  733658, 733659, 733660, # 733665 FAILED open or read 
#  733675, 
#  733678, 733679, # 733680 FAILED 
#  733683, 733686, 733688, # 733689 FAILED 
 # 733692, 733693, 733696, 733698, 733699, 733700, 733701, 733702, 733705, 733707, 733708, 733710, 733711, # 733718 FAILED 
 # 733719, 733720, 733722, 733723, 733724, 733725, 733728, 733738, 733740, 733741, 733742, 733743, 733748,
#  733750,
 # 733754, 733756, 733757, 733758, 733759]
# runList=[733758, 733759]

# SPS_06_2016
# runList = [732695, 732773, 732777, 732778, 732725, 732728, 732729, 732774, 732867, 732817, 732818, 732820, 732826, 732872, 732875, 732877, 732878, 732865, 732882, 732889]


# SPS_10_2016
# runList = [733724, 733728, 733742, 733743, 733748, 733750, 733754, 733756, 733711, 733718, 733654, 733655, 733656, 733660, 733665, 733675, 733680, 733683, 733686, 733689, 733693, 733696]

# runList=[174,176,180,183,197,199,204,223,
     # 401,403,405,406,408,410,412]
par = []
args = [[str(i), runPeriod] for i in runList]
a = len(args)
par = [args[i] for i in range(a)]
par
s = ArgSplitter()
s.args = par
j.splitter = s
queues.add(j.submit)
