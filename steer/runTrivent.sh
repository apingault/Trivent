#!/bin/sh 
runLocal=false

# period="SPS_10_2016"
run=$1
period=$2
maxEvt=0
# geometryFile='m3_bonneteau.xml'
geometryFile='m3_bonneteau_avril2015.xml'
CerenkovTimeWindow=30
CerenkovDifId=3


if [ "$#" -lt 2 ]; then
  echo "Wring Arguments, usage : runTrivent.sh runNumber tbPeriod"
  return 111
fi

if [ "$runLocal" = true ]; then
  echo " Running on my Pee Aye Aime Pee"
  export MARLIN_DLL=/eos/user/a/apingaul/CALICE/Software/Trivent/lib/libTrivent.dylib
  # dataFolder="/eos/user/a/apingaul/CALICE/Data/$period"
  dataFolder="/Volumes/PagosDisk/CALICE/data/$period"
  sourceilc
else 
  echo " Running on Lx+" 
  # export MARLIN_DLL=/eos/user/a/apingaul/CALICE/Software/sdhcal_analysis/lib/libTrivent.so
  # source ~/init_ilcsoft.sh
  # dataFolder="/eos/user/a/apingaul/CALICE/Data/$period"
  export MARLIN_DLL=libTrivent.so
  dataFolder="/grid/calice/SDHCAL/pingault/TB/CERN/$period"
  source init_ilcsoft.sh
fi




# dataFolder="/eos/user/a/apingaul/CALICE/CALICE/Data/SPS_10_2016/"

echo " .... Running Trivent process, local File : $(ls -al)"

export LFC_HOST=grid-lfc.desy.de

echo "LFC_HOST = '${LFC_HOST}'"

INPUTFOLDER="${dataFolder}/Streamout"

# FILELIST=`ls ${INPUTFOLDER} | grep TDHCAL_${run}.slcio`
FILELIST=$(ls | grep "DHCAL_${run}_SO_Antoine.slcio")
if
    [[ $FILELIST == "DHCAL_${run}_SO_Antoine.slcio" ]]
    then
    echo ${FILELIST}
else
    echo "Looking on grid for files, run = '${run}' - inputGridFolder = '${INPUTFOLDER}' "
    echo "lfc-ls ${INPUTFOLDER} | grep ${run} | grep slcio "
    echo "$(lfc-ls ${INPUTFOLDER} | grep ${run} | grep slcio)"
    
    FILELIST=$(lfc-ls ${INPUTFOLDER} | grep ${run} | grep slcio)
    echo " Found : ${FILELIST}"
    if [ -z "$FILELIST" ]; then
        exit 111
    fi
    for file in ${FILELIST};
      do source carefulDownload.sh ${file} ${INPUTFOLDER}
    done
fi

echo "FILELIST = " ${FILELIST}
echo " --- Found" ${FILELIST} " in " ${INPUTFOLDER}

OUTPUTFOLDER="${dataFolder}/Trivent/"
OUTPUTFILENAME="TDHCAL_${run}_Antoine" 


  # <parameter name="LCIOInputFiles"> ${INPUTFOLDER}/${FILELIST}  </parameter>
cat > TriventLCIO.xml <<EOF
<marlin>
  <execute>
    <processor name="MyTriventProcessor"/>
  </execute>
  <global>
    <parameter name="LCIOInputFiles">
    ${FILELIST} 
    </parameter>
    <parameter name="MaxRecordNumber">${maxEvt}</parameter>
    <parameter name="SkipNEvents">0</parameter>
    <parameter name="Verbosity" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT">DEBUG0</parameter>
  </global>
  <processor name="MyTriventProcessor" type="TriventProc">
    <parameter name="Verbosity" type="string" options="DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT">DEBUG0</parameter>
    <parameter name="InputCollectionNames" type="StringVec">DHCALRawHits</parameter>
    <parameter name="OutputCollectionName" type="string">SDHCAL_HIT</parameter>
    <parameter name="LCIOOutputFile" type="string">${OUTPUTFILENAME}.slcio</parameter>
    <parameter name="ROOTOutputFile" type="string">${OUTPUTFILENAME}.root</parameter>
    <parameter name="PlotFolder" type="string">"./"</parameter>
    <parameter name="SetupGeometry" type="string"> ${geometryFile}  </parameter>
    <parameter name="GainCorrectionMode" type="bool">False</parameter>
    <parameter name="ElectronicNoiseCut" type="int">350000</parameter>
    <parameter name="LayerGap" type="double">2.8</parameter>
    <parameter name="LayerCut" type="int">7</parameter>
    <parameter name="NoiseCut" type="int">10</parameter>
    <parameter name="TimeWin" type="int">2</parameter>
    <parameter name="HasCerenkovDif" type="bool">True</parameter>
    <parameter name="CerenkovDifId" type="int">${CerenkovDifId}</parameter>
    <parameter name="CerenkovTimeWindow" type="int">${CerenkovTimeWindow}</parameter>
  </processor>
</marlin>
EOF

Marlin TriventLCIO.xml
rm TriventLCIO.xml
for file in ${FILELIST}; do
    rm file
done
echo " Marlin Done...List of local files : $(ls -al)"
echo " Uploading Files (slcio) : $(ls -al | grep slcio)"
source carefulUpload.sh ${OUTPUTFILENAME}.slcio ${OUTPUTFOLDER};
echo " Uploading Files (root) : $(ls -al | grep root)"
source carefulUpload.sh ${OUTPUTFILENAME}.root ${OUTPUTFOLDER};
