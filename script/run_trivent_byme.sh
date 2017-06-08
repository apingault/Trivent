#!/bin/bash

echo "==================================================="

date 
uname -a

echo "============= Trivent parameters =================="

#runPeriod="October2015"
runPeriod="June2016"

RUN=$1
if [ "$2" == "" ]
then
    Sync=""
else 
    Sync="_$2"
fi

# verbosity="MESSAGE"
verbosity="DEBUG2"

MaxRecordNumber=0
skipNEvents=0
eventPrint=10

GAIN_CORRECTION_MODE=false
electronic_noise_cut=175000
layerGap=2.8
LayerCut=7
noiseCut=10
timeWin=2
cerenkovID=3
cerenkovOffset=0
cerenkovLength=1


if [ ${RUN} -lt 732000 ]; then
  dataFolder="/Users/antoine/CALICE/DataAnalysis/data" # 2015
  echo -e "RunNumber : ${RUN} < 732000 -> 2015 data..."
else 
  dataFolder="/Volumes/PagosDisk/CALICE/data/SPS_06_2016" #2016
  echo -e "RunNumber : ${RUN} > 732000 -> 2016 data..."
fi
echo -e "\t -> looking for slcio files in $dataFolder"

# dataFolder="/Users/antoine/CALICE/DataAnalysis/data"
# dataFolder="/Volumes/PagosDisk/CALICE/data/SPS_06_2016"

triventFolder="/Users/antoine/CALICE/DataAnalysis/triventArnaud"
geomFile="${triventFolder}/DifGeom/m3_bonneteau_avril2015.xml"

outputFolder="/Users/antoine/CALICE/DataAnalysis/data"
outputFile="${outputFolder}/TDHCAL_${RUN}${Sync}_Arnaud" #root/slcio extension added in the xml file


echo -e "MaxRecordNumber : "$MaxRecordNumber
echo -e "GAIN_CORRECTION_MODE : "$GAIN_CORRECTION_MODE
echo -e "electronic_noise_cut : "$electronic_noise_cut
echo -e "layerGap : "$layerGap
echo -e "LayerCut : "$LayerCut
echo -e "noiseCut : "$noiseCut
echo -e "timeWin : "$timeWin
echo -e "cerenkovID : "$cerenkovID
echo -e "cerenkovOffset : "$cerenkovOffset
echo -e "cerenkovLength : "$cerenkovLength
echo -e "dataFolder : "$dataFolder
echo -e "triventFolder : "$triventFolder
echo -e "Sync = " $Sync

echo "================  ILCSoft - ======================="

sourceilc
export MARLIN_DLL=${triventFolder}/lib/libTriventArnaud.dylib
#export MARLIN_DLL=${triventFolder}/lib/libTrivent.so

echo "========  steering file for Marlin ============="

FILELIST=`ls ${dataFolder} | grep DHCAL_${RUN}_SO${Sync}.slcio`
if
    [[ $FILELIST == "DHCAL_${RUN}_SO${Sync}.slcio" ]]
then
    echo $FILELIST
    echo -e "input files : ${dataFolder}/DHCAL_${RUN}_SO${Sync}.slcio"
else 
  echo "No Files Found"    
# else 
#     scp -r lyosdhcal10:/data/NAS/${runPeriod}/STREAMOUT_tmp/DHCAL_${RUN}_SO${Sync}.slcio ${dataFolder}/Streamout/DHCAL_${RUN}_SO${Sync}.slcio
#     FILELIST=DHCAL_${RUN}_SO${Sync}.slcio
#     #FILELIST=`ls  /data/NAS/December2014/STREAMOUT_tmp/DHCAL_${RUN}_I0_0s.slcio`
#     echo $FILELIST
fi
#FILELIST=`ls q| grep ${RUN}`
#echo $FILELIST

cat >LCIO.xml<<EOF
<marlin>
  <!--##########################################
      #                                        #
      #     Steering file for Trivent_beta     #
      #                                        #
      ##########################################-->
  <execute>
    <processor name="MyTriventProc"/>
  </execute>
  
  
  <global>
    <parameter name="LCIOInputFiles">  ${dataFolder}/DHCAL_${RUN}_SO${Sync}.slcio  </parameter>
    <parameter name="MaxRecordNumber">  ${MaxRecordNumber}  </parameter>
    <parameter name="SkipNEvents" value="${skipNEvents}"/>
    <parameter name="Verbosity" type="string"> ${verbosity} </parameter>
    <parameter name="SupressCheck" type="bool"> ${GAIN_CORRECTION_MODE} </parameter>
  </global>
  

  <processor name="MyTriventProc" type="TriventProc">
    <parameter name="HCALCollections" type="StringVec">DHCALRawHits </parameter>
    <parameter name="DIFMapping" type="string">
      ./DifGeom/Slot1_39_Geom.txt
    </parameter>
    <parameter name="SetupGeometry" type="string"> ${geomFile} </parameter>
    <parameter name="GAIN_CORRECTION_MODE" type="bool">${GAIN_CORRECTION_MODE} </parameter>
    <parameter name="ElectronicNoiseCut" type="int"> ${electronic_noise_cut} </parameter>
    <parameter name="LayerGap" type="double"> ${layerGap} </parameter>
    <parameter name="LayerCut" type="int"> ${layerCut} </parameter>
    <parameter name="HitPerLayerCut" type="int"> ${hitPerLayerCut} </parameter>
    <parameter name="NoiseCut" type="int"> ${noiseCut} </parameter>
    <parameter name="TimeWin" type="int"> ${timeWin} </parameter>
    <parameter name="LCIOOutputFile" value="${outputFile}.slcio"/>
    <parameter name="ROOTOutputFile" value="${outputFile}.root"/>
  </processor>
  
  
  
</marlin>

EOF
echo "========  RUNNING TRIVENT RECO ================="
Marlin LCIO.xml
echo " === RootFile: ${outputFile}.root"
# rm LCIO.xml
