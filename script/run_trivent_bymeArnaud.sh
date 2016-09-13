#!/bin/bash

echo "==================================================="

date 
uname -a

echo "============= Trivent parameters =================="

RUN=$1
MaxRecordNumber=5000000
GAIN_CORRECTION_MODE=false
electronic_noise_cut=100000
layerGap=2.8
LayerCut=7
noiseCut=7
timeWin=2
dataFolder="/Users/antoine/CALICE/DataAnalysis/data"

echo $GAIN_CORRECTION_MODE

echo "================  ILCSoft - ======================="

source /opt/ilcsoft/v01-17-06/init_ilcsoft.sh
export MARLIN_DLL=/Users/antoine/CALICE/DataAnalysis/triventArnaud/lib/libTriventArnaud.dylib

echo "========  steering file for Marlin ============="

FILELIST=`ls ${dataFolder}/Streamout/ | grep DHCAL_${RUN}_SO.slcio`
if
    [[ $FILELIST == "DHCAL_${RUN}_SO.slcio" ]]
then
    echo $FILELIST
# else 
#     scp acqilc@lyosdhcal10.cern.ch:/data/NAS/Avril2015/STREAMOUT_tmp/DHCAL_${RUN}_SO.slcio DHCAL_${RUN}_SO.slcio
#     FILELIST=DHCAL_${RUN}_SO.slcio
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
    <parameter name="LCIOInputFiles">
     ${dataFolder}/Streamout/DHCAL_${RUN}_SO.slcio
    </parameter>
    <parameter name="MaxRecordNumber"  value="${MaxRecordNumber}" />
    <parameter name="SkipNEvents" value="0" />
    <parameter name="Verbosity" type="string"> MESSAGE </parameter>
    <parameter name="SupressCheck" type="bool"> ${GAIN_CORRECTION_MODE} </parameter>
  </global>
  

  <processor name="MyTriventProc" type="TriventProc">
    <parameter name="HitCollectionName" type="StringVec">DHCALRawHits </parameter>
    <parameter name="DIFMapping" type="string">
      ./DifGeom/Slot1_39_Geom.txt
    </parameter>
    <parameter name="setup_geometry" type="string"> /Users/antoine/CALICE/DataAnalysis/triventArnaud/DifGeom/m3_bonneteau_avril2015.xml </parameter>
    <parameter name="GAIN_CORRECTION_MODE" type="bool">false </parameter>
    <parameter name="electronic_noise_cut" type="int">500000 </parameter>
    <parameter name="layerGap" type="double"> 2.8 </parameter>
    <parameter name="LayerCut" type="int">5 </parameter>
    <parameter name="noiseCut" type="int">7 </parameter>
    <parameter name="timeWin" type="int">2 </parameter>
    <parameter name="LCIOOutputFile" value="${dataFolder}/Trivent/slcio/TDHCAL_${RUN}.slcio"/>
  </processor>
  
  
  
</marlin>

EOF
echo "========  RUNNING TRIVENT RECO ================="

Marlin LCIO.xml
# rm LCIO.xml