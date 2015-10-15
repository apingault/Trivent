#!/bin/bash

echo "==================================================="

date 
uname -a

#echo "============= Trivent parameters =================="

RUN=$1
GAIN_CORRECTION_MODE=false
electronic_noise_cut=100000
layerGap=2.8
LayerCut=7
noiseCut=7
timeWin=2

echo "================  ILCSoft - ======================="

source /home/steen/ilcsoft/v01-17-05/init_ilcsoft.sh
export MARLIN_DLL=/home/steen/trivent/trunk/lib/libTrivent.so

echo "========  steering file for Marlin ============="

cat >LCIO.xml<<EOF
<marlin>
  <!--##########################################
      #                                        #
      #     Steering file for Trivent_beta     #
      #                                        #
      ##########################################-->
  <execute>
    <processor name="MySDHCAL_RawData_Processor"/>
  </execute>
  
  
  <global>
    <parameter name="LCIOInputFiles">
     /data/NAS/December2014/DHCAL_${RUN}_I0_0.slcio
    </parameter>
    <parameter name="MaxRecordNumber" value="100"/>
    <parameter name="SkipNEvents" value="0"/>
    <parameter name="Verbosity" type="string"> MESSAGE </parameter>
    <parameter name="SupressCheck" value="false"/>
  </global>
  
  <processor name="MySDHCAL_RawData_Processor" type="SDHCAL_RawData_Processor">
   <!--SDHCAL_RawData_Processor prints info on the Raw data-->
   <!--verbosity level of this processor ("DEBUG0-4,MESSAGE0-4,WARNING0-4,ERROR0-4,SILENT")-->
   <parameter name="Verbosity" type="string">DEBUG </parameter>
   <!--parameter name="Verbosity" type="string">MESSAGE </parameter-->
   <!--XDAQ produced collection name-->
   <parameter name="XDAQCollectionName" type="string" lcioInType="LCGenericObject">RU_XDAQ </parameter>
   <!--Name of output collection containing raw calorimeter hits-->
   <parameter name="OutputRawCaloHitCollectionName" type="string" lcioOutType="RawCalorimeterHit">DHCALRawHits </parameter>
   <!--Turn ON/OFF debug mode : Warning Debug mode uses assert and may crash the application-->
   <parameter name="DebugMode" type="bool"> true </parameter>
  </processor>
  
</marlin>

EOF
echo "========  RUNNING TRIVENT RECO ================="

Marlin LCIO.xml
rm LCIO.xml
