#!/usr/bin/env bash
export PATH=/cvmfs/sft.cern.ch/lcg/contrib/gcc/4.9.3/x86_64-slc6/bin:/cvmfs/sft.cern.ch/lcg/releases/LCG_87/Python/2.7.10/x86_64-slc6-gcc49-opt/bin:${PATH}
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_87/Python/2.7.10/x86_64-slc6-gcc49-opt/lib/:${LD_LIBRARY_PATH}
wget https://bootstrap.pypa.io/get-pip.py && python get-pip.py --user
cd ../.local/bin
./pip install pyyaml --user
cd -
echo " ================================================== "
echo " --- Building sources "
tar -xvzf processor.tgz
source init_ilcsoft.sh # ilcsoft script
mkdir build; cd build; cmake -C $ILCSOFT/ILCSoft.cmake ..; make install -j4; cd -
echo " --- Done building sources\n"
echo " ================================================== "
echo " --- Runnning script " $1 " with option " $2
ls -alhtr
python $1 $2
echo " --- Done running script"
echo " ================================================== "
ls -alhtr
echo " ================ "
echo " --- finished --- "
