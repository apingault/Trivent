#!/bin/bash


if [ -z "$ILCSOFT" ]; then
   echo " You should source an init_ilcsoft.sh script"
   exit
fi



generateRootDict()
{
    echo "============================================================"
    echo "=                 Generate root dictionnary                ="
    echo "============================================================"
    ORIGDIR=${PWD}
    cd ${PWD}/SDHCALReco/src
    if [ -e rootDict.cc ]; then
	echo "The file rootDict.cc is already exist .."
	rm rootDict.*
	rootcint rootDict.cc -c -I../../EventBuilder/include/ ../include/TCaloHitRoot.hh ../include/TRawCaloHitRoot.hh
    else
	echo "Create a new dictionary rootDict.cc .."	
	rootcint rootDict.cc -c -I../../EventBuilder/include/ ../include/TCaloHitRoot.hh ../include/TRawCaloHitRoot.hh
    fi
    ls -lthr

    cd $ORIGDIR
}


CreateBuildDir()
{
    if [ -d build ]; then
	rm -r build
    fi
    mkdir build
}

CreateMakefile()
{
    echo "============================================================"
    echo "=                 Build Makefile for compilation           ="
    echo "============================================================"    
    echo "  Create a new Makefile  .. "
    cmake -C ${ILCSOFT}/ILCSoft.cmake ..
    ls -lthr  
}

Compile()
{
    echo "============================================================"
    echo "=                      Compilation                         ="
    echo "============================================================"
    make install
}

#The script

WHAT=$1

if [ "$WHAT" == "Full" ]; then
    generateRootDict
    CreateBuildDir
    cd build
    CreateMakefile
    cd ..
fi

cd build
Compile
echo " Compilation done"


