# Trivent - SDHCAL timing event builder

## Dependencies

You need cmake>3.5, gcc>5.0 (for c++14 features)
From ilcsoft you also need, Marlin, lcio, ilcutil, Root6

For ease of use you can use the clicdp installation of ilcsoft available on cvmfs at  `/cvmfs/clicdp.cern.ch/iLCSoft/builds/current/CI_gcc`
If for whatever reason this build doesn't work anymore you can use `/cvmfs/clicdp.cern.ch/iLCSoft/builds/2019-07-09/x86_64-slc6-gcc7-opt/`.

## Installation

```bash
export ilcsoftScript=/cvmfs/clicdp.cern.ch/iLCSoft/builds/current/CI_gcc/init_ilcsoft.sh # adapt to your needs
git clone https://github.com/apingault/Trivent
source $ilcsoftScript
cd Trivent; mkdir build; cd build; cmake -C $ILCSOFT/ILCSoft.cmake ..
make install
cd -
```

## Running

This is a classic Marlin processor so you can run it as usual with

```bash
export MARLIN_DLL=lib/marlin_dll/libTriventProc.so  # or .dylib if running on mac
Marlin config/TriventProcessor.xml
```

or with changes to Marlin options on the command line:

```bash
Marlin  --global.MaxRecordNumber="1000" --MyTriventProc.NoiseCut="10" config/TriventProcessor.xml
```

For more flexibility in the configuration process you can use [marlinpy_sdhcal](https://gitlab.cern.ch/apingaul/marlinpy_sdhcal)

## Geometry

### Convention in the PAD/Layer frame:

* I is the side perpendicular to the DIF goes from 1->96. Origin is at the top of the layer
* J is the side with the DIF goes from 1->96. Origin is at the left of the layer.
* K is the layer number goes from 0->nLayer. Origin is at the front of the calorimeter (first layer)

When facing the calorimeter the origin is on the top left corner of the first layer:

* I goes from top to bottom from 1 to 96
* J goes from left to right from 1 to 96
* K Goes from first layer to last from 0 to nLayer

### Convention in the X/Y/Z frame:

* X is the side with the DIF/parallel to the ground. Origin is on the left side of the layer.
* Y is the side perpendicular to DIF/ground. Origin is on the bottom of the layer.
* Z is the position. Origin is at the front of the first layer + eventual shift.

When facing the calorimeter the origin is on the bottom left corner of the first layer(+- shift in the z direction):

* X goes from left to right from 5.2mm to 10004,4mm
* Y goes from bottom to top from 5.2mm to 10004,4mm
* Z Goes from first layer to last

Hit is placed in the center of a pad and start of the layer.

### Conversion from PAD/Layer geometry to x/y/z:

* Shift on the z axis  is set to +225mm by default (corresponds to convention defined in 09_2018 for ecal common testbeam, Z=0 being the front of the first ecal plate)

* X = (J - 1) * pad_dimension + pad_dimension / 2
* Y = (96 - I) * pad_dimension + pad_dimension / 2
* Z = K * layer_thickness + pos_shift
