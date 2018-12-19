# Trivent - SDHCAL timing event builder

## Dependencies

You need cmake>3.5, gcc>5.0 (for c++14 features) and python>2.7
From ilcsoft you also need, Marlin, lcio, ilcutil, Root6

For ease of use you can use the clicdp installation of ilcsoft available on cvmfs at  `/cvmfs/clicdp.cern.ch/iLCSoft/builds/current/CI_gcc`
If for whatever reason this build doesn't work anymore you can use `/cvmfs/clicdp.cern.ch/iLCSoft/builds/2018-10-11/x86_64-slc6-gcc62-opt/`.

## Installation

```bash
    export ilcsoftScript=/cvmfs/clicdp.cern.ch/iLCSoft/builds/current/CI_gcc/init_ilcsoft.sh  # adapt to your needs
    git clone --branch CommonEcalSdhcal https://github.com/apingault/Trivent/
    source $ilcsoftScript
    cd Trivent; mkdir build; cd build; cmake -C $ILCSOFT/ILCSoft.cmake ..
    make install
    cd -
```

## Configuration
Edit `config/configTriventLocal.py` with the paths to data/ilcsoft script etc.


## Running

To ensure proper environment I would suggest using a virtualenv: (These steps are needed only the first time around, it should take ~2min on lxplus)

``` bash
    # From the main dir
    source $ilcsoftScript
    curl -LO https://bootstrap.pypa.io/get-pip.py
    python get-pip.py --user
    export PY_USER_BIN=$(python -c 'import site; print(site.USER_BASE + "/bin")')
    export PATH=$PY_USER_BIN:$PATH
    python -m pip install --user virtualenv
    virtualenv --no-site-packages -p $(which python) triventEnv  # Ensure you use the correct Python version not the ols system one(2.6.6)
    source triventEnv/bin/activate
    # Update python installation
    pip install pyyaml pymysql
```

Once this is set-up you can just:

```bash
    export PY_USER_BIN=$(python -c 'import site; print(site.USER_BASE + "/bin")')
    export PATH=$PY_USER_BIN:$PATH
    source triventEnv/bin/activate  # Use the virtual env
    python/pyMarlin.py config/config_triventLocalEcal  # without the .py extension at the end of the configFile
```

github sometimes messes with my files permissions, if you have a 

```bash
    bash: ./python/pyMarlin.py: Permission denied
```

Just do:

```bash
    chmod +x python/pyMarlin.py
```

Once you are done you can leave the virtualenv by calling `deactivate`


## Running Marlin directly

You can also directly run the marlin processor from the command line as usual. Edit the config file (default is in `config/TriventProcessor.xml`) and run with

```bash
    export MARLIN_DLL=lib/libTriventProc.so  # or .dylib if running on mac
    Marlin config/TriventProcessor.xml  # 
```

You can also directly specify Marlin option on the command line, this will override the value from the config file.

```bash
    Marlin --global.MaxRecordNumber="1000" config/TriventProcessor.xml
```
