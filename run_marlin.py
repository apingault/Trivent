#!/usr/bin/python

from marlin import Marlin
import sys

if __name__ == '__main__':
    try :
        cfgFile = sys.argv[1]
    except IndexError:
        sys.exit('[run_marlin.py] --- ERROR: No configuration file provided')
    marlin = Marlin()

    marlin.readConfigFile(cfgFile)
    marlin.run(cfgFile)
