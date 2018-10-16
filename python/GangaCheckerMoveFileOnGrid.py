#!/usr/bin/env python2
from Ganga.Utility.GridShell import getShell
from Ganga.Utility.Config import getConfig
from Ganga.Core.exceptions import ProtectedAttributeError
import time
import os

'''
lcgse TDHCAL_Antoine_736572POS.slcio TDHCAL_Antoine_736572POS.slcio ->guid:382bf599-77b1-49e7-8ead-ffba57bea9e1
\nlcgse TDHCAL_Antoine_736572POS.root TDHCAL_Antoine_736572POS.root ->guid:4a7697ce-f654-4c6c-8f9f-27d7f12f5235
\n
'''


def fileExists(f, fName):
    cmd = 'lcg-ls --vo {} lfn:{}'.format(f.credential_requirements.vo, fName)
    # print 'lsCmd: ', cmd
    (exitcode, output, m) = getShell(f.credential_requirements).cmd1(cmd, capture_stderr=True)
    if exitcode == 0:  # File exists
        return True
    elif 'No such file or directory' in output:
        return False
    else:
        f.failureReason += '\n' + " The file can't be uploaded because of lsCmd output : {}".format(output)
        return True


def makeGridAlias(f, fName):
    if not f.locations:
        output = "The file can't be uploaded because no guid found for '{}'".format(fName)
        f.failureReason += '\n' + output
        return 1, output
    cmd = 'lcg-aa --vo {} {} lfn:{}'.format(f.credential_requirements.vo, f.locations, fName)
    # print 'aliasCmd: ', cmd
    (exitcode, output, m) = getShell(f.credential_requirements).cmd1(cmd, capture_stderr=True)
    if exitcode != 0:
        f.failureReason += "\nThe file can't be uploaded because of {}".format(output)
    return exitcode, output


def fixPostProcessFile(job):
    postprocessLocationsPath = os.path.join(job.outputdir, getConfig('Output')['PostProcessLocationsFileName'])
    postprocesslocations = None
    with open(postprocessLocationsPath, 'r') as postprocesslocations:
        all_lines = postprocesslocations.readlines()
    all_lines = [line.replace('\\n', '\n') for line in all_lines]
    with open(postprocessLocationsPath, 'w') as postprocesslocations:
        postprocesslocations.write(''.join(all_lines))


def check(j):
    ''' postProcessor for gangaJob. Used to move file to their proper location on the grid
        They go by default to /grid/calice/generated/date/file-SomeRandomNumber
        Function takes the guid of the file, make an alias (lcg-aa) to its proper location
        on the grid
    '''
    outputFiles = [file for file in j.outputfiles]
    good = True
    for f in outputFiles:
        fName = '/grid/{}/{}{}'.format(f.credential_requirements.vo, f.se_rpath, f.namePattern)
        # TODO: check rules in ganga DiracFile for renaming!
        # print 'fName: ', fName
        if fileExists(f, fName):
            fName += time.strftime("-%Y-%m-%d_%H-%M-%S")

        # lcgse not able to process multiple output files due to typo in postprocesslocations file
        # Try fixing the typo to get around
        if not f.locations:
            fixPostProcessFile(j)
            f.setLocation()
        rc, output = makeGridAlias(f, fName)
        if rc != 0:
            good = False
    return good
