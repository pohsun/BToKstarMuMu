#!/usr/bin/env python

### Remark optparse is replaced with argparth in python 2.7.

import sys, os, shutil, re, subprocess, time
from   optparse import OptionParser
# import ROOT

pbsBatchTemplate = """#!/usr/bin/env bash

#PBS -V
#PBS -j oe
#PBS -q cms
#PBS -N JOBNAME
#PBS -d JOBDIR
#PBS -o JOBSTDOUT
#PBS -e JOBSTDERR
#PBS -l walltime=00:20:00

if [ ! -e JOBDIR/pbsDone ] && [ ! -e JOBDIR/pbsRuntime ]; then
    cd JOBDIR
    touch pbsRuntime
    CMD
    /bin/rm pbsRuntime
fi

touch JOBDIR/pbsDone

"""

def main():
    # parse all options
    parser = OptionParser()
    parser.add_option('-f', '--fit', dest='fit', help='The pathname to fit function', default='../fit', type='string')
    parser.add_option('-b', '--bin', dest='q2bin', help='The q2 bin to be processed', default=None, type='int')
    # parser.add_option('-o', '--output',    dest='outpath' ,      help='output path',                         default='./',    type='string')
    (opt, args) = parser.parse_args()

    workBin = opt.q2bin
    if workBin == None :
        print "ERROR   : q2 bin must be given using '-b' option."
        return 1

    fitFunc         =os.path.abspath(opt.fit)
    binfolder       =os.path.abspath("../limit/bin"+str(workBin))
    # binfolder       =os.path.abspath("../limit/validation/bin"+str(workBin))
    iwspacepath     =os.path.abspath("../wspace_ANv16/")
    idatacardpath   =os.path.abspath("../datacard_ANv16/")

    # Create the list of headfolders
    headfolders = []
    headfolderPatternAfb = re.compile("^afb....$")
    headfolderPatternFl  = re.compile("^fl....$")
    for idx in os.listdir(binfolder):
        if headfolderPatternAfb.match(idx) or headfolderPatternFl.match(idx):
            headfolders.append(idx)

    # Copy the wspaces and datacards, no overwritten
    wspacePattern = re.compile("^wspace_.*.root$")
    for idx in os.listdir(iwspacepath):
        if wspacePattern.match(idx):
            for headfolder in headfolders:
                if idx not in os.listdir(os.path.join(binfolder,headfolder)):
                    shutil.copy2(os.path.join(iwspacepath,idx),os.path.join(binfolder,headfolder))

    datacardPattern = re.compile("fitParameters"+str(workBin)+"\.txt$")
    for idx in os.listdir(idatacardpath):
        if datacardPattern.match(idx):
            for headfolder in headfolders:
                if idx not in os.listdir(os.path.join(binfolder,headfolder)):
                    shutil.copy2(os.path.join(idatacardpath,idx),os.path.join(binfolder,headfolder))

    # Start loop, collect all jobs in binfolder/batchJobs
    subfolderPattern = re.compile("^set[0-9]...$")
    if not os.path.exists(os.path.join(binfolder,"batchJobs")):
        os.makedirs(os.path.join(binfolder,"batchJobs"), 0755)
    for headfolder in headfolders:
        print "Processing headfolder {0}".format(headfolder)
        for idx in os.listdir(os.path.join(binfolder,headfolder)):
            if not subfolderPattern.match(idx):
                continue
            fullpath = os.path.join(binfolder,headfolder,idx)
            cmd = """{0} angular3D_bins "{2}/*.root" --iallpath {1} --iCombBkgWspacepath {2} --oallpath {2} --keepparam > {2}/runtime.log 2>&1""".format(fitFunc,os.path.join(binfolder,headfolder),fullpath)
            # cmd = """{0} wideQ "{2}/*.root" --iallpath {1} --iCombBkgWspacepath {2} --oallpath {2} --keepparam --keeplog > {2}/runtime.log 2>&1""".format(fitFunc,os.path.join(binfolder,headfolder),fullpath)
            batchScriptFile = open(os.path.join(binfolder,'batchJobs',idx+'_'+headfolder+'.sh'),'w')
            batchScript = re.sub('JOBDIR',fullpath,pbsBatchTemplate)
            batchScript = re.sub('JOBNAME',idx+'_'+headfolder,batchScript)
            batchScript = re.sub('JOBSTDOUT',os.path.join(fullpath,'stdout.log'),batchScript)
            batchScript = re.sub('JOBSTDERR',os.path.join(fullpath,'stderr.log'),batchScript)
            batchScript = re.sub('CMD',cmd,batchScript)
            batchScriptFile.write(batchScript)
            batchScriptFile.close()
    return 0

if __name__ == "__main__":
    sys.exit(main())

