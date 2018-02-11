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
#PBS -l walltime=01:30:00

source ~pchen/local/root/bin/thisroot.sh
cd /wk_cms/pchen/work/BuToKstarMuMu/devel/CMSSW_5_3_20/src/BphAna/BToKstarMuMu/plugins/macro
CMD

"""

def main():
    # parse all options
    parser = OptionParser()
    parser.add_option('-b', '--bin', dest='binfolder'           , help='The q2 bin to be processed', default=None, type='str')
    parser.add_option('-p', '--pat', dest='headfolderPattern'   , help='The pattern of headfolder to be processed', default=None, type='str')
    parser.add_option('-n', '--set', dest='nSets'               , help='The maximum number of toy sets', default=500 , type='int')
    (opt, args) = parser.parse_args()

    if opt.binfolder is None :
        print "ERROR   : bin must be given using '-b' option."
        return 1
    binfolder = opt.binfolder

    # Start loop, collect all jobs in binfolder/batchJobs
    headfolders = []
    if opt.headfolderPattern is None:
        headfolderPatterns = [ re.compile("^afb....$"), re.compile("^fl\+...$") ]
    elif opt.headfolderPattern != "":
        headfolderPatterns.append(re.compile(opt.headfolderPattern))
    else:
        headfolderPatterns = []
        headfolders.append("")

    for headfolder in os.listdir("../limit/{0}".format(binfolder)):
        if any( [ p.match(headfolder) for p in headfolderPatterns ] ):
            headfolders.append(headfolder)

    if not os.path.exists(os.path.join("../limit",binfolder,"batchJobs")):
        os.makedirs(os.path.join("../limit",binfolder,"batchJobs"), 0755)
    for headfolder in headfolders:
        print "Processing headfolder {0}".format(headfolder)
        fullpath = os.path.join("../limit",binfolder,headfolder)
        for iPart in range(0,opt.nSets/100+1):
            cmd = "./fitFCToys.py {0} \"{1}\" \"set{2:02d}..\"".format(binfolder,headfolder.replace("+","\+"),iPart)
            batchScriptFile = open(os.path.join("../limit",binfolder,'batchJobs',"{0}part{1:02d}.sh".format(headfolder+"_" if headfolder != "" else "",iPart+1)),'w')
            batchScript = re.sub('JOBDIR',os.path.abspath(fullpath),pbsBatchTemplate)
            batchScript = re.sub('JOBNAME',str(binfolder)+'_'+headfolder+str(iPart),batchScript)
            batchScript = re.sub('CMD',cmd,batchScript)
            batchScriptFile.write(batchScript)
            batchScriptFile.close()
    return 0

if __name__ == "__main__":
    sys.exit(main())

