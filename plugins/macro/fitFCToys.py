#!/usr/bin/env python

import os
import subprocess
import time
import sys
import re
import shutil
from threading import Thread

number_of_threads = 11

class transfer_thread(Thread):
    def __init__ (self,headfolder,subfolder):
        Thread.__init__(self)
        self.headfolder = headfolder
        self.subfolder  = subfolder
        self.fullpath   = headfolder+"/"+subfolder
        self.status=-1

    def run(self):
        print time.ctime()
        print "Processing folder",self.headfolder,self.subfolder
        sys.stdout.flush()
        cmd=["../fit", "angular3D_bins"  ,self.fullpath+"/*.root", "--iallpath="+self.headfolder, "--oallpath="+self.fullpath, "--iCombBkgWspacepath="+self.fullpath, "--keeplog", "--keepparam"]
        # cmd=["echo", "Dry run in "+self.fullpath]
        # if "UnderProcessing" in os.listdir(self.fullpath):
            # cmd=["echo", "Skip UnderProcessing tag in "+self.fullpath]
        # print cmd
        time.sleep(3)
        operation = subprocess.Popen(cmd)
        operation.wait()
        self.status=0

proclist = []

binfolder="../limit/bin"
if len(sys.argv) == 2 :
    theBin = sys.argv[1]
    binfolder = binfolder+sys.argv[1]
    print "Running with "+binfolder
else:
    exit(1)


headfolders = []
headfolderPattern = re.compile("^afb...._fl....$")
for idx in os.listdir(binfolder):
    if headfolderPattern.match(idx):
        headfolders.append(idx)

# Prepare wspaces
iwspacepath="../wspace_ANv11_v3p5/"
wspacePattern = re.compile("^wspace_.*.root$")
for idx in os.listdir(iwspacepath):
    if wspacePattern.match(idx):
        for headfolder in headfolders:
            if idx not in os.listdir(binfolder+'/'+headfolder):
                shutil.copy2(iwspacepath+idx,binfolder+'/'+headfolder)

idatacardpath="../datacard_ANv11_v3p5/"
datacardPattern = re.compile("^fitParameters"+str(theBin)+"\.txt$")
for idx in os.listdir(idatacardpath):
    if datacardPattern.match(idx):
        for headfolder in headfolders:
            if idx not in os.listdir(binfolder+'/'+headfolder):
                shutil.copy2(idatacardpath+idx,binfolder+'/'+headfolder)

# Start loop
subfolderPattern = re.compile("^set[0-9]...$")
fitResultPattern = "wspace_angular3D_bin"+str(theBin)+".root"
for headfolder in headfolders:
    print "Processing headfolder {0}".format(headfolder)
    for idx in os.listdir(binfolder+'/'+headfolder):
        fullpath = binfolder+'/'+headfolder+'/'+idx
        if not subfolderPattern.match(idx):
            continue

        if fitResultPattern in os.listdir(fullpath):
            # print "Old fitting result "+fullpath+'/'+fitResultPattern+" is found. Skip!"
            if "UnderProcessing" in os.listdir(fullpath):
                os.remove(fullpath+"/UnderProcessing")

        elif "UnderProcessing" in os.listdir(fullpath):
            # Clean the tags during development
            # if fitResultPattern not in os.listdir(fullpath):
                # os.remove(fullpath+"/UnderProcessing")

            print "Tag for running job is found in "+fullpath+". Skip!"
        else:
            print "Processing subfolder {0}".format(idx)
            file=open(fullpath+"/UnderProcessing",'w')
            file.close()

            n_active_proc = number_of_threads
            while n_active_proc>=number_of_threads:
                n_active_proc = 0
                for proc in proclist:
                    if proc.status<0: n_active_proc = n_active_proc + 1
                time.sleep(1)

            current = transfer_thread(binfolder+'/'+headfolder,idx)
            proclist.append(current)
            current.start()
            time.sleep(1)

