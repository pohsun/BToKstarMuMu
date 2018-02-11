#!/usr/bin/env python

# Usage: ./fitFCToys.py binfolder [headfolderPattern] [subfolderPattern]

import os
import subprocess
import time
import sys
import re
import shutil
from threading import Thread

number_of_threads = 11

# Prepare wspaces
iwspacepath="../wspace_ANv17/"
idatacardpath="../datacard_ANv17/"

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
        ### This is for nominal bins
        cmd=["../fit", "angular3D_bins"  ,self.fullpath+"/*.root",
                "--iallpath="+self.headfolder,
                "--oallpath="+self.fullpath,
                # "--iCombBkgWspacepath="+self.fullpath,
                "--keeplog", "--keepparam",
                ">/dev/null" , "2>&1"]
        # cmd=["echo", "Dry run in "+self.fullpath]
        if "pbsRuntime" in os.listdir(self.fullpath):
            cmd=["echo", "Skip pbsRuntime tag in "+self.fullpath]
        print cmd
        os.system("touch "+self.fullpath+"/pbsRuntime")
        time.sleep(3)
        operation = subprocess.Popen(cmd)
        operation.wait()
        self.status=0

proclist = []

if len(sys.argv) >= 1 :
    if sys.argv[1] not in os.listdir("../limit/"):
        exit(1)
    binfolder = "../limit/"+sys.argv[1]
    print "Running with "+binfolder
else:
    exit(1)

headfolders = []
headfolderPatterns = []
if len(sys.argv) > 2:
    print "Appending headfolderPattern: {0}".format(sys.argv[2])
    if len(sys.argv[2]) != 0:
        headfolderPatterns.append(re.compile(sys.argv[2]))
    else:
        headfolders.append("")
else:
    headfolderPatterns.append(re.compile("^afb....$"))
    headfolderPatterns.append(re.compile("^fl\+...$") )

if len(sys.argv) > 3:
    print "Appending subfolderPattern: {0}".format(sys.argv[3])
    subfolderPattern = re.compile(sys.argv[3])
else:
    subfolderPattern = re.compile("^set[0-9]...$")

for pat in headfolderPatterns:
    for idx in os.listdir(binfolder):
        if pat.match(idx):
            headfolders.append(idx)

wspacePattern = re.compile("^wspace_.*.root$")
for idx in os.listdir(iwspacepath):
    if wspacePattern.match(idx):
        for headfolder in headfolders:
            if idx not in os.listdir(binfolder+'/'+headfolder):
                shutil.copy2(iwspacepath+idx,binfolder+'/'+headfolder)

datacardPattern = re.compile("^fitParameters.*\.txt$")
for idx in os.listdir(idatacardpath):
    if datacardPattern.match(idx):
        for headfolder in headfolders:
            if idx not in os.listdir(binfolder+'/'+headfolder):
                shutil.copy2(idatacardpath+idx,binfolder+'/'+headfolder)

# Start loop
fitResultPattern = "wspace_angular3D_bin.*\.root"
for headfolder in headfolders:
    print "Processing headfolder {0}".format(headfolder)
    for idx in os.listdir(binfolder+'/'+headfolder):
        if not subfolderPattern.match(idx):
            continue
        fullpath = binfolder+'/'+headfolder+'/'+idx
        if fitResultPattern in os.listdir(fullpath):
            print "Old fitting status "+fullpath+'/'+fitResultPattern+" is found. Try skip!"
            if "pbsRuntime" in os.listdir(fullpath):
                os.remove(os.path.join(fullpath,"pbsRuntime"))
            continue
        elif "pbsRuntime" in os.listdir(fullpath):
            print "Tag for running job is found in "+fullpath+". Skip!"
            continue

        print "Processing subfolder {0}".format(idx)
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

