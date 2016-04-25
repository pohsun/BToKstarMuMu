#!/usr/bin/env python

import os
import subprocess
import time
import sys
import re
import shutil
from threading import Thread

number_of_threads = 8

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
        time.sleep(5)
        cmd=["../fit", "angular3D_bins"  ,self.fullpath+"/*.root", "--iallpath="+self.headfolder, "--oallpath="+self.fullpath, "--iCombBkgWspacepath="+self.fullpath]
        # print cmd
        operation = subprocess.Popen(cmd)
        operation.wait()
        self.status=0

proclist = []

headfolder="../limit"
if len(sys.argv) == 2 :
    headfolder = sys.argv[1]
else:
    exit(1)

# if os.path.isfile(headfiler+"/FCScanResult.pdf"):
    # exit(0)

# Prepare wspaces
iwspacepath="../"
wspacePattern = re.compile("^wspace_.*.root$")
for idx in os.listdir(iwspacepath):
    if wspacePattern.match(idx):
        shutil.copy2(iwspacepath+idx,headfolder)

idatacardpath="../"
datacardPattern = re.compile("^fitParameters.\.txt$")
for idx in os.listdir(idatacardpath):
    if datacardPattern.match(idx):
        shutil.copy2(idatacardpath+idx,headfolder)

# Start loop
subfolderPattern = re.compile("^set[0-9]...$")
for idx in os.listdir(headfolder):
    if not subfolderPattern.match(idx):
        continue

    n_active_proc = number_of_threads
    while n_active_proc>=number_of_threads:
        n_active_proc = 0
        for proc in proclist:
            if proc.status<0: n_active_proc = n_active_proc + 1
        time.sleep(1)

    current = transfer_thread(headfolder,idx)
    proclist.append(current)
    current.start()
    time.sleep(1)

