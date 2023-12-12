#!/usr/bin/env python3
# 
#
from __future__ import print_function

import os
import glob
import time
import re
import sys
from optparse import OptionParser

if __name__ == '__main__':

    parser = OptionParser()
    parser.usage = """
    %prog <dir>
    Scan chunk logs in <dir>, and print job CPU statistics for each sample
    """

    (options,args) = parser.parse_args()
    if len(args)>=1:
        path = args[0]
    else :
        path="AAAOK/Chunks"

    print("Scanning logs in", path)
    files = glob.glob(path+"/*Chunk*/log/*.log")

    if len(files)==0:
        print ("No log found.")
        sys.exit(1)
    
    datasets = {}
    for file in files:
        chunkname = os.path.basename(os.path.dirname((os.path.dirname(file))))
        dataset = re.sub("_Chunk.*$","",chunkname)  
        with open(file, 'r') as f:
            for line in f.readlines():
                if "Total Remote Usage" in line:
                    cpustring=line.split()[2]
                    cpustring=cpustring.replace(",","")
                    tt=time.strptime(cpustring,"%H:%M:%S")
                    cpu=tt.tm_hour*3600 + tt.tm_min*60 + tt.tm_sec
#                    print(chunkname, dataset,  line.split()[2], cpu)
                    if dataset in datasets:
                        datasets[dataset] = [min(cpu,datasets[dataset][0]),cpu+datasets[dataset][1],max(cpu,datasets[dataset][2]),datasets[dataset][3]+1]
                    else:
                        datasets[dataset] = [cpu, cpu, cpu, 1]
    
    print("CPU usage: njobs, min/avg/max (s):")
    for dataset in datasets:
        print(dataset,  datasets[dataset][3], "{:.1f}/{:.1f}/{:.1f}".format(datasets[dataset][0],datasets[dataset][1]/datasets[dataset][3],datasets[dataset][2]))
