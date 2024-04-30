#!/usr/bin/env python3
##
# Split all Chunks in the current dir in 1-file chunks.
# Note: current implementation supports only NanoAnalysis productions.
##

import os, shutil, glob
from ZZAnalysis.NanoAnalysis.tools import setConf, getConf


def splitChunk(chunk_) :
    chunk = chunk_.strip('/')
    with open(chunk+'/run_cfg.py', 'r') as fp:
        lines = fp.readlines()
        filesline = -1
        for line in lines:
            if line.find('setConf("fileNames"') != -1:
                exec(line)
                filesline = lines.index(line)
                break

        if filesline == -1 :
            raise ValueError(chunk+' is apparently not a nanoAOD job')
        
        fileNames = getConf('fileNames')
        if len(fileNames) > 1 :
            print('Splitting', chunk, ':')
            for ifile in range(len(fileNames)) :
                if ifile == 0 :
                    newchunk=chunk
                    shutil.copy2(chunk+'/run_cfg.py',newchunk+'/run_cfg.py.orig')
                else :
                    newchunk=chunk+'-'+str(ifile)
                    os.makedirs(newchunk+'/log', exist_ok=True)
                    shutil.copy2(chunk+'/batchScript.sh',newchunk)
                    
                print('  ', newchunk)
                pyfile = open(newchunk+"/run_cfg.py", "w")                
                for iline, line in enumerate(lines) :
                    if iline == filesline :
                        pyfile.write('setConf("fileNames",["'+fileNames[ifile]+'"])\n')
                    else: 
                        pyfile.write(line)
                pyfile.close()
        else :
            print(chunk, "is already a single-file job.")


chunks = glob.glob('*Chunk*/')
for chunk in chunks:
    splitChunk(chunk)
            
