#!/usr/bin/env python

'''
Check the list of samples in a .py file like analyzer_2011.py
and write a black list of the samples with broken path:

chmod u+x checkSample.py
./checkSample.py -f analyzer_2011.py

N.B.: Require CMGTools/Production to be installed!

'''

import subprocess
import os
from CMGTools.Production.eostools import *
from optparse import OptionParser


cmgToolBase = '/store/cmst3/user/cmgtools/CMG'
cmgGroupBase = '/store/cmst3/group/cmgtools/CMG'

def main():
    parser = OptionParser()
    parser.add_option("-f", "--file", type="string", dest="fileName", metavar="SAMPLEFILE", default="analyzer_2012.py",
                      help="name of the file name")
    (options, args) = parser.parse_args()

    if os.path.exists(options.fileName) == False:
        print 'Sorry, file ',options.fileName,' does not exist'
        exit(1)
    else:
        print 'Parsing sample names from ',options.fileName
        if os.path.exists('blackList.txt') == True:
            print 'A file called blackList.txt already exists, remove it first'
            exit(1)
        else:
            outList = open('blackList.txt','w')
            parseListOfSamples(options.fileName,outList)
            outList.close()



def parseListOfSamples(file, blackList):
    inputFile = open(file,'r')
    for aLine in inputFile:
        if 'cmgTuple.*root' in aLine:

            '''Get relevant informations from parsing the line'''
            splitALine = aLine.split("'")
            sampleName = splitALine[1]
            cmgTag = splitALine[3]
            samplePath = splitALine[5]

            '''Skip commented samples (since they are not relevant)'''
            if '#' in splitALine[0]:
                continue

            '''Translate info in a useful way'''
            if cmgTag == 'cmgtools':
                samplePath = cmgToolBase+samplePath
            elif cmgTag == 'cmgtools_group':
                samplePath = cmgGroupBase+samplePath
            else:
                print 'Uncorrect CMG directory tag for sample ',sampleName,': ',cmgTag
                exit(1)

            ''' Check samples and write down the black list of broken one'''
            print  'Checking sample ',sampleName,' in path: ',samplePath

            lsOut = runEOSCommand(samplePath, 'ls')[0]
            if 'cmgTuple' in lsOut:
                print '---> OK'
                continue
            else:
                print 'No cmgTuple in path: ',samplePath,' for sample ',sampleName
                blackList.write('Check path: '+samplePath+' for sample '+sampleName)
                

            



if __name__ == '__main__':
    main()
