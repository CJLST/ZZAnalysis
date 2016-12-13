#!/bin/env python

import sys
import imp
import copy
import os
import shutil
import pickle
import math
import pprint
import subprocess
from optparse import OptionParser


class TPLtoPY:
   def __init__(self):
      # define options and arguments ====================================
      self.parser = OptionParser()
      self.parser.add_option("--outdir", dest="outdir", type="string",
                                 help="Name of the output directory",
                                 default=None)

      self.parser.add_option("--outname", dest="outname", type="string",
                                 help="Name of the output fragment",
                                 default=None)

      self.parser.add_option("--template", dest="tplFragment", type="string",
                                help="Template fragment file",
                                default=None)

      self.parser.add_option("--indicators", dest="indicators", type="string",
                                 help="Comma-separated indicators to pick. Format: tplkey:replacevalue,tplkey:replacevalue",
                                 default=None)


      (self.opt,self.args) = self.parser.parse_args()
      self.mkdir(self.opt.outdir)
      self.tplarray = self.split(self.opt.indicators)
      self.outfile = "{}/{}".format(self.opt.outdir,self.opt.outname)
      self.replaceStrings()

   def replaceStrings(self):
      os.system("cp {} {}".format(self.opt.tplFragment,self.outfile))
      for pair in self.tplarray:
         os.system("sed -i \'s/{}/{}/g\' {}".format(pair[0],pair[1],self.outfile))

   def split(self,string):
      if string is None:
         sys.exit("TPLtoPY::split: No string provided")
      theList=[]
      x = string.split(',')
      for y in x:
         yy = y.split(':')
         theList.append(yy)
      return theList

   def mkdir(self,dirname):
      if( os.system('mkdir -p {}'.format(dirname)) != 0 ):
         sys.exit("TPLtoPY::mkdir: Please remove or rename directory: {}".format(dirname))


if __name__ == '__main__':
   theClass = TPLtoPY()
