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

      self.parser.add_option("--indicators", dest="indicators", type="string", action="append",
                                 help="Indicators to pick; can specify multiple indicators. Format: tplkey:replacevalue",
                                 default=None)


      (self.opt,self.args) = self.parser.parse_args()
      self.mkdir(self.opt.outdir)
      self.outfile = "{}/{}".format(self.opt.outdir,self.opt.outname)
      self.replaceStrings()

   def replaceStrings(self):
      os.system("cp {} {}".format(self.opt.tplFragment,self.outfile))
      for indicator in self.opt.indicators:
         pair = indicator.split(':')
         if len(pair)!=2:
            sys.exit("Pair size is not 2, original indicator was {}".format(indicator))
         #print "Replacing {} with {}".format(pair[0],pair[1])
         os.system("sed -i \'s/{}/{}/g\' {}".format(pair[0],pair[1],self.outfile))

   def mkdir(self,dirname):
      if( os.system('mkdir -p {}'.format(dirname)) != 0 ):
         sys.exit("TPLtoPY::mkdir: Please remove or rename directory: {}".format(dirname))


if __name__ == '__main__':
   theClass = TPLtoPY()
