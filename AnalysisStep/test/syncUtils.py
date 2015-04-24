import os, sys, ctypes
from ctypes import *

class Event:

    def __init__(self,iBC,run,lumi,event):

        self.iBC   = iBC
        self.run   = run
        self.lumi  = lumi
        self.event = event

    def printOut(self):

        line   = ""
        line  += str(int(self.run))
        line  += ":" + str(int(self.lumi))
        line  += ":" + str(int(self.event))
        
        return line

class KDs:

    def __init__(self,p0plus_VAJHU,p0minus_VAJHU,p0hplus_VAJHU,p1plus_VAJHU,p1_VAJHU,p2_VAJHU,p2qqb_VAJHU,bkg_VAMCFM):

        self.p0plus_VAJHU  = p0plus_VAJHU
        self.p0minus_VAJHU = p0minus_VAJHU
        self.p0hplus_VAJHU = p0hplus_VAJHU
        self.p1plus_VAJHU  = p1plus_VAJHU 
        self.p1_VAJHU      = p1_VAJHU     
        self.p2_VAJHU      = p2_VAJHU     
        self.p2qqb_VAJHU   = p2qqb_VAJHU              
        self.bkg_VAMCFM    = bkg_VAMCFM
        self.KD            = -1.
        self.KD_pseudo     = -1.
        self.KD_highdim    = -1.
        self.KD_vec        = -1.
        self.KD_psvec      = -1.
        self.KD_gggrav     = -1.
        self.KD_qqgrav     = -1.
        self.computeKDs()


    def computeKDs(self):

        self.KD         = self.p0plus_VAJHU/(self.p0plus_VAJHU + self.bkg_VAMCFM) 	 
        self.KD_pseudo  = self.p0plus_VAJHU/(self.p0plus_VAJHU + self.p0minus_VAJHU)
        self.KD_highdim = self.p0plus_VAJHU/(self.p0plus_VAJHU + self.p0hplus_VAJHU)
        self.KD_vec     = self.p0plus_VAJHU/(self.p0plus_VAJHU + self.p1plus_VAJHU) 	 
        self.KD_psvec   = self.p0plus_VAJHU/(self.p0plus_VAJHU + self.p1_VAJHU) 	 
        self.KD_gggrav  = self.p0plus_VAJHU/(self.p0plus_VAJHU + self.p2_VAJHU) 	 
        self.KD_qqgrav  = self.p0plus_VAJHU/(self.p0plus_VAJHU + self.p2qqb_VAJHU) 


class Candidate:

    def __init__(self,event,m,mZ1,mZ2,mErr,mErrCorr,pt,nExtraLep,jets30pt,jets30eta,njets30Btag,mjj,detajj,kds):

        self.eventInfo   = event
        self.mass4l      = m
        self.mZ1         = mZ1
        self.mZ2         = mZ2
        self.massErrRaw  = mErr
        self.massErrCorr = mErrCorr
        self.kds         = kds
        self.pt4l        = pt
        self.nExtraLep   = nExtraLep
        self.jets30pt    = jets30pt
        self.jets30eta   = jets30eta
        self.njets30     = len(jets30pt)
        self.njets30Btag = njets30Btag
        self.mjj         = mjj
        self.detajj      = detajj
        self.fishjj      = -1.
        self.isDiJet     = False
        self.jet1pt      = -1.
        self.jet2pt      = -1.
        self.jet1eta     = -99.
        self.jet2eta     = -99.
        self.fillJetInfo()
        self.category    = ctypes.CDLL('libZZAnalysisAnalysisStep.so').category(
            c_int(nExtraLep),
            c_float(self.pt4l),
            c_float(self.mass4l),
            c_int(self.njets30),
            c_int(self.njets30Btag),
            c_float(self.jet1pt),
            c_float(self.jet2pt),
            c_float(self.jet1eta),
            c_float(self.jet2eta),
            c_float(self.mjj),
            c_float(self.fishjj),
            )

    def fillJetInfo(self):
        
        if self.njets30==1:
            self.jet1pt = self.jets30pt[0]
            self.jet1eta = self.jets30eta[0]
            self.mjj = -1.
            self.detajj = -1.
        elif self.njets30>=2:
            self.jet1pt = self.jets30pt[0]            
            self.jet2pt = self.jets30pt[1]
            self.jet1eta = self.jets30eta[0]            
            self.jet2eta = self.jets30eta[1]                    
            self.fishjj = 0.18*abs(self.detajj) + 1.92e-04*self.mjj
        else:
            self.mjj = -1.
            self.detajj = -1.
            

    def printOut(self):

        line = ""
        line  += "{0:.2f}".format(self.mass4l)
        line  += ":" + "{0:.2f}".format(self.mZ1)
        line  += ":" + "{0:.2f}".format(self.mZ2)
#        line  += ":" + "{0:.2f}".format(self.massErrRaw)
#        line  += ":" + "{0:.2f}".format(self.massErrCorr)
        line  += ":" + "{0:.3f}".format(self.kds.KD)
#        line  += ":" + "{0:.2f}".format(self.pt4l)
        line  += ":" + "{0:d}".format(self.njets30)
        line  += ":" + "{0:.2f}".format(self.jet1pt)
        line  += ":" + "{0:.2f}".format(self.jet2pt)
#        line  += ":" + "{0:.2f}".format(self.mjj)
#        line  += ":" + "{0:.3f}".format(self.detajj)
#        line  += ":" + "{0:.3f}".format(self.fishjj)
#        line  += ":" + "{0:.3f}".format(self.kds.KD_pseudo)
#        line  += ":" + "{0:.3f}".format(self.kds.KD_highdim)
#        line  += ":" + "{0:.3f}".format(self.kds.KD_vec)
#        line  += ":" + "{0:.3f}".format(self.kds.KD_psvec)
#        line  += ":" + "{0:.3f}".format(self.kds.KD_gggrav)
#        line  += ":" + "{0:.3f}".format(self.kds.KD_qqgrav) 
#        line  += ":" + "{0:d}".format(self.category)       

        return line
