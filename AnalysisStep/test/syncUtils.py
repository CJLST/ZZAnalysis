import math, os, sys, ctypes
from ROOT import *
from ctypes import *

def sign(x):
    if x > 0:
        return 1.
    elif x < 0:
        return -1.
    elif x == 0:
        return 0.
    else:
        return x

def cubicroot(x):
    if x >= 0:
        return x**(1./3.)
    elif x < 0:
        return float('nan')

class Event:

    def __init__(self,run,lumi,event):

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

    def __init__(self,p_GG_SIG_ghg2_1_ghz1_1_JHUGen,p_GG_SIG_ghg2_1_ghz4_1_JHUGen,p_GG_SIG_ghg2_1_ghz2_1_JHUGen,p_QQB_SIG_ZPqqLR_1_gZPz2_1_JHUGen,p_QQB_SIG_ZPqqLR_1_gZPz1_1_JHUGen,p_GG_SIG_gXg1_1_gXz1_1_gXz5_1_JHUGen,p_QQB_SIG_XqqLR_1_gXz1_1_gXz5_1_JHUGen,p_QQB_BKG_MCFM,p_m4l_SIG,p_m4l_BKG,p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,p_HadWH_SIG_ghv1_1_JHUGen_JECNominal,p_HadZH_SIG_ghz1_1_JHUGen_JECNominal,njets30,jets30QGLikelihood,jets30phi,ZZFlav,ZZMass):

        self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen  = p_GG_SIG_ghg2_1_ghz1_1_JHUGen
        self.p_GG_SIG_ghg2_1_ghz4_1_JHUGen = p_GG_SIG_ghg2_1_ghz4_1_JHUGen
        self.p_GG_SIG_ghg2_1_ghz2_1_JHUGen = p_GG_SIG_ghg2_1_ghz2_1_JHUGen
        self.p_QQB_SIG_ZPqqLR_1_gZPz2_1_JHUGen  = p_QQB_SIG_ZPqqLR_1_gZPz2_1_JHUGen
        self.p_QQB_SIG_ZPqqLR_1_gZPz1_1_JHUGen      = p_QQB_SIG_ZPqqLR_1_gZPz1_1_JHUGen
        self.p_GG_SIG_gXg1_1_gXz1_1_gXz5_1_JHUGen      = p_GG_SIG_gXg1_1_gXz1_1_gXz5_1_JHUGen
        self.p_QQB_SIG_XqqLR_1_gXz1_1_gXz5_1_JHUGen   = p_QQB_SIG_XqqLR_1_gXz1_1_gXz5_1_JHUGen
        self.p_QQB_BKG_MCFM    = p_QQB_BKG_MCFM
        self.p_m4l_SIG    = p_m4l_SIG
        self.p_m4l_BKG       = p_m4l_BKG
        self.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal    = p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal
        self.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal    = p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal
        self.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal     = p_JQCD_SIG_ghg2_1_JHUGen_JECNominal
        self.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal     = p_JVBF_SIG_ghv1_1_JHUGen_JECNominal
        self.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal     = pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal
        self.p_HadWH_SIG_ghv1_1_JHUGen_JECNominal = p_HadWH_SIG_ghv1_1_JHUGen_JECNominal
        self.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal = p_HadZH_SIG_ghz1_1_JHUGen_JECNominal
        self.njets30       = njets30
        self.jets30QGL     = jets30QGLikelihood
        self.jets30phi     = jets30phi
        self.ZZFlav        = ZZFlav
        self.ZZMass        = ZZMass
        self.D_bkg_kin     = -1.
        self.D_g4          = -1.
        self.Dbkg          = -1.
        self.KD_pseudo     = -1.
        self.KD_highdim    = -1.
        self.KD_vec        = -1.
        self.KD_psvec      = -1.
        self.KD_gggrav     = -1.
        self.KD_qqgrav     = -1.
        self.Djet_VAJHU    = -1.
        self.D_WHh_VAJHU   = -1.
        self.D_ZHh_VAJHU   = -1.
        self.D_VBF1j_VAJHU = -1.
        self.Dfull_VBF2j   = -1.
        self.Dfull_WHh     = -1.
        self.Dfull_ZHh     = -1.
        self.Dfull_VBF1j   = -1.
        self.computeKDs()


    def computeKDs(self):

        lib = ctypes.CDLL('libZZAnalysisAnalysisStep.so')
        lib.getDbkgkinConstant.restype = ctypes.c_float
        lib.getDbkgConstant.restype = ctypes.c_float
        lib.getDVBF2jetsConstant.restype = ctypes.c_float
        lib.getDVBF1jetConstant.restype = ctypes.c_float

        self.D_bkg_kin  = self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + self.p_QQB_BKG_MCFM*lib.getDbkgkinConstant(c_int(int(self.ZZFlav)),c_float(self.ZZMass)))
        self.D_bkg      = self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen*self.p_m4l_SIG/(self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen*self.p_m4l_SIG+self.p_QQB_BKG_MCFM*self.p_m4l_BKG*lib.getDbkgConstant(c_int(int(self.ZZFlav)),c_float(self.ZZMass)))
        self.D_g4       = self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + pow(2.521, 2)*self.p_GG_SIG_ghg2_1_ghz4_1_JHUGen) # D_0-
        self.KD_highdim = self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + self.p_GG_SIG_ghg2_1_ghz2_1_JHUGen)
        self.KD_vec     = self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + self.p_QQB_SIG_ZPqqLR_1_gZPz2_1_JHUGen)
        self.KD_psvec   = self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + self.p_QQB_SIG_ZPqqLR_1_gZPz1_1_JHUGen)
        self.KD_gggrav  = self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + self.p_GG_SIG_gXg1_1_gXz1_1_gXz5_1_JHUGen)
        self.KD_qqgrav  = self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(self.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + self.p_QQB_SIG_XqqLR_1_gXz1_1_gXz5_1_JHUGen)
        ##MELA-only production discriminants:
        if self.njets30 >= 2 :
            self.Djet_VAJHU = self.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal/(self.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal+self.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal*lib.getDVBF2jetsConstant(c_float(self.ZZMass))) # VBF(2j) vs. gg->H+2j
            self.D_WHh_VAJHU   = self.p_HadWH_SIG_ghv1_1_JHUGen_JECNominal/(self.p_HadWH_SIG_ghv1_1_JHUGen_JECNominal+100000.*self.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal) # W(->2j)H vs. gg->H+2j
            self.D_ZHh_VAJHU   = self.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal/(self.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal+10000.*self.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal) # Z(->2j)H vs. gg->H+2j
        if self.njets30 == 1 :
            self.D_VBF1j_VAJHU = self.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*self.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal/(self.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal*self.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal+self.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal*lib.getDVBF1jetConstant(c_float(self.ZZMass))) # VBF(1j) vs. gg->H+1j
        ##MELA+q/g production discriminants:
        jets30PgOverPq = []
        for i in range(self.njets30):
            if self.jets30QGL[i]<0. and i<2 :
                rand = TRandom3()
                rand.SetSeed(abs(int(math.sin(self.jets30phi[i])*100000)))
                jets30PgOverPq.append( 1./rand.Uniform() - 1. )
            else:
                jets30PgOverPq.append( 1./self.jets30QGL[i] - 1. )
        if self.njets30 >= 2 :
            if self.jets30QGL[0] == 0. or self.jets30QGL[1] == 0.:
                self.Dfull_VBF2j = 0.
                self.Dfull_WHh = 0.
                self.Dfull_ZHh = 0.
            else:
                self.Dfull_VBF2j = 1/(1+ (1./self.Djet_VAJHU-1.) * cubicroot(jets30PgOverPq[0]*jets30PgOverPq[1]) ) ; # VBF(2j) vs. gg->H+2j
                self.Dfull_WHh = 1/(1+ (1./self.D_WHh_VAJHU-1.) * cubicroot(jets30PgOverPq[0]*jets30PgOverPq[1]) ) ; # W(->2j)H vs. gg->H+2j
                self.Dfull_ZHh = 1/(1+ (1./self.D_ZHh_VAJHU-1.) * cubicroot(jets30PgOverPq[0]*jets30PgOverPq[1]) ) ;  # VBF(2j) vs. gg->H+2j
        if self.njets30 == 1 :
            if self.jets30QGL[0] == 0.:
                self.Dfull_VBF1j = 0.
            else:
                self.Dfull_VBF1j = 1/(1+ (1./self.D_VBF1j_VAJHU-1.) * cubicroot(jets30PgOverPq[0]) ) ; # VBF(1j) vs. gg->H+1j


class Candidate:

    def __init__(self,event,m,mZ1,mZ2,mErr,mErrCorr,m4lRefit,m4lRefitErr,pt,nExtraLep,nExtraZ,jets30pt,jets30eta,jets30phi,jets30mass,njets30,njets30Btag,mjj,detajj,kds,weight,
                 jets30QGLikelihood,
                 p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal,p_JQCD_SIG_ghg2_1_JHUGen_JECNominal,p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal,p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,p_HadWH_SIG_ghv1_1_JHUGen_JECNominal,p_HadZH_SIG_ghz1_1_JHUGen_JECNominal
                 ):

        self.eventInfo   = event
        self.weight      = weight
        self.mass4l      = m
        self.mZ1         = mZ1
        self.mZ2         = mZ2
        self.massErrRaw  = mErr
        self.massErrCorr = mErrCorr
        self.m4lRefit    = m4lRefit
        self.m4lRefitErr = m4lRefitErr
        self.kds         = kds
        self.pt4l        = pt
        self.nExtraLep   = nExtraLep
        self.nExtraZ     = nExtraZ
        self.jets30pt    = jets30pt
        self.jets30eta   = jets30eta
        self.jets30phi   = jets30phi
        self.jets30mass  = jets30mass
        self.jets30QGLikelihood = jets30QGLikelihood
        self.njets30     = njets30
        self.njets30Btag = njets30Btag
        self.mjj         = mjj
        self.detajj      = detajj
        self.fishjj      = -1.
        self.isDiJet     = False
        self.jet1pt      = -1.
        self.jet2pt      = -1.
        self.jet1qgl     = -1.
        self.jet2qgl     = -1.
        self.fillJetInfo()

        # Winter 2015 version
#         self.category    = ctypes.CDLL('libZZAnalysisAnalysisStep.so').category(
#             c_int(nExtraLep),
#             c_float(self.pt4l),
#             c_float(self.mass4l),
#             c_int(self.njets30),
#             c_int(self.njets30Btag),
#             (ctypes.c_float * len(self.jets30pt  ))(*self.jets30pt  ),
#             (ctypes.c_float * len(self.jets30eta ))(*self.jets30eta ),
#             (ctypes.c_float * len(self.jets30phi ))(*self.jets30phi ),
#             (ctypes.c_float * len(self.jets30mass))(*self.jets30mass),
#             c_float(self.fishjj),
#             )
        # Summer 2016 version
        self.category    = ctypes.CDLL('libZZAnalysisAnalysisStep.so').categoryIchep16(
            c_int(nExtraLep),
            c_int(nExtraZ),
            c_int(njets30),
            c_int(njets30Btag),
            (ctypes.c_float * len(jets30QGLikelihood))(*jets30QGLikelihood),
            c_float(p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal),
            c_float(p_JQCD_SIG_ghg2_1_JHUGen_JECNominal),
            c_float(p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal),
            c_float(p_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
            c_float(pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal,
            c_float(p_HadWH_SIG_ghv1_1_JHUGen_JECNominal),
            c_float(p_HadZH_SIG_ghz1_1_JHUGen_JECNominal),
            (ctypes.c_float * len(jets30phi))(*jets30phi),
            c_float(self.mass4l),
            c_bool(False)
            )


    def fillJetInfo(self):

        if self.njets30==1:
            self.jet1pt = self.jets30pt[0]
            self.jet1qgl = self.jets30QGLikelihood[0]
            self.mjj = -1.
            self.detajj = -1.
        elif self.njets30>=2:
            self.jet1pt = self.jets30pt[0]
            self.jet2pt = self.jets30pt[1]
            self.jet1qgl = self.jets30QGLikelihood[0]
            self.jet2qgl = self.jets30QGLikelihood[1]
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
        line  += ":" + "{0:.3f}".format(self.kds.D_bkg_kin)
        line  += ":" + "{0:.3f}".format(self.kds.D_bkg)
        line  += ":" + "{0:.3f}".format(self.kds.Djet_VAJHU)
        line  += ":" + "{0:.3f}".format(self.kds.D_g4)
        line  += ":" + "{0:.3f}".format(self.kds.D_VBF1j_VAJHU)
        line  += ":" + "{0:.3f}".format(self.kds.D_WHh_VAJHU)
        line  += ":" + "{0:.3f}".format(self.kds.D_ZHh_VAJHU)
#        line  += ":" + "{0:.2f}".format(self.pt4l)
        line  += ":" + "{0:d}".format(self.njets30)
        line  += ":" + "{0:.2f}".format(self.jet1pt)
        line  += ":" + "{0:.2f}".format(self.jet2pt)
        line  += ":" + "{0:.3f}".format(self.jet1qgl)
        line  += ":" + "{0:.3f}".format(self.jet2qgl)
        line  += ":" + "{0:.3f}".format(self.kds.Dfull_VBF2j)
        line  += ":" + "{0:.3f}".format(self.kds.Dfull_VBF1j)
        line  += ":" + "{0:.3f}".format(self.kds.Dfull_WHh)
        line  += ":" + "{0:.3f}".format(self.kds.Dfull_ZHh)
#        line  += ":" + "{0:.2f}".format(self.mjj)
#        line  += ":" + "{0:.3f}".format(self.detajj)
#        line  += ":" + "{0:.3f}".format(self.fishjj)
#        line  += ":" + "{0:.3f}".format(self.kds.KD_highdim)
#        line  += ":" + "{0:.3f}".format(self.kds.KD_vec)
#        line  += ":" + "{0:.3f}".format(self.kds.KD_psvec)
#        line  += ":" + "{0:.3f}".format(self.kds.KD_gggrav)
#        line  += ":" + "{0:.3f}".format(self.kds.KD_qqgrav)
        line  += ":" + "{0:d}".format(self.category)
        if self.m4lRefit>=0:
            line  += ":" + "{0:.2f}".format(self.m4lRefit)
            line  += ":" + "{0:.2f}".format(self.m4lRefitErr)
            line  += ":" + "{0:.3f}".format(self.weight)

        return line
