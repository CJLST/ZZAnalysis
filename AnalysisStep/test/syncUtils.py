import math, os, sys#, ctypes
from ROOT import *
#from import c_int, c_float
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


class Candidate:

    def __init__(self, treeEntry, options) :

	isMC = False
        if treeEntry.GetBranch("genHEPMCweight") :
            isMC = True

        self.Z1Flav           = treeEntry.Z1Flav
        self.Z2Flav           = treeEntry.Z2Flav
        self.ZZFlav           = self.Z1Flav * self.Z2Flav

        self.Z1Mass           = treeEntry.Z1Mass
        self.Z2Mass           = treeEntry.Z2Mass
        self.massErrRaw    = treeEntry.ZZMassErr
        self.massErrCorr   = treeEntry.ZZMassErrCorr
        self.m4lRefit      = -1
        self.m4lRefitErr   = -1
        if options.longOutput:
            self.m4lRefit      = treeEntry.ZZMassRefit
            self.m4lRefitErr   = treeEntry.ZZMassRefitErr
 

        self.run   = treeEntry.RunNumber
        self.lumi  = treeEntry.LumiNumber
        self.event = treeEntry.EventNumber
        self.weight      = 1
        self.ZZMass      = treeEntry.ZZMass

        self.pt4l          = treeEntry.ZZPt
        self.nExtraLep     = treeEntry.nExtraLep
        self.nExtraZ       = treeEntry.nExtraZ
        self.jetpt         = treeEntry.JetPt
        self.jeteta        = treeEntry.JetEta
        self.jetphi        = treeEntry.JetPhi
        self.jetmass       = treeEntry.JetMass
        self.jetQGLikelihood = treeEntry.JetQGLikelihood
        self.njets30       = treeEntry.nCleanedJetsPt30
        self.njets30Btag   = treeEntry.nCleanedJetsPt30BTagged
        self.mjj           = treeEntry.DiJetMass
        self.detajj        = treeEntry.DiJetDEta
        self.pfMet         = treeEntry.PFMET
        self.weight        = 1.
        if (isMC) : self.weight = sign(treeEntry.genHEPMCweight) * treeEntry.PUWeight * treeEntry.dataMCWeight

        self.jets30pt = []
        self.jets30eta = []
        self.jets30phi = []
        self.jets30mass = []
        self.jets30QGLikelihood = []

        for i in range(len(treeEntry.JetPt)):
            if treeEntry.JetPt[i]>30.:
                self.jets30pt.append(treeEntry.JetPt[i])
                self.jets30eta.append(treeEntry.JetEta[i])
                self.jets30phi.append(treeEntry.JetPhi[i])
                self.jets30mass.append(treeEntry.JetMass[i])
                self.jets30QGLikelihood.append(treeEntry.JetQGLikelihood[i])



        self.kds         = None
        self.fishjj      = -1.
        self.isDiJet     = False
        self.jet1pt      = -1.
        self.jet2pt      = -1.
        self.jet1qgl     = -1.
        self.jet2qgl     = -1.
        self.fillJetInfo()

	if options.synchMode == 'HZZ' :
	    self.computeKDs(treeEntry)
		
#	    # ICHEP2016 categories
#	    self.category    = CDLL('libZZAnalysisAnalysisStep.so').categoryIchep16(
#	        c_int(self.nExtraLep),
#	        c_int(self.nExtraZ),
#	        c_int(self.njets30),
#	        c_int(self.njets30Btag),
#	        (c_float * len(self.jets30QGLikelihood))(*self.jets30QGLikelihood),
#	        c_float(treeEntry.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal),
#	        c_float(treeEntry.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal),
#	        c_float(treeEntry.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal),
#	        c_float(treeEntry.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal),
#	        c_float(treeEntry.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal),
#	        c_float(treeEntry.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal),
#	        c_float(treeEntry.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal),
#	        (c_float * len(self.jets30phi))(*self.jets30phi),
#	        c_float(self.ZZMass),
#	        c_bool(False)
#	        )
		
	    # Moriond2017 categories
	    self.category    = CDLL('libZZAnalysisAnalysisStep.so').categoryMor17(
	        c_int(self.nExtraLep),
	        c_int(self.nExtraZ),
	        c_int(self.njets30),
	        c_int(self.njets30Btag),
	        (c_float * len(self.jets30QGLikelihood))(*self.jets30QGLikelihood),
	        c_float(treeEntry.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal),
	        c_float(treeEntry.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal),
	        c_float(treeEntry.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal),
	        c_float(treeEntry.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal),
	        c_float(treeEntry.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal),
	        c_float(treeEntry.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal),
	        c_float(treeEntry.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal),
	        (c_float * len(self.jets30phi))(*self.jets30phi),
	        c_float(self.ZZMass),
	        c_float(self.pfMet),
	        c_bool(True), #useVHMETTagged
	        c_bool(False) #useQGTagging
	        )

    def computeKDs(self, treeEntry):

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


        lib = CDLL('libZZAnalysisAnalysisStep.so')
        lib.getDbkgkinConstant.restype = c_float
        lib.getDbkgConstant.restype = c_float
        lib.DVBF2j_ME.restype = c_float
        lib.DVBF1j_ME.restype = c_float
        lib.DWHh_ME.restype = c_float
        lib.DZHh_ME.restype = c_float
        lib.DVBF2j_ME_QG.restype = c_float
        lib.DVBF1j_ME_QG.restype = c_float
        lib.DWHh_ME_QG.restype = c_float
        lib.DZHh_ME_QG.restype = c_float

        self.D_bkg_kin  = treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + treeEntry.p_QQB_BKG_MCFM*lib.getDbkgkinConstant(c_int(int(self.ZZFlav)),c_float(self.ZZMass)))
        self.D_bkg      = treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen*treeEntry.p_m4l_SIG/(treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen*treeEntry.p_m4l_SIG+treeEntry.p_QQB_BKG_MCFM*treeEntry.p_m4l_BKG*lib.getDbkgConstant(c_int(int(self.ZZFlav)),c_float(self.ZZMass)))
        self.D_g4       = treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + pow(2.521, 2)*treeEntry.p_GG_SIG_ghg2_1_ghz4_1_JHUGen) # D_0-
        self.KD_highdim = treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + treeEntry.p_GG_SIG_ghg2_1_ghz2_1_JHUGen)
        self.KD_vec     = treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + treeEntry.p_QQB_SIG_ZPqqLR_1_gZPz2_1_JHUGen)
        self.KD_psvec   = treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + treeEntry.p_QQB_SIG_ZPqqLR_1_gZPz1_1_JHUGen)
        self.KD_gggrav  = treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + treeEntry.p_GG_SIG_gXg1_1_gXz1_1_gXz5_1_JHUGen)
        self.KD_qqgrav  = treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen/(treeEntry.p_GG_SIG_ghg2_1_ghz1_1_JHUGen + treeEntry.p_QQB_SIG_XqqLR_1_gXz1_1_gXz5_1_JHUGen)
        ##MELA-only production discriminants:
        if self.njets30 >= 2 :
            self.Djet_VAJHU = lib.DVBF2j_ME(
                c_float(treeEntry.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal),
                c_float(treeEntry.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal),
                c_float(self.ZZMass))
            self.D_WHh_VAJHU = lib.DWHh_ME(
                c_float(treeEntry.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal),
                c_float(treeEntry.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal),
                c_float(self.ZZMass))
            self.D_ZHh_VAJHU = lib.DZHh_ME(
                c_float(treeEntry.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal),
                c_float(treeEntry.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal),
                c_float(self.ZZMass))
        if self.njets30 == 1 :
            self.D_VBF1j_VAJHU = lib.DVBF1j_ME(
                c_float(treeEntry.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal),
                c_float(treeEntry.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal),
                c_float(treeEntry.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal),
                c_float(self.ZZMass))
        ##MELA+q/g production discriminants:
        if self.njets30 >= 2 :
            self.Dfull_VBF2j = lib.DVBF2j_ME_QG(
                c_float(treeEntry.p_JJVBF_SIG_ghv1_1_JHUGen_JECNominal),
                c_float(treeEntry.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal),
                c_float(self.ZZMass),
                (c_float * len(self.jets30QGLikelihood))(*self.jets30QGLikelihood),
                (c_float * len(self.jets30phi))(*self.jets30phi))
            self.Dfull_WHh = lib.DWHh_ME_QG(
                c_float(treeEntry.p_HadWH_SIG_ghw1_1_JHUGen_JECNominal),
                c_float(treeEntry.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal),
                c_float(self.ZZMass),
                (c_float * len(self.jets30QGLikelihood))(*self.jets30QGLikelihood),
                (c_float * len(self.jets30phi))(*self.jets30phi))
            self.Dfull_ZHh = lib.DZHh_ME_QG(
                c_float(treeEntry.p_HadZH_SIG_ghz1_1_JHUGen_JECNominal),
                c_float(treeEntry.p_JJQCD_SIG_ghg2_1_JHUGen_JECNominal),
                c_float(self.ZZMass),
                (c_float * len(self.jets30QGLikelihood))(*self.jets30QGLikelihood),
                (c_float * len(self.jets30phi))(*self.jets30phi))
        if self.njets30 == 1 :
            self.Dfull_VBF1j = lib.DVBF1j_ME_QG(
                c_float(treeEntry.p_JVBF_SIG_ghv1_1_JHUGen_JECNominal),
                c_float(treeEntry.pAux_JVBF_SIG_ghv1_1_JHUGen_JECNominal),
                c_float(treeEntry.p_JQCD_SIG_ghg2_1_JHUGen_JECNominal),
                c_float(self.ZZMass),
                (c_float * len(self.jets30QGLikelihood))(*self.jets30QGLikelihood),
                (c_float * len(self.jets30phi))(*self.jets30phi))


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


    def printOut(self, options):
        line = ""
	if options.synchMode == 'HZZ' :
            line  += str(int(self.run))
            line  += ":" + str(int(self.lumi))
            line  += ":" + str(int(self.event))
            line  += ":{0:.2f}".format(self.ZZMass)
            line  += ":" + "{0:.2f}".format(self.Z1Mass)
            line  += ":" + "{0:.2f}".format(self.Z2Mass)
#            line  += ":" + "{0:.2f}".format(self.massErrRaw)
#            line  += ":" + "{0:.2f}".format(self.massErrCorr)
            line  += ":" + "{0:.3f}".format(self.D_bkg_kin)
            line  += ":" + "{0:.3f}".format(self.D_bkg)
            line  += ":" + "{0:.3f}".format(self.Djet_VAJHU)
            line  += ":" + "{0:.3f}".format(self.D_g4)
            line  += ":" + "{0:.3f}".format(self.D_VBF1j_VAJHU)
            line  += ":" + "{0:.3f}".format(self.D_WHh_VAJHU)
            line  += ":" + "{0:.3f}".format(self.D_ZHh_VAJHU)
#            line  += ":" + "{0:.2f}".format(self.pt4l)
            line  += ":" + "{0:d}".format(self.njets30)
            line  += ":" + "{0:.2f}".format(self.jet1pt)
            line  += ":" + "{0:.2f}".format(self.jet2pt)
            line  += ":" + "{0:.3f}".format(self.jet1qgl)
            line  += ":" + "{0:.3f}".format(self.jet2qgl)
            line  += ":" + "{0:.3f}".format(self.Dfull_VBF2j)
            line  += ":" + "{0:.3f}".format(self.Dfull_VBF1j)
            line  += ":" + "{0:.3f}".format(self.Dfull_WHh)
            line  += ":" + "{0:.3f}".format(self.Dfull_ZHh)
#            line  += ":" + "{0:.2f}".format(self.mjj)
#            line  += ":" + "{0:.3f}".format(self.detajj)
#            line  += ":" + "{0:.3f}".format(self.fishjj)
#            line  += ":" + "{0:.3f}".format(self.kds.KD_highdim)
#            line  += ":" + "{0:.3f}".format(self.kds.KD_vec)
#            line  += ":" + "{0:.3f}".format(self.kds.KD_psvec)
#            line  += ":" + "{0:.3f}".format(self.kds.KD_gggrav)
#            line  += ":" + "{0:.3f}".format(self.kds.KD_qqgrav)
            line  += ":" + "{0:.3f}".format(self.pfMet)
            line  += ":" + "{0:d}".format(self.category)
            if self.m4lRefit>=0:
                line  += ":" + "{0:.2f}".format(self.m4lRefit)
                line  += ":" + "{0:.2f}".format(self.m4lRefitErr)
                line  += ":" + "{0:.3f}".format(self.weight)

	if options.synchMode == 'VBS' :
            line  += str(int(self.run))
            line  += ":" + str(int(self.lumi))
            line  += ":" + str(int(self.event))
	    channel = ''
	    if self.ZZFlav == 11**4 : channel = 'eeee'
	    if self.ZZFlav == 11**2*13**2 : channel = 'eemm'
	    if self.ZZFlav == 13**4 : channel = 'mmmm'
	    line  += ":%s"%channel
            line  += ":{0:.2f}".format(self.ZZMass)
            line  += ":" + "{0:.2f}".format(self.Z1Mass)
            line  += ":" + "{0:.2f}".format(self.Z2Mass)
            line  += ":" + "{0:d}".format(self.njets30)
            line  += ":" + "{0:.2f}".format(self.jet1pt)
            line  += ":" + "{0:.2f}".format(self.jet2pt)
            line  += ":" + "{0:.2f}".format(self.mjj)

        return line
