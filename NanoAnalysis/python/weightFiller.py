##
# Compute per-event MC weights.
##
from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
import os
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
import h5py
import scipy.interpolate as interp

class weightFiller(Module):
    def __init__(self, XS, APPLY_K_NNLOQCD_ZZGG, APPLY_K_NNLOQCD_NLOEW_ZZQQB, APPLY_QCD_GGF_UNCERT, LEPTON_SETUP):
        print(
            "***weightFiller: XS:", XS, 
            "APPLY_K_NNLOQCD_ZZGG:", APPLY_K_NNLOQCD_ZZGG,
            "APPLY_K_NNLOQCD_NLOEW_ZZQQB:", APPLY_K_NNLOQCD_NLOEW_ZZQQB, 
            "APPLY_QCD_GGF_UNCERT:", APPLY_QCD_GGF_UNCERT, 
            "LEPTON_SETUP:", LEPTON_SETUP
            flush=True
        )
        self.writeHistFile = False
        self.XS = XS
        self.APPLY_K_NNLOQCD_ZZGG = APPLY_K_NNLOQCD_ZZGG
        self.APPLY_K_NNLOQCD_NLOEW_ZZQQB = APPLY_K_NNLOQCD_NLOEW_ZZQQB
        self.APPLY_QCD_GGF_UNCERT = APPLY_QCD_GGF_UNCERT
        self.LEPTON_SETUP = LEPTON_SETUP

        # Just use the lepton setup as an analog for the run number
        if LEPTON_SETUP in (2016, 2017, 2018):
            run_number = 2
        elif LEPTON_SETUP in (2022, 2023, 2024, 2025):
            run_number = 3

        basePath_nano=f'{os.environ['CMSSW_BASE']}/src/ZZAnalysis/NanoAnalysis/data/kFactors'

        if self.APPLY_K_NNLOQCD_ZZGG in (1,2):
            ## ggZZ QCD k-factors
            self.spkfactor_ggzz = [None]*4

            #Apart from the nominal case, the other 3 are simply up and down factors
            strZZGGKFVar = ["nominal", "PDF", "PDF_aS", "QCD_mu"]
            if self.APPLY_K_NNLOQCD_ZZGG == 1:
                raw_dat = h5py.File(f"{basePath_nano}/gluonFusion/Run{run_number}-NNLO.h5")
            elif self.APPLY_K_NNLOQCD_ZZGG == 2:
                raw_dat = h5py.File(f"{basePath_nano}/gluonFusion/Run{run_number}-NLO-NNLO.h5")

            for i, var in enumerate(strZZGGKFVar):
                self.spkfactor_ggzz[i] = interp.make_interp_spline(
                    raw_dat['hmass'], raw_dat[var], k=1
                )
            raw_dat.close()
            del raw_dat
        elif self.APPLY_K_NNLOQCD_ZZGG != 0:
            print ("Unsupported: APPLY_K_NNLOQCD_ZZGG=", self.APPLY_K_NNLOQCD_ZZGG)
            print ("Supported values are: 0 (do not apply), 1 (LO to NNLO), 2 (NLO to NNLO)")
            exit(1)

        if self.APPLY_K_NNLOQCD_NLOEW_ZZQQB :
            strQQZZKFVar = ["nominal", "QCD_up", "QCD_dn", "EW_factor", "smoothing_factor"]
            self.spkfactor_qqzz = [[None]*len(strQQZZKFVar)]*4
            for i in range(4):
                #cos(theta^*) is symmetric around 0, but
                #the actual files go from -1 to 0, so just reverse them for aesthetic reasons
                with h5py.File(f"{basePath_nano}/qqBarToZZ/Run{run_number}/kfac_m4l__cos{3-i}.h5") as raw_dat:
                    for j, var in enumerate(strQQZZKFVar):
                        self.spkfactor_qqzz[i][j] = interp.make_interp_spline(
                            raw_dat['zzmass'], raw_dat[var], k=1
                        )
                del raw_dat

        # ggH NNLOPS weights
        if self.APPLY_QCD_GGF_UNCERT :
            basePath = f'{os.environ['CMSSW_BASE']}/src/ZZAnalysis/AnalysisStep/'
            if self.LEPTON_SETUP >= 2022:
                NNLOPS_weight_file = ROOT.TFile.Open(basePath+'data/ggH_NNLOPS_Weights/NNLOPS_reweight_13p6.root')
            else:
                NNLOPS_weight_file = ROOT.TFile.Open(basePath+'data/ggH_NNLOPS_Weights/NNLOPS_reweight.root')
            self.gr_NNLOPSratio_pt_powheg_0jet = NNLOPS_weight_file.Get("gr_NNLOPSratio_pt_powheg_0jet")
            self.gr_NNLOPSratio_pt_powheg_1jet = NNLOPS_weight_file.Get("gr_NNLOPSratio_pt_powheg_1jet")
            self.gr_NNLOPSratio_pt_powheg_2jet = NNLOPS_weight_file.Get("gr_NNLOPSratio_pt_powheg_2jet")
            self.gr_NNLOPSratio_pt_powheg_3jet = NNLOPS_weight_file.Get("gr_NNLOPSratio_pt_powheg_3jet")


    def evalSpline(self, cStar, zzMass, j):
        cStar = abs(cStar)
        if cStar < 0.25:
            i = 0
        elif cStar < 0.5:
            i = 1
        elif cStar < 0.75:
            i = 2
        else:
            i = 3
        return self.spkfactor_qqzz[i][j](zzMass)


    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.APPLY_K_NNLOQCD_ZZGG > 0 :
            self.out.branch("KFactor_QCD_ggZZ_Nominal_Weight", "F", title="QCD k-factor for ggZZ")
            self.out.branch("KFactor_QCD_ggZZ_PDF_Factor", "F", title="Multiplicative factor for PDF variation")
            self.out.branch("KFactor_QCD_ggZZ_QCD_Factor", "F", title="Multiplicative factor for QCD variation")
            self.out.branch("KFactor_QCD_ggZZ_aS_Factor", "F", title="Multiplicative factor for aS variation")
        if self.APPLY_K_NNLOQCD_NLOEW_ZZQQB :
            self.out.branch("KFactor_qqZZ_Nominal_Weight", "F", title="Combined EW/QCD k-factor for qqZZ, with proper treatment of EW=1 below 2mZ")
            self.out.branch("KFactor_qqZZ_QCDup_Factor", "F", title="Multiplicative factor for QCD_up variation")
            self.out.branch("KFactor_qqZZ_QCDdn_Factor", "F", title="Multiplicative factor for QCD_down variation")
            self.out.branch("KFactor_qqZZ_EW_Factor", "F", title="Multiplicative factor for EW factorization variation")
            self.out.branch("KFactor_qqZZ_smooth_Factor", "F", title="Multiplicative factor for EW smoothing variation")
        if self.APPLY_QCD_GGF_UNCERT :
            self.out.branch("ggH_NNLOPS_Weight", "F", title="Reweighting for ggH as a function of njets and pT")

            
        self.out.branch("overallEventWeight", "F", title="Event weight: Generator_weight*XS*puWeight*(relevant k-factors where applicable). Must be normalized by sum of genEventSumw in the Runs tree")

    def analyze(self, event):
        KFactor_ZZQQB_Nominal = 1.
        KFactor_ZZQQB_QCD_up = 1.
        KFactor_ZZQQB_QCD_dn = 1.
        KFactor_ZZQQB_EW_factor = 1.
        KFactor_ZZQQB_smooth_factor = 1.
    ############ GLUON FUSION KFACTOR VALUES ##########
        KFactor_QCD_ggZZ_Nominal = 1.
        KFactor_QCD_ggZZ_PDF = 1.
        KFactor_QCD_ggZZ_aS = 1.
        KFactor_QCD_ggZZ_QCD = 1.

        ggH_NNLOPS_Weight = 1.

        ### QCD weights for ggH, ggZZ
        if self.APPLY_K_NNLOQCD_ZZGG :
            KFactor_QCD_ggZZ_Nominal = self.spkfactor_ggzz[0](event.GenZZ_mass)

            KFactor_QCD_ggZZ_PDF = self.spkfactor_ggzz[1](event.GenZZ_mass)
            KFactor_QCD_ggZZ_aS = self.spkfactor_ggzz[2](event.GenZZ_mass)
            KFactor_QCD_ggZZ_QCD = self.spkfactor_ggzz[3](event.GenZZ_mass)

        if self.APPLY_K_NNLOQCD_NLOEW_ZZQQB:
            KFactor_ZZQQB_Nominal = self.evalSpline(
                event.LHEMela_costhetastar,event.GenZZ_mass, 0
            )
            KFactor_ZZQQB_QCD_up = self.evalSpline(
                event.LHEMela_costhetastar,event.GenZZ_mass, 1
            )
            KFactor_ZZQQB_QCD_dn = self.evalSpline(
                event.LHEMela_costhetastar,event.GenZZ_mass, 2
            )
            KFactor_ZZQQB_EW_factor = self.evalSpline(
                event.LHEMela_costhetastar,event.GenZZ_mass, 3
            )
            KFactor_ZZQQB_smooth_factor = self.evalSpline(
                event.LHEMela_costhetastar,event.GenZZ_mass, 4
            )


        if self.APPLY_QCD_GGF_UNCERT :
            htxsNJets = event.HTXS_njets30
            htxsHPt = event.HTXS_Higgs_pt
            if htxsNJets==0 :
                ggH_NNLOPS_Weight = self.gr_NNLOPSratio_pt_powheg_0jet.Eval(min(htxsHPt,125.0))
            elif htxsNJets==1 :
                ggH_NNLOPS_Weight = self.gr_NNLOPSratio_pt_powheg_1jet.Eval(min(htxsHPt,625.0))
            elif htxsNJets==2 :
                ggH_NNLOPS_Weight = self.gr_NNLOPSratio_pt_powheg_2jet.Eval(min(htxsHPt,800.0))
            elif htxsNJets>=3 :
                ggH_NNLOPS_Weight = self.gr_NNLOPSratio_pt_powheg_3jet.Eval(min(htxsHPt,925.0))


        # L1 pre-firing weights. FIXME: not available in Nano02Apr2020
        # L1prefiringWeight = event.L1PreFiringWeight_Nom
        # L1prefiringWeightUp = event.L1PreFiringWeight_Up
        # L1prefiringWeightDn = event.L1PreFiringWeight_Dn
        
        #FIXME: event.ZZ_dataMCWeight is not included, since that can be stored per-candidate if storeAllCands=True.
        w_total = self.XS * event.Generator_weight * event.puWeight * KFactor_ZZQQB_Nominal * KFactor_QCD_ggZZ_Nominal * ggH_NNLOPS_Weight

        if self.APPLY_K_NNLOQCD_ZZGG : 
            self.out.fillBranch("KFactor_QCD_ggZZ_Nominal_Weight", KFactor_QCD_ggZZ_Nominal)
            #For systematic variations, multiply w_total by either (1+factor) or (1-factor)
            self.out.fillBranch("KFactor_QCD_ggZZ_PDF_Factor", KFactor_QCD_ggZZ_PDF)
            self.out.fillBranch("KFactor_QCD_ggZZ_aS_Factor", KFactor_QCD_ggZZ_aS)
            self.out.fillBranch("KFactor_QCD_ggZZ_QCD_Factor", KFactor_QCD_ggZZ_QCD)
        if self.APPLY_K_NNLOQCD_NLOEW_ZZQQB :
            self.out.fillBranch("KFactor_qqZZ_Nominal_Weight", KFactor_ZZQQB_Nominal)
            #for QCD variations multiply w_total by factor up or down
            self.out.fillBranch("KFactor_qqZZ_QCDup_Factor", KFactor_ZZQQB_QCD_up)
            self.out.fillBranch("KFactor_qqZZ_QCDdn_Factor", KFactor_ZZQQB_QCD_dn)
            #for these variations multiply w_total by (1+factor) or (1-factor)
            self.out.fillBranch("KFactor_qqZZ_EW_Factor", KFactor_ZZQQB_EW_factor)
            self.out.fillBranch("KFactor_qqZZ_smooth_Factor", KFactor_ZZQQB_smooth_factor)
        if self.APPLY_QCD_GGF_UNCERT :
            self.out.fillBranch("ggH_NNLOPS_Weight", ggH_NNLOPS_Weight)

        self.out.fillBranch("overallEventWeight", w_total)

        return True


