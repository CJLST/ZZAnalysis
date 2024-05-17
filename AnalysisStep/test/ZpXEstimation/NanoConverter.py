## python3 NanoConverter.py
import ROOT
import os,sys
import optparse
from ZZAnalysis.NanoAnalysis.tools import get_genEventSumw

usage = ('usage: %prog [options]\n'
             + '%prog -h for help')
parser = optparse.OptionParser(usage)
parser.add_option('',   '--mc', action='store_true', dest='MC', default=False, help='MC samples')
parser.add_option('',   '--input',dest='INPUT',type='string',default='', help='Name and path to the input file')
parser.add_option('',   '--output',dest='OUTPUT',type='string',default='', help='Name and path to the output file')
parser.add_option('',   '--skipZL', action='store_true', dest='SKIPZL', default=False, help='Skip the ZL tree (e.g., necessary for ZZto4l samples)')
(opt, args) = parser.parse_args()

def makeCR(_df, _flag):

    # CRZLLss 21 --> 2097152
    # CRZLLos_2P2F 22 --> 4194304
    # CRZLLos_3P1F 23 --> 8388608
    if _flag == '3P1F':
        bit = '8388608'
    elif _flag == '2P2F':
        bit = '4194304'
    elif _flag == 'SS':
        bit = '2097152'
    elif _flag == 'SIPCR':
        bit = '21'
    else:
        raise Exception("The CR "+_flag+" is not known")

    df_out = ( _df.Filter('ZLLbest'+_flag+'Idx>-1').Define("ZZMass", "ZLLCand_mass[ZLLbest"+_flag+"Idx]")
                                                   .Define("RunNumber", "run")
                                                   .Define("EventNumber", "event")
                                                   .Define("LumiNumber", "luminosityBlock")
                                                   .Define("CRflag", bit)
                                                   .Define("Z1Mass", "ZLLCand_Z1mass[ZLLbest"+_flag+"Idx]")
                                                   .Define("Z2Mass", "ZLLCand_Z2mass[ZLLbest"+_flag+"Idx]")
                                                   .Define("Z1Flav", "ZLLCand_Z1flav[ZLLbest"+_flag+"Idx]")
                                                   .Define("Z2Flav", "ZLLCand_Z2flav[ZLLbest"+_flag+"Idx]")
                                                   .Define('Leptons_pt', "concatenate(Electron_pt,Muon_pt)")
                                                   .Define('Leptons_eta', "concatenate(Electron_eta,Muon_eta)")
                                                   .Define('Leptons_phi', "concatenate(Electron_phi,Muon_phi)")
                                                   .Define('Leptons_dxy', "concatenate(Electron_dxy,Muon_dxy)")
                                                   .Define('Leptons_dz', "concatenate(Electron_dz,Muon_dz)")
                                                   .Define('Leptons_id', "concatenate(Electron_pdgId,Muon_pdgId)")
                                                   .Define('Leptons_sip', "concatenate(Electron_sip3d,Muon_sip3d)")
                                                   .Define('Leptons_iso', "concatenate(Electron_pfRelIso03FsrCorr,Muon_pfRelIso03FsrCorr)")
                                                   .Define('Leptons_isid', "concatenate(Electron_passBDT,Muon_passID)")
                                                   ## Need to add the LepMissingHit branch for SS FR method
                                                   ## First create a dummy branch for muons filled with zeroes
                                                   .Define('Muon_lostHits', "addDummyBranch(Muon_pt)")
                                                   .Define('Leptons_missinghit', "concatenate(Electron_lostHits, Muon_lostHits)")
                                                   ## Variable miniAOD-style
                                                   .Define('LepPt', "std::vector<float> LepPt{Leptons_pt[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_pt[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_pt[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_pt[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return LepPt")
                                                   .Define('LepEta', "std::vector<float> LepEta{Leptons_eta[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_eta[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_eta[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_eta[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return LepEta")
                                                   .Define('LepPhi', "std::vector<float> LepPhi{Leptons_phi[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_phi[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_phi[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_phi[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return LepPhi")
                                                   .Define('Lepdxy', "std::vector<float> Lepdxy{Leptons_dxy[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_dxy[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_dxy[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_dxy[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return Lepdxy")
                                                   .Define('Lepdz', "std::vector<float> Lepdz{Leptons_dz[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_dz[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_dz[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_dz[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return Lepdz")
                                                   .Define('LepLepId', "std::vector<short> LepLepId{Leptons_id[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_id[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_id[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_id[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return LepLepId")
                                                   .Define('LepSIP', "std::vector<float> LepSIP{Leptons_sip[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_sip[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_sip[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_sip[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return LepSIP")
                                                   .Define('LepCombRelIsoPF', "std::vector<float> LepCombRelIsoPF{Leptons_iso[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_iso[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_iso[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_iso[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return LepCombRelIsoPF")
                                                   .Define('LepisID', "std::vector<bool> LepisID{Leptons_isid[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_isid[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_isid[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_isid[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return Leptons_isid")
                                                   .Define('LepMissingHit', "std::vector<unsigned char> LepMissingHit{Leptons_missinghit[ZLLCand_Z1l1Idx[ZLLbest"+_flag+"Idx]], Leptons_missinghit[ZLLCand_Z1l2Idx[ZLLbest"+_flag+"Idx]], Leptons_missinghit[ZLLCand_Z2l1Idx[ZLLbest"+_flag+"Idx]], Leptons_missinghit[ZLLCand_Z2l2Idx[ZLLbest"+_flag+"Idx]]}; return LepMissingHit")
                                                   .Define('PFMET', "MET_pt")
                                                   ## overallEventWeight contains everything in NanoAODs
                                                   .Define('L1prefiringWeight', "1") ## Dummy
                                                   .Define('KFactor_EW_qqZZ', "1") ## Dummy
                                                   .Define('KFactor_QCD_qqZZ_M', "1") ## Dummy
                                                   .Define('KFactor_QCD_ggZZ_Nominal', '1') ## Dummy
                                                   .Define('xsec', '1') ## Dummy
                                                   )
    return df_out

ROOT.gInterpreter.Declare("""
ROOT::RVec<float> concatenate(ROOT::RVec<float> &A, ROOT::RVec<float> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<short> concatenate(ROOT::RVec<int> &A, ROOT::RVec<int> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<bool> concatenate(ROOT::RVec<bool> &A, ROOT::RVec<bool> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")

ROOT.gInterpreter.Declare("""
ROOT::RVec<unsigned char> concatenate(ROOT::RVec<unsigned char> &A, ROOT::RVec<unsigned char> &B){
    int sizeA = A.size();
    int sizeB = B.size();
    int sizeC = sizeA + sizeB;
    ROOT::RVec<float> C(sizeC);

    for (int i = 0; i < sizeA; ++i) {
        C[i] = A[i];
    }
    for (int i = 0; i < sizeB; ++i) {
        C[sizeA+i] = B[i];
    }

    return C;
}
""")


ROOT.gInterpreter.Declare("""
ROOT::RVec<unsigned char> addDummyBranch(ROOT::RVec<float> &A){
    int sizeA = A.size();
    ROOT::RVec<unsigned char> B(sizeA);

    for (int i = 0; i < sizeA; ++i) {
        B[i] = 0;
    }

    return B;
}
""")



##################################### MAIN #####################################
MC = opt.MC

inFileName = opt.INPUT
outFileName = opt.OUTPUT

df = ROOT.RDataFrame('Events', inFileName)

df_3P1F = makeCR(df, "3P1F")
df_2P2F = makeCR(df, "2P2F")
df_2P2Lss = makeCR(df, "SS")
df_SIP = makeCR(df, "SIPCR")

opts = ROOT.RDF.RSnapshotOptions()
opts.fMode = 'RECREATE'
## Variables to store in the output root file
vars = {'RunNumber',
        'EventNumber',
        'LumiNumber',
        'ZZMass',
        'CRflag',
        'Z1Flav',
        'Z2Flav',
        'Z1Mass',
        'Z2Mass',
        'LepisID',
        'PFMET',
        'LepPt',
        'LepEta',
        'LepPhi',
        'Lepdxy',
        'Lepdz',
        'LepLepId',
        'LepSIP',
        'LepCombRelIsoPF',
        'LepMissingHit',
        }
if MC:
    vars.add('overallEventWeight')
    vars.add('xsec')
    vars.add('L1prefiringWeight')
    vars.add('KFactor_EW_qqZZ')
    vars.add('KFactor_QCD_qqZZ_M')
    vars.add('KFactor_QCD_ggZZ_Nominal')


df_3P1F.Snapshot('CRZLLTree/candTree', "test_3P1F.root", vars, opts)
df_2P2F.Snapshot('CRZLLTree/candTree', "test_2P2F.root", vars, opts)
df_2P2Lss.Snapshot('CRZLLTree/candTree', "test_2P2Lss.root", vars, opts)
df_SIP.Snapshot('CRZLLTree/candTree', "test_SIP.root", vars, opts)

## RooDataFrames cannot be concatenated.
## Solution: save each CR in a different tree and then merge them through another RooDataFrame
## Not fancy, but it works
df_bis = ROOT.RDataFrame('CRZLLTree/candTree', "test_*.root")
skip_ZLL = True
if df_bis.Count().GetValue() != 0:
    ## If the number of entries in df_bis is zero (when debugging on few events)
    ## you get an error when trying to take the snapshot
    df_bis.Snapshot('CRZLLTree/candTree', outFileName, vars, opts)
    skip_ZLL = False
## Remove the intermediate files, we don't need them anymore
os.system('rm test_*.root')

## SR
df_SR = ( df.Filter('bestCandIdx>=0').Define("ZZMass", "ZZCand_mass[bestCandIdx]") ## Dummy
                                     .Define("CRflag", "0") ## Dummy
                                     .Define("Z1Flav", "ZZCand_Z1flav[bestCandIdx]")
                                     .Define("Z2Flav", "ZZCand_Z2flav[bestCandIdx]") ## Dummy
                                     .Define("Z1Mass", "ZZCand_Z1mass[bestCandIdx]")
                                     .Define("Z2Mass", "ZZCand_Z2mass[bestCandIdx]") ## Dummy
                                     .Define("RunNumber", "run")
                                     .Define("EventNumber", "event")
                                     .Define("LumiNumber", "luminosityBlock")
                                     .Define('Leptons_pt', "concatenate(Electron_pt,Muon_pt)")
                                     .Define('Leptons_eta', "concatenate(Electron_eta,Muon_eta)")
                                     .Define('Leptons_phi', "concatenate(Electron_phi,Muon_phi)")
                                     .Define('Leptons_dxy', "concatenate(Electron_dxy,Muon_dxy)")
                                     .Define('Leptons_dz', "concatenate(Electron_dz,Muon_dz)")
                                     .Define('Leptons_id', "concatenate(Electron_pdgId,Muon_pdgId)")
                                     .Define('Leptons_sip', "concatenate(Electron_sip3d,Muon_sip3d)")
                                     .Define('Leptons_iso', "concatenate(Electron_pfRelIso03FsrCorr,Muon_pfRelIso03FsrCorr)")
                                     .Define('Leptons_isid', "concatenate(Electron_passBDT,Muon_passID)")
                                     ## Need to add the LepMissingHit branch for SS FR method
                                     ## First create a dummy branch for muons filled with zeroes
                                     .Define('Muon_lostHits', "addDummyBranch(Muon_pt)")
                                     .Define('Leptons_missinghit', "concatenate(Electron_lostHits, Muon_lostHits)")
                                     ## Variable miniAOD-style
                                     .Define('LepPt', "std::vector<float> LepPt{Leptons_pt[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_pt[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_pt[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_pt[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepPt")
                                     .Define('LepEta', "std::vector<float> LepEta{Leptons_eta[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_eta[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_eta[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_eta[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepEta")
                                     .Define('LepPhi', "std::vector<float> LepPhi{Leptons_phi[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_phi[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_phi[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_phi[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepPhi")
                                     .Define('Lepdxy', "std::vector<float> Lepdxy{Leptons_dxy[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_dxy[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_dxy[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_dxy[ZZCand_Z2l2Idx[bestCandIdx]]}; return Lepdxy")
                                     .Define('Lepdz', "std::vector<float> Lepdz{Leptons_dz[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_dz[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_dz[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_dz[ZZCand_Z2l2Idx[bestCandIdx]]}; return Lepdz")
                                     .Define('LepLepId', "std::vector<short> LepLepId{Leptons_id[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_id[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_id[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_id[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepLepId")
                                     .Define('LepSIP', "std::vector<float> LepSIP{Leptons_sip[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_sip[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_sip[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_sip[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepSIP")
                                     .Define('LepCombRelIsoPF', "std::vector<float> LepCombRelIsoPF{Leptons_iso[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_iso[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_iso[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_iso[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepCombRelIsoPF")
                                     .Define('LepisID', "std::vector<bool> LepisID{Leptons_isid[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_isid[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_isid[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_isid[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepisID")
                                     .Define('LepMissingHit', "std::vector<unsigned char> LepMissingHit{Leptons_missinghit[ZZCand_Z1l1Idx[bestCandIdx]], Leptons_missinghit[ZZCand_Z1l2Idx[bestCandIdx]], Leptons_missinghit[ZZCand_Z2l1Idx[bestCandIdx]], Leptons_missinghit[ZZCand_Z2l2Idx[bestCandIdx]]}; return LepMissingHit")
                                     .Define('PFMET', "MET_pt")
                                     ## overallEventWeight contains everything in NanoAODs
                                     .Define('L1prefiringWeight', "1") ## Dummy
                                     .Define('KFactor_EW_qqZZ', "1") ## Dummy
                                     .Define('KFactor_QCD_qqZZ_M', "1") ## Dummy
                                     .Define('KFactor_QCD_ggZZ_Nominal', '1') ## Dummy
                                     .Define('xsec', '1') ## Dummy
                                     )

opts.fMode = 'UPDATE'
df_SR.Snapshot('ZZTree/candTree', outFileName, vars, opts)



## ZL CR for the computation of fake rates
if not opt.SKIPZL:
    df_ZL = ( df.Filter('ZLCand_lepIdx>-1').Define("ZZMass", "0") ## Dummy
                                           .Define("CRflag", "0") ## Dummy
                                           .Define("Z1Flav", "ZCand_flav[bestZIdx]")
                                           .Define("Z2Flav", "0") ## Dummy
                                           .Define("Z1Mass", "ZCand_mass[bestZIdx]")
                                           .Define("Z2Mass", "0") ## Dummy
                                           .Define("RunNumber", "run")
                                           .Define("EventNumber", "event")
                                           .Define("LumiNumber", "luminosityBlock")
                                           .Define('Leptons_pt', "concatenate(Electron_pt,Muon_pt)")
                                           .Define('Leptons_eta', "concatenate(Electron_eta,Muon_eta)")
                                           .Define('Leptons_phi', "concatenate(Electron_phi,Muon_phi)")
                                           .Define('Leptons_dxy', "concatenate(Electron_dxy,Muon_dxy)")
                                           .Define('Leptons_dz', "concatenate(Electron_dz,Muon_dz)")
                                           .Define('Leptons_id', "concatenate(Electron_pdgId,Muon_pdgId)")
                                           .Define('Leptons_sip', "concatenate(Electron_sip3d,Muon_sip3d)")
                                           .Define('Leptons_iso', "concatenate(Electron_pfRelIso03FsrCorr,Muon_pfRelIso03FsrCorr)")
                                           .Define('Leptons_isid', "concatenate(Electron_passBDT,Muon_passID)")
                                           ## Need to add the LepMissingHit branch for SS FR method
                                           ## First create a dummy branch for muons filled with zeroes
                                           .Define('Muon_lostHits', "addDummyBranch(Muon_pt)")
                                           .Define('Leptons_missinghit', "concatenate(Electron_lostHits, Muon_lostHits)")
                                           ## Variable miniAOD-style
                                           .Define('LepPt', "std::vector<float> LepPt{Leptons_pt[ZCand_l1Idx[bestZIdx]], Leptons_pt[ZCand_l2Idx[bestZIdx]], Leptons_pt[ZLCand_lepIdx]}; return LepPt")
                                           .Define('LepEta', "std::vector<float> LepEta{Leptons_eta[ZCand_l1Idx[bestZIdx]], Leptons_eta[ZCand_l2Idx[bestZIdx]], Leptons_eta[ZLCand_lepIdx]}; return LepEta")
                                           .Define('LepPhi', "std::vector<float> LepPhi{Leptons_phi[ZCand_l1Idx[bestZIdx]], Leptons_phi[ZCand_l2Idx[bestZIdx]], Leptons_phi[ZLCand_lepIdx]}; return LepPhi")
                                           .Define('Lepdxy', "std::vector<float> Lepdxy{Leptons_dxy[ZCand_l1Idx[bestZIdx]], Leptons_dxy[ZCand_l2Idx[bestZIdx]], Leptons_dxy[ZLCand_lepIdx]}; return Lepdxy")
                                           .Define('Lepdz', "std::vector<float> Lepdz{Leptons_dz[ZCand_l1Idx[bestZIdx]], Leptons_dz[ZCand_l2Idx[bestZIdx]], Leptons_dz[ZLCand_lepIdx]}; return Lepdz")
                                           .Define('LepLepId', "std::vector<short> LepLepId{Leptons_id[ZCand_l1Idx[bestZIdx]], Leptons_id[ZCand_l2Idx[bestZIdx]], Leptons_id[ZLCand_lepIdx]}; return LepLepId")
                                           .Define('LepSIP', "std::vector<float> LepSIP{Leptons_sip[ZCand_l1Idx[bestZIdx]], Leptons_sip[ZCand_l2Idx[bestZIdx]], Leptons_sip[ZLCand_lepIdx]}; return LepSIP")
                                           .Define('LepCombRelIsoPF', "std::vector<float> LepCombRelIsoPF{Leptons_iso[ZCand_l1Idx[bestZIdx]], Leptons_iso[ZCand_l2Idx[bestZIdx]], Leptons_iso[ZLCand_lepIdx]}; return LepCombRelIsoPF")
                                           .Define('LepisID', "std::vector<bool> LepisID{Leptons_isid[ZCand_l1Idx[bestZIdx]], Leptons_isid[ZCand_l2Idx[bestZIdx]], Leptons_isid[ZLCand_lepIdx]}; return LepisID")
                                           .Define('LepMissingHit', "std::vector<unsigned char> LepMissingHit{Leptons_missinghit[ZCand_l1Idx[bestZIdx]], Leptons_missinghit[ZCand_l2Idx[bestZIdx]], Leptons_missinghit[ZLCand_lepIdx]}; return LepMissingHit")
                                           .Define('PFMET', "MET_pt")
                                           ## overallEventWeight contains everything in NanoAODs
                                           .Define('L1prefiringWeight', "1") ## Dummy
                                           .Define('KFactor_EW_qqZZ', "1") ## Dummy
                                           .Define('KFactor_QCD_qqZZ_M', "1") ## Dummy
                                           .Define('KFactor_QCD_ggZZ_Nominal', '1') ## Dummy
                                           .Define('xsec', '1') ## Dummy
                                           )

    opts.fMode = 'UPDATE'
    df_ZL.Snapshot('CRZLTree/candTree', outFileName, vars, opts)

## Add counter only with the 40th entry
counters = ROOT.TH1F("Counters", "Counters", 50, 0, 100)
if opt.MC:
    root = ROOT.TFile.Open(inFileName)
    genEventSumw = get_genEventSumw(root)
    counters.SetBinContent(40, genEventSumw)
    root.Close()

root_file = ROOT.TFile(outFileName, "UPDATE")
if not skip_ZLL:
    sub_dir = root_file.Get("CRZLLTree")
    sub_dir.cd()
    counters.Write()

if not opt.SKIPZL:
    sub_dir = root_file.Get("CRZLTree")
    sub_dir.cd()
    counters.Write()
root_file.Close()
