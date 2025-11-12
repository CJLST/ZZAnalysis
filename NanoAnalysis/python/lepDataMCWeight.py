'''Add data/MC weights to leptons.
'''

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ROOT import LeptonSFHelper


class lepDataMCWeight(Module):
    def __init__(self, year, data_tag, muonIdByMVA = False):
        '''Add data/MC weights to leptons.'''
        
        print("***lepDataMCWeight: year:", year, "data_tag:", data_tag, "muonIdByMVA:", muonIdByMVA, flush=True)
        self.year = year
        self.lepSFHelper = LeptonSFHelper(year, data_tag+("MUON_ID_BYMVA" if muonIdByMVA else "")) # Squeeze muonIdByMVA in data_tag to avoid changing the C++ interface

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Muon_dataMC", "F", lenVar="nMuon", title="data/MC correction", limitedPrecision=12)
        self.out.branch("Muon_dataMCUnc", "F", lenVar="nMuon", title="data/MC correction relative uncertainty", limitedPrecision=12)
        self.out.branch("Electron_dataMC", "F", lenVar="nElectron", title="data/MC correction", limitedPrecision=12)
        self.out.branch("Electron_dataMCUnc", "F", lenVar="nElectron", title="data/MC correction relative uncertainty", limitedPrecision=12)

    def analyze(self, event):
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")

        e_SFs = [1.]*event.nElectron
        e_SFsUnc = [1.]*event.nElectron
        for ie, ele in enumerate(electrons):
            e_SFs[ie], e_SFsUnc[ie] = self.getLepSF(ele)

        m_SFs = [1.]*event.nMuon
        m_SFsUnc = [1.]*event.nMuon
        for im, mu in enumerate(muons):
            m_SFs[im], m_SFsUnc[im] = self.getLepSF(mu)
        
        self.out.fillBranch("Electron_dataMC", e_SFs)    
        self.out.fillBranch("Electron_dataMCUnc", e_SFsUnc)
        self.out.fillBranch("Muon_dataMC", m_SFs)
        self.out.fillBranch("Muon_dataMCUnc", m_SFsUnc)

        return True

    
    def getLepSF(self, lep):
        '''Return lepton efficiency scale factor'''

        myLepID = abs(lep.pdgId)
        mySCeta = lep.eta
        isCrack = False # FIXME: isGap() is not available in nanoAODs, and cannot be recomputed easily based on eta, phi. We thus use the non-gap SFs for all electrons.
        isHoleBPix = False  # default

        if myLepID==11 :
            mySCeta = lep.eta + lep.deltaEtaSC # Use the SC eta and not the electron eta

        # Deal with very rare cases when SCeta is out of 2.5 bounds
        mySCeta = min(mySCeta,2.49)
        mySCeta = max(mySCeta,-2.49)

        pair = self.lepSFHelper.getSF(myLepID, lep.pt, lep.eta, mySCeta, lep.phi, isCrack)
        SF = pair.first
        SFerror = pair.second

        # Add a protection for leptons outside standard acceptance (pt<5/7 for mu/ele or |eta|>2.4 for mu) which get SF=0 and SFError = nan, since they may be still
        # be used for dedicated studies
        if SF==0 :
            SF, SFerror = 1., 0.5
        return SF, SFerror
      
