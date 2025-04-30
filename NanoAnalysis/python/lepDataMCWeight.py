'''Add data/MC weights to leptons.
'''

from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from ROOT import LeptonSFHelper


class lepDataMCWeight(Module):
    def __init__(self, year, data_tag):
        '''Add data/MC weights to leptons.'''
        
        print("***lepDataMCWeight: year:", year, "data_tag:", data_tag, flush=True)
        self.year = year
        self.lepSFHelper = LeptonSFHelper(data_tag)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("Muon_dataMC", "F", lenVar="nMuon", title="data/MC correction", limitedPrecision=12)
        self.out.branch("Electron_dataMC", "F", lenVar="nElectron", title="data/MC correction", limitedPrecision=12)

    def analyze(self, event):
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")

        e_SFs = [1.]*event.nElectron
        for ie, ele in enumerate(electrons):
            e_SFs[ie] = self.getLepSF(ele)

        m_SFs = [1.]*event.nMuon
        for im, mu in enumerate(muons):
            m_SFs[im] = self.getLepSF(mu)
        
        self.out.fillBranch("Electron_dataMC", e_SFs)    
        self.out.fillBranch("Muon_dataMC", m_SFs)    

        return True

    
    def getLepSF(self, lep):
        '''Return lepton efficiency scale factor'''

        if self.year > 2023 : #FIXME: 2023/2024 SFs not yet implemented!
            return 1.
        myLepID = abs(lep.pdgId)
        mySCeta = lep.eta
        isCrack = False # FIXME: isGap() is not available in nanoAODs, and cannot be recomputed easily based on eta, phi. We thus use the non-gap SFs for all electrons.
        if myLepID==11 :
            mySCeta = lep.eta + lep.deltaEtaSC # Use the SC eta and not the electron eta

        # Deal with very rare cases when SCeta is out of 2.5 bounds
        mySCeta = min(mySCeta,2.49)
        mySCeta = max(mySCeta,-2.49)

        SF = self.lepSFHelper.getSF(self.year, myLepID, lep.pt, lep.eta, mySCeta, isCrack)
        # SF_Unc = self.lepSFHelper.getSFError(year, myLepID, lep.pt, lep.eta, mySCeta, isCrack)
        return SF
      
