from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.HeppyCore.utils.deltar import deltaR


class LHEFiller(Module):
    """Categorize LHEParticles for the purposes of calculating MELA probabilities. 
    """
    def __init__(self):
        print("***LHEFiller", flush=True)
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("LHEPart_MELAStatus", "S", lenVar="nLHEPart", title="Classification for LHEMothers (MELASTATUS=1), LHEDaughters (MELASTATUS=2), LHEAssociatedParticles (MELASTATUS=3) formerly used in the miniAOD format.")
    def analyze(self, event):
        LHEPart = Collection(event, 'LHEPart')
        MELA_Status = [0.]*len(LHEPart)

        for i, lp in enumerate(LHEPart):
            if lp.status == -1:
                # intermediate or initial particle
                MELA_Status[i] = 1

            elif lp.status == 1:
                parentIdx = lp.firstMotherIdx

                # guard against invalid parent index
                if parentIdx < 0 :
                    parentPdg = 0
                    gParentPdg = 0
                    MELA_Status[i] = -1  # or 0 if you prefer a neutral code
                else:
                    parentPdg = LHEPart[parentIdx].pdgId
                    gParentIdx = LHEPart[parentIdx].firstMotherIdx

                    # guard against invalid grandparent index
                    if gParentIdx < 0 :
                        gParentPdg = 0
                    else:
                        gParentPdg = LHEPart[gParentIdx].pdgId

                    if parentPdg == 23 and gParentPdg == 25:
                        MELA_Status[i] = 2  # Z from Higgs
                    else:
                        MELA_Status[i] = 3  # other final-state particle

            else:
                MELA_Status[i] = -1

        self.out.fillBranch("LHEPart_MELAStatus", MELA_Status)
        return True
