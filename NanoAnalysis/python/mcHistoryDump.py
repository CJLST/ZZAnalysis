"""Print Gen history and LHE particles on output.
See also: https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/common/hepmcDump.py
"""

from __future__ import print_function
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from prettytable import PrettyTable

class mcHistoryDump(Module):
    def __init__(self, printGen=True, printLHE=False, nanoversion=12) :
        """Print Gen history and LHE particles.
        """ 
        self.printGen=printGen
        self.printLHE=printLHE
        self.nanoversion = nanoversion

    def analyze(self, event):
        if self.printGen :
            print(f"\n--------mcHistoryDump for event {event.run}:{event.luminosityBlock}:{event.event}" )
            print("---GenPart---") 
            genpart=Collection(event,"GenPart")
            self.gen_logger(genpart)
        if self.printLHE :
            LHE=Collection(event,"LHEPart")
            self.lhe_logger(LHE)
        return True

    def format_value(self, val):
        if isinstance(val, float):
            return f"{val:.3f}"
        return str(val)

    def gen_logger(self, genpart):
        '''Print genpart history on output
        '''
        table = PrettyTable(['i', 'pdgId', 'mother', 'pT', 'eta', 'phi', 'mass', 'status'])
        for i, gp in enumerate(genpart) :
            # motherId=-1
            # gmotherId=-1
            # if gp.genPartIdxMother >= 0 :
            #     motherId = genpart[gp.genPartIdxMother].pdgId
            #     if genpart[gp.genPartIdxMother].genPartIdxMother >= 0 :
            #         gmotherId = genpart[genpart[gp.genPartIdxMother].genPartIdxMother].pdgId
            table.add_row([self.format_value(val) for val in [i, gp.pdgId, gp.genPartIdxMother, gp.pt, gp.eta, gp.phi, gp.p4().M(), gp.status]])
        for col in table.field_names:
            table.align[col] = "l"
        print(table)

    def lhe_logger(self, LHEPart):
        '''Print LHE particle history on output
        '''

        print("---LHEPart---")
        if self.nanoversion<15:
            table = PrettyTable(['i', 'pdgId', 'pT', 'eta', 'phi', 'status', 'incomingpz'])
        else:
            table = PrettyTable(['i', 'pdgId', 'mother', 'pT', 'eta', 'phi', 'status', 'incomingpz'])
        for i, Lp in enumerate(LHEPart):
            if self.nanoversion<15:
                table.add_row([self.format_value(val) for val in [i, Lp.pdgId, Lp.pt, Lp.eta, Lp.mass, Lp.status, Lp.incomingpz]])
            else:
                table.add_row([self.format_value(val) for val in [i, Lp.pdgId, Lp.firstMotherIdx, Lp.pt, Lp.eta, Lp.mass, Lp.status, Lp.incomingpz]])
        for col in table.field_names:
            table.align[col] = "l"
        print(table)


