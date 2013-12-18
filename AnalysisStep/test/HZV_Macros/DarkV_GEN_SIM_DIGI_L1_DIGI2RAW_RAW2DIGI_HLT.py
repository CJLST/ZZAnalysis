# Auto generated configuration file
# using: 
# Revision: 1.381.2.7 
# Source: /local/reps/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: reco -s GEN,SIM,DIGI,L1,DIGI2RAW,RAW2DIGI,HLT --conditions START53_V10::All --filein=file:/data2/p/pellicci/DarkBoson/Prod/Gen/HZV_mH126_mV10.root --fileout=HZV_mH126_mV10.root --mc --no_exec --scenario=pp
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.GeometrySimDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic8TeVCollision_cfi')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)

process.RandomNumberGeneratorService.generator = cms.PSet(
        initialSeed = cms.untracked.uint32(123456789),
        VtxSmeared = cms.untracked.uint32(223458),
        engineName = cms.untracked.string('TRandom3')
)

from Configuration.Generator.PythiaUESettings_cfi import *
process.generator = cms.EDFilter("Pythia6GeneratorFilter",
    pythiaPylistVerbosity = cms.untracked.int32(1),
    displayPythiaCards = cms.untracked.bool(False),                                 
    # put here the efficiency of your filter (1. if no filter)
    filterEfficiency = cms.untracked.double(1.0),
    pythiaHepMCVerbosity = cms.untracked.bool(True),
    # put here the cross section of your process (in pb)
    #crossSection = cms.untracked.double(0.5),
    comEnergy = cms.double(8000.0),
    maxEventsToPrint = cms.untracked.int32(2),
    PythiaParameters = cms.PSet(
        pythiaUESettingsBlock,
        processParameters = cms.vstring('PMAS(36,1)=126.0      !mass of H0',
            'PMAS(25,1)=30.0        !mass of h0',                                        
            'MSEL=0                  ! user selection for process', 
            'MSUB(157)=1             !ggH', 
            'MDME(420,1)=0           !Higgs decay into dd', 
            'MDME(421,1)=0           !Higgs decay into uu', 
            'MDME(422,1)=0           !Higgs decay into ss', 
            'MDME(423,1)=0           !Higgs decay into cc', 
            'MDME(424,1)=0           !Higgs decay into bb', 
            'MDME(425,1)=0           !Higgs decay into tt', 
            'MDME(426,1)=0           !Higgs decay into', 
            'MDME(427,1)=0           !Higgs decay into Higgs decay', 
            'MDME(428,1)=0           !Higgs decay into e nu e', 
            'MDME(429,1)=0           !Higgs decay into mu nu mu', 
            'MDME(430,1)=0           !Higgs decay into tau nu tau', 
            'MDME(431,1)=0           !Higgs decay into Higgs decay', 
            'MDME(432,1)=0           !Higgs decay into g g', 
            'MDME(433,1)=0           !Higgs decay into gam gam', 
            'MDME(434,1)=0           !Higgs decay into gam Z', 
            'MDME(435,1)=0           !Higgs decay into Z Z',
            'MDME(436,1)=0           !Higgs decay into W W',
            'MDME(437,1)=1           !Higgs decay into Z h0',
            'MDME(438,1)=0           !Higgs decay into h0 h0',
            'MDME(439,1)=0           !Higgs decay into A0 A0',                                                                                
            'MDME(210,1)=0           !h0 decay into dd', 
            'MDME(211,1)=0           !h0 decay into uu', 
            'MDME(212,1)=0           !h0 decay into ss', 
            'MDME(213,1)=0           !h0 decay into cc', 
            'MDME(214,1)=0           !h0 decay into bb', 
            'MDME(215,1)=0           !h0 decay into tt', 
            'MDME(218,1)=1           !h0 decay into e e', 
            'MDME(219,1)=1           !h0 decay into mu mu', 
            'MDME(220,1)=0           !h0 decay into tau tau', 
            'MDME(222,1)=0           !h0 decay into g g', 
            'MDME(223,1)=0           !h0 decay into gam gam', 
            'MDME(224,1)=0           !h0 decay into gam Z', 
            'MDME(225,1)=0           !h0 decay into Z Z', 
            'MDME(226,1)=0           !h0 decay into W W',
            'MDME(174,1)=0    !Z decay into d dbar', 
            'MDME(175,1)=0    !Z decay into u ubar', 
            'MDME(176,1)=0    !Z decay into s sbar', 
            'MDME(177,1)=0    !Z decay into c cbar', 
            'MDME(178,1)=0    !Z decay into b bbar', 
            'MDME(179,1)=0    !Z decay into t tbar', 
            'MDME(182,1)=1    !Z decay into e- e+', 
            'MDME(183,1)=0    !Z decay into nu_e nu_ebar', 
            'MDME(184,1)=1    !Z decay into mu- mu+', 
            'MDME(185,1)=0    !Z decay into nu_mu nu_mubar', 
            'MDME(186,1)=0    !Z decay into tau- tau+', 
            'MDME(187,1)=0    !Z decay into nu_tau nu_taubar',
            'BRAT(218) = 0.5    ',
            'BRAT(219) = 0.5    ',
            'MWID(25) = 2'),
        # This is a vector of ParameterSet names to be read, in this order
        parameterSets = cms.vstring('pythiaUESettings', 
            'processParameters')
    )
)

process.load("Configuration.StandardSequences.Generator_cff")

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.1 $'),
    annotation = cms.untracked.string('reco nevts:1'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    fileName = cms.untracked.string('GenSimDigiHLT.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'START53_V10::All', '')

# Path and EndPath definitions
process.generation_step = cms.Path(process.generator * process.pgen)
process.simulation_step = cms.Path(process.psim)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.simulation_step,process.digitisation_step,process.L1simulation_step,process.digi2raw_step,process.raw2digi_step)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.RECOSIMoutput_step])

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# End of customisation functions
