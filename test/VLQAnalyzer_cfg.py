import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')
options.register('outFilename', 'test.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name"
                 )
options.register('skim', True,
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool,
                 "skim events with one lepton and 2 jets"
                 )
options.register('updateJEC', '',
                 VarParsing.multiplicity.list,
                 VarParsing.varType.string,
                 "Name of the SQLite file (with path and extension) used to update the jet collection to the latest JEC and the era of the new JEC"
                )
options.parseArguments()

process = cms.Process("TpAna")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
       #'/store/mc/PhaseIITDRFall17MiniAOD/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v2-v3/00000/0497DC9A-D3C1-E711-80BB-008CFA1C645C.root',
       '/store/mc/PhaseIITDRFall17MiniAOD/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/MINIAODSIM/noPU_93X_upgrade2023_realistic_v2-v1/150000/04B58981-6FBA-E711-B496-002590E7DF2A.root', 
       # '/store/mc/PhaseIITDRFall17MiniAOD/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v2-v1/150000/0E6628CF-AABE-E711-99B8-9CDC714AC590.root',
       # '/store/mc/PhaseIITDRFall17MiniAOD/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v2-v1/150000/56C4EC43-71BE-E711-8111-0025905B8594.root',
       # '/store/mc/PhaseIITDRFall17MiniAOD/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v2-v1/150000/8446A966-8CBE-E711-ABCA-34E6D7E05F0E.root',
       # '/store/mc/PhaseIITDRFall17MiniAOD/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v2-v1/150000/8E0420F9-1EBF-E711-B6E1-001C23C0B671.root',
       # '/store/mc/PhaseIITDRFall17MiniAOD/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v2-v1/150000/B06DB791-5ABE-E711-BBAB-002590E3A286.root',
       # '/store/mc/PhaseIITDRFall17MiniAOD/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v2-v1/150000/C276E7A3-43BE-E711-98F0-008CFAFBF1EC.root',
       # '/store/mc/PhaseIITDRFall17MiniAOD/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v2-v1/150000/C86BA597-55BE-E711-870C-003048CDBB94.root',
        ),
                            secondaryFileNames = cms.untracked.vstring(
       #'/store/mc/PhaseIITDRFall17DR/TT_TuneCUETP8M2T4_14TeV-powheg-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/30000/F24662E9-C2B0-E711-B263-0242AC130002.root',
       'file://AEE9EBD2-A1B9-E711-B422-0CC47A1E046E.root',
        #'/store/mc/PhaseIITDRFall17DR/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/150000/6430FCA9-91B8-E711-B58E-24BE05CE6D61.root',
        #'/store/mc/PhaseIITDRFall17DR/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/150000/8C648B49-42B9-E711-9F52-24BE05CE3EA1.root',
        #'/store/mc/PhaseIITDRFall17DR/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/150000/9A17A2BA-86B8-E711-8F16-24BE05C618F1.root',
        #'/store/mc/PhaseIITDRFall17DR/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/150000/9C6B9264-87B8-E711-A399-E0071B741D70.root',
        #'/store/mc/PhaseIITDRFall17DR/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/150000/A0C0DA74-8AB8-E711-A1A9-E0071B6CAD10.root',
        #'/store/mc/PhaseIITDRFall17DR/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/150000/B87037C6-8CB8-E711-A21A-24BE05C68671.root',
        #'/store/mc/PhaseIITDRFall17DR/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/150000/BA99979D-7AB8-E711-877B-E0071B73B6B0.root',
        #'/store/mc/PhaseIITDRFall17DR/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/150000/CE4A6F51-52B8-E711-9098-5065F3815241.root',
        )
)
        
# Geometry, GT, and other standard sequences
process.load('Configuration.Geometry.GeometryExtended2023D17Reco_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '93X_upgrade2023_realistic_v2', '')

# HGCAL EGamma ID
process.load("RecoEgamma.Phase2InterimID.phase2EgammaPAT_cff")

# load the input tags
process.ana = cms.EDAnalyzer("VLQAnalyzer")
process.load("Upgrades.VLQAnalyzer.CfiFile_cfi")

# load the output root file service routine
process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFilename))

# sequence the process
process.p = cms.Path(process.phase2Egamma * process.ana)

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 100
