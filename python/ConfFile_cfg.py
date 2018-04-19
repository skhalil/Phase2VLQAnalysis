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
        '/store/mc/PhaseIITDRFall17MiniAOD/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/MINIAODSIM/noPU_93X_upgrade2023_realistic_v2-v1/150000/04B58981-6FBA-E711-B496-002590E7DF2A.root',        
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17MiniAOD/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/MINIAODSIM/PU200_93X_upgrade2023_realistic_v2-v2/20000/1E43B80D-3CC2-E711-9953-008CFA197D10.root', 
        ),
                            secondaryFileNames = cms.untracked.vstring(
        #'/store/mc/PhaseIITDRFall17DR/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/GEN-SIM-RECO/noPU_93X_upgrade2023_realistic_v2-v1/150000/AEE9EBD2-A1B9-E711-B422-0CC47A1E046E.root',
        'file:/uscms_data/d2/skhalil/AEE9EBD2-A1B9-E711-B422-0CC47A1E046E.root',
        #'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17DR/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/GEN-SIM-RECO/PU200_93X_upgrade2023_realistic_v2-v1/150000/3E8063E8-5DB8-E711-94A0-A4BF0112E310.root',
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
