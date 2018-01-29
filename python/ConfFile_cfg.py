import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')
options.register('outFilename', 'VLQTree.root',
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
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
       'root://cms-xrd-global.cern.ch//store/mc/PhaseIITDRFall17MiniAOD/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/MINIAODSIM/noPU_93X_upgrade2023_realistic_v2-v1/150000/04B58981-6FBA-E711-B496-002590E7DF2A.root', 
     #  'file:myfile.root'
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

process.ana = cms.EDAnalyzer("VLQAnalyzer",
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    puInfo   = cms.InputTag("slimmedAddPileupInfo"),
    #pileup = cms.uint32(200),
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outFilename))

process.p = cms.Path(process.ana)

process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    allowUnscheduled = cms.untracked.bool(True)
)

process.MessageLogger.cerr.FwkReport.reportEvery = 100
