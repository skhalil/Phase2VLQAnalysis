import FWCore.ParameterSet.Config as cms

process = cms.Process("Gen")


process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.MessageLogger.cerr.FwkJob.limit=1
process.MessageLogger.cerr.ERROR = cms.untracked.PSet( limit = cms.untracked.int32(0) )


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
       '/store/mc/PhaseIITDRFall17MiniAOD/TprimeBToTH_M-1000_Width-10p_LH_TuneCUETP8M2T4_14TeV-madgraph-pythia8/MINIAODSIM/noPU_93X_upgrade2023_realistic_v2-v1/150000/04B58981-6FBA-E711-B496-002590E7DF2A.root', 
        ),
                            secondaryFileNames = cms.untracked.vstring(
       'file://AEE9EBD2-A1B9-E711-B422-0CC47A1E046E.root',
        )
)

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.printTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("prunedGenParticles"),
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(False),
                                   printStatus = cms.untracked.bool(False),
                                   printIndex = cms.untracked.bool(True),
                                   status = cms.untracked.vint32( 3 )
                                   )

process.printList = cms.EDAnalyzer("ParticleListDrawer",
  maxEventsToPrint = cms.untracked.int32(1),
  printVertex = cms.untracked.bool(False),
  src = cms.InputTag("prunedGenParticles")
)

process.p0 = cms.Path(process.printTree * process.printList)

#process.out = cms.OutputModule("PoolOutputModule",
#                               fileName = cms.untracked.string('genOutputFile.root'),
#                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p0') ),
#                               outputCommands = cms.untracked.vstring('drop *')
#                               )

#process.e = cms.EndPath(process.out)


