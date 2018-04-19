import FWCore.ParameterSet.Config as cms

ana = cms.EDAnalyzer('VLQAnalyzer',
    genParts    = cms.InputTag("packedGenParticles"),    
    vertices    = cms.InputTag("offlineSlimmedPrimaryVertices"),
    puInfo      = cms.InputTag("slimmedAddPileupInfo"),
    electrons   = cms.InputTag("phase2Electrons"),
    muons       = cms.InputTag("slimmedMuons"),
    beamspot    = cms.InputTag("offlineBeamSpot"),
    conversions = cms.InputTag("reducedEgamma", "reducedConversions", "PAT"),
    genJets     = cms.InputTag("slimmedGenJets"),
    jets        = cms.InputTag("slimmedJets"),#slimmedJetsPuppi
    jets_ak8    = cms.InputTag("slimmedJetsAK8"),#slimmedJetsAK8Puppi
    mets        = cms.InputTag("slimmedMETs"),#slimmedMETsPuppi
    subjets_ak8  = cms.InputTag("slimmedJets"),
    ak4ptmin     = cms.double(20.), 
    ak4etamax    = cms.double(5.), 
    ak8ptmin     = cms.double(200.), 
    ak8etamax = cms.double(3.)
)
