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
    jets        = cms.InputTag("slimmedJetsPuppi"),#slimmedJets
    jets_ak8    = cms.InputTag("slimmedJetsAK8"),
    mets        = cms.InputTag("slimmedMETsPuppi"),#slimmedMETs
    subjets_ak8 = cms.InputTag("slimmedJets"),#slimmedJetsAK8PFPuppiSoftDropPacked_Subjets, this varaible is not used
    usePuppi    = cms.bool(True),# turn it off when not using Puppi jets/MET
    ak4ptmin    = cms.double(20.), 
    ak4etamax   = cms.double(5.), 
    ak8ptmin    = cms.double(200.), 
    ak8etamax   = cms.double(3.)
)
