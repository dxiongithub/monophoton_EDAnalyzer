import FWCore.ParameterSet.Config as cms

process = cms.Process("demo")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = -1
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.source = cms.Source("PoolSource",
# replace 'myfile.root' with the source file you want to use
fileNames = cms.untracked.vstring(
                                'file:FB1DF528-86EF-374F-8074-660836ACAC19.root'
)
)

process.demo = cms.EDAnalyzer('mvaAnalyzer',
photons = cms.InputTag("slimmedPhotons"),
OOTphotons = cms.InputTag("slimmedOOTPhotons"),
ebReducedRecHitCollection = cms.InputTag("ecalRecalibRecHit", "reducedEBRecHits"),
eeReducedRecHitCollection = cms.InputTag("ecalRecalibRecHit", "reducedEERecHits"),
esReducedRecHitCollection = cms.InputTag("ecalRecalibRecHit", "reducedESRecHits"),
triggerResults                      = cms.InputTag("TriggerResults", "", "HLT"),
muons = cms.InputTag("slimmedMuons"),
electrons = cms.InputTag("slimmedElectrons"),
mets = cms.InputTag("slimmedMETs"),
ecalBadCalibReducedMINIAODFilter    = cms.InputTag("ecalBadCalibReducedMINIAODFilter"),
vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
tracks = cms.InputTag("packedPFCandidates"),
genparticles = cms.bool(False),
miniaodinfo = cms.bool(True),
year     = cms.int32(2017)   #2017 and 2018 share the same trigger info 

)

process.TFileService = cms.Service("TFileService",
fileName = cms.string('monopho17.root')
)

process.p = cms.Path(process.demo)