import FWCore.ParameterSet.Config as cms

process = cms.Process("demo")
process.load('FWCore.MessageService.MessageLogger_cfi')

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff" )
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v11', '')

process.MessageLogger.cerr.FwkReport.reportEvery = -1
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.source = cms.Source("PoolSource",
# replace 'myfile.root' with the source file you want to use
fileNames = cms.untracked.vstring(
                                'file:703D04DA-BF2C-2B4A-906B-0EF17E153839.root'
                                ))

process.load( "PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff" )
process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
process.load( "PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff" )

### fix a bug in the ECAL-Tracker momentum combination when applying the scale and smearing
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runVID=True, 
                       era='2017-UL',
                       eleIDModules=['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Fall17_94X_V2_cff',
                                     'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
                                     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff',
                                     'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_noIso_V2_cff'],
                       phoIDModules=['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Fall17_94X_V2_cff',
                                     'RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff']
                       )

from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process,  names=['Photons', 'Electrons','Muons','Taus','Jets'], outputModules = [] )


from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD (
        process,
        isData = True, # false for MC    
        fixEE2017 = True,
        fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
        postfix = "ModifiedMET"
)

process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
                                   src = cms.InputTag("slimmedMuons"),
                                   preselection = cms.string("track.isNonnull"),
                                   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
                                   fractionOfSharedSegments = cms.double(0.499))

process.demo = cms.EDAnalyzer('photonAnalyzer',
photons                             = cms.InputTag("slimmedPhotons"),
OOTphotons                          = cms.InputTag("slimmedOOTPhotons"),
ebReducedRecHitCollection           = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
eeReducedRecHitCollection           = cms.InputTag("reducedEgamma", "reducedEERecHits"),
esReducedRecHitCollection           = cms.InputTag("reducedEgamma", "reducedESRecHits"),
triggerResults                      = cms.InputTag("TriggerResults", "", "HLT"),
patTriggerResults                   = cms.InputTag("TriggerResults", "", "RECO"),
muons                               = cms.InputTag("slimmedMuons"),
electrons                           = cms.InputTag("slimmedElectrons"),
mets                                = cms.InputTag("slimmedMETs"),
ecalBadCalibReducedMINIAODFilter    = cms.InputTag("ecalBadCalibReducedMINIAODFilter"),
vertices                            = cms.InputTag("offlineSlimmedPrimaryVertices"),
rhoLabel                            = cms.InputTag("fixedGridRhoFastjetAll"),
tracks                              = cms.InputTag("packedPFCandidates"),
genparticles                        = cms.bool(False),
miniaodinfo                         = cms.bool(True),
newParticles                        = cms.vint32(1000006, 1000021, 1000022, 1000024, 1000025, 1000039, 3000001, 3000002, 35),
runOnParticleGun                    = cms.bool(False),
genParticleSrc                      = cms.InputTag("prunedGenParticles"),
year                                = cms.int32(2017),   #2017 and 2018 share the same trigger info 


)

process.TFileService = cms.Service("TFileService",
fileName = cms.string('monopho17_genphostudy.root')
)

process.p = cms.Path(
    process.egammaPostRecoSeq *
    process.cleanedMu *
    process.demo)