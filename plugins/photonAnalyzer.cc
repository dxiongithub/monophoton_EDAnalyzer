// -*- C++ -*-
//
// Package:    monophoAnalysis/photonAnalyzer
// Class:      photonAnalyzer
//
/**\class photonAnalyzer photonAnalyzer.cc monophoAnalysis/photonAnalyzer/plugins/photonAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  De-Lin Xiong
//         Created:  Fri, 26 Mar 2021 16:07:56 GMT
//
//

/*
==================== 03262021 ====================
my own Ntuplizer
test in photon
... - great!

then track
...


vertex, electron, photon, muon, met


now getting everything
*/


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/Framework/interface/ESHandle.h"

//getting the object info (from .h)
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/DetId/interface/DetId.h"
#include <map>

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/METReco/interface/BeamHaloSummary.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TTree.h"




using namespace std;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


void setbit (UShort_t& x, UShort_t bit);

class photonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit photonAnalyzer(const edm::ParameterSet&);
      ~photonAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      

   private:
        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        void branchPhotons        (TTree*);
        void branchOOTPhotons     (TTree*);
        void branchCrystals       (TTree*);
        void branchTriggerinfo    (TTree*);
        void branchMuons          (TTree*);
        void branchElectrons      (TTree*);
        void branchMETs           (TTree*);
        void branchVertices       (TTree*);
        //void branchEvents         (TTree*);
        void branchTracks         (TTree*);
        //void branchGenInfo        (TTree*);
        void branchBeamHaloMuon   (TTree*);
        


        void fillPhotons      (const edm::Event&, const edm::EventSetup&);
        void fillOOTPhotons   (const edm::Event&, const edm::EventSetup&);
        void fillCrystals     (const edm::Event&, const edm::EventSetup&);
        void fillTriggerinfo  (const edm::Event&, const edm::EventSetup&);
        void fillMuons        (const edm::Event&, const edm::EventSetup&);
        void fillElectrons    (const edm::Event&, const edm::EventSetup&);
        void fillMETs         (const edm::Event&, const edm::EventSetup&);
        void fillVertices     (const edm::Event&, const edm::EventSetup&);
        //void fillEvents       (const edm::Event&, const edm::EventSetup&);
        void fillTracks       (const edm::Event&, const edm::EventSetup&);
        //void fillGenInfo      (const edm::Event&, const edm::EventSetup&);
        void fillBeamHaloMuon (const edm::Event&, const edm::EventSetup&);

        // ----------member data ---------------------------
        edm::EDGetTokenT<edm::View<pat::Photon> >           photonCollection_;
        edm::EDGetTokenT<edm::View<pat::Photon> >           OOTphotonCollection_;
        edm::EDGetTokenT<EcalRecHitCollection>              ebReducedRecHitCollection_;
        edm::EDGetTokenT<EcalRecHitCollection>              eeReducedRecHitCollection_;
        edm::EDGetTokenT<EcalRecHitCollection>              esReducedRecHitCollection_; 
        edm::EDGetTokenT<edm::TriggerResults>               trgResultsLabel_;
        edm::EDGetTokenT<edm::TriggerResults>               patTrgResultsLabel_;
        edm::EDGetTokenT<edm::View<pat::Muon> >             muonCollection_;
        edm::EDGetTokenT<edm::View<pat::Electron> >         electronCollection_;
        edm::EDGetTokenT<edm::View<pat::MET> >              pfMETlabel_;
        edm::EDGetTokenT<bool>                              ecalBadCalibFilterUpdate_;
        edm::EDGetTokenT<std::vector<reco::Vertex>>         vtxLabel_;
        edm::EDGetTokenT<double>                            rhoLabel_;
        edm::EDGetTokenT<std::vector<pat::PackedCandidate>> tracklabel_;
        bool doGenParticles_;
        bool addFilterInfoMINIAOD_;
        //edm::EDGetTokenT<std::vector<reco::GenParticle> >   genParticlesCollection_;
        edm::EDGetTokenT<reco::BeamHaloSummary>             beamHaloSummaryToken_;
        

        static bool vtxSort( const reco::Vertex &  a, const reco::Vertex & b );
        

        TTree   *tree_;
        
        //*********** photons ***********
        Int_t               nPho_;
        vector<float>       phoEt_;
        vector<float>       phoEta_;
        vector<float>       phoPhi_;
        vector<UShort_t>    phoIDbit_;
        vector<float>       phoSeedIEta_;
        vector<float>       phoSeedIPhi_;
        

        vector<int>         phohasPixelSeed_;
        vector<float>       phoR9_;
        vector<float>       phoSeedTime_;
        vector<float>       phoSeedEnergy_;
        vector<float>       phoMIPTotEnergy_;
        vector<float>       phoCalibEt_;
        

        //photon isolation
        vector<float>       phoPFChIso_;
        vector<float>       phoPFChPVIso_;
        vector<float>       phoPFPhoIso_;
        vector<float>       phoPFNeuIso_;
        vector<float>       phoPFChWorstVetoIso_;
        vector<float>       phoPFChWorstIso_;

        //photon mva variables
        vector<float>       phoIDMVA_;
        vector<float>       phoHoverE_;
        vector<float>       phoTrkSumPtHollowConeDR03_;
        vector<float>       phoPFClusEcalIso_;
        vector<float>       phoPFClusHcalIso_;

        vector<float>       phoE2x2Full5x5_;
        vector<float>       phoE1x3Full5x5_;
        vector<float>       phoE2ndFull5x5_;
        vector<float>       phoE2x5Full5x5_;
        vector<float>       phoMaxEnergyXtal_;
        vector<float>       phoR9Full5x5_;
        vector<float>       phoE5x5Full5x5_;
        vector<float>       phoSigmaIEtaIEtaFull5x5_;
        vector<float>       phoSigmaIEtaIPhiFull5x5_;
        vector<float>       phoSigmaIPhiIPhiFull5x5_;
        vector<float>       phoSCEtaWidth_;
        vector<float>       phoSCPhiWidth_;
        vector<float>       phoSCRawE_;


        //*********** ootphotons ***********
        Int_t               onPho_;
        vector<float>       ophoEt_;
        vector<float>       ophoEta_;
        vector<float>       ophoPhi_;
        vector<float>       ophoSeedIEta_;
        vector<float>       ophoSeedIPhi_;

        vector<int>         ophohasPixelSeed_;
        vector<float>       ophoR9_;
        vector<float>       ophoSeedTime_;
        vector<float>       ophoSeedEnergy_;
        vector<float>       ophoMIPTotEnergy_;
        vector<float>       ophoCalibEt_;
        

        //ophoton mva variables
        vector<float>       ophoHoverE_;
        vector<float>       ophoTrkSumPtHollowConeDR03_;
        vector<float>       ophoPFClusEcalIso_;
        vector<float>       ophoPFClusHcalIso_;

        vector<float>       ophoE2x2Full5x5_;
        vector<float>       ophoE1x3Full5x5_;
        vector<float>       ophoE2ndFull5x5_;
        vector<float>       ophoE2x5Full5x5_;
        vector<float>       ophoE5x5Full5x5_;
        vector<float>       ophoMaxEnergyXtal_;
        vector<float>       ophoR9Full5x5_;
        vector<float>       ophoSigmaIEtaIEtaFull5x5_;
        vector<float>       ophoSigmaIEtaIPhiFull5x5_;
        vector<float>       ophoSigmaIPhiIPhiFull5x5_;
        vector<float>       ophoSCEtaWidth_;
        vector<float>       ophoSCPhiWidth_;
        vector<float>       ophoSCRawE_;
        

        //photon crystals
        Int_t               nPFPhoton;
        Int_t               nAllCellsEB;
        Int_t               AllCellsIEtaEB[30000];
        Int_t               AllCellsIPhiEB[30000];
        Float_t             AllCellsE_EB[30000];
        Float_t             AllTimeEB[30000];
        Int_t               AllClusteredEB[30000][30];

        //triggerinfo
        ULong64_t           HLTPho_;

        //muons
        Int_t               nMu_;
        vector<float>       muPt_;
        vector<float>       muEta_;
        vector<float>       muPhi_;
        vector<int>         muIDbit_;

        //Electrons
        Int_t               nEle_;
        vector<float>       elePt_;
        vector<float>       eleEta_;
        vector<float>       elePhi_;
        vector<UShort_t>    eleIDbit_;


        
        //METs
        Int_t metFilters_;

        float genMET_;
        float genMETPhi_;
        float pfMET_;
        float pfMETPhi_;



        //vertexs
        Int_t nVtx_;
        vector<float> vtx_x_;
        vector<float> vtx_y_;
        vector<float> vtx_z_;
        vector<float> vtx_xError_;
        vector<float> vtx_yError_;
        vector<float> vtx_zError_;
        int           year_;


        //events variables
        Int_t       run_;
        Long64_t    event_;
        Int_t       lumis_;
        float       rho_;


        //tracks
        Int_t nTracks_;
        vector<float> trackPt_;
        vector<float> trackEta_;
        vector<float> trackPhi_;


        //beamhalomuon
        Int_t            nBHmuon_;
        vector<UShort_t> BHIDbit_;
/*
        //genInfo
        Int_t            nMC_;
        vector<int>      mcPID;
        vector<float>    mcPt;
        vector<float>    mcMass;
        vector<float>    mcEta;
        vector<float>    mcPhi;
        vector<float>    mcEt;
        vector<int>      mcStatus;


        Int_t            ngenPho_;
        vector<float>    genPhoPt;
        vector<float>    genPhoMass;
        vector<float>    genPhoEta;
        vector<float>    genPhoPhi;
        vector<float>    genPhoEt;
        vector<int>      genPhoStatus;


        vector<int>     newparticles_;
        bool runOnParticleGun_;
*/

        int N_vtx;
        double vtx_x, vtx_y, vtx_z, vtx_xError, vtx_yError, vtx_zError;


};

void setbit(UShort_t& x, UShort_t bit) 
{
    UShort_t a = 1;
    x |= (a << bit);
}

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
photonAnalyzer::photonAnalyzer(const edm::ParameterSet& ps):
    photonCollection_(consumes<edm::View<pat::Photon> >         (ps.getParameter<edm::InputTag>("photons"))),
    OOTphotonCollection_(consumes<edm::View<pat::Photon> >      (ps.getParameter<edm::InputTag>("OOTphotons"))),
    ebReducedRecHitCollection_(consumes<EcalRecHitCollection>   (ps.getParameter<edm::InputTag>("ebReducedRecHitCollection"))),
    eeReducedRecHitCollection_(consumes<EcalRecHitCollection>   (ps.getParameter<edm::InputTag>("eeReducedRecHitCollection"))),
    esReducedRecHitCollection_(consumes<EcalRecHitCollection>   (ps.getParameter<edm::InputTag>("esReducedRecHitCollection"))),
    trgResultsLabel_(consumes<edm::TriggerResults>              (ps.getParameter<edm::InputTag>("triggerResults"))),
    patTrgResultsLabel_(consumes<edm::TriggerResults>           (ps.getParameter<edm::InputTag>("patTriggerResults"))),
    muonCollection_(consumes<edm::View<pat::Muon> >             (ps.getParameter<edm::InputTag>("muons"))),
    electronCollection_(consumes<edm::View<pat::Electron> >     (ps.getParameter<edm::InputTag>("electrons"))),
    pfMETlabel_(consumes<edm::View<pat::MET> >                  (ps.getParameter<edm::InputTag>("mets"))),
    ecalBadCalibFilterUpdate_(consumes<bool>                    (ps.getParameter<edm::InputTag>("ecalBadCalibReducedMINIAODFilter"))), //MET filters 
    vtxLabel_(consumes<std::vector<reco::Vertex>>               (ps.getParameter<edm::InputTag>("vertices"))),
    rhoLabel_(consumes<double>                                  (ps.getParameter<edm::InputTag>("rhoLabel"))),
    tracklabel_(consumes<std::vector<pat::PackedCandidate>>     (ps.getParameter<edm::InputTag>("tracks"))),
    //doGenParticles_                                             (ps.getParameter<bool>("genparticles")),
    addFilterInfoMINIAOD_                                       (ps.getParameter<bool>("miniaodinfo")),
    beamHaloSummaryToken_(consumes<reco::BeamHaloSummary>       (ps.getParameter<edm::InputTag>("beamhalomuon"))),
    //genParticlesCollection_(consumes<std::vector<reco::GenParticle>> (ps.getParameter<edm::InputTag>("genParticleSrc"))),
    year_                                                       (ps.getParameter<int>("year"))
    //newparticles_                                               (ps.getParameter<std::vector<int>>("newParticles")),
    //runOnParticleGun_                                           (ps.getParameter<bool>("runOnParticleGun"))

{
    //now do what ever initialization is needed
    
    edm::Service<TFileService> fs;
    tree_    = fs->make<TTree>("EventTree", "eventree");

    branchPhotons(tree_);
    branchOOTPhotons(tree_);
    branchCrystals(tree_);
    branchTriggerinfo(tree_);
    branchMuons(tree_);
    branchElectrons(tree_);
    branchMETs(tree_);
    branchVertices(tree_);
    branchTracks(tree_);
    //branchGenInfo(tree_);
    branchBeamHaloMuon(tree_);

}


photonAnalyzer::~photonAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


bool photonAnalyzer::vtxSort( const reco::Vertex &  a, const reco::Vertex & b ){
  if( a.tracksSize() != b.tracksSize() )
    return  a.tracksSize() > b.tracksSize() ? true : false ;
  else
    return  a.chi2() < b.chi2() ? true : false ;
}

//Necessary for the Photon Footprint removal
template <class T, class U>
bool isInFootprint(const T& thefootprint, const U& theCandidate) {
  for ( auto itr = thefootprint.begin(); itr != thefootprint.end(); ++itr ) {

    if( itr.key() == theCandidate.key() ) return true;
    
  }
  return false;
}

void photonAnalyzer::branchPhotons(TTree* tree)
{
    tree->Branch("nPho",                    &nPho_);
    tree->Branch("phoEt",                   &phoEt_);
    tree->Branch("phoEta",                  &phoEta_);
    tree->Branch("phoPhi",                  &phoPhi_);
    tree->Branch("phoIDbit",                &phoIDbit_);
    tree->Branch("phoSeedIEta",             &phoSeedIEta_);
    tree->Branch("phoSeedIPhi",             &phoSeedIPhi_);
    

    tree->Branch("phohasPixelSeed",         &phohasPixelSeed_);
    tree->Branch("phoR9",                   &phoR9_);
    tree->Branch("phoSeedTime",             &phoSeedTime_);
    tree->Branch("phoSeedEnergy",           &phoSeedEnergy_);
    tree->Branch("phoMIPTotEnergy",         &phoMIPTotEnergy_);
    tree->Branch("phoCalibEt",              &phoCalibEt_);


    tree->Branch("phoPFChIso",              &phoPFChIso_);
    tree->Branch("phoPFChPVIso",            &phoPFChPVIso_);
    tree->Branch("phoPFPhoIso",             &phoPFPhoIso_);
    tree->Branch("phoPFNeuIso",             &phoPFNeuIso_);
    tree->Branch("phoPFChWorstIso",         &phoPFChWorstIso_);
    tree->Branch("phoPFChWorstVetoIso",     &phoPFChWorstVetoIso_);
    
    tree->Branch("phoIDMVA",                &phoIDMVA_);
    tree->Branch("phoHoverE",               &phoHoverE_);
    tree->Branch("phoTrkSumPtHollowConeDR03", &phoTrkSumPtHollowConeDR03_);
    tree->Branch("phoPFClusEcalIso", &phoPFClusEcalIso_);
    tree->Branch("phoPFClusHcalIso", &phoPFClusHcalIso_);

    tree->Branch("phoE2x2Full5x5", &phoE2x2Full5x5_);
    tree->Branch("phoE2ndFull5x5", &phoE2ndFull5x5_);
    tree->Branch("phoE2x5Full5x5", &phoE2x5Full5x5_);
    tree->Branch("phoMaxEnergyXtal", &phoMaxEnergyXtal_);
    tree->Branch("phoR9Full5x5", &phoR9Full5x5_); 
    tree->Branch("phoE5x5Full5x5",          &phoE5x5Full5x5_);
    tree->Branch("phoE1x3Full5x5", &phoE1x3Full5x5_);
    tree->Branch("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5_);
    tree->Branch("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5_);
    tree->Branch("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5_);
    tree->Branch("phoSCEtaWidth", &phoSCEtaWidth_);
    tree->Branch("phoSCPhiWidth", &phoSCPhiWidth_);
    tree->Branch("phoSCRawE", &phoSCRawE_);
}


void photonAnalyzer::fillPhotons(const edm::Event& e, const edm::EventSetup& es) 
{

    using namespace edm;
    using namespace std;

    phoEt_          .clear();
    phoEta_         .clear();
    phoPhi_         .clear();
    phoIDbit_       .clear();
    phoSeedIEta_    .clear();
    phoSeedIPhi_   .clear();

    phohasPixelSeed_    .clear();
    phoR9_              .clear();
    phoSeedTime_        .clear();
    phoSeedEnergy_      .clear();
    phoMIPTotEnergy_    .clear(); 
    phoCalibEt_         .clear();


    phoPFChIso_         .clear();
    phoPFChPVIso_       .clear();
    phoPFPhoIso_        .clear();
    phoPFNeuIso_        .clear();
    phoPFChWorstVetoIso_.clear();
    phoPFChWorstIso_    .clear();
    

   
    phoIDMVA_               .clear();
    phoHoverE_              .clear();
    phoTrkSumPtHollowConeDR03_ .clear();
    phoPFClusEcalIso_       .clear();
    phoPFClusHcalIso_       .clear();

    phoE2x2Full5x5_         .clear();
    phoE1x3Full5x5_         .clear();
    phoE2ndFull5x5_         .clear();
    phoE2x5Full5x5_         .clear();
    phoMaxEnergyXtal_       .clear();
    phoR9Full5x5_           .clear();
    phoE5x5Full5x5_         .clear();
    phoSigmaIEtaIEtaFull5x5_ .clear();
    phoSigmaIEtaIPhiFull5x5_ .clear();
    phoSigmaIPhiIPhiFull5x5_ .clear();
    phoSCEtaWidth_          .clear();
    phoSCPhiWidth_          .clear();
    phoSCRawE_              .clear();

    nPho_ = 0;

    //photon info
    edm::Handle<edm::View<pat::Photon> > photonHandle;
    e.getByToken(photonCollection_, photonHandle);
    
    EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
    noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

    for (edm::View<pat::Photon>::const_iterator iPho = photonHandle->begin(); iPho != photonHandle->end(); ++iPho) 
    {
        phoEt_              .push_back(iPho->et());
        phoEta_             .push_back(iPho->eta());
        phoPhi_             .push_back(iPho->phi());

        UShort_t tmpphoIDbit = 0;        
        bool isPassLoose  = iPho->photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
        if (isPassLoose)  setbit(tmpphoIDbit, 0);   
        bool isPassMedium = iPho->photonID("cutBasedPhotonID-Fall17-94X-V2-medium");
        if (isPassMedium) setbit(tmpphoIDbit, 1);    
        bool isPassTight  = iPho->photonID("cutBasedPhotonID-Fall17-94X-V2-tight");
        if (isPassTight)  setbit(tmpphoIDbit, 2);
        
        phoIDbit_.push_back(tmpphoIDbit);      


        phohasPixelSeed_    .push_back((Int_t)iPho->hasPixelSeed());
        phoR9_              .push_back(iPho->r9());

        //seed info
        DetId seed = (iPho->superCluster()->seed()->hitsAndFractions())[0].first;
        bool isBarrel = seed.subdetId() == EcalBarrel;
        const EcalRecHitCollection * rechits = (isBarrel?lazyTool.getEcalEBRecHitCollection():lazyTool.getEcalEERecHitCollection());
                
        EcalRecHitCollection::const_iterator theSeedHit = rechits->find(seed);
        if (theSeedHit != rechits->end()) 
        {

            phoSeedTime_  .push_back((*theSeedHit).time());
            phoSeedEnergy_.push_back((*theSeedHit).energy());

            EBDetId dit   = theSeedHit->detid();
            phoSeedIEta_  .push_back(dit.ieta());
            phoSeedIPhi_  .push_back(dit.iphi());
        } 
        else
        {
            phoSeedTime_  .push_back(-99.);
            phoSeedEnergy_.push_back(-99.);

            phoSeedIEta_  .push_back(-99.);
            phoSeedIPhi_  .push_back(-99.);
        }

        phoMIPTotEnergy_         .push_back(iPho->mipTotEnergy());

        phoCalibEt_       .push_back(iPho->et()*iPho->userFloat("ecalEnergyPostCorr")/iPho->energy());
        phoIDMVA_           .push_back(iPho->userFloat("PhotonMVAEstimatorRunIIFall17v2Values"));  
        phoHoverE_          .push_back(iPho->hadTowOverEm());
        phoTrkSumPtHollowConeDR03_ .push_back(iPho->trkSumPtHollowConeDR03());  
        phoPFClusEcalIso_ .push_back(iPho->ecalPFClusterIso());
        phoPFClusHcalIso_ .push_back(iPho->hcalPFClusterIso());
        
        phoPFChIso_       .push_back(iPho->chargedHadronIso()); 
        phoPFChPVIso_       .push_back(iPho->chargedHadronPFPVIso()); 
        phoPFPhoIso_        .push_back(iPho->photonIso());
        phoPFNeuIso_        .push_back(iPho->neutralHadronIso());
        phoPFChWorstIso_    .push_back(iPho->chargedHadronWorstVtxIso()); 
        phoPFChWorstVetoIso_.push_back(iPho->chargedHadronWorstVtxGeomVetoIso());

        phoE2x2Full5x5_ .push_back(iPho->full5x5_showerShapeVariables().e2x2);
        phoE1x3Full5x5_ .push_back(iPho->full5x5_showerShapeVariables().e1x3);
        phoE2ndFull5x5_ .push_back(iPho->full5x5_showerShapeVariables().e2nd);
        phoE2x5Full5x5_ .push_back(iPho->full5x5_showerShapeVariables().e2x5);
        phoE5x5Full5x5_ .push_back(iPho->full5x5_e5x5());
        phoMaxEnergyXtal_ .push_back(iPho->full5x5_maxEnergyXtal());
        phoR9Full5x5_ .push_back(iPho->full5x5_r9());
        phoSigmaIEtaIEtaFull5x5_ .push_back(iPho->full5x5_sigmaIetaIeta());
        phoSigmaIEtaIPhiFull5x5_ .push_back(iPho->full5x5_showerShapeVariables().sigmaIetaIphi);
        phoSigmaIPhiIPhiFull5x5_ .push_back(iPho->full5x5_showerShapeVariables().sigmaIphiIphi);
        phoSCEtaWidth_ .push_back((*iPho).superCluster()->etaWidth());
        phoSCPhiWidth_ .push_back((*iPho).superCluster()->phiWidth());
        phoSCRawE_ .push_back((*iPho).superCluster()->rawEnergy());  
        
        nPho_++;
    }
}


void photonAnalyzer::branchOOTPhotons(TTree* tree)
{
    tree->Branch("onPho",        &onPho_);
    tree->Branch("ophoEt",       &ophoEt_);
    tree->Branch("ophoEta",      &ophoEta_);
    tree->Branch("ophoPhi",      &ophoPhi_);

    tree->Branch("ophoSeedIEta",            &ophoSeedIEta_);
    tree->Branch("ophoSeedIPhi",            &ophoSeedIPhi_);

    tree->Branch("ophohasPixelSeed",         &ophohasPixelSeed_);
    tree->Branch("ophoR9",                   &ophoR9_);
    tree->Branch("ophoSeedTime",             &ophoSeedTime_);
    tree->Branch("ophoSeedEnergy",           &ophoSeedEnergy_);
    tree->Branch("ophoMIPTotEnergy",         &ophoMIPTotEnergy_);
    tree->Branch("ophoCalibEt",              &ophoCalibEt_);


    tree->Branch("ophoHoverE",               &ophoHoverE_);
    tree->Branch("ophoTrkSumPtHollowConeDR03", &ophoTrkSumPtHollowConeDR03_);
    tree->Branch("ophoPFClusEcalIso", &ophoPFClusEcalIso_);
    tree->Branch("ophoPFClusHcalIso", &ophoPFClusHcalIso_);

    tree->Branch("ophoE2x2Full5x5", &ophoE2x2Full5x5_);
    tree->Branch("ophoE2ndFull5x5", &ophoE2ndFull5x5_);
    tree->Branch("ophoE2x5Full5x5", &ophoE2x5Full5x5_);
    tree->Branch("ophoMaxEnergyXtal", &ophoMaxEnergyXtal_);
    tree->Branch("ophoR9Full5x5", &ophoR9Full5x5_); 
    tree->Branch("ophoE1x3Full5x5", &ophoE1x3Full5x5_);
    tree->Branch("ophoE5x5Full5x5",          &ophoE5x5Full5x5_);
    tree->Branch("ophoSigmaIEtaIEtaFull5x5", &ophoSigmaIEtaIEtaFull5x5_);
    tree->Branch("ophoSigmaIEtaIPhiFull5x5", &ophoSigmaIEtaIPhiFull5x5_);
    tree->Branch("ophoSigmaIPhiIPhiFull5x5", &ophoSigmaIPhiIPhiFull5x5_);
    tree->Branch("ophoSCEtaWidth", &ophoSCEtaWidth_);
    tree->Branch("ophoSCPhiWidth", &ophoSCPhiWidth_);
    tree->Branch("ophoSCRawE", &ophoSCRawE_);
}

void photonAnalyzer::fillOOTPhotons(const edm::Event& e, const edm::EventSetup& es) 
{

    using namespace edm;
    using namespace std;

    ophoEt_          .clear();
    ophoEta_         .clear();
    ophoPhi_         .clear();

    ophoSeedIEta_    .clear();
    ophoSeedIPhi_   .clear();

    ophohasPixelSeed_    .clear();
    ophoR9_              .clear();
    ophoSeedTime_        .clear();
    ophoSeedEnergy_      .clear();
    ophoMIPTotEnergy_    .clear();  
    ophoCalibEt_         .clear();
   

    ophoHoverE_                 .clear();
    ophoTrkSumPtHollowConeDR03_ .clear();
    ophoPFClusEcalIso_          .clear();
    ophoPFClusHcalIso_          .clear();

    ophoE2x2Full5x5_         .clear();
    ophoE1x3Full5x5_         .clear();
    ophoE2ndFull5x5_         .clear();
    ophoE2x5Full5x5_         .clear();
    ophoMaxEnergyXtal_       .clear();
    ophoR9Full5x5_           .clear();
    ophoE5x5Full5x5_         .clear();
    ophoSigmaIEtaIEtaFull5x5_ .clear();
    ophoSigmaIEtaIPhiFull5x5_ .clear();
    ophoSigmaIPhiIPhiFull5x5_ .clear();
    ophoSCEtaWidth_          .clear();
    ophoSCPhiWidth_          .clear();
    ophoSCRawE_              .clear();

    onPho_ = 0;

    //photon info
    edm::Handle<edm::View<pat::Photon> > OOTPhotonsH;
    e.getByToken(OOTphotonCollection_, OOTPhotonsH);

    EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
    noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
    

    for (edm::View<pat::Photon>::const_iterator oPho = OOTPhotonsH->begin(); oPho != OOTPhotonsH->end(); ++oPho) 
    {
        ophoEt_              .push_back(oPho->et());
        ophoEta_             .push_back(oPho->eta());
        ophoPhi_             .push_back(oPho->phi());

        ophohasPixelSeed_    .push_back((Int_t)oPho->hasPixelSeed());
        ophoR9_              .push_back(oPho->r9());

        //seed info
        DetId seed = (oPho->superCluster()->seed()->hitsAndFractions())[0].first;
        bool isBarrel = seed.subdetId() == EcalBarrel;
        const EcalRecHitCollection * rechits = (isBarrel?lazyTool.getEcalEBRecHitCollection():lazyTool.getEcalEERecHitCollection());
                
        EcalRecHitCollection::const_iterator theSeedHit = rechits->find(seed);
        if (theSeedHit != rechits->end()) 
        {
            ophoSeedTime_  .push_back((*theSeedHit).time());
            ophoSeedEnergy_.push_back((*theSeedHit).energy());

            EBDetId dit   = theSeedHit->detid();
            ophoSeedIEta_  .push_back(dit.ieta());
            ophoSeedIPhi_  .push_back(dit.iphi());
        } 
        else
        {
            ophoSeedTime_  .push_back(-99.);
            ophoSeedEnergy_.push_back(-99.);

            ophoSeedIEta_  .push_back(-99.);
            ophoSeedIPhi_  .push_back(-99.);
        }

        ophoMIPTotEnergy_         .push_back(oPho->mipTotEnergy());
        //ophoCalibEt_              .push_back(oPho->et()*oPho->userFloat("ecalEnergyPostCorr")/oPho->energy());

        ophoHoverE_          .push_back(oPho->hadTowOverEm());
        ophoTrkSumPtHollowConeDR03_ .push_back(oPho->trkSumPtHollowConeDR03());  
        ophoPFClusEcalIso_ .push_back(oPho->ecalPFClusterIso());
        ophoPFClusHcalIso_ .push_back(oPho->hcalPFClusterIso());

        ophoE2x2Full5x5_ .push_back(oPho->full5x5_showerShapeVariables().e2x2);
        ophoE1x3Full5x5_ .push_back(oPho->full5x5_showerShapeVariables().e1x3);
        ophoE2ndFull5x5_ .push_back(oPho->full5x5_showerShapeVariables().e2nd);
        ophoE2x5Full5x5_ .push_back(oPho->full5x5_showerShapeVariables().e2x5);
        ophoE5x5Full5x5_  .push_back(oPho->full5x5_e5x5());
        ophoMaxEnergyXtal_ .push_back(oPho->full5x5_maxEnergyXtal());
        ophoR9Full5x5_ .push_back(oPho->full5x5_r9());
        ophoSigmaIEtaIEtaFull5x5_ .push_back(oPho->full5x5_sigmaIetaIeta());
        ophoSigmaIEtaIPhiFull5x5_ .push_back(oPho->full5x5_showerShapeVariables().sigmaIetaIphi);
        ophoSigmaIPhiIPhiFull5x5_ .push_back(oPho->full5x5_showerShapeVariables().sigmaIphiIphi);
        ophoSCEtaWidth_ .push_back((*oPho).superCluster()->etaWidth());
        ophoSCPhiWidth_ .push_back((*oPho).superCluster()->phiWidth());
        ophoSCRawE_ .push_back((*oPho).superCluster()->rawEnergy());  
        
        onPho_++;
    }
}


void photonAnalyzer::branchCrystals(TTree* tree)
{
    tree->Branch("nPFPhoton", &nPFPhoton, "nPFPhoton/I");
    tree->Branch("nAllCellsEB", &nAllCellsEB, "nAllCellsEB/I");
    tree->Branch("AllCellsIEtaEB", &AllCellsIEtaEB, "AllCellsIEtaEB[nAllCellsEB]/I");
    tree->Branch("AllCellsIPhiEB", &AllCellsIPhiEB, "AllCellsIPhiEB[nAllCellsEB]/I");  
    tree->Branch("AllCellsE_EB", &AllCellsE_EB, "AllCellsE_EB[nAllCellsEB]/F");
    tree->Branch("AllTimeEB", &AllTimeEB, "AllFracEB[nAllCellsEB]/F");
    tree->Branch("AllClusteredEB", AllClusteredEB, "AllClusteredEB[nAllCellsEB][30]/I");
}

void photonAnalyzer::fillCrystals(const edm::Event& e, const edm::EventSetup& es) 
{
    edm::Handle<edm::View<pat::Photon> > photonHandle;
    e.getByToken(photonCollection_, photonHandle);
    const edm::View<pat::Photon> *mphotons = photonHandle.product();
    edm::View<pat::Photon>::const_iterator photonit = mphotons->begin();

    std::map<DetId, vector<std::pair<Int_t, Float_t> > > crysclusEB;
    nPFPhoton = 0;
    for(photonit = mphotons->begin(); photonit!= mphotons->end()&&nPFPhoton<100; ++photonit)
    {
        std::vector<std::pair<DetId, float> > clusdet = photonit->superCluster()->hitsAndFractions();
        for (int iii=0;iii<int(clusdet.size());++iii)
        {
            DetId blarg = clusdet[iii].first;
            
            std::map<DetId, vector< pair<Int_t, Float_t> > >::iterator selcrys = crysclusEB.find(blarg); 
            
            if (selcrys == crysclusEB.end())
            {
                vector< std::pair<Int_t, Float_t> > blerg;
                blerg.push_back( make_pair(nPFPhoton, clusdet[iii].second));
                crysclusEB.insert(make_pair(blarg,blerg));
            }
            else
            {
                selcrys->second.push_back(make_pair(nPFPhoton, clusdet[iii].second));
            }
        }
    }

    //getting the info of detid
    edm::Handle<EcalRecHitCollection> ecalhitsCollHEB;
    e.getByToken(ebReducedRecHitCollection_, ecalhitsCollHEB);
    const EcalRecHitCollection* rechitsCollectionEB_ = ecalhitsCollHEB.product();

    
    //AllClusteredEB[30000][30] = {0};
    nAllCellsEB=0;

    for (EcalRecHitCollection::const_iterator it = rechitsCollectionEB_->begin();it!=rechitsCollectionEB_->end() 
    && nAllCellsEB < 30000; ++it)
    {
        DetId blarg = it->detid();
        //Here's where the map is used to put what photon this crystal goes with. 
        std::map<DetId, vector< pair<Int_t, Float_t> > >::const_iterator selcrys = crysclusEB.find(blarg);
        if (selcrys==crysclusEB.end())
        {
            AllClusteredEB[nAllCellsEB][0]=-1;
        }
        else
        {
            vector< pair<Int_t, Float_t> > blerg = selcrys->second;
            for (int ikk=0; ikk<int(blerg.size()) && ikk<30; ikk++)
            {
                AllClusteredEB[nAllCellsEB][ikk]=blerg[ikk].first;
            }
        }
        
        //array
        EBDetId dit = it->detid();
        AllCellsIEtaEB[nAllCellsEB]=dit.ieta();
        AllCellsIPhiEB[nAllCellsEB]=dit.iphi();
        AllCellsE_EB[nAllCellsEB]=it->energy();
        AllTimeEB[nAllCellsEB]=it->time();
        

        nAllCellsEB++;

    }
}


void photonAnalyzer::branchTriggerinfo(TTree* tree)
{
    tree->Branch("HLTPho",      &HLTPho_);
}

void photonAnalyzer::fillTriggerinfo(const edm::Event& e, const edm::EventSetup& es) 
{
    HLTPho_ = 0;

    edm::Handle<edm::TriggerResults> trgResultsHandle;
    e.getByToken(trgResultsLabel_, trgResultsHandle);

    const edm::TriggerNames &trgNames = e.triggerNames(*trgResultsHandle);

    for (size_t i = 0; i < trgNames.size(); ++i) 
    {
        const std::string &name = trgNames.triggerName(i);
        int bitPho    = -1;

        if (year_ == 2017) 
        {
            if      (name.find("HLT_Photon175_v")                   != std::string::npos) bitPho =  7; // 2017
            else if (name.find("HLT_Photon250_v")                   != std::string::npos) bitPho =  8; // 2017
            else if (name.find("HLT_Photon300_NoHE_v")              != std::string::npos) bitPho =  9; // 2017
            else if (name.find("HLT_Photon200_v")                   != std::string::npos) bitPho = 10; // 2017
            else if (name.find("HLT_Photon165_HE10_v")              != std::string::npos) bitPho = 12; 
        }
        
        ULong64_t isFired     = (trgResultsHandle->accept(i)) ? 1 : 0;

        if (bitPho >= 0) 
        {
            HLTPho_            |= (isFired << bitPho);
        }

    }
}


void photonAnalyzer::branchMuons(TTree* tree)
{
    tree->Branch("nMu",     &nMu_);
    tree->Branch("muPt",    &muPt_);
    tree->Branch("muEta",   &muEta_);
    tree->Branch("muPhi",   &muPhi_);
    tree->Branch("muIDbit", &muIDbit_);
}


void photonAnalyzer::fillMuons(const edm::Event& e, const edm::EventSetup& es) 
{
    muPt_                  .clear();
    muEta_                 .clear();
    muPhi_                 .clear();
    muIDbit_               .clear();

    nMu_ = 0;

    edm::Handle<edm::View<pat::Muon> > muonHandle;
    e.getByToken(muonCollection_, muonHandle);

    for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) 
    {
        muPt_    .push_back(iMu->pt());
        muEta_   .push_back(iMu->eta());
        muPhi_   .push_back(iMu->phi());

        int tmpmuIDbit = 0;
        if (iMu->passed(reco::Muon::CutBasedIdLoose))        tmpmuIDbit += pow(2,  0);
        if (iMu->passed(reco::Muon::CutBasedIdMedium))       tmpmuIDbit += pow(2,  1);
        if (iMu->passed(reco::Muon::CutBasedIdMediumPrompt)) tmpmuIDbit += pow(2,  2);
        if (iMu->passed(reco::Muon::CutBasedIdTight))        tmpmuIDbit += pow(2,  3);
        if (iMu->passed(reco::Muon::CutBasedIdGlobalHighPt)) tmpmuIDbit += pow(2,  4);
        if (iMu->passed(reco::Muon::CutBasedIdTrkHighPt))    tmpmuIDbit += pow(2,  5);
        if (iMu->passed(reco::Muon::PFIsoVeryLoose))         tmpmuIDbit += pow(2,  6);
        if (iMu->passed(reco::Muon::PFIsoLoose))             tmpmuIDbit += pow(2,  7);
        if (iMu->passed(reco::Muon::PFIsoMedium))            tmpmuIDbit += pow(2,  8);
        if (iMu->passed(reco::Muon::PFIsoTight))             tmpmuIDbit += pow(2,  9);
        if (iMu->passed(reco::Muon::PFIsoVeryTight))         tmpmuIDbit += pow(2, 10);
        if (iMu->passed(reco::Muon::TkIsoLoose))             tmpmuIDbit += pow(2, 11);
        if (iMu->passed(reco::Muon::TkIsoTight))             tmpmuIDbit += pow(2, 12);
        if (iMu->passed(reco::Muon::SoftCutBasedId))         tmpmuIDbit += pow(2, 13);
        if (iMu->passed(reco::Muon::SoftMvaId))              tmpmuIDbit += pow(2, 14);
        if (iMu->passed(reco::Muon::MvaLoose))               tmpmuIDbit += pow(2, 15);
        if (iMu->passed(reco::Muon::MvaMedium))              tmpmuIDbit += pow(2, 16);
        if (iMu->passed(reco::Muon::MvaTight))               tmpmuIDbit += pow(2, 17);
        if (iMu->passed(reco::Muon::MiniIsoLoose))           tmpmuIDbit += pow(2, 18);
        if (iMu->passed(reco::Muon::MiniIsoMedium))          tmpmuIDbit += pow(2, 19);
        if (iMu->passed(reco::Muon::MiniIsoTight))           tmpmuIDbit += pow(2, 20);
        if (iMu->passed(reco::Muon::MiniIsoVeryTight))       tmpmuIDbit += pow(2, 21);
        muIDbit_.push_back(tmpmuIDbit);


        nMu_++;
    }

}



void photonAnalyzer::branchElectrons(TTree* tree)
{
    tree->Branch("nEle",        &nEle_);
    tree->Branch("elePt",       &elePt_);
    tree->Branch("eleEta",      &eleEta_);
    tree->Branch("elePhi",      &elePhi_);
    tree->Branch("eleIDbit",    &eleIDbit_);
}

void photonAnalyzer::fillElectrons(const edm::Event& e, const edm::EventSetup& es) 
{
    elePt_                      .clear();
    eleEta_                     .clear();
    elePhi_                     .clear();
    eleIDbit_                   .clear();

    nEle_ = 0;


    edm::Handle<edm::View<pat::Electron> > electronHandle;
    e.getByToken(electronCollection_, electronHandle);

    for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) 
    {
        elePt_              .push_back(iEle->pt());
        eleEta_             .push_back(iEle->eta());
        elePhi_             .push_back(iEle->phi());


        //EleID
        UShort_t tmpeleIDbit = 0;   
        bool isPassVeto   = iEle->electronID("cutBasedElectronID-Fall17-94X-V2-veto");
        if (isPassVeto)   setbit(tmpeleIDbit, 0);    
        bool isPassLoose  = iEle->electronID("cutBasedElectronID-Fall17-94X-V2-loose");
        if (isPassLoose)  setbit(tmpeleIDbit, 1);   
        bool isPassMedium = iEle->electronID("cutBasedElectronID-Fall17-94X-V2-medium");
        if (isPassMedium) setbit(tmpeleIDbit, 2);    
        bool isPassTight  = iEle->electronID("cutBasedElectronID-Fall17-94X-V2-tight");
        if (isPassTight)  setbit(tmpeleIDbit, 3);    
        bool isPassHEEP   = iEle->electronID("heepElectronID-HEEPV70");
        if (isPassHEEP)   setbit(tmpeleIDbit, 4);

        eleIDbit_.push_back(tmpeleIDbit);

        nEle_++;
    }

}

void photonAnalyzer::branchMETs(TTree* tree)
{
    
    if (doGenParticles_) 
    {
        tree->Branch("genMET",      &genMET_);
        tree->Branch("genMETPhi",   &genMETPhi_);
    }
    if (addFilterInfoMINIAOD_)
    {
        tree->Branch("metFilters",      &metFilters_);

        tree->Branch("pfMET",            &pfMET_);
        tree->Branch("pfMETPhi",         &pfMETPhi_);
    }
}


void photonAnalyzer::fillMETs(const edm::Event& e, const edm::EventSetup& es) 
{
    metFilters_ = 0;

    if (addFilterInfoMINIAOD_) 
    {
        std::string filterNamesToCheck[9] = 
        {
            "Flag_HBHENoiseFilter",
            "Flag_HBHENoiseIsoFilter", 
            "Flag_globalSuperTightHalo2016Filter",
            "Flag_goodVertices",
            "Flag_eeBadScFilter",
            "Flag_EcalDeadCellTriggerPrimitiveFilter",
            "Flag_BadPFMuonFilter",
            "Flag_ecalBadCalibReducedMINIAODFilter",
            "Flag_BadChargedCandidateFilter"
        };


        edm::Handle<edm::TriggerResults> patFilterResultsHandle;
        e.getByToken(patTrgResultsLabel_, patFilterResultsHandle);
        edm::TriggerResults const& patFilterResults = *patFilterResultsHandle;
        
        auto&& filterNames = e.triggerNames(patFilterResults);

        for (unsigned iF = 0; iF < 9; ++iF) 
        {
            unsigned index = filterNames.triggerIndex(filterNamesToCheck[iF]);
            if ( index == filterNames.size() ) 
            {
                edm::Handle<bool> passecalBadCalibFilterUpdate;
                e.getByToken(ecalBadCalibFilterUpdate_, passecalBadCalibFilterUpdate);
                if (passecalBadCalibFilterUpdate.isValid()) 
                {
                    bool passecalBadCalibFilterUpdate_ = (*passecalBadCalibFilterUpdate);	
                    if (passecalBadCalibFilterUpdate_) metFilters_ += pow(2, iF+1);
                }
            } 
            else 
            {
                if ( !patFilterResults.accept(index) ) 
                {
                    metFilters_ += pow(2, iF+1);
                }
            }
        }
    }

    edm::Handle<edm::View<pat::MET> > pfMETHandle;
    e.getByToken(pfMETlabel_, pfMETHandle);

    genMET_    = -99;
    genMETPhi_ = -99;

    pfMET_     = -99;
    pfMETPhi_  = -99;

    if (pfMETHandle.isValid()) 
    {
        const pat::MET *pfMET = 0;
        pfMET     = &(pfMETHandle->front());
        pfMET_    = pfMET->et();
        pfMETPhi_ = pfMET->phi();

        if (!e.isRealData()) 
        {
            genMET_    = pfMET->genMET()->et();
            genMETPhi_ = pfMET->genMET()->phi();
        }
    }
}


void photonAnalyzer::branchVertices(TTree* tree)
{
    tree->Branch("nVtx_", &nVtx_);
    tree->Branch("vtx_x", &vtx_x_);
    tree->Branch("vtx_y", &vtx_y_);
    tree->Branch("vtx_z", &vtx_z_);
    tree->Branch("vtx_xError", &vtx_xError_);
    tree->Branch("vtx_yError", &vtx_yError_);
    tree->Branch("vtx_zError", &vtx_zError_);

    tree->Branch("run",     &run_);
    tree->Branch("event",   &event_);
    tree->Branch("lumis",   &lumis_);
    tree->Branch("rho",     &rho_);
}

void photonAnalyzer::fillVertices(const edm::Event& e, const edm::EventSetup& es) 
{
    using namespace edm;
    using namespace std;
  

    vtx_x_      .clear();
    vtx_y_      .clear();
    vtx_z_      .clear();
    vtx_xError_ .clear();
    vtx_yError_ .clear();
    vtx_zError_ .clear();
    nVtx_ = 0;


    edm::Handle<std::vector<reco::Vertex>> vtxHandle;
    e.getByToken(vtxLabel_, vtxHandle);
    std::vector<reco::Vertex> vtx_sorted = *vtxHandle;
    std::sort( vtx_sorted.begin(), vtx_sorted.end(), photonAnalyzer::vtxSort );
    if(vtx_sorted.size() == 0) return;
    if(fabs(vtx_sorted.begin()->position().z()) > 15.0) return; //default
    
    vtx_x_.push_back(vtx_sorted.begin()->position().x());
    vtx_y_.push_back(vtx_sorted.begin()->position().y()); 
    vtx_z_.push_back(vtx_sorted.begin()->position().z());
    vtx_xError_.push_back(vtx_sorted.begin()->xError());
    vtx_yError_.push_back(vtx_sorted.begin()->yError());
    vtx_zError_.push_back(vtx_sorted.begin()->zError());
    nVtx_ = vtx_sorted.size();


    edm::Handle<double> rhoHandle;
    e.getByToken(rhoLabel_, rhoHandle);

    run_    = e.id().run();
    event_  = e.id().event();
    lumis_  = e.luminosityBlock();

    rho_    = *(rhoHandle.product());
}


void photonAnalyzer::branchTracks(TTree* tree)
{
    tree->Branch("nTracks", &nTracks_);
    tree->Branch("trackPt", &trackPt_);
    tree->Branch("trackEta", &trackEta_);
    tree->Branch("trackPhi", &trackPhi_);
}


void photonAnalyzer::fillTracks(const edm::Event& e, const edm::EventSetup& es) 
{
    using namespace edm;
    using namespace std;

    trackPt_    .clear();
    trackEta_   .clear();
    trackPhi_   .clear();
    nTracks_ = 0;


    edm::Handle<std::vector<reco::Vertex>> vtxHandle;
    e.getByToken(vtxLabel_, vtxHandle);
    std::vector<reco::Vertex> vtx_sorted = *vtxHandle;
    std::sort( vtx_sorted.begin(), vtx_sorted.end(), photonAnalyzer::vtxSort );
    if(vtx_sorted.size() == 0) return;
    if(fabs(vtx_sorted.begin()->position().z()) > 15.0) return; //default
    
    vtx_x = (double)vtx_sorted.begin()->position().x(); 
    vtx_y = (double)vtx_sorted.begin()->position().y(); 
    vtx_z = (double)vtx_sorted.begin()->position().z(); 
    vtx_xError = (double)vtx_sorted.begin()->xError();
    vtx_yError = (double)vtx_sorted.begin()->yError();
    vtx_zError = (double)vtx_sorted.begin()->zError();
    N_vtx = vtx_sorted.size();


    //track info
    edm::Handle<std::vector<pat::PackedCandidate> > trackCollection;
    e.getByToken(tracklabel_, trackCollection);

    math::XYZPoint vtx(vtx_x,vtx_y,vtx_z);
    for(std::vector<pat::PackedCandidate>::const_iterator itrk = trackCollection->begin(); itrk != trackCollection->end(); ++itrk)
    {
        if (itrk->hasTrackDetails())
        {
            double aux_tk_dz_vtx = (double)itrk->dz(vtx);
            double aux_tk_dzError_vtx  = (double)sqrt(itrk->dzError()*itrk->dzError()+vtx_zError*vtx_zError);
            double aux_tk_dxy_vtx = (double)itrk->dxy(vtx);
            double aux_tk_dxyError_vtx  = (double)sqrt(itrk->dxyError()*itrk->dxyError()+vtx_xError*vtx_yError);
        

            if(itrk->trackHighPurity() && fabs(aux_tk_dz_vtx/aux_tk_dzError_vtx) > 3.0 
            && fabs(aux_tk_dxy_vtx/aux_tk_dxyError_vtx) > 3.0 && itrk->pixelLayersWithMeasurement() != 0)
            {
                //trackPt->Fill(itrk->pt());
                trackPt_ .push_back(itrk->pt());
                trackEta_ .push_back(itrk->eta());
                trackPhi_ .push_back(itrk->phi());
                nTracks_++;
            }
        }
    }
}


void photonAnalyzer::branchBeamHaloMuon(TTree* tree)
{
    tree->Branch("nBHmuon",       &nBHmuon_);
    tree->Branch("BHIDbit",       &BHIDbit_);
}

void photonAnalyzer::fillBeamHaloMuon(const edm::Event& e, const edm::EventSetup& es)
{
    BHIDbit_ .clear();

    nBHmuon_ = 0;

    edm::Handle<reco::BeamHaloSummary> beamHaloSummaryHandle;
    e.getByToken(beamHaloSummaryToken_, beamHaloSummaryHandle);

    UShort_t tmpBHIDbit = 0; 

    if(beamHaloSummaryHandle.isValid())
    {
        if(beamHaloSummaryHandle->CSCLooseHaloId()) setbit(tmpBHIDbit, 0);
        if(beamHaloSummaryHandle->CSCTightHaloId()) setbit(tmpBHIDbit, 1);
        if(beamHaloSummaryHandle->CSCTightHaloId2015()) setbit(tmpBHIDbit, 2);
        if(beamHaloSummaryHandle->CSCTightHaloIdTrkMuUnveto()) setbit(tmpBHIDbit, 3);
        if(beamHaloSummaryHandle->EcalLooseHaloId()) setbit(tmpBHIDbit, 4);
        if(beamHaloSummaryHandle->EcalTightHaloId()) setbit(tmpBHIDbit, 5);
        if(beamHaloSummaryHandle->EventSmellsLikeHalo()) setbit(tmpBHIDbit, 6);
        if(beamHaloSummaryHandle->ExtremeTightId()) setbit(tmpBHIDbit, 7);
        if(beamHaloSummaryHandle->GlobalLooseHaloId()) setbit(tmpBHIDbit, 8);
        if(beamHaloSummaryHandle->GlobalSuperTightHaloId2016()) setbit(tmpBHIDbit, 9);
        if(beamHaloSummaryHandle->GlobalTightHaloId()) setbit(tmpBHIDbit, 10);
        if(beamHaloSummaryHandle->GlobalTightHaloId2016()) setbit(tmpBHIDbit, 11);
        if(beamHaloSummaryHandle->HcalLooseHaloId()) setbit(tmpBHIDbit, 12);
        if(beamHaloSummaryHandle->HcalTightHaloId()) setbit(tmpBHIDbit, 13);
        if(beamHaloSummaryHandle->LooseId()) setbit(tmpBHIDbit, 14);
        if(beamHaloSummaryHandle->TightId()) setbit(tmpBHIDbit, 15);

        BHIDbit_.push_back(tmpBHIDbit);

        nBHmuon_++;
    }
}
//
// member functions
//

// ------------ method called for each event  ------------
void
photonAnalyzer::analyze(const edm::Event& e, const edm::EventSetup& es)
{
    fillPhotons(e, es);
    fillOOTPhotons(e, es); 
    fillCrystals(e, es);
    fillTriggerinfo(e, es);
    fillMuons(e, es);
    fillElectrons(e, es);
    fillMETs(e, es);
    fillVertices(e, es);
    //fillEvents(e, es);
    fillTracks(e, es);
    //fillGenInfo(e, es);
    fillBeamHaloMuon(e,es);

    tree_->Fill();
    
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void
photonAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
photonAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
photonAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(photonAnalyzer);



/******** log *********
06232021 added beamhalomuon information for getting beamhalo from muon info
reference: https://mattermost.web.cern.ch/cms-ntgc-metg/messages/@shilpi

*/
