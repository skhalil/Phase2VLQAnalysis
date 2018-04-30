// -*- C++ -*-
//
// Package:    Upgrades/VLQAnalyzer
// Class:      VLQAnalyzer
// 
/**\class VLQAnalyzer VLQAnalyzer.cc Upgrades/VLQAnalyzer/plugins/VLQAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Sadia Khalil
//         Created:  Mon, 15 Jan 2018 19:19:16 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"

#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "Upgrades/VLQAnalyzer/interface/EventInfoTree.h"
#include "Upgrades/VLQAnalyzer/interface/METTree.h"
#include "Upgrades/VLQAnalyzer/interface/ElectronTree.h"
#include "Upgrades/VLQAnalyzer/interface/MuonTree.h"
#include "Upgrades/VLQAnalyzer/interface/JetTree.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class VLQAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
   explicit VLQAnalyzer(const edm::ParameterSet&);
   ~VLQAnalyzer();
   
   static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
   
   
private:
   virtual void beginRun(edm::Run const&, edm::EventSetup const&); //override;
   virtual void beginJob() override;
   virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
   virtual void endRun(edm::Run const&, edm::EventSetup const&); //override;
   virtual void endJob() override;
   
   bool isME0MuonSelNew(const reco::Muon&, double, double, double, edm::EventSetup const& ); //copy paste from Jan & Maria 
  
   // ----------member data ---------------------------
   edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genPartsToken_;
   edm::EDGetTokenT<GenEventInfoProduct>         genToken_;
   edm::EDGetTokenT<LHEEventProduct>             genlheToken_;
   edm::InputTag                                 puInfo_;
   edm::EDGetTokenT<std::vector<reco::Vertex>>   vtxToken_;
   edm::EDGetTokenT<std::vector<pat::Electron>>  elecsToken_;
   edm::EDGetTokenT<std::vector<pat::Muon>>      muonsToken_;
   edm::EDGetTokenT<std::vector<pat::MET>>       metsToken_;
   PFJetIDSelectionFunctor                       jetIDLoose_;
   PFJetIDSelectionFunctor                       jetIDTight_;
   edm::EDGetTokenT<edm::View<pat::Jet> >        ak4jetsToken_;
   edm::EDGetTokenT<edm::View<pat::Jet> >        ak8jetsToken_;
   edm::EDGetTokenT<edm::View<pat::Jet> >        subak8jetsToken_;
   edm::EDGetTokenT<std::vector<reco::GenJet>>   genJetsToken_;

   const ME0Geometry* ME0Geometry_; 
   bool usePuppi_;
   float ak4ptmin_, ak4etamax_, ak8ptmin_, ak8etamax_;

   edm::Service<TFileService> fs_;
   TTree* tree_;

   GenInfoTree   genevt_;  
   EventInfoTree evt_; 
   METTree       met_;
   ElectronTree  ele_;
   MuonTree      mu_;
   JetTree       ak4jet_;
   JetTree       ak8jet_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
VLQAnalyzer::VLQAnalyzer(const edm::ParameterSet& iConfig):   
   genPartsToken_  (consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
   genToken_       (consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
   genlheToken_    (consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer",""))),
   puInfo_         (iConfig.getParameter<edm::InputTag>("puInfo")),
   vtxToken_       (consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
   elecsToken_     (consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
   muonsToken_     (consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
   metsToken_      (consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
   jetIDLoose_     (PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE), 
   jetIDTight_     (PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT), 
   ak4jetsToken_   (consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
   ak8jetsToken_   (consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets_ak8"))),
   subak8jetsToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("subjets_ak8"))),
   genJetsToken_   (consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
   usePuppi_       (iConfig.getParameter<bool>("usePuppi")),
   ak4ptmin_       (iConfig.getParameter<double>("ak4ptmin")),
   ak4etamax_      (iConfig.getParameter<double>("ak4etamax")),
   ak8ptmin_       (iConfig.getParameter<double>("ak8ptmin")),
   ak8etamax_      (iConfig.getParameter<double>("ak8etamax"))
{
   consumes<std::vector<PileupSummaryInfo>>(puInfo_);
   //now do what ever initialization is needed
   usesResource("TFileService");
   
}


VLQAnalyzer::~VLQAnalyzer()
{
   
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
VLQAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   genevt_.clearTreeVectors();  
   evt_.clearTreeVectors();  
   ele_.clearTreeVectors();
   mu_.clearTreeVectors();
   ak4jet_.clearTreeVectors();
   ak8jet_.clearTreeVectors();

   // Basic event info
   evt_.runno = iEvent.eventAuxiliary().run();
   evt_.lumisec = iEvent.eventAuxiliary().luminosityBlock();
   evt_.evtno = iEvent.eventAuxiliary().event();

   //Generator weights
   float genwt = 1.0;
   edm::Handle<GenEventInfoProduct> genInfo;
   iEvent.getByToken( genToken_,genInfo);
   
   if (genInfo.isValid()){
      genwt = genInfo->weight();
      genwt /= std::abs(genwt);   
   }
   genevt_.genWt = genwt;
   
   //LHE weight
   edm::Handle<LHEEventProduct> lheInfo;
   iEvent.getByToken(genlheToken_, lheInfo);

   if(lheInfo.isValid()) {
      double asdd = lheInfo->weights()[0].wgt;
      int size = lheInfo->weights().size();
      genevt_.lheWtIDs.reserve(size);
      genevt_.lheWts.reserve(size);
      for(unsigned int i=0; i<lheInfo->weights().size(); ++i) {         
         double asdde =lheInfo->weights()[i].wgt;
         int asddeID = std::stoi(lheInfo->weights()[i].id);
         genevt_.lheWtIDs.push_back(asddeID);
         genevt_.lheWts.push_back(asdde/asdd);
      }
   }
   
   //Pileup
   edm::Handle<std::vector< PileupSummaryInfo > > puInfo;
   iEvent.getByLabel(puInfo_, puInfo);

   std::vector<PileupSummaryInfo>::const_iterator pvi;
   for(pvi = puInfo->begin(); pvi != puInfo->end(); ++pvi) {
      evt_.npuTrue = pvi->getTrueNumInteractions(); 
      evt_.npuInt = pvi->getBunchCrossing(); 
      evt_.puBX = pvi->getPU_NumInteractions();
   }
  
   edm::Handle<reco::VertexCollection> vertices;
   iEvent.getByToken(vtxToken_, vertices);

   if (vertices->empty()) return;
   evt_.npv=0; evt_.nGoodVtx = 0;
   evt_.npv = vertices->size();
   
   int prVtx = -1;
   for (size_t i = 0; i < vertices->size(); i++) {
      if (vertices->at(i).isFake()) continue;
      if (vertices->at(i).ndof() <= 4) continue; 
      prVtx = i;  
      evt_.vPt2.push_back(vertices->at(i).p4().pt());
      evt_.nGoodVtx++;
   }
   auto primaryVertex=vertices->at(prVtx);
   
   // MET
   edm::Handle<std::vector<pat::MET>> mets;
   iEvent.getByToken(metsToken_, mets);
   met_.px=0.;  met_.py=0.; met_.pt=0.; met_.eta=0.; met_.phi=0.;
   if (mets->size() != 0 ){
         met_.px = mets->at(0).px();
         met_.py = mets->at(0).py();
         met_.pt = mets->at(0).pt();
         met_.eta = mets->at(0).eta();
         met_.phi = mets->at(0).phi();
      }

   //Electrons
   edm::Handle<std::vector<pat::Electron>> elecs;
   iEvent.getByToken(elecsToken_, elecs);
    
   //ele_.nLoose=0; ele_.nMedium=0.; ele_.nTight=0.;
   for (const pat::Electron & i : *elecs) {
      if (i.pt() < 20.) continue;
      if (fabs(i.eta()) > 4.) continue;
      float mvaValue = i.userFloat("mvaValue");
      bool isEB = i.isEB();
      bool isLoose (0), isMedium (0), isTight (0);
      
       if( isEB ) {
          if (i.pt() < 20.) {
            isLoose  = (mvaValue > -0.661);
            isMedium = (mvaValue > 0.855);
            isTight  = (mvaValue > 0.986);
         }
         else {
            isLoose  = (mvaValue > -0.797);
            isMedium = (mvaValue > 0.723);
            isTight  = (mvaValue > 0.988);
         }
      }
      else {
         if (not (i.userFloat("hgcElectronID:ecEnergy") > 0)) continue;
         if (not (i.userFloat("hgcElectronID:sigmaUU") > 0)) continue;
         if (not (i.fbrem() > -1)) continue;
         if (not (i.userFloat("hgcElectronID:measuredDepth") < 40)) continue;
         if (not (i.userFloat("hgcElectronID:nLayers") > 20)) continue;
         if (i.pt() < 20.) {
            isLoose = (mvaValue > -0.320);
            isMedium = mvaValue > 0.777;
            isTight = (mvaValue > 0.969);
         }
         else {
            isLoose = (mvaValue > -0.919);
            isMedium = mvaValue > 0.591;
            isTight = (mvaValue > 0.983);
         }
      }

      float dxy=0;
      float dz=0;
      if(i.gsfTrack().isNonnull()){
         dxy=std::abs(i.gsfTrack()->dxy(primaryVertex.position()));
         dz=std::abs(i.gsfTrack()->dz(primaryVertex.position()));
      }
      float iso =   i.puppiNoLeptonsChargedHadronIso() + i.puppiNoLeptonsNeutralHadronIso() + i.puppiNoLeptonsPhotonIso();

      if (!isLoose) continue;
      
      int eleType = 1*isLoose + 2*isMedium + 4*isTight;
      ele_.eleWP.   push_back(eleType);
      
      //if ( (eleType & 1) == 1) {ele_.nLoose++ ;}
      //if ( (eleType & 2) == 2) {ele_.nMedium++ ;}
      //if ( (eleType & 4) == 4) {ele_.nTight++ ;}
          
      ele_.pt.     push_back(i.pt()) ;
      ele_.eta.    push_back(i.eta());
      ele_.phi.    push_back(i.phi());
      ele_.mass.   push_back(i.mass());
      ele_.energy. push_back(i.energy());
      ele_.charge. push_back(i.charge());
      ele_.dz.     push_back(dz);
      ele_.dxy.    push_back(dxy); 
      ele_.mva.    push_back(mvaValue);
      ele_.relIso. push_back(iso/i.pt());
   }

   //Muons
   // muon isolation comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2MuonBarrelRecipes#Muon_isolation
   // muon ID comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2MuonBarrelRecipes#Muon_identification
   Handle<std::vector<pat::Muon>> mus;
   iEvent.getByToken(muonsToken_, mus);
   
   //mu_.nLoose=0; mu_.nMedium=0.; mu_.nTight=0.;
   for (const pat::Muon & i : *mus) {
      if (i.pt() < 20.) continue;
      if (fabs(i.eta()) > 4.) continue;
      
      bool ipxy = false, ipz = false, validPxlHit = false, highPurity = false;
      float dz=0, dxy=0;
      if (i.innerTrack().isNonnull()){
         dxy=std::abs(i.muonBestTrack()->dxy(vertices->at(prVtx).position()));
         dz= std::abs(i.muonBestTrack()->dz(vertices->at(prVtx).position()));
         ipxy = dxy < 0.2;
         ipz = dz < 0.5;
         validPxlHit = i.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
         highPurity = i.innerTrack()->quality(reco::Track::highPurity);
      } 
      float iso = i.trackIso(); 

      //Loose
      double dPhiCut = std::min(std::max(1.2/i.p(),1.2/100),0.056);
      double dPhiBendCut = std::min(std::max(0.2/i.p(),0.2/100),0.0096);       
      bool isLoose = (fabs(i.eta()) < 2.4 && muon::isLooseMuon(i)) || 
         (fabs(i.eta()) > 2.4 && isME0MuonSelNew(i, 0.077, dPhiCut, dPhiBendCut,iSetup));
      
      // Medium ID -- needs to be updated
      bool isMedium = (fabs(i.eta()) < 2.4 && muon::isMediumMuon(i)) || 
         (fabs(i.eta()) > 2.4 && isME0MuonSelNew(i, 0.077, dPhiCut, dPhiBendCut, iSetup) && ipxy && ipz && validPxlHit && highPurity);
      
      // Tight ID
      dPhiCut = std::min(std::max(1.2/i.p(),1.2/100),0.032);
      dPhiBendCut = std::min(std::max(0.2/i.p(),0.2/100),0.0041);
      bool isTight = (fabs(i.eta()) < 2.4 && vertices->size() > 0 && muon::isTightMuon(i,vertices->at(prVtx))) || 
         (fabs(i.eta()) > 2.4 && isME0MuonSelNew(i, 0.048, dPhiCut, dPhiBendCut,iSetup) && ipxy && ipz && validPxlHit && highPurity);
      
      //
      int muType = 1*isLoose + 2*isMedium + 4*isTight;
      mu_.muWP.   push_back(muType);
      mu_.charge. push_back(i.charge());
      mu_.pt.     push_back(i.pt());
      mu_.phi.    push_back(i.phi());
      mu_.eta.    push_back(i.eta());
      mu_.mass.   push_back(i.mass());
      mu_.dxy.    push_back(dxy);
      mu_.dz.     push_back(dz);
      mu_.relIso. push_back(iso/i.pt());
   }

   //Gen Jets
   std::vector<reco::GenJet> cleanGenJets;
   
   edm::Handle<std::vector<pat::PackedGenParticle>> genParts;
   iEvent.getByToken(genPartsToken_, genParts);
   
   //Jets
   edm::Handle<edm::View<pat::Jet> > ak4jets;
   iEvent.getByToken(ak4jetsToken_, ak4jets);

   for (const pat::Jet & j : *ak4jets) {
      if (j.pt() < ak4ptmin_ || abs(j.eta()) > ak4etamax_) continue; 
      
      const reco::GenJet* genjet = j.genJet() ;
      bool overlaps = false;
      
      for (const pat::PackedGenParticle & gp : *genParts) {
         if (abs(gp.pdgId()) != 11 && abs(gp.pdgId()) != 13) continue;
         if (genjet != nullptr && 
             fabs(j.genJet()->pt()-gp.pt()) < 0.01*gp.pt() && 
             ROOT::Math::VectorUtil::DeltaR(gp.p4(),j.genJet()->p4()) < 0.01) 
         {overlaps = true; break;}
      }
      
      if (genjet != nullptr && !overlaps) {
         ak4jet_.genjetpt.push_back(j.genJet()->pt());
         ak4jet_.genjeteta.push_back(j.genJet()->eta());
         ak4jet_.genjetphi.push_back(j.genJet()->phi());
         ak4jet_.genjetenergy.push_back(j.genJet()->energy());
         ak4jet_.genjetmass.push_back(j.genJet()->mass()) ;
         cleanGenJets.push_back(*genjet);
      }
      else {
         ak4jet_.genjetpt.push_back(-9999)    ;
         ak4jet_.genjeteta.push_back(-9999)   ;
         ak4jet_.genjetphi.push_back(-9999)   ;
         ak4jet_.genjetenergy.push_back(-9999);  
         ak4jet_.genjetmass.push_back(-9999)  ;
      }

      ak4jet_.pt           .push_back(j.pt())  ;
      ak4jet_.eta          .push_back(j.eta()) ;
      ak4jet_.phi          .push_back(j.phi()) ;
      ak4jet_.energy       .push_back(j.energy());
      ak4jet_.mass         .push_back(j.mass());
      ak4jet_.partonFlavour.push_back(j.partonFlavour()); 
      ak4jet_.hadronFlavour.push_back(j.hadronFlavour()); 
      ak4jet_.csvv2        .push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      ak4jet_.deepcsv      .push_back(j.bDiscriminator("pfDeepCSVJetTags:probb") 
                                       + j.bDiscriminator("pfDeepCSVJetTags:probbb"));
      if (!usePuppi_){
         ak4jet_.pujetid   .push_back(j.userFloat("pileupJetId:fullDiscriminant")); 
      }
   } 
   
   edm::Handle<edm::View<pat::Jet> > ak8jets;
   iEvent.getByToken(ak8jetsToken_, ak8jets);
   for (const pat::Jet & j : *ak8jets) {
      if (j.pt() < ak8ptmin_ || abs(j.eta()) > ak8etamax_) continue; 
      const reco::GenJet* genjet = j.genJet() ; 
      if (genjet != nullptr) {
         ak8jet_.genjetpt.push_back(j.genJet()->pt())  ;
      }
      else {
         ak8jet_.genjetpt.push_back(-9999)  ;
      }
      ak8jet_.pt               .push_back(j.pt())  ;
      ak8jet_.eta              .push_back(j.eta()) ;
      ak8jet_.phi              .push_back(j.phi()) ;
      ak8jet_.energy           .push_back(j.energy());
      ak8jet_.mass             .push_back(j.mass());
      ak8jet_.ptCHS            .push_back(j.userFloat("ak8PFJetsCHSValueMap:pt"))  ;
      ak8jet_.etaCHS           .push_back(j.userFloat("ak8PFJetsCHSValueMap:eta"))  ;
      ak8jet_.phiCHS           .push_back(j.userFloat("ak8PFJetsCHSValueMap:phi"))  ;
      ak8jet_.massCHS          .push_back(j.userFloat("ak8PFJetsCHSValueMap:mass"))  ;
      ak8jet_.softDropMassCHS  .push_back(j.userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass")) ;  
      ak8jet_.prunedMassCHS    .push_back(j.userFloat("ak8PFJetsCHSValueMap:ak8PFJetsCHSPrunedMass")) ;  
      ak8jet_.tau1CHS          .push_back(j.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau1"));
      ak8jet_.tau2CHS          .push_back(j.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau2")); 
      ak8jet_.tau3CHS          .push_back(j.userFloat("ak8PFJetsCHSValueMap:NjettinessAK8CHSTau3"));
      ak8jet_.softDropMassPuppi.push_back(j.userFloat("ak8PFJetsPuppiSoftDropMass")) ;  
      ak8jet_.tau1Puppi        .push_back(j.userFloat("NjettinessAK8Puppi:tau1"));
      ak8jet_.tau2Puppi        .push_back(j.userFloat("NjettinessAK8Puppi:tau2")); 
      ak8jet_.tau3Puppi        .push_back(j.userFloat("NjettinessAK8Puppi:tau3"));
      ak8jet_.csvv2            .push_back(j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      ak8jet_.deepcsv          .push_back(j.bDiscriminator("pfDeepCSVJetTags:probb") 
                                           + j.bDiscriminator("pfDeepCSVJetTags:probbb"));
      ak8jet_.partonFlavour    .push_back(j.partonFlavour()); 
      ak8jet_.hadronFlavour    .push_back(j.hadronFlavour()); 
      
      std::vector<edm::Ptr<pat::Jet> > const& sdsubjets = j.subjets("SoftDropPuppi") ;
      if (sdsubjets.size() < 2) continue ;
      ak8jet_.sj0pt           .push_back(sdsubjets.at(0)->pt()) ; 
      ak8jet_.sj1pt           .push_back(sdsubjets.at(1)->pt()) ; 
      ak8jet_.sj0eta          .push_back(sdsubjets.at(0)->eta()) ; 
      ak8jet_.sj1eta          .push_back(sdsubjets.at(1)->eta()) ; 
      ak8jet_.sj0phi          .push_back(sdsubjets.at(0)->phi()) ; 
      ak8jet_.sj1phi          .push_back(sdsubjets.at(1)->phi()) ; 
      ak8jet_.sj0partonFlavour.push_back(sdsubjets.at(0)->partonFlavour()); 
      ak8jet_.sj0hadronFlavour.push_back(sdsubjets.at(0)->hadronFlavour()); 
      ak8jet_.sj1partonFlavour.push_back(sdsubjets.at(1)->partonFlavour()); 
      ak8jet_.sj1hadronFlavour.push_back(sdsubjets.at(1)->hadronFlavour()); 
      ak8jet_.sj0csvv2        .push_back(sdsubjets.at(0)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      ak8jet_.sj1csvv2        .push_back(sdsubjets.at(1)->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      ak8jet_.sj0deepcsv      .push_back(sdsubjets.at(0)->bDiscriminator("pfDeepCSVJetTags:probb") 
                                          + sdsubjets.at(0)->bDiscriminator("pfDeepCSVJetTags:probbb"));
      ak8jet_.sj1deepcsv      .push_back(sdsubjets.at(1)->bDiscriminator("pfDeepCSVJetTags:probb") 
                                          + sdsubjets.at(1)->bDiscriminator("pfDeepCSVJetTags:probbb"));
   }
   
   //std::cout << "to this end : " << std::endl;
   tree_->Fill();
}

// ------------ private functions --------------
bool
VLQAnalyzer::isME0MuonSelNew(const reco::Muon& muon, double dEtaCut, double dPhiCut, double dPhiBendCut, edm::EventSetup const& iSetup)
{
   bool result = false;
   bool isME0 = muon.isME0Muon();

  if(isME0){
     double deltaEta = 999;
     double deltaPhi = 999;
     double deltaPhiBend = 999;
     
     if(!ME0Geometry_){
       edm::ESHandle<ME0Geometry> hGeom;
       iSetup.get<MuonGeometryRecord>().get(hGeom);
       ME0Geometry_ =( &*hGeom);
       if(!ME0Geometry_)
          return false;
    }

     const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
     for( std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber ){       
        if (chamber->detector() == 5){
           
           for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); 
                 segment != chamber->me0Matches.end(); ++segment ){

             LocalPoint trk_loc_coord(chamber->x, chamber->y, 0);
             LocalPoint seg_loc_coord(segment->x, segment->y, 0);
             LocalVector trk_loc_vec(chamber->dXdZ, chamber->dYdZ, 1);
             LocalVector seg_loc_vec(segment->dXdZ, segment->dYdZ, 1);
             
             const ME0Chamber * me0chamber = ME0Geometry_->chamber(chamber->id);
             if(!me0chamber)continue;
             
             GlobalPoint trk_glb_coord = me0chamber->toGlobal(trk_loc_coord);
             GlobalPoint seg_glb_coord = me0chamber->toGlobal(seg_loc_coord);
             
             //double segDPhi = segment->me0SegmentRef->deltaPhi();
             // need to check if this works
             double segDPhi = me0chamber->computeDeltaPhi(seg_loc_coord, seg_loc_vec);
             double trackDPhi = me0chamber->computeDeltaPhi(trk_loc_coord, trk_loc_vec);
             
             deltaEta = std::abs(trk_glb_coord.eta() - seg_glb_coord.eta() );
             deltaPhi = std::abs(trk_glb_coord.phi() - seg_glb_coord.phi() );
             deltaPhiBend = std::abs(segDPhi - trackDPhi);
             
             if (deltaEta < dEtaCut && deltaPhi < dPhiCut && deltaPhiBend < dPhiBendCut) result = true;
          }
       }
    }//for loop
  }//if
  return result;
}

// ------------ method called once each run ----------------
void
VLQAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
  edm::ESHandle<ME0Geometry> hGeom;
  iSetup.get<MuonGeometryRecord>().get(hGeom);
  ME0Geometry_ =( &*hGeom);
}

// ------------ method called when ending the processing of a run  ------------
  void
VLQAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called once each job just before starting event loop  ------------
void 
VLQAnalyzer::beginJob()
{
  tree_ = fs_->make<TTree>("anatree", "anatree") ;
  evt_.RegisterTree(tree_, "SelectedEvt") ;
  genevt_.RegisterTree(tree_, "GenEvt") ;
  ele_.RegisterTree(tree_, "Electons") ;
  met_.RegisterTree(tree_, "MET");
  ak4jet_.RegisterTree(tree_, "AK4JetsCHS") ; 
  ak8jet_.RegisterTree(tree_, "AK8Jets") ; 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VLQAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VLQAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(VLQAnalyzer);
