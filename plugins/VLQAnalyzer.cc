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

#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include "Upgrades/VLQAnalyzer/interface/EventInfoTree.h"
#include "Upgrades/VLQAnalyzer/interface/MET.h"
#include "Upgrades/VLQAnalyzer/interface/ElectronTree.h"
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
   virtual void beginJob() override;
   virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
   virtual void endJob() override;
   
   // ----------member data ---------------------------
   edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genPartsToken_;
   edm::EDGetTokenT<GenEventInfoProduct> genToken_;
   edm::EDGetTokenT<LHEEventProduct>     genlheToken_;
   //unsigned int pileup_;
   //edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puInfo_;
   edm::InputTag puInfo_;
   //edm::EDGetTokenT<reco::VertexCollection >   vtxToken_;
   edm::EDGetTokenT<std::vector<reco::Vertex>>     vtxToken_;
   edm::EDGetTokenT<std::vector<pat::Electron>>    elecsToken_;
   edm::EDGetTokenT<reco::BeamSpot>                bsToken_;
   edm::EDGetTokenT<std::vector<reco::Conversion>> convToken_;
   
   edm::EDGetTokenT<std::vector<pat::MET>> metsToken_;
   PFJetIDSelectionFunctor                 jetIDLoose_;
   PFJetIDSelectionFunctor                 jetIDTight_;
   edm::EDGetTokenT<edm::View<pat::Jet> >  ak4jetsToken_;
   edm::EDGetTokenT<edm::View<pat::Jet> >  ak8jetsToken_;
   edm::EDGetTokenT<edm::View<pat::Jet> >  subak8jetsToken_;
   //edm::EDGetTokenT<reco::GenParticleCollection> genparToken_;
   edm::EDGetTokenT<std::vector<reco::GenJet>> genJetsToken_;
   float ak4ptmin_, ak4etamax_, ak8ptmin_, ak8etamax_;

   edm::Service<TFileService> fs_;
   TTree* tree_;

   GenInfoTree   genevt_;  
   EventInfoTree evt_; 
   METTree       met_;
   ElectronTree  ele_;
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
   //pileup_         (iConfig.getParameter<unsigned int>("pileup"))
   //vtxToken_       (consumes<reco::VertexCollection> (iConfig.getParameter<edm::InputTag>("vertices"))),
   vtxToken_       (consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
   elecsToken_     (consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
   bsToken_        (consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
   convToken_      (consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions"))),
   metsToken_      (consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
   jetIDLoose_     (PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE), 
   jetIDTight_     (PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT), 
   ak4jetsToken_   (consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets"))),
   ak8jetsToken_   (consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("jets_ak8"))),
   subak8jetsToken_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("subjets_ak8"))),
   genJetsToken_   (consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
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
   evt_.npv = vertices->size();
   evt_.nGoodVtx = 0;
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
   //Handle<reco::ConversionCollection> conversions;
   //iEvent.getByToken(convToken_, conversions);
   //Handle<reco::BeamSpot> bsHandle;
   //iEvent.getByToken(bsToken_, bsHandle);
   //const reco::BeamSpot &beamspot = *bsHandle.product(); 

   //int nLooseEle(0), nMediumEle(0), nTightEle(0);
   ele_.nLoose=0; ele_.nMedium=0.; ele_.nTight=0.;
   for (const pat::Electron & i : *elecs) {

      if (i.pt() < 10.) continue;
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
      
      if ( (eleType & 1) == 1) {ele_.nLoose++ ;}
      if ( (eleType & 2) == 4) {ele_.nMedium++ ;}
      if ( (eleType & 4) == 4) {ele_.nTight++ ;}
      //nLooseEle++;
     
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
/*
      if (!isMedium) continue;
      nMediumEle++;

      ele_.ptM.     push_back(i.pt()) ;
      ele_.etaM.    push_back(i.eta());
      ele_.phiM.    push_back(i.phi());
      ele_.massM.   push_back(i.mass());
      ele_.energyM. push_back(i.energy());
      ele_.chargeM. push_back(i.charge());
      ele_.dzM.     push_back(dz);
      ele_.dxyM.    push_back(dxy); 
      ele_.mvaM.    push_back(mvaValue);
      ele_.relIsoM. push_back(iso/i.pt());

      if (!isTight) continue;
      nTightEle++;   

      ele_.ptT.     push_back(i.pt()) ;
      ele_.etaT.    push_back(i.eta());
      ele_.phiT.    push_back(i.phi());
      ele_.massT.   push_back(i.mass());
      ele_.energyT. push_back(i.energy());
      ele_.chargeT. push_back(i.charge());
      ele_.dzT.     push_back(dz);
      ele_.dxyT.    push_back(dxy); 
      ele_.mvaT.    push_back(mvaValue);
      ele_.relIsoT. push_back(iso/i.pt());   
   }
   ele_.nL.     push_back(nLooseEle);
   ele_.nM.     push_back(nMediumEle);
   ele_.nT.     push_back(nTightEle);
*/
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
      
      if (genjet != nullptr && !overlaps) {// || !overlaps
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
         //std::cout << ", overlaps : " << overlaps << std::endl;
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
      ak4jet_.pujetid      .push_back(j.userFloat("pileupJetId:fullDiscriminant")); 
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
