#include <string>
#include <vector>
#include <TTree.h>
#include <TLorentzVector.h>

class JetTree {
public:
   std::vector<int>  idx; 
   std::vector<double> genjetpt;
   std::vector<double> pt;
   std::vector<double> eta;
   std::vector<double> phi;
   std::vector<double> energy;
   std::vector<double> mass;
   std::vector<double> ptCHS;
   std::vector<double> etaCHS;
   std::vector<double> phiCHS;
   std::vector<double> massCHS;
   std::vector<double> softDropMassCHS;
   std::vector<double> prunedMassCHS;
   std::vector<double> tau1CHS;
   std::vector<double> tau2CHS;
   std::vector<double> tau3CHS;
   std::vector<double> softDropMassPuppi;
   std::vector<double> tau1Puppi;
   std::vector<double> tau2Puppi;
   std::vector<double> tau3Puppi;
   std::vector<double> csvv2;
   std::vector<double> deepcsv;
   std::vector<double> pujetid;
   std::vector<double> partonFlavour; 
   std::vector<double> hadronFlavour; 
   std::vector<double> sj0csvv2;
   std::vector<double> sj1csvv2;
   std::vector<double> sj0deepcsv;
   std::vector<double> sj1deepcsv;
   std::vector<double> sj0partonFlavour; 
   std::vector<double> sj1partonFlavour; 
   std::vector<double> sj0hadronFlavour; 
   std::vector<double> sj1hadronFlavour; 
   std::vector<double> sj0pt;
   std::vector<double> sj1pt;
   std::vector<double> sj0eta;
   std::vector<double> sj1eta;
   std::vector<double> sj0phi;
   std::vector<double> sj1phi;
  
   void clearTreeVectors() {
      idx               .clear() ;     
      genjetpt          .clear() ;    
      pt                .clear() ;    
      eta               .clear() ;    
      phi               .clear() ;    
      energy            .clear() ;    
      mass              .clear() ;    
      ptCHS             .clear() ;    
      etaCHS            .clear() ;    
      phiCHS            .clear() ;    
      massCHS           .clear() ;    
      softDropMassCHS   .clear() ;    
      prunedMassCHS     .clear() ;    
      tau1CHS           .clear() ;    
      tau2CHS           .clear() ;    
      tau3CHS           .clear() ;    
      softDropMassPuppi .clear() ;    
      tau1Puppi         .clear() ;    
      tau2Puppi         .clear() ;    
      tau3Puppi         .clear() ;    
      csvv2             .clear() ;    
      deepcsv           .clear() ;    
      pujetid           .clear() ;    
      partonFlavour     .clear() ;     
      hadronFlavour     .clear() ;     
      sj0csvv2          .clear() ;    
      sj1csvv2          .clear() ;    
      sj0deepcsv        .clear() ;    
      sj1deepcsv        .clear() ;    
      sj0partonFlavour  .clear() ;     
      sj1partonFlavour  .clear() ;     
      sj0hadronFlavour  .clear() ;     
      sj1hadronFlavour  .clear() ;     
      sj0pt             .clear() ;    
      sj1pt             .clear() ;    
      sj0eta            .clear() ;    
      sj1eta            .clear() ;    
      sj0phi            .clear() ;    
      sj1phi            .clear() ;    
    }

   void RegisterTree(TTree* tree, std::string name="Jets") {
      tree->Branch((name+"_idx").c_str(), &idx); 
      tree->Branch((name+"_genjetpt").c_str(), &genjetpt); 
      tree->Branch((name+"_pt").c_str(), &pt); 
      tree->Branch((name+"_eta").c_str(), &eta);
      tree->Branch((name+"_phi").c_str(), &phi);
      tree->Branch((name+"_energy").c_str(), &energy);
      tree->Branch((name+"_mass").c_str(), &mass);
      tree->Branch((name+"_csvv2").c_str(), &csvv2);
      tree->Branch((name+"_deepcsv").c_str(), &deepcsv);
      tree->Branch((name+"_pujetid").c_str(), &pujetid);
      
      if ( name.find("AK8Jets") != std::string::npos ) {
        tree->Branch((name+"_ptCHS").c_str(), &ptCHS);
        tree->Branch((name+"_etaCHS").c_str(), &etaCHS);
        tree->Branch((name+"_phiCHS").c_str(), &phiCHS);
        tree->Branch((name+"_massCHS").c_str(), &massCHS);
        tree->Branch((name+"_softDropMassCHS").c_str(), &softDropMassCHS);
        tree->Branch((name+"_prunedCHS").c_str(), &prunedMassCHS);
        tree->Branch((name+"_tau1CHS").c_str(), &tau1CHS);
        tree->Branch((name+"_tau2CHS").c_str(), &tau2CHS);
        tree->Branch((name+"_tau3CHS").c_str(), &tau3CHS);
        tree->Branch((name+"_softDropMassPuppi").c_str(), &softDropMassPuppi);
        tree->Branch((name+"_tau1Puppi").c_str(), &tau1Puppi);
        tree->Branch((name+"_tau2Puppi").c_str(), &tau2Puppi);
        tree->Branch((name+"_tau3Puppi").c_str(), &tau3Puppi);
        tree->Branch((name+"_partonFlavour").c_str(), &partonFlavour);
        tree->Branch((name+"_hadronFlavour").c_str(), &hadronFlavour);
        tree->Branch((name+"_sj0csvv2").c_str(),&sj0csvv2);
        tree->Branch((name+"_sj1csvv2").c_str(),&sj1csvv2);
        tree->Branch((name+"_sj0deepcsv").c_str(),&sj0deepcsv);
        tree->Branch((name+"_sj1deepcsv").c_str(),&sj1deepcsv);
        tree->Branch((name+"_sj0partonflavour").c_str(),&sj0partonFlavour);
        tree->Branch((name+"_sj1partonflavour").c_str(),&sj1partonFlavour);
        tree->Branch((name+"_sj0hadronflavour").c_str(),&sj0hadronFlavour);
        tree->Branch((name+"_sj1hadronflavour").c_str(),&sj1hadronFlavour);
        tree->Branch((name+"_sj0pt").c_str(),&sj0pt);
        tree->Branch((name+"_sj1pt").c_str(),&sj1pt);
        tree->Branch((name+"_sj0eta").c_str(),&sj0eta);
        tree->Branch((name+"_sj1eta").c_str(),&sj1eta);
        tree->Branch((name+"_sj0phi").c_str(),&sj0phi);
        tree->Branch((name+"_sj1phi").c_str(),&sj1phi);
      }
      
   }

};
