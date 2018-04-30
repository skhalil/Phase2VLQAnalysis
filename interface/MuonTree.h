#include <string>
#include <vector>
#include <TTree.h>
#include <TLorentzVector.h>

class MuonTree {
public:
   //int nLoose, nMedium, nTight;
   std::vector<int>   muWP;   
   std::vector<float> pt; 
   std::vector<float> eta;
   std::vector<float> phi;
   std::vector<float> mass;
   std::vector<float> energy; 
   std::vector<float> charge; 
   std::vector<float> dz; 
   std::vector<float> dxy;
   //std::vector<float> mva;
   std::vector<float> relIso; 
   
   void clearTreeVectors() {
      muWP. clear();
      pt.    clear();
      eta.   clear();
      phi.   clear();
      mass.  clear();
      energy.clear();
      charge.clear();
      dz.    clear();
      dxy.   clear();
      //mva.   clear();
      relIso.clear();
   }

   void RegisterTree(TTree* tree, std::string name="Muons") {
      //tree->Branch((name+"_nLoose") .c_str(), &nLoose);
      //tree->Branch((name+"_nMedium").c_str(), &nMedium);   
      //tree->Branch((name+"_nTight") .c_str(), &nTight);
      tree->Branch((name+"_muWP")   .c_str(), &muWP);
      tree->Branch((name+"_pt")     .c_str(), &pt);
      tree->Branch((name+"_eta")    .c_str(), &eta); 
      tree->Branch((name+"_phi")    .c_str(), &phi);
      tree->Branch((name+"_mass")   .c_str(), &mass);
      tree->Branch((name+"_energy") .c_str(), &energy); 
      tree->Branch((name+"_charge") .c_str(), &charge);
      tree->Branch((name+"_dz")     .c_str(), &dz);
      tree->Branch((name+"_dxy")    .c_str(), &dxy); 
      //tree->Branch((name+"_mva")    .c_str(), &mva);
      tree->Branch((name+"_relIso") .c_str(), &relIso); 
   }
};
