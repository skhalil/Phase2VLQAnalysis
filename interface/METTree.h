#include <string>
#include <vector>
#include <TTree.h>
#include <TLorentzVector.h>

class METTree {
  public:

   float  px, py, pz, et, pt, phi;
   
    void RegisterTree(TTree* tree, std::string name="MET") {
      tree->Branch((name+"_px").c_str(), &py, (name+"_px/F").c_str());
      tree->Branch((name+"_py").c_str(), &py, (name+"_py/F").c_str());
      tree->Branch((name+"_pz").c_str(), &pz, (name+"_pz/F").c_str());
      tree->Branch((name+"_et").c_str(), &et, (name+"_et/F").c_str());
      tree->Branch((name+"_pt").c_str(), &pt, (name+"_pt/F").c_str());
      tree->Branch((name+"_phi").c_str(), &phi, (name+"_phi/F").c_str());
    }
};
