#include <string>
#include <vector>
#include <TTree.h>
#include <TLorentzVector.h>

class METTree {
  public:

   float  px, py, pt, eta, phi;
   
    void RegisterTree(TTree* tree, std::string name="MET") {
      tree->Branch((name+"_px").c_str(), &py, (name+"_px/F").c_str());
      tree->Branch((name+"_py").c_str(), &py, (name+"_py/F").c_str());
      tree->Branch((name+"_pt").c_str(), &pt, (name+"_pt/F").c_str());
      tree->Branch((name+"_eta").c_str(), &eta, (name+"_eta/F").c_str());
      tree->Branch((name+"_phi").c_str(), &phi, (name+"_phi/F").c_str());
    }
};
