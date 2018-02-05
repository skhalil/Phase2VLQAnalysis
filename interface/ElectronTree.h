#include <string>
#include <vector>
#include <TTree.h>
#include <TLorentzVector.h>

class ElectronTree {
public:
   int   size ;  
   
   void RegisterTree(TTree* tree, std::string name="ElectronLoose") {
   tree->Branch((name+"_size").c_str(), &size, (name+"_size/I").c_str());
      
   }
};
