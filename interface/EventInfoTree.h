#include <string>
#include <vector>
#include <TTree.h>
#include <TLorentzVector.h>

class EventInfoTree {
  public:
    int    countEvents ;
   //int    runno ;
    int    lumisec ;
    int    evtno;
    int    npv;   
    int    npuTrue;
    int    puBX;
    int    npuInt;
    int    nGoodVtx;
    std::vector<int>   ndofVtx;
    std::vector<float> chi2Vtx;
    std::vector<float> zVtx;
    std::vector<float> rhoVtx;
    std::vector<float> ptVtx;

    void clearTreeVectors() {
       ndofVtx.clear();
       chi2Vtx.clear();
       zVtx.clear();
       rhoVtx.clear();
       ptVtx.clear();
    }

    void RegisterTree(TTree* tree, std::string name="SelectedEvents") {
      tree->Branch((name+"_countEvents").c_str(), &countEvents, (name+"_countEvents/I").c_str());
      //tree->Branch((name+"_runno").c_str(), &runno, (name+"_runno/I").c_str());
      tree->Branch((name+"_lumisec").c_str(), &lumisec, (name+"_lumisec/I").c_str());
      tree->Branch((name+"_evtno").c_str(), &evtno, (name+"_evtno/I").c_str());
      tree->Branch((name+"_npv").c_str(), &npv, "npv/I");
      tree->Branch((name+"_npuTrue").c_str(), &npuTrue, "npuTrue/I");
      tree->Branch((name+"_puBX").c_str(), &puBX, "puBX/I");
      tree->Branch((name+"_npuInt").c_str(), &npuInt, "npuInt/I");      
      tree->Branch((name+"_nGoodVtx").c_str(), &nGoodVtx, (name+"_nGoodVtx/I").c_str());
      tree->Branch((name+"_ndofVtx").c_str(), &ndofVtx, (name+"_ndofVtx/I").c_str());
      tree->Branch((name+"_chi2Vtx").c_str(), &chi2Vtx, (name+"_chi2Vtx/F").c_str());
      tree->Branch((name+"_zVtx").c_str(), &zVtx, (name+"_zVtx/F").c_str());
      tree->Branch((name+"_rhoVtx").c_str(), &rhoVtx, (name+"_rhoVtx/F").c_str()); 
      tree->Branch((name+"_ptVtx").c_str(), &ptVtx, (name+"_ptVtx/F").c_str());
    }
};

class GenInfoTree {
public:
   float genWt;
   std::vector<float> lheWts;
   std::vector<int> lheWtIDs;
   void clearTreeVectors() {
      lheWtIDs.clear();
      lheWts.clear();
   }
   void RegisterTree(TTree* tree, std::string name="EventWeights"){
    tree->Branch((name+"_genWt").c_str(), &genWt, (name+"_genWt/F").c_str());  
    tree->Branch((name+"_lheWts").c_str(), &lheWts, (name+"_lheWts/F").c_str());
    tree->Branch((name+"_lheWtIDs").c_str(), &lheWtIDs, (name+"_lheWtIDs/F").c_str());
   } 

};

