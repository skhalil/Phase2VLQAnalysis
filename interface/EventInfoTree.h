#include <string>
#include <vector>
#include <TTree.h>
#include <TLorentzVector.h>

class EventInfoTree {
  public:

    int    runno ;
    int    lumisec ;
    int    evtno;
    int    npv;   
    int    npuTrue;
    int    puBX;
    int    npuInt;
    std::vector<double> vPt2;
    int    nvtx;

    void clearTreeVectors() {
         vPt2.clear();
    }

    void RegisterTree(TTree* tree, std::string name="SelectedEvents") {
      tree->Branch((name+"_runno").c_str(), &runno, (name+"_runno/I").c_str());
      tree->Branch((name+"_lumisec").c_str(), &lumisec, (name+"_lumisec/I").c_str());
      tree->Branch((name+"_evtno").c_str(), &evtno, (name+"_evtno/I").c_str());
      tree->Branch((name+"_npv").c_str(), &npv, "npv/I");
      tree->Branch((name+"_npuTrue").c_str(), &npuTrue, "npuTrue/I");
      tree->Branch((name+"_puBX").c_str(), &puBX, "puBX/I");
      tree->Branch((name+"_npuInt").c_str(), &npuInt, "npuInt/I");
      tree->Branch((name+"_vPt2").c_str(), &vPt2, (name+"_vPt2/D").c_str());
      tree->Branch((name+"_nvtx").c_str(), &nvtx, (name+"_nvtx/D").c_str());
    }
};
