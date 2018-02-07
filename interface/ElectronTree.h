#include <string>
#include <vector>
#include <TTree.h>
#include <TLorentzVector.h>

class ElectronTree {
public:
   std::vector<int>   nL, nM, nT;  
   std::vector<float> ptL, ptM, ptT;
   std::vector<float> etaL, etaM, etaT;
   std::vector<float> phiL, phiM, phiT;
   std::vector<float> massL, massM, massT;
   std::vector<float> energyL, energyM, energyT;
   std::vector<float> chargeL, chargeM, chargeT;
   std::vector<float> dzL, dzM, dzT;
   std::vector<float> dxyL, dxyM, dxyT;
   std::vector<float> mvaL, mvaM, mvaT;
   std::vector<float> relIsoL, relIsoM, relIsoT;
   
   void clearTreeVectors() {
      nL.     clear();
      ptL.    clear();
      etaL.   clear();
      phiL.   clear();
      massL.  clear();
      energyL.clear();
      dzL.    clear();
      dxyL.   clear();
      mvaL.   clear();
      relIsoL.clear();

      nM.     clear();
      ptM.    clear();
      etaM.   clear();
      phiM.   clear();
      massM.  clear();
      energyM.clear();
      dzM.    clear();
      dxyM.   clear();
      mvaM.   clear();
      relIsoM.clear();

      nT.     clear();
      ptT.    clear();
      etaT.   clear();
      phiT.   clear();
      massT.  clear();
      energyT.clear();
      dzT.    clear();
      dxyT.   clear();
      mvaT.   clear();
      relIsoT.clear();
   }

   void RegisterTree(TTree* tree, std::string name="Electons") {
      tree->Branch((name+"_nL")     .c_str(), &nL);
      tree->Branch((name+"_ptL")    .c_str(), &ptL);
      tree->Branch((name+"_etaL")   .c_str(), &etaL); 
      tree->Branch((name+"_phiL")   .c_str(), &phiL);
      tree->Branch((name+"_massL")  .c_str(), &massL);
      tree->Branch((name+"_energyL").c_str(), &energyL); 
      tree->Branch((name+"_chargeL").c_str(), &chargeL);
      tree->Branch((name+"_dzL")    .c_str(), &dzL);
      tree->Branch((name+"_dxyL")   .c_str(), &dxyL); 
      tree->Branch((name+"_mvaL")   .c_str(), &mvaL);
      tree->Branch((name+"_relIsoL").c_str(), &relIsoL);

      tree->Branch((name+"_nM")     .c_str(), &nM);
      tree->Branch((name+"_ptM")    .c_str(), &ptM);
      tree->Branch((name+"_etaM")   .c_str(), &etaM); 
      tree->Branch((name+"_phiM")   .c_str(), &phiM);
      tree->Branch((name+"_massM")  .c_str(), &massM);
      tree->Branch((name+"_energyM").c_str(), &energyM); 
      tree->Branch((name+"_chargeM").c_str(), &chargeM);
      tree->Branch((name+"_dzM")    .c_str(), &dzM);
      tree->Branch((name+"_dxyM")   .c_str(), &dxyM); 
      tree->Branch((name+"_mvaM")   .c_str(), &mvaM);
      tree->Branch((name+"_relIsoM").c_str(), &relIsoM);

      tree->Branch((name+"_nT")     .c_str(), &nT);
      tree->Branch((name+"_ptT")    .c_str(), &ptT);
      tree->Branch((name+"_etaT")   .c_str(), &etaT); 
      tree->Branch((name+"_phiT")   .c_str(), &phiT);
      tree->Branch((name+"_massT")  .c_str(), &massT);
      tree->Branch((name+"_energyT").c_str(), &energyT); 
      tree->Branch((name+"_chargeT").c_str(), &chargeT);
      tree->Branch((name+"_dzT")    .c_str(), &dzT);
      tree->Branch((name+"_dxyT")   .c_str(), &dxyT); 
      tree->Branch((name+"_mvaT")   .c_str(), &mvaT);
      tree->Branch((name+"_relIsoT").c_str(), &relIsoT);  
   }
};
