#!/usr/bin/env python
import os, sys
import numpy as np
from ROOT import gROOT,std,ROOT,TFile,TTree,TH1D,TH2D,TStopwatch,TMatrix,TLorentzVector,TMath,TVector
gROOT.Macro("~/rootlogon.C")

def Overlaps2D(jetsP4, lepP4, drMax, ptrelMin):
   overlaps = False
   for j in range(0, len(jetsP4)): 
      #print j    
      dr = jetsP4[j].DeltaR(lepP4)
      #print 'loop over P4 jet and dr = ', dr
      jp3 = jetsP4[j].Vect()
      lp3 = lepP4.Vect()
      ptRel = jetsP4[j].Perp(lp3)
      dPtRel = (jp3.Cross( lp3 )).Mag()/ jp3.Mag()
      #print 'in function: dr = ', dr, 'dPtRel = ', dPtRel
      if(dr < drMax and ptRel < ptrelMin) : 
      #if not (dr > drMax or dPtRel > ptrelMin) : 
         overlaps = True
         #print 'in function: dr = ', dr, 'dPtRel = ', dPtRel
         return overlaps
   print 'overlap ?? ', overlaps
   return overlaps 
      
from optparse import OptionParser

# Create a command line option parser
options  = OptionParser()

options.add_option('--inDir', metavar='T', type='string', action='store',
                  default='inDir', 
                  dest='input directory', 
                  help='input data directory name')
options.add_option('-f', '--files',  
                   dest="files", 
                   default="T_M1000_W10.txt",
                   type="string")
options.add_option('-n', '--maxEvts',  
                   dest="maxEvts", 
                   default=-1,
                   type="int")

(options,args) = options.parse_args()
# ==========end: options =============

##options.add_option_group(evtsel)

#opt, remainder = options.parse_args()

print options
maxEvts = options.maxEvts

# Define the output histograms
fout = TFile(options.files.rstrip().replace('.txt', '_out.root'), 'RECREATE')
fout.cd()

hCutflow = TH1D("hCutflow" ,";;Events;" ,10, 0.5, 10.5)
cutsName = ['Total', '== 1 lep', '2D lep Iso', 'N(fjet) #geq 1', 'N(jet) #geq 2', 'leading jet pt > 200', '2nd jet pt > 80', 'N(b jet) #geq 1', 'MET #geq 20', 'N(Higgs) #geq 1']
ibin = 0
for n in cutsName:
   ibin = ibin+1
   hCutflow.GetXaxis().SetBinLabel(ibin, n)

 
hlepPt     = TH1D("hlepPt",  "Lepton p_{T}; p_{T}(GeV); Events/60 GeV;", 50, 0, 300)
hlepEta    = TH1D("hlepEta", "Lepton #eta; #eta; Events/10 bins;", 80, -4.0, 4.0)
hlepIso    = TH1D("hlepIso", "Lepton isolation; Isolation; Events;", 210, -100.0, 5.0)
hlepIso_sig    = TH1D("hlepIso_sig", "Lepton isolation; Isolation; Events;", 210, -100.0, 5.0)
hjetsPt    = TH1D("hJetsPt", "Jets Pt; Jets Pt (GeV); Events/15 GeV", 80, 0, 2000)
hjetsEta   = TH1D("hJetsEta", "Jets #eta; Jets #eta; Events", 50, -5.0, 5.0)
hDRMin     = TH1D("hDRMin", "#Delta R_{MIN}(l, jet); #Delta R_{MIN}(l,jet); Events", 30, 0.0,3.0)
hDR        = TH1D("hDR", "#Delta R(l, jet); #Delta R (l,jet); Events", 30, 0.0,3.0) 
hPtRel     = TH1D("hPtRel", "p_{T,rel}; p_{T,rel} (GeV); Events/20 GeV", 50, 0, 100); 
hDPtRel    = TH1D("hDPtRel", "#Delta p_{T}^{REL}; #Delta p_{T}^{REL} (GeV); Events/20 GeV",50, 0, 100)   
hDelPtRel  = TH1D("hDelPtRel", "#Delta p_{T}^{REL}; #Delta p_{T}^{REL} [GeV]; Events/20 GeV;", 50, 0, 100) 
h2DdPtReldR = TH2D("h2DdPtReldR", ";#Delta R(l,j); #Delta p_{T}^{REL} [GeV]", 50, 0.0, 1.0, 20, 0., 200.)
h2DdPtRelDRMin = TH2D("h2DdPtRelDRMin", ";#Delta R_{MIN}(l,j); min #Delta p_{T}^{REL} (GeV)", 50, 0.0, 1.0, 20, 0., 200.)
h2DPtRelDRMin = TH2D("h2DPtRelDRMin",";#Delta R_{MIN}(l,j); p_{T,rel} (GeV)", 50, 0.0, 1.0, 20, 0., 200.) 

lepP4        = TLorentzVector(0.0, 0.0, 0.0, 0.0)
#jetP4        = TLorentzVector(0.0, 0.0, 0.0, 0.0)
nearestJetP4 = TLorentzVector(0.0, 0.0, 0.0, 0.0)

# Open the input ntuples
fnames = [line.strip() for line in open(options.files, 'r')]

# Begin running over all trees
ievt = 0
for fname in fnames:
  if maxEvts > 0 and ievt > maxEvts: break
  #if ievt%100 == 0: print " Processing evt %i" % ievt

  print 'Opening file %s' % fname
  f = TFile.Open(fname)
  print f.ls()

  tree = f.Get("ana/anatree")
  entries = tree.GetEntriesFast()
  
  for t in tree:
    ncut = 0
    if maxEvts > 0 and ievt > maxEvts: break
    if ievt%100 == 0: print " Processing evt %i" % ievt
    
    ievt += 1
    #ncut += 1
    # call the Gen weights
    evtwt = t.GenEvt_genWt
    #hCutflow.Fill(ncut, evtwt)
    
    # require at least one good primary vertix
    vertices = t.SelectedEvt_nGoodVtx
    if vertices <= 0: continue;
    ncut += 1
    hCutflow.Fill(ncut, evtwt)
    
  # require exactly one medium WP electron or muon in event


  # loop on electrons, selecting medium WP electrons
    nMedium=0
    ele_eta = t.Electrons_eta
    ele_phi = t.Electrons_phi
    ele_pt  = t.Electrons_pt
    ele_m   = t.Electrons_mass
    ele_mva = t.Electrons_mva
    ele_iso = t.Electrons_relIso
    #print 'size of electrons :' , len(t.Electrons_pt)
    for i in range(0, len(t.Electrons_pt)):
      # choose the WP: 1 = loose, 2 = medium, 4 = tight
      eWP = t.Electrons_eleWP[i]    
      if (eWP & 1) != 1: continue
      if ele_pt[i] < 40.: continue 
      if abs(ele_eta[i]) < 2.8: continue
      nMedium += 1
      #print 'mva = ',  ele_mva[i], 'eta = ', abs(ele_eta[i]), 'pt = ', ele_pt[i], 'eWP = ', eWP 

    #require exactly one lepton   
    if nMedium != 1: continue 
    ncut += 1
    hCutflow.Fill(ncut, evtwt)
    
    # store lep variables:
    hlepPt.Fill(ele_pt[0], evtwt)
    hlepEta.Fill(ele_eta[0], evtwt)
    hlepIso.Fill(ele_iso[0], evtwt) 
    lepP4.SetPtEtaPhiM(ele_pt[0], ele_eta[0], ele_phi[0], ele_m[0])
    lep_p3 = lepP4.Vect();

    # loop over jets 
    ak4jet_pt  = t.AK4JetsCHS_pt 
    ak4jet_eta = t.AK4JetsCHS_eta #barrel: eta< 1.479
    ak4jet_phi = t.AK4JetsCHS_phi
    ak4jet_m   = t.AK4JetsCHS_mass 
    nak4jet = len(ak4jet_pt)

    jetsP4 = [] #define a list to store P4 of all good jets
    dR = 900.0; dRMin = 999.0; delPtRel = 999.0
    
    for j in range(0,nak4jet):
      #funny enough: if I define the jetP4 in this line, then the address of object is never changed, hence its content
      hjetsPt.Fill(ak4jet_pt[j], evtwt)
      hjetsEta.Fill(ak4jet_eta[j], evtwt)
      if ak4jet_pt[j] < 30. : continue
      if ak4jet_eta[j] > 5.0 : continue
      jetP4 = TLorentzVector(0.0, 0.0, 0.0, 0.0)
      jetP4.SetPtEtaPhiM(ak4jet_pt[j], ak4jet_eta[j], ak4jet_phi[j], ak4jet_m[j])     
      #print 'what is going inside: ', jetP4.Pt()
      jetsP4.append(jetP4)
      dR = jetP4.DeltaR(lepP4)      
      jet_p3 = jetP4.Vect()       
      delPtRel = (jet_p3.Cross( lep_p3 )).Mag()/ jet_p3.Mag()
      hDR.Fill(dR, evtwt)
      hDelPtRel.Fill(delPtRel, evtwt)
      h2DdPtReldR.Fill(dR, delPtRel, evtwt)
      #print 'in main jet loop: dr = ', dR, 'dPtRel = ', delPtRel  
      if dR < dRMin:
        nearestJetP4 = jetP4
        dRMin = dR
      #print 'what is stored inside: ', jetsP4, 'with j :', j

    # Store extra variables 
    ptRel = nearestJetP4.Perp( lepP4.Vect() )
    hPtRel.Fill(ptRel, evtwt)
    
    hDRMin.Fill(dRMin, evtwt)

    nearestJet_p3 = nearestJetP4.Vect()
    dPtRel = (nearestJet_p3.Cross( lep_p3 )).Mag()/ nearestJet_p3.Mag()   
    hDPtRel.Fill(dPtRel, evtwt)

    h2DPtRelDRMin.Fill(dRMin, ptRel, evtwt)

    h2DdPtRelDRMin.Fill(dRMin, dPtRel, evtwt)

    pass2D = ptRel > 20. or dRMin > 0.4
    if not pass2D: continue
    #if Overlaps2D(jetsP4, lepP4, 0.4, 20): continue
    
    ncut += 1
    hCutflow.Fill(ncut, evtwt)
    # separate the jets into central and forward jet collections  
    
    hlepIso_sig.Fill(ele_iso[0], evtwt) 
      
    del jetsP4[:] 
    #print ievt
fout.Write()
fout.Close()
