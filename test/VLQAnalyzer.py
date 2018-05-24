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

def SolveNuPz(vlep, vnu, wmass, nuz1, nuz2):
    solutionBool = True
    x = vlep.X()*vnu.X() + vlep.Y()*vnu.Y() + wmass*wmass/2
    a = vlep.Z()*vlep.Z() - vlep.E()*vlep.E()
    b = 2*x*vlep.Z()
    c = x*x - vnu.Perp2() * vlep.E()*vlep.E()
    d = b*b - 4*a*c
    #print 'x: ', x
    #print 'a: ', a
    #print 'b: ', b
    #print 'c: ', c
    #print 'd: ', d
    if d < 0:
        d = 0
        solutionBool = False
    nuz1 = (-b + np.sqrt(d))/2/a
    nuz2 = (-b -np.sqrt(d))/2/a
    #print 'nuz1: ', nuz1
    #print 'nuz2: ', nuz2
    if abs(nuz1) > abs(nuz2):
        lowsol = nuz2
        highsol = nuz1
    else:
        lowsol = nuz1
	highsol = nuz2
    return solutionBool, lowsol, highsol;

def AdjustEnergyForMass(v, mass):
    v.SetE(np.sqrt(v.Vect().Mag2() + mass*mass))
    return;

def nextPerm(lis):
    n = len(lis)

    i = n-2
    while i >= 0 and lis[i] >= lis[i+1]:
	i -= 1

    if i == -1:
	return False;

    j = i + 1
    while j < n and lis[j] > lis[i]:
	j += 1
    j -= 1

    lis[i], lis[j] = lis[j], lis[i]

    left = i + 1
    right = n -1

    while left < right:
	lis[left], lis[right] = lis[right], lis[left]
	left += 1
	right -= 1

    return True;

def GetChi2(JetsP4, LeptonP4, NuP4, topMass, higgsMass, topP4, higgsP4, dR):
    top = 0.0; higgs = 0.0; top_chi2 = 0.0; higgs_chi2 = 0.0; dR_topH = 0.0; dR_topH_chi2 = 0.0
    top = abs( (JetsP4[0] + LeptonP4 + NuP4).M() - topMass)
    top_chi2 = top*top / (14.5*14.5) #where did 14.5 come from?
    topP4 = JetsP4[0] + LeptonP4 + NuP4

    higgs = abs((JetsP4[1] + JetsP4[2]).M() - higgsMass)
    higgs_chi2 = higgs*higgs / (14.5*14.5) #again, where did 14.5 come from?
    higgsP4 = JetsP4[1] + JetsP4[2]

    dR_topH = abs( topP4.DeltaR(higgsP4) - 3.15) #3.15
    dR_topH_chi2 = dR_topH*dR_topH / (0.196*0.196) #0.196
    dR = topP4.DeltaR(higgsP4)

    if JetsP4[0].DeltaR(higgsP4) < 1.0:
	return 100000.0, 100000.0;
    else:
	return top_chi2 + higgs_chi2 + dR_topH_chi2, dR
    
def GetChi2Boosted(ak4jetsP4, higgsJet, higgsSoftDropMass, leptonP4, nuP4, topMass, higgsMass, topP4, higgsP4, dR):
    top = 0.0; higgs = 0.0; top_chi2 = 0.0; higgs_chi2 = 0.0; dR_topH = 0.0; dR_topH_chi2 = 0.0
    
    top = abs( (ak4jetsP4[0] + leptonP4 + nuP4).M() - topMass)
    top_chi2 = top*top / (14.5*14.5)
    topP4 = ak4jetsP4[0] + leptonP4 + nuP4

    higgs = abs(higgsSoftDropMass[0] - higgsMass)
    higgs_chi2 = higgs*higgs / (14.5*14.5)
    higgsP4 = higgsJet[0]

    dR_topH = abs(topP4.DeltaR(higgsP4) - 3.15)
    dR_topH_chi2 = dR_topH*dR_topH / (0.196*0.196)
    #print 'infuncdRBefore: ', dR
    dR = topP4.DeltaR(higgsP4)
    #print 'infuncdRAfter: ', dR

    if ak4jetsP4[0].DeltaR(higgsP4) < 1.0:
	return 100000.0, 100000.0, TLorentzVector(100000.0, 100000.0, 100000.0, 100000.0), TLorentzVector(100000.0, 100000.0, 100000.0, 100000.0)
    else:
	return top_chi2 + higgs_chi2 + dR_topH_chi2, dR, topP4, higgsP4
 
def DoMassReco(jetColl, LeptonP4, NuP4, higgsMass, topMass, chi2_dR1, chi2_dR2, chi2_higgs1, chi2_higgs2, chi2_top1, chi2_top2):
    chi2_higgs1 = 100000.0; chi2_top1 = 100000.0
    JetsP4 = []
    index_list = [0, 1, 2, 3]
    index_list1 = [0, 1, 2]

    chi2 = 100000.0; dR = 100000.0; minChi2 = 100000.0

    topP4 = TLorentzVector(0.0, 0.0, 0.0, 0.0)
    higgsP4 = TLorentzVector(0.0, 0.0, 0.0, 0.0)

    if len(jetColl) > 3:
	while True:
	    i0 = index_list[0]
	    i1 = index_list[1]
	    i2 = index_list[2]
    	    i3 = index_list[3]

	    jetsP4PassToChi2 = []

	    if len(jetsP4PassToChi2) != 0:
	        jetsP4PassToChi2.clear()
	    topP4.Clear()
	    higgsP4.Clear()
	    JetsP4 = []

	    JetsP4.append(jetColl[i0])
	    JetsP4.append(jetColl[i1])
	    JetsP4.append(jetColl[i2])
	    JetsP4.append(jetColl[i3])
	    jetsP4PassToChi2.append(JetsP4[0])
	    jetsP4PassToChi2.append(JetsP4[1])
	    jetsP4PassToChi2.append(JetsP4[2])
	    jetsP4PassToChi2.append(JetsP4[3])

	    chi2, dR = GetChi2(jetsP4PassToChi2, LeptonP4, NuP4, topMass, higgsMass, topP4, higgsP4, dR)

	    if (chi2 < minChi2):
	        minChi2 = chi2
	    
	        chi2_higgs1 = minChi2
	        chi2_higgs2 = higgsP4

	        chi2_top1 = minChi2
	        chi2_top2 = topP4

	        chi2_dR1 = minChi2
	        chi2_dR2 = dR
	    if not nextPerm(index_list):
	        break
    elif len(jetColl) == 3:
	while True:
	    i0 = index_list[0]
	    i1 = index_list[1]
	    i2 = index_list[2]

	    jetsP4PassToChi2 = []

	    if len(jetsP4PassToChi2) != 0:
		jetsP4PassToChi2.clear()
	    topP4.Clear()
	    higgsP4.Clear()
	    JetsP4 = []

	    JetsP4.append(jetColl[i0])
            JetsP4.append(jetColl[i1])
            JetsP4.append(jetColl[i2])

	    jetsP4PassToChi2.append(JetsP4[0])
	    jetsP4PassToChi2.append(JetsP4[1])
	    jetsP4PassToChi2.append(JetsP4[2])

	    chi2, dR = GetChi2(jetsP4PassToChi2, LeptonP4, NuP4, topMass, higgsMass, topP4, higgsP4, dR)

	    if chi2 < minChi2:
		minChi2 = chi2
		
		chi2_higgs1 = minChi2
		chi2_higgs2 = higgsP4

		chi2_top1 = minChi2
		chi2_top2 = topP4

		chi2_dR1 = minChi2
		chi2_dR2 = dR
	    if not nextPerm(index_list1):
                break

    return chi2_higgs1, chi2_higgs2, chi2_top2, chi2_dR2;

def DoMassRecoBoost(ak4Jets, higgsJet, higgsSoftDropMass, leptonP4, nuP4, higgsMass, topMass, chi2_dR1, chi2_dR2, chi2_higgs1, chi2_higgs2, chi2_top1, chi2_top2):
    ak4JetsP4 = [TLorentzVector(0.0, 0.0, 0.0, 0.0), TLorentzVector(0.0, 0.0, 0.0, 0.0)]
    index_list = [0, 1]
    chi2 = 100000.0; dR = 100000.0; minChi2 = 100000.0
    topP4 = TLorentzVector(0.0, 0.0, 0.0, 0.0)
    higgsP4 = TLorentzVector(0.0, 0.0, 0.0, 0.0)

    if len(higgsJet) > 0:
	while True:
	    i0 = index_list[0]
	    i1 = index_list[1]
	    jetsP4PassToChi2 = []

	    if len(jetsP4PassToChi2) != 0:
                jetsP4PassToChi2.clear()
	    topP4.Clear()
	    higgsP4.Clear()

	    ak4JetsP4[0] = ak4Jets[i0]
	    if len(ak4Jets) > 1:
		ak4JetsP4[1] = ak4Jets[i1]
	    else:
		ak4JetsP4[1].SetPtEtaPhiM(0.0, 0.0, 0.0, 0.0)

	    jetsP4PassToChi2.append(ak4Jets[i0])
	    jetsP4PassToChi2.append(ak4Jets[i1])

            chi2, dR, topP4, higgsP4 = GetChi2Boosted(jetsP4PassToChi2, higgsJet, higgsSoftDropMass, leptonP4, nuP4, topMass, higgsMass, topP4, higgsP4, dR)

	    if chi2 < minChi2:
		minChi2 = chi2

		chi2_higgs1 = minChi2
		chi2_higgs2 = higgsP4

		chi2_top1 = minChi2
		chi2_top2 = topP4

		chi2_dR1 = minChi2
		chi2_dR2 = dR		
	    if not nextPerm(index_list):
                break

    return chi2_higgs1, chi2_higgs2, chi2_top2, chi2_dR2

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
h2DPtRelDRMin = TH2D("h2DPtRelDRMin", ";#Delta R_{MIN}(l,j); p_{T,rel} (GeV)", 50, 0.0, 1.0, 20, 0., 200.) 
hNForwardJets = TH1D("hNForwardJets", "Number of Forward Jets; Number of Forward Jets; Events;", 40, 0, 40)
hLeadingJetPt = TH1D("hLeadingJetPt", "Leading Jet p_{T}; p_{T} {GeV}; Events;", 50, 0, 1000)
hLeadingJetEta = TH1D("hLeadingJetEta", "Leading Jet #eta; #eta; Events;", 80, -4.0, 4.0)
hNCentJets = TH1D("hNCentJets", "Number of Central Jets; Number of Central Jets; Events;", 40, 0, 40)
hSecJetPt = TH1D("hSecJetPt", "Second Leading Jet p_{T}; p_{T} {GeV}; Events;", 50, 0, 1000)
hSecJetEta = TH1D("hSecJetEta", "Second Leading Jet #eta; #eta; Events;", 80, -4.0, 4.0)
hLargestDeepCSVB = TH1D("hLargestDeepCSVB", "Largest B-Tag Value from DeepCSV; B-Tag Value; Events;", 500, 0, 1.5)
hMET = TH1D("hMET", "Missing E_{T}; MET {GeV}; Events;", 50, 0, 1000)
hak8JetPt = TH1D("hak8JetPt", "ak8Jet p_{T}; p_{T} {GeV}; Events;", 50, 0, 1000)
hak8JetEta = TH1D("hak8JetEta", "ak8Jet #eta; #eta; Events;", 80, -4.0, 4.0)
hhiggsJetPt = TH1D("hhiggsJetPt", "higgsJet p_{T}; p_{T} {GeV}; Events;", 50, 0, 1000)
hhiggsJetEta = TH1D("hhiggsJetEta", "higgsJet #eta; #eta; Events;", 80, -4.0, 4.0)

lepP4        = TLorentzVector(0.0, 0.0, 0.0, 0.0)
jetP4        = TLorentzVector(0.0, 0.0, 0.0, 0.0)
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
      if abs(ele_eta[i]) > 2.8: continue
      nMedium += 1
      #print 'mva = ',  ele_mva[i], 'eta = ', abs(ele_eta[i]), 'pt = ', ele_pt[i], 'eWP = ', eWP 

  # loop on muons, selecting medium WP muons
    nmuMedium=0
    mu_eta = t.Muons_eta
    mu_phi = t.Muons_phi
    mu_pt  = t.Muons_pt
    mu_m   = t.Muons_mass
    mu_iso = t.Muons_relIso
    #print 'size of muons :' , len(t.Muons_pt)
    for i in range(0, len(t.Muons_pt)):
      # choose the WP: 1 = loose, 2 = medium, 4 = tight
      muWP = t.Muons_muWP[i]
      if (muWP & 1) != 1: continue
      if mu_pt[i] < 40.: continue
      if abs(mu_eta[i]) > 4.0: continue
      nmuMedium += 1
      #print 'mva = ',  mu_mva[i], 'eta = ', abs(mu_eta[i]), 'pt = ', mu_pt[i], 'muWP = ', muWP 

    #require exactly one lepton   
    if nMedium + nmuMedium != 1: continue 
    ncut += 1
    hCutflow.Fill(ncut, evtwt)
    
    # store lep variables:
    if len(ele_pt) > 0:
    	hlepPt.Fill(ele_pt[0], evtwt)
    	hlepEta.Fill(ele_eta[0], evtwt)
    	hlepIso.Fill(ele_iso[0], evtwt) 
    	lepP4.SetPtEtaPhiM(ele_pt[0], ele_eta[0], ele_phi[0], ele_m[0])
    else:
	hlepPt.Fill(mu_pt[0], evtwt)
        hlepEta.Fill(mu_eta[0], evtwt)
        hlepIso.Fill(mu_iso[0], evtwt)
        lepP4.SetPtEtaPhiM(mu_pt[0], mu_eta[0], mu_phi[0], mu_m[0])
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
      if abs(ak4jet_eta[j]) > 5.0 : continue
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


    if len(ele_iso) > 0:
    	hlepIso_sig.Fill(ele_iso[0], evtwt)
    else:
	hlepIso_sig.Fill(mu_iso[0], evtwt) 

    # separate the jets into central and forward jet collections    
    centjets = []
    fjets = []
    for j in jetsP4:
	jet1 = j
	eta1 = j.Eta()
	if abs(eta1) < 2.4:
	    centjets.append(jet1)
	else:
	    fjets.append(jet1)
 
    #print 'Length of good jets: ', len(jetsP4)
    #print 'Length of Btag list: ', len(jetsBtag)

    # histograms showing number of forward and central jets before jet cuts
    hNForwardJets.Fill(len(fjets), evtwt)
    hNCentJets.Fill(len(centjets), evtwt)

    # 1 or more forward jets
    if len(fjets) < 1: continue
    ncut += 1
    hCutflow.Fill(ncut, evtwt)
    # 3 or more central jets
    if len(centjets) < 3: continue
    ncut += 1
    hCutflow.Fill(ncut, evtwt)

    if centjets[0].Pt() < centjets[1].Pt(): print 'not ordered properly'

    # leading jet pt > 200
    leadjet = centjets[0]
    leadjetpt = leadjet.Pt()
    leadjeteta = leadjet.Eta()
    if leadjetpt <= 200: continue
    ncut += 1
    hCutflow.Fill(ncut, evtwt)

    # second leading jet pt > 80
    secondjet = centjets[1]
    secondjetpt = secondjet.Pt()
    secondjeteta = secondjet.Eta()
    if secondjetpt <= 80: continue
    ncut += 1
    hCutflow.Fill(ncut, evtwt)

    # other variables and histograms with jets
    hLeadingJetPt.Fill(leadjetpt, evtwt)
    hLeadingJetEta.Fill(leadjeteta, evtwt)
    hSecJetPt.Fill(secondjetpt, evtwt)
    hSecJetEta.Fill(secondjeteta, evtwt)
   
    # cut out events without 1 B jet (using deepcsvm value) 
    ak4jet_deepcsv = t.AK4JetsCHS_deepcsv
    ngoodjets = len(jetsP4)
    #print 'ngoodjets: ', ngoodjets
    jetsBtag = []
    deepcsvm = 0.4941
    nMedBjets = 0
    largestB = 0
    for j in range(0, ngoodjets):
	jetsBtag.append(ak4jet_deepcsv[j])
    nbseljets = len(jetsBtag)
    #print 'nbseljets: ', nbseljets
    for j in jetsBtag:
        #print 'Btag values before cut: ', j
        if j < deepcsvm: continue
        nMedBjets += 1
	#print 'Btag values after cut: ', j
        if j > largestB:
	    largestB = j
    if nMedBjets == 0: continue
    ncut += 1
    hCutflow.Fill(ncut, evtwt)
    #print 'largestB: ', largestB
    hLargestDeepCSVB.Fill(largestB, evtwt)

    # Met > 20 GeV
    met_pt = t.MET_pt
    met_px = t.MET_px
    met_py = t.MET_py
    met_pz = t.MET_pz
    if met_pt <= 20: continue
    ncut += 1
    hCutflow.Fill(ncut, evtwt)
    hMET.Fill(met_pt, evtwt)

    #print 'met_px: ', met_px
    #print 'met_py: ', met_py
    #print 'met_pz: ', met_pz
 
    # Higgs Tagging
    ak8jet_pt = t.AK8Jets_ptCHS
    ak8jet_eta = t.AK8Jets_etaCHS
    ak8jet_phi = t.AK8Jets_phiCHS
    ak8jet_m = t.AK8Jets_massCHS
    ak8jet_tau1 = t.AK8Jets_tau1CHS
    ak8jet_tau2 = t.AK8Jets_tau2CHS
    ak8jet_sj1pt = t.AK8Jets_sj0pt
    ak8jet_sj2pt = t.AK8Jets_sj1pt
    ak8jet_sdmass = t.AK8Jets_softDropMassCHS
    ak8jet_sj1deepcsv = t.AK8Jets_sj0deepcsv
    ak8jet_sj2deepcsv = t.AK8Jets_sj1deepcsv
    nak8jet = len(ak8jet_pt)
    nHiggs = 0

    ak8jetsP4 = []
    higgsjets = []
    higgsSoftDropM = []

    for j in range(0, nak8jet):
	# fill histos with ak8pt and ak8eta
	hak8JetPt.Fill(ak8jet_pt[j], evtwt)
	hak8JetEta.Fill(ak8jet_eta[j], evtwt)
	#print 'ak8 tau2: ', ak8jet_tau2[j]
	#print 'ak8 tau1: ', ak8jet_tau1[j]
	#print 'subjettiness: ', ak8jet_tau2[j]/ak8jet_tau1[j]
	#print 'length of subjet1: ', len(ak8jet_sj1pt)
	#print 'length of subjet2: ', len(ak8jet_sj2pt)
	if ak8jet_pt[j] < 300.0: continue
	if abs(ak8jet_eta[j]) > 2.4: continue
	if (ak8jet_tau2[j]/ak8jet_tau1[j]) > 0.6: continue
	if len(ak8jet_sj1pt) == 0 and len(ak8jet_sj2pt) == 0: continue
	ak8jetP4 = TLorentzVector(0.0, 0.0, 0.0, 0.0)
	ak8jetP4.SetPtEtaPhiM(ak8jet_pt[j], ak8jet_eta[j], ak8jet_phi[j], ak8jet_m[j])
	#print 'still pt here: ', ak8jetP4.Pt()
	#print 'deltaR: ', ak8jetP4.DeltaR(lepP4)
	#print 'softdropmass: ', ak8jet_sdmass[j]
	if ak8jetP4.DeltaR(lepP4) <= 1.0: continue
	if ak8jet_sdmass[j] > 160 or ak8jet_sdmass[j] < 90: continue
       
	# b tagging of subjets
	if ak8jet_sj1deepcsv[j] < deepcsvm: continue
	if ak8jet_sj2deepcsv[j] < deepcsvm: continue
	#print 'subjet1 deepcsv: ', ak8jet_sj1deepcsv[j]
	#print 'subjet2 deepcsv: ', ak8jet_sj2deepcsv[j]
 
	nHiggs += 1

	higgsjets.append(ak8jetP4)
	higgsSoftDropM.append(ak8jet_sdmass[j])
	hhiggsJetPt.Fill(ak8jetP4.Pt(), evtwt)
	hhiggsJetEta.Fill(ak8jetP4.Eta(), evtwt)

    if nHiggs == 0: continue
    ncut += 1
    hCutflow.Fill(ncut, evtwt)

    # Top mass reconstruction
    nuP4 = TLorentzVector(met_px, met_py, met_pz, 0.0)
    lowersol = 0.0; highersol = 0.0
    isNuPz, lowersol, highersol = SolveNuPz(lepP4, nuP4, 80.4, lowersol, highersol)

    #print 'sollow: ', lowersol
    #print 'solhigh: ', highersol
    #print 'isNuPz: ', isNuPz

    nuP4.SetPz(lowersol)
    #print 'nuP4before: ', nuP4.Px(), nuP4.Py(), nuP4.Pz(), nuP4.E()
    AdjustEnergyForMass(nuP4, 0.0)
    #print 'nuP4after: ', nuP4.Px(), nuP4.Py(), nuP4.Pz(), nuP4.E()

    topMass = 174.0; higgsMass = 125.0; chi2_dR1 = 0.0; chi2_dR2 = 100000.0; chi2_higgs1 = 0.0; chi2_top1 = 0.0; chi2_dR_boost1 = 0.0; chi2_dR_boost2 = 0.0; chi2_higgs_boost1 = 0.0; chi2_top_boost1 = 0.0
    chi2_higgs2 = TLorentzVector(0.0, 0.0, 0.0, 0.0)
    chi2_top2 = TLorentzVector(0.0, 0.0, 0.0, 0.0)
    chi2_higgs_boost2 = TLorentzVector(0.0, 0.0, 0.0, 0.0)
    chi2_top_boost2 = TLorentzVector(0.0, 0.0, 0.0, 0.0)
    higgsP4 = TLorentzVector(0.0, 0.0, 0.0, 0.0)
    topP4 = TLorentzVector(0.0, 0.0, 0.0, 0.0)
    chi2 = 100000.0

    print 'beforechi2: ', chi2
    print 'beforedR: ', dR
    print 'beforetopP4: ', topP4.Px(), topP4.Py(), topP4.Pz(), topP4.E()
    print 'beforehiggsP4: ', higgsP4.Px(), higgsP4.Py(), higgsP4.Pz(), higgsP4.E()
 
    chi2, higgsP4, topP4, dR = DoMassRecoBoost(centjets, higgsjets, higgsSoftDropM, lepP4, nuP4, higgsMass, topMass, chi2_dR_boost1, chi2_dR_boost2, chi2_higgs_boost1, chi2_higgs_boost2, chi2_top_boost1, chi2_top_boost2)

    print 'afterchi2: ', chi2
    print 'afterdR: ', dR
    print 'aftertopP4: ', topP4.Px(), topP4.Py(), topP4.Pz(), topP4.E()
    print 'afterhiggsP4: ', higgsP4.Px(), higgsP4.Py(), higgsP4.Pz(), higgsP4.E()

    if dR < 2.0: continue # we want to get rid of these results, because it is difficult to discern between the top and higgs in these cases

    # Now have the topP4 and higgsP4. So we can reconstruct the Tprime quark and calculate a mass

    

    del jetsP4[:]
    del fjets[:]
    del centjets[:]
    #print ievt

fout.Write()
fout.Close()

# Create Mass Reconstruction Function
# Create Boosted Mass Reconstruction Function
# Create Chi2 Function
# Create Boosted Chi2 Function
