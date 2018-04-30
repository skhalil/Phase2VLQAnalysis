#!/usr/bin/env python
import os, sys
import numpy as np
from ROOT import gROOT,std,ROOT,TFile,TTree,TH1D,TStopwatch,TMatrix,TLorentzVector,TMath,TVector
gROOT.Macro("~/rootlogon.C")

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
                   type="string",
                   )
options.add_option('-n', '--maxEvts',  
                   dest="maxEvts", 
                   default=-1,
                   type="int",
                   )

(options,args) = options.parse_args()
# ==========end: options =============

##options.add_option_group(evtsel)

#opt, remainder = options.parse_args()

print options
maxEvts = options.maxEvts

# Define the output histograms
# This assumes that the input and output files are in 
# different directories, otherwise the input will be overwritten 
# by the output directory name
fout = TFile(options.files.rstrip().replace('.txt', '_out.root'), 'RECREATE')
fout.cd()

hCutflow = TH1D("hCutflow" ,";;Events;" ,10 ,0.5 ,10.5 )
hCutflow.GetXaxis().SetBinLabel(1 ,"All evts") ;
hCutflow.GetXaxis().SetBinLabel(2 ,"good PV") ;
hCutflow.GetXaxis().SetBinLabel(3 ,"== 1 lepton") ;

 
hlepPt  = TH1D("hlepPt" ,";;Lepton p_{T}(GeV);" ,50 ,0 ,300  )
hlepEta = TH1D("hlepEta" ,";;Lepton #eta;" ,80 ,-4 ,4  )

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
    
    if maxEvts > 0 and ievt > maxEvts: break
    if ievt%100 == 0: print " Processing evt %i" % ievt
    
    ievt += 1

    # call the Gen weights
    evtwt = t.GenEvt_genWt
    hCutflow.Fill(1, evtwt)
    
    # require at least one good primary vertix
    vertices = t.SelectedEvt_nGoodVtx
    if vertices > 0: hCutflow.Fill(2, evtwt)
    else: continue
    
  # require exactly one medium WP electron or muon in event


  # loop on electrons, selecting medium WP electrons
    nMedium=0
    ele_eta = t.Electons_eta
    ele_pt = t.Electons_pt

    for i in range(0, len(t.Electons_pt)):
      # choose the WP: 1 = loose, 2 = medium, 4 = tight
      eWP = t.Electons_eleWP[i]    
      if (eWP & 2) != 2: continue
      # basic kinematics
      if ele_pt[i] < 20: continue 
      if abs(ele_eta[i]) < 2.8: continue
      hlepPt.Fill(ele_pt[i])
      hlepEta.Fill(ele_eta[i])
      nMedium += 1
       
    # require exactly one lepton
    
    if nMedium == 1: 
      hCutflow.Fill(3, evtwt) 
    else: continue 
    
    
    #print ievt
    
fout.Write()
fout.Close()
