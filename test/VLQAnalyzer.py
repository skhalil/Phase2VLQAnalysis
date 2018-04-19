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
fout = TFile(options.files.rstrip().replace('txt', 'root'), 'RECREATE')
fout.cd()

hCutflow = TH1D("hCutflow" ,";;Events;" ,10 ,0.5 ,10.5 )
hCutflow.GetXaxis().SetBinLabel(1 ,"All evts") ; 

# Open the input ntuples
fnames = [line.strip() for line in open(options.files, 'r')]

# Begin running over all trees
ievt = 0
for fname in fnames:
  if maxEvts > 0 and ievt > maxEvts: break
  if ievt%100 == 0: print " Processing evt %i" % ievt

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
    if vertices == 0: continue
    hCutflow.Fill(2, evtwt)

  # require exactly one medium WP electron or muon in event


  # loop on electrons, selecting medium WP electrons
    nmedium=0
    for i in range(0, len(t.Electons_pt)):
      eWP = t.Electons_eleWP[i]
      print 'eWP =', eWP
      # choose the medium WP: 1 = loose, 2 = medium, 4 = tight
      if (eWP & 2) != 2: continue
      nmedium += 1
       
      #if (eWP & 1) == 1 : loose_ele = 1 
      #if (eWP & 2) == 2 : medium_ele = 1
      #if (eWP & 4) == 4 : tight_ele = 1
      
      #print 'loose = ', loose_ele, 'medium = ', medium_ele, 'tight = ', tight_ele
    #if not tight_ele: 
    
    nET = t.Electons_nMedium
    print 'nMedium = ', ntight, 'nMedium = ', nET
    #print ievt
