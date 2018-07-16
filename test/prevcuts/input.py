import sys, os
from ROOT import TH1D, TFile, TLegend

# Inputs

gSF = 0.97 #trigger SF
path = '/home/t3-ku/z184o935/CMSSW_9_3_2/src/Upgrades/VLQAnalyzer/test/prevcuts/'

f_tt_M2T4 = TFile(path+'TT_M2T4_out.root')
f_WToLNu_0J = TFile(path+'WToLNu_0J_out.root')
f_WToLNu_1J = TFile(path+'WToLNu_1J_out.root')
f_WToLNu_2J = TFile(path+'WToLNu_2J_out.root')
f_WToLNu_3J = TFile(path+'WToLNu_3J_out.root')
f_DYToLL_0J = TFile(path+'DYToLL_0J_out.root')
f_DYToLL_1J = TFile(path+'DYToLL_1J_out.root')
f_DYToLL_2J = TFile(path+'DYToLL_2J_out.root')
f_DYToLL_3J = TFile(path+'DYToLL_3J_out.root')
f_st_tch_antitop = TFile(path+'ST_tch_antitop_out.root')
f_st_tch_top = TFile(path+'ST_tch_top_out.root')
f_st_tW_DR_antitop = TFile(path+'ST_tW_antitop_out.root')
f_st_tW_DR_top = TFile(path+'ST_tW_top_out.root')
f_Tbj_M1000_W10 = TFile(path+'T_M1000_W10_out.root')
f_Tbj_M1000_W20 = TFile(path+'T_M1000_W20_out.root')
f_Tbj_M1000_W30 = TFile(path+'T_M1000_W30_out.root')
f_Tbj_M1500_W10 = TFile(path+'T_M1500_W10_out.root')
f_Tbj_M1500_W20 = TFile(path+'T_M1500_W20_out.root')
f_Tbj_M1500_W30 = TFile(path+'T_M1500_W30_out.root')
f_Tbj_M2000_W10 = TFile(path+'T_M2000_W10_out.root')
f_Tbj_M2000_W20 = TFile(path+'T_M2000_W20_out.root')
f_Tbj_M2000_W30 = TFile(path+'T_M2000_W30_out.root')
f_Tbj_M2500_W10 = TFile(path+'T_M2500_W10_out.root')
f_Tbj_M2500_W20 = TFile(path+'T_M2500_W20_out.root')
f_Tbj_M2500_W30 = TFile(path+'T_M2500_W30_out.root')
f_Tbj_M3000_W10 = TFile(path+'T_M3000_W10_out.root')
f_Tbj_M3000_W20 = TFile(path+'T_M3000_W20_out.root')
f_Tbj_M3000_W30 = TFile(path+'T_M3000_W30_out.root')

# cross sections

tt_M2T4_xs = 864.4 * gSF
WToLNu_0J_xs = 38870.0 * gSF
WToLNu_1J_xs = 10330.0 * gSF
WToLNu_2J_xs = 3314.0 * gSF
WToLNu_3J_xs = 1891.0 * gSF
DYToLL_0J_xs = 3668.0 * gSF
DYToLL_1J_xs = 1094.0 * gSF
DYToLL_2J_xs = 369.7 * gSF
DYToLL_3J_xs = 190.2 * gSF
st_tch_antitop_xs = 93.28
st_tch_top_xs = 154.76
st_tW_DR_antitop_xs = 42.2
st_tW_DR_top_xs = 42.2
Tbj_M1000_W10_xs = 1.0 * gSF * 0.58
Tbj_M1000_W20_xs = 1.0 * gSF * 0.58
Tbj_M1000_W30_xs = 1.0 * gSF * 0.58
Tbj_M1500_W10_xs = 1.0 * gSF * 0.58
Tbj_M1500_W20_xs = 1.0 * gSF * 0.58
Tbj_M1500_W30_xs = 1.0 * gSF * 0.58
Tbj_M2000_W10_xs = 1.0 * gSF * 0.58
Tbj_M2000_W20_xs = 1.0 * gSF * 0.58
Tbj_M2000_W30_xs = 1.0 * gSF * 0.58
Tbj_M2500_W10_xs = 1.0 * gSF * 0.58
Tbj_M2500_W20_xs = 1.0 * gSF * 0.58
Tbj_M2500_W30_xs = 1.0 * gSF * 0.58
Tbj_M3000_W10_xs = 1.0 * gSF * 0.58
Tbj_M3000_W20_xs = 1.0 * gSF * 0.58
Tbj_M3000_W30_xs = 1.0 * gSF * 0.58

# num events

tt_M2T4_num = f_tt_M2T4.Get('hNGenEvents').GetBinContent(1)
WToLNu_0J_num = f_WToLNu_0J.Get('hNGenEvents').GetBinContent(1)
WToLNu_1J_num = f_WToLNu_1J.Get('hNGenEvents').GetBinContent(1)
WToLNu_2J_num = f_WToLNu_2J.Get('hNGenEvents').GetBinContent(1)
WToLNu_3J_num = f_WToLNu_3J.Get('hNGenEvents').GetBinContent(1)
DYToLL_0J_num = f_DYToLL_0J.Get('hNGenEvents').GetBinContent(1)
DYToLL_1J_num = f_DYToLL_1J.Get('hNGenEvents').GetBinContent(1)
DYToLL_2J_num = f_DYToLL_2J.Get('hNGenEvents').GetBinContent(1)
DYToLL_3J_num = f_DYToLL_3J.Get('hNGenEvents').GetBinContent(1)
st_tch_antitop_num = f_st_tch_antitop.Get('hNGenEvents').GetBinContent(1)
st_tch_top_num = f_st_tch_top.Get('hNGenEvents').GetBinContent(1)
st_tW_DR_antitop_num = f_st_tW_DR_antitop.Get('hNGenEvents').GetBinContent(1)
st_tW_DR_top_num = f_st_tW_DR_top.Get('hNGenEvents').GetBinContent(1)
Tbj_M1000_W10_num = f_Tbj_M1000_W10.Get('hNGenEvents').GetBinContent(1)
Tbj_M1000_W20_num = f_Tbj_M1000_W20.Get('hNGenEvents').GetBinContent(1)
Tbj_M1000_W30_num = f_Tbj_M1000_W30.Get('hNGenEvents').GetBinContent(1)
Tbj_M1500_W10_num = f_Tbj_M1500_W10.Get('hNGenEvents').GetBinContent(1)
Tbj_M1500_W20_num = f_Tbj_M1500_W20.Get('hNGenEvents').GetBinContent(1)
Tbj_M1500_W30_num = f_Tbj_M1500_W30.Get('hNGenEvents').GetBinContent(1)
Tbj_M2000_W10_num = f_Tbj_M2000_W10.Get('hNGenEvents').GetBinContent(1)
Tbj_M2000_W20_num = f_Tbj_M2000_W20.Get('hNGenEvents').GetBinContent(1)
Tbj_M2000_W30_num = f_Tbj_M2000_W30.Get('hNGenEvents').GetBinContent(1)
Tbj_M2500_W10_num = f_Tbj_M2500_W10.Get('hNGenEvents').GetBinContent(1)
Tbj_M2500_W20_num = f_Tbj_M2500_W20.Get('hNGenEvents').GetBinContent(1)
Tbj_M2500_W30_num = f_Tbj_M2500_W30.Get('hNGenEvents').GetBinContent(1)
Tbj_M3000_W10_num = f_Tbj_M3000_W10.Get('hNGenEvents').GetBinContent(1)
Tbj_M3000_W20_num = f_Tbj_M3000_W20.Get('hNGenEvents').GetBinContent(1)
Tbj_M3000_W30_num = f_Tbj_M3000_W30.Get('hNGenEvents').GetBinContent(1)

tt_M2T4nums = []
WToLNunums = []
DYToLLnums = []
stnums = []
T_M1000_W10nums = []
T_M1500_W10nums = []
T_M2000_W10nums = []
T_M2500_W10nums = []
T_M3000_W10nums = []

tt_M2T4cutflow = f_tt_M2T4.Get('hCutflow')
WToLNu0cutflow = f_WToLNu_0J.Get('hCutflow')
WToLNu1cutflow = f_WToLNu_1J.Get('hCutflow')
WToLNu2cutflow = f_WToLNu_2J.Get('hCutflow')
WToLNu3cutflow = f_WToLNu_3J.Get('hCutflow')
DYToLL0cutflow = f_DYToLL_0J.Get('hCutflow')
DYToLL1cutflow = f_DYToLL_1J.Get('hCutflow')
DYToLL2cutflow = f_DYToLL_2J.Get('hCutflow')
DYToLL3cutflow = f_DYToLL_3J.Get('hCutflow')
st_tch_antitopcutflow = f_st_tch_antitop.Get('hCutflow')
st_tch_topcutflow = f_st_tch_top.Get('hCutflow')
st_tW_DR_antitopcutflow = f_st_tW_DR_antitop.Get('hCutflow')
st_tW_DR_topcutflow = f_st_tW_DR_top.Get('hCutflow')
T_M1000_W10cutflow = f_Tbj_M1000_W10.Get('hCutflow')
T_M1500_W10cutflow = f_Tbj_M1500_W10.Get('hCutflow')
T_M2000_W10cutflow = f_Tbj_M2000_W10.Get('hCutflow')
T_M2500_W10cutflow = f_Tbj_M2500_W10.Get('hCutflow')
T_M3000_W10cutflow = f_Tbj_M3000_W10.Get('hCutflow')

for i in range(0,tt_M2T4cutflow.GetSize()):
    tt_M2T4nums.append(tt_M2T4cutflow.GetBinContent(i))
for i in range(0,WToLNu0cutflow.GetSize()):
    WToLNunums.append(WToLNu0cutflow.GetBinContent(i) + WToLNu1cutflow.GetBinContent(i) + WToLNu2cutflow.GetBinContent(i) + WToLNu3cutflow.GetBinContent(i))
for i in range(0,DYToLL0cutflow.GetSize()):
    DYToLLnums.append(DYToLL0cutflow.GetBinContent(i) + DYToLL1cutflow.GetBinContent(i) + DYToLL2cutflow.GetBinContent(i) + DYToLL3cutflow.GetBinContent(i))
for i in range(0,st_tch_antitopcutflow.GetSize()):
    stnums.append(st_tch_antitopcutflow.GetBinContent(i) + st_tch_topcutflow.GetBinContent(i) + st_tW_DR_antitopcutflow.GetBinContent(i) + st_tW_DR_topcutflow.GetBinContent(i))
for i in range(0,T_M1000_W10cutflow.GetSize()):
    T_M1000_W10nums.append(T_M1000_W10cutflow.GetBinContent(i))
for i in range(0,T_M1500_W10cutflow.GetSize()):
    T_M1500_W10nums.append(T_M1500_W10cutflow.GetBinContent(i))
for i in range(0,T_M2000_W10cutflow.GetSize()):
    T_M2000_W10nums.append(T_M2000_W10cutflow.GetBinContent(i))
for i in range(0,T_M2500_W10cutflow.GetSize()):
    T_M2500_W10nums.append(T_M2500_W10cutflow.GetBinContent(i))
for i in range(0, T_M3000_W10cutflow.GetSize()):
    T_M3000_W10nums.append(T_M3000_W10cutflow.GetBinContent(i))

# variables

leg = TLegend(0.80,0.99,0.99,0.8) 
leg.SetBorderSize(0)
leg.SetFillColor(10)
leg.SetLineColor(10)
leg.SetLineWidth(0)

# functions

def overUnderFlow(hist):
    xbins = hist.GetNbinsX()
    hist.SetBinContent(xbins, hist.GetBinContent(xbins)+hist.GetBinContent(xbins+1))
    hist.SetBinContent(1, hist.GetBinContent(0)+hist.GetBinContent(1))
    hist.SetBinError(xbins, TMath.Sqrt(TMath.Power(hist.GetBinError(xbins),2)+TMath.Power(hist.GetBinError(xbins+1),2)))
    hist.SetBinError(1, TMath.Sqrt(TMath.Power(hist.GetBinError(0),2)+TMath.Power(hist.GetBinError(1),2)))
    hist.SetBinContent(xbins+1, 0.)
    hist.SetBinContent(0, 0.)
    hist.SetBinError(xbins+1, 0.)
    hist.SetBinError(0, 0.)

def setCosmetics(hist, legname, hname, color):
    hist.Rebin(rebinS)
    hist.SetLineColor(color)
    hist.SetName(hname)
    hist.SetTitle("")
    if 'Tbj' in hname:
        hist.SetLineWidth(2)
        leg.AddEntry(hist, legname, 'l')
    else:
        hist.SetFillColor(color)
        leg.AddEntry(hist, legname, 'f')

def setTitle(hs,xTitle,yTitle):
    y = hs.GetYaxis()
    x = hs.GetXaxis()
    #y.SetTitle("Events / Bin")
    y.SetTitle(yTitle)
    x.SetTitle(xTitle)
    y.SetLabelSize(0.04)
    y.SetTitleSize(0.06)
    y.SetTitleOffset(0.75)
    y.SetTitleFont(42)
    x.SetTitleSize(0.04)
    x.SetTitleFont(42)

def getHisto( label, leg, var, Samples, color, verbose) :
    histos = []
    for iSample in Samples :
        ifile = iSample[0]
        xs = iSample[1]
        nevt = iSample[2]
        lumi = iSample[3]
        readname = var
        hist  = ifile.Get( readname ).Clone()
        if verbose:
            print 'file: {0:<20}, histo:{1:<10}, integral before weight:{2:<3.3f}, nEntries:{3:<3.0f}, weight:{4:<2.3f}'.format(
                ifile.GetName(),
                hist.GetName(),
                hist.Integral(), hist.GetEntries(), xs * lumi /nevt
                )
        hist.Sumw2()
        hist.Scale( xs * lumi /nevt)
        #hist.SetName(leg+hist.GetName())
        histos.append( hist )

    histo = histos[0]
    setCosmetics(histo, leg, label+var, color)
    for ihisto in range(1, len(histos) ):
        #print 'ihisto =', ihisto, 'integral', histos[ihisto].Integral(), ', entries', histos[ihisto].GetEntries()
        histo.Add( histos[ihisto] )
        #print 'after addition', histo.Integral()
    if verbose:
        print 'newName: {0:<5}, Entries:{1:5.2f},  newIntegral: {2:5.2f}'.format(label+var, histo.GetEntries(), histo.Integral() )
    return histo
