#!/usr/bin/env python2

from ROOT import gROOT, TFile, TCanvas, TFractionFitter, TObjArray, TLegend, TPaveText, kRed, kBlue, kWhite, kBlack
import ROOT
import numpy as np
import ctypes

filename = "plots.root"

# bins = [x for x in xrange(0, 103, 2.5)]
npBins = np.linspace(0, 100, 40)
bins = [round(x * 2) / 2 for x in npBins]
# region = "Bar_Pt50To75"
# region = "Bar_2Jet_Pt50To75"
# region = "End1_Pt350To400"
# region = "End1_2Jet_Pt350To400"
region = "Bar_Pt120To140"
# region = "Bar_Pt175To200"
# region = "Bar_Pt250To300"
# region = "Bar_Pt300To350"
# region = "Bar_Pt350To400"
# region = "End1_Pt120To140"

gROOT.SetBatch(True)

tfile = TFile.Open(filename)
dataHist = tfile.Get("TrkIso_Data/DataTrkIso_"+region)
eleHist = tfile.Get("TrkIso_Electrons/ElesTrkIso_"+region)
jetsHist = tfile.Get("TrkIso_Jets/JetsTrkIso_"+region)
dyElesHist = tfile.Get("TrkIso_DYElectrons/DYElesTrkIso_"+region)

c = TCanvas()
c.cd()
c.SetLogy()

dataHist = dataHist.Rebin(len(bins)-1, dataHist.GetName(), np.array(bins))
eleHist = eleHist.Rebin(len(bins)-1, eleHist.GetName(), np.array(bins))
jetsHist = jetsHist.Rebin(len(bins)-1, jetsHist.GetName(), np.array(bins))
dyElesHist = dyElesHist.Rebin(len(bins)-1, dyElesHist.GetName(), np.array(bins))
dataHist.Draw()
dataHist.SetStats(0)
dataHist.GetXaxis().SetRangeUser(0, 50)
dataHist.GetXaxis().SetNdivisions(515)
dataHist.GetXaxis().SetTitle("Trk p_{T} Iso. [GeV]")
dataHist.SetMinimum(10)
dataHist.SetMaximum(10*dataHist.GetMaximum())
dataHist.Draw()
jetsHist.SetLineColor(kRed+1)
jetsHist.SetLineStyle(3)
jetsHist.SetLineWidth(2)
jetsHist.SetMarkerColor(kRed+1)
jetsHist.SetMarkerSize(1e-100)
eleHist.SetLineColor(kBlue)
eleHist.SetMarkerColor(kBlue)
eleHist.Draw("same")

# MC amc@NLO can have negative bins
for iBin in xrange(0, dyElesHist.GetNbinsX()+1):
    if dyElesHist.GetBinContent(iBin) < 0:
        print "INFO: Found negative value in dyElesHist bin:", iBin, "; set to zero"
        dyElesHist.SetBinContent(iBin, 0)
dyElesHist.SetLineColor(kBlue)
dyElesHist.SetMarkerColor(kBlue)
dyElesHist.SetMarkerSize(1e-100)
dyElesHist.SetLineStyle(3)
dyElesHist.SetLineWidth(2)

# fraction fitting part
myObjArr = TObjArray(2)
myObjArr.Add(dyElesHist)
myObjArr.Add(jetsHist)
fracFit = TFractionFitter(eleHist, myObjArr)
fracFit.Constrain(0, 0.0, 1.0)
fracFit.Constrain(1, 0.0, 1.0)
fracFit.GetFitter().Config().ParSettings(0).Set("DYElectrons", 0.6, 0.01, 0, 1.0)
fracFit.GetFitter().Config().ParSettings(1).Set("Jets", 0.4, 0.01, 0, 1.0)
fracFit.SetRangeX(1, 20)
status = fracFit.Fit()
print "fit status=", status
fitPlot = fracFit.GetPlot()
fitPlot.SetLineColor(kBlack)
fitPlot.SetMarkerColor(kBlack)
fitPlot.SetLineStyle(3)
fitPlot.SetLineWidth(2)
fitPlot.Draw("esamehist")
ROOT.SetOwnership(fracFit, False)
eleFrac = ctypes.c_double()
eleFracErr = ctypes.c_double()
fracFit.GetResult(0, eleFrac, eleFracErr)
dyElesHist.Scale(eleFrac.value*eleHist.Integral()/dyElesHist.Integral())
jetsFrac = ctypes.c_double()
jetsFracErr = ctypes.c_double()
fracFit.GetResult(1, jetsFrac, jetsFracErr)
jetsHist.Scale(jetsFrac.value*eleHist.Integral()/jetsHist.Integral())

dyElesHist.Draw("esamehist")
jetsHist.Draw("esamehist")

leg = TLegend(0.62, 0.73, 0.87, 0.90)
leg.SetFillColor(kWhite)
leg.SetFillStyle(1001)
leg.SetBorderSize(0)
leg.SetShadowColor(10)
leg.SetMargin(0.2)
leg.SetTextFont(132)
leg.AddEntry(dataHist, "denominator (loose electrons)", "lp")
leg.AddEntry(eleHist, "data (HEEP' electrons)", "lp")
leg.AddEntry(jetsHist, "jet template", "l")
leg.AddEntry(dyElesHist, "ele template (DY MC)", "l")
leg.AddEntry(fitPlot, "ele+jets template", "l")
leg.Draw("same")

text = TPaveText(0.47, 0.78, 0.57, 0.87, "NDC")
text.SetBorderSize(0)
text.SetTextColor(kBlack)
text.SetTextSize(3e-2)
regSplit = region.split("_")
det = regSplit[0]
if det == "Bar":
    det = "barrel"
ptSplit = regSplit[1].split("To")
ptMin = ptSplit[0][2:]
ptMax = ptSplit[1]
text.AddText(ptMin+" < P_{T} < "+ptMax+" GeV")
text.AddText(det)
text.Draw()

c.Print("trkIsoOverlay_"+region+".pdf")

tfile.Close()
