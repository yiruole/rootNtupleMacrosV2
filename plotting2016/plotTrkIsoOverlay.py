#!/usr/bin/env python2

from ROOT import gROOT, gPad, gStyle, TFile, TCanvas, TPad, TFractionFitter, TObjArray, TLegend, TPaveText, kRed, kBlue, kWhite, kBlack
import ROOT
import numpy as np
import ctypes
import copy
import os

# configurables
filename = "plots.root"
makeRatioPlot = True
plotsDir = "trkIsoPlots"

# start running
npBins = np.linspace(0, 100, 40)
bins = [round(x * 2) / 2 for x in npBins]

regList = ["Bar_Pt120To140", "Bar_Pt140To175", "Bar_Pt175To200", "Bar_Pt350To400"]
regList.extend(["End1_Pt120To140", "End1_Pt140To175", "End1_Pt175To200", "End1_Pt350To400"])
regList.extend(["End2_Pt120To140", "End2_Pt140To175", "End2_Pt175To200", "End2_Pt350To400"])

gROOT.SetBatch(True)

# gStyle.SetPadTopMargin(0.08)
# gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadRightMargin(0.02)
gStyle.SetPadLeftMargin(0.1)

if not os.path.isdir(plotsDir):
    os.mkdir(plotsDir)

tfile = TFile.Open(filename)

for region in regList:
    denomHist = tfile.Get("TrkIso_Data/DataTrkIso_"+region)
    heepPrimeHist = tfile.Get("TrkIso_Electrons/ElesTrkIso_"+region)
    jetsHist = tfile.Get("TrkIso_Jets/JetsTrkIso_"+region)
    dyElesHist = tfile.Get("TrkIso_DYElectrons/DYElesTrkIso_"+region)
    # rebin to 2.5 GeV bins
    denomHist = denomHist.Rebin(len(bins)-1, denomHist.GetName(), np.array(bins))
    heepPrimeHist = heepPrimeHist.Rebin(len(bins)-1, heepPrimeHist.GetName(), np.array(bins))
    jetsHist = jetsHist.Rebin(len(bins)-1, jetsHist.GetName(), np.array(bins))
    dyElesHist = dyElesHist.Rebin(len(bins)-1, dyElesHist.GetName(), np.array(bins))

    canvas = TCanvas()
    if makeRatioPlot:
        ratioPadFraction = 0.3
        fPads1 = TPad("pad1", "", 0.00, ratioPadFraction, 0.99, 0.99)
        fPads2 = TPad("pad2", "", 0.00, 0.00, 0.99, ratioPadFraction)
        fPads1.SetFillColor(0)
        fPads1.SetLineColor(0)
        fPads2.SetFillColor(0)
        fPads2.SetLineColor(0)
        fPads1.SetBottomMargin(1e-2)
        fPads2.SetTopMargin(3e-2)
        fPads2.SetBottomMargin(3e-1)
        fPads1.Draw()
        fPads2.Draw()
    else:
        fPads1 = TPad("pad1", "", 0.00, 0.0, 0.99, 0.99)
        fPads1.SetFillColor(0)
        fPads1.SetLineColor(0)
        fPads1.Draw()

    fPads1.cd()
    fPads1.SetLogy()
    # styles
    denomHist.Draw()
    denomHist.SetStats(0)
    denomHist.GetXaxis().SetRangeUser(0, 50)
    denomHist.GetXaxis().SetNdivisions(515)
    denomHist.GetXaxis().SetTitle("Trk p_{T} Iso. [GeV]")
    denomHist.SetMinimum(10)
    denomHist.SetMaximum(10*denomHist.GetMaximum())
    denomHist.Draw()
    # jets hist style
    jetsHist.SetLineColor(kRed+1)
    jetsHist.SetLineStyle(3)
    jetsHist.SetLineWidth(2)
    jetsHist.SetMarkerColor(kRed+1)
    jetsHist.SetMarkerSize(1e-100)
    # eles hist style
    heepPrimeHist.SetLineColor(kBlue)
    heepPrimeHist.SetMarkerColor(kBlue)
    heepPrimeHist.Draw("same")

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
    fracFit = TFractionFitter(heepPrimeHist, myObjArr)
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
    dyElesHist.Scale(eleFrac.value*heepPrimeHist.Integral()/dyElesHist.Integral())
    jetsFrac = ctypes.c_double()
    jetsFracErr = ctypes.c_double()
    fracFit.GetResult(1, jetsFrac, jetsFracErr)
    jetsHist.Scale(jetsFrac.value*heepPrimeHist.Integral()/jetsHist.Integral())

    dyElesHist.Draw("esamehist")
    jetsHist.Draw("esamehist")

    leg = TLegend(0.75, 0.7, 0.98, 0.87)
    leg.SetFillColor(kWhite)
    leg.SetFillStyle(1001)
    leg.SetBorderSize(0)
    leg.SetShadowColor(10)
    leg.SetMargin(0.2)
    leg.SetTextFont(132)
    leg.AddEntry(denomHist, "denominator (loose electrons)", "lp")
    leg.AddEntry(heepPrimeHist, "HEEP' electrons", "lp")
    leg.AddEntry(jetsHist, "jet template", "l")
    leg.AddEntry(dyElesHist, "ele template (DY MC)", "l")
    leg.AddEntry(fitPlot, "ele+jets template", "l")
    leg.Draw("same")

    text = TPaveText(0.57, 0.78, 0.67, 0.87, "NDC")
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

    canvas.Update()
    gPad.RedrawAxis()
    gPad.Modified()

    # ratio plot
    h_ratio1 = copy.deepcopy(heepPrimeHist)
    h_ratio1.SetStats(0)
    if makeRatioPlot == 1:
        denomHist.GetXaxis().SetTitle()
        fPads2.cd()
        fPads2.SetGridy()
        h_ratio1.Divide(fitPlot)

        h_ratio1.GetXaxis().SetTitle("Trk p_{T} Iso. [GeV]")
        h_ratio1.GetXaxis().SetNdivisions(515)
        h_ratio1.GetXaxis().SetTitleSize(0.12)
        h_ratio1.GetXaxis().SetLabelSize(0.12)
        h_ratio1.GetXaxis().SetRangeUser(0.0, 50)
        h_ratio1.GetYaxis().SetRangeUser(0.0, 2)
        h_ratio1.GetYaxis().SetTitle("HEEP'/tot. templ.")
        h_ratio1.GetYaxis().SetLabelSize(0.1)
        h_ratio1.GetYaxis().SetTitleSize(0.1)
        h_ratio1.GetYaxis().SetTitleOffset(0.3)
        h_ratio1.Draw("e0p")

        canvas.Update()
        gPad.RedrawAxis()
        gPad.Modified()

    canvas.Print(plotsDir+"/trkIsoOverlay_"+region+".pdf")

tfile.Close()
