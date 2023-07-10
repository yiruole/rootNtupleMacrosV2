#!/usr/bin/env python2

from __future__ import print_function

from ROOT import (
    gROOT,
    gPad,
    gStyle,
    TFile,
    TCanvas,
    TPad,
    TFractionFitter,
    TObjArray,
    TLegend,
    TPaveText,
    kRed,
    kBlue,
    kWhite,
    kBlack,
    kGreen,
    kAzure,
    kGray,
    TH1D,
    TEfficiency,
    TGraphAsymmErrors,
)
import ROOT
import numpy as np
import ctypes
import copy
import os
import math

from numpy.core.numerictypes import LOWER_TABLE


def GetCPErrorLow(passing, total):
    # use unweighted pois calc from TGraphAsymmErrors::Divide code
    # https://root.cern.ch/doc/master/TGraphAsymmErrors_8cxx_source.html#l00582
    low = TEfficiency.ClopperPearson(passing + total, passing, 0.683, False)
    eff = passing / (passing + total)
    ratio = eff / (1.0 - eff)
    low = low / (1.0 - low)
    eff = ratio
    return eff - low


# configurables
filename = "$LQDATA/2016/qcdFakeRateCalc/calcFR_2016HEEPpre/fakeRate_plots.root"
makeRatioPlot = False
doHEEPFR = True
if doHEEPFR:
    varName = "TrkIsoHEEP7"
    varTitle = "Trk p_{T} Iso. [GeV]"
    varXRangeMax = 50
    varSBMinLowPt = 10
    varSBMaxLowPt = 20
    varSBMinPtGt500 = varSBMinLowPt
    varSBMaxPtGt500 = varSBMaxLowPt
    varSRMax = 5
    fractionFitXRangeMax = varXRangeMax
else:
    varName = "PFRelIsoAll"
    varTitle = "PFRelIso"
    varXRangeMax = 1.5
    varSBMinLowPt = 0.25
    varSBMaxLowPt = 0.45
    # varSBMin = 0.3
    # varSBMax = 0.5
    varSBMinPtGt500 = 0.15
    varSBMaxPtGt500 = 0.35
    varSRMax = 0.112
    # varName = "Full5x5SigmaIEtaIEta"
    # varTitle = "Full5x5SigmaIEtaIEta"
    # varXRangeMax = 0.1
    fractionFitXRangeMax = 0.8
plotsDir = "$LQDATA/2016/qcdFakeRateCalc/calcFR_2016HEEPpre/ExpandedRangeIsoPlots"
doFracFit = False
do2016 = True

# start running
npBins = np.linspace(0, 100, 40)
bins = [round(x * 2) / 2 for x in npBins]

if do2016:
    # ptBins = ["Pt120To140", "Pt140To175", "Pt175To200", "Pt350To400"]
    ptBinsEndcap = [
        36,
        50,
        75,
        90,
        120,
        140,
        175,
        200,
        225,
        250,
        300,
        350,
        400,
        500,
        600,
        1000,
    ]
    ptBinsHighEndcap = [
        36,
        50,
        75,
        90,
        120,
        140,
        175,
        200,
        225,
        250,
        300,
        350,
        400,
        500,
        1000,
    ]
else:
    # ptBins = ["Pt120To150", "Pt150To175", "Pt175To200", "Pt350To400"]
    ptBinsEndcap = [
        36,
        50,
        75,
        90,
        120,
        150,
        175,
        200,
        225,
        250,
        300,
        350,
        400,
        500,
        600,
        1000,
    ]
    ptBinsHighEndcap = [
        36,
        50,
        75,
        90,
        120,
        150,
        175,
        200,
        225,
        250,
        300,
        350,
        400,
        500,
        1000,
    ]
# TODO HEM1516Only
ptBinsEndcapHEM1516Only = [
    36,
    50,
    75,
    90,
    120,
    150,
    175,
    200,
    225,
    250,
    300,
    350,
    400,
    500,
    1000,
]
ptBins = [
    "Pt" + str(lowEnd) + "To" + str(ptBinsEndcap[ptBinsEndcap.index(lowEnd) + 1])
    for lowEnd in ptBinsEndcap[:-1]
]

ptBinsEnd2 = [
    "Pt" + str(lowEnd) + "To" + str(ptBinsHighEndcap[ptBinsHighEndcap.index(lowEnd) + 1])
    for lowEnd in ptBinsHighEndcap[:-1]
]
etaRegions = ["Bar", "End1", "End2"]
regList = []
for reg in etaRegions:
    if reg == "End2":
        regList.extend([reg + "_" + x for x in ptBinsEnd2])
    else:
        regList.extend([reg + "_" + x for x in ptBins])
dyEleFracInSBByReg = {}
dyEleFracErrInSBByReg = {}
dyEleInSBHistsByEtaReg = {}
dyEleTotalHistsByEtaReg = {}
mcEleInSBHistsByEtaReg = {}
mcEleTotalHistsByEtaReg = {}
for etaReg in etaRegions:
    dyEleInSBHistsByEtaReg[etaReg] = TH1D(
        "dyEleInSB_" + etaReg,
        "",
        len(ptBinsEndcap) - 1,
        np.array(ptBinsEndcap, dtype=float),
    )
    dyEleTotalHistsByEtaReg[etaReg] = TH1D(
        "dyEleTotal_" + etaReg,
        "",
        len(ptBinsEndcap) - 1,
        np.array(ptBinsEndcap, dtype=float),
    )
    mcEleInSBHistsByEtaReg[etaReg] = TH1D(
        "mcEleInSB_" + etaReg,
        "",
        len(ptBinsEndcap) - 1,
        np.array(ptBinsEndcap, dtype=float),
    )
    mcEleTotalHistsByEtaReg[etaReg] = TH1D(
        "mcEleTotal_" + etaReg,
        "",
        len(ptBinsEndcap) - 1,
        np.array(ptBinsEndcap, dtype=float),
    )

gROOT.SetBatch(True)

# gStyle.SetPadTopMargin(0.08)
# gStyle.SetPadBottomMargin(0.12)
gStyle.SetPadRightMargin(0.02)
gStyle.SetPadLeftMargin(0.1)

if not os.path.isdir(plotsDir):
    os.mkdir(plotsDir)

tfile = TFile.Open(filename)

for region in regList:
    denomHist = tfile.Get(varName + "_Data/Data" + varName + "_" + region)
    heepPrimeHist = tfile.Get(varName + "_Electrons/Eles" + varName + "_" + region)
    jetsHist = tfile.Get(varName + "_Jets/Jets" + varName + "_" + region)
    dyElesHist = tfile.Get(varName + "_DYElectrons/DYEles" + varName + "_" + region)
    mcElesHist = tfile.Get(varName + "_MCElectrons/" + varName + "MC_" + region)

    if denomHist == 0 or denomHist is None:
        raise RuntimeError(
            "Did not find hist named {} in file {}".format(
                varName + "_Data/Data" + varName + "_" + region, tfile.GetName()
            )
        )
    if dyElesHist == 0 or dyElesHist is None:
        raise RuntimeError(
            "Did not find hist named {} in file {}".format(
                varName + "_DYElectrons/DYEles" + varName + "_" + region,
                tfile.GetName(),
            )
        )
    if mcElesHist == 0 or mcElesHist is None:
        raise RuntimeError(
            "Did not find hist named {} in file {}".format(
                varName + "_MCElectrons/Eles" + varName + "_" + region,
                tfile.GetName(),
            )
        )
    print()
    print(
        "Run for this region: found hist named {} in file {}".format(
            varName + "_Data/Data" + varName + "_" + region, tfile.GetName()
        )
    )
    # rebin to 2.5 GeV bins
    if "TrkIso" in varName:
        denomHist = denomHist.Rebin(len(bins) - 1, denomHist.GetName(), np.array(bins))
        heepPrimeHist = heepPrimeHist.Rebin(
            len(bins) - 1, heepPrimeHist.GetName(), np.array(bins)
        )
        jetsHist = jetsHist.Rebin(len(bins) - 1, jetsHist.GetName(), np.array(bins))
        dyElesHist = dyElesHist.Rebin(
            len(bins) - 1, dyElesHist.GetName(), np.array(bins)
        )
        mcElesHist = mcElesHist.Rebin(
            len(bins) - 1, mcElesHist.GetName(), np.array(bins)
        )
    # else:
    #     denomHist = denomHist.Rebin(4)
    #     heepPrimeHist = heepPrimeHist.Rebin(4)
    #     jetsHist = jetsHist.Rebin(4)
    #     dyElesHist = dyElesHist.Rebin(4)
    #     mcElesHist = mcElesHist.Rebin(4)

    canvas = TCanvas()
    if makeRatioPlot:
        ratioPadFraction = 0.3
        fPads1 = TPad("pad1", "", 0.00, ratioPadFraction, 0.99, 0.99)
        fPads2 = TPad("pad2", "", 0.00, 0.00, 0.99, ratioPadFraction + 0.001)
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
    denomHist.SetLineColor(kGreen + 2)
    denomHist.SetMarkerColor(kGreen + 2)
    denomHist.SetStats(0)
    denomHist.GetXaxis().SetRangeUser(0, varXRangeMax)
    denomHist.GetXaxis().SetNdivisions(515)
    denomHist.GetXaxis().SetTitle(varTitle)
    denomHist.SetMinimum(0.01)
    denomHist.SetMaximum(10 * denomHist.GetMaximum())
    denomHist.Draw()
    # jets hist style
    jetsHist.SetLineColor(kRed + 1)
    jetsHist.SetLineStyle(3)
    jetsHist.SetLineWidth(2)
    jetsHist.SetMarkerColor(kRed + 1)
    jetsHist.SetMarkerSize(1e-100)
    # eles hist style
    heepPrimeHist.SetLineColor(kAzure + 10)
    heepPrimeHist.SetMarkerColor(kAzure + 10)
    heepPrimeHist.Draw("same")

    # MC amc@NLO can have negative bins
    # for iBin in range(0, dyElesHist.GetNbinsX() + 1):
    #     if dyElesHist.GetBinContent(iBin) < 0:
    #         # print("INFO: Found negative value in dyElesHist bin:", iBin, "; set to zero")
    #         dyElesHist.SetBinContent(iBin, 0)
    # dyElesHist.SetLineColor(kBlue)
    # dyElesHist.SetMarkerColor(kBlue)
    # dyElesHist.SetMarkerSize(1e-100)
    # dyElesHist.SetLineStyle(3)
    # dyElesHist.SetLineWidth(2)
    # if dyElesHist.GetEntries() < 10:
    #     doFracFit = False
    #     print(
    #         "WARN: DY electrons hist had",
    #         dyElesHist.GetEntries(),
    #         "entries; skipping TFractionFitter part.",
    #     )
    for iBin in range(0, mcElesHist.GetNbinsX() + 1):
        if mcElesHist.GetBinContent(iBin) < 0:
            # print("INFO: Found negative value in mcElesHist bin:", iBin, "; set to zero")
            mcElesHist.SetBinContent(iBin, 0)
    mcElesHist.SetLineColor(kBlue)
    mcElesHist.SetMarkerColor(kBlue)
    mcElesHist.SetMarkerSize(1e-100)
    mcElesHist.SetLineStyle(3)
    mcElesHist.SetLineWidth(2)
    if mcElesHist.GetEntries() < 10:
        doFracFit = False
        print(
            "WARN: mc electrons hist had",
            mcElesHist.GetEntries(),
            "entries; skipping TFractionFitter part.",
        )
    etaReg = region.split("_")[0]
    regSplit = region.split("_")
    det = regSplit[0]
    if det == "Bar":
        det = "barrel"
    ptSplit = regSplit[1].split("To")
    ptMin = ptSplit[0][2:]
    ptMax = ptSplit[1]
    ptAvg = (float(ptMin) + float(ptMax)) / 2.0
    if float(ptMin) >= 500:
        varSBMin = 0.15  # varSBMinPtGt500
        varSBMax = 0.2  # varSBMaxPtGt500
        if float(ptMin) >= 600:
            varSBMin = 0.15
            varSBMax = 0.2
    else:
        varSBMin = varSBMinLowPt
        varSBMax = varSBMaxLowPt
    if "bar" in det.lower():
        varSRMax = 0.112 + 0.506 / float(ptMin)
    else:
        varSRMax = 0.108 + 0.963 / float(ptMin)
    # print "First bin contents (before fit):"
    # print "\theepPrime=", heepPrimeHist.GetBinContent(1)
    # print "\tDYEles=", dyElesHist.GetBinContent(1)
    # print "\tJets=", jetsHist.GetBinContent(1)
    # print "\tdenom (FR loose eles)=", denomHist.GetBinContent(1)
    elePrimeSBErr = ctypes.c_double()
    elePrimeSB = heepPrimeHist.IntegralAndError(
        heepPrimeHist.FindFixBin(varSBMin),
        heepPrimeHist.FindFixBin(varSBMax) - 1,
        elePrimeSBErr,
    )
    # print("elePrime: bin for varSBMax={} is {}, binLowEdge={}, binUpEdge={}".format(varSBMax, heepPrimeHist.FindFixBin(varSBMax),
    #       heepPrimeHist.GetXaxis().GetBinLowEdge(heepPrimeHist.FindFixBin(varSBMax)), heepPrimeHist.GetXaxis().GetBinUpEdge(heepPrimeHist.FindFixBin(varSBMax))))
    elePrimeSR = heepPrimeHist.Integral(1, heepPrimeHist.FindFixBin(varSRMax) - 1)
    jetsSBErr = ctypes.c_double()
    jetsSB = jetsHist.IntegralAndError(
        jetsHist.FindFixBin(varSBMin), jetsHist.FindFixBin(varSBMax) - 1, jetsSBErr
    )
    jetsSRErr = ctypes.c_double()
    jetsSR = jetsHist.IntegralAndError(1, jetsHist.FindFixBin(varSRMax) - 1, jetsSRErr)
    dyElesSBErr = ctypes.c_double()
    dyElesSB = dyElesHist.IntegralAndError(
        dyElesHist.FindFixBin(varSBMin), dyElesHist.FindFixBin(varSBMax) - 1, dyElesSBErr
    )
    mcElesSBErr = ctypes.c_double()
    mcElesSB = mcElesHist.IntegralAndError(
        mcElesHist.FindFixBin(varSBMin), mcElesHist.FindFixBin(varSBMax) - 1, mcElesSBErr
    )
    rJets = jetsSR / jetsSB
    nJetsTemplate = rJets * elePrimeSB
    # rJetsErr = (jetsSR/jetsSB)*math.sqrt((jetsSRErr.value/jetsSR)**2+(jetsSBErr.value/jetsSB)**2)
    rJetsErr = GetCPErrorLow(jetsSR, jetsSB)
    nJetsTemplateErr = nJetsTemplate * math.sqrt(
        (rJetsErr / rJets) ** 2 + (elePrimeSBErr.value / elePrimeSB) ** 2
    )
    looseElesErr = ctypes.c_double()
    looseEles = denomHist.IntegralAndError(0, -1, looseElesErr)
    frErr = GetCPErrorLow(nJetsTemplate, looseEles)
    print("Using SB {}-{} and SR 0-{}".format(varSBMin, varSBMax, varSRMax))
    # print(
    #     "dyEles(SB)/dyEles(total)={}/{}={}".format(
    #         elesSB, dyElesHist.Integral(), elesSB / dyElesHist.Integral()
    #     )
    # )
    print(
        "mcEles(SB)/mcEles(total)={}/{}={}".format(
            mcElesSB, mcElesHist.Integral(), mcElesSB / mcElesHist.Integral()
        )
    )
    print("ele'SB = {}".format(elePrimeSB))
    print("ele'SR = {}".format(elePrimeSR))
    print("jetsSB = {}".format(jetsSB))
    print("jetsSR = {}".format(jetsSR))
    print("nJets  = {} +/- {}".format(nJetsTemplate, nJetsTemplateErr))
    print("FR  = {} +/- {}".format(nJetsTemplate / looseEles, frErr))
    dyElesIntegralErr = ctypes.c_double()
    dyElesIntegral = dyElesHist.IntegralAndError(0, -1, dyElesIntegralErr)
    dyEleTotalHistsByEtaReg[etaReg].SetBinContent(
        dyEleTotalHistsByEtaReg[etaReg].FindBin(ptAvg), dyElesIntegral
    )
    dyEleTotalHistsByEtaReg[etaReg].SetBinError(
        dyEleTotalHistsByEtaReg[etaReg].FindBin(ptAvg), dyElesIntegralErr.value
    )
    dyEleInSBHistsByEtaReg[etaReg].SetBinContent(
        dyEleInSBHistsByEtaReg[etaReg].FindBin(ptAvg), dyElesSB
    )
    dyEleInSBHistsByEtaReg[etaReg].SetBinError(
        dyEleInSBHistsByEtaReg[etaReg].FindBin(ptAvg), dyElesSBErr.value
    )
    mcElesIntegralErr = ctypes.c_double()
    mcElesIntegral = mcElesHist.IntegralAndError(0, -1, mcElesIntegralErr)
    mcEleTotalHistsByEtaReg[etaReg].SetBinContent(
        mcEleTotalHistsByEtaReg[etaReg].FindBin(ptAvg), mcElesIntegral
    )
    mcEleTotalHistsByEtaReg[etaReg].SetBinError(
        mcEleTotalHistsByEtaReg[etaReg].FindBin(ptAvg), mcElesIntegralErr.value
    )
    mcEleInSBHistsByEtaReg[etaReg].SetBinContent(
        mcEleInSBHistsByEtaReg[etaReg].FindBin(ptAvg), mcElesSB
    )
    mcEleInSBHistsByEtaReg[etaReg].SetBinError(
        mcEleInSBHistsByEtaReg[etaReg].FindBin(ptAvg), mcElesSBErr.value
    )
    print("etaReg={}, ptAvg={}, mcElesInSB={} +/- {}, mcElesTotal={} +/- {}".format(etaReg, ptAvg, mcElesSB, mcElesSBErr.value, dyElesIntegral, dyElesIntegralErr.value))
    # fraction fitting part
    if doFracFit:
        myObjArr = TObjArray(2)
        # myObjArr.Add(dyElesHist)
        myObjArr.Add(mcElesHist)
        myObjArr.Add(jetsHist)
        fracFit = TFractionFitter(heepPrimeHist, myObjArr)
        fracFit.Constrain(0, 0.0, 1.0)
        fracFit.Constrain(1, 0.0, 1.0)
        fracFit.GetFitter().Config().ParSettings(0).Set(
            "DYElectrons", 0.7, 0.01, 0, 1.0
        )
        fracFit.GetFitter().Config().ParSettings(1).Set("Jets", 0.3, 0.01, 0, 1.0)
        fracFit.SetRangeX(1, heepPrimeHist.FindFixBin(fractionFitXRangeMax))
        status = fracFit.Fit()
        print("fit status=", status)
        fitPlot = fracFit.GetPlot()
        fitPlot.SetLineColor(kGray + 2)
        fitPlot.SetMarkerColor(kGray + 2)
        fitPlot.SetLineStyle(3)
        fitPlot.SetLineWidth(2)
        fitPlot.Draw("esamehist")
        ROOT.SetOwnership(fracFit, False)
        eleFrac = ctypes.c_double()
        eleFracErr = ctypes.c_double()
        fracFit.GetResult(0, eleFrac, eleFracErr)
        # print "scale DYEles by {}*{}/{}={}".format(eleFrac.value, heepPrimeHist.Integral(), dyElesHist.Integral(), eleFrac.value*heepPrimeHist.Integral()/dyElesHist.Integral())
        # dyElesHist.Scale(
        #     eleFrac.value
        #     * heepPrimeHist.Integral(1, heepPrimeHist.FindFixBin(fractionFitXRangeMax))
        #     / dyElesHist.Integral(1, dyElesHist.FindFixBin(fractionFitXRangeMax))
        # )
        mcElesHist.Scale(
            eleFrac.value * heepPrimeHist.Integral(1, heepPrimeHist.FindFixBin(fractionFitXRangeMax)) / mcElesHist.Integral(1, mcElesHist.FindFixBin(fractionFitXRangeMax))
        )
        jetsFrac = ctypes.c_double()
        jetsFracErr = ctypes.c_double()
        fracFit.GetResult(1, jetsFrac, jetsFracErr)
        # print "scale jets by {}*{}/{}={}".format(jetsFrac.value, heepPrimeHist.Integral(), jetsHist.Integral(), jetsFrac.value*heepPrimeHist.Integral()/jetsHist.Integral())
        jetsHist.Scale(
            jetsFrac.value
            * heepPrimeHist.Integral(1, heepPrimeHist.FindFixBin(fractionFitXRangeMax))
            / jetsHist.Integral(1, jetsHist.FindFixBin(fractionFitXRangeMax))
        )

    # dyElesHist.Draw("esamehist")
    mcElesHist.Draw("esamehist")
    jetsHist.Draw("esamehist")
    # print "First bin contents (after fit):"
    # print "\theepPrime=", heepPrimeHist.GetBinContent(1)
    # print "\tDYEles=", dyElesHist.GetBinContent(1)
    # print "\tJets=", jetsHist.GetBinContent(1)
    # print "\tdenom (FR loose eles)=", denomHist.GetBinContent(1)
    jetsSR = jetsHist.IntegralAndError(1, jetsHist.FindBin(varSRMax) - 1, jetsSRErr)
    print("nJets = {} +/- {}".format(jetsSR, jetsSRErr.value))
    frErr = GetCPErrorLow(jetsSR, looseEles)
    print("FR  = {} +/- {}".format(jetsSR / looseEles, frErr))
    # dyElesSR = dyElesHist.Integral(1, dyElesHist.FindBin(varSRMax) - 1)
    # print("dyElesSR = {}".format(dyElesSR))
    # print("jetsSR/ele'SR = {}".format(jetsSR / elePrimeSR))
    # print("dyElesSR/ele'SR = {}".format(dyElesSR / elePrimeSR))
    mcElesSR = mcElesHist.Integral(1, mcElesHist.FindBin(varSRMax) - 1)
    print("mcElesSR = {}".format(mcElesSR))
    print("jetsSR/ele'SR = {}".format(jetsSR / elePrimeSR))
    print("mcElesSR/ele'SR = {}".format(mcElesSR / elePrimeSR))
    frRatioErr = GetCPErrorLow(nJetsTemplate, jetsSR)
    print(
        "NjetsRatio/NjetsFracFit = {} +/- {}".format(nJetsTemplate / jetsSR, frRatioErr)
    )

    leg = TLegend(0.75, 0.7, 0.97, 0.87)
    leg.SetFillColor(kWhite)
    leg.SetFillStyle(1001)
    leg.SetBorderSize(0)
    # leg.SetShadowColor(10)
    leg.SetMargin(0.2)
    leg.SetTextFont(132)
    leg.AddEntry(denomHist, "denominator (loose electrons)", "lp")
    leg.AddEntry(heepPrimeHist, "ele' electrons", "lp")
    leg.AddEntry(jetsHist, "jet template", "l")
    # leg.AddEntry(dyElesHist, "ele template (DY MC)", "l")
    leg.AddEntry(mcElesHist, "ele template (MC)", "l")
    if doFracFit:
        leg.AddEntry(fitPlot, "ele+jets template", "l")
    leg.Draw("same")

    text = TPaveText(0.57, 0.78, 0.67, 0.87, "NDC")
    text.SetBorderSize(0)
    text.SetTextColor(kBlack)
    text.SetTextSize(3e-2)
    text.SetFillStyle(0)
    text.AddText(ptMin + " < P_{T} < " + ptMax + " GeV")
    text.AddText(det)
    text.Draw()

    canvas.Update()
    gPad.RedrawAxis()
    gPad.Modified()

    # ratio plot
    h_ratio1 = copy.deepcopy(heepPrimeHist)
    h_ratio1.SetStats(0)
    h_ratio1.SetTitle("")
    if makeRatioPlot == 1 and doFracFit:
        denomHist.GetXaxis().SetTitle()
        fPads2.cd()
        fPads2.SetGridy()
        h_ratio1.Divide(fitPlot)

        h_ratio1.GetXaxis().SetTitle(varTitle)
        h_ratio1.GetXaxis().SetNdivisions(515)
        h_ratio1.GetXaxis().SetTitleSize(0.12)
        h_ratio1.GetXaxis().SetLabelSize(0.12)
        h_ratio1.GetXaxis().SetRangeUser(0.0, varXRangeMax)
        h_ratio1.GetYaxis().SetRangeUser(0.0, 2)
        h_ratio1.GetYaxis().SetTitle("ele'/tot. templ.")
        h_ratio1.GetYaxis().SetLabelSize(0.1)
        h_ratio1.GetYaxis().SetTitleSize(0.1)
        h_ratio1.GetYaxis().SetTitleOffset(0.3)
        h_ratio1.Draw("e0p")

        canvas.Update()
        gPad.RedrawAxis()
        gPad.Modified()

    canvas.Print(plotsDir + "/" + varName + "Overlay_" + region + ".pdf")

tfile.Close()

# make plots of real ele frac. in SB
gStyle.SetPadRightMargin(0.04)
mcEleFracGraphs = {}
for etaReg in etaRegions:
    mcEleFracGraphs[etaReg] = TGraphAsymmErrors()
    mcEleFracGraph = mcEleFracGraphs[etaReg]
    mcEleFracGraph.SetName("mcEleFrac_" + etaReg)
    mcEleFracGraph.Divide(
        mcEleInSBHistsByEtaReg[etaReg],
        mcEleTotalHistsByEtaReg[etaReg],
        "cl=0.683 b(1,1) mode",
    )
for etaReg in etaRegions:
    can = TCanvas()
    can.cd()
    can.SetGridy()
    mcEleFracGraphs[etaReg].Draw()
    mcEleFracGraphs[etaReg].GetXaxis().SetTitle("P_{T} [GeV]")
    mcEleFracGraphs[etaReg].GetYaxis().SetTitle("mcEle(SB)/mcEle(tot)")
    mcEleFracGraphs[etaReg].GetYaxis().SetNdivisions(515)
    mcEleFracGraphs[etaReg].SetMaximum(0.05)
    # mcEleFracGraphs[etaReg].SetStats(0)
    mcEleFracGraphs[etaReg].Draw("ap")
    text = TPaveText(0.57, 0.78, 0.67, 0.87, "NDC")
    text.SetBorderSize(0)
    text.SetTextColor(kBlack)
    text.SetTextSize(3e-2)
    text.SetFillStyle(0)
    text.AddText(etaReg)
    text.Draw()
    can.Print(plotsDir + "/mcEleFracInSB" + etaReg + ".pdf")
