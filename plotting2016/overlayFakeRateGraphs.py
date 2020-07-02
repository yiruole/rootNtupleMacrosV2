#!/usr/bin/env python

import copy

import ROOT as r

r.gROOT.SetBatch(True)

fileNames = [
        "$LQANA/versionsOfFakeRate/2016/may25/plots.root",
        "$LQANA/versionsOfFakeRate/2017/apr17/plots.root",
        "$LQANA/versionsOfFakeRate/2018/may14/plots.root"
        ]
tfiles = [r.TFile(filename) for filename in fileNames]
legLabels = ['2016', '2017', '2018']
colorList = [r.kSpring-1, r.kAzure+1, r.kBlue, r.kGreen]
markerStyles = [25, 23, 22]

#nBins = 40050
#minRange = -5
#maxRange = 4000
#histNomOrig = r.TH1F("histNomOrig", "histNomOrig", nBins, minRange, maxRange)
#histNanoOrig = r.TH1F("histNanoOrig", "histNanoOrig", nBins, minRange, maxRange)
#histNanoOrig.SetLineColor(2)
#histNanoOrig.SetLineStyle(2)
#
#histNames = [
#    "Ele1_PtHeep",
#    "Ele1_TrkIsoHEEP7",
#    "Ele1_SCEt",
#    "Ele1_Eta",
#    "Ele1_Phi",
#]
#xaxisTitles = [
#    "Ele1_PtHeep [GeV]",
#    "Ele1_TrkIso [GeV]",
#    "Ele1_SCEt [GeV]",
#    "Ele1_Eta",
#    "Ele1_Phi",
#]
#minRangesX = [
#    0,
#    0,
#    0,
#    -4,
#    -2,
#]
#maxRangesX = [
#    2500,
#    7,
#    2000,
#    4,
#    2,
#]
#histNames.extend(
#    ["Ele2_PtHeep", "Ele2_TrkIsoHEEP7", "Ele2_SCEt", "Ele2_Eta", "Ele2_Phi",]
#)
#xaxisTitles.extend(
#    [
#        "Ele2_PtHeep [GeV]",
#        "Ele2_TrkIso [GeV]",
#        "Ele2_SCEt [GeV]",
#        "Ele2_Eta",
#        "Ele2_Phi",
#    ]
#)
#minRangesX.extend(
#    [0, 0, 0, -4, -2,]
#)
#maxRangesX.extend(
#    [2000, 7, 2000, 4, 2,]
#)
#histNames.extend(["sT_eejj", "M_e1e2"])
#xaxisTitles.extend(["S_{T} [GeV]", "M_{ee} [GeV]"])
#minRangesX.extend([0, 0])
#maxRangesX.extend([3000, 2000])
#histNames.extend(
#    ["Jet1_Pt", "Jet1_Eta", "Jet1_Phi",]
#)
#xaxisTitles.extend(
#    ["Jet1_Pt [GeV]", "Jet1_Eta", "Jet1_Phi",]
#)
#minRangesX.extend(
#    [0, -4, -2,]
#)
#maxRangesX.extend(
#    [2000, 4, 2,]
#)
#histNames.extend(
#    ["Jet2_Pt", "Jet2_Eta", "Jet2_Phi",]
#)
#xaxisTitles.extend(
#    ["Jet2_Pt [GeV]", "Jet2_Eta", "Jet2_Phi",]
#)
#minRangesX.extend(
#    [0, -4, -2,]
#)
#maxRangesX.extend(
#    [2000, 4, 2,]
#)

histNames = [
        "frBar_2Jet_template",
        "frEnd1_2Jet_template",
        "frEnd2_2Jet_template"
        ]

xaxisTitle = "Pt [GeV]"

for i, histName in enumerate(histNames):
    print "histName=", histName
    hists = []
    # maxRangeX = maxRangesX[i]
    # minRangeX = minRangesX[i]
    # xaxisTitle = xaxisTitles[i]
    can = r.TCanvas()
    # for ratio
    # pad1 = r.TPad("pad1", "pad1", 0, 0.2, 1, 1)
    # pad1.SetBottomMargin(0.1)
    pad1 = r.TPad("pad1", "pad1", 0, 0, 1, 1)
    # pad1.SetLogy()
    pad1.SetFillColor(0)
    pad1.SetLineColor(0)
    pad1.Draw()
    pad1.cd()
    for j, tfile in enumerate(tfiles):
        # print "file:", tfile.GetName()
        histOrig = tfile.Get(histName)
        histOrig.SetName(histOrig.GetName()+"_"+str(j))
        hist = copy.deepcopy(histOrig)
        hist.SetLineColor(colorList[j])
        hist.SetMarkerColor(colorList[j])
        hist.SetMarkerStyle(markerStyles[j])
        # these are actually graphs
        if j == 0:
            hist.Draw("ap")
        else:
            hist.Draw("p")
        legPos = "top-right"
        if "sT" in histName or "Pt" in histName or "Et" in histName or "e1e2" in histName:
            if "Eta" not in histName:
                rebinFactor = 801  # 534 #445 #267 #178
                hist = hist.Rebin(rebinFactor)
                legPos = "bottom-left"
        hist.GetXaxis().SetTitle(xaxisTitle)
        hist.GetXaxis().SetTitleSize(0.06)
        # hist.GetXaxis().SetRangeUser(minRangeX, maxRangeX)
        hist.SetTitle("")
        if j == 0:
            hist.Draw("ap")
        else:
            hist.Draw("p")
        hists.append(hist)

    # reposition stats box
    pad1.Draw()
    pad1.Modified()
    can.Update()
    pad1.Update()
    #st = histNano.GetListOfFunctions().FindObject("stats")
    #st.SetX1NDC(0.596)
    #st.SetX2NDC(0.796)
    #st.SetY1NDC(0.662)
    #st.SetY2NDC(0.996)
    #st.Draw()
    ##
    #stNom = histNom.GetListOfFunctions().FindObject("stats")
    #stNom.SetX1NDC(0.798)
    #stNom.SetX2NDC(0.997)
    #stNom.SetY1NDC(0.662)
    #stNom.SetY2NDC(0.996)
    #stNom.Draw()
    # legend
    if legPos == "top-left":
        leg = r.TLegend(0.2, 0.7, 0.4, 0.85)
    elif legPos == "bottom-left":
        leg = r.TLegend(0.2, 0.2, 0.4, 0.35)
    elif legPos == "top-right":
        leg = r.TLegend(0.75, 0.7, 0.95, 0.85)
    leg.SetBorderSize(0)
    for j, hist in enumerate(hists):
        leg.AddEntry(hist, legLabels[j], "lp")
    leg.Draw()
    can.Update()
    pad1.Update()

    # # ratio
    # can.cd()
    # pad2 = r.TPad("pad2", "pad2", 0, 0, 1, 0.2)
    # pad2.SetTopMargin(0.075)
    # pad2.SetFillColor(0)
    # pad2.SetLineColor(0)
    # pad2.Draw()
    # #
    # pad2.cd()
    # # histRatio = histNano.Clone()
    # # FIXME convert to BayesDivide
    # # histRatio.Divide(histNom)
    # histRatio = r.TGraphAsymmErrors()
    # histRatio.Divide(histNano, histNom, "pois")
    # histRatio.SetLineStyle(1)
    # histRatio.SetLineColor(1)
    # histRatio.Draw("ap")
    # pad2.Draw()
    # pad2.Update()
    # can.Update()
    # histRatio.GetYaxis().SetRangeUser(0, 2)
    # histRatio.GetXaxis().SetLimits(minRangeX, maxRangeX)
    # histRatio.GetXaxis().SetRangeUser(minRangeX, maxRangeX)
    # histRatio.GetXaxis().SetTitle("")
    # histRatio.SetTitle("")
    # histRatio.GetYaxis().SetTitle("nano/mini")
    # histRatio.GetYaxis().SetTitleSize(0.2)
    # histRatio.GetYaxis().SetTitleOffset(0.2)
    # histRatio.GetXaxis().SetLabelSize(0.1)
    # # histRatio.GetXaxis().SetLabelSize(0.)
    # histRatio.GetYaxis().SetLabelSize(0.15)
    # lineAtOne = r.TLine(
    #     histRatio.GetXaxis().GetXmin(), 1, histRatio.GetXaxis().GetXmax(), 1
    # )
    # lineAtOne.SetLineColor(2)
    # lineAtOne.Draw()

    r.gPad.Modified()
    r.gPad.Update()
    can.Print(histName + "_overlay.png")
    can.Print(histName + "_overlay.pdf")


for tfile in tfiles:
    tfile.Close()

