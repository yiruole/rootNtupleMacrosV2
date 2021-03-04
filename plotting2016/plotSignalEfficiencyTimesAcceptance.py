#!/usr/bin/env python

from ROOT import TH1D, TGraphAsymmErrors, TCanvas, gROOT, TLatex, TFile
from plot_class import GetFile, GetHisto

gROOT.SetBatch(False)

# filePath = "/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/nanoV7/2016/analysis/eejj_20oct2020_optFinalSels/condor/analysisClass_lq_eejj___{}/output/analysisClass_lq_eejj___{}_0.root"
# filePath = "/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/nanoV7/2016/analysis/eejj_26oct2020_optFinalSels_noPdfReweight/condor/analysisClass_lq_eejj___{}/output/analysisClass_lq_eejj___{}_0.root"
# filePath = "$LQDATA/nanoV7/2017/analysis/prefire_eejj_23oct2020_optFinalSels/condor/analysisClass_lq_eejj___{}/output/analysisClass_lq_eejj___{}_0.root"
filePath = "$LQDATA/nanoV7/2016/analysis/eejj_9nov2020_optFinalSelsOld/condor/analysisClass_lq_eejj___{}/output/analysisClass_lq_eejj___{}_0.root"
# LQToDEle
signalNameTemplate = "LQToDEle_M-{}_pair_pythia8"
mass_points = [i for i in range(300, 3100, 100)]  # go from 300-3000 in 100 GeV steps
mass_points.extend([3500, 4000])
if "2016" in filePath:
    mass_points.remove(2500)  # FIXME 2016
elif "2017" in filePath:
    mass_points.remove(3000)  # FIXME 2017
# LQToUE
# signalNameTemplate = "LQToUE_M-{}_BetaOne_pythia8"
# mass_points = [i for i in range(300, 2100, 100)]  # go from 300-2000 in 100 GeV steps
#
if "2016" in filePath:
    year = 2016
elif "2017" in filePath:
    year = 2017
if "UE" in signalNameTemplate:
    signalNameShort = "eeuu"
elif "DEle" in signalNameTemplate:
    signalNameShort = "eedd"

# histNameBase = "histo1D__{}__EventsPassingCuts"
histName = "EventsPassingCuts"
totalEventsByMass = []
eventsAtFinalSelByMass = []
firstFinalSelBin = 39

histTotal = TH1D("total", "total", 38, 250, 4050)
histTotal.Sumw2()
histPass = TH1D("pass", "pass", 38, 250, 4050)
histPass.Sumw2()

for i, mass in enumerate(mass_points):
    sampleName = signalNameTemplate.format(mass)
    filename = filePath.format(sampleName, sampleName)
    tfile = GetFile(filename)
    # histName = histNameBase.format(sampleName)
    hist = GetHisto(histName, tfile)
    # noCutEntries = hist.GetBinContent(1)
    sumOfWeightsHist = GetHisto("SumOfWeights", tfile)
    sumWeights = sumOfWeightsHist.GetBinContent(1)
    lhePdfWeightsHist = tfile.Get("LHEPdfSumw")
    lhePdfWeightSumw = lhePdfWeightsHist.GetBinContent(1)  # sum[genWeight*pdfWeight_0]
    # print hist.GetXaxis().GetBinLabel(2)
    # finalSelName = "min_M_ej_LQ{}".format(mass)
    # finalSelBin = hist.GetXaxis().FindBin(finalSelName)
    # merged hists have no bin labels; have to use super ugly hack
    finalSelBin = firstFinalSelBin+(i*3)
    # print "mass={}, finalSelBin ={}".format(mass, finalSelBin)
    finalSelEntries = hist.GetBinContent(finalSelBin)
    # totalEventsByMass.append(round(noCutEntries, 3))
    totalEventsByMass.append(round(sumWeights, 3))
    print "filename={}".format(filename)
    if "2016" in filename:
        if "LQToBEle" in filename or "LQToDEle" in filename:
            totalEventsByMass.pop()
            totalEventsByMass.append(round(lhePdfWeightSumw, 3))
            print "\tapplying LHEPdfWeight={} to dataset={}".format(lhePdfWeightSumw, filename)+"[instead of original sumWeights={}]".format(sumWeights)

    eventsAtFinalSelByMass.append(round(finalSelEntries, 4))
    histTotal.SetBinContent(histTotal.FindBin(mass), totalEventsByMass[i])
    histTotal.SetBinError(histTotal.FindBin(mass), hist.GetBinError(1))
    histPass.SetBinContent(histPass.FindBin(mass), eventsAtFinalSelByMass[i])
    histPass.SetBinError(histTotal.FindBin(mass), hist.GetBinError(finalSelBin))
    tfile.Close()


print "masses:", mass_points
print "total:", totalEventsByMass
print "pass:", eventsAtFinalSelByMass
print "effAcc:", [n/d for n, d in zip(eventsAtFinalSelByMass, totalEventsByMass)]

# tcan2 = TCanvas()
# tcan2.cd()
# histTotal.Draw()

# tcan3 = TCanvas()
# tcan3.cd()
# histPass.Draw()

tcan = TCanvas()
tcan.cd()
graph = TGraphAsymmErrors()
graph.Divide(histPass, histTotal, "cl=0.683 b(1,1) modev")
#graph.Divide(histPass, histTotal, "cp")
graph.Draw("ap")
graph.GetYaxis().SetRangeUser(0, 0.7)
graph.GetYaxis().SetTitle("Acceptance #times efficiency")
graph.GetYaxis().SetTitleFont(42)
graph.GetYaxis().SetLabelFont(42)
# graph.GetYaxis().SetLabelOffset(0.007)
# graph.GetYaxis().SetLabelSize(0.06)
# graph.GetYaxis().SetTitleOffset(1.)
# graph.GetYaxis().SetTitleSize(0.07)
graph.GetYaxis().CenterTitle(1)
#
graph.GetXaxis().SetTitle('#it{m}_{LQ} [GeV]')
graph.GetXaxis().SetTitleFont(42)
graph.GetXaxis().SetLabelFont(42)
# graph.GetXaxis().SetLabelOffset(0.01)
graph.GetXaxis().SetTitleOffset(1.)
# graph.GetXaxis().SetLabelSize(0.06)
# graph.GetXaxis().SetTitleSize(0.07)
graph.GetXaxis().SetNdivisions(505)
graph.SetName("accTimesEffGraph")
tcan.Update()

# -- draw label
labelOffset = 0.03
ystart = 0.45
xstart = 217
hsize = 0.21
vsize = 0.25
l = TLatex()
l.SetTextAlign(12)
l.SetTextFont(132)
# l.SetTextSize(0.065)
l.SetTextSize(0.055)
l.DrawLatex(
    xstart - hsize + 0, ystart + vsize - 0.05, "CMS Preliminary {}".format(year)
)
# l.DrawLatex(
#     xstart - hsize + 0, ystart + vsize - 0.10, "#it{Simulation}"
# )
l.DrawLatex(
    xstart - hsize + 0, ystart + vsize - 0.10, "Scalar LQ #bar{LQ} #rightarrow "+signalNameShort
)
tcan.Update()
tcan.Print("accTimesEff.pdf")

outFile = TFile("accTimesEff.root", "recreate")
tcan.Write()
graph.Write()
outFile.Close()
