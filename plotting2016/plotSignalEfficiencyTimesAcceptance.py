#!/usr/bin/env python

from ROOT import TH1D, TGraphAsymmErrors, TCanvas, gROOT
from plot_class import GetFile, GetHisto

gROOT.SetBatch(False)

# filePath = "/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/nanoV7/2016/analysis/eejj_20oct2020_optFinalSels/condor/analysisClass_lq_eejj___{}/output/analysisClass_lq_eejj___{}_0.root"
filePath = "$LQDATA/nanoV7/2017/analysis/prefire_eejj_23oct2020_optFinalSels/condor/analysisClass_lq_eejj___{}/output/analysisClass_lq_eejj___{}_0.root"
signalNameTemplate = "LQToDEle_M-{}_pair_pythia8"
mass_points = [i for i in range(300, 3100, 100)]  # go from 300-3000 in 100 GeV steps
mass_points.extend([3500, 4000])
# mass_points.remove(2500)  # FIXME 2016
mass_points.remove(3000)  # FIXME 2017
#
#signalNameTemplate = "LQToUE_M-{}_BetaOne_pythia8"
#mass_points = [i for i in range(300, 2100, 100)]  # go from 300-2000 in 100 GeV steps

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
    noCutEntries = hist.GetBinContent(1)
    # print hist.GetXaxis().GetBinLabel(2)
    # finalSelName = "min_M_ej_LQ{}".format(mass)
    # finalSelBin = hist.GetXaxis().FindBin(finalSelName)
    # merged hists have no bin labels; have to use super ugly hack
    finalSelBin = firstFinalSelBin+(i*3)
    # print "mass={}, finalSelBin ={}".format(mass, finalSelBin)
    finalSelEntries = hist.GetBinContent(finalSelBin)
    totalEventsByMass.append(round(noCutEntries, 3))
    eventsAtFinalSelByMass.append(round(finalSelEntries, 3))
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
tcan.Update()

# rep = ''
# while rep not in ['c', 'C']:
#     rep = raw_input('enter "c" to continue: ')
#     if 1 < len(rep):
#         rep = rep[0]
