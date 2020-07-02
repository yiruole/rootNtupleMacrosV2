#!/usr/bin/env python

##############################################################################
# # USER CODE IS TOWARD THE END OF THE FILE
##############################################################################

##############################################################################
# ############ DON'T NEED TO MODIFY ANYTHING HERE - BEGIN #####################

# ---Import
import sys
import string

# from optparse import OptionParser
import os.path
from ROOT import kTRUE, gROOT, gStyle, TFile, TCanvas, TRandom3, kWhite
import re

# from array import array
import copy
import math

# --- ROOT general options
gROOT.SetBatch(kTRUE)
gStyle.SetOptStat(0)
gStyle.SetPalette(1)
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
# --- TODO: WHY IT DOES NOT LOAD THE DEFAULT ROOTLOGON.C ? ---#


def GetFile(filename):
    file = TFile(filename)
    if not file:
        print "ERROR: file " + filename + " not found"
        print "exiting..."
        sys.exit()
    return file


def GetHisto(histoName, file, scale=1):
    histo = file.Get(histoName)
    if not histo:
        print "ERROR: histo " + histoName + " not found in " + file.GetName()
        print "exiting..."
        sys.exit()
    new = copy.deepcopy(histo)
    if scale != 1:
        new.Scale(scale)
    return new


def GetIntegralTH1(histo, xmin, xmax, verbose=False):
    # get integral
    axis = histo.GetXaxis()
    bmin = axis.FindBin(xmin)
    bmax = axis.FindBin(xmax)
    bminResidual = (
        histo.GetBinContent(bmin)
        * (xmin - axis.GetBinLowEdge(bmin))
        / axis.GetBinWidth(bmin)
    )
    bmaxResidual = (
        histo.GetBinContent(bmax)
        * (axis.GetBinUpEdge(bmax) - xmax)
        / axis.GetBinWidth(bmax)
    )
    integral = histo.Integral(bmin, bmax) - bminResidual - bmaxResidual
    if verbose:
        print "GetIntegralTH1(" + histo.GetName(), xmin, xmax, ")=", integral
    return integral


def GetErrorIntegralTH1(histo, xmin, xmax, verbose=False):
    if verbose:
        print "## calculating error for integral of histo " + str(histo)
        print "## in the x range [" + str(xmin) + "," + str(xmax) + "]"
    # get integral
    axis = histo.GetXaxis()
    bmin = axis.FindBin(xmin)
    bmax = axis.FindBin(xmax)
    bminResidual = (
        histo.GetBinContent(bmin)
        * (xmin - axis.GetBinLowEdge(bmin))
        / axis.GetBinWidth(bmin)
    )
    bmaxResidual = (
        histo.GetBinContent(bmax)
        * (axis.GetBinUpEdge(bmax) - xmax)
        / axis.GetBinWidth(bmax)
    )
    integral = histo.Integral(bmin, bmax) - bminResidual - bmaxResidual
    error = 0
    for bin in range(bmin, bmax + 1):
        # print "bin: " +str(bin)
        if bin == bmax and bmaxResidual == histo.GetBinContent(
            bmax
        ):  # skip last bin if out of range
            if verbose:
                print "     --> skip bin: " + str(bin)
        else:
            error = error + histo.GetBinError(bin) ** 2
            # print "error**2 : " + str(error)

    error = math.sqrt(error)
    if verbose:
        print " "
    return error


# The Plot class: add members if needed
class Plot:
    histoDATA = ""  # DATA
    histoTTbar = ""  # MCTTbar
    histoMCall = ""  # MCall
    histoQCD = ""  # QCD
    histoZJet = ""
    histoWJet = ""
    histoSingleTop = ""
    histoPhotonJets = ""
    histoDiboson = ""
    xtit = ""  # xtitle
    ytit = ""  # ytitle
    xmin = ""  # set xmin to calculate rescaling factor (-- please take into account the bin size of histograms --)
    xmax = ""  # # set xmax to calculate rescaling factor (-- please take into account the bin size of histograms --)
    xminplot = ""  # min x axis range (need to set both min and max. Leave it as is for full range)
    xmaxplot = ""  # max x axis range (need to set both min and max. Leave it as is for full range)
    yminplot = ""  # min y axis range (need to set both min and max. Leave it as is for full range)
    ymaxplot = ""  # max y axis range (need to set both min and max. Leave it as is for full range)
    lpos = (
        ""  # legend position (default = top-right, option="bottom-center", "top-left")
    )
    #    xlog         = "" # log scale of X axis (default = no, option="yes") ### IT SEEMS IT DOES NOT WORK
    ylog = ""  # log scale of Y axis (default = no, option="yes")
    # rebin       = "" # rebin x axis (default = 1, option = set it to whatever you want )
    name = ""  # name of the final plots
    lint = "2.6 fb^{-1}"  # integrated luminosity of the sample ( example "10 pb^{-1}" )
    fileXsectionNoRescale = ""  # cross section file (with no rescale
    datasetName = ""  # string for pattern recognition of dataset name (rescaling will be done only on matched datasets)

    def CheckMCDataConsistency(self):
        # checks
        if self.histoMCall.GetNbinsX() != self.histoDATA.GetNbinsX():
            print "WARNING! number of bins is different between DATA and MC"
            print "exiting..."
            sys.exit()
        if self.histoMCall.GetBinWidth(1) != self.histoDATA.GetBinWidth(1):
            print "WARNING! bin width is different between DATA and MC"
            print "exiting..."
            sys.exit()


def GetRTTBarDYJets(plotObjTTBar, plotObjDYJets, randomize=False, verbose=False):

    # integrals: ttbar
    integralDATA_ttbar = GetIntegralTH1(
        plotObjTTBar.histoDATA, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralDATA_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoDATA, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    integralMCall_ttbar = GetIntegralTH1(
        plotObjTTBar.histoMCall, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralMCall_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoMCall, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    integralTTbar_ttbar = GetIntegralTH1(
        plotObjTTBar.histoTTbar, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralTTbar_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoTTbar, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    integralDYJets_ttbar = GetIntegralTH1(
        plotObjTTBar.histoZJet, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralDYJets_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoZJet, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    integralQCD_ttbar = GetIntegralTH1(
        plotObjTTBar.histoQCD, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralQCD_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoQCD, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    # contamination from other backgrounds (except TTbar and DYJets) in the integral range [QCD is not in MCall]
    integralMCothers_ttbar = (
        integralMCall_ttbar - integralTTbar_ttbar - integralDYJets_ttbar
    )
    ERRintegralMCothers_ttbar = math.sqrt(
        ERRintegralMCall_ttbar ** 2 + ERRintegralTTbar_ttbar ** 2
    )

    # integrals: wjets
    integralDATA_dyjets = GetIntegralTH1(
        plotObjDYJets.histoDATA, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralDATA_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoDATA, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    integralMCall_dyjets = GetIntegralTH1(
        plotObjDYJets.histoMCall, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralMCall_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoMCall, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    integralTTbar_dyjets = GetIntegralTH1(
        plotObjDYJets.histoTTbar, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralTTbar_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoTTbar, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    integralDYJets_dyjets = GetIntegralTH1(
        plotObjDYJets.histoZJet, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralDYJets_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoZJet, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    integralQCD_dyjets = GetIntegralTH1(
        plotObjDYJets.histoQCD, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralQCD_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoQCD, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    # contamination from other backgrounds (except dyjets and TTBar) in the integral range [QCD is not in MCall]
    integralMCothers_dyjets = (
        integralMCall_dyjets - integralDYJets_dyjets - integralTTbar_dyjets
    )
    ERRintegralMCothers_dyjets = math.sqrt(
        ERRintegralMCall_dyjets ** 2 + ERRintegralDYJets_dyjets ** 2
    )
    contamination_dyjets = (
        integralMCothers_dyjets + integralTTbar_dyjets + integralQCD_dyjets
    ) / (integralMCall_dyjets + integralQCD_dyjets)

    # randomize
    trand = TRandom3(0)
    integralDATA_ttbar = trand.Gaus(integralDATA_ttbar, ERRintegralDATA_ttbar)
    integralMCall_ttbar = trand.Gaus(integralMCall_ttbar, ERRintegralMCall_ttbar)
    integralTTbar_ttbar = trand.Gaus(integralTTbar_ttbar, ERRintegralTTbar_ttbar)
    integralDYJets_ttbar = trand.Gaus(integralDYJets_ttbar, ERRintegralDYJets_ttbar)
    integralQCD_ttbar = trand.Gaus(integralQCD_ttbar, ERRintegralQCD_ttbar)
    integralMCothers_ttbar = (
        integralMCall_ttbar - integralTTbar_ttbar - integralDYJets_ttbar
    )
    #
    integralDATA_dyjets = trand.Gaus(integralDATA_dyjets, ERRintegralDATA_dyjets)
    integralMCall_dyjets = trand.Gaus(integralMCall_dyjets, ERRintegralMCall_dyjets)
    integralTTbar_dyjets = trand.Gaus(integralTTbar_dyjets, ERRintegralTTbar_dyjets)
    integraldyjets_dyjets = trand.Gaus(integralDYJets_dyjets, ERRintegralDYJets_dyjets)
    integralQCD_dyjets = trand.Gaus(integralQCD_dyjets, ERRintegralQCD_dyjets)

    # solve the system of equations
    # (1) --> dyjets
    # (2) --> ttbar
    rTTBar = (
        integralDATA_dyjets * integralDYJets_ttbar
        - (integralMCothers_dyjets + integralQCD_dyjets) * integralDYJets_ttbar
    )
    rTTBar += (
        integralMCothers_ttbar + integralQCD_ttbar - integralDATA_ttbar
    ) * integraldyjets_dyjets
    if verbose:
        print integralTTbar_dyjets, integralDYJets_ttbar, integralTTbar_ttbar, integralDYJets_dyjets
    try:
        rTTBar /= (
            integralTTbar_dyjets * integralDYJets_ttbar
            - integralTTbar_ttbar * integralDYJets_dyjets
        )
    except ZeroDivisionError:
        print "ERROR: ZeroDivisionError: one or more of the integrals above are zero"
        return -1, -1
    rDYJets = (
        integralDATA_dyjets * integralTTbar_ttbar
        - (integralMCothers_dyjets + integralQCD_dyjets) * integralTTbar_ttbar
    )
    rDYJets += (
        integralMCothers_ttbar + integralQCD_ttbar - integralDATA_ttbar
    ) * integralTTbar_dyjets
    rDYJets /= (
        integralTTbar_ttbar * integralDYJets_dyjets
        - integralTTbar_dyjets * integralDYJets_ttbar
    )

    return rTTBar, rDYJets


def CalculateRescaleFactor(plotObjTTBar, plotObjDYJets, fileps):
    # calculate rescaling factor for Z/gamma+jet background and create new cross section file
    canvas = TCanvas()

    plotObjTTBar.CheckMCDataConsistency()
    plotObjDYJets.CheckMCDataConsistency()

    # do the following N times, so we can calculate the stat. uncertainty
    N = 100
    # NB: the commented code below makes a nice progress bar but causes the dict to be undefined...
    # steps = N
    print "Randomizing histos and calculating scale factors:",
    # progressString = "0% [" + " " * steps + "] 100%"
    # print progressString,
    # print "\b" * (len(progressString) - 3),
    sys.stdout.flush()
    rTTBarList = []
    rDYJetsList = []
    for i in range(0, N):
        # print ""
        rTTBar, rDYJets = GetRTTBarDYJets(plotObjTTBar, plotObjDYJets, True)
        rTTBarList.append(rTTBar)
        rDYJetsList.append(rDYJets)
        # print "\b.",
        # sys.stdout.flush()
    # print "\b] 100%"
    print "Done."

    ttMean = 0
    for rtt in rTTBarList:
        ttMean += rtt
    ttMean /= N
    #
    wMean = 0
    for rw in rDYJetsList:
        wMean += rw
    wMean /= N

    rTTBarSigma = 0
    for rtt in rTTBarList:
        rTTBarSigma += pow(rtt - ttMean, 2)
    rTTBarSigma /= N
    rTTBarSigma = math.sqrt(rTTBarSigma)
    #
    rDYJetsSigma = 0
    for rdy in rDYJetsList:
        rDYJetsSigma += pow(rdy - wMean, 2)
    rDYJetsSigma /= N
    rDYJetsSigma = math.sqrt(rDYJetsSigma)

    # integrals: ttbar
    integralDATA_ttbar = GetIntegralTH1(
        plotObjTTBar.histoDATA, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralDATA_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoDATA, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    integralMCall_ttbar = GetIntegralTH1(
        plotObjTTBar.histoMCall, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralMCall_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoMCall, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    integralTTbar_ttbar = GetIntegralTH1(
        plotObjTTBar.histoTTbar, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    print "plotObjTTBar.histoTTbar=", plotObjTTBar.histoTTbar
    ERRintegralTTbar_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoTTbar, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    integralDYJets_ttbar = GetIntegralTH1(
        plotObjTTBar.histoZJet, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralDYJets_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoZJet, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    integralQCD_ttbar = GetIntegralTH1(
        plotObjTTBar.histoQCD, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralQCD_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoQCD, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    # contamination from other backgrounds (except TTbar and DYJets) in the integral range [QCD is not in MCall]
    integralMCothers_ttbar = (
        integralMCall_ttbar - integralTTbar_ttbar - integralDYJets_ttbar
    )
    ERRintegralMCothers_ttbar = math.sqrt(
        ERRintegralMCall_ttbar ** 2 + ERRintegralTTbar_ttbar ** 2
    )
    try:
        contamination_ttbar = (
            integralMCothers_ttbar + integralDYJets_ttbar + integralQCD_ttbar
        ) / (integralMCall_ttbar + integralQCD_ttbar)
    except ZeroDivisionError:
        print "ERROR: ZeroDivisionError: integralMCall_ttbar+integralQCD_ttbar is zero; integralMCall_ttbar=", integralMCall_ttbar, "integralQCD_ttbar=", integralQCD_ttbar
        contamination_ttbar = -1

    # integrals: dyjets
    integralDATA_dyjets = GetIntegralTH1(
        plotObjDYJets.histoDATA, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralDATA_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoDATA, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    integralMCall_dyjets = GetIntegralTH1(
        plotObjDYJets.histoMCall, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralMCall_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoMCall, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    integralTTbar_dyjets = GetIntegralTH1(
        plotObjDYJets.histoTTbar, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    print "plotObjDYJets.histoTTbar=", plotObjDYJets.histoTTbar
    ERRintegralTTbar_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoTTbar, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    integralDYJets_dyjets = GetIntegralTH1(
        plotObjDYJets.histoZJet, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralDYJets_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoZJet, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    integralQCD_dyjets = GetIntegralTH1(
        plotObjDYJets.histoQCD, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralQCD_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoQCD, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    # contamination from other backgrounds (except WJets and TTBar) in the integral range [QCD is not in MCall]
    integralMCothers_dyjets = (
        integralMCall_dyjets - integralDYJets_dyjets - integralTTbar_dyjets
    )
    ERRintegralMCothers_dyjets = math.sqrt(
        ERRintegralMCall_dyjets ** 2 + ERRintegralDYJets_dyjets ** 2
    )
    contamination_dyjets = (
        integralMCothers_dyjets + integralTTbar_dyjets + integralQCD_dyjets
    ) / (integralMCall_dyjets + integralQCD_dyjets)

    # solve the system of equations
    # (1) --> dyjets
    # (2) --> ttbar
    rTTBar = (
        integralDATA_dyjets * integralDYJets_ttbar
        - (integralMCothers_dyjets + integralQCD_dyjets) * integralDYJets_ttbar
    )
    rTTBar += (
        integralMCothers_ttbar + integralQCD_ttbar - integralDATA_ttbar
    ) * integralDYJets_dyjets
    try:
        rTTBar /= (
            integralTTbar_dyjets * integralDYJets_ttbar
            - integralTTbar_ttbar * integralDYJets_dyjets
        )
    except ZeroDivisionError:
        print "ERROR in rTTBar: ZeroDivisionError: one or more of the integrals is zero:"
        print integralTTbar_dyjets, integralDYJets_ttbar, "-", integralTTbar_ttbar, integralDYJets_dyjets
    rDYJets = (
        integralDATA_dyjets * integralTTbar_ttbar
        - (integralMCothers_dyjets + integralQCD_dyjets) * integralTTbar_ttbar
    )
    rDYJets += (
        integralMCothers_ttbar + integralQCD_ttbar - integralDATA_ttbar
    ) * integralTTbar_dyjets
    try:
        rDYJets /= (
            integralTTbar_ttbar * integralDYJets_dyjets
            - integralTTbar_dyjets * integralDYJets_ttbar
        )
    except ZeroDivisionError:
        print "ERROR in rWJets: ZeroDivisionError: one or more of the integrals is zero:"
        print integralTTbar_ttbar, integralDYJets_dyjets, "-", integralTTbar_dyjets, integralDYJets_ttbar

    # FIXME
    ##draw histo
    # self.histoMCall.SetFillColor(kBlue)
    # self.histoDATA.SetMarkerStyle(20)

    # self.histoMCall.Draw("HIST")
    # self.histoDATA.Draw("psame")
    # self.histoMCall.GetXaxis().SetRangeUser(self.xminplot,self.xmaxplot)
    # self.histoMCall.GetYaxis().SetRangeUser(self.yminplot,self.ymaxplot)

    # canvas.Update()
    # gPad.RedrawAxis()
    # gPad.Modified()
    ##canvas.SaveAs(self.name + ".eps","eps")
    ##canvas.SaveAs(self.name + ".pdf","pdf")
    # canvas.Print(fileps)
    # canvas.Print(self.name + ".C")
    ## make root file
    # tfile = TFile(self.name+'.root','recreate')
    # tfile.cd()
    # self.histoDATA.Write()
    # self.histoTTbar.Write()
    # self.histoMCall.Write()
    # self.histoQCD.Write()
    # self.histoZJet.Write()
    # self.histoWJet.Write()
    # self.histoSingleTop.Write()
    # self.histoPhotonJets.Write()
    # self.histoDiboson.Write()
    # tfile.Close()

    # printout
    print
    print " TTBar "
    print "######################################## "
    print "name:                         " + plotObjTTBar.name
    print "integral range:               " + str(
        plotObjTTBar.xmin
    ) + " < Mee < " + str(plotObjTTBar.xmax) + " GeV/c2"
    print "integral MC All:              " + str(integralMCall_ttbar) + " +/- " + str(
        ERRintegralMCall_ttbar
    )
    print "integral QCD:                 " + str(integralQCD_ttbar) + " +/- " + str(
        ERRintegralQCD_ttbar
    )
    print "integral MC TTbar:            " + str(integralTTbar_ttbar) + " +/- " + str(
        ERRintegralTTbar_ttbar
    )
    print "integral MC DYJets:            " + str(integralDYJets_ttbar) + " +/- " + str(
        ERRintegralDYJets_ttbar
    )
    print "integral MC other:            " + str(
        integralMCothers_ttbar
    ) + " +/- " + str(ERRintegralMCothers_ttbar)
    print "rescaled integral MC TTbar:   " + str(
        rTTBar * integralTTbar_ttbar
    ) + " +/- " + str(rTTBar * ERRintegralTTbar_ttbar)
    print "rescaled integral MC DYJets:   " + str(
        rDYJets * integralDYJets_ttbar
    ) + " +/- " + str(rDYJets * ERRintegralDYJets_ttbar)
    print "integral DATA:                " + str(integralDATA_ttbar) + " +/- " + str(
        ERRintegralDATA_ttbar
    )
    print "contribution from other bkgs (except TTbar): " + str(
        contamination_ttbar * 100
    ) + "%"
    # print "integral DATA (corrected for contribution from other bkgs): "  + str( integralDATAcorr_ttbar ) + " +/- " + str( ERRintegralDATAcorr_ttbar )
    print "rescale factor for TTbar background: " + str(rTTBar) + " +/- " + str(
        rTTBarSigma
    )
    # print "systematical uncertainty of TTbar background modeling: " + str(relERRrescale*100) + "%"
    print "######################################## "
    print

    print
    print " DYJets "
    print "######################################## "
    print "name:                         " + plotObjDYJets.name
    print "integral range:               " + str(
        plotObjDYJets.xmin
    ) + " < Mee < " + str(plotObjDYJets.xmax) + " GeV/c2"
    print "integral MC All:              " + str(integralMCall_dyjets) + " +/- " + str(
        ERRintegralMCall_dyjets
    )
    print "integral QCD:                 " + str(integralQCD_dyjets) + " +/- " + str(
        ERRintegralQCD_dyjets
    )
    print "integral MC TTbar:            " + str(integralTTbar_dyjets) + " +/- " + str(
        ERRintegralTTbar_dyjets
    )
    print "integral MC DYJets:            " + str(integralDYJets_dyjets) + " +/- " + str(
        ERRintegralDYJets_dyjets
    )
    print "integral MC other:            " + str(
        integralMCothers_dyjets
    ) + " +/- " + str(ERRintegralMCothers_dyjets)
    print "rescaled integral MC TTbar:   " + str(
        rTTBar * integralTTbar_dyjets
    ) + " +/- " + str(rTTBar * ERRintegralTTbar_dyjets)
    print "rescaled integral MC DYJets:   " + str(
        rDYJets * integralDYJets_dyjets
    ) + " +/- " + str(rDYJets * ERRintegralDYJets_dyjets)
    print "integral DATA:                " + str(integralDATA_dyjets) + " +/- " + str(
        ERRintegralDATA_dyjets
    )
    print "contribution from other bkgs (except wjets): " + str(
        contamination_dyjets * 100
    ) + "%"
    # print "integral DATA (corrected for contribution from other bkgs): "  + str( integralDATAcorr_wjets ) + " +/- " + str( ERRintegralDATAcorr_wjets )
    print "rescale factor for DYJets background: " + str(rDYJets) + " +/- " + str(
        rDYJetsSigma
    )
    print "######################################## "
    print

    # create new cross section file -- ttbar
    originalFileName = string.split(
        string.split(plotTTbar.fileXsectionNoRescale, "/")[-1], "."
    )[0]
    ttbarFileName = originalFileName + "_" + plotObjTTBar.name + ".txt"
    os.system("rm -f " + ttbarFileName)
    outputFile = open(ttbarFileName, "w")

    for line in open(plotTTbar.fileXsectionNoRescale):
        line = string.strip(line, "\n")
        lineNoComments = line.split("#")[
            0
        ]  # strip off anything after any '#' if present
        # ignore empty lines
        if len(lineNoComments) <= 0:
            print >> outputFile, line
            continue
        # FIXME this doesn't support keeping comments at the end of lines

        if re.search(plotTTbar.datasetName, line):
            list = re.split("\s+", line)
            newline = (
                str(list[0]) + "    " + str("%.6f" % (float(list[1]) * float(rTTBar)))
            )
            print >> outputFile, newline
        else:
            print >> outputFile, line

    outputFile.close()
    print "New xsection file (after TTbar rescaling) is: " + ttbarFileName
    print " "

    # create new cross section file -- DYJets
    originalFileName = string.split(string.split(ttbarFileName, "/")[-1], ".")[0]
    newFileName = originalFileName + "_" + plotObjDYJets.name + ".txt"
    os.system("rm -f " + newFileName)
    outputFile = open(newFileName, "w")

    for line in open(ttbarFileName):
        line = string.strip(line, "\n")
        lineNoComments = line.split("#")[
            0
        ]  # strip off anything after any '#' if present
        # ignore empty lines
        if len(lineNoComments) <= 0:
            print >> outputFile, line
            continue

        if re.search(plotDYJets.datasetName, line):
            list = re.split("\s+", line)
            newline = (
                str(list[0]) + "    " + str("%.6f" % (float(list[1]) * float(rDYJets)))
            )
            print >> outputFile, newline
        else:
            print >> outputFile, line

    outputFile.close()
    print "New xsection file (after WJets rescaling) is: " + newFileName
    print " "


# ############ DON'T NEED TO MODIFY ANYTHING HERE - END #######################
##############################################################################


##############################################################################
# ############ USER CODE - BEGIN ##############################################

# --- Input files
File_QCD_preselection = GetFile(
        # "$LQDATA/nanoV6/2016/analysis/qcdYield_24jun2020/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root"
        # "$LQDATA/nanoV6/2017/analysis/qcdYield_25jun2020/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root"
        "$LQDATA/nanoV6/2018/analysis/qcdYield_25jun2020/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root"
)

File_preselection = GetFile(
    # "$LQDATA/nanoV6/2016/analysis/prefire_19may2020/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root"
    # "$LQDATA/nanoV6/2017/analysis/noPrefire_22may2020/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root"
    # "$LQDATA/nanoV6/2017/analysis/prefire_22may2020/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root"
    "$LQDATA/nanoV6/2018/analysis/eejj_5may2020/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root"
)

xsectionFile = "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleAnalyzerV2/config/xsection_13TeV_2015.txt"

histNameDefault = "histo1D__SAMPLE__"
# histNameReleaseMee = "histo1D__SAMPLE__cutHisto_allOtherCuts___________"
histNameAllPrevCuts = "histo1D__SAMPLE__cutHisto_allPreviousCuts________"
#
# allBkg = "ALLBKG_powhegTTBar_ZJetAMCJetPtBinnedWJetAMCJetBinned_DibosonPyth"
# allBkg = "ALLBKG_powhegTTBar_ZJetPtWJetInc_NLODiboson_triboson"
if "2016" in File_preselection.GetName():
    # wjet = "WJet_amcatnlo_Inc"
    wjet = "WJet_amcatnlo_jetBinned"
    zjet = "ZJet_amcatnlo_ptBinned"
    zjetDatasetName = "DYJetsToLL_Pt.+Tune"
    diboson = "DIBOSON_nlo"
    allBkg = "ALLBKG_powhegTTBar_ZJetPtWJetAMCJetBinned_NLODiboson_triboson"
else:
    wjet = "WJet_amcatnlo_jetBinned"
    zjet = "ZJet_jetAndPtBinned"
    zjetDatasetName = "DY.+JetsToLL_M-50_LHEZpT_.+Tune"
    diboson = "DIBOSON_nlo"
    allBkg = "ALLBKG_powhegTTBar_ZJetAMCJetPtBinnedWJetAMCJetBinned_NLODiboson_triboson"
# wjet = "WJet_amcatnlo_ptBinned"
# diboson = "DIBOSON"
ttbar = "TTbar_powheg"
ttbarDatasetName = "TT"
singletop = "SingleTop"
photonjets = "PhotonJets_Madgraph"
# qcd = "QCD_EMEnriched"
qcd = "QCDFakes_DATA"
# File_QCD_preselection = File_preselection

print "INFO: using file: " + File_preselection.GetName()
print "INFO: using QCD file: " + File_QCD_preselection.GetName()
print "INFO: using samples:"
print "\t allBkg ------>", allBkg
print "\t DY ---------->", zjet, "; datasetname =", zjetDatasetName
print "\t W ----------->", wjet
print "\t ttbar ------->", ttbar, "; datasetname =", ttbarDatasetName
print "\t diboson ----->", diboson
print "\t QCD --------->", qcd
print "\t SingleTop --->", singletop
print "\t PhotonJets -->", photonjets
# --- Rescaling of W + jet and ttbar+jets background

histBaseNames = []
# nominal
histBaseNames.append("Mee_PAS")
histBaseNames.append("Mee_80_100_Preselection")
histBaseNames.append("Mee_70_110_Preselection")
histBaseNames.append("Mee_EBEB_PAS")
histBaseNames.append("Mee_EEEE_PAS")
##histBaseNames.append('Mee_NJetEq2_PAS')
##histBaseNames.append('Mee_NJetEq3_PAS')
##histBaseNames.append('Mee_NJetEq4_PAS')
##histBaseNames.append('Mee_NJetEq5_PAS')
##histBaseNames.append('Mee_NJetEq6_PAS')
##histBaseNames.append('Mee_NJetEq7_PAS')
##histBaseNames.append('Mee_NJetGeq3_PAS')
##histBaseNames.append('Mee_NJetGeq4_PAS')
histBaseNames.append("Mee_sT300To500_PAS")
histBaseNames.append("Mee_sT500To750_PAS")
histBaseNames.append("Mee_sT750To1250_PAS")
histBaseNames.append("Mee_sT1250ToInf_PAS")
histBaseNames.append("Mee_MejMin100To200_PAS")
histBaseNames.append("Mee_MejMin200To300_PAS")
histBaseNames.append("Mee_MejMin300To400_PAS")
histBaseNames.append("Mee_MejMin400To500_PAS")
histBaseNames.append("Mee_MejMin500To650_PAS")
histBaseNames.append("Mee_MejMin650ToInf_PAS")
##histBaseNames.append( "Mee_sT340_PAS")
##histBaseNames.append( "Mee_sT405_PAS")
##histBaseNames.append( "Mee_sT470_PAS")
##histBaseNames.append( "Mee_sT535_PAS")
##histBaseNames.append( "Mee_sT595_PAS")
##histBaseNames.append( "Mee_sT660_PAS")
##histBaseNames.append( "Mee_sT720_PAS")
##histBaseNames.append( "Mee_sT780_PAS")
##histBaseNames.append( "Mee_sT840_PAS")
##histBaseNames.append( "Mee_sT900_PAS")
##histBaseNames.append( "Mee_sT960_PAS")
##histBaseNames.append( "Mee_sT1015_PAS")
##histBaseNames.append( "Mee_sT1075_PAS")
##histBaseNames.append( "Mee_sT1130_PAS")
##histBaseNames.append( "Mee_sT1190_PAS")
##histBaseNames.append( "Mee_sT1245_PAS")
##histBaseNames.append( "Mee_sT1300_PAS")
##histBaseNames.append( "Mee_sT1355_PAS")
##histBaseNames.append( "Mee_sT1410_PAS")
##histBaseNames.append( "Mee_sT1460_PAS")
##histBaseNames.append( "Mee_sT1515_PAS")
##histBaseNames.append( "Mee_sT1565_PAS")
##histBaseNames.append( "Mee_sT1615_PAS")
##histBaseNames.append( "Mee_sT1670_PAS")
##histBaseNames.append( "Mee_sT1720_PAS")
##histBaseNames.append( "Mee_sT1770_PAS")
##histBaseNames.append( "Mee_sT1815_PAS")
histBaseNames.append("Mee_70_110_LQ300")
histBaseNames.append("Mee_70_110_LQ600")
histBaseNames.append("Mee_70_110_LQ800")
histBaseNames.append("Mee_70_110_LQ900")

meeMinTTBar = 100
meeMaxTTBar = 200 # TODO check range
meeMinDYJets = 80
meeMaxDYJets = 100

plotsTTBar = []
plotsDYJets = []

for idx, histBaseName in enumerate(histBaseNames):
    # -----------------------------------------
    # for ttbar-enriched region
    # plotBaseName = 'MTenu_50_110_Njet_gte4'
    # plotBaseName = 'MTenu_50_110_gteOneBtaggedJet'
    # plotBaseName = 'MTenu_70_150_gteOneBtaggedJet'
    # plotBaseName = 'MTenu_110_190_gteOneBtaggedJet'
    # plotBaseName= 'MTenu_50_110_gteOneBtaggedJet_LQ800'
    thisHistName = histNameDefault

    # plotBaseName = histBaseName.replace("BJETBIN", "gteOneBtaggedJet")
    plotBaseName = histBaseName
    print "for TTBar, using plotBaseName:", plotBaseName

    h_ALLBKG_powheg_ttbar = GetHisto(
        thisHistName.replace("SAMPLE", allBkg)
        + plotBaseName,
        File_preselection,
    )  # MC all

    # powheg ttbar
    h_TTbar_powheg_ttbar = GetHisto(
        thisHistName.replace("SAMPLE", ttbar) + plotBaseName, File_preselection
    )  # MC TTbar
    h_ZJets_amcatnlo_ttbar = GetHisto(
        thisHistName.replace("SAMPLE", zjet) + plotBaseName,
        File_preselection,
    )
    h_WJets_amcatnlo_ttbar = GetHisto(
        thisHistName.replace("SAMPLE", wjet) + plotBaseName,
        File_preselection,
    )
    # h_WJets_amcatnlo_ttbar = GetHisto("histo1D__WJet_amcatnlo_Inc__"+plotBaseName, File_preselection)
    h_SingleTop_ttbar = GetHisto(
        thisHistName.replace("SAMPLE", singletop) + plotBaseName, File_preselection
    )
    h_PhotonJets_ttbar = GetHisto(
        thisHistName.replace("SAMPLE", photonjets) + plotBaseName,
        File_preselection,
    )
    # h_Diboson_ttbar = GetHisto("histo1D__DIBOSON__"+plotBaseName, File_preselection)
    h_Diboson_ttbar = GetHisto(
        thisHistName.replace("SAMPLE", diboson) + plotBaseName,
        File_preselection,
    )

    # DATA
    h_DATA_ttbar = GetHisto(
        thisHistName.replace("SAMPLE", "DATA") + plotBaseName, File_preselection
    )  # DATA
    # QCD
    h_QCD = GetHisto(
        thisHistName.replace("SAMPLE", qcd)
        #+ plotBaseName.replace("_btagSFDownShift", "").replace("_btagSFUpShift", ""),
        + plotBaseName,
        File_QCD_preselection,
    )
    # h_QCD_ttbar = GetHisto("histo1D__QCD_EMEnriched__"+plotBaseName,File_QCD_preselection)

    plotTTbar = Plot()
    plotTTbar.histoDATA = h_DATA_ttbar
    plotTTbar.histoMCall = h_ALLBKG_powheg_ttbar
    plotTTbar.histoTTbar = h_TTbar_powheg_ttbar
    # plotTTbar.histoQCD = h_QCD_ttbar
    plotTTbar.histoQCD = h_QCD
    plotTTbar.histoZJet = h_ZJets_amcatnlo_ttbar
    plotTTbar.histoWJet = h_WJets_amcatnlo_ttbar

    plotTTbar.histoSingleTop = h_SingleTop_ttbar
    plotTTbar.histoPhotonJets = h_PhotonJets_ttbar
    plotTTbar.histoDiboson = h_Diboson_ttbar
    plotTTbar.xmin = meeMinTTBar
    # plotTTbar.xmax = h_TTbar_amcatnlo_ttbar.GetXaxis().GetXmax()
    plotTTbar.xmax = meeMaxTTBar
    # plotTTbar.name = "TTbarRescale"
    plotTTbar.name = plotBaseName + "_TTbar"
    plotTTbar.fileXsectionNoRescale = xsectionFile
    plotTTbar.xminplot = 0
    plotTTbar.xmaxplot = 2000
    plotTTbar.yminplot = 0
    plotTTbar.ymaxplot = 2000
    plotTTbar.datasetName = ttbarDatasetName
    # plotTTbar.datasetName = "TT"
    # plot0.datasetName = "DYJetsToLL_M-50_HT.+Tune"
    # plot0.datasetName = "Z.+Jets_Pt.+alpgen"
    # example: this match with /Z3Jets_Pt300to800-alpgen/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
    plotsTTBar.append(plotTTbar)

    # -----------------------------------------
    # for wjets-enriched region
    # plotBaseName = 'MTenu_50_110_Njet_lte3'
    # plotBaseName = 'MTenu_50_110_noBtaggedJets'
    # plotBaseName = 'MTenu_70_150_noBtaggedJets'
    # plotBaseName = 'MTenu_110_190_noBtaggedJets'
    # plotBaseName= 'MTenu_50_110_noBtaggedJets_LQ800'

    # plotBaseName = histBaseName.replace("BJETBIN", "noBtaggedJets")
    print "for DYJets, using plotBaseName:", plotBaseName

    h_ALLBKG_powheg_dyjets = GetHisto(
        thisHistName.replace("SAMPLE", allBkg)
        + plotBaseName,
        File_preselection,
    )  # MC all
    h_TTbar_powheg_dyjets = GetHisto(
        thisHistName.replace("SAMPLE", ttbar) + plotBaseName, File_preselection
    )  # MC TTbar
    h_ZJets_amcatnlo_dyjets = GetHisto(
        thisHistName.replace("SAMPLE", zjet) + plotBaseName,
        File_preselection,
    )
    h_WJets_amcatnlo_dyjets = GetHisto(
        thisHistName.replace("SAMPLE", wjet) + plotBaseName,
        File_preselection,
    )
    # h_dyjets_amcatnlo_dyjets = GetHisto("histo1D__WJet_amcatnlo_Inc__"+plotBaseName, File_preselection)
    h_SingleTop_dyjets = GetHisto(
        thisHistName.replace("SAMPLE", singletop) + plotBaseName, File_preselection
    )
    h_PhotonJets_dyjets = GetHisto(
        thisHistName.replace("SAMPLE", photonjets) + plotBaseName,
        File_preselection,
    )
    # h_Diboson_dyjets = GetHisto("histo1D__DIBOSON__"+plotBaseName, File_preselection)
    h_Diboson_dyjets = GetHisto(
        thisHistName.replace("SAMPLE", diboson) + plotBaseName,
        File_preselection,
    )

    # DATA
    h_DATA_dyjets = GetHisto(
        thisHistName.replace("SAMPLE", "DATA") + plotBaseName, File_preselection
    )  # DATA
    # QCD
    h_QCD = GetHisto(
        thisHistName.replace("SAMPLE", qcd)
        # + plotBaseName.replace("_btagSFDownShift", "").replace("_btagSFUpShift", ""),
        + plotBaseName,
        File_QCD_preselection,
    )
    # h_QCD_wjets = GetHisto("histo1D__QCD_EMEnriched__"+plotBaseName,File_QCD_preselection)

    plotDYJets = Plot()
    plotDYJets.histoDATA = h_DATA_dyjets
    plotDYJets.histoMCall = h_ALLBKG_powheg_dyjets
    plotDYJets.histoTTbar = h_TTbar_powheg_dyjets
    # plotDYJets.histoQCD = h_QCD_dyjets
    plotDYJets.histoQCD = h_QCD
    plotDYJets.histoZJet = h_ZJets_amcatnlo_dyjets
    plotDYJets.histoWJet = h_WJets_amcatnlo_dyjets

    plotDYJets.histoSingleTop = h_SingleTop_dyjets
    plotDYJets.histoPhotonJets = h_PhotonJets_dyjets
    plotDYJets.histoDiboson = h_Diboson_dyjets
    plotDYJets.xmin = meeMinDYJets
    # plotDYJets.xmax = h_TTbar_amcatnlo_wjets.GetXaxis().GetXmax()
    plotDYJets.xmax = meeMaxDYJets
    # plotDYJets.name = "DYJetsRescale"
    plotDYJets.name = plotBaseName + "_DYJets"
    plotDYJets.fileXsectionNoRescale = xsectionFile
    plotDYJets.xminplot = 0
    plotDYJets.xmaxplot = 2000
    plotDYJets.yminplot = 0
    plotDYJets.ymaxplot = 2000
    plotDYJets.datasetName = zjetDatasetName
    # plotDYJets.datasetName = "DYJetsToLNu_Pt.+Tune"
    # plot0.datasetName = "Z.+Jets_Pt.+alpgen"
    # example: this match with /Z3Jets_Pt300to800-alpgen/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
    plotsDYJets.append(plotDYJets)

# -----------------------------------------------------------------------------------


# ############ USER CODE - END ################################################
##############################################################################


# --- Generate and print the plots from the list 'plots' define above

# --- Output files
fileps = "allPlots_calc_dyJetsAndTTBarRescale_And_xsecFile.ps"

# --- Generate and print the plots from the list 'plots' define above
# c = TCanvas()
# c.Print(fileps+"[")
for idx, plot in enumerate(plotsTTBar):
    print "plot:",plot.name
    CalculateRescaleFactor(plot, plotsDYJets[idx], fileps)
# c.Print(fileps+"]")
# os.system('ps2pdf '+fileps)
