#!/usr/bin/env python

##############################################################################
## USER CODE IS TOWARD THE END OF THE FILE
##############################################################################

##############################################################################
############# DON'T NEED TO MODIFY ANYTHING HERE - BEGIN #####################

#---Import
import sys
import string
from optparse import OptionParser
import os.path
from ROOT import *
import re
from array import array
import copy
import math

#--- ROOT general options
gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetPalette(1)
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1);
gStyle.SetPadTickY(1);
#--- TODO: WHY IT DOES NOT LOAD THE DEFAULT ROOTLOGON.C ? ---#

def GetFile(filename):
    file = TFile(filename)
    if( not file):
        print "ERROR: file " + filename + " not found"
        print "exiting..."
        sys.exit()
    return file


def GetHisto( histoName , file, scale = 1 ):
    histo = file.Get( histoName )
    if( not histo):
        print "ERROR: histo " + histoName + " not found in " + file.GetName()
        print "exiting..."
        sys.exit()
    new = copy.deepcopy(histo)
    if(scale!=1):
        new.Scale(scale)
    return new

def GetIntegralTH1( histo, xmin, xmax):
    #get integral
    axis = histo.GetXaxis()
    bmin = axis.FindBin(xmin)
    bmax = axis.FindBin(xmax)
    bminResidual = histo.GetBinContent(bmin)*(xmin-axis.GetBinLowEdge(bmin)) / axis.GetBinWidth(bmin)
    bmaxResidual = histo.GetBinContent(bmax)*(axis.GetBinUpEdge(bmax)-xmax) / axis.GetBinWidth(bmax)
    integral = histo.Integral(bmin,bmax) - bminResidual - bmaxResidual
    print 'GetIntegralTH1('+histo.GetName(),xmin,xmax,')=',integral
    return integral

def GetErrorIntegralTH1( histo, xmin, xmax):
    print "## calculating error for integral of histo " + str(histo)
    print "## in the x range [" + str(xmin) + "," + str(xmax) + "]"
    #get integral
    axis = histo.GetXaxis()
    bmin = axis.FindBin(xmin)
    bmax = axis.FindBin(xmax)
    bminResidual = histo.GetBinContent(bmin)*(xmin-axis.GetBinLowEdge(bmin)) / axis.GetBinWidth(bmin)
    bmaxResidual = histo.GetBinContent(bmax)*(axis.GetBinUpEdge(bmax)-xmax) / axis.GetBinWidth(bmax)
    integral = histo.Integral(bmin,bmax) - bminResidual - bmaxResidual
    error = 0
    for bin in range(bmin, bmax+1):
	#print "bin: " +str(bin)
        if(bin==bmax and bmaxResidual==histo.GetBinContent(bmax)): # skip last bin if out of range
            print "     --> skip bin: " + str(bin)
        else:
            error = error + histo.GetBinError(bin)**2
            #print "error**2 : " + str(error)

    error = math.sqrt(error)
    print  " "
    return error


## The Plot class: add members if needed
class Plot:
    histoDATA    = "" # DATA
    histoTTbar   = "" # MCTTbar
    histoMCall   = "" # MCall
    histoQCD     = "" # QCD
    histoZJet    = ""
    histoWJet    = ""
    histoSingleTop = ""
    histoPhotonJets = ""
    histoDiboson = ""
    xtit         = "" # xtitle
    ytit         = "" # ytitle
    xmin         = "" # set xmin to calculate rescaling factor (-- please take into account the bin size of histograms --)
    xmax         = "" # # set xmax to calculate rescaling factor (-- please take into account the bin size of histograms --)
    xminplot     = "" # min x axis range (need to set both min and max. Leave it as is for full range)
    xmaxplot     = "" # max x axis range (need to set both min and max. Leave it as is for full range)
    yminplot     = "" # min y axis range (need to set both min and max. Leave it as is for full range)
    ymaxplot     = "" # max y axis range (need to set both min and max. Leave it as is for full range)
    lpos         = "" # legend position (default = top-right, option="bottom-center", "top-left")
    #    xlog         = "" # log scale of X axis (default = no, option="yes") ### IT SEEMS IT DOES NOT WORK
    ylog         = "" # log scale of Y axis (default = no, option="yes")
    #rebin       = "" # rebin x axis (default = 1, option = set it to whatever you want )
    name         = "" # name of the final plots
    lint         = "2.6 fb^{-1}" # integrated luminosity of the sample ( example "10 pb^{-1}" )
    fileXsectionNoRescale = "" #cross section file (with no rescale
    datasetName = "" # string for pattern recognition of dataset name (rescaling will be done only on matched datasets)

    def CheckMCDataConsistency(self):
        #checks
        if(self.histoMCall.GetNbinsX()!=self.histoDATA.GetNbinsX()):
            print "WARNING! number of bins is different between DATA and MC"
            print "exiting..."
            sys.exit()
        if(self.histoMCall.GetBinWidth(1)!=self.histoDATA.GetBinWidth(1)):
            print "WARNING! bin width is different between DATA and MC"
            print "exiting..."
            sys.exit()


def GetRTTBarWJets(plotObjTTBar, plotObjWJets, randomize=False):

    #integrals: ttbar
    integralDATA_ttbar = GetIntegralTH1(plotObjTTBar.histoDATA,plotObjTTBar.xmin,plotObjTTBar.xmax)
    ERRintegralDATA_ttbar = GetErrorIntegralTH1(plotObjTTBar.histoDATA,plotObjTTBar.xmin,plotObjTTBar.xmax)
    integralMCall_ttbar = GetIntegralTH1(plotObjTTBar.histoMCall,plotObjTTBar.xmin,plotObjTTBar.xmax)
    ERRintegralMCall_ttbar = GetErrorIntegralTH1(plotObjTTBar.histoMCall,plotObjTTBar.xmin,plotObjTTBar.xmax)
    integralTTbar_ttbar = GetIntegralTH1(plotObjTTBar.histoTTbar,plotObjTTBar.xmin,plotObjTTBar.xmax)
    ERRintegralTTbar_ttbar = GetErrorIntegralTH1(plotObjTTBar.histoTTbar,plotObjTTBar.xmin,plotObjTTBar.xmax)
    integralWJets_ttbar = GetIntegralTH1(plotObjTTBar.histoWJet,plotObjTTBar.xmin,plotObjTTBar.xmax)
    ERRintegralWJets_ttbar = GetErrorIntegralTH1(plotObjTTBar.histoWJet,plotObjTTBar.xmin,plotObjTTBar.xmax)
    integralQCD_ttbar = GetIntegralTH1(plotObjTTBar.histoQCD,plotObjTTBar.xmin,plotObjTTBar.xmax)
    ERRintegralQCD_ttbar = GetErrorIntegralTH1(plotObjTTBar.histoQCD,plotObjTTBar.xmin,plotObjTTBar.xmax)
    #contamination from other backgrounds (except TTbar and WJets) in the integral range [QCD is not in MCall]
    integralMCothers_ttbar = integralMCall_ttbar - integralTTbar_ttbar - integralWJets_ttbar
    ERRintegralMCothers_ttbar = math.sqrt(ERRintegralMCall_ttbar**2 + ERRintegralTTbar_ttbar**2)

    # integrals: wjets
    integralDATA_wjets = GetIntegralTH1(plotObjWJets.histoDATA,plotObjWJets.xmin,plotObjWJets.xmax)
    ERRintegralDATA_wjets = GetErrorIntegralTH1(plotObjWJets.histoDATA,plotObjWJets.xmin,plotObjWJets.xmax)
    integralMCall_wjets = GetIntegralTH1(plotObjWJets.histoMCall,plotObjWJets.xmin,plotObjWJets.xmax)
    ERRintegralMCall_wjets = GetErrorIntegralTH1(plotObjWJets.histoMCall,plotObjWJets.xmin,plotObjWJets.xmax)
    integralTTbar_wjets = GetIntegralTH1(plotObjWJets.histoTTbar,plotObjWJets.xmin,plotObjWJets.xmax)
    ERRintegralTTbar_wjets = GetErrorIntegralTH1(plotObjWJets.histoTTbar,plotObjWJets.xmin,plotObjWJets.xmax)
    integralWJets_wjets = GetIntegralTH1(plotObjWJets.histoWJet,plotObjWJets.xmin,plotObjWJets.xmax)
    ERRintegralWJets_wjets = GetErrorIntegralTH1(plotObjWJets.histoWJet,plotObjWJets.xmin,plotObjWJets.xmax)
    integralQCD_wjets = GetIntegralTH1(plotObjWJets.histoQCD,plotObjWJets.xmin,plotObjWJets.xmax)
    ERRintegralQCD_wjets = GetErrorIntegralTH1(plotObjWJets.histoQCD,plotObjWJets.xmin,plotObjWJets.xmax)
    #contamination from other backgrounds (except WJets and TTBar) in the integral range [QCD is not in MCall]
    integralMCothers_wjets = integralMCall_wjets - integralWJets_wjets - integralTTbar_wjets
    ERRintegralMCothers_wjets = math.sqrt(ERRintegralMCall_wjets**2 + ERRintegralWJets_wjets**2)
    contamination_wjets = (integralMCothers_wjets+integralTTbar_wjets+integralQCD_wjets) / (integralMCall_wjets+integralQCD_wjets)

    # randomize
    trand = TRandom3(0)
    integralDATA_ttbar = trand.Gaus(integralDATA_ttbar,ERRintegralDATA_ttbar)
    integralMCall_ttbar = trand.Gaus(integralMCall_ttbar,ERRintegralMCall_ttbar)
    integralTTbar_ttbar = trand.Gaus(integralTTbar_ttbar,ERRintegralTTbar_ttbar)
    integralWJets_ttbar = trand.Gaus(integralWJets_ttbar,ERRintegralWJets_ttbar)
    integralQCD_ttbar = trand.Gaus(integralQCD_ttbar,ERRintegralQCD_ttbar)
    integralMCothers_ttbar = integralMCall_ttbar - integralTTbar_ttbar - integralWJets_ttbar
    #
    integralDATA_wjets = trand.Gaus(integralDATA_wjets,ERRintegralDATA_wjets)
    integralMCall_wjets = trand.Gaus(integralMCall_wjets,ERRintegralMCall_wjets)
    integralTTbar_wjets = trand.Gaus(integralTTbar_wjets,ERRintegralTTbar_wjets)
    integralWJets_wjets = trand.Gaus(integralWJets_wjets,ERRintegralWJets_wjets)
    integralQCD_wjets = trand.Gaus(integralQCD_wjets,ERRintegralQCD_wjets)

    # solve the system of equations
    # (1) --> wjets
    # (2) --> ttbar
    rTTBar = integralDATA_wjets*integralWJets_ttbar - (integralMCothers_wjets+integralQCD_wjets)*integralWJets_ttbar
    rTTBar+= (integralMCothers_ttbar+integralQCD_ttbar-integralDATA_ttbar)*integralWJets_wjets
    print integralTTbar_wjets,integralWJets_ttbar,integralTTbar_ttbar,integralWJets_wjets
    try:
      rTTBar/= (integralTTbar_wjets*integralWJets_ttbar - integralTTbar_ttbar*integralWJets_wjets)
    except ZeroDivisionError:
      print 'ERROR: ZeroDivisionError: one or more of the integrals above are zero'
      return -1,-1
    rWJets = integralDATA_wjets*integralTTbar_ttbar - (integralMCothers_wjets+integralQCD_wjets)*integralTTbar_ttbar
    rWJets+= (integralMCothers_ttbar+integralQCD_ttbar-integralDATA_ttbar)*integralTTbar_wjets
    rWJets/= (integralTTbar_ttbar*integralWJets_wjets-integralTTbar_wjets*integralWJets_ttbar)

    return rTTBar,rWJets

def CalculateRescaleFactor(plotObjTTBar, plotObjWJets, fileps):
    #calculate rescaling factor for Z/gamma+jet background and create new cross section file
    canvas = TCanvas()

    plotObjTTBar.CheckMCDataConsistency()
    plotObjWJets.CheckMCDataConsistency()

    # do the following N times, so we can calculate the stat. uncertainty
    N=100
    # NB: the commented code below makes a nice progress bar but causes the dict to be undefined...
    steps = N
    print 'Randomizing histos and calculating scale factors:'
    progressString = '0% ['+' '*steps+'] 100%'
    print progressString,
    print '\b'*(len(progressString)-3),
    sys.stdout.flush()
    rTTBarList = []
    rWJetsList = []
    for i in range(0,N):
        print ''
        rTTBar,rWJets = GetRTTBarWJets(plotObjTTBar,plotObjWJets,True)
        rTTBarList.append(rTTBar)
        rWJetsList.append(rWJets)
        print '\b.',
        sys.stdout.flush()
    print '\b] 100%'

    ttMean = 0
    for rtt in rTTBarList:
      ttMean+=rtt
    ttMean/=N
    #
    wMean = 0
    for rw in rWJetsList:
      wMean+=rw
    wMean/=N
    
    rTTBarSigma = 0
    for rtt in rTTBarList:
      rTTBarSigma+=pow(rtt-ttMean,2)
    rTTBarSigma/=N
    rTTBarSigma = math.sqrt(rTTBarSigma)
    #
    rWJetsSigma = 0
    for rw in rWJetsList:
      rWJetsSigma+=pow(rw-wMean,2)
    rWJetsSigma/=N
    rWJetsSigma = math.sqrt(rWJetsSigma)
    
    #integrals: ttbar
    integralDATA_ttbar = GetIntegralTH1(plotObjTTBar.histoDATA,plotObjTTBar.xmin,plotObjTTBar.xmax)
    ERRintegralDATA_ttbar = GetErrorIntegralTH1(plotObjTTBar.histoDATA,plotObjTTBar.xmin,plotObjTTBar.xmax)
    integralMCall_ttbar = GetIntegralTH1(plotObjTTBar.histoMCall,plotObjTTBar.xmin,plotObjTTBar.xmax)
    ERRintegralMCall_ttbar = GetErrorIntegralTH1(plotObjTTBar.histoMCall,plotObjTTBar.xmin,plotObjTTBar.xmax)
    integralTTbar_ttbar = GetIntegralTH1(plotObjTTBar.histoTTbar,plotObjTTBar.xmin,plotObjTTBar.xmax)
    print 'plotObjTTBar.histoTTbar=',plotObjTTBar.histoTTbar
    ERRintegralTTbar_ttbar = GetErrorIntegralTH1(plotObjTTBar.histoTTbar,plotObjTTBar.xmin,plotObjTTBar.xmax)
    integralWJets_ttbar = GetIntegralTH1(plotObjTTBar.histoWJet,plotObjTTBar.xmin,plotObjTTBar.xmax)
    ERRintegralWJets_ttbar = GetErrorIntegralTH1(plotObjTTBar.histoWJet,plotObjTTBar.xmin,plotObjTTBar.xmax)
    integralQCD_ttbar = GetIntegralTH1(plotObjTTBar.histoQCD,plotObjTTBar.xmin,plotObjTTBar.xmax)
    ERRintegralQCD_ttbar = GetErrorIntegralTH1(plotObjTTBar.histoQCD,plotObjTTBar.xmin,plotObjTTBar.xmax)
    #contamination from other backgrounds (except TTbar and WJets) in the integral range [QCD is not in MCall]
    integralMCothers_ttbar = integralMCall_ttbar - integralTTbar_ttbar - integralWJets_ttbar
    ERRintegralMCothers_ttbar = math.sqrt(ERRintegralMCall_ttbar**2 + ERRintegralTTbar_ttbar**2)
    try:
      contamination_ttbar = (integralMCothers_ttbar+integralWJets_ttbar+integralQCD_ttbar) / (integralMCall_ttbar+integralQCD_ttbar)
    except ZeroDivisionError:
      print 'ERROR: ZeroDivisionError: integralMCall_ttbar+integralQCD_ttbar is zero; integralMCall_ttbar=',integralMCall_ttbar,'integralQCD_ttbar=',integralQCD_ttbar
      contamination_ttbar = -1

    # integrals: wjets
    integralDATA_wjets = GetIntegralTH1(plotObjWJets.histoDATA,plotObjWJets.xmin,plotObjWJets.xmax)
    ERRintegralDATA_wjets = GetErrorIntegralTH1(plotObjWJets.histoDATA,plotObjWJets.xmin,plotObjWJets.xmax)
    integralMCall_wjets = GetIntegralTH1(plotObjWJets.histoMCall,plotObjWJets.xmin,plotObjWJets.xmax)
    ERRintegralMCall_wjets = GetErrorIntegralTH1(plotObjWJets.histoMCall,plotObjWJets.xmin,plotObjWJets.xmax)
    integralTTbar_wjets = GetIntegralTH1(plotObjWJets.histoTTbar,plotObjWJets.xmin,plotObjWJets.xmax)
    print 'plotObjWJets.histoTTbar=',plotObjWJets.histoTTbar
    ERRintegralTTbar_wjets = GetErrorIntegralTH1(plotObjWJets.histoTTbar,plotObjWJets.xmin,plotObjWJets.xmax)
    integralWJets_wjets = GetIntegralTH1(plotObjWJets.histoWJet,plotObjWJets.xmin,plotObjWJets.xmax)
    ERRintegralWJets_wjets = GetErrorIntegralTH1(plotObjWJets.histoWJet,plotObjWJets.xmin,plotObjWJets.xmax)
    integralQCD_wjets = GetIntegralTH1(plotObjWJets.histoQCD,plotObjWJets.xmin,plotObjWJets.xmax)
    ERRintegralQCD_wjets = GetErrorIntegralTH1(plotObjWJets.histoQCD,plotObjWJets.xmin,plotObjWJets.xmax)
    #contamination from other backgrounds (except WJets and TTBar) in the integral range [QCD is not in MCall]
    integralMCothers_wjets = integralMCall_wjets - integralWJets_wjets - integralTTbar_wjets
    ERRintegralMCothers_wjets = math.sqrt(ERRintegralMCall_wjets**2 + ERRintegralWJets_wjets**2)
    contamination_wjets = (integralMCothers_wjets+integralTTbar_wjets+integralQCD_wjets) / (integralMCall_wjets+integralQCD_wjets)

    # solve the system of equations
    # (1) --> wjets
    # (2) --> ttbar
    rTTBar = integralDATA_wjets*integralWJets_ttbar - (integralMCothers_wjets+integralQCD_wjets)*integralWJets_ttbar
    rTTBar+= (integralMCothers_ttbar+integralQCD_ttbar-integralDATA_ttbar)*integralWJets_wjets
    try:
      rTTBar/= (integralTTbar_wjets*integralWJets_ttbar - integralTTbar_ttbar*integralWJets_wjets)
    except ZeroDivisionError:
      print 'ERROR in rTTBar: ZeroDivisionError: one or more of the integrals is zero:'
      print integralTTbar_wjets,integralWJets_ttbar,'-',integralTTbar_ttbar,integralWJets_wjets
    rWJets = integralDATA_wjets*integralTTbar_ttbar - (integralMCothers_wjets+integralQCD_wjets)*integralTTbar_ttbar
    rWJets+= (integralMCothers_ttbar+integralQCD_ttbar-integralDATA_ttbar)*integralTTbar_wjets
    try:
      rWJets/= (integralTTbar_ttbar*integralWJets_wjets-integralTTbar_wjets*integralWJets_ttbar)
    except ZeroDivisionError:
      print 'ERROR in rWJets: ZeroDivisionError: one or more of the integrals is zero:'
      print integralTTbar_ttbar,integralWJets_wjets,'-',integralTTbar_wjets,integralWJets_ttbar


# FIXME
    ##draw histo
    #self.histoMCall.SetFillColor(kBlue)
    #self.histoDATA.SetMarkerStyle(20)

    #self.histoMCall.Draw("HIST")
    #self.histoDATA.Draw("psame")
    #self.histoMCall.GetXaxis().SetRangeUser(self.xminplot,self.xmaxplot)
    #self.histoMCall.GetYaxis().SetRangeUser(self.yminplot,self.ymaxplot)

    #canvas.Update()
    #gPad.RedrawAxis()
    #gPad.Modified()
    ##canvas.SaveAs(self.name + ".eps","eps")
    ##canvas.SaveAs(self.name + ".pdf","pdf")
    #canvas.Print(fileps)
    #canvas.Print(self.name + ".C")
    ## make root file
    #tfile = TFile(self.name+'.root','recreate')
    #tfile.cd()
    #self.histoDATA.Write()
    #self.histoTTbar.Write()
    #self.histoMCall.Write()
    #self.histoQCD.Write()
    #self.histoZJet.Write()
    #self.histoWJet.Write()
    #self.histoSingleTop.Write()
    #self.histoPhotonJets.Write()
    #self.histoDiboson.Write()
    #tfile.Close()

    #printout
    print
    print " TTBar "
    print "######################################## "
    print "name:                         " + plotObjTTBar.name
    print "integral range:               " + str(plotObjTTBar.xmin) + " < MTenu < " + str(plotObjTTBar.xmax) + " GeV/c2"
    print "integral MC All:              "   + str( integralMCall_ttbar ) + " +/- " + str( ERRintegralMCall_ttbar )
    print "integral QCD:                 "   + str( integralQCD_ttbar ) + " +/- " + str( ERRintegralQCD_ttbar )
    print "integral MC TTbar:            "   + str( integralTTbar_ttbar) + " +/- " + str( ERRintegralTTbar_ttbar )
    print "integral MC WJets:            "   + str( integralWJets_ttbar) + " +/- " + str( ERRintegralWJets_ttbar )
    print "integral MC other:            "   + str( integralMCothers_ttbar) + " +/- " + str( ERRintegralMCothers_ttbar )
    print "rescaled integral MC TTbar:   "   + str( rTTBar*integralTTbar_ttbar) + " +/- " + str( rTTBar*ERRintegralTTbar_ttbar )
    print "rescaled integral MC WJets:   "   + str( rWJets*integralWJets_ttbar) + " +/- " + str( rWJets*ERRintegralWJets_ttbar )
    print "integral DATA:                "   + str( integralDATA_ttbar ) + " +/- " + str( ERRintegralDATA_ttbar )
    print "contribution from other bkgs (except TTbar): " + str(contamination_ttbar*100) + "%"
    #print "integral DATA (corrected for contribution from other bkgs): "  + str( integralDATAcorr_ttbar ) + " +/- " + str( ERRintegralDATAcorr_ttbar )
    print "rescale factor for TTbar background: " + str(rTTBar) + " +\- " + str(rTTBarSigma)
    #print "systematical uncertainty of TTbar background modeling: " + str(relERRrescale*100) + "%"
    print "######################################## "
    print

    print
    print " WJets "
    print "######################################## "
    print "name:                         " + plotObjWJets.name
    print "integral range:               " + str(plotObjWJets.xmin) + " < MTenu < " + str(plotObjWJets.xmax) + " GeV/c2"
    print "integral MC All:              "   + str( integralMCall_wjets ) + " +/- " + str( ERRintegralMCall_wjets )
    print "integral QCD:                 "   + str( integralQCD_wjets ) + " +/- " + str( ERRintegralQCD_wjets )
    print "integral MC TTbar:            "   + str( integralTTbar_wjets) + " +/- " + str( ERRintegralTTbar_wjets )
    print "integral MC WJets:            "   + str( integralWJets_wjets) + " +/- " + str( ERRintegralWJets_wjets )
    print "integral MC other:            "   + str( integralMCothers_wjets) + " +/- " + str( ERRintegralMCothers_wjets )
    print "rescaled integral MC TTbar:   "   + str( rTTBar*integralTTbar_wjets) + " +/- " + str( rTTBar*ERRintegralTTbar_wjets )
    print "rescaled integral MC WJets:   "   + str( rWJets*integralWJets_wjets) + " +/- " + str( rWJets*ERRintegralWJets_wjets )
    print "integral DATA:                "   + str( integralDATA_wjets ) + " +/- " + str( ERRintegralDATA_wjets )
    print "contribution from other bkgs (except wjets): " + str(contamination_wjets*100) + "%"
    #print "integral DATA (corrected for contribution from other bkgs): "  + str( integralDATAcorr_wjets ) + " +/- " + str( ERRintegralDATAcorr_wjets )
    print "rescale factor for WJets background: " + str(rWJets) + " +\- " + str(rWJetsSigma)
    print "######################################## "
    print

    #create new cross section file -- ttbar
    originalFileName = string.split( string.split(plotTTbar.fileXsectionNoRescale, "/" )[-1], "." ) [0]
    ttbarFileName = originalFileName + "_" + plotTTbar.name +".txt"
    os.system('rm -f '+ ttbarFileName)
    outputFile = open(ttbarFileName,'w')

    for line in open( plotTTbar.fileXsectionNoRescale ):
        line = string.strip(line,"\n")
        lineNoComments = line.split("#")[0] # strip off anything after any '#' if present
        # ignore empty lines
        if len(lineNoComments) <= 0:
          print >> outputFile, line
          continue
        #FIXME this doesn't support keeping comments at the end of lines

        if( re.search(plotTTbar.datasetName, line) ):
            list = re.split( '\s+' , line  )
            newline = str(list[0]) + "    "  + str("%.6f" % (float(list[1])*float(rTTBar)) )
            print >> outputFile, newline
        else:
            print >> outputFile, line

    outputFile.close()
    print "New xsection file (after TTbar rescaling) is: " + ttbarFileName
    print " "

    #create new cross section file -- WJets
    originalFileName = string.split( string.split(ttbarFileName, "/" )[-1], "." ) [0]
    newFileName = originalFileName + "_" + plotWJets.name +".txt"
    os.system('rm -f '+ newFileName)
    outputFile = open(newFileName,'w')

    for line in open( ttbarFileName ):
        line = string.strip(line,"\n")
        lineNoComments = line.split("#")[0] # strip off anything after any '#' if present
        # ignore empty lines
        if len(lineNoComments) <= 0:
          print >> outputFile, line
          continue

        if( re.search(plotWJets.datasetName, line) ):
            list = re.split( '\s+' , line  )
            newline = str(list[0]) + "    "  + str("%.6f" % (float(list[1])*float(rWJets)) )
            print >> outputFile, newline
        else:
            print >> outputFile, line

    outputFile.close()
    print "New xsection file (after WJets rescaling) is: " + newFileName
    print " "


############# DON'T NEED TO MODIFY ANYTHING HERE - END #######################
##############################################################################


##############################################################################
############# USER CODE - BEGIN ##############################################

#--- Input files
#preselection
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_onRSK_local_nov1_addStSFplots_allDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_orig.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_onRSK_local_nov1_addStSFplots_allDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_QCD_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_local_nov18_addStSFplots_allDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")

#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_local_nov18_addStSFplots_ICHEPDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_local_nov18_addStSFplots_ICHEPDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_QCD_plots.root")

#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_dec3_ICHEPDataExcludeEarlyRunsAndMC_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_dec3_ICHEPDataExcludeEarlyRunsAndMC_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_QCD_plots.root")

#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_dec14_ICHEPDataExcludeEarlyRunsAndMC_rereco_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_dec14_ICHEPDataExcludeEarlyRunsAndMC_rereco_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_QCD_plots.root")

#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jan20_rereco_stitch120_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
##File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jan20_rereco_stitch120_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_QCD_jan22_rereco_enujj2012FinSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")

#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb28_recoHeepSF_rereco_stitchDYW_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_QCD_feb27_rereco_enujj2012FinSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")

# unscaled
# no topPt reweight
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_mar28_recoHeepSF_rereco_stitchDYW_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
# with reweight
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_mar30_topPtWeight_recoHeepSF_rereco_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_apr11_ele27wptightOREle115_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_feb28_recoHeepSF_rereco_stitchDYW_ele27wptightEta2p1CurveMC/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_jun1_ele27wptightOREle115ORPhoton175/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_sep25_btags__enujjOptFinalSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_oct16_newFR_updatedFinalSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_oct19_reminiAODFR_updatedFinalSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_oct23_updateFR_updatedFinalSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_oct24_updateFRAndOct21QCDRSK_updatedFinalSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_rsk_nov10_failHEEP_updateFRSCEt_preselOnly/output_cutTable_lq_enujj_MT_QCD_preselOnly/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk_nov17_highMejPlots_preselOnly/output_cutTable_lq_enujj_MT_QCD_preselOnly/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk237_jan17_preselOnly/output_cutTable_lq_enujj_MT_QCD_preselOnly/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk237_jan24_gsfEtaCheck_preselOnly/output_cutTable_lq_enujj_MT_QCD_preselOnly/analysisClass_lq_enujj_QCD_plots.root")
# MET 100
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk237_jan31_gsfEtaCheck_MET100/output_cutTable_lq_enujj_MT_QCD_preselOnly_MET100/analysisClass_lq_enujj_QCD_plots.root")
# with Pt(e,met)
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk237_feb4_gsfEtaCheck_MET100_PtEMET70/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk237_feb5_gsfEtaCheck_MET100_PtEMET70_prevCutHists/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk237_feb7_gsfEtaCheck_MET100_PtEMET70_addHists/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_feb10_bugfix/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_feb14_dPhiEleMET0p8/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_mar16_fixMuons/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#
File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_mar26_addMTPlots/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")


# unscaled
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_may8_ele27wptightOREle115_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
# unscaled ele27 from Bibhu
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_may11_ele27wptight_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jun1_ele27wptightOREle115ORPhoton175_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jul4_ele27wptightOREle115ORPhoton175_enujjOptFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root.jul4")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jul4_ele27wptightOREle115ORPhoton175_enujjOptFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_sep26_bTagFix_enujjPowhegOptFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")

#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_oct3_bTagFix_enujjPowhegOptFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")

#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_oct5_noTrigEffMC_enujjPowhegOptFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_oct6_finerBinnedTrigEff_enujjPowhegOptFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_oct15_extendMTRange_updatedFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_nov10_reprodNominal/output_cutTable_lq_enujj_MT_preselOnly/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_nov22_fixTrigEff_usePtHeep_preselOnly/output_cutTable_lq_enujj_MT_preselOnly/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_dec13_fixTrigEff_usePtHeep_preselOnly/output_cutTable_lq_enujj_MT_preselOnly/analysisClass_lq_enujj_MT_plots_unscaled.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jan24_v237_preselOnly/output_cutTable_lq_enujj_MT_preselOnly/analysisClass_lq_enujj_MT_plots_unscaled.root")
# MET100
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb1_v237_MET100/output_cutTable_lq_enujj_MT_preselOnly_MET100/analysisClass_lq_enujj_MT_plots_unscaled.root")
# with Pt(e,MET) > 70
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb4_v237_MET100_PtEMET70/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb4_v237_MET100_PtEMET70_prevCutHists/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb7_v237_MET100_PtEMET70_addHists/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb10_v237_bugfix/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb14_dPhiEleMet0p8/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb19_dPhiEleMet0p8_btagSysts/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb20_dPhiEleMet0p8_newSingTop/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb23_dPhiEleMet0p8_fixMTPlots/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_mar5_removeTopPtReweight/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_mar16_fixMuons/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#
File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_mar26_addMTPlots/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")

histBaseNames = []
mtRanges = []
histBaseNames.append('MTenu_50_110_BJETBIN') # nominal
mtRanges.append([50,110])
#
#histBaseNames.append('MTenu_50_110_BJETBIN_btagSFDownShift')
#mtRanges.append([50,110])
#histBaseNames.append('MTenu_50_110_BJETBIN_btagSFUpShift')
#mtRanges.append([50,110])
##
#histBaseNames.append('MTenu_50_110_BJETBIN_LQ300')
#mtRanges.append([50,110])
#histBaseNames.append('MTenu_50_110_BJETBIN_LQ400')
#mtRanges.append([50,110])
#histBaseNames.append('MTenu_50_110_BJETBIN_LQ500')
#mtRanges.append([50,110])
#histBaseNames.append('MTenu_50_110_BJETBIN_LQ600')
#mtRanges.append([50,110])
#histBaseNames.append('MTenu_50_110_BJETBIN_LQ700')
#mtRanges.append([50,110])
#histBaseNames.append('MTenu_50_110_BJETBIN_LQ800') # no stats
#mtRanges.append([50,110])
#
#histBaseNames.append('MTenu_50_110_BJETBIN_sT300To500_PAS')
#mtRanges.append([50,110])
#histBaseNames.append('MTenu_50_110_BJETBIN_sT500To750_PAS')
#mtRanges.append([50,110])
#histBaseNames.append('MTenu_50_110_BJETBIN_sT750To1250_PAS')
#mtRanges.append([50,110])
#histBaseNames.append('MTenu_50_110_BJETBIN_sT1250ToInf_PAS')
#mtRanges.append([50,110])
#
# To do the MT variations, need to change the xmin/xmax (i.e., MT) integration ranges
histBaseNames.append('MTenu_BJETBIN_MT200To400_PAS')
mtRanges.append([200,400])
histBaseNames.append('MTenu_BJETBIN_MT400To600_PAS')
mtRanges.append([400,600])
histBaseNames.append('MTenu_BJETBIN_MT600To900_PAS')
mtRanges.append([600,900])
histBaseNames.append('MTenu_BJETBIN_MT900ToInf_PAS')
mtRanges.append([900,2000])
#histBaseNames.append('MTenu_70_150_BJETBIN')
#mtRanges.append([70,150])
#histBaseNames.append('MTenu_110_190_BJETBIN')
#mtRanges.append([110,190])

plotsTTBar = []
plotsWJets = []
useMGHT = False
useAMCAtNLOTTBar = False # if False use powheg

histNameDefault = "histo1D__SAMPLE__"
#histNameReleaseMee = "histo1D__SAMPLE__cutHisto_allOtherCuts___________"
histNameAllPrevCuts = "histo1D__SAMPLE__cutHisto_allPreviousCuts________"

#--- Rescaling of W + jet and ttbar+jets background

for idx,histBaseName in enumerate(histBaseNames):
  #-----------------------------------------
  # for ttbar-enriched region
  #plotBaseName = 'MTenu_50_110_Njet_gte4'
  #plotBaseName = 'MTenu_50_110_gteOneBtaggedJet'
  #plotBaseName = 'MTenu_70_150_gteOneBtaggedJet'
  #plotBaseName = 'MTenu_110_190_gteOneBtaggedJet'
  #plotBaseName= 'MTenu_50_110_gteOneBtaggedJet_LQ800'
  thisHistName = histNameDefault

  mtMin = mtRanges[idx][0]
  mtMax = mtRanges[idx][1]

  plotBaseName = histBaseName.replace('BJETBIN','gteOneBtaggedJet')
  print 'for TTBar, using plotBaseName:',plotBaseName
  
  ## MG HT BKG
  #h_ALLBKG_HT_ttbar = GetHisto("histo1D__ALLBKG_MG_HT__"+plotBaseName, File_preselection) # MC all
  # amc@NLO BKG
  #h_ALLBKG_amcatnlo_ttbar = GetHisto("histo1D__ALLBKG_amcAtNLOIncTTBar_ZJetWJetPt__"+plotBaseName, File_preselection) # MC all
  #h_ALLBKG_amcatnlo_ttbar = GetHisto(thisHistName.replace('SAMPLE','ALLBKG_amcAtNLOIncTTBar_ZJetWJetPt')+plotBaseName, File_preselection) # MC all
  # powheg ttbar
  #h_ALLBKG_powheg_ttbar = GetHisto("histo1D__ALLBKG_powhegTTBar_ZJetWJetPt__"+plotBaseName, File_preselection) # MC all
  h_ALLBKG_powheg_ttbar = GetHisto(thisHistName.replace('SAMPLE',"ALLBKG_powhegTTBar_ZJetWJetPt_amcAtNLODiboson")+plotBaseName, File_preselection) # MC all
  
  
  #h_TTbar_MG_ttbar = GetHisto("histo1D__TTbar_Madgraph__"+plotBaseName, File_preselection) # MC TTbar
  #h_ZJets_MGHT_ttbar = GetHisto("histo1D__ZJet_Madgraph_HT__"+plotBaseName, File_preselection)
  #h_WJets_MGHT_ttbar = GetHisto("histo1D__WJet_Madgraph_HT__"+plotBaseName, File_preselection)
  
  #h_TTbar_amcatnlo_ttbar = GetHisto(thisHistName.replace('SAMPLE',"TTbar_amcatnlo_Inc")+plotBaseName, File_preselection) # MC TTbar
  # powheg ttbar
  h_TTbar_powheg_ttbar = GetHisto(thisHistName.replace('SAMPLE',"TTbar_powheg")+plotBaseName, File_preselection) # MC TTbar
  h_ZJets_amcatnlo_ttbar = GetHisto(thisHistName.replace('SAMPLE',"ZJet_amcatnlo_ptBinned")+plotBaseName, File_preselection)
  h_WJets_amcatnlo_ttbar = GetHisto(thisHistName.replace('SAMPLE',"WJet_amcatnlo_ptBinned")+plotBaseName, File_preselection)
  #h_WJets_amcatnlo_ttbar = GetHisto("histo1D__WJet_amcatnlo_Inc__"+plotBaseName, File_preselection)
  h_SingleTop_ttbar = GetHisto(thisHistName.replace('SAMPLE',"SingleTop")+plotBaseName, File_preselection)
  h_PhotonJets_ttbar = GetHisto(thisHistName.replace('SAMPLE',"PhotonJets_Madgraph")+plotBaseName, File_preselection)
  #h_Diboson_ttbar = GetHisto("histo1D__DIBOSON__"+plotBaseName, File_preselection)
  h_Diboson_ttbar = GetHisto(thisHistName.replace('SAMPLE',"DIBOSON_amcatnlo")+plotBaseName, File_preselection)
  
  # DATA
  h_DATA_ttbar = GetHisto(thisHistName.replace('SAMPLE',"DATA")+plotBaseName, File_preselection) #DATA
  # QCD
  h_QCD_DataDriven = GetHisto(thisHistName.replace('SAMPLE',"QCDFakes_DATA")+plotBaseName.replace('_btagSFDownShift','').replace('_btagSFUpShift',''),File_QCD_preselection)
  #h_QCD_ttbar = GetHisto("histo1D__QCD_EMEnriched__"+plotBaseName,File_QCD_preselection)
  
  plotTTbar = Plot()
  plotTTbar.histoDATA = h_DATA_ttbar
  # amc@NLO
  if not useMGHT:
    if useAMCAtNLOTTBar:
      plotTTbar.histoMCall = h_ALLBKG_amcatnlo_ttbar
      plotTTbar.histoTTbar = h_TTbar_amcatnlo_ttbar
    else:
      plotTTbar.histoMCall = h_ALLBKG_powheg_ttbar
      plotTTbar.histoTTbar = h_TTbar_powheg_ttbar
    #plotTTbar.histoQCD = h_QCD_ttbar
    plotTTbar.histoQCD = h_QCD_DataDriven
    plotTTbar.histoZJet = h_ZJets_amcatnlo_ttbar
    plotTTbar.histoWJet = h_WJets_amcatnlo_ttbar
  else:
    # MG HT
    plotTTbar.histoMCall = h_ALLBKG_HT_ttbar
    plotTTbar.histoTTbar = h_TTbar_MG_ttbar
    plotTTbar.histoQCD = h_QCD_ttbar
    plotTTbar.histoZJet = h_ZJets_MGHT_ttbar
    plotTTbar.histoWJet = h_WJets_MGHT_ttbar
  
  plotTTbar.histoSingleTop = h_SingleTop_ttbar
  plotTTbar.histoPhotonJets = h_PhotonJets_ttbar
  plotTTbar.histoDiboson = h_Diboson_ttbar
  plotTTbar.xmin = mtMin
  #plotTTbar.xmax = h_TTbar_amcatnlo_ttbar.GetXaxis().GetXmax()
  plotTTbar.xmax = mtMax
  #plotTTbar.name = "TTbarRescale"
  plotTTbar.name = plotBaseName+'_TTbar'
  plotTTbar.fileXsectionNoRescale = "/afs/cern.ch/user/s/scooper/work/private/cmssw/8011/TestRootNTuplizerRecipe/src/Leptoquarks/analyzer/rootNtupleAnalyzerV2/config/xsection_13TeV_2015.txt"
  plotTTbar.xminplot = 0
  plotTTbar.xmaxplot = 2000
  plotTTbar.yminplot = 0
  plotTTbar.ymaxplot = 2000
  plotTTbar.datasetName = "TT"
  #plot0.datasetName = "DYJetsToLL_M-50_HT.+Tune"
  #plot0.datasetName = "Z.+Jets_Pt.+alpgen"
  # example: this match with /Z3Jets_Pt300to800-alpgen/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
  plotsTTBar.append(plotTTbar)
  
  #-----------------------------------------
  # for wjets-enriched region
  #plotBaseName = 'MTenu_50_110_Njet_lte3'
  #plotBaseName = 'MTenu_50_110_noBtaggedJets'
  #plotBaseName = 'MTenu_70_150_noBtaggedJets'
  #plotBaseName = 'MTenu_110_190_noBtaggedJets'
  #plotBaseName= 'MTenu_50_110_noBtaggedJets_LQ800'
  
  plotBaseName = histBaseName.replace('BJETBIN','noBtaggedJets')
  print 'for WJets, using plotBaseName:',plotBaseName
  
  ## MG HT BKG
  #h_ALLBKG_HT_wjets = GetHisto("histo1D__ALLBKG_MG_HT__"+plotBaseName, File_preselection) # MC all
  # amc@NLO BKG
  #h_ALLBKG_amcatnlo_wjets = GetHisto(thisHistName.replace('SAMPLE',"ALLBKG_amcAtNLOIncTTBar_ZJetWJetPt")+plotBaseName, File_preselection) # MC all
  # powheg
  #h_ALLBKG_powheg_wjets = GetHisto("histo1D__ALLBKG_powhegTTBar_ZJetWJetPt__"+plotBaseName, File_preselection) # MC all
  h_ALLBKG_powheg_wjets = GetHisto(thisHistName.replace('SAMPLE',"ALLBKG_powhegTTBar_ZJetWJetPt_amcAtNLODiboson")+plotBaseName, File_preselection) # MC all
  
  #h_TTbar_MG_wjets = GetHisto("histo1D__TTbar_Madgraph__"+plotBaseName, File_preselection) # MC TTbar
  #h_ZJets_MGHT_wjets = GetHisto("histo1D__ZJet_Madgraph_HT__"+plotBaseName, File_preselection)
  #h_WJets_MGHT_wjets = GetHisto("histo1D__WJet_Madgraph_HT__"+plotBaseName, File_preselection)
  #h_TTbar_amcatnlo_wjets = GetHisto(thisHistName.replace('SAMPLE',"TTbar_amcatnlo_Inc")+plotBaseName, File_preselection) # MC TTbar
  h_TTbar_powheg_wjets = GetHisto(thisHistName.replace('SAMPLE',"TTbar_powheg")+plotBaseName, File_preselection) # MC TTbar
  h_ZJets_amcatnlo_wjets = GetHisto(thisHistName.replace('SAMPLE',"ZJet_amcatnlo_ptBinned")+plotBaseName, File_preselection)
  h_WJets_amcatnlo_wjets = GetHisto(thisHistName.replace('SAMPLE',"WJet_amcatnlo_ptBinned")+plotBaseName, File_preselection)
  #h_WJets_amcatnlo_wjets = GetHisto("histo1D__WJet_amcatnlo_Inc__"+plotBaseName, File_preselection)
  h_SingleTop_wjets = GetHisto(thisHistName.replace('SAMPLE',"SingleTop")+plotBaseName, File_preselection)
  h_PhotonJets_wjets = GetHisto(thisHistName.replace('SAMPLE',"PhotonJets_Madgraph")+plotBaseName, File_preselection)
  #h_Diboson_wjets = GetHisto("histo1D__DIBOSON__"+plotBaseName, File_preselection)
  h_Diboson_wjets = GetHisto(thisHistName.replace('SAMPLE',"DIBOSON_amcatnlo")+plotBaseName, File_preselection)
  
  # DATA
  h_DATA_wjets = GetHisto(thisHistName.replace('SAMPLE',"DATA")+plotBaseName, File_preselection) #DATA
  # QCD
  h_QCD_DataDriven = GetHisto(thisHistName.replace('SAMPLE',"QCDFakes_DATA")+plotBaseName.replace('_btagSFDownShift','').replace('_btagSFUpShift',''),File_QCD_preselection)
  #h_QCD_wjets = GetHisto("histo1D__QCD_EMEnriched__"+plotBaseName,File_QCD_preselection)
  
  plotWJets = Plot()
  plotWJets.histoDATA = h_DATA_wjets
  if not useMGHT:
    # amc@NLO
    if useAMCAtNLOTTBar:
      plotWJets.histoMCall = h_ALLBKG_amcatnlo_wjets
      plotWJets.histoTTbar = h_TTbar_amcatnlo_wjets
    else:
      plotWJets.histoMCall = h_ALLBKG_powheg_wjets
      plotWJets.histoTTbar = h_TTbar_powheg_wjets
    #plotWJets.histoQCD = h_QCD_wjets
    plotWJets.histoQCD = h_QCD_DataDriven
    plotWJets.histoZJet = h_ZJets_amcatnlo_wjets
    plotWJets.histoWJet = h_WJets_amcatnlo_wjets
  else:
    # MG HT
    plotWJets.histoMCall = h_ALLBKG_HT_wjets
    plotWJets.histoTTbar = h_TTbar_MG_wjets
    plotWJets.histoQCD = h_QCD_wjets
    plotWJets.histoZJet = h_ZJets_MGHT_wjets
    plotWJets.histoWJet = h_WJets_MGHT_wjets
    
  plotWJets.histoSingleTop = h_SingleTop_wjets
  plotWJets.histoPhotonJets = h_PhotonJets_wjets
  plotWJets.histoDiboson = h_Diboson_wjets
  plotWJets.xmin = mtMin
  #plotWJets.xmax = h_TTbar_amcatnlo_wjets.GetXaxis().GetXmax()
  plotWJets.xmax = mtMax
  #plotWJets.name = "WJetsRescale"
  plotWJets.name = plotBaseName+'_WJets'
  plotWJets.fileXsectionNoRescale = "/afs/cern.ch/user/s/scooper/work/private/cmssw/8011/TestRootNTuplizerRecipe/src/Leptoquarks/analyzer/rootNtupleAnalyzerV2/config/xsection_13TeV_2015.txt"
  plotWJets.xminplot = 0
  plotWJets.xmaxplot = 2000
  plotWJets.yminplot = 0
  plotWJets.ymaxplot = 2000
  plotWJets.datasetName = "WJetsToLNu_Pt.+Tune"
  #plot0.datasetName = "DYJetsToLL_M-50_HT.+Tune"
  #plot0.datasetName = "Z.+Jets_Pt.+alpgen"
  # example: this match with /Z3Jets_Pt300to800-alpgen/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
  plotsWJets.append(plotWJets)

#-----------------------------------------------------------------------------------


############# USER CODE - END ################################################
##############################################################################



#--- Generate and print the plots from the list 'plots' define above

#--- Output files
fileps = "allPlots_calc_WJetsAndTTBarRescale_And_xsecFile.ps"

#--- Generate and print the plots from the list 'plots' define above
#c = TCanvas()
#c.Print(fileps+"[")
for idx,plot in enumerate(plotsTTBar):
  CalculateRescaleFactor(plot,plotsWJets[idx],fileps)
#c.Print(fileps+"]")
#os.system('ps2pdf '+fileps)

