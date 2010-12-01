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
    return integral

def GetErrorIntegralTH1( histo, xmin, xmax):
    #print "## calculating error for integral of histo " + str(histo)
    #print "## in the x range [" + str(xmin) + "," + str(xmax) + "]"
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
            #print "     --> skip bin: " + str(bin)
            pippo = 0
        else:
            error = error + histo.GetBinError(bin)**2
            #print "error**2 : " + str(error)

    error = sqrt(error)
    #print  " "
    return error



#--- Input files
File_Data_and_MC = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/33.2pb-1_sT_presel_250_Zrescale1.20_Wrescale1.06/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
File_QCD = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/6.1pb-1_QCD_HLT30_sT_presel_250/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")

#--- Histograms

##MTenu_PAS
#h_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MTenu_PAS", File_Data_and_MC) #DATA
#h_ALLBKG = GetHisto("histo1D__ALLBKG__cutHisto_allPreviousCuts________MTenu_PAS", File_Data_and_MC) # MC ALLBKG
#h_QCD = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MTenu_PAS", File_QCD, float(33.2/6.1)) #QCD (data-driven) scaled to the correct integrated lumi
##h_WJetAlpgen_MTenu = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection) # MC W
##h_ZJetAlpgen_MTenu = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection) # MC Z
##h_TTbar_Madgraph_MTenu = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection) # MC TTbar
##h_SingleTop_MTenu = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection) # MC SingleTop
##h_WW_MTenu = GetHisto("histo1D__WW__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection) # MC WW
##h_WZ_MTenu = GetHisto("histo1D__WZ__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection) # MC WZ
##h_ZZ_MTenu = GetHisto("histo1D__ZZ__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection) # MC ZZ

##Pt1stEle_PAS
h_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_Data_and_MC) #DATA
h_ALLBKG = GetHisto("histo1D__ALLBKG__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_Data_and_MC) # MC ALLBKG
h_QCD = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_QCD, float(33.2/6.1)) #QCD (data-driven) scaled to the correct integrated lumi


#--- Range
#variable = "MTenu"
variable = "Pt1stEle"
xmin = 250
xmax = 1000

#check
if(h_ALLBKG.GetNbinsX()!=h_DATA.GetNbinsX()):
    print "WARNING! number of bins is different between DATA and MC"
    print "exiting..."
    sys.exit()
if(h_ALLBKG.GetBinWidth(1)!=h_DATA.GetBinWidth(1)):
    print "WARNING! bin width is different between DATA and MC"
    print "exiting..."
    sys.exit()
    
#integrals
integralDATA = GetIntegralTH1(h_DATA,xmin,xmax)
ERRintegralDATA = GetErrorIntegralTH1(h_DATA,xmin,xmax)
integralMCall = GetIntegralTH1(h_ALLBKG,xmin,xmax)
ERRintegralMCall = GetErrorIntegralTH1(h_ALLBKG,xmin,xmax)
integralQCD = GetIntegralTH1(h_QCD,xmin,xmax)
ERRintegralQCD = GetErrorIntegralTH1(h_QCD,xmin,xmax)
#MC all + QCD
integralMCall_plus_QCD = integralMCall + integralQCD
ERRintegralMCall_plus_QCD = sqrt(ERRintegralMCall**2 + ERRintegralQCD**2)

##Difference between data and background prediction
diff_DATA_minus_Bkg = integralDATA - integralMCall_plus_QCD
ERRdiff_DATA_minus_Bkg = sqrt(ERRintegralDATA**2 + ERRintegralMCall_plus_QCD**2)

#printout
print " "
print "######################################## "
print "integral range: "          + str(xmin) + " < " + str(variable) + " < " + str(xmax) + " GeV"
print "integral DATA: "           + str( integralDATA ) + " +/- " + str( ERRintegralDATA )
print "integral MC All: "         + str( integralMCall ) + " +/- " + str( ERRintegralMCall )
print "integral QCD: "            + str( integralQCD ) + " +/- " + str( ERRintegralQCD )
print "integral MC All + QCD: "   + str( integralMCall_plus_QCD ) + " +/- " + str( ERRintegralMCall_plus_QCD )
print "difference DATA - TotBkg: "    + str( diff_DATA_minus_Bkg ) + " +/- " + str( ERRdiff_DATA_minus_Bkg )
print "######################################## "
print " "
    

