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
File_QCD_fakeRate = GetFile("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/35.84pb-1_QCD_UseHLTPrescales_sT_presel_250/analysisClass_enujjSample_QCD_plots.root")
File_QCD_invertIso = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_QCD_InvertEleIso_sT_presel_250_Zrescale1.20_Wrescale1_fullntuples/output_cutTable_enujjSample_QCD_invertIso/analysisClass_enujjSample_QCD_invertIso_plots.root")

#--- Histograms

##MTenu_PAS
h_fakeRate_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MTenu_PAS", File_QCD_fakeRate) # DATA - QCD fake rate
h_invertIso_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MTenu_PAS", File_QCD_invertIso) # DATA - inverted ele isolation
h_invertIso_MC = GetHisto("histo1D__ALLBKG__cutHisto_allPreviousCuts________MTenu_PAS", File_QCD_invertIso) # MC - inverted ele isolation

#--- Range
variable = "MTenu"
xmin = 0
xmax = 100

#check
if(h_fakeRate_DATA.GetNbinsX()!=h_invertIso_DATA.GetNbinsX()):
    print "WARNING! number of bins is different between DATA and MC"
    print "exiting..."
    sys.exit()
if(h_fakeRate_DATA.GetBinWidth(1)!=h_invertIso_DATA.GetBinWidth(1)):
    print "WARNING! bin width is different between DATA and MC"
    print "exiting..."
    sys.exit()
    
#integrals
integral_fakeRate_DATA = GetIntegralTH1(h_fakeRate_DATA,xmin,xmax)
ERRintegral_fakeRate_DATA = GetErrorIntegralTH1(h_fakeRate_DATA,xmin,xmax)
integral_invertIso_DATA = GetIntegralTH1(h_invertIso_DATA,xmin,xmax)
ERRintegral_invertIso_DATA = GetErrorIntegralTH1(h_invertIso_DATA,xmin,xmax)
integral_invertIso_MC = GetIntegralTH1(h_invertIso_MC,xmin,xmax)
ERRintegral_invertIso_MC = GetErrorIntegralTH1(h_invertIso_MC,xmin,xmax)

scaleFactor = integral_fakeRate_DATA/integral_invertIso_DATA
ERRscaleFactor = scaleFactor * sqrt( (ERRintegral_fakeRate_DATA/integral_fakeRate_DATA)**2
                                     + (ERRintegral_invertIso_DATA/integral_invertIso_DATA)**2 )
h_invertIso_DATA.Scale(scaleFactor)

#plots
canvas = TCanvas()
canvas.SetLogy()

h_fakeRate_DATA.SetLineColor(3)
h_fakeRate_DATA.SetMarkerColor(3)
h_fakeRate_DATA.SetMarkerStyle(20)
h_fakeRate_DATA.Rebin(2)

h_invertIso_DATA.SetLineColor(4)
h_invertIso_DATA.SetMarkerColor(4)
h_invertIso_DATA.SetMarkerStyle(20)
h_invertIso_DATA.Rebin(2)

h_fakeRate_DATA.GetXaxis().SetRangeUser(0,600)

h_fakeRate_DATA.Draw();
h_invertIso_DATA.Draw("HISTEsames")

#legend
legend = TLegend( 0.36,0.63,0.74,0.85)
legend.AddEntry(h_fakeRate_DATA, "QCD c#nujj (fake rate)", "pl")
legend.AddEntry(h_invertIso_DATA, "QCD e#nujj (inverted iso)", "pl")
legend.Draw()

canvas.SaveAs("calc_QCDInvertedEleIso_Rescale.root")

#printout
print " "
print "######################################## "
print "integral range: "          + str(xmin) + " < " + str(variable) + " < " + str(xmax) + " GeV"
print "integral fakeRate_DATA: "           + str( integral_fakeRate_DATA ) + " +/- " + str( ERRintegral_fakeRate_DATA )
print "integral invertIso_DATA: "         + str( integral_invertIso_DATA ) + " +/- " + str( ERRintegral_invertIso_DATA )
print "integral invertIso_MC: "            + str( integral_invertIso_MC ) + " +/- " + str( ERRintegral_invertIso_MC )
print "integral scaleFactor: "            + str( scaleFactor ) + " +/- " + str( ERRscaleFactor )
print "######################################## "
print " "
    

