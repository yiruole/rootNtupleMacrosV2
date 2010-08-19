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
import ROOT
from array import array

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


def GetHisto( histoName , file ):
    histo = file.Get( histoName )
    if( not histo):
        print "ERROR: histo " + histoName + " not found in " + file.GetName()
        print "exiting..."
        sys.exit()
    return histo

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
	print "bin: " +str(bin)
        if(bin==bmax and bmaxResidual==histo.GetBinContent(bmax)): # skip last bin if out of range
            print "     --> skip bin: " + str(bin)
        else:           
            error = error + histo.GetBinError(bin)**2    
            print "error**2 : " + str(error)

    error = sqrt(error)
    print  " "
    return error


## The Plot class: add members if needed
class Plot:
    histoDATA   = "" # DATA
    histoMCZ    = "" # MCZ
    histoMCall  = "" # MCall
    xtit        = "" # xtitle
    ytit        = "" # ytitle 
    xmin        = "" # set xmin to calculate rescaling factor (-- please take into account the bin size of histograms --)
    xmax        = "" # # set xmax to calculate rescaling factor (-- please take into account the bin size of histograms --)
    xminplot    = "" # min x axis range (need to set both min and max. Leave it as is for full range)
    xmaxplot    = "" # max x axis range (need to set both min and max. Leave it as is for full range)
    yminplot    = "" # min y axis range (need to set both min and max. Leave it as is for full range)
    ymaxplot    = "" # max y axis range (need to set both min and max. Leave it as is for full range)
    lpos        = "" # legend position (default = top-right, option="bottom-center", "top-left")
    #    xlog        = "" # log scale of X axis (default = no, option="yes") ### IT SEEMS IT DOES NOT WORK
    ylog        = "" # log scale of Y axis (default = no, option="yes")
    #rebin      = "" # rebin x axis (default = 1, option = set it to whatever you want )
    name        = "" # name of the final plots
    lint        = "828 nb^{-1}" # integrated luminosity of the sample ( example "10 pb^{-1}" )

    def CalculateRescaleFactor(self, fileps):
        #calculate rescaling factor for Z/gamma+jet background and create new cross section file
        canvas = TCanvas()

        #check
        if(plot.histoMCall.GetNbinsX()!=plot.histoDATA.GetNbinsX()):
            print "WARNING! number of bins is different between DATA and MC"
            print "existing..."
            sys.exit()
        if(plot.histoMCall.GetBinWidth(1)!=plot.histoDATA.GetBinWidth(1)):
            print "WARNING! bin width is different between DATA and MC"
            print "existing..."
            sys.exit()

        #integrals
        integralMCall = GetIntegralTH1(plot.histoMCall,plot.xmin,plot.xmax)
        ERRintegralMCall = GetErrorIntegralTH1(plot.histoMCall,plot.xmin,plot.xmax)
        integralDATA = GetIntegralTH1(plot.histoDATA,plot.xmin,plot.xmax)
        ERRintegralDATA = GetErrorIntegralTH1(plot.histoDATA,plot.xmin,plot.xmax)
        integralMCZ = GetIntegralTH1(plot.histoMCZ,plot.xmin,plot.xmax)
        ERRintegralMCZ = GetErrorIntegralTH1(plot.histoMCZ,plot.xmin,plot.xmax)

        #contamination
        contamination = (integralMCall - integralMCZ) / integralMCall

        #rescale factor
        rescale = integralMCall / integralDATA
        relERRintegralMCall = ERRintegralMCall / integralMCall
        relERRintegralDATA = ERRintegralDATA / integralDATA
        ERRrescale = sqrt(relERRintegralMCall**2 + relERRintegralDATA**2)


        #draw histo
        plot.histoMCall.SetFillColor(kBlue)
        plot.histoDATA.SetMarkerStyle(20)
        
        plot.histoMCall.Draw("HIST")
        plot.histoDATA.Draw("psame")
        plot.histoMCall.GetXaxis().SetRangeUser(plot.xminplot,plot.xmaxplot)
        plot.histoMCall.GetYaxis().SetRangeUser(plot.yminplot,plot.ymaxplot)
        
        canvas.Update()
        gPad.RedrawAxis()
        gPad.Modified()
        #canvas.SaveAs(plot.name + ".eps","eps")
        #canvas.SaveAs(plot.name + ".pdf","pdf")
        canvas.Print(fileps)

        #printout
        print " "
        print "######################################## "
        print "integralMCall: "   + str( integralMCall ) + " +/- " + str( ERRintegralMCall )
        print "integralDATA: "   + str( integralDATA ) + " +/- " + str( ERRintegralDATA )
        print "contribution from all other backgrounds (except Z+jet): " + str(contamination*100) + "%"
        print "rescale factor for Z background: " + str(rescale) + " +\- " + str(ERRrescale*rescale)
        print "systematical uncertainty of Z+jet background modeling: " + str(ERRrescale*100) + "%"
        print "######################################## "
        print " "

############# DON'T NEED TO MODIFY ANYTHING HERE - END #######################
##############################################################################


##############################################################################
############# USER CODE - BEGING #############################################

#--- Input files
#preselection
File_preselection = GetFile("$LQDATA/eejj_analysis/eejj_825nb-1_preSelJet20GeV_noDeltaEta/output_cutTable_eejjSample/analysisClass_eejjSample_plots.root")

#xsection file
File_xsection_noRescale = "$LQANA/config/xsection_7TeV.txt"


#--- Rescaling of Z/gamma + jet background

#-----------------------------------------
h_ALLBKG_Mee = GetHisto("histo1D__ALLBKG__cutHisto_allPreviousCuts________Mee", File_preselection).Clone() # MC all
h_ZJetAlpgen_Mee = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mee", File_preselection).Clone() # MC Z
h_DATA_Mee = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mee", File_preselection).Clone() #DATA

plot0 = Plot()
plot0.histoDATA = h_DATA_Mee
plot0.histoMCall = h_ALLBKG_Mee
plot0.histoMCZ = h_ZJetAlpgen_Mee
plot0.xmin = 80  
plot0.xmax = 100 
plot0.name = "h_Mee_compare"
plot0.xminplot = 0
plot0.xmaxplot = 200
plot0.yminplot = 0
plot0.ymaxplot = 15
#datasetName = "Z0Jets_Pt0to100-alpgen"# string for pattern recognition of dataset name
#                                      # (rescaling will be applied only to those datasets)

plots = [plot0]

#-----------------------------------------------------------------------------------


############# USER CODE - END ################################################
##############################################################################



#--- Generate and print the plots from the list 'plots' define above

#--- Output files
fileps = "allPlots_calc_MCrescale_AND_xsecFile.ps"

#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()
c.Print(fileps+"[")
for plot in plots:
    plot.CalculateRescaleFactor(fileps)
c.Print(fileps+"]")
os.system('ps2pdf '+fileps)

