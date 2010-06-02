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


## The Plot class: add members if needed
class Plot:
    histos      = [] # list of histograms to be plotted in this plot
    keys        = [] # list of keys to be put in the legend (1 key per histo)
    histosStack = [] # list of histograms to be plotted in this plot -- stack histo format
    keysStack   = [] # list of keys to be put in the legend (1 key per histo) -- stack histo format
    xtit        = "" # xtitle
    ytit        = "" # ytitle 
    xmin        = "" # min x axis range (need to set both min and max. Leave it as is for full range)
    xmax        = "" # max x axis range (need to set both min and max. Leave it as is for full range)
    ymin        = "" # min y axis range (need to set both min and max. Leave it as is for full range)
    ymax        = "" # max y axis range (need to set both min and max. Leave it as is for full range)
    lpos        = "" # legend position (default = top-right, option="bottom-center")
    logscale    = "" # log scale of Y axis (default = no, option="yes")
    rebin       = "" # rebin x axis (default = 1, option = set it to whatever you want )
    name        = "" # name of the final plots
    
    def Draw(self, fileps):

        #-- create canvas
        canvas = TCanvas()
        stack = {}

        #-- log scale
        if (plot.logscale == "yes"):
            canvas.SetLogy();

        #-- legend
        if (plot.lpos=="bottom-center"):
            legend = TLegend(0.35, 0.25, 0.35+0.30, 0.25+0.20)
        else:
            legend = TLegend(0.60, 0.70, 0.89, 0.90)
        legend.SetFillColor(kWhite)
        legend.SetMargin(0.2)

        #-- loop over histograms (stacked)
        Nstacked = len(plot.histosStack)
        for iter in range(0, Nstacked):
            #make this stack
            stack[iter] = TH1F()
            Nloop = Nstacked - iter
            for iter1 in range(0,Nloop):            
                histo = plot.histosStack[iter1]
                if(iter1==0):
                    stack[iter].SetBins( histo.GetNbinsX(), histo.GetXaxis().GetXmin(), histo.GetXaxis().GetXmax() )
                    #stack[iter].SetName( plot.keysStack[iter] )
                stack[iter].Add(histo)
            #define style
            if(plot.rebin!=""):
                stack[iter].Rebin(plot.rebin)
            stack[iter].SetMarkerStyle(20+2*iter)
            stack[iter].SetMarkerColor(1+iter)
            stack[iter].SetLineColor(1+iter)
            stack[iter].SetFillColor(1+iter)
            legend.AddEntry(stack[iter], plot.keysStack[Nstacked - iter - 1],"lf")            
            #draw stack
            if iter==0:
                thisMin = stack[iter].GetXaxis().GetXmin()
                thisMax = stack[iter].GetXaxis().GetXmax()
                thisNbins = stack[iter].GetNbinsX()
                newBinning = (thisMax - thisMin) / thisNbins
                stack[iter].SetTitle("")
                stack[iter].GetXaxis().SetTitle(plot.xtit)
                stack[iter].GetYaxis().SetTitle(plot.ytit + " / ( "+ str(newBinning) + " )")
                if (plot.xmin!="" and plot.xmax!=""):
                    stack[iter].GetXaxis().SetRangeUser(plot.xmin,plot.xmax)
                if (plot.ymin!="" and plot.ymax!=""):
                    stack[iter].GetYaxis().SetLimits(plot.ymin,plot.ymax)
                    stack[iter].GetYaxis().SetRangeUser(plot.ymin,plot.ymax)
                #search for maximum of histograms
                #maxHisto = stack[iter].GetMaximum()
                #print maxHisto
                #for hh in plot.histos:
                #    if(plot.rebin!=""):
                #        if(hh.GetMaximum()*plot.rebin > maxHisto):
                #            maxHisto = hh.GetMaximum()*plot.rebin
                #    else:
                #        if(hh.GetMaximum() > maxHisto):
                #            maxHisto = hh.GetMaximum()
                #stack[iter].GetYaxis().SetLimits(0.,maxHisto*1.2)
                #stack[iter].GetYaxis().SetRangeUser(0.001,maxHisto*1.2)                    
                #draw first histo
                stack[iter].Draw("HIST")
            else:
                stack[iter].Draw("HISTsame")

        #-- loop over histograms (overlaid)
        ih=0 # index of histo within a plot
        for histo in plot.histos:
            if(plot.rebin!=""):
                histo.Rebin(plot.rebin)
            histo.SetMarkerStyle(20+2*ih)
            histo.SetMarkerColor(1+ih)
            histo.SetLineColor(1+ih)
            legend.AddEntry(histo, plot.keys[ih],"l")
            histo.Draw("HISTsame")
            ih=ih+1

        #-- draw label
        l = TLatex()
        l.SetTextAlign(12)
        l.SetTextSize(0.04)
        l.SetTextFont(62)
        l.SetNDC()
        if (plot.lpos=="bottom-center"):
            l.DrawLatex(0.35,0.20,"CMS 2010 Preliminary")
            l.DrawLatex(0.35,0.15,"L_{int} = 10 pb^{-1}")
        else:
            l.DrawLatex(0.60,0.65,"CMS 2010 Preliminary")
            l.DrawLatex(0.60,0.60,"L_{int} = 10 pb^{-1}")

        #-- end
        legend.Draw()
        canvas.Update()
        gPad.RedrawAxis()
        gPad.Modified()
        canvas.SaveAs(plot.name + ".eps")
        canvas.Print(fileps)



############# DON'T NEED TO MODIFY ANYTHING HERE - END #######################
##############################################################################


##############################################################################
############# USER CODE - BEGING #############################################

#--- Input root file

File1 = GetFile("/afs/cern.ch/user/s/santanas/scratch0/Releases/CMSSW_3_5_7_LQ_May19/src/Leptoquarks/rootNtupleAnalyzerV2/LQ_PAS_June2010/output/analysisClass_eejjSample_plots.root")

                
#--- Final plots are defined here

# Simply add or remove plots and update the list plots = [plot0, plot1, ...] below
    
#--- sT ---
h_sT_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allOtherCuts___________sT", File1)
h_sT_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allOtherCuts___________sT", File1)
h_sT_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allOtherCuts___________sT", File1)
h_sT_TTbar = GetHisto("histo1D__TTbar__cutHisto_allOtherCuts___________sT", File1)
h_sT_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allOtherCuts___________sT", File1)
h_sT_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allOtherCuts___________sT", File1)
h_sT_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allOtherCuts___________sT", File1)
h_sT_VVjets = GetHisto("histo1D__VVjets__cutHisto_allOtherCuts___________sT", File1)
h_sT_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allOtherCuts___________sT", File1)

plot0 = Plot() 
plot0.histosStack     = [h_sT_TTbar, h_sT_ZJetAlpgen, h_sT_QCD_Madgraph,
                         h_sT_SingleTop, h_sT_VVjets, h_sT_WJetAlpgen]
plot0.keysStack       = ["ttbar", "Z/#gamma/Z* + jets", "QCD multi-jets",
                         "single top", "di-bosons + jets", "W/W* + jets" ]
plot0.histos          = [h_sT_LQeejj_M300, h_sT_LQeejj_M400, h_sT_LQeejj_M500]
plot0.keys            = ["LQ eejj M300","LQ eejj M400","LQ eejj M500"]
plot0.xtit            = "sT (GeV)"
plot0.ytit            = "Number of events"
plot0.logscale        = "yes"
plot0.rebin           = 4
plot0.ymin            = 0.001
plot0.ymax            = 2
#plot0.lpos = "bottom-center"
plot0.name            = "sT_allOtherCuts"


#--- Mej ---
h_Mej_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allOtherCuts___________Mej", File1)
h_Mej_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allOtherCuts___________Mej", File1)
h_Mej_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allOtherCuts___________Mej", File1)
h_Mej_TTbar = GetHisto("histo1D__TTbar__cutHisto_allOtherCuts___________Mej", File1)
h_Mej_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allOtherCuts___________Mej", File1)
h_Mej_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allOtherCuts___________Mej", File1)
h_Mej_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allOtherCuts___________Mej", File1)
h_Mej_VVjets = GetHisto("histo1D__VVjets__cutHisto_allOtherCuts___________Mej", File1)
h_Mej_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allOtherCuts___________Mej", File1)

plot1 = Plot() 
plot1.histosStack     = [h_Mej_TTbar, h_Mej_ZJetAlpgen, h_Mej_QCD_Madgraph,
                         h_Mej_SingleTop, h_Mej_VVjets, h_Mej_WJetAlpgen]
plot1.keysStack       = ["ttbar", "Z/#gamma/Z* + jets", "QCD multi-jets",
                         "single top", "di-bosons + jets", "W/W* + jets" ]
plot1.histos          = [h_Mej_LQeejj_M300, h_Mej_LQeejj_M400, h_Mej_LQeejj_M500]
plot1.keys            = ["LQ eejj M300","LQ eejj M400","LQ eejj M500"]
plot1.xtit            = "Mej (GeV)"
plot1.ytit            = "Number of events"
plot1.logscale        = "yes"
plot1.rebin           = 4
plot1.ymin            = 0.001
plot1.ymax            = 3
plot1.name            = "Mej_allOtherCuts"


# List of plots to be plotted
plots = [plot0, plot1]


############# USER CODE - END ################################################
##############################################################################


#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()
c.Print("allPlots.ps[")
for plot in plots:
    plot.Draw("allPlots.ps")
c.Print("allPlots.ps]")
