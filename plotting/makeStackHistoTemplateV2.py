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
    lpos        = "" # legend position (default = top-right, option="bottom-center", "top-left")
    logscale    = "" # log scale of Y axis (default = no, option="yes")
    rebin       = "" # rebin x axis (default = 1, option = set it to whatever you want )
    name        = "" # name of the final plots
    lint        = "11.81 nb^{-1}" # integrated luminosity of the sample ( example "10 pb^{-1}" )
    histodata   = "" # data histogram
    
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
        elif(plot.lpos=="top-left"):
            legend = TLegend(0.12, 0.70, 0.15+0.30, 0.70+0.20)            
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

        #-- plot data
        if(plot.rebin!=""):
            plot.histodata.Rebin(plot.rebin)
        plot.histodata.SetMarkerStyle(20)
        legend.AddEntry(plot.histodata, "data","p")
        plot.histodata.Draw("psame")

        #-- draw label
        l = TLatex()
        l.SetTextAlign(12)
        l.SetTextSize(0.04)
        l.SetTextFont(62)
        l.SetNDC()
        if (plot.lpos=="bottom-center"):
            l.DrawLatex(0.35,0.20,"CMS 2010 Preliminary")
            l.DrawLatex(0.35,0.15,"L_{int} = " + plot.lint)
        if (plot.lpos=="top-left"):
            l.DrawLatex(0.12,0.65,"CMS 2010 Preliminary")
            l.DrawLatex(0.12,0.60,"L_{int} = " + plot.lint)
        else:
            l.DrawLatex(0.60,0.65,"CMS 2010 Preliminary")
            l.DrawLatex(0.60,0.60,"L_{int} = " + plot.lint)

        #-- end
        legend.Draw()
        canvas.Update()
        gPad.RedrawAxis()
        gPad.Modified()
        canvas.SaveAs(plot.name + ".gif","gif")
        canvas.Print(fileps)



############# DON'T NEED TO MODIFY ANYTHING HERE - END #######################
##############################################################################


##############################################################################
############# USER CODE - BEGING #############################################

#--- Input root file

FileMC = GetFile("/afs/cern.ch/user/p/prumerio/scratch0/lq/collisions/data/output/analysisClass_eejjSample_plots.root")
#FileDATA = GetFile("/afs/cern.ch/user/p/prumerio/scratch0/lq/collisions/data/output/analysisClass_eejjSample___EG__Run2010A-PromptReco-v1__RECO.root")
#FileDATA = GetFile("/afs/cern.ch/user/p/prumerio/scratch0/lq/collisions/data/output/analysisClass_eejjSample_plots.root")
FileDATA = FileMC
               
#--- Final plots are defined here

# Simply add or remove plots and update the list plots = [plot0, plot1, ...] below
    
#--- Mee ---

h_Mee_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allPreviousCuts________Mee_PAS", FileMC)
h_Mee_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allPreviousCuts________Mee_PAS", FileMC)
h_Mee_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allPreviousCuts________Mee_PAS", FileMC)
h_Mee_TTbar = GetHisto("histo1D__TTbar__cutHisto_allPreviousCuts________Mee_PAS", FileMC)
h_Mee_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mee_PAS", FileMC)
#h_Mee_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Mee_PAS", FileMC)
h_Mee_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Mee_PAS", FileMC)
h_Mee_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Mee_PAS", FileMC)
h_Mee_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Mee_PAS", FileMC)
h_Mee_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Mee_PAS", FileMC)
h_Mee_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mee_PAS", FileDATA)


plot0 = Plot() 
## inputs for stacked histograms
## it created h_Mee_TTbar, h_Mee_TTbar+h_Mee_ZJetAlpgen , h_Mee_TTbar+h_Mee_ZJetAlpgen+h_Mee_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
plot0.histosStack     = [h_Mee_TTbar, h_Mee_ZJetAlpgen, h_Mee_QCDPt15,
                         h_Mee_SingleTop, h_Mee_VVjets, h_Mee_WJetAlpgen]
plot0.keysStack       = ["ttbar", "Z/#gamma/Z* + jets", "QCD multi-jets",
                         "single top", "di-bosons + jets", "W/W* + jets" ]
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot0.histos          = [h_Mee_LQeejj_M300, h_Mee_LQeejj_M400, h_Mee_LQeejj_M500]
plot0.keys            = ["LQ eejj M300","LQ eejj M400","LQ eejj M500"]
plot0.xtit            = "M(ee) (GeV)"
plot0.ytit            = "Number of events"
plot0.logscale        = "yes"
plot0.rebin           = 1
plot0.ymin            = 0.00000001
plot0.ymax            = 5
#plot0.lpos = "bottom-center"
plot0.name            = "Mee_allPreviousCuts"
plot0.histodata       = h_Mee_DATA

#--- nEle_PtCut_IDISO_noOvrlp ---

h_nEle_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", FileMC)
h_nEle_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", FileMC)
h_nEle_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", FileMC)
h_nEle_TTbar = GetHisto("histo1D__TTbar__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", FileMC)
h_nEle_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", FileMC)
#h_nEle_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", FileMC)
h_nEle_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", FileMC)
h_nEle_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", FileMC)
h_nEle_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", FileMC)
h_nEle_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", FileMC)
h_nEle_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", FileDATA)

plot1 = Plot() 
## inputs for stacked histograms
plot1.histosStack     = [h_nEle_TTbar, h_nEle_ZJetAlpgen, h_nEle_QCDPt15,
                         h_nEle_SingleTop, h_nEle_VVjets, h_nEle_WJetAlpgen]
plot1.keysStack       = ["ttbar", "Z/#gamma/Z* + jets", "QCD multi-jets",
                         "single top", "di-bosons + jets", "W/W* + jets" ]
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot1.histos          = [h_nEle_LQeejj_M300, h_nEle_LQeejj_M400, h_nEle_LQeejj_M500]
plot1.keys            = ["LQ eejj M300","LQ eejj M400","LQ eejj M500"]
plot1.xtit            = "Number of Heep electrons (pT>25 GeV)"
plot1.ytit            = "Number of events"
plot1.logscale        = "yes"
plot1.rebin           = 1
plot1.ymin            = 0.0001
plot1.ymax            = 1000000
#plot1.lpos = "bottom-center"
plot1.name            = "nEle_allPreviousCuts"
plot1.histodata       = h_nEle_DATA




#--- Pt1stEle_IDISO_NoOvrlp  ---

h_pT1stEle_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allPreviousCuts________Pt1stEle_IDISO_NoOvrlp", FileMC)
h_pT1stEle_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allPreviousCuts________Pt1stEle_IDISO_NoOvrlp", FileMC)
h_pT1stEle_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allPreviousCuts________Pt1stEle_IDISO_NoOvrlp", FileMC)
h_pT1stEle_TTbar = GetHisto("histo1D__TTbar__cutHisto_allPreviousCuts________Pt1stEle_IDISO_NoOvrlp", FileMC)
h_pT1stEle_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt1stEle_IDISO_NoOvrlp", FileMC)
#h_pT1stEle_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt1stEle_IDISO_NoOvrlp", FileMC)
h_pT1stEle_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt1stEle_IDISO_NoOvrlp", FileMC)
h_pT1stEle_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt1stEle_IDISO_NoOvrlp", FileMC)
h_pT1stEle_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt1stEle_IDISO_NoOvrlp", FileMC)
h_pT1stEle_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt1stEle_IDISO_NoOvrlp", FileMC)
h_pT1stEle_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt1stEle_IDISO_NoOvrlp", FileDATA)

plot2 = Plot() 
## inputs for stacked histograms
plot2.histosStack     = [h_pT1stEle_TTbar, h_pT1stEle_ZJetAlpgen, h_pT1stEle_QCDPt15,
                         h_pT1stEle_SingleTop, h_pT1stEle_VVjets, h_pT1stEle_WJetAlpgen]
plot2.keysStack       = ["ttbar", "Z/#gamma/Z* + jets", "QCD multi-jets",
                         "single top", "di-bosons + jets", "W/W* + jets" ]
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot2.histos          = [h_pT1stEle_LQeejj_M300, h_pT1stEle_LQeejj_M400, h_pT1stEle_LQeejj_M500]
plot2.keys            = ["LQ eejj M300","LQ eejj M400","LQ eejj M500"]
plot2.xtit            = "pT 1st HEEP electron (GeV)"
plot2.ytit            = "Number of events"
plot2.logscale        = "yes"
plot2.rebin           = 1
plot2.ymin            = 0.00000001
plot2.ymax            = 5
#plot2.lpos = "bottom-center"
plot2.name            = "pT1stEle_allPreviousCuts"
plot2.histodata       = h_pT1stEle_DATA


#--- Eta1stEle_IDISO_NoOvrlp  ---

h_Eta1stEle_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allPreviousCuts________Eta1stEle_IDISO_NoOvrlp", FileMC)
h_Eta1stEle_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allPreviousCuts________Eta1stEle_IDISO_NoOvrlp", FileMC)
h_Eta1stEle_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allPreviousCuts________Eta1stEle_IDISO_NoOvrlp", FileMC)
h_Eta1stEle_TTbar = GetHisto("histo1D__TTbar__cutHisto_allPreviousCuts________Eta1stEle_IDISO_NoOvrlp", FileMC)
h_Eta1stEle_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta1stEle_IDISO_NoOvrlp", FileMC)
#h_Eta1stEle_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta1stEle_IDISO_NoOvrlp", FileMC)
h_Eta1stEle_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta1stEle_IDISO_NoOvrlp", FileMC)
h_Eta1stEle_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta1stEle_IDISO_NoOvrlp", FileMC)
h_Eta1stEle_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta1stEle_IDISO_NoOvrlp", FileMC)
h_Eta1stEle_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta1stEle_IDISO_NoOvrlp", FileMC)
h_Eta1stEle_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta1stEle_IDISO_NoOvrlp", FileDATA)

plot3 = Plot() 
## inputs for stacked histograms
plot3.histosStack     = [h_Eta1stEle_TTbar, h_Eta1stEle_ZJetAlpgen, h_Eta1stEle_QCDPt15,
                         h_Eta1stEle_SingleTop, h_Eta1stEle_VVjets, h_Eta1stEle_WJetAlpgen]
plot3.keysStack       = ["ttbar", "Z/#gamma/Z* + jets", "QCD multi-jets",
                         "single top", "di-bosons + jets", "W/W* + jets" ]
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot3.histos          = [h_Eta1stEle_LQeejj_M300, h_Eta1stEle_LQeejj_M400, h_Eta1stEle_LQeejj_M500]
plot3.keys            = ["LQ eejj M300","LQ eejj M400","LQ eejj M500"]
plot3.xtit            = "#eta 1st HEEP electron"
plot3.ytit            = "Number of events"
plot3.logscale        = "yes"
plot3.rebin           = 10
plot3.ymin            = 0.0001
plot3.ymax            = 20
plot3.lpos = "top-left"
plot3.name            = "Eta1stEle_allPreviousCuts"
plot3.histodata       = h_Eta1stEle_DATA



#--- Pt2ndEle_IDISO_NoOvrlp  ---

h_pT2ndEle_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allPreviousCuts________Pt2ndEle_IDISO_NoOvrlp", FileMC)
h_pT2ndEle_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allPreviousCuts________Pt2ndEle_IDISO_NoOvrlp", FileMC)
h_pT2ndEle_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allPreviousCuts________Pt2ndEle_IDISO_NoOvrlp", FileMC)
h_pT2ndEle_TTbar = GetHisto("histo1D__TTbar__cutHisto_allPreviousCuts________Pt2ndEle_IDISO_NoOvrlp", FileMC)
h_pT2ndEle_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt2ndEle_IDISO_NoOvrlp", FileMC)
#h_pT2ndEle_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt2ndEle_IDISO_NoOvrlp", FileMC)
h_pT2ndEle_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt2ndEle_IDISO_NoOvrlp", FileMC)
h_pT2ndEle_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt2ndEle_IDISO_NoOvrlp", FileMC)
h_pT2ndEle_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt2ndEle_IDISO_NoOvrlp", FileMC)
h_pT2ndEle_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt2ndEle_IDISO_NoOvrlp", FileMC)
h_pT2ndEle_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt2ndEle_IDISO_NoOvrlp", FileDATA)

plot4 = Plot() 
## inputs for stacked histograms
plot4.histosStack     = [h_pT2ndEle_TTbar, h_pT2ndEle_ZJetAlpgen, h_pT2ndEle_QCDPt15,
                         h_pT2ndEle_SingleTop, h_pT2ndEle_VVjets, h_pT2ndEle_WJetAlpgen]
plot4.keysStack       = ["ttbar", "Z/#gamma/Z* + jets", "QCD multi-jets",
                         "single top", "di-bosons + jets", "W/W* + jets" ]
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot4.histos          = [h_pT2ndEle_LQeejj_M300, h_pT2ndEle_LQeejj_M400, h_pT2ndEle_LQeejj_M500]
plot4.keys            = ["LQ eejj M300","LQ eejj M400","LQ eejj M500"]
plot4.xtit            = "pT 2nd HEEP electron (GeV)"
plot4.ytit            = "Number of events"
plot4.logscale        = "yes"
plot4.rebin           = 1
plot4.ymin            = 0.00000001
plot4.ymax            = 5
#plot4.lpos = "bottom-center"
plot4.name            = "pT2ndEle_allPreviousCuts"
plot4.histodata       = h_pT2ndEle_DATA


#--- Eta2ndEle_IDISO_NoOvrlp  ---

h_Eta2ndEle_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allPreviousCuts________Eta2ndEle_IDISO_NoOvrlp", FileMC)
h_Eta2ndEle_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allPreviousCuts________Eta2ndEle_IDISO_NoOvrlp", FileMC)
h_Eta2ndEle_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allPreviousCuts________Eta2ndEle_IDISO_NoOvrlp", FileMC)
h_Eta2ndEle_TTbar = GetHisto("histo1D__TTbar__cutHisto_allPreviousCuts________Eta2ndEle_IDISO_NoOvrlp", FileMC)
h_Eta2ndEle_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta2ndEle_IDISO_NoOvrlp", FileMC)
#h_Eta2ndEle_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta2ndEle_IDISO_NoOvrlp", FileMC)
h_Eta2ndEle_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta2ndEle_IDISO_NoOvrlp", FileMC)
h_Eta2ndEle_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta2ndEle_IDISO_NoOvrlp", FileMC)
h_Eta2ndEle_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta2ndEle_IDISO_NoOvrlp", FileMC)
h_Eta2ndEle_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta2ndEle_IDISO_NoOvrlp", FileMC)
h_Eta2ndEle_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta2ndEle_IDISO_NoOvrlp", FileDATA)

plot5 = Plot() 
## inputs for stacked histograms
plot5.histosStack     = [h_Eta2ndEle_TTbar, h_Eta2ndEle_ZJetAlpgen, h_Eta2ndEle_QCDPt15,
                         h_Eta2ndEle_SingleTop, h_Eta2ndEle_VVjets, h_Eta2ndEle_WJetAlpgen]
plot5.keysStack       = ["ttbar", "Z/#gamma/Z* + jets", "QCD multi-jets",
                         "single top", "di-bosons + jets", "W/W* + jets" ]
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot5.histos          = [h_Eta2ndEle_LQeejj_M300, h_Eta2ndEle_LQeejj_M400, h_Eta2ndEle_LQeejj_M500]
plot5.keys            = ["LQ eejj M300","LQ eejj M400","LQ eejj M500"]
plot5.xtit            = "#eta 2nd HEEP electron"
plot5.ytit            = "Number of events"
plot5.logscale        = "yes"
plot5.rebin           = 10
plot5.ymin            = 0.0001
plot5.ymax            = 20
plot5.lpos = "top-left"
plot5.name            = "Eta2ndEle_allPreviousCuts"
plot5.histodata       = h_Eta2ndEle_DATA


#--- nJet_PAS ---

h_nJet_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allPreviousCuts________nJet_PAS", FileMC)
h_nJet_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allPreviousCuts________nJet_PAS", FileMC)
h_nJet_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allPreviousCuts________nJet_PAS", FileMC)
h_nJet_TTbar = GetHisto("histo1D__TTbar__cutHisto_allPreviousCuts________nJet_PAS", FileMC)
h_nJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________nJet_PAS", FileMC)
#h_nJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________nJet_PAS", FileMC)
h_nJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________nJet_PAS", FileMC)
h_nJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________nJet_PAS", FileMC)
h_nJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________nJet_PAS", FileMC)
h_nJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________nJet_PAS", FileMC)
h_nJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________nJet_PAS", FileDATA)


plot6 = Plot() 
## inputs for stacked histograms
plot6.histosStack     = [h_nJet_TTbar, h_nJet_ZJetAlpgen, h_nJet_QCDPt15,
                         h_nJet_SingleTop, h_nJet_VVjets, h_nJet_WJetAlpgen]
plot6.keysStack       = ["ttbar", "Z/#gamma/Z* + jets", "QCD multi-jets",
                         "single top", "di-bosons + jets", "W/W* + jets" ]
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot6.histos          = [h_nJet_LQeejj_M300, h_nJet_LQeejj_M400, h_nJet_LQeejj_M500]
plot6.keys            = ["LQ eejj M300","LQ eejj M400","LQ eejj M500"]
plot6.xtit            = "Number of jets (pT>10GeV)"
plot6.ytit            = "Number of events"
plot6.logscale        = "yes"
plot6.rebin           = 1
plot6.ymin            = 0.0001
plot6.ymax            = 5
#plot6.lpos = "bottom-center"
plot6.name            = "nJet_allPreviousCuts"
plot6.histodata       = h_nJet_DATA



#--- Pt1stJet_PAS ---

h_Pt1stJet_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allPreviousCuts________Pt1stJet_PAS", FileMC)
h_Pt1stJet_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allPreviousCuts________Pt1stJet_PAS", FileMC)
h_Pt1stJet_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allPreviousCuts________Pt1stJet_PAS", FileMC)
h_Pt1stJet_TTbar = GetHisto("histo1D__TTbar__cutHisto_allPreviousCuts________Pt1stJet_PAS", FileMC)
h_Pt1stJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt1stJet_PAS", FileMC)
#h_Pt1stJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt1stJet_PAS", FileMC)
h_Pt1stJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt1stJet_PAS", FileMC)
h_Pt1stJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt1stJet_PAS", FileMC)
h_Pt1stJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt1stJet_PAS", FileMC)
h_Pt1stJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt1stJet_PAS", FileMC)
h_Pt1stJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt1stJet_PAS", FileDATA)


plot7 = Plot() 
## inputs for stacked histograms
plot7.histosStack     = [h_Pt1stJet_TTbar, h_Pt1stJet_ZJetAlpgen, h_Pt1stJet_QCDPt15,
                         h_Pt1stJet_SingleTop, h_Pt1stJet_VVjets, h_Pt1stJet_WJetAlpgen]
plot7.keysStack       = ["ttbar", "Z/#gamma/Z* + jets", "QCD multi-jets",
                         "single top", "di-bosons + jets", "W/W* + jets" ]
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot7.histos          = [h_Pt1stJet_LQeejj_M300, h_Pt1stJet_LQeejj_M400, h_Pt1stJet_LQeejj_M500]
plot7.keys            = ["LQ eejj M300","LQ eejj M400","LQ eejj M500"]
plot7.xtit            = "pT 1st jet (GeV)"
plot7.ytit            = "Number of events"
plot7.logscale        = "yes"
plot7.rebin           = 1
plot7.ymin            = 0.00000001
plot7.ymax            = 5
#plot7.lpos = "bottom-center"
plot7.name            = "Pt1stJet_allPreviousCuts"
plot7.histodata       = h_Pt1stJet_DATA


#--- Eta1stJet_PAS ---

h_Eta1stJet_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allPreviousCuts________Eta1stJet_PAS", FileMC)
h_Eta1stJet_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allPreviousCuts________Eta1stJet_PAS", FileMC)
h_Eta1stJet_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allPreviousCuts________Eta1stJet_PAS", FileMC)
h_Eta1stJet_TTbar = GetHisto("histo1D__TTbar__cutHisto_allPreviousCuts________Eta1stJet_PAS", FileMC)
h_Eta1stJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta1stJet_PAS", FileMC)
#h_Eta1stJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta1stJet_PAS", FileMC)
h_Eta1stJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta1stJet_PAS", FileMC)
h_Eta1stJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta1stJet_PAS", FileMC)
h_Eta1stJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta1stJet_PAS", FileMC)
h_Eta1stJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta1stJet_PAS", FileMC)
h_Eta1stJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta1stJet_PAS", FileDATA)


plot8 = Plot() 
## inputs for stacked histograms
plot8.histosStack     = [h_Eta1stJet_TTbar, h_Eta1stJet_ZJetAlpgen, h_Eta1stJet_QCDPt15,
                         h_Eta1stJet_SingleTop, h_Eta1stJet_VVjets, h_Eta1stJet_WJetAlpgen]
plot8.keysStack       = ["ttbar", "Z/#gamma/Z* + jets", "QCD multi-jets",
                         "single top", "di-bosons + jets", "W/W* + jets" ]
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot8.histos          = [h_Eta1stJet_LQeejj_M300, h_Eta1stJet_LQeejj_M400, h_Eta1stJet_LQeejj_M500]
plot8.keys            = ["LQ eejj M300","LQ eejj M400","LQ eejj M500"]
plot8.xtit            = "#eta 1st jet (GeV)"
plot8.ytit            = "Number of events"
plot8.logscale        = "yes"
plot8.rebin           = 20
plot8.ymin            = 0.0001
plot8.ymax            = 20
plot8.lpos = "top-left"
plot8.name            = "Eta1stJet_allPreviousCuts"
plot8.histodata       = h_Eta1stJet_DATA


#--- Pt2ndJet_PAS ---

h_Pt2ndJet_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allPreviousCuts________Pt2ndJet_PAS", FileMC)
h_Pt2ndJet_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allPreviousCuts________Pt2ndJet_PAS", FileMC)
h_Pt2ndJet_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allPreviousCuts________Pt2ndJet_PAS", FileMC)
h_Pt2ndJet_TTbar = GetHisto("histo1D__TTbar__cutHisto_allPreviousCuts________Pt2ndJet_PAS", FileMC)
h_Pt2ndJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt2ndJet_PAS", FileMC)
#h_Pt2ndJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt2ndJet_PAS", FileMC)
h_Pt2ndJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt2ndJet_PAS", FileMC)
h_Pt2ndJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt2ndJet_PAS", FileMC)
h_Pt2ndJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt2ndJet_PAS", FileMC)
h_Pt2ndJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt2ndJet_PAS", FileMC)
h_Pt2ndJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt2ndJet_PAS", FileDATA)


plot9 = Plot() 
## inputs for stacked histograms
plot9.histosStack     = [h_Pt2ndJet_TTbar, h_Pt2ndJet_ZJetAlpgen, h_Pt2ndJet_QCDPt15,
                         h_Pt2ndJet_SingleTop, h_Pt2ndJet_VVjets, h_Pt2ndJet_WJetAlpgen]
plot9.keysStack       = ["ttbar", "Z/#gamma/Z* + jets", "QCD multi-jets",
                         "single top", "di-bosons + jets", "W/W* + jets" ]
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot9.histos          = [h_Pt2ndJet_LQeejj_M300, h_Pt2ndJet_LQeejj_M400, h_Pt2ndJet_LQeejj_M500]
plot9.keys            = ["LQ eejj M300","LQ eejj M400","LQ eejj M500"]
plot9.xtit            = "pT 2nd jet (GeV)"
plot9.ytit            = "Number of events"
plot9.logscale        = "yes"
plot9.rebin           = 1
plot9.ymin            = 0.00000001
plot9.ymax            = 5
#plot9.lpos = "bottom-center"
plot9.name            = "Pt2ndJet_allPreviousCuts"
plot9.histodata       = h_Pt2ndJet_DATA


#--- Eta2ndJet_PAS ---

h_Eta2ndJet_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allPreviousCuts________Eta2ndJet_PAS", FileMC)
h_Eta2ndJet_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allPreviousCuts________Eta2ndJet_PAS", FileMC)
h_Eta2ndJet_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allPreviousCuts________Eta2ndJet_PAS", FileMC)
h_Eta2ndJet_TTbar = GetHisto("histo1D__TTbar__cutHisto_allPreviousCuts________Eta2ndJet_PAS", FileMC)
h_Eta2ndJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta2ndJet_PAS", FileMC)
#h_Eta2ndJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta2ndJet_PAS", FileMC)
h_Eta2ndJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta2ndJet_PAS", FileMC)
h_Eta2ndJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta2ndJet_PAS", FileMC)
h_Eta2ndJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta2ndJet_PAS", FileMC)
h_Eta2ndJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta2ndJet_PAS", FileMC)
h_Eta2ndJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta2ndJet_PAS", FileDATA)


plot10 = Plot() 
## inputs for stacked histograms
plot10.histosStack     = [h_Eta2ndJet_TTbar, h_Eta2ndJet_ZJetAlpgen, h_Eta2ndJet_QCDPt15,
                         h_Eta2ndJet_SingleTop, h_Eta2ndJet_VVjets, h_Eta2ndJet_WJetAlpgen]
plot10.keysStack       = ["ttbar", "Z/#gamma/Z* + jets", "QCD multi-jets",
                         "single top", "di-bosons + jets", "W/W* + jets" ]
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot10.histos          = [h_Eta2ndJet_LQeejj_M300, h_Eta2ndJet_LQeejj_M400, h_Eta2ndJet_LQeejj_M500]
plot10.keys            = ["LQ eejj M300","LQ eejj M400","LQ eejj M500"]
plot10.xtit            = "#eta 2nd jet (GeV)"
plot10.ytit            = "Number of events"
plot10.logscale        = "yes"
plot10.rebin           = 10
plot10.ymin            = 0.0001
plot10.ymax            = 20
plot10.lpos = "top-left"
plot10.name            = "Eta2ndJet_allPreviousCuts"
plot10.histodata       = h_Eta2ndJet_DATA

#--- sT ---

h_sT_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allPreviousCuts________sT_PAS", FileMC)
h_sT_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allPreviousCuts________sT_PAS", FileMC)
h_sT_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allPreviousCuts________sT_PAS", FileMC)
h_sT_TTbar = GetHisto("histo1D__TTbar__cutHisto_allPreviousCuts________sT_PAS", FileMC)
h_sT_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________sT_PAS", FileMC)
#h_sT_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________sT_PAS", FileMC)
h_sT_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________sT_PAS", FileMC)
h_sT_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________sT_PAS", FileMC)
h_sT_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________sT_PAS", FileMC)
h_sT_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________sT_PAS", FileMC)
h_sT_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________sT_PAS", FileDATA)


plot11 = Plot() 
## inputs for stacked histograms
## it created h_sT_TTbar, h_sT_TTbar+h_sT_ZJetAlpgen , h_sT_TTbar+h_sT_ZJetAlpgen+h_sT_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
plot11.histosStack     = [h_sT_TTbar, h_sT_ZJetAlpgen, h_sT_QCDPt15,
                         h_sT_SingleTop, h_sT_VVjets, h_sT_WJetAlpgen]
plot11.keysStack       = ["ttbar", "Z/#gamma/Z* + jets", "QCD multi-jets",
                         "single top", "di-bosons + jets", "W/W* + jets" ]
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot11.histos          = [h_sT_LQeejj_M300, h_sT_LQeejj_M400, h_sT_LQeejj_M500]
plot11.keys            = ["LQ eejj M300","LQ eejj M400","LQ eejj M500"]
plot11.xtit            = "St (GeV)"
plot11.ytit            = "Number of events"
plot11.logscale        = "yes"
plot11.rebin           = 1
plot11.ymin            = 0.00000001
plot11.ymax            = 5
#plot11.lpos = "bottom-center"
plot11.name            = "sT_allPreviousCuts"
plot11.histodata       = h_sT_DATA




# List of plots to be plotted
plots = [plot0, plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9, plot10, plot11]


############# USER CODE - END ################################################
##############################################################################


#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()
c.Print("allPlots.ps[")
for plot in plots:
    plot.Draw("allPlots.ps")
c.Print("allPlots.ps]")
