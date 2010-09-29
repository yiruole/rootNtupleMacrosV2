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


def GetHisto( histoName , file , scale = 1 ):
    file.cd()
    histo = file.Get( histoName )
    if(scale!=1):
        histo.Scale(scale)
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
    #    xlog        = "" # log scale of X axis (default = no, option="yes") ### IT SEEMS IT DOES NOT WORK
    ylog        = "" # log scale of Y axis (default = no, option="yes")
    rebin       = "" # rebin x axis (default = 1, option = set it to whatever you want )
    name        = "" # name of the final plots
    lint        = "2.9 pb^{-1}" # integrated luminosity of the sample ( example "10 pb^{-1}" )
    addZUncBand = "no" # add an uncertainty band coming from the data-MC Z+jets rescaling (default = "no", option="yes")
    ZUncKey     = "Z/#gamma/Z* + jets unc." # key to be put in the legend for the Z+jets uncertainty band
    ZPlotIndex  = 1 # index of the Z+jets plots in the histosStack list (default = 1)
    ZScaleUnc   = 0.20 # uncertainty of the data-MC Z+jets scale factor
    histodata   = "" # data histogram

    def Draw(self, fileps):

        #-- create canvas
        canvas = TCanvas()
        stack = {}

        #-- log scale
#             xlog may npot work         if (plot.xlog     == "yes"):
#             canvas.SetLogx();
        if (plot.ylog     == "yes"):
            canvas.SetLogy();

        #-- legend
        hsize=0.20
        vsize=0.25
        if (plot.lpos=="bottom-center"):
            xstart=0.35
            ystart=0.25
        elif(plot.lpos=="top-left"):
            xstart=0.12
            ystart=0.63
        else:
            xstart=0.68
            ystart=0.63
        legend = TLegend(xstart, ystart, xstart+hsize, ystart+vsize)
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
            stack[iter].SetMarkerColor(15+10*iter)
            stack[iter].SetLineColor(  15+10*iter)
            stack[iter].SetFillColor(  15+10*iter)
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

        #-- Z+jets uncertainty band
        if(plot.addZUncBand == "yes"):
            Zhisto = plot.histosStack[plot.ZPlotIndex].Clone()
            if(plot.rebin!=""):
                Zhisto.Rebin(plot.rebin)
            zUncHisto = stack[0].Clone()
            for bin in range(0,Zhisto.GetNbinsX()):
              zUncHisto.SetBinError(bin+1,plot.ZScaleUnc*Zhisto.GetBinContent(bin+1))
            zUncHisto.SetMarkerStyle(0)
            zUncHisto.SetLineColor(0)
            zUncHisto.SetFillColor(5)
            zUncHisto.SetFillStyle(3154)
            zUncHisto.Draw("E2same")
            legend.AddEntry(zUncHisto, plot.ZUncKey,"f")

        #-- loop over histograms (overlaid)
        ih=0 # index of histo within a plot
        for histo in plot.histos:
            if(plot.rebin!=""):
                histo.Rebin(plot.rebin)
            histo.SetMarkerStyle(20+2*ih)
            histo.SetMarkerColor(2+2*ih)
            histo.SetLineColor(  2+2*ih)
            legend.AddEntry(histo, plot.keys[ih],"l")
            histo.Draw("HISTsame")
            ih=ih+1

        #-- plot data
        if(plot.histodata!=""):
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
#        l.DrawLatex(xstart,ystart-0.05,"CMS Preliminary 2010")
#        l.DrawLatex(xstart,ystart-0.10,"L_{int} = " + plot.lint)
        if (plot.lpos=="bottom-center"):
            l.DrawLatex(0.35,0.20,"CMS Preliminary 2010")
            l.DrawLatex(0.35,0.15,"L_{int} = " + plot.lint)
        if (plot.lpos=="top-left"):
            l.DrawLatex(xstart+hsize+0.02,ystart+vsize-0.03,"CMS Preliminary 2010")
            l.DrawLatex(xstart+hsize+0.02,ystart+vsize-0.08,"L_{int} = " + plot.lint)
        else:
            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.03,"CMS Preliminary 2010")
            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.08,"L_{int} = " + plot.lint)

        #-- end
        legend.Draw()
        canvas.Update()
        gPad.RedrawAxis()
        gPad.Modified()
        canvas.SaveAs(plot.name + ".eps","eps")
        canvas.SaveAs(plot.name + ".pdf","pdf")
        canvas.Print(fileps)



############# DON'T NEED TO MODIFY ANYTHING HERE - END #######################
##############################################################################


##############################################################################
############# USER CODE - BEGING #############################################

#--- Input root file


#usually preselection and selection file are the same actual file (but we keep the name separate anyway)
File_preselection = GetFile("$LQDATA/enujj_analysis/2.9pb-1_v7/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
File_selection    = GetFile("$LQDATA/enujj_analysis/2.9pb-1_v7/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")

UseQCDFromData    = 1 #set to zero if you don't use QCD from data
#always put an existing file under File_QCD (otherwise the code will crash)
File_QCD          = GetFile("$LQDATA/enujj_analysis/2.9pb-1_v7_QCD_HLT/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")
QCDscaleFactor    = 1.26 # ratio between integrated lumi of the signal sample (i.e. 2.9 pb-1) / integrated lumi of the QCD sample (i.e. 2.3 pb-1 from HLT Photon20)

                  
#### Common values for plots:
#otherBkgsKey="single top, VV+jets, Z/Z*/gamma+jets, QCD?"
otherBkgsKey="Other Bkgs"
zUncBand="no"

pt_xmin=0
pt_xmax=800
pt_ymin=0.001
pt_ymax=500

eta_rebin=4
eta_ymin=0
eta_ymax=60


#--- Final plots are defined here

# Simply add or remove plots and update the list plots = [plot0, plot1, ...] below


#--- nEle_PtCut_IDISO_noOvrlp ---

h_nEle_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
h_nEle_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
h_nEle_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
h_nEle_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
h_nEle_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
h_nEle_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
h_nEle_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
h_nEle_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
#h_nEle_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
#h_nEle_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
h_nEle_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
h_nEle_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
h_nEle_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
h_nEle_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
if(UseQCDFromData):
    h_nEle_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_QCD, QCDscaleFactor).Clone()

plot0 = Plot()
## inputs for stacked histograms
if(UseQCDFromData):
    plot0.histosStack     = [h_nEle_QCDFromData, h_nEle_TTbar, h_nEle_WJetAlpgen, h_nEle_OTHERBKG]
    plot0.keysStack       = ["QCD", "ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot0.histosStack     = [h_nEle_TTbar, h_nEle_WJetAlpgen, h_nEle_OTHERBKG]
    plot0.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot0.histos          = [h_nEle_LQenujj_M200, h_nEle_LQenujj_M300]
plot0.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot0.xtit            = "Number of electrons"
plot0.ytit            = "Number of events"
plot0.ylog            = "yes"
plot0.rebin           = 1
plot0.ymin            = 0.0001
plot0.ymax            = 100000000
#plot0.lpos = "bottom-center"
plot0.name            = "nEle_allPreviousCuts"
plot0.addZUncBand     = zUncBand
plot0.histodata       = h_nEle_DATA


#--- Pt1stEle_PAS  ---

h_pT1stEle_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_preselection).Clone()
h_pT1stEle_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_preselection).Clone()
h_pT1stEle_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_preselection).Clone()
h_pT1stEle_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_preselection).Clone()
h_pT1stEle_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_preselection).Clone()
h_pT1stEle_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_preselection).Clone()
h_pT1stEle_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_preselection).Clone()
h_pT1stEle_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_preselection).Clone()
#h_pT1stEle_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_preselection).Clone()
#h_pT1stEle_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_preselection).Clone()
h_pT1stEle_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_preselection).Clone()
h_pT1stEle_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_preselection).Clone()
h_pT1stEle_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_preselection).Clone()
h_pT1stEle_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_pT1stEle_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt1stEle_PAS", File_QCD, QCDscaleFactor).Clone()

plot1 = Plot()
## inputs for stacked histograms
if(UseQCDFromData):
    plot1.histosStack     = [h_pT1stEle_QCDFromData, h_pT1stEle_TTbar, h_pT1stEle_WJetAlpgen, h_pT1stEle_OTHERBKG]
    plot1.keysStack       = ["QCD", "ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot1.histosStack     = [h_pT1stEle_TTbar, h_pT1stEle_WJetAlpgen, h_pT1stEle_OTHERBKG]
    plot1.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot1.histos          = [h_pT1stEle_LQenujj_M200, h_pT1stEle_LQenujj_M300]
plot1.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot1.xtit            = "pT 1st electron (GeV/c)"
plot1.ytit            = "Number of events"
plot1.ylog            = "yes"
plot1.rebin           = 1
plot1.xmin            = pt_xmin
plot1.xmax            = pt_xmax
plot1.ymin            = pt_ymin
plot1.ymax            = pt_ymax
#plot1.lpos = "bottom-center"
plot1.name            = "pT1stEle_allPreviousCuts"
plot1.addZUncBand     = zUncBand
plot1.histodata       = h_pT1stEle_DATA


#--- Eta1stEle_PAS  ---

h_Eta1stEle_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________Eta1stEle_PAS", File_preselection).Clone()
h_Eta1stEle_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________Eta1stEle_PAS", File_preselection).Clone()
h_Eta1stEle_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________Eta1stEle_PAS", File_preselection).Clone()
h_Eta1stEle_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________Eta1stEle_PAS", File_preselection).Clone()
h_Eta1stEle_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________Eta1stEle_PAS", File_preselection).Clone()
h_Eta1stEle_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Eta1stEle_PAS", File_preselection).Clone()
h_Eta1stEle_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta1stEle_PAS", File_preselection).Clone()
h_Eta1stEle_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Eta1stEle_PAS", File_preselection).Clone()
#h_Eta1stEle_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta1stEle_PAS", File_preselection).Clone()
#h_Eta1stEle_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta1stEle_PAS", File_preselection).Clone()
h_Eta1stEle_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta1stEle_PAS", File_preselection).Clone()
h_Eta1stEle_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta1stEle_PAS", File_preselection).Clone()
h_Eta1stEle_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta1stEle_PAS", File_preselection).Clone()
h_Eta1stEle_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta1stEle_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_Eta1stEle_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta1stEle_PAS", File_QCD, QCDscaleFactor).Clone()

plot2 = Plot()
## inputs for stacked histograms
if(UseQCDFromData):
    plot2.histosStack     = [h_Eta1stEle_QCDFromData, h_Eta1stEle_TTbar, h_Eta1stEle_WJetAlpgen, h_Eta1stEle_OTHERBKG]
    plot2.keysStack       = ["QCD","ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot2.histosStack     = [h_Eta1stEle_TTbar, h_Eta1stEle_WJetAlpgen, h_Eta1stEle_OTHERBKG]
    plot2.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot2.histos          = [h_Eta1stEle_LQenujj_M200, h_Eta1stEle_LQenujj_M300]
plot2.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot2.xtit            = "#eta 1st electron"
plot2.ytit            = "Number of events"
plot2.rebin           = eta_rebin
plot2.ymin            = eta_ymin
plot2.ymax            = eta_ymax
plot2.lpos = "top-left"
plot2.name            = "Eta1stEle_allPreviousCuts"
plot2.addZUncBand     = zUncBand
plot2.histodata       = h_Eta1stEle_DATA


#--- MET ---

h_MET_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________MET_PAS", File_preselection).Clone()
h_MET_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________MET_PAS", File_preselection).Clone()
h_MET_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________MET_PAS", File_preselection).Clone()
h_MET_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________MET_PAS", File_preselection).Clone()
h_MET_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________MET_PAS", File_preselection).Clone()
h_MET_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________MET_PAS", File_preselection).Clone()
h_MET_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________MET_PAS", File_preselection).Clone()
h_MET_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________MET_PAS", File_preselection).Clone()
#h_MET_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________MET_PAS", File_preselection).Clone()
#h_MET_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________MET_PAS", File_preselection).Clone()
h_MET_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________MET_PAS", File_preselection).Clone()
h_MET_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________MET_PAS", File_preselection).Clone()
h_MET_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________MET_PAS", File_preselection).Clone()
h_MET_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MET_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_MET_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MET_PAS", File_QCD, QCDscaleFactor).Clone()

plot3 = Plot()
## inputs for stacked histograms
## it created h_MET_TTbar, h_MET_TTbar+h_MET_ZJetAlpgen , h_MET_TTbar+h_MET_ZJetAlpgen+h_MET_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
if(UseQCDFromData):
    plot3.histosStack     = [h_MET_QCDFromData, h_MET_TTbar, h_MET_WJetAlpgen, h_MET_OTHERBKG]
    plot3.keysStack       = ["QCD", "ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot3.histosStack     = [h_MET_TTbar, h_MET_WJetAlpgen, h_MET_OTHERBKG]
    plot3.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot3.histos          = [h_MET_LQenujj_M200, h_MET_LQenujj_M300, h_MET_LQenujj_M400, h_MET_LQenujj_M500]
plot3.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300", "LQ e\\nujj M400", "LQ e\\nujj M500"]
plot3.xtit            = "pfMET (GeV/c)"
plot3.ytit            = "Number of events"
#plot3.xlog            = "yes"
plot3.ylog            = "yes"
plot3.rebin           = 1
plot3.xmin            = 0
plot3.xmax            = 800
plot3.ymin            = 0.001
plot3.ymax            = 500
#plot3.lpos = "bottom-center"
plot3.name            = "MET_allPreviousCuts"
plot3.addZUncBand     = zUncBand
plot3.histodata       = h_MET_DATA


#--- minMETPt1stEle_PAS ---

h_minMETPt1stEle_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________minMETPt1stEle_PAS", File_preselection).Clone()
h_minMETPt1stEle_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________minMETPt1stEle_PAS", File_preselection).Clone()
h_minMETPt1stEle_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________minMETPt1stEle_PAS", File_preselection).Clone()
h_minMETPt1stEle_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________minMETPt1stEle_PAS", File_preselection).Clone()
h_minMETPt1stEle_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________minMETPt1stEle_PAS", File_preselection).Clone()
h_minMETPt1stEle_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________minMETPt1stEle_PAS", File_preselection).Clone()
h_minMETPt1stEle_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________minMETPt1stEle_PAS", File_preselection).Clone()
h_minMETPt1stEle_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________minMETPt1stEle_PAS", File_preselection).Clone()
#h_minMETPt1stEle_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________minMETPt1stEle_PAS", File_preselection).Clone()
#h_minMETPt1stEle_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________minMETPt1stEle_PAS", File_preselection).Clone()
h_minMETPt1stEle_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________minMETPt1stEle_PAS", File_preselection).Clone()
h_minMETPt1stEle_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________minMETPt1stEle_PAS", File_preselection).Clone()
h_minMETPt1stEle_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________minMETPt1stEle_PAS", File_preselection).Clone()
h_minMETPt1stEle_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________minMETPt1stEle_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_minMETPt1stEle_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________minMETPt1stEle_PAS", File_QCD, QCDscaleFactor).Clone()

plot4 = Plot()
## inputs for stacked histograms
## it created h_minMETPt1stEle_TTbar, h_minMETPt1stEle_TTbar+h_minMETPt1stEle_ZJetAlpgen , h_minMETPt1stEle_TTbar+h_minMETPt1stEle_ZJetAlpgen+h_minMETPt1stEle_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
if(UseQCDFromData):
    plot4.histosStack     = [h_minMETPt1stEle_QCDFromData, h_minMETPt1stEle_TTbar, h_minMETPt1stEle_WJetAlpgen, h_minMETPt1stEle_OTHERBKG]
    plot4.keysStack       = ["QCD", "ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot4.histosStack     = [h_minMETPt1stEle_TTbar, h_minMETPt1stEle_WJetAlpgen, h_minMETPt1stEle_OTHERBKG]
    plot4.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot4.histos          = [h_minMETPt1stEle_LQenujj_M200, h_minMETPt1stEle_LQenujj_M300, h_minMETPt1stEle_LQenujj_M400, h_minMETPt1stEle_LQenujj_M500]
plot4.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300", "LQ e\\nujj M400", "LQ e\\nujj M500"]
plot4.xtit            = "min(pT 1st electron,pfMET) (GeV/c)"
plot4.ytit            = "Number of events"
#plot4.xlog            = "yes"
plot4.ylog            = "yes"
plot4.rebin           = 1
plot4.xmin            = 0
plot4.xmax            = 500
plot4.ymin            = 0.001
plot4.ymax            = 500
#plot4.lpos = "bottom-center"
plot4.name            = "minMETPt1stEle_allPreviousCuts"
plot4.addZUncBand     = zUncBand
plot4.histodata       = h_minMETPt1stEle_DATA


#--- nJet_WithJetEtaCut ---

h_nJet_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________nJet_WithJetEtaCut", File_preselection).Clone()
h_nJet_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________nJet_WithJetEtaCut", File_preselection).Clone()
h_nJet_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________nJet_WithJetEtaCut", File_preselection).Clone()
h_nJet_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________nJet_WithJetEtaCut", File_preselection).Clone()
h_nJet_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________nJet_WithJetEtaCut", File_preselection).Clone()
h_nJet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________nJet_WithJetEtaCut", File_preselection).Clone()
h_nJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________nJet_WithJetEtaCut", File_preselection).Clone()
h_nJet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________nJet_WithJetEtaCut", File_preselection).Clone()
#h_nJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________nJet_WithJetEtaCut", File_preselection).Clone()
#h_nJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________nJet_WithJetEtaCut", File_preselection).Clone()
h_nJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________nJet_WithJetEtaCut", File_preselection).Clone()
h_nJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________nJet_WithJetEtaCut", File_preselection).Clone()
h_nJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________nJet_WithJetEtaCut", File_preselection).Clone()
h_nJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________nJet_WithJetEtaCut", File_preselection).Clone()
if(UseQCDFromData):
    h_nJet_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________nJet_WithJetEtaCut", File_QCD, QCDscaleFactor).Clone()

plot5 = Plot()
## inputs for stacked histograms
if(UseQCDFromData):
    plot5.histosStack     = [h_nJet_QCDFromData, h_nJet_TTbar, h_nJet_WJetAlpgen, h_nJet_OTHERBKG]
    plot5.keysStack       = ["QCD", "ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot5.histosStack     = [h_nJet_TTbar, h_nJet_WJetAlpgen, h_nJet_OTHERBKG]
    plot5.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot5.histos          = [h_nJet_LQenujj_M200, h_nJet_LQenujj_M300]
plot5.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot5.xtit            = "Number of jets"
plot5.ytit            = "Number of events"
plot5.ylog            = "yes"
plot5.rebin           = 1
plot5.xmin            = 0
plot5.xmax            = 12
plot5.ymin            = 0.01
plot5.ymax            = 100000000
#plot5.lpos = "bottom-center"
plot5.name            = "nJet_allPreviousCuts"
plot5.addZUncBand     = zUncBand
plot5.histodata       = h_nJet_DATA



#--- Pt1stJet_PAS ---

h_Pt1stJet_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_Pt1stJet_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_Pt1stJet_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_Pt1stJet_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_Pt1stJet_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_Pt1stJet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_Pt1stJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_Pt1stJet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
#h_Pt1stJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
#h_Pt1stJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_Pt1stJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_Pt1stJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_Pt1stJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_Pt1stJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_Pt1stJet_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_QCD, QCDscaleFactor).Clone()

plot6 = Plot()
## inputs for stacked histograms
if(UseQCDFromData):
    plot6.histosStack     = [h_Pt1stJet_QCDFromData,h_Pt1stJet_TTbar, h_Pt1stJet_WJetAlpgen, h_Pt1stJet_OTHERBKG]
    plot6.keysStack       = ["QCD","ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot6.histosStack     = [h_Pt1stJet_TTbar, h_Pt1stJet_WJetAlpgen, h_Pt1stJet_OTHERBKG]
    plot6.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot6.histos          = [h_Pt1stJet_LQenujj_M200, h_Pt1stJet_LQenujj_M300]
plot6.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot6.xtit            = "pT 1st jet (GeV/c)"
plot6.ytit            = "Number of events"
plot6.ylog            = "yes"
plot6.rebin           = 1
plot6.xmin            = pt_xmin
plot6.xmax            = pt_xmax
plot6.ymin            = pt_ymin
plot6.ymax            = pt_ymax
#plot6.lpos = "bottom-center"
plot6.name            = "Pt1stJet_allPreviousCuts"
plot6.addZUncBand     = zUncBand
plot6.histodata       = h_Pt1stJet_DATA


#--- Eta1stJet_PAS ---

h_Eta1stJet_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_Eta1stJet_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_Eta1stJet_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_Eta1stJet_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_Eta1stJet_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_Eta1stJet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_Eta1stJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_Eta1stJet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
#h_Eta1stJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
#h_Eta1stJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_Eta1stJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_Eta1stJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_Eta1stJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_Eta1stJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_Eta1stJet_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_QCD, QCDscaleFactor).Clone()

plot7 = Plot()
## inputs for stacked histograms
if(UseQCDFromData):
    plot7.histosStack     = [h_Eta1stJet_QCDFromData, h_Eta1stJet_TTbar, h_Eta1stJet_WJetAlpgen, h_Eta1stJet_OTHERBKG]
    plot7.keysStack       = ["QCD", "ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot7.histosStack     = [h_Eta1stJet_TTbar, h_Eta1stJet_WJetAlpgen, h_Eta1stJet_OTHERBKG]
    plot7.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot7.histos          = [h_Eta1stJet_LQenujj_M200, h_Eta1stJet_LQenujj_M300]
plot7.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot7.xtit            = "#eta 1st jet"
plot7.ytit            = "Number of events"
plot7.rebin           = eta_rebin
plot7.ymin            = eta_ymin
plot7.ymax            = eta_ymax
plot7.lpos = "top-left"
plot7.name            = "Eta1stJet_allPreviousCuts"
plot7.addZUncBand     = zUncBand
plot7.histodata       = h_Eta1stJet_DATA


#--- Pt2ndJet_PAS ---

h_Pt2ndJet_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
h_Pt2ndJet_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
h_Pt2ndJet_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
h_Pt2ndJet_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
h_Pt2ndJet_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
h_Pt2ndJet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
h_Pt2ndJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
h_Pt2ndJet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
#h_Pt2ndJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
#h_Pt2ndJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
h_Pt2ndJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
h_Pt2ndJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
h_Pt2ndJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
h_Pt2ndJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_Pt2ndJet_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_QCD, QCDscaleFactor).Clone()

plot8 = Plot()
## inputs for stacked histograms
if(UseQCDFromData):
    plot8.histosStack     = [h_Pt2ndJet_QCDFromData, h_Pt2ndJet_TTbar, h_Pt2ndJet_WJetAlpgen, h_Pt2ndJet_OTHERBKG]
    plot8.keysStack       = ["QCD","ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot8.histosStack     = [h_Pt2ndJet_TTbar, h_Pt2ndJet_WJetAlpgen, h_Pt2ndJet_OTHERBKG]
    plot8.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot8.histos          = [h_Pt2ndJet_LQenujj_M200, h_Pt2ndJet_LQenujj_M300]
plot8.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot8.xtit            = "pT 2nd jet (GeV/c)"
plot8.ytit            = "Number of events"
plot8.ylog            = "yes"
plot8.xmin            = pt_xmin
plot8.xmax            = pt_xmax
plot8.ymin            = pt_ymin
plot8.ymax            = pt_ymax
#plot8.lpos = "bottom-center"
plot8.name            = "Pt2ndJet_allPreviousCuts"
plot8.addZUncBand     = zUncBand
plot8.histodata       = h_Pt2ndJet_DATA
    
#--- Eta2ndJet_PAS ---

h_Eta2ndJet_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
h_Eta2ndJet_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
h_Eta2ndJet_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
h_Eta2ndJet_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
h_Eta2ndJet_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
h_Eta2ndJet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
h_Eta2ndJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
h_Eta2ndJet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
#h_Eta2ndJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
#h_Eta2ndJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
h_Eta2ndJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
h_Eta2ndJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
h_Eta2ndJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
h_Eta2ndJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_Eta2ndJet_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_QCD, QCDscaleFactor).Clone()

plot9 = Plot()
## inputs for stacked histograms
if(UseQCDFromData):
    plot9.histosStack     = [h_Eta2ndJet_QCDFromData,h_Eta2ndJet_TTbar, h_Eta2ndJet_WJetAlpgen, h_Eta2ndJet_OTHERBKG]
    plot9.keysStack       = ["QCD","ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot9.histosStack     = [h_Eta2ndJet_TTbar, h_Eta2ndJet_WJetAlpgen, h_Eta2ndJet_OTHERBKG]
    plot9.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot9.histos          = [h_Eta2ndJet_LQenujj_M200, h_Eta2ndJet_LQenujj_M300]
plot9.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot9.xtit            = "#eta 2nd jet"
plot9.ytit            = "Number of events"
plot9.rebin           = eta_rebin
plot9.ymin            = eta_ymin
plot9.ymax            = eta_ymax
plot9.lpos = "top-left"
plot9.name            = "Eta2ndJet_allPreviousCuts"
plot9.addZUncBand     = zUncBand
plot9.histodata       = h_Eta2ndJet_DATA


##--- Pt Jets AllPreviousCuts ---

h_pTJets_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_pTJets_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_pTJets_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_pTJets_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_pTJets_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_pTJets_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_pTJets_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_pTJets_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
#h_pTJets_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
#h_pTJets_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_pTJets_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_pTJets_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_pTJets_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_pTJets_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_pTJets_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_QCD, QCDscaleFactor).Clone()

h_pTJets_LQenujj_M100.Add(GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
h_pTJets_LQenujj_M200.Add(GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
h_pTJets_LQenujj_M300.Add(GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
h_pTJets_LQenujj_M400.Add(GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
h_pTJets_LQenujj_M500.Add(GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
h_pTJets_TTbar.Add(GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
h_pTJets_ZJetAlpgen.Add(GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
h_pTJets_OTHERBKG.Add(GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
#h_pTJets_QCD_Madgraph.Add(GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
#h_pTJets_QCDPt15.Add(GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
h_pTJets_SingleTop.Add(GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
h_pTJets_VVjets.Add(GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
h_pTJets_WJetAlpgen.Add(GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
h_pTJets_DATA.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
if(UseQCDFromData):
    h_pTJets_QCDFromData.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_QCD, QCDscaleFactor))

plot6and8 = Plot()
## inputs for stacked histograms
if(UseQCDFromData):
    plot6and8.histosStack     = [h_pTJets_QCDFromData, h_pTJets_TTbar, h_pTJets_WJetAlpgen, h_pTJets_OTHERBKG]
    plot6and8.keysStack       = ["QCD","ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot6and8.histosStack     = [h_pTJets_TTbar, h_pTJets_WJetAlpgen, h_pTJets_OTHERBKG]
    plot6and8.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot6and8.histos          = [h_pTJets_LQenujj_M200, h_pTJets_LQenujj_M300]
plot6and8.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot6and8.xtit            = "pT jets (GeV/c)"
plot6and8.ytit            = "Number of events x 2"
plot6and8.ylog            = "yes"
plot6and8.rebin           = 1
plot6and8.xmin            = pt_xmin
plot6and8.xmax            = pt_xmax
plot6and8.ymin            = pt_ymin
plot6and8.ymax            = pt_ymax
#plot6and8.lpos = "bottom-center"
plot6and8.name            = "pTJets_allPreviousCuts"
plot6and8.addZUncBand     = zUncBand
plot6and8.histodata       = h_pTJets_DATA

##--- Eta Eles AllPreviousCuts ---

h_etaJets_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_etaJets_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_etaJets_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_etaJets_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_etaJets_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_etaJets_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_etaJets_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_etaJets_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
#h_etaJets_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
#h_etaJets_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_etaJets_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_etaJets_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_etaJets_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_etaJets_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_etaJets_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_QCD, QCDscaleFactor).Clone()

h_etaJets_LQenujj_M100.Add(GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
h_etaJets_LQenujj_M200.Add(GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
h_etaJets_LQenujj_M300.Add(GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
h_etaJets_LQenujj_M400.Add(GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
h_etaJets_LQenujj_M500.Add(GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
h_etaJets_TTbar.Add(GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
h_etaJets_ZJetAlpgen.Add(GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
h_etaJets_OTHERBKG.Add(GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
#h_etaJets_QCD_Madgraph.Add(GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
#h_etaJets_QCDPt15.Add(GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
h_etaJets_SingleTop.Add(GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
h_etaJets_VVjets.Add(GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
h_etaJets_WJetAlpgen.Add(GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
h_etaJets_DATA.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
if(UseQCDFromData):
    h_etaJets_QCDFromData.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_QCD, QCDscaleFactor))

plot7and9 = Plot()
## inputs for stacked histograms
if(UseQCDFromData):
    plot7and9.histosStack     = [h_etaJets_QCDFromData,h_etaJets_TTbar, h_etaJets_WJetAlpgen, h_etaJets_OTHERBKG]
    plot7and9.keysStack       = ["QCD","ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot7and9.histosStack     = [h_etaJets_TTbar, h_etaJets_WJetAlpgen, h_etaJets_OTHERBKG]
    plot7and9.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot7and9.histos          = [h_etaJets_LQenujj_M200, h_etaJets_LQenujj_M300]
plot7and9.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot7and9.xtit            = "#eta jets"
plot7and9.ytit            = "Number of events x 2"
plot7and9.rebin           = eta_rebin/2
plot7and9.ymin            = eta_ymin
plot7and9.ymax            = eta_ymax*1.3
plot7and9.lpos            = "top-left"
#plot7and9.lpos = "bottom-center"
plot7and9.name            = "etaJets_allPreviousCuts"
plot7and9.addZUncBand     = zUncBand
plot7and9.histodata       = h_etaJets_DATA


#--- nMuon_PtCut_IDISO_PAS ---

h_nMuon_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________nMuon_PtCut_IDISO_PAS", File_preselection).Clone()
h_nMuon_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________nMuon_PtCut_IDISO_PAS", File_preselection).Clone()
h_nMuon_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________nMuon_PtCut_IDISO_PAS", File_preselection).Clone()
h_nMuon_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________nMuon_PtCut_IDISO_PAS", File_preselection).Clone()
h_nMuon_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________nMuon_PtCut_IDISO_PAS", File_preselection).Clone()
h_nMuon_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________nMuon_PtCut_IDISO_PAS", File_preselection).Clone()
h_nMuon_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________nMuon_PtCut_IDISO_PAS", File_preselection).Clone()
h_nMuon_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________nMuon_PtCut_IDISO_PAS", File_preselection).Clone()
#h_nMuon_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________nMuon_PtCut_IDISO_PAS", File_preselection).Clone()
#h_nMuon_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________nMuon_PtCut_IDISO_PAS", File_preselection).Clone()
h_nMuon_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________nMuon_PtCut_IDISO_PAS", File_preselection).Clone()
h_nMuon_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________nMuon_PtCut_IDISO_PAS", File_preselection).Clone()
h_nMuon_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________nMuon_PtCut_IDISO_PAS", File_preselection).Clone()
h_nMuon_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________nMuon_PtCut_IDISO_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_nMuon_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________nMuon_PtCut_IDISO_PAS", File_QCD, QCDscaleFactor).Clone()

plot10 = Plot()
## inputs for stacked histograms
if(UseQCDFromData):
    plot10.histosStack     = [h_nMuon_QCDFromData,h_nMuon_TTbar, h_nMuon_WJetAlpgen, h_nMuon_OTHERBKG]
    plot10.keysStack       = ["QCD","ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot10.histosStack     = [h_nMuon_TTbar, h_nMuon_WJetAlpgen, h_nMuon_OTHERBKG]
    plot10.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot10.histos          = [h_nMuon_LQenujj_M200, h_nMuon_LQenujj_M300]
plot10.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot10.xtit            = "Number of muons"
plot10.ytit            = "Number of events"
plot10.ylog            = "yes"
plot10.rebin           = 1
plot10.xmin            = 0
plot10.xmax            = 12
plot10.ymin            = 0.01
plot10.ymax            = 10000
#plot10.lpos = "bottom-center"
plot10.name            = "nMuon_allPreviousCuts"
plot10.addZUncBand     = zUncBand
plot10.histodata       = h_nMuon_DATA


#--- mDeltaPhiMETEle_PAS ---

h_mDeltaPhiMETEle_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________mDeltaPhiMETEle_PAS", File_preselection).Clone()
h_mDeltaPhiMETEle_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________mDeltaPhiMETEle_PAS", File_preselection).Clone()
h_mDeltaPhiMETEle_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________mDeltaPhiMETEle_PAS", File_preselection).Clone()
h_mDeltaPhiMETEle_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________mDeltaPhiMETEle_PAS", File_preselection).Clone()
h_mDeltaPhiMETEle_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________mDeltaPhiMETEle_PAS", File_preselection).Clone()
h_mDeltaPhiMETEle_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________mDeltaPhiMETEle_PAS", File_preselection).Clone()
h_mDeltaPhiMETEle_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________mDeltaPhiMETEle_PAS", File_preselection).Clone()
h_mDeltaPhiMETEle_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________mDeltaPhiMETEle_PAS", File_preselection).Clone()
#h_mDeltaPhiMETEle_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________mDeltaPhiMETEle_PAS", File_preselection).Clone()
#h_mDeltaPhiMETEle_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________mDeltaPhiMETEle_PAS", File_preselection).Clone()
h_mDeltaPhiMETEle_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________mDeltaPhiMETEle_PAS", File_preselection).Clone()
h_mDeltaPhiMETEle_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________mDeltaPhiMETEle_PAS", File_preselection).Clone()
h_mDeltaPhiMETEle_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________mDeltaPhiMETEle_PAS", File_preselection).Clone()
h_mDeltaPhiMETEle_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________mDeltaPhiMETEle_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_mDeltaPhiMETEle_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________mDeltaPhiMETEle_PAS", File_QCD, QCDscaleFactor).Clone()

plot11 = Plot()
## inputs for stacked histograms
## it created h_sT_TTbar, h_sT_TTbar+h_sT_ZJetAlpgen , h_sT_TTbar+h_sT_ZJetAlpgen+h_sT_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
if(UseQCDFromData):
    plot11.histosStack     = [h_mDeltaPhiMETEle_QCDFromData, h_mDeltaPhiMETEle_TTbar, h_mDeltaPhiMETEle_WJetAlpgen, h_mDeltaPhiMETEle_OTHERBKG]
    plot11.keysStack       = ["QCD","ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot11.histosStack     = [h_mDeltaPhiMETEle_TTbar, h_mDeltaPhiMETEle_WJetAlpgen, h_mDeltaPhiMETEle_OTHERBKG]
    plot11.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot11.histos          = [h_mDeltaPhiMETEle_LQenujj_M200, h_mDeltaPhiMETEle_LQenujj_M300]
plot11.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot11.xtit            = "#Delta#phi(MET,e) (rad.)"
plot11.ytit            = "Number of events"
#plot11.xlog            = "yes"
plot11.ylog            = "yes"
plot11.rebin           = 1
plot11.xmin            = 0
plot11.xmax            = 5
plot11.ymin            = 0.001
plot11.ymax            = 10000
#plot11.lpos = "bottom-center"
plot11.name            = "mDeltaPhiMETEle_allPreviousCuts"
plot11.addZUncBand     = zUncBand
plot11.histodata       = h_mDeltaPhiMETEle_DATA


#--- mDeltaPhiMET1stJet_PAS ---

h_mDeltaPhiMET1stJet_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________mDeltaPhiMET1stJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET1stJet_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________mDeltaPhiMET1stJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET1stJet_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________mDeltaPhiMET1stJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET1stJet_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________mDeltaPhiMET1stJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET1stJet_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________mDeltaPhiMET1stJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET1stJet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________mDeltaPhiMET1stJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET1stJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________mDeltaPhiMET1stJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET1stJet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________mDeltaPhiMET1stJet_PAS", File_preselection).Clone()
#h_mDeltaPhiMET1stJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________mDeltaPhiMET1stJet_PAS", File_preselection).Clone()
#h_mDeltaPhiMET1stJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________mDeltaPhiMET1stJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET1stJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________mDeltaPhiMET1stJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET1stJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________mDeltaPhiMET1stJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET1stJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________mDeltaPhiMET1stJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET1stJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________mDeltaPhiMET1stJet_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_mDeltaPhiMET1stJet_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________mDeltaPhiMET1stJet_PAS", File_QCD, QCDscaleFactor).Clone()

plot12 = Plot()
## inputs for stacked histograms
## it created h_sT_TTbar, h_sT_TTbar+h_sT_ZJetAlpgen , h_sT_TTbar+h_sT_ZJetAlpgen+h_sT_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
if(UseQCDFromData):
    plot12.histosStack     = [h_mDeltaPhiMET1stJet_QCDFromData,h_mDeltaPhiMET1stJet_TTbar, h_mDeltaPhiMET1stJet_WJetAlpgen, h_mDeltaPhiMET1stJet_OTHERBKG]
    plot12.keysStack       = ["QCD","ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot12.histosStack     = [h_mDeltaPhiMET1stJet_TTbar, h_mDeltaPhiMET1stJet_WJetAlpgen, h_mDeltaPhiMET1stJet_OTHERBKG]
    plot12.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot12.histos          = [h_mDeltaPhiMET1stJet_LQenujj_M200, h_mDeltaPhiMET1stJet_LQenujj_M300]
plot12.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot12.xtit            = "#Delta#phi(MET,1^{st} jet) (rad.)"
plot12.ytit            = "Number of events"
#plot12.xlog            = "yes"
plot12.ylog            = "yes"
plot12.rebin           = 1
plot12.xmin            = 0
plot12.xmax            = 5
plot12.ymin            = 0.001
plot12.ymax            = 10000
#plot12.lpos = "bottom-center"
plot12.name            = "mDeltaPhiMET1stJet_allPreviousCuts"
plot12.addZUncBand     = zUncBand
plot12.histodata       = h_mDeltaPhiMET1stJet_DATA


#--- mDeltaPhiMET2ndJet_PAS ---

h_mDeltaPhiMET2ndJet_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________mDeltaPhiMET2ndJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET2ndJet_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________mDeltaPhiMET2ndJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET2ndJet_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________mDeltaPhiMET2ndJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET2ndJet_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________mDeltaPhiMET2ndJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET2ndJet_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________mDeltaPhiMET2ndJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET2ndJet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________mDeltaPhiMET2ndJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET2ndJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________mDeltaPhiMET2ndJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET2ndJet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________mDeltaPhiMET2ndJet_PAS", File_preselection).Clone()
#h_mDeltaPhiMET2ndJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________mDeltaPhiMET2ndJet_PAS", File_preselection).Clone()
#h_mDeltaPhiMET2ndJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________mDeltaPhiMET2ndJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET2ndJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________mDeltaPhiMET2ndJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET2ndJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________mDeltaPhiMET2ndJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET2ndJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________mDeltaPhiMET2ndJet_PAS", File_preselection).Clone()
h_mDeltaPhiMET2ndJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________mDeltaPhiMET2ndJet_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_mDeltaPhiMET2ndJet_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________mDeltaPhiMET2ndJet_PAS", File_QCD, QCDscaleFactor).Clone()

plot13 = Plot()
## inputs for stacked histograms
## it created h_sT_TTbar, h_sT_TTbar+h_sT_ZJetAlpgen , h_sT_TTbar+h_sT_ZJetAlpgen+h_sT_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
if(UseQCDFromData):
    plot13.histosStack     = [h_mDeltaPhiMET2ndJet_QCDFromData, h_mDeltaPhiMET2ndJet_TTbar, h_mDeltaPhiMET2ndJet_WJetAlpgen, h_mDeltaPhiMET2ndJet_OTHERBKG]
    plot13.keysStack       = ["QCD","ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot13.histosStack     = [h_mDeltaPhiMET2ndJet_TTbar, h_mDeltaPhiMET2ndJet_WJetAlpgen, h_mDeltaPhiMET2ndJet_OTHERBKG]
    plot13.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot13.histos          = [h_mDeltaPhiMET2ndJet_LQenujj_M200, h_mDeltaPhiMET2ndJet_LQenujj_M300]
plot13.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot13.xtit            = "#Delta#phi(MET,2^{nd} jet) (rad.)"
plot13.ytit            = "Number of events"
#plot13.xlog            = "yes"
plot13.ylog            = "yes"
plot13.rebin           = 1
plot13.xmin            = 0
plot13.xmax            = 5
plot13.ymin            = 0.001
plot13.ymax            = 10000
#plot13.lpos = "bottom-center"
plot13.name            = "mDeltaPhiMET2ndJet_allPreviousCuts"
plot13.addZUncBand     = zUncBand
plot13.histodata       = h_mDeltaPhiMET2ndJet_DATA


#--- MTenu_PAS (after preselection) ---

h_MTenu_FullPreSel_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection).Clone()
h_MTenu_FullPreSel_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection).Clone()
h_MTenu_FullPreSel_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection).Clone()
h_MTenu_FullPreSel_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection).Clone()
h_MTenu_FullPreSel_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection).Clone()
h_MTenu_FullPreSel_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection).Clone()
h_MTenu_FullPreSel_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection).Clone()
h_MTenu_FullPreSel_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection).Clone()
#h_MTenu_FullPreSel_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection).Clone()
#h_MTenu_FullPreSel_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection).Clone()
h_MTenu_FullPreSel_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection).Clone()
h_MTenu_FullPreSel_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection).Clone()
h_MTenu_FullPreSel_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection).Clone()
h_MTenu_FullPreSel_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_MTenu_FullPreSel_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MTenu_PAS", File_QCD, QCDscaleFactor).Clone()

plot14 = Plot()
## inputs for stacked histograms
## it created h_MTenu_FullPreSel_TTbar, h_MTenu_FullPreSel_TTbar+h_MTenu_FullPreSel_ZJetAlpgen , h_MTenu_FullPreSel_TTbar+h_MTenu_FullPreSel_ZJetAlpgen+h_MTenu_FullPreSel_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
if(UseQCDFromData):
    plot14.histosStack     = [h_MTenu_FullPreSel_QCDFromData, h_MTenu_FullPreSel_TTbar, h_MTenu_FullPreSel_WJetAlpgen, h_MTenu_FullPreSel_OTHERBKG]
    plot14.keysStack       = ["QCD","ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot14.histosStack     = [h_MTenu_FullPreSel_TTbar, h_MTenu_FullPreSel_WJetAlpgen, h_MTenu_FullPreSel_OTHERBKG]
    plot14.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
#plot14.histos          = [h_MTenu_FullPreSel_LQenujj_M100, h_MTenu_FullPreSel_LQenujj_M200, h_MTenu_FullPreSel_LQenujj_M300]
#plot14.keys            = ["LQ e\\nujj M100","LQ e\\nujj M200","LQ e\\nujj M300"]
plot14.histos          = [h_MTenu_FullPreSel_LQenujj_M200]
plot14.keys            = ["LQ e\\nujj M200"]
plot14.xtit            = "M_{T}(e\\nu) (GeV/c^{2})"
plot14.ytit            = "Number of events"
# plot14.ylog            = "yes"
# plot14.rebin           = 1
# plot14.ymin            = 0.00000001
# plot14.ymax            = 20
plot14.ylog            = "no"
plot14.rebin           = 1
plot14.xmin            = 0
plot14.xmax            = 400
plot14.ymin            = 0
plot14.ymax            = 60
#plot14.lpos = "bottom-center"
plot14.name            = "MTenu_allPreviousCuts_ylin"
plot14.addZUncBand     = zUncBand
plot14.histodata       = h_MTenu_FullPreSel_DATA

plot14_ylog = Plot()
## inputs for stacked histograms
## it created h_MTenu_FullPreSel_TTbar, h_MTenu_FullPreSel_TTbar+h_MTenu_FullPreSel_ZJetAlpgen , h_MTenu_FullPreSel_TTbar+h_MTenu_FullPreSel_ZJetAlpgen+h_MTenu_FullPreSel_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
if(UseQCDFromData):
    plot14_ylog.histosStack     = [h_MTenu_FullPreSel_QCDFromData, h_MTenu_FullPreSel_TTbar, h_MTenu_FullPreSel_WJetAlpgen, h_MTenu_FullPreSel_OTHERBKG]
    plot14_ylog.keysStack       = ["QCD", "ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot14_ylog.histosStack     = [h_MTenu_FullPreSel_TTbar, h_MTenu_FullPreSel_WJetAlpgen, h_MTenu_FullPreSel_OTHERBKG]
    plot14_ylog.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]
    

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot14_ylog.histos          = [h_MTenu_FullPreSel_LQenujj_M200, h_MTenu_FullPreSel_LQenujj_M300]
plot14_ylog.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot14_ylog.xtit            = "M_{T}(e\\nu) (GeV/c^{2})"
plot14_ylog.ytit            = "Number of events"
plot14_ylog.ylog            = "yes"
plot14_ylog.rebin           = 1 # don't change it (since a rebinning is already applied above on the same histo)
plot14_ylog.xmin            = 0
plot14_ylog.xmax            = 500
plot14_ylog.ymin            = 0.001
plot14_ylog.ymax            = 100
#plot14_ylog.lpos = "bottom-center"
plot14_ylog.name            = "MTenu_allPreviousCuts"
plot14_ylog.addZUncBand     = zUncBand
plot14_ylog.histodata       = h_MTenu_FullPreSel_DATA



#--- sT ---

h_sT_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
h_sT_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
h_sT_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
h_sT_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
h_sT_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
h_sT_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
h_sT_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
h_sT_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
#h_sT_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
#h_sT_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
h_sT_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
h_sT_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
h_sT_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
h_sT_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_sT_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________sT_PAS", File_QCD, QCDscaleFactor).Clone()

plot15 = Plot()
## inputs for stacked histograms
## it created h_sT_TTbar, h_sT_TTbar+h_sT_ZJetAlpgen , h_sT_TTbar+h_sT_ZJetAlpgen+h_sT_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
if(UseQCDFromData):
    plot15.histosStack     = [h_sT_QCDFromData, h_sT_TTbar, h_sT_WJetAlpgen, h_sT_OTHERBKG]
    plot15.keysStack       = ["QCD","ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot15.histosStack     = [h_sT_TTbar, h_sT_WJetAlpgen, h_sT_OTHERBKG]
    plot15.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot15.histos          = [h_sT_LQenujj_M200, h_sT_LQenujj_M300]
plot15.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot15.xtit            = "St (GeV/c)"
plot15.ytit            = "Number of events"
#plot15.xlog            = "yes"
plot15.ylog            = "yes"
plot15.rebin           = 2
plot15.xmin            = 50
plot15.xmax            = 2000
plot15.ymin            = 0.001
plot15.ymax            = 1000
#plot15.lpos = "bottom-center"
plot15.name            = "sT_allPreviousCuts"
plot15.addZUncBand     = zUncBand
plot15.histodata       = h_sT_DATA

#--- Mjj_PAS (after preselection) ---

h_Mjj_FullPreSel_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
h_Mjj_FullPreSel_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
h_Mjj_FullPreSel_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
h_Mjj_FullPreSel_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
h_Mjj_FullPreSel_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
h_Mjj_FullPreSel_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
h_Mjj_FullPreSel_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
h_Mjj_FullPreSel_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
#h_Mjj_FullPreSel_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
#h_Mjj_FullPreSel_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
h_Mjj_FullPreSel_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
h_Mjj_FullPreSel_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
h_Mjj_FullPreSel_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
h_Mjj_FullPreSel_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_Mjj_FullPreSel_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mjj_PAS", File_QCD, QCDscaleFactor).Clone()

plot16 = Plot()
## inputs for stacked histograms
## it created h_Mjj_FullPreSel_TTbar, h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen , h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen+h_Mjj_FullPreSel_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
if(UseQCDFromData):
    plot16.histosStack     = [h_Mjj_FullPreSel_QCDFromData,h_Mjj_FullPreSel_TTbar, h_Mjj_FullPreSel_WJetAlpgen, h_Mjj_FullPreSel_OTHERBKG]
    plot16.keysStack       = ["QCD","ttbar", "W//W* + jets", otherBkgsKey]
else:
    plot16.histosStack     = [h_Mjj_FullPreSel_TTbar, h_Mjj_FullPreSel_WJetAlpgen, h_Mjj_FullPreSel_OTHERBKG]
    plot16.keysStack       = ["ttbar", "W//W* + jets", otherBkgsKey]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
#plot16.histos          = [h_Mjj_FullPreSel_LQenujj_M100, h_Mjj_FullPreSel_LQenujj_M200, h_Mjj_FullPreSel_LQenujj_M300]
#plot16.keys            = ["LQ enujj M100","LQ enujj M200","LQ enujj M300"]
plot16.histos          = [h_Mjj_FullPreSel_LQenujj_M200, h_Mjj_FullPreSel_LQenujj_M300]
plot16.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot16.xtit            = "M(jj) (GeV/c^{2})"
plot16.ytit            = "Number of events"
# plot16.ylog            = "yes"
# plot16.rebin           = 1
# plot16.ymin            = 0.00000001
# plot16.ymax            = 20
plot16.ylog            = "yes"
plot16.rebin           = 2
plot16.ymin            = 0.001
plot16.ymax            = 500
plot16.xmin            = 0
plot16.xmax            = 2000
#plot16.lpos = "bottom-center"
plot16.name            = "Mjj_FullPreSel_allPreviousCuts_ylin"
plot16.addZUncBand     = zUncBand
plot16.histodata       = h_Mjj_FullPreSel_DATA

# plot16_ylog = Plot()
# ## inputs for stacked histograms
# ## it created h_Mjj_FullPreSel_TTbar, h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen , h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen+h_Mjj_FullPreSel_QCD_Madgraph etc..
# ## and plot them one on top of each other to effectly create a stacked histogram
# if(UseQCDFromData):
#     plot16_ylog.histosStack     = [h_Mjj_FullPreSel_QCDFromData, h_Mjj_FullPreSel_TTbar, h_Mjj_FullPreSel_WJetAlpgen, h_Mjj_FullPreSel_OTHERBKG]
#     plot16_ylog.keysStack       = ["QCD","ttbar", "W/W* + jets", otherBkgsKey]
# else:
#     plot16_ylog.histosStack     = [h_Mjj_FullPreSel_TTbar, h_Mjj_FullPreSel_WJetAlpgen, h_Mjj_FullPreSel_OTHERBKG]
#     plot16_ylog.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]


# ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot16_ylog.histos          = [h_Mjj_FullPreSel_LQenujj_M200, h_Mjj_FullPreSel_LQenujj_M300]
# plot16_ylog.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
# plot16_ylog.xtit            = "M(jj) (GeV/c^{2})"
# plot16_ylog.ytit            = "Number of events"
# plot16_ylog.ylog            = "yes"
# plot16_ylog.rebin           = 1 # don't change it (since a rebinning is already applied above on the same histo)
# plot16_ylog.ymin            = 0.001
# plot16_ylog.ymax            = 500
# plot16_ylog.xmin            = 0
# plot16_ylog.xmax            = 2000
# #plot16_ylog.lpos = "bottom-center"
# plot16_ylog.name            = "Mjj_FullPreSel_allPreviousCuts"
# plot16_ylog.addZUncBand     = zUncBand
# plot16_ylog.histodata       = h_Mjj_FullPreSel_DATA


##--- Mej preselection

h_Mej_presel_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________Mej_1stPair_PAS", File_preselection).Clone()
h_Mej_presel_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________Mej_1stPair_PAS", File_preselection).Clone()
h_Mej_presel_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________Mej_1stPair_PAS", File_preselection).Clone()
h_Mej_presel_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________Mej_1stPair_PAS", File_preselection).Clone()
h_Mej_presel_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________Mej_1stPair_PAS", File_preselection).Clone()
h_Mej_presel_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Mej_1stPair_PAS", File_preselection).Clone()
h_Mej_presel_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mej_1stPair_PAS", File_preselection).Clone()
h_Mej_presel_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Mej_1stPair_PAS", File_preselection).Clone()
#h_Mej_presel_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Mej_1stPair_PAS", File_preselection).Clone()
#h_Mej_presel_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Mej_1stPair_PAS", File_preselection).Clone()
h_Mej_presel_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Mej_1stPair_PAS", File_preselection).Clone()
h_Mej_presel_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Mej_1stPair_PAS", File_preselection).Clone()
h_Mej_presel_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Mej_1stPair_PAS", File_preselection).Clone()
h_Mej_presel_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mej_1stPair_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_Mej_presel_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mej_1stPair_PAS", File_QCD, QCDscaleFactor).Clone()

h_Mej_presel_LQenujj_M100.Add(GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________Mej_2ndPair_PAS", File_preselection))
h_Mej_presel_LQenujj_M200.Add(GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________Mej_2ndPair_PAS", File_preselection))
h_Mej_presel_LQenujj_M300.Add(GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________Mej_2ndPair_PAS", File_preselection))
h_Mej_presel_LQenujj_M400.Add(GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________Mej_2ndPair_PAS", File_preselection))
h_Mej_presel_LQenujj_M500.Add(GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________Mej_2ndPair_PAS", File_preselection))
h_Mej_presel_TTbar.Add(GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Mej_2ndPair_PAS", File_preselection))
h_Mej_presel_ZJetAlpgen.Add(GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mej_2ndPair_PAS", File_preselection))
h_Mej_presel_OTHERBKG.Add(GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Mej_2ndPair_PAS", File_preselection))
#h_Mej_presel_QCD_Madgraph.Add(GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Mej_2ndPair_PAS", File_preselection))
#h_Mej_presel_QCDPt15.Add(GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Mej_2ndPair_PAS", File_preselection))
h_Mej_presel_SingleTop.Add(GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Mej_2ndPair_PAS", File_preselection))
h_Mej_presel_VVjets.Add(GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Mej_2ndPair_PAS", File_preselection))
h_Mej_presel_WJetAlpgen.Add(GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Mej_2ndPair_PAS", File_preselection))
h_Mej_presel_DATA.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mej_2ndPair_PAS", File_preselection))
if(UseQCDFromData):
    h_Mej_presel_QCDFromData.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mej_2ndPair_PAS", File_QCD, QCDscaleFactor))

plot17 = Plot()
if(UseQCDFromData):
    plot17.histosStack     = [h_Mej_presel_QCDFromData,h_Mej_presel_TTbar, h_Mej_presel_WJetAlpgen, h_Mej_presel_OTHERBKG]
    plot17.keysStack       = ["QCD","ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot17.histosStack     = [h_Mej_presel_TTbar, h_Mej_presel_WJetAlpgen, h_Mej_presel_OTHERBKG]
    plot17.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot17.histos          = [h_Mej_presel_LQenujj_M200, h_Mej_presel_LQenujj_M300]
plot17.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot17.xtit            = "M(ej) (GeV/c^{2})"
plot17.ytit            = "Number of events x 2"
plot17.ylog            = "yes"
plot17.rebin           = 2
plot17.xmin            = 0
plot17.xmax            = 1000
plot17.ymin            = 0.001
plot17.ymax            = 500
#plot17.lpos = "bottom-center"
plot17.name            = "Mej_allPreviousCuts"
plot17.addZUncBand     = zUncBand
plot17.histodata       = h_Mej_presel_DATA


##--- MTnuj preselection

h_MTnuj_presel_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________MTnuj_1stPair_PAS", File_preselection).Clone()
h_MTnuj_presel_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________MTnuj_1stPair_PAS", File_preselection).Clone()
h_MTnuj_presel_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________MTnuj_1stPair_PAS", File_preselection).Clone()
h_MTnuj_presel_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________MTnuj_1stPair_PAS", File_preselection).Clone()
h_MTnuj_presel_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________MTnuj_1stPair_PAS", File_preselection).Clone()
h_MTnuj_presel_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________MTnuj_1stPair_PAS", File_preselection).Clone()
h_MTnuj_presel_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________MTnuj_1stPair_PAS", File_preselection).Clone()
h_MTnuj_presel_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________MTnuj_1stPair_PAS", File_preselection).Clone()
#h_MTnuj_presel_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________MTnuj_1stPair_PAS", File_preselection).Clone()
#h_MTnuj_presel_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________MTnuj_1stPair_PAS", File_preselection).Clone()
h_MTnuj_presel_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________MTnuj_1stPair_PAS", File_preselection).Clone()
h_MTnuj_presel_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________MTnuj_1stPair_PAS", File_preselection).Clone()
h_MTnuj_presel_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________MTnuj_1stPair_PAS", File_preselection).Clone()
h_MTnuj_presel_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MTnuj_1stPair_PAS", File_preselection).Clone()
if(UseQCDFromData):
    h_MTnuj_presel_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MTnuj_1stPair_PAS", File_QCD, QCDscaleFactor).Clone()

h_MTnuj_presel_LQenujj_M100.Add(GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________MTnuj_2ndPair_PAS", File_preselection))
h_MTnuj_presel_LQenujj_M200.Add(GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________MTnuj_2ndPair_PAS", File_preselection))
h_MTnuj_presel_LQenujj_M300.Add(GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________MTnuj_2ndPair_PAS", File_preselection))
h_MTnuj_presel_LQenujj_M400.Add(GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________MTnuj_2ndPair_PAS", File_preselection))
h_MTnuj_presel_LQenujj_M500.Add(GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________MTnuj_2ndPair_PAS", File_preselection))
h_MTnuj_presel_TTbar.Add(GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________MTnuj_2ndPair_PAS", File_preselection))
h_MTnuj_presel_ZJetAlpgen.Add(GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________MTnuj_2ndPair_PAS", File_preselection))
h_MTnuj_presel_OTHERBKG.Add(GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________MTnuj_2ndPair_PAS", File_preselection))
#h_MTnuj_presel_QCD_Madgraph.Add(GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________MTnuj_2ndPair_PAS", File_preselection))
#h_MTnuj_presel_QCDPt15.Add(GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________MTnuj_2ndPair_PAS", File_preselection))
h_MTnuj_presel_SingleTop.Add(GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________MTnuj_2ndPair_PAS", File_preselection))
h_MTnuj_presel_VVjets.Add(GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________MTnuj_2ndPair_PAS", File_preselection))
h_MTnuj_presel_WJetAlpgen.Add(GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________MTnuj_2ndPair_PAS", File_preselection))
h_MTnuj_presel_DATA.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MTnuj_2ndPair_PAS", File_preselection))
if(UseQCDFromData):
    h_MTnuj_presel_QCDFromData.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MTnuj_2ndPair_PAS", File_QCD, QCDscaleFactor))

plot18 = Plot()
if(UseQCDFromData):
    plot18.histosStack     = [h_MTnuj_presel_QCDFromData, h_MTnuj_presel_TTbar, h_MTnuj_presel_WJetAlpgen, h_MTnuj_presel_OTHERBKG]
    plot18.keysStack       = ["QCD", "ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot18.histosStack     = [h_MTnuj_presel_TTbar, h_MTnuj_presel_WJetAlpgen, h_MTnuj_presel_OTHERBKG]
    plot18.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot18.histos          = [h_MTnuj_presel_LQenujj_M200, h_MTnuj_presel_LQenujj_M300]
plot18.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot18.xtit            = "M_{T}(\\nuj) (GeV/c^{2})"
plot18.ytit            = "Number of events x 2"
plot18.ylog            = "yes"
plot18.rebin           = 2
plot18.xmin            = 0
plot18.xmax            = 1000
plot18.ymin            = 0.001
plot18.ymax            = 500
#plot18.lpos = "bottom-center"
plot18.name            = "MTnuj_allPreviousCuts"
plot18.addZUncBand     = zUncBand
plot18.histodata       = h_MTnuj_presel_DATA



# ############################ Plots below to be done after full selection ######################

##--- sT full selection ---

plot20 = Plot()

if(UseQCDFromData):
    plot20.histosStack     = [
        GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________sT", File_QCD, QCDscaleFactor).Clone(),
        GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
        GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
        GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
        #     GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
        #     GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
        #     GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
        #     GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________sT", File_selection).Clone()
        ]
    plot20.keysStack       = [
        "QCD","ttbar", "W/W* + jets", otherBkgsKey]
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
else:
    plot20.histosStack     = [
        GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
        GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
        GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
        #     GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
        #     GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
        #     GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
        #     GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________sT", File_selection).Clone()
        ]
    plot20.keysStack       = [
        "ttbar", "W/W* + jets", otherBkgsKey]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot20.histos          = [
    GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
    GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________sT", File_selection).Clone()
    ]
plot20.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot20.xtit            = "St (GeV/c)"
plot20.ytit            = "Number of events"
plot20.ylog            = "yes"
plot20.rebin           = 10
plot20.xmin            = 100
plot20.xmax            = 1000
plot20.ymin            = 0.001
plot20.ymax            = 100
#plot20.lpos = "bottom-center"
plot20.name            = "sT_fullSelection"
plot20.addZUncBand     = zUncBand
plot20.histodata       = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________sT", File_selection).Clone()


##--- Mej fullselection --

h_Mej_fullsel_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________Mej_1stPair", File_selection).Clone()
h_Mej_fullsel_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________Mej_1stPair", File_selection).Clone()
h_Mej_fullsel_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________Mej_1stPair", File_selection).Clone()
h_Mej_fullsel_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________Mej_1stPair", File_selection).Clone()
h_Mej_fullsel_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________Mej_1stPair", File_selection).Clone()
h_Mej_fullsel_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Mej_1stPair", File_selection).Clone()
h_Mej_fullsel_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mej_1stPair", File_selection).Clone()
h_Mej_fullsel_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Mej_1stPair", File_selection).Clone()
#h_Mej_fullsel_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Mej_1stPair", File_selection).Clone()
#h_Mej_fullsel_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Mej_1stPair", File_selection).Clone()
h_Mej_fullsel_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Mej_1stPair", File_selection).Clone()
h_Mej_fullsel_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Mej_1stPair", File_selection).Clone()
h_Mej_fullsel_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Mej_1stPair", File_selection).Clone()
h_Mej_fullsel_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mej_1stPair", File_selection).Clone()
if(UseQCDFromData):
    h_Mej_fullsel_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mej_1stPair", File_QCD, QCDscaleFactor).Clone()

h_Mej_fullsel_LQenujj_M100.Add(GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________Mej_2ndPair", File_selection))
h_Mej_fullsel_LQenujj_M200.Add(GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________Mej_2ndPair", File_selection))
h_Mej_fullsel_LQenujj_M300.Add(GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________Mej_2ndPair", File_selection))
h_Mej_fullsel_LQenujj_M400.Add(GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________Mej_2ndPair", File_selection))
h_Mej_fullsel_LQenujj_M500.Add(GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________Mej_2ndPair", File_selection))
h_Mej_fullsel_TTbar.Add(GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Mej_2ndPair", File_selection))
h_Mej_fullsel_ZJetAlpgen.Add(GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mej_2ndPair", File_selection))
h_Mej_fullsel_OTHERBKG.Add(GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Mej_2ndPair", File_selection))
#h_Mej_fullsel_QCD_Madgraph.Add(GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Mej_2ndPair", File_selection))
#h_Mej_fullsel_QCDPt15.Add(GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Mej_2ndPair", File_selection))
h_Mej_fullsel_SingleTop.Add(GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Mej_2ndPair", File_selection))
h_Mej_fullsel_VVjets.Add(GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Mej_2ndPair", File_selection))
h_Mej_fullsel_WJetAlpgen.Add(GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Mej_2ndPair", File_selection))
h_Mej_fullsel_DATA.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mej_2ndPair", File_selection))
if(UseQCDFromData):
    h_Mej_fullsel_QCDFromData.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mej_2ndPair", File_QCD, QCDscaleFactor))

plot21 = Plot()
if(UseQCDFromData):
    plot21.histosStack     = [h_Mej_fullsel_QCDFromData, h_Mej_fullsel_TTbar, h_Mej_fullsel_WJetAlpgen, h_Mej_fullsel_OTHERBKG]
    plot21.keysStack       = ["QCD","ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot21.histosStack     = [h_Mej_fullsel_TTbar, h_Mej_fullsel_WJetAlpgen, h_Mej_fullsel_OTHERBKG]
    plot21.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot21.histos          = [h_Mej_fullsel_LQenujj_M200, h_Mej_fullsel_LQenujj_M300]
plot21.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot21.xtit            = "M(ej) (GeV/c^{2})"
plot21.ytit            = "Number of events x 2"
plot21.ylog            = "yes"
plot21.rebin           = 10
plot21.xmin            = 0
plot21.xmax            = 1000
plot21.ymin            = 0.001
plot21.ymax            = 500
#plot21.lpos = "bottom-center"
plot21.name            = "Mej_fullSelection"
plot21.addZUncBand     = zUncBand
plot21.histodata       = h_Mej_fullsel_DATA


##--- MTnuj fullselection --

h_MTnuj_fullsel_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________MTnuj_1stPair", File_selection).Clone()
h_MTnuj_fullsel_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________MTnuj_1stPair", File_selection).Clone()
h_MTnuj_fullsel_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________MTnuj_1stPair", File_selection).Clone()
h_MTnuj_fullsel_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________MTnuj_1stPair", File_selection).Clone()
h_MTnuj_fullsel_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________MTnuj_1stPair", File_selection).Clone()
h_MTnuj_fullsel_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________MTnuj_1stPair", File_selection).Clone()
h_MTnuj_fullsel_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________MTnuj_1stPair", File_selection).Clone()
h_MTnuj_fullsel_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________MTnuj_1stPair", File_selection).Clone()
#h_MTnuj_fullsel_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________MTnuj_1stPair", File_selection).Clone()
#h_MTnuj_fullsel_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________MTnuj_1stPair", File_selection).Clone()
h_MTnuj_fullsel_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________MTnuj_1stPair", File_selection).Clone()
h_MTnuj_fullsel_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________MTnuj_1stPair", File_selection).Clone()
h_MTnuj_fullsel_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________MTnuj_1stPair", File_selection).Clone()
h_MTnuj_fullsel_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MTnuj_1stPair", File_selection).Clone()
if(UseQCDFromData):
    h_MTnuj_fullsel_QCDFromData = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MTnuj_1stPair", File_QCD, QCDscaleFactor).Clone()

h_MTnuj_fullsel_LQenujj_M100.Add(GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________MTnuj_2ndPair", File_selection))
h_MTnuj_fullsel_LQenujj_M200.Add(GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________MTnuj_2ndPair", File_selection))
h_MTnuj_fullsel_LQenujj_M300.Add(GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________MTnuj_2ndPair", File_selection))
h_MTnuj_fullsel_LQenujj_M400.Add(GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________MTnuj_2ndPair", File_selection))
h_MTnuj_fullsel_LQenujj_M500.Add(GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________MTnuj_2ndPair", File_selection))
h_MTnuj_fullsel_TTbar.Add(GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________MTnuj_2ndPair", File_selection))
h_MTnuj_fullsel_ZJetAlpgen.Add(GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________MTnuj_2ndPair", File_selection))
h_MTnuj_fullsel_OTHERBKG.Add(GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________MTnuj_2ndPair", File_selection))
#h_MTnuj_fullsel_QCD_Madgraph.Add(GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________MTnuj_2ndPair", File_selection))
#h_MTnuj_fullsel_QCDPt15.Add(GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________MTnuj_2ndPair", File_selection))
h_MTnuj_fullsel_SingleTop.Add(GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________MTnuj_2ndPair", File_selection))
h_MTnuj_fullsel_VVjets.Add(GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________MTnuj_2ndPair", File_selection))
h_MTnuj_fullsel_WJetAlpgen.Add(GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________MTnuj_2ndPair", File_selection))
h_MTnuj_fullsel_DATA.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MTnuj_2ndPair", File_selection))
if(UseQCDFromData):
    h_MTnuj_fullsel_QCDFromData.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MTnuj_2ndPair", File_QCD, QCDscaleFactor))

plot22 = Plot()
if(UseQCDFromData):
    plot22.histosStack     = [h_MTnuj_fullsel_QCDFromData, h_MTnuj_fullsel_TTbar, h_MTnuj_fullsel_WJetAlpgen, h_MTnuj_fullsel_OTHERBKG]
    plot22.keysStack       = ["QCD", "ttbar", "W/W* + jets", otherBkgsKey]
else:
    plot22.histosStack     = [h_MTnuj_fullsel_TTbar, h_MTnuj_fullsel_WJetAlpgen, h_MTnuj_fullsel_OTHERBKG]
    plot22.keysStack       = ["ttbar", "W/W* + jets", otherBkgsKey]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot22.histos          = [h_MTnuj_fullsel_LQenujj_M200, h_MTnuj_fullsel_LQenujj_M300]
plot22.keys            = ["LQ e\\nujj M200","LQ e\\nujj M300"]
plot22.xtit            = "M_{T}(\\nuj) (GeV/c^{2})"
plot22.ytit            = "Number of events x 2"
plot22.ylog            = "yes"
plot22.rebin           = 10
plot22.xmin            = 0
plot22.xmax            = 1000
plot22.ymin            = 0.001
plot22.ymax            = 500
#plot22.lpos = "bottom-center"
plot22.name            = "MTnuj_fullSelection"
plot22.addZUncBand     = zUncBand
plot22.histodata       = h_MTnuj_fullsel_DATA

#-----------------------------------------------------------------------------------


# List of plots to be plotted
plots = [plot0, plot1, plot2, plot3, plot4, plot6, plot7, plot8, plot9, plot6and8, plot7and9, 
         plot10, plot11, plot12, plot13,
         plot14, plot14_ylog, plot15, plot16, #plot16_ylog,
         plot17, plot18, # produced using preselection root file
         plot20, plot21, plot22] # produced using full selection root file



############# USER CODE - END ################################################
##############################################################################


#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()
c.Print("allPlots.ps[")
for plot in plots:
    plot.Draw("allPlots.ps")
c.Print("allPlots.ps]")
