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
    lint        = "3.04 pb^{-1}" # integrated luminosity of the sample ( example "10 pb^{-1}" )
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

File_preselection = GetFile("$LQDATA/muons/3.04pb-1/output_cutTable_mmjjSample/analysisClass_mmjjSample_plots.root")
File_selection    = File_preselection

#### Common values for plots:
#otherBkgsKey="QCD, single top, VV+jets, W/W*+jets"
otherBkgsKey="Other Bkgs"
zUncBand="no"

pt_xmin=0
pt_xmax=400
pt_ymin=0.001
pt_ymax=300

eta_rebin=10
eta_ymin=0
eta_ymax=30



#--- Final plots are defined here

# Simply add or remove plots and update the list plots = [plot0, plot1, ...] below

#--- Mmm_TwoMuOnly ---

h_Mmm_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Mmm_TwoMuOnly", File_preselection).Clone()
h_Mmm_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mmm_TwoMuOnly", File_preselection).Clone()
h_Mmm_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mmm_TwoMuOnly", File_preselection).Clone()

plot0 = Plot()
plot0.histosStack     = [h_Mmm_TTbar, h_Mmm_ZJetAlpgen]
plot0.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot0.xtit            = "M(mm) (GeV/c^{2})"
plot0.ytit            = "Number of events"
# plot0.ylog            = "yes"
# plot0.rebin           = 1
# plot0.ymin            = 0.00000001
# plot0.ymax            = 20
plot0.ylog            = "no"
plot0.rebin           = 1
plot0.ymin            = 0
plot0.ymax            = 700
plot0.xmin            = 0
plot0.xmax            = 200
#plot0.lpos = "bottom-center"
plot0.name            = "Mmm_allPreviousCuts_ylin"
plot0.addZUncBand     = zUncBand
plot0.histodata       = h_Mmm_DATA

plot0_ylog = Plot()
plot0_ylog.histosStack     = [h_Mmm_TTbar, h_Mmm_ZJetAlpgen]
plot0_ylog.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot0_ylog.xtit            = "M(mm) (GeV/c^{2})"
plot0_ylog.ytit            = "Number of events"
plot0_ylog.ylog            = "yes"
plot0_ylog.rebin           = 1
plot0_ylog.ymin            = 0.001
plot0_ylog.ymax            = 1000
plot0_ylog.xmin            = 0
plot0_ylog.xmax            = 1000
#plot0_ylog.lpos = "bottom-center"
plot0_ylog.name            = "Mmm_allPreviousCuts"
plot0_ylog.addZUncBand     = zUncBand
plot0_ylog.histodata       = h_Mmm_DATA


# #--- nMu_PtCut_IDISO_noOvrlp ---

h_nMu_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________nMu_PtCut_IDISO", File_preselection).Clone()
h_nMu_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________nMu_PtCut_IDISO", File_preselection).Clone()
h_nMu_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________nMu_PtCut_IDISO", File_preselection).Clone()

plot1 = Plot()
## inputs for stacked histograms
plot1.histosStack     = [h_nMu_TTbar, h_nMu_ZJetAlpgen]
plot1.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot1.xtit            = "Number of muons"
plot1.ytit            = "Number of events"
plot1.ylog            = "yes"
plot1.rebin           = 1
plot1.ymin            = 0.0001
plot1.ymax            = 60000000
#plot1.lpos = "bottom-center"
plot1.name            = "nMu_allPreviousCuts"
plot1.addZUncBand     = zUncBand
plot1.histodata       = h_nMu_DATA




#--- Pt1stMu_PAS  ---

h_pT1stMu_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()
h_pT1stMu_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()
h_pT1stMu_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()

plot2 = Plot()
## inputs for stacked histograms
plot2.histosStack     = [h_pT1stMu_TTbar, h_pT1stMu_ZJetAlpgen]
plot2.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot2.xtit            = "pT 1st muon (GeV/c)"
plot2.ytit            = "Number of events"
plot2.ylog            = "yes"
plot2.rebin           = 1
plot2.xmin            = pt_xmin
plot2.xmax            = pt_xmax
plot2.ymin            = pt_ymin
plot2.ymax            = pt_ymax
#plot2.lpos = "bottom-center"
plot2.name            = "pT1stMu_allPreviousCuts"
plot2.addZUncBand     = zUncBand
plot2.histodata       = h_pT1stMu_DATA


#--- Eta1stMu_PAS  ---

h_Eta1stMu_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()
h_Eta1stMu_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()
h_Eta1stMu_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()

plot3 = Plot()
## inputs for stacked histograms
plot3.histosStack     = [h_Eta1stMu_TTbar, h_Eta1stMu_ZJetAlpgen]
plot3.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot3.xtit            = "#eta 1st muon"
plot3.ytit            = "Number of events"
plot3.rebin           = eta_rebin
plot3.ymin            = eta_ymin
plot3.ymax            = eta_ymax
plot3.lpos = "top-left"
plot3.name            = "Eta1stMu_allPreviousCuts"
plot3.addZUncBand     = zUncBand
plot3.histodata       = h_Eta1stMu_DATA



#--- Pt2ndMu_PAS  ---

h_pT2ndMu_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection).Clone()
h_pT2ndMu_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection).Clone()
h_pT2ndMu_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection).Clone()

plot4 = Plot()
## inputs for stacked histograms
plot4.histosStack     = [h_pT2ndMu_TTbar, h_pT2ndMu_ZJetAlpgen]
plot4.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot4.xtit            = "pT 2nd muon (GeV/c)"
plot4.ytit            = "Number of events"
plot4.ylog            = "yes"
plot4.rebin           = 1
plot4.xmin            = pt_xmin
plot4.xmax            = pt_xmax
plot4.ymin            = pt_ymin
plot4.ymax            = pt_ymax
#plot4.lpos = "bottom-center"
plot4.name            = "pT2ndMu_allPreviousCuts"
plot4.addZUncBand     = zUncBand
plot4.histodata       = h_pT2ndMu_DATA


#--- Eta2ndMu_PAS  ---

h_Eta2ndMu_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection).Clone()
h_Eta2ndMu_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection).Clone()
h_Eta2ndMu_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection).Clone()

plot5 = Plot()
## inputs for stacked histograms
plot5.histosStack     = [h_Eta2ndMu_TTbar, h_Eta2ndMu_ZJetAlpgen]
plot5.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot5.xtit            = "#eta 2nd muon"
plot5.ytit            = "Number of events"
plot5.rebin           = eta_rebin
plot5.ymin            = eta_ymin
plot5.ymax            = eta_ymax
plot5.lpos = "top-left"
plot5.name            = "Eta2ndMu_allPreviousCuts"
plot5.addZUncBand     = zUncBand
plot5.histodata       = h_Eta2ndMu_DATA


# #--- nJet_TwoMuOnly_EtaCut ---

# h_nJet_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________nJet_TwoMuOnly_EtaCut", File_preselection).Clone()
# h_nJet_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________nJet_TwoMuOnly_EtaCut", File_preselection).Clone()
# h_nJet_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________nJet_TwoMuOnly_EtaCut", File_preselection).Clone()
# h_nJet_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________nJet_TwoMuOnly_EtaCut", File_preselection).Clone()
# h_nJet_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________nJet_TwoMuOnly_EtaCut", File_preselection).Clone()
h_nJet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________nJet_TwoMuOnly_EtaCut", File_preselection).Clone()
h_nJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________nJet_TwoMuOnly_EtaCut", File_preselection).Clone()
# h_nJet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________nJet_TwoMuOnly_EtaCut", File_preselection).Clone()
# #h_nJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________nJet_TwoMuOnly_EtaCut", File_preselection).Clone()
# #h_nJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________nJet_TwoMuOnly_EtaCut", File_preselection).Clone()
# h_nJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________nJet_TwoMuOnly_EtaCut", File_preselection).Clone()
# h_nJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________nJet_TwoMuOnly_EtaCut", File_preselection).Clone()
# h_nJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________nJet_TwoMuOnly_EtaCut", File_preselection).Clone()
h_nJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________nJet_TwoMuOnly_EtaCut", File_preselection).Clone()


plot6 = Plot()
## inputs for stacked histograms
plot6.histosStack     = [h_nJet_TTbar, h_nJet_ZJetAlpgen]
plot6.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot6.histos          = [h_nJet_LQmmjj_M200, h_nJet_LQmmjj_M300]
# plot6.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot6.xtit            = "Number of jets"
plot6.ytit            = "Number of events"
plot6.ylog            = "yes"
plot6.rebin           = 1
plot6.xmin            = 0
plot6.xmax            = 12
plot6.ymin            = 0.01
plot6.ymax            = 2000
#plot6.lpos = "bottom-center"
plot6.name            = "nJet_allPreviousCuts"
plot6.addZUncBand     = zUncBand
plot6.histodata       = h_nJet_DATA



# #--- Pt1stJet_PAS ---

# h_Pt1stJet_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_Pt1stJet_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_Pt1stJet_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_Pt1stJet_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_Pt1stJet_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_Pt1stJet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_Pt1stJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_Pt1stJet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# #h_Pt1stJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# #h_Pt1stJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_Pt1stJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_Pt1stJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_Pt1stJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_Pt1stJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()


plot7 = Plot()
## inputs for stacked histograms
plot7.histosStack     = [h_Pt1stJet_TTbar, h_Pt1stJet_ZJetAlpgen]
plot7.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot7.histos          = [h_Pt1stJet_LQmmjj_M200, h_Pt1stJet_LQmmjj_M300]
# plot7.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot7.xtit            = "pT 1st jet (GeV/c)"
plot7.ytit            = "Number of events"
plot7.ylog            = "yes"
plot7.rebin           = 1
plot7.xmin            = pt_xmin
plot7.xmax            = pt_xmax
plot7.ymin            = pt_ymin
plot7.ymax            = pt_ymax
#plot7.lpos = "bottom-center"
plot7.name            = "Pt1stJet_allPreviousCuts"
plot7.addZUncBand     = zUncBand
plot7.histodata       = h_Pt1stJet_DATA


# #--- Eta1stJet_PAS ---

# h_Eta1stJet_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_Eta1stJet_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_Eta1stJet_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_Eta1stJet_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_Eta1stJet_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_Eta1stJet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_Eta1stJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_Eta1stJet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# #h_Eta1stJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# #h_Eta1stJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_Eta1stJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_Eta1stJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_Eta1stJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_Eta1stJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()


plot8 = Plot()
## inputs for stacked histograms
plot8.histosStack     = [h_Eta1stJet_TTbar, h_Eta1stJet_ZJetAlpgen]
plot8.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot8.histos          = [h_Eta1stJet_LQmmjj_M200, h_Eta1stJet_LQmmjj_M300]
# plot8.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot8.xtit            = "#eta 1st jet"
plot8.ytit            = "Number of events"
plot8.rebin           = eta_rebin
plot8.ymin            = eta_ymin
plot8.ymax            = eta_ymax
plot8.lpos = "top-left"
plot8.name            = "Eta1stJet_allPreviousCuts"
plot8.addZUncBand     = zUncBand
plot8.histodata       = h_Eta1stJet_DATA


# #--- Pt2ndJet_PAS ---

# h_Pt2ndJet_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
# h_Pt2ndJet_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
# h_Pt2ndJet_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
# h_Pt2ndJet_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
# h_Pt2ndJet_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
h_Pt2ndJet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
h_Pt2ndJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
# h_Pt2ndJet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
# #h_Pt2ndJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
# #h_Pt2ndJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
# h_Pt2ndJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
# h_Pt2ndJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
# h_Pt2ndJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()
h_Pt2ndJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection).Clone()


plot9 = Plot()
## inputs for stacked histograms
plot9.histosStack     = [h_Pt2ndJet_TTbar, h_Pt2ndJet_ZJetAlpgen]
plot9.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot9.histos          = [h_Pt2ndJet_LQmmjj_M200, h_Pt2ndJet_LQmmjj_M300]
# plot9.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot9.xtit            = "pT 2nd jet (GeV/c)"
plot9.ytit            = "Number of events"
plot9.ylog            = "yes"
plot9.xmin            = pt_xmin
plot9.xmax            = pt_xmax
plot9.ymin            = pt_ymin
plot9.ymax            = pt_ymax
#plot9.lpos = "bottom-center"
plot9.name            = "Pt2ndJet_allPreviousCuts"
plot9.addZUncBand     = zUncBand
plot9.histodata       = h_Pt2ndJet_DATA


# #--- Eta2ndJet_PAS ---

# h_Eta2ndJet_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
# h_Eta2ndJet_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
# h_Eta2ndJet_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
# h_Eta2ndJet_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
# h_Eta2ndJet_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
h_Eta2ndJet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
h_Eta2ndJet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
# h_Eta2ndJet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
# #h_Eta2ndJet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
# #h_Eta2ndJet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
# h_Eta2ndJet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
# h_Eta2ndJet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
# h_Eta2ndJet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()
h_Eta2ndJet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection).Clone()


plot10 = Plot()
## inputs for stacked histograms
plot10.histosStack     = [h_Eta2ndJet_TTbar, h_Eta2ndJet_ZJetAlpgen]
plot10.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot10.histos          = [h_Eta2ndJet_LQmmjj_M200, h_Eta2ndJet_LQmmjj_M300]
# plot10.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot10.xtit            = "#eta 2nd jet"
plot10.ytit            = "Number of events"
plot10.rebin           = eta_rebin
plot10.ymin            = eta_ymin
plot10.ymax            = eta_ymax
plot10.lpos = "top-left"
plot10.name            = "Eta2ndJet_allPreviousCuts"
plot10.addZUncBand     = zUncBand
plot10.histodata       = h_Eta2ndJet_DATA

# #--- sT ---

# h_sT_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
# h_sT_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
# h_sT_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
# h_sT_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
# h_sT_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
h_sT_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
h_sT_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
# h_sT_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
# #h_sT_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
# #h_sT_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
# h_sT_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
# h_sT_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
# h_sT_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()
h_sT_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________sT_PAS", File_preselection).Clone()


plot11 = Plot()
## inputs for stacked histograms
## it created h_sT_TTbar, h_sT_TTbar+h_sT_ZJetAlpgen , h_sT_TTbar+h_sT_ZJetAlpgen+h_sT_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
plot11.histosStack     = [h_sT_TTbar, h_sT_ZJetAlpgen]
plot11.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot11.histos          = [h_sT_LQmmjj_M200, h_sT_LQmmjj_M300]
# plot11.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot11.xtit            = "St (GeV/c)"
plot11.ytit            = "Number of events"
#plot11.xlog            = "yes"
plot11.ylog            = "yes"
plot11.rebin           = 2
plot11.xmin            = 50
plot11.xmax            = 1000
plot11.ymin            = 0.001
plot11.ymax            = 100
#plot11.lpos = "bottom-center"
plot11.name            = "sT_allPreviousCuts"
plot11.addZUncBand     = zUncBand
plot11.histodata       = h_sT_DATA


# #--- sTmuon ---

# h_sT_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________sTmuon_PAS", File_preselection).Clone()
# h_sT_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________sTmuon_PAS", File_preselection).Clone()
# h_sT_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________sTmuon_PAS", File_preselection).Clone()
# h_sT_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________sTmuon_PAS", File_preselection).Clone()
# h_sT_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________sTmuon_PAS", File_preselection).Clone()
# h_sT_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________sTmuon_PAS", File_preselection).Clone()
h_sT_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________sTmuon_PAS", File_preselection).Clone()
h_sT_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________sTmuon_PAS", File_preselection).Clone()
# #h_sT_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________sTmuon_PAS", File_preselection).Clone()
# #h_sT_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________sTmuon_PAS", File_preselection).Clone()
# h_sT_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________sTmuon_PAS", File_preselection).Clone()
# h_sT_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________sTmuon_PAS", File_preselection).Clone()
# h_sT_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________sTmuon_PAS", File_preselection).Clone()
h_sT_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________sTmuon_PAS", File_preselection).Clone()


plot11_ele = Plot()
## inputs for stacked histograms
## it created h_sT_TTbar, h_sT_TTbar+h_sT_ZJetAlpgen , h_sT_TTbar+h_sT_ZJetAlpgen+h_sT_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
plot11_ele.histosStack     = [h_sT_TTbar, h_sT_ZJetAlpgen]
plot11_ele.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot11_ele.histos          = [h_sT_LQmmjj_M200, h_sT_LQmmjj_M300]
# plot11_ele.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot11_ele.xtit            = "St muons (GeV/c)"
plot11_ele.ytit            = "Number of events"
#plot11_ele.xlog            = "yes"
plot11_ele.ylog            = "yes"
plot11_ele.rebin           = 2
plot11_ele.xmin            = 50
plot11_ele.xmax            = 1000
plot11_ele.ymin            = 0.001
plot11_ele.ymax            = 100
#plot11_ele.lpos = "bottom-center"
plot11_ele.name            = "sTmuon_allPreviousCuts"
plot11_ele.addZUncBand     = zUncBand
plot11_ele.histodata       = h_sT_DATA


# #--- sTjet ---

# h_sT_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________sTjet_PAS", File_preselection).Clone()
# h_sT_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________sTjet_PAS", File_preselection).Clone()
# h_sT_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________sTjet_PAS", File_preselection).Clone()
# h_sT_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________sTjet_PAS", File_preselection).Clone()
# h_sT_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________sTjet_PAS", File_preselection).Clone()
h_sT_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________sTjet_PAS", File_preselection).Clone()
h_sT_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________sTjet_PAS", File_preselection).Clone()
# h_sT_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________sTjet_PAS", File_preselection).Clone()
# #h_sT_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________sTjet_PAS", File_preselection).Clone()
# #h_sT_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________sTjet_PAS", File_preselection).Clone()
# h_sT_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________sTjet_PAS", File_preselection).Clone()
# h_sT_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________sTjet_PAS", File_preselection).Clone()
# h_sT_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________sTjet_PAS", File_preselection).Clone()
h_sT_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________sTjet_PAS", File_preselection).Clone()


plot11_jet = Plot()
## inputs for stacked histograms
## it created h_sT_TTbar, h_sT_TTbar+h_sT_ZJetAlpgen , h_sT_TTbar+h_sT_ZJetAlpgen+h_sT_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
plot11_jet.histosStack     = [h_sT_TTbar, h_sT_ZJetAlpgen]
plot11_jet.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot11_jet.histos          = [h_sT_LQmmjj_M200, h_sT_LQmmjj_M300]
# plot11_jet.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot11_jet.xtit            = "St jets (GeV/c)"
plot11_jet.ytit            = "Number of events"
#plot11_jet.xlog            = "yes"
plot11_jet.ylog            = "yes"
plot11_jet.rebin           = 2
plot11_jet.xmin            = 50
plot11_jet.xmax            = 1000
plot11_jet.ymin            = 0.001
plot11_jet.ymax            = 100
#plot11_jet.lpos = "bottom-center"
plot11_jet.name            = "sTjet_allPreviousCuts"
plot11_jet.addZUncBand     = zUncBand
plot11_jet.histodata       = h_sT_DATA




# ##--- Mmj preselection

# h_Mmj_presel_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Mmj_1stPair_PAS", File_preselection).Clone()
# h_Mmj_presel_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Mmj_1stPair_PAS", File_preselection).Clone()
# h_Mmj_presel_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Mmj_1stPair_PAS", File_preselection).Clone()
# h_Mmj_presel_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Mmj_1stPair_PAS", File_preselection).Clone()
# h_Mmj_presel_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Mmj_1stPair_PAS", File_preselection).Clone()
h_Mmj_presel_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Mmj_1stPair_PAS", File_preselection).Clone()
h_Mmj_presel_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mmj_1stPair_PAS", File_preselection).Clone()
# h_Mmj_presel_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Mmj_1stPair_PAS", File_preselection).Clone()
# #h_Mmj_presel_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Mmj_1stPair_PAS", File_preselection).Clone()
# #h_Mmj_presel_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Mmj_1stPair_PAS", File_preselection).Clone()
# h_Mmj_presel_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Mmj_1stPair_PAS", File_preselection).Clone()
# h_Mmj_presel_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Mmj_1stPair_PAS", File_preselection).Clone()
# h_Mmj_presel_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Mmj_1stPair_PAS", File_preselection).Clone()
h_Mmj_presel_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mmj_1stPair_PAS", File_preselection).Clone()

# h_Mmj_presel_LQmmjj_M100.Add(GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Mmj_2ndPair_PAS", File_preselection))
# h_Mmj_presel_LQmmjj_M200.Add(GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Mmj_2ndPair_PAS", File_preselection))
# h_Mmj_presel_LQmmjj_M300.Add(GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Mmj_2ndPair_PAS", File_preselection))
# h_Mmj_presel_LQmmjj_M400.Add(GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Mmj_2ndPair_PAS", File_preselection))
# h_Mmj_presel_LQmmjj_M500.Add(GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Mmj_2ndPair_PAS", File_preselection))
h_Mmj_presel_TTbar.Add(GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Mmj_2ndPair_PAS", File_preselection))
h_Mmj_presel_ZJetAlpgen.Add(GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mmj_2ndPair_PAS", File_preselection))
# h_Mmj_presel_OTHERBKG.Add(GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Mmj_2ndPair_PAS", File_preselection))
# #h_Mmj_presel_QCD_Madgraph.Add(GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Mmj_2ndPair_PAS", File_preselection))
# #h_Mmj_presel_QCDPt15.Add(GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Mmj_2ndPair_PAS", File_preselection))
# h_Mmj_presel_SingleTop.Add(GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Mmj_2ndPair_PAS", File_preselection))
# h_Mmj_presel_VVjets.Add(GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Mmj_2ndPair_PAS", File_preselection))
# h_Mmj_presel_WJetAlpgen.Add(GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Mmj_2ndPair_PAS", File_preselection))
h_Mmj_presel_DATA.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mmj_2ndPair_PAS", File_preselection))

plot12 = Plot()
plot12.histosStack     = [h_Mmj_presel_TTbar, h_Mmj_presel_ZJetAlpgen]
plot12.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot12.histos          = [h_Mmj_presel_LQmmjj_M200, h_Mmj_presel_LQmmjj_M300]
# plot12.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot12.xtit            = "Mmj (GeV/c^{2})"
plot12.ytit            = "Number of events x 2"
plot12.ylog            = "yes"
plot12.rebin           = 2
plot12.xmin            = 0
plot12.xmax            = 1000
plot12.ymin            = 0.001
plot12.ymax            = 500
#plot12.lpos = "bottom-center"
plot12.name            = "Mmj_allPreviousCuts"
plot12.addZUncBand     = zUncBand
plot12.histodata       = h_Mmj_presel_DATA





# #--- Mmm_PAS (after preselection) ---

# h_Mmm_FullPreSel_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Mmm_PAS", File_preselection).Clone()
# h_Mmm_FullPreSel_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Mmm_PAS", File_preselection).Clone()
# h_Mmm_FullPreSel_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Mmm_PAS", File_preselection).Clone()
# h_Mmm_FullPreSel_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Mmm_PAS", File_preselection).Clone()
# h_Mmm_FullPreSel_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Mmm_PAS", File_preselection).Clone()
h_Mmm_FullPreSel_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Mmm_PAS", File_preselection).Clone()
h_Mmm_FullPreSel_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mmm_PAS", File_preselection).Clone()
# h_Mmm_FullPreSel_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Mmm_PAS", File_preselection).Clone()
# #h_Mmm_FullPreSel_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Mmm_PAS", File_preselection).Clone()
# #h_Mmm_FullPreSel_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Mmm_PAS", File_preselection).Clone()
# h_Mmm_FullPreSel_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Mmm_PAS", File_preselection).Clone()
# h_Mmm_FullPreSel_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Mmm_PAS", File_preselection).Clone()
# h_Mmm_FullPreSel_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Mmm_PAS", File_preselection).Clone()
h_Mmm_FullPreSel_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mmm_PAS", File_preselection).Clone()

plot13 = Plot()
## inputs for stacked histograms
## it created h_Mmm_FullPreSel_TTbar, h_Mmm_FullPreSel_TTbar+h_Mmm_FullPreSel_ZJetAlpgen , h_Mmm_FullPreSel_TTbar+h_Mmm_FullPreSel_ZJetAlpgen+h_Mmm_FullPreSel_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
plot13.histosStack     = [h_Mmm_FullPreSel_TTbar, h_Mmm_FullPreSel_ZJetAlpgen]
plot13.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot13.histos          = [h_Mmm_FullPreSel_LQmmjj_M200]
# plot13.keys            = ["LQ mmjj M200"]
plot13.xtit            = "M(mm) (GeV/c^{2})"
plot13.ytit            = "Number of events"
# plot13.ylog            = "yes"
# plot13.rebin           = 1
# plot13.ymin            = 0.00000001
# plot13.ymax            = 20
plot13.ylog            = "no"
plot13.rebin           = 1
plot13.ymin            = 0
plot13.ymax            = 40
plot13.xmin            = 0
plot13.xmax            = 200
#plot13.lpos = "bottom-center"
plot13.name            = "Mmm_FullPreSel_allPreviousCuts_ylin"
plot13.addZUncBand     = zUncBand
plot13.histodata       = h_Mmm_FullPreSel_DATA

plot13_ylog = Plot()
## inputs for stacked histograms
plot13_ylog.histosStack     = [h_Mmm_FullPreSel_TTbar, h_Mmm_FullPreSel_ZJetAlpgen]
plot13_ylog.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot13_ylog.histos          = [h_Mmm_FullPreSel_LQmmjj_M200, h_Mmm_FullPreSel_LQmmjj_M300]
# plot13_ylog.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot13_ylog.xtit            = "M(mm) (GeV/c^{2})"
plot13_ylog.ytit            = "Number of events"
plot13_ylog.ylog            = "yes"
plot13_ylog.rebin           = 1 # don't change it (since a rebinning is already applied above on the same histo)
plot13_ylog.ymin            = 0.001
plot13_ylog.ymax            = 100
plot13_ylog.xmin            = 0
plot13_ylog.xmax            = 1000
#plot13_ylog.lpos = "bottom-center"
plot13_ylog.name            = "Mmm_FullPreSel_allPreviousCuts"
plot13_ylog.addZUncBand     = zUncBand
plot13_ylog.histodata       = h_Mmm_FullPreSel_DATA



# #--- Mjj_PAS (after preselection) ---

# h_Mjj_FullPreSel_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
# h_Mjj_FullPreSel_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
# h_Mjj_FullPreSel_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
# h_Mjj_FullPreSel_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
# h_Mjj_FullPreSel_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
h_Mjj_FullPreSel_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
h_Mjj_FullPreSel_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
# h_Mjj_FullPreSel_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
# #h_Mjj_FullPreSel_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
# #h_Mjj_FullPreSel_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
# h_Mjj_FullPreSel_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
# h_Mjj_FullPreSel_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
# h_Mjj_FullPreSel_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()
h_Mjj_FullPreSel_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mjj_PAS", File_preselection).Clone()

plot14 = Plot()
## inputs for stacked histograms
## it created h_Mjj_FullPreSel_TTbar, h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen , h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen+h_Mjj_FullPreSel_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
plot14.histosStack     = [h_Mjj_FullPreSel_TTbar, h_Mjj_FullPreSel_ZJetAlpgen]
plot14.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]


## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot14.histos          = [h_Mjj_FullPreSel_LQmmjj_M200, h_Mjj_FullPreSel_LQmmjj_M300]
# plot14.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot14.xtit            = "M(jj) (GeV/c^{2})"
plot14.ytit            = "Number of events"
# plot14.ylog            = "yes"
# plot14.rebin           = 1
# plot14.ymin            = 0.00000001
# plot14.ymax            = 20
plot14.ylog            = "no"
plot14.rebin           = 1
plot14.ymin            = 0
plot14.ymax            = 10
plot14.xmin            = 0
plot14.xmax            = 1000
#plot14.lpos = "bottom-center"
plot14.name            = "Mjj_FullPreSel_allPreviousCuts_ylin"
plot14.addZUncBand     = zUncBand
plot14.histodata       = h_Mjj_FullPreSel_DATA

plot14_ylog = Plot()
## inputs for stacked histograms
## it created h_Mjj_FullPreSel_TTbar, h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen , h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen+h_Mjj_FullPreSel_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
plot14_ylog.histosStack     = [h_Mjj_FullPreSel_TTbar, h_Mjj_FullPreSel_ZJetAlpgen]
plot14_ylog.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot14_ylog.histos          = [h_Mjj_FullPreSel_LQmmjj_M200, h_Mjj_FullPreSel_LQmmjj_M300]
# plot14_ylog.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot14_ylog.xtit            = "M(jj) (GeV/c^{2})"
plot14_ylog.ytit            = "Number of events"
plot14_ylog.ylog            = "yes"
plot14_ylog.rebin           = 1 # don't change it (since a rebinning is already applied above on the same histo)
plot14_ylog.ymin            = 0.001
plot14_ylog.ymax            = 100
plot14_ylog.xmin            = 0
plot14_ylog.xmax            = 1000
#plot14_ylog.lpos = "bottom-center"
plot14_ylog.name            = "Mjj_FullPreSel_allPreviousCuts"
plot14_ylog.addZUncBand     = zUncBand
plot14_ylog.histodata       = h_Mjj_FullPreSel_DATA





# ##--- Pt Mus AllPreviousCuts ---

# h_pTMus_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()
# h_pTMus_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()
# h_pTMus_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()
# h_pTMus_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()
# h_pTMus_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()
h_pTMus_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()
h_pTMus_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()
# h_pTMus_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()
# #h_pTMus_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()
# #h_pTMus_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()
# h_pTMus_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()
# h_pTMus_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()
# h_pTMus_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()
h_pTMus_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt1stMu_PAS", File_preselection).Clone()

# h_pTMus_LQmmjj_M100.Add(GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection))
# h_pTMus_LQmmjj_M200.Add(GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection))
# h_pTMus_LQmmjj_M300.Add(GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection))
# h_pTMus_LQmmjj_M400.Add(GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection))
# h_pTMus_LQmmjj_M500.Add(GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection))
h_pTMus_TTbar.Add(GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection))
h_pTMus_ZJetAlpgen.Add(GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection))
# h_pTMus_OTHERBKG.Add(GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection))
# #h_pTMus_QCD_Madgraph.Add(GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection))
# #h_pTMus_QCDPt15.Add(GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection))
# h_pTMus_SingleTop.Add(GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection))
# h_pTMus_VVjets.Add(GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection))
# h_pTMus_WJetAlpgen.Add(GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection))
h_pTMus_DATA.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt2ndMu_PAS", File_preselection))

plot2and4 = Plot()
## inputs for stacked histograms
plot2and4.histosStack     = [h_pTMus_TTbar, h_pTMus_ZJetAlpgen]
plot2and4.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot2and4.histos          = [h_pTMus_LQmmjj_M200, h_pTMus_LQmmjj_M300]
# plot2and4.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot2and4.xtit            = "pT muons (GeV/c)"
plot2and4.ytit            = "Number of events x 2"
plot2and4.ylog            = "yes"
plot2and4.rebin           = 1
plot2and4.xmin            = pt_xmin
plot2and4.xmax            = pt_xmax
plot2and4.ymin            = pt_ymin
plot2and4.ymax            = pt_ymax
#plot2and4.lpos = "bottom-center"
plot2and4.name            = "pTMus_allPreviousCuts"
plot2and4.addZUncBand     = zUncBand
plot2and4.histodata       = h_pTMus_DATA

# ##--- Eta Mus AllPreviousCuts ---

# h_etaMus_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()
# h_etaMus_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()
# h_etaMus_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()
# h_etaMus_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()
# h_etaMus_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()
h_etaMus_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()
h_etaMus_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()
# h_etaMus_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()
# #h_etaMus_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()
# #h_etaMus_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()
# h_etaMus_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()
# h_etaMus_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()
# h_etaMus_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()
h_etaMus_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta1stMu_PAS", File_preselection).Clone()

# h_etaMus_LQmmjj_M100.Add(GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection))
# h_etaMus_LQmmjj_M200.Add(GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection))
# h_etaMus_LQmmjj_M300.Add(GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection))
# h_etaMus_LQmmjj_M400.Add(GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection))
# h_etaMus_LQmmjj_M500.Add(GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection))
h_etaMus_TTbar.Add(GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection))
h_etaMus_ZJetAlpgen.Add(GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection))
# h_etaMus_OTHERBKG.Add(GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection))
# #h_etaMus_QCD_Madgraph.Add(GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection))
# #h_etaMus_QCDPt15.Add(GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection))
# h_etaMus_SingleTop.Add(GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection))
# h_etaMus_VVjets.Add(GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection))
# h_etaMus_WJetAlpgen.Add(GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection))
h_etaMus_DATA.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta2ndMu_PAS", File_preselection))

plot3and5 = Plot()
## inputs for stacked histograms
plot3and5.histosStack     = [h_etaMus_TTbar, h_etaMus_ZJetAlpgen]
plot3and5.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot3and5.histos          = [h_etaMus_LQmmjj_M200, h_etaMus_LQmmjj_M300]
# plot3and5.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot3and5.xtit            = "#eta muons"
plot3and5.ytit            = "Number of events x 2"
plot3and5.rebin           = eta_rebin/2
plot3and5.ymin            = eta_ymin
plot3and5.ymax            = eta_ymax
plot3and5.lpos            = "top-left"
#plot3and5.lpos = "bottom-center"
plot3and5.name            = "etaMus_allPreviousCuts"
plot3and5.addZUncBand     = zUncBand
plot3and5.histodata       = h_etaMus_DATA


# ##--- Pt Jets AllPreviousCuts ---

# h_pTJets_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_pTJets_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_pTJets_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_pTJets_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_pTJets_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_pTJets_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_pTJets_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_pTJets_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# #h_pTJets_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# #h_pTJets_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_pTJets_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_pTJets_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
# h_pTJets_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()
h_pTJets_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt1stJet_PAS", File_preselection).Clone()

# h_pTJets_LQmmjj_M100.Add(GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
# h_pTJets_LQmmjj_M200.Add(GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
# h_pTJets_LQmmjj_M300.Add(GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
# h_pTJets_LQmmjj_M400.Add(GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
# h_pTJets_LQmmjj_M500.Add(GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
h_pTJets_TTbar.Add(GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
h_pTJets_ZJetAlpgen.Add(GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
# h_pTJets_OTHERBKG.Add(GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
# #h_pTJets_QCD_Madgraph.Add(GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
# #h_pTJets_QCDPt15.Add(GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
# h_pTJets_SingleTop.Add(GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
# h_pTJets_VVjets.Add(GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
# h_pTJets_WJetAlpgen.Add(GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))
h_pTJets_DATA.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt2ndJet_PAS", File_preselection))

plot7and9 = Plot()
## inputs for stacked histograms
plot7and9.histosStack     = [h_pTJets_TTbar, h_pTJets_ZJetAlpgen]
plot7and9.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot7and9.histos          = [h_pTJets_LQmmjj_M200, h_pTJets_LQmmjj_M300]
# plot7and9.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot7and9.xtit            = "pT jets (GeV/c)"
plot7and9.ytit            = "Number of events x 2"
plot7and9.ylog            = "yes"
plot7and9.rebin           = 1
plot7and9.xmin            = pt_xmin
plot7and9.xmax            = pt_xmax
plot7and9.ymin            = pt_ymin
plot7and9.ymax            = pt_ymax
#plot7and9.lpos = "bottom-center"
plot7and9.name            = "pTJets_allPreviousCuts"
plot7and9.addZUncBand     = zUncBand
plot7and9.histodata       = h_pTJets_DATA

# ##--- Eta Mus AllPreviousCuts ---

# h_etaJets_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_etaJets_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_etaJets_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_etaJets_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_etaJets_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_etaJets_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_etaJets_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_etaJets_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# #h_etaJets_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# #h_etaJets_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_etaJets_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_etaJets_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
# h_etaJets_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()
h_etaJets_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta1stJet_PAS", File_preselection).Clone()

# h_etaJets_LQmmjj_M100.Add(GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
# h_etaJets_LQmmjj_M200.Add(GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
# h_etaJets_LQmmjj_M300.Add(GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
# h_etaJets_LQmmjj_M400.Add(GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
# h_etaJets_LQmmjj_M500.Add(GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
h_etaJets_TTbar.Add(GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
h_etaJets_ZJetAlpgen.Add(GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
# h_etaJets_OTHERBKG.Add(GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
# #h_etaJets_QCD_Madgraph.Add(GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
# #h_etaJets_QCDPt15.Add(GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
# h_etaJets_SingleTop.Add(GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
# h_etaJets_VVjets.Add(GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
# h_etaJets_WJetAlpgen.Add(GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))
h_etaJets_DATA.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta2ndJet_PAS", File_preselection))

plot8and10 = Plot()
## inputs for stacked histograms
plot8and10.histosStack     = [h_etaJets_TTbar, h_etaJets_ZJetAlpgen]
plot8and10.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot8and10.histos          = [h_etaJets_LQmmjj_M200, h_etaJets_LQmmjj_M300]
# plot8and10.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot8and10.xtit            = "#eta jets"
plot8and10.ytit            = "Number of events x 2"
plot8and10.rebin           = eta_rebin/2
plot8and10.ymin            = eta_ymin
plot8and10.ymax            = eta_ymax
plot8and10.lpos            = "top-left"
#plot8and10.lpos = "bottom-center"
plot8and10.name            = "etaJets_allPreviousCuts"
plot8and10.addZUncBand     = zUncBand
plot8and10.histodata       = h_etaJets_DATA


# #--- pfMET ---

# h_pfMET_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________pfMET_PAS", File_preselection).Clone()
# h_pfMET_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________pfMET_PAS", File_preselection).Clone()
# h_pfMET_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________pfMET_PAS", File_preselection).Clone()
# h_pfMET_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________pfMET_PAS", File_preselection).Clone()
# h_pfMET_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________pfMET_PAS", File_preselection).Clone()
h_pfMET_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________pfMET_PAS", File_preselection).Clone()
h_pfMET_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________pfMET_PAS", File_preselection).Clone()
# h_pfMET_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________pfMET_PAS", File_preselection).Clone()
# #h_pfMET_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________pfMET_PAS", File_preselection).Clone()
# #h_pfMET_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________pfMET_PAS", File_preselection).Clone()
# h_pfMET_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________pfMET_PAS", File_preselection).Clone()
# h_pfMET_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________pfMET_PAS", File_preselection).Clone()
# h_pfMET_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________pfMET_PAS", File_preselection).Clone()
h_pfMET_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________pfMET_PAS", File_preselection).Clone()


plot15 = Plot()
## inputs for stacked histograms
## it created h_pfMET_TTbar, h_pfMET_TTbar+h_pfMET_ZJetAlpgen , h_pfMET_TTbar+h_pfMET_ZJetAlpgen+h_pfMET_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
plot15.histosStack     = [h_pfMET_TTbar, h_pfMET_ZJetAlpgen]
plot15.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot15.histos          = [h_pfMET_LQmmjj_M200, h_pfMET_LQmmjj_M300, h_pfMET_LQmmjj_M400, h_pfMET_LQmmjj_M500]
# plot15.keys            = ["LQ mmjj M200","LQ mmjj M300", "LQ mmjj M400", "LQ mmjj M500"]
plot15.xtit            = "pfMET (GeV/c)"
plot15.ytit            = "Number of events"
#plot15.xlog            = "yes"
plot15.ylog            = "yes"
plot15.rebin           = 1
plot15.xmin            = 0
plot15.xmax            = 300
plot15.ymin            = 0.001
plot15.ymax            = 400
#plot15.lpos = "bottom-center"
plot15.name            = "pfMET_allPreviousCuts"
plot15.addZUncBand     = zUncBand
plot15.histodata       = h_pfMET_DATA


# #--- caloMET ---

# h_caloMET_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________caloMET_PAS", File_preselection).Clone()
# h_caloMET_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________caloMET_PAS", File_preselection).Clone()
# h_caloMET_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________caloMET_PAS", File_preselection).Clone()
# h_caloMET_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________caloMET_PAS", File_preselection).Clone()
# h_caloMET_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________caloMET_PAS", File_preselection).Clone()
h_caloMET_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________caloMET_PAS", File_preselection).Clone()
h_caloMET_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________caloMET_PAS", File_preselection).Clone()
# h_caloMET_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________caloMET_PAS", File_preselection).Clone()
# #h_caloMET_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________caloMET_PAS", File_preselection).Clone()
# #h_caloMET_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________caloMET_PAS", File_preselection).Clone()
# h_caloMET_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________caloMET_PAS", File_preselection).Clone()
# h_caloMET_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________caloMET_PAS", File_preselection).Clone()
# h_caloMET_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________caloMET_PAS", File_preselection).Clone()
h_caloMET_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________caloMET_PAS", File_preselection).Clone()


plot16 = Plot()
## inputs for stacked histograms
## it created h_caloMET_TTbar, h_caloMET_TTbar+h_caloMET_ZJetAlpgen , h_caloMET_TTbar+h_caloMET_ZJetAlpgen+h_caloMET_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
plot16.histosStack     = [h_caloMET_TTbar, h_caloMET_ZJetAlpgen]
plot16.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot16.histos          = [h_caloMET_LQmmjj_M200, h_caloMET_LQmmjj_M300, h_caloMET_LQmmjj_M400, h_caloMET_LQmmjj_M500]
# plot16.keys            = ["LQ mmjj M200","LQ mmjj M300", "LQ mmjj M400", "LQ mmjj M500"]
plot16.xtit            = "caloMET (GeV/c)"
plot16.ytit            = "Number of events"
#plot16.xlog            = "yes"
plot16.ylog            = "yes"
plot16.rebin           = 1
plot16.xmin            = 0
plot16.xmax            = 300
plot16.ymin            = 0.001
plot16.ymax            = 500
#plot16.lpos = "bottom-center"
plot16.name            = "caloMET_allPreviousCuts"
plot16.addZUncBand     = zUncBand
plot16.histodata       = h_caloMET_DATA




# #--- Pt1stMu_IDISO  ---

# h_pT1stMu_nojet_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Pt1stMu_IDISO", File_preselection).Clone()
# h_pT1stMu_nojet_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Pt1stMu_IDISO", File_preselection).Clone()
# h_pT1stMu_nojet_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Pt1stMu_IDISO", File_preselection).Clone()
# h_pT1stMu_nojet_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Pt1stMu_IDISO", File_preselection).Clone()
# h_pT1stMu_nojet_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Pt1stMu_IDISO", File_preselection).Clone()
h_pT1stMu_nojet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Pt1stMu_IDISO", File_preselection).Clone()
h_pT1stMu_nojet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt1stMu_IDISO", File_preselection).Clone()
# h_pT1stMu_nojet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Pt1stMu_IDISO", File_preselection).Clone()
# #h_pT1stMu_nojet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt1stMu_IDISO", File_preselection).Clone()
# #h_pT1stMu_nojet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt1stMu_IDISO", File_preselection).Clone()
# h_pT1stMu_nojet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt1stMu_IDISO", File_preselection).Clone()
# h_pT1stMu_nojet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt1stMu_IDISO", File_preselection).Clone()
# h_pT1stMu_nojet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt1stMu_IDISO", File_preselection).Clone()
h_pT1stMu_nojet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt1stMu_IDISO", File_preselection).Clone()

plot2_nojet = Plot()
## inputs for stacked histograms
plot2_nojet.histosStack     = [h_pT1stMu_nojet_TTbar, h_pT1stMu_nojet_ZJetAlpgen]
plot2_nojet.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot2_nojet.histos          = [h_pT1stMu_nojet_LQmmjj_M200, h_pT1stMu_nojet_LQmmjj_M300]
# plot2_nojet.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot2_nojet.xtit            = "pT 1st muon (GeV/c)"
plot2_nojet.ytit            = "Number of events"
plot2_nojet.ylog            = "yes"
plot2_nojet.rebin           = 1
plot2_nojet.xmin            = pt_xmin
plot2_nojet.xmax            = pt_xmax
plot2_nojet.ymin            = pt_ymin
plot2_nojet.ymax            = pt_ymax*10
#plot2_nojet.lpos = "bottom-center"
plot2_nojet.name            = "pT1stMu_nojet_allPreviousCuts"
plot2_nojet.addZUncBand     = zUncBand
plot2_nojet.histodata       = h_pT1stMu_nojet_DATA


# #--- Eta1stMu_IDISO  ---

# h_Eta1stMu_nojet_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Eta1stMu_IDISO", File_preselection).Clone()
# h_Eta1stMu_nojet_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Eta1stMu_IDISO", File_preselection).Clone()
# h_Eta1stMu_nojet_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Eta1stMu_IDISO", File_preselection).Clone()
# h_Eta1stMu_nojet_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Eta1stMu_IDISO", File_preselection).Clone()
# h_Eta1stMu_nojet_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Eta1stMu_IDISO", File_preselection).Clone()
h_Eta1stMu_nojet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Eta1stMu_IDISO", File_preselection).Clone()
h_Eta1stMu_nojet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta1stMu_IDISO", File_preselection).Clone()
# h_Eta1stMu_nojet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Eta1stMu_IDISO", File_preselection).Clone()
# #h_Eta1stMu_nojet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta1stMu_IDISO", File_preselection).Clone()
# #h_Eta1stMu_nojet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta1stMu_IDISO", File_preselection).Clone()
# h_Eta1stMu_nojet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta1stMu_IDISO", File_preselection).Clone()
# h_Eta1stMu_nojet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta1stMu_IDISO", File_preselection).Clone()
# h_Eta1stMu_nojet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta1stMu_IDISO", File_preselection).Clone()
h_Eta1stMu_nojet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta1stMu_IDISO", File_preselection).Clone()

plot3_nojet = Plot()
## inputs for stacked histograms
plot3_nojet.histosStack     = [h_Eta1stMu_nojet_TTbar, h_Eta1stMu_nojet_ZJetAlpgen]
plot3_nojet.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot3_nojet.histos          = [h_Eta1stMu_nojet_LQmmjj_M200, h_Eta1stMu_nojet_LQmmjj_M300]
# plot3_nojet.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot3_nojet.xtit            = "#eta 1st muon"
plot3_nojet.ytit            = "Number of events"
plot3_nojet.rebin           = eta_rebin
plot3_nojet.ymin            = eta_ymin
plot3_nojet.ymax            = eta_ymax*10
plot3_nojet.lpos = "top-left"
plot3_nojet.name            = "Eta1stMu_nojet_allPreviousCuts"
plot3_nojet.addZUncBand     = zUncBand
plot3_nojet.histodata       = h_Eta1stMu_nojet_DATA



# #--- Pt2ndMu_IDISO  ---

# h_pT2ndMu_nojet_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Pt2ndMu_IDISO", File_preselection).Clone()
# h_pT2ndMu_nojet_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Pt2ndMu_IDISO", File_preselection).Clone()
# h_pT2ndMu_nojet_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Pt2ndMu_IDISO", File_preselection).Clone()
# h_pT2ndMu_nojet_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Pt2ndMu_IDISO", File_preselection).Clone()
# h_pT2ndMu_nojet_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Pt2ndMu_IDISO", File_preselection).Clone()
h_pT2ndMu_nojet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Pt2ndMu_IDISO", File_preselection).Clone()
h_pT2ndMu_nojet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Pt2ndMu_IDISO", File_preselection).Clone()
# h_pT2ndMu_nojet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Pt2ndMu_IDISO", File_preselection).Clone()
# #h_pT2ndMu_nojet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Pt2ndMu_IDISO", File_preselection).Clone()
# #h_pT2ndMu_nojet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Pt2ndMu_IDISO", File_preselection).Clone()
# h_pT2ndMu_nojet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Pt2ndMu_IDISO", File_preselection).Clone()
# h_pT2ndMu_nojet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Pt2ndMu_IDISO", File_preselection).Clone()
# h_pT2ndMu_nojet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Pt2ndMu_IDISO", File_preselection).Clone()
h_pT2ndMu_nojet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Pt2ndMu_IDISO", File_preselection).Clone()

plot4_nojet = Plot()
## inputs for stacked histograms
plot4_nojet.histosStack     = [h_pT2ndMu_nojet_TTbar, h_pT2ndMu_nojet_ZJetAlpgen]
plot4_nojet.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot4_nojet.histos          = [h_pT2ndMu_nojet_LQmmjj_M200, h_pT2ndMu_nojet_LQmmjj_M300]
# plot4_nojet.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot4_nojet.xtit            = "pT 2nd muon (GeV/c)"
plot4_nojet.ytit            = "Number of events"
plot4_nojet.ylog            = "yes"
plot4_nojet.rebin           = 1
plot4_nojet.xmin            = pt_xmin
plot4_nojet.xmax            = pt_xmax
plot4_nojet.ymin            = pt_ymin
plot4_nojet.ymax            = pt_ymax*10
#plot4_nojet.lpos = "bottom-center"
plot4_nojet.name            = "pT2ndMu_nojet_allPreviousCuts"
plot4_nojet.addZUncBand     = zUncBand
plot4_nojet.histodata       = h_pT2ndMu_nojet_DATA


# #--- Eta2ndMu_IDISO  ---

# h_Eta2ndMu_nojet_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Eta2ndMu_IDISO", File_preselection).Clone()
# h_Eta2ndMu_nojet_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Eta2ndMu_IDISO", File_preselection).Clone()
# h_Eta2ndMu_nojet_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Eta2ndMu_IDISO", File_preselection).Clone()
# h_Eta2ndMu_nojet_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Eta2ndMu_IDISO", File_preselection).Clone()
# h_Eta2ndMu_nojet_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Eta2ndMu_IDISO", File_preselection).Clone()
h_Eta2ndMu_nojet_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Eta2ndMu_IDISO", File_preselection).Clone()
h_Eta2ndMu_nojet_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Eta2ndMu_IDISO", File_preselection).Clone()
# h_Eta2ndMu_nojet_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Eta2ndMu_IDISO", File_preselection).Clone()
# #h_Eta2ndMu_nojet_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Eta2ndMu_IDISO", File_preselection).Clone()
# #h_Eta2ndMu_nojet_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Eta2ndMu_IDISO", File_preselection).Clone()
# h_Eta2ndMu_nojet_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Eta2ndMu_IDISO", File_preselection).Clone()
# h_Eta2ndMu_nojet_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Eta2ndMu_IDISO", File_preselection).Clone()
# h_Eta2ndMu_nojet_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Eta2ndMu_IDISO", File_preselection).Clone()
h_Eta2ndMu_nojet_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Eta2ndMu_IDISO", File_preselection).Clone()

plot5_nojet = Plot()
## inputs for stacked histograms
plot5_nojet.histosStack     = [h_Eta2ndMu_nojet_TTbar, h_Eta2ndMu_nojet_ZJetAlpgen]
plot5_nojet.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot5_nojet.histos          = [h_Eta2ndMu_nojet_LQmmjj_M200, h_Eta2ndMu_nojet_LQmmjj_M300]
# plot5_nojet.keys            = ["LQ mmjj M200","LQ mmjj M300"]
plot5_nojet.xtit            = "#eta 2nd muon"
plot5_nojet.ytit            = "Number of events"
plot5_nojet.rebin           = eta_rebin
plot5_nojet.ymin            = eta_ymin
plot5_nojet.ymax            = eta_ymax*10
plot5_nojet.lpos = "top-left"
plot5_nojet.name            = "Eta2ndMu_nojet_allPreviousCuts"
plot5_nojet.addZUncBand     = zUncBand
plot5_nojet.histodata       = h_Eta2ndMu_nojet_DATA



# # ############################ Plots below to be done after full selection ######################

# ##--- sT AllOtherCuts ---

# plot20 = Plot()
# plot20.histosStack     = [
#     GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
#     GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
#     GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
# #     GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
# #     GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
# #     GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
# #     GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________sT", File_selection).Clone()
#     ]
# plot20.keysStack       = [
#     "ttbar", "Z/#gamma/Z* + jets"]
# ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot20.histos          = [
#     GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________sT", File_selection).Clone(),
#     GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________sT", File_selection).Clone()
#     ]
# plot20.keys            = ["LQ mmjj M200","LQ mmjj M300"]
# plot20.xtit            = "St (GeV/c)"
# plot20.ytit            = "Number of events"
# plot20.ylog            = "yes"
# plot20.rebin           = 10
# plot20.xmin            = 0
# plot20.xmax            = 1000
# plot20.ymin            = 0.001
# plot20.ymax            = 100
# #plot20.lpos = "bottom-center"
# plot20.name            = "sT_allOtherCuts"
# plot20.addZUncBand     = zUncBand
# plot20.histodata       = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________sT", File_selection).Clone()


# ##--- Mmj AllOtherCuts ---

# h_Mmj_LQmmjj_M100 = GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Mmj_1stPair", File_selection).Clone()
# h_Mmj_LQmmjj_M200 = GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Mmj_1stPair", File_selection).Clone()
# h_Mmj_LQmmjj_M300 = GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Mmj_1stPair", File_selection).Clone()
# h_Mmj_LQmmjj_M400 = GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Mmj_1stPair", File_selection).Clone()
# h_Mmj_LQmmjj_M500 = GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Mmj_1stPair", File_selection).Clone()
# h_Mmj_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Mmj_1stPair", File_selection).Clone()
# h_Mmj_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mmj_1stPair", File_selection).Clone()
# h_Mmj_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Mmj_1stPair", File_selection).Clone()
# #h_Mmj_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Mmj_1stPair", File_selection).Clone()
# #h_Mmj_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Mmj_1stPair", File_selection).Clone()
# h_Mmj_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Mmj_1stPair", File_selection).Clone()
# h_Mmj_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Mmj_1stPair", File_selection).Clone()
# h_Mmj_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Mmj_1stPair", File_selection).Clone()
# h_Mmj_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mmj_1stPair", File_selection).Clone()

# h_Mmj_LQmmjj_M100.Add(GetHisto("histo1D__LQmmjj_M100__cutHisto_allPreviousCuts________Mmj_2ndPair", File_selection))
# h_Mmj_LQmmjj_M200.Add(GetHisto("histo1D__LQmmjj_M200__cutHisto_allPreviousCuts________Mmj_2ndPair", File_selection))
# h_Mmj_LQmmjj_M300.Add(GetHisto("histo1D__LQmmjj_M300__cutHisto_allPreviousCuts________Mmj_2ndPair", File_selection))
# h_Mmj_LQmmjj_M400.Add(GetHisto("histo1D__LQmmjj_M400__cutHisto_allPreviousCuts________Mmj_2ndPair", File_selection))
# h_Mmj_LQmmjj_M500.Add(GetHisto("histo1D__LQmmjj_M500__cutHisto_allPreviousCuts________Mmj_2ndPair", File_selection))
# h_Mmj_TTbar.Add(GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Mmj_2ndPair", File_selection))
# h_Mmj_ZJetAlpgen.Add(GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mmj_2ndPair", File_selection))
# h_Mmj_OTHERBKG.Add(GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Mmj_2ndPair", File_selection))
# #h_Mmj_QCD_Madgraph.Add(GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Mmj_2ndPair", File_selection))
# #h_Mmj_QCDPt15.Add(GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Mmj_2ndPair", File_selection))
# h_Mmj_SingleTop.Add(GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Mmj_2ndPair", File_selection))
# h_Mmj_VVjets.Add(GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Mmj_2ndPair", File_selection))
# h_Mmj_WJetAlpgen.Add(GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Mmj_2ndPair", File_selection))
# h_Mmj_DATA.Add(GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mmj_2ndPair", File_selection))

# plot21 = Plot()
# plot21.histosStack     = [h_Mmj_TTbar, h_Mmj_ZJetAlpgen, h_Mmj_OTHERBKG]
# plot21.keysStack       = ["ttbar", "Z/#gamma/Z* + jets"]
# ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot21.histos          = [h_Mmj_LQmmjj_M200, h_Mmj_LQmmjj_M300]
# plot21.keys            = ["LQ mmjj M200","LQ mmjj M300"]
# plot21.xtit            = "Mmj (GeV/c^{2})"
# plot21.ytit            = "Number of events x 2"
# plot21.ylog            = "yes"
# plot21.rebin           = 1
# plot21.xmin            = 0
# plot21.xmax            = 600
# plot21.ymin            = 0.001
# plot21.ymax            = 20
# #plot21.lpos = "bottom-center"
# plot21.name            = "Mmj_allOtherCuts"
# plot21.addZUncBand     = zUncBand
# plot21.histodata       = h_Mmj_DATA


#-----------------------------------------------------------------------------------


# List of plots to be plotted
plots = [plot0, plot0_ylog, plot2_nojet, plot3_nojet, plot4_nojet, plot5_nojet, plot1, plot2, plot3, plot4, plot5,
         plot2and4, plot3and5, plot6, plot7, plot8,
         plot9, plot10, plot7and9, plot8and10, plot11, plot11_ele, plot11_jet, plot12, plot13, plot13_ylog, plot14, plot14_ylog,
         plot15, plot16]  # produced using preselection root file
#          plot20, plot21] # produced using full selection root file



############# USER CODE - END ################################################
##############################################################################


#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()
c.Print("allPlots.ps[")
for plot in plots:
    plot.Draw("allPlots.ps")
c.Print("allPlots.ps]")
