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

def generateHistoList( histoBaseName , samples , variableName, fileName ):
    histolist = []
    for sample in samples:
        hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
        histolist.append(GetHisto(hname, fileName).Clone())
    return histolist
                                                    
def generateAndAddHistoList( histoBaseName , samples , variableNames, fileName ):
    histolist = []
    for sample in samples:
        iv=0
        for variableName in variableNames:
            hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
            if (iv==0):
                histo = GetHisto(hname, fileName).Clone()
            else:
                histo.Add(GetHisto(hname, fileName).Clone())
            iv=iv+1
        histolist.append(histo)
    return histolist

def generateHisto( histoBaseName , sample , variableName, fileName ):
    hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
    histo = GetHisto(hname, fileName).Clone()
    return histo

def generateAndAddHisto( histoBaseName , sample , variableNames, fileName ):
    iv=0
    for variableName in variableNames:
        hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
        if (iv==0):
            histo = GetHisto(hname, fileName).Clone()
        else:
            histo.Add(GetHisto(hname, fileName).Clone())
        iv=iv+1
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
            legend.AddEntry(plot.histodata, "True ecj","p")
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

File_true = GetFile("$LQDATA/qcd/2.47pb-1/output_cutTable_ecj_ClosureTest_true/analysisClass_ecj_ClosureTest_plots.root")
File_pred = GetFile("$LQDATA/qcd/2.47pb-1/output_cutTable_ecj_ClosureTest_predicted/analysisClass_ecj_ClosureTest_plots.root")

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

histoBaseName = "histo1D__SAMPLE__cutHisto_allPreviousCuts________VARIABLE"

samplesForStackHistos = ["DATA"]
keysStack =             ["Predicted ecj"]

samplesForHistos = []
keys             = []

sampleForDataHisto = "DATA"


#--- Mee ---
variableName = "Mee"
plot1 = Plot()
## inputs for stacked histograms
plot1.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_pred)
plot1.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
#plot1.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_pred)
plot1.keys            = keys
plot1.xtit            = "Mec (GeV/c^{2})"
plot1.ytit            = "Number of events"
#plot1.xlog            = "yes"
plot1.ylog            = "yes"
plot1.rebin           = 2
plot1.xmin            = 100
plot1.xmax            = 1000
plot1.ymin            = 0.01
plot1.ymax            = 100
#plot1.lpos = "bottom-center"
plot1.name            = "Mec"
plot1.addZUncBand     = zUncBand
plot1.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_true)

#--- sT ---
variableName = "sT"
plot2 = Plot()
## inputs for stacked histograms
plot2.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_pred)
plot2.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
#plot2.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_pred)
plot2.keys            = keys
plot2.xtit            = "St (GeV/c)"
plot2.ytit            = "Number of events"
#plot2.xlog            = "yes"
plot2.ylog            = "yes"
plot2.rebin           = 2
plot2.xmin            = 50
plot2.xmax            = 1000
plot2.ymin            = 0.01
plot2.ymax            = 100
#plot2.lpos = "bottom-center"
plot2.name            = "sT"
plot2.addZUncBand     = zUncBand
plot2.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_true)

# #--- Pt1stEle_PAS ---
# variableName = "Pt1stEle_PAS"
# plot2 = Plot()
# ## inputs for stacked histograms
# plot2.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_pred)
# plot2.keysStack       = keysStack
# ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# #plot2.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_pred)
# plot2.keys            = keys
# plot2.xtit            = "Pt leading supercluster (GeV/c)"
# plot2.ytit            = "Number of events"
# #plot2.xlog            = "yes"
# plot2.ylog            = "yes"
# plot2.rebin           = 2
# plot2.xmin            = 0
# plot2.xmax            = 500
# plot2.ymin            = 0.01
# plot2.ymax            = 10000
# #plot2.lpos = "bottom-center"
# plot2.name            = "Pt1stEle_PAS"
# plot2.addZUncBand     = zUncBand
# plot2.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_true)

# #--- Pt2ndEle_PAS ---
# variableName = "Pt2ndEle_PAS"
# plot3 = Plot()
# ## inputs for stacked histograms
# plot3.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_pred)
# plot3.keysStack       = keysStack
# ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# #plot3.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_pred)
# plot3.keys            = keys
# plot3.xtit            = "Pt second leading supercluster (GeV/c)"
# plot3.ytit            = "Number of events"
# #plot3.xlog            = "yes"
# plot3.ylog            = "yes"
# plot3.rebin           = 2
# plot3.xmin            = 0
# plot3.xmax            = 500
# plot3.ymin            = 0.01
# plot3.ymax            = 10000
# #plot3.lpos = "bottom-center"
# plot3.name            = "Pt2ndEle_PAS"
# plot3.addZUncBand     = zUncBand
# plot3.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_true)

# #--- Pt1stJet_PAS ---
# variableName = "Pt1stJet_PAS"
# plot4 = Plot()
# ## inputs for stacked histograms
# plot4.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_pred)
# plot4.keysStack       = keysStack
# ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# #plot4.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_pred)
# plot4.keys            = keys
# plot4.xtit            = "Pt leading supercluster (GeV/c)"
# plot4.ytit            = "Number of events"
# #plot4.xlog            = "yes"
# plot4.ylog            = "yes"
# plot4.rebin           = 2
# plot4.xmin            = 0
# plot4.xmax            = 500
# plot4.ymin            = 0.01
# plot4.ymax            = 10000
# #plot4.lpos = "bottom-center"
# plot4.name            = "Pt1stJet_PAS"
# plot4.addZUncBand     = zUncBand
# plot4.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_true)




#-----------------------------------------------------------------------------------


# List of plots to be plotted
plots = [plot1, plot2] 



############# USER CODE - END ################################################
##############################################################################


#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()
c.Print("allPlots.ps[")
for plot in plots:
    plot.Draw("allPlots.ps")
c.Print("allPlots.ps]")
