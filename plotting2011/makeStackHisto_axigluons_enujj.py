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
gStyle.SetPadTopMargin(0.08);
gStyle.SetPadBottomMargin(0.12);
#gStyle.SetTitleSize(0.05, "XYZ");
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
    if( not histo ):
        print "ERROR: histo " + histoName + " not found in " + file.GetName()
        print "exiting..."
        sys.exit()
    new = copy.deepcopy(histo)
    if(scale!=1):
        new.Scale(scale)
    return new

def generateHistoList( histoBaseName , samples , variableName, fileName , scale = 1):
    histolist = []
    for sample in samples:
        hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
        histolist.append(GetHisto(hname, fileName, scale))
    return histolist

def generateAndAddHistoList( histoBaseName , samples , variableNames, fileName , scale = 1):
    histolist = []
    for sample in samples:
        iv=0
        histo = TH1F()
        for variableName in variableNames:
            hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
            if (iv==0):
                histo = GetHisto(hname, fileName, scale)
            else:
                histo.Add(GetHisto(hname, fileName, scale))
            iv=iv+1
        histolist.append(histo)
    return histolist

def generateHisto( histoBaseName , sample , variableName, fileName , scale = 1):
    hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
    histo = GetHisto(hname, fileName)
    return histo

def generateAndAddHisto( histoBaseName , sample , variableNames, fileName , scale = 1):
    iv=0
    histo = TH1F()
    for variableName in variableNames:
        hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
        if (iv==0):
            histo = GetHisto(hname, fileName, scale)
        else:
            histo.Add(GetHisto(hname, fileName, scale))
        iv=iv+1
    return histo

def rebinHisto( histo, xmin, xmax, rebin, xbins, addOvfl ):
    new_histo = TH1F()
    minBinWidth = 0
    if( xmin!="" and xmax!="" and rebin!="var" ):
        if(rebin!=""):
            histo.Rebin(rebin)
        minBinWidth = histo.GetBinWidth(1)
        xbinmin = histo.GetXaxis().FindBin(xmin)
        xbinmax = histo.GetXaxis().FindBin(xmax-0.000001)
        underflowBinContent = 0
        underflowBinError2 = 0
        for iter in range(0,xbinmin):
            underflowBinContent = underflowBinContent + histo.GetBinContent(iter)
            underflowBinError2 = underflowBinError2 + histo.GetBinError(iter)**2
        overflowBinContent = 0
        overflowBinError2 = 0
        for iter in range(xbinmax+1,histo.GetNbinsX()+2):
            overflowBinContent = overflowBinContent + histo.GetBinContent(iter)
            overflowBinError2 = overflowBinError2 + histo.GetBinError(iter)**2
        nbins = (xbinmax-xbinmin+1)
        xmin = histo.GetXaxis().GetBinLowEdge(xbinmin)
        xmax = histo.GetXaxis().GetBinUpEdge(xbinmax)
        new_histo.SetBins( nbins, xmin, xmax )
        for iter in range(1,nbins+1):
            new_histo.SetBinContent( iter, histo.GetBinContent(xbinmin+iter-1) )
            new_histo.SetBinError( iter, histo.GetBinError(xbinmin+iter-1) )
        new_histo.SetBinContent( 0, underflowBinContent )
        new_histo.SetBinError( 0, sqrt(underflowBinError2) )
        if( addOvfl=="yes"):
            new_histo.SetBinContent( nbins, new_histo.GetBinContent(nbins) + overflowBinContent )
            new_histo.SetBinError( nbins, sqrt(new_histo.GetBinError(nbins)**2 + overflowBinError2) )
    elif( xbins!="" and rebin=="var" ):
        binWidths = []
        for iter in range(0,len(xbins)-1):
            binWidths.append(float(xbins[iter+1] - xbins[iter]))
        minBinWidth = min(binWidths)
        xbinmin = histo.GetXaxis().FindBin(xbins[0])
        xbinmax = histo.GetXaxis().FindBin(xbins[-1]-0.000001)
        underflowBinContent = 0
        underflowBinError2 = 0
        for iter in range(0,xbinmin):
            underflowBinContent = underflowBinContent + histo.GetBinContent(iter)
            underflowBinError2 = underflowBinError2 + histo.GetBinError(iter)**2
        overflowBinContent = 0
        overflowBinError2 = 0
        for iter in range(xbinmax+1,histo.GetNbinsX()+2):
            overflowBinContent = overflowBinContent + histo.GetBinContent(iter)
            overflowBinError2 = overflowBinError2 + histo.GetBinError(iter)**2
        xbins[0] = histo.GetXaxis().GetBinLowEdge(xbinmin)
        xbins[-1] = histo.GetXaxis().GetBinUpEdge(xbinmax)
        xbinsFinal = array( 'd', xbins )
        nbins = len(xbinsFinal)-1
        new_histo = histo.Rebin( nbins , "new_histo", xbinsFinal )
        new_histo.SetBinContent( 0, underflowBinContent )
        new_histo.SetBinError( 0, sqrt(underflowBinError2) )
        for iter in range(1,nbins):
            new_histo.SetBinContent( iter, new_histo.GetBinContent(iter)*(minBinWidth/binWidths[iter-1]) )
            new_histo.SetBinError( iter, new_histo.GetBinError(iter)*(minBinWidth/binWidths[iter-1]) )
        if( addOvfl=="yes"):
            new_histo.SetBinContent( nbins, (new_histo.GetBinContent(nbins) + overflowBinContent)*(minBinWidth/binWidths[nbins-1]) )
            new_histo.SetBinError( nbins, sqrt(new_histo.GetBinError(nbins)**2 + overflowBinError2)*(minBinWidth/binWidths[nbins-1]) )
        else:
            new_histo.SetBinContent( nbins, new_histo.GetBinContent(nbins)*(minBinWidth/binWidths[nbins-1]) )
            new_histo.SetBinError( nbins, new_histo.GetBinError(nbins)*(minBinWidth/binWidths[nbins-1]) )
    else:
        if(rebin!=""):
            histo.Rebin(rebin)
        new_histo = histo
        minBinWidth = histo.GetBinWidth(1)
    return [new_histo, minBinWidth]

def rebinHistos( histos, xmin, xmax, rebin, xbins, addOvfl ):
    new_histos = []
    for histo in histos:
        new_histo = TH1F()
        new_histo = rebinHisto( histo, xmin, xmax, rebin, xbins, addOvfl )[0]
        new_histos.append(new_histo)
    return new_histos


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
    addOvfl     = "yes" # add the overflow bin to the last visible bin (default = "yes", option="no")
    name        = "" # name of the final plots
    lint        = "330 pb^{-1}" # integrated luminosity of the sample ( example "10 pb^{-1}" )
    addZUncBand = "no" # add an uncertainty band coming from the data-MC Z+jets rescaling (default = "no", option="yes")
    ZUncKey     = "Z/#gamma/Z* + jets unc." # key to be put in the legend for the Z+jets uncertainty band
    ZPlotIndex  = 1 # index of the Z+jets plots in the histosStack list (default = 1)
    ZScaleUnc   = 0.20 # uncertainty of the data-MC Z+jets scale factor
    makeRatio   = "" # 1=simple ratio, 2=ratio of cumulative histograms
    xbins       = "" #array with variable bin structure
    histodata   = "" # data histogram

    def Draw(self, fileps):

        self.histos = rebinHistos( self.histos, self.xmin, self.xmax, self.rebin, self.xbins, self.addOvfl )
        self.histosStack = rebinHistos( self.histosStack, self.xmin, self.xmax, self.rebin, self.xbins, self.addOvfl )
        resultArray = rebinHisto( self.histodata, self.xmin, self.xmax, self.rebin, self.xbins, self.addOvfl )
        self.histodata = resultArray[0]
        minBinW = resultArray[1]

        #-- create canvas
        canvas = TCanvas()
        stack = {}

        if(self.makeRatio==1):
            fPads1 = TPad("pad1", "", 0.00, 0.20, 0.99, 0.99)
            fPads2 = TPad("pad2", "", 0.00, 0.00, 0.99, 0.20)
            fPads1.SetFillColor(0)
            fPads1.SetLineColor(0)
            fPads2.SetFillColor(0)
            fPads2.SetLineColor(0)
            fPads1.Draw()
            fPads2.Draw()
        else:
            fPads1 = TPad("pad1", "", 0.00, 0.0, 0.99, 0.99)
            fPads1.SetFillColor(0)
            fPads1.SetLineColor(0)
            fPads1.Draw()

        #-- 1st pad
        fPads1.cd()
        #-- log scale
        # xlog may not work if (self.xlog == "yes"):
        # fPads1.SetLogx()
        if (self.ylog     == "yes"):
            fPads1.SetLogy()

        #-- legend
#        hsize=0.22
#        vsize=0.26
        hsize=0.35
        vsize=0.35
        if (self.lpos=="bottom-center"):
            xstart=0.35
            ystart=0.25
        elif(self.lpos=="top-left"):
            xstart=0.12
#            ystart=0.63
            ystart=0.54
        else:
            xstart=0.55
            ystart=0.52
#            xstart=0.65
#            ystart=0.63
        legend = TLegend(xstart, ystart, xstart+hsize, ystart+vsize)
        legend.SetFillColor(kWhite)
        legend.SetBorderSize(0)
        legend.SetShadowColor(10)
        legend.SetMargin(0.2)
        legend.SetTextFont(132)
        legend.AddEntry(plot.histodata, "Data, 330 pb^{-1}","lp")

        #-- loop over histograms (stacked)
        Nstacked = len(self.histosStack)
        #stackColorIndexes = [20,38,14,45,20,38,14,45]
        #stackColorIndexes = [20,38,12,14,20,38,12,14]
        #stackColorIndexes = [2,4,3,14,2,4,3,14]
        stackColorIndexes = [2,4,3,14,92,6,2,4]
        stackFillStyleIds = [3354,3345,3395,3344,3354,3345,3395,3344]
        stkcp = []
        for iter in range(0, Nstacked):
            #make this stack
            stack[iter] = TH1F()
            Nloop = Nstacked - iter
            for iter1 in range(0,Nloop):
                histo = copy.deepcopy(self.histosStack[iter1])
                if(iter1==0):
                    stack[iter] = histo
                    #stack[iter].SetName( self.keysStack[iter] )
                else:
                    stack[iter].Add(histo)
            #define style
            stack[iter].SetMarkerStyle(20+2*iter)
            stack[iter].SetMarkerColor(stackColorIndexes[iter])
            stack[iter].SetLineColor(  stackColorIndexes[iter])
            stack[iter].SetLineWidth(  2 )
            stack[iter].SetFillColor(  stackColorIndexes[iter])
            stack[iter].SetFillStyle(  stackFillStyleIds[iter])
            legend.AddEntry(stack[iter], self.keysStack[Nstacked - iter - 1],"lf")
            #draw stack
            if iter==0:
                stack[iter].SetTitle("")
                stack[iter].GetXaxis().SetTitle(self.xtit)
                stack[iter].GetXaxis().SetTitleFont(132)
                stack[iter].GetXaxis().SetTitleOffset(0.8)
                stack[iter].GetXaxis().SetLabelOffset(0.0)
                stack[iter].GetXaxis().SetTitleSize(0.065)
                stack[iter].GetXaxis().SetLabelSize(0.055)
                stack[iter].GetXaxis().SetLabelFont(132)
                stack[iter].GetYaxis().SetTitleFont(132)
                stack[iter].GetYaxis().SetTitleOffset(0.7)
                stack[iter].GetYaxis().SetTitleSize(0.065)
                stack[iter].GetYaxis().SetLabelSize(0.055)
                stack[iter].GetYaxis().SetLabelOffset(0.0)
                stack[iter].GetYaxis().SetLabelFont(132)
#                stack[iter].GetXaxis().SetTitleFont(132)
#                stack[iter].GetXaxis().SetTitleOffset(0.62)
#                stack[iter].GetXaxis().SetTitleSize(0.05)
#                stack[iter].GetXaxis().SetLabelFont(132)
#                stack[iter].GetXaxis().SetTitleSize(0.05)
#                stack[iter].GetXaxis().SetLabelSize(0.045)
#                stack[iter].GetYaxis().SetTitleFont(132)
#                stack[iter].GetYaxis().SetLabelFont(132)
#                stack[iter].GetYaxis().SetTitleOffset(0.8)
#                stack[iter].GetYaxis().SetTitleSize(0.05)
#                stack[iter].GetYaxis().SetLabelSize(0.045)
                stack[iter].GetYaxis().SetTitle(self.ytit + " #times ("+ str(minBinW) + ")/(bin width)") # units omitted or no units for x-axis
                #stack[iter].GetYaxis().SetTitle((self.ytit + " #times (%.0f GeV)/(bin width)")%(minBinW)) # for x-axis in units of GeV
                if (self.ymin!="" and self.ymax!=""):
                    #stack[iter].GetYaxis().SetLimits(self.ymin,self.ymax)
                    stack[iter].GetYaxis().SetRangeUser(self.ymin,self.ymax)
                #search for maximum of histograms
                #maxHisto = stack[iter].GetMaximum()
                #print maxHisto
                #for hh in self.histos:
                #    if(hh.GetMaximum() > maxHisto):
                #        maxHisto = hh.GetMaximum()
                #stack[iter].GetYaxis().SetLimits(0.,maxHisto*1.2)
                #stack[iter].GetYaxis().SetRangeUser(0.001,maxHisto*1.2)
                #draw first histo
                stack[iter].Draw("HIST")
                stkcp.append(copy.deepcopy(stack[iter]))
            else:
                stkcp.append(copy.deepcopy(stack[iter])) # this is the only way that I figured out to cover the previous histogram!
                stkcp[iter].SetFillStyle(1001)
                stkcp[iter].SetFillColor(10)
                stkcp[iter].Draw("HISTsame")
                stack[iter].Draw("HISTsame")

        #-- Z+jets uncertainty band
        if(self.addZUncBand == "yes"):
            Zhisto = copy.deepcopy(self.histosStack[self.ZPlotIndex])
            zUncHisto = copy.deepcopy(stack[0])
            for bin in range(0,Zhisto.GetNbinsX()):
              zUncHisto.SetBinError(bin+1,self.ZScaleUnc*Zhisto.GetBinContent(bin+1))
            zUncHisto.SetMarkerStyle(0)
            zUncHisto.SetLineColor(0)
            zUncHisto.SetFillColor(5)
            zUncHisto.SetFillStyle(3154)
            zUncHisto.Draw("E2same")
            legend.AddEntry(zUncHisto, self.ZUncKey,"f")

        #-- loop over histograms (overlaid)
        ih=0 # index of histo within a plot
        dataColorIndexes = [1,4,1,1,4,1]
        #dataLineIndexes = [1,2,3,1,2,3]
        dataLineIndexes = [2,1,3,1,2,3]
        for histo in self.histos:
            histo.SetMarkerStyle(dataColorIndexes[ih])
            histo.SetMarkerColor(dataColorIndexes[ih])
            histo.SetLineColor(  dataColorIndexes[ih])
            histo.SetLineStyle(  dataLineIndexes[ih])
            #            histo.SetMarkerStyle(20+2*ih)
            #            histo.SetMarkerColor(2+2*ih)
            #            histo.SetLineColor(  2+2*ih)
            histo.SetLineWidth(2)
            legend.AddEntry(histo, self.keys[ih],"l")
            histo.Draw("HISTsame")
            ih=ih+1

        #-- plot data
        if(self.histodata!=""):
            self.histodata.SetMarkerStyle(20)
            self.histodata.SetLineWidth(2)
            #legend.AddEntry(self.histodata, "Data","lp")
            self.histodata.Draw("psame")

        #-- draw label
        l = TLatex()
        l.SetTextAlign(12)
        l.SetTextFont(132)
        l.SetTextSize(0.065)
        l.SetNDC()
#        l.DrawLatex(xstart,ystart-0.05,"CMS Preliminary 2010")
#        l.DrawLatex(xstart,ystart-0.10,"L_{int} = " + self.lint)
        if (self.lpos=="bottom-center"):
            l.DrawLatex(0.35,0.20,"CMS")
            #            l.DrawLatex(0.35,0.20,"CMS Preliminary 2010")
            #            l.DrawLatex(0.35,0.10,"#intLdt = " + self.lint)
        if (self.lpos=="top-left"):
            l.DrawLatex(xstart+hsize+0.02,ystart+vsize-0.03,"CMS")
#            l.DrawLatex(xstart+hsize+0.02,ystart+vsize-0.03,"CMS Preliminary 2010")
#            l.DrawLatex(xstart+hsize+0.02,ystart+vsize-0.13,"#intLdt = " + self.lint)
        else:
            #l.DrawLatex(xstart-hsize+0.12,ystart+vsize-0.05,"CMS")
            l.DrawLatex(xstart-hsize+0.04,ystart+vsize-0.05,"CMS Preliminary 2011")
            l.DrawLatex(xstart-hsize+0.12,ystart+vsize-0.15,"#sqrt{s} = 7 TeV")
#            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.03,"CMS Preliminary 2010")
#            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.13,"#intLdt = " + self.lint)

        legend.Draw()
        canvas.Update()
        gPad.RedrawAxis()
        gPad.Modified()

        #-- 2nd pad (ratio)
        if(self.makeRatio==1):
            fPads2.cd()
            #fPads2.SetLogy()
            h_bkgTot = copy.deepcopy(stack[0])
            h_ratio = copy.deepcopy(self.histodata)
            h_bkgTot1 = TH1F()
            h_ratio1 = TH1F()

            if( self.xbins!="" and self.rebin!="var" ): ## Variable binning
                xbinsFinal = array( 'd', self.xbins )
                length = len(xbinsFinal)-1
                h_bkgTot1 = h_bkgTot.Rebin( length , "h_bkgTot1", xbinsFinal)
                h_ratio1 = h_ratio.Rebin( length , "h_ratio1" , xbinsFinal)
            else:
                h_bkgTot1 = h_bkgTot
                h_ratio1 = h_ratio

            h_ratio1.SetStats(0)
            if( self.xmin!="" and self.xmax!="" and self.rebin!="var" ):
                h_bkgTot1.GetXaxis().SetRangeUser(self.xmin,self.xmax-0.000001)
                h_ratio1.GetXaxis().SetRangeUser(self.xmin,self.xmax-0.000001)
            h_ratio1.Divide(h_bkgTot1)
            h_ratio1.GetXaxis().SetTitle("")
            h_ratio1.GetXaxis().SetTitleSize(0.06)
            h_ratio1.GetXaxis().SetLabelSize(0.1)
            h_ratio1.GetYaxis().SetLimits(0.,2)
            h_ratio1.GetYaxis().SetRangeUser(0.,2)
            #h_ratio1.GetYaxis().SetRangeUser(0.1,10)
            h_ratio1.GetYaxis().SetTitle("Data/MC")
            h_ratio1.GetYaxis().SetLabelSize(0.1)
            h_ratio1.GetYaxis().SetTitleSize(0.13)
            h_ratio1.GetYaxis().SetTitleOffset(0.3)

            h_ratio1.Draw("p")

            lineAtOne = TLine(h_ratio.GetXaxis().GetXmin(),1,h_ratio.GetXaxis().GetXmax(),1)
            lineAtOne.SetLineColor(2)
            lineAtOne.Draw()


        #-- end
        canvas.SaveAs(self.name + ".eps","eps")
        #canvas.SaveAs(self.name + ".png","png")
        #canvas.SaveAs(self.name + ".root","root")
        #canvas.SaveAs(self.name + ".pdf","pdf") # do not use this line because root creates rotated pdf plot - see end of the file instead
        canvas.Print(fileps)



############# DON'T NEED TO MODIFY ANYTHING HERE - END #######################
##############################################################################


##############################################################################
############# USER CODE - BEGING #############################################

#--- Input root file


File_preselection = GetFile("/home/santanas/Axigluons/data/output_fromAFS/axigluons_enujj/330pb-1_21072011/output_cutTable_axigluons_enujj/analysisClass_axigluons_enujj_plots.root")
#File_preselection = GetFile("/home/santanas/Axigluons/data/output_fromAFS/axigluons_enujj/330pb-1_18072011/output_cutTable_axigluons_enujj/analysisClass_axigluons_enujj_plots.root")
File_selection    = File_preselection

#UseQCDFromData    = 1 # always put an existing file under File_QCD (otherwise the code will crash)
#File_QCD = GetFile("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/35.8pb-1_QCD_UseHLTPrescales_presel_MET45_presel_sT250_Feb112011/analysisClass_enujjSample_QCD_plots.root")
#QCDscaleFactor    = 1 # no need to rescale anymore since we are using the HLT prescales (36/35.84 can be ignored)

#### Common values for plots:
#otherBkgsKey="QCD, single top, VV+jets, W/W*+jets"
otherBkgsKey="Other backgrounds"
zUncBand="no"
makeRatio=1
doExtraPlots = False

counter_ymax=100000000

pt_xmin=0
pt_xmax_lin=200
pt_xmax_log=1000
pt_xmax_log_high=2000 #for Mjj, sT, Mej, etc..
pt_ymin=0.01
pt_ymax_lin=2000
pt_ymax_log=10000

eta_rebin=2
eta_ymin=0
eta_ymax=1000

#--- Final plots are defined here

# Simply add or remove plots and update the list plots = [plot0, plot1, ...] below

histoBaseName = "histo1D__SAMPLE__cutHisto_allPreviousCuts________VARIABLE"
histoBaseName_userDef = "histo1D__SAMPLE__VARIABLE"

#samplesForStackHistosQCD = ["DATA"]
#samplesForStackHistos = ["OTHERBKG","DIBOSON","TTbar_Madgraph","WJet_Madgraph"]
#keysStack =             [otherBkgsKey,"WW + WZ + ZZ","t#bar{t}", "W/W* + jets"]
#keysStack =             ["QCD multijet",otherBkgsKey,"t#bar{t}", "W/W* + jets"]
samplesForStackHistos = ["ZJet_Madgraph","PhotonJets","SingleTop","DIBOSON","TTbar_Madgraph","WJet_Madgraph"]
keysStack =             ["Z/Z* + jets","#gamma + jets","single top","WW + WZ + ZZ","t#bar{t}", "W/W* + jets"]

#samplesForHistos = ["LQenujj_M250", "LQenujj_M300","LQenujj_M340"]
#keys             = ["LQ, M = 250 GeV","LQ, M = 300 GeV","LQ, M = 340 GeV"]
samplesForHistos = []
keys             = []

sampleForDataHisto = "DATA"

#--- nEle ---
variableName = "nEle"

plot0 = Plot()
## inputs for stacked histograms
#plot0.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot0.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot0.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot0.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot0.keys            = keys
plot0.xtit            = "Number of electrons (Preselection)"
plot0.ytit            = "Events"
plot0.ylog            = "yes"
plot0.rebin           = 1
plot0.xmin            = -0.5
plot0.xmax            = 6.5
plot0.ymin            = 0.1
plot0.ymax            = counter_ymax
#plot0.lpos = "bottom-center"
plot0.name            = "nEle_preselection"
plot0.addZUncBand     = zUncBand
plot0.makeRatio       = makeRatio
plot0.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

#--- Ele1_Pt  ---
variableName = "Ele1_Pt"

plot1_ylin = Plot()
## inputs for stacked histograms
#plot1_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot1_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot1_ylin.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot1_ylin.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot1_ylin.keys            = keys
plot1_ylin.xtit            = "p_{T} electron [GeV] (Preselection)"
plot1_ylin.ytit            = "Events"
plot1_ylin.ylog            = "no"
#plot1_ylin.rebin           = "var"
plot1_ylin.xmin            = pt_xmin
plot1_ylin.xmax            = pt_xmax_lin
plot1_ylin.ymin            = pt_ymin
plot1_ylin.ymax            = pt_ymax_lin
#plot1_ylin.lpos = "bottom-center"
plot1_ylin.name            = "Ele1_Pt_preselection_ylin"
plot1_ylin.addZUncBand     = zUncBand
plot1_ylin.makeRatio       = makeRatio
plot1_ylin.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,pt_xmax_lin]
plot1_ylin.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

plot1_ylog = Plot()
## inputs for stacked histograms
#plot1_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot1_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot1_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot1_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot1_ylog.keys            = keys
plot1_ylog.xtit            = "p_{T} electron [GeV] (Preselection)"
plot1_ylog.ytit            = "Events"
plot1_ylog.ylog            = "yes"
#plot1_ylog.rebin           = "var"
plot1_ylog.xmin            = pt_xmin
plot1_ylog.xmax            = pt_xmax_log
plot1_ylog.ymin            = pt_ymin
plot1_ylog.ymax            = pt_ymax_log
#plot1_ylog.lpos = "bottom-center"
plot1_ylog.name            = "Ele1_Pt_preselection_ylog"
plot1_ylog.addZUncBand     = zUncBand
plot1_ylog.makeRatio       = makeRatio
plot1_ylog.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,200,300,400,600,800,pt_xmax_log]
plot1_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

#--- Ele1_Eta  ---
variableName = "Ele1_Eta"

plot2 = Plot()
## inputs for stacked histograms
#plot2.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot2.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot2.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot2.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot2.keys            = keys
plot2.xtit            = "#eta electron (Preselection)"
plot2.ytit            = "Events"
plot2.rebin           = eta_rebin
plot2.xmin            = -3
plot2.xmax            = 3
plot2.ymin            = eta_ymin
plot2.ymax            = eta_ymax
#plot2.lpos            = "top-left"
plot2.name            = "Ele1_Eta_preselection_ylin"
plot2.addZUncBand     = zUncBand
plot2.makeRatio       = makeRatio
plot2.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

#--- Ele1_Phi  ---
variableName = "Ele1_Phi"

plot3 = Plot()
## inputs for stacked histograms
#plot3.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot3.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot3.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot3.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot3.keys            = keys
plot3.xtit            = "#phi electron [rad.] (Preselection)"
plot3.ytit            = "Events"
plot3.rebin           = eta_rebin*2
#plot3.xmin            = -3.15
#plot3.xmax            = 3.15
plot3.ymin            = eta_ymin
plot3.ymax            = eta_ymax
#plot3.lpos            = "top-left"
plot3.name            = "Ele1_Phi_preselection_ylin"
plot3.addZUncBand     = zUncBand
plot3.makeRatio       = makeRatio
plot3.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- Ele1_Charge  ---
variableName = "Ele1_Charge"

plot4 = Plot()
## inputs for stacked histograms
#plot4.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot4.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot4.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot4.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot4.keys            = keys
plot4.xtit            = "Charge electron (Preselection)"
plot4.ytit            = "Events"
plot4.ylog            = "no"
plot4.rebin           = 1
#plot4.xmin            = -1.0001
#plot4.xmax            = 1.0001
plot4.ymin            = 1
plot4.ymax            = pt_ymax_log
#plot4.lpos = "bottom-center"
plot4.name            = "Ele1_Eta_preselection_ylin"
plot4.addZUncBand     = zUncBand
plot4.makeRatio       = makeRatio
plot4.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- DR_Ele1Jets ---
variableNames = ["DR_Ele1Jet1","DR_Ele1Jet2"]

plot5 = Plot()
## inputs for stacked histograms
#plot5.histosStack     = generateAndAddHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateAndAddHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot5.histosStack     = generateAndAddHistoList( histoBaseName_userDef, samplesForStackHistos, variableNames, File_preselection)
plot5.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot5.histos          = generateAndAddHistoList( histoBaseName_userDef, samplesForHistos, variableNames, File_preselection)
plot5.keys            = keys
plot5.xtit            = "DeltaR(e,1^{st} and 2^{nd} jets)"
plot5.ytit            = "2 X Events"
plot5.ylog            = "yes"
plot5.rebin           = 2
plot5.xmin            = 0
plot5.xmax            = 7
plot5.ymin            = pt_ymin
plot5.ymax            = pt_ymax_log*1000
#plot5.lpos = "bottom-center"
plot5.name            = "DR_Ele1Jets_preselection_ylog"
plot5.addZUncBand     = zUncBand
plot5.makeRatio       = makeRatio
plot5.histodata       = generateAndAddHisto( histoBaseName_userDef, sampleForDataHisto, variableNames, File_preselection)



#--- MET_Pt  ---
variableName = "MET_Pt"

plot20_ylin = Plot()
## inputs for stacked histograms
#plot20_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot20_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot20_ylin.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot20_ylin.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot20_ylin.keys            = keys
plot20_ylin.xtit            = "MET [GeV] (Preselection)"
plot20_ylin.ytit            = "Events"
plot20_ylin.ylog            = "no"
#plot20_ylin.rebin           = "var"
plot20_ylin.xmin            = 0
plot20_ylin.xmax            = pt_xmax_lin
plot20_ylin.ymin            = 0.01
plot20_ylin.ymax            = pt_ymax_lin
#plot20_ylin.lpos = "bottom-center"
plot20_ylin.name            = "MET_Pt_preselection_ylin"
plot20_ylin.addZUncBand     = zUncBand
plot20_ylin.makeRatio       = makeRatio
plot20_ylin.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,pt_xmax_lin]
plot20_ylin.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

plot20_ylog = Plot()
## inputs for stacked histograms
#plot20_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot20_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot20_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot20_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot20_ylog.keys            = keys
plot20_ylog.xtit            = "MET [GeV] (Preselection)"
plot20_ylog.ytit            = "Events"
plot20_ylog.ylog            = "yes"
#plot20_ylog.rebin           = "var"
plot20_ylog.xmin            = 0
plot20_ylog.xmax            = pt_xmax_log
plot20_ylog.ymin            = 0.01
plot20_ylog.ymax            = pt_ymax_log
#plot20_ylog.lpos = "bottom-center"
plot20_ylog.name            = "MET_Pt_preselection_ylog"
plot20_ylog.addZUncBand     = zUncBand
plot20_ylog.makeRatio       = makeRatio
plot20_ylog.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,200,300,400,500,600,800,pt_xmax_log]
plot20_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- MET_Phi  ---
variableName = "MET_Phi"

plot21 = Plot()
## inputs for stacked histograms
#plot21.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot21.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot21.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot21.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot21.keys            = keys
plot21.xtit            = "MET #phi [rad.] (Preselection)"
plot21.ytit            = "Events"
plot21.rebin           = eta_rebin*2
plot21.ymin            = eta_ymin
plot21.ymax            = eta_ymax
#plot21.lpos = "bottom-center"
plot21.name            = "MET_Phi_preselection_ylin"
plot21.addZUncBand     = zUncBand
plot21.makeRatio       = makeRatio
plot21.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- MT_Ele1MET ---
variableName = "MT_Ele1MET"

plot22_ylin = Plot()
## inputs for stacked histograms
#plot22_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot22_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot22_ylin.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot22_ylin.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot22_ylin.keys            = keys
plot22_ylin.xtit            = "m_{T,e#nu} [GeV] (Preselection)"
plot22_ylin.ytit            = "Events"
plot22_ylin.ylog            = "no"
plot22_ylin.rebin           = 2
plot22_ylin.xmin            = 0
plot22_ylin.xmax            = pt_xmax_lin 
plot22_ylin.ymin            = 0
plot22_ylin.ymax            = pt_ymax_lin
#plot22_ylin.lpos = "bottom-center"
plot22_ylin.name            = "MT_Ele1MET_preselection_ylin"
plot22_ylin.addZUncBand     = zUncBand
plot22_ylin.makeRatio       = makeRatio
plot22_ylin.xbins           = [0,10,20,30,40,50,60,70,80,90,100,120,130,140,150,160,170,180,190,pt_xmax_lin]
plot22_ylin.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

plot22_ylog = Plot()
## inputs for stacked histograms
#plot22_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot22_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot22_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot22_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot22_ylog.keys            = keys
plot22_ylog.xtit            = "m_{T,e#nu} [GeV] (Preselection)"
plot22_ylog.ytit            = "Events"
plot22_ylog.ylog            = "yes"
plot22_ylog.rebin           = 2
plot22_ylog.xmin            = 0
plot22_ylog.xmax            = pt_xmax_log
plot22_ylog.ymin            = 0.01
plot22_ylog.ymax            = pt_ymax_log
#plot22_ylog.lpos = "bottom-center"
plot22_ylog.name            = "MT_Ele1MET_preselection_ylog"
plot22_ylog.addZUncBand     = zUncBand
plot22_ylog.makeRatio       = makeRatio                                   
plot22_ylog.xbins           = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,350,400,500,600,800,pt_xmax_log]
plot22_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- mDPhi_METEle1 ---
variableName = "mDPhi_METEle1"

plot23_ylin = Plot()
## inputs for stacked histograms
#plot23_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot23_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot23_ylin.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot23_ylin.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot23_ylin.keys            = keys
plot23_ylin.xtit            = "#Delta#phi(MET,e) [rad.] (Preselection)"
plot23_ylin.ytit            = "Events"
#plot23_ylin.xlog            = "yes"
plot23_ylin.ylog            = "no"
plot23_ylin.rebin           = eta_rebin
plot23_ylin.ymin            = eta_ymin
plot23_ylin.ymax            = eta_ymax
#plot23_ylin.lpos = "bottom-center"
plot23_ylin.name            = "mDPhi_METEle1_preselection_ylin"
plot23_ylin.addZUncBand     = zUncBand
plot23_ylin.makeRatio       = makeRatio
plot23_ylin.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

plot23_ylog = Plot()
## inputs for stacked histograms
#plot23_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot23_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot23_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot23_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot23_ylog.keys            = keys
plot23_ylog.xtit            = "#Delta#phi(MET,e) [rad.] (Preselection)"
plot23_ylog.ytit            = "Events"
#plot23_ylog.xlog            = "yes"
plot23_ylog.ylog            = "yes"
plot23_ylog.rebin           = eta_rebin
plot23_ylog.ymin            = 0.1
plot23_ylog.ymax            = eta_ymax*10000
#plot23_ylog.lpos = "bottom-center"
plot23_ylog.name            = "mDPhi_METEle1_preselection_ylog"
plot23_ylog.addZUncBand     = zUncBand
plot23_ylog.makeRatio       = makeRatio
plot23_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

#--- mDPhi_METJet1 ---
variableName = "mDPhi_METJet1"

plot24_ylin = Plot()
## inputs for stacked histograms
#plot24_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot24_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot24_ylin.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot24_ylin.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot24_ylin.keys            = keys
plot24_ylin.xtit            = "#Delta#phi(MET,1^{st} jet) [rad.] (Preselection)"
plot24_ylin.ytit            = "Events"
#plot24_ylin.xlog            = "yes"
plot24_ylin.ylog            = "no"
plot24_ylin.rebin           = eta_rebin
plot24_ylin.ymin            = eta_ymin
plot24_ylin.ymax            = eta_ymax + eta_ymax/2 
#plot24_ylin.lpos = "bottom-center"
plot24_ylin.name            = "mDPhi_METJet1_preselection_ylin"
plot24_ylin.addZUncBand     = zUncBand
plot24_ylin.makeRatio       = makeRatio
plot24_ylin.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

plot24_ylog = Plot()
## inputs for stacked histograms
#plot24_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot24_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot24_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot24_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot24_ylog.keys            = keys
plot24_ylog.xtit            = "#Delta#phi(MET,1^{st} jet) [rad.] (Preselection)"
plot24_ylog.ytit            = "Events"
#plot24_ylog.xlog            = "yes"
plot24_ylog.ylog            = "yes"
plot24_ylog.rebin           = eta_rebin
plot24_ylog.ymin            = 0.1
plot24_ylog.ymax            = eta_ymax*10000
#plot24_ylog.lpos = "bottom-center"
plot24_ylog.name            = "mDPhi_METJet1_preselection_ylog"
plot24_ylog.addZUncBand     = zUncBand
plot24_ylog.makeRatio       = makeRatio
plot24_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

#--- mDPhi_METJet2 ---
variableName = "mDPhi_METJet2"

plot25_ylin = Plot()
## inputs for stacked histograms
#plot25_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot25_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot25_ylin.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot25_ylin.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot25_ylin.keys            = keys
plot25_ylin.xtit            = "#Delta#phi(MET,2^{nd} jet) [rad.] (Preselection)"
plot25_ylin.ytit            = "Events"
#plot25_ylin.xlog            = "yes"
plot25_ylin.ylog            = "no"
plot25_ylin.rebin           = eta_rebin
plot25_ylin.ymin            = eta_ymin
plot25_ylin.ymax            = eta_ymax
#plot25_ylin.lpos = "bottom-center"
plot25_ylin.name            = "mDPhi_METJet2_preselection_ylin"
plot25_ylin.addZUncBand     = zUncBand
plot25_ylin.makeRatio       = makeRatio
plot25_ylin.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

plot25_ylog = Plot()
## inputs for stacked histograms
#plot25_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot25_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot25_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot25_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot25_ylog.keys            = keys
plot25_ylog.xtit            = "#Delta#phi(MET,2^{nd} jet) [rad.] (Preselection)"
plot25_ylog.ytit            = "Events"
#plot25_ylog.xlog            = "yes"
plot25_ylog.ylog            = "yes"
plot25_ylog.rebin           = eta_rebin
plot25_ylog.ymin            = 0.1
plot25_ylog.ymax            = eta_ymax*10000
#plot25_ylog.lpos = "bottom-center"
plot25_ylog.name            = "mDPhi_METJet2_preselection_ylog"
plot25_ylog.addZUncBand     = zUncBand
plot25_ylog.makeRatio       = makeRatio
plot25_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

#--- nJet ---
variableName = "nJet"

plot40 = Plot()
## inputs for stacked histograms
#plot40.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot40.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot40.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot40.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot40.keys            = keys
plot40.xtit            = "Number of jets (Preselection)"
plot40.ytit            = "Events"
plot40.ylog            = "yes"
plot40.rebin           = 1
plot40.xmin            = -0.5
plot40.xmax            = 15.5
plot40.ymin            = 0.1
plot40.ymax            = counter_ymax
#plot40.lpos = "bottom-center"
plot40.name            = "nJet_preselection"
plot40.addZUncBand     = zUncBand
plot40.makeRatio       = makeRatio
plot40.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

#--- nJet_btagTCHE ---
variableName = "nJet_btagTCHE"

plot41 = Plot()
## inputs for stacked histograms
#plot41.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot41.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot41.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot41.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot41.keys            = keys
plot41.xtit            = "Number of b-tagged jets [TCHEL] (Preselection)"
plot41.ytit            = "Events"
plot41.ylog            = "yes"
plot41.rebin           = 1
plot41.xmin            = -0.5
plot41.xmax            = 15.5
plot41.ymin            = 0.1
plot41.ymax            = counter_ymax
#plot41.lpos = "bottom-center"
plot41.name            = "nJet_btagTCHE_preselection"
plot41.addZUncBand     = zUncBand
plot41.makeRatio       = makeRatio
plot41.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

#--- Jet1_Pt  ---
variableName = "Jet1_Pt"

plot42_ylin = Plot()
## inputs for stacked histograms
#plot42_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot42_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot42_ylin.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot42_ylin.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot42_ylin.keys            = keys
plot42_ylin.xtit            = "p_{T} 1^{st} jet [GeV] (Preselection)"
plot42_ylin.ytit            = "Events"
plot42_ylin.ylog            = "no"
#plot42_ylin.rebin           = "var"
plot42_ylin.xmin            = pt_xmin
plot42_ylin.xmax            = pt_xmax_lin
plot42_ylin.ymin            = pt_ymin
plot42_ylin.ymax            = pt_ymax_lin
#plot42_ylin.lpos = "bottom-center"
plot42_ylin.name            = "Jet1_Pt_preselection_ylin"
plot42_ylin.addZUncBand     = zUncBand
plot42_ylin.makeRatio       = makeRatio
plot42_ylin.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,pt_xmax_lin]
plot42_ylin.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

plot42_ylog = Plot()
## inputs for stacked histograms
#plot42_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot42_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot42_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot42_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot42_ylog.keys            = keys
plot42_ylog.xtit            = "p_{T} 1^{st} jet [GeV] (Preselection)"
plot42_ylog.ytit            = "Events"
plot42_ylog.ylog            = "yes"
#plot42_ylog.rebin           = "var"
plot42_ylog.xmin            = pt_xmin
plot42_ylog.xmax            = pt_xmax_log
plot42_ylog.ymin            = pt_ymin
plot42_ylog.ymax            = pt_ymax_log
#plot42_ylog.lpos = "bottom-center"
plot42_ylog.name            = "Jet1_Pt_preselection_ylog"
plot42_ylog.addZUncBand     = zUncBand
plot42_ylog.makeRatio       = makeRatio
plot42_ylog.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,200,300,400,600,800,pt_xmax_log]
plot42_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

#--- Jet1_Eta  ---
variableName = "Jet1_Eta"

plot43 = Plot()
## inputs for stacked histograms
#plot43.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot43.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot43.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot43.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot43.keys            = keys
plot43.xtit            = "#eta 1^{st} jet (Preselection)"
plot43.ytit            = "Events"
plot43.rebin           = eta_rebin
plot43.xmin            = -5
plot43.xmax            = 5
plot43.ymin            = eta_ymin
plot43.ymax            = eta_ymax
#plot43.lpos            = "top-left"
plot43.name            = "Jet1_Eta_preselection_ylin"
plot43.addZUncBand     = zUncBand
plot43.makeRatio       = makeRatio
plot43.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

#--- Jet1_Phi  ---
variableName = "Jet1_Phi"

plot44 = Plot()
## inputs for stacked histograms
#plot44.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot44.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot44.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot44.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot44.keys            = keys
plot44.xtit            = "#phi 1^{st} jet [rad.] (Preselection)"
plot44.ytit            = "Events"
plot44.rebin           = eta_rebin*2
#plot44.xmin            = -3.15
#plot44.xmax            = 3.15
plot44.ymin            = eta_ymin
plot44.ymax            = eta_ymax
#plot44.lpos            = "top-left"
plot44.name            = "Jet1_Phi_preselection_ylin"
plot44.addZUncBand     = zUncBand
plot44.makeRatio       = makeRatio
plot44.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- Jet1_btagTCHE ---
variableName = "Jet1_btagTCHE"

plot45_ylin = Plot()
## inputs for stacked histograms
#plot45_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot45_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot45_ylin.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot45_ylin.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot45_ylin.keys            = keys
plot45_ylin.xtit            = "B-tag discriminator [TCHE] 1^{st} jet (Preselection)"
plot45_ylin.ytit            = "Events"
plot45_ylin.ylog            = "no"
#plot45_ylin.rebin           = "var"
plot45_ylin.ymin            = pt_ymin
plot45_ylin.ymax            = pt_ymax_lin / 2
plot45_ylin.xmin            = 0
plot45_ylin.xmax            = 10
#plot45_ylin.lpos            = "top-left"
plot45_ylin.name            = "Jet1_btagTCHE_preselection_ylin"
plot45_ylin.addZUncBand     = zUncBand
plot45_ylin.makeRatio       = makeRatio
#plot45_ylin.xbins           = [0,0.4,0.8,1.2,1.6,2.0,2.4,2.8,4,5,7,10,15,20,30,40,50]
plot45_ylin.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

plot45_ylog = Plot()
## inputs for stacked histograms
#plot45_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot45_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot45_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot45_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot45_ylog.keys            = keys
plot45_ylog.xtit            = "B-tag discriminator [TCHE] 1^{st} jet (Preselection)"
plot45_ylog.ytit            = "Events"
plot45_ylog.ylog            = "yes"
#plot45_ylog.rebin           = "var"
plot45_ylog.ymin            = pt_ymin
plot45_ylog.ymax            = pt_ymax_log*10
#plot45_ylog.lpos            = "top-left"
plot45_ylog.name            = "Jet1_btagTCHE_preselection_ylog"
plot45_ylog.addZUncBand     = zUncBand
plot45_ylog.makeRatio       = makeRatio
plot45_ylog.xbins           = [0,0.4,0.8,1.2,1.6,2.0,2.4,2.8,4,5,7,10,15,20,30,40,50]
plot45_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)



#--- Jet2_Pt  ---
variableName = "Jet2_Pt"

plot46_ylin = Plot()
## inputs for stacked histograms
#plot46_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot46_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot46_ylin.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot46_ylin.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot46_ylin.keys            = keys
plot46_ylin.xtit            = "p_{T} 2^{nd} jet [GeV] (Preselection)"
plot46_ylin.ytit            = "Events"
plot46_ylin.ylog            = "no"
#plot46_ylin.rebin           = "var"
plot46_ylin.xmin            = pt_xmin
plot46_ylin.xmax            = pt_xmax_lin
plot46_ylin.ymin            = pt_ymin
plot46_ylin.ymax            = pt_ymax_lin
#plot46_ylin.lpos = "bottom-center"
plot46_ylin.name            = "Jet2_Pt_preselection_ylin"
plot46_ylin.addZUncBand     = zUncBand
plot46_ylin.makeRatio       = makeRatio
plot46_ylin.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,pt_xmax_lin]
plot46_ylin.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

plot46_ylog = Plot()
## inputs for stacked histograms
#plot46_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot46_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot46_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot46_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot46_ylog.keys            = keys
plot46_ylog.xtit            = "p_{T} 2^{nd} jet [GeV] (Preselection)"
plot46_ylog.ytit            = "Events"
plot46_ylog.ylog            = "yes"
#plot46_ylog.rebin           = "var"
plot46_ylog.xmin            = pt_xmin
plot46_ylog.xmax            = pt_xmax_log
plot46_ylog.ymin            = pt_ymin
plot46_ylog.ymax            = pt_ymax_log
#plot46_ylog.lpos = "bottom-center"
plot46_ylog.name            = "Jet2_Pt_preselection_ylog"
plot46_ylog.addZUncBand     = zUncBand
plot46_ylog.makeRatio       = makeRatio
plot46_ylog.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,200,300,400,600,800,pt_xmax_log]
plot46_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

#--- Jet2_Eta  ---
variableName = "Jet2_Eta"

plot47 = Plot()
## inputs for stacked histograms
#plot47.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot47.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot47.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot47.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot47.keys            = keys
plot47.xtit            = "#eta 2^{nd} jet (Preselection)"
plot47.ytit            = "Events"
plot47.rebin           = eta_rebin
plot47.xmin            = -5
plot47.xmax            = 5
plot47.ymin            = eta_ymin
plot47.ymax            = eta_ymax
#plot47.lpos            = "top-left"
plot47.name            = "Jet2_Eta_preselection_ylin"
plot47.addZUncBand     = zUncBand
plot47.makeRatio       = makeRatio
plot47.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

#--- Jet2_Phi  ---
variableName = "Jet2_Phi"

plot48 = Plot()
## inputs for stacked histograms
#plot48.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot48.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot48.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot48.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot48.keys            = keys
plot48.xtit            = "#phi 2^{nd} jet [rad.] (Preselection)"
plot48.ytit            = "Events"
plot48.rebin           = eta_rebin*2
#plot48.xmin            = -3.15
#plot48.xmax            = 3.15
plot48.ymin            = eta_ymin
plot48.ymax            = eta_ymax
#plot48.lpos            = "top-left"
plot48.name            = "Jet2_Phi_preselection_ylin"
plot48.addZUncBand     = zUncBand
plot48.makeRatio       = makeRatio
plot48.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- Jet2_btagTCHE ---
variableName = "Jet2_btagTCHE"

plot49_ylin = Plot()
## inputs for stacked histograms
#plot49_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot49_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot49_ylin.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot49_ylin.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot49_ylin.keys            = keys
plot49_ylin.xtit            = "B-tag discriminator [TCHE] 2^{nd} jet (Preselection)"
plot49_ylin.ytit            = "Events"
plot49_ylin.ylog            = "no"
#plot49_ylin.rebin           = "var"
plot49_ylin.ymin            = pt_ymin
plot49_ylin.ymax            = pt_ymax_lin / 2
plot49_ylin.xmin            = 0
plot49_ylin.xmax            = 10
#plot49_ylin.lpos            = "top-left"
plot49_ylin.name            = "Jet2_btagTCHE_preselection_ylin"
plot49_ylin.addZUncBand     = zUncBand
plot49_ylin.makeRatio       = makeRatio
#plot49_ylin.xbins           = [0,0.4,0.8,1.2,1.6,2.0,2.4,2.8,4,5,7,10,15,20,30,40,50]
plot49_ylin.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

plot49_ylog = Plot()
## inputs for stacked histograms
#plot49_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot49_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot49_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot49_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot49_ylog.keys            = keys
plot49_ylog.xtit            = "B-tag discriminator [TCHE] 2^{nd} jet (Preselection)"
plot49_ylog.ytit            = "Events"
plot49_ylog.ylog            = "yes"
#plot49_ylog.rebin           = "var"
plot49_ylog.ymin            = pt_ymin
plot49_ylog.ymax            = pt_ymax_log*10
#plot49_ylog.lpos            = "top-left"
plot49_ylog.name            = "Jet2_btagTCHE_preselection_ylog"
plot49_ylog.addZUncBand     = zUncBand
plot49_ylog.makeRatio       = makeRatio
plot49_ylog.xbins           = [0,0.4,0.8,1.2,1.6,2.0,2.4,2.8,4,5,7,10,15,20,30,40,50]
plot49_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- Ptjets ---
variableNames = ["Jet1_Pt","Jet2_Pt"]

plot50_ylin = Plot()
## inputs for stacked histograms
#plot50_ylin.histosStack     = generateAndAddHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableNames, File_QCD, QCDscaleFactor) + generateAndAddHistoList( histoBaseName_userDef, samplesForStackHistos, variableNames, File_preselection)
plot50_ylin.histosStack     = generateAndAddHistoList( histoBaseName_userDef, samplesForStackHistos, variableNames, File_preselection)
plot50_ylin.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot50_ylin.histos          = generateAndAddHistoList( histoBaseName_userDef, samplesForHistos, variableNames, File_preselection)
plot50_ylin.keys            = keys
plot50_ylin.xtit            = "p_{T} 1^{st} and 2^{nd} jets [GeV] (Preselection)"
plot50_ylin.ytit            = "2 x Events"
plot50_ylin.ylog            = "no"
#plot50_ylin.rebin           = "var"
plot50_ylin.xmin            = pt_xmin
plot50_ylin.xmax            = pt_xmax_lin
plot50_ylin.ymin            = pt_ymin
plot50_ylin.ymax            = pt_ymax_lin*1.5
#plot50_ylin.lpos = "bottom-center"
plot50_ylin.name            = "Jets1and2_Pt_preselection_ylin"
plot50_ylin.addZUncBand     = zUncBand
plot50_ylin.makeRatio       = makeRatio
plot50_ylin.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,pt_xmax_lin]
plot50_ylin.histodata       = generateAndAddHisto( histoBaseName_userDef, sampleForDataHisto, variableNames, File_preselection)

plot50_ylog = Plot()
## inputs for stacked histograms
#plot50_ylog.histosStack     = generateAndAddHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableNames, File_QCD, QCDscaleFactor) + generateAndAddHistoList( histoBaseName_userDef, samplesForStackHistos, variableNames, File_preselection)
plot50_ylog.histosStack     = generateAndAddHistoList( histoBaseName_userDef, samplesForStackHistos, variableNames, File_preselection)
plot50_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot50_ylog.histos          = generateAndAddHistoList( histoBaseName_userDef, samplesForHistos, variableNames, File_preselection)
plot50_ylog.keys            = keys
plot50_ylog.xtit            = "p_{T} 1^{st} and 2^{nd} jets [GeV] (Preselection)"
plot50_ylog.ytit            = "2 x Events"
plot50_ylog.ylog            = "yes"
#plot50_ylog.rebin           = "var"
plot50_ylog.xmin            = pt_xmin
plot50_ylog.xmax            = pt_xmax_log
plot50_ylog.ymin            = pt_ymin
plot50_ylog.ymax            = pt_ymax_log*2
#plot50_ylog.lpos = "bottom-center"
plot50_ylog.name            = "Jets1and2_Pt_preselection_ylog"
plot50_ylog.addZUncBand     = zUncBand
plot50_ylog.makeRatio       = makeRatio
plot50_ylog.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,200,300,400,600,800,pt_xmax_log]
plot50_ylog.histodata       = generateAndAddHisto( histoBaseName_userDef, sampleForDataHisto, variableNames, File_preselection)


#--- Etajets  ---
variableNames = ["Jet1_Eta","Jet2_Eta"]

plot51 = Plot()
## inputs for stacked histograms
#plot51.histosStack     = generateAndAddHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableNames, File_QCD, QCDscaleFactor) + generateAndAddHistoList( histoBaseName_userDef, samplesForStackHistos, variableNames, File_preselection)
plot51.histosStack     = generateAndAddHistoList( histoBaseName_userDef, samplesForStackHistos, variableNames, File_preselection)
plot51.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot51.histos          = generateAndAddHistoList( histoBaseName_userDef, samplesForHistos, variableNames, File_preselection)
plot51.keys            = keys
plot51.xtit            = "#eta 1^{st}and 2^{nd} jets (Preselection)"
plot51.ytit            = "2 x Events"
plot51.rebin           = eta_rebin
plot51.xmin            = -5
plot51.xmax            = 5
plot51.ymin            = eta_ymin
plot51.ymax            = eta_ymax*1.5
#plot51.lpos            = "top-left"
plot51.name            = "Jets1and2_Eta_preselection_ylin"
plot51.addZUncBand     = zUncBand
plot51.makeRatio       = makeRatio
plot51.histodata       = generateAndAddHisto( histoBaseName_userDef, sampleForDataHisto, variableNames, File_preselection)


#--- M_j1j2 ---
variableName = "M_j1j2"

plot52_ylin = Plot()
## inputs for stacked histograms
#plot52_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableNames, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableNames, File_preselection)
plot52_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot52_ylin.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot52_ylin.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot52_ylin.keys            = keys
plot52_ylin.xtit            = "M(j1j2) [GeV]"
plot52_ylin.ytit            = "Events"
plot52_ylin.ylog            = "no"
#plot52_ylin.rebin           = "var"
plot52_ylin.ymin            = pt_ymin
plot52_ylin.ymax            = pt_ymax_lin/4
plot52_ylin.xmin            = pt_xmin
#plot52_ylin.xmax            = pt_xmax_lin*2.5  ==> FIXME THIS DOS NOT SHOW THE OVERFLOW BIN CONTNET IN THE LAST BIN OF THE PLOT
plot52_ylin.xmax            = pt_xmax_lin*3  
#plot52_ylin.lpos = "bottom-center"
plot52_ylin.name            = "M_j1j2_preselection_ylin"
plot52_ylin.addZUncBand     = zUncBand
plot52_ylin.makeRatio       = makeRatio
#plot52_ylin.xbins           = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,160,180,200,220,240,260,280,300,350,pt_xmax_lin*2]
plot52_ylin.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


plot52_ylog = Plot()
## inputs for stacked histograms
#plot52_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableNames, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableNames, File_preselection)
plot52_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot52_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot52_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot52_ylog.keys            = keys
plot52_ylog.xtit            = "M(j1j2) [GeV]"
plot52_ylog.ytit            = "Events"
plot52_ylog.ylog            = "yes"
#plot52_ylog.rebin           = "var"
plot52_ylog.ymin            = pt_ymin
plot52_ylog.ymax            = pt_ymax_log
#plot52_ylog.xmin            = pt_xmin
plot52_ylog.xmax            = pt_xmax_log_high
#plot52_ylog.lpos = "bottom-center"
plot52_ylog.name            = "M_j1j2_preselection_ylog"
plot52_ylog.addZUncBand     = zUncBand
plot52_ylog.makeRatio       = makeRatio
plot52_ylog.xbins           = [0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,160,180,200,220,240,260,280,320,340,360,380,400,420,460,500,550,
                               600,650,700,750,800,850,900,950,1000,1050,1100,1200,1300,1400,1500,1600,1700,1800,1900,pt_xmax_log_high]
plot52_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)



#--- nMuon ---
variableName = "nMuon"

plot100 = Plot()
## inputs for stacked histograms
#plot100.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot100.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot100.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot100.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot100.keys            = keys
plot100.xtit            = "Number of muons (Preselection)"
plot100.ytit            = "Events"
plot100.ylog            = "yes"
plot100.rebin           = 1
plot100.xmin            = -0.5
plot100.xmax            = 6.5
plot100.ymin            = 0.1
plot100.ymax            = counter_ymax
#plot100.lpos = "bottom-center"
plot100.name            = "nMuon_preselection"
plot100.addZUncBand     = zUncBand
plot100.makeRatio       = makeRatio
plot100.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- Pt_Ele1MET  ---
variableName = "Pt_Ele1MET"

plot120_ylin = Plot()
## inputs for stacked histograms
#plot120_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot120_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot120_ylin.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot120_ylin.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot120_ylin.keys            = keys
plot120_ylin.xtit            = "p_{T}(e#nu) [GeV] (Preselection)"
plot120_ylin.ytit            = "Events"
plot120_ylin.ylog            = "no"
#plot120_ylin.rebin           = "var"
plot120_ylin.xmin            = pt_xmin
plot120_ylin.xmax            = pt_xmax_lin + pt_xmax_lin/2
plot120_ylin.ymin            = pt_ymin
plot120_ylin.ymax            = pt_ymax_lin/2
#plot120_ylin.lpos = "bottom-center"
plot120_ylin.name            = "Pt_Ele1MET_preselection_ylin"
plot120_ylin.addZUncBand     = zUncBand
plot120_ylin.makeRatio       = makeRatio
plot120_ylin.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,pt_xmax_lin + pt_xmax_lin/2]
plot120_ylin.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

plot120_ylog = Plot()
## inputs for stacked histograms
#plot120_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot120_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot120_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot120_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot120_ylog.keys            = keys
plot120_ylog.xtit            = "p_{T}(e#nu) [GeV] (Preselection)"
plot120_ylog.ytit            = "Events"
plot120_ylog.ylog            = "yes"
#plot120_ylog.rebin           = "var"
plot120_ylog.xmin            = pt_xmin
plot120_ylog.xmax            = pt_xmax_log
plot120_ylog.ymin            = pt_ymin
plot120_ylog.ymax            = pt_ymax_log
#plot120_ylog.lpos = "bottom-center"
plot120_ylog.name            = "Pt_Ele1MET_preselection_ylog"
plot120_ylog.addZUncBand     = zUncBand
plot120_ylog.makeRatio       = makeRatio
plot120_ylog.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,200,300,400,600,800,pt_xmax_log]
plot120_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- Pt_j1j2  ---
variableName = "Pt_j1j2"

plot121_ylin = Plot()
## inputs for stacked histograms
#plot121_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot121_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot121_ylin.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot121_ylin.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot121_ylin.keys            = keys
plot121_ylin.xtit            = "p_{T}(j1j2) [GeV] (Preselection)"
plot121_ylin.ytit            = "Events"
plot121_ylin.ylog            = "no"
#plot121_ylin.rebin           = "var"
plot121_ylin.xmin            = pt_xmin
plot121_ylin.xmax            = pt_xmax_lin + pt_xmax_lin/2
plot121_ylin.ymin            = pt_ymin
plot121_ylin.ymax            = pt_ymax_lin/2
#plot121_ylin.lpos = "bottom-center"
plot121_ylin.name            = "Pt_j1j2_preselection_ylin"
plot121_ylin.addZUncBand     = zUncBand
plot121_ylin.makeRatio       = makeRatio
plot121_ylin.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,pt_xmax_lin + pt_xmax_lin/2]
plot121_ylin.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

plot121_ylog = Plot()
## inputs for stacked histograms
#plot121_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot121_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot121_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot121_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot121_ylog.keys            = keys
plot121_ylog.xtit            = "p_{T}(j1j2) [GeV] (Preselection)"
plot121_ylog.ytit            = "Events"
plot121_ylog.ylog            = "yes"
#plot121_ylog.rebin           = "var"
plot121_ylog.xmin            = pt_xmin
plot121_ylog.xmax            = pt_xmax_log
plot121_ylog.ymin            = pt_ymin
plot121_ylog.ymax            = pt_ymax_log
#plot121_ylog.lpos = "bottom-center"
plot121_ylog.name            = "Pt_j1j2_preselection_ylog"
plot121_ylog.addZUncBand     = zUncBand
plot121_ylog.makeRatio       = makeRatio
plot121_ylog.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,200,300,400,600,800,pt_xmax_log]
plot121_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- sT_enujj  ---
variableName = "sT_enujj"

plot122_ylin = Plot()
## inputs for stacked histograms
#plot122_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot122_ylin.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot122_ylin.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot122_ylin.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot122_ylin.keys            = keys
plot122_ylin.xtit            = "S_{T} [GeV] (Preselection)"
plot122_ylin.ytit            = "Events"
plot122_ylin.ylog            = "no"
#plot122_ylin.rebin           = "var"
plot122_ylin.xmin            = 200
plot122_ylin.xmax            = pt_xmax_lin*3
plot122_ylin.ymin            = pt_ymin
plot122_ylin.ymax            = pt_ymax_lin/2
#plot122_ylin.lpos = "bottom-center"
plot122_ylin.name            = "sT_enujj_preselection_ylin"
plot122_ylin.addZUncBand     = zUncBand
plot122_ylin.makeRatio       = makeRatio
plot122_ylin.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

plot122_ylog = Plot()
## inputs for stacked histograms
#plot122_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot122_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot122_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot122_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot122_ylog.keys            = keys
plot122_ylog.xtit            = "S_{T} [GeV] (Preselection)"
plot122_ylog.ytit            = "Events"
plot122_ylog.ylog            = "yes"
#plot122_ylog.rebin           = "var"
plot122_ylog.xmin            = 200
plot122_ylog.xmax            = pt_xmax_log_high
plot122_ylog.ymin            = pt_ymin
plot122_ylog.ymax            = pt_ymax_log
#plot122_ylog.lpos = "bottom-center"
plot122_ylog.name            = "sT_enujj_preselection_ylog"
plot122_ylog.addZUncBand     = zUncBand
plot122_ylog.makeRatio       = makeRatio
plot122_ylog.xbins           = [200,220,240,260,280,320,340,360,380,400,420,460,500,550,
                               600,650,700,750,800,850,900,950,1000,1050,1100,1200,1300,1400,1500,1600,1700,1800,1900,pt_xmax_log_high]
plot122_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- nVertex ---
variableName = "nVertex"

plot123 = Plot()
## inputs for stacked histograms
#plot123.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot123.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot123.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot123.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot123.keys            = keys
plot123.xtit            = "Number of vertices (Preselection)"
plot123.ytit            = "Events"
plot123.ylog            = "yes"
plot123.rebin           = 1
plot123.xmin            = -0.5
plot123.xmax            = 30.5
plot123.ymin            = 0.1
plot123.ymax            = counter_ymax
#plot123.lpos = "bottom-center"
plot123.name            = "nVertex_preselection"
plot123.addZUncBand     = zUncBand
plot123.makeRatio       = makeRatio
plot123.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- nVertex_good ---
variableName = "nVertex_good"

plot124 = Plot()
## inputs for stacked histograms
#plot124.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot124.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot124.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot124.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot124.keys            = keys
plot124.xtit            = "Number of good-quality vertices (Preselection)"
plot124.ytit            = "Events"
plot124.ylog            = "yes"
plot124.rebin           = 1
plot124.xmin            = -0.5
plot124.xmax            = 30.5
plot124.ymin            = 0.1
plot124.ymax            = counter_ymax
#plot124.lpos = "bottom-center"
plot124.name            = "nVertex_good_preselection"
plot124.addZUncBand     = zUncBand
plot124.makeRatio       = makeRatio
plot124.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)



# ############################ Plots below to be done after full selection ######################

# ############# Extra plots #############
#if doExtraPlots:

#-----------------------------------------------------------------------------------


# list of plots to be plotted
plots  = [
          # Electrons
          plot0, plot1_ylin, plot1_ylog, plot2, plot3, plot4, plot5,
          # MET
          plot20_ylin, plot20_ylog, plot21, plot22_ylin, plot22_ylog, plot23_ylin,
          plot23_ylog, plot24_ylin, plot24_ylog, plot25_ylin, plot25_ylog,
          # Jets
          plot40, plot41, plot42_ylin, plot42_ylog, plot43, plot44,
          plot45_ylin, plot45_ylog, plot46_ylin, plot46_ylog, plot47, plot48, plot49_ylin, plot49_ylog,
          plot50_ylin, plot50_ylog, plot51, plot52_ylin, plot52_ylog,
          #Muons
          plot100,
          #Others
          plot120_ylin, plot120_ylog, plot121_ylin, plot121_ylog, plot122_ylin, plot122_ylog, plot123, plot124
          ]
#if doExtraPlots:
#  extra_plots = [plot_Vtxd0, plot_MissHits, plot_Dist, plot_DCotTheta
#                ]
#plots = plots + extra_plots

############# USER CODE - END ################################################
##############################################################################


#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()
c.Print("allPlots.ps[")
for plot in plots:
    plot.Draw("allPlots.ps")
c.Print("allPlots.ps]")

# Uncomment the 3 lines below to create (not rotated) pdf files,
# but keep them commented out in cvs since may slow down the execution
#print "Converting eps files into pdf files ..."
#for plot in plots:
#    system("convert "+plot.name+".eps "+plot.name+".pdf") # instead, uncomment this line to create pdf plots (a bit slow)

# create a file with list of eps and pdf files, to facilitate copying them to svn area for AN/PAS/PAPER
print "Creating file listOfEpsFiles.txt and listOfPdfFiles.txt ..."
system("rm listOfEpsFiles.txt listOfPdfFiles.txt")
for plot in plots:
    system("echo "+plot.name+".eps >> listOfEpsFiles.txt; echo "+plot.name+".pdf >> listOfPdfFiles.txt")
print "Use for example: scp `cat listOfPdfFiles.txt` pcuscms46:/home/santanas/Documents/CMSNotes/notes/AN-10-361/trunk/plots"



