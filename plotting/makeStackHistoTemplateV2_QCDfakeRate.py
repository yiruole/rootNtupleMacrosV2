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
            new_histo.SetBinContent( nbins, (new_histo.GetBinContent(nbins) + overflowBinContent)*(minBinWidth/binWidths[iter-1]) )
            new_histo.SetBinError( nbins, sqrt(new_histo.GetBinError(nbins)**2 + overflowBinError2)*(minBinWidth/binWidths[iter-1]) )
        else:
            new_histo.SetBinContent( nbins, new_histo.GetBinContent(nbins)*(minBinWidth/binWidths[iter-1]) )
            new_histo.SetBinError( nbins, new_histo.GetBinError(nbins)*(minBinWidth/binWidths[iter-1]) )
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
    #lint        = "3.0 pb^{-1}" # integrated luminosity of the sample ( example "10 pb^{-1}" )
    lint        = "33.0 pb^{-1}" # integrated luminosity of the sample ( example "10 pb^{-1}" )
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
        #             xlog may npot work         if (self.xlog     == "yes"):
        #             fPads1.SetLogx()
        if (self.ylog     == "yes"):
            fPads1.SetLogy()

        #-- legend
        hsize=0.20
        vsize=0.25
        if (self.lpos=="bottom-center"):
            xstart=0.35
            ystart=0.25
        elif(self.lpos=="top-left"):
            xstart=0.12
            ystart=0.63
        else:
            xstart=0.68
            ystart=0.63
        legend = TLegend(xstart, ystart, xstart+hsize, ystart+vsize)
        legend.SetFillColor(kWhite)
        legend.SetMargin(0.2)

        #-- loop over histograms (stacked)
        Nstacked = len(self.histosStack)
        #stackColorIndexes = [20,38,14,45,20,38,14,45]
        stackColorIndexes = [20,38,12,14,20,38,12,14]
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
            stack[iter].SetFillColor(  stackColorIndexes[iter])
#            stack[iter].SetMarkerColor(15+10*iter)
#            stack[iter].SetLineColor(  15+10*iter)
#            stack[iter].SetFillColor(  15+10*iter)
            legend.AddEntry(stack[iter], self.keysStack[Nstacked - iter - 1],"lf")
            #draw stack
            if iter==0:
                stack[iter].SetTitle("")
                stack[iter].GetXaxis().SetTitle(self.xtit)
                stack[iter].GetYaxis().SetTitle(self.ytit + " #times ("+ str(minBinW) + ")/(bin width)")
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
            else:
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
        dataLineIndexes = [1,2,3,1,2,3]
        for histo in self.histos:
            histo.SetMarkerStyle(dataColorIndexes[ih])
            histo.SetMarkerColor(dataColorIndexes[ih])
            histo.SetLineColor(  dataColorIndexes[ih])
            histo.SetLineStyle(  dataLineIndexes[ih])
            #            histo.SetMarkerStyle(20+2*ih)
            #            histo.SetMarkerColor(2+2*ih)
            #            histo.SetLineColor(  2+2*ih)
            legend.AddEntry(histo, self.keys[ih],"l")
            histo.Draw("HISTsame")
            ih=ih+1

        #-- plot data
        if(self.histodata!=""):
            self.histodata.SetMarkerStyle(20)
            legend.AddEntry(self.histodata, "data","p")
            self.histodata.Draw("psame")

        #-- draw label
        l = TLatex()
        l.SetTextAlign(12)
        l.SetTextSize(0.04)
        l.SetTextFont(62)
        l.SetNDC()
#        l.DrawLatex(xstart,ystart-0.05,"CMS Preliminary 2010")
#        l.DrawLatex(xstart,ystart-0.10,"L_{int} = " + self.lint)
        if (self.lpos=="bottom-center"):
            l.DrawLatex(0.35,0.20,"CMS 2010")
            l.DrawLatex(0.35,0.15,"L_{int} = " + self.lint)
        if (self.lpos=="top-left"):
            l.DrawLatex(xstart+hsize+0.02,ystart+vsize-0.03,"CMS 2010")
            l.DrawLatex(xstart+hsize+0.02,ystart+vsize-0.08,"L_{int} = " + self.lint)
        else:
            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.03,"CMS 2010")
            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.08,"L_{int} = " + self.lint)

        legend.Draw()
        canvas.Update()
        gPad.RedrawAxis()
        gPad.Modified()

        #-- 2nd pad (ratio)
        if(self.makeRatio==1):
            fPads2.cd()
            fPads2.SetLogy()
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
            #h_ratio1.GetYaxis().SetLimits(0.,2)
            h_ratio1.GetYaxis().SetRangeUser(0.1,100)
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
        #canvas.SaveAs(self.name + ".root","root")
        #canvas.SaveAs(self.name + ".pdf","pdf") # do not use this line because root creates rotated pdf plot - see end of the file instead
        canvas.Print(fileps)



############# DON'T NEED TO MODIFY ANYTHING HERE - END #######################
##############################################################################


##############################################################################
############# USER CODE - BEGING #############################################

#--- Input root file


#QCD fake rate: 3.0pb-1_QCD_fakeRate_run_lt144114_njet_More1_MET_lt35
#File_selection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/3.0pb-1_QCD_fakeRate_run_lt144114_njet_More1_MET_lt35/output_cutTable_QCD_fakeRate/analysisClass_QCD_fakeRate_plots.root")

#QCD fake rate: 33.0pb-1_QCD_fakeRate_run_gt144114_njet_More1_MET_lt35
File_selection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/33.0pb-1_QCD_fakeRate_run_gt144114_njet_More1_MET_lt35/output_cutTable_QCD_fakeRate/analysisClass_QCD_fakeRate_plots.root")

#QCD fake rate: 33.0pb-1_QCD_fakeRate_run_gt144114_njet_More1
#File_selection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/33.0pb-1_QCD_fakeRate_run_gt144114_njet_More1/output_cutTable_QCD_fakeRate/analysisClass_QCD_fakeRate_plots.root")

#### Common values for plots:
otherBkgsKey="Other Bkgs"
zUncBand="no"
makeRatio=1

pt_xmin=35
pt_xmax=500
pt_ymin=0.01
pt_ymax=1000000

eta_rebin=1
eta_ymin=0
eta_ymax=15000

#--- Final plots are defined here

# Simply add or remove plots and update the list plots = [plot0, plot1, ...] below

#histoBaseName = "histo1D__SAMPLE__cutHisto_allPreviousCuts________VARIABLE"
histoBaseName = "histo1D__SAMPLE__cutHisto_allCuts________________VARIABLE"
histoBaseName_userDef = "histo1D__SAMPLE__VARIABLE"

samplesForStackHistos = ["OTHERBKG","TTbar_Madgraph","PhotJetPt30","WJetAlpgen"]
keysStack =             [otherBkgsKey,"ttbar", "\\gamma + jet" ,"W/W* + jets"]

samplesForHistos = ["LQenujj_M200", "LQenujj_M250","LQenujj_M300"]
keys             = ["LQ e\\nujj M200","LQ e\\nujj M250","LQ e\\nujj M300"]

sampleForDataHisto = "DATA"

#--- nEle_PtCut_IDISO_noOvrlp ---
variableName = "nEle_PtCut_IDISO_noOvrlp"

plot0 = Plot()
## inputs for stacked histograms
plot0.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_selection)
plot0.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot0.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_selection)
plot0.keys            = keys
plot0.xtit            = "Number of electrons"
plot0.ytit            = "Number of events"
plot0.ylog            = "yes"
plot0.rebin           = 1
plot0.xmin            = -0.5
plot0.xmax            = 6.5
plot0.ymin            = 0.0001
plot0.ymax            = 100000000
#plot0.lpos = "bottom-center"
plot0.name            = "nEle_QCDfakeRate"
plot0.addZUncBand     = zUncBand
plot0.makeRatio       = makeRatio
plot0.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_selection)


#--- h1_ElePt  ---
variableName = "h1_ElePt"

plot1 = Plot()
## inputs for stacked histograms
plot1.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_selection)
plot1.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot1.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_selection)
plot1.keys            = keys
plot1.xtit            = "pT 1st electron (GeV) - barrel + endcap"
plot1.ytit            = "Number of events"
plot1.ylog            = "yes"
plot1.rebin           = "var"
plot1.xmin            = pt_xmin
plot1.xmax            = pt_xmax
plot1.ymin            = pt_ymin
plot1.ymax            = pt_ymax
#plot1.lpos = "bottom-center"
plot1.name            = "pT1stEle_QCDfakeRate_bar_end"
plot1.addZUncBand     = zUncBand
plot1.makeRatio       = makeRatio
plot1.xbins           = [35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,200,220,240,260,280,300,350,400,500]
plot1.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_selection)

#--- h1_ElePt_barrel1  ---
variableName = "h1_ElePt_barrel1"

plot2 = Plot()
## inputs for stacked histograms
plot2.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_selection)
plot2.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot2.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_selection)
plot2.keys            = keys
plot2.xtit            = "pT 1st electron (GeV) - barrel 1"
plot2.ytit            = "Number of events"
plot2.ylog            = "yes"
plot2.rebin           = "var"
plot2.xmin            = pt_xmin
plot2.xmax            = pt_xmax
plot2.ymin            = pt_ymin
plot2.ymax            = pt_ymax
#plot2.lpos = "bottom-center"
plot2.name            = "pT1stEle_QCDfakeRate_bar1"
plot2.addZUncBand     = zUncBand
plot2.makeRatio       = makeRatio
plot2.xbins           = [35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,200,220,240,260,280,300,350,400,500]
plot2.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_selection)

#--- h1_ElePt_barrel2  ---
variableName = "h1_ElePt_barrel2"

plot3 = Plot()
## inputs for stacked histograms
plot3.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_selection)
plot3.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot3.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_selection)
plot3.keys            = keys
plot3.xtit            = "pT 1st electron (GeV) - barrel 2"
plot3.ytit            = "Number of events"
plot3.ylog            = "yes"
plot3.rebin           = "var"
plot3.xmin            = pt_xmin
plot3.xmax            = pt_xmax
plot3.ymin            = pt_ymin
plot3.ymax            = pt_ymax
#plot3.lpos = "bottom-center"
plot3.name            = "pT1stEle_QCDfakeRate_bar2"
plot3.addZUncBand     = zUncBand
plot3.makeRatio       = makeRatio
plot3.xbins           = [35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,200,220,240,260,280,300,350,400,500]
plot3.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_selection)

#--- h1_ElePt_endcap1  ---
variableName = "h1_ElePt_endcap1"

plot4 = Plot()
## inputs for stacked histograms
plot4.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_selection)
plot4.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot4.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_selection)
plot4.keys            = keys
plot4.xtit            = "pT 1st electron (GeV) - endcap 1"
plot4.ytit            = "Number of events"
plot4.ylog            = "yes"
plot4.rebin           = "var"
plot4.xmin            = pt_xmin
plot4.xmax            = pt_xmax
plot4.ymin            = pt_ymin
plot4.ymax            = pt_ymax
#plot4.lpos = "bottom-center"
plot4.name            = "pT1stEle_QCDfakeRate_end1"
plot4.addZUncBand     = zUncBand
plot4.makeRatio       = makeRatio
plot4.xbins           = [35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,200,220,240,260,280,300,350,400,500]
plot4.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_selection)

#--- h1_ElePt_endcap2  ---
variableName = "h1_ElePt_endcap2"

plot5 = Plot()
## inputs for stacked histograms
plot5.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_selection)
plot5.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot5.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_selection)
plot5.keys            = keys
plot5.xtit            = "pT 1st electron (GeV) - endcap 2"
plot5.ytit            = "Number of events"
plot5.ylog            = "yes"
plot5.rebin           = "var"
plot5.xmin            = pt_xmin
plot5.xmax            = pt_xmax
plot5.ymin            = pt_ymin
plot5.ymax            = pt_ymax
#plot5.lpos = "bottom-center"
plot5.name            = "pT1stEle_QCDfakeRate_end2"
plot5.addZUncBand     = zUncBand
plot5.makeRatio       = makeRatio
plot5.xbins           = [35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,200,220,240,260,280,300,350,400,500]
plot5.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_selection)

#--- h1_EleEta ---
variableName = "h1_EleEta"

plot6 = Plot()
## inputs for stacked histograms
plot6.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_selection)
plot6.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot6.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_selection)
plot6.keys            = keys
plot6.xtit            = "#eta 1st electron"
plot6.ytit            = "Number of events"
plot6.rebin           = eta_rebin
plot6.xmin            = -3
plot6.xmax            = 3
plot6.ymin            = eta_ymin
plot6.ymax            = eta_ymax
#plot6.lpos            = "top-left"
plot6.name            = "Eta1stEle_QCDfakeRate"
plot6.addZUncBand     = zUncBand
plot6.makeRatio       = makeRatio
plot6.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_selection)

#--- Phi1stEle ---
variableName = "Phi1stEle"

plot7 = Plot()
## inputs for stacked histograms
plot7.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_selection)
plot7.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot7.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_selection)
plot7.keys            = keys
plot7.xtit            = "#phi 1st electron (rad.)"
plot7.ytit            = "Number of events"
plot7.rebin           = eta_rebin
#plot7.xmin            = -3.15
#plot7.xmax            = 3.15
plot7.ymin            = eta_ymin
plot7.ymax            = eta_ymax
#plot7.lpos            = "top-left"
plot7.name            = "Phi1stEle_QCDfakeRate"
plot7.addZUncBand     = zUncBand
plot7.makeRatio       = makeRatio
plot7.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_selection)


#--- Charge1stEle_PAS  ---
variableName = "Charge1stEle"

plot8 = Plot()
## inputs for stacked histograms
plot8.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_selection)
plot8.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot8.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_selection)
plot8.keys            = keys
plot8.xtit            = "Charge of reconstructed electron"
plot8.ytit            = "Number of events"
plot8.ylog            = "yes"
plot8.rebin           = 1
#plot8.xmin            = -1.0001
#plot8.xmax            = 1.0001
plot8.ymin            = 1
plot8.ymax            = pt_ymax*10
#plot8.lpos = "bottom-center"
plot8.name            = "charge1stEle_QCDfakeRate"
plot8.addZUncBand     = zUncBand
plot8.makeRatio       = makeRatio
plot8.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_selection)


#--- nJet_WithJetEtaCut ---
variableName = "nJet_PtCut_noOvrlp_ID"

plot9 = Plot()
## inputs for stacked histograms
plot9.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_selection)
plot9.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot9.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_selection)
plot9.keys            = keys
plot9.xtit            = "Number of jets"
plot9.ytit            = "Number of events"
plot9.ylog            = "yes"
plot9.rebin           = 1
plot9.xmin            = -0.5
plot9.xmax            = 11.5
plot9.ymin            = 0.01
plot9.ymax            = 1000000
#plot9.lpos = "bottom-center"
plot9.name            = "nJet_QCDfakeRate"
plot9.addZUncBand     = zUncBand
plot9.makeRatio       = makeRatio
plot9.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_selection)


#--- MET_PAS  ---
variableName = "MET"

plot10 = Plot()
## inputs for stacked histograms
plot10.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_selection)
plot10.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot10.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_selection)
plot10.keys            = keys
plot10.xtit            = "pfMET (GeV)"
plot10.ytit            = "Number of events"
plot10.ylog            = "yes"
plot10.rebin           = "var"
plot10.xmin            = 0
plot10.xmax            = 300
plot10.ymin            = 0.01
plot10.ymax            = 100000
#plot10.lpos = "bottom-center"
plot10.name            = "MET_QCDfakeRate"
plot10.addZUncBand     = zUncBand
plot10.makeRatio       = makeRatio
plot10.xbins           = [0,10,20,30,40,50,70,90,110,130,150,200,250,300]
plot10.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_selection)


#--- MTenu ---
variableName = "MTenu"

plot11 = Plot()
## inputs for stacked histograms
plot11.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_selection)
plot11.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot11.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_selection)
plot11.keys            = keys
plot11.xtit            = "M_{T}(e\\nu) (GeV)"
plot11.ytit            = "Number of events"
# plot11.ylog            = "yes"
# plot11.rebin           = 1
# plot11.ymin            = 0.00000001
# plot11.ymax            = 20
plot11.ylog            = "no"
plot11.rebin           = 1
plot11.xmin            = 0
plot11.xmax            = 400
plot11.ymin            = 0
plot11.ymax            = 100000
#plot11.lpos = "bottom-center"
plot11.name            = "MTenu_QCDfakeRate_ylin"
plot11.addZUncBand     = zUncBand
plot11.makeRatio       = makeRatio
plot11.xbins           = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,200,300,400,500]
plot11.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_selection)

plot11_ylog = Plot()
## inputs for stacked histograms
plot11_ylog.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_selection)
plot11_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot11_ylog.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_selection)
plot11_ylog.keys            = keys
plot11_ylog.xtit            = "M_{T}(e\\nu) (GeV)"
plot11_ylog.ytit            = "Number of events"
plot11_ylog.ylog            = "yes"
plot11_ylog.rebin           = "var" # don't change it (since a rebinning is already applied above on the same histo)
plot11_ylog.xmin            = 0
plot11_ylog.xmax            = 1000
plot11_ylog.ymin            = 0.01
plot11_ylog.ymax            = 100000
#plot11_ylog.lpos = "bottom-center"
plot11_ylog.name            = "MTenu_QCDfakeRate"
plot11_ylog.addZUncBand     = zUncBand
plot11_ylog.makeRatio       = makeRatio
plot11_ylog.xbins           = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,200,300,500]
plot11_ylog.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_selection)


#--- h1_Ele_EtOverPt ---
variableName = "h1_Ele_EtOverPt"

plot12 = Plot()
## inputs for stacked histograms
plot12.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_selection)
plot12.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot12.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_selection)
plot12.keys            = keys
plot12.xtit            = "Et/Pt of reconstructed electron - barrel + endcap"
plot12.ytit            = "Number of events"
plot12.ylog            = "no"
plot12.rebin           = 10
plot12.xmin            = 0
plot12.xmax            = 10
plot12.ymin            = 0.
plot12.ymax            = 15000
#plot12.lpos            = "top-left"
plot12.name            = "EtOverPtEle_QCDfakeRate_ylin"
plot12.addZUncBand     = zUncBand
plot12.makeRatio       = makeRatio
plot12.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_selection)

plot12_ylog = Plot()
## inputs for stacked histograms
plot12_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_selection)
plot12_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot12_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_selection)
plot12_ylog.keys            = keys
plot12_ylog.xtit            = "Et/Pt of reconstructed electron - barrel + endcap"
plot12_ylog.ytit            = "Number of events"
plot12_ylog.ylog            = "yes"
plot12_ylog.rebin           = 10
plot12_ylog.xmin            = 0
plot12_ylog.xmax            = 10
plot12_ylog.ymin            = 0.1
plot12_ylog.ymax            = 100000
#plot12_ylog.lpos            = "top-left"
plot12_ylog.name            = "EtOverPtEle_QCDfakeRate"
plot12_ylog.addZUncBand     = zUncBand
plot12_ylog.makeRatio       = makeRatio
plot12_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_selection)


#--- h1_Ele_EtOverPt_barrel1 ---
variableName = "h1_Ele_EtOverPt_barrel1"

plot13 = Plot()
## inputs for stacked histograms
plot13.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_selection)
plot13.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot13.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_selection)
plot13.keys            = keys
plot13.xtit            = "Et/Pt of reconstructed electron - barrel1"
plot13.ytit            = "Number of events"
plot13.ylog            = "no"
plot13.rebin           = 10
plot13.xmin            = 0
plot13.xmax            = 10
plot13.ymin            = 0.
plot13.ymax            = 10000
#plot13.lpos            = "top-left"
plot13.name            = "EtOverPtEle_QCDfakeRate_ylin_bar1"
plot13.addZUncBand     = zUncBand
plot13.makeRatio       = makeRatio
plot13.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_selection)

#--- h1_Ele_EtOverPt_barrel2 ---
variableName = "h1_Ele_EtOverPt_barrel2"

plot14 = Plot()
## inputs for stacked histograms
plot14.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_selection)
plot14.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot14.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_selection)
plot14.keys            = keys
plot14.xtit            = "Et/Pt of reconstructed electron - barrel2"
plot14.ytit            = "Number of events"
plot14.ylog            = "no"
plot14.rebin           = 10
plot14.xmin            = 0
plot14.xmax            = 10
plot14.ymin            = 0.
plot14.ymax            = 5000
#plot14.lpos            = "top-left"
plot14.name            = "EtOverPtEle_QCDfakeRate_ylin_bar2"
plot14.addZUncBand     = zUncBand
plot14.makeRatio       = makeRatio
plot14.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_selection)

#--- h1_Ele_EtOverPt_endcap1 ---
variableName = "h1_Ele_EtOverPt_endcap1"

plot15 = Plot()
## inputs for stacked histograms
plot15.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_selection)
plot15.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot15.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_selection)
plot15.keys            = keys
plot15.xtit            = "Et/Pt of reconstructed electron - endcap1"
plot15.ytit            = "Number of events"
plot15.ylog            = "no"
plot15.rebin           = 10
plot15.xmin            = 0
plot15.xmax            = 10
plot15.ymin            = 0.
plot15.ymax            = 5000
#plot15.lpos            = "top-left"
plot15.name            = "EtOverPtEle_QCDfakeRate_ylin_end1"
plot15.addZUncBand     = zUncBand
plot15.makeRatio       = makeRatio
plot15.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_selection)

#--- h1_Ele_EtOverPt_endcap2 ---
variableName = "h1_Ele_EtOverPt_endcap2"

plot16 = Plot()
## inputs for stacked histograms
plot16.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_selection)
plot16.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot16.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_selection)
plot16.keys            = keys
plot16.xtit            = "Et/Pt of reconstructed electron - endcap2"
plot16.ytit            = "Number of events"
plot16.ylog            = "no"
plot16.rebin           = 10
plot16.xmin            = 0
plot16.xmax            = 10
plot16.ymin            = 0.
plot16.ymax            = 5000
#plot16.lpos            = "top-left"
plot16.name            = "EtOverPtEle_QCDfakeRate_ylin_end2"
plot16.addZUncBand     = zUncBand
plot16.makeRatio       = makeRatio
plot16.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_selection)

#--- mDeltaPhiMETEle ---
variableName = "mDeltaPhiMETEle"

plot17 = Plot()
## inputs for stacked histograms
plot17.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_selection)
plot17.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot17.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_selection)
plot17.keys            = keys
plot17.xtit            = "#Delta#phi(MET,e) (rad.)"
plot17.ytit            = "Number of events"
#plot17.xlog            = "yes"
plot17.ylog            = "yes"
plot17.rebin           = 5
#plot17.xmin            = 0
#plot17.xmax            = 3.14
plot17.ymin            = 0.1
plot17.ymax            = 10000000
#plot17.lpos = "bottom-center"
plot17.name            = "mDeltaPhiMETEle_QCDfakeRate"
plot17.addZUncBand     = zUncBand
plot17.makeRatio       = makeRatio
plot17.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_selection)


#--- minDRscjets ---
variableName = "minDRscjets"

plot18 = Plot()
## inputs for stacked histograms
plot18.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_selection)
plot18.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot18.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_selection)
plot18.keys            = keys
plot18.xtit            = "min#Delta R(e,jets)"
plot18.ytit            = "Number of events"
plot18.ylog            = "yes"
plot18.rebin           = 2
plot18.xmin            = 0
plot18.xmax            = 7
plot18.ymin            = 0.01
plot18.ymax            = 100000
#plot18.lpos = "bottom-center"
plot18.name            = "minDRej_QCDfakeRate"
plot18.addZUncBand     = zUncBand
plot18.makeRatio       = makeRatio
#plot18.xbins           = [30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot18.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_selection)


#--- Pt1stJet ---
variableName = "Pt1stJet"

plot19 = Plot()
## inputs for stacked histograms
plot19.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_selection)
plot19.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot19.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_selection)
plot19.keys            = keys
plot19.xtit            = "pT 1st jet (GeV)"
plot19.ytit            = "Number of events"
plot19.ylog            = "yes"
plot19.rebin           = 1
plot19.xmin            = pt_xmin
plot19.xmax            = pt_xmax
plot19.ymin            = pt_ymin
plot19.ymax            = pt_ymax
#plot19.lpos = "bottom-center"
plot19.name            = "Pt1stJet_QCDfakeRate"
plot19.addZUncBand     = zUncBand
plot19.makeRatio       = makeRatio
plot19.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,200,300,400,800]
plot19.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_selection)

#--- Eta1stJet ---
variableName = "Eta1stJet"

plot20 = Plot()
## inputs for stacked histograms
plot20.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_selection)
plot20.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot20.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_selection)
plot20.keys            = keys
plot20.xtit            = "#eta 1st jet"
plot20.ytit            = "Number of events"
plot20.rebin           = eta_rebin
plot20.ymin            = eta_ymin
plot20.ymax            = eta_ymax
#plot20.lpos            = "top-left"
plot20.name            = "Eta1stJet_QCDfakeRate"
plot20.addZUncBand     = zUncBand
plot20.makeRatio       = makeRatio
plot20.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_selection)



# #--- mDeltaPhiMET1stJet_PAS ---
# variableName = "mDeltaPhiMET1stJet_PAS"

# plot12 = Plot()
# ## inputs for stacked histograms
# plot12.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_selection)
# plot12.keysStack       = keysStack
# ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot12.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_selection)
# plot12.keys            = keys
# plot12.xtit            = "#Delta#phi(MET,1^{st} jet) (rad.)"
# plot12.ytit            = "Number of events"
# #plot12.xlog            = "yes"
# plot12.ylog            = "yes"
# plot12.rebin           = 6
# #plot12.xmin            = 0
# #plot12.xmax            = 3.146
# plot12.ymin            = 0.1
# plot12.ymax            = 1000000
# #plot12.lpos = "bottom-center"
# plot12.name            = "mDeltaPhiMET1stJet_QCDfakeRate"
# plot12.addZUncBand     = zUncBand
# plot12.makeRatio       = makeRatio
# plot12.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_selection)




# #--- METPhi_PAS  ---
# variableName = "METPhi_PAS"

# plot22 = Plot()
# ## inputs for stacked histograms
# plot22.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_selection)
# plot22.keysStack       = keysStack
# ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot22.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_selection)
# plot22.keys            = keys
# plot22.xtit            = "pfMET #phi (rad.)"
# plot22.ytit            = "Number of events"
# plot22.rebin           = eta_rebin*2
# plot22.ymin            = eta_ymin
# plot22.ymax            = eta_ymax
# #plot22.lpos = "bottom-center"
# plot22.name            = "METPhi_QCDfakeRate"
# plot22.addZUncBand     = zUncBand
# plot22.makeRatio       = makeRatio
# plot22.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_selection)




#-----------------------------------------------------------------------------------


# List of plots to be plotted
plots  = [plot0, plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8,
          plot9, plot10, plot11, plot11_ylog,
          plot12, plot12_ylog, plot13, plot14, plot15, plot16, plot17,
          plot18, plot19, plot20
          ]



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

