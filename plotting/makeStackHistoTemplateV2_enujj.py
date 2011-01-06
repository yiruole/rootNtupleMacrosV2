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
    lint        = "36.0 pb^{-1}" # integrated luminosity of the sample ( example "10 pb^{-1}" )
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
            h_ratio1.GetYaxis().SetRangeUser(0.1,10)
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

File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_fullntuples_MET45/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/36.05pb-1_sT_presel_250_Zrescale1.20/analysisClass_enujjSample_plots.root")
##File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/33.2pb-1_sT_presel_250_Zrescale1.20_Wrescale1.06_extraPlotsDec9/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/f/ferencek/scratch0/LQ/CMSSW_3_5_7/test/Leptoquarks/rootNtupleAnalyzerV2/data/output/36.05pb-1_sT_presel_250_Zrescale1.20/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/f/ferencek/scratch0/LQ/CMSSW_3_5_7/test/Leptoquarks/rootNtupleAnalyzerV2/data/output/36.05pb-1_sT_presel_250_MET_presel_45_Zrescale1.20_Wrescale1.19/analysisClass_enujjSample_plots.root")

File_selection    = File_preselection

UseQCDFromData    = 1 # always put an existing file under File_QCD (otherwise the code will crash)

File_QCD          = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_fullntuples_MET45/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")
#File_QCD          = GetFile("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/35.84pb-1_QCD_UseHLTPrescales_sT_presel_250/analysisClass_enujjSample_QCD_plots.root")
##File_QCD          = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/6.1pb-1_QCD_HLT30_sT_presel_250_extraPlotsDec9/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")
#File_QCD          = GetFile("/afs/cern.ch/user/f/ferencek/scratch0/LQ/CMSSW_3_5_7/test/Leptoquarks/rootNtupleAnalyzerV2/data/output/35.84pb-1_QCD_UseHLTPrescales_sT_presel_250/analysisClass_enujjSample_QCD_plots.root")
#File_QCD          = GetFile("/afs/cern.ch/user/f/ferencek/scratch0/LQ/CMSSW_3_5_7/test/Leptoquarks/rootNtupleAnalyzerV2/data/output/35.84pb-1_QCD_UseHLTPrescales_sT_presel_250_MET_presel_45/analysisClass_enujjSample_QCD_plots.root")

QCDscaleFactor    = 1 # no need to rescale anymore since we are using the HLT prescales (36/35.84 can be ignored)

#### Common values for plots:
#otherBkgsKey="QCD, single top, VV+jets, W/W*+jets"
otherBkgsKey="Other Bkgs"
zUncBand="no"
makeRatio=1

pt_xmin=0
pt_xmax=800
pt_ymin=0.01
pt_ymax=1000

eta_rebin=2
eta_ymin=0
eta_ymax=100

#--- Final plots are defined here

# Simply add or remove plots and update the list plots = [plot0, plot1, ...] below

histoBaseName = "histo1D__SAMPLE__cutHisto_allPreviousCuts________VARIABLE"
histoBaseName_userDef = "histo1D__SAMPLE__VARIABLE"

samplesForStackHistosQCD = ["DATA"]
samplesForStackHistos = ["OTHERBKG","TTbar_Madgraph","WJetAlpgen"]
keysStack =             ["QCD",otherBkgsKey,"ttbar", "W/W* + jets"]

samplesForHistos = ["LQenujj_M200", "LQenujj_M250","LQenujj_M300"]
keys             = ["LQ e\\nujj M200","LQ e\\nujj M250","LQ e\\nujj M300"]

sampleForDataHisto = "DATA"

#--- nEle_PtCut_IDISO_noOvrlp ---
variableName = "nEle_PtCut_IDISO_noOvrlp"

#h_nEle_LQenujj_M100 = GetHisto("histo1D__LQenujj_M100__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
#h_nEle_LQenujj_M200 = GetHisto("histo1D__LQenujj_M200__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
#h_nEle_LQenujj_M300 = GetHisto("histo1D__LQenujj_M300__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
#h_nEle_LQenujj_M400 = GetHisto("histo1D__LQenujj_M400__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
#h_nEle_LQenujj_M500 = GetHisto("histo1D__LQenujj_M500__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
#h_nEle_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
#h_nEle_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
#h_nEle_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
#h_nEle_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
#h_nEle_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
#h_nEle_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
#h_nEle_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
#h_nEle_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()
#h_nEle_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________nEle_PtCut_IDISO_noOvrlp", File_preselection).Clone()

plot0 = Plot()
## inputs for stacked histograms
plot0.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot0.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot0.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
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
plot0.name            = "nEle_allPreviousCuts"
plot0.addZUncBand     = zUncBand
plot0.makeRatio       = makeRatio
plot0.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Pt1stEle_PAS  ---
variableName = "Pt1stEle_PAS"

plot1 = Plot()
## inputs for stacked histograms
plot1.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot1.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot1.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot1.keys            = keys
plot1.xtit            = "pT 1st electron (GeV)"
plot1.ytit            = "Number of events"
plot1.ylog            = "yes"
plot1.rebin           = "var"
plot1.xmin            = pt_xmin
plot1.xmax            = pt_xmax
plot1.ymin            = pt_ymin
plot1.ymax            = pt_ymax
#plot1.lpos = "bottom-center"
plot1.name            = "pT1stEle_allPreviousCuts"
plot1.addZUncBand     = zUncBand
plot1.makeRatio       = makeRatio
plot1.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,200,300,400,800]
plot1.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- Eta1stEle_PAS  ---
variableName = "Eta1stEle_PAS"

plot2 = Plot()
## inputs for stacked histograms
plot2.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot2.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot2.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot2.keys            = keys
plot2.xtit            = "#eta 1st electron"
plot2.ytit            = "Number of events"
plot2.rebin           = eta_rebin
plot2.xmin            = -3
plot2.xmax            = 3
plot2.ymin            = eta_ymin
plot2.ymax            = eta_ymax
#plot2.lpos            = "top-left"
plot2.name            = "Eta1stEle_allPreviousCuts"
plot2.addZUncBand     = zUncBand
plot2.makeRatio       = makeRatio
plot2.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- Charge1stEle_PAS  ---
variableName = "Charge1stEle_PAS"

plot_after2 = Plot()
## inputs for stacked histograms
plot_after2.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot_after2.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot_after2.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot_after2.keys            = keys
plot_after2.xtit            = "Charge of reconstructed electron"
plot_after2.ytit            = "Number of events"
plot_after2.ylog            = "yes"
plot_after2.rebin           = 1
#plot_after2.xmin            = -1.0001
#plot_after2.xmax            = 1.0001
plot_after2.ymin            = 1
plot_after2.ymax            = pt_ymax*50
#plot_after2.lpos = "bottom-center"
plot_after2.name            = "charge1stEle_allPreviousCuts"
plot_after2.addZUncBand     = zUncBand
plot_after2.makeRatio       = makeRatio
plot_after2.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- MET_PAS  ---
variableName = "MET_PAS"

plot3 = Plot()
## inputs for stacked histograms
plot3.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot3.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot3.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot3.keys            = keys
plot3.xtit            = "pfMET (GeV)"
plot3.ytit            = "Number of events"
plot3.ylog            = "yes"
plot3.rebin           = "var"
plot3.xmin            = 0
plot3.xmax            = 500
plot3.ymin            = 0.01
plot3.ymax            = 5000
#plot3.lpos = "bottom-center"
plot3.name            = "MET_allPreviousCuts"
plot3.addZUncBand     = zUncBand
plot3.makeRatio       = makeRatio
plot3.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,200,300,400,500]
plot3.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- minMETPt1stEle_PAS ---
variableName = "minMETPt1stEle_PAS"

plot4 = Plot()
## inputs for stacked histograms
plot4.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot4.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot4.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot4.keys            = keys
plot4.xtit            = "min(pT 1st electron,pfMET) (GeV)"
plot4.ytit            = "Number of events"
plot4.ylog            = "yes"
plot4.rebin           = "var"
plot4.xmin            = 0
plot4.xmax            = 300
plot4.ymin            = 0.01
plot4.ymax            = 5000
#plot4.lpos = "bottom-center"
plot4.name            = "minMETPt1stEle_allPreviousCuts"
plot4.addZUncBand     = zUncBand
plot4.makeRatio       = makeRatio
plot4.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,200,300]
plot4.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- nJet_WithJetEtaCut ---
variableName = "nJet_WithJetEtaCut"

plot5 = Plot()
## inputs for stacked histograms
plot5.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot5.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot5.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot5.keys            = keys
plot5.xtit            = "Number of jets"
plot5.ytit            = "Number of events"
plot5.ylog            = "yes"
plot5.rebin           = 1
plot5.xmin            = -0.5
plot5.xmax            = 11.5
plot5.ymin            = 0.01
plot5.ymax            = 1000000
#plot5.lpos = "bottom-center"
plot5.name            = "nJet_allPreviousCuts"
plot5.addZUncBand     = zUncBand
plot5.makeRatio       = makeRatio
plot5.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- Pt1stJet_PAS ---
variableName = "Pt1stJet_PAS"

plot6 = Plot()
## inputs for stacked histograms
plot6.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot6.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot6.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot6.keys            = keys
plot6.xtit            = "pT 1st jet (GeV)"
plot6.ytit            = "Number of events"
plot6.ylog            = "yes"
plot6.rebin           = "var"
plot6.xmin            = pt_xmin
plot6.xmax            = pt_xmax
plot6.ymin            = pt_ymin
plot6.ymax            = pt_ymax
#plot6.lpos = "bottom-center"
plot6.name            = "Pt1stJet_allPreviousCuts"
plot6.addZUncBand     = zUncBand
plot6.makeRatio       = makeRatio
plot6.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,200,300,400,800]
plot6.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- Eta1stJet_PAS ---
variableName = "Eta1stJet_PAS"

plot7 = Plot()
## inputs for stacked histograms
plot7.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot7.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot7.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot7.keys            = keys
plot7.xtit            = "#eta 1st jet"
plot7.ytit            = "Number of events"
plot7.rebin           = eta_rebin
plot7.ymin            = eta_ymin
plot7.ymax            = eta_ymax
#plot7.lpos            = "top-left"
plot7.name            = "Eta1stJet_allPreviousCuts"
plot7.addZUncBand     = zUncBand
plot7.makeRatio       = makeRatio
plot7.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- Pt2ndJet_PAS ---
variableName = "Pt2ndJet_PAS"

plot8 = Plot()
## inputs for stacked histograms
plot8.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot8.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot8.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot8.keys            = keys
plot8.xtit            = "pT 2nd jet (GeV)"
plot8.ytit            = "Number of events"
plot8.ylog            = "yes"
plot8.rebin           = "var"
plot8.xmin            = pt_xmin
plot8.xmax            = pt_xmax
plot8.ymin            = pt_ymin
plot8.ymax            = pt_ymax
#plot8.lpos = "bottom-center"
plot8.name            = "Pt2ndJet_allPreviousCuts"
plot8.addZUncBand     = zUncBand
plot8.makeRatio       = makeRatio
plot8.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,200,300,400,800]
plot8.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- Eta2ndJet_PAS ---
variableName = "Eta2ndJet_PAS"

plot9 = Plot()
## inputs for stacked histograms
plot9.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot9.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot9.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot9.keys            = keys
plot9.xtit            = "#eta 2nd jet"
plot9.ytit            = "Number of events"
plot9.rebin           = eta_rebin
plot9.ymin            = eta_ymin
plot9.ymax            = eta_ymax
#plot9.lpos            = "top-left"
plot9.name            = "Eta2ndJet_allPreviousCuts"
plot9.addZUncBand     = zUncBand
plot9.makeRatio       = makeRatio
plot9.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Ptjets ---
variableNames = ["Pt1stJet_PAS","Pt2ndJet_PAS"]

plot6and8 = Plot()
## inputs for stacked histograms
plot6and8.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistosQCD, variableNames, File_QCD, QCDscaleFactor) + generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot6and8.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot6and8.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot6and8.keys            = keys
plot6and8.xtit            = "pT jets (GeV)"
plot6and8.ytit            = "Number of events"
plot6and8.ylog            = "yes"
plot6and8.rebin           = "var"
plot6and8.xmin            = pt_xmin
plot6and8.xmax            = pt_xmax
plot6and8.ymin            = pt_ymin
plot6and8.ymax            = pt_ymax
#plot6and8.lpos = "bottom-center"
plot6and8.name            = "PtJets_allPreviousCuts"
plot6and8.addZUncBand     = zUncBand
plot6and8.makeRatio       = makeRatio
plot6and8.xbins           = [0,10,20,30,40,50,60,70,80,90,100,125,150,175,200,300,400,800]
plot6and8.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)


#--- Etajets ---
variableNames = ["Eta1stJet_PAS","Eta2ndJet_PAS"]

plot7and9 = Plot()
## inputs for stacked histograms
plot7and9.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistosQCD, variableNames, File_QCD, QCDscaleFactor) + generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot7and9.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot7and9.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot7and9.keys            = keys
plot7and9.xtit            = "#eta jets"
plot7and9.ytit            = "Number of events"
plot7and9.rebin           = eta_rebin
plot7and9.ymin            = eta_ymin
plot7and9.ymax            = eta_ymax*2
#plot7and9.lpos            = "top-left"
plot7and9.name            = "EtaJets_allPreviousCuts"
plot7and9.addZUncBand     = zUncBand
plot7and9.makeRatio       = makeRatio
plot7and9.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)


#--- nMuon_PtCut_IDISO_PAS ---
variableName = "nMuon_PtCut_IDISO_PAS"

plot10 = Plot()
## inputs for stacked histograms
plot10.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot10.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot10.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot10.keys            = keys
plot10.xtit            = "Number of muons"
plot10.ytit            = "Number of events"
plot10.ylog            = "yes"
plot10.rebin           = 1
plot10.xmin            = -0.5
plot10.xmax            = 4.5
plot10.ymin            = 0.01
plot10.ymax            = 10000
#plot10.lpos = "bottom-center"
plot10.name            = "nMuon_allPreviousCuts"
plot10.addZUncBand     = zUncBand
plot10.makeRatio       = makeRatio
plot10.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- mDeltaPhiMETEle_PAS ---
variableName = "mDeltaPhiMETEle_PAS"

plot11 = Plot()
## inputs for stacked histograms
plot11.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot11.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot11.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot11.keys            = keys
plot11.xtit            = "#Delta#phi(MET,e) (rad.)"
plot11.ytit            = "Number of events"
#plot11.xlog            = "yes"
plot11.ylog            = "yes"
plot11.rebin           = 5
#plot11.xmin            = 0
#plot11.xmax            = 3.14
plot11.ymin            = 0.1
plot11.ymax            = 1000000
#plot11.lpos = "bottom-center"
plot11.name            = "mDeltaPhiMETEle_allPreviousCuts"
plot11.addZUncBand     = zUncBand
plot11.makeRatio       = makeRatio
plot11.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- mDeltaPhiMET1stJet_PAS ---
variableName = "mDeltaPhiMET1stJet_PAS"

plot12 = Plot()
## inputs for stacked histograms
plot12.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot12.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot12.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot12.keys            = keys
plot12.xtit            = "#Delta#phi(MET,1^{st} jet) (rad.)"
plot12.ytit            = "Number of events"
#plot12.xlog            = "yes"
plot12.ylog            = "yes"
plot12.rebin           = 5
#plot12.xmin            = 0
#plot12.xmax            = 3.146
plot12.ymin            = 0.1
plot12.ymax            = 1000000
#plot12.lpos = "bottom-center"
plot12.name            = "mDeltaPhiMET1stJet_allPreviousCuts"
plot12.addZUncBand     = zUncBand
plot12.makeRatio       = makeRatio
plot12.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- mDeltaPhiMET2ndJet_PAS ---
variableName = "mDeltaPhiMET2ndJet_PAS"

plot13 = Plot()
## inputs for stacked histograms
plot13.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot13.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot13.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot13.keys            = keys
plot13.xtit            = "#Delta#phi(MET,2^{nd} jet) (rad.)"
plot13.ytit            = "Number of events"
#plot13.xlog            = "yes"
plot13.ylog            = "yes"
plot13.rebin           = 5
#plot13.xmin            = 0
#plot13.xmax            = 3.146
plot13.ymin            = 0.1
plot13.ymax            = 1000000
#plot13.lpos = "bottom-center"
plot13.name            = "mDeltaPhiMET2ndJet_allPreviousCuts"
plot13.addZUncBand     = zUncBand
plot13.makeRatio       = makeRatio
plot13.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- mDeltaPhiMET2ndJet_PAS ---
variableName = "mDeltaPhiMET2ndJet"

plot13_afterOtherDfCuts = Plot()
## inputs for stacked histograms
plot13_afterOtherDfCuts.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot13_afterOtherDfCuts.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot13_afterOtherDfCuts.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot13_afterOtherDfCuts.keys            = keys
plot13_afterOtherDfCuts.xtit            = "#Delta#phi(MET,2^{nd} jet) (rad.)"
plot13_afterOtherDfCuts.ytit            = "Number of events"
#plot13_afterOtherDfCuts.xlog            = "yes"
plot13_afterOtherDfCuts.ylog            = "yes"
plot13_afterOtherDfCuts.rebin           = 5
#plot13_afterOtherDfCuts.xmin            = 0
#plot13_afterOtherDfCuts.xmax            = 3.146
plot13_afterOtherDfCuts.ymin            = 0.1
plot13_afterOtherDfCuts.ymax            = 1000000
#plot13_afterOtherDfCuts.lpos = "bottom-center"
plot13_afterOtherDfCuts.name            = "mDeltaPhiMET2ndJet_afterOtherDeltaPhiCuts"
plot13_afterOtherDfCuts.addZUncBand     = zUncBand
plot13_afterOtherDfCuts.makeRatio       = makeRatio
plot13_afterOtherDfCuts.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)



#--- MTenu_PAS ---
variableName = "MTenu_PAS"

plot14 = Plot()
## inputs for stacked histograms
plot14.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot14.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot14.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot14.keys            = keys
plot14.xtit            = "M_{T}(e\\nu) (GeV)"
plot14.ytit            = "Number of events"
# plot14.ylog            = "yes"
# plot14.rebin           = 1
# plot14.ymin            = 0.00000001
# plot14.ymax            = 20
plot14.ylog            = "no"
plot14.rebin           = "var"
plot14.xmin            = 0
plot14.xmax            = 400
plot14.ymin            = 0
plot14.ymax            = 150
#plot14.lpos = "bottom-center"
plot14.name            = "MTenu_allPreviousCuts_ylin"
plot14.addZUncBand     = zUncBand
plot14.makeRatio       = makeRatio
plot14.xbins           = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,200,300,400,500]
plot14.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

plot14_ylog = Plot()
## inputs for stacked histograms
plot14_ylog.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot14_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot14_ylog.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot14_ylog.keys            = keys
plot14_ylog.xtit            = "M_{T}(e\\nu) (GeV)"
plot14_ylog.ytit            = "Number of events"
plot14_ylog.ylog            = "yes"
plot14_ylog.rebin           = "var" # don't change it (since a rebinning is already applied above on the same histo)
plot14_ylog.xmin            = 0
plot14_ylog.xmax            = 500
plot14_ylog.ymin            = 0.01
plot14_ylog.ymax            = 1000
#plot14_ylog.lpos = "bottom-center"
plot14_ylog.name            = "MTenu_allPreviousCuts"
plot14_ylog.addZUncBand     = zUncBand
plot14_ylog.makeRatio       = makeRatio
plot14_ylog.xbins           = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,200,300,400,500]
plot14_ylog.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- MTenu_PAS DeltaPhi (0,1) ---
variableName = "h1_MTenu_PAS_DeltaPhiMET2ndJet_0_1"

plot14_0_1 = Plot()
## inputs for stacked histograms
plot14_0_1.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot14_0_1.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot14_0_1.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot14_0_1.keys            = keys
plot14_0_1.xtit            = "M_{T}(e\\nu) (GeV) - #Delta#phi(MET,2nd jet)<=1"
plot14_0_1.ytit            = "Number of events"
# plot14_0_1.ylog            = "yes"
# plot14_0_1.rebin           = 1
# plot14_0_1.ymin            = 0.00000001
# plot14_0_1.ymax            = 20
plot14_0_1.ylog            = "no"
plot14_0_1.rebin           = "var"
plot14_0_1.xmin            = 0
plot14_0_1.xmax            = 400
plot14_0_1.ymin            = 0
plot14_0_1.ymax            = 100
#plot14_0_1.lpos = "bottom-center"
plot14_0_1.name            = "MTenu_allPreviousCuts_ylin_Low_DfMETe"
plot14_0_1.addZUncBand     = zUncBand
plot14_0_1.makeRatio       = makeRatio
plot14_0_1.xbins           = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,200,300,400,500]
plot14_0_1.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

#--- MTenu_PAS DeltaPhi 1_2 ---
variableName = "h1_MTenu_PAS_DeltaPhiMET2ndJet_1_2"

plot14_1_2 = Plot()
## inputs for stacked histograms
plot14_1_2.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot14_1_2.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot14_1_2.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot14_1_2.keys            = keys
plot14_1_2.xtit            = "M_{T}(e\\nu) (GeV) - 1<#Delta#phi(MET,2nd jet)<2"
plot14_1_2.ytit            = "Number of events"
# plot14_1_2.ylog            = "yes"
# plot14_1_2.rebin           = 1
# plot14_1_2.ymin            = 0.00000001
# plot14_1_2.ymax            = 20
plot14_1_2.ylog            = "no"
plot14_1_2.rebin           = "var"
plot14_1_2.xmin            = 0
plot14_1_2.xmax            = 400
plot14_1_2.ymin            = 0
plot14_1_2.ymax            = 100
#plot14_1_2.lpos = "bottom-center"
plot14_1_2.name            = "MTenu_allPreviousCuts_ylin_Mid_DfMETe"
plot14_1_2.addZUncBand     = zUncBand
plot14_1_2.makeRatio       = makeRatio
plot14_1_2.xbins           = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,200,300,400,500]
plot14_1_2.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- MTenu_PAS DeltaPhi 2_pi ---
variableName = "h1_MTenu_PAS_DeltaPhiMET2ndJet_2_pi"

plot14_2_pi = Plot()
## inputs for stacked histograms
plot14_2_pi.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot14_2_pi.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot14_2_pi.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot14_2_pi.keys            = keys
plot14_2_pi.xtit            = "M_{T}(e\\nu) (GeV) - #Delta#phi(MET,2nd jet)>=2"
plot14_2_pi.ytit            = "Number of events"
# plot14_2_pi.ylog            = "yes"
# plot14_2_pi.rebin           = 1
# plot14_2_pi.ymin            = 0.00000001
# plot14_2_pi.ymax            = 20
plot14_2_pi.ylog            = "no"
plot14_2_pi.rebin           = "var"
plot14_2_pi.xmin            = 0
plot14_2_pi.xmax            = 400
plot14_2_pi.ymin            = 0
plot14_2_pi.ymax            = 100
#plot14_2_pi.lpos = "bottom-center"
plot14_2_pi.name            = "MTenu_allPreviousCuts_ylin_Large_DfMETe"
plot14_2_pi.addZUncBand     = zUncBand
plot14_2_pi.makeRatio       = makeRatio
plot14_2_pi.xbins           = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,200,300,400,500]
plot14_2_pi.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- MTenu_PAS_plus ---
variableName = "h1_MTenu_PAS_plus"

plot14_plus = Plot()
## inputs for stacked histograms
plot14_plus.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot14_plus.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot14_plus.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot14_plus.keys            = keys
plot14_plus.xtit            = "M_{T}(e+\\nu) (GeV)"
plot14_plus.ytit            = "Number of events"
# plot14_plus.ylog            = "yes"
# plot14_plus.rebin           = 1
# plot14_plus.ymin            = 0.00000001
# plot14_plus.ymax            = 20
plot14_plus.ylog            = "no"
plot14_plus.rebin           = "var"
plot14_plus.xmin            = 0
plot14_plus.xmax            = 400
plot14_plus.ymin            = 0
plot14_plus.ymax            = 70
#plot14_plus.lpos = "bottom-center"
plot14_plus.name            = "MTenu_plus_allPreviousCuts_ylin"
plot14_plus.addZUncBand     = zUncBand
plot14_plus.makeRatio       = makeRatio
plot14_plus.xbins           = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,200,300,400,500]
plot14_plus.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

plot14_plus_ylog = Plot()
## inputs for stacked histograms
plot14_plus_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot14_plus_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot14_plus_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot14_plus_ylog.keys            = keys
plot14_plus_ylog.xtit            = "M_{T}(e+\\nu) (GeV)"
plot14_plus_ylog.ytit            = "Number of events"
plot14_plus_ylog.ylog            = "yes"
plot14_plus_ylog.rebin           = "var" # don't change it (since a rebinning is already applied above on the same histo)
plot14_plus_ylog.xmin            = 0
plot14_plus_ylog.xmax            = 500
plot14_plus_ylog.ymin            = 0.01
plot14_plus_ylog.ymax            = 1000
#plot14_plus_ylog.lpos = "bottom-center"
plot14_plus_ylog.name            = "MTenu_plus_allPreviousCuts"
plot14_plus_ylog.addZUncBand     = zUncBand
plot14_plus_ylog.makeRatio       = makeRatio
plot14_plus_ylog.xbins           = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,200,300,400,500]
plot14_plus_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- MTenu_PAS_minus ---
variableName = "h1_MTenu_PAS_minus"

plot14_minus = Plot()
## inputs for stacked histograms
plot14_minus.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot14_minus.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot14_minus.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot14_minus.keys            = keys
plot14_minus.xtit            = "M_{T}(e-\\nu) (GeV)"
plot14_minus.ytit            = "Number of events"
# plot14_minus.ylog            = "yes"
# plot14_minus.rebin           = 1
# plot14_minus.ymin            = 0.00000001
# plot14_minus.ymax            = 20
plot14_minus.ylog            = "no"
plot14_minus.rebin           = "var"
plot14_minus.xmin            = 0
plot14_minus.xmax            = 400
plot14_minus.ymin            = 0
plot14_minus.ymax            = 70
#plot14_minus.lpos = "bottom-center"
plot14_minus.name            = "MTenu_minus_allPreviousCuts_ylin"
plot14_minus.addZUncBand     = zUncBand
plot14_minus.makeRatio       = makeRatio
plot14_minus.xbins           = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,200,300,400,500]
plot14_minus.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

plot14_minus_ylog = Plot()
## inputs for stacked histograms
plot14_minus_ylog.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot14_minus_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot14_minus_ylog.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot14_minus_ylog.keys            = keys
plot14_minus_ylog.xtit            = "M_{T}(e-\\nu) (GeV)"
plot14_minus_ylog.ytit            = "Number of events"
plot14_minus_ylog.ylog            = "yes"
plot14_minus_ylog.rebin           = "var" # don't change it (since a rebinning is already applied above on the same histo)
plot14_minus_ylog.xmin            = 0
plot14_minus_ylog.xmax            = 500
plot14_minus_ylog.ymin            = 0.01
plot14_minus_ylog.ymax            = 1000
#plot14_minus_ylog.lpos = "bottom-center"
plot14_minus_ylog.name            = "MTenu_minus_allPreviousCuts"
plot14_minus_ylog.addZUncBand     = zUncBand
plot14_minus_ylog.makeRatio       = makeRatio
plot14_minus_ylog.xbins           = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,200,300,400,500]
plot14_minus_ylog.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)



#--- sT_PAS ---
variableName = "sT_PAS"

plot15 = Plot()
## inputs for stacked histograms
plot15.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot15.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot15.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot15.keys            = keys
plot15.xtit            = "St (GeV)"
plot15.ytit            = "Number of events"
plot15.ylog            = "yes"
plot15.rebin           = "var"
plot15.xmin            = 50
plot15.xmax            = 2000
plot15.ymin            = 0.01
plot15.ymax            = 1000
#plot15.lpos = "bottom-center"
plot15.name            = "sT_allPreviousCuts"
plot15.addZUncBand     = zUncBand
plot15.makeRatio       = makeRatio
plot15.xbins           = [50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000,1500,2000]
plot15.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- sTlep_PAS ---
variableName = "sTlep_PAS"

plot15_lep = Plot()
## inputs for stacked histograms
plot15_lep.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot15_lep.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot15_lep.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot15_lep.keys            = keys
plot15_lep.xtit            = "St leptons (GeV)"
plot15_lep.ytit            = "Number of events"
plot15_lep.ylog            = "yes"
plot15_lep.rebin           = "var"
plot15_lep.xmin            = 50
plot15_lep.xmax            = 1000
plot15_lep.ymin            = 0.01
plot15_lep.ymax            = 1000
#plot15_lep.lpos = "bottom-center"
plot15_lep.name            = "sTlep_allPreviousCuts"
plot15_lep.addZUncBand     = zUncBand
plot15_lep.makeRatio       = makeRatio
plot15_lep.xbins           = [50,70,90,110,130,150,170,190,210,230,250,300,400,500,600,700,800,1000]
plot15_lep.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- sTjet_PAS ---
variableName = "sTjet_PAS"

plot15_jet = Plot()
## inputs for stacked histograms
plot15_jet.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot15_jet.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot15_jet.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot15_jet.keys            = keys
plot15_jet.xtit            = "St jets (GeV)"
plot15_jet.ytit            = "Number of events"
plot15_jet.ylog            = "yes"
plot15_jet.rebin           = "var"
plot15_jet.xmin            = 50
plot15_jet.xmax            = 2000
plot15_jet.ymin            = 0.01
plot15_jet.ymax            = 1000
#plot15_jet.lpos = "bottom-center"
plot15_jet.name            = "sTjet_allPreviousCuts"
plot15_jet.addZUncBand     = zUncBand
plot15_jet.makeRatio       = makeRatio
plot15_jet.xbins           = [50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000,1500,2000]
plot15_jet.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)



#--- Mjj_PAS (after preselection) ---
variableName = "Mjj_PAS"

plot16 = Plot()
## inputs for stacked histograms
plot16.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot16.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot16.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot16.keys            = keys
plot16.xtit            = "M(jj) (GeV)"
plot16.ytit            = "Number of events"
plot16.ylog            = "yes"
plot16.rebin           = "var"
plot16.ymin            = 0.01
plot16.ymax            = 1000
plot16.xmin            = 0
plot16.xmax            = 2000
#plot16.lpos = "bottom-center"
plot16.name            = "Mjj_allPreviousCuts"
plot16.addZUncBand     = zUncBand
plot16.makeRatio       = makeRatio
plot16.xbins           = [0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,350,400,500,800,1000,2000]
plot16.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- Mej (after preselection) ---
variableNames = ["Mej_1stPair_PAS","Mej_2ndPair_PAS"]

plot17 = Plot()
## inputs for stacked histograms
plot17.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistosQCD, variableNames, File_QCD, QCDscaleFactor) + generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot17.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot17.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot17.keys            = keys
plot17.xtit            = "M(ej) (GeV)"
plot17.ytit            = "Number of events x 2"
plot17.ylog            = "no"
plot17.rebin           = "var"
plot17.xmin            = 0
plot17.xmax            = 1500
plot17.ymin            = 0.
plot17.ymax            = 150
#plot17.lpos = "bottom-center"
plot17.name            = "Mej_allPreviousCuts_ylin"
plot17.addZUncBand     = zUncBand
plot17.makeRatio       = makeRatio
plot17.xbins           = [0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,325,350,375,400,425,450,475,500,550,600,800,1000,1500]
plot17.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)

plot17_ylog = Plot()
## inputs for stacked histograms
plot17_ylog.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistosQCD, variableNames, File_QCD, QCDscaleFactor) + generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot17_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot17_ylog.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot17_ylog.keys            = keys
plot17_ylog.xtit            = "M(ej) (GeV)"
plot17_ylog.ytit            = "Number of events x 2"
plot17_ylog.ylog            = "yes"
plot17_ylog.rebin           = "var"
plot17_ylog.xmin            = 0
plot17_ylog.xmax            = 2000
plot17_ylog.ymin            = 0.01
plot17_ylog.ymax            = 1000
#plot17_ylog.lpos = "bottom-center"
plot17_ylog.name            = "Mej_allPreviousCuts"
plot17_ylog.addZUncBand     = zUncBand
plot17_ylog.makeRatio       = makeRatio
plot17_ylog.xbins           = [0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,325,350,375,400,425,450,475,500,550,600,800,1500,2000]
plot17_ylog.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)


##--- MTnuj (after preselection) ---
variableNames = ["MTnuj_1stPair_PAS","MTnuj_2ndPair_PAS"]

plot18 = Plot()
## inputs for stacked histograms
plot18.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistosQCD, variableNames, File_QCD, QCDscaleFactor) + generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot18.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot18.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot18.keys            = keys
plot18.xtit            = "M_{T}(\\nuj) (GeV)"
plot18.ytit            = "Number of events x 2"
plot18.ylog            = "yes"
plot18.rebin           = "var"
plot18.xmin            = 0
plot18.xmax            = 1000
plot18.ymin            = 0.01
plot18.ymax            = 1000
#plot18.lpos = "bottom-center"
plot18.name            = "MTnuj_allPreviousCuts"
plot18.addZUncBand     = zUncBand
plot18.makeRatio       = makeRatio
plot18.xbins           = [0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,500,600,800,1000]
plot18.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)


#--- Phi1stEle_PAS ---
variableName = "Phi1stEle_PAS"

plot19 = Plot()
## inputs for stacked histograms
plot19.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot19.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot19.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot19.keys            = keys
plot19.xtit            = "#phi 1st electron (rad.)"
plot19.ytit            = "Number of events"
plot19.rebin           = eta_rebin*2
#plot19.xmin            = -3.15
#plot19.xmax            = 3.15
plot19.ymin            = eta_ymin
plot19.ymax            = eta_ymax
#plot19.lpos            = "top-left"
plot19.name            = "Phi1stEle_allPreviousCuts"
plot19.addZUncBand     = zUncBand
plot19.makeRatio       = makeRatio
plot19.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- Phi1stEle_PAS_EleBarrel  ---
variableName = "h1_Phi1stEle_PAS_EleBarrel"

plot19_EleBarrel = Plot()
## inputs for stacked histograms
plot19_EleBarrel.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot19_EleBarrel.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot19_EleBarrel.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot19_EleBarrel.keys            = keys
plot19_EleBarrel.xtit            = "#phi 1st electron in barrel (rad.)"
plot19_EleBarrel.ytit            = "Number of events"
plot19_EleBarrel.rebin           = eta_rebin*2
#plot19_EleBarrel.xmin            = -3.15
#plot19_EleBarrel.xmax            = 3.15
plot19_EleBarrel.ymin            = eta_ymin
plot19_EleBarrel.ymax            = eta_ymax
#plot19_EleBarrel.lpos            = "top-left"
plot19_EleBarrel.name            = "Phi1stEle_EleBarrel_allPreviousCuts"
plot19_EleBarrel.addZUncBand     = zUncBand
plot19_EleBarrel.makeRatio       = makeRatio
plot19_EleBarrel.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

#--- Phi1stEle_PAS_EleEndcap  ---
variableName = "h1_Phi1stEle_PAS_EleEndcap"

plot19_EleEndcap = Plot()
## inputs for stacked histograms
plot19_EleEndcap.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot19_EleEndcap.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot19_EleEndcap.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot19_EleEndcap.keys            = keys
plot19_EleEndcap.xtit            = "#phi 1st electron in endcap (rad.)"
plot19_EleEndcap.ytit            = "Number of events"
plot19_EleEndcap.rebin           = eta_rebin*2
#plot19_EleEndcap.xmin            = -3.15
#plot19_EleEndcap.xmax            = 3.15
plot19_EleEndcap.ymin            = eta_ymin
plot19_EleEndcap.ymax            = eta_ymax
#plot19_EleEndcap.lpos            = "top-left"
plot19_EleEndcap.name            = "Phi1stEle_EleEndcap_allPreviousCuts"
plot19_EleEndcap.addZUncBand     = zUncBand
plot19_EleEndcap.makeRatio       = makeRatio
plot19_EleEndcap.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- Phi1stJet_PAS ---
variableName = "Phi1stJet_PAS"

plot20 = Plot()
## inputs for stacked histograms
plot20.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot20.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot20.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot20.keys            = keys
plot20.xtit            = "#phi 1st jet (rad.)"
plot20.ytit            = "Number of events"
plot20.rebin           = eta_rebin*2
plot20.ymin            = eta_ymin
plot20.ymax            = eta_ymax
#plot20.lpos            = "top-left"
plot20.name            = "Phi1stJet_allPreviousCuts"
plot20.addZUncBand     = zUncBand
plot20.makeRatio       = makeRatio
plot20.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Phi2ndJet_PAS ---
variableName = "Phi2ndJet_PAS"

plot21 = Plot()
## inputs for stacked histograms
plot21.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot21.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot21.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot21.keys            = keys
plot21.xtit            = "#phi 2nd jet (rad.)"
plot21.ytit            = "Number of events"
plot21.rebin           = eta_rebin*2
plot21.ymin            = eta_ymin
plot21.ymax            = eta_ymax
#plot21.lpos            = "top-left"
plot21.name            = "Phi2ndJet_allPreviousCuts"
plot21.addZUncBand     = zUncBand
plot21.makeRatio       = makeRatio
plot21.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- METPhi_PAS  ---
variableName = "METPhi_PAS"

plot22 = Plot()
## inputs for stacked histograms
plot22.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot22.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot22.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot22.keys            = keys
plot22.xtit            = "pfMET #phi (rad.)"
plot22.ytit            = "Number of events"
plot22.rebin           = eta_rebin*2
plot22.ymin            = eta_ymin
plot22.ymax            = eta_ymax
#plot22.lpos = "bottom-center"
plot22.name            = "METPhi_allPreviousCuts"
plot22.addZUncBand     = zUncBand
plot22.makeRatio       = makeRatio
plot22.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- METPhi_PAS_EleBarrel  ---
variableName = "h1_METPhi_PAS_EleBarrel"

plot22_EleBarrel = Plot()
## inputs for stacked histograms
plot22_EleBarrel.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot22_EleBarrel.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot22_EleBarrel.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot22_EleBarrel.keys            = keys
plot22_EleBarrel.xtit            = "pfMET #phi for electron in barrel (rad.)"
plot22_EleBarrel.ytit            = "Number of events"
plot22_EleBarrel.rebin           = eta_rebin*2
#plot22_EleBarrel.xmin            = -3.15
#plot22_EleBarrel.xmax            = 3.15
plot22_EleBarrel.ymin            = eta_ymin
plot22_EleBarrel.ymax            = eta_ymax
#plot22_EleBarrel.lpos            = "top-left"
plot22_EleBarrel.name            = "METPhi_EleBarrel_allPreviousCuts"
plot22_EleBarrel.addZUncBand     = zUncBand
plot22_EleBarrel.makeRatio       = makeRatio
plot22_EleBarrel.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

#--- METPhi_PAS_EleEndcap  ---
variableName = "h1_METPhi_PAS_EleEndcap"

plot22_EleEndcap = Plot()
## inputs for stacked histograms
plot22_EleEndcap.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot22_EleEndcap.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot22_EleEndcap.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot22_EleEndcap.keys            = keys
plot22_EleEndcap.xtit            = "pfMET #phi for electron in endcap (rad.)"
plot22_EleEndcap.ytit            = "Number of events"
plot22_EleEndcap.rebin           = eta_rebin*2
#plot22_EleEndcap.xmin            = -3.15
#plot22_EleEndcap.xmax            = 3.15
plot22_EleEndcap.ymin            = eta_ymin
plot22_EleEndcap.ymax            = eta_ymax
#plot22_EleEndcap.lpos            = "top-left"
plot22_EleEndcap.name            = "METPhi_EleEndcap_allPreviousCuts"
plot22_EleEndcap.addZUncBand     = zUncBand
plot22_EleEndcap.makeRatio       = makeRatio
plot22_EleEndcap.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- Ptenu_PAS ---
variableName = "Ptenu_PAS"

plot23 = Plot()
## inputs for stacked histograms
plot23.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot23.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot23.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot23.keys            = keys
plot23.xtit            = "Pt(e\\nu) (GeV)"
plot23.ytit            = "Number of events"
plot23.ylog            = "yes"
plot23.rebin           = "var"
plot23.xmin            = 0
plot23.xmax            = 600
plot23.ymin            = 0.01
plot23.ymax            = 1000
#plot23.lpos = "bottom-center"
plot23.name            = "Ptenu_allPreviousCuts"
plot23.addZUncBand     = zUncBand
plot23.makeRatio       = makeRatio
plot23.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,300,400,500,600,700,800,1000]
plot23.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- minDRej ---
variableName = "minDRej"

plot24 = Plot()
## inputs for stacked histograms
plot24.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot24.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot24.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot24.keys            = keys
plot24.xtit            = "min#Delta R(e,jets)"
plot24.ytit            = "Number of events"
plot24.ylog            = "yes"
plot24.rebin           = 2
plot24.xmin            = 0
plot24.xmax            = 7
plot24.ymin            = 0.01
plot24.ymax            = 5000
#plot24.lpos = "bottom-center"
plot24.name            = "minDRej_allPreviousCuts"
plot24.addZUncBand     = zUncBand
plot24.makeRatio       = makeRatio
#plot24.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot24.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- DeltaRjets_PAS ---
variableName = "DeltaRjets_PAS"

plot25 = Plot()
## inputs for stacked histograms
plot25.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot25.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot25.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot25.keys            = keys
plot25.xtit            = "#Delta R(j1,j2)"
plot25.ytit            = "Number of events"
plot25.ylog            = "yes"
plot25.rebin           = 2
plot25.xmin            = 0
plot25.xmax            = 7
plot25.ymin            = 0.01
plot25.ymax            = 5000
#plot25.lpos = "bottom-center"
plot25.name            = "DRjets_allPreviousCuts"
plot25.addZUncBand     = zUncBand
plot25.makeRatio       = makeRatio
#plot25.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot25.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)





# ############################ Plots below to be done after full selection ######################


##--- sT (after full selection) ---
variableName = "sT"

plot30 = Plot()
## inputs for stacked histograms
plot30.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot30.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot30.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot30.keys            = keys
plot30.xtit            = "St (GeV)"
plot30.ytit            = "Number of events"
plot30.ylog            = "yes"
plot30.rebin           = 10
plot30.xmin            = 100
plot30.xmax            = 2000
plot30.ymin            = 0.01
plot30.ymax            = 100
#plot30.lpos = "bottom-center"
plot30.name            = "sT_fullSelection"
plot30.addZUncBand     = zUncBand
plot30.makeRatio       = makeRatio
plot30.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Mej (after full selection) ---
variableNames = ["Mej_1stPair","Mej_2ndPair"]

plot31 = Plot()
## inputs for stacked histograms
plot31.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistosQCD, variableNames, File_QCD, QCDscaleFactor) + generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot31.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot31.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot31.keys            = keys
plot31.xtit            = "M(ej) (GeV)"
plot31.ytit            = "Number of events x 2"
plot31.ylog            = "yes"
plot31.rebin           = 10
plot31.xmin            = 0
plot31.xmax            = 1000
plot31.ymin            = 0.01
plot31.ymax            = 100
#plot31.lpos = "bottom-center"
plot31.name            = "Mej_fullSelection"
plot31.addZUncBand     = zUncBand
plot31.makeRatio       = makeRatio
plot31.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)


##--- MTnuj (after full selection) ---
variableNames = ["MTnuj_1stPair","MTnuj_2ndPair"]

plot32 = Plot()
## inputs for stacked histograms
plot32.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistosQCD, variableNames, File_QCD, QCDscaleFactor) + generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot32.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot32.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot32.keys            = keys
plot32.xtit            = "M_{T}(\\nuj) (GeV)"
plot32.ytit            = "Number of events x 2"
plot32.ylog            = "yes"
plot32.rebin           = 10
plot32.xmin            = 0
plot32.xmax            = 2000
plot32.ymin            = 0.01
plot32.ymax            = 100
#plot32.lpos = "bottom-center"
plot32.name            = "MTnuj_fullSelection"
plot32.addZUncBand     = zUncBand
plot32.makeRatio       = makeRatio
plot32.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)

#--- Ptenu (after full selection) ---
variableName = "Ptenu"

plot33 = Plot()
## inputs for stacked histograms
plot33.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot33.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot33.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot33.keys            = keys
plot33.xtit            = "Pt(e\\nu) (GeV)"
plot33.ytit            = "Number of events"
plot33.ylog            = "yes"
plot33.rebin           = 1
plot33.xmin            = 0
plot33.xmax            = 1000
plot33.ymin            = 0.01
plot33.ymax            = 100
#plot33.lpos = "bottom-center"
plot33.name            = "Ptenu_fullSelection"
plot33.addZUncBand     = zUncBand
plot33.makeRatio       = makeRatio
plot33.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot33.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Charge1stEle (after full selection) ---
variableName = "h1_Charge1stEle_PAS__sT"

plot34 = Plot()
## inputs for stacked histograms
plot34.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot34.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot34.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot34.keys            = keys
plot34.xtit            = "Charge of reconstructed electron"
plot34.ytit            = "Number of events"
plot34.ylog            = "yes"
plot34.rebin           = 1
#plot34.xmin            = -1.0001
#plot34.xmax            = 1.0001
plot34.ymin            = 0.01
plot34.ymax            = 1000
#plot34.lpos = "bottom-center"
plot34.name            = "charge1stEle_fullSelection"
plot34.addZUncBand     = zUncBand
plot34.makeRatio       = makeRatio
plot34.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- Eta1stEle (full selection)  ---
variableName = "h1_Eta1stEle_PAS__sT"

plot35 = Plot()
## inputs for stacked histograms
plot35.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot35.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot35.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot35.keys            = keys
plot35.xtit            = "#eta 1st electron"
plot35.ytit            = "Number of events"
plot35.rebin           = eta_rebin
plot35.xmin            = -3
plot35.xmax            = 3
plot35.ymin            = 0
plot35.ymax            = 5
#plot35.lpos            = "top-left"
plot35.name            = "Eta1stEle_fullSelection"
plot35.addZUncBand     = zUncBand
plot35.makeRatio       = makeRatio
plot35.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)



#--- h1_mDeltaPhiMETEle_fullSel ---
variableName = "h1_mDeltaPhiMETEle_fullSel"

plot110 = Plot()
## inputs for stacked histograms
plot110.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot110.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot110.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot110.keys            = keys
plot110.xtit            = "#Delta#phi(MET,e) (rad.) - fullSel"
plot110.ytit            = "Number of events"
#plot110.xlog            = "yes"
plot110.ylog            = "no"
plot110.rebin           = 20
#plot110.xmin            = 0
#plot110.xmax            = 3.14
plot110.ymin            = 0
plot110.ymax            = 20
#plot110.lpos = "bottom-center"
plot110.name            = "mDeltaPhiMETEle_presel_fullSel"
plot110.addZUncBand     = zUncBand
plot110.makeRatio       = makeRatio
plot110.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_mDeltaPhiMET1stJet_fullSel ---
variableName = "h1_mDeltaPhiMET1stJet_fullSel"

plot111 = Plot()
## inputs for stacked histograms
plot111.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot111.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot111.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot111.keys            = keys
plot111.xtit            = "#Delta#phi(MET,1^{st} jet) (rad.) - fullSel"
plot111.ytit            = "Number of events"
#plot111.xlog            = "yes"
plot111.ylog            = "no"
plot111.rebin           = 20
#plot111.xmin            = 0
#plot111.xmax            = 3.146
plot111.ymin            = 0.
plot111.ymax            = 20
#plot111.lpos = "bottom-center"
plot111.name            = "mDeltaPhiMET1stJet_presel_fullSel"
plot111.addZUncBand     = zUncBand
plot111.makeRatio       = makeRatio
plot111.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_mDeltaPhiMET2ndJet_fullSel ---
variableName = "h1_mDeltaPhiMET2ndJet_fullSel"

plot112 = Plot()
## inputs for stacked histograms
plot112.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot112.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot112.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot112.keys            = keys
plot112.xtit            = "#Delta#phi(MET,2^{nd} jet) (rad.) - fullSel"
plot112.ytit            = "Number of events"
#plot112.xlog            = "yes"
plot112.ylog            = "no"
plot112.rebin           = 20
#plot112.xmin            = 0
#plot112.xmax            = 3.146
plot112.ymin            = 0.
plot112.ymax            = 20
#plot112.lpos = "bottom-center"
plot112.name            = "mDeltaPhiMET2ndJet_presel_fullSel"
plot112.addZUncBand     = zUncBand
plot112.makeRatio       = makeRatio
plot112.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_minDRej_fullSel ---
variableName = "h1_minDRej_fullSel"

plot113 = Plot()
## inputs for stacked histograms
plot113.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot113.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot113.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot113.keys            = keys
plot113.xtit            = "min#Delta R(e,jets) - fullSel"
plot113.ytit            = "Number of events"
plot113.ylog            = "no"
plot113.rebin           = 10
plot113.xmin            = 0
plot113.xmax            = 7
plot113.ymin            = 0.
plot113.ymax            = 25
#plot113.lpos = "bottom-center"
plot113.name            = "minDRej_presel_fullSel"
plot113.addZUncBand     = zUncBand
plot113.makeRatio       = makeRatio
#plot113.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot113.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_DeltaRjets_PAS_fullSel ---
variableName = "h1_DeltaRjets_PAS_fullSel"

plot114 = Plot()
## inputs for stacked histograms
plot114.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot114.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot114.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot114.keys            = keys
plot114.xtit            = "#Delta R(j1,j2) - fullSel"
plot114.ytit            = "Number of events"
plot114.ylog            = "no"
plot114.rebin           = 10
plot114.xmin            = 0
plot114.xmax            = 7
plot114.ymin            = 0.
plot114.ymax            = 20
#plot114.lpos = "bottom-center"
plot114.name            = "DRjets_presel_fullSel"
plot114.addZUncBand     = zUncBand
plot114.makeRatio       = makeRatio
#plot114.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot114.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Njet_fullSel ---
variableName = "h1_Njet_fullSel"

plot115 = Plot()
## inputs for stacked histograms
plot115.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot115.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot115.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot115.keys            = keys
plot115.xtit            = "Number of jets - fullSel"
plot115.ytit            = "Number of events"
plot115.ylog            = "no"
plot115.rebin           = 1
plot115.xmin            = -0.5
plot115.xmax            = 11.5
plot115.ymin            = 0.
plot115.ymax            = 20
#plot115.lpos = "bottom-center"
plot115.name            = "nJet_presel_fullSel"
plot115.addZUncBand     = zUncBand
plot115.makeRatio       = makeRatio
plot115.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)






# ############################ Plots below to be done after pre-selection + high MT ######################

#--- h1_mDeltaPhiMETEle_highMT ---
variableName = "h1_mDeltaPhiMETEle_highMT"

plot90 = Plot()
## inputs for stacked histograms
plot90.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot90.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot90.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot90.keys            = keys
plot90.xtit            = "#Delta#phi(MET,e) (rad.) - highMT"
plot90.ytit            = "Number of events"
#plot90.xlog            = "yes"
plot90.ylog            = "no"
plot90.rebin           = 20
#plot90.xmin            = 0
#plot90.xmax            = 3.14
plot90.ymin            = 0
plot90.ymax            = 20
#plot90.lpos = "bottom-center"
plot90.name            = "mDeltaPhiMETEle_presel_highMT"
plot90.addZUncBand     = zUncBand
plot90.makeRatio       = makeRatio
plot90.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_mDeltaPhiMET1stJet_highMT ---
variableName = "h1_mDeltaPhiMET1stJet_highMT"

plot91 = Plot()
## inputs for stacked histograms
plot91.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot91.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot91.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot91.keys            = keys
plot91.xtit            = "#Delta#phi(MET,1^{st} jet) (rad.) - highMT"
plot91.ytit            = "Number of events"
#plot91.xlog            = "yes"
plot91.ylog            = "no"
plot91.rebin           = 20
#plot91.xmin            = 0
#plot91.xmax            = 3.146
plot91.ymin            = 0.
plot91.ymax            = 20
#plot91.lpos = "bottom-center"
plot91.name            = "mDeltaPhiMET1stJet_presel_highMT"
plot91.addZUncBand     = zUncBand
plot91.makeRatio       = makeRatio
plot91.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_mDeltaPhiMET2ndJet_highMT ---
variableName = "h1_mDeltaPhiMET2ndJet_highMT"

plot92 = Plot()
## inputs for stacked histograms
plot92.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot92.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot92.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot92.keys            = keys
plot92.xtit            = "#Delta#phi(MET,2^{nd} jet) (rad.) - highMT"
plot92.ytit            = "Number of events"
#plot92.xlog            = "yes"
plot92.ylog            = "no"
plot92.rebin           = 20
#plot92.xmin            = 0
#plot92.xmax            = 3.146
plot92.ymin            = 0.
plot92.ymax            = 20
#plot92.lpos = "bottom-center"
plot92.name            = "mDeltaPhiMET2ndJet_presel_highMT"
plot92.addZUncBand     = zUncBand
plot92.makeRatio       = makeRatio
plot92.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_minDRej_highMT ---
variableName = "h1_minDRej_highMT"

plot93 = Plot()
## inputs for stacked histograms
plot93.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot93.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot93.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot93.keys            = keys
plot93.xtit            = "min#Delta R(e,jets) - highMT"
plot93.ytit            = "Number of events"
plot93.ylog            = "no"
plot93.rebin           = 10
plot93.xmin            = 0
plot93.xmax            = 7
plot93.ymin            = 0.
plot93.ymax            = 25
#plot93.lpos = "bottom-center"
plot93.name            = "minDRej_presel_highMT"
plot93.addZUncBand     = zUncBand
plot93.makeRatio       = makeRatio
#plot93.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot93.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_DeltaRjets_PAS_highMT ---
variableName = "h1_DeltaRjets_PAS_highMT"

plot94 = Plot()
## inputs for stacked histograms
plot94.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot94.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot94.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot94.keys            = keys
plot94.xtit            = "#Delta R(j1,j2) - highMT"
plot94.ytit            = "Number of events"
plot94.ylog            = "no"
plot94.rebin           = 10
plot94.xmin            = 0
plot94.xmax            = 7
plot94.ymin            = 0.
plot94.ymax            = 20
#plot94.lpos = "bottom-center"
plot94.name            = "DRjets_presel_highMT"
plot94.addZUncBand     = zUncBand
plot94.makeRatio       = makeRatio
#plot94.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot94.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Njet_highMT ---
variableName = "h1_Njet_highMT"

plot95 = Plot()
## inputs for stacked histograms
plot95.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot95.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot95.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot95.keys            = keys
plot95.xtit            = "Number of jets - highMT"
plot95.ytit            = "Number of events"
plot95.ylog            = "no"
plot95.rebin           = 1
plot95.xmin            = -0.5
plot95.xmax            = 11.5
plot95.ymin            = 0.
plot95.ymax            = 20
#plot95.lpos = "bottom-center"
plot95.name            = "nJet_presel_highMT"
plot95.addZUncBand     = zUncBand
plot95.makeRatio       = makeRatio
plot95.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)




#-----------------------------------------------------------------------------------


# List of plots to be plotted
plots  = [plot0, plot1, plot2, plot_after2, plot3, plot4, plot5, plot6, plot7, plot8, plot9
          , plot6and8, plot7and9
          , plot10, plot11, plot12, plot13, plot13_afterOtherDfCuts, plot14, plot14_ylog
          , plot14_plus, plot14_plus_ylog, plot14_minus, plot14_minus_ylog
          , plot14_0_1, plot14_1_2, plot14_2_pi
          , plot15, plot15_lep, plot15_jet, plot16
          , plot17, plot17_ylog, plot18
          , plot19, plot19_EleBarrel, plot19_EleEndcap, plot20, plot21, plot22, plot22_EleBarrel, plot22_EleEndcap, plot23, plot24, plot25
          , plot30, plot31, plot32, plot33, plot34, plot35
          , plot90, plot91, plot92, plot93, plot94, plot95
          , plot110, plot111, plot112, plot113, plot114, plot115
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

