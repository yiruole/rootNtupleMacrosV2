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


File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_enujjskim_MET45_Jan11Prod/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_fullntuples_MET45_Jan11Prod/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_fullntuples_MET45_Jan11Prod_type1PFMET/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_fullntuples_MET45/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/36.05pb-1_sT_presel_250_Zrescale1.20/analysisClass_enujjSample_plots.root")
##File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/33.2pb-1_sT_presel_250_Zrescale1.20_Wrescale1.06_extraPlotsDec9/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/f/ferencek/scratch0/LQ/CMSSW_3_8_6/test/Leptoquarks/rootNtupleAnalyzerV2/data/output/analysisClass_enujjSample_plots.root")

File_selection    = File_preselection

UseQCDFromData    = 1 # always put an existing file under File_QCD (otherwise the code will crash)

File_QCD          = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_enujjskim_MET45/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")
#File_QCD          = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_fullntuples_MET45/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")
#File_QCD          = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_fullntuples_MET45_type1PFMET/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")
#File_QCD          = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_fullntuples_MET45_more1SC/output_cutTable_enujjSample_QCD_more1SC/analysisClass_enujjSample_QCD_more1SC_plots.root")
#File_QCD          = GetFile("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/35.84pb-1_QCD_UseHLTPrescales_sT_presel_250/analysisClass_enujjSample_QCD_plots.root")
#File_QCD = GetFile("/afs/cern.ch/user/f/ferencek/scratch0/LQ/CMSSW_3_8_6/test/Leptoquarks/rootNtupleAnalyzerV2/data/output/analysisClass_enujjSample_QCD_plots.root")

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


#--- Vtxd01stEle_PAS ---
variableName = "Vtxd01stEle_PAS"

plot_Vtxd0 = Plot()
## inputs for stacked histograms
plot_Vtxd0.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot_Vtxd0.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot_Vtxd0.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot_Vtxd0.keys            = keys
plot_Vtxd0.xtit            = "Vtx d0 1st electron [cm]"
plot_Vtxd0.ytit            = "Number of events"
plot_Vtxd0.rebin           = 5
#plot_Vtxd0.xmin            = -0.01
#plot_Vtxd0.xmax            = 0.01
plot_Vtxd0.ymin            = 0.
plot_Vtxd0.ymax            = 70.
#plot_Vtxd0.lpos            = "top-left"
plot_Vtxd0.name            = "Vtxd01stEle_allPreviousCuts"
plot_Vtxd0.addZUncBand     = zUncBand
plot_Vtxd0.makeRatio       = makeRatio
plot_Vtxd0.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- MissingHits1stEle_PAS ---
variableName = "MissingHits1stEle_PAS"

plot_MissHits = Plot()
## inputs for stacked histograms
plot_MissHits.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot_MissHits.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot_MissHits.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot_MissHits.keys            = keys
plot_MissHits.xtit            = "Missing Hits 1st electron"
plot_MissHits.ytit            = "Number of events"
plot_MissHits.ylog            = "yes"
plot_MissHits.rebin           = 1
#plot_MissHits.xmin            = 0
#plot_MissHits.xmax            = 2.0
plot_MissHits.ymin            = 0.1
plot_MissHits.ymax            = 1000
#plot_MissHits.lpos            = "top-left"
plot_MissHits.name            = "MissingHits1stEle_allPreviousCuts"
plot_MissHits.addZUncBand     = zUncBand
plot_MissHits.makeRatio       = makeRatio
plot_MissHits.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Dist1stEle_PAS ---
variableName = "Dist1stEle_PAS"

plot_Dist = Plot()
## inputs for stacked histograms
plot_Dist.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot_Dist.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot_Dist.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot_Dist.keys            = keys
plot_Dist.xtit            = "|Dist| 1st electron"
plot_Dist.ytit            = "Number of events"
plot_Dist.ylog            = "yes"
plot_Dist.rebin           = "var"
#plot_Dist.xmin            = 0
#plot_Dist.xmax            = 2.0
plot_Dist.ymin            = 0.1
plot_Dist.ymax            = 110.
#plot_Dist.lpos            = "top-left"
plot_Dist.name            = "Dist1stEle_allPreviousCuts"
plot_Dist.addOvfl         = "no"
plot_Dist.addZUncBand     = zUncBand
plot_Dist.makeRatio       = makeRatio
plot_Dist.xbins           = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
plot_Dist.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- DCotTheta1stEle_PAS ---
variableName = "DCotTheta1stEle_PAS"

plot_DCotTheta = Plot()
## inputs for stacked histograms
plot_DCotTheta.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot_DCotTheta.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot_DCotTheta.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot_DCotTheta.keys            = keys
plot_DCotTheta.xtit            = "|#DeltaCot(#theta)| 1st electron"
plot_DCotTheta.ytit            = "Number of events"
plot_DCotTheta.ylog            = "yes"
plot_DCotTheta.rebin           = "var"
#plot_DCotTheta.xmin            = 0
#plot_DCotTheta.xmax            = 2.0
plot_DCotTheta.ymin            = 0.1
plot_DCotTheta.ymax            = 200.
#plot_DCotTheta.lpos            = "top-left"
plot_DCotTheta.name            = "DCotTheta1stEle_allPreviousCuts"
plot_DCotTheta.addOvfl         = "no"
plot_DCotTheta.addZUncBand     = zUncBand
plot_DCotTheta.makeRatio       = makeRatio
plot_DCotTheta.xbins           = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
plot_DCotTheta.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


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


#--- nJet_TCHELBTag ---
variableName = "nJet_TCHELBTag"

plot_TCHELBTag = Plot()
## inputs for stacked histograms
plot_TCHELBTag.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot_TCHELBTag.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot_TCHELBTag.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot_TCHELBTag.keys            = keys
plot_TCHELBTag.xtit            = "Number of b-tagged jets (TCHEL)"
plot_TCHELBTag.ytit            = "Number of events"
plot_TCHELBTag.ylog            = "yes"
plot_TCHELBTag.rebin           = 1
plot_TCHELBTag.xmin            = -0.5
plot_TCHELBTag.xmax            = 11.5
plot_TCHELBTag.ymin            = 0.01
plot_TCHELBTag.ymax            = 1000000
#plot_TCHELBTag.lpos = "bottom-center"
plot_TCHELBTag.name            = "nJetTCHELBtag_allPreviousCuts"
plot_TCHELBTag.addZUncBand     = zUncBand
plot_TCHELBTag.makeRatio       = makeRatio
plot_TCHELBTag.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


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


#--- TCHE1stJet_PAS ---
variableName = "TCHE1stJet_PAS"

plot_TCHE1 = Plot()
## inputs for stacked histograms
plot_TCHE1.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot_TCHE1.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot_TCHE1.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot_TCHE1.keys            = keys
plot_TCHE1.xtit            = "TCHE discriminator 1st jet"
plot_TCHE1.ytit            = "Number of events"
plot_TCHE1.ylog            = "yes"
plot_TCHE1.rebin           = "var"
plot_TCHE1.ymin            = 0.1
plot_TCHE1.ymax            = 150
#plot_TCHE1.lpos            = "top-left"
plot_TCHE1.name            = "TCHE1stJet_allPreviousCuts"
plot_TCHE1.addZUncBand     = zUncBand
plot_TCHE1.makeRatio       = makeRatio
plot_TCHE1.xbins           = [0,0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.0,4,5,7,10,15,20]
plot_TCHE1.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- TCHE2ndJet_PAS ---
variableName = "TCHE2ndJet_PAS"

plot_TCHE2 = Plot()
## inputs for stacked histograms
plot_TCHE2.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot_TCHE2.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot_TCHE2.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot_TCHE2.keys            = keys
plot_TCHE2.xtit            = "TCHE discriminator 2nd jet"
plot_TCHE2.ytit            = "Number of events"
plot_TCHE2.ylog            = "yes"
plot_TCHE2.rebin           = "var"
plot_TCHE2.ymin            = 0.1
plot_TCHE2.ymax            = 150
#plot_TCHE2.lpos            = "top-left"
plot_TCHE2.name            = "TCHE2ndJet_allPreviousCuts"
plot_TCHE2.addZUncBand     = zUncBand
plot_TCHE2.makeRatio       = makeRatio
plot_TCHE2.xbins           = [0,0.4,0.8,1.2,1.6,2.0,2.4,2.8,3.0,4,5,7,10,15,20]
plot_TCHE2.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


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
plot24.xtit            = "min#DeltaR(e,jets)"
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
plot25.xtit            = "#DeltaR(j1,j2)"
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
plot111.ymin            = 0
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
plot112.ymin            = 0
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
plot113.xtit            = "min#DeltaR(e,jets) - fullSel"
plot113.ytit            = "Number of events"
plot113.ylog            = "no"
plot113.rebin           = 10
plot113.xmin            = 0
plot113.xmax            = 7
plot113.ymin            = 0
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
plot114.xtit            = "#DeltaR(j1,j2) - fullSel"
plot114.ytit            = "Number of events"
plot114.ylog            = "no"
plot114.rebin           = 10
plot114.xmin            = 0
plot114.xmax            = 7
plot114.ymin            = 0
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
plot115.ymin            = 0
plot115.ymax            = 20
#plot115.lpos = "bottom-center"
plot115.name            = "nJet_presel_fullSel"
plot115.addZUncBand     = zUncBand
plot115.makeRatio       = makeRatio
plot115.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_NjetTCHELBTag_fullSel ---
variableName = "h1_NjetTCHELBTag_fullSel"

plot116 = Plot()
## inputs for stacked histograms
plot116.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot116.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot116.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot116.keys            = keys
plot116.xtit            = "Number of b-tagged jets (TCHEL) - fullSel"
plot116.ytit            = "Number of events"
plot116.ylog            = "no"
plot116.rebin           = 1
plot116.xmin            = -0.5
plot116.xmax            = 11.5
plot116.ymin            = 0
plot116.ymax            = 20
#plot116.lpos = "bottom-center"
plot116.name            = "nJetTCHELBTag_presel_fullSel"
plot116.addZUncBand     = zUncBand
plot116.makeRatio       = makeRatio
plot116.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Vtxd01stEle_PAS_fullSel ---
variableName = "h1_Vtxd01stEle_PAS_fullSel"

plot117 = Plot()
## inputs for stacked histograms
plot117.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot117.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot117.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot117.keys            = keys
plot117.xtit            = "Vtx d0 1st electron [cm] - fullSel"
plot117.ytit            = "Number of events"
plot117.ylog            = "no"
plot117.rebin           = 10
#plot117.xmin            = -0.01
#plot117.xmax            = 0.01
plot117.ymin            = 0
plot117.ymax            = 8
#plot117.lpos = "bottom-center"
plot117.name            = "Vtxd01stEle_presel_fullSel"
plot117.addZUncBand     = zUncBand
plot117.makeRatio       = makeRatio
plot117.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_MissingHits1stEle_PAS_fullSel ---
variableName = "h1_MissingHits1stEle_PAS_fullSel"

plot118 = Plot()
## inputs for stacked histograms
plot118.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot118.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot118.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot118.keys            = keys
plot118.xtit            = "Missing Hits 1st electron - fullSel"
plot118.ytit            = "Number of events"
plot118.ylog            = "no"
plot118.rebin           = 1
#plot118.xmin            = -0.01
#plot118.xmax            = 0.01
plot118.ymin            = 0
plot118.ymax            = 20
#plot118.lpos = "bottom-center"
plot118.name            = "MissingHits1stEle_presel_fullSel"
plot118.addZUncBand     = zUncBand
plot118.makeRatio       = makeRatio
plot118.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Dist1stEle_PAS_fullSel ---
variableName = "h1_Dist1stEle_PAS_fullSel"

plot119 = Plot()
## inputs for stacked histograms
plot119.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot119.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot119.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot119.keys            = keys
plot119.xtit            = "|Dist| 1st electron - fullSel"
plot119.ytit            = "Number of events"
plot119.ylog            = "no"
plot119.rebin           = 2
plot119.xmin            = 0.
plot119.xmax            = 0.6
plot119.ymin            = 0
plot119.ymax            = 10
#plot119.lpos = "bottom-center"
plot119.name            = "Dist1stEle_presel_fullSel"
plot119.addZUncBand     = zUncBand
plot119.makeRatio       = makeRatio
plot119.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_DCotTheta1stEle_PAS_fullSel ---
variableName = "h1_DCotTheta1stEle_PAS_fullSel"

plot120 = Plot()
## inputs for stacked histograms
plot120.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot120.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot120.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot120.keys            = keys
plot120.xtit            = "|#DeltaCot(#theta)| 1st electron - fullSel"
plot120.ytit            = "Number of events"
plot120.ylog            = "no"
plot120.rebin           = 2
plot120.xmin            = 0.
plot120.xmax            = 0.8
plot120.ymin            = 0
plot120.ymax            = 10
#plot120.lpos = "bottom-center"
plot120.name            = "DCotTheta1stEle_presel_fullSel"
plot120.addZUncBand     = zUncBand
plot120.makeRatio       = makeRatio
plot120.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


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
plot90.ymax            = 25
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
plot91.ymin            = 0
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
#plot92.xmax            = 3.1416
plot92.ymin            = 0
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
plot93.xtit            = "min#DeltaR(e,jets) - highMT"
plot93.ytit            = "Number of events"
plot93.ylog            = "no"
plot93.rebin           = 10
plot93.xmin            = 0
plot93.xmax            = 7
plot93.ymin            = 0
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
plot94.xtit            = "#DeltaR(j1,j2) - highMT"
plot94.ytit            = "Number of events"
plot94.ylog            = "no"
plot94.rebin           = 10
plot94.xmin            = 0
plot94.xmax            = 7
plot94.ymin            = 0
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
plot95.ymin            = 0
plot95.ymax            = 20
#plot95.lpos = "bottom-center"
plot95.name            = "nJet_presel_highMT"
plot95.addZUncBand     = zUncBand
plot95.makeRatio       = makeRatio
plot95.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_NjetTCHELBTag_highMT ---
variableName = "h1_NjetTCHELBTag_highMT"

plot96 = Plot()
## inputs for stacked histograms
plot96.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot96.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot96.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot96.keys            = keys
plot96.xtit            = "Number of b-tagged jets (TCHEL) - highMT"
plot96.ytit            = "Number of events"
plot96.ylog            = "no"
plot96.rebin           = 1
plot96.xmin            = -0.5
plot96.xmax            = 11.5
plot96.ymin            = 0
plot96.ymax            = 20
#plot96.lpos = "bottom-center"
plot96.name            = "nJetTCHELBTag_presel_highMT"
plot96.addZUncBand     = zUncBand
plot96.makeRatio       = makeRatio
plot96.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Vtxd01stEle_PAS_highMT ---
variableName = "h1_Vtxd01stEle_PAS_highMT"

plot97 = Plot()
## inputs for stacked histograms
plot97.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot97.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot97.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot97.keys            = keys
plot97.xtit            = "Vtx d0 1st electron [cm] - highMT"
plot97.ytit            = "Number of events"
plot97.ylog            = "no"
plot97.rebin           = 10
#plot97.xmin            = -0.01
#plot97.xmax            = 0.01
plot97.ymin            = 0
plot97.ymax            = 8
#plot97.lpos = "bottom-center"
plot97.name            = "Vtxd01stEle_presel_highMT"
plot97.addZUncBand     = zUncBand
plot97.makeRatio       = makeRatio
plot97.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_MissingHits1stEle_PAS_highMT ---
variableName = "h1_MissingHits1stEle_PAS_highMT"

plot98 = Plot()
## inputs for stacked histograms
plot98.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot98.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot98.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot98.keys            = keys
plot98.xtit            = "Missing Hits 1st electron - highMT"
plot98.ytit            = "Number of events"
plot98.ylog            = "no"
plot98.rebin           = 1
#plot98.xmin            = -0.01
#plot98.xmax            = 0.01
plot98.ymin            = 0
plot98.ymax            = 20
#plot98.lpos = "bottom-center"
plot98.name            = "MissingHits1stEle_presel_highMT"
plot98.addZUncBand     = zUncBand
plot98.makeRatio       = makeRatio
plot98.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Dist1stEle_PAS_highMT ---
variableName = "h1_Dist1stEle_PAS_highMT"

plot99 = Plot()
## inputs for stacked histograms
plot99.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot99.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot99.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot99.keys            = keys
plot99.xtit            = "|Dist| 1st electron - highMT"
plot99.ytit            = "Number of events"
plot99.ylog            = "no"
plot99.rebin           = 2
plot99.xmin            = 0.
plot99.xmax            = 0.6
plot99.ymin            = 0
plot99.ymax            = 10
#plot99.lpos = "bottom-center"
plot99.name            = "Dist1stEle_presel_highMT"
plot99.addZUncBand     = zUncBand
plot99.makeRatio       = makeRatio
plot99.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_DCotTheta1stEle_PAS_highMT ---
variableName = "h1_DCotTheta1stEle_PAS_highMT"

plot100 = Plot()
## inputs for stacked histograms
plot100.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot100.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot100.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot100.keys            = keys
plot100.xtit            = "|#DeltaCot(#theta)| 1st electron - highMT"
plot100.ytit            = "Number of events"
plot100.ylog            = "no"
plot100.rebin           = 2
plot100.xmin            = 0.
plot100.xmax            = 0.8
plot100.ymin            = 0
plot100.ymax            = 10
#plot100.lpos = "bottom-center"
plot100.name            = "DCotTheta1stEle_presel_highMT"
plot100.addZUncBand     = zUncBand
plot100.makeRatio       = makeRatio
plot100.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


# ############################ Plots below to be done after pre-selection + high Pt1stEle ######################

#--- h1_mDeltaPhiMETEle_highPt1stEle ---
variableName = "h1_mDeltaPhiMETEle_highPt1stEle"

plot190 = Plot()
## inputs for stacked histograms
plot190.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot190.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot190.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot190.keys            = keys
plot190.xtit            = "#Delta#phi(MET,e) (rad.) - highPt1stEle"
plot190.ytit            = "Number of events"
#plot190.xlog            = "yes"
plot190.ylog            = "no"
plot190.rebin           = 20
#plot190.xmin            = 0
#plot190.xmax            = 3.14
plot190.ymin            = 0
plot190.ymax            = 10
#plot190.lpos = "bottom-center"
plot190.name            = "mDeltaPhiMETEle_presel_highPt1stEle"
plot190.addZUncBand     = zUncBand
plot190.makeRatio       = makeRatio
plot190.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_mDeltaPhiMET1stJet_highPt1stEle ---
variableName = "h1_mDeltaPhiMET1stJet_highPt1stEle"

plot191 = Plot()
## inputs for stacked histograms
plot191.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot191.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot191.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot191.keys            = keys
plot191.xtit            = "#Delta#phi(MET,1^{st} jet) (rad.) - highPt1stEle"
plot191.ytit            = "Number of events"
#plot191.xlog            = "yes"
plot191.ylog            = "no"
plot191.rebin           = 20
#plot191.xmin            = 0
#plot191.xmax            = 3.146
plot191.ymin            = 0
plot191.ymax            = 10
#plot191.lpos = "bottom-center"
plot191.name            = "mDeltaPhiMET1stJet_presel_highPt1stEle"
plot191.addZUncBand     = zUncBand
plot191.makeRatio       = makeRatio
plot191.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_mDeltaPhiMET2ndJet_highPt1stEle ---
variableName = "h1_mDeltaPhiMET2ndJet_highPt1stEle"

plot192 = Plot()
## inputs for stacked histograms
plot192.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot192.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot192.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot192.keys            = keys
plot192.xtit            = "#Delta#phi(MET,2^{nd} jet) (rad.) - highPt1stEle"
plot192.ytit            = "Number of events"
#plot192.xlog            = "yes"
plot192.ylog            = "no"
plot192.rebin           = 20
#plot192.xmin            = 0
#plot192.xmax            = 3.146
plot192.ymin            = 0
plot192.ymax            = 10
#plot192.lpos = "bottom-center"
plot192.name            = "mDeltaPhiMET2ndJet_presel_highPt1stEle"
plot192.addZUncBand     = zUncBand
plot192.makeRatio       = makeRatio
plot192.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_minDRej_highPt1stEle ---
variableName = "h1_minDRej_highPt1stEle"

plot193 = Plot()
## inputs for stacked histograms
plot193.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot193.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot193.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot193.keys            = keys
plot193.xtit            = "min#DeltaR(e,jets) - highPt1stEle"
plot193.ytit            = "Number of events"
plot193.ylog            = "no"
plot193.rebin           = 10
plot193.xmin            = 0
plot193.xmax            = 7
plot193.ymin            = 0
plot193.ymax            = 10
#plot193.lpos = "bottom-center"
plot193.name            = "minDRej_presel_highPt1stEle"
plot193.addZUncBand     = zUncBand
plot193.makeRatio       = makeRatio
#plot193.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot193.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_DeltaRjets_PAS_highPt1stEle ---
variableName = "h1_DeltaRjets_PAS_highPt1stEle"

plot194 = Plot()
## inputs for stacked histograms
plot194.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot194.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot194.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot194.keys            = keys
plot194.xtit            = "#DeltaR(j1,j2) - highPt1stEle"
plot194.ytit            = "Number of events"
plot194.ylog            = "no"
plot194.rebin           = 10
plot194.xmin            = 0
plot194.xmax            = 7
plot194.ymin            = 0
plot194.ymax            = 10
#plot194.lpos = "bottom-center"
plot194.name            = "DRjets_presel_highPt1stEle"
plot194.addZUncBand     = zUncBand
plot194.makeRatio       = makeRatio
#plot194.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot194.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Njet_highPt1stEle ---
variableName = "h1_Njet_highPt1stEle"

plot195 = Plot()
## inputs for stacked histograms
plot195.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot195.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot195.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot195.keys            = keys
plot195.xtit            = "Number of jets - highPt1stEle"
plot195.ytit            = "Number of events"
plot195.ylog            = "no"
plot195.rebin           = 1
plot195.xmin            = -0.5
plot195.xmax            = 11.5
plot195.ymin            = 0
plot195.ymax            = 10
#plot195.lpos = "bottom-center"
plot195.name            = "nJet_presel_highPt1stEle"
plot195.addZUncBand     = zUncBand
plot195.makeRatio       = makeRatio
plot195.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_NjetTCHELBTag_highPt1stEle ---
variableName = "h1_NjetTCHELBTag_highPt1stEle"

plot196 = Plot()
## inputs for stacked histograms
plot196.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot196.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot196.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot196.keys            = keys
plot196.xtit            = "Number of b-tagged jets (TCHEL) - highPt1stEle"
plot196.ytit            = "Number of events"
plot196.ylog            = "no"
plot196.rebin           = 1
plot196.xmin            = -0.5
plot196.xmax            = 11.5
plot196.ymin            = 0
plot196.ymax            = 10
#plot196.lpos = "bottom-center"
plot196.name            = "nJetTCHELBTag_presel_highPt1stEle"
plot196.addZUncBand     = zUncBand
plot196.makeRatio       = makeRatio
plot196.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Vtxd01stEle_PAS_highPt1stEle ---
variableName = "h1_Vtxd01stEle_PAS_highPt1stEle"

plot197 = Plot()
## inputs for stacked histograms
plot197.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot197.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot197.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot197.keys            = keys
plot197.xtit            = "Vtx d0 1st electron [cm] - highPt1stEle"
plot197.ytit            = "Number of events"
plot197.ylog            = "no"
plot197.rebin           = 10
#plot197.xmin            = -0.01
#plot197.xmax            = 0.01
plot197.ymin            = 0
plot197.ymax            = 10
#plot197.lpos = "bottom-center"
plot197.name            = "Vtxd01stEle_presel_highPt1stEle"
plot197.addZUncBand     = zUncBand
plot197.makeRatio       = makeRatio
plot197.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_MissingHits1stEle_PAS_highPt1stEle ---
variableName = "h1_MissingHits1stEle_PAS_highPt1stEle"

plot198 = Plot()
## inputs for stacked histograms
plot198.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot198.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot198.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot198.keys            = keys
plot198.xtit            = "Missing Hits 1st electron - highPt1stEle"
plot198.ytit            = "Number of events"
plot198.ylog            = "no"
plot198.rebin           = 1
#plot198.xmin            = -0.01
#plot198.xmax            = 0.01
plot198.ymin            = 0
plot198.ymax            = 10
#plot198.lpos = "bottom-center"
plot198.name            = "MissingHits1stEle_presel_highPt1stEle"
plot198.addZUncBand     = zUncBand
plot198.makeRatio       = makeRatio
plot198.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Dist1stEle_PAS_highPt1stEle ---
variableName = "h1_Dist1stEle_PAS_highPt1stEle"

plot199 = Plot()
## inputs for stacked histograms
plot199.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot199.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot199.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot199.keys            = keys
plot199.xtit            = "|Dist| 1st electron - highPt1stEle"
plot199.ytit            = "Number of events"
plot199.ylog            = "no"
plot199.rebin           = 2
plot199.xmin            = 0.
plot199.xmax            = 0.6
plot199.ymin            = 0
plot199.ymax            = 5
#plot199.lpos = "bottom-center"
plot199.name            = "Dist1stEle_presel_highPt1stEle"
plot199.addZUncBand     = zUncBand
plot199.makeRatio       = makeRatio
plot199.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_DCotTheta1stEle_PAS_highPt1stEle ---
variableName = "h1_DCotTheta1stEle_PAS_highPt1stEle"

plot200 = Plot()
## inputs for stacked histograms
plot200.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot200.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot200.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot200.keys            = keys
plot200.xtit            = "|#DeltaCot(#theta)| 1st electron - highPt1stEle"
plot200.ytit            = "Number of events"
plot200.ylog            = "no"
plot200.rebin           = 2
plot200.xmin            = 0.
plot200.xmax            = 0.8
plot200.ymin            = 0
plot200.ymax            = 5
#plot200.lpos = "bottom-center"
plot200.name            = "DCotTheta1stEle_presel_highPt1stEle"
plot200.addZUncBand     = zUncBand
plot200.makeRatio       = makeRatio
plot200.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


# ############################ Plots below to be done after pre-selection + high Pt1stJet ######################

#--- h1_mDeltaPhiMETEle_highPt1stEle ---
variableName = "h1_mDeltaPhiMETEle_highPt1stJet"

plot290 = Plot()
## inputs for stacked histograms
plot290.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot290.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot290.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot290.keys            = keys
plot290.xtit            = "#Delta#phi(MET,e) (rad.) - highPt1stJet"
plot290.ytit            = "Number of events"
#plot290.xlog            = "yes"
plot290.ylog            = "no"
plot290.rebin           = 20
#plot290.xmin            = 0
#plot290.xmax            = 3.14
plot290.ymin            = 0
plot290.ymax            = 10
#plot290.lpos = "bottom-center"
plot290.name            = "mDeltaPhiMETEle_presel_highPt1stJet"
plot290.addZUncBand     = zUncBand
plot290.makeRatio       = makeRatio
plot290.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_mDeltaPhiMET1stJet_highPt1stJet ---
variableName = "h1_mDeltaPhiMET1stJet_highPt1stJet"

plot291 = Plot()
## inputs for stacked histograms
plot291.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot291.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot291.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot291.keys            = keys
plot291.xtit            = "#Delta#phi(MET,1^{st} jet) (rad.) - highPt1stJet"
plot291.ytit            = "Number of events"
#plot291.xlog            = "yes"
plot291.ylog            = "no"
plot291.rebin           = 20
#plot291.xmin            = 0
#plot291.xmax            = 3.146
plot291.ymin            = 0
plot291.ymax            = 15
#plot291.lpos = "bottom-center"
plot291.name            = "mDeltaPhiMET1stJet_presel_highPt1stJet"
plot291.addZUncBand     = zUncBand
plot291.makeRatio       = makeRatio
plot291.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_mDeltaPhiMET2ndJet_highPt1stJet ---
variableName = "h1_mDeltaPhiMET2ndJet_highPt1stJet"

plot292 = Plot()
## inputs for stacked histograms
plot292.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot292.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot292.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot292.keys            = keys
plot292.xtit            = "#Delta#phi(MET,2^{nd} jet) (rad.) - highPt1stJet"
plot292.ytit            = "Number of events"
#plot292.xlog            = "yes"
plot292.ylog            = "no"
plot292.rebin           = 20
#plot292.xmin            = 0
#plot292.xmax            = 3.146
plot292.ymin            = 0
plot292.ymax            = 10
#plot292.lpos = "bottom-center"
plot292.name            = "mDeltaPhiMET2ndJet_presel_highPt1stJet"
plot292.addZUncBand     = zUncBand
plot292.makeRatio       = makeRatio
plot292.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_minDRej_highPt1stJet ---
variableName = "h1_minDRej_highPt1stJet"

plot293 = Plot()
## inputs for stacked histograms
plot293.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot293.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot293.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot293.keys            = keys
plot293.xtit            = "min#DeltaR(e,jets) - highPt1stJet"
plot293.ytit            = "Number of events"
plot293.ylog            = "no"
plot293.rebin           = 10
plot293.xmin            = 0
plot293.xmax            = 7
plot293.ymin            = 0
plot293.ymax            = 10
#plot293.lpos = "bottom-center"
plot293.name            = "minDRej_presel_highPt1stJet"
plot293.addZUncBand     = zUncBand
plot293.makeRatio       = makeRatio
#plot293.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot293.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_DeltaRjets_PAS_highPt1stJet ---
variableName = "h1_DeltaRjets_PAS_highPt1stJet"

plot294 = Plot()
## inputs for stacked histograms
plot294.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot294.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot294.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot294.keys            = keys
plot294.xtit            = "#DeltaR(j1,j2) - highPt1stJet"
plot294.ytit            = "Number of events"
plot294.ylog            = "no"
plot294.rebin           = 10
plot294.xmin            = 0
plot294.xmax            = 7
plot294.ymin            = 0
plot294.ymax            = 10
#plot294.lpos = "bottom-center"
plot294.name            = "DRjets_presel_highPt1stJet"
plot294.addZUncBand     = zUncBand
plot294.makeRatio       = makeRatio
#plot294.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot294.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Njet_highPt1stJet ---
variableName = "h1_Njet_highPt1stJet"

plot295 = Plot()
## inputs for stacked histograms
plot295.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot295.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot295.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot295.keys            = keys
plot295.xtit            = "Number of jets - highPt1stJet"
plot295.ytit            = "Number of events"
plot295.ylog            = "no"
plot295.rebin           = 1
plot295.xmin            = -0.5
plot295.xmax            = 11.5
plot295.ymin            = 0
plot295.ymax            = 10
#plot295.lpos = "bottom-center"
plot295.name            = "nJet_presel_highPt1stJet"
plot295.addZUncBand     = zUncBand
plot295.makeRatio       = makeRatio
plot295.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_NjetTCHELBTag_highPt1stJet ---
variableName = "h1_NjetTCHELBTag_highPt1stJet"

plot296 = Plot()
## inputs for stacked histograms
plot296.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot296.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot296.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot296.keys            = keys
plot296.xtit            = "Number of b-tagged jets (TCHEL) - highPt1stJet"
plot296.ytit            = "Number of events"
plot296.ylog            = "no"
plot296.rebin           = 1
plot296.xmin            = -0.5
plot296.xmax            = 11.5
plot296.ymin            = 0
plot296.ymax            = 10
#plot296.lpos = "bottom-center"
plot296.name            = "nJetTCHELBTag_presel_highPt1stJet"
plot296.addZUncBand     = zUncBand
plot296.makeRatio       = makeRatio
plot296.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Vtxd01stEle_PAS_highPt1stJet ---
variableName = "h1_Vtxd01stEle_PAS_highPt1stJet"

plot297 = Plot()
## inputs for stacked histograms
plot297.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot297.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot297.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot297.keys            = keys
plot297.xtit            = "Vtx d0 1st electron [cm] - highPt1stJet"
plot297.ytit            = "Number of events"
plot297.ylog            = "no"
plot297.rebin           = 10
#plot297.xmin            = -0.01
#plot297.xmax            = 0.01
plot297.ymin            = 0
plot297.ymax            = 8
#plot297.lpos = "bottom-center"
plot297.name            = "Vtxd01stEle_presel_highPt1stJet"
plot297.addZUncBand     = zUncBand
plot297.makeRatio       = makeRatio
plot297.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_MissingHits1stEle_PAS_highPt1stJet ---
variableName = "h1_MissingHits1stEle_PAS_highPt1stJet"

plot298 = Plot()
## inputs for stacked histograms
plot298.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot298.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot298.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot298.keys            = keys
plot298.xtit            = "Missing Hits 1st electron - highPt1stJet"
plot298.ytit            = "Number of events"
plot298.ylog            = "no"
plot298.rebin           = 1
#plot298.xmin            = -0.01
#plot298.xmax            = 0.01
plot298.ymin            = 0
plot298.ymax            = 20
#plot298.lpos = "bottom-center"
plot298.name            = "MissingHits1stEle_presel_highPt1stJet"
plot298.addZUncBand     = zUncBand
plot298.makeRatio       = makeRatio
plot298.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Dist1stEle_PAS_highPt1stJet ---
variableName = "h1_Dist1stEle_PAS_highPt1stJet"

plot299 = Plot()
## inputs for stacked histograms
plot299.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot299.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot299.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot299.keys            = keys
plot299.xtit            = "|Dist| 1st electron - highPt1stJet"
plot299.ytit            = "Number of events"
plot299.ylog            = "no"
plot299.rebin           = 2
plot299.xmin            = 0.
plot299.xmax            = 0.6
plot299.ymin            = 0
plot299.ymax            = 15
#plot299.lpos = "bottom-center"
plot299.name            = "Dist1stEle_presel_highPt1stJet"
plot299.addZUncBand     = zUncBand
plot299.makeRatio       = makeRatio
plot299.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_DCotTheta1stEle_PAS_highPt1stJet ---
variableName = "h1_DCotTheta1stEle_PAS_highPt1stJet"

plot300 = Plot()
## inputs for stacked histograms
plot300.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot300.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot300.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot300.keys            = keys
plot300.xtit            = "|#DeltaCot(#theta)| 1st electron - highPt1stJet"
plot300.ytit            = "Number of events"
plot300.ylog            = "no"
plot300.rebin           = 2
plot300.xmin            = 0.
plot300.xmax            = 0.8
plot300.ymin            = 0
plot300.ymax            = 10
#plot300.lpos = "bottom-center"
plot300.name            = "DCotTheta1stEle_presel_highPt1stJet"
plot300.addZUncBand     = zUncBand
plot300.makeRatio       = makeRatio
plot300.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


# ############################ Plots below to be done after pre-selection + high Pt2ndJet ######################

#--- h1_mDeltaPhiMETEle_highPt1stEle ---
variableName = "h1_mDeltaPhiMETEle_highPt2ndJet"

plot390 = Plot()
## inputs for stacked histograms
plot390.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot390.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot390.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot390.keys            = keys
plot390.xtit            = "#Delta#phi(MET,e) (rad.) - highPt2ndJet"
plot390.ytit            = "Number of events"
#plot390.xlog            = "yes"
plot390.ylog            = "no"
plot390.rebin           = 20
#plot390.xmin            = 0
#plot390.xmax            = 3.14
plot390.ymin            = 0
plot390.ymax            = 10
#plot390.lpos = "bottom-center"
plot390.name            = "mDeltaPhiMETEle_presel_highPt2ndJet"
plot390.addZUncBand     = zUncBand
plot390.makeRatio       = makeRatio
plot390.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_mDeltaPhiMET1stJet_highPt2ndJet ---
variableName = "h1_mDeltaPhiMET1stJet_highPt2ndJet"

plot391 = Plot()
## inputs for stacked histograms
plot391.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot391.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot391.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot391.keys            = keys
plot391.xtit            = "#Delta#phi(MET,1^{st} jet) (rad.) - highPt2ndJet"
plot391.ytit            = "Number of events"
#plot391.xlog            = "yes"
plot391.ylog            = "no"
plot391.rebin           = 20
#plot391.xmin            = 0
#plot391.xmax            = 3.146
plot391.ymin            = 0
plot391.ymax            = 15
#plot391.lpos = "bottom-center"
plot391.name            = "mDeltaPhiMET1stJet_presel_highPt2ndJet"
plot391.addZUncBand     = zUncBand
plot391.makeRatio       = makeRatio
plot391.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_mDeltaPhiMET2ndJet_highPt2ndJet ---
variableName = "h1_mDeltaPhiMET2ndJet_highPt2ndJet"

plot392 = Plot()
## inputs for stacked histograms
plot392.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot392.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot392.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot392.keys            = keys
plot392.xtit            = "#Delta#phi(MET,2^{nd} jet) (rad.) - highPt2ndJet"
plot392.ytit            = "Number of events"
#plot392.xlog            = "yes"
plot392.ylog            = "no"
plot392.rebin           = 20
#plot392.xmin            = 0
#plot392.xmax            = 3.146
plot392.ymin            = 0
plot392.ymax            = 10
#plot392.lpos = "bottom-center"
plot392.name            = "mDeltaPhiMET2ndJet_presel_highPt2ndJet"
plot392.addZUncBand     = zUncBand
plot392.makeRatio       = makeRatio
plot392.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_minDRej_highPt2ndJet ---
variableName = "h1_minDRej_highPt2ndJet"

plot393 = Plot()
## inputs for stacked histograms
plot393.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot393.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot393.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot393.keys            = keys
plot393.xtit            = "min#DeltaR(e,jets) - highPt2ndJet"
plot393.ytit            = "Number of events"
plot393.ylog            = "no"
plot393.rebin           = 10
plot393.xmin            = 0
plot393.xmax            = 7
plot393.ymin            = 0
plot393.ymax            = 10
#plot393.lpos = "bottom-center"
plot393.name            = "minDRej_presel_highPt2ndJet"
plot393.addZUncBand     = zUncBand
plot393.makeRatio       = makeRatio
#plot393.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot393.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_DeltaRjets_PAS_highPt2ndJet ---
variableName = "h1_DeltaRjets_PAS_highPt2ndJet"

plot394 = Plot()
## inputs for stacked histograms
plot394.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot394.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot394.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot394.keys            = keys
plot394.xtit            = "#DeltaR(j1,j2) - highPt2ndJet"
plot394.ytit            = "Number of events"
plot394.ylog            = "no"
plot394.rebin           = 10
plot394.xmin            = 0
plot394.xmax            = 7
plot394.ymin            = 0
plot394.ymax            = 10
#plot394.lpos = "bottom-center"
plot394.name            = "DRjets_presel_highPt2ndJet"
plot394.addZUncBand     = zUncBand
plot394.makeRatio       = makeRatio
#plot394.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot394.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Njet_highPt2ndJet ---
variableName = "h1_Njet_highPt2ndJet"

plot395 = Plot()
## inputs for stacked histograms
plot395.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot395.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot395.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot395.keys            = keys
plot395.xtit            = "Number of jets - highPt2ndJet"
plot395.ytit            = "Number of events"
plot395.ylog            = "no"
plot395.rebin           = 1
plot395.xmin            = -0.5
plot395.xmax            = 11.5
plot395.ymin            = 0
plot395.ymax            = 10
#plot395.lpos = "bottom-center"
plot395.name            = "nJet_presel_highPt2ndJet"
plot395.addZUncBand     = zUncBand
plot395.makeRatio       = makeRatio
plot395.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_NjetTCHELBTag_highPt2ndJet ---
variableName = "h1_NjetTCHELBTag_highPt2ndJet"

plot396 = Plot()
## inputs for stacked histograms
plot396.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot396.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot396.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot396.keys            = keys
plot396.xtit            = "Number of b-tagged jets (TCHEL) - highPt2ndJet"
plot396.ytit            = "Number of events"
plot396.ylog            = "no"
plot396.rebin           = 1
plot396.xmin            = -0.5
plot396.xmax            = 11.5
plot396.ymin            = 0
plot396.ymax            = 10
#plot396.lpos = "bottom-center"
plot396.name            = "nJetTCHELBTag_presel_highPt2ndJet"
plot396.addZUncBand     = zUncBand
plot396.makeRatio       = makeRatio
plot396.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Vtxd01stEle_PAS_highPt2ndJet ---
variableName = "h1_Vtxd01stEle_PAS_highPt2ndJet"

plot397 = Plot()
## inputs for stacked histograms
plot397.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot397.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot397.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot397.keys            = keys
plot397.xtit            = "Vtx d0 1st electron [cm] - highPt2ndJet"
plot397.ytit            = "Number of events"
plot397.ylog            = "no"
plot397.rebin           = 10
#plot397.xmin            = -0.01
#plot397.xmax            = 0.01
plot397.ymin            = 0
plot397.ymax            = 8
#plot397.lpos = "bottom-center"
plot397.name            = "Vtxd01stEle_presel_highPt2ndJet"
plot397.addZUncBand     = zUncBand
plot397.makeRatio       = makeRatio
plot397.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_MissingHits1stEle_PAS_highPt2ndJet ---
variableName = "h1_MissingHits1stEle_PAS_highPt2ndJet"

plot398 = Plot()
## inputs for stacked histograms
plot398.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot398.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot398.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot398.keys            = keys
plot398.xtit            = "Missing Hits 1st electron - highPt2ndJet"
plot398.ytit            = "Number of events"
plot398.ylog            = "no"
plot398.rebin           = 1
#plot398.xmin            = -0.01
#plot398.xmax            = 0.01
plot398.ymin            = 0
plot398.ymax            = 20
#plot398.lpos = "bottom-center"
plot398.name            = "MissingHits1stEle_presel_highPt2ndJet"
plot398.addZUncBand     = zUncBand
plot398.makeRatio       = makeRatio
plot398.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Dist1stEle_PAS_highPt2ndJet ---
variableName = "h1_Dist1stEle_PAS_highPt2ndJet"

plot399 = Plot()
## inputs for stacked histograms
plot399.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot399.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot399.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot399.keys            = keys
plot399.xtit            = "|Dist| 1st electron - highPt2ndJet"
plot399.ytit            = "Number of events"
plot399.ylog            = "no"
plot399.rebin           = 2
plot399.xmin            = 0.
plot399.xmax            = 0.6
plot399.ymin            = 0
plot399.ymax            = 5
#plot399.lpos = "bottom-center"
plot399.name            = "Dist1stEle_presel_highPt2ndJet"
plot399.addZUncBand     = zUncBand
plot399.makeRatio       = makeRatio
plot399.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_DCotTheta1stEle_PAS_highPt2ndJet ---
variableName = "h1_DCotTheta1stEle_PAS_highPt2ndJet"

plot400 = Plot()
## inputs for stacked histograms
plot400.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot400.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot400.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot400.keys            = keys
plot400.xtit            = "|#DeltaCot(#theta)| 1st electron - highPt2ndJet"
plot400.ytit            = "Number of events"
plot400.ylog            = "no"
plot400.rebin           = 2
plot400.xmin            = 0.
plot400.xmax            = 0.8
plot400.ymin            = 0
plot400.ymax            = 5
#plot400.lpos = "bottom-center"
plot400.name            = "DCotTheta1stEle_presel_highPt2ndJet"
plot400.addZUncBand     = zUncBand
plot400.makeRatio       = makeRatio
plot400.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


# ############################ Plots below to be done after pre-selection + high Mej ######################

#--- h1_mDeltaPhiMETEle_highPt1stEle ---
variableName = "h1_mDeltaPhiMETEle_highMej"

plot490 = Plot()
## inputs for stacked histograms
plot490.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot490.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot490.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot490.keys            = keys
plot490.xtit            = "#Delta#phi(MET,e) (rad.) - highMej"
plot490.ytit            = "Number of events"
#plot490.xlog            = "yes"
plot490.ylog            = "no"
plot490.rebin           = 20
#plot490.xmin            = 0
#plot490.xmax            = 3.14
plot490.ymin            = 0
plot490.ymax            = 50
#plot490.lpos = "bottom-center"
plot490.name            = "mDeltaPhiMETEle_presel_highMej"
plot490.addZUncBand     = zUncBand
plot490.makeRatio       = makeRatio
plot490.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_mDeltaPhiMET1stJet_highMej ---
variableName = "h1_mDeltaPhiMET1stJet_highMej"

plot491 = Plot()
## inputs for stacked histograms
plot491.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot491.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot491.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot491.keys            = keys
plot491.xtit            = "#Delta#phi(MET,1^{st} jet) (rad.) - highMej"
plot491.ytit            = "Number of events"
#plot491.xlog            = "yes"
plot491.ylog            = "no"
plot491.rebin           = 20
#plot491.xmin            = 0
#plot491.xmax            = 3.146
plot491.ymin            = 0
plot491.ymax            = 50
#plot491.lpos = "bottom-center"
plot491.name            = "mDeltaPhiMET1stJet_presel_highMej"
plot491.addZUncBand     = zUncBand
plot491.makeRatio       = makeRatio
plot491.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_mDeltaPhiMET2ndJet_highMej ---
variableName = "h1_mDeltaPhiMET2ndJet_highMej"

plot492 = Plot()
## inputs for stacked histograms
plot492.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot492.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot492.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot492.keys            = keys
plot492.xtit            = "#Delta#phi(MET,2^{nd} jet) (rad.) - highMej"
plot492.ytit            = "Number of events"
#plot492.xlog            = "yes"
plot492.ylog            = "no"
plot492.rebin           = 20
#plot492.xmin            = 0
#plot492.xmax            = 3.146
plot492.ymin            = 0
plot492.ymax            = 50
#plot492.lpos = "bottom-center"
plot492.name            = "mDeltaPhiMET2ndJet_presel_highMej"
plot492.addZUncBand     = zUncBand
plot492.makeRatio       = makeRatio
plot492.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_minDRej_highMej ---
variableName = "h1_minDRej_highMej"

plot493 = Plot()
## inputs for stacked histograms
plot493.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot493.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot493.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot493.keys            = keys
plot493.xtit            = "min#DeltaR(e,jets) - highMej"
plot493.ytit            = "Number of events"
plot493.ylog            = "no"
plot493.rebin           = 10
plot493.xmin            = 0
plot493.xmax            = 7
plot493.ymin            = 0
plot493.ymax            = 50
#plot493.lpos = "bottom-center"
plot493.name            = "minDRej_presel_highMej"
plot493.addZUncBand     = zUncBand
plot493.makeRatio       = makeRatio
#plot493.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot493.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_DeltaRjets_PAS_highMej ---
variableName = "h1_DeltaRjets_PAS_highMej"

plot494 = Plot()
## inputs for stacked histograms
plot494.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot494.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot494.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot494.keys            = keys
plot494.xtit            = "#DeltaR(j1,j2) - highMej"
plot494.ytit            = "Number of events"
plot494.ylog            = "no"
plot494.rebin           = 10
plot494.xmin            = 0
plot494.xmax            = 7
plot494.ymin            = 0
plot494.ymax            = 50
#plot494.lpos = "bottom-center"
plot494.name            = "DRjets_presel_highMej"
plot494.addZUncBand     = zUncBand
plot494.makeRatio       = makeRatio
#plot494.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
plot494.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Njet_highMej ---
variableName = "h1_Njet_highMej"

plot495 = Plot()
## inputs for stacked histograms
plot495.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot495.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot495.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot495.keys            = keys
plot495.xtit            = "Number of jets - highMej"
plot495.ytit            = "Number of events"
plot495.ylog            = "no"
plot495.rebin           = 1
plot495.xmin            = -0.5
plot495.xmax            = 11.5
plot495.ymin            = 0
plot495.ymax            = 50
#plot495.lpos = "bottom-center"
plot495.name            = "nJet_presel_highMej"
plot495.addZUncBand     = zUncBand
plot495.makeRatio       = makeRatio
plot495.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_NjetTCHELBTag_highMej ---
variableName = "h1_NjetTCHELBTag_highMej"

plot496 = Plot()
## inputs for stacked histograms
plot496.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
plot496.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot496.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot496.keys            = keys
plot496.xtit            = "Number of b-tagged jets (TCHEL) - highMej"
plot496.ytit            = "Number of events"
plot496.ylog            = "no"
plot496.rebin           = 1
plot496.xmin            = -0.5
plot496.xmax            = 11.5
plot496.ymin            = 0
plot496.ymax            = 50
#plot496.lpos = "bottom-center"
plot496.name            = "nJetTCHELBTag_presel_highMej"
plot496.addZUncBand     = zUncBand
plot496.makeRatio       = makeRatio
plot496.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Vtxd01stEle_PAS_highMej ---
variableName = "h1_Vtxd01stEle_PAS_highMej"

plot497 = Plot()
## inputs for stacked histograms
plot497.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot497.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot497.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot497.keys            = keys
plot497.xtit            = "Vtx d0 1st electron [cm] - highMej"
plot497.ytit            = "Number of events"
plot497.ylog            = "no"
plot497.rebin           = 10
#plot497.xmin            = -0.01
#plot497.xmax            = 0.01
plot497.ymin            = 0
plot497.ymax            = 15
#plot497.lpos = "bottom-center"
plot497.name            = "Vtxd01stEle_presel_highMej"
plot497.addZUncBand     = zUncBand
plot497.makeRatio       = makeRatio
plot497.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_MissingHits1stEle_PAS_highMej ---
variableName = "h1_MissingHits1stEle_PAS_highMej"

plot498 = Plot()
## inputs for stacked histograms
plot498.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot498.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot498.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot498.keys            = keys
plot498.xtit            = "Missing Hits 1st electron - highMej"
plot498.ytit            = "Number of events"
plot498.ylog            = "no"
plot498.rebin           = 1
#plot498.xmin            = -0.01
#plot498.xmax            = 0.01
plot498.ymin            = 0
plot498.ymax            = 80
#plot498.lpos = "bottom-center"
plot498.name            = "MissingHits1stEle_presel_highMej"
plot498.addZUncBand     = zUncBand
plot498.makeRatio       = makeRatio
plot498.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_Dist1stEle_PAS_highMej ---
variableName = "h1_Dist1stEle_PAS_highMej"

plot499 = Plot()
## inputs for stacked histograms
plot499.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot499.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot499.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot499.keys            = keys
plot499.xtit            = "|Dist| 1st electron - highMej"
plot499.ytit            = "Number of events"
plot499.ylog            = "no"
plot499.rebin           = 2
plot499.xmin            = 0.
plot499.xmax            = 0.6
plot499.ymin            = 0
plot499.ymax            = 25
#plot499.lpos = "bottom-center"
plot499.name            = "Dist1stEle_presel_highMej"
plot499.addZUncBand     = zUncBand
plot499.makeRatio       = makeRatio
plot499.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#--- h1_DCotTheta1stEle_PAS_highMej ---
variableName = "h1_DCotTheta1stEle_PAS_highMej"

plot500 = Plot()
## inputs for stacked histograms
plot500.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
keysStacknoQCD = copy.deepcopy(keysStack)
keysStacknoQCD.remove("QCD")
plot500.keysStack       = keysStacknoQCD
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot500.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
plot500.keys            = keys
plot500.xtit            = "|#DeltaCot(#theta)| 1st electron - highMej"
plot500.ytit            = "Number of events"
plot500.ylog            = "no"
plot500.rebin           = 2
plot500.xmin            = 0.
plot500.xmax            = 0.8
plot500.ymin            = 0
plot500.ymax            = 25
#plot500.lpos = "bottom-center"
plot500.name            = "DCotTheta1stEle_presel_highMej"
plot500.addZUncBand     = zUncBand
plot500.makeRatio       = makeRatio
plot500.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#-----------------------------------------------------------------------------------


# List of plots to be plotted
plots  = [plot0, plot1, plot2, plot_after2, plot_Vtxd0, plot_MissHits, plot_Dist, plot_DCotTheta, plot3, plot4, plot5, plot_TCHELBTag, plot6, plot7, plot8, plot9, plot_TCHE1, plot_TCHE2
          , plot6and8, plot7and9
          , plot10, plot11, plot12, plot13, plot13_afterOtherDfCuts, plot14, plot14_ylog
          , plot14_plus, plot14_plus_ylog, plot14_minus, plot14_minus_ylog
          , plot14_0_1, plot14_1_2, plot14_2_pi
          , plot15, plot15_lep, plot15_jet, plot16
          , plot17, plot17_ylog, plot18
          , plot19, plot19_EleBarrel, plot19_EleEndcap, plot20, plot21, plot22, plot22_EleBarrel, plot22_EleEndcap, plot23, plot24, plot25
          , plot30, plot31, plot32, plot33, plot34, plot35
          , plot90, plot91, plot92, plot93, plot94, plot95, plot96, plot97, plot98, plot99, plot100
          , plot190, plot191, plot192, plot193, plot194, plot195, plot196, plot197, plot198, plot199, plot200
          , plot290, plot291, plot292, plot293, plot294, plot295, plot296, plot297, plot298, plot299, plot300
          , plot390, plot391, plot392, plot393, plot394, plot395, plot396, plot397, plot398, plot399, plot400
          , plot490, plot491, plot492, plot493, plot494, plot495, plot496, plot497, plot498, plot499, plot500
          , plot110, plot111, plot112, plot113, plot114, plot115, plot116, plot117, plot118, plot119, plot120

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

