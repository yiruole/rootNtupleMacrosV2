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


#File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_enujjskim_MET45_Jan11Prod/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_fullntuples_MET45_Jan11Prod/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_fullntuples_MET45_Jan11Prod_type1PFMET/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_fullntuples_MET45/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/36.05pb-1_sT_presel_250_Zrescale1.20/analysisClass_enujjSample_plots.root")
##File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/33.2pb-1_sT_presel_250_Zrescale1.20_Wrescale1.06_extraPlotsDec9/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
File_preselection = GetFile("/afs/cern.ch/user/f/ferencek/scratch0/LQ/CMSSW_3_8_6/test/Leptoquarks/rootNtupleAnalyzerV2/data/output/analysisClass_enujjSample_plots.root")

File_selection    = File_preselection

UseQCDFromData    = 1 # always put an existing file under File_QCD (otherwise the code will crash)

#File_QCD          = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_enujjskim_MET45/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")
#File_QCD          = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_fullntuples_MET45/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")
#File_QCD          = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_fullntuples_MET45_type1PFMET/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")
#File_QCD          = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_fullntuples_MET45_more1SC/output_cutTable_enujjSample_QCD_more1SC/analysisClass_enujjSample_QCD_more1SC_plots.root")
#File_QCD          = GetFile("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/35.84pb-1_QCD_UseHLTPrescales_sT_presel_250/analysisClass_enujjSample_QCD_plots.root")
File_QCD = GetFile("/afs/cern.ch/user/f/ferencek/scratch0/LQ/CMSSW_3_8_6/test/Leptoquarks/rootNtupleAnalyzerV2/data/output/analysisClass_enujjSample_QCD_plots.root")

QCDscaleFactor    = 1 # no need to rescale anymore since we are using the HLT prescales (36/35.84 can be ignored)

#### Common values for plots:
#otherBkgsKey="QCD, single top, VV+jets, W/W*+jets"
otherBkgsKey="Other Bkgs"
zUncBand="no"
makeRatio=1
doExtraPlots = False

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
plot7and9.ymax            = eta_ymax*1.5
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



# ############# Extra plots #############
if doExtraPlots:


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
  #plot_Dist.addOvfl         = "no"
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
  #plot_DCotTheta.addOvfl         = "no"
  plot_DCotTheta.addZUncBand     = zUncBand
  plot_DCotTheta.makeRatio       = makeRatio
  plot_DCotTheta.xbins           = [0,0.05,0.1,0.15,0.2,0.3,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
  plot_DCotTheta.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMETEle_fullSel ---
  variableName = "h1_mDeltaPhiMETEle_fullSel"

  plot100 = Plot()
  ## inputs for stacked histograms
  plot100.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot100.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot100.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot100.keys            = keys
  plot100.xtit            = "#Delta#phi(MET,e) (rad.) - fullSel"
  plot100.ytit            = "Number of events"
  #plot100.xlog            = "yes"
  plot100.ylog            = "no"
  plot100.rebin           = 20
  #plot100.xmin            = 0
  #plot100.xmax            = 3.14
  plot100.ymin            = 0
  plot100.ymax            = 20
  #plot100.lpos = "bottom-center"
  plot100.name            = "mDeltaPhiMETEle_presel_fullSel"
  plot100.addZUncBand     = zUncBand
  plot100.makeRatio       = makeRatio
  plot100.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET1stJet_fullSel ---
  variableName = "h1_mDeltaPhiMET1stJet_fullSel"

  plot101 = Plot()
  ## inputs for stacked histograms
  plot101.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot101.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot101.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot101.keys            = keys
  plot101.xtit            = "#Delta#phi(MET,1^{st} jet) (rad.) - fullSel"
  plot101.ytit            = "Number of events"
  #plot101.xlog            = "yes"
  plot101.ylog            = "no"
  plot101.rebin           = 20
  #plot101.xmin            = 0
  #plot101.xmax            = 3.146
  plot101.ymin            = 0
  plot101.ymax            = 20
  #plot101.lpos = "bottom-center"
  plot101.name            = "mDeltaPhiMET1stJet_presel_fullSel"
  plot101.addZUncBand     = zUncBand
  plot101.makeRatio       = makeRatio
  plot101.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_fullSel ---
  variableName = "h1_mDeltaPhiMET2ndJet_fullSel"

  plot102 = Plot()
  ## inputs for stacked histograms
  plot102.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot102.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot102.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot102.keys            = keys
  plot102.xtit            = "#Delta#phi(MET,2^{nd} jet) (rad.) - fullSel"
  plot102.ytit            = "Number of events"
  #plot102.xlog            = "yes"
  plot102.ylog            = "no"
  plot102.rebin           = 20
  #plot102.xmin            = 0
  #plot102.xmax            = 3.146
  plot102.ymin            = 0
  plot102.ymax            = 20
  #plot102.lpos = "bottom-center"
  plot102.name            = "mDeltaPhiMET2ndJet_presel_fullSel"
  plot102.addZUncBand     = zUncBand
  plot102.makeRatio       = makeRatio
  plot102.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_minDRej_fullSel ---
  variableName = "h1_minDRej_fullSel"

  plot103 = Plot()
  ## inputs for stacked histograms
  plot103.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot103.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot103.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot103.keys            = keys
  plot103.xtit            = "min#DeltaR(e,jets) - fullSel"
  plot103.ytit            = "Number of events"
  plot103.ylog            = "no"
  plot103.rebin           = 10
  plot103.xmin            = 0
  plot103.xmax            = 7
  plot103.ymin            = 0
  plot103.ymax            = 25
  #plot103.lpos = "bottom-center"
  plot103.name            = "minDRej_presel_fullSel"
  plot103.addZUncBand     = zUncBand
  plot103.makeRatio       = makeRatio
  #plot103.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
  plot103.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DeltaRjets_PAS_fullSel ---
  variableName = "h1_DeltaRjets_PAS_fullSel"

  plot104 = Plot()
  ## inputs for stacked histograms
  plot104.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot104.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot104.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot104.keys            = keys
  plot104.xtit            = "#DeltaR(j1,j2) - fullSel"
  plot104.ytit            = "Number of events"
  plot104.ylog            = "no"
  plot104.rebin           = 10
  plot104.xmin            = 0
  plot104.xmax            = 7
  plot104.ymin            = 0
  plot104.ymax            = 20
  #plot104.lpos = "bottom-center"
  plot104.name            = "DRjets_presel_fullSel"
  plot104.addZUncBand     = zUncBand
  plot104.makeRatio       = makeRatio
  #plot104.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
  plot104.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_fullSel ---
  variableName = "h1_Njet_fullSel"

  plot105 = Plot()
  ## inputs for stacked histograms
  plot105.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot105.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot105.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot105.keys            = keys
  plot105.xtit            = "Number of jets - fullSel"
  plot105.ytit            = "Number of events"
  plot105.ylog            = "no"
  plot105.rebin           = 1
  plot105.xmin            = -0.5
  plot105.xmax            = 11.5
  plot105.ymin            = 0
  plot105.ymax            = 20
  #plot105.lpos = "bottom-center"
  plot105.name            = "nJet_presel_fullSel"
  plot105.addZUncBand     = zUncBand
  plot105.makeRatio       = makeRatio
  plot105.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_fullSel ---
  variableName = "h1_NjetTCHELBTag_fullSel"

  plot106 = Plot()
  ## inputs for stacked histograms
  plot106.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot106.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot106.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot106.keys            = keys
  plot106.xtit            = "Number of b-tagged jets (TCHEL) - fullSel"
  plot106.ytit            = "Number of events"
  plot106.ylog            = "no"
  plot106.rebin           = 1
  plot106.xmin            = -0.5
  plot106.xmax            = 11.5
  plot106.ymin            = 0
  plot106.ymax            = 20
  #plot106.lpos = "bottom-center"
  plot106.name            = "nJetTCHELBTag_presel_fullSel"
  plot106.addZUncBand     = zUncBand
  plot106.makeRatio       = makeRatio
  plot106.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Vtxd01stEle_PAS_fullSel ---
  variableName = "h1_Vtxd01stEle_PAS_fullSel"

  plot107 = Plot()
  ## inputs for stacked histograms
  plot107.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot107.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot107.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot107.keys            = keys
  plot107.xtit            = "Vtx d0 1st electron [cm] - fullSel"
  plot107.ytit            = "Number of events"
  plot107.ylog            = "no"
  plot107.rebin           = 10
  #plot107.xmin            = -0.01
  #plot107.xmax            = 0.01
  plot107.ymin            = 0
  plot107.ymax            = 8
  #plot107.lpos = "bottom-center"
  plot107.name            = "Vtxd01stEle_presel_fullSel"
  plot107.addZUncBand     = zUncBand
  plot107.makeRatio       = makeRatio
  plot107.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MissingHits1stEle_PAS_fullSel ---
  variableName = "h1_MissingHits1stEle_PAS_fullSel"

  plot108 = Plot()
  ## inputs for stacked histograms
  plot108.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot108.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot108.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot108.keys            = keys
  plot108.xtit            = "Missing Hits 1st electron - fullSel"
  plot108.ytit            = "Number of events"
  plot108.ylog            = "no"
  plot108.rebin           = 1
  #plot108.xmin            = -0.01
  #plot108.xmax            = 0.01
  plot108.ymin            = 0
  plot108.ymax            = 20
  #plot108.lpos = "bottom-center"
  plot108.name            = "MissingHits1stEle_presel_fullSel"
  plot108.addZUncBand     = zUncBand
  plot108.makeRatio       = makeRatio
  plot108.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Dist1stEle_PAS_fullSel ---
  variableName = "h1_Dist1stEle_PAS_fullSel"

  plot109 = Plot()
  ## inputs for stacked histograms
  plot109.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot109.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot109.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot109.keys            = keys
  plot109.xtit            = "|Dist| 1st electron - fullSel"
  plot109.ytit            = "Number of events"
  plot109.ylog            = "no"
  plot109.rebin           = 2
  plot109.xmin            = 0.
  plot109.xmax            = 0.6
  plot109.ymin            = 0
  plot109.ymax            = 10
  #plot109.lpos = "bottom-center"
  plot109.name            = "Dist1stEle_presel_fullSel"
  plot109.addZUncBand     = zUncBand
  plot109.makeRatio       = makeRatio
  plot109.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DCotTheta1stEle_PAS_fullSel ---
  variableName = "h1_DCotTheta1stEle_PAS_fullSel"

  plot110 = Plot()
  ## inputs for stacked histograms
  plot110.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot110.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot110.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot110.keys            = keys
  plot110.xtit            = "|#DeltaCot(#theta)| 1st electron - fullSel"
  plot110.ytit            = "Number of events"
  plot110.ylog            = "no"
  plot110.rebin           = 2
  plot110.xmin            = 0.
  plot110.xmax            = 0.8
  plot110.ymin            = 0
  plot110.ymax            = 10
  #plot110.lpos = "bottom-center"
  plot110.name            = "DCotTheta1stEle_presel_fullSel"
  plot110.addZUncBand     = zUncBand
  plot110.makeRatio       = makeRatio
  plot110.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  # ############################ Plots below to be done after pre-selection + high MT ######################

  #--- h1_mDeltaPhiMETEle_highMT ---
  variableName = "h1_mDeltaPhiMETEle_highMT"

  plot200 = Plot()
  ## inputs for stacked histograms
  plot200.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot200.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot200.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot200.keys            = keys
  plot200.xtit            = "#Delta#phi(MET,e) (rad.) - highMT"
  plot200.ytit            = "Number of events"
  #plot200.xlog            = "yes"
  plot200.ylog            = "no"
  plot200.rebin           = 20
  #plot200.xmin            = 0
  #plot200.xmax            = 3.14
  plot200.ymin            = 0
  plot200.ymax            = 25
  #plot200.lpos = "bottom-center"
  plot200.name            = "mDeltaPhiMETEle_presel_highMT"
  plot200.addZUncBand     = zUncBand
  plot200.makeRatio       = makeRatio
  plot200.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET1stJet_highMT ---
  variableName = "h1_mDeltaPhiMET1stJet_highMT"

  plot201 = Plot()
  ## inputs for stacked histograms
  plot201.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot201.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot201.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot201.keys            = keys
  plot201.xtit            = "#Delta#phi(MET,1^{st} jet) (rad.) - highMT"
  plot201.ytit            = "Number of events"
  #plot201.xlog            = "yes"
  plot201.ylog            = "no"
  plot201.rebin           = 20
  #plot201.xmin            = 0
  #plot201.xmax            = 3.146
  plot201.ymin            = 0
  plot201.ymax            = 20
  #plot201.lpos = "bottom-center"
  plot201.name            = "mDeltaPhiMET1stJet_presel_highMT"
  plot201.addZUncBand     = zUncBand
  plot201.makeRatio       = makeRatio
  plot201.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_highMT ---
  variableName = "h1_mDeltaPhiMET2ndJet_highMT"

  plot202 = Plot()
  ## inputs for stacked histograms
  plot202.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot202.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot202.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot202.keys            = keys
  plot202.xtit            = "#Delta#phi(MET,2^{nd} jet) (rad.) - highMT"
  plot202.ytit            = "Number of events"
  #plot202.xlog            = "yes"
  plot202.ylog            = "no"
  plot202.rebin           = 20
  #plot202.xmin            = 0
  #plot202.xmax            = 3.1416
  plot202.ymin            = 0
  plot202.ymax            = 20
  #plot202.lpos = "bottom-center"
  plot202.name            = "mDeltaPhiMET2ndJet_presel_highMT"
  plot202.addZUncBand     = zUncBand
  plot202.makeRatio       = makeRatio
  plot202.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_minDRej_highMT ---
  variableName = "h1_minDRej_highMT"

  plot203 = Plot()
  ## inputs for stacked histograms
  plot203.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot203.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot203.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot203.keys            = keys
  plot203.xtit            = "min#DeltaR(e,jets) - highMT"
  plot203.ytit            = "Number of events"
  plot203.ylog            = "no"
  plot203.rebin           = 10
  plot203.xmin            = 0
  plot203.xmax            = 7
  plot203.ymin            = 0
  plot203.ymax            = 25
  #plot203.lpos = "bottom-center"
  plot203.name            = "minDRej_presel_highMT"
  plot203.addZUncBand     = zUncBand
  plot203.makeRatio       = makeRatio
  #plot203.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
  plot203.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DeltaRjets_PAS_highMT ---
  variableName = "h1_DeltaRjets_PAS_highMT"

  plot204 = Plot()
  ## inputs for stacked histograms
  plot204.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot204.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot204.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot204.keys            = keys
  plot204.xtit            = "#DeltaR(j1,j2) - highMT"
  plot204.ytit            = "Number of events"
  plot204.ylog            = "no"
  plot204.rebin           = 10
  plot204.xmin            = 0
  plot204.xmax            = 7
  plot204.ymin            = 0
  plot204.ymax            = 20
  #plot204.lpos = "bottom-center"
  plot204.name            = "DRjets_presel_highMT"
  plot204.addZUncBand     = zUncBand
  plot204.makeRatio       = makeRatio
  #plot204.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
  plot204.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_highMT ---
  variableName = "h1_Njet_highMT"

  plot205 = Plot()
  ## inputs for stacked histograms
  plot205.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot205.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot205.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot205.keys            = keys
  plot205.xtit            = "Number of jets - highMT"
  plot205.ytit            = "Number of events"
  plot205.ylog            = "no"
  plot205.rebin           = 1
  plot205.xmin            = -0.5
  plot205.xmax            = 11.5
  plot205.ymin            = 0
  plot205.ymax            = 20
  #plot205.lpos = "bottom-center"
  plot205.name            = "nJet_presel_highMT"
  plot205.addZUncBand     = zUncBand
  plot205.makeRatio       = makeRatio
  plot205.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_highMT ---
  variableName = "h1_NjetTCHELBTag_highMT"

  plot206 = Plot()
  ## inputs for stacked histograms
  plot206.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot206.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot206.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot206.keys            = keys
  plot206.xtit            = "Number of b-tagged jets (TCHEL) - highMT"
  plot206.ytit            = "Number of events"
  plot206.ylog            = "no"
  plot206.rebin           = 1
  plot206.xmin            = -0.5
  plot206.xmax            = 11.5
  plot206.ymin            = 0
  plot206.ymax            = 20
  #plot206.lpos = "bottom-center"
  plot206.name            = "nJetTCHELBTag_presel_highMT"
  plot206.addZUncBand     = zUncBand
  plot206.makeRatio       = makeRatio
  plot206.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Vtxd01stEle_PAS_highMT ---
  variableName = "h1_Vtxd01stEle_PAS_highMT"

  plot207 = Plot()
  ## inputs for stacked histograms
  plot207.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot207.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot207.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot207.keys            = keys
  plot207.xtit            = "Vtx d0 1st electron [cm] - highMT"
  plot207.ytit            = "Number of events"
  plot207.ylog            = "no"
  plot207.rebin           = 10
  #plot207.xmin            = -0.01
  #plot207.xmax            = 0.01
  plot207.ymin            = 0
  plot207.ymax            = 8
  #plot207.lpos = "bottom-center"
  plot207.name            = "Vtxd01stEle_presel_highMT"
  plot207.addZUncBand     = zUncBand
  plot207.makeRatio       = makeRatio
  plot207.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MissingHits1stEle_PAS_highMT ---
  variableName = "h1_MissingHits1stEle_PAS_highMT"

  plot208 = Plot()
  ## inputs for stacked histograms
  plot208.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot208.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot208.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot208.keys            = keys
  plot208.xtit            = "Missing Hits 1st electron - highMT"
  plot208.ytit            = "Number of events"
  plot208.ylog            = "no"
  plot208.rebin           = 1
  #plot208.xmin            = -0.01
  #plot208.xmax            = 0.01
  plot208.ymin            = 0
  plot208.ymax            = 20
  #plot208.lpos = "bottom-center"
  plot208.name            = "MissingHits1stEle_presel_highMT"
  plot208.addZUncBand     = zUncBand
  plot208.makeRatio       = makeRatio
  plot208.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Dist1stEle_PAS_highMT ---
  variableName = "h1_Dist1stEle_PAS_highMT"

  plot209 = Plot()
  ## inputs for stacked histograms
  plot209.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot209.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot209.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot209.keys            = keys
  plot209.xtit            = "|Dist| 1st electron - highMT"
  plot209.ytit            = "Number of events"
  plot209.ylog            = "no"
  plot209.rebin           = 2
  plot209.xmin            = 0.
  plot209.xmax            = 0.6
  plot209.ymin            = 0
  plot209.ymax            = 10
  #plot209.lpos = "bottom-center"
  plot209.name            = "Dist1stEle_presel_highMT"
  plot209.addZUncBand     = zUncBand
  plot209.makeRatio       = makeRatio
  plot209.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DCotTheta1stEle_PAS_highMT ---
  variableName = "h1_DCotTheta1stEle_PAS_highMT"

  plot210 = Plot()
  ## inputs for stacked histograms
  plot210.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot210.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot210.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot210.keys            = keys
  plot210.xtit            = "|#DeltaCot(#theta)| 1st electron - highMT"
  plot210.ytit            = "Number of events"
  plot210.ylog            = "no"
  plot210.rebin           = 2
  plot210.xmin            = 0.
  plot210.xmax            = 0.8
  plot210.ymin            = 0
  plot210.ymax            = 10
  #plot210.lpos = "bottom-center"
  plot210.name            = "DCotTheta1stEle_presel_highMT"
  plot210.addZUncBand     = zUncBand
  plot210.makeRatio       = makeRatio
  plot210.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  # ############################ Plots below to be done after pre-selection + high Pt1stEle ######################

  #--- h1_mDeltaPhiMETEle_highPt1stEle ---
  variableName = "h1_mDeltaPhiMETEle_highPt1stEle"

  plot300 = Plot()
  ## inputs for stacked histograms
  plot300.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot300.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot300.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot300.keys            = keys
  plot300.xtit            = "#Delta#phi(MET,e) (rad.) - highPt1stEle"
  plot300.ytit            = "Number of events"
  #plot300.xlog            = "yes"
  plot300.ylog            = "no"
  plot300.rebin           = 20
  #plot300.xmin            = 0
  #plot300.xmax            = 3.14
  plot300.ymin            = 0
  plot300.ymax            = 10
  #plot300.lpos = "bottom-center"
  plot300.name            = "mDeltaPhiMETEle_presel_highPt1stEle"
  plot300.addZUncBand     = zUncBand
  plot300.makeRatio       = makeRatio
  plot300.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET1stJet_highPt1stEle ---
  variableName = "h1_mDeltaPhiMET1stJet_highPt1stEle"

  plot301 = Plot()
  ## inputs for stacked histograms
  plot301.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot301.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot301.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot301.keys            = keys
  plot301.xtit            = "#Delta#phi(MET,1^{st} jet) (rad.) - highPt1stEle"
  plot301.ytit            = "Number of events"
  #plot301.xlog            = "yes"
  plot301.ylog            = "no"
  plot301.rebin           = 20
  #plot301.xmin            = 0
  #plot301.xmax            = 3.146
  plot301.ymin            = 0
  plot301.ymax            = 10
  #plot301.lpos = "bottom-center"
  plot301.name            = "mDeltaPhiMET1stJet_presel_highPt1stEle"
  plot301.addZUncBand     = zUncBand
  plot301.makeRatio       = makeRatio
  plot301.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_highPt1stEle ---
  variableName = "h1_mDeltaPhiMET2ndJet_highPt1stEle"

  plot302 = Plot()
  ## inputs for stacked histograms
  plot302.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot302.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot302.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot302.keys            = keys
  plot302.xtit            = "#Delta#phi(MET,2^{nd} jet) (rad.) - highPt1stEle"
  plot302.ytit            = "Number of events"
  #plot302.xlog            = "yes"
  plot302.ylog            = "no"
  plot302.rebin           = 20
  #plot302.xmin            = 0
  #plot302.xmax            = 3.146
  plot302.ymin            = 0
  plot302.ymax            = 10
  #plot302.lpos = "bottom-center"
  plot302.name            = "mDeltaPhiMET2ndJet_presel_highPt1stEle"
  plot302.addZUncBand     = zUncBand
  plot302.makeRatio       = makeRatio
  plot302.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_minDRej_highPt1stEle ---
  variableName = "h1_minDRej_highPt1stEle"

  plot303 = Plot()
  ## inputs for stacked histograms
  plot303.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot303.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot303.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot303.keys            = keys
  plot303.xtit            = "min#DeltaR(e,jets) - highPt1stEle"
  plot303.ytit            = "Number of events"
  plot303.ylog            = "no"
  plot303.rebin           = 10
  plot303.xmin            = 0
  plot303.xmax            = 7
  plot303.ymin            = 0
  plot303.ymax            = 10
  #plot303.lpos = "bottom-center"
  plot303.name            = "minDRej_presel_highPt1stEle"
  plot303.addZUncBand     = zUncBand
  plot303.makeRatio       = makeRatio
  #plot303.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
  plot303.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DeltaRjets_PAS_highPt1stEle ---
  variableName = "h1_DeltaRjets_PAS_highPt1stEle"

  plot304 = Plot()
  ## inputs for stacked histograms
  plot304.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot304.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot304.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot304.keys            = keys
  plot304.xtit            = "#DeltaR(j1,j2) - highPt1stEle"
  plot304.ytit            = "Number of events"
  plot304.ylog            = "no"
  plot304.rebin           = 10
  plot304.xmin            = 0
  plot304.xmax            = 7
  plot304.ymin            = 0
  plot304.ymax            = 10
  #plot304.lpos = "bottom-center"
  plot304.name            = "DRjets_presel_highPt1stEle"
  plot304.addZUncBand     = zUncBand
  plot304.makeRatio       = makeRatio
  #plot304.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
  plot304.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_highPt1stEle ---
  variableName = "h1_Njet_highPt1stEle"

  plot305 = Plot()
  ## inputs for stacked histograms
  plot305.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot305.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot305.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot305.keys            = keys
  plot305.xtit            = "Number of jets - highPt1stEle"
  plot305.ytit            = "Number of events"
  plot305.ylog            = "no"
  plot305.rebin           = 1
  plot305.xmin            = -0.5
  plot305.xmax            = 11.5
  plot305.ymin            = 0
  plot305.ymax            = 10
  #plot305.lpos = "bottom-center"
  plot305.name            = "nJet_presel_highPt1stEle"
  plot305.addZUncBand     = zUncBand
  plot305.makeRatio       = makeRatio
  plot305.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_highPt1stEle ---
  variableName = "h1_NjetTCHELBTag_highPt1stEle"

  plot306 = Plot()
  ## inputs for stacked histograms
  plot306.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot306.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot306.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot306.keys            = keys
  plot306.xtit            = "Number of b-tagged jets (TCHEL) - highPt1stEle"
  plot306.ytit            = "Number of events"
  plot306.ylog            = "no"
  plot306.rebin           = 1
  plot306.xmin            = -0.5
  plot306.xmax            = 11.5
  plot306.ymin            = 0
  plot306.ymax            = 10
  #plot306.lpos = "bottom-center"
  plot306.name            = "nJetTCHELBTag_presel_highPt1stEle"
  plot306.addZUncBand     = zUncBand
  plot306.makeRatio       = makeRatio
  plot306.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Vtxd01stEle_PAS_highPt1stEle ---
  variableName = "h1_Vtxd01stEle_PAS_highPt1stEle"

  plot307 = Plot()
  ## inputs for stacked histograms
  plot307.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot307.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot307.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot307.keys            = keys
  plot307.xtit            = "Vtx d0 1st electron [cm] - highPt1stEle"
  plot307.ytit            = "Number of events"
  plot307.ylog            = "no"
  plot307.rebin           = 10
  #plot307.xmin            = -0.01
  #plot307.xmax            = 0.01
  plot307.ymin            = 0
  plot307.ymax            = 10
  #plot307.lpos = "bottom-center"
  plot307.name            = "Vtxd01stEle_presel_highPt1stEle"
  plot307.addZUncBand     = zUncBand
  plot307.makeRatio       = makeRatio
  plot307.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MissingHits1stEle_PAS_highPt1stEle ---
  variableName = "h1_MissingHits1stEle_PAS_highPt1stEle"

  plot308 = Plot()
  ## inputs for stacked histograms
  plot308.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot308.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot308.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot308.keys            = keys
  plot308.xtit            = "Missing Hits 1st electron - highPt1stEle"
  plot308.ytit            = "Number of events"
  plot308.ylog            = "no"
  plot308.rebin           = 1
  #plot308.xmin            = -0.01
  #plot308.xmax            = 0.01
  plot308.ymin            = 0
  plot308.ymax            = 10
  #plot308.lpos = "bottom-center"
  plot308.name            = "MissingHits1stEle_presel_highPt1stEle"
  plot308.addZUncBand     = zUncBand
  plot308.makeRatio       = makeRatio
  plot308.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Dist1stEle_PAS_highPt1stEle ---
  variableName = "h1_Dist1stEle_PAS_highPt1stEle"

  plot309 = Plot()
  ## inputs for stacked histograms
  plot309.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot309.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot309.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot309.keys            = keys
  plot309.xtit            = "|Dist| 1st electron - highPt1stEle"
  plot309.ytit            = "Number of events"
  plot309.ylog            = "no"
  plot309.rebin           = 2
  plot309.xmin            = 0.
  plot309.xmax            = 0.6
  plot309.ymin            = 0
  plot309.ymax            = 5
  #plot309.lpos = "bottom-center"
  plot309.name            = "Dist1stEle_presel_highPt1stEle"
  plot309.addZUncBand     = zUncBand
  plot309.makeRatio       = makeRatio
  plot309.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DCotTheta1stEle_PAS_highPt1stEle ---
  variableName = "h1_DCotTheta1stEle_PAS_highPt1stEle"

  plot310 = Plot()
  ## inputs for stacked histograms
  plot310.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot310.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot310.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot310.keys            = keys
  plot310.xtit            = "|#DeltaCot(#theta)| 1st electron - highPt1stEle"
  plot310.ytit            = "Number of events"
  plot310.ylog            = "no"
  plot310.rebin           = 2
  plot310.xmin            = 0.
  plot310.xmax            = 0.8
  plot310.ymin            = 0
  plot310.ymax            = 5
  #plot310.lpos = "bottom-center"
  plot310.name            = "DCotTheta1stEle_presel_highPt1stEle"
  plot310.addZUncBand     = zUncBand
  plot310.makeRatio       = makeRatio
  plot310.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  # ############################ Plots below to be done after pre-selection + high Pt1stJet ######################

  #--- h1_mDeltaPhiMETEle_highPt1stEle ---
  variableName = "h1_mDeltaPhiMETEle_highPt1stJet"

  plot400 = Plot()
  ## inputs for stacked histograms
  plot400.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot400.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot400.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot400.keys            = keys
  plot400.xtit            = "#Delta#phi(MET,e) (rad.) - highPt1stJet"
  plot400.ytit            = "Number of events"
  #plot400.xlog            = "yes"
  plot400.ylog            = "no"
  plot400.rebin           = 20
  #plot400.xmin            = 0
  #plot400.xmax            = 3.14
  plot400.ymin            = 0
  plot400.ymax            = 10
  #plot400.lpos = "bottom-center"
  plot400.name            = "mDeltaPhiMETEle_presel_highPt1stJet"
  plot400.addZUncBand     = zUncBand
  plot400.makeRatio       = makeRatio
  plot400.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET1stJet_highPt1stJet ---
  variableName = "h1_mDeltaPhiMET1stJet_highPt1stJet"

  plot401 = Plot()
  ## inputs for stacked histograms
  plot401.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot401.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot401.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot401.keys            = keys
  plot401.xtit            = "#Delta#phi(MET,1^{st} jet) (rad.) - highPt1stJet"
  plot401.ytit            = "Number of events"
  #plot401.xlog            = "yes"
  plot401.ylog            = "no"
  plot401.rebin           = 20
  #plot401.xmin            = 0
  #plot401.xmax            = 3.146
  plot401.ymin            = 0
  plot401.ymax            = 15
  #plot401.lpos = "bottom-center"
  plot401.name            = "mDeltaPhiMET1stJet_presel_highPt1stJet"
  plot401.addZUncBand     = zUncBand
  plot401.makeRatio       = makeRatio
  plot401.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_highPt1stJet ---
  variableName = "h1_mDeltaPhiMET2ndJet_highPt1stJet"

  plot402 = Plot()
  ## inputs for stacked histograms
  plot402.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot402.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot402.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot402.keys            = keys
  plot402.xtit            = "#Delta#phi(MET,2^{nd} jet) (rad.) - highPt1stJet"
  plot402.ytit            = "Number of events"
  #plot402.xlog            = "yes"
  plot402.ylog            = "no"
  plot402.rebin           = 20
  #plot402.xmin            = 0
  #plot402.xmax            = 3.146
  plot402.ymin            = 0
  plot402.ymax            = 10
  #plot402.lpos = "bottom-center"
  plot402.name            = "mDeltaPhiMET2ndJet_presel_highPt1stJet"
  plot402.addZUncBand     = zUncBand
  plot402.makeRatio       = makeRatio
  plot402.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_minDRej_highPt1stJet ---
  variableName = "h1_minDRej_highPt1stJet"

  plot403 = Plot()
  ## inputs for stacked histograms
  plot403.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot403.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot403.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot403.keys            = keys
  plot403.xtit            = "min#DeltaR(e,jets) - highPt1stJet"
  plot403.ytit            = "Number of events"
  plot403.ylog            = "no"
  plot403.rebin           = 10
  plot403.xmin            = 0
  plot403.xmax            = 7
  plot403.ymin            = 0
  plot403.ymax            = 10
  #plot403.lpos = "bottom-center"
  plot403.name            = "minDRej_presel_highPt1stJet"
  plot403.addZUncBand     = zUncBand
  plot403.makeRatio       = makeRatio
  #plot403.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
  plot403.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DeltaRjets_PAS_highPt1stJet ---
  variableName = "h1_DeltaRjets_PAS_highPt1stJet"

  plot404 = Plot()
  ## inputs for stacked histograms
  plot404.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot404.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot404.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot404.keys            = keys
  plot404.xtit            = "#DeltaR(j1,j2) - highPt1stJet"
  plot404.ytit            = "Number of events"
  plot404.ylog            = "no"
  plot404.rebin           = 10
  plot404.xmin            = 0
  plot404.xmax            = 7
  plot404.ymin            = 0
  plot404.ymax            = 10
  #plot404.lpos = "bottom-center"
  plot404.name            = "DRjets_presel_highPt1stJet"
  plot404.addZUncBand     = zUncBand
  plot404.makeRatio       = makeRatio
  #plot404.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
  plot404.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_highPt1stJet ---
  variableName = "h1_Njet_highPt1stJet"

  plot405 = Plot()
  ## inputs for stacked histograms
  plot405.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot405.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot405.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot405.keys            = keys
  plot405.xtit            = "Number of jets - highPt1stJet"
  plot405.ytit            = "Number of events"
  plot405.ylog            = "no"
  plot405.rebin           = 1
  plot405.xmin            = -0.5
  plot405.xmax            = 11.5
  plot405.ymin            = 0
  plot405.ymax            = 10
  #plot405.lpos = "bottom-center"
  plot405.name            = "nJet_presel_highPt1stJet"
  plot405.addZUncBand     = zUncBand
  plot405.makeRatio       = makeRatio
  plot405.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_highPt1stJet ---
  variableName = "h1_NjetTCHELBTag_highPt1stJet"

  plot406 = Plot()
  ## inputs for stacked histograms
  plot406.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot406.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot406.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot406.keys            = keys
  plot406.xtit            = "Number of b-tagged jets (TCHEL) - highPt1stJet"
  plot406.ytit            = "Number of events"
  plot406.ylog            = "no"
  plot406.rebin           = 1
  plot406.xmin            = -0.5
  plot406.xmax            = 11.5
  plot406.ymin            = 0
  plot406.ymax            = 10
  #plot406.lpos = "bottom-center"
  plot406.name            = "nJetTCHELBTag_presel_highPt1stJet"
  plot406.addZUncBand     = zUncBand
  plot406.makeRatio       = makeRatio
  plot406.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Vtxd01stEle_PAS_highPt1stJet ---
  variableName = "h1_Vtxd01stEle_PAS_highPt1stJet"

  plot407 = Plot()
  ## inputs for stacked histograms
  plot407.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot407.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot407.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot407.keys            = keys
  plot407.xtit            = "Vtx d0 1st electron [cm] - highPt1stJet"
  plot407.ytit            = "Number of events"
  plot407.ylog            = "no"
  plot407.rebin           = 10
  #plot407.xmin            = -0.01
  #plot407.xmax            = 0.01
  plot407.ymin            = 0
  plot407.ymax            = 8
  #plot407.lpos = "bottom-center"
  plot407.name            = "Vtxd01stEle_presel_highPt1stJet"
  plot407.addZUncBand     = zUncBand
  plot407.makeRatio       = makeRatio
  plot407.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MissingHits1stEle_PAS_highPt1stJet ---
  variableName = "h1_MissingHits1stEle_PAS_highPt1stJet"

  plot408 = Plot()
  ## inputs for stacked histograms
  plot408.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot408.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot408.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot408.keys            = keys
  plot408.xtit            = "Missing Hits 1st electron - highPt1stJet"
  plot408.ytit            = "Number of events"
  plot408.ylog            = "no"
  plot408.rebin           = 1
  #plot408.xmin            = -0.01
  #plot408.xmax            = 0.01
  plot408.ymin            = 0
  plot408.ymax            = 20
  #plot408.lpos = "bottom-center"
  plot408.name            = "MissingHits1stEle_presel_highPt1stJet"
  plot408.addZUncBand     = zUncBand
  plot408.makeRatio       = makeRatio
  plot408.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Dist1stEle_PAS_highPt1stJet ---
  variableName = "h1_Dist1stEle_PAS_highPt1stJet"

  plot409 = Plot()
  ## inputs for stacked histograms
  plot409.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot409.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot409.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot409.keys            = keys
  plot409.xtit            = "|Dist| 1st electron - highPt1stJet"
  plot409.ytit            = "Number of events"
  plot409.ylog            = "no"
  plot409.rebin           = 2
  plot409.xmin            = 0.
  plot409.xmax            = 0.6
  plot409.ymin            = 0
  plot409.ymax            = 15
  #plot409.lpos = "bottom-center"
  plot409.name            = "Dist1stEle_presel_highPt1stJet"
  plot409.addZUncBand     = zUncBand
  plot409.makeRatio       = makeRatio
  plot409.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DCotTheta1stEle_PAS_highPt1stJet ---
  variableName = "h1_DCotTheta1stEle_PAS_highPt1stJet"

  plot410 = Plot()
  ## inputs for stacked histograms
  plot410.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot410.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot410.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot410.keys            = keys
  plot410.xtit            = "|#DeltaCot(#theta)| 1st electron - highPt1stJet"
  plot410.ytit            = "Number of events"
  plot410.ylog            = "no"
  plot410.rebin           = 2
  plot410.xmin            = 0.
  plot410.xmax            = 0.8
  plot410.ymin            = 0
  plot410.ymax            = 10
  #plot410.lpos = "bottom-center"
  plot410.name            = "DCotTheta1stEle_presel_highPt1stJet"
  plot410.addZUncBand     = zUncBand
  plot410.makeRatio       = makeRatio
  plot410.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  # ############################ Plots below to be done after pre-selection + high Pt2ndJet ######################

  #--- h1_mDeltaPhiMETEle_highPt1stEle ---
  variableName = "h1_mDeltaPhiMETEle_highPt2ndJet"

  plot500 = Plot()
  ## inputs for stacked histograms
  plot500.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot500.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot500.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot500.keys            = keys
  plot500.xtit            = "#Delta#phi(MET,e) (rad.) - highPt2ndJet"
  plot500.ytit            = "Number of events"
  #plot500.xlog            = "yes"
  plot500.ylog            = "no"
  plot500.rebin           = 20
  #plot500.xmin            = 0
  #plot500.xmax            = 3.14
  plot500.ymin            = 0
  plot500.ymax            = 10
  #plot500.lpos = "bottom-center"
  plot500.name            = "mDeltaPhiMETEle_presel_highPt2ndJet"
  plot500.addZUncBand     = zUncBand
  plot500.makeRatio       = makeRatio
  plot500.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET1stJet_highPt2ndJet ---
  variableName = "h1_mDeltaPhiMET1stJet_highPt2ndJet"

  plot501 = Plot()
  ## inputs for stacked histograms
  plot501.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot501.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot501.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot501.keys            = keys
  plot501.xtit            = "#Delta#phi(MET,1^{st} jet) (rad.) - highPt2ndJet"
  plot501.ytit            = "Number of events"
  #plot501.xlog            = "yes"
  plot501.ylog            = "no"
  plot501.rebin           = 20
  #plot501.xmin            = 0
  #plot501.xmax            = 3.146
  plot501.ymin            = 0
  plot501.ymax            = 15
  #plot501.lpos = "bottom-center"
  plot501.name            = "mDeltaPhiMET1stJet_presel_highPt2ndJet"
  plot501.addZUncBand     = zUncBand
  plot501.makeRatio       = makeRatio
  plot501.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_highPt2ndJet ---
  variableName = "h1_mDeltaPhiMET2ndJet_highPt2ndJet"

  plot502 = Plot()
  ## inputs for stacked histograms
  plot502.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot502.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot502.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot502.keys            = keys
  plot502.xtit            = "#Delta#phi(MET,2^{nd} jet) (rad.) - highPt2ndJet"
  plot502.ytit            = "Number of events"
  #plot502.xlog            = "yes"
  plot502.ylog            = "no"
  plot502.rebin           = 20
  #plot502.xmin            = 0
  #plot502.xmax            = 3.146
  plot502.ymin            = 0
  plot502.ymax            = 10
  #plot502.lpos = "bottom-center"
  plot502.name            = "mDeltaPhiMET2ndJet_presel_highPt2ndJet"
  plot502.addZUncBand     = zUncBand
  plot502.makeRatio       = makeRatio
  plot502.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_minDRej_highPt2ndJet ---
  variableName = "h1_minDRej_highPt2ndJet"

  plot503 = Plot()
  ## inputs for stacked histograms
  plot503.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot503.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot503.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot503.keys            = keys
  plot503.xtit            = "min#DeltaR(e,jets) - highPt2ndJet"
  plot503.ytit            = "Number of events"
  plot503.ylog            = "no"
  plot503.rebin           = 10
  plot503.xmin            = 0
  plot503.xmax            = 7
  plot503.ymin            = 0
  plot503.ymax            = 10
  #plot503.lpos = "bottom-center"
  plot503.name            = "minDRej_presel_highPt2ndJet"
  plot503.addZUncBand     = zUncBand
  plot503.makeRatio       = makeRatio
  #plot503.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
  plot503.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DeltaRjets_PAS_highPt2ndJet ---
  variableName = "h1_DeltaRjets_PAS_highPt2ndJet"

  plot504 = Plot()
  ## inputs for stacked histograms
  plot504.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot504.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot504.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot504.keys            = keys
  plot504.xtit            = "#DeltaR(j1,j2) - highPt2ndJet"
  plot504.ytit            = "Number of events"
  plot504.ylog            = "no"
  plot504.rebin           = 10
  plot504.xmin            = 0
  plot504.xmax            = 7
  plot504.ymin            = 0
  plot504.ymax            = 10
  #plot504.lpos = "bottom-center"
  plot504.name            = "DRjets_presel_highPt2ndJet"
  plot504.addZUncBand     = zUncBand
  plot504.makeRatio       = makeRatio
  #plot504.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
  plot504.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_highPt2ndJet ---
  variableName = "h1_Njet_highPt2ndJet"

  plot505 = Plot()
  ## inputs for stacked histograms
  plot505.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot505.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot505.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot505.keys            = keys
  plot505.xtit            = "Number of jets - highPt2ndJet"
  plot505.ytit            = "Number of events"
  plot505.ylog            = "no"
  plot505.rebin           = 1
  plot505.xmin            = -0.5
  plot505.xmax            = 11.5
  plot505.ymin            = 0
  plot505.ymax            = 10
  #plot505.lpos = "bottom-center"
  plot505.name            = "nJet_presel_highPt2ndJet"
  plot505.addZUncBand     = zUncBand
  plot505.makeRatio       = makeRatio
  plot505.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_highPt2ndJet ---
  variableName = "h1_NjetTCHELBTag_highPt2ndJet"

  plot506 = Plot()
  ## inputs for stacked histograms
  plot506.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot506.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot506.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot506.keys            = keys
  plot506.xtit            = "Number of b-tagged jets (TCHEL) - highPt2ndJet"
  plot506.ytit            = "Number of events"
  plot506.ylog            = "no"
  plot506.rebin           = 1
  plot506.xmin            = -0.5
  plot506.xmax            = 11.5
  plot506.ymin            = 0
  plot506.ymax            = 10
  #plot506.lpos = "bottom-center"
  plot506.name            = "nJetTCHELBTag_presel_highPt2ndJet"
  plot506.addZUncBand     = zUncBand
  plot506.makeRatio       = makeRatio
  plot506.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Vtxd01stEle_PAS_highPt2ndJet ---
  variableName = "h1_Vtxd01stEle_PAS_highPt2ndJet"

  plot507 = Plot()
  ## inputs for stacked histograms
  plot507.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot507.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot507.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot507.keys            = keys
  plot507.xtit            = "Vtx d0 1st electron [cm] - highPt2ndJet"
  plot507.ytit            = "Number of events"
  plot507.ylog            = "no"
  plot507.rebin           = 10
  #plot507.xmin            = -0.01
  #plot507.xmax            = 0.01
  plot507.ymin            = 0
  plot507.ymax            = 8
  #plot507.lpos = "bottom-center"
  plot507.name            = "Vtxd01stEle_presel_highPt2ndJet"
  plot507.addZUncBand     = zUncBand
  plot507.makeRatio       = makeRatio
  plot507.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MissingHits1stEle_PAS_highPt2ndJet ---
  variableName = "h1_MissingHits1stEle_PAS_highPt2ndJet"

  plot508 = Plot()
  ## inputs for stacked histograms
  plot508.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot508.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot508.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot508.keys            = keys
  plot508.xtit            = "Missing Hits 1st electron - highPt2ndJet"
  plot508.ytit            = "Number of events"
  plot508.ylog            = "no"
  plot508.rebin           = 1
  #plot508.xmin            = -0.01
  #plot508.xmax            = 0.01
  plot508.ymin            = 0
  plot508.ymax            = 20
  #plot508.lpos = "bottom-center"
  plot508.name            = "MissingHits1stEle_presel_highPt2ndJet"
  plot508.addZUncBand     = zUncBand
  plot508.makeRatio       = makeRatio
  plot508.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Dist1stEle_PAS_highPt2ndJet ---
  variableName = "h1_Dist1stEle_PAS_highPt2ndJet"

  plot509 = Plot()
  ## inputs for stacked histograms
  plot509.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot509.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot509.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot509.keys            = keys
  plot509.xtit            = "|Dist| 1st electron - highPt2ndJet"
  plot509.ytit            = "Number of events"
  plot509.ylog            = "no"
  plot509.rebin           = 2
  plot509.xmin            = 0.
  plot509.xmax            = 0.6
  plot509.ymin            = 0
  plot509.ymax            = 5
  #plot509.lpos = "bottom-center"
  plot509.name            = "Dist1stEle_presel_highPt2ndJet"
  plot509.addZUncBand     = zUncBand
  plot509.makeRatio       = makeRatio
  plot509.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DCotTheta1stEle_PAS_highPt2ndJet ---
  variableName = "h1_DCotTheta1stEle_PAS_highPt2ndJet"

  plot510 = Plot()
  ## inputs for stacked histograms
  plot510.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot510.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot510.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot510.keys            = keys
  plot510.xtit            = "|#DeltaCot(#theta)| 1st electron - highPt2ndJet"
  plot510.ytit            = "Number of events"
  plot510.ylog            = "no"
  plot510.rebin           = 2
  plot510.xmin            = 0.
  plot510.xmax            = 0.8
  plot510.ymin            = 0
  plot510.ymax            = 5
  #plot510.lpos = "bottom-center"
  plot510.name            = "DCotTheta1stEle_presel_highPt2ndJet"
  plot510.addZUncBand     = zUncBand
  plot510.makeRatio       = makeRatio
  plot510.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  # ############################ Plots below to be done after pre-selection + high Mej ######################

  #--- h1_mDeltaPhiMETEle_highPt1stEle ---
  variableName = "h1_mDeltaPhiMETEle_highMej"

  plot600 = Plot()
  ## inputs for stacked histograms
  plot600.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot600.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot600.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot600.keys            = keys
  plot600.xtit            = "#Delta#phi(MET,e) (rad.) - highMej"
  plot600.ytit            = "Number of events"
  #plot600.xlog            = "yes"
  plot600.ylog            = "no"
  plot600.rebin           = 20
  #plot600.xmin            = 0
  #plot600.xmax            = 3.14
  plot600.ymin            = 0
  plot600.ymax            = 50
  #plot600.lpos = "bottom-center"
  plot600.name            = "mDeltaPhiMETEle_presel_highMej"
  plot600.addZUncBand     = zUncBand
  plot600.makeRatio       = makeRatio
  plot600.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET1stJet_highMej ---
  variableName = "h1_mDeltaPhiMET1stJet_highMej"

  plot601 = Plot()
  ## inputs for stacked histograms
  plot601.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot601.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot601.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot601.keys            = keys
  plot601.xtit            = "#Delta#phi(MET,1^{st} jet) (rad.) - highMej"
  plot601.ytit            = "Number of events"
  #plot601.xlog            = "yes"
  plot601.ylog            = "no"
  plot601.rebin           = 5
  #plot601.xmin            = 0
  #plot601.xmax            = 3.146
  plot601.ymin            = 0
  plot601.ymax            = 25
  #plot601.lpos = "bottom-center"
  plot601.name            = "mDeltaPhiMET1stJet_presel_highMej"
  plot601.addZUncBand     = zUncBand
  plot601.makeRatio       = makeRatio
  plot601.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_highMej ---
  variableName = "h1_mDeltaPhiMET2ndJet_highMej"

  plot602 = Plot()
  ## inputs for stacked histograms
  plot602.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot602.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot602.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot602.keys            = keys
  plot602.xtit            = "#Delta#phi(MET,2^{nd} jet) (rad.) - highMej"
  plot602.ytit            = "Number of events"
  #plot602.xlog            = "yes"
  plot602.ylog            = "no"
  plot602.rebin           = 20
  #plot602.xmin            = 0
  #plot602.xmax            = 3.146
  plot602.ymin            = 0
  plot602.ymax            = 50
  #plot602.lpos = "bottom-center"
  plot602.name            = "mDeltaPhiMET2ndJet_presel_highMej"
  plot602.addZUncBand     = zUncBand
  plot602.makeRatio       = makeRatio
  plot602.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_minDRej_highMej ---
  variableName = "h1_minDRej_highMej"

  plot603 = Plot()
  ## inputs for stacked histograms
  plot603.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot603.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot603.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot603.keys            = keys
  plot603.xtit            = "min#DeltaR(e,jets) - highMej"
  plot603.ytit            = "Number of events"
  plot603.ylog            = "no"
  plot603.rebin           = 10
  plot603.xmin            = 0
  plot603.xmax            = 7
  plot603.ymin            = 0
  plot603.ymax            = 50
  #plot603.lpos = "bottom-center"
  plot603.name            = "minDRej_presel_highMej"
  plot603.addZUncBand     = zUncBand
  plot603.makeRatio       = makeRatio
  #plot603.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
  plot603.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DeltaRjets_PAS_highMej ---
  variableName = "h1_DeltaRjets_PAS_highMej"

  plot604 = Plot()
  ## inputs for stacked histograms
  plot604.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot604.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot604.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot604.keys            = keys
  plot604.xtit            = "#DeltaR(j1,j2) - highMej"
  plot604.ytit            = "Number of events"
  plot604.ylog            = "no"
  plot604.rebin           = 10
  plot604.xmin            = 0
  plot604.xmax            = 7
  plot604.ymin            = 0
  plot604.ymax            = 50
  #plot604.lpos = "bottom-center"
  plot604.name            = "DRjets_presel_highMej"
  plot604.addZUncBand     = zUncBand
  plot604.makeRatio       = makeRatio
  #plot604.xbins           = [0,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000]
  plot604.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_highMej ---
  variableName = "h1_Njet_highMej"

  plot605 = Plot()
  ## inputs for stacked histograms
  plot605.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot605.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot605.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot605.keys            = keys
  plot605.xtit            = "Number of jets - highMej"
  plot605.ytit            = "Number of events"
  plot605.ylog            = "no"
  plot605.rebin           = 1
  plot605.xmin            = -0.5
  plot605.xmax            = 11.5
  plot605.ymin            = 0
  plot605.ymax            = 50
  #plot605.lpos = "bottom-center"
  plot605.name            = "nJet_presel_highMej"
  plot605.addZUncBand     = zUncBand
  plot605.makeRatio       = makeRatio
  plot605.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_highMej ---
  variableName = "h1_NjetTCHELBTag_highMej"

  plot606 = Plot()
  ## inputs for stacked histograms
  plot606.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot606.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot606.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot606.keys            = keys
  plot606.xtit            = "Number of b-tagged jets (TCHEL) - highMej"
  plot606.ytit            = "Number of events"
  plot606.ylog            = "no"
  plot606.rebin           = 1
  plot606.xmin            = -0.5
  plot606.xmax            = 11.5
  plot606.ymin            = 0
  plot606.ymax            = 50
  #plot606.lpos = "bottom-center"
  plot606.name            = "nJetTCHELBTag_presel_highMej"
  plot606.addZUncBand     = zUncBand
  plot606.makeRatio       = makeRatio
  plot606.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Vtxd01stEle_PAS_highMej ---
  variableName = "h1_Vtxd01stEle_PAS_highMej"

  plot607 = Plot()
  ## inputs for stacked histograms
  plot607.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot607.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot607.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot607.keys            = keys
  plot607.xtit            = "Vtx d0 1st electron [cm] - highMej"
  plot607.ytit            = "Number of events"
  plot607.ylog            = "no"
  plot607.rebin           = 10
  #plot607.xmin            = -0.01
  #plot607.xmax            = 0.01
  plot607.ymin            = 0
  plot607.ymax            = 20
  #plot607.lpos = "bottom-center"
  plot607.name            = "Vtxd01stEle_presel_highMej"
  plot607.addZUncBand     = zUncBand
  plot607.makeRatio       = makeRatio
  plot607.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MissingHits1stEle_PAS_highMej ---
  variableName = "h1_MissingHits1stEle_PAS_highMej"

  plot608 = Plot()
  ## inputs for stacked histograms
  plot608.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot608.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot608.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot608.keys            = keys
  plot608.xtit            = "Missing Hits 1st electron - highMej"
  plot608.ytit            = "Number of events"
  plot608.ylog            = "no"
  plot608.rebin           = 1
  #plot608.xmin            = -0.01
  #plot608.xmax            = 0.01
  plot608.ymin            = 0
  plot608.ymax            = 80
  #plot608.lpos = "bottom-center"
  plot608.name            = "MissingHits1stEle_presel_highMej"
  plot608.addZUncBand     = zUncBand
  plot608.makeRatio       = makeRatio
  plot608.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Dist1stEle_PAS_highMej ---
  variableName = "h1_Dist1stEle_PAS_highMej"

  plot609 = Plot()
  ## inputs for stacked histograms
  plot609.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot609.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot609.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot609.keys            = keys
  plot609.xtit            = "|Dist| 1st electron - highMej"
  plot609.ytit            = "Number of events"
  plot609.ylog            = "no"
  plot609.rebin           = 2
  plot609.xmin            = 0.
  plot609.xmax            = 0.6
  plot609.ymin            = 0
  plot609.ymax            = 25
  #plot609.lpos = "bottom-center"
  plot609.name            = "Dist1stEle_presel_highMej"
  plot609.addZUncBand     = zUncBand
  plot609.makeRatio       = makeRatio
  plot609.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DCotTheta1stEle_PAS_highMej ---
  variableName = "h1_DCotTheta1stEle_PAS_highMej"

  plot610 = Plot()
  ## inputs for stacked histograms
  plot610.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot610.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot610.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot610.keys            = keys
  plot610.xtit            = "|#DeltaCot(#theta)| 1st electron - highMej"
  plot610.ytit            = "Number of events"
  plot610.ylog            = "no"
  plot610.rebin           = 2
  plot610.xmin            = 0.
  plot610.xmax            = 0.8
  plot610.ymin            = 0
  plot610.ymax            = 25
  #plot610.lpos = "bottom-center"
  plot610.name            = "DCotTheta1stEle_presel_highMej"
  plot610.addZUncBand     = zUncBand
  plot610.makeRatio       = makeRatio
  plot610.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Conversion1stEle_highMej ---
  variableName = "h1_Conversion1stEle_highMej"

  plot611 = Plot()
  ## inputs for stacked histograms
  plot611.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot611.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot611.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot611.keys            = keys
  plot611.xtit            = "Conversion flag 1st electron - highMej"
  plot611.ytit            = "Number of events"
  plot611.ylog            = "no"
  #plot611.rebin           = 1
  #plot611.xmin            = -0.5
  #plot611.xmax            = 1.5
  plot611.ymin            = 0
  plot611.ymax            = 90
  #plot611.lpos = "bottom-center"
  plot611.name            = "Conversion1stEle_presel_highMej"
  plot611.addZUncBand     = zUncBand
  plot611.makeRatio       = makeRatio
  plot611.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

  # ## Additional plots to investigate excess around DeltaPhiMET1stJet=3

  # mDeltaPhiMET1stJet>2.5

  #--- h1_Pt1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot700 = Plot()
  ## inputs for stacked histograms
  plot700.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot700.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot700.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot700.keys            = keys
  plot700.xtit            = "p_{T} 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot700.ytit            = "Number of events"
  plot700.ylog            = "no"
  plot700.rebin           = "var"
  #plot700.xmin            = 0
  #plot700.xmax            = 1000
  plot700.ymin            = 0
  plot700.ymax            = 15
  #plot700.lpos = "bottom-center"
  plot700.name            = "Pt1stJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot700.addZUncBand     = zUncBand
  plot700.makeRatio       = makeRatio
  plot700.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot700.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot701 = Plot()
  ## inputs for stacked histograms
  plot701.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot701.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot701.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot701.keys            = keys
  plot701.xtit            = "#eta 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot701.ytit            = "Number of events"
  plot701.ylog            = "no"
  plot701.rebin           = 5
  #plot701.xmin            = -5
  #plot701.xmax            = 5
  plot701.ymin            = 0
  plot701.ymax            = 15
  #plot701.lpos = "bottom-center"
  plot701.name            = "Eta1stJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot701.addZUncBand     = zUncBand
  plot701.makeRatio       = makeRatio
  plot701.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot702 = Plot()
  ## inputs for stacked histograms
  plot702.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot702.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot702.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot702.keys            = keys
  plot702.xtit            = "#phi 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot702.ytit            = "Number of events"
  plot702.ylog            = "no"
  plot702.rebin           = 10
  #plot702.xmin            = -3.1416
  #plot702.xmax            = 3.1416
  plot702.ymin            = 0
  plot702.ymax            = 20
  #plot702.lpos = "bottom-center"
  plot702.name            = "Phi1stJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot702.addZUncBand     = zUncBand
  plot702.makeRatio       = makeRatio
  plot702.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot703 = Plot()
  ## inputs for stacked histograms
  plot703.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot703.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot703.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot703.keys            = keys
  plot703.xtit            = "p_{T} 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot703.ytit            = "Number of events"
  plot703.ylog            = "no"
  plot703.rebin           = "var"
  #plot703.xmin            = 0
  #plot703.xmax            = 1000
  plot703.ymin            = 0
  plot703.ymax            = 20
  #plot703.lpos = "bottom-center"
  plot703.name            = "Pt2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot703.addZUncBand     = zUncBand
  plot703.makeRatio       = makeRatio
  plot703.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot703.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot704 = Plot()
  ## inputs for stacked histograms
  plot704.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot704.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot704.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot704.keys            = keys
  plot704.xtit            = "#eta 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot704.ytit            = "Number of events"
  plot704.ylog            = "no"
  plot704.rebin           = 5
  #plot704.xmin            = -5
  #plot704.xmax            = 5
  plot704.ymin            = 0
  plot704.ymax            = 15
  #plot704.lpos = "bottom-center"
  plot704.name            = "Eta2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot704.addZUncBand     = zUncBand
  plot704.makeRatio       = makeRatio
  plot704.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot705 = Plot()
  ## inputs for stacked histograms
  plot705.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot705.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot705.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot705.keys            = keys
  plot705.xtit            = "#phi 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot705.ytit            = "Number of events"
  plot705.ylog            = "no"
  plot705.rebin           = 10
  #plot705.xmin            = -3.1416
  #plot705.xmax            = 3.1416
  plot705.ymin            = 0
  plot705.ymax            = 20
  #plot705.lpos = "bottom-center"
  plot705.name            = "Phi2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot705.addZUncBand     = zUncBand
  plot705.makeRatio       = makeRatio
  plot705.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_E1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_E1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot706 = Plot()
  ## inputs for stacked histograms
  plot706.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot706.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot706.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot706.keys            = keys
  plot706.xtit            = "Energy 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot706.ytit            = "Number of events"
  plot706.ylog            = "no"
  plot706.rebin           = "var"
  #plot706.xmin            = 0
  #plot706.xmax            = 1000
  plot706.ymin            = 0
  plot706.ymax            = 20
  #plot706.lpos = "bottom-center"
  plot706.name            = "E1stEle_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot706.addZUncBand     = zUncBand
  plot706.makeRatio       = makeRatio
  plot706.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot706.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot707 = Plot()
  ## inputs for stacked histograms
  plot707.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot707.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot707.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot707.keys            = keys
  plot707.xtit            = "p_{T} 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot707.ytit            = "Number of events"
  plot707.ylog            = "no"
  plot707.rebin           = "var"
  #plot707.xmin            = 0
  #plot707.xmax            = 1000
  plot707.ymin            = 0
  plot707.ymax            = 20
  #plot707.lpos = "bottom-center"
  plot707.name            = "Pt1stEle_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot707.addZUncBand     = zUncBand
  plot707.makeRatio       = makeRatio
  plot707.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot707.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot708 = Plot()
  ## inputs for stacked histograms
  plot708.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot708.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot708.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot708.keys            = keys
  plot708.xtit            = "#eta 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot708.ytit            = "Number of events"
  plot708.ylog            = "no"
  plot708.rebin           = 5
  #plot708.xmin            = -5
  #plot708.xmax            = 5
  plot708.ymin            = 0
  plot708.ymax            = 15
  #plot708.lpos = "bottom-center"
  plot708.name            = "Eta1stEle_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot708.addZUncBand     = zUncBand
  plot708.makeRatio       = makeRatio
  plot708.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot709 = Plot()
  ## inputs for stacked histograms
  plot709.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot709.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot709.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot709.keys            = keys
  plot709.xtit            = "#phi 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot709.ytit            = "Number of events"
  plot709.ylog            = "no"
  plot709.rebin           = 10
  #plot709.xmin            = -3.1416
  #plot709.xmax            = 3.1416
  plot709.ymin            = 0
  plot709.ymax            = 20
  #plot709.lpos = "bottom-center"
  plot709.name            = "Phi1stEle_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot709.addZUncBand     = zUncBand
  plot709.makeRatio       = makeRatio
  plot709.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Charge1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Charge1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot710 = Plot()
  ## inputs for stacked histograms
  plot710.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot710.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot710.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot710.keys            = keys
  plot710.xtit            = "Charge 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot710.ytit            = "Number of events"
  plot710.ylog            = "no"
  #plot710.rebin           = 1
  #plot710.xmin            = -1.001
  #plot710.xmax            = 1.001
  plot710.ymin            = 0
  plot710.ymax            = 50
  #plot710.lpos = "bottom-center"
  plot710.name            = "Charge1stEle_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot710.addZUncBand     = zUncBand
  plot710.makeRatio       = makeRatio
  plot710.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Conversion1stEle_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Conversion1stEle_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot711 = Plot()
  ## inputs for stacked histograms
  plot711.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot711.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot711.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot711.keys            = keys
  plot711.xtit            = "Conversion flag 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot711.ytit            = "Number of events"
  plot711.ylog            = "no"
  #plot711.rebin           = 1
  #plot711.xmin            = -0.5
  #plot711.xmax            = 1.5
  plot711.ymin            = 0
  plot711.ymax            = 60
  #plot711.lpos = "bottom-center"
  plot711.name            = "Conversion1stEle_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot711.addZUncBand     = zUncBand
  plot711.makeRatio       = makeRatio
  plot711.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMETEle_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiMETEle_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot712 = Plot()
  ## inputs for stacked histograms
  plot712.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot712.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot712.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot712.keys            = keys
  plot712.xtit            = "#Delta#phi(MET,1^{st} ele) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot712.ytit            = "Number of events"
  plot712.ylog            = "no"
  plot712.rebin           = 10
  #plot712.xmin            = 0
  #plot712.xmax            = 3.1416
  plot712.ymin            = 0
  plot712.ymax            = 20
  #plot712.lpos = "bottom-center"
  plot712.name            = "mDeltaPhiMETEle_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot712.addZUncBand     = zUncBand
  plot712.makeRatio       = makeRatio
  plot712.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiEle2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot713 = Plot()
  ## inputs for stacked histograms
  plot713.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot713.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot713.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot713.keys            = keys
  plot713.xtit            = "#Delta#phi(1^{st} ele,2^{nd} jet) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot713.ytit            = "Number of events"
  plot713.ylog            = "no"
  plot713.rebin           = 10
  #plot713.xmin            = 0
  #plot713.xmax            = 3.1416
  plot713.ymin            = 0
  plot713.ymax            = 20
  #plot713.lpos = "bottom-center"
  plot713.name            = "mDeltaPhiEle2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot713.addZUncBand     = zUncBand
  plot713.makeRatio       = makeRatio
  plot713.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiMET2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot714 = Plot()
  ## inputs for stacked histograms
  plot714.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot714.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot714.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot714.keys            = keys
  plot714.xtit            = "#Delta#phi(MET,2^{nd} jet) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot714.ytit            = "Number of events"
  plot714.ylog            = "no"
  plot714.rebin           = 10
  #plot714.xmin            = 0
  #plot714.xmax            = 3.1416
  plot714.ymin            = 0
  plot714.ymax            = 20
  #plot714.lpos = "bottom-center"
  plot714.name            = "mDeltaPhiMET2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot714.addZUncBand     = zUncBand
  plot714.makeRatio       = makeRatio
  plot714.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Ptenu_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Ptenu_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot715 = Plot()
  ## inputs for stacked histograms
  plot715.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot715.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot715.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot715.keys            = keys
  plot715.xtit            = "p_{T}(e#nu) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot715.ytit            = "Number of events"
  plot715.ylog            = "no"
  plot715.rebin           = "var"
  #plot715.xmin            = 0
  #plot715.xmax            = 2000
  plot715.ymin            = 0
  plot715.ymax            = 20
  #plot715.lpos = "bottom-center"
  plot715.name            = "Ptenu_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot715.addZUncBand     = zUncBand
  plot715.makeRatio       = makeRatio
  plot715.xbins           = [0,40,80,120,160,200,300,600]
  plot715.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_1stJet_PTOverPTPlusMET_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_1stJet_PTOverPTPlusMET_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot716 = Plot()
  ## inputs for stacked histograms
  plot716.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot716.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot716.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot716.keys            = keys
  plot716.xtit            = "Pt1stJet/(Pt1stJet+MET) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot716.ytit            = "Number of events"
  plot716.ylog            = "no"
  plot716.rebin           = 5
  #plot716.xmin            = 0
  #plot716.xmax            = 1
  plot716.ymin            = 0
  plot716.ymax            = 20
  #plot716.lpos = "bottom-center"
  plot716.name            = "1stJet_PTOverPTPlusMET_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot716.addZUncBand     = zUncBand
  plot716.makeRatio       = makeRatio
  plot716.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MET_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_MET_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot717 = Plot()
  ## inputs for stacked histograms
  plot717.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot717.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot717.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot717.keys            = keys
  plot717.xtit            = "pfMET - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot717.ytit            = "Number of events"
  plot717.ylog            = "no"
  plot717.rebin           = "var"
  #plot717.xmin            = 0
  #plot717.xmax            = 1000
  plot717.ymin            = 0
  plot717.ymax            = 20
  #plot717.lpos = "bottom-center"
  plot717.name            = "MET_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot717.addZUncBand     = zUncBand
  plot717.makeRatio       = makeRatio
  plot717.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot717.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Njet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot718 = Plot()
  ## inputs for stacked histograms
  plot718.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot718.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot718.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot718.keys            = keys
  plot718.xtit            = "Number of jets - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot718.ytit            = "Number of events"
  plot718.ylog            = "no"
  plot718.rebin           = 1
  plot718.xmin            = -0.5
  plot718.xmax            = 11.5
  plot718.ymin            = 0
  plot718.ymax            = 25
  #plot718.lpos = "bottom-center"
  plot718.name            = "Njet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot718.addZUncBand     = zUncBand
  plot718.makeRatio       = makeRatio
  plot718.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NjetTCHELBTag_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot719 = Plot()
  ## inputs for stacked histograms
  plot719.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot719.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot719.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot719.keys            = keys
  plot719.xtit            = "Number of b-tagged jets (TCHEL) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot719.ytit            = "Number of events"
  plot719.ylog            = "no"
  plot719.rebin           = 1
  plot719.xmin            = -0.5
  plot719.xmax            = 11.5
  plot719.ymin            = 0
  plot719.ymax            = 25
  #plot719.lpos = "bottom-center"
  plot719.name            = "NjetTCHELBTag_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot719.addZUncBand     = zUncBand
  plot719.makeRatio       = makeRatio
  plot719.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  # mDeltaPhiMET1stJet>2.5, e+ only

  #--- h1_Pt1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1000 = Plot()
  ## inputs for stacked histograms
  plot1000.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1000.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1000.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1000.keys            = keys
  plot1000.xtit            = "p_{T} 1^{st} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1000.ytit            = "Number of events"
  plot1000.ylog            = "no"
  plot1000.rebin           = "var"
  #plot1000.xmin            = 0
  #plot1000.xmax            = 1000
  plot1000.ymin            = 0
  plot1000.ymax            = 15
  #plot1000.lpos = "bottom-center"
  plot1000.name            = "Pt1stJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1000.addZUncBand     = zUncBand
  plot1000.makeRatio       = makeRatio
  plot1000.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot1000.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1001 = Plot()
  ## inputs for stacked histograms
  plot1001.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1001.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1001.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1001.keys            = keys
  plot1001.xtit            = "#eta 1^{st} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1001.ytit            = "Number of events"
  plot1001.ylog            = "no"
  plot1001.rebin           = 5
  #plot1001.xmin            = -5
  #plot1001.xmax            = 5
  plot1001.ymin            = 0
  plot1001.ymax            = 15
  #plot1001.lpos = "bottom-center"
  plot1001.name            = "Eta1stJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1001.addZUncBand     = zUncBand
  plot1001.makeRatio       = makeRatio
  plot1001.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1002 = Plot()
  ## inputs for stacked histograms
  plot1002.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1002.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1002.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1002.keys            = keys
  plot1002.xtit            = "#phi 1^{st} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1002.ytit            = "Number of events"
  plot1002.ylog            = "no"
  plot1002.rebin           = 10
  #plot1002.xmin            = -3.1416
  #plot1002.xmax            = 3.1416
  plot1002.ymin            = 0
  plot1002.ymax            = 20
  #plot1002.lpos = "bottom-center"
  plot1002.name            = "Phi1stJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1002.addZUncBand     = zUncBand
  plot1002.makeRatio       = makeRatio
  plot1002.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1003 = Plot()
  ## inputs for stacked histograms
  plot1003.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1003.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1003.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1003.keys            = keys
  plot1003.xtit            = "p_{T} 2^{nd} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1003.ytit            = "Number of events"
  plot1003.ylog            = "no"
  plot1003.rebin           = "var"
  #plot1003.xmin            = 0
  #plot1003.xmax            = 1000
  plot1003.ymin            = 0
  plot1003.ymax            = 20
  #plot1003.lpos = "bottom-center"
  plot1003.name            = "Pt2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1003.addZUncBand     = zUncBand
  plot1003.makeRatio       = makeRatio
  plot1003.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot1003.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1004 = Plot()
  ## inputs for stacked histograms
  plot1004.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1004.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1004.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1004.keys            = keys
  plot1004.xtit            = "#eta 2^{nd} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1004.ytit            = "Number of events"
  plot1004.ylog            = "no"
  plot1004.rebin           = 5
  #plot1004.xmin            = -5
  #plot1004.xmax            = 5
  plot1004.ymin            = 0
  plot1004.ymax            = 15
  #plot1004.lpos = "bottom-center"
  plot1004.name            = "Eta2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1004.addZUncBand     = zUncBand
  plot1004.makeRatio       = makeRatio
  plot1004.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1005 = Plot()
  ## inputs for stacked histograms
  plot1005.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1005.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1005.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1005.keys            = keys
  plot1005.xtit            = "#phi 2^{nd} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1005.ytit            = "Number of events"
  plot1005.ylog            = "no"
  plot1005.rebin           = 10
  #plot1005.xmin            = -3.1416
  #plot1005.xmax            = 3.1416
  plot1005.ymin            = 0
  plot1005.ymax            = 20
  #plot1005.lpos = "bottom-center"
  plot1005.name            = "Phi2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1005.addZUncBand     = zUncBand
  plot1005.makeRatio       = makeRatio
  plot1005.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_E1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_E1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1006 = Plot()
  ## inputs for stacked histograms
  plot1006.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1006.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1006.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1006.keys            = keys
  plot1006.xtit            = "Energy 1^{st} electron - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1006.ytit            = "Number of events"
  plot1006.ylog            = "no"
  plot1006.rebin           = "var"
  #plot1006.xmin            = 0
  #plot1006.xmax            = 1000
  plot1006.ymin            = 0
  plot1006.ymax            = 20
  #plot1006.lpos = "bottom-center"
  plot1006.name            = "E1stEle_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1006.addZUncBand     = zUncBand
  plot1006.makeRatio       = makeRatio
  plot1006.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot1006.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1007 = Plot()
  ## inputs for stacked histograms
  plot1007.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1007.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1007.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1007.keys            = keys
  plot1007.xtit            = "p_{T} 1^{st} electron - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1007.ytit            = "Number of events"
  plot1007.ylog            = "no"
  plot1007.rebin           = "var"
  #plot1007.xmin            = 0
  #plot1007.xmax            = 1000
  plot1007.ymin            = 0
  plot1007.ymax            = 20
  #plot1007.lpos = "bottom-center"
  plot1007.name            = "Pt1stEle_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1007.addZUncBand     = zUncBand
  plot1007.makeRatio       = makeRatio
  plot1007.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot1007.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1008 = Plot()
  ## inputs for stacked histograms
  plot1008.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1008.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1008.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1008.keys            = keys
  plot1008.xtit            = "#eta 1^{st} electron - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1008.ytit            = "Number of events"
  plot1008.ylog            = "no"
  plot1008.rebin           = 5
  #plot1008.xmin            = -5
  #plot1008.xmax            = 5
  plot1008.ymin            = 0
  plot1008.ymax            = 15
  #plot1008.lpos = "bottom-center"
  plot1008.name            = "Eta1stEle_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1008.addZUncBand     = zUncBand
  plot1008.makeRatio       = makeRatio
  plot1008.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1009 = Plot()
  ## inputs for stacked histograms
  plot1009.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1009.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1009.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1009.keys            = keys
  plot1009.xtit            = "#phi 1^{st} electron - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1009.ytit            = "Number of events"
  plot1009.ylog            = "no"
  plot1009.rebin           = 10
  #plot1009.xmin            = -3.1416
  #plot1009.xmax            = 3.1416
  plot1009.ymin            = 0
  plot1009.ymax            = 20
  #plot1009.lpos = "bottom-center"
  plot1009.name            = "Phi1stEle_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1009.addZUncBand     = zUncBand
  plot1009.makeRatio       = makeRatio
  plot1009.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Charge1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Charge1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1010 = Plot()
  ## inputs for stacked histograms
  plot1010.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1010.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1010.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1010.keys            = keys
  plot1010.xtit            = "Charge 1^{st} electron - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1010.ytit            = "Number of events"
  plot1010.ylog            = "no"
  #plot1010.rebin           = 1
  #plot1010.xmin            = -1.001
  #plot1010.xmax            = 1.001
  plot1010.ymin            = 0
  plot1010.ymax            = 50
  #plot1010.lpos = "bottom-center"
  plot1010.name            = "Charge1stEle_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1010.addZUncBand     = zUncBand
  plot1010.makeRatio       = makeRatio
  plot1010.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Conversion1stEle_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Conversion1stEle_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1011 = Plot()
  ## inputs for stacked histograms
  plot1011.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot1011.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1011.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1011.keys            = keys
  plot1011.xtit            = "Conversion flag 1^{st} electron - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1011.ytit            = "Number of events"
  plot1011.ylog            = "no"
  #plot1011.rebin           = 1
  #plot1011.xmin            = -0.5
  #plot1011.xmax            = 1.5
  plot1011.ymin            = 0
  plot1011.ymax            = 60
  #plot1011.lpos = "bottom-center"
  plot1011.name            = "Conversion1stEle_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1011.addZUncBand     = zUncBand
  plot1011.makeRatio       = makeRatio
  plot1011.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMETEle_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiMETEle_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1012 = Plot()
  ## inputs for stacked histograms
  plot1012.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1012.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1012.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1012.keys            = keys
  plot1012.xtit            = "#Delta#phi(MET,1^{st} ele) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1012.ytit            = "Number of events"
  plot1012.ylog            = "no"
  plot1012.rebin           = 10
  #plot1012.xmin            = 0
  #plot1012.xmax            = 3.1416
  plot1012.ymin            = 0
  plot1012.ymax            = 20
  #plot1012.lpos = "bottom-center"
  plot1012.name            = "mDeltaPhiMETEle_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1012.addZUncBand     = zUncBand
  plot1012.makeRatio       = makeRatio
  plot1012.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiEle2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1013 = Plot()
  ## inputs for stacked histograms
  plot1013.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1013.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1013.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1013.keys            = keys
  plot1013.xtit            = "#Delta#phi(1^{st} ele,2^{nd} jet) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1013.ytit            = "Number of events"
  plot1013.ylog            = "no"
  plot1013.rebin           = 10
  #plot1013.xmin            = 0
  #plot1013.xmax            = 3.1416
  plot1013.ymin            = 0
  plot1013.ymax            = 20
  #plot1013.lpos = "bottom-center"
  plot1013.name            = "mDeltaPhiEle2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1013.addZUncBand     = zUncBand
  plot1013.makeRatio       = makeRatio
  plot1013.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiMET2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1014 = Plot()
  ## inputs for stacked histograms
  plot1014.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1014.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1014.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1014.keys            = keys
  plot1014.xtit            = "#Delta#phi(MET,2^{nd} jet) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1014.ytit            = "Number of events"
  plot1014.ylog            = "no"
  plot1014.rebin           = 10
  #plot1014.xmin            = 0
  #plot1014.xmax            = 3.1416
  plot1014.ymin            = 0
  plot1014.ymax            = 20
  #plot1014.lpos = "bottom-center"
  plot1014.name            = "mDeltaPhiMET2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1014.addZUncBand     = zUncBand
  plot1014.makeRatio       = makeRatio
  plot1014.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Ptenu_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Ptenu_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1015 = Plot()
  ## inputs for stacked histograms
  plot1015.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1015.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1015.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1015.keys            = keys
  plot1015.xtit            = "p_{T}(e#nu) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1015.ytit            = "Number of events"
  plot1015.ylog            = "no"
  plot1015.rebin           = "var"
  #plot1015.xmin            = 0
  #plot1015.xmax            = 2000
  plot1015.ymin            = 0
  plot1015.ymax            = 20
  #plot1015.lpos = "bottom-center"
  plot1015.name            = "Ptenu_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1015.addZUncBand     = zUncBand
  plot1015.makeRatio       = makeRatio
  plot1015.xbins           = [0,40,80,120,160,200,300,600]
  plot1015.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_1stJet_PTOverPTPlusMET_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_1stJet_PTOverPTPlusMET_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1016 = Plot()
  ## inputs for stacked histograms
  plot1016.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1016.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1016.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1016.keys            = keys
  plot1016.xtit            = "Pt1stJet/(Pt1stJet+MET) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1016.ytit            = "Number of events"
  plot1016.ylog            = "no"
  plot1016.rebin           = 5
  #plot1016.xmin            = 0
  #plot1016.xmax            = 1
  plot1016.ymin            = 0
  plot1016.ymax            = 20
  #plot1016.lpos = "bottom-center"
  plot1016.name            = "1stJet_PTOverPTPlusMET_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1016.addZUncBand     = zUncBand
  plot1016.makeRatio       = makeRatio
  plot1016.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MET_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_MET_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1017 = Plot()
  ## inputs for stacked histograms
  plot1017.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1017.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1017.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1017.keys            = keys
  plot1017.xtit            = "pfMET - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1017.ytit            = "Number of events"
  plot1017.ylog            = "no"
  plot1017.rebin           = "var"
  #plot1017.xmin            = 0
  #plot1017.xmax            = 1000
  plot1017.ymin            = 0
  plot1017.ymax            = 20
  #plot1017.lpos = "bottom-center"
  plot1017.name            = "MET_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1017.addZUncBand     = zUncBand
  plot1017.makeRatio       = makeRatio
  plot1017.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot1017.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Njet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1018 = Plot()
  ## inputs for stacked histograms
  plot1018.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1018.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1018.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1018.keys            = keys
  plot1018.xtit            = "Number of jets - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1018.ytit            = "Number of events"
  plot1018.ylog            = "no"
  plot1018.rebin           = 1
  plot1018.xmin            = -0.5
  plot1018.xmax            = 11.5
  plot1018.ymin            = 0
  plot1018.ymax            = 25
  #plot1018.lpos = "bottom-center"
  plot1018.name            = "Njet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1018.addZUncBand     = zUncBand
  plot1018.makeRatio       = makeRatio
  plot1018.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NjetTCHELBTag_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1019 = Plot()
  ## inputs for stacked histograms
  plot1019.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1019.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1019.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1019.keys            = keys
  plot1019.xtit            = "Number of b-tagged jets (TCHEL) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1019.ytit            = "Number of events"
  plot1019.ylog            = "no"
  plot1019.rebin           = 1
  plot1019.xmin            = -0.5
  plot1019.xmax            = 11.5
  plot1019.ymin            = 0
  plot1019.ymax            = 25
  #plot1019.lpos = "bottom-center"
  plot1019.name            = "NjetTCHELBTag_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1019.addZUncBand     = zUncBand
  plot1019.makeRatio       = makeRatio
  plot1019.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  # mDeltaPhiMET1stJet>2.5, e- only

  #--- h1_Pt1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1100 = Plot()
  ## inputs for stacked histograms
  plot1100.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1100.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1100.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1100.keys            = keys
  plot1100.xtit            = "p_{T} 1^{st} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1100.ytit            = "Number of events"
  plot1100.ylog            = "no"
  plot1100.rebin           = "var"
  #plot1100.xmin            = 0
  #plot1100.xmax            = 1000
  plot1100.ymin            = 0
  plot1100.ymax            = 15
  #plot1100.lpos = "bottom-center"
  plot1100.name            = "Pt1stJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1100.addZUncBand     = zUncBand
  plot1100.makeRatio       = makeRatio
  plot1100.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot1100.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1101 = Plot()
  ## inputs for stacked histograms
  plot1101.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1101.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1101.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1101.keys            = keys
  plot1101.xtit            = "#eta 1^{st} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1101.ytit            = "Number of events"
  plot1101.ylog            = "no"
  plot1101.rebin           = 5
  #plot1101.xmin            = -5
  #plot1101.xmax            = 5
  plot1101.ymin            = 0
  plot1101.ymax            = 15
  #plot1101.lpos = "bottom-center"
  plot1101.name            = "Eta1stJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1101.addZUncBand     = zUncBand
  plot1101.makeRatio       = makeRatio
  plot1101.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1102 = Plot()
  ## inputs for stacked histograms
  plot1102.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1102.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1102.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1102.keys            = keys
  plot1102.xtit            = "#phi 1^{st} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1102.ytit            = "Number of events"
  plot1102.ylog            = "no"
  plot1102.rebin           = 10
  #plot1102.xmin            = -3.1416
  #plot1102.xmax            = 3.1416
  plot1102.ymin            = 0
  plot1102.ymax            = 20
  #plot1102.lpos = "bottom-center"
  plot1102.name            = "Phi1stJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1102.addZUncBand     = zUncBand
  plot1102.makeRatio       = makeRatio
  plot1102.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1103 = Plot()
  ## inputs for stacked histograms
  plot1103.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1103.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1103.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1103.keys            = keys
  plot1103.xtit            = "p_{T} 2^{nd} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1103.ytit            = "Number of events"
  plot1103.ylog            = "no"
  plot1103.rebin           = "var"
  #plot1103.xmin            = 0
  #plot1103.xmax            = 1000
  plot1103.ymin            = 0
  plot1103.ymax            = 20
  #plot1103.lpos = "bottom-center"
  plot1103.name            = "Pt2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1103.addZUncBand     = zUncBand
  plot1103.makeRatio       = makeRatio
  plot1103.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot1103.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1104 = Plot()
  ## inputs for stacked histograms
  plot1104.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1104.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1104.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1104.keys            = keys
  plot1104.xtit            = "#eta 2^{nd} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1104.ytit            = "Number of events"
  plot1104.ylog            = "no"
  plot1104.rebin           = 5
  #plot1104.xmin            = -5
  #plot1104.xmax            = 5
  plot1104.ymin            = 0
  plot1104.ymax            = 15
  #plot1104.lpos = "bottom-center"
  plot1104.name            = "Eta2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1104.addZUncBand     = zUncBand
  plot1104.makeRatio       = makeRatio
  plot1104.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1105 = Plot()
  ## inputs for stacked histograms
  plot1105.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1105.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1105.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1105.keys            = keys
  plot1105.xtit            = "#phi 2^{nd} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1105.ytit            = "Number of events"
  plot1105.ylog            = "no"
  plot1105.rebin           = 10
  #plot1105.xmin            = -3.1416
  #plot1105.xmax            = 3.1416
  plot1105.ymin            = 0
  plot1105.ymax            = 20
  #plot1105.lpos = "bottom-center"
  plot1105.name            = "Phi2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1105.addZUncBand     = zUncBand
  plot1105.makeRatio       = makeRatio
  plot1105.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_E1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_E1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1106 = Plot()
  ## inputs for stacked histograms
  plot1106.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1106.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1106.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1106.keys            = keys
  plot1106.xtit            = "Energy 1^{st} electron - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1106.ytit            = "Number of events"
  plot1106.ylog            = "no"
  plot1106.rebin           = "var"
  #plot1106.xmin            = 0
  #plot1106.xmax            = 1000
  plot1106.ymin            = 0
  plot1106.ymax            = 20
  #plot1106.lpos = "bottom-center"
  plot1106.name            = "E1stEle_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1106.addZUncBand     = zUncBand
  plot1106.makeRatio       = makeRatio
  plot1106.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot1106.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1107 = Plot()
  ## inputs for stacked histograms
  plot1107.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1107.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1107.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1107.keys            = keys
  plot1107.xtit            = "p_{T} 1^{st} electron - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1107.ytit            = "Number of events"
  plot1107.ylog            = "no"
  plot1107.rebin           = "var"
  #plot1107.xmin            = 0
  #plot1107.xmax            = 1000
  plot1107.ymin            = 0
  plot1107.ymax            = 20
  #plot1107.lpos = "bottom-center"
  plot1107.name            = "Pt1stEle_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1107.addZUncBand     = zUncBand
  plot1107.makeRatio       = makeRatio
  plot1107.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot1107.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1108 = Plot()
  ## inputs for stacked histograms
  plot1108.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1108.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1108.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1108.keys            = keys
  plot1108.xtit            = "#eta 1^{st} electron - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1108.ytit            = "Number of events"
  plot1108.ylog            = "no"
  plot1108.rebin           = 5
  #plot1108.xmin            = -5
  #plot1108.xmax            = 5
  plot1108.ymin            = 0
  plot1108.ymax            = 15
  #plot1108.lpos = "bottom-center"
  plot1108.name            = "Eta1stEle_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1108.addZUncBand     = zUncBand
  plot1108.makeRatio       = makeRatio
  plot1108.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1109 = Plot()
  ## inputs for stacked histograms
  plot1109.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1109.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1109.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1109.keys            = keys
  plot1109.xtit            = "#phi 1^{st} electron - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1109.ytit            = "Number of events"
  plot1109.ylog            = "no"
  plot1109.rebin           = 10
  #plot1109.xmin            = -3.1416
  #plot1109.xmax            = 3.1416
  plot1109.ymin            = 0
  plot1109.ymax            = 20
  #plot1109.lpos = "bottom-center"
  plot1109.name            = "Phi1stEle_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1109.addZUncBand     = zUncBand
  plot1109.makeRatio       = makeRatio
  plot1109.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Charge1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Charge1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1110 = Plot()
  ## inputs for stacked histograms
  plot1110.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1110.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1110.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1110.keys            = keys
  plot1110.xtit            = "Charge 1^{st} electron - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1110.ytit            = "Number of events"
  plot1110.ylog            = "no"
  #plot1110.rebin           = 1
  #plot1110.xmin            = -1.001
  #plot1110.xmax            = 1.001
  plot1110.ymin            = 0
  plot1110.ymax            = 50
  #plot1110.lpos = "bottom-center"
  plot1110.name            = "Charge1stEle_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1110.addZUncBand     = zUncBand
  plot1110.makeRatio       = makeRatio
  plot1110.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Conversion1stEle_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Conversion1stEle_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1111 = Plot()
  ## inputs for stacked histograms
  plot1111.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot1111.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1111.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1111.keys            = keys
  plot1111.xtit            = "Conversion flag 1^{st} electron - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1111.ytit            = "Number of events"
  plot1111.ylog            = "no"
  #plot1111.rebin           = 1
  #plot1111.xmin            = -0.5
  #plot1111.xmax            = 1.5
  plot1111.ymin            = 0
  plot1111.ymax            = 60
  #plot1111.lpos = "bottom-center"
  plot1111.name            = "Conversion1stEle_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1111.addZUncBand     = zUncBand
  plot1111.makeRatio       = makeRatio
  plot1111.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMETEle_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiMETEle_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1112 = Plot()
  ## inputs for stacked histograms
  plot1112.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1112.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1112.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1112.keys            = keys
  plot1112.xtit            = "#Delta#phi(MET,1^{st} ele) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1112.ytit            = "Number of events"
  plot1112.ylog            = "no"
  plot1112.rebin           = 10
  #plot1112.xmin            = 0
  #plot1112.xmax            = 3.1416
  plot1112.ymin            = 0
  plot1112.ymax            = 20
  #plot1112.lpos = "bottom-center"
  plot1112.name            = "mDeltaPhiMETEle_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1112.addZUncBand     = zUncBand
  plot1112.makeRatio       = makeRatio
  plot1112.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiEle2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1113 = Plot()
  ## inputs for stacked histograms
  plot1113.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1113.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1113.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1113.keys            = keys
  plot1113.xtit            = "#Delta#phi(1^{st} ele,2^{nd} jet) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1113.ytit            = "Number of events"
  plot1113.ylog            = "no"
  plot1113.rebin           = 10
  #plot1113.xmin            = 0
  #plot1113.xmax            = 3.1416
  plot1113.ymin            = 0
  plot1113.ymax            = 20
  #plot1113.lpos = "bottom-center"
  plot1113.name            = "mDeltaPhiEle2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1113.addZUncBand     = zUncBand
  plot1113.makeRatio       = makeRatio
  plot1113.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiMET2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1114 = Plot()
  ## inputs for stacked histograms
  plot1114.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1114.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1114.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1114.keys            = keys
  plot1114.xtit            = "#Delta#phi(MET,2^{nd} jet) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1114.ytit            = "Number of events"
  plot1114.ylog            = "no"
  plot1114.rebin           = 10
  #plot1114.xmin            = 0
  #plot1114.xmax            = 3.1416
  plot1114.ymin            = 0
  plot1114.ymax            = 20
  #plot1114.lpos = "bottom-center"
  plot1114.name            = "mDeltaPhiMET2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1114.addZUncBand     = zUncBand
  plot1114.makeRatio       = makeRatio
  plot1114.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Ptenu_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Ptenu_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1115 = Plot()
  ## inputs for stacked histograms
  plot1115.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1115.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1115.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1115.keys            = keys
  plot1115.xtit            = "p_{T}(e#nu) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1115.ytit            = "Number of events"
  plot1115.ylog            = "no"
  plot1115.rebin           = "var"
  #plot1115.xmin            = 0
  #plot1115.xmax            = 2000
  plot1115.ymin            = 0
  plot1115.ymax            = 20
  #plot1115.lpos = "bottom-center"
  plot1115.name            = "Ptenu_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1115.addZUncBand     = zUncBand
  plot1115.makeRatio       = makeRatio
  plot1115.xbins           = [0,40,80,120,160,200,300,600]
  plot1115.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_1stJet_PTOverPTPlusMET_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_1stJet_PTOverPTPlusMET_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1116 = Plot()
  ## inputs for stacked histograms
  plot1116.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1116.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1116.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1116.keys            = keys
  plot1116.xtit            = "Pt1stJet/(Pt1stJet+MET) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1116.ytit            = "Number of events"
  plot1116.ylog            = "no"
  plot1116.rebin           = 5
  #plot1116.xmin            = 0
  #plot1116.xmax            = 1
  plot1116.ymin            = 0
  plot1116.ymax            = 20
  #plot1116.lpos = "bottom-center"
  plot1116.name            = "1stJet_PTOverPTPlusMET_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1116.addZUncBand     = zUncBand
  plot1116.makeRatio       = makeRatio
  plot1116.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MET_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_MET_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1117 = Plot()
  ## inputs for stacked histograms
  plot1117.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1117.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1117.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1117.keys            = keys
  plot1117.xtit            = "pfMET - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1117.ytit            = "Number of events"
  plot1117.ylog            = "no"
  plot1117.rebin           = "var"
  #plot1117.xmin            = 0
  #plot1117.xmax            = 1000
  plot1117.ymin            = 0
  plot1117.ymax            = 20
  #plot1117.lpos = "bottom-center"
  plot1117.name            = "MET_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1117.addZUncBand     = zUncBand
  plot1117.makeRatio       = makeRatio
  plot1117.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot1117.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Njet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1118 = Plot()
  ## inputs for stacked histograms
  plot1118.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1118.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1118.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1118.keys            = keys
  plot1118.xtit            = "Number of jets - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1118.ytit            = "Number of events"
  plot1118.ylog            = "no"
  plot1118.rebin           = 1
  plot1118.xmin            = -0.5
  plot1118.xmax            = 11.5
  plot1118.ymin            = 0
  plot1118.ymax            = 25
  #plot1118.lpos = "bottom-center"
  plot1118.name            = "Njet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1118.addZUncBand     = zUncBand
  plot1118.makeRatio       = makeRatio
  plot1118.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NjetTCHELBTag_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot1119 = Plot()
  ## inputs for stacked histograms
  plot1119.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1119.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1119.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1119.keys            = keys
  plot1119.xtit            = "Number of b-tagged jets (TCHEL) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot1119.ytit            = "Number of events"
  plot1119.ylog            = "no"
  plot1119.rebin           = 1
  plot1119.xmin            = -0.5
  plot1119.xmax            = 11.5
  plot1119.ymin            = 0
  plot1119.ymax            = 25
  #plot1119.lpos = "bottom-center"
  plot1119.name            = "NjetTCHELBTag_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot1119.addZUncBand     = zUncBand
  plot1119.makeRatio       = makeRatio
  plot1119.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  # mDeltaPhiMET1stJet<2.5

  #--- h1_Pt1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Pt1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot800 = Plot()
  ## inputs for stacked histograms
  plot800.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot800.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot800.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot800.keys            = keys
  plot800.xtit            = "p_{T} 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot800.ytit            = "Number of events"
  plot800.ylog            = "no"
  plot800.rebin           = "var"
  #plot800.xmin            = 0
  #plot800.xmax            = 1000
  plot800.ymin            = 0
  plot800.ymax            = 15
  #plot800.lpos = "bottom-center"
  plot800.name            = "Pt1stJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot800.addZUncBand     = zUncBand
  plot800.makeRatio       = makeRatio
  plot800.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot800.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Eta1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot801 = Plot()
  ## inputs for stacked histograms
  plot801.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot801.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot801.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot801.keys            = keys
  plot801.xtit            = "#eta 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot801.ytit            = "Number of events"
  plot801.ylog            = "no"
  plot801.rebin           = 5
  #plot801.xmin            = -5
  #plot801.xmax            = 5
  plot801.ymin            = 0
  plot801.ymax            = 15
  #plot801.lpos = "bottom-center"
  plot801.name            = "Eta1stJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot801.addZUncBand     = zUncBand
  plot801.makeRatio       = makeRatio
  plot801.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Phi1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot802 = Plot()
  ## inputs for stacked histograms
  plot802.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot802.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot802.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot802.keys            = keys
  plot802.xtit            = "#phi 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot802.ytit            = "Number of events"
  plot802.ylog            = "no"
  plot802.rebin           = 10
  #plot802.xmin            = -3.1416
  #plot802.xmax            = 3.1416
  plot802.ymin            = 0
  plot802.ymax            = 20
  #plot802.lpos = "bottom-center"
  plot802.name            = "Phi1stJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot802.addZUncBand     = zUncBand
  plot802.makeRatio       = makeRatio
  plot802.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Pt2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot803 = Plot()
  ## inputs for stacked histograms
  plot803.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot803.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot803.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot803.keys            = keys
  plot803.xtit            = "p_{T} 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot803.ytit            = "Number of events"
  plot803.ylog            = "no"
  plot803.rebin           = "var"
  #plot803.xmin            = 0
  #plot803.xmax            = 1000
  plot803.ymin            = 0
  plot803.ymax            = 20
  #plot803.lpos = "bottom-center"
  plot803.name            = "Pt2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot803.addZUncBand     = zUncBand
  plot803.makeRatio       = makeRatio
  plot803.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot803.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Eta2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot804 = Plot()
  ## inputs for stacked histograms
  plot804.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot804.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot804.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot804.keys            = keys
  plot804.xtit            = "#eta 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot804.ytit            = "Number of events"
  plot804.ylog            = "no"
  plot804.rebin           = 5
  #plot804.xmin            = -5
  #plot804.xmax            = 5
  plot804.ymin            = 0
  plot804.ymax            = 15
  #plot804.lpos = "bottom-center"
  plot804.name            = "Eta2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot804.addZUncBand     = zUncBand
  plot804.makeRatio       = makeRatio
  plot804.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Phi2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot805 = Plot()
  ## inputs for stacked histograms
  plot805.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot805.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot805.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot805.keys            = keys
  plot805.xtit            = "#phi 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot805.ytit            = "Number of events"
  plot805.ylog            = "no"
  plot805.rebin           = 10
  #plot805.xmin            = -3.1416
  #plot805.xmax            = 3.1416
  plot805.ymin            = 0
  plot805.ymax            = 20
  #plot805.lpos = "bottom-center"
  plot805.name            = "Phi2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot805.addZUncBand     = zUncBand
  plot805.makeRatio       = makeRatio
  plot805.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_E1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_E1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot806 = Plot()
  ## inputs for stacked histograms
  plot806.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot806.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot806.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot806.keys            = keys
  plot806.xtit            = "Energy 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot806.ytit            = "Number of events"
  plot806.ylog            = "no"
  plot806.rebin           = "var"
  #plot806.xmin            = 0
  #plot806.xmax            = 1000
  plot806.ymin            = 0
  plot806.ymax            = 20
  #plot806.lpos = "bottom-center"
  plot806.name            = "E1stEle_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot806.addZUncBand     = zUncBand
  plot806.makeRatio       = makeRatio
  plot806.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot806.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Pt1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot807 = Plot()
  ## inputs for stacked histograms
  plot807.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot807.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot807.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot807.keys            = keys
  plot807.xtit            = "p_{T} 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot807.ytit            = "Number of events"
  plot807.ylog            = "no"
  plot807.rebin           = "var"
  #plot807.xmin            = 0
  #plot807.xmax            = 1000
  plot807.ymin            = 0
  plot807.ymax            = 20
  #plot807.lpos = "bottom-center"
  plot807.name            = "Pt1stEle_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot807.addZUncBand     = zUncBand
  plot807.makeRatio       = makeRatio
  plot807.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot807.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Eta1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot808 = Plot()
  ## inputs for stacked histograms
  plot808.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot808.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot808.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot808.keys            = keys
  plot808.xtit            = "#eta 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot808.ytit            = "Number of events"
  plot808.ylog            = "no"
  plot808.rebin           = 5
  #plot808.xmin            = -5
  #plot808.xmax            = 5
  plot808.ymin            = 0
  plot808.ymax            = 15
  #plot808.lpos = "bottom-center"
  plot808.name            = "Eta1stEle_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot808.addZUncBand     = zUncBand
  plot808.makeRatio       = makeRatio
  plot808.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Phi1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot809 = Plot()
  ## inputs for stacked histograms
  plot809.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot809.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot809.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot809.keys            = keys
  plot809.xtit            = "#phi 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot809.ytit            = "Number of events"
  plot809.ylog            = "no"
  plot809.rebin           = 10
  #plot809.xmin            = -3.1416
  #plot809.xmax            = 3.1416
  plot809.ymin            = 0
  plot809.ymax            = 20
  #plot809.lpos = "bottom-center"
  plot809.name            = "Phi1stEle_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot809.addZUncBand     = zUncBand
  plot809.makeRatio       = makeRatio
  plot809.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Charge1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Charge1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot810 = Plot()
  ## inputs for stacked histograms
  plot810.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot810.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot810.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot810.keys            = keys
  plot810.xtit            = "Charge 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot810.ytit            = "Number of events"
  plot810.ylog            = "no"
  #plot810.rebin           = 1
  #plot810.xmin            = -1.001
  #plot810.xmax            = 1.001
  plot810.ymin            = 0
  plot810.ymax            = 50
  #plot810.lpos = "bottom-center"
  plot810.name            = "Charge1stEle_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot810.addZUncBand     = zUncBand
  plot810.makeRatio       = makeRatio
  plot810.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Conversion1stEle_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Conversion1stEle_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot811 = Plot()
  ## inputs for stacked histograms
  plot811.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot811.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot811.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot811.keys            = keys
  plot811.xtit            = "Conversion flag 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot811.ytit            = "Number of events"
  plot811.ylog            = "no"
  #plot811.rebin           = 1
  #plot811.xmin            = -0.5
  #plot811.xmax            = 1.5
  plot811.ymin            = 0
  plot811.ymax            = 60
  #plot811.lpos = "bottom-center"
  plot811.name            = "Conversion1stEle_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot811.addZUncBand     = zUncBand
  plot811.makeRatio       = makeRatio
  plot811.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMETEle_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_mDeltaPhiMETEle_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot812 = Plot()
  ## inputs for stacked histograms
  plot812.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot812.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot812.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot812.keys            = keys
  plot812.xtit            = "#Delta#phi(MET,1^{st} ele) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot812.ytit            = "Number of events"
  plot812.ylog            = "no"
  plot812.rebin           = 10
  #plot812.xmin            = 0
  #plot812.xmax            = 3.1416
  plot812.ymin            = 0
  plot812.ymax            = 20
  #plot812.lpos = "bottom-center"
  plot812.name            = "mDeltaPhiMETEle_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot812.addZUncBand     = zUncBand
  plot812.makeRatio       = makeRatio
  plot812.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_mDeltaPhiEle2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot813 = Plot()
  ## inputs for stacked histograms
  plot813.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot813.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot813.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot813.keys            = keys
  plot813.xtit            = "#Delta#phi(1^{st} ele,2^{nd} jet) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot813.ytit            = "Number of events"
  plot813.ylog            = "no"
  plot813.rebin           = 10
  #plot813.xmin            = 0
  #plot813.xmax            = 3.1416
  plot813.ymin            = 0
  plot813.ymax            = 20
  #plot813.lpos = "bottom-center"
  plot813.name            = "mDeltaPhiEle2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot813.addZUncBand     = zUncBand
  plot813.makeRatio       = makeRatio
  plot813.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_mDeltaPhiMET2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot814 = Plot()
  ## inputs for stacked histograms
  plot814.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot814.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot814.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot814.keys            = keys
  plot814.xtit            = "#Delta#phi(MET,2^{nd} jet) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot814.ytit            = "Number of events"
  plot814.ylog            = "no"
  plot814.rebin           = 10
  #plot814.xmin            = 0
  #plot814.xmax            = 3.1416
  plot814.ymin            = 0
  plot814.ymax            = 20
  #plot814.lpos = "bottom-center"
  plot814.name            = "mDeltaPhiMET2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot814.addZUncBand     = zUncBand
  plot814.makeRatio       = makeRatio
  plot814.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Ptenu_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Ptenu_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot815 = Plot()
  ## inputs for stacked histograms
  plot815.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot815.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot815.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot815.keys            = keys
  plot815.xtit            = "p_{T}(e#nu) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot815.ytit            = "Number of events"
  plot815.ylog            = "no"
  plot815.rebin           = "var"
  #plot815.xmin            = 0
  #plot815.xmax            = 2000
  plot815.ymin            = 0
  plot815.ymax            = 20
  #plot815.lpos = "bottom-center"
  plot815.name            = "Ptenu_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot815.addZUncBand     = zUncBand
  plot815.makeRatio       = makeRatio
  plot815.xbins           = [0,40,80,120,160,200,300,600]
  plot815.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_1stJet_PTOverPTPlusMET_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_1stJet_PTOverPTPlusMET_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot816 = Plot()
  ## inputs for stacked histograms
  plot816.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot816.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot816.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot816.keys            = keys
  plot816.xtit            = "Pt1stJet/(Pt1stJet+MET) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot816.ytit            = "Number of events"
  plot816.ylog            = "no"
  plot816.rebin           = 5
  #plot816.xmin            = 0
  #plot816.xmax            = 1
  plot816.ymin            = 0
  plot816.ymax            = 20
  #plot816.lpos = "bottom-center"
  plot816.name            = "1stJet_PTOverPTPlusMET_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot816.addZUncBand     = zUncBand
  plot816.makeRatio       = makeRatio
  plot816.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MET_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_MET_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot817 = Plot()
  ## inputs for stacked histograms
  plot817.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot817.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot817.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot817.keys            = keys
  plot817.xtit            = "pfMET - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot817.ytit            = "Number of events"
  plot817.ylog            = "no"
  plot817.rebin           = "var"
  #plot817.xmin            = 0
  #plot817.xmax            = 1000
  plot817.ymin            = 0
  plot817.ymax            = 20
  #plot817.lpos = "bottom-center"
  plot817.name            = "MET_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot817.addZUncBand     = zUncBand
  plot817.makeRatio       = makeRatio
  plot817.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot817.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Njet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot818 = Plot()
  ## inputs for stacked histograms
  plot818.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot818.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot818.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot818.keys            = keys
  plot818.xtit            = "Number of jets - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot818.ytit            = "Number of events"
  plot818.ylog            = "no"
  plot818.rebin           = 1
  plot818.xmin            = -0.5
  plot818.xmax            = 11.5
  plot818.ymin            = 0
  plot818.ymax            = 25
  #plot818.lpos = "bottom-center"
  plot818.name            = "Njet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot818.addZUncBand     = zUncBand
  plot818.makeRatio       = makeRatio
  plot818.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_NjetTCHELBTag_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot819 = Plot()
  ## inputs for stacked histograms
  plot819.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot819.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot819.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot819.keys            = keys
  plot819.xtit            = "Number of b-tagged jets (TCHEL) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot819.ytit            = "Number of events"
  plot819.ylog            = "no"
  plot819.rebin           = 1
  plot819.xmin            = -0.5
  plot819.xmax            = 11.5
  plot819.ymin            = 0
  plot819.ymax            = 25
  #plot819.lpos = "bottom-center"
  plot819.name            = "NjetTCHELBTag_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot819.addZUncBand     = zUncBand
  plot819.makeRatio       = makeRatio
  plot819.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  # ## Additional plots to investigate bump in Eta1stJet distribution

  #--- h1_Pt1stJet_PAS_Eta1stJetBump ---
  variableName = "h1_Pt1stJet_PAS_Eta1stJetBump"

  plot900 = Plot()
  ## inputs for stacked histograms
  plot900.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot900.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot900.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot900.keys            = keys
  plot900.xtit            = "p_{T} 1^{st} jet - Eta1stJetBump"
  plot900.ytit            = "Number of events"
  plot900.ylog            = "no"
  plot900.rebin           = "var"
  #plot900.xmin            = 0
  #plot900.xmax            = 1000
  plot900.ymin            = 0
  plot900.ymax            = 20
  #plot900.lpos = "bottom-center"
  plot900.name            = "Pt1stJet_presel_Eta1stJetBump"
  plot900.addZUncBand     = zUncBand
  plot900.makeRatio       = makeRatio
  plot900.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot900.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  ##--- h1_Eta1stJet_PAS_Eta1stJetBump ---
  #variableName = "h1_Eta1stJet_PAS_Eta1stJetBump"

  #plot901 = Plot()
  ### inputs for stacked histograms
  #plot901.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  #plot901.keysStack       = keysStack
  ### this is the list of histograms that should be simply overlaid on top of the stacked histogram
  #plot901.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  #plot901.keys            = keys
  #plot901.xtit            = "#eta 1^{st} jet - Eta1stJetBump"
  #plot901.ytit            = "Number of events"
  #plot901.ylog            = "no"
  #plot901.rebin           = 5
  ##plot901.xmin            = -5
  ##plot901.xmax            = 5
  #plot901.ymin            = 0
  #plot901.ymax            = 15
  ##plot901.lpos = "bottom-center"
  #plot901.name            = "Eta1stJet_presel_Eta1stJetBump"
  #plot901.addZUncBand     = zUncBand
  #plot901.makeRatio       = makeRatio
  #plot901.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stJet_PAS_Eta1stJetBump ---
  variableName = "h1_Phi1stJet_PAS_Eta1stJetBump"

  plot902 = Plot()
  ## inputs for stacked histograms
  plot902.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot902.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot902.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot902.keys            = keys
  plot902.xtit            = "#phi 1^{st} jet - Eta1stJetBump"
  plot902.ytit            = "Number of events"
  plot902.ylog            = "no"
  plot902.rebin           = 10
  #plot902.xmin            = -3.1416
  #plot902.xmax            = 3.1416
  plot902.ymin            = 0
  plot902.ymax            = 20
  #plot902.lpos = "bottom-center"
  plot902.name            = "Phi1stJet_presel_Eta1stJetBump"
  plot902.addZUncBand     = zUncBand
  plot902.makeRatio       = makeRatio
  plot902.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt2ndJet_PAS_Eta1stJetBump ---
  variableName = "h1_Pt2ndJet_PAS_Eta1stJetBump"

  plot903 = Plot()
  ## inputs for stacked histograms
  plot903.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot903.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot903.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot903.keys            = keys
  plot903.xtit            = "p_{T} 2^{nd} jet - Eta1stJetBump"
  plot903.ytit            = "Number of events"
  plot903.ylog            = "no"
  plot903.rebin           = "var"
  #plot903.xmin            = 0
  #plot903.xmax            = 1000
  plot903.ymin            = 0
  plot903.ymax            = 25
  #plot903.lpos = "bottom-center"
  plot903.name            = "Pt2ndJet_presel_Eta1stJetBump"
  plot903.addZUncBand     = zUncBand
  plot903.makeRatio       = makeRatio
  plot903.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot903.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta2ndJet_PAS_Eta1stJetBump ---
  variableName = "h1_Eta2ndJet_PAS_Eta1stJetBump"

  plot904 = Plot()
  ## inputs for stacked histograms
  plot904.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot904.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot904.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot904.keys            = keys
  plot904.xtit            = "#eta 2^{nd} jet - Eta1stJetBump"
  plot904.ytit            = "Number of events"
  plot904.ylog            = "no"
  plot904.rebin           = 5
  #plot904.xmin            = -5
  #plot904.xmax            = 5
  plot904.ymin            = 0
  plot904.ymax            = 15
  #plot904.lpos = "bottom-center"
  plot904.name            = "Eta2ndJet_presel_Eta1stJetBump"
  plot904.addZUncBand     = zUncBand
  plot904.makeRatio       = makeRatio
  plot904.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi2ndJet_PAS_Eta1stJetBump ---
  variableName = "h1_Phi2ndJet_PAS_Eta1stJetBump"

  plot905 = Plot()
  ## inputs for stacked histograms
  plot905.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot905.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot905.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot905.keys            = keys
  plot905.xtit            = "#phi 2^{nd} jet - Eta1stJetBump"
  plot905.ytit            = "Number of events"
  plot905.ylog            = "no"
  plot905.rebin           = 10
  #plot905.xmin            = -3.1416
  #plot905.xmax            = 3.1416
  plot905.ymin            = 0
  plot905.ymax            = 20
  #plot905.lpos = "bottom-center"
  plot905.name            = "Phi2ndJet_presel_Eta1stJetBump"
  plot905.addZUncBand     = zUncBand
  plot905.makeRatio       = makeRatio
  plot905.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_E1stEle_PAS_Eta1stJetBump ---
  variableName = "h1_E1stEle_PAS_Eta1stJetBump"

  plot906 = Plot()
  ## inputs for stacked histograms
  plot906.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot906.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot906.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot906.keys            = keys
  plot906.xtit            = "Energy 1^{st} electron - Eta1stJetBump"
  plot906.ytit            = "Number of events"
  plot906.ylog            = "no"
  plot906.rebin           = "var"
  #plot906.xmin            = 0
  #plot906.xmax            = 1000
  plot906.ymin            = 0
  plot906.ymax            = 20
  #plot906.lpos = "bottom-center"
  plot906.name            = "E1stEle_presel_Eta1stJetBump"
  plot906.addZUncBand     = zUncBand
  plot906.makeRatio       = makeRatio
  plot906.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot906.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt1stEle_PAS_Eta1stJetBump ---
  variableName = "h1_Pt1stEle_PAS_Eta1stJetBump"

  plot907 = Plot()
  ## inputs for stacked histograms
  plot907.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot907.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot907.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot907.keys            = keys
  plot907.xtit            = "p_{T} 1^{st} electron - Eta1stJetBump"
  plot907.ytit            = "Number of events"
  plot907.ylog            = "no"
  plot907.rebin           = "var"
  #plot907.xmin            = 0
  #plot907.xmax            = 1000
  plot907.ymin            = 0
  plot907.ymax            = 25
  #plot907.lpos = "bottom-center"
  plot907.name            = "Pt1stEle_presel_Eta1stJetBump"
  plot907.addZUncBand     = zUncBand
  plot907.makeRatio       = makeRatio
  plot907.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot907.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stEle_PAS_Eta1stJetBump ---
  variableName = "h1_Eta1stEle_PAS_Eta1stJetBump"

  plot908 = Plot()
  ## inputs for stacked histograms
  plot908.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot908.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot908.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot908.keys            = keys
  plot908.xtit            = "#eta 1^{st} electron - Eta1stJetBump"
  plot908.ytit            = "Number of events"
  plot908.ylog            = "no"
  plot908.rebin           = 5
  #plot908.xmin            = -5
  #plot908.xmax            = 5
  plot908.ymin            = 0
  plot908.ymax            = 15
  #plot908.lpos = "bottom-center"
  plot908.name            = "Eta1stEle_presel_Eta1stJetBump"
  plot908.addZUncBand     = zUncBand
  plot908.makeRatio       = makeRatio
  plot908.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stEle_PAS_Eta1stJetBump ---
  variableName = "h1_Phi1stEle_PAS_Eta1stJetBump"

  plot909 = Plot()
  ## inputs for stacked histograms
  plot909.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot909.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot909.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot909.keys            = keys
  plot909.xtit            = "#phi 1^{st} electron - Eta1stJetBump"
  plot909.ytit            = "Number of events"
  plot909.ylog            = "no"
  plot909.rebin           = 10
  #plot909.xmin            = -3.1416
  #plot909.xmax            = 3.1416
  plot909.ymin            = 0
  plot909.ymax            = 20
  #plot909.lpos = "bottom-center"
  plot909.name            = "Phi1stEle_presel_Eta1stJetBump"
  plot909.addZUncBand     = zUncBand
  plot909.makeRatio       = makeRatio
  plot909.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Charge1stEle_PAS_Eta1stJetBump ---
  variableName = "h1_Charge1stEle_PAS_Eta1stJetBump"

  plot910 = Plot()
  ## inputs for stacked histograms
  plot910.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot910.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot910.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot910.keys            = keys
  plot910.xtit            = "Charge 1^{st} electron - Eta1stJetBump"
  plot910.ytit            = "Number of events"
  plot910.ylog            = "no"
  #plot910.rebin           = 1
  #plot910.xmin            = -1.001
  #plot910.xmax            = 1.001
  plot910.ymin            = 0
  plot910.ymax            = 50
  #plot910.lpos = "bottom-center"
  plot910.name            = "Charge1stEle_presel_Eta1stJetBump"
  plot910.addZUncBand     = zUncBand
  plot910.makeRatio       = makeRatio
  plot910.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Conversion1stEle_Eta1stJetBump ---
  variableName = "h1_Conversion1stEle_Eta1stJetBump"

  plot911 = Plot()
  ## inputs for stacked histograms
  plot911.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot911.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot911.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot911.keys            = keys
  plot911.xtit            = "Conversion flag 1^{st} electron - Eta1stJetBump"
  plot911.ytit            = "Number of events"
  plot911.ylog            = "no"
  #plot911.rebin           = 1
  #plot911.xmin            = -0.5
  #plot911.xmax            = 1.5
  plot911.ymin            = 0
  plot911.ymax            = 60
  #plot911.lpos = "bottom-center"
  plot911.name            = "Conversion1stEle_presel_Eta1stJetBump"
  plot911.addZUncBand     = zUncBand
  plot911.makeRatio       = makeRatio
  plot911.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMETEle_Eta1stJetBump ---
  variableName = "h1_mDeltaPhiMETEle_Eta1stJetBump"

  plot912 = Plot()
  ## inputs for stacked histograms
  plot912.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot912.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot912.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot912.keys            = keys
  plot912.xtit            = "#Delta#phi(MET,1^{st} ele) - Eta1stJetBump"
  plot912.ytit            = "Number of events"
  plot912.ylog            = "no"
  plot912.rebin           = 10
  #plot912.xmin            = 0
  #plot912.xmax            = 3.1416
  plot912.ymin            = 0
  plot912.ymax            = 20
  #plot912.lpos = "bottom-center"
  plot912.name            = "mDeltaPhiMETEle_presel_Eta1stJetBump"
  plot912.addZUncBand     = zUncBand
  plot912.makeRatio       = makeRatio
  plot912.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle2ndJet_PAS_Eta1stJetBump ---
  variableName = "h1_mDeltaPhiEle2ndJet_PAS_Eta1stJetBump"

  plot913 = Plot()
  ## inputs for stacked histograms
  plot913.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot913.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot913.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot913.keys            = keys
  plot913.xtit            = "#Delta#phi(1^{st} ele,2^{nd} jet) - Eta1stJetBump"
  plot913.ytit            = "Number of events"
  plot913.ylog            = "no"
  plot913.rebin           = 10
  #plot913.xmin            = 0
  #plot913.xmax            = 3.1416
  plot913.ymin            = 0
  plot913.ymax            = 20
  #plot913.lpos = "bottom-center"
  plot913.name            = "mDeltaPhiEle2ndJet_presel_Eta1stJetBump"
  plot913.addZUncBand     = zUncBand
  plot913.makeRatio       = makeRatio
  plot913.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_Eta1stJetBump ---
  variableName = "h1_mDeltaPhiMET2ndJet_Eta1stJetBump"

  plot914 = Plot()
  ## inputs for stacked histograms
  plot914.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot914.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot914.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot914.keys            = keys
  plot914.xtit            = "#Delta#phi(MET,2^{nd} jet) - Eta1stJetBump"
  plot914.ytit            = "Number of events"
  plot914.ylog            = "no"
  plot914.rebin           = 10
  #plot914.xmin            = 0
  #plot914.xmax            = 3.1416
  plot914.ymin            = 0
  plot914.ymax            = 20
  #plot914.lpos = "bottom-center"
  plot914.name            = "mDeltaPhiMET2ndJet_presel_Eta1stJetBump"
  plot914.addZUncBand     = zUncBand
  plot914.makeRatio       = makeRatio
  plot914.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Ptenu_PAS_Eta1stJetBump ---
  variableName = "h1_Ptenu_PAS_Eta1stJetBump"

  plot915 = Plot()
  ## inputs for stacked histograms
  plot915.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot915.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot915.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot915.keys            = keys
  plot915.xtit            = "p_{T}(e#nu) - Eta1stJetBump"
  plot915.ytit            = "Number of events"
  plot915.ylog            = "no"
  plot915.rebin           = "var"
  #plot915.xmin            = 0
  #plot915.xmax            = 2000
  plot915.ymin            = 0
  plot915.ymax            = 20
  #plot915.lpos = "bottom-center"
  plot915.name            = "Ptenu_presel_Eta1stJetBump"
  plot915.addZUncBand     = zUncBand
  plot915.makeRatio       = makeRatio
  plot915.xbins           = [0,40,80,120,160,200,300,600]
  plot915.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_1stJet_PTOverPTPlusMET_Eta1stJetBump ---
  variableName = "h1_1stJet_PTOverPTPlusMET_Eta1stJetBump"

  plot916 = Plot()
  ## inputs for stacked histograms
  plot916.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot916.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot916.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot916.keys            = keys
  plot916.xtit            = "Pt1stJet/(Pt1stJet+MET) - Eta1stJetBump"
  plot916.ytit            = "Number of events"
  plot916.ylog            = "no"
  plot916.rebin           = 5
  #plot916.xmin            = 0
  #plot916.xmax            = 1
  plot916.ymin            = 0
  plot916.ymax            = 20
  #plot916.lpos = "bottom-center"
  plot916.name            = "1stJet_PTOverPTPlusMET_presel_Eta1stJetBump"
  plot916.addZUncBand     = zUncBand
  plot916.makeRatio       = makeRatio
  plot916.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MET_PAS_Eta1stJetBump ---
  variableName = "h1_MET_PAS_Eta1stJetBump"

  plot917 = Plot()
  ## inputs for stacked histograms
  plot917.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot917.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot917.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot917.keys            = keys
  plot917.xtit            = "pfMET - Eta1stJetBump"
  plot917.ytit            = "Number of events"
  plot917.ylog            = "no"
  plot917.rebin           = "var"
  #plot917.xmin            = 0
  #plot917.xmax            = 1000
  plot917.ymin            = 0
  plot917.ymax            = 20
  #plot917.lpos = "bottom-center"
  plot917.name            = "MET_presel_Eta1stJetBump"
  plot917.addZUncBand     = zUncBand
  plot917.makeRatio       = makeRatio
  plot917.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot917.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_Eta1stJetBump ---
  variableName = "h1_Njet_Eta1stJetBump"

  plot918 = Plot()
  ## inputs for stacked histograms
  plot918.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot918.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot918.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot918.keys            = keys
  plot918.xtit            = "Number of jets - Eta1stJetBump"
  plot918.ytit            = "Number of events"
  plot918.ylog            = "no"
  plot918.rebin           = 1
  plot918.xmin            = -0.5
  plot918.xmax            = 11.5
  plot918.ymin            = 0
  plot918.ymax            = 25
  #plot918.lpos = "bottom-center"
  plot918.name            = "Njet_presel_Eta1stJetBump"
  plot918.addZUncBand     = zUncBand
  plot918.makeRatio       = makeRatio
  plot918.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_Eta1stJetBump ---
  variableName = "h1_NjetTCHELBTag_Eta1stJetBump"

  plot919 = Plot()
  ## inputs for stacked histograms
  plot919.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot919.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot919.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot919.keys            = keys
  plot919.xtit            = "Number of b-tagged jets (TCHEL) - Eta1stJetBump"
  plot919.ytit            = "Number of events"
  plot919.ylog            = "no"
  plot919.rebin           = 1
  plot919.xmin            = -0.5
  plot919.xmax            = 11.5
  plot919.ymin            = 0
  plot919.ymax            = 25
  #plot919.lpos = "bottom-center"
  plot919.name            = "NjetTCHELBTag_presel_Eta1stJetBump"
  plot919.addZUncBand     = zUncBand
  plot919.makeRatio       = makeRatio
  plot919.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt1stJet_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Pt1stJet_PAS_OutsideEta1stJetBump"

  plot1200 = Plot()
  ## inputs for stacked histograms
  plot1200.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1200.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1200.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1200.keys            = keys
  plot1200.xtit            = "p_{T} 1^{st} jet - OutsideEta1stJetBump"
  plot1200.ytit            = "Number of events"
  plot1200.ylog            = "no"
  plot1200.rebin           = "var"
  #plot1200.xmin            = 0
  #plot1200.xmax            = 1000
  plot1200.ymin            = 0
  plot1200.ymax            = 50
  #plot1200.lpos = "bottom-center"
  plot1200.name            = "Pt1stJet_presel_OutsideEta1stJetBump"
  plot1200.addZUncBand     = zUncBand
  plot1200.makeRatio       = makeRatio
  plot1200.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot1200.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  ##--- h1_Eta1stJet_PAS_OutsideEta1stJetBump ---
  #variableName = "h1_Eta1stJet_PAS_OutsideEta1stJetBump"

  #plot1201 = Plot()
  ### inputs for stacked histograms
  #plot1201.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  #plot1201.keysStack       = keysStack
  ### this is the list of histograms that should be simply overlaid on top of the stacked histogram
  #plot1201.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  #plot1201.keys            = keys
  #plot1201.xtit            = "#eta 1^{st} jet - OutsideEta1stJetBump"
  #plot1201.ytit            = "Number of events"
  #plot1201.ylog            = "no"
  #plot1201.rebin           = 5
  ##plot1201.xmin            = -5
  ##plot1201.xmax            = 5
  #plot1201.ymin            = 0
  #plot1201.ymax            = 15
  ##plot1201.lpos = "bottom-center"
  #plot1201.name            = "Eta1stJet_presel_OutsideEta1stJetBump"
  #plot1201.addZUncBand     = zUncBand
  #plot1201.makeRatio       = makeRatio
  #plot1201.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stJet_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Phi1stJet_PAS_OutsideEta1stJetBump"

  plot1202 = Plot()
  ## inputs for stacked histograms
  plot1202.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1202.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1202.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1202.keys            = keys
  plot1202.xtit            = "#phi 1^{st} jet - OutsideEta1stJetBump"
  plot1202.ytit            = "Number of events"
  plot1202.ylog            = "no"
  plot1202.rebin           = 10
  #plot1202.xmin            = -3.1416
  #plot1202.xmax            = 3.1416
  plot1202.ymin            = 0
  plot1202.ymax            = 50
  #plot1202.lpos = "bottom-center"
  plot1202.name            = "Phi1stJet_presel_OutsideEta1stJetBump"
  plot1202.addZUncBand     = zUncBand
  plot1202.makeRatio       = makeRatio
  plot1202.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt2ndJet_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Pt2ndJet_PAS_OutsideEta1stJetBump"

  plot1203 = Plot()
  ## inputs for stacked histograms
  plot1203.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1203.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1203.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1203.keys            = keys
  plot1203.xtit            = "p_{T} 2^{nd} jet - OutsideEta1stJetBump"
  plot1203.ytit            = "Number of events"
  plot1203.ylog            = "no"
  plot1203.rebin           = "var"
  #plot1203.xmin            = 0
  #plot1203.xmax            = 1000
  plot1203.ymin            = 0
  plot1203.ymax            = 60
  #plot1203.lpos = "bottom-center"
  plot1203.name            = "Pt2ndJet_presel_OutsideEta1stJetBump"
  plot1203.addZUncBand     = zUncBand
  plot1203.makeRatio       = makeRatio
  plot1203.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot1203.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta2ndJet_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Eta2ndJet_PAS_OutsideEta1stJetBump"

  plot1204 = Plot()
  ## inputs for stacked histograms
  plot1204.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1204.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1204.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1204.keys            = keys
  plot1204.xtit            = "#eta 2^{nd} jet - OutsideEta1stJetBump"
  plot1204.ytit            = "Number of events"
  plot1204.ylog            = "no"
  plot1204.rebin           = 5
  #plot1204.xmin            = -5
  #plot1204.xmax            = 5
  plot1204.ymin            = 0
  plot1204.ymax            = 50
  #plot1204.lpos = "bottom-center"
  plot1204.name            = "Eta2ndJet_presel_OutsideEta1stJetBump"
  plot1204.addZUncBand     = zUncBand
  plot1204.makeRatio       = makeRatio
  plot1204.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi2ndJet_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Phi2ndJet_PAS_OutsideEta1stJetBump"

  plot1205 = Plot()
  ## inputs for stacked histograms
  plot1205.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1205.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1205.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1205.keys            = keys
  plot1205.xtit            = "#phi 2^{nd} jet - OutsideEta1stJetBump"
  plot1205.ytit            = "Number of events"
  plot1205.ylog            = "no"
  plot1205.rebin           = 10
  #plot1205.xmin            = -3.1416
  #plot1205.xmax            = 3.1416
  plot1205.ymin            = 0
  plot1205.ymax            = 50
  #plot1205.lpos = "bottom-center"
  plot1205.name            = "Phi2ndJet_presel_OutsideEta1stJetBump"
  plot1205.addZUncBand     = zUncBand
  plot1205.makeRatio       = makeRatio
  plot1205.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_E1stEle_PAS_OutsideEta1stJetBump ---
  variableName = "h1_E1stEle_PAS_OutsideEta1stJetBump"

  plot1206 = Plot()
  ## inputs for stacked histograms
  plot1206.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1206.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1206.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1206.keys            = keys
  plot1206.xtit            = "Energy 1^{st} electron - OutsideEta1stJetBump"
  plot1206.ytit            = "Number of events"
  plot1206.ylog            = "no"
  plot1206.rebin           = "var"
  #plot1206.xmin            = 0
  #plot1206.xmax            = 1000
  plot1206.ymin            = 0
  plot1206.ymax            = 60
  #plot1206.lpos = "bottom-center"
  plot1206.name            = "E1stEle_presel_OutsideEta1stJetBump"
  plot1206.addZUncBand     = zUncBand
  plot1206.makeRatio       = makeRatio
  plot1206.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot1206.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt1stEle_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Pt1stEle_PAS_OutsideEta1stJetBump"

  plot1207 = Plot()
  ## inputs for stacked histograms
  plot1207.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1207.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1207.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1207.keys            = keys
  plot1207.xtit            = "p_{T} 1^{st} electron - OutsideEta1stJetBump"
  plot1207.ytit            = "Number of events"
  plot1207.ylog            = "no"
  plot1207.rebin           = "var"
  #plot1207.xmin            = 0
  #plot1207.xmax            = 1000
  plot1207.ymin            = 0
  plot1207.ymax            = 80
  #plot1207.lpos = "bottom-center"
  plot1207.name            = "Pt1stEle_presel_OutsideEta1stJetBump"
  plot1207.addZUncBand     = zUncBand
  plot1207.makeRatio       = makeRatio
  plot1207.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot1207.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stEle_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Eta1stEle_PAS_OutsideEta1stJetBump"

  plot1208 = Plot()
  ## inputs for stacked histograms
  plot1208.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1208.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1208.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1208.keys            = keys
  plot1208.xtit            = "#eta 1^{st} electron - OutsideEta1stJetBump"
  plot1208.ytit            = "Number of events"
  plot1208.ylog            = "no"
  plot1208.rebin           = 5
  #plot1208.xmin            = -5
  #plot1208.xmax            = 5
  plot1208.ymin            = 0
  plot1208.ymax            = 50
  #plot1208.lpos = "bottom-center"
  plot1208.name            = "Eta1stEle_presel_OutsideEta1stJetBump"
  plot1208.addZUncBand     = zUncBand
  plot1208.makeRatio       = makeRatio
  plot1208.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stEle_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Phi1stEle_PAS_OutsideEta1stJetBump"

  plot1209 = Plot()
  ## inputs for stacked histograms
  plot1209.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1209.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1209.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1209.keys            = keys
  plot1209.xtit            = "#phi 1^{st} electron - OutsideEta1stJetBump"
  plot1209.ytit            = "Number of events"
  plot1209.ylog            = "no"
  plot1209.rebin           = 10
  #plot1209.xmin            = -3.1416
  #plot1209.xmax            = 3.1416
  plot1209.ymin            = 0
  plot1209.ymax            = 60
  #plot1209.lpos = "bottom-center"
  plot1209.name            = "Phi1stEle_presel_OutsideEta1stJetBump"
  plot1209.addZUncBand     = zUncBand
  plot1209.makeRatio       = makeRatio
  plot1209.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Charge1stEle_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Charge1stEle_PAS_OutsideEta1stJetBump"

  plot1210 = Plot()
  ## inputs for stacked histograms
  plot1210.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1210.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1210.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1210.keys            = keys
  plot1210.xtit            = "Charge 1^{st} electron - OutsideEta1stJetBump"
  plot1210.ytit            = "Number of events"
  plot1210.ylog            = "no"
  #plot1210.rebin           = 1
  #plot1210.xmin            = -1.001
  #plot1210.xmax            = 1.001
  plot1210.ymin            = 0
  plot1210.ymax            = 140
  #plot1210.lpos = "bottom-center"
  plot1210.name            = "Charge1stEle_presel_OutsideEta1stJetBump"
  plot1210.addZUncBand     = zUncBand
  plot1210.makeRatio       = makeRatio
  plot1210.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Conversion1stEle_OutsideEta1stJetBump ---
  variableName = "h1_Conversion1stEle_OutsideEta1stJetBump"

  plot1211 = Plot()
  ## inputs for stacked histograms
  plot1211.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot1211.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1211.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1211.keys            = keys
  plot1211.xtit            = "Conversion flag 1^{st} electron - OutsideEta1stJetBump"
  plot1211.ytit            = "Number of events"
  plot1211.ylog            = "no"
  #plot1211.rebin           = 1
  #plot1211.xmin            = -0.5
  #plot1211.xmax            = 1.5
  plot1211.ymin            = 0
  plot1211.ymax            = 200
  #plot1211.lpos = "bottom-center"
  plot1211.name            = "Conversion1stEle_presel_OutsideEta1stJetBump"
  plot1211.addZUncBand     = zUncBand
  plot1211.makeRatio       = makeRatio
  plot1211.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMETEle_OutsideEta1stJetBump ---
  variableName = "h1_mDeltaPhiMETEle_OutsideEta1stJetBump"

  plot1212 = Plot()
  ## inputs for stacked histograms
  plot1212.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1212.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1212.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1212.keys            = keys
  plot1212.xtit            = "#Delta#phi(MET,1^{st} ele) - OutsideEta1stJetBump"
  plot1212.ytit            = "Number of events"
  plot1212.ylog            = "no"
  plot1212.rebin           = 10
  #plot1212.xmin            = 0
  #plot1212.xmax            = 3.1416
  plot1212.ymin            = 0
  plot1212.ymax            = 60
  #plot1212.lpos = "bottom-center"
  plot1212.name            = "mDeltaPhiMETEle_presel_OutsideEta1stJetBump"
  plot1212.addZUncBand     = zUncBand
  plot1212.makeRatio       = makeRatio
  plot1212.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle2ndJet_PAS_OutsideEta1stJetBump ---
  variableName = "h1_mDeltaPhiEle2ndJet_PAS_OutsideEta1stJetBump"

  plot1213 = Plot()
  ## inputs for stacked histograms
  plot1213.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1213.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1213.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1213.keys            = keys
  plot1213.xtit            = "#Delta#phi(1^{st} ele,2^{nd} jet) - OutsideEta1stJetBump"
  plot1213.ytit            = "Number of events"
  plot1213.ylog            = "no"
  plot1213.rebin           = 10
  #plot1213.xmin            = 0
  #plot1213.xmax            = 3.1416
  plot1213.ymin            = 0
  plot1213.ymax            = 50
  #plot1213.lpos = "bottom-center"
  plot1213.name            = "mDeltaPhiEle2ndJet_presel_OutsideEta1stJetBump"
  plot1213.addZUncBand     = zUncBand
  plot1213.makeRatio       = makeRatio
  plot1213.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_OutsideEta1stJetBump ---
  variableName = "h1_mDeltaPhiMET2ndJet_OutsideEta1stJetBump"

  plot1214 = Plot()
  ## inputs for stacked histograms
  plot1214.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1214.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1214.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1214.keys            = keys
  plot1214.xtit            = "#Delta#phi(MET,2^{nd} jet) - OutsideEta1stJetBump"
  plot1214.ytit            = "Number of events"
  plot1214.ylog            = "no"
  plot1214.rebin           = 10
  #plot1214.xmin            = 0
  #plot1214.xmax            = 3.1416
  plot1214.ymin            = 0
  plot1214.ymax            = 50
  #plot1214.lpos = "bottom-center"
  plot1214.name            = "mDeltaPhiMET2ndJet_presel_OutsideEta1stJetBump"
  plot1214.addZUncBand     = zUncBand
  plot1214.makeRatio       = makeRatio
  plot1214.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Ptenu_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Ptenu_PAS_OutsideEta1stJetBump"

  plot1215 = Plot()
  ## inputs for stacked histograms
  plot1215.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1215.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1215.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1215.keys            = keys
  plot1215.xtit            = "p_{T}(e#nu) - OutsideEta1stJetBump"
  plot1215.ytit            = "Number of events"
  plot1215.ylog            = "no"
  plot1215.rebin           = "var"
  #plot1215.xmin            = 0
  #plot1215.xmax            = 2000
  plot1215.ymin            = 0
  plot1215.ymax            = 70
  #plot1215.lpos = "bottom-center"
  plot1215.name            = "Ptenu_presel_OutsideEta1stJetBump"
  plot1215.addZUncBand     = zUncBand
  plot1215.makeRatio       = makeRatio
  plot1215.xbins           = [0,40,80,120,160,200,300,600]
  plot1215.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_1stJet_PTOverPTPlusMET_OutsideEta1stJetBump ---
  variableName = "h1_1stJet_PTOverPTPlusMET_OutsideEta1stJetBump"

  plot1216 = Plot()
  ## inputs for stacked histograms
  plot1216.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1216.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1216.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1216.keys            = keys
  plot1216.xtit            = "Pt1stJet/(Pt1stJet+MET) - OutsideEta1stJetBump"
  plot1216.ytit            = "Number of events"
  plot1216.ylog            = "no"
  plot1216.rebin           = 5
  #plot1216.xmin            = 0
  #plot1216.xmax            = 1
  plot1216.ymin            = 0
  plot1216.ymax            = 60
  #plot1216.lpos = "bottom-center"
  plot1216.name            = "1stJet_PTOverPTPlusMET_presel_OutsideEta1stJetBump"
  plot1216.addZUncBand     = zUncBand
  plot1216.makeRatio       = makeRatio
  plot1216.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MET_PAS_OutsideEta1stJetBump ---
  variableName = "h1_MET_PAS_OutsideEta1stJetBump"

  plot1217 = Plot()
  ## inputs for stacked histograms
  plot1217.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1217.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1217.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1217.keys            = keys
  plot1217.xtit            = "pfMET - OutsideEta1stJetBump"
  plot1217.ytit            = "Number of events"
  plot1217.ylog            = "no"
  plot1217.rebin           = "var"
  #plot1217.xmin            = 0
  #plot1217.xmax            = 1000
  plot1217.ymin            = 0
  plot1217.ymax            = 70
  #plot1217.lpos = "bottom-center"
  plot1217.name            = "MET_presel_OutsideEta1stJetBump"
  plot1217.addZUncBand     = zUncBand
  plot1217.makeRatio       = makeRatio
  plot1217.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot1217.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_OutsideEta1stJetBump ---
  variableName = "h1_Njet_OutsideEta1stJetBump"

  plot1218 = Plot()
  ## inputs for stacked histograms
  plot1218.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1218.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1218.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1218.keys            = keys
  plot1218.xtit            = "Number of jets - OutsideEta1stJetBump"
  plot1218.ytit            = "Number of events"
  plot1218.ylog            = "no"
  plot1218.rebin           = 1
  plot1218.xmin            = -0.5
  plot1218.xmax            = 11.5
  plot1218.ymin            = 0
  plot1218.ymax            = 100
  #plot1218.lpos = "bottom-center"
  plot1218.name            = "Njet_presel_OutsideEta1stJetBump"
  plot1218.addZUncBand     = zUncBand
  plot1218.makeRatio       = makeRatio
  plot1218.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_OutsideEta1stJetBump ---
  variableName = "h1_NjetTCHELBTag_OutsideEta1stJetBump"

  plot1219 = Plot()
  ## inputs for stacked histograms
  plot1219.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1219.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1219.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1219.keys            = keys
  plot1219.xtit            = "Number of b-tagged jets (TCHEL) - OutsideEta1stJetBump"
  plot1219.ytit            = "Number of events"
  plot1219.ylog            = "no"
  plot1219.rebin           = 1
  plot1219.xmin            = -0.5
  plot1219.xmax            = 11.5
  plot1219.ymin            = 0
  plot1219.ymax            = 80
  #plot1219.lpos = "bottom-center"
  plot1219.name            = "NjetTCHELBTag_presel_OutsideEta1stJetBump"
  plot1219.addZUncBand     = zUncBand
  plot1219.makeRatio       = makeRatio
  plot1219.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


#-----------------------------------------------------------------------------------


# List of plots to be plotted
plots  = [plot0, plot1, plot2, plot_after2, plot3, plot4, plot5, plot_TCHELBTag, plot6, plot7, plot8, plot9, plot_TCHE1, plot_TCHE2
          , plot6and8, plot7and9
          , plot10, plot11, plot12, plot13, plot13_afterOtherDfCuts, plot14, plot14_ylog
          , plot14_plus, plot14_plus_ylog, plot14_minus, plot14_minus_ylog
          , plot14_0_1, plot14_1_2, plot14_2_pi
          , plot15, plot15_lep, plot15_jet, plot16
          , plot17, plot17_ylog, plot18
          , plot19, plot19_EleBarrel, plot19_EleEndcap, plot20, plot21, plot22, plot22_EleBarrel, plot22_EleEndcap, plot23, plot24, plot25
          , plot30, plot31, plot32, plot33, plot34, plot35
          ]

if doExtraPlots:
  extra_plots = [plot_Vtxd0, plot_MissHits, plot_Dist, plot_DCotTheta
                ,plot100, plot101, plot102, plot103, plot104, plot105, plot106, plot107, plot108, plot109, plot110
                ,plot200, plot201, plot202, plot203, plot204, plot205, plot206, plot207, plot208, plot209, plot210
                ,plot300, plot301, plot302, plot303, plot304, plot305, plot306, plot307, plot308, plot309, plot310
                ,plot400, plot401, plot402, plot403, plot404, plot405, plot406, plot407, plot408, plot409, plot410
                ,plot500, plot501, plot502, plot503, plot504, plot505, plot506, plot507, plot508, plot509, plot510
                ,plot600, plot601, plot602, plot603, plot604, plot605, plot606, plot607, plot608, plot609, plot610, plot611
                ,plot700, plot701, plot702, plot703, plot704, plot705, plot706, plot707, plot708, plot709, plot710, plot711, plot712, plot713, plot714, plot715, plot716, plot717, plot718, plot719
                ,plot1000, plot1001, plot1002, plot1003, plot1004, plot1005, plot1006, plot1007, plot1008, plot1009, plot1010, plot1011, plot1012, plot1013, plot1014, plot1015, plot1016, plot1017, plot1018, plot1019
                ,plot1100, plot1101, plot1102, plot1103, plot1104, plot1105, plot1106, plot1107, plot1108, plot1109, plot1110, plot1111, plot1112, plot1113, plot1114, plot1115, plot1116, plot1117, plot1118, plot1119
                ,plot800, plot801, plot802, plot803, plot804, plot805, plot806, plot807, plot808, plot809, plot810, plot811, plot812, plot813, plot814, plot815, plot816, plot817, plot818, plot819
                ,plot900, plot902, plot903, plot904, plot905, plot906, plot907, plot908, plot909, plot910, plot911, plot912, plot913, plot914, plot915, plot916, plot917, plot918, plot919
                ,plot1200, plot1202, plot1203, plot1204, plot1205, plot1206, plot1207, plot1208, plot1209, plot1210, plot1211, plot1212, plot1213, plot1214, plot1215, plot1216, plot1217, plot1218, plot1219
                ]

  plots = plots + extra_plots



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

