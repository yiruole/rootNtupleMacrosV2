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
        # xlog may not work if (self.xlog == "yes"):
        # fPads1.SetLogx()
        if (self.ylog     == "yes"):
            fPads1.SetLogy()

        #-- legend
        hsize=0.22
        vsize=0.26
        if (self.lpos=="bottom-center"):
            xstart=0.35
            ystart=0.25
        elif(self.lpos=="top-left"):
            xstart=0.12
            ystart=0.63
        else:
            xstart=0.65
            ystart=0.63
        legend = TLegend(xstart, ystart, xstart+hsize, ystart+vsize)
        legend.SetFillColor(kWhite)
        legend.SetMargin(0.2)
        legend.SetTextFont(132)

        #-- loop over histograms (stacked)
        Nstacked = len(self.histosStack)
        #stackColorIndexes = [20,38,14,45,20,38,14,45]
        #stackColorIndexes = [20,38,12,14,20,38,12,14]
        stackColorIndexes = [2,4,3,14,2,4,3,14]
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
                stack[iter].GetXaxis().SetLabelFont(132)
                stack[iter].GetXaxis().SetTitleOffset(1.0)
                stack[iter].GetXaxis().SetTitleSize(0.05)
                stack[iter].GetXaxis().SetLabelSize(0.045)
                stack[iter].GetYaxis().SetTitleFont(132)
                stack[iter].GetYaxis().SetLabelFont(132)
                stack[iter].GetYaxis().SetTitleOffset(0.8)
                stack[iter].GetYaxis().SetTitleSize(0.05)
                stack[iter].GetYaxis().SetLabelSize(0.045)
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
        dataLineIndexes = [1,2,3,1,2,3]
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
            legend.AddEntry(self.histodata, "Data","p")
            self.histodata.Draw("psame")

        #-- draw label
        l = TLatex()
        l.SetTextAlign(12)
        l.SetTextFont(132)
        l.SetTextSize(0.04)
        l.SetNDC()
#        l.DrawLatex(xstart,ystart-0.05,"CMS Preliminary 2010")
#        l.DrawLatex(xstart,ystart-0.10,"L_{int} = " + self.lint)
        if (self.lpos=="bottom-center"):
            l.DrawLatex(0.35,0.20,"CMS Preliminary")
            l.DrawLatex(0.35,0.09,"#intLdt = " + self.lint)
        if (self.lpos=="top-left"):
            l.DrawLatex(xstart+hsize+0.02,ystart+vsize-0.03,"CMS Preliminary")
            l.DrawLatex(xstart+hsize+0.02,ystart+vsize-0.11,"#intLdt = " + self.lint)
        else:
            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.03,"CMS Preliminary")
            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.11,"#intLdt = " + self.lint)

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

File_preselection = GetFile("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/36.0pb-1_presel_MET45_presel_sT250_Wrescale1.18_Feb112011/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/36.0pb-1_presel_MET45_presel_sT250_Wrescale1.16_PuMETSmearing_FullNtuples_Feb182011/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_enujjskim_MejStudies_v3/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_enujjskim_MET45_Jan11Prod/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_fullntuples_MET45_Jan11Prod/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_fullntuples_MET45_Jan11Prod_type1PFMET/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_fullntuples_MET45/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")
#File_preselection = GetFile("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/36.05pb-1_sT_presel_250_Zrescale1.20/analysisClass_enujjSample_plots.root")
##File_preselection = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/33.2pb-1_sT_presel_250_Zrescale1.20_Wrescale1.06_extraPlotsDec9/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")

File_selection    = File_preselection

UseQCDFromData    = 1 # always put an existing file under File_QCD (otherwise the code will crash)

File_QCD = GetFile("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/35.8pb-1_QCD_UseHLTPrescales_presel_MET45_presel_sT250_Feb112011/analysisClass_enujjSample_QCD_plots.root")
#File_QCD = GetFile("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/35.8pb-1_QCD_UseHLTPrescales_presel_MET45_presel_sT250_FullNtuples_Feb182011/analysisClass_enujjSample_QCD_plots.root")
#File_QCD  = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_enujjskim_MejStudies_v3/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")
#File_QCD  = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_enujjskim_MET45/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")
#File_QCD  = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_fullntuples_MET45/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")
#File_QCD  = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_fullntuples_MET45_type1PFMET/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")
#File_QCD  = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_fullntuples_MET45_more1SC/output_cutTable_enujjSample_QCD_more1SC/analysisClass_enujjSample_QCD_more1SC_plots.root")
#File_QCD  = GetFile("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/35.84pb-1_QCD_UseHLTPrescales_sT_presel_250/analysisClass_enujjSample_QCD_plots.root")

QCDscaleFactor    = 1 # no need to rescale anymore since we are using the HLT prescales (36/35.84 can be ignored)

#### Common values for plots:
#otherBkgsKey="QCD, single top, VV+jets, W/W*+jets"
otherBkgsKey="Other backgrounds"
zUncBand="no"
makeRatio=0
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
keysStack =             ["QCD",otherBkgsKey,"t#bar{t}", "W/W* + jets"]

samplesForHistos = ["LQenujj_M250", "LQenujj_M300","LQenujj_M340"]
keys             = ["LQ, M=250 GeV","LQ, M=300 GeV","LQ, M=340 GeV"]

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
plot1.xtit            = "p_{T} 1^{st} electron [GeV]"
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
plot2.xtit            = "#eta 1^{st} electron"
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
plot3.xtit            = "pfMET [GeV]"
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
plot4.xtit            = "min(p_{T} 1^{st} electron,pfMET) [GeV]"
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
plot6.xtit            = "p_{T} 1^{st} jet [GeV]"
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
plot7.xtit            = "#eta 1^{st} jet"
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
plot8.xtit            = "p_{T} 2^{nd} jet [GeV]"
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
plot9.xtit            = "#eta 2^{nd} jet"
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
plot_TCHE1.xtit            = "TCHE discriminator 1^{st} jet"
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
plot_TCHE2.xtit            = "TCHE discriminator 2^{nd} jet"
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
plot6and8.xtit            = "p_{T} jets [GeV]"
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
variableName = "nMuon_PtCut_IDISO"

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
plot11.xtit            = "#Delta#phi(MET,e) [rad.]"
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
plot12.xtit            = "#Delta#phi(MET,1^{st} jet) [rad.]"
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
plot13.xtit            = "#Delta#phi(MET,2^{nd} jet) [rad.]"
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
plot13_afterOtherDfCuts.xtit            = "#Delta#phi(MET,2^{nd} jet) [rad.]"
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
plot14.xtit            = "M_{T}(e#nu) [GeV]"
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
plot14_ylog.xtit            = "M_{T}(e#nu) [GeV]"
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
plot14_0_1.xtit            = "M_{T}(e#nu) [GeV] - #Delta#phi(MET,2nd jet)<=1"
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
plot14_1_2.xtit            = "M_{T}(e#nu) [GeV] - 1<#Delta#phi(MET,2nd jet)<2"
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
plot14_2_pi.xtit            = "M_{T}(e#nu) [GeV] - #Delta#phi(MET,2nd jet)>=2"
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
plot14_plus.xtit            = "M_{T}(e+#nu) [GeV]"
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
plot14_plus_ylog.xtit            = "M_{T}(e+#nu) [GeV]"
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
plot14_minus.xtit            = "M_{T}(e-#nu) [GeV]"
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
plot14_minus_ylog.xtit            = "M_{T}(e-#nu) [GeV]"
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
plot15.xtit            = "S_{T} [GeV]"
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
plot15_lep.xtit            = "S_{T} leptons [GeV]"
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
plot15_jet.xtit            = "S_{T} jets [GeV]"
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
plot16.xtit            = "M(jj) [GeV]"
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
plot17.xtit            = "M(ej) [GeV]"
plot17.ytit            = "Number of entries"
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
plot17_ylog.xtit            = "M(ej) [GeV]"
plot17_ylog.ytit            = "Number of entries"
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


#--- max Mej (after preselection) ---
variableName = "Mej_1stPair_PAS"

plot17_1 = Plot()
## inputs for stacked histograms
plot17_1.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot17_1.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot17_1.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot17_1.keys            = keys
plot17_1.xtit            = "max M(ej) [GeV]"
plot17_1.ytit            = "Number of entries"
plot17_1.ylog            = "no"
plot17_1.rebin           = "var"
plot17_1.xmin            = 0
plot17_1.xmax            = 1500
plot17_1.ymin            = 0.
plot17_1.ymax            = 150
#plot17_1.lpos = "bottom-center"
plot17_1.name            = "Mej_max_allPreviousCuts_ylin"
plot17_1.addZUncBand     = zUncBand
plot17_1.makeRatio       = makeRatio
plot17_1.xbins           = [0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,325,350,375,400,425,450,475,500,550,600,800,1000,1500]
plot17_1.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

plot17_1_ylog = Plot()
## inputs for stacked histograms
plot17_1_ylog.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot17_1_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot17_1_ylog.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot17_1_ylog.keys            = keys
plot17_1_ylog.xtit            = "max M(ej) [GeV]"
plot17_1_ylog.ytit            = "Number of entries"
plot17_1_ylog.ylog            = "yes"
plot17_1_ylog.rebin           = "var"
plot17_1_ylog.xmin            = 0
plot17_1_ylog.xmax            = 2000
plot17_1_ylog.ymin            = 0.01
plot17_1_ylog.ymax            = 1000
#plot17_1_ylog.lpos = "bottom-center"
plot17_1_ylog.name            = "Mej_max_allPreviousCuts"
plot17_1_ylog.addZUncBand     = zUncBand
plot17_1_ylog.makeRatio       = makeRatio
plot17_1_ylog.xbins           = [0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,325,350,375,400,425,450,475,500,550,600,800,1500,2000]
plot17_1_ylog.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- min Mej (after preselection) ---
variableName = "Mej_2ndPair_PAS"

plot17_2 = Plot()
## inputs for stacked histograms
plot17_2.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot17_2.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot17_2.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot17_2.keys            = keys
plot17_2.xtit            = "min M(ej) [GeV]"
plot17_2.ytit            = "Number of entries"
plot17_2.ylog            = "no"
plot17_2.rebin           = "var"
plot17_2.xmin            = 0
plot17_2.xmax            = 1500
plot17_2.ymin            = 0.
plot17_2.ymax            = 150
#plot17_2.lpos = "bottom-center"
plot17_2.name            = "Mej_min_allPreviousCuts_ylin"
plot17_2.addZUncBand     = zUncBand
plot17_2.makeRatio       = makeRatio
plot17_2.xbins           = [0,20,40,60,80,100,120,140,160,200,240,300,350,400,600,1000]
plot17_2.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

plot17_2_ylog = Plot()
## inputs for stacked histograms
plot17_2_ylog.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot17_2_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot17_2_ylog.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot17_2_ylog.keys            = keys
plot17_2_ylog.xtit            = "min M(ej) [GeV]"
plot17_2_ylog.ytit            = "Number of entries"
plot17_2_ylog.ylog            = "yes"
plot17_2_ylog.rebin           = "var"
plot17_2_ylog.xmin            = 0
plot17_2_ylog.xmax            = 2000
plot17_2_ylog.ymin            = 0.01
plot17_2_ylog.ymax            = 1000
#plot17_2_ylog.lpos = "bottom-center"
plot17_2_ylog.name            = "Mej_min_allPreviousCuts"
plot17_2_ylog.addZUncBand     = zUncBand
plot17_2_ylog.makeRatio       = makeRatio
plot17_2_ylog.xbins           = [0,20,40,60,80,100,120,140,160,200,240,300,350,400,600,1000]
plot17_2_ylog.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


##--- MTnuj (after preselection) ---
variableNames = ["MTnuj_1stPair_PAS","MTnuj_2ndPair_PAS"]

plot18 = Plot()
## inputs for stacked histograms
plot18.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistosQCD, variableNames, File_QCD, QCDscaleFactor) + generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot18.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot18.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot18.keys            = keys
plot18.xtit            = "M_{T}(#nuj) [GeV]"
plot18.ytit            = "Number of entries"
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
plot19.xtit            = "#phi 1^{st} electron [rad.]"
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
plot19_EleBarrel.xtit            = "#phi 1^{st} electron in barrel [rad.]"
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
plot19_EleEndcap.xtit            = "#phi 1^{st} electron in endcap [rad.]"
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
plot20.xtit            = "#phi 1^{st} jet [rad.]"
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
plot21.xtit            = "#phi 2^{nd} jet [rad.]"
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
plot22.xtit            = "pfMET #phi [rad.]"
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
plot22_EleBarrel.xtit            = "pfMET #phi for electron in barrel [rad.]"
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
plot22_EleEndcap.xtit            = "pfMET #phi for electron in endcap [rad.]"
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
plot23.xtit            = "p_{T}(e#nu) [GeV]"
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
#variableName = "sT"
variableName = "sT_MLQ250"

plot30 = Plot()
## inputs for stacked histograms
plot30.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot30.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot30.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot30.keys            = keys
plot30.xtit            = "S_{T} [GeV]"
plot30.ytit            = "Number of events"
plot30.ylog            = "yes"
plot30.rebin           = 10
plot30.xmin            = 100
plot30.xmax            = 2000
plot30.ymin            = 0.01
plot30.ymax            = 100
#plot30.lpos = "bottom-center"
#plot30.name            = "sT_fullSelection"
plot30.name            = "sT_fullSelection_M200"
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
plot31.xtit            = "M(ej) [GeV]"
plot31.ytit            = "Number of entries"
plot31.ylog            = "yes"
plot31.rebin           = 2
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
plot32.xtit            = "M_{T}(#nuj) [GeV]"
plot32.ytit            = "Number of entries"
plot32.ylog            = "yes"
plot32.rebin           = 2
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
plot33.xtit            = "p_{T}(e#nu) [GeV]"
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
plot35.xtit            = "#eta 1^{st} electron"
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
  plot_Vtxd0.xtit            = "Vtx d0 1^{st} electron [cm]"
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
  plot_MissHits.xtit            = "Missing Hits 1^{st} electron"
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
  plot_Dist.xtit            = "|Dist| 1^{st} electron"
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
  plot_DCotTheta.xtit            = "|#DeltaCot(#theta)| 1^{st} electron"
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
  plot100.xtit            = "#Delta#phi(MET,e) [rad.] - fullSel"
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
  plot101.xtit            = "#Delta#phi(MET,1^{st} jet) [rad.] - fullSel"
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
  plot102.xtit            = "#Delta#phi(MET,2^{nd} jet) [rad.] - fullSel"
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
  plot107.xtit            = "Vtx d0 1^{st} electron [cm] - fullSel"
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
  plot108.xtit            = "Missing Hits 1^{st} electron - fullSel"
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
  plot109.xtit            = "|Dist| 1^{st} electron - fullSel"
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
  plot110.xtit            = "|#DeltaCot(#theta)| 1^{st} electron - fullSel"
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
  plot200.xtit            = "#Delta#phi(MET,e) [rad.] - highMT"
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
  plot201.xtit            = "#Delta#phi(MET,1^{st} jet) [rad.] - highMT"
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
  plot202.xtit            = "#Delta#phi(MET,2^{nd} jet) [rad.] - highMT"
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
  plot207.xtit            = "Vtx d0 1^{st} electron [cm] - highMT"
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
  plot208.xtit            = "Missing Hits 1^{st} electron - highMT"
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
  plot209.xtit            = "|Dist| 1^{st} electron - highMT"
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
  plot210.xtit            = "|#DeltaCot(#theta)| 1^{st} electron - highMT"
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
  plot300.xtit            = "#Delta#phi(MET,e) [rad.] - highPt1stEle"
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
  plot301.xtit            = "#Delta#phi(MET,1^{st} jet) [rad.] - highPt1stEle"
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
  plot302.xtit            = "#Delta#phi(MET,2^{nd} jet) [rad.] - highPt1stEle"
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
  plot307.xtit            = "Vtx d0 1^{st} electron [cm] - highPt1stEle"
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
  plot308.xtit            = "Missing Hits 1^{st} electron - highPt1stEle"
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
  plot309.xtit            = "|Dist| 1^{st} electron - highPt1stEle"
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
  plot310.xtit            = "|#DeltaCot(#theta)| 1^{st} electron - highPt1stEle"
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
  plot400.xtit            = "#Delta#phi(MET,e) [rad.] - highPt1stJet"
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
  plot401.xtit            = "#Delta#phi(MET,1^{st} jet) [rad.] - highPt1stJet"
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
  plot402.xtit            = "#Delta#phi(MET,2^{nd} jet) [rad.] - highPt1stJet"
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
  plot407.xtit            = "Vtx d0 1^{st} electron [cm] - highPt1stJet"
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
  plot408.xtit            = "Missing Hits 1^{st} electron - highPt1stJet"
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
  plot409.xtit            = "|Dist| 1^{st} electron - highPt1stJet"
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
  plot410.xtit            = "|#DeltaCot(#theta)| 1^{st} electron - highPt1stJet"
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
  plot500.xtit            = "#Delta#phi(MET,e) [rad.] - highPt2ndJet"
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
  plot501.xtit            = "#Delta#phi(MET,1^{st} jet) [rad.] - highPt2ndJet"
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
  plot502.xtit            = "#Delta#phi(MET,2^{nd} jet) [rad.] - highPt2ndJet"
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
  plot507.xtit            = "Vtx d0 1^{st} electron [cm] - highPt2ndJet"
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
  plot508.xtit            = "Missing Hits 1^{st} electron - highPt2ndJet"
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
  plot509.xtit            = "|Dist| 1^{st} electron - highPt2ndJet"
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
  plot510.xtit            = "|#DeltaCot(#theta)| 1^{st} electron - highPt2ndJet"
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
  plot600.xtit            = "#Delta#phi(MET,e) [rad.] - highMej"
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
  plot601.xtit            = "#Delta#phi(MET,1^{st} jet) [rad.] - highMej"
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
  plot602.xtit            = "#Delta#phi(MET,2^{nd} jet) [rad.] - highMej"
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
  plot607.xtit            = "Vtx d0 1^{st} electron [cm] - highMej"
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
  plot608.xtit            = "Missing Hits 1^{st} electron - highMej"
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
  plot609.xtit            = "|Dist| 1^{st} electron - highMej"
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
  plot610.xtit            = "|#DeltaCot(#theta)| 1^{st} electron - highMej"
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
  plot611.xtit            = "Conversion flag 1^{st} electron - highMej"
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


  #--- h1_CHF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_CHF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot703 = Plot()
  ## inputs for stacked histograms
  plot703.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot703.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot703.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot703.keys            = keys
  plot703.xtit            = "Charged hadron energy fraction 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot703.ytit            = "Number of events"
  plot703.ylog            = "no"
  plot703.rebin           = 10
  #plot703.xmin            = 0
  #plot703.xmax            = 1
  plot703.ymin            = 0
  plot703.ymax            = 20
  #plot703.lpos = "bottom-center"
  plot703.name            = "CHF1stJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot703.addZUncBand     = zUncBand
  plot703.makeRatio       = makeRatio
  plot703.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NHF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NHF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot704 = Plot()
  ## inputs for stacked histograms
  plot704.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot704.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot704.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot704.keys            = keys
  plot704.xtit            = "Neutral hadron energy fraction 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot704.ytit            = "Number of events"
  plot704.ylog            = "no"
  plot704.rebin           = 10
  #plot704.xmin            = 0
  #plot704.xmax            = 1
  plot704.ymin            = 0
  plot704.ymax            = 40
  #plot704.lpos = "bottom-center"
  plot704.name            = "NHF1stJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot704.addZUncBand     = zUncBand
  plot704.makeRatio       = makeRatio
  plot704.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CEF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_CEF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot705 = Plot()
  ## inputs for stacked histograms
  plot705.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot705.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot705.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot705.keys            = keys
  plot705.xtit            = "Charged EM energy fraction 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot705.ytit            = "Number of events"
  plot705.ylog            = "no"
  plot705.rebin           = 10
  #plot705.xmin            = 0
  #plot705.xmax            = 1
  plot705.ymin            = 0
  plot705.ymax            = 60
  #plot705.lpos = "bottom-center"
  plot705.name            = "CEF1stJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot705.addZUncBand     = zUncBand
  plot705.makeRatio       = makeRatio
  plot705.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NEF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NEF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot706 = Plot()
  ## inputs for stacked histograms
  plot706.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot706.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot706.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot706.keys            = keys
  plot706.xtit            = "Neutral EM energy fraction 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot706.ytit            = "Number of events"
  plot706.ylog            = "no"
  plot706.rebin           = 10
  #plot706.xmin            = 0
  #plot706.xmax            = 1
  plot706.ymin            = 0
  plot706.ymax            = 20
  #plot706.lpos = "bottom-center"
  plot706.name            = "NEF1stJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot706.addZUncBand     = zUncBand
  plot706.makeRatio       = makeRatio
  plot706.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NCH1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NCH1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot707 = Plot()
  ## inputs for stacked histograms
  plot707.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot707.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot707.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot707.keys            = keys
  plot707.xtit            = "Charged multiplicity 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot707.ytit            = "Number of events"
  plot707.ylog            = "no"
  plot707.rebin           = 10
  #plot707.xmin            = -0.5
  #plot707.xmax            = 99.5
  plot707.ymin            = 0
  plot707.ymax            = 40
  #plot707.lpos = "bottom-center"
  plot707.name            = "NCH1stJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot707.addZUncBand     = zUncBand
  plot707.makeRatio       = makeRatio
  plot707.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NN1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NN1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot708 = Plot()
  ## inputs for stacked histograms
  plot708.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot708.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot708.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot708.keys            = keys
  plot708.xtit            = "Neutral multiplicity 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot708.ytit            = "Number of events"
  plot708.ylog            = "no"
  plot708.rebin           = 10
  #plot708.xmin            = -0.5
  #plot708.xmax            = 99.5
  plot708.ymin            = 0
  plot708.ymax            = 30
  #plot708.lpos = "bottom-center"
  plot708.name            = "NN1stJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot708.addZUncBand     = zUncBand
  plot708.makeRatio       = makeRatio
  plot708.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NC1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NC1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot709 = Plot()
  ## inputs for stacked histograms
  plot709.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot709.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot709.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot709.keys            = keys
  plot709.xtit            = "Number of constituents 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot709.ytit            = "Number of events"
  plot709.ylog            = "no"
  plot709.rebin           = 10
  #plot709.xmin            = -0.5
  #plot709.xmax            = 99.5
  plot709.ymin            = 0
  plot709.ymax            = 40
  #plot709.lpos = "bottom-center"
  plot709.name            = "NC1stJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot709.addZUncBand     = zUncBand
  plot709.makeRatio       = makeRatio
  plot709.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot710 = Plot()
  ## inputs for stacked histograms
  plot710.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot710.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot710.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot710.keys            = keys
  plot710.xtit            = "p_{T} 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot710.ytit            = "Number of events"
  plot710.ylog            = "no"
  plot710.rebin           = "var"
  #plot710.xmin            = 0
  #plot710.xmax            = 1000
  plot710.ymin            = 0
  plot710.ymax            = 20
  #plot710.lpos = "bottom-center"
  plot710.name            = "Pt2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot710.addZUncBand     = zUncBand
  plot710.makeRatio       = makeRatio
  plot710.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot710.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot711 = Plot()
  ## inputs for stacked histograms
  plot711.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot711.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot711.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot711.keys            = keys
  plot711.xtit            = "#eta 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot711.ytit            = "Number of events"
  plot711.ylog            = "no"
  plot711.rebin           = 5
  #plot711.xmin            = -5
  #plot711.xmax            = 5
  plot711.ymin            = 0
  plot711.ymax            = 15
  #plot711.lpos = "bottom-center"
  plot711.name            = "Eta2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot711.addZUncBand     = zUncBand
  plot711.makeRatio       = makeRatio
  plot711.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot712 = Plot()
  ## inputs for stacked histograms
  plot712.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot712.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot712.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot712.keys            = keys
  plot712.xtit            = "#phi 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot712.ytit            = "Number of events"
  plot712.ylog            = "no"
  plot712.rebin           = 10
  #plot712.xmin            = -3.1416
  #plot712.xmax            = 3.1416
  plot712.ymin            = 0
  plot712.ymax            = 20
  #plot712.lpos = "bottom-center"
  plot712.name            = "Phi2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot712.addZUncBand     = zUncBand
  plot712.makeRatio       = makeRatio
  plot712.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CHF2ndJet_highMej_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_CHF2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot713 = Plot()
  ## inputs for stacked histograms
  plot713.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot713.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot713.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot713.keys            = keys
  plot713.xtit            = "Charged hadron energy fraction 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot713.ytit            = "Number of events"
  plot713.ylog            = "no"
  plot713.rebin           = 10
  #plot713.xmin            = 0
  #plot713.xmax            = 1
  plot713.ymin            = 0
  plot713.ymax            = 20
  #plot713.lpos = "bottom-center"
  plot713.name            = "CHF2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot713.addZUncBand     = zUncBand
  plot713.makeRatio       = makeRatio
  plot713.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NHF2ndJet_highMej_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_NHF2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot714 = Plot()
  ## inputs for stacked histograms
  plot714.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot714.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot714.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot714.keys            = keys
  plot714.xtit            = "Neutral hadron energy fraction 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot714.ytit            = "Number of events"
  plot714.ylog            = "no"
  plot714.rebin           = 10
  #plot714.xmin            = 0
  #plot714.xmax            = 1
  plot714.ymin            = 0
  plot714.ymax            = 40
  #plot714.lpos = "bottom-center"
  plot714.name            = "NHF2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot714.addZUncBand     = zUncBand
  plot714.makeRatio       = makeRatio
  plot714.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CEF2ndJet_highMej_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_CEF2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot715 = Plot()
  ## inputs for stacked histograms
  plot715.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot715.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot715.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot715.keys            = keys
  plot715.xtit            = "Charged EM energy fraction 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot715.ytit            = "Number of events"
  plot715.ylog            = "no"
  plot715.rebin           = 10
  #plot715.xmin            = 0
  #plot715.xmax            = 1
  plot715.ymin            = 0
  plot715.ymax            = 60
  #plot715.lpos = "bottom-center"
  plot715.name            = "CEF2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot715.addZUncBand     = zUncBand
  plot715.makeRatio       = makeRatio
  plot715.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NEF2ndJet_highMej_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_NEF2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot716 = Plot()
  ## inputs for stacked histograms
  plot716.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot716.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot716.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot716.keys            = keys
  plot716.xtit            = "Neutral EM energy fraction 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot716.ytit            = "Number of events"
  plot716.ylog            = "no"
  plot716.rebin           = 10
  #plot716.xmin            = 0
  #plot716.xmax            = 1
  plot716.ymin            = 0
  plot716.ymax            = 20
  #plot716.lpos = "bottom-center"
  plot716.name            = "NEF2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot716.addZUncBand     = zUncBand
  plot716.makeRatio       = makeRatio
  plot716.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NCH2ndJet_highMej_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_NCH2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot717 = Plot()
  ## inputs for stacked histograms
  plot717.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot717.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot717.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot717.keys            = keys
  plot717.xtit            = "Charged multiplicity 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot717.ytit            = "Number of events"
  plot717.ylog            = "no"
  plot717.rebin           = 10
  #plot717.xmin            = -0.5
  #plot717.xmax            = 99.5
  plot717.ymin            = 0
  plot717.ymax            = 40
  #plot717.lpos = "bottom-center"
  plot717.name            = "NCH2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot717.addZUncBand     = zUncBand
  plot717.makeRatio       = makeRatio
  plot717.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NN2ndJet_highMej_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_NN2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot718 = Plot()
  ## inputs for stacked histograms
  plot718.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot718.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot718.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot718.keys            = keys
  plot718.xtit            = "Neutral multiplicity 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot718.ytit            = "Number of events"
  plot718.ylog            = "no"
  plot718.rebin           = 10
  #plot718.xmin            = -0.5
  #plot718.xmax            = 99.5
  plot718.ymin            = 0
  plot718.ymax            = 30
  #plot718.lpos = "bottom-center"
  plot718.name            = "NN2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot718.addZUncBand     = zUncBand
  plot718.makeRatio       = makeRatio
  plot718.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NC2ndJet_highMej_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_NC2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot719 = Plot()
  ## inputs for stacked histograms
  plot719.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot719.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot719.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot719.keys            = keys
  plot719.xtit            = "Number of constituents 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot719.ytit            = "Number of events"
  plot719.ylog            = "no"
  plot719.rebin           = 10
  #plot719.xmin            = -0.5
  #plot719.xmax            = 99.5
  plot719.ymin            = 0
  plot719.ymax            = 30
  #plot719.lpos = "bottom-center"
  plot719.name            = "NC2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot719.addZUncBand     = zUncBand
  plot719.makeRatio       = makeRatio
  plot719.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_E1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_E1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot720 = Plot()
  ## inputs for stacked histograms
  plot720.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot720.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot720.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot720.keys            = keys
  plot720.xtit            = "Energy 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot720.ytit            = "Number of events"
  plot720.ylog            = "no"
  plot720.rebin           = "var"
  #plot720.xmin            = 0
  #plot720.xmax            = 1000
  plot720.ymin            = 0
  plot720.ymax            = 20
  #plot720.lpos = "bottom-center"
  plot720.name            = "E1stEle_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot720.addZUncBand     = zUncBand
  plot720.makeRatio       = makeRatio
  plot720.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot720.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot721 = Plot()
  ## inputs for stacked histograms
  plot721.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot721.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot721.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot721.keys            = keys
  plot721.xtit            = "p_{T} 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot721.ytit            = "Number of events"
  plot721.ylog            = "no"
  plot721.rebin           = "var"
  #plot721.xmin            = 0
  #plot721.xmax            = 1000
  plot721.ymin            = 0
  plot721.ymax            = 20
  #plot721.lpos = "bottom-center"
  plot721.name            = "Pt1stEle_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot721.addZUncBand     = zUncBand
  plot721.makeRatio       = makeRatio
  plot721.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot721.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot722 = Plot()
  ## inputs for stacked histograms
  plot722.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot722.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot722.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot722.keys            = keys
  plot722.xtit            = "#eta 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot722.ytit            = "Number of events"
  plot722.ylog            = "no"
  plot722.rebin           = 5
  #plot722.xmin            = -5
  #plot722.xmax            = 5
  plot722.ymin            = 0
  plot722.ymax            = 15
  #plot722.lpos = "bottom-center"
  plot722.name            = "Eta1stEle_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot722.addZUncBand     = zUncBand
  plot722.makeRatio       = makeRatio
  plot722.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot723 = Plot()
  ## inputs for stacked histograms
  plot723.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot723.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot723.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot723.keys            = keys
  plot723.xtit            = "#phi 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot723.ytit            = "Number of events"
  plot723.ylog            = "no"
  plot723.rebin           = 10
  #plot723.xmin            = -3.1416
  #plot723.xmax            = 3.1416
  plot723.ymin            = 0
  plot723.ymax            = 20
  #plot723.lpos = "bottom-center"
  plot723.name            = "Phi1stEle_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot723.addZUncBand     = zUncBand
  plot723.makeRatio       = makeRatio
  plot723.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Charge1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Charge1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot724 = Plot()
  ## inputs for stacked histograms
  plot724.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot724.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot724.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot724.keys            = keys
  plot724.xtit            = "Charge 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot724.ytit            = "Number of events"
  plot724.ylog            = "no"
  #plot724.rebin           = 1
  #plot724.xmin            = -1.001
  #plot724.xmax            = 1.001
  plot724.ymin            = 0
  plot724.ymax            = 50
  #plot724.lpos = "bottom-center"
  plot724.name            = "Charge1stEle_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot724.addZUncBand     = zUncBand
  plot724.makeRatio       = makeRatio
  plot724.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Conversion1stEle_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Conversion1stEle_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot725 = Plot()
  ## inputs for stacked histograms
  plot725.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot725.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot725.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot725.keys            = keys
  plot725.xtit            = "Conversion flag 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot725.ytit            = "Number of events"
  plot725.ylog            = "no"
  #plot725.rebin           = 1
  #plot725.xmin            = -0.5
  #plot725.xmax            = 1.5
  plot725.ymin            = 0
  plot725.ymax            = 60
  #plot725.lpos = "bottom-center"
  plot725.name            = "Conversion1stEle_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot725.addZUncBand     = zUncBand
  plot725.makeRatio       = makeRatio
  plot725.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMETEle_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiMETEle_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot726 = Plot()
  ## inputs for stacked histograms
  plot726.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot726.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot726.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot726.keys            = keys
  plot726.xtit            = "#Delta#phi(MET,1^{st} ele) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot726.ytit            = "Number of events"
  plot726.ylog            = "no"
  plot726.rebin           = 10
  #plot726.xmin            = 0
  #plot726.xmax            = 3.1416
  plot726.ymin            = 0
  plot726.ymax            = 20
  #plot726.lpos = "bottom-center"
  plot726.name            = "mDeltaPhiMETEle_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot726.addZUncBand     = zUncBand
  plot726.makeRatio       = makeRatio
  plot726.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiEle1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot727 = Plot()
  ## inputs for stacked histograms
  plot727.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot727.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot727.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot727.keys            = keys
  plot727.xtit            = "#Delta#phi(1^{st} ele,1^{st} jet) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot727.ytit            = "Number of events"
  plot727.ylog            = "no"
  plot727.rebin           = 10
  #plot727.xmin            = 0
  #plot727.xmax            = 3.1416
  plot727.ymin            = 0
  plot727.ymax            = 20
  #plot727.lpos = "bottom-center"
  plot727.name            = "mDeltaPhiEle1stJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot727.addZUncBand     = zUncBand
  plot727.makeRatio       = makeRatio
  plot727.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiEle2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot728 = Plot()
  ## inputs for stacked histograms
  plot728.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot728.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot728.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot728.keys            = keys
  plot728.xtit            = "#Delta#phi(1^{st} ele,2^{nd} jet) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot728.ytit            = "Number of events"
  plot728.ylog            = "no"
  plot728.rebin           = 10
  #plot728.xmin            = 0
  #plot728.xmax            = 3.1416
  plot728.ymin            = 0
  plot728.ymax            = 20
  #plot728.lpos = "bottom-center"
  plot728.name            = "mDeltaPhiEle2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot728.addZUncBand     = zUncBand
  plot728.makeRatio       = makeRatio
  plot728.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiMET2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot729 = Plot()
  ## inputs for stacked histograms
  plot729.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot729.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot729.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot729.keys            = keys
  plot729.xtit            = "#Delta#phi(MET,2^{nd} jet) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot729.ytit            = "Number of events"
  plot729.ylog            = "no"
  plot729.rebin           = 10
  #plot729.xmin            = 0
  #plot729.xmax            = 3.1416
  plot729.ymin            = 0
  plot729.ymax            = 20
  #plot729.lpos = "bottom-center"
  plot729.name            = "mDeltaPhiMET2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot729.addZUncBand     = zUncBand
  plot729.makeRatio       = makeRatio
  plot729.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Ptenu_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Ptenu_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot730 = Plot()
  ## inputs for stacked histograms
  plot730.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot730.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot730.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot730.keys            = keys
  plot730.xtit            = "p_{T}(e#nu) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot730.ytit            = "Number of events"
  plot730.ylog            = "no"
  plot730.rebin           = "var"
  #plot730.xmin            = 0
  #plot730.xmax            = 2000
  plot730.ymin            = 0
  plot730.ymax            = 20
  #plot730.lpos = "bottom-center"
  plot730.name            = "Ptenu_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot730.addZUncBand     = zUncBand
  plot730.makeRatio       = makeRatio
  plot730.xbins           = [0,40,80,120,160,200,300,600]
  plot730.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_1stJet_PTOverPTPlusMET_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_1stJet_PTOverPTPlusMET_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot731 = Plot()
  ## inputs for stacked histograms
  plot731.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot731.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot731.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot731.keys            = keys
  plot731.xtit            = "Pt1stJet/(Pt1stJet+MET) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot731.ytit            = "Number of events"
  plot731.ylog            = "no"
  plot731.rebin           = 5
  #plot731.xmin            = 0
  #plot731.xmax            = 1
  plot731.ymin            = 0
  plot731.ymax            = 20
  #plot731.lpos = "bottom-center"
  plot731.name            = "1stJet_PTOverPTPlusMET_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot731.addZUncBand     = zUncBand
  plot731.makeRatio       = makeRatio
  plot731.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MET_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_MET_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot732 = Plot()
  ## inputs for stacked histograms
  plot732.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot732.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot732.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot732.keys            = keys
  plot732.xtit            = "pfMET - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot732.ytit            = "Number of events"
  plot732.ylog            = "no"
  plot732.rebin           = "var"
  #plot732.xmin            = 0
  #plot732.xmax            = 1000
  plot732.ymin            = 0
  plot732.ymax            = 20
  #plot732.lpos = "bottom-center"
  plot732.name            = "MET_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot732.addZUncBand     = zUncBand
  plot732.makeRatio       = makeRatio
  plot732.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot732.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Njet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot733 = Plot()
  ## inputs for stacked histograms
  plot733.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot733.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot733.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot733.keys            = keys
  plot733.xtit            = "Number of jets - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot733.ytit            = "Number of events"
  plot733.ylog            = "no"
  plot733.rebin           = 1
  plot733.xmin            = -0.5
  plot733.xmax            = 11.5
  plot733.ymin            = 0
  plot733.ymax            = 25
  #plot733.lpos = "bottom-center"
  plot733.name            = "Njet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot733.addZUncBand     = zUncBand
  plot733.makeRatio       = makeRatio
  plot733.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NjetTCHELBTag_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot734 = Plot()
  ## inputs for stacked histograms
  plot734.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot734.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot734.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot734.keys            = keys
  plot734.xtit            = "Number of b-tagged jets (TCHEL) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot734.ytit            = "Number of events"
  plot734.ylog            = "no"
  plot734.rebin           = 1
  plot734.xmin            = -0.5
  plot734.xmax            = 11.5
  plot734.ymin            = 0
  plot734.ymax            = 25
  #plot734.lpos = "bottom-center"
  plot734.name            = "NjetTCHELBTag_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot734.addZUncBand     = zUncBand
  plot734.makeRatio       = makeRatio
  plot734.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_minDRej_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_minDRej_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot735 = Plot()
  ## inputs for stacked histograms
  plot735.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot735.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot735.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot735.keys            = keys
  plot735.xtit            = "min#DeltaR(e,j) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot735.ytit            = "Number of events"
  plot735.ylog            = "no"
  plot735.rebin           = 10
  #plot735.xmin            = 0
  #plot735.xmax            = 10
  plot735.ymin            = 0
  plot735.ymax            = 25
  #plot735.lpos = "bottom-center"
  plot735.name            = "minDRej_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot735.addZUncBand     = zUncBand
  plot735.makeRatio       = makeRatio
  plot735.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_maxDRej_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_maxDRej_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot736 = Plot()
  ## inputs for stacked histograms
  plot736.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot736.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot736.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot736.keys            = keys
  plot736.xtit            = "max#DeltaR(e,j) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot736.ytit            = "Number of events"
  plot736.ylog            = "no"
  plot736.rebin           = 10
  #plot736.xmin            = 0
  #plot736.xmax            = 10
  plot736.ymin            = 0
  plot736.ymax            = 30
  #plot736.lpos = "bottom-center"
  plot736.name            = "maxDRej_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot736.addZUncBand     = zUncBand
  plot736.makeRatio       = makeRatio
  plot736.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DRjets_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_DRjets_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot737 = Plot()
  ## inputs for stacked histograms
  plot737.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot737.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot737.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot737.keys            = keys
  plot737.xtit            = "#DeltaR(j1,j2) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot737.ytit            = "Number of events"
  plot737.ylog            = "no"
  plot737.rebin           = 10
  #plot737.xmin            = 0
  #plot737.xmax            = 10
  plot737.ymin            = 0
  plot737.ymax            = 30
  #plot737.lpos = "bottom-center"
  plot737.name            = "DRjets_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot737.addZUncBand     = zUncBand
  plot737.makeRatio       = makeRatio
  plot737.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Mass1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Mass1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot738 = Plot()
  ## inputs for stacked histograms
  plot738.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot738.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot738.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot738.keys            = keys
  plot738.xtit            = "Mass 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot738.ytit            = "Number of events"
  plot738.ylog            = "no"
  plot738.rebin           = 5
  #plot738.xmin            = -0.5
  #plot738.xmax            = 99.5
  plot738.ymin            = 0
  plot738.ymax            = 20
  #plot738.lpos = "bottom-center"
  plot738.name            = "Mass1stJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot738.addZUncBand     = zUncBand
  plot738.makeRatio       = makeRatio
  plot738.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

  #--- h1_minDR_1stPFjet_CaloJet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_minDR_1stPFjet_CaloJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot739 = Plot()
  ## inputs for stacked histograms
  plot739.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot739.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot739.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot739.keys            = keys
  plot739.xtit            = "min#DeltaR(PFj1,CaloJet) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot739.ytit            = "Number of events"
  plot739.ylog            = "no"
  plot739.rebin           = 1
  plot739.xmin            = 0
  plot739.xmax            = 0.5
  plot739.ymin            = 0
  plot739.ymax            = 20
  #plot739.lpos = "bottom-center"
  plot739.name            = "minDR_1stPFjet_CaloJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot739.addZUncBand     = zUncBand
  plot739.makeRatio       = makeRatio
  plot739.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

  #--- h1_PT1stJet_over_PTCaloJet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_PT1stJet_over_PTCaloJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot740 = Plot()
  ## inputs for stacked histograms
  plot740.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot740.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot740.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot740.keys            = keys
  plot740.xtit            = "PT(1^{st} PF jet)/PT(Closest Calo jet) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot740.ytit            = "Number of events"
  plot740.ylog            = "no"
  plot740.rebin           = 5
  #plot740.xmin            = 0
  #plot740.xmax            = 10
  plot740.ymin            = 0
  plot740.ymax            = 20
  #plot740.lpos = "bottom-center"
  plot740.name            = "PT1stJet_over_PTCaloJet_CaloJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot740.addZUncBand     = zUncBand
  plot740.makeRatio       = makeRatio
  plot740.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

  #--- h1_NC1stJet_over_n90HitsCaloJet_1stPFjet_CaloJet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NC1stJet_over_n90HitsCaloJet_1stPFjet_CaloJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot741 = Plot()
  ## inputs for stacked histograms
  plot741.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot741.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot741.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot741.keys            = keys
  plot741.xtit            = "NC(1^{st} PF jet)/n90Hits(Closest Calo jet) - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot741.ytit            = "Number of events"
  plot741.ylog            = "no"
  plot741.rebin           = 5
  #plot741.xmin            = 0
  #plot741.xmax            = 10
  plot741.ymin            = 0
  plot741.ymax            = 20
  #plot741.lpos = "bottom-center"
  plot741.name            = "NC1stJet_over_n90HitsCaloJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot741.addZUncBand     = zUncBand
  plot741.makeRatio       = makeRatio
  plot741.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Mass2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Mass2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot742 = Plot()
  ## inputs for stacked histograms
  plot742.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot742.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot742.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot742.keys            = keys
  plot742.xtit            = "Mass 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot742.ytit            = "Number of events"
  plot742.ylog            = "no"
  plot742.rebin           = 5
  #plot742.xmin            = -0.5
  #plot742.xmax            = 99.5
  plot742.ymin            = 0
  plot742.ymax            = 20
  #plot742.lpos = "bottom-center"
  plot742.name            = "Mass2ndJet_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot742.addZUncBand     = zUncBand
  plot742.makeRatio       = makeRatio
  plot742.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CaloMET_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_CaloMET_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot743 = Plot()
  ## inputs for stacked histograms
  plot743.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot743.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot743.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot743.keys            = keys
  plot743.xtit            = "caloMET - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot743.ytit            = "Number of events"
  plot743.ylog            = "no"
  plot743.rebin           = "var"
  #plot743.xmin            = 0
  #plot743.xmax            = 1000
  plot743.ymin            = 0
  plot743.ymax            = 20
  #plot743.lpos = "bottom-center"
  plot743.name            = "CaloMET_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot743.addZUncBand     = zUncBand
  plot743.makeRatio       = makeRatio
  plot743.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot743.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_TCMET_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_TCMET_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot744 = Plot()
  ## inputs for stacked histograms
  plot744.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot744.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot744.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot744.keys            = keys
  plot744.xtit            = "tcMET - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot744.ytit            = "Number of events"
  plot744.ylog            = "no"
  plot744.rebin           = "var"
  #plot744.xmin            = 0
  #plot744.xmax            = 1000
  plot744.ymin            = 0
  plot744.ymax            = 20
  #plot744.lpos = "bottom-center"
  plot744.name            = "TCMET_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot744.addZUncBand     = zUncBand
  plot744.makeRatio       = makeRatio
  plot744.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot744.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

  #--- h1_MET_over_CaloMET_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_MET_over_CaloMET_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot745 = Plot()
  ## inputs for stacked histograms
  plot745.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot745.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot745.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot745.keys            = keys
  plot745.xtit            = "pfMET/CaloMET - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot745.ytit            = "Number of events"
  plot745.ylog            = "no"
  plot745.rebin           = 5
  #plot745.xmin            = 0
  #plot745.xmax            = 10
  plot745.ymin            = 0
  plot745.ymax            = 30
  #plot745.lpos = "bottom-center"
  plot745.name            = "MET_over_CaloMET_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot745.addZUncBand     = zUncBand
  plot745.makeRatio       = makeRatio
  plot745.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

  #--- h1_MTenu_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_MTenu_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5"

  plot746 = Plot()
  ## inputs for stacked histograms
  plot746.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot746.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot746.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot746.keys            = keys
  plot746.xtit            = "MTenu - highMej, #Delta#phi(MET,1^{st} jet)>2.5"
  plot746.ytit            = "Number of events"
  plot746.ylog            = "no"
  plot746.rebin           = "var"
  #plot746.xmin            = 0
  #plot746.xmax            = 1000
  plot746.ymin            = 0
  plot746.ymax            = 20
  #plot746.lpos = "bottom-center"
  plot746.name            = "MTenu_presel_highMej_mDeltaPhiMET1stJet_gt_2.5"
  plot746.addZUncBand     = zUncBand
  plot746.makeRatio       = makeRatio
  plot746.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot746.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)




  # mDeltaPhiMET1stJet>2.5, e+ only

  #--- h1_Pt1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot800 = Plot()
  ## inputs for stacked histograms
  plot800.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot800.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot800.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot800.keys            = keys
  plot800.xtit            = "p_{T} 1^{st} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot800.ytit            = "Number of events"
  plot800.ylog            = "no"
  plot800.rebin           = "var"
  #plot800.xmin            = 0
  #plot800.xmax            = 1000
  plot800.ymin            = 0
  plot800.ymax            = 15
  #plot800.lpos = "bottom-center"
  plot800.name            = "Pt1stJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot800.addZUncBand     = zUncBand
  plot800.makeRatio       = makeRatio
  plot800.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot800.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot801 = Plot()
  ## inputs for stacked histograms
  plot801.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot801.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot801.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot801.keys            = keys
  plot801.xtit            = "#eta 1^{st} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot801.ytit            = "Number of events"
  plot801.ylog            = "no"
  plot801.rebin           = 5
  #plot801.xmin            = -5
  #plot801.xmax            = 5
  plot801.ymin            = 0
  plot801.ymax            = 15
  #plot801.lpos = "bottom-center"
  plot801.name            = "Eta1stJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot801.addZUncBand     = zUncBand
  plot801.makeRatio       = makeRatio
  plot801.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot802 = Plot()
  ## inputs for stacked histograms
  plot802.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot802.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot802.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot802.keys            = keys
  plot802.xtit            = "#phi 1^{st} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot802.ytit            = "Number of events"
  plot802.ylog            = "no"
  plot802.rebin           = 10
  #plot802.xmin            = -3.1416
  #plot802.xmax            = 3.1416
  plot802.ymin            = 0
  plot802.ymax            = 20
  #plot802.lpos = "bottom-center"
  plot802.name            = "Phi1stJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot802.addZUncBand     = zUncBand
  plot802.makeRatio       = makeRatio
  plot802.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CHF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_CHF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot803 = Plot()
  ## inputs for stacked histograms
  plot803.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot803.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot803.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot803.keys            = keys
  plot803.xtit            = "Charged hadron energy fraction 1^{st} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot803.ytit            = "Number of events"
  plot803.ylog            = "no"
  plot803.rebin           = 10
  #plot803.xmin            = 0
  #plot803.xmax            = 1
  plot803.ymin            = 0
  plot803.ymax            = 20
  #plot803.lpos = "bottom-center"
  plot803.name            = "CHF1stJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot803.addZUncBand     = zUncBand
  plot803.makeRatio       = makeRatio
  plot803.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NHF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NHF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot804 = Plot()
  ## inputs for stacked histograms
  plot804.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot804.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot804.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot804.keys            = keys
  plot804.xtit            = "Neutral hadron energy fraction 1^{st} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot804.ytit            = "Number of events"
  plot804.ylog            = "no"
  plot804.rebin           = 10
  #plot804.xmin            = 0
  #plot804.xmax            = 1
  plot804.ymin            = 0
  plot804.ymax            = 40
  #plot804.lpos = "bottom-center"
  plot804.name            = "NHF1stJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot804.addZUncBand     = zUncBand
  plot804.makeRatio       = makeRatio
  plot804.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CEF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_CEF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot805 = Plot()
  ## inputs for stacked histograms
  plot805.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot805.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot805.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot805.keys            = keys
  plot805.xtit            = "Charged EM energy fraction 1^{st} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot805.ytit            = "Number of events"
  plot805.ylog            = "no"
  plot805.rebin           = 10
  #plot805.xmin            = 0
  #plot805.xmax            = 1
  plot805.ymin            = 0
  plot805.ymax            = 60
  #plot805.lpos = "bottom-center"
  plot805.name            = "CEF1stJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot805.addZUncBand     = zUncBand
  plot805.makeRatio       = makeRatio
  plot805.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NEF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NEF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot806 = Plot()
  ## inputs for stacked histograms
  plot806.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot806.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot806.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot806.keys            = keys
  plot806.xtit            = "Neutral EM energy fraction 1^{st} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot806.ytit            = "Number of events"
  plot806.ylog            = "no"
  plot806.rebin           = 10
  #plot806.xmin            = 0
  #plot806.xmax            = 1
  plot806.ymin            = 0
  plot806.ymax            = 20
  #plot806.lpos = "bottom-center"
  plot806.name            = "NEF1stJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot806.addZUncBand     = zUncBand
  plot806.makeRatio       = makeRatio
  plot806.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NCH1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NCH1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot807 = Plot()
  ## inputs for stacked histograms
  plot807.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot807.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot807.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot807.keys            = keys
  plot807.xtit            = "Charged multiplicity 1^{st} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot807.ytit            = "Number of events"
  plot807.ylog            = "no"
  plot807.rebin           = 10
  #plot807.xmin            = -0.5
  #plot807.xmax            = 99.5
  plot807.ymin            = 0
  plot807.ymax            = 40
  #plot807.lpos = "bottom-center"
  plot807.name            = "NCH1stJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot807.addZUncBand     = zUncBand
  plot807.makeRatio       = makeRatio
  plot807.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NN1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NN1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot808 = Plot()
  ## inputs for stacked histograms
  plot808.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot808.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot808.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot808.keys            = keys
  plot808.xtit            = "Neutral multiplicity 1^{st} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot808.ytit            = "Number of events"
  plot808.ylog            = "no"
  plot808.rebin           = 10
  #plot808.xmin            = -0.5
  #plot808.xmax            = 99.5
  plot808.ymin            = 0
  plot808.ymax            = 20
  #plot808.lpos = "bottom-center"
  plot808.name            = "NN1stJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot808.addZUncBand     = zUncBand
  plot808.makeRatio       = makeRatio
  plot808.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NC1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NC1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot809 = Plot()
  ## inputs for stacked histograms
  plot809.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot809.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot809.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot809.keys            = keys
  plot809.xtit            = "Number of constituents 1^{st} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot809.ytit            = "Number of events"
  plot809.ylog            = "no"
  plot809.rebin           = 10
  #plot809.xmin            = -0.5
  #plot809.xmax            = 99.5
  plot809.ymin            = 0
  plot809.ymax            = 40
  #plot809.lpos = "bottom-center"
  plot809.name            = "NC1stJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot809.addZUncBand     = zUncBand
  plot809.makeRatio       = makeRatio
  plot809.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot810 = Plot()
  ## inputs for stacked histograms
  plot810.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot810.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot810.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot810.keys            = keys
  plot810.xtit            = "p_{T} 2^{nd} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot810.ytit            = "Number of events"
  plot810.ylog            = "no"
  plot810.rebin           = "var"
  #plot810.xmin            = 0
  #plot810.xmax            = 1000
  plot810.ymin            = 0
  plot810.ymax            = 20
  #plot810.lpos = "bottom-center"
  plot810.name            = "Pt2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot810.addZUncBand     = zUncBand
  plot810.makeRatio       = makeRatio
  plot810.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot810.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot811 = Plot()
  ## inputs for stacked histograms
  plot811.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot811.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot811.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot811.keys            = keys
  plot811.xtit            = "#eta 2^{nd} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot811.ytit            = "Number of events"
  plot811.ylog            = "no"
  plot811.rebin           = 5
  #plot811.xmin            = -5
  #plot811.xmax            = 5
  plot811.ymin            = 0
  plot811.ymax            = 15
  #plot811.lpos = "bottom-center"
  plot811.name            = "Eta2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot811.addZUncBand     = zUncBand
  plot811.makeRatio       = makeRatio
  plot811.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot812 = Plot()
  ## inputs for stacked histograms
  plot812.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot812.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot812.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot812.keys            = keys
  plot812.xtit            = "#phi 2^{nd} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot812.ytit            = "Number of events"
  plot812.ylog            = "no"
  plot812.rebin           = 10
  #plot812.xmin            = -3.1416
  #plot812.xmax            = 3.1416
  plot812.ymin            = 0
  plot812.ymax            = 20
  #plot812.lpos = "bottom-center"
  plot812.name            = "Phi2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot812.addZUncBand     = zUncBand
  plot812.makeRatio       = makeRatio
  plot812.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CHF2ndJet_highMePlusj_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_CHF2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot813 = Plot()
  ## inputs for stacked histograms
  plot813.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot813.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot813.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot813.keys            = keys
  plot813.xtit            = "Charged hadron energy fraction 2^{nd} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot813.ytit            = "Number of events"
  plot813.ylog            = "no"
  plot813.rebin           = 10
  #plot813.xmin            = 0
  #plot813.xmax            = 1
  plot813.ymin            = 0
  plot813.ymax            = 20
  #plot813.lpos = "bottom-center"
  plot813.name            = "CHF2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot813.addZUncBand     = zUncBand
  plot813.makeRatio       = makeRatio
  plot813.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NHF2ndJet_highMePlusj_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_NHF2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot814 = Plot()
  ## inputs for stacked histograms
  plot814.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot814.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot814.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot814.keys            = keys
  plot814.xtit            = "Neutral hadron energy fraction 2^{nd} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot814.ytit            = "Number of events"
  plot814.ylog            = "no"
  plot814.rebin           = 10
  #plot814.xmin            = 0
  #plot814.xmax            = 1
  plot814.ymin            = 0
  plot814.ymax            = 40
  #plot814.lpos = "bottom-center"
  plot814.name            = "NHF2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot814.addZUncBand     = zUncBand
  plot814.makeRatio       = makeRatio
  plot814.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CEF2ndJet_highMePlusj_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_CEF2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot815 = Plot()
  ## inputs for stacked histograms
  plot815.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot815.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot815.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot815.keys            = keys
  plot815.xtit            = "Charged EM energy fraction 2^{nd} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot815.ytit            = "Number of events"
  plot815.ylog            = "no"
  plot815.rebin           = 10
  #plot815.xmin            = 0
  #plot815.xmax            = 1
  plot815.ymin            = 0
  plot815.ymax            = 60
  #plot815.lpos = "bottom-center"
  plot815.name            = "CEF2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot815.addZUncBand     = zUncBand
  plot815.makeRatio       = makeRatio
  plot815.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NEF2ndJet_highMePlusj_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_NEF2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot816 = Plot()
  ## inputs for stacked histograms
  plot816.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot816.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot816.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot816.keys            = keys
  plot816.xtit            = "Neutral EM energy fraction 2^{nd} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot816.ytit            = "Number of events"
  plot816.ylog            = "no"
  plot816.rebin           = 10
  #plot816.xmin            = 0
  #plot816.xmax            = 1
  plot816.ymin            = 0
  plot816.ymax            = 20
  #plot816.lpos = "bottom-center"
  plot816.name            = "NEF2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot816.addZUncBand     = zUncBand
  plot816.makeRatio       = makeRatio
  plot816.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NCH2ndJet_highMePlusj_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_NCH2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot817 = Plot()
  ## inputs for stacked histograms
  plot817.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot817.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot817.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot817.keys            = keys
  plot817.xtit            = "Charged multiplicity 2^{nd} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot817.ytit            = "Number of events"
  plot817.ylog            = "no"
  plot817.rebin           = 10
  #plot817.xmin            = -0.5
  #plot817.xmax            = 99.5
  plot817.ymin            = 0
  plot817.ymax            = 40
  #plot817.lpos = "bottom-center"
  plot817.name            = "NCH2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot817.addZUncBand     = zUncBand
  plot817.makeRatio       = makeRatio
  plot817.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NN2ndJet_highMePlusj_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_NN2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot818 = Plot()
  ## inputs for stacked histograms
  plot818.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot818.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot818.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot818.keys            = keys
  plot818.xtit            = "Neutral multiplicity 2^{nd} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot818.ytit            = "Number of events"
  plot818.ylog            = "no"
  plot818.rebin           = 10
  #plot818.xmin            = -0.5
  #plot818.xmax            = 99.5
  plot818.ymin            = 0
  plot818.ymax            = 20
  #plot818.lpos = "bottom-center"
  plot818.name            = "NN2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot818.addZUncBand     = zUncBand
  plot818.makeRatio       = makeRatio
  plot818.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NC2ndJet_highMePlusj_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_NC2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot819 = Plot()
  ## inputs for stacked histograms
  plot819.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot819.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot819.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot819.keys            = keys
  plot819.xtit            = "Number of constituents 2^{nd} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot819.ytit            = "Number of events"
  plot819.ylog            = "no"
  plot819.rebin           = 10
  #plot819.xmin            = -0.5
  #plot819.xmax            = 99.5
  plot819.ymin            = 0
  plot819.ymax            = 30
  #plot819.lpos = "bottom-center"
  plot819.name            = "NC2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot819.addZUncBand     = zUncBand
  plot819.makeRatio       = makeRatio
  plot819.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_E1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_E1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot820 = Plot()
  ## inputs for stacked histograms
  plot820.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot820.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot820.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot820.keys            = keys
  plot820.xtit            = "Energy 1^{st} electron - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot820.ytit            = "Number of events"
  plot820.ylog            = "no"
  plot820.rebin           = "var"
  #plot820.xmin            = 0
  #plot820.xmax            = 1000
  plot820.ymin            = 0
  plot820.ymax            = 20
  #plot820.lpos = "bottom-center"
  plot820.name            = "E1stEle_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot820.addZUncBand     = zUncBand
  plot820.makeRatio       = makeRatio
  plot820.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot820.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot821 = Plot()
  ## inputs for stacked histograms
  plot821.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot821.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot821.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot821.keys            = keys
  plot821.xtit            = "p_{T} 1^{st} electron - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot821.ytit            = "Number of events"
  plot821.ylog            = "no"
  plot821.rebin           = "var"
  #plot821.xmin            = 0
  #plot821.xmax            = 1000
  plot821.ymin            = 0
  plot821.ymax            = 20
  #plot821.lpos = "bottom-center"
  plot821.name            = "Pt1stEle_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot821.addZUncBand     = zUncBand
  plot821.makeRatio       = makeRatio
  plot821.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot821.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot822 = Plot()
  ## inputs for stacked histograms
  plot822.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot822.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot822.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot822.keys            = keys
  plot822.xtit            = "#eta 1^{st} electron - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot822.ytit            = "Number of events"
  plot822.ylog            = "no"
  plot822.rebin           = 5
  #plot822.xmin            = -5
  #plot822.xmax            = 5
  plot822.ymin            = 0
  plot822.ymax            = 15
  #plot822.lpos = "bottom-center"
  plot822.name            = "Eta1stEle_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot822.addZUncBand     = zUncBand
  plot822.makeRatio       = makeRatio
  plot822.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot823 = Plot()
  ## inputs for stacked histograms
  plot823.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot823.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot823.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot823.keys            = keys
  plot823.xtit            = "#phi 1^{st} electron - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot823.ytit            = "Number of events"
  plot823.ylog            = "no"
  plot823.rebin           = 10
  #plot823.xmin            = -3.1416
  #plot823.xmax            = 3.1416
  plot823.ymin            = 0
  plot823.ymax            = 20
  #plot823.lpos = "bottom-center"
  plot823.name            = "Phi1stEle_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot823.addZUncBand     = zUncBand
  plot823.makeRatio       = makeRatio
  plot823.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Charge1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Charge1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot824 = Plot()
  ## inputs for stacked histograms
  plot824.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot824.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot824.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot824.keys            = keys
  plot824.xtit            = "Charge 1^{st} electron - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot824.ytit            = "Number of events"
  plot824.ylog            = "no"
  #plot824.rebin           = 1
  #plot824.xmin            = -1.001
  #plot824.xmax            = 1.001
  plot824.ymin            = 0
  plot824.ymax            = 50
  #plot824.lpos = "bottom-center"
  plot824.name            = "Charge1stEle_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot824.addZUncBand     = zUncBand
  plot824.makeRatio       = makeRatio
  plot824.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Conversion1stEle_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Conversion1stEle_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot825 = Plot()
  ## inputs for stacked histograms
  plot825.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot825.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot825.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot825.keys            = keys
  plot825.xtit            = "Conversion flag 1^{st} electron - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot825.ytit            = "Number of events"
  plot825.ylog            = "no"
  #plot825.rebin           = 1
  #plot825.xmin            = -0.5
  #plot825.xmax            = 1.5
  plot825.ymin            = 0
  plot825.ymax            = 60
  #plot825.lpos = "bottom-center"
  plot825.name            = "Conversion1stEle_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot825.addZUncBand     = zUncBand
  plot825.makeRatio       = makeRatio
  plot825.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMETEle_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiMETEle_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot826 = Plot()
  ## inputs for stacked histograms
  plot826.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot826.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot826.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot826.keys            = keys
  plot826.xtit            = "#Delta#phi(MET,1^{st} ele) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot826.ytit            = "Number of events"
  plot826.ylog            = "no"
  plot826.rebin           = 10
  #plot826.xmin            = 0
  #plot826.xmax            = 3.1416
  plot826.ymin            = 0
  plot826.ymax            = 20
  #plot826.lpos = "bottom-center"
  plot826.name            = "mDeltaPhiMETEle_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot826.addZUncBand     = zUncBand
  plot826.makeRatio       = makeRatio
  plot826.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiEle1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot827 = Plot()
  ## inputs for stacked histograms
  plot827.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot827.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot827.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot827.keys            = keys
  plot827.xtit            = "#Delta#phi(1^{st} ele,1^{st} jet) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot827.ytit            = "Number of events"
  plot827.ylog            = "no"
  plot827.rebin           = 10
  #plot827.xmin            = 0
  #plot827.xmax            = 3.1416
  plot827.ymin            = 0
  plot827.ymax            = 20
  #plot827.lpos = "bottom-center"
  plot827.name            = "mDeltaPhiEle1stJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot827.addZUncBand     = zUncBand
  plot827.makeRatio       = makeRatio
  plot827.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiEle2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot828 = Plot()
  ## inputs for stacked histograms
  plot828.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot828.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot828.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot828.keys            = keys
  plot828.xtit            = "#Delta#phi(1^{st} ele,2^{nd} jet) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot828.ytit            = "Number of events"
  plot828.ylog            = "no"
  plot828.rebin           = 10
  #plot828.xmin            = 0
  #plot828.xmax            = 3.1416
  plot828.ymin            = 0
  plot828.ymax            = 20
  #plot828.lpos = "bottom-center"
  plot828.name            = "mDeltaPhiEle2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot828.addZUncBand     = zUncBand
  plot828.makeRatio       = makeRatio
  plot828.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiMET2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot829 = Plot()
  ## inputs for stacked histograms
  plot829.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot829.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot829.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot829.keys            = keys
  plot829.xtit            = "#Delta#phi(MET,2^{nd} jet) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot829.ytit            = "Number of events"
  plot829.ylog            = "no"
  plot829.rebin           = 10
  #plot829.xmin            = 0
  #plot829.xmax            = 3.1416
  plot829.ymin            = 0
  plot829.ymax            = 20
  #plot829.lpos = "bottom-center"
  plot829.name            = "mDeltaPhiMET2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot829.addZUncBand     = zUncBand
  plot829.makeRatio       = makeRatio
  plot829.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Ptenu_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Ptenu_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot830 = Plot()
  ## inputs for stacked histograms
  plot830.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot830.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot830.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot830.keys            = keys
  plot830.xtit            = "p_{T}(e#nu) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot830.ytit            = "Number of events"
  plot830.ylog            = "no"
  plot830.rebin           = "var"
  #plot830.xmin            = 0
  #plot830.xmax            = 2000
  plot830.ymin            = 0
  plot830.ymax            = 20
  #plot830.lpos = "bottom-center"
  plot830.name            = "Ptenu_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot830.addZUncBand     = zUncBand
  plot830.makeRatio       = makeRatio
  plot830.xbins           = [0,40,80,120,160,200,300,600]
  plot830.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_1stJet_PTOverPTPlusMET_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_1stJet_PTOverPTPlusMET_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot831 = Plot()
  ## inputs for stacked histograms
  plot831.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot831.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot831.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot831.keys            = keys
  plot831.xtit            = "Pt1stJet/(Pt1stJet+MET) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot831.ytit            = "Number of events"
  plot831.ylog            = "no"
  plot831.rebin           = 5
  #plot831.xmin            = 0
  #plot831.xmax            = 1
  plot831.ymin            = 0
  plot831.ymax            = 20
  #plot831.lpos = "bottom-center"
  plot831.name            = "1stJet_PTOverPTPlusMET_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot831.addZUncBand     = zUncBand
  plot831.makeRatio       = makeRatio
  plot831.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MET_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_MET_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot832 = Plot()
  ## inputs for stacked histograms
  plot832.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot832.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot832.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot832.keys            = keys
  plot832.xtit            = "pfMET - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot832.ytit            = "Number of events"
  plot832.ylog            = "no"
  plot832.rebin           = "var"
  #plot832.xmin            = 0
  #plot832.xmax            = 1000
  plot832.ymin            = 0
  plot832.ymax            = 20
  #plot832.lpos = "bottom-center"
  plot832.name            = "MET_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot832.addZUncBand     = zUncBand
  plot832.makeRatio       = makeRatio
  plot832.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot832.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Njet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot833 = Plot()
  ## inputs for stacked histograms
  plot833.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot833.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot833.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot833.keys            = keys
  plot833.xtit            = "Number of jets - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot833.ytit            = "Number of events"
  plot833.ylog            = "no"
  plot833.rebin           = 1
  plot833.xmin            = -0.5
  plot833.xmax            = 11.5
  plot833.ymin            = 0
  plot833.ymax            = 25
  #plot833.lpos = "bottom-center"
  plot833.name            = "Njet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot833.addZUncBand     = zUncBand
  plot833.makeRatio       = makeRatio
  plot833.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NjetTCHELBTag_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot834 = Plot()
  ## inputs for stacked histograms
  plot834.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot834.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot834.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot834.keys            = keys
  plot834.xtit            = "Number of b-tagged jets (TCHEL) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot834.ytit            = "Number of events"
  plot834.ylog            = "no"
  plot834.rebin           = 1
  plot834.xmin            = -0.5
  plot834.xmax            = 11.5
  plot834.ymin            = 0
  plot834.ymax            = 25
  #plot834.lpos = "bottom-center"
  plot834.name            = "NjetTCHELBTag_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot834.addZUncBand     = zUncBand
  plot834.makeRatio       = makeRatio
  plot834.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_minDRej_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_minDRej_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot835 = Plot()
  ## inputs for stacked histograms
  plot835.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot835.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot835.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot835.keys            = keys
  plot835.xtit            = "min#DeltaR(e,j) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot835.ytit            = "Number of events"
  plot835.ylog            = "no"
  plot835.rebin           = 10
  #plot835.xmin            = 0
  #plot835.xmax            = 10
  plot835.ymin            = 0
  plot835.ymax            = 25
  #plot835.lpos = "bottom-center"
  plot835.name            = "minDRej_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot835.addZUncBand     = zUncBand
  plot835.makeRatio       = makeRatio
  plot835.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_maxDRej_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_maxDRej_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot836 = Plot()
  ## inputs for stacked histograms
  plot836.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot836.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot836.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot836.keys            = keys
  plot836.xtit            = "max#DeltaR(e,j) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot836.ytit            = "Number of events"
  plot836.ylog            = "no"
  plot836.rebin           = 10
  #plot836.xmin            = 0
  #plot836.xmax            = 10
  plot836.ymin            = 0
  plot836.ymax            = 30
  #plot836.lpos = "bottom-center"
  plot836.name            = "maxDRej_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot836.addZUncBand     = zUncBand
  plot836.makeRatio       = makeRatio
  plot836.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DRjets_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_DRjets_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot837 = Plot()
  ## inputs for stacked histograms
  plot837.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot837.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot837.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot837.keys            = keys
  plot837.xtit            = "#DeltaR(j1,j2) - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot837.ytit            = "Number of events"
  plot837.ylog            = "no"
  plot837.rebin           = 10
  #plot837.xmin            = 0
  #plot837.xmax            = 10
  plot837.ymin            = 0
  plot837.ymax            = 30
  #plot837.lpos = "bottom-center"
  plot837.name            = "DRjets_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot837.addZUncBand     = zUncBand
  plot837.makeRatio       = makeRatio
  plot837.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Mass1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Mass1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot838 = Plot()
  ## inputs for stacked histograms
  plot838.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot838.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot838.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot838.keys            = keys
  plot838.xtit            = "Mass 1^{st} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot838.ytit            = "Number of events"
  plot838.ylog            = "no"
  plot838.rebin           = 5
  #plot838.xmin            = -0.5
  #plot838.xmax            = 99.5
  plot838.ymin            = 0
  plot838.ymax            = 20
  #plot838.lpos = "bottom-center"
  plot838.name            = "Mass1stJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot838.addZUncBand     = zUncBand
  plot838.makeRatio       = makeRatio
  plot838.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Mass2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Mass2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot842 = Plot()
  ## inputs for stacked histograms
  plot842.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot842.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot842.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot842.keys            = keys
  plot842.xtit            = "Mass 2^{nd} jet - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot842.ytit            = "Number of events"
  plot842.ylog            = "no"
  plot842.rebin           = 5
  #plot842.xmin            = -0.5
  #plot842.xmax            = 99.5
  plot842.ymin            = 0
  plot842.ymax            = 20
  #plot842.lpos = "bottom-center"
  plot842.name            = "Mass2ndJet_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot842.addZUncBand     = zUncBand
  plot842.makeRatio       = makeRatio
  plot842.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CaloMET_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_CaloMET_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot843 = Plot()
  ## inputs for stacked histograms
  plot843.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot843.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot843.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot843.keys            = keys
  plot843.xtit            = "caloMET - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot843.ytit            = "Number of events"
  plot843.ylog            = "no"
  plot843.rebin           = "var"
  #plot843.xmin            = 0
  #plot843.xmax            = 1000
  plot843.ymin            = 0
  plot843.ymax            = 20
  #plot843.lpos = "bottom-center"
  plot843.name            = "CaloMET_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot843.addZUncBand     = zUncBand
  plot843.makeRatio       = makeRatio
  plot843.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot843.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_TCMET_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_TCMET_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot844 = Plot()
  ## inputs for stacked histograms
  plot844.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot844.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot844.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot844.keys            = keys
  plot844.xtit            = "tcMET - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot844.ytit            = "Number of events"
  plot844.ylog            = "no"
  plot844.rebin           = "var"
  #plot844.xmin            = 0
  #plot844.xmax            = 1000
  plot844.ymin            = 0
  plot844.ymax            = 20
  #plot844.lpos = "bottom-center"
  plot844.name            = "TCMET_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot844.addZUncBand     = zUncBand
  plot844.makeRatio       = makeRatio
  plot844.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot844.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

  #--- h1_MET_over_CaloMET_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_MET_over_CaloMET_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot845 = Plot()
  ## inputs for stacked histograms
  plot845.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot845.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot845.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot845.keys            = keys
  plot845.xtit            = "pfMET/CaloMET - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot845.ytit            = "Number of events"
  plot845.ylog            = "no"
  plot845.rebin           = 5
  #plot845.xmin            = 0
  #plot845.xmax            = 10
  plot845.ymin            = 0
  plot845.ymax            = 20
  #plot845.lpos = "bottom-center"
  plot845.name            = "MET_over_CaloMET_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot845.addZUncBand     = zUncBand
  plot845.makeRatio       = makeRatio
  plot845.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MTenu_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_MTenu_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"

  plot846 = Plot()
  ## inputs for stacked histograms
  plot846.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot846.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot846.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot846.keys            = keys
  plot846.xtit            = "MTenu - highMePlusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot846.ytit            = "Number of events"
  plot846.ylog            = "no"
  plot846.rebin           = "var"
  #plot846.xmin            = 0
  #plot846.xmax            = 1000
  plot846.ymin            = 0
  plot846.ymax            = 20
  #plot846.lpos = "bottom-center"
  plot846.name            = "MTenu_presel_highMePlusj_mDeltaPhiMET1stJet_gt_2.5"
  plot846.addZUncBand     = zUncBand
  plot846.makeRatio       = makeRatio
  plot846.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot846.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)



  # mDeltaPhiMET1stJet>2.5, e- only

  #--- h1_Pt1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot900 = Plot()
  ## inputs for stacked histograms
  plot900.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot900.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot900.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot900.keys            = keys
  plot900.xtit            = "p_{T} 1^{st} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot900.ytit            = "Number of events"
  plot900.ylog            = "no"
  plot900.rebin           = "var"
  #plot900.xmin            = 0
  #plot900.xmax            = 1000
  plot900.ymin            = 0
  plot900.ymax            = 15
  #plot900.lpos = "bottom-center"
  plot900.name            = "Pt1stJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot900.addZUncBand     = zUncBand
  plot900.makeRatio       = makeRatio
  plot900.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot900.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot901 = Plot()
  ## inputs for stacked histograms
  plot901.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot901.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot901.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot901.keys            = keys
  plot901.xtit            = "#eta 1^{st} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot901.ytit            = "Number of events"
  plot901.ylog            = "no"
  plot901.rebin           = 5
  #plot901.xmin            = -5
  #plot901.xmax            = 5
  plot901.ymin            = 0
  plot901.ymax            = 15
  #plot901.lpos = "bottom-center"
  plot901.name            = "Eta1stJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot901.addZUncBand     = zUncBand
  plot901.makeRatio       = makeRatio
  plot901.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot902 = Plot()
  ## inputs for stacked histograms
  plot902.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot902.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot902.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot902.keys            = keys
  plot902.xtit            = "#phi 1^{st} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot902.ytit            = "Number of events"
  plot902.ylog            = "no"
  plot902.rebin           = 10
  #plot902.xmin            = -3.1416
  #plot902.xmax            = 3.1416
  plot902.ymin            = 0
  plot902.ymax            = 20
  #plot902.lpos = "bottom-center"
  plot902.name            = "Phi1stJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot902.addZUncBand     = zUncBand
  plot902.makeRatio       = makeRatio
  plot902.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CHF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_CHF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot903 = Plot()
  ## inputs for stacked histograms
  plot903.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot903.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot903.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot903.keys            = keys
  plot903.xtit            = "Charged hadron energy fraction 1^{st} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot903.ytit            = "Number of events"
  plot903.ylog            = "no"
  plot903.rebin           = 10
  #plot903.xmin            = 0
  #plot903.xmax            = 1
  plot903.ymin            = 0
  plot903.ymax            = 20
  #plot903.lpos = "bottom-center"
  plot903.name            = "CHF1stJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot903.addZUncBand     = zUncBand
  plot903.makeRatio       = makeRatio
  plot903.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NHF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NHF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot904 = Plot()
  ## inputs for stacked histograms
  plot904.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot904.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot904.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot904.keys            = keys
  plot904.xtit            = "Neutral hadron energy fraction 1^{st} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot904.ytit            = "Number of events"
  plot904.ylog            = "no"
  plot904.rebin           = 10
  #plot904.xmin            = 0
  #plot904.xmax            = 1
  plot904.ymin            = 0
  plot904.ymax            = 40
  #plot904.lpos = "bottom-center"
  plot904.name            = "NHF1stJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot904.addZUncBand     = zUncBand
  plot904.makeRatio       = makeRatio
  plot904.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CEF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_CEF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot905 = Plot()
  ## inputs for stacked histograms
  plot905.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot905.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot905.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot905.keys            = keys
  plot905.xtit            = "Charged EM energy fraction 1^{st} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot905.ytit            = "Number of events"
  plot905.ylog            = "no"
  plot905.rebin           = 10
  #plot905.xmin            = 0
  #plot905.xmax            = 1
  plot905.ymin            = 0
  plot905.ymax            = 60
  #plot905.lpos = "bottom-center"
  plot905.name            = "CEF1stJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot905.addZUncBand     = zUncBand
  plot905.makeRatio       = makeRatio
  plot905.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NEF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NEF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot906 = Plot()
  ## inputs for stacked histograms
  plot906.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot906.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot906.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot906.keys            = keys
  plot906.xtit            = "Neutral EM energy fraction 1^{st} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot906.ytit            = "Number of events"
  plot906.ylog            = "no"
  plot906.rebin           = 10
  #plot906.xmin            = 0
  #plot906.xmax            = 1
  plot906.ymin            = 0
  plot906.ymax            = 20
  #plot906.lpos = "bottom-center"
  plot906.name            = "NEF1stJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot906.addZUncBand     = zUncBand
  plot906.makeRatio       = makeRatio
  plot906.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NCH1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NCH1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot907 = Plot()
  ## inputs for stacked histograms
  plot907.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot907.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot907.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot907.keys            = keys
  plot907.xtit            = "Charged multiplicity 1^{st} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot907.ytit            = "Number of events"
  plot907.ylog            = "no"
  plot907.rebin           = 10
  #plot907.xmin            = -0.5
  #plot907.xmax            = 99.5
  plot907.ymin            = 0
  plot907.ymax            = 40
  #plot907.lpos = "bottom-center"
  plot907.name            = "NCH1stJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot907.addZUncBand     = zUncBand
  plot907.makeRatio       = makeRatio
  plot907.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NN1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NN1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot908 = Plot()
  ## inputs for stacked histograms
  plot908.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot908.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot908.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot908.keys            = keys
  plot908.xtit            = "Neutral multiplicity 1^{st} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot908.ytit            = "Number of events"
  plot908.ylog            = "no"
  plot908.rebin           = 10
  #plot908.xmin            = -0.5
  #plot908.xmax            = 99.5
  plot908.ymin            = 0
  plot908.ymax            = 20
  #plot908.lpos = "bottom-center"
  plot908.name            = "NN1stJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot908.addZUncBand     = zUncBand
  plot908.makeRatio       = makeRatio
  plot908.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NC1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NC1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot909 = Plot()
  ## inputs for stacked histograms
  plot909.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot909.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot909.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot909.keys            = keys
  plot909.xtit            = "Number of constituents 1^{st} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot909.ytit            = "Number of events"
  plot909.ylog            = "no"
  plot909.rebin           = 10
  #plot909.xmin            = -0.5
  #plot909.xmax            = 99.5
  plot909.ymin            = 0
  plot909.ymax            = 40
  #plot909.lpos = "bottom-center"
  plot909.name            = "NC1stJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot909.addZUncBand     = zUncBand
  plot909.makeRatio       = makeRatio
  plot909.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot910 = Plot()
  ## inputs for stacked histograms
  plot910.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot910.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot910.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot910.keys            = keys
  plot910.xtit            = "p_{T} 2^{nd} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot910.ytit            = "Number of events"
  plot910.ylog            = "no"
  plot910.rebin           = "var"
  #plot910.xmin            = 0
  #plot910.xmax            = 1000
  plot910.ymin            = 0
  plot910.ymax            = 20
  #plot910.lpos = "bottom-center"
  plot910.name            = "Pt2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot910.addZUncBand     = zUncBand
  plot910.makeRatio       = makeRatio
  plot910.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot910.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot911 = Plot()
  ## inputs for stacked histograms
  plot911.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot911.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot911.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot911.keys            = keys
  plot911.xtit            = "#eta 2^{nd} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot911.ytit            = "Number of events"
  plot911.ylog            = "no"
  plot911.rebin           = 5
  #plot911.xmin            = -5
  #plot911.xmax            = 5
  plot911.ymin            = 0
  plot911.ymax            = 15
  #plot911.lpos = "bottom-center"
  plot911.name            = "Eta2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot911.addZUncBand     = zUncBand
  plot911.makeRatio       = makeRatio
  plot911.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot912 = Plot()
  ## inputs for stacked histograms
  plot912.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot912.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot912.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot912.keys            = keys
  plot912.xtit            = "#phi 2^{nd} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot912.ytit            = "Number of events"
  plot912.ylog            = "no"
  plot912.rebin           = 10
  #plot912.xmin            = -3.1416
  #plot912.xmax            = 3.1416
  plot912.ymin            = 0
  plot912.ymax            = 20
  #plot912.lpos = "bottom-center"
  plot912.name            = "Phi2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot912.addZUncBand     = zUncBand
  plot912.makeRatio       = makeRatio
  plot912.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CHF2ndJet_highMeMinusj_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_CHF2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot913 = Plot()
  ## inputs for stacked histograms
  plot913.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot913.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot913.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot913.keys            = keys
  plot913.xtit            = "Charged hadron energy fraction 2^{nd} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot913.ytit            = "Number of events"
  plot913.ylog            = "no"
  plot913.rebin           = 10
  #plot913.xmin            = 0
  #plot913.xmax            = 1
  plot913.ymin            = 0
  plot913.ymax            = 20
  #plot913.lpos = "bottom-center"
  plot913.name            = "CHF2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot913.addZUncBand     = zUncBand
  plot913.makeRatio       = makeRatio
  plot913.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NHF2ndJet_highMeMinusj_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_NHF2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot914 = Plot()
  ## inputs for stacked histograms
  plot914.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot914.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot914.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot914.keys            = keys
  plot914.xtit            = "Neutral hadron energy fraction 2^{nd} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot914.ytit            = "Number of events"
  plot914.ylog            = "no"
  plot914.rebin           = 10
  #plot914.xmin            = 0
  #plot914.xmax            = 1
  plot914.ymin            = 0
  plot914.ymax            = 40
  #plot914.lpos = "bottom-center"
  plot914.name            = "NHF2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot914.addZUncBand     = zUncBand
  plot914.makeRatio       = makeRatio
  plot914.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CEF2ndJet_highMeMinusj_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_CEF2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot915 = Plot()
  ## inputs for stacked histograms
  plot915.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot915.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot915.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot915.keys            = keys
  plot915.xtit            = "Charged EM energy fraction 2^{nd} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot915.ytit            = "Number of events"
  plot915.ylog            = "no"
  plot915.rebin           = 10
  #plot915.xmin            = 0
  #plot915.xmax            = 1
  plot915.ymin            = 0
  plot915.ymax            = 60
  #plot915.lpos = "bottom-center"
  plot915.name            = "CEF2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot915.addZUncBand     = zUncBand
  plot915.makeRatio       = makeRatio
  plot915.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NEF2ndJet_highMeMinusj_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_NEF2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot916 = Plot()
  ## inputs for stacked histograms
  plot916.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot916.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot916.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot916.keys            = keys
  plot916.xtit            = "Neutral EM energy fraction 2^{nd} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot916.ytit            = "Number of events"
  plot916.ylog            = "no"
  plot916.rebin           = 10
  #plot916.xmin            = 0
  #plot916.xmax            = 1
  plot916.ymin            = 0
  plot916.ymax            = 20
  #plot916.lpos = "bottom-center"
  plot916.name            = "NEF2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot916.addZUncBand     = zUncBand
  plot916.makeRatio       = makeRatio
  plot916.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NCH2ndJet_highMeMinusj_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_NCH2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot917 = Plot()
  ## inputs for stacked histograms
  plot917.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot917.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot917.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot917.keys            = keys
  plot917.xtit            = "Charged multiplicity 2^{nd} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot917.ytit            = "Number of events"
  plot917.ylog            = "no"
  plot917.rebin           = 10
  #plot917.xmin            = -0.5
  #plot917.xmax            = 99.5
  plot917.ymin            = 0
  plot917.ymax            = 40
  #plot917.lpos = "bottom-center"
  plot917.name            = "NCH2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot917.addZUncBand     = zUncBand
  plot917.makeRatio       = makeRatio
  plot917.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NN2ndJet_highMeMinusj_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_NN2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot918 = Plot()
  ## inputs for stacked histograms
  plot918.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot918.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot918.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot918.keys            = keys
  plot918.xtit            = "Neutral multiplicity 2^{nd} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot918.ytit            = "Number of events"
  plot918.ylog            = "no"
  plot918.rebin           = 10
  #plot918.xmin            = -0.5
  #plot918.xmax            = 99.5
  plot918.ymin            = 0
  plot918.ymax            = 20
  #plot918.lpos = "bottom-center"
  plot918.name            = "NN2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot918.addZUncBand     = zUncBand
  plot918.makeRatio       = makeRatio
  plot918.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NC2ndJet_highMeMinusj_mDeltaPhiMET2ndJet_gt_2.5 ---
  variableName = "h1_NC2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot919 = Plot()
  ## inputs for stacked histograms
  plot919.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot919.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot919.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot919.keys            = keys
  plot919.xtit            = "Number of constituents 2^{nd} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot919.ytit            = "Number of events"
  plot919.ylog            = "no"
  plot919.rebin           = 10
  #plot919.xmin            = -0.5
  #plot919.xmax            = 99.5
  plot919.ymin            = 0
  plot919.ymax            = 30
  #plot919.lpos = "bottom-center"
  plot919.name            = "NC2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot919.addZUncBand     = zUncBand
  plot919.makeRatio       = makeRatio
  plot919.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_E1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_E1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot920 = Plot()
  ## inputs for stacked histograms
  plot920.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot920.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot920.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot920.keys            = keys
  plot920.xtit            = "Energy 1^{st} electron - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot920.ytit            = "Number of events"
  plot920.ylog            = "no"
  plot920.rebin           = "var"
  #plot920.xmin            = 0
  #plot920.xmax            = 1000
  plot920.ymin            = 0
  plot920.ymax            = 20
  #plot920.lpos = "bottom-center"
  plot920.name            = "E1stEle_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot920.addZUncBand     = zUncBand
  plot920.makeRatio       = makeRatio
  plot920.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot920.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Pt1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot921 = Plot()
  ## inputs for stacked histograms
  plot921.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot921.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot921.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot921.keys            = keys
  plot921.xtit            = "p_{T} 1^{st} electron - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot921.ytit            = "Number of events"
  plot921.ylog            = "no"
  plot921.rebin           = "var"
  #plot921.xmin            = 0
  #plot921.xmax            = 1000
  plot921.ymin            = 0
  plot921.ymax            = 20
  #plot921.lpos = "bottom-center"
  plot921.name            = "Pt1stEle_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot921.addZUncBand     = zUncBand
  plot921.makeRatio       = makeRatio
  plot921.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot921.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Eta1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot922 = Plot()
  ## inputs for stacked histograms
  plot922.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot922.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot922.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot922.keys            = keys
  plot922.xtit            = "#eta 1^{st} electron - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot922.ytit            = "Number of events"
  plot922.ylog            = "no"
  plot922.rebin           = 5
  #plot922.xmin            = -5
  #plot922.xmax            = 5
  plot922.ymin            = 0
  plot922.ymax            = 15
  #plot922.lpos = "bottom-center"
  plot922.name            = "Eta1stEle_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot922.addZUncBand     = zUncBand
  plot922.makeRatio       = makeRatio
  plot922.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Phi1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot923 = Plot()
  ## inputs for stacked histograms
  plot923.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot923.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot923.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot923.keys            = keys
  plot923.xtit            = "#phi 1^{st} electron - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot923.ytit            = "Number of events"
  plot923.ylog            = "no"
  plot923.rebin           = 10
  #plot923.xmin            = -3.1416
  #plot923.xmax            = 3.1416
  plot923.ymin            = 0
  plot923.ymax            = 20
  #plot923.lpos = "bottom-center"
  plot923.name            = "Phi1stEle_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot923.addZUncBand     = zUncBand
  plot923.makeRatio       = makeRatio
  plot923.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Charge1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Charge1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot924 = Plot()
  ## inputs for stacked histograms
  plot924.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot924.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot924.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot924.keys            = keys
  plot924.xtit            = "Charge 1^{st} electron - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot924.ytit            = "Number of events"
  plot924.ylog            = "no"
  #plot924.rebin           = 1
  #plot924.xmin            = -1.001
  #plot924.xmax            = 1.001
  plot924.ymin            = 0
  plot924.ymax            = 50
  #plot924.lpos = "bottom-center"
  plot924.name            = "Charge1stEle_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot924.addZUncBand     = zUncBand
  plot924.makeRatio       = makeRatio
  plot924.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Conversion1stEle_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Conversion1stEle_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot925 = Plot()
  ## inputs for stacked histograms
  plot925.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot925.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot925.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot925.keys            = keys
  plot925.xtit            = "Conversion flag 1^{st} electron - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot925.ytit            = "Number of events"
  plot925.ylog            = "no"
  #plot925.rebin           = 1
  #plot925.xmin            = -0.5
  #plot925.xmax            = 1.5
  plot925.ymin            = 0
  plot925.ymax            = 60
  #plot925.lpos = "bottom-center"
  plot925.name            = "Conversion1stEle_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot925.addZUncBand     = zUncBand
  plot925.makeRatio       = makeRatio
  plot925.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMETEle_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiMETEle_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot926 = Plot()
  ## inputs for stacked histograms
  plot926.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot926.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot926.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot926.keys            = keys
  plot926.xtit            = "#Delta#phi(MET,1^{st} ele) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot926.ytit            = "Number of events"
  plot926.ylog            = "no"
  plot926.rebin           = 10
  #plot926.xmin            = 0
  #plot926.xmax            = 3.1416
  plot926.ymin            = 0
  plot926.ymax            = 20
  #plot926.lpos = "bottom-center"
  plot926.name            = "mDeltaPhiMETEle_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot926.addZUncBand     = zUncBand
  plot926.makeRatio       = makeRatio
  plot926.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiEle1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot927 = Plot()
  ## inputs for stacked histograms
  plot927.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot927.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot927.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot927.keys            = keys
  plot927.xtit            = "#Delta#phi(1^{st} ele,1^{st} jet) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot927.ytit            = "Number of events"
  plot927.ylog            = "no"
  plot927.rebin           = 10
  #plot927.xmin            = 0
  #plot927.xmax            = 3.1416
  plot927.ymin            = 0
  plot927.ymax            = 20
  #plot927.lpos = "bottom-center"
  plot927.name            = "mDeltaPhiEle1stJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot927.addZUncBand     = zUncBand
  plot927.makeRatio       = makeRatio
  plot927.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiEle2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot928 = Plot()
  ## inputs for stacked histograms
  plot928.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot928.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot928.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot928.keys            = keys
  plot928.xtit            = "#Delta#phi(1^{st} ele,2^{nd} jet) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot928.ytit            = "Number of events"
  plot928.ylog            = "no"
  plot928.rebin           = 10
  #plot928.xmin            = 0
  #plot928.xmax            = 3.1416
  plot928.ymin            = 0
  plot928.ymax            = 20
  #plot928.lpos = "bottom-center"
  plot928.name            = "mDeltaPhiEle2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot928.addZUncBand     = zUncBand
  plot928.makeRatio       = makeRatio
  plot928.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_mDeltaPhiMET2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot929 = Plot()
  ## inputs for stacked histograms
  plot929.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot929.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot929.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot929.keys            = keys
  plot929.xtit            = "#Delta#phi(MET,2^{nd} jet) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot929.ytit            = "Number of events"
  plot929.ylog            = "no"
  plot929.rebin           = 10
  #plot929.xmin            = 0
  #plot929.xmax            = 3.1416
  plot929.ymin            = 0
  plot929.ymax            = 20
  #plot929.lpos = "bottom-center"
  plot929.name            = "mDeltaPhiMET2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot929.addZUncBand     = zUncBand
  plot929.makeRatio       = makeRatio
  plot929.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Ptenu_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Ptenu_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot930 = Plot()
  ## inputs for stacked histograms
  plot930.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot930.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot930.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot930.keys            = keys
  plot930.xtit            = "p_{T}(e#nu) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot930.ytit            = "Number of events"
  plot930.ylog            = "no"
  plot930.rebin           = "var"
  #plot930.xmin            = 0
  #plot930.xmax            = 2000
  plot930.ymin            = 0
  plot930.ymax            = 20
  #plot930.lpos = "bottom-center"
  plot930.name            = "Ptenu_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot930.addZUncBand     = zUncBand
  plot930.makeRatio       = makeRatio
  plot930.xbins           = [0,40,80,120,160,200,300,600]
  plot930.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_1stJet_PTOverPTPlusMET_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_1stJet_PTOverPTPlusMET_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot931 = Plot()
  ## inputs for stacked histograms
  plot931.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot931.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot931.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot931.keys            = keys
  plot931.xtit            = "Pt1stJet/(Pt1stJet+MET) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot931.ytit            = "Number of events"
  plot931.ylog            = "no"
  plot931.rebin           = 5
  #plot931.xmin            = 0
  #plot931.xmax            = 1
  plot931.ymin            = 0
  plot931.ymax            = 20
  #plot931.lpos = "bottom-center"
  plot931.name            = "1stJet_PTOverPTPlusMET_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot931.addZUncBand     = zUncBand
  plot931.makeRatio       = makeRatio
  plot931.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MET_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_MET_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot932 = Plot()
  ## inputs for stacked histograms
  plot932.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot932.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot932.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot932.keys            = keys
  plot932.xtit            = "pfMET - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot932.ytit            = "Number of events"
  plot932.ylog            = "no"
  plot932.rebin           = "var"
  #plot932.xmin            = 0
  #plot932.xmax            = 1000
  plot932.ymin            = 0
  plot932.ymax            = 20
  #plot932.lpos = "bottom-center"
  plot932.name            = "MET_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot932.addZUncBand     = zUncBand
  plot932.makeRatio       = makeRatio
  plot932.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot932.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Njet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot933 = Plot()
  ## inputs for stacked histograms
  plot933.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot933.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot933.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot933.keys            = keys
  plot933.xtit            = "Number of jets - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot933.ytit            = "Number of events"
  plot933.ylog            = "no"
  plot933.rebin           = 1
  plot933.xmin            = -0.5
  plot933.xmax            = 11.5
  plot933.ymin            = 0
  plot933.ymax            = 25
  #plot933.lpos = "bottom-center"
  plot933.name            = "Njet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot933.addZUncBand     = zUncBand
  plot933.makeRatio       = makeRatio
  plot933.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_NjetTCHELBTag_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot934 = Plot()
  ## inputs for stacked histograms
  plot934.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot934.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot934.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot934.keys            = keys
  plot934.xtit            = "Number of b-tagged jets (TCHEL) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot934.ytit            = "Number of events"
  plot934.ylog            = "no"
  plot934.rebin           = 1
  plot934.xmin            = -0.5
  plot934.xmax            = 11.5
  plot934.ymin            = 0
  plot934.ymax            = 25
  #plot934.lpos = "bottom-center"
  plot934.name            = "NjetTCHELBTag_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot934.addZUncBand     = zUncBand
  plot934.makeRatio       = makeRatio
  plot934.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_minDRej_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_minDRej_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot935 = Plot()
  ## inputs for stacked histograms
  plot935.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot935.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot935.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot935.keys            = keys
  plot935.xtit            = "min#DeltaR(e,j) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot935.ytit            = "Number of events"
  plot935.ylog            = "no"
  plot935.rebin           = 10
  #plot935.xmin            = 0
  #plot935.xmax            = 10
  plot935.ymin            = 0
  plot935.ymax            = 25
  #plot935.lpos = "bottom-center"
  plot935.name            = "minDRej_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot935.addZUncBand     = zUncBand
  plot935.makeRatio       = makeRatio
  plot935.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_maxDRej_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_maxDRej_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot936 = Plot()
  ## inputs for stacked histograms
  plot936.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot936.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot936.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot936.keys            = keys
  plot936.xtit            = "max#DeltaR(e,j) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot936.ytit            = "Number of events"
  plot936.ylog            = "no"
  plot936.rebin           = 10
  #plot936.xmin            = 0
  #plot936.xmax            = 10
  plot936.ymin            = 0
  plot936.ymax            = 30
  #plot936.lpos = "bottom-center"
  plot936.name            = "maxDRej_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot936.addZUncBand     = zUncBand
  plot936.makeRatio       = makeRatio
  plot936.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DRjets_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_DRjets_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot937 = Plot()
  ## inputs for stacked histograms
  plot937.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot937.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot937.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot937.keys            = keys
  plot937.xtit            = "#DeltaR(j1,j2) - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot937.ytit            = "Number of events"
  plot937.ylog            = "no"
  plot937.rebin           = 10
  #plot937.xmin            = 0
  #plot937.xmax            = 10
  plot937.ymin            = 0
  plot937.ymax            = 30
  #plot937.lpos = "bottom-center"
  plot937.name            = "DRjets_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot937.addZUncBand     = zUncBand
  plot937.makeRatio       = makeRatio
  plot937.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)



  #--- h1_Mass1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Mass1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot938 = Plot()
  ## inputs for stacked histograms
  plot938.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot938.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot938.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot938.keys            = keys
  plot938.xtit            = "Mass 1^{st} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot938.ytit            = "Number of events"
  plot938.ylog            = "no"
  plot938.rebin           = 5
  #plot938.xmin            = -0.5
  #plot938.xmax            = 99.5
  plot938.ymin            = 0
  plot938.ymax            = 20
  #plot938.lpos = "bottom-center"
  plot938.name            = "Mass1stJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot938.addZUncBand     = zUncBand
  plot938.makeRatio       = makeRatio
  plot938.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Mass2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_Mass2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot942 = Plot()
  ## inputs for stacked histograms
  plot942.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot942.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot942.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot942.keys            = keys
  plot942.xtit            = "Mass 2^{nd} jet - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot942.ytit            = "Number of events"
  plot942.ylog            = "no"
  plot942.rebin           = 5
  #plot942.xmin            = -0.5
  #plot942.xmax            = 99.5
  plot942.ymin            = 0
  plot942.ymax            = 20
  #plot942.lpos = "bottom-center"
  plot942.name            = "Mass2ndJet_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot942.addZUncBand     = zUncBand
  plot942.makeRatio       = makeRatio
  plot942.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CaloMET_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_CaloMET_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot943 = Plot()
  ## inputs for stacked histograms
  plot943.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot943.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot943.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot943.keys            = keys
  plot943.xtit            = "caloMET - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot943.ytit            = "Number of events"
  plot943.ylog            = "no"
  plot943.rebin           = "var"
  #plot943.xmin            = 0
  #plot943.xmax            = 1000
  plot943.ymin            = 0
  plot943.ymax            = 20
  #plot943.lpos = "bottom-center"
  plot943.name            = "CaloMET_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot943.addZUncBand     = zUncBand
  plot943.makeRatio       = makeRatio
  plot943.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot943.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_TCMET_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_TCMET_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot944 = Plot()
  ## inputs for stacked histograms
  plot944.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot944.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot944.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot944.keys            = keys
  plot944.xtit            = "tcMET - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot944.ytit            = "Number of events"
  plot944.ylog            = "no"
  plot944.rebin           = "var"
  #plot944.xmin            = 0
  #plot944.xmax            = 1000
  plot944.ymin            = 0
  plot944.ymax            = 20
  #plot944.lpos = "bottom-center"
  plot944.name            = "TCMET_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot944.addZUncBand     = zUncBand
  plot944.makeRatio       = makeRatio
  plot944.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot944.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

  #--- h1_MET_over_CaloMET_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_MET_over_CaloMET_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot945 = Plot()
  ## inputs for stacked histograms
  plot945.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot945.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot945.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot945.keys            = keys
  plot945.xtit            = "pfMET/CaloMET - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot945.ytit            = "Number of events"
  plot945.ylog            = "no"
  plot945.rebin           = 5
  #plot945.xmin            = 0
  #plot945.xmax            = 10
  plot945.ymin            = 0
  plot945.ymax            = 20
  #plot945.lpos = "bottom-center"
  plot945.name            = "MET_over_CaloMET_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot945.addZUncBand     = zUncBand
  plot945.makeRatio       = makeRatio
  plot945.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

  #--- h1_MTenu_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5 ---
  variableName = "h1_MTenu_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"

  plot946 = Plot()
  ## inputs for stacked histograms
  plot946.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot946.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot946.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot946.keys            = keys
  plot946.xtit            = "MTenu - highMeMinusj, #Delta#phi(MET,1^{st} jet)>2.5"
  plot946.ytit            = "Number of events"
  plot946.ylog            = "no"
  plot946.rebin           = "var"
  #plot946.xmin            = 0
  #plot946.xmax            = 1000
  plot946.ymin            = 0
  plot946.ymax            = 20
  #plot946.lpos = "bottom-center"
  plot946.name            = "MTenu_presel_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5"
  plot946.addZUncBand     = zUncBand
  plot946.makeRatio       = makeRatio
  plot946.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot946.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)



  # mDeltaPhiMET1stJet<2.5

  #--- h1_Pt1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Pt1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1000 = Plot()
  ## inputs for stacked histograms
  plot1000.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1000.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1000.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1000.keys            = keys
  plot1000.xtit            = "p_{T} 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1000.ytit            = "Number of events"
  plot1000.ylog            = "no"
  plot1000.rebin           = "var"
  #plot1000.xmin            = 0
  #plot1000.xmax            = 1000
  plot1000.ymin            = 0
  plot1000.ymax            = 15
  #plot1000.lpos = "bottom-center"
  plot1000.name            = "Pt1stJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1000.addZUncBand     = zUncBand
  plot1000.makeRatio       = makeRatio
  plot1000.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot1000.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Eta1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1001 = Plot()
  ## inputs for stacked histograms
  plot1001.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1001.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1001.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1001.keys            = keys
  plot1001.xtit            = "#eta 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1001.ytit            = "Number of events"
  plot1001.ylog            = "no"
  plot1001.rebin           = 5
  #plot1001.xmin            = -5
  #plot1001.xmax            = 5
  plot1001.ymin            = 0
  plot1001.ymax            = 15
  #plot1001.lpos = "bottom-center"
  plot1001.name            = "Eta1stJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1001.addZUncBand     = zUncBand
  plot1001.makeRatio       = makeRatio
  plot1001.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Phi1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1002 = Plot()
  ## inputs for stacked histograms
  plot1002.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1002.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1002.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1002.keys            = keys
  plot1002.xtit            = "#phi 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1002.ytit            = "Number of events"
  plot1002.ylog            = "no"
  plot1002.rebin           = 10
  #plot1002.xmin            = -3.1416
  #plot1002.xmax            = 3.1416
  plot1002.ymin            = 0
  plot1002.ymax            = 20
  #plot1002.lpos = "bottom-center"
  plot1002.name            = "Phi1stJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1002.addZUncBand     = zUncBand
  plot1002.makeRatio       = makeRatio
  plot1002.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CHF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_CHF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1003 = Plot()
  ## inputs for stacked histograms
  plot1003.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1003.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1003.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1003.keys            = keys
  plot1003.xtit            = "Charged hadron energy fraction 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1003.ytit            = "Number of events"
  plot1003.ylog            = "no"
  plot1003.rebin           = 10
  #plot1003.xmin            = 0
  #plot1003.xmax            = 1
  plot1003.ymin            = 0
  plot1003.ymax            = 20
  #plot1003.lpos = "bottom-center"
  plot1003.name            = "CHF1stJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1003.addZUncBand     = zUncBand
  plot1003.makeRatio       = makeRatio
  plot1003.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NHF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_NHF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1004 = Plot()
  ## inputs for stacked histograms
  plot1004.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1004.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1004.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1004.keys            = keys
  plot1004.xtit            = "Neutral hadron energy fraction 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1004.ytit            = "Number of events"
  plot1004.ylog            = "no"
  plot1004.rebin           = 10
  #plot1004.xmin            = 0
  #plot1004.xmax            = 1
  plot1004.ymin            = 0
  plot1004.ymax            = 40
  #plot1004.lpos = "bottom-center"
  plot1004.name            = "NHF1stJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1004.addZUncBand     = zUncBand
  plot1004.makeRatio       = makeRatio
  plot1004.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CEF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_CEF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1005 = Plot()
  ## inputs for stacked histograms
  plot1005.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1005.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1005.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1005.keys            = keys
  plot1005.xtit            = "Charged EM energy fraction 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1005.ytit            = "Number of events"
  plot1005.ylog            = "no"
  plot1005.rebin           = 10
  #plot1005.xmin            = 0
  #plot1005.xmax            = 1
  plot1005.ymin            = 0
  plot1005.ymax            = 60
  #plot1005.lpos = "bottom-center"
  plot1005.name            = "CEF1stJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1005.addZUncBand     = zUncBand
  plot1005.makeRatio       = makeRatio
  plot1005.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NEF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_NEF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1006 = Plot()
  ## inputs for stacked histograms
  plot1006.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1006.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1006.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1006.keys            = keys
  plot1006.xtit            = "Neutral EM energy fraction 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1006.ytit            = "Number of events"
  plot1006.ylog            = "no"
  plot1006.rebin           = 10
  #plot1006.xmin            = 0
  #plot1006.xmax            = 1
  plot1006.ymin            = 0
  plot1006.ymax            = 20
  #plot1006.lpos = "bottom-center"
  plot1006.name            = "NEF1stJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1006.addZUncBand     = zUncBand
  plot1006.makeRatio       = makeRatio
  plot1006.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NCH1stJet_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_NCH1stJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1007 = Plot()
  ## inputs for stacked histograms
  plot1007.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1007.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1007.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1007.keys            = keys
  plot1007.xtit            = "Charged multiplicity 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1007.ytit            = "Number of events"
  plot1007.ylog            = "no"
  plot1007.rebin           = 10
  #plot1007.xmin            = -0.5
  #plot1007.xmax            = 99.5
  plot1007.ymin            = 0
  plot1007.ymax            = 40
  #plot1007.lpos = "bottom-center"
  plot1007.name            = "NCH1stJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1007.addZUncBand     = zUncBand
  plot1007.makeRatio       = makeRatio
  plot1007.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NN1stJet_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_NN1stJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1008 = Plot()
  ## inputs for stacked histograms
  plot1008.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1008.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1008.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1008.keys            = keys
  plot1008.xtit            = "Neutral multiplicity 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1008.ytit            = "Number of events"
  plot1008.ylog            = "no"
  plot1008.rebin           = 10
  #plot1008.xmin            = -0.5
  #plot1008.xmax            = 99.5
  plot1008.ymin            = 0
  plot1008.ymax            = 30
  #plot1008.lpos = "bottom-center"
  plot1008.name            = "NN1stJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1008.addZUncBand     = zUncBand
  plot1008.makeRatio       = makeRatio
  plot1008.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NC1stJet_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_NC1stJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1009 = Plot()
  ## inputs for stacked histograms
  plot1009.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1009.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1009.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1009.keys            = keys
  plot1009.xtit            = "Number of constituents 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1009.ytit            = "Number of events"
  plot1009.ylog            = "no"
  plot1009.rebin           = 10
  #plot1009.xmin            = -0.5
  #plot1009.xmax            = 99.5
  plot1009.ymin            = 0
  plot1009.ymax            = 40
  #plot1009.lpos = "bottom-center"
  plot1009.name            = "NC1stJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1009.addZUncBand     = zUncBand
  plot1009.makeRatio       = makeRatio
  plot1009.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Pt2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1010 = Plot()
  ## inputs for stacked histograms
  plot1010.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1010.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1010.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1010.keys            = keys
  plot1010.xtit            = "p_{T} 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1010.ytit            = "Number of events"
  plot1010.ylog            = "no"
  plot1010.rebin           = "var"
  #plot1010.xmin            = 0
  #plot1010.xmax            = 1000
  plot1010.ymin            = 0
  plot1010.ymax            = 20
  #plot1010.lpos = "bottom-center"
  plot1010.name            = "Pt2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1010.addZUncBand     = zUncBand
  plot1010.makeRatio       = makeRatio
  plot1010.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot1010.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Eta2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1011 = Plot()
  ## inputs for stacked histograms
  plot1011.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1011.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1011.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1011.keys            = keys
  plot1011.xtit            = "#eta 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1011.ytit            = "Number of events"
  plot1011.ylog            = "no"
  plot1011.rebin           = 5
  #plot1011.xmin            = -5
  #plot1011.xmax            = 5
  plot1011.ymin            = 0
  plot1011.ymax            = 15
  #plot1011.lpos = "bottom-center"
  plot1011.name            = "Eta2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1011.addZUncBand     = zUncBand
  plot1011.makeRatio       = makeRatio
  plot1011.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Phi2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1012 = Plot()
  ## inputs for stacked histograms
  plot1012.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1012.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1012.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1012.keys            = keys
  plot1012.xtit            = "#phi 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1012.ytit            = "Number of events"
  plot1012.ylog            = "no"
  plot1012.rebin           = 10
  #plot1012.xmin            = -3.1416
  #plot1012.xmax            = 3.1416
  plot1012.ymin            = 0
  plot1012.ymax            = 20
  #plot1012.lpos = "bottom-center"
  plot1012.name            = "Phi2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1012.addZUncBand     = zUncBand
  plot1012.makeRatio       = makeRatio
  plot1012.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CHF2ndJet_highMej_mDeltaPhiMET2ndJet_le_2.5 ---
  variableName = "h1_CHF2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1013 = Plot()
  ## inputs for stacked histograms
  plot1013.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1013.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1013.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1013.keys            = keys
  plot1013.xtit            = "Charged hadron energy fraction 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1013.ytit            = "Number of events"
  plot1013.ylog            = "no"
  plot1013.rebin           = 10
  #plot1013.xmin            = 0
  #plot1013.xmax            = 1
  plot1013.ymin            = 0
  plot1013.ymax            = 20
  #plot1013.lpos = "bottom-center"
  plot1013.name            = "CHF2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1013.addZUncBand     = zUncBand
  plot1013.makeRatio       = makeRatio
  plot1013.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NHF2ndJet_highMej_mDeltaPhiMET2ndJet_le_2.5 ---
  variableName = "h1_NHF2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1014 = Plot()
  ## inputs for stacked histograms
  plot1014.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1014.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1014.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1014.keys            = keys
  plot1014.xtit            = "Neutral hadron energy fraction 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1014.ytit            = "Number of events"
  plot1014.ylog            = "no"
  plot1014.rebin           = 10
  #plot1014.xmin            = 0
  #plot1014.xmax            = 1
  plot1014.ymin            = 0
  plot1014.ymax            = 40
  #plot1014.lpos = "bottom-center"
  plot1014.name            = "NHF2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1014.addZUncBand     = zUncBand
  plot1014.makeRatio       = makeRatio
  plot1014.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CEF2ndJet_highMej_mDeltaPhiMET2ndJet_le_2.5 ---
  variableName = "h1_CEF2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1015 = Plot()
  ## inputs for stacked histograms
  plot1015.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1015.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1015.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1015.keys            = keys
  plot1015.xtit            = "Charged EM energy fraction 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1015.ytit            = "Number of events"
  plot1015.ylog            = "no"
  plot1015.rebin           = 10
  #plot1015.xmin            = 0
  #plot1015.xmax            = 1
  plot1015.ymin            = 0
  plot1015.ymax            = 60
  #plot1015.lpos = "bottom-center"
  plot1015.name            = "CEF2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1015.addZUncBand     = zUncBand
  plot1015.makeRatio       = makeRatio
  plot1015.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NEF2ndJet_highMej_mDeltaPhiMET2ndJet_le_2.5 ---
  variableName = "h1_NEF2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1016 = Plot()
  ## inputs for stacked histograms
  plot1016.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1016.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1016.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1016.keys            = keys
  plot1016.xtit            = "Neutral EM energy fraction 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1016.ytit            = "Number of events"
  plot1016.ylog            = "no"
  plot1016.rebin           = 10
  #plot1016.xmin            = 0
  #plot1016.xmax            = 1
  plot1016.ymin            = 0
  plot1016.ymax            = 20
  #plot1016.lpos = "bottom-center"
  plot1016.name            = "NEF2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1016.addZUncBand     = zUncBand
  plot1016.makeRatio       = makeRatio
  plot1016.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NCH2ndJet_highMej_mDeltaPhiMET2ndJet_le_2.5 ---
  variableName = "h1_NCH2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1017 = Plot()
  ## inputs for stacked histograms
  plot1017.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1017.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1017.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1017.keys            = keys
  plot1017.xtit            = "Charged multiplicity 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1017.ytit            = "Number of events"
  plot1017.ylog            = "no"
  plot1017.rebin           = 10
  #plot1017.xmin            = -0.5
  #plot1017.xmax            = 99.5
  plot1017.ymin            = 0
  plot1017.ymax            = 40
  #plot1017.lpos = "bottom-center"
  plot1017.name            = "NCH2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1017.addZUncBand     = zUncBand
  plot1017.makeRatio       = makeRatio
  plot1017.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NN2ndJet_highMej_mDeltaPhiMET2ndJet_le_2.5 ---
  variableName = "h1_NN2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1018 = Plot()
  ## inputs for stacked histograms
  plot1018.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1018.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1018.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1018.keys            = keys
  plot1018.xtit            = "Neutral multiplicity 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1018.ytit            = "Number of events"
  plot1018.ylog            = "no"
  plot1018.rebin           = 10
  #plot1018.xmin            = -0.5
  #plot1018.xmax            = 99.5
  plot1018.ymin            = 0
  plot1018.ymax            = 30
  #plot1018.lpos = "bottom-center"
  plot1018.name            = "NN2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1018.addZUncBand     = zUncBand
  plot1018.makeRatio       = makeRatio
  plot1018.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NC2ndJet_highMej_mDeltaPhiMET2ndJet_le_2.5 ---
  variableName = "h1_NC2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1019 = Plot()
  ## inputs for stacked histograms
  plot1019.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1019.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1019.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1019.keys            = keys
  plot1019.xtit            = "Number of constituents 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1019.ytit            = "Number of events"
  plot1019.ylog            = "no"
  plot1019.rebin           = 10
  #plot1019.xmin            = -0.5
  #plot1019.xmax            = 99.5
  plot1019.ymin            = 0
  plot1019.ymax            = 30
  #plot1019.lpos = "bottom-center"
  plot1019.name            = "NC2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1019.addZUncBand     = zUncBand
  plot1019.makeRatio       = makeRatio
  plot1019.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_E1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_E1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1020 = Plot()
  ## inputs for stacked histograms
  plot1020.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1020.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1020.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1020.keys            = keys
  plot1020.xtit            = "Energy 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1020.ytit            = "Number of events"
  plot1020.ylog            = "no"
  plot1020.rebin           = "var"
  #plot1020.xmin            = 0
  #plot1020.xmax            = 1000
  plot1020.ymin            = 0
  plot1020.ymax            = 20
  #plot1020.lpos = "bottom-center"
  plot1020.name            = "E1stEle_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1020.addZUncBand     = zUncBand
  plot1020.makeRatio       = makeRatio
  plot1020.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot1020.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Pt1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1021 = Plot()
  ## inputs for stacked histograms
  plot1021.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1021.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1021.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1021.keys            = keys
  plot1021.xtit            = "p_{T} 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1021.ytit            = "Number of events"
  plot1021.ylog            = "no"
  plot1021.rebin           = "var"
  #plot1021.xmin            = 0
  #plot1021.xmax            = 1000
  plot1021.ymin            = 0
  plot1021.ymax            = 20
  #plot1021.lpos = "bottom-center"
  plot1021.name            = "Pt1stEle_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1021.addZUncBand     = zUncBand
  plot1021.makeRatio       = makeRatio
  plot1021.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot1021.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Eta1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1022 = Plot()
  ## inputs for stacked histograms
  plot1022.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1022.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1022.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1022.keys            = keys
  plot1022.xtit            = "#eta 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1022.ytit            = "Number of events"
  plot1022.ylog            = "no"
  plot1022.rebin           = 5
  #plot1022.xmin            = -5
  #plot1022.xmax            = 5
  plot1022.ymin            = 0
  plot1022.ymax            = 15
  #plot1022.lpos = "bottom-center"
  plot1022.name            = "Eta1stEle_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1022.addZUncBand     = zUncBand
  plot1022.makeRatio       = makeRatio
  plot1022.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Phi1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1023 = Plot()
  ## inputs for stacked histograms
  plot1023.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1023.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1023.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1023.keys            = keys
  plot1023.xtit            = "#phi 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1023.ytit            = "Number of events"
  plot1023.ylog            = "no"
  plot1023.rebin           = 10
  #plot1023.xmin            = -3.1416
  #plot1023.xmax            = 3.1416
  plot1023.ymin            = 0
  plot1023.ymax            = 20
  #plot1023.lpos = "bottom-center"
  plot1023.name            = "Phi1stEle_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1023.addZUncBand     = zUncBand
  plot1023.makeRatio       = makeRatio
  plot1023.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Charge1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Charge1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1024 = Plot()
  ## inputs for stacked histograms
  plot1024.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1024.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1024.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1024.keys            = keys
  plot1024.xtit            = "Charge 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1024.ytit            = "Number of events"
  plot1024.ylog            = "no"
  #plot1024.rebin           = 1
  #plot1024.xmin            = -1.001
  #plot1024.xmax            = 1.001
  plot1024.ymin            = 0
  plot1024.ymax            = 50
  #plot1024.lpos = "bottom-center"
  plot1024.name            = "Charge1stEle_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1024.addZUncBand     = zUncBand
  plot1024.makeRatio       = makeRatio
  plot1024.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Conversion1stEle_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Conversion1stEle_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1025 = Plot()
  ## inputs for stacked histograms
  plot1025.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot1025.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1025.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1025.keys            = keys
  plot1025.xtit            = "Conversion flag 1^{st} electron - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1025.ytit            = "Number of events"
  plot1025.ylog            = "no"
  #plot1025.rebin           = 1
  #plot1025.xmin            = -0.5
  #plot1025.xmax            = 1.5
  plot1025.ymin            = 0
  plot1025.ymax            = 60
  #plot1025.lpos = "bottom-center"
  plot1025.name            = "Conversion1stEle_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1025.addZUncBand     = zUncBand
  plot1025.makeRatio       = makeRatio
  plot1025.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMETEle_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_mDeltaPhiMETEle_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1026 = Plot()
  ## inputs for stacked histograms
  plot1026.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1026.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1026.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1026.keys            = keys
  plot1026.xtit            = "#Delta#phi(MET,1^{st} ele) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1026.ytit            = "Number of events"
  plot1026.ylog            = "no"
  plot1026.rebin           = 10
  #plot1026.xmin            = 0
  #plot1026.xmax            = 3.1416
  plot1026.ymin            = 0
  plot1026.ymax            = 20
  #plot1026.lpos = "bottom-center"
  plot1026.name            = "mDeltaPhiMETEle_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1026.addZUncBand     = zUncBand
  plot1026.makeRatio       = makeRatio
  plot1026.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_mDeltaPhiEle1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1027 = Plot()
  ## inputs for stacked histograms
  plot1027.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1027.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1027.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1027.keys            = keys
  plot1027.xtit            = "#Delta#phi(1^{st} ele,1^{st} jet) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1027.ytit            = "Number of events"
  plot1027.ylog            = "no"
  plot1027.rebin           = 10
  #plot1027.xmin            = 0
  #plot1027.xmax            = 3.1416
  plot1027.ymin            = 0
  plot1027.ymax            = 20
  #plot1027.lpos = "bottom-center"
  plot1027.name            = "mDeltaPhiEle1stJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1027.addZUncBand     = zUncBand
  plot1027.makeRatio       = makeRatio
  plot1027.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_mDeltaPhiEle2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1028 = Plot()
  ## inputs for stacked histograms
  plot1028.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1028.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1028.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1028.keys            = keys
  plot1028.xtit            = "#Delta#phi(1^{st} ele,2^{nd} jet) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1028.ytit            = "Number of events"
  plot1028.ylog            = "no"
  plot1028.rebin           = 10
  #plot1028.xmin            = 0
  #plot1028.xmax            = 3.1416
  plot1028.ymin            = 0
  plot1028.ymax            = 20
  #plot1028.lpos = "bottom-center"
  plot1028.name            = "mDeltaPhiEle2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1028.addZUncBand     = zUncBand
  plot1028.makeRatio       = makeRatio
  plot1028.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_mDeltaPhiMET2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1029 = Plot()
  ## inputs for stacked histograms
  plot1029.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1029.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1029.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1029.keys            = keys
  plot1029.xtit            = "#Delta#phi(MET,2^{nd} jet) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1029.ytit            = "Number of events"
  plot1029.ylog            = "no"
  plot1029.rebin           = 10
  #plot1029.xmin            = 0
  #plot1029.xmax            = 3.1416
  plot1029.ymin            = 0
  plot1029.ymax            = 20
  #plot1029.lpos = "bottom-center"
  plot1029.name            = "mDeltaPhiMET2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1029.addZUncBand     = zUncBand
  plot1029.makeRatio       = makeRatio
  plot1029.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Ptenu_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Ptenu_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1030 = Plot()
  ## inputs for stacked histograms
  plot1030.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1030.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1030.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1030.keys            = keys
  plot1030.xtit            = "p_{T}(e#nu) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1030.ytit            = "Number of events"
  plot1030.ylog            = "no"
  plot1030.rebin           = "var"
  #plot1030.xmin            = 0
  #plot1030.xmax            = 2000
  plot1030.ymin            = 0
  plot1030.ymax            = 20
  #plot1030.lpos = "bottom-center"
  plot1030.name            = "Ptenu_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1030.addZUncBand     = zUncBand
  plot1030.makeRatio       = makeRatio
  plot1030.xbins           = [0,40,80,120,160,200,300,600]
  plot1030.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_1stJet_PTOverPTPlusMET_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_1stJet_PTOverPTPlusMET_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1031 = Plot()
  ## inputs for stacked histograms
  plot1031.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1031.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1031.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1031.keys            = keys
  plot1031.xtit            = "Pt1stJet/(Pt1stJet+MET) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1031.ytit            = "Number of events"
  plot1031.ylog            = "no"
  plot1031.rebin           = 5
  #plot1031.xmin            = 0
  #plot1031.xmax            = 1
  plot1031.ymin            = 0
  plot1031.ymax            = 20
  #plot1031.lpos = "bottom-center"
  plot1031.name            = "1stJet_PTOverPTPlusMET_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1031.addZUncBand     = zUncBand
  plot1031.makeRatio       = makeRatio
  plot1031.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MET_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_MET_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1032 = Plot()
  ## inputs for stacked histograms
  plot1032.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1032.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1032.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1032.keys            = keys
  plot1032.xtit            = "pfMET - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1032.ytit            = "Number of events"
  plot1032.ylog            = "no"
  plot1032.rebin           = "var"
  #plot1032.xmin            = 0
  #plot1032.xmax            = 1000
  plot1032.ymin            = 0
  plot1032.ymax            = 20
  #plot1032.lpos = "bottom-center"
  plot1032.name            = "MET_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1032.addZUncBand     = zUncBand
  plot1032.makeRatio       = makeRatio
  plot1032.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot1032.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Njet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1033 = Plot()
  ## inputs for stacked histograms
  plot1033.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1033.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1033.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1033.keys            = keys
  plot1033.xtit            = "Number of jets - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1033.ytit            = "Number of events"
  plot1033.ylog            = "no"
  plot1033.rebin           = 1
  plot1033.xmin            = -0.5
  plot1033.xmax            = 11.5
  plot1033.ymin            = 0
  plot1033.ymax            = 25
  #plot1033.lpos = "bottom-center"
  plot1033.name            = "Njet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1033.addZUncBand     = zUncBand
  plot1033.makeRatio       = makeRatio
  plot1033.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_NjetTCHELBTag_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1034 = Plot()
  ## inputs for stacked histograms
  plot1034.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1034.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1034.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1034.keys            = keys
  plot1034.xtit            = "Number of b-tagged jets (TCHEL) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1034.ytit            = "Number of events"
  plot1034.ylog            = "no"
  plot1034.rebin           = 1
  plot1034.xmin            = -0.5
  plot1034.xmax            = 11.5
  plot1034.ymin            = 0
  plot1034.ymax            = 25
  #plot1034.lpos = "bottom-center"
  plot1034.name            = "NjetTCHELBTag_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1034.addZUncBand     = zUncBand
  plot1034.makeRatio       = makeRatio
  plot1034.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_minDRej_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_minDRej_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1035 = Plot()
  ## inputs for stacked histograms
  plot1035.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1035.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1035.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1035.keys            = keys
  plot1035.xtit            = "min#DeltaR(e,j) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1035.ytit            = "Number of events"
  plot1035.ylog            = "no"
  plot1035.rebin           = 10
  #plot1035.xmin            = 0
  #plot1035.xmax            = 10
  plot1035.ymin            = 0
  plot1035.ymax            = 25
  #plot1035.lpos = "bottom-center"
  plot1035.name            = "minDRej_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1035.addZUncBand     = zUncBand
  plot1035.makeRatio       = makeRatio
  plot1035.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_maxDRej_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_maxDRej_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1036 = Plot()
  ## inputs for stacked histograms
  plot1036.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1036.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1036.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1036.keys            = keys
  plot1036.xtit            = "max#DeltaR(e,j) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1036.ytit            = "Number of events"
  plot1036.ylog            = "no"
  plot1036.rebin           = 10
  #plot1036.xmin            = 0
  #plot1036.xmax            = 10
  plot1036.ymin            = 0
  plot1036.ymax            = 30
  #plot1036.lpos = "bottom-center"
  plot1036.name            = "maxDRej_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1036.addZUncBand     = zUncBand
  plot1036.makeRatio       = makeRatio
  plot1036.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DRjets_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_DRjets_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1037 = Plot()
  ## inputs for stacked histograms
  plot1037.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1037.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1037.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1037.keys            = keys
  plot1037.xtit            = "#DeltaR(j1,j2) - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1037.ytit            = "Number of events"
  plot1037.ylog            = "no"
  plot1037.rebin           = 10
  #plot1037.xmin            = 0
  #plot1037.xmax            = 10
  plot1037.ymin            = 0
  plot1037.ymax            = 30
  #plot1037.lpos = "bottom-center"
  plot1037.name            = "DRjets_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1037.addZUncBand     = zUncBand
  plot1037.makeRatio       = makeRatio
  plot1037.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Mass1stJet_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Mass1stJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1038 = Plot()
  ## inputs for stacked histograms
  plot1038.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1038.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1038.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1038.keys            = keys
  plot1038.xtit            = "Mass 1^{st} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1038.ytit            = "Number of events"
  plot1038.ylog            = "no"
  plot1038.rebin           = 5
  #plot1038.xmin            = -0.5
  #plot1038.xmax            = 99.5
  plot1038.ymin            = 0
  plot1038.ymax            = 20
  #plot1038.lpos = "bottom-center"
  plot1038.name            = "Mass1stJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1038.addZUncBand     = zUncBand
  plot1038.makeRatio       = makeRatio
  plot1038.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Mass2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_Mass2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1042 = Plot()
  ## inputs for stacked histograms
  plot1042.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1042.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1042.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1042.keys            = keys
  plot1042.xtit            = "Mass 2^{nd} jet - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1042.ytit            = "Number of events"
  plot1042.ylog            = "no"
  plot1042.rebin           = 5
  #plot1042.xmin            = -0.5
  #plot1042.xmax            = 99.5
  plot1042.ymin            = 0
  plot1042.ymax            = 20
  #plot1042.lpos = "bottom-center"
  plot1042.name            = "Mass2ndJet_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1042.addZUncBand     = zUncBand
  plot1042.makeRatio       = makeRatio
  plot1042.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CaloMET_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_CaloMET_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1043 = Plot()
  ## inputs for stacked histograms
  plot1043.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1043.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1043.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1043.keys            = keys
  plot1043.xtit            = "caloMET - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1043.ytit            = "Number of events"
  plot1043.ylog            = "no"
  plot1043.rebin           = "var"
  #plot1043.xmin            = 0
  #plot1043.xmax            = 1000
  plot1043.ymin            = 0
  plot1043.ymax            = 20
  #plot1043.lpos = "bottom-center"
  plot1043.name            = "CaloMET_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1043.addZUncBand     = zUncBand
  plot1043.makeRatio       = makeRatio
  plot1043.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot1043.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_TCMET_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_TCMET_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1044 = Plot()
  ## inputs for stacked histograms
  plot1044.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1044.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1044.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1044.keys            = keys
  plot1044.xtit            = "tcMET - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1044.ytit            = "Number of events"
  plot1044.ylog            = "no"
  plot1044.rebin           = "var"
  #plot1044.xmin            = 0
  #plot1044.xmax            = 1000
  plot1044.ymin            = 0
  plot1044.ymax            = 20
  #plot1044.lpos = "bottom-center"
  plot1044.name            = "TCMET_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1044.addZUncBand     = zUncBand
  plot1044.makeRatio       = makeRatio
  plot1044.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot1044.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

  #--- h1_MET_over_CaloMET_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_MET_over_CaloMET_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1045 = Plot()
  ## inputs for stacked histograms
  plot1045.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1045.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1045.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1045.keys            = keys
  plot1045.xtit            = "pfMET/CaloMET - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1045.ytit            = "Number of events"
  plot1045.ylog            = "no"
  plot1045.rebin           = 5
  #plot1045.xmin            = 0
  #plot1045.xmax            = 10
  plot1045.ymin            = 0
  plot1045.ymax            = 20
  #plot1045.lpos = "bottom-center"
  plot1045.name            = "MET_over_CaloMET_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1045.addZUncBand     = zUncBand
  plot1045.makeRatio       = makeRatio
  plot1045.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MTenu_PAS_highMej_mDeltaPhiMET1stJet_le_2.5 ---
  variableName = "h1_MTenu_PAS_highMej_mDeltaPhiMET1stJet_le_2.5"

  plot1046 = Plot()
  ## inputs for stacked histograms
  plot1046.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1046.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1046.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1046.keys            = keys
  plot1046.xtit            = "MTenu - highMej, #Delta#phi(MET,1^{st} jet)<2.5"
  plot1046.ytit            = "Number of events"
  plot1046.ylog            = "no"
  plot1046.rebin           = "var"
  #plot1046.xmin            = 0
  #plot1046.xmax            = 1000
  plot1046.ymin            = 0
  plot1046.ymax            = 20
  #plot1046.lpos = "bottom-center"
  plot1046.name            = "MTenu_presel_highMej_mDeltaPhiMET1stJet_le_2.5"
  plot1046.addZUncBand     = zUncBand
  plot1046.makeRatio       = makeRatio
  plot1046.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot1046.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  # ## Additional plots to investigate bump in Eta1stJet distribution

  #--- h1_Pt1stJet_PAS_Eta1stJetBump ---
  variableName = "h1_Pt1stJet_PAS_Eta1stJetBump"

  plot1100 = Plot()
  ## inputs for stacked histograms
  plot1100.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1100.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1100.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1100.keys            = keys
  plot1100.xtit            = "p_{T} 1^{st} jet - Eta1stJetBump"
  plot1100.ytit            = "Number of events"
  plot1100.ylog            = "no"
  plot1100.rebin           = "var"
  #plot1100.xmin            = 0
  #plot1100.xmax            = 1000
  plot1100.ymin            = 0
  plot1100.ymax            = 20
  #plot1100.lpos = "bottom-center"
  plot1100.name            = "Pt1stJet_presel_Eta1stJetBump"
  plot1100.addZUncBand     = zUncBand
  plot1100.makeRatio       = makeRatio
  plot1100.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot1100.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  ##--- h1_Eta1stJet_PAS_Eta1stJetBump ---
  #variableName = "h1_Eta1stJet_PAS_Eta1stJetBump"

  #plot1101 = Plot()
  ### inputs for stacked histograms
  #plot1101.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  #plot1101.keysStack       = keysStack
  ### this is the list of histograms that should be simply overlaid on top of the stacked histogram
  #plot1101.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  #plot1101.keys            = keys
  #plot1101.xtit            = "#eta 1^{st} jet - Eta1stJetBump"
  #plot1101.ytit            = "Number of events"
  #plot1101.ylog            = "no"
  #plot1101.rebin           = 5
  ##plot1101.xmin            = -5
  ##plot1101.xmax            = 5
  #plot1101.ymin            = 0
  #plot1101.ymax            = 15
  ##plot1101.lpos = "bottom-center"
  #plot1101.name            = "Eta1stJet_presel_Eta1stJetBump"
  #plot1101.addZUncBand     = zUncBand
  #plot1101.makeRatio       = makeRatio
  #plot1101.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stJet_PAS_Eta1stJetBump ---
  variableName = "h1_Phi1stJet_PAS_Eta1stJetBump"

  plot1102 = Plot()
  ## inputs for stacked histograms
  plot1102.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1102.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1102.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1102.keys            = keys
  plot1102.xtit            = "#phi 1^{st} jet - Eta1stJetBump"
  plot1102.ytit            = "Number of events"
  plot1102.ylog            = "no"
  plot1102.rebin           = 10
  #plot1102.xmin            = -3.1416
  #plot1102.xmax            = 3.1416
  plot1102.ymin            = 0
  plot1102.ymax            = 20
  #plot1102.lpos = "bottom-center"
  plot1102.name            = "Phi1stJet_presel_Eta1stJetBump"
  plot1102.addZUncBand     = zUncBand
  plot1102.makeRatio       = makeRatio
  plot1102.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CHF1stJet_Eta1stJetBump ---
  variableName = "h1_CHF1stJet_Eta1stJetBump"

  plot1103 = Plot()
  ## inputs for stacked histograms
  plot1103.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1103.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1103.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1103.keys            = keys
  plot1103.xtit            = "Charged hadron energy fraction 1^{st} jet - Eta1stJetBump"
  plot1103.ytit            = "Number of events"
  plot1103.ylog            = "no"
  plot1103.rebin           = 10
  #plot1103.xmin            = 0
  #plot1103.xmax            = 1
  plot1103.ymin            = 0
  plot1103.ymax            = 20
  #plot1103.lpos = "bottom-center"
  plot1103.name            = "CHF1stJet_presel_Eta1stJetBump"
  plot1103.addZUncBand     = zUncBand
  plot1103.makeRatio       = makeRatio
  plot1103.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NHF1stJet_Eta1stJetBump ---
  variableName = "h1_NHF1stJet_Eta1stJetBump"

  plot1104 = Plot()
  ## inputs for stacked histograms
  plot1104.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1104.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1104.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1104.keys            = keys
  plot1104.xtit            = "Neutral hadron energy fraction 1^{st} jet - Eta1stJetBump"
  plot1104.ytit            = "Number of events"
  plot1104.ylog            = "no"
  plot1104.rebin           = 10
  #plot1104.xmin            = 0
  #plot1104.xmax            = 1
  plot1104.ymin            = 0
  plot1104.ymax            = 40
  #plot1104.lpos = "bottom-center"
  plot1104.name            = "NHF1stJet_presel_Eta1stJetBump"
  plot1104.addZUncBand     = zUncBand
  plot1104.makeRatio       = makeRatio
  plot1104.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CEF1stJet_Eta1stJetBump ---
  variableName = "h1_CEF1stJet_Eta1stJetBump"

  plot1105 = Plot()
  ## inputs for stacked histograms
  plot1105.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1105.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1105.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1105.keys            = keys
  plot1105.xtit            = "Charged EM energy fraction 1^{st} jet - Eta1stJetBump"
  plot1105.ytit            = "Number of events"
  plot1105.ylog            = "no"
  plot1105.rebin           = 10
  #plot1105.xmin            = 0
  #plot1105.xmax            = 1
  plot1105.ymin            = 0
  plot1105.ymax            = 60
  #plot1105.lpos = "bottom-center"
  plot1105.name            = "CEF1stJet_presel_Eta1stJetBump"
  plot1105.addZUncBand     = zUncBand
  plot1105.makeRatio       = makeRatio
  plot1105.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NEF1stJet_Eta1stJetBump ---
  variableName = "h1_NEF1stJet_Eta1stJetBump"

  plot1106 = Plot()
  ## inputs for stacked histograms
  plot1106.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1106.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1106.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1106.keys            = keys
  plot1106.xtit            = "Neutral EM energy fraction 1^{st} jet - Eta1stJetBump"
  plot1106.ytit            = "Number of events"
  plot1106.ylog            = "no"
  plot1106.rebin           = 10
  #plot1106.xmin            = 0
  #plot1106.xmax            = 1
  plot1106.ymin            = 0
  plot1106.ymax            = 20
  #plot1106.lpos = "bottom-center"
  plot1106.name            = "NEF1stJet_presel_Eta1stJetBump"
  plot1106.addZUncBand     = zUncBand
  plot1106.makeRatio       = makeRatio
  plot1106.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NCH1stJet_Eta1stJetBump ---
  variableName = "h1_NCH1stJet_Eta1stJetBump"

  plot1107 = Plot()
  ## inputs for stacked histograms
  plot1107.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1107.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1107.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1107.keys            = keys
  plot1107.xtit            = "Charged multiplicity 1^{st} jet - Eta1stJetBump"
  plot1107.ytit            = "Number of events"
  plot1107.ylog            = "no"
  plot1107.rebin           = 10
  #plot1107.xmin            = -0.5
  #plot1107.xmax            = 99.5
  plot1107.ymin            = 0
  plot1107.ymax            = 40
  #plot1107.lpos = "bottom-center"
  plot1107.name            = "NCH1stJet_presel_Eta1stJetBump"
  plot1107.addZUncBand     = zUncBand
  plot1107.makeRatio       = makeRatio
  plot1107.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NN1stJet_Eta1stJetBump ---
  variableName = "h1_NN1stJet_Eta1stJetBump"

  plot1108 = Plot()
  ## inputs for stacked histograms
  plot1108.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1108.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1108.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1108.keys            = keys
  plot1108.xtit            = "Neutral multiplicity 1^{st} jet - Eta1stJetBump"
  plot1108.ytit            = "Number of events"
  plot1108.ylog            = "no"
  plot1108.rebin           = 10
  #plot1108.xmin            = -0.5
  #plot1108.xmax            = 99.5
  plot1108.ymin            = 0
  plot1108.ymax            = 40
  #plot1108.lpos = "bottom-center"
  plot1108.name            = "NN1stJet_presel_Eta1stJetBump"
  plot1108.addZUncBand     = zUncBand
  plot1108.makeRatio       = makeRatio
  plot1108.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NC1stJet_Eta1stJetBump ---
  variableName = "h1_NC1stJet_Eta1stJetBump"

  plot1109 = Plot()
  ## inputs for stacked histograms
  plot1109.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1109.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1109.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1109.keys            = keys
  plot1109.xtit            = "Number of constituents 1^{st} jet - Eta1stJetBump"
  plot1109.ytit            = "Number of events"
  plot1109.ylog            = "no"
  plot1109.rebin           = 10
  #plot1109.xmin            = -0.5
  #plot1109.xmax            = 99.5
  plot1109.ymin            = 0
  plot1109.ymax            = 40
  #plot1109.lpos = "bottom-center"
  plot1109.name            = "NC1stJet_presel_Eta1stJetBump"
  plot1109.addZUncBand     = zUncBand
  plot1109.makeRatio       = makeRatio
  plot1109.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt2ndJet_PAS_Eta1stJetBump ---
  variableName = "h1_Pt2ndJet_PAS_Eta1stJetBump"

  plot1110 = Plot()
  ## inputs for stacked histograms
  plot1110.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1110.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1110.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1110.keys            = keys
  plot1110.xtit            = "p_{T} 2^{nd} jet - Eta1stJetBump"
  plot1110.ytit            = "Number of events"
  plot1110.ylog            = "no"
  plot1110.rebin           = "var"
  #plot1110.xmin            = 0
  #plot1110.xmax            = 1000
  plot1110.ymin            = 0
  plot1110.ymax            = 25
  #plot1110.lpos = "bottom-center"
  plot1110.name            = "Pt2ndJet_presel_Eta1stJetBump"
  plot1110.addZUncBand     = zUncBand
  plot1110.makeRatio       = makeRatio
  plot1110.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot1110.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta2ndJet_PAS_Eta1stJetBump ---
  variableName = "h1_Eta2ndJet_PAS_Eta1stJetBump"

  plot1111 = Plot()
  ## inputs for stacked histograms
  plot1111.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1111.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1111.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1111.keys            = keys
  plot1111.xtit            = "#eta 2^{nd} jet - Eta1stJetBump"
  plot1111.ytit            = "Number of events"
  plot1111.ylog            = "no"
  plot1111.rebin           = 5
  #plot1111.xmin            = -5
  #plot1111.xmax            = 5
  plot1111.ymin            = 0
  plot1111.ymax            = 15
  #plot1111.lpos = "bottom-center"
  plot1111.name            = "Eta2ndJet_presel_Eta1stJetBump"
  plot1111.addZUncBand     = zUncBand
  plot1111.makeRatio       = makeRatio
  plot1111.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi2ndJet_PAS_Eta1stJetBump ---
  variableName = "h1_Phi2ndJet_PAS_Eta1stJetBump"

  plot1112 = Plot()
  ## inputs for stacked histograms
  plot1112.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1112.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1112.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1112.keys            = keys
  plot1112.xtit            = "#phi 2^{nd} jet - Eta1stJetBump"
  plot1112.ytit            = "Number of events"
  plot1112.ylog            = "no"
  plot1112.rebin           = 10
  #plot1112.xmin            = -3.1416
  #plot1112.xmax            = 3.1416
  plot1112.ymin            = 0
  plot1112.ymax            = 20
  #plot1112.lpos = "bottom-center"
  plot1112.name            = "Phi2ndJet_presel_Eta1stJetBump"
  plot1112.addZUncBand     = zUncBand
  plot1112.makeRatio       = makeRatio
  plot1112.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CHF2ndJet_highMej_Eta1stJetBump ---
  variableName = "h1_CHF2ndJet_Eta1stJetBump"

  plot1113 = Plot()
  ## inputs for stacked histograms
  plot1113.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1113.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1113.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1113.keys            = keys
  plot1113.xtit            = "Charged hadron energy fraction 2^{nd} jet - Eta1stJetBump"
  plot1113.ytit            = "Number of events"
  plot1113.ylog            = "no"
  plot1113.rebin           = 10
  #plot1113.xmin            = 0
  #plot1113.xmax            = 1
  plot1113.ymin            = 0
  plot1113.ymax            = 20
  #plot1113.lpos = "bottom-center"
  plot1113.name            = "CHF2ndJet_presel_Eta1stJetBump"
  plot1113.addZUncBand     = zUncBand
  plot1113.makeRatio       = makeRatio
  plot1113.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NHF2ndJet_highMej_Eta1stJetBump ---
  variableName = "h1_NHF2ndJet_Eta1stJetBump"

  plot1114 = Plot()
  ## inputs for stacked histograms
  plot1114.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1114.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1114.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1114.keys            = keys
  plot1114.xtit            = "Neutral hadron energy fraction 2^{nd} jet - Eta1stJetBump"
  plot1114.ytit            = "Number of events"
  plot1114.ylog            = "no"
  plot1114.rebin           = 10
  #plot1114.xmin            = 0
  #plot1114.xmax            = 1
  plot1114.ymin            = 0
  plot1114.ymax            = 40
  #plot1114.lpos = "bottom-center"
  plot1114.name            = "NHF2ndJet_presel_Eta1stJetBump"
  plot1114.addZUncBand     = zUncBand
  plot1114.makeRatio       = makeRatio
  plot1114.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CEF2ndJet_highMej_Eta1stJetBump ---
  variableName = "h1_CEF2ndJet_Eta1stJetBump"

  plot1115 = Plot()
  ## inputs for stacked histograms
  plot1115.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1115.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1115.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1115.keys            = keys
  plot1115.xtit            = "Charged EM energy fraction 2^{nd} jet - Eta1stJetBump"
  plot1115.ytit            = "Number of events"
  plot1115.ylog            = "no"
  plot1115.rebin           = 10
  #plot1115.xmin            = 0
  #plot1115.xmax            = 1
  plot1115.ymin            = 0
  plot1115.ymax            = 60
  #plot1115.lpos = "bottom-center"
  plot1115.name            = "CEF2ndJet_presel_Eta1stJetBump"
  plot1115.addZUncBand     = zUncBand
  plot1115.makeRatio       = makeRatio
  plot1115.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NEF2ndJet_highMej_Eta1stJetBump ---
  variableName = "h1_NEF2ndJet_Eta1stJetBump"

  plot1116 = Plot()
  ## inputs for stacked histograms
  plot1116.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1116.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1116.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1116.keys            = keys
  plot1116.xtit            = "Neutral EM energy fraction 2^{nd} jet - Eta1stJetBump"
  plot1116.ytit            = "Number of events"
  plot1116.ylog            = "no"
  plot1116.rebin           = 10
  #plot1116.xmin            = 0
  #plot1116.xmax            = 1
  plot1116.ymin            = 0
  plot1116.ymax            = 20
  #plot1116.lpos = "bottom-center"
  plot1116.name            = "NEF2ndJet_presel_Eta1stJetBump"
  plot1116.addZUncBand     = zUncBand
  plot1116.makeRatio       = makeRatio
  plot1116.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NCH2ndJet_highMej_Eta1stJetBump ---
  variableName = "h1_NCH2ndJet_Eta1stJetBump"

  plot1117 = Plot()
  ## inputs for stacked histograms
  plot1117.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1117.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1117.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1117.keys            = keys
  plot1117.xtit            = "Charged multiplicity 2^{nd} jet - Eta1stJetBump"
  plot1117.ytit            = "Number of events"
  plot1117.ylog            = "no"
  plot1117.rebin           = 10
  #plot1117.xmin            = -0.5
  #plot1117.xmax            = 99.5
  plot1117.ymin            = 0
  plot1117.ymax            = 40
  #plot1117.lpos = "bottom-center"
  plot1117.name            = "NCH2ndJet_presel_Eta1stJetBump"
  plot1117.addZUncBand     = zUncBand
  plot1117.makeRatio       = makeRatio
  plot1117.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NN2ndJet_highMej_Eta1stJetBump ---
  variableName = "h1_NN2ndJet_Eta1stJetBump"

  plot1118 = Plot()
  ## inputs for stacked histograms
  plot1118.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1118.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1118.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1118.keys            = keys
  plot1118.xtit            = "Neutral multiplicity 2^{nd} jet - Eta1stJetBump"
  plot1118.ytit            = "Number of events"
  plot1118.ylog            = "no"
  plot1118.rebin           = 10
  #plot1118.xmin            = -0.5
  #plot1118.xmax            = 99.5
  plot1118.ymin            = 0
  plot1118.ymax            = 40
  #plot1118.lpos = "bottom-center"
  plot1118.name            = "NN2ndJet_presel_Eta1stJetBump"
  plot1118.addZUncBand     = zUncBand
  plot1118.makeRatio       = makeRatio
  plot1118.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NC2ndJet_highMej_Eta1stJetBump ---
  variableName = "h1_NC2ndJet_Eta1stJetBump"

  plot1119 = Plot()
  ## inputs for stacked histograms
  plot1119.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1119.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1119.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1119.keys            = keys
  plot1119.xtit            = "Number of constituents 2^{nd} jet - Eta1stJetBump"
  plot1119.ytit            = "Number of events"
  plot1119.ylog            = "no"
  plot1119.rebin           = 10
  #plot1119.xmin            = -0.5
  #plot1119.xmax            = 99.5
  plot1119.ymin            = 0
  plot1119.ymax            = 30
  #plot1119.lpos = "bottom-center"
  plot1119.name            = "NC2ndJet_presel_Eta1stJetBump"
  plot1119.addZUncBand     = zUncBand
  plot1119.makeRatio       = makeRatio
  plot1119.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_E1stEle_PAS_Eta1stJetBump ---
  variableName = "h1_E1stEle_PAS_Eta1stJetBump"

  plot1120 = Plot()
  ## inputs for stacked histograms
  plot1120.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1120.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1120.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1120.keys            = keys
  plot1120.xtit            = "Energy 1^{st} electron - Eta1stJetBump"
  plot1120.ytit            = "Number of events"
  plot1120.ylog            = "no"
  plot1120.rebin           = "var"
  #plot1120.xmin            = 0
  #plot1120.xmax            = 1000
  plot1120.ymin            = 0
  plot1120.ymax            = 20
  #plot1120.lpos = "bottom-center"
  plot1120.name            = "E1stEle_presel_Eta1stJetBump"
  plot1120.addZUncBand     = zUncBand
  plot1120.makeRatio       = makeRatio
  plot1120.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot1120.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt1stEle_PAS_Eta1stJetBump ---
  variableName = "h1_Pt1stEle_PAS_Eta1stJetBump"

  plot1121 = Plot()
  ## inputs for stacked histograms
  plot1121.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1121.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1121.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1121.keys            = keys
  plot1121.xtit            = "p_{T} 1^{st} electron - Eta1stJetBump"
  plot1121.ytit            = "Number of events"
  plot1121.ylog            = "no"
  plot1121.rebin           = "var"
  #plot1121.xmin            = 0
  #plot1121.xmax            = 1000
  plot1121.ymin            = 0
  plot1121.ymax            = 30
  #plot1121.lpos = "bottom-center"
  plot1121.name            = "Pt1stEle_presel_Eta1stJetBump"
  plot1121.addZUncBand     = zUncBand
  plot1121.makeRatio       = makeRatio
  plot1121.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot1121.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stEle_PAS_Eta1stJetBump ---
  variableName = "h1_Eta1stEle_PAS_Eta1stJetBump"

  plot1122 = Plot()
  ## inputs for stacked histograms
  plot1122.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1122.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1122.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1122.keys            = keys
  plot1122.xtit            = "#eta 1^{st} electron - Eta1stJetBump"
  plot1122.ytit            = "Number of events"
  plot1122.ylog            = "no"
  plot1122.rebin           = 5
  #plot1122.xmin            = -5
  #plot1122.xmax            = 5
  plot1122.ymin            = 0
  plot1122.ymax            = 15
  #plot1122.lpos = "bottom-center"
  plot1122.name            = "Eta1stEle_presel_Eta1stJetBump"
  plot1122.addZUncBand     = zUncBand
  plot1122.makeRatio       = makeRatio
  plot1122.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stEle_PAS_Eta1stJetBump ---
  variableName = "h1_Phi1stEle_PAS_Eta1stJetBump"

  plot1123 = Plot()
  ## inputs for stacked histograms
  plot1123.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1123.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1123.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1123.keys            = keys
  plot1123.xtit            = "#phi 1^{st} electron - Eta1stJetBump"
  plot1123.ytit            = "Number of events"
  plot1123.ylog            = "no"
  plot1123.rebin           = 10
  #plot1123.xmin            = -3.1416
  #plot1123.xmax            = 3.1416
  plot1123.ymin            = 0
  plot1123.ymax            = 20
  #plot1123.lpos = "bottom-center"
  plot1123.name            = "Phi1stEle_presel_Eta1stJetBump"
  plot1123.addZUncBand     = zUncBand
  plot1123.makeRatio       = makeRatio
  plot1123.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Charge1stEle_PAS_Eta1stJetBump ---
  variableName = "h1_Charge1stEle_PAS_Eta1stJetBump"

  plot1124 = Plot()
  ## inputs for stacked histograms
  plot1124.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1124.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1124.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1124.keys            = keys
  plot1124.xtit            = "Charge 1^{st} electron - Eta1stJetBump"
  plot1124.ytit            = "Number of events"
  plot1124.ylog            = "no"
  #plot1124.rebin           = 1
  #plot1124.xmin            = -1.001
  #plot1124.xmax            = 1.001
  plot1124.ymin            = 0
  plot1124.ymax            = 50
  #plot1124.lpos = "bottom-center"
  plot1124.name            = "Charge1stEle_presel_Eta1stJetBump"
  plot1124.addZUncBand     = zUncBand
  plot1124.makeRatio       = makeRatio
  plot1124.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Conversion1stEle_Eta1stJetBump ---
  variableName = "h1_Conversion1stEle_Eta1stJetBump"

  plot1125 = Plot()
  ## inputs for stacked histograms
  plot1125.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot1125.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1125.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1125.keys            = keys
  plot1125.xtit            = "Conversion flag 1^{st} electron - Eta1stJetBump"
  plot1125.ytit            = "Number of events"
  plot1125.ylog            = "no"
  #plot1125.rebin           = 1
  #plot1125.xmin            = -0.5
  #plot1125.xmax            = 1.5
  plot1125.ymin            = 0
  plot1125.ymax            = 60
  #plot1125.lpos = "bottom-center"
  plot1125.name            = "Conversion1stEle_presel_Eta1stJetBump"
  plot1125.addZUncBand     = zUncBand
  plot1125.makeRatio       = makeRatio
  plot1125.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMETEle_Eta1stJetBump ---
  variableName = "h1_mDeltaPhiMETEle_Eta1stJetBump"

  plot1126 = Plot()
  ## inputs for stacked histograms
  plot1126.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1126.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1126.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1126.keys            = keys
  plot1126.xtit            = "#Delta#phi(MET,1^{st} ele) - Eta1stJetBump"
  plot1126.ytit            = "Number of events"
  plot1126.ylog            = "no"
  plot1126.rebin           = 10
  #plot1126.xmin            = 0
  #plot1126.xmax            = 3.1416
  plot1126.ymin            = 0
  plot1126.ymax            = 25
  #plot1126.lpos = "bottom-center"
  plot1126.name            = "mDeltaPhiMETEle_presel_Eta1stJetBump"
  plot1126.addZUncBand     = zUncBand
  plot1126.makeRatio       = makeRatio
  plot1126.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle1stJet_PAS_Eta1stJetBump ---
  variableName = "h1_mDeltaPhiEle1stJet_PAS_Eta1stJetBump"

  plot1127 = Plot()
  ## inputs for stacked histograms
  plot1127.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1127.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1127.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1127.keys            = keys
  plot1127.xtit            = "#Delta#phi(1^{st} ele,1^{st} jet) - Eta1stJetBump"
  plot1127.ytit            = "Number of events"
  plot1127.ylog            = "no"
  plot1127.rebin           = 10
  #plot1127.xmin            = 0
  #plot1127.xmax            = 3.1416
  plot1127.ymin            = 0
  plot1127.ymax            = 27
  #plot1127.lpos = "bottom-center"
  plot1127.name            = "mDeltaPhiEle1stJet_presel_Eta1stJetBump"
  plot1127.addZUncBand     = zUncBand
  plot1127.makeRatio       = makeRatio
  plot1127.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle2ndJet_PAS_Eta1stJetBump ---
  variableName = "h1_mDeltaPhiEle2ndJet_PAS_Eta1stJetBump"

  plot1128 = Plot()
  ## inputs for stacked histograms
  plot1128.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1128.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1128.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1128.keys            = keys
  plot1128.xtit            = "#Delta#phi(1^{st} ele,2^{nd} jet) - Eta1stJetBump"
  plot1128.ytit            = "Number of events"
  plot1128.ylog            = "no"
  plot1128.rebin           = 10
  #plot1128.xmin            = 0
  #plot1128.xmax            = 3.1416
  plot1128.ymin            = 0
  plot1128.ymax            = 20
  #plot1128.lpos = "bottom-center"
  plot1128.name            = "mDeltaPhiEle2ndJet_presel_Eta1stJetBump"
  plot1128.addZUncBand     = zUncBand
  plot1128.makeRatio       = makeRatio
  plot1128.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_Eta1stJetBump ---
  variableName = "h1_mDeltaPhiMET2ndJet_Eta1stJetBump"

  plot1129 = Plot()
  ## inputs for stacked histograms
  plot1129.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1129.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1129.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1129.keys            = keys
  plot1129.xtit            = "#Delta#phi(MET,2^{nd} jet) - Eta1stJetBump"
  plot1129.ytit            = "Number of events"
  plot1129.ylog            = "no"
  plot1129.rebin           = 10
  #plot1129.xmin            = 0
  #plot1129.xmax            = 3.1416
  plot1129.ymin            = 0
  plot1129.ymax            = 20
  #plot1129.lpos = "bottom-center"
  plot1129.name            = "mDeltaPhiMET2ndJet_presel_Eta1stJetBump"
  plot1129.addZUncBand     = zUncBand
  plot1129.makeRatio       = makeRatio
  plot1129.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Ptenu_PAS_Eta1stJetBump ---
  variableName = "h1_Ptenu_PAS_Eta1stJetBump"

  plot1130 = Plot()
  ## inputs for stacked histograms
  plot1130.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1130.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1130.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1130.keys            = keys
  plot1130.xtit            = "p_{T}(e#nu) - Eta1stJetBump"
  plot1130.ytit            = "Number of events"
  plot1130.ylog            = "no"
  plot1130.rebin           = "var"
  #plot1130.xmin            = 0
  #plot1130.xmax            = 2000
  plot1130.ymin            = 0
  plot1130.ymax            = 25
  #plot1130.lpos = "bottom-center"
  plot1130.name            = "Ptenu_presel_Eta1stJetBump"
  plot1130.addZUncBand     = zUncBand
  plot1130.makeRatio       = makeRatio
  plot1130.xbins           = [0,40,80,120,160,200,300,600]
  plot1130.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_1stJet_PTOverPTPlusMET_Eta1stJetBump ---
  variableName = "h1_1stJet_PTOverPTPlusMET_Eta1stJetBump"

  plot1131 = Plot()
  ## inputs for stacked histograms
  plot1131.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1131.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1131.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1131.keys            = keys
  plot1131.xtit            = "Pt1stJet/(Pt1stJet+MET) - Eta1stJetBump"
  plot1131.ytit            = "Number of events"
  plot1131.ylog            = "no"
  plot1131.rebin           = 5
  #plot1131.xmin            = 0
  #plot1131.xmax            = 1
  plot1131.ymin            = 0
  plot1131.ymax            = 20
  #plot1131.lpos = "bottom-center"
  plot1131.name            = "1stJet_PTOverPTPlusMET_presel_Eta1stJetBump"
  plot1131.addZUncBand     = zUncBand
  plot1131.makeRatio       = makeRatio
  plot1131.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MET_PAS_Eta1stJetBump ---
  variableName = "h1_MET_PAS_Eta1stJetBump"

  plot1132 = Plot()
  ## inputs for stacked histograms
  plot1132.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1132.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1132.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1132.keys            = keys
  plot1132.xtit            = "pfMET - Eta1stJetBump"
  plot1132.ytit            = "Number of events"
  plot1132.ylog            = "no"
  plot1132.rebin           = "var"
  #plot1132.xmin            = 0
  #plot1132.xmax            = 1000
  plot1132.ymin            = 0
  plot1132.ymax            = 20
  #plot1132.lpos = "bottom-center"
  plot1132.name            = "MET_presel_Eta1stJetBump"
  plot1132.addZUncBand     = zUncBand
  plot1132.makeRatio       = makeRatio
  plot1132.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot1132.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_Eta1stJetBump ---
  variableName = "h1_Njet_Eta1stJetBump"

  plot1133 = Plot()
  ## inputs for stacked histograms
  plot1133.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1133.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1133.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1133.keys            = keys
  plot1133.xtit            = "Number of jets - Eta1stJetBump"
  plot1133.ytit            = "Number of events"
  plot1133.ylog            = "no"
  plot1133.rebin           = 1
  plot1133.xmin            = -0.5
  plot1133.xmax            = 11.5
  plot1133.ymin            = 0
  plot1133.ymax            = 30
  #plot1133.lpos = "bottom-center"
  plot1133.name            = "Njet_presel_Eta1stJetBump"
  plot1133.addZUncBand     = zUncBand
  plot1133.makeRatio       = makeRatio
  plot1133.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_Eta1stJetBump ---
  variableName = "h1_NjetTCHELBTag_Eta1stJetBump"

  plot1134 = Plot()
  ## inputs for stacked histograms
  plot1134.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1134.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1134.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1134.keys            = keys
  plot1134.xtit            = "Number of b-tagged jets (TCHEL) - Eta1stJetBump"
  plot1134.ytit            = "Number of events"
  plot1134.ylog            = "no"
  plot1134.rebin           = 1
  plot1134.xmin            = -0.5
  plot1134.xmax            = 11.5
  plot1134.ymin            = 0
  plot1134.ymax            = 25
  #plot1134.lpos = "bottom-center"
  plot1134.name            = "NjetTCHELBTag_presel_Eta1stJetBump"
  plot1134.addZUncBand     = zUncBand
  plot1134.makeRatio       = makeRatio
  plot1134.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_minDRej_Eta1stJetBump ---
  variableName = "h1_minDRej_Eta1stJetBump"

  plot1135 = Plot()
  ## inputs for stacked histograms
  plot1135.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1135.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1135.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1135.keys            = keys
  plot1135.xtit            = "min#DeltaR(e,j) - Eta1stJetBump"
  plot1135.ytit            = "Number of events"
  plot1135.ylog            = "no"
  plot1135.rebin           = 10
  #plot1135.xmin            = 0
  #plot1135.xmax            = 10
  plot1135.ymin            = 0
  plot1135.ymax            = 35
  #plot1135.lpos = "bottom-center"
  plot1135.name            = "minDRej_presel_Eta1stJetBump"
  plot1135.addZUncBand     = zUncBand
  plot1135.makeRatio       = makeRatio
  plot1135.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_maxDRej_Eta1stJetBump ---
  variableName = "h1_maxDRej_Eta1stJetBump"

  plot1136 = Plot()
  ## inputs for stacked histograms
  plot1136.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1136.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1136.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1136.keys            = keys
  plot1136.xtit            = "max#DeltaR(e,j) - Eta1stJetBump"
  plot1136.ytit            = "Number of events"
  plot1136.ylog            = "no"
  plot1136.rebin           = 10
  #plot1136.xmin            = 0
  #plot1136.xmax            = 10
  plot1136.ymin            = 0
  plot1136.ymax            = 50
  #plot1136.lpos = "bottom-center"
  plot1136.name            = "maxDRej_presel_Eta1stJetBump"
  plot1136.addZUncBand     = zUncBand
  plot1136.makeRatio       = makeRatio
  plot1136.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DRjets_Eta1stJetBump ---
  variableName = "h1_DRjets_Eta1stJetBump"

  plot1137 = Plot()
  ## inputs for stacked histograms
  plot1137.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1137.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1137.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1137.keys            = keys
  plot1137.xtit            = "#DeltaR(j1,j2) - Eta1stJetBump"
  plot1137.ytit            = "Number of events"
  plot1137.ylog            = "no"
  plot1137.rebin           = 10
  #plot1137.xmin            = 0
  #plot1137.xmax            = 10
  plot1137.ymin            = 0
  plot1137.ymax            = 30
  #plot1137.lpos = "bottom-center"
  plot1137.name            = "DRjets_presel_Eta1stJetBump"
  plot1137.addZUncBand     = zUncBand
  plot1137.makeRatio       = makeRatio
  plot1137.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Mass1stJet_Eta1stJetBump ---
  variableName = "h1_Mass1stJet_Eta1stJetBump"

  plot1138 = Plot()
  ## inputs for stacked histograms
  plot1138.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1138.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1138.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1138.keys            = keys
  plot1138.xtit            = "Mass 1^{st} jet - Eta1stJetBump"
  plot1138.ytit            = "Number of events"
  plot1138.ylog            = "no"
  plot1138.rebin           = 5
  #plot1138.xmin            = -0.5
  #plot1138.xmax            = 99.5
  plot1138.ymin            = 0
  plot1138.ymax            = 40
  #plot1138.lpos = "bottom-center"
  plot1138.name            = "Mass1stJet_presel_Eta1stJetBump"
  plot1138.addZUncBand     = zUncBand
  plot1138.makeRatio       = makeRatio
  plot1138.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Mass2ndJet_Eta1stJetBump ---
  variableName = "h1_Mass2ndJet_Eta1stJetBump"

  plot1142 = Plot()
  ## inputs for stacked histograms
  plot1142.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1142.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1142.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1142.keys            = keys
  plot1142.xtit            = "Mass 2^{nd} jet - Eta1stJetBump"
  plot1142.ytit            = "Number of events"
  plot1142.ylog            = "no"
  plot1142.rebin           = 5
  #plot1142.xmin            = -0.5
  #plot1142.xmax            = 99.5
  plot1142.ymin            = 0
  plot1142.ymax            = 40
  #plot1142.lpos = "bottom-center"
  plot1142.name            = "Mass2ndJet_presel_Eta1stJetBump"
  plot1142.addZUncBand     = zUncBand
  plot1142.makeRatio       = makeRatio
  plot1142.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CaloMET_PAS_Eta1stJetBump ---
  variableName = "h1_CaloMET_PAS_Eta1stJetBump"

  plot1143 = Plot()
  ## inputs for stacked histograms
  plot1143.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1143.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1143.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1143.keys            = keys
  plot1143.xtit            = "caloMET - Eta1stJetBump"
  plot1143.ytit            = "Number of events"
  plot1143.ylog            = "no"
  plot1143.rebin           = "var"
  #plot1143.xmin            = 0
  #plot1143.xmax            = 1000
  plot1143.ymin            = 0
  plot1143.ymax            = 20
  #plot1143.lpos = "bottom-center"
  plot1143.name            = "CaloMET_presel_Eta1stJetBump"
  plot1143.addZUncBand     = zUncBand
  plot1143.makeRatio       = makeRatio
  plot1143.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot1143.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_TCMET_PAS_Eta1stJetBump ---
  variableName = "h1_TCMET_PAS_Eta1stJetBump"

  plot1144 = Plot()
  ## inputs for stacked histograms
  plot1144.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1144.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1144.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1144.keys            = keys
  plot1144.xtit            = "tcMET - Eta1stJetBump"
  plot1144.ytit            = "Number of events"
  plot1144.ylog            = "no"
  plot1144.rebin           = "var"
  #plot1144.xmin            = 0
  #plot1144.xmax            = 1000
  plot1144.ymin            = 0
  plot1144.ymax            = 20
  #plot1144.lpos = "bottom-center"
  plot1144.name            = "TCMET_presel_Eta1stJetBump"
  plot1144.addZUncBand     = zUncBand
  plot1144.makeRatio       = makeRatio
  plot1144.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot1144.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

  #--- h1_MET_over_CaloMET_Eta1stJetBump ---
  variableName = "h1_MET_over_CaloMET_Eta1stJetBump"

  plot1145 = Plot()
  ## inputs for stacked histograms
  plot1145.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1145.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1145.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1145.keys            = keys
  plot1145.xtit            = "pfMET/CaloMET - Eta1stJetBump"
  plot1145.ytit            = "Number of events"
  plot1145.ylog            = "no"
  plot1145.rebin           = 5
  #plot1145.xmin            = 0
  #plot1145.xmax            = 10
  plot1145.ymin            = 0
  plot1145.ymax            = 50
  #plot1145.lpos = "bottom-center"
  plot1145.name            = "MET_over_CaloMET_presel_Eta1stJetBump"
  plot1145.addZUncBand     = zUncBand
  plot1145.makeRatio       = makeRatio
  plot1145.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

  #--- h1_MTenu_PAS_Eta1stJetBump ---
  variableName = "h1_MTenu_PAS_Eta1stJetBump"

  plot1146 = Plot()
  ## inputs for stacked histograms
  plot1146.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1146.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1146.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1146.keys            = keys
  plot1146.xtit            = "MTenu - Eta1stJetBump"
  plot1146.ytit            = "Number of events"
  plot1146.ylog            = "no"
  plot1146.rebin           = "var"
  #plot1146.xmin            = 0
  #plot1146.xmax            = 1000
  plot1146.ymin            = 0
  plot1146.ymax            = 20
  #plot1146.lpos = "bottom-center"
  plot1146.name            = "MTenu_presel_Eta1stJetBump"
  plot1146.addZUncBand     = zUncBand
  plot1146.makeRatio       = makeRatio
  plot1146.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot1146.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)




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


  #--- h1_CHF1stJet_OutsideEta1stJetBump ---
  variableName = "h1_CHF1stJet_OutsideEta1stJetBump"

  plot1203 = Plot()
  ## inputs for stacked histograms
  plot1203.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1203.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1203.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1203.keys            = keys
  plot1203.xtit            = "Charged hadron energy fraction 1^{st} jet - OutsideEta1stJetBump"
  plot1203.ytit            = "Number of events"
  plot1203.ylog            = "no"
  plot1203.rebin           = 10
  #plot1203.xmin            = 0
  #plot1203.xmax            = 1
  plot1203.ymin            = 0
  plot1203.ymax            = 60
  #plot1203.lpos = "bottom-center"
  plot1203.name            = "CHF1stJet_presel_OutsideEta1stJetBump"
  plot1203.addZUncBand     = zUncBand
  plot1203.makeRatio       = makeRatio
  plot1203.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NHF1stJet_OutsideEta1stJetBump ---
  variableName = "h1_NHF1stJet_OutsideEta1stJetBump"

  plot1204 = Plot()
  ## inputs for stacked histograms
  plot1204.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1204.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1204.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1204.keys            = keys
  plot1204.xtit            = "Neutral hadron energy fraction 1^{st} jet - OutsideEta1stJetBump"
  plot1204.ytit            = "Number of events"
  plot1204.ylog            = "no"
  plot1204.rebin           = 10
  #plot1204.xmin            = 0
  #plot1204.xmax            = 1
  plot1204.ymin            = 0
  plot1204.ymax            = 150
  #plot1204.lpos = "bottom-center"
  plot1204.name            = "NHF1stJet_presel_OutsideEta1stJetBump"
  plot1204.addZUncBand     = zUncBand
  plot1204.makeRatio       = makeRatio
  plot1204.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CEF1stJet_OutsideEta1stJetBump ---
  variableName = "h1_CEF1stJet_OutsideEta1stJetBump"

  plot1205 = Plot()
  ## inputs for stacked histograms
  plot1205.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1205.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1205.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1205.keys            = keys
  plot1205.xtit            = "Charged EM energy fraction 1^{st} jet - OutsideEta1stJetBump"
  plot1205.ytit            = "Number of events"
  plot1205.ylog            = "no"
  plot1205.rebin           = 10
  #plot1205.xmin            = 0
  #plot1205.xmax            = 1
  plot1205.ymin            = 0
  plot1205.ymax            = 200
  #plot1205.lpos = "bottom-center"
  plot1205.name            = "CEF1stJet_presel_OutsideEta1stJetBump"
  plot1205.addZUncBand     = zUncBand
  plot1205.makeRatio       = makeRatio
  plot1205.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NEF1stJet_OutsideEta1stJetBump ---
  variableName = "h1_NEF1stJet_OutsideEta1stJetBump"

  plot1206 = Plot()
  ## inputs for stacked histograms
  plot1206.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1206.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1206.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1206.keys            = keys
  plot1206.xtit            = "Neutral EM energy fraction 1^{st} jet - OutsideEta1stJetBump"
  plot1206.ytit            = "Number of events"
  plot1206.ylog            = "no"
  plot1206.rebin           = 10
  #plot1206.xmin            = 0
  #plot1206.xmax            = 1
  plot1206.ymin            = 0
  plot1206.ymax            = 60
  #plot1206.lpos = "bottom-center"
  plot1206.name            = "NEF1stJet_presel_OutsideEta1stJetBump"
  plot1206.addZUncBand     = zUncBand
  plot1206.makeRatio       = makeRatio
  plot1206.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NCH1stJet_OutsideEta1stJetBump ---
  variableName = "h1_NCH1stJet_OutsideEta1stJetBump"

  plot1207 = Plot()
  ## inputs for stacked histograms
  plot1207.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1207.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1207.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1207.keys            = keys
  plot1207.xtit            = "Charged multiplicity 1^{st} jet - OutsideEta1stJetBump"
  plot1207.ytit            = "Number of events"
  plot1207.ylog            = "no"
  plot1207.rebin           = 10
  #plot1207.xmin            = -0.5
  #plot1207.xmax            = 99.5
  plot1207.ymin            = 0
  plot1207.ymax            = 150
  #plot1207.lpos = "bottom-center"
  plot1207.name            = "NCH1stJet_presel_OutsideEta1stJetBump"
  plot1207.addZUncBand     = zUncBand
  plot1207.makeRatio       = makeRatio
  plot1207.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NN1stJet_OutsideEta1stJetBump ---
  variableName = "h1_NN1stJet_OutsideEta1stJetBump"

  plot1208 = Plot()
  ## inputs for stacked histograms
  plot1208.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1208.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1208.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1208.keys            = keys
  plot1208.xtit            = "Neutral multiplicity 1^{st} jet - OutsideEta1stJetBump"
  plot1208.ytit            = "Number of events"
  plot1208.ylog            = "no"
  plot1208.rebin           = 10
  #plot1208.xmin            = -0.5
  #plot1208.xmax            = 99.5
  plot1208.ymin            = 0
  plot1208.ymax            = 100
  #plot1208.lpos = "bottom-center"
  plot1208.name            = "NN1stJet_presel_OutsideEta1stJetBump"
  plot1208.addZUncBand     = zUncBand
  plot1208.makeRatio       = makeRatio
  plot1208.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NC1stJet_OutsideEta1stJetBump ---
  variableName = "h1_NC1stJet_OutsideEta1stJetBump"

  plot1209 = Plot()
  ## inputs for stacked histograms
  plot1209.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1209.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1209.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1209.keys            = keys
  plot1209.xtit            = "Number of constituents 1^{st} jet - OutsideEta1stJetBump"
  plot1209.ytit            = "Number of events"
  plot1209.ylog            = "no"
  plot1209.rebin           = 10
  #plot1209.xmin            = -0.5
  #plot1209.xmax            = 99.5
  plot1209.ymin            = 0
  plot1209.ymax            = 80
  #plot1209.lpos = "bottom-center"
  plot1209.name            = "NC1stJet_presel_OutsideEta1stJetBump"
  plot1209.addZUncBand     = zUncBand
  plot1209.makeRatio       = makeRatio
  plot1209.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt2ndJet_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Pt2ndJet_PAS_OutsideEta1stJetBump"

  plot1210 = Plot()
  ## inputs for stacked histograms
  plot1210.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1210.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1210.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1210.keys            = keys
  plot1210.xtit            = "p_{T} 2^{nd} jet - OutsideEta1stJetBump"
  plot1210.ytit            = "Number of events"
  plot1210.ylog            = "no"
  plot1210.rebin           = "var"
  #plot1210.xmin            = 0
  #plot1210.xmax            = 1000
  plot1210.ymin            = 0
  plot1210.ymax            = 60
  #plot1210.lpos = "bottom-center"
  plot1210.name            = "Pt2ndJet_presel_OutsideEta1stJetBump"
  plot1210.addZUncBand     = zUncBand
  plot1210.makeRatio       = makeRatio
  plot1210.xbins           = [0,20,40,60,80,100,120,160,200,300,400,600,1000]
  plot1210.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta2ndJet_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Eta2ndJet_PAS_OutsideEta1stJetBump"

  plot1211 = Plot()
  ## inputs for stacked histograms
  plot1211.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1211.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1211.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1211.keys            = keys
  plot1211.xtit            = "#eta 2^{nd} jet - OutsideEta1stJetBump"
  plot1211.ytit            = "Number of events"
  plot1211.ylog            = "no"
  plot1211.rebin           = 5
  #plot1211.xmin            = -5
  #plot1211.xmax            = 5
  plot1211.ymin            = 0
  plot1211.ymax            = 50
  #plot1211.lpos = "bottom-center"
  plot1211.name            = "Eta2ndJet_presel_OutsideEta1stJetBump"
  plot1211.addZUncBand     = zUncBand
  plot1211.makeRatio       = makeRatio
  plot1211.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi2ndJet_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Phi2ndJet_PAS_OutsideEta1stJetBump"

  plot1212 = Plot()
  ## inputs for stacked histograms
  plot1212.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1212.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1212.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1212.keys            = keys
  plot1212.xtit            = "#phi 2^{nd} jet - OutsideEta1stJetBump"
  plot1212.ytit            = "Number of events"
  plot1212.ylog            = "no"
  plot1212.rebin           = 10
  #plot1212.xmin            = -3.1416
  #plot1212.xmax            = 3.1416
  plot1212.ymin            = 0
  plot1212.ymax            = 50
  #plot1212.lpos = "bottom-center"
  plot1212.name            = "Phi2ndJet_presel_OutsideEta1stJetBump"
  plot1212.addZUncBand     = zUncBand
  plot1212.makeRatio       = makeRatio
  plot1212.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CHF2ndJet_highMej_OutsideEta1stJetBump ---
  variableName = "h1_CHF2ndJet_OutsideEta1stJetBump"

  plot1213 = Plot()
  ## inputs for stacked histograms
  plot1213.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1213.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1213.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1213.keys            = keys
  plot1213.xtit            = "Charged hadron energy fraction 2^{nd} jet - OutsideEta1stJetBump"
  plot1213.ytit            = "Number of events"
  plot1213.ylog            = "no"
  plot1213.rebin           = 10
  #plot1213.xmin            = 0
  #plot1213.xmax            = 1
  plot1213.ymin            = 0
  plot1213.ymax            = 80
  #plot1213.lpos = "bottom-center"
  plot1213.name            = "CHF2ndJet_presel_OutsideEta1stJetBump"
  plot1213.addZUncBand     = zUncBand
  plot1213.makeRatio       = makeRatio
  plot1213.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NHF2ndJet_highMej_OutsideEta1stJetBump ---
  variableName = "h1_NHF2ndJet_OutsideEta1stJetBump"

  plot1214 = Plot()
  ## inputs for stacked histograms
  plot1214.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1214.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1214.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1214.keys            = keys
  plot1214.xtit            = "Neutral hadron energy fraction 2^{nd} jet - OutsideEta1stJetBump"
  plot1214.ytit            = "Number of events"
  plot1214.ylog            = "no"
  plot1214.rebin           = 10
  #plot1214.xmin            = 0
  #plot1214.xmax            = 1
  plot1214.ymin            = 0
  plot1214.ymax            = 140
  #plot1214.lpos = "bottom-center"
  plot1214.name            = "NHF2ndJet_presel_OutsideEta1stJetBump"
  plot1214.addZUncBand     = zUncBand
  plot1214.makeRatio       = makeRatio
  plot1214.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CEF2ndJet_highMej_OutsideEta1stJetBump ---
  variableName = "h1_CEF2ndJet_OutsideEta1stJetBump"

  plot1215 = Plot()
  ## inputs for stacked histograms
  plot1215.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1215.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1215.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1215.keys            = keys
  plot1215.xtit            = "Charged EM energy fraction 2^{nd} jet - OutsideEta1stJetBump"
  plot1215.ytit            = "Number of events"
  plot1215.ylog            = "no"
  plot1215.rebin           = 10
  #plot1215.xmin            = 0
  #plot1215.xmax            = 1
  plot1215.ymin            = 0
  plot1215.ymax            = 180
  #plot1215.lpos = "bottom-center"
  plot1215.name            = "CEF2ndJet_presel_OutsideEta1stJetBump"
  plot1215.addZUncBand     = zUncBand
  plot1215.makeRatio       = makeRatio
  plot1215.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NEF2ndJet_highMej_OutsideEta1stJetBump ---
  variableName = "h1_NEF2ndJet_OutsideEta1stJetBump"

  plot1216 = Plot()
  ## inputs for stacked histograms
  plot1216.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1216.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1216.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1216.keys            = keys
  plot1216.xtit            = "Neutral EM energy fraction 2^{nd} jet - OutsideEta1stJetBump"
  plot1216.ytit            = "Number of events"
  plot1216.ylog            = "no"
  plot1216.rebin           = 10
  #plot1216.xmin            = 0
  #plot1216.xmax            = 1
  plot1216.ymin            = 0
  plot1216.ymax            = 60
  #plot1216.lpos = "bottom-center"
  plot1216.name            = "NEF2ndJet_presel_OutsideEta1stJetBump"
  plot1216.addZUncBand     = zUncBand
  plot1216.makeRatio       = makeRatio
  plot1216.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NCH2ndJet_highMej_OutsideEta1stJetBump ---
  variableName = "h1_NCH2ndJet_OutsideEta1stJetBump"

  plot1217 = Plot()
  ## inputs for stacked histograms
  plot1217.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1217.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1217.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1217.keys            = keys
  plot1217.xtit            = "Charged multiplicity 2^{nd} jet - OutsideEta1stJetBump"
  plot1217.ytit            = "Number of events"
  plot1217.ylog            = "no"
  plot1217.rebin           = 10
  #plot1217.xmin            = -0.5
  #plot1217.xmax            = 99.5
  plot1217.ymin            = 0
  plot1217.ymax            = 100
  #plot1217.lpos = "bottom-center"
  plot1217.name            = "NCH2ndJet_presel_OutsideEta1stJetBump"
  plot1217.addZUncBand     = zUncBand
  plot1217.makeRatio       = makeRatio
  plot1217.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NN2ndJet_highMej_OutsideEta1stJetBump ---
  variableName = "h1_NN2ndJet_OutsideEta1stJetBump"

  plot1218 = Plot()
  ## inputs for stacked histograms
  plot1218.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1218.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1218.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1218.keys            = keys
  plot1218.xtit            = "Neutral multiplicity 2^{nd} jet - OutsideEta1stJetBump"
  plot1218.ytit            = "Number of events"
  plot1218.ylog            = "no"
  plot1218.rebin           = 10
  #plot1218.xmin            = -0.5
  #plot1218.xmax            = 99.5
  plot1218.ymin            = 0
  plot1218.ymax            = 100
  #plot1218.lpos = "bottom-center"
  plot1218.name            = "NN2ndJet_presel_OutsideEta1stJetBump"
  plot1218.addZUncBand     = zUncBand
  plot1218.makeRatio       = makeRatio
  plot1218.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NC2ndJet_highMej_OutsideEta1stJetBump ---
  variableName = "h1_NC2ndJet_OutsideEta1stJetBump"

  plot1219 = Plot()
  ## inputs for stacked histograms
  plot1219.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1219.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1219.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1219.keys            = keys
  plot1219.xtit            = "Number of constituents 2^{nd} jet - OutsideEta1stJetBump"
  plot1219.ytit            = "Number of events"
  plot1219.ylog            = "no"
  plot1219.rebin           = 10
  #plot1219.xmin            = -0.5
  #plot1219.xmax            = 99.5
  plot1219.ymin            = 0
  plot1219.ymax            = 100
  #plot1219.lpos = "bottom-center"
  plot1219.name            = "NC2ndJet_presel_OutsideEta1stJetBump"
  plot1219.addZUncBand     = zUncBand
  plot1219.makeRatio       = makeRatio
  plot1219.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_E1stEle_PAS_OutsideEta1stJetBump ---
  variableName = "h1_E1stEle_PAS_OutsideEta1stJetBump"

  plot1220 = Plot()
  ## inputs for stacked histograms
  plot1220.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1220.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1220.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1220.keys            = keys
  plot1220.xtit            = "Energy 1^{st} electron - OutsideEta1stJetBump"
  plot1220.ytit            = "Number of events"
  plot1220.ylog            = "no"
  plot1220.rebin           = "var"
  #plot1220.xmin            = 0
  #plot1220.xmax            = 1000
  plot1220.ymin            = 0
  plot1220.ymax            = 60
  #plot1220.lpos = "bottom-center"
  plot1220.name            = "E1stEle_presel_OutsideEta1stJetBump"
  plot1220.addZUncBand     = zUncBand
  plot1220.makeRatio       = makeRatio
  plot1220.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot1220.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Pt1stEle_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Pt1stEle_PAS_OutsideEta1stJetBump"

  plot1221 = Plot()
  ## inputs for stacked histograms
  plot1221.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1221.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1221.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1221.keys            = keys
  plot1221.xtit            = "p_{T} 1^{st} electron - OutsideEta1stJetBump"
  plot1221.ytit            = "Number of events"
  plot1221.ylog            = "no"
  plot1221.rebin           = "var"
  #plot1221.xmin            = 0
  #plot1221.xmax            = 1000
  plot1221.ymin            = 0
  plot1221.ymax            = 80
  #plot1221.lpos = "bottom-center"
  plot1221.name            = "Pt1stEle_presel_OutsideEta1stJetBump"
  plot1221.addZUncBand     = zUncBand
  plot1221.makeRatio       = makeRatio
  plot1221.xbins           = [0,20,40,60,80,100,120,160,200,300,500]
  plot1221.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Eta1stEle_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Eta1stEle_PAS_OutsideEta1stJetBump"

  plot1222 = Plot()
  ## inputs for stacked histograms
  plot1222.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1222.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1222.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1222.keys            = keys
  plot1222.xtit            = "#eta 1^{st} electron - OutsideEta1stJetBump"
  plot1222.ytit            = "Number of events"
  plot1222.ylog            = "no"
  plot1222.rebin           = 5
  #plot1222.xmin            = -5
  #plot1222.xmax            = 5
  plot1222.ymin            = 0
  plot1222.ymax            = 50
  #plot1222.lpos = "bottom-center"
  plot1222.name            = "Eta1stEle_presel_OutsideEta1stJetBump"
  plot1222.addZUncBand     = zUncBand
  plot1222.makeRatio       = makeRatio
  plot1222.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Phi1stEle_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Phi1stEle_PAS_OutsideEta1stJetBump"

  plot1223 = Plot()
  ## inputs for stacked histograms
  plot1223.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1223.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1223.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1223.keys            = keys
  plot1223.xtit            = "#phi 1^{st} electron - OutsideEta1stJetBump"
  plot1223.ytit            = "Number of events"
  plot1223.ylog            = "no"
  plot1223.rebin           = 10
  #plot1223.xmin            = -3.1416
  #plot1223.xmax            = 3.1416
  plot1223.ymin            = 0
  plot1223.ymax            = 60
  #plot1223.lpos = "bottom-center"
  plot1223.name            = "Phi1stEle_presel_OutsideEta1stJetBump"
  plot1223.addZUncBand     = zUncBand
  plot1223.makeRatio       = makeRatio
  plot1223.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Charge1stEle_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Charge1stEle_PAS_OutsideEta1stJetBump"

  plot1224 = Plot()
  ## inputs for stacked histograms
  plot1224.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1224.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1224.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1224.keys            = keys
  plot1224.xtit            = "Charge 1^{st} electron - OutsideEta1stJetBump"
  plot1224.ytit            = "Number of events"
  plot1224.ylog            = "no"
  #plot1224.rebin           = 1
  #plot1224.xmin            = -1.001
  #plot1224.xmax            = 1.001
  plot1224.ymin            = 0
  plot1224.ymax            = 140
  #plot1224.lpos = "bottom-center"
  plot1224.name            = "Charge1stEle_presel_OutsideEta1stJetBump"
  plot1224.addZUncBand     = zUncBand
  plot1224.makeRatio       = makeRatio
  plot1224.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Conversion1stEle_OutsideEta1stJetBump ---
  variableName = "h1_Conversion1stEle_OutsideEta1stJetBump"

  plot1225 = Plot()
  ## inputs for stacked histograms
  plot1225.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  keysStacknoQCD = copy.deepcopy(keysStack)
  keysStacknoQCD.remove("QCD")
  plot1225.keysStack       = keysStacknoQCD
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1225.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1225.keys            = keys
  plot1225.xtit            = "Conversion flag 1^{st} electron - OutsideEta1stJetBump"
  plot1225.ytit            = "Number of events"
  plot1225.ylog            = "no"
  #plot1225.rebin           = 1
  #plot1225.xmin            = -0.5
  #plot1225.xmax            = 1.5
  plot1225.ymin            = 0
  plot1225.ymax            = 200
  #plot1225.lpos = "bottom-center"
  plot1225.name            = "Conversion1stEle_presel_OutsideEta1stJetBump"
  plot1225.addZUncBand     = zUncBand
  plot1225.makeRatio       = makeRatio
  plot1225.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMETEle_OutsideEta1stJetBump ---
  variableName = "h1_mDeltaPhiMETEle_OutsideEta1stJetBump"

  plot1226 = Plot()
  ## inputs for stacked histograms
  plot1226.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1226.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1226.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1226.keys            = keys
  plot1226.xtit            = "#Delta#phi(MET,1^{st} ele) - OutsideEta1stJetBump"
  plot1226.ytit            = "Number of events"
  plot1226.ylog            = "no"
  plot1226.rebin           = 10
  #plot1226.xmin            = 0
  #plot1226.xmax            = 3.1416
  plot1226.ymin            = 0
  plot1226.ymax            = 60
  #plot1226.lpos = "bottom-center"
  plot1226.name            = "mDeltaPhiMETEle_presel_OutsideEta1stJetBump"
  plot1226.addZUncBand     = zUncBand
  plot1226.makeRatio       = makeRatio
  plot1226.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle1stJet_PAS_OutsideEta1stJetBump ---
  variableName = "h1_mDeltaPhiEle1stJet_PAS_OutsideEta1stJetBump"

  plot1227 = Plot()
  ## inputs for stacked histograms
  plot1227.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1227.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1227.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1227.keys            = keys
  plot1227.xtit            = "#Delta#phi(1^{st} ele,1^{st} jet) - OutsideEta1stJetBump"
  plot1227.ytit            = "Number of events"
  plot1227.ylog            = "no"
  plot1227.rebin           = 10
  #plot1227.xmin            = 0
  #plot1227.xmax            = 3.1416
  plot1227.ymin            = 0
  plot1227.ymax            = 60
  #plot1227.lpos = "bottom-center"
  plot1227.name            = "mDeltaPhiEle1stJet_presel_OutsideEta1stJetBump"
  plot1227.addZUncBand     = zUncBand
  plot1227.makeRatio       = makeRatio
  plot1227.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiEle2ndJet_PAS_OutsideEta1stJetBump ---
  variableName = "h1_mDeltaPhiEle2ndJet_PAS_OutsideEta1stJetBump"

  plot1228 = Plot()
  ## inputs for stacked histograms
  plot1228.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1228.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1228.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1228.keys            = keys
  plot1228.xtit            = "#Delta#phi(1^{st} ele,2^{nd} jet) - OutsideEta1stJetBump"
  plot1228.ytit            = "Number of events"
  plot1228.ylog            = "no"
  plot1228.rebin           = 10
  #plot1228.xmin            = 0
  #plot1228.xmax            = 3.1416
  plot1228.ymin            = 0
  plot1228.ymax            = 50
  #plot1228.lpos = "bottom-center"
  plot1228.name            = "mDeltaPhiEle2ndJet_presel_OutsideEta1stJetBump"
  plot1228.addZUncBand     = zUncBand
  plot1228.makeRatio       = makeRatio
  plot1228.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_mDeltaPhiMET2ndJet_OutsideEta1stJetBump ---
  variableName = "h1_mDeltaPhiMET2ndJet_OutsideEta1stJetBump"

  plot1229 = Plot()
  ## inputs for stacked histograms
  plot1229.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1229.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1229.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1229.keys            = keys
  plot1229.xtit            = "#Delta#phi(MET,2^{nd} jet) - OutsideEta1stJetBump"
  plot1229.ytit            = "Number of events"
  plot1229.ylog            = "no"
  plot1229.rebin           = 10
  #plot1229.xmin            = 0
  #plot1229.xmax            = 3.1416
  plot1229.ymin            = 0
  plot1229.ymax            = 50
  #plot1229.lpos = "bottom-center"
  plot1229.name            = "mDeltaPhiMET2ndJet_presel_OutsideEta1stJetBump"
  plot1229.addZUncBand     = zUncBand
  plot1229.makeRatio       = makeRatio
  plot1229.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Ptenu_PAS_OutsideEta1stJetBump ---
  variableName = "h1_Ptenu_PAS_OutsideEta1stJetBump"

  plot1230 = Plot()
  ## inputs for stacked histograms
  plot1230.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1230.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1230.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1230.keys            = keys
  plot1230.xtit            = "p_{T}(e#nu) - OutsideEta1stJetBump"
  plot1230.ytit            = "Number of events"
  plot1230.ylog            = "no"
  plot1230.rebin           = "var"
  #plot1230.xmin            = 0
  #plot1230.xmax            = 2000
  plot1230.ymin            = 0
  plot1230.ymax            = 70
  #plot1230.lpos = "bottom-center"
  plot1230.name            = "Ptenu_presel_OutsideEta1stJetBump"
  plot1230.addZUncBand     = zUncBand
  plot1230.makeRatio       = makeRatio
  plot1230.xbins           = [0,40,80,120,160,200,300,600]
  plot1230.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_1stJet_PTOverPTPlusMET_OutsideEta1stJetBump ---
  variableName = "h1_1stJet_PTOverPTPlusMET_OutsideEta1stJetBump"

  plot1231 = Plot()
  ## inputs for stacked histograms
  plot1231.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1231.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1231.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1231.keys            = keys
  plot1231.xtit            = "Pt1stJet/(Pt1stJet+MET) - OutsideEta1stJetBump"
  plot1231.ytit            = "Number of events"
  plot1231.ylog            = "no"
  plot1231.rebin           = 5
  #plot1231.xmin            = 0
  #plot1231.xmax            = 1
  plot1231.ymin            = 0
  plot1231.ymax            = 60
  #plot1231.lpos = "bottom-center"
  plot1231.name            = "1stJet_PTOverPTPlusMET_presel_OutsideEta1stJetBump"
  plot1231.addZUncBand     = zUncBand
  plot1231.makeRatio       = makeRatio
  plot1231.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_MET_PAS_OutsideEta1stJetBump ---
  variableName = "h1_MET_PAS_OutsideEta1stJetBump"

  plot1232 = Plot()
  ## inputs for stacked histograms
  plot1232.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1232.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1232.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1232.keys            = keys
  plot1232.xtit            = "pfMET - OutsideEta1stJetBump"
  plot1232.ytit            = "Number of events"
  plot1232.ylog            = "no"
  plot1232.rebin           = "var"
  #plot1232.xmin            = 0
  #plot1232.xmax            = 1000
  plot1232.ymin            = 0
  plot1232.ymax            = 70
  #plot1232.lpos = "bottom-center"
  plot1232.name            = "MET_presel_OutsideEta1stJetBump"
  plot1232.addZUncBand     = zUncBand
  plot1232.makeRatio       = makeRatio
  plot1232.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot1232.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Njet_OutsideEta1stJetBump ---
  variableName = "h1_Njet_OutsideEta1stJetBump"

  plot1233 = Plot()
  ## inputs for stacked histograms
  plot1233.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1233.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1233.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1233.keys            = keys
  plot1233.xtit            = "Number of jets - OutsideEta1stJetBump"
  plot1233.ytit            = "Number of events"
  plot1233.ylog            = "no"
  plot1233.rebin           = 1
  plot1233.xmin            = -0.5
  plot1233.xmax            = 11.5
  plot1233.ymin            = 0
  plot1233.ymax            = 100
  #plot1233.lpos = "bottom-center"
  plot1233.name            = "Njet_presel_OutsideEta1stJetBump"
  plot1233.addZUncBand     = zUncBand
  plot1233.makeRatio       = makeRatio
  plot1233.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_NjetTCHELBTag_OutsideEta1stJetBump ---
  variableName = "h1_NjetTCHELBTag_OutsideEta1stJetBump"

  plot1234 = Plot()
  ## inputs for stacked histograms
  plot1234.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1234.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1234.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1234.keys            = keys
  plot1234.xtit            = "Number of b-tagged jets (TCHEL) - OutsideEta1stJetBump"
  plot1234.ytit            = "Number of events"
  plot1234.ylog            = "no"
  plot1234.rebin           = 1
  plot1234.xmin            = -0.5
  plot1234.xmax            = 11.5
  plot1234.ymin            = 0
  plot1234.ymax            = 80
  #plot1234.lpos = "bottom-center"
  plot1234.name            = "NjetTCHELBTag_presel_OutsideEta1stJetBump"
  plot1234.addZUncBand     = zUncBand
  plot1234.makeRatio       = makeRatio
  plot1234.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_minDRej_OutsideEta1stJetBump ---
  variableName = "h1_minDRej_OutsideEta1stJetBump"

  plot1235 = Plot()
  ## inputs for stacked histograms
  plot1235.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1235.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1235.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1235.keys            = keys
  plot1235.xtit            = "min#DeltaR(e,j) - OutsideEta1stJetBump"
  plot1235.ytit            = "Number of events"
  plot1235.ylog            = "no"
  plot1235.rebin           = 10
  #plot1235.xmin            = 0
  #plot1235.xmax            = 10
  plot1235.ymin            = 0
  plot1235.ymax            = 100
  #plot1235.lpos = "bottom-center"
  plot1235.name            = "minDRej_presel_OutsideEta1stJetBump"
  plot1235.addZUncBand     = zUncBand
  plot1235.makeRatio       = makeRatio
  plot1235.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_maxDRej_OutsideEta1stJetBump ---
  variableName = "h1_maxDRej_OutsideEta1stJetBump"

  plot1236 = Plot()
  ## inputs for stacked histograms
  plot1236.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1236.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1236.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1236.keys            = keys
  plot1236.xtit            = "max#DeltaR(e,j) - OutsideEta1stJetBump"
  plot1236.ytit            = "Number of events"
  plot1236.ylog            = "no"
  plot1236.rebin           = 10
  #plot1236.xmin            = 0
  #plot1236.xmax            = 10
  plot1236.ymin            = 0
  plot1236.ymax            = 100
  #plot1236.lpos = "bottom-center"
  plot1236.name            = "maxDRej_presel_OutsideEta1stJetBump"
  plot1236.addZUncBand     = zUncBand
  plot1236.makeRatio       = makeRatio
  plot1236.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_DRjets_OutsideEta1stJetBump ---
  variableName = "h1_DRjets_OutsideEta1stJetBump"

  plot1237 = Plot()
  ## inputs for stacked histograms
  plot1237.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1237.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1237.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1237.keys            = keys
  plot1237.xtit            = "#DeltaR(j1,j2) - OutsideEta1stJetBump"
  plot1237.ytit            = "Number of events"
  plot1237.ylog            = "no"
  plot1237.rebin           = 10
  #plot1237.xmin            = 0
  #plot1237.xmax            = 10
  plot1237.ymin            = 0
  plot1237.ymax            = 100
  #plot1237.lpos = "bottom-center"
  plot1237.name            = "DRjets_presel_OutsideEta1stJetBump"
  plot1237.addZUncBand     = zUncBand
  plot1237.makeRatio       = makeRatio
  plot1237.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Mass1stJet_OutsideEta1stJetBump ---
  variableName = "h1_Mass1stJet_OutsideEta1stJetBump"

  plot1238 = Plot()
  ## inputs for stacked histograms
  plot1238.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1238.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1238.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1238.keys            = keys
  plot1238.xtit            = "Mass 1^{st} jet - OutsideEta1stJetBump"
  plot1238.ytit            = "Number of events"
  plot1238.ylog            = "no"
  plot1238.rebin           = 5
  #plot1238.xmin            = -0.5
  #plot1238.xmax            = 99.5
  plot1238.ymin            = 0
  plot1238.ymax            = 100
  #plot1238.lpos = "bottom-center"
  plot1238.name            = "Mass1stJet_presel_OutsideEta1stJetBump"
  plot1238.addZUncBand     = zUncBand
  plot1238.makeRatio       = makeRatio
  plot1238.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_Mass2ndJet_OutsideEta1stJetBump ---
  variableName = "h1_Mass2ndJet_OutsideEta1stJetBump"

  plot1242 = Plot()
  ## inputs for stacked histograms
  plot1242.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1242.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1242.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1242.keys            = keys
  plot1242.xtit            = "Mass 2^{nd} jet - OutsideEta1stJetBump"
  plot1242.ytit            = "Number of events"
  plot1242.ylog            = "no"
  plot1242.rebin           = 5
  #plot1242.xmin            = -0.5
  #plot1242.xmax            = 99.5
  plot1242.ymin            = 0
  plot1242.ymax            = 100
  #plot1242.lpos = "bottom-center"
  plot1242.name            = "Mass2ndJet_presel_OutsideEta1stJetBump"
  plot1242.addZUncBand     = zUncBand
  plot1242.makeRatio       = makeRatio
  plot1242.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_CaloMET_PAS_OutsideEta1stJetBump ---
  variableName = "h1_CaloMET_PAS_OutsideEta1stJetBump"

  plot1243 = Plot()
  ## inputs for stacked histograms
  plot1243.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1243.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1243.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1243.keys            = keys
  plot1243.xtit            = "caloMET - OutsideEta1stJetBump"
  plot1243.ytit            = "Number of events"
  plot1243.ylog            = "no"
  plot1243.rebin           = "var"
  #plot1243.xmin            = 0
  #plot1243.xmax            = 1000
  plot1243.ymin            = 0
  plot1243.ymax            = 70
  #plot1243.lpos = "bottom-center"
  plot1243.name            = "CaloMET_presel_OutsideEta1stJetBump"
  plot1243.addZUncBand     = zUncBand
  plot1243.makeRatio       = makeRatio
  plot1243.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot1243.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)


  #--- h1_TCMET_PAS_OutsideEta1stJetBump ---
  variableName = "h1_TCMET_PAS_OutsideEta1stJetBump"

  plot1244 = Plot()
  ## inputs for stacked histograms
  plot1244.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1244.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1244.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1244.keys            = keys
  plot1244.xtit            = "tcMET - OutsideEta1stJetBump"
  plot1244.ytit            = "Number of events"
  plot1244.ylog            = "no"
  plot1244.rebin           = "var"
  #plot1244.xmin            = 0
  #plot1244.xmax            = 1000
  plot1244.ymin            = 0
  plot1244.ymax            = 70
  #plot1244.lpos = "bottom-center"
  plot1244.name            = "TCMET_presel_OutsideEta1stJetBump"
  plot1244.addZUncBand     = zUncBand
  plot1244.makeRatio       = makeRatio
  plot1244.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot1244.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

  #--- h1_MET_over_CaloMET_OutsideEta1stJetBump ---
  variableName = "h1_MET_over_CaloMET_OutsideEta1stJetBump"

  plot1245 = Plot()
  ## inputs for stacked histograms
  plot1245.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1245.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1245.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1245.keys            = keys
  plot1245.xtit            = "pfMET/CaloMET - OutsideEta1stJetBump"
  plot1245.ytit            = "Number of events"
  plot1245.ylog            = "no"
  plot1245.rebin           = 5
  #plot1245.xmin            = 0
  #plot1245.xmax            = 10
  plot1245.ymin            = 0
  plot1245.ymax            = 100
  #plot1245.lpos = "bottom-center"
  plot1245.name            = "MET_over_CaloMET_presel_OutsideEta1stJetBump"
  plot1245.addZUncBand     = zUncBand
  plot1245.makeRatio       = makeRatio
  plot1245.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)

  #--- h1_MTenu_PAS_OutsideEta1stJetBump ---
  variableName = "h1_MTenu_PAS_OutsideEta1stJetBump"

  plot1246 = Plot()
  ## inputs for stacked histograms
  plot1246.histosStack     = generateHistoList( histoBaseName_userDef, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName_userDef, samplesForStackHistos, variableName, File_preselection)
  plot1246.keysStack       = keysStack
  ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
  plot1246.histos          = generateHistoList( histoBaseName_userDef, samplesForHistos, variableName, File_preselection)
  plot1246.keys            = keys
  plot1246.xtit            = "MTenu - OutsideEta1stJetBump"
  plot1246.ytit            = "Number of events"
  plot1246.ylog            = "no"
  plot1246.rebin           = "var"
  #plot1246.xmin            = 0
  #plot1246.xmax            = 1000
  plot1246.ymin            = 0
  plot1246.ymax            = 100
  #plot1246.lpos = "bottom-center"
  plot1246.name            = "MTenu_presel_OutsideEta1stJetBump"
  plot1246.addZUncBand     = zUncBand
  plot1246.makeRatio       = makeRatio
  plot1246.xbins           = [0,20,40,60,80,100,120,160,200,300,600]
  plot1246.histodata       = generateHisto( histoBaseName_userDef, sampleForDataHisto, variableName, File_preselection)



#-----------------------------------------------------------------------------------


# list of plots to be plotted
plots  = [plot0, plot1, plot2, plot_after2, plot3, plot4, plot5, plot_TCHELBTag, plot6, plot7, plot8, plot9, plot_TCHE1, plot_TCHE2
          , plot6and8, plot7and9
          , plot10, plot11, plot12, plot13, plot13_afterOtherDfCuts, plot14, plot14_ylog
          , plot14_plus, plot14_plus_ylog, plot14_minus, plot14_minus_ylog
          , plot14_0_1, plot14_1_2, plot14_2_pi
          , plot15, plot15_lep, plot15_jet, plot16
          , plot17, plot17_ylog, plot17_1, plot17_1_ylog, plot17_2, plot17_2_ylog, plot18
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

                ,plot700, plot701, plot702, plot703, plot704, plot705, plot706, plot707, plot708, plot709,     plot738, plot739, plot740, plot741,     plot710, plot711, plot712, plot713, plot714, plot715, plot716, plot717, plot718, plot719,    plot742,    plot720, plot721, plot722, plot723, plot724, plot725, plot726, plot727, plot728, plot729, plot730, plot731, plot732,     plot743, plot744, plot745,    plot733, plot734, plot735, plot736, plot737, plot746

                ,plot800, plot801, plot802, plot803, plot804, plot805, plot806, plot807, plot808, plot809,     plot838,                                plot810, plot811, plot812, plot813, plot814, plot815, plot816, plot817, plot818, plot819,    plot842,    plot820, plot821, plot822, plot823, plot824, plot825, plot826, plot827, plot828, plot829, plot830, plot831, plot832,     plot843, plot844, plot845,    plot833, plot834, plot835, plot836, plot837, plot846

                ,plot900, plot901, plot902, plot903, plot904, plot905, plot906, plot907, plot908, plot909,     plot938,                                plot910, plot911, plot912, plot913, plot914, plot915, plot916, plot917, plot918, plot919,    plot942,    plot920, plot921, plot922, plot923, plot924, plot925, plot926, plot927, plot928, plot929, plot930, plot931, plot932,     plot943, plot944, plot945,    plot933, plot934, plot935, plot936, plot937, plot946

                ,plot1000, plot1001, plot1002, plot1003, plot1004, plot1005, plot1006, plot1007, plot1008, plot1009,   plot1038,   plot1010, plot1011, plot1012, plot1013, plot1014, plot1015, plot1016, plot1017, plot1018, plot1019,   plot1042,   plot1020, plot1021, plot1022, plot1023, plot1024, plot1025, plot1026, plot1027, plot1028, plot1029, plot1030, plot1031, plot1032,    plot1043, plot1044, plot1045,    plot1033, plot1034, plot1035, plot1036, plot1037, plot1046

                ,plot1100, plot1102, plot1103, plot1104, plot1105, plot1106, plot1107, plot1108, plot1109,             plot1138,   plot1110, plot1111, plot1112, plot1113, plot1114, plot1115, plot1116, plot1117, plot1118, plot1119,   plot1142,   plot1120, plot1121, plot1122, plot1123, plot1124, plot1125, plot1126, plot1127, plot1128, plot1129, plot1130, plot1131, plot1132,    plot1143, plot1144, plot1145,    plot1133, plot1134, plot1135, plot1136, plot1137, plot1146

                ,plot1200, plot1202, plot1203, plot1204, plot1205, plot1206, plot1207, plot1208, plot1209,             plot1238,   plot1210, plot1211, plot1212, plot1213, plot1214, plot1215, plot1216, plot1217, plot1218, plot1219,   plot1242,   plot1220, plot1221, plot1222, plot1223, plot1224, plot1225, plot1226, plot1227, plot1228, plot1229, plot1230, plot1231, plot1232,    plot1243, plot1244, plot1245,    plot1233, plot1234, plot1235, plot1236, plot1237, plot1246
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

