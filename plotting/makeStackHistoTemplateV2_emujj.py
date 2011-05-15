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

def rebinHisto( histo, xmin, xmax, rebin, xbins ):
    new_histo = TH1F()
    if( xmin!="" and xmax!="" and rebin!="var" ):
        if(rebin!=""):
            histo.Rebin(rebin)
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
        #new_histo.SetBinContent( nbins, new_histo.GetBinContent(nbins) + overflowBinContent )
        #new_histo.SetBinError( nbins, sqrt( new_histo.GetBinError(nbins)**2 + overflowBinError2 ) )
    elif( xbins!="" and rebin=="var" ):
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
        #new_histo.SetBinContent( nbins, new_histo.GetBinContent(nbins) + overflowBinContent )
        #new_histo.SetBinError( nbins, sqrt( new_histo.GetBinError(nbins)**2 + overflowBinError2 ) )
    else:
        new_histo = histo
    return new_histo

def rebinHistos( histos, xmin, xmax, rebin, xbins ):
    new_histos = []
    for histo in histos:
        new_histo = TH1F()
        new_histo = rebinHisto( histo, xmin, xmax, rebin, xbins )
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

        self.histos = rebinHistos( self.histos, self.xmin, self.xmax, self.rebin, self.xbins )
        self.histosStack = rebinHistos( self.histosStack, self.xmin, self.xmax, self.rebin, self.xbins )
        self.histodata = rebinHisto( self.histodata, self.xmin, self.xmax, self.rebin, self.xbins )

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
            fPads1 = TPad("pad1", "", 0.00, 0.00, 0.99, 0.99)
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
        hsize=0.33
        vsize=0.33
        if (self.lpos=="bottom-center"):
            xstart=0.35
            ystart=0.25
        elif(self.lpos=="top-left"):
            xstart=0.12
            ystart=0.54
        else:
            xstart=0.54
            ystart=0.54
        legend = TLegend(xstart, ystart, xstart+hsize, ystart+vsize)
        legend.SetFillColor(kWhite)
        legend.SetMargin(0.2)
        legend.SetTextFont(132)
        legend.AddEntry(plot.histodata, "Data, 36.0 pb^{-1}","lp")

        #-- loop over histograms (stacked)
        Nstacked = len(self.histosStack)
        #stackColorIndexes = [20,38,14,45,20,38,14,45]
        stackColorIndexes = [2,4,3]
        stackFillStyleIds = [3354,3345,3395]
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
            legend.AddEntry(stack[iter], plot.keysStack[Nstacked - iter - 1],"lf")
            #draw stack
            if iter==0:
                stack[iter].SetTitle("")
                stack[iter].GetXaxis().SetTitle(self.xtit)
                stack[iter].GetXaxis().SetTitleFont(132)
                stack[iter].GetXaxis().SetTitleOffset(0.62)
                stack[iter].GetXaxis().SetLabelOffset(0.0)
                stack[iter].GetXaxis().SetTitleSize(0.07)
                stack[iter].GetXaxis().SetLabelSize(0.045)
                stack[iter].GetXaxis().SetLabelFont(132)
                stack[iter].GetYaxis().SetTitleFont(132)
                stack[iter].GetYaxis().SetTitleOffset(0.65)
                stack[iter].GetYaxis().SetTitleSize(0.07)
                stack[iter].GetYaxis().SetLabelSize(0.045)
                stack[iter].GetYaxis().SetLabelFont(132)
                #thisMin = stack[iter].GetXaxis().GetXmin()
                #thisMax = stack[iter].GetXaxis().GetXmax()
                #thisNbins = stack[iter].GetNbinsX()
                #newBinning = (thisMax - thisMin) / thisNbins
                #stack[iter].GetYaxis().SetTitle(self.ytit + " / ( "+ str(newBinning) + " )")
                stack[iter].GetYaxis().SetTitle(self.ytit + " / bin")
                if (self.ymin!="" and self.ymax!=""):
                    stack[iter].GetYaxis().SetLimits(self.ymin,self.ymax)
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
                stkcp.append(stack[iter].Clone())
            else:
                stkcp.append(stack[iter].Clone()) # this is the only way that I figured out to cover the previous histogram!
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
        dataColorIndexes = [1,1,1]
        dataLineIndexes = [2,1,3]
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
            #legend.AddEntry(self.histodata, "data","p")
            self.histodata.Draw("psame")

        #-- draw label
        l = TLatex()
        l.SetTextAlign(12)
        l.SetTextFont(132)
        l.SetTextSize(0.09)
        l.SetNDC()
#        l.DrawLatex(xstart,ystart-0.05,"CMS Preliminary 2010")
#        l.DrawLatex(xstart,ystart-0.10,"L_{int} = " + self.lint)
        if (self.lpos=="bottom-center"):
            l.DrawLatex(0.35,0.20,"CMS")
            #l.DrawLatex(0.35,0.15,"L_{int} = " + self.lint)
        if (self.lpos=="top-left"):
            l.DrawLatex(xstart+hsize+0.02,ystart+vsize-0.03,"CMS")
            #l.DrawLatex(xstart+hsize+0.02,ystart+vsize-0.08,"L_{int} = " + self.lint)
        else:
            l.DrawLatex(xstart-hsize+0.17,ystart+vsize-0.05,"CMS")
            #l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.08,"L_{int} = " + self.lint)

        legend.Draw("SAME")
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

#File_preselection = GetFile("$LQDATA/emujj/10.9pb-1/output_cutTable_emujjSample/analysisClass_emujjSample_plots.root")
File_preselection = GetFile("$LQDATA/emujj/36.0pb-1/output_cutTable_emujjSample/analysisClass_emujjSample_plots.root")

File_selection    = File_preselection

#### Common values for plots:
#otherBkgsKey="QCD, single top, VV+jets, W/W*+jets"
otherBkgsKey="Other Bkgs"
zUncBand="no"

pt_xmin=0
pt_xmax=500
pt_ymin=0.01
pt_ymax=50

eta_rebin=4
eta_ymin=0
eta_ymax=20



#--- Final plots are defined here

# Simply add or remove plots and update the list plots = [plot0, plot1, ...] below

histoBaseName = "histo1D__SAMPLE__cutHisto_allPreviousCuts________VARIABLE"

samplesForStackHistos = ["TTbar_Madgraph","ZJetAlpgen","OTHERBKG"]
keysStack =             ["ttbar", "Z/#gamma/Z* + jets", otherBkgsKey]

# samplesForHistos = ["LQenujj_M400"]
# keys             = ["LQ enujj M400"]
samplesForHistos = []
keys             = []

sampleForDataHisto = "DATA"


#--- Mee_TwoEleOnly ---
variableName = "Memu_OneEleOneMu"

# h_Mee_LQeejj_M100 = GetHisto("histo1D__LQeejj_M100__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_LQeejj_M200 = GetHisto("histo1D__LQeejj_M200__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# #h_Mee_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# #h_Mee_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()

plot0 = Plot()
## inputs for stacked histograms
plot0.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot0.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot0.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot0.keys            = keys
plot0.xtit            = "M(e#mu) (GeV/c^{2})"
plot0.ytit            = "Number of events"
# plot0.ylog            = "yes"
# plot0.rebin           = 1
# plot0.ymin            = 0.00000001
# plot0.ymax            = 20
plot0.ylog            = "no"
plot0.rebin           = 2
plot0.ymin            = 0
plot0.ymax            = 20
plot0.xmin            = 0
plot0.xmax            = 400
#plot0.lpos = "bottom-center"
plot0.name            = "Memu_allPreviousCuts_ylin"
plot0.addZUncBand     = zUncBand
plot0.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

plot0_ylog = Plot()
## inputs for stacked histograms
plot0_ylog.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot0_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot0_ylog.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot0_ylog.keys            = keys
plot0_ylog.xtit            = "M(e#mu) (GeV/c^{2})"
plot0_ylog.ytit            = "Number of events"
plot0_ylog.ylog            = "yes"
plot0_ylog.rebin           = 2
plot0_ylog.ymin            = 0.01
plot0_ylog.ymax            = 100
plot0_ylog.xmin            = 0
plot0_ylog.xmax            = 1000
#plot0_ylog.lpos = "bottom-center"
plot0_ylog.name            = "Memu_allPreviousCuts"
plot0_ylog.addZUncBand     = zUncBand
plot0_ylog.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- nEle_PtCut_IDISO_noOvrlp ---
variableName = "nEle_PtCut_IDISO_noOvrlp"

plot1 = Plot()
## inputs for stacked histograms
plot1.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot1.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot1.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot1.keys            = keys
plot1.xtit            = "Number of electrons"
plot1.ytit            = "Number of events"
plot1.ylog            = "yes"
plot1.rebin           = 1
plot1.ymin            = 0.0001
plot1.ymax            = 60000000
#plot1.lpos = "bottom-center"
plot1.name            = "nEle_allPreviousCuts"
plot1.addZUncBand     = zUncBand
plot1.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- nMu_PtCut_IDISO ---
variableName = "nMu_PtCut_IDISO"

plot1b = Plot()
## inputs for stacked histograms
plot1b.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot1b.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot1b.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot1b.keys            = keys
plot1b.xtit            = "Number of muons"
plot1b.ytit            = "Number of events"
plot1b.ylog            = "yes"
plot1b.rebin           = 1
plot1b.ymin            = 0.0001
plot1b.ymax            = 60000000
#plot1b.lpos = "bottom-center"
plot1b.name            = "nMu_allPreviousCuts"
plot1b.addZUncBand     = zUncBand
plot1b.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- nLept_PtCut_IDISO ---
variableName = "nLept_PtCut_IDISO"

plot1c = Plot()
## inputs for stacked histograms
plot1c.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot1c.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot1c.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot1c.keys            = keys
plot1c.xtit            = "Number of leptons"
plot1c.ytit            = "Number of events"
plot1c.ylog            = "yes"
plot1c.rebin           = 1
plot1c.ymin            = 0.0001
plot1c.ymax            = 1000
#plot1c.lpos = "bottom-center"
plot1c.name            = "nMu_allPreviousCuts"
plot1c.addZUncBand     = zUncBand
plot1c.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)




#--- Pt1stEle_PAS  ---
variableName = "Pt1stEle_PAS"

plot2 = Plot()
## inputs for stacked histograms
plot2.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot2.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot2.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot2.keys            = keys
plot2.xtit            = "pT 1st electron (GeV/c)"
plot2.ytit            = "Number of events"
plot2.ylog            = "yes"
plot2.rebin           = 1
plot2.xmin            = pt_xmin
plot2.xmax            = pt_xmax
plot2.ymin            = pt_ymin
plot2.ymax            = pt_ymax
#plot2.lpos = "bottom-center"
plot2.name            = "pT1stEle_allPreviousCuts"
plot2.addZUncBand     = zUncBand
plot2.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Eta1stEle_PAS  ---
variableName = "Eta1stEle_PAS"

plot3 = Plot()
## inputs for stacked histograms
plot3.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot3.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot3.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot3.keys            = keys
plot3.xtit            = "#eta 1st electron"
plot3.ytit            = "Number of events"
plot3.xmin            = -5
plot3.xmax            = 5
plot3.rebin           = eta_rebin
plot3.ymin            = eta_ymin
plot3.ymax            = eta_ymax
plot3.lpos = "top-left"
plot3.name            = "Eta1stEle_allPreviousCuts"
plot3.addZUncBand     = zUncBand
plot3.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)



#--- Pt2ndEle_PAS  ---
variableName = "Pt1stMu_PAS"

plot4 = Plot()
## inputs for stacked histograms
plot4.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot4.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot4.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot4.keys            = keys
plot4.xtit            = "pT 1st muon (GeV/c)"
plot4.ytit            = "Number of events"
plot4.ylog            = "yes"
plot4.rebin           = 1
plot4.xmin            = pt_xmin
plot4.xmax            = pt_xmax
plot4.ymin            = pt_ymin
plot4.ymax            = pt_ymax
#plot4.lpos = "bottom-center"
plot4.name            = "pT1stMu_allPreviousCuts"
plot4.addZUncBand     = zUncBand
plot4.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Eta2ndEle_PAS  ---
variableName = "Eta1stMu_PAS"

plot5 = Plot()
## inputs for stacked histograms
plot5.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot5.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot5.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot5.keys            = keys
plot5.xtit            = "#eta 1st muon"
plot5.ytit            = "Number of events"
plot5.xmin            = -5
plot5.xmax            = 5
plot5.rebin           = eta_rebin
plot5.ymin            = eta_ymin
plot5.ymax            = eta_ymax
plot5.lpos = "top-left"
plot5.name            = "Eta1stMu_allPreviousCuts"
plot5.addZUncBand     = zUncBand
plot5.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- nJet_TwoEleOnly_EtaCut ---
variableName = "nJet_OneEleOneMu_EtaCut"

plot6 = Plot()
## inputs for stacked histograms
plot6.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot6.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot6.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot6.keys            = keys
plot6.xtit            = "Number of jets"
plot6.ytit            = "Number of events"
plot6.ylog            = "yes"
plot6.rebin           = 1
plot6.xmin            = 0
plot6.xmax            = 12
plot6.ymin            = 0.01
plot6.ymax            = 100
#plot6.lpos = "bottom-center"
plot6.name            = "nJet_allPreviousCuts"
plot6.addZUncBand     = zUncBand
plot6.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)



#--- Pt1stJet_PAS ---
variableName = "Pt1stJet_PAS"

plot7 = Plot()
## inputs for stacked histograms
plot7.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot7.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot7.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot7.keys            = keys
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
plot7.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Eta1stJet_PAS ---
variableName = "Eta1stJet_PAS"

plot8 = Plot()
## inputs for stacked histograms
plot8.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot8.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot8.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot8.keys            = keys
plot8.xtit            = "#eta 1st jet"
plot8.ytit            = "Number of events"
plot8.xmin            = -5
plot8.xmax            = 5
plot8.rebin           = eta_rebin
plot8.ymin            = eta_ymin
plot8.ymax            = eta_ymax
plot8.lpos = "top-left"
plot8.name            = "Eta1stJet_allPreviousCuts"
plot8.addZUncBand     = zUncBand
plot8.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Pt2ndJet_PAS ---
variableName = "Pt2ndJet_PAS"

plot9 = Plot()
## inputs for stacked histograms
plot9.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot9.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot9.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot9.keys            = keys
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
plot9.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Eta2ndJet_PAS ---
variableName = "Eta2ndJet_PAS"

plot10 = Plot()
## inputs for stacked histograms
plot10.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot10.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot10.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot10.keys            = keys
plot10.xtit            = "#eta 2nd jet"
plot10.ytit            = "Number of events"
plot10.xmin            = -5
plot10.xmax            = 5
plot10.rebin           = eta_rebin
plot10.ymin            = eta_ymin
plot10.ymax            = eta_ymax
plot10.lpos = "top-left"
plot10.name            = "Eta2ndJet_allPreviousCuts"
plot10.addZUncBand     = zUncBand
plot10.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- sT ---
variableName = "sT_PAS"

plot11 = Plot()
## inputs for stacked histograms
plot11.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot11.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot11.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot11.keys            = keys
plot11.xtit            = "St (GeV/c)"
plot11.ytit            = "Number of events"
#plot11.xlog            = "yes"
plot11.ylog            = "yes"
plot11.rebin           = 2
plot11.xmin            = 0
plot11.xmax            = 1000
plot11.ymin            = 0.01
plot11.ymax            = 100
#plot11.lpos = "bottom-center"
plot11.name            = "sT_allPreviousCuts"
plot11.addZUncBand     = zUncBand
plot11.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- sTele ---
variableName = "sTemu_PAS"

plot11_ele = Plot()
## inputs for stacked histograms
plot11_ele.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot11_ele.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot11_ele.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot11_ele.keys            = keys
plot11_ele.xtit            = "St leptons (GeV/c)"
plot11_ele.ytit            = "Number of events"
#plot11_ele.xlog            = "yes"
plot11_ele.ylog            = "yes"
plot11_ele.rebin           = 2
plot11_ele.xmin            = 0
plot11_ele.xmax            = 1000
plot11_ele.ymin            = 0.01
plot11_ele.ymax            = 100
#plot11_ele.lpos = "bottom-center"
plot11_ele.name            = "sTemu_allPreviousCuts"
plot11_ele.addZUncBand     = zUncBand
plot11_ele.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- sTjet ---
variableName = "sTjet_PAS"

plot11_jet = Plot()
## inputs for stacked histograms
plot11_jet.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot11_jet.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot11_jet.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot11_jet.keys            = keys
plot11_jet.xtit            = "St jets (GeV/c)"
plot11_jet.ytit            = "Number of events"
#plot11_jet.xlog            = "yes"
plot11_jet.ylog            = "yes"
plot11_jet.rebin           = 2
plot11_jet.xmin            = 0
plot11_jet.xmax            = 1000
plot11_jet.ymin            = 0.01
plot11_jet.ymax            = 100
#plot11_jet.lpos = "bottom-center"
plot11_jet.name            = "sTjet_allPreviousCuts"
plot11_jet.addZUncBand     = zUncBand
plot11_jet.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)




##--- Mej preselection
variableNames = ["Mlj_1stPair_PAS", "Mlj_2ndPair_PAS"]

plot12 = Plot()
#plot12.histosStack     = [h_Mej_presel_TTbar, h_Mej_presel_ZJetAlpgen, h_Mej_presel_OTHERBKG]
plot12.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot12.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot12.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot12.keys            = keys
plot12.xtit            = "Mlj (GeV/c^{2})"
plot12.ytit            = "Number of events x 2"
plot12.ylog            = "yes"
plot12.rebin           = 2
plot12.xmin            = 0
plot12.xmax            = 1000
plot12.ymin            = 0.01
plot12.ymax            = 100
#plot12.lpos = "bottom-center"
plot12.name            = "Mlj_allPreviousCuts"
plot12.addZUncBand     = zUncBand
plot12.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)


#--- Mee_PAS (after preselection) ---
variableName = "Memu_PAS"

plot13 = Plot()
plot13.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot13.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot13.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot13.keys            = keys
plot13.xtit            = "M(e#mu) (GeV/c^{2})"
plot13.ytit            = "Number of events"
# plot13.ylog            = "yes"
# plot13.rebin           = 1
# plot13.ymin            = 0.00000001
# plot13.ymax            = 20
plot13.ylog            = "no"
plot13.rebin           = 2
plot13.ymin            = 0
plot13.ymax            = 10
plot13.xmin            = 0
plot13.xmax            = 500
#plot13.lpos = "bottom-center"
plot13.name            = "Memu_FullPreSel_allPreviousCuts_ylin"
plot13.addZUncBand     = zUncBand
plot13.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

plot13_ylog = Plot()
## inputs for stacked histograms
plot13_ylog.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot13_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot13_ylog.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot13_ylog.keys            = keys
plot13_ylog.xtit            = "M(e#mu) (GeV/c^{2})"
plot13_ylog.ytit            = "Number of events"
plot13_ylog.ylog            = "yes"
plot13_ylog.rebin           = 2
plot13_ylog.ymin            = 0.01
plot13_ylog.ymax            = 100
plot13_ylog.xmin            = 0
plot13_ylog.xmax            = 1000
#plot13_ylog.lpos = "bottom-center"
plot13_ylog.name            = "Memu_FullPreSel_allPreviousCuts"
plot13_ylog.addZUncBand     = zUncBand
plot13_ylog.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)



#--- Mjj_PAS (after preselection) ---
variableName = "Mjj_PAS"

plot14 = Plot()
## inputs for stacked histograms
## it created h_Mjj_FullPreSel_TTbar, h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen , h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen+h_Mjj_FullPreSel_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
plot14.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot14.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot14.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot14.keys            = keys
plot14.xtit            = "M(jj) (GeV/c^{2})"
plot14.ytit            = "Number of events"
# plot14.ylog            = "yes"
# plot14.rebin           = 1
# plot14.ymin            = 0.00000001
# plot14.ymax            = 20
plot14.ylog            = "no"
plot14.rebin           = 2
plot14.ymin            = 0
plot14.ymax            = 10
plot14.xmin            = 0
plot14.xmax            = 500
#plot14.lpos = "bottom-center"
plot14.name            = "Mjj_FullPreSel_allPreviousCuts_ylin"
plot14.addZUncBand     = zUncBand
plot14.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

plot14_ylog = Plot()
## inputs for stacked histograms
## it created h_Mjj_FullPreSel_TTbar, h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen , h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen+h_Mjj_FullPreSel_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
plot14_ylog.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot14_ylog.keysStack       = keysStack

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot14_ylog.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot14_ylog.keys            = keys
plot14_ylog.xtit            = "M(jj) (GeV/c^{2})"
plot14_ylog.ytit            = "Number of events"
plot14_ylog.ylog            = "yes"
plot14_ylog.rebin           = 2
plot14_ylog.ymin            = 0.01
plot14_ylog.ymax            = 100
plot14_ylog.xmin            = 0
plot14_ylog.xmax            = 1000
#plot14_ylog.lpos = "bottom-center"
plot14_ylog.name            = "Mjj_FullPreSel_allPreviousCuts"
plot14_ylog.addZUncBand     = zUncBand
plot14_ylog.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


##--- Pt leptons AllPreviousCuts ---
variableNames = ["Pt1stEle_PAS","Pt1stMu_PAS"]

plot2and4 = Plot()
## inputs for stacked histograms
plot2and4.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot2and4.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot2and4.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot2and4.keys            = keys
plot2and4.xtit            = "pT leptons (GeV/c)"
plot2and4.ytit            = "Number of events x 2"
plot2and4.ylog            = "yes"
plot2and4.rebin           = 1
plot2and4.xmin            = pt_xmin
plot2and4.xmax            = pt_xmax
plot2and4.ymin            = pt_ymin
plot2and4.ymax            = pt_ymax
#plot2and4.lpos = "bottom-center"
plot2and4.name            = "pTEleMu_allPreviousCuts"
plot2and4.addZUncBand     = zUncBand
plot2and4.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)


##--- Eta Leptons AllPreviousCuts ---
variableNames = ["Eta1stEle_PAS","Eta1stMu_PAS"]

plot3and5 = Plot()
## inputs for stacked histograms
plot3and5.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot3and5.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot3and5.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot3and5.keys            = keys
plot3and5.xtit            = "#eta leptons"
plot3and5.ytit            = "Number of events x 2"
plot3and5.xmin            = -5
plot3and5.xmax            = 5
plot3and5.rebin           = eta_rebin/2
plot3and5.ymin            = eta_ymin
plot3and5.ymax            = eta_ymax
plot3and5.lpos            = "top-left"
#plot3and5.lpos = "bottom-center"
plot3and5.name            = "etaEleMu_allPreviousCuts"
plot3and5.addZUncBand     = zUncBand
plot3and5.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)


##--- Pt Jets AllPreviousCuts ---
variableNames = ["Pt1stJet_PAS","Pt2ndJet_PAS"]

plot7and9 = Plot()
## inputs for stacked histograms
plot7and9.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot7and9.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot7and9.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot7and9.keys            = keys
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
plot7and9.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)

##--- Eta Jets AllPreviousCuts ---
variableNames = ["Eta1stJet_PAS","Eta2ndJet_PAS"]

plot8and10 = Plot()
## inputs for stacked histograms
plot8and10.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot8and10.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot8and10.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot8and10.keys            = keys
plot8and10.xtit            = "#eta jets"
plot8and10.ytit            = "Number of events x 2"
plot8and10.xmin            = -5
plot8and10.xmax            = 5
plot8and10.rebin           = eta_rebin/2
plot8and10.ymin            = eta_ymin
plot8and10.ymax            = eta_ymax
plot8and10.lpos            = "top-left"
#plot8and10.lpos = "bottom-center"
plot8and10.name            = "etaJets_allPreviousCuts"
plot8and10.addZUncBand     = zUncBand
plot8and10.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)


#--- pfMET ---
variableName = "pfMET_PAS"

plot15 = Plot()
## inputs for stacked histograms
plot15.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot15.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot15.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot15.keys            = keys
plot15.xtit            = "pfMET (GeV/c)"
plot15.ytit            = "Number of events"
#plot15.xlog            = "yes"
plot15.ylog            = "yes"
plot15.rebin           = 2
plot15.xmin            = 0
plot15.xmax            = 500
plot15.ymin            = 0.01
plot15.ymax            = 100
#plot15.lpos = "bottom-center"
plot15.name            = "pfMET_allPreviousCuts"
plot15.addZUncBand     = zUncBand
plot15.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- caloMET ---
variableName = "caloMET_PAS"

plot16 = Plot()
## inputs for stacked histograms
plot16.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot16.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot16.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot16.keys            = keys
plot16.xtit            = "caloMET (GeV/c)"
plot16.ytit            = "Number of events"
#plot16.xlog            = "yes"
plot16.ylog            = "yes"
plot16.rebin           = 2
plot16.xmin            = 0
plot16.xmax            = 500
plot16.ymin            = 0.01
plot16.ymax            = 100
#plot16.lpos = "bottom-center"
plot16.name            = "caloMET_allPreviousCuts"
plot16.addZUncBand     = zUncBand
plot16.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Pt1stEle_IDISO_NoOvrlp  ---
variableName = "Pt1stEle_IDISO_NoOvrlp"

plot2_nojet = Plot()
## inputs for stacked histograms
plot2_nojet.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot2_nojet.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot2_nojet.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot2_nojet.keys            = keys
plot2_nojet.xtit            = "pT 1st electron (GeV/c)"
plot2_nojet.ytit            = "Number of events"
plot2_nojet.ylog            = "yes"
plot2_nojet.rebin           = 1
plot2_nojet.xmin            = pt_xmin
plot2_nojet.xmax            = pt_xmax
plot2_nojet.ymin            = pt_ymin
plot2_nojet.ymax            = pt_ymax*10
#plot2_nojet.lpos = "bottom-center"
plot2_nojet.name            = "pT1stEle_nojet_allPreviousCuts"
plot2_nojet.addZUncBand     = zUncBand
plot2_nojet.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Eta1stEle_IDISO_NoOvrlp  ---
variableName = "Eta1stEle_IDISO_NoOvrlp"

plot3_nojet = Plot()
## inputs for stacked histograms
plot3_nojet.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot3_nojet.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot3_nojet.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot3_nojet.keys            = keys
plot3_nojet.xtit            = "#eta 1st electron"
plot3_nojet.ytit            = "Number of events"
plot3_nojet.xmin            = -5
plot3_nojet.xmax            = 5
plot3_nojet.rebin           = eta_rebin
plot3_nojet.ymin            = eta_ymin
plot3_nojet.ymax            = eta_ymax
plot3_nojet.lpos = "top-left"
plot3_nojet.name            = "Eta1stEle_nojet_allPreviousCuts"
plot3_nojet.addZUncBand     = zUncBand
plot3_nojet.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Pt2ndEle_IDISO_NoOvrlp  ---
variableName = "Pt1stMu_IDISO"

plot4_nojet = Plot()
## inputs for stacked histograms
plot4_nojet.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot4_nojet.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot4_nojet.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot4_nojet.keys            = keys
plot4_nojet.xtit            = "pT 1st muon (GeV/c)"
plot4_nojet.ytit            = "Number of events"
plot4_nojet.ylog            = "yes"
plot4_nojet.rebin           = 1
plot4_nojet.xmin            = pt_xmin
plot4_nojet.xmax            = pt_xmax
plot4_nojet.ymin            = pt_ymin
plot4_nojet.ymax            = pt_ymax*10
#plot4_nojet.lpos = "bottom-center"
plot4_nojet.name            = "pT1stMu_nojet_allPreviousCuts"
plot4_nojet.addZUncBand     = zUncBand
plot4_nojet.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Eta2ndEle_IDISO_NoOvrlp  ---
variableName = "Eta1stMu_IDISO"

plot5_nojet = Plot()
## inputs for stacked histograms
plot5_nojet.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot5_nojet.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot5_nojet.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot5_nojet.keys            = keys
plot5_nojet.xtit            = "#eta 1st muon"
plot5_nojet.ytit            = "Number of events"
plot5_nojet.xmin            = -5
plot5_nojet.xmax            = 5
plot5_nojet.rebin           = eta_rebin
plot5_nojet.ymin            = eta_ymin
plot5_nojet.ymax            = eta_ymax
plot5_nojet.lpos = "top-left"
plot5_nojet.name            = "Eta1stMu_nojet_allPreviousCuts"
plot5_nojet.addZUncBand     = zUncBand
plot5_nojet.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)



# ############################ Plots below to be done after full selection ######################

# ##--- sT AllOtherCuts ---
# variableName = "sT"

# plot20 = Plot()
# plot20.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_selection)
# plot20.keysStack       = keysStack
# ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot20.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_selection)
# plot20.keys            = keys
# plot20.xtit            = "St (GeV/c)"
# plot20.ytit            = "Number of events"
# plot20.ylog            = "yes"
# plot20.rebin           = 5
# plot20.xmin            = 0
# plot20.xmax            = 1000
# plot20.ymin            = 0.01
# plot20.ymax            = 100
# #plot20.lpos = "bottom-center"
# plot20.name            = "sT_allOtherCuts"
# plot20.addZUncBand     = zUncBand
# plot20.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_selection)


# ##--- Mlj AllOtherCuts ---
# variableNames = ["Mlj_1stPair","Mlj_2ndPair"]

# plot21 = Plot()
# plot21.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_selection)
# plot21.keysStack       = keysStack
# ## this is the list of histograms that should be simply overlaid on top of the stacked histogram
# plot21.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_selection)
# plot21.keys            = keys
# plot21.xtit            = "Mlj (GeV/c^{2})"
# plot21.ytit            = "Number of events x 2"
# plot21.ylog            = "yes"
# plot21.rebin           = 5
# plot21.xmin            = 0
# plot21.xmax            = 1000
# plot21.ymin            = 0.01
# plot21.ymax            = 100
# #plot21.lpos = "bottom-center"
# plot21.name            = "Mlj_allOtherCuts"
# plot21.addZUncBand     = zUncBand
# plot21.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_selection)


#-----------------------------------------------------------------------------------


# List of plots to be plotted
plots = [plot0, plot0_ylog, plot1, plot1b, plot1c, plot2_nojet, plot3_nojet, plot4_nojet, plot5_nojet,
         plot6, plot2, plot3, plot4, plot5, 
         plot2and4, plot3and5, plot7, plot8, plot9, plot10, plot7and9, plot8and10,
         plot11, plot11_ele, plot11_jet,
         plot12, plot13, plot13_ylog, plot14, plot14_ylog,
         plot15, plot16]  # produced using preselection root file
#         plot20, plot21] # produced using full selection root file



############# USER CODE - END ################################################
##############################################################################


#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()
c.Print("allPlots.ps[")
for plot in plots:
    plot.Draw("allPlots.ps")
c.Print("allPlots.ps]")
