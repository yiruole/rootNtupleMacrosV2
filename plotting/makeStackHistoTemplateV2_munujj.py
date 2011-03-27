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
    lint        = "36 pb^{-1}" # integrated luminosity of the sample ( example "10 pb^{-1}" )
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
        hsize=0.33
        vsize=0.33
        if (self.lpos=="bottom-center"):
            xstart=0.35
            ystart=0.25
        elif(self.lpos=="top-left"):
            xstart=0.12
#            ystart=0.63
            ystart=0.54
        else:
            xstart=0.58
            ystart=0.54
#            xstart=0.65
#            ystart=0.63
        legend = TLegend(xstart, ystart, xstart+hsize, ystart+vsize)
        legend.SetFillColor(kWhite)
        legend.SetBorderSize(0)
        legend.SetShadowColor(10)
        legend.SetMargin(0.2)
        legend.SetTextFont(132)
        legend.AddEntry(plot.histodata, "Data, 36 pb^{-1}","lp")

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
                stack[iter].GetXaxis().SetTitleOffset(0.8)
                stack[iter].GetXaxis().SetLabelOffset(0.0)
                stack[iter].GetXaxis().SetTitleSize(0.055)
                stack[iter].GetXaxis().SetLabelSize(0.045)
                stack[iter].GetXaxis().SetLabelFont(132)
                stack[iter].GetYaxis().SetTitleFont(132)
                stack[iter].GetYaxis().SetTitleOffset(0.75)
                stack[iter].GetYaxis().SetTitleSize(0.055)
                stack[iter].GetYaxis().SetLabelSize(0.045)
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
            legend.AddEntry(plot.histodata, "","")
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
            l.DrawLatex(xstart-hsize+0.17,ystart+vsize-0.05,"CMS")
            l.DrawLatex(xstart-hsize+0.15,ystart+vsize-0.15,"#sqrt{s}=7 TeV")
#            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.03,"CMS Preliminary 2010")
#            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.13,"#intLdt = " + self.lint)

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
        canvas.SaveAs(self.name + ".png","png")
        #canvas.SaveAs(self.name + ".root","root")
        #canvas.SaveAs(self.name + ".pdf","pdf") # do not use this line because root creates rotated pdf plot - see end of the file instead
        canvas.Print(fileps)


############# DON'T NEED TO MODIFY ANYTHING HERE - END #######################
##############################################################################


##############################################################################
############# USER CODE - BEGING #############################################

#--- Input root file

File_preselection = GetFile("presel_signal300.root")

File_fullselection    = GetFile("fullselection_300.root")

WscaleFactor    = 1.20 # no need to rescale anymore since we are using the HLT prescales (36/35.84 can be ignored)

#### Common values for plots:
otherBkgsKey="Other backgrounds"
zUncBand="no"
makeRatio=0
doExtraPlots = False

#--- Final plots are defined here

# Simply add or remove plots and update the list plots = [plot0, plot1, ...] below

histoBaseName = "VARIABLE_SAMPLE"
#histoBaseName = "histo1D__SAMPLE__cutHisto_allPreviousCuts________VARIABLE"
#histoBaseName_userDef = "histo1D__SAMPLE__VARIABLE"

samplesForStackHistos = ["others","tt"]
samplesForStackHistosW = ["wj"]
keysStack =             [otherBkgsKey,"t#bar{t}", "W/W* + jets"]

samplesForHistos = ["lq"]
keys             = ["LQ, M=300 GeV"]

sampleForDataHisto = "data"


#--- MTmunu_PAS ---
variableName = "wtm"

plot1 = Plot()
## inputs for stacked histograms
plot1.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistosW, variableName, File_preselection, WscaleFactor)
plot1.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot1.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot1.keys            = keys
plot1.xtit            = "m_{T,#mu#nu} [GeV]"
plot1.ytit            = "Number of events"
plot1.ylog            = "yes"
plot1.rebin           = "var" # don't change it (since a rebinning is already applied above on the same histo)
plot1.xmin            = 0
plot1.xmax            = 500
plot1.ymin            = 0.01
plot1.ymax            = 1000
#plot1.lpos = "bottom-center"
plot1.name            = "MTmunu_allPreviousCuts"
plot1.addZUncBand     = zUncBand
plot1.makeRatio       = makeRatio
plot1.xbins           = [0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,110,120,130,140,150,200,300,400,500]
plot1.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

# ############################ Plots below to be done after full selection ######################


##--- sT (after full selection) ---
variableName = "st"

plot2 = Plot()
## inputs for stacked histograms
plot2.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_fullselection) + generateHistoList( histoBaseName, samplesForStackHistosW, variableName, File_fullselection, WscaleFactor)
plot2.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot2.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_fullselection)
plot2.keys            = keys
plot2.xtit            = "S_{T} [GeV]"
plot2.ytit            = "Number of events"
plot2.ylog            = "yes"
plot2.rebin           = "var"#10
plot2.xmin            = 100
plot2.xmax            = 2000
plot2.ymin            = 0.01
plot2.ymax            = 500
#plot2.lpos = "bottom-center"
#plot2.name            = "sT_fullSelection"
plot2.name            = "sT_fullSelection_M200_munujj"
plot2.addZUncBand     = zUncBand
plot2.makeRatio       = makeRatio
plot2.xbins           = [200,300,400,500,600,700,1200]
plot2.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_fullselection)


#-----------------------------------------------------------------------------------


# List of plots to be plotted
plots  = [plot1, plot2]

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
print "Converting eps files into pdf files ..."
for plot in plots:
    system("convert "+plot.name+".eps "+plot.name+".pdf") # instead, uncomment this line to create pdf plots (a bit slow)

# create a file with list of eps and pdf files, to facilitate copying them to svn area for AN/PAS/PAPER
print "Creating file listOfEpsFiles.txt and listOfPdfFiles.txt ..."
system("rm listOfEpsFiles.txt listOfPdfFiles.txt")
for plot in plots:
    system("echo "+plot.name+".eps >> listOfEpsFiles.txt; echo "+plot.name+".pdf >> listOfPdfFiles.txt")
print "Use for example: scp `cat listOfPdfFiles.txt` pcuscms46:/home/santanas/Documents/CMSNotes/notes/AN-10-361/trunk/plots"

