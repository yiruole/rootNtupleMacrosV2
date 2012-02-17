
##############################################################################
## USER CODE IS TOWARD THE END OF THE FILE
##############################################################################

##############################################################################
############# DON'T NEED TO MODIFY ANYTHING HERE - BEGIN #####################

#---Import
import sys, math
import string
from optparse import OptionParser
import os.path
from ROOT import *
import re
from array import array
import copy
import numpy

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
    new = copy.deepcopy(histo)
    if(scale!=1):
        new.Scale(scale)
    return new

def generateIntegralHisto ( original_histo ) :
    integral_histo = copy.deepcopy ( original_histo )
    
    n_bins = original_histo.GetNbinsX()

    for ibin in range (0, n_bins + 1 ):
        error = Double () 
        integral = original_histo.IntegralAndError ( ibin, n_bins + 1, error ) 
        integral_histo.SetBinContent ( ibin, integral ) 
        integral_histo.SetBinError   ( ibin, error    ) 
    
    return integral_histo


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
    data_file_name = ""
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
    lint        = "2130 pb^{-1}" # integrated luminosity of the sample ( example "10 pb^{-1}" )
    addZUncBand = "no" # add an uncertainty band coming from the data-MC Z+jets rescaling (default = "no", option="yes")
    ZUncKey     = "Z/#gamma/Z* + jets unc." # key to be put in the legend for the Z+jets uncertainty band
    ZPlotIndex  = 1 # index of the Z+jets plots in the histosStack list (default = 1)
    ZScaleUnc   = 0.20 # uncertainty of the data-MC Z+jets scale factor
    makeRatio   = ""  # 1=simple ratio, 2=ratio of cumulative histograms
    makeNSigma  = "" # 1= yes, 0 = no
    xbins       = "" #array with variable bin structure
    histodata   = "" # data histogram
    gif_folder  = "/tmp/"
    eps_folder  = "/tmp/"
    lumi_pb = "0.0"
    suffix = ""
    stackColorIndexes = []
    stackFillStyleIds = []
    is_integral = False

    def Draw(self, fileps):


        self.histos = rebinHistos( self.histos, self.xmin, self.xmax, self.rebin, self.xbins, self.addOvfl )
        self.histosStack = rebinHistos( self.histosStack, self.xmin, self.xmax, self.rebin, self.xbins, self.addOvfl )
        resultArray = rebinHisto( self.histodata, self.xmin, self.xmax, self.rebin, self.xbins, self.addOvfl )
        self.histodata = resultArray[0]
        minBinW = resultArray[1]

        if self.is_integral : 

            integral_histoStack = []
            integral_histos     = []
            integral_histodata  = generateIntegralHisto ( self.histodata ) 

            for histo in self.histosStack:
                integral_histo = generateIntegralHisto ( histo ) 
                integral_histoStack.append ( integral_histo ) 

            for histo in self.histos:
                integral_histo = generateIntegralHisto ( histo ) 
                integral_histos.append ( integral_histo ) 
                
            self.histodata = integral_histodata
            self.histosStack = integral_histoStack
            self.histos = integral_histos
            
        #-- create canvas
        canvas = TCanvas()
        stack = {}

        if(self.makeNSigma==0):
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
        
        if(self.makeNSigma==1):
            if(self.makeRatio==1):
                fPads1 = TPad("pad1", "", 0.00, 0.40, 0.99, 0.99)
                fPads3 = TPad("pad3", "", 0.00, 0.20, 0.99, 0.40)
                fPads2 = TPad("pad2", "", 0.00, 0.00, 0.99, 0.20)
                fPads1.SetFillColor(0)
                fPads1.SetLineColor(0)
                fPads2.SetFillColor(0)
                fPads2.SetLineColor(0)
                fPads3.SetFillColor(0)
                fPads3.SetLineColor(0)
                fPads1.Draw()
                fPads2.Draw()
                fPads3.Draw()
            else:
                fPads1 = TPad("pad1", "", 0.00, 0.20, 0.99, 0.99)
                fPads3 = TPad("pad3", "", 0.00, 0.00, 0.99, 0.20)
                fPads1.SetFillColor(0)
                fPads1.SetLineColor(0)
                fPads3.SetFillColor(0)
                fPads3.SetLineColor(0)
                fPads1.Draw()
                fPads3.Draw()
                

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
        hsize=0.20
        vsize=0.35
        if (self.lpos=="bottom-center"):
            xstart=0.35
            ystart=0.25
        elif(self.lpos=="top-left"):
            xstart=0.12
#            ystart=0.63
            ystart=0.54
        else:
            xstart=0.75
            ystart=0.52
#            xstart=0.65
#            ystart=0.63
        legend = TLegend(xstart, ystart, xstart+hsize, ystart+vsize)
        legend.SetFillColor(kWhite)
        legend.SetBorderSize(0)
        legend.SetShadowColor(10)
        legend.SetMargin(0.2)
        legend.SetTextFont(132)
        legend.AddEntry(self.histodata, "Data, "+self.lumi_pb+" pb^{-1}","lp")

        #-- loop over histograms (stacked)
        Nstacked = len(self.histosStack)
        #stackColorIndexes = [20,38,14,45,20,38,14,45]
        #stackColorIndexes = [20,38,12,14,20,38,12,14]
        
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
            stack[iter].SetMarkerColor(self.stackColorIndexes[iter])
            stack[iter].SetLineColor(  self.stackColorIndexes[iter])
            stack[iter].SetLineWidth(  2 )
            stack[iter].SetFillColor(  self.stackColorIndexes[iter])
            stack[iter].SetFillStyle(  self.stackFillStyleIds[iter])
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
        dataLineIndexes = [1,2,3,1,2,3]
        # dataLineIndexes = [2,1,3,1,2,3]
        
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
            self.histodata.SetMarkerStyle(1)
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
            l.DrawLatex(xstart-hsize+0,ystart+vsize-0.05,"CMS")
            l.DrawLatex(xstart-hsize+0,ystart+vsize-0.15,"#sqrt{s} = 7 TeV")
#            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.03,"CMS Preliminary 2010")
#            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.13,"#intLdt = " + self.lint)
        
        legend.Draw()
        canvas.Update()
        gPad.RedrawAxis()
        gPad.Modified()

        #-- 2nd pad (ratio)
        if(self.makeRatio==1 or self.makeNSigma==1):

            h_bkgTot = copy.deepcopy(stack[0])
            h_ratio = copy.deepcopy(self.histodata)
            h_nsigma = copy.deepcopy(self.histodata)
            h_bkgTot1 = TH1F()
            h_ratio1 = TH1F()
            h_nsigma1 = TH1F()

            if( self.xbins!="" and self.rebin!="var" ): ## Variable binning
                xbinsFinal = array( 'd', self.xbins )
                length = len(xbinsFinal)-1
                h_bkgTot1 = h_bkgTot.Rebin( length , "h_bkgTot1", xbinsFinal)
                h_ratio1 = h_ratio.Rebin( length , "h_ratio1" , xbinsFinal)
                h_nsigma1 = h_nsigma.Rebin( length, "h_nsigma1", xbinsFinal)
            else:
                h_bkgTot1 = h_bkgTot
                h_ratio1 = h_ratio
                h_nsigma1 = h_nsigma

            h_ratio1.SetStats(0)
            if( self.xmin!="" and self.xmax!="" and self.rebin!="var" ):
                h_bkgTot1.GetXaxis().SetRangeUser(self.xmin,self.xmax-0.000001)
                h_ratio1.GetXaxis().SetRangeUser(self.xmin,self.xmax-0.000001)
                h_nsigma1.GetXaxis().SetRangeUser(self.xmin,self.xmax-0.000001)

            if ( self.makeRatio == 1 ):
                fPads2.cd()
                # fPads2.SetLogy()
                h_ratio1.Divide(h_bkgTot1)
                
                h_ratio1.GetXaxis().SetTitle("")
                h_ratio1.GetXaxis().SetTitleSize(0.06)
                h_ratio1.GetXaxis().SetLabelSize(0.1)
                h_ratio1.GetYaxis().SetRangeUser(0.,2)
                h_ratio1.GetYaxis().SetTitle("Data/MC")
                h_ratio1.GetYaxis().SetLabelSize(0.1)
                h_ratio1.GetYaxis().SetTitleSize(0.13)
                h_ratio1.GetYaxis().SetTitleOffset(0.3)
                h_ratio1.SetMarkerStyle ( 1 )
                
                h_ratio1.Draw("p")

                lineAtOne = TLine(h_ratio.GetXaxis().GetXmin(),1,h_ratio.GetXaxis().GetXmax(),1)
                lineAtOne.SetLineColor(2)
                lineAtOne.Draw()

            if ( self.makeNSigma == 1 ):

                fPads3.cd()

                nsigma_x = []
                nsigma_y = []

                for ibin in range (0, h_nsigma1.GetNbinsX() + 1):

                    data  = h_nsigma1.GetBinContent (ibin)
                    bkgd  = h_bkgTot1.GetBinContent (ibin)
                    eData = h_nsigma1.GetBinError   (ibin)
                    eBkgd = h_bkgTot1.GetBinError   (ibin)

                    x = h_ratio1.GetBinCenter ( ibin )

                    diff   = data - bkgd
                    sigma  = math.sqrt ( ( eData * eData ) + ( eBkgd * eBkgd ))

                    if ( sigma != 0.0 and data != 0.0 ):
                        nsigma_x.append ( float (x ) )
                        nsigma_y.append ( float (diff / sigma ) )
                
                    if len ( nsigma_x ) != 0 : 
               		
                        nsigma1_x = numpy.array ( nsigma_x ) 
                        nsigma1_y = numpy.array ( nsigma_y ) 

                        xmin = h_ratio.GetXaxis().GetXmin()
                        xmax = h_ratio.GetXaxis().GetXmax()

                        # if ( len ( nsigma_x ) == 0 ) continue

                        g_nsigma = TGraph ( len ( nsigma_x ) , nsigma1_x, nsigma1_y ) 

                        g_nsigma.Draw("AP")
                        g_nsigma.SetMarkerStyle ( 8 )

                        g_nsigma.SetMarkerSize ( 0.4 )

                        g_nsigma.GetHistogram().GetXaxis().SetTitle("")
                        g_nsigma.GetHistogram().GetXaxis().SetTitleSize(0.06)
                        g_nsigma.GetHistogram().GetXaxis().SetLabelSize(0.1)
                        g_nsigma.GetHistogram().GetYaxis().SetRangeUser(-8.,8)
                        g_nsigma.GetHistogram().GetXaxis().SetLimits ( self.histodata.GetXaxis().GetXmin(), self.histodata.GetXaxis().GetXmax() )

                        g_nsigma.GetHistogram().GetYaxis().SetTitle("N(#sigma) Diff")
                        g_nsigma.GetHistogram().GetYaxis().SetLabelSize(0.1)
                        g_nsigma.GetHistogram().GetYaxis().SetTitleSize(0.13)
                        g_nsigma.GetHistogram().GetYaxis().SetTitleOffset(0.3)

                        lineAtZero     = TLine(h_ratio.GetXaxis().GetXmin(),0 ,h_ratio.GetXaxis().GetXmax(),0 )
                        lineAtPlusTwo  = TLine(h_ratio.GetXaxis().GetXmin(),2 ,h_ratio.GetXaxis().GetXmax(),2 )
                        lineAtMinusTwo = TLine(h_ratio.GetXaxis().GetXmin(),-2,h_ratio.GetXaxis().GetXmax(),-2)

                        lineAtZero    .SetLineColor(2) 
                        lineAtPlusTwo .SetLineColor(1) 
                        lineAtMinusTwo.SetLineColor(1) 

                        g_nsigma.Draw("AP")
                        lineAtZero.Draw("SAME")
                        lineAtPlusTwo .Draw("SAME")
                        lineAtMinusTwo.Draw("SAME")

        #-- end
        if self.suffix == "" : 
            canvas.SaveAs(self.eps_folder + "/" + self.name + ".eps","eps")
            canvas.SaveAs(self.gif_folder + "/" + self.name + ".gif","gif")
        else:
            canvas.SaveAs(self.eps_folder + "/" + self.name + "_" + self.suffix +  ".eps","eps")
            canvas.SaveAs(self.gif_folder + "/" + self.name + "_" + self.suffix +  ".gif","gif")
        #canvas.SaveAs(self.name + ".png","png")
        #canvas.SaveAs(self.name + ".root","root")
        #canvas.SaveAs(self.name + ".pdf","pdf") # do not use this line because root creates rotated pdf plot - see end of the file instead
        canvas.Print(fileps)


