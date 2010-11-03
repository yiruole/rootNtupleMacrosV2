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
    new = histo.Clone()
    #ROOT.SetOwnership( new, True )    
    if(scale!=1):
        new.Scale(scale)
    if( not new):
        print "ERROR: histo " + histoName + " not found in " + file.GetName()
        print "exiting..."
        sys.exit()
    return new

def generateHistoList( histoBaseName , samples , variableName, fileName , scale = 1):
    histolist = []
    for sample in samples:
        hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
        histolist.append(GetHisto(hname, fileName, scale).Clone())
    return histolist
                                                    
def generateAndAddHistoList( histoBaseName , samples , variableNames, fileName , scale = 1):
    histolist = []
    for sample in samples:
        iv=0
        for variableName in variableNames:
            hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
            if (iv==0):
                histo = GetHisto(hname, fileName, scale).Clone()
            else:
                histo.Add(GetHisto(hname, fileName, scale).Clone())
            iv=iv+1
        histolist.append(histo)
    return histolist

def generateHisto( histoBaseName , sample , variableName, fileName , scale = 1):
    hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
    histo = GetHisto(hname, fileName).Clone()
    return histo

def generateAndAddHisto( histoBaseName , sample , variableNames, fileName , scale = 1):
    iv=0
    for variableName in variableNames:
        hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
        if (iv==0):
            histo = GetHisto(hname, fileName, scale).Clone()
        else:
            histo.Add(GetHisto(hname, fileName, scale).Clone())
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
    lint        = "21.9 pb^{-1}" # integrated luminosity of the sample ( example "10 pb^{-1}" )
    addZUncBand = "no" # add an uncertainty band coming from the data-MC Z+jets rescaling (default = "no", option="yes")
    ZUncKey     = "Z/#gamma/Z* + jets unc." # key to be put in the legend for the Z+jets uncertainty band
    ZPlotIndex  = 1 # index of the Z+jets plots in the histosStack list (default = 1)
    ZScaleUnc   = 0.20 # uncertainty of the data-MC Z+jets scale factor
    makeRatio   = "" # 1=simple ratio, 2=ratio of cumulative histograms
    xbins       = "" #array with variable bin structure
    histodata   = "" # data histogram

    def Draw(self, fileps):

        #-- create canvas
        canvas = TCanvas()
        stack = {}

        if(plot.makeRatio==1):    
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
        #             xlog may npot work         if (plot.xlog     == "yes"):
        #             fPads1.SetLogx()
        if (plot.ylog     == "yes"):
            fPads1.SetLogy()

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
                    stack[iter].GetXaxis().SetRangeUser(plot.xmin,plot.xmax-0.000001)
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
            legend.AddEntry(plot.histodata, "data","p")
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

        legend.Draw()
        canvas.Update()
        gPad.RedrawAxis()
        gPad.Modified()

        #-- 2nd pad (ratio)
        if(plot.makeRatio==1):    
            fPads2.cd()
            #fPads2.SetLogy()
            h_bkgTot = stack[0].Clone() 
            h_ratio = plot.histodata.Clone()

            if( plot.xbins!="" ): ## Variable binning
                xbinsFinal = array( 'd', plot.xbins )
                lenght = len(xbinsFinal)-1            
                h_bkgTot1 = h_bkgTot.Rebin( lenght , "h_bkgTot1", xbinsFinal) 
                h_ratio1 = h_ratio.Rebin( lenght , "h_ratio1" , xbinsFinal)
            else:                 ## Fixed binning
                h_bkgTot1 = h_bkgTot.Rebin( 1 , "h_bkgTot1" ) 
                h_ratio1 = h_ratio.Rebin( 1 , "h_ratio1" )            

            h_ratio1.SetStats(0)
            if (plot.xmin!="" and plot.xmax!=""):
                h_bkgTot1.GetXaxis().SetRangeUser(plot.xmin,plot.xmax-0.000001)
                h_ratio1.GetXaxis().SetRangeUser(plot.xmin,plot.xmax-0.000001)
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
            
            if (plot.xmin!="" and plot.xmax!=""):
                lineAtOne = TLine(plot.xmin,1,plot.xmax,1)
            else:
                lineAtOne = TLine(h_ratio.GetXaxis().GetXmin(),1,h_ratio.GetXaxis().GetXmax(),1)
            lineAtOne.SetLineColor(2)
            lineAtOne.Draw()


        #-- end
        canvas.SaveAs(plot.name + ".eps","eps")
        #canvas.SaveAs(plot.name + ".pdf","pdf") # do not use this line because root creates rotated pdf plot - see end of the file instead
        canvas.Print(fileps)



############# DON'T NEED TO MODIFY ANYTHING HERE - END #######################
##############################################################################


##############################################################################
############# USER CODE - BEGING #############################################

#--- Input root file

File_preselection = GetFile("$LQDATA/enujj_analysis/21.9pb-1_v2/output_cutTable_enujjSample/analysisClass_enujjSample_plots.root")

File_selection    = File_preselection

UseQCDFromData    = 1 # always put an existing file under File_QCD (otherwise the code will crash)

File_QCD          = GetFile("$LQDATA/enujj_analysis/7.4pb-1_v5_QCD_HLT30_AllDeltaPhi/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")
QCDscaleFactor    = 2.96 # ratio between integrated lumi of the signal
                        # sample (i.e. 21.9 pb-1) / integrated lumi of the QCD sample (i.e. 7.4 pb-1 from HLT Photon20)

##File_QCD          = GetFile("$LQDATA/enujj_analysis/2.5pb-1_v3_usePF_QCD_HLT20/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_plots.root")
##QCDscaleFactor    = 5.964 # ratio between integrated lumi of the signal
                          # sample (i.e. 15.1 pb-1) / integrated lumi of the QCD sample (i.e. 2.532 pb-1 from HLT Photon20)

#### Common values for plots:
#otherBkgsKey="QCD, single top, VV+jets, W/W*+jets"
otherBkgsKey="Other Bkgs"
zUncBand="no"
makeRatio=1

pt_xmin=0
pt_xmax=800
pt_ymin=0.001
pt_ymax=1000

eta_rebin=4
eta_ymin=0
eta_ymax=220

#--- Final plots are defined here

# Simply add or remove plots and update the list plots = [plot0, plot1, ...] below

histoBaseName = "histo1D__SAMPLE__cutHisto_allPreviousCuts________VARIABLE"
histoBaseName_userDef = "histo1D__SAMPLE__VARIABLE"

samplesForStackHistosQCD = ["DATA"]
samplesForStackHistos = ["TTbar_Madgraph","WJetAlpgen","OTHERBKG"]
keysStack =             ["QCD","ttbar", "W/W* + jets", otherBkgsKey]

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
plot1.xtit            = "pT 1st electron (GeV/c)"
plot1.ytit            = "Number of events"
plot1.ylog            = "yes"
plot1.rebin           = 1
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

#--- MET_PAS  ---
variableName = "MET_PAS"

plot3 = Plot()
## inputs for stacked histograms
plot3.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot3.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot3.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot3.keys            = keys
plot3.xtit            = "pfMET (GeV/c)"
plot3.ytit            = "Number of events"
plot3.ylog            = "yes"
plot3.rebin           = 2
plot3.xmin            = 0
plot3.xmax            = 500
plot3.ymin            = 0.001
plot3.ymax            = 1000
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
plot4.xtit            = "min(pT 1st electron,pfMET) (GeV/c)"
plot4.ytit            = "Number of events"
plot4.ylog            = "yes"
plot4.rebin           = 1
plot4.xmin            = 0
plot4.xmax            = 300
plot4.ymin            = 0.001
plot4.ymax            = 1000
#plot4.lpos = "bottom-center"
plot4.name            = "minMETPt1stEle_allPreviousCuts"
plot4.addZUncBand     = zUncBand
plot4.makeRatio       = makeRatio
plot4.xbins           = [0,5,10,20,25,30,35,40,45,50,55,60,65,70,75,80,90,100,125,150,200,300]
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
plot6.xtit            = "pT 1st jet (GeV/c)"
plot6.ytit            = "Number of events"
plot6.ylog            = "yes"
plot6.rebin           = 1
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
plot8.xtit            = "pT 2nd jet (GeV/c)"
plot8.ytit            = "Number of events"
plot8.ylog            = "yes"
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
plot6and8.xtit            = "pT jets (GeV/c)"
plot6and8.ytit            = "Number of events"
plot6and8.ylog            = "yes"
plot6and8.rebin           = 1
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
plot11.rebin           = 6
#plot11.xmin            = 0
#plot11.xmax            = 3.14
plot11.ymin            = 0.001
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
plot12.rebin           = 6
#plot12.xmin            = 0
#plot12.xmax            = 3.146
plot12.ymin            = 0.001
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
plot13.rebin           = 6
#plot13.xmin            = 0
#plot13.xmax            = 3.146
plot13.ymin            = 0.001
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
plot13_afterOtherDfCuts.rebin           = 6
#plot13_afterOtherDfCuts.xmin            = 0
#plot13_afterOtherDfCuts.xmax            = 3.146
plot13_afterOtherDfCuts.ymin            = 0.001
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
plot14.xtit            = "M_{T}(e\\nu) (GeV/c^{2})"
plot14.ytit            = "Number of events"
# plot14.ylog            = "yes"
# plot14.rebin           = 1
# plot14.ymin            = 0.00000001
# plot14.ymax            = 20
plot14.ylog            = "no"
plot14.rebin           = 1
plot14.xmin            = 0
plot14.xmax            = 400
plot14.ymin            = 0
plot14.ymax            = 250
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
plot14_ylog.xtit            = "M_{T}(e\\nu) (GeV/c^{2})"
plot14_ylog.ytit            = "Number of events"
plot14_ylog.ylog            = "yes"
plot14_ylog.rebin           = 1 # don't change it (since a rebinning is already applied above on the same histo)
plot14_ylog.xmin            = 0
plot14_ylog.xmax            = 500
plot14_ylog.ymin            = 0.001
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
plot14_0_1.xtit            = "M_{T}(e\\nu) (GeV/c^{2}) - #Delta#phi(MET,2nd jet)<=1"
plot14_0_1.ytit            = "Number of events"
# plot14_0_1.ylog            = "yes"
# plot14_0_1.rebin           = 1
# plot14_0_1.ymin            = 0.00000001
# plot14_0_1.ymax            = 20
plot14_0_1.ylog            = "no"
plot14_0_1.rebin           = 1
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
plot14_1_2.xtit            = "M_{T}(e\\nu) (GeV/c^{2}) - 1<#Delta#phi(MET,2nd jet)<2"
plot14_1_2.ytit            = "Number of events"
# plot14_1_2.ylog            = "yes"
# plot14_1_2.rebin           = 1
# plot14_1_2.ymin            = 0.00000001
# plot14_1_2.ymax            = 20
plot14_1_2.ylog            = "no"
plot14_1_2.rebin           = 1
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
plot14_2_pi.xtit            = "M_{T}(e\\nu) (GeV/c^{2}) - #Delta#phi(MET,2nd jet)>=2"
plot14_2_pi.ytit            = "Number of events"
# plot14_2_pi.ylog            = "yes"
# plot14_2_pi.rebin           = 1
# plot14_2_pi.ymin            = 0.00000001
# plot14_2_pi.ymax            = 20
plot14_2_pi.ylog            = "no"
plot14_2_pi.rebin           = 1
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



#--- sT_PAS ---
variableName = "sT_PAS"

plot15 = Plot()
## inputs for stacked histograms
plot15.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot15.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot15.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot15.keys            = keys
plot15.xtit            = "St (GeV/c)"
plot15.ytit            = "Number of events"
plot15.ylog            = "yes"
plot15.rebin           = 2
plot15.xmin            = 50
plot15.xmax            = 2000
plot15.ymin            = 0.001
plot15.ymax            = 1000
#plot15.lpos = "bottom-center"
plot15.name            = "sT_allPreviousCuts"
plot15.addZUncBand     = zUncBand
plot15.makeRatio       = makeRatio
plot15.xbins           = [50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370,400,500,600,700,800,1000,1500,2000]
plot15.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Mjj_PAS (after preselection) ---
variableName = "Mjj_PAS"

plot16 = Plot()
## inputs for stacked histograms
plot16.histosStack     = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDscaleFactor) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot16.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot16.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot16.keys            = keys
plot16.xtit            = "M(jj) (GeV/c^{2})"
plot16.ytit            = "Number of events"
plot16.ylog            = "yes"
plot16.rebin           = 2
plot16.ymin            = 0.001
plot16.ymax            = 1000
plot16.xmin            = 0
plot16.xmax            = 2000
#plot16.lpos = "bottom-center"
plot16.name            = "Mjj_FullPreSel_allPreviousCuts_ylin"
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
plot17.xtit            = "M(ej) (GeV/c^{2})"
plot17.ytit            = "Number of events x 2"
plot17.ylog            = "no"
plot17.rebin           = 1
plot17.xmin            = 0
plot17.xmax            = 1500
plot17.ymin            = 0.
plot17.ymax            = 300 #450
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
plot17_ylog.xtit            = "M(ej) (GeV/c^{2})"
plot17_ylog.ytit            = "Number of events x 2"
plot17_ylog.ylog            = "yes"
plot17_ylog.rebin           = 2
plot17_ylog.xmin            = 0
plot17_ylog.xmax            = 2000
plot17_ylog.ymin            = 0.001
plot17_ylog.ymax            = 1000
#plot17_ylog.lpos = "bottom-center"
plot17_ylog.name            = "Mej_allPreviousCuts"
plot17_ylog.addZUncBand     = zUncBand
plot17_ylog.makeRatio       = makeRatio
plot17_ylog.xbins           = [0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,325,350,375,400,425,450,475,500,550,600,800,1000,1500,2000]
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
plot18.xtit            = "M_{T}(\\nuj) (GeV/c^{2})"
plot18.ytit            = "Number of events x 2"
plot18.ylog            = "yes"
plot18.rebin           = 2
plot18.xmin            = 0
plot18.xmax            = 1000
plot18.ymin            = 0.001
plot18.ymax            = 1000
#plot18.lpos = "bottom-center"
plot18.name            = "MTnuj_allPreviousCuts"
plot18.addZUncBand     = zUncBand
plot18.makeRatio       = makeRatio
plot18.xbins           = [0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,425,450,475,500,550,600,800,1000]
plot18.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)


#--- Phi1stEle_PAS  ---
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
plot19.rebin           = eta_rebin
#plot19.xmin            = -3.15
#plot19.xmax            = 3.15
plot19.ymin            = eta_ymin
plot19.ymax            = eta_ymax
#plot19.lpos            = "top-left"
plot19.name            = "Phi1stEle_allPreviousCuts"
plot19.addZUncBand     = zUncBand
plot19.makeRatio       = makeRatio
plot19.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


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
plot20.rebin           = eta_rebin
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
plot21.rebin           = eta_rebin
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
plot22.rebin           = eta_rebin
plot22.ymin            = eta_ymin
plot22.ymax            = eta_ymax
#plot22.lpos = "bottom-center"
plot22.name            = "METPhi_allPreviousCuts"
plot22.addZUncBand     = zUncBand
plot22.makeRatio       = makeRatio
plot22.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)



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
plot30.xtit            = "St (GeV/c)"
plot30.ytit            = "Number of events"
plot30.ylog            = "yes"
plot30.rebin           = 10
plot30.xmin            = 100
plot30.xmax            = 2000
plot30.ymin            = 0.001
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
plot31.xtit            = "M(ej) (GeV/c^{2})"
plot31.ytit            = "Number of events x 2"
plot31.ylog            = "yes"
plot31.rebin           = 10
plot31.xmin            = 0
plot31.xmax            = 2000
plot31.ymin            = 0.001
plot31.ymax            = 500
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
plot32.xtit            = "M_{T}(\\nuj) (GeV/c^{2})"
plot32.ytit            = "Number of events x 2"
plot32.ylog            = "yes"
plot32.rebin           = 10
plot32.xmin            = 0
plot32.xmax            = 2000
plot32.ymin            = 0.001
plot32.ymax            = 500
#plot32.lpos = "bottom-center"
plot32.name            = "MTnuj_fullSelection"
plot32.addZUncBand     = zUncBand
plot32.makeRatio       = makeRatio
plot32.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)



#-----------------------------------------------------------------------------------


# List of plots to be plotted
plots  = [plot0, plot1, plot2, plot3, plot4, plot5, plot6, plot7, plot8, plot9
          , plot6and8, plot7and9
          , plot10, plot11, plot12, plot13, plot13_afterOtherDfCuts, plot14, plot14_ylog
          , plot14_0_1, plot14_1_2, plot14_2_pi  
          , plot15, plot16
          , plot17, plot17_ylog, plot18
          , plot19, plot20, plot21, plot22 
          , plot30, plot31, plot32]



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
# print "Converting eps files into pdf files ..."
# for plot in plots:
#     system("convert "+plot.name+".eps "+plot.name+".pdf") # instead, uncomment this line to create pdf plots (a bit slow)

# create a file with list of eps and pdf files, to facilitate copying them to svn area for AN/PAS/PAPER
print "Creating file listOfEpsFiles.txt and listOfPdfFiles.txt ..."
system("rm listOfEpsFiles.txt listOfPdfFiles.txt")
for plot in plots:
    system("echo "+plot.name+".eps >> listOfEpsFiles.txt; echo "+plot.name+".pdf >> listOfPdfFiles.txt") 
print "Use for example: scp `cat listOfPdfFiles.txt` pcuscms41:/home/prumerio/doc/cms/phys/leptoquark/ANPAS2010/papers/EXO-10-005/trunk/plots"

