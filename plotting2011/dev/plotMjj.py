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

def Plot4Histo(h1, h2, h3, h4  , rebin , xmin , xmax , ymin , ymax , autoyaxis , xtitle , ytitle, algolabel , masslabel, first, last, outputfile, outputplot):
    c = TCanvas()
    c.SetGridx()
    c.SetGridy()

    hAll = copy.deepcopy(h1) 
    hNoTop = copy.deepcopy(h2) 
    hTop = copy.deepcopy(h3) 
    hGen = copy.deepcopy(h4) 
    
    hGen.Scale( 1 / hGen.GetEntries() )
    hAll.Scale( 1 /  hAll.GetEntries() )
    hNoTop.Scale( 1 /  hAll.GetEntries() )
    hTop.Scale( 1 /  hAll.GetEntries() )
    
    hGen.Rebin(rebin)
    hAll.Rebin(rebin)
    hNoTop.Rebin(rebin)
    hTop.Rebin(rebin)
    
    hGen.SetLineColor(0)
    hAll.SetLineColor(1)
    hNoTop.SetLineColor(2)
    hTop.SetLineColor(4)
    
    hGen.SetFillColor(1)
    hGen.SetFillStyle(3004)
    hTop.SetFillColor(4)
    
    max = 0
    if ( hGen.GetMaximum() > max ):
        max = hGen.GetMaximum()
    if (  hAll.GetMaximum() > max ):
        max =  hAll.GetMaximum()
    if ( hNoTop.GetMaximum() > max ):
        max = hNoTop.GetMaximum()
    if ( hTop.GetMaximum() > max ):
        max = hTop.GetMaximum()
                    
    hAll.GetXaxis().SetRangeUser(xmin,xmax)
    if(autoyaxis):
        hAll.GetYaxis().SetRangeUser(0,max+0.15*max)
    else:
        hAll.GetYaxis().SetRangeUser(ymin,ymax)
    
    hAll.GetXaxis().SetTitle(xtitle)
    hAll.GetYaxis().SetTitle(ytitle)
    hAll.SetTitle("")
    
    hAll.Draw()
    hNoTop.Draw("sames")
    hTop.Draw("samesHISTE")
    hGen.Draw("samesHISTE")
    
    legend = TLegend(0.114943,0.652542,0.46408,0.891949)
    legend.SetFillColor(kWhite)
    legend.SetMargin(0.3)
    legend.AddEntry( hAll,"All decays","l")
    legend.AddEntry(hNoTop,"Axigluon \\rightarrow qqbar (q=u,d,c,s,b)","l")
    legend.AddEntry(hTop,"Axigluon \\rightarrow ttbar","f")
    legend.AddEntry(hGen,"Generated axigluon mass","f")
    legend.Draw()
    
    l = TLatex()
    l.SetTextAlign(12)
    l.SetTextFont(132)
    l.SetTextSize(0.065)
    l.SetNDC()
    l.DrawLatex(0.15,0.60,algolabel)
    l.DrawLatex(0.15,0.40,masslabel)
    
    c.Update()
    gPad.RedrawAxis()
    gPad.Modified()

    if(first and last):    
        c.Print(outputfile)
    if(first and not last):
        c.Print(outputfile+"(")
    if(last and not first):
        c.Print(outputfile+")")
    if(not first and not last):    
        c.Print(outputfile)

    c.SaveAs(outputplot)

def Plot3Histo(h1, h2, h3, rebin , xmin , xmax , ymin , ymax , autoyaxis , xtitle , ytitle, algolabel , masslabel, first, last, outputfile, outputplot):
    c = TCanvas()
    c.SetGridx()
    c.SetGridy()

    hAll = copy.deepcopy(h1) 
    hNoTop = copy.deepcopy(h2) 
    hTop = copy.deepcopy(h3) 
    
    hAll.Scale( 1 /  hAll.GetEntries() )
    hNoTop.Scale( 1 /  hAll.GetEntries() )
    hTop.Scale( 1 /  hAll.GetEntries() )
    
    hAll.Rebin(rebin)
    hNoTop.Rebin(rebin)
    hTop.Rebin(rebin)
    
    hAll.SetLineColor(1)
    hNoTop.SetLineColor(2)
    hTop.SetLineColor(4)

    hTop.SetFillColor(4)
    
    max = 0
    if (  hAll.GetMaximum() > max ):
        max =  hAll.GetMaximum()
    if ( hNoTop.GetMaximum() > max ):
        max = hNoTop.GetMaximum()
    if ( hTop.GetMaximum() > max ):
        max = hTop.GetMaximum()
                    
    hAll.GetXaxis().SetRangeUser(xmin,xmax)
    if(autoyaxis):
        hAll.GetYaxis().SetRangeUser(0,max+0.15*max)
    else:
        hAll.GetYaxis().SetRangeUser(ymin,ymax)
    
    hAll.GetXaxis().SetTitle(xtitle)
    hAll.GetYaxis().SetTitle(ytitle)
    hAll.SetTitle("")
    
    hAll.Draw("HISTE")
    hNoTop.Draw("samesHISTE")
    hTop.Draw("samesHISTE")
    
    legend = TLegend(0.114943,0.652542,0.46408,0.891949)
    legend.SetFillColor(kWhite)
    legend.SetMargin(0.3)
    legend.AddEntry( hAll,"All decays","l")
    legend.AddEntry(hNoTop,"Axigluon \\rightarrow qqbar (q=u,d,c,s,b)","l")
    legend.AddEntry(hTop,"Axigluon \\rightarrow ttbar","f")
    legend.Draw()
    
    l = TLatex()
    l.SetTextAlign(12)
    l.SetTextFont(132)
    l.SetTextSize(0.065)
    l.SetNDC()
    l.DrawLatex(0.15,0.60,algolabel)
    l.DrawLatex(0.15,0.40,masslabel)
    
    c.Update()
    gPad.RedrawAxis()
    gPad.Modified()

    if(first and last):    
        c.Print(outputfile)
    if(first and not last):
        c.Print(outputfile+"(")
    if(last and not first):
        c.Print(outputfile+")")
    if(not first and not last):    
        c.Print(outputfile)

    c.SaveAs(outputplot)

#--- Files

File500  = GetFile("/afs/cern.ch/user/s/santanas/scratch0/Axigluons/rootNtupleAnalyzerV2/AGW500Gen.root");
File1000 = GetFile("/afs/cern.ch/user/s/santanas/scratch0/Axigluons/rootNtupleAnalyzerV2/AGW1000Gen.root");

outputFile = "massPlots.ps"

#--- Histograms

h_Mjj_genlevel_AGW500 = GetHisto("h1_Mass_AG" , File500)

h_Mjj_genjets_AGW500 = GetHisto("h1_Mjj_genJet_1st2nd" , File500)
h_Mjj_genjets_noTop_AGW500 = GetHisto("h1_Mjj_genJet_1st2nd_noTop" , File500)
h_Mjj_genjets_Top_AGW500 = GetHisto("h1_Mjj_genJet_1st2nd_Top" , File500)
h_Mjj_fatgenjets_AGW500 = GetHisto("h1_Mjj_fatgenJet_1st2nd" , File500)
h_Mjj_fatgenjets_noTop_AGW500 = GetHisto("h1_Mjj_fatgenJet_1st2nd_noTop" , File500)
h_Mjj_fatgenjets_Top_AGW500 = GetHisto("h1_Mjj_fatgenJet_1st2nd_Top" , File500)
h_Mjj_recojets_AGW500 = GetHisto("h1_Mjj_recoJet_1st2nd" , File500)
h_Mjj_recojets_noTop_AGW500 = GetHisto("h1_Mjj_recoJet_1st2nd_noTop" , File500)
h_Mjj_recojets_Top_AGW500 = GetHisto("h1_Mjj_recoJet_1st2nd_Top" , File500)
h_Mjj_fatrecojets_AGW500 = GetHisto("h1_Mjj_fatrecoJet_1st2nd" , File500)
h_Mjj_fatrecojets_noTop_AGW500 = GetHisto("h1_Mjj_fatrecoJet_1st2nd_noTop" , File500)
h_Mjj_fatrecojets_Top_AGW500 = GetHisto("h1_Mjj_fatrecoJet_1st2nd_Top" , File500)

h_Mjj_over_massAG_genjets_AGW500 = GetHisto("h1_Mjj_genJet_1st2nd_over_massAG" , File500)
h_Mjj_over_massAG_genjets_noTop_AGW500 = GetHisto("h1_Mjj_genJet_1st2nd_over_massAG_noTop" , File500)
h_Mjj_over_massAG_genjets_Top_AGW500 = GetHisto("h1_Mjj_genJet_1st2nd_over_massAG_Top" , File500)
h_Mjj_over_massAG_fatgenjets_AGW500 = GetHisto("h1_Mjj_fatgenJet_1st2nd_over_massAG" , File500)
h_Mjj_over_massAG_fatgenjets_noTop_AGW500 = GetHisto("h1_Mjj_fatgenJet_1st2nd_over_massAG_noTop" , File500)
h_Mjj_over_massAG_fatgenjets_Top_AGW500 = GetHisto("h1_Mjj_fatgenJet_1st2nd_over_massAG_Top" , File500)
h_Mjj_over_massAG_recojets_AGW500 = GetHisto("h1_Mjj_recoJet_1st2nd_over_massAG" , File500)
h_Mjj_over_massAG_recojets_noTop_AGW500 = GetHisto("h1_Mjj_recoJet_1st2nd_over_massAG_noTop" , File500)
h_Mjj_over_massAG_recojets_Top_AGW500 = GetHisto("h1_Mjj_recoJet_1st2nd_over_massAG_Top" , File500)
h_Mjj_over_massAG_fatrecojets_AGW500 = GetHisto("h1_Mjj_fatrecoJet_1st2nd_over_massAG" , File500)
h_Mjj_over_massAG_fatrecojets_noTop_AGW500 = GetHisto("h1_Mjj_fatrecoJet_1st2nd_over_massAG_noTop" , File500)
h_Mjj_over_massAG_fatrecojets_Top_AGW500 = GetHisto("h1_Mjj_fatrecoJet_1st2nd_over_massAG_Top" , File500)

h_Mjj_genlevel_AGW1000 = GetHisto("h1_Mass_AG" , File1000)

h_Mjj_genjets_AGW1000 = GetHisto("h1_Mjj_genJet_1st2nd" , File1000)
h_Mjj_genjets_noTop_AGW1000 = GetHisto("h1_Mjj_genJet_1st2nd_noTop" , File1000)
h_Mjj_genjets_Top_AGW1000 = GetHisto("h1_Mjj_genJet_1st2nd_Top" , File1000)
h_Mjj_fatgenjets_AGW1000 = GetHisto("h1_Mjj_fatgenJet_1st2nd" , File1000)
h_Mjj_fatgenjets_noTop_AGW1000 = GetHisto("h1_Mjj_fatgenJet_1st2nd_noTop" , File1000)
h_Mjj_fatgenjets_Top_AGW1000 = GetHisto("h1_Mjj_fatgenJet_1st2nd_Top" , File1000)
h_Mjj_recojets_AGW1000 = GetHisto("h1_Mjj_recoJet_1st2nd" , File1000)
h_Mjj_recojets_noTop_AGW1000 = GetHisto("h1_Mjj_recoJet_1st2nd_noTop" , File1000)
h_Mjj_recojets_Top_AGW1000 = GetHisto("h1_Mjj_recoJet_1st2nd_Top" , File1000)
h_Mjj_fatrecojets_AGW1000 = GetHisto("h1_Mjj_fatrecoJet_1st2nd" , File1000)
h_Mjj_fatrecojets_noTop_AGW1000 = GetHisto("h1_Mjj_fatrecoJet_1st2nd_noTop" , File1000)
h_Mjj_fatrecojets_Top_AGW1000 = GetHisto("h1_Mjj_fatrecoJet_1st2nd_Top" , File1000)

h_Mjj_over_massAG_genjets_AGW1000 = GetHisto("h1_Mjj_genJet_1st2nd_over_massAG" , File1000)
h_Mjj_over_massAG_genjets_noTop_AGW1000 = GetHisto("h1_Mjj_genJet_1st2nd_over_massAG_noTop" , File1000)
h_Mjj_over_massAG_genjets_Top_AGW1000 = GetHisto("h1_Mjj_genJet_1st2nd_over_massAG_Top" , File1000)
h_Mjj_over_massAG_fatgenjets_AGW1000 = GetHisto("h1_Mjj_fatgenJet_1st2nd_over_massAG" , File1000)
h_Mjj_over_massAG_fatgenjets_noTop_AGW1000 = GetHisto("h1_Mjj_fatgenJet_1st2nd_over_massAG_noTop" , File1000)
h_Mjj_over_massAG_fatgenjets_Top_AGW1000 = GetHisto("h1_Mjj_fatgenJet_1st2nd_over_massAG_Top" , File1000)
h_Mjj_over_massAG_recojets_AGW1000 = GetHisto("h1_Mjj_recoJet_1st2nd_over_massAG" , File1000)
h_Mjj_over_massAG_recojets_noTop_AGW1000 = GetHisto("h1_Mjj_recoJet_1st2nd_over_massAG_noTop" , File1000)
h_Mjj_over_massAG_recojets_Top_AGW1000 = GetHisto("h1_Mjj_recoJet_1st2nd_over_massAG_Top" , File1000)
h_Mjj_over_massAG_fatrecojets_AGW1000 = GetHisto("h1_Mjj_fatrecoJet_1st2nd_over_massAG" , File1000)
h_Mjj_over_massAG_fatrecojets_noTop_AGW1000 = GetHisto("h1_Mjj_fatrecoJet_1st2nd_over_massAG_noTop" , File1000)
h_Mjj_over_massAG_fatrecojets_Top_AGW1000 = GetHisto("h1_Mjj_fatrecoJet_1st2nd_over_massAG_Top" , File1000)


#--- Plot

Plot4Histo( h_Mjj_genjets_AGW500 , h_Mjj_genjets_noTop_AGW500 , h_Mjj_genjets_Top_AGW500 , h_Mjj_genlevel_AGW500, 
            10 , 0, 799, 0, 0.3, 0, "M(j1,j2) [GeV]", "a.u.", "GenJets AK5", "M500", 1, 0, outputFile,"mjj_genjets_m500.png")
Plot4Histo( h_Mjj_fatgenjets_AGW500 , h_Mjj_fatgenjets_noTop_AGW500 , h_Mjj_fatgenjets_Top_AGW500 , h_Mjj_genlevel_AGW500, 
            10 , 0, 799, 0, 0.3, 0, "M(j1,j2) [GeV]", "a.u.", "FatGenJets (R=1.1)", "M500", 0, 0, outputFile,"mjj_fatgenjets_m500.png")
Plot4Histo( h_Mjj_recojets_AGW500 , h_Mjj_recojets_noTop_AGW500 , h_Mjj_recojets_Top_AGW500 , h_Mjj_genlevel_AGW500, 
            10 , 0, 799, 0, 0.3, 0, "M(j1,j2) [GeV]", "a.u.", "PFJets AK5", "M500", 0, 0, outputFile,"mjj_recojets_m500.png")
Plot4Histo( h_Mjj_fatrecojets_AGW500 , h_Mjj_fatrecojets_noTop_AGW500 , h_Mjj_fatrecojets_Top_AGW500 , h_Mjj_genlevel_AGW500, 
            10 , 0, 799, 0, 0.3, 0, "M(j1,j2) [GeV]", "a.u.", "FatPFJets (R=1.1)", "M500", 0, 0, outputFile,"mjj_fatrecojets_m500.png")

Plot4Histo( h_Mjj_genjets_AGW1000 , h_Mjj_genjets_noTop_AGW1000 , h_Mjj_genjets_Top_AGW1000 , h_Mjj_genlevel_AGW1000, 
            10 , 0, 1399, 0, 0.3, 0, "M(j1,j2) [GeV]", "a.u.", "GenJets AK5", "M1000", 0, 0, outputFile,"mjj_genjets_m1000.png")
Plot4Histo( h_Mjj_fatgenjets_AGW1000 , h_Mjj_fatgenjets_noTop_AGW1000 , h_Mjj_fatgenjets_Top_AGW1000 , h_Mjj_genlevel_AGW1000, 
            10 , 0, 1399, 0, 0.3, 0, "M(j1,j2) [GeV]", "a.u.", "FatGenJets (R=1.1)", "M1000", 0, 0, outputFile, "mjj_fatgenjets_m1000.png")
Plot4Histo( h_Mjj_recojets_AGW1000 , h_Mjj_recojets_noTop_AGW1000 , h_Mjj_recojets_Top_AGW1000 , h_Mjj_genlevel_AGW1000, 
            10 , 0, 1399, 0, 0.3, 0, "M(j1,j2) [GeV]", "a.u.", "PFJets AK5", "M1000",0, 0, outputFile, "mjj_recojets_m1000.png")
Plot4Histo( h_Mjj_fatrecojets_AGW1000 , h_Mjj_fatrecojets_noTop_AGW1000 , h_Mjj_fatrecojets_Top_AGW1000 , h_Mjj_genlevel_AGW1000, 
            10 , 0, 1399, 0, 0.3, 0, "M(j1,j2) [GeV]", "a.u.", "FatPFJets (R=1.1)", "M1000", 0, 0, outputFile, "mjj_fatrecojets_m1000.png")

Plot3Histo( h_Mjj_over_massAG_genjets_AGW500 , h_Mjj_over_massAG_genjets_noTop_AGW500 , h_Mjj_over_massAG_genjets_Top_AGW500 ,
            1 , 0, 1.8, 0, 0.5, 0, "M(j1,j2) / M_{Axigluon}", "a.u.", "GenJets AK5", "M500",0, 0, outputFile,"mjjRatio_genjets_m500.png")
Plot3Histo( h_Mjj_over_massAG_fatgenjets_AGW500 , h_Mjj_over_massAG_fatgenjets_noTop_AGW500 , h_Mjj_over_massAG_fatgenjets_Top_AGW500 ,
            1 , 0, 1.8, 0, 0.5, 0, "M(j1,j2) / M_{Axigluon}", "a.u.", "FatGenJets (R=1.1)", "M500",0, 0, outputFile,"mjjRatio_fatgenjets_m500.png")
Plot3Histo( h_Mjj_over_massAG_recojets_AGW500 , h_Mjj_over_massAG_recojets_noTop_AGW500 , h_Mjj_over_massAG_recojets_Top_AGW500 ,  
            1 , 0, 1.8, 0, 0.5, 0, "M(j1,j2) / M_{Axigluon}", "a.u.", "PFJets AK5", "M500",0, 0, outputFile,"mjjRatio_recojets_m500.png")
Plot3Histo( h_Mjj_over_massAG_fatrecojets_AGW500 , h_Mjj_over_massAG_fatrecojets_noTop_AGW500 , h_Mjj_over_massAG_fatrecojets_Top_AGW500 , 
            1 , 0, 1.8, 0, 0.5, 0, "M(j1,j2) / M_{Axigluon}", "a.u.", "FatPFJets (R=1.1)", "M500",0, 0, outputFile,"mjjRatio_fatrecojets_m500.png")

Plot3Histo( h_Mjj_over_massAG_genjets_AGW1000 , h_Mjj_over_massAG_genjets_noTop_AGW1000 , h_Mjj_over_massAG_genjets_Top_AGW1000 ,
            1 , 0, 1.8, 0, 0.5, 0, "M(j1,j2) / M_{Axigluon}", "a.u.", "GenJets AK5", "M1000",0, 0, outputFile,"mjjRatio_genjets_m1000.png")
Plot3Histo( h_Mjj_over_massAG_fatgenjets_AGW1000 , h_Mjj_over_massAG_fatgenjets_noTop_AGW1000 , h_Mjj_over_massAG_fatgenjets_Top_AGW1000 , 
            1 , 0, 1.8, 0, 0.5, 0, "M(j1,j2) / M_{Axigluon}", "a.u.", "FatGenJets (R=1.1)", "M1000",0, 0, outputFile,"mjjRatio_fatgenjets_m1000.png")
Plot3Histo( h_Mjj_over_massAG_recojets_AGW1000 , h_Mjj_over_massAG_recojets_noTop_AGW1000 , h_Mjj_over_massAG_recojets_Top_AGW1000 , 
            1 , 0, 1.8, 0, 0.5, 0, "M(j1,j2) / M_{Axigluon}", "a.u.", "PFJets AK5", "M1000",0, 0, outputFile,"mjjRatio_recojets_m1000.png")
Plot3Histo( h_Mjj_over_massAG_fatrecojets_AGW1000 , h_Mjj_over_massAG_fatrecojets_noTop_AGW1000 , h_Mjj_over_massAG_fatrecojets_Top_AGW1000 , 
            1 , 0, 1.8, 0, 0.5, 0, "M(j1,j2) / M_{Axigluon}", "a.u.", "FatPFJets (R=1.1)", "M1000",0, 1, outputFile,"mjjRatio_fatrecojets_m1000.png")


## Terminate the program
print "Press ENTER to terminate"
wait=raw_input()
