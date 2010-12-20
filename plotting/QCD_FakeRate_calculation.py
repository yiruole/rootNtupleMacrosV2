#!/usr/bin/env python

## Instructions and description:
##     This code uses defined methods to divide histograms and compute the proper error bars (only for the fixed bin case, not yet for the var bin case)
##     for efficiencies from weighted samples.
##     Change the name of the input files to be the output of analysisClass_QCD_fakeRate.C
##     Towards the bottom of the macro change the names of the histograms passed to GetEffFixBinning to be those for the appropriate dataset.
##     Adjust plotting option at the very end of the macro.


#---Import
import sys
import string
from optparse import OptionParser
import os.path
from ROOT import *
import re
import ROOT
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

def GetHisto( histoName , file, scale = 1 ):
    histo = file.Get( histoName )
    if( not histo):
        print "ERROR: histo " + histoName + " not found in " + file.GetName()
        print "exiting..."
        sys.exit()
    new = copy.deepcopy(histo)
    if(scale!=1):
        new.Scale(scale)
    return new

def GetIntegralTH1( histo, xmin, xmax):
    #get integral
    axis = histo.GetXaxis()
    bmin = axis.FindBin(xmin)
    bmax = axis.FindBin(xmax)
    bminResidual = histo.GetBinContent(bmin)*(xmin-axis.GetBinLowEdge(bmin)) / axis.GetBinWidth(bmin)
    bmaxResidual = histo.GetBinContent(bmax)*(axis.GetBinUpEdge(bmax)-xmax) / axis.GetBinWidth(bmax)
    integral = histo.Integral(bmin,bmax) - bminResidual - bmaxResidual
    return integral

def GetErrorIntegralTH1( histo, xmin, xmax):
    print "## calculating error for integral of histo " + str(histo)
    print "## in the x range [" + str(xmin) + "," + str(xmax) + "]"
    #get integral
    axis = histo.GetXaxis()
    bmin = axis.FindBin(xmin)
    bmax = axis.FindBin(xmax)
    bminResidual = histo.GetBinContent(bmin)*(xmin-axis.GetBinLowEdge(bmin)) / axis.GetBinWidth(bmin)
    bmaxResidual = histo.GetBinContent(bmax)*(axis.GetBinUpEdge(bmax)-xmax) / axis.GetBinWidth(bmax)
    integral = histo.Integral(bmin,bmax) - bminResidual - bmaxResidual
    error = 0
    for bin in range(bmin, bmax+1):
	print "bin: " +str(bin)
        if(bin==bmax and bmaxResidual==histo.GetBinContent(bmax)): # skip last bin if out of range
            print "     --> skip bin: " + str(bin)
        else:
            error = error + histo.GetBinError(bin)**2
            print "error**2 : " + str(error)

    error = sqrt(error)
    print  " "
    return error

def GetEffFixBinning( num , den , m_size , m_style , m_color , xtitle , ytitle , min , max, rebin):
    #print num;
    num.Rebin(rebin)
    den.Rebin(rebin)
    nBins_num = num.GetXaxis().GetNbins()
    nBins_den = den.GetXaxis().GetNbins()
    Xmin = num.GetXaxis().GetXmin()
    Xmax = num.GetXaxis().GetXmax()

##     for i in range (1,nBins_num):
##         print num.GetBinError(i)

    Tmp = TH1F("","",nBins_num,Xmin,Xmax)
    Tmp.Divide(num,den,1,1)
    
    GraphEff = TGraphAsymmErrors()
    GraphEff.SetMarkerSize( m_size )
    GraphEff.SetMarkerStyle( m_style )
    GraphEff.SetMarkerColor( m_color )
    GraphEff.SetLineColor( m_color )
    GraphEff.GetXaxis().SetTitle(xtitle)
    GraphEff.GetYaxis().SetTitle(ytitle)
    GraphEff.GetXaxis().SetRangeUser(min,max);
    for bin in range (1,nBins_num):
        if (Tmp.GetXaxis().GetBinUpEdge(bin)<=min): continue
        if (Tmp.GetXaxis().GetBinLowEdge(bin)>max): continue
        GraphEff.SetPoint(bin-1,Tmp.GetBinCenter(bin),Tmp.GetBinContent(bin))
        #print Tmp.GetBinCenter(bin)
        eff = Tmp.GetBinContent(bin)
        dNum = num.GetBinError(bin)
        dDenom = den.GetBinError(bin)
        Denom = den.GetBinContent(bin)
        if (Denom!=0):
            sigma = sqrt((dNum*dNum)+(eff*eff*dDenom*dDenom))/Denom
            GraphEff.SetPointEYhigh(bin-1,sigma) # 0 is underflow bins for histo, but first point in TGraph
            GraphEff.SetPointEYlow(bin-1,sigma)

    f1 = TF1("f1","pol0",30,500);
    #GraphEff.Fit("f1","R")

    return GraphEff

def GetEffVarBinning( num , den , m_size , m_style , m_color , xtitle , ytitle , min , max, Bins):

    NbinTot = len(Bins) 
    BinsFinal = array( 'f', Bins ) 

    num_varBin = TH1F("num_varBin",
                      "num_varBin",
                      NbinTot - 1 , BinsFinal );
    nBins_num = num.GetXaxis().GetNbins()
    getxmax_num = num.GetXaxis().GetXmax();
    getxmin_num = num.GetXaxis().GetXmin();
    step_num = (getxmax_num - getxmin_num) / nBins_num;
    for bin in range( 1 , nBins_num+1 ):
        #for entry in range( 0 , int(num.GetBinContent(bin)) ):
            num_varBin.Fill(step_num * bin-1, num.GetBinContent(bin))

    den_varBin = TH1F("den_varBin",
                      "den_varBin",
                      NbinTot - 1 , BinsFinal );
    nBins_den = den.GetXaxis().GetNbins()
    getxmax_den = den.GetXaxis().GetXmax();
    getxmin_den = den.GetXaxis().GetXmin();
    step_den = (getxmax_den - getxmin_den) / nBins_den;
    for bin in range( 1 , nBins_den+1 ):
        #for entry in range( 0 , int(den.GetBinContent(bin)) ):
            den_varBin.Fill(step_den * bin-1, den.GetBinContent(bin) )

    #--- create final graph
    Tmp = TH1F("","",NbinTot - 1 , BinsFinal)
    Tmp.Divide(num_varBin,den_varBin,1,1)
    GraphEff = TGraphAsymmErrors(Tmp)

##     for i in range (1,NbinTot):
##         print num_varBin.GetBinError(i))

    GraphEff.SetMarkerSize( m_size )
    GraphEff.SetMarkerStyle( m_style )
    GraphEff.SetMarkerColor( m_color )
    GraphEff.SetLineColor( m_color )
    GraphEff.GetXaxis().SetTitle(xtitle)
    GraphEff.GetYaxis().SetTitle(ytitle)
    GraphEff.GetXaxis().SetRangeUser(min,max);

##     f1 = TF1("f1","pol1",30,100);
##     GraphEff.Fit("f1","R")

    for bin in range (1,nBins_num):
        eff = Tmp.GetBinContent(bin)
        dNum = num_varBin.GetBinError(bin)
        dDenom = den_varBin.GetBinError(bin)
        Denom = den_varBin.GetBinContent(bin)
        if (Denom!=0):
            sigma = sqrt((dNum*dNum)+(eff*eff*dDenom*dDenom))/Denom
            GraphEff.SetPointEYhigh(bin-1,sigma) # 0 is underflow bins for histo, but first point in TGraph
            GraphEff.SetPointEYlow(bin-1,sigma)

    return GraphEff

def GetRatioEff( num , den , m_size , m_style , m_color , xtitle , ytitle ):

    ratio = TGraphAsymmErrors()
    npoint = 0;

    if(num.GetN() != den.GetN()):
        print "ERROR in GetRatioEff: num.GetN() != den.GetN()"
        print "exiting..."
        sys.exit()
        
    for point in range( 0 , den.GetN() ):
        
        x1, y1, x2, y2 = ROOT.Double(1), ROOT.Double(1), ROOT.Double(1), ROOT.Double(1)

        
        num.GetPoint(point, x1, y1)
        ehx1 = num.GetErrorXhigh(point)  
        ehy1 = num.GetErrorYhigh(point)  
        elx1 = num.GetErrorXlow(point)  
        ely1 = num.GetErrorYlow(point)  
        
        den.GetPoint(point, x2, y2)
        ehx2 = den.GetErrorXhigh(point)  
        ehy2 = den.GetErrorYhigh(point)  
        elx2 = den.GetErrorXlow(point)  
        ely2 = den.GetErrorYlow(point)  
        
        if( y2 == 0 or y1==0 ):
            ratio.SetPoint(npoint,x1, 0)
            ratio.SetPointEXhigh(npoint, ehx1)
            ratio.SetPointEXlow(npoint, elx1)
            ratio.SetPointEYhigh(npoint, 0)
            ratio.SetPointEYlow(npoint, 0)
        else:
            r = y1/y2
            erelhy1 = ehy1 / y1
            erelhy2 = ehy2 / y2
            erelly1 = ely1 / y1
            erelly2 = ely2 / y2
            erelhr = sqrt(erelhy1 * erelhy1 + erelhy2 * erelhy2)
            erellr = sqrt(erelly1 * erelly1 + erelly2 * erelly2)
            ehr=erelhr*r
            elr=erellr*r
            
            ratio.SetPoint(npoint, x1, r)
            ratio.SetPointEXhigh(npoint, ehx1)
            ratio.SetPointEXlow(npoint, elx1)
            ratio.SetPointEYhigh(npoint, ehr)
            ratio.SetPointEYlow(npoint, elr)
        
        npoint = npoint + 1

    ratio.SetMarkerSize( m_size )
    ratio.SetMarkerStyle( m_style )
    ratio.SetMarkerColor( m_color )
    ratio.GetXaxis().SetTitle(xtitle)
    ratio.GetYaxis().SetTitle(ytitle)

    return ratio


# %%%%%%% END %%%%%%% 

#QCD fake rate: 33.0pb-1_QCD_fakeRate_run_gt144114_njet_More1_MET_lt35
File1 = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/33.0pb-1_QCD_fakeRate_run_gt144114_njet_More1_MET_lt35/output_cutTable_QCD_fakeRate/analysisClass_QCD_fakeRate_plots.root")

#QCD fake rate: 33.0pb-1_QCD_fakeRate_run_gt144114_njet_More1
#File1 = GetFile("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/33.0pb-1_QCD_fakeRate_run_gt144114_njet_More1/output_cutTable_QCD_fakeRate/analysisClass_QCD_fakeRate_plots.root")


#--- Define all the histograms
#def GetEffFixBinning( num , den , m_size , m_style , m_color , xtitle , ytitle , min , max, rebin):

h1_SuperClusterPt_DATA_all = GetHisto( "histo1D__DATA__h1_SuperClusterPt" , File1)
h1_SuperClusterPt_DATA_barrel1 = GetHisto( "histo1D__DATA__h1_SuperClusterPt_barrel1" , File1)
h1_SuperClusterPt_DATA_barrel2 = GetHisto( "histo1D__DATA__h1_SuperClusterPt_barrel2" , File1)
h1_SuperClusterPt_DATA_endcap1 = GetHisto( "histo1D__DATA__h1_SuperClusterPt_endcap1" , File1)
h1_SuperClusterPt_DATA_endcap2 = GetHisto( "histo1D__DATA__h1_SuperClusterPt_endcap2" , File1)
#h1_SuperClusterEta_DATA_all = GetHisto( "histo1D__DATA__h1_SuperClusterEta" , File1)
#sums
h1_SuperClusterPt_DATA_barrel  = GetHisto( "histo1D__DATA__h1_SuperClusterPt_barrel1" , File1)
h1_SuperClusterPt_DATA_barrel.Add(h1_SuperClusterPt_DATA_barrel2)
h1_SuperClusterPt_DATA_endcap  = GetHisto( "histo1D__DATA__h1_SuperClusterPt_endcap1" , File1)
h1_SuperClusterPt_DATA_endcap.Add(h1_SuperClusterPt_DATA_endcap2)


h1_ElePt_DATA_all = GetHisto( "histo1D__DATA__h1_ElePt" , File1)
h1_ElePt_DATA_barrel1 = GetHisto( "histo1D__DATA__h1_ElePt_barrel1" , File1)
h1_ElePt_DATA_barrel2 = GetHisto( "histo1D__DATA__h1_ElePt_barrel2" , File1)
h1_ElePt_DATA_endcap1 = GetHisto( "histo1D__DATA__h1_ElePt_endcap1" , File1)
h1_ElePt_DATA_endcap2 = GetHisto( "histo1D__DATA__h1_ElePt_endcap2" , File1)
#h1_EleEta_DATA_all = GetHisto( "histo1D__DATA__h1_EleEta" , File1)
#sums
h1_ElePt_DATA_barrel  = GetHisto( "histo1D__DATA__h1_ElePt_barrel1" , File1)
h1_ElePt_DATA_barrel.Add(h1_ElePt_DATA_barrel2)
h1_ElePt_DATA_endcap  = GetHisto( "histo1D__DATA__h1_ElePt_endcap1" , File1)
h1_ElePt_DATA_endcap.Add(h1_ElePt_DATA_endcap2)

# print "-------------------------------------------------------------------------"
# print "-------------------------------------------------------------------------"

# print "h1_SuperClusterPt_DATA_all.GetEntries(): " + str(h1_SuperClusterPt_DATA_all.GetEntries())
# print "-------------------------------------------------------------------------"
# print "h1_SuperClusterPt_DATA_barrel.GetEntries(): " + str(h1_SuperClusterPt_DATA_barrel.GetEntries())
# print "h1_SuperClusterPt_DATA_barrel1.GetEntries(): " + str(h1_SuperClusterPt_DATA_barrel1.GetEntries())
# print "h1_SuperClusterPt_DATA_barrel2.GetEntries(): " + str(h1_SuperClusterPt_DATA_barrel2.GetEntries())
# print "-------------------------------------------------------------------------"
# print "h1_SuperClusterPt_DATA_endcap.GetEntries(): " + str(h1_SuperClusterPt_DATA_endcap.GetEntries())
# print "h1_SuperClusterPt_DATA_endcap1.GetEntries(): " + str(h1_SuperClusterPt_DATA_endcap1.GetEntries())
# print "h1_SuperClusterPt_DATA_endcap2.GetEntries(): " + str(h1_SuperClusterPt_DATA_endcap2.GetEntries())

# print "-------------------------------------------------------------------------"
# print "-------------------------------------------------------------------------"

# print "h1_ElePt_DATA_all.GetEntries(): " + str(h1_ElePt_DATA_all.GetEntries())
# print "-------------------------------------------------------------------------"
# print "h1_ElePt_DATA_barrel.GetEntries(): " + str(h1_ElePt_DATA_barrel.GetEntries())
# print "h1_ElePt_DATA_barrel1.GetEntries(): " + str(h1_ElePt_DATA_barrel1.GetEntries())
# print "h1_ElePt_DATA_barrel2.GetEntries(): " + str(h1_ElePt_DATA_barrel2.GetEntries())
# print "-------------------------------------------------------------------------"
# print "h1_ElePt_DATA_endcap.GetEntries(): " + str(h1_ElePt_DATA_endcap.GetEntries())
# print "h1_ElePt_DATA_endcap1.GetEntries(): " + str(h1_ElePt_DATA_endcap1.GetEntries())
# print "h1_ElePt_DATA_endcap2.GetEntries(): " + str(h1_ElePt_DATA_endcap2.GetEntries())
# print "-------------------------------------------------------------------------"

# print "-------------------------------------------------------------------------"
# print "-------------------------------------------------------------------------"

#Calculate fake rate
fakeRate_vs_Pt_DATA_all = GetEffFixBinning( h1_ElePt_DATA_all , h1_SuperClusterPt_DATA_all
                                            , 1.2 , 20 , 1
                                            , "SuperCluster Et (GeV)" , "Fake Probability - All"
                                            , 30 , 500, 10)
# fakeRate_vs_Eta_DATA_all = GetEffFixBinning( h1_EleEta_DATA_all , h1_SuperClusterEta_DATA_all
#                                              , 1.2 , 20 , 1
#                                              , "SuperCluster Eta" , "Fake Probability - All"
#                                              , -5 , 5, 1)

fakeRate_vs_Pt_DATA_barrel = GetEffFixBinning( h1_ElePt_DATA_barrel , h1_SuperClusterPt_DATA_barrel
                                               , 1.2 , 20 , 1
                                               , "SuperCluster Et (GeV)" , "Fake Probability - Barrel"
                                               , 30 , 500, 10)
fakeRate_vs_Pt_DATA_barrel1 = GetEffFixBinning( h1_ElePt_DATA_barrel1 , h1_SuperClusterPt_DATA_barrel1
                                                , 1.2 , 20 , 1
                                                , "SuperCluster Et (GeV)" , "Fake Probability - Barrel1"
                                                , 30 , 500, 10)
fakeRate_vs_Pt_DATA_barrel2 = GetEffFixBinning( h1_ElePt_DATA_barrel2 , h1_SuperClusterPt_DATA_barrel2
                                                , 1.2 , 20 , 1
                                                , "SuperCluster Et (GeV)" , "Fake Probability - Barrel2"
                                                , 30 , 500, 10)

fakeRate_vs_Pt_DATA_endcap = GetEffFixBinning( h1_ElePt_DATA_endcap , h1_SuperClusterPt_DATA_endcap
                                               , 1.2 , 20 , 1
                                               , "SuperCluster Et (GeV)" , "Fake Probability - Endcap"
                                               , 30 , 500, 10)
fakeRate_vs_Pt_DATA_endcap1 = GetEffFixBinning( h1_ElePt_DATA_endcap1 , h1_SuperClusterPt_DATA_endcap1
                                                , 1.2 , 20 , 1
                                                , "SuperCluster Et (GeV)" , "Fake Probability - Endcap1"
                                                , 30 , 500, 10)
fakeRate_vs_Pt_DATA_endcap2 = GetEffFixBinning( h1_ElePt_DATA_endcap2 , h1_SuperClusterPt_DATA_endcap2
                                                , 1.2 , 20 , 1
                                                , "SuperCluster Et (GeV)" , "Fake Probability - Endcap2"
                                                , 30 , 500, 10)


#--- Final plots
cTot = TCanvas()
cTot.Print("FakeRatePlots.ps[")

#fake rate pt (all)
call = TCanvas()
call.SetGridy()
call.SetGridx()
fakeRate_vs_Pt_DATA_all.GetYaxis().SetRangeUser(0,0.2)
fakeRate_vs_Pt_DATA_all.GetXaxis().SetRangeUser(30,1000)
fakeRate_vs_Pt_DATA_all.Draw("ap")
fakeRate_vs_Pt_DATA_all.SetLineWidth(2)
fakeRate_vs_Pt_DATA_all.GetXaxis().SetTitle("SuperCluster Et (GeV)")
fakeRate_vs_Pt_DATA_all.GetYaxis().SetTitle("Fake Probability - All")
call.Print("FakeRate_all.root","root")
call.Print("FakeRatePlots.ps")

#fake rate (barrel)
cbarrel = TCanvas()
cbarrel.SetGridy()
cbarrel.SetGridx()
fakeRate_vs_Pt_DATA_barrel.GetYaxis().SetRangeUser(0,0.2)
fakeRate_vs_Pt_DATA_barrel.GetXaxis().SetRangeUser(30,1000)
fakeRate_vs_Pt_DATA_barrel.Draw("ap")
fakeRate_vs_Pt_DATA_barrel.SetLineWidth(2)
fakeRate_vs_Pt_DATA_barrel.GetXaxis().SetTitle("SuperCluster Et (GeV)")
fakeRate_vs_Pt_DATA_barrel.GetYaxis().SetTitle("Fake Probability - Barrel")
cbarrel.Print("FakeRate_barrel.root","root")
cbarrel.Print("FakeRatePlots.ps")

#fake rate (barrel1)
cbarrel1 = TCanvas()
cbarrel1.SetGridy()
cbarrel1.SetGridx()
fakeRate_vs_Pt_DATA_barrel1.GetYaxis().SetRangeUser(0,0.2)
fakeRate_vs_Pt_DATA_barrel1.GetXaxis().SetRangeUser(30,1000)
fakeRate_vs_Pt_DATA_barrel1.Draw("ap")
fakeRate_vs_Pt_DATA_barrel1.SetLineWidth(2)
fakeRate_vs_Pt_DATA_barrel1.GetXaxis().SetTitle("SuperCluster Et (GeV)")
fakeRate_vs_Pt_DATA_barrel1.GetYaxis().SetTitle("Fake Probability - Barrel1")
cbarrel1.Print("FakeRate_barrel1.root","root")
cbarrel1.Print("FakeRatePlots.ps")

#fake rate (barrel2)
cbarrel2 = TCanvas()
cbarrel2.SetGridy()
cbarrel2.SetGridx()
fakeRate_vs_Pt_DATA_barrel2.GetYaxis().SetRangeUser(0,0.2)
fakeRate_vs_Pt_DATA_barrel2.GetXaxis().SetRangeUser(30,1000)
fakeRate_vs_Pt_DATA_barrel2.Draw("ap")
fakeRate_vs_Pt_DATA_barrel2.SetLineWidth(2)
fakeRate_vs_Pt_DATA_barrel2.GetXaxis().SetTitle("SuperCluster Et (GeV)")
fakeRate_vs_Pt_DATA_barrel2.GetYaxis().SetTitle("Fake Probability - Barrel2")
cbarrel2.Print("FakeRate_barrel2.root","root")
cbarrel2.Print("FakeRatePlots.ps")

#fake rate (endcap)
cendcap = TCanvas()
cendcap.SetGridy()
cendcap.SetGridx()
fakeRate_vs_Pt_DATA_endcap.GetYaxis().SetRangeUser(0,0.2)
fakeRate_vs_Pt_DATA_endcap.GetXaxis().SetRangeUser(30,1000)
fakeRate_vs_Pt_DATA_endcap.Draw("ap")
fakeRate_vs_Pt_DATA_endcap.SetLineWidth(2)
fakeRate_vs_Pt_DATA_endcap.GetXaxis().SetTitle("SuperCluster Et (GeV)")
fakeRate_vs_Pt_DATA_endcap.GetYaxis().SetTitle("Fake Probability - Endcap")
cendcap.Print("FakeRate_endcap.root","root")
cendcap.Print("FakeRatePlots.ps")

#fake rate (endcap1)
cendcap1 = TCanvas()
cendcap1.SetGridy()
cendcap1.SetGridx()
fakeRate_vs_Pt_DATA_endcap1.GetYaxis().SetRangeUser(0,0.2)
fakeRate_vs_Pt_DATA_endcap1.GetXaxis().SetRangeUser(30,1000)
fakeRate_vs_Pt_DATA_endcap1.Draw("ap")
fakeRate_vs_Pt_DATA_endcap1.SetLineWidth(2)
fakeRate_vs_Pt_DATA_endcap1.GetXaxis().SetTitle("SuperCluster Et (GeV)")
fakeRate_vs_Pt_DATA_endcap1.GetYaxis().SetTitle("Fake Probability - Endcap1")
cendcap1.Print("FakeRate_endcap1.root","root")
cendcap1.Print("FakeRatePlots.ps")

#fake rate (endcap2)
cendcap2 = TCanvas()
cendcap2.SetGridy()
cendcap2.SetGridx()
fakeRate_vs_Pt_DATA_endcap2.GetYaxis().SetRangeUser(0,0.2)
fakeRate_vs_Pt_DATA_endcap2.GetXaxis().SetRangeUser(30,1000)
fakeRate_vs_Pt_DATA_endcap2.Draw("ap")
fakeRate_vs_Pt_DATA_endcap2.SetLineWidth(2)
fakeRate_vs_Pt_DATA_endcap2.GetXaxis().SetTitle("SuperCluster Et (GeV)")
fakeRate_vs_Pt_DATA_endcap2.GetYaxis().SetTitle("Fake Probability - Endcap2")
cendcap2.Print("FakeRate_endcap2.root","root")
cendcap2.Print("FakeRatePlots.ps")

# #fake rate eta (all)
# ceta = TCanvas()
# ceta.SetGridy()
# ceta.SetGridx()
# fakeRate_vs_Eta_DATA_all.GetYaxis().SetRangeUser(0,10)
# fakeRate_vs_Eta_DATA_all.GetXaxis().SetRangeUser(0,1)
# fakeRate_vs_Eta_DATA_all.Draw("ap")
# fakeRate_vs_Eta_DATA_all.SetLineWidth(2)
# fakeRate_vs_Eta_DATA_all.GetXaxis().SetTitle("SuperCluster Eta")
# fakeRate_vs_Eta_DATA_all.GetYaxis().SetTitle("Fake Probability - All")
# ceta.Print("FakeRate_all_eta.root","root")
# ceta.Print("FakeRatePlots.ps")


#---
cTot.Print("FakeRatePlots.ps]")

## Terminate the program
#print "Press ENTER to terminate"
#wait=raw_input()
