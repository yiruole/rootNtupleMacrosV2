#!/usr/bin/env python

## Instructions and description:
##     This code uses defined methods to divide histograms and compute the proper error bars (only for the fixed bin case, not yet for the var bin case)
##     for efficiencies from weighted samples.
##     Change the name of the input files to be the output of analysisClass_SCFakeRate.C
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

#--- ROOT general options
gStyle.SetOptStat(0)
gStyle.SetPalette(1)
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1);
gStyle.SetPadTickY(1);
#--- TODO: WHY IT DOES NOT LOAD THE DEFAULT ROOTLOGON.C ? ---#


#--- Root files

# %%%%%%% BEGIN %%%%%%%
def GetFile(filename):
    file = TFile(filename)
    if( not file):
        print "ERROR: file " + filename + " not found"
        print "exiting..."
        sys.exit()
    return file
# %%%%%%% END %%%%%%%


#--- Define functions

# %%%%%%% BEGIN %%%%%%%     

def GetHisto( histoName , file ):
    histo = file.Get( histoName )
    if( not histo):
        print "ERROR: histo " + histoName + " not found in " + file.GetName()
        print "exiting..."
        sys.exit()
    return histo



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
    GraphEff.Fit("f1","R")

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


def GetHistoRescaled( histo , function ):

    Bins = []
    BinsContent = []
    BinsContentErrMax = []
    for point in range( 0 , function.GetN() ):
        x, y = ROOT.Double(1), ROOT.Double(1)
        function.GetPoint(point, x, y)
        ehx = function.GetErrorXhigh(point)  
        elx = function.GetErrorXlow(point)
        ehy = function.GetErrorYhigh(point)
        ely = function.GetErrorYlow(point)  
        up = x + ehx
        low = x - elx
        if(point == 0):
            Bins.append( low )
        Bins.append( up )
        BinsContent.append( y )
        BinsContentErrMax.append( max(ehy,ely) )

    BinsFinal = array( 'f', Bins ) 

    ##define the tmp-histo
    histoFinal = TH1F("histoFinal",
                      histo.GetName()+"_rescaled",
                      len(Bins) - 1 , BinsFinal )
    histoFinal.Sumw2()

    ##fill the tmp-histo with the content of the histo passed as argument to the function
    for bin1 in range( 1 , histoFinal.GetNbinsX()+1 ):
        center = histoFinal.GetBinCenter(bin1)
        halfwidth = histoFinal.GetBinWidth(bin1) / 2
        #print center
        #print halfwidth
        newValue = 0
        newErr = 0
        for bin2 in range( 1 , histo.GetNbinsX()+1 ):
            if(histo.GetBinCenter(bin2) > (center - halfwidth)
               and
               histo.GetBinCenter(bin2) < (center + halfwidth) ):
                newValue = newValue + histo.GetBinContent(bin2)
                newErr = newErr + histo.GetBinError(bin2)*histo.GetBinError(bin2)

        histoFinal.SetBinContent(bin1,newValue)
        histoFinal.SetBinError(bin1,sqrt(newErr))

#    for bin in range( 1 , histo.GetNbinsX()+1 ):
#        for entry in range( 0 , int(histo.GetBinContent(bin)) ):
#            histoFinal.SetBinContent( histo.GetBinCenter(bin) , float(histo.GetBinContent(bin) )

    ##rescale tmp-histo with the function
    for index, rescaleFactor in enumerate( BinsContent ):
        #print index
        #print "histoFinal.GetBinContent(index+1): " + str(histoFinal.GetBinContent(index+1))
        #print "histoFinal.GetBinError(index+1): " + str(histoFinal.GetBinError(index+1))
        #print "rescaleFactor: " + str(rescaleFactor)
        #print "BinsContentErrMax[index]: " + str(BinsContentErrMax[index])
        rescaledBinContent = rescaleFactor * histoFinal.GetBinContent(index+1)
        if( histoFinal.GetBinContent(index+1) >0 and rescaleFactor >0 ):
            erelh = histoFinal.GetBinError(index+1) / histoFinal.GetBinContent(index+1)
            erelf = BinsContentErrMax[index] / rescaleFactor
            erelhfin = sqrt(erelh*erelh + erelf*erelf)
            rescaledBinError = erelhfin * rescaledBinContent
            #print "rescaledBinContent: " + str(rescaledBinContent)
            #print "rescaledBinError: " + str(rescaledBinError)
            histoFinal.SetBinContent(index+1 , rescaledBinContent )
            histoFinal.SetBinError(index+1 , rescaledBinError )
        else:
            histoFinal.SetBinContent(index+1 , 0 )
            histoFinal.SetBinError(index+1 , 0 )
        #print ""

    return histoFinal

def GetIntegral(histo):
    return histo.Integral()

def GetIntegralError(histo):
    integralError = 0
    for bin in range( 1 , histo.GetNbinsX()+1 ):
        integralError = integralError + (histoFinal.GetBinError(bin)*histoFinal.GetBinError(bin))
    integralError = sqrt(integralError)
    return integralError
        
# %%%%%%% END %%%%%%% 

                #output from /rootNtupleMacros/src/analysisClass_elecStudies2.C

#File2 = GetFile("$LQDATA/qcd/2.9pb-1/output_cutTable_SCFakeRate_fromEllie_cutCvs1.10_macroCvs1.26/analysisClass_SCFakeRate_plots.root")
#File2 = GetFile("$LQDATA/qcd/2.9pb-1/output_cutTable_SCFakeRate_Pt1stSCgt25/analysisClass_SCFakeRate_plots.root")
# File2 = GetFile("$LQDATA/qcd/2.9pb-1/output_cutTable_SCFakeRate_Pt1stSCgt25_nJetPtCut4orMore_MET20/analysisClass_SCFakeRate_plots.root")
# File2 = GetFile("$LQDATA/qcd/2.9pb-1/output_cutTable_SCFakeRate_Pt1stSCgt25_nJetPtCut3orMore_MET20/analysisClass_SCFakeRate_plots.root")
# File2 = GetFile("$LQDATA/qcd/2.9pb-1/output_cutTable_SCFakeRate_Pt1stSCgt25_nJetPtCut2orMore_MET20/analysisClass_SCFakeRate_plots.root")
# File2 = GetFile("$LQDATA/qcd/2.9pb-1/output_cutTable_SCFakeRate_ele_PtCut25/analysisClass_SCFakeRate_plots.root")
# File2 = GetFile("$LQDATA/qcd/2.47pb-1/output_cutTable_SCFakeRate_Pt1stSCgt30_nJetPtCut2orMore_MET20_noDEtaInEE/analysisClass_SCFakeRate_plots.root")
# File2 = GetFile("$LQDATA/qcd/2.47pb-1/output_cutTable_SCFakeRate_Pt1stSCgt30_nJetPtCut2orMore_MET20_withDEtaInEE/analysisClass_SCFakeRate_plots.root")
# File2 = GetFile("$LQDATA/qcd/2.47pb-1/output_cutTable_SCFakeRate_Pt1stSCgt30_nJetPtCut3orMore_MET20_withDEtaInEE/analysisClass_SCFakeRate_plots.root")
#File2 = GetFile("$LQDATA/qcd/2010dataset/output_cutTable_SCFakeRate_allruns/analysisClass_SCFakeRate_plots.root")
File2 = GetFile("$LQDATA/qcd/2010dataset/output_cutTable_SCFakeRate_run_le_144114/analysisClass_SCFakeRate_plots.root")
#File2 = GetFile("$LQDATA/qcd/2010dataset/output_cutTable_SCFakeRate_run_ge144114_lt148058/analysisClass_SCFakeRate_plots.root")
#File2 = GetFile("$LQDATA/qcd/2010dataset/output_cutTable_SCFakeRate_run_gt148058/analysisClass_SCFakeRate_plots.root")


#--- Define all the histograms
MyBins = [10,15,20,25,30,40,60,80,100,200,300,400,500]



#PT ratios
h_ele_pt_bottom_NewData_tight_barrel = GetHisto( "histo1D__DATA__goodSCPt_Barrel" , File2)
h_ele_pt_heep_NewData_tight_barrel = GetHisto( "histo1D__DATA__goodEleSCPt_Barrel" , File2)
fakeRate_NewData_tight_barrel = GetEffFixBinning( h_ele_pt_heep_NewData_tight_barrel , h_ele_pt_bottom_NewData_tight_barrel
                                 , 1.2 , 20 , 1
                                 , "Supercluster Pt (GeV)" , "fake probability - Barrel"
                                 , 30 , 500, 5)


h_ele_pt_bottom_NewData_tight_endcap = GetHisto( "histo1D__DATA__goodSCPt_Endcap" , File2)
h_ele_pt_heep_NewData_tight_endcap = GetHisto( "histo1D__DATA__goodEleSCPt_Endcap" , File2)
fakeRate_NewData_tight_endcap = GetEffFixBinning( h_ele_pt_heep_NewData_tight_endcap , h_ele_pt_bottom_NewData_tight_endcap
                                 , 1.2 , 20 , 1
                                 , "Supercluster Pt (GeV)" , "fake probability - Endcap"
                                 , 30 , 500, 5)

# #PT ratios
# h_ele_pt_bottom_NewData_tight_barrel = GetHisto( "histo1D__PhotonJet_Pt15__goodSCPt_Barrel" , File2)
# h_ele_pt_heep_NewData_tight_barrel = GetHisto( "histo1D__PhotonJet_Pt15__goodEleSCPt_Barrel" , File2)
# fakeRate_NewData_tight_barrel = GetEffFixBinning( h_ele_pt_heep_NewData_tight_barrel , h_ele_pt_bottom_NewData_tight_barrel
#                                  , 1.2 , 20 , 1
#                                  , "Supercluster Pt (GeV)" , "fake probability - Barrel"
#                                  , 30 , 100, 10)


# h_ele_pt_bottom_NewData_tight_endcap = GetHisto( "histo1D__PhotonJet_Pt15__goodSCPt_Endcap" , File2)
# h_ele_pt_heep_NewData_tight_endcap = GetHisto( "histo1D__PhotonJet_Pt15__goodEleSCPt_Endcap" , File2)
# fakeRate_NewData_tight_endcap = GetEffFixBinning( h_ele_pt_heep_NewData_tight_endcap , h_ele_pt_bottom_NewData_tight_endcap
#                                  , 1.2 , 20 , 1
#                                  , "Supercluster Pt (GeV)" , "fake probability - Endcap"
#                                  , 30 , 100, 10)


#--- Final plots
cAll = TCanvas();
cAll.Print("FakeRatePlots_fit.ps[");

#all fake rate

c2 = TCanvas()
c2.SetGridy();
c2.SetGridx();
#fakeRate_NewData_tight_barrel.GetYaxis().SetRangeUser(0,0.2)
fakeRate_NewData_tight_barrel.GetYaxis().SetRangeUser(0,0.05)
fakeRate_NewData_tight_barrel.GetXaxis().SetRangeUser(30,500)
fakeRate_NewData_tight_barrel.Draw("ap")
fakeRate_NewData_tight_barrel.SetLineWidth(2)
fakeRate_NewData_tight_barrel.GetXaxis().SetTitle("Supercluster Pt (GeV)")
fakeRate_NewData_tight_barrel.GetYaxis().SetTitle("fake probability TIGHT - Barrel")

## c2.Print("FakeRate_tight_barrel_final.pdf","pdf");
c2.Print("FakeRate_tight_barrel_fit.eps","eps");
c2.Print("FakeRatePlots_fit.ps");

c3 = TCanvas()
c3.SetGridy();
c3.SetGridx();
#fakeRate_NewData_tight_endcap.GetYaxis().SetRangeUser(0,0.2)
fakeRate_NewData_tight_endcap.GetYaxis().SetRangeUser(0,0.20)
fakeRate_NewData_tight_endcap.GetXaxis().SetRangeUser(30,500)
fakeRate_NewData_tight_endcap.Draw("ap")
fakeRate_NewData_tight_endcap.SetLineWidth(2)
fakeRate_NewData_tight_endcap.GetXaxis().SetTitle("Supercluster Pt (GeV)")
fakeRate_NewData_tight_endcap.GetYaxis().SetTitle("fake probability TIGHT - Endcap")

## c3.Print("FakeRate_endcap_final.gif","gif");
c3.Print("FakeRate_tight_endcap_fit.eps","eps");
c3.Print("FakeRatePlots_fit.ps");


## slope = 0.0
## intercept = 0.00715202
## x_line=array("d",[20.0,100.0])
## y_line= array("d",[slope*x_line[0]+intercept,slope*x_line[1]+intercept])
## line = TPolyLine(2,x_line,y_line,"L")
## line.SetLineWidth(2)
## line.Draw()

## slope = 0.0
## intercept = 0.006245
## x_line_sigPlus=array("d",[30.0,100.0])
## y_line_sigPlus= array("d",[slope*x_line_sigPlus[0]+intercept,slope*x_line_sigPlus[1]+intercept])
## line_sigPlus = TPolyLine(2,x_line_sigPlus,y_line_sigPlus,"L")
## line_sigPlus.SetLineColor(2)
## line_sigPlus.SetLineWidth(2)
## line_sigPlus.SetLineStyle(5)
## line_sigPlus.Draw()

## slope = 0.0
## intercept = 0.0054311
## x_line_sigMinus=array("d",[30.0,100.0])
## y_line_sigMinus= array("d",[slope*x_line_sigMinus[0]+intercept,slope*x_line_sigMinus[1]+intercept])
## line_sigMinus = TPolyLine(2,x_line_sigMinus,y_line_sigMinus,"L")
## line_sigMinus.SetLineColor(2)
## line_sigMinus.SetLineWidth(2)
## line_sigMinus.SetLineStyle(5)
## line_sigMinus.Draw()


## c2 = TCanvas()
## c2.SetGridy();
## c2.SetGridx();
## fakeRate_qcd_eta.GetYaxis().SetRangeUser(0,0.08);
## fakeRate_qcd_eta.Draw("ap")
## fakeRate_data_eta.Draw("pSAME")
## ## legend.Draw()
## c2.Print("FakeRate_eta_data.gif","gif");
## c2.Print("FakeRate_eta_data.pdf","pdf");
## c2.Print("FakeRate_eta_data.eps","eps");
## c2.Print("FakeRatePlots.ps");

## c1 = TCanvas()
## c1.SetGridy();
## c1.SetGridx();
## fakeRate_eta.Draw("ap")
## ## fakeRate_eta_2.Draw("pSAME")
## c1.Print("FakeRate_eta.gif","gif");
## c1.Print("FakeRatePlots.ps");

## c2 = TCanvas()
## c2.SetGridy();
## c2.SetGridx();
## fakeRate_ST.Draw("ap")
## c2.Print("FakeRate_ST.gif","gif");
## c2.Print("FakeRatePlots.ps");

## ### 1_ele fake rate
## line1 = TLine();
## line1.SetLineStyle(7);
## line1.DrawLine(50,0.04437,300,0.135173);





cAll.Print("FakeRatePlots_fit.ps]");


# lasthisto = File5.Get( "histo1D__QCD_MdGrph__cutHisto_allOtherCuts___________sT" )
# lasthisto.GetXaxis().SetTitle("St")
# c4 = TCanvas()
# c4.SetGridy();
# c4.SetGridx();
# lasthisto.Draw("ap")


## Terminate the program
print "Press ENTER to terminate"
wait=raw_input()
