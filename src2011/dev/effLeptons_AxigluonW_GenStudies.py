#!/usr/bin/env python

##############################################################################
## DONT'T MODIFY WITHIN "# %%%%%%% BEGIN %%%%%%%"  and "# %%%%%%% END %%%%%%%"
##############################################################################

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
gStyle.SetPadTopMargin(0.08);
gStyle.SetPadBottomMargin(0.12);
#gStyle.SetTitleSize(0.05, "XYZ");
#--- TODO: WHY IT DOES NOT LOAD THE DEFAULT ROOTLOGON.C ? ---#

#--- Define functions

# %%%%%%% BEGIN %%%%%%%     

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

def GetEffFixBinning( num , den , m_size , m_style , m_color , xtitle , ytitle , min , max, EffMin, EffMax):
    GraphEffTemp = TGraphAsymmErrors( num , den )
    GraphEff = TGraphAsymmErrors()
    counter=0
    for point in range( 0 , GraphEffTemp.GetN() ):
        x = ROOT.Double(1)
        y = ROOT.Double(1)
        GraphEffTemp.GetPoint( point , x , y) 
        if( x > max or x < min ):
            continue
        GraphEff.SetPoint(counter,x,y);
        GraphEff.SetPointError(counter,
                               GraphEffTemp.GetErrorXlow(point),
                               GraphEffTemp.GetErrorXhigh(point),
                               GraphEffTemp.GetErrorYlow(point),
                               GraphEffTemp.GetErrorYhigh(point) )
        counter = counter + 1

    GraphEff.SetMarkerSize( m_size )
    GraphEff.SetMarkerStyle( m_style )
    GraphEff.SetMarkerColor( m_color )
    GraphEff.GetXaxis().SetTitle(xtitle)
    GraphEff.GetYaxis().SetTitle(ytitle)

    GraphEff.GetYaxis().SetRangeUser(EffMin,EffMax)

    return GraphEff



def GetEffVarBinning( num , den , m_size , m_style , m_color , xtitle , ytitle , min , max, Bins, EffMin, EffMax):

    NbinTot = len(Bins) 
    print Bins
    BinsFinal = array( 'f', Bins ) 

    num_varBin = TH1F("num_varBin",
                      "num_varBin",
                      NbinTot - 1 , BinsFinal );
    nBins_num = num.GetXaxis().GetNbins()
    getxmax_num = num.GetXaxis().GetXmax();
    getxmin_num = num.GetXaxis().GetXmin();
    step_num = (getxmax_num - getxmin_num) / nBins_num;
    for bin in range( 1 , nBins_num+1 ):
        for entry in range( 0 , int(num.GetBinContent(bin)) ):            
            num_varBin.Fill(step_num * bin)

    den_varBin = TH1F("den_varBin",
                      "den_varBin",
                      NbinTot - 1 , BinsFinal );
    nBins_den = den.GetXaxis().GetNbins()
    getxmax_den = den.GetXaxis().GetXmax();
    getxmin_den = den.GetXaxis().GetXmin();
    step_den = (getxmax_den - getxmin_den) / nBins_den;
    for bin in range( 1 , nBins_den+1 ):
        for entry in range( 0 , int(den.GetBinContent(bin)) ):
            den_varBin.Fill(step_den * bin)


    #--- create final graph
    #    print num_varBin.GetXaxis().GetNbins()
    #    print den_varBin.GetXaxis().GetNbins()
    #    print num_varBin.GetXaxis().GetXmax()
    #    print den_varBin.GetXaxis().GetXmax()
    #    print num_varBin.GetXaxis().GetXmin()
    #    print den_varBin.GetXaxis().GetXmin()
    #    num_varBin.SaveAs("num.root")
    #    den_varBin.SaveAs("den.root")
    GraphEffTemp = TGraphAsymmErrors( num_varBin , den_varBin )
    GraphEff = TGraphAsymmErrors()
    counter=0
    for point in range( 0 , GraphEffTemp.GetN() ):
        x = ROOT.Double(1)
        y = ROOT.Double(1)
        GraphEffTemp.GetPoint( point , x , y) 
        if( x > max or x < min ):
            continue
        GraphEff.SetPoint(counter,x,y);
        GraphEff.SetPointError(counter,
                               GraphEffTemp.GetErrorXlow(point),
                               GraphEffTemp.GetErrorXhigh(point),
                               GraphEffTemp.GetErrorYlow(point),
                               GraphEffTemp.GetErrorYhigh(point) )
        counter = counter + 1

    GraphEff.SetMarkerSize( m_size )
    GraphEff.SetMarkerStyle( m_style )
    GraphEff.SetMarkerColor( m_color )
    GraphEff.GetXaxis().SetTitle(xtitle)
    GraphEff.GetYaxis().SetTitle(ytitle)

    GraphEff.GetYaxis().SetRangeUser(EffMin,EffMax)

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


#--- Root files

File1 = GetFile("AG150WGen.root")
#File1 = GetFile("AG1000WGen.root")

#--- Define all the histograms

#--------- Electrons ----------
# outputFile = "effEle_AG150WGen.ps"
# #outputFile = "effEle_AG1000WGen.ps"
# leptonType = "electron"
# ## Eta histograms
# h_lepton_eta_gen_all = GetHisto( "h1_Eta_W_e" , File1)
# h_lepton_eta_gen_cuts = GetHisto( "h1_Eta_W_e_gencuts" , File1)
# h_lepton_eta_gen_RECO = GetHisto( "h1_Eta_W_e_recomatched" , File1)
# h_lepton_eta_gen_RECO_ID_ISO = GetHisto( "h1_Eta_W_e_recomatched_IDISO" , File1)
# ## Pt histograms
# h_lepton_pt_gen_all = GetHisto( "h1_Pt_W_e" , File1)
# h_lepton_pt_gen_cuts = GetHisto( "h1_Pt_W_e_gencuts" , File1)
# h_lepton_pt_gen_RECO = GetHisto( "h1_Pt_W_e_recomatched" , File1)
# h_lepton_pt_gen_RECO_ID_ISO = GetHisto( "h1_Pt_W_e_recomatched_IDISO" , File1)
# ## DR histograms
# h_lepton_DRleptonquark_all = GetHisto( "h1_minDR_elequark" , File1)
# h_lepton_DRleptonquark_cuts = GetHisto( "h1_minDR_elequark_gencuts" , File1)
# h_lepton_DRleptonquark_RECO = GetHisto( "h1_minDR_elequark_recomatched" , File1)
# h_lepton_DRleptonquark_RECO_ID_ISO= GetHisto( "h1_minDR_elequark_recomatched_IDISO" , File1)
# ## nVtx histograms
# h_lepton_nVtx_gen_all = GetHisto( "h1_nVtx_e" , File1)
# h_lepton_nVtx_gen_cuts = GetHisto( "h1_nVtx_e_gencuts" , File1)
# h_lepton_nVtx_gen_RECO = GetHisto( "h1_nVtx_e_recomatched" , File1)
# h_lepton_nVtx_gen_RECO_ID_ISO = GetHisto( "h1_nVtx_e_recomatched_IDISO" , File1)


#--------- Muons ----------
outputFile = "effMuon_AG150WGen.ps"
#outputFile = "effMuon_AG1000WGen.ps"
leptonType = "muon"
## Eta histograms
h_lepton_eta_gen_all = GetHisto( "h1_Eta_W_mu" , File1)
h_lepton_eta_gen_cuts = GetHisto( "h1_Eta_W_mu_gencuts" , File1)
h_lepton_eta_gen_RECO = GetHisto( "h1_Eta_W_mu_recomatched" , File1)
h_lepton_eta_gen_RECO_ID_ISO = GetHisto( "h1_Eta_W_mu_recomatched_IDISO" , File1)
## Pt histograms
h_lepton_pt_gen_all = GetHisto( "h1_Pt_W_mu" , File1)
h_lepton_pt_gen_cuts = GetHisto( "h1_Pt_W_mu_gencuts" , File1)
h_lepton_pt_gen_RECO = GetHisto( "h1_Pt_W_mu_recomatched" , File1)
h_lepton_pt_gen_RECO_ID_ISO = GetHisto( "h1_Pt_W_mu_recomatched_IDISO" , File1)
## DR histograms
h_lepton_DRleptonquark_all = GetHisto( "h1_minDR_muonquark" , File1)
h_lepton_DRleptonquark_cuts = GetHisto( "h1_minDR_muonquark_gencuts" , File1)
h_lepton_DRleptonquark_RECO = GetHisto( "h1_minDR_muonquark_recomatched" , File1)
h_lepton_DRleptonquark_RECO_ID_ISO= GetHisto( "h1_minDR_muonquark_recomatched_IDISO" , File1)
## nVtx histograms
h_lepton_nVtx_gen_all = GetHisto( "h1_nVtx_mu" , File1)
h_lepton_nVtx_gen_cuts = GetHisto( "h1_nVtx_mu_gencuts" , File1)
h_lepton_nVtx_gen_RECO = GetHisto( "h1_nVtx_mu_recomatched" , File1)
h_lepton_nVtx_gen_RECO_ID_ISO = GetHisto( "h1_nVtx_mu_recomatched_IDISO" , File1)


############################
#variable binning
MyBins = [0,10,20,30,40,50,60,70,80,90,100,150,200,300]
############################


#--- Calculate Acceptance
accept_lepton_eta_gen = GetEffFixBinning( h_lepton_eta_gen_cuts , h_lepton_eta_gen_all
                                       , 0.9 , 20 , 1
                                       , leptonType+" \\eta_{gen}" , "acceptance"
                                       , -2.5 , 2.5, 0, 1)

accept_lepton_pt_gen = GetEffVarBinning( h_lepton_pt_gen_cuts , h_lepton_pt_gen_all
                                      , 0.9 , 20 , 1
                                      , leptonType+" pT_{gen} (GeV)" , "acceptance"
                                      , 0 , 300
                                      , MyBins, 0.5, 1)

accept_lepton_DRleptonquark_gen = GetEffFixBinning( h_lepton_DRleptonquark_cuts , h_lepton_DRleptonquark_all
                                              , 0.9 , 20 , 1
                                              , "min\\Delta R (lepton,quark)" , "acceptance"
                                              , 0 , 10, 0, 1)


#--- Calculate Absolute Efficiency
eff_lepton_eta_gen_RECO = GetEffFixBinning( h_lepton_eta_gen_RECO , h_lepton_eta_gen_cuts
                                           , 0.9 , 20 , 1
                                           , leptonType+" \\eta_{gen}" , "relative efficiency wrt to acceptance"
                                           , -2.5 , 2.5, 0, 1)

eff_lepton_eta_gen_RECO_ID_ISO = GetEffFixBinning( h_lepton_eta_gen_RECO_ID_ISO , h_lepton_eta_gen_cuts
                                                , 0.9 , 20 , 4
                                                , leptonType+" \\eta_{gen}" , "relative efficiency wrt to acceptance"
                                                , -2.5 , 2.5, 0, 1)

eff_lepton_pt_gen_RECO = GetEffVarBinning( h_lepton_pt_gen_RECO , h_lepton_pt_gen_cuts
                                        , 0.9 , 20 , 1
                                        , leptonType+" pT_{gen} (GeV)" , "relative efficiency wrt to acceptance"
                                        , 0 , 300
                                        , MyBins, 0.5, 1)

eff_lepton_pt_gen_RECO_ID_ISO = GetEffVarBinning( h_lepton_pt_gen_RECO_ID_ISO , h_lepton_pt_gen_cuts
                                               , 0.9 , 20 , 4
                                               , leptonType+" pT_{gen} (GeV)" , "relative efficiency wrt to acceptance"
                                               , 0 , 300
                                               , MyBins, 0.5, 1)

eff_lepton_DRleptonquark_gen_RECO = GetEffFixBinning( h_lepton_DRleptonquark_RECO , h_lepton_DRleptonquark_cuts
                                                , 0.9 , 20 , 1
                                                , "min\\Delta R ("+leptonType+",quark)" , "relative efficiency wrt to acceptance"
                                                , 0 , 10, 0, 1)

eff_lepton_DRleptonquark_gen_RECO_ID_ISO = GetEffFixBinning( h_lepton_DRleptonquark_RECO_ID_ISO , h_lepton_DRleptonquark_cuts
                                                       , 0.9 , 20 , 4
                                                       , "min\\Delta R ("+leptonType+",quark)" , "relative efficiency wrt to acceptance"
                                                       , 0 , 10, 0, 1)


#--- Calculate Relative Efficiency


eff_lepton_pt_gen_RECO_ID_ISO_rel = GetEffVarBinning( h_lepton_pt_gen_RECO_ID_ISO , h_lepton_pt_gen_RECO
                                                   , 0.9 , 20 , 4
                                                   , leptonType+" pT_{gen} (GeV)" , "relative efficiency"
                                                   , 0 , 300
                                                   , MyBins, 0.5, 1)

eff_lepton_DRleptonquark_gen_RECO_ID_ISO_rel = GetEffFixBinning( h_lepton_DRleptonquark_RECO_ID_ISO , h_lepton_DRleptonquark_RECO
                                                           , 0.9 , 20 , 4
                                                           , "min\\Delta R ("+leptonType+",quark)" , "relative efficiency"
                                                           , 0 , 10 , 0., 1)


#--- Calculate total acceptance x efficiency

eff_times_acc_lepton_nVtx_lepton_gen_RECO_ID_ISO = GetEffFixBinning( h_lepton_nVtx_gen_RECO_ID_ISO , h_lepton_nVtx_gen_all
                                                                     , 0.9 , 20 , 4
                                                                     , "Number of reco vertices (in events with "+leptonType+")" , "acceptance x efficiency"
                                                                     , 0 , 30 , 0., 1)

#--- Final plots

#------------------------------------------------------

## Eta distributions - lepton
c0 = TCanvas()
h_lepton_eta_gen_all.SetMarkerStyle(20)
h_lepton_eta_gen_all.SetMarkerColor(kBlack)
h_lepton_eta_gen_all.SetTitle("")
h_lepton_eta_gen_all.GetXaxis().SetTitle(leptonType+" \\eta_{gen}")
h_lepton_eta_gen_all.GetYaxis().SetTitle("arbitrary unit")
h_lepton_eta_gen_all.Draw("p")
h_lepton_eta_gen_cuts.SetMarkerStyle(22)
h_lepton_eta_gen_cuts.SetMarkerColor(kRed)
h_lepton_eta_gen_cuts.Draw("psame")

#---Create legend
globals()['legend0'] = TLegend(0.564648,0.6618,0.869067,0.851582)
legend0.SetFillColor(kWhite)
legend0.SetMargin(0.3)
legend0.AddEntry(h_lepton_eta_gen_all,"no cuts","p")
legend0.AddEntry(h_lepton_eta_gen_cuts,"acceptance cuts","p")
legend0.Draw()
c0.Update()
gPad.RedrawAxis()
gPad.Modified()

c0.Print(outputFile+"(")

#------------------------------------------------------

## pt distributions - lepton
c00 = TCanvas()
h_lepton_pt_gen_all.SetMarkerStyle(20)
h_lepton_pt_gen_all.SetMarkerColor(kBlack)
h_lepton_pt_gen_all.SetTitle("")
h_lepton_pt_gen_all.GetXaxis().SetTitle(leptonType+" pt_{gen} (GeV)")
h_lepton_pt_gen_all.GetYaxis().SetTitle("arbitrary unit")
h_lepton_pt_gen_all.Draw("p")
h_lepton_pt_gen_cuts.SetMarkerStyle(22)
h_lepton_pt_gen_cuts.SetMarkerColor(kRed)
h_lepton_pt_gen_cuts.Draw("psame")

#---Create legend
globals()['legend00'] = TLegend(0.564648,0.6618,0.869067,0.851582)
legend00.SetFillColor(kWhite)
legend00.SetMargin(0.3)
legend00.AddEntry(h_lepton_pt_gen_all,"no cuts","p")
legend00.AddEntry(h_lepton_pt_gen_cuts,"acceptance cuts","p")
legend00.Draw()
c00.Update()
gPad.RedrawAxis()
gPad.Modified()

c00.Print(outputFile)

#------------------------------------------------------

## DR distributions - lepton
c01 = TCanvas()

h_lepton_DRleptonquark_all.SetMarkerStyle(20)
h_lepton_DRleptonquark_all.SetMarkerColor(kBlack)
h_lepton_DRleptonquark_all.SetTitle("")
h_lepton_DRleptonquark_all.GetXaxis().SetTitle("min\\Delta R ("+leptonType+",quark)")
h_lepton_DRleptonquark_all.GetYaxis().SetTitle("arbitrary unit")
h_lepton_DRleptonquark_all.Draw("p")
h_lepton_DRleptonquark_cuts.SetMarkerStyle(22)
h_lepton_DRleptonquark_cuts.SetMarkerColor(kRed)
h_lepton_DRleptonquark_cuts.Draw("psame")

#---Create legend
globals()['legend01'] = TLegend(0.564648,0.6618,0.869067,0.851582)
legend01.SetFillColor(kWhite)
legend01.SetMargin(0.3)
legend01.AddEntry(h_lepton_DRleptonquark_all,"no cuts","p")
legend01.AddEntry(h_lepton_DRleptonquark_cuts,"acceptance cuts","p")
legend01.Draw()
c01.Update()
gPad.RedrawAxis()
gPad.Modified()

c01.Print(outputFile)

#------------------------------------------------------

## Acceptance vs eta - lepton
c1 = TCanvas()
c1.SetGridy();
c1.SetGridx();
accept_lepton_eta_gen.Draw("ap")

c1.Print(outputFile)

## Acceptance vs pt - lepton
c2 = TCanvas()
c2.SetGridy();
c2.SetGridx();
accept_lepton_pt_gen.Draw("ap")

c2.Print(outputFile)

## Acceptance vs DR - lepton
c21 = TCanvas()
c21.SetGridy();
c21.SetGridx();
accept_lepton_DRleptonquark_gen.Draw("ap")

c21.Print(outputFile)

#------------------------------------------------------

## Eff vs eta - lepton
c3 = TCanvas()
c3.SetGridy();
c3.SetGridx();
eff_lepton_eta_gen_RECO.Draw("ap")
eff_lepton_eta_gen_RECO_ID_ISO.Draw("psame")

#---Create legend
globals()['legend3'] = TLegend(0.348609,0.184915,0.653028,0.377129)
legend3.SetFillColor(kWhite)
legend3.SetMargin(0.3)
legend3.AddEntry(eff_lepton_eta_gen_RECO,"only reco","p")
legend3.AddEntry(eff_lepton_eta_gen_RECO_ID_ISO,"reco+ID+ISO","p")
legend3.Draw()
c3.Update()
gPad.RedrawAxis()
gPad.Modified()

c3.Print(outputFile)

#------------------------------------------------------

## Eff vs pt - lepton
c4 = TCanvas()
c4.SetGridy();
c4.SetGridx();
eff_lepton_pt_gen_RECO.Draw("ap")
eff_lepton_pt_gen_RECO_ID_ISO.Draw("psame")

#---Create legend
globals()['legend4'] = TLegend(0.348609,0.184915,0.653028,0.377129)
legend4.SetFillColor(kWhite)
legend4.SetMargin(0.3)
legend4.AddEntry(eff_lepton_pt_gen_RECO,"only reco","p")
legend4.AddEntry(eff_lepton_pt_gen_RECO_ID_ISO,"reco+ID+ISO","p")
legend4.Draw()
c4.Update()
gPad.RedrawAxis()
gPad.Modified()

c4.Print(outputFile)

#------------------------------------------------------

## Eff vs DR - lepton
c41 = TCanvas()
c41.SetGridy();
c41.SetGridx();
eff_lepton_DRleptonquark_gen_RECO.Draw("ap")
eff_lepton_DRleptonquark_gen_RECO_ID_ISO.Draw("psame")

#---Create legend
globals()['legend41'] = TLegend(0.348609,0.184915,0.653028,0.377129)
legend41.SetFillColor(kWhite)
legend41.SetMargin(0.3)
legend41.AddEntry(eff_lepton_DRleptonquark_gen_RECO,"only reco","p")
legend41.AddEntry(eff_lepton_DRleptonquark_gen_RECO_ID_ISO,"reco+ID+ISO","p")
legend41.Draw()
c41.Update()
gPad.RedrawAxis()
gPad.Modified()

c41.Print(outputFile)

#------------------------------------------------------

## Relative Eff vs pt - lepton
c5 = TCanvas()
c5.SetGridy();
c5.SetGridx();
eff_lepton_pt_gen_RECO_ID_ISO_rel.Draw("ap")

#---Create legend
globals()['legend5'] = TLegend(0.348609,0.184915,0.653028,0.377129)
legend5.SetFillColor(kWhite)
legend5.SetMargin(0.3)
legend5.AddEntry(eff_lepton_pt_gen_RECO_ID_ISO_rel,"ID+ISO eff wrt to RECO","p")
legend5.Draw()
c5.Update()
gPad.RedrawAxis()
gPad.Modified()

c5.Print(outputFile)

#------------------------------------------------------

## Relative Eff vs DR - lepton
c51 = TCanvas()
c51.SetGridy();
c51.SetGridx();
eff_lepton_DRleptonquark_gen_RECO_ID_ISO_rel.Draw("ap")

#---Create legend
globals()['legend51'] = TLegend(0.348609,0.184915,0.653028,0.377129)
legend51.SetFillColor(kWhite)
legend51.SetMargin(0.3)
legend51.AddEntry(eff_lepton_DRleptonquark_gen_RECO_ID_ISO_rel,"ID+ISO eff wrt to RECO","p")
legend51.Draw()
c51.Update()
gPad.RedrawAxis()
gPad.Modified()

c51.Print(outputFile)

#------------------------------------------------------

## Eff x Acceptance vs nVtx - lepton
c52 = TCanvas()
c52.SetGridy();
c52.SetGridx();
eff_times_acc_lepton_nVtx_lepton_gen_RECO_ID_ISO.Draw("ap")

#---Create legend
#globals()['legend52'] = TLegend(0.348609,0.584915,0.653028,0.777129)
globals()['legend52'] = TLegend(0.348609,0.184915,0.653028,0.377129)
legend52.SetFillColor(kWhite)
legend52.SetMargin(0.3)
legend52.AddEntry(eff_times_acc_lepton_nVtx_lepton_gen_RECO_ID_ISO,"Acc. x RECO+ID+ISO Eff.","p")
legend52.Draw()
c52.Update()
gPad.RedrawAxis()
gPad.Modified()

c52.Print(outputFile+")")

#------------------------------------------------------


#Print out
print "Acceptance = " + str(h_lepton_eta_gen_cuts.Integral()/h_lepton_eta_gen_all.Integral())
print "Reco eff wrt to acceptance = " + str(h_lepton_eta_gen_RECO.Integral()/h_lepton_eta_gen_cuts.Integral())
print "Reco+ID+ISO eff wrt to acceptance = " + str(h_lepton_eta_gen_RECO_ID_ISO.Integral()/h_lepton_eta_gen_cuts.Integral())
print "Total Reco+ID+ISO eff = " + str(h_lepton_eta_gen_RECO_ID_ISO.Integral()/h_lepton_eta_gen_all.Integral())

## Terminate the program
print "Press ENTER to terminate"
wait=raw_input()
