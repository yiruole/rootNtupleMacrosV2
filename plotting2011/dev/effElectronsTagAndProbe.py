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
gROOT.SetBatch(kFALSE);
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

File1 = GetFile("$LQDATA/TagAndProbe_1/TagAndProbeTest.root")

#--- Define all the histograms

#--------- TagAndProbe ----------
outputFile = "effElectrons.ps"
leptonType = "electron"
## Eta histograms
h_eta_recoEleMatchProbe_PassEleOffline = GetHisto( "eta_recoEleMatchProbe_PassEleOffline" , File1)
h_eta_recoEleMatchProbe_PassEleOfflineAndWP80 = GetHisto( "eta_recoEleMatchProbe_PassEleOfflineAndWP80" , File1)
h_eta_recoEleMatchProbe_PassEleOfflineAndTag = GetHisto( "eta_recoEleMatchProbe_PassEleOfflineAndTag" , File1)

## Pt histograms
h_pt_recoEleMatchProbe_PassEleOffline = GetHisto( "pt_recoEleMatchProbe_PassEleOffline" , File1)
h_pt_recoEleMatchProbe_PassEleOfflineAndWP80 = GetHisto( "pt_recoEleMatchProbe_PassEleOfflineAndWP80" , File1)
h_pt_recoEleMatchProbe_PassEleOfflineAndTag = GetHisto( "pt_recoEleMatchProbe_PassEleOfflineAndTag" , File1)

## nVtx histograms
h_NPV_recoEleMatchProbe_PassEleOffline = GetHisto( "NPV_recoEleMatchProbe_PassEleOffline" , File1)
h_NPV_recoEleMatchProbe_PassEleOfflineAndWP80 = GetHisto( "NPV_recoEleMatchProbe_PassEleOfflineAndWP80" , File1)
h_NPV_recoEleMatchProbe_PassEleOfflineAndTag = GetHisto( "NPV_recoEleMatchProbe_PassEleOfflineAndTag" , File1)

h_NPV_recoEleMatchProbe_PassEleOffline_bar = GetHisto( "NPV_recoEleMatchProbe_PassEleOffline_bar" , File1)
h_NPV_recoEleMatchProbe_PassEleOfflineAndWP80_bar = GetHisto( "NPV_recoEleMatchProbe_PassEleOfflineAndWP80_bar" , File1)
h_NPV_recoEleMatchProbe_PassEleOfflineAndTag_bar = GetHisto( "NPV_recoEleMatchProbe_PassEleOfflineAndTag_bar" , File1)

h_NPV_recoEleMatchProbe_PassEleOffline_end = GetHisto( "NPV_recoEleMatchProbe_PassEleOffline_end" , File1)
h_NPV_recoEleMatchProbe_PassEleOfflineAndWP80_end = GetHisto( "NPV_recoEleMatchProbe_PassEleOfflineAndWP80_end" , File1)
h_NPV_recoEleMatchProbe_PassEleOfflineAndTag_end = GetHisto( "NPV_recoEleMatchProbe_PassEleOfflineAndTag_end" , File1)

############################
#variable binning
MyBins = [0,10,20,30,40,50,60,70,80,90,100,150,200]
############################


#--- WP80 Eff
eff_lepton_eta_wp80 = GetEffFixBinning( h_eta_recoEleMatchProbe_PassEleOfflineAndWP80 , h_eta_recoEleMatchProbe_PassEleOffline
                                        , 0.9 , 20 , 1
                                        , leptonType+" \\eta_{reco}" , "trigger eff WP80 wrt to good reco electron"
                                        , -2.5 , 2.5, 0, 1)

eff_lepton_pt_wp80 = GetEffVarBinning( h_pt_recoEleMatchProbe_PassEleOfflineAndWP80 , h_pt_recoEleMatchProbe_PassEleOffline
                                       , 0.9 , 20 , 1
                                       , leptonType+" pT_{reco} (GeV)" , "trigger eff WP80 wrt to good reco electron"
                                       , 0 , 300
                                       , MyBins, 0.5, 1)

eff_lepton_NPV_wp80 = GetEffFixBinning( h_NPV_recoEleMatchProbe_PassEleOfflineAndWP80 , h_NPV_recoEleMatchProbe_PassEleOffline
                                        , 0.9 , 20 , 1
                                        , "NPV" , "trigger eff WP80 wrt to good reco electron"
                                        , 0 , 30 , 0. , 1)

eff_lepton_NPV_wp80_bar = GetEffFixBinning( h_NPV_recoEleMatchProbe_PassEleOfflineAndWP80_bar , h_NPV_recoEleMatchProbe_PassEleOffline_bar
                                            , 0.9 , 20 , 1
                                            , "NPV" , "trigger eff WP80 wrt to good reco electron (barrel)"
                                            , 0 , 30 , 0. , 1)

eff_lepton_NPV_wp80_end = GetEffFixBinning( h_NPV_recoEleMatchProbe_PassEleOfflineAndWP80_end , h_NPV_recoEleMatchProbe_PassEleOffline_end
                                            , 0.9 , 20 , 1
                                            , "NPV" , "trigger eff WP80 wrt to good reco electron (endcap)"
                                            , 0 , 30 , 0. , 1)

#--- TAG Eff
eff_lepton_eta_tag = GetEffFixBinning( h_eta_recoEleMatchProbe_PassEleOfflineAndTag , h_eta_recoEleMatchProbe_PassEleOffline
                                        , 0.9 , 20 , 1
                                        , leptonType+" \\eta_{reco}" , "trigger eff TAG wrt to good reco electron"
                                        , -2.5 , 2.5, 0, 1)

eff_lepton_pt_tag = GetEffVarBinning( h_pt_recoEleMatchProbe_PassEleOfflineAndTag , h_pt_recoEleMatchProbe_PassEleOffline
                                       , 0.9 , 20 , 1
                                       , leptonType+" pT_{reco} (GeV)" , "trigger eff TAG wrt to good reco electron"
                                       , 0 , 300
                                       , MyBins, 0.5, 1)

eff_lepton_NPV_tag = GetEffFixBinning( h_NPV_recoEleMatchProbe_PassEleOfflineAndTag , h_NPV_recoEleMatchProbe_PassEleOffline
                                       , 0.9 , 20 , 1
                                       , "NPV" , "trigger eff TAG wrt to good reco electron"
                                       , 0 , 30 , 0. , 1)

eff_lepton_NPV_tag_bar = GetEffFixBinning( h_NPV_recoEleMatchProbe_PassEleOfflineAndTag_bar , h_NPV_recoEleMatchProbe_PassEleOffline_bar
                                           , 0.9 , 20 , 1
                                           , "NPV" , "trigger eff TAG wrt to good reco electron (barrel)"
                                           , 0 , 30 , 0. , 1)

eff_lepton_NPV_tag_end = GetEffFixBinning( h_NPV_recoEleMatchProbe_PassEleOfflineAndTag_end , h_NPV_recoEleMatchProbe_PassEleOffline_end
                                           , 0.9 , 20 , 1
                                           , "NPV" , "trigger eff TAG wrt to good reco electron (endcap)"
                                           , 0 , 30 , 0. , 1)


#--- Calculate total acceptance x efficiency

#eff_times_acc_lepton_nVtx_lepton_gen_RECO_ID_ISO = GetEffFixBinning( h_lepton_nVtx_gen_RECO_ID_ISO , h_lepton_nVtx_gen_all
#                                                                     , 0.9 , 20 , 4
#                                                                     , "Number of reco vertices (in events with "+leptonType+")" , "acceptance x efficiency"
#                                                                     , 0 , 30 , 0., 1)

#--- Final plots

#------------------------------------------------------

#--- WP80 Eff

## Acceptance vs eta - lepton
c1 = TCanvas()
c1.SetGridy();
c1.SetGridx();
eff_lepton_eta_wp80.Draw("ap")

c1.Print(outputFile+"(")

## Acceptance vs pt - lepton
c2 = TCanvas()
c2.SetGridy();
c2.SetGridx();
eff_lepton_pt_wp80.Draw("ap")

c2.Print(outputFile)

## Acceptance vs NPV - lepton
c3 = TCanvas()
c3.SetGridy();
c3.SetGridx();
eff_lepton_NPV_wp80.Draw("ap")

c3.Print(outputFile)

## Acceptance vs NPV (barrel) - lepton
c4 = TCanvas()
c4.SetGridy();
c4.SetGridx();
eff_lepton_NPV_wp80_bar.Draw("ap")

c4.Print(outputFile)

## Acceptance vs NPV (endcap) - lepton
c5 = TCanvas()
c5.SetGridy();
c5.SetGridx();
eff_lepton_NPV_wp80_end.Draw("ap")

c5.Print(outputFile)

#--- TAG Eff

## Acceptance vs eta - lepton
c11 = TCanvas()
c11.SetGridy();
c11.SetGridx();
eff_lepton_eta_tag.Draw("ap")

c11.Print(outputFile)

## Acceptance vs pt - lepton
c21 = TCanvas()
c21.SetGridy();
c21.SetGridx();
eff_lepton_pt_tag.Draw("ap")

c21.Print(outputFile)

## Acceptance vs NPV - lepton
c31 = TCanvas()
c31.SetGridy();
c31.SetGridx();
eff_lepton_NPV_tag.Draw("ap")

c31.Print(outputFile)

## Acceptance vs NPV (barrel) - lepton
c41 = TCanvas()
c41.SetGridy();
c41.SetGridx();
eff_lepton_NPV_tag_bar.Draw("ap")

c41.Print(outputFile)

## Acceptance vs NPV (endcap) - lepton
c51 = TCanvas()
c51.SetGridy();
c51.SetGridx();
eff_lepton_NPV_tag_end.Draw("ap")

c51.Print(outputFile+")")

#---Create legend
#globals()['legend3'] = TLegend(0.348609,0.184915,0.653028,0.377129)
#legend3.SetFillColor(kWhite)
#legend3.SetMargin(0.3)
#legend3.AddEntry(eff_lepton_eta_gen_RECO,"only reco","p")
#legend3.AddEntry(eff_lepton_eta_gen_RECO_ID_ISO,"reco+ID+ISO","p")
#legend3.Draw()
#c3.Update()
#gPad.RedrawAxis()
#gPad.Modified()
#c3.Print(outputFile)

#------------------------------------------------------

h_NPV_recoEleMatchProbe_PassEleOffline = GetHisto( "NPV_recoEleMatchProbe_PassEleOffline" , File1)
h_NPV_recoEleMatchProbe_PassEleOfflineAndWP80 = GetHisto( "NPV_recoEleMatchProbe_PassEleOfflineAndWP80" , File1)
h_NPV_recoEleMatchProbe_PassEleOfflineAndTag = GetHisto( "NPV_recoEleMatchProbe_PassEleOfflineAndTag" , File1)

h_NPV_recoEleMatchProbe_PassEleOffline_bar = GetHisto( "NPV_recoEleMatchProbe_PassEleOffline_bar" , File1)
h_NPV_recoEleMatchProbe_PassEleOfflineAndWP80_bar = GetHisto( "NPV_recoEleMatchProbe_PassEleOfflineAndWP80_bar" , File1)
h_NPV_recoEleMatchProbe_PassEleOfflineAndTag_bar = GetHisto( "NPV_recoEleMatchProbe_PassEleOfflineAndTag_bar" , File1)

h_NPV_recoEleMatchProbe_PassEleOffline_end = GetHisto( "NPV_recoEleMatchProbe_PassEleOffline_end" , File1)
h_NPV_recoEleMatchProbe_PassEleOfflineAndWP80_end = GetHisto( "NPV_recoEleMatchProbe_PassEleOfflineAndWP80_end" , File1)
h_NPV_recoEleMatchProbe_PassEleOfflineAndTag_end = GetHisto( "NPV_recoEleMatchProbe_PassEleOfflineAndTag_end" , File1)


#Print out
print "Eff WP80 (Bar+End) = " + str(h_NPV_recoEleMatchProbe_PassEleOfflineAndWP80.Integral()/h_NPV_recoEleMatchProbe_PassEleOffline.Integral())
print "Eff WP80 (Bar) = " + str(h_NPV_recoEleMatchProbe_PassEleOfflineAndWP80_bar.Integral()/h_NPV_recoEleMatchProbe_PassEleOffline_bar.Integral())
print "Eff WP80 (End) = " + str(h_NPV_recoEleMatchProbe_PassEleOfflineAndWP80_end.Integral()/h_NPV_recoEleMatchProbe_PassEleOffline_end.Integral())

print "Eff TAG (Bar+End) = " + str(h_NPV_recoEleMatchProbe_PassEleOfflineAndTag.Integral()/h_NPV_recoEleMatchProbe_PassEleOffline.Integral())
print "Eff TAG (Bar) = " + str(h_NPV_recoEleMatchProbe_PassEleOfflineAndTag_bar.Integral()/h_NPV_recoEleMatchProbe_PassEleOffline_bar.Integral())
print "Eff TAG (End) = " + str(h_NPV_recoEleMatchProbe_PassEleOfflineAndTag_end.Integral()/h_NPV_recoEleMatchProbe_PassEleOffline_end.Integral())

LintWP80 = 0.75 #fb-1
LintTAG = 3.9   #fb-1

eff_2011_default = ( h_NPV_recoEleMatchProbe_PassEleOfflineAndWP80.Integral()/h_NPV_recoEleMatchProbe_PassEleOffline.Integral() * LintWP80 + h_NPV_recoEleMatchProbe_PassEleOfflineAndTag.Integral()/h_NPV_recoEleMatchProbe_PassEleOffline.Integral() * LintTAG ) / (LintWP80 + LintTAG)
eff_2011_bar = ( h_NPV_recoEleMatchProbe_PassEleOfflineAndWP80_bar.Integral()/h_NPV_recoEleMatchProbe_PassEleOffline_bar.Integral() * LintWP80 + h_NPV_recoEleMatchProbe_PassEleOfflineAndTag_bar.Integral()/h_NPV_recoEleMatchProbe_PassEleOffline_bar.Integral() * LintTAG ) / (LintWP80 + LintTAG)
eff_2011_end = ( h_NPV_recoEleMatchProbe_PassEleOfflineAndWP80_end.Integral()/h_NPV_recoEleMatchProbe_PassEleOffline_end.Integral() * LintWP80 + h_NPV_recoEleMatchProbe_PassEleOfflineAndTag_end.Integral()/h_NPV_recoEleMatchProbe_PassEleOffline_end.Integral() * LintTAG ) / (LintWP80 + LintTAG)


print "Eff 2011 = " + str(eff_2011_default) + " + "  + str(eff_2011_bar - eff_2011_default) + " - " + str(eff_2011_default - eff_2011_end)


## Terminate the program
print "Press ENTER to terminate"
wait=raw_input()
