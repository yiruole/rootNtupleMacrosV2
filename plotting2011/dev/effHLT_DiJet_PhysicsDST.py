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

File1 = GetFile("lowHTSkim.root")

#--------- Muons ----------
outputFile = "effHLT.ps"
## Eta histograms
h_M_FatPFJet1FatPFJet2_HT150 = GetHisto( "M_FatPFJet1FatPFJet2_HT150" , File1)
h_M_FatPFJet1FatPFJet2_HT150_M400 = GetHisto( "M_FatPFJet1FatPFJet2_HT150_M400" , File1)
h_M_FatPFJet1FatPFJet2_HT150_M400ORHT350 = GetHisto( "M_FatPFJet1FatPFJet2_HT150_M400ORHT350" , File1)

h_M_FatPFJet1FatPFJet2_HT200 = GetHisto( "M_FatPFJet1FatPFJet2_HT200" , File1)
h_M_FatPFJet1FatPFJet2_HT200_M400 = GetHisto( "M_FatPFJet1FatPFJet2_HT200_M400" , File1)
h_M_FatPFJet1FatPFJet2_HT200_M400ORHT350 = GetHisto( "M_FatPFJet1FatPFJet2_HT200_M400ORHT350" , File1)

h_M_FatPFJet1FatPFJet2_HT250 = GetHisto( "M_FatPFJet1FatPFJet2_HT250" , File1)
h_M_FatPFJet1FatPFJet2_HT250_M400 = GetHisto( "M_FatPFJet1FatPFJet2_HT250_M400" , File1)
h_M_FatPFJet1FatPFJet2_HT250_M400ORHT350 = GetHisto( "M_FatPFJet1FatPFJet2_HT250_M400ORHT350" , File1)

############################
#variable binning
#MyBins = [200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840,860,880,900,920,940,960,980,1000]
#MyBins = [200,240,280,320,360,400,440,480,520,560,600,640,680,720,760,800,840,880,920,960,1000]
MyBins = [200,240,280,320,360,400,440,480,520,560,600,640,680,1000]
############################

#--- Calculate HLT Efficiency wrt HT150
eff_HLTM400_wrt_HT150 = GetEffVarBinning( h_M_FatPFJet1FatPFJet2_HT150_M400 , h_M_FatPFJet1FatPFJet2_HT150
                                          , 0.9 , 20 , 1
                                          , "FatPFDiJetMass@HLT [GeV]" , "HLT Efficiency wrt HT150"
                                          , 0, 1000, MyBins, 0, 1.1)

eff_HLTM400ORHT350_wrt_HT150 = GetEffVarBinning( h_M_FatPFJet1FatPFJet2_HT150_M400ORHT350 , h_M_FatPFJet1FatPFJet2_HT150
                                                 , 0.9 , 20 , 2
                                                 , "FatPFDiJetMass@HLT [GeV]" , "HLT Efficiency wrt HT150"
                                                 , 0, 1000, MyBins, 0, 1.1)

eff_HLTM400_wrt_HT200 = GetEffVarBinning( h_M_FatPFJet1FatPFJet2_HT200_M400 , h_M_FatPFJet1FatPFJet2_HT200
                                          , 0.9 , 20 , 1
                                          , "FatPFDiJetMass@HLT [GeV]" , "HLT Efficiency wrt HT200"
                                          , 0, 1000, MyBins, 0, 1.1)

eff_HLTM400ORHT350_wrt_HT200 = GetEffVarBinning( h_M_FatPFJet1FatPFJet2_HT200_M400ORHT350 , h_M_FatPFJet1FatPFJet2_HT200
                                                 , 0.9 , 20 , 2
                                                 , "FatPFDiJetMass@HLT [GeV]" , "HLT Efficiency wrt HT200"
                                                 , 0, 1000, MyBins, 0, 1.1)

eff_HLTM400_wrt_HT250 = GetEffVarBinning( h_M_FatPFJet1FatPFJet2_HT250_M400 , h_M_FatPFJet1FatPFJet2_HT250
                                          , 0.9 , 20 , 1
                                          , "FatPFDiJetMass@HLT [GeV]" , "HLT Efficiency wrt HT250"
                                          , 0, 1000, MyBins, 0, 1.1)

eff_HLTM400ORHT350_wrt_HT250 = GetEffVarBinning( h_M_FatPFJet1FatPFJet2_HT250_M400ORHT350 , h_M_FatPFJet1FatPFJet2_HT250
                                                 , 0.9 , 20 , 2
                                                 , "FatPFDiJetMass@HLT [GeV]" , "HLT Efficiency wrt HT250"
                                                 , 0, 1000, MyBins, 0, 1.1)

#--- Final plots

#------------------------------------------------------


#------------------------------------------------------

## Eff wrt HT150
c1 = TCanvas()
c1.SetGridy();
c1.SetGridx();
eff_HLTM400_wrt_HT150.Draw("ap")
eff_HLTM400ORHT350_wrt_HT150.Draw("psame")

#---Create legend
globals()['legend1'] = TLegend(0.534483,0.184322,0.840517,0.375)
legend1.SetFillColor(kWhite)
legend1.SetMargin(0.2)
legend1.AddEntry(eff_HLTM400_wrt_HT150,"FatJetMass400","p")
legend1.AddEntry(eff_HLTM400ORHT350_wrt_HT150,"FatJetMass400 OR HT350","p")
legend1.Draw()
c1.Update()
gPad.RedrawAxis()
gPad.Modified()
c1.Print(outputFile+"(")

## Eff wrt HT200
c2 = TCanvas()
c2.SetGridy();
c2.SetGridx();
eff_HLTM400_wrt_HT200.Draw("ap")
eff_HLTM400ORHT350_wrt_HT200.Draw("psame")

#---Create legend
globals()['legend2'] = TLegend(0.534483,0.184322,0.840517,0.375) 
legend2.SetFillColor(kWhite)
legend2.SetMargin(0.2)
legend2.AddEntry(eff_HLTM400_wrt_HT200,"FatJetMass400","p")
legend2.AddEntry(eff_HLTM400ORHT350_wrt_HT200,"FatJetMass400 OR HT350","p")
legend2.Draw()
c2.Update()
gPad.RedrawAxis()
gPad.Modified()
c2.Print(outputFile)

## Eff wrt HT250
c3 = TCanvas()
c3.SetGridy();
c3.SetGridx();
eff_HLTM400_wrt_HT250.Draw("ap")
eff_HLTM400ORHT350_wrt_HT250.Draw("psame")

#---Create legend
globals()['legend3'] = TLegend(0.534483,0.184322,0.840517,0.375) 
legend3.SetFillColor(kWhite)
legend3.SetMargin(0.2)
legend3.AddEntry(eff_HLTM400_wrt_HT250,"FatJetMass400","p")
legend3.AddEntry(eff_HLTM400ORHT350_wrt_HT250,"FatJetMass400 OR HT350","p")
legend3.Draw()
c3.Update()
gPad.RedrawAxis()
gPad.Modified()
c3.Print(outputFile+")")

#------------------------------------------------------

## Terminate the program
print "Press ENTER to terminate"
wait=raw_input()
