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
gStyle.SetOptStat(0)
gStyle.SetPalette(1)
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
#gStyle.SetPadTickX(1);
#gStyle.SetPadTickY(1);
#--- TODO: WHY IT DOES NOT LOAD THE DEFAULT ROOTLOGON.C ? ---#

def GetFile(filename):
    file = TFile(filename)
    if( not file):
        print "ERROR: file " + filename + " not found"
        print "exiting..."
        sys.exit()
    return file


def GetHisto( histoName , file ):
    histo = file.Get( histoName )
    if( not histo):
        print "ERROR: histo " + histoName + " not found in " + file.GetName()
        print "exiting..."
        sys.exit()
    return histo



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
        
############# DON'T NEED TO MODIFY ANYTHING HERE - END #######################
##############################################################################


##############################################################################
############# USER CODE - BEGING #############################################

#--- Input root file

#File1 = GetFile("/home/prumerio/cms/phys/lq/2010/heepEffic/rootNtupleAnalyzerV2/data/output/MLQ300.root")      #MLQ300
File1 = GetFile("/home/prumerio/cms/phys/lq/2010/heepEffic/rootNtupleAnalyzerV2/data/output/relValZee.root")  #RelValZee

                
#--- Get and define all the histograms

# Number of electrons
h_N_ele_Gen = GetHisto("N_ele_Gen", File1)
h_N_ele_Gen_etaPtCut = GetHisto("N_ele_Gen_etaPtCut", File1)
h_N_ele_Gen_matched = GetHisto("N_ele_Gen_matched", File1)
#h_N_ele_Pt_Gen_matched_ID = GetHisto("N_ele_Pt_Gen_matched_ID", File1)
h_N_ele_Pt_Gen_matched_ID_ISO = GetHisto("N_ele_Pt_Gen_matched_ID_ISO", File1)
          
#-- reconstruction efficiencies

## Eta histograms
h_ele_Eta_Gen = GetHisto( "ele_Eta_Gen" , File1)
h_ele_Eta_Gen_etaPtCut = GetHisto( "ele_Eta_Gen_etaPtCut" , File1)
h_ele_Eta_Gen_matched = GetHisto( "ele_Eta_Gen_matched" , File1)
###h_ele_Eta_Gen_matched_ID = GetHisto( "ele_Eta_Gen_matched_ID" , File1)
h_ele_Eta_Gen_matched_ID_ISO = GetHisto( "ele_Eta_Gen_matched_ID_ISO" , File1)

## Pt histograms
h_ele_Pt_Gen = GetHisto( "ele_Pt_Gen" , File1)
h_ele_Pt_Gen_etaPtCut = GetHisto( "ele_Pt_Gen_etaPtCut" , File1)
h_ele_Pt_Gen_matched = GetHisto( "ele_Pt_Gen_matched" , File1)
###h_ele_Pt_Gen_matched_ID = GetHisto( "ele_Pt_Gen_matched_ID" , File1)
h_ele_Pt_Gen_matched_ID_ISO = GetHisto( "ele_Pt_Gen_matched_ID_ISO" , File1)


############################
#variable binning
#MyBins = [25, 100, 150, 200, 250, 300, 350, 400, 500, 600, 800, 1200] # MLQ300
MyBins = [20, 30, 40, 50, 60, 80, 100, 150, 400]                   # RelValZee
############################


#--- Calculate Acceptance
accept_ele_eta_gen = GetEffFixBinning( h_ele_Eta_Gen_etaPtCut , h_ele_Eta_Gen
                                       , 0.9 , 20 , 1
                                       , "electron \\eta_{gen}" , "acceptance"
                                       , -2.5 , 2.5, 0, 1)

# accept_ele_pt_gen = GetEffVarBinning( h_ele_Pt_Gen_etaPtCut , h_ele_Pt_Gen
#                                       , 0.9 , 20 , 1
#                                       , "electron pT_{gen} (GeV)" , "acceptance"
#                                       , 0 , 1200
#                                       , MyBins, 0 , 1)
accept_ele_pt_gen = GetEffFixBinning( h_ele_Pt_Gen_etaPtCut , h_ele_Pt_Gen
                                      , 0.9 , 20 , 1
                                      , "electron pT_{gen} (GeV)" , "acceptance"
                                      , 0 , 1200, 0 , 1)

#--- Calculate Absolute Efficiency
eff_ele_eta_gen_RECO = GetEffFixBinning( h_ele_Eta_Gen_matched , h_ele_Eta_Gen_etaPtCut
                                           , 0.9 , 20 , 1
                                           , "electron \\eta_{gen}" , "relative efficiency wrt to acceptance"
                                           , -2.5 , 2.5, 0, 1)

# eff_ele_eta_gen_RECO_ID = GetEffFixBinning( h_ele_Eta_Gen_matched_ID , h_ele_Eta_Gen_etaPtCut
#                                             , 0.9 , 20 , 2
#                                             , "electron \\eta_{gen}" , "relative efficiency wrt to acceptance"
#                                             , -2.5 , 2.5, 0, 1)

eff_ele_eta_gen_RECO_ID_ISO = GetEffFixBinning( h_ele_Eta_Gen_matched_ID_ISO , h_ele_Eta_Gen_etaPtCut
                                                , 0.9 , 20 , 4
                                                , "electron \\eta_{gen}" , "relative efficiency wrt to acceptance"
                                                , -2.5 , 2.5, 0, 1)


# eff_ele_pt_gen_RECO = GetEffVarBinning( h_ele_Pt_Gen_matched , h_ele_Pt_Gen_etaPtCut
#                                         , 0.9 , 20 , 1
#                                         , "electron pT_{gen} (GeV)" , "relative efficiency wrt to acceptance"
#                                         , 0 , 1200
#                                         , MyBins, 0 , 1)
eff_ele_pt_gen_RECO = GetEffFixBinning( h_ele_Pt_Gen_matched , h_ele_Pt_Gen_etaPtCut
                                        , 0.9 , 20 , 1
                                        , "electron pT_{gen} (GeV)" , "relative efficiency wrt to acceptance"
                                        , 0 , 1200, 0 , 1)

# eff_ele_pt_gen_RECO_ID = GetEffVarBinning( h_ele_Pt_Gen_matched_ID , h_ele_Pt_Gen_etaPtCut
#                                            , 0.9 , 20 , 2
#                                            , "electron pT_{gen} (GeV)" , "relative efficiency wrt to acceptance"
#                                            , 0 , 1200
#                                            , MyBins, 0.5, 1)

# eff_ele_pt_gen_RECO_ID_ISO = GetEffVarBinning( h_ele_Pt_Gen_matched_ID_ISO , h_ele_Pt_Gen_etaPtCut
#                                                , 0.9 , 20 , 4
#                                                , "electron pT_{gen} (GeV)" , "relative efficiency wrt to acceptance"
#                                                , 0 , 1200
#                                                , MyBins, 0 , 1)
eff_ele_pt_gen_RECO_ID_ISO = GetEffFixBinning( h_ele_Pt_Gen_matched_ID_ISO , h_ele_Pt_Gen_etaPtCut
                                               , 0.9 , 20 , 4
                                               , "electron pT_{gen} (GeV)" , "relative efficiency wrt to acceptance"
                                               , 0 , 1200, 0 , 1)


#--- Calculate Relative Efficiency


# eff_ele_pt_gen_RECO_ID_rel = GetEffVarBinning( h_ele_Pt_Gen_matched_ID , h_ele_Pt_Gen_matched
#                                                , 0.9 , 20 , 2
#                                                , "electron pT_{gen} (GeV)" , "relative efficiency"
#                                                , 0 , 1200
#                                                , MyBins, 0.5, 1)

# eff_ele_pt_gen_RECO_ID_ISO_rel = GetEffVarBinning( h_ele_Pt_Gen_matched_ID_ISO , h_ele_Pt_Gen_matched_ID
#                                                , 0.9 , 20 , 4
#                                                , "electron pT_{gen} (GeV)" , "relative efficiency"
#                                                , 0 , 1200
#                                                , MyBins, 0.5, 1)

eff_ele_eta_gen_HEEP_RECO = GetEffFixBinning( h_ele_Eta_Gen_matched_ID_ISO , h_ele_Eta_Gen_matched
                                            , 0.9 , 20 , 2
                                            , "electron \\eta_{gen}" , "relative efficiency"
                                            , -2.5 , 2.5, 0, 1)

# eff_ele_pt_gen_HEEP_RECO = GetEffVarBinning( h_ele_Pt_Gen_matched_ID_ISO , h_ele_Pt_Gen_matched
#                                             , 0.9 , 20 , 2
#                                             , "electron pT_{gen} (GeV)" , "relative efficiency"
#                                             , 0 , 1200
#                                             , MyBins, 0.5, 1)
eff_ele_pt_gen_HEEP_RECO = GetEffFixBinning( h_ele_Pt_Gen_matched_ID_ISO , h_ele_Pt_Gen_matched
                                            , 0.9 , 20 , 2
                                            , "electron pT_{gen} (GeV)" , "relative efficiency"
                                            , 0 , 1200, 0, 1)


#--- Final plots

########################################################################
## The Plot class: add members if needed
class Plot:
    histos  = [] # list of hisotgrams to be plotted in this plot
    keys    = [] # list of keys to be put in the legend (1 key per histo)
    htit    = "" # histo title
    xtit    = "" # xtitle
    ytit    = "" # ytitle 
    xmin    = "" # min x axis range (need to set both min and max. Leave it as is for full range)
    xmax    = "" # max x axis range (need to set both min and max. Leave it as is for full range)
    ymin    = "" # min y axis range (need to set both min and max. Leave it as is for full range)
    ymax    = "" # max y axis range (need to set both min and max. Leave it as is for full range)
    lpos    = "" # legend position (default = top-right, options="bottom-center","bottom-right")
    lxpos   = 0.60 # legend x size (do not set lpos if want to use lxpos)
    lypos   = 0.70 # legend y size (do not set lpos if want to use lypos)
    lxsiz   = 0.29 # legend x size
    lysiz   = 0.20 # legend y size
    def Draw(self, fileps):
        ih=0 # index of histo within a plot
        canvas = TCanvas()
        canvas.SetGridx()
        canvas.SetGridy()
        if (plot.lpos=="bottom-center"):
            plot.lxpos=0.35
            plot.lypos=0.15
        if (plot.lpos=="bottom-right"):
            plot.lxpos=0.60
            plot.lypos=0.15
        legend = TLegend(plot.lxpos, plot.lypos, plot.lxpos+plot.lxsiz, plot.lypos+plot.lysiz)
        legend.SetFillColor(kWhite)
        legend.SetMargin(0.2)
        for histo in plot.histos:
            histo.SetMarkerStyle(20+2*ih)
            histo.SetMarkerColor(1+ih)
            legend.AddEntry(histo, plot.keys[ih],"p")
            legend.SetTextSize(0.03)
            legend.SetMargin(0.1)
            if ih==0:
                histo.SetTitle(plot.htit)
                histo.GetXaxis().SetTitle(plot.xtit)
                histo.GetYaxis().SetTitle(plot.ytit)
                if (plot.xmin!="" and plot.xmax!=""):
                    histo.GetXaxis().SetRangeUser(plot.xmin,plot.xmax)
                if (plot.ymin!="" and plot.ymax!=""):
                    histo.GetYaxis().SetRangeUser(plot.ymin,plot.ymax)
                histo.Draw("ap")
            else:
                histo.Draw("psame")
            ih=ih+1
        legend.Draw()
        canvas.Update()
        gPad.RedrawAxis()
        gPad.Modified()
        canvas.Print(fileps)
#########################################################################

## Final plots are defined here.
## Simply add or remove plots and update the list plots = [plot0, plot1, ...] below
    
## Eta distributions - electron
plot0 = Plot() 
plot0.htit   = "N. of electrons vs. \\eta_{gen}"
plot0.histos = [h_ele_Eta_Gen, h_ele_Eta_Gen_etaPtCut, h_ele_Eta_Gen_matched, h_ele_Eta_Gen_matched_ID_ISO]
plot0.keys   = ["No cuts"    , "Acceptance cuts"     ,  "RECO"              , "HEEP"]
plot0.xtit   = "electron \\eta_{gen}"
plot0.xmin   = -6
plot0.xmax   = 6

## pt distributions - electron
plot1 = Plot()
plot1.htit   = "N. of electrons vs. pt_{gen}"
plot1.histos = [h_ele_Pt_Gen, h_ele_Pt_Gen_etaPtCut, h_ele_Pt_Gen_matched, h_ele_Pt_Gen_matched_ID_ISO]
plot1.keys   = ["No cuts"   , "Acceptance cuts"    , "RECO"              , "HEEP"]
plot1.xtit   = "electron pt_{gen} (GeV)"
plot1.xmin   = 0
plot1.xmax   = 150 # RelValZee 
#plot1.xmax   = 600 # MLQ300

## Acceptance vs eta - electron
plot2 = Plot()
plot2.htit   = "Acceptance vs. \\eta_{gen}"
plot2.histos = [accept_ele_eta_gen]
plot2.keys   = ["Acceptance cuts: Pt>25 GeV AND (|\\eta|<1.442 OR 1.560<|\\eta|<2.5)"]
plot2.xtit   = "electron \\eta_{gen}"
plot2.ytit   = "Acceptance"
plot2.ymin   = 0
plot2.ymax   = 1.05
plot2.lxpos  = 0.20
plot2.lypos  = 0.15
plot2.lxsiz  = 0.65
plot2.lysiz  = 0.10

## Acceptance vs pt - electron
plot3 = Plot()
plot3.htit   = "Acceptance vs. pt_{gen}"
plot3.histos = [accept_ele_pt_gen]
plot3.keys   = ["Acceptance cuts: Pt>25 GeV AND (|\\eta|<1.442 OR 1.560<|\\eta|<2.5)"]
plot3.xtit   = "electron pt_{gen} (GeV)"
plot3.ytit   = "Acceptance"
plot3.ymin   = 0
plot3.ymax   = 1.05
plot3.lxpos  = 0.20
plot3.lypos  = 0.15
plot3.lxsiz  = 0.65
plot3.lysiz  = 0.10

## Eff vs eta - electron
plot4 = Plot()
plot4.htit   = "Efficiency vs. \\eta_{gen} within acceptance"
plot4.histos = [eff_ele_eta_gen_RECO, eff_ele_eta_gen_RECO_ID_ISO] # eff_ele_eta_gen_RECO_ID, 
plot4.keys   = ["RECO wrt Acceptance", "HEEP wrt Acceptance"] # "reco+ID"
plot4.xtit   = "electron \\eta_{gen}"
plot4.ytit   = "Efficiency"
plot4.ymin   = 0.5
plot4.ymax   = 1.01
plot4.lpos = "bottom-center"

## Eff vs pt - electron
plot5 = Plot()
plot5.htit   = "Efficiency vs. pt_{gen} within acceptance"
plot5.histos = [eff_ele_pt_gen_RECO, eff_ele_pt_gen_RECO_ID_ISO] # eff_ele_pt_gen_RECO_ID
plot5.keys   = ["RECO wrt Acceptance", "HEEP wrt Acceptance"] # "reco+ID"
plot5.xtit   = "electron pt_{gen} (GeV)"
plot5.ytit   = "Efficiency"
plot5.ymin   = 0.5
plot5.ymax   = 1.01
plot5.lpos = "bottom-right"

# ## Relative Eff vs eta - electron
# plot6 = Plot()
# plot6.histos = [eff_ele_eta_gen_RECO_ID_ISO_rel] # eff_ele_eta_gen_RECO_ID_rel, 
# plot6.keys   = ["HEEP eff wrt to RECO"] # "ID eff wrt to RECO",  
# plot6.xtit   = "electron \\eta_{gen}"
# plot6.ytit   = "Relative Efficiency"
# plot6.xmin   = ""
# plot6.xmax   = ""
# plot6.ymin   = 0
# plot6.ymax   = 1.05

# ## Relative Eff vs pt - electron
# plot7 = Plot()
# plot7.histos = [eff_ele_pt_gen_RECO_ID_ISO_rel] # eff_ele_pt_gen_RECO_ID_rel, 
# plot7.keys   = ["HEEP eff wrt to RECO"] # "ID eff wrt to RECO",  
# plot7.xtit   = "electron pt_{gen} (GeV)"
# plot7.ytit   = "Relative Efficiency"
# plot7.xmin   = ""
# plot7.xmax   = ""
# plot7.ymin   = 0
# plot7.ymax   = 1.05

## Relative Eff vs eta - electron
plot8 = Plot()
plot8.htit   = "Relative Efficiency vs. \\eta_{gen} (within acceptance)"
plot8.histos = [eff_ele_eta_gen_HEEP_RECO]  
plot8.keys   = ["HEEP eff wrt to RECO"]  
plot8.xtit   = "electron \\eta_{gen}"
plot8.ytit   = "Relative Efficiency"
plot8.ymin   = 0.8
plot8.ymax   = 1.0
plot8.lpos   = "bottom-center"

## Relative Eff vs pt - electron
plot9 = Plot()
plot9.htit   = "Relative Efficiency vs. pt_{gen} (within acceptance)"
plot9.histos = [eff_ele_pt_gen_HEEP_RECO]  
plot9.keys   = ["HEEP eff wrt to RECO"]  
plot9.xtit   = "electron pt_{gen} (GeV)"
plot9.ytit   = "Relative Efficiency"
plot9.ymin   = 0.8
plot9.ymax   = 1.0
plot9.lpos = "bottom-right"

## Number of electrons
plot10 = Plot()
plot10.histos = [h_N_ele_Gen, h_N_ele_Gen_etaPtCut, h_N_ele_Gen_matched, h_N_ele_Pt_Gen_matched_ID_ISO] # h_N_ele_Pt_Gen_matched_ID
plot10.keys   = ["","","",""] 
plot10.xtit   = "Number of electrons"


## List of plots to be plotted
plots = [plot0, plot1, plot2, plot3, plot4, plot5, plot8, plot9, plot10]

## Print out efficiencies
print "Acceptance = " + str(h_ele_Eta_Gen_etaPtCut.Integral()/h_ele_Eta_Gen.Integral())
print "Reco eff wrt to acceptance = " + str(h_ele_Eta_Gen_matched.Integral()/h_ele_Eta_Gen_etaPtCut.Integral())
#print "Reco+ID eff wrt to acceptance = " + str(h_ele_Eta_Gen_matched_ID.Integral()/h_ele_Eta_Gen_etaPtCut.Integral())
print "Reco+ID+ISO eff wrt to acceptance = " + str(h_ele_Eta_Gen_matched_ID_ISO.Integral()/h_ele_Eta_Gen_etaPtCut.Integral())
print "Reco+ID+ISO eff wrt generated = " + str(h_ele_Eta_Gen_matched_ID_ISO.Integral()/h_ele_Eta_Gen.Integral())
print "Reco+ID+ISO eff wrt RECO = " + str(h_ele_Eta_Gen_matched_ID_ISO.Integral()/h_ele_Pt_Gen_matched.Integral())


############# USER CODE - END ################################################
##############################################################################

## Generate and print the plots from the list 'plots' define above
#fileps="MLQ300.ps"
fileps="RelValZee.ps"
c = TCanvas()
c.Print(fileps+"[")
for plot in plots:
    plot.Draw(fileps)
c.Print(fileps+"]")


## Terminate the program
#print "Press ENTER to terminate"
#wait=raw_input()
