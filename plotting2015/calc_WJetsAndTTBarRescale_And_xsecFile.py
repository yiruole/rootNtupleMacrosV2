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
import math

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
	#print "bin: " +str(bin)
        if(bin==bmax and bmaxResidual==histo.GetBinContent(bmax)): # skip last bin if out of range
            print "     --> skip bin: " + str(bin)
        else:
            error = error + histo.GetBinError(bin)**2
            #print "error**2 : " + str(error)

    error = math.sqrt(error)
    print  " "
    return error


## The Plot class: add members if needed
class Plot:
    histoDATA    = "" # DATA
    histoTTbar   = "" # MCTTbar
    histoMCall   = "" # MCall
    histoQCD     = "" # QCD
    histoZJet    = ""
    histoWJet    = ""
    histoSingleTop = ""
    histoPhotonJets = ""
    histoDiboson = ""
    xtit         = "" # xtitle
    ytit         = "" # ytitle
    xmin         = "" # set xmin to calculate rescaling factor (-- please take into account the bin size of histograms --)
    xmax         = "" # # set xmax to calculate rescaling factor (-- please take into account the bin size of histograms --)
    xminplot     = "" # min x axis range (need to set both min and max. Leave it as is for full range)
    xmaxplot     = "" # max x axis range (need to set both min and max. Leave it as is for full range)
    yminplot     = "" # min y axis range (need to set both min and max. Leave it as is for full range)
    ymaxplot     = "" # max y axis range (need to set both min and max. Leave it as is for full range)
    lpos         = "" # legend position (default = top-right, option="bottom-center", "top-left")
    #    xlog         = "" # log scale of X axis (default = no, option="yes") ### IT SEEMS IT DOES NOT WORK
    ylog         = "" # log scale of Y axis (default = no, option="yes")
    #rebin       = "" # rebin x axis (default = 1, option = set it to whatever you want )
    name         = "" # name of the final plots
    lint         = "2.6 fb^{-1}" # integrated luminosity of the sample ( example "10 pb^{-1}" )
    fileXsectionNoRescale = "" #cross section file (with no rescale
    datasetName = "" # string for pattern recognition of dataset name (rescaling will be done only on matched datasets)

    def CheckMCDataConsistency(self):
        #checks
        if(self.histoMCall.GetNbinsX()!=self.histoDATA.GetNbinsX()):
            print "WARNING! number of bins is different between DATA and MC"
            print "exiting..."
            sys.exit()
        if(self.histoMCall.GetBinWidth(1)!=self.histoDATA.GetBinWidth(1)):
            print "WARNING! bin width is different between DATA and MC"
            print "exiting..."
            sys.exit()


def CalculateRescaleFactor(plotObjTTBar, plotObjWJets, fileps):
    #calculate rescaling factor for Z/gamma+jet background and create new cross section file
    canvas = TCanvas()

    plotObjTTBar.CheckMCDataConsistency()
    plotObjWJets.CheckMCDataConsistency()

    #integrals: ttbar
    integralDATA_ttbar = GetIntegralTH1(plotObjTTBar.histoDATA,plotObjTTBar.xmin,plotObjTTBar.xmax)
    ERRintegralDATA_ttbar = GetErrorIntegralTH1(plotObjTTBar.histoDATA,plotObjTTBar.xmin,plotObjTTBar.xmax)
    integralMCall_ttbar = GetIntegralTH1(plotObjTTBar.histoMCall,plotObjTTBar.xmin,plotObjTTBar.xmax)
    ERRintegralMCall_ttbar = GetErrorIntegralTH1(plotObjTTBar.histoMCall,plotObjTTBar.xmin,plotObjTTBar.xmax)
    integralTTbar_ttbar = GetIntegralTH1(plotObjTTBar.histoTTbar,plotObjTTBar.xmin,plotObjTTBar.xmax)
    ERRintegralTTbar_ttbar = GetErrorIntegralTH1(plotObjTTBar.histoTTbar,plotObjTTBar.xmin,plotObjTTBar.xmax)
    integralWJets_ttbar = GetIntegralTH1(plotObjTTBar.histoWJet,plotObjTTBar.xmin,plotObjTTBar.xmax)
    ERRintegralWJets_ttbar = GetErrorIntegralTH1(plotObjTTBar.histoWJet,plotObjTTBar.xmin,plotObjTTBar.xmax)
    integralQCD_ttbar = GetIntegralTH1(plotObjTTBar.histoQCD,plotObjTTBar.xmin,plotObjTTBar.xmax)
    ERRintegralQCD_ttbar = GetErrorIntegralTH1(plotObjTTBar.histoQCD,plotObjTTBar.xmin,plotObjTTBar.xmax)
    #contamination from other backgrounds (except TTbar and WJets) in the integral range [QCD is not in MCall]
    integralMCothers_ttbar = integralMCall_ttbar - integralTTbar_ttbar - integralWJets_ttbar
    ERRintegralMCothers_ttbar = math.sqrt(ERRintegralMCall_ttbar**2 + ERRintegralTTbar_ttbar**2)
    contamination_ttbar = (integralMCothers_ttbar+integralWJets_ttbar+integralQCD_ttbar) / (integralMCall_ttbar+integralQCD_ttbar)

    # integrals: wjets
    integralDATA_wjets = GetIntegralTH1(plotObjWJets.histoDATA,plotObjWJets.xmin,plotObjWJets.xmax)
    ERRintegralDATA_wjets = GetErrorIntegralTH1(plotObjWJets.histoDATA,plotObjWJets.xmin,plotObjWJets.xmax)
    integralMCall_wjets = GetIntegralTH1(plotObjWJets.histoMCall,plotObjWJets.xmin,plotObjWJets.xmax)
    ERRintegralMCall_wjets = GetErrorIntegralTH1(plotObjWJets.histoMCall,plotObjWJets.xmin,plotObjWJets.xmax)
    integralTTbar_wjets = GetIntegralTH1(plotObjWJets.histoTTbar,plotObjWJets.xmin,plotObjWJets.xmax)
    ERRintegralTTbar_wjets = GetErrorIntegralTH1(plotObjWJets.histoTTbar,plotObjWJets.xmin,plotObjWJets.xmax)
    integralWJets_wjets = GetIntegralTH1(plotObjWJets.histoWJet,plotObjWJets.xmin,plotObjWJets.xmax)
    ERRintegralWJets_wjets = GetErrorIntegralTH1(plotObjWJets.histoWJet,plotObjWJets.xmin,plotObjWJets.xmax)
    integralQCD_wjets = GetIntegralTH1(plotObjWJets.histoQCD,plotObjWJets.xmin,plotObjWJets.xmax)
    ERRintegralQCD_wjets = GetErrorIntegralTH1(plotObjWJets.histoQCD,plotObjWJets.xmin,plotObjWJets.xmax)
    #contamination from other backgrounds (except WJets and TTBar) in the integral range [QCD is not in MCall]
    integralMCothers_wjets = integralMCall_wjets - integralWJets_wjets - integralTTbar_wjets
    ERRintegralMCothers_wjets = math.sqrt(ERRintegralMCall_wjets**2 + ERRintegralWJets_wjets**2)
    contamination_wjets = (integralMCothers_wjets+integralTTbar_wjets+integralQCD_wjets) / (integralMCall_wjets+integralQCD_wjets)

    # solve the system of equations
    # (1) --> wjets
    # (2) --> ttbar
    rTTBar = integralDATA_wjets*integralWJets_ttbar - (integralMCothers_wjets+integralQCD_wjets)*integralWJets_ttbar
    rTTBar+= (integralMCothers_ttbar+integralQCD_ttbar-integralDATA_ttbar)*integralWJets_wjets
    rTTBar/= (integralTTbar_wjets*integralWJets_ttbar - integralTTbar_ttbar*integralWJets_wjets)
    rWJets = integralDATA_wjets*integralTTbar_ttbar - (integralMCothers_wjets+integralQCD_wjets)*integralTTbar_ttbar
    rWJets+= (integralMCothers_ttbar+integralQCD_ttbar-integralDATA_ttbar)*integralTTbar_wjets
    rWJets/= (integralTTbar_ttbar*integralWJets_wjets-integralTTbar_wjets*integralWJets_ttbar)


    ##DATA corrected for other bkg contamination --> best estimate of DATA (due to Z only)
    #integralDATAcorr = (integralDATA - integralMCothers)
    #ERRintegralDATAcorr = math.sqrt(ERRintegralDATA**2 + ERRintegralMCothers**2)

    ##rescale factor
    #rescale = integralDATAcorr / integralTTbar
    #relERRintegralDATAcorr = ERRintegralDATAcorr / integralDATAcorr
    #relERRintegralTTbar = ERRintegralTTbar / integralTTbar
    #relERRrescale = math.sqrt(relERRintegralDATAcorr**2 + relERRintegralTTbar**2)

# FIXME
    ##draw histo
    #self.histoMCall.SetFillColor(kBlue)
    #self.histoDATA.SetMarkerStyle(20)

    #self.histoMCall.Draw("HIST")
    #self.histoDATA.Draw("psame")
    #self.histoMCall.GetXaxis().SetRangeUser(self.xminplot,self.xmaxplot)
    #self.histoMCall.GetYaxis().SetRangeUser(self.yminplot,self.ymaxplot)

    #canvas.Update()
    #gPad.RedrawAxis()
    #gPad.Modified()
    ##canvas.SaveAs(self.name + ".eps","eps")
    ##canvas.SaveAs(self.name + ".pdf","pdf")
    #canvas.Print(fileps)
    #canvas.Print(self.name + ".C")
    ## make root file
    #tfile = TFile(self.name+'.root','recreate')
    #tfile.cd()
    #self.histoDATA.Write()
    #self.histoTTbar.Write()
    #self.histoMCall.Write()
    #self.histoQCD.Write()
    #self.histoZJet.Write()
    #self.histoWJet.Write()
    #self.histoSingleTop.Write()
    #self.histoPhotonJets.Write()
    #self.histoDiboson.Write()
    #tfile.Close()

    #printout
    print
    print " TTBar "
    print "######################################## "
    print "integral range:               " + str(plotObjTTBar.xmin) + " < MTenu < " + str(plotObjTTBar.xmax) + " GeV/c2"
    print "integral MC All:              "   + str( integralMCall_ttbar ) + " +/- " + str( ERRintegralMCall_ttbar )
    print "integral QCD:                 "   + str( integralQCD_ttbar ) + " +/- " + str( ERRintegralQCD_ttbar )
    print "integral MC TTbar:            "   + str( integralTTbar_ttbar) + " +/- " + str( ERRintegralTTbar_ttbar )
    print "integral MC WJets:            "   + str( integralWJets_ttbar) + " +/- " + str( ERRintegralWJets_ttbar )
    print "integral MC other:            "   + str( integralMCothers_ttbar) + " +/- " + str( ERRintegralMCothers_ttbar )
    print "rescaled integral MC TTbar:   "   + str( rTTBar*integralTTbar_ttbar) + " +/- " + str( rTTBar*ERRintegralTTbar_ttbar )
    print "rescaled integral MC WJets:   "   + str( rWJets*integralWJets_ttbar) + " +/- " + str( rWJets*ERRintegralWJets_ttbar )
    print "integral DATA:                "   + str( integralDATA_ttbar ) + " +/- " + str( ERRintegralDATA_ttbar )
    print "contribution from other bkgs (except TTbar): " + str(contamination_ttbar*100) + "%"
    #print "integral DATA (corrected for contribution from other bkgs): "  + str( integralDATAcorr_ttbar ) + " +/- " + str( ERRintegralDATAcorr_ttbar )
    print "rescale factor for TTbar background: " + str(rTTBar) #+ " +\- " + str(relERRrescale*rescale)
    # FIXME
    #print "systematical uncertainty of TTbar background modeling: " + str(relERRrescale*100) + "%"
    print "######################################## "
    print

    print
    print " WJets "
    print "######################################## "
    print "integral range:               " + str(plotObjWJets.xmin) + " < MTenu < " + str(plotObjWJets.xmax) + " GeV/c2"
    print "integral MC All:              "   + str( integralMCall_wjets ) + " +/- " + str( ERRintegralMCall_wjets )
    print "integral QCD:                 "   + str( integralQCD_wjets ) + " +/- " + str( ERRintegralQCD_wjets )
    print "integral MC TTbar:            "   + str( integralTTbar_wjets) + " +/- " + str( ERRintegralTTbar_wjets )
    print "integral MC WJets:            "   + str( integralWJets_wjets) + " +/- " + str( ERRintegralWJets_wjets )
    print "integral MC other:            "   + str( integralMCothers_wjets) + " +/- " + str( ERRintegralMCothers_wjets )
    print "rescaled integral MC TTbar:   "   + str( rTTBar*integralTTbar_wjets) + " +/- " + str( rTTBar*ERRintegralTTbar_wjets )
    print "rescaled integral MC WJets:   "   + str( rWJets*integralWJets_wjets) + " +/- " + str( rWJets*ERRintegralWJets_wjets )
    print "integral DATA:                "   + str( integralDATA_wjets ) + " +/- " + str( ERRintegralDATA_wjets )
    print "contribution from other bkgs (except wjets): " + str(contamination_wjets*100) + "%"
    #print "integral DATA (corrected for contribution from other bkgs): "  + str( integralDATAcorr_wjets ) + " +/- " + str( ERRintegralDATAcorr_wjets )
    print "rescale factor for WJets background: " + str(rWJets) #+ " +\- " + str(relERRrescale*rescale)
    # FIXME
    #print "systematical uncertainty of WJets background modeling: " + str(relERRrescale*100) + "%"
    print "######################################## "
    print

#TODO FIXME
    ##create new cross section file
    #originalFileName = string.split( string.split(self.fileXsectionNoRescale, "/" )[-1], "." ) [0]
    #newFileName = originalFileName + "_" + self.name +".txt"
    #os.system('rm -f '+ newFileName)
    #outputFile = open(newFileName,'w')

    #for line in open( self.fileXsectionNoRescale ):
    #    line = string.strip(line,"\n")

    #    if( re.search(self.datasetName, line) ):
    #        list = re.split( '\s+' , line  )
    #        newline = str(list[0]) + "    "  + str("%.6f" % (float(list[1])*float(rescale)) )
    #        print >> outputFile, newline
    #    else:
    #        print >> outputFile, line

    #outputFile.close
    #print "New xsection file (after TTbar rescaling) is: " + newFileName
    #print " "


############# DON'T NEED TO MODIFY ANYTHING HERE - END #######################
##############################################################################


##############################################################################
############# USER CODE - BEGIN ##############################################

#--- Input files
#preselection
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_onRSK_local_nov1_addStSFplots_allDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_orig.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_onRSK_local_nov1_addStSFplots_allDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_QCD_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_local_nov18_addStSFplots_allDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")

#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_local_nov18_addStSFplots_ICHEPDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_local_nov18_addStSFplots_ICHEPDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_QCD_plots.root")

File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_dec3_ICHEPDataExcludeEarlyRunsAndMC_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_dec3_ICHEPDataExcludeEarlyRunsAndMC_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_QCD_plots.root")

#--- Rescaling of W + jet and ttbar+jets background

plots = []
#-----------------------------------------
# for ttbar-enriched region
plotBaseName = 'MTenu_50_110_Njet_gte4'

# MG HT BKG
h_ALLBKG_HT_ttbar = GetHisto("histo1D__ALLBKG_MG_HT__"+plotBaseName, File_preselection) # MC all
# amc@NLO BKG
h_ALLBKG_amcatnlo_ttbar = GetHisto("histo1D__ALLBKG_amcAtNLOIncTTBar_ZJetWJetPt__"+plotBaseName, File_preselection) # MC all

h_TTbar_MG_ttbar = GetHisto("histo1D__TTbar_Madgraph__"+plotBaseName, File_preselection) # MC TTbar
h_ZJets_MGHT_ttbar = GetHisto("histo1D__ZJet_Madgraph_HT__"+plotBaseName, File_preselection)
h_WJets_MGHT_ttbar = GetHisto("histo1D__WJet_Madgraph_HT__"+plotBaseName, File_preselection)
h_TTbar_amcatnlo_ttbar = GetHisto("histo1D__TTbar_amcatnlo_Inc__"+plotBaseName, File_preselection) # MC TTbar
h_ZJets_amcatnlo_ttbar = GetHisto("histo1D__ZJet_amcatnlo_ptBinned__"+plotBaseName, File_preselection)
h_WJets_amcatnlo_ttbar = GetHisto("histo1D__WJet_amcatnlo_ptBinned__"+plotBaseName, File_preselection)
h_SingleTop_ttbar = GetHisto("histo1D__SingleTop__"+plotBaseName, File_preselection)
h_PhotonJets_ttbar = GetHisto("histo1D__PhotonJets_Madgraph__"+plotBaseName, File_preselection)
h_Diboson_ttbar = GetHisto("histo1D__DIBOSON__"+plotBaseName, File_preselection)

# DATA
h_DATA_ttbar = GetHisto("histo1D__DATA__"+plotBaseName, File_preselection) #DATA
# QCD
#h_QCD_DataDriven = GetHisto("histo1D__QCDFakes_DATA__"+plotBaseName,File_QCD_preselection)
h_QCD_ttbar = GetHisto("histo1D__QCD_EMEnriched__"+plotBaseName,File_QCD_preselection)

plotTTbar = Plot()
plotTTbar.histoDATA = h_DATA_ttbar
useMGHT = False
# amc@NLO
if not useMGHT:
  plotTTbar.histoMCall = h_ALLBKG_amcatnlo_ttbar
  plotTTbar.histoTTbar = h_TTbar_amcatnlo_ttbar
  plotTTbar.histoQCD = h_QCD_ttbar
  plotTTbar.histoZJet = h_ZJets_amcatnlo_ttbar
  plotTTbar.histoWJet = h_WJets_amcatnlo_ttbar
else:
  # MG HT
  plotTTbar.histoMCall = h_ALLBKG_HT_ttbar
  plotTTbar.histoTTbar = h_TTbar_MG_ttbar
  plotTTbar.histoQCD = h_QCD_ttbar
  plotTTbar.histoZJet = h_ZJets_MGHT_ttbar
  plotTTbar.histoWJet = h_WJets_MGHT_ttbar

plotTTbar.histoSingleTop = h_SingleTop_ttbar
plotTTbar.histoPhotonJets = h_PhotonJets_ttbar
plotTTbar.histoDiboson = h_Diboson_ttbar
plotTTbar.xmin = 50
#plotTTbar.xmax = h_TTbar_amcatnlo_ttbar.GetXaxis().GetXmax()
plotTTbar.xmax = 110
plotTTbar.name = "TTbarRescale"
plotTTbar.fileXsectionNoRescale = "/afs/cern.ch/user/s/scooper/work/private/cmssw/8011/TestRootNTuplizerRecipe/src/Leptoquarks/analyzer/rootNtupleAnalyzerV2/config/xsection_13TeV_2015.txt"
plotTTbar.xminplot = 0
plotTTbar.xmaxplot = 2000
plotTTbar.yminplot = 0
plotTTbar.ymaxplot = 2000
plotTTbar.datasetName = "TTJets_.+Tune"
#plot0.datasetName = "DYJetsToLL_M-50_HT.+Tune"
#plot0.datasetName = "Z.+Jets_Pt.+alpgen"
# example: this match with /Z3Jets_Pt300to800-alpgen/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
plots.append(plotTTbar)

#-----------------------------------------
# for wjets-enriched region
plotBaseName = 'MTenu_50_110_Njet_lte3'

# MG HT BKG
h_ALLBKG_HT_wjets = GetHisto("histo1D__ALLBKG_MG_HT__"+plotBaseName, File_preselection) # MC all
# amc@NLO BKG
h_ALLBKG_amcatnlo_wjets = GetHisto("histo1D__ALLBKG_amcAtNLOIncTTBar_ZJetWJetPt__"+plotBaseName, File_preselection) # MC all

h_TTbar_MG_wjets = GetHisto("histo1D__TTbar_Madgraph__"+plotBaseName, File_preselection) # MC TTbar
h_ZJets_MGHT_wjets = GetHisto("histo1D__ZJet_Madgraph_HT__"+plotBaseName, File_preselection)
h_WJets_MGHT_wjets = GetHisto("histo1D__WJet_Madgraph_HT__"+plotBaseName, File_preselection)
h_TTbar_amcatnlo_wjets = GetHisto("histo1D__TTbar_amcatnlo_Inc__"+plotBaseName, File_preselection) # MC TTbar
h_ZJets_amcatnlo_wjets = GetHisto("histo1D__ZJet_amcatnlo_ptBinned__"+plotBaseName, File_preselection)
h_WJets_amcatnlo_wjets = GetHisto("histo1D__WJet_amcatnlo_ptBinned__"+plotBaseName, File_preselection)
h_SingleTop_wjets = GetHisto("histo1D__SingleTop__"+plotBaseName, File_preselection)
h_PhotonJets_wjets = GetHisto("histo1D__PhotonJets_Madgraph__"+plotBaseName, File_preselection)
h_Diboson_wjets = GetHisto("histo1D__DIBOSON__"+plotBaseName, File_preselection)

# DATA
h_DATA_wjets = GetHisto("histo1D__DATA__"+plotBaseName, File_preselection) #DATA
# QCD
#h_QCD_DataDriven = GetHisto("histo1D__QCDFakes_DATA__"+plotBaseName,File_QCD_preselection)
h_QCD_wjets = GetHisto("histo1D__QCD_EMEnriched__"+plotBaseName,File_QCD_preselection)

plotWJets = Plot()
plotWJets.histoDATA = h_DATA_wjets
if not useMGHT:
  # amc@NLO
  plotWJets.histoMCall = h_ALLBKG_amcatnlo_wjets
  plotWJets.histoTTbar = h_TTbar_amcatnlo_wjets
  plotWJets.histoQCD = h_QCD_wjets
  plotWJets.histoZJet = h_ZJets_amcatnlo_wjets
  plotWJets.histoWJet = h_WJets_amcatnlo_wjets
else:
  # MG HT
  plotWJets.histoMCall = h_ALLBKG_HT_wjets
  plotWJets.histoTTbar = h_TTbar_MG_wjets
  plotWJets.histoQCD = h_QCD_wjets
  plotWJets.histoZJet = h_ZJets_MGHT_wjets
  plotWJets.histoWJet = h_WJets_MGHT_wjets
  
plotWJets.histoSingleTop = h_SingleTop_wjets
plotWJets.histoPhotonJets = h_PhotonJets_wjets
plotWJets.histoDiboson = h_Diboson_wjets
plotWJets.xmin = 50
#plotWJets.xmax = h_TTbar_amcatnlo_wjets.GetXaxis().GetXmax()
plotWJets.xmax = 110
plotWJets.name = "WJetsRescale"
plotWJets.fileXsectionNoRescale = "/afs/cern.ch/user/s/scooper/work/private/cmssw/8011/TestRootNTuplizerRecipe/src/Leptoquarks/analyzer/rootNtupleAnalyzerV2/config/xsection_13TeV_2015.txt"
plotWJets.xminplot = 0
plotWJets.xmaxplot = 2000
plotWJets.yminplot = 0
plotWJets.ymaxplot = 2000
plotWJets.datasetName = "WJetsToLNu_.+Tune"
#plot0.datasetName = "DYJetsToLL_M-50_HT.+Tune"
#plot0.datasetName = "Z.+Jets_Pt.+alpgen"
# example: this match with /Z3Jets_Pt300to800-alpgen/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO
plots.append(plotWJets)

#-----------------------------------------------------------------------------------


############# USER CODE - END ################################################
##############################################################################



#--- Generate and print the plots from the list 'plots' define above

#--- Output files
fileps = "allPlots_calc_WJetsAndTTBarRescale_And_xsecFile.ps"

#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()
c.Print(fileps+"[")
CalculateRescaleFactor(plotTTbar,plotWJets,fileps)
c.Print(fileps+"]")
os.system('ps2pdf '+fileps)

