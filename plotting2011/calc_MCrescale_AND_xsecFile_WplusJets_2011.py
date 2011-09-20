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


## The Plot class: add members if needed
class Plot:
    histoDATA   = "" # DATA
    histoMCW    = "" # MC W
    histoMCTTbar = "" #MC TTbar
    histoMCZ    = "" # MC Z
    histoMCSingleTop = "" # MC SingleTop
    histoMCWW   = "" # MC WW
    histoMCWZ   = "" # MC WZ
    histoMCZZ   = "" # MC ZZ
    histoMCPhotonJets   = "" # MC PhotonJets
    histoMCBJets   = "" # MC BJets
    # histoQCD    = "" # QCD (data-driven)
    MCSyst     = 0  # Systematic uncertainty on MC samples other than TTbar (0.1 = 10%)
    MCTTbarSyst = 0  # Systematic uncertainty on MC TTbar (0.1 = 10%)
    QCDSyst     = 0  # Systematic uncertainty on data-driven QCD (0.1 = 10%)
    xtit        = "" # xtitle
    ytit        = "" # ytitle
    xmin        = "" # set xmin to calculate rescaling factor (-- please take into account the bin size of histograms --)
    xmax        = "" # # set xmax to calculate rescaling factor (-- please take into account the bin size of histograms --)
    xminplot    = "" # min x axis range (need to set both min and max. Leave it as is for full range)
    xmaxplot    = "" # max x axis range (need to set both min and max. Leave it as is for full range)
    yminplot    = "" # min y axis range (need to set both min and max. Leave it as is for full range)
    ymaxplot    = "" # max y axis range (need to set both min and max. Leave it as is for full range)
    lpos        = "" # legend position (default = top-right, option="bottom-center", "top-left")
    #    xlog        = "" # log scale of X axis (default = no, option="yes") ### IT SEEMS IT DOES NOT WORK
    ylog        = "" # log scale of Y axis (default = no, option="yes")
    #rebin      = "" # rebin x axis (default = 1, option = set it to whatever you want )
    name        = "" # name of the final plots
    lint        = "36.0 pb^{-1}" # integrated luminosity of the sample ( example "10 pb^{-1}" )
    fileXsectionNoRescale = "" #cross section file (with no rescale
    datasetName = "" # string for pattern recognition of dataset name (rescaling will be done only on matched datasets)

    def CalculateRescaleFactor(self, fileps):
        #calculate rescaling factor for W+jets background and create new cross section file
        canvas = TCanvas()

        #check
        if(self.histoMCW.GetNbinsX()!=self.histoDATA.GetNbinsX()):
            print "WARNING! number of bins is different between DATA and MC"
            print "exiting..."
            sys.exit()
        if(self.histoMCW.GetBinWidth(1)!=self.histoDATA.GetBinWidth(1)):
            print "WARNING! bin width is different between DATA and MC"
            print "exiting..."
            sys.exit()

        #integrals
        integralDATA = GetIntegralTH1(self.histoDATA,self.xmin,self.xmax)
        ERRintegralDATA = GetErrorIntegralTH1(self.histoDATA,self.xmin,self.xmax)
        integralMCW = GetIntegralTH1(self.histoMCW,self.xmin,self.xmax)
        ERRintegralMCW = GetErrorIntegralTH1(self.histoMCW,self.xmin,self.xmax)
        integralMCTTbar = GetIntegralTH1(self.histoMCTTbar,self.xmin,self.xmax)
        ERRintegralMCTTbar = sqrt( GetErrorIntegralTH1(self.histoMCTTbar,self.xmin,self.xmax)**2 + (self.MCTTbarSyst*integralMCTTbar)**2 )
        integralMCZ = GetIntegralTH1(self.histoMCZ,self.xmin,self.xmax)
        ERRintegralMCZ = sqrt( GetErrorIntegralTH1(self.histoMCZ,self.xmin,self.xmax)**2 + (self.MCSyst*integralMCZ)**2 )
        integralMCSingleTop = GetIntegralTH1(self.histoMCSingleTop,self.xmin,self.xmax)
        ERRintegralMCSingleTop = sqrt( GetErrorIntegralTH1(self.histoMCSingleTop,self.xmin,self.xmax)**2 + (self.MCSyst*integralMCSingleTop)**2 )
        integralMCWW = GetIntegralTH1(self.histoMCWW,self.xmin,self.xmax)
        ERRintegralMCWW = sqrt( GetErrorIntegralTH1(self.histoMCWW,self.xmin,self.xmax)**2 + (self.MCSyst*integralMCWW)**2 )
        integralMCWZ = GetIntegralTH1(self.histoMCWZ,self.xmin,self.xmax)
        ERRintegralMCWZ = sqrt( GetErrorIntegralTH1(self.histoMCWZ,self.xmin,self.xmax)**2 + (self.MCSyst*integralMCWZ)**2 )
        integralMCZZ = GetIntegralTH1(self.histoMCZZ,self.xmin,self.xmax)
        ERRintegralMCZZ = sqrt( GetErrorIntegralTH1(self.histoMCZZ,self.xmin,self.xmax)**2 + (self.MCSyst*integralMCZZ)**2 )
        integralMCPhotonJets = GetIntegralTH1(self.histoMCPhotonJets,self.xmin,self.xmax)
        ERRintegralMCPhotonJets = sqrt( GetErrorIntegralTH1(self.histoMCPhotonJets,self.xmin,self.xmax)**2 + (self.MCSyst*integralMCPhotonJets)**2 )
        # integralMCBJets = GetIntegralTH1(self.histoMCBJets,self.xmin,self.xmax)
        # ERRintegralMCBJets = sqrt( GetErrorIntegralTH1(self.histoMCBJets,self.xmin,self.xmax)**2 + (self.MCSyst*integralMCBJets)**2 )
        integralMCBJets = 0.0
        ERRintegralMCBJets = 0.0
        # integralQCD = GetIntegralTH1(self.histoQCD,self.xmin,self.xmax)
        # ERRintegralQCD = sqrt( GetErrorIntegralTH1(self.histoQCD,self.xmin,self.xmax)**2 + (self.QCDSyst*integralQCD)**2 )
        integralQCD = 0.0
        ERRintegralQCD = 0.0

        #contamination from other backgrounds (other than W+jets) in the integral range
        integralOthers =  integralMCTTbar + integralMCZ + integralMCSingleTop + integralMCWW + integralMCWZ + integralMCZZ + integralMCPhotonJets + integralMCBJets + integralQCD
        ERRintegralOthers = sqrt( ERRintegralMCTTbar**2 + ERRintegralMCZ**2 + ERRintegralMCSingleTop**2 + ERRintegralMCWW**2 + ERRintegralMCWZ**2 + ERRintegralMCZZ**2 + ERRintegralMCPhotonJets**2 + ERRintegralMCBJets**2 + ERRintegralQCD**2 )
        integralMCall = integralMCW+ integralMCTTbar + integralMCZ  + integralMCSingleTop + integralMCWW + integralMCWZ + integralMCZZ + integralMCPhotonJets + integralMCBJets
        ERRintegralMCall = sqrt( ERRintegralMCW**2 + ERRintegralMCTTbar**2 + ERRintegralMCZ**2 + ERRintegralMCSingleTop**2 + ERRintegralMCWW**2 + ERRintegralMCWZ**2 + ERRintegralMCZZ**2 + ERRintegralMCPhotonJets**2 + ERRintegralMCBJets**2 )
        contamination = integralOthers / (integralMCW + integralOthers)

        #DATA corrected for other bkg contamination --> best estimate of DATA (due to W only)
        integralDATAcorr = (integralDATA - integralOthers)
        ERRintegralDATAcorr = sqrt(ERRintegralDATA**2 + ERRintegralOthers**2)

        #rescale factor
        rescale = integralDATAcorr / integralMCW
        relERRintegralDATAcorr = ERRintegralDATAcorr / integralDATAcorr
        relERRintegralMCW = ERRintegralMCW / integralMCW
        relERRrescale = sqrt(relERRintegralDATAcorr**2 + relERRintegralMCW**2)

        #draw histo
        histoAll = copy.deepcopy(self.histoMCW)
        histoAll.Add(self.histoMCTTbar)
        histoAll.Add(self.histoMCZ)
        histoAll.Add(self.histoMCSingleTop)
        histoAll.Add(self.histoMCWW)
        histoAll.Add(self.histoMCWZ)
        histoAll.Add(self.histoMCZZ)
        histoAll.Add(self.histoMCPhotonJets)
        # histoAll.Add(self.histoMCBJets)
        # histoAll.Add(self.histoQCD)
        histoAll.SetFillColor(kBlue)
        self.histoDATA.SetMarkerStyle(20)

        histoAll.Draw("HIST")
        self.histoDATA.Draw("psame")
        histoAll.GetXaxis().SetRangeUser(self.xminplot,self.xmaxplot)
        histoAll.GetYaxis().SetRangeUser(self.yminplot,self.ymaxplot)

        canvas.Update()
        gPad.RedrawAxis()
        gPad.Modified()
        #canvas.SaveAs(self.name + ".eps","eps")
        #canvas.SaveAs(self.name + ".pdf","pdf")
        canvas.Print(fileps)

        #printout
        print " "
        print "######################################## "
        print "integral range: " + str(self.xmin) + " < MTenu < " + str(self.xmax) + " GeV/c2"
        print "integral DATA: "   + str( integralDATA ) + " +/- " + str( ERRintegralDATA )
        print "integral MC All: "   + str( integralMCall ) + " +/- " + str( ERRintegralMCall )
        print "integral QCD: "   + str( integralQCD ) + " +/- " + str( ERRintegralQCD )
        print "integral MC W: "   + str( integralMCW ) + " +/- " + str( ERRintegralMCW )
        print "contribution from other bkgs (other than W+jets): " + str(contamination*100) + "%"
        print "integral DATA (corrected for contribution from other bkgs): "  + str( integralDATAcorr ) + " +/- " + str( ERRintegralDATAcorr )
        print "rescale factor for W background: " + str(rescale) + " +/- " + str(relERRrescale*rescale)
        print "systematical uncertainty of W+jets background modeling: " + str(relERRrescale*100) + "%"
        print "######################################## "
        print " "

        #create new cross section file
        originalFileName = string.split( string.split(self.fileXsectionNoRescale, "/" )[-1], "." ) [0]
        newFileName = originalFileName + "_" + self.name +".txt"
        os.system('rm -f '+ newFileName)
        outputFile = open(newFileName,'w')

        for line in open( self.fileXsectionNoRescale ):
            line = string.strip(line,"\n")

            if( re.search(self.datasetName, line) ):
                list = re.split( '\s+' , line  )
                newline = str(list[0]) + "    "  + str("%.6f" % (float(list[1])*float(rescale)) )
                print >> outputFile, newline
            else:
                print >> outputFile, line

        outputFile.close
        print "New xsection file (after W rescaling) is: " + newFileName
        print " "


############# DON'T NEED TO MODIFY ANYTHING HERE - END #######################
##############################################################################


##############################################################################
############# USER CODE - BEGIN ##############################################

#--- Input files
#preselection
File_preselection = GetFile("/Users/eberry/Code/ROOT/LQ/output/triggerStudy/analysisClass_lq_enujj_forTriggerStudy_plots.root")
# File_preselection_QCD = GetFile("/afs/cern.ch/user/f/ferencek/scratch0/LQ/CMSSW_3_8_6/test/Leptoquarks/rootNtupleAnalyzerV2/data/output/analysisClass_enujjSample_QCD_plots.root")

#--- Rescaling of W+jets background

#-----------------------------------------
h_DATA_MTenu             = GetHisto("histo1D__DATA__MTenu_PAS", File_preselection) #DATA
h_WJetAlpgen_MTenu       = GetHisto("histo1D__WJet_Madgraph__MTenu_PAS", File_preselection) # MC W
h_TTbar_Madgraph_MTenu   = GetHisto("histo1D__TTbar_Madgraph__MTenu_PAS", File_preselection) # MC TTbar
h_ZJetAlpgen_MTenu       = GetHisto("histo1D__ZJet_Madgraph__MTenu_PAS", File_preselection) # MC Z
h_SingleTop_MTenu        = GetHisto("histo1D__SingleTop__MTenu_PAS", File_preselection) # MC SingleTop
h_WW_MTenu               = GetHisto("histo1D__WW__MTenu_PAS", File_preselection) # MC WW
h_WZ_MTenu               = GetHisto("histo1D__WZ__MTenu_PAS", File_preselection) # MC WZ
h_ZZ_MTenu               = GetHisto("histo1D__ZZ__MTenu_PAS", File_preselection) # MC ZZ
h_PhotonJets_MTenu       = GetHisto("histo1D__PhotonJets__MTenu_PAS", File_preselection) # MC PhotonJets
# h_BJets_MTenu            = GetHisto("histo1D__BJets__MTenu_PAS", File_preselection) # MC BJets
# h_QCD_MTenu = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________MTenu_PAS", File_preselection_QCD, 1. ) #QCD (data-driven) scaled to the correct integrated luminosity (36.05/35.84~=1)

##-- W- only
#h_DATA_MTenu = GetHisto("histo1D__DATA__h1_MTenu_PAS_minus", File_preselection) #DATA
#h_WJetAlpgen_MTenu = GetHisto("histo1D__WJetAlpgen__h1_MTenu_PAS_minus", File_preselection) # MC W
#h_ZJetAlpgen_MTenu = GetHisto("histo1D__ZJetAlpgen__h1_MTenu_PAS_minus", File_preselection) # MC Z
#h_TTbar_Madgraph_MTenu = GetHisto("histo1D__TTbar_Madgraph__h1_MTenu_PAS_minus", File_preselection) # MC TTbar
#h_SingleTop_MTenu = GetHisto("histo1D__SingleTop__h1_MTenu_PAS_minus", File_preselection) # MC SingleTop
#h_WW_MTenu = GetHisto("histo1D__WW__h1_MTenu_PAS_minus", File_preselection) # MC WW
#h_WZ_MTenu = GetHisto("histo1D__WZ__h1_MTenu_PAS_minus", File_preselection) # MC WZ
#h_ZZ_MTenu = GetHisto("histo1D__ZZ__h1_MTenu_PAS_minus", File_preselection) # MC ZZ
#h_PhotonJets_MTenu = GetHisto("histo1D__PhotonJets__h1_MTenu_PAS_minus", File_preselection) # MC PhotonJets
#h_BJets_MTenu = GetHisto("histo1D__BJets__h1_MTenu_PAS_minus", File_preselection) # MC BJets
#h_QCD_MTenu = GetHisto("histo1D__DATA__h1_MTenu_PAS_minus", File_preselection_QCD, 1. ) #QCD (data-driven) scaled to the correct integrated luminosity (36.05/35.84~=1)

##-- W+ only
#h_DATA_MTenu = GetHisto("histo1D__DATA__h1_MTenu_PAS_plus", File_preselection) #DATA
#h_WJetAlpgen_MTenu = GetHisto("histo1D__WJetAlpgen__h1_MTenu_PAS_plus", File_preselection) # MC W
#h_ZJetAlpgen_MTenu = GetHisto("histo1D__ZJetAlpgen__h1_MTenu_PAS_plus", File_preselection) # MC Z
#h_TTbar_Madgraph_MTenu = GetHisto("histo1D__TTbar_Madgraph__h1_MTenu_PAS_plus", File_preselection) # MC TTbar
#h_SingleTop_MTenu = GetHisto("histo1D__SingleTop__h1_MTenu_PAS_plus", File_preselection) # MC SingleTop
#h_WW_MTenu = GetHisto("histo1D__WW__h1_MTenu_PAS_plus", File_preselection) # MC WW
#h_WZ_MTenu = GetHisto("histo1D__WZ__h1_MTenu_PAS_plus", File_preselection) # MC WZ
#h_ZZ_MTenu = GetHisto("histo1D__ZZ__h1_MTenu_PAS_plus", File_preselection) # MC ZZ
#h_PhotonJets_MTenu = GetHisto("histo1D__PhotonJets__h1_MTenu_PAS_plus", File_preselection) # MC PhotonJets
#h_BJets_MTenu = GetHisto("histo1D__BJets__h1_MTenu_PAS_plus", File_preselection) # MC BJets
#h_QCD_MTenu = GetHisto("histo1D__DATA__h1_MTenu_PAS_plus", File_preselection_QCD, 1. ) #QCD (data-driven) scaled to the correct integrated luminosity (36.05/35.84~=1)


plot0 = Plot()
plot0.histoDATA = h_DATA_MTenu
plot0.histoMCW = h_WJetAlpgen_MTenu
plot0.histoMCTTbar = h_TTbar_Madgraph_MTenu
plot0.histoMCZ = h_ZJetAlpgen_MTenu
plot0.histoMCSingleTop = h_SingleTop_MTenu
plot0.histoMCWW = h_WW_MTenu
plot0.histoMCWZ = h_WZ_MTenu
plot0.histoMCZZ = h_ZZ_MTenu
plot0.histoMCPhotonJets = h_PhotonJets_MTenu
# plot0.histoMCBJets = h_BJets_MTenu
# plot0.histoQCD = h_QCD_MTenu
plot0.MCSyst = 0.11
plot0.MCTTbarSyst = 0.17
plot0.QCDSyst = 0.25
plot0.xmin = 50
plot0.xmax = 110
plot0.name = "Wrescale"
plot0.fileXsectionNoRescale = "/Users/eberry/Code/ROOT/LQ/dev_LQANA/rootNtupleAnalyzerV2/config/xsection_7TeV_2011.txt"
plot0.xminplot = 0
plot0.xmaxplot = 200
plot0.yminplot = 0
plot0.ymaxplot = 120
plot0.datasetName = "WJetsToLNu_TuneZ2_7TeV"
# example: this match with /W3Jets_Pt300to800-alpgen/Spring10-START3X_V26_S09-v1/GEN-SIM-RECO

plots = [plot0]

#-----------------------------------------------------------------------------------


############# USER CODE - END ################################################
##############################################################################



#--- Generate and print the plots from the list 'plots' define above

#--- Output files
fileps = "WplusJets_MCrescale.ps"

#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()
c.Print(fileps+"[")
for plot in plots:
    plot.CalculateRescaleFactor(fileps)
c.Print(fileps+"]")
os.system('ps2pdf '+fileps)

