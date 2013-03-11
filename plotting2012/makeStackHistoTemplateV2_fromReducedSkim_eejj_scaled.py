#!/usr/bin/env python

from plot_class import *
from ROOT import *

File_preselection = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj//scaled_output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root" )
File_QCD          = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj_qcd//output_cutTable_lq_eejj/analysisClass_lq_eejj_QCD_plots.root")

zUncBand="no"
makeRatio=1
makeNSigma=1

pt_rebin = 2

histoBaseName = "histo1D__SAMPLE__cutHisto_allOtherCuts___________VARIABLE"
histoBaseName_userDef = "histo1D__SAMPLE__VARIABLE"

histoBaseName2D = "histo2D__SAMPLE__cutHisto_allOtherCuts___________VARIABLE"
histoBaseName2D_userDef = "histo2D__SAMPLE__VARIABLE"

samplesForStackHistos_other = [ "OTHERBKG" ]
samplesForStackHistosQCD     = ["DATA"]
samplesForStackHistos_ZJets  = [ "TTbar_Madgraph", "ZJet_Madgraph" ]
samplesForStackHistos = samplesForStackHistos_other + samplesForStackHistos_ZJets

keysStack             = [ "QCD multijets", "Other backgrounds", "t#bar{t}"  ,  "Z/#gamma* + jets"  ]
stackColorIndexes     = [ 432            , 9                  , 600         ,  kRed           ]
stackFillStyleIds     = [ 3005           , 3008               , 3004        ,  3345           ]

stackColorIndexes.reverse()
stackFillStyleIds.reverse()

samplesForHistos = ["LQ_M900"      ]
keys             = ["LQ, M=900 GeV"]


samplesForHistos_blank = []
keys_blank             = []

sampleForDataHisto = "DATA"

QCDScale   = 1.0

def makeDefaultPlot ( variableName, histoBaseName, 
                      samplesForHistos, keys,
                      samplesForStackHistos, keysStack,
                      sampleForDataHisto,
                      zUncBand, makeRatio ) :
    plot                   = Plot() 
    plot.histosStack       =  ( generateHistoList( histoBaseName, samplesForStackHistosQCD   , variableName, File_QCD, QCDScale) + 
                                generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection ) + 
                                generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection )  ) 
    plot.keysStack         = keysStack
    plot.histos            = generateHistoList( histoBaseName, samplesForHistos_blank, variableName, File_preselection)
    plot.keys              = keys_blank
    plot.addZUncBand       = zUncBand
    plot.makeRatio         = makeRatio
    plot.makeNSigma        = makeNSigma
    plot.histodata         = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.ytit              = "N(Events)"
    plot.ylog              = "no"
    plot.name              = variableName
    plot.stackColorIndexes = stackColorIndexes
    plot.stackFillStyleIds = stackFillStyleIds 
    plot.gif_folder        = "gif_eejj_scaled/"
    plot.eps_folder        = "eps_eejj_scaled/"
    plot.pdf_folder        = "pdf_eejj_scaled/"
    plot.png_folder        = "png_eejj_scaled/"
    plot.suffix            = "eejj"
    plot.lumi_fb           = "19.5"
    
    return plot

def makeDefaultPlot2D ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2D()
    plot.hasData        = True
    plot.name           = variableName
    plot.histosStack    = ( generateHistoList( histoBaseName, samplesForStackHistosQCD   , variableName, File_QCD, QCDScale) + 
                            generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection ) + 
                            generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection )  ) 
    plot.histodata      =   generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_eejj_scaled/"
    plot.eps_folder     = "eps_eejj_scaled/"
    plot.suffix         = "eejj"
    
    return plot


def makeDefaultPlot2D_NoData ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2D()
    plot.hasData        = False
    plot.name           = variableName 
    plot.histosStack    = ( generateHistoList ( histoBaseName, samplesForStackHistosQCD   , variableName, File_QCD, QCDScale) + 
                            generateHistoList ( histoBaseName, samplesForStackHistos_other, variableName, File_preselection ) + 
                            generateHistoList ( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection )  ) 
    plot.histodata      =   generateHistoBlank( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_eejj_scaled/"
    plot.eps_folder     = "eps_eejj_scaled/"
    plot.suffix         = "eejj"
    
    return plot


plots = []

plots.append ( makeDefaultPlot ( "nElectron_PAS"      , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax  = 10000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = -0.5
plots[-1].xmax  = 6.5
plots[-1].ylog  = "yes"
plots[-1].xtit  = "Number of electrons [Preselection]"

plots.append ( makeDefaultPlot ( "nMuon_PAS"      , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax  = 10000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = -0.5
plots[-1].xmax  = 6.5
plots[-1].ylog  = "yes"
plots[-1].xtit  = "Number of muons [Preselection]"

plots.append ( makeDefaultPlot ( "nJet_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit  = "Number of jets [Preselection]"
plots[-1].ymax  = 200000000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin  = -0.5
plots[-1].xmax  = 10.5


plots.append ( makeDefaultPlot ( "nJet_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit  = "Number of jets [Preselection]"
plots[-1].ymax  = 10000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "no"
plots[-1].xmin  = -0.5
plots[-1].xmax  = 10.5


plots.append ( makeDefaultPlot ( "nJet_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit  = "Number of jets [Preselection + M(ee) > 100]"
plots[-1].ymax  = 200000000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin  = -0.5
plots[-1].xmax  = 10.5


plots.append ( makeDefaultPlot ( "nJet_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit  = "Number of jets [Preselection + M(ee) > 100]"
plots[-1].ymax  = 2000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "no"
plots[-1].xmin  = -0.5
plots[-1].xmax  = 10.5


plots.append ( makeDefaultPlot ( "nJet_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit  = "Number of jets [Preselection + region of interest]"
plots[-1].ymax  = 200000000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin  = -0.5
plots[-1].xmax  = 10.5



plots.append ( makeDefaultPlot ( "nJet_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit  = "Number of jets [Preselection + region of interest]"
plots[-1].ymax  = 40
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "no"
plots[-1].xmin  = -0.5
plots[-1].xmax  = 10.5




plots.append ( makeDefaultPlot ( "EleChargeSum_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit  = "Electron 1 charge + Electron 2 charge [Preselection]"
plots[-1].ymax  = 200000000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin  = -1.5
plots[-1].xmax  = 1.5


plots.append ( makeDefaultPlot ( "EleChargeSum_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit  = "Electron 1 charge + Electron 2 charge [Preselection + M(ee) > 100]"
plots[-1].ymax  = 200000000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin  = -1.5
plots[-1].xmax  = 1.5


plots.append ( makeDefaultPlot ( "EleChargeSum_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit  = "Electron 1 charge + Electron 2 charge [Preselection + region of interest]"
plots[-1].ymax  = 200000000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin  = -1.5
plots[-1].xmax  = 1.5



plots.append ( makeDefaultPlot ( "EleChargeSum_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit  = "Electron 1 charge + Electron 2 charge [Preselection]"
plots[-1].ymax  = 20000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = -1.5
plots[-1].xmax  = 1.5


plots.append ( makeDefaultPlot ( "EleChargeSum_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit  = "Electron 1 charge + Electron 2 charge [Preselection + M(ee) > 100]"
plots[-1].ymax  = 4000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = -1.5
plots[-1].xmax  = 1.5


plots.append ( makeDefaultPlot ( "EleChargeSum_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit  = "Electron 1 charge + Electron 2 charge [Preselection + region of interest]"
plots[-1].ymax  = 100
plots[-1].ymin  = 1e-1
plots[-1].xmin  = -1.5
plots[-1].xmax  = 1.5


plots.append ( makeDefaultPlot ( "Pt1stEle_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "1st Electron p_{T} (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Pt2ndEle_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "2nd Electron p_{T} (GeV) [Preselection]"


plots.append ( makeDefaultPlot ( "Pt1stEle_PASandMee100"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "1st Electron p_{T} (GeV) [Preselection + M(ee) > 100]"

plots.append ( makeDefaultPlot ( "Pt2ndEle_PASandMee100"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "2nd Electron p_{T} (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "Pt1stEle_ROI"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "1st Electron p_{T} (GeV) [Preselection + region of interest]"

plots.append ( makeDefaultPlot ( "Pt2ndEle_ROI"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "2nd Electron p_{T} (GeV) [Preselection + region of interest]"


plots.append ( makeDefaultPlot ( "Eta1stEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Electron #eta [Preselection]"   
plots[-1].rebin = 2
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = -3
plots[-1].xmax  = 3
plots[-1].ylog  = "yes"


# plots.append ( makeDefaultPlot ( "Eta1stEle_PASandMee100"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
# plots[-1].xtit = "1st Electron #eta [Preselection + M(ee) > 100]"   
# plots[-1].rebin = 2
# plots[-1].ymax = 200000000
# plots[-1].ymin = 1e-1
# plots[-1].xmin  = -3
# plots[-1].xmax  = 3
# plots[-1].ylog  = "yes"



plots.append ( makeDefaultPlot ( "Eta1stEle_ROI"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Electron #eta [Preselection + region of interest]"   
plots[-1].rebin = 2
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = -3
plots[-1].xmax  = 3
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "Eta1stEle_ROI"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Electron #eta [Preselection + region of interest]"   
plots[-1].rebin = 2
plots[-1].ymax = 30
plots[-1].ymin = 1e-1
plots[-1].xmin  = -3
plots[-1].xmax  = 3



plots.append ( makeDefaultPlot ( "Phi1stEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )   
plots[-1].xtit = "1st Electron #phi [Preselection]"
plots[-1].rebin = 4
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"


# plots.append ( makeDefaultPlot ( "Phi1stEle_PASandMee100"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )   
# plots[-1].xtit = "1st Electron #phi [Preselection + M(ee) > 100 ]"
# plots[-1].rebin = 4
# plots[-1].ymax = 20000000
# plots[-1].ymin = 1e-1
# plots[-1].ylog  = "yes"



plots.append ( makeDefaultPlot ( "Phi1stEle_ROI"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )   
plots[-1].xtit = "1st Electron #phi [Preselection + region of interest]"
plots[-1].rebin = 4
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"



plots.append ( makeDefaultPlot ( "Phi1stEle_ROI"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )   
plots[-1].xtit = "1st Electron #phi [Preselection + region of interest]"
plots[-1].rebin = 1
plots[-1].ymax = 20
plots[-1].ymin = 1e-1


plots.append ( makeDefaultPlot ( "Phi2ndEle_ROI"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )   
plots[-1].xtit = "2nd Electron #phi [Preselection + region of interest]"
plots[-1].rebin = 1
plots[-1].ymax = 20
plots[-1].ymin = 1e-1


plots.append ( makeDefaultPlot ( "Eta2ndEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "2nd Electron #eta [Preselection]"   
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = -3
plots[-1].xmax  = 3
plots[-1].rebin = 2
plots[-1].ylog  = "yes"


# plots.append ( makeDefaultPlot ( "Eta2ndEle_PASandMee100"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
# plots[-1].xtit = "2nd Electron #eta [Preselection + M(ee) > 100]"   
# plots[-1].rebin = 2
# plots[-1].ymax = 200000000
# plots[-1].ymin = 1e-1
# plots[-1].xmin  = -3
# plots[-1].xmax  = 3
# plots[-1].ylog  = "yes"



plots.append ( makeDefaultPlot ( "Eta2ndEle_ROI"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "2nd Electron #eta [Preselection + region of interest]"   
plots[-1].rebin = 2
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = -3
plots[-1].xmax  = 3
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "Eta2ndEle_ROI"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "2nd Electron #eta [Preselection + region of interest]"   
plots[-1].rebin = 2
plots[-1].ymax = 30
plots[-1].ymin = 1e-1
plots[-1].xmin  = -3
plots[-1].xmax  = 3


plots.append ( makeDefaultPlot ( "Phi2ndEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )   
plots[-1].xtit = "2nd Electron #phi [Preselection]"
plots[-1].rebin = 4
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"



# plots.append ( makeDefaultPlot ( "Phi2ndEle_PASandMee100"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )   
# plots[-1].xtit = "2nd Electron #phi [Preselection + M(ee) > 100 ]"
# plots[-1].rebin = 4
# plots[-1].ymax = 20000000
# plots[-1].ymin = 1e-1
# plots[-1].ylog  = "yes"



plots.append ( makeDefaultPlot ( "Phi2ndEle_ROI"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )   
plots[-1].xtit = "2nd Electron #phi [Preselection + region of interest]"
plots[-1].rebin = 4
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"



plots.append ( makeDefaultPlot ( "Phi2ndEle_ROI"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )   
plots[-1].xtit = "2nd Electron #phi [Preselection + region of interest]"
plots[-1].rebin = 4
plots[-1].ymax = 30
plots[-1].ymin = 1e-1



plots.append ( makeDefaultPlot ( "Charge1stEle_PAS"      ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit = "1st Electron Charge [Preselection]"
plots[-1].ymin = 0.0
plots[-1].ymax = 20000.

plots.append ( makeDefaultPlot ( "Charge2ndEle_PAS"      ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit = "2nd Electron Charge [Preselection]"
plots[-1].ymin = 0.0
plots[-1].ymax = 20000.



plots.append ( makeDefaultPlot ( "MTenu_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T}(e_{1}m PFMET [Preselection]) (GeV)"
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 500
plots[-1].xmin = 0
plots[-1].rebin = 2
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "MET_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET (GeV) [Preselection]"
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 500
plots[-1].xmin = 0
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "METPhi_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET #phi [Preselection]"
plots[-1].rebin = 4
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Pt1stJet_PAS"          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Jet p_{T} (GeV) [Preselection]"
plots[-1].rebin = pt_rebin
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmax = 1000
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Pt2ndJet_PAS"          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmax = 1000
plots[-1].xmin = 0
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet p_{T} (GeV) [Preselection]"


plots.append ( makeDefaultPlot ( "Pt1stJet_PASandMee100"          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Jet p_{T} (GeV) [Preselection + M(ee) > 100]"
plots[-1].rebin = pt_rebin
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmax = 1000
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

# plots.append ( makeDefaultPlot ( "Pt2ndJet_PASandMee100"          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
# plots[-1].rebin = pt_rebin
# plots[-1].ymax = 2000000
# plots[-1].ymin = 1e-1
# plots[-1].xmax = 1000
# plots[-1].xmin = 0
# plots[-1].ylog  = "yes"
# plots[-1].xtit = "2nd Jet p_{T} (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "Pt1stJet_ROI"          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Jet p_{T} (GeV) [Preselection + region of interest]"
plots[-1].rebin = pt_rebin
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmax = 1000
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Pt2ndJet_ROI"          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmax = 1000
plots[-1].xmin = 0
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet p_{T} (GeV) [Preselection + region of interest]"


plots.append ( makeDefaultPlot ( "Eta1stJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = -3
plots[-1].xmax  = 3
plots[-1].ylog  = "yes"
plots[-1].xtit = "1st Jet #eta [Preselection]"


# plots.append ( makeDefaultPlot ( "Eta1stJet_PASandMee100"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
# plots[-1].rebin = 2
# plots[-1].ymax = 200000000
# plots[-1].ymin = 1e-1
# plots[-1].xmin  = -3
# plots[-1].xmax  = 3
# plots[-1].ylog  = "yes"
# plots[-1].xtit = "1st Jet #eta [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "Eta1stJet_ROI"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = -3
plots[-1].xmax  = 3
plots[-1].ylog  = "yes"
plots[-1].xtit = "1st Jet #eta [Preselection + region of interest]"



plots.append ( makeDefaultPlot ( "Eta2ndJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = -3
plots[-1].xmax  = 3
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet #eta [Preselection]"


# plots.append ( makeDefaultPlot ( "Eta2ndJet_PASandMee100"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
# plots[-1].rebin = 2
# plots[-1].ymax = 200000000
# plots[-1].ymin = 1e-1
# plots[-1].xmin  = -3
# plots[-1].xmax  = 3
# plots[-1].ylog  = "yes"
# plots[-1].xtit = "2nd Jet #eta [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "Eta2ndJet_ROI"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = -3
plots[-1].xmax  = 3
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet #eta [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "Phi1stJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "1st Jet #phi [Preselection]"


# plots.append ( makeDefaultPlot ( "Phi1stJet_PASandMee100"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
# plots[-1].rebin = 1
# plots[-1].ymax = 2000000
# plots[-1].ymin = 1e-1
# plots[-1].ylog  = "yes"
# plots[-1].xtit = "1st Jet #phi [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "Phi1stJet_ROI"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "1st Jet #phi [Preselection + region of interest]"


plots.append ( makeDefaultPlot ( "Phi2ndJet_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet #phi [Preselection]"


# plots.append ( makeDefaultPlot ( "Phi2ndJet_PASandMee100"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
# plots[-1].rebin = 1
# plots[-1].ymax = 2000000
# plots[-1].ymin = 1e-1
# plots[-1].ylog  = "yes"
# plots[-1].xtit = "2nd Jet #phi [Preselection]"


plots.append ( makeDefaultPlot ( "Phi2ndJet_ROI"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet #phi [Preselection + region of interest]"



plots.append ( makeDefaultPlot ( "sTlep_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (1st Electron, 2nd Electron) (GeV) [Preselection]"


plots.append ( makeDefaultPlot ( "sTjet_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (1st Jet, 2nd Jet) (GeV) [Preselection]"


plots.append ( makeDefaultPlot ( "sTlep_PASandMee100"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (1st Electron, 2nd Electron) (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "sTjet_PASandMee100"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (1st Jet, 2nd Jet) (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "sTlep_ROI"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (1st Electron, 2nd Electron) (GeV) [Preselection + region of interest]"


plots.append ( makeDefaultPlot ( "sTjet_ROI"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (1st Jet, 2nd Jet) (GeV) [Preselection + region of interest]"


plots.append ( makeDefaultPlot ( "sT_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "sT_PASandMee100"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "sT_PASandMee110"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 110]"


plots.append ( makeDefaultPlot ( "sT_PASandMee120"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 120]"


plots.append ( makeDefaultPlot ( "sT_PASandMee130"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 130]"



plots.append ( makeDefaultPlot ( "sT_PASandMee140"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 140]"


plots.append ( makeDefaultPlot ( "sT_PASandMee150"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 150]"



plots.append ( makeDefaultPlot ( "sT_PASandMee160"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 160]"

plots.append ( makeDefaultPlot ( "sT_PASandMee170"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 170]"




plots.append ( makeDefaultPlot ( "sT_PASandMee180"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 180]"


plots.append ( makeDefaultPlot ( "sT_PASandMee190"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 190]"



plots.append ( makeDefaultPlot ( "sT_PASandMee200"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 200]"



plots.append ( makeDefaultPlot ( "sT_PASandMee100"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ylog  = "yes"
plots[-1].rebin = "var"
plots[-1].xbins = [ 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1110, 1140, 1180, 1220, 1260, 1300, 1400, 1500, 1600, 1700, 1800, 2000 ]
plots[-1].ymin = 1e-1
plots[-1].ymax = 20000
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 100, var. binning]"



plots.append ( makeDefaultPlot ( "sT_PASandMee110"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ylog  = "yes"
plots[-1].rebin = "var"
plots[-1].xbins = [ 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1110, 1140, 1180, 1220, 1260, 1300, 1400, 1500, 1600, 1700, 1800, 2000 ]
plots[-1].ymin = 1e-1
plots[-1].ymax = 20000
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 110, var. binning]"


plots.append ( makeDefaultPlot ( "sT_PASandMee120"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ylog  = "yes"
plots[-1].rebin = "var"
plots[-1].xbins = [ 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1110, 1140, 1180, 1220, 1260, 1300, 1400, 1500, 1600, 1700, 1800, 2000 ]
plots[-1].ymin = 1e-1
plots[-1].ymax = 20000
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 120, var. binning]"


plots.append ( makeDefaultPlot ( "sT_PASandMee130"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ylog  = "yes"
plots[-1].rebin = "var"
plots[-1].xbins = [ 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1110, 1140, 1180, 1220, 1260, 1300, 1400, 1500, 1600, 1700, 1800, 2000 ]
plots[-1].ymin = 1e-1
plots[-1].ymax = 20000
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 130, var. binning]"



plots.append ( makeDefaultPlot ( "sT_PASandMee140"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ylog  = "yes"
plots[-1].rebin = "var"
plots[-1].xbins = [ 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1110, 1140, 1180, 1220, 1260, 1300, 1400, 1500, 1600, 1700, 1800, 2000 ]
plots[-1].ymin = 1e-1
plots[-1].ymax = 20000
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 140, var. binning]"

plots.append ( makeDefaultPlot ( "sT_PASandMee150"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ylog  = "yes"
plots[-1].rebin = "var"
plots[-1].xbins = [ 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1110, 1140, 1180, 1220, 1260, 1300, 1400, 1500, 1600, 1700, 1800, 2000 ]
plots[-1].ymin = 1e-1
plots[-1].ymax = 20000
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 150, var. binning]"



plots.append ( makeDefaultPlot ( "sT_PASandMee160"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ylog  = "yes"
plots[-1].rebin = "var"
plots[-1].xbins = [ 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1110, 1140, 1180, 1220, 1260, 1300, 1400, 1500, 1600, 1700, 1800, 2000 ]
plots[-1].ymin = 1e-1
plots[-1].ymax = 20000
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 160, var. binning]"


plots.append ( makeDefaultPlot ( "sT_PASandMee170"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ylog  = "yes"
plots[-1].rebin = "var"
plots[-1].xbins = [ 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1110, 1140, 1180, 1220, 1260, 1300, 1400, 1500, 1600, 1700, 1800, 2000 ]
plots[-1].ymin = 1e-1
plots[-1].ymax = 20000
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 170, var. binning]"


plots.append ( makeDefaultPlot ( "sT_PASandMee180"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ylog  = "yes"
plots[-1].rebin = "var"
plots[-1].xbins = [ 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1110, 1140, 1180, 1220, 1260, 1300, 1400, 1500, 1600, 1700, 1800, 2000 ]
plots[-1].ymin = 1e-1
plots[-1].ymax = 20000
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 180, var. binning]"


plots.append ( makeDefaultPlot ( "sT_PASandMee190"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ylog  = "yes"
plots[-1].rebin = "var"
plots[-1].xbins = [ 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1110, 1140, 1180, 1220, 1260, 1300, 1400, 1500, 1600, 1700, 1800, 2000 ]
plots[-1].ymin = 1e-1
plots[-1].ymax = 20000
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 190, var. binning]"



plots.append ( makeDefaultPlot ( "sT_PASandMee200"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ylog  = "yes"
plots[-1].rebin = "var"
plots[-1].xbins = [ 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 840, 880, 920, 960, 1000, 1040, 1080, 1110, 1140, 1180, 1220, 1260, 1300, 1400, 1500, 1600, 1700, 1800, 2000 ]
plots[-1].ymin = 1e-1
plots[-1].ymax = 20000
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 200, var. binning]"




plots.append ( makeDefaultPlot ( "sT_ROI"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection + region of interest]"



plots.append ( makeDefaultPlot ( "sT_zjj_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (Z, 2jet) (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "sT_zjj_PASandMee100"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (Z, 2jet) (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "sT_zjj_ROI"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (Z, 2jet) (GeV) [Preselection + region of interest]"



plots.append ( makeDefaultPlot ( "sTfrac_Jet1_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from jet 1 (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "sTfrac_Jet1_PASandMee100"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from jet 1 (GeV) [Preselection + M(ee) > 100]"

plots.append ( makeDefaultPlot ( "sTfrac_Jet1_ROI"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from jet 1 (GeV) [Preselection + region of interest]"




plots.append ( makeDefaultPlot ( "sTfrac_Jet2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from jet 2 (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "sTfrac_Jet2_PASandMee100"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from jet 2 (GeV) [Preselection + M(ee) > 100]"

plots.append ( makeDefaultPlot ( "sTfrac_Jet2_ROI"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from jet 2 (GeV) [Preselection + region of interest]"






plots.append ( makeDefaultPlot ( "sTfrac_Ele1_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from ele 1 (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "sTfrac_Ele1_PASandMee100"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from ele 1 (GeV) [Preselection + M(ee) > 100]"

plots.append ( makeDefaultPlot ( "sTfrac_Ele1_ROI"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from ele 1 (GeV) [Preselection + region of interest]"




plots.append ( makeDefaultPlot ( "sTfrac_Ele2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from ele 2 (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "sTfrac_Ele2_PASandMee100"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from ele 2 (GeV) [Preselection + M(ee) > 100]"

plots.append ( makeDefaultPlot ( "sTfrac_Ele2_ROI"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from ele 2 (GeV) [Preselection + region of interest]"





plots.append ( makeDefaultPlot ( "sTfrac_Ele_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from electrons (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "sTfrac_Ele_PASandMee100"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from electrons (GeV) [Preselection + M(ee) > 100]"

plots.append ( makeDefaultPlot ( "sTfrac_Ele_ROI"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from electrons (GeV) [Preselection + region of interest]"



plots.append ( makeDefaultPlot ( "sTfrac_Jet_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from jets (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "sTfrac_Jet_PASandMee100"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from jets (GeV) [Preselection + M(ee) > 100]"

plots.append ( makeDefaultPlot ( "sTfrac_Jet_ROI"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.0
plots[-1].xmax = 1.0
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Fraction S_{T} from jets (GeV) [Preselection + region of interest]"





plots.append ( makeDefaultPlot ( "Mjj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Dijet Mass (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Mjj_PASandMee100"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Dijet Mass (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "Mjj_ROI"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Dijet Mass (GeV) [Preselection + region of interest]"





plots.append ( makeDefaultPlot ( "M_j2j3_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(j2,j3) Mass (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "M_j2j3_PASandMee100"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(j2,j3) Mass (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "M_j2j3_ROI"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(j2,j3) Mass (GeV) [Preselection + region of interest]"





plots.append ( makeDefaultPlot ( "M_j1j3_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(j1,j3) Mass (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "M_j1j3_PASandMee100"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(j1,j3) Mass (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "M_j1j3_ROI"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(j1,j3) Mass (GeV) [Preselection + region of interest]"


plots.append ( makeDefaultPlot ( "Me1j1_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e1,j1) Mass (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Me1j1_PASandMee100"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e1,j1) Mass (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "Me1j1_ROI"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e1,j1) Mass (GeV) [Preselection + region of interest]"




plots.append ( makeDefaultPlot ( "Me1j1_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 2000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].xtit = "M(e1,j1) Mass (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Me1j1_PASandMee100"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 500
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].xtit = "M(e1,j1) Mass (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "Me1j1_ROI"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].xtit = "M(e1,j1) Mass (GeV) [Preselection + region of interest]"





plots.append ( makeDefaultPlot ( "M_e1j3_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e1,j3) Mass (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "M_e1j3_PASandMee100"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e1,j3) Mass (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "M_e1j3_ROI"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e1,j3) Mass (GeV) [Preselection + region of interest]"




plots.append ( makeDefaultPlot ( "M_e2j3_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e2,j3) Mass (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "M_e2j3_PASandMee100"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e2,j3) Mass (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "M_e2j3_ROI"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e2,j3) Mass (GeV) [Preselection + region of interest]"


plots.append ( makeDefaultPlot ( "M_eejjj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e1, e2, j1, j2, j3) (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "M_eejjj_PASandMee100"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e1, e2, j1, j2, j3) (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "M_eejjj_ROI"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e1, e2, j1, j2, j3) (GeV) [Preselection + region of interest]"









plots.append ( makeDefaultPlot ( "Meejj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 5
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "Mass_{eejj} (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Mee_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ee) (GeV) [Preselection]"



plots.append ( makeDefaultPlot ( "Mee_ROI"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ee) (GeV) [Preselection + region of interest]"


plots.append ( makeDefaultPlot ( "Mee_PASandST445"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ee) (GeV) [Preselection + S_{T} > 445]"

plots.append ( makeDefaultPlot ( "Mee_EBEB_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 3000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) (preselection, EB-EB) (GeV)"

plots.append ( makeDefaultPlot ( "Mee_EBEE_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 600
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) (preselection, EB-EE) (GeV)"


plots.append ( makeDefaultPlot ( "Mee_EEEE_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 300
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) (preselection, EE-EE) (GeV)"

plots.append ( makeDefaultPlot ( "Mee_EB_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 3000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) (preselection, EB-EE and EB-EB) (GeV)"


plots.append ( makeDefaultPlot ( "Mee_EBEB_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 3000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [80, 100] (preselection, EB-EB) (GeV)"

plots.append ( makeDefaultPlot ( "Mee_EBEE_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 600
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [80, 100] (preselection, EB-EE) (GeV)"


plots.append ( makeDefaultPlot ( "Mee_EEEE_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 300
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [80, 100] (preselection, EE-EE) (GeV)"

plots.append ( makeDefaultPlot ( "Mee_EB_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 3000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [80, 100] (preselection, EB-EE and EB-EB) (GeV)"

plots.append ( makeDefaultPlot ( "Mee_80_100_Preselection" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].xmin = 70.
plots[-1].xmax = 110.
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ee) [80, 100] (GeV) [Preselection]"


plots.append ( makeDefaultPlot ( "Mee_EBEB_70_110_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 3000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [70, 110] (preselection, EB-EB) (GeV)"

plots.append ( makeDefaultPlot ( "Mee_EBEE_70_110_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 600
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [70, 110] (preselection, EB-EE) (GeV)"


plots.append ( makeDefaultPlot ( "Mee_EEEE_70_110_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 300
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [70, 110] (preselection, EE-EE) (GeV)"

plots.append ( makeDefaultPlot ( "Mee_EB_70_110_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 3000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [70, 110] (preselection, EB-EE and EB-EB) (GeV)"

plots.append ( makeDefaultPlot ( "Mee_70_110_Preselection" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].xmin = 70.
plots[-1].xmax = 110.
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ee) [70, 110] (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Ptee_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "P_{T}(ee) (GeV) [Preselection]"


plots.append ( makeDefaultPlot ( "Ptee_PASandMee100"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "P_{T}(ee) (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "Ptee_ROI"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "P_{T}(ee) (GeV) [Preselection + region of interest]"


plots.append ( makeDefaultPlot ( "Ptj1j2_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "P_{T}(j_{1}, j_{2}) (GeV) [Preselection]"
			     

plots.append ( makeDefaultPlot ( "Ptj2j3_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "P_{T}(jet_{2}, j_{3}) (GeV) [Preselection]"
			     

plots.append ( makeDefaultPlot ( "Ptj1j3_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "P_{T}(jet_{1}, j_{3}) (GeV) [Preselection]"
			     

plots.append ( makeDefaultPlot ( "Ptj1j2j3_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "P_{T}(jet_{1}, jet_{2}, j_{3}) (GeV) [Preselection]"
			     


plots.append ( makeDefaultPlot ( "Ptee_Minus_Ptj1j2_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = -500
plots[-1].xmax  = 500
plots[-1].ylog  = "yes"
plots[-1].xtit  = "P_{T}(ee) - P_{T}(j_{1}, j_{2}) (GeV) [Preselection]"
			     


plots.append ( makeDefaultPlot ( "Ptee_Minus_Ptj1j2_PASandMee100"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = -500
plots[-1].xmax  = 500
plots[-1].ylog  = "yes"
plots[-1].xtit  = "P_{T}(ee) - P_{T}(j_{1}, j_{2}) (GeV) [Preselection + M(ee) > 100]"
			     

plots.append ( makeDefaultPlot ( "Ptee_Minus_Ptj1j2_ROI"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = -500
plots[-1].xmax  = 500
plots[-1].ylog  = "yes"
plots[-1].xtit  = "P_{T}(ee) - P_{T}(j_{1}, j_{2}) (GeV) [region of interest]"
			     


plots.append ( makeDefaultPlot ( "Ptee_Minus_Ptj1j2j3_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = -500
plots[-1].xmax  = 500
plots[-1].ylog  = "yes"
plots[-1].xtit  = "P_{T}(ee) - P_{T}(j_{1}, j_{2}, j_{3}) (GeV) [Preselection]"


plots.append ( makeDefaultPlot ( "Ptee_Minus_Ptj1j2j3_PASandMee100"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = -500
plots[-1].xmax  = 500
plots[-1].ylog  = "yes"
plots[-1].xtit  = "P_{T}(ee) - P_{T}(j_{1}, j_{2}, j_{3}) (GeV) [Preselection + M(ee) > 100]"



plots.append ( makeDefaultPlot ( "Ptee_Minus_Ptj1j2j3_ROI"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = -500
plots[-1].xmax  = 500
plots[-1].ylog  = "yes"
plots[-1].xtit  = "P_{T}(ee) - P_{T}(j_{1}, j_{2}, j_{3}) (GeV) [region of interest]"



plots.append ( makeDefaultPlot ( "Me1j1_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 2000
plots[-1].rebin = 5
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e_{1}j_{1}) (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Me1j2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 2000
plots[-1].rebin = 5

plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e_{1}j_{2}) (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Me2j1_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 2000
plots[-1].rebin = 5

plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e_{2}j_{1}) (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Me2j2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e_{2}j_{2}) (GeV) [Preselection]"
plots[-1].xmin  = 0
plots[-1].xmax  = 2000
plots[-1].rebin = 5


plots.append ( makeDefaultPlot ( "Mej_selected_avg_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "no"
plots[-1].xtit = "M(ej) average (GeV) [Preselection]"
plots[-1].xmin  = 0
plots[-1].xmax  = 2000
plots[-1].rebin = 5

plots.append ( makeDefaultPlot ( "Mej_selected_avg_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ej) average (GeV) [Preselection]"
plots[-1].xmin  = 0
plots[-1].xmax  = 2000
plots[-1].rebin = 5


plots.append ( makeDefaultPlot ( "Mej_selected_avg_PASandMee100"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ej) average (GeV) [Preselection + M(ee) > 100]"
plots[-1].xmin  = 0
plots[-1].xmax  = 2000
plots[-1].rebin = 5


plots.append ( makeDefaultPlot ( "Mej_selected_avg_ROI"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ej) average (GeV) [Preselection + region of interest]"
plots[-1].xmin  = 0
plots[-1].xmax  = 2000
plots[-1].rebin = 5

plots.append ( makeDefaultPlot ( "nVertex_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xmin = -0.5
plots[-1].xmax = 40.5
plots[-1].ymin = 1e-1
plots[-1].ymax = 100000
plots[-1].ylog  = "yes"
plots[-1].xtit = "n(vertexes) [Preselection]"


plots.append ( makeDefaultPlot ( "nVertex_PASandMee100"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xmin = -0.5
plots[-1].xmax = 40.5
plots[-1].ymin = 1e-1
plots[-1].ymax = 100000
plots[-1].ylog  = "yes"
plots[-1].xtit = "n(vertexes) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "nVertex_ROI"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xmin = -0.5
plots[-1].xmax = 40.5
plots[-1].ymin = 1e-1
plots[-1].ymax = 100000
plots[-1].ylog  = "yes"
plots[-1].xtit = "n(vertexes) [Preselection + region of interest]"


plots.append ( makeDefaultPlot ( "nVertex_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xmin = -0.5
plots[-1].xmax = 40.5
plots[-1].ymin = 1e-1
plots[-1].ymax = 1000
plots[-1].xtit = "n(vertexes) [Preselection]"


plots.append ( makeDefaultPlot ( "nVertex_PASandMee100"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xmin = -0.5
plots[-1].xmax = 40.5
plots[-1].ymin = 1e-1
plots[-1].ymax = 300
plots[-1].xtit = "n(vertexes) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot ( "nVertex_ROI"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xmin = -0.5
plots[-1].xmax = 40.5
plots[-1].ymin = 1e-1
plots[-1].ymax = 20
plots[-1].xtit = "n(vertexes) [Preselection + region of interest]"

plots.append ( makeDefaultPlot ( "DR_Ele1Jet1_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xtit = "#DeltaR(e_{1},j_{1})"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0
plots[-1].xmax = 6
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "DR_Ele1Jet2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xtit = "#DeltaR(e_{1},j_{2})"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0
plots[-1].xmax = 6
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "DR_Ele2Jet1_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xtit = "#DeltaR(e_{2},j_{1})"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0
plots[-1].xmax = 6
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "DR_Ele2Jet2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xtit = "#DeltaR(e_{2},j_{2})"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0
plots[-1].xmax = 6
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "DR_Ele1Ele2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 4
plots[-1].xtit = "#DeltaR(e_{1},e_{2})"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0
plots[-1].xmax = 6
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "minDR_EleJet_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xtit = "Minimum #DeltaR(e_{1},(j_{1}, j_{2}, j_{3}))"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin = 0
plots[-1].xmax = 6

plots.append ( makeDefaultPlot ( "minDR_ZJet_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 4
plots[-1].xtit = "Minimum #DeltaR(ee, (j_{1}, j_{2})) [Preselection]"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin = 0
plots[-1].xmax = 6


plots.append ( makeDefaultPlot ( "minDR_ZJet_ROI"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 4
plots[-1].xtit = "Minimum #DeltaR(ee, (j_{1}, j_{2})) [Preselection + region of interest]"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin = 0
plots[-1].xmax = 6


plots.append ( makeDefaultPlot ( "DR_ZJet1_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 4
plots[-1].xtit = "#DeltaR(ee, j_{1}) [Preselection]"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin = 0
plots[-1].xmax = 6


plots.append ( makeDefaultPlot ( "DR_ZJet1_ROI"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 4
plots[-1].xtit = "#DeltaR(ee, j_{1}) [Preselection + region of interest]"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin = 0
plots[-1].xmax = 6



plots.append ( makeDefaultPlot ( "DR_ZJet2_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 4
plots[-1].xtit = "#DeltaR(ee, j_{2}) [Preselection]"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin = 0
plots[-1].xmax = 6


plots.append ( makeDefaultPlot ( "DR_ZJet2_ROI"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 4
plots[-1].xtit = "#DeltaR(ee, j_{2}) [Preselection + region of interest]"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin = 0
plots[-1].xmax = 6

plots.append ( makeDefaultPlot ( "Meejj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 10
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "Mass_{eejj} (GeV) [Preselection]"


plots.append ( makeDefaultPlot ( "Meejj_ROI"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 5
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "Mass_{eejj} (GeV) [Preselection + region of interest]"


plots.append ( makeDefaultPlot ( "Meej_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 10
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "Mass_{eej} (GeV) [Preselection]"


plots.append ( makeDefaultPlot ( "Meej_ROI"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 10
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "Mass_{eej} (GeV) [Preselection + region of interest]"

plots.append ( makeDefaultPlot ( "Mejj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 10
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "Mass_{ejj} (GeV) [Preselection]"


plots.append ( makeDefaultPlot ( "Mejj_ROI"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 10
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "Mass_{ejj} (GeV) [Preselection + region of interest]"


plots.append ( makeDefaultPlot ( "MET_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET (GeV) [Preselection]"
plots[-1].rebin = 4
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 500
plots[-1].xmin = 0
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "MET_ROI"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET (GeV) [Preselection + region of interest]"
plots[-1].rebin = 4
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 500
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot2D_NoData ( "MeeVsST_PAS" , histoBaseName2D_userDef, samplesForStackHistos, sampleForDataHisto ) )
plots[-1].xrebin = 2
plots[-1].yrebin = 2
plots[-1].xtit = "M(ee) [GeV]"
plots[-1].ytit = "S_{T}(eejj) [GeV]"
plots[-1].zlog = "yes"

plots.append ( makeDefaultPlot2D_NoData ( "MeeVsST_PAS" , histoBaseName2D_userDef, samplesForStackHistos, sampleForDataHisto ) )
plots[-1].xrebin = 2
plots[-1].yrebin = 2
plots[-1].xtit = "M(ee) [GeV]"
plots[-1].ytit = "S_{T}(eejj) [GeV]"

plots.append ( makeDefaultPlot2D_NoData ( "MeeVsST_PASandMee100" , histoBaseName2D_userDef, samplesForStackHistos, sampleForDataHisto ) )
plots[-1].xrebin = 2
plots[-1].yrebin = 2
plots[-1].xtit = "M(ee) [GeV]"
plots[-1].ytit = "S_{T}(eejj) [GeV]"
plots[-1].zlog = "yes"

plots.append ( makeDefaultPlot2D_NoData ( "MeeVsST_PASandMee100" , histoBaseName2D_userDef, samplesForStackHistos, sampleForDataHisto ) )
plots[-1].xrebin = 2
plots[-1].yrebin = 2
plots[-1].xtit = "M(ee) [GeV]"
plots[-1].ytit = "S_{T}(eejj) [GeV]"

plots.append ( makeDefaultPlot2D_NoData ( "MeeVsST_ROI" , histoBaseName2D_userDef, samplesForStackHistos, sampleForDataHisto ) )
plots[-1].xrebin = 2
plots[-1].yrebin = 2
plots[-1].xtit = "M(ee) [GeV]"
plots[-1].ytit = "S_{T}(eejj) [GeV]"
plots[-1].zlog = "yes"

plots.append ( makeDefaultPlot2D_NoData ( "MeeVsST_ROI" , histoBaseName2D_userDef, samplesForStackHistos, sampleForDataHisto ) )
plots[-1].xrebin = 2
plots[-1].yrebin = 2
plots[-1].xtit = "M(ee) [GeV]"
plots[-1].ytit = "S_{T}(eejj) [GeV]"



#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()

fileps = "allPlots_eejj_scaled_analysis.pdf"

c.Print(fileps + "[")
for i_plot, plot in enumerate(plots):
    plot.Draw(fileps, i_plot + 1)
c.Print(fileps+"]")
# os.system('ps2pdf '+fileps)
# os.system('rm '+fileps)

makeTOC ( "allPlots_eejj_scaled_analysis_toc.tex" , fileps, plots ) 

