#!/usr/bin/env python

from plot_class import *
from ROOT import *

File_preselection = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj//scaled_output_cutTable_lq_eejj_2012C/analysisClass_lq_eejj_plots.root" )
File_QCD          = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj_qcd//output_cutTable_lq_eejj_2012C/analysisClass_lq_eejj_QCD_plots.root")

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
    plot.ylog              = "yes"
    plot.name              = variableName
    plot.stackColorIndexes = stackColorIndexes
    plot.stackFillStyleIds = stackFillStyleIds 
    plot.gif_folder        = "gif_eejj_scaled_2012C/"
    plot.eps_folder        = "eps_eejj_scaled_2012C/"
    plot.pdf_folder        = "pdf_eejj_scaled_2012C/"
    plot.png_folder        = "png_eejj_scaled_2012C/"
    plot.suffix            = "eejj"
    plot.lumi_fb           = "6.88"
    plot.ymax  = 2000000
    plot.ymin  = 1e-1
    
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

plots.append ( makeDefaultPlot ("Classif_1stEle_PAS"               , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron classification [Preselection]"
plots[-1].ymax  = 2e9

plots.append ( makeDefaultPlot ("Classif_1stEle_PASandMee100"               , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron classification [Preselection + M(ee) > 100]"
plots[-1].ymax  = 2e9


plots.append ( makeDefaultPlot ("Classif_1stEle_ROI"               , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron classification [Preselection + region of interest]"
plots[-1].ymax  = 2e9


plots.append ( makeDefaultPlot ("Classif_2ndEle_PAS"               , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron classification [Preselection]"
plots[-1].ymax  = 2e9

plots.append ( makeDefaultPlot ("Classif_2ndEle_PASandMee100"               , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron classification [Preselection + M(ee) > 100]"
plots[-1].ymax  = 2e9


plots.append ( makeDefaultPlot ("Classif_2ndEle_ROI"               , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron classification [Preselection + region of interest]"
plots[-1].ymax  = 2e9

plots.append ( makeDefaultPlot ("CorrIsolation_1stEle_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron corrected HEEP isolation [Preselection]"
plots[-1].xmax  = 1.0
plots[-1].xmin  = -25
plots[-1].ymax  = 2e8
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("CorrIsolation_1stEle_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron corrected HEEP isolation [Preselection + M(ee) > 100]"
plots[-1].xmax  = 1.0
plots[-1].xmin  = -25
plots[-1].ymax  = 2e8
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("CorrIsolation_1stEle_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron corrected HEEP isolation [Preselection + region of interest]"
plots[-1].xmax  = 1.0
plots[-1].xmin  = -25
plots[-1].ymax  = 2e8
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("CorrIsolation_2ndEle_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron corrected HEEP isolation [Preselection]"
plots[-1].xmax  = 1.0
plots[-1].xmin  = -25
plots[-1].ymax  = 2e8
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("CorrIsolation_2ndEle_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron corrected HEEP isolation [Preselection + M(ee) > 100]"
plots[-1].xmax  = 1.0
plots[-1].xmin  = -25
plots[-1].ymax  = 2e8
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("CorrIsolation_2ndEle_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron corrected HEEP isolation [Preselection + region of interest]"
plots[-1].xmax  = 1.0
plots[-1].xmin  = -25
plots[-1].ymax  = 2e8
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("DeltaEtaTrkSC_1stEle_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron #Delta#eta(track, SC) [Preselection]"
plots[-1].xmin  = -0.007
plots[-1].xmax  =  0.007
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("DeltaEtaTrkSC_1stEle_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron #Delta#eta(track, SC) [Preselection + M(ee) > 100]"
plots[-1].xmin  = -0.007
plots[-1].xmax  =  0.007
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("DeltaEtaTrkSC_1stEle_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron #Delta#eta(track, SC) [Preselection + region of interest]"
plots[-1].xmin  = -0.007
plots[-1].xmax  =  0.007
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("DeltaEtaTrkSC_2ndEle_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron #Delta#eta(track, SC) [Preselection]"
plots[-1].xmin  = -0.007
plots[-1].xmax  =  0.007
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("DeltaEtaTrkSC_2ndEle_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron #Delta#eta(track, SC) [Preselection + M(ee) > 100]"
plots[-1].xmin  = -0.007
plots[-1].xmax  =  0.007
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("DeltaEtaTrkSC_2ndEle_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron #Delta#eta(track, SC) [Preselection + region of interest]"
plots[-1].xmin  = -0.007
plots[-1].xmax  =  0.007
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("DeltaPhiTrkSC_1stEle_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron #Delta#phi(track, SC) [Preselection]"
plots[-1].xmin  = -0.06
plots[-1].xmax  =  0.06
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("DeltaPhiTrkSC_1stEle_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron #Delta#phi(track, SC) [Preselection + M(ee) > 100]"
plots[-1].xmin  = -0.06
plots[-1].xmax  =  0.06
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("DeltaPhiTrkSC_1stEle_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron #Delta#phi(track, SC) [Preselection + region of interest]"
plots[-1].xmin  = -0.06
plots[-1].xmax  =  0.06
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("DeltaPhiTrkSC_2ndEle_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron #Delta#phi(track, SC) [Preselection]"
plots[-1].xmin  = -0.06
plots[-1].xmax  =  0.06
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("DeltaPhiTrkSC_2ndEle_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron #Delta#phi(track, SC) [Preselection + M(ee) > 100]"
plots[-1].xmin  = -0.06
plots[-1].xmax  =  0.06
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("DeltaPhiTrkSC_2ndEle_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron #Delta#phi(track, SC) [Preselection + region of interest]"
plots[-1].xmin  = -0.06
plots[-1].xmax  =  0.06
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("E1x5OverE5x5_1stEle_PAS"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron E_{1x5} / E_{5x5} [Preselection]"
plots[-1].xmin  =  0.20
plots[-1].xmax  =  1.40
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("E1x5OverE5x5_1stEle_PASandMee100"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron E_{1x5} / E_{5x5} [Preselection + M(ee) > 100]"
plots[-1].xmin  =  0.20
plots[-1].xmax  =  1.40
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("E1x5OverE5x5_1stEle_ROI"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron E_{1x5} / E_{5x5} [Preselection + region of interest]"
plots[-1].xmin  =  0.20
plots[-1].xmax  =  1.40
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("E1x5OverE5x5_2ndEle_PAS"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron E_{1x5} / E_{5x5} [Preselection]"
plots[-1].xmin  =  0.20
plots[-1].xmax  =  1.40
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("E1x5OverE5x5_2ndEle_PASandMee100"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron E_{1x5} / E_{5x5} [Preselection + M(ee) > 100]"
plots[-1].xmin  =  0.20
plots[-1].xmax  =  1.40
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("E1x5OverE5x5_2ndEle_ROI"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron E_{1x5} / E_{5x5} [Preselection + region of interest]"
plots[-1].xmin  =  0.20
plots[-1].xmax  =  1.40
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("E2x5OverE5x5_1stEle_PAS"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron E_{2x5} / E_{5x5} [Preselection]"
plots[-1].xmin  =  0.20
plots[-1].xmax  =  1.40
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("E2x5OverE5x5_1stEle_PASandMee100"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron E_{2x5} / E_{5x5} [Preselection + M(ee) > 100]"
plots[-1].xmin  =  0.20
plots[-1].xmax  =  1.40
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("E2x5OverE5x5_1stEle_ROI"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron E_{2x5} / E_{5x5} [Preselection + region of interest]"
plots[-1].xmin  =  0.20
plots[-1].xmax  =  1.40
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("E2x5OverE5x5_2ndEle_PAS"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron E_{2x5} / E_{5x5} [Preselection]"
plots[-1].xmin  =  0.20
plots[-1].xmax  =  1.40
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("E2x5OverE5x5_2ndEle_PASandMee100"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron E_{2x5} / E_{5x5} [Preselection + M(ee) > 100]"
plots[-1].xmin  =  0.20
plots[-1].xmax  =  1.40
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("E2x5OverE5x5_2ndEle_ROI"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron E_{2x5} / E_{5x5} [Preselection + region of interest]"
plots[-1].xmin  =  0.20
plots[-1].xmax  =  1.40
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("EcalIsolation_1stEle_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron HEEP ECAL isolation [Preselection]"
plots[-1].xmin  =   0.00
plots[-1].xmax  =  15.00
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("EcalIsolation_1stEle_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron HEEP ECAL isolation [Preselection + M(ee) > 100]"
plots[-1].xmin  =   0.00
plots[-1].xmax  =  15.00
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("EcalIsolation_1stEle_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron HEEP ECAL isolation [Preselection + region of interest]"
plots[-1].xmin  =   0.00
plots[-1].xmax  =  15.00
plots[-1].rebin = 4



plots.append ( makeDefaultPlot ("EcalIsolation_2ndEle_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron HEEP ECAL isolation [Preselection]"
plots[-1].xmin  =   0.00
plots[-1].xmax  =  15.00
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("EcalIsolation_2ndEle_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron HEEP ECAL isolation [Preselection + M(ee) > 100]"
plots[-1].xmin  =   0.00
plots[-1].xmax  =  15.00
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("EcalIsolation_2ndEle_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron HEEP ECAL isolation [Preselection + region of interest]"
plots[-1].xmin  =   0.00
plots[-1].xmax  =  15.00
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("HcalIsolation_1stEle_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron HEEP HCAL isolation [Preselection]"
plots[-1].xmin  =   0.00
plots[-1].xmax  =  15.00
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("HcalIsolation_1stEle_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron HEEP HCAL isolation [Preselection + M(ee) > 100]"
plots[-1].xmin  =   0.00
plots[-1].xmax  =  15.00
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("HcalIsolation_1stEle_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron HEEP HCAL isolation [Preselection + region of interest]"
plots[-1].xmin  =   0.00
plots[-1].xmax  =  15.00
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("HcalIsolation_2ndEle_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron HEEP HCAL isolation [Preselection]"
plots[-1].xmin  =   0.00
plots[-1].xmax  =  15.00
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("HcalIsolation_2ndEle_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron HEEP HCAL isolation [Preselection + M(ee) > 100]"
plots[-1].xmin  =   0.00
plots[-1].xmax  =  15.00
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("HcalIsolation_2ndEle_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron HEEP HCAL isolation [Preselection + region of interest]"
plots[-1].xmin  =   0.00
plots[-1].xmax  =  15.00
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("TrkIsolation_1stEle_PAS"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron HEEP tracker isolation [Preselection]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("TrkIsolation_1stEle_PASandMee100"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron HEEP tracker isolation [Preselection + M(ee) > 100]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("TrkIsolation_1stEle_ROI"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "1st Electron HEEP tracker isolation [Preselection + region of interest]"
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("TrkIsolation_2ndEle_PAS"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron HEEP tracker isolation [Preselection]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("TrkIsolation_2ndEle_PASandMee100"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron HEEP tracker isolation [Preselection + M(ee) > 100]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("TrkIsolation_2ndEle_ROI"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
plots[-1].xtit  = "2nd Electron HEEP tracker isolation [Preselection + region of interest]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("FBrem_1stEle_PAS"                 , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron brem fraction [Preselection]"
plots[-1].xmin  =   -5.0
plots[-1].xmax  =    4.0

plots.append ( makeDefaultPlot ("FBrem_1stEle_PASandMee100"                 , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron brem fraction [Preselection + M(ee) > 100]"
plots[-1].xmin  =   -5.0
plots[-1].xmax  =    4.0

plots.append ( makeDefaultPlot ("FBrem_1stEle_ROI"                 , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron brem fraction [Preselection + region of interest]"
plots[-1].xmin  =   -5.0
plots[-1].xmax  =    4.0


plots.append ( makeDefaultPlot ("FBrem_2ndEle_PAS"                 , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron brem fraction [Preselection]"
plots[-1].xmin  =   -5.0
plots[-1].xmax  =    4.0

plots.append ( makeDefaultPlot ("FBrem_2ndEle_PASandMee100"                 , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron brem fraction [Preselection + M(ee) > 100]"
plots[-1].xmin  =   -5.0
plots[-1].xmax  =    4.0

plots.append ( makeDefaultPlot ("FBrem_2ndEle_ROI"                 , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron brem fraction [Preselection + region of interest]"
plots[-1].xmin  =   -5.0
plots[-1].xmax  =    4.0


plots.append ( makeDefaultPlot ("GsfCtfCharge_1stEle_PAS"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron GSF CTF charge [Preselection]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10

plots.append ( makeDefaultPlot ("GsfCtfCharge_1stEle_PASandMee100"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron GSF CTF charge [Preselection + M(ee) > 100]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10

plots.append ( makeDefaultPlot ("GsfCtfCharge_1stEle_ROI"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron GSF CTF charge [Preselection + region of interest]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10


plots.append ( makeDefaultPlot ("GsfCtfCharge_2ndEle_PAS"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron GSF CTF charge [Preselection]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10

plots.append ( makeDefaultPlot ("GsfCtfCharge_2ndEle_PASandMee100"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron GSF CTF charge [Preselection + M(ee) > 100]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10

plots.append ( makeDefaultPlot ("GsfCtfCharge_2ndEle_ROI"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron GSF CTF charge [Preselection + region of interest]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10

plots.append ( makeDefaultPlot ("GsfCtfScPixCharge_1stEle_PAS"     , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron GSF CTF SC Pixel charge [Preselection]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10

plots.append ( makeDefaultPlot ("GsfCtfScPixCharge_1stEle_PASandMee100"     , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron GSF CTF SC Pixel charge [Preselection + M(ee) > 100]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10

plots.append ( makeDefaultPlot ("GsfCtfScPixCharge_1stEle_ROI"     , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron GSF CTF SC Pixel charge [Preselection + region of interest]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10


plots.append ( makeDefaultPlot ("GsfCtfScPixCharge_2ndEle_PAS"     , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron GSF CTF SC Pixel charge [Preselection]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10

plots.append ( makeDefaultPlot ("GsfCtfScPixCharge_2ndEle_PASandMee100"     , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron GSF CTF SC Pixel charge [Preselection + M(ee) > 100]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10

plots.append ( makeDefaultPlot ("GsfCtfScPixCharge_2ndEle_ROI"     , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron GSF CTF SC Pixel charge [Preselection + region of interest]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10

plots.append ( makeDefaultPlot ("GsfScPixCharge_1stEle_PAS"        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron GSF SC Pixel charge [Preselection]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10

plots.append ( makeDefaultPlot ("GsfScPixCharge_1stEle_PASandMee100"        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron GSF SC Pixel charge [Preselection + M(ee) > 100]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10

plots.append ( makeDefaultPlot ("GsfScPixCharge_1stEle_ROI"        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron GSF SC Pixel charge [Preselection + region of interest]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10


plots.append ( makeDefaultPlot ("GsfScPixCharge_2ndEle_PAS"        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron GSF SC Pixel charge [Preselection]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10

plots.append ( makeDefaultPlot ("GsfScPixCharge_2ndEle_PASandMee100"        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron GSF SC Pixel charge [Preselection + M(ee) > 100]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10

plots.append ( makeDefaultPlot ("GsfScPixCharge_2ndEle_ROI"        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron GSF SC Pixel charge [Preselection + region of interest]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e10

plots.append ( makeDefaultPlot ("HasMatchedPhot_1stEle_PAS"        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron has matched photon [Preselection]"

plots.append ( makeDefaultPlot ("HasMatchedPhot_1stEle_PASandMee100"        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron has matched photon [Preselection + M(ee) > 100]"

plots.append ( makeDefaultPlot ("HasMatchedPhot_1stEle_ROI"        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron has matched photon [Preselection + region of interest]"

plots.append ( makeDefaultPlot ("HasMatchedPhot_2ndEle_PAS"        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron has matched photon [Preselection]"

plots.append ( makeDefaultPlot ("HasMatchedPhot_2ndEle_PASandMee100"        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron has matched photon [Preselection + M(ee) > 100]"

plots.append ( makeDefaultPlot ("HasMatchedPhot_2ndEle_ROI"        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron has matched photon [Preselection + region of interest]"


plots.append ( makeDefaultPlot ("HoE_1stEle_PAS"                   , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron H/E [Preselection]"
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("HoE_1stEle_PASandMee100"                   , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron H/E [Preselection + M(ee) > 100]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("HoE_1stEle_ROI"                   , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron H/E [Preselection + region of interest]"
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("HoE_2ndEle_PAS"                   , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron H/E [Preselection]"
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("HoE_2ndEle_PASandMee100"                   , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron H/E [Preselection + M(ee) > 100]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("HoE_2ndEle_ROI"                   , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron H/E [Preselection + region of interest]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("LeadVtxDistXY_1stEle_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron D_{XY} vs leading vertex [Preselection]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("LeadVtxDistXY_1stEle_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron D_{XY} vs leading vertex [Preselection + M(ee) > 100]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("LeadVtxDistXY_1stEle_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron D_{XY} vs leading vertex [Preselection + region of interest]"
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("LeadVtxDistXY_2ndEle_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron D_{XY} vs leading vertex [Preselection]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("LeadVtxDistXY_2ndEle_PASandMee100"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron D_{XY} vs leading vertex [Preselection + M(ee) > 100]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("LeadVtxDistXY_2ndEle_ROI"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron D_{XY} vs leading vertex [Preselection + region of interest]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("LeadVtxDistZ_1stEle_PAS"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron D_{Z} vs leading vertex [Preselection]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("LeadVtxDistZ_1stEle_PASandMee100"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron D_{Z} vs leading vertex [Preselection + M(ee) > 100]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("LeadVtxDistZ_1stEle_ROI"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron D_{Z} vs leading vertex [Preselection + region of interest]"
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("LeadVtxDistZ_2ndEle_PAS"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron D_{Z} vs leading vertex [Preselection]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("LeadVtxDistZ_2ndEle_PASandMee100"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron D_{Z} vs leading vertex [Preselection + M(ee) > 100]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("LeadVtxDistZ_2ndEle_ROI"          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron D_{Z} vs leading vertex [Preselection + region of interest]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("MissingHits_1stEle_PAS"           , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron N(missing hits) [Preselection]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e7

plots.append ( makeDefaultPlot ("MissingHits_1stEle_PASandMee100"           , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron N(missing hits) [Preselection + M(ee) > 100]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e7

plots.append ( makeDefaultPlot ("MissingHits_1stEle_ROI"           , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron N(missing hits) [Preselection + region of interest]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e7


plots.append ( makeDefaultPlot ("MissingHits_2ndEle_PAS"           , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron N(missing hits) [Preselection]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e7

plots.append ( makeDefaultPlot ("MissingHits_2ndEle_PASandMee100"           , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron N(missing hits) [Preselection + M(ee) > 100]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e7

plots.append ( makeDefaultPlot ("MissingHits_2ndEle_ROI"           , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron N(missing hits) [Preselection + region of interest]"
plots[-1].ymin  = 1e-1
plots[-1].ymax  = 1e7

plots.append ( makeDefaultPlot ("NBrems_1stEle_PAS"                , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron N(brems) [Preselection]"

plots.append ( makeDefaultPlot ("NBrems_1stEle_PASandMee100"                , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron N(brems) [Preselection + M(ee) > 100]"

plots.append ( makeDefaultPlot ("NBrems_1stEle_ROI"                , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron N(brems) [Preselection + region of interest]"


plots.append ( makeDefaultPlot ("NBrems_2ndEle_PAS"                , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron N(brems) [Preselection]"

plots.append ( makeDefaultPlot ("NBrems_2ndEle_PASandMee100"                , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron N(brems) [Preselection + M(ee) > 100]"

plots.append ( makeDefaultPlot ("NBrems_2ndEle_ROI"                , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron N(brems) [Preselection + region of interest]"


plots.append ( makeDefaultPlot ("EnergyORawEnergy_1stEle_PAS"      , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron energy correction factor [Preselection]"
plots[-1].xmin  = 0.95
plots[-1].xmax  = 1.3
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("EnergyORawEnergy_1stEle_PASandMee100"      , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron energy correction factor [Preselection + M(ee) > 100]"
plots[-1].xmin  = 0.95
plots[-1].xmax  = 1.3
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("EnergyORawEnergy_1stEle_ROI"      , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron energy correction factor [Preselection + region of interest]"
plots[-1].xmin  = 0.95
plots[-1].xmax  = 1.3
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("EnergyORawEnergy_2ndEle_PAS"      , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron energy correction factor [Preselection]"
plots[-1].xmin  = 0.95
plots[-1].xmax  = 1.3
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("EnergyORawEnergy_2ndEle_PASandMee100"      , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron energy correction factor [Preselection + M(ee) > 100]"
plots[-1].xmin  = 0.95
plots[-1].xmax  = 1.3
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("EnergyORawEnergy_2ndEle_ROI"      , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron energy correction factor [Preselection + region of interest]"
plots[-1].xmin  = 0.95
plots[-1].xmax  = 1.3
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("SigmaEtaEta_Barrel_1stEle_PAS"    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron #sigma_{#eta#eta} [Preselection, Barrel]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.015
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaEtaEta_Barrel_1stEle_PASandMee100"    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron #sigma_{#eta#eta} [Preselection, Barrel + M(ee) > 100]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.015
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaEtaEta_Barrel_1stEle_ROI"    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron #sigma_{#eta#eta} [Preselection, Barrel + region of interest]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.015
plots[-1].rebin = 2


plots.append ( makeDefaultPlot ("SigmaEtaEta_Barrel_2ndEle_PAS"    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron #sigma_{#eta#eta} [Preselection, Barrel]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.015
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaEtaEta_Barrel_2ndEle_PASandMee100"    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron #sigma_{#eta#eta} [Preselection, Barrel + M(ee) > 100]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.015
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaEtaEta_Barrel_2ndEle_ROI"    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron #sigma_{#eta#eta} [Preselection, Barrel + region of interest]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.015
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaEtaEta_Endcap_1stEle_PAS"    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron #sigma_{#eta#eta} [Preselection, Endcap]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.04
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaEtaEta_Endcap_1stEle_PASandMee100"    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron #sigma_{#eta#eta} [Preselection, Endcap + M(ee) > 100]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.04
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaEtaEta_Endcap_1stEle_ROI"    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron #sigma_{#eta#eta} [Preselection, Endcap + region of interest]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.04
plots[-1].rebin = 2


plots.append ( makeDefaultPlot ("SigmaEtaEta_Endcap_2ndEle_PAS"    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron #sigma_{#eta#eta} [Preselection, Endcap]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.04
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaEtaEta_Endcap_2ndEle_PASandMee100"    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron #sigma_{#eta#eta} [Preselection, Endcap + M(ee) > 100]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.04
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaEtaEta_Endcap_2ndEle_ROI"    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron #sigma_{#eta#eta} [Preselection, Endcap + region of interest]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.04
plots[-1].rebin = 2


plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Barrel_1stEle_PAS"  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron #sigma_{i#etai#eta} [Preselection, Barrel]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.015
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Barrel_1stEle_PASandMee100"  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron #sigma_{i#etai#eta} [Preselection, Barrel + M(ee) > 100]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.015
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Barrel_1stEle_ROI"  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron #sigma_{i#etai#eta} [Preselection, Barrel + region of interest]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.015
plots[-1].rebin = 2



plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Barrel_2ndEle_PAS"  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron #sigma_{i#etai#eta} [Preselection, Barrel]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.015
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Barrel_2ndEle_PASandMee100"  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron #sigma_{i#etai#eta} [Preselection, Barrel + M(ee) > 100]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.015
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Barrel_2ndEle_ROI"  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron #sigma_{i#etai#eta} [Preselection, Barrel + region of interest]"
plots[-1].xmin  = 0.005
plots[-1].xmax  = 0.015
plots[-1].rebin = 2


plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Endcap_1stEle_PAS"  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron #sigma_{i#etai#eta} [Preselection, Endcap]"
plots[-1].xmin  = 0.015
plots[-1].xmax  = 0.035
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Endcap_1stEle_PASandMee100"  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron #sigma_{i#etai#eta} [Preselection, Endcap + M(ee) > 100]"
plots[-1].xmin  = 0.015
plots[-1].xmax  = 0.035
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Endcap_1stEle_ROI"  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron #sigma_{i#etai#eta} [Preselection, Endcap + region of interest]"
plots[-1].xmin  = 0.015
plots[-1].xmax  = 0.035
plots[-1].rebin = 2


plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Endcap_2ndEle_PAS"  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron #sigma_{i#etai#eta} [Preselection, Endcap]"
plots[-1].xmin  = 0.015
plots[-1].xmax  = 0.035
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Endcap_2ndEle_PASandMee100"  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron #sigma_{i#etai#eta} [Preselection, Endcap + M(ee) > 100]"
plots[-1].xmin  = 0.015
plots[-1].xmax  = 0.035
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Endcap_2ndEle_ROI"  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron #sigma_{i#etai#eta} [Preselection, Endcap + region of interest]"
plots[-1].xmin  = 0.015
plots[-1].xmax  = 0.035
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ("TrkPtOPt_1stEle_PAS"              , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron track pt / SC pt [Preselection]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("TrkPtOPt_1stEle_PASandMee100"              , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron track pt / SC pt [Preselection + M(ee) > 100]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("TrkPtOPt_1stEle_ROI"              , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron track pt / SC pt [Preselection + region of interest]"
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("TrkPtOPt_2ndEle_PAS"              , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron track pt / SC pt [Preselection]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("TrkPtOPt_2ndEle_PASandMee100"              , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron track pt / SC pt [Preselection + M(ee) > 100]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("TrkPtOPt_2ndEle_ROI"              , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron track pt / SC pt [Preselection + region of interest]"
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("ValidFrac_1stEle_PAS"             , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron valid fraction of hits [Preselection]"
plots[-1].xmin  = 0.00
plots[-1].xmax  = 1.50
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("ValidFrac_1stEle_PASandMee100"             , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron valid fraction of hits [Preselection + M(ee) > 100]"
plots[-1].xmin  = 0.00
plots[-1].xmax  = 1.50
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("ValidFrac_1stEle_ROI"             , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "1st Electron valid fraction of hits [Preselection + region of interest]" 
plots[-1].xmin  = 0.00
plots[-1].xmax  = 1.50
plots[-1].rebin = 4


plots.append ( makeDefaultPlot ("ValidFrac_2ndEle_PAS"             , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron valid fraction of hits [Preselection]"
plots[-1].xmin  = 0.00
plots[-1].xmax  = 1.50
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("ValidFrac_2ndEle_PASandMee100"             , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron valid fraction of hits [Preselection + M(ee) > 100]"
plots[-1].xmin  = 0.00
plots[-1].xmax  = 1.50
plots[-1].rebin = 4

plots.append ( makeDefaultPlot ("ValidFrac_2ndEle_ROI"             , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit  = "2nd Electron valid fraction of hits [Preselection + region of interest]" 
plots[-1].xmin  = 0.00
plots[-1].xmax  = 1.50
plots[-1].rebin = 4


#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()

fileps = "allPlots_eejj_2012C_scaled_electronQuality_analysis.pdf"

c.Print(fileps + "[")
for i_plot, plot in enumerate(plots):
    plot.Draw(fileps, i_plot + 1)
c.Print(fileps+"]")
# os.system('ps2pdf '+fileps)
# os.system('rm '+fileps)


makeTOC ( "allPlots_eejj_2012C_scaled_electronQuality_analysis_toc.tex" , fileps, plots ) 
