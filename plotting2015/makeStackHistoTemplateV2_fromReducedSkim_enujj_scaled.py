#!/usr/bin/env python

from plot_class import *
from ROOT import *

#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_local_nov18_addStSFplots_allDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_noIncWStitch.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_local_nov18_addStSFplots_allDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_QCD_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_local_nov18_addStSFplots_ICHEPDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_local_nov18_addStSFplots_ICHEPDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_QCD_plots.root")

File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jan20_rereco_stitch120_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jan20_rereco_stitch120_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_QCD_plots.root")
File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_QCD_jan22_rereco_enujj2012FinSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")

#### Common values for plots:
#otherBkgsKey="QCD, single top, VV+jets, W+jets"
zUncBand="no"
makeRatio=1
makeNSigma=1
doExtraPlots = False

pt_rebin   = 2
eta_rebin  = 2
st_rebin   = 2
mt_rebin   = 2
mass_rebin = 2
dphi_rebin = 2
dr_rebin   = 2

ymin = 1e-2

QCDScale = 1.0

#--- Final plots are defined here

# Simply add or remove plots and update the list plots = [plot0, plot1, ...] below

histoBaseName = "histo1D__SAMPLE__cutHisto_allOtherCuts___________VARIABLE"
histoBaseName_userDef = "histo1D__SAMPLE__VARIABLE"


#samplesForStackHistos_ZJets  = [ "TTbar_Madgraph", "ZJet_amcAtNLO_Inc" ]
#samplesForStackHistos_other = [ "WJet_amcAtNLO_Inc" , "SingleTop", "QCD_EMEnriched", "DIBOSON"]
#samplesForStackHistos_other = [ "WJet_amcAtNLO_Inc" , "SingleTop", "PhotonJets_Madgraph", "QCD_EMEnriched", "DIBOSON"]

## amc@NLO
#samplesForStackHistos_ZJets  = [ "TTbar_amcatnlo_Inc", "ZJet_amcatnlo_Inc" ]
#samplesForStackHistos_other = [ "OTHERBKG_amcAtNLOInc" ]
#keysStack             = [ "Other backgrounds", "t#bar{t} (amc@NLO)"  ,  "Z/#gamma* + jets (amc@NLO)"  ]

## MG Inc
#samplesForStackHistos_other = [ "OTHERBKG_MGInc" ]
#samplesForStackHistos_ZJets  = [ "TTbar_Madgraph", "ZJet_Madgraph_Inc" ]
#keysStack             = [ "Other backgrounds", "t#bar{t} (Madgraph)"  ,  "Z/#gamma* + jets (MG Inc)"  ]

### MG HT
#samplesForStackHistos_other = [ "OTHERBKG_MG_HT" ]
##samplesForStackHistos_ZJets  = [ "TTbar_Madgraph", "ZJet_Madgraph_HT" ]
#samplesForStackHistos_WJets  = [ "TTbar_Madgraph", "WJet_Madgraph_HT" ]
##keysStack             = [ "Other backgrounds", "t#bar{t} (Madgraph)"  ,  "Z/#gamma* + jets (MG HT)"  ]
#keysStack             = [ "Other backgrounds", "QCD multijet", "t#bar{t} (Madgraph)"  ,  "W + jets (MG HT)"  ]

## amc@NLO WJets PtBinned
#samplesForStackHistos_other = [ "OTHERBKG_ZJetWJetPt" ]
samplesForStackHistos_other = [ "OTHERBKG_amcAtNLOIncTTBar_ZJetWJetPt" ]
#samplesForStackHistos_WJets  = [ "TTbar_Madgraph", "WJet_amcatnlo_ptBinned" ]
##keysStack             = [ "Other backgrounds", "t#bar{t} (Madgraph)"  ,  "Z/#gamma* + jets (MG HT)"  ]
#keysStack             = [ "Other backgrounds", "QCD multijet", "t#bar{t} (Madgraph)"  ,  "W + jets (amc@NLO Pt)"  ]
samplesForStackHistos_WJets  = [ "TTbar_amcatnlo_Inc", "WJet_amcatnlo_ptBinned" ]
keysStack             = [ "Other backgrounds", "QCD multijet (data)", "t#bar{t} (amc@NLO)"  ,  "W + jets (amc@NLO Pt)"  ]

# QCD
samplesForStackHistos_QCD = ["QCDFakes_DATA"]
#samplesForStackHistos_QCD = ["QCD_EMEnriched"]
keysForStackHistos_QCD = ["QCD multijet (data)"]


#samplesForStackHistos_ZJets  = [ "TTbar_FromData", "ZJet_Madgraph" ]
# older
#samplesForStackHistos = samplesForStackHistos_other + samplesForStackHistos_ZJets
samplesForStackHistos = samplesForStackHistos_other + samplesForStackHistos_QCD + samplesForStackHistos_WJets
#keysStack             = [ "Other backgrounds", "QCD multijet", "t#bar{t} (Madgraph)"  ,  "W + jets (MG HT)"  ]
stackColorIndexes     = [ 9                  , kCyan         ,         600            ,  kRed           ]
stackFillStyleIds     = [ 1001               , 1001          ,  1001                  , 1001   ]

##keysStack             = [ "Other backgrounds", "t#bar{t} (Madgraph)"  ,  "Z/#gamma* + jets (MG Inc)"  ]
##stackColorIndexes     = [ 9                  , 600         ,  kRed           ]
###stackFillStyleIds     = [ 3008               , 3004        ,  3345           ]
##stackFillStyleIds     = [ 1001               , 1001        ,  1001           ]
#
##keysStack             = ["t#bar{t} (Madgraph)"     , "SingleTop"  ,"DIBOSON"   ,  "WJet(amc@NLO)",   "PhotonJets_Madgraph",  "QCD_EMEnriched",     "Z/#gamma* + jets"  ]
##stackColorIndexes     = [   600                    , kMagenta     , kGreen+3   ,  kCyan-3        ,           7            ,    kMagenta+3    ,     kRed+1           ]
##stackFillStyleIds     = [   1001                   , 1001         ,    1001    ,  1001           ,           1001         ,        1001      ,     3345           ]
##keysStack             = [ "WJet(amc@NLO)"    , "SingleTop"  ,  "QCD_EMEnriched",  "DIBOSON"   , "t#bar{t} (Madgraph)"  ,  "Z/#gamma* + jets"  ]
##stackColorIndexes     = [ 9                  , kGreen+3     ,        6         ,        3     ,    600                 ,  kRed           ]
##stackFillStyleIds     = [ 1001               , 1001         ,        1001      ,      1001    ,    3004                ,  3345           ]
#
##stackColorIndexes.reverse()
##stackFillStyleIds.reverse()

#samplesForHistos = ["LQ_M600"      ]
#keys             = ["LQ, M=600 GeV"]
samplesForHistos = []
keys             = []


samplesForHistos_blank = []
keys_blank             = []

sampleForDataHisto = "DATA"

def makeDefaultPlot ( variableName, histoBaseName, 
                      samplesForHistos, keys,
                      samplesForStackHistos, keysStack,
                      sampleForDataHisto,
                      zUncBand, makeRatio ) :
    plot                   = Plot() 
    plot.histosStack       =  ( generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection ) + 
                                generateHistoList( histoBaseName, samplesForStackHistos_QCD, variableName, File_QCD_preselection ) +
                                generateHistoList( histoBaseName, samplesForStackHistos_WJets, variableName, File_preselection )  ) 
    plot.keysStack         = keysStack
    plot.histos            = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
    plot.keys              = keys
    plot.addZUncBand       = zUncBand
    plot.makeRatio         = makeRatio
    plot.makeNSigma        = makeNSigma
    if sampleForDataHisto != '':
      plot.histodata         = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.ytit              = "N(Events)"
    plot.ylog              = "no"
    plot.name              = variableName
    plot.stackColorIndexes = stackColorIndexes
    plot.stackFillStyleIds = stackFillStyleIds 
    #plot.gif_folder        = "gif_enujj_scaled_preselectionOnly/"
    #plot.eps_folder        = "eps_enujj_scaled_preselectionOnly/"
    plot.pdf_folder        = "pdf_enujj_scaled_preselectionOnly/"
    plot.png_folder        = "png_enujj_scaled_preselectionOnly/"
    plot.suffix            = "enujj"
    #plot.lumi_fb           = "12.9"
    plot.lumi_fb           = "36.8"
    
    return plot

def makeDefaultPlot2D ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2D()
    plot.hasData        = True
    plot.name           = variableName
    plot.histosStack    = ( generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection ) + 
                            generateHistoList( histoBaseName, samplesForStackHistos_WJets, variableName, File_preselection )  ) 
    plot.histodata      =   generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_enujj_scaled_preselectionOnly/"
    plot.eps_folder     = "eps_enujj_scaled_preselectionOnly/"
    plot.suffix         = "enujj"
    
    return plot


def makeDefaultPlot2D_NoData ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2D()
    plot.hasData        = False
    plot.name           = variableName 
    plot.histosStack    = ( generateHistoList ( histoBaseName, samplesForStackHistos_other, variableName, File_preselection ) + 
                            generateHistoList ( histoBaseName, samplesForStackHistos_WJets, variableName, File_preselection )  ) 
    plot.histodata      =   generateHistoBlank( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_enujj_scaled_preselectionOnly/"
    plot.eps_folder     = "eps_enujj_scaled_preselectionOnly/"
    plot.suffix         = "enujj"
    
    return plot

def makeDefaultPlot2D_NSigma ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2DNSigma()
    plot.histosStack =   generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistos_WJets, variableName, File_preselection) 
    plot.histodata      = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_enujj_scaled_finalOnly_2012A/"
    plot.eps_folder     = "eps_enujj_scaled_finalOnly_2012A/"
    plot.suffix         = "enujj_2DNSigma_finalOnly"


    return plot

def makeDefaultPlot2D_Ratio ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2DRatio()
    plot.histosStack = generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistos_WJets, variableName, File_preselection) 
    plot.histodata      = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_enujj_scaled_finalOnly_2012A/"
    plot.eps_folder     = "eps_enujj_scaled_finalOnly_2012A/"
    plot.suffix         = "enujj_2DRatio_finalOnly"
    
    return plot

plots = []

plots.append ( makeDefaultPlot ( "GeneratorWeight", histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit            = "Generator weight"
plots[-1].ylog            = "yes"
plots[-1].rebin           = 1
plots[-1].ymin            = 0.0001
plots[-1].ymax            = 10000000000

plots.append ( makeDefaultPlot ( "PileupWeight", histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit            = "Pileup weight"
plots[-1].ylog            = "yes"
plots[-1].rebin           = 1
plots[-1].ymin            = 0.0001
plots[-1].ymax            = 10000000000


plots.append ( makeDefaultPlot ( "nJet_PAS", histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit            = "Number of jets [Preselection]"
plots[-1].ylog            = "yes"
plots[-1].rebin           = 1
plots[-1].xmin            = -0.5
plots[-1].xmax            = 10.5
plots[-1].ymin            = 0.01
plots[-1].ymax            = 10000000000

plots.append ( makeDefaultPlot ( "nElectron_PAS", histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit            = "Number of electrons [Preselection]"
plots[-1].ylog            = "yes"
plots[-1].rebin           = 1
plots[-1].xmin            = -0.5
plots[-1].xmax            = 6.5
plots[-1].ymin            = 0.0001
plots[-1].ymax            = 10000000000

plots.append ( makeDefaultPlot ( "nMuon_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "Number of muons [Preselection]"
plots[-1].ylog            = "yes"
plots[-1].rebin           = 1
plots[-1].xmin            = -0.5
plots[-1].xmax            = 6.5
plots[-1].ymin            = 0.0001
plots[-1].ymax            = 10000000000

plots.append ( makeDefaultPlot ( "Pt1stEle_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Electron p_{T} (GeV) [Preselection]"
plots[-1].xmin = 0.
plots[-1].xmax = 600.
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].rebin = pt_rebin

plots.append ( makeDefaultPlot ( "Eta1stEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Electron #eta [Preselection]"   
plots[-1].ymax = 200000000
plots[-1].rebin = eta_rebin
plots[-1].ymin = 1e-1
plots[-1].xmin = -3.
plots[-1].xmax =  3.
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Phi1stEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )   
plots[-1].xtit = "1st Electron #phi [Preselection]"
plots[-1].rebin = 1
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Charge1stEle_PAS"      ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit = "1st Electron Charge [Preselection]"
plots[-1].ylog  = "yes"
plots[-1].ymin = 1e-1
plots[-1].ymax = 10000000000

plots.append ( makeDefaultPlot ( "MET_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET (GeV) [Preselection]"
plots[-1].rebin = pt_rebin
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "METPhi_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET #phi [Preselection]"
plots[-1].rebin = 1
plots[-1].ymax = 200000000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "MET_Type01_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET0+1 (GeV) [Preselection]"
plots[-1].rebin = pt_rebin
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "MET_Type01_Phi_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET0+1 #phi [Preselection]"
plots[-1].rebin = 1
plots[-1].ymax = 200000000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "minMETPt1stEle_PAS"    ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 500
plots[-1].xmin = 0
plots[-1].ylog  = "yes"
plots[-1].xtit = "Min (PFMET, 1st Electron p_{T}) (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Pt1stJet_PAS"          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Jet p_{T} (GeV) [Preselection]"
plots[-1].rebin = pt_rebin
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Pt2ndJet_PAS"          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet p_{T} (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Eta1stJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = eta_rebin
plots[-1].ymax = 2000000000
plots[-1].ymin = 1e-1
plots[-1].xmin = -3.
plots[-1].xmax =  3.
plots[-1].ylog  = "yes"
plots[-1].xtit = "1st Jet #eta [Preselection]"

plots.append ( makeDefaultPlot ( "Eta2ndJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = eta_rebin
plots[-1].ymax = 2000000000
plots[-1].ymin = 1e-1
plots[-1].xmin = -3.
plots[-1].xmax =  3.
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet #eta [Preselection]"

plots.append ( makeDefaultPlot ( "Phi1stJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "1st Jet #phi [Preselection]"

plots.append ( makeDefaultPlot ( "Phi2ndJet_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet #phi [Preselection]"

plots.append ( makeDefaultPlot ( "CSV1stJet_PAS"        ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 5
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "1st Jet CSV [Preselection]"

plots.append ( makeDefaultPlot ( "CSV2ndJet_PAS"        ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 5
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet CSV [Preselection]"

plots.append ( makeDefaultPlot ( "MTenu_PAS"             ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (e, PFMET) (GeV) [Preselection]"
plots[-1].rebin = mt_rebin
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "MTenu_Type01_PAS"             ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (e, PFMET0+1) (GeV) [Preselection]"
plots[-1].rebin = mt_rebin
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "MTenu_50_110"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "50 < M_{T} (e, PFMET) < 110 (GeV) [Preselection] "
plots[-1].rebin = mt_rebin
plots[-1].ymax = 2e4
plots[-1].ymin = 200
plots[-1].xmax = 145
plots[-1].xmin = 40
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "MTenu_50_110_Njet_gte5"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "50 < M_{T} (e, PFMET) < 110 (GeV) [Preselection + N(Jet) #geq 5]"
plots[-1].rebin = mt_rebin
plots[-1].ymax = 3e3
plots[-1].ymin = 1
plots[-1].xmax = 145
plots[-1].xmin = 40
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "MTenu_50_110_Njet_lte4"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "50 < M_{T} (e, PFMET) < 110 (GeV) [Preselection + N(Jet) #leq 4]"
plots[-1].rebin = mt_rebin
plots[-1].ymax = 2e4
plots[-1].ymin = 200
plots[-1].xmax = 145
plots[-1].xmin = 40
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "MTenu_50_110_Njet_gte4"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "50 < M_{T} (e, PFMET) < 110 (GeV) [Preselection + N(Jet) #geq 4]"
plots[-1].rebin = mt_rebin
plots[-1].ymax = 3e3
plots[-1].ymin = 20
plots[-1].xmax = 145
plots[-1].xmin = 40
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "MTenu_50_110_Njet_lte3"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "50 < M_{T} (e, PFMET) < 110 (GeV) [Preselection + N(Jet) #leq 3]"
plots[-1].rebin = mt_rebin
plots[-1].ymax = 4e4
plots[-1].ymin = 200
plots[-1].xmax = 145
plots[-1].xmin = 40
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Ptenu_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"
plots[-1].xtit = "p_{T} (e, PFMET) (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "sTlep_PAS"             ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = st_rebin
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (e, PFMET) (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "sTjet_PAS"             ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = st_rebin
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (1st Jet, 2nd Jet) (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "sT_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = "var"
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xbins = [ 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 980, 1000, 1020, 1040, 1060, 1080, 1100, 1140, 1180, 1220, 1260, 1300, 1400, 1500, 1600, 1700, 1800, 2000 ]
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Mjj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = mass_rebin
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(jj) (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Mej1_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = mass_rebin
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 2000
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "Mass (1st Electron, 1st Jet) (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Mej2_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = mass_rebin
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 2000
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "Mass (1st Electron, 2nd Jet) (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Mej_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = mass_rebin
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 2000
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ej) (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "MTjnu_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = mass_rebin
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "M_{T}(j, PFMET) (GeV) [Preselection]"

plots.append ( makeDefaultPlot (  "DCotTheta1stEle_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Electron D #times Cot(theta) [cm] [Preselection]" 
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot (  "Dist1stEle_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Electron Distance [cm] [Preselection] " 
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot (  "mDPhi1stEleMET_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "#Delta#phi( 1st Electron, PFMET ) [Preselection]"
plots[-1].rebin = dphi_rebin
plots[-1].ymax = 5000
plots[-1].ymin = 0

plots.append ( makeDefaultPlot (  "mDPhi1stJetMET_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "#Delta#phi( 1st Jet, PFMET ) [Preselection]"
plots[-1].rebin = dphi_rebin
plots[-1].ymax = 10000
plots[-1].ymin = 0

plots.append ( makeDefaultPlot (  "mDPhi2ndJetMET_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "#Delta#phi( 2nd Jet, PFMET ) [Preselection]"
plots[-1].rebin = dphi_rebin
plots[-1].ymax = 4000
plots[-1].ymin = 0

plots.append ( makeDefaultPlot ( "minDR_EleJet_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = dr_rebin
plots[-1].xtit = "Minimum #DeltaR(e_{1},(j_{1}, j_{2}, j_{3})) [Preselection]"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "nVertex_PAS",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "n(vertex) [Preselection]"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmax = 60.5
plots[-1].xmin = -0.5
plots[-1].rebin = 1
plots[-1].ylog  = "yes"


#-----------------------------------------------------------------------------------

if doExtraPlots:
  extra_plots = []
  plots = plots + extra_plots

############# USER CODE - END ################################################
##############################################################################


#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()
fileps = "allPlots_enujj_scaled_analysis.pdf"

c.Print(fileps + "[")
for i_plot, plot in enumerate(plots):
    #print 'draw plot:',plot
    plot.Draw(fileps, i_plot + 1)
c.Print(fileps+"]")

makeTOC ( "allPlots_enujj_analysis_toc.tex" , fileps, plots ) 

