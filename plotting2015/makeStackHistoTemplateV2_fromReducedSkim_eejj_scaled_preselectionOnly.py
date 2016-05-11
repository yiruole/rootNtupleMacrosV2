#!/usr/bin/env python

from plot_class import *
from ROOT import *

#File_preselection = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj//scaled_output_cutTable_lq_eejj_preselectionOnly/analysisClass_lq_eejj_preselectionOnly_plots.root" )
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_analysis_2Nov2015_MiniAODV2_DataSigAndSomeMoreBackgroundMC_preselOnly_withTrig//output_cutTable_lq_eejj_preselectionOnly/analysisClass_lq_eejj_preselectionOnly_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_analysis_3Nov2015_MiniAODV2_DataSigAndSomeMoreBackgroundMC_preselOnly_withTrig_withMETFilters_NoEEBadSC//output_cutTable_lq_eejj_preselectionOnly/analysisClass_lq_eejj_preselectionOnly_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_analysis_6Nov2015_1547invPb_MiniAODV2_presel_withTrig_withMETFilters_NoEEBadSC//output_cutTable_lq_eejj_2012preselectionOnly/analysisClass_lq_eejj_preselectionOnly_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_analysis_6Nov2015_1547invPb_MiniAODV2_presel_withTrig_withMETFilters_NoEEBadSC/output_cutTable_lq_eejj_preselectionOnly/analysisClass_lq_eejj_preselectionOnly_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_analysis_10Nov2015_1547invPb_MiniAODV2_presels_hltDataOnly_withMETFilters_NoEEBadSC//output_cutTable_lq_eejj_preselectionOnly_tightenEleCuts/analysisClass_lq_eejj_preselectionOnly_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_analysis_24Nov2015_AK4CHS_1547invPb_MiniAODV2_presels_hltDataOnly_withMETFilters_NoEEBadSC//output_cutTable_lq_eejj_preselectionOnly_tightenEleCuts/analysisClass_lq_eejj_preselectionOnly_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_analysis_7Dec2015_AK5_2094invPb_MiniAODV2_presels_hltDataOnly_withOldMETFilters//output_cutTable_lq_eejj_preselectionOnly/analysisClass_lq_eejj_preselectionOnly_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_29Dec2015_AK4CHS_v1-4-3_1JetOrLess/output_cutTable_lq_eejj_1jetOrLess/analysisClass_lq_eejj_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_20jan2016_v1-4-3_updateEcalCondsRun2015D/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_16Dec2015_AK4CHS_v1-4-3_Few2012LQFinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_1feb2016_v1-5-2/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_ele27wplooseData_eles35gev_noJetCuts_6feb2016_v1-5-2/output_cutTable_lq_eejj_loosenEleRequirements_noJetRequirement/analysisClass_lq_eejj_noJets_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_ele27wplooseData_eles35GeV_1Jet50GeV_7feb2016_v1-5-2/output_cutTable_lq_eejj_loosenEleRequirements_1Jet/analysisClass_lq_eejj_1Jet_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_ele27wplooseData_eles35gev_noJetCuts_6feb2016_v1-5-2/output_cutTable_lq_eejj_loosenEleRequirements_noJetRequirement/analysisClass_lq_eejj_noJets_plots.root")
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_ele27wplooseData_eles35GeV_noJets_noSt_16feb2016_v1-5-2/output_cutTable_lq_eejj_loosenEleRequirements_noJetRequirement/analysisClass_lq_eejj_noJets_plots.root")
File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_ele27wplooseData_eles45GeV_2jets_18feb2016_v1-5-2/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
# 2012
#File_preselection = GetFile("/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/analysisClass_lq_eejj_plots.root")

zUncBand="no"
makeRatio=1
makeNSigma=1

pt_rebin = 2

histoBaseName = "histo1D__SAMPLE__cutHisto_allOtherCuts___________VARIABLE"
histoBaseName_userDef = "histo1D__SAMPLE__VARIABLE"

histoBaseName2D = "histo2D__SAMPLE__cutHisto_allOtherCuts___________VARIABLE"
histoBaseName2D_userDef = "histo2D__SAMPLE__VARIABLE"

samplesForStackHistos_other = [ "OTHERBKG" ]
samplesForStackHistos_ZJets  = [ "TTbar_Madgraph", "ZJet_Madgraph_Inc" ]
#samplesForStackHistos_ZJets  = [ "TTbar_Madgraph", "ZJet_Madgraph_HT" ]
#samplesForStackHistos_ZJets  = [ "TTbar_FromData", "ZJet_Madgraph" ]
samplesForStackHistos = samplesForStackHistos_other + samplesForStackHistos_ZJets

keysStack             = [ "Other backgrounds", "t#bar{t} (Madgraph)"  ,  "Z/#gamma* + jets"  ]
stackColorIndexes     = [ 9                  , 600         ,  kRed           ]
stackFillStyleIds     = [ 3008               , 3004        ,  3345           ]

stackColorIndexes.reverse()
stackFillStyleIds.reverse()

samplesForHistos = ["LQ_M600"      ]
keys             = ["LQ, M=600 GeV"]


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
                                generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection )  ) 
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
    plot.gif_folder        = "gif_eejj_scaled_preselectionOnly/"
    plot.eps_folder        = "eps_eejj_scaled_preselectionOnly/"
    plot.pdf_folder        = "pdf_eejj_scaled_preselectionOnly/"
    plot.png_folder        = "png_eejj_scaled_preselectionOnly/"
    plot.suffix            = "eejj"
    plot.lumi_fb           = "2.1"
    #plot.lumi_fb           = "1.9"
    #plot.lumi_fb           = "2.138"
    #plot.lumi_fb           = "2.094"
    #plot.lumi_fb           = "1.547"
    #plot.lumi_fb           = "0.553"
    
    return plot

def makeDefaultPlot2D ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2D()
    plot.hasData        = True
    plot.name           = variableName
    plot.histosStack    = ( generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection ) + 
                            generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection )  ) 
    plot.histodata      =   generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_eejj_scaled_preselectionOnly/"
    plot.eps_folder     = "eps_eejj_scaled_preselectionOnly/"
    plot.suffix         = "eejj"
    
    return plot


def makeDefaultPlot2D_NoData ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2D()
    plot.hasData        = False
    plot.name           = variableName 
    plot.histosStack    = ( generateHistoList ( histoBaseName, samplesForStackHistos_other, variableName, File_preselection ) + 
                            generateHistoList ( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection )  ) 
    plot.histodata      =   generateHistoBlank( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_eejj_scaled_preselectionOnly/"
    plot.eps_folder     = "eps_eejj_scaled_preselectionOnly/"
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

#plots.append ( makeDefaultPlot ( "EleChargeSum_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].xtit  = "Electron 1 charge + Electron 2 charge [Preselection]"
#plots[-1].ymax  = 200000000
#plots[-1].ymin  = 1e-1
#plots[-1].ylog  = "yes"
#plots[-1].xmin  = -1.5
#plots[-1].xmax  = 1.5
#
#
#
#plots.append ( makeDefaultPlot ( "EleChargeSum_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].xtit  = "Electron 1 charge + Electron 2 charge [Preselection]"
#plots[-1].ymax  = 20000
#plots[-1].ymin  = 1e-1
#plots[-1].xmin  = -1.5
#plots[-1].xmax  = 1.5


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


plots.append ( makeDefaultPlot ( "Eta1stEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Electron #eta [Preselection]"   
plots[-1].rebin = 2
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = -3
plots[-1].xmax  = 3
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "Phi1stEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )   
plots[-1].xtit = "1st Electron #phi [Preselection]"
plots[-1].rebin = 4
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Eta2ndEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "2nd Electron #eta [Preselection]"   
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = -3
plots[-1].xmax  = 3
plots[-1].rebin = 2
plots[-1].ylog  = "yes"



plots.append ( makeDefaultPlot ( "Phi2ndEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )   
plots[-1].xtit = "2nd Electron #phi [Preselection]"
plots[-1].rebin = 4
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

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



plots.append ( makeDefaultPlot ( "Eta1stJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = -3
plots[-1].xmax  = 3
plots[-1].ylog  = "yes"
plots[-1].xtit = "1st Jet #eta [Preselection]"


plots.append ( makeDefaultPlot ( "Eta2ndJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = -3
plots[-1].xmax  = 3
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet #eta [Preselection]"


plots.append ( makeDefaultPlot ( "Phi1stJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "1st Jet #phi [Preselection]"

plots.append ( makeDefaultPlot ( "Phi2ndJet_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet #phi [Preselection]"



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


plots.append ( makeDefaultPlot ( "sT_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 20
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "sT_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 5
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection]"

#plots.append ( makeDefaultPlot ( "sT_zjj_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].rebin = 1
#plots[-1].ymax = 20000
#plots[-1].ymin = 1e-1
#plots[-1].xmin = 200.
#plots[-1].xmax = 2000.
#plots[-1].rebin = 2
#plots[-1].ylog  = "yes"
#plots[-1].xtit = "S_{T} (Z, 2jet) (GeV) [Preselection]"
#
#
#plots.append ( makeDefaultPlot ( "sTfrac_Jet1_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].ymax = 2000000
#plots[-1].ymin = 1e-1
#plots[-1].xmin = 0.0
#plots[-1].xmax = 1.0
#plots[-1].rebin = 4
#plots[-1].ylog  = "yes"
#plots[-1].xtit = "Fraction S_{T} from jet 1 (GeV) [Preselection]"
#
#
#plots.append ( makeDefaultPlot ( "sTfrac_Jet2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].ymax = 2000000
#plots[-1].ymin = 1e-1
#plots[-1].xmin = 0.0
#plots[-1].xmax = 1.0
#plots[-1].rebin = 4
#plots[-1].ylog  = "yes"
#plots[-1].xtit = "Fraction S_{T} from jet 2 (GeV) [Preselection]"
#
#
#
#plots.append ( makeDefaultPlot ( "sTfrac_Ele1_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].ymax = 2000000
#plots[-1].ymin = 1e-1
#plots[-1].xmin = 0.0
#plots[-1].xmax = 1.0
#plots[-1].rebin = 4
#plots[-1].ylog  = "yes"
#plots[-1].xtit = "Fraction S_{T} from ele 1 (GeV) [Preselection]"
#
#
#plots.append ( makeDefaultPlot ( "sTfrac_Ele2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].ymax = 2000000
#plots[-1].ymin = 1e-1
#plots[-1].xmin = 0.0
#plots[-1].xmax = 1.0
#plots[-1].rebin = 4
#plots[-1].ylog  = "yes"
#plots[-1].xtit = "Fraction S_{T} from ele 2 (GeV) [Preselection]"
#
#
#plots.append ( makeDefaultPlot ( "sTfrac_Ele_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].ymax = 2000000
#plots[-1].ymin = 1e-1
#plots[-1].xmin = 0.0
#plots[-1].xmax = 1.0
#plots[-1].rebin = 4
#plots[-1].ylog  = "yes"
#plots[-1].xtit = "Fraction S_{T} from electrons (GeV) [Preselection]"
#
#
#plots.append ( makeDefaultPlot ( "sTfrac_Jet_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].ymax = 2000000
#plots[-1].ymin = 1e-1
#plots[-1].xmin = 0.0
#plots[-1].xmax = 1.0
#plots[-1].rebin = 4
#plots[-1].ylog  = "yes"
#plots[-1].xtit = "Fraction S_{T} from jets (GeV) [Preselection]"


plots.append ( makeDefaultPlot ( "Mjj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Dijet Mass (GeV) [Preselection]"


#plots.append ( makeDefaultPlot ( "M_j2j3_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
#plots[-1].rebin = 1
#plots[-1].ymax = 20000
#plots[-1].ymin = 1e-1
#plots[-1].rebin = 4
#plots[-1].ylog  = "yes"
#plots[-1].xtit = "M(j2,j3) Mass (GeV) [Preselection]"
#
#plots.append ( makeDefaultPlot ( "M_j1j3_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
#plots[-1].rebin = 1
#plots[-1].ymax = 20000
#plots[-1].ymin = 1e-1
#plots[-1].rebin = 4
#plots[-1].ylog  = "yes"
#plots[-1].xtit = "M(j1,j3) Mass (GeV) [Preselection]"


plots.append ( makeDefaultPlot ( "Me1j1_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 2000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e1,j1) Mass (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Me1j1_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 200
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].xtit = "M(e1,j1) Mass (GeV) [Preselection]"

#plots.append ( makeDefaultPlot ( "M_e1j3_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
#plots[-1].rebin = 1
#plots[-1].ymax = 20000
#plots[-1].ymin = 1e-1
#plots[-1].rebin = 4
#plots[-1].ylog  = "yes"
#plots[-1].xtit = "M(e1,j3) Mass (GeV) [Preselection]"
#
#plots.append ( makeDefaultPlot ( "M_e2j3_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
#plots[-1].rebin = 1
#plots[-1].ymax = 20000
#plots[-1].ymin = 1e-1
#plots[-1].rebin = 4
#plots[-1].ylog  = "yes"
#plots[-1].xtit = "M(e2,j3) Mass (GeV) [Preselection]"
#
#plots.append ( makeDefaultPlot ( "M_eejjj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
#plots[-1].rebin = 1
#plots[-1].ymax = 20000
#plots[-1].ymin = 1e-1
#plots[-1].rebin = 4
#plots[-1].ylog  = "yes"
#plots[-1].xtit = "M(e1, e2, j1, j2, j3) (GeV) [Preselection]"


plots.append ( makeDefaultPlot ( "Meejj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 5
plots[-1].ymax = 2000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "Mass_{eejj} (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Mee_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ee) (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Mee_EBEB_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 100
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) (preselection, EB-EB) (GeV)"

plots.append ( makeDefaultPlot ( "Mee_EBEE_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) (preselection, EB-EE) (GeV)"


plots.append ( makeDefaultPlot ( "Mee_EEEE_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) (preselection, EE-EE) (GeV)"

plots.append ( makeDefaultPlot ( "Mee_EB_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 100
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) (preselection, EB-EE and EB-EB) (GeV)"


plots.append ( makeDefaultPlot ( "Mee_EBEB_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 100
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [80, 100] (preselection, EB-EB) (GeV)"

plots.append ( makeDefaultPlot ( "Mee_EBEE_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 25
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [80, 100] (preselection, EB-EE) (GeV)"


plots.append ( makeDefaultPlot ( "Mee_EEEE_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 25
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [80, 100] (preselection, EE-EE) (GeV)"

plots.append ( makeDefaultPlot ( "Mee_EB_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 100
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
plots[-1].ymax = 100
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [70, 110] (preselection, EB-EB) (GeV)"

plots.append ( makeDefaultPlot ( "Mee_EBEE_70_110_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 30
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [70, 110] (preselection, EB-EE) (GeV)"


plots.append ( makeDefaultPlot ( "Mee_EEEE_70_110_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 30
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [70, 110] (preselection, EE-EE) (GeV)"

plots.append ( makeDefaultPlot ( "Mee_EB_70_110_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 150
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [70, 110] (preselection, EB-EE and EB-EB) (GeV)"

plots.append ( makeDefaultPlot ( "Mee_70_110_Preselection" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2e4
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].xmin = 70.
plots[-1].xmax = 110.
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ee) [70, 110] (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Ptee_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2e3
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "P_{T}(ee) (GeV) [Preselection]"

plots.append ( makeDefaultPlot ( "Mee_70_110_Preselection" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 4
plots[-1].ymax = 200000000
plots[-1].ymin = 1e-1
plots[-1].rebin = 8
plots[-1].xmin = 70.
plots[-1].xmax = 110.
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ee) [70, 110] (GeV) [Preselection]"

#plots.append ( makeDefaultPlot ( "Ptj1j2_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].rebin = 2
#plots[-1].ymax  = 2000000
#plots[-1].ymin  = 1e-1
#plots[-1].xmin  = 0
#plots[-1].xmax  = 1000
#plots[-1].ylog  = "yes"
#plots[-1].xtit  = "P_{T}(j_{1}, j_{2}) (GeV) [Preselection]"
#			     
#
#plots.append ( makeDefaultPlot ( "Ptj2j3_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].rebin = 2
#plots[-1].ymax  = 2000000
#plots[-1].ymin  = 1e-1
#plots[-1].xmin  = 0
#plots[-1].xmax  = 1000
#plots[-1].ylog  = "yes"
#plots[-1].xtit  = "P_{T}(jet_{2}, j_{3}) (GeV) [Preselection]"
#			     
#
#plots.append ( makeDefaultPlot ( "Ptj1j3_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].rebin = 2
#plots[-1].ymax  = 2000000
#plots[-1].ymin  = 1e-1
#plots[-1].xmin  = 0
#plots[-1].xmax  = 1000
#plots[-1].ylog  = "yes"
#plots[-1].xtit  = "P_{T}(jet_{1}, j_{3}) (GeV) [Preselection]"
#			     
#
#plots.append ( makeDefaultPlot ( "Ptj1j2j3_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].rebin = 2
#plots[-1].ymax  = 2000000
#plots[-1].ymin  = 1e-1
#plots[-1].xmin  = 0
#plots[-1].xmax  = 1000
#plots[-1].ylog  = "yes"
#plots[-1].xtit  = "P_{T}(jet_{1}, jet_{2}, j_{3}) (GeV) [Preselection]"
#			     
#
#
#plots.append ( makeDefaultPlot ( "Ptee_Minus_Ptj1j2_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].rebin = 2
#plots[-1].ymax  = 2000000
#plots[-1].ymin  = 1e-1
#plots[-1].xmin  = -500
#plots[-1].xmax  = 500
#plots[-1].ylog  = "yes"
#plots[-1].xtit  = "P_{T}(ee) - P_{T}(j_{1}, j_{2}) (GeV) [Preselection]"
			     

#plots.append ( makeDefaultPlot ( "Ptee_Minus_Ptj1j2j3_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].rebin = 2
#plots[-1].ymax  = 2000000
#plots[-1].ymin  = 1e-1
#plots[-1].xmin  = -500
#plots[-1].xmax  = 500
#plots[-1].ylog  = "yes"
#plots[-1].xtit  = "P_{T}(ee) - P_{T}(j_{1}, j_{2}, j_{3}) (GeV) [Preselection]"


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
plots[-1].ymax = 500
plots[-1].ymin = 1e-1
plots[-1].ylog  = "no"
plots[-1].xtit = "M(ej) average (GeV) [Preselection]"
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].rebin = 5


plots.append ( makeDefaultPlot ( "Mej_selected_min_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 500
plots[-1].ymin = 1e-1
plots[-1].ylog  = "no"
plots[-1].xtit = "M(ej) minimum (GeV) [Preselection]"
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].rebin = 5


plots.append ( makeDefaultPlot ( "Mej_selected_max_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 500
plots[-1].ymin = 1e-1
plots[-1].ylog  = "no"
plots[-1].xtit = "M(ej) maximum (GeV) [Preselection]"
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].rebin = 5

plots.append ( makeDefaultPlot ( "Mej_selected_avg_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ej) average (GeV) [Preselection]"
plots[-1].xmin  = 0
plots[-1].xmax  = 1200
plots[-1].rebin = 5


plots.append ( makeDefaultPlot ( "Mej_selected_min_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ej) minimum (GeV) [Preselection]"
plots[-1].xmin  = 0
plots[-1].xmax  = 1200
plots[-1].rebin = 5


plots.append ( makeDefaultPlot ( "Mej_selected_max_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ej) maximum (GeV) [Preselection]"
plots[-1].xmin  = 0
plots[-1].xmax  = 1200
plots[-1].rebin = 5


#plots.append ( makeDefaultPlot ( "Mej_selected_diff_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].rebin = 1
#plots[-1].ymax = 2000000
#plots[-1].ymin = 1e-1
#plots[-1].ylog  = "yes"
#plots[-1].xtit = "M(ej) diff (GeV) [Preselection]"
#plots[-1].xmin  = 0
#plots[-1].xmax  = 500
#plots[-1].rebin = 5


plots.append ( makeDefaultPlot ( "nVertex_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xmin = -0.5
plots[-1].xmax = 40.5
plots[-1].ymin = 1e-1
plots[-1].ymax = 100000
plots[-1].ylog  = "yes"
plots[-1].xtit = "n(vertexes) [Preselection]"



plots.append ( makeDefaultPlot ( "nVertex_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xmin = -0.5
plots[-1].xmax = 40.5
plots[-1].ymin = 1e-1
plots[-1].ymax = 1200
plots[-1].xtit = "n(vertexes) [Preselection]"

# XXX REMOVE FOR NO JET REQUIREMENTS
##plots.append ( makeDefaultPlot ( "DR_Ele1Jet1_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
##plots[-1].rebin = 1
##plots[-1].xtit = "#DeltaR(e_{1},j_{1})"
##plots[-1].ymax = 2000000
##plots[-1].ymin = 1e-1
##plots[-1].xmin = 0
##plots[-1].xmax = 6
##plots[-1].ylog  = "yes"
##
##plots.append ( makeDefaultPlot ( "DR_Ele1Jet2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
##plots[-1].rebin = 1
##plots[-1].xtit = "#DeltaR(e_{1},j_{2})"
##plots[-1].ymax = 2000000
##plots[-1].ymin = 1e-1
##plots[-1].xmin = 0
##plots[-1].xmax = 6
##plots[-1].ylog  = "yes"
##
##plots.append ( makeDefaultPlot ( "DR_Ele2Jet1_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
##plots[-1].rebin = 1
##plots[-1].xtit = "#DeltaR(e_{2},j_{1})"
##plots[-1].ymax = 2000000
##plots[-1].ymin = 1e-1
##plots[-1].xmin = 0
##plots[-1].xmax = 6
##plots[-1].ylog  = "yes"
##
##plots.append ( makeDefaultPlot ( "DR_Ele2Jet2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
##plots[-1].rebin = 1
##plots[-1].xtit = "#DeltaR(e_{2},j_{2})"
##plots[-1].ymax = 2000000
##plots[-1].ymin = 1e-1
##plots[-1].xmin = 0
##plots[-1].xmax = 6
##plots[-1].ylog  = "yes"
##
##plots.append ( makeDefaultPlot ( "DR_Ele1Ele2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
##plots[-1].rebin = 4
##plots[-1].xtit = "#DeltaR(e_{1},e_{2})"
##plots[-1].ymax = 2000000
##plots[-1].ymin = 1e-1
##plots[-1].xmin = 0
##plots[-1].xmax = 6
##plots[-1].ylog  = "yes"
##
##plots.append ( makeDefaultPlot ( "minDR_EleJet_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
##plots[-1].rebin = 1
##plots[-1].xtit = "Minimum #DeltaR(e_{1},(j_{1}, j_{2}, j_{3}))"
##plots[-1].ymax = 2000000
##plots[-1].ymin = 1e-1
##plots[-1].ylog  = "yes"
##plots[-1].xmin = 0
##plots[-1].xmax = 6

#plots.append ( makeDefaultPlot ( "minDR_ZJet_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].rebin = 4
#plots[-1].xtit = "Minimum #DeltaR(ee, (j_{1}, j_{2})) [Preselection]"
#plots[-1].ymax = 2000000
#plots[-1].ymin = 1e-1
#plots[-1].ylog  = "yes"
#plots[-1].xmin = 0
#plots[-1].xmax = 6
#
#plots.append ( makeDefaultPlot ( "DR_ZJet1_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].rebin = 4
#plots[-1].xtit = "#DeltaR(ee, j_{1}) [Preselection]"
#plots[-1].ymax = 2000000
#plots[-1].ymin = 1e-1
#plots[-1].ylog  = "yes"
#plots[-1].xmin = 0
#plots[-1].xmax = 6
#
#plots.append ( makeDefaultPlot ( "DR_ZJet2_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
#plots[-1].rebin = 4
#plots[-1].xtit = "#DeltaR(ee, j_{2}) [Preselection]"
#plots[-1].ymax = 2000000
#plots[-1].ymin = 1e-1
#plots[-1].ylog  = "yes"
#plots[-1].xmin = 0
#plots[-1].xmax = 6


plots.append ( makeDefaultPlot ( "Meejj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 10
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "Mass_{eejj} (GeV) [Preselection]"


#plots.append ( makeDefaultPlot ( "Meej_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
#plots[-1].rebin = 10
#plots[-1].ymax = 200000
#plots[-1].ymin = 1e-1
#plots[-1].ylog  = "yes"
#plots[-1].xtit = "Mass_{eej} (GeV) [Preselection]"

# plots.append ( makeDefaultPlot ( "Mejj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
# plots[-1].rebin = 10
# plots[-1].ymax = 200000
# plots[-1].ymin = 1e-1
# plots[-1].ylog  = "yes"
# plots[-1].xtit = "Mass_{ejj} (GeV) [Preselection]"


plots.append ( makeDefaultPlot ( "MET_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET (GeV) [Preselection]"
plots[-1].rebin = 4
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 500
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

#plots.append ( makeDefaultPlot2D_NoData ( "MeeVsST_PAS" , histoBaseName2D_userDef, samplesForStackHistos, sampleForDataHisto ) )
#plots[-1].xrebin = 2
#plots[-1].yrebin = 2
#plots[-1].xtit = "M(ee) [GeV]"
#plots[-1].ytit = "S_{T}(eejj) [GeV]"
#plots[-1].zlog = "yes"
#
#plots.append ( makeDefaultPlot2D_NoData ( "MeeVsST_PAS" , histoBaseName2D_userDef, samplesForStackHistos, sampleForDataHisto ) )
#plots[-1].xrebin = 2
#plots[-1].yrebin = 2
#plots[-1].xtit = "M(ee) [GeV]"
#plots[-1].ytit = "S_{T}(eejj) [GeV]"
#
#
#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()

fileps = "allPlots_eejj_scaled_analysis_preselectionOnly.pdf"

c.Print(fileps + "[")
for i_plot, plot in enumerate(plots):
    #print 'draw plot:',plot
    plot.Draw(fileps, i_plot + 1)
c.Print(fileps+"]")

makeTOC ( "allPlots_eejj_scaled_analysis_preselectionOnly_toc.tex" , fileps, plots ) 

