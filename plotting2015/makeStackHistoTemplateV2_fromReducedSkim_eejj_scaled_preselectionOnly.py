#!/usr/bin/env python

from plot_class import *
from ROOT import *

#File_preselection     = GetFile("$LQDATA/2016analysis/dec13_onPSK_addStSFplots_ICHEPDataExcludeEarlyRuns_rereco_ele27wplooseEta2p1Data2016CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#File_QCD_preselection = GetFile("$LQDATA/2016analysis/dec13_onPSK_addStSFplots_ICHEPDataExcludeEarlyRuns_rereco_ele27wplooseEta2p1Data2016CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_QCD_plots.root")
#File_preselection     = GetFile("$LQDATA/2016analysis/nov20_addStSFplots_allDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#File_QCD_preselection = GetFile("$LQDATA/2016analysis/nov20_addStSFplots_allDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_QCD_plots.root")

#File_preselection     = GetFile("$LQDATA/2016analysis/nov26_addStSFplots_ICHEPDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#File_QCD_preselection = GetFile("$LQDATA/2016analysis/nov26_addStSFplots_ICHEPDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_QCD_plots.root")
#File_preselection =     GetFile("$LQDATA/2016analysis/nov28_noJets_ele27wptightOrPhoton175Data2015CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj_noJets/analysisClass_lq_eejj_noJets_plots.root")
#File_QCD_preselection = GetFile("$LQDATA/2016analysis/nov28_noJets_ele27wptightOrPhoton175Data2015CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj_noJets/analysisClass_lq_eejj_noJets_QCD_hack_plots.root")
#File_preselection =     GetFile("$LQDATA/2016analysis/nov26_addStSFplots_ICHEPDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#File_QCD_preselection = GetFile("$LQDATA/2016analysis/nov26_addStSFplots_ICHEPDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_QCD_plots.root")
#File_preselection     = GetFile("$LQDATA/2016analysis/nov28_onRSK_addStSFplots_ICHEPDataExcludeEarlyRuns_ele27wplooseEta2p1Data2016CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#File_QCD_preselection = GetFile("$LQDATA/2016analysis/nov28_onRSK_addStSFplots_ICHEPDataExcludeEarlyRuns_ele27wplooseEta2p1Data2016CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_QCD_plots.root")


#File_preselection     = GetFile("$LQDATA/2016analysis/jan20_onPSK_rereco_DYWStitch120GeV_ele27wptightEta2p1Data2016CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#File_preselection     = GetFile("$LQDATA/2016analysis/jan20_onPSK_rereco_DYWStitch120GeV_ele27wptightEta2p1Data2016CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_scaled.root")
##File_QCD_preselection = GetFile("$LQDATA/2016analysis/eejj_psk_QCD_jan22_rereco_eejj2015FinSels//output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
#File_QCD_preselection = GetFile("$LQDATA/2016analysis/eejj_psk_QCD_jan24_rereco_eejj2015FinSels//output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")

##File_QCD_preselection = GetFile("$LQDATA/2016analysis/feb22_onPsk_QCD_jan24_rereco/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
#File_QCD_preselection = GetFile("$LQDATA/2016qcd/eejj_psk_mar26_recoHeepSFs_reMiniAOD_ele27wptightEta2p1Data2016CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
## unscaled
##File_preselection     = GetFile("$LQDATA/2016analysis/eejj_feb28_recoHeepSFs_onPSK_rereco_ele27wptightEta2p1Data2016CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root")
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_mar28_recoHeepSFs_onPSK_reminiAOD_ele27wptightEta2p1Data2016CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
## scaled
##File_preselection     = GetFile("$LQDATA/2016analysis/eejj_feb28_recoHeepSFs_onPSK_rereco_ele27wptightEta2p1Data2016CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_noTTBarDataDriven.root")
## scaled+TTBarData
##File_preselection     = GetFile("$LQDATA/2016analysis/eejj_feb28_recoHeepSFs_onPSK_rereco_ele27wptightEta2p1Data2016CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")

#File_QCD_preselection = GetFile("$LQDATA/2016qcd/eejj_psk_mar26_recoHeepSFs_reMiniAOD_ele27wptightEta2p1Data2016CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
#File_QCD_preselection = GetFile("$LQDATA/2016analysis/eejj_QCD_psk_may22_ele27wptightOREle115_eejjOptFinalSels//output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
#File_QCD_preselection = GetFile("$LQDATA/2016analysis/eejj_QCD_psk_may29_ele27wptightOREle115ORPhoton175_eejjOptFinalSels/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
#File_QCD_preselection = GetFile("$LQDATA/2016analysis/eejj_QCD_psk_jul2_ele27wptightOREle115ORPhoton175_eejjOptFinalSels/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
#File_QCD_preselection = GetFile("$LQDATA/2016analysis/eejj_QCD_psk_sep29_ptEECut/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
#File_QCD_preselection = GetFile("$LQDATA/2016qcd/eejj_QCD_psk_oct6_ptEECut_updateFinalSels/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
#File_QCD_preselection = GetFile("$LQDATA/2016qcd/eejj_psk_oct6_ptEECut_actualUpdateFinalSels/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
#
#File_QCD_preselection = GetFile("$LQDATA/2016qcd/eejj_psk_oct28_newFR_ptEECut_noMuonReq/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
#File_QCD_preselection = GetFile("$LQDATA/2016qcd/eejj_psk_nov19_finalSels_noMuonVeto_nEleGte2/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
#
#File_QCD_preselection = GetFile("$LQDATA/2016qcd/eejj_psk_nov27_finalSels_muonVeto35GeV_nEleGte2/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
#
#File_QCD_preselection = GetFile("$LQDATA/2016qcd/eejj_QCD_jan26_gsfEtaCheck_finalSels/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
#
#File_QCD_preselection = GetFile("$LQDATA/2016qcd/eejj_QCD_feb10_bugfix/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
#
#File_QCD_preselection = GetFile("$LQDATA/2016qcd/eejj_QCD_feb13_addPlot/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
#
#File_QCD_preselection = GetFile("$LQDATA/2016qcd/eejj_QCD_mar16_fixMuons/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")
File_QCD_preselection = GetFile("$LQDATA/2016qcd/eejj_QCD_mar20_fixPlots/output_cutTable_lq_eejj_QCD/analysisClass_lq_eejj_QCD_plots.root")

# unscaled
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_mar30_topPtWeight_recoHeepSFs_ele27wptightEta2p1Data2016CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root")
# scaled
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_mar30_topPtWeight_recoHeepSFs_ele27wptightEta2p1Data2016CurveMC_eejj2015FinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
# ele27||ele115
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_apr11_ele27wptightOREle115_eejjOptFinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root")
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_apr11_ele27wptightOREle115_eejjOptFinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
# unscaled
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_may8_lowWZPt_ele27wptightOREle115_eejjBadOptFinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root")
# scaled
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_may8_lowWZPt_ele27wptightOREle115_eejjBadOptFinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
# unscaled
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_may22_lowWZPt_ele27wptightOREle115_eejjBadOptFinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root")
# scaled
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_may22_lowWZPt_ele27wptightOREle115_eejjBadOptFinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_may30_properEle27wptightOREle115ORPhoton175_eejjBadOptFinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root")
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_may30_properEle27wptightOREle115ORPhoton175_eejjBadOptFinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_jul4_properEle27wptightOREle115ORPhoton175_eejjOptFinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root")
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_jul4_properEle27wptightOREle115ORPhoton175_eejjOptFinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_sep17_ptEECut_properEle27wptightOREle115ORPhoton175_eejjOptFinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root")
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_sep17_ptEECut_properEle27wptightOREle115ORPhoton175_eejjOptFinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_oct6_ptEECut_updateFinalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_oct28_ptEECut_noMuonReq/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_nov19_updateTrigEff_finalSels_noMuonVeto_nEleGte2/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_nov24_fixTrigEff_finalSels_muonVetoDef35GeV_nEleGte2/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root")
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_nov24_fixTrigEff_finalSels_muonVetoDef35GeV_nEleGte2/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_jan26_gsfEtaCheck_finalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root")
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_jan26_gsfEtaCheck_finalSels/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_feb10_bugfix/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root")
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_feb10_bugfix/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_feb20_newSingTop/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root")
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_feb20_newSingTop/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
#
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_mar16_fixMuons/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots_unscaled.root")
#File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_mar16_fixMuons/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")
File_preselection      = GetFile("$LQDATA/2016analysis/eejj_psk_mar20_fixPlots/output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root")


#File_ttbar_preselection = GetFile("$LQDATA/2016ttbar/nov19_emujj/output_cutTable_lq_ttbar_emujj_correctTrig/analysisClass_lq_ttbarEst_plots.root")
#File_ttbar_preselection = GetFile("$LQDATA/2016ttbar/feb2_newSkim_emujj_correctTrig_finalSelections/output_cutTable_lq_ttbar_emujj_correctTrig/analysisClass_lq_ttbarEst_plots.root")
#File_ttbar_preselection = GetFile("$LQDATA/2016ttbar/feb11_emujj_correctTrig/output_cutTable_lq_ttbar_emujj_correctTrig/analysisClass_lq_ttbarEst_plots.root")
#File_ttbar_preselection = GetFile("$LQDATA/2016ttbar/feb28_emujj_RTrigBugFix_correctTrig/output_cutTable_lq_ttbar_emujj_correctTrig/analysisClass_lq_ttbarEst_plots.root")
# fixed ttbar
#File_ttbar_preselection = GetFile("$LQDATA/2016ttbar/mar1_emujj_RedoRTrig/output_cutTable_lq_ttbar_emujj_correctTrig/analysisClass_lq_ttbarEst_plots.root")
#
#File_ttbar_preselection = GetFile("$LQDATA/2016ttbar/mar17_emujj_fixMuons/output_cutTable_lq_ttbar_emujj_correctTrig/analysisClass_lq_ttbarEst_plots.root")
File_ttbar_preselection = GetFile("$LQDATA/2016ttbar/mar20_emujj_fixPlots/output_cutTable_lq_ttbar_emujj_correctTrig/analysisClass_lq_ttbarEst_plots.root")


LQmasses = [650,1500]
LQmassesFinalSelection = [400,650,1000]
doPreselPlots = True
doFinalSelectionPlots = True

zUncBand="no"
bkgUncBand=True
## FIXME parse these from dat files or hists
#preselAllBkg = 50632.42
#backgroundSyst = GetBackgroundSyst(preselAllBkg,42597.24,6262.92,36.85)
#backgroundSyst /= preselAllBkg
makeRatio=1
makeNSigma=1


pt_rebin = 2

histoBaseName = "histo1D__SAMPLE__cutHisto_allOtherCuts___________VARIABLE"
histoBaseName_userDef = "histo1D__SAMPLE__VARIABLE"

histoBaseName2D = "histo2D__SAMPLE__cutHisto_allOtherCuts___________VARIABLE"
histoBaseName2D_userDef = "histo2D__SAMPLE__VARIABLE"

#samplesForStackHistos_ZJets  = [ "TTbar_Madgraph", "ZJet_amcAtNLO_Inc" ]
#samplesForStackHistos_other = [ "WJet_amcAtNLO_Inc" , "SingleTop", "QCD_EMEnriched", "DIBOSON"]
#samplesForStackHistos_other = [ "WJet_amcAtNLO_Inc" , "SingleTop", "PhotonJets_Madgraph", "QCD_EMEnriched", "DIBOSON"]

## amc@NLO
#samplesForStackHistos_ZJets  = [ "TTbar_amcatnlo_Inc", "ZJet_amcatnlo_Inc" ]
#samplesForStackHistos_other = [ "OTHERBKG_amcAtNLOInc" ]
#keysStack             = [ "Other backgrounds", "t#bar{t} (MG5_aMC)"  ,  "Z/#gamma* + jets (MG5_aMC)"  ]

## MG Inc
#samplesForStackHistos_other = [ "OTHERBKG_MGInc" ]
#samplesForStackHistos_ZJets  = [ "TTbar_Madgraph", "ZJet_Madgraph_Inc" ]
#keysStack             = [ "Other backgrounds", "t#bar{t} (Madgraph)"  ,  "Z/#gamma* + jets (MG Inc)"  ]

useMGHT=False
if useMGHT:
    # MG HT
    samplesForStackHistos_other = [ "OTHERBKG_MG_HT" ]
    samplesForStackHistos_ZJets  = [ "TTbar_Madgraph", "ZJet_Madgraph_HT" ]
    #keysStack             = [ "Other backgrounds", "t#bar{t} (Madgraph)"  ,  "Z/#gamma* + jets (MG HT)"  ]
    keysStack             = [ "Other backgrounds", "QCD multijet", "t#bar{t} (Madgraph)"  ,  "Z/#gamma* + jets (MG HT)"  ]

## amc@NLO Pt Z, TTBar MG
#samplesForStackHistos_ZJets  = [ "TTbar_Madgraph", "ZJet_amcatnlo_ptBinned" ]
#samplesForStackHistos_other = [ "OTHERBKG_MG_ZJetPt" ]
##keysStack             = [ "Other backgrounds", "t#bar{t} (MG)"  ,  "Z/#gamma* + jets (MG5_aMC Pt)"  ]
#keysStack             = [ "Other backgrounds", "QCD multijet", "t#bar{t} (Madgraph)"  ,  "Z/#gamma* + jets (MG5_aMC Pt)"  ]

else:
    # amc@NLO Pt ZJets and TTBar
    #samplesForStackHistos_other = [ "OTHERBKG_WJetPt" ]
    samplesForStackHistos_other = [ "OTHERBKG_WJetPt_amcAtNLODiboson" ]
    # MC ttbar
    #samplesForStackHistos_ttbar = [ "TTbar_amcatnlo_Inc" ]
    #keysStack             = ["QCD multijet (data)", "Other backgrounds", "t#bar{t} (MG5_aMC)"  ,  "Z/#gamma* + jets (MG5_aMC Pt)"  ]
    # data-driven ttbar
    samplesForStackHistos_ttbar  = [ "TTBarFromDATA" ]
    keysStack             = ["QCD multijet (data)", "Other backgrounds", "t#bar{t} (data)", "Z/#gamma* + jets (MG5_aMC Pt)" ,   ]
    samplesForStackHistos_ZJets  = [ "ZJet_amcatnlo_ptBinned" ]
    systTypes             = ['qcd', 'mc', 'ttbarfromdata', 'zjets']

# QCD
samplesForStackHistos_QCD = ["QCDFakes_DATA"]
#samplesForStackHistos_QCD = ["QCD_EMEnriched"]
#keysForStackHistos_QCD = ["QCD multijet (data)"]


#samplesForStackHistos_ZJets  = [ "TTbar_FromData", "ZJet_Madgraph" ]
# older
#samplesForStackHistos = samplesForStackHistos_other + samplesForStackHistos_ZJets
samplesForStackHistos =  samplesForStackHistos_QCD + samplesForStackHistos_other + samplesForStackHistos_ttbar + samplesForStackHistos_ZJets
#print 'samplesForStackHistos',samplesForStackHistos
#keysStack             = [ "Other backgrounds", "QCD multijet", "t#bar{t} (Madgraph)"  ,  "Z/#gamma* + jets (MG HT)"  ]
#keysStack             = [ "Other backgrounds", "QCD multijet", "t#bar{t} (Madgraph)"  ,  "Z/#gamma* + jets (amc@NLO Pt)"  ]
stackColorIndexes     = [kCyan         , 9                  ,    600                  ,  kRed  ]
stackFillStyleIds     = [1001          , 1001               ,   1001                  , 1001   ]

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

samplesForHistos = ["LQ_M{}".format(lqmass) for lqmass in LQmasses      ]
keys             = ["LQ, M={} GeV".format(lqmass) for lqmass in LQmasses ]
# no signal
#samplesForHistos = []
#keys             = []


samplesForHistos_blank = []
keys_blank             = []

sampleForDataHisto = "DATA"
#dataBlindAbovePt1 = 800 # GeV; used for ele Pt1, Mee, Mej
#dataBlindAbovePt2 = 400 # for ele pt2
#dataBlindAboveSt = 1500 # for St plots
dataBlindAboveSt = -1 # for St plots
dataBlindAboveAll = -1
dataBlindAbovePt1 = -1 # GeV; used for ele Pt1, Mee, Mej
dataBlindAbovePt2 = -1 # for ele pt2


def makeDefaultPlot ( variableName, histoBaseName, 
                      samplesForHistos, keys,
                      samplesForStackHistos, keysStack,
                      sampleForDataHisto,
                      zUncBand, makeRatio, dataBlindAbove = -1):
    plot                   = Plot() 
    plot.histosStack       =  ( generateHistoList( histoBaseName, samplesForStackHistos_QCD, variableName, File_QCD_preselection ) +
                                generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection ) + 
                                generateHistoList( histoBaseName, samplesForStackHistos_ttbar, variableName, File_ttbar_preselection ) + 
                                generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection )  ) 
    plot.keysStack         = keysStack
    plot.systs             = [GetBackgroundSyst(systType,True) for systType in systTypes]
    plot.histos            = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
    plot.keys              = keys
    plot.addZUncBand       = zUncBand
    plot.makeRatio         = makeRatio
    plot.makeNSigma        = makeNSigma
    if sampleForDataHisto != '':
      scale = 1.0
      plot.histodata         = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection, scale, dataBlindAbove)
      plot.histodataBlindAbove = dataBlindAbove
    plot.ytit              = "N(Events)"
    plot.ylog              = "no"
    plot.name              = variableName
    plot.stackColorIndexes = stackColorIndexes
    plot.stackFillStyleIds = stackFillStyleIds 
    #plot.gif_folder        = "gif_eejj_scaled_preselectionOnly/"
    #plot.eps_folder        = "eps_eejj_scaled_preselectionOnly/"
    plot.pdf_folder        = "pdf_eejj_scaled_preselectionOnly/"
    plot.png_folder        = "png_eejj_scaled_preselectionOnly/"
    plot.root_folder          = "root_eejj_scaled_preselectionOnly/"
    plot.suffix            = "eejj"
    plot.lumi_fb           = "35.9"
    plot.addBkgUncBand     = bkgUncBand
    
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

def makeDefaultPlot2D_NSigma ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2DNSigma()
    plot.histosStack =   generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection) 
    plot.histodata      = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_eejj_scaled_finalOnly_2012A/"
    plot.eps_folder     = "eps_eejj_scaled_finalOnly_2012A/"
    plot.suffix         = "eejj_2DNSigma_finalOnly"


    return plot

def makeDefaultPlot2D_Ratio ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2DRatio()
    plot.histosStack = generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection) 
    plot.histodata      = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_eejj_scaled_finalOnly_2012A/"
    plot.eps_folder     = "eps_eejj_scaled_finalOnly_2012A/"
    plot.suffix         = "eejj_2DRatio_finalOnly"
    
    return plot


plots = []

####################################################################################################
#  PRESELECTION
####################################################################################################
if doPreselPlots:
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
    plots[-1].ymax  = 5e4
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
    
    
    plots.append ( makeDefaultPlot ( "PtHeep1stEle_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    plots[-1].rebin = 2
    plots[-1].ymax  = 1e4
    plots[-1].ymin  = 1e-1
    plots[-1].xmin  = 0
    plots[-1].xmax  = 1500
    plots[-1].ylog  = "yes"
    plots[-1].xtit  = "1st Electron p_{T} (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "PtHeep2ndEle_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt2) )
    plots[-1].rebin = 2
    plots[-1].ymax  = 1e5
    plots[-1].ymin  = 1e-1
    plots[-1].xmin  = 0
    plots[-1].xmax  = 1000
    plots[-1].ylog  = "yes"
    plots[-1].xtit  = "2nd Electron p_{T} (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "Pt1stMuon_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    #plots[-1].rebin = 2
    plots[-1].ymax  = 1e4
    plots[-1].ymin  = 1e-1
    plots[-1].xmin  = 0
    plots[-1].xmax  = 500
    plots[-1].ylog  = "yes"
    plots[-1].xtit  = "1st Muon p_{T} (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "Pt2ndMuon_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt2) )
    #plots[-1].rebin = 2
    plots[-1].ymax  = 1e5
    plots[-1].ymin  = 1e-1
    plots[-1].xmin  = 0
    plots[-1].xmax  = 200
    plots[-1].ylog  = "yes"
    plots[-1].xtit  = "2nd Muon p_{T} (GeV) [Preselection]"
    
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
    plots[-1].ymax = 5e7
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
    plots[-1].ymax = 5e7
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "Charge1stEle_PAS"      ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
    plots[-1].xtit = "1st Electron Charge [Preselection]"
    plots[-1].ymin = 0.0
    plots[-1].ymax = 8e4
    
    plots.append ( makeDefaultPlot ( "Charge2ndEle_PAS"      ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
    plots[-1].xtit = "2nd Electron Charge [Preselection]"
    plots[-1].ymin = 0.0
    plots[-1].ymax = 8e4
    
    
    
    plots.append ( makeDefaultPlot ( "MTenu_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "M_{T}(e_{1}m PFMET [Preselection]) (GeV)"
    plots[-1].rebin = 1
    plots[-1].ymax = 300000
    plots[-1].ymin = 1e-1
    plots[-1].xmax = 500
    plots[-1].xmin = 0
    plots[-1].rebin = 2
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MET_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "PFMET (GeV) [Preselection]"
    plots[-1].rebin = 1
    plots[-1].ymax = 1e4
    plots[-1].ymin = 1e-1
    plots[-1].xmax = 800
    plots[-1].xmin = 0
    plots[-1].ylog  = "yes"
    
    
    plots.append ( makeDefaultPlot ( "METPhi_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "PFMET #phi [Preselection]"
    plots[-1].rebin = 4
    plots[-1].ymax = 5e7
    plots[-1].ymin = 1e-1
    #plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "Pt1stJet_PAS"          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "1st Jet p_{T} (GeV) [Preselection]"
    plots[-1].rebin = pt_rebin
    plots[-1].ymax = 1e4
    plots[-1].ymin = 1e-1
    plots[-1].xmax = 2000
    plots[-1].xmin = 0
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "Pt2ndJet_PAS"          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = pt_rebin
    plots[-1].ymax = 1e4
    plots[-1].ymin = 1e-1
    plots[-1].xmax = 2000
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
    plots[-1].ymax = 4e7
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "1st Jet #phi [Preselection]"
    
    plots.append ( makeDefaultPlot ( "Phi2ndJet_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
    plots[-1].rebin = 1
    plots[-1].ymax = 4e7
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "2nd Jet #phi [Preselection]"
    
    
    
    plots.append ( makeDefaultPlot ( "sTlep_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 1e6
    plots[-1].ymin = 1e-1
    plots[-1].xmax  = 3000
    plots[-1].xmin  = 200
    plots[-1].rebin = 4
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "S_{T} (1st Electron, 2nd Electron) (GeV) [Preselection]"
    
    
    plots.append ( makeDefaultPlot ( "sTjet_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 1e6
    plots[-1].ymin = 1e-1
    plots[-1].xmin  = 200
    plots[-1].xmax  = 3000
    plots[-1].rebin = 4
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "S_{T} (1st Jet, 2nd Jet) (GeV) [Preselection]"
    
    
    plots.append ( makeDefaultPlot ( "sT_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAboveSt) )
    plots[-1].rebin = 1
    plots[-1].ymax = 100000
    plots[-1].ymin = 1e-1
    plots[-1].xmin = 200.
    plots[-1].xmax = 3000.
    plots[-1].rebin = 20
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "S_{T} (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "sT_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAboveSt) )
    plots[-1].rebin = 1
    plots[-1].ymax = 1e4
    plots[-1].ymin = 1e-1
    plots[-1].xmin = 200.
    plots[-1].xmax = 3000.
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
    plots[-1].ymax = 1e6
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
    plots[-1].ymax = 1e5
    plots[-1].ymin = 1e-1
    plots[-1].rebin = 4
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "M(e1,j1) Mass (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "Me1j1_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
    plots[-1].rebin = 1
    plots[-1].ymax = 1e5
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
    
    
    plots.append ( makeDefaultPlot ( "Meejj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAboveSt) )    
    plots[-1].rebin = 5
    plots[-1].ymax = 6e5
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "Mass_{eejj} (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "Mee_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].rebin = 1
    plots[-1].rebin = "var"
    plots[-1].xbins = range(0,410,10)+[420,440,460,480,500,550,600,650,700,800,900,1000,1200,1400,1600,1800,2000]
    plots[-1].ymax = 2e4
    plots[-1].ymin = 1e-1
    plots[-1].xmin = 0.
    plots[-1].xmax = 2000.
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "M(ee) (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "Mee_EBEB_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 1e4
    plots[-1].ymin = 1e-1
    plots[-1].xmin = 0.
    plots[-1].xmax = 1000.
    plots[-1].xtit = "M(ee) (preselection, EB-EB) (GeV)"
    
    plots.append ( makeDefaultPlot ( "Mee_EBEE_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 2e3
    plots[-1].ymin = 1e-1
    plots[-1].xmin = 0.
    plots[-1].xmax = 1000.
    plots[-1].xtit = "M(ee) (preselection, EB-EE) (GeV)"
    
    
    plots.append ( makeDefaultPlot ( "Mee_EEEE_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 1e3
    plots[-1].ymin = 1e-1
    plots[-1].xmin = 0.
    plots[-1].xmax = 1000.
    plots[-1].xtit = "M(ee) (preselection, EE-EE) (GeV)"
    
    plots.append ( makeDefaultPlot ( "Mee_EB_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 8e3
    plots[-1].ymin = 1e-1
    plots[-1].xmin = 0.
    plots[-1].xmax = 1000.
    plots[-1].xtit = "M(ee) (preselection, EB-EE and EB-EB) (GeV)"
    
    
    plots.append ( makeDefaultPlot ( "Mee_EBEB_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 8e3
    plots[-1].ymin = 1e-1
    plots[-1].xmin = 0.
    plots[-1].xmax = 1000.
    plots[-1].xtit = "M(ee) [80, 100] (preselection, EB-EB) (GeV)"
    
    plots.append ( makeDefaultPlot ( "Mee_EBEE_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 1e3
    plots[-1].ymin = 1e-1
    plots[-1].xmin = 0.
    plots[-1].xmax = 1000.
    plots[-1].xtit = "M(ee) [80, 100] (preselection, EB-EE) (GeV)"
    
    
    plots.append ( makeDefaultPlot ( "Mee_EEEE_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 1e3
    plots[-1].ymin = 1e-1
    plots[-1].xmin = 0.
    plots[-1].xmax = 1000.
    plots[-1].xtit = "M(ee) [80, 100] (preselection, EE-EE) (GeV)"
    
    plots.append ( makeDefaultPlot ( "Mee_EB_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 6e3
    plots[-1].ymin = 1e-1
    plots[-1].xmin = 0.
    plots[-1].xmax = 1000.
    plots[-1].xtit = "M(ee) [80, 100] (preselection, EB-EE and EB-EB) (GeV)"
    
    plots.append ( makeDefaultPlot ( "Mee_80_100_Preselection" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 1e8
    plots[-1].ymin = 1e-1
    plots[-1].rebin = 4
    plots[-1].xmin = 70.
    plots[-1].xmax = 110.
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "M(ee) [80, 100] (GeV) [Preselection]"
    
    
    plots.append ( makeDefaultPlot ( "Mee_EBEB_70_110_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 1e4
    plots[-1].ymin = 1e-1
    plots[-1].xmin = 0.
    plots[-1].xmax = 1000.
    plots[-1].xtit = "M(ee) [70, 110] (preselection, EB-EB) (GeV)"
    
    plots.append ( makeDefaultPlot ( "Mee_EBEE_70_110_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 2e3
    plots[-1].ymin = 1e-1
    plots[-1].xmin = 0.
    plots[-1].xmax = 1000.
    plots[-1].xtit = "M(ee) [70, 110] (preselection, EB-EE) (GeV)"
    
    
    plots.append ( makeDefaultPlot ( "Mee_EEEE_70_110_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 2e3
    plots[-1].ymin = 1e-1
    plots[-1].xmin = 0.
    plots[-1].xmax = 1000.
    plots[-1].xtit = "M(ee) [70, 110] (preselection, EE-EE) (GeV)"
    
    plots.append ( makeDefaultPlot ( "Mee_EB_70_110_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 6e3
    plots[-1].ymin = 1e-1
    plots[-1].xmin = 0.
    plots[-1].xmax = 1000.
    plots[-1].xtit = "M(ee) [70, 110] (preselection, EB-EE and EB-EB) (GeV)"
    
    plots.append ( makeDefaultPlot ( "Mee_70_110_Preselection" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 2e8
    plots[-1].ymin = 1e-1
    plots[-1].rebin = 4
    plots[-1].xmin = 70.
    plots[-1].xmax = 110.
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "M(ee) [70, 110] (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "Ptee_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 2
    plots[-1].ymax  = 4e6
    plots[-1].ymin  = 1e-1
    plots[-1].xmin  = 0
    plots[-1].xmax  = 1000
    plots[-1].ylog  = "yes"
    plots[-1].xtit  = "P_{T}(ee) (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "Mee_70_110_Preselection" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 4
    plots[-1].ymax = 2e8
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
    plots[-1].ymax = 1e5
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "no"
    plots[-1].xtit = "M(ej) average (GeV) [Preselection]"
    plots[-1].xmin  = 0
    plots[-1].xmax  = 1000
    plots[-1].rebin = 5
    
    
    plots.append ( makeDefaultPlot ( "Mej_selected_min_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 1e5
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "no"
    plots[-1].xtit = "M(ej) minimum (GeV) [Preselection]"
    plots[-1].xmin  = 0
    plots[-1].xmax  = 1000
    plots[-1].rebin = 5
    
    
    plots.append ( makeDefaultPlot ( "Mej_selected_max_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 1e5
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "no"
    plots[-1].xtit = "M(ej) maximum (GeV) [Preselection]"
    plots[-1].xmin  = 0
    plots[-1].xmax  = 1000
    plots[-1].rebin = 5
    
    # USED IN AN
    plots.append ( makeDefaultPlot ( "Mej_selected_avg_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 4e4
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "M(ej) average (GeV) [Preselection]"
    plots[-1].xmin  = 0
    plots[-1].xmax  = 1200
    plots[-1].rebin = 5
    
    
    plots.append ( makeDefaultPlot ( "Mej_selected_min_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 1e6
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "M(ej) minimum (GeV) [Preselection]"
    plots[-1].xmin  = 0
    plots[-1].xmax  = 1200
    plots[-1].rebin = 5
    
    
    plots.append ( makeDefaultPlot ( "Mej_selected_max_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 1
    plots[-1].ymax = 1000000
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
    plots[-1].xmax = 60.5
    plots[-1].ymin = 1e-1
    plots[-1].ymax = 5e5
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "n(vertices) [Preselection]"
    
    
    
    #plots.append ( makeDefaultPlot ( "nVertex_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].rebin = 1
    #plots[-1].xmin = -0.5
    #plots[-1].xmax = 40.5
    #plots[-1].ymin = 1e-1
    #plots[-1].ymax = 6e2
    #plots[-1].xtit = "n(vertices) [Preselection]"
    
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
    plots[-1].ymax = 2000000
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
    plots[-1].ymax = 1e6
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
    
####################################################################################################
#  FINAL SELECTION
####################################################################################################
if doFinalSelectionPlots:
    for mass_point in LQmassesFinalSelection:
        samplesForHistos = ["LQ_M{}".format(mass_point) ]
        keys             = ["LQ, M={} GeV".format(mass_point) ]

        plots.append ( makeDefaultPlot ( "Pt1stEle_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "1st Electron p_{T} (GeV) [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 2
    
        plots.append ( makeDefaultPlot ( "Pt2ndEle_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "2nd Electron p_{T} (GeV) [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 2
    
        plots.append ( makeDefaultPlot ( "Pt1stMuon_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "1st Muon p_{T} (GeV) [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 2
    
        plots.append ( makeDefaultPlot ( "Pt1stMuon_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "2nd Muon p_{T} (GeV) [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 2
    
        plots.append ( makeDefaultPlot ( "Eta1stEle_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "1st Electron #eta [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 2
    
        plots.append ( makeDefaultPlot ( "Eta2ndEle_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "2nd Electron #eta [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 2
    
        plots.append ( makeDefaultPlot ( "Phi1stEle_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "1st Electron #phi [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 4
    
        plots.append ( makeDefaultPlot ( "Phi2ndEle_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "2nd Electron #phi [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 4
    
        plots.append ( makeDefaultPlot ( "Pt1stJet_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "1st Jet p_{T} (GeV) [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 2
    
        plots.append ( makeDefaultPlot ( "Pt2ndJet_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "2nd Jet p_{T} (GeV) [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 2
    
        plots.append ( makeDefaultPlot ( "Eta1stJet_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "1st Jet #eta [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 2
    
        plots.append ( makeDefaultPlot ( "Eta2ndJet_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "2nd Jet #eta [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 2
    
        plots.append ( makeDefaultPlot ( "Phi1stJet_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "1st Jet #phi [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 4
    
        plots.append ( makeDefaultPlot ( "Phi2ndJet_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "2nd Jet #phi [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 4
    
        
        
        mej_rebin = 1
        mej_xmin = 0
        mej_xmax = 2300
        mee_rebin = 1
        st_rebin  = 1
        dr_rebin  = 2
        plots.append ( makeDefaultPlot ( "Me1j1_LQ" + str ( mass_point ) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].rebin = mej_rebin * 5
        plots[-1].xtit = "M_{e1j1} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin = mej_xmin
        plots[-1].xmax = mej_xmax
    
    
        plots.append ( makeDefaultPlot ( "Me1j2_LQ" + str ( mass_point ) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].rebin = mej_rebin * 5
        plots[-1].xtit = "M_{e1j2} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin = mej_xmin
        plots[-1].xmax = mej_xmax
    
    
        plots.append ( makeDefaultPlot ( "Me2j1_LQ" + str ( mass_point ) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].rebin = mej_rebin * 5
        plots[-1].xtit = "M_{e2j1} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin = mej_xmin
        plots[-1].xmax = mej_xmax
    
    
        plots.append ( makeDefaultPlot ( "Me2j2_LQ" + str ( mass_point ) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].rebin = mej_rebin * 5
        plots[-1].xtit = "M_{e2j2} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin = mej_xmin
        plots[-1].xmax = mej_xmax
    
        plots.append ( makeDefaultPlot ( "Mej_selected_avg_LQ" + str ( mass_point) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].rebin = mej_rebin
        plots[-1].xtit = "M_{ej}^{avg} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin = mej_xmin
        plots[-1].xmax = mej_xmax
    
        #plots.append ( makeDefaultPlot ( "Mej_selected_avg_LQ" + str ( mass_point) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        #plots[-1].rebin = mej_rebin
        #plots[-1].xtit = "M_{ej}^{avg} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].xmin = mej_xmin
        #plots[-1].xmax = mej_xmax
    
        plots.append ( makeDefaultPlot ( "Mej_selected_min_LQ" + str( mass_point ) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].rebin = mej_rebin
        plots[-1].xtit = "M_{ej}^{min} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin = mej_xmin
        plots[-1].xmax = mej_xmax
    
        plots.append ( makeDefaultPlot ( "Mej_selected_max_LQ" + str( mass_point ) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].rebin = mej_rebin
        plots[-1].xtit = "M_{ej}^{max} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin = mej_xmin
        plots[-1].xmax = mej_xmax
    
        plots.append ( makeDefaultPlot ( "Mej_minmax_LQ" + str( mass_point ), histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].rebin = mej_rebin
        plots[-1].xtit = "M_{ej} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].ytit = "2 #times Entries"
    
        plots.append ( makeDefaultPlot ( "sT_eejj_LQ" + str( mass_point ), histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit = "S_{T} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        
    
        plots.append ( makeDefaultPlot ( "sTlep_LQ" + str( mass_point ), histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit = "S_{T} (1st Electron, 2nd Electron) (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].rebin = 10
    
        plots.append ( makeDefaultPlot ( "sTjet_LQ" + str( mass_point ), histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit = "S_{T} (1st Jet, 2nd Jet) (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].rebin = 10
    
        plots.append ( makeDefaultPlot ( "sTfrac_Jet1_LQ" + str( mass_point )                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xmin = 0.0
        plots[-1].xmax = 1.0
        plots[-1].rebin = 4
        plots[-1].xtit = "Fraction S_{T} from jet 1 (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ( "sTfrac_Jet2_LQ" + str( mass_point )               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xmin = 0.0
        plots[-1].xmax = 1.0
        plots[-1].rebin = 4
        plots[-1].xtit = "Fraction S_{T} from jet 2 (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ( "sTfrac_Ele1_LQ" + str( mass_point )                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xmin = 0.0
        plots[-1].xmax = 1.0
        plots[-1].rebin = 4
        plots[-1].xtit = "Fraction S_{T} from ele 1 (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ( "sTfrac_Ele2_LQ" + str( mass_point )                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xmin = 0.0
        plots[-1].xmax = 1.0
        plots[-1].rebin = 4
        plots[-1].xtit = "Fraction S_{T} from ele 2 (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ( "sTfrac_Ele_LQ" + str( mass_point )                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xmin = 0.0
        plots[-1].xmax = 1.0
        plots[-1].rebin = 4
        plots[-1].xtit = "Fraction S_{T} from electrons (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ( "sTfrac_Jet_LQ" + str (mass_point) ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xmin = 0.0
        plots[-1].xmax = 1.0
        plots[-1].rebin = 4
        plots[-1].xtit = "Fraction S_{T} from jets (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ( "Ptee_LQ" + str( mass_point ), histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].rebin = 2
        plots[-1].xmin  = 0
        plots[-1].xmax  = 1000
        plots[-1].xtit  = "P_{T}(ee), (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ( "Mjj_LQ" + str( mass_point )          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
        plots[-1].rebin = 4
        plots[-1].xtit = "Dijet Mass (GeV) (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ( "Ptj1j2_LQ" + str( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].rebin = 2
        plots[-1].xmin  = 0
        plots[-1].xmax  = 1000
        plots[-1].xtit  = "P_{T}(j_{1}, j_{2}) (GeV), (LQ M = " + str ( mass_point ) + " selection)"
    
        plots.append ( makeDefaultPlot ( "nVertex_LQ" + str( mass_point )                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xmin = -0.5
        plots[-1].xmax = 40.5
        plots[-1].xtit = "n(vertexes) [Preselection]"
        
        plots.append ( makeDefaultPlot ( "Mee_LQ" + str ( mass_point ) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].rebin = mee_rebin
        plots[-1].xtit = "M(ee) (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ( "DR_Ele1Jet1_LQ" + str ( mass_point ) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].rebin = dr_rebin
        plots[-1].xtit = "#DeltaR(e1,j1) (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot2D ( "Mej_selected_min_vs_max_LQ" + str ( mass_point ) , histoBaseName2D_userDef, samplesForStackHistos, sampleForDataHisto ) )
        plots[-1].xrebin = 2
        plots[-1].yrebin = 2
        plots[-1].xtit = "M(ej)^{min} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].ytit = "M(ej)^{max} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
    
        plots.append ( makeDefaultPlot2D_Ratio ( "Mej_selected_min_vs_max_LQ" + str ( mass_point ) , histoBaseName2D_userDef, samplesForStackHistos, sampleForDataHisto ) )
        plots[-1].xrebin = 2
        plots[-1].yrebin = 2
        plots[-1].xtit = "M(ej)^{min} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].ytit = "M(ej)^{max} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot2D_NSigma ( "Mej_selected_min_vs_max_LQ" + str ( mass_point ) , histoBaseName2D_userDef, samplesForStackHistos, sampleForDataHisto ) )
        plots[-1].xrebin = 2
        plots[-1].yrebin = 2
        plots[-1].xtit = "M(ej)^{min} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].ytit = "M(ej)^{max} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
    
        
        plots.append ( makeDefaultPlot ("Classif_1stEle_LQ" + str ( int ( mass_point )  )               , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron classification, (LQ M = " + str ( mass_point ) + " selection)"
        
        
        plots.append ( makeDefaultPlot ("Classif_2ndEle_LQ" + str ( int ( mass_point )  )               , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron classification, (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ("CorrIsolation_1stEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron corrected HEEP isolation, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmax  = 1.0
        plots[-1].xmin  = -25
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("CorrIsolation_2ndEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron corrected HEEP isolation, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmax  = 1.0
        plots[-1].xmin  = -25
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("DeltaEtaTrkSC_1stEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        plots[-1].xtit  = "1st Electron #Delta#eta(track, SC), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  = -0.007
        plots[-1].xmax  =  0.007
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("DeltaEtaTrkSC_2ndEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        plots[-1].xtit  = "2nd Electron #Delta#eta(track, SC), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  = -0.007
        plots[-1].xmax  =  0.007
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("DeltaPhiTrkSC_1stEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        plots[-1].xtit  = "1st Electron #Delta#phi(track, SC), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  = -0.06
        plots[-1].xmax  =  0.06
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("DeltaPhiTrkSC_2ndEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        plots[-1].xtit  = "2nd Electron #Delta#phi(track, SC), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  = -0.06
        plots[-1].xmax  =  0.06
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("Full5x5E1x5OverE5x5_1stEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        plots[-1].xtit  = "1st Electron E_{1x5} / E_{5x5}, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  =  0.20
        plots[-1].xmax  =  1.40
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("Full5x5E1x5OverE5x5_2ndEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        plots[-1].xtit  = "2nd Electron E_{1x5} / E_{5x5}, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  =  0.20
        plots[-1].xmax  =  1.40
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("Full5x5E2x5OverE5x5_1stEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        plots[-1].xtit  = "1st Electron E_{2x5} / E_{5x5}, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  =  0.20
        plots[-1].xmax  =  1.40
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("Full5x5E2x5OverE5x5_2ndEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        plots[-1].xtit  = "2nd Electron E_{2x5} / E_{5x5}, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  =  0.20
        plots[-1].xmax  =  1.40
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("EcalIsolation_1stEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        plots[-1].xtit  = "1st Electron HEEP ECAL isolation, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  =   0.00
        plots[-1].xmax  =  15.00
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("EcalIsolation_2ndEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        plots[-1].xtit  = "2nd Electron HEEP ECAL isolation, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  =   0.00
        plots[-1].xmax  =  15.00
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("HcalIsolation_1stEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        plots[-1].xtit  = "1st Electron HEEP HCAL isolation, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  =   0.00
        plots[-1].xmax  =  15.00
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("HcalIsolation_2ndEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        plots[-1].xtit  = "2nd Electron HEEP HCAL isolation, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  =   0.00
        plots[-1].xmax  =  15.00
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("TrkIsolation_1stEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        plots[-1].xtit  = "1st Electron HEEP tracker isolation, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("TrkIsolation_2ndEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        plots[-1].xtit  = "2nd Electron HEEP tracker isolation, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("FBrem_1stEle_LQ" + str ( int ( mass_point )  )                 , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron brem fraction, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  =   -5.0
        plots[-1].xmax  =    4.0
        
        plots.append ( makeDefaultPlot ("FBrem_2ndEle_LQ" + str ( int ( mass_point )  )                 , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron brem fraction, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  =   -5.0
        plots[-1].xmax  =    4.0
        
        plots.append ( makeDefaultPlot ("GsfCtfCharge_1stEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron GSF CTF charge, (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ("GsfCtfCharge_2ndEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron GSF CTF charge, (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ("GsfCtfScPixCharge_1stEle_LQ" + str ( int ( mass_point )  )     , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron GSF CTF SC Pixel charge, (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ("GsfCtfScPixCharge_2ndEle_LQ" + str ( int ( mass_point )  )     , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron GSF CTF SC Pixel charge, (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ("GsfScPixCharge_1stEle_LQ" + str ( int ( mass_point )  )        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron GSF SC Pixel charge, (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ("GsfScPixCharge_2ndEle_LQ" + str ( int ( mass_point )  )        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron GSF SC Pixel charge, (LQ M = " + str ( mass_point ) + " selection)"
    
        plots.append ( makeDefaultPlot ("HasMatchedPhot_1stEle_LQ" + str ( int ( mass_point )  )        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron has matched photon, (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ("HasMatchedPhot_2ndEle_LQ" + str ( int ( mass_point )  )        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron has matched photon, (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ("HoE_1stEle_LQ" + str ( int ( mass_point )  )                   , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron H/E, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("HoE_2ndEle_LQ" + str ( int ( mass_point )  )                   , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron H/E, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("LeadVtxDistXY_1stEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron D_{XY} vs leading vertex, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("LeadVtxDistXY_2ndEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron D_{XY} vs leading vertex, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("LeadVtxDistZ_1stEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron D_{Z} vs leading vertex, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("LeadVtxDistZ_2ndEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron D_{Z} vs leading vertex, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("MissingHits_1stEle_LQ" + str ( int ( mass_point )  )           , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron N(missing hits), (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ("MissingHits_2ndEle_LQ" + str ( int ( mass_point )  )           , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron N(missing hits), (LQ M = " + str ( mass_point ) + " selection)"
    
        plots.append ( makeDefaultPlot ("NBrems_1stEle_LQ" + str ( int ( mass_point )  )                , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron N(brems), (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ("NBrems_2ndEle_LQ" + str ( int ( mass_point )  )                , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron N(brems), (LQ M = " + str ( mass_point ) + " selection)"
        
        plots.append ( makeDefaultPlot ("EnergyORawEnergy_1stEle_LQ" + str ( int ( mass_point )  )      , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron energy correction factor, (LQ M = " + str ( mass_point ) + " selection)" 
        plots[-1].xmin  = 0.95
        plots[-1].xmax  = 1.3
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("EnergyORawEnergy_2ndEle_LQ" + str ( int ( mass_point )  )      , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron energy correction factor, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  = 0.95
        plots[-1].xmax  = 1.3
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("SigmaEtaEta_Barrel_1stEle_LQ" + str ( int ( mass_point )  )    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron #sigma_{#eta#eta} [Preselection, Barrel]"
        plots[-1].xmin  = 0.005
        plots[-1].xmax  = 0.015
        plots[-1].rebin = 2
        
        plots.append ( makeDefaultPlot ("SigmaEtaEta_Barrel_2ndEle_LQ" + str ( int ( mass_point )  )    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron #sigma_{#eta#eta} [Preselection, Barrel]"
        plots[-1].xmin  = 0.005
        plots[-1].xmax  = 0.015
        plots[-1].rebin = 2
        
        plots.append ( makeDefaultPlot ("SigmaEtaEta_Endcap_1stEle_LQ" + str ( int ( mass_point )  )    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron #sigma_{#eta#eta} [Preselection, Endcap]"
        plots[-1].xmin  = 0.005
        plots[-1].xmax  = 0.04
        plots[-1].rebin = 2
        
        plots.append ( makeDefaultPlot ("SigmaEtaEta_Endcap_2ndEle_LQ" + str ( int ( mass_point )  )    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron #sigma_{#eta#eta} [Preselection, Endcap]"
        plots[-1].xmin  = 0.005
        plots[-1].xmax  = 0.04
        plots[-1].rebin = 2
        
        #plots.append ( makeDefaultPlot ("Full5x5SigmaIEtaIEta_Barrel_1stEle_LQ" + str ( int ( mass_point )  )  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron #sigma_{i#etai#eta} [Preselection, Barrel]"
        #plots[-1].xmin  = 0.005
        #plots[-1].xmax  = 0.015
        #plots[-1].rebin = 2
        #
        #plots.append ( makeDefaultPlot ("Full5x5SigmaIEtaIEta_Barrel_2ndEle_LQ" + str ( int ( mass_point )  )  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "2nd Electron #sigma_{i#etai#eta} [Preselection, Barrel]"
        #plots[-1].xmin  = 0.005
        #plots[-1].xmax  = 0.015
        #plots[-1].rebin = 2
        #
        #plots.append ( makeDefaultPlot ("Full5x5SigmaIEtaIEta_Endcap_1stEle_LQ" + str ( int ( mass_point )  )  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron #sigma_{i#etai#eta} [Preselection, Endcap]"
        #plots[-1].xmin  = 0.015
        #plots[-1].xmax  = 0.035
        #plots[-1].rebin = 2
        #
        #plots.append ( makeDefaultPlot ("Full5x5SigmaIEtaIEta_Endcap_2ndEle_LQ" + str ( int ( mass_point )  )  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "2nd Electron #sigma_{i#etai#eta} [Preselection, Endcap]"
        #plots[-1].xmin  = 0.015
        #plots[-1].xmax  = 0.035
        #plots[-1].rebin = 2
        
        plots.append ( makeDefaultPlot ("TrkPtOPt_1stEle_LQ" + str ( int ( mass_point )  )              , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron track pt / SC pt, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("TrkPtOPt_2ndEle_LQ" + str ( int ( mass_point )  )              , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron track pt / SC pt, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("ValidFrac_1stEle_LQ" + str ( int ( mass_point )  )             , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "1st Electron valid fraction of hits, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  = 0.00
        plots[-1].xmax  = 1.50
        plots[-1].rebin = 4
        
        plots.append ( makeDefaultPlot ("ValidFrac_2ndEle_LQ" + str ( int ( mass_point )  )              , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        plots[-1].xtit  = "2nd Electron valid fraction of hits, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  = 0.00
        plots[-1].xmax  = 1.50
        plots[-1].rebin = 4
    
        plots.append ( makeDefaultPlot ( "EleChargeSum_LQ" + str( int (mass_point))        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "Electron 1 charge + Electron 2 charge, (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin  = -1.5
        plots[-1].xmax  = 1.5
    
        plots.append ( makeDefaultPlot ( "MET_LQ"+ str( int (mass_point)) ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit = "PFMET (GeV), (LQ M = " + str ( mass_point ) + " selection)" 
        plots[-1].rebin = 4
        plots[-1].xmax = 500
        plots[-1].xmin = 0

#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()

fileps = "allPlots_eejj_analysis.pdf"

c.Print(fileps + "[")
for i_plot, plot in enumerate(plots):
    #print 'draw plot:',plot
    plot.Draw(fileps, i_plot + 1)
c.Print(fileps+"]")

makeTOC ( "allPlots_eejj_analysis_toc.tex" , fileps, plots ) 

