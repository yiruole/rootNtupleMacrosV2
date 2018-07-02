#!/usr/bin/env python

from plot_class import *
from ROOT import *

#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_local_nov18_addStSFplots_allDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_noIncWStitch.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_local_nov18_addStSFplots_allDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_QCD_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_local_nov18_addStSFplots_ICHEPDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_local_nov18_addStSFplots_ICHEPDataAndMC_ele27wptightOrPhoton175Data2015CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_QCD_plots.root")

#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jan20_rereco_stitch120_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jan20_rereco_stitch120_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_scaled.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jan20_rereco_stitch120_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_QCD_jan22_rereco_enujj2012FinSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")

#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jan27_rereco_stitch120_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_scaled.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_QCD_jan25_rereco_enujj2012FinSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")

# unscaled
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb28_recoHeepSF_rereco_stitchDYW_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb28_recoHeepSF_rereco_stitchDYW_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_QCD_feb27_rereco_enujj2012FinSels//output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")

# unscaled
# no top pt weighting
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_mar28_recoHeepSF_rereco_stitchDYW_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
# with top
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_mar30_topPtWeight_recoHeepSF_rereco_ele27wptightEta2p1CurveMC_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
# new trigger
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_apr11_ele27wptightOREle115_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_feb28_recoHeepSF_rereco_stitchDYW_ele27wptightEta2p1CurveMC/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
# unscaled
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_may8_ele27wptightOREle115_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
# scaled
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_may8_ele27wptightOREle115_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
# unscaled, old wp27 trigger curve
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_may10_ele27wptightEta2p1_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
# scaled, old wp27 trigger curve
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_may10_ele27wptightEta2p1_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
# unscaled ele27 from Bibhu
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_may11_ele27wptight_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
# scaled ele27 from Bibhu
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_may11_ele27wptight_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")

#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_jun1_ele27wptightOREle115ORPhoton175/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_sep25_btags__enujjOptFinalSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_oct6_updateFinalSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_oct8_updatedFinalSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_oct16_newFR_updatedFinalSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_oct18_mejEtaPlots_updatedFinalSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_oct19_reminiAODFR_updatedFinalSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_oct23_updateFR_updatedFinalSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_psk_oct24_updateFRAndOct21QCDRSK_updatedFinalSels/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_rsk_oct28_updateFRAndOct21QCDRSK_noMuonReq/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_rsk_oct31_updateFRSCEtAndOct21QCDRSK/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_rsk_nov10_failHEEP_updateFRSCEt_preselOnly/output_cutTable_lq_enujj_MT_QCD_preselOnly/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_rsk_nov14_attemptMjjControlRegions_preselOnly/output_cutTable_lq_enujj_MT_QCD_preselOnly/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk_nov15_attemptMjjControlRegions_preselOnly/output_cutTable_lq_enujj_MT_QCD_preselOnly/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk_nov16_preselOnly/output_cutTable_lq_enujj_MT_QCD_preselOnly/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk_nov16_addDiElectronPlot_preselOnly/output_cutTable_lq_enujj_MT_QCD_preselOnly/analysisClass_lq_enujj_QCD_plots.root")
#
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk_nov17_highMejPlots_preselOnly/output_cutTable_lq_enujj_MT_QCD_preselOnly/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk_dec10_withAndWithoutSt_preselOnly/output_cutTable_lq_enujj_MT_QCD_preselOnly_fakeRateETHLTWithStMET/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk_dec13_preselOnly/output_cutTable_lq_enujj_MT_QCD_preselOnly/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_jan8_preselOnly/output_cutTable_lq_enujj_MT_QCD_preselOnly/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk237_jan17_preselOnly/output_cutTable_lq_enujj_MT_QCD_preselOnly/analysisClass_lq_enujj_QCD_plots.root")
# MET 100
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk237_jan31_gsfEtaCheck_MET100/output_cutTable_lq_enujj_MT_QCD_preselOnly_MET100/analysisClass_lq_enujj_QCD_plots.root")
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk237_feb2_gsfEtaCheck_MET100_addPlots/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
# with Pt(e,met)
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_newRsk237_feb4_gsfEtaCheck_MET100_PtEMET70/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_feb10_bugfix/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_feb10_addPlot/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_feb14_dPhiEleMET0p8/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_feb27_addMoreFinalSelPlots/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_mar16_fixMuons/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#
#File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_mar20_addPlots/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")
#
File_QCD_preselection = GetFile(os.environ["LQDATA"] + "/2016qcd/enujj_mar26_addMTPlots/output_cutTable_lq_enujj_MT_QCD/analysisClass_lq_enujj_QCD_plots.root")

#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jun1_ele27wptightOREle115ORPhoton175_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jun1_ele27wptightOREle115ORPhoton175_enujj2012FinSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jul4_ele27wptightOREle115ORPhoton175_enujjOptFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root.jul4")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jul4_ele27wptightOREle115ORPhoton175_enujjOptFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")

#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jul4_ele27wptightOREle115ORPhoton175_enujjOptFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_sep13_ele27wptightOREle115ORPhoton175_enujjPowhegOptFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")

#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_oct3_bTagFix_enujjPowhegOptFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_oct5_noTrigEffMC_enujjPowhegOptFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_oct6_finerBinnedTrigEff_enujjPowhegOptFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_oct6_finerBinnedTrigEff_updateFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_oct8_finerBinnedTrigEff_updatedFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_oct15_extendMTRange_updatedFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_scaled.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_oct18_mejEtaPlots_updatedFinalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
##
##File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_rsk_oct28_noMuonReq/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_nov10_reprodNominal/output_cutTable_lq_enujj_MT_preselOnly/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_nov14_attemptMjjControlRegions_preselOnly/output_cutTable_lq_enujj_MT_preselOnly/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_nov16_preselOnly/output_cutTable_lq_enujj_MT_preselOnly/analysisClass_lq_enujj_MT_plots.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_nov16_addDiElectronPlot_preselOnly/output_cutTable_lq_enujj_MT_preselOnly/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_nov17_newTrigEff_preselOnly/output_cutTable_lq_enujj_MT_preselOnly/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_nov17_highMejPlots_preselOnly/output_cutTable_lq_enujj_MT_preselOnly/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_nov22_fixTrigEff_usePtHeep_preselOnly/output_cutTable_lq_enujj_MT_preselOnly/analysisClass_lq_enujj_MT_plots.root")
# scaled
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_dec13_fixTrigEff_usePtHeep_preselOnly/output_cutTable_lq_enujj_MT_preselOnly/analysisClass_lq_enujj_MT_plots.root")
# unscaled
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_dec13_fixTrigEff_usePtHeep_preselOnly/output_cutTable_lq_enujj_MT_preselOnly/analysisClass_lq_enujj_MT_plots_unscaled.root")
# scaled
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jan20_usePtHeep_finalSels/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")

# MET100
# unscaled
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jan31_v237_MET100/output_cutTable_lq_enujj_MT_preselOnly_MET100/analysisClass_lq_enujj_MT_plots_unscaled.root")
# scaled
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_jan31_v237_MET100/output_cutTable_lq_enujj_MT_preselOnly_MET100/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb2_v237_MET100_addMTPlots/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
# with Pt(e,MET) > 70
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb4_v237_MET100_PtEMET70/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb4_v237_MET100_PtEMET70/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb10_v237_bugfix/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb10_v237_bugfix/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb12_v237_fixPlot/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb14_dPhiEleMet0p8/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb14_dPhiEleMet0p8/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb20_dPhiEleMet0p8_newSingTop/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb20_dPhiEleMet0p8_newSingTop/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_feb27_dPhiEleMet0p8_fixMETPlot/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_mar5_removeTopPtReweight/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_mar5_removeTopPtReweight/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_mar16_fixMuons/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_mar16_fixMuons/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_mar20_addPlots/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
#
#File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_mar26_addMTPlots/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots_unscaled.root")
File_preselection     = GetFile(os.environ["LQDATA"] + "/2016analysis/enujj_psk_mar26_addMTPlots/output_cutTable_lq_enujj_MT/analysisClass_lq_enujj_MT_plots.root")


LQmasses = [650,1200]
LQmassesFinalSelection = [400, 650, 1000]

#### Common values for plots:
#otherBkgsKey="QCD, single top, VV+jets, W+jets"
zUncBand="no"
bkgUncBand=True
makeRatio=1
makeNSigma=1
doExtraPlots = False

doPreselectionPlots = True
doFinalSelectionPlots = False

pt_rebin   = 2
eta_rebin  = 2
st_rebin   = 2
mt_rebin   = 2
mass_rebin = 2
dphi_rebin = 2
dr_rebin   = 2
mej_rebin = 8

#met_bins = [0,5,10,15,20,25,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,400,500,700,1000,1400,2000]
met_bins = range(0,105,5)
met_bins.extend(range(110,210,10))
met_bins.extend(range(220,320,20))
met_bins.extend([340,380,420,480,540,600,650,700,800,900,1000,1200,1500,2000])
print(met_bins)

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
#keysStack             = [ "Other backgrounds", "t#bar{t} (MG5_aMC)"  ,  "Z/#gamma* + jets (MG5_aMC)"  ]

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
#samplesForStackHistos_other = [ "OTHERBKG_ZJetPt" ]
samplesForStackHistos_other = [ "OTHERBKG_ZJetPt_amcAtNLODiboson" ]
#samplesForStackHistos_WJets  = [ "TTbar_Madgraph", "WJet_amcatnlo_ptBinned" ]
##keysStack             = [ "Other backgrounds", "t#bar{t} (Madgraph)"  ,  "Z/#gamma* + jets (MG HT)"  ]
#keysStack             = [ "Other backgrounds", "QCD multijet", "t#bar{t} (Madgraph)"  ,  "W + jets (MG5_aMC Pt)"  ]

## amcatnlo TTBar and WJets amcAtNLO Pt
#samplesForStackHistos_WJets  = [ "TTbar_amcatnlo_Inc", "WJet_amcatnlo_ptBinned" ]
#keysStack             = ["QCD multijet (data)",  "Other backgrounds", "t#bar{t} (MG5_aMC)"  ,  "W + jets (MG5_aMC Pt)"  ]
# powheg TTBar and WJets amcAtNLO Pt
samplesForStackHistos_WJets  = [ "TTbar_powheg", "WJet_amcatnlo_ptBinned" ]
keysStack             = ["QCD multijet (data)",  "Other backgrounds", "t#bar{t} (powheg)"  ,  "W + jets (MG5_aMC Pt)"  ]
systTypes             = ['qcd', 'mc', 'ttbar', 'wjets']
# WJets inclusive samples
## amcatnlo TTBar and WJets amcAtNLO Inc
#samplesForStackHistos_WJets  = [ "TTbar_amcatnlo_Inc", "WJet_amcatnlo_Inc" ]
#keysStack             = ["QCD multijet (data)",  "Other backgrounds", "t#bar{t} (MG5_aMC)"  ,  "W + jets (MG5_aMC Inc.)"  ]
## powheg TTBar and WJets amcAtNLO Inc
#samplesForStackHistos_WJets  = [ "TTbar_powheg", "WJet_amcatnlo_Inc" ]
#keysStack             = ["QCD multijet (data)",  "Other backgrounds", "t#bar{t} (powheg)"  ,  "W + jets (MG5_aMC Inc.)"  ]

# QCD
samplesForStackHistos_QCD = ["QCDFakes_DATA"]
#samplesForStackHistos_QCD = ["QCD_EMEnriched"]
#keysForStackHistos_QCD = ["QCD multijet (data)"]


#samplesForStackHistos_ZJets  = [ "TTbar_FromData", "ZJet_Madgraph" ]
# older
#samplesForStackHistos = samplesForStackHistos_other + samplesForStackHistos_ZJets
samplesForStackHistos = samplesForStackHistos_QCD + samplesForStackHistos_other + samplesForStackHistos_WJets
#keysStack             = [ "Other backgrounds", "QCD multijet", "t#bar{t} (Madgraph)"  ,  "W + jets (MG HT)"  ]
stackColorIndexes     = [ kCyan         , 9                  ,   600            , kRed  ]
stackFillStyleIds     = [ 1001          , 1001               ,  1001            , 1001  ]

##keysStack             = [ "Other backgrounds", "t#bar{t} (Madgraph)"  ,  "Z/#gamma* + jets (MG Inc)"  ]
##stackColorIndexes     = [ 9                  , 600         ,  kRed           ]
###stackFillStyleIds     = [ 3008               , 3004        ,  3345           ]
##stackFillStyleIds     = [ 1001               , 1001        ,  1001           ]
#
##keysStack             = ["t#bar{t} (Madgraph)"     , "SingleTop"  ,"DIBOSON"   ,  "WJet(MG5_aMC)",   "PhotonJets_Madgraph",  "QCD_EMEnriched",     "Z/#gamma* + jets"  ]
##stackColorIndexes     = [   600                    , kMagenta     , kGreen+3   ,  kCyan-3        ,           7            ,    kMagenta+3    ,     kRed+1           ]
##stackFillStyleIds     = [   1001                   , 1001         ,    1001    ,  1001           ,           1001         ,        1001      ,     3345           ]
##keysStack             = [ "WJet(MG5_aMC)"    , "SingleTop"  ,  "QCD_EMEnriched",  "DIBOSON"   , "t#bar{t} (Madgraph)"  ,  "Z/#gamma* + jets"  ]
##stackColorIndexes     = [ 9                  , kGreen+3     ,        6         ,        3     ,    600                 ,  kRed           ]
##stackFillStyleIds     = [ 1001               , 1001         ,        1001      ,      1001    ,    3004                ,  3345           ]
#
##stackColorIndexes.reverse()
##stackFillStyleIds.reverse()

#samplesForHistos = ["LQ_M1200"      ]
#keys             = ["LQ, M=1200 GeV"]
#samplesForHistos = ["LQ_M600", "LQ_M1200"]
#keys             = ["LQ, M=600 GeV", "LQ, M=1200 GeV"]
samplesForHistos = ["LQ_M{}".format(lqmass) for lqmass in LQmasses      ]
keys             = ["LQ, M={} GeV".format(lqmass) for lqmass in LQmasses ]
# no signal
#samplesForHistos = []
#keys             = []


samplesForHistos_blank = []
keys_blank             = []

sampleForDataHisto = "DATA"
#dataBlindAbovePt1 = 800 # GeV; used for ele Pt1, MET, Mej, MT
#dataBlindAboveSt = 1500 # for St plots
dataBlindAbovePt1 = -1 # GeV; used for ele Pt1, MET, Mej, MT
dataBlindAboveSt = -1 # for St plots

def makeDefaultPlot ( variableName, histoBaseName, 
                      samplesForHistos, keys,
                      samplesForStackHistos, keysStack,
                      sampleForDataHisto,
                      zUncBand, makeRatio, dataBlindAbove = -1):
    ## special handling for Pt1st ele for QCD
    #if variableName=='Pt1stEle_PAS':
    #  variableNameQCD='SCEt1stEle_Presel'
    ##elif variableName=='Ele1_PtHeep':
    ##  variableNameQCD='Ele1_SCEt'
    #elif variableName=='PtHeep1stEle_Presel':
    #  variableNameQCD='SCEt1stEle_Presel'
    #else:
    #  variableNameQCD=variableName
    variableNameQCD=variableName
    
    plot                   = Plot() 
    plot.histosStack       =  ( 
                                generateHistoList( histoBaseName, samplesForStackHistos_QCD, variableNameQCD, File_QCD_preselection ) +
                                generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection   ) +
                                generateHistoList( histoBaseName, samplesForStackHistos_WJets, variableName, File_preselection   ) ) 
    plot.keysStack         = keysStack
    plot.systs             = [GetBackgroundSyst(systType,False) for systType in systTypes]
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
    #plot.gif_folder        = "gif_enujj_scaled_preselectionOnly/"
    #plot.eps_folder        = "eps_enujj_scaled_preselectionOnly/"
    plot.pdf_folder        = "pdf_enujj_scaled_preselectionOnly/"
    plot.png_folder        = "png_enujj_scaled_preselectionOnly/"
    #plot.c_folder          = "c_enujj_scaled_preselectionOnly/"
    plot.root_folder       = "root_enujj_scaled_preselectionOnly/"
    plot.suffix            = "enujj"
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

if doPreselectionPlots:
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
    
    ##plots.append ( makeDefaultPlot ( "nEle", histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots.append ( makeDefaultPlot ( "nEle", histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].xtit            = "Number of electrons [Preselection]"
    #plots[-1].ylog            = "yes"
    #plots[-1].rebin           = 1
    #plots[-1].xmin            = -0.5
    #plots[-1].xmax            = 6.5
    #plots[-1].ymin            = 0.0001
    #plots[-1].ymax            = 10000000000
    
    plots.append ( makeDefaultPlot ( "nElectron_PAS", histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit            = "Number of electrons [Preselection]"
    plots[-1].ylog            = "yes"
    plots[-1].rebin           = 1
    plots[-1].xmin            = -0.5
    plots[-1].xmax            = 6.5
    plots[-1].ymin            = 0.0001
    plots[-1].ymax            = 10000000000
    
    plots.append ( makeDefaultPlot ( "nMuon_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots.append ( makeDefaultPlot ( "nMuon" ,  histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "Number of muons [Preselection]"
    plots[-1].ylog            = "yes"
    plots[-1].rebin           = 1
    plots[-1].xmin            = -0.5
    plots[-1].xmax            = 6.5
    plots[-1].ymin            = 0.0001
    plots[-1].ymax            = 10000000000
    
    #plots.append ( makeDefaultPlot ( "PtHeep1stEle_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    ##plots.append ( makeDefaultPlot ( "Ele1_PtHeep"	         ,  histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    #plots[-1].xtit = "1st Electron p_{T} (GeV) [Preselection]"
    #plots[-1].xmin = 0.
    #plots[-1].xmax = 1500.
    #plots[-1].ymax = 200000
    #plots[-1].ymin = 1e-1
    #plots[-1].ylog  = "yes"
    #plots[-1].rebin = pt_rebin
    
    plots.append ( makeDefaultPlot ( "PtHeep1stEle_Presel"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    #plots.append ( makeDefaultPlot ( "Ele1_PtHeep"	         ,  histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    plots[-1].xtit = "1st Electron p_{T} (GeV) [Preselection]"
    plots[-1].xmin = 0.
    plots[-1].xmax = 1500.
    plots[-1].ymax = 200000
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    plots[-1].rebin = pt_rebin
    
    plots.append ( makeDefaultPlot ( "SCEt1stEle_Presel"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    #plots.append ( makeDefaultPlot ( "Ele1_PtHeep"	         ,  histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    plots[-1].xtit = "1st Electron SC p_{T} (GeV) [Preselection]"
    plots[-1].xmin = 0.
    plots[-1].xmax = 1500.
    plots[-1].ymax = 200000
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    plots[-1].rebin = pt_rebin
    
    #plots.append ( makeDefaultPlot ( "Pt1stMuon_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    #plots[-1].xtit = "1st Muon p_{T} (GeV) [Preselection]"
    #plots[-1].xmin = 0.
    #plots[-1].xmax = 500.
    #plots[-1].ymax = 200000
    #plots[-1].ymin = 1e-1
    #plots[-1].ylog  = "yes"
    ##plots[-1].rebin = pt_rebin
    #plots.append ( makeDefaultPlot ( "Pt2ndMuon_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    #plots[-1].xtit = "2nd Muon p_{T} (GeV) [Preselection]"
    #plots[-1].xmin = 0.
    #plots[-1].xmax = 200.
    #plots[-1].ymax = 200000
    #plots[-1].ymin = 1e-1
    #plots[-1].ylog  = "yes"
    ##plots[-1].rebin = pt_rebin
    
    plots.append ( makeDefaultPlot ( "Eta1stEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "1st Electron #eta [Preselection]"   
    plots[-1].ymax = 5e10
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
    
    plots.append ( makeDefaultPlot ( "MET_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio ) )
    plots[-1].xtit = "PFMET (GeV) [Preselection]"
    #plots[-1].rebin = pt_rebin
    plots[-1].rebin = 'var'
    plots[-1].xbins = met_bins
    plots[-1].ymax = 200000
    plots[-1].ymin = 1e-1
    plots[-1].xmax = 2500
    plots[-1].xmin = 0
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "METPhi_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "PFMET #phi [Preselection]"
    plots[-1].rebin = 1
    plots[-1].ymax = 6000
    plots[-1].ymin = 1e-1
    #plots[-1].ylog  = "yes"
    
    
    #plots.append ( makeDefaultPlot ( "MET_Type01_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    #plots[-1].xtit = "PFMET0+1 (GeV) [Preselection]"
    #plots[-1].rebin = pt_rebin
    #plots[-1].ymax = 200000000
    #plots[-1].ymin = 1e-1
    #plots[-1].xmax = 600
    #plots[-1].xmin = 0
    #plots[-1].ylog  = "yes"
    #
    #plots.append ( makeDefaultPlot ( "MET_Type01_Phi_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].xtit = "PFMET0+1 #phi [Preselection]"
    #plots[-1].rebin = 1
    #plots[-1].ymax = 200000000000
    #plots[-1].ymin = 1e-1
    #plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "minMETPt1stEle_PAS"    ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = pt_rebin
    plots[-1].ymax = 1e7
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
    plots[-1].xmax = 2000
    plots[-1].xmin = 0
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "Pt2ndJet_PAS"          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = pt_rebin
    plots[-1].ymax = 200000
    plots[-1].ymin = 1e-1
    plots[-1].xmax = 2000
    plots[-1].xmin = 0
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "2nd Jet p_{T} (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "Eta1stJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = eta_rebin
    plots[-1].ymax = 5e10
    plots[-1].ymin = 1e-1
    plots[-1].xmin = -3.
    plots[-1].xmax =  3.
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "1st Jet #eta [Preselection]"
    
    plots.append ( makeDefaultPlot ( "Eta2ndJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = eta_rebin
    plots[-1].ymax = 5e10
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
    
    plots.append ( makeDefaultPlot ( "CISV1stJet_PAS"        ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 5
    plots[-1].ymax = 20000000
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "1st Jet CISV [Preselection]"
    
    plots.append ( makeDefaultPlot ( "CISV2ndJet_PAS"        ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = 5
    plots[-1].ymax = 200000000
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "2nd Jet CISV [Preselection]"
    
    plots.append ( makeDefaultPlot ( "MTenu_PAS"             ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    plots[-1].xtit = "M_{T} (e, PFMET) (GeV) [Preselection]"
    #plots[-1].rebin = mt_rebin
    plots[-1].rebin = 'var'
    plots[-1].xbins = met_bins
    plots[-1].ymax = 20000000
    plots[-1].ymin = 1e-1
    plots[-1].xmax = 2000
    plots[-1].xmin = 0
    plots[-1].ylog  = "yes"
    
    #plots.append ( makeDefaultPlot ( "MTenu_PAS"             ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    plots.append ( makeDefaultPlot ( "MTenu_PAS"             ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "M_{T} (e, PFMET) (GeV) [Preselection]"
    #plots[-1].rebin = mt_rebin
    plots[-1].rebin = 'var'
    plots[-1].xbins = met_bins
    plots[-1].ymax = 20000000
    plots[-1].ymin = 1e-1
    plots[-1].xmax = 2000
    plots[-1].xmin = 0
    plots[-1].ylog  = "yes"
    
    ##plots.append ( makeDefaultPlot ( "MTenu_Type01_PAS"             ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    #plots.append ( makeDefaultPlot ( "MTenu_Type01_PAS"             ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].xtit = "M_{T} (e, PFMET0+1) (GeV) [Preselection]"
    #plots[-1].rebin = mt_rebin
    #plots[-1].ymax = 20000000
    #plots[-1].ymin = 1e-1
    #plots[-1].xmax = 1000
    #plots[-1].xmin = 0
    #plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MTenu_50_110"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "50 < M_{T} (e, PFMET) < 110 (GeV) [Preselection] "
    plots[-1].rebin = mt_rebin
    plots[-1].ymax = 2e4
    plots[-1].ymin = 200
    plots[-1].xmax = 145
    plots[-1].xmin = 40
    plots[-1].ylog  = "yes"
    
    
    #plots.append ( makeDefaultPlot ( "MTenu_50_110_Njet_gte5"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].xtit = "50 < M_{T} (e, PFMET) < 110 (GeV) [Preselection + N(Jet) #geq 5]"
    #plots[-1].rebin = mt_rebin
    #plots[-1].ymax = 3e3
    #plots[-1].ymin = 1
    #plots[-1].xmax = 145
    #plots[-1].xmin = 40
    #plots[-1].ylog  = "yes"
    
    
    #plots.append ( makeDefaultPlot ( "MTenu_50_110_Njet_lte4"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].xtit = "50 < M_{T} (e, PFMET) < 110 (GeV) [Preselection + N(Jet) #leq 4]"
    #plots[-1].rebin = mt_rebin
    #plots[-1].ymax = 2e4
    #plots[-1].ymin = 200
    #plots[-1].xmax = 145
    #plots[-1].xmin = 40
    #plots[-1].ylog  = "yes"
    
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
    
    plots.append ( makeDefaultPlot ( "Pt1stEle_MTenu_50_110_Njet_gte4"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "Pt 1st electron [50 < M_{T} (e, PFMET) < 110 (GeV) + Preselection + N(Jet) #geq 4]"
    plots[-1].xmin = 0.
    plots[-1].xmax = 1500.
    plots[-1].ymax = 200000
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    plots[-1].rebin = pt_rebin
    
    plots.append ( makeDefaultPlot ( "Pt1stJet_MTenu_50_110_Njet_gte4"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "Pt 1st jet [50 < M_{T} (e, PFMET) < 110 (GeV) + Preselection + N(Jet) #geq 4]"
    plots[-1].rebin = pt_rebin
    plots[-1].ymax = 20000000
    plots[-1].ymin = 1e-1
    plots[-1].xmax = 1500
    plots[-1].xmin = 0
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "Pt2ndJet_MTenu_50_110_Njet_gte4"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "Pt 2nd jet [50 < M_{T} (e, PFMET) < 110 (GeV) + Preselection + N(Jet) #geq 4]"
    plots[-1].rebin = pt_rebin
    plots[-1].ymax = 20000000
    plots[-1].ymin = 1e-1
    plots[-1].xmax = 1500
    plots[-1].xmin = 0
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MET_MTenu_50_110_Njet_gte4"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "PFMET [50 < M_{T} (e, PFMET) < 110 (GeV) + Preselection + N(Jet) #geq 4]"
    plots[-1].rebin = pt_rebin
    plots[-1].ymax = 20000000
    plots[-1].ymin = 1e-1
    plots[-1].xmax = 1500
    plots[-1].xmin = 0
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "Pt1stEle_MTenu_50_110_Njet_lte3"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "Pt 1st electron [50 < M_{T} (e, PFMET) < 110 (GeV) + Preselection + N(Jet) #leq 3]"
    plots[-1].xmin = 0.
    plots[-1].xmax = 1500.
    plots[-1].ymax = 200000
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    plots[-1].rebin = pt_rebin
    
    plots.append ( makeDefaultPlot ( "Pt1stJet_MTenu_50_110_Njet_lte3"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "Pt 1st jet [50 < M_{T} (e, PFMET) < 110 (GeV) + Preselection + N(Jet) #leq 3]"
    plots[-1].rebin = pt_rebin
    plots[-1].ymax = 20000000
    plots[-1].ymin = 1e-1
    plots[-1].xmax = 1500
    plots[-1].xmin = 0
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "Pt2ndJet_MTenu_50_110_Njet_lte3"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "Pt 2nd jet [50 < M_{T} (e, PFMET) < 110 (GeV) + Preselection + N(Jet) #leq 3]"
    plots[-1].rebin = pt_rebin
    plots[-1].ymax = 20000000
    plots[-1].ymin = 1e-1
    plots[-1].xmax = 1500
    plots[-1].xmin = 0
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MET_MTenu_50_110_Njet_lte3"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "PFMET [50 < M_{T} (e, PFMET) < 110 (GeV) + Preselection + N(Jet) #leq 3]"
    plots[-1].rebin = pt_rebin
    plots[-1].ymax = 20000000
    plots[-1].ymin = 1e-1
    plots[-1].xmax = 1500
    plots[-1].xmin = 0
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MTenu_50_110_gteOneBtaggedJet"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "50 < M_{T} (e, PFMET) < 110 (GeV) [Preselection + N(B-tag jets) #geq 1]"
    plots[-1].rebin = mt_rebin
    plots[-1].ymax = 4e4
    plots[-1].ymin = 200
    plots[-1].xmax = 145
    plots[-1].xmin = 40
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MTenu_50_110_noBtaggedJets"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "50 < M_{T} (e, PFMET) < 110 (GeV) [Preselection + N(B-tag jets) = 0]"
    plots[-1].rebin = mt_rebin
    plots[-1].ymax = 4e4
    plots[-1].ymin = 200
    plots[-1].xmax = 145
    plots[-1].xmin = 40
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MTenu_70_150_gteOneBtaggedJet"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "70 < M_{T} (e, PFMET) < 150 (GeV) [Preselection + N(B-tag jets) #geq 1]"
    plots[-1].rebin = mt_rebin
    plots[-1].ymax = 4e4
    plots[-1].ymin = 200
    plots[-1].xmax = 160
    plots[-1].xmin = 60
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MTenu_70_150_noBtaggedJets"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "70 < M_{T} (e, PFMET) < 150 (GeV) [Preselection + N(B-tag jets) = 0]"
    plots[-1].rebin = mt_rebin
    plots[-1].ymax = 4e4
    plots[-1].ymin = 200
    plots[-1].xmax = 160
    plots[-1].xmin = 60
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MTenu_110_190_gteOneBtaggedJet"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "110 < M_{T} (e, PFMET) < 190 (GeV) [Preselection + N(B-tag jets) #geq 1]"
    plots[-1].rebin = mt_rebin
    plots[-1].ymax = 4e4
    plots[-1].ymin = 200
    plots[-1].xmax = 145
    plots[-1].xmin = 40
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MTenu_110_190_noBtaggedJets"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "110 < M_{T} (e, PFMET) < 190 (GeV) [Preselection + N(B-tag jets) = 0]"
    plots[-1].rebin = mt_rebin
    plots[-1].ymax = 4e4
    plots[-1].ymin = 200
    plots[-1].xmax = 190
    plots[-1].xmin = 40
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MTenu_noBtaggedJets_MT200To400_PAS"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "200 < M_{T} (e, PFMET) < 400 (GeV) [Preselection + N(B-tag jets)=0]"
    plots[-1].rebin = mt_rebin
    plots[-1].ymax = 4e4
    plots[-1].ymin = 1
    plots[-1].xmax = 400
    plots[-1].xmin = 200
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MTenu_gteOneBtaggedJet_MT200To400_PAS"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "200 < M_{T} (e, PFMET) < 400 (GeV) [Preselection + N(B-tag jets) #geq 1]"
    plots[-1].rebin = mt_rebin
    plots[-1].ymax = 4e4
    plots[-1].ymin = 1
    plots[-1].xmax = 400
    plots[-1].xmin = 200
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MTenu_50_110_Mjj50to110"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "50 < M_{T} (e, PFMET) < 110 (GeV) [Preselection + 50 < M_{jj} < 110 (GeV)]"
    plots[-1].rebin = mt_rebin
    plots[-1].ymax = 4e4
    plots[-1].ymin = 200
    plots[-1].xmax = 190
    plots[-1].xmin = 40
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MTenu_50_110_Mjj50to110_addBtagJet"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "50 < M_{T} (e, PFMET) < 110 (GeV) [Preselection + 50 < M_{jj} < 110 (GeV) + add b-tagged jet]"
    plots[-1].rebin = mt_rebin
    plots[-1].ymax = 4e4
    plots[-1].ymin = 200
    plots[-1].xmax = 145
    plots[-1].xmin = 40
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MTenu_50_110_Mjj50to110_noAddBtagJets"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "50 < M_{T} (e, PFMET) < 110 (GeV) [Preselection + 50 < M_{jj} < 110 (GeV) + no addl b-tag jets]"
    plots[-1].rebin = mt_rebin
    plots[-1].ymax = 4e4
    plots[-1].ymin = 200
    plots[-1].xmax = 145
    plots[-1].xmin = 40
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "MTenu_50_110_MjjGte110"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "50 < M_{T} (e, PFMET) < 110 (GeV) [Preselection + M_{jj} > 110 (GeV)]"
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
    plots[-1].ymax = 2e5
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "S_{T} (e, PFMET) (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "sTjet_PAS"             ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = st_rebin
    plots[-1].ymax = 2e5
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "S_{T} (1st Jet, 2nd Jet) (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "sT_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAboveSt) )
    #plots.append ( makeDefaultPlot ( "ST"                ,  histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = "var"
    plots[-1].ymax = 1e6
    plots[-1].ymin = 1e-1
    plots[-1].xbins = [ 0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400, 420, 440, 460, 480, 500, 520, 540, 560, 580, 600, 620, 640, 660, 680, 700, 720, 740, 760, 780, 800, 820, 840, 860, 880, 900, 920, 940, 960, 980, 1000, 1020, 1040, 1060, 1080, 1100, 1140, 1180, 1220, 1260, 1300, 1400, 1500, 1600, 1700, 1800, 2000 ]
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "S_{T} (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "Mjj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
    plots[-1].rebin = mass_rebin
    plots[-1].ymax = 2e5
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "M(jj) (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "Mej1_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    plots[-1].rebin = mass_rebin
    plots[-1].ymax = 2000000
    plots[-1].ymin = 1e-1
    plots[-1].xmin  = 0
    plots[-1].xmax  = 2000
    plots[-1].rebin = 2
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "Mass (1st Electron, 1st Jet) (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "Mej2_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    plots[-1].rebin = mass_rebin
    plots[-1].ymax = 2000000
    plots[-1].ymin = 1e-1
    plots[-1].xmin  = 0
    plots[-1].xmax  = 2000
    plots[-1].rebin = 2
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "Mass (1st Electron, 2nd Jet) (GeV) [Preselection]"
    
    #plots.append ( makeDefaultPlot ( "Mej_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    plots.append ( makeDefaultPlot ( "Mej_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = mass_rebin
    plots[-1].ymax = 2000000
    plots[-1].ymin = 1e-1
    plots[-1].xmin  = 0
    plots[-1].xmax  = 2000
    plots[-1].rebin = mej_rebin
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "M(ej) (GeV) [Preselection]"
    
    plots.append ( makeDefaultPlot ( "Mej_Barrel_Presel"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = mass_rebin
    plots[-1].ymax = 2000000
    plots[-1].ymin = 1e-1
    plots[-1].xmin  = 0
    plots[-1].xmax  = 2000
    plots[-1].rebin = mej_rebin
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "M(ej) (GeV) [Barrel, Preselection]"
    
    plots.append ( makeDefaultPlot ( "Mej_Endcap1_Presel"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = mass_rebin
    plots[-1].ymax = 2000000
    plots[-1].ymin = 1e-1
    plots[-1].xmin  = 0
    plots[-1].xmax  = 2000
    plots[-1].rebin = mej_rebin
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "M(ej) (GeV) [Endcap1, Preselection]"
    
    plots.append ( makeDefaultPlot ( "Mej_Endcap2_Presel"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = mass_rebin
    plots[-1].ymax = 2000000
    plots[-1].ymin = 1e-1
    plots[-1].xmin  = 0
    plots[-1].xmax  = 2000
    plots[-1].rebin = mej_rebin
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "M(ej) (GeV) [Endcap2, Preselection]"
    
    #plots.append ( makeDefaultPlot ( "MejGte1500_Pt1stEle_Barrel_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].rebin = pt_rebin
    #plots[-1].ymax = 2000000
    #plots[-1].ymin = 1e-1
    #plots[-1].xmin  = 0
    #plots[-1].xmax  = 2000
    #plots[-1].rebin = pt_rebin
    #plots[-1].ylog  = "yes"
    #plots[-1].xtit = "P_{T} Ele (GeV) [M_{ej} > 1500 GeV, barrel, Preselection]"
    #
    #plots.append ( makeDefaultPlot ( "MejGte1500_Pt1stEle_Endcap1_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].rebin = pt_rebin
    #plots[-1].ymax = 2000000
    #plots[-1].ymin = 1e-1
    #plots[-1].xmin  = 0
    #plots[-1].xmax  = 2000
    #plots[-1].rebin = pt_rebin
    #plots[-1].ylog  = "yes"
    #plots[-1].xtit = "P_{T} Ele (GeV) [M_{ej} > 1500 GeV, endcap1, Preselection]"
    #
    #plots.append ( makeDefaultPlot ( "MejGte1500_Pt1stEle_Endcap2_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].rebin = pt_rebin
    #plots[-1].ymax = 2000000
    #plots[-1].ymin = 1e-1
    #plots[-1].xmin  = 0
    #plots[-1].xmax  = 2000
    #plots[-1].rebin = pt_rebin
    #plots[-1].ylog  = "yes"
    #plots[-1].xtit = "P_{T} Ele (GeV) [M_{ej} > 1500 GeV, endcap2, Preselection]"
    #
    #plots.append ( makeDefaultPlot ( "MejGte1500_PtJet_Barrel_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].rebin = pt_rebin
    #plots[-1].ymax = 2000000
    #plots[-1].ymin = 1e-1
    #plots[-1].xmin  = 0
    #plots[-1].xmax  = 2000
    #plots[-1].rebin = pt_rebin
    #plots[-1].ylog  = "yes"
    #plots[-1].xtit = "P_{T} Jet (GeV) [M_{ej} > 1500 GeV, barrel, Preselection]"
    #
    #plots.append ( makeDefaultPlot ( "MejGte1500_PtJet_Endcap1_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].rebin = pt_rebin
    #plots[-1].ymax = 2000000
    #plots[-1].ymin = 1e-1
    #plots[-1].xmin  = 0
    #plots[-1].xmax  = 2000
    #plots[-1].rebin = pt_rebin
    #plots[-1].ylog  = "yes"
    #plots[-1].xtit = "P_{T} Jet (GeV) [M_{ej} > 1500 GeV, endcap1, Preselection]"
    #
    #plots.append ( makeDefaultPlot ( "MejGte1500_PtJet_Endcap2_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].rebin = pt_rebin
    #plots[-1].ymax = 2000000
    #plots[-1].ymin = 1e-1
    #plots[-1].xmin  = 0
    #plots[-1].xmax  = 2000
    #plots[-1].rebin = pt_rebin
    #plots[-1].ylog  = "yes"
    #plots[-1].xtit = "P_{T} Jet (GeV) [M_{ej} > 1500 GeV, endcap2, Preselection]"
    #
    #plots.append ( makeDefaultPlot (  "MejGte1500_DeltaPhiEleMET_Barrel_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].xtit = "#Delta#phi( 1st Electron, PFMET ) [M_{ej} > 1500 GeV, barrel, Preselection]"
    #plots[-1].rebin = dphi_rebin
    #plots[-1].ymax = 20
    #plots[-1].ymin = 0
    #
    #plots.append ( makeDefaultPlot (  "MejGte1500_DeltaPhiEleMET_Endcap1_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].xtit = "#Delta#phi( 1st Electron, PFMET ) [M_{ej} > 1500 GeV, endcap1, Preselection]"
    #plots[-1].rebin = dphi_rebin
    #plots[-1].ymax = 20
    #plots[-1].ymin = 0
    #
    #plots.append ( makeDefaultPlot (  "MejGte1500_DeltaPhiEleMET_Endcap2_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].xtit = "#Delta#phi( 1st Electron, PFMET ) [M_{ej} > 1500 GeV, endcap2, Preselection]"
    #plots[-1].rebin = dphi_rebin
    #plots[-1].ymax = 20
    #plots[-1].ymin = 0
    #
    #plots.append ( makeDefaultPlot ( "MejGte1500_DREleJet_Barrel_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].rebin = dr_rebin
    #plots[-1].xtit = "#DeltaR(e_{1},j) [M_{ej} > 1500 GeV, barrel, Preselection]"
    #plots[-1].ymax = 1e5
    #plots[-1].ymin = 1e-1
    #plots[-1].ylog  = "yes"
    #
    #plots.append ( makeDefaultPlot ( "MejGte1500_DREleJet_Endcap1_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].rebin = dr_rebin
    #plots[-1].xtit = "#DeltaR(e_{1},j) [M_{ej} > 1500 GeV, endcap1, Preselection]"
    #plots[-1].ymax = 1e5
    #plots[-1].ymin = 1e-1
    #plots[-1].ylog  = "yes"
    #
    #plots.append ( makeDefaultPlot ( "MejGte1500_DREleJet_Endcap2_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].rebin = dr_rebin
    #plots[-1].xtit = "#DeltaR(e_{1},j) [M_{ej} > 1500 GeV, endcap2, Preselection]"
    #plots[-1].ymax = 1e5
    #plots[-1].ymin = 1e-1
    #plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "Mee_allElectrons_3EleEvents_Presel"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = mass_rebin
    plots[-1].ymax = 2000000
    plots[-1].ymin = 1e-1
    plots[-1].xmin  = 0
    plots[-1].xmax  = 200
    plots[-1].rebin = mass_rebin
    plots[-1].ylog  = "yes"
    plots[-1].xtit = "M(ee) (GeV) [3ele events, Preselection]"
    
    plots.append ( makeDefaultPlot ( "MTjnu_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio, dataBlindAbovePt1) )
    #plots[-1].rebin = mass_rebin
    plots[-1].rebin = 2
    #plots[-1].rebin = 'var'
    #plots[-1].xbins = met_bins
    plots[-1].ymax = 2000000
    plots[-1].ymin = 1e-1
    plots[-1].xmin  = 0
    plots[-1].xmax  = 1000
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
    #plots[-1].ymax = 2.5e4 # linear
    plots[-1].ymax = 1e8
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot (  "mDPhi1stJetMET_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "#Delta#phi( 1st Jet, PFMET ) [Preselection]"
    plots[-1].rebin = dphi_rebin
    #plots[-1].ymax = 3e4 # linear
    plots[-1].ymax = 1e10
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot (  "mDPhi2ndJetMET_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "#Delta#phi( 2nd Jet, PFMET ) [Preselection]"
    plots[-1].rebin = dphi_rebin
    #plots[-1].ymax = 1.5e4 # linear
    plots[-1].ymax = 1e10
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    
    plots.append ( makeDefaultPlot ( "minDR_EleJet_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = dr_rebin
    plots[-1].xtit = "Minimum #DeltaR(e_{1},(j_{1}, j_{2}, j_{3})) [Preselection]"
    plots[-1].ymax = 1e5
    plots[-1].ymin = 1e-1
    plots[-1].ylog  = "yes"
    
    #plots.append ( makeDefaultPlot ( "MejGte1500_minDR_EleJet_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].rebin = dr_rebin
    #plots[-1].xtit = "Minimum #DeltaR(e_{1},(j_{1}, j_{2}, j_{3})) [M_{ej}>1500+Preselection]"
    #plots[-1].ymax = 1e5
    #plots[-1].ymin = 1e-1
    #plots[-1].ylog  = "yes"
    #
    #plots.append ( makeDefaultPlot ( "MejGte1500_minDR_EleJet_Barrel_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].rebin = dr_rebin
    #plots[-1].xtit = "Minimum #DeltaR(e_{1},(j_{1}, j_{2}, j_{3})) [barrel,M_{ej}>1500+Preselection]"
    #plots[-1].ymax = 1e5
    #plots[-1].ymin = 1e-1
    #plots[-1].ylog  = "yes"
    #
    #plots.append ( makeDefaultPlot ( "MejGte1500_minDR_EleJet_Endcap1_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].rebin = dr_rebin
    #plots[-1].xtit = "Minimum #DeltaR(e_{1},(j_{1}, j_{2}, j_{3})) [endcap1,M_{ej}>1500+Preselection]"
    #plots[-1].ymax = 1e5
    #plots[-1].ymin = 1e-1
    #plots[-1].ylog  = "yes"
    #
    #plots.append ( makeDefaultPlot ( "MejGte1500_minDR_EleJet_Endcap2_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    #plots[-1].rebin = dr_rebin
    #plots[-1].xtit = "Minimum #DeltaR(e_{1},(j_{1}, j_{2}, j_{3})) [endcap2,M_{ej}>1500+Preselection]"
    #plots[-1].ymax = 1e5
    #plots[-1].ymin = 1e-1
    #plots[-1].ylog  = "yes"
    
    
    plots.append ( makeDefaultPlot ( "nVertex_PAS",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit = "n(vertex) [Preselection]"
    plots[-1].ymax = 2e8
    plots[-1].ymin = 1e-1
    plots[-1].xmax = 60.5
    plots[-1].xmin = -0.5
    plots[-1].rebin = 1
    plots[-1].ylog  = "yes"


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
    
        plots.append ( makeDefaultPlot ( "Pt1stMuon_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "1st Muon p_{T} (GeV) [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 2
    
        plots.append ( makeDefaultPlot ( "Eta1stEle_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "1st Electron #eta [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 2
    
        plots.append ( makeDefaultPlot ( "Phi1stEle_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "1st Electron #phi [LQ M = " + str ( mass_point ) + " selection]"
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
    
        plots.append ( makeDefaultPlot ( "MET_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit  = "PFMET [LQ M = " + str ( mass_point ) + " selection]"
        plots[-1].rebin = 10
        
        
        mej_rebin = 4
        mej_xmin = 0
        mej_xmax = 2300
        mee_rebin = 12
        st_rebin  = 10
        dr_rebin  = 2
        plots.append ( makeDefaultPlot ( "Mej_LQ" + str ( mass_point ) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].rebin = mej_rebin
        plots[-1].xtit = "M_{ej} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].xmin = mej_xmin
        plots[-1].xmax = mej_xmax
    
        plots.append ( makeDefaultPlot ( "ST_LQ" + str( mass_point ), histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].xtit = "S_{T} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        plots[-1].rebin = st_rebin
        
        plots.append ( makeDefaultPlot ( "MTenu_LQ" + str ( mass_point ) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        plots[-1].rebin = mee_rebin
        plots[-1].xtit = "M_{T}(e#nu) (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        
    
    #FIXME ALL
        #plots.append ( makeDefaultPlot ( "sTlep_LQ" + str( mass_point ), histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        #plots[-1].xtit = "S_{T} (1st Electron, 2nd Electron) (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].rebin = 10
    
        #plots.append ( makeDefaultPlot ( "sTjet_LQ" + str( mass_point ), histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        #plots[-1].xtit = "S_{T} (1st Jet, 2nd Jet) (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].rebin = 10
    
        #plots.append ( makeDefaultPlot ( "sTfrac_Jet1_LQ" + str( mass_point )                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        #plots[-1].xmin = 0.0
        #plots[-1].xmax = 1.0
        #plots[-1].rebin = 4
        #plots[-1].xtit = "Fraction S_{T} from jet 1 (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        #
        #plots.append ( makeDefaultPlot ( "sTfrac_Jet2_LQ" + str( mass_point )               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        #plots[-1].xmin = 0.0
        #plots[-1].xmax = 1.0
        #plots[-1].rebin = 4
        #plots[-1].xtit = "Fraction S_{T} from jet 2 (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        #
        #plots.append ( makeDefaultPlot ( "sTfrac_Ele1_LQ" + str( mass_point )                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        #plots[-1].xmin = 0.0
        #plots[-1].xmax = 1.0
        #plots[-1].rebin = 4
        #plots[-1].xtit = "Fraction S_{T} from ele 1 (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        #
        #plots.append ( makeDefaultPlot ( "sTfrac_Ele_LQ" + str( mass_point )                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        #plots[-1].xmin = 0.0
        #plots[-1].xmax = 1.0
        #plots[-1].rebin = 4
        #plots[-1].xtit = "Fraction S_{T} from electrons (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        #
        #plots.append ( makeDefaultPlot ( "sTfrac_Jet_LQ" + str (mass_point) ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        #plots[-1].xmin = 0.0
        #plots[-1].xmax = 1.0
        #plots[-1].rebin = 4
        #plots[-1].xtit = "Fraction S_{T} from jets (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        #
        #plots.append ( makeDefaultPlot ( "sTfrac_MET_LQ" + str (mass_point) ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        #plots[-1].xmin = 0.0
        #plots[-1].xmax = 1.0
        #plots[-1].rebin = 4
        #plots[-1].xtit = "Fraction S_{T} from MET (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        
        #plots.append ( makeDefaultPlot ( "Mjj_LQ" + str( mass_point )          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
        #plots[-1].rebin = 4
        #plots[-1].xtit = "Dijet Mass (GeV) (LQ M = " + str ( mass_point ) + " selection)"
        #
        #plots.append ( makeDefaultPlot ( "Ptj1j2_LQ" + str( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        #plots[-1].rebin = 2
        #plots[-1].xmin  = 0
        #plots[-1].xmax  = 1000
        #plots[-1].xtit  = "P_{T}(j_{1}, j_{2}) (GeV), (LQ M = " + str ( mass_point ) + " selection)"
    
        #plots.append ( makeDefaultPlot ( "nVertex_LQ" + str( mass_point )                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        #plots[-1].xmin = -0.5
        #plots[-1].xmax = 40.5
        #plots[-1].xtit = "n(vertexes) [Preselection]"
        
        #plots.append ( makeDefaultPlot ( "DR_Ele1Jet1_LQ" + str ( mass_point ) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        #plots[-1].rebin = dr_rebin
        #plots[-1].xtit = "#DeltaR(e1,j1) (GeV), (LQ M = " + str ( mass_point ) + " selection)"
        #
        
        #plots.append ( makeDefaultPlot ("Classif_1stEle_LQ" + str ( int ( mass_point )  )               , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron classification, (LQ M = " + str ( mass_point ) + " selection)"
        #
        #
        #plots.append ( makeDefaultPlot ("CorrIsolation_1stEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron corrected HEEP isolation, (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].xmax  = 1.0
        #plots[-1].xmin  = -25
        #plots[-1].rebin = 4
        #
        #plots.append ( makeDefaultPlot ("DeltaEtaTrkSC_1stEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        #plots[-1].xtit  = "1st Electron #Delta#eta(track, SC), (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].xmin  = -0.007
        #plots[-1].xmax  =  0.007
        #plots[-1].rebin = 4
        #
        #plots.append ( makeDefaultPlot ("DeltaPhiTrkSC_1stEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        #plots[-1].xtit  = "1st Electron #Delta#phi(track, SC), (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].xmin  = -0.06
        #plots[-1].xmax  =  0.06
        #plots[-1].rebin = 4
        #
        #plots.append ( makeDefaultPlot ("Full5x5E1x5OverE5x5_1stEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        #plots[-1].xtit  = "1st Electron E_{1x5} / E_{5x5}, (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].xmin  =  0.20
        #plots[-1].xmax  =  1.40
        #plots[-1].rebin = 4
        #
        #plots.append ( makeDefaultPlot ("Full5x5E2x5OverE5x5_1stEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        #plots[-1].xtit  = "1st Electron E_{2x5} / E_{5x5}, (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].xmin  =  0.20
        #plots[-1].xmax  =  1.40
        #plots[-1].rebin = 4
        #
        #plots.append ( makeDefaultPlot ("EcalIsolation_1stEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        #plots[-1].xtit  = "1st Electron HEEP ECAL isolation, (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].xmin  =   0.00
        #plots[-1].xmax  =  15.00
        #plots[-1].rebin = 4
        #
        #plots.append ( makeDefaultPlot ("HcalIsolation_1stEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        #plots[-1].xtit  = "1st Electron HEEP HCAL isolation, (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].xmin  =   0.00
        #plots[-1].xmax  =  15.00
        #plots[-1].rebin = 4
        
        #plots.append ( makeDefaultPlot ("TrkIsolation_1stEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
        #plots[-1].xtit  = "1st Electron HEEP tracker isolation, (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].rebin = 4
        #
        #plots.append ( makeDefaultPlot ("FBrem_1stEle_LQ" + str ( int ( mass_point )  )                 , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron brem fraction, (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].xmin  =   -5.0
        #plots[-1].xmax  =    4.0
        #
        #plots.append ( makeDefaultPlot ("GsfCtfCharge_1stEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron GSF CTF charge, (LQ M = " + str ( mass_point ) + " selection)"
        #
        #plots.append ( makeDefaultPlot ("GsfCtfScPixCharge_1stEle_LQ" + str ( int ( mass_point )  )     , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron GSF CTF SC Pixel charge, (LQ M = " + str ( mass_point ) + " selection)"
        #
        #plots.append ( makeDefaultPlot ("GsfScPixCharge_1stEle_LQ" + str ( int ( mass_point )  )        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron GSF SC Pixel charge, (LQ M = " + str ( mass_point ) + " selection)"
        #
        #plots.append ( makeDefaultPlot ("HasMatchedPhot_1stEle_LQ" + str ( int ( mass_point )  )        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron has matched photon, (LQ M = " + str ( mass_point ) + " selection)"
        #
        #plots.append ( makeDefaultPlot ("HoE_1stEle_LQ" + str ( int ( mass_point )  )                   , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron H/E, (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].rebin = 4
        #
        #plots.append ( makeDefaultPlot ("LeadVtxDistXY_1stEle_LQ" + str ( int ( mass_point )  )         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron D_{XY} vs leading vertex, (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].rebin = 4
        #
        #plots.append ( makeDefaultPlot ("LeadVtxDistZ_1stEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron D_{Z} vs leading vertex, (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].rebin = 4
        #
        #plots.append ( makeDefaultPlot ("MissingHits_1stEle_LQ" + str ( int ( mass_point )  )           , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron N(missing hits), (LQ M = " + str ( mass_point ) + " selection)"
        #
        #plots.append ( makeDefaultPlot ("NBrems_1stEle_LQ" + str ( int ( mass_point )  )                , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron N(brems), (LQ M = " + str ( mass_point ) + " selection)"
        #
        #plots.append ( makeDefaultPlot ("EnergyORawEnergy_1stEle_LQ" + str ( int ( mass_point )  )      , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron energy correction factor, (LQ M = " + str ( mass_point ) + " selection)" 
        #plots[-1].xmin  = 0.95
        #plots[-1].xmax  = 1.3
        #plots[-1].rebin = 4
        #
        #plots.append ( makeDefaultPlot ("SigmaEtaEta_Barrel_1stEle_LQ" + str ( int ( mass_point )  )    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron #sigma_{#eta#eta} [Preselection, Barrel]"
        #plots[-1].xmin  = 0.005
        #plots[-1].xmax  = 0.015
        #plots[-1].rebin = 2
        #
        #plots.append ( makeDefaultPlot ("SigmaEtaEta_Endcap_1stEle_LQ" + str ( int ( mass_point )  )    , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron #sigma_{#eta#eta} [Preselection, Endcap]"
        #plots[-1].xmin  = 0.005
        #plots[-1].xmax  = 0.04
        #plots[-1].rebin = 2
        #
        ##plots.append ( makeDefaultPlot ("Full5x5SigmaIEtaIEta_Barrel_1stEle_LQ" + str ( int ( mass_point )  )  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        ##plots[-1].xtit  = "1st Electron #sigma_{i#etai#eta} [Preselection, Barrel]"
        ##plots[-1].xmin  = 0.005
        ##plots[-1].xmax  = 0.015
        ##plots[-1].rebin = 2
        ##
        ##plots.append ( makeDefaultPlot ("Full5x5SigmaIEtaIEta_Endcap_1stEle_LQ" + str ( int ( mass_point )  )  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        ##plots[-1].xtit  = "1st Electron #sigma_{i#etai#eta} [Preselection, Endcap]"
        ##plots[-1].xmin  = 0.015
        ##plots[-1].xmax  = 0.035
        ##plots[-1].rebin = 2
        ##
        #plots.append ( makeDefaultPlot ("TrkPtOPt_1stEle_LQ" + str ( int ( mass_point )  )              , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron track pt / SC pt, (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].rebin = 4
        #
        #plots.append ( makeDefaultPlot ("ValidFrac_1stEle_LQ" + str ( int ( mass_point )  )             , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
        #plots[-1].xtit  = "1st Electron valid fraction of hits, (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].xmin  = 0.00
        #plots[-1].xmax  = 1.50
        #plots[-1].rebin = 4
        #
        #plots.append ( makeDefaultPlot ( "EleChargeSum_LQ" + str( int (mass_point))        , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
        #plots[-1].xtit  = "Electron 1 charge + Electron 2 charge, (LQ M = " + str ( mass_point ) + " selection)"
        #plots[-1].xmin  = -1.5
        #plots[-1].xmax  = 1.5

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

