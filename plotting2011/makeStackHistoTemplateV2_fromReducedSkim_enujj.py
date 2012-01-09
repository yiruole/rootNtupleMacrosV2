#!/usr/bin/env python

from plot_class import *
from ROOT import *

File_preselection = GetFile(os.environ['LQDATA']+"/enujj_analysis/enujj/output_cutTable_lq_enujj/analysisClass_lq_enujj_plots.root")

File_QCD = GetFile(os.environ['LQDATA']+"/enujj_analysis/enujj_qcd/output_cutTable_lq_enujj/analysisClass_lq_enujj_QCD_plots.root")

#### Common values for plots:
#otherBkgsKey="QCD, single top, VV+jets, W/W*+jets"
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

QCDScale = 1.0

#--- Final plots are defined here

# Simply add or remove plots and update the list plots = [plot0, plot1, ...] below

histoBaseName = "histo1D__SAMPLE__cutHisto_allOtherCuts___________VARIABLE"
histoBaseName_userDef = "histo1D__SAMPLE__VARIABLE"


samplesForStackHistosQCD = ["DATA"]

samplesForStackHistos = ["ZJet_Madgraph","PhotonJets","SingleTop","DIBOSON","TTbar_Madgraph","WJet_Madgraph"]
keysStack =             ["QCD multijets", "Z/Z* + jets"  ,"#gamma + jets","single top","WW + WZ + ZZ","t#bar{t}", "W/W* + jets"]
stackColorIndexes     = [7   , 6              , 14           ,  4            , 3             , 92             , 2             ]
stackFillStyleIds     = [3345, 3345           , 3344         ,  3345         , 3354          , 3354           , 3395          ]

stackColorIndexes.reverse()
stackFillStyleIds.reverse()

samplesForHistos = ["LQ_M450"      , "LQ_M550"      , "LQ_M650"      ]
keys             = ["LQ, M=450 GeV", "LQ, M=550 GeV", "LQ, M=650 GeV"]

sampleForDataHisto = "DATA"

def makeDefaultPlot ( variableName, histoBaseName, 
                      samplesForHistos, keys,
                      samplesForStackHistos, keysStack,
                      sampleForDataHisto,
                      zUncBand, makeRatio ) : 
    plot                = Plot() 
    plot.histosStack    = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDScale) + generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection) 
    plot.keysStack      = keysStack
    plot.histos         = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
    plot.keys           = keys
    plot.addZUncBand    = zUncBand
    plot.makeRatio      = makeRatio
    plot.makeNSigma     = makeNSigma
    plot.histodata      = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.ytit           = "Events"
    plot.ylog           = "no"
    plot.name           = variableName
    plot.stackColorIndexes = stackColorIndexes
    plot.stackFillStyleIds = stackFillStyleIds 
    plot.gif_folder     = "gif_enujj/"
    plot.eps_folder     = "eps_enujj/"
    plot.suffix        = "enujj"
    plot.lumi_pb        = "4623"

    return plot

plots = []

plots.append ( makeDefaultPlot ( "nEle", histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit            = "Number of electrons (cut)"
plots[-1].ylog            = "yes"
plots[-1].xmin            = -0.5
plots[-1].xmax            = 6.5
plots[-1].ymin            = 0.0001
plots[-1].ymax            = 10000000000

plots.append ( makeDefaultPlot ( "Ele1_Pt"  , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax  = 200000000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit  = "1st Electron p_{T} [GeV] (cut)"

plots.append ( makeDefaultPlot ( "Ele1_Eta" , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = eta_rebin
plots[-1].ymax  = 20000000000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit  = "1st Electron #eta (cut)"

plots.append ( makeDefaultPlot ( "MET"    , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax  = 10000000
plots[-1].ymin  = 1e-1
plots[-1].xmax  = 500.

plots[-1].ylog  = "yes"
plots[-1].xtit  = "PFMET [GeV] (cut)"

plots.append ( makeDefaultPlot ( "nJet" , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit  = "Number of jets"
plots[-1].ymax  = 200000000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin  = -0.5
plots[-1].xmax  = 10.5

plots.append ( makeDefaultPlot ( "Jet1_Pt"     , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit  = "1st Jet p_{T} [GeV] (cut)"

plots.append ( makeDefaultPlot ( "Jet1_Eta"    , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = eta_rebin
plots[-1].ymax  = 4000000000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit  = "1st Jet #eta (cut)"

plots.append ( makeDefaultPlot ( "Jet2_Pt"     , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit  = "1st Jet p_{T} [GeV] (cut)"

plots.append ( makeDefaultPlot ( "Jet2_Eta"    , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = eta_rebin
plots[-1].ymax  = 4000000000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit  = "1st Jet #eta (cut)"

plots.append ( makeDefaultPlot ( "ST"    , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = st_rebin
plots[-1].ymax  = 100000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 200.0
plots[-1].xmax  = 1700.
plots[-1].ylog  = "yes"
plots[-1].xtit  = "ST (1st Electron, MET, 1st Jet, 2nd Jet) [GeV] (cut)"

plots.append ( makeDefaultPlot ( "nMuon"      , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax  = 10000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = -0.5
plots[-1].xmax  = 6.5
plots[-1].ylog  = "yes"
plots[-1].xtit  = "Number of muons (cut)"

plots.append ( makeDefaultPlot ( "DR_Ele1Jet1"                ,  histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = dr_rebin
plots[-1].xtit = "#DeltaR(e_{1},j_{1}) (cut)"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "DR_Ele1Jet2"                ,  histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = dr_rebin
plots[-1].xtit = "#DeltaR(e_{1},j_{2}) (cut)"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "mDeltaPhiMETEle"   , 	  histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "#Delta#phi(e_{1},MET) (cut)"
plots[-1].rebin = dphi_rebin
plots[-1].ymax = 2000
plots[-1].ymin = 1e-1

plots.append ( makeDefaultPlot ( "mDeltaPhiMET1stJet",    histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "#Delta#phi(j_{1},MET) (cut)"
plots[-1].rebin = dphi_rebin
plots[-1].ymax = 2000
plots[-1].ymin = 1e-1

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

plots.append ( makeDefaultPlot ( "nElectron_PAS", histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit            = "Number of electrons (preselection)"
plots[-1].ylog            = "yes"
plots[-1].rebin           = 1
plots[-1].xmin            = -0.5
plots[-1].xmax            = 6.5
plots[-1].ymin            = 0.0001
plots[-1].ymax            = 10000000000

plots.append ( makeDefaultPlot ( "nMuon_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "Number of muons (Preselection)"
plots[-1].ylog            = "yes"
plots[-1].rebin           = 1
plots[-1].xmin            = -0.5
plots[-1].xmax            = 6.5
plots[-1].ymin            = 0.0001
plots[-1].ymax            = 10000000000

plots.append ( makeDefaultPlot ( "Pt1stEle_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Electron p_{T} (Preselection) [GeV]"
plots[-1].xmin = 0.
plots[-1].xmax = 600.
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].rebin = pt_rebin

plots.append ( makeDefaultPlot ( "Eta1stEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Electron #eta (Preselection)"   
plots[-1].ymax = 20000000
plots[-1].rebin = eta_rebin
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Phi1stEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )   
plots[-1].xtit = "1st Electron #phi (Preselection)"
plots[-1].rebin = 1
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Charge1stEle_PAS"      ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit = "1st Electron Charge (Preselection)"
plots[-1].ylog  = "yes"
plots[-1].ymin = 1e-1
plots[-1].ymax = 10000000000

plots.append ( makeDefaultPlot ( "MatchPhotonConv1stEle_PAS"      ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit = "1st Electron matched to a photon (Preselection)"
plots[-1].ylog  = "yes"
plots[-1].ymin = 1e-1
plots[-1].ymax = 10000000000

plots.append ( makeDefaultPlot ( "MET_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET (Preselection) [GeV]"
plots[-1].rebin = pt_rebin
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "METSig_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET Significance (Preselection) [GeV]"
plots[-1].rebin = pt_rebin
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 700
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "METPhi_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET #phi (Preselection)"
plots[-1].rebin = 1
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "METCharged_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET Charged (Preselection) [GeV]"
plots[-1].rebin = pt_rebin
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "METChargedPhi_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET Charged #phi (Preselection)"
plots[-1].rebin = dphi_rebin
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "METType1_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET Type1-Corrected (Preselection) [GeV]"
plots[-1].rebin = pt_rebin
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "METType1Phi_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET Type1-Corrected #phi (Preselection)"
plots[-1].rebin = 1
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "minMETPt1stEle_PAS"    ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 500
plots[-1].xmin = 0
plots[-1].ylog  = "yes"
plots[-1].xtit = "Min (PFMET, 1st Electron p_{T}) (Preselection) [GeV]"

plots.append ( makeDefaultPlot ( "Pt1stJet_PAS"          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Jet p_{T} (Preselection) [GeV]"
plots[-1].rebin = pt_rebin
plots[-1].ymax = 200000
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
plots[-1].xtit = "2nd Jet p_{T} (Preselection) [GeV]"

plots.append ( makeDefaultPlot ( "Eta1stJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = eta_rebin
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "1st Jet #eta (Preselection)"

plots.append ( makeDefaultPlot ( "Eta2ndJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = eta_rebin
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet #eta (Preselection)"

plots.append ( makeDefaultPlot ( "Phi1stJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "1st Jet #phi (Preselection)"

plots.append ( makeDefaultPlot ( "Phi2ndJet_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet #phi (Preselection)"

plots.append ( makeDefaultPlot ( "TCHE1stJet_PAS"        ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "1st Jet TCHE (Preselection)"
plots.append ( makeDefaultPlot ( "TCHE2ndJet_PAS"        ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet TCHE (Preselection)"

plots.append ( makeDefaultPlot ( "MT_charged_enu_PAS",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (1st Electron, Charged PFMET) (Preselection) [GeV]"
plots[-1].rebin = mt_rebin
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "MT_type1_enu_PAS"  ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (1st Electron, Type 1 PFMET) (Preselection) [GeV]"
plots[-1].rebin = mt_rebin
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "MTenu_PAS"             ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (1st Electron, PFMET) (Preselection) [GeV]"
plots[-1].rebin = mt_rebin
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "MTenu_50_110"            ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "50 < M_{T} (1st Electron, PFMET) (Preselection) < 100 [GeV]"
plots[-1].rebin = mt_rebin
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].xmax = 140
plots[-1].xmin = 40
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Ptenu_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"
plots[-1].xtit = "p_{T} (1st Electron, PFMET) (Preselection) [GeV]"

plots.append ( makeDefaultPlot ( "sTlep_PAS"             ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = st_rebin
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (1st Electron, PFMET) (Preselection) [GeV]"

plots.append ( makeDefaultPlot ( "sTjet_PAS"             ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = st_rebin
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (1st Jet, 2nd Jet) (Preselection) [GeV]"

plots.append ( makeDefaultPlot ( "sT_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = st_rebin
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (1st Electron, 1st Jet, 2nd Jet, PFMET) (Preselection) [GeV]"

plots.append ( makeDefaultPlot ( "Mjj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = mass_rebin
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "Dijet Mass (Preselection) [GeV]"

plots.append ( makeDefaultPlot ( "Mej1_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = mass_rebin
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "Mass (1st Electron, 1st Jet) (Preselection) [GeV]"

plots.append ( makeDefaultPlot ( "Mej2_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = mass_rebin
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "Mass (1st Electron, 2nd Jet) (Preselection) [GeV]"

plots.append ( makeDefaultPlot (  "DCotTheta1stEle_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Electron D #times Cot(theta) [cm] (Preselection)" 
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot (  "Dist1stEle_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Electron Distance [cm] (Preselection) " 
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot (  "mDPhi1stEleMET_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "#Delta#phi( 1st Electron, PFMET ) (Preselection)"
plots[-1].rebin = dphi_rebin
plots[-1].ymax = 2000
plots[-1].ymin = 0

plots.append ( makeDefaultPlot (  "mDPhi1stJetMET_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "#Delta#phi( 1st Jet, PFMET ) (Preselection)"
plots[-1].rebin = dphi_rebin
plots[-1].ymax = 2000
plots[-1].ymin = 0

plots.append ( makeDefaultPlot (  "mDPhi2ndJetMET_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "#Delta#phi( 2nd Jet, PFMET ) (Preselection)"
plots[-1].rebin = dphi_rebin
plots[-1].ymax = 2000
plots[-1].ymin = 0

plots.append ( makeDefaultPlot (  "MT_GoodVtxLTE3_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (1st Electron, PFMET) (Preselection + Good vertices #leq 3) [GeV]"
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].rebin = mt_rebin
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot (  "MT_GoodVtxGTE4_LTE8_PAS"  ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (1st Electron, PFMET) (Preselection + 4 #geq Good vertices #leq 8) [GeV]"
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].rebin = mt_rebin
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot (  "MT_GoodVtxGTE9_LTE15_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (1st Electron, PFMET) (Preselection + 9 #geq Good vertices #leq 15) [GeV]"
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].rebin = mt_rebin
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot (  "MT_GoodVtxGTE16_PAS"      ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (1st Electron, PFMET) (Preselection + Good vertices #geq 16) [GeV]"
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].rebin = mt_rebin
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot (  "MTCharged_GoodVtxLTE3_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (1st Electron, Charged PFMET) (Preselection + Good vertices #leq 3) [GeV]"
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].rebin = mt_rebin
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot (  "MTCharged_GoodVtxGTE4_LTE8_PAS"  ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (1st Electron, Charged PFMET) (Preselection + 4 #geq Good vertices #leq 8) [GeV]"
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].rebin = mt_rebin
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot (  "MTCharged_GoodVtxGTE9_LTE15_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (1st Electron, Charged PFMET) (Preselection + 9 #geq Good vertices #leq 15) [GeV]"
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].rebin = mt_rebin
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot (  "MTCharged_GoodVtxGTE16_PAS"      ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (1st Electron, Charged PFMET) (Preselection + Good vertices #geq 16) [GeV]"
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].rebin = mt_rebin
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot (  "MTType1_GoodVtxLTE3_PAS"       ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (1st Electron, Type1 PFMET) (Preselection + Good vertices #leq 3) [GeV]"
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].rebin = mt_rebin
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot (  "MTType1_GoodVtxGTE4_LTE8_PAS"  ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (1st Electron, Type1 PFMET) (Preselection + 4 #geq Good vertices #leq 8) [GeV]"
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].rebin = mt_rebin
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot (  "MTType1_GoodVtxGTE9_LTE15_PAS" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (1st Electron, Type1 PFMET) (Preselection + 9 #geq Good vertices #leq 15) [GeV]"
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].rebin = mt_rebin
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot (  "MTType1_GoodVtxGTE16_PAS"      ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T} (1st Electron, Type1 PFMET) (Preselection + Good vertices #geq 16) [GeV]"
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].rebin = mt_rebin
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "minDR_EleJet_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = dr_rebin
plots[-1].xtit = "Minimum #DeltaR(e_{1},(j_{1}, j_{2}, j_{3})) (cut)"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "nVertex_good_PAS",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "n(good vertexes)"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmax = 30.5
plots[-1].xmin = -0.5
plots[-1].rebin = 1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "nVertex_PAS",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "n(vertex)"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmax = 30.5
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
fileps = "allPlots_enujj_analysis.ps" 

c.Print(fileps + "[")
for plot in plots:
    plot.Draw(fileps)
c.Print(fileps+"]")
os.system('ps2pdf '+fileps)
os.system('rm '+fileps)
