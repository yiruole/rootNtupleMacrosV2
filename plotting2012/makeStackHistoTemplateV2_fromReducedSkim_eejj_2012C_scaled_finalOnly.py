#!/usr/bin/env python

from plot_class import *
from ROOT import *

# mass_points = [ 300,  350,  400, 450, 500, 550,  600,  650, 700,  750,  800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 ]
mass_points = [ 300, 500, 650, 900 ]

File_preselection = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj//scaled_output_cutTable_lq_eejj_2012C/analysisClass_lq_eejj_plots.root" )
File_TTBar        = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj//scaled_output_cutTable_lq_eejj_2012C/analysisClass_lq_eejj_plots.root" )
File_QCD          = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj_qcd//output_cutTable_lq_eejj_2012C/analysisClass_lq_eejj_QCD_plots.root")

#### Common values for plots:
zUncBand="no"
makeRatio=1
makeNSigma=1

pt_rebin = 2

histoBaseName = "histo1D__SAMPLE__cutHisto_allOtherCuts___________VARIABLE"
histoBaseName_userDef = "histo1D__SAMPLE__VARIABLE"

histoBaseName2D = "histo2D__SAMPLE__cutHisto_allOtherCuts___________VARIABLE"
histoBaseName2D_userDef = "histo2D__SAMPLE__VARIABLE"


samplesForStackHistos_other  = [ "OTHERBKG"    ]
samplesForStackHistosQCD     = [ "DATA"        ]
samplesForStackHistosTTBar   = [ "TTbar_Madgraph"]
samplesForStackHistos_ZJets  = [ "ZJet_Madgraph" ]
samplesForStackHistos = samplesForStackHistos_other + samplesForStackHistos_ZJets

keysStack             = [ "QCD multijets", "Other backgrounds", "t#bar{t}"  ,  "Z/#gamma* + jets"  ]
stackColorIndexes     = [ 432            , 9                  , 600         ,  kRed           ]
stackFillStyleIds     = [ 3005           , 3008               , 3005        ,  3004           ]

# keysStack             = [ "QCD multijets", "Other backgrounds", "t#bar{t}"  ,  "Z/#gamma* + jets" ]
# stackColorIndexes     = [ 432            , 9                  , 600         ,  kRed           ]
# stackFillStyleIds     = [ 3005           , 3008               , 3005        ,  3004           ]

stackColorIndexes.reverse()
stackFillStyleIds.reverse()

samplesForHistos = []
keys             = []

sampleForDataHisto = "DATA"

QCDScale   = 1.0
TTBarScale = 1.0

#--- nEle_PtCut_IDISO_noOvrlp ---

def makeDefaultPlot ( variableName, histoBaseName, 
                      samplesForHistos, keys,
                      samplesForStackHistos, keysStack,
                      sampleForDataHisto,
                      zUncBand, makeRatio ) :
    plot                   = Plot() 
    plot.histosStack       = generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDScale) + generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistosTTBar, variableName, File_TTBar, TTBarScale) + generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection) 
    # plot.histosStack       = generateHistoList( histoBaseName, samplesForStackHistosTTBar, variableName, File_TTBar, TTBarScale) + generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistos_Signal, variableName, File_preselection) 
    # plot.histosStack       = generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDScale) +  generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection) 
    plot.keysStack         = keysStack
    plot.histos            = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
    plot.keys              = keys
    plot.addZUncBand       = zUncBand
    plot.makeRatio         = makeRatio
    plot.makeNSigma        = makeNSigma
    plot.histodata         = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.ytit              = "Events"
    plot.ylog              = "no"
    plot.name              = variableName
    plot.stackColorIndexes = stackColorIndexes
    plot.stackFillStyleIds = stackFillStyleIds 
    plot.gif_folder        = "gif_eejj_scaled_finalOnly_2012C/"
    plot.eps_folder        = "eps_eejj_scaled_finalOnly_2012C/"
    plot.png_folder        = "png_eejj_scaled_finalOnly_2012C/"
    plot.pdf_folder        = "pdf_eejj_scaled_finalOnly_2012C/"
    plot.suffix            = "eejj_finalOnly"
    plot.lumi_fb           = "6.88"
    
    return plot

def makeDefaultPlot2D ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2D()
    plot.histosStack = generateHistoList ( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDScale) + generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistosTTBar, variableName, File_TTBar, TTBarScale) + generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection) 
    plot.histodata      = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_eejj_scaled_finalOnly_2012C/"
    plot.eps_folder     = "eps_eejj_scaled_finalOnly_2012C/"
    plot.suffix         = "eejj_2D_finalOnly"


    return plot


def makeDefaultPlot2D_NSigma ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2DNSigma()
    plot.histosStack = generateHistoList ( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDScale) + generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistosTTBar, variableName, File_TTBar, TTBarScale) + generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection) 
    plot.histodata      = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_eejj_scaled_finalOnly_2012C/"
    plot.eps_folder     = "eps_eejj_scaled_finalOnly_2012C/"
    plot.suffix         = "eejj_2DNSigma_finalOnly"


    return plot


def makeDefaultPlot2D_Ratio ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2DRatio()
    plot.histosStack = generateHistoList ( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDScale) + generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistosTTBar, variableName, File_TTBar, TTBarScale) + generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection) 
    plot.histodata      = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_eejj_scaled_finalOnly_2012C/"
    plot.eps_folder     = "eps_eejj_scaled_finalOnly_2012C/"
    plot.suffix         = "eejj_2DRatio_finalOnly"
    
    return plot


plots = []

Mej_selected_avg_max = 20000
Mej_selected_min_max = 20000
Mej_minmax_max = 20000
sT_eejj_max          = 20000000
Mee_max              = 20000
ymin_final_selection = 1e-2

mej_xmin = 0
mej_xmax = 2300

st_xmin = 400
st_xmax = 2500

mej_rebin = 2
mee_rebin = 1
st_rebin  = 1
dr_rebin  = 2


for mass_point in mass_points:
    
    plots.append ( makeDefaultPlot ( "Pt1stEle_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit  = "1st Electron p_{T} (GeV) [LQ M = " + str ( mass_point ) + " selection]"
    plots[-1].rebin = 2

    plots.append ( makeDefaultPlot ( "Pt2ndEle_LQ" + str ( mass_point )  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].xtit  = "2nd Electron p_{T} (GeV) [LQ M = " + str ( mass_point ) + " selection]"
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

    plots.append ( makeDefaultPlot ( "Mej_selected_avg_LQ" + str ( mass_point) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = mej_rebin
    plots[-1].xtit = "M_{ej}^{avg} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
    plots[-1].xmin = mej_xmin
    plots[-1].xmax = mej_xmax

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
    
    plots.append ( makeDefaultPlot ( "Mjj_LQ" + str( mass_point )	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
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
    
    plots.append ( makeDefaultPlot ("E1x5OverE5x5_1stEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
    plots[-1].xtit  = "1st Electron E_{1x5} / E_{5x5}, (LQ M = " + str ( mass_point ) + " selection)"
    plots[-1].xmin  =  0.20
    plots[-1].xmax  =  1.40
    plots[-1].rebin = 4
    
    plots.append ( makeDefaultPlot ("E1x5OverE5x5_2ndEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
    plots[-1].xtit  = "2nd Electron E_{1x5} / E_{5x5}, (LQ M = " + str ( mass_point ) + " selection)"
    plots[-1].xmin  =  0.20
    plots[-1].xmax  =  1.40
    plots[-1].rebin = 4
    
    plots.append ( makeDefaultPlot ("E2x5OverE5x5_1stEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
    plots[-1].xtit  = "1st Electron E_{2x5} / E_{5x5}, (LQ M = " + str ( mass_point ) + " selection)"
    plots[-1].xmin  =  0.20
    plots[-1].xmax  =  1.40
    plots[-1].rebin = 4
    
    plots.append ( makeDefaultPlot ("E2x5OverE5x5_2ndEle_LQ" + str ( int ( mass_point )  )          , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )  
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
    
    plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Barrel_1stEle_LQ" + str ( int ( mass_point )  )  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
    plots[-1].xtit  = "1st Electron #sigma_{i#etai#eta} [Preselection, Barrel]"
    plots[-1].xmin  = 0.005
    plots[-1].xmax  = 0.015
    plots[-1].rebin = 2
    
    plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Barrel_2ndEle_LQ" + str ( int ( mass_point )  )  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
    plots[-1].xtit  = "2nd Electron #sigma_{i#etai#eta} [Preselection, Barrel]"
    plots[-1].xmin  = 0.005
    plots[-1].xmax  = 0.015
    plots[-1].rebin = 2
    
    plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Endcap_1stEle_LQ" + str ( int ( mass_point )  )  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
    plots[-1].xtit  = "1st Electron #sigma_{i#etai#eta} [Preselection, Endcap]"
    plots[-1].xmin  = 0.015
    plots[-1].xmax  = 0.035
    plots[-1].rebin = 2
    
    plots.append ( makeDefaultPlot ("SigmaIEtaIEta_Endcap_2ndEle_LQ" + str ( int ( mass_point )  )  , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
    plots[-1].xtit  = "2nd Electron #sigma_{i#etai#eta} [Preselection, Endcap]"
    plots[-1].xmin  = 0.015
    plots[-1].xmax  = 0.035
    plots[-1].rebin = 2
    
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

fileps = "allPlots_eejj_scaled_2012C_analysis_finalOnly.pdf"

c.Print(fileps + "[")
for i_plot, plot in enumerate(plots):
    plot.Draw(fileps, i_plot + 1 )
c.Print(fileps+"]")


