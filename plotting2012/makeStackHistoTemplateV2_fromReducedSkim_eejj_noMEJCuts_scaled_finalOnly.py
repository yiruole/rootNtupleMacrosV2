#!/usr/bin/env python

from plot_class import *
from ROOT import *

mass_points = [ 300,  350,  400, 450, 500, 550,  600,  650, 700,  750,  800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 ]

File_preselection = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj//scaled_output_cutTable_lq_eejj_noMEJCuts/analysisClass_lq_eejj_plots.root"     )
File_TTBar        = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj//scaled_output_cutTable_lq_eejj_noMEJCuts/analysisClass_lq_eejj_plots.root"     )
File_QCD          = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj_qcd//output_cutTable_lq_eejj_noMEJCuts/analysisClass_lq_eejj_QCD_plots.root"    )

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
samplesForStackHistos_Signal = [ "LQ_M700"     ] 
samplesForStackHistos = samplesForStackHistos_other + samplesForStackHistos_ZJets

keysStack             = [ "QCD multijets", "Other backgrounds", "t#bar{t}"  ,  "Z/#gamma* + jets"  ]
stackColorIndexes     = [ 432            , 9                  , 600         ,  kRed  ]
stackFillStyleIds     = [ 3005           , 3008               , 3005        ,  3004  ]

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
    plot.gif_folder        = "gif_eejj_noMEJCuts_scaled_finalOnly/"
    plot.eps_folder        = "eps_eejj_noMEJCuts_scaled_finalOnly/"
    plot.png_folder        = "png_eejj_noMEJCuts_scaled_finalOnly/"
    plot.pdf_folder        = "pdf_eejj_noMEJCuts_scaled_finalOnly/"
    plot.suffix            = "eejj_noMEJCuts_finalOnly"
    plot.lumi_fb           = "19.5"
    
    return plot

def makeDefaultPlot2D ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2D()
    plot.histosStack = generateHistoList ( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDScale) + generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistosTTBar, variableName, File_TTBar, TTBarScale) + generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection) 
    plot.histodata      = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_eejj_noMEJCuts_scaled_finalOnly/"
    plot.eps_folder     = "eps_eejj_noMEJCuts_scaled_finalOnly/"
    plot.suffix         = "eejj_2D_noMEJCuts_finalOnly"


    return plot


def makeDefaultPlot2D_NSigma ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2DNSigma()
    plot.histosStack = generateHistoList ( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDScale) + generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistosTTBar, variableName, File_TTBar, TTBarScale) + generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection) 
    plot.histodata      = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_eejj_noMEJCuts_scaled_finalOnly/"
    plot.eps_folder     = "eps_eejj_noMEJCuts_scaled_finalOnly/"
    plot.suffix         = "eejj_2DNSigma_noMEJCuts_finalOnly"


    return plot


def makeDefaultPlot2D_Ratio ( variableName, 
                        histoBaseName, 
                        samplesForStackHistos, 
                        sampleForDataHisto ):

    plot = Plot2DRatio()
    plot.histosStack = generateHistoList ( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDScale) + generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistosTTBar, variableName, File_TTBar, TTBarScale) + generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection) 
    plot.histodata      = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.gif_folder     = "gif_eejj_noMEJCuts_scaled_finalOnly/"
    plot.eps_folder     = "eps_eejj_noMEJCuts_scaled_finalOnly/"
    plot.suffix         = "eejj_2DRatio_noMEJCuts_finalOnly"
    
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


for mass_point in mass_points:

    plots.append ( makeDefaultPlot ( "Mej_selected_avg_LQ" + str ( mass_point) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = mej_rebin
    plots[-1].ymin = ymin_final_selection
    plots[-1].ymax = Mej_selected_avg_max
    plots[-1].xtit = "M_{ej}^{avg} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
    plots[-1].ylog  = "yes"
    plots[-1].xmin = mej_xmin
    plots[-1].xmax = mej_xmax

    plots.append ( makeDefaultPlot ( "Mej_selected_avg_LQ" + str ( mass_point) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = mej_rebin
    plots[-1].ymin = ymin_final_selection
    # plots[-1].ymax = 1000
    plots[-1].xtit = "M_{ej}^{avg} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
    plots[-1].ylog  = "no"
    plots[-1].xmin = mej_xmin
    plots[-1].xmax = mej_xmax

    plots.append ( makeDefaultPlot ( "Mej_selected_min_LQ" + str( mass_point ) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = mej_rebin
    plots[-1].ymin = ymin_final_selection
    plots[-1].ymax = Mej_selected_min_max
    plots[-1].xtit = "M_{ej}^{min} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
    plots[-1].ylog  = "yes"
    plots[-1].xmin = mej_xmin
    plots[-1].xmax = mej_xmax

    plots.append ( makeDefaultPlot ( "Mej_selected_max_LQ" + str( mass_point ) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = mej_rebin
    plots[-1].ymin = ymin_final_selection
    plots[-1].ymax = Mej_selected_min_max
    plots[-1].xtit = "M_{ej}^{max} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
    plots[-1].ylog  = "yes"
    plots[-1].xmin = mej_xmin
    plots[-1].xmax = mej_xmax

    plots.append ( makeDefaultPlot ( "Mej_minmax_LQ" + str( mass_point ), histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = mej_rebin
    plots[-1].ymin = ymin_final_selection
    plots[-1].ymax = Mej_minmax_max
    plots[-1].xtit = "M_{ej} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
    plots[-1].ytit = "2 #times Entries"
    plots[-1].ylog  = "yes"

    plots.append ( makeDefaultPlot ( "sT_eejj_LQ" + str( mass_point ), histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].ymin = ymin_final_selection
    plots[-1].ymax = sT_eejj_max
    plots[-1].xtit = "S_{T}^{ee} (GeV), (LQ M = " + str ( mass_point ) + " selection)"
    plots[-1].ylog = "yes"

    plots.append ( makeDefaultPlot ( "Mee_LQ" + str ( mass_point ) , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
    plots[-1].rebin = mee_rebin
    plots[-1].ymin = ymin_final_selection
    plots[-1].ymax = Mee_max
    plots[-1].xtit = "M(ee) (GeV), (LQ M = " + str ( mass_point ) + " selection)"
    plots[-1].ylog  = "yes"
    
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
    
#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()

fileps = "allPlots_eejj_noMEJCuts_scaled_analysis_finalOnly.pdf"

c.Print(fileps + "[")
for i_plot, plot in enumerate(plots):
    plot.Draw(fileps, i_plot + 1 )
c.Print(fileps+"]")
