#!/usr/bin/env python

from plot_class import *
from ROOT import *

File_EEJ          = GetFile(os.environ["LQDATA"] + "/eej_analysis/eej/output_cutTable_lq_eej/analysisClass_lq_eej_plots.root" ) 
File_preselection = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj//scaled_output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root" )
File_QCD          = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj_qcd//output_cutTable_lq_eejj/analysisClass_lq_eejj_QCD_plots.root")

zUncBand="no"
makeRatio=1
makeNSigma=1

histoBaseName = "histo1D__SAMPLE__cutHisto_allOtherCuts___________VARIABLE"
histoBaseName_userDef = "histo1D__SAMPLE__VARIABLE"

samplesForStackHistos_other = [ "OTHERBKG" ]
samplesForStackHistosQCD     = ["DATA"]
samplesForStackHistosEEJ     = ["DATA"]
samplesForStackHistos_ZJets  = [ "TTbar_Madgraph", "ZJet_Madgraph" ]
samplesForStackHistos = samplesForStackHistos_other + samplesForStackHistos_ZJets

keysStack             = [ "QCD multijets", "Other backgrounds", "t#bar{t}"  ,  "Z/#gamma* + jets"  ]
stackColorIndexes     = [ 432            , 9                  , 600         ,  kRed           ]
stackFillStyleIds     = [ 3005           , 3008               , 3004        ,  3345           ]

eejKeysStack             = [ "Data (eej)" ]
eejStackColorIndexes     = [ 1 ] 
eejStackFillStyleIds     = [ 3005 ] 


samplesForHistos = ["LQ_M900"      ]
keys             = ["LQ, M=900 GeV"]


stackColorIndexes.reverse()
stackFillStyleIds.reverse()

samplesForHistos_blank = []
keys_blank             = []

sampleForDataHisto = "DATA"

QCDScale   = 1.0
EEJScale   = 2.945
# EEJScale   = 1.0

def makeDefaultPlot_eejj_vs_MC ( variableName, histoBaseName, 
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
    plot.gif_folder        = "gif_eej_comparison/"
    plot.eps_folder        = "eps_eej_comparison/"
    plot.pdf_folder        = "pdf_eej_comparison/"
    plot.png_folder        = "png_eej_comparison/"
    plot.suffix            = "eejjVSmc"
    plot.lumi_fb           = "19.5"
    
    return plot


def makeDefaultPlot_eejj_vs_eej ( variableName, histoBaseName, 
                                 samplesForHistos, keys,
                                 samplesForStackHistos, keysStack,
                                 sampleForDataHisto,
                                 zUncBand, makeRatio ) :
    plot                   = Plot() 
    plot.histosStack       = generateHistoList( histoBaseName, samplesForStackHistosEEJ  , variableName, File_EEJ, EEJScale) 
    plot.keysStack         = eejKeysStack
    plot.histos            = generateHistoList( histoBaseName, samplesForHistos_blank, variableName, File_preselection)
    plot.keys              = keys_blank
    plot.addZUncBand       = zUncBand
    plot.makeRatio         = makeRatio
    plot.makeNSigma        = makeNSigma
    plot.histodata         = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.ytit              = "N(Events)"
    plot.ylog              = "no"
    plot.name              = variableName
    plot.stackColorIndexes = eejStackColorIndexes
    plot.stackFillStyleIds = eejStackFillStyleIds 
    plot.gif_folder        = "gif_eej_comparison/"
    plot.eps_folder        = "eps_eej_comparison/"
    plot.pdf_folder        = "pdf_eej_comparison/"
    plot.png_folder        = "png_eej_comparison/"
    plot.suffix            = "eejjVSeej"
    plot.lumi_fb           = "19.5"
    
    return plot

plots = []


plots.append ( makeDefaultPlot_eejj_vs_MC ( "sT_PAS"   ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection]"


plots.append ( makeDefaultPlot_eejj_vs_eej ( "sT_PAS"   ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection]"


plots.append ( makeDefaultPlot_eejj_vs_MC ( "sT_PASandMee100"   ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 100]"


plots.append ( makeDefaultPlot_eejj_vs_eej ( "sT_PASandMee100"   ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 200.
plots[-1].xmax = 2000.
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (GeV) [Preselection + M(ee) > 100]"


#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()

fileps = "allPlots_eej_scaled_analysis.pdf"

c.Print(fileps + "[")
for i_plot, plot in enumerate(plots):
    plot.Draw(fileps, i_plot + 1)
c.Print(fileps+"]")
# os.system('ps2pdf '+fileps)
# os.system('rm '+fileps)
