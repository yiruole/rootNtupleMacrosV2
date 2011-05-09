#!/usr/bin/env python

import pprint # for pretty printing
import math


## 2011 pre-approval ##
# data files (V00-01-04 to V00-01-06 MC ntuples, enujj skim, jetid flags + QCD fake rate for Njet>=2 and full dataset)
f1 = open("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/36.0pb-1_presel_MET45_presel_sT250_Wrescale1.18_Feb112011/analysisClass_enujjSample_tables.dat")
f2 = open("/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/35.8pb-1_QCD_UseHLTPrescales_presel_MET45_presel_sT250_Feb112011/analysisClass_enujjSample_QCD_tables.dat")
# data files (V00-01-04 to V00-01-06 MC ntuples, enujj skim, jetid flags)
#f1 = open("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_enujjskim_MET45_Jan11Prod/output_cutTable_enujjSample/analysisClass_enujjSample_tables.dat")
#f2 = open("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_enujjskim_MET45/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_tables.dat")
# data files (V00-01-04 to V00-01-06 MC ntuples)
#f1 = open("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_fullntuples_MET45_Jan11Prod/output_cutTable_enujjSample/analysisClass_enujjSample_tables.dat")
#f2 = open("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_fullntuples_MET45/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_tables.dat")
# data files (V00-01-04 to V00-01-06 MC ntuples + type1 PFMET)
#f1 = open("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_fullntuples_MET45_Jan11Prod_type1PFMET/output_cutTable_enujjSample/analysisClass_enujjSample_tables.dat")
#f2 = open("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_fullntuples_MET45_type1PFMET/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_tables.dat")
# data files (V00-00-XX MC ntuples)
#f1 = open("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/36.0pb-1_sT_presel_250_Zrescale1.20_Wrescale1.19_fullntuples_MET45/output_cutTable_enujjSample/analysisClass_enujjSample_tables.dat")
#f2 = open("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/35.8pb-1_QCD_sT_presel_250_UseHLTPrescales_fullntuples_MET45/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_tables.dat")


## Dec 2010 pre-approval ##
#f1 = open("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/33.2pb-1_sT_presel_250_Zrescale1.20_Wrescale1.06_extraPlotsDec9/output_cutTable_enujjSample/analysisClass_enujjSample_tables.dat")
#f2 = open("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/6.1pb-1_QCD_HLT30_sT_presel_250_extraPlotsDec9/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_tables.dat")
#f1 = open("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/33.2pb-1_sT_presel_250_Zrescale1.20_Wrescale1.06/output_cutTable_enujjSample/analysisClass_enujjSample_tables.dat")
#f2 = open("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/6.1pb-1_QCD_HLT30_sT_presel_250/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_tables.dat")


#dict containing all values
d = {}

###########
# User data
###########

## List of cuts
cutNames = [ "sT_MLQ200",
             "sT_MLQ250",
             "sT_MLQ280",
             "sT_MLQ300",
             "sT_MLQ320",
             "sT_MLQ340",
             "sT_MLQ370",
             "sT_MLQ400",
             "sT_MLQ450",
             "sT_MLQ500",
           ]
cutLabels = [ r"200 (\st$>350$)",
              r"250 (\st$>410$)",
              r"280 (\st$>460$)",
              r"300 (\st$>490$)",
              r"320 (\st$>520$)",
              r"340 (\st$>540$)",
              r"370 (\st$>570$)",
              r"400 (\st$>600$)",
              r"450 (\st$>640$)",
              r"500 (\st$>670$)",
            ]
LQsampleForEachCut = [ "LQenujj_M200",
                       "LQenujj_M250",
                       "LQenujj_M280",
                       "LQenujj_M300",
                       "LQenujj_M320",
                       "LQenujj_M340",
                       "LQenujj_M370",
                       "LQenujj_M400",
                       "LQenujj_M450",
                       "LQenujj_M500",
                       ]

systUncert = { 'data': {'TTbar_Madgraph': 14, 'WJetAlpgen': 10, 'QCD': 25, 'OTHERBKG': 0, 'LQ': 0, 'magnitude': '-'} ,
               'Wshape': {'TTbar_Madgraph': 0, 'WJetAlpgen': 44, 'QCD': 0, 'OTHERBKG': 0, 'LQ': 0, 'magnitude': '44'} ,
               'JES':  {'TTbar_Madgraph': 9, 'WJetAlpgen': 6, 'QCD': 0, 'OTHERBKG': 9, 'LQ': 5, 'magnitude': '5'} ,
               'EES':  {'TTbar_Madgraph': 4, 'WJetAlpgen': 2, 'QCD': 0, 'OTHERBKG': 1, 'LQ': 1, 'magnitude': '1(3)'} ,
               'lumi': {'TTbar_Madgraph': 0, 'WJetAlpgen': 0, 'QCD': 0, 'OTHERBKG': 4, 'LQ': 0, 'magnitude': '4'} ,
               'ele':  {'TTbar_Madgraph': 0, 'WJetAlpgen': 0, 'QCD': 0, 'OTHERBKG': 6, 'LQ': 6, 'magnitude': '6'} ,
               }
CutForSystematics = 'sT_MLQ300'
SampleForSystematics = 'LQenujj_M300'

## List of samples
blocks = { 'all' : {"ALLBKG":         {"rescale": 0.001, "label":  "All Bkgs"},
                    "TTbar_Madgraph": {"rescale": 0.001, "label": r"$t\bar{t}$ + jets"},
                    "WJetAlpgen":     {"rescale": 0.001, "label": r"$W$ + jets"},
                    "OTHERBKG":       {"rescale": 0.001, "label":  "Other Bkgs"},
                    "DATA":           {"rescale": 1,     "label":  "DATA"},
                    "LQenujj_M200":   {"rescale": 0.001, "label":  "LQenujj200"},
                    "LQenujj_M250":   {"rescale": 0.001, "label":  "LQenujj250"},
                    "LQenujj_M280":   {"rescale": 0.001, "label":  "LQenujj280"},
                    "LQenujj_M300":   {"rescale": 0.001, "label":  "LQenujj300"},
                    "LQenujj_M320":   {"rescale": 0.001, "label":  "LQenujj320"},
                    "LQenujj_M340":   {"rescale": 0.001, "label":  "LQenujj340"},
                    "LQenujj_M370":   {"rescale": 0.001, "label":  "LQenujj370"},
                    "LQenujj_M400":   {"rescale": 0.001, "label":  "LQenujj400"},
                    "LQenujj_M450":   {"rescale": 0.001, "label":  "LQenujj450"},
                    "LQenujj_M500":   {"rescale": 0.001, "label":  "LQenujj500"},
                    },
           'QCD' : {"DATA":           {"rescale": 1,     "label":  "QCD"},
                    }
#           'QCD' : {"DATA":           {"rescale": 5.4, "label": "QCD"},
#                    }

         }


################################################
### HOW TO SORT A DICT ACCORDINGLY WITH THE ORIGINAL ORDERING ??
#cutNames = { "Pt1stEle_PAS": r"$e\nu jj$ pre-sel.",
#             "nMuon_PtCut_IDISO": r"\nmuon$=0$",
#             "minMETPt1stEle": r"\minptmet$>85$ \GeVmom",
#             "MTenu": r"\mt$>125$ \GeVmass",
#             "sT_MLQ280": r"\st$>460$ \GeVmom",
#           }
#cutNames = { "Pt1stEle_PAS":        {"label": r"$e\nu jj$ pre-sel.", "level": 1},
#             "nMuon_PtCut_IDISO":   {"label": r"\nmuon$=0$"        , "level": 2},
#             "minMETPt1stEle":      {"label": r"\minptmet$>85$ \GeVmom", "level": 3},
#             "MTenu":               {"label": r"\mt$>125$ \GeVmass", "level": 4},
#             "sT_MLQ280":           {"label": r"\st$>460$ \GeVmom", "level": 5},
#           }
################################################



###########
# Algorithm
###########
for sample in blocks:
  if sample == "all": ff = f1
  if sample == "QCD": ff = f2
  yesBlock = False
  for ln, line in enumerate(ff): #.readlines()):
        #print "line no.", ln
        line = line.strip()
        if line in blocks[sample]:
            yesBlock = True
            block = line
            continue
        if yesBlock:
            if line == "":
              yesBlock = False
              continue
            col = line.split()
            cutName = col[0].strip()
            #print cutName
            if cutName in cutNames:

                cut = float(col[1].strip())
                Npass = float(col[5].strip())
                errNpass = float(col[6].strip())
                EffAbs = float(col[9].strip())
                errEffAbs = float(col[10].strip())

                if not cutName in d:
                    d[cutName] = {}
                if not sample in d[cutName]:
                    d[cutName][sample] = {}
                if not block in d[cutName][sample]:
                    d[cutName][sample][block] = {}
                d[cutName][sample][block]['cut'] = cut
                d[cutName][sample][block]['Npass'] = Npass * blocks[sample][block]['rescale']
                d[cutName][sample][block]['errNpass'] = errNpass * blocks[sample][block]['rescale']
                d[cutName][sample][block]['EffAbs'] = EffAbs
                d[cutName][sample][block]['errEffAbs'] = errEffAbs



######################
# Operations on values
######################

for cut in cutNames:
  if not 'ALLBKG+QCD' in d[cut]['all']:
    d[cut]['all']['ALLBKG+QCD'] = {}
  d[cut]['all']['ALLBKG+QCD']['Npass']=d[cut]['all']['ALLBKG']['Npass'] + d[cut]['QCD']['DATA']['Npass']
  d[cut]['all']['ALLBKG+QCD']['errNpass']=math.sqrt( d[cut]['all']['ALLBKG']['errNpass']**2
                                                     + (d[cut]['QCD']['DATA']['errNpass'])**2 )



######################################
# Print the dictionary with all values
######################################
pprint.pprint(d)


############################
# Output on LaTeX table
############################
fout = open("table_FinalSel_enujj.tex", "w")
#fout.write(r"\begin{table}[]"+"\n")
for idx, cutName in enumerate(cutNames):
  print cutName
  fout.write(r" %s & %.3f$\pm$%.3f & %.3f & %.2f$\pm$%.2f & %.2f$\pm$%.2f & %.2f$\pm$%.2f & %.2f$\pm$%.2f & %.2f$\pm$%.2f & %i \\" %
             (
               #cutNames[cutName]['label'],
               cutLabels[idx],
               d[cutName]['all'][LQsampleForEachCut[idx]]['Npass'],
               d[cutName]['all'][LQsampleForEachCut[idx]]['errNpass'],
               d[cutName]['all'][LQsampleForEachCut[idx]]['EffAbs'],
               #d[cutName]['all'][LQsampleForEachCut[idx]]['errEffAbs'],
               d[cutName]['all']["TTbar_Madgraph"]['Npass'],
               d[cutName]['all']["TTbar_Madgraph"]['errNpass'],
               d[cutName]['all']["WJetAlpgen"]['Npass'],
               d[cutName]['all']["WJetAlpgen"]['errNpass'],
               d[cutName]['all']["OTHERBKG"]['Npass'],
               d[cutName]['all']["OTHERBKG"]['errNpass'],
               d[cutName]['QCD']["DATA"]['Npass'],
               d[cutName]['QCD']["DATA"]['errNpass'],
               d[cutName]['all']["ALLBKG+QCD"]['Npass'],
               d[cutName]['all']["ALLBKG+QCD"]['errNpass'],
               d[cutName]['all']["DATA"]['Npass'],
               ) + "\n"
             )

fout.close()


#############
# systematics
#############


errDataALL = math.sqrt(
                        (
                          (systUncert['data']['WJetAlpgen']     * d[CutForSystematics]['all']["WJetAlpgen"]['Npass'] )
                          + (systUncert['data']['TTbar_Madgraph'] * d[CutForSystematics]['all']["TTbar_Madgraph"]['Npass'] )
                        )**2
                          + (systUncert['data']['QCD']            * d[CutForSystematics]['QCD']["DATA"]['Npass'] )**2
                          + (systUncert['data']['OTHERBKG']       * d[CutForSystematics]['all']["OTHERBKG"]['Npass'] )**2
                          + ( d[CutForSystematics]['QCD']["DATA"]['errNpass'] )**2
                      )        /       d[CutForSystematics]['all']["ALLBKG+QCD"]['Npass']

errWshapeALL = (
  (systUncert['Wshape']['WJetAlpgen']     * d[CutForSystematics]['all']["WJetAlpgen"]['Npass'] )
  )         /       d[CutForSystematics]['all']["ALLBKG+QCD"]['Npass']

errJESALL = (
  (systUncert['JES']['WJetAlpgen']     * d[CutForSystematics]['all']["WJetAlpgen"]['Npass'] )
  + (systUncert['JES']['TTbar_Madgraph'] * d[CutForSystematics]['all']["TTbar_Madgraph"]['Npass'] )
  + (systUncert['JES']['QCD']            * d[CutForSystematics]['QCD']["DATA"]['Npass'] )
  + (systUncert['JES']['OTHERBKG']       * d[CutForSystematics]['all']["OTHERBKG"]['Npass'] )
  )         /       d[CutForSystematics]['all']["ALLBKG+QCD"]['Npass']

errEESALL = (
  (systUncert['EES']['WJetAlpgen']     * d[CutForSystematics]['all']["WJetAlpgen"]['Npass'] )
  + (systUncert['EES']['TTbar_Madgraph'] * d[CutForSystematics]['all']["TTbar_Madgraph"]['Npass'] )
  + (systUncert['EES']['QCD']            * d[CutForSystematics]['QCD']["DATA"]['Npass'] )
  + (systUncert['EES']['OTHERBKG']       * d[CutForSystematics]['all']["OTHERBKG"]['Npass'] )
  )         /       d[CutForSystematics]['all']["ALLBKG+QCD"]['Npass']

errlumiALL = (
  (systUncert['lumi']['WJetAlpgen']     * d[CutForSystematics]['all']["WJetAlpgen"]['Npass'] )
  + (systUncert['lumi']['TTbar_Madgraph'] * d[CutForSystematics]['all']["TTbar_Madgraph"]['Npass'] )
  + (systUncert['lumi']['QCD']            * d[CutForSystematics]['QCD']["DATA"]['Npass'] )
  + (systUncert['lumi']['OTHERBKG']       * d[CutForSystematics]['all']["OTHERBKG"]['Npass'] )
  )         /       d[CutForSystematics]['all']["ALLBKG+QCD"]['Npass']

erreleALL = (
  (systUncert['ele']['WJetAlpgen']     * d[CutForSystematics]['all']["WJetAlpgen"]['Npass'] )
  + (systUncert['ele']['TTbar_Madgraph'] * d[CutForSystematics]['all']["TTbar_Madgraph"]['Npass'] )
  + (systUncert['ele']['QCD']            * d[CutForSystematics]['QCD']["DATA"]['Npass'] )
  + (systUncert['ele']['OTHERBKG']       * d[CutForSystematics]['all']["OTHERBKG"]['Npass'] )
  )         /       d[CutForSystematics]['all']["ALLBKG+QCD"]['Npass']




errALL = math.sqrt(errDataALL**2 + errWshapeALL**2 + errJESALL**2 + errEESALL**2 + errlumiALL**2 + erreleALL**2
                   + (100 * d[CutForSystematics]['all']["ALLBKG+QCD"]['errNpass'] / d[CutForSystematics]['all']["ALLBKG+QCD"]['Npass'])**2 )


errLQ = math.sqrt ( systUncert['data']['LQ']**2
                    + systUncert['JES']['LQ']**2
                    + systUncert['EES']['LQ']**2
                    + systUncert['lumi']['LQ']**2
                    + (100 * d[CutForSystematics]['all'][SampleForSystematics]['errNpass'] / d[CutForSystematics]['all'][SampleForSystematics]['Npass'])**2
                    + systUncert['ele']['LQ']**2
                    )

errW = math.sqrt ( systUncert['data']['WJetAlpgen']**2
                   + systUncert['Wshape']['WJetAlpgen']**2
                   + systUncert['JES']['WJetAlpgen']**2
                   + systUncert['EES']['WJetAlpgen']**2
                   + systUncert['lumi']['WJetAlpgen']**2
                   + (100 * d[CutForSystematics]['all']["WJetAlpgen"]['errNpass'] / d[CutForSystematics]['all']["WJetAlpgen"]['Npass'])**2
                   + systUncert['ele']['WJetAlpgen']**2
                   )

errTTBAR = math.sqrt ( systUncert['data']['TTbar_Madgraph']**2
                       + systUncert['JES']['TTbar_Madgraph']**2
                       + systUncert['EES']['TTbar_Madgraph']**2
                       + systUncert['lumi']['TTbar_Madgraph']**2
                       + (100 * d[CutForSystematics]['all']["TTbar_Madgraph"]['errNpass'] / d[CutForSystematics]['all']["TTbar_Madgraph"]['Npass'])**2
                       + systUncert['ele']['TTbar_Madgraph']**2
                       )

errQCD = math.sqrt ( systUncert['data']['QCD']**2
                     + systUncert['JES']['QCD']**2
                     + systUncert['EES']['QCD']**2
                     + systUncert['lumi']['QCD']**2
                     + (100 * d[CutForSystematics]['QCD']["DATA"]['errNpass'] / d[CutForSystematics]['QCD']["DATA"]['Npass'])**2
                     + systUncert['ele']['QCD']**2
                     )

errOTHER = math.sqrt ( systUncert['data']['OTHERBKG']**2
                       + systUncert['JES']['OTHERBKG']**2
                       + systUncert['EES']['OTHERBKG']**2
                       + systUncert['lumi']['OTHERBKG']**2
                       + (100 * d[CutForSystematics]['all']["OTHERBKG"]['errNpass'] / d[CutForSystematics]['all']["OTHERBKG"]['Npass'])**2
                       + systUncert['ele']["OTHERBKG"]**2
                       )

fout2 = open("table_systematics_enujj.tex", "w")

# Data-Driven Uncert.
fout2.write(r"   %s   &   %s   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   \\" %
              (
              'Data-Driven Uncert.' ,
              systUncert['data']['magnitude'] ,
              systUncert['data']['LQ'] ,
              systUncert['data']['WJetAlpgen'] ,
              systUncert['data']['TTbar_Madgraph'] ,
              math.sqrt( (systUncert['data']['QCD'])**2 + (100 * d[CutForSystematics]['QCD']["DATA"]['errNpass'] / d[CutForSystematics]['QCD']["DATA"]['Npass'])**2 ),
              systUncert['data']['OTHERBKG'] ,
              errDataALL,
              ) + "\n"
            )

#W+jets Shape
fout2.write(r"   %s   &   %s   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   \\" %
              (
              'W+jets Background Shape' ,
              systUncert['Wshape']['magnitude'] ,
              systUncert['Wshape']['LQ'] ,
              systUncert['Wshape']['WJetAlpgen'] ,
              systUncert['Wshape']['TTbar_Madgraph'] ,
              systUncert['Wshape']['QCD'] ,
              systUncert['Wshape']['OTHERBKG'] ,
              errWshapeALL,
              ) + "\n"
            )

#JES
fout2.write(r"   %s   &   %s   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   \\" %
              (
              'Jet/MET Energy Scale' ,
              systUncert['JES']['magnitude'] ,
              systUncert['JES']['LQ'] ,
              systUncert['JES']['WJetAlpgen'] ,
              systUncert['JES']['TTbar_Madgraph'] ,
              systUncert['JES']['QCD'] ,
              systUncert['JES']['OTHERBKG'] ,
              errJESALL,
              ) + "\n"
            )


#EES
fout2.write(r"   %s   &   %s   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   \\" %
              (
              'Electron Energy Scale EB(EE)' ,
              systUncert['EES']['magnitude'] ,
              systUncert['EES']['LQ'] ,
              systUncert['EES']['WJetAlpgen'] ,
              systUncert['EES']['TTbar_Madgraph'] ,
              systUncert['EES']['QCD'] ,
              systUncert['EES']['OTHERBKG'] ,
              errEESALL,
              ) + "\n"
            )

#Integrated luminosity
fout2.write(r"   %s   &   %s   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   \\" %
              (
              'Integrated Luminosity' ,
              systUncert['lumi']['magnitude'] ,
              systUncert['lumi']['LQ'] ,
              systUncert['lumi']['WJetAlpgen'] ,
              systUncert['lumi']['TTbar_Madgraph'] ,
              systUncert['lumi']['QCD'] ,
              systUncert['lumi']['OTHERBKG'] ,
              errlumiALL,
              ) + "\n"
            )

#MC Statistics
fout2.write(r"   %s   &   %s   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   \\" %
              (
              'MC Statistics' ,
              '-' ,
              100 * d[CutForSystematics]['all'][SampleForSystematics]['errNpass'] / d[CutForSystematics]['all'][SampleForSystematics]['Npass'] ,
              100 * d[CutForSystematics]['all']["WJetAlpgen"]['errNpass'] / d[CutForSystematics]['all']["WJetAlpgen"]['Npass'] ,
              100 * d[CutForSystematics]['all']["TTbar_Madgraph"]['errNpass'] / d[CutForSystematics]['all']["TTbar_Madgraph"]['Npass'] ,
              0,
              #              100 * d[CutForSystematics]['QCD']["DATA"]['errNpass'] / d[CutForSystematics]['QCD']["DATA"]['Npass'] ,
              100 * d[CutForSystematics]['all']["OTHERBKG"]['errNpass'] / d[CutForSystematics]['all']["OTHERBKG"]['Npass'] ,
              100 * d[CutForSystematics]['all']["ALLBKG+QCD"]['errNpass'] / d[CutForSystematics]['all']["ALLBKG+QCD"]['Npass'],
              ) + "\n"
            )

#Electron HLT/Reco/ID/Iso
fout2.write(r"   %s   &   %s   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   \\" %
              (
              'Electron HLT/Reco/ID/Iso' ,
              systUncert['ele']['magnitude'] ,
              systUncert['ele']['LQ'] ,
              systUncert['ele']['WJetAlpgen'] ,
              systUncert['ele']['TTbar_Madgraph'] ,
              systUncert['ele']['QCD'] ,
              systUncert['ele']['OTHERBKG'] ,
              erreleALL,
              ) + "\n"
            )

#Total
fout2.write(r"\hline\hline" + "\n" )
fout2.write(r"   %s   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   &   %.1f   \\" %
              (
              r'\multicolumn{2}{|l||}{Total}' ,
              errLQ ,
              errW ,
              errTTBAR ,
              errQCD ,
              errOTHER ,
              errALL ,
              ) + "\n"
            )


fout2.close()


# examples
# ##uncorrelated
#               math.sqrt(
#                            (systUncert['data']['WJetAlpgen']     * d[CutForSystematics]['all']["WJetAlpgen"]['Npass'] ) **2
#                          + (systUncert['data']['TTbar_Madgraph'] * d[CutForSystematics]['all']["TTbar_Madgraph"]['Npass'] ) **2
#                          + (systUncert['data']['QCD']            * d[CutForSystematics]['all']["QCD"]['Npass'] ) **2
#                          + (systUncert['data']['OTHERBKG']       * d[CutForSystematics]['all']["OTHERBKG"]['Npass'] ) **2
#                          )         /       d[CutForSystematics]['all']["ALLBKG+QCD"]['Npass'],
# #fully correlated
#               systUncert['data']['OTHERBKG'] ,
#                          (
#                            (systUncert['data']['WJetAlpgen']     * d[CutForSystematics]['all']["WJetAlpgen"]['Npass'] )
#                          + (systUncert['data']['TTbar_Madgraph'] * d[CutForSystematics]['all']["TTbar_Madgraph"]['Npass'] )
#                          + (systUncert['data']['QCD']            * d[CutForSystematics]['all']["QCD"]['Npass'] )
#                          + (systUncert['data']['OTHERBKG']       * d[CutForSystematics]['all']["OTHERBKG"]['Npass'] )
#                          )         /       d[CutForSystematics]['all']["ALLBKG+QCD"]['Npass'],

