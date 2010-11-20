
import pprint # for pretty printing
import math

# data files
f1 = open("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/34.7pb-1_v7_EleEtaCut2.2/output_cutTable_enujjSample/analysisClass_enujjSample_tables.dat")
f2 = open("/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/7.4pb-1_v7_QCD_HLT30_EleEtaCut2.2/output_cutTable_enujjSample_QCD/analysisClass_enujjSample_QCD_tables.dat")

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
           ]
cutLabels = [ r"200 (\st$>350$)",
              r"250 (\st$>410$)",
              r"280 (\st$>460$)",
              r"300 (\st$>490$)",
              r"320 (\st$>520$)",
            ]
LQsampleForEachCut = [ "LQenujj_M200",
                       "LQenujj_M250",
                       "LQenujj_M280",
                       "LQenujj_M300",
                       "LQenujj_M320",
                       ]


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
                    },
           'QCD' : {"DATA":           {"rescale": 4.689, "label": "QCD"},
                    }
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
  fout.write(r" %s && %.2f\pm%.2f && %.3f && %.2f\pm%.2f && %.2f\pm%.2f && %.2f\pm%.2f && %.2f\pm%.2f && %.2f\pm%.2f && %i \\" % 
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


        

