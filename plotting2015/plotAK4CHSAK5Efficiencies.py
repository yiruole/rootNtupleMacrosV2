#!/usr/bin/env python

from ROOT import *
import os
import math
import array
from prettytable import PrettyTable
from collections import OrderedDict
 

class VarInfo:
  def __init__(self,initVarName='',initMin1=0.0,initMax1=0.0,initMin2=0.0,initMax2=0.0,initNpass=0.0,initErrNpass=0.0,initEffRel=0.0,initErrEffRel=0.0,initEffAbs=0.0,initErrEffAbs=0.0):
    self.varName = initVarName
    self.min1 = initMin1
    self.max1 = initMax1
    self.min2 = initMin2
    self.max2 = initMax2
    self.nPass = initNpass
    self.errNpass = initErrNpass
    self.effRel = initEffRel
    self.errEffRel = initErrEffRel
    self.effAbs = initEffAbs
    self.errEffAbs = initErrEffAbs

class SampleData:
  def __init__(self,initSampleName=''):
    self.sampleName = initSampleName
    self.varInfoDict = OrderedDict()

def GetVal(item):
  if item != '-':
    return float(item)
  else:
    return -1.0

def FillSampleDataDict(tableFile,sampleDataDict):
  with open(tableFile) as tableFileOpened:
    for line in tableFileOpened:
      #if len(line) <= 0:
      #  continue
      split = line.split()
      if len(split) == 1:
        sampleName = split[0]
        currentSample = sampleName
        sampleData = SampleData(sampleName)
        sampleDataDict[sampleName] = sampleData
      elif len(split) > 1:
        if 'variableName' in line:
          continue
        try:
          varName = split[0]
        except IndexError:
          print 'Could not get varName from split[0]; line looks like:',line
          print '    split looks like:',split
          print 'Quitting'
          exit(-1)
        min1 = GetVal(split[1])
        max1 = GetVal(split[2])
        min2 = GetVal(split[3])
        max2 = GetVal(split[4])
        nPass = GetVal(split[5])
        errNpass = GetVal(split[6])
        effRel = GetVal(split[7])
        errEffRel = GetVal(split[8])
        effAbs = GetVal(split[9])
        errEffAbs = GetVal(split[10])
        sampleDataDict[currentSample].varInfoDict[varName] = VarInfo(varName,min1,max1,min2,max2,nPass,errNpass,effRel,errEffRel,effAbs,errEffAbs)


####################################################################################################
# user config options
####################################################################################################
fileAK4CHS = os.getenv('LQDATA')+'/RunII/'+'eejj_analysis_7Dec2015_AK4CHS_2094invPb_MiniAODV2_presels_hltDataOnly_withOldMETFilters//output_cutTable_lq_eejj_preselectionOnly/analysisClass_lq_eejj_preselectionOnly_tables.dat'
fileAK5 = os.getenv('LQDATA')+'/RunII/'+'eejj_analysis_7Dec2015_AK5_2094invPb_MiniAODV2_presels_hltDataOnly_withOldMETFilters/output_cutTable_lq_eejj_preselectionOnly/analysisClass_lq_eejj_preselectionOnly_tables.dat'
doPrint = True
####################################################################################################
# end of user config options
####################################################################################################


####################################################################################################
# Run
####################################################################################################
# fill sample data objects
sampleDataDictAK4CHS = {}
FillSampleDataDict(fileAK4CHS,sampleDataDictAK4CHS)
sampleDataDictAK5 = {}
FillSampleDataDict(fileAK5,sampleDataDictAK5)

samplePreselectionAK4FinalAbsEffDict = {}
samplePreselectionAK4FinalAbsEffErrDict = {}
samplePreselectionAK5FinalAbsEffDict = {}
samplePreselectionAK5FinalAbsEffErrDict = {}

# loop over them and make tables, plots, etc.
for keyAK4,sampleDataAK4 in sampleDataDictAK4CHS.iteritems():
  sampleDataAK5 = sampleDataDictAK5[keyAK4]
  t = PrettyTable(['VarName', 'EffRelAK5','EffRelAK4','deltaEffRel/effRel','EffAbsAK5','EffAbsAK4','deltaEffAvs/effAbs'])
  t.float_format = "4.3"
  #t.align['VarName'] = 'l'
  t.align = 'l'
  maxAbsRelEffRelDiff = 0.0
  maxAbsRelEffRelVar = ''
  if len(sampleDataAK4.varInfoDict) < 1:
    continue
  for var,varInfoAK4 in sampleDataAK4.varInfoDict.iteritems():
    #print 'variable:',var,'effRelAK4',varInfoAK4.effRel,'effRelAK5:',sampleDataAK5.varInfoDict[var].effRel
    try:
      relEffRelDiff = (varInfoAK4.effRel-sampleDataAK5.varInfoDict[var].effRel)/sampleDataAK5.varInfoDict[var].effRel
      absEffRelDiff = (varInfoAK4.effAbs-sampleDataAK5.varInfoDict[var].effAbs)/sampleDataAK5.varInfoDict[var].effAbs
    except ZeroDivisionError:
      print 'WARN: AK5 had 0 passing events for sample:',keyAK4,'and var:',var
      continue
    t.add_row([var,sampleDataAK5.varInfoDict[var].effRel,varInfoAK4.effRel,relEffRelDiff,sampleDataAK5.varInfoDict[var].effAbs,varInfoAK4.effAbs,absEffRelDiff])
    if math.fabs(relEffRelDiff) > maxAbsRelEffRelDiff:
      maxAbsRelEffRelDiff = math.fabs(relEffRelDiff)
      maxAbsRelEffRelVar = var
  samplePreselectionAK4FinalAbsEffDict[keyAK4] = sampleDataAK4.varInfoDict[next(reversed(sampleDataAK4.varInfoDict))].effAbs
  samplePreselectionAK4FinalAbsEffErrDict[keyAK4] = sampleDataAK4.varInfoDict[next(reversed(sampleDataAK4.varInfoDict))].errEffAbs
  samplePreselectionAK5FinalAbsEffDict[keyAK4] = sampleDataAK5.varInfoDict[next(reversed(sampleDataAK5.varInfoDict))].effAbs
  samplePreselectionAK5FinalAbsEffErrDict[keyAK4] = sampleDataAK5.varInfoDict[next(reversed(sampleDataAK5.varInfoDict))].errEffAbs
  if doPrint:
    print 'sample:',keyAK4,'MaxRelEffRelDiff for var:',maxAbsRelEffRelVar
    print t
    print

# make graph
hist = TH1F("hist","Final absolute efficiency",len(samplePreselectionAK4FinalAbsEffDict),0,len(samplePreselectionAK4FinalAbsEffDict)-1)
indexList = []
ak4effList = []
ak4effErrList = []
ak5effList = []
ak5effErrList = []
for index,key in enumerate(sorted(samplePreselectionAK4FinalAbsEffDict)):
  hist.GetXaxis().SetBinLabel(index+1,key)
  indexList.append(index+0.5)
  ak4effList.append(samplePreselectionAK4FinalAbsEffDict[key])
  ak4effErrList.append(samplePreselectionAK4FinalAbsEffErrDict[key])
  ak5effList.append(samplePreselectionAK5FinalAbsEffDict[key])
  ak5effErrList.append(samplePreselectionAK5FinalAbsEffErrDict[key])

c = TCanvas()
c.SetLogy()
c.SetGridx()
hist.SetStats(False)
hist.SetMaximum(1)
hist.SetMinimum(5e-9)
hist.GetXaxis().SetLabelSize(0.02)
hist.GetXaxis().SetNdivisions(hist.GetNbinsX())
hist.Draw()
#gr = TGraph(len(samplePreselectionAK4FinalAbsEffDict),array.array('f',indexList),array.array('f',ak4effList))
gr4 = TGraphErrors(len(samplePreselectionAK4FinalAbsEffDict),array.array('f',indexList),array.array('f',ak4effList),array.array('f',[0]*len(ak4effList)),array.array('f',ak4effErrList))
gr5 = TGraphErrors(len(samplePreselectionAK5FinalAbsEffDict),array.array('f',indexList),array.array('f',ak5effList),array.array('f',[0]*len(ak5effList)),array.array('f',ak5effErrList))
gr4.Draw('p')
gr5.SetLineColor(2)
gr5.SetMarkerColor(2)
gr5.Draw('p')

t = TLegend(0.2,0.2,0.3,0.3)
t.SetBorderSize(0)
t.AddEntry(gr5,'AK5','lp')
t.AddEntry(gr4,'AK4','lp')
t.Draw()

## wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
if __name__ == '__main__':
   rep = ''
   while not rep in [ 'q', 'Q' ]:
      rep = raw_input( 'enter "q" to quit: ' )
      if 1 < len(rep):
         rep = rep[0]

exit()
