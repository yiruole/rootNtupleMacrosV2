#!/usr/bin/env python

#
# To combine the root files (inside the LQDATA/blah directory):
# for i in `find . -maxdepth 1 -type d`; do cd $i; echo $PWD; rm combined.root; hadd combined.root *.root; cd $LQDATA/trigEffStudyEejjSignals; done;
#
from ROOT import *
import os
import copy


def SetGraphStyle(graph,counter):
  if counter==0:
    graph.SetLineColor(4)
    graph.SetMarkerColor(4)
    graph.SetMarkerStyle(22)
  elif counter==1:
    graph.SetLineColor(2)
    graph.SetMarkerColor(2)
    graph.SetMarkerStyle(23)
  elif counter==2:
    graph.SetLineColor(8)
    graph.SetMarkerColor(8)
  elif counter==3:
    graph.SetLineColor(9)
    graph.SetMarkerColor(9)
  elif counter==4:
    graph.SetLineColor(kMagenta-4)
    graph.SetMarkerColor(kMagenta-4)

  return graph


def MakeTrigEffGraphs(inputFilesTrigs, IsEEJJ):
  trigEffGraphs = []
  for fileName in inputFilesTrigs:
    tfile = TFile.Open(fileName)
    print 'MakeTrigEffGraphs() looking in file:',fileName
  
    denominator = tfile.Get('NEventsTwoGenEleECAL') if IsEEJJ else tfile.Get('NEventsOneGenEleECAL')
    numerator = tfile.Get('NEventsTwoGenEleECALPlusTrig') if IsEEJJ else tfile.Get('NEventsOneGenEleECALPlusTrig')
    print '\tdenominator bins:',denominator.GetNbinsX()
    print '\tnumerator bins:',numerator.GetNbinsX()
  
    trigEff = TGraphAsymmErrors()
    trigEff.BayesDivide(numerator,denominator)
    trigEff.GetXaxis().SetTitle('LQ mass [GeV]')
    trigEff.GetXaxis().SetNdivisions(512)
    if IsEEJJ:
      trigEff.GetYaxis().SetTitle('TrigEff. [eejj]')
    else:
      trigEff.GetYaxis().SetTitle('TrigEff. [e#nujj]')
    trigEff.GetYaxis().SetTitleOffset(0.7)
    trigEff.GetYaxis().SetNdivisions(514)
    trigNum = fileName[fileName.find('trigger')+7:fileName.find('trigger')+8]
    trigEff.SetName(GetTrigNameFromIndex(trigNum)+'_trigEffGraph')
    #t3 = TCanvas()
    #t3.cd()
    #trigEff.Draw('ap')
    trigEffGraphs.append(trigEff)
  
    tfile.Close()
  return trigEffGraphs


def GetTrigNameFromIndex(index):
  index=int(index)
  if index==1:
    return 'HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1'
  elif index==2:
    return 'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v1'
  elif index==3:
    return 'HLT_Ele27_WP85_Gsf_v1'
  elif index==4:
    return 'HLT_Ele27_WPTight_Gsf_v2'
  elif index==5:
    return 'HLT_Ele27_WPLoose_Gsf_v2'
  elif index==6:
    return 'HLT_Ele27_WPTight_Gsf || HLT_Photon175'
  elif index==7:
    return 'HLT_Ele27_WPTight_Gsf || HLT_Ele105_CaloIdVT_GsfTrkIdT'


def GetIDNameFromIndex(index):
  index=int(index)
  if index==0:
    return 'HEEPv6.1'
  elif index==1:
    return 'HEEPv6.0'
  elif index==2:
    return 'HEEPv5.1'
  elif index==3:
    return 'EGamma Tight'
  elif index==4:
    return 'EGamma Medium'
  elif index==5:
    return 'EGamma Loose'



def MakeIDEffGraphs(inputFilesIDs, IsEEJJ):
  trigNames = []
  idNames = []
  idEffGraphs = []
  trigAndIdEffGraphs = []
  idEffVsPt1stGenEleGraphs = []
  for i,fileName in enumerate(inputFilesIDs):
    tfile = TFile.Open(fileName)
    print 'MakeIDEffGraphs(): looking in file:',fileName
  
    #nevtsTwoGenEleEcalWithTrigHist = tfile.Get('NEventsTwoGenEleECALPlusTrig')
    #nevtsTwoGenEleEcalWithTrigWithTwoIDHist = tfile.Get('NEventsTwoGenEleECALPlusTrigPlusTwoID')
    denominator = tfile.Get('NEventsTwoGenEleECALPlusTrig') if IsEEJJ else tfile.Get('NEventsOneGenEleECALPlusTrig')
    numerator = tfile.Get('NEventsTwoGenEleECALPlusTrigPlusTwoID') if IsEEJJ else tfile.Get('NEventsOneGenEleECALPlusTrigPlusOneID')
    numeratorTrigAndID = tfile.Get('NEventsTwoGenEleECALPlusTrigPlusTwoID') if IsEEJJ else tfile.Get('NEventsOneGenEleECALPlusTrigPlusOneID')
    denominatorTrigAndID = tfile.Get('NEventsTwoGenEleECAL') if IsEEJJ else tfile.Get('NEventsOneGenEleECAL')
    
    idEff = TGraphAsymmErrors()
    idEff.BayesDivide(numerator,denominator)
    idEff.GetXaxis().SetTitle('LQ mass [GeV]')
    idEff.GetXaxis().SetNdivisions(512)
    if IsEEJJ:
      idEff.GetYaxis().SetTitle('Eff. [eejj]')
    else:
      idEff.GetYaxis().SetTitle('Eff. [e#nujj]')
    idEff.GetYaxis().SetTitleOffset(0.9)
    idEff.GetYaxis().SetNdivisions(514)
    #t4 = TCanvas()
    #t4.cd()
    #idEff.Draw('ap')
    idEffGraphs.append(idEff)
    # trig and ID
    trigAndIdEff = TGraphAsymmErrors()
    trigAndIdEff.BayesDivide(numeratorTrigAndID,denominatorTrigAndID)
    trigAndIdEff.GetXaxis().SetTitle('LQ mass [GeV]')
    trigAndIdEff.GetXaxis().SetNdivisions(512)
    if IsEEJJ:
      trigAndIdEff.GetYaxis().SetTitle('Eff. [eejj]')
    else:
      trigAndIdEff.GetYaxis().SetTitle('Eff. [e#nujj]')
    trigAndIdEff.GetYaxis().SetTitleOffset(0.9)
    trigAndIdEff.GetYaxis().SetNdivisions(514)
    trigAndIdEffGraphs.append(trigAndIdEff)
    #
    ptFirstGenElePlusTrig = tfile.Get('Pt1stGenEle_passTrigger')
    ptFirstGenElePlusTrig.Rebin(2)
    ptFirstGenElePlusTrigPlusID = tfile.Get('Pt1stGenEle_passTrigger_passID')
    ptFirstGenElePlusTrigPlusID.Rebin(2)
    idEffPt = TGraphAsymmErrors()
    idEffPt.BayesDivide(ptFirstGenElePlusTrigPlusID,ptFirstGenElePlusTrig)
    idEffPt.GetXaxis().SetTitle('1st e GenPt [GeV]')
    idEffPt.GetXaxis().SetNdivisions(512)
    idEffPt.GetYaxis().SetTitle('Eff. [eejj]')
    idEffPt.GetYaxis().SetTitleOffset(0.9)
    idEffPt.GetYaxis().SetNdivisions(514)
    idEffVsPt1stGenEleGraphs.append(idEffPt)
    trigNum = fileName[fileName.find('trigger')+7:fileName.find('trigger')+8]
    trigName = GetTrigNameFromIndex(trigNum)
    trigNames.append(trigName)
    #print 'fileName=',fileName
    #print 'trigNum=',trigNum
    #print 'GetTrigNameFromIndex(trigNum)=',GetTrigNameFromIndex(trigNum)
    idNum = fileName[fileName.find('id')+2:fileName.find('id')+3]
    idNames.append(GetIDNameFromIndex(idNum))
    #print 'idNum=',idNum
    #print 'GetIDNameFromIndex(idNum)=',GetIDNameFromIndex(idNum)
    #if i%2==0:
    #  trigNames.append('HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1')
    #else:
    #  trigNames.append('HLT_Ele27_WP85_Gsf_v1')
    #if i < 2:
    #  idNames.append('HEEPv6.0')
    #elif i < 4:
    #  idNames.append('HEEPv5.1')
    #elif i < 6:
    #  idNames.append('Egamma Tight')
    #elif i < 8:
    #  idNames.append('Egamma Medium')
    #else:
    #  idNames.append('Egamma Loose')
  
    tfile.Close()
  return trigNames, idNames, idEffGraphs, trigAndIdEffGraphs, idEffVsPt1stGenEleGraphs


def MakeTrigEffMultigraph(trigEffGraphs,IsEEJJ):
  mg = TMultiGraph()
  leg = TLegend(0.375,0.178,0.852,0.451)
  for i,graph in enumerate(trigEffGraphs):
    if i==0:
      graph.SetLineColor(2)
      graph.SetMarkerColor(2)
    elif i==1:
      graph.SetLineColor(4)
      graph.SetMarkerColor(4)
    elif i==2:
      graph.SetLineColor(8)
      graph.SetMarkerColor(8)
    elif i==3:
      graph.SetLineColor(6)
      graph.SetMarkerColor(6)
    graphName = graph.GetName()
    leg.AddEntry(graph,graphName[:graphName.find('_trigEff')-7],'lp')
    mg.Add(graph)
  tc = TCanvas()
  tc.cd()
  mg.Draw('ap')
  mg.GetXaxis().SetTitle('LQ mass [GeV]')
  mg.GetXaxis().SetNdivisions(512)
  if IsEEJJ:
    mg.GetYaxis().SetTitle('Efficiency [eejj]')
  else:
    mg.GetYaxis().SetTitle('Efficiency [e#nujj]')
  mg.GetYaxis().SetTitleOffset(0.7)
  mg.GetYaxis().SetNdivisions(514)
  mg.Draw('ap')
  leg.SetBorderSize(0)
  leg.Draw()
  imageName='triggerEff_'+('eejj' if IsEEJJ else 'enujj')+'.png'
  tc.Print(imageName)


def MakeIDEffMultigraph(trigger,trigNames,idNames,idEffGraphs,IsEEJJ):
  mg2 = TMultiGraph()
  leg2 = TLegend(0.375,0.178,0.852,0.451)
  leg2.SetHeader(trigger)
  counter = 0
  if not trigger in trigNames:
    print 'ERROR: requested trigger',trigger,'not in trigNames=',trigNames
    exit(-1)
  for i,graph in enumerate(idEffGraphs):
    if trigNames[i] != trigger:
      continue
    graph = SetGraphStyle(graph,counter)
    leg2.AddEntry(graph,idNames[i],'lp')
    mg2.Add(graph)
    counter+=1
  t6 = TCanvas()
  t6.cd()
  mg2.Draw('ap')
  mg2.GetXaxis().SetTitle('LQ mass [GeV]')
  mg2.GetXaxis().SetNdivisions(512)
  if IsEEJJ:
    mg2.GetYaxis().SetTitle('Efficiency [eejj]')
  else:
    mg2.GetYaxis().SetTitle('Efficiency [e#nujj]')
  mg2.GetYaxis().SetTitleOffset(0.8)
  mg2.GetYaxis().SetNdivisions(514)
  mg2.GetYaxis().SetRangeUser(0.6,1.0)
  if IsEEJJ:
    mg2.GetYaxis().SetRangeUser(0.3,1.0)
  mg2.Draw('ap')
  leg2.SetBorderSize(0)
  leg2.Draw()
  imageName = 'idEff_'+trigger+('_eejj' if IsEEJJ else '_enujj')+'.png'
  t6.Print(imageName)


def MakeTrigAndIDEffMultigraph(trigger,trigNames,idNames,trigAndIdEffGraphsOrig,IsEEJJ):
  # work with a copy, or else there are seg faults later
  trigAndIdEffGraphs = copy.deepcopy(trigAndIdEffGraphsOrig)
  mg3 = TMultiGraph()
  leg3 = TLegend(0.375,0.178,0.852,0.451)
  leg3.SetHeader(trigger)
  counter = 0
  if not trigger in trigNames:
    print 'ERROR: requested trigger',trigger,'not in trigNames=',trigNames
    exit(-1)
  for i,graph in enumerate(trigAndIdEffGraphs):
    if trigNames[i] != trigger:
      continue
    graph = SetGraphStyle(graph,counter)
    leg3.AddEntry(graph,idNames[i],'lp')
    mg3.Add(graph)
    counter+=1
  t7 = TCanvas()
  t7.cd()
  t7.SetGridy()
  mg3.Draw('ap')
  mg3.GetXaxis().SetTitle('LQ mass [GeV]')
  mg3.GetXaxis().SetNdivisions(512)
  if IsEEJJ:
    mg3.GetYaxis().SetTitle('Efficiency [eejj]')
  else:
    mg3.GetYaxis().SetTitle('Efficiency [e#nujj]')
  mg3.GetYaxis().SetTitleOffset(0.8)
  mg3.GetYaxis().SetNdivisions(514)
  mg3.GetYaxis().SetRangeUser(0.6,1.0)
  if IsEEJJ:
    mg3.GetYaxis().SetRangeUser(0.3,1.0)
  mg3.Draw('ap')
  leg3.SetBorderSize(0)
  leg3.Draw()
  imageName = 'trigAndIdEff_'+trigger+('_eejj' if IsEEJJ else '_enujj')+'.png'
  t7.Print(imageName)


def MakeTrigAndIDEffMultigraphWithID(idname,trigNames,idNames,trigAndIdEffGraphsOrig,IsEEJJ):
  # work with a copy, or else there are seg faults later
  trigAndIdEffGraphs = copy.deepcopy(trigAndIdEffGraphsOrig)
  mg3 = TMultiGraph()
  if IsEEJJ:
    leg3 = TLegend(0.375,0.6,0.852,0.9)
  else:
    leg3 = TLegend(0.375,0.178,0.852,0.451)
  leg3.SetHeader(idname)
  counter = 0
  if not idname in idNames:
    print 'ERROR: requested ID',idname,'not in idNames=',idNames
    exit(-1)
  for i,graph in enumerate(trigAndIdEffGraphs):
    if idNames[i] != idname:
      continue
    graph = SetGraphStyle(graph,counter)
    leg3.AddEntry(graph,trigNames[i],'lp')
    mg3.Add(graph)
    counter+=1
  t7 = TCanvas()
  t7.cd()
  t7.SetGridy()
  mg3.Draw('ap')
  mg3.GetXaxis().SetTitle('LQ mass [GeV]')
  mg3.GetXaxis().SetNdivisions(512)
  if IsEEJJ:
    mg3.GetYaxis().SetTitle('Efficiency [eejj]')
  else:
    mg3.GetYaxis().SetTitle('Efficiency [e#nujj]')
  mg3.GetYaxis().SetTitleOffset(0.8)
  mg3.GetYaxis().SetNdivisions(514)
  mg3.GetYaxis().SetRangeUser(0.5,0.9)
  if IsEEJJ:
    mg3.GetYaxis().SetRangeUser(0.5,0.9)
  mg3.Draw('ap')
  leg3.SetBorderSize(0)
  leg3.Draw()
  imageName = 'trigAndIdEff_'+idname.replace(' ','_')+('_eejj' if IsEEJJ else '_enujj')+'.png'
  t7.Print(imageName)


def MakeFirstEleGenPtMultigraph(trigger,trigNames,idNames,idEffVsPt1stGenEleGraphs,IsEEJJ):
  mg4 = TMultiGraph()
  leg4 = TLegend(0.375,0.178,0.852,0.451)
  leg4.SetHeader(trigger)
  counter = 0
  for i,graph in enumerate(idEffVsPt1stGenEleGraphs):
    if trigNames[i] != trigger:
      continue
    graph = SetGraphStyle(graph,counter)
    leg4.AddEntry(graph,idNames[i],'lp')
    mg4.Add(graph)
    counter+=1
  t8 = TCanvas()
  t8.cd()
  mg4.Draw('ap')
  mg4.GetXaxis().SetTitle('Lead e GenPt [GeV]')
  mg4.GetXaxis().SetNdivisions(512)
  if IsEEJJ:
    mg4.GetYaxis().SetTitle('Efficiency [eejj]')
  else:
    mg4.GetYaxis().SetTitle('Efficiency [e#nujj]')
  mg4.GetYaxis().SetTitleOffset(0.8)
  mg4.GetYaxis().SetNdivisions(514)
  mg4.GetYaxis().SetRangeUser(0.3,0.9)
  mg4.Draw('ap')
  leg4.SetBorderSize(0)
  leg4.Draw()
  imageName = 'idEffVs1stEleGenPt_'+trigger+('_eejj' if IsEEJJ else '_enujj')+'.png'
  t8.Print(imageName)



#####################################################
# RUN
#####################################################

gROOT.SetBatch(True)

lqdata = os.environ['LQDATA']
#eejjAnalysisDir = lqdata+'/2016analysis/trigEffStudies_2016sep13/'
eejjAnalysisDir = lqdata+'/2016analysis/trigEffStudies_2016sep19/'

# diff trigs
eejjInputFilesTrigs = []
eejjInputFilesTrigs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger4_id1_effiStudy/combined_eejj.root')
eejjInputFilesTrigs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger5_id1_effiStudy/combined_eejj.root')
eejjInputFilesTrigs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger4_id0_effiStudy/combined_eejj.root')
eejjInputFilesTrigs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger5_id0_effiStudy/combined_eejj.root')
eejjInputFilesTrigs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger6_id0_effiStudy/combined_eejj.root')
eejjInputFilesTrigs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger7_id0_effiStudy/combined_eejj.root')
#eejjInputFilesTrigs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger3_id2_effiStudy/combined.root')
#eejjInputFilesTrigs.append('/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_trigEff_9Dec2015_AK5_reHLTSignalsFromRaw/output_cutTable_lq1_eejj_trigger5_id1_effiStudy/combined.root')
# diff IDs
eejjInputFilesIDs = []
#eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger1_id1_effiStudy/combined.root')
#eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger3_id1_effiStudy/combined.root')
#eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger1_id2_effiStudy/combined.root')
#eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger3_id2_effiStudy/combined.root')
#eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger1_id3_effiStudy/combined.root')
#eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger3_id3_effiStudy/combined.root')
#eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger1_id4_effiStudy/combined.root')
#eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger3_id4_effiStudy/combined.root')
##eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger1_id5_effiStudy/combined.root')
##eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger3_id5_effiStudy/combined.root')
#eejjInputFilesIDs.append('/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_trigEff_9Dec2015_AK5_reHLTSignalsFromRaw/output_cutTable_lq1_eejj_trigger5_id1_effiStudy/combined.root')
#eejjInputFilesIDs.append('/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/eejj_trigEff_9Dec2015_AK5_reHLTSignalsFromRaw/output_cutTable_lq1_eejj_trigger5_id4_effiStudy/combined.root')
eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger4_id1_effiStudy/combined_eejj.root')
eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger5_id1_effiStudy/combined_eejj.root')
eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger4_id0_effiStudy/combined_eejj.root')
eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger5_id0_effiStudy/combined_eejj.root')
eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger6_id0_effiStudy/combined_eejj.root')
eejjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger7_id0_effiStudy/combined_eejj.root')

## ENUJJ
#enujjAnalysisDir = lqdata+'/2016analysis/trigEffStudies_2016sep13/'
enujjAnalysisDir = lqdata+'/2016analysis/trigEffStudies_2016sep19/'
# diff trigs
enujjInputFilesTrigs = []
enujjInputFilesTrigs.append(enujjAnalysisDir+'output_cutTable_lq1_eejj_trigger4_id1_effiStudy/combined_enujj.root')
enujjInputFilesTrigs.append(enujjAnalysisDir+'output_cutTable_lq1_eejj_trigger5_id1_effiStudy/combined_enujj.root')
enujjInputFilesTrigs.append(enujjAnalysisDir+'output_cutTable_lq1_eejj_trigger4_id0_effiStudy/combined_enujj.root')
enujjInputFilesTrigs.append(enujjAnalysisDir+'output_cutTable_lq1_eejj_trigger5_id0_effiStudy/combined_enujj.root')
enujjInputFilesTrigs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger6_id0_effiStudy/combined_enujj.root')
enujjInputFilesTrigs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger7_id0_effiStudy/combined_enujj.root')
#enujjInputFilesTrigs.append(enujjAnalysisDir+'output_cutTable_lq1_enujj_trigger4_id1_effiStudy/combined.root')
#enujjInputFilesTrigs.append('/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/enujj_trigEff_10Dec2015_AK5_reHLT/output_cutTable_lq1_enujj_trigger5_id1_effiStudy/combined.root')
# diff IDs
enujjInputFilesIDs = []
#enujjInputFilesIDs.append(enujjAnalysisDir+'output_cutTable_lq1_enujj_trigger1_id1_effiStudy/combined.root')
#enujjInputFilesIDs.append(enujjAnalysisDir+'output_cutTable_lq1_enujj_trigger1_id2_effiStudy/combined.root')
#enujjInputFilesIDs.append(enujjAnalysisDir+'output_cutTable_lq1_enujj_trigger1_id3_effiStudy/combined.root')
#enujjInputFilesIDs.append(enujjAnalysisDir+'output_cutTable_lq1_enujj_trigger1_id4_effiStudy/combined.root')
##enujjInputFilesIDs.append(enujjAnalysisDir+'output_cutTable_lq1_enujj_trigger1_id5_effiStudy/combined.root')
#enujjInputFilesIDs.append(enujjAnalysisDir+'output_cutTable_lq1_enujj_trigger2_id1_effiStudy/combined.root')
#enujjInputFilesIDs.append(enujjAnalysisDir+'output_cutTable_lq1_enujj_trigger3_id1_effiStudy/combined.root')
#enujjInputFilesIDs.append('/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/enujj_trigEff_10Dec2015_AK5_reHLT/output_cutTable_lq1_enujj_trigger5_id1_effiStudy/combined.root')
#enujjInputFilesIDs.append('/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/RunII/enujj_trigEff_10Dec2015_AK5_reHLT/output_cutTable_lq1_enujj_trigger5_id4_effiStudy/combined.root')
enujjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger4_id1_effiStudy/combined_enujj.root')
enujjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger5_id1_effiStudy/combined_enujj.root')
enujjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger4_id0_effiStudy/combined_enujj.root')
enujjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger5_id0_effiStudy/combined_enujj.root')
enujjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger6_id0_effiStudy/combined_enujj.root')
enujjInputFilesIDs.append(eejjAnalysisDir+'output_cutTable_lq1_eejj_trigger7_id0_effiStudy/combined_enujj.root')


eejjTrigEffGraphs = MakeTrigEffGraphs(eejjInputFilesTrigs,True)
eejjTrigNames, eejjIdNames, eejjIdEffGraphs, eejjTrigAndIdEffGraphs, eejjIdEffVsPt1stGenEleGraphs = MakeIDEffGraphs(eejjInputFilesIDs,True)
# EEJJ trig eff multigraph
MakeTrigEffMultigraph(eejjTrigEffGraphs,True)
## EEJJ ID eff multigraph: trig 1 = HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1
#MakeIDEffMultigraph('HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1',eejjTrigNames,eejjIdNames,eejjIdEffGraphs,True)
MakeIDEffMultigraph('HLT_Ele27_WPTight_Gsf_v2',eejjTrigNames,eejjIdNames,eejjIdEffGraphs,True)
MakeIDEffMultigraph('HLT_Ele27_WPLoose_Gsf_v2',eejjTrigNames,eejjIdNames,eejjIdEffGraphs,True)
#MakeFirstEleGenPtMultigraph('HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1',eejjTrigNames,eejjIdNames,eejjIdEffVsPt1stGenEleGraphs,True)
MakeFirstEleGenPtMultigraph('HLT_Ele27_WPTight_Gsf_v2',eejjTrigNames,eejjIdNames,eejjIdEffVsPt1stGenEleGraphs,True)
MakeFirstEleGenPtMultigraph('HLT_Ele27_WPLoose_Gsf_v2',eejjTrigNames,eejjIdNames,eejjIdEffVsPt1stGenEleGraphs,True)
#MakeTrigAndIDEffMultigraph('HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1',eejjTrigNames,eejjIdNames,eejjTrigAndIdEffGraphs,True)
MakeTrigAndIDEffMultigraph('HLT_Ele27_WPTight_Gsf_v2',eejjTrigNames,eejjIdNames,eejjTrigAndIdEffGraphs,True)
MakeTrigAndIDEffMultigraph('HLT_Ele27_WPLoose_Gsf_v2',eejjTrigNames,eejjIdNames,eejjTrigAndIdEffGraphs,True)
#
MakeTrigAndIDEffMultigraphWithID('HEEPv6.0',eejjTrigNames,eejjIdNames,eejjTrigAndIdEffGraphs,True)
MakeTrigAndIDEffMultigraphWithID('HEEPv6.1',eejjTrigNames,eejjIdNames,eejjTrigAndIdEffGraphs,True)
#MakeTrigAndIDEffMultigraphWithID('EGamma Medium',eejjTrigNames,eejjIdNames,eejjTrigAndIdEffGraphs,True)

enujjTrigEffGraphs = MakeTrigEffGraphs(enujjInputFilesTrigs,False)
enujjTrigNames, enujjIdNames, enujjIdEffGraphs, enujjTrigAndIdEffGraphs, enujjIdEffVsPt1stGenEleGraphs = MakeIDEffGraphs(enujjInputFilesIDs,False)
# ENUJJ trig eff multigraph
MakeTrigEffMultigraph(enujjTrigEffGraphs,False)
# ENUJJ ID eff multigraph: trig 1 = HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1
#MakeIDEffMultigraph('HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1',enujjTrigNames,enujjIdNames,enujjIdEffGraphs,False)
MakeIDEffMultigraph('HLT_Ele27_WPTight_Gsf_v2',enujjTrigNames,enujjIdNames,enujjIdEffGraphs,False)
MakeIDEffMultigraph('HLT_Ele27_WPLoose_Gsf_v2',enujjTrigNames,enujjIdNames,enujjIdEffGraphs,False)
##MakeFirstEleGenPtMultigraph('HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1',enujjTrigNames,enujjIdNames,enujjIdEffVsPt1stGenEleGraphs,False)
#MakeTrigAndIDEffMultigraph('HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1',enujjTrigNames,enujjIdNames,enujjTrigAndIdEffGraphs,False)
MakeTrigAndIDEffMultigraph('HLT_Ele27_WPTight_Gsf_v2',enujjTrigNames,enujjIdNames,enujjTrigAndIdEffGraphs,False)
MakeTrigAndIDEffMultigraph('HLT_Ele27_WPLoose_Gsf_v2',enujjTrigNames,enujjIdNames,enujjTrigAndIdEffGraphs,False)
##
MakeTrigAndIDEffMultigraphWithID('HEEPv6.0',enujjTrigNames,enujjIdNames,enujjTrigAndIdEffGraphs,False)
MakeTrigAndIDEffMultigraphWithID('HEEPv6.1',enujjTrigNames,enujjIdNames,enujjTrigAndIdEffGraphs,False)
##MakeTrigAndIDEffMultigraphWithID('EGamma Medium',enujjTrigNames,enujjIdNames,enujjTrigAndIdEffGraphs,False)


### wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
#if __name__ == '__main__':
#   rep = ''
#   while not rep in [ 'q', 'Q' ]:
#      rep = raw_input( 'enter "q" to quit: ' )
#      if 1 < len(rep):
#         rep = rep[0]

exit()
