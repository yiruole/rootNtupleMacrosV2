#!/usr/bin/env python

#
# To combine the root files (inside the LQDATA/blah directory):
# for i in `find . -maxdepth 1 -type d`; do cd $i; echo $PWD; rm combined.root; hadd combined.root *.root; cd $LQDATA/trigEffStudyEejjSignals; done;
#
from ROOT import *
import os


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



#####################################################
# RUN
#####################################################

lqdata = os.environ['LQDATA']
analysisDir = lqdata+'/trigIDSpring15EffStudyEejjSignals/'

# diff trigs
inputFilesTrigs = []
inputFilesTrigs.append(analysisDir+'output_cutTable_lq1_eejj_trigger1_id1_effiStudy/combined.root')
inputFilesTrigs.append(analysisDir+'output_cutTable_lq1_eejj_trigger2_id1_effiStudy/combined.root')
inputFilesTrigs.append(analysisDir+'output_cutTable_lq1_eejj_trigger3_id1_effiStudy/combined.root')
inputFilesIDs = []
# diff IDs
inputFilesIDs.append(analysisDir+'output_cutTable_lq1_eejj_trigger1_id1_effiStudy/combined.root')
inputFilesIDs.append(analysisDir+'output_cutTable_lq1_eejj_trigger3_id1_effiStudy/combined.root')
inputFilesIDs.append(analysisDir+'output_cutTable_lq1_eejj_trigger1_id2_effiStudy/combined.root')
inputFilesIDs.append(analysisDir+'output_cutTable_lq1_eejj_trigger3_id2_effiStudy/combined.root')
inputFilesIDs.append(analysisDir+'output_cutTable_lq1_eejj_trigger1_id3_effiStudy/combined.root')
inputFilesIDs.append(analysisDir+'output_cutTable_lq1_eejj_trigger3_id3_effiStudy/combined.root')
inputFilesIDs.append(analysisDir+'output_cutTable_lq1_eejj_trigger1_id4_effiStudy/combined.root')
inputFilesIDs.append(analysisDir+'output_cutTable_lq1_eejj_trigger3_id4_effiStudy/combined.root')
inputFilesIDs.append(analysisDir+'output_cutTable_lq1_eejj_trigger1_id5_effiStudy/combined.root')
inputFilesIDs.append(analysisDir+'output_cutTable_lq1_eejj_trigger3_id5_effiStudy/combined.root')


trigEffGraphs = []
for fileName in inputFilesTrigs:
  tfile = TFile.Open(fileName)
  print '[1] looking in file:',fileName

  nevtsTwoGenEleEcalHist = tfile.Get('NEventsTwoGenEleECAL')
  nevtsTwoGenEleEcalWithTrigHist = tfile.Get('NEventsTwoGenEleECALPlusTrig')

  trigEff = TGraphAsymmErrors()
  trigEff.BayesDivide(nevtsTwoGenEleEcalWithTrigHist,nevtsTwoGenEleEcalHist)
  trigEff.GetXaxis().SetTitle('LQ mass [GeV]')
  trigEff.GetXaxis().SetNdivisions(512)
  trigEff.GetYaxis().SetTitle('TrigEff.')
  trigEff.GetYaxis().SetTitleOffset(0.7)
  trigEff.GetYaxis().SetNdivisions(514)
  #t3 = TCanvas()
  #t3.cd()
  #trigEff.Draw('ap')
  trigEffGraphs.append(trigEff)

  tfile.Close()


trigNames = []
idNames = []
idEffGraphs = []
idEffVsPt1stGenEleGraphs = []
for i,fileName in enumerate(inputFilesIDs):
  tfile = TFile.Open(fileName)
  print '[2] looking in file:',fileName

  nevtsTwoGenEleEcalWithTrigHist = tfile.Get('NEventsTwoGenEleECALPlusTrig')
  nevtsTwoGenEleEcalWithTrigWithTwoIDHist = tfile.Get('NEventsTwoGenEleECALPlusTrigPlusTwoID')
  
  idEff = TGraphAsymmErrors()
  idEff.BayesDivide(nevtsTwoGenEleEcalWithTrigWithTwoIDHist,nevtsTwoGenEleEcalWithTrigHist)
  idEff.GetXaxis().SetTitle('LQ mass [GeV]')
  idEff.GetXaxis().SetNdivisions(512)
  idEff.GetYaxis().SetTitle('Eff. [eejj]')
  idEff.GetYaxis().SetTitleOffset(0.9)
  idEff.GetYaxis().SetNdivisions(514)
  #t4 = TCanvas()
  #t4.cd()
  #idEff.Draw('ap')
  idEffGraphs.append(idEff)
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
  if i%2==0:
    trigNames.append('HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1')
  else:
    trigNames.append('HLT_Ele27_WP85_Gsf_v1')
  if i < 2:
    idNames.append('HEEPv6.0')
  elif i < 4:
    idNames.append('HEEPv5.1')
  elif i < 6:
    idNames.append('Egamma Tight')
  elif i < 8:
    idNames.append('Egamma Medium')
  else:
    idNames.append('Egamma Loose')

  tfile.Close()


# trig eff multigraph
mg = TMultiGraph()
leg = TLegend(0.375,0.178,0.852,0.451)
for i,graph in enumerate(trigEffGraphs):
  if i==0:
    graph.SetLineColor(2)
    graph.SetMarkerColor(2)
    leg.AddEntry(graph,'HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1','lp')
  elif i==1:
    graph.SetLineColor(4)
    graph.SetMarkerColor(4)
    leg.AddEntry(graph,'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v1','lp')
  elif i==2:
    graph.SetLineColor(8)
    graph.SetMarkerColor(8)
    leg.AddEntry(graph,'HLT_Ele27_WP85_Gsf_v1','lp')
  mg.Add(graph)
t5 = TCanvas()
t5.cd()
mg.Draw('ap')
mg.GetXaxis().SetTitle('LQ mass [GeV]')
mg.GetXaxis().SetNdivisions(512)
mg.GetYaxis().SetTitle('Efficiency [eejj]')
mg.GetYaxis().SetTitleOffset(0.7)
mg.GetYaxis().SetNdivisions(514)
mg.Draw('ap')
leg.SetBorderSize(0)
leg.Draw()
t5.Print('triggerEff.png')


# ID eff multigraph: trig 1
mg2 = TMultiGraph()
leg2 = TLegend(0.375,0.178,0.852,0.451)
leg2.SetHeader('HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1')
counter = 0
for i,graph in enumerate(idEffGraphs):
  if trigNames[i] != 'HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1':
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
mg2.GetYaxis().SetTitle('Efficiency [eejj]')
mg2.GetYaxis().SetTitleOffset(0.8)
mg2.GetYaxis().SetNdivisions(514)
mg2.GetYaxis().SetRangeUser(0.3,0.9)
mg2.Draw('ap')
leg2.SetBorderSize(0)
leg2.Draw()
t6.Print('idEff_trig1.png')


# ID eff multigraph: trig 3
mg3 = TMultiGraph()
leg3 = TLegend(0.375,0.178,0.852,0.451)
leg3.SetHeader('HLT_Ele27_WP85_Gsf_v1')
counter = 0
for i,graph in enumerate(idEffGraphs):
  if trigNames[i] != 'HLT_Ele27_WP85_Gsf_v1':
    continue
  graph = SetGraphStyle(graph,counter)
  leg3.AddEntry(graph,idNames[i],'lp')
  mg3.Add(graph)
  counter+=1
t7 = TCanvas()
t7.cd()
mg3.Draw('ap')
mg3.GetXaxis().SetTitle('LQ mass [GeV]')
mg3.GetXaxis().SetNdivisions(512)
mg3.GetYaxis().SetTitle('Efficiency [eejj]')
mg3.GetYaxis().SetTitleOffset(0.8)
mg3.GetYaxis().SetNdivisions(514)
mg2.GetYaxis().SetRangeUser(0.3,0.9)
mg3.Draw('ap')
leg3.SetTextSize(0.025)
leg3.SetBorderSize(0)
leg3.Draw()
t7.Print('idEff_trig3.png')


# ID eff 1st ele genpt multigraph: trig 1
mg4 = TMultiGraph()
leg4 = TLegend(0.375,0.178,0.852,0.451)
leg4.SetHeader('HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1')
counter = 0
for i,graph in enumerate(idEffVsPt1stGenEleGraphs):
  if trigNames[i] != 'HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1':
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
mg4.GetYaxis().SetTitle('Efficiency [eejj]')
mg4.GetYaxis().SetTitleOffset(0.8)
mg4.GetYaxis().SetNdivisions(514)
mg4.GetYaxis().SetRangeUser(0.3,0.9)
mg4.Draw('ap')
leg4.SetBorderSize(0)
leg4.Draw()
t8.Print('idEffVs1stEleGenPt_trig1.png')


## wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
if __name__ == '__main__':
   rep = ''
   while not rep in [ 'q', 'Q' ]:
      rep = raw_input( 'enter "q" to quit: ' )
      if 1 < len(rep):
         rep = rep[0]

exit()
