#!/usr/bin/env python

from ROOT import TGraphErrors,TCanvas,TLegend,TLine,kBlue,TGraph
import numpy as numpy
import math
import os

#txtFilename = os.getenv('LQANA')+'/versionsOfAnalysis_eejj/nov24_muonVeto35GeV/unscaled/rescale.log'
#txtFilename = os.getenv('LQANA')+'/versionsOfAnalysis_eejj/feb1/unscaled/rescale.log'
#txtFilename = os.getenv('LQANA')+'/versionsOfAnalysis_eejj/feb6/unscaled/test/rescale.log'
#txtFilename = os.getenv('LQANA')+'/versionsOfAnalysis_eejj/feb20/unscaled/rescale.log'
txtFilename = os.getenv('LQANA')+'/versionsOfAnalysis_eejj/mar17/unscaled/rescale.log'

#txtFilename = os.getenv('LQANA')+'/versionsOfAnalysis_enujj/feb5_SF_finalSels/unscaled/rescale.log'
#txtFilename = os.getenv('LQANA')+'/versionsOfAnalysis_enujj/feb20/unscaled/rescaleFinalSels.log'
#txtFilename = os.getenv('LQANA')+'/versionsOfAnalysis_enujj/mar6_noTopPtReweight/unscaled/rescale_finalSelections.log'
#txtFilename = os.getenv('LQANA')+'/versionsOfAnalysis_enujj/mar17/unscaled/rescale_lqBjetVariation.log'
#txtFilename = os.getenv('LQANA')+'/versionsOfAnalysis_enujj/mar17/unscaled/newSingleTop/rescaleVariations.log'
txtFilename = os.getenv('LQANA')+'/versionsOfAnalysis_enujj/mar17/unscaled/newSingleTop/rescale_stVariation.log'

bkgName = 'zjet'
#varList = ['sT','MejMin']
#txtFilename = 'ttbar/logAllVariations.log'
#bkgName = 'ttbar'
# use sT and MejMin, one at a time
#varList = ['sT']
#varList = ['MejMin']
varList = ['LQ']

# for enujj
#bkgName = 'wjet'
bkgName = 'ttbar'
#varList = ['LQ']
varList = ['sT']

eejjMode = True
if 'enujj' in txtFilename:
  eejjMode=False
  # in enujj mode, default x/y is ttbar

for var in varList:
  xPoints = []
  xPointErrs = []
  yPoints = []
  yPointErrs = []
  if not eejjMode:
    xPointsWJ = []
    xPointErrsWJ = []
    yPointsWJ = []
    yPointErrsWJ = []
  with open(txtFilename,'r') as thisFile:
    isVarPoint = False
    for line in thisFile:
      line = line.strip()
      middle = 0.0
      delta = 0.0
      rescale = 0.0
      rescaleErr = 0.0
      if 'name' in line:
        line = line.split(':')
        name = line[1].strip()
        if not var in name:
          continue
        isVarPoint = True
        if not 'LQ' in name:
          print 'name=',name
          if not 'To' in name:
            # doesn't quite work for Menu; use separate script for plotting MT variations
            split = name.split('_')
            lowerVarBound = float(split[1])
            upperVarBound = float(split[2])
            middle = float((lowerVarBound+upperVarBound)/2.0)
            delta = upperVarBound-middle
            varBin = name[name.find(var)+len(var):name.rfind('_')]
            print 'lowerVarBound=',lowerVarBound,'upperVarBound=',upperVarBound
          else:
            #print 'varBin=',varBin
            varBin = name[name.find(var)+len(var):name.rfind('_')]
            lowerVarBound = float(varBin[0:varBin.find('To')])
            upperVarBound = float(varBin[varBin.find('To')+2:varBin.rfind('_')])
            # above is for bins; below is for threshold plots
            #lowerVarBound = float(varBin)-5
            #upperVarBound = float(varBin)+5
            if upperVarBound > 2000:
              upperVarBound = 2000
            print 'lowerVarBound=',lowerVarBound,'upperVarBound=',upperVarBound
            middle = float((lowerVarBound+upperVarBound)/2.0)
            delta = upperVarBound-middle
        else: # for LQ 'final selection' plots
          #print 'name=',name
          varBin = name[name.find('LQ')+2:name.rfind('_')]
          if len(varBin)<=0:
            varBin = name[name.find('LQ')+2:]
          #print 'varBin=',varBin
          middle = float(varBin)
          delta = 0
        if eejjMode or not 'WJet' in name:
          xPoints.append(middle)
          xPointErrs.append(delta)
        elif not eejjMode and 'WJet' in name:
          xPointsWJ.append(middle)
          xPointErrsWJ.append(delta)
      if 'rescale factor' in line and isVarPoint:
        line = line.split(':')
        numbers = line[1]
        #print 'numbers="'+numbers+'"'
        rescale = float(numbers[0:numbers.find('+\-')])
        rescaleErr = float(numbers[numbers.find('+\-')+4:])
        if eejjMode or not 'WJet' in name:
          print 'yPoints append: rescale =',rescale,'rescaleErr=',rescaleErr
          yPoints.append(rescale)
          if rescaleErr < 0:
            rescaleErr = math.fabs(rescaleErr)
          yPointErrs.append(rescaleErr)
        elif not eejjMode and 'WJet' in name:
          print 'yPointsWJ append: rescale =',rescale,'rescaleErr=',rescaleErr
          yPointsWJ.append(rescale)
          if rescaleErr < 0:
            rescaleErr = math.fabs(rescaleErr)
          yPointErrsWJ.append(rescaleErr)
        isVarPoint = False
 

  if len(xPoints) <= 0:
    print 'No data in logfile found for var:',var
    continue
  # make graph
  #print numpy.array(xPoints)
  #print numpy.array(xPointErrs)
  #print numpy.array(yPoints)
  #print numpy.array(yPointErrs)
  if bkgName=='wjet' and not eejjMode:
    xPoints = xPointsWJ
    xPointErrs = xPointErrsWJ
    yPoints = yPointsWJ
    yPointErrs = yPointErrsWJ
  canvas = TCanvas()
  canvas.cd()
  graph = TGraphErrors(len(xPoints),numpy.array(xPoints),numpy.array(yPoints),numpy.array(xPointErrs),numpy.array(yPointErrs))
  graph.SetTitle('')
  graph.SetMarkerColor(kBlue)
  graph.SetLineColor(kBlue)
  graph.Draw('ap')
  if var=='sT':
    graph.GetXaxis().SetTitle('S_{T} [GeV]')
  elif var=='MejMin':
    graph.GetXaxis().SetTitle('M_{ej}^{min} [GeV]')
  elif 'LQ' in var:
    graph.GetXaxis().SetTitle('M_{LQ} [GeV]')
  graph.GetXaxis().SetNdivisions(516)
  if eejjMode:
    if bkgName=='ttbar':
      graph.GetYaxis().SetTitle('t#bar{t} scale factor')
      sfNom = 0.81479
      sfNomErr = 0.036965
    elif bkgName=='zjet':
      graph.GetYaxis().SetTitle('Z+jets scale factor')
      #sfNom = 0.94 # MGHT
      #sfNomErr = 0.01
      #sfNom = 1.03 # amc@NLO PtBinned
      #sfNomErr = 0.02
      #sfNom = 1.05 # amc@NLO PtBinned
      #sfNomErr = 0.01
      #sfNom = 0.98 # amc@NLO PtBinned, Sep. 29 PtEE
      sfNom = 0.97 # amc@NLO PtBinned, Nov. 19, amcAtNLO diboson
      sfNomErr = 0.01
  else:
    if bkgName=='ttbar':
      graph.GetYaxis().SetTitle('t#bar{t} scale factor')
      graph.GetYaxis().SetRangeUser(0,2)
      #sfNom = 0.95
      sfNom = 0.92 # no top pt reweight
      sfNomErr = 0.01
    elif bkgName=='wjet':
      graph.GetYaxis().SetTitle('W+jets scale factor')
      graph.GetYaxis().SetRangeUser(0,2)
      sfNom = 0.87
      sfNomErr = 0.012
  uncertaintyXpoints = [xPoints[0]-xPointErrs[0],xPoints[-1]+xPointErrs[-1]]
  uncertaintyYpoints = [sfNom,sfNom]
  uncertaintyXpointsErrs = [0,0]
  uncertaintyYpointsErrs = [sfNomErr,sfNomErr]
  uncertaintyRegionGraph = TGraphErrors(len(uncertaintyXpoints),numpy.array(uncertaintyXpoints),numpy.array(uncertaintyYpoints),numpy.array(uncertaintyXpointsErrs),numpy.array(uncertaintyYpointsErrs))
  uncertaintyRegionGraph.SetFillColor(15)
  uncertaintyRegionGraph.SetFillStyle(3001)
  uncertaintyRegionGraph.SetLineColor(1)
  uncertaintyRegionGraph.Draw('3l')
  #lineSFNom = TLine(xPoints[0]-xPointErrs[0],sfNom,xPoints[0]-xPointErrs[0],sfNom)
  #lineSFNom.SetLineColor(2)
  #lineSFNom.Draw()
  # draw the polygon for shading
  #uncertaintyXpoints = [xPoints[0]-xPointErrs[0],xPoints[-1]+xPointErrs[-1],xPoints[-1]+xPointErrs[-1],xPoints[0]-xPointErrs[0]]
  #uncertaintyYpoints = [sfNom+sfNomErr,sfNom+sfNomErr,sfNom-sfNomErr,sfNom-sfNomErr]
  #uncertaintyRegionGraph = TGraph(len(uncertaintyXpoints),numpy.array(uncertaintyXpoints),numpy.array(uncertaintyYpoints))
  #uncertaintyRegionGraph.SetFillColor(15)
  #uncertaintyRegionGraph.SetFillStyle(3001)
  #uncertaintyRegionGraph.Draw('f')
  #lineUp = TLine(graph.GetXaxis().GetXmin(),sfNom*1.1,graph.GetXaxis().GetXmax(),sfNom*1.1)
  #lineUp.SetLineStyle(2)
  #lineUp.Draw()
  #lineDown = TLine(graph.GetXaxis().GetXmin(),sfNom*0.9,graph.GetXaxis().GetXmax(),sfNom*0.9)
  #lineDown.SetLineStyle(2)
  #lineDown.Draw()
  if eejjMode:
    if bkgName=='ttbar':
      if var=='sT':
        leg = TLegend(0.19,0.68,0.51,0.89)
      elif var=='MejMin':
        leg = TLegend(0.19,0.19,0.51,0.40)
      leg.AddEntry(graph,'t#bar{t} scale factor variations','lp')
    elif bkgName=='zjet':
      if var=='sT':
        #leg = TLegend(0.19,0.79,0.51,0.89)
        leg = TLegend(0.19,0.19,0.51,0.40)
      elif var=='MejMin':
        #leg = TLegend(0.19,0.79,0.51,0.89)
        leg = TLegend(0.19,0.19,0.51,0.40)
      else:
        leg = TLegend(0.19,0.19,0.51,0.40)
      leg.AddEntry(graph,'Z+jets scale factor variations','lp')
  else:
    if bkgName=='ttbar':
      if var=='sT':
        leg = TLegend(0.19,0.68,0.51,0.89)
      elif var=='MejMin':
        leg = TLegend(0.19,0.19,0.51,0.40)
      else:
        leg = TLegend(0.19,0.19,0.51,0.40)
      leg.AddEntry(graph,'t#bar{t} scale factor variations','lp')
    elif bkgName=='wjet':
      if var=='sT':
        #leg = TLegend(0.19,0.79,0.51,0.89)
        leg = TLegend(0.19,0.19,0.51,0.40)
      elif var=='MejMin':
        #leg = TLegend(0.19,0.79,0.51,0.89)
        leg = TLegend(0.19,0.19,0.51,0.40)
      else:
        leg = TLegend(0.19,0.68,0.51,0.89)
      leg.AddEntry(graph,'W+jets scale factor variations','lp')
  leg.AddEntry(uncertaintyRegionGraph,'Nominal scale factor = '+str(round(sfNom,3))+' #pm '+str(round(sfNomErr,3)),'fl')
  #leg.AddEntry(lineUp,'Nominal scale factor #pm 10%','l')
  leg.SetBorderSize(0)
  leg.Draw()
  canvas.Modified()

  # print
  if var=='sT':
    baseName = bkgName+'_scaleFactorVariation_sTBins' 
  elif var=='MejMin':
    baseName = bkgName+'_scaleFactorVariation_mejMinBins'
  elif 'LQ' in var:
    baseName = bkgName+'_scaleFactorVariation_LQBins'
  elif 'MTenu' in var:
    baseName = bkgName+'_scaleFactorVariation_MTBins'
  canvas.Print(baseName+'.png')
  canvas.Print(baseName+'.pdf')

  ## wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
  if __name__ == '__main__':
    rep = ''
    while not rep in [ 'c', 'C' ]:
     rep = raw_input( 'enter "c" to continue: ' )
     if 1 < len(rep):
       rep = rep[0]

exit(0)

