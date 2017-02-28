#!/usr/bin/env python

from ROOT import TGraphErrors,TCanvas,TLegend,TLine,kBlue,TGraph
import numpy as numpy
import math
import os

#txtFilename = os.getenv('LQANA')+'/versionsOfAnalysis_eejj/oct27_tupleV211/zjet_ht/log.log'
txtFilename = os.getenv('LQANA')+'/versionsOfAnalysis_eejj/feb28/logDYrescale.log'
bkgName = 'zjet'
#varList = ['sT','MejMin']
#txtFilename = 'ttbar/logAllVariations.log'
#bkgName = 'ttbar'
varList = ['sT']

for var in varList:
  xPoints = []
  xPointErrs = []
  yPoints = []
  yPointErrs = []
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
        #print 'name=',name
        varBin = name[name.find(var)+len(var):name.rfind('_')]
        #print 'varBin=',varBin
        lowerVarBound = float(varBin[0:varBin.find('To')])
        upperVarBound = float(varBin[varBin.find('To')+2:])
        # above is for bins; below is for threshold plots
        #lowerVarBound = float(varBin)-5
        #upperVarBound = float(varBin)+5
        if upperVarBound > 2000:
          upperVarBound = 2000
        #print 'lowerVarBound=',lowerVarBound,'upperVarBound=',upperVarBound
        middle = float((lowerVarBound+upperVarBound)/2.0)
        delta = upperVarBound-middle
        xPoints.append(middle)
        xPointErrs.append(delta)
      if 'rescale' in line and isVarPoint:
        line = line.split(':')
        numbers = line[1]
        #print 'numbers="'+numbers+'"'
        rescale = float(numbers[0:numbers.find('+\-')])
        rescaleErr = float(numbers[numbers.find('+\-')+4:])
        print 'rescale =',rescale,'rescaleErr=',rescaleErr
        yPoints.append(rescale)
        if rescaleErr < 0:
          rescaleErr = math.fabs(rescaleErr)
        yPointErrs.append(rescaleErr)
        isVarPoint = False
 
  if len(xPoints) <= 0:
    print 'No data in logfile found for var:',var
    continue
  # make graph
  #print numpy.array(xPoints)
  #print numpy.array(xPointErrs)
  #print numpy.array(yPoints)
  #print numpy.array(yPointErrs)
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
  graph.GetXaxis().SetNdivisions(516)
  if bkgName=='ttbar':
    graph.GetYaxis().SetTitle('t#bar{t} scale factor')
    sfNom = 0.81479
    sfNomErr = 0.036965
  elif bkgName=='zjet':
    graph.GetYaxis().SetTitle('Z+jets scale factor')
    #sfNom = 0.94 # MGHT
    #sfNomErr = 0.01
    sfNom = 0.95 # amc@NLO PtBinned
    sfNomErr = 0.02
  lineSFNom = TLine(graph.GetXaxis().GetXmin(),sfNom,graph.GetXaxis().GetXmax(),sfNom)
  lineSFNom.Draw()
  uncertaintyXpoints = [xPoints[0]-xPointErrs[0],xPoints[-1]+xPointErrs[-1],xPoints[-1]+xPointErrs[-1],xPoints[0]-xPointErrs[0]]
  uncertaintyYpoints = [sfNom+sfNomErr,sfNom+sfNomErr,sfNom-sfNomErr,sfNom-sfNomErr]
  uncertaintyRegionGraph = TGraph(len(uncertaintyXpoints),numpy.array(uncertaintyXpoints),numpy.array(uncertaintyYpoints))
  uncertaintyRegionGraph.SetFillColor(15)
  uncertaintyRegionGraph.SetFillStyle(3001)
  uncertaintyRegionGraph.Draw('f')
  #lineUp = TLine(graph.GetXaxis().GetXmin(),sfNom*1.1,graph.GetXaxis().GetXmax(),sfNom*1.1)
  #lineUp.SetLineStyle(2)
  #lineUp.Draw()
  #lineDown = TLine(graph.GetXaxis().GetXmin(),sfNom*0.9,graph.GetXaxis().GetXmax(),sfNom*0.9)
  #lineDown.SetLineStyle(2)
  #lineDown.Draw()
  if bkgName=='ttbar':
    if var=='sT':
      leg = TLegend(0.19,0.68,0.51,0.89)
    elif var=='MejMin':
      leg = TLegend(0.19,0.19,0.51,0.40)
    leg.AddEntry(graph,'t#bar{t} scale factor variations','lp')
  elif bkgName=='zjet':
    if var=='sT':
      leg = TLegend(0.19,0.19,0.51,0.40)
    elif var=='MejMin':
      leg = TLegend(0.19,0.79,0.51,0.89)
    leg.AddEntry(graph,'Z+jets scale factor variations','lp')
  leg.AddEntry(lineSFNom,'Nominal scale factor = '+str(round(sfNom,3))+' #pm '+str(round(sfNomErr,3)),'l')
  #leg.AddEntry(lineUp,'Nominal scale factor #pm 10%','l')
  leg.SetBorderSize(0)
  leg.Draw()
  canvas.Modified()

  # print
  baseName = bkgName+'_scaleFactorVariation_sTBins' if var=='sT' else bkgName+'_scaleFactorVariation_mejMinBins'
  canvas.Print(baseName+'.png')
  #canvas.Print(baseName+'.pdf')

  ## wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
  if __name__ == '__main__':
    rep = ''
    while not rep in [ 'c', 'C' ]:
     rep = raw_input( 'enter "c" to continue: ' )
     if 1 < len(rep):
       rep = rep[0]

exit(0)

