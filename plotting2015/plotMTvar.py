#!/usr/bin/env python

from ROOT import *
import numpy

mtLow = [70.0,110.0,150.0,200,400,]
mtHigh = [110.0,150.0,190.0,400,600,]
mtMiddle = [90.0,130.0,170.0,300,500,]
mtDelta = [20.0,20.0,20.0,100,100,]

#wJetsSF = [0.87,0.77,0.68]
#wJetsSFErr = [0.01,0.03,0.08]
# no top PtReweight Mar6 2018
#wJetsSF = [0.87,0.77,0.68]
#wJetsSFErr = [0.013,0.03,0.07]
# fix muons Mar17 2018
#wJetsSF = [0.870390454678,0.771572269699,0.707961636976]
#wJetsSFErr = [0.0113683134421,0.0285582518664,0.0805988823919]
wJetsSF = [0.869075214034,0.77094051958,0.703717693798,0.793739236228,0.784642078141]
wJetsSFErr = [0.0117150147512,0.0297768124288,0.0807847584029,0.0589450751429,0.0631985710737]

#ttSF = [0.96,0.93,0.95]
#ttSFErr = [0.01,0.01,0.02]
# no top PtReweight Mar6 2018
#ttSF = [0.93,0.91,0.92]
#ttSFErr = [0.01,0.014,0.02]
# fix muons Mar17 2018
#ttSF = [0.927925925574,0.905547940866,0.915299946222]
#ttSFErr = [0.00990741202651,0.0130330847157,0.0181914543556]
ttSF = [0.928068644323,0.903234385133,0.915210492272,0.883567955986,1.03237742053,]
ttSFErr = [0.00865079550222,0.013136462854,0.0200283014143,0.0159271814383,0.178092382728]

#wjetsSFNom = 0.87
#wjetsSFNomErr = 0.01
# no top PtReweight Mar6 2018
wjetsSFNom = 0.87
wjetsSFNomErr = 0.014
# add from Mar17
ttbarSFNom = 0.92
ttbarSFNomErr = 0.0096

canvas = TCanvas()
canvas.cd()

#print len(wJetsSF),numpy.array(mtMiddle),numpy.array(wJetsSF),numpy.array(mtDelta),numpy.array(wJetsSFErr)
graphWJets = TGraphErrors(len(wJetsSF),numpy.array(mtMiddle),numpy.array(wJetsSF),numpy.array(mtDelta),numpy.array(wJetsSFErr))
graphWJets.SetTitle('')
graphWJets.SetMarkerColor(kBlue)
graphWJets.SetLineColor(kBlue)
graphWJets.Draw('ap')
graphWJets.GetXaxis().SetTitle('M_{T}(e#nu) [GeV]')

uncertaintyXpoints = [mtMiddle[0]-mtDelta[0],mtMiddle[-1]+mtDelta[-1]]
uncertaintyYpoints = [wjetsSFNom,wjetsSFNom]
uncertaintyXpointsErrs = [0,0]
uncertaintyYpointsErrs = [wjetsSFNomErr,wjetsSFNomErr]
#print len(uncertaintyXpoints),numpy.array(uncertaintyXpoints),numpy.array(uncertaintyYpoints),numpy.array(uncertaintyXpointsErrs),numpy.array(uncertaintyYpointsErrs)
uncertaintyRegionGraph = TGraphErrors(len(uncertaintyXpoints),numpy.array(uncertaintyXpoints),numpy.array(uncertaintyYpoints),numpy.array(uncertaintyXpointsErrs),numpy.array(uncertaintyYpointsErrs))
uncertaintyRegionGraph.SetFillColor(15)
uncertaintyRegionGraph.SetFillStyle(3001)
uncertaintyRegionGraph.SetLineColor(1)
uncertaintyRegionGraph.Draw('3l')
graphWJets.Draw('psame')

leg = TLegend(0.19,0.19,0.51,0.40)
leg.AddEntry(graphWJets,'W+jets scale factor variations','lp')
leg.AddEntry(uncertaintyRegionGraph,'Nominal scale factor = '+str(round(wjetsSFNom,3))+' #pm '+str(round(wjetsSFNomErr,3)),'fl')
#leg.AddEntry(lineUp,'Nominal scale factor #pm 10%','l')
leg.SetBorderSize(0)
leg.Draw()
canvas.Modified()
canvas.Print('wjet_scaleFactorVariation_MTBins.pdf')


canvas2 = TCanvas()
canvas2.cd()
graphTT = TGraphErrors(len(ttSF),numpy.array(mtMiddle),numpy.array(ttSF),numpy.array(mtDelta),numpy.array(ttSFErr))
graphTT.SetTitle('')
graphTT.SetMarkerColor(kBlue)
graphTT.SetLineColor(kBlue)
graphTT.Draw('ap')
graphTT.GetXaxis().SetTitle('M_{T}(e#nu) [GeV]')

uncertaintyTTXpoints = [mtMiddle[0]-mtDelta[0],mtMiddle[-1]+mtDelta[-1]]
uncertaintyTTYpoints = [ttbarSFNom,ttbarSFNom]
uncertaintyTTXpointsErrs = [0,0]
uncertaintyTTYpointsErrs = [ttbarSFNomErr,ttbarSFNomErr]
uncertaintyRegionTTGraph = TGraphErrors(len(uncertaintyTTXpoints),numpy.array(uncertaintyTTXpoints),numpy.array(uncertaintyTTYpoints),numpy.array(uncertaintyTTXpointsErrs),numpy.array(uncertaintyTTYpointsErrs))
uncertaintyRegionTTGraph.SetFillColor(15)
uncertaintyRegionTTGraph.SetFillStyle(3001)
uncertaintyRegionTTGraph.SetLineColor(1)
uncertaintyRegionTTGraph.Draw('3l')
graphTT.Draw('psame')

leg2 = TLegend(0.19,0.19,0.51,0.40)
leg2.AddEntry(graphTT,'t#bar{t} scale factor variations','lp')
leg2.AddEntry(uncertaintyRegionGraph,'Nominal scale factor = '+str(round(ttbarSFNom,3))+' #pm '+str(round(ttbarSFNomErr,3)),'fl')
#leg2.AddEntry(lineUp,'Nominal scale factor #pm 10%','l')
leg2.SetBorderSize(0)
leg2.Draw()
canvas2.Modified()
#canvas2.Print('ttbar_scaleFactorVariation_MTBins.pdf')
