import os
import sys
import math

import ROOT
# import subprocess as sp
import array

ROOT.gStyle.SetOptStat(0)
ROOT.gErrorIgnoreLevel = 1
# ROOT.gROOT.SetBatch(True);


def fitf(x, par):

    xx = x[0]

    if xx < par[0]:
        fitval = par[1] + par[2] * xx
    else:
        fitval = par[1] + par[2] * par[0]

    return fitval


def getNumeratorHist(root_file, variable, ptMin, ptMax, rebin):
    # data_hist = root_file.Get("histo1D__DATA__Pass_"   + variable)
    # data_hist.Rebin ( rebin )
    # mc_hist   = root_file.Get("histo1D__ALLBKG__Pass_" + variable)
    # mc_hist.Rebin ( rebin )
    # print "Numerator: "
    # print "\tData integral =", data_hist.Integral(0, data_hist.GetNbinsX() + 1 )
    # print "\tMC   integral =", mc_hist  .Integral(0, mc_hist  .GetNbinsX() + 1 )
    # print "\tcontamination =", 100. * mc_hist  .Integral(0, mc_hist  .GetNbinsX() + 1 ) / ( mc_hist  .Integral(0, mc_hist  .GetNbinsX() + 1 ) + data_hist.Integral(0, data_hist.GetNbinsX() + 1 ) )
    # hist = data_hist.Clone()
    # hist.Add ( mc_hist, -1 )

    jets_hist2D = root_file.Get("histo2D__QCDFakes_DATA__Jets_" + variable)
    jets_hist1D = jets_hist2D.ProjectionY(
        "Jets_" + variable + "_1D",
        jets_hist2D.GetXaxis().FindBin(ptMin),
        jets_hist2D.GetXaxis().FindBin(ptMax),
    )
    hist = jets_hist1D.Clone()
    return hist


def getDenominatorHist(root_file, variable, ptMin, ptMax, rebin):
    total_hist2D = root_file.Get("histo2D__QCDFakes_DATA__Total_" + variable)
    total_hist1D = total_hist2D.ProjectionY(
        "Total_" + variable + "_1D",
        total_hist2D.GetXaxis().FindBin(ptMin),
        total_hist2D.GetXaxis().FindBin(ptMax),
    )
    hist = total_hist1D.Clone()
    return hist


def getElectronHist(root_file, variable, ptMin, ptMax, rebin):
    total_hist2D = root_file.Get("histo2D__QCDFakes_DATA__Electrons_" + variable)
    total_hist1D = total_hist2D.ProjectionY(
        "Electrons_" + variable + "_1D",
        total_hist2D.GetXaxis().FindBin(ptMin),
        total_hist2D.GetXaxis().FindBin(ptMax),
    )
    hist = total_hist1D.Clone()
    return hist


def makePlots(
    psfile, canvas, root_file, variable, xtitle, xmin, xmax, ymin, ymax, rebin
):

    canvas.Clear()
    canvas.cd()

    ptBins = [
        35,
        45,
        55,
        65,
        75,
        85,
        95,
        105,
        125,
        145,
        165,
        185,
        205,
        250,
        300,
        400,
        500,
        600,
        700,
        800,
        900,
        1000,
    ]
    ptBinsArr = array.array("f", ptBins)
    # ratio = TGraphErrors()
    ratio = ROOT.TH1F("fakeRate", "", len(ptBins) - 1, ptBinsArr)
    for index, pt in enumerate(ptBins):
        if index >= len(ptBins) - 1:
            break
        # print 'len=',len(ptBins),'index=',index
        ptMin = pt
        ptMax = ptBins[index + 1]
        n_hist = getNumeratorHist(root_file, variable, ptMin, ptMax, rebin)
        e_hist = getElectronHist(root_file, variable, ptMin, ptMax, rebin)
        # TrkIso sideband: 10-20 GeV
        # TrkIso signal: < 5 GeV
        # take ratio of jets_SR/jets_SB
        jets_SR_err = ROOT.Double(0)
        jets_SB_err = ROOT.Double(0)
        jets_SR = n_hist.IntegralAndError(
            n_hist.FindBin(0), n_hist.FindBin(5), jets_SR_err
        )
        jets_SB = n_hist.IntegralAndError(
            n_hist.FindBin(10), n_hist.FindBin(20), jets_SB_err
        )
        ratioJets_SRSB = jets_SR / jets_SB
        ratioJets_SRSB_err = ratioJets_SRSB * math.sqrt(
            pow(jets_SB_err / jets_SB, 2) + pow(jets_SR_err / jets_SR, 2)
        )
        # multiply by HEEP electrons in SB
        nHEEP_SB_err = ROOT.Double(0)
        nHEEP_SB = e_hist.IntegralAndError(
            n_hist.FindBin(10), n_hist.FindBin(20), nHEEP_SB_err
        )
        nJets_sig = ratioJets_SRSB * nHEEP_SB
        nJets_sig_err = nJets_sig * math.sqrt(
            pow(ratioJets_SRSB_err / ratioJets_SRSB, 2)
            + pow(nHEEP_SB_err / nHEEP_SB, 2)
        )
        d_hist = getDenominatorHist(root_file, variable, ptMin, ptMax, rebin)
        d_int_err = ROOT.Double(0)
        d_int = d_hist.IntegralAndError(0, d_hist.GetNbinsX(), d_int_err)
        fake_rate = nJets_sig / d_int
        fake_rate_err = fake_rate * math.sqrt(
            pow(nJets_sig_err / nJets_sig, 2) + pow(d_int_err / d_int, 2)
        )
        print "pt:", ptMin, "-", ptMax, "\tfake_rate=", fake_rate, "+/-", fake_rate_err
        ratio.SetBinContent(index + 1, fake_rate)
        ratio.SetBinError(index + 1, fake_rate_err)

    ratio.Draw("p")
    ratio.GetYaxis().SetRangeUser(0, 0.05)
    ratio.GetYaxis().SetNdivisions(515)
    ratio.GetXaxis().SetNdivisions(515)
    canvas.SetGridy()
    canvas.Update()
    canvas.Print("fakeRate_Jets0_barrel.png")
    # eventually use the loop to fill the fake rate hist or graph
    # r_hist = n_hist.Clone()
    # r_hist.Divide ( d_hist )

    # d_hist.GetXaxis().SetRangeUser ( xmin, xmax )
    # n_hist.GetXaxis().SetRangeUser ( xmin, xmax )
    # r_hist.GetXaxis().SetRangeUser ( xmin, xmax )
    #
    # d_hist.SetMaximum ( d_hist.GetMaximum() * 10000 )
    # n_hist.SetMaximum ( d_hist.GetMaximum() * 10000 )
    # d_hist.SetMinimum ( 0.1 )
    # n_hist.SetMinimum ( 0.1 )

    # r_hist.SetMaximum ( ymax )
    # r_hist.SetMinimum ( ymin )

    # for ibin in range ( 0, r_hist.GetNbinsX() + 2):
    #    error_bar_length = r_hist.GetBinError ( ibin )
    #    mean_value = r_hist.GetBinContent ( ibin )
    #    if error_bar_length < 0.02: continue

    #    bin_min = r_hist.GetXaxis().GetBinLowEdge(ibin)
    #    bin_max = r_hist.GetXaxis().GetBinUpEdge (ibin)

    #    if bin_max > xmax : continue

    # d_hist.GetXaxis().SetTitle ( xtitle )
    # d_hist.GetYaxis().SetTitleOffset ( 1.5 )
    ##d_hist.GetYaxis().SetTitle ( "Data - non-QCD MC, width = " + str ( d_hist.GetBinWidth(1)) + " GeV" )

    # d_hist.SetLineColor  (ROOT.kRed)
    # d_hist.SetMarkerColor(ROOT.kRed)
    # n_hist.SetLineColor  (ROOT.kBlue)
    # n_hist.SetMarkerColor(ROOT.kBlue)

    # hsize=0.60
    # vsize=0.15
    # xstart=0.30
    # ystart=0.75

    # legend = ROOT.TLegend(xstart, ystart, xstart+hsize, ystart+vsize)
    # legend.SetFillColor(ROOT.kWhite)
    # legend.SetBorderSize(0)
    # legend.SetShadowColor(10)
    # legend.SetMargin(0.2)
    # legend.SetTextFont(132)
    ##legend.AddEntry ( d_hist,"Data - non-QCD MC", "pl" )
    ##legend.AddEntry ( n_hist,"Data - non-QCD MC, passing HEEP v4.1", "pl" )

    # canvas.Divide ( 2, 1)
    # pad1 = canvas.GetPad ( 1 )
    # pad2 = canvas.GetPad ( 2 )
    # pad1.cd()
    # pad1.SetLogy()
    # d_hist.Draw()
    # n_hist.Draw("SAME")
    # legend.Draw()

    # pad2.cd()
    # r_hist.GetXaxis().SetTitle ( xtitle )
    # r_hist.GetYaxis().SetTitle ( "Fake rate" )
    # r_hist.GetYaxis().SetTitleOffset ( 1.5 )
    # r_hist.Draw()

    # ROOT.gErrorIgnoreLevel = 1
    # canvas.Print (psfile)
    # canvas.SaveAs ( "fakeRate_eps/" + variable + ".eps"  );
    # canvas.SaveAs ( "fakeRate_gif/" + variable + ".gif"  );

    # return parameters, par_errors


# -----------------------------------------------------------------
# Get the ROOT file
# -----------------------------------------------------------------

if len(sys.argv) != 2:
    print "Usage: python calculateFakeRate.py <root file>"
    sys.exit()

if not os.path.isfile(sys.argv[1]):
    print "Usage: python calculateFakeRate.py <root file>"
    print sys.argv[1] + " is not a file "
    sys.exit()

root_file_path = sys.argv[1]

root_file = ROOT.TFile(root_file_path)

# -----------------------------------------------------------------
# Make the plots
# -----------------------------------------------------------------

plot_file_name = "plots.ps"
plot_file_pdf = plot_file_name.replace(".ps", ".pdf")

canvas = ROOT.TCanvas("canvas", "", 800, 500)

ROOT.gErrorIgnoreLevel = 1
# canvas.Print (plot_file_name + "[")

makePlots(
    plot_file_name,
    canvas,
    root_file,
    "Bar_TrkIsoHEEP7vsPt_PAS",
    "1st Ele p_{T} (Barrel) [GeV]",
    35,
    500,
    0.0,
    0.03,
    10,
)
# par_inc_bar , err_inc_bar  = makePlots ( plot_file_name, canvas, root_file, "Bar_Pt1stEle_PAS" , "1st Ele p_{T} (Barrel) [GeV]"              , 35, 500, 0.0, 0.03, 10 )
# par_inc_end1, err_inc_end1 = makePlots ( plot_file_name, canvas, root_file, "End1_Pt1stEle_PAS", "1st Ele p_{T} (Endcap, |#eta| < 2.0) [GeV]", 35, 500, 0.0, 0.12, 10 )
# par_inc_end2, err_inc_end2 = makePlots ( plot_file_name, canvas, root_file, "End2_Pt1stEle_PAS", "1st Ele p_{T} (Endcap, |#eta| > 2.0) [GeV]", 35, 500, 0.0, 0.12, 10 )
#
# par_1jet_bar , err_1jet_bar  = makePlots ( plot_file_name, canvas, root_file, "Bar_1Jet_Pt1stEle_PAS" , "1st Ele p_{T} (Barrel, 1 Jet) [GeV]"                   , 35, 500, 0.0, 0.12, 10 )
# par_1jet_end1, err_1jet_end1 = makePlots ( plot_file_name, canvas, root_file, "End1_1Jet_Pt1stEle_PAS", "1st Ele p_{T} (Endcap, |#eta| < 2.0, #geq 1 Jet) [GeV]", 35, 500, 0.0, 0.12, 10 )
# par_1jet_end2, err_1jet_end2 = makePlots ( plot_file_name, canvas, root_file, "End2_1Jet_Pt1stEle_PAS", "1st Ele p_{T} (Endcap, |#eta| > 2.0, #geq 1 Jet) [GeV]", 35, 500, 0.0, 0.12, 10 )
#
# par_2jet_bar , err_2jet_bar  = makePlots ( plot_file_name, canvas, root_file, "Bar_2Jet_Pt1stEle_PAS" , "1st Ele p_{T} (Barrel, 2 Jets) [GeV]"                    , 35, 500, 0.0, 0.03, 10 )
# par_2jet_end1, err_2jet_end1 = makePlots ( plot_file_name, canvas, root_file, "End1_2Jet_Pt1stEle_PAS", "1st Ele p_{T} (Endcap, |#eta| < 2.0, #geq 2 Jets ) [GeV]", 35, 500, 0.0, 0.12, 10 )
# par_2jet_end2, err_2jet_end2 = makePlots ( plot_file_name, canvas, root_file, "End2_2Jet_Pt1stEle_PAS", "1st Ele p_{T} (Endcap, |#eta| > 2.0, #geq 2 Jets ) [GeV]", 35, 500, 0.0, 0.12, 10 )
#
# par_3jet_bar , err_3jet_bar  = makePlots ( plot_file_name, canvas, root_file, "Bar_3Jet_Pt1stEle_PAS" , "1st Ele p_{T} (Barrel, 3 Jets) [GeV]"                    , 35, 500, 0.0, 0.12, 10 )
# par_3jet_end1, err_3jet_end1 = makePlots ( plot_file_name, canvas, root_file, "End1_3Jet_Pt1stEle_PAS", "1st Ele p_{T} (Endcap, |#eta| < 2.0, #geq 3 Jets ) [GeV]", 35, 500, 0.0, 0.12, 10 )
# par_3jet_end2, err_3jet_end2 = makePlots ( plot_file_name, canvas, root_file, "End2_3Jet_Pt1stEle_PAS", "1st Ele p_{T} (Endcap, |#eta| > 2.0, #geq 3 Jets ) [GeV]", 35, 500, 0.0, 0.12, 10 )

# ROOT.gErrorIgnoreLevel = 1
# canvas.Print (plot_file_name + "]")
#
# print "Converting file"
#
# os.system ( "ps2pdf " + plot_file_name )
# os.system ( "rm " + plot_file_name )
#
# print "...Done"

# TODO
##-----------------------------------------------------------------
## Make the table
##-----------------------------------------------------------------
#
# print "Writing latex"
#
# table_file_name = "table.tex"
# table_file = open ( table_file_name, "w" )
# table_file_pdf = table_file_name.replace (".tex",".pdf")
#
# table_file.write("\documentclass{article}\n")
# table_file.write("\usepackage{multirow} \n")
# table_file.write("\usepackage{amsmath} \n")
# table_file.write("\usepackage[landscape, top=1cm, bottom=1cm, left=1cm, right=1cm]{geometry}\n")
# table_file.write("\pagestyle{empty} \n")
# table_file.write("\\begin{document}\n")
# table_file.write("\\begin{table}[h]\n")
# table_file.write("\\begin{tabular}{| c | c | l |}\n")
# table_file.write("  \hline \n")
# table_file.write("jet requirement & $\eta$ region & functional form \\\ \n")
#
# table_file.write("  \hline \n")
#
# if par_inc_bar[1] < 0: par2_string = " - "
# else                  : par2_string = " + "
#
# table_file.write("& $|\eta| < $ 1.442 &  \\text{for} $E_{T} < $" + " %3.1f:" % par_inc_bar[0] + " %1.3f" % par_inc_bar[1] + par2_string + " %1.5f" % abs(par_inc_bar[2]) + " $\\times E_{T} \\text{, else }$ %1.4f" % ( par_inc_bar[1] + ( par_inc_bar[2] * par_inc_bar[0] )) +  "\\\ \n")
#
# if par_inc_end1[1] < 0: par2_string = " - "
# else                   : par2_string = " + "
#
# table_file.write("$\\text{N(jet)} \geq 0 $ & 1.56 $< |\eta| < $ 2.0 &  \\text{for} $E_{T} < $" + " %3.1f:" % par_inc_end1[0] +  " %1.3f" % par_inc_end1[1] + par2_string + " %1.5f" % abs(par_inc_end1[2]) + " $\\times E_{T} \\text{, else }$ %1.4f" % ( par_inc_end1[1] + ( par_inc_end1[2] * par_inc_end1[0] )) +  "\\\ \n")
#
# if par_inc_end2[1] < 0: par2_string = " - "
# else                   : par2_string = " + "
#
# table_file.write("& 2.0 $< |\eta| < $ 2.5 &  \\text{for} $E_{T} < $" + " %3.1f:" % par_inc_end2[0] +  " %1.3f" % par_inc_end2[1] + par2_string + " %1.5f" % abs(par_inc_end2[2]) + " $\\times E_{T} \\text{, else }$ %1.4f" % ( par_inc_end2[1] + ( par_inc_end2[2] * par_inc_end2[0] )) +  "\\\ \n")
#
# table_file.write("  \hline \n")
#
# if par_1jet_bar[1] < 0: par2_string = " - "
# else                  : par2_string = " + "
#
# table_file.write("& $|\eta| < $ 1.442 &  \\text{for} $E_{T} < $" + " %3.1f:" % par_1jet_bar[0] + " %1.3f" % par_1jet_bar[1] + par2_string + " %1.5f" % abs(par_1jet_bar[2]) + " $\\times E_{T} \\text{, else }$ %1.4f" % ( par_1jet_bar[1] + ( par_1jet_bar[2] * par_1jet_bar[0] )) +  "\\\ \n")
#
# if par_1jet_end1[1] < 0: par2_string = " - "
# else                   : par2_string = " + "
#
# table_file.write("$\\text{N(jet)} \geq 1$ & 1.56 $< |\eta| < $ 2.0 &  \\text{for} $E_{T} < $" + " %3.1f:" % par_1jet_end1[0] +  " %1.3f" % par_1jet_end1[1] + par2_string + " %1.5f" % abs(par_1jet_end1[2]) + " $\\times E_{T} \\text{, else }$ %1.4f" % ( par_1jet_end1[1] + ( par_1jet_end1[2] * par_1jet_end1[0] )) +  "\\\ \n")
#
# if par_1jet_end2[1] < 0: par2_string = " - "
# else                   : par2_string = " + "
#
# table_file.write("& 2.0 $< |\eta| < $ 2.5 &  \\text{for} $E_{T} < $" + " %3.1f:" % par_1jet_end2[0] +  " %1.3f" % par_1jet_end2[1] + par2_string + " %1.5f" % abs(par_1jet_end2[2]) + " $\\times E_{T} \\text{, else }$ %1.4f" % ( par_1jet_end2[1] + ( par_1jet_end2[2] * par_1jet_end2[0] )) +  "\\\ \n")
#
# table_file.write("  \hline \n")
#
# if par_2jet_bar[1] < 0: par2_string = " - "
# else                  : par2_string = " + "
#
# table_file.write("& $|\eta| < $ 1.442 &  \\text{for} $E_{T} < $" + " %3.1f:" % par_2jet_bar[0] + " %1.3f" % par_2jet_bar[1] + par2_string + " %1.5f" % abs(par_2jet_bar[2]) + " $\\times E_{T} \\text{, else }$ %1.4f" % ( par_2jet_bar[1] + ( par_2jet_bar[2] * par_2jet_bar[0] )) +  "\\\ \n")
#
# if par_2jet_end1[1] < 0: par2_string = " - "
# else                   : par2_string = " + "
#
# table_file.write("$\\text{N(jet)} \geq 2$ & 1.56 $< |\eta| < $ 2.0 &  \\text{for} $E_{T} < $" + " %3.1f:" % par_2jet_end1[0] +  " %1.3f" % par_2jet_end1[1] + par2_string + " %1.5f" % abs(par_2jet_end1[2]) + " $\\times E_{T} \\text{, else }$ %1.4f" % ( par_2jet_end1[1] + ( par_2jet_end1[2] * par_2jet_end1[0] )) +  "\\\ \n")
#
# if par_2jet_end2[1] < 0: par2_string = " - "
# else                   : par2_string = " + "
#
# table_file.write("& 2.0 $< |\eta| < $ 2.5 &  \\text{for} $E_{T} < $" + " %3.1f:" % par_2jet_end2[0] +  " %1.3f" % par_2jet_end2[1] + par2_string + " %1.5f" % abs(par_2jet_end2[2]) + " $\\times E_{T} \\text{, else }$ %1.4f" % ( par_2jet_end2[1] + ( par_2jet_end2[2] * par_2jet_end2[0] )) +  "\\\ \n")
#
# table_file.write("  \hline \n")
#
# if par_3jet_bar[1] < 0: par2_string = " - "
# else                  : par2_string = " + "
#
# table_file.write("& $|\eta| < $ 1.442 &  \\text{for} $E_{T} < $" + " %3.1f:" % par_3jet_bar[0] + " %1.3f" % par_3jet_bar[1] + par2_string + " %1.5f" % abs(par_3jet_bar[2]) + " $\\times E_{T} \\text{, else }$ %1.4f" % ( par_3jet_bar[1] + ( par_3jet_bar[2] * par_3jet_bar[0] )) +  "\\\ \n")
#
# if par_3jet_end1[1] < 0: par2_string = " - "
# else                   : par2_string = " + "
#
# table_file.write("$\\text{N(jet)} \geq 3$ & 1.56 $< |\eta| < $ 2.0 &  \\text{for} $E_{T} < $" + " %3.1f:" % par_3jet_end1[0] +  " %1.3f" % par_3jet_end1[1] + par2_string + " %1.5f" % abs(par_3jet_end1[2]) + " $\\times E_{T} \\text{, else }$ %1.4f" % ( par_3jet_end1[1] + ( par_3jet_end1[2] * par_3jet_end1[0] )) +  "\\\ \n")
#
# if par_3jet_end2[1] < 0: par2_string = " - "
# else                   : par2_string = " + "
#
# table_file.write("& 2.0 $< |\eta| < $ 2.5 &  \\text{for} $E_{T} < $" + " %3.1f:" % par_3jet_end2[0] +  " %1.3f" % par_3jet_end2[1] + par2_string + " %1.5f" % abs(par_3jet_end2[2]) + " $\\times E_{T} \\text{, else }$ %1.4f" % ( par_3jet_end2[1] + ( par_3jet_end2[2] * par_3jet_end2[0] )) +  "\\\ \n")
#
# table_file.write("  \hline \n")
# table_file.write("\end{tabular}\n")
# table_file.write("\end{table}\n")
# table_file.write("\end{document}\n")
#
# table_file.close()
#
# print "...Done"
#
# print "Making latex pdf"
#
# os.system ("pdflatex " + table_file_name +  " > /dev/null ")
#
# print "...Done"
#
##-----------------------------------------------------------------
## Combine the PDFs
##-----------------------------------------------------------------
#
# command = "gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=fakeRate.pdf " + table_file_pdf + " " + plot_file_pdf
#
# os.system ( command )
#
# os.system ( "rm " + table_file_pdf )
# os.system ( "rm " + plot_file_pdf )
# os.system ( "rm *.aux *.log")
#
##-----------------------------------------------------------------
## Tell the user what to do:
##-----------------------------------------------------------------
#
# bar_high = par_2jet_bar[1] + ( par_2jet_bar[2] * par_2jet_bar[0] )
# e_bar_high1 = par_2jet_bar[2] * par_2jet_bar[0] * math.sqrt ( ( ( err_2jet_bar[2] / par_2jet_bar[2] ) * ( err_2jet_bar[2] / par_2jet_bar[2] ) ) +
#                                                              ( ( err_2jet_bar[0] / par_2jet_bar[0] ) * ( err_2jet_bar[0] / par_2jet_bar[0] ) ) )
# e_bar_high = math.sqrt ( ( err_2jet_bar[1] * err_2jet_bar[1] ) + ( e_bar_high1 * e_bar_high1 ) )
#
# end1_high = par_2jet_end1[1] + ( par_2jet_end1[2] * par_2jet_end1[0] )
# e_end1_high1 = par_2jet_end1[2] * par_2jet_end1[0] * math.sqrt ( ( ( err_2jet_end1[2] / par_2jet_end1[2] ) * ( err_2jet_end1[2] / par_2jet_end1[2] ) ) +
#                                                              ( ( err_2jet_end1[0] / par_2jet_end1[0] ) * ( err_2jet_end1[0] / par_2jet_end1[0] ) ) )
# e_end1_high = math.sqrt ( ( err_2jet_end1[1] * err_2jet_end1[1] ) + ( e_end1_high1 * e_end1_high1 ) )
#
# end2_high = par_2jet_end2[1] + ( par_2jet_end2[2] * par_2jet_end2[0] )
# e_end2_high1 = par_2jet_end2[2] * par_2jet_end2[0] * math.sqrt ( ( ( err_2jet_end2[2] / par_2jet_end2[2] ) * ( err_2jet_end2[2] / par_2jet_end2[2] ) ) +
#                                                              ( ( err_2jet_end2[0] / par_2jet_end2[0] ) * ( err_2jet_end2[0] / par_2jet_end2[0] ) ) )
# e_end2_high = math.sqrt ( ( err_2jet_end2[1] * err_2jet_end2[1] ) + ( e_end2_high1 * e_end2_high1 ) )
#
# print "\n\nPut this in your cut file (for analysis): \n\n"
# print "fakeRate_bar  \t"  + "%.3e" % par_2jet_bar [1] + "\t" + "%.3e" % par_2jet_bar [2] + "\t" + "%.3e" % bar_high    + "\t" +"%.3f" % par_2jet_bar [0] +  "\t\t-1"
# print "fakeRate_end1 \t"  + "%.3e" % par_2jet_end1[1] + "\t" + "%.3e" % par_2jet_end1[2] + "\t" + "%.3e" % end1_high   + "\t" +"%.3f" % par_2jet_end1[0] +  "\t\t-1"
# print "fakeRate_end2 \t"  + "%.3e" % par_2jet_end2[1] + "\t" + "%.3e" % par_2jet_end2[2] + "\t" + "%.3e" % end2_high   + "\t" +"%.3f" % par_2jet_end2[0] +  "\t\t-1"
# print "\n\n"
# print "eFakeRate_bar  \t" + "%.3e" % err_2jet_bar [1] + "\t" + "%.3e" % err_2jet_bar [2] + "\t" + "%.3e" % e_bar_high  + "\t" +"%.3f" % err_2jet_bar [0] +  "\t\t-1"
# print "eFakeRate_end1 \t" + "%.3e" % err_2jet_end1[1] + "\t" + "%.3e" % err_2jet_end1[2] + "\t" + "%.3e" % e_end1_high + "\t" +"%.3f" % err_2jet_end1[0] +  "\t\t-1"
# print "eFakeRate_end2 \t" + "%.3e" % err_2jet_end2[1] + "\t" + "%.3e" % err_2jet_end2[2] + "\t" + "%.3e" % e_end2_high + "\t" +"%.3f" % err_2jet_end2[0] +  "\t\t-1"
# print "\n\n"
#
#
# print "\n\nPut this in your cut file (for closure test): \n\n"
# print "fakeRate_bar  \t"  + "%.3e" % par_1jet_bar [1] + "\t" + "%.3e" % par_1jet_bar [2] + "\t" + "%.3e" % bar_high    + "\t" +"%.3f" % par_1jet_bar [0] +  "\t\t-1"
# print "fakeRate_end1 \t"  + "%.3e" % par_1jet_end1[1] + "\t" + "%.3e" % par_1jet_end1[2] + "\t" + "%.3e" % end1_high   + "\t" +"%.3f" % par_1jet_end1[0] +  "\t\t-1"
# print "fakeRate_end2 \t"  + "%.3e" % par_1jet_end2[1] + "\t" + "%.3e" % par_1jet_end2[2] + "\t" + "%.3e" % end2_high   + "\t" +"%.3f" % par_1jet_end2[0] +  "\t\t-1"
# print "\n\n"
# print "eFakeRate_bar  \t" + "%.3e" % err_1jet_bar [1] + "\t" + "%.3e" % err_1jet_bar [2] + "\t" + "%.3e" % e_bar_high  + "\t" +"%.3f" % err_1jet_bar [0] +  "\t\t-1"
# print "eFakeRate_end1 \t" + "%.3e" % err_1jet_end1[1] + "\t" + "%.3e" % err_1jet_end1[2] + "\t" + "%.3e" % e_end1_high + "\t" +"%.3f" % err_1jet_end1[0] +  "\t\t-1"
# print "eFakeRate_end2 \t" + "%.3e" % err_1jet_end2[1] + "\t" + "%.3e" % err_1jet_end2[2] + "\t" + "%.3e" % e_end2_high + "\t" +"%.3f" % err_1jet_end2[0] +  "\t\t-1"
# print "\n\n"
#
