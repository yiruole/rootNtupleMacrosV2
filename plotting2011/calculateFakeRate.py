import ROOT
import os, sys
import subprocess as sp

ROOT.gErrorIgnoreLevel = 1

def getNumeratorHist   (root_file, variable ) :
    data_hist = ROOT.TH1F()
    data_hist = root_file.Get("histo1D__DATA__Pass_"   + variable)
    mc_hist   = root_file.Get("histo1D__ALLBKG__Pass_" + variable)
    hist = data_hist.Clone() 
    hist.Add ( mc_hist, -1 ) 
    return hist
    
def getDenominatorHist (root_file, variable ) :
    data_hist = ROOT.TH1F()
    data_hist = root_file.Get("histo1D__DATA__Total_" + variable)
    mc_hist   = root_file.Get("histo1D__ALLBKG__Total_" + variable)
    hist = data_hist.Clone() 
    hist.Add ( mc_hist, -1 ) 
    return hist

def makePlots ( psfile, canvas, root_file, variable, xtitle, xmin, xmax, ymin, ymax ) :

    canvas.Clear()

    n_hist = getNumeratorHist   (root_file, variable ) 
    d_hist = getDenominatorHist (root_file, variable ) 
    r_hist = n_hist.Clone()
    r_hist.Divide ( d_hist ) 
    
    d_hist.GetXaxis().SetRangeUser ( xmin, xmax ) 
    n_hist.GetXaxis().SetRangeUser ( xmin, xmax ) 
    r_hist.GetXaxis().SetRangeUser ( xmin, xmax ) 

    d_hist.SetMaximum ( d_hist.GetMaximum() * 10000 ) 
    n_hist.SetMaximum ( d_hist.GetMaximum() * 10000 ) 
    d_hist.SetMinimum ( 0.1 ) 
    n_hist.SetMinimum ( 0.1 ) 

    r_hist.SetMaximum ( ymax ) 
    r_hist.SetMinimum ( ymin ) 

    d_hist.GetXaxis().SetTitle ( xtitle ) 
    d_hist.GetYaxis().SetTitleOffset ( 1.5 ) 
    d_hist.GetYaxis().SetTitle ( "Data - non-QCD MC" )

    d_hist.SetLineColor  (ROOT.kRed)
    d_hist.SetMarkerColor(ROOT.kRed)
    n_hist.SetLineColor  (ROOT.kBlue)
    n_hist.SetMarkerColor(ROOT.kBlue)

    hsize=0.60
    vsize=0.15
    xstart=0.30
    ystart=0.75

    legend = ROOT.TLegend(xstart, ystart, xstart+hsize, ystart+vsize)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetBorderSize(0)
    legend.SetShadowColor(10)
    legend.SetMargin(0.2)
    legend.SetTextFont(132)
    legend.AddEntry ( d_hist,"Data - non-QCD MC", "pl" ) 
    legend.AddEntry ( n_hist,"Data - non-QCD MC, passing HEEP v3.1", "pl" ) 

    canvas.Divide ( 2, 1)
    pad1 = canvas.GetPad ( 1 ) 
    pad2 = canvas.GetPad ( 2 ) 
    pad1.cd() 
    pad1.SetLogy()
    d_hist.Draw()
    n_hist.Draw("SAME")
    legend.Draw()

    pad2.cd()
    r_hist.GetXaxis().SetTitle ( xtitle ) 
    r_hist.GetYaxis().SetTitle ( "Fake rate" ) 
    r_hist.GetYaxis().SetTitleOffset ( 1.5 ) 
    r_hist.Draw()

    f1 = ROOT.TF1 ("f1","pol1[0]", xmin,  100 )
    f2 = ROOT.TF1 ("f2","pol0[0]", 100  , xmax)
    r_hist.Fit ("f1","QR")
    pad2.Update()
    stats1        = r_hist.GetListOfFunctions().FindObject("stats")
    cloned_stats1 = stats1.Clone();
    x1 = cloned_stats1.GetX1NDC();
    x2 = cloned_stats1.GetX2NDC();
    length = x2 - x1
    cloned_stats1.SetX1NDC ( x1 - ( 1.1 * length ) )
    cloned_stats1.SetX2NDC ( cloned_stats1.GetX1NDC() + length )
    cloned_stats1.SetBorderSize(2)

    r_hist.Fit ("f2","QR")
    f1.Draw("SAME")
    cloned_stats1.Draw("SAME")

    ROOT.gErrorIgnoreLevel = 1
    canvas.Print (psfile)
    canvas.SaveAs ( "fakeRate_eps/" + variable + ".eps"  );
    canvas.SaveAs ( "fakeRate_gif/" + variable + ".gif"  );

    parameters = [] 
    par_errors = [] 

    parameters.append ( f1.GetParameter ( 0 ) )
    parameters.append ( f1.GetParameter ( 1 ) )
    parameters.append ( f2.GetParameter ( 0 ) )

    par_errors.append ( f1.GetParError  ( 0 ) )
    par_errors.append ( f1.GetParError  ( 1 ) )
    par_errors.append ( f2.GetParError  ( 0 ) )

    return parameters, par_errors
    
#-----------------------------------------------------------------
# Get the ROOT file
#-----------------------------------------------------------------

if len ( sys.argv ) != 2 : 
    print "Usage: python calculateFakeRate.py <root file>"
    sys.exit() 

if not os.path.isfile ( sys.argv[1] ) : 
    print "Usage: python calculateFakeRate.py <root file>"
    print sys.argv[1] + " is not a file " 
    sys.exit()

root_file_path = sys.argv[1]

root_file = ROOT.TFile ( root_file_path ) 

#-----------------------------------------------------------------
# Make the plots
#-----------------------------------------------------------------

plot_file_name = "plots.ps"
plot_file_pdf = plot_file_name.replace(".ps",".pdf")

canvas = ROOT.TCanvas("canvas","",1200,500)

ROOT.gErrorIgnoreLevel = 1
canvas.Print (plot_file_name + "[")

par_inc_bar , err_inc_bar  = makePlots ( plot_file_name, canvas, root_file, "Bar_Pt1stEle_PAS" , "1st Electron p_{T} (Barrel) [GeV]"              , 35, 250, 0.0, 0.15 )
par_inc_end1, err_inc_end1 = makePlots ( plot_file_name, canvas, root_file, "End1_Pt1stEle_PAS", "1st Electron p_{T} (Endcap, |#eta| < 2.0) [GeV]", 35, 250, 0.0, 0.15 )
par_inc_end2, err_inc_end2 = makePlots ( plot_file_name, canvas, root_file, "End2_Pt1stEle_PAS", "1st Electron p_{T} (Endcap, |#eta| > 2.0) [GeV]", 35, 250, 0.0, 0.15 )

par_1jet_bar , err_1jet_bar  = makePlots ( plot_file_name, canvas, root_file, "Bar_1Jet_Pt1stEle_PAS" , "1st Electron p_{T} (Barrel, 1 Jet) [GeV]"                   , 35, 250, 0.0, 0.15 )
par_1jet_end1, err_1jet_end1 = makePlots ( plot_file_name, canvas, root_file, "End1_1Jet_Pt1stEle_PAS", "1st Electron p_{T} (Endcap, |#eta| < 2.0, #geq 1 Jet) [GeV]", 35, 250, 0.0, 0.15 )
par_1jet_end2, err_1jet_end2 = makePlots ( plot_file_name, canvas, root_file, "End2_1Jet_Pt1stEle_PAS", "1st Electron p_{T} (Endcap, |#eta| > 2.0, #geq 1 Jet) [GeV]", 35, 250, 0.0, 0.15 )

par_2jet_bar , err_2jet_bar  = makePlots ( plot_file_name, canvas, root_file, "Bar_2Jet_Pt1stEle_PAS" , "1st Electron p_{T} (Barrel, 2 Jets) [GeV]"                    , 35, 250, 0.0, 0.15 )
par_2jet_end1, err_2jet_end1 = makePlots ( plot_file_name, canvas, root_file, "End1_2Jet_Pt1stEle_PAS", "1st Electron p_{T} (Endcap, |#eta| < 2.0, #geq 2 Jets ) [GeV]", 35, 250, 0.0, 0.15 )
par_2jet_end2, err_2jet_end2 = makePlots ( plot_file_name, canvas, root_file, "End2_2Jet_Pt1stEle_PAS", "1st Electron p_{T} (Endcap, |#eta| > 2.0, #geq 2 Jets ) [GeV]", 35, 250, 0.0, 0.15 )

par_3jet_bar , err_3jet_bar  = makePlots ( plot_file_name, canvas, root_file, "Bar_3Jet_Pt1stEle_PAS" , "1st Electron p_{T} (Barrel, 3 Jets) [GeV]"                    , 35, 250, 0.0, 0.15 )
par_3jet_end1, err_3jet_end1 = makePlots ( plot_file_name, canvas, root_file, "End1_3Jet_Pt1stEle_PAS", "1st Electron p_{T} (Endcap, |#eta| < 2.0, #geq 3 Jets ) [GeV]", 35, 250, 0.0, 0.15 )
par_3jet_end2, err_3jet_end2 = makePlots ( plot_file_name, canvas, root_file, "End2_3Jet_Pt1stEle_PAS", "1st Electron p_{T} (Endcap, |#eta| > 2.0, #geq 3 Jets ) [GeV]", 35, 250, 0.0, 0.15 )

ROOT.gErrorIgnoreLevel = 1
canvas.Print (plot_file_name + "]")

os.system ( "ps2pdf " + plot_file_name ) 
os.system ( "rm " + plot_file_name ) 

#-----------------------------------------------------------------
# Make the table
#-----------------------------------------------------------------



table_file_name = "table.tex"
table_file = open ( table_file_name, "w" )  
table_file_pdf = table_file_name.replace (".tex",".pdf")

table_file.write("\documentclass{article}\n")
table_file.write("\usepackage{multirow} \n")
table_file.write("\usepackage{amsmath} \n")
table_file.write("\usepackage[landscape, top=1cm, bottom=1cm, left=1cm, right=1cm]{geometry}\n")
table_file.write("\pagestyle{empty} \n")
table_file.write("\\begin{document}\n")
table_file.write("\[\n")
table_file.write("fit(x) =\n")
table_file.write("\\begin{cases}\n")
table_file.write("p_{0} + p_{1} x & \\text{if } x < 100 \\\ \n")
table_file.write("p_{2}       & \\text{if } x \geq 100 \n")
table_file.write("\end{cases} \n")
table_file.write("\\] \n")
table_file.write("\\begin{table}[h]\n")
table_file.write("\\begin{tabular}{l c c c}\n")
table_file.write(" Selection & $p_{0}$ & $p_{1}$ & $p_{2}$ \\\ \n" )
table_file.write("  \hline \n")
table_file.write("  \hline \n")
table_file.write("inclusive, barrel & "                   + "%.3e" % par_inc_bar  [0] + " $\pm$ " + "%.3e" % err_inc_bar  [0]+  " & " + "%.3e" % par_inc_bar  [1] + " $\pm$ " + "%.3e" % err_inc_bar  [1]+ " & " +"%.3e" % par_inc_bar  [2]  + " $\pm$ " + "%.3e" % err_inc_bar  [2]  + " \\\ \n")
table_file.write("inclusive, endcap ($|\eta| < 2.0$ ) & " + "%.3e" % par_inc_end1 [0] + " $\pm$ " + "%.3e" % err_inc_end1 [0]+  " & " + "%.3e" % par_inc_end1 [1] + " $\pm$ " + "%.3e" % err_inc_end1 [1]+ " & " +"%.3e" % par_inc_end1 [2]  + " $\pm$ " + "%.3e" % err_inc_end1 [2]  + " \\\ \n")
table_file.write("inclusive, endcap ($|\eta| > 2.0$ ) & " + "%.3e" % par_inc_end2 [0] + " $\pm$ " + "%.3e" % err_inc_end2 [0]+  " & " + "%.3e" % par_inc_end2 [1] + " $\pm$ " + "%.3e" % err_inc_end2 [1]+ " & " +"%.3e" % par_inc_end2 [2]  + " $\pm$ " + "%.3e" % err_inc_end2 [2]  + " \\\ \n")
table_file.write("  \hline \n") 
table_file.write("1 jet, barrel & "                       + "%.3e" % par_1jet_bar [0] + " $\pm$ " + "%.3e" % err_1jet_bar [0]+  " & " + "%.3e" % par_1jet_bar [1] + " $\pm$ " + "%.3e" % err_1jet_bar [1]+ " & " +"%.3e" % par_1jet_bar [2]  + " $\pm$ " + "%.3e" % err_1jet_bar [2]  + " \\\ \n")
table_file.write("1 jet, endcap ($|\eta| < 2.0$ ) & "     + "%.3e" % par_1jet_end1[0] + " $\pm$ " + "%.3e" % err_1jet_end1[0]+  " & " + "%.3e" % par_1jet_end1[1] + " $\pm$ " + "%.3e" % err_1jet_end1[1]+ " & " +"%.3e" % par_1jet_end1[2]  + " $\pm$ " + "%.3e" % err_1jet_end1[2]  + " \\\ \n")
table_file.write("1 jet, endcap ($|\eta| > 2.0$ ) & "     + "%.3e" % par_1jet_end2[0] + " $\pm$ " + "%.3e" % err_1jet_end2[0]+  " & " + "%.3e" % par_1jet_end2[1] + " $\pm$ " + "%.3e" % err_1jet_end2[1]+ " & " +"%.3e" % par_1jet_end2[2]  + " $\pm$ " + "%.3e" % err_1jet_end2[2]  + " \\\ \n")
table_file.write("  \hline \n") 
table_file.write("2 jets, barrel & "                      + "%.3e" % par_2jet_bar [0] + " $\pm$ " + "%.3e" % err_2jet_bar [0]+  " & " + "%.3e" % par_2jet_bar [1] + " $\pm$ " + "%.3e" % err_2jet_bar [1]+ " & " +"%.3e" % par_2jet_bar [2]  + " $\pm$ " + "%.3e" % err_2jet_bar [2]  + " \\\ \n")
table_file.write("2 jets, endcap ($|\eta| < 2.0$ ) & "    + "%.3e" % par_2jet_end1[0] + " $\pm$ " + "%.3e" % err_2jet_end1[0]+  " & " + "%.3e" % par_2jet_end1[1] + " $\pm$ " + "%.3e" % err_2jet_end1[1]+ " & " +"%.3e" % par_2jet_end1[2]  + " $\pm$ " + "%.3e" % err_2jet_end1[2]  + " \\\ \n")
table_file.write("2 jets, endcap ($|\eta| > 2.0$ ) & "    + "%.3e" % par_2jet_end2[0] + " $\pm$ " + "%.3e" % err_2jet_end2[0]+  " & " + "%.3e" % par_2jet_end2[1] + " $\pm$ " + "%.3e" % err_2jet_end2[1]+ " & " +"%.3e" % par_2jet_end2[2]  + " $\pm$ " + "%.3e" % err_2jet_end2[2]  + " \\\ \n")
table_file.write("  \hline \n")  
table_file.write("3 jets, barrel & "                      + "%.3e" % par_3jet_bar [0] + " $\pm$ " + "%.3e" % err_3jet_bar [0]+  " & " + "%.3e" % par_3jet_bar [1] + " $\pm$ " + "%.3e" % err_3jet_bar [1]+ " & " +"%.3e" % par_3jet_bar [2]  + " $\pm$ " + "%.3e" % err_3jet_bar [2]  + " \\\ \n")
table_file.write("3 jets, endcap ($|\eta| < 2.0$ ) & "    + "%.3e" % par_3jet_end1[0] + " $\pm$ " + "%.3e" % err_3jet_end1[0]+  " & " + "%.3e" % par_3jet_end1[1] + " $\pm$ " + "%.3e" % err_3jet_end1[1]+ " & " +"%.3e" % par_3jet_end1[2]  + " $\pm$ " + "%.3e" % err_3jet_end1[2]  + " \\\ \n")
table_file.write("3 jets, endcap ($|\eta| > 2.0$ ) & "    + "%.3e" % par_3jet_end2[0] + " $\pm$ " + "%.3e" % err_3jet_end2[0]+  " & " + "%.3e" % par_3jet_end2[1] + " $\pm$ " + "%.3e" % err_3jet_end2[1]+ " & " +"%.3e" % par_3jet_end2[2]  + " $\pm$ " + "%.3e" % err_3jet_end2[2]  + " \\\ \n")
table_file.write("  \hline \n")
table_file.write("\end{tabular}\n")
table_file.write("\end{table}\n")
table_file.write("\end{document}\n")

table_file.close()

os.system ("pdflatex " + table_file_name + " > /dev/null ") 

#-----------------------------------------------------------------
# Combine the PDFs
#-----------------------------------------------------------------

command = "gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=fakeRate.pdf " + table_file_pdf + " " + plot_file_pdf

os.system ( command ) 

os.system ( "rm " + table_file_pdf ) 
os.system ( "rm " + plot_file_pdf ) 
os.system ( "rm *.aux *.log") 

#-----------------------------------------------------------------
# Tell the user what to do:
#-----------------------------------------------------------------

print "\n\nPut this in your cut file: \n\n"
print "fakeRate_bar  \t" + "%.3e" % par_2jet_bar [0] + "\t" + "%.3e" % par_2jet_bar [1] + "\t" + "%.3e" % par_2jet_bar [2] + "\t-\t-1"
print "fakeRate_end1 \t" + "%.3e" % par_2jet_end1[0] + "\t" + "%.3e" % par_2jet_end1[1] + "\t" + "%.3e" % par_2jet_end1[2] + "\t-\t-1"
print "fakeRate_end2 \t" + "%.3e" % par_2jet_end2[0] + "\t" + "%.3e" % par_2jet_end2[1] + "\t" + "%.3e" % par_2jet_end2[2] + "\t-\t-1"
print "\n\n"
print "eFakeRate_bar  \t" + "%.3e" % err_2jet_bar [0] + "\t" + "%.3e" % err_2jet_bar [1] + "\t" + "%.3e" % err_2jet_bar [2] + "\t-\t-1"
print "eFakeRate_end1 \t" + "%.3e" % err_2jet_end1[0] + "\t" + "%.3e" % err_2jet_end1[1] + "\t" + "%.3e" % err_2jet_end1[2] + "\t-\t-1"
print "eFakeRate_end2 \t" + "%.3e" % err_2jet_end2[0] + "\t" + "%.3e" % err_2jet_end2[1] + "\t" + "%.3e" % err_2jet_end2[2] + "\t-\t-1"
print "\n\n"
