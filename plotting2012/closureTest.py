
import sys,os, math
import subprocess as sp
from ROOT import *


ROOT.gStyle.SetOptStat(0) 

do_plots = True

cej_data_dat_file_name = os.environ["LQDATA"] + "//eejj_analysis/QCDClosureTest_cej_DataAndMC/output_cutTable_lq_QCD_FakeRateClosureTest_cej_DataAndMC/analysisClass_lq_QCD_FakeRateClosureTest_tables.dat"

ccj_data_dat_file_name = os.environ["LQDATA"] + "/eejj_analysis/QCDClosureTest_ccj_DataOnly/output_cutTable_lq_QCD_FakeRateClosureTest_ccj_DataOnly/analysisClass_lq_QCD_FakeRateClosureTest_tables.dat"

ccjp_data_dat_file_name = os.environ["LQDATA"] + "/eejj_analysis/QCDClosureTest_ccj_DataOnly_Plus1Sigma/output_cutTable_lq_QCD_FakeRateClosureTest_ccj_DataOnly/analysisClass_lq_QCD_FakeRateClosureTest_tables.dat"

ccjm_data_dat_file_name = os.environ["LQDATA"] + "/eejj_analysis/QCDClosureTest_ccj_DataOnly_Minus1Sigma/output_cutTable_lq_QCD_FakeRateClosureTest_ccj_DataOnly/analysisClass_lq_QCD_FakeRateClosureTest_tables.dat"

ccj_data_file_name = os.environ["LQDATA"] + "/eejj_analysis/QCDClosureTest_ccj_DataOnly/output_cutTable_lq_QCD_FakeRateClosureTest_ccj_DataOnly/analysisClass_lq_QCD_FakeRateClosureTest_plots.root"

cej_mc_file_name =  os.environ["LQDATA"] + "//eejj_analysis/QCDClosureTest_cej_DataAndMC/output_cutTable_lq_QCD_FakeRateClosureTest_cej_DataAndMC/analysisClass_lq_QCD_FakeRateClosureTest_plots.root"

print "c-e-j dat file:" 
os.system ( "ls -l " + cej_data_dat_file_name )
print "\n"

print "c-c-j dat file:" 
os.system ( "ls -l " + ccj_data_dat_file_name )
print "\n"

print "c-c-j dat file (+1 sigma):"
os.system ("ls -l " + ccjp_data_dat_file_name )
print "\n"

print "c-c-j dat file (-1 sigma):"
os.system ("ls -l " + ccjm_data_dat_file_name )
print "\n"

print "c-c-j data plots:"
os.system ("ls -l " + ccj_data_file_name )
print "\n"

print "c-e-j data and mc plots:"
os.system ("ls -l " + cej_mc_file_name )
print "\n"

def getFRHist ( ccj_data_file, variable ) :
    hist = ccj_data_file.Get("histo1D__DATA__" + variable)
    return hist

def getMCHist ( cej_mc_file, variable ) :
    data_hist = TH1F()
    data_hist = cej_mc_file.Get("histo1D__DATA__"   + variable)
    mc_hist   = cej_mc_file.Get("histo1D__ALLBKG__" + variable)
    hist = data_hist.Clone() 
    hist.Add ( mc_hist, -1 ) 
    return hist

def makePlots ( psfile, canvas, ccj_data_file, cej_mc_file , variable, xtitle, xmin, xmax, ymin, ymax, rebin, log ) :
    
    canvas.Clear()
    canvas.cd()
    if log == True: canvas.SetLogy()
    else          : canvas.SetLogy(0)

    ccj_data_file = TFile ( ccj_data_file_name ) 
    cej_mc_file = TFile ( cej_mc_file_name ) 

    fr_hist = getFRHist ( ccj_data_file, variable ) 
    mc_hist = getMCHist ( cej_mc_file, variable ) 

    fr_hist.GetYaxis().SetTitle("Events") 
    fr_hist.GetXaxis().SetTitle(xtitle)
    
    fr_hist.SetMaximum(ymax)
    mc_hist.SetMaximum(ymax)
    fr_hist.SetMinimum(ymin)
    mc_hist.SetMinimum(ymin)

    fr_hist_int = float ( fr_hist.Integral ( 0 , fr_hist.GetNbinsX() + 1 ))
    mc_hist_int = float ( mc_hist.Integral ( 0 , mc_hist.GetNbinsX() + 1 ))

    print "fr hist int = " + str ( fr_hist_int ) 
    print "mc hist int = " + str ( mc_hist_int ) 
    print "fr / mc     = " + str ( fr_hist_int / mc_hist_int ) 
    
    fr_hist.GetXaxis().SetRangeUser ( xmin, xmax ) 
    mc_hist.GetXaxis().SetRangeUser ( xmin, xmax ) 

    fr_hist.Rebin ( rebin ) 
    mc_hist.Rebin ( rebin ) 

    hsize=0.60
    vsize=0.15
    xstart=0.30
    ystart=0.75
    
    fr_hist.SetLineColor(kBlue)
    fr_hist.SetMarkerColor(kBlue)
    fr_hist.SetMarkerStyle(20)
    mc_hist.SetLineColor(kRed)
    mc_hist.SetMarkerColor(kRed)
    mc_hist.SetMarkerStyle(23)

    legend = TLegend(xstart, ystart, xstart+hsize, ystart+vsize)
    legend.SetFillColor(kWhite)
    legend.SetBorderSize(0)
    legend.SetShadowColor(10)
    legend.SetMargin(0.2)
    legend.SetTextFont(132)
    legend.AddEntry ( fr_hist,"QCD: Fake rate prediction"  , "pl" ) 
    legend.AddEntry ( mc_hist,"QCD: Data - non QCD MC truth", "pl" ) 

    fr_hist.Draw()
    mc_hist.Draw("SAME")
    legend.Draw()

    canvas.Print (psfile)
    canvas.SaveAs ( "closureTest_eps/" + variable + ".eps"  );
    canvas.SaveAs ( "closureTest_gif/" + variable + ".gif"  );
    
    ccj_data_file.Close()
    cej_mc_file.Close()


def getDict ( file_name, data_name ) :
    file = open ( file_name, "r" ) 
    file_info = file.read() 
    data_info = file_info.split(data_name)[1].split("\n\n")[0]
    
    list = []
    dict = {} 

    data_info_lines = data_info.split("\n")
    for line in data_info_lines : 
        split_line = line.split()
        print line, len ( split_line ) 
        if len ( split_line ) < 11 : continue
        if split_line[0] == "variableName" : continue
        
        variable = split_line [0]
        min1     = split_line [1]
        max1     = split_line [2]
        npass    = split_line [5]
        eNpass   = split_line [6] 
        list.append ( variable ) 
        dict [ variable ] = [npass, eNpass, min1, max1 ]
    return list, dict 
        

variable_list , cej_allbkg_dict = getDict ( cej_data_dat_file_name , "ALLBKG" ) 
print "1",  cej_allbkg_dict["sT_eej_200"][2], cej_allbkg_dict["sT_eej_200"][3]

variable_list , cej_data_dict   = getDict ( cej_data_dat_file_name , "DATA"   ) 
print "2",  cej_data_dict["sT_eej_200"][2], cej_data_dict["sT_eej_200"][3]

variable_list , ccj_data_dict   = getDict ( ccj_data_dat_file_name  , "DATA"   ) 
print "3",  ccj_data_dict["sT_eej_200"][2], ccj_data_dict["sT_eej_200"][3]

variable_list , ccjp_data_dict  = getDict ( ccjp_data_dat_file_name, "DATA"   ) 
print "4",  ccjp_data_dict["sT_eej_200"][2], ccjp_data_dict["sT_eej_200"][3]

variable_list , ccjm_data_dict  = getDict ( ccjm_data_dat_file_name, "DATA"   ) 
print "5",  ccjm_data_dict["sT_eej_200"][2], ccjm_data_dict["sT_eej_200"][3]

tab_length = 7

latex_file_name = "closureTest_tmp.tex"
latex_file_pdf = latex_file_name.replace (".tex", ".pdf")

latex_file = open (latex_file_name,"w")

latex_file.write("\documentclass{article}\n")
latex_file.write("\usepackage{multirow} \n")
latex_file.write("\usepackage[landscape, top=1cm, bottom=1cm, left=1cm, right=1cm]{geometry}\n")
latex_file.write("\pagestyle{empty} \n")
latex_file.write("\\begin{document}\n")
latex_file.write("\\begin{table}\n")
latex_file.write("\\begin{tabular}{| l | c | c | c |} \n")
latex_file.write("  \hline \n")
latex_file.write("Cut & Prediction from fake rate & Prediction from Data - MC & Ratio \\\ \n" )
latex_file.write("  \hline \n")
for variable in variable_list[6:] : 
    
    predicted        = ccj_data_dict   [variable][0]

    predicted_p      = ccjp_data_dict  [variable][0]
    predicted_m      = ccjm_data_dict  [variable][0]

    ePredicted = max ( abs ( float (predicted) - float(predicted_m) ) , abs ( float (predicted) - float(predicted_p) ) ) 
                                       
    observed_data    = cej_data_dict   [variable][0]
    eObserved_data   = cej_data_dict   [variable][1]
                                       
    observed_allbkg  = cej_allbkg_dict [variable][0]
    eObserved_allbkg = cej_allbkg_dict [variable][1]

    cut_min1 = cej_data_dict[variable][2]
    cut_max1 = cej_data_dict[variable][3]

    line = ""
    
    variable_tex = variable.replace ("_", "\_")

    if variable_tex == "nocut" or variable_tex == "skim":
        line = variable_tex
    elif cut_min1.strip() == "-inf" and cut_max1.strip() == "+inf":
        continue
    elif cut_min1.strip() != "-inf" and cut_max1.strip() == "+inf":
        line = variable_tex + " $>$ " + cut_min1 
    elif cut_min1.strip() == "-inf" and cut_max1.strip() != "+inf":
        line = variable_tex + " $\leq$ " + cut_max1 
    elif cut_min1.strip() != "-inf" and cut_max1.strip() != "+inf":
        line = cut_min1 + " $<$ " + variable_tex + " $\leq$ " + cut_max1 

    observed         = float(observed_data) - (float(observed_allbkg) / 1000.) 
    eObserved        = math.sqrt ( ( (float(eObserved_data  )      ) * ( float(eObserved_data  )        ) ) + 
                                   ( (float(eObserved_allbkg)/1000.) * ( float(eObserved_allbkg) / 1000.) ) )
    
    ratio  = float ( predicted ) / float ( observed ) 
    eRatio = ratio * math.sqrt ( ( ( float(ePredicted) / float(predicted ) ) * ( float(ePredicted) / float(predicted ) ) ) + 
                                 ( ( float(eObserved ) / float(observed  ) ) * ( float(eObserved ) / float(observed  ) ) ) )

    space = "\t\t\t"
    if len (variable.strip() ) < tab_length:     space = space + "\t"
    if len (variable.strip() ) > tab_length * 2: space = space[:-1]
    if len (variable.strip() ) > tab_length * 3: space = space[:-1]

    line = line + " & " + "%.3e" % float(predicted) +" $\pm$ " + "%.3e" % ePredicted + " & " + "%.3e" % observed +" $\pm$ " + "%.3e" % eObserved + " & " +  "%.3e" % ratio +" $\pm$ " + "%.3e" % eRatio + " \\\ \n" 
    
    latex_file.write( line ) 

latex_file.write("  \hline \n")
latex_file.write("\end{tabular}\n")
latex_file.write("\end{table}\n")
latex_file.write("\end{document}\n")

latex_file.close()

os.system ( "pdflatex " + latex_file_name ) 

if do_plots:

    plot_file_name = "plots.ps"
    plot_file_pdf = plot_file_name.replace(".ps",".pdf")
    
    canvas = TCanvas("canvas","")
    
    canvas.Print (plot_file_name + "[")
    
    makePlots ( plot_file_name, canvas, ccj_data_file_name, cej_mc_file_name , "sT_PAS" , "s_{T} (e_{1} p_{T} + e_{2} p_{T} + j_{1} p_{T}) [GeV] (preselection)"  , 0   , 2000 , 0.005, 10000000, 10, True ) 
    makePlots ( plot_file_name, canvas, ccj_data_file_name, cej_mc_file_name , "Mee_PAS", "M(e_{1}, e_{2}) [GeV] (preselection)"                                  , 0   , 1000 , 0.005, 1000000 , 10, True ) 
    makePlots ( plot_file_name, canvas, ccj_data_file_name, cej_mc_file_name , "Pt1stEle_PAS" , "1st electron p_{T} [GeV] (preselection)"                         , 0   , 1000 , 0.005, 1000000 , 5, True ) 
    makePlots ( plot_file_name, canvas, ccj_data_file_name, cej_mc_file_name , "Eta1stEle_PAS", "1st electron #eta (preselection)"                                , -3.0, 3.0  , 0.005, 10000000, 5, True ) 
    makePlots ( plot_file_name, canvas, ccj_data_file_name, cej_mc_file_name , "Eta1stEle_PAS", "1st electron #eta (preselection)"                                , -3.0, 3.0  , 0.005, 10000   , 5, False ) 
    makePlots ( plot_file_name, canvas, ccj_data_file_name, cej_mc_file_name , "Pt2ndEle_PAS" , "2nd electron p_{T} [GeV] (preselection)"                         , 0   , 1000 , 0.005, 200000  , 5, True ) 
    makePlots ( plot_file_name, canvas, ccj_data_file_name, cej_mc_file_name , "Eta2ndEle_PAS", "2nd electron #eta (preselection)"                                , -3.0, 3.0  , 0.005, 10000000, 5, True ) 
    makePlots ( plot_file_name, canvas, ccj_data_file_name, cej_mc_file_name , "Eta2ndEle_PAS", "2nd electron #eta (preselection)"                                , -3.0, 3.0  , 0.005, 20000  , 5, False ) 
    makePlots ( plot_file_name, canvas, ccj_data_file_name, cej_mc_file_name , "Pt1stJet_PAS" , "1st jet p_{T} [GeV] (preselection)"                              , 0   , 1000 , 0.005, 1000000 , 5, True ) 
    makePlots ( plot_file_name, canvas, ccj_data_file_name, cej_mc_file_name , "Eta1stJet_PAS", "1st jet #eta (preselection)"                                     , -3.0, 3.0  , 0.005, 10000000, 5, True ) 
    makePlots ( plot_file_name, canvas, ccj_data_file_name, cej_mc_file_name , "Me1j1_PAS"    , "M(e_{1}, j_{1}) [GeV] (preselection)"                            , 0   , 1000 , 0.005, 1000000  , 5, True ) 
    makePlots ( plot_file_name, canvas, ccj_data_file_name, cej_mc_file_name , "Me2j1_PAS"    , "M(e_{2}, j_{1}) [GeV] (preselection)"                            , 0   , 1000 , 0.005, 1000000 , 5, True ) 
    
    canvas.Print (plot_file_name + "]")

    os.system ( "ps2pdf " + plot_file_name ) 
    os.system ( "rm " + plot_file_name ) 


#-----------------------------------------------------------------
# Combine the PDFs
#-----------------------------------------------------------------

    command = "gs -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=closureTest.pdf " + latex_file_pdf + " " + plot_file_pdf

    os.system ( command ) 

if do_plots: os.system ( "rm " + latex_file_pdf ) 
os.system ( "rm *.aux *.tex *.log") 

if do_plots: os.system ( "rm " + plot_file_pdf ) 
