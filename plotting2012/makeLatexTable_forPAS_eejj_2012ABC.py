import os, sys, math
import subprocess 
from ROOT import *

def getDataDrivenDict ( file_name ) :
    
    file = open ( file_name, "r" )
    file_contents = file.read()
    file.close() 
    
    raw_table = file_contents.split("DATA")[1].split("\n\n")[0].strip()
    column_names = raw_table.split("\n")[0].split()
    rows = raw_table.split("\n")
    
    d_cutName_cutData = {} 

    for row in rows[1:]:
        row = row.strip()
        row_fields = row.split() 
        cut_data = dict ( zip ( column_names, row_fields ) )
        cut_name = row_fields[0]
        d_cutName_cutData [ cut_name ] = cut_data 

    return d_cutName_cutData

def getScale (sample_name, lumi, n_events):
    xsection_file = open ("config/xsection_8TeV_2012ABC_ZScaled.txt","r")
    xsection_data = xsection_file.read() 
    xsection_file.close() 
    
    xsection = float(xsection_data.split(sample_name)[1].split("\n")[0].strip())
               
    scale = xsection * lumi / n_events
    print sample_name, xsection, scale


    return scale

TTBarScale = 0.49
QCDScale = 1.0

dat_file_name     = os.environ['LQDATA'] + "//eejj_analysis/eejj//scaled_2012ABC_output_cutTable_lq_eejj/analysisClass_lq_eejj_tables.dat"
qcd_dat_file_name = os.environ['LQDATA'] + "//eejj_analysis/eejj_qcd//output_cutTable_lq_eejj/analysisClass_lq_eejj_QCD_tables.dat"

lumi = float ( 12200. )

LQ_masses = [ 300,  350,  400, 450, 500, 550,  600,  650, 700,  750,  800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 ]

sampleNames = [
    "ZJet_Madgraph" ,
    "TTbar_Madgraph",
    "QCD_fromDATA"  ,
    "OTHERBKG"      ,
    "2012ABC"          
]

d_sampleName_sampleTitle = {
    "DIBOSON"       :"WW + WZ + ZZ",
    "TTbar_Madgraph":"$t\\bar{t}$" ,
    "TTbar_fromDATA":"$t\\bar{t}$" ,
    "QCD_fromDATA"  :"QCD"         ,
    "WJet_Madgraph"   :"W+Jets"      ,
    "ZJet_Madgraph"   :"Z+Jets"      ,
    "ZJet_tarball_Madgraph"   :"Z+Jets"      ,
    "PhotonJets"    :"$\gamma$ + jets",
    "SingleTop"     :"Single Top"  ,
    "2012ABC"          :"Data", 
    "OTHERBKG"      :"Other"
}

final_selection_cut_prefix = "M_e1e2_LQ"
pre_selection_cut_name     = "M_e1e2"

signalSampleNames = []
for mass in LQ_masses:
    signalSampleNames.append ( "LQ_M" + str ( mass ) )

d_sampleName_cutData = {} 

dat_file = open ( dat_file_name,"r") 
dat_file_contents = dat_file.read() 
dat_file.close()

for sample_name in sampleNames + signalSampleNames:
    if "fromDATA" in sample_name: continue
    raw_table = dat_file_contents.split(sample_name+"\n")[1].split("\n\n")[0].strip()
    column_names = raw_table.split("\n")[0].split()
    rows = raw_table.split("\n")
    
    d_cutName_cutData = {} 

    for row in rows[1:]:
        row = row.strip()
        row_fields = row.split() 
        cut_data = dict ( zip ( column_names, row_fields ) )
        cut_name = row_fields[0]
        d_cutName_cutData [ cut_name ] = cut_data 

    d_sampleName_cutData [ sample_name ] = d_cutName_cutData 

d_sampleName_cutData [ "QCD_fromDATA"  ] =  getDataDrivenDict ( qcd_dat_file_name ) 

latex_file_name = "table_finalSelection_eejj_2012ABC.tex"
latex_file = open ( latex_file_name, "w" )  

latex_file.write("\documentclass{article}\n")
latex_file.write("\usepackage{multirow} \n")
latex_file.write("\usepackage[landscape, top=1cm, bottom=1cm, left=1cm, right=1cm]{geometry}\n")
latex_file.write("\pagestyle{empty} \n")
latex_file.write("\\begin{document}\n")
latex_file.write("\\begin{table} \n")
latex_file.write("\\small \n")
line = "\\begin{tabular}{| l | c |"
for sampleName in sampleNames:
    line = line + " c |" 
line = line + " c |} \n"
latex_file.write ( line ) 
latex_file.write("  \hline \n")
line = "$M_{LQ}$ & LQ Signal & "
for sampleName in sampleNames:
    line = line + d_sampleName_sampleTitle[sampleName] + " & "
line = line + " Total BG \\\ \n"
latex_file.write ( line ) 
latex_file.write("  \hline \n")
latex_file.write("  \hline \n")
line = "Presel & - & "

total_bkg = 0
total_eBkg2 = 0

for sampleName in sampleNames:
    npass  = float ( d_sampleName_cutData [ sampleName ][ pre_selection_cut_name ][ "Npass"    ] )
    eNpass = float ( d_sampleName_cutData [ sampleName ][ pre_selection_cut_name ][ "errNpass" ] )

    if "fromdata" in sampleName.lower() :
        if "qcd" in sampleName.lower() : 
            npass  = npass  * QCDScale 
            eNpass = eNpass * QCDScale 
            print "qcd npass =", npass , "+/-", eNpass

        if "ttb" in sampleName.lower() : 
            npass  = npass  * TTBarScale 
            eNpass = eNpass * TTBarScale 
            
        total_bkg   = total_bkg + npass
        total_eBkg2 = total_eBkg2 + ( eNpass * eNpass ) 
        if ( eNpass >= 0.01 ):
            # line = line + " $ " + "%.2f" % npass + " \pm " +"%.2f" % eNpass + " $ &"
            line = line + " $ " + "%.3f" % npass + " \pm " +"%.3f" % eNpass + " $ &"
        else:
            # line = line + " $ " + "%.3f" % npass + " \pm " +"%.3f" % eNpass + " $ &"
            line = line + " $ " + "%.4f" % npass + " \pm " +"%.4f" % eNpass + " $ &"
    elif ( sampleName.lower() != "2012abc" ):
        npass = npass / 1000.
        eNpass = eNpass / 1000.
        total_bkg   = total_bkg + npass
        total_eBkg2 = total_eBkg2 + ( eNpass * eNpass ) 
        line = line + " $ " + "%.1f" % npass + " \pm " +"%.1f" % eNpass + " $ &"
    else:
        line = line + "%d" % npass + " &"
total_eBkg = math.sqrt ( total_eBkg2 ) 
line = line + " $ %.1f \pm " % total_bkg + "%.1f" % total_eBkg  + " $ \\\ \n"
latex_file.write ( line ) 
latex_file.write("  \hline \n")
    
for LQ_mass in LQ_masses:
    line = str ( LQ_mass ) + " & " 
    nLQPass  = float ( d_sampleName_cutData [ "LQ_M" + str ( LQ_mass ) ][ final_selection_cut_prefix + str ( LQ_mass )  ][ "Npass"    ] ) / 1000
    eNLQPass = float ( d_sampleName_cutData [ "LQ_M" + str ( LQ_mass ) ][ final_selection_cut_prefix + str ( LQ_mass )  ][ "errNpass" ] ) / 1000
    if ( eNLQPass >= 0.1 ): 
        # line = line + " $ %.1f " % nLQPass + "\pm %.1f $ & " % eNLQPass 
        line = line + " $ %.2f " % nLQPass + "\pm %.2f $ & " % eNLQPass 
    else:
        # line = line + " $ %.2f " % nLQPass + "\pm %.2f $ & " % eNLQPass 
        line = line + " $ %.4f " % nLQPass + "\pm %.4f $ & " % eNLQPass 
    total_bkg = 0
    total_eBkg2_lo = 0
    total_eBkg2_hi = 0
    for sampleName in sampleNames:
        npass     = float ( d_sampleName_cutData [ sampleName ][ final_selection_cut_prefix + str ( LQ_mass )  ][ "Npass"    ] )
        eNpass_lo = float ( d_sampleName_cutData [ sampleName ][ final_selection_cut_prefix + str ( LQ_mass )  ][ "errNpass" ] )
        eNpass_hi = eNpass_lo
	
        if "fromdata" in sampleName.lower() :

            if "qcd" in sampleName.lower() : 
                npass     = npass     * QCDScale 
                eNpass_lo = eNpass_lo * QCDScale 
                eNpass_hi = eNpass_hi * QCDScale 

            if "ttb" in sampleName.lower() : 
                npass     = npass     * TTBarScale 
                eNpass_lo = eNpass_lo * TTBarScale 
                eNpass_hi = eNpass_hi * TTBarScale 

            if npass != 0.0 and eNpass_lo >= 0.01:
                # line = line + " $ " + "%.2f" % npass + " \pm " +"%.2f" % eNpass_lo + " $ &"
                line = line + " $ " + "%.3f" % npass + " \pm " +"%.3f" % eNpass_lo + " $ &"
            elif npass != 0.0 and eNpass_lo < 0.01:
                # line = line + " $ " + "%.3f" % npass + " \pm " +"%.4f" % eNpass_lo + " $ &"
                line = line + " $ " + "%.5f" % npass + " \pm " +"%.5f" % eNpass_lo + " $ &"
            else           : 
                eNpass_lo = 0.0
                eNpass_hi = TTBarScale * 1.14
                # line = line + " $ " + "%.2f" % npass + "_{-0.00}^{+%.2f}" % ( eNpass_hi ) + "$ & "
                line = line + " $ " + "%.5f" % npass + "_{-0.00}^{+%.5f}" % ( eNpass_hi ) + "$ & "

            total_bkg   = total_bkg + npass
            total_eBkg2_lo = total_eBkg2_lo + ( eNpass_lo * eNpass_lo ) 
            total_eBkg2_hi = total_eBkg2_hi + ( eNpass_hi * eNpass_hi ) 

        elif sampleName.lower() != "2012abc" : 
            npass     = npass     / 1000.
            eNpass_lo = eNpass_lo / 1000.
            eNpass_hi = eNpass_hi / 1000.
            
            if npass != 0.0: 
                if ( eNpass_lo >= 0.01 ):
                    # line = line + " $ " + "%.2f" % npass + " \pm " +"%.2f" % eNpass_lo + " $ &"
                    line = line + " $ " + "%.3f" % npass + " \pm " +"%.3f" % eNpass_lo + " $ &"
                elif ( eNpass_lo >= 0.001 ):
                    # line = line + " $ " + "%.3f" % npass + " \pm " +"%.3f" % eNpass_lo + " $ &"
                    line = line + " $ " + "%.5f" % npass + " \pm " +"%.5f" % eNpass_lo + " $ &"
                else:
                    # line = line + " $ " + "%.4f" % npass + " \pm " +"%.4f" % eNpass_lo + " $ &"
                    line = line + " $ " + "%.5f" % npass + " \pm " +"%.5f" % eNpass_lo + " $ &"
            elif "WJet" in sampleName and npass < 1:
                
                eNpass_lo = 0.0
                eNpass_hi = 0.935
                line = line + " $ " + "%.2f" % npass + "_{-0.00}^{+%.2f}" % ( eNpass_hi ) + "$ & "
            else           : 
                if   "WJet" in sampleName: 
                    eNpass_lo = 0.0
                    eNpass_hi = 0.935
                    # line = line + " $ " + "%.2f" % npass + "_{-0.00}^{+%.2f}" % ( eNpass_hi ) + "$ & "
                    line = line + " $ " + "%.3f" % npass + "_{-0.00}^{+%.3f}" % ( eNpass_hi ) + "$ & "
                elif "ZJet" in sampleName: 
                    eNpass_lo = 0.0
                    eNpass_hi = zscale * 1.14
                    # line = line + " $ " + "%.2f" % npass + "_{-0.00}^{+%.2f}" % ( eNpass_hi ) + "$ & "
                    line = line + " $ " + "%.3f" % npass + "_{-0.00}^{+%.3f}" % ( eNpass_hi ) + "$ & "
                else                     : 
                    # line = line + " $ " + "%.2f" % npass + "_{-0.00}^{NO INFO} $ & " 
                    line = line + " $ " + "%.3f" % npass + "_{-0.00}^{NO INFO} $ & " 
            
            total_bkg      = total_bkg + npass
            total_eBkg2_lo = total_eBkg2_lo + ( eNpass_lo * eNpass_lo ) 
            total_eBkg2_hi = total_eBkg2_hi + ( eNpass_hi * eNpass_hi ) 
            
        else:
            line = line + "%d" % npass + " &"
    total_eBkg_lo = math.sqrt ( total_eBkg2_lo ) 
    total_eBkg_hi = math.sqrt ( total_eBkg2_hi ) 
    # if ( total_eBkg_lo == total_eBkg_hi ) : line = line + " $ %.1f \pm " % total_bkg + "%.1f" % total_eBkg_lo  + " $ \\\ \n"
    if ( total_eBkg_lo == total_eBkg_hi ) : line = line + " $ %.2f \pm " % total_bkg + "%.2f" % total_eBkg_lo  + " $ \\\ \n"
    # else                                  : line = line + " $ " + "%.2f" % total_bkg + "_{-%.2f}" % total_eBkg_lo + "^{+%.2f}" % total_eBkg_hi + "$ \\\ \n"
    else                                  : line = line + " $ " + "%.3f" % total_bkg + "_{-%.3f}" % total_eBkg_lo + "^{+%.3f}" % total_eBkg_hi + "$ \\\ \n"

    latex_file.write ( line ) 
    
latex_file.write("  \hline \n")
latex_file.write("\end{tabular}\n")
latex_file.write("\end{table}\n")
latex_file.write("\end{document}\n")
latex_file.close()

os.system ( "pdflatex " + latex_file_name ) 

pdf_file_name = os.environ["PWD"] + "/" + latex_file_name
pdf_file_name = pdf_file_name.replace (".tex", ".pdf" )

print "ls " + pdf_file_name
print "open " + pdf_file_name.replace("/afs/cern.ch/user/e/eberry/", "/mnt/lxplus/")
