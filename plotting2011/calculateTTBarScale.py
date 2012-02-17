from ROOT import *
import os, math, sys
from optparse import OptionParser

def GetIntegralTH1( histo, xmin, xmax):
    #get integral
    axis = histo.GetXaxis()
    bmin = axis.FindBin(xmin)
    bmax = axis.FindBin(xmax)
    bminResidual = histo.GetBinContent(bmin)*(xmin-axis.GetBinLowEdge(bmin)) / axis.GetBinWidth(bmin)
    bmaxResidual = histo.GetBinContent(bmax)*(axis.GetBinUpEdge(bmax)-xmax) / axis.GetBinWidth(bmax)
    integral = histo.Integral(bmin,bmax) - bminResidual - bmaxResidual
    return integral

def GetErrorIntegralTH1( histo, xmin, xmax):
    print "## calculating error for integral of histo " + str(histo)
    print "## in the x range [" + str(xmin) + "," + str(xmax) + "]"
    #get integral
    axis = histo.GetXaxis()
    bmin = axis.FindBin(xmin)
    bmax = axis.FindBin(xmax)
    bminResidual = histo.GetBinContent(bmin)*(xmin-axis.GetBinLowEdge(bmin)) / axis.GetBinWidth(bmin)
    bmaxResidual = histo.GetBinContent(bmax)*(axis.GetBinUpEdge(bmax)-xmax) / axis.GetBinWidth(bmax)
    integral = histo.Integral(bmin,bmax) - bminResidual - bmaxResidual
    error = 0
    for bin in range(bmin, bmax+1):
	print "bin: " +str(bin)
        if(bin==bmax and bmaxResidual==histo.GetBinContent(bmax)): # skip last bin if out of range
            print "     --> skip bin: " + str(bin)
        else:
            error = error + histo.GetBinError(bin)**2
            print "error**2 : " + str(error)

    error = sqrt(error)
    print  " "
    return error

# def GetIntegralAndError ( histo, xmin, xmax ) : 
#     integral = GetErrorIntegralTH1( histo, xmin, xmax)
#     error    = GetIntegralTH1     ( histo, xmin, xmax)
#     return integral, error

def getIntAndError ( hist ) :
    error = Double(0.0)
    integral = hist.IntegralAndError ( 0, hist.GetNbinsX() + 1, error )
    return integral, error 

usage = "usage: %prog [options] \nExample: python calculateTTBarScale.py -p " + os.environ['LQDATA'] + "/enujj_analysis/enujj/output_cutTable_lq_enujj/analysisClass_lq_enujj_plots.root -q " + os.environ['LQDATA'] + "/enujj_analysis/enujj_qcd/output_cutTable_lq_enujj/analysisClass_lq_enujj_QCD_plots.root -x config/xsection_7TeV_Summer2011.txt -s TTJets_TuneZ2_7TeV-madgraph-tauola -n MTenu_50_110_4jet -l 1.0"

parser = OptionParser(usage=usage)

parser.add_option("-p", "--preselectionFile", dest="preselection_file_name",
                  help="Monte Carlo and Data preselection file name",
                  metavar="PRESELECTION")

parser.add_option("-q", "--qcdFile", dest="qcd_file_name",
                  help="Monte Carlo and Data preselection file name",
                  metavar="PRESELECTION")

parser.add_option("-x", "--xsectionFile", dest="xsection_file_name",
                  help="Original, unscaled cross section file",
                  metavar="XSECTION")

parser.add_option("-s", "--sampleName", dest="sample_name",
                  help="Name of your Z sample",
                  metavar="SAMPLENAME" )

parser.add_option("-n","--plotName", dest="plot_name",
                  help="Name of the plot in your root file to be used for normalizing",
                  metavar="PLOTNAME")

parser.add_option("-l", "--qcdScale", dest="qcd_scale",
                  help="Scale factor for QCD histograms",
                  metavar="QCDSCALE")

parser.add_option("-r", "--useSherpa", dest="use_sherpa",
                  help="Whether to use SHERPA for W/Z samples (default is 'no')",
                  metavar="USESHERPA")

parser.add_option("-f","--newXSectionFile", dest="new_file",
                  help="Name of the new cross section file",
                  metavar="NEWFILE" )


(options, args) = parser.parse_args()

if len ( sys.argv ) < 12 : 
    print usage 
    sys.exit()

preselection_file_name = options.preselection_file_name
qcd_file_name          = options.qcd_file_name         
xsection_file_name     = options.xsection_file_name    
sample_name            = options.sample_name           
plot_name              = options.plot_name    
use_sherpa             = options.use_sherpa
qcd_scale              = float ( options.qcd_scale )

if options.use_sherpa:
    use_sherpa = options.use_sherpa.lower()

other_background_samples = ["ZJet_Madgraph"     ,"PhotonJets"     ,"SingleTop" ,"DIBOSON"     , "WJet_Madgraph"     ]
other_background_titles  = ["$Z/Z^{*}$ + jets"  ,"$\gamma$ + jets","single top","WW + WZ + ZZ", "$W/W^{*}$ + jets"  ]

if options.use_sherpa:
    if use_sherpa == "yes" :
        other_background_samples = ["ZJet_Sherpa"       ,"PhotonJets"     ,"SingleTop" ,"DIBOSON"     , "WJet_Sherpa"       ]
        other_background_titles  = ["$Z/Z^{*}$ + jets"  ,"$\gamma$ + jets","single top","WW + WZ + ZZ", "$W/W^{*}$ + jets"  ]

scale_background_sample  = "TTbar_Madgraph"
scale_background_title   = "$t\\bar{t}$"   

data_sample              = ["DATA"]

other_background_hists = []
other_background_ints  = []
other_background_eInts  = []

xmin = 40
xmax = 120

preselection_file = TFile ( preselection_file_name ) 
qcd_file = TFile ( qcd_file_name ) 

data_hist = preselection_file.Get("histo1D__DATA__" + plot_name )
data_int, data_eInt = getIntAndError ( data_hist ) 

qcd_background_hist = qcd_file.Get("histo1D__DATA__" + plot_name )
qcd_background_hist.Scale ( qcd_scale ) 
qcd_background_int, qcd_background_eInt = getIntAndError ( qcd_background_hist )

scale_hist = preselection_file.Get("histo1D__" + scale_background_sample + "__" + plot_name )
scale_int, scale_eInt = getIntAndError ( scale_hist )

other_backgrounds_int = 0.0
other_backgrounds_eInt = 0.0

for sample in other_background_samples:
    hist_name = "histo1D__"+ sample  +"__" + plot_name
    hist = TH1D()
    hist = preselection_file.Get( hist_name ) 
    other_background_hists.append ( hist ) 

    hist_int, hist_eInt = getIntAndError ( hist )

    other_background_ints.append ( hist_int ) 
    other_background_eInts.append ( hist_eInt ) 

    other_backgrounds_int = other_backgrounds_int + hist_int
    other_backgrounds_eInt = other_backgrounds_eInt + ( hist_eInt * hist_eInt ) 

other_backgrounds_int  = other_backgrounds_int + qcd_background_int 
other_backgrounds_eInt = other_backgrounds_eInt + (  qcd_background_eInt *  qcd_background_eInt )

other_backgrounds_eInt = math.sqrt ( other_backgrounds_eInt ) 

Rnum = data_int - other_backgrounds_int 
eRnum = math.sqrt ( ( data_eInt * data_eInt ) + ( other_backgrounds_eInt * other_backgrounds_eInt ) )
                           
Rden = scale_int 
eRden = scale_eInt

R = Rnum / Rden 
eR = R * math.sqrt ( ( ( eRnum / Rnum ) * ( eRnum / Rnum ) ) + 
                     ( ( eRden / Rden ) * ( eRden / Rden ) ) )

latex_file_name = "table_TTBarScale.tex"

if options.use_sherpa:
    if use_sherpa == "yes" :
        latex_file_name = "table_TTBarScale_WZSherpa.tex"

latex_file = open ( latex_file_name, "w" )  

latex_file.write("\documentclass{article}\n")
latex_file.write("\usepackage{multirow} \n")
latex_file.write("\usepackage[landscape, top=1cm, bottom=1cm, left=1cm, right=1cm]{geometry}\n")
latex_file.write("\pagestyle{empty} \n")
latex_file.write("\\begin{document}\n")
latex_file.write("\\begin{table} \n")
latex_file.write("\\begin{tabular}{ l | c  } \n")
latex_file.write ("Title & Value \\\ \n")
latex_file.write ("\\hline \n")
latex_file.write ("\\hline \n")
for i, title in enumerate ( other_background_titles ) :
    latex_file.write( title + " & " + "%.3f" %  other_background_ints[i] + "\t$\pm$\t" + "%.3f" % other_background_eInts[i] + " \\\ \n" )
latex_file.write ( "QCD & " + "%.3f" % qcd_background_int + "\t$\pm$\t"  + "%.3f" % qcd_background_eInt + " \\\ \n" )

latex_file.write ("\\hline \n")
latex_file.write ("Data " + "&" + str ( int ( data_int ) ) + "\t$\pm$\t"  + str(int( data_eInt) ) + " \\\ \n" )
latex_file.write ( scale_background_title + " & "  + "%.3f" % scale_int + "\t$\pm$\t"  + "%.3f" % scale_eInt + " \\\ \n" )
latex_file.write ( "Other Backgrounds & " + "%.3f" % other_backgrounds_int + "\t$\pm$\t"  + "%.3f" % other_backgrounds_eInt + " \\\ \n" )
latex_file.write ("\\hline \n")
latex_file.write ( "$R_{t\\bar{t}}$ & " + "%.3f" % R + "$\pm$" + "%.3f" % eR + " \\\ \n" )

latex_file.write("\end{tabular}\n")
latex_file.write("\end{table}\n")
latex_file.write("\end{document}\n")

latex_file.close()

os.system ( "pdflatex " + latex_file_name )

command = "ls "+latex_file_name[:-4]+".* | egrep -v \"tex|pdf\" | xargs rm " 

os.system ( command ) 

print "\n\n I made this latex table: "
print " open " + os.getcwd() + "/" + latex_file_name.replace(".tex", ".pdf")

print "\n\n According to preselection file: "
os.system ("ls -l " + preselection_file_name )

print "\n\n And QCD file: "
os.system ("ls -l " + qcd_file_name )

print "\n\n I calculate: "
print " R_TTBar = " + "%.3f" % R + " +/- " + "%.3f" % eR 

if not options.new_file:
    new_xsection_file_name = xsection_file_name.replace(".txt", "_TTBarScaled.txt")
else :
    new_xsection_file_name = options.new_file

xsection_file = open ( xsection_file_name, "r" )
new_xsection_file = open (new_xsection_file_name, "w" ) 

for line in xsection_file:
    
    if len ( line.split() ) < 2 : continue

    this_sample_name = line.split()[0]
    this_xsection    = line.split()[1]

    if sample_name not in this_sample_name:
        new_xsection_file.write( line.strip() + "\n" )
    else :
        new_line = this_sample_name + " \t " + str ( float ( this_xsection )  * R )
        new_xsection_file.write( new_line.strip() + "\n" )

xsection_file.close()
new_xsection_file.close()

print "\n\n From this xsection file: "
print " " + xsection_file_name

print "\n\n I made a new xsection file: " 
print " " + new_xsection_file_name

print "\n\n With this sample name scaled:"
print " " + sample_name

print "\n\n"
