#!/usr/bin/python

import sys
import os
import csv
from math import sqrt 

#--------------------------------------------------------------------------------------
# Define functions
#--------------------------------------------------------------------------------------

def printUsage() : 
    print "Usage: "
    print "  python makeLatexTable.py <.dat file> " 
    sys.exit() 

def printChannelInfo ( channel_name, row_names, column_names, channel_data ):


    column_name_header = ""
    for column_name in column_names:
        column_name_header = column_name_header + "\t" + column_name
    for row_name in row_names:
        row = "" 
        for column_name in column_names:
            data = channel_data[channel_name][row_name][column_name]
            row = row + "\t" + str(data) 
        print row

#-------------------------------------------------------
# There should be only one argument (input list)
#-------------------------------------------------------

if len (sys.argv) != 2 and len (sys.argv) != 3 : printUsage()

txt_file_path = sys.argv[1]

qcd_file_path = ""
if len (sys.argv) == 3 : qcd_file_path = sys.argv[2]

#-------------------------------------------------------
# Does the file exist?
#-------------------------------------------------------

if not os.path.exists ( txt_file_path ) :
    print "ERROR: file \""+txt_file_path+"\" does not exist"
    printUsage()

if len ( sys.argv ) == 3 and not os.path.exists ( qcd_file_path ) :
    print "ERROR: file \""+qcd_file_path+"\" does not exist"
    printUsage()

channel_names_to_include  = ["ZJet_Madgraph","PhotonJets"     ,"SingleTop" ,"DIBOSON"     ,"TTbar_Madgraph"   ,"WJet_Madgraph"]
channel_titles_to_include = ["Z/Z* + jets"  ,"$\gamma$ + jets","single top","WW + WZ + ZZ","$t\\bar{t}$"      , "W/W* + jets" ]

# channel_names_to_include  = ["TTbar_Madgraph"]
# channel_titles_to_include = ["$t\\bar{t}$"   ]

#-------------------------------------------------------
# open file and scan contents
#-------------------------------------------------------

txt_file = open ( txt_file_path, "r" ) 
txt_file_start = txt_file.tell() 
channel_data  = {}
all_channel_names = []
channel_names = []
column_names  = []
row_names     = []

empty_channel_names = []

for line in txt_file : 
    if len ( line.split() ) == 1 : all_channel_names.append ( line.strip () ) 

txt_file.seek ( txt_file_start ) 
txt_file_data = txt_file.read()

for i,channel_name in enumerate(all_channel_names):
    if i != ( len ( all_channel_names ) - 1 ): 
        channel_raw_data = txt_file_data.split ( all_channel_names[i] )[1].split(all_channel_names[i+1])[0]
    else :
        channel_raw_data = txt_file_data.split ( all_channel_names[i] )[1]    
    
    channel_raw_data_rows = channel_raw_data.strip().split("\n")
    n_channel_raw_data_rows = len ( channel_raw_data_rows ) 

    if n_channel_raw_data_rows == 1 : 
        print "WARNING: " + channel_name + " is empty, and I will not consider it." 
        empty_channel_names.append( channel_name ) 
        continue
    channel_names.append ( channel_name ) 
    channel_data[channel_name] = {}

    column_names = channel_raw_data_rows[0].split()

    for irow, channel_raw_data_row in enumerate(channel_raw_data_rows[1:]):
        fields = channel_raw_data_row.split()
        row_name = fields[0]
        if row_name not in row_names : row_names.append ( row_name )
        channel_raw_data_row_tuple = zip ( column_names, fields ) 
        channel_raw_data_row_dict = {} 
        for i,(name, value) in enumerate(channel_raw_data_row_tuple):
            if i != 0 and value.strip() != "-":
                new_value = float ( value )
                channel_raw_data_row_dict [str(name)] = new_value 
            else :
                new_value = value.replace ( "\n", " " ) 
                channel_raw_data_row_dict [str(name)] = str(new_value.strip())

        channel_data [channel_name][row_name] = channel_raw_data_row_dict

#-------------------------------------------------------
# QCD file:
#-------------------------------------------------------

if len ( sys.argv ) == 3 : 
    qcd_file = open ( qcd_file_path, "r" ) 
    qcd_file_data = qcd_file.read() 
    
    qcd_raw_data = qcd_file_data.split("DATA")[1].split("\n\n")[0].strip()
    
    qcd_raw_data_rows = qcd_raw_data.split("\n")
    n_qcd_raw_data_rows = len ( qcd_raw_data_rows ) 
    
    if n_qcd_raw_data_rows == 1 : 
        print "WARNING: DATA is empty in " + qcd_file_path
        sys.exit()

    channel_names.append ( "QCD" )
    channel_data["QCD"] = {}

    
    channel_names_to_include .append ("QCD")
    channel_titles_to_include.append ("QCD")
    
    for irow, qcd_raw_data_row in enumerate(qcd_raw_data_rows[1:]):
        fields = qcd_raw_data_row.split()
        row_name = fields[0]
        qcd_raw_data_row_tuple = zip ( column_names, fields ) 
        qcd_raw_data_row_dict = {} 
        for i,(name, value) in enumerate(qcd_raw_data_row_tuple):
            if i != 0 and value.strip() != "-":
                new_value = float ( value )
                qcd_raw_data_row_dict [str(name)] = new_value 
            else :
                new_value = value.replace ( "\n", " " ) 
                qcd_raw_data_row_dict [str(name)] = str(new_value.strip())

        channel_data["QCD"][row_name] = qcd_raw_data_row_dict

for j, empty_channel_name in enumerate ( empty_channel_names ) :
    for i, channel_name_to_include in enumerate ( channel_names_to_include ) :
        if empty_channel_name == channel_name_to_include:
            channel_names_to_include.remove(channel_name_to_include)
            channel_titles_to_include.remove(channel_titles_to_include[i])

latex_file_name = "table.tex"
latex_file = open ( latex_file_name, "w" )  

latex_file.write("\documentclass{article}\n")
latex_file.write("\usepackage{multirow} \n")
latex_file.write("\usepackage[landscape, top=1cm, bottom=1cm, left=1cm, right=1cm]{geometry}\n")
latex_file.write("\pagestyle{empty} \n")
latex_file.write("\\begin{document}\n")
latex_file.write("\\begin{table} \n")
latex_file.write("\\small \n")
line = "\\begin{tabular}{| l | "
for channel_title_to_include in channel_titles_to_include:
    line = line + " c" 
line = line + " c | c |} \n"
latex_file.write ( line ) 
latex_file.write("  \hline \n")
if len ( channel_names_to_include ) != 0 :
    latex_file.write("  Cut & \multicolumn{"+str ( len (channel_names_to_include ) + 1 )  + "}{c|}{MC and QCD Background Samples} & Events \\\ \n" )
    latex_file.write("      & \multicolumn{"+str ( len (channel_names_to_include ) + 1 )  + "}{c|}{Selected Events in}            & in     \\\ \n" )
else : 
    latex_file.write("  Cut & Events \\\ \n" )
    latex_file.write("      & in     \\\ \n" )
line = "     "
for channel_title_to_include in channel_titles_to_include:
    line = line + " & " + channel_title_to_include 
if len ( channel_names_to_include ) != 0 :
    line = line + " & All Bkgds & Data \\\ \n"
else : 
    line = line + " & Data \\\ \n"
latex_file.write ( line ) 
latex_file.write("  \hline \n")
for i,row_name in enumerate(row_names):
    if row_name == "variableName" : continue
    variable = row_name.replace("_","\_")
    
    if len ( channel_names_to_include ) > 0 : 
        cut_min1 = str ( channel_data[channel_names_to_include[0]][row_name]["min1"] )
        cut_max1 = str ( channel_data[channel_names_to_include[0]][row_name]["max1"] )
    else :
        cut_min1 = str ( channel_data["DATA"][row_name]["min1"] )
        cut_max1 = str ( channel_data["DATA"][row_name]["max1"] )

    line = ""

    variable_tex = variable

    if variable_tex == "nocut" or variable_tex == "skim":
        line = variable_tex
    elif cut_min1.strip() == "-inf" and cut_max1.strip() == "+inf" and variable != "Reweighting":
        continue
    elif cut_min1.strip() != "-inf" and cut_max1.strip() == "+inf":
        line = variable_tex + " $>$ " + cut_min1  
    elif cut_min1.strip() == "-inf" and cut_max1.strip() != "+inf":
        line = variable_tex + " $\leq$ " + cut_max1 
    elif cut_min1.strip() != "-inf" and cut_max1.strip() != "+inf":
        line = cut_min1 + " $<$ " + variable_tex + " $\leq$ " + cut_max1 

    total_data = 0
    total_eDataSqr = 0
    for i, channel_name in enumerate ( channel_names_to_include ) :
        if channel_name != "QCD" : 
            data     = int   ( channel_data[channel_name][row_name]["Npass"   ] / 1000.0 )
            eData    = int   ( channel_data[channel_name][row_name]["errNpass"] / 1000.0 )
        if channel_name == "QCD" : 
            data     = int   ( channel_data[channel_name][row_name]["Npass"   ]  )
            eData    = int   ( channel_data[channel_name][row_name]["errNpass"]  )
        total_data = data + total_data
        total_eDataSqr = total_eDataSqr + ( eData * eData ) 
        if ( eData != 0 ) : line = line + " & $" + str( (data )) + "\pm"+ str( eData ) + "$"
        if ( eData == 0 ) : line = line + " & $" + str( (data )) + "$"

    total_eData = int ( sqrt ( total_eDataSqr )  )

    if len ( channel_names_to_include ) > 0 : 
        if ( total_eData != 0 ) : line = line + " & $" + str(total_data) + "\pm" + str ( total_eData ) + "$"
        if ( total_eData == 0 ) : line = line + " & $" + str(total_data) + "$"
    if ( "DATA" not in empty_channel_names ) : 
        line = line + " & " + str(int (channel_data["DATA"][row_name]["Npass"])) + " \\\ \n"
    else :
        line = line + " & NO DATA \\\ \n"
    latex_file.write ( line ) 
latex_file.write("  \hline \n")
latex_file.write("\end{tabular}\n")
latex_file.write("\end{table}\n")
latex_file.write("\end{document}\n")

latex_file.close()

os.system ( "pdflatex " + latex_file_name ) 
os.system ( "open " + latex_file_name.replace (".tex",".pdf") )
