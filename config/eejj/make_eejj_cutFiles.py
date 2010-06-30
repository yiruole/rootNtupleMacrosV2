#!/usr/bin/env python

import os

#### USER"S INPUTS HERE ####
inputFileName="../cutTable_eejjSample.txt" # Will be used as a template and left untouched
## St and Mee values - One cut file will be created in this dir for each St,Mee pair
StValues =["140", "240", "300"] # Opt values for MLQ=100,200,300 GeV at ILum=100 nb-1
MeeValues=["100", "100", "100"] # Optimized Mee was 95, but we set it to 100 GeV
#### End OF USER"S INPUTS ####


def replaceCutValue(filename, variable, newvalue) :
    f = file(filename)
    newlines = []
    for line in f:
        replaceline=False
        newline = ""
        for i,field in enumerate(line.split()):
            #print "i=", i , "  field= " , field                
            if (field == variable) :
                replaceline=True
            if (i==1 and replaceline) :
                field=newvalue
            if (replaceline) :
                if (i==0) :
                    newline = str(field)
                elif (i==1) :
                    newline += "\t\t\t\t"+str(field)
                elif (i==6) :
                    newline += "\t"+str(field)
                elif (i==7 or i==8) :
                    newline += " "+str(field)
                else :
                    newline += "\t\t"+str(field)
        if (replaceline) :
            line=newline+"\n"
        newlines.append(line)
    outfile = file(filename, 'w')
    outfile.writelines(newlines)

for i, unused in enumerate(StValues) :
    outFileName=os.path.basename(inputFileName)
    outFileName=outFileName.strip(".txt")
    outFileName += "_Mee"+MeeValues[i]+"_St"+StValues[i]+".txt"
    os.system("cp "+inputFileName+" "+outFileName)
    replaceCutValue(outFileName,"sT",StValues[i])
    replaceCutValue(outFileName,"Mee",MeeValues[i])


