#---Import
import sys
import string
from optparse import OptionParser
import os.path
#from ROOT import *
import re
import pprint # for pretty printing
import math

#==============================================================

#--- input file(s)
#inputFile = "*.log"
inputFile = "*Zrescale1.20_Wrescale1.06_W0jets_W1jets_HighPtElePrintOut.log"

#--- matching string for grep
#matchstring_sel = "PassFullSelection"
#matchstring_sel = "PassMETThreshold"
#matchstring_sel = "PassMinMETPt1stEleThreshold"
#matchstring_sel = "PassPt1stEleThreshold"
#matchstring_sel = "PassPt1stJetThreshold"
#matchstring_sel = "PassPt2ndJetThreshold"
#matchstring_sel = "PassMTenuThreshold"
#matchstring_sel = "PasssTThreshold"
#matchstring_sel = "PassPtEleW0jetThreshold"
matchstring_sel = "PassPtEleW1jetThreshold"

#--- matching string to create list of run, ls, event number
matchstring_run_lS_evNum = "Run"

#--- exclude columns in this range
min = 0
max = 7

#==============================================================

#--- create file with event details
outputFile = matchstring_sel+"_details.txt"
os.system("rm -f "+outputFile)
os.system("grep " + matchstring_sel + " " +  inputFile + " | awk -v f="+str(min)+" -v t="+str(max)+" \'{ for (i=1; i<=NF;i++) if( i>=f && i<=t) continue; else printf(\"%s%s\", $i,(i!=NF) ? OFS : ORS) }\' >> "+outputFile)

#--- read run, ls, event number --> put it in a dict --> sort it by run+ls
d = {}
index = 0
for ln, line in enumerate( open(outputFile) ):
    line = line.strip()
    if(  re.search(matchstring_run_lS_evNum, line) ):
        ##print index
        ##index=index+1
        #print ln
        line = line.split()
        run = string.split(line[4],",")[0]
        ls = string.split(line[5],",")[0]
        ev = string.split(line[6],",")[0]
        #print run+":"+ls+":"+ev
        indice = int(100000*int(run))+int(ls)
        d[indice]=run+":"+ls+":"+ev
        #d[10000*int(run)+int(ls)]=run+":"+ls+":"+ev

sum = d.keys()
sum.sort()
#print sum

#--- create file with list of events to be used with "edmPickEvents.py"
fout = open(matchstring_sel+"_list.txt", "w")
for v, value in enumerate(sum):
    #print value    
    print str(v+1) + " " + d[value]
    fout.write( r"%s" % d[value] + "\n" ) 
fout.close()
