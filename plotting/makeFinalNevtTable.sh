#!/bin/bash

##########  USER'S INPUTS BEGIN HERE ############################################################

#FILE=$LQDATA/eejj/10.9pb-1/output_cutTable_eejjSample/analysisClass_eejjSample_tables.dat
FILE=$LQDATA/eejj/34.7pb-1/output_cutTable_eejjSample_noPreStCut_ZjetsRescaled/analysisClass_eejjSample_tables.dat

MASSES="250 280 300 320 340 370 400"

# taken from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Leptoquarks/CommonToolsV2/test/exclusion_plot_beta_vs_m_2010Note.C
# Number and order of the limit values must correspond to number and order of MASSES above !!! 
xsUp_observed=( 0.364746 0.290527 0.277832 0.268799 0.263184 0.252686 0.244141 )
xsUp_expected=( 0.424484 0.34855  0.310529 0.285282 0.262771 0.239716 0.219592 )



##########  USER'S INPUTS END HERE ###############################################################


ii=0  
for MASS in $MASSES
do
  CUT="sT_MLQ$MASS"
  LQSAMPLE="LQeejj_M$MASS"
  SAMPLES="$LQSAMPLE TTbar_Madgraph ZJetAlpgen OTHERBKG ALLBKG DATA"
  
  ALLSTRINGS=""
  for SAMPLE in $SAMPLES
  do
    STRING=`awk -v factor=0.001 -v cut=$CUT '{if( NF==1 ) name=$1; if( $1==cut ) printf("%20s %20s  %i   %f %f        %f %f \n",name,$1, $2, $6*factor,$7*factor,$10,$11) }' $FILE | egrep $SAMPLE`
    if [ $SAMPLE == $LQSAMPLE ]; then
	THISSTRING=`echo $STRING | awk -v mass=$MASS '{printf("%i ($S_T>$%i) & %.3g$pm$%.1g & %.3g$pm$%.1g ",mass,$3,$4,$5,$6,$7)}'`
    elif [ $SAMPLE == "DATA" ]; then
	THISSTRING=`echo $STRING | awk '{printf("& %i ", $4*1000)}'`
    else
	THISSTRING=`echo $STRING | awk '{printf("& %.2g$pm$%.1g ",$4,$5)}'`
    fi
    ALLSTRINGS="$ALLSTRINGS $THISSTRING"
  done
  LIMITSTRING=`echo ${xsUp_observed[$ii]}" "${xsUp_expected[$ii]} | awk '{printf("& %.3g / %.3g \n",$1,$2)}' `
  FINALSTRING=$ALLSTRINGS" "$LIMITSTRING
  let "ii += 1"
  #echo $ALLSTRINGS'\\' 
  echo $FINALSTRING'\\' 
done


