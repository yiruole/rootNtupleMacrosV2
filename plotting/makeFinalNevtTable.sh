#!/bin/sh

##########  USER'S INPUTS BEGIN HERE ############################################################

FILE=$LQDATA/collisions/10.9pb-1/output_cutTable_eejjSample/analysisClass_eejjSample_tables.dat

MASSES="250 300 320 340 400"

##########  USER'S INPUTS END HERE ###############################################################


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
	THISSTRING=`echo $STRING | awk -v mass=$MASS '{printf("%i ($S_T>$%i) & %.3g$pm$%.1g & %.3g$pm$%.1g \n",mass,$3,$4,$5,$6,$7)}'`
    elif [ $SAMPLE == "DATA" ]; then
	THISSTRING=`echo $STRING | awk '{printf("& %i \n", $4*1000)}'`
    else
	THISSTRING=`echo $STRING | awk '{printf("& %.3g$pm$%.1g \n",$4,$5)}'`
    fi
    ALLSTRINGS="$ALLSTRINGS $THISSTRING"
  done
  echo $ALLSTRINGS'\\'
done



