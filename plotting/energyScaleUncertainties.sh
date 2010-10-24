#!/bin/sh

files=`ls $LQDATA/collisions/10.9pb-1/output_cutTable_eejjSample/analysisClass_eejjSample_tables.dat $LQDATA/collisions/10.9pb-1/output_cutTable_eejjSample_[E,J]ES*/analysisClass_eejjSample_tables.dat`

for file in $files
do
  echo $file
  echo "              Sample               cut        PassedEvents+Error       Efficiency+Error"
  awk -v factor=0.001 '{if( NF==1 ) name=$1; if( $1=="sT_MLQ300" ) printf("%20s %20s     %f %f        %f %f \n",name,$1,$6*factor,$7*factor,$10,$11) }' $file | egrep 'LQeejj_M300 | TTbar_Madgraph | ZJetAlpgen'
done