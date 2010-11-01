#!/bin/sh

##########  USER'S INPUTS BEGIN HERE ############################################################

FILE________STD=$LQDATA/eejj/10.9pb-1/output_cutTable_eejjSample_stdEneScale/analysisClass_eejjSample_tables.dat
FILE_EES_EB_POS=$LQDATA/eejj/10.9pb-1/output_cutTable_eejjSample_EES_EB_1.01/analysisClass_eejjSample_tables.dat
FILE_EES_EB_NEG=$LQDATA/eejj/10.9pb-1/output_cutTable_eejjSample_EES_EB_0.99/analysisClass_eejjSample_tables.dat
FILE_EES_EE_POS=$LQDATA/eejj/10.9pb-1/output_cutTable_eejjSample_EES_EE_1.03/analysisClass_eejjSample_tables.dat
FILE_EES_EE_NEG=$LQDATA/eejj/10.9pb-1/output_cutTable_eejjSample_EES_EE_0.97/analysisClass_eejjSample_tables.dat
FILE_JES____POS=$LQDATA/eejj/10.9pb-1/output_cutTable_eejjSample_JES_1.10/analysisClass_eejjSample_tables.dat
FILE_JES____NEG=$LQDATA/eejj/10.9pb-1/output_cutTable_eejjSample_JES_0.90/analysisClass_eejjSample_tables.dat

SAMPLES='LQeejj_M300 TTbar_Madgraph ZJetAlpgen ALLBKG'

##########  USER'S INPUTS END HERE ###############################################################

#files=`ls $LQDATA/eejj/10.9pb-1/output_cutTable_eejjSample_stdEneScale/analysisClass_eejjSample_tables.dat $LQDATA/eejj/10.9pb-1/output_cutTable_eejjSample_[E,J]ES*/analysisClass_eejjSample_tables.dat
files="$FILE________STD $FILE_EES_EB_POS $FILE_EES_EB_NEG $FILE_EES_EE_POS $FILE_EES_EE_NEG $FILE_JES____POS $FILE_JES____NEG"


# for file in $files
# do
#   echo $file
#   echo "              Sample               cut        PassedEvents+Error       Efficiency+Error"
#   awk -v factor=0.001 '{if( NF==1 ) name=$1; if( $1=="sT_MLQ300" ) printf("%20s %20s     %f %f        %f %f \n",name,$1,$6*factor,$7*factor,$10,$11) }' $file | egrep 'LQeejj_M300 | TTbar_Madgraph | ZJetAlpgen | ALLBKG'
# done


echo "The numbers below are the relative change of each energy scale variation w.r.t the default energy scale"
echo "    SAMPLES            WHAT      DEFAULT        EES_EB_POS      EES_EB_NEG      EES_EE_POS      EES_EE_NEG      JES____POS      JES____NEG"
for SAMPLE in $SAMPLES
do
  NEVTS=""
  for file in $files
  do
    #echo $file
    #echo "              Sample               cut        PassedEvents+Error       Efficiency+Error"
    STRING=`awk -v factor=0.001 '{if( NF==1 ) name=$1; if( $1=="sT_MLQ300" ) printf("%20s %20s     %f %f \n",name,$1,$6*factor,$7*factor) }' $file | egrep $SAMPLE`
    THISNEVTS=`echo $STRING | awk '{print $3}'`
    NEVTS="$NEVTS $THISNEVTS"
  done
  echo $NEVTS | awk -v sample=$SAMPLE '{printf("%20s  RelChange\t%f\t%f\t%f\t%f\t%f\t%f\t%f \n",sample, $1/$1-1, $2/$1-1, $3/$1-1, $4/$1-1, $5/$1-1, $6/$1-1, $7/$1-1) }'
#  echo $NEVTS | awk -v sample=$SAMPLE '{printf("%20s  nEvtsPass\t%f\t%f\t%f\t%f\t%f\t%f\t%f \n",sample, $1     , $2     , $3     , $4     , $5     , $6     , $7     ) }'
done




