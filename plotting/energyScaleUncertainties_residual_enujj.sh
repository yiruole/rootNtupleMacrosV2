#!/bin/sh

# It calculates the residual effect due to JES uncertainties between 2 cuts (set in the variable CUTS).
# Important: if you do not have al samples required below, set the missing files to be the same as FILE________STD
# so you'll get zeros as results for the missing files, but will not mess up the results for the available files.

##########  USER'S INPUTS BEGIN HERE ############################################################

FILE________STD=/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/36.0pb-1_presel_MET45_presel_sT250_Wrescale1.18_Feb112011/analysisClass_enujjSample_tables.dat
FILE_EES_EB_POS=/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/36.0pb-1_presel_MET45_presel_sT250_Wrescale1.18_EES_EB1.01_FullNtuples_Mar012011/analysisClass_enujjSample_tables.dat
FILE_EES_EB_NEG=/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/36.0pb-1_presel_MET45_presel_sT250_Wrescale1.18_EES_EB0.99_FullNtuples_Mar012011/analysisClass_enujjSample_tables.dat
FILE_EES_EE_POS=/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/36.0pb-1_presel_MET45_presel_sT250_Wrescale1.18_EES_EE1.03_FullNtuples_Mar012011/analysisClass_enujjSample_tables.dat
FILE_EES_EE_NEG=/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/36.0pb-1_presel_MET45_presel_sT250_Wrescale1.18_EES_EE0.97_FullNtuples_Mar012011/analysisClass_enujjSample_tables.dat
FILE_JES____POS=/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/36.0pb-1_presel_MET45_presel_sT250_Wrescale1.18_JES1.05_FullNtuples_Mar012011/analysisClass_enujjSample_tables.dat
FILE_JES____NEG=/home/ferencek/work/Leptoquarks/output_fromAFS/enujj_analysis/36.0pb-1_presel_MET45_presel_sT250_Wrescale1.18_JES0.95_FullNtuples_Mar012011/analysisClass_enujjSample_tables.dat

## LQ M300 signal
# SAMPLES='LQenujj_M300'
# CUTS='nocut sT_MLQ300' # first cut is the one where data-MC comparison is done, second cut is the final cut

## TTbar inputs (preselection):
# SAMPLES='TTbar_Madgraph'
# CUTS='maxMEtaJets_noOvrlp_ID sT_MLQ300' # first cut is the one where data-MC comparison is done, second cut is the final cut

## W+jets inputs (preselection + preliminary sT cut)
# SAMPLES='WJetAlpgen'
# CUTS='nMuon_PtCut_IDISO sT_MLQ300' # first cut is the one where data-MC comparison is done, second cut is the final cut

## Other backgrounds
SAMPLES='OTHERBKG'
CUTS='nocut sT_MLQ300' # first cut is the one where data-MC comparison is done, second cut is the final cut

##########  USER'S INPUTS END HERE ###############################################################


files="$FILE________STD $FILE_EES_EB_POS $FILE_EES_EB_NEG $FILE_EES_EE_POS $FILE_EES_EE_NEG $FILE_JES____POS $FILE_JES____NEG"

outfile='temp.txt'
echo "" > $outfile

for CUT in $CUTS
do
  #echo "The numbers below are the relative change of each energy scale variation w.r.t the default energy scale" > $outfile
  echo "    SAMPLES            WHAT      DEFAULT        EES_EB_POS      EES_EB_NEG      EES_EE_POS      EES_EE_NEG      JES____POS      JES____NEG" >> $outfile
  for SAMPLE in $SAMPLES
    do
    NEVTS=""
    for file in $files
      do
      #echo $file
      #echo "              Sample               cut        PassedEvents+Error       Efficiency+Error">> $outfile
      STRING=`awk -v factor=0.001 -v cut=$CUT '{if( NF==1 ) name=$1; if( $1==cut ) printf("%20s %20s     %f %f \n",name,$1,$6*factor,$7*factor) }' $file | egrep $SAMPLE`
      THISNEVTS=`echo $STRING | awk '{print $3}'`
      NEVTS="$NEVTS $THISNEVTS"
    done
    #echo $NEVTS | awk -v sample=$SAMPLE '{printf("%20s  RelChange\t%f\t%f\t%f\t%f\t%f\t%f\t%f \n",sample, $1/$1-1, $2/$1-1, $3/$1-1, $4/$1-1, $5/$1-1, $6/$1-1, $7/$1-1) }' >> $outfile
    echo $NEVTS | awk -v sample=$SAMPLE '{printf("%20s  nEvtsPass\t%f\t%f\t%f\t%f\t%f\t%f\t%f \n",sample, $1     , $2     , $3     , $4     , $5     , $6     , $7     ) }' >> $outfile
  done
done

echo ""
echo "Events passed for after these 2 cuts $CUTS are, respectively:"
cat temp.txt

echo ""
echo "The variation in the relative efficiencies are:"
echo ""
echo "    SAMPLES           WHAT                      DEFAULT        EES_EB_POS      EES_EB_NEG      EES_EE_POS      EES_EE_NEG      JES____POS      JES____NEG"
for SAMP in $SAMPLES
  do
  TWOLINES=`cat $outfile | grep $SAMP`
  #echo $TWOLINES
  EFFRATIOS=`echo $TWOLINES | awk -v sample=$SAMP '{printf("%20s  effRatio\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e \n",sample, $12/$3, $13/$4, $14/$5, $15/$6, $16/$7, $17/$8, $18/$9  ) }'`
  #echo $EFFRATIOS
  echo $EFFRATIOS | awk -v sample=$SAMP '{printf("%20s  EffRatio_variation\t%f\t%f\t%f\t%f\t%f\t%f\t%f \n",sample, $3/$3-1, $4/$3-1, $5/$3-1, $6/$3-1, $7/$3-1, $8/$3-1, $9/$3-1) }'
done

