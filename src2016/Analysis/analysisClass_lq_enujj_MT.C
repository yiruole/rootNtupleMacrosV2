#define analysisClass_cxx
#include "analysisClass.h"
#include "likelihoodGetter.h"
#include <TH2.h>
#include <TH1F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
// for scale factors
#include "ElectronScaleFactors.C"
#include "MuonScaleFactors.C"
// 2016 trigger efficiency
#include "TriggerEfficiency2016.h"


analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile){}

analysisClass::~analysisClass(){}

void analysisClass::Loop()
{
   std::cout << "analysisClass::Loop() begins" <<std::endl;   
   bool do_eleQuality_plots = false;
   bool do_extra_finalSelection_plots = true;

   char cut_name[100];
   char st_cut_name [100];
   char mt_cut_name [100];
   char mej_cut_name[100];
   char met_cut_name[100];   
   
   //--------------------------------------------------------------------------
   // Final selection mass points
   //--------------------------------------------------------------------------

   const int n_lq_mass = 37;
   int LQ_MASS[n_lq_mass] = { 
     200,  250,
     300,  350,  400, 450, 500, 550,  600,  650,
     700,  750,  800, 850, 900, 950, 1000, 1050,
     1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450,
     1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850,
     1900, 1950, 2000
   };

   // SIC only look at LQ650 selection for now
   // LQ650 only 2012
   //const int n_lq_mass = 1;
   //int LQ_MASS[n_lq_mass] = { 
   //   650
   //};

   // turn off totally for optimization
   //const int n_lq_mass = 0;
   //int LQ_MASS[n_lq_mass];


   std::vector<bool> passed_vector;
   
   std::vector<bool> passed_st_vector;
   std::vector<bool> passed_mt_vector;
   std::vector<bool> passed_met_vector;
   std::vector<bool> passed_mej_vector;

   std::vector<bool> passed_stAndMtAndMet_vector ;
   std::vector<bool> passed_stAndMtAndMej_vector ;
   std::vector<bool> passed_stAndMetAndMej_vector;
   std::vector<bool> passed_mtAndMetAndMej_vector;
   
   std::vector<bool> passed_stAndMej_vector;
   std::vector<bool> passed_mtAndMej_vector;
   std::vector<bool> passed_metAndMej_vector;
   std::vector<bool> passed_stAndMet_vector;
   std::vector<bool> passed_mtAndMet_vector;
   std::vector<bool> passed_stAndMt_vector;
   
   //--------------------------------------------------------------------------
   // Setup likelihood
   //--------------------------------------------------------------------------
   
   std::vector<std::string> variables;
   variables.push_back ( std::string("Mej_PDF"));
   variables.push_back ( std::string("sT_PDF" ));
   variables.push_back ( std::string("MET_PDF"));
   variables.push_back ( std::string("MT_PDF"));
   
   //std::vector<std::string> signals;
   //signals.push_back ( std::string("LQ_M300"));
   //signals.push_back ( std::string("LQ_M350"));
   //signals.push_back ( std::string("LQ_M400"));
   //signals.push_back ( std::string("LQ_M450"));
   //signals.push_back ( std::string("LQ_M500"));
   //signals.push_back ( std::string("LQ_M550"));
   //signals.push_back ( std::string("LQ_M600"));
   //signals.push_back ( std::string("LQ_M650"));
   //signals.push_back ( std::string("LQ_M700"));
   //signals.push_back ( std::string("LQ_M750"));
   //signals.push_back ( std::string("LQ_M800"));
   //signals.push_back ( std::string("LQ_M850"));
   //signals.push_back ( std::string("LQ_M900"));
   //signals.push_back ( std::string("LQ_M950"));
   //signals.push_back ( std::string("LQ_M1000"));
   //signals.push_back ( std::string("LQ_M1050"));
   //signals.push_back ( std::string("LQ_M1100"));
   //signals.push_back ( std::string("LQ_M1150"));
   //signals.push_back ( std::string("LQ_M1200"));

   //std::string pdf_file_name("/afs/cern.ch/user/e/eberry/public/LQPDF/enujj_pdfs_withMT.root");
   //likelihoodGetter * myLikelihoodGetter = new likelihoodGetter ( pdf_file_name, variables, signals ) ;
   //std::vector<double> likelihood_values;
   
   //--------------------------------------------------------------------------
   // Decide which plots to save (default is to save everything)
   //--------------------------------------------------------------------------
   
   fillSkim                         ( !true  ) ;
   fillAllPreviousCuts              ( !true  ) ;
   fillAllOtherCuts                 ( !true  ) ;
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;
   
   //--------------------------------------------------------------------------
   // Get pre-cut values
   //--------------------------------------------------------------------------

   // eta boundaries
   double eleEta_bar            = getPreCutValue1("eleEta_bar");
   double eleEta_end1_min       = getPreCutValue1("eleEta_end1");
   double eleEta_end1_max       = getPreCutValue2("eleEta_end1");
   double eleEta_end2_min       = getPreCutValue1("eleEta_end2");
   double eleEta_end2_max       = getPreCutValue2("eleEta_end2");

   // std::string test_string = getPreCutString1("test");
   // std::cout << "test string = \"" << test_string << "\"" << std::endl;

   //--------------------------------------------------------------------------
   // 2016 trigger efficiency
   //--------------------------------------------------------------------------
   std::string trigEffFileName = getPreCutString1("TriggerEfficiencyFileName");
   //std::string graphName = getPreCutString1("TriggerEfficiencyGraphName");
   // this needs to match the path we use for data below
   //TriggerEfficiency triggerEfficiency(trigEffFileName, "hEff_Ele27");
   //TriggerEfficiency triggerEfficiency(trigEffFileName, "hEff_Ele27OR115");
   TriggerEfficiency triggerEfficiency(trigEffFileName, "hEff_Ele27OR115OR175");

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------

   CreateUserTH1D( "ProcessID"             ,    21  , -0.5    , 20.5     );
   CreateUserTH1D( "ProcessID_PAS"         ,    21  , -0.5    , 20.5     );
   CreateUserTH1D( "ProcessID_WWindow"     ,    21  , -0.5    , 20.5     );
   

   //CreateUserTH1D( "split_PAS"                ,    38    , -0.5    , 37.5     );
   //CreateUserTH1D( "split_1fb_PAS"            ,    38    , -0.5    , 37.5     );

   CreateUserTH1D( "MET_PDF"                  , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PDF"                   , 500 , 0       , 5000	 ); 
   CreateUserTH1D( "MT_PDF"                   , 400 , 0       , 4000	 ); 
   CreateUserTH1D( "Mej_PDF"                  , 400 , 0       , 4000	 ); 

   CreateUserTH1D( "nElectron_PAS"            , 5   , -0.5    , 4.5      );
   CreateUserTH1D( "nMuon_PAS"                , 5   , -0.5    , 4.5      );
   CreateUserTH1D( "nJet_PAS"                 , 11  , -0.5    , 10.5     );
   CreateUserTH1D( "nJet_PASandFrancesco"     , 11  , -0.5    , 10.5     );
   CreateUserTH1D( "Pt1stEle_PAS"	      , 200 , 0       , 2000     ); 
   CreateUserTH1D( "Energy1stEle_PAS"	      , 200 , 0       , 2000     ); 
   CreateUserTH1D( "Eta1stEle_PAS"	      , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stEle_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Charge1stEle_PAS"	      , 2   , -1.0001 , 1.0001	 ); 
   CreateUserTH1D( "MET_PAS"                  , 600 , 0       , 3000	 ); 
   CreateUserTH1D( "Pt1stEle_Barrel_PAS"	      , 200 , 0       , 2000     ); 
   CreateUserTH1D( "Pt1stEle_Endcap1_PAS"	      , 200 , 0       , 2000     ); 
   CreateUserTH1D( "Pt1stEle_Endcap2_PAS"	      , 200 , 0       , 2000     ); 
   CreateUserTH1D( "Mee_allElectrons_Presel"                  , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mee_allElectrons_3EleEvents_Presel"                  , 200 , 0       , 2000	 ); 

   CreateUserTH1D( "METPhi_PAS"		      , 60  , -3.1416 , +3.1416	 ); 
   //CreateUserTH1D( "MET_Type1_PAS"           , 200 , 0       , 1000	 ); 
   //CreateUserTH1D( "MET_Type1_Phi_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "minMETPt1stEle_PAS"       , 200 , 0       , 1000	 ); 
   //CreateUserTH1D( "minMET1Pt1stEle_PAS"     , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "Pt1stJet_PAS"             , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Pt2ndJet_PAS"             , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Eta1stJet_PAS"            , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Eta2ndJet_PAS"            , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stJet_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Phi2ndJet_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   // muon kinematics
   CreateUserTH1D( "Pt1stMuon_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta1stMuon_PAS"	      , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stMuon_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Pt2ndMuon_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta2ndMuon_PAS"	      , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi2ndMuon_PAS"	      , 60  , -3.1416 , +3.1416	 ); 

   CreateUserTH1D( "Mass1stJet_PAS"           , 100 , 0       , 500      );
   CreateUserTH1D( "Mass2ndJet_PAS"           , 100 , 0       , 500      );
   CreateUserTH1D( "CISV1stJet_PAS"            , 200 , 0       , 1.0	 ); 
   CreateUserTH1D( "CISV2ndJet_PAS"            , 200 , 0       , 1.0	 ); 
   CreateUserTH1D( "nMuon_PtCut_IDISO_PAS"    , 16  , -0.5    , 15.5	 ); 
   CreateUserTH1D( "MTenu_PAS"                , 400 , 0       , 2000	 ); 
   CreateUserTH1D( "Ptenu_PAS"		      , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTlep_PAS"                , 200 , 0       , 2000	 ); 
   //CreateUserTH1D( "sTlep_Type1_PAS"         , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTjet_PAS"                , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PAS"                   , 300 , 0       , 3000	 ); 
   CreateUserTH1D( "sT_ControlRegion_Njet_gte4", 300 , 0       , 3000	 ); 
   CreateUserTH1D( "sT_ControlRegion_Njet_lte3", 300 , 0       , 3000	 ); 
   //CreateUserTH1D( "sT_Type1_PAS"            , 300 , 0       , 3000	 ); 
   CreateUserTH1D( "Mjj_PAS"		      , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mej1_PAS"                 , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mej2_PAS"                 , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mej_PAS"                  , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mejj_PAS"                 , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "MTjnu_PAS"                , 200 , 0       , 1000     );
   CreateUserTH1D( "MTjjnu_PAS"               , 200 , 0       , 1000     );
   CreateUserTH1D( "DCotTheta1stEle_PAS"      , 100 , 0.0     , 1.0      );
   CreateUserTH1D( "Dist1stEle_PAS"           , 100 , 0.0     , 1.0      );  
   CreateUserTH1D( "DR_Ele1Jet1_PAS"	      , 100 , 0       , 10       ); 
   CreateUserTH1D( "DR_Ele1Jet2_PAS"	      , 100 , 0       , 10       ); 
   CreateUserTH1D( "DR_Jet1Jet2_PAS"	      , 100 , 0       , 10       );
   CreateUserTH1D( "minDR_EleJet_PAS"         , 100 , 0       , 10       );
   CreateUserTH1D( "mDPhiJet1Jet2_PAS"        , 100 , 0.      , 3.14159  );
   CreateUserTH1D( "mDPhi1stEleMET_PAS"       , 100 , 0.      , 3.14159  );
   CreateUserTH1D( "mDPhi1stJetMET_PAS"       , 100 , 0.      , 3.14159  );
   CreateUserTH1D( "mDPhi2ndJetMET_PAS"       , 100 , 0.      , 3.14159  );
   CreateUserTH1D( "GeneratorWeight"          , 200 , -2.0    , 2.0      );
   CreateUserTH1D( "PileupWeight"             , 200 , -2.0    , 2.0      );
   CreateUserTH1D( "EOverP_1stEle_PAS"        , 20000 , 0.      , 100.     );

   CreateUserTH1D( "Pt1stEle_Pt40to45_EtaGT2p1", 150, 35., 50. );

   CreateUserTH1D( "MT_GoodVtxLTE3_PAS"       , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MT_GoodVtxGTE4_LTE8_PAS"  , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MT_GoodVtxGTE9_LTE15_PAS" , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MT_GoodVtxGTE16_PAS"      , 200 , 0.      ,  1000    );

   CreateUserTH1D( "nVertex_PAS"           ,    101   , -0.5   , 100.5	 ) ; 
   
   // MT plots
   CreateUserTH1D( "MTenu_50_110", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte3", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte5", 240, 40, 160 );

   CreateUserTH1D( "MTenu_70_110", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_70_110"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte3", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte5", 240, 40, 160 );

   CreateUserTH1D( "MTenu_70_150", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_70_150"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_70_150_Njet_gte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_lte3", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_lte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_gte5", 240, 40, 160 );

   CreateUserTH1D( "MTenu_110_190", 200, 100, 200 );
   CreateUserTH1D( "nJets_MTenu_110_190"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_110_190_Njet_gte4", 200, 100, 200 );
   CreateUserTH1D( "MTenu_110_190_Njet_lte3", 200, 100, 200 );
   CreateUserTH1D( "MTenu_110_190_Njet_lte4", 200, 100, 200 );
   CreateUserTH1D( "MTenu_110_190_Njet_gte5", 200, 100, 200 );

   // without btags
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_btagSFUpShift", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_btagSFDownShift", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_noBtaggedJets"    , 20 , -0.5, 19.5 );
   // for scale factor dependence studies
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_sT300To500_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_sT500To750_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_sT750To1250_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_sT1250ToInf_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_Mej100To200_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_Mej200To300_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_Mej300To400_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_Mej400To500_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_Mej500To650_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_Mej650ToInf_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_noBtaggedJets_MT200To400_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_noBtaggedJets_MT400To600_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_noBtaggedJets_MT600To900_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_noBtaggedJets_MT900ToInf_PAS"		             ,    200   , 0       , 2000	  ); 
   //
   CreateUserTH1D( "MTenu_50_110_Njet_gte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte3_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte5_noBtaggedJets", 240, 40, 160 );

   CreateUserTH1D( "MTenu_70_110_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_70_110_noBtaggedJets"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte3_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte5_noBtaggedJets", 240, 40, 160 );

   CreateUserTH1D( "MTenu_70_150_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_70_150_noBtaggedJets"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_70_150_Njet_gte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_lte3_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_lte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_gte5_noBtaggedJets", 240, 40, 160 );

   CreateUserTH1D( "MTenu_110_190_noBtaggedJets", 200, 100, 200 );
   CreateUserTH1D( "nJets_MTenu_110_190_noBtaggedJets"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_110_190_Njet_gte4_noBtaggedJets", 200, 100, 200 );
   CreateUserTH1D( "MTenu_110_190_Njet_lte3_noBtaggedJets", 200, 100, 200 );
   CreateUserTH1D( "MTenu_110_190_Njet_lte4_noBtaggedJets", 200, 100, 200 );
   CreateUserTH1D( "MTenu_110_190_Njet_gte5_noBtaggedJets", 200, 100, 200 );

   // with a btag
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_btagSFUpShift", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_btagSFDownShift", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet"    , 20 , -0.5, 19.5 );
   // for scale factor dependence studies
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_sT300To500_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_sT500To750_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_sT750To1250_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_sT1250ToInf_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_Mej100To200_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_Mej200To300_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_Mej300To400_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_Mej400To500_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_Mej500To650_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_Mej650ToInf_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_gteOneBtaggedJet_MT200To400_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_gteOneBtaggedJet_MT400To600_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_gteOneBtaggedJet_MT600To900_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_gteOneBtaggedJet_MT900ToInf_PAS"		             ,    200   , 0       , 2000	  ); 
   //
   CreateUserTH1D( "MTenu_50_110_Njet_gte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte3_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte5_gteOneBtaggedJet", 240, 40, 160 );

   CreateUserTH1D( "MTenu_70_110_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_70_110_gteOneBtaggedJet"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte3_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte5_gteOneBtaggedJet", 240, 40, 160 );

   CreateUserTH1D( "MTenu_70_150_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_70_150_gteOneBtaggedJet"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_70_150_Njet_gte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_lte3_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_lte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_gte5_gteOneBtaggedJet", 240, 40, 160 );

   CreateUserTH1D( "MTenu_110_190_gteOneBtaggedJet", 200, 100, 200 );
   CreateUserTH1D( "nJets_MTenu_110_190_gteOneBtaggedJet"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_110_190_Njet_gte4_gteOneBtaggedJet", 200, 100, 200 );
   CreateUserTH1D( "MTenu_110_190_Njet_lte3_gteOneBtaggedJet", 200, 100, 200 );
   CreateUserTH1D( "MTenu_110_190_Njet_lte4_gteOneBtaggedJet", 200, 100, 200 );
   CreateUserTH1D( "MTenu_110_190_Njet_gte5_gteOneBtaggedJet", 200, 100, 200 );

   // with M(jj)~M(W)
   CreateUserTH1D( "MTenu_50_110_Mjj50to110", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_Mjj50to110"    , 20 , -0.5, 19.5 );
   // with M(jj)~M(W) and at least one additional b-tagged jet
   CreateUserTH1D( "MTenu_50_110_Mjj50to110_addBtagJet", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_Mjj50to110_addBtagJet"    , 20 , -0.5, 19.5 );
   // with M(jj)~M(W) and no additional b-tagged jets
   CreateUserTH1D( "MTenu_50_110_Mjj50to110_noAddBtagJets", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_Mjj50to110_noAddBtagJets"    , 20 , -0.5, 19.5 );
   // with M(jj) !~ M(W)
   CreateUserTH1D( "MTenu_50_110_MjjGte110", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_MjjGte110"    , 20 , -0.5, 19.5 );


   CreateUserTH1D( "run_PAS"               ,    20000 , 160300  , 180300 );
   CreateUserTH1D( "run_HLT"               ,    20000 , 160300  , 180300 );

   CreateUserTH1D( "Eta1stJet_PASand2Jet"  , 100 , -5 , 5 ); 
   CreateUserTH1D( "Eta1stJet_PASand3Jet"  , 100 , -5 , 5 ); 
   CreateUserTH1D( "Eta1stJet_PASand4Jet"  , 100 , -5 , 5 ); 
   CreateUserTH1D( "Eta1stJet_PASand5Jet"  , 100 , -5 , 5 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_50_110_Njet_lte3", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_50_110_Njet_lte3", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_50_110_Njet_lte3", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_50_110_Njet_lte3"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_50_110_Njet_lte3"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_50_110_Njet_lte3"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_50_110_Njet_lte4", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_50_110_Njet_lte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_50_110_Njet_lte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_50_110_Njet_lte4"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_50_110_Njet_lte4"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_50_110_Njet_lte4"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_50_110_Njet_gte4", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_50_110_Njet_gte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_50_110_Njet_gte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_50_110_Njet_gte4"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_50_110_Njet_gte4"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_50_110_Njet_gte4"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_50_110_Njet_gte5", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_50_110_Njet_gte5", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_50_110_Njet_gte5", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_50_110_Njet_gte5"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_50_110_Njet_gte5"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_50_110_Njet_gte5"     , 200 , 0 , 2000 ); 
   
   // 70-150
   CreateUserTH1D( "Pt1stEle_MTenu_70_150_Njet_lte3", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_150_Njet_lte3", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_150_Njet_lte3", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_150_Njet_lte3"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_150_Njet_lte3"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_150_Njet_lte3"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_150_Njet_lte4", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_150_Njet_lte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_150_Njet_lte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_150_Njet_lte4"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_150_Njet_lte4"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_150_Njet_lte4"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_150_Njet_gte4", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_150_Njet_gte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_150_Njet_gte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_150_Njet_gte4"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_150_Njet_gte4"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_150_Njet_gte4"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_150_Njet_gte5", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_150_Njet_gte5", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_150_Njet_gte5", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_150_Njet_gte5"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_150_Njet_gte5"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_150_Njet_gte5"     , 200 , 0 , 2000 ); 
   
   // 110-190
   CreateUserTH1D( "Pt1stEle_MTenu_110_190_Njet_lte3", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_110_190_Njet_lte3", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_110_190_Njet_lte3", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_110_190_Njet_lte3"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_110_190_Njet_lte3"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_110_190_Njet_lte3"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_110_190_Njet_lte4", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_110_190_Njet_lte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_110_190_Njet_lte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_110_190_Njet_lte4"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_110_190_Njet_lte4"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_110_190_Njet_lte4"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_110_190_Njet_gte4", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_110_190_Njet_gte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_110_190_Njet_gte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_110_190_Njet_gte4"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_110_190_Njet_gte4"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_110_190_Njet_gte4"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_110_190_Njet_gte5", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_110_190_Njet_gte5", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_110_190_Njet_gte5", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_110_190_Njet_gte5"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_110_190_Njet_gte5"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_110_190_Njet_gte5"     , 200 , 0 , 2000 ); 

   // no btags
   CreateUserTH1D( "Pt1stEle_MTenu_50_110_Njet_lte3_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_50_110_Njet_lte3_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_50_110_Njet_lte3_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_50_110_Njet_lte3_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_50_110_Njet_lte3_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_50_110_Njet_lte3_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_50_110_Njet_lte4_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_50_110_Njet_lte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_50_110_Njet_lte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_50_110_Njet_lte4_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_50_110_Njet_lte4_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_50_110_Njet_lte4_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_50_110_Njet_gte4_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_50_110_Njet_gte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_50_110_Njet_gte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_50_110_Njet_gte4_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_50_110_Njet_gte4_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_50_110_Njet_gte4_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_50_110_Njet_gte5_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_50_110_Njet_gte5_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_50_110_Njet_gte5_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_50_110_Njet_gte5_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_50_110_Njet_gte5_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_50_110_Njet_gte5_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   // 70-150
   CreateUserTH1D( "Pt1stEle_MTenu_70_150_Njet_lte3_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_150_Njet_lte3_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_150_Njet_lte3_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_150_Njet_lte3_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_150_Njet_lte3_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_150_Njet_lte3_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_150_Njet_lte4_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_150_Njet_lte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_150_Njet_lte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_150_Njet_lte4_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_150_Njet_lte4_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_150_Njet_lte4_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_150_Njet_gte4_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_150_Njet_gte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_150_Njet_gte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_150_Njet_gte4_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_150_Njet_gte4_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_150_Njet_gte4_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_150_Njet_gte5_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_150_Njet_gte5_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_150_Njet_gte5_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_150_Njet_gte5_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_150_Njet_gte5_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_150_Njet_gte5_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   // 110-190
   CreateUserTH1D( "Pt1stEle_MTenu_110_190_Njet_lte3_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_110_190_Njet_lte3_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_110_190_Njet_lte3_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_110_190_Njet_lte3_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_110_190_Njet_lte3_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_110_190_Njet_lte3_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_110_190_Njet_lte4_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_110_190_Njet_lte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_110_190_Njet_lte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_110_190_Njet_lte4_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_110_190_Njet_lte4_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_110_190_Njet_lte4_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_110_190_Njet_gte4_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_110_190_Njet_gte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_110_190_Njet_gte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_110_190_Njet_gte4_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_110_190_Njet_gte4_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_110_190_Njet_gte4_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_110_190_Njet_gte5_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_110_190_Njet_gte5_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_110_190_Njet_gte5_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_110_190_Njet_gte5_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_110_190_Njet_gte5_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_110_190_Njet_gte5_noBtaggedJets"     , 200 , 0 , 2000 ); 

   // at least one btag
   CreateUserTH1D( "Pt1stEle_MTenu_50_110_Njet_lte3_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_50_110_Njet_lte3_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_50_110_Njet_lte3_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_50_110_Njet_lte3_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_50_110_Njet_lte3_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_50_110_Njet_lte3_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_50_110_Njet_lte4_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_50_110_Njet_lte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_50_110_Njet_lte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_50_110_Njet_lte4_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_50_110_Njet_lte4_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_50_110_Njet_lte4_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_50_110_Njet_gte4_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_50_110_Njet_gte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_50_110_Njet_gte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_50_110_Njet_gte4_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_50_110_Njet_gte4_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_50_110_Njet_gte4_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_50_110_Njet_gte5_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_50_110_Njet_gte5_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_50_110_Njet_gte5_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_50_110_Njet_gte5_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_50_110_Njet_gte5_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_50_110_Njet_gte5_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   // 70-150
   CreateUserTH1D( "Pt1stEle_MTenu_70_150_Njet_lte3_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_150_Njet_lte3_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_150_Njet_lte3_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_150_Njet_lte3_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_150_Njet_lte3_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_150_Njet_lte3_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_150_Njet_lte4_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_150_Njet_lte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_150_Njet_lte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_150_Njet_lte4_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_150_Njet_lte4_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_150_Njet_lte4_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_150_Njet_gte4_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_150_Njet_gte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_150_Njet_gte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_150_Njet_gte4_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_150_Njet_gte4_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_150_Njet_gte4_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_150_Njet_gte5_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_150_Njet_gte5_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_150_Njet_gte5_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_150_Njet_gte5_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_150_Njet_gte5_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_150_Njet_gte5_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   // 110-190
   CreateUserTH1D( "Pt1stEle_MTenu_110_190_Njet_lte3_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_110_190_Njet_lte3_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_110_190_Njet_lte3_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_110_190_Njet_lte3_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_110_190_Njet_lte3_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_110_190_Njet_lte3_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_110_190_Njet_lte4_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_110_190_Njet_lte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_110_190_Njet_lte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_110_190_Njet_lte4_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_110_190_Njet_lte4_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_110_190_Njet_lte4_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_110_190_Njet_gte4_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_110_190_Njet_gte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_110_190_Njet_gte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_110_190_Njet_gte4_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_110_190_Njet_gte4_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_110_190_Njet_gte4_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_110_190_Njet_gte5_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_110_190_Njet_gte5_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_110_190_Njet_gte5_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_110_190_Njet_gte5_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_110_190_Njet_gte5_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_110_190_Njet_gte5_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 

   // Ele vars
   CreateUserTH1D("DeltaPhi1stEle_Presel", 200, -0.1, 0.1);
   CreateUserTH1D("SCEta1stEle_Presel", 100, -5, 5);
   CreateUserTH1D("PtHeep1stEle_Presel",200,0,2000);
   CreateUserTH1D("SCEt1stEle_Presel",200,0,2000);
   CreateUserTH1D("HoE1stEle_Presel",400,0,0.1);
   CreateUserTH1D("EcalIso1stEle_Presel",200,0,200);
   CreateUserTH1D("HcalIso1stEle_Presel",200,0,200);
   CreateUserTH1D("TrkIsoHEEP701stEle_Presel",200,0,200);
   CreateUserTH1D("LeadVtxDistXY1stEle_Presel", 200, -0.05,   0.05 );
   CreateUserTH1D("DeltaEtaEleTrk1stEle_Presel", 400, -0.5,   0.5 );

   // lower Pt regions
   CreateUserTH1D("DeltaPhi1stEle_Pt50to100", 200, -0.1, 0.1);
   CreateUserTH1D("SCEta1stEle_Pt50to100", 100, -5, 5);
   CreateUserTH1D("PtHeep1stEle_Pt50to100",200,0,2000);
   CreateUserTH1D("HoE1stEle_Pt50to100",400,0,0.1);
   CreateUserTH1D("EcalIso1stEle_Pt50to100",200,0,200);
   CreateUserTH1D("HcalIso1stEle_Pt50to100",200,0,200);
   CreateUserTH1D("TrkIsoHEEP701stEle_Pt50to100",200,0,200);
   CreateUserTH1D("LeadVtxDistXY1stEle_Pt50to100", 200, -0.05,   0.05 );
   //
   CreateUserTH1D("DeltaPhi1stEle_Pt130to200", 200, -0.1, 0.1);
   CreateUserTH1D("SCEta1stEle_Pt130to200", 100, -5, 5);
   CreateUserTH1D("PtHeep1stEle_Pt130to200",200,0,2000);
   CreateUserTH1D("HoE1stEle_Pt130to200",400,0,0.1);
   CreateUserTH1D("EcalIso1stEle_Pt130to200",200,0,200);
   CreateUserTH1D("HcalIso1stEle_Pt130to200",200,0,200);
   CreateUserTH1D("TrkIsoHEEP701stEle_Pt130to200",200,0,200);
   CreateUserTH1D("LeadVtxDistXY1stEle_Pt130to200", 200, -0.05,   0.05 );

   CreateUserTH1D( "Mej_Barrel_Presel"                  , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mej_Endcap1_Presel"                  , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mej_Endcap2_Presel"                  , 200 , 0       , 2000	 ); 
   //

   char plot_name[200];
   
   //for (int i_lq_mass = 0; i_lq_mass < n_lq_mass ; ++i_lq_mass ) { 
   //  int lq_mass = LQ_MASS[i_lq_mass];
   //  sprintf(plot_name, "likelihood_LQ%d"   , lq_mass ); CreateUserTH1D ( plot_name, 200, 0.0, 1.0 );
   //  sprintf(plot_name, "logLikelihood_LQ%d", lq_mass ); CreateUserTH1D ( plot_name, 100, 0.0, 100.0 );
   //}
   
   //--------------------------------------------------------------------------
   // Final selection plots
   //--------------------------------------------------------------------------
   bool doFinalSelections = false;
   // check if there is a final Mej specific in cutfile for any LQ mass
   for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
     int lq_mass = LQ_MASS[i_lq_mass];
     sprintf(cut_name, "Mej_LQ%d"   , lq_mass );
     if(hasCut(cut_name)) {
       doFinalSelections = true;
       break;
     }
   }
   // now, we must have an Mej cut and optimization must be off to have final selections enabled
   doFinalSelections = doFinalSelections && !isOptimizationEnabled();
   if ( doFinalSelections ) { 
     for (int i_lq_mass = 0; i_lq_mass < n_lq_mass ; ++i_lq_mass ) { 
       int lq_mass = LQ_MASS[i_lq_mass];
       //sprintf(plot_name, "split_1fb_LQ%d", lq_mass ); CreateUserTH1D ( plot_name, 21  , -0.5, 20.5);
       sprintf(plot_name, "MET_LQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 600 , 0 , 3000 ); 
       //sprintf(plot_name, "MET_UnclUp_LQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 
       //sprintf(plot_name, "MET_UnclDown_LQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 
       sprintf(plot_name, "Mej_LQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "Mejj_LQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "Mjj_LQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "Mej3_LQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "ST_LQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "MTenu_LQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 );
       sprintf(plot_name, "MTjnu_LQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MTjjnu_LQ%d"   , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MTejjnu_LQ%d"  , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "nJets_LQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 21  , -0.5, 20.5 );
       sprintf(plot_name, "mDPhi_Jet1Jet2_LQ%d", lq_mass ); CreateUserTH1D( plot_name, 100 , 0.      , 3.14159  );
       sprintf(plot_name, "mDPhi_Ele1MET_LQ%d" , lq_mass ); CreateUserTH1D( plot_name, 100 , 0.      , 3.14159  );
       sprintf(plot_name, "nBJet_loose_LQ%d" , lq_mass ); CreateUserTH1D( plot_name ,     10 , -0.5    , 9.5      );
       sprintf(plot_name, "nBJet_medium_LQ%d", lq_mass ); CreateUserTH1D( plot_name ,     10 , -0.5    , 9.5      );
       sprintf(plot_name, "nBJet_tight_LQ%d" , lq_mass ); CreateUserTH1D( plot_name ,     10 , -0.5    , 9.5      );

       sprintf(plot_name, "Ele_EtaVsPhi_LQ%d" , lq_mass ); CreateUserTH2D( plot_name ,    100 , -2.5,   2.5, 100, -3.15, 3.15 );
       sprintf(plot_name, "Jet_EtaVsPhi_LQ%d" , lq_mass ); CreateUserTH2D( plot_name ,    100 , -2.5,   2.5, 100, -3.15, 3.15 );
       sprintf(plot_name, "Jet1_EtaVsPhi_LQ%d", lq_mass ); CreateUserTH2D( plot_name ,    100 , -2.5,   2.5, 100, -3.15, 3.15 );
       sprintf(plot_name, "Jet2_EtaVsPhi_LQ%d", lq_mass ); CreateUserTH2D( plot_name ,    100 , -2.5,   2.5, 100, -3.15, 3.15 );

       sprintf(plot_name, "ST_StLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_StLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_StLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 );
       sprintf(plot_name, "MET_StLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_MtLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_MtLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_MtLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 );
       sprintf(plot_name, "MET_MtLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_MejLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_MejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_MejLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 );
       sprintf(plot_name, "MET_MejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_MetLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_MetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_MetLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 );
       sprintf(plot_name, "MET_MetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_stAndMtAndMetLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_stAndMtAndMetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_stAndMtAndMetLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 );
       sprintf(plot_name, "MET_stAndMtAndMetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_stAndMtAndMejLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_stAndMtAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_stAndMtAndMejLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 );
       sprintf(plot_name, "MET_stAndMtAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_stAndMetAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_stAndMetAndMejLQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_stAndMetAndMejLQ%d"   , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 );
       sprintf(plot_name, "MET_stAndMetAndMejLQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_mtAndMetAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_mtAndMetAndMejLQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_mtAndMetAndMejLQ%d"   , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 );
       sprintf(plot_name, "MET_mtAndMetAndMejLQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 
       
       sprintf(plot_name, "ST_StAndMejLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_StAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_StAndMejLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 );
       sprintf(plot_name, "MET_StAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_MtAndMejLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_MtAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_MtAndMejLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 );
       sprintf(plot_name, "MET_MtAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_MetAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_MetAndMejLQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_MetAndMejLQ%d"   , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 );
       sprintf(plot_name, "MET_MetAndMejLQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_StAndMetLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_StAndMetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_StAndMetLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 );
       sprintf(plot_name, "MET_StAndMetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_MtAndMetLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_MtAndMetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_MtAndMetLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 );
       sprintf(plot_name, "MET_MtAndMetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_StAndMtLQ%d"        , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_StAndMtLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_StAndMtLQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 );
       sprintf(plot_name, "MET_StAndMtLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "Mej_Barrel_LQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000	); 
       sprintf(plot_name, "Mej_Endcap1_LQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000	); 
       sprintf(plot_name, "Mej_Endcap2_LQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000	); 

       if ( do_extra_finalSelection_plots ) { 
         sprintf(plot_name, "sTfrac_Jet1_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
         sprintf(plot_name, "sTfrac_Jet2_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
         sprintf(plot_name, "sTfrac_Ele1_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
         sprintf(plot_name, "sTfrac_MET_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
         sprintf(plot_name, "sTfrac_Jet_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
         sprintf(plot_name, "sTfrac_Lep_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
         sprintf(plot_name, "Charge1stEle_LQ%d"	       , lq_mass ); CreateUserTH1D( plot_name ,      3 ,  -1.5   , 1.5      );
         sprintf(plot_name, "Energy1stEle_LQ%d"	       , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
         sprintf(plot_name, "Pt1stEle_LQ%d"	       , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
         sprintf(plot_name, "PtHeep1stEle_LQ%d"	       , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
         sprintf(plot_name, "SCEt1stEle_LQ%d"	       , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
         sprintf(plot_name, "Eta1stEle_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
         sprintf(plot_name, "SCEta1stEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
         sprintf(plot_name, "Phi1stEle_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , +3.1416  ); 
         sprintf(plot_name, "Pt1stJet_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
         sprintf(plot_name, "Pt2ndJet_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
         sprintf(plot_name, "Eta1stJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
         sprintf(plot_name, "Eta2ndJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
         sprintf(plot_name, "Phi1stJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , 3.1416   ); 
         sprintf(plot_name, "Phi2ndJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , 3.1416   ); 
         sprintf(plot_name, "sTlep_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
         sprintf(plot_name, "sTjet_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
         sprintf(plot_name, "Me1j1_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
         sprintf(plot_name, "Me1j2_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
         sprintf(plot_name, "nVertex_LQ%d"             , lq_mass ); CreateUserTH1D( plot_name ,    101 , -0.5    ,  100.5   ); 
         sprintf(plot_name, "MeeVsST_LQ%d"             , lq_mass ); CreateUserTH2D( plot_name ,    200 ,  0.0, 2000, 200, 0, 2000);
         // muon kinematics
         sprintf(plot_name, "Pt1stMuon_LQ%d"	       , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
         sprintf(plot_name, "Eta1stMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
         sprintf(plot_name, "Phi1stMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , +3.1416  ); 
         sprintf(plot_name, "Pt2ndMuon_LQ%d"	       , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
         sprintf(plot_name, "Eta2ndMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
         sprintf(plot_name, "Phi2ndMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , +3.1416  ); 
       }

       if ( do_eleQuality_plots ) { 
         sprintf(plot_name, "BeamSpotDXY_1stEle_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   0.5  );
         sprintf(plot_name, "Classif_1stEle_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name , 5  , -0.5 ,   4.5  );
         sprintf(plot_name, "CorrIsolation_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200,-25.0 ,  25.0  );
         sprintf(plot_name, "DeltaEtaTrkSC_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.01,   0.01 );
         sprintf(plot_name, "DeltaPhiTrkSC_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.1 ,   0.1  );
         sprintf(plot_name, "Full5x5E1x5OverE5x5_1stEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );
         sprintf(plot_name, "Full5x5E2x5OverE5x5_1stEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );
         sprintf(plot_name, "EcalIsolation_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,  20.0  );
         sprintf(plot_name, "HcalIsolation_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,  20.0  );
         sprintf(plot_name, "TrkIsolation_1stEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    5.0  );
         sprintf(plot_name, "Energy_1stEle_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,3000.0  );
         sprintf(plot_name, "FBrem_1stEle_LQ%d"                , lq_mass ); CreateUserTH1D( plot_name , 200,-10.0 ,  10.0  );
         sprintf(plot_name, "GsfCtfCharge_1stEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
         sprintf(plot_name, "GsfCtfScPixCharge_1stEle_LQ%d"    , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
         sprintf(plot_name, "GsfScPixCharge_1stEle_LQ%d"       , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
         sprintf(plot_name, "HasMatchedPhot_1stEle_LQ%d"       , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
         sprintf(plot_name, "HoE_1stEle_LQ%d"                  , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   0.05 );
         sprintf(plot_name, "LeadVtxDistXY_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.05,   0.05 );
         sprintf(plot_name, "LeadVtxDistZ_1stEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 200, -0.2 ,   0.2  );
         sprintf(plot_name, "MissingHits_1stEle_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name , 2  , -0.5,    1.5  );
         sprintf(plot_name, "NBrems_1stEle_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name , 11 , -0.5,   10.5  );
         sprintf(plot_name, "EnergyORawEnergy_1stEle_LQ%d"     , lq_mass ); CreateUserTH1D( plot_name , 200,  0.9,    1.4  );
         sprintf(plot_name, "SigmaEtaEta_Barrel_1stEle_LQ%d"   , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.02 );
         sprintf(plot_name, "SigmaEtaEta_Endcap_1stEle_LQ%d"   , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.1  );
         sprintf(plot_name, "Full5x5SigmaIEtaIEta_Barrel_1stEle_LQ%d" , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.02 );
         sprintf(plot_name, "Full5x5SigmaIEtaIEta_Endcap_1stEle_LQ%d" , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.1  );
         sprintf(plot_name, "TrkPtOPt_1stEle_LQ%d"             , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,  100.0  );
         sprintf(plot_name, "ValidFrac_1stEle_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );
         sprintf(plot_name, "EOverP_1stEle_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name , 20000,  0.0 , 100.0  );
       }
     }
   }
   
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_LQ300", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ300"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ300", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ300"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_LQ400", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ400"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ400", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ400"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_LQ500", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ500"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ500", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ500"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_LQ600", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ600"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ600", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ600"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_LQ700", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ700"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ700", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ700"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_LQ800", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ800"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ800", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ800"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_LQ900", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ900"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ900", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ900"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets_LQ1000", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ1000"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ1000", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ1000"    , 20 , -0.5, 19.5 );
   //
   CreateUserTH1D( "MET_LQ300_NoMETCut", 200 , 0 , 1000 );
   CreateUserTH1D( "Mej_LQ300_NoMETCut", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MT_LQ300_NoMETCut" , 200 , 0 , 1000 ); 
   CreateUserTH1D( "ST_LQ300_NoMETCut" , 300 , 0 , 3000 ); 

   //--------------------------------------------------------------------------
   // Loop over the chain
   //--------------------------------------------------------------------------

   if (fChain == 0) return;
   //std::cout << "fChain = " << fChain << std::endl;

   Long64_t nentries = fChain->GetEntries();
   // Long64_t nentries = 1;
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0)
     {
       std::cout << "ERROR: Could not read from TTree; exiting." << std::endl;
       exit(-1);
     }
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if (nb < 0)
     {
       std::cout << "ERROR: Could not read entry from TTree: read " << nb << "bytes; exiting." << std::endl;
       exit(-2);
     }

     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << "/" << nentries << std::endl;   

     //--------------------------------------------------------------------------
     // Reset the cuts
     //--------------------------------------------------------------------------

     resetCuts();

     //--------------------------------------------------------------------------
     // Check good run list
     //--------------------------------------------------------------------------
     
     int passedJSON = passJSON ( run, ls , isData ) ;

     //--------------------------------------------------------------------------
     // Do pileup re-weighting
     //--------------------------------------------------------------------------
     
     double pileup_weight = getPileupWeight ( nPileUpInt_True, isData ) ;
     
     //--------------------------------------------------------------------------
     // Get information about gen-level reweighting (should be for Sherpa only)
     //--------------------------------------------------------------------------

     double gen_weight = Weight;
     if ( isData ) gen_weight = 1.0;

     std::string current_file_name ( fChain->GetCurrentFile()->GetName());

     //// TopPt reweight
     //// only valid for powheg
     //if(current_file_name.find("TT_") != std::string::npos) {
     //  gen_weight*=TopPtWeight;
     //}
     // removed 2018 Mar. 2

     // Electron scale factors for MC only
     if(!isData)
     {
       float recoSFEle1 = ElectronScaleFactors2016::LookupRecoSF(Ele1_SCEta);
       float heepSFEle1 = ElectronScaleFactors2016::LookupHeepSF(Ele1_SCEta);
       float totalScaleFactor = recoSFEle1*heepSFEle1;
       gen_weight*=totalScaleFactor;
     }

     //if ( Ele2_ValidFrac > 998. && Ele1_ValidFrac > 998. ) nMuon_ptCut = 0;

     //--------------------------------------------------------------------------
     // Is this a barrel electron?
     //--------------------------------------------------------------------------

     bool Ele1_IsBarrel = bool ( fabs (Ele1_Eta) < 1.442 );

     //--------------------------------------------------------------------------
     // First variable to fill just shows the "reweighting".  Always passes.
     //--------------------------------------------------------------------------

     fillVariableWithValue ( "Reweighting", 1, gen_weight * pileup_weight  );

     ////--------------------------------------------------------------------------
     //// Special treatment of inclusive W/Z
     ////--------------------------------------------------------------------------
     //bool passGenWZPt = true;
     //// inclusive
     //if(current_file_name.find("WJetsToLNu_ext1_amcatnloFXFX") != std::string::npos 
     //    || current_file_name.find("WJetsToLNu_amcatnloFXFX") != std::string::npos) {
     //  if(GenW1_Pt > 120) passGenWZPt = false; // if W Pt > 120 GeV, cut it out
     //}
     //if(current_file_name.find("DYJetsToLL_M-50") != std::string::npos) {
     //  if(GenZGamma1_Pt > 70) passGenWZPt = false; // if Z/gamma Pt > 70 GeV, cut it out
     //}
     //// first pt bin
     //if(current_file_name.find("WJetsToLNu_Pt-100") != std::string::npos) {
     //  if(GenW1_Pt <= 120) passGenWZPt = false;
     //}
     //if(current_file_name.find("DYJetsToLL_Pt-50") != std::string::npos) {
     //  if(GenZGamma1_Pt <= 70) passGenWZPt = false;
     //}
     ////// testing
     ////if(current_file_name.find("WJetsToLNu_ext1_amcatnloFXFX") != std::string::npos 
     ////    || current_file_name.find("WJetsToLNu_amcatnloFXFX") != std::string::npos) {
     ////  if(GenW1_Pt <= 100) passGenWZPt = false; // if W Pt > 100 GeV, cut it out
     ////}
     ////if(current_file_name.find("DYJetsToLL_M-50_amcatnloFXFX") != std::string::npos) {
     ////  if(GenZGamma1_Pt <= 100) passGenWZPt = false; // if Z/gamma Pt > 100 GeV, cut it out
     ////}
     //fillVariableWithValue("PassGenWZPt",passGenWZPt,gen_weight*pileup_weight);

     //--------------------------------------------------------------------------
     // Fill JSON variable
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue ("PassJSON", passedJSON, gen_weight * pileup_weight); 

     //--------------------------------------------------------------------------
     // Fill noise filters
     //--------------------------------------------------------------------------
     // see: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
     // we filled these at skim time
     fillVariableWithValue( "PassGlobalTightHalo2016Filter" , PassGlobalTightHalo2016Filter  , gen_weight * pileup_weight );
     fillVariableWithValue( "PassGoodVertices"	            , PassGoodVertices               , gen_weight * pileup_weight );
     fillVariableWithValue( "PassHBHENoiseFilter"	          , PassHBHENoiseFilter            , gen_weight * pileup_weight );
     fillVariableWithValue( "PassHBHENoiseIsoFilter"	      , PassHBHENoiseIsoFilter         , gen_weight * pileup_weight );
     fillVariableWithValue( "PassBadEESupercrystalFilter"   , PassBadEESupercrystalFilter    , gen_weight * pileup_weight );
     fillVariableWithValue( "PassEcalDeadCellTrigPrim"      , PassEcalDeadCellTrigPrim       , gen_weight * pileup_weight );
     fillVariableWithValue( "PassChargedCandidateFilter"    , PassChargedCandidateFilter     , gen_weight * pileup_weight );
     fillVariableWithValue( "PassBadPFMuonFilter"           , PassBadPFMuonFilter            , gen_weight * pileup_weight );
     //

     
     //--------------------------------------------------------------------------
     // Fill HLT
     //--------------------------------------------------------------------------

     int passHLT = 0;
     bool verboseTrigEff = false;
     if (isData) {
       if(current_file_name.find("SinglePhoton") != std::string::npos) {
         if (H_Photon175 == 1) // take events triggered by Photon175 only plus those triggered by Photon175 AND Ele27/Ele115
           passHLT = 1;
       }
       else if(current_file_name.find("SingleElectron") != std::string::npos) {
         if (H_Photon175 != 1 && (H_Ele27_WPTight == 1 || H_Ele115_CIdVT_GsfIdT == 1) ) // take events triggered only by Ele27 OR Ele115
           passHLT = 1;
       }
     }
     else {
       passHLT = triggerEfficiency.PassTrigger(Ele1_SCEta,Ele1_SCEnergy/cosh(Ele1_SCEta),verboseTrigEff) ? 1 : 0;
       // for enujj, we only get one chance to trigger
     }

     fillVariableWithValue ( "PassHLT", passHLT, gen_weight * pileup_weight  ) ;     

     //--------------------------------------------------------------------------
     // Fill more variables
     //--------------------------------------------------------------------------
     // Muon variables ( for veto ) 					      
     fillVariableWithValue(   "nMuon"                    , nMuon_ptCut             , gen_weight * pileup_weight );
     // now check if there are no muons but gen muons passing pt/eta which should have passed ID but did not
     // correct by a factor of (1-effData)/(1-effMC) to evaluate uncertainty
     if(nMuon_ptCut==0) // no ID'ed muons in event
     {
       // check for gen muons which should have passed ID but did not
       // pt threshold is defined in highPt muon ID in MuonIDs
       if(GenMu1_Pt > 35.0 && fabs(GenMu1_Eta) <= 2.4)
         gen_weight*=MuonScaleFactors::GetVetoMCDataEffRatio(GenMu1_Eta);
       if(GenMu2_Pt > 35.0 && fabs(GenMu2_Eta) <= 2.4)
         gen_weight*=MuonScaleFactors::GetVetoMCDataEffRatio(GenMu2_Eta);
       if(GenMu3_Pt > 35.0 && fabs(GenMu3_Eta) <= 2.4)
         gen_weight*=MuonScaleFactors::GetVetoMCDataEffRatio(GenMu3_Eta);
     }
			                                      		                
     // 1st Electron variables
     fillVariableWithValue(   "nEle"                     , nEle_ptCut              , gen_weight * pileup_weight ); 
     //fillVariableWithValue(   "Ele1_SCEt"                , Ele1_SCEnergy/cosh(Ele1_SCEta)             , gen_weight * pileup_weight );
     fillVariableWithValue(   "Ele1_PtHeep"              , Ele1_PtHeep             , gen_weight * pileup_weight );
     fillVariableWithValue(   "Ele1_Eta"                 , Ele1_Eta                , gen_weight * pileup_weight );
     fillVariableWithValue(   "Ele1_IsBarrel"            , Ele1_IsBarrel           , gen_weight * pileup_weight );
     fillVariableWithValue(   "MTenu"                    , MT_Ele1MET              , gen_weight * pileup_weight );
     fillVariableWithValue(   "AbsDeltaEtaEleTrk"        , fabs(Ele1_Eta-Ele1_TrkEta),gen_weight * pileup_weight );

     // MET variables	                                      		           
     fillVariableWithValue(   "MET"                      , PFMET_Type1XY_Pt       , gen_weight * pileup_weight );
     fillVariableWithValue(   "Pt_EMET"                  , Pt_Ele1MET             , gen_weight * pileup_weight );
     fillVariableWithValue(   "mDeltaPhiMETEle"          , mDPhi_METEle1           , gen_weight * pileup_weight );
     									           
     // 1st JET variables                                     		           
     fillVariableWithValue(   "nJet"                     , nJet_ptCut              , gen_weight * pileup_weight );
     
     double MT_Jet1MET, MT_Jet2MET, MT_Ele1Jet1, MT_Ele1Jet2, MT_Ele1MET_Type01;
     double mDPhi_METType01_Ele1, mDPhi_METType01_Jet1, mDPhi_METType01_Jet2;
     // SIC add
     double MT_Jet1MET_Type01XY_EleEnUp, MT_Jet2MET_Type01XY_EleEnUp, MT_Ele1MET_Type01XY_EleEnUp;
     double MT_Jet1MET_Type01XY_EleEnDown, MT_Jet2MET_Type01XY_EleEnDown, MT_Ele1MET_Type01XY_EleEnDown;
     double MT_Jet1MET_Type01XY_JetEnUp, MT_Jet2MET_Type01XY_JetEnUp, MT_Ele1MET_Type01XY_JetEnUp;
     double MT_Jet1MET_Type01XY_JetEnDown, MT_Jet2MET_Type01XY_JetEnDown, MT_Ele1MET_Type01XY_JetEnDown;
     double MT_Jet1MET_Type01XY_JetResUp, MT_Jet2MET_Type01XY_JetResUp, MT_Ele1MET_Type01XY_JetResUp;
     double MT_Jet1MET_Type01XY_JetResDown, MT_Jet2MET_Type01XY_JetResDown, MT_Ele1MET_Type01XY_JetResDown;
     double mDPhi_METType01XYEleEnUp_Ele1, mDPhi_METType01XYEleEnDown_Ele1;
     double mDPhi_METType01XYEleEnUp_Jet1, mDPhi_METType01XYEleEnDown_Jet1, mDPhi_METType01XYEleEnUp_Jet2, mDPhi_METType01XYEleEnDown_Jet2;
     double mDPhi_METType01XYJetEnUp_Ele1, mDPhi_METType01XYJetEnDown_Ele1;
     double mDPhi_METType01XYJetEnUp_Jet1, mDPhi_METType01XYJetEnDown_Jet1, mDPhi_METType01XYJetEnUp_Jet2, mDPhi_METType01XYJetEnDown_Jet2;
     double mDPhi_METType01XYJetResUp_Ele1, mDPhi_METType01XYJetResDown_Ele1;
     double mDPhi_METType01XYJetResUp_Jet1, mDPhi_METType01XYJetResDown_Jet1, mDPhi_METType01XYJetResUp_Jet2, mDPhi_METType01XYJetResDown_Jet2;

     // alternate METs

     if ( nEle_store > 0 ) {
       TVector2 v_ele;
       TVector2 v_MET_Type01;
       TLorentzVector v_ele_lorentz, v_met_lorentz, v_ele_met_lorentz;
       v_ele_lorentz.SetPtEtaPhiM( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
       v_met_lorentz.SetPtEtaPhiM( PFMET_Type1XY_Pt, 0.0, PFMET_Type1XY_Phi, 0.0 );
       v_ele_met_lorentz = v_ele_lorentz + v_met_lorentz;
       

       v_ele.SetMagPhi( Ele1_Pt, Ele1_Phi );
       mDPhi_METType01_Ele1 = fabs(v_MET_Type01.DeltaPhi ( v_ele ));
       float deltaphi = v_MET_Type01.DeltaPhi(v_ele);
     }
     
     // 1st JET variables                                     		           
     if ( nJet_store > 0 ) { 						           
       fillVariableWithValue( "Jet1_Pt"                  , Jet1_Pt                 , gen_weight * pileup_weight );
       fillVariableWithValue( "Jet1_Eta"                 , Jet1_Eta                , gen_weight * pileup_weight );
       fillVariableWithValue( "mDeltaPhiMET1stJet"       , mDPhi_METJet1           , gen_weight * pileup_weight );

       TVector2 v_MET;
       TVector2 v_jet;
       TVector2 v_MET_Type01;
       v_MET.SetMagPhi( PFMET_Type1XY_Pt , PFMET_Type1XY_Phi  );
       v_jet.SetMagPhi( Jet1_Pt, Jet1_Phi );
       mDPhi_METType01_Jet1 = fabs(v_MET_Type01.DeltaPhi ( v_jet ));
       float deltaphi = v_MET.DeltaPhi(v_jet);
       MT_Jet1MET = sqrt ( 2 * Jet1_Pt * PFMET_Type1XY_Pt * ( 1 - cos ( deltaphi ) ) );
     }									           
     
     // 2nd JET variables                                     		           
     if ( nJet_store > 1 ) { 	                                      	           
       fillVariableWithValue( "Jet2_Pt"                  , Jet2_Pt                 , gen_weight * pileup_weight );
       fillVariableWithValue( "Jet2_Eta"                 , Jet2_Eta                , gen_weight * pileup_weight );
       fillVariableWithValue( "ST"                       , sT_enujj                , gen_weight * pileup_weight );
       
       TVector2 v_MET;
       TVector2 v_jet;
       TVector2 v_MET_Type01;
       v_MET.SetMagPhi( PFMET_Type1XY_Pt , PFMET_Type1XY_Phi );
       v_jet.SetMagPhi( Jet2_Pt, Jet2_Phi );
       mDPhi_METType01_Jet2 = fabs(v_MET_Type01.DeltaPhi ( v_jet ));
       float deltaphi = v_MET.DeltaPhi(v_jet);
       MT_Jet2MET = sqrt ( 2 * Jet2_Pt * PFMET_Type1XY_Pt * ( 1 - cos ( deltaphi ) ) );
     }

     // 2nd JET variables  
     double MT_ejjnu = -999;
     double MT_jjnu = -999;
     double Mejj = -999;
     if ( nJet_store > 1 ) { 	                                      	           
       TVector2 v_MET, v_dijet, v_edijet;
       TLorentzVector v_ele, v_jet1, v_jet2, v_jet1p2, v_ejj;

       v_ele .SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
       v_jet1.SetPtEtaPhiM ( Jet1_Pt, Jet1_Eta, Jet1_Phi, 0.0 );
       v_jet2.SetPtEtaPhiM ( Jet2_Pt, Jet2_Eta, Jet2_Phi, 0.0 );
       v_jet1p2 = v_jet1 + v_jet2;
       v_ejj    = v_ele + v_jet1 + v_jet2;       
       
       v_edijet.SetMagPhi ( v_ejj.Pt()     , v_ejj.Phi()     );
       v_dijet.SetMagPhi ( v_jet1p2.Pt()     , v_jet1p2.Phi()     );
       v_MET  .SetMagPhi ( PFMET_Type1XY_Pt , PFMET_Type1XY_Phi );
       
       float deltaphi = v_MET.DeltaPhi(v_dijet);
       float deltaphi2 = v_MET.DeltaPhi ( v_edijet );
       MT_jjnu  = sqrt ( 2 * v_jet1p2.Pt() * PFMET_Type1XY_Pt * ( 1 - cos ( deltaphi  ) ) );
       MT_ejjnu = sqrt ( 2 * v_ejj.Pt()    * PFMET_Type1XY_Pt * ( 1 - cos ( deltaphi2 ) ) );
       Mejj = v_ejj.M();
     }
     
     // 3rd JET variables 
     // if ( nJet_store > 2 ) {
     //   fillVariableWithValue( "Jet3_Pt"                  , Jet3_Pt                 , gen_weight * pileup_weight );
     //   fillVariableWithValue( "Jet3_Eta"                 , Jet3_Eta                , gen_weight * pileup_weight );
     // }

     // 1 electron, 1 jet variables 
     if ( nEle_ptCut > 0 && nJet_ptCut > 0 ) { 
       fillVariableWithValue ( "DR_Ele1Jet1"             , DR_Ele1Jet1             , gen_weight * pileup_weight );

       TVector2 v_ele;
       TVector2 v_jet1;
       v_ele .SetMagPhi ( Ele1_Pt, Ele1_Phi );
       v_jet1.SetMagPhi ( Jet1_Pt, Jet1_Phi );
       float deltaphi = v_ele.DeltaPhi ( v_jet1 );
       MT_Ele1Jet1 = sqrt ( 2 * Jet1_Pt * Ele1_Pt * ( 1 - cos ( deltaphi ) ) );
     }

     // 1 electron, 2 jet variables 
     if ( nEle_store > 0 && nJet_store > 1 ) { 
       fillVariableWithValue ( "DR_Ele1Jet2"             , DR_Ele1Jet2             , gen_weight * pileup_weight );

       TVector2 v_ele;
       TVector2 v_jet2;
       v_ele .SetMagPhi ( Ele1_Pt, Ele1_Phi );
       v_jet2.SetMagPhi ( Jet2_Pt, Jet2_Phi );
       float deltaphi = v_ele.DeltaPhi ( v_jet2 );
       MT_Ele1Jet2 = sqrt ( 2 * Jet2_Pt * Ele1_Pt * ( 1 - cos ( deltaphi ) ) );
     }


     // 1 electron, 3 jet variables 
     double Mej3 = 0.0;
     if ( nEle_store > 0 && nJet_store > 2 ) { 
       TLorentzVector v_ele;
       TLorentzVector v_jet3;
       TLorentzVector v_ej3;
       v_ele .SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
       v_jet3.SetPtEtaPhiM ( Jet3_Pt, Jet3_Eta, Jet3_Phi, 0.0 );
       v_ej3 = v_ele + v_jet3 ;
       Mej3 = v_ej3.M();
     }


     double M_e1e2 = 0.0;
     if ( nEle_store > 1) {
       TLorentzVector v_ele1, v_ele2, v_sum;
       v_ele1.SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
       v_ele2.SetPtEtaPhiM ( Ele2_Pt, Ele2_Eta, Ele2_Phi, 0.0 );
       v_sum = v_ele1+v_ele2;
       M_e1e2 = v_sum.M();
     }
     double M_e1e3 = 0.0;
     double M_e2e3 = 0.0;
     if ( nEle_store > 2) {
       TLorentzVector v_ele1, v_ele2, v_ele3, v_sum;
       v_ele1.SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
       v_ele2.SetPtEtaPhiM ( Ele2_Pt, Ele2_Eta, Ele2_Phi, 0.0 );
       v_ele3.SetPtEtaPhiM ( Ele3_Pt, Ele3_Eta, Ele3_Phi, 0.0 );
       v_sum = v_ele1+v_ele3;
       M_e1e3 = v_sum.M();
       v_sum = v_ele1+v_ele2;
       M_e2e3 = v_sum.M();
     }

     double MT_JetMET;
     double MT_EleJet;
     double Mej;
     bool mejSelectedJet1 = true;
     
     if ( fabs ( MT_Jet1MET - MT_Ele1Jet2 ) < fabs( MT_Jet2MET - MT_Ele1Jet1 )){
       MT_JetMET = MT_Jet1MET;
       MT_EleJet = MT_Ele1Jet2;
       Mej = M_e1j2;
       mejSelectedJet1 = false;
     } else { 
       MT_JetMET = MT_Jet2MET;
       MT_EleJet = MT_Ele1Jet1;
       Mej = M_e1j1;
     }	 
     

       double Jet1_Mass, Jet2_Mass;

       TLorentzVector jet1, jet2;
       jet1.SetPtEtaPhiE ( Jet1_Pt, Jet1_Eta, Jet1_Phi, Jet1_Energy );
       jet2.SetPtEtaPhiE ( Jet2_Pt, Jet2_Eta, Jet2_Phi, Jet2_Energy );
       Jet1_Mass = jet1.M();
       Jet2_Mass = jet2.M();

       double mDPhi_Jet1Jet2 = fabs(jet1.DeltaPhi ( jet2 ));

       double min_DR_EleJet = 999.0;
       double DR_Ele1Jet3 = 999.0;
       if ( nJet_store > 2 ) {
         TLorentzVector ele1,  jet3;
         ele1.SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
         jet3.SetPtEtaPhiM ( Jet3_Pt, Jet3_Eta, Jet3_Phi, 0.0 );	 
         DR_Ele1Jet3 = ele1.DeltaR ( jet3 ) ;
       }

       if ( DR_Ele1Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet1;
       if ( DR_Ele1Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet2;
       if ( nJet_store > 2 ) {
         if ( DR_Ele1Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet3;
       }

     // Optimization variables
     fillVariableWithValue( "ST_opt"   , sT_enujj   , gen_weight * pileup_weight );
     fillVariableWithValue( "Mej_min_opt"  , Mej        , gen_weight * pileup_weight );
     fillVariableWithValue( "MET_opt"  , PFMET_Type1XY_Pt , gen_weight * pileup_weight );
     fillVariableWithValue( "MT_opt"   , MT_Ele1MET , gen_weight * pileup_weight );
     //std::cout << "Mej_opt=" << Mej << std::endl;
     //std::cout << "ST_opt=" << sT_enujj << std::endl;
     //std::cout << "MET_opt=" << PFMET_Type1XY_Pt << std::endl;
     //std::cout << "MT_opt=" << MT_Ele1MET << std::endl;
     
     // Dummy variables
     fillVariableWithValue ("preselection",1, gen_weight * pileup_weight ); 
	
     //--------------------------------------------------------------------------
     // Fill bjet variables
     //--------------------------------------------------------------------------
     // from Dave
     //munu_wjets = '(MT_uv>70)*(MT_uv<110)*(((CISV_jet1>0.5426)+(CISV_jet2>0.5426))<1) * (2-0.887973*((1.+(0.0523821*Pt_jet1))/(1.+(0.0460876*Pt_jet1))))'
     //munu_ttbar = '(MT_uv>70)*(MT_uv<110)*(((CISV_jet1>0.8484)+(CISV_jet2>0.8484))>=1) * (0.561694*((1.+(0.31439*Pt_jet1))/(1.+(0.17756*Pt_jet1))))' 
     // 0.5426 --> CISVv2L "loose"
     // 0.8484 --> CISVv2M "medium"
     // require at least 1 medium b-tagged jet for TTBar control region
     // veto loose b-tagged jets for WJets control region
     // and then apply the b-tag scale factors as a function of Pt (2-SF in the veto case)
     // this is based on JetN_btagCISV in the ntuple
     // see: https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
     // --> CSVv2_Moriond17_B_H.csv 
     // specification of CSV file: https://twiki.cern.ch/twiki/bin/view/CMS/BTagCalibration
     // for event weights, we follow: https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1c_Event_reweighting_using_scale

     double weightZeroBJets = 1.0;
     double weightAtLeastOneBJet = 1.0;
     double weightZeroBJetsBeyondLeadingTwo = 1.0;
     double weightAtLeastOneBJetBeyondLeadingTwo = 1.0;
     double weightZeroBJetsUpShift = 1.0;
     double weightAtLeastOneBJetUpShift = 1.0;
     double weightZeroBJetsDownShift = 1.0;
     double weightAtLeastOneBJetDownShift = 1.0;

     double btagCISV_loose_cut  = 0.5426;
     double btagCISV_medium_cut = 0.8484;
     double btagCISV_tight_cut  = 0.9535;
       
     if(!isData)
     {
       // calculate and apply scale factors to MC only
       TFormula btagScaleFactorFormula("btagScaleFactorFormula","[0]*((1+([1]*x))/(1+([2]*x)))");
       // where x is jet Pt
       // comb, central, medium WP, jet flavor = 0 (B jets)
       btagScaleFactorFormula.SetParameters(0.561694,0.31439,0.17756);
       // up shift
       // we only need to worry about >= 50 GeV Pt, since that's the cut threshold
       TFormula btagScaleFactorUpFormula("btagScaleFactorUpFormula","[0]*((1+([1]*x))/(1+([2]*x)))+(0.037118)*(x>=50&&x<70)+(0.036822)*(x>=70&&x<100)+(0.034398)*(x>=100&&x<140)+(0.036239)*(x>=140&&x<200)+(0.044986)*(x>=200&&x<300)+(0.064243)*(x>=300&&x<600)+(0.097131)*(x>=600)");
       btagScaleFactorUpFormula.SetParameters(0.561694,0.31439,0.17756);
       // down shift
       // we only need to worry about >= 50 GeV Pt, since that's the cut threshold
       TFormula btagScaleFactorDownFormula("btagScaleFactorDownFormula","[0]*((1+([1]*x))/(1+([2]*x)))-(0.037118)*(x>=50&&x<70)-(0.036822)*(x>=70&&x<100)-(0.034398)*(x>=100&&x<140)-(0.036239)*(x>=140&&x<200)-(0.044986)*(x>=200&&x<300)-(0.064243)*(x>=300&&x<600)-(0.097131)*(x>=600)");
       btagScaleFactorDownFormula.SetParameters(0.561694,0.31439,0.17756);
       //
       if ( Jet1_btagCISV > btagCISV_medium_cut ) weightZeroBJets*=(1-btagScaleFactorFormula.Eval(Jet1_Pt));
       if ( Jet2_btagCISV > btagCISV_medium_cut ) weightZeroBJets*=(1-btagScaleFactorFormula.Eval(Jet2_Pt));
       if ( Jet3_btagCISV > btagCISV_medium_cut ) weightZeroBJets*=(1-btagScaleFactorFormula.Eval(Jet3_Pt));
       if ( Jet4_btagCISV > btagCISV_medium_cut ) weightZeroBJets*=(1-btagScaleFactorFormula.Eval(Jet4_Pt));
       if ( Jet5_btagCISV > btagCISV_medium_cut ) weightZeroBJets*=(1-btagScaleFactorFormula.Eval(Jet5_Pt));
       weightAtLeastOneBJet = 1 - weightZeroBJets;
       //
       if ( Jet3_btagCISV > btagCISV_medium_cut ) weightZeroBJetsBeyondLeadingTwo*=(1-btagScaleFactorFormula.Eval(Jet3_Pt));
       if ( Jet4_btagCISV > btagCISV_medium_cut ) weightZeroBJetsBeyondLeadingTwo*=(1-btagScaleFactorFormula.Eval(Jet4_Pt));
       if ( Jet5_btagCISV > btagCISV_medium_cut ) weightZeroBJetsBeyondLeadingTwo*=(1-btagScaleFactorFormula.Eval(Jet5_Pt));
       weightAtLeastOneBJetBeyondLeadingTwo = 1 - weightZeroBJetsBeyondLeadingTwo;
       //
       if ( Jet1_btagCISV > btagCISV_medium_cut ) weightZeroBJetsUpShift*=(1-btagScaleFactorUpFormula.Eval(Jet1_Pt));
       if ( Jet2_btagCISV > btagCISV_medium_cut ) weightZeroBJetsUpShift*=(1-btagScaleFactorUpFormula.Eval(Jet2_Pt));
       if ( Jet3_btagCISV > btagCISV_medium_cut ) weightZeroBJetsUpShift*=(1-btagScaleFactorUpFormula.Eval(Jet3_Pt));
       if ( Jet4_btagCISV > btagCISV_medium_cut ) weightZeroBJetsUpShift*=(1-btagScaleFactorUpFormula.Eval(Jet4_Pt));
       if ( Jet5_btagCISV > btagCISV_medium_cut ) weightZeroBJetsUpShift*=(1-btagScaleFactorUpFormula.Eval(Jet5_Pt));
       weightAtLeastOneBJetUpShift = 1 - weightZeroBJetsUpShift;
       //
       if ( Jet1_btagCISV > btagCISV_medium_cut ) weightZeroBJetsDownShift*=(1-btagScaleFactorDownFormula.Eval(Jet1_Pt));
       if ( Jet2_btagCISV > btagCISV_medium_cut ) weightZeroBJetsDownShift*=(1-btagScaleFactorDownFormula.Eval(Jet2_Pt));
       if ( Jet3_btagCISV > btagCISV_medium_cut ) weightZeroBJetsDownShift*=(1-btagScaleFactorDownFormula.Eval(Jet3_Pt));
       if ( Jet4_btagCISV > btagCISV_medium_cut ) weightZeroBJetsDownShift*=(1-btagScaleFactorDownFormula.Eval(Jet4_Pt));
       if ( Jet5_btagCISV > btagCISV_medium_cut ) weightZeroBJetsDownShift*=(1-btagScaleFactorDownFormula.Eval(Jet5_Pt));
       weightAtLeastOneBJetDownShift = 1 - weightZeroBJetsDownShift;
     }
     int nBJet_loose_ptCut  = 0;
     int nBJet_medium_ptCut = 0;
     int nBJet_tight_ptCut  = 0;
     int nBJet_medium_ptCut_beyondLeadingTwo = 0;
     
     if ( Jet1_btagCISV > btagCISV_loose_cut  ) nBJet_loose_ptCut++;
     if ( Jet2_btagCISV > btagCISV_loose_cut  ) nBJet_loose_ptCut++;
     if ( Jet3_btagCISV > btagCISV_loose_cut  ) nBJet_loose_ptCut++;
     if ( Jet4_btagCISV > btagCISV_loose_cut  ) nBJet_loose_ptCut++;
     if ( Jet5_btagCISV > btagCISV_loose_cut  ) nBJet_loose_ptCut++;
     
     // medium used below for control regions
     if ( Jet1_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut++;
     if ( Jet2_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut++;
     if ( Jet3_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut++;
     if ( Jet4_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut++;
     if ( Jet5_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut++;
     //
     if ( Jet3_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut_beyondLeadingTwo++;
     if ( Jet4_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut_beyondLeadingTwo++;
     if ( Jet5_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut_beyondLeadingTwo++;
     
     if ( Jet1_btagCISV > btagCISV_tight_cut  ) nBJet_tight_ptCut++;
     if ( Jet2_btagCISV > btagCISV_tight_cut  ) nBJet_tight_ptCut++;
     if ( Jet3_btagCISV > btagCISV_tight_cut  ) nBJet_tight_ptCut++;
     if ( Jet4_btagCISV > btagCISV_tight_cut  ) nBJet_tight_ptCut++;
     if ( Jet5_btagCISV > btagCISV_tight_cut  ) nBJet_tight_ptCut++;
     	
     //std::cout << "INFO: weightAtLeastOneBJet=" << weightAtLeastOneBJet << "; while weightZeroBJets=" << weightZeroBJets << "; this event has " << nJet_store << " stored jets, with " << 
     // nBJet_medium_ptCut << " passing the medium Btag cut." << std::endl;

     //--------------------------------------------------------------------------
     // Fill final selection cuts
     //--------------------------------------------------------------------------

     if (doFinalSelections) { 
       for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
         // SIC note: This appears to be where the final selection variables are defined (in terms of quantites calculated/known to the analysisClass)
         int lq_mass = LQ_MASS[i_lq_mass];
         // SIC edit
         //sprintf(cut_name, "MET_LQ%d"   , lq_mass ); fillVariableWithValue( cut_name , PFMET_Type1XY_Pt , gen_weight * pileup_weight );
         //sprintf(cut_name, "Mej_LQ%d"   , lq_mass ); fillVariableWithValue( cut_name , Mej        , gen_weight * pileup_weight );
         //sprintf(cut_name, "ST_LQ%d"    , lq_mass ); fillVariableWithValue( cut_name , sT_enujj   , gen_weight * pileup_weight );
         //sprintf(cut_name, "MT_LQ%d"    , lq_mass ); fillVariableWithValue( cut_name , MT_Ele1MET , gen_weight * pileup_weight );
         sprintf(cut_name, "Mej_LQ%d"   , lq_mass ); fillVariableWithValue( cut_name , Mej        , gen_weight * pileup_weight );
         sprintf(cut_name, "ST_LQ%d"    , lq_mass ); fillVariableWithValue( cut_name , sT_enujj   , gen_weight * pileup_weight );
         // SIC XXX: also change sT_enujj as defined above
         // Nominal
         sprintf(cut_name, "MT_LQ%d"    , lq_mass ); fillVariableWithValue( cut_name , MT_Ele1MET , gen_weight * pileup_weight );
         sprintf(cut_name, "MET_LQ%d"   , lq_mass ); fillVariableWithValue( cut_name , PFMET_Type1XY_Pt , gen_weight * pileup_weight );
         // UnclEnUp
         //sprintf(cut_name, "MT_LQ%d"    , lq_mass ); fillVariableWithValue( cut_name , MT_Ele1MET_Type01XY_UnclEnUp , gen_weight * pileup_weight );
         //sprintf(cut_name, "MET_LQ%d"   , lq_mass ); fillVariableWithValue( cut_name , PFMET_Type01XY_UnclUp_Pt , gen_weight * pileup_weight );
         // UnclEnDown
         //sprintf(cut_name, "MT_LQ%d"    , lq_mass ); fillVariableWithValue( cut_name , MT_Ele1MET_Type01XY_UnclEnDown , gen_weight * pileup_weight );
         //sprintf(cut_name, "MET_LQ%d"   , lq_mass ); fillVariableWithValue( cut_name , PFMET_Type01XY_UnclDown_Pt , gen_weight * pileup_weight );
         // end SIC edit
       }
     }
     
     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------

     evaluateCuts();

     if (!isData && !passedCut ("PassJSON")){
       std::cout << "ERROR: This event did not pass the JSON file!" << std::endl;
       std::cout << "  isData = " << isData << std::endl;
       std::cout << "  passedJSON = " << passedJSON << std::endl;
     }
     
     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------
     
     bool passed_preselection = passedAllPreviousCuts("preselection");
     //bool passed_minimum      = ( passedAllPreviousCuts("PassTrackingFailure") && passedCut ("PassTrackingFailure"));
     bool passed_minimum      = ( passedAllPreviousCuts("PassBadPFMuonFilter") && passedCut ("PassBadPFMuonFilter"));
    
     if ( passed_minimum && isData ){ 
       FillUserTH1D ("run_HLT", run );
     }

     //--------------------------------------------------------------------------
     // Did we pass any final selection cuts?
     //--------------------------------------------------------------------------

     if ( doFinalSelections ) { 
       passed_vector.clear();

       passed_st_vector .clear();;
       passed_mt_vector .clear();;
       passed_met_vector.clear();;
       passed_mej_vector.clear();;

       passed_stAndMej_vector .clear();
       passed_mtAndMej_vector .clear();
       passed_metAndMej_vector.clear();
       passed_stAndMet_vector .clear();
       passed_mtAndMet_vector .clear();
       passed_stAndMt_vector  .clear();
       
       passed_stAndMtAndMet_vector .clear();
       passed_stAndMtAndMej_vector .clear();
       passed_stAndMetAndMej_vector.clear();
       passed_mtAndMetAndMej_vector.clear();

       for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
         int lq_mass = LQ_MASS[i_lq_mass];
         sprintf(cut_name    ,"Mej_LQ%d", lq_mass );
         sprintf(mt_cut_name ,"MT_LQ%d" , lq_mass );
         sprintf(st_cut_name ,"ST_LQ%d" , lq_mass );
         sprintf(met_cut_name,"MET_LQ%d", lq_mass );
         sprintf(mej_cut_name,"Mej_LQ%d", lq_mass );

         bool decision     = bool ( passedAllPreviousCuts(      cut_name) && passedCut (    cut_name));
         bool st_decision  = bool ( passedAllPreviousCuts("preselection") && passedCut ( st_cut_name));
         bool mt_decision  = bool ( passedAllPreviousCuts("preselection") && passedCut ( mt_cut_name));
         bool met_decision = bool ( passedAllPreviousCuts("preselection") && passedCut (met_cut_name));
         bool mej_decision = bool ( passedAllPreviousCuts("preselection") && passedCut (mej_cut_name));

         passed_vector.push_back (decision);

         passed_st_vector .push_back ( st_decision  );
         passed_mt_vector .push_back ( mt_decision  );
         passed_met_vector.push_back ( met_decision );
         passed_mej_vector.push_back ( mej_decision );

         passed_stAndMej_vector .push_back( st_decision  && mej_decision );
         passed_mtAndMej_vector .push_back( mt_decision  && mej_decision );
         passed_metAndMej_vector.push_back( met_decision && mej_decision );
         passed_stAndMet_vector .push_back( st_decision  && met_decision );
         passed_mtAndMet_vector .push_back( mt_decision  && met_decision );
         passed_stAndMt_vector  .push_back( st_decision  && mt_decision  );

         passed_stAndMtAndMet_vector .push_back(st_decision && mt_decision  && met_decision ); 
         passed_stAndMtAndMej_vector .push_back(st_decision && mt_decision  && mej_decision ); 
         passed_stAndMetAndMej_vector.push_back(st_decision && met_decision && mej_decision ); 
         passed_mtAndMetAndMej_vector.push_back(mt_decision && met_decision && mej_decision ); 

       }
     }
     
     FillUserTH1D( "ProcessID"      , ProcessID, pileup_weight * gen_weight ) ;
     
     // for 3 electron event dists
     // passed preselection without nEle cut
     if(passedAllOtherCuts("nEle")) {
       if(nEle_store > 1)
         FillUserTH1D("Mee_allElectrons_Presel", M_e1e2, pileup_weight * gen_weight);
       if(nEle_store > 2) {
         FillUserTH1D("Mee_allElectrons_Presel", M_e1e3, pileup_weight * gen_weight);
         FillUserTH1D("Mee_allElectrons_Presel", M_e2e3, pileup_weight * gen_weight);
         FillUserTH1D("Mee_allElectrons_3EleEvents_Presel", M_e1e3, pileup_weight * gen_weight);
         FillUserTH1D("Mee_allElectrons_3EleEvents_Presel", M_e2e3, pileup_weight * gen_weight);
         FillUserTH1D("Mee_allElectrons_3EleEvents_Presel", M_e1e2, pileup_weight * gen_weight);
       }
     }

     if ( passed_preselection ) { 

       //std::cout.precision(0);
       //if(MT_Ele1MET > 1000 || PFMET_Type1XY_Pt > 1000) {
       //  // run ls event
       //  if(current_file_name.find("SingleElectron") != std::string::npos) {
       //    std::cout << fixed << "[Preselection SingleElectron] passing run/ls/event: " << static_cast<int>(run) << " " << static_cast<int>(ls) << " " << ((unsigned int)event) << std::endl;
       //  }
       //  else if(current_file_name.find("SinglePhoton") != std::string::npos) {
       //    std::cout << fixed << "[Preselection SinglePhoton] passing run/ls/event: " << static_cast<int>(run) << " " << static_cast<int>(ls) << " " << ((unsigned int)event) << std::endl;
       //  }
       //  std::cout << "in filename: " << current_file_name << std::endl;
       //  // kinematics
       //  std::cout.precision(3);
       //  std::cout << fixed << "Mej      = " << Mej << std::endl;
       //  std::cout << fixed << "sT       = " << sT_enujj << std::endl;
       //  std::cout << fixed << "MTenu = " << MT_Ele1MET << std::endl;
       //  std::cout << fixed << "MET Pt  = " << PFMET_Type1XY_Pt << "\t, Phi = " << PFMET_Type1XY_Phi << std::endl;
       //  //std::cout << fixed << "Mee      = " << M_e1e2 << std::endl;
       //  std::cout << fixed << "Ele1 Pt  = " << Ele1_Pt << "\t, Eta = " << Ele1_Eta << "\t, Phi = " << Ele1_Phi << std::endl;
       //  std::cout << fixed << "Ele1 dPhi  = " << Ele1_DeltaPhiTrkSC << "\t, HoE = " << Ele1_HoE << "\t, sIetaIeta = " << Ele1_Full5x5SigmaIEtaIEta << std::endl;
       //  //std::cout << fixed << "  Ele2 Pt  = " << Ele2_Pt << "\t, Eta = " << Ele2_Eta << "\t, Phi = " << Ele2_Phi << std::endl;
       //  //std::cout << fixed << "  Ele2 dPhi  = " << Ele2_DeltaPhiTrkSC << "\t, HoE = " << Ele2_HoE << "\t, sIetaIeta = " << Ele2_SigmaIEtaIEta << std::endl;
       //  std::cout << fixed << "Jet1 Pt  = " << Jet1_Pt << "\t, Eta = " << Jet1_Eta << "\t, Phi = " << Jet1_Phi << std::endl;
       //  std::cout << fixed << "Jet2 Pt  = " << Jet2_Pt << "\t, Eta = " << Jet2_Eta << "\t, Phi = " << Jet2_Phi << std::endl;
       //  std::cout << fixed << "deltaPhi(jet1,MET) = " << mDPhi_METJet1 << "\t, deltaPhi(jet2,MET) = " << mDPhi_METJet2 << std::endl;
       //  std::cout << fixed << "deltaPhi(ele,MET) = " << mDPhi_METEle1 << std::endl;
       //  std::cout << fixed << "MT(j,MET) = " << MT_JetMET << std::endl;
       //  std::cout << fixed << "min[dR(ele,jets)] = " << min_DR_EleJet << std::endl;
       //  std::cout << std::endl;
       //}


       FillUserTH1D( "ProcessID_PAS"      , ProcessID, pileup_weight * gen_weight ) ;

       //--------------------------------------------------------------------------
       // Fill control region
       //--------------------------------------------------------------------------
       
       
       if ( MT_Ele1MET  > 110. && 
           sT_enujj > 300. && 
           sT_enujj < 495. ){
         if ( nJet_ptCut >= 4 ){ 
           FillUserTH1D("sT_ControlRegion_Njet_gte4", sT_enujj, pileup_weight * gen_weight );
         }
         if ( nJet_ptCut <= 3 ){
           FillUserTH1D("sT_ControlRegion_Njet_lte3", sT_enujj, pileup_weight * gen_weight );
         }
       }

       //--------------------------------------------------------------------------
       // Fill skim tree, if necessary
       //--------------------------------------------------------------------------

       if ( isData ) {
         FillUserTH1D("run_PAS", run ) ;
         //FillUserTH1D("split_PAS"    , get_split     ( run ), pileup_weight * gen_weight ) ;
         //FillUserTH1D("split_1fb_PAS", get_split_1fb ( run ), pileup_weight * gen_weight ) ;
       }

       FillUserTH1D( "nElectron_PAS"              , nEle_ptCut                                   , pileup_weight * gen_weight); 
       FillUserTH1D( "nMuon_PAS"                  , nMuon_ptCut                                  , pileup_weight * gen_weight); 
       FillUserTH1D( "EOverP_1stEle_PAS"          , Ele1_EOverP                                  , pileup_weight * gen_weight); 
       FillUserTH1D( "Pt1stEle_PAS"	          , Ele1_Pt                                      , pileup_weight * gen_weight); 
       FillUserTH1D( "Energy1stEle_PAS"	          , Ele1_Energy                                  , pileup_weight * gen_weight); 
       FillUserTH1D( "Eta1stEle_PAS"	          , Ele1_Eta                                     , pileup_weight * gen_weight);
       FillUserTH1D( "Phi1stEle_PAS"	          , Ele1_Phi                                     , pileup_weight * gen_weight);
       FillUserTH1D( "Charge1stEle_PAS"           , Ele1_Charge                                  , pileup_weight * gen_weight);   
       FillUserTH1D( "MET_PAS"                    , PFMET_Type1XY_Pt                            , pileup_weight * gen_weight);
       // muon kinematics
       FillUserTH1D( "Pt1stMuon_PAS"	          , Muon1_Pt                                      , pileup_weight * gen_weight); 
       FillUserTH1D( "Eta1stMuon_PAS"	          , Muon1_Eta                                     , pileup_weight * gen_weight);
       FillUserTH1D( "Phi1stMuon_PAS"	          , Muon1_Phi                                     , pileup_weight * gen_weight);
       FillUserTH1D( "Pt2ndMuon_PAS"	          , Muon2_Pt                                      , pileup_weight * gen_weight); 
       FillUserTH1D( "Eta2ndMuon_PAS"	          , Muon2_Eta                                     , pileup_weight * gen_weight);
       FillUserTH1D( "Phi2ndMuon_PAS"	          , Muon2_Phi                                     , pileup_weight * gen_weight);
       //// SIC ADD
       //FillUserTH1D( "MET_PAS_UnclUp"                    , PFMET_Type01XY_UnclUp_Pt                            , pileup_weight * gen_weight);
       //FillUserTH1D( "MET_PAS_UnclDown"                    , PFMET_Type01XY_UnclDown_Pt                            , pileup_weight * gen_weight);
       //// END SIC ADD
       FillUserTH1D( "MET_PDF"                    , PFMET_Type1XY_Pt                            , pileup_weight * gen_weight);
       FillUserTH1D( "METPhi_PAS"	          , PFMET_Type1XY_Phi                           , pileup_weight * gen_weight);   
       //FillUserTH1D( "MET_Type01_PAS"             , PFMET_Type01_Pt                              , pileup_weight * gen_weight);
       //FillUserTH1D( "MET_Type01_Phi_PAS"	  , PFMET_Type01_Phi                             , pileup_weight * gen_weight);   
       //FillUserTH1D( "MET_Type1_PAS"             , PFMET_Type1_Pt                              , pileup_weight * gen_weight);
       //FillUserTH1D( "MET_Type1_Phi_PAS"	  , PFMET_Type1_Phi                             , pileup_weight * gen_weight);   
       FillUserTH1D( "minMETPt1stEle_PAS"         , TMath::Min ( Ele1_Pt, PFMET_Type1XY_Pt  )   , pileup_weight * gen_weight);
       //FillUserTH1D( "minMET01Pt1stEle_PAS"       , TMath::Min ( Ele1_Pt, PFMET_Type01_Pt    )   , pileup_weight * gen_weight);
       //FillUserTH1D( "minMET1Pt1stEle_PAS"       , TMath::Min ( Ele1_Pt, PFMET_Type1_Pt    )   , pileup_weight * gen_weight);
       FillUserTH1D( "Pt1stJet_PAS"               , Jet1_Pt                                      , pileup_weight * gen_weight);
       FillUserTH1D( "Pt2ndJet_PAS"               , Jet2_Pt                                      , pileup_weight * gen_weight);
       FillUserTH1D( "Eta1stJet_PAS"              , Jet1_Eta                                     , pileup_weight * gen_weight);
       FillUserTH1D( "Eta2ndJet_PAS"              , Jet2_Eta                                     , pileup_weight * gen_weight);
       FillUserTH1D( "Phi1stJet_PAS"              , Jet1_Phi                                     , pileup_weight * gen_weight);
       FillUserTH1D( "Phi2ndJet_PAS"	          , Jet2_Phi                                     , pileup_weight * gen_weight);
       FillUserTH1D( "Mass1stJet_PAS"             , Jet1_Mass                                    , pileup_weight * gen_weight);
       FillUserTH1D( "Mass2ndJet_PAS"             , Jet2_Mass                                    , pileup_weight * gen_weight);
       FillUserTH1D( "CISV1stJet_PAS"              , Jet1_btagCISV                                 , pileup_weight * gen_weight);
       FillUserTH1D( "CISV2ndJet_PAS"              , Jet2_btagCISV                                 , pileup_weight * gen_weight);
       FillUserTH1D( "nMuon_PtCut_IDISO_PAS"      , nMuon_ptCut                                  , pileup_weight * gen_weight); 
       FillUserTH1D( "MTenu_PAS"                  , MT_Ele1MET                                   , pileup_weight * gen_weight);
       FillUserTH1D( "MT_PDF"                     , MT_Ele1MET                                   , pileup_weight * gen_weight);
       //FillUserTH1D( "MTenu_Type01_PAS"           , MT_Ele1MET_Type01                            , pileup_weight * gen_weight);
       FillUserTH1D( "Ptenu_PAS"	          , Pt_Ele1MET                                   , pileup_weight * gen_weight);
       FillUserTH1D( "sTlep_PAS"                  , Ele1_Pt + PFMET_Type1XY_Pt                  , pileup_weight * gen_weight);
       //FillUserTH1D( "sTlep_Type01_PAS"           , Ele1_Pt + PFMET_Type01_Pt                    , pileup_weight * gen_weight);
       //FillUserTH1D( "sTlep_Type1_PAS"           , Ele1_Pt + PFMET_Type1_Pt                    , pileup_weight * gen_weight);
       FillUserTH1D( "sT_PAS"                     , sT_enujj                                     , pileup_weight * gen_weight);
       FillUserTH1D( "sT_PDF"                     , sT_enujj                                     , pileup_weight * gen_weight);
       //FillUserTH1D( "sT_Type01_PAS"              , Ele1_Pt + PFMET_Type01_Pt + Jet1_Pt + Jet2_Pt, pileup_weight * gen_weight);
       //FillUserTH1D( "sT_Type1_PAS"              , Ele1_Pt + PFMET_Type1_Pt + Jet1_Pt + Jet2_Pt, pileup_weight * gen_weight);
       FillUserTH1D( "sTjet_PAS"                  , Jet1_Pt + Jet2_Pt                            , pileup_weight * gen_weight);
       FillUserTH1D( "Mjj_PAS"	                  , M_j1j2                                       , pileup_weight * gen_weight);   
       FillUserTH1D( "DCotTheta1stEle_PAS"        , Ele1_DCotTheta                               , pileup_weight * gen_weight);
       FillUserTH1D( "Dist1stEle_PAS"             , Ele1_Dist                                    , pileup_weight * gen_weight);
       FillUserTH1D( "mDPhiJet1Jet2_PAS"          , mDPhi_Jet1Jet2                               , pileup_weight * gen_weight);
       FillUserTH1D( "mDPhi1stEleMET_PAS"         , mDPhi_METEle1                                , pileup_weight * gen_weight);
       FillUserTH1D( "mDPhi1stJetMET_PAS"         , mDPhi_METJet1                                , pileup_weight * gen_weight);
       FillUserTH1D( "mDPhi2ndJetMET_PAS"         , mDPhi_METJet2                                , pileup_weight * gen_weight); 
       //FillUserTH1D( "mDPhi1stEleMET_Type01_PAS"  , mDPhi_METType01_Ele1                         , pileup_weight * gen_weight);
       //FillUserTH1D( "mDPhi1stJetMET_Type01_PAS"  , mDPhi_METType01_Jet1                         , pileup_weight * gen_weight);
       //FillUserTH1D( "mDPhi2ndJetMET_Type01_PAS"  , mDPhi_METType01_Jet2                         , pileup_weight * gen_weight); 
       FillUserTH1D( "Mej1_PAS"                   , M_e1j1                                       , pileup_weight * gen_weight);
       FillUserTH1D( "Mej2_PAS"                   , M_e1j2                                       , pileup_weight * gen_weight);
       FillUserTH1D( "Mej_PAS"                    , Mej                                          , pileup_weight * gen_weight);
       FillUserTH1D( "Mej_PDF"                    , Mej                                          , pileup_weight * gen_weight);
       FillUserTH1D( "Mejj_PAS"                   , Mejj                                         , pileup_weight * gen_weight);
       FillUserTH1D( "MTjnu_PAS"                  , MT_JetMET                                    , pileup_weight * gen_weight);
       FillUserTH1D( "MTjjnu_PAS"                 , MT_jjnu                                      , pileup_weight * gen_weight);
       FillUserTH1D( "DR_Ele1Jet1_PAS"	          , DR_Ele1Jet1                                  , pileup_weight * gen_weight);
       FillUserTH1D( "DR_Ele1Jet2_PAS"	          , DR_Ele1Jet2                                  , pileup_weight * gen_weight);
       FillUserTH1D( "DR_Jet1Jet2_PAS"	          , DR_Jet1Jet2                                  , pileup_weight * gen_weight);
       FillUserTH1D( "minDR_EleJet_PAS"           , min_DR_EleJet                                , pileup_weight * gen_weight);
       FillUserTH1D( "nVertex_PAS"                , nVertex                                      , pileup_weight * gen_weight);
       FillUserTH1D( "nJet_PAS"                   , nJet_ptCut                                   , pileup_weight * gen_weight);
       FillUserTH1D( "GeneratorWeight"            , gen_weight                                                               );
       FillUserTH1D( "PileupWeight"               , pileup_weight                                                            );
       //
       FillUserTH1D( "DeltaPhi1stEle_Presel"             , Ele1_DeltaPhiTrkSC                           , pileup_weight * gen_weight);
       FillUserTH1D( "SCEta1stEle_Presel"                , Ele1_SCEta                                   , pileup_weight * gen_weight);
       FillUserTH1D( "PtHeep1stEle_Presel"               , Ele1_PtHeep                                  , pileup_weight * gen_weight);
       FillUserTH1D( "SCEt1stEle_Presel"                 , Ele1_SCEnergy/cosh(Ele1_SCEta)               , pileup_weight * gen_weight);
       FillUserTH1D( "HoE1stEle_Presel"                  , Ele1_HoE                                     , pileup_weight * gen_weight);
       FillUserTH1D( "EcalIso1stEle_Presel"              , Ele1_EcalIsolation                           , pileup_weight * gen_weight);
       FillUserTH1D( "HcalIso1stEle_Presel"              , Ele1_HcalIsolation                           , pileup_weight * gen_weight);
       FillUserTH1D( "TrkIsoHEEP701stEle_Presel"         , Ele1_TrkIsoHEEP7                             , pileup_weight * gen_weight);
       FillUserTH1D( "LeadVtxDistXY1stEle_Presel"        , Ele1_LeadVtxDistXY                           , pileup_weight * gen_weight);
       FillUserTH1D( "DeltaEtaEleTrk1stEle_Presel"       , fabs(Ele1_Eta-Ele1_TrkEta)                   , pileup_weight * gen_weight);

       if(Ele1_Pt >= 50 && Ele1_Pt <= 100)
       {
         FillUserTH1D( "DeltaPhi1stEle_Pt50to100"             , Ele1_DeltaPhiTrkSC                           , pileup_weight * gen_weight);
         FillUserTH1D( "SCEta1stEle_Pt50to100"                , Ele1_SCEta                                   , pileup_weight * gen_weight);
         FillUserTH1D( "PtHeep1stEle_Pt50to100"               , Ele1_PtHeep                                  , pileup_weight * gen_weight);
         FillUserTH1D( "HoE1stEle_Pt50to100"                  , Ele1_HoE                                     , pileup_weight * gen_weight);
         FillUserTH1D( "EcalIso1stEle_Pt50to100"              , Ele1_EcalIsolation                           , pileup_weight * gen_weight);
         FillUserTH1D( "HcalIso1stEle_Pt50to100"              , Ele1_HcalIsolation                           , pileup_weight * gen_weight);
         FillUserTH1D( "TrkIsoHEEP701stEle_Pt50to100"         , Ele1_TrkIsoHEEP7                             , pileup_weight * gen_weight);
         FillUserTH1D( "LeadVtxDistXY1stEle_Pt50to100"        , Ele1_LeadVtxDistXY                           , pileup_weight * gen_weight);
       }
       else if(Ele1_Pt >= 130 && Ele1_Pt <= 200)
       {
         FillUserTH1D( "DeltaPhi1stEle_Pt130to200"             , Ele1_DeltaPhiTrkSC                           , pileup_weight * gen_weight);
         FillUserTH1D( "SCEta1stEle_Pt130to200"                , Ele1_SCEta                                   , pileup_weight * gen_weight);
         FillUserTH1D( "PtHeep1stEle_Pt130to200"               , Ele1_PtHeep                                  , pileup_weight * gen_weight);
         FillUserTH1D( "HoE1stEle_Pt130to200"                  , Ele1_HoE                                     , pileup_weight * gen_weight);
         FillUserTH1D( "EcalIso1stEle_Pt130to200"              , Ele1_EcalIsolation                           , pileup_weight * gen_weight);
         FillUserTH1D( "HcalIso1stEle_Pt130to200"              , Ele1_HcalIsolation                           , pileup_weight * gen_weight);
         FillUserTH1D( "TrkIsoHEEP701stEle_Pt130to200"         , Ele1_TrkIsoHEEP7                             , pileup_weight * gen_weight);
         FillUserTH1D( "LeadVtxDistXY1stEle_Pt130to200"        , Ele1_LeadVtxDistXY                           , pileup_weight * gen_weight);
       }

       if ( fabs(Ele1_Eta) <= eleEta_bar ) { 
         FillUserTH1D( "Mej_Barrel_Presel"                     , Mej                                          , pileup_weight * gen_weight);
         FillUserTH1D( "Pt1stEle_Barrel_PAS"                  , Ele1_Pt                                       , pileup_weight * gen_weight);
       }
       else if ( fabs(Ele1_Eta) >= eleEta_end1_min && fabs(Ele1_Eta) < eleEta_end1_max) {
         FillUserTH1D( "Mej_Endcap1_Presel"                    , Mej                                          , pileup_weight * gen_weight);
         FillUserTH1D( "Pt1stEle_Endcap1_PAS"                  , Ele1_Pt                                      , pileup_weight * gen_weight);
       }
       else if ( fabs(Ele1_Eta) >= eleEta_end2_min && fabs(Ele1_Eta) < eleEta_end2_max) {
         FillUserTH1D( "Mej_Endcap2_Presel"                    , Mej                                          , pileup_weight * gen_weight);
         FillUserTH1D( "Pt1stEle_Endcap2_PAS"                  , Ele1_Pt                                      , pileup_weight * gen_weight);
       }

       //likelihood_values.clear();
       //likelihood_values.push_back ( Mej );  
       //likelihood_values.push_back ( sT_enujj  );
       //likelihood_values.push_back ( PFMET_Type1XY_Pt   );
       //likelihood_values.push_back ( MT_Ele1MET );

       // possibly FIXME
       //for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
       //  int  lq_mass = LQ_MASS      [i_lq_mass];
       //  char signal_name[100]; sprintf( signal_name, "LQ_M%d", lq_mass );
       //  double likelihood     = myLikelihoodGetter -> getLikelihood    ( signal_name, likelihood_values ) ;
       //  double log_likelihood = myLikelihoodGetter -> getLogLikelihood ( signal_name, likelihood_values ) ;
       //  sprintf(plot_name, "likelihood_LQ%d"   , lq_mass ); FillUserTH1D ( plot_name, likelihood    , pileup_weight * gen_weight );	 
       //  sprintf(plot_name, "logLikelihood_LQ%d", lq_mass ); FillUserTH1D ( plot_name, log_likelihood, pileup_weight * gen_weight );	 
       //}

       if ( Pt_Ele1MET > 200. && 
           Jet1_Pt    > 200. && 
           Jet1_Mass  > 50.  ) {
         FillUserTH1D("nJet_PASandFrancesco", nJet_ptCut, pileup_weight * gen_weight);
       }

       if ( nJet_ptCut == 2 ) FillUserTH1D( "Eta1stJet_PASand2Jet", Jet1_Eta, pileup_weight * gen_weight);
       if ( nJet_ptCut == 3 ) FillUserTH1D( "Eta1stJet_PASand3Jet", Jet1_Eta, pileup_weight * gen_weight);
       if ( nJet_ptCut == 4 ) FillUserTH1D( "Eta1stJet_PASand4Jet", Jet1_Eta, pileup_weight * gen_weight);
       if ( nJet_ptCut == 5 ) FillUserTH1D( "Eta1stJet_PASand5Jet", Jet1_Eta, pileup_weight * gen_weight);

       if ( Ele1_Pt > 40 && Ele1_Pt < 45 && fabs ( Ele1_Eta ) > 2.1 ){
         FillUserTH1D( "Pt1stEle_Pt40to45_EtaGT2p1", Ele1_Pt, pileup_weight * gen_weight);
       }

       if ( MT_Ele1MET > 50 && MT_Ele1MET < 110 ){

         FillUserTH1D( "ProcessID_WWindow", ProcessID, pileup_weight * gen_weight );
         FillUserTH1D( "MTenu_50_110"      , MT_Ele1MET,  pileup_weight * gen_weight ) ;
         FillUserTH1D( "nJets_MTenu_50_110", nJet_ptCut,  pileup_weight * gen_weight ) ;

         if ( nJet_ptCut <= 3 ){ 
           FillUserTH1D(   "MTenu_50_110_Njet_lte3"         , MT_Ele1MET       ,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_lte3", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_lte3", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_lte3", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_50_110_Njet_lte3"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_50_110_Njet_lte3"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_50_110_Njet_lte3"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut <= 4 ){ 
           FillUserTH1D(   "MTenu_50_110_Njet_lte4", MT_Ele1MET,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_lte4", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_lte4", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_lte4", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_50_110_Njet_lte4"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_50_110_Njet_lte4"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_50_110_Njet_lte4"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut >= 4 ){ 
           FillUserTH1D(   "MTenu_50_110_Njet_gte4", MT_Ele1MET,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_gte4", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_gte4", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_gte4", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_50_110_Njet_gte4"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_50_110_Njet_gte4"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_50_110_Njet_gte4"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut >= 5 ){ 
           FillUserTH1D(   "MTenu_50_110_Njet_gte5", MT_Ele1MET,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_gte5", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_gte5", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_gte5", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_50_110_Njet_gte5"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_50_110_Njet_gte5"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_50_110_Njet_gte5"     , Mej              ,  pileup_weight * gen_weight ) ;
         }
       }

       if ( MT_Ele1MET > 70 && MT_Ele1MET < 150 ){

         FillUserTH1D( "ProcessID_WWindow", ProcessID, pileup_weight * gen_weight );
         FillUserTH1D( "MTenu_70_150"      , MT_Ele1MET,  pileup_weight * gen_weight ) ;
         FillUserTH1D( "nJets_MTenu_70_150", nJet_ptCut,  pileup_weight * gen_weight ) ;

         if ( nJet_ptCut <= 3 ){ 
           FillUserTH1D(   "MTenu_70_150_Njet_lte3"         , MT_Ele1MET       ,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_lte3", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_lte3", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_lte3", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_70_150_Njet_lte3"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_70_150_Njet_lte3"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_70_150_Njet_lte3"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut <= 4 ){ 
           FillUserTH1D(   "MTenu_70_150_Njet_lte4", MT_Ele1MET,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_lte4", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_lte4", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_lte4", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_70_150_Njet_lte4"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_70_150_Njet_lte4"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_70_150_Njet_lte4"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut >= 4 ){ 
           FillUserTH1D(   "MTenu_70_150_Njet_gte4", MT_Ele1MET,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_gte4", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_gte4", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_gte4", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_70_150_Njet_gte4"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_70_150_Njet_gte4"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_70_150_Njet_gte4"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut >= 5 ){ 
           FillUserTH1D(   "MTenu_70_150_Njet_gte5", MT_Ele1MET,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_gte5", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_gte5", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_gte5", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_70_150_Njet_gte5"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_70_150_Njet_gte5"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_70_150_Njet_gte5"     , Mej              ,  pileup_weight * gen_weight ) ;
         }
       }

       if ( MT_Ele1MET > 110 && MT_Ele1MET < 190 ){

         FillUserTH1D( "MTenu_110_190"      , MT_Ele1MET,  pileup_weight * gen_weight ) ;
         FillUserTH1D( "nJets_MTenu_110_190", nJet_ptCut,  pileup_weight * gen_weight ) ;

         if ( nJet_ptCut <= 3 ){ 
           FillUserTH1D(   "MTenu_110_190_Njet_lte3"         , MT_Ele1MET       ,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_110_190_Njet_lte3", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_110_190_Njet_lte3", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_110_190_Njet_lte3", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_110_190_Njet_lte3"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_110_190_Njet_lte3"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_110_190_Njet_lte3"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut <= 4 ){ 
           FillUserTH1D(   "MTenu_110_190_Njet_lte4", MT_Ele1MET,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_110_190_Njet_lte4", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_110_190_Njet_lte4", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_110_190_Njet_lte4", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_110_190_Njet_lte4"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_110_190_Njet_lte4"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_110_190_Njet_lte4"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut >= 4 ){ 
           FillUserTH1D(   "MTenu_110_190_Njet_gte4", MT_Ele1MET,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_110_190_Njet_gte4", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_110_190_Njet_gte4", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_110_190_Njet_gte4", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_110_190_Njet_gte4"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_110_190_Njet_gte4"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_110_190_Njet_gte4"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut >= 5 ){ 
           FillUserTH1D(   "MTenu_110_190_Njet_gte5", MT_Ele1MET,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_110_190_Njet_gte5", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_110_190_Njet_gte5", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_110_190_Njet_gte5", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_110_190_Njet_gte5"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_110_190_Njet_gte5"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_110_190_Njet_gte5"     , Mej              ,  pileup_weight * gen_weight ) ;
         }
       }


       //-------------------------------------------------------------------------- 
       // no b tags
       //-------------------------------------------------------------------------- 
       if((isData && nBJet_medium_ptCut==0) || !isData) {
         if ( MT_Ele1MET > 70 && MT_Ele1MET < 110 ){

         //  FillUserTH1D( "ProcessID_WWindow", ProcessID, pileup_weight * gen_weight * weightZeroBJets );
         FillUserTH1D( "MTenu_70_110_noBtaggedJets"      , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //  FillUserTH1D( "nJets_MTenu_70_110_noBtaggedJets", nJet_ptCut,  pileup_weight * gen_weight * weightZeroBJets ) ;

         //  if ( nJet_ptCut <= 3 ){ 
         //    FillUserTH1D(   "MTenu_70_110_Njet_lte3_noBtaggedJets"         , MT_Ele1MET       ,  pileup_weight * gen_weight * weightZeroBJets ) ;

         //    FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_lte3_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_lte3_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_lte3_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "MET_MTenu_70_110_Njet_lte3_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "sT_MTenu_70_110_Njet_lte3_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "Mej_MTenu_70_110_Njet_lte3_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //  }

         //  if ( nJet_ptCut <= 4 ){ 
         //    FillUserTH1D(   "MTenu_70_110_Njet_lte4_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;

         //    FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_lte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_lte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_lte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "MET_MTenu_70_110_Njet_lte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "sT_MTenu_70_110_Njet_lte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "Mej_MTenu_70_110_Njet_lte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //  }

         //  if ( nJet_ptCut >= 4 ){ 
         //    FillUserTH1D(   "MTenu_70_110_Njet_gte4_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;

         //    FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_gte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_gte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_gte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "MET_MTenu_70_110_Njet_gte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "sT_MTenu_70_110_Njet_gte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "Mej_MTenu_70_110_Njet_gte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //  }

         //  if ( nJet_ptCut >= 5 ){ 
         //    FillUserTH1D(   "MTenu_70_110_Njet_gte5_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;

         //    FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_gte5_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_gte5_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_gte5_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "MET_MTenu_70_110_Njet_gte5_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "sT_MTenu_70_110_Njet_gte5_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //    FillUserTH1D(   "Mej_MTenu_70_110_Njet_gte5_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
         //  }
         }

         if ( MT_Ele1MET > 50 && MT_Ele1MET < 110 ){

           FillUserTH1D( "ProcessID_WWindow", ProcessID, pileup_weight * gen_weight * weightZeroBJets );
           FillUserTH1D( "MTenu_50_110_noBtaggedJets"      , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;
           FillUserTH1D( "MTenu_50_110_noBtaggedJets_btagSFUpShift"      , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJetsUpShift ) ;
           FillUserTH1D( "MTenu_50_110_noBtaggedJets_btagSFDownShift"      , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJetsDownShift ) ;
           FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets", nJet_ptCut,  pileup_weight * gen_weight * weightZeroBJets ) ;
           // scale factor dependence plots
           if(sT_enujj > 300 && sT_enujj < 500)
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_sT300To500_PAS"		         , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ); 
           else if(sT_enujj > 500 && sT_enujj < 750)
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_sT500To750_PAS"		         , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ); 
           else if(sT_enujj > 750 && sT_enujj < 1250)
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_sT750To1250_PAS"		         , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ); 
           else if(sT_enujj > 1250)
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_sT1250ToInf_PAS"		         , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ); 
           if(Mej > 100 && Mej < 200)
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_Mej100To200_PAS"		     , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ); 
           else if(Mej > 200 && Mej < 300)
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_Mej200To300_PAS"		     , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ); 
           else if(Mej > 300 && Mej < 400)
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_Mej300To400_PAS"		     , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ); 
           else if(Mej > 400 && Mej < 500)
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_Mej400To500_PAS"		     , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ); 
           else if(Mej > 500 && Mej < 650)
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_Mej500To650_PAS"		     , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ); 
           else if(Mej > 650)
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_Mej650ToInf_PAS"		     , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ); 

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_lte3_noBtaggedJets"         , MT_Ele1MET       ,  pileup_weight * gen_weight * weightZeroBJets ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_lte3_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_lte3_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_lte3_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_lte3_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_lte3_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_lte3_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_lte4_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_lte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_lte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_lte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_lte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_lte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_lte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_gte4_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_gte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_gte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_gte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_gte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_gte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_gte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_gte5_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_gte5_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_gte5_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_gte5_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_gte5_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_gte5_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_gte5_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
           }
         }

         if ( MT_Ele1MET > 70 && MT_Ele1MET < 150 ){

           FillUserTH1D( "ProcessID_WWindow", ProcessID, pileup_weight * gen_weight * weightZeroBJets );
           FillUserTH1D( "MTenu_70_150_noBtaggedJets"      , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;
           FillUserTH1D( "nJets_MTenu_70_150_noBtaggedJets", nJet_ptCut,  pileup_weight * gen_weight * weightZeroBJets ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_lte3_noBtaggedJets"         , MT_Ele1MET       ,  pileup_weight * gen_weight * weightZeroBJets ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_lte3_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_lte3_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_lte3_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_lte3_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_lte3_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_lte3_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_lte4_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_lte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_lte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_lte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_lte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_lte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_lte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_gte4_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_gte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_gte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_gte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_gte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_gte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_gte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_gte5_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_gte5_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_gte5_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_gte5_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_gte5_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_gte5_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_gte5_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
           }
         }

         if ( MT_Ele1MET > 110 && MT_Ele1MET < 190 ){

           FillUserTH1D( "MTenu_110_190_noBtaggedJets"      , MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;
           FillUserTH1D( "nJets_MTenu_110_190_noBtaggedJets", nJet_ptCut,  pileup_weight * gen_weight * weightZeroBJets ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_110_190_Njet_lte3_noBtaggedJets"         , MT_Ele1MET       ,  pileup_weight * gen_weight * weightZeroBJets ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_110_190_Njet_lte3_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_110_190_Njet_lte3_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_110_190_Njet_lte3_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "MET_MTenu_110_190_Njet_lte3_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "sT_MTenu_110_190_Njet_lte3_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Mej_MTenu_110_190_Njet_lte3_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_110_190_Njet_lte4_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_110_190_Njet_lte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_110_190_Njet_lte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_110_190_Njet_lte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "MET_MTenu_110_190_Njet_lte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "sT_MTenu_110_190_Njet_lte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Mej_MTenu_110_190_Njet_lte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_110_190_Njet_gte4_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_110_190_Njet_gte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_110_190_Njet_gte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_110_190_Njet_gte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "MET_MTenu_110_190_Njet_gte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "sT_MTenu_110_190_Njet_gte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Mej_MTenu_110_190_Njet_gte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_110_190_Njet_gte5_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight * weightZeroBJets ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_110_190_Njet_gte5_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_110_190_Njet_gte5_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_110_190_Njet_gte5_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "MET_MTenu_110_190_Njet_gte5_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "sT_MTenu_110_190_Njet_gte5_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight * weightZeroBJets ) ;
             FillUserTH1D(   "Mej_MTenu_110_190_Njet_gte5_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight * weightZeroBJets ) ;
           }
         }

         if(MT_Ele1MET >= 200 && MT_Ele1MET < 400)
           FillUserTH1D( "MTenu_noBtaggedJets_MT200To400_PAS" , MT_Ele1MET, pileup_weight * gen_weight * weightZeroBJets ) ;
         else if(MT_Ele1MET >= 400 && MT_Ele1MET < 600)
           FillUserTH1D( "MTenu_noBtaggedJets_MT400To600_PAS" , MT_Ele1MET, pileup_weight * gen_weight * weightZeroBJets ) ;
         else if(MT_Ele1MET >= 600 && MT_Ele1MET < 900)
           FillUserTH1D( "MTenu_noBtaggedJets_MT600To900_PAS" , MT_Ele1MET, pileup_weight * gen_weight * weightZeroBJets ) ;
         else if(MT_Ele1MET >= 900)
           FillUserTH1D( "MTenu_noBtaggedJets_MT900ToInf_PAS" , MT_Ele1MET, pileup_weight * gen_weight * weightZeroBJets ) ;
       }

       //-------------------------------------------------------------------------- 
       // at least 1 b tag
       //-------------------------------------------------------------------------- 
       if((isData && nBJet_medium_ptCut>=1) || !isData) {
         if ( MT_Ele1MET > 70 && MT_Ele1MET < 110 ){

         //  FillUserTH1D( "ProcessID_WWindow", ProcessID, pileup_weight * gen_weight * weightAtLeastOneBJet );
         FillUserTH1D( "MTenu_70_110_gteOneBtaggedJet"      , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //  FillUserTH1D( "nJets_MTenu_70_110_gteOneBtaggedJet", nJet_ptCut,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

         //  if ( nJet_ptCut <= 3 ){ 
         //    FillUserTH1D(   "MTenu_70_110_Njet_lte3_gteOneBtaggedJet"         , MT_Ele1MET       ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

         //    FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_lte3_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_lte3_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_lte3_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "MET_MTenu_70_110_Njet_lte3_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "sT_MTenu_70_110_Njet_lte3_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "Mej_MTenu_70_110_Njet_lte3_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //  }

         //  if ( nJet_ptCut <= 4 ){ 
         //    FillUserTH1D(   "MTenu_70_110_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

         //    FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_lte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_lte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_lte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "MET_MTenu_70_110_Njet_lte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "sT_MTenu_70_110_Njet_lte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "Mej_MTenu_70_110_Njet_lte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //  }

         //  if ( nJet_ptCut >= 4 ){ 
         //    FillUserTH1D(   "MTenu_70_110_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

         //    FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_gte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_gte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_gte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "MET_MTenu_70_110_Njet_gte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "sT_MTenu_70_110_Njet_gte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "Mej_MTenu_70_110_Njet_gte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //  }

         //  if ( nJet_ptCut >= 5 ){ 
         //    FillUserTH1D(   "MTenu_70_110_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

         //    FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_gte5_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_gte5_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_gte5_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "MET_MTenu_70_110_Njet_gte5_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "sT_MTenu_70_110_Njet_gte5_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //    FillUserTH1D(   "Mej_MTenu_70_110_Njet_gte5_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         //  }
         }

         if ( MT_Ele1MET > 50 && MT_Ele1MET < 110 ){

           FillUserTH1D( "ProcessID_WWindow", ProcessID, pileup_weight * gen_weight * weightAtLeastOneBJet );
           FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet"      , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_btagSFUpShift"        , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJetUpShift ) ;
           FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_btagSFDownShift"      , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJetDownShift ) ;
           FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet", nJet_ptCut,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           // scale factor dependence plots
           if(sT_enujj > 300 && sT_enujj < 500)
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_sT300To500_PAS"		         , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ); 
           else if(sT_enujj > 500 && sT_enujj < 750)
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_sT500To750_PAS"		         , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ); 
           else if(sT_enujj > 750 && sT_enujj < 1250)
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_sT750To1250_PAS"		         , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ); 
           else if(sT_enujj > 1250)
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_sT1250ToInf_PAS"		         , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ); 
           if(Mej > 100 && Mej < 200)
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_Mej100To200_PAS"		     , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ); 
           else if(Mej > 200 && Mej < 300)
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_Mej200To300_PAS"		     , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ); 
           else if(Mej > 300 && Mej < 400)
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_Mej300To400_PAS"		     , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ); 
           else if(Mej > 400 && Mej < 500)
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_Mej400To500_PAS"		     , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ); 
           else if(Mej > 500 && Mej < 650)
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_Mej500To650_PAS"		     , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ); 
           else if(Mej > 650)
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_Mej650ToInf_PAS"		     , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ); 

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_lte3_gteOneBtaggedJet"         , MT_Ele1MET       ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_lte3_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_lte3_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_lte3_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_lte3_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_lte3_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_lte3_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_lte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_lte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_lte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_lte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_lte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_lte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_gte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_gte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_gte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_gte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_gte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_gte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_gte5_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_gte5_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_gte5_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_gte5_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_gte5_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_gte5_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           }
         }

         if ( MT_Ele1MET > 70 && MT_Ele1MET < 150 ){

           FillUserTH1D( "ProcessID_WWindow", ProcessID, pileup_weight * gen_weight * weightAtLeastOneBJet );
           FillUserTH1D( "MTenu_70_150_gteOneBtaggedJet"      , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           FillUserTH1D( "nJets_MTenu_70_150_gteOneBtaggedJet", nJet_ptCut,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_lte3_gteOneBtaggedJet"         , MT_Ele1MET       ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_lte3_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_lte3_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_lte3_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_lte3_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_lte3_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_lte3_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_lte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_lte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_lte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_lte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_lte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_lte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_gte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_gte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_gte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_gte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_gte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_gte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_gte5_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_gte5_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_gte5_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_gte5_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_gte5_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_gte5_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           }
         }

         if ( MT_Ele1MET > 110 && MT_Ele1MET < 190 ){

           FillUserTH1D( "MTenu_110_190_gteOneBtaggedJet"      , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           FillUserTH1D( "nJets_MTenu_110_190_gteOneBtaggedJet", nJet_ptCut,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_110_190_Njet_lte3_gteOneBtaggedJet"         , MT_Ele1MET       ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_110_190_Njet_lte3_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_110_190_Njet_lte3_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_110_190_Njet_lte3_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "MET_MTenu_110_190_Njet_lte3_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "sT_MTenu_110_190_Njet_lte3_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Mej_MTenu_110_190_Njet_lte3_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_110_190_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_110_190_Njet_lte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_110_190_Njet_lte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_110_190_Njet_lte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "MET_MTenu_110_190_Njet_lte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "sT_MTenu_110_190_Njet_lte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Mej_MTenu_110_190_Njet_lte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_110_190_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_110_190_Njet_gte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_110_190_Njet_gte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_110_190_Njet_gte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "MET_MTenu_110_190_Njet_gte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "sT_MTenu_110_190_Njet_gte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Mej_MTenu_110_190_Njet_gte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_110_190_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_110_190_Njet_gte5_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_110_190_Njet_gte5_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_110_190_Njet_gte5_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "MET_MTenu_110_190_Njet_gte5_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "sT_MTenu_110_190_Njet_gte5_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
             FillUserTH1D(   "Mej_MTenu_110_190_Njet_gte5_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
           }
         }

         if(MT_Ele1MET >= 200 && MT_Ele1MET < 400)
           FillUserTH1D( "MTenu_gteOneBtaggedJet_MT200To400_PAS" , MT_Ele1MET, pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         else if(MT_Ele1MET >= 400 && MT_Ele1MET < 600)
           FillUserTH1D( "MTenu_gteOneBtaggedJet_MT400To600_PAS" , MT_Ele1MET, pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         else if(MT_Ele1MET >= 600 && MT_Ele1MET < 900)
           FillUserTH1D( "MTenu_gteOneBtaggedJet_MT600To900_PAS" , MT_Ele1MET, pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
         else if(MT_Ele1MET >= 900)
           FillUserTH1D( "MTenu_gteOneBtaggedJet_MT900ToInf_PAS" , MT_Ele1MET, pileup_weight * gen_weight * weightAtLeastOneBJet ) ;
       }

       //-------------------------------------------------------------------------- 
       // dijet mass control regions
       //-------------------------------------------------------------------------- 
       if(MT_Ele1MET > 50 && MT_Ele1MET < 110 ) {
         if(M_j1j2 > 50 && M_j1j2 < 110) {
           FillUserTH1D( "MTenu_50_110_Mjj50to110", MT_Ele1MET, pileup_weight * gen_weight);
           FillUserTH1D( "nJets_MTenu_50_110_Mjj50to110", nJet_ptCut, pileup_weight * gen_weight);
           if((isData && nBJet_medium_ptCut_beyondLeadingTwo>=1) || !isData) {
             FillUserTH1D( "MTenu_50_110_Mjj50to110_addBtagJet", MT_Ele1MET, pileup_weight * gen_weight * weightAtLeastOneBJetBeyondLeadingTwo );
             FillUserTH1D( "nJets_MTenu_50_110_Mjj50to110_addBtagJet", nJet_ptCut, pileup_weight * gen_weight * weightAtLeastOneBJetBeyondLeadingTwo );
           }
           if((isData && nBJet_medium_ptCut_beyondLeadingTwo<1) || !isData) {
             FillUserTH1D( "MTenu_50_110_Mjj50to110_noAddBtagJets", MT_Ele1MET, pileup_weight * gen_weight * weightZeroBJetsBeyondLeadingTwo );
             FillUserTH1D( "nJets_MTenu_50_110_Mjj50to110_noAddBtagJets", nJet_ptCut, pileup_weight * gen_weight * weightZeroBJetsBeyondLeadingTwo );
           }
         }
         else if(M_j1j2 > 110) {
           FillUserTH1D( "MTenu_50_110_MjjGte110", MT_Ele1MET, pileup_weight * gen_weight);
           FillUserTH1D( "nJets_MTenu_50_110_MjjGte110", nJet_ptCut, pileup_weight * gen_weight);
         }
       }

       //-------------------------------------------------------------------------- 
       // Final selection plots
       //-------------------------------------------------------------------------- 
       if ( doFinalSelections ) { 

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS         [i_lq_mass];
           bool pass    = passed_st_vector[i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_StLQ%d"       , lq_mass ); FillUserTH1D ( plot_name, sT_enujj         ,  pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_StLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, Mej              ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MTenu_StLQ%d"    , lq_mass ); FillUserTH1D ( plot_name, MT_Ele1MET       ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MET_StLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, PFMET_Type1XY_Pt,  pileup_weight * gen_weight ); 
         }

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS         [i_lq_mass];
           bool pass    = passed_mt_vector[i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_MtLQ%d"       , lq_mass ); FillUserTH1D ( plot_name, sT_enujj         ,  pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_MtLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, Mej              ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MTenu_MtLQ%d"    , lq_mass ); FillUserTH1D ( plot_name, MT_Ele1MET       ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MET_MtLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, PFMET_Type1XY_Pt,  pileup_weight * gen_weight ); 
         }

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS         [i_lq_mass];
           bool pass    = passed_mej_vector[i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_MejLQ%d"       , lq_mass ); FillUserTH1D ( plot_name, sT_enujj         ,  pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_MejLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, Mej              ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MTenu_MejLQ%d"    , lq_mass ); FillUserTH1D ( plot_name, MT_Ele1MET       ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MET_MejLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, PFMET_Type1XY_Pt,  pileup_weight * gen_weight ); 
         }

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS         [i_lq_mass];
           bool pass    = passed_met_vector[i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_MetLQ%d"       , lq_mass ); FillUserTH1D ( plot_name, sT_enujj         ,  pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_MetLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, Mej              ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MTenu_MetLQ%d"    , lq_mass ); FillUserTH1D ( plot_name, MT_Ele1MET       ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MET_MetLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, PFMET_Type1XY_Pt,  pileup_weight * gen_weight ); 
         }

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS         [i_lq_mass];
           bool pass    = passed_st_vector[i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_StLQ%d"       , lq_mass ); FillUserTH1D ( plot_name, sT_enujj         ,  pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_StLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, Mej              ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MTenu_StLQ%d"    , lq_mass ); FillUserTH1D ( plot_name, MT_Ele1MET       ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MET_StLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, PFMET_Type1XY_Pt,  pileup_weight * gen_weight ); 
         }

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS                     [i_lq_mass];
           bool pass    = passed_stAndMtAndMet_vector [i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_stAndMtAndMetLQ%d"       , lq_mass ); FillUserTH1D ( plot_name, sT_enujj         ,  pileup_weight * gen_weight);
           sprintf(plot_name, "Mej_stAndMtAndMetLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, Mej              ,  pileup_weight * gen_weight);
           sprintf(plot_name, "MTenu_stAndMtAndMetLQ%d"    , lq_mass ); FillUserTH1D ( plot_name, MT_Ele1MET       ,  pileup_weight * gen_weight);
           sprintf(plot_name, "MET_stAndMtAndMetLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, PFMET_Type1XY_Pt,  pileup_weight * gen_weight); 
         }

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS                     [i_lq_mass];
           bool pass    = passed_stAndMtAndMej_vector [i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_stAndMtAndMejLQ%d"       , lq_mass ); FillUserTH1D ( plot_name, sT_enujj         ,  pileup_weight * gen_weight);
           sprintf(plot_name, "Mej_stAndMtAndMejLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, Mej              ,  pileup_weight * gen_weight);
           sprintf(plot_name, "MTenu_stAndMtAndMejLQ%d"    , lq_mass ); FillUserTH1D ( plot_name, MT_Ele1MET       ,  pileup_weight * gen_weight);
           sprintf(plot_name, "MET_stAndMtAndMejLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, PFMET_Type1XY_Pt,  pileup_weight * gen_weight); 
         }

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS                     [i_lq_mass];
           bool pass    = passed_stAndMetAndMej_vector[i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_stAndMetAndMejLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, sT_enujj         ,  pileup_weight * gen_weight);
           sprintf(plot_name, "Mej_stAndMetAndMejLQ%d"     , lq_mass ); FillUserTH1D ( plot_name, Mej              ,  pileup_weight * gen_weight);
           sprintf(plot_name, "MTenu_stAndMetAndMejLQ%d"   , lq_mass ); FillUserTH1D ( plot_name, MT_Ele1MET       ,  pileup_weight * gen_weight);
           sprintf(plot_name, "MET_stAndMetAndMejLQ%d"     , lq_mass ); FillUserTH1D ( plot_name, PFMET_Type1XY_Pt,  pileup_weight * gen_weight); 
         }

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS                     [i_lq_mass];
           bool pass    = passed_mtAndMetAndMej_vector[i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_mtAndMetAndMejLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, sT_enujj         ,  pileup_weight * gen_weight);
           sprintf(plot_name, "Mej_mtAndMetAndMejLQ%d"     , lq_mass ); FillUserTH1D ( plot_name, Mej              ,  pileup_weight * gen_weight);
           sprintf(plot_name, "MTenu_mtAndMetAndMejLQ%d"   , lq_mass ); FillUserTH1D ( plot_name, MT_Ele1MET       ,  pileup_weight * gen_weight);
           sprintf(plot_name, "MET_mtAndMetAndMejLQ%d"     , lq_mass ); FillUserTH1D ( plot_name, PFMET_Type1XY_Pt,  pileup_weight * gen_weight); 
         }

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS      [i_lq_mass];
           bool pass = passed_stAndMej_vector [i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_StAndMejLQ%d"       , lq_mass ); FillUserTH1D ( plot_name, sT_enujj         ,  pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_StAndMejLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, Mej              ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MTenu_StAndMejLQ%d"    , lq_mass ); FillUserTH1D ( plot_name, MT_Ele1MET       ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MET_StAndMejLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, PFMET_Type1XY_Pt,  pileup_weight * gen_weight ); 
         }

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS      [i_lq_mass];
           bool pass = passed_mtAndMej_vector [i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_MtAndMejLQ%d"       , lq_mass ); FillUserTH1D ( plot_name, sT_enujj         ,  pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_MtAndMejLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, Mej              ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MTenu_MtAndMejLQ%d"    , lq_mass ); FillUserTH1D ( plot_name, MT_Ele1MET       ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MET_MtAndMejLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, PFMET_Type1XY_Pt,  pileup_weight * gen_weight ); 
         }

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS      [i_lq_mass];
           bool pass = passed_metAndMej_vector[i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_MetAndMejLQ%d"       , lq_mass ); FillUserTH1D ( plot_name, sT_enujj         ,  pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_MetAndMejLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, Mej              ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MTenu_MetAndMejLQ%d"    , lq_mass ); FillUserTH1D ( plot_name, MT_Ele1MET       ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MET_MetAndMejLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, PFMET_Type1XY_Pt,  pileup_weight * gen_weight ); 
         }

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS      [i_lq_mass];
           bool pass = passed_stAndMet_vector [i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_StAndMetLQ%d"       , lq_mass ); FillUserTH1D ( plot_name, sT_enujj         ,  pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_StAndMetLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, Mej              ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MTenu_StAndMetLQ%d"    , lq_mass ); FillUserTH1D ( plot_name, MT_Ele1MET       ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MET_StAndMetLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, PFMET_Type1XY_Pt,  pileup_weight * gen_weight ); 
         }

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS      [i_lq_mass];
           bool pass = passed_mtAndMet_vector [i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_MtAndMetLQ%d"       , lq_mass ); FillUserTH1D ( plot_name, sT_enujj         ,  pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_MtAndMetLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, Mej              ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MTenu_MtAndMetLQ%d"    , lq_mass ); FillUserTH1D ( plot_name, MT_Ele1MET       ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MET_MtAndMetLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, PFMET_Type1XY_Pt,  pileup_weight * gen_weight ); 
         }

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS      [i_lq_mass];
           bool pass = passed_stAndMt_vector  [i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_StAndMtLQ%d"       , lq_mass ); FillUserTH1D ( plot_name, sT_enujj         ,  pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_StAndMtLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, Mej              ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MTenu_StAndMtLQ%d"    , lq_mass ); FillUserTH1D ( plot_name, MT_Ele1MET       ,  pileup_weight * gen_weight );
           sprintf(plot_name, "MET_StAndMtLQ%d"      , lq_mass ); FillUserTH1D ( plot_name, PFMET_Type1XY_Pt,  pileup_weight * gen_weight ); 
         }

         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS      [i_lq_mass];
           bool pass    = passed_vector[i_lq_mass];
           if ( !pass ) continue;

           // if ( lq_mass == 650 && isData > 0.5 ) { 
           //   std::string root_file_name( fChain -> GetCurrentFile() -> GetName() );
           //   std::size_t found = root_file_name.find ( std::string ("ElectronHad"));
           //   if ( found != std::string::npos ) {
           //     ofstream myfile;
           //     myfile.open ( file_name.c_str() ,  ios::out | ios::app );
           //     myfile << fixed << int(run) << ":" << int(ls) << ":" << int(event) << "\n";
           //     myfile.close();
           //   }
           // }
           //if(lq_mass==1150)
           //  std::cout << "[LQ1150] passing run/ls/event: " << static_cast<int>(run) << " " << static_cast<int>(ls) << " " << ((unsigned int)event) << 
           //    ": sT=" << sT_enujj << "; PtEle1=" << Ele1_Pt << "; EtSCEle1=" << Ele1_SCEt << "; MET=" << PFMET_Type1XY_Pt << "; PT(e,MET)=" << Pt_Ele1MET <<
           //    std::endl;

           sprintf(plot_name, "Ele_EtaVsPhi_LQ%d" , lq_mass ); FillUserTH2D( plot_name , Ele1_Eta, Ele1_Phi, gen_weight * pileup_weight );
           sprintf(plot_name, "Jet_EtaVsPhi_LQ%d" , lq_mass ); FillUserTH2D( plot_name , Jet1_Eta, Jet1_Phi, gen_weight * pileup_weight );
           sprintf(plot_name, "Jet_EtaVsPhi_LQ%d" , lq_mass ); FillUserTH2D( plot_name , Jet2_Eta, Jet2_Phi, gen_weight * pileup_weight );
           sprintf(plot_name, "Jet1_EtaVsPhi_LQ%d", lq_mass ); FillUserTH2D( plot_name , Jet1_Eta, Jet1_Phi, gen_weight * pileup_weight );
           sprintf(plot_name, "Jet2_EtaVsPhi_LQ%d", lq_mass ); FillUserTH2D( plot_name , Jet2_Eta, Jet2_Phi, gen_weight * pileup_weight );

           sprintf(plot_name, "MET_LQ%d"    , lq_mass ); FillUserTH1D( plot_name , PFMET_Type1XY_Pt ,gen_weight * pileup_weight );
           //// SIC ADD
           //sprintf(plot_name, "MET_UnclUp_LQ%d"    , lq_mass ); FillUserTH1D( plot_name , PFMET_Type1XY_UnclUp_Pt ,gen_weight * pileup_weight );
           //sprintf(plot_name, "MET_UnclDown_LQ%d"    , lq_mass ); FillUserTH1D( plot_name , PFMET_Type1XY_UnclDown_Pt ,gen_weight * pileup_weight );
           sprintf(plot_name, "Mej_LQ%d"    , lq_mass ); FillUserTH1D( plot_name , Mej               ,gen_weight * pileup_weight );
           sprintf(plot_name, "Mejj_LQ%d"   , lq_mass ); FillUserTH1D( plot_name , Mejj              ,gen_weight * pileup_weight );
           sprintf(plot_name, "Mjj_LQ%d"    , lq_mass ); FillUserTH1D( plot_name , M_j1j2            ,gen_weight * pileup_weight );
           sprintf(plot_name, "ST_LQ%d"     , lq_mass ); FillUserTH1D( plot_name , sT_enujj          ,gen_weight * pileup_weight );
           sprintf(plot_name, "MTjnu_LQ%d"  , lq_mass ); FillUserTH1D( plot_name , MT_JetMET         ,gen_weight * pileup_weight );
           sprintf(plot_name, "MTenu_LQ%d"  , lq_mass ); FillUserTH1D( plot_name , MT_Ele1MET        ,gen_weight * pileup_weight );
           sprintf(plot_name, "MTjjnu_LQ%d" , lq_mass ); FillUserTH1D( plot_name , MT_jjnu           ,gen_weight * pileup_weight );
           sprintf(plot_name, "MTejjnu_LQ%d", lq_mass ); FillUserTH1D( plot_name , MT_ejjnu          ,gen_weight * pileup_weight );
           sprintf(plot_name, "nJets_LQ%d"  , lq_mass ); FillUserTH1D( plot_name , nJet_ptCut        ,gen_weight * pileup_weight );
           sprintf(plot_name, "mDPhi_Jet1Jet2_LQ%d", lq_mass ); FillUserTH1D( plot_name, mDPhi_Jet1Jet2,gen_weight * pileup_weight );
           sprintf(plot_name, "mDPhi_Ele1MET_LQ%d" , lq_mass ); FillUserTH1D( plot_name, mDPhi_METEle1 ,gen_weight * pileup_weight );
           sprintf(plot_name, "nBJet_loose_LQ%d" , lq_mass ); FillUserTH1D( plot_name , nBJet_loose_ptCut , pileup_weight * gen_weight );
           sprintf(plot_name, "nBJet_medium_LQ%d", lq_mass ); FillUserTH1D( plot_name , nBJet_medium_ptCut, pileup_weight * gen_weight );
           sprintf(plot_name, "nBJet_tight_LQ%d" , lq_mass ); FillUserTH1D( plot_name , nBJet_tight_ptCut , pileup_weight * gen_weight );
           if ( nJet_store > 2 ) 
             sprintf(plot_name, "Mej3_LQ%d"  , lq_mass ); FillUserTH1D( plot_name , Mej3              , gen_weight * pileup_weight );


           if ( do_extra_finalSelection_plots ) { 
             sprintf(plot_name, "Me1j1_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , M_e1j1                         , pileup_weight * gen_weight );
             sprintf(plot_name, "Me1j2_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , M_e1j2                         , pileup_weight * gen_weight );
             sprintf(plot_name, "Eta1stJet_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Jet1_Eta                       , pileup_weight * gen_weight );
             sprintf(plot_name, "Eta2ndJet_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Jet2_Eta                       , pileup_weight * gen_weight );
             sprintf(plot_name, "Eta1stEle_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Ele1_Eta                       , pileup_weight * gen_weight );
             sprintf(plot_name, "Charge1stEle_LQ%d"      , lq_mass ); FillUserTH1D( plot_name , Ele1_Charge                    , pileup_weight * gen_weight );
             sprintf(plot_name, "Phi1stEle_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Ele1_Phi                       , pileup_weight * gen_weight );
             sprintf(plot_name, "Phi1stJet_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Jet1_Phi                       , pileup_weight * gen_weight );
             sprintf(plot_name, "Phi2ndJet_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Jet2_Phi                       , pileup_weight * gen_weight );
             sprintf(plot_name, "nVertex_LQ%d"           , lq_mass ); FillUserTH1D( plot_name , nVertex                        , pileup_weight * gen_weight );
             sprintf(plot_name, "sTlep_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , Ele1_Pt + PFMET_Type1XY_Pt     , pileup_weight * gen_weight );
             sprintf(plot_name, "sTjet_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , Jet1_Pt + Jet2_Pt              , pileup_weight * gen_weight );
             sprintf(plot_name, "Energy1stEle_LQ%d"      , lq_mass ); FillUserTH1D( plot_name , Ele1_Energy                    , pileup_weight * gen_weight );
             sprintf(plot_name, "Pt1stEle_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Ele1_Pt                        , pileup_weight * gen_weight );
             sprintf(plot_name, "PtHeep1stEle_LQ%d"      , lq_mass ); FillUserTH1D( plot_name , Ele1_PtHeep                    , pileup_weight * gen_weight );
             sprintf(plot_name, "SCEt1stEle_LQ%d"        , lq_mass ); FillUserTH1D( plot_name , Ele1_SCEnergy/cosh(Ele1_SCEta) , pileup_weight * gen_weight );
             sprintf(plot_name, "SCEta1stEle_LQ%d"       , lq_mass ); FillUserTH1D( plot_name , Ele1_SCEta                     , pileup_weight * gen_weight );
             sprintf(plot_name, "Pt1stJet_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Jet1_Pt                        , pileup_weight * gen_weight );
             sprintf(plot_name, "Pt2ndJet_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Jet2_Pt                        , pileup_weight * gen_weight );
             sprintf(plot_name, "sTfrac_Jet1_LQ%d"       , lq_mass ); FillUserTH1D( plot_name , Jet1_Pt / sT_enujj              , pileup_weight * gen_weight );
             sprintf(plot_name, "sTfrac_Jet2_LQ%d"       , lq_mass ); FillUserTH1D( plot_name , Jet2_Pt / sT_enujj              , pileup_weight * gen_weight );
             sprintf(plot_name, "sTfrac_Ele1_LQ%d"       , lq_mass ); FillUserTH1D( plot_name , Ele1_Pt / sT_enujj              , pileup_weight * gen_weight );
             sprintf(plot_name, "sTfrac_MET_LQ%d"        , lq_mass ); FillUserTH1D( plot_name , PFMET_Type1XY_Pt / sT_enujj              , pileup_weight * gen_weight );
             sprintf(plot_name, "sTfrac_Jet_LQ%d"        , lq_mass ); FillUserTH1D( plot_name , ( Jet1_Pt + Jet2_Pt ) / sT_enujj, pileup_weight * gen_weight );
             sprintf(plot_name, "sTfrac_Lep_LQ%d"        , lq_mass ); FillUserTH1D( plot_name , ( Ele1_Pt + PFMET_Type1XY_Pt ) / sT_enujj, pileup_weight * gen_weight );
             // muon kinematics
             sprintf(plot_name, "Eta1stMuon_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Muon1_Eta                       , pileup_weight * gen_weight );
             sprintf(plot_name, "Phi1stMuon_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Muon1_Phi                       , pileup_weight * gen_weight );
             sprintf(plot_name, "Pt1stMuon_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Muon1_Pt                        , pileup_weight * gen_weight );
             sprintf(plot_name, "Eta2ndMuon_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Muon2_Eta                       , pileup_weight * gen_weight );
             sprintf(plot_name, "Phi2ndMuon_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Muon2_Phi                       , pileup_weight * gen_weight );
             sprintf(plot_name, "Pt2ndMuon_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Muon2_Pt                        , pileup_weight * gen_weight );
           }

           if ( do_eleQuality_plots ) { 
             sprintf(plot_name, "BeamSpotDXY_1stEle_LQ%d"        , lq_mass );   FillUserTH1D(plot_name,  Ele1_BeamSpotDXY               , pileup_weight * gen_weight ); 
             sprintf(plot_name, "Classif_1stEle_LQ%d"            , lq_mass );   FillUserTH1D(plot_name,  Ele1_Classif                   , pileup_weight * gen_weight ); 
             sprintf(plot_name, "CorrIsolation_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele1_CorrIsolation             , pileup_weight * gen_weight ); 
             sprintf(plot_name, "DeltaEtaTrkSC_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele1_DeltaEtaTrkSC             , pileup_weight * gen_weight ); 
             sprintf(plot_name, "DeltaPhiTrkSC_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele1_DeltaPhiTrkSC             , pileup_weight * gen_weight ); 
             sprintf(plot_name, "Full5x5E1x5OverE5x5_1stEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  Ele1_Full5x5E1x5OverE5x5              , pileup_weight * gen_weight ); 
             sprintf(plot_name, "Full5x5E2x5OverE5x5_1stEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  Ele1_Full5x5E2x5OverE5x5              , pileup_weight * gen_weight ); 
             sprintf(plot_name, "EcalIsolation_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele1_EcalIsolation             , pileup_weight * gen_weight ); 
             sprintf(plot_name, "HcalIsolation_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele1_HcalIsolation             , pileup_weight * gen_weight ); 
             sprintf(plot_name, "TrkIsolation_1stEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  Ele1_TrkIsolation              , pileup_weight * gen_weight ); 
             sprintf(plot_name, "Energy_1stEle_LQ%d"             , lq_mass );   FillUserTH1D(plot_name,  Ele1_Energy                    , pileup_weight * gen_weight ); 
             sprintf(plot_name, "FBrem_1stEle_LQ%d"              , lq_mass );   FillUserTH1D(plot_name,  Ele1_FBrem                     , pileup_weight * gen_weight ); 
             sprintf(plot_name, "GsfCtfCharge_1stEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  Ele1_GsfCtfCharge              , pileup_weight * gen_weight ); 
             sprintf(plot_name, "GsfCtfScPixCharge_1stEle_LQ%d"  , lq_mass );   FillUserTH1D(plot_name,  Ele1_GsfCtfScPixCharge         , pileup_weight * gen_weight ); 
             sprintf(plot_name, "GsfScPixCharge_1stEle_LQ%d"     , lq_mass );   FillUserTH1D(plot_name,  Ele1_GsfScPixCharge            , pileup_weight * gen_weight ); 
             sprintf(plot_name, "HasMatchedPhot_1stEle_LQ%d"     , lq_mass );   FillUserTH1D(plot_name,  Ele1_HasMatchedPhot            , pileup_weight * gen_weight ); 
             sprintf(plot_name, "HoE_1stEle_LQ%d"                , lq_mass );   FillUserTH1D(plot_name,  Ele1_HoE                       , pileup_weight * gen_weight ); 
             sprintf(plot_name, "LeadVtxDistXY_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele1_LeadVtxDistXY             , pileup_weight * gen_weight ); 
             sprintf(plot_name, "LeadVtxDistZ_1stEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  Ele1_LeadVtxDistZ              , pileup_weight * gen_weight ); 
             sprintf(plot_name, "MissingHits_1stEle_LQ%d"        , lq_mass );   FillUserTH1D(plot_name,  Ele1_MissingHits               , pileup_weight * gen_weight ); 
             sprintf(plot_name, "NBrems_1stEle_LQ%d"             , lq_mass );   FillUserTH1D(plot_name,  Ele1_NBrems                    , pileup_weight * gen_weight ); 
             sprintf(plot_name, "ValidFrac_1stEle_LQ%d"          , lq_mass );   FillUserTH1D(plot_name,  Ele1_ValidFrac                 , pileup_weight * gen_weight ); 
             sprintf(plot_name, "EnergyORawEnergy_1stEle_LQ%d"   , lq_mass );   FillUserTH1D(plot_name,  Ele1_Energy / Ele1_RawEnergy   , pileup_weight * gen_weight ); 
             sprintf(plot_name, "TrkPtOPt_1stEle_LQ%d"           , lq_mass );   FillUserTH1D(plot_name,  Ele1_TrkPt  / Ele1_Pt          , pileup_weight * gen_weight ); 
             sprintf(plot_name, "EOverP_1stEle_LQ%d"             , lq_mass );   FillUserTH1D(plot_name,  Ele1_EOverP                    , pileup_weight * gen_weight ); 

             if ( fabs(Ele1_Eta) < eleEta_bar ) { 
               sprintf(plot_name, "SigmaEtaEta_Barrel_1stEle_LQ%d"  , lq_mass ); FillUserTH1D( plot_name , Ele1_SigmaEtaEta   , pileup_weight * gen_weight    ); 
               sprintf(plot_name, "Full5x5SigmaIEtaIEta_Barrel_1stEle_LQ%d", lq_mass ); FillUserTH1D( plot_name , Ele1_Full5x5SigmaIEtaIEta , pileup_weight * gen_weight    ); 
             }
             else if ( fabs(Ele1_Eta) > eleEta_end1_min && fabs(Ele1_Eta) < eleEta_end2_max ){
               sprintf(plot_name, "SigmaEtaEta_Endcap_1stEle_LQ%d"  , lq_mass ); FillUserTH1D( plot_name , Ele1_SigmaEtaEta   , pileup_weight * gen_weight    ); 
               sprintf(plot_name, "Full5x5SigmaIEtaIEta_Endcap_1stEle_LQ%d", lq_mass ); FillUserTH1D( plot_name , Ele1_Full5x5SigmaIEtaIEta , pileup_weight * gen_weight    ); 
             }

             if ( fabs(Ele1_Eta) <= eleEta_bar ) { 
               sprintf(plot_name,"Mej_Barrel_LQ%d", lq_mass);  FillUserTH1D( plot_name, Mej, pileup_weight * gen_weight);
             }
             else if ( fabs(Ele1_Eta) >= eleEta_end1_min && fabs(Ele1_Eta) < eleEta_end1_max) {
               sprintf(plot_name,"Mej_Endcap1_LQ%d", lq_mass);  FillUserTH1D( plot_name, Mej, pileup_weight * gen_weight);
             }
             else if ( fabs(Ele1_Eta) >= eleEta_end2_min && fabs(Ele1_Eta) < eleEta_end2_max) {
               sprintf(plot_name,"Mej_Endcap2_LQ%d", lq_mass);  FillUserTH1D( plot_name, Mej, pileup_weight * gen_weight);
             }
           }

           //if ( isData == 1 ) {
           //  sprintf(plot_name, "split_1fb_LQ%d", lq_mass ); 
           //  int split_1fb = get_split_1fb ( run );
           //  FillUserTH1D( plot_name, split_1fb , pileup_weight * gen_weight );
           //  if ( lq_mass == 600 || lq_mass == 650 ) { 
           //    std::cout << std::fixed << lq_mass << ", " << int(run) << ":" << int(ls) << ":" << int(event) << std::endl;
           //  }
           //}
         }

         //if ( passedCut ("ST_LQ300") && passedCut ("Mej_LQ300") ){
         //  FillUserTH1D("MET_LQ300_NoMETCut", PFMET_Type1XY_Pt , gen_weight * pileup_weight );
         //  FillUserTH1D("Mej_LQ300_NoMETCut", Mej               , gen_weight * pileup_weight );
         //  FillUserTH1D("ST_LQ300_NoMETCut" , sT_enujj          , gen_weight * pileup_weight ); 
         //  FillUserTH1D("MT_LQ300_NoMETCut" , MT_Ele1MET        , gen_weight * pileup_weight ); 
         //}
         // for scale factor at "final selection" studies
         if((isData && nBJet_medium_ptCut==0) || !isData) {
           if ( hasCut("ST_LQ300") && passedCut("ST_LQ300") && passedCut("Mej_LQ300") && passedCut("MET_LQ300") ){
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ300"       , MT_Ele1MET, pileup_weight * gen_weight * weightZeroBJets);
             FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ300" , nJet_ptCut, pileup_weight * gen_weight * weightZeroBJets);
           }
           if ( hasCut("ST_LQ400") && passedCut("ST_LQ400") && passedCut("Mej_LQ400") && passedCut("MET_LQ400") ){
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ400"       , MT_Ele1MET, pileup_weight * gen_weight * weightZeroBJets);
             FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ400" , nJet_ptCut, pileup_weight * gen_weight * weightZeroBJets);
           }
           if ( hasCut("ST_LQ500") && passedCut("ST_LQ500") && passedCut("Mej_LQ500") && passedCut("MET_LQ500") ){
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ500"       , MT_Ele1MET, pileup_weight * gen_weight * weightZeroBJets);
             FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ500" , nJet_ptCut, pileup_weight * gen_weight * weightZeroBJets);
           }
           if ( hasCut("ST_LQ600") && passedCut("ST_LQ600") && passedCut("Mej_LQ600") && passedCut("MET_LQ600") ){
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ600"       , MT_Ele1MET, pileup_weight * gen_weight * weightZeroBJets);
             FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ600" , nJet_ptCut, pileup_weight * gen_weight * weightZeroBJets);
           }
           if ( hasCut("ST_LQ700") && passedCut("ST_LQ700") && passedCut("Mej_LQ700") && passedCut("MET_LQ700") ){
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ700"       , MT_Ele1MET, pileup_weight * gen_weight * weightZeroBJets);
             FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ700" , nJet_ptCut, pileup_weight * gen_weight * weightZeroBJets);
           }
           if ( hasCut("ST_LQ800") && passedCut("ST_LQ800") && passedCut("Mej_LQ800") && passedCut("MET_LQ800") ){
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ800"          , MT_Ele1MET, pileup_weight * gen_weight * weightZeroBJets);
             FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ800"    , nJet_ptCut, pileup_weight * gen_weight * weightZeroBJets);
           }
           if ( hasCut("ST_LQ900") && passedCut("ST_LQ900") && passedCut("Mej_LQ900") && passedCut("MET_LQ900") ){
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ900"          , MT_Ele1MET, pileup_weight * gen_weight * weightZeroBJets);
             FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ900"    , nJet_ptCut, pileup_weight * gen_weight * weightZeroBJets);
           }
           if ( hasCut("ST_LQ1000") && passedCut("ST_LQ1000") && passedCut("Mej_LQ1000") && passedCut("MET_LQ1000") ){
           }
           FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ1000"       ,MT_Ele1MET, pileup_weight * gen_weight * weightZeroBJets); 
           FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ1000" ,nJet_ptCut, pileup_weight * gen_weight * weightZeroBJets);
         }
         if((isData && nBJet_medium_ptCut>=1) || !isData) {
           if ( hasCut("ST_LQ300") && passedCut("ST_LQ300") && passedCut("Mej_LQ300") && passedCut("MET_LQ300") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ300"       , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ300" , nJet_ptCut,  pileup_weight * gen_weight * weightAtLeastOneBJet );
           }
           if ( hasCut("ST_LQ400") && passedCut("ST_LQ400") && passedCut("Mej_LQ400") && passedCut("MET_LQ400") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ400"       , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ400" , nJet_ptCut,  pileup_weight * gen_weight * weightAtLeastOneBJet );
           }
           if ( hasCut("ST_LQ500") && passedCut("ST_LQ500") && passedCut("Mej_LQ500") && passedCut("MET_LQ500") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ500"       , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ500" , nJet_ptCut,  pileup_weight * gen_weight * weightAtLeastOneBJet );
           }
           if ( hasCut("ST_LQ600") && passedCut("ST_LQ600") && passedCut("Mej_LQ600") && passedCut("MET_LQ600") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ600"       , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ600" , nJet_ptCut,  pileup_weight * gen_weight * weightAtLeastOneBJet );
           }
           if ( hasCut("ST_LQ700") && passedCut("ST_LQ700") && passedCut("Mej_LQ700") && passedCut("MET_LQ700") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ700"       , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ700" , nJet_ptCut,  pileup_weight * gen_weight * weightAtLeastOneBJet );
           }
           if ( hasCut("ST_LQ800") && passedCut("ST_LQ800") && passedCut("Mej_LQ800") && passedCut("MET_LQ800") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ800"       , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ800" , nJet_ptCut,  pileup_weight * gen_weight * weightAtLeastOneBJet );
           }
           if ( hasCut("ST_LQ900") && passedCut("ST_LQ900") && passedCut("Mej_LQ900") && passedCut("MET_LQ900") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ900"       , MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ900" , nJet_ptCut,  pileup_weight * gen_weight * weightAtLeastOneBJet );
           }
           if ( hasCut("ST_LQ1000") && passedCut("ST_LQ1000") && passedCut("Mej_LQ1000") && passedCut("MET_LQ1000") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ1000"       ,MT_Ele1MET,  pileup_weight * gen_weight * weightAtLeastOneBJet );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ1000" ,nJet_ptCut,  pileup_weight * gen_weight * weightAtLeastOneBJet );
           }
         }
       } // End do final selection
     } // End do preselection
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
