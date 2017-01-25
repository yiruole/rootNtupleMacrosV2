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
// for trigger turn-on
#include "Ele27WPLooseTrigTurnOn.C"
// for fake rate
#include "include/QCDFakeRate.h"


analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile){}

analysisClass::~analysisClass(){}

void analysisClass::Loop()
{
   std::cout << "analysisClass::Loop() begins" <<std::endl;   
   bool do_final_selection = true;
   bool do_eleQuality_plots = true;
   bool do_extra_finalSelection_plots = true;

   char cut_name[100];
   char st_cut_name [100];
   char mt_cut_name [100];
   char mej_cut_name[100];
   char met_cut_name[100];   

   //--------------------------------------------------------------------------
   // MET weight function
   //--------------------------------------------------------------------------
   
   TF1 * met_weight_func = new TF1( "met_func", "pol1", 0, 5000 );
   met_weight_func -> SetParameter(0, 0.988719);
   met_weight_func -> SetParameter(1, -0.000967365);
   
   TF1 * mt_weight_func = new TF1 ("mt_func", "gaus(0)+pol1(3)", 0, 5000 );
   mt_weight_func -> SetParameter(0, 0.105337    );
   mt_weight_func -> SetParameter(1, 162.649     );
   mt_weight_func -> SetParameter(2, 38.3520     );
   mt_weight_func -> SetParameter(3, 0.942490    );
   mt_weight_func -> SetParameter(4, 0.000374290 );
   
   TF1 * mt_weight_func_first = new TF1 ("mt_func_first", "gaus(0)+pol1(3)", 0, 5000 );
   mt_weight_func_first -> SetParameter(0, 0.136393     );
   mt_weight_func_first -> SetParameter(1, 158.679	);
   mt_weight_func_first -> SetParameter(2, 40.6906	);
   mt_weight_func_first -> SetParameter(3, 0.854018	);
   mt_weight_func_first -> SetParameter(4, 0.0000162925 );

   TF1 * mt_weight_func_muon_background = new TF1 ("mt_func_muon_background", "gaus(0)+pol1(3)", 0, 5000 );
   mt_weight_func_muon_background -> SetParameter(0, .363729);
   mt_weight_func_muon_background -> SetParameter(1, 139.421);
   mt_weight_func_muon_background -> SetParameter(2, 47.7152);
   mt_weight_func_muon_background -> SetParameter(3, .772203);
   mt_weight_func_muon_background -> SetParameter(4, .000326885);
   
   //--------------------------------------------------------------------------
   // Event list
   //--------------------------------------------------------------------------

   std::string file_name ("enujj_event_list.txt");
   
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

   //// SIC only look at LQ650 selection for now
   //// LQ650 only 2012
   //const int n_lq_mass = 1;
   //int LQ_MASS[n_lq_mass] = { 
   //   6502012
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
   
   std::vector<std::string> signals;
   signals.push_back ( std::string("LQ_M300"));
   signals.push_back ( std::string("LQ_M350"));
   signals.push_back ( std::string("LQ_M400"));
   signals.push_back ( std::string("LQ_M450"));
   signals.push_back ( std::string("LQ_M500"));
   signals.push_back ( std::string("LQ_M550"));
   signals.push_back ( std::string("LQ_M600"));
   signals.push_back ( std::string("LQ_M650"));
   signals.push_back ( std::string("LQ_M700"));
   signals.push_back ( std::string("LQ_M750"));
   signals.push_back ( std::string("LQ_M800"));
   signals.push_back ( std::string("LQ_M850"));
   signals.push_back ( std::string("LQ_M900"));
   signals.push_back ( std::string("LQ_M950"));
   signals.push_back ( std::string("LQ_M1000"));
   signals.push_back ( std::string("LQ_M1050"));
   signals.push_back ( std::string("LQ_M1100"));
   signals.push_back ( std::string("LQ_M1150"));
   signals.push_back ( std::string("LQ_M1200"));

   std::string pdf_file_name("/afs/cern.ch/user/e/eberry/public/LQPDF/enujj_pdfs_withMT.root");
   
   likelihoodGetter * myLikelihoodGetter = new likelihoodGetter ( pdf_file_name, variables, signals ) ;
   
   std::vector<double> likelihood_values;
   
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

   double eleEta_bar_max = getPreCutValue1("eleEta_bar");
   double eleEta_end_min = getPreCutValue1("eleEta_end1");
   double eleEta_end_max = getPreCutValue2("eleEta_end2");


   // eta boundaries

   double eleEta_bar            = getPreCutValue1("eleEta_bar");
   double eleEta_end1_min       = getPreCutValue1("eleEta_end1");
   double eleEta_end1_max       = getPreCutValue2("eleEta_end1");
   double eleEta_end2_min       = getPreCutValue1("eleEta_end2");
   double eleEta_end2_max       = getPreCutValue2("eleEta_end2");

   // std::string test_string = getPreCutString1("test");
   // std::cout << "test string = \"" << test_string << "\"" << std::endl;

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------

   CreateUserTH1D( "ProcessID"             ,    21  , -0.5    , 20.5     );
   CreateUserTH1D( "ProcessID_PAS"         ,    21  , -0.5    , 20.5     );
   CreateUserTH1D( "ProcessID_WWindow"     ,    21  , -0.5    , 20.5     );
   

   CreateUserTH1D( "split_PAS"                ,    38    , -0.5    , 37.5     );
   CreateUserTH1D( "split_1fb_PAS"            ,    38    , -0.5    , 37.5     );

   CreateUserTH1D( "MET_PDF"                  , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PDF"                   , 500 , 0       , 5000	 ); 
   CreateUserTH1D( "MT_PDF"                   , 400 , 0       , 4000	 ); 
   CreateUserTH1D( "Mej_PDF"                  , 400 , 0       , 4000	 ); 

   CreateUserTH1D( "nElectron_PAS"            , 5   , -0.5    , 4.5      );
   CreateUserTH1D( "nMuon_PAS"                , 5   , -0.5    , 4.5      );
   CreateUserTH1D( "nJet_PAS"                 , 11  , -0.5    , 10.5     );
   CreateUserTH1D( "nJet_PASandFrancesco"     , 11  , -0.5    , 10.5     );
   CreateUserTH1D( "Pt1stEle_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "Energy1stEle_PAS"	      , 200 , 0       , 2000     ); 
   CreateUserTH1D( "Eta1stEle_PAS"	      , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stEle_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Charge1stEle_PAS"	      , 2   , -1.0001 , 1.0001	 ); 
   CreateUserTH1D( "MET_PAS"                  , 200 , 0       , 1000	 ); 
   // SIC ADD
   //CreateUserTH1D( "MET_PAS_UnclUp"                  , 200 , 0       , 1000	 ); 
   //CreateUserTH1D( "MET_PAS_UnclDown"                  , 200 , 0       , 1000	 ); 
   // END SIC ADD

   CreateUserTH1D( "METPhi_PAS"		      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "MET_Type01_PAS"           , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "MET_Type01_Phi_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "MET_Type1_PAS"           , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "MET_Type1_Phi_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "minMETPt1stEle_PAS"       , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "minMET01Pt1stEle_PAS"     , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "minMET1Pt1stEle_PAS"     , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "Pt1stJet_PAS"             , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Pt2ndJet_PAS"             , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Eta1stJet_PAS"            , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Eta2ndJet_PAS"            , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stJet_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Phi2ndJet_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Mass1stJet_PAS"           , 100 , 0       , 500      );
   CreateUserTH1D( "Mass2ndJet_PAS"           , 100 , 0       , 500      );
   CreateUserTH1D( "CSV1stJet_PAS"            , 200 , 0       , 1.0	 ); 
   CreateUserTH1D( "CSV2ndJet_PAS"            , 200 , 0       , 1.0	 ); 
   CreateUserTH1D( "nMuon_PtCut_IDISO_PAS"    , 16  , -0.5    , 15.5	 ); 
   CreateUserTH1D( "MTenu_PAS"                , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "MTenu_Type01_PAS"         , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "Ptenu_PAS"		      , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTlep_PAS"                , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTlep_Type01_PAS"         , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTlep_Type1_PAS"         , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTjet_PAS"                , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PAS"                   , 300 , 0       , 3000	 ); 
   CreateUserTH1D( "sT_ControlRegion_Njet_gte4", 300 , 0       , 3000	 ); 
   CreateUserTH1D( "sT_ControlRegion_Njet_lte3", 300 , 0       , 3000	 ); 
   CreateUserTH1D( "sT_Type01_PAS"            , 300 , 0       , 3000	 ); 
   CreateUserTH1D( "sT_Type1_PAS"            , 300 , 0       , 3000	 ); 
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
   CreateUserTH1D( "mDPhi1stEleMET_Type01_PAS", 100 , 0.      , 3.14159  );
   CreateUserTH1D( "mDPhi1stJetMET_Type01_PAS", 100 , 0.      , 3.14159  );
   CreateUserTH1D( "mDPhi2ndJetMET_Type01_PAS", 100 , 0.      , 3.14159  );
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

   CreateUserTH1D( "MTenu_Type01_50_110", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_Type01_50_110"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_Type01_50_110_Njet_gte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_50_110_Njet_lte3", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_50_110_Njet_lte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_50_110_Njet_gte5", 240, 40, 160 );

   CreateUserTH1D( "MTenu_70_110", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_70_110"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte3", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte5", 240, 40, 160 );

   CreateUserTH1D( "MTenu_Type01_70_110", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_Type01_70_110"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_Type01_70_110_Njet_gte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_110_Njet_lte3", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_110_Njet_lte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_110_Njet_gte5", 240, 40, 160 );

   CreateUserTH1D( "MTenu_70_150", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_70_150"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_70_150_Njet_gte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_lte3", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_lte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_gte5", 240, 40, 160 );

   CreateUserTH1D( "MTenu_Type01_70_150", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_Type01_70_150"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_Type01_70_150_Njet_gte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_150_Njet_lte3", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_150_Njet_lte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_150_Njet_gte5", 240, 40, 160 );

   // without btags
   CreateUserTH1D( "MTenu_50_110_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_noBtaggedJets"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte3_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte5_noBtaggedJets", 240, 40, 160 );

   CreateUserTH1D( "MTenu_Type01_50_110_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_Type01_50_110_noBtaggedJets"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_Type01_50_110_Njet_gte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_50_110_Njet_lte3_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_50_110_Njet_lte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_50_110_Njet_gte5_noBtaggedJets", 240, 40, 160 );

   CreateUserTH1D( "MTenu_70_110_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_70_110_noBtaggedJets"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte3_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte5_noBtaggedJets", 240, 40, 160 );

   CreateUserTH1D( "MTenu_Type01_70_110_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_Type01_70_110_noBtaggedJets"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_Type01_70_110_Njet_gte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_110_Njet_lte3_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_110_Njet_lte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_110_Njet_gte5_noBtaggedJets", 240, 40, 160 );

   CreateUserTH1D( "MTenu_70_150_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_70_150_noBtaggedJets"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_70_150_Njet_gte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_lte3_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_lte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_gte5_noBtaggedJets", 240, 40, 160 );

   CreateUserTH1D( "MTenu_Type01_70_150_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_Type01_70_150_noBtaggedJets"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_Type01_70_150_Njet_gte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_150_Njet_lte3_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_150_Njet_lte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_150_Njet_gte5_noBtaggedJets", 240, 40, 160 );

   // with a btag
   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte3_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte5_gteOneBtaggedJet", 240, 40, 160 );

   CreateUserTH1D( "MTenu_Type01_50_110_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_Type01_50_110_gteOneBtaggedJet"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet", 240, 40, 160 );

   CreateUserTH1D( "MTenu_70_110_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_70_110_gteOneBtaggedJet"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte3_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte5_gteOneBtaggedJet", 240, 40, 160 );

   CreateUserTH1D( "MTenu_Type01_70_110_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_Type01_70_110_gteOneBtaggedJet"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet", 240, 40, 160 );

   CreateUserTH1D( "MTenu_70_150_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_70_150_gteOneBtaggedJet"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_70_150_Njet_gte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_lte3_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_lte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_150_Njet_gte5_gteOneBtaggedJet", 240, 40, 160 );

   CreateUserTH1D( "MTenu_Type01_70_150_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_Type01_70_150_gteOneBtaggedJet"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_Type01_70_150_Njet_gte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_150_Njet_lte3_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_150_Njet_lte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_Type01_70_150_Njet_gte5_gteOneBtaggedJet", 240, 40, 160 );


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
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_50_110_Njet_lte3", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_50_110_Njet_lte3", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_50_110_Njet_lte3", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_50_110_Njet_lte3"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_50_110_Njet_lte3"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_50_110_Njet_lte3"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_50_110_Njet_lte4", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_50_110_Njet_lte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_50_110_Njet_lte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_50_110_Njet_lte4"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_50_110_Njet_lte4"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_50_110_Njet_lte4"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_50_110_Njet_gte4", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_50_110_Njet_gte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_50_110_Njet_gte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_50_110_Njet_gte4"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_50_110_Njet_gte4"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_50_110_Njet_gte4"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_50_110_Njet_gte5", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_50_110_Njet_gte5", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_50_110_Njet_gte5", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_50_110_Njet_gte5"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_50_110_Njet_gte5"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_50_110_Njet_gte5"     , 200 , 0 , 2000 );

   // 70-110
   CreateUserTH1D( "Pt1stEle_MTenu_70_110_Njet_lte3", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_110_Njet_lte3", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_110_Njet_lte3", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_110_Njet_lte3"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_110_Njet_lte3"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_110_Njet_lte3"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_110_Njet_lte4", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_110_Njet_lte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_110_Njet_lte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_110_Njet_lte4"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_110_Njet_lte4"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_110_Njet_lte4"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_110_Njet_gte4", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_110_Njet_gte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_110_Njet_gte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_110_Njet_gte4"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_110_Njet_gte4"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_110_Njet_gte4"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_110_Njet_gte5", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_110_Njet_gte5", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_110_Njet_gte5", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_110_Njet_gte5"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_110_Njet_gte5"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_110_Njet_gte5"     , 200 , 0 , 2000 ); 
   
   // 70-110 type01
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_110_Njet_lte3", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_110_Njet_lte3", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_110_Njet_lte3", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_110_Njet_lte3"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_110_Njet_lte3"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_110_Njet_lte3"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_110_Njet_lte4", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_110_Njet_lte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_110_Njet_lte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_110_Njet_lte4"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_110_Njet_lte4"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_110_Njet_lte4"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_110_Njet_gte4", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_110_Njet_gte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_110_Njet_gte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_110_Njet_gte4"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_110_Njet_gte4"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_110_Njet_gte4"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_110_Njet_gte5", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_110_Njet_gte5", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_110_Njet_gte5", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_110_Njet_gte5"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_110_Njet_gte5"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_110_Njet_gte5"     , 200 , 0 , 2000 );

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
   
   // 70-150 type01
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_150_Njet_lte3", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_150_Njet_lte3", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_150_Njet_lte3", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_150_Njet_lte3"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_150_Njet_lte3"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_150_Njet_lte3"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_150_Njet_lte4", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_150_Njet_lte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_150_Njet_lte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_150_Njet_lte4"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_150_Njet_lte4"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_150_Njet_lte4"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_150_Njet_gte4", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_150_Njet_gte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_150_Njet_gte4", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_150_Njet_gte4"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_150_Njet_gte4"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_150_Njet_gte4"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_150_Njet_gte5", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_150_Njet_gte5", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_150_Njet_gte5", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_150_Njet_gte5"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_150_Njet_gte5"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_150_Njet_gte5"     , 200 , 0 , 2000 );

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
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_50_110_Njet_lte3_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_50_110_Njet_lte3_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_50_110_Njet_lte3_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_50_110_Njet_lte3_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_50_110_Njet_lte3_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_50_110_Njet_lte3_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_50_110_Njet_lte4_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_50_110_Njet_lte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_50_110_Njet_lte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_50_110_Njet_lte4_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_50_110_Njet_lte4_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_50_110_Njet_lte4_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_50_110_Njet_gte4_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_50_110_Njet_gte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_50_110_Njet_gte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_50_110_Njet_gte4_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_50_110_Njet_gte4_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_50_110_Njet_gte4_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_50_110_Njet_gte5_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_50_110_Njet_gte5_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_50_110_Njet_gte5_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_50_110_Njet_gte5_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_50_110_Njet_gte5_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_50_110_Njet_gte5_noBtaggedJets"     , 200 , 0 , 2000 );

   // 70-110
   CreateUserTH1D( "Pt1stEle_MTenu_70_110_Njet_lte3_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_110_Njet_lte3_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_110_Njet_lte3_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_110_Njet_lte3_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_110_Njet_lte3_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_110_Njet_lte3_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_110_Njet_lte4_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_110_Njet_lte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_110_Njet_lte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_110_Njet_lte4_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_110_Njet_lte4_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_110_Njet_lte4_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_110_Njet_gte4_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_110_Njet_gte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_110_Njet_gte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_110_Njet_gte4_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_110_Njet_gte4_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_110_Njet_gte4_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_110_Njet_gte5_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_110_Njet_gte5_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_110_Njet_gte5_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_110_Njet_gte5_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_110_Njet_gte5_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_110_Njet_gte5_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   // 70-110 type01
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_110_Njet_lte3_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_110_Njet_lte3_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_110_Njet_lte3_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_110_Njet_lte3_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_110_Njet_lte3_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_110_Njet_lte3_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_110_Njet_lte4_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_110_Njet_lte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_110_Njet_lte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_110_Njet_lte4_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_110_Njet_lte4_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_110_Njet_lte4_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_110_Njet_gte4_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_110_Njet_gte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_110_Njet_gte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_110_Njet_gte4_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_110_Njet_gte4_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_110_Njet_gte4_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_110_Njet_gte5_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_110_Njet_gte5_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_110_Njet_gte5_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_110_Njet_gte5_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_110_Njet_gte5_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_110_Njet_gte5_noBtaggedJets"     , 200 , 0 , 2000 );

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
   
   // 70-150 type01
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_150_Njet_lte3_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_150_Njet_lte3_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_150_Njet_lte3_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_150_Njet_lte3_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_150_Njet_lte3_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_150_Njet_lte3_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_150_Njet_lte4_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_150_Njet_lte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_150_Njet_lte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_150_Njet_lte4_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_150_Njet_lte4_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_150_Njet_lte4_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_150_Njet_gte4_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_150_Njet_gte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_150_Njet_gte4_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_150_Njet_gte4_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_150_Njet_gte4_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_150_Njet_gte4_noBtaggedJets"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_150_Njet_gte5_noBtaggedJets", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_150_Njet_gte5_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_150_Njet_gte5_noBtaggedJets", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_150_Njet_gte5_noBtaggedJets"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_150_Njet_gte5_noBtaggedJets"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_150_Njet_gte5_noBtaggedJets"     , 200 , 0 , 2000 );

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
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet"     , 200 , 0 , 2000 );

   // 70-110
   CreateUserTH1D( "Pt1stEle_MTenu_70_110_Njet_lte3_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_110_Njet_lte3_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_110_Njet_lte3_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_110_Njet_lte3_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_110_Njet_lte3_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_110_Njet_lte3_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_110_Njet_lte4_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_110_Njet_lte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_110_Njet_lte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_110_Njet_lte4_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_110_Njet_lte4_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_110_Njet_lte4_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_110_Njet_gte4_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_110_Njet_gte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_110_Njet_gte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_110_Njet_gte4_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_110_Njet_gte4_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_110_Njet_gte4_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_70_110_Njet_gte5_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_70_110_Njet_gte5_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_70_110_Njet_gte5_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_70_110_Njet_gte5_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_70_110_Njet_gte5_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_70_110_Njet_gte5_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   // 70-110 type01
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet"     , 200 , 0 , 2000 );

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
   
   // 70-150 type01
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_150_Njet_lte3_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_150_Njet_lte3_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_150_Njet_lte3_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_150_Njet_lte3_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_150_Njet_lte3_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_150_Njet_lte3_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_150_Njet_lte4_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_150_Njet_lte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_150_Njet_lte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_150_Njet_lte4_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_150_Njet_lte4_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_150_Njet_lte4_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_150_Njet_gte4_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_150_Njet_gte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_150_Njet_gte4_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_150_Njet_gte4_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_150_Njet_gte4_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_150_Njet_gte4_gteOneBtaggedJet"     , 200 , 0 , 2000 ); 
   
   CreateUserTH1D( "Pt1stEle_MTenu_Type01_70_150_Njet_gte5_gteOneBtaggedJet", 100 , 0 , 1000 ); 
   CreateUserTH1D( "Pt1stJet_MTenu_Type01_70_150_Njet_gte5_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "Pt2ndJet_MTenu_Type01_70_150_Njet_gte5_gteOneBtaggedJet", 200 , 0 , 2000 ); 
   CreateUserTH1D( "MET_MTenu_Type01_70_150_Njet_gte5_gteOneBtaggedJet"     , 200 , 0 , 1000 ); 
   CreateUserTH1D( "sT_MTenu_Type01_70_150_Njet_gte5_gteOneBtaggedJet"      , 300 , 0 , 3000 ); 
   CreateUserTH1D( "Mej_MTenu_Type01_70_150_Njet_gte5_gteOneBtaggedJet"     , 200 , 0 , 2000 );

   char plot_name[200];
   
   for (int i_lq_mass = 0; i_lq_mass < n_lq_mass ; ++i_lq_mass ) { 
     int lq_mass = LQ_MASS[i_lq_mass];
     sprintf(plot_name, "likelihood_LQ%d"   , lq_mass ); CreateUserTH1D ( plot_name, 200, 0.0, 1.0 );
     sprintf(plot_name, "logLikelihood_LQ%d", lq_mass ); CreateUserTH1D ( plot_name, 100, 0.0, 100.0 );
   }
   
   //--------------------------------------------------------------------------
   // Final selection plots
   //--------------------------------------------------------------------------
   
   if ( do_final_selection ) { 
     for (int i_lq_mass = 0; i_lq_mass < n_lq_mass ; ++i_lq_mass ) { 
       int lq_mass = LQ_MASS[i_lq_mass];
       sprintf(plot_name, "split_1fb_LQ%d", lq_mass ); CreateUserTH1D ( plot_name, 21  , -0.5, 20.5);
       sprintf(plot_name, "MET_LQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 
       //sprintf(plot_name, "MET_UnclUp_LQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 
       //sprintf(plot_name, "MET_UnclDown_LQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 
       sprintf(plot_name, "Mej_LQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "Mejj_LQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "Mjj_LQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "Mej3_LQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "ST_LQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "MTenu_LQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
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
       sprintf(plot_name, "MTenu_StLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MET_StLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_MtLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_MtLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_MtLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MET_MtLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_MejLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_MejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_MejLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MET_MejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_MetLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_MetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_MetLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MET_MetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_stAndMtAndMetLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_stAndMtAndMetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_stAndMtAndMetLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MET_stAndMtAndMetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_stAndMtAndMejLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_stAndMtAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_stAndMtAndMejLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MET_stAndMtAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_stAndMetAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_stAndMetAndMejLQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_stAndMetAndMejLQ%d"   , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MET_stAndMetAndMejLQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_mtAndMetAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_mtAndMetAndMejLQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_mtAndMetAndMejLQ%d"   , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MET_mtAndMetAndMejLQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 
       
       sprintf(plot_name, "ST_StAndMejLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_StAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_StAndMejLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MET_StAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_MtAndMejLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_MtAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_MtAndMejLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MET_MtAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_MetAndMejLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_MetAndMejLQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_MetAndMejLQ%d"   , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MET_MetAndMejLQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_StAndMetLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_StAndMetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_StAndMetLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MET_StAndMetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_MtAndMetLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_MtAndMetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_MtAndMetLQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MET_MtAndMetLQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 

       sprintf(plot_name, "ST_StAndMtLQ%d"        , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_StAndMtLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "MTenu_StAndMtLQ%d"     , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 );
       sprintf(plot_name, "MET_StAndMtLQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 1000 ); 
       
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
	 sprintf(plot_name, "Eta1stEle_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
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
       }
       
       if ( do_eleQuality_plots ) { 
	 sprintf(plot_name, "BeamSpotDXY_1stEle_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   0.5  );
	 sprintf(plot_name, "Classif_1stEle_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name , 5  , -0.5 ,   4.5  );
	 sprintf(plot_name, "CorrIsolation_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200,-25.0 ,  25.0  );
	 sprintf(plot_name, "DeltaEtaTrkSC_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.01,   0.01 );
	 sprintf(plot_name, "DeltaPhiTrkSC_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.1 ,   0.1  );
	 sprintf(plot_name, "E1x5OverE5x5_1stEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );
	 sprintf(plot_name, "E2x5OverE5x5_1stEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );
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
	 sprintf(plot_name, "SigmaIEtaIEta_Barrel_1stEle_LQ%d" , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.02 );
	 sprintf(plot_name, "SigmaIEtaIEta_Endcap_1stEle_LQ%d" , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.1  );
	 sprintf(plot_name, "TrkPtOPt_1stEle_LQ%d"             , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,  100.0  );
	 sprintf(plot_name, "ValidFrac_1stEle_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );
	 sprintf(plot_name, "EOverP_1stEle_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name , 20000,  0.0 , 100.0  );
       }
     }
   }
   
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
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
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

     //if ( isData && Ele2_ValidFrac > 998. && Ele1_ValidFrac > 998. ){
     //  double mu_trigger_effi_1 = 0.0;
     //  if      ( fabs(Muon1_Eta) > 0.0 && fabs(Muon1_Eta) <= 0.9 ) mu_trigger_effi_1 = 0.94;
     //  else if ( fabs(Muon1_Eta) > 0.9 && fabs(Muon1_Eta) <= 1.2 ) mu_trigger_effi_1 = 0.84;
     //  else if ( fabs(Muon1_Eta) > 1.2 && fabs(Muon1_Eta) <= 2.1 ) mu_trigger_effi_1 = 0.82;

     //  double ele_trigger_effi_1 = 0.0;
     //  if      ( fabs(Muon1_Eta) > 0.000 && fabs(Muon1_Eta) <= 1.442 ) ele_trigger_effi_1 = 0.974;
     //  else if ( fabs(Muon1_Eta) > 1.442 && fabs(Muon1_Eta) <= 2.500 ) ele_trigger_effi_1 = 0.958;

     //  double mt_resolution_weight = mt_weight_func_muon_background -> Eval ( MT_Ele1MET );
     //  // double mt_resolution_weight = 1.0;

     //  double enujj_over_munujj_reco_weight = 0.972;
     //  double total_weight = mt_resolution_weight * enujj_over_munujj_reco_weight * ele_trigger_effi_1 / mu_trigger_effi_1;

     //  gen_weight = total_weight;
     //  
     //}

     //--------------------------------------------------------------------------
     // Reweighting MET & MT
     //--------------------------------------------------------------------------

     /*
     if ( isData == 0 && Ele2_ValidFrac < 998. && Ele1_ValidFrac < 998. ){
       double met_weight      = met_weight_func      -> Eval ( PFMET_Type01XY_Pt );
       double mt_weight       = mt_weight_func       -> Eval ( MT_Ele1MET        );
       double mt_weight_first = mt_weight_func_first -> Eval ( MT_Ele1MET        );
       // gen_weight = met_weight * mt_weight;
       gen_weight = mt_weight_first;
     }
     */
	  
     //if ( Ele2_ValidFrac > 998. && Ele1_ValidFrac > 998. ) nMuon_ptCut = 0;

     //--------------------------------------------------------------------------
     // Reweight W+Jets binned samples
     //--------------------------------------------------------------------------

     /*
     std::string current_file_name ( fChain->GetCurrentFile()->GetName() );
     
     if      ( ProcessID == 0 && current_file_name.find ( std::string ("W0Jet")) != std::string::npos ) gen_weight *= 1.0; // Don't touch, not enough stats
     else if ( ProcessID == 0 && current_file_name.find ( std::string ("W1Jet")) != std::string::npos ) gen_weight *= 1.0; // Don't touch, not enough stats
     else if ( ProcessID == 0 && current_file_name.find ( std::string ("W2Jet")) != std::string::npos ) gen_weight *= 0.995806;
     else if ( ProcessID == 0 && current_file_name.find ( std::string ("W3Jet")) != std::string::npos ) gen_weight *= 1.023070;
     else if ( ProcessID == 0 && current_file_name.find ( std::string ("W4Jet")) != std::string::npos ) gen_weight *= 0.966947;
     else if ( ProcessID == 1 && current_file_name.find ( std::string ("W0Jet")) != std::string::npos ) gen_weight *= 1.0; // Don't touch, not enough stats
     else if ( ProcessID == 1 && current_file_name.find ( std::string ("W1Jet")) != std::string::npos ) gen_weight *= 1.0; // Don't touch, not enough stats
     else if ( ProcessID == 1 && current_file_name.find ( std::string ("W2Jet")) != std::string::npos ) gen_weight *= 0.869381;
     else if ( ProcessID == 1 && current_file_name.find ( std::string ("W3Jet")) != std::string::npos ) gen_weight *= 0.972974;
     else if ( ProcessID == 1 && current_file_name.find ( std::string ("W4Jet")) != std::string::npos ) gen_weight *= 0.954368;
     */
     
     //--------------------------------------------------------------------------
     // Is this a barrel electron?
     //--------------------------------------------------------------------------

     bool Ele1_IsBarrel = bool ( fabs (Ele1_Eta) < 1.442 );

     //--------------------------------------------------------------------------
     // First variable to fill just shows the "reweighting".  Always passes.
     //--------------------------------------------------------------------------

     fillVariableWithValue ( "Reweighting", 1, gen_weight * pileup_weight  );

     //--------------------------------------------------------------------------
     // Special treatment of inclusive W/Z
     //--------------------------------------------------------------------------
     bool passGenWZPt = true;
     std::string current_file_name ( fChain->GetCurrentFile()->GetName());
     // inclusive
     if(current_file_name.find("WJetsToLNu_ext1_amcatnloFXFX") != std::string::npos 
         || current_file_name.find("WJetsToLNu_amcatnloFXFX") != std::string::npos) {
       if(GenW1_Pt > 120) passGenWZPt = false; // if W Pt > 120 GeV, cut it out
     }
     if(current_file_name.find("DYJetsToLL_M-50_amcatnloFXFX") != std::string::npos) {
       if(GenZGamma1_Pt > 120) passGenWZPt = false; // if Z/gamma Pt > 120 GeV, cut it out
     }
     // first pt bin
     if(current_file_name.find("WJetsToLNu_Pt-100") != std::string::npos) {
       if(GenW1_Pt <= 120) passGenWZPt = false;
     }
     if(current_file_name.find("DYJetsToLL_Pt-100") != std::string::npos) {
       if(GenZGamma1_Pt <= 120) passGenWZPt = false;
     }
     //// testing
     //if(current_file_name.find("WJetsToLNu_ext1_amcatnloFXFX") != std::string::npos 
     //    || current_file_name.find("WJetsToLNu_amcatnloFXFX") != std::string::npos) {
     //  if(GenW1_Pt <= 100) passGenWZPt = false; // if W Pt > 100 GeV, cut it out
     //}
     //if(current_file_name.find("DYJetsToLL_M-50_amcatnloFXFX") != std::string::npos) {
     //  if(GenZGamma1_Pt <= 100) passGenWZPt = false; // if Z/gamma Pt > 100 GeV, cut it out
     //}
     fillVariableWithValue("PassGenWZPt",passGenWZPt,gen_weight*pileup_weight);

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

     int passHLT = 1;
     //if ( isData ) { 
     //  passHLT = 0;
     //  if ( H_Ele27_WPTight == 1 || H_Photon175 == 1)
     //    passHLT = 1;
     //}
     //// FIXME: Update to new 2016 curve?
     //else // using the turn-on in the MC
     //{
     //  // a la Z', throw a random number and if it's below the efficiency at this pt/eta, pass the event
     //  //   we get two chances to pass since we may have two electrons in the event
     //  // UPDATE TO SCETA/ptheep
     //  passHLT = trigEle27::passTrig(Ele1_PtHeep,Ele1_SCEta) ? 1 : 0;
     //  if(!passHLT) // if the first one doesn't pass, try the second one
     //    passHLT = trigEle27::passTrig(Ele2_PtHeep,Ele2_SCEta) ? 1 : 0;
     //}
     if (isData) {
       passHLT = 0;
       if (H_Ele27_WPTight_eta2p1 == 1)
         passHLT = 1;
       //if ( H_Ele27_WPLoose == 1)
       //if ( H_Ele27_WPTight == 1)
       //if ( H_Ele27_WPTight == 1 || H_Photon175 == 1)
       //if ( H_Ele27_WPLoose_eta2p1 == 1)
       //  passHLT = 1;
       //if(run < 273726) // bad endcap alignment affecting deltaEtaIn cut
       //  passHLT = 0;
       //if(run < 275676) // L1 EGM efficiency going to zero
       //  passHLT = 0;
     }
     else {
       passHLT = trigEle27::passTrig(Ele1_PtHeep,Ele1_SCEta) ? 1 : 0;
       // for enujj, we only get one chance to trigger
     }

     fillVariableWithValue ( "PassHLT", passHLT, gen_weight * pileup_weight  ) ;     

     //--------------------------------------------------------------------------
     // Fill more variables
     //--------------------------------------------------------------------------
     // Muon variables ( for veto ) 					      
     fillVariableWithValue(   "nMuon"                    , nMuon_ptCut             , gen_weight * pileup_weight );
			                                      		                
     // 1st Electron variables
     fillVariableWithValue(   "nEle"                     , nEle_ptCut              , gen_weight * pileup_weight ); 
     fillVariableWithValue(   "Ele1_PtHeep"              , Ele1_PtHeep             , gen_weight * pileup_weight );
     fillVariableWithValue(   "Ele1_Eta"                 , Ele1_Eta                , gen_weight * pileup_weight );
     fillVariableWithValue(   "Ele1_IsBarrel"            , Ele1_IsBarrel           , gen_weight * pileup_weight );
     fillVariableWithValue(   "MTenu"                    , MT_Ele1MET              , gen_weight * pileup_weight );

     // MET variables	                                      		           
     fillVariableWithValue(   "MET"                      , PFMET_Type1XY_Pt       , gen_weight * pileup_weight );
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
     //double MT_Jet1MET_Type01XY_UnclEnUp, MT_Jet2MET_Type01XY_UnclEnUp, MT_Ele1MET_Type01XY_UnclEnUp;
     //double MT_Jet1MET_Type01XY_UnclEnDown, MT_Jet2MET_Type01XY_UnclEnDown, MT_Ele1MET_Type01XY_UnclEnDown;
     double mDPhi_METType01XYEleEnUp_Ele1, mDPhi_METType01XYEleEnDown_Ele1;
     double mDPhi_METType01XYEleEnUp_Jet1, mDPhi_METType01XYEleEnDown_Jet1, mDPhi_METType01XYEleEnUp_Jet2, mDPhi_METType01XYEleEnDown_Jet2;
     double mDPhi_METType01XYJetEnUp_Ele1, mDPhi_METType01XYJetEnDown_Ele1;
     double mDPhi_METType01XYJetEnUp_Jet1, mDPhi_METType01XYJetEnDown_Jet1, mDPhi_METType01XYJetEnUp_Jet2, mDPhi_METType01XYJetEnDown_Jet2;
     double mDPhi_METType01XYJetResUp_Ele1, mDPhi_METType01XYJetResDown_Ele1;
     double mDPhi_METType01XYJetResUp_Jet1, mDPhi_METType01XYJetResDown_Jet1, mDPhi_METType01XYJetResUp_Jet2, mDPhi_METType01XYJetResDown_Jet2;
     //double mDPhi_METType01XYUnclEnUp_Ele1, mDPhi_METType01XYUnclEnDown_Ele1;
     //double mDPhi_METType01XYUnclEnUp_Jet1, mDPhi_METType01XYUnclEnDown_Jet1, mDPhi_METType01XYUnclEnUp_Jet2, mDPhi_METType01XYUnclEnDown_Jet2;
     //// XXX I shall redefine the ST here
     //if ( nJet_store >= 2 ) 
     //  sT_enujj = Jet1_Pt+Jet2_Pt+Ele1_Pt+PFMET_Type01XY_Pt;
     //  //sT_enujj = Jet1_Pt+Jet2_Pt+Ele1_Pt+PFMET_Type01XY_UnclUp_Pt;
     //  //sT_enujj = Jet1_Pt+Jet2_Pt+Ele1_Pt+PFMET_Type01XY_UnclDown_Pt;
     //// SIC end add

     // alternate METs

     if ( nEle_store > 0 ) {
       TVector2 v_ele;
       TVector2 v_MET_Type01;
       TLorentzVector v_ele_lorentz, v_met_lorentz, v_ele_met_lorentz;
       v_ele_lorentz.SetPtEtaPhiM( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
       v_met_lorentz.SetPtEtaPhiM( PFMET_Type1XY_Pt, 0.0, PFMET_Type1XY_Phi, 0.0 );
       v_ele_met_lorentz = v_ele_lorentz + v_met_lorentz;
       

       v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       v_ele.SetMagPhi( Ele1_Pt, Ele1_Phi );
       mDPhi_METType01_Ele1 = fabs(v_MET_Type01.DeltaPhi ( v_ele ));
       float deltaphi = v_MET_Type01.DeltaPhi(v_ele);
       MT_Ele1MET_Type01 = sqrt ( 2 * Ele1_Pt * PFMET_Type01_Pt * ( 1 - cos ( deltaphi ) ) );
       //// SIC add
       //// Here I think we can calculate the MT with the shifted METs too and use them later in the final selection
       //// define them above with the others
       //// Ele EnUp/Down
       //TVector2 v_MET_Type01XYEleEnUp;
       //v_MET_Type01XYEleEnUp.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       //v_ele.SetMagPhi( Ele1_Pt, Ele1_Phi );
       //mDPhi_METType01XYEleEnUp_Ele1 = fabs(v_MET_Type01XYEleEnUp.DeltaPhi ( v_ele ));
       //deltaphi = v_MET_Type01XYEleEnUp.DeltaPhi(v_ele);
       //MT_Ele1MET_Type01XY_EleEnUp = sqrt ( 2 * Ele1_Pt * PFMET_Type01XY_EleEnUp_Pt * ( 1 - cos ( deltaphi ) ) );
       //TVector2 v_MET_Type01XYEleEnDown;
       //v_MET_Type01XYEleEnDown.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       //v_ele.SetMagPhi( Ele1_Pt, Ele1_Phi );
       //mDPhi_METType01XYEleEnDown_Ele1 = fabs(v_MET_Type01XYEleEnDown.DeltaPhi ( v_ele ));
       //deltaphi = v_MET_Type01XYEleEnDown.DeltaPhi(v_ele);
       //MT_Ele1MET_Type01XY_EleEnDown = sqrt ( 2 * Ele1_Pt * PFMET_Type01XY_EleEnDown_Pt * ( 1 - cos ( deltaphi ) ) );
       //// Jet EnUp/Down
       //TVector2 v_MET_Type01XYJetEnUp;
       //v_MET_Type01XYJetEnUp.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       //v_ele.SetMagPhi( Jet1_Pt, Jet1_Phi );
       //mDPhi_METType01XYJetEnUp_Ele1 = fabs(v_MET_Type01XYJetEnUp.DeltaPhi ( v_ele ));
       //deltaphi = v_MET_Type01XYJetEnUp.DeltaPhi(v_ele);
       //MT_Jet1MET_Type01XY_JetEnUp = sqrt ( 2 * Jet1_Pt * PFMET_Type01XY_JetEnUp_Pt * ( 1 - cos ( deltaphi ) ) );
       //TVector2 v_MET_Type01XYJetEnDown;
       //v_MET_Type01XYJetEnDown.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       //v_ele.SetMagPhi( Jet1_Pt, Jet1_Phi );
       //mDPhi_METType01XYJetEnDown_Ele1 = fabs(v_MET_Type01XYJetEnDown.DeltaPhi ( v_ele ));
       //deltaphi = v_MET_Type01XYJetEnDown.DeltaPhi(v_ele);
       //MT_Jet1MET_Type01XY_JetEnDown = sqrt ( 2 * Jet1_Pt * PFMET_Type01XY_JetEnDown_Pt * ( 1 - cos ( deltaphi ) ) );
       //// Jet ResUp/Down
       //TVector2 v_MET_Type01XYJetResUp;
       //v_MET_Type01XYJetResUp.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       //v_ele.SetMagPhi( Ele1_Pt, Ele1_Phi );
       //mDPhi_METType01XYJetResUp_Ele1 = fabs(v_MET_Type01XYJetResUp.DeltaPhi ( v_ele ));
       //deltaphi = v_MET_Type01XYJetResUp.DeltaPhi(v_ele);
       //MT_Ele1MET_Type01XY_JetResUp = sqrt ( 2 * Ele1_Pt * PFMET_Type01XY_JetResUp_Pt * ( 1 - cos ( deltaphi ) ) );
       //TVector2 v_MET_Type01XYJetResDown;
       //v_MET_Type01XYJetResDown.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       //v_ele.SetMagPhi( Ele1_Pt, Ele1_Phi );
       //mDPhi_METType01XYJetResDown_Ele1 = fabs(v_MET_Type01XYJetResDown.DeltaPhi ( v_ele ));
       //deltaphi = v_MET_Type01XYJetResDown.DeltaPhi(v_ele);
       //MT_Ele1MET_Type01XY_JetResDown = sqrt ( 2 * Ele1_Pt * PFMET_Type01XY_JetResDown_Pt * ( 1 - cos ( deltaphi ) ) );
       //// UnclEnUp/Down
       //TVector2 v_MET_Type01XYUnclEnUp;
       //v_MET_Type01XYUnclEnUp.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       //v_ele.SetMagPhi( Ele1_Pt, Ele1_Phi );
       //mDPhi_METType01XYUnclEnUp_Ele1 = fabs(v_MET_Type01XYUnclEnUp.DeltaPhi ( v_ele ));
       //deltaphi = v_MET_Type01XYUnclEnUp.DeltaPhi(v_ele);
       //MT_Ele1MET_Type01XY_UnclEnUp = sqrt ( 2 * Ele1_Pt * PFMET_Type01XY_UnclUp_Pt * ( 1 - cos ( deltaphi ) ) );
       //TVector2 v_MET_Type01XYUnclEnDown;
       //v_MET_Type01XYUnclEnDown.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       //v_ele.SetMagPhi( Ele1_Pt, Ele1_Phi );
       //mDPhi_METType01XYUnclEnDown_Ele1 = fabs(v_MET_Type01XYUnclEnDown.DeltaPhi ( v_ele ));
       //deltaphi = v_MET_Type01XYUnclEnDown.DeltaPhi(v_ele);
       //MT_Ele1MET_Type01XY_UnclEnDown = sqrt ( 2 * Ele1_Pt * PFMET_Type01XY_UnclDown_Pt * ( 1 - cos ( deltaphi ) ) );
       //// SIC end add
     }
     
     // 1st JET variables                                     		           
     if ( nJet_store > 0 ) { 						           
       fillVariableWithValue( "Jet1_Pt"                  , Jet1_Pt                 , gen_weight * pileup_weight );
       fillVariableWithValue( "Jet1_Eta"                 , Jet1_Eta                , gen_weight * pileup_weight );
       fillVariableWithValue( "mDeltaPhiMET1stJet"       , mDPhi_METJet1           , gen_weight * pileup_weight );

       TVector2 v_MET;
       TVector2 v_jet;
       TVector2 v_MET_Type01;
       v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       v_MET.SetMagPhi( PFMET_Type1XY_Pt , PFMET_Type1XY_Phi  );
       v_jet.SetMagPhi( Jet1_Pt, Jet1_Phi );
       mDPhi_METType01_Jet1 = fabs(v_MET_Type01.DeltaPhi ( v_jet ));
       float deltaphi = v_MET.DeltaPhi(v_jet);
       MT_Jet1MET = sqrt ( 2 * Jet1_Pt * PFMET_Type1XY_Pt * ( 1 - cos ( deltaphi ) ) );
       //// SIC add
       //TVector2 v_MET_Type01XYEleEnUp, v_MET_Type01XYEleEnDown, v_MET_Type01XYJetEnUp, v_MET_Type01XYJetEnDown;
       //TVector2 v_MET_Type01XYJetResUp, v_MET_Type01XYJetResDown, v_MET_Type01XYUnclEnUp, v_MET_Type01XYUnclEnDown;
       //// EleEnUp/Down
       //v_MET_Type01XYEleEnUp.SetMagPhi( PFMET_Type01XY_EleEnUp_Pt , PFMET_Type01XY_EleEnUp_Phi  );
       //mDPhi_METType01XYEleEnUp_Jet1 = fabs(v_MET_Type01XYEleEnUp.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYEleEnUp.DeltaPhi(v_jet);
       //MT_Jet1MET_Type01XY_EleEnUp = sqrt ( 2 * Jet1_Pt * PFMET_Type01XY_EleEnUp_Pt * ( 1 - cos ( deltaphi ) ) );
       //v_MET_Type01XYEleEnDown.SetMagPhi( PFMET_Type01XY_EleEnDown_Pt , PFMET_Type01XY_EleEnDown_Phi  );
       //mDPhi_METType01XYEleEnDown_Jet1 = fabs(v_MET_Type01XYEleEnDown.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYEleEnDown.DeltaPhi(v_jet);
       //MT_Jet1MET_Type01XY_EleEnDown = sqrt ( 2 * Jet1_Pt * PFMET_Type01XY_EleEnDown_Pt * ( 1 - cos ( deltaphi ) ) );
       //// JetEnUp/Down
       //v_MET_Type01XYJetEnUp.SetMagPhi( PFMET_Type01XY_JetEnUp_Pt , PFMET_Type01XY_JetEnUp_Phi  );
       //mDPhi_METType01XYJetEnUp_Jet1 = fabs(v_MET_Type01XYJetEnUp.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYJetEnUp.DeltaPhi(v_jet);
       //MT_Jet1MET_Type01XY_JetEnUp = sqrt ( 2 * Jet1_Pt * PFMET_Type01XY_JetEnUp_Pt * ( 1 - cos ( deltaphi ) ) );
       //v_MET_Type01XYJetEnDown.SetMagPhi( PFMET_Type01XY_JetEnDown_Pt , PFMET_Type01XY_JetEnDown_Phi  );
       //mDPhi_METType01XYJetEnDown_Jet1 = fabs(v_MET_Type01XYJetEnDown.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYJetEnDown.DeltaPhi(v_jet);
       //MT_Jet1MET_Type01XY_JetEnDown = sqrt ( 2 * Jet1_Pt * PFMET_Type01XY_JetEnDown_Pt * ( 1 - cos ( deltaphi ) ) );
       //// JetResUp/Down
       //v_MET_Type01XYJetResUp.SetMagPhi( PFMET_Type01XY_JetResUp_Pt , PFMET_Type01XY_JetResUp_Phi  );
       //mDPhi_METType01XYJetResUp_Jet1 = fabs(v_MET_Type01XYJetResUp.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYJetResUp.DeltaPhi(v_jet);
       //MT_Jet1MET_Type01XY_JetResUp = sqrt ( 2 * Jet1_Pt * PFMET_Type01XY_JetResUp_Pt * ( 1 - cos ( deltaphi ) ) );
       //v_MET_Type01XYJetResDown.SetMagPhi( PFMET_Type01XY_JetResDown_Pt , PFMET_Type01XY_JetResDown_Phi  );
       //mDPhi_METType01XYJetResDown_Jet1 = fabs(v_MET_Type01XYJetResDown.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYJetResDown.DeltaPhi(v_jet);
       //MT_Jet1MET_Type01XY_JetResDown = sqrt ( 2 * Jet1_Pt * PFMET_Type01XY_JetResDown_Pt * ( 1 - cos ( deltaphi ) ) );
       //// UnclEnUp/Down
       //v_MET_Type01XYUnclEnUp.SetMagPhi( PFMET_Type01XY_UnclUp_Pt , PFMET_Type01XY_UnclUp_Phi  );
       //mDPhi_METType01XYUnclEnUp_Jet1 = fabs(v_MET_Type01XYUnclEnUp.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYUnclEnUp.DeltaPhi(v_jet);
       //MT_Jet1MET_Type01XY_UnclEnUp = sqrt ( 2 * Jet1_Pt * PFMET_Type01XY_UnclUp_Pt * ( 1 - cos ( deltaphi ) ) );
       //v_MET_Type01XYUnclEnDown.SetMagPhi( PFMET_Type01XY_UnclDown_Pt , PFMET_Type01XY_UnclDown_Phi  );
       //mDPhi_METType01XYUnclEnDown_Jet1 = fabs(v_MET_Type01XYUnclEnDown.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYUnclEnDown.DeltaPhi(v_jet);
       //MT_Jet1MET_Type01XY_UnclEnDown = sqrt ( 2 * Jet1_Pt * PFMET_Type01XY_UnclDown_Pt * ( 1 - cos ( deltaphi ) ) );
       //// SIC end add
     }									           
     
     // 2nd JET variables                                     		           
     if ( nJet_store > 1 ) { 	                                      	           
       fillVariableWithValue( "Jet2_Pt"                  , Jet2_Pt                 , gen_weight * pileup_weight );
       fillVariableWithValue( "Jet2_Eta"                 , Jet2_Eta                , gen_weight * pileup_weight );
       fillVariableWithValue( "ST"                       , sT_enujj                , gen_weight * pileup_weight );
       
       TVector2 v_MET;
       TVector2 v_jet;
       TVector2 v_MET_Type01;
       v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       v_MET.SetMagPhi( PFMET_Type1XY_Pt , PFMET_Type1XY_Phi );
       v_jet.SetMagPhi( Jet2_Pt, Jet2_Phi );
       mDPhi_METType01_Jet2 = fabs(v_MET_Type01.DeltaPhi ( v_jet ));
       float deltaphi = v_MET.DeltaPhi(v_jet);
       MT_Jet2MET = sqrt ( 2 * Jet2_Pt * PFMET_Type1XY_Pt * ( 1 - cos ( deltaphi ) ) );
       //// SIC add
       //TVector2 v_MET_Type01XYEleEnUp, v_MET_Type01XYEleEnDown, v_MET_Type01XYJetEnUp, v_MET_Type01XYJetEnDown;
       //TVector2 v_MET_Type01XYJetResUp, v_MET_Type01XYJetResDown, v_MET_Type01XYUnclEnUp, v_MET_Type01XYUnclEnDown;
       //// EleEnUp/Down
       //v_MET_Type01XYEleEnUp.SetMagPhi( PFMET_Type01XY_EleEnUp_Pt , PFMET_Type01XY_EleEnUp_Phi  );
       //mDPhi_METType01XYEleEnUp_Jet2 = fabs(v_MET_Type01XYEleEnUp.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYEleEnUp.DeltaPhi(v_jet);
       //MT_Jet2MET_Type01XY_EleEnUp = sqrt ( 2 * Jet2_Pt * PFMET_Type01XY_EleEnUp_Pt * ( 1 - cos ( deltaphi ) ) );
       //v_MET_Type01XYEleEnDown.SetMagPhi( PFMET_Type01XY_EleEnDown_Pt , PFMET_Type01XY_EleEnDown_Phi  );
       //mDPhi_METType01XYEleEnDown_Jet2 = fabs(v_MET_Type01XYEleEnDown.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYEleEnDown.DeltaPhi(v_jet);
       //MT_Jet2MET_Type01XY_EleEnDown = sqrt ( 2 * Jet2_Pt * PFMET_Type01XY_EleEnDown_Pt * ( 1 - cos ( deltaphi ) ) );
       //// JetEnUp/Down
       //v_MET_Type01XYJetEnUp.SetMagPhi( PFMET_Type01XY_JetEnUp_Pt , PFMET_Type01XY_JetEnUp_Phi  );
       //mDPhi_METType01XYJetEnUp_Jet2 = fabs(v_MET_Type01XYJetEnUp.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYJetEnUp.DeltaPhi(v_jet);
       //MT_Jet2MET_Type01XY_JetEnUp = sqrt ( 2 * Jet2_Pt * PFMET_Type01XY_JetEnUp_Pt * ( 1 - cos ( deltaphi ) ) );
       //v_MET_Type01XYJetEnDown.SetMagPhi( PFMET_Type01XY_JetEnDown_Pt , PFMET_Type01XY_JetEnDown_Phi  );
       //mDPhi_METType01XYJetEnDown_Jet2 = fabs(v_MET_Type01XYJetEnDown.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYJetEnDown.DeltaPhi(v_jet);
       //MT_Jet2MET_Type01XY_JetEnDown = sqrt ( 2 * Jet2_Pt * PFMET_Type01XY_JetEnDown_Pt * ( 1 - cos ( deltaphi ) ) );
       //// JetResUp/Down
       //v_MET_Type01XYJetResUp.SetMagPhi( PFMET_Type01XY_JetResUp_Pt , PFMET_Type01XY_JetResUp_Phi  );
       //mDPhi_METType01XYJetResUp_Jet2 = fabs(v_MET_Type01XYJetResUp.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYJetResUp.DeltaPhi(v_jet);
       //MT_Jet2MET_Type01XY_JetResUp = sqrt ( 2 * Jet2_Pt * PFMET_Type01XY_JetResUp_Pt * ( 1 - cos ( deltaphi ) ) );
       //v_MET_Type01XYJetResDown.SetMagPhi( PFMET_Type01XY_JetResDown_Pt , PFMET_Type01XY_JetResDown_Phi  );
       //mDPhi_METType01XYJetResDown_Jet2 = fabs(v_MET_Type01XYJetResDown.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYJetResDown.DeltaPhi(v_jet);
       //MT_Jet2MET_Type01XY_JetResDown = sqrt ( 2 * Jet2_Pt * PFMET_Type01XY_JetResDown_Pt * ( 1 - cos ( deltaphi ) ) );
       //// UnclEnUp/Down
       //v_MET_Type01XYUnclEnUp.SetMagPhi( PFMET_Type01XY_UnclUp_Pt , PFMET_Type01XY_UnclUp_Phi  );
       //mDPhi_METType01XYUnclEnUp_Jet2 = fabs(v_MET_Type01XYUnclEnUp.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYUnclEnUp.DeltaPhi(v_jet);
       //MT_Jet2MET_Type01XY_UnclEnUp = sqrt ( 2 * Jet2_Pt * PFMET_Type01XY_UnclUp_Pt * ( 1 - cos ( deltaphi ) ) );
       //v_MET_Type01XYUnclEnDown.SetMagPhi( PFMET_Type01XY_UnclDown_Pt , PFMET_Type01XY_UnclDown_Phi  );
       //mDPhi_METType01XYUnclEnDown_Jet2 = fabs(v_MET_Type01XYUnclEnDown.DeltaPhi ( v_jet ));
       //deltaphi = v_MET_Type01XYUnclEnDown.DeltaPhi(v_jet);
       //MT_Jet2MET_Type01XY_UnclEnDown = sqrt ( 2 * Jet2_Pt * PFMET_Type01XY_UnclDown_Pt * ( 1 - cos ( deltaphi ) ) );
       //// SIC end add
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


     double MT_JetMET;
     double MT_EleJet;
     double Mej;
     
     if ( fabs ( MT_Jet1MET - MT_Ele1Jet2 ) < fabs( MT_Jet2MET - MT_Ele1Jet1 )){
       MT_JetMET = MT_Jet1MET;
       MT_EleJet = MT_Ele1Jet2;
       Mej = M_e1j2;
     } else { 
       MT_JetMET = MT_Jet2MET;
       MT_EleJet = MT_Ele1Jet1;
       Mej = M_e1j1;
     }	 
     
     // Optimization variables
     fillVariableWithValue( "ST_opt"   , sT_enujj   , gen_weight * pileup_weight );
     fillVariableWithValue( "Mej_opt"  , Mej        , gen_weight * pileup_weight );
     fillVariableWithValue( "MET_opt"  , PFMET_Type1XY_Pt , gen_weight * pileup_weight );
     fillVariableWithValue( "MT_opt"   , MT_Ele1MET , gen_weight * pileup_weight );
     
     // Dummy variables
     fillVariableWithValue ("preselection",1, gen_weight * pileup_weight ); 
	
     //--------------------------------------------------------------------------
     // Fill bjet variables
     //--------------------------------------------------------------------------
     
     double btagCSV_loose_cut  = 0.460;
     double btagCSV_medium_cut = 0.800;
     double btagCSV_tight_cut  = 0.935;
       
     int nBJet_loose_ptCut  = 0;
     int nBJet_medium_ptCut = 0;
     int nBJet_tight_ptCut  = 0;
     
     if ( Jet1_btagCSV > btagCSV_loose_cut  ) nBJet_loose_ptCut++;
     if ( Jet2_btagCSV > btagCSV_loose_cut  ) nBJet_loose_ptCut++;
     if ( Jet3_btagCSV > btagCSV_loose_cut  ) nBJet_loose_ptCut++;
     if ( Jet4_btagCSV > btagCSV_loose_cut  ) nBJet_loose_ptCut++;
     if ( Jet5_btagCSV > btagCSV_loose_cut  ) nBJet_loose_ptCut++;
     
     if ( Jet1_btagCSV > btagCSV_medium_cut ) nBJet_medium_ptCut++;
     if ( Jet2_btagCSV > btagCSV_medium_cut ) nBJet_medium_ptCut++;
     if ( Jet3_btagCSV > btagCSV_medium_cut ) nBJet_medium_ptCut++;
     if ( Jet4_btagCSV > btagCSV_medium_cut ) nBJet_medium_ptCut++;
     if ( Jet5_btagCSV > btagCSV_medium_cut ) nBJet_medium_ptCut++;
     
     if ( Jet1_btagCSV > btagCSV_tight_cut  ) nBJet_tight_ptCut++;
     if ( Jet2_btagCSV > btagCSV_tight_cut  ) nBJet_tight_ptCut++;
     if ( Jet3_btagCSV > btagCSV_tight_cut  ) nBJet_tight_ptCut++;
     if ( Jet4_btagCSV > btagCSV_tight_cut  ) nBJet_tight_ptCut++;
     if ( Jet5_btagCSV > btagCSV_tight_cut  ) nBJet_tight_ptCut++;
     	
     //--------------------------------------------------------------------------
     // Fill final selection cuts
     //--------------------------------------------------------------------------

     if (!isOptimizationEnabled()) { 
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

     if ( do_final_selection ) { 
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
     
     if ( passed_preselection ) { 

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
       //// SIC ADD
       //FillUserTH1D( "MET_PAS_UnclUp"                    , PFMET_Type01XY_UnclUp_Pt                            , pileup_weight * gen_weight);
       //FillUserTH1D( "MET_PAS_UnclDown"                    , PFMET_Type01XY_UnclDown_Pt                            , pileup_weight * gen_weight);
       //// END SIC ADD
       FillUserTH1D( "MET_PDF"                    , PFMET_Type1XY_Pt                            , pileup_weight * gen_weight);
       FillUserTH1D( "METPhi_PAS"	          , PFMET_Type1XY_Phi                           , pileup_weight * gen_weight);   
       FillUserTH1D( "MET_Type01_PAS"             , PFMET_Type01_Pt                              , pileup_weight * gen_weight);
       FillUserTH1D( "MET_Type01_Phi_PAS"	  , PFMET_Type01_Phi                             , pileup_weight * gen_weight);   
       FillUserTH1D( "MET_Type1_PAS"             , PFMET_Type1_Pt                              , pileup_weight * gen_weight);
       FillUserTH1D( "MET_Type1_Phi_PAS"	  , PFMET_Type1_Phi                             , pileup_weight * gen_weight);   
       FillUserTH1D( "minMETPt1stEle_PAS"         , TMath::Min ( Ele1_Pt, PFMET_Type1XY_Pt  )   , pileup_weight * gen_weight);
       FillUserTH1D( "minMET01Pt1stEle_PAS"       , TMath::Min ( Ele1_Pt, PFMET_Type01_Pt    )   , pileup_weight * gen_weight);
       FillUserTH1D( "minMET1Pt1stEle_PAS"       , TMath::Min ( Ele1_Pt, PFMET_Type1_Pt    )   , pileup_weight * gen_weight);
       FillUserTH1D( "Pt1stJet_PAS"               , Jet1_Pt                                      , pileup_weight * gen_weight);
       FillUserTH1D( "Pt2ndJet_PAS"               , Jet2_Pt                                      , pileup_weight * gen_weight);
       FillUserTH1D( "Eta1stJet_PAS"              , Jet1_Eta                                     , pileup_weight * gen_weight);
       FillUserTH1D( "Eta2ndJet_PAS"              , Jet2_Eta                                     , pileup_weight * gen_weight);
       FillUserTH1D( "Phi1stJet_PAS"              , Jet1_Phi                                     , pileup_weight * gen_weight);
       FillUserTH1D( "Phi2ndJet_PAS"	          , Jet2_Phi                                     , pileup_weight * gen_weight);
       FillUserTH1D( "Mass1stJet_PAS"             , Jet1_Mass                                    , pileup_weight * gen_weight);
       FillUserTH1D( "Mass2ndJet_PAS"             , Jet2_Mass                                    , pileup_weight * gen_weight);
       FillUserTH1D( "CSV1stJet_PAS"              , Jet1_btagCSV                                 , pileup_weight * gen_weight);
       FillUserTH1D( "CSV2ndJet_PAS"              , Jet2_btagCSV                                 , pileup_weight * gen_weight);
       FillUserTH1D( "nMuon_PtCut_IDISO_PAS"      , nMuon_ptCut                                  , pileup_weight * gen_weight); 
       FillUserTH1D( "MTenu_PAS"                  , MT_Ele1MET                                   , pileup_weight * gen_weight);
       FillUserTH1D( "MT_PDF"                     , MT_Ele1MET                                   , pileup_weight * gen_weight);
       FillUserTH1D( "MTenu_Type01_PAS"           , MT_Ele1MET_Type01                            , pileup_weight * gen_weight);
       FillUserTH1D( "Ptenu_PAS"	          , Pt_Ele1MET                                   , pileup_weight * gen_weight);
       FillUserTH1D( "sTlep_PAS"                  , Ele1_Pt + PFMET_Type1XY_Pt                  , pileup_weight * gen_weight);
       FillUserTH1D( "sTlep_Type01_PAS"           , Ele1_Pt + PFMET_Type01_Pt                    , pileup_weight * gen_weight);
       FillUserTH1D( "sTlep_Type1_PAS"           , Ele1_Pt + PFMET_Type1_Pt                    , pileup_weight * gen_weight);
       FillUserTH1D( "sT_PAS"                     , sT_enujj                                     , pileup_weight * gen_weight);
       FillUserTH1D( "sT_PDF"                     , sT_enujj                                     , pileup_weight * gen_weight);
       FillUserTH1D( "sT_Type01_PAS"              , Ele1_Pt + PFMET_Type01_Pt + Jet1_Pt + Jet2_Pt, pileup_weight * gen_weight);
       FillUserTH1D( "sT_Type1_PAS"              , Ele1_Pt + PFMET_Type1_Pt + Jet1_Pt + Jet2_Pt, pileup_weight * gen_weight);
       FillUserTH1D( "sTjet_PAS"                  , Jet1_Pt + Jet2_Pt                            , pileup_weight * gen_weight);
       FillUserTH1D( "Mjj_PAS"	                  , M_j1j2                                       , pileup_weight * gen_weight);   
       FillUserTH1D( "DCotTheta1stEle_PAS"        , Ele1_DCotTheta                               , pileup_weight * gen_weight);
       FillUserTH1D( "Dist1stEle_PAS"             , Ele1_Dist                                    , pileup_weight * gen_weight);
       FillUserTH1D( "mDPhiJet1Jet2_PAS"          , mDPhi_Jet1Jet2                               , pileup_weight * gen_weight);
       FillUserTH1D( "mDPhi1stEleMET_PAS"         , mDPhi_METEle1                                , pileup_weight * gen_weight);
       FillUserTH1D( "mDPhi1stJetMET_PAS"         , mDPhi_METJet1                                , pileup_weight * gen_weight);
       FillUserTH1D( "mDPhi2ndJetMET_PAS"         , mDPhi_METJet2                                , pileup_weight * gen_weight); 
       FillUserTH1D( "mDPhi1stEleMET_Type01_PAS"  , mDPhi_METType01_Ele1                         , pileup_weight * gen_weight);
       FillUserTH1D( "mDPhi1stJetMET_Type01_PAS"  , mDPhi_METType01_Jet1                         , pileup_weight * gen_weight);
       FillUserTH1D( "mDPhi2ndJetMET_Type01_PAS"  , mDPhi_METType01_Jet2                         , pileup_weight * gen_weight); 
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

       likelihood_values.clear();
       likelihood_values.push_back ( Mej );  
       likelihood_values.push_back ( sT_enujj  );
       likelihood_values.push_back ( PFMET_Type1XY_Pt   );
       likelihood_values.push_back ( MT_Ele1MET );

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


       if ( MT_Ele1MET > 70 && MT_Ele1MET < 110 ){

         FillUserTH1D( "ProcessID_WWindow", ProcessID, pileup_weight * gen_weight );
         FillUserTH1D( "MTenu_70_110"      , MT_Ele1MET,  pileup_weight * gen_weight ) ;
         FillUserTH1D( "nJets_MTenu_70_110", nJet_ptCut,  pileup_weight * gen_weight ) ;

         if ( nJet_ptCut <= 3 ){ 
           FillUserTH1D(   "MTenu_70_110_Njet_lte3"         , MT_Ele1MET       ,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_lte3", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_lte3", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_lte3", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_70_110_Njet_lte3"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_70_110_Njet_lte3"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_70_110_Njet_lte3"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut <= 4 ){ 
           FillUserTH1D(   "MTenu_70_110_Njet_lte4", MT_Ele1MET,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_lte4", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_lte4", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_lte4", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_70_110_Njet_lte4"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_70_110_Njet_lte4"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_70_110_Njet_lte4"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut >= 4 ){ 
           FillUserTH1D(   "MTenu_70_110_Njet_gte4", MT_Ele1MET,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_gte4", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_gte4", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_gte4", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_70_110_Njet_gte4"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_70_110_Njet_gte4"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_70_110_Njet_gte4"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut >= 5 ){ 
           FillUserTH1D(   "MTenu_70_110_Njet_gte5", MT_Ele1MET,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_gte5", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_gte5", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_gte5", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_70_110_Njet_gte5"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_70_110_Njet_gte5"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_70_110_Njet_gte5"     , Mej              ,  pileup_weight * gen_weight ) ;
         }
       }

       if ( MT_Ele1MET_Type01 > 70 && MT_Ele1MET_Type01 < 110 ){

         FillUserTH1D( "MTenu_Type01_70_110"      , MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;
         FillUserTH1D( "nJets_MTenu_Type01_70_110", nJet_ptCut,  pileup_weight * gen_weight ) ;

         if ( nJet_ptCut <= 3 ){ 
           FillUserTH1D(   "MTenu_Type01_70_110_Njet_lte3", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_110_Njet_lte3", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_110_Njet_lte3", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_110_Njet_lte3", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_Type01_70_110_Njet_lte3"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_Type01_70_110_Njet_lte3"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_Type01_70_110_Njet_lte3"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut <= 4 ){ 
           FillUserTH1D(   "MTenu_Type01_70_110_Njet_lte4", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_110_Njet_lte4", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_110_Njet_lte4", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_110_Njet_lte4", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_Type01_70_110_Njet_lte4"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_Type01_70_110_Njet_lte4"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_Type01_70_110_Njet_lte4"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut >= 4 ){ 
           FillUserTH1D(   "MTenu_Type01_70_110_Njet_gte4", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_110_Njet_gte4", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_110_Njet_gte4", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_110_Njet_gte4", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_Type01_70_110_Njet_gte4"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_Type01_70_110_Njet_gte4"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_Type01_70_110_Njet_gte4"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut >= 5 ){ 
           FillUserTH1D(   "MTenu_Type01_70_110_Njet_gte5", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_110_Njet_gte5", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_110_Njet_gte5", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_110_Njet_gte5", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_Type01_70_110_Njet_gte5"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_Type01_70_110_Njet_gte5"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_Type01_70_110_Njet_gte5"     , Mej              ,  pileup_weight * gen_weight ) ;
         }
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

       if ( MT_Ele1MET_Type01 > 50 && MT_Ele1MET_Type01 < 110 ){

         FillUserTH1D( "MTenu_Type01_50_110"      , MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;
         FillUserTH1D( "nJets_MTenu_Type01_50_110", nJet_ptCut,  pileup_weight * gen_weight ) ;

         if ( nJet_ptCut <= 3 ){ 
           FillUserTH1D(   "MTenu_Type01_50_110_Njet_lte3", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_Type01_50_110_Njet_lte3", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_Type01_50_110_Njet_lte3", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_Type01_50_110_Njet_lte3", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_Type01_50_110_Njet_lte3"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_Type01_50_110_Njet_lte3"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_Type01_50_110_Njet_lte3"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut <= 4 ){ 
           FillUserTH1D(   "MTenu_Type01_50_110_Njet_lte4", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_Type01_50_110_Njet_lte4", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_Type01_50_110_Njet_lte4", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_Type01_50_110_Njet_lte4", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_Type01_50_110_Njet_lte4"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_Type01_50_110_Njet_lte4"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_Type01_50_110_Njet_lte4"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut >= 4 ){ 
           FillUserTH1D(   "MTenu_Type01_50_110_Njet_gte4", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_Type01_50_110_Njet_gte4", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_Type01_50_110_Njet_gte4", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_Type01_50_110_Njet_gte4", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_Type01_50_110_Njet_gte4"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_Type01_50_110_Njet_gte4"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_Type01_50_110_Njet_gte4"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut >= 5 ){ 
           FillUserTH1D(   "MTenu_Type01_50_110_Njet_gte5", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_Type01_50_110_Njet_gte5", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_Type01_50_110_Njet_gte5", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_Type01_50_110_Njet_gte5", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_Type01_50_110_Njet_gte5"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_Type01_50_110_Njet_gte5"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_Type01_50_110_Njet_gte5"     , Mej              ,  pileup_weight * gen_weight ) ;
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

       if ( MT_Ele1MET_Type01 > 70 && MT_Ele1MET_Type01 < 150 ){

         FillUserTH1D( "MTenu_Type01_70_150"      , MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;
         FillUserTH1D( "nJets_MTenu_Type01_70_150", nJet_ptCut,  pileup_weight * gen_weight ) ;

         if ( nJet_ptCut <= 3 ){ 
           FillUserTH1D(   "MTenu_Type01_70_150_Njet_lte3", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_150_Njet_lte3", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_150_Njet_lte3", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_150_Njet_lte3", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_Type01_70_150_Njet_lte3"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_Type01_70_150_Njet_lte3"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_Type01_70_150_Njet_lte3"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut <= 4 ){ 
           FillUserTH1D(   "MTenu_Type01_70_150_Njet_lte4", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_150_Njet_lte4", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_150_Njet_lte4", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_150_Njet_lte4", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_Type01_70_150_Njet_lte4"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_Type01_70_150_Njet_lte4"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_Type01_70_150_Njet_lte4"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut >= 4 ){ 
           FillUserTH1D(   "MTenu_Type01_70_150_Njet_gte4", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_150_Njet_gte4", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_150_Njet_gte4", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_150_Njet_gte4", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_Type01_70_150_Njet_gte4"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_Type01_70_150_Njet_gte4"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_Type01_70_150_Njet_gte4"     , Mej              ,  pileup_weight * gen_weight ) ;
         }

         if ( nJet_ptCut >= 5 ){ 
           FillUserTH1D(   "MTenu_Type01_70_150_Njet_gte5", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

           FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_150_Njet_gte5", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_150_Njet_gte5", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_150_Njet_gte5", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "MET_MTenu_Type01_70_150_Njet_gte5"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "sT_MTenu_Type01_70_150_Njet_gte5"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
           FillUserTH1D(   "Mej_MTenu_Type01_70_150_Njet_gte5"     , Mej              ,  pileup_weight * gen_weight ) ;
         }
       }

       // no b tags
       //-------------------------------------------------------------------------- 
       if (nBJet_medium_ptCut==0) {

         if ( MT_Ele1MET > 70 && MT_Ele1MET < 110 ){

           FillUserTH1D( "ProcessID_WWindow", ProcessID, pileup_weight * gen_weight );
           FillUserTH1D( "MTenu_70_110_noBtaggedJets"      , MT_Ele1MET,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "nJets_MTenu_70_110_noBtaggedJets", nJet_ptCut,  pileup_weight * gen_weight ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_70_110_Njet_lte3_noBtaggedJets"         , MT_Ele1MET       ,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_lte3_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_lte3_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_lte3_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_110_Njet_lte3_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_110_Njet_lte3_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_110_Njet_lte3_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_70_110_Njet_lte4_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_lte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_lte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_lte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_110_Njet_lte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_110_Njet_lte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_110_Njet_lte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_70_110_Njet_gte4_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_gte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_gte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_gte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_110_Njet_gte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_110_Njet_gte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_110_Njet_gte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_70_110_Njet_gte5_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_gte5_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_gte5_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_gte5_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_110_Njet_gte5_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_110_Njet_gte5_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_110_Njet_gte5_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }
         }

         if ( MT_Ele1MET_Type01 > 70 && MT_Ele1MET_Type01 < 110 ){

           FillUserTH1D( "MTenu_Type01_70_110_noBtaggedJets"      , MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "nJets_MTenu_Type01_70_110_noBtaggedJets", nJet_ptCut,  pileup_weight * gen_weight ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_Type01_70_110_Njet_lte3_noBtaggedJets", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_110_Njet_lte3_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_110_Njet_lte3_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_110_Njet_lte3_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_110_Njet_lte3_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_110_Njet_lte3_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_110_Njet_lte3_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_Type01_70_110_Njet_lte4_noBtaggedJets", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_110_Njet_lte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_110_Njet_lte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_110_Njet_lte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_110_Njet_lte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_110_Njet_lte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_110_Njet_lte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_Type01_70_110_Njet_gte4_noBtaggedJets", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_110_Njet_gte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_110_Njet_gte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_110_Njet_gte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_110_Njet_gte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_110_Njet_gte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_110_Njet_gte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_Type01_70_110_Njet_gte5_noBtaggedJets", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_110_Njet_gte5_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_110_Njet_gte5_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_110_Njet_gte5_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_110_Njet_gte5_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_110_Njet_gte5_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_110_Njet_gte5_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }
         }

         if ( MT_Ele1MET > 50 && MT_Ele1MET < 110 ){

           FillUserTH1D( "ProcessID_WWindow", ProcessID, pileup_weight * gen_weight );
           FillUserTH1D( "MTenu_50_110_noBtaggedJets"      , MT_Ele1MET,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets", nJet_ptCut,  pileup_weight * gen_weight ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_lte3_noBtaggedJets"         , MT_Ele1MET       ,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_lte3_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_lte3_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_lte3_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_lte3_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_lte3_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_lte3_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_lte4_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_lte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_lte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_lte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_lte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_lte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_lte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_gte4_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_gte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_gte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_gte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_gte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_gte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_gte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_gte5_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_gte5_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_gte5_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_gte5_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_gte5_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_gte5_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_gte5_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }
         }

         if ( MT_Ele1MET_Type01 > 50 && MT_Ele1MET_Type01 < 110 ){

           FillUserTH1D( "MTenu_Type01_50_110_noBtaggedJets"      , MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "nJets_MTenu_Type01_50_110_noBtaggedJets", nJet_ptCut,  pileup_weight * gen_weight ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_Type01_50_110_Njet_lte3_noBtaggedJets", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_50_110_Njet_lte3_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_50_110_Njet_lte3_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_50_110_Njet_lte3_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_50_110_Njet_lte3_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_50_110_Njet_lte3_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_50_110_Njet_lte3_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_Type01_50_110_Njet_lte4_noBtaggedJets", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_50_110_Njet_lte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_50_110_Njet_lte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_50_110_Njet_lte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_50_110_Njet_lte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_50_110_Njet_lte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_50_110_Njet_lte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_Type01_50_110_Njet_gte4_noBtaggedJets", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_50_110_Njet_gte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_50_110_Njet_gte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_50_110_Njet_gte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_50_110_Njet_gte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_50_110_Njet_gte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_50_110_Njet_gte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_Type01_50_110_Njet_gte5_noBtaggedJets", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_50_110_Njet_gte5_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_50_110_Njet_gte5_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_50_110_Njet_gte5_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_50_110_Njet_gte5_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_50_110_Njet_gte5_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_50_110_Njet_gte5_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }
         }

         if ( MT_Ele1MET > 70 && MT_Ele1MET < 150 ){

           FillUserTH1D( "ProcessID_WWindow", ProcessID, pileup_weight * gen_weight );
           FillUserTH1D( "MTenu_70_150_noBtaggedJets"      , MT_Ele1MET,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "nJets_MTenu_70_150_noBtaggedJets", nJet_ptCut,  pileup_weight * gen_weight ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_lte3_noBtaggedJets"         , MT_Ele1MET       ,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_lte3_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_lte3_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_lte3_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_lte3_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_lte3_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_lte3_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_lte4_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_lte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_lte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_lte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_lte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_lte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_lte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_gte4_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_gte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_gte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_gte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_gte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_gte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_gte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_gte5_noBtaggedJets", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_gte5_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_gte5_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_gte5_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_gte5_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_gte5_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_gte5_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }
         }

         if ( MT_Ele1MET_Type01 > 70 && MT_Ele1MET_Type01 < 150 ){

           FillUserTH1D( "MTenu_Type01_70_150_noBtaggedJets"      , MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "nJets_MTenu_Type01_70_150_noBtaggedJets", nJet_ptCut,  pileup_weight * gen_weight ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_Type01_70_150_Njet_lte3_noBtaggedJets", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_150_Njet_lte3_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_150_Njet_lte3_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_150_Njet_lte3_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_150_Njet_lte3_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_150_Njet_lte3_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_150_Njet_lte3_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_Type01_70_150_Njet_lte4_noBtaggedJets", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_150_Njet_lte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_150_Njet_lte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_150_Njet_lte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_150_Njet_lte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_150_Njet_lte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_150_Njet_lte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_Type01_70_150_Njet_gte4_noBtaggedJets", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_150_Njet_gte4_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_150_Njet_gte4_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_150_Njet_gte4_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_150_Njet_gte4_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_150_Njet_gte4_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_150_Njet_gte4_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_Type01_70_150_Njet_gte5_noBtaggedJets", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_150_Njet_gte5_noBtaggedJets", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_150_Njet_gte5_noBtaggedJets", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_150_Njet_gte5_noBtaggedJets", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_150_Njet_gte5_noBtaggedJets"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_150_Njet_gte5_noBtaggedJets"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_150_Njet_gte5_noBtaggedJets"     , Mej              ,  pileup_weight * gen_weight ) ;
           }
         }
       }

       // at least 1 b tag
       //-------------------------------------------------------------------------- 
       if (nBJet_medium_ptCut>=1) {

         if ( MT_Ele1MET > 70 && MT_Ele1MET < 110 ){

           FillUserTH1D( "ProcessID_WWindow", ProcessID, pileup_weight * gen_weight );
           FillUserTH1D( "MTenu_70_110_gteOneBtaggedJet"      , MT_Ele1MET,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "nJets_MTenu_70_110_gteOneBtaggedJet", nJet_ptCut,  pileup_weight * gen_weight ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_70_110_Njet_lte3_gteOneBtaggedJet"         , MT_Ele1MET       ,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_lte3_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_lte3_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_lte3_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_110_Njet_lte3_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_110_Njet_lte3_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_110_Njet_lte3_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_70_110_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_lte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_lte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_lte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_110_Njet_lte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_110_Njet_lte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_110_Njet_lte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_70_110_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_gte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_gte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_gte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_110_Njet_gte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_110_Njet_gte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_110_Njet_gte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_70_110_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_110_Njet_gte5_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_110_Njet_gte5_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_110_Njet_gte5_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_110_Njet_gte5_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_110_Njet_gte5_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_110_Njet_gte5_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }
         }

         if ( MT_Ele1MET_Type01 > 70 && MT_Ele1MET_Type01 < 110 ){

           FillUserTH1D( "MTenu_Type01_70_110_gteOneBtaggedJet"      , MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "nJets_MTenu_Type01_70_110_gteOneBtaggedJet", nJet_ptCut,  pileup_weight * gen_weight ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }
         }

         if ( MT_Ele1MET > 50 && MT_Ele1MET < 110 ){

           FillUserTH1D( "ProcessID_WWindow", ProcessID, pileup_weight * gen_weight );
           FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet"      , MT_Ele1MET,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet", nJet_ptCut,  pileup_weight * gen_weight ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_lte3_gteOneBtaggedJet"         , MT_Ele1MET       ,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_lte3_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_lte3_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_lte3_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_lte3_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_lte3_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_lte3_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_lte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_lte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_lte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_lte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_lte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_lte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_gte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_gte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_gte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_gte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_gte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_gte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_gte5_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_gte5_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_gte5_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_50_110_Njet_gte5_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_50_110_Njet_gte5_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_50_110_Njet_gte5_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }
         }

         if ( MT_Ele1MET_Type01 > 50 && MT_Ele1MET_Type01 < 110 ){

           FillUserTH1D( "MTenu_Type01_50_110_gteOneBtaggedJet"      , MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "nJets_MTenu_Type01_50_110_gteOneBtaggedJet", nJet_ptCut,  pileup_weight * gen_weight ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }
         }

         if ( MT_Ele1MET > 70 && MT_Ele1MET < 150 ){

           FillUserTH1D( "ProcessID_WWindow", ProcessID, pileup_weight * gen_weight );
           FillUserTH1D( "MTenu_70_150_gteOneBtaggedJet"      , MT_Ele1MET,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "nJets_MTenu_70_150_gteOneBtaggedJet", nJet_ptCut,  pileup_weight * gen_weight ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_lte3_gteOneBtaggedJet"         , MT_Ele1MET       ,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_lte3_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_lte3_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_lte3_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_lte3_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_lte3_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_lte3_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_lte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_lte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_lte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_lte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_lte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_lte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_gte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_gte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_gte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_gte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_gte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_gte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_70_150_Njet_gte5_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_70_150_Njet_gte5_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_70_150_Njet_gte5_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_70_150_Njet_gte5_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_70_150_Njet_gte5_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_70_150_Njet_gte5_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }
         }

         if ( MT_Ele1MET_Type01 > 70 && MT_Ele1MET_Type01 < 150 ){

           FillUserTH1D( "MTenu_Type01_70_150_gteOneBtaggedJet"      , MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "nJets_MTenu_Type01_70_150_gteOneBtaggedJet", nJet_ptCut,  pileup_weight * gen_weight ) ;

           if ( nJet_ptCut <= 3 ){ 
             FillUserTH1D(   "MTenu_Type01_70_150_Njet_lte3_gteOneBtaggedJet", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_150_Njet_lte3_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_150_Njet_lte3_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_150_Njet_lte3_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_150_Njet_lte3_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_150_Njet_lte3_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_150_Njet_lte3_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut <= 4 ){ 
             FillUserTH1D(   "MTenu_Type01_70_150_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_150_Njet_lte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_150_Njet_lte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_150_Njet_lte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_150_Njet_lte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_150_Njet_lte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_150_Njet_lte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 4 ){ 
             FillUserTH1D(   "MTenu_Type01_70_150_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_150_Njet_gte4_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_150_Njet_gte4_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_150_Njet_gte4_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_150_Njet_gte4_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_150_Njet_gte4_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_150_Njet_gte4_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }

           if ( nJet_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_Type01_70_150_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET_Type01,  pileup_weight * gen_weight ) ;

             FillUserTH1D(   "Pt1stEle_MTenu_Type01_70_150_Njet_gte5_gteOneBtaggedJet", Ele1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt1stJet_MTenu_Type01_70_150_Njet_gte5_gteOneBtaggedJet", Jet1_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Pt2ndJet_MTenu_Type01_70_150_Njet_gte5_gteOneBtaggedJet", Jet2_Pt          ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "MET_MTenu_Type01_70_150_Njet_gte5_gteOneBtaggedJet"     , PFMET_Type1XY_Pt,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "sT_MTenu_Type01_70_150_Njet_gte5_gteOneBtaggedJet"      , sT_enujj         ,  pileup_weight * gen_weight ) ;
             FillUserTH1D(   "Mej_MTenu_Type01_70_150_Njet_gte5_gteOneBtaggedJet"     , Mej              ,  pileup_weight * gen_weight ) ;
           }
         }
       }

       //-------------------------------------------------------------------------- 
       // Final selection plots
       //-------------------------------------------------------------------------- 
       if ( !isOptimizationEnabled() ) { 

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
             sprintf(plot_name, "sTlep_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , Ele1_Pt + PFMET_Type1XY_Pt    , pileup_weight * gen_weight );
             sprintf(plot_name, "sTjet_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , Jet1_Pt + Jet2_Pt              , pileup_weight * gen_weight );
             sprintf(plot_name, "Energy1stEle_LQ%d"      , lq_mass ); FillUserTH1D( plot_name , Ele1_Energy                    , pileup_weight * gen_weight );
             sprintf(plot_name, "Pt1stEle_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Ele1_Pt                        , pileup_weight * gen_weight );
             sprintf(plot_name, "Pt1stJet_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Jet1_Pt                        , pileup_weight * gen_weight );
             sprintf(plot_name, "Pt2ndJet_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Jet2_Pt                        , pileup_weight * gen_weight );
             sprintf(plot_name, "sTfrac_Jet1_LQ%d"       , lq_mass ); FillUserTH1D( plot_name , Jet1_Pt / sT_enujj              , pileup_weight * gen_weight );
             sprintf(plot_name, "sTfrac_Jet2_LQ%d"       , lq_mass ); FillUserTH1D( plot_name , Jet2_Pt / sT_enujj              , pileup_weight * gen_weight );
             sprintf(plot_name, "sTfrac_Ele1_LQ%d"       , lq_mass ); FillUserTH1D( plot_name , Ele1_Pt / sT_enujj              , pileup_weight * gen_weight );
             sprintf(plot_name, "sTfrac_MET_LQ%d"        , lq_mass ); FillUserTH1D( plot_name , PFMET_Type1XY_Pt / sT_enujj              , pileup_weight * gen_weight );
             sprintf(plot_name, "sTfrac_Jet_LQ%d"        , lq_mass ); FillUserTH1D( plot_name , ( Jet1_Pt + Jet2_Pt ) / sT_enujj, pileup_weight * gen_weight );
             sprintf(plot_name, "sTfrac_Lep_LQ%d"        , lq_mass ); FillUserTH1D( plot_name , ( Ele1_Pt + PFMET_Type1XY_Pt ) / sT_enujj, pileup_weight * gen_weight );
           }

           if ( do_eleQuality_plots ) { 
             sprintf(plot_name, "BeamSpotDXY_1stEle_LQ%d"        , lq_mass );   FillUserTH1D(plot_name,  Ele1_BeamSpotDXY               , pileup_weight * gen_weight ); 
             sprintf(plot_name, "Classif_1stEle_LQ%d"            , lq_mass );   FillUserTH1D(plot_name,  Ele1_Classif                   , pileup_weight * gen_weight ); 
             sprintf(plot_name, "CorrIsolation_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele1_CorrIsolation             , pileup_weight * gen_weight ); 
             sprintf(plot_name, "DeltaEtaTrkSC_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele1_DeltaEtaTrkSC             , pileup_weight * gen_weight ); 
             sprintf(plot_name, "DeltaPhiTrkSC_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele1_DeltaPhiTrkSC             , pileup_weight * gen_weight ); 
             sprintf(plot_name, "E1x5OverE5x5_1stEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  Ele1_E1x5OverE5x5              , pileup_weight * gen_weight ); 
             sprintf(plot_name, "E2x5OverE5x5_1stEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  Ele1_E2x5OverE5x5              , pileup_weight * gen_weight ); 
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
               sprintf(plot_name, "SigmaIEtaIEta_Barrel_1stEle_LQ%d", lq_mass ); FillUserTH1D( plot_name , Ele1_SigmaIEtaIEta , pileup_weight * gen_weight    ); 
             }
             else if ( fabs(Ele1_Eta) > eleEta_end1_min && fabs(Ele2_Eta) < eleEta_end2_max ){
               sprintf(plot_name, "SigmaEtaEta_Endcap_1stEle_LQ%d"  , lq_mass ); FillUserTH1D( plot_name , Ele1_SigmaEtaEta   , pileup_weight * gen_weight    ); 
               sprintf(plot_name, "SigmaIEtaIEta_Endcap_1stEle_LQ%d", lq_mass ); FillUserTH1D( plot_name , Ele1_SigmaIEtaIEta , pileup_weight * gen_weight    ); 
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

         if ( passedCut ("ST_LQ300") && passedCut ("Mej_LQ300") ){
           FillUserTH1D("MET_LQ300_NoMETCut", PFMET_Type1XY_Pt , gen_weight * pileup_weight );
           FillUserTH1D("Mej_LQ300_NoMETCut", Mej               , gen_weight * pileup_weight );
           FillUserTH1D("ST_LQ300_NoMETCut" , sT_enujj          , gen_weight * pileup_weight ); 
           FillUserTH1D("MT_LQ300_NoMETCut" , MT_Ele1MET        , gen_weight * pileup_weight ); 
         }
       } // End do final selection
     } // End do preselection
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
