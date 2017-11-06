#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
// for fake rate
#include "include/QCDFakeRate.h"

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile){}

analysisClass::~analysisClass(){}

void analysisClass::Loop()
{
   std::cout << "analysisClass::Loop() begins" <<std::endl;   
   bool do_extra_finalSelection_plots = true;

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
   
   //--------------------------------------------------------------------------
   // Decide which plots to save (default is to save everything)
   //--------------------------------------------------------------------------
   
   fillSkim                         ( !true  ) ;
   fillAllPreviousCuts              ( !true  ) ;
   fillAllOtherCuts                 (  true  ) ;
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

  //--------------------------------------------------------------------------
  // QCD Fake Rate loading part
  //--------------------------------------------------------------------------
  std::string qcdFileName = getPreCutString1("QCDFakeRateFileName");
  //FIXME TODO eventually remove hardcoding of graph name
  QCDFakeRate qcdFakeRate(qcdFileName);
   
   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------

   CreateUserTH1D( "nElectron_PAS"            , 5   , -0.5    , 4.5      );
   CreateUserTH1D( "nMuon_PAS"                , 5   , -0.5    , 4.5      );
   CreateUserTH1D( "nJet_PAS"                 , 11  , -0.5    , 10.5     );
   CreateUserTH1D( "nJet_PASandFrancesco"     , 11  , -0.5    , 10.5     );
   CreateUserTH1D( "Pt1stEle_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "PtHeep1stEle_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "SCEt1stEle_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta1stEle_PAS"	      , 100 , -5      , 5	 ); 
   CreateUserTH1D( "SCEta1stEle_PAS"	      , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stEle_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Pt1stEle_Barrel_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "Pt1stEle_Endcap1_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "Pt1stEle_Endcap2_PAS"	      , 100 , 0       , 1000     ); 
   // muon kinematics
   CreateUserTH1D( "Pt1stMuon_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta1stMuon_PAS"	      , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stMuon_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Pt2ndMuon_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta2ndMuon_PAS"	      , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi2ndMuon_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   //
   CreateUserTH1D( "Charge1stEle_PAS"	      , 2   , -1.0001 , 1.0001	 ); 
   CreateUserTH1D( "MET_PAS"                  , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "METPhi_PAS"		      , 60  , -3.1416 , +3.1416	 ); 
   //CreateUserTH1D( "MET_Type01_PAS"           , 200 , 0       , 1000	 ); 
   //CreateUserTH1D( "MET_Type01_Phi_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "minMETPt1stEle_PAS"       , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "Pt1stJet_PAS"             , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Pt2ndJet_PAS"             , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Eta1stJet_PAS"            , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Eta2ndJet_PAS"            , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stJet_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Phi2ndJet_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Mass1stJet_PAS"           , 100 , 0       , 500      );
   CreateUserTH1D( "Mass2ndJet_PAS"           , 100 , 0       , 500      );
   CreateUserTH1D( "CISV1stJet_PAS"            , 200 , 0       , 1.0	 ); 
   CreateUserTH1D( "CISV2ndJet_PAS"            , 200 , 0       , 1.0	 ); 
   CreateUserTH1D( "nMuon_PtCut_IDISO_PAS"    , 16  , -0.5    , 15.5	 ); 
   CreateUserTH1D( "MTenu_PAS"                , 400 , 0       , 2000	 ); 
   //CreateUserTH1D( "MTenu_Type01_PAS"         , 400 , 0       , 2000	 ); 
   CreateUserTH1D( "MT_charged_enu_PAS"       , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "MT_type1_enu_PAS"         , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "Ptenu_PAS"		      , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTlep_PAS"                , 200 , 0       , 2000	 ); 
   //CreateUserTH1D( "sTlep_Type01_PAS"         , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTjet_PAS"                , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PAS"                   , 300 , 0       , 3000	 ); 
   //CreateUserTH1D( "sT_Type01_PAS"            , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mjj_PAS"		      , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mej1_PAS"                 , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mej2_PAS"                 , 200 , 0       , 2000	 );
   CreateUserTH1D( "Mej_PAS"                  , 200 , 0       , 2000	 );  
   CreateUserTH1D( "MTjnu_PAS"                , 200 , 0       , 1000     );
   CreateUserTH1D( "DCotTheta1stEle_PAS"      , 100 , 0.0     , 1.0      );
   CreateUserTH1D( "Dist1stEle_PAS"           , 100 , 0.0     , 1.0      );  
   CreateUserTH1D( "DR_Ele1Jet1_PAS"	      , 100 , 0       , 10       ); 
   CreateUserTH1D( "DR_Ele1Jet2_PAS"	      , 100 , 0       , 10       ); 
   CreateUserTH1D( "DR_Jet1Jet2_PAS"	      , 100 , 0       , 10       ); 
   CreateUserTH1D( "minDR_EleJet_PAS"         , 100 , 0       , 10       ); 
   CreateUserTH1D( "mDPhi1stEleMET_PAS"       , 100 , 0.      ,  3.14159 );
   CreateUserTH1D( "mDPhi1stJetMET_PAS"       , 100 , 0.      ,  3.14159 );
   CreateUserTH1D( "mDPhi2ndJetMET_PAS"       , 100 , 0.      ,  3.14159 );
   //CreateUserTH1D( "mDPhi1stEleMET_Type01_PAS", 100 , 0.      , 3.14159  );
   //CreateUserTH1D( "mDPhi1stJetMET_Type01_PAS", 100 , 0.      , 3.14159  );
   //CreateUserTH1D( "mDPhi2ndJetMET_Type01_PAS", 100 , 0.      , 3.14159  );

   CreateUserTH1D( "Pt1stEle_Pt40to45_EtaGT2p1", 150, 35., 50. );

   CreateUserTH1D( "MT_GoodVtxLTE3_PAS"       , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MT_GoodVtxGTE4_LTE8_PAS"  , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MT_GoodVtxGTE9_LTE15_PAS" , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MT_GoodVtxGTE16_PAS"      , 200 , 0.      ,  1000    );

   CreateUserTH1D( "GeneratorWeight"       , 200 , -2.0    , 2.0      );
   CreateUserTH1D( "PileupWeight"          , 200 , -2.0    , 2.0      );

   CreateUserTH1D( "nVertex_PAS"           ,    101   , -0.5   , 100.5	 ) ; 

   CreateUserTH1D( "MTCharged_GoodVtxLTE3_PAS"       , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTCharged_GoodVtxGTE4_LTE8_PAS"  , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTCharged_GoodVtxGTE9_LTE15_PAS" , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTCharged_GoodVtxGTE16_PAS"      , 200 , 0.      ,  1000    );

   CreateUserTH1D( "MTType1_GoodVtxLTE3_PAS"       , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTType1_GoodVtxGTE4_LTE8_PAS"  , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTType1_GoodVtxGTE9_LTE15_PAS" , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTType1_GoodVtxGTE16_PAS"      , 200 , 0.      ,  1000    );
   
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
   
   CreateUserTH1D( "MTenu_50_110", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte3", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte5", 240, 40, 160 );

   //CreateUserTH1D( "MTenu_Type01_50_110", 200, 40, 140 );
   //CreateUserTH1D( "nJets_MTenu_Type01_50_110"    , 20 , -0.5, 19.5 );
   //CreateUserTH1D( "MTenu_Type01_50_110_Njet_gte4", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_50_110_Njet_lte3", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_50_110_Njet_lte4", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_50_110_Njet_gte5", 200, 40, 140 );

   CreateUserTH1D( "MTenu_70_110", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_70_110"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte3", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte4", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte5", 240, 40, 160 );
   //CreateUserTH1D( "MTenu_Type01_70_110", 200, 40, 140 );
   //CreateUserTH1D( "nJets_MTenu_Type01_70_110"    , 20 , -0.5, 19.5 );
   //CreateUserTH1D( "MTenu_Type01_70_110_Njet_gte4", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_70_110_Njet_lte3", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_70_110_Njet_lte4", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_70_110_Njet_gte5", 200, 40, 140 );

   CreateUserTH1D( "MTenu_50_110_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_noBtaggedJets"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte3_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte5_noBtaggedJets", 240, 40, 160 );
   //CreateUserTH1D( "MTenu_Type01_50_110_noBtaggedJets", 200, 40, 140 );
   //CreateUserTH1D( "nJets_MTenu_Type01_50_110_noBtaggedJets"    , 20 , -0.5, 19.5 );
   //CreateUserTH1D( "MTenu_Type01_50_110_Njet_gte4_noBtaggedJets", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_50_110_Njet_lte3_noBtaggedJets", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_50_110_Njet_lte4_noBtaggedJets", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_50_110_Njet_gte5_noBtaggedJets", 200, 40, 140 );

   CreateUserTH1D( "MTenu_70_110_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_70_110_noBtaggedJets"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte3_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte5_noBtaggedJets", 240, 40, 160 );
   //CreateUserTH1D( "MTenu_Type01_70_110_noBtaggedJets", 200, 40, 140 );
   //CreateUserTH1D( "nJets_MTenu_Type01_70_110_noBtaggedJets"    , 20 , -0.5, 19.5 );
   //CreateUserTH1D( "MTenu_Type01_70_110_Njet_gte4_noBtaggedJets", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_70_110_Njet_lte3_noBtaggedJets", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_70_110_Njet_lte4_noBtaggedJets", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_70_110_Njet_gte5_noBtaggedJets", 200, 40, 140 );

   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte3_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte5_gteOneBtaggedJet", 240, 40, 160 );
   //CreateUserTH1D( "MTenu_Type01_50_110_gteOneBtaggedJet", 200, 40, 140 );
   //CreateUserTH1D( "nJets_MTenu_Type01_50_110_gteOneBtaggedJet"    , 20 , -0.5, 19.5 );
   //CreateUserTH1D( "MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet", 200, 40, 140 );

   CreateUserTH1D( "MTenu_70_110_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "nJets_MTenu_70_110_gteOneBtaggedJet"    , 20 , -0.5, 19.5 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte3_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_lte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_70_110_Njet_gte5_gteOneBtaggedJet", 240, 40, 160 );
   //CreateUserTH1D( "MTenu_Type01_70_110_gteOneBtaggedJet", 200, 40, 140 );
   //CreateUserTH1D( "nJets_MTenu_Type01_70_110_gteOneBtaggedJet"    , 20 , -0.5, 19.5 );
   //CreateUserTH1D( "MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet", 200, 40, 140 );
   //CreateUserTH1D( "MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet", 200, 40, 140 );

   CreateUserTH1D( "Eta1stJet_PASand2Jet"  , 100 , -5 , 5 ); 
   CreateUserTH1D( "Eta1stJet_PASand3Jet"  , 100 , -5 , 5 ); 
   CreateUserTH1D( "Eta1stJet_PASand4Jet"  , 100 , -5 , 5 ); 
   CreateUserTH1D( "Eta1stJet_PASand5Jet"  , 100 , -5 , 5 ); 

   CreateUserTH1D( "Mej_Barrel_Presel"                  , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mej_Endcap1_Presel"                  , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mej_Endcap2_Presel"                  , 200 , 0       , 2000	 ); 

   //--------------------------------------------------------------------------
   // Final selection plots
   //--------------------------------------------------------------------------
   
   char plot_name[100];
  
   if (!isOptimizationEnabled() ) { 
     for (int i_lq_mass = 0; i_lq_mass < n_lq_mass ; ++i_lq_mass ) { 
       int lq_mass = LQ_MASS[i_lq_mass];
       sprintf(plot_name, "MTenu_LQ%d" , lq_mass ); CreateUserTH1D ( plot_name, 400 , 0 , 2000 ); 
       sprintf(plot_name, "Mej_LQ%d"   , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000 );
       sprintf(plot_name, "ST_LQ%d"    , lq_mass ); CreateUserTH1D ( plot_name, 300 , 0 , 3000 );
       sprintf(plot_name, "Mej_Barrel_LQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000	); 
       sprintf(plot_name, "Mej_Endcap1_LQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000	); 
       sprintf(plot_name, "Mej_Endcap2_LQ%d"      , lq_mass ); CreateUserTH1D ( plot_name, 200 , 0 , 2000	); 
       if ( do_extra_finalSelection_plots ) { 
         //sprintf(plot_name, "sTfrac_Jet1_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
         //sprintf(plot_name, "sTfrac_Jet2_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
         //sprintf(plot_name, "sTfrac_Ele1_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
         //sprintf(plot_name, "sTfrac_MET_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
         //sprintf(plot_name, "sTfrac_Jet_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
         //sprintf(plot_name, "sTfrac_Lep_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
         //sprintf(plot_name, "Charge1stEle_LQ%d"	       , lq_mass ); CreateUserTH1D( plot_name ,      3 ,  -1.5   , 1.5      );
         //sprintf(plot_name, "Energy1stEle_LQ%d"	       , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
         sprintf(plot_name, "Pt1stEle_LQ%d"	       , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
         sprintf(plot_name, "PtHeep1stEle_LQ%d"	       , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
         sprintf(plot_name, "SCEt1stEle_LQ%d"	       , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
         sprintf(plot_name, "Eta1stEle_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
         sprintf(plot_name, "SCEta1stEle_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
         sprintf(plot_name, "Phi1stEle_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , +3.1416  ); 
         sprintf(plot_name, "Pt1stJet_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
         sprintf(plot_name, "Pt2ndJet_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
         sprintf(plot_name, "Eta1stJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
         sprintf(plot_name, "Eta2ndJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
         sprintf(plot_name, "Phi1stJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , 3.1416   ); 
         sprintf(plot_name, "Phi2ndJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , 3.1416   ); 
         //sprintf(plot_name, "sTlep_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
         //sprintf(plot_name, "sTjet_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
         //sprintf(plot_name, "Me1j1_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
         //sprintf(plot_name, "Me1j2_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
         //sprintf(plot_name, "nVertex_LQ%d"             , lq_mass ); CreateUserTH1D( plot_name ,    101 , -0.5    ,  100.5   ); 
         //sprintf(plot_name, "MeeVsST_LQ%d"             , lq_mass ); CreateUserTH2D( plot_name ,    200 ,  0.0, 2000, 200, 0, 2000);
         // muon kinematics
         sprintf(plot_name, "Pt1stMuon_LQ%d"	       , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
         sprintf(plot_name, "Eta1stMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
         sprintf(plot_name, "Phi1stMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , +3.1416  ); 
         sprintf(plot_name, "Pt2ndMuon_LQ%d"	       , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
         sprintf(plot_name, "Eta2ndMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
         sprintf(plot_name, "Phi2ndMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , +3.1416  ); 
       }
     }
   }
   
   //--------------------------------------------------------------------------
   // Loop over the chain
   //--------------------------------------------------------------------------

   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntries();
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
     // Find the right prescale for this event
     //--------------------------------------------------------------------------
     
     int min_prescale = 0;
     int passTrigger  = 0;

     if ( LooseEle1_hltPhotonPt > 0.0 ) { 
       if ( H_Photon22   > 0.1 && LooseEle1_hltPhotonPt >= 22.  && LooseEle1_hltPhotonPt < 30. ) { passTrigger = 1; min_prescale = H_Photon22  ; } 
       if ( H_Photon30   > 0.1 && LooseEle1_hltPhotonPt >= 30.  && LooseEle1_hltPhotonPt < 36. ) { passTrigger = 1; min_prescale = H_Photon30  ; } 
       if ( H_Photon36   > 0.1 && LooseEle1_hltPhotonPt >= 36.  && LooseEle1_hltPhotonPt < 50. ) { passTrigger = 1; min_prescale = H_Photon36  ; } 
       if ( H_Photon50   > 0.1 && LooseEle1_hltPhotonPt >= 50.  && LooseEle1_hltPhotonPt < 75. ) { passTrigger = 1; min_prescale = H_Photon50  ; } 
       if ( H_Photon75   > 0.1 && LooseEle1_hltPhotonPt >= 75.  && LooseEle1_hltPhotonPt < 90. ) { passTrigger = 1; min_prescale = H_Photon75  ; } 
       if ( H_Photon90   > 0.1 && LooseEle1_hltPhotonPt >= 90.  && LooseEle1_hltPhotonPt < 120.) { passTrigger = 1; min_prescale = H_Photon90  ; } 
       if ( H_Photon120  > 0.1 && LooseEle1_hltPhotonPt >= 120. && LooseEle1_hltPhotonPt < 175.) { passTrigger = 1; min_prescale = H_Photon120 ; } 
       if ( H_Photon175  > 0.1 && LooseEle1_hltPhotonPt >= 175.) { passTrigger = 1; min_prescale = H_Photon175 ; } 
     }

     //--------------------------------------------------------------------------
     // What kind of event is this?
     //   - Barrel
     //   - Endcap 1 (eta < 2.0)
     //   - Endcap 2 (eta > 2.0) 
     //--------------------------------------------------------------------------
     
     bool ele1_isBarrel  = false;
     bool ele1_isEndcap1 = false;
     bool ele1_isEndcap2 = false;
     bool ele2_isBarrel  = false;
     bool ele2_isEndcap1 = false;
     bool ele2_isEndcap2 = false;
     
     if( fabs( LooseEle1_Eta  ) < eleEta_bar )        ele1_isBarrel  = true;
     if( fabs( LooseEle1_Eta  ) > eleEta_end1_min &&
	 fabs( LooseEle1_Eta  ) < eleEta_end1_max )   ele1_isEndcap1 = true;
     if( fabs( LooseEle1_Eta  ) > eleEta_end2_min &&
	 fabs( LooseEle1_Eta  ) < eleEta_end2_max )   ele1_isEndcap2 = true;

     if( fabs( LooseEle2_Eta  ) < eleEta_bar )        ele2_isBarrel  = true;
     if( fabs( LooseEle2_Eta  ) > eleEta_end1_min &&
	 fabs( LooseEle2_Eta  ) < eleEta_end1_max )   ele2_isEndcap1 = true;
     if( fabs( LooseEle2_Eta  ) > eleEta_end2_min &&
	 fabs( LooseEle2_Eta  ) < eleEta_end2_max )   ele2_isEndcap2 = true;

     //--------------------------------------------------------------------------
     // Make this a QCD fake rate calculation
     //--------------------------------------------------------------------------
     //float fakeRate1 = qcdFakeRate.GetFakeRate(LooseEle1_SCEta,LooseEle1_PtHeep);
     float fakeRate1 = qcdFakeRate.GetFakeRate(LooseEle1_SCEta,LooseEle1_SCEnergy/cosh(LooseEle1_SCEta));
     
     //--------------------------------------------------------------------------
     // Finally have the effective fake rate
     //--------------------------------------------------------------------------

     double fakeRateEffective  = fakeRate1;
     double eFakeRateEffective = 0.0; //FIXME eFakeRate1;
     
     //--------------------------------------------------------------------------
     // Calculate some variables:
     //--------------------------------------------------------------------------
     
     TLorentzVector loose_ele1, loose_ele2, jet1, jet2, met;
     loose_ele1.SetPtEtaPhiM ( LooseEle1_Pt , LooseEle1_Eta , LooseEle1_Phi , 0.0 );
     loose_ele2.SetPtEtaPhiM ( LooseEle2_Pt , LooseEle2_Eta , LooseEle2_Phi , 0.0 );
     jet1.SetPtEtaPhiM       ( JetLooseEle1_Pt, JetLooseEle1_Eta, JetLooseEle1_Phi, 0.0 );
     jet2.SetPtEtaPhiM       ( JetLooseEle2_Pt, JetLooseEle2_Eta, JetLooseEle2_Phi, 0.0 );
     met.SetPtEtaPhiM        ( PFMET_Type1XY_Pt         , 0.0             , PFMET_Type1XY_Phi         , 0.0 );

     TLorentzVector e1met = loose_ele1 + met;
     TLorentzVector j1j2 = jet1 + jet2;
     TLorentzVector e1j1 = loose_ele1 + jet1;
     TLorentzVector e1j2 = loose_ele1 + jet2;

     M_e1j1 = e1j1.M();
     M_e1j2 = e1j2.M();

     Pt_Ele1MET = e1met.Pt();
     MT_Ele1MET = sqrt(2 * LooseEle1_Pt * PFMET_Type1XY_Pt  * (1 - cos(loose_ele1.DeltaPhi ( met) ) ) );

     mDPhi_METEle1= fabs(loose_ele1.DeltaPhi ( met ));
     mDPhi_METJet1= fabs(jet1.DeltaPhi ( met ));
     mDPhi_METJet2= fabs(jet2.DeltaPhi ( met ));

     sT_enujj = LooseEle1_Pt + PFMET_Type1XY_Pt + JetLooseEle1_Pt + JetLooseEle2_Pt ;
     
     DR_Ele1Jet1 = loose_ele1.DeltaR ( jet1 ) ;
     DR_Ele1Jet2 = loose_ele1.DeltaR ( jet2 ) ;

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "Reweighting"              , 1                       , min_prescale * fakeRateEffective ); 
     fillVariableWithValue(   "PassJSON"                 , passedJSON              , min_prescale * fakeRateEffective ); 
     									          
     // HLT variable							           
     fillVariableWithValue(   "PassHLT"                  , passTrigger             , min_prescale * fakeRateEffective );
     
     				      
     //--------------------------------------------------------------------------
     // Fill noise filters
     //--------------------------------------------------------------------------
     // see: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
     // we filled these at skim time
     fillVariableWithValue( "PassGlobalTightHalo2016Filter" , PassGlobalTightHalo2016Filter  , fakeRateEffective * min_prescale );
     fillVariableWithValue( "PassGoodVertices"	             , PassGoodVertices               , fakeRateEffective * min_prescale );
     fillVariableWithValue( "PassHBHENoiseFilter"	         , PassHBHENoiseFilter            , fakeRateEffective * min_prescale );
     fillVariableWithValue( "PassHBHENoiseIsoFilter"	       , PassHBHENoiseIsoFilter         , fakeRateEffective * min_prescale );
     fillVariableWithValue( "PassBadEESupercrystalFilter"   , PassBadEESupercrystalFilter    , fakeRateEffective * min_prescale );
     fillVariableWithValue( "PassEcalDeadCellTrigPrim"      , PassEcalDeadCellTrigPrim       , fakeRateEffective * min_prescale );
     fillVariableWithValue( "PassChargedCandidateFilter"    , PassChargedCandidateFilter     , fakeRateEffective * min_prescale );
     fillVariableWithValue( "PassBadPFMuonFilter"           , PassBadPFMuonFilter            , fakeRateEffective * min_prescale );
     //

     // Muon variables ( for veto ) 					      
     fillVariableWithValue(   "nMuon"                    , nMuon_ptCut             , min_prescale * fakeRateEffective );
			                                      		                
     // 1st Electron variables				      		                
     fillVariableWithValue(   "nEle"                     , nLooseEle_ptCut       , min_prescale * fakeRateEffective ); 
     fillVariableWithValue(   "Ele1_SCEt"              , LooseEle1_SCEnergy/cosh(LooseEle1_SCEta)      , min_prescale * fakeRateEffective );
     fillVariableWithValue(   "Ele1_Eta"                 , LooseEle1_Eta         , min_prescale * fakeRateEffective );
     fillVariableWithValue(   "Ele1_IsBarrel"            , ele1_isBarrel         , min_prescale * fakeRateEffective );
     fillVariableWithValue(   "MTenu"                    , MT_Ele1MET            , min_prescale * fakeRateEffective );
									           
     // MET variables	                                      		           
     fillVariableWithValue(   "MET"                      , PFMET_Type1XY_Pt                  , min_prescale * fakeRateEffective );
     fillVariableWithValue(   "mDeltaPhiMETEle"          , mDPhi_METEle1           , min_prescale * fakeRateEffective );
     									           
     // 1st JET variables                                     		           
     fillVariableWithValue(   "nJet"                     , nJetLooseEle_ptCut      , min_prescale * fakeRateEffective );
			
     double MT_Jet1MET, MT_Jet2MET, MT_Ele1Jet1, MT_Ele1Jet2, MT_Ele1MET_Type01;
     double mDPhi_METType01_Ele1, mDPhi_METType01_Jet1, mDPhi_METType01_Jet2;

     // alternate METs
     if ( nLooseEle_store > 0 ) {
       TVector2 v_ele;
       TVector2 v_MET_Type01;
       v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       v_ele.SetMagPhi( LooseEle1_Pt, LooseEle1_Phi );
       mDPhi_METType01_Ele1 = fabs(v_MET_Type01.DeltaPhi ( v_ele ));
       float deltaphi = v_MET_Type01.DeltaPhi(v_ele);
       MT_Ele1MET_Type01 = sqrt ( 2 * LooseEle1_Pt * PFMET_Type01_Pt * ( 1 - cos ( deltaphi ) ) );
     }
				           
     // 1st JET variables                                     		           
     if ( nJetLooseEle_store > 0 ) { 						           
       fillVariableWithValue( "Jet1_Pt"                  , JetLooseEle1_Pt         , min_prescale * fakeRateEffective );
       fillVariableWithValue( "Jet1_Eta"                 , JetLooseEle1_Eta        , min_prescale * fakeRateEffective );
       fillVariableWithValue( "mDeltaPhiMET1stJet"       , mDPhi_METJet1           , min_prescale * fakeRateEffective );

       TVector2 v_MET;
       TVector2 v_jet;
       TVector2 v_MET_Type01;
       v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       v_MET.SetMagPhi( PFMET_Type1XY_Pt , PFMET_Type1XY_Phi  );
       v_jet.SetMagPhi( JetLooseEle1_Pt, JetLooseEle1_Phi );
       mDPhi_METType01_Jet1 = fabs(v_MET_Type01.DeltaPhi ( v_jet ));
       float deltaphi = v_MET.DeltaPhi(v_jet);
       MT_Jet1MET = sqrt ( 2 * JetLooseEle1_Pt * PFMET_Type1XY_Pt * ( 1 - cos ( deltaphi ) ) );
     }									           
     									           
     // 2nd JET variables                                     		           
     if ( nJetLooseEle_store > 1 ) { 	                                      	           
       fillVariableWithValue( "Jet2_Pt"                  , JetLooseEle2_Pt         , min_prescale * fakeRateEffective );
       fillVariableWithValue( "Jet2_Eta"                 , JetLooseEle2_Eta        , min_prescale * fakeRateEffective );
       fillVariableWithValue( "ST"                       , sT_enujj                , min_prescale * fakeRateEffective );

       TVector2 v_MET;
       TVector2 v_jet;
       TVector2 v_MET_Type01;
       v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       v_MET.SetMagPhi( PFMET_Type1XY_Pt , PFMET_Type1XY_Phi  );
       v_jet.SetMagPhi( JetLooseEle2_Pt, JetLooseEle2_Phi );
       mDPhi_METType01_Jet2 = fabs(v_MET_Type01.DeltaPhi ( v_jet ));
       float deltaphi = v_MET.DeltaPhi(v_jet);
       MT_Jet2MET = sqrt ( 2 * JetLooseEle2_Pt * PFMET_Type1XY_Pt * ( 1 - cos ( deltaphi ) ) );
     }

     // 3rd JET variables 
     // if ( nJetLooseEle_store > 2 ) {
     //   fillVariableWithValue( "Jet3_Pt"                  , JetLooseEle3_Pt         , min_prescale * fakeRateEffective );
     //   fillVariableWithValue( "Jet3_Eta"                 , JetLooseEle3_Eta        , min_prescale * fakeRateEffective );
     // }
     
     // 1 electron, 1 jet variables 
     if ( nLooseEle_store > 0 && nJetLooseEle_store > 0 ) { 
       fillVariableWithValue ( "DR_Ele1Jet1"             , DR_Ele1Jet1             , min_prescale * fakeRateEffective );


       TVector2 v_ele;
       TVector2 v_jet1;
       TVector2 v_MET_Type01;
       v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       v_ele .SetMagPhi ( LooseEle1_Pt, LooseEle1_Phi );
       v_jet1.SetMagPhi ( JetLooseEle1_Pt, JetLooseEle1_Phi );
       float deltaphi = v_ele.DeltaPhi ( v_jet1 );
       MT_Ele1Jet2 = sqrt ( 2 * JetLooseEle1_Pt * LooseEle1_Pt * ( 1 - cos ( deltaphi ) ) );

     }

     // 1 electron, 2 jet variables 
     if ( nLooseEle_store > 0 && nJetLooseEle_store > 1 ) { 
       fillVariableWithValue ( "DR_Ele1Jet2"             , DR_Ele1Jet2             , min_prescale * fakeRateEffective );
       
       TVector2 v_ele;
       TVector2 v_jet2;
       TVector2 v_MET_Type01;
       v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       v_ele .SetMagPhi ( LooseEle1_Pt, LooseEle1_Phi );
       v_jet2.SetMagPhi ( JetLooseEle2_Pt, JetLooseEle2_Phi );
       float deltaphi = v_ele.DeltaPhi ( v_jet2 );
       MT_Ele1Jet2 = sqrt ( 2 * JetLooseEle2_Pt * LooseEle1_Pt * ( 1 - cos ( deltaphi ) ) );

     }
     
     double MT_JetMET;
     double Mej;
     
     if ( fabs ( MT_Jet1MET - MT_Ele1Jet2 ) < fabs( MT_Jet2MET - MT_Ele1Jet1 )){
       MT_JetMET = MT_Jet1MET;
       Mej = M_e1j2;
     } else { 
       MT_JetMET = MT_Jet2MET;
       Mej = M_e1j1;
     }	 

     
     // Optimization variables
     fillVariableWithValue( "ST_opt"   , sT_enujj   , min_prescale * fakeRateEffective );
     fillVariableWithValue( "Mej_min_opt"  , Mej        , min_prescale * fakeRateEffective );
     fillVariableWithValue( "MT_opt", MT_Ele1MET , min_prescale * fakeRateEffective );
     fillVariableWithValue( "MET_opt" , PFMET_Type1XY_Pt   , min_prescale * fakeRateEffective );

     // Dummy variables
     fillVariableWithValue ("preselection",1, min_prescale * fakeRateEffective );

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
     // for event weights, we follow: https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods#1c_Event_reweighting_using_scale

     double weightZeroBJets = 1.0;
     double weightAtLeastOneBJet = 1.0;
     double btagCISV_loose_cut  = 0.5426;
     double btagCISV_medium_cut = 0.8484;
     double btagCISV_tight_cut  = 0.9535;
       
     int nBJet_loose_ptCut  = 0;
     int nBJet_medium_ptCut = 0;
     int nBJet_tight_ptCut  = 0;
     
     if ( JetLooseEle1_btagCISV > btagCISV_loose_cut  ) nBJet_loose_ptCut++;
     if ( JetLooseEle2_btagCISV > btagCISV_loose_cut  ) nBJet_loose_ptCut++;
     if ( JetLooseEle3_btagCISV > btagCISV_loose_cut  ) nBJet_loose_ptCut++;
     if ( JetLooseEle4_btagCISV > btagCISV_loose_cut  ) nBJet_loose_ptCut++;
     if ( JetLooseEle5_btagCISV > btagCISV_loose_cut  ) nBJet_loose_ptCut++;
     
     // medium used below for control regions
     if ( JetLooseEle1_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut++;
     if ( JetLooseEle2_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut++;
     if ( JetLooseEle3_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut++;
     if ( JetLooseEle4_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut++;
     if ( JetLooseEle5_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut++;
     
     if ( JetLooseEle1_btagCISV > btagCISV_tight_cut  ) nBJet_tight_ptCut++;
     if ( JetLooseEle2_btagCISV > btagCISV_tight_cut  ) nBJet_tight_ptCut++;
     if ( JetLooseEle3_btagCISV > btagCISV_tight_cut  ) nBJet_tight_ptCut++;
     if ( JetLooseEle4_btagCISV > btagCISV_tight_cut  ) nBJet_tight_ptCut++;
     if ( JetLooseEle5_btagCISV > btagCISV_tight_cut  ) nBJet_tight_ptCut++;
     	
     //--------------------------------------------------------------------------
     // Fill final selection cuts
     //--------------------------------------------------------------------------

     char cut_name[100];
     if(!isOptimizationEnabled()) {
       for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
         int lq_mass = LQ_MASS[i_lq_mass];
         sprintf(cut_name, "MT_LQ%d" , lq_mass ); fillVariableWithValue( cut_name , MT_Ele1MET , min_prescale * fakeRateEffective );
         sprintf(cut_name, "Mej_LQ%d"   , lq_mass ); fillVariableWithValue( cut_name , Mej        , min_prescale * fakeRateEffective );
         sprintf(cut_name, "ST_LQ%d"    , lq_mass ); fillVariableWithValue( cut_name , sT_enujj   , min_prescale * fakeRateEffective );
         sprintf(cut_name, "MET_LQ%d"  , lq_mass ); fillVariableWithValue( cut_name , PFMET_Type1XY_Pt   , min_prescale * fakeRateEffective );
       }
     }
     
     
     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------
     
     bool passed_preselection = passedAllPreviousCuts("preselection");


     //--------------------------------------------------------------------------
     // Did we pass any final selection cuts?
     //--------------------------------------------------------------------------

     if(!isOptimizationEnabled() ) {
       passed_vector.clear();
       for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
         int lq_mass = LQ_MASS[i_lq_mass];
         sprintf(cut_name, "Mej_LQ%d", lq_mass );
         bool decision = bool ( passedAllPreviousCuts(cut_name) && passedCut (cut_name));
         passed_vector.push_back (decision);
       }
     }
     
     if ( passed_preselection ) { 

       //--------------------------------------------------------------------------
       // Fill skim tree, if necessary
       //--------------------------------------------------------------------------
       
       double JetLooseEle1_Mass, JetLooseEle2_Mass;
       TLorentzVector jetm1, jetm2;
       jetm1.SetPtEtaPhiE       ( JetLooseEle1_Pt, JetLooseEle1_Eta, JetLooseEle1_Phi, JetLooseEle1_Energy );
       jetm2.SetPtEtaPhiE       ( JetLooseEle2_Pt, JetLooseEle2_Eta, JetLooseEle2_Phi, JetLooseEle2_Energy );
       JetLooseEle1_Mass = jetm1.M();
       JetLooseEle2_Mass = jetm2.M();

       double min_DR_EleJet = 999.0;
       double DR_Ele1Jet3 = 999.0;
       if ( nJetLooseEle_store > 2 ) {
         TLorentzVector ele1, jet3;
         ele1.SetPtEtaPhiM ( LooseEle1_Pt, LooseEle1_Eta, LooseEle1_Phi, 0.0 );
         jet3.SetPtEtaPhiM ( JetLooseEle3_Pt, JetLooseEle3_Eta, JetLooseEle3_Phi, 0.0 );
         DR_Ele1Jet3 = ele1.DeltaR ( jet3 ) ;
       }

       if ( DR_Ele1Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet1;
       if ( DR_Ele1Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet2;
       if ( nJetLooseEle_store > 2 ) {
         if ( DR_Ele1Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet3;
       }

       double sT_enujj_Type01 = LooseEle1_Pt + PFMET_Type01_Pt + JetLooseEle1_Pt + JetLooseEle2_Pt;

       FillUserTH1D( "nElectron_PAS"              , nLooseEle_ptCut                                  , min_prescale * fakeRateEffective); 
       FillUserTH1D( "nMuon_PAS"                  , nMuon_ptCut                                      , min_prescale * fakeRateEffective); 
       FillUserTH1D( "Pt1stEle_PAS"	              , LooseEle1_Pt                                     , min_prescale * fakeRateEffective); 
       FillUserTH1D( "PtHeep1stEle_PAS"	          , LooseEle1_PtHeep                                 , min_prescale * fakeRateEffective); 
       FillUserTH1D( "SCEta1stEle_PAS"	          , LooseEle1_SCEta                                  , min_prescale * fakeRateEffective); 
       FillUserTH1D( "SCEt1stEle_PAS"	            , LooseEle1_SCEnergy/cosh(LooseEle1_SCEta)         , min_prescale * fakeRateEffective); 
       FillUserTH1D( "Eta1stEle_PAS"	            , LooseEle1_Eta                                    , min_prescale * fakeRateEffective);
       FillUserTH1D( "Phi1stEle_PAS"	            , LooseEle1_Phi                                    , min_prescale * fakeRateEffective);
       // muon kinematics
       FillUserTH1D( "Pt1stMuon_PAS"	            , Muon1_Pt                                     , min_prescale * fakeRateEffective); 
       FillUserTH1D( "Eta1stMuon_PAS"	            , Muon1_Eta                                    , min_prescale * fakeRateEffective);
       FillUserTH1D( "Phi1stMuon_PAS"	            , Muon1_Phi                                    , min_prescale * fakeRateEffective);
       FillUserTH1D( "Pt2ndMuon_PAS"	            , Muon2_Pt                                     , min_prescale * fakeRateEffective); 
       FillUserTH1D( "Eta2ndMuon_PAS"	            , Muon2_Eta                                    , min_prescale * fakeRateEffective);
       FillUserTH1D( "Phi2ndMuon_PAS"	            , Muon2_Phi                                    , min_prescale * fakeRateEffective);
       //
       FillUserTH1D( "Charge1stEle_PAS"           , LooseEle1_Charge                                 , min_prescale * fakeRateEffective);   
       FillUserTH1D( "MET_PAS"                    , PFMET_Type1XY_Pt                                , min_prescale * fakeRateEffective);
       FillUserTH1D( "METPhi_PAS"	                , PFMET_Type1XY_Phi                               , min_prescale * fakeRateEffective);   
       //FillUserTH1D( "MET_Type01_PAS"             , PFMET_Type01_Pt                                  , min_prescale * fakeRateEffective);
       //FillUserTH1D( "MET_Type01_Phi_PAS"	        , PFMET_Type01_Phi                                 , min_prescale * fakeRateEffective);   
       FillUserTH1D( "minMETPt1stEle_PAS"         , TMath::Min ( LooseEle1_Pt, PFMET_Type1XY_Pt  )  , min_prescale * fakeRateEffective);
       FillUserTH1D( "Pt1stJet_PAS"               , JetLooseEle1_Pt                                  , min_prescale * fakeRateEffective);
       FillUserTH1D( "Pt2ndJet_PAS"               , JetLooseEle2_Pt                                  , min_prescale * fakeRateEffective);
       FillUserTH1D( "Eta1stJet_PAS"              , JetLooseEle1_Eta                                 , min_prescale * fakeRateEffective);
       FillUserTH1D( "Eta2ndJet_PAS"              , JetLooseEle2_Eta                                 , min_prescale * fakeRateEffective);
       FillUserTH1D( "Phi1stJet_PAS"              , JetLooseEle1_Phi                                 , min_prescale * fakeRateEffective);
       FillUserTH1D( "Phi2ndJet_PAS"	          , JetLooseEle2_Phi                                 , min_prescale * fakeRateEffective);
       FillUserTH1D( "Mass1stJet_PAS"             , JetLooseEle1_Mass                                , min_prescale * fakeRateEffective);
       FillUserTH1D( "Mass2ndJet_PAS"             , JetLooseEle2_Mass                                , min_prescale * fakeRateEffective);
       FillUserTH1D( "CISV1stJet_PAS"              , JetLooseEle1_btagCISV                             , min_prescale * fakeRateEffective);
       FillUserTH1D( "CISV2ndJet_PAS"              , JetLooseEle2_btagCISV                             , min_prescale * fakeRateEffective);
       FillUserTH1D( "nMuon_PtCut_IDISO_PAS"      , nMuon_ptCut                                      , min_prescale * fakeRateEffective); 
       FillUserTH1D( "MTenu_PAS"                  , MT_Ele1MET                                       , min_prescale * fakeRateEffective);
       //FillUserTH1D( "MTenu_Type01_PAS"           , MT_Ele1MET_Type01                                , min_prescale * fakeRateEffective);
       FillUserTH1D( "Ptenu_PAS"	          , Pt_Ele1MET                                       , min_prescale * fakeRateEffective);
       FillUserTH1D( "sTlep_PAS"                  , LooseEle1_Pt + PFMET_Type1XY_Pt                 , min_prescale * fakeRateEffective);
       //FillUserTH1D( "sTlep_Type01_PAS"           , LooseEle1_Pt + PFMET_Type01_Pt                   , min_prescale * fakeRateEffective);
       FillUserTH1D( "sTjet_PAS"                  , JetLooseEle1_Pt + JetLooseEle2_Pt                , min_prescale * fakeRateEffective);
       FillUserTH1D( "sT_PAS"                     , sT_enujj                                         , min_prescale * fakeRateEffective);
       //FillUserTH1D( "sT_Type01_PAS"              , sT_enujj_Type01                                  , min_prescale * fakeRateEffective);
       FillUserTH1D( "Mjj_PAS"	                  , M_j1j2                                           , min_prescale * fakeRateEffective);   
       FillUserTH1D( "DCotTheta1stEle_PAS"        , LooseEle1_DCotTheta                              , min_prescale * fakeRateEffective);
       FillUserTH1D( "Dist1stEle_PAS"             , LooseEle1_Dist                                   , min_prescale * fakeRateEffective);
       FillUserTH1D( "mDPhi1stEleMET_PAS"         , mDPhi_METEle1                                    , min_prescale * fakeRateEffective);
       FillUserTH1D( "mDPhi1stJetMET_PAS"         , mDPhi_METJet1                                    , min_prescale * fakeRateEffective);
       FillUserTH1D( "mDPhi2ndJetMET_PAS"         , mDPhi_METJet2                                    , min_prescale * fakeRateEffective); 
       //FillUserTH1D( "mDPhi1stEleMET_Type01_PAS"  , mDPhi_METType01_Ele1                             , min_prescale * fakeRateEffective);
       //FillUserTH1D( "mDPhi1stJetMET_Type01_PAS"  , mDPhi_METType01_Jet1                             , min_prescale * fakeRateEffective);
       //FillUserTH1D( "mDPhi2ndJetMET_Type01_PAS"  , mDPhi_METType01_Jet2                             , min_prescale * fakeRateEffective); 
       FillUserTH1D( "Mej1_PAS"                   , M_e1j1                                           , min_prescale * fakeRateEffective);
       FillUserTH1D( "Mej2_PAS"                   , M_e1j2                                           , min_prescale * fakeRateEffective);
       FillUserTH1D( "Mej_PAS"                    , Mej                                              , min_prescale * fakeRateEffective);
       FillUserTH1D( "MTjnu_PAS"                  , MT_JetMET                                        , min_prescale * fakeRateEffective);
       FillUserTH1D( "DR_Ele1Jet1_PAS"	          , DR_Ele1Jet1                                      , min_prescale * fakeRateEffective);
       FillUserTH1D( "DR_Ele1Jet2_PAS"	          , DR_Ele1Jet2                                      , min_prescale * fakeRateEffective);
       FillUserTH1D( "DR_Jet1Jet2_PAS"	          , DR_Jet1Jet2                                      , min_prescale * fakeRateEffective);
       FillUserTH1D( "minDR_EleJet_PAS"           , min_DR_EleJet                                    , min_prescale * fakeRateEffective);
       FillUserTH1D( "nVertex_PAS"                , nVertex                                          , min_prescale * fakeRateEffective);
       FillUserTH1D( "nJet_PAS"                   , nJetLooseEle_ptCut                               , min_prescale * fakeRateEffective);
       FillUserTH1D( "GeneratorWeight"            , -1.0             );
       FillUserTH1D( "PileupWeight"               , -1.0             );
       
       if ( Pt_Ele1MET         > 200. && 
           JetLooseEle1_Pt    > 200. && 
           JetLooseEle1_Mass  > 50.  ) {
         FillUserTH1D("nJet_PASandFrancesco", nJetLooseEle_ptCut, min_prescale * fakeRateEffective);
       }

       if ( nJetLooseEle_ptCut == 2 ) FillUserTH1D( "Eta1stJet_PASand2Jet", JetLooseEle1_Eta, min_prescale * fakeRateEffective);
       if ( nJetLooseEle_ptCut == 3 ) FillUserTH1D( "Eta1stJet_PASand3Jet", JetLooseEle1_Eta, min_prescale * fakeRateEffective);
       if ( nJetLooseEle_ptCut == 4 ) FillUserTH1D( "Eta1stJet_PASand4Jet", JetLooseEle1_Eta, min_prescale * fakeRateEffective);
       if ( nJetLooseEle_ptCut == 5 ) FillUserTH1D( "Eta1stJet_PASand5Jet", JetLooseEle1_Eta, min_prescale * fakeRateEffective);

       if ( LooseEle1_Pt > 40 && LooseEle1_Pt < 45 && fabs ( LooseEle1_Eta ) > 2.1 ){
         FillUserTH1D( "Pt1stEle_Pt40to45_EtaGT2p1", LooseEle1_Pt, min_prescale * fakeRateEffective);
       }

       if ( MT_Ele1MET > 50 && MT_Ele1MET < 110 ){

         FillUserTH1D( "MTenu_50_110"      , MT_Ele1MET, min_prescale * fakeRateEffective );
         FillUserTH1D( "nJets_MTenu_50_110", nJetLooseEle_ptCut, min_prescale * fakeRateEffective );

         if ( nJetLooseEle_ptCut <= 3 ){
           FillUserTH1D(   "MTenu_50_110_Njet_lte3", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_lte3", LooseEle1_Pt          ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_lte3", JetLooseEle1_Pt          ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_lte3", JetLooseEle2_Pt          ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "MET_MTenu_50_110_Njet_lte3"     , PFMET_Type1XY_Pt,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "sT_MTenu_50_110_Njet_lte3"      , sT_enujj         ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Mej_MTenu_50_110_Njet_lte3"     , Mej              ,  min_prescale * fakeRateEffective ) ;
         }

         if ( nJetLooseEle_ptCut <= 4 ){
           FillUserTH1D(   "MTenu_50_110_Njet_lte4", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_lte4", LooseEle1_Pt          ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_lte4", JetLooseEle1_Pt          ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_lte4", JetLooseEle2_Pt          ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "MET_MTenu_50_110_Njet_lte4"     , PFMET_Type1XY_Pt,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "sT_MTenu_50_110_Njet_lte4"      , sT_enujj         ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Mej_MTenu_50_110_Njet_lte4"     , Mej              ,  min_prescale * fakeRateEffective ) ;
         }

         if ( nJetLooseEle_ptCut >= 4 ){
           FillUserTH1D(   "MTenu_50_110_Njet_gte4", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_gte4", LooseEle1_Pt          ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_gte4", JetLooseEle1_Pt          ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_gte4", JetLooseEle2_Pt          ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "MET_MTenu_50_110_Njet_gte4"     , PFMET_Type1XY_Pt,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "sT_MTenu_50_110_Njet_gte4"      , sT_enujj         ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Mej_MTenu_50_110_Njet_gte4"     , Mej              ,  min_prescale * fakeRateEffective ) ;
         }

         if ( nJetLooseEle_ptCut >= 5 ){ 
           FillUserTH1D(   "MTenu_50_110_Njet_gte5", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Pt1stEle_MTenu_50_110_Njet_gte5", LooseEle1_Pt          ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Pt1stJet_MTenu_50_110_Njet_gte5", JetLooseEle1_Pt          ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Pt2ndJet_MTenu_50_110_Njet_gte5", JetLooseEle2_Pt          ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "MET_MTenu_50_110_Njet_gte5"     , PFMET_Type1XY_Pt,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "sT_MTenu_50_110_Njet_gte5"      , sT_enujj         ,  min_prescale * fakeRateEffective ) ;
           FillUserTH1D(   "Mej_MTenu_50_110_Njet_gte5"     , Mej              ,  min_prescale * fakeRateEffective ) ;
         }
       }

       //if ( MT_Ele1MET_Type01 > 50 && MT_Ele1MET_Type01 < 110 ){

       //  FillUserTH1D( "MTenu_Type01_50_110"      , MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
       //  FillUserTH1D( "nJets_MTenu_Type01_50_110", nJetLooseEle_ptCut,  min_prescale * fakeRateEffective ) ;

       //  if ( nJetLooseEle_ptCut <= 3 ){ 
       //    FillUserTH1D(   "MTenu_Type01_50_110_Njet_lte3", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
       //  }

       //  if ( nJetLooseEle_ptCut <= 4 ){ 
       //    FillUserTH1D(   "MTenu_Type01_50_110_Njet_lte4", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
       //  }

       //  if ( nJetLooseEle_ptCut >= 4 ){ 
       //    FillUserTH1D(   "MTenu_Type01_50_110_Njet_gte4", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
       //  }

       //  if ( nJetLooseEle_ptCut >= 5 ){ 
       //    FillUserTH1D(   "MTenu_Type01_50_110_Njet_gte5", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
       //  }
       //}
       
       if ( MT_Ele1MET > 70 && MT_Ele1MET < 110 ){

         FillUserTH1D( "MTenu_70_110"      , MT_Ele1MET, min_prescale * fakeRateEffective );
         FillUserTH1D( "nJets_MTenu_70_110", nJetLooseEle_ptCut, min_prescale * fakeRateEffective );

         if ( nJetLooseEle_ptCut <= 3 ){
           FillUserTH1D(   "MTenu_70_110_Njet_lte3", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
         }

         if ( nJetLooseEle_ptCut <= 4 ){
           FillUserTH1D(   "MTenu_70_110_Njet_lte4", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
         }

         if ( nJetLooseEle_ptCut >= 4 ){
           FillUserTH1D(   "MTenu_70_110_Njet_gte4", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
         }

         if ( nJetLooseEle_ptCut >= 5 ){ 
           FillUserTH1D(   "MTenu_70_110_Njet_gte5", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
         }
       }

       //if ( MT_Ele1MET_Type01 > 70 && MT_Ele1MET_Type01 < 110 ){

       //  FillUserTH1D( "MTenu_Type01_70_110"      , MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
       //  FillUserTH1D( "nJets_MTenu_Type01_70_110", nJetLooseEle_ptCut,  min_prescale * fakeRateEffective ) ;

       //  if ( nJetLooseEle_ptCut <= 3 ){ 
       //    FillUserTH1D(   "MTenu_Type01_70_110_Njet_lte3", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
       //  }

       //  if ( nJetLooseEle_ptCut <= 4 ){ 
       //    FillUserTH1D(   "MTenu_Type01_70_110_Njet_lte4", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
       //  }

       //  if ( nJetLooseEle_ptCut >= 4 ){ 
       //    FillUserTH1D(   "MTenu_Type01_70_110_Njet_gte4", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
       //  }

       //  if ( nJetLooseEle_ptCut >= 5 ){ 
       //    FillUserTH1D(   "MTenu_Type01_70_110_Njet_gte5", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
       //  }
       //}


       // no b-tagged jets
       if (nBJet_medium_ptCut==0) {

         if ( MT_Ele1MET > 50 && MT_Ele1MET < 110 ){

           FillUserTH1D( "MTenu_50_110_noBtaggedJets"      , MT_Ele1MET, min_prescale * fakeRateEffective );
           FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets", nJetLooseEle_ptCut, min_prescale * fakeRateEffective );

           if ( nJetLooseEle_ptCut <= 3 ){
             FillUserTH1D(   "MTenu_50_110_Njet_lte3_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut <= 4 ){
             FillUserTH1D(   "MTenu_50_110_Njet_lte4_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 4 ){
             FillUserTH1D(   "MTenu_50_110_Njet_gte4_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_gte5_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }
         }

         //if ( MT_Ele1MET_Type01 > 50 && MT_Ele1MET_Type01 < 110 ){

         //  FillUserTH1D( "MTenu_Type01_50_110_noBtaggedJets"      , MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  FillUserTH1D( "nJets_MTenu_Type01_50_110_noBtaggedJets", nJetLooseEle_ptCut,  min_prescale * fakeRateEffective ) ;

         //  if ( nJetLooseEle_ptCut <= 3 ){ 
         //    FillUserTH1D(   "MTenu_Type01_50_110_Njet_lte3_noBtaggedJets", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }

         //  if ( nJetLooseEle_ptCut <= 4 ){ 
         //    FillUserTH1D(   "MTenu_Type01_50_110_Njet_lte4_noBtaggedJets", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }

         //  if ( nJetLooseEle_ptCut >= 4 ){ 
         //    FillUserTH1D(   "MTenu_Type01_50_110_Njet_gte4_noBtaggedJets", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }

         //  if ( nJetLooseEle_ptCut >= 5 ){ 
         //    FillUserTH1D(   "MTenu_Type01_50_110_Njet_gte5_noBtaggedJets", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }
         //}

         if ( MT_Ele1MET > 70 && MT_Ele1MET < 110 ){

           FillUserTH1D( "MTenu_70_110_noBtaggedJets"      , MT_Ele1MET, min_prescale * fakeRateEffective );
           FillUserTH1D( "nJets_MTenu_70_110_noBtaggedJets", nJetLooseEle_ptCut, min_prescale * fakeRateEffective );

           if ( nJetLooseEle_ptCut <= 3 ){
             FillUserTH1D(   "MTenu_70_110_Njet_lte3_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut <= 4 ){
             FillUserTH1D(   "MTenu_70_110_Njet_lte4_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 4 ){
             FillUserTH1D(   "MTenu_70_110_Njet_gte4_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_70_110_Njet_gte5_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }
         }

         //if ( MT_Ele1MET_Type01 > 70 && MT_Ele1MET_Type01 < 110 ){

         //  FillUserTH1D( "MTenu_Type01_70_110_noBtaggedJets"      , MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  FillUserTH1D( "nJets_MTenu_Type01_70_110_noBtaggedJets", nJetLooseEle_ptCut,  min_prescale * fakeRateEffective ) ;

         //  if ( nJetLooseEle_ptCut <= 3 ){ 
         //    FillUserTH1D(   "MTenu_Type01_70_110_Njet_lte3_noBtaggedJets", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }

         //  if ( nJetLooseEle_ptCut <= 4 ){ 
         //    FillUserTH1D(   "MTenu_Type01_70_110_Njet_lte4_noBtaggedJets", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }

         //  if ( nJetLooseEle_ptCut >= 4 ){ 
         //    FillUserTH1D(   "MTenu_Type01_70_110_Njet_gte4_noBtaggedJets", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }

         //  if ( nJetLooseEle_ptCut >= 5 ){ 
         //    FillUserTH1D(   "MTenu_Type01_70_110_Njet_gte5_noBtaggedJets", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }
         //}
       }

       // at least one b-tagged jet
       if (nBJet_medium_ptCut>=1) {

         if ( MT_Ele1MET > 50 && MT_Ele1MET < 110 ){

           FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet"      , MT_Ele1MET, min_prescale * fakeRateEffective );
           FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet", nJetLooseEle_ptCut, min_prescale * fakeRateEffective );

           if ( nJetLooseEle_ptCut <= 3 ){
             FillUserTH1D(   "MTenu_50_110_Njet_lte3_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut <= 4 ){
             FillUserTH1D(   "MTenu_50_110_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 4 ){
             FillUserTH1D(   "MTenu_50_110_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_50_110_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }
         }

         //if ( MT_Ele1MET_Type01 > 50 && MT_Ele1MET_Type01 < 110 ){

         //  FillUserTH1D( "MTenu_Type01_50_110_gteOneBtaggedJet"      , MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  FillUserTH1D( "nJets_MTenu_Type01_50_110_gteOneBtaggedJet", nJetLooseEle_ptCut,  min_prescale * fakeRateEffective ) ;

         //  if ( nJetLooseEle_ptCut <= 3 ){ 
         //    FillUserTH1D(   "MTenu_Type01_50_110_Njet_lte3_gteOneBtaggedJet", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }

         //  if ( nJetLooseEle_ptCut <= 4 ){ 
         //    FillUserTH1D(   "MTenu_Type01_50_110_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }

         //  if ( nJetLooseEle_ptCut >= 4 ){ 
         //    FillUserTH1D(   "MTenu_Type01_50_110_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }

         //  if ( nJetLooseEle_ptCut >= 5 ){ 
         //    FillUserTH1D(   "MTenu_Type01_50_110_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }
         //}

         if ( MT_Ele1MET > 70 && MT_Ele1MET < 110 ){

           FillUserTH1D( "MTenu_70_110_gteOneBtaggedJet"      , MT_Ele1MET, min_prescale * fakeRateEffective );
           FillUserTH1D( "nJets_MTenu_70_110_gteOneBtaggedJet", nJetLooseEle_ptCut, min_prescale * fakeRateEffective );

           if ( nJetLooseEle_ptCut <= 3 ){
             FillUserTH1D(   "MTenu_70_110_Njet_lte3_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut <= 4 ){
             FillUserTH1D(   "MTenu_70_110_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 4 ){
             FillUserTH1D(   "MTenu_70_110_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_70_110_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }
         }

         //if ( MT_Ele1MET_Type01 > 70 && MT_Ele1MET_Type01 < 110 ){

         //  FillUserTH1D( "MTenu_Type01_70_110_gteOneBtaggedJet"      , MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  FillUserTH1D( "nJets_MTenu_Type01_70_110_gteOneBtaggedJet", nJetLooseEle_ptCut,  min_prescale * fakeRateEffective ) ;

         //  if ( nJetLooseEle_ptCut <= 3 ){ 
         //    FillUserTH1D(   "MTenu_Type01_70_110_Njet_lte3_gteOneBtaggedJet", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }

         //  if ( nJetLooseEle_ptCut <= 4 ){ 
         //    FillUserTH1D(   "MTenu_Type01_70_110_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }

         //  if ( nJetLooseEle_ptCut >= 4 ){ 
         //    FillUserTH1D(   "MTenu_Type01_70_110_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }

         //  if ( nJetLooseEle_ptCut >= 5 ){ 
         //    FillUserTH1D(   "MTenu_Type01_70_110_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET_Type01,  min_prescale * fakeRateEffective ) ;
         //  }
         //}
       }

       if ( fabs(LooseEle1_SCEta) <= eleEta_bar ) { 
         sprintf(plot_name,"Mej_Barrel_Presel");  FillUserTH1D( plot_name, Mej, min_prescale * fakeRateEffective);
         sprintf(plot_name,"Pt1stEle_Barrel_PAS");  FillUserTH1D( plot_name, LooseEle1_Pt, min_prescale * fakeRateEffective);
       }
       else if ( fabs(LooseEle1_SCEta) >= eleEta_end1_min && fabs(LooseEle1_SCEta) < eleEta_end1_max) {
         sprintf(plot_name,"Mej_Endcap1_Presel");  FillUserTH1D( plot_name, Mej, min_prescale * fakeRateEffective);
         sprintf(plot_name,"Pt1stEle_Endcap1_PAS");  FillUserTH1D( plot_name, LooseEle1_Pt, min_prescale * fakeRateEffective);
       }
       else if ( fabs(LooseEle1_SCEta) >= eleEta_end2_min && fabs(LooseEle1_SCEta) < eleEta_end2_max) {
         sprintf(plot_name,"Mej_Endcap2_Presel");  FillUserTH1D( plot_name, Mej, min_prescale * fakeRateEffective);
         sprintf(plot_name,"Pt1stEle_Endcap2_PAS");  FillUserTH1D( plot_name, LooseEle1_Pt, min_prescale * fakeRateEffective);
       }

       //-------------------------------------------------------------------------- 
       // Final selection plots
       //-------------------------------------------------------------------------- 

       if(!isOptimizationEnabled()) {
         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS      [i_lq_mass];
           bool pass    = passed_vector[i_lq_mass];
           if ( !pass ) continue;
           sprintf(plot_name, "ST_LQ%d"    , lq_mass ); FillUserTH1D( plot_name , sT_enujj   , min_prescale * fakeRateEffective ) ;
           sprintf(plot_name, "MTenu_LQ%d" , lq_mass ); FillUserTH1D( plot_name , MT_Ele1MET , min_prescale * fakeRateEffective ) ;
           sprintf(plot_name, "Mej_LQ%d"   , lq_mass ); FillUserTH1D( plot_name , Mej        , min_prescale * fakeRateEffective ) ;
           if ( fabs(LooseEle1_SCEta) <= eleEta_bar ) { 
             sprintf(plot_name,"Mej_Barrel_LQ%d", lq_mass);  FillUserTH1D( plot_name, Mej, min_prescale * fakeRateEffective);
           }
           else if ( fabs(LooseEle1_SCEta) >= eleEta_end1_min && fabs(LooseEle1_SCEta) < eleEta_end1_max) {
             sprintf(plot_name,"Mej_Endcap1_LQ%d", lq_mass);  FillUserTH1D( plot_name, Mej, min_prescale * fakeRateEffective);
           }
           else if ( fabs(LooseEle1_SCEta) >= eleEta_end2_min && fabs(LooseEle1_SCEta) < eleEta_end2_max) {
             sprintf(plot_name,"Mej_Endcap2_LQ%d", lq_mass);  FillUserTH1D( plot_name, Mej, min_prescale * fakeRateEffective);
           }
           if ( do_extra_finalSelection_plots ) {
             sprintf(plot_name, "Pt1stEle_LQ%d"	           , lq_mass ); FillUserTH1D( plot_name ,LooseEle1_Pt,          min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "PtHeep1stEle_LQ%d"	       , lq_mass ); FillUserTH1D( plot_name ,LooseEle1_PtHeep,      min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "Eta1stEle_LQ%d"           , lq_mass ); FillUserTH1D( plot_name ,LooseEle1_Eta,         min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "SCEta1stEle_LQ%d"         , lq_mass ); FillUserTH1D( plot_name ,LooseEle1_SCEta,       min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "SCEt1stEle_LQ%d"          , lq_mass ); FillUserTH1D( plot_name ,LooseEle1_SCEnergy/cosh(LooseEle1_SCEta),min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "Phi1stEle_LQ%d"           , lq_mass ); FillUserTH1D( plot_name ,LooseEle1_Phi,         min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "Pt1stJet_LQ%d"            , lq_mass ); FillUserTH1D( plot_name ,JetLooseEle1_Pt,       min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "Pt2ndJet_LQ%d"            , lq_mass ); FillUserTH1D( plot_name ,JetLooseEle2_Pt,       min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "Eta1stJet_LQ%d"           , lq_mass ); FillUserTH1D( plot_name ,JetLooseEle1_Eta,      min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "Eta2ndJet_LQ%d"           , lq_mass ); FillUserTH1D( plot_name ,JetLooseEle2_Eta,      min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "Phi1stJet_LQ%d"           , lq_mass ); FillUserTH1D( plot_name ,JetLooseEle1_Phi,      min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "Phi2ndJet_LQ%d"           , lq_mass ); FillUserTH1D( plot_name ,JetLooseEle2_Phi,      min_prescale * fakeRateEffective ); 
             // muon kinematics
             sprintf(plot_name, "Pt1stMuon_LQ%d"	         , lq_mass ); FillUserTH1D( plot_name ,Muon1_Pt,  min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "Eta1stMuon_LQ%d"          , lq_mass ); FillUserTH1D( plot_name ,Muon1_Eta, min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "Phi1stMuon_LQ%d"          , lq_mass ); FillUserTH1D( plot_name ,Muon1_Phi, min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "Pt2ndMuon_LQ%d"	         , lq_mass ); FillUserTH1D( plot_name ,Muon2_Pt,  min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "Eta2ndMuon_LQ%d"         , lq_mass );  FillUserTH1D( plot_name ,Muon2_Eta, min_prescale * fakeRateEffective ); 
             sprintf(plot_name, "Phi2ndMuon_LQ%d"         , lq_mass );  FillUserTH1D( plot_name ,Muon2_Phi, min_prescale * fakeRateEffective ); 
           }
         }
       }

     } // if passed_preselection
   } // end loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
