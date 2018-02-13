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
// for scale factors
#include "ElectronScaleFactors.C"

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
   
   char cut_name[100];

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
   CreateUserTH1D( "Pt1stEle_PAS"	      , 200 , 0       , 2000     ); 
   CreateUserTH1D( "PtHeep1stEle_Presel"	      , 200 , 0       , 2000     ); 
   CreateUserTH1D( "SCEt1stEle_Presel"	      , 200 , 0       , 2000     ); 
   CreateUserTH1D( "Eta1stEle_PAS"	      , 100 , -5      , 5	 ); 
   CreateUserTH1D( "SCEta1stEle_Presel"	      , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stEle_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Pt1stEle_Barrel_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "Pt1stEle_Endcap1_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "Pt1stEle_Endcap2_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "DeltaEtaEleTrk1stEle_Presel", 400, -0.5,   0.5 );
   // muon kinematics
   CreateUserTH1D( "Pt1stMuon_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta1stMuon_PAS"	      , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stMuon_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Pt2ndMuon_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta2ndMuon_PAS"	      , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi2ndMuon_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   //
   CreateUserTH1D( "Charge1stEle_PAS"	      , 2   , -1.0001 , 1.0001	 ); 
   CreateUserTH1D( "MET_PAS"                  , 600 , 0       , 3000	 ); 
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
   CreateUserTH1D( "Mee_allElectrons_Presel"                  , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mee_allElectrons_3EleEvents_Presel"                  , 200 , 0       , 2000	 ); 
   //
   CreateUserTH2D( "Mej_vs_EleSCEt" ,    200 ,  0, 2000, 200, 0, 2000);
   CreateUserTH2D( "Mej_vs_Jet1Pt" ,    200 ,  0, 2000, 200, 0, 2000);
   CreateUserTH2D( "Mej_vs_Jet2Pt" ,    200 ,  0, 2000, 200, 0, 2000);
   CreateUserTH2D( "Mej_vs_SelJetPt" ,    200 ,  0, 2000, 200, 0, 2000);

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

   CreateUserTH1D( "MTenu_50_110_noBtaggedJets", 240, 40, 160 );
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
   //
   CreateUserTH1D( "MTenu_50_110_Njet_gte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte3_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte4_noBtaggedJets", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte5_noBtaggedJets", 240, 40, 160 );

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

   CreateUserTH1D( "MTenu_50_110_gteOneBtaggedJet", 240, 40, 160 );
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
   //
   CreateUserTH1D( "MTenu_50_110_Njet_gte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte3_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_lte4_gteOneBtaggedJet", 240, 40, 160 );
   CreateUserTH1D( "MTenu_50_110_Njet_gte5_gteOneBtaggedJet", 240, 40, 160 );

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

   CreateUserTH1D( "Eta1stJet_PASand2Jet"  , 100 , -5 , 5 ); 
   CreateUserTH1D( "Eta1stJet_PASand3Jet"  , 100 , -5 , 5 ); 
   CreateUserTH1D( "Eta1stJet_PASand4Jet"  , 100 , -5 , 5 ); 
   CreateUserTH1D( "Eta1stJet_PASand5Jet"  , 100 , -5 , 5 ); 

   CreateUserTH1D( "Mej_Barrel_Presel"                  , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mej_Endcap1_Presel"                  , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mej_Endcap2_Presel"                  , 200 , 0       , 2000	 ); 
   //
   CreateUserTH1D("MejGte1500_Pt1stEle_Barrel_PAS"      ,200,0,2000);
   CreateUserTH1D("MejGte1500_PtJet_Barrel_PAS"         ,200,0,2000);
   CreateUserTH1D("MejGte1500_DREleJet_Barrel_PAS"      ,100,0,10);
   CreateUserTH1D("MejGte1500_DeltaPhiEleMET_Barrel_PAS",100,0,3.14159);
   CreateUserTH1D("MejGte1500_Pt1stEle_Endcap1_PAS"      ,200,0,2000);
   CreateUserTH1D("MejGte1500_PtJet_Endcap1_PAS"         ,200,0,2000);
   CreateUserTH1D("MejGte1500_DREleJet_Endcap1_PAS"      ,100,0,10);
   CreateUserTH1D("MejGte1500_DeltaPhiEleMET_Endcap1_PAS",100,0,3.14159);
   CreateUserTH1D("MejGte1500_Pt1stEle_Endcap2_PAS"      ,200,0,2000);
   CreateUserTH1D("MejGte1500_PtJet_Endcap2_PAS"         ,200,0,2000);
   CreateUserTH1D("MejGte1500_DREleJet_Endcap2_PAS"      ,100,0,10);
   CreateUserTH1D("MejGte1500_DeltaPhiEleMET_Endcap2_PAS",100,0,3.14159);
   CreateUserTH1D("MejGte1500_minDR_EleJet_PAS"         , 100 , 0       , 10       ); 
   CreateUserTH1D("MejGte1500_minDR_EleJet_Barrel_PAS"         , 100 , 0       , 10       ); 
   CreateUserTH1D("MejGte1500_minDR_EleJet_Endcap1_PAS"         , 100 , 0       , 10       ); 
   CreateUserTH1D("MejGte1500_minDR_EleJet_Endcap2_PAS"         , 100 , 0       , 10       ); 

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
   
   char plot_name[100];
  
   if (doFinalSelections ) { 
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

     if(jentry < 10 || jentry%10000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << "/" << nentries << std::endl;   

     //XXX SIC TEST: only run over 275376
     //if(static_cast<int>(run)!=275376) continue;
     //if(lumi!=1539) continue;
     //if(event!=2376331393) continue;
     //XXX end SIC TEST

     //// run ls event
     ////std::cout << "[Preselection] passing run/ls/event: " << static_cast<int>(run) << " " << static_cast<int>(ls) << " " << ((unsigned int)event) << std::endl;
     //std::cout << static_cast<int>(run) << " " << static_cast<int>(ls) << " " << ((unsigned int)event) << std::endl;
     ////std::cout << "\tMET: " << PFMET_Type1XY_Pt << std::endl;
     ////std::cout << "\tMETphi: " << PFMET_Type1XY_Phi << std::endl;
     //std::cout << "\tSCPt: " << LooseEle1_SCEt << std::endl;
     //std::cout << "\tSCEta: " << LooseEle1_Eta << std::endl;
     //std::cout << "\tFailHEEP: " << !LooseEle1_PassHEEPID << std::endl;
     //std::cout << "\tEcalDriven: " << LooseEle1_EcalDriven << std::endl;
     //std::cout << "\tSigmaIetaIeta: " << LooseEle1_Full5x5SigmaIEtaIEta << std::endl;
     //std::cout << "\tdxy: " << LooseEle1_LeadVtxDistXY << std::endl;
     //float hoe2 = LooseEle1_HoE * LooseEle1_PtHeep/LooseEle1_SCEt;
     //std::cout << "\tHoE: " << hoe2 << std::endl;
     //std::cout << "\tMissingHits: " << LooseEle1_MissingHits << std::endl;


     //--------------------------------------------------------------------------
     // Reset the cuts
     //--------------------------------------------------------------------------

     resetCuts();

     //--------------------------------------------------------------------------
     // Check good run list
     //--------------------------------------------------------------------------
     
     int passedJSON = passJSON ( run, ls , isData ) ;

     //--------------------------------------------------------------------------
     // Do pileup re-weighting for MC
     //--------------------------------------------------------------------------
     double pileup_weight = getPileupWeight ( nPileUpInt_True, isData ) ;
     
     //--------------------------------------------------------------------------
     // Get information about gen-level reweighting (should be for Sherpa only)
     //--------------------------------------------------------------------------
     double gen_weight = Weight;
     if ( isData ) gen_weight = 1.0;

     // TopPt reweight
     // only valid for powheg
     std::string current_file_name ( fChain->GetCurrentFile()->GetName());
     if(current_file_name.find("TT_") != std::string::npos) {
       gen_weight*=TopPtWeight;
     }

     // Electron scale factors for MC only
     if(!isData)
     {
       float recoSFEle1 = ElectronScaleFactors2016::LookupRecoSF(LooseEle1_SCEta);
       float heepSFEle1 = ElectronScaleFactors2016::LookupHeepSF(LooseEle1_SCEta);
       float totalScaleFactor = recoSFEle1*heepSFEle1;
       gen_weight*=totalScaleFactor;
     }

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
     // for MC, need to set prescale to use to gen/pileup weight
     // we can still check trigger and trigger matching as done above
     // should apply scale factors to cover any data/MC trigger efficiency differences, but this is probably better than ignoring the trigger decision in the MC
     if(!isData) {
       min_prescale = gen_weight*pileup_weight;
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
     // PtHeep has bugged energy correction (for loose electrons); don't use it
     float fakeRate1 = qcdFakeRate.GetFakeRate(LooseEle1_SCEta,LooseEle1_SCEnergy/cosh(LooseEle1_SCEta));
     
     //--------------------------------------------------------------------------
     // Finally have the effective fake rate
     //--------------------------------------------------------------------------

     //double fakeRateEffective  = fakeRate1;
     double fakeRateEffective  = fakeRate1/(1-fakeRate1); // require loose electron to fail HEEP ID
     //double fakeRateEffective = 1.0; // turn off fake rate
     double eFakeRateEffective = 0.0; //FIXME eFakeRate1;
     
     //--------------------------------------------------------------------------
     // Calculate some variables:
     //--------------------------------------------------------------------------
     
     TLorentzVector loose_ele1, loose_ele2, jet1, jet2, met;
     //loose_ele1.SetPtEtaPhiM ( LooseEle1_Pt , LooseEle1_Eta , LooseEle1_Phi , 0.0 );
     //loose_ele2.SetPtEtaPhiM ( LooseEle2_Pt , LooseEle2_Eta , LooseEle2_Phi , 0.0 );
     // need to use uncorrected Pt
     //loose_ele1.SetPtEtaPhiM ( LooseEle1_SCEnergy/cosh(LooseEle1_SCEta) , LooseEle1_Eta , LooseEle1_Phi , 0.0 );
     //loose_ele2.SetPtEtaPhiM ( LooseEle2_SCEnergy/cosh(LooseEle2_SCEta) , LooseEle2_Eta , LooseEle2_Phi , 0.0 );
     loose_ele1.SetPtEtaPhiM ( LooseEle1_SCEnergy/cosh(LooseEle1_SCEta) , LooseEle1_Eta , LooseEle1_Phi , 0.0 );
     loose_ele2.SetPtEtaPhiM ( LooseEle2_SCEnergy/cosh(LooseEle2_SCEta) , LooseEle2_Eta , LooseEle2_Phi , 0.0 );
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
     //MT_Ele1MET = sqrt(2 * LooseEle1_Pt * PFMET_Type1XY_Pt  * (1 - cos(loose_ele1.DeltaPhi ( met) ) ) );
     MT_Ele1MET = sqrt(2 * LooseEle1_SCEnergy/cosh(LooseEle1_SCEta) * PFMET_Type1XY_Pt  * (1 - cos(loose_ele1.DeltaPhi ( met) ) ) );

     mDPhi_METEle1= fabs(loose_ele1.DeltaPhi ( met ));
     mDPhi_METJet1= fabs(jet1.DeltaPhi ( met ));
     mDPhi_METJet2= fabs(jet2.DeltaPhi ( met ));

     //sT_enujj = LooseEle1_Pt + PFMET_Type1XY_Pt + JetLooseEle1_Pt + JetLooseEle2_Pt ;
     sT_enujj = LooseEle1_SCEnergy/cosh(LooseEle1_SCEta) + PFMET_Type1XY_Pt + JetLooseEle1_Pt + JetLooseEle2_Pt ;
     
     DR_Ele1Jet1 = loose_ele1.DeltaR ( jet1 ) ;
     DR_Ele1Jet2 = loose_ele1.DeltaR ( jet2 ) ;

     //XXX SIC REMOVE BAD ENERGY CORRECTIONS TO HoE
     float hoe = LooseEle1_HoE * LooseEle1_PtHeep/LooseEle1_SCEt;
     float hoeThresh = ele1_isBarrel ? 0.15 : 0.10;
     int passHoE = 0;
     if(hoe < hoeThresh)
       passHoE = 1;
     //std::cout << "correct HoE from: " << LooseEle1_HoE << " to " << hoe << "; corr=ptheep/scEt=" << LooseEle1_PtHeep << "/" << LooseEle1_SCEt << " = " << LooseEle1_PtHeep/LooseEle1_SCEt << endl;

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
     fillVariableWithValue(   "Ele1_SCEt"                , LooseEle1_SCEnergy/cosh(LooseEle1_SCEta)      , min_prescale * fakeRateEffective );
     fillVariableWithValue(   "Ele1_Eta"                 , LooseEle1_Eta         , min_prescale * fakeRateEffective );
     fillVariableWithValue(   "Ele1_IsBarrel"            , ele1_isBarrel         , min_prescale * fakeRateEffective );
     fillVariableWithValue(   "MTenu"                    , MT_Ele1MET            , min_prescale * fakeRateEffective );
     fillVariableWithValue(   "Ele1_PassHEEPID"          , LooseEle1_PassHEEPID  , min_prescale * fakeRateEffective );
     fillVariableWithValue(   "Ele1_PassHoE"             , passHoE  , min_prescale * fakeRateEffective );
     fillVariableWithValue(   "AbsDeltaEtaEleTrk"        , fabs(LooseEle1_Eta-LooseEle1_TrkEta),min_prescale * fakeRateEffective  );
									           
     // MET variables	                                      		           
     fillVariableWithValue(   "MET"                      , PFMET_Type1XY_Pt                  , min_prescale * fakeRateEffective );
     fillVariableWithValue(   "Pt_EMET"                  , Pt_Ele1MET             , min_prescale * fakeRateEffective );
     fillVariableWithValue(   "mDeltaPhiMETEle"          , mDPhi_METEle1           , min_prescale * fakeRateEffective );
     									           
     // 1st JET variables                                     		           
     fillVariableWithValue(   "nJet"                     , nJetLooseEle_ptCut      , min_prescale * fakeRateEffective );
			
     double MT_Jet1MET, MT_Jet2MET, MT_Ele1Jet1, MT_Ele1Jet2;//, MT_Ele1MET_Type01;
     //double mDPhi_METType01_Ele1, mDPhi_METType01_Jet1, mDPhi_METType01_Jet2;

     //// alternate METs
     //if ( nLooseEle_store > 0 ) {
     //  TVector2 v_ele;
     //  TVector2 v_MET_Type01;
     //  v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
     //  v_ele.SetMagPhi( LooseEle1_Pt, LooseEle1_Phi );
     //  mDPhi_METType01_Ele1 = fabs(v_MET_Type01.DeltaPhi ( v_ele ));
     //  float deltaphi = v_MET_Type01.DeltaPhi(v_ele);
     //  MT_Ele1MET_Type01 = sqrt ( 2 * LooseEle1_Pt * PFMET_Type01_Pt * ( 1 - cos ( deltaphi ) ) );
     //}
				           
     // 1st JET variables                                     		           
     if ( nJetLooseEle_store > 0 ) { 						           
       fillVariableWithValue( "Jet1_Pt"                  , JetLooseEle1_Pt         , min_prescale * fakeRateEffective );
       fillVariableWithValue( "Jet1_Eta"                 , JetLooseEle1_Eta        , min_prescale * fakeRateEffective );
       fillVariableWithValue( "mDeltaPhiMET1stJet"       , mDPhi_METJet1           , min_prescale * fakeRateEffective );

       TVector2 v_MET;
       TVector2 v_jet;
       //TVector2 v_MET_Type01;
       //v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       v_MET.SetMagPhi( PFMET_Type1XY_Pt , PFMET_Type1XY_Phi  );
       v_jet.SetMagPhi( JetLooseEle1_Pt, JetLooseEle1_Phi );
       //mDPhi_METType01_Jet1 = fabs(v_MET_Type01.DeltaPhi ( v_jet ));
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
       //TVector2 v_MET_Type01;
       //v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       v_MET.SetMagPhi( PFMET_Type1XY_Pt , PFMET_Type1XY_Phi  );
       v_jet.SetMagPhi( JetLooseEle2_Pt, JetLooseEle2_Phi );
       //mDPhi_METType01_Jet2 = fabs(v_MET_Type01.DeltaPhi ( v_jet ));
       float deltaphi = v_MET.DeltaPhi(v_jet);
       MT_Jet2MET = sqrt ( 2 * JetLooseEle2_Pt * PFMET_Type1XY_Pt * ( 1 - cos ( deltaphi ) ) );
     }
     else
     {
       fillVariableWithValue( "Jet2_Pt"                  , -1.0        , min_prescale * fakeRateEffective );
       fillVariableWithValue( "Jet2_Eta"                 , -1.0        , min_prescale * fakeRateEffective );
       fillVariableWithValue( "ST"                       , sT_enujj                , min_prescale * fakeRateEffective );
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
       //TVector2 v_MET_Type01;
       //v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       //v_ele .SetMagPhi ( LooseEle1_Pt, LooseEle1_Phi );
       // need to use uncorrected Pt
       v_ele .SetMagPhi ( LooseEle1_SCEnergy/cosh(LooseEle1_SCEta), LooseEle1_Phi );
       v_jet1.SetMagPhi ( JetLooseEle1_Pt, JetLooseEle1_Phi );
       float deltaphi = v_ele.DeltaPhi ( v_jet1 );
       //MT_Ele1Jet1 = sqrt ( 2 * JetLooseEle1_Pt * LooseEle1_Pt * ( 1 - cos ( deltaphi ) ) );
       MT_Ele1Jet1 = sqrt ( 2 * JetLooseEle1_Pt * LooseEle1_SCEnergy/cosh(LooseEle1_SCEta) * ( 1 - cos ( deltaphi ) ) );

     }

     // 1 electron, 2 jet variables 
     if ( nLooseEle_store > 0 && nJetLooseEle_store > 1 ) { 
       fillVariableWithValue ( "DR_Ele1Jet2"             , DR_Ele1Jet2             , min_prescale * fakeRateEffective );
       
       TVector2 v_ele;
       TVector2 v_jet2;
       //TVector2 v_MET_Type01;
       //v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
       //v_ele .SetMagPhi ( LooseEle1_Pt, LooseEle1_Phi );
       // need to use uncorrected Pt
       v_ele .SetMagPhi ( LooseEle1_SCEnergy/cosh(LooseEle1_SCEta), LooseEle1_Phi );
       v_jet2.SetMagPhi ( JetLooseEle2_Pt, JetLooseEle2_Phi );
       float deltaphi = v_ele.DeltaPhi ( v_jet2 );
       //MT_Ele1Jet2 = sqrt ( 2 * JetLooseEle2_Pt * LooseEle1_Pt * ( 1 - cos ( deltaphi ) ) );
       MT_Ele1Jet2 = sqrt ( 2 * JetLooseEle2_Pt * LooseEle1_SCEnergy/cosh(LooseEle1_SCEta) * ( 1 - cos ( deltaphi ) ) );

     }
     else
       fillVariableWithValue ( "DR_Ele1Jet2"             , -1.0             , min_prescale * fakeRateEffective );
     
     double MT_JetMET;
     double Mej;
     bool mejSelectedJet1 = true;
     
     if ( fabs ( MT_Jet1MET - MT_Ele1Jet2 ) < fabs( MT_Jet2MET - MT_Ele1Jet1 )){
       MT_JetMET = MT_Jet1MET;
       Mej = M_e1j2;
       mejSelectedJet1 = false;
     } else { 
       MT_JetMET = MT_Jet2MET;
       Mej = M_e1j1;
     }	 

     double M_e1e2 = 0.0;
     if ( nLooseEle_store > 1) {
       TLorentzVector v_ele1, v_ele2, v_sum;
       v_ele1.SetPtEtaPhiM ( LooseEle1_SCEnergy/cosh(LooseEle1_SCEta), LooseEle1_Eta, LooseEle1_Phi, 0.0 );
       v_ele2.SetPtEtaPhiM ( LooseEle2_SCEnergy/cosh(LooseEle2_SCEta), LooseEle2_Eta, LooseEle2_Phi, 0.0 );
       v_sum = v_ele1+v_ele2;
       M_e1e2 = v_sum.M();
     }
     double M_e1e3 = 0.0;
     double M_e2e3 = 0.0;
     if ( nLooseEle_store > 2) {
       TLorentzVector v_ele1, v_ele2, v_ele3, v_sum;
       v_ele1.SetPtEtaPhiM ( LooseEle1_SCEnergy/cosh(LooseEle1_SCEta), LooseEle1_Eta, LooseEle1_Phi, 0.0 );
       v_ele2.SetPtEtaPhiM ( LooseEle2_SCEnergy/cosh(LooseEle2_SCEta), LooseEle2_Eta, LooseEle2_Phi, 0.0 );
       v_ele3.SetPtEtaPhiM ( LooseEle3_SCEnergy/cosh(LooseEle3_SCEta), LooseEle3_Eta, LooseEle3_Phi, 0.0 );
       v_sum = v_ele1+v_ele3;
       M_e1e3 = v_sum.M();
       v_sum = v_ele2+v_ele3;
       M_e2e3 = v_sum.M();
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
     int nBJet_medium_ptCut_beyondLeadingTwo = 0;
     
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
     //
     if ( JetLooseEle3_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut_beyondLeadingTwo++;
     if ( JetLooseEle4_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut_beyondLeadingTwo++;
     if ( JetLooseEle5_btagCISV > btagCISV_medium_cut ) nBJet_medium_ptCut_beyondLeadingTwo++;
     
     if ( JetLooseEle1_btagCISV > btagCISV_tight_cut  ) nBJet_tight_ptCut++;
     if ( JetLooseEle2_btagCISV > btagCISV_tight_cut  ) nBJet_tight_ptCut++;
     if ( JetLooseEle3_btagCISV > btagCISV_tight_cut  ) nBJet_tight_ptCut++;
     if ( JetLooseEle4_btagCISV > btagCISV_tight_cut  ) nBJet_tight_ptCut++;
     if ( JetLooseEle5_btagCISV > btagCISV_tight_cut  ) nBJet_tight_ptCut++;
     	
     //--------------------------------------------------------------------------
     // Fill final selection cuts
     //--------------------------------------------------------------------------

     char cut_name[100];
     if(doFinalSelections) {
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

     if(doFinalSelections ) {
       passed_vector.clear();
       for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
         int lq_mass = LQ_MASS[i_lq_mass];
         sprintf(cut_name, "Mej_LQ%d", lq_mass );
         bool decision = bool ( passedAllPreviousCuts(cut_name) && passedCut (cut_name));
         passed_vector.push_back (decision);
       }
     }
     
     // for 3 electron event dists
     // passed preselection without nEle cut
     if(passedAllOtherCuts("nEle")) {
       if(nLooseEle_store > 1)
         FillUserTH1D("Mee_allElectrons_Presel", M_e1e2, min_prescale * fakeRateEffective);
       if(nLooseEle_store > 2) {
         FillUserTH1D("Mee_allElectrons_Presel", M_e1e3, min_prescale * fakeRateEffective);
         FillUserTH1D("Mee_allElectrons_Presel", M_e2e3, min_prescale * fakeRateEffective);
         FillUserTH1D("Mee_allElectrons_3EleEvents_Presel", M_e1e3, min_prescale * fakeRateEffective);
         FillUserTH1D("Mee_allElectrons_3EleEvents_Presel", M_e2e3, min_prescale * fakeRateEffective);
         FillUserTH1D("Mee_allElectrons_3EleEvents_Presel", M_e1e2, min_prescale * fakeRateEffective);
       }
     }

     if ( passed_preselection ) { 

       ////// run ls event
       //std::cout << "\t[Preselection] passing run/ls/event: " << static_cast<int>(run) << " " << static_cast<int>(ls) << " " << ((unsigned int)event) << std::endl;
       ////std::cout << static_cast<int>(run) << " " << static_cast<int>(ls) << " " << ((unsigned int)event) << std::endl;
       //////std::cout << "\tMET: " << PFMET_Type1XY_Pt << std::endl;
       //////std::cout << "\tMETphi: " << PFMET_Type1XY_Phi << std::endl;
       ////std::cout << "\tFailHEEP: " << !LooseEle1_PassHEEPID << std::endl;
       ////std::cout << "\tEcalDriven: " << LooseEle1_EcalDriven << std::endl;
       ////std::cout << "\tSigmaIetaIeta: " << LooseEle1_Full5x5SigmaIEtaIEta << std::endl;
       ////std::cout << "\tdxy: " << LooseEle1_LeadVtxDistXY << std::endl;
       ////std::cout << "\tHoE: " << hoe << std::endl;
       ////std::cout << "\tMissingHits: " << LooseEle1_MissingHits << std::endl;

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
         //ele1.SetPtEtaPhiM ( LooseEle1_Pt, LooseEle1_Eta, LooseEle1_Phi, 0.0 );
         // need to use uncorrected Pt
         ele1.SetPtEtaPhiM ( LooseEle1_SCEnergy/cosh(LooseEle1_Eta), LooseEle1_Eta, LooseEle1_Phi, 0.0 );
         jet3.SetPtEtaPhiM ( JetLooseEle3_Pt, JetLooseEle3_Eta, JetLooseEle3_Phi, 0.0 );
         DR_Ele1Jet3 = ele1.DeltaR ( jet3 ) ;
       }

       if ( DR_Ele1Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet1;
       if ( DR_Ele1Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet2;
       if ( nJetLooseEle_store > 2 ) {
         if ( DR_Ele1Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet3;
       }
       //if(min_DR_EleJet < 0.3) {
       //  cout << "WARNING: FOUND AN EVENT WITH MinDR(ele,jet) < 0.3; it's: " << min_DR_EleJet << endl;
       //  cout << "DR_Ele1Jet1 = " << DR_Ele1Jet1 << " DR_Ele1Jet2 = " << DR_Ele1Jet2;
       //  if ( nJetLooseEle_store > 2 )
       //    cout << "DR_Ele1Jet3 = " << DR_Ele1Jet3 << endl;
       //  else
       //    cout << endl;
       //  std::cout << static_cast<int>(run) << " " << static_cast<int>(ls) << " " << ((unsigned int)event) << std::endl;
       //}

       //double sT_enujj_Type01 = LooseEle1_Pt + PFMET_Type01_Pt + JetLooseEle1_Pt + JetLooseEle2_Pt;

       FillUserTH1D( "nElectron_PAS"              , nLooseEle_ptCut                                  , min_prescale * fakeRateEffective); 
       FillUserTH1D( "nMuon_PAS"                  , nMuon_ptCut                                      , min_prescale * fakeRateEffective); 
       FillUserTH1D( "Pt1stEle_PAS"	              , LooseEle1_Pt                                     , min_prescale * fakeRateEffective); 
       FillUserTH1D( "PtHeep1stEle_Presel"	          , LooseEle1_PtHeep                                 , min_prescale * fakeRateEffective); 
       FillUserTH1D( "SCEta1stEle_Presel"	          , LooseEle1_SCEta                                  , min_prescale * fakeRateEffective); 
       FillUserTH1D( "SCEt1stEle_Presel"	            , LooseEle1_SCEnergy/cosh(LooseEle1_SCEta)         , min_prescale * fakeRateEffective); 
       FillUserTH1D( "Eta1stEle_PAS"	            , LooseEle1_Eta                                    , min_prescale * fakeRateEffective);
       FillUserTH1D( "Phi1stEle_PAS"	            , LooseEle1_Phi                                    , min_prescale * fakeRateEffective);
       FillUserTH1D( "DeltaEtaEleTrk1stEle_Presel"       , fabs(LooseEle1_Eta-LooseEle1_TrkEta)  , min_prescale * fakeRateEffective);
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
       //FillUserTH1D( "minMETPt1stEle_PAS"         , TMath::Min ( LooseEle1_Pt, PFMET_Type1XY_Pt  )  , min_prescale * fakeRateEffective);
       // need to use uncorrected Pt
       FillUserTH1D( "minMETPt1stEle_PAS"         , TMath::Min ( LooseEle1_SCEnergy/cosh(LooseEle1_SCEta), PFMET_Type1XY_Pt  )  , min_prescale * fakeRateEffective);
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
       // need to use uncorrected Pt
       //FillUserTH1D( "sTlep_PAS"                  , LooseEle1_Pt + PFMET_Type1XY_Pt                 , min_prescale * fakeRateEffective);
       FillUserTH1D( "sTlep_PAS"                  , LooseEle1_SCEnergy/cosh(LooseEle1_SCEta) + PFMET_Type1XY_Pt                 , min_prescale * fakeRateEffective);
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
       FillUserTH2D("Mej_vs_EleSCEt",LooseEle1_SCEnergy/cosh(LooseEle1_SCEta),Mej, min_prescale * fakeRateEffective);
       FillUserTH2D("Mej_vs_Jet1Pt",JetLooseEle1_Pt,Mej, min_prescale * fakeRateEffective);
       FillUserTH2D("Mej_vs_Jet2Pt",JetLooseEle1_Pt,Mej, min_prescale * fakeRateEffective);
       if(mejSelectedJet1)
         FillUserTH2D("Mej_vs_SelJetPt",JetLooseEle1_Pt,Mej, min_prescale * fakeRateEffective);
       else
         FillUserTH2D("Mej_vs_SelJetPt",JetLooseEle2_Pt,Mej, min_prescale * fakeRateEffective);
       
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
       
       if ( MT_Ele1MET > 70 && MT_Ele1MET < 150 ){

         FillUserTH1D( "MTenu_70_150"      , MT_Ele1MET, min_prescale * fakeRateEffective );
         FillUserTH1D( "nJets_MTenu_70_150", nJetLooseEle_ptCut, min_prescale * fakeRateEffective );

         if ( nJetLooseEle_ptCut <= 3 ){
           FillUserTH1D(   "MTenu_70_150_Njet_lte3", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
         }

         if ( nJetLooseEle_ptCut <= 4 ){
           FillUserTH1D(   "MTenu_70_150_Njet_lte4", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
         }

         if ( nJetLooseEle_ptCut >= 4 ){
           FillUserTH1D(   "MTenu_70_150_Njet_gte4", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
         }

         if ( nJetLooseEle_ptCut >= 5 ){ 
           FillUserTH1D(   "MTenu_70_150_Njet_gte5", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
         }
       }

       if ( MT_Ele1MET > 110 && MT_Ele1MET < 190 ){

         FillUserTH1D( "MTenu_110_190"      , MT_Ele1MET, min_prescale * fakeRateEffective );
         FillUserTH1D( "nJets_MTenu_110_190", nJetLooseEle_ptCut, min_prescale * fakeRateEffective );

         if ( nJetLooseEle_ptCut <= 3 ){
           FillUserTH1D(   "MTenu_110_190_Njet_lte3", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
         }

         if ( nJetLooseEle_ptCut <= 4 ){
           FillUserTH1D(   "MTenu_110_190_Njet_lte4", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
         }

         if ( nJetLooseEle_ptCut >= 4 ){
           FillUserTH1D(   "MTenu_110_190_Njet_gte4", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
         }

         if ( nJetLooseEle_ptCut >= 5 ){ 
           FillUserTH1D(   "MTenu_110_190_Njet_gte5", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
         }
       }

       // no b-tagged jets
       if (nBJet_medium_ptCut==0) {

         if ( MT_Ele1MET > 50 && MT_Ele1MET < 110 ){

           FillUserTH1D( "MTenu_50_110_noBtaggedJets"      , MT_Ele1MET, min_prescale * fakeRateEffective );
           FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets", nJetLooseEle_ptCut, min_prescale * fakeRateEffective );
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


         if ( MT_Ele1MET > 70 && MT_Ele1MET < 150 ){

           FillUserTH1D( "MTenu_70_150_noBtaggedJets"      , MT_Ele1MET, min_prescale * fakeRateEffective );
           FillUserTH1D( "nJets_MTenu_70_150_noBtaggedJets", nJetLooseEle_ptCut, min_prescale * fakeRateEffective );

           if ( nJetLooseEle_ptCut <= 3 ){
             FillUserTH1D(   "MTenu_70_150_Njet_lte3_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut <= 4 ){
             FillUserTH1D(   "MTenu_70_150_Njet_lte4_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 4 ){
             FillUserTH1D(   "MTenu_70_150_Njet_gte4_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_gte5_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }
         }

         if ( MT_Ele1MET > 110 && MT_Ele1MET < 190 ){

           FillUserTH1D( "MTenu_110_190_noBtaggedJets"      , MT_Ele1MET, min_prescale * fakeRateEffective );
           FillUserTH1D( "nJets_MTenu_110_190_noBtaggedJets", nJetLooseEle_ptCut, min_prescale * fakeRateEffective );

           if ( nJetLooseEle_ptCut <= 3 ){
             FillUserTH1D(   "MTenu_110_190_Njet_lte3_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut <= 4 ){
             FillUserTH1D(   "MTenu_110_190_Njet_lte4_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 4 ){
             FillUserTH1D(   "MTenu_110_190_Njet_gte4_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_110_190_Njet_gte5_noBtaggedJets", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }
         }
       }

       // at least one b-tagged jet
       if (nBJet_medium_ptCut>=1) {

         if ( MT_Ele1MET > 50 && MT_Ele1MET < 110 ){

           FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet"      , MT_Ele1MET, min_prescale * fakeRateEffective );
           FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet", nJetLooseEle_ptCut, min_prescale * fakeRateEffective );
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

         if ( MT_Ele1MET > 70 && MT_Ele1MET < 150 ){

           FillUserTH1D( "MTenu_70_150_gteOneBtaggedJet"      , MT_Ele1MET, min_prescale * fakeRateEffective );
           FillUserTH1D( "nJets_MTenu_70_150_gteOneBtaggedJet", nJetLooseEle_ptCut, min_prescale * fakeRateEffective );

           if ( nJetLooseEle_ptCut <= 3 ){
             FillUserTH1D(   "MTenu_70_150_Njet_lte3_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut <= 4 ){
             FillUserTH1D(   "MTenu_70_150_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 4 ){
             FillUserTH1D(   "MTenu_70_150_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_70_150_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }
         }

         if ( MT_Ele1MET > 110 && MT_Ele1MET < 190 ){

           FillUserTH1D( "MTenu_110_190_gteOneBtaggedJet"      , MT_Ele1MET, min_prescale * fakeRateEffective );
           FillUserTH1D( "nJets_MTenu_110_190_gteOneBtaggedJet", nJetLooseEle_ptCut, min_prescale * fakeRateEffective );

           if ( nJetLooseEle_ptCut <= 3 ){
             FillUserTH1D(   "MTenu_110_190_Njet_lte3_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut <= 4 ){
             FillUserTH1D(   "MTenu_110_190_Njet_lte4_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 4 ){
             FillUserTH1D(   "MTenu_110_190_Njet_gte4_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }

           if ( nJetLooseEle_ptCut >= 5 ){ 
             FillUserTH1D(   "MTenu_110_190_Njet_gte5_gteOneBtaggedJet", MT_Ele1MET,  min_prescale * fakeRateEffective ) ;
           }
         }
       }

       if ( fabs(LooseEle1_SCEta) <= eleEta_bar ) { 
         sprintf(plot_name,"Mej_Barrel_Presel");  FillUserTH1D( plot_name, Mej, min_prescale * fakeRateEffective);
         sprintf(plot_name,"Pt1stEle_Barrel_PAS");  FillUserTH1D( plot_name, LooseEle1_SCEnergy/cosh(LooseEle1_SCEta), min_prescale * fakeRateEffective);
       }
       else if ( fabs(LooseEle1_SCEta) >= eleEta_end1_min && fabs(LooseEle1_SCEta) < eleEta_end1_max) {
         sprintf(plot_name,"Mej_Endcap1_Presel");  FillUserTH1D( plot_name, Mej, min_prescale * fakeRateEffective);
         sprintf(plot_name,"Pt1stEle_Endcap1_PAS");  FillUserTH1D( plot_name, LooseEle1_SCEnergy/cosh(LooseEle1_SCEta), min_prescale * fakeRateEffective);
       }
       else if ( fabs(LooseEle1_SCEta) >= eleEta_end2_min && fabs(LooseEle1_SCEta) < eleEta_end2_max) {
         sprintf(plot_name,"Mej_Endcap2_Presel");  FillUserTH1D( plot_name, Mej, min_prescale * fakeRateEffective);
         sprintf(plot_name,"Pt1stEle_Endcap2_PAS");  FillUserTH1D( plot_name, LooseEle1_SCEnergy/cosh(LooseEle1_SCEta), min_prescale * fakeRateEffective);
       }
       //if(MT_JetMET < 100) {
       //  //// run ls event
       //  std::cout << "\t[Preselection MT_JetMET < 100] passing run/ls/event: " << static_cast<int>(run) << " " << static_cast<int>(ls) << " " << ((unsigned int)event) << std::endl;
       //  //std::cout << static_cast<int>(run) << " " << static_cast<int>(ls) << " " << ((unsigned int)event) << std::endl;
       //  //std::cout << "\thltPhotonPt: " << LooseEle1_hltPhotonPt << std::endl;
       //  //std::cout << "\tmin_prescale: " << min_prescale << std::endl;
       //  //std::cout << "\tFailHEEP: " << !LooseEle1_PassHEEPID << std::endl;
       //  //std::cout << "\tEcalDriven: " << LooseEle1_EcalDriven << std::endl;
       //  //std::cout << "\tSigmaIetaIeta: " << LooseEle1_Full5x5SigmaIEtaIEta << std::endl;
       //  //std::cout << "\tdxy: " << LooseEle1_LeadVtxDistXY << std::endl;
       //  //std::cout << "\tHoE: " << hoe << std::endl;
       //  //std::cout << "\tMissingHits: " << LooseEle1_MissingHits << std::endl;
       //  std::cout << "\tMT_JetMET: " << MT_JetMET << std::endl;
       //  std::cout << "\tMT_Jet1MET: " << MT_Jet1MET << std::endl;
       //  std::cout << "\tMT_Jet2MET: " << MT_Jet2MET << std::endl;
       //  std::cout << "\tMT_Ele1Jet1: " << MT_Ele1Jet1 << std::endl;
       //  std::cout << "\tMT_Ele1Jet2: " << MT_Ele1Jet2 << std::endl;
       //  std::cout << "\tM_e1j1: " << M_e1j1 << std::endl;
       //  std::cout << "\tM_e1j2: " << M_e1j2 << std::endl;
       //  std::cout << "\tfabs(MT_Jet1MET - MT_Ele1Jet2) = " << fabs(MT_Jet1MET - MT_Ele1Jet2) << std::endl;
       //  std::cout << "\tfabs(MT_Jet2MET - MT_Ele1Jet1) = " << fabs(MT_Jet2MET - MT_Ele1Jet1) << std::endl;
       //  if ( fabs ( MT_Jet1MET - MT_Ele1Jet2 ) < fabs( MT_Jet2MET - MT_Ele1Jet1 ))
       //    std::cout << "\t\tfabs ( MT_Jet1MET - MT_Ele1Jet2 ) < fabs( MT_Jet2MET - MT_Ele1Jet1 ) --> M_e1j2 selected" << std::endl;
       //  else
       //    std::cout << "\t\tfabs ( MT_Jet2MET - MT_Ele1Jet1 ) < fabs( MT_Jet1MET - MT_Ele1Jet2 ) --> M_e1j1 selected" << std::endl;
       //  std::cout << "\tSelected M_ej: " << Mej << std::endl;
       //  std::cout << "\tE1: pt=" << LooseEle1_SCEt << ", [check Pt]=" << LooseEle1_SCEnergy/cosh(LooseEle1_SCEta) << ", eta=" << LooseEle1_Eta << ", SCEta=" << LooseEle1_SCEta << ", phi=" << LooseEle1_Phi << std::endl;
       //  std::cout << "\tJ1: pt=" << JetLooseEle1_Pt << ", eta=" << JetLooseEle1_Eta << ", phi=" << JetLooseEle1_Phi << std::endl;
       //  std::cout << "\tJ2: pt=" << JetLooseEle2_Pt << ", eta=" << JetLooseEle2_Eta << ", phi=" << JetLooseEle2_Phi << std::endl;
       //  std::cout << "\tMET: " << PFMET_Type1XY_Pt << std::endl;
       //  std::cout << "\tMETphi: " << PFMET_Type1XY_Phi << std::endl;
       //}
       //// high Mej plots
       //if(Mej >= 1500) {
       //  ////// run ls event
       //  //std::cout << "\t[Preselection Mej>1500] passing run/ls/event: " << static_cast<int>(run) << " " << static_cast<int>(ls) << " " << ((unsigned int)event) << std::endl;
       //  ////std::cout << static_cast<int>(run) << " " << static_cast<int>(ls) << " " << ((unsigned int)event) << std::endl;
       //  //std::cout << "\thltPhotonPt: " << LooseEle1_hltPhotonPt << std::endl;
       //  //std::cout << "\tmin_prescale: " << min_prescale << std::endl;
       //  //std::cout << "\tFailHEEP: " << !LooseEle1_PassHEEPID << std::endl;
       //  //std::cout << "\tEcalDriven: " << LooseEle1_EcalDriven << std::endl;
       //  //std::cout << "\tSigmaIetaIeta: " << LooseEle1_Full5x5SigmaIEtaIEta << std::endl;
       //  //std::cout << "\tdxy: " << LooseEle1_LeadVtxDistXY << std::endl;
       //  //std::cout << "\tHoE: " << hoe << std::endl;
       //  //std::cout << "\tMissingHits: " << LooseEle1_MissingHits << std::endl;
       //  //std::cout << "\tMT_Jet1MET: " << MT_Jet1MET << std::endl;
       //  //std::cout << "\tMT_Jet2MET: " << MT_Jet2MET << std::endl;
       //  //std::cout << "\tMT_Ele1Jet1: " << MT_Ele1Jet1 << std::endl;
       //  //std::cout << "\tMT_Ele1Jet2: " << MT_Ele1Jet2 << std::endl;
       //  //std::cout << "\tM_e1j1: " << M_e1j1 << std::endl;
       //  //std::cout << "\tM_e1j2: " << M_e1j2 << std::endl;
       //  //std::cout << "\tfabs(MT_Jet1MET - MT_Ele1Jet2) = " << fabs(MT_Jet1MET - MT_Ele1Jet2) << std::endl;
       //  //std::cout << "\tfabs(MT_Jet2MET - MT_Ele1Jet1) = " << fabs(MT_Jet2MET - MT_Ele1Jet1) << std::endl;
       //  //if ( fabs ( MT_Jet1MET - MT_Ele1Jet2 ) < fabs( MT_Jet2MET - MT_Ele1Jet1 ))
       //  //  std::cout << "\t\tfabs ( MT_Jet1MET - MT_Ele1Jet2 ) < fabs( MT_Jet2MET - MT_Ele1Jet1 ) --> M_e1j2 selected" << std::endl;
       //  //else
       //  //  std::cout << "\t\tfabs ( MT_Jet2MET - MT_Ele1Jet1 ) < fabs( MT_Jet1MET - MT_Ele1Jet2 ) --> M_e1j1 selected" << std::endl;
       //  //std::cout << "\tSelected M_ej: " << Mej << std::endl;
       //  //std::cout << "\tE1: pt=" << LooseEle1_SCEt << ", [check Pt]=" << LooseEle1_SCEnergy/cosh(LooseEle1_SCEta) << ", eta=" << LooseEle1_Eta << ", SCEta=" << LooseEle1_SCEta << ", phi=" << LooseEle1_Phi << std::endl;
       //  //std::cout << "\tJ1: pt=" << JetLooseEle1_Pt << ", eta=" << JetLooseEle1_Eta << ", phi=" << JetLooseEle1_Phi << std::endl;
       //  //std::cout << "\tJ2: pt=" << JetLooseEle2_Pt << ", eta=" << JetLooseEle2_Eta << ", phi=" << JetLooseEle2_Phi << std::endl;
       //  //std::cout << "\tMET: " << PFMET_Type1XY_Pt << std::endl;
       //  //std::cout << "\tMETphi: " << PFMET_Type1XY_Phi << std::endl;

       //  FillUserTH1D( "MejGte1500_minDR_EleJet_PAS"           , min_DR_EleJet                                    , min_prescale * fakeRateEffective);
       //  if ( fabs(LooseEle1_SCEta) <= eleEta_bar ) { 
       //    FillUserTH1D( "MejGte1500_minDR_EleJet_Barrel_PAS"           , min_DR_EleJet                                    , min_prescale * fakeRateEffective);
       //    sprintf(plot_name,"MejGte1500_Pt1stEle_Barrel_PAS");  FillUserTH1D( plot_name, LooseEle1_SCEnergy/cosh(LooseEle1_SCEta), min_prescale * fakeRateEffective);
       //    sprintf(plot_name,"MejGte1500_DeltaPhiEleMET_Barrel_PAS");  FillUserTH1D( plot_name, mDPhi_METEle1, min_prescale * fakeRateEffective);
       //    if(mejSelectedJet1) {
       //      sprintf(plot_name,"MejGte1500_PtJet_Barrel_PAS");  FillUserTH1D( plot_name, JetLooseEle1_Pt, min_prescale * fakeRateEffective);
       //      sprintf(plot_name,"MejGte1500_DREleJet_Barrel_PAS");  FillUserTH1D( plot_name, DR_Ele1Jet1, min_prescale * fakeRateEffective);
       //    }
       //    else {
       //      sprintf(plot_name,"MejGte1500_PtJet_Barrel_PAS");  FillUserTH1D( plot_name, JetLooseEle2_Pt, min_prescale * fakeRateEffective);
       //      sprintf(plot_name,"MejGte1500_DREleJet_Barrel_PAS");  FillUserTH1D( plot_name, DR_Ele1Jet2, min_prescale * fakeRateEffective);
       //    }
       //  }
       //  else if ( fabs(LooseEle1_SCEta) >= eleEta_end1_min && fabs(LooseEle1_SCEta) < eleEta_end1_max) {
       //    FillUserTH1D( "MejGte1500_minDR_EleJet_Endcap1_PAS"           , min_DR_EleJet                                    , min_prescale * fakeRateEffective);
       //    sprintf(plot_name,"MejGte1500_Pt1stEle_Endcap1_PAS");  FillUserTH1D( plot_name, LooseEle1_SCEnergy/cosh(LooseEle1_SCEta), min_prescale * fakeRateEffective);
       //    sprintf(plot_name,"MejGte1500_DeltaPhiEleMET_Endcap1_PAS");  FillUserTH1D( plot_name, mDPhi_METEle1, min_prescale * fakeRateEffective);
       //    if(mejSelectedJet1) {
       //      sprintf(plot_name,"MejGte1500_PtJet_Endcap1_PAS");  FillUserTH1D( plot_name, JetLooseEle1_Pt, min_prescale * fakeRateEffective);
       //      sprintf(plot_name,"MejGte1500_DREleJet_Endcap1_PAS");  FillUserTH1D( plot_name, DR_Ele1Jet1, min_prescale * fakeRateEffective);
       //    }
       //    else {
       //      sprintf(plot_name,"MejGte1500_PtJet_Endcap1_PAS");  FillUserTH1D( plot_name, JetLooseEle2_Pt, min_prescale * fakeRateEffective);
       //      sprintf(plot_name,"MejGte1500_DREleJet_Endcap1_PAS");  FillUserTH1D( plot_name, DR_Ele1Jet2, min_prescale * fakeRateEffective);
       //    }
       //  }
       //  else if ( fabs(LooseEle1_SCEta) >= eleEta_end2_min && fabs(LooseEle1_SCEta) < eleEta_end2_max) {
       //    FillUserTH1D( "MejGte1500_minDR_EleJet_Endcap2_PAS"           , min_DR_EleJet                                    , min_prescale * fakeRateEffective);
       //    sprintf(plot_name,"MejGte1500_Pt1stEle_Endcap2_PAS");  FillUserTH1D( plot_name, LooseEle1_SCEnergy/cosh(LooseEle1_SCEta), min_prescale * fakeRateEffective);
       //    sprintf(plot_name,"MejGte1500_DeltaPhiEleMET_Endcap2_PAS");  FillUserTH1D( plot_name, mDPhi_METEle1, min_prescale * fakeRateEffective);
       //    if(mejSelectedJet1) {
       //      sprintf(plot_name,"MejGte1500_PtJet_Endcap2_PAS");  FillUserTH1D( plot_name, JetLooseEle1_Pt, min_prescale * fakeRateEffective);
       //      sprintf(plot_name,"MejGte1500_DREleJet_Endcap2_PAS");  FillUserTH1D( plot_name, DR_Ele1Jet1, min_prescale * fakeRateEffective);
       //    }
       //    else {
       //      sprintf(plot_name,"MejGte1500_PtJet_Endcap2_PAS");  FillUserTH1D( plot_name, JetLooseEle2_Pt, min_prescale * fakeRateEffective);
       //      sprintf(plot_name,"MejGte1500_DREleJet_Endcap2_PAS");  FillUserTH1D( plot_name, DR_Ele1Jet2, min_prescale * fakeRateEffective);
       //    }
       //  }
       //}

       //-------------------------------------------------------------------------- 
       // dijet mass control regions
       //-------------------------------------------------------------------------- 
       if(MT_Ele1MET > 50 && MT_Ele1MET < 110 ) {
         if(M_j1j2 > 50 && M_j1j2 < 110) {
           FillUserTH1D( "MTenu_50_110_Mjj50to110", MT_Ele1MET, min_prescale * fakeRateEffective);
           FillUserTH1D( "nJets_MTenu_50_110_Mjj50to110", nJetLooseEle_ptCut, min_prescale * fakeRateEffective);
           if(nBJet_medium_ptCut_beyondLeadingTwo>=1) {
             FillUserTH1D( "MTenu_50_110_Mjj50to110_addBtagJet", MT_Ele1MET, min_prescale * fakeRateEffective);
             FillUserTH1D( "nJets_MTenu_50_110_Mjj50to110_addBtagJet", nJetLooseEle_ptCut, min_prescale * fakeRateEffective);
           }
           if(nBJet_medium_ptCut_beyondLeadingTwo<1) {
             FillUserTH1D( "MTenu_50_110_Mjj50to110_noAddBtagJets", MT_Ele1MET, min_prescale * fakeRateEffective);
             FillUserTH1D( "nJets_MTenu_50_110_Mjj50to110_noAddBtagJets", nJetLooseEle_ptCut, min_prescale * fakeRateEffective);
           }
         }
         else if(M_j1j2 > 110) {
           FillUserTH1D( "MTenu_50_110_MjjGte110", MT_Ele1MET, min_prescale * fakeRateEffective);
           FillUserTH1D( "nJets_MTenu_50_110_MjjGte110", nJetLooseEle_ptCut, min_prescale * fakeRateEffective);
         }
       }

       //-------------------------------------------------------------------------- 
       // Final selection plots
       //-------------------------------------------------------------------------- 

       if(doFinalSelections) {
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

         // for scale factor at "final selection" studies
         if(nBJet_medium_ptCut==0) {
           if ( passedCut("ST_LQ300") && passedCut("Mej_LQ300") && passedCut("MET_LQ300") ){
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ300"       , MT_Ele1MET, min_prescale * fakeRateEffective);
             FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ300" , nJetLooseEle_ptCut, min_prescale * fakeRateEffective);
           }
           if ( passedCut("ST_LQ400") && passedCut("Mej_LQ400") && passedCut("MET_LQ400") ){
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ400"       , MT_Ele1MET, min_prescale * fakeRateEffective);
             FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ400" , nJetLooseEle_ptCut, min_prescale * fakeRateEffective);
           }
           if ( passedCut("ST_LQ500") && passedCut("Mej_LQ500") && passedCut("MET_LQ500") ){
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ500"       , MT_Ele1MET, min_prescale * fakeRateEffective);
             FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ500" , nJetLooseEle_ptCut, min_prescale * fakeRateEffective);
           }
           if ( passedCut("ST_LQ600") && passedCut("Mej_LQ600") && passedCut("MET_LQ600") ){
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ600"       , MT_Ele1MET, min_prescale * fakeRateEffective);
             FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ600" , nJetLooseEle_ptCut, min_prescale * fakeRateEffective);
           }
           if ( passedCut("ST_LQ700") && passedCut("Mej_LQ700") && passedCut("MET_LQ700") ){
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ700"       , MT_Ele1MET, min_prescale * fakeRateEffective);
             FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ700" , nJetLooseEle_ptCut, min_prescale * fakeRateEffective);
           }
           if ( passedCut("ST_LQ800") && passedCut("Mej_LQ800") && passedCut("MET_LQ800") ){
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ800"          , MT_Ele1MET, min_prescale * fakeRateEffective);
             FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ800"    , nJetLooseEle_ptCut, min_prescale * fakeRateEffective);
           }
           if ( passedCut("ST_LQ900") && passedCut("Mej_LQ900") && passedCut("MET_LQ900") ){
             FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ900"          , MT_Ele1MET, min_prescale * fakeRateEffective);
             FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ900"    , nJetLooseEle_ptCut, min_prescale * fakeRateEffective);
           }
           if ( passedCut("ST_LQ1000") && passedCut("Mej_LQ1000") && passedCut("MET_LQ1000") ){
           }
           FillUserTH1D( "MTenu_50_110_noBtaggedJets_LQ1000"       ,MT_Ele1MET, min_prescale * fakeRateEffective); 
           FillUserTH1D( "nJets_MTenu_50_110_noBtaggedJets_LQ1000" ,nJetLooseEle_ptCut, min_prescale * fakeRateEffective);
         }
         if(nBJet_medium_ptCut>=1) {
           if ( passedCut("ST_LQ300") && passedCut("Mej_LQ300") && passedCut("MET_LQ300") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ300"       , MT_Ele1MET,  min_prescale * fakeRateEffective );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ300" , nJetLooseEle_ptCut,  min_prescale * fakeRateEffective );
           }
           if ( passedCut("ST_LQ400") && passedCut("Mej_LQ400") && passedCut("MET_LQ400") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ400"       , MT_Ele1MET,  min_prescale * fakeRateEffective );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ400" , nJetLooseEle_ptCut,  min_prescale * fakeRateEffective );
           }
           if ( passedCut("ST_LQ500") && passedCut("Mej_LQ500") && passedCut("MET_LQ500") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ500"       , MT_Ele1MET,  min_prescale * fakeRateEffective );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ500" , nJetLooseEle_ptCut,  min_prescale * fakeRateEffective );
           }
           if ( passedCut("ST_LQ600") && passedCut("Mej_LQ600") && passedCut("MET_LQ600") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ600"       , MT_Ele1MET,  min_prescale * fakeRateEffective );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ600" , nJetLooseEle_ptCut,  min_prescale * fakeRateEffective );
           }
           if ( passedCut("ST_LQ700") && passedCut("Mej_LQ700") && passedCut("MET_LQ700") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ700"       , MT_Ele1MET,  min_prescale * fakeRateEffective );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ700" , nJetLooseEle_ptCut,  min_prescale * fakeRateEffective );
           }
           if ( passedCut("ST_LQ800") && passedCut("Mej_LQ800") && passedCut("MET_LQ800") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ800"       , MT_Ele1MET,  min_prescale * fakeRateEffective );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ800" , nJetLooseEle_ptCut,  min_prescale * fakeRateEffective );
           }
           if ( passedCut("ST_LQ900") && passedCut("Mej_LQ900") && passedCut("MET_LQ900") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ900"       , MT_Ele1MET,  min_prescale * fakeRateEffective );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ900" , nJetLooseEle_ptCut,  min_prescale * fakeRateEffective );
           }
           if ( passedCut("ST_LQ1000") && passedCut("Mej_LQ1000") && passedCut("MET_LQ1000") ){
             FillUserTH1D( "MTenu_50_110_gteOneBtaggedJet_LQ1000"       ,MT_Ele1MET,  min_prescale * fakeRateEffective );
             FillUserTH1D( "nJets_MTenu_50_110_gteOneBtaggedJet_LQ1000" ,nJetLooseEle_ptCut,  min_prescale * fakeRateEffective );
           }
         }

       } // end do final selections
     } // end passed_preselection
   } // end loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
