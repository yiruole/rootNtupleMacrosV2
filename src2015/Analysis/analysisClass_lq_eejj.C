#define analysisClass_cxx
#define USE_SINGLE_ELE_REDUCED_NTUPLE
#include "analysisClass.h"
#include <assert.h>
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TProfile.h>
// for trigger turn-on
#include "Ele27WPLooseTrigTurnOn.C"
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
   
   char cut_name[100];

   //--------------------------------------------------------------------------
   // Decide which plots to save (default is to save everything)
   //--------------------------------------------------------------------------
   
   fillSkim                         ( true  ) ;
   fillAllPreviousCuts              ( true  ) ;
   fillAllOtherCuts                 ( true  ) ;
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;

   bool do_roi_plots = false;
   //--------------------------------------------------------------------------
   // Any extra features
   //--------------------------------------------------------------------------
   
   TProfile * profile_run_vs_nvtx_HLT = new TProfile("run_vs_nvtx_HLT", "", 20000 , 160300  , 180300 );
   TProfile * profile_run_vs_nvtx_PAS = new TProfile("run_vs_nvtx_PAS", "", 20000 , 160300  , 180300 );
   
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

   //--------------------------------------------------------------------------
   // 2016 trigger efficiency
   //--------------------------------------------------------------------------
   std::string trigEffFileName = getPreCutString1("TriggerEfficiencyFileName");
   //std::string graphName = getPreCutString1("TriggerEfficiencyGraphName");
   // this needs to match the path we use for data below
   //TriggerEfficiency triggerEfficiency(trigEffFileName, "hEff_Ele27OR115");
   TriggerEfficiency triggerEfficiency(trigEffFileName, "hEff_Ele27OR115OR175");

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------
   
   CreateUserTH1D( "sTfrac_Jet1_PAS"                 ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Jet2_PAS"                 ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele1_PAS"                 ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele2_PAS"                 ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Jet_PAS"                  ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele_PAS"                  ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Jet1_PASandMee100"        ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Jet2_PASandMee100"        ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele1_PASandMee100"        ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele2_PASandMee100"        ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Jet_PASandMee100"         ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele_PASandMee100"         ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "Mej_selected_avg_PASandMee100"   ,    200   , 0       , 2000     );
   CreateUserTH1D( "EleChargeSum_PAS"                ,    3     , -2.5    , 2.5      );
   CreateUserTH1D( "EleChargeSum_PASandMee100"       ,    3     , -2.5    , 2.5      );
   //CreateUserTH1D( "ProcessID"                       ,    21    , -0.5    , 20.5     );
   //CreateUserTH1D( "ProcessID_PAS"                   ,    21    , -0.5    , 20.5     );
   //CreateUserTH1D( "ProcessID_ZWindow"               ,    21    , -0.5    , 20.5     );
   CreateUserTH1D( "nElectron_PAS"                   ,    5     , -0.5    , 4.5      );
   CreateUserTH1D( "nMuon_PAS"                       ,    5     , -0.5    , 4.5      );
   CreateUserTH1D( "nJet_PAS"                        ,    10    , -0.5    , 9.5      );
   CreateUserTH1D( "nJet_PASandMee100"               ,    10    , -0.5    , 9.5      );
   CreateUserTH1D( "Pt1stEle_PAS"	             , 	100    , 0       , 1000     ); 
   CreateUserTH1D( "PtHeep1stEle_PAS"	             , 	100    , 0       , 1000     ); 
   CreateUserTH1D( "Pt1stEle_PASandMee100"           , 	100    , 0       , 1000     ); 
   CreateUserTH1D( "Eta1stEle_PAS"	             , 	100    , -5      , 5	  ); 
   CreateUserTH1D( "SCEta1stEle_PAS"	             , 	100    , -5      , 5	  ); 
   CreateUserTH1D( "DeltaEtaEleTrk1stEle_Presel", 400, -0.5,   0.5 );
   CreateUserTH1D( "Phi1stEle_PAS"	             , 	60     , -3.1416 , +3.1416  ); 
   CreateUserTH1D( "Pt2ndEle_PAS"	             , 	300    , 0       , 3000     ); 
   CreateUserTH1D( "PtHeep2ndEle_PAS"	             , 	300    , 0       , 3000     ); 
   CreateUserTH1D( "Pt2ndEle_PASandMee100"           , 	300    , 0       , 3000     ); 
   CreateUserTH1D( "Eta2ndEle_PAS"	             , 	100    , -5      , 5	  ); 
   CreateUserTH1D( "SCEta2ndEle_PAS"	             , 	100    , -5      , 5	  ); 
   CreateUserTH1D( "DeltaEtaEleTrk2ndEle_Presel", 400, -0.5,   0.5 );
   CreateUserTH1D( "Phi2ndEle_PAS"	             , 	60     , -3.1416 , +3.1416  ); 
   CreateUserTH1D( "Charge1stEle_PAS"	             , 	2      , -1.0001 , 1.0001	  ); 
   CreateUserTH1D( "Charge2ndEle_PAS"	             , 	2      , -1.0001 , 1.0001	  ); 
   CreateUserTH1D( "MET_PAS"                         ,    200   , 0       , 1000	  ); 
   CreateUserTH1D( "METPhi_PAS"		             , 	60     , -3.1416 , +3.1416  ); 
   CreateUserTH1D( "Pt1stJet_PAS"                    ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Pt2ndJet_PAS"                    ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Pt1stJet_PASandMee100"           ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Pt2ndJet_PASandMee100"           ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Eta1stJet_PAS"                   ,    100   , -5      , 5	  ); 
   CreateUserTH1D( "Eta2ndJet_PAS"                   ,    100   , -5      , 5	  ); 
   CreateUserTH1D( "Phi1stJet_PAS"	             , 	 60    , -3.1416 , +3.1416  ); 
   CreateUserTH1D( "Phi2ndJet_PAS"	             , 	 60    , -3.1416 , +3.1416  ); 
   CreateUserTH1D( "sTlep_PAS"                       ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sTlep_PASandMee100"              ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sTjet_PAS"                       ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sTjet_PASandMee100"              ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sT_PAS"                          ,    300   , 0       , 3000	  ); 
   CreateUserTH1D( "sT_zjj_PAS"                      ,    300   , 0       , 3000	  ); 
   CreateUserTH1D( "sT_zjj_PASandMee100"             ,    300   , 0       , 3000	  ); 
   CreateUserTH1D( "sT_PASandMee100"                 ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sT_PASandMee110"                 ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sT_PASandMee120"                 ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sT_PASandMee130"                 ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sT_PASandMee140"                 ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sT_PASandMee150"                 ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sT_PASandMee160"                 ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sT_PASandMee170"                 ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sT_PASandMee180"                 ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sT_PASandMee190"                 ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sT_PASandMee200"                 ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mjj_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mjj_PASandMee100"	             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_PASandST445"                 ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_PAS"                       ,    200   , 0       , 1000	  ); 
   CreateUserTH1D( "Me1j1_PAS"                       ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Me1j2_PAS"                       ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Me2j1_PAS"                       ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Me2j2_PAS"                       ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Me1j_selected_PAS"               ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Me2j_selected_PAS"               ,    200   , 0       , 2000   );
   CreateUserTH1D( "Mej_selected_avg_PAS"            ,    200   , 0       , 2000   ); 
   CreateUserTH1D( "Mej_selected_min_PAS"            ,    200   , 0       , 2000   ); 
   CreateUserTH1D( "Mej_selected_max_PAS"            ,    200   , 0       , 2000   ); 
   CreateUserTH1D( "Mej_minmax_PAS"                  ,    200   , 0       , 2000   ); 
   CreateUserTH1D( "Meejj_PAS"                       ,    400   , 0       , 4000   );
   CreateUserTH1D( "Mejj_PAS"                        ,    400   , 0       , 4000   );
   CreateUserTH1D( "Meej_PAS"                        ,    400   , 0       , 4000   );
   CreateUserTH1D( "run_PAS"                         ,    20000 , 270000  , 290000 );
   CreateUserTH1D( "run_HLT"                         ,    20000 , 270000  , 290000 );
						     
   CreateUserTH1D( "Ptj1j2j3_PAS"                    ,    400 , 0       , 4000     );
   CreateUserTH1D( "Ptj1j2_PAS"                      ,    400 , 0       , 4000     );
   CreateUserTH1D( "Ptj2j3_PAS"                      ,    400 , 0       , 4000     );
   CreateUserTH1D( "Ptj1j3_PAS"                      ,    400 , 0       , 4000     );
						     
   CreateUserTH1D( "Ptee_Minus_Ptj1j2_PAS"           ,    200 , -500    , 500      );
   CreateUserTH1D( "Ptee_Minus_Ptj1j2j3_PAS"         ,    200 , -500    , 500      );
   CreateUserTH1D( "Ptj1j2j3_PASandMee100"           ,    200 , 0       , 2000     );
   CreateUserTH1D( "Ptj1j2_PASandMee100"             ,    200 , 0       , 2000     );
   CreateUserTH1D( "Ptj2j3_PASandMee100"             ,    200 , 0       , 2000     );
   CreateUserTH1D( "Ptj1j3_PASandMee100"             ,    200 , 0       , 2000     );
   CreateUserTH1D( "Ptee_Minus_Ptj1j2_PASandMee100"  ,    200 , -500    , 500      );
   CreateUserTH1D( "Ptee_Minus_Ptj1j2j3_PASandMee100",    200 , -500    , 500      );
   CreateUserTH1D( "Ptee_PAS"                        ,    200 , 0       , 2000     );
   CreateUserTH1D( "Ptee_PASandMee100"               ,    200 , 0       , 2000     );
					             
   CreateUserTH1D( "M_j1j3_PAS"                      ,    200 , 0       , 2000	 );    
   CreateUserTH1D( "M_j2j3_PAS"                      ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "M_e1j3_PAS"                      ,    200 , 0       , 2000	 );    
   CreateUserTH1D( "M_e2j3_PAS"                      ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "M_eejjj_PAS"                     ,    500 , 0       , 5000	 ); 
					             
   CreateUserTH1D( "Me1j1_PASandMee100"             ,    200 , 0       , 2000   );
   CreateUserTH1D( "Me1j2_PASandMee100"             ,    200 , 0       , 2000   );
   CreateUserTH1D( "Me2j1_PASandMee100"             ,    200 , 0       , 2000   );
   CreateUserTH1D( "Me2j2_PASandMee100"             ,    200 , 0       , 2000   );

   CreateUserTH1D( "M_j1j3_PASandMee100"             ,    200 , 0       , 2000	 );    
   CreateUserTH1D( "M_j2j3_PASandMee100"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "M_e1j3_PASandMee100"             ,    200 , 0       , 2000	 );    
   CreateUserTH1D( "M_e2j3_PASandMee100"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "M_eejjj_PASandMee100"            ,    500 , 0       , 5000	 ); 
   
   CreateUserTH1D( "DCotTheta1stEle_PAS"             ,    100 , 0.0, 1.0);
   CreateUserTH1D( "Dist1stEle_PAS"                  ,    100 , 0.0, 1.0);  
   CreateUserTH1D( "DCotTheta2ndEle_PAS"             ,    100 , 0.0, 1.0);
   CreateUserTH1D( "Dist2ndEle_PAS"                  ,    100 , 0.0, 1.0);  
		                                     
   CreateUserTH1D( "nVertex_PAS"                     ,    101   , -0.5   , 100.5	 ) ; 
   CreateUserTH1D( "nVertex_PASandMee100"            ,    101   , -0.5   , 100.5	 ) ; 
   
   CreateUserTH1D( "DR_Ele1Jet1_PAS"	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele1Jet2_PAS"	   , 	getHistoNBins("DR_Ele1Jet2"), getHistoMin("DR_Ele1Jet2"), getHistoMax("DR_Ele1Jet2")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet1_PAS"	   , 	getHistoNBins("DR_Ele2Jet1"), getHistoMin("DR_Ele2Jet1"), getHistoMax("DR_Ele2Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet2_PAS"	   , 	getHistoNBins("DR_Ele2Jet2"), getHistoMin("DR_Ele2Jet2"), getHistoMax("DR_Ele2Jet2")     ) ; 
   CreateUserTH1D( "DR_Jet1Jet2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
   CreateUserTH1D( "DR_Ele1Ele2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
   CreateUserTH1D( "minDR_EleJet_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
   CreateUserTH1D( "minDR_ZJet_PAS"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
   
   CreateUserTH1D( "DR_ZJet1_PAS"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
   CreateUserTH1D( "DR_ZJet2_PAS"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 

   CreateUserTH2D( "MeeVsST_PAS"                 ,     400, 0, 4000, 400, 0, 4000) ;
   CreateUserTH2D( "MeeVsST_PASandMee100"        ,     400, 0, 4000, 400, 0, 4000) ;
   CreateUserTH2D( "MeeVsPtee_PAS"                 ,     400, 0, 4000, 400, 0, 4000) ;

   CreateUserTH2D( "Me1jVsMe2j_selected",     200, 0, 2000, 200, 0, 2000) ;
   CreateUserTH2D( "Me1jVsMe2j_rejected",     200, 0, 2000, 200, 0, 2000) ;

   CreateUserTH1D( "MTeemunu_PAS"          ,    200 , 0       , 1000	 ); 

   CreateUserTH1D( "Mee_80_100_Preselection", 200, 60, 120 );
   CreateUserTH1D( "Mee_70_110_Preselection", 200, 60, 120 );
   CreateUserTH1D( "Mee_70_110_ST600_Preselection", 200, 60, 120 );

   CreateUserTH1D( "Mee_70_110_Preselection_Process0", 200, 60, 120 );
   CreateUserTH1D( "Mee_70_110_Preselection_Process1", 200, 60, 120 );
   CreateUserTH1D( "Mee_70_110_Preselection_Process2", 200, 60, 120 );
   CreateUserTH1D( "Mee_70_110_Preselection_Process3", 200, 60, 120 );
   CreateUserTH1D( "Mee_70_110_Preselection_Process4", 200, 60, 120 );

   CreateUserTH1D( "Mee_EBEB_PAS"		   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EBEE_PAS"		   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EEEE_PAS"		   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EB_PAS" 		   ,    60 , 60       , 120	 ); 

   CreateUserTH1D( "Mee_EBEB_80_100_PAS"	   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EBEE_80_100_PAS"	   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EEEE_80_100_PAS"	   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EB_80_100_PAS" 	     	   ,    60 , 60       , 120	 ); 
   
   CreateUserTH1D( "Mee_EBEB_70_110_PAS"	   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EBEE_70_110_PAS"	   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EEEE_70_110_PAS"	   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EB_70_110_PAS" 	     	   ,    60 , 60       , 120	 ); 

   CreateUserTH1D( "PileupWeight"   , 100, -10, 10 );
   CreateUserTH1D( "GeneratorWeight", 100, -2.0, 2.0);

   // basic muon kinematics
   CreateUserTH1D( "Pt1stMuon_PAS"	             , 	100    , 0       , 1000     ); 
   CreateUserTH1D( "Eta1stMuon_PAS"	             , 	100    , -5      , 5	  ); 
   CreateUserTH1D( "Phi1stMuon_PAS"	             , 	60     , -3.1416 , +3.1416  ); 
   CreateUserTH1D( "Pt2ndMuon_PAS"	             , 	100    , 0       , 1000     ); 
   CreateUserTH1D( "Eta2ndMuon_PAS"	             , 	100    , -5      , 5	  ); 
   CreateUserTH1D( "Phi2ndMuon_PAS"	             , 	60     , -3.1416 , +3.1416  ); 


   CreateUserTH1D("BeamSpotDXY_1stEle_PAS"                   , 200,  0.0 ,   0.5  ); CreateUserTH1D("BeamSpotDXY_2ndEle_PAS"                   , 200,  0.0 ,   0.5  );
   CreateUserTH1D("Classif_1stEle_PAS"                       , 5  , -0.5 ,   4.5  ); CreateUserTH1D("Classif_2ndEle_PAS"                       , 5  , -0.5 ,   4.5  );
   CreateUserTH1D("CorrIsolation_1stEle_PAS"                 , 200,-25.0 ,  25.0  ); CreateUserTH1D("CorrIsolation_2ndEle_PAS"                 , 200,-25.0 ,  25.0  );
   CreateUserTH1D("DeltaEtaTrkSC_1stEle_PAS"                 , 200, -0.01,   0.01 ); CreateUserTH1D("DeltaEtaTrkSC_2ndEle_PAS"                 , 200, -0.01,   0.01 );
   CreateUserTH1D("DeltaPhiTrkSC_1stEle_PAS"                 , 200, -0.1 ,   0.1  ); CreateUserTH1D("DeltaPhiTrkSC_2ndEle_PAS"                 , 200, -0.1 ,   0.1  );
   CreateUserTH1D("Full5x5E1x5OverE5x5_1stEle_PAS"           , 200,  0.0 ,   2.0  ); CreateUserTH1D("Full5x5E1x5OverE5x5_2ndEle_PAS"           , 200,  0.0 ,   2.0  );
   CreateUserTH1D("Full5x5E2x5OverE5x5_1stEle_PAS"           , 200,  0.0 ,   2.0  ); CreateUserTH1D("Full5x5E2x5OverE5x5_2ndEle_PAS"           , 200,  0.0 ,   2.0  );
   CreateUserTH1D("EcalIsolation_1stEle_PAS"                 , 200,  0.0 ,  20.0  ); CreateUserTH1D("EcalIsolation_2ndEle_PAS"                 , 200,  0.0 ,  20.0  );
   CreateUserTH1D("HcalIsolation_1stEle_PAS"                 , 200,  0.0 ,  20.0  ); CreateUserTH1D("HcalIsolation_2ndEle_PAS"                 , 200,  0.0 ,  20.0  );
   CreateUserTH1D("TrkIsolation_1stEle_PAS"                  , 200,  0.0,    5.0  ); CreateUserTH1D("TrkIsolation_2ndEle_PAS"                  , 200,  0.0,    5.0  );
   CreateUserTH1D("Energy_1stEle_PAS"                        , 200,  0.0 ,3000.0  ); CreateUserTH1D("Energy_2ndEle_PAS"                        , 200,  0.0 ,3000.0  );
   CreateUserTH1D("FBrem_1stEle_PAS"                         , 200,-10.0 ,  10.0  ); CreateUserTH1D("FBrem_2ndEle_PAS"                         , 200,-10.0 ,  10.0  );
   CreateUserTH1D("GsfCtfCharge_1stEle_PAS"                  , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfCtfCharge_2ndEle_PAS"                  , 2,   -0.5 ,   1.5  );
   CreateUserTH1D("GsfCtfScPixCharge_1stEle_PAS"             , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfCtfScPixCharge_2ndEle_PAS"             , 2,   -0.5 ,   1.5  );
   CreateUserTH1D("GsfScPixCharge_1stEle_PAS"                , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfScPixCharge_2ndEle_PAS"                , 2,   -0.5 ,   1.5  );
   CreateUserTH1D("HasMatchedPhot_1stEle_PAS"                , 2,   -0.5 ,   1.5  ); CreateUserTH1D("HasMatchedPhot_2ndEle_PAS"                , 2,   -0.5 ,   1.5  );
   CreateUserTH1D("HoE_1stEle_PAS"                           , 200,  0.0 ,   0.05 ); CreateUserTH1D("HoE_2ndEle_PAS"                           , 200,  0.0 ,   0.05 );
   CreateUserTH1D("LeadVtxDistXY_1stEle_PAS"                 , 200, -0.05,   0.05 ); CreateUserTH1D("LeadVtxDistXY_2ndEle_PAS"                 , 200, -0.05,   0.05 );
   CreateUserTH1D("LeadVtxDistZ_1stEle_PAS"                  , 200, -0.2 ,   0.2  ); CreateUserTH1D("LeadVtxDistZ_2ndEle_PAS"                  , 200, -0.2 ,   0.2  );
   CreateUserTH1D("MissingHits_1stEle_PAS"                   , 2  , -0.5,    1.5  ); CreateUserTH1D("MissingHits_2ndEle_PAS"                   , 2  , -0.5,    1.5  );
   CreateUserTH1D("NBrems_1stEle_PAS"                        , 11 , -0.5,   10.5  ); CreateUserTH1D("NBrems_2ndEle_PAS"                        , 11 , -0.5,   10.5  );
   CreateUserTH1D("EnergyORawEnergy_1stEle_PAS"              , 200,  0.9,    1.4  ); CreateUserTH1D("EnergyORawEnergy_2ndEle_PAS"              , 200,  0.9,    1.4  );
   CreateUserTH1D("SigmaEtaEta_Barrel_1stEle_PAS"            , 200,  0.0,    0.02 ); CreateUserTH1D("SigmaEtaEta_Barrel_2ndEle_PAS"            , 200,  0.0,    0.02 );
   CreateUserTH1D("SigmaEtaEta_Endcap_1stEle_PAS"            , 200,  0.0,    0.1  ); CreateUserTH1D("SigmaEtaEta_Endcap_2ndEle_PAS"            , 200,  0.0,    0.1  );
   CreateUserTH1D("Full5x5SigmaIEtaIEta_Barrel_1stEle_PAS"          , 200,  0.0,    0.04 ); CreateUserTH1D("Full5x5SigmaIEtaIEta_Barrel_2ndEle_PAS"          , 200,  0.0,    0.04 );
   CreateUserTH1D("Full5x5SigmaIEtaIEta_Endcap_1stEle_PAS"          , 200,  0.0,    0.1  ); CreateUserTH1D("Full5x5SigmaIEtaIEta_Endcap_2ndEle_PAS"          , 200,  0.0,    0.1  );
   CreateUserTH1D("TrkPtOPt_1stEle_PAS"                      , 200,  0.0,  100.0  ); CreateUserTH1D("TrkPtOPt_2ndEle_PAS"                      , 200,  0.0,  100.0  );
   CreateUserTH1D("ValidFrac_1stEle_PAS"                     , 200,  0.0 ,   2.0  ); CreateUserTH1D("ValidFrac_2ndEle_PAS"                     , 200,  0.0 ,   2.0  );
                                                                                                                                                                      
   CreateUserTH1D("BeamSpotDXY_1stEle_PASandMee100"          , 200,  0.0 ,   0.5  ); CreateUserTH1D("BeamSpotDXY_2ndEle_PASandMee100"          , 200,  0.0 ,   0.5  );
   CreateUserTH1D("Classif_1stEle_PASandMee100"              , 5  , -0.5 ,   4.5  ); CreateUserTH1D("Classif_2ndEle_PASandMee100"              , 5  , -0.5 ,   4.5  );
   CreateUserTH1D("CorrIsolation_1stEle_PASandMee100"        , 200,-25.0 ,  25.0  ); CreateUserTH1D("CorrIsolation_2ndEle_PASandMee100"        , 200,-25.0 ,  25.0  );
   CreateUserTH1D("DeltaEtaTrkSC_1stEle_PASandMee100"        , 200, -0.01,   0.01 ); CreateUserTH1D("DeltaEtaTrkSC_2ndEle_PASandMee100"        , 200, -0.01,   0.01 );
   CreateUserTH1D("DeltaPhiTrkSC_1stEle_PASandMee100"        , 200, -0.1 ,   0.1  ); CreateUserTH1D("DeltaPhiTrkSC_2ndEle_PASandMee100"        , 200, -0.1 ,   0.1  );
   CreateUserTH1D("Full5x5E1x5OverE5x5_1stEle_PASandMee100"  , 200,  0.0 ,   2.0  ); CreateUserTH1D("Full5x5E1x5OverE5x5_2ndEle_PASandMee100"  , 200,  0.0 ,   2.0  );
   CreateUserTH1D("Full5x5E2x5OverE5x5_1stEle_PASandMee100"  , 200,  0.0 ,   2.0  ); CreateUserTH1D("Full5x5E2x5OverE5x5_2ndEle_PASandMee100"  , 200,  0.0 ,   2.0  );
   CreateUserTH1D("EcalIsolation_1stEle_PASandMee100"        , 200,  0.0 ,  20.0  ); CreateUserTH1D("EcalIsolation_2ndEle_PASandMee100"        , 200,  0.0 ,  20.0  );
   CreateUserTH1D("HcalIsolation_1stEle_PASandMee100"        , 200,  0.0 ,  20.0  ); CreateUserTH1D("HcalIsolation_2ndEle_PASandMee100"        , 200,  0.0 ,  20.0  );
   CreateUserTH1D("TrkIsolation_1stEle_PASandMee100"         , 200,  0.0,    5.0  ); CreateUserTH1D("TrkIsolation_2ndEle_PASandMee100"         , 200,  0.0,    5.0  );
   CreateUserTH1D("Energy_1stEle_PASandMee100"               , 200,  0.0 ,3000.0  ); CreateUserTH1D("Energy_2ndEle_PASandMee100"               , 200,  0.0 ,3000.0  );
   CreateUserTH1D("FBrem_1stEle_PASandMee100"                , 200,-10.0 ,  10.0  ); CreateUserTH1D("FBrem_2ndEle_PASandMee100"                , 200,-10.0 ,  10.0  );
   CreateUserTH1D("GsfCtfCharge_1stEle_PASandMee100"         , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfCtfCharge_2ndEle_PASandMee100"         , 2,   -0.5 ,   1.5  );
   CreateUserTH1D("GsfCtfScPixCharge_1stEle_PASandMee100"    , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfCtfScPixCharge_2ndEle_PASandMee100"    , 2,   -0.5 ,   1.5  );
   CreateUserTH1D("GsfScPixCharge_1stEle_PASandMee100"       , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfScPixCharge_2ndEle_PASandMee100"       , 2,   -0.5 ,   1.5  );
   CreateUserTH1D("HasMatchedPhot_1stEle_PASandMee100"       , 2,   -0.5 ,   1.5  ); CreateUserTH1D("HasMatchedPhot_2ndEle_PASandMee100"       , 2,   -0.5 ,   1.5  );
   CreateUserTH1D("HoE_1stEle_PASandMee100"                  , 200,  0.0 ,   0.05 ); CreateUserTH1D("HoE_2ndEle_PASandMee100"                  , 200,  0.0 ,   0.05 );
   CreateUserTH1D("LeadVtxDistXY_1stEle_PASandMee100"        , 200, -0.05,   0.05 ); CreateUserTH1D("LeadVtxDistXY_2ndEle_PASandMee100"        , 200, -0.05,   0.05 );
   CreateUserTH1D("LeadVtxDistZ_1stEle_PASandMee100"         , 200, -0.2 ,   0.2  ); CreateUserTH1D("LeadVtxDistZ_2ndEle_PASandMee100"         , 200, -0.2 ,   0.2  );
   CreateUserTH1D("MissingHits_1stEle_PASandMee100"          , 2  , -0.5,    1.5  ); CreateUserTH1D("MissingHits_2ndEle_PASandMee100"          , 2  , -0.5,    1.5  );
   CreateUserTH1D("NBrems_1stEle_PASandMee100"               , 11 , -0.5,   10.5  ); CreateUserTH1D("NBrems_2ndEle_PASandMee100"               , 11 , -0.5,   10.5  );
   CreateUserTH1D("EnergyORawEnergy_1stEle_PASandMee100"     , 200,  0.9,    1.4  ); CreateUserTH1D("EnergyORawEnergy_2ndEle_PASandMee100"     , 200,  0.9,    1.4  );
   CreateUserTH1D("SigmaEtaEta_Barrel_1stEle_PASandMee100"   , 200,  0.0,    0.02 ); CreateUserTH1D("SigmaEtaEta_Barrel_2ndEle_PASandMee100"   , 200,  0.0,    0.02 );
   CreateUserTH1D("SigmaEtaEta_Endcap_1stEle_PASandMee100"   , 200,  0.0,    0.1  ); CreateUserTH1D("SigmaEtaEta_Endcap_2ndEle_PASandMee100"   , 200,  0.0,    0.1  );
   CreateUserTH1D("Full5x5SigmaIEtaIEta_Barrel_1stEle_PASandMee100" , 200,  0.0,    0.02 ); CreateUserTH1D("Full5x5SigmaIEtaIEta_Barrel_2ndEle_PASandMee100" , 200,  0.0,    0.02 );
   CreateUserTH1D("Full5x5SigmaIEtaIEta_Endcap_1stEle_PASandMee100" , 200,  0.0,    0.1  ); CreateUserTH1D("Full5x5SigmaIEtaIEta_Endcap_2ndEle_PASandMee100" , 200,  0.0,    0.1  );
   CreateUserTH1D("TrkPtOPt_1stEle_PASandMee100"             , 200,  0.0,  100.0  ); CreateUserTH1D("TrkPtOPt_2ndEle_PASandMee100"             , 200,  0.0,  100.0  );
   CreateUserTH1D("ValidFrac_1stEle_PASandMee100"            , 200,  0.0 ,   2.0  ); CreateUserTH1D("ValidFrac_2ndEle_PASandMee100"            , 200,  0.0 ,   2.0  );
                                                                                                                                                                      
   if(do_roi_plots) {										     										 
     CreateUserTH1D( "EleChargeSum_ROI"                ,    3     , -2.5    , 2.5      );
     CreateUserTH1D( "Mej_selected_avg_ROI"            ,    200   , 0       , 2000     );
     CreateUserTH1D( "sTfrac_Jet1_ROI"                 ,   100    ,  0.0    , 1.0      );
     CreateUserTH1D( "sTfrac_Jet2_ROI"                 ,   100    ,  0.0    , 1.0      );
     CreateUserTH1D( "sTfrac_Ele1_ROI"                 ,   100    ,  0.0    , 1.0      );
     CreateUserTH1D( "sTfrac_Ele2_ROI"                 ,   100    ,  0.0    , 1.0      );
     CreateUserTH1D( "sTfrac_Jet_ROI"                  ,   100    ,  0.0    , 1.0      );
     CreateUserTH1D( "sTfrac_Ele_ROI"                  ,   100    ,  0.0    , 1.0      );
     CreateUserTH1D( "nJet_ROI"                        ,    10    , -0.5    , 9.5      );
     CreateUserTH1D( "Pt1stEle_ROI"	             , 	100    , 0       , 1000     ); 
     CreateUserTH1D( "Eta1stEle_ROI"	             , 	100    , -5      , 5	  ); 
     CreateUserTH1D( "Phi1stEle_ROI"	             , 	60     , -3.1416 , +3.1416  ); 
     CreateUserTH1D( "Pt2ndEle_ROI"	             , 	300    , 0       , 3000     ); 
     CreateUserTH1D( "Eta2ndEle_ROI"	             , 	100    , -5      , 5	  ); 
     CreateUserTH1D( "Phi2ndEle_ROI"	             , 	60     , -3.1416 , +3.1416  ); 
     CreateUserTH1D( "MET_ROI"                         ,    200   , 0       , 1000	  ); 
     CreateUserTH1D( "Pt1stJet_ROI"                    ,    200   , 0       , 2000	  ); 
     CreateUserTH1D( "Pt2ndJet_ROI"                    ,    200   , 0       , 2000	  ); 
     CreateUserTH1D( "Eta1stJet_ROI"                   ,    100   , -5      , 5	  ); 
     CreateUserTH1D( "Eta2ndJet_ROI"                   ,    100   , -5      , 5	  ); 
     CreateUserTH1D( "Phi1stJet_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
     CreateUserTH1D( "Phi2ndJet_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
     CreateUserTH1D( "sTlep_ROI"                       ,    200   , 0       , 2000	  ); 
     CreateUserTH1D( "sTjet_ROI"                       ,    200   , 0       , 2000	  ); 
     CreateUserTH1D( "sT_zjj_ROI"                      ,    200   , 0       , 2000	  ); 
     CreateUserTH1D( "sT_ROI"                          ,    200   , 0       , 2000	  ); 
     CreateUserTH1D( "Mjj_ROI"		             ,    200   , 0       , 2000	  ); 
     CreateUserTH1D( "Mee_ROI"		             ,    200   , 0       , 2000	  ); 
     CreateUserTH1D( "Meejj_ROI"                       ,    400   , 0       , 4000   );
     CreateUserTH1D( "Mejj_ROI"                        ,    400   , 0       , 4000   );
     CreateUserTH1D( "Meej_ROI"                        ,    400   , 0       , 4000   );
     CreateUserTH1D( "Ptj1j2j3_ROI"                    ,    200 , 0       , 2000     );
     CreateUserTH1D( "Ptj1j2_ROI"                      ,    200 , 0       , 2000     );
     CreateUserTH1D( "Ptj2j3_ROI"                      ,    200 , 0       , 2000     );
     CreateUserTH1D( "Ptj1j3_ROI"                      ,    200 , 0       , 2000     );
     CreateUserTH1D( "Ptee_Minus_Ptj1j2_ROI"           ,    200 , -500    , 500      );
     CreateUserTH1D( "Ptee_Minus_Ptj1j2j3_ROI"         ,    200 , -500    , 500      );
     CreateUserTH1D( "Ptee_ROI"                        ,    200 , 0       , 2000     );
     CreateUserTH1D( "Me1j1_ROI"                      ,    200 , 0       , 2000   );
     CreateUserTH1D( "Me1j2_ROI"                      ,    200 , 0       , 2000   );
     CreateUserTH1D( "Me2j1_ROI"                      ,    200 , 0       , 2000   );
     CreateUserTH1D( "Me2j2_ROI"                      ,    200 , 0       , 2000   );
     CreateUserTH1D( "M_j1j3_ROI"                      ,    200 , 0       , 2000	 );    
     CreateUserTH1D( "M_j2j3_ROI"                      ,    200 , 0       , 2000	 ); 
     CreateUserTH1D( "M_e1j3_ROI"                      ,    200 , 0       , 2000	 );    
     CreateUserTH1D( "M_e2j3_ROI"                      ,    200 , 0       , 2000	 ); 
     CreateUserTH1D( "M_eejjj_ROI"                     ,    500 , 0       , 5000	 ); 
     CreateUserTH1D( "nVertex_ROI"                     ,    101   , -0.5   , 100.5	 ) ; 
     CreateUserTH1D( "minDR_ZJet_ROI"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
     CreateUserTH1D( "DR_ZJet1_ROI"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
     CreateUserTH1D( "DR_ZJet2_ROI"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
     CreateUserTH2D( "MeeVsST_ROI"                 ,     400, 0, 4000, 400, 0, 4000) ;
     CreateUserTH1D("BeamSpotDXY_1stEle_ROI"                   , 200,  0.0 ,   0.5  ); CreateUserTH1D("BeamSpotDXY_2ndEle_ROI"                   , 200,  0.0 ,   0.5  );
     CreateUserTH1D("Classif_1stEle_ROI"                       , 5  , -0.5 ,   4.5  ); CreateUserTH1D("Classif_2ndEle_ROI"                       , 5  , -0.5 ,   4.5  );
     CreateUserTH1D("CorrIsolation_1stEle_ROI"                 , 200,-25.0 ,  25.0  ); CreateUserTH1D("CorrIsolation_2ndEle_ROI"                 , 200,-25.0 ,  25.0  );
     CreateUserTH1D("DeltaEtaTrkSC_1stEle_ROI"                 , 200, -0.01,   0.01 ); CreateUserTH1D("DeltaEtaTrkSC_2ndEle_ROI"                 , 200, -0.01,   0.01 );
     CreateUserTH1D("DeltaPhiTrkSC_1stEle_ROI"                 , 200, -0.1 ,   0.1  ); CreateUserTH1D("DeltaPhiTrkSC_2ndEle_ROI"                 , 200, -0.1 ,   0.1  );
     CreateUserTH1D("E1x5OverE5x5_1stEle_ROI"                  , 200,  0.0 ,   2.0  ); CreateUserTH1D("E1x5OverE5x5_2ndEle_ROI"                  , 200,  0.0 ,   2.0  );
     CreateUserTH1D("E2x5OverE5x5_1stEle_ROI"                  , 200,  0.0 ,   2.0  ); CreateUserTH1D("E2x5OverE5x5_2ndEle_ROI"                  , 200,  0.0 ,   2.0  );
     CreateUserTH1D("EcalIsolation_1stEle_ROI"                 , 200,  0.0 ,  20.0  ); CreateUserTH1D("EcalIsolation_2ndEle_ROI"                 , 200,  0.0 ,  20.0  );
     CreateUserTH1D("HcalIsolation_1stEle_ROI"                 , 200,  0.0 ,  20.0  ); CreateUserTH1D("HcalIsolation_2ndEle_ROI"                 , 200,  0.0 ,  20.0  );
     CreateUserTH1D("TrkIsolation_1stEle_ROI"                  , 200,  0.0,    5.0  ); CreateUserTH1D("TrkIsolation_2ndEle_ROI"                  , 200,  0.0,    5.0  );
     CreateUserTH1D("Energy_1stEle_ROI"                        , 200,  0.0 ,3000.0  ); CreateUserTH1D("Energy_2ndEle_ROI"                        , 200,  0.0 ,3000.0  );
     CreateUserTH1D("FBrem_1stEle_ROI"                         , 200,-10.0 ,  10.0  ); CreateUserTH1D("FBrem_2ndEle_ROI"                         , 200,-10.0 ,  10.0  );
     CreateUserTH1D("GsfCtfCharge_1stEle_ROI"                  , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfCtfCharge_2ndEle_ROI"                  , 2,   -0.5 ,   1.5  );
     CreateUserTH1D("GsfCtfScPixCharge_1stEle_ROI"             , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfCtfScPixCharge_2ndEle_ROI"             , 2,   -0.5 ,   1.5  );
     CreateUserTH1D("GsfScPixCharge_1stEle_ROI"                , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfScPixCharge_2ndEle_ROI"                , 2,   -0.5 ,   1.5  );
     CreateUserTH1D("HasMatchedPhot_1stEle_ROI"                , 2,   -0.5 ,   1.5  ); CreateUserTH1D("HasMatchedPhot_2ndEle_ROI"                , 2,   -0.5 ,   1.5  );
     CreateUserTH1D("HoE_1stEle_ROI"                           , 200,  0.0 ,   0.05 ); CreateUserTH1D("HoE_2ndEle_ROI"                           , 200,  0.0 ,   0.05 );
     CreateUserTH1D("LeadVtxDistXY_1stEle_ROI"                 , 200, -0.05,   0.05 ); CreateUserTH1D("LeadVtxDistXY_2ndEle_ROI"                 , 200, -0.05,   0.05 );
     CreateUserTH1D("LeadVtxDistZ_1stEle_ROI"                  , 200, -0.2 ,   0.2  ); CreateUserTH1D("LeadVtxDistZ_2ndEle_ROI"                  , 200, -0.2 ,   0.2  );
     CreateUserTH1D("MissingHits_1stEle_ROI"                   , 2  , -0.5,    1.5  ); CreateUserTH1D("MissingHits_2ndEle_ROI"                   , 2  , -0.5,    1.5  );
     CreateUserTH1D("NBrems_1stEle_ROI"                        , 11 , -0.5,   10.5  ); CreateUserTH1D("NBrems_2ndEle_ROI"                        , 11 , -0.5,   20.5  );
     CreateUserTH1D("EnergyORawEnergy_1stEle_ROI"              , 200,  0.9,    1.4  ); CreateUserTH1D("EnergyORawEnergy_2ndEle_ROI"              , 200,  0.9,    1.4  );
     CreateUserTH1D("SigmaEtaEta_Barrel_1stEle_ROI"            , 200,  0.0,    0.02 ); CreateUserTH1D("SigmaEtaEta_Barrel_2ndEle_ROI"            , 200,  0.0,    0.02 );
     CreateUserTH1D("SigmaEtaEta_Endcap_1stEle_ROI"            , 200,  0.0,    0.1  ); CreateUserTH1D("SigmaEtaEta_Endcap_2ndEle_ROI"            , 200,  0.0,    0.1  );
     CreateUserTH1D("SigmaIEtaIEta_Barrel_1stEle_ROI"          , 200,  0.0,    0.02 ); CreateUserTH1D("SigmaIEtaIEta_Barrel_2ndEle_ROI"          , 200,  0.0,    0.02 );
     CreateUserTH1D("SigmaIEtaIEta_Endcap_1stEle_ROI"          , 200,  0.0,    0.1  ); CreateUserTH1D("SigmaIEtaIEta_Endcap_2ndEle_ROI"          , 200,  0.0,    0.1  );
     CreateUserTH1D("TrkPtOPt_1stEle_ROI"                      , 200,  0.0,  100.0  ); CreateUserTH1D("TrkPtOPt_2ndEle_ROI"                      , 200,  0.0,  100.0  );
     CreateUserTH1D("ValidFrac_1stEle_ROI"                     , 200,  0.0 ,   2.0  ); CreateUserTH1D("ValidFrac_2ndEle_ROI"                     , 200,  0.0 ,   2.0  );
   }
   
   // for scale factor dependence studies
   CreateUserTH1D( "Mee_NJetEq2_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_NJetEq3_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_NJetEq4_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_NJetEq5_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_NJetEq6_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_NJetEq7_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_NJetGeq3_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_NJetGeq4_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT300To500_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT500To750_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT750To1250_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1250ToInf_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_MejMin100To200_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_MejMin200To300_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_MejMin300To400_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_MejMin400To500_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_MejMin500To650_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_MejMin650ToInf_PAS"		             ,    200   , 0       , 2000	  ); 
   // for Ptee bins
   CreateUserTH1D( "Mee_Ptee0To100_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_Ptee100To150_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_Ptee150To200_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_Ptee200To250_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_Ptee250To300_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_Ptee300To350_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_Ptee350To400_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_Ptee400ToInf_PAS"		             ,    200   , 0       , 2000	  ); 
   // more plots with sT cuts for final selection thresholds
   CreateUserTH1D( "Mee_sT340_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT405_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT470_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT535_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT595_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT660_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT720_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT780_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT840_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT900_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT960_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1015_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1075_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1130_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1190_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1245_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1300_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1355_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1410_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1460_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1515_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1565_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1615_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1670_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1720_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1770_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_sT1815_PAS"		             ,    200   , 0       , 2000	  ); 
   //// test opt
   //CreateUserTH1D( "Mee_sT2000_PAS"		             ,    200   , 0       , 2000	  ); 
   //CreateUserTH1D( "OptBinLQ600", 200, 0, 2000);
   //CreateUserTH1D( "PtZforOptBin600",1000,0,2000);
   //CreateUserTH1D( "PtEEforOptBin600",1000,0,2000);
   //CreateUserTH2D( "PtEEVsZPtForOptBin600",1000,0,2000,1000,0,2000);
   //CreateUserTH1D( "OptBinLQ650", 200, 0, 2000);
   //CreateUserTH1D( "OptBinLQ700", 200, 0, 2000);
   //CreateUserTH1D( "OptBinLQ600_noWeight", 200, 0, 2000);
   //CreateUserTH1D( "OptBinLQ650_noWeight", 200, 0, 2000);
   //CreateUserTH1D( "OptBinLQ700_noWeight", 200, 0, 2000);
   // 3D opt cut space
   CreateUserTH3D( "OptimizationCutSpace", 200, 0, 2000, 200, 0, 2000, 200, 0, 2000);

   //CreateUserTH1D( "WZ_system_Pt" , 200,0,2000);
   //CreateUserTH1D( "WZ_system_Pt_GenEle" , 200,0,2000);
   //CreateUserTH1D( "WZ_system_Pt_WZPtCut" , 200,0,2000);
   // checking electrons
   CreateUserTH1D( "SCEta_1stEle_Presel",4000,0,10);
   CreateUserTH1D( "EleEta_1stEle_Presel",4000,0,10);
   CreateUserTH1D( "DeltaEtaTrkSC_1stEle_Presel",4000,0,10);
   CreateUserTH1D( "SCEtaMinusEleEta_1stEle_Presel",4000,0,10);
   CreateUserTH1D( "DeltaEtaTrkSC_Minus_SCEtaMinusEleEta_1stEle_Presel",4000,0,10);
   CreateUserTH1D( "SCEta_2ndEle_Presel",4000,0,10);
   CreateUserTH1D( "EleEta_2ndEle_Presel",4000,0,10);
   CreateUserTH1D( "DeltaEtaTrkSC_2ndEle_Presel",4000,0,10);
   CreateUserTH1D( "SCEtaMinusEleEta_2ndEle_Presel",4000,0,10);
   CreateUserTH1D( "DeltaEtaTrkSC_Minus_SCEtaMinusEleEta_2ndEle_Presel",4000,0,10);

   //--------------------------------------------------------------------------
   // Final selection plots
   //--------------------------------------------------------------------------
   bool doFinalSelections = false;
   // check if there is a final Mej specific in cutfile for any LQ mass
   for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
     int lq_mass = LQ_MASS[i_lq_mass];
     sprintf(cut_name, "min_M_ej_LQ%d"   , lq_mass );
     if(hasCut(cut_name)) {
       doFinalSelections = true;
       break;
     }
   }
   // now, we must have an Mej cut and optimization must be off to have final selections enabled
   doFinalSelections = doFinalSelections && !isOptimizationEnabled();
   
   char plot_name[100];
   
   if(doFinalSelections)
   {
     for (int i_lq_mass = 0; i_lq_mass < n_lq_mass ; ++i_lq_mass ) { 
       int lq_mass = LQ_MASS[i_lq_mass];
       sprintf(plot_name, "Mej_selected_avg_LQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 60  , 0 , 3000 );
       sprintf(plot_name, "Mej_selected_min_LQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 60  , 0 , 3000 );
       sprintf(plot_name, "Mej_selected_max_LQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 60  , 0 , 3000 );
       sprintf(plot_name, "Mej_minmax_LQ%d"             , lq_mass ); CreateUserTH1D ( plot_name, 60  , 0 , 3000 );
       sprintf(plot_name, "sT_eejj_LQ%d"                , lq_mass ); CreateUserTH1D ( plot_name, 30  , 0 , 3000 );
       sprintf(plot_name, "Mee_LQ%d"                    , lq_mass ); CreateUserTH1D ( plot_name, 40  , 0 , 2000 );
       sprintf(plot_name, "Mej_selected_min_vs_max_LQ%d", lq_mass ); CreateUserTH2D ( plot_name, 150  , 0 , 3000, 150  , 0 , 3000 );
       sprintf(plot_name, "DR_Ele1Jet1_LQ%d"            , lq_mass ); CreateUserTH1D ( plot_name, 
           getHistoNBins("DR_Ele1Jet1"), 
           getHistoMin  ("DR_Ele1Jet1"), 
           getHistoMax  ("DR_Ele1Jet1"));

       sprintf(plot_name, "nElectron_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name , 5,  -0.5 ,   4.5  );
       sprintf(plot_name, "BeamSpotDXY_1stEle_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   0.5  );
       sprintf(plot_name, "Classif_1stEle_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name , 5  , -0.5 ,   4.5  );
       sprintf(plot_name, "CorrIsolation_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200,-25.0 ,  25.0  );
       sprintf(plot_name, "DeltaEtaTrkSC_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.01,   0.01 );
       sprintf(plot_name, "DeltaPhiTrkSC_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.1 ,   0.1  );
       sprintf(plot_name, "Full5x5E1x5OverE5x5_1stEle_LQ%d"  , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );
       sprintf(plot_name, "Full5x5E2x5OverE5x5_1stEle_LQ%d"  , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );
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


       sprintf(plot_name, "BeamSpotDXY_2ndEle_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   0.5  );
       sprintf(plot_name, "Classif_2ndEle_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name , 5  , -0.5 ,   4.5  );
       sprintf(plot_name, "CorrIsolation_2ndEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200,-25.0 ,  25.0  );
       sprintf(plot_name, "DeltaEtaTrkSC_2ndEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.01,   0.01 );
       sprintf(plot_name, "DeltaPhiTrkSC_2ndEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.1 ,   0.1  );
       sprintf(plot_name, "Full5x5E1x5OverE5x5_2ndEle_LQ%d"  , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );
       sprintf(plot_name, "Full5x5E2x5OverE5x5_2ndEle_LQ%d"  , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );
       sprintf(plot_name, "EcalIsolation_2ndEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,  20.0  );
       sprintf(plot_name, "HcalIsolation_2ndEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,  20.0  );
       sprintf(plot_name, "TrkIsolation_2ndEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    5.0  );
       sprintf(plot_name, "Energy_2ndEle_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,3000.0  );
       sprintf(plot_name, "FBrem_2ndEle_LQ%d"                , lq_mass ); CreateUserTH1D( plot_name , 200,-10.0 ,  10.0  );
       sprintf(plot_name, "GsfCtfCharge_2ndEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
       sprintf(plot_name, "GsfCtfScPixCharge_2ndEle_LQ%d"    , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
       sprintf(plot_name, "GsfScPixCharge_2ndEle_LQ%d"       , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
       sprintf(plot_name, "HasMatchedPhot_2ndEle_LQ%d"       , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
       sprintf(plot_name, "HoE_2ndEle_LQ%d"                  , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   0.05 );
       sprintf(plot_name, "LeadVtxDistXY_2ndEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.05,   0.05 );
       sprintf(plot_name, "LeadVtxDistZ_2ndEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 200, -0.2 ,   0.2  );
       sprintf(plot_name, "MissingHits_2ndEle_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name , 2  , -0.5,    1.5  );
       sprintf(plot_name, "NBrems_2ndEle_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name , 11 , -0.5,   10.5  );
       sprintf(plot_name, "EnergyORawEnergy_2ndEle_LQ%d"     , lq_mass ); CreateUserTH1D( plot_name , 200,  0.9,    1.4  );
       sprintf(plot_name, "SigmaEtaEta_Barrel_2ndEle_LQ%d"   , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.02 );
       sprintf(plot_name, "SigmaEtaEta_Endcap_2ndEle_LQ%d"   , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.1  );
       sprintf(plot_name, "Full5x5SigmaIEtaIEta_Barrel_2ndEle_LQ%d" , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.02 );
       sprintf(plot_name, "Full5x5SigmaIEtaIEta_Endcap_2ndEle_LQ%d" , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.1  );
       sprintf(plot_name, "TrkPtOPt_2ndEle_LQ%d"             , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,  100.0  );
       sprintf(plot_name, "ValidFrac_2ndEle_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );


       sprintf(plot_name, "EleChargeSum_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name ,      3 , -2.5    , 2.5      );
       sprintf(plot_name, "sTfrac_Jet1_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
       sprintf(plot_name, "sTfrac_Jet2_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
       sprintf(plot_name, "sTfrac_Ele1_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
       sprintf(plot_name, "sTfrac_Ele2_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
       sprintf(plot_name, "sTfrac_Jet_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
       sprintf(plot_name, "sTfrac_Ele_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
       sprintf(plot_name, "nJet_LQ%d"                , lq_mass ); CreateUserTH1D( plot_name ,     10 , -0.5    , 9.5      );
       sprintf(plot_name, "Pt1stEle_LQ%d"	           , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
       sprintf(plot_name, "Eta1stEle_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
       sprintf(plot_name, "Phi1stEle_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , +3.1416  ); 
       sprintf(plot_name, "Pt2ndEle_LQ%d"	           , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
       sprintf(plot_name, "Eta2ndEle_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
       sprintf(plot_name, "Phi2ndEle_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , 3.1416   ); 
       sprintf(plot_name, "MET_LQ%d"                 , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 1000     ); 
       sprintf(plot_name, "Pt1stJet_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
       sprintf(plot_name, "Pt2ndJet_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
       sprintf(plot_name, "Eta1stJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
       sprintf(plot_name, "Eta2ndJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
       sprintf(plot_name, "Phi1stJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , 3.1416   ); 
       sprintf(plot_name, "Phi2ndJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , 3.1416   ); 
       sprintf(plot_name, "sTlep_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
       sprintf(plot_name, "sTjet_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
       sprintf(plot_name, "sT_zjj_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
       sprintf(plot_name, "Mjj_LQ%d"		   , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
       sprintf(plot_name, "Meejj_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    400 ,  0.0    , 4000     );
       sprintf(plot_name, "Mejj_LQ%d"                , lq_mass ); CreateUserTH1D( plot_name ,    400 ,  0.0    , 4000     );
       sprintf(plot_name, "Meej_LQ%d"                , lq_mass ); CreateUserTH1D( plot_name ,    400 ,  0.0    , 4000     );
       sprintf(plot_name, "Ptj1j2j3_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
       sprintf(plot_name, "Ptj1j2_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
       sprintf(plot_name, "Ptj2j3_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
       sprintf(plot_name, "Ptj1j3_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
       sprintf(plot_name, "Ptee_Minus_Ptj1j2_LQ%d"   , lq_mass ); CreateUserTH1D( plot_name ,    200 , -500    ,  500     );
       sprintf(plot_name, "Ptee_Minus_Ptj1j2j3_LQ%d" , lq_mass ); CreateUserTH1D( plot_name ,    200 , -500    ,  500     );
       sprintf(plot_name, "Ptee_LQ%d"                , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
       sprintf(plot_name, "Me1j1_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
       sprintf(plot_name, "Me1j2_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
       sprintf(plot_name, "Me2j1_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
       sprintf(plot_name, "Me2j2_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
       sprintf(plot_name, "M_j1j3_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );    
       sprintf(plot_name, "M_j2j3_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
       sprintf(plot_name, "M_e1j3_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );    
       sprintf(plot_name, "M_e2j3_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
       sprintf(plot_name, "M_eejjj_LQ%d"             , lq_mass ); CreateUserTH1D( plot_name ,    500 ,  0.0    , 5000     ); 
       sprintf(plot_name, "nVertex_LQ%d"             , lq_mass ); CreateUserTH1D( plot_name ,    101 , -0.5    ,  100.5   ); 
       sprintf(plot_name, "MeeVsST_LQ%d"             , lq_mass ); CreateUserTH2D( plot_name ,    200 ,  0.0, 2000, 200, 0, 2000) ;
       sprintf(plot_name, "minDR_ZJet_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
       sprintf(plot_name, "DR_ZJet1_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
       sprintf(plot_name, "DR_ZJet2_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
       // checking electrons
       sprintf(plot_name, "SCEta_1stEle_LQ%d", lq_mass ); CreateUserTH1D( plot_name , 4000,0,10);
       sprintf(plot_name, "EleEta_1stEle_LQ%d", lq_mass ); CreateUserTH1D( plot_name , 4000,0,10);
       sprintf(plot_name, "SCEtaMinusEleEta_1stEle_LQ%d", lq_mass ); CreateUserTH1D( plot_name , 4000,0,10);
       sprintf(plot_name, "DeltaEtaTrkSC_Minus_SCEtaMinusEleEta_1stEle_LQ%d", lq_mass ); CreateUserTH1D( plot_name , 4000,0,10);
       sprintf(plot_name, "SCEta_2ndEle_LQ%d", lq_mass ); CreateUserTH1D( plot_name , 4000,0,10);
       sprintf(plot_name, "EleEta_2ndEle_LQ%d", lq_mass ); CreateUserTH1D( plot_name , 4000,0,10);
       sprintf(plot_name, "SCEtaMinusEleEta_2ndEle_LQ%d", lq_mass ); CreateUserTH1D( plot_name , 4000,0,10);
       sprintf(plot_name, "DeltaEtaTrkSC_Minus_SCEtaMinusEleEta_2ndEle_LQ%d", lq_mass ); CreateUserTH1D( plot_name , 4000,0,10);
       // muon kinematics
       sprintf(plot_name, "Pt1stMuon_LQ%d"	           , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
       sprintf(plot_name, "Eta1stMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
       sprintf(plot_name, "Phi1stMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , +3.1416  ); 
       sprintf(plot_name, "Pt2ndMuon_LQ%d"	           , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
       sprintf(plot_name, "Eta2ndMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
       sprintf(plot_name, "Phi2ndMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , 3.1416   ); 
     }
     // for SF at final selections
     CreateUserTH1D( "Mee_70_110_LQ300", 200, 60, 120 );
     CreateUserTH1D( "Mee_70_110_LQ600", 200, 60, 120 );
     CreateUserTH1D( "Mee_70_110_LQ800", 200, 60, 120 );
     CreateUserTH1D( "Mee_70_110_LQ900", 200, 60, 120 );
     CreateUserTH1D( "Mee_70_110_LQ1000", 200, 60, 120 );
   }

   //--------------------------------------------------------------------------
   // Tell the user how many entries we'll look at
   //--------------------------------------------------------------------------

   Long64_t nentries = GetTreeEntries();
   std::cout << "analysisClass::analysisClass(): nentries = " << nentries << std::endl;

   //--------------------------------------------------------------------------
   // Loop over the chain
   //--------------------------------------------------------------------------
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     readerTools_->LoadEntry(jentry);
     //-----------------------------------------------------------------
     // Print progress
     //-----------------------------------------------------------------
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass:Loop(): jentry = " << jentry << "/" << nentries << std::endl;
     //// run ls event
     //std::cout << static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(ls) << " " << static_cast<unsigned int>(event) << std::endl;


     //--------------------------------------------------------------------------
     // Reset the cuts
     //--------------------------------------------------------------------------

     resetCuts();

     //--------------------------------------------------------------------------
     // Check good run list
     //--------------------------------------------------------------------------
     
     double run = readerTools_->ReadValueBranch<Double_t>("run");
     int passedJSON = passJSON ( run,
         readerTools_->ReadValueBranch<Double_t>("ls"),
         isData() ) ;

     //--------------------------------------------------------------------------
     // Do pileup re-weighting
     //--------------------------------------------------------------------------
     
     double pileup_weight = readerTools_->ReadValueBranch<Double_t>("puWeight");
     
     //--------------------------------------------------------------------------
     // Get information about gen-level reweighting (should be for Sherpa only)
     //--------------------------------------------------------------------------

     double gen_weight = readerTools_->ReadValueBranch<Double_t>("Weight");
     if ( isData() ) gen_weight = 1.0;
     //if ( isData && Ele2_ValidFrac > 998. ){
     //  gen_weight = 0.0;
     //  if      (  60.0 < M_e1e2 < 120. ) gen_weight = 0.61;
     //  else if ( 120.0 < M_e1e2 < 200. ) gen_weight = 0.42;
     //  else if ( 200.0 < M_e1e2        ) gen_weight = 0.42;
     //}

     // std::cout << "Gen weight = " << int ( 1.0 / gen_weight ) << std::endl;
     //std::cout << "Gen weight = " << gen_weight << std::endl;

     std::string current_file_name ( readerTools_->GetTree()->GetCurrentFile()->GetName());

     // SIC remove March 2018
     //// TopPt reweight
     //// only valid for powheg
     //if(current_file_name.find("TT_") != std::string::npos) {
     //  gen_weight*=TopPtWeight;
     //}
    
     // Electron scale factors for MC only
     if(!isData())
     {
       float recoSFEle1 = ElectronScaleFactors2016::LookupRecoSF(readerTools_->ReadValueBranch<Double_t>("Ele1_SCEta"));
       float recoSFEle2 = ElectronScaleFactors2016::LookupRecoSF(readerTools_->ReadValueBranch<Double_t>("Ele2_SCEta"));
       float heepSFEle1 = ElectronScaleFactors2016::LookupHeepSF(readerTools_->ReadValueBranch<Double_t>("Ele1_SCEta"));
       float heepSFEle2 = ElectronScaleFactors2016::LookupHeepSF(readerTools_->ReadValueBranch<Double_t>("Ele2_SCEta"));
       float totalScaleFactor = recoSFEle1*recoSFEle2*heepSFEle1*heepSFEle2;
       gen_weight*=totalScaleFactor;
     }

     //--------------------------------------------------------------------------
     // First variable to fill just shows the "reweighting".  Always passes.
     //--------------------------------------------------------------------------

     fillVariableWithValue ( "Reweighting", 1, gen_weight * pileup_weight  );

     ////--------------------------------------------------------------------------
     //// Special treatment of inclusive W/Z
     ////--------------------------------------------------------------------------
     //bool passGenWZPt = true;
     //// inclusive
     ////if(current_file_name.find("WJetsToLNu_ext1_amcatnloFXFX") != std::string::npos 
     ////    || current_file_name.find("WJetsToLNu_amcatnloFXFX") != std::string::npos) {
     ////  if(GenW1_Pt > 120) passGenWZPt = false; // if W Pt > 120 GeV, cut it out
     ////}
     //if(current_file_name.find("DYJetsToLL_M-50_amcatnloFXFX") != std::string::npos) {
     //  if(GenZGamma1_Pt > 70) passGenWZPt = false; // if Z/gamma Pt > 70 GeV, cut it out
     //}
     //// first pt bin
     ////if(current_file_name.find("WJetsToLNu_Pt-100") != std::string::npos) {
     ////  if(GenW1_Pt <= 120) passGenWZPt = false;
     ////}
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
     fillVariableWithValue("PassGlobalTightHalo2016Filter" , int(readerTools_->ReadValueBranch<Double_t>("PassGlobalTightHalo2016Filter")          == 1));
     fillVariableWithValue("PassGoodVertices"              , int(readerTools_->ReadValueBranch<Double_t>("PassGoodVertices")                       == 1));
     fillVariableWithValue("PassHBHENoiseFilter"           , int(readerTools_->ReadValueBranch<Double_t>("PassHBHENoiseFilter")                    == 1));
     fillVariableWithValue("PassHBHENoiseIsoFilter"        , int(readerTools_->ReadValueBranch<Double_t>("PassHBHENoiseIsoFilter")                 == 1));
     fillVariableWithValue("PassBadEESupercrystalFilter"   , int(readerTools_->ReadValueBranch<Double_t>("PassBadEESupercrystalFilter")            == 1));
     fillVariableWithValue("PassEcalDeadCellTrigPrim"      , int(readerTools_->ReadValueBranch<Double_t>("PassEcalDeadCellTrigPrim")               == 1));
     fillVariableWithValue("PassChargedCandidateFilter"    , int(readerTools_->ReadValueBranch<Double_t>("PassChargedCandidateFilter")             == 1));
     fillVariableWithValue("PassBadPFMuonFilter"           , int(readerTools_->ReadValueBranch<Double_t>("PassBadPFMuonFilter")                    == 1));

     //--------------------------------------------------------------------------
     // Fill HLT
     //--------------------------------------------------------------------------

     //int passHLT = 1;
     //if ( isData ) { 
     //  passHLT = 0;
     //  //if ( H_Ele27_WPTight == 1 || H_Photon175 == 1)
     //  //  passHLT = 1;
     //  if (H_Ele27_WPTight_eta2p1 == 1)
     //    passHLT = 1;
     //  ////XXX SIC FIXME TEST
     //  //if(run < 275676) // L1 EGM efficiency going to zero
     //  //  passHLT = 0;
     //}
     //// FIXME: Update to new 2016 curve
     //else // using the turn-on in the MC
     //{
     //  // a la Z', throw a random number and if it's below the efficiency at this pt/eta, pass the event
     //  //   we get two chances to pass since we may have two electrons in the event
     //  passHLT = trigEle27::passTrig(Ele1_PtHeep,Ele1_SCEta) ? 1 : 0;
     //  //passHLT = trigEle27::passTrig(Ele1_Pt,Ele1_Eta) ? 1 : 0;
     //  if(!passHLT) // if the first one doesn't pass, try the second one
     //    passHLT = trigEle27::passTrig(Ele2_PtHeep,Ele2_SCEta) ? 1 : 0;
     //    //passHLT = trigEle27::passTrig(Ele2_Pt,Ele2_Eta) ? 1 : 0;
     //}
     //////XXX SIC FIXME TEST
     ////if (isData) {
     ////  passHLT = 0;
     ////  //if ( H_Ele27_WPLoose == 1)
     ////  //if ( H_Ele27_WPTight == 1)
     ////  //if ( H_Ele27_WPTight == 1 || H_Photon175 == 1)
     ////  if ( H_Ele27_WPLoose_eta2p1 == 1)
     ////    passHLT = 1;
     ////  if(run < 273726) // bad endcap alignment affecting deltaEtaIn cut
     ////    passHLT = 0;
     ////  if(run < 275676) // L1 EGM efficiency going to zero
     ////    passHLT = 0;
     ////}
     ////else {
     ////  passHLT = trigEle27::passTrig(Ele1_PtHeep,Ele1_SCEta) ? 1 : 0;
     ////  if(!passHLT) // if the first one doesn't pass, try the second one
     ////    passHLT = trigEle27::passTrig(Ele2_PtHeep,Ele2_SCEta) ? 1 : 0;
     ////}

     int passHLT = 0;
     bool verboseTrigEff = false;
     if ( isData() ) { 
       if(current_file_name.find("SinglePhoton") != std::string::npos) {
         if (readerTools_->ReadValueBranch<Double_t>("H_Photon175") == 1) // take events triggered by Photon175 only plus those triggered by Photon175 AND Ele27/Ele115
           passHLT = 1;
       }
       else if(current_file_name.find("SingleElectron") != std::string::npos) {
         if (readerTools_->ReadValueBranch<Double_t>("H_Photon175") != 1 && 
             (readerTools_->ReadValueBranch<Double_t>("H_Ele27_WPTight == 1") || readerTools_->ReadValueBranch<Double_t>("H_Ele115_CIdVT_GsfIdT") == 1) ) // take events triggered only by Ele27 OR Ele115
           passHLT = 1;
       }
       //if (H_Ele27_WPTight == 1 || H_Ele115_CIdVT_GsfIdT == 1)
       //  passHLT = 1;
     }
     else // using the turn-on in the MC
     {
       // a la Z', throw a random number and if it's below the efficiency at this pt/eta, pass the event
       //   we get two chances to pass since we may have two electrons in the event
       // trigger efficiency is binned in SCEta and SCEt (uncorrected)
       passHLT = triggerEfficiency.PassTrigger(readerTools_->ReadValueBranch<Double_t>("Ele1_SCEta"),readerTools_->ReadValueBranch<Double_t>("Ele1_SCEnergy")/cosh(readerTools_->ReadValueBranch<Double_t>("Ele1_SCEta")),verboseTrigEff) ? 1 : 0;
       if(!passHLT) // if the first one doesn't pass, try the second one
         passHLT = triggerEfficiency.PassTrigger(readerTools_->ReadValueBranch<Double_t>("Ele2_SCEta"),readerTools_->ReadValueBranch<Double_t>("Ele2_SCEnergy")/cosh(readerTools_->ReadValueBranch<Double_t>("Ele2_SCEta")),verboseTrigEff) ? 1 : 0;
     }
     fillVariableWithValue ( "PassHLT", passHLT, gen_weight * pileup_weight  ) ;     

     //--------------------------------------------------------------------------
     // Calculate variables for trigger matching 
     //--------------------------------------------------------------------------

     int nEle_hltMatched = 0.0;
     //FIXME reduced skim
     //if ( Ele1_hltEleSignalPt > 0.0 ) nEle_hltMatched++;
     //if ( Ele2_hltEleSignalPt > 0.0 ) nEle_hltMatched++;

     int nJet_hltMatched = 0.0;
     if ( readerTools_->ReadValueBranch<Double_t>("Jet1_hltJetPt") > 0.0 ) nJet_hltMatched++;
     if ( readerTools_->ReadValueBranch<Double_t>("Jet2_hltJetPt") > 0.0 ) nJet_hltMatched++;

     fillVariableWithValue("nEle_hltMatched",nEle_hltMatched, gen_weight * pileup_weight  );
     fillVariableWithValue("nJet_hltMatched",nJet_hltMatched, gen_weight * pileup_weight  );

     //--------------------------------------------------------------------------
     // Pass number of muons & electrons 
     // --> Special consideration if ttbar is derived from data
     //--------------------------------------------------------------------------

     // Muons and electrons
     bool is_ttbar_from_data = false;
     if ( readerTools_->ReadValueBranch<Double_t>("Ele2_ValidFrac") > 998. ) is_ttbar_from_data = true;
     // NB: we are doing data-driven ttbar separately so this shouldn't be the case here

     int PassNEle = 0;
     double nEle_ptCut = readerTools_->ReadValueBranch<Double_t>("nEle_ptCut");
     // nEle_ptCut are HEEP ID'ed electrons passing the Pt cut in the skim
     //if ( !is_ttbar_from_data && nEle_ptCut == 2 ) PassNEle = 1;
     //if (  is_ttbar_from_data && nEle_ptCut == 2 ) PassNEle = 1;
     if ( !is_ttbar_from_data && nEle_ptCut >= 2 )
       PassNEle = 1;
     if (  is_ttbar_from_data && nEle_ptCut >= 2 )
       PassNEle = 1;

     int PassNMuon = 0;
     double nMuon_ptCut = readerTools_->ReadValueBranch<Double_t>("nMuon_ptCut");
     if ( !is_ttbar_from_data && nMuon_ptCut == 0 )
       PassNMuon = 1;
     if (  is_ttbar_from_data && nMuon_ptCut >  0 )
       PassNMuon = 1;
     // now check if there are no muons but gen muons passing pt/eta which should have passed ID but did not
     // correct by a factor of (1-effData)/(1-effMC) to evaluate uncertainty
     if(PassNMuon) // no ID'ed muons in event
     {
       double GenMu1_Pt =  readerTools_->ReadValueBranch<Double_t>("GenMu1_Pt");
       double GenMu1_Eta = readerTools_->ReadValueBranch<Double_t>("GenMu1_Eta");
       double GenMu2_Pt =  readerTools_->ReadValueBranch<Double_t>("GenMu2_Pt");
       double GenMu2_Eta = readerTools_->ReadValueBranch<Double_t>("GenMu2_Eta");
       double GenMu3_Pt =  readerTools_->ReadValueBranch<Double_t>("GenMu3_Pt");
       double GenMu3_Eta = readerTools_->ReadValueBranch<Double_t>("GenMu3_Eta");
       // check for gen muons which should have passed ID but did not
       // pt threshold is defined in highPt muon ID in MuonIDs
       if(GenMu1_Pt > 35.0 && fabs(GenMu1_Eta) <= 2.4) {
         gen_weight*=MuonScaleFactors::GetVetoMCDataEffRatio(GenMu1_Eta);
         //cout << "Found gen muon 1 with Pt: " << GenMu1_Pt << " and eta: " << GenMu1_Eta << endl;
         //cout << "\tCorrect by: " << MuonScaleFactors::GetVetoMCDataEffRatio(GenMu1_Eta) << endl;
       }
       if(GenMu2_Pt > 35.0 && fabs(GenMu2_Eta) <= 2.4) {
         gen_weight*=MuonScaleFactors::GetVetoMCDataEffRatio(GenMu2_Eta);
         //cout << "Found gen muon 2 with Pt: " << GenMu2_Pt << " and eta: " << GenMu2_Eta << endl;
         //cout << "\tCorrect by: " << MuonScaleFactors::GetVetoMCDataEffRatio(GenMu2_Eta) << endl;
       }
       if(GenMu3_Pt > 35.0 && fabs(GenMu3_Eta) <= 2.4) {
         gen_weight*=MuonScaleFactors::GetVetoMCDataEffRatio(GenMu3_Eta);
         //cout << "Found gen muon 3 with Pt: " << GenMu3_Pt << " and eta: " << GenMu3_Eta << endl;
         //cout << "\tCorrect by: " << MuonScaleFactors::GetVetoMCDataEffRatio(GenMu3_Eta) << endl;
       }
     }


     fillVariableWithValue("PassNEle" , PassNEle , gen_weight * pileup_weight);
     fillVariableWithValue("PassNMuon", PassNMuon, gen_weight * pileup_weight);

     //--------------------------------------------------------------------------
     // Calculate electron-jet pair mass values
     //--------------------------------------------------------------------------

     double M_ej_avg, M_ej_min, M_ej_max;

     double nEle_store = readerTools_->ReadValueBranch<Double_t>("nEle_store");
     double nJet_store = readerTools_->ReadValueBranch<Double_t>("nJet_store");
     double M_e1j1 = readerTools_->ReadValueBranch<Double_t>("M_e1j1");
     double M_e1j2 = readerTools_->ReadValueBranch<Double_t>("M_e1j2");
     double M_e2j1 = readerTools_->ReadValueBranch<Double_t>("M_e2j1");
     double M_e2j2 = readerTools_->ReadValueBranch<Double_t>("M_e2j2");
     double Pt_e1e2 = readerTools_->ReadValueBranch<Double_t>("Pt_e1e2");
     if ( nEle_store >= 2 && nJet_store >= 2) {
       if ( fabs(M_e1j1-M_e2j2) < fabs(M_e1j2-M_e2j1) )  {
         M_ej_avg = (M_e1j1 + M_e2j2) / 2.0;
         if    ( M_e1j1 < M_e2j2 ) { M_ej_min = M_e1j1; M_ej_max = M_e2j2; }
         else                      { M_ej_min = M_e2j2; M_ej_max = M_e1j1; }
       }
       else { 
         M_ej_avg = (M_e1j2 + M_e2j1) / 2.0;
         if    ( M_e1j2 < M_e2j1 ) { M_ej_min = M_e1j2; M_ej_max = M_e2j1; }
         else                      { M_ej_min = M_e2j1; M_ej_max = M_e1j2; }
       }
     }

     double sT_zjj = Pt_e1e2 + readerTools_->ReadValueBranch<Double_t>("Jet1_Pt") + readerTools_->ReadValueBranch<Double_t>("Jet2_Pt");

     //--------------------------------------------------------------------------
     // Fill electron variables 
     //--------------------------------------------------------------------------
     double Ele1_PtHeep = readerTools_->ReadValueBranch<Double_t>("Ele1_PtHeep");
     double Ele2_PtHeep = readerTools_->ReadValueBranch<Double_t>("Ele2_PtHeep");
     double Ele1_Eta = readerTools_->ReadValueBranch<Double_t>("Ele1_Eta");
     double Ele2_Eta = readerTools_->ReadValueBranch<Double_t>("Ele2_Eta");
     double Ele1_TrkEta = readerTools_->ReadValueBranch<Double_t>("Ele1_TrkEta");
     double Ele2_TrkEta = readerTools_->ReadValueBranch<Double_t>("Ele2_TrkEta");
     
     if ( nEle_store >= 1 ) {
       fillVariableWithValue( "Ele1_PtHeep",            Ele1_PtHeep, gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "Ele1_AbsDeltaEtaEleTrk",
           fabs(Ele1_Eta-Ele1_TrkEta), gen_weight * pileup_weight );
     }
     if ( nEle_store >= 2 ) {
       fillVariableWithValue( "Ele2_PtHeep",            Ele2_PtHeep, gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "Ele2_AbsDeltaEtaEleTrk",
           fabs(Ele2_Eta-Ele2_TrkEta), gen_weight * pileup_weight );
     }
			
     //--------------------------------------------------------------------------
     // Fill jet variables 
     //--------------------------------------------------------------------------

     double nJet_ptCut = readerTools_->ReadValueBranch<Double_t>("nJet_ptCut");
     double DR_Jet1Jet2 = readerTools_->ReadValueBranch<Double_t>("DR_Jet1Jet2");
     // Jets								    
     fillVariableWithValue("nJet", nJet_ptCut, gen_weight * pileup_weight );
     if ( nJet_store >= 1 ) { 						    
       fillVariableWithValue( "Jet1_Pt"    , readerTools_->ReadValueBranch<Double_t>("Jet1_Pt")     , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "Jet1_Eta"   , readerTools_->ReadValueBranch<Double_t>("Jet1_Eta")    , gen_weight * pileup_weight  ) ;
     }
     if ( nJet_store >= 2 ) { 
       fillVariableWithValue( "Jet2_Pt"    , readerTools_->ReadValueBranch<Double_t>("Jet2_Pt")     , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "Jet2_Eta"   , readerTools_->ReadValueBranch<Double_t>("Jet2_Eta")    , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "DR_Jet1Jet2", DR_Jet1Jet2 , gen_weight * pileup_weight  ) ;
     }

     //--------------------------------------------------------------------------
     // Fill DeltaR variables
     //--------------------------------------------------------------------------

     double DR_Ele1Jet1 = readerTools_->ReadValueBranch<Double_t>("DR_Ele1Jet1");
     double DR_Ele2Jet1 = readerTools_->ReadValueBranch<Double_t>("DR_Ele2Jet1");
     double DR_Ele1Jet2 = readerTools_->ReadValueBranch<Double_t>("DR_Ele1Jet2");
     double DR_Ele2Jet2 = readerTools_->ReadValueBranch<Double_t>("DR_Ele2Jet2");
     if ( nEle_store >= 2 && nJet_store >= 1) {
       fillVariableWithValue( "DR_Ele1Jet1"  , DR_Ele1Jet1, gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "DR_Ele2Jet1"  , DR_Ele2Jet1 , gen_weight * pileup_weight  ) ;
       if(nJet_store >= 2) {
         fillVariableWithValue( "DR_Ele1Jet2", DR_Ele1Jet2 , gen_weight * pileup_weight  ) ;
         fillVariableWithValue( "DR_Ele2Jet2", DR_Ele2Jet2 , gen_weight * pileup_weight  ) ;
       }
     }


     //--------------------------------------------------------------------------
     // Multi-object variables
     //--------------------------------------------------------------------------
     double sT_eejj = readerTools_->ReadValueBranch<Double_t>("sT_eejj");
     double M_e1e2 =  readerTools_->ReadValueBranch<Double_t>("M_e1e2");

     if ( nEle_store >= 2 ) { 						    
       fillVariableWithValue( "M_e1e2"     , M_e1e2 , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "M_e1e2_opt" , M_e1e2 , gen_weight * pileup_weight  ) ;
       // for Ptee cut
       fillVariableWithValue( "Pt_e1e2"    , Pt_e1e2 , gen_weight * pileup_weight  ) ;

       if ( nJet_store >= 2 ) { 
         // SIC recompute sT using PtHeep. FIXME: this is now being done in skims
         //sT_eejj = Ele1_PtHeep+Ele2_PtHeep+Jet1_Pt+Jet2_Pt;
         fillVariableWithValue( "sT_eejj"    , sT_eejj , gen_weight * pileup_weight  ) ;
         fillVariableWithValue( "sT_eejj_opt", sT_eejj , gen_weight * pileup_weight  ) ;
         fillVariableWithValue( "Mej_min_opt", M_ej_min, gen_weight * pileup_weight  ) ;
       }      
     }

     //--------------------------------------------------------------------------
     // Fill final selection cuts
     //--------------------------------------------------------------------------

     if(doFinalSelections)
     {
       for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
         int lq_mass = LQ_MASS[i_lq_mass];
         sprintf(cut_name, "M_e1e2_LQ%d"  , lq_mass );
         fillVariableWithValue ( cut_name, M_e1e2  , gen_weight * pileup_weight  ) ;
         sprintf(cut_name, "sT_eejj_LQ%d" , lq_mass );
         fillVariableWithValue ( cut_name, sT_eejj , gen_weight * pileup_weight  ) ;
         sprintf(cut_name, "min_M_ej_LQ%d", lq_mass );
         fillVariableWithValue ( cut_name, M_ej_min, gen_weight * pileup_weight  ) ;
       }
     }
     
     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //std::cout << "getVariableValue(\"PassNEle\")=" << getVariableValue("PassNEle") << "; passedCut(\"PassNEle\")=" << passedCut("PassNEle") << std::endl;
     //--------------------------------------------------------------------------
     // Did we at least pass the noise filtes?
     //--------------------------------------------------------------------------
     
     bool passed_minimum = ( passedAllPreviousCuts("PassBadPFMuonFilter") && passedCut ("PassBadPFMuonFilter"));
     
     //--------------------------------------------------------------------------
     // Did we pass preselection?
     //--------------------------------------------------------------------------

     bool passed_preselection = ( passedAllPreviousCuts("M_e1e2") && passedCut ("M_e1e2") );

     //--------------------------------------------------------------------------
     // Are we in the region of interest?
     //--------------------------------------------------------------------------

     bool passed_region_of_interest = bool (
         passed_preselection && M_e1e2 > 170. && sT_eejj > 900.0
         );

     //--------------------------------------------------------------------------
     // Did we pass any final selection cuts?
     //--------------------------------------------------------------------------

     passed_vector.clear();
     if(doFinalSelections)
     {
       for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
         int lq_mass = LQ_MASS[i_lq_mass];
         //sprintf(cut_name, "M_e1e2_LQ%d", lq_mass );
         sprintf(cut_name, "min_M_ej_LQ%d", lq_mass ); // this is actually the last cut in the cut file...!
         bool decision = bool ( passedAllPreviousCuts(cut_name) && passedCut (cut_name));
         passed_vector.push_back (decision);
       }
     }

     //--------------------------------------------------------------------------
     // Fill plots with no selection applied
     //--------------------------------------------------------------------------

     FillUserTH1D( "PileupWeight"   , pileup_weight );
     FillUserTH1D( "GeneratorWeight", gen_weight ) ;

     //--------------------------------------------------------------------------
     // Fill noise filter level plots
     //--------------------------------------------------------------------------

     double nVertex = readerTools_->ReadValueBranch<Double_t>("nVertex");
     if ( passed_minimum && isData() ){ 
       FillUserTH1D ("run_HLT", run );
       profile_run_vs_nvtx_HLT -> Fill (run, nVertex, 1);
     }

     //--------------------------------------------------------------------------
     // Print if desired
     //--------------------------------------------------------------------------

     /*
     if ( isData ) {
       std::cout.precision(0);
       std::cout << fixed <<  "Run = " << run << ", event = " << event << ", ls = " << ls << std::endl;
       std::cout.precision(3);
       std::cout << fixed <<  "  Mej      = " << M_ej_avg << std::endl;
       std::cout << fixed <<  "  Mee      = " << M_e1e2 << std::endl;
       std::cout << fixed <<  "  sT       = " << sT_enujj << std::endl;
       std::cout << fixed <<  "  Ele1 Pt  = " << Ele1_Pt << "\t, Eta = " << Ele1_Eta << "\t, Phi = " << Ele1_Phi << std::endl;
       std::cout << fixed <<  "  Ele2 Pt  = " << Ele2_Pt << "\t, Eta = " << Ele2_Eta << "\t, Phi = " << Ele2_Phi << std::endl;
       std::cout << fixed <<  "  Jet1 Pt  = " << Jet1_Pt << "\t, Eta = " << Jet1_Eta << "\t, Phi = " << Jet1_Phi << std::endl;
       std::cout << fixed <<  "  Jet2 Pt  = " << Jet2_Pt << "\t, Eta = " << Jet2_Eta << "\t, Phi = " << Jet2_Phi << std::endl;
     }
     */

     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------

     if ( passed_preselection ) {

       //--------------------------------------------------------------------------
       // Recalculate some variables
       //--------------------------------------------------------------------------

       TLorentzVector e1, j1, e2, j2,j3, mu, met;
       TLorentzVector eejj, e1e2mu;
       TLorentzVector eej, ejj, ee;
       TLorentzVector e1j3, e2j3, j1j3, j2j3, j1j2, j1j2j3, eejjj;
       double Ele1_Pt = readerTools_->ReadValueBranch<Double_t>("Ele1_Pt");
       double Ele1_Eta = readerTools_->ReadValueBranch<Double_t>("Ele1_Eta");
       double Ele1_Phi = readerTools_->ReadValueBranch<Double_t>("Ele1_Phi");
       double Ele2_Pt = readerTools_->ReadValueBranch<Double_t>("Ele2_Pt");
       double Ele2_Eta = readerTools_->ReadValueBranch<Double_t>("Ele2_Eta");
       double Ele2_Phi = readerTools_->ReadValueBranch<Double_t>("Ele2_Phi");
       double Jet1_Pt = readerTools_->ReadValueBranch<Double_t>("Jet1_Pt");
       double Jet1_Eta = readerTools_->ReadValueBranch<Double_t>("Jet1_Eta");
       double Jet1_Phi = readerTools_->ReadValueBranch<Double_t>("Jet1_Phi");
       double Jet2_Pt = readerTools_->ReadValueBranch<Double_t>("Jet2_Pt");
       double Jet2_Phi = readerTools_->ReadValueBranch<Double_t>("Jet2_Phi");
       double Jet2_Eta = readerTools_->ReadValueBranch<Double_t>("Jet2_Eta");
       double Muon1_Pt = readerTools_->ReadValueBranch<Double_t>("Muon1_Pt");
       double Muon1_Eta = readerTools_->ReadValueBranch<Double_t>("Muon1_Eta");
       double Muon1_Phi = readerTools_->ReadValueBranch<Double_t>("Muon1_Phi");
       double Muon2_Pt = readerTools_->ReadValueBranch<Double_t>("Muon2_Pt");
       double Muon2_Eta = readerTools_->ReadValueBranch<Double_t>("Muon2_Eta");
       double Muon2_Phi = readerTools_->ReadValueBranch<Double_t>("Muon2_Phi");
       double PFMET_Type1_Pt = readerTools_->ReadValueBranch<Double_t>("PFMET_Type1_Pt");
       double PFMET_Type1_Phi = readerTools_->ReadValueBranch<Double_t>("PFMET_Type1_Phi");

       e1.SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
       e2.SetPtEtaPhiM ( Ele2_Pt, Ele2_Eta, Ele2_Phi, 0.0 );
       j1.SetPtEtaPhiM ( Jet1_Pt, Jet1_Eta, Jet1_Phi, 0.0 );
       j2.SetPtEtaPhiM ( Jet2_Pt, Jet2_Eta, Jet2_Phi, 0.0 );
       mu.SetPtEtaPhiM ( Muon1_Pt,Muon1_Eta,Muon1_Phi, 0.0 );
       met.SetPtEtaPhiM ( PFMET_Type1_Pt, 0.0, PFMET_Type1_Phi, 0.0 );

       eejj = e1 + e2 + j1 + j2 ; 
       eej  = e1 + e2 + j1;
       ejj  = e1 + j1 + j2;
       ee   = e1 + e2;
       j1j2 = j1 + j2;

       double min_DR_EleJet;
       double DR_Ele1Jet3;
       double DR_Ele2Jet3;
       double DR_Ele1Ele2 = e1.DeltaR( e2 );

       double M_eejj = eejj.M();
       double M_eej  = eej.M();
       double M_ejj  = ejj.M();

       double Pt_j1j2 = j1j2.Pt();
       double Pt_j1j3;
       double Pt_j2j3;
       double Pt_j1j2j3;

       double M_j1j3, M_j2j3;
       double M_e1j3, M_e2j3;
       double M_eejjj;

       double min_DeltaR_Zj = 999.;
       if ( ee.DeltaR ( j1 ) < min_DeltaR_Zj ) min_DeltaR_Zj = ee.DeltaR ( j1 );
       if ( ee.DeltaR ( j2 ) < min_DeltaR_Zj ) min_DeltaR_Zj = ee.DeltaR ( j2 );
       double DR_ZJ1 = ee.DeltaR ( j1 );
       double DR_ZJ2 = ee.DeltaR ( j2 );

       if ( nJet_ptCut > 2 ) { 
         j3.SetPtEtaPhiM ( readerTools_->ReadValueBranch<Double_t>("Jet3_Pt"), readerTools_->ReadValueBranch<Double_t>("Jet3_Eta"), readerTools_->ReadValueBranch<Double_t>("Jet3_Phi"), 0.0 );

         e1j3 = e1 + j3;
         e2j3 = e2 + j3;
         j1j3 = j1 + j3;
         j2j3 = j2 + j3;
         j1j2j3 = j1 + j2 + j3;
         eejjj = e1 + e2 + j1 + j2 + j3;

         DR_Ele1Jet3 = e1.DeltaR( j3 );
         DR_Ele2Jet3 = e2.DeltaR( j3 );
         M_e1j3 = e1j3.M();
         M_e2j3 = e2j3.M();
         M_j1j3 = j1j3.M();
         M_j2j3 = j2j3.M();
         Pt_j1j3 = j1j3.Pt();
         Pt_j2j3 = j2j3.Pt();
         Pt_j1j2j3 = j1j2j3.Pt();
         M_eejjj = eejjj.M();
       }

       if ( DR_Ele1Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet1;
       if ( DR_Ele1Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet2;
       if ( DR_Ele2Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet1;
       if ( DR_Ele2Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet2;
       if (nJet_ptCut > 2 ) {
         if ( DR_Ele1Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet3;
         if ( DR_Ele2Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet3;
       }

       //--------------------------------------------------------------------------
       // Determine what kind of event this is:
       // - EB-EB
       // - EB-EE
       // - EE-EE
       //--------------------------------------------------------------------------

       bool isEB1 = ( fabs(Ele1_Eta) < eleEta_bar_max  ) ;
       bool isEE1 = ( fabs(Ele1_Eta) > eleEta_end_min &&
           fabs(Ele1_Eta) < eleEta_end_max ) ;

       bool isEB2 = ( fabs(Ele2_Eta) < eleEta_bar_max  ) ;
       bool isEE2 = ( fabs(Ele2_Eta) > eleEta_end_min &&
           fabs(Ele2_Eta) < eleEta_end_max ) ;

       bool isEBEB = ( isEB1 && isEB2 ) ;
       bool isEBEE = ( ( isEB1 && isEE2 ) ||
           ( isEE1 && isEB2 ) );
       bool isEEEE = ( isEE1  && isEE2  );
       bool isEB   = ( isEBEB || isEBEE );


       //--------------------------------------------------------------------------
       // Electron quality histograms (preselection)
       //--------------------------------------------------------------------------
       double Ele1_BeamSpotDXY          = readerTools_->ReadValueBranch<Double_t>("Ele1_BeamSpotDXY")        ; 
       double Ele1_Classif              = readerTools_->ReadValueBranch<Double_t>("Ele1_Classif")            ; 
       double Ele1_CorrIsolation        = readerTools_->ReadValueBranch<Double_t>("Ele1_CorrIsolation")      ; 
       double Ele1_DeltaEtaTrkSC        = readerTools_->ReadValueBranch<Double_t>("Ele1_DeltaEtaTrkSC")      ; 
       double Ele1_DeltaPhiTrkSC        = readerTools_->ReadValueBranch<Double_t>("Ele1_DeltaPhiTrkSC")      ; 
       double Ele1_Full5x5E1x5OverE5x5  = readerTools_->ReadValueBranch<Double_t>("Ele1_Full5x5E1x5OverE5x5"); 
       double Ele1_Full5x5E2x5OverE5x5  = readerTools_->ReadValueBranch<Double_t>("Ele1_Full5x5E2x5OverE5x5"); 
       double Ele1_EcalIsolation        = readerTools_->ReadValueBranch<Double_t>("Ele1_EcalIsolation")      ; 
       double Ele1_HcalIsolation        = readerTools_->ReadValueBranch<Double_t>("Ele1_HcalIsolation")      ; 
       double Ele1_TrkIsolation         = readerTools_->ReadValueBranch<Double_t>("Ele1_TrkIsolation")       ; 
       double Ele1_Energy               = readerTools_->ReadValueBranch<Double_t>("Ele1_Energy")             ; 
       double Ele1_FBrem                = readerTools_->ReadValueBranch<Double_t>("Ele1_FBrem")              ; 
       double Ele1_GsfCtfCharge         = readerTools_->ReadValueBranch<Double_t>("Ele1_GsfCtfCharge")       ; 
       double Ele1_GsfCtfScPixCharge    = readerTools_->ReadValueBranch<Double_t>("Ele1_GsfCtfScPixCharge")  ; 
       double Ele1_GsfScPixCharge       = readerTools_->ReadValueBranch<Double_t>("Ele1_GsfScPixCharge")     ; 
       double Ele1_HasMatchedPhot       = readerTools_->ReadValueBranch<Double_t>("Ele1_HasMatchedPhot")     ; 
       double Ele1_HoE                  = readerTools_->ReadValueBranch<Double_t>("Ele1_HoE")                ; 
       double Ele1_LeadVtxDistXY        = readerTools_->ReadValueBranch<Double_t>("Ele1_LeadVtxDistXY")      ; 
       double Ele1_LeadVtxDistZ         = readerTools_->ReadValueBranch<Double_t>("Ele1_LeadVtxDistZ")       ; 
       double Ele1_MissingHits          = readerTools_->ReadValueBranch<Double_t>("Ele1_MissingHits")        ; 
       double Ele1_NBrems               = readerTools_->ReadValueBranch<Double_t>("Ele1_NBrems")             ; 
       double Ele1_ValidFrac            = readerTools_->ReadValueBranch<Double_t>("Ele1_ValidFrac")          ; 
       double Ele1_RawEnergy            = readerTools_->ReadValueBranch<Double_t>("Ele1_RawEnergy");
       double Ele1_SCEta            = readerTools_->ReadValueBranch<Double_t>("Ele1_SCEta");
       double Ele1_TrkPt                = readerTools_->ReadValueBranch<Double_t>("Ele1_TrkPt"); 
       double Ele1_SigmaEtaEta          = readerTools_->ReadValueBranch<Double_t>("Ele1_SigmaEtaEta");
       double Ele1_Full5x5SigmaIEtaIEta = readerTools_->ReadValueBranch<Double_t>("Ele1_Full5x5SigmaIEtaIEta");
       double Ele1_Charge               = readerTools_->ReadValueBranch<Double_t>("Ele1_Charge");
       double Ele1_PtHeep               = readerTools_->ReadValueBranch<Double_t>("Ele1_PtHeep");

       FillUserTH1D("BeamSpotDXY_1stEle_PAS"         ,Ele1_BeamSpotDXY         , pileup_weight * gen_weight    ); 
       FillUserTH1D("Classif_1stEle_PAS"             ,Ele1_Classif             , pileup_weight * gen_weight    ); 
       FillUserTH1D("CorrIsolation_1stEle_PAS"       ,Ele1_CorrIsolation       , pileup_weight * gen_weight    ); 
       FillUserTH1D("DeltaEtaTrkSC_1stEle_PAS"       ,Ele1_DeltaEtaTrkSC       , pileup_weight * gen_weight    ); 
       FillUserTH1D("DeltaPhiTrkSC_1stEle_PAS"       ,Ele1_DeltaPhiTrkSC       , pileup_weight * gen_weight    ); 
       FillUserTH1D("Full5x5E1x5OverE5x5_1stEle_PAS" ,Ele1_Full5x5E1x5OverE5x5 , pileup_weight * gen_weight    ); 
       FillUserTH1D("Full5x5E2x5OverE5x5_1stEle_PAS" ,Ele1_Full5x5E2x5OverE5x5 , pileup_weight * gen_weight    ); 
       FillUserTH1D("EcalIsolation_1stEle_PAS"       ,Ele1_EcalIsolation       , pileup_weight * gen_weight    ); 
       FillUserTH1D("HcalIsolation_1stEle_PAS"       ,Ele1_HcalIsolation       , pileup_weight * gen_weight    ); 
       FillUserTH1D("TrkIsolation_1stEle_PAS"        ,Ele1_TrkIsolation        , pileup_weight * gen_weight    ); 
       FillUserTH1D("Energy_1stEle_PAS"              ,Ele1_Energy              , pileup_weight * gen_weight    ); 
       FillUserTH1D("FBrem_1stEle_PAS"               ,Ele1_FBrem               , pileup_weight * gen_weight    ); 
       FillUserTH1D("GsfCtfCharge_1stEle_PAS"        ,Ele1_GsfCtfCharge        , pileup_weight * gen_weight    ); 
       FillUserTH1D("GsfCtfScPixCharge_1stEle_PAS"   ,Ele1_GsfCtfScPixCharge   , pileup_weight * gen_weight    ); 
       FillUserTH1D("GsfScPixCharge_1stEle_PAS"      ,Ele1_GsfScPixCharge      , pileup_weight * gen_weight    ); 
       FillUserTH1D("HasMatchedPhot_1stEle_PAS"      ,Ele1_HasMatchedPhot      , pileup_weight * gen_weight    ); 
       FillUserTH1D("HoE_1stEle_PAS"                 ,Ele1_HoE                 , pileup_weight * gen_weight    ); 
       FillUserTH1D("LeadVtxDistXY_1stEle_PAS"       ,Ele1_LeadVtxDistXY       , pileup_weight * gen_weight    ); 
       FillUserTH1D("LeadVtxDistZ_1stEle_PAS"        ,Ele1_LeadVtxDistZ        , pileup_weight * gen_weight    ); 
       FillUserTH1D("MissingHits_1stEle_PAS"         ,Ele1_MissingHits         , pileup_weight * gen_weight    ); 
       FillUserTH1D("NBrems_1stEle_PAS"              ,Ele1_NBrems              , pileup_weight * gen_weight    ); 
       FillUserTH1D("ValidFrac_1stEle_PAS"           ,Ele1_ValidFrac           , pileup_weight * gen_weight    ); 
       FillUserTH1D("EnergyORawEnergy_1stEle_PAS"    ,Ele1_Energy / Ele1_RawEnergy, pileup_weight * gen_weight    ); 
       FillUserTH1D("TrkPtOPt_1stEle_PAS"            ,Ele1_TrkPt  / Ele1_Pt   , pileup_weight * gen_weight    ); 
       if ( fabs(Ele1_Eta) < eleEta_bar ) { 
         FillUserTH1D("SigmaEtaEta_Barrel_1stEle_PAS"  , Ele1_SigmaEtaEta                    , pileup_weight * gen_weight    ); 
         FillUserTH1D("Full5x5SigmaIEtaIEta_Barrel_1stEle_PAS", Ele1_Full5x5SigmaIEtaIEta                  , pileup_weight * gen_weight    ); 
       }
       else if ( fabs(Ele1_Eta) > eleEta_end1_min && fabs(Ele2_Eta) < eleEta_end2_max ){
         FillUserTH1D("SigmaEtaEta_Endcap_1stEle_PAS"  , Ele1_SigmaEtaEta                    , pileup_weight * gen_weight    ); 
         FillUserTH1D("Full5x5SigmaIEtaIEta_Endcap_1stEle_PAS", Ele1_Full5x5SigmaIEtaIEta                  , pileup_weight * gen_weight    ); 
       }

      double Ele2_BeamSpotDXY          = readerTools_->ReadValueBranch<Double_t>("Ele2_BeamSpotDXY")        ; 
      double Ele2_Classif              = readerTools_->ReadValueBranch<Double_t>("Ele2_Classif")            ; 
      double Ele2_CorrIsolation        = readerTools_->ReadValueBranch<Double_t>("Ele2_CorrIsolation")      ; 
      double Ele2_DeltaEtaTrkSC        = readerTools_->ReadValueBranch<Double_t>("Ele2_DeltaEtaTrkSC")      ; 
      double Ele2_DeltaPhiTrkSC        = readerTools_->ReadValueBranch<Double_t>("Ele2_DeltaPhiTrkSC")      ; 
      double Ele2_Full5x5E1x5OverE5x5  = readerTools_->ReadValueBranch<Double_t>("Ele2_Full5x5E1x5OverE5x5"); 
      double Ele2_Full5x5E2x5OverE5x5  = readerTools_->ReadValueBranch<Double_t>("Ele2_Full5x5E2x5OverE5x5"); 
      double Ele2_EcalIsolation        = readerTools_->ReadValueBranch<Double_t>("Ele2_EcalIsolation")      ; 
      double Ele2_HcalIsolation        = readerTools_->ReadValueBranch<Double_t>("Ele2_HcalIsolation")      ; 
      double Ele2_TrkIsolation         = readerTools_->ReadValueBranch<Double_t>("Ele2_TrkIsolation")       ; 
      double Ele2_Energy               = readerTools_->ReadValueBranch<Double_t>("Ele2_Energy")             ; 
      double Ele2_FBrem                = readerTools_->ReadValueBranch<Double_t>("Ele2_FBrem")              ; 
      double Ele2_GsfCtfCharge         = readerTools_->ReadValueBranch<Double_t>("Ele2_GsfCtfCharge")       ; 
      double Ele2_GsfCtfScPixCharge    = readerTools_->ReadValueBranch<Double_t>("Ele2_GsfCtfScPixCharge")  ; 
      double Ele2_GsfScPixCharge       = readerTools_->ReadValueBranch<Double_t>("Ele2_GsfScPixCharge")     ; 
      double Ele2_HasMatchedPhot       = readerTools_->ReadValueBranch<Double_t>("Ele2_HasMatchedPhot")     ; 
      double Ele2_HoE                  = readerTools_->ReadValueBranch<Double_t>("Ele2_HoE")                ; 
      double Ele2_LeadVtxDistXY        = readerTools_->ReadValueBranch<Double_t>("Ele2_LeadVtxDistXY")      ; 
      double Ele2_LeadVtxDistZ         = readerTools_->ReadValueBranch<Double_t>("Ele2_LeadVtxDistZ")       ; 
      double Ele2_MissingHits          = readerTools_->ReadValueBranch<Double_t>("Ele2_MissingHits")        ; 
      double Ele2_NBrems               = readerTools_->ReadValueBranch<Double_t>("Ele2_NBrems")             ; 
      double Ele2_ValidFrac            = readerTools_->ReadValueBranch<Double_t>("Ele2_ValidFrac")          ; 
      double Ele2_RawEnergy            = readerTools_->ReadValueBranch<Double_t>("Ele2_RawEnergy");
      double Ele2_SCEta            = readerTools_->ReadValueBranch<Double_t>("Ele2_SCEta");
      double Ele2_TrkPt                = readerTools_->ReadValueBranch<Double_t>("Ele2_TrkPt"); 
      double Ele2_SigmaEtaEta          = readerTools_->ReadValueBranch<Double_t>("Ele2_SigmaEtaEta");
      double Ele2_Full5x5SigmaIEtaIEta = readerTools_->ReadValueBranch<Double_t>("Ele2_Full5x5SigmaIEtaIEta");
      double Ele2_Charge               = readerTools_->ReadValueBranch<Double_t>("Ele2_Charge");
      double Ele2_PtHeep               = readerTools_->ReadValueBranch<Double_t>("Ele2_PtHeep");

       FillUserTH1D("BeamSpotDXY_2ndEle_PAS"         ,Ele2_BeamSpotDXY              , pileup_weight * gen_weight    ); 
       FillUserTH1D("Classif_2ndEle_PAS"             ,Ele2_Classif                  , pileup_weight * gen_weight    ); 
       FillUserTH1D("CorrIsolation_2ndEle_PAS"       ,Ele2_CorrIsolation            , pileup_weight * gen_weight    ); 
       FillUserTH1D("DeltaEtaTrkSC_2ndEle_PAS"       ,Ele2_DeltaEtaTrkSC            , pileup_weight * gen_weight    ); 
       FillUserTH1D("DeltaPhiTrkSC_2ndEle_PAS"       ,Ele2_DeltaPhiTrkSC            , pileup_weight * gen_weight    ); 
       FillUserTH1D("Full5x5E1x5OverE5x5_2ndEle_PAS" ,Ele2_Full5x5E1x5OverE5x5      , pileup_weight * gen_weight    ); 
       FillUserTH1D("Full5x5E2x5OverE5x5_2ndEle_PAS" ,Ele2_Full5x5E2x5OverE5x5      , pileup_weight * gen_weight    ); 
       FillUserTH1D("EcalIsolation_2ndEle_PAS"       ,Ele2_EcalIsolation            , pileup_weight * gen_weight    ); 
       FillUserTH1D("HcalIsolation_2ndEle_PAS"       ,Ele2_HcalIsolation            , pileup_weight * gen_weight    ); 
       FillUserTH1D("TrkIsolation_2ndEle_PAS"        ,Ele2_TrkIsolation             , pileup_weight * gen_weight    ); 
       FillUserTH1D("Energy_2ndEle_PAS"              ,Ele2_Energy                   , pileup_weight * gen_weight    ); 
       FillUserTH1D("FBrem_2ndEle_PAS"               ,Ele2_FBrem                    , pileup_weight * gen_weight    ); 
       FillUserTH1D("GsfCtfCharge_2ndEle_PAS"        ,Ele2_GsfCtfCharge             , pileup_weight * gen_weight    ); 
       FillUserTH1D("GsfCtfScPixCharge_2ndEle_PAS"   ,Ele2_GsfCtfScPixCharge        , pileup_weight * gen_weight    ); 
       FillUserTH1D("GsfScPixCharge_2ndEle_PAS"      ,Ele2_GsfScPixCharge           , pileup_weight * gen_weight    ); 
       FillUserTH1D("HasMatchedPhot_2ndEle_PAS"      ,Ele2_HasMatchedPhot           , pileup_weight * gen_weight    ); 
       FillUserTH1D("HoE_2ndEle_PAS"                 ,Ele2_HoE                      , pileup_weight * gen_weight    ); 
       FillUserTH1D("LeadVtxDistXY_2ndEle_PAS"       ,Ele2_LeadVtxDistXY            , pileup_weight * gen_weight    ); 
       FillUserTH1D("LeadVtxDistZ_2ndEle_PAS"        ,Ele2_LeadVtxDistZ             , pileup_weight * gen_weight    ); 
       FillUserTH1D("MissingHits_2ndEle_PAS"         ,Ele2_MissingHits              , pileup_weight * gen_weight    ); 
       FillUserTH1D("NBrems_2ndEle_PAS"              ,Ele2_NBrems                   , pileup_weight * gen_weight    ); 
       FillUserTH1D("ValidFrac_2ndEle_PAS"           ,Ele2_ValidFrac                , pileup_weight * gen_weight    ); 
       FillUserTH1D("EnergyORawEnergy_2ndEle_PAS"    ,Ele2_Energy / Ele2_RawEnergy , pileup_weight * gen_weight    ); 
       FillUserTH1D("TrkPtOPt_2ndEle_PAS"            ,Ele2_TrkPt  / Ele2_Pt        , pileup_weight * gen_weight    ); 
       if ( fabs(Ele2_Eta) < eleEta_bar ) { 
         FillUserTH1D("SigmaEtaEta_Barrel_2ndEle_PAS"  , Ele2_SigmaEtaEta                    , pileup_weight * gen_weight    ); 
         FillUserTH1D("Full5x5SigmaIEtaIEta_Barrel_2ndEle_PAS", Ele2_Full5x5SigmaIEtaIEta                  , pileup_weight * gen_weight    ); 
       }
       else if ( fabs(Ele2_Eta) > eleEta_end1_min && fabs(Ele2_Eta) < eleEta_end2_max ){
         FillUserTH1D("SigmaEtaEta_Endcap_2ndEle_PAS"  , Ele2_SigmaEtaEta                    , pileup_weight * gen_weight    ); 
         FillUserTH1D("Full5x5SigmaIEtaIEta_Endcap_2ndEle_PAS", Ele2_Full5x5SigmaIEtaIEta                  , pileup_weight * gen_weight    ); 
       }

       //--------------------------------------------------------------------------
       // Preselection histograms
       //--------------------------------------------------------------------------
       double M_j1j2 = readerTools_->ReadValueBranch<Double_t>("M_j1j2");

       FillUserTH1D( "Ptj1j2_PAS"           , Pt_j1j2                        , pileup_weight * gen_weight );
       FillUserTH1D( "Ptee_Minus_Ptj1j2_PAS", Pt_e1e2 - Pt_j1j2              , pileup_weight * gen_weight );
       //FillUserTH1D("ProcessID_PAS"         , ProcessID                      , pileup_weight * gen_weight );
       FillUserTH1D("minDR_EleJet_PAS"      , min_DR_EleJet                  , pileup_weight * gen_weight );
       FillUserTH1D("DR_Ele1Ele2_PAS"	    , DR_Ele1Ele2                    , pileup_weight * gen_weight );
       FillUserTH1D("EleChargeSum_PAS"      , Ele1_Charge + Ele2_Charge      , pileup_weight * gen_weight );
       FillUserTH1D("nElectron_PAS"         , nEle_ptCut                     , pileup_weight * gen_weight );
       FillUserTH1D("nMuon_PAS"             , nMuon_ptCut                    , pileup_weight * gen_weight );
       FillUserTH1D("nJet_PAS"              , nJet_ptCut                     , pileup_weight * gen_weight );
       FillUserTH1D("Pt1stEle_PAS"	        , Ele1_Pt                        , pileup_weight * gen_weight );
       FillUserTH1D("PtHeep1stEle_PAS"	    , Ele1_PtHeep                    , pileup_weight * gen_weight );
       FillUserTH1D("Eta1stEle_PAS"	        , Ele1_Eta                       , pileup_weight * gen_weight );
       FillUserTH1D("SCEta1stEle_PAS"	      , Ele1_SCEta                     , pileup_weight * gen_weight );
       FillUserTH1D("DeltaEtaEleTrk1stEle_Presel"       , fabs(Ele1_Eta-Ele1_TrkEta)                   , pileup_weight * gen_weight );
       FillUserTH1D("Phi1stEle_PAS"	        , Ele1_Phi                       , pileup_weight * gen_weight );
       FillUserTH1D("Pt2ndEle_PAS"	        , Ele2_Pt                        , pileup_weight * gen_weight );
       FillUserTH1D("PtHeep2ndEle_PAS"	    , Ele2_PtHeep                    , pileup_weight * gen_weight );
       FillUserTH1D("Eta2ndEle_PAS"	        , Ele2_Eta                       , pileup_weight * gen_weight );
       FillUserTH1D("SCEta2ndEle_PAS"	      , Ele2_SCEta                     , pileup_weight * gen_weight );
       FillUserTH1D("DeltaEtaEleTrk2ndEle_Presel"       , fabs(Ele2_Eta-Ele2_TrkEta)                   , pileup_weight * gen_weight );
       FillUserTH1D("Phi2ndEle_PAS"	    , Ele2_Phi                       , pileup_weight * gen_weight );
       FillUserTH1D("Charge1stEle_PAS"	    , Ele1_Charge                    , pileup_weight * gen_weight );
       FillUserTH1D("Charge2ndEle_PAS"	    , Ele2_Charge                    , pileup_weight * gen_weight );
       FillUserTH1D("MET_PAS"               , PFMET_Type1_Pt              , pileup_weight * gen_weight );
       FillUserTH1D("METPhi_PAS"	    , PFMET_Type1_Phi             , pileup_weight * gen_weight );
       FillUserTH1D("Pt1stJet_PAS"          , Jet1_Pt                        , pileup_weight * gen_weight );
       FillUserTH1D("Pt2ndJet_PAS"          , Jet2_Pt                        , pileup_weight * gen_weight );
       FillUserTH1D("Eta1stJet_PAS"         , Jet1_Eta                       , pileup_weight * gen_weight );
       FillUserTH1D("Eta2ndJet_PAS"         , Jet2_Eta                       , pileup_weight * gen_weight );
       FillUserTH1D("Phi1stJet_PAS"	    , Jet1_Phi                       , pileup_weight * gen_weight );
       FillUserTH1D("Phi2ndJet_PAS"	    , Jet2_Phi                       , pileup_weight * gen_weight );
       FillUserTH1D("sTlep_PAS"             , Ele1_Pt + Ele2_Pt              , pileup_weight * gen_weight );
       FillUserTH1D("sTjet_PAS"             , Jet1_Pt + Jet2_Pt              , pileup_weight * gen_weight );
       FillUserTH1D("sT_PAS"                , sT_eejj                        , pileup_weight * gen_weight );
       FillUserTH1D("sT_zjj_PAS"            , sT_zjj                         , pileup_weight * gen_weight );
       FillUserTH1D("Mjj_PAS"		    , M_j1j2                         , pileup_weight * gen_weight );
       FillUserTH1D("Mee_PAS"		    , M_e1e2                         , pileup_weight * gen_weight );
       FillUserTH1D( "MTenu_PAS"            , readerTools_->ReadValueBranch<Double_t>("MT_Ele1MET")                     , pileup_weight * gen_weight );
       FillUserTH1D("Me1j1_PAS"             , M_e1j1                         , pileup_weight * gen_weight );
       // muon kinematics
       FillUserTH1D("Pt1stMuon_PAS"	      , Muon1_Pt                        , pileup_weight * gen_weight );
       FillUserTH1D("Eta1stMuon_PAS"	    , Muon1_Eta                       , pileup_weight * gen_weight );
       FillUserTH1D("Phi1stMuon_PAS"	    , Muon1_Phi                       , pileup_weight * gen_weight );
       FillUserTH1D("Pt2ndMuon_PAS"	      , Muon2_Pt                        , pileup_weight * gen_weight );
       FillUserTH1D("Eta2ndMuon_PAS"	    , Muon2_Eta                       , pileup_weight * gen_weight );
       FillUserTH1D("Phi2ndMuon_PAS"	    , Muon2_Phi                       , pileup_weight * gen_weight );
       // scale factor dependence histos
       if ( nJet_ptCut == 2 )
         FillUserTH1D("Mee_NJetEq2_PAS", M_e1e2                         , pileup_weight * gen_weight );
       else if( nJet_ptCut == 3 )
         FillUserTH1D("Mee_NJetEq3_PAS", M_e1e2                         , pileup_weight * gen_weight );
       else if( nJet_ptCut == 4 )
         FillUserTH1D("Mee_NJetEq4_PAS", M_e1e2                         , pileup_weight * gen_weight );
       else if( nJet_ptCut == 5 )
         FillUserTH1D("Mee_NJetEq5_PAS", M_e1e2                         , pileup_weight * gen_weight );
       else if( nJet_ptCut == 6 )
         FillUserTH1D("Mee_NJetEq6_PAS", M_e1e2                         , pileup_weight * gen_weight );
       else if( nJet_ptCut == 7 )
         FillUserTH1D("Mee_NJetEq7_PAS", M_e1e2                         , pileup_weight * gen_weight );
       //
       if ( nJet_ptCut >= 3 )
         FillUserTH1D("Mee_NJetGeq3_PAS", M_e1e2                         , pileup_weight * gen_weight );
       if ( nJet_ptCut >= 4 )
         FillUserTH1D("Mee_NJetGeq4_PAS", M_e1e2                         , pileup_weight * gen_weight );
       //
       if (sT_eejj >= 300 && sT_eejj < 500)
         FillUserTH1D("Mee_sT300To500_PAS", M_e1e2                         , pileup_weight * gen_weight );
       else if (sT_eejj >= 500 && sT_eejj < 750)
         FillUserTH1D("Mee_sT500To750_PAS", M_e1e2                         , pileup_weight * gen_weight );
       else if (sT_eejj >= 750 && sT_eejj < 1250)
         FillUserTH1D("Mee_sT750To1250_PAS", M_e1e2                         , pileup_weight * gen_weight );
       else if (sT_eejj >= 1250)
         FillUserTH1D("Mee_sT1250ToInf_PAS", M_e1e2                         , pileup_weight * gen_weight );
       if (sT_eejj > 340)
         FillUserTH1D( "Mee_sT340_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 405)
         FillUserTH1D( "Mee_sT405_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 470)
         FillUserTH1D( "Mee_sT470_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 535)
         FillUserTH1D( "Mee_sT535_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 595)
         FillUserTH1D( "Mee_sT595_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 660)
         FillUserTH1D( "Mee_sT660_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 720)
         FillUserTH1D( "Mee_sT720_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 780)
         FillUserTH1D( "Mee_sT780_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 840)
         FillUserTH1D( "Mee_sT840_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 900)
         FillUserTH1D( "Mee_sT900_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 960)
         FillUserTH1D( "Mee_sT960_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1015)
         FillUserTH1D( "Mee_sT1015_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1075)
         FillUserTH1D( "Mee_sT1075_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1130)
         FillUserTH1D( "Mee_sT1130_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1190)
         FillUserTH1D( "Mee_sT1190_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1245)
         FillUserTH1D( "Mee_sT1245_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1300)
         FillUserTH1D( "Mee_sT1300_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1355)
         FillUserTH1D( "Mee_sT1355_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1410)
         FillUserTH1D( "Mee_sT1410_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1460)
         FillUserTH1D( "Mee_sT1460_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1515)
         FillUserTH1D( "Mee_sT1515_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1565)
         FillUserTH1D( "Mee_sT1565_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1615)
         FillUserTH1D( "Mee_sT1615_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1670)
         FillUserTH1D( "Mee_sT1670_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1720)
         FillUserTH1D( "Mee_sT1720_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1770)
         FillUserTH1D( "Mee_sT1770_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       if (sT_eejj > 1815)
         FillUserTH1D( "Mee_sT1815_PAS"		             ,M_e1e2                         , pileup_weight * gen_weight ); 
       ////
       //if (sT_eejj>2000)
       //  FillUserTH1D("Mee_sT2000_PAS"         ,M_e1e2                         , pileup_weight * gen_weight );
       //if (sT_eejj>926.316 && M_ej_min > 413.158 && M_e1e2 > 310.526) {
       //  FillUserTH1D("OptBinLQ600", M_e1e2, pileup_weight*gen_weight);
       //  FillUserTH1D("PtZforOptBin600",GenZGamma1_Pt, pileup_weight*gen_weight);
       //  FillUserTH1D("PtEEforOptBin600",Pt_e1e2, pileup_weight*gen_weight);
       //  FillUserTH1D("OptBinLQ600_noWeight", M_e1e2, pileup_weight*gen_weight);
       //  FillUserTH2D("PtEEVsZPtForOptBin600", Pt_e1e2, GenZGamma1_Pt, pileup_weight*gen_weight);
       //  //// printing
       //  //std::cout.precision(0);
       //  //std::cout << fixed <<  "Run = " << run << ", event = " << event << ", ls = " << ls << std::endl;
       //  //std::cout.precision(3);
       //  //std::cout << fixed <<  "  Mej      = " << M_ej_avg << std::endl;
       //  //std::cout << fixed <<  "  Mee      = " << M_e1e2 << std::endl;
       //  //std::cout << fixed <<  "  sT       = " << sT_enujj << std::endl;
       //  //std::cout << fixed <<  "  Ele1 Pt  = " << Ele1_Pt << "\t, Eta = " << Ele1_Eta << "\t, Phi = " << Ele1_Phi << std::endl;
       //  //std::cout << fixed <<  "  Ele1 dPhi  = " << Ele1_DeltaPhiTrkSC << "\t, HoE = " << Ele1_HoE << "\t, sIetaIeta = " << Ele1_SigmaIEtaIEta << std::endl;
       //  //std::cout << fixed <<  "  Ele2 Pt  = " << Ele2_Pt << "\t, Eta = " << Ele2_Eta << "\t, Phi = " << Ele2_Phi << std::endl;
       //  //std::cout << fixed <<  "  Jet1 Pt  = " << Jet1_Pt << "\t, Eta = " << Jet1_Eta << "\t, Phi = " << Jet1_Phi << std::endl;
       //  //std::cout << fixed <<  "  Jet2 Pt  = " << Jet2_Pt << "\t, Eta = " << Jet2_Eta << "\t, Phi = " << Jet2_Phi << std::endl;
       //}
       //if (sT_eejj>926.316 && M_ej_min > 534.211 && M_e1e2 > 239.474) {
       //  FillUserTH1D("OptBinLQ650", M_e1e2, pileup_weight*gen_weight);
       //  FillUserTH1D("OptBinLQ650_noWeight", M_e1e2);
       //}
       //if (sT_eejj>1105.26 && M_ej_min > 594.737 && M_e1e2 > 239.474) {
       //  FillUserTH1D("OptBinLQ700", M_e1e2, pileup_weight*gen_weight);
       //  FillUserTH1D("OptBinLQ700_noWeight", M_e1e2);
       //}
       //
       FillUserTH3D("OptimizationCutSpace",sT_eejj,M_ej_min,M_e1e2, pileup_weight*gen_weight);
       if (M_ej_min >= 100 && M_ej_min < 200)
         FillUserTH1D("Mee_MejMin100To200_PAS", M_e1e2                         , pileup_weight * gen_weight );
       else if (M_ej_min >= 200 && M_ej_min < 300)
         FillUserTH1D("Mee_MejMin200To300_PAS", M_e1e2                         , pileup_weight * gen_weight );
       else if (M_ej_min >= 300 && M_ej_min < 400)
         FillUserTH1D("Mee_MejMin300To400_PAS", M_e1e2                         , pileup_weight * gen_weight );
       else if (M_ej_min >= 400 && M_ej_min < 500)
         FillUserTH1D("Mee_MejMin400To500_PAS", M_e1e2                         , pileup_weight * gen_weight );
       else if (M_ej_min >= 500 && M_ej_min < 650)
         FillUserTH1D("Mee_MejMin500To650_PAS", M_e1e2                         , pileup_weight * gen_weight );
       else if (M_ej_min >= 650)
         FillUserTH1D("Mee_MejMin650ToInf_PAS", M_e1e2                         , pileup_weight * gen_weight );
       //
       if (Pt_e1e2 < 100)
         FillUserTH1D( "Mee_Ptee0To100_PAS"		               ,  Pt_e1e2 , pileup_weight * gen_weight ); 
       else if(Pt_e1e2 < 150 && Pt_e1e2 >= 100)
         FillUserTH1D( "Mee_Ptee100To150_PAS"		             ,  Pt_e1e2 , pileup_weight * gen_weight ); 
       else if(Pt_e1e2 < 200 && Pt_e1e2 >= 150)
         FillUserTH1D( "Mee_Ptee150To200_PAS"		             ,  Pt_e1e2 , pileup_weight * gen_weight ); 
       else if(Pt_e1e2 < 250 && Pt_e1e2 >= 200)
         FillUserTH1D( "Mee_Ptee200To250_PAS"		             ,  Pt_e1e2 , pileup_weight * gen_weight ); 
       else if(Pt_e1e2 < 300 && Pt_e1e2 >= 250)
         FillUserTH1D( "Mee_Ptee250To300_PAS"		             ,  Pt_e1e2 , pileup_weight * gen_weight ); 
       else if(Pt_e1e2 < 350 && Pt_e1e2 >= 300)
         FillUserTH1D( "Mee_Ptee300To350_PAS"		             ,  Pt_e1e2 , pileup_weight * gen_weight ); 
       else if(Pt_e1e2 < 400 && Pt_e1e2 >= 350)
         FillUserTH1D( "Mee_Ptee350To400_PAS"		             ,  Pt_e1e2 , pileup_weight * gen_weight ); 
       else
         FillUserTH1D( "Mee_Ptee400ToInf_PAS"		             ,  Pt_e1e2 , pileup_weight * gen_weight ); 
       FillUserTH1D("Me1j2_PAS"             , M_e1j2                         , pileup_weight * gen_weight );
       FillUserTH1D("Me2j1_PAS"             , M_e2j1                         , pileup_weight * gen_weight );
       FillUserTH1D("Me2j2_PAS"             , M_e2j2                         , pileup_weight * gen_weight );
       FillUserTH1D("Ptee_PAS"              , Pt_e1e2                        , pileup_weight * gen_weight );
       FillUserTH1D("DCotTheta1stEle_PAS"   , readerTools_->ReadValueBranch<Double_t>("Ele1_DCotTheta")                 , pileup_weight * gen_weight );
       FillUserTH1D("Dist1stEle_PAS"        , readerTools_->ReadValueBranch<Double_t>("Ele1_Dist")                      , pileup_weight * gen_weight );
       FillUserTH1D("DCotTheta2ndEle_PAS"   , readerTools_->ReadValueBranch<Double_t>("Ele2_DCotTheta")                 , pileup_weight * gen_weight );
       FillUserTH1D("Dist2ndEle_PAS"        , readerTools_->ReadValueBranch<Double_t>("Ele2_Dist")                      , pileup_weight * gen_weight );
       FillUserTH1D("nVertex_PAS"           , readerTools_->ReadValueBranch<Double_t>("nVertex")                        , pileup_weight * gen_weight );
       FillUserTH1D("DR_Ele1Jet1_PAS"	    , DR_Ele1Jet1                    , pileup_weight * gen_weight );
       FillUserTH1D("DR_Ele1Jet2_PAS"	    , DR_Ele1Jet2                    , pileup_weight * gen_weight );
       FillUserTH1D("DR_Ele2Jet1_PAS"	    , DR_Ele2Jet1                    , pileup_weight * gen_weight );
       FillUserTH1D("DR_Ele2Jet2_PAS"	    , DR_Ele2Jet2                    , pileup_weight * gen_weight );
       FillUserTH1D("DR_Jet1Jet2_PAS"	    , DR_Jet1Jet2                    , pileup_weight * gen_weight );
       FillUserTH1D("Meejj_PAS"             , M_eejj                         , pileup_weight * gen_weight );
       FillUserTH1D("Meej_PAS"              , M_eej                          , pileup_weight * gen_weight );
       FillUserTH1D("Mejj_PAS"              , M_ejj                          , pileup_weight * gen_weight );
       FillUserTH1D("minDR_ZJet_PAS"        , min_DeltaR_Zj                  , pileup_weight * gen_weight );
       FillUserTH1D("DR_ZJet1_PAS"          , DR_ZJ1                         , pileup_weight * gen_weight );
       FillUserTH1D("DR_ZJet2_PAS"          , DR_ZJ2                         , pileup_weight * gen_weight );
       FillUserTH1D("sTfrac_Jet1_PAS"       , Jet1_Pt / sT_eejj              , pileup_weight * gen_weight );
       FillUserTH1D("sTfrac_Jet2_PAS"       , Jet2_Pt / sT_eejj              , pileup_weight * gen_weight );
       FillUserTH1D("sTfrac_Ele1_PAS"       , Ele1_Pt / sT_eejj              , pileup_weight * gen_weight );
       FillUserTH1D("sTfrac_Ele2_PAS"       , Ele2_Pt / sT_eejj              , pileup_weight * gen_weight );
       FillUserTH1D("sTfrac_Jet_PAS"        , ( Jet1_Pt + Jet2_Pt ) / sT_eejj, pileup_weight * gen_weight );
       FillUserTH1D("sTfrac_Ele_PAS"        , ( Ele1_Pt + Ele2_Pt ) / sT_eejj, pileup_weight * gen_weight );
       FillUserTH1D("Mej_selected_avg_PAS"  , M_ej_avg                       , pileup_weight * gen_weight );	   
       FillUserTH1D("Mej_selected_min_PAS"  , M_ej_min                       , pileup_weight * gen_weight );	   
       FillUserTH1D("Mej_selected_max_PAS"  , M_ej_max                       , pileup_weight * gen_weight );	   
       FillUserTH1D("Mej_minmax_PAS"        , M_ej_min                       , pileup_weight * gen_weight );	   
       FillUserTH1D("Mej_minmax_PAS"        , M_ej_max                       , pileup_weight * gen_weight );	   

       FillUserTH2D("MeeVsST_PAS" , M_e1e2, sT_eejj, pileup_weight * gen_weight );	   
       FillUserTH2D("MeeVsPtee_PAS" , M_e1e2, Pt_e1e2, pileup_weight * gen_weight );	   

       // checking electrons
       FillUserTH1D( "SCEta_1stEle_Presel",Ele1_SCEta, pileup_weight*gen_weight);
       FillUserTH1D( "EleEta_1stEle_Presel",Ele1_Eta,pileup_weight*gen_weight);
       FillUserTH1D( "DeltaEtaTrkSC_1stEle_Presel",Ele1_DeltaEtaTrkSC,pileup_weight*gen_weight);
       FillUserTH1D( "SCEtaMinusEleEta_1stEle_Presel",Ele1_SCEta-Ele1_Eta,pileup_weight*gen_weight);
       FillUserTH1D( "DeltaEtaTrkSC_Minus_SCEtaMinusEleEta_1stEle_Presel",Ele1_DeltaEtaTrkSC-(Ele1_SCEta-Ele1_Eta),pileup_weight*gen_weight);
       FillUserTH1D( "SCEta_2ndEle_Presel",Ele2_SCEta, pileup_weight*gen_weight);
       FillUserTH1D( "EleEta_2ndEle_Presel",Ele2_Eta,pileup_weight*gen_weight);
       FillUserTH1D( "DeltaEtaTrkSC_2ndEle_Presel",Ele2_DeltaEtaTrkSC,pileup_weight*gen_weight);
       FillUserTH1D( "SCEtaMinusEleEta_2ndEle_Presel",Ele2_SCEta-Ele2_Eta,pileup_weight*gen_weight);
       FillUserTH1D( "DeltaEtaTrkSC_Minus_SCEtaMinusEleEta_2ndEle_Presel",Ele2_DeltaEtaTrkSC-(Ele2_SCEta-Ele2_Eta),pileup_weight*gen_weight);
       //--------------------------------------------------------------------------
       // Mass-pairing histograms at preselection
       //--------------------------------------------------------------------------

       if ( fabs(M_e1j1-M_e2j2) < fabs(M_e1j2-M_e2j1) ) {
         FillUserTH1D("Me1j_selected_PAS"   , M_e1j1          , pileup_weight * gen_weight);	   
         FillUserTH1D("Me2j_selected_PAS"   , M_e2j2          , pileup_weight * gen_weight);	   
         FillUserTH2D("Me1jVsMe2j_selected" , M_e1j1  , M_e2j2, pileup_weight * gen_weight);
         FillUserTH2D("Me1jVsMe2j_rejected" , M_e1j2  , M_e2j1, pileup_weight * gen_weight);
       }
       else {
         FillUserTH1D("Me1j_selected_PAS"   , M_e1j2          , pileup_weight * gen_weight);	   
         FillUserTH1D("Me2j_selected_PAS"   , M_e2j1          , pileup_weight * gen_weight);	   
         FillUserTH2D("Me1jVsMe2j_selected" , M_e1j2  , M_e2j1, pileup_weight * gen_weight);
         FillUserTH2D("Me1jVsMe2j_rejected" , M_e1j1  , M_e2j2, pileup_weight * gen_weight);
       }

       //--------------------------------------------------------------------------
       // Preselection + data-only
       //--------------------------------------------------------------------------

       if ( isData() ) { 
         FillUserTH1D("run_PAS"  , run );
         profile_run_vs_nvtx_PAS -> Fill ( run, nVertex, 1 );
       }

       //--------------------------------------------------------------------------
       // Preselection + N(Jet) > 2 
       //--------------------------------------------------------------------------

       if ( nJet_ptCut > 2 ){ 
         FillUserTH1D( "M_e1j3_PAS"  , M_e1j3, pileup_weight * gen_weight ) ;
         FillUserTH1D( "M_e2j3_PAS"  , M_e2j3, pileup_weight * gen_weight ) ;
         FillUserTH1D( "M_j1j3_PAS"  , M_j1j3, pileup_weight * gen_weight ) ;
         FillUserTH1D( "M_j2j3_PAS"  , M_j2j3, pileup_weight * gen_weight ) ;
         FillUserTH1D( "M_eejjj_PAS" , M_eejjj,pileup_weight * gen_weight ) ;

         FillUserTH1D( "Ptj1j2j3_PAS"            , Pt_j1j2j3           , pileup_weight * gen_weight );
         FillUserTH1D( "Ptj2j3_PAS"              , Pt_j2j3             , pileup_weight * gen_weight );
         FillUserTH1D( "Ptj1j3_PAS"              , Pt_j1j3             , pileup_weight * gen_weight );
         FillUserTH1D( "Ptee_Minus_Ptj1j2j3_PAS" , Pt_e1e2 - Pt_j1j2j3 , pileup_weight * gen_weight ); 
       }

       //--------------------------------------------------------------------------
       // Preselection + event type (EBEB, EEEB, EEEE, etc)
       //--------------------------------------------------------------------------

       if      ( isEB   ) FillUserTH1D( "Mee_EB_PAS"  , M_e1e2, pileup_weight * gen_weight ); 
       if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_PAS", M_e1e2, pileup_weight * gen_weight ); 
       else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_PAS", M_e1e2, pileup_weight * gen_weight ); 
       else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_PAS", M_e1e2, pileup_weight * gen_weight ); 

       //--------------------------------------------------------------------------
       // Preselection + high ST plot
       //--------------------------------------------------------------------------

       if ( sT_eejj > 445. ) FillUserTH1D("Mee_PASandST445"   , M_e1e2  , pileup_weight * gen_weight ) ;

       //--------------------------------------------------------------------------
       // High M(ee) plots
       //--------------------------------------------------------------------------

       if ( M_e1e2  > 100. ) FillUserTH1D("sT_PASandMee100"   , sT_eejj , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 110. ) FillUserTH1D("sT_PASandMee110"   , sT_eejj , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 120. ) FillUserTH1D("sT_PASandMee120"   , sT_eejj , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 130. ) FillUserTH1D("sT_PASandMee130"   , sT_eejj , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 140. ) FillUserTH1D("sT_PASandMee140"   , sT_eejj , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 150. ) FillUserTH1D("sT_PASandMee150"   , sT_eejj , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 160. ) FillUserTH1D("sT_PASandMee160"   , sT_eejj , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 170. ) FillUserTH1D("sT_PASandMee170"   , sT_eejj , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 180. ) FillUserTH1D("sT_PASandMee180"   , sT_eejj , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 190. ) FillUserTH1D("sT_PASandMee190"   , sT_eejj , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 200. ) FillUserTH1D("sT_PASandMee200"   , sT_eejj , pileup_weight * gen_weight ) ;

       if ( M_e1e2 > 100. ) { 

         FillUserTH1D("BeamSpotDXY_1stEle_PASandMee100"           , Ele1_BeamSpotDXY                    , pileup_weight * gen_weight    ); 
         FillUserTH1D("Classif_1stEle_PASandMee100"               , Ele1_Classif                        , pileup_weight * gen_weight    ); 
         FillUserTH1D("CorrIsolation_1stEle_PASandMee100"         , Ele1_CorrIsolation                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("DeltaEtaTrkSC_1stEle_PASandMee100"         , Ele1_DeltaEtaTrkSC                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("DeltaPhiTrkSC_1stEle_PASandMee100"         , Ele1_DeltaPhiTrkSC                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("Full5x5E1x5OverE5x5_1stEle_PASandMee100"   , Ele1_Full5x5E1x5OverE5x5            , pileup_weight * gen_weight    ); 
         FillUserTH1D("Full5x5E2x5OverE5x5_1stEle_PASandMee100"   , Ele1_Full5x5E2x5OverE5x5            , pileup_weight * gen_weight    ); 
         FillUserTH1D("EcalIsolation_1stEle_PASandMee100"         , Ele1_EcalIsolation                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("HcalIsolation_1stEle_PASandMee100"         , Ele1_HcalIsolation                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("TrkIsolation_1stEle_PASandMee100"          , Ele1_TrkIsolation                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("Energy_1stEle_PASandMee100"                , Ele1_Energy                         , pileup_weight * gen_weight    ); 
         FillUserTH1D("FBrem_1stEle_PASandMee100"                 , Ele1_FBrem                          , pileup_weight * gen_weight    ); 
         FillUserTH1D("GsfCtfCharge_1stEle_PASandMee100"          , Ele1_GsfCtfCharge                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("GsfCtfScPixCharge_1stEle_PASandMee100"     , Ele1_GsfCtfScPixCharge              , pileup_weight * gen_weight    ); 
         FillUserTH1D("GsfScPixCharge_1stEle_PASandMee100"        , Ele1_GsfScPixCharge                 , pileup_weight * gen_weight    ); 
         FillUserTH1D("HasMatchedPhot_1stEle_PASandMee100"        , Ele1_HasMatchedPhot                 , pileup_weight * gen_weight    ); 
         FillUserTH1D("HoE_1stEle_PASandMee100"                   , Ele1_HoE                            , pileup_weight * gen_weight    ); 
         FillUserTH1D("LeadVtxDistXY_1stEle_PASandMee100"         , Ele1_LeadVtxDistXY                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("LeadVtxDistZ_1stEle_PASandMee100"          , Ele1_LeadVtxDistZ                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("MissingHits_1stEle_PASandMee100"           , Ele1_MissingHits                    , pileup_weight * gen_weight    ); 
         FillUserTH1D("NBrems_1stEle_PASandMee100"                , Ele1_NBrems                         , pileup_weight * gen_weight    ); 
         FillUserTH1D("ValidFrac_1stEle_PASandMee100"             , Ele1_ValidFrac                      , pileup_weight * gen_weight    ); 
         FillUserTH1D("EnergyORawEnergy_1stEle_PASandMee100"      , Ele1_Energy / Ele1_RawEnergy        , pileup_weight * gen_weight    ); 
         FillUserTH1D("TrkPtOPt_1stEle_PASandMee100"              , Ele1_TrkPt  / Ele1_Pt               , pileup_weight * gen_weight    ); 
         if ( fabs(Ele1_Eta) < eleEta_bar ) { 
           FillUserTH1D("SigmaEtaEta_Barrel_1stEle_PASandMee100"  , Ele1_SigmaEtaEta                    , pileup_weight * gen_weight    ); 
           FillUserTH1D("Full5x5SigmaIEtaIEta_Barrel_1stEle_PASandMee100", Ele1_Full5x5SigmaIEtaIEta                  , pileup_weight * gen_weight    ); 
         }
         else if ( fabs(Ele1_Eta) > eleEta_end1_min && fabs(Ele2_Eta) < eleEta_end2_max ){
           FillUserTH1D("SigmaEtaEta_Endcap_1stEle_PASandMee100"  , Ele1_SigmaEtaEta                    , pileup_weight * gen_weight    ); 
           FillUserTH1D("Full5x5SigmaIEtaIEta_Endcap_1stEle_PASandMee100", Ele1_Full5x5SigmaIEtaIEta                  , pileup_weight * gen_weight    ); 
         }

         FillUserTH1D("BeamSpotDXY_2ndEle_PASandMee100"           , Ele2_BeamSpotDXY                    , pileup_weight * gen_weight    ); 
         FillUserTH1D("Classif_2ndEle_PASandMee100"               , Ele2_Classif                        , pileup_weight * gen_weight    ); 
         FillUserTH1D("CorrIsolation_2ndEle_PASandMee100"         , Ele2_CorrIsolation                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("DeltaEtaTrkSC_2ndEle_PASandMee100"         , Ele2_DeltaEtaTrkSC                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("DeltaPhiTrkSC_2ndEle_PASandMee100"         , Ele2_DeltaPhiTrkSC                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("Full5x5E1x5OverE5x5_2ndEle_PASandMee100"   , Ele2_Full5x5E1x5OverE5x5            , pileup_weight * gen_weight    ); 
         FillUserTH1D("Full5x5E2x5OverE5x5_2ndEle_PASandMee100"   , Ele2_Full5x5E2x5OverE5x5            , pileup_weight * gen_weight    ); 
         FillUserTH1D("EcalIsolation_2ndEle_PASandMee100"         , Ele2_EcalIsolation                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("HcalIsolation_2ndEle_PASandMee100"         , Ele2_HcalIsolation                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("TrkIsolation_2ndEle_PASandMee100"          , Ele2_TrkIsolation                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("Energy_2ndEle_PASandMee100"                , Ele2_Energy                         , pileup_weight * gen_weight    ); 
         FillUserTH1D("FBrem_2ndEle_PASandMee100"                 , Ele2_FBrem                          , pileup_weight * gen_weight    ); 
         FillUserTH1D("GsfCtfCharge_2ndEle_PASandMee100"          , Ele2_GsfCtfCharge                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("GsfCtfScPixCharge_2ndEle_PASandMee100"     , Ele2_GsfCtfScPixCharge              , pileup_weight * gen_weight    ); 
         FillUserTH1D("GsfScPixCharge_2ndEle_PASandMee100"        , Ele2_GsfScPixCharge                 , pileup_weight * gen_weight    ); 
         FillUserTH1D("HasMatchedPhot_2ndEle_PASandMee100"        , Ele2_HasMatchedPhot                 , pileup_weight * gen_weight    ); 
         FillUserTH1D("HoE_2ndEle_PASandMee100"                   , Ele2_HoE                            , pileup_weight * gen_weight    ); 
         FillUserTH1D("LeadVtxDistXY_2ndEle_PASandMee100"         , Ele2_LeadVtxDistXY                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("LeadVtxDistZ_2ndEle_PASandMee100"          , Ele2_LeadVtxDistZ                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("MissingHits_2ndEle_PASandMee100"           , Ele2_MissingHits                    , pileup_weight * gen_weight    ); 
         FillUserTH1D("NBrems_2ndEle_PASandMee100"                , Ele2_NBrems                         , pileup_weight * gen_weight    ); 
         FillUserTH1D("ValidFrac_2ndEle_PASandMee100"             , Ele2_ValidFrac                      , pileup_weight * gen_weight    ); 
         FillUserTH1D("EnergyORawEnergy_2ndEle_PASandMee100"      , Ele2_Energy / Ele2_RawEnergy        , pileup_weight * gen_weight    ); 
         FillUserTH1D("TrkPtOPt_2ndEle_PASandMee100"              , Ele2_TrkPt  / Ele2_Pt               , pileup_weight * gen_weight    ); 
         if ( fabs(Ele2_Eta) < eleEta_bar ) { 
           FillUserTH1D("SigmaEtaEta_Barrel_2ndEle_PASandMee100"  , Ele2_SigmaEtaEta                    , pileup_weight * gen_weight    ); 
           FillUserTH1D("Full5x5SigmaIEtaIEta_Barrel_2ndEle_PASandMee100", Ele2_Full5x5SigmaIEtaIEta                  , pileup_weight * gen_weight    ); 
         }
         else if ( fabs(Ele2_Eta) > eleEta_end1_min && fabs(Ele2_Eta) < eleEta_end2_max ){
           FillUserTH1D("SigmaEtaEta_Endcap_2ndEle_PASandMee100"  , Ele2_SigmaEtaEta                    , pileup_weight * gen_weight    ); 
           FillUserTH1D("Full5x5SigmaIEtaIEta_Endcap_2ndEle_PASandMee100", Ele2_Full5x5SigmaIEtaIEta                  , pileup_weight * gen_weight    ); 
         }

         FillUserTH1D("Ptee_PASandMee100"              , Pt_e1e2                        , pileup_weight * gen_weight );
         FillUserTH2D("MeeVsST_PASandMee100" , M_e1e2, sT_eejj, pileup_weight * gen_weight );	   
         FillUserTH1D("sT_zjj_PASandMee100"            , sT_zjj                         , pileup_weight * gen_weight );
         FillUserTH1D("nVertex_PASandMee100"           , nVertex                        , pileup_weight * gen_weight );
         FillUserTH1D("EleChargeSum_PASandMee100"      , Ele1_Charge + Ele2_Charge      , pileup_weight * gen_weight ) ;
         FillUserTH1D("nJet_PASandMee100"              , nJet_ptCut                     , pileup_weight * gen_weight ) ;
         FillUserTH1D("sTlep_PASandMee100"             , Ele1_Pt + Ele2_Pt              , pileup_weight * gen_weight ) ;
         FillUserTH1D("sTjet_PASandMee100"             , Jet1_Pt + Jet2_Pt              , pileup_weight * gen_weight ) ;
         FillUserTH1D("Mjj_PASandMee100"               , M_j1j2                         , pileup_weight * gen_weight ) ;
         FillUserTH1D("Me1j1_PASandMee100"             , M_e1j1                         , pileup_weight * gen_weight ) ;
         FillUserTH1D("Me1j2_PASandMee100"             , M_e1j2                         , pileup_weight * gen_weight ) ;
         FillUserTH1D("Me2j1_PASandMee100"             , M_e1j1                         , pileup_weight * gen_weight ) ;
         FillUserTH1D("Me2j2_PASandMee100"             , M_e1j2                         , pileup_weight * gen_weight ) ;
         FillUserTH1D("Pt1stEle_PASandMee100"          , Ele1_Pt                        , pileup_weight * gen_weight ) ;
         FillUserTH1D("Pt2ndEle_PASandMee100"          , Ele2_Pt                        , pileup_weight * gen_weight ) ;
         FillUserTH1D("Pt1stJet_PASandMee100"          , Jet1_Pt                        , pileup_weight * gen_weight ) ;
         FillUserTH1D("Pt2ndJet_PASandMee100"          , Jet2_Pt                        , pileup_weight * gen_weight ) ;
         FillUserTH1D("Mej_selected_avg_PASandMee100"  , M_ej_avg                       , pileup_weight * gen_weight ) ;
         FillUserTH1D("sTfrac_Jet1_PASandMee100"       , Jet1_Pt / sT_eejj              , pileup_weight * gen_weight ) ;
         FillUserTH1D("sTfrac_Jet2_PASandMee100"       , Jet2_Pt / sT_eejj              , pileup_weight * gen_weight ) ;
         FillUserTH1D("sTfrac_Ele1_PASandMee100"       , Ele1_Pt / sT_eejj              , pileup_weight * gen_weight ) ;
         FillUserTH1D("sTfrac_Ele2_PASandMee100"       , Ele2_Pt / sT_eejj              , pileup_weight * gen_weight ) ;
         FillUserTH1D("sTfrac_Jet_PASandMee100"        , ( Jet1_Pt + Jet2_Pt ) / sT_eejj, pileup_weight * gen_weight ) ;
         FillUserTH1D("sTfrac_Ele_PASandMee100"        , ( Ele1_Pt + Ele2_Pt ) / sT_eejj, pileup_weight * gen_weight ) ;
         FillUserTH1D("Ptj1j2_PASandMee100"            , Pt_j1j2                        , pileup_weight * gen_weight ) ;
         FillUserTH1D("Ptee_Minus_Ptj1j2_PASandMee100" , Pt_e1e2 - Pt_j1j2              , pileup_weight * gen_weight ) ;

         if ( nJet_ptCut > 2 ) { 
           FillUserTH1D( "M_e1j3_PASandMee100" , M_e1j3, pileup_weight * gen_weight ) ;
           FillUserTH1D( "M_e2j3_PASandMee100" , M_e2j3, pileup_weight * gen_weight ) ;
           FillUserTH1D( "M_j1j3_PASandMee100" , M_j1j3, pileup_weight * gen_weight ) ;
           FillUserTH1D( "M_j2j3_PASandMee100" , M_j2j3, pileup_weight * gen_weight ) ;
           FillUserTH1D( "M_eejjj_PASandMee100", M_eejjj,pileup_weight * gen_weight ) ;
           FillUserTH1D( "Ptj1j2j3_PASandMee100"            , Pt_j1j2j3           , pileup_weight * gen_weight );
           FillUserTH1D( "Ptj2j3_PASandMee100"              , Pt_j2j3             , pileup_weight * gen_weight );
           FillUserTH1D( "Ptj1j3_PASandMee100"              , Pt_j1j3             , pileup_weight * gen_weight );
           FillUserTH1D( "Ptee_Minus_Ptj1j2j3_PASandMee100" , Pt_e1e2 - Pt_j1j2j3 , pileup_weight * gen_weight ); 
         }
       }

       //--------------------------------------------------------------------------
       // Preselection + M(ee) normalization region plots
       //--------------------------------------------------------------------------

       if ( M_e1e2 > 80.0 && M_e1e2 < 100.0 ){
         FillUserTH1D("Mee_80_100_Preselection", M_e1e2, pileup_weight * gen_weight ) ;
         if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_80_100_PAS", M_e1e2, pileup_weight * gen_weight ); 
         else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_80_100_PAS", M_e1e2, pileup_weight * gen_weight ); 
         else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_80_100_PAS", M_e1e2, pileup_weight * gen_weight ); 
         if      ( isEB   ) FillUserTH1D( "Mee_EB_80_100_PAS"  , M_e1e2, pileup_weight * gen_weight ); 
       } 

       if ( M_e1e2 > 70.0 && M_e1e2 < 110.0 ){
         FillUserTH1D("Mee_70_110_Preselection", M_e1e2, pileup_weight * gen_weight ) ;
         if ( sT_eejj > 600 ) FillUserTH1D("Mee_70_110_ST600_Preselection", M_e1e2, pileup_weight * gen_weight ) ;
         if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_70_110_PAS", M_e1e2, pileup_weight * gen_weight ); 
         else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_70_110_PAS", M_e1e2, pileup_weight * gen_weight ); 
         else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_70_110_PAS", M_e1e2, pileup_weight * gen_weight ); 
         if      ( isEB   ) FillUserTH1D( "Mee_EB_70_110_PAS"  , M_e1e2, pileup_weight * gen_weight ); 

         //FillUserTH1D( "ProcessID_ZWindow", ProcessID, pileup_weight * gen_weight );
         //if ( ProcessID == 0 ) FillUserTH1D ( "Mee_70_110_Preselection_Process0", M_e1e2, pileup_weight * gen_weight );
         //if ( ProcessID == 1 ) FillUserTH1D ( "Mee_70_110_Preselection_Process1", M_e1e2, pileup_weight * gen_weight );
         //if ( ProcessID == 2 ) FillUserTH1D ( "Mee_70_110_Preselection_Process2", M_e1e2, pileup_weight * gen_weight );
         //if ( ProcessID == 3 ) FillUserTH1D ( "Mee_70_110_Preselection_Process3", M_e1e2, pileup_weight * gen_weight );
         //if ( ProcessID == 4 ) FillUserTH1D ( "Mee_70_110_Preselection_Process4", M_e1e2, pileup_weight * gen_weight );
       }

       //--------------------------------------------------------------------------
       // Region of interest plots
       //-------------------------------------------------------------------------- 

       if ( do_roi_plots && passed_region_of_interest ) { 


         FillUserTH1D("BeamSpotDXY_1stEle_ROI"           , Ele1_BeamSpotDXY                    , pileup_weight * gen_weight    ); 
         FillUserTH1D("Classif_1stEle_ROI"               , Ele1_Classif                        , pileup_weight * gen_weight    ); 
         FillUserTH1D("CorrIsolation_1stEle_ROI"         , Ele1_CorrIsolation                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("DeltaEtaTrkSC_1stEle_ROI"         , Ele1_DeltaEtaTrkSC                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("DeltaPhiTrkSC_1stEle_ROI"         , Ele1_DeltaPhiTrkSC                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("E1x5OverE5x5_1stEle_ROI"          , Ele1_Full5x5E1x5OverE5x5                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("E2x5OverE5x5_1stEle_ROI"          , Ele1_Full5x5E2x5OverE5x5                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("EcalIsolation_1stEle_ROI"         , Ele1_EcalIsolation                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("HcalIsolation_1stEle_ROI"         , Ele1_HcalIsolation                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("TrkIsolation_1stEle_ROI"          , Ele1_TrkIsolation                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("Energy_1stEle_ROI"                , Ele1_Energy                         , pileup_weight * gen_weight    ); 
         FillUserTH1D("FBrem_1stEle_ROI"                 , Ele1_FBrem                          , pileup_weight * gen_weight    ); 
         FillUserTH1D("GsfCtfCharge_1stEle_ROI"          , Ele1_GsfCtfCharge                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("GsfCtfScPixCharge_1stEle_ROI"     , Ele1_GsfCtfScPixCharge              , pileup_weight * gen_weight    ); 
         FillUserTH1D("GsfScPixCharge_1stEle_ROI"        , Ele1_GsfScPixCharge                 , pileup_weight * gen_weight    ); 
         FillUserTH1D("HasMatchedPhot_1stEle_ROI"        , Ele1_HasMatchedPhot                 , pileup_weight * gen_weight    ); 
         FillUserTH1D("HoE_1stEle_ROI"                   , Ele1_HoE                            , pileup_weight * gen_weight    ); 
         FillUserTH1D("LeadVtxDistXY_1stEle_ROI"         , Ele1_LeadVtxDistXY                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("LeadVtxDistZ_1stEle_ROI"          , Ele1_LeadVtxDistZ                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("MissingHits_1stEle_ROI"           , Ele1_MissingHits                    , pileup_weight * gen_weight    ); 
         FillUserTH1D("NBrems_1stEle_ROI"                , Ele1_NBrems                         , pileup_weight * gen_weight    ); 
         FillUserTH1D("ValidFrac_1stEle_ROI"             , Ele1_ValidFrac                      , pileup_weight * gen_weight    ); 
         FillUserTH1D("EnergyORawEnergy_1stEle_ROI"      , Ele1_Energy / Ele1_RawEnergy        , pileup_weight * gen_weight    ); 
         FillUserTH1D("TrkPtOPt_1stEle_ROI"              , Ele1_TrkPt  / Ele1_Pt               , pileup_weight * gen_weight    ); 
         if ( fabs(Ele1_Eta) < eleEta_bar ) { 
           FillUserTH1D("SigmaEtaEta_Barrel_1stEle_ROI"  , Ele1_SigmaEtaEta                    , pileup_weight * gen_weight    ); 
           FillUserTH1D("SigmaIEtaIEta_Barrel_1stEle_ROI", Ele1_Full5x5SigmaIEtaIEta                  , pileup_weight * gen_weight    ); 
         }
         else if ( fabs(Ele1_Eta) > eleEta_end1_min && fabs(Ele2_Eta) < eleEta_end2_max ){
           FillUserTH1D("SigmaEtaEta_Endcap_1stEle_ROI"  , Ele1_SigmaEtaEta                    , pileup_weight * gen_weight    ); 
           FillUserTH1D("SigmaIEtaIEta_Endcap_1stEle_ROI", Ele1_Full5x5SigmaIEtaIEta                  , pileup_weight * gen_weight    ); 
         }

         FillUserTH1D("BeamSpotDXY_2ndEle_ROI"           , Ele2_BeamSpotDXY                    , pileup_weight * gen_weight    ); 
         FillUserTH1D("Classif_2ndEle_ROI"               , Ele2_Classif                        , pileup_weight * gen_weight    ); 
         FillUserTH1D("CorrIsolation_2ndEle_ROI"         , Ele2_CorrIsolation                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("DeltaEtaTrkSC_2ndEle_ROI"         , Ele2_DeltaEtaTrkSC                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("DeltaPhiTrkSC_2ndEle_ROI"         , Ele2_DeltaPhiTrkSC                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("E1x5OverE5x5_2ndEle_ROI"          , Ele2_Full5x5E1x5OverE5x5                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("E2x5OverE5x5_2ndEle_ROI"          , Ele2_Full5x5E2x5OverE5x5                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("EcalIsolation_2ndEle_ROI"         , Ele2_EcalIsolation                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("HcalIsolation_2ndEle_ROI"         , Ele2_HcalIsolation                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("TrkIsolation_2ndEle_ROI"          , Ele2_TrkIsolation                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("Energy_2ndEle_ROI"                , Ele2_Energy                         , pileup_weight * gen_weight    ); 
         FillUserTH1D("FBrem_2ndEle_ROI"                 , Ele2_FBrem                          , pileup_weight * gen_weight    ); 
         FillUserTH1D("GsfCtfCharge_2ndEle_ROI"          , Ele2_GsfCtfCharge                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("GsfCtfScPixCharge_2ndEle_ROI"     , Ele2_GsfCtfScPixCharge              , pileup_weight * gen_weight    ); 
         FillUserTH1D("GsfScPixCharge_2ndEle_ROI"        , Ele2_GsfScPixCharge                 , pileup_weight * gen_weight    ); 
         FillUserTH1D("HasMatchedPhot_2ndEle_ROI"        , Ele2_HasMatchedPhot                 , pileup_weight * gen_weight    ); 
         FillUserTH1D("HoE_2ndEle_ROI"                   , Ele2_HoE                            , pileup_weight * gen_weight    ); 
         FillUserTH1D("LeadVtxDistXY_2ndEle_ROI"         , Ele2_LeadVtxDistXY                  , pileup_weight * gen_weight    ); 
         FillUserTH1D("LeadVtxDistZ_2ndEle_ROI"          , Ele2_LeadVtxDistZ                   , pileup_weight * gen_weight    ); 
         FillUserTH1D("MissingHits_2ndEle_ROI"           , Ele2_MissingHits                    , pileup_weight * gen_weight    ); 
         FillUserTH1D("NBrems_2ndEle_ROI"                , Ele2_NBrems                         , pileup_weight * gen_weight    ); 
         FillUserTH1D("ValidFrac_2ndEle_ROI"             , Ele2_ValidFrac                      , pileup_weight * gen_weight    ); 
         FillUserTH1D("EnergyORawEnergy_2ndEle_ROI"      , Ele2_Energy / Ele2_RawEnergy        , pileup_weight * gen_weight    ); 
         FillUserTH1D("TrkPtOPt_2ndEle_ROI"              , Ele2_TrkPt  / Ele2_Pt               , pileup_weight * gen_weight    ); 
         if ( fabs(Ele2_Eta) < eleEta_bar ) { 
           FillUserTH1D("SigmaEtaEta_Barrel_2ndEle_ROI"  , Ele2_SigmaEtaEta                    , pileup_weight * gen_weight    ); 
           FillUserTH1D("SigmaIEtaIEta_Barrel_2ndEle_ROI", Ele2_Full5x5SigmaIEtaIEta                  , pileup_weight * gen_weight    ); 
         }
         else if ( fabs(Ele2_Eta) > eleEta_end1_min && fabs(Ele2_Eta) < eleEta_end2_max ){
           FillUserTH1D("SigmaEtaEta_Endcap_2ndEle_ROI"  , Ele2_SigmaEtaEta                    , pileup_weight * gen_weight    ); 
           FillUserTH1D("SigmaIEtaIEta_Endcap_2ndEle_ROI", Ele2_Full5x5SigmaIEtaIEta                  , pileup_weight * gen_weight    ); 
         }


         FillUserTH1D("Me1j1_ROI"           , M_e1j1                         , pileup_weight * gen_weight );
         FillUserTH1D("Me1j2_ROI"           , M_e1j2                         , pileup_weight * gen_weight );
         FillUserTH1D("Me2j1_ROI"           , M_e2j1                         , pileup_weight * gen_weight );
         FillUserTH1D("Me2j2_ROI"           , M_e2j2                         , pileup_weight * gen_weight );
         FillUserTH1D("Ptee_ROI"            , Pt_e1e2                        , pileup_weight * gen_weight );
         FillUserTH1D("Eta1stJet_ROI"       , Jet1_Eta                       , pileup_weight * gen_weight );
         FillUserTH1D("Eta2ndJet_ROI"       , Jet2_Eta                       , pileup_weight * gen_weight );
         FillUserTH1D("Eta1stEle_ROI"	    , Ele1_Eta                       , pileup_weight * gen_weight );
         FillUserTH1D("Eta2ndEle_ROI"	    , Ele2_Eta                       , pileup_weight * gen_weight );
         FillUserTH1D("Phi1stJet_ROI"       , Jet1_Phi                       , pileup_weight * gen_weight );
         FillUserTH1D("Phi2ndJet_ROI"       , Jet2_Phi                       , pileup_weight * gen_weight );
         FillUserTH1D("Phi1stEle_ROI"	    , Ele1_Phi                       , pileup_weight * gen_weight );
         FillUserTH1D("Phi2ndEle_ROI"	    , Ele2_Phi                       , pileup_weight * gen_weight );
         FillUserTH2D("MeeVsST_ROI" , M_e1e2, sT_eejj, pileup_weight * gen_weight );	   
         FillUserTH1D("Mee_ROI"		    , M_e1e2                         , pileup_weight * gen_weight );
         FillUserTH1D("sT_zjj_ROI"          , sT_zjj                         , pileup_weight * gen_weight );
         FillUserTH1D("nVertex_ROI"         , nVertex                        , pileup_weight * gen_weight );
         FillUserTH1D("nJet_ROI"            , nJet_ptCut                     , pileup_weight * gen_weight );
         FillUserTH1D("EleChargeSum_ROI"    , Ele1_Charge + Ele2_Charge      , pileup_weight * gen_weight );
         FillUserTH1D("Meejj_ROI"           , M_eejj                         , pileup_weight * gen_weight );
         FillUserTH1D("Meej_ROI"            , M_eej                          , pileup_weight * gen_weight );
         FillUserTH1D("Mejj_ROI"            , M_ejj                          , pileup_weight * gen_weight );
         FillUserTH1D("Mjj_ROI"             , M_j1j2                         , pileup_weight * gen_weight );
         FillUserTH1D("Mej_selected_avg_ROI", M_ej_avg                       , pileup_weight * gen_weight );
         FillUserTH1D("minDR_ZJet_ROI"      , min_DeltaR_Zj                  , pileup_weight * gen_weight );
         FillUserTH1D("DR_ZJet1_ROI"        , DR_ZJ1                         , pileup_weight * gen_weight );
         FillUserTH1D("DR_ZJet2_ROI"        , DR_ZJ2                         , pileup_weight * gen_weight );
         FillUserTH1D("MET_ROI"             , PFMET_Type1_Pt              , pileup_weight * gen_weight );
         FillUserTH1D("sT_ROI"              , sT_eejj                        , pileup_weight * gen_weight );
         FillUserTH1D("sTlep_ROI"           , Ele1_Pt + Ele2_Pt              , pileup_weight * gen_weight );
         FillUserTH1D("sTjet_ROI"           , Jet1_Pt + Jet2_Pt              , pileup_weight * gen_weight );
         FillUserTH1D("Pt1stEle_ROI"        , Ele1_Pt                        , pileup_weight * gen_weight );
         FillUserTH1D("Pt2ndEle_ROI"        , Ele2_Pt                        , pileup_weight * gen_weight );
         FillUserTH1D("Pt1stJet_ROI"        , Jet1_Pt                        , pileup_weight * gen_weight );
         FillUserTH1D("Pt2ndJet_ROI"        , Jet2_Pt                        , pileup_weight * gen_weight );
         FillUserTH1D( "sTfrac_Jet1_ROI"    , Jet1_Pt / sT_eejj              , pileup_weight * gen_weight );
         FillUserTH1D( "sTfrac_Jet2_ROI"    , Jet2_Pt / sT_eejj              , pileup_weight * gen_weight );
         FillUserTH1D( "sTfrac_Ele1_ROI"    , Ele1_Pt / sT_eejj              , pileup_weight * gen_weight );
         FillUserTH1D( "sTfrac_Ele2_ROI"    , Ele2_Pt / sT_eejj              , pileup_weight * gen_weight );
         FillUserTH1D( "sTfrac_Jet_ROI"     , ( Jet1_Pt + Jet2_Pt ) / sT_eejj, pileup_weight * gen_weight );
         FillUserTH1D( "sTfrac_Ele_ROI"     , ( Ele1_Pt + Ele2_Pt ) / sT_eejj, pileup_weight * gen_weight );
         FillUserTH1D("Ptj1j2_ROI"            , Pt_j1j2                        , pileup_weight * gen_weight ) ;
         FillUserTH1D("Ptee_Minus_Ptj1j2_ROI" , Pt_e1e2 - Pt_j1j2              , pileup_weight * gen_weight ) ;


         if ( nJet_ptCut > 2 ) { 
           FillUserTH1D( "M_e1j3_ROI"  , M_e1j3,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "M_e2j3_ROI"  , M_e2j3,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "M_j1j3_ROI"  , M_j1j3,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "M_j2j3_ROI"  , M_j2j3,  pileup_weight * gen_weight ) ;
           FillUserTH1D( "M_eejjj_ROI" , M_eejjj, pileup_weight * gen_weight ) ;
           FillUserTH1D( "Ptj1j2j3_ROI"            , Pt_j1j2j3           , pileup_weight * gen_weight );
           FillUserTH1D( "Ptj2j3_ROI"              , Pt_j2j3             , pileup_weight * gen_weight );
           FillUserTH1D( "Ptj1j3_ROI"              , Pt_j1j3             , pileup_weight * gen_weight );
           FillUserTH1D( "Ptee_Minus_Ptj1j2j3_ROI" , Pt_e1e2 - Pt_j1j2j3 , pileup_weight * gen_weight ); 
         }
       }

       //-------------------------------------------------------------------------- 
       // Final selection plots
       //-------------------------------------------------------------------------- 

       if(doFinalSelections)
       {
         for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
           int  lq_mass = LQ_MASS      [i_lq_mass];
           bool pass    = passed_vector[i_lq_mass];
           if ( !pass ) continue;

           sprintf(plot_name, "Mej_selected_avg_LQ%d"       , lq_mass ); FillUserTH1D ( plot_name, M_ej_avg          , pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_selected_min_LQ%d"       , lq_mass ); FillUserTH1D ( plot_name, M_ej_min          , pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_selected_max_LQ%d"       , lq_mass ); FillUserTH1D ( plot_name, M_ej_max          , pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_minmax_LQ%d"             , lq_mass ); FillUserTH1D ( plot_name, M_ej_min          , pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_minmax_LQ%d"             , lq_mass ); FillUserTH1D ( plot_name, M_ej_max          , pileup_weight * gen_weight );
           sprintf(plot_name, "sT_eejj_LQ%d"                , lq_mass ); FillUserTH1D ( plot_name, sT_eejj           , pileup_weight * gen_weight );
           sprintf(plot_name, "Mee_LQ%d"                    , lq_mass ); FillUserTH1D ( plot_name, M_e1e2            , pileup_weight * gen_weight );
           sprintf(plot_name, "DR_Ele1Jet1_LQ%d"            , lq_mass ); FillUserTH1D ( plot_name, DR_Ele1Jet1       , pileup_weight * gen_weight );
           sprintf(plot_name, "Mej_selected_min_vs_max_LQ%d", lq_mass ); FillUserTH2D ( plot_name, M_ej_min, M_ej_max, pileup_weight * gen_weight );


           sprintf(plot_name, "BeamSpotDXY_1stEle_LQ%d"        , lq_mass );   FillUserTH1D(plot_name,  Ele1_BeamSpotDXY               , pileup_weight * gen_weight ); 
           sprintf(plot_name, "Classif_1stEle_LQ%d"            , lq_mass );   FillUserTH1D(plot_name,  Ele1_Classif                   , pileup_weight * gen_weight ); 
           sprintf(plot_name, "CorrIsolation_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele1_CorrIsolation             , pileup_weight * gen_weight ); 
           sprintf(plot_name, "DeltaEtaTrkSC_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele1_DeltaEtaTrkSC             , pileup_weight * gen_weight ); 
           sprintf(plot_name, "DeltaPhiTrkSC_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele1_DeltaPhiTrkSC             , pileup_weight * gen_weight ); 
           sprintf(plot_name, "Full5x5E1x5OverE5x5_1stEle_LQ%d", lq_mass );   FillUserTH1D(plot_name,  Ele1_Full5x5E1x5OverE5x5       , pileup_weight * gen_weight ); 
           sprintf(plot_name, "Full5x5E2x5OverE5x5_1stEle_LQ%d", lq_mass );   FillUserTH1D(plot_name,  Ele1_Full5x5E2x5OverE5x5       , pileup_weight * gen_weight ); 
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

           if ( fabs(Ele1_Eta) < eleEta_bar ) { 
             sprintf(plot_name, "SigmaEtaEta_Barrel_1stEle_LQ%d"  , lq_mass ); FillUserTH1D( plot_name , Ele1_SigmaEtaEta   , pileup_weight * gen_weight    ); 
             sprintf(plot_name, "Full5x5SigmaIEtaIEta_Barrel_1stEle_LQ%d", lq_mass ); FillUserTH1D( plot_name , Ele1_Full5x5SigmaIEtaIEta , pileup_weight * gen_weight    ); 
           }
           else if ( fabs(Ele1_Eta) > eleEta_end1_min && fabs(Ele2_Eta) < eleEta_end2_max ){
             sprintf(plot_name, "SigmaEtaEta_Endcap_1stEle_LQ%d"  , lq_mass ); FillUserTH1D( plot_name , Ele1_SigmaEtaEta   , pileup_weight * gen_weight    ); 
             sprintf(plot_name, "Full5x5SigmaIEtaIEta_Endcap_1stEle_LQ%d", lq_mass ); FillUserTH1D( plot_name , Ele1_Full5x5SigmaIEtaIEta , pileup_weight * gen_weight    ); 
           }

           sprintf(plot_name, "BeamSpotDXY_2ndEle_LQ%d"        , lq_mass );   FillUserTH1D(plot_name,  Ele2_BeamSpotDXY               , pileup_weight * gen_weight ); 
           sprintf(plot_name, "Classif_2ndEle_LQ%d"            , lq_mass );   FillUserTH1D(plot_name,  Ele2_Classif                   , pileup_weight * gen_weight ); 
           sprintf(plot_name, "CorrIsolation_2ndEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele2_CorrIsolation             , pileup_weight * gen_weight ); 
           sprintf(plot_name, "DeltaEtaTrkSC_2ndEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele2_DeltaEtaTrkSC             , pileup_weight * gen_weight ); 
           sprintf(plot_name, "DeltaPhiTrkSC_2ndEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele2_DeltaPhiTrkSC             , pileup_weight * gen_weight ); 
           sprintf(plot_name, "Full5x5E1x5OverE5x5_2ndEle_LQ%d", lq_mass );   FillUserTH1D(plot_name,  Ele2_Full5x5E1x5OverE5x5       , pileup_weight * gen_weight ); 
           sprintf(plot_name, "Full5x5E2x5OverE5x5_2ndEle_LQ%d", lq_mass );   FillUserTH1D(plot_name,  Ele2_Full5x5E2x5OverE5x5       , pileup_weight * gen_weight ); 
           sprintf(plot_name, "EcalIsolation_2ndEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele2_EcalIsolation             , pileup_weight * gen_weight ); 
           sprintf(plot_name, "HcalIsolation_2ndEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele2_HcalIsolation             , pileup_weight * gen_weight ); 
           sprintf(plot_name, "TrkIsolation_2ndEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  Ele2_TrkIsolation              , pileup_weight * gen_weight ); 
           sprintf(plot_name, "Energy_2ndEle_LQ%d"             , lq_mass );   FillUserTH1D(plot_name,  Ele2_Energy                    , pileup_weight * gen_weight ); 
           sprintf(plot_name, "FBrem_2ndEle_LQ%d"              , lq_mass );   FillUserTH1D(plot_name,  Ele2_FBrem                     , pileup_weight * gen_weight ); 
           sprintf(plot_name, "GsfCtfCharge_2ndEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  Ele2_GsfCtfCharge              , pileup_weight * gen_weight ); 
           sprintf(plot_name, "GsfCtfScPixCharge_2ndEle_LQ%d"  , lq_mass );   FillUserTH1D(plot_name,  Ele2_GsfCtfScPixCharge         , pileup_weight * gen_weight ); 
           sprintf(plot_name, "GsfScPixCharge_2ndEle_LQ%d"     , lq_mass );   FillUserTH1D(plot_name,  Ele2_GsfScPixCharge            , pileup_weight * gen_weight ); 
           sprintf(plot_name, "HasMatchedPhot_2ndEle_LQ%d"     , lq_mass );   FillUserTH1D(plot_name,  Ele2_HasMatchedPhot            , pileup_weight * gen_weight ); 
           sprintf(plot_name, "HoE_2ndEle_LQ%d"                , lq_mass );   FillUserTH1D(plot_name,  Ele2_HoE                       , pileup_weight * gen_weight ); 
           sprintf(plot_name, "LeadVtxDistXY_2ndEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  Ele2_LeadVtxDistXY             , pileup_weight * gen_weight ); 
           sprintf(plot_name, "LeadVtxDistZ_2ndEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  Ele2_LeadVtxDistZ              , pileup_weight * gen_weight ); 
           sprintf(plot_name, "MissingHits_2ndEle_LQ%d"        , lq_mass );   FillUserTH1D(plot_name,  Ele2_MissingHits               , pileup_weight * gen_weight ); 
           sprintf(plot_name, "NBrems_2ndEle_LQ%d"             , lq_mass );   FillUserTH1D(plot_name,  Ele2_NBrems                    , pileup_weight * gen_weight ); 
           sprintf(plot_name, "ValidFrac_2ndEle_LQ%d"          , lq_mass );   FillUserTH1D(plot_name,  Ele2_ValidFrac                 , pileup_weight * gen_weight ); 
           sprintf(plot_name, "EnergyORawEnergy_2ndEle_LQ%d"   , lq_mass );   FillUserTH1D(plot_name,  Ele2_Energy / Ele2_RawEnergy   , pileup_weight * gen_weight ); 
           sprintf(plot_name, "TrkPtOPt_2ndEle_LQ%d"           , lq_mass );   FillUserTH1D(plot_name,  Ele2_TrkPt  / Ele2_Pt          , pileup_weight * gen_weight ); 

           if ( fabs(Ele2_Eta) < eleEta_bar ) { 
             sprintf(plot_name, "SigmaEtaEta_Barrel_2ndEle_LQ%d"  , lq_mass ); FillUserTH1D( plot_name , Ele2_SigmaEtaEta   , pileup_weight * gen_weight    ); 
             sprintf(plot_name, "Full5x5SigmaIEtaIEta_Barrel_2ndEle_LQ%d", lq_mass ); FillUserTH1D( plot_name , Ele2_Full5x5SigmaIEtaIEta , pileup_weight * gen_weight    ); 
           }
           else if ( fabs(Ele2_Eta) > eleEta_end2_min && fabs(Ele2_Eta) < eleEta_end2_max ){
             sprintf(plot_name, "SigmaEtaEta_Endcap_2ndEle_LQ%d"  , lq_mass ); FillUserTH1D( plot_name , Ele2_SigmaEtaEta   , pileup_weight * gen_weight    ); 
             sprintf(plot_name, "Full5x5SigmaIEtaIEta_Endcap_2ndEle_LQ%d", lq_mass ); FillUserTH1D( plot_name , Ele2_Full5x5SigmaIEtaIEta , pileup_weight * gen_weight    ); 
           }

           sprintf(plot_name, "Me1j1_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , M_e1j1                         , pileup_weight * gen_weight );
           sprintf(plot_name, "Me1j2_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , M_e1j2                         , pileup_weight * gen_weight );
           sprintf(plot_name, "Me2j1_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , M_e2j1                         , pileup_weight * gen_weight );
           sprintf(plot_name, "Me2j2_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , M_e2j2                         , pileup_weight * gen_weight );
           sprintf(plot_name, "Ptee_LQ%d"              , lq_mass ); FillUserTH1D( plot_name , Pt_e1e2                        , pileup_weight * gen_weight );
           sprintf(plot_name, "Eta1stJet_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Jet1_Eta                       , pileup_weight * gen_weight );
           sprintf(plot_name, "Eta2ndJet_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Jet2_Eta                       , pileup_weight * gen_weight );
           sprintf(plot_name, "Eta1stEle_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Ele1_Eta                       , pileup_weight * gen_weight );
           sprintf(plot_name, "Eta2ndEle_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Ele2_Eta                       , pileup_weight * gen_weight );
           sprintf(plot_name, "Phi1stJet_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Jet1_Phi                       , pileup_weight * gen_weight );
           sprintf(plot_name, "Phi2ndJet_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Jet2_Phi                       , pileup_weight * gen_weight );
           sprintf(plot_name, "Phi1stEle_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Ele1_Phi                       , pileup_weight * gen_weight );
           sprintf(plot_name, "Phi2ndEle_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Ele2_Phi                       , pileup_weight * gen_weight );
           sprintf(plot_name, "MeeVsST_LQ%d"           , lq_mass ); FillUserTH2D( plot_name , M_e1e2, sT_eejj                , pileup_weight * gen_weight );	   
           sprintf(plot_name, "sT_zjj_LQ%d"            , lq_mass ); FillUserTH1D( plot_name , sT_zjj                         , pileup_weight * gen_weight );
           sprintf(plot_name, "nVertex_LQ%d"           , lq_mass ); FillUserTH1D( plot_name , nVertex                        , pileup_weight * gen_weight );
           sprintf(plot_name, "nJet_LQ%d"              , lq_mass ); FillUserTH1D( plot_name , nJet_ptCut                     , pileup_weight * gen_weight );
           sprintf(plot_name, "EleChargeSum_LQ%d"      , lq_mass ); FillUserTH1D( plot_name , Ele1_Charge + Ele2_Charge      , pileup_weight * gen_weight );
           sprintf(plot_name, "Meejj_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , M_eejj                         , pileup_weight * gen_weight );
           sprintf(plot_name, "Meej_LQ%d"              , lq_mass ); FillUserTH1D( plot_name , M_eej                          , pileup_weight * gen_weight );
           sprintf(plot_name, "Mejj_LQ%d"              , lq_mass ); FillUserTH1D( plot_name , M_ejj                          , pileup_weight * gen_weight );
           sprintf(plot_name, "Mjj_LQ%d"               , lq_mass ); FillUserTH1D( plot_name , M_j1j2                         , pileup_weight * gen_weight );
           sprintf(plot_name, "minDR_ZJet_LQ%d"        , lq_mass ); FillUserTH1D( plot_name , min_DeltaR_Zj                  , pileup_weight * gen_weight );
           sprintf(plot_name, "DR_ZJet1_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , DR_ZJ1                         , pileup_weight * gen_weight );
           sprintf(plot_name, "DR_ZJet2_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , DR_ZJ2                         , pileup_weight * gen_weight );
           sprintf(plot_name, "MET_LQ%d"               , lq_mass ); FillUserTH1D( plot_name , PFMET_Type1_Pt              , pileup_weight * gen_weight );
           sprintf(plot_name, "sTlep_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , Ele1_Pt + Ele2_Pt              , pileup_weight * gen_weight );
           sprintf(plot_name, "sTjet_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , Jet1_Pt + Jet2_Pt              , pileup_weight * gen_weight );
           sprintf(plot_name, "Pt1stEle_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Ele1_Pt                        , pileup_weight * gen_weight );
           sprintf(plot_name, "Pt2ndEle_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Ele2_Pt                        , pileup_weight * gen_weight );
           sprintf(plot_name, "Pt1stJet_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Jet1_Pt                        , pileup_weight * gen_weight );
           sprintf(plot_name, "Pt2ndJet_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Jet2_Pt                        , pileup_weight * gen_weight );
           sprintf(plot_name, "sTfrac_Jet1_LQ%d"       , lq_mass ); FillUserTH1D( plot_name , Jet1_Pt / sT_eejj              , pileup_weight * gen_weight );
           sprintf(plot_name, "sTfrac_Jet2_LQ%d"       , lq_mass ); FillUserTH1D( plot_name , Jet2_Pt / sT_eejj              , pileup_weight * gen_weight );
           sprintf(plot_name, "sTfrac_Ele1_LQ%d"       , lq_mass ); FillUserTH1D( plot_name , Ele1_Pt / sT_eejj              , pileup_weight * gen_weight );
           sprintf(plot_name, "sTfrac_Ele2_LQ%d"       , lq_mass ); FillUserTH1D( plot_name , Ele2_Pt / sT_eejj              , pileup_weight * gen_weight );
           sprintf(plot_name, "sTfrac_Jet_LQ%d"        , lq_mass ); FillUserTH1D( plot_name , ( Jet1_Pt + Jet2_Pt ) / sT_eejj, pileup_weight * gen_weight );
           sprintf(plot_name, "sTfrac_Ele_LQ%d"        , lq_mass ); FillUserTH1D( plot_name , ( Ele1_Pt + Ele2_Pt ) / sT_eejj, pileup_weight * gen_weight );
           sprintf(plot_name, "Ptj1j2_LQ%d"            , lq_mass ); FillUserTH1D( plot_name , Pt_j1j2                        , pileup_weight * gen_weight );
           sprintf(plot_name, "Ptee_Minus_Ptj1j2_LQ%d" , lq_mass ); FillUserTH1D( plot_name , Pt_e1e2 - Pt_j1j2              , pileup_weight * gen_weight );
           // checking electrons
           sprintf(plot_name,"SCEta_1stEle_LQ%d" , lq_mass );                                FillUserTH1D(plot_name,Ele1_SCEta, pileup_weight*gen_weight);
           sprintf(plot_name,"EleEta_1stEle_LQ%d" , lq_mass );                               FillUserTH1D(plot_name,Ele1_Eta,pileup_weight*gen_weight);
           sprintf(plot_name,"DeltaEtaTrkSC_1stEle_LQ%d" , lq_mass );                        FillUserTH1D(plot_name,Ele1_DeltaEtaTrkSC,pileup_weight*gen_weight);
           sprintf(plot_name,"SCEtaMinusEleEta_1stEle_LQ%d" , lq_mass );                     FillUserTH1D(plot_name,Ele1_SCEta-Ele1_Eta,pileup_weight*gen_weight);
           sprintf(plot_name,"DeltaEtaTrkSC_Minus_SCEtaMinusEleEta_1stEle_LQ%d" , lq_mass ); FillUserTH1D(plot_name,Ele1_DeltaEtaTrkSC-(Ele1_SCEta-Ele1_Eta),pileup_weight*gen_weight);
           sprintf(plot_name,"SCEta_2ndEle_LQ%d" , lq_mass );                                FillUserTH1D(plot_name,Ele2_SCEta, pileup_weight*gen_weight);
           sprintf(plot_name,"EleEta_2ndEle_LQ%d" , lq_mass );                               FillUserTH1D(plot_name,Ele2_Eta,pileup_weight*gen_weight);
           sprintf(plot_name,"DeltaEtaTrkSC_2ndEle_LQ%d" , lq_mass );                        FillUserTH1D(plot_name,Ele2_DeltaEtaTrkSC,pileup_weight*gen_weight);
           sprintf(plot_name,"SCEtaMinusEleEta_2ndEle_LQ%d", lq_mass );                      FillUserTH1D(plot_name,Ele2_SCEta-Ele2_Eta,pileup_weight*gen_weight);
           sprintf(plot_name,"DeltaEtaTrkSC_Minus_SCEtaMinusEleEta_2ndEle_LQ%d" , lq_mass ); FillUserTH1D(plot_name,Ele2_DeltaEtaTrkSC-(Ele2_SCEta-Ele2_Eta),pileup_weight*gen_weight);
           // muon kinematics
           sprintf(plot_name, "Eta1stMuon_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Muon1_Eta                       , pileup_weight * gen_weight );
           sprintf(plot_name, "Eta2ndMuon_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Muon2_Eta                       , pileup_weight * gen_weight );
           sprintf(plot_name, "Phi1stMuon_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Muon1_Phi                       , pileup_weight * gen_weight );
           sprintf(plot_name, "Pt1stMuon_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Muon1_Pt                        , pileup_weight * gen_weight );
           sprintf(plot_name, "Pt2ndMuon_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Muon2_Pt                        , pileup_weight * gen_weight );

         } // End final selection

         if( hasCut("sT_eejj_LQ300") && passedCut("sT_eejj_LQ300") && passedCut("min_M_ej_LQ300"))
           FillUserTH1D("Mee_70_110_LQ300", M_e1e2 , pileup_weight * gen_weight );
         if( hasCut("sT_eejj_LQ600") && passedCut("sT_eejj_LQ600") && passedCut("min_M_ej_LQ600"))
           FillUserTH1D("Mee_70_110_LQ600", M_e1e2 , pileup_weight * gen_weight );
         if( hasCut("sT_eejj_LQ800") && passedCut("sT_eejj_LQ800") && passedCut("min_M_ej_LQ800"))
           FillUserTH1D("Mee_70_110_LQ800", M_e1e2 , pileup_weight * gen_weight );
         if( hasCut("sT_eejj_LQ900") && passedCut("sT_eejj_LQ900") && passedCut("min_M_ej_LQ900"))
           FillUserTH1D("Mee_70_110_LQ900", M_e1e2 , pileup_weight * gen_weight );
         if( hasCut("sT_eejj_LQ1000") && passedCut("sT_eejj_LQ1000") && passedCut("min_M_ej_LQ1000"))
           FillUserTH1D("Mee_70_110_LQ1000", M_e1e2 , pileup_weight * gen_weight );

       }

       //--------------------------------------------------------------------------
       // Fill skim tree
       //--------------------------------------------------------------------------
       fillSkimTree();

     } // End preselection 
   } // End loop over events
   
   output_root_ -> cd();
   profile_run_vs_nvtx_HLT -> Write();
   profile_run_vs_nvtx_PAS -> Write();
   

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

