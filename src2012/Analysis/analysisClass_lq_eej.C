#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TProfile.h>

int run_max[38] = { 191718, 193621, 194223, 194533, 194912, 195304, 195552, 195930, 196239, 196453, 
                    196531, 198941, 199319, 199571, 199812, 200042, 200369, 200991, 201278, 201705, 
                    202045, 202209, 202477, 203002, 204250, 204601, 205339, 205694, 206187, 206389, 
                    206512, 206745, 207214, 207372, 207515, 208307, 208487, 208686 };



int run_min[38] = { 190645, 191720, 193834, 194224, 194619, 194914, 195378, 195633, 195937, 196249, 
                    196495, 198049, 198954, 199336, 199572, 199833, 200049, 200381, 200992, 201554, 
                    201706, 202054, 202237, 202478, 203894, 204511, 205086, 205344, 205718, 206188, 
                    206391, 206513, 206859, 207217, 207397, 207517, 208339, 208509 };


int run_max_1fb[20] = { 193621, 194533, 195304, 195915, 196452, 196531, 199319, 199804, 200244, 201196, 
			202014, 202305, 203002, 204577, 205666, 206246, 206598, 207273, 207905, 208686 };



int run_min_1fb[20] = { 190645, 193834, 194619, 195378, 195916, 196453, 198049, 199336, 199812, 200245, 
			201197, 202016, 202314, 203894, 204599, 205667, 206257, 206605, 207279, 207920 };


int get_split ( int run ) { 
  
  int n_split = 38;
  int split = -999;
  
  for (int i_split = 0; i_split < n_split; ++i_split ){ 
    int tmp_run_min = run_min [i_split];
    int tmp_run_max = run_max [i_split];
    if ( run >= tmp_run_min &&
	 run <= tmp_run_max ) {
      split = i_split;
    }
  }

  if ( split < 0 ) std::cout << "ERROR: could not find a split for run = " << run << std::endl;
  
  return split;
}



int get_split_1fb ( int run ) { 
  
  int n_split = 20;
  int split = -999;
  
  for (int i_split = 0; i_split < n_split; ++i_split ){ 
    int tmp_run_min = run_min_1fb [i_split];
    int tmp_run_max = run_max_1fb [i_split];
    if ( run >= tmp_run_min &&
	 run <= tmp_run_max ) {
      split = i_split;
    }
  }

  if ( split < 0 ) std::cout << "ERROR: could not find a split (1 fb-1) for run = " << run << std::endl;
  
  return split;
}


analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile){}

analysisClass::~analysisClass(){}

void analysisClass::Loop()
{
   std::cout << "analysisClass::Loop() begins" <<std::endl;   

   //--------------------------------------------------------------------------
   // Final selection mass points
   //--------------------------------------------------------------------------

   const int n_lq_mass = 19;

   int LQ_MASS[n_lq_mass] = { 
     300,  350,  400, 450, 500, 550,  600,  650,
     700,  750,  800, 850, 900, 950, 1000, 1050,
     1100, 1150, 1200
   };

   std::vector<bool> passed_vector;
   
   //--------------------------------------------------------------------------
   // Decide which plots to save (default is to save everything)
   //--------------------------------------------------------------------------
   
   fillSkim                         ( !true  ) ;
   fillAllPreviousCuts              ( !true  ) ;
   fillAllOtherCuts                 ( !true  ) ;
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;

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

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------
   
   CreateUserTH1D( "split_PAS"                       ,    38    , -0.5    , 37.5     );
   CreateUserTH1D( "split_1fb_PAS"                   ,    38    , -0.5    , 37.5     );
   CreateUserTH1D( "split_ROI"                       ,    38    , -0.5    , 37.5     );
   CreateUserTH1D( "split_PASandMee100"              ,    38    , -0.5    , 37.5     );
   CreateUserTH1D( "sTfrac_Jet1_PAS"                 ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele1_PAS"                 ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele2_PAS"                 ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Jet_PAS"                  ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele_PAS"                  ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Jet1_PASandMee100"        ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele1_PASandMee100"        ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele2_PASandMee100"        ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Jet_PASandMee100"         ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele_PASandMee100"         ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "EleChargeSum_PAS"                ,    3     , -2.5    , 2.5      );
   CreateUserTH1D( "EleChargeSum_PASandMee100"       ,    3     , -2.5    , 2.5      );
   CreateUserTH1D( "EleChargeSum_ROI"                ,    3     , -2.5    , 2.5      );
   CreateUserTH1D( "sTfrac_Jet1_ROI"                 ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele1_ROI"                 ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele2_ROI"                 ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Jet_ROI"                  ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele_ROI"                  ,   100    ,  0.0    , 1.0      );
   CreateUserTH1D( "ProcessID"                       ,    21    , -0.5    , 20.5     );
   CreateUserTH1D( "ProcessID_PAS"                   ,    21    , -0.5    , 20.5     );
   CreateUserTH1D( "ProcessID_ZWindow"               ,    21    , -0.5    , 20.5     );
   CreateUserTH1D( "nElectron_PAS"                   ,    5     , -0.5    , 4.5      );
   CreateUserTH1D( "nMuon_PAS"                       ,    5     , -0.5    , 4.5      );
   CreateUserTH1D( "nJet_PAS"                        ,    10    , -0.5    , 9.5      );
   CreateUserTH1D( "nJet_PASandMee100"               ,    10    , -0.5    , 9.5      );
   CreateUserTH1D( "nJet_ROI"                        ,    10    , -0.5    , 9.5      );
   CreateUserTH1D( "Pt1stEle_PAS"	             , 	100    , 0       , 1000     ); 
   CreateUserTH1D( "Pt1stEle_PASandMee100"           , 	100    , 0       , 1000     ); 
   CreateUserTH1D( "Pt1stEle_ROI"	             , 	100    , 0       , 1000     ); 
   CreateUserTH1D( "Eta1stEle_PAS"	             , 	100    , -5      , 5	  ); 
   CreateUserTH1D( "Eta1stEle_ROI"	             , 	100    , -5      , 5	  ); 
   CreateUserTH1D( "Phi1stEle_PAS"	             , 	60     , -3.1416 , +3.1416  ); 
   CreateUserTH1D( "Phi1stEle_ROI"	             , 	60     , -3.1416 , +3.1416  ); 
   CreateUserTH1D( "Pt2ndEle_PAS"	             , 	100    , 0       , 1000     ); 
   CreateUserTH1D( "Pt2ndEle_PASandMee100"           , 	100    , 0       , 1000     ); 
   CreateUserTH1D( "Pt2ndEle_ROI"	             , 	100    , 0       , 1000     ); 
   CreateUserTH1D( "Eta2ndEle_PAS"	             , 	100    , -5      , 5	  ); 
   CreateUserTH1D( "Eta2ndEle_ROI"	             , 	100    , -5      , 5	  ); 
   CreateUserTH1D( "Phi2ndEle_PAS"	             , 	60     , -3.1416 , +3.1416  ); 
   CreateUserTH1D( "Phi2ndEle_ROI"	             , 	60     , -3.1416 , +3.1416  ); 
   CreateUserTH1D( "Charge1stEle_PAS"	             , 	2      , -1.0001 , 1.0001	  ); 
   CreateUserTH1D( "Charge2ndEle_PAS"	             , 	2      , -1.0001 , 1.0001	  ); 
   CreateUserTH1D( "MET_PAS"                         ,    200   , 0       , 1000	  ); 
   CreateUserTH1D( "MET_ROI"                         ,    200   , 0       , 1000	  ); 
   CreateUserTH1D( "METPhi_PAS"		             , 	60     , -3.1416 , +3.1416  ); 
   CreateUserTH1D( "Pt1stJet_PAS"                    ,    100   , 0       , 1000	  ); 
   CreateUserTH1D( "Pt1stJet_PASandMee100"           ,    100   , 0       , 1000	  ); 
   CreateUserTH1D( "Pt1stJet_ROI"                    ,    100   , 0       , 1000	  ); 
   CreateUserTH1D( "Eta1stJet_PAS"                   ,    100   , -5      , 5	  ); 
   CreateUserTH1D( "Eta1stJet_ROI"                   ,    100   , -5      , 5	  ); 
   CreateUserTH1D( "Phi1stJet_PAS"	             , 	 60    , -3.1416 , +3.1416  ); 
   CreateUserTH1D( "Phi1stJet_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
   CreateUserTH1D( "sTlep_PAS"                       ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sTlep_PASandMee100"              ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sTlep_ROI"                       ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sTjet_PAS"                       ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sTjet_PASandMee100"              ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sTjet_ROI"                       ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "sT_PAS"                          ,    200   , 0       , 2000	  ); 
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
   CreateUserTH1D( "sT_ROI"                          ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mjj_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mjj_PASandMee100"	             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mjj_ROI"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_PAS"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_ROI"		             ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Mee_PASandST445"                 ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "MTenu_PAS"                       ,    200   , 0       , 1000	  ); 
   CreateUserTH1D( "Me1j1_PAS"                       ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Me2j1_PAS"                       ,    200   , 0       , 2000	  ); 
   CreateUserTH1D( "Meej_PAS"                        ,    400   , 0       , 4000   );
   CreateUserTH1D( "Meej_ROI"                        ,    400   , 0       , 4000   );
   CreateUserTH1D( "run_PAS"                         ,    20000 , 160300  , 180300 );
   CreateUserTH1D( "run_HLT"                         ,    20000 , 160300  , 180300 );
			
   CreateUserTH1D( "Ptee_PAS"                        ,    200 , 0       , 2000     );
   CreateUserTH1D( "Ptee_PASandMee100"               ,    200 , 0       , 2000     );
   CreateUserTH1D( "Ptee_ROI"                        ,    200 , 0       , 2000     );
					             
   CreateUserTH1D( "Me1j1_PASandMee100"             ,    200 , 0       , 2000   );   
   CreateUserTH1D( "Me1j1_ROI"                      ,    200 , 0       , 2000   );
   					             
   CreateUserTH1D( "DCotTheta1stEle_PAS"             ,    100 , 0.0, 1.0);
   CreateUserTH1D( "Dist1stEle_PAS"                  ,    100 , 0.0, 1.0);  
   CreateUserTH1D( "DCotTheta2ndEle_PAS"             ,    100 , 0.0, 1.0);
   CreateUserTH1D( "Dist2ndEle_PAS"                  ,    100 , 0.0, 1.0);  
		                                     
   CreateUserTH1D( "nVertex_PAS"                     ,    101   , -0.5   , 100.5	 ) ; 
   CreateUserTH1D( "nVertex_PASandMee100"            ,    101   , -0.5   , 100.5	 ) ; 
   CreateUserTH1D( "nVertex_ROI"                     ,    101   , -0.5   , 100.5	 ) ; 
   
   CreateUserTH1D( "DR_Ele1Jet1_PAS"	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet1_PAS"	   , 	getHistoNBins("DR_Ele2Jet1"), getHistoMin("DR_Ele2Jet1"), getHistoMax("DR_Ele2Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele1Ele2_PAS"	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
   
   CreateUserTH1D( "DR_ZJet1_PAS"        ,    getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
   CreateUserTH1D( "DR_ZJet1_ROI"        ,    getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 


   CreateUserTH2D( "MeeVsST_PAS"                 ,     200, 0, 2000, 200, 0, 2000) ;
   CreateUserTH2D( "MeeVsST_PASandMee100"        ,     200, 0, 2000, 200, 0, 2000) ;
   CreateUserTH2D( "MeeVsST_ROI"                 ,     200, 0, 2000, 200, 0, 2000) ;

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

   //--------------------------------------------------------------------------
   // Loop over the chain
   //--------------------------------------------------------------------------

   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntries();
   //Long64_t nentries = 5;
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
     
     int    passedJSON = passJSON ( run, ls , isData ) ;

     //--------------------------------------------------------------------------
     // Do pileup re-weighting
     //--------------------------------------------------------------------------
     
     double pileup_weight = getPileupWeight ( nPileUpInt_True, isData ) ;
     
     //--------------------------------------------------------------------------
     // Get information about gen-level reweighting (should be for Sherpa only)
     //--------------------------------------------------------------------------

     double gen_weight = Weight;
     if ( isData ) gen_weight = 1.0;
     if ( isData && Ele2_ValidFrac > 998. ){
       gen_weight = 0.0;
       if      (  60.0 < M_e1e2 < 120. ) gen_weight = 0.61;
       else if ( 120.0 < M_e1e2 < 200. ) gen_weight = 0.42;
       else if ( 200.0 < M_e1e2        ) gen_weight = 0.42;
     }

     // std::cout << "Gen weight = " << int ( 1.0 / gen_weight ) << std::endl;

     //--------------------------------------------------------------------------
     // First variable to fill just shows the "reweighting".  Always passes.
     //--------------------------------------------------------------------------

     fillVariableWithValue ( "Reweighting", 1, gen_weight * pileup_weight  );

     //--------------------------------------------------------------------------
     // Fill JSON variable
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue ("PassJSON", passedJSON, gen_weight * pileup_weight); 

     //--------------------------------------------------------------------------
     // Fill noise filters
     //--------------------------------------------------------------------------

     // Noise filters
     fillVariableWithValue(   "PassHBHENoiseFilter"	      , PassHBHENoiseFilter                              , gen_weight * pileup_weight );
     fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight                          , gen_weight * pileup_weight );
     fillVariableWithValue(   "PassBadEESupercrystalFilter"   , ( isData == 1 ) ? PassBadEESupercrystalFilter : 1, gen_weight * pileup_weight );
     fillVariableWithValue(   "PassBeamScraping"	      , ( isData == 1 ) ? PassBeamScraping	      : 1, gen_weight * pileup_weight );
     fillVariableWithValue(   "PassEcalDeadCellBoundEnergy"   , PassEcalDeadCellBoundEnergy                      , gen_weight * pileup_weight );
     fillVariableWithValue(   "PassEcalDeadCellTrigPrim"      , PassEcalDeadCellTrigPrim                         , gen_weight * pileup_weight );
     fillVariableWithValue(   "PassEcalLaserCorrFilter"       , ( isData == 1 ) ? PassEcalLaserCorrFilter     : 1, gen_weight * pileup_weight );
     fillVariableWithValue(   "PassHcalLaserEventFilter"      , ( isData == 1 ) ? PassHcalLaserEventFilter    : 1, gen_weight * pileup_weight );
     fillVariableWithValue(   "PassPhysDecl"		      , ( isData == 1 ) ? PassPhysDecl		      : 1, gen_weight * pileup_weight );
     fillVariableWithValue(   "PassPrimaryVertex"	      , PassPrimaryVertex                                , gen_weight * pileup_weight );
     fillVariableWithValue(   "PassTrackingFailure"	      , ( isData == 1 ) ? PassTrackingFailure	      : 1, gen_weight * pileup_weight );
     
     //--------------------------------------------------------------------------
     // Fill HLT
     //--------------------------------------------------------------------------

     int passHLT = 1;
     if ( isData ) { 
       passHLT = 0;
       if ( H_DoubleEle33_CIdL_GsfIdVL == 1 ) { 
       	 passHLT = 1;
       }
     }
     
     fillVariableWithValue ( "PassHLT", passHLT, gen_weight * pileup_weight  ) ;     

     //--------------------------------------------------------------------------
     // Calculate variables for trigger matching 
     //--------------------------------------------------------------------------

     int nEle_hltMatched = 0.0;
     if ( Ele1_hltEleSignalPt > 0.0 ) nEle_hltMatched++;
     if ( Ele2_hltEleSignalPt > 0.0 ) nEle_hltMatched++;
     
     int nJet_hltMatched = 0.0;
     if ( Jet1_hltNoPUJetPt > 0.0 || Jet1_hltJetPt > 0.0 ) nJet_hltMatched++;
     if ( Jet2_hltNoPUJetPt > 0.0 || Jet2_hltJetPt > 0.0 ) nJet_hltMatched++;

     fillVariableWithValue("nEle_hltMatched",nEle_hltMatched, gen_weight * pileup_weight  );
     fillVariableWithValue("nJet_hltMatched",nJet_hltMatched, gen_weight * pileup_weight  );
     
     //--------------------------------------------------------------------------
     // Pass number of muons & electrons 
     // --> Special consideration if ttbar is derived from data
     //--------------------------------------------------------------------------

     // Muons and electrons
     bool is_ttbar_from_data = false;
     if ( Ele2_ValidFrac > 998. ) is_ttbar_from_data = true;
     
     int PassNEle = 0;
     if ( !is_ttbar_from_data && nEle_ptCut == 2 ) PassNEle = 1;
     if (  is_ttbar_from_data && nEle_ptCut == 2 ) PassNEle = 1;

     int PassNMuon = 0;
     if ( !is_ttbar_from_data && nMuon_ptCut == 0 ) PassNMuon = 1;
     if (  is_ttbar_from_data && nMuon_ptCut >  0 ) PassNMuon = 1;

     fillVariableWithValue("PassNEle" , PassNEle , gen_weight * pileup_weight);
     fillVariableWithValue("PassNMuon", PassNMuon, gen_weight * pileup_weight);

     //--------------------------------------------------------------------------
     // Calculate electron-jet pair mass values
     //--------------------------------------------------------------------------
     
     double M_ej_avg, M_ej_min, M_ej_max;
     
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
     
     double sT_zjj = Pt_e1e2 + Jet1_Pt + Jet2_Pt;

     //--------------------------------------------------------------------------
     // Fill electron variables 
     //--------------------------------------------------------------------------
     
     if ( nEle_store >= 1 ) fillVariableWithValue( "Ele1_Pt", Ele1_Pt, gen_weight * pileup_weight  ) ;
     if ( nEle_store >= 2 ) fillVariableWithValue( "Ele2_Pt", Ele2_Pt, gen_weight * pileup_weight  ) ;
			
     //--------------------------------------------------------------------------
     // Fill jet variables 
     //--------------------------------------------------------------------------
					    
     // Jets								    
     fillVariableWithValue("nJet", nJet_ptCut, gen_weight * pileup_weight );
     if ( nJet_store >= 1 ) { 						    
       fillVariableWithValue( "Jet1_Pt"    , Jet1_Pt     , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "Jet1_Eta"   , Jet1_Eta    , gen_weight * pileup_weight  ) ;
     }

     //--------------------------------------------------------------------------
     // Fill DeltaR variables
     //--------------------------------------------------------------------------

     if ( nEle_store >= 2 && nJet_store >= 1) {
       fillVariableWithValue( "DR_Ele1Jet1"  , DR_Ele1Jet1 , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "DR_Ele2Jet1"  , DR_Ele2Jet1 , gen_weight * pileup_weight  ) ;
     }
     
     //--------------------------------------------------------------------------
     // Multi-object variables
     //--------------------------------------------------------------------------

     double sT_eej = Ele1_Pt + Ele2_Pt + Jet1_Pt;

     if ( nEle_store >= 2 ) { 						    
       fillVariableWithValue( "M_e1e2"     , M_e1e2 , gen_weight * pileup_weight  ) ;
       if ( nJet_store >= 1 ) { 
	 fillVariableWithValue( "sT_eej"     , sT_eej  , gen_weight * pileup_weight  ) ;
       }      
     }

     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Did we at least pass the noise filtes?
     //--------------------------------------------------------------------------
     
     bool passed_minimum = ( passedAllPreviousCuts("PassTrackingFailure") && passedCut ("PassTrackingFailure"));
     
     //--------------------------------------------------------------------------
     // Did we pass preselection?
     //--------------------------------------------------------------------------
     
     bool passed_preselection = ( passedAllPreviousCuts("M_e1e2") && passedCut ("M_e1e2") );

     //--------------------------------------------------------------------------
     // Are we in the region of interest?
     //--------------------------------------------------------------------------

     bool passed_region_of_interest = bool ( passed_preselection && M_e1e2 > 170. && sT_eej > 600.0 );

     //--------------------------------------------------------------------------
     // Fill plots with no selection applied
     //--------------------------------------------------------------------------
     
     FillUserTH1D( "PileupWeight"   , pileup_weight );
     FillUserTH1D( "GeneratorWeight", gen_weight ) ;
     
     //--------------------------------------------------------------------------
     // Fill noise filter level plots
     //--------------------------------------------------------------------------
     
     if ( passed_minimum && isData ){ 
       FillUserTH1D ("run_HLT", run );
       profile_run_vs_nvtx_HLT -> Fill ( run, nVertex, 1 ) ;
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

       TLorentzVector e1, j1, e2, mu, met;
       TLorentzVector e1e2mu;
       TLorentzVector eej, ee;
       
       e1.SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
       e2.SetPtEtaPhiM ( Ele2_Pt, Ele2_Eta, Ele2_Phi, 0.0 );
       j1.SetPtEtaPhiM ( Jet1_Pt, Jet1_Eta, Jet1_Phi, 0.0 );
       mu.SetPtEtaPhiM ( Muon1_Pt, Muon1_Eta, Muon1_Phi, 0.0 );
       met.SetPtEtaPhiM ( PFMET_Type01XY_Pt, 0.0, PFMET_Type01XY_Phi, 0.0 );

       eej  = e1 + e2 + j1;
       ee   = e1 + e2;
       
       double DR_Ele1Ele2 = e1.DeltaR( e2 );
       double M_eej  = eej.M();
       double DR_ZJ1 = ee.DeltaR ( j1 );
       
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
       // Preselection histograms
       //--------------------------------------------------------------------------
       
       FillUserTH1D("ProcessID_PAS"         , ProcessID                      , pileup_weight * gen_weight );
       FillUserTH1D("DR_Ele1Ele2_PAS"	    , DR_Ele1Ele2                    , pileup_weight * gen_weight );
       FillUserTH1D("EleChargeSum_PAS"      , Ele1_Charge + Ele2_Charge      , pileup_weight * gen_weight );
       FillUserTH1D("nElectron_PAS"         , nEle_ptCut                     , pileup_weight * gen_weight );
       FillUserTH1D("nMuon_PAS"             , nMuon_ptCut                    , pileup_weight * gen_weight );
       FillUserTH1D("nJet_PAS"              , nJet_ptCut                     , pileup_weight * gen_weight );
       FillUserTH1D("Pt1stEle_PAS"	    , Ele1_Pt                        , pileup_weight * gen_weight );
       FillUserTH1D("Eta1stEle_PAS"	    , Ele1_Eta                       , pileup_weight * gen_weight );
       FillUserTH1D("Phi1stEle_PAS"	    , Ele1_Phi                       , pileup_weight * gen_weight );
       FillUserTH1D("Pt2ndEle_PAS"	    , Ele2_Pt                        , pileup_weight * gen_weight );
       FillUserTH1D("Eta2ndEle_PAS"	    , Ele2_Eta                       , pileup_weight * gen_weight );
       FillUserTH1D("Phi2ndEle_PAS"	    , Ele2_Phi                       , pileup_weight * gen_weight );
       FillUserTH1D("Charge1stEle_PAS"	    , Ele1_Charge                    , pileup_weight * gen_weight );
       FillUserTH1D("Charge2ndEle_PAS"	    , Ele2_Charge                    , pileup_weight * gen_weight );
       FillUserTH1D("MET_PAS"               , PFMET_Type01XY_Pt              , pileup_weight * gen_weight );
       FillUserTH1D("METPhi_PAS"	    , PFMET_Type01XY_Phi             , pileup_weight * gen_weight );
       FillUserTH1D("Pt1stJet_PAS"          , Jet1_Pt                        , pileup_weight * gen_weight );
       FillUserTH1D("Eta1stJet_PAS"         , Jet1_Eta                       , pileup_weight * gen_weight );
       FillUserTH1D("Phi1stJet_PAS"	    , Jet1_Phi                       , pileup_weight * gen_weight );
       FillUserTH1D("sTlep_PAS"             , Ele1_Pt + Ele2_Pt              , pileup_weight * gen_weight );
       FillUserTH1D("sT_PAS"                , sT_eej                         , pileup_weight * gen_weight );
       FillUserTH1D("Mee_PAS"		    , M_e1e2                         , pileup_weight * gen_weight );
       FillUserTH1D( "MTenu_PAS"            , MT_Ele1MET                     , pileup_weight * gen_weight );
       FillUserTH1D("Me1j1_PAS"             , M_e1j1                         , pileup_weight * gen_weight );
       FillUserTH1D("Me2j1_PAS"             , M_e2j1                         , pileup_weight * gen_weight );
       FillUserTH1D("Ptee_PAS"              , Pt_e1e2                        , pileup_weight * gen_weight );
       FillUserTH1D("DCotTheta1stEle_PAS"   , Ele1_DCotTheta                 , pileup_weight * gen_weight );
       FillUserTH1D("Dist1stEle_PAS"        , Ele1_Dist                      , pileup_weight * gen_weight );
       FillUserTH1D("DCotTheta2ndEle_PAS"   , Ele2_DCotTheta                 , pileup_weight * gen_weight );
       FillUserTH1D("Dist2ndEle_PAS"        , Ele2_Dist                      , pileup_weight * gen_weight );
       FillUserTH1D("nVertex_PAS"           , nVertex                        , pileup_weight * gen_weight );
       FillUserTH1D("DR_Ele1Jet1_PAS"	    , DR_Ele1Jet1                    , pileup_weight * gen_weight );
       FillUserTH1D("DR_Ele2Jet1_PAS"	    , DR_Ele2Jet1                    , pileup_weight * gen_weight );
       FillUserTH1D("Meej_PAS"              , M_eej                          , pileup_weight * gen_weight );
       FillUserTH1D("DR_ZJet1_PAS"          , DR_ZJ1                         , pileup_weight * gen_weight );
       FillUserTH1D("sTfrac_Jet1_PAS"       , Jet1_Pt / sT_eej               , pileup_weight * gen_weight );
       FillUserTH1D("sTfrac_Ele1_PAS"       , Ele1_Pt / sT_eej               , pileup_weight * gen_weight );
       FillUserTH1D("sTfrac_Ele2_PAS"       , Ele2_Pt / sT_eej               , pileup_weight * gen_weight );
       FillUserTH1D("sTfrac_Jet_PAS"        , ( Jet1_Pt ) / sT_eej           , pileup_weight * gen_weight );
       FillUserTH1D("sTfrac_Ele_PAS"        , ( Ele1_Pt + Ele2_Pt ) / sT_eej , pileup_weight * gen_weight );

       FillUserTH2D("MeeVsST_PAS" , M_e1e2, sT_eej, pileup_weight * gen_weight );	   
       
       //--------------------------------------------------------------------------
       // Preselection + data-only
       //--------------------------------------------------------------------------
       
       if ( isData == 1 ) { 
	 FillUserTH1D("split_PAS"    , get_split     ( run ), pileup_weight * gen_weight ) ;
	 FillUserTH1D("split_1fb_PAS", get_split_1fb ( run ), pileup_weight * gen_weight ) ;
	 FillUserTH1D("run_PAS"  , run );
	 profile_run_vs_nvtx_PAS -> Fill ( run, nVertex, 1 );
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

       if ( sT_eej  > 445. ) FillUserTH1D("Mee_PASandST445"   , M_e1e2  , pileup_weight * gen_weight ) ;

       //--------------------------------------------------------------------------
       // High M(ee) plots
       //--------------------------------------------------------------------------

       if ( M_e1e2  > 100. ) FillUserTH1D("sT_PASandMee100"   , sT_eej , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 110. ) FillUserTH1D("sT_PASandMee110"   , sT_eej , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 120. ) FillUserTH1D("sT_PASandMee120"   , sT_eej , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 130. ) FillUserTH1D("sT_PASandMee130"   , sT_eej , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 140. ) FillUserTH1D("sT_PASandMee140"   , sT_eej , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 150. ) FillUserTH1D("sT_PASandMee150"   , sT_eej , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 160. ) FillUserTH1D("sT_PASandMee160"   , sT_eej , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 170. ) FillUserTH1D("sT_PASandMee170"   , sT_eej , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 180. ) FillUserTH1D("sT_PASandMee180"   , sT_eej , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 190. ) FillUserTH1D("sT_PASandMee190"   , sT_eej , pileup_weight * gen_weight ) ;
       if ( M_e1e2  > 200. ) FillUserTH1D("sT_PASandMee200"   , sT_eej , pileup_weight * gen_weight ) ;
       
       if ( M_e1e2 > 100. ) { 
	 if ( isData == 1 ) FillUserTH1D("split_PASandMee100", get_split ( run ), pileup_weight * gen_weight ) ;
	 FillUserTH1D("Ptee_PASandMee100"              , Pt_e1e2                        , pileup_weight * gen_weight );
	 FillUserTH2D("MeeVsST_PASandMee100"           , M_e1e2, sT_eej                 , pileup_weight * gen_weight );	   
	 FillUserTH1D("nVertex_PASandMee100"           , nVertex                        , pileup_weight * gen_weight );
	 FillUserTH1D("EleChargeSum_PASandMee100"      , Ele1_Charge + Ele2_Charge      , pileup_weight * gen_weight ) ;
	 FillUserTH1D("nJet_PASandMee100"              , nJet_ptCut                     , pileup_weight * gen_weight ) ;
	 FillUserTH1D("sTlep_PASandMee100"             , Ele1_Pt + Ele2_Pt              , pileup_weight * gen_weight ) ;
	 FillUserTH1D("Me1j1_PASandMee100"             , M_e1j1                         , pileup_weight * gen_weight ) ;
	 FillUserTH1D("Pt1stEle_PASandMee100"          , Ele1_Pt                        , pileup_weight * gen_weight ) ;
	 FillUserTH1D("Pt2ndEle_PASandMee100"          , Ele2_Pt                        , pileup_weight * gen_weight ) ;
	 FillUserTH1D("Pt1stJet_PASandMee100"          , Jet1_Pt                        , pileup_weight * gen_weight ) ;
	 FillUserTH1D("sTfrac_Jet1_PASandMee100"       , Jet1_Pt / sT_eej               , pileup_weight * gen_weight ) ;
	 FillUserTH1D("sTfrac_Ele1_PASandMee100"       , Ele1_Pt / sT_eej               , pileup_weight * gen_weight ) ;
	 FillUserTH1D("sTfrac_Ele2_PASandMee100"       , Ele2_Pt / sT_eej               , pileup_weight * gen_weight ) ;
	 FillUserTH1D("sTfrac_Jet_PASandMee100"        , ( Jet1_Pt )           / sT_eej , pileup_weight * gen_weight ) ;
	 FillUserTH1D("sTfrac_Ele_PASandMee100"        , ( Ele1_Pt + Ele2_Pt ) / sT_eej , pileup_weight * gen_weight ) ;
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
	 if ( sT_eej > 600 ) 	 FillUserTH1D("Mee_70_110_ST600_Preselection", M_e1e2, pileup_weight * gen_weight ); 
	 if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_70_110_PAS", M_e1e2, pileup_weight * gen_weight ); 
	 else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_70_110_PAS", M_e1e2, pileup_weight * gen_weight ); 
	 else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_70_110_PAS", M_e1e2, pileup_weight * gen_weight ); 
	 if      ( isEB   ) FillUserTH1D( "Mee_EB_70_110_PAS"  , M_e1e2, pileup_weight * gen_weight ); 

	 FillUserTH1D( "ProcessID_ZWindow", ProcessID, pileup_weight * gen_weight );
	 if ( ProcessID == 0 ) FillUserTH1D ( "Mee_70_110_Preselection_Process0", M_e1e2, pileup_weight * gen_weight );
	 if ( ProcessID == 1 ) FillUserTH1D ( "Mee_70_110_Preselection_Process1", M_e1e2, pileup_weight * gen_weight );
	 if ( ProcessID == 2 ) FillUserTH1D ( "Mee_70_110_Preselection_Process2", M_e1e2, pileup_weight * gen_weight );
	 if ( ProcessID == 3 ) FillUserTH1D ( "Mee_70_110_Preselection_Process3", M_e1e2, pileup_weight * gen_weight );
	 if ( ProcessID == 4 ) FillUserTH1D ( "Mee_70_110_Preselection_Process4", M_e1e2, pileup_weight * gen_weight );
       }

       //--------------------------------------------------------------------------
       // Region of interest plots
       //-------------------------------------------------------------------------- 

       if ( passed_region_of_interest ) { 

	 std::cout << "*** " << int(run) << ":" << int(ls) << ":" << int (event) << std::endl;

	 if ( isData == 1 ) FillUserTH1D( "split_ROI" , get_split_1fb ( run ) , pileup_weight * gen_weight );
	 
	 FillUserTH1D("Me1j1_ROI"           , M_e1j1                         , pileup_weight * gen_weight );
	 FillUserTH1D("Ptee_ROI"            , Pt_e1e2                        , pileup_weight * gen_weight );
	 FillUserTH1D("Eta1stJet_ROI"       , Jet1_Eta                       , pileup_weight * gen_weight );
	 FillUserTH1D("Eta1stEle_ROI"	    , Ele1_Eta                       , pileup_weight * gen_weight );
	 FillUserTH1D("Eta2ndEle_ROI"	    , Ele2_Eta                       , pileup_weight * gen_weight );
	 FillUserTH1D("Phi1stJet_ROI"       , Jet1_Phi                       , pileup_weight * gen_weight );
	 FillUserTH1D("Phi1stEle_ROI"	    , Ele1_Phi                       , pileup_weight * gen_weight );
	 FillUserTH1D("Phi2ndEle_ROI"	    , Ele2_Phi                       , pileup_weight * gen_weight );
	 FillUserTH2D("MeeVsST_ROI"         , M_e1e2, sT_eej                 , pileup_weight * gen_weight );	   
	 FillUserTH1D("Mee_ROI"		    , M_e1e2                         , pileup_weight * gen_weight );
	 FillUserTH1D("nVertex_ROI"         , nVertex                        , pileup_weight * gen_weight );
	 FillUserTH1D("nJet_ROI"            , nJet_ptCut                     , pileup_weight * gen_weight );
	 FillUserTH1D("EleChargeSum_ROI"    , Ele1_Charge + Ele2_Charge      , pileup_weight * gen_weight );
	 FillUserTH1D("Meej_ROI"            , M_eej                          , pileup_weight * gen_weight );
	 FillUserTH1D("DR_ZJet1_ROI"        , DR_ZJ1                         , pileup_weight * gen_weight );
	 FillUserTH1D("MET_ROI"             , PFMET_Type01XY_Pt              , pileup_weight * gen_weight );
	 FillUserTH1D("sT_ROI"              , sT_eej                         , pileup_weight * gen_weight );
	 FillUserTH1D("sTlep_ROI"           , Ele1_Pt + Ele2_Pt              , pileup_weight * gen_weight );
	 FillUserTH1D("Pt1stEle_ROI"        , Ele1_Pt                        , pileup_weight * gen_weight );
	 FillUserTH1D("Pt2ndEle_ROI"        , Ele2_Pt                        , pileup_weight * gen_weight );
	 FillUserTH1D("Pt1stJet_ROI"        , Jet1_Pt                        , pileup_weight * gen_weight );
	 FillUserTH1D( "sTfrac_Jet1_ROI"    , Jet1_Pt / sT_eej               , pileup_weight * gen_weight );
	 FillUserTH1D( "sTfrac_Ele1_ROI"    , Ele1_Pt / sT_eej               , pileup_weight * gen_weight );
	 FillUserTH1D( "sTfrac_Ele2_ROI"    , Ele2_Pt / sT_eej               , pileup_weight * gen_weight );
	 FillUserTH1D( "sTfrac_Jet_ROI"     , ( Jet1_Pt )           / sT_eej , pileup_weight * gen_weight );
	 FillUserTH1D( "sTfrac_Ele_ROI"     , ( Ele1_Pt + Ele2_Pt ) / sT_eej , pileup_weight * gen_weight );

       }

     } // End preselection 
   } // End loop over events
   
   output_root_ -> cd();
   profile_run_vs_nvtx_HLT -> Write();
   profile_run_vs_nvtx_PAS -> Write();
   
   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

