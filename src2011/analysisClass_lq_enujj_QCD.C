#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile){}

analysisClass::~analysisClass(){}

void analysisClass::Loop()
{
   std::cout << "analysisClass::Loop() begins" <<std::endl;   

   //--------------------------------------------------------------------------
   // Decide which plots to save (default is to save everything)
   //--------------------------------------------------------------------------
   
   fillSkim                         (  true  ) ;
   fillAllPreviousCuts              (  true  ) ;
   fillAllOtherCuts                 ( !true  ) ; 
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------

   CreateUserTH1D( "Total_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_nMuon_PAS"             ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_nMuon_PAS"             ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_Pt1stEle_PAS"	         ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_Pt1stEle_PAS"	        ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_Eta1stEle_PAS"	 ,    100 , -5      , 5	       );  CreateUserTH1D( "Pass_Eta1stEle_PAS"	        ,    100 , -5      , 5	      );
   CreateUserTH1D( "Total_Phi1stEle_PAS"	 ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_Phi1stEle_PAS"	        ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_Charge1stEle_PAS"	 ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_Charge1stEle_PAS"	,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_MET_PAS"               ,    200 , 0       , 1000     );  CreateUserTH1D( "Pass_MET_PAS"               ,    200 , 0       , 1000     );
   CreateUserTH1D( "Total_METPhi_PAS"		 ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_METPhi_PAS"		,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_minMETPt1stEle_PAS"    ,    200 , 0       , 1000     );  CreateUserTH1D( "Pass_minMETPt1stEle_PAS"    ,    200 , 0       , 1000     );
   CreateUserTH1D( "Total_Pt1stJet_PAS"          ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_Pt1stJet_PAS"          ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_Pt2ndJet_PAS"          ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_Pt2ndJet_PAS"          ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_Eta1stJet_PAS"         ,    100 , -5      , 5	       );  CreateUserTH1D( "Pass_Eta1stJet_PAS"         ,    100 , -5      , 5	      );
   CreateUserTH1D( "Total_Eta2ndJet_PAS"         ,    100 , -5      , 5	       );  CreateUserTH1D( "Pass_Eta2ndJet_PAS"         ,    100 , -5      , 5	      );
   CreateUserTH1D( "Total_Phi1stJet_PAS"	 ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_Phi1stJet_PAS"	        ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_Phi2ndJet_PAS"	 ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_Phi2ndJet_PAS"	        ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_TCHE1stJet_PAS"        ,    100 , 0       , 20       );  CreateUserTH1D( "Pass_TCHE1stJet_PAS"        ,    100 , 0       , 20       );
   CreateUserTH1D( "Total_TCHE2ndJet_PAS"        ,    100 , 0       , 20       );  CreateUserTH1D( "Pass_TCHE2ndJet_PAS"        ,    100 , 0       , 20       );
   CreateUserTH1D( "Total_nMuon_PtCut_IDISO_PAS" ,    16  , -0.5    , 15.5     );  CreateUserTH1D( "Pass_nMuon_PtCut_IDISO_PAS" ,    16  , -0.5    , 15.5     );
   CreateUserTH1D( "Total_MTenu_PAS"             ,    200 , 0       , 1000     );  CreateUserTH1D( "Pass_MTenu_PAS"             ,    200 , 0       , 1000     );
   CreateUserTH1D( "Total_Ptenu_PAS"		 ,    200 , 0       , 2000     );  CreateUserTH1D( "Pass_Ptenu_PAS"		,    200 , 0       , 2000     );
   CreateUserTH1D( "Total_sTlep_PAS"             ,    200 , 0       , 2000     );  CreateUserTH1D( "Pass_sTlep_PAS"             ,    200 , 0       , 2000     );
   CreateUserTH1D( "Total_sTjet_PAS"             ,    200 , 0       , 2000     );  CreateUserTH1D( "Pass_sTjet_PAS"             ,    200 , 0       , 2000     );
   CreateUserTH1D( "Total_sT_PAS"                ,    200 , 0       , 2000     );  CreateUserTH1D( "Pass_sT_PAS"                ,    200 , 0       , 2000     );
   CreateUserTH1D( "Total_Mjj_PAS"		 ,    200 , 0       , 2000     );  CreateUserTH1D( "Pass_Mjj_PAS"		,    200 , 0       , 2000     );
   CreateUserTH1D( "Total_Mej_1stPair_PAS"       ,    200 , 0       , 2000     );  CreateUserTH1D( "Pass_Mej_1stPair_PAS"       ,    200 , 0       , 2000     );
   CreateUserTH1D( "Total_Mej_2ndPair_PAS"       ,    200 , 0       , 2000     );  CreateUserTH1D( "Pass_Mej_2ndPair_PAS"       ,    200 , 0       , 2000     );
   CreateUserTH1D( "Total_HcalIso1stEle_PAS"     ,    200 , 0       , 20       );  CreateUserTH1D( "Pass_HcalIso1stEle_PAS"     ,    200 , 0       , 20       );
   CreateUserTH1D( "Total_EcalIso1stEle_PAS"     ,    200 , 0       , 20       );  CreateUserTH1D( "Pass_EcalIso1stEle_PAS"     ,    200 , 0       , 20       );
   CreateUserTH1D( "Total_RelIso1stEle_PAS"      ,    200 , 0       , 1.0      );  CreateUserTH1D( "Pass_RelIso1stEle_PAS"      ,    200 , 0       , 1.0      );
   CreateUserTH1D( "Total_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_mDPhi1stEleMET"        ,    100 , 0.      , 3.14159  );  CreateUserTH1D( "Pass_mDPhi1stEleMET"        ,    100 , 0.      , 3.14159  );
   CreateUserTH1D( "Total_mDPhi1stJetMET"        ,    100 , 0.      , 3.14159  );  CreateUserTH1D( "Pass_mDPhi1stJetMET"        ,    100 , 0.      , 3.14159  );
   CreateUserTH1D( "Total_mDPhi2ndJetMET"        ,    100 , 0.      , 3.14159  );  CreateUserTH1D( "Pass_mDPhi2ndJetMET"        ,    100 , 0.      , 3.14159  );
   CreateUserTH1D( "Total_MT_GoodVtxLTE5"        ,    200 , 0.      , 1000     );  CreateUserTH1D( "Pass_MT_GoodVtxLTE5"        ,    200 , 0.      , 1000     );
   CreateUserTH1D( "Total_MT_GoodVtxGT5"         ,    200 , 0.      , 1000     );  CreateUserTH1D( "Pass_MT_GoodVtxGT5"         ,    200 , 0.      , 1000     );

   //--------------------------------------------------------------------------
   // Get preselection values
   //--------------------------------------------------------------------------
   
   //--------------------------------------------------------------------------
   // Loop over the chain
   //--------------------------------------------------------------------------

   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntries();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%10000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << "/" << nentries << std::endl;   

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
     
     double weight     = getPileupWeight ( nPileUpInteractions, isData ) ;

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "PassJSON"                      , passedJSON    ); 
     
     // Filters
     fillVariableWithValue(   "PassHBHENoiseFilter"           , PassHBHENoiseFilter ) ;
     fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight ) ;
									      
     // Muon variables ( for veto ) 					      
     fillVariableWithValue(   "nMuon_PtCut_ID_ISO"            , nMuon_Ana     );
			                                      		      
     // 1st Electron variables				      		      
     fillVariableWithValue(   "nEle_PtCut_IDISO_noOvrlp"      , nEle_Ana      ); 
     fillVariableWithValue(   "Pt1stEle_PtCut_IDISO_noOvrlp"  , Ele1_Pt       );
     fillVariableWithValue(   "Eta1stEle_PtCut_IDISO_noOvrlp" , Ele1_Eta      );

     // 1st JET variables                                     
     fillVariableWithValue(   "nJet_PtCut_ID_noOvrlp"         , nJet_Ana      );

     // Dummy variables
     fillVariableWithValue ("denominator",1);

     // Trigger 
     
     int passTrigger = 1;
     double prescale = 1.0;

     fillVariableWithValue ( "Photon1Pt" , Photon1_Pt  );
     fillVariableWithValue ( "Photon1HoE", Photon1_HoE );
     
     if ( isData ) { 
       std::vector <std::string> pass_names;
       std::vector <int>         pass_prescales;
       std::vector <int>         pass_thresholds;
       
       if ( HLT_Photon30_CaloIdVL_v1 > -500. ) { pass_names.push_back ( "HLT_Photon30_CaloIdVL_v1" ); pass_prescales.push_back ( HLT_Photon30_CaloIdVL_v1); pass_thresholds.push_back (30); } 
       // if ( HLT_Photon30_CaloIdVL_v2 > -500. ) { pass_names.push_back ( "HLT_Photon30_CaloIdVL_v2" ); pass_prescales.push_back ( HLT_Photon30_CaloIdVL_v2); pass_thresholds.push_back (30); } 
       // if ( HLT_Photon30_CaloIdVL_v3 > -500. ) { pass_names.push_back ( "HLT_Photon30_CaloIdVL_v3" ); pass_prescales.push_back ( HLT_Photon30_CaloIdVL_v3); pass_thresholds.push_back (30); } 
       // if ( HLT_Photon30_CaloIdVL_v4 > -500. ) { pass_names.push_back ( "HLT_Photon30_CaloIdVL_v4" ); pass_prescales.push_back ( HLT_Photon30_CaloIdVL_v4); pass_thresholds.push_back (30); } 
       // if ( HLT_Photon30_CaloIdVL_v5 > -500. ) { pass_names.push_back ( "HLT_Photon30_CaloIdVL_v5" ); pass_prescales.push_back ( HLT_Photon30_CaloIdVL_v5); pass_thresholds.push_back (30); } 
       if ( HLT_Photon50_CaloIdVL_v1 > -500. ) { pass_names.push_back ( "HLT_Photon50_CaloIdVL_v1" ); pass_prescales.push_back ( HLT_Photon50_CaloIdVL_v1); pass_thresholds.push_back (50); } 
       // if ( HLT_Photon50_CaloIdVL_v2 > -500. ) { pass_names.push_back ( "HLT_Photon50_CaloIdVL_v2" ); pass_prescales.push_back ( HLT_Photon50_CaloIdVL_v2); pass_thresholds.push_back (50); } 
       // if ( HLT_Photon75_CaloIdVL_v1 > -500. ) { pass_names.push_back ( "HLT_Photon75_CaloIdVL_v1" ); pass_prescales.push_back ( HLT_Photon75_CaloIdVL_v1); pass_thresholds.push_back (75); } 
       // if ( HLT_Photon75_CaloIdVL_v2 > -500. ) { pass_names.push_back ( "HLT_Photon75_CaloIdVL_v2" ); pass_prescales.push_back ( HLT_Photon75_CaloIdVL_v2); pass_thresholds.push_back (75); } 
       // if ( HLT_Photon75_CaloIdVL_v3 > -500. ) { pass_names.push_back ( "HLT_Photon75_CaloIdVL_v3" ); pass_prescales.push_back ( HLT_Photon75_CaloIdVL_v3); pass_thresholds.push_back (75); } 
       // if ( HLT_Photon75_CaloIdVL_v4 > -500. ) { pass_names.push_back ( "HLT_Photon75_CaloIdVL_v4" ); pass_prescales.push_back ( HLT_Photon75_CaloIdVL_v4); pass_thresholds.push_back (75); } 
       // if ( HLT_Photon75_CaloIdVL_v5 > -500. ) { pass_names.push_back ( "HLT_Photon75_CaloIdVL_v5" ); pass_prescales.push_back ( HLT_Photon75_CaloIdVL_v5); pass_thresholds.push_back (75); } 
       // if ( HLT_Photon90_CaloIdVL_v1 > -500. ) { pass_names.push_back ( "HLT_Photon90_CaloIdVL_v1" ); pass_prescales.push_back ( HLT_Photon90_CaloIdVL_v1); pass_thresholds.push_back (90); } 
       // if ( HLT_Photon90_CaloIdVL_v2 > -500. ) { pass_names.push_back ( "HLT_Photon90_CaloIdVL_v2" ); pass_prescales.push_back ( HLT_Photon90_CaloIdVL_v2); pass_thresholds.push_back (90); } 

       int n_pass = (int) pass_names.size();

       passTrigger = ( n_pass >= 1 ) ? 1 : 0;
       
       if ( n_pass > 1 ) {
	 
	 int n_prescaled = 0;
	 int n_unprescaled = 0;
	 for (int i = 0; i < n_pass ; ++i){
	   if ( pass_prescales[i] ==  1 ) n_unprescaled++;
	   if ( pass_prescales[i] >   1 ) n_prescaled++;
	 }

	 if ( n_unprescaled  > 0 ) prescale = 1.0;
	 
	 if ( n_unprescaled == 0 ) { 
	   
	   double pt = Photon1_Pt;
	   int nPhoton = (int) nPhoton_Ana;
	   if ( nPhoton_Ana == 0 ) pt = Ele1_Pt;
	   
	   int i_max_trigger_under_pt = -1;
	   int i_min_trigger_over_pt  = -1;

	   double max_threshold_under_pt = -1;
	   double min_threshold_over_pt = 9999;

	   for (int i = 0; i < n_pass ; ++i){
	     double threshold = pass_thresholds[i];
	     if ( pt < (double) threshold ) continue;
	     if ( threshold > max_threshold_under_pt ) {
	       max_threshold_under_pt = threshold;
	       i_max_trigger_under_pt = i;
	     }
	   }

	   for (int i = 0; i < n_pass ; ++i){
	     double threshold = pass_thresholds[i];
	     if ( pt > (double) threshold ) continue;
	     if ( threshold < min_threshold_over_pt ) {
	       min_threshold_over_pt = threshold;
	       i_min_trigger_over_pt = i;
	       prescale = pass_prescales[i];
	     }
	   }

	   if ( i_max_trigger_under_pt < 0 ) {
	     i_max_trigger_under_pt = 0;
	     max_threshold_under_pt = pass_thresholds[0];
	     prescale = (double) pass_prescales[0];
	   }

	   /*
	   std::cout << "------------------------------------------------------" << std::endl;
	   std::cout << "NPhoton          = " << int(nPhoton_Stored) << std::endl;
	   std::cout << "Lead photon   pt = " << Photon1_Pt << std::endl;
	   std::cout << "Lead electron pt = " << Ele1_Pt    << std::endl;
	   std::cout << "PT               = " << pt << std::endl;
	   std::cout << "Triggers with the following thresholds fired:";
	   for (int i = 0; i < n_pass ; ++i){
	     std::cout << " " << pass_thresholds[i];
	   }
	   std::cout << std::endl;
	   std::cout << "Triggers with the following prescales  fired:";
	   for (int i = 0; i < n_pass ; ++i){
	     std::cout << " " << pass_prescales[i];
	   }
	   std::cout << std::endl;
	   std::cout << "Max trigger      = " << pass_names[i_max_trigger_under_pt] << ", " << pass_prescales [i_max_trigger_under_pt] << std::endl;
	   */
	 }
       }
       else if ( n_pass == 1 ) { 
	 prescale = pass_prescales[0];
       }
     }

     weight *= prescale;

     fillVariableWithValue ("PassTrigger", passTrigger ) ;

     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------
     
     bool passed_denominator = passedAllPreviousCuts("denominator");
       
     if ( passed_denominator ) { 
       
        FillUserTH1D( "Total_nElectron_PAS"         , nEle_Ana                       , weight);
        FillUserTH1D( "Total_nMuon_PAS"             , nMuon_Ana                      , weight);
        FillUserTH1D( "Total_Pt1stEle_PAS"	    , Ele1_Pt                        , weight);
       	FillUserTH1D( "Total_Eta1stEle_PAS"	    , Ele1_Eta                       , weight);
       	FillUserTH1D( "Total_Phi1stEle_PAS"	    , Ele1_Phi                       , weight);
        FillUserTH1D( "Total_Charge1stEle_PAS"      , Ele1_Charge                    , weight);  
       	FillUserTH1D( "Total_MET_PAS"               , MET_Pt                         , weight);
        FillUserTH1D( "Total_METPhi_PAS"	    , MET_Phi                        , weight);  
       	FillUserTH1D( "Total_minMETPt1stEle_PAS"    , TMath::Min ( Ele1_Pt, MET_Pt  ), weight);
       	FillUserTH1D( "Total_Pt1stJet_PAS"          , Jet1_Pt                        , weight);
       	FillUserTH1D( "Total_Pt2ndJet_PAS"          , Jet2_Pt                        , weight);
       	FillUserTH1D( "Total_Eta1stJet_PAS"         , Jet1_Eta                       , weight);
       	FillUserTH1D( "Total_Eta2ndJet_PAS"         , Jet2_Eta                       , weight);
       	FillUserTH1D( "Total_Phi1stJet_PAS"         , Jet1_Phi                       , weight);
       	FillUserTH1D( "Total_Phi2ndJet_PAS"	    , Jet2_Phi                       , weight);
       	FillUserTH1D( "Total_TCHE1stJet_PAS"        , Jet1_btagTCHE                  , weight);
       	FillUserTH1D( "Total_TCHE2ndJet_PAS"        , Jet2_btagTCHE                  , weight);
        FillUserTH1D( "Total_nMuon_PtCut_IDISO_PAS" , nMuon_Ana                      , weight);
       	FillUserTH1D( "Total_MTenu_PAS"             , MT_Ele1MET                     , weight);
       	FillUserTH1D( "Total_Ptenu_PAS"	            , Pt_Ele1MET                     , weight);
       	FillUserTH1D( "Total_sTlep_PAS"             , Ele1_Pt + MET_Pt               , weight);
       	FillUserTH1D( "Total_sTjet_PAS"             , Jet1_Pt + Jet2_Pt              , weight);
       	FillUserTH1D( "Total_sT_PAS"                , sT_enujj                       , weight);
        FillUserTH1D( "Total_Mjj_PAS"	            , M_j1j2                         , weight);  
       	FillUserTH1D( "Total_Mej_1stPair_PAS"       , M_ej1                          , weight);
       	FillUserTH1D( "Total_Mej_2ndPair_PAS"       , M_ej2                          , weight);
       	FillUserTH1D( "Total_HcalIso1stEle_PAS"     , Ele1_HcalIso                   , weight);
       	FillUserTH1D( "Total_EcalIso1stEle_PAS"     , Ele1_EcalIso                   , weight);
       	FillUserTH1D( "Total_RelIso1stEle_PAS"      , Ele1_RelIso                    , weight);
       	FillUserTH1D( "Total_DCotTheta1stEle_PAS"   , Ele1_DCotTheta                 , weight);
       	FillUserTH1D( "Total_Dist1stEle_PAS"        , Ele1_Dist                      , weight);
       	FillUserTH1D( "Total_mDPhi1stEleMET"        , mDPhi_METEle1                  , weight);
       	FillUserTH1D( "Total_mDPhi1stJetMET"        , mDPhi_METJet1                  , weight);
       	FillUserTH1D( "Total_mDPhi2ndJetMET"        , mDPhi_METJet2                  , weight);
	
	if ( nVertex_good <= 5 ) FillUserTH1D ( "Total_MT_GoodVtxLTE5", MT_Ele1MET , weight);
	if ( nVertex_good >  5 ) FillUserTH1D ( "Total_MT_GoodVtxGT5" , MT_Ele1MET , weight);
	
	if ( Ele1_PassHEEP ) { 
	  
	  FillUserTH1D( "Pass_nElectron_PAS"         , nEle_Ana                       , weight);
	  FillUserTH1D( "Pass_nMuon_PAS"             , nMuon_Ana                      , weight);
	  FillUserTH1D( "Pass_Pt1stEle_PAS"	     , Ele1_Pt                        , weight);
	  FillUserTH1D( "Pass_Eta1stEle_PAS"	     , Ele1_Eta                       , weight);
	  FillUserTH1D( "Pass_Phi1stEle_PAS"	     , Ele1_Phi                       , weight);
	  FillUserTH1D( "Pass_Charge1stEle_PAS"      , Ele1_Charge                    , weight);
	  FillUserTH1D( "Pass_MET_PAS"               , MET_Pt                         , weight);
	  FillUserTH1D( "Pass_METPhi_PAS"	     , MET_Phi                        , weight);
	  FillUserTH1D( "Pass_minMETPt1stEle_PAS"    , TMath::Min ( Ele1_Pt, MET_Pt  ), weight);
	  FillUserTH1D( "Pass_Pt1stJet_PAS"          , Jet1_Pt                        , weight);
	  FillUserTH1D( "Pass_Pt2ndJet_PAS"          , Jet2_Pt                        , weight);
	  FillUserTH1D( "Pass_Eta1stJet_PAS"         , Jet1_Eta                       , weight);
	  FillUserTH1D( "Pass_Eta2ndJet_PAS"         , Jet2_Eta                       , weight);
	  FillUserTH1D( "Pass_Phi1stJet_PAS"         , Jet1_Phi                       , weight);
	  FillUserTH1D( "Pass_Phi2ndJet_PAS"	     , Jet2_Phi                       , weight);
	  FillUserTH1D( "Pass_TCHE1stJet_PAS"        , Jet1_btagTCHE                  , weight);
	  FillUserTH1D( "Pass_TCHE2ndJet_PAS"        , Jet2_btagTCHE                  , weight);
	  FillUserTH1D( "Pass_nMuon_PtCut_IDISO_PAS" , nMuon_Ana                      , weight);
	  FillUserTH1D( "Pass_MTenu_PAS"             , MT_Ele1MET                     , weight);
	  FillUserTH1D( "Pass_Ptenu_PAS"	     , Pt_Ele1MET                     , weight);
	  FillUserTH1D( "Pass_sTlep_PAS"             , Ele1_Pt + MET_Pt               , weight);
	  FillUserTH1D( "Pass_sTjet_PAS"             , Jet1_Pt + Jet2_Pt              , weight);
	  FillUserTH1D( "Pass_sT_PAS"                , sT_enujj                       , weight);
	  FillUserTH1D( "Pass_Mjj_PAS"	             , M_j1j2                         , weight);
	  FillUserTH1D( "Pass_Mej_1stPair_PAS"       , M_ej1                          , weight);
	  FillUserTH1D( "Pass_Mej_2ndPair_PAS"       , M_ej2                          , weight);
	  FillUserTH1D( "Pass_HcalIso1stEle_PAS"     , Ele1_HcalIso                   , weight);
	  FillUserTH1D( "Pass_EcalIso1stEle_PAS"     , Ele1_EcalIso                   , weight);
	  FillUserTH1D( "Pass_RelIso1stEle_PAS"      , Ele1_RelIso                    , weight);
	  FillUserTH1D( "Pass_DCotTheta1stEle_PAS"   , Ele1_DCotTheta                 , weight);
	  FillUserTH1D( "Pass_Dist1stEle_PAS"        , Ele1_Dist                      , weight);
	  FillUserTH1D( "Pass_mDPhi1stEleMET"        , mDPhi_METEle1                  , weight);
	  FillUserTH1D( "Pass_mDPhi1stJetMET"        , mDPhi_METJet1                  , weight);
	  FillUserTH1D( "Pass_mDPhi2ndJetMET"        , mDPhi_METJet2                  , weight);

	  if ( nVertex_good <= 5 ) FillUserTH1D ( "Pass_MT_GoodVtxLTE5", MT_Ele1MET , weight);
	  if ( nVertex_good >  5 ) FillUserTH1D ( "Pass_MT_GoodVtxGT5" , MT_Ele1MET , weight);
	 
	}
     }
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
