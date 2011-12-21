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

   // fake rates
   
   double fakeRate_low_bar_p0    = getPreCutValue1 ( "fakeRate_bar"  );
   double fakeRate_low_bar_p1    = getPreCutValue2 ( "fakeRate_bar"  );
   double fakeRate_high_bar_p0   = getPreCutValue3 ( "fakeRate_bar"  );
				 
   double fakeRate_low_end1_p0   = getPreCutValue1 ( "fakeRate_end1" );
   double fakeRate_low_end1_p1   = getPreCutValue2 ( "fakeRate_end1" );
   double fakeRate_high_end1_p0  = getPreCutValue3 ( "fakeRate_end1" );
				 
   double fakeRate_low_end2_p0   = getPreCutValue1 ( "fakeRate_end2" );
   double fakeRate_low_end2_p1   = getPreCutValue2 ( "fakeRate_end2" );
   double fakeRate_high_end2_p0  = getPreCutValue3 ( "fakeRate_end2" );
   
   double eFakeRate_low_bar_p0   = getPreCutValue1 ( "eFakeRate_bar" );
   double eFakeRate_low_bar_p1   = getPreCutValue2 ( "eFakeRate_bar" );
   double eFakeRate_high_bar_p0  = getPreCutValue3 ( "eFakeRate_bar" );

   double eFakeRate_low_end1_p0  = getPreCutValue1 ( "eFakeRate_end1");
   double eFakeRate_low_end1_p1  = getPreCutValue2 ( "eFakeRate_end1");
   double eFakeRate_high_end1_p0 = getPreCutValue3 ( "eFakeRate_end1");

   double eFakeRate_low_end2_p0  = getPreCutValue1 ( "eFakeRate_end2");
   double eFakeRate_low_end2_p1  = getPreCutValue2 ( "eFakeRate_end2");
   double eFakeRate_high_end2_p0 = getPreCutValue3 ( "eFakeRate_end2");

   // trigger requirements
    
   double trigger_tolerance = getPreCutValue1("trigger_tolerance"); 
   
   //--------------------------------------------------------------------------
   // Set global variables
   //--------------------------------------------------------------------------
   
   TVector2 v_METCharged, v_METType1, v_ele;

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------

   CreateUserTH1D( "nElectron_PAS"            , 5   , -0.5    , 4.5      );
   CreateUserTH1D( "nMuon_PAS"                , 5   , -0.5    , 4.5      );
   CreateUserTH1D( "Pt1stEle_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta1stEle_PAS"	      , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stEle_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Charge1stEle_PAS"	      , 2   , -1.0001 , 1.0001	 ); 
   CreateUserTH1D( "MatchPhotonConv1stEle_PAS",	2   , -0.5    , 1.5      );
   CreateUserTH1D( "MatchPhotonConv2ndEle_PAS",	2   , -0.5    , 1.5      );
   CreateUserTH1D( "MET_PAS"                  , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "METPhi_PAS"		      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "METCharged_PAS"           , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "METChargedPhi_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "METType1_PAS"             , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "METType1Phi_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "METSig_PAS"               , 100 , 0       , 800      );
   CreateUserTH1D( "minMETPt1stEle_PAS"       , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "Pt1stJet_PAS"             , 100 , 0       , 1000	 ); 
   CreateUserTH1D( "Pt2ndJet_PAS"             , 100 , 0       , 1000	 ); 
   CreateUserTH1D( "Eta1stJet_PAS"            , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Eta2ndJet_PAS"            , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stJet_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Phi2ndJet_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "TCHE1stJet_PAS"           , 100 , 0       , 20	 ); 
   CreateUserTH1D( "TCHE2ndJet_PAS"           , 100 , 0       , 20	 ); 
   CreateUserTH1D( "nMuon_PtCut_IDISO_PAS"    , 16  , -0.5    , 15.5	 ); 
   CreateUserTH1D( "MTenu_PAS"                , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "MT_charged_enu_PAS"       , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "MT_type1_enu_PAS"         , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "Ptenu_PAS"		      , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTlep_PAS"                , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTjet_PAS"                , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PAS"                   , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mjj_PAS"		      , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mej1_PAS"                 , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mej2_PAS"                 , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "DCotTheta1stEle_PAS"      , 100 , 0.0     , 1.0      );
   CreateUserTH1D( "Dist1stEle_PAS"           , 100 , 0.0     , 1.0      );  
   CreateUserTH1D( "DR_Ele1Jet1_PAS"	      , 100 , 0       , 10       ); 
   CreateUserTH1D( "DR_Ele1Jet2_PAS"	      , 100 , 0       , 10       ); 
   CreateUserTH1D( "DR_Jet1Jet2_PAS"	      , 100 , 0       , 10       ); 
   CreateUserTH1D( "minDR_EleJet_PAS"         , 100 , 0       , 10       ); 
   CreateUserTH1D( "mDPhi1stEleMET_PAS"       , 100 , 0.      ,  3.14159 );
   CreateUserTH1D( "mDPhi1stJetMET_PAS"       , 100 , 0.      ,  3.14159 );
   CreateUserTH1D( "mDPhi2ndJetMET_PAS"       , 100 , 0.      ,  3.14159 );

   CreateUserTH1D( "MT_GoodVtxLTE3_PAS"       , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MT_GoodVtxGTE4_LTE8_PAS"  , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MT_GoodVtxGTE9_LTE15_PAS" , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MT_GoodVtxGTE16_PAS"      , 200 , 0.      ,  1000    );

   CreateUserTH1D( "GeneratorWeight"       , 200 , -2.0    , 2.0      );
   CreateUserTH1D( "PileupWeight"          , 200 , -2.0    , 2.0      );

   CreateUserTH1D( "nVertex_PAS"           ,    31   , -0.5   , 30.5	 ) ; 
   CreateUserTH1D( "nVertex_good_PAS"      ,    31   , -0.5   , 30.5	 ) ; 

   CreateUserTH1D( "MTCharged_GoodVtxLTE3_PAS"       , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTCharged_GoodVtxGTE4_LTE8_PAS"  , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTCharged_GoodVtxGTE9_LTE15_PAS" , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTCharged_GoodVtxGTE16_PAS"      , 200 , 0.      ,  1000    );

   CreateUserTH1D( "MTType1_GoodVtxLTE3_PAS"       , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTType1_GoodVtxGTE4_LTE8_PAS"  , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTType1_GoodVtxGTE9_LTE15_PAS" , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTType1_GoodVtxGTE16_PAS"      , 200 , 0.      ,  1000    );
   
   CreateUserTH1D( "MTenu_50_110", 200, 40, 140 );
   
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
     // Do pileup re-weighting
     //--------------------------------------------------------------------------
     
     // int NPILEUP_AVE = int( (nPileUpInt_BXminus1 + nPileUpInt_BX0 + nPileUpInt_BXplus1)/3 );
     int NPILEUP_AVE = int( nPileUpInt_BX0 );
     int NPILEUP_FINAL = min ( NPILEUP_AVE , 25 );
     double pileup_weight = getPileupWeight ( NPILEUP_FINAL, isData ) ;

     //--------------------------------------------------------------------------
     // Check good run list
     //--------------------------------------------------------------------------
     
     int passedJSON = passJSON ( run, ls , isData ) ;

     //--------------------------------------------------------------------------
     // 7 trigger paths (some with multiple versions)
     //--------------------------------------------------------------------------
     
     // Number of times a path fired per event ( should be 0 or 1 )
     
     int N_Photon30_CIdVL  = 0;
     int N_Photon50_CIdVL  = 0;
     int N_Photon75_CIdVL  = 0;
     int N_Photon90_CIdVL  = 0;
     int N_Photon125       = 0;
     int N_Photon135       = 0;
     int N_Photon400       = 0;
     
     // Trigger prescale in an event
     
     int PS_Photon30_CIdVL = 0;
     int PS_Photon50_CIdVL = 0;
     int PS_Photon75_CIdVL = 0;
     int PS_Photon90_CIdVL = 0;
     int PS_Photon125      = 0;
     int PS_Photon135      = 0;
     int PS_Photon400      = 0;
     
     //--------------------------------------------------------------------------
     // Find the right prescale for this event
     //--------------------------------------------------------------------------
     
     // Did the HLT_Photon30_CaloIdVL trigger fire?
     
     if ( H_Photon30_CIdVL_1 > 0 && H_Photon30_CIdVL_1 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_1; } 
     if ( H_Photon30_CIdVL_2 > 0 && H_Photon30_CIdVL_2 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_2; } 
     if ( H_Photon30_CIdVL_3 > 0 && H_Photon30_CIdVL_3 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_3; } 
     if ( H_Photon30_CIdVL_4 > 0 && H_Photon30_CIdVL_4 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_4; } 
     if ( H_Photon30_CIdVL_5 > 0 && H_Photon30_CIdVL_5 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_5; } 
     if ( H_Photon30_CIdVL_6 > 0 && H_Photon30_CIdVL_6 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_6; } 
     if ( H_Photon30_CIdVL_7 > 0 && H_Photon30_CIdVL_7 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_7; } 
     
     // Did the HLT_Photon50_CaloIdVL trigger fire?
     
     if ( H_Photon50_CIdVL_1 > 0 && H_Photon50_CIdVL_1 != 999 ) { N_Photon50_CIdVL++; PS_Photon50_CIdVL = H_Photon50_CIdVL_1; } 
     if ( H_Photon50_CIdVL_2 > 0 && H_Photon50_CIdVL_2 != 999 ) { N_Photon50_CIdVL++; PS_Photon50_CIdVL = H_Photon50_CIdVL_2; } 
     if ( H_Photon50_CIdVL_3 > 0 && H_Photon50_CIdVL_3 != 999 ) { N_Photon50_CIdVL++; PS_Photon50_CIdVL = H_Photon50_CIdVL_3; } 
     if ( H_Photon50_CIdVL_4 > 0 && H_Photon50_CIdVL_4 != 999 ) { N_Photon50_CIdVL++; PS_Photon50_CIdVL = H_Photon50_CIdVL_4; } 
     
     // Did the HLT_Photon75_CaloIdVL trigger fire?
     
     if ( H_Photon75_CIdVL_1 > 0 && H_Photon75_CIdVL_1 != 999 ) { N_Photon75_CIdVL++; PS_Photon75_CIdVL = H_Photon75_CIdVL_1; } 
     if ( H_Photon75_CIdVL_2 > 0 && H_Photon75_CIdVL_2 != 999 ) { N_Photon75_CIdVL++; PS_Photon75_CIdVL = H_Photon75_CIdVL_2; } 
     if ( H_Photon75_CIdVL_3 > 0 && H_Photon75_CIdVL_3 != 999 ) { N_Photon75_CIdVL++; PS_Photon75_CIdVL = H_Photon75_CIdVL_3; } 
     if ( H_Photon75_CIdVL_4 > 0 && H_Photon75_CIdVL_4 != 999 ) { N_Photon75_CIdVL++; PS_Photon75_CIdVL = H_Photon75_CIdVL_4; } 
     if ( H_Photon75_CIdVL_5 > 0 && H_Photon75_CIdVL_5 != 999 ) { N_Photon75_CIdVL++; PS_Photon75_CIdVL = H_Photon75_CIdVL_5; } 
     if ( H_Photon75_CIdVL_6 > 0 && H_Photon75_CIdVL_6 != 999 ) { N_Photon75_CIdVL++; PS_Photon75_CIdVL = H_Photon75_CIdVL_6; } 
     if ( H_Photon75_CIdVL_7 > 0 && H_Photon75_CIdVL_7 != 999 ) { N_Photon75_CIdVL++; PS_Photon75_CIdVL = H_Photon75_CIdVL_7; }
     
     // Did the HLT_Photon90_CaloIdVL trigger fire?
     
     if ( H_Photon90_CIdVL_1 > 0 && H_Photon90_CIdVL_1 != 999 ) { N_Photon90_CIdVL++; PS_Photon90_CIdVL = H_Photon90_CIdVL_1; } 
     if ( H_Photon90_CIdVL_2 > 0 && H_Photon90_CIdVL_2 != 999 ) { N_Photon90_CIdVL++; PS_Photon90_CIdVL = H_Photon90_CIdVL_2; } 
     if ( H_Photon90_CIdVL_3 > 0 && H_Photon90_CIdVL_3 != 999 ) { N_Photon90_CIdVL++; PS_Photon90_CIdVL = H_Photon90_CIdVL_3; } 
     if ( H_Photon90_CIdVL_4 > 0 && H_Photon90_CIdVL_4 != 999 ) { N_Photon90_CIdVL++; PS_Photon90_CIdVL = H_Photon90_CIdVL_4; } 
     
     // Did the HLT_Photon125 trigger fire?
     
     if ( H_Photon125_1      > 0 && H_Photon125_1      != 999 ) { N_Photon125     ++; PS_Photon125      = H_Photon125_1     ; } 
     if ( H_Photon125_2      > 0 && H_Photon125_2      != 999 ) { N_Photon125     ++; PS_Photon125      = H_Photon125_2     ; } 
     
     // Did the HLT_Photon135 trigger fire?
     
     if ( H_Photon135_1      > 0 && H_Photon135_1      != 999 ) { N_Photon135     ++; PS_Photon135      = H_Photon135_1     ; } 
     if ( H_Photon135_2      > 0 && H_Photon135_2      != 999 ) { N_Photon135     ++; PS_Photon135      = H_Photon135_2     ; } 
     
     // Did the HLT_Photon400 trigger fire?
     
     if ( H_Photon400_1      > 0 && H_Photon400_1      != 999 ) { N_Photon400     ++; PS_Photon400      = H_Photon400_1     ; } 
     
     // Sanity check: make sure two versions of the same trigger didn't fire in the same event (impossible)
     
     if ( N_Photon30_CIdVL > 1 ) { std::cout << "ERROR: trigger overlap in N_Photon30_CIdVL" << std::endl; exit (0); }
     if ( N_Photon50_CIdVL > 1 ) { std::cout << "ERROR: trigger overlap in N_Photon50_CIdVL" << std::endl; exit (0); }
     if ( N_Photon75_CIdVL > 1 ) { std::cout << "ERROR: trigger overlap in N_Photon75_CIdVL" << std::endl; exit (0); }
     if ( N_Photon90_CIdVL > 1 ) { std::cout << "ERROR: trigger overlap in N_Photon90_CIdVL" << std::endl; exit (0); }
     if ( N_Photon125      > 1 ) { std::cout << "ERROR: trigger overlap in N_Photon125"      << std::endl; exit (0); }
     if ( N_Photon135      > 1 ) { std::cout << "ERROR: trigger overlap in N_Photon135"      << std::endl; exit (0); }
     if ( N_Photon400      > 1 ) { std::cout << "ERROR: trigger overlap in N_Photon400"      << std::endl; exit (0); }
     
     // What is the lowest-prescale trigger that this electron could have fired?
     
     int min_prescale      = 999999;
     std::string min_prescale_name("");

     if ( N_Photon30_CIdVL != 0 && QCDFakeEle1_Pt > 30. * trigger_tolerance  && PS_Photon30_CIdVL <= min_prescale ) { min_prescale = PS_Photon30_CIdVL; min_prescale_name = std::string("PS_Photon30_CIdVL"); }
     if ( N_Photon50_CIdVL != 0 && QCDFakeEle1_Pt > 50. * trigger_tolerance  && PS_Photon50_CIdVL <= min_prescale ) { min_prescale = PS_Photon50_CIdVL; min_prescale_name = std::string("PS_Photon50_CIdVL"); }
     if ( N_Photon75_CIdVL != 0 && QCDFakeEle1_Pt > 75. * trigger_tolerance  && PS_Photon75_CIdVL <= min_prescale ) { min_prescale = PS_Photon75_CIdVL; min_prescale_name = std::string("PS_Photon75_CIdVL"); }
     if ( N_Photon90_CIdVL != 0 && QCDFakeEle1_Pt > 90. * trigger_tolerance  && PS_Photon90_CIdVL <= min_prescale ) { min_prescale = PS_Photon90_CIdVL; min_prescale_name = std::string("PS_Photon90_CIdVL"); }
     if ( N_Photon125      != 0 && QCDFakeEle1_Pt > 125.* trigger_tolerance  && PS_Photon125      <= min_prescale ) { min_prescale = PS_Photon125     ; min_prescale_name = std::string("PS_Photon125"     ); }
     if ( N_Photon135      != 0 && QCDFakeEle1_Pt > 135.* trigger_tolerance  && PS_Photon135      <= min_prescale ) { min_prescale = PS_Photon135     ; min_prescale_name = std::string("PS_Photon135"     ); }
     if ( N_Photon400      != 0 && QCDFakeEle1_Pt > 400.* trigger_tolerance  && PS_Photon400      <= min_prescale ) { min_prescale = PS_Photon400     ; min_prescale_name = std::string("PS_Photon400"     ); }
     
     // If we find a suitable trigger, scale this event by that trigger's prescale
     
     int passedHLT = 0;
     if ( min_prescale != 999999 ) { // if I found some suitable trigger that fired
       passedHLT = 1;     
     } else {                        // if I was not able to find a suitable trigger that fired
       min_prescale = 0;            
     }

     if ( !isData ) {
       min_prescale = 1;
       passedHLT    = 1;
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
     
     if( fabs( QCDFakeEle1_Eta  ) < eleEta_bar )        ele1_isBarrel  = true;
     if( fabs( QCDFakeEle1_Eta  ) > eleEta_end1_min &&
	 fabs( QCDFakeEle1_Eta  ) < eleEta_end1_max )   ele1_isEndcap1 = true;
     if( fabs( QCDFakeEle1_Eta  ) > eleEta_end2_min &&
	 fabs( QCDFakeEle1_Eta  ) < eleEta_end2_max )   ele1_isEndcap2 = true;

     if( fabs( QCDFakeEle2_Eta  ) < eleEta_bar )        ele2_isBarrel  = true;
     if( fabs( QCDFakeEle2_Eta  ) > eleEta_end1_min &&
	 fabs( QCDFakeEle2_Eta  ) < eleEta_end1_max )   ele2_isEndcap1 = true;
     if( fabs( QCDFakeEle2_Eta  ) > eleEta_end2_min &&
	 fabs( QCDFakeEle2_Eta  ) < eleEta_end2_max )   ele2_isEndcap2 = true;

     
     //--------------------------------------------------------------------------
     // Determine which fake rates to use
     //--------------------------------------------------------------------------
     
     double fakeRate1  = 0.0;
     double eFakeRate1 = 0.0 ;
     
     if ( QCDFakeEle1_Pt < 100 ){
       if ( ele1_isBarrel  ) {
	 fakeRate1  = fakeRate_low_bar_p0  + fakeRate_low_bar_p1  * QCDFakeEle1_Pt;
	 eFakeRate1 = sqrt ((( eFakeRate_low_bar_p1  * QCDFakeEle1_Pt       ) * 
			     ( eFakeRate_low_bar_p1  * QCDFakeEle1_Pt       )) + 
			    (( eFakeRate_low_bar_p0  * eFakeRate_low_bar_p0 ) * 
			     ( eFakeRate_low_bar_p0  * eFakeRate_low_bar_p0 )));
       }
       if ( ele1_isEndcap1 ) {
	 fakeRate1  = fakeRate_low_end1_p0 + fakeRate_low_end1_p1 * QCDFakeEle1_Pt;
	 eFakeRate1 = sqrt ((( eFakeRate_low_end1_p1  * QCDFakeEle1_Pt       ) * 
			     ( eFakeRate_low_end1_p1  * QCDFakeEle1_Pt       )) + 
			    (( eFakeRate_low_end1_p0  * eFakeRate_low_end1_p0 ) * 
			     ( eFakeRate_low_end1_p0  * eFakeRate_low_end1_p0 )));
	 
       }
       if ( ele1_isEndcap2 ) {
	 fakeRate1 = fakeRate_low_end2_p0 + fakeRate_low_end2_p1 * QCDFakeEle1_Pt;
	 eFakeRate1 = sqrt ((( eFakeRate_low_end2_p1  * QCDFakeEle1_Pt       ) * 
			     ( eFakeRate_low_end2_p1  * QCDFakeEle1_Pt       )) + 
			    (( eFakeRate_low_end2_p0  * eFakeRate_low_end2_p0 ) * 
			     ( eFakeRate_low_end2_p0  * eFakeRate_low_end2_p0 )));
       } 
     }
     else if  ( QCDFakeEle1_Pt >= 100 ){
       if ( ele1_isBarrel  ) {
	 fakeRate1  = fakeRate_high_bar_p0  ;
	 eFakeRate1 = eFakeRate_high_bar_p0 ;
       }
       if ( ele1_isEndcap1 ) { 
	 fakeRate1 = fakeRate_high_end1_p0 ;
	 eFakeRate1 = eFakeRate_high_end1_p0 ;
       }
       if ( ele1_isEndcap2 ) {
	 fakeRate1 = fakeRate_high_end2_p0 ;
	 eFakeRate1 = eFakeRate_high_end2_p0 ;
       }
     }
     
     //--------------------------------------------------------------------------
     // Finally have the effective fake rate
     //--------------------------------------------------------------------------
     
     double fakeRate  = fakeRate1;
     double eFakeRate = eFakeRate1;


     //--------------------------------------------------------------------------
     // Calculate some variables:
     //--------------------------------------------------------------------------
     
     TLorentzVector loose_ele1, loose_ele2, jet1, jet2, met;
     loose_ele1.SetPtEtaPhiM ( QCDFakeEle1_Pt , QCDFakeEle1_Eta , QCDFakeEle1_Phi , 0.0 );
     loose_ele2.SetPtEtaPhiM ( QCDFakeEle2_Pt , QCDFakeEle2_Eta , QCDFakeEle2_Phi , 0.0 );
     jet1.SetPtEtaPhiM       ( JetLooseEle1_Pt, JetLooseEle1_Eta, JetLooseEle1_Phi, 0.0 );
     jet2.SetPtEtaPhiM       ( JetLooseEle2_Pt, JetLooseEle2_Eta, JetLooseEle2_Phi, 0.0 );
     met.SetPtEtaPhiM        ( MET_Pt         , 0.0             , MET_Phi         , 0.0 );

     TLorentzVector e1met = loose_ele1 + met;
     TLorentzVector j1j2 = jet1 + jet2;
     TLorentzVector e1j1 = loose_ele1 + jet1;
     TLorentzVector e1j2 = loose_ele1 + jet2;

     M_e1j1 = e1j1.M();
     M_e1j2 = e1j2.M();

     Pt_Ele1MET = e1met.Pt();
     MT_Ele1MET = sqrt(2 * QCDFakeEle1_Pt * MET_Pt  * (1 - cos(loose_ele1.DeltaPhi ( met) ) ) );

     mDPhi_METEle1= fabs(loose_ele1.DeltaPhi ( met ));
     mDPhi_METJet1= fabs(jet1.DeltaPhi ( met ));
     mDPhi_METJet2= fabs(jet2.DeltaPhi ( met ));

     sT_enujj = QCDFakeEle1_Pt + MET_Pt + JetLooseEle1_Pt + JetLooseEle2_Pt ;
     
     int nEle_QCDFake_Ana = 0;
     if      ( QCDFakeEle3_Pt > 35. ) nEle_QCDFake_Ana = 3;
     else if ( QCDFakeEle2_Pt > 35. ) nEle_QCDFake_Ana = 2;
     else if ( QCDFakeEle1_Pt > 35. ) nEle_QCDFake_Ana = 1;
     else                             nEle_QCDFake_Ana = 0;
     
     DR_Ele1Jet1 = loose_ele1.DeltaR ( jet1 ) ;
     DR_Ele1Jet2 = loose_ele1.DeltaR ( jet2 ) ;

     // bug: don't have QCD btag values
     double JetLooseEle1_btagTCHE = 0.0;
     double JetLooseEle2_btagTCHE = 0.0;

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "Reweighting"              , 1                       , min_prescale * fakeRate ); 
     fillVariableWithValue(   "PassJSON"                 , passedJSON              , min_prescale * fakeRate ); 
     									          
     // HLT variable							           
     fillVariableWithValue(   "PassHLT"                  , passedHLT               , min_prescale * fakeRate );
     
     // Dataset variable 
     fillVariableWithValue(   "PassDataset"                   , 1                  , min_prescale * fakeRate );

     // Filters
     fillVariableWithValue(   "PassHBHENoiseFilter"      , PassHBHENoiseFilter     , min_prescale * fakeRate );
     fillVariableWithValue(   "PassBeamHaloFilterTight"  , PassBeamHaloFilterTight , min_prescale * fakeRate );
									      
     // Muon variables ( for veto ) 					      
     fillVariableWithValue(   "nMuon"                    , nMuon_Ana               , min_prescale * fakeRate );
			                                      		                
     // 1st Electron variables				      		                
     fillVariableWithValue(   "nEle"                     , nEle_QCDFake_Ana        , min_prescale * fakeRate ); 
     fillVariableWithValue(   "Ele1_Pt"                  , QCDFakeEle1_Pt          , min_prescale * fakeRate );
     fillVariableWithValue(   "Ele1_Eta"                 , QCDFakeEle1_Eta         , min_prescale * fakeRate );
									           
     // MET variables	                                      		           
     fillVariableWithValue(   "MET"                      , MET_Pt                  , min_prescale * fakeRate );
     fillVariableWithValue(   "mDeltaPhiMETEle"          , mDPhi_METEle1           , min_prescale * fakeRate );
     									           
     // 1st JET variables                                     		           
     fillVariableWithValue(   "nJet"                     , nJetLooseEle_Ana        , min_prescale * fakeRate );
									           
     // 1st JET variables                                     		           
     if ( nJet_Stored > 0 ) { 						           
       fillVariableWithValue( "Jet1_Pt"                  , JetLooseEle1_Pt         , min_prescale * fakeRate );
       fillVariableWithValue( "Jet1_Eta"                 , JetLooseEle1_Eta        , min_prescale * fakeRate );
       fillVariableWithValue( "mDeltaPhiMET1stJet"       , mDPhi_METJet1           , min_prescale * fakeRate );
     }									           
     									           
     // 2nd JET variables                                     		           
     if ( nJet_Stored > 1 ) { 	                                      	           
       fillVariableWithValue( "Jet2_Pt"                  , JetLooseEle2_Pt         , min_prescale * fakeRate );
       fillVariableWithValue( "Jet2_Eta"                 , JetLooseEle2_Eta        , min_prescale * fakeRate );
       fillVariableWithValue( "ST"                       , sT_enujj                , min_prescale * fakeRate );
     }

     // 1 electron, 1 jet variables 
     if ( nEle_QCDFake > 0 && nJet_Ana > 0 ) { 
       fillVariableWithValue ( "DR_Ele1Jet1"             , DR_Ele1Jet1             , min_prescale * fakeRate );
     }

     // 1 electron, 2 jet variables 
     if ( nEle_QCDFake > 0 && nJet_Ana > 1 ) { 
       fillVariableWithValue ( "DR_Ele1Jet2"             , DR_Ele1Jet2             , min_prescale * fakeRate );
     }
     
     // Dummy variables
     fillVariableWithValue ("preselection",1, min_prescale * fakeRate );

     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------
     
     bool passed_preselection = passedAllPreviousCuts("preselection");
       
     if ( passed_preselection ) { 

       bool use_charged_met = (PFMETCharged < MET_Pt);

       if ( use_charged_met ) v_METCharged.SetMagPhi(PFMETCharged , PFMETChargedPhi );
       else                   v_METCharged.SetMagPhi(MET_Pt       , MET_Phi         );
       
       v_METType1.SetMagPhi  (PFMETType1Cor , PFMETPhiType1Cor);
       v_ele.SetMagPhi       (QCDFakeEle1_Pt, QCDFakeEle1_Phi );
       
       double deltaphi_charged = v_METCharged.DeltaPhi(v_ele);
       double deltaphi_type1   = v_METType1  .DeltaPhi(v_ele);
       double MTCharged = sqrt(2 * QCDFakeEle1_Pt * PFMETCharged  * (1 - cos(deltaphi_charged)) );
       double MTType1   = sqrt(2 * QCDFakeEle1_Pt * PFMETType1Cor * (1 - cos(deltaphi_type1  )) );

       double min_DR_EleJet = 999.0;
       double DR_Ele1Jet3 = 999.0;
       if ( nJetLooseEle_Ana > 2 ) {
	 TLorentzVector ele1, jet3;
	 ele1.SetPtEtaPhiM ( QCDFakeEle1_Pt, QCDFakeEle1_Eta, QCDFakeEle1_Phi, 0.0 );
	 jet3.SetPtEtaPhiM ( JetLooseEle3_Pt, JetLooseEle3_Eta, JetLooseEle3_Phi, 0.0 );
	 DR_Ele1Jet3 = ele1.DeltaR ( jet3 ) ;
       }

       if ( DR_Ele1Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet1;
       if ( DR_Ele1Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet2;
       if ( nJet_Ana > 2 ) {
	 if ( DR_Ele1Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet3;
       }

       FillUserTH1D( "MT_charged_enu_PAS"         , MTCharged                             , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "MT_type1_enu_PAS"           , MTType1                               , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "nElectron_PAS"              , nEle_QCDFake_Ana                      , pileup_weight * min_prescale * fakeRate); 
       FillUserTH1D( "nMuon_PAS"                  , nMuon_Ana                             , pileup_weight * min_prescale * fakeRate); 
       FillUserTH1D( "Pt1stEle_PAS"	          , QCDFakeEle1_Pt                        , pileup_weight * min_prescale * fakeRate); 
       FillUserTH1D( "Eta1stEle_PAS"	          , QCDFakeEle1_Eta                       , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "Phi1stEle_PAS"	          , QCDFakeEle1_Phi                       , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "Charge1stEle_PAS"           , QCDFakeEle1_Charge                    , pileup_weight * min_prescale * fakeRate);   
       FillUserTH1D( "METSig_PAS"	          , PFMETSig                              , pileup_weight * min_prescale * fakeRate);   
       FillUserTH1D( "MET_PAS"                    , MET_Pt                                , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "METPhi_PAS"	          , MET_Phi                               , pileup_weight * min_prescale * fakeRate);   
       FillUserTH1D( "METCharged_PAS"             , PFMETCharged                          , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "METChargedPhi_PAS"          , PFMETChargedPhi                       , pileup_weight * min_prescale * fakeRate);   
       FillUserTH1D( "METType1_PAS"               , PFMETType1Cor                         , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "METType1Phi_PAS"            , PFMETPhiType1Cor                      , pileup_weight * min_prescale * fakeRate);   
       FillUserTH1D( "minMETPt1stEle_PAS"         , TMath::Min ( QCDFakeEle1_Pt, MET_Pt  ), pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "Pt1stJet_PAS"               , JetLooseEle1_Pt                       , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "Pt2ndJet_PAS"               , JetLooseEle2_Pt                       , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "Eta1stJet_PAS"              , JetLooseEle1_Eta                      , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "Eta2ndJet_PAS"              , JetLooseEle2_Eta                      , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "Phi1stJet_PAS"              , JetLooseEle1_Phi                      , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "Phi2ndJet_PAS"	          , JetLooseEle2_Phi                      , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "TCHE1stJet_PAS"             , JetLooseEle1_btagTCHE                 , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "TCHE2ndJet_PAS"             , JetLooseEle2_btagTCHE                 , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "nMuon_PtCut_IDISO_PAS"      , nMuon_Ana                             , pileup_weight * min_prescale * fakeRate); 
       FillUserTH1D( "MTenu_PAS"                  , MT_Ele1MET                            , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "Ptenu_PAS"	          , Pt_Ele1MET                            , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "sTlep_PAS"                  , QCDFakeEle1_Pt + MET_Pt               , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "sTjet_PAS"                  , JetLooseEle1_Pt + JetLooseEle2_Pt     , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "sT_PAS"                     , sT_enujj                              , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "Mjj_PAS"	                  , M_j1j2                                , pileup_weight * min_prescale * fakeRate);   
       FillUserTH1D( "DCotTheta1stEle_PAS"        , QCDFakeEle1_DCotTheta                 , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "Dist1stEle_PAS"             , QCDFakeEle1_Dist                      , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "mDPhi1stEleMET_PAS"         , mDPhi_METEle1                         , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "mDPhi1stJetMET_PAS"         , mDPhi_METJet1                         , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "mDPhi2ndJetMET_PAS"         , mDPhi_METJet2                         , pileup_weight * min_prescale * fakeRate); 
       FillUserTH1D( "Mej1_PAS"                   , M_e1j1                                , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "Mej2_PAS"                   , M_e1j2                                , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "DR_Ele1Jet1_PAS"	          , DR_Ele1Jet1                           , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "DR_Ele1Jet2_PAS"	          , DR_Ele1Jet2                           , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "DR_Jet1Jet2_PAS"	          , DR_Jet1Jet2                           , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "minDR_EleJet_PAS"           , min_DR_EleJet                         , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "nVertex_PAS"                , nVertex                               , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "nVertex_good_PAS"           , nVertex_good                          , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "MatchPhotonConv1stEle_PAS"  , Ele1_MatchPhotConv                    , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "MatchPhotonConv2ndEle_PAS"  , Ele2_MatchPhotConv                    , pileup_weight * min_prescale * fakeRate);
       FillUserTH1D( "GeneratorWeight"            , -1.0             );
       FillUserTH1D( "PileupWeight"               , -1.0             );

       if ( MT_Ele1MET > 50 && MT_Ele1MET < 110 ){
	 FillUserTH1D( "MTenu_50_110", MT_Ele1MET, pileup_weight * min_prescale * fakeRate ) ;
       }

       if ( nVertex_good >= 0 && nVertex_good <= 3 ) {
	 FillUserTH1D( "MT_GoodVtxLTE3_PAS"              , MT_Ele1MET, pileup_weight * min_prescale * fakeRate ) ;
	 FillUserTH1D( "MTCharged_GoodVtxLTE3_PAS"       , MTCharged , pileup_weight * min_prescale * fakeRate ) ;
	 FillUserTH1D( "MTType1_GoodVtxLTE3_PAS"         , MTType1   , pileup_weight * min_prescale * fakeRate ) ;
       }						 
       							 
       if ( nVertex_good >= 4 && nVertex_good <= 8 ) {	 
	 FillUserTH1D( "MT_GoodVtxGTE4_LTE8_PAS"         , MT_Ele1MET, pileup_weight * min_prescale * fakeRate ) ;
	 FillUserTH1D( "MTCharged_GoodVtxGTE4_LTE8_PAS"  , MTCharged , pileup_weight * min_prescale * fakeRate ) ;
	 FillUserTH1D( "MTType1_GoodVtxGTE4_LTE8_PAS"    , MTType1   , pileup_weight * min_prescale * fakeRate ) ;
       }
       
       if ( nVertex_good >= 9 && nVertex_good <= 15) {
	 FillUserTH1D( "MT_GoodVtxGTE9_LTE15_PAS"        , MT_Ele1MET, pileup_weight * min_prescale * fakeRate ) ;
	 FillUserTH1D( "MTCharged_GoodVtxGTE9_LTE15_PAS" , MTCharged , pileup_weight * min_prescale * fakeRate ) ;
	 FillUserTH1D( "MTType1_GoodVtxGTE9_LTE15_PAS"   , MTType1   , pileup_weight * min_prescale * fakeRate ) ;
       }
       
       if ( nVertex_good >= 16                     ) {
	 FillUserTH1D( "MT_GoodVtxGTE16_PAS"             , MT_Ele1MET, pileup_weight * min_prescale * fakeRate ) ;
	 FillUserTH1D( "MTCharged_GoodVtxGTE16_PAS"      , MTCharged , pileup_weight * min_prescale * fakeRate ) ;
	 FillUserTH1D( "MTType1_GoodVtxGTE16_PAS"        , MTType1   , pileup_weight * min_prescale * fakeRate ) ;
       }
              
     }
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
