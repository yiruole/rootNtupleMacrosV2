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

   // override the fake rate?
   double fakeRate_override = getPreCutValue1("fakeRate_override");
   bool override_fakeRate = ( fakeRate_override > 0.0 );

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------
   
   CreateUserTH1D( "nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "nMuon_PAS"             ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "nJet_PAS"              ,    10  , -0.5    , 9.5      );
   CreateUserTH1D( "Pt1stEle_PAS"	   , 	100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta1stEle_PAS"	   , 	100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stEle_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Pt2ndEle_PAS"	   , 	100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta2ndEle_PAS"	   , 	100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi2ndEle_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Charge1stEle_PAS"	   , 	2   , -1.0001 , 1.0001	 ); 
   CreateUserTH1D( "Charge2ndEle_PAS"	   , 	2   , -1.0001 , 1.0001	 ); 
   CreateUserTH1D( "MET_PAS"               ,    200 , 0       , 1000	 ); 
   CreateUserTH1D( "METPhi_PAS"		   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "METCharged_PAS"        ,    200 , 0       , 1000	 ); 
   CreateUserTH1D( "METChargedPhi_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "METType1_PAS"          ,    200 , 0       , 1000	 ); 
   CreateUserTH1D( "METType1Phi_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Pt1stJet_PAS"          ,    100 , 0       , 1000	 ); 
   CreateUserTH1D( "Eta1stJet_PAS"         ,    100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "sT_PAS"                ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mee_PAS"		   ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me1j1_PAS"		   ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me2j1_PAS"		   ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Meejj_PAS"             ,    200 , 0       , 2000     );
   CreateUserTH1D( "Ptee_PAS"              ,    200 , 0       , 2000     );
   		                           
		                           
   CreateUserTH1D( "nVertex_PAS"           ,    31   , -0.5   , 30.5	 ) ; 
   CreateUserTH1D( "nVertex_good_PAS"      ,    31   , -0.5   , 30.5	 ) ; 
		                           
   CreateUserTH1D( "DR_Ele1Jet1_PAS"	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet1_PAS"	   , 	getHistoNBins("DR_Ele2Jet1"), getHistoMin("DR_Ele2Jet1"), getHistoMax("DR_Ele2Jet1")     ) ; 

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
     
     int NPILEUP_AVE = int( nPileUpInt_BX0 );
     int NPILEUP_FINAL = min( NPILEUP_AVE , 25 );
     double pileup_weight = getPileupWeight ( NPILEUP_FINAL, isData ) ;

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "PassJSON"                      , passedJSON ) ; 

     // Noise filters
     fillVariableWithValue(   "PassHBHENoiseFilter"           , PassHBHENoiseFilter ) ; 
     fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight ) ; 

     //--------------------------------------------------------------------------
     // Trigger
     //--------------------------------------------------------------------------

     int min_prescale;
     int passTrigger;
     std::string min_prescale_name;

     if ( isData ) {

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
 
       min_prescale      = 999999;
       min_prescale_name = std::string("");

       if ( N_Photon30_CIdVL != 0 && QCDFakeEle1_Pt > 30. * trigger_tolerance  && PS_Photon30_CIdVL <= min_prescale ) { min_prescale = PS_Photon30_CIdVL; min_prescale_name = std::string("PS_Photon30_CIdVL"); }
       if ( N_Photon50_CIdVL != 0 && QCDFakeEle1_Pt > 50. * trigger_tolerance  && PS_Photon50_CIdVL <= min_prescale ) { min_prescale = PS_Photon50_CIdVL; min_prescale_name = std::string("PS_Photon50_CIdVL"); }
       if ( N_Photon75_CIdVL != 0 && QCDFakeEle1_Pt > 75. * trigger_tolerance  && PS_Photon75_CIdVL <= min_prescale ) { min_prescale = PS_Photon75_CIdVL; min_prescale_name = std::string("PS_Photon75_CIdVL"); }
       if ( N_Photon90_CIdVL != 0 && QCDFakeEle1_Pt > 90. * trigger_tolerance  && PS_Photon90_CIdVL <= min_prescale ) { min_prescale = PS_Photon90_CIdVL; min_prescale_name = std::string("PS_Photon90_CIdVL"); }
       if ( N_Photon125      != 0 && QCDFakeEle1_Pt > 125.* trigger_tolerance  && PS_Photon125      <= min_prescale ) { min_prescale = PS_Photon125     ; min_prescale_name = std::string("PS_Photon125"     ); }
       if ( N_Photon135      != 0 && QCDFakeEle1_Pt > 135.* trigger_tolerance  && PS_Photon135      <= min_prescale ) { min_prescale = PS_Photon135     ; min_prescale_name = std::string("PS_Photon135"     ); }
       if ( N_Photon400      != 0 && QCDFakeEle1_Pt > 400.* trigger_tolerance  && PS_Photon400      <= min_prescale ) { min_prescale = PS_Photon400     ; min_prescale_name = std::string("PS_Photon400"     ); }

       // If we find a suitable trigger, scale this event by that trigger's prescale

       passTrigger = 0;
       if ( min_prescale != 999999 ) { // if I found some suitable trigger that fired
	 passTrigger = 1;     
       } else {                        // if I was not able to find a suitable trigger that fired
	 min_prescale = 0;            
       }

     }  // end if (isData) 
					       
     else { // i.e., if this is Monte Carlo
       min_prescale = 1;
       passTrigger = 1 ;
     }


     //--------------------------------------------------------------------------
     // What kind of event is this?
     //  - EB-EB
     //  - EB-EE
     //  - EE-EE
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

     bool isEBEB   = ( ele1_isBarrel  && ele2_isBarrel  );
     bool isEBEE1  = ( ele1_isBarrel  && ele2_isEndcap1 );
     bool isEBEE2  = ( ele1_isBarrel  && ele2_isEndcap2 );

     bool isEE1EB  = ( ele1_isEndcap1 && ele2_isBarrel  );
     bool isEE1EE1 = ( ele1_isEndcap1 && ele2_isEndcap1 );
     bool isEE1EE2 = ( ele1_isEndcap1 && ele2_isEndcap2 );

     bool isEE2EB  = ( ele1_isEndcap2 && ele2_isBarrel  );
     bool isEE2EE1 = ( ele1_isEndcap2 && ele2_isEndcap1 );
     bool isEE2EE2 = ( ele1_isEndcap2 && ele2_isEndcap2 );

     bool isEBEE   = ( isEBEE1  || isEBEE2  || isEE1EB  || isEE2EB  );
     bool isEEEE   = ( isEE1EE1 || isEE1EE2 || isEE2EE1 || isEE2EE2 );

     //--------------------------------------------------------------------------
     // Determine which fake rates to use
     //--------------------------------------------------------------------------

     double fakeRate1  = 0.0;
     double fakeRate2  = 0.0;
     double eFakeRate1 = 0.0 ;
     double eFakeRate2 = 0.0;
     
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

     if ( QCDFakeEle2_Pt < 100 ){
       if ( ele1_isBarrel  ) {
	 fakeRate2  = fakeRate_low_bar_p0  + fakeRate_low_bar_p1  * QCDFakeEle2_Pt;
	 eFakeRate2 = sqrt ((( eFakeRate_low_bar_p1  * QCDFakeEle2_Pt       ) * 
			     ( eFakeRate_low_bar_p1  * QCDFakeEle2_Pt       )) + 
			    (( eFakeRate_low_bar_p0  * eFakeRate_low_bar_p0 ) * 
			     ( eFakeRate_low_bar_p0  * eFakeRate_low_bar_p0 )));
       }
       if ( ele1_isEndcap1 ) {
	 fakeRate2 = fakeRate_low_end1_p0 + fakeRate_low_end1_p1 * QCDFakeEle2_Pt;
	 eFakeRate2 = sqrt ((( eFakeRate_low_end1_p1  * QCDFakeEle2_Pt       ) * 
			     ( eFakeRate_low_end1_p1  * QCDFakeEle2_Pt       )) + 
			    (( eFakeRate_low_end1_p0  * eFakeRate_low_end1_p0 ) * 
			     ( eFakeRate_low_end1_p0  * eFakeRate_low_end1_p0 )));
       }
       if ( ele1_isEndcap2 ) {
	 fakeRate2 = fakeRate_low_end2_p0 + fakeRate_low_end2_p1 * QCDFakeEle2_Pt;
	 eFakeRate2 = sqrt ((( eFakeRate_low_end2_p1  * QCDFakeEle2_Pt       ) * 
			     ( eFakeRate_low_end2_p1  * QCDFakeEle2_Pt       )) + 
			    (( eFakeRate_low_end2_p0  * eFakeRate_low_end2_p0 ) * 
			     ( eFakeRate_low_end2_p0  * eFakeRate_low_end2_p0 )));
       }
     }

     else if  ( QCDFakeEle2_Pt >= 100 ){
       if ( ele1_isBarrel  ) {
	 fakeRate2  = fakeRate_high_bar_p0  ;
	 eFakeRate2 = eFakeRate_high_bar_p0;
       }
       if ( ele1_isEndcap1 ) { 
	 fakeRate2 = fakeRate_high_end1_p0 ;
	 eFakeRate2 = eFakeRate_high_end1_p0;
       }
       if ( ele1_isEndcap2 ) {
	 fakeRate2 = fakeRate_high_end2_p0 ;
	 eFakeRate2 = eFakeRate_high_end2_p0;
       }
     }

     //--------------------------------------------------------------------------
     // Finally have the effective fake rate
     //--------------------------------------------------------------------------

     double fakeRateEffective  = fakeRate1 + fakeRate2;
     double eFakeRateEffective = sqrt ( ( eFakeRate1 * eFakeRate1 ) +
					( eFakeRate2 * eFakeRate2 ) );

     // fakeRateEffective += eFakeRateEffective;

     //--------------------------------------------------------------------------
     // User has the option to use a flat fake rate (e.g. 1.0 = no fake rate)
     //--------------------------------------------------------------------------
     
     if ( override_fakeRate ) fakeRateEffective = fakeRate_override;

     //--------------------------------------------------------------------------
     // Bug: we don't have the number for number of loose jets stored... 
     // Have to derive it
     //--------------------------------------------------------------------------

     int nJetLooseEle_Stored = 0;
     if      ( JetLooseEle3_Pt > 10.0 ) nJetLooseEle_Stored = 3;
     else if ( JetLooseEle2_Pt > 10.0 ) nJetLooseEle_Stored = 2;
     else if ( JetLooseEle1_Pt > 10.0 ) nJetLooseEle_Stored = 1;
     else                               nJetLooseEle_Stored = 0;

     //--------------------------------------------------------------------------
     // How many loose electrons have HEEP ID?
     //--------------------------------------------------------------------------

     int nPass = 0;
     if ( QCDFakeEle1_PassID == 1 ) nPass ++;
     if ( QCDFakeEle2_PassID == 1 ) nPass ++;

     //--------------------------------------------------------------------------
     // Calculate a few missing variables
     //--------------------------------------------------------------------------

     TLorentzVector loose_ele1, loose_ele2 , jet1, jet2;
     loose_ele1.SetPtEtaPhiM ( QCDFakeEle1_Pt , QCDFakeEle1_Eta , QCDFakeEle1_Phi , 0.0 );
     loose_ele2.SetPtEtaPhiM ( QCDFakeEle2_Pt , QCDFakeEle2_Eta , QCDFakeEle2_Phi , 0.0 );
     jet1.SetPtEtaPhiM       ( JetLooseEle1_Pt, JetLooseEle1_Eta, JetLooseEle1_Phi, 0.0 );
     jet2.SetPtEtaPhiM       ( JetLooseEle2_Pt, JetLooseEle2_Eta, JetLooseEle2_Phi, 0.0 );

     TLorentzVector loose_e1e2 = loose_ele1 + loose_ele2;
     TLorentzVector j1j2 = jet1 + jet2;

     TLorentzVector e1j1 = loose_ele1 + jet1;
     TLorentzVector e2j1 = loose_ele2 + jet1;

     //--------------------------------------------------------------------------
     // Now fill cut values
     // DON'T use the pileup weight ... it's included by default
     // DO    use the prescale.  It's already 1.0 for MC.
     // DO    use the effective fake rate.  It will only mean anything when you run
     //       over data, though, so be sure to override it to 1.0 
     //       if you run over Monte Carlo
     //--------------------------------------------------------------------------
     
     fillVariableWithValue ( "PassHLT", passTrigger, min_prescale * fakeRateEffective ) ;
     
     // Electrons
     fillVariableWithValue(   "nEleLoose"                     , nEle_QCDFake     , min_prescale  * fakeRateEffective );
     fillVariableWithValue(   "nEleTight"                     , nPass            , min_prescale  * fakeRateEffective );
     if ( nEle_QCDFake >= 1 ) { 
       fillVariableWithValue( "Ele1_Pt"                       , QCDFakeEle1_Pt   , min_prescale  * fakeRateEffective ) ;
       fillVariableWithValue( "Ele1_Eta"                      , QCDFakeEle1_Eta  , min_prescale  * fakeRateEffective ) ;
     }
     if ( nEle_QCDFake >= 2 ) { 
       fillVariableWithValue( "Ele2_Pt"                       , QCDFakeEle2_Pt   , min_prescale  * fakeRateEffective ) ;
       fillVariableWithValue( "Ele2_Eta"                      , QCDFakeEle2_Eta  , min_prescale  * fakeRateEffective ) ;
       fillVariableWithValue( "M_e1e2"                        , loose_e1e2.M()   , min_prescale  * fakeRateEffective );
       fillVariableWithValue( "Pt_e1e2"                       , loose_e1e2.Pt()  , min_prescale  * fakeRateEffective );
     }

     // Jets
     fillVariableWithValue(   "nJet"                          , nJetLooseEle_Stored , min_prescale  * fakeRateEffective );
     if ( nJetLooseEle_Stored >= 1 ) {
       fillVariableWithValue( "Jet1_Pt"                       , JetLooseEle1_Pt   , min_prescale  * fakeRateEffective ) ;
       fillVariableWithValue( "Jet1_Eta"                      , JetLooseEle1_Eta , min_prescale  * fakeRateEffective ) ;
     }

     // Muons
     fillVariableWithValue(   "nMuon"                         , nMuon_Ana           , min_prescale  * fakeRateEffective ); 
										      			      
     // DeltaR
     if ( nEle_QCDFake >= 2 && nJetLooseEle_Stored >= 1) {
       double st = QCDFakeEle1_Pt + QCDFakeEle2_Pt + JetLooseEle1_Pt ;
       fillVariableWithValue( "DR_Ele1Jet1"                   , loose_ele1.DeltaR ( jet1 ), min_prescale  * fakeRateEffective );
       fillVariableWithValue( "DR_Ele2Jet1"                   , loose_ele2.DeltaR ( jet1 ), min_prescale  * fakeRateEffective );
       fillVariableWithValue( "sT_eej"                        , st                        , min_prescale  * fakeRateEffective );
     }      

     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Fill preselection plots
     // DO    use the pileup weight.  It's equal to 1.0 for data.  
     // DO    use the min_prescale.  It's equal to 1.0 for Monte Carlo
     // DO    use the effective fake rate.  It will only mean anything when you run
     //       over data, though, so be sure to override it to 1.0 
     //       if you run over Monte Carlo
     //--------------------------------------------------------------------------     

     bool passed_preselection = ( passedAllPreviousCuts("sT_eej") && passedCut ("sT_eej") );
     
     if ( passed_preselection ) {

       double sT_eej = QCDFakeEle1_Pt + QCDFakeEle2_Pt + JetLooseEle1_Pt ;

       FillUserTH1D("nElectron_PAS"        , nEle_QCDFake              , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("nMuon_PAS"            , nMuon_Stored              , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("nJet_PAS"             , nJetLooseEle_Stored       , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Pt1stEle_PAS"	   , QCDFakeEle1_Pt            , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Eta1stEle_PAS"	   , QCDFakeEle1_Eta           , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Phi1stEle_PAS"	   , QCDFakeEle1_Phi           , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Pt2ndEle_PAS"	   , QCDFakeEle2_Pt            , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Eta2ndEle_PAS"	   , QCDFakeEle2_Eta           , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Phi2ndEle_PAS"	   , QCDFakeEle2_Phi           , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Charge1stEle_PAS"	   , QCDFakeEle1_Charge        , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Charge2ndEle_PAS"	   , QCDFakeEle2_Charge        , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("MET_PAS"              , MET_Pt                    , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("METPhi_PAS"	   , MET_Phi                   , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("METCharged_PAS"       , PFMETCharged              , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("METChargedPhi_PAS"    , PFMETChargedPhi           , pileup_weight * min_prescale * fakeRateEffective );   
       FillUserTH1D("METType1_PAS"         , PFMETType1Cor             , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("METType1Phi_PAS"      , PFMETPhiType1Cor          , pileup_weight * min_prescale * fakeRateEffective );   
       FillUserTH1D("Pt1stJet_PAS"         , JetLooseEle1_Pt           , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Eta1stJet_PAS"        , JetLooseEle1_Eta          , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Phi1stJet_PAS"	   , JetLooseEle1_Phi          , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("sT_PAS"               , sT_eej                    , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Mee_PAS"		   , loose_e1e2.M()            , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Ptee_PAS"             , loose_e1e2.Pt()           , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("nVertex_PAS"          , nVertex                   , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("nVertex_good_PAS"     , nVertex_good              , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("DR_Ele1Jet1_PAS"	   , loose_ele1.DeltaR ( jet1 ), pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("DR_Ele2Jet1_PAS"	   , loose_ele2.DeltaR ( jet1 ), pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Me1j1_PAS"            , e1j1.M()                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Me2j1_PAS"            , e2j1.M()                  , pileup_weight * min_prescale * fakeRateEffective );
     }
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

