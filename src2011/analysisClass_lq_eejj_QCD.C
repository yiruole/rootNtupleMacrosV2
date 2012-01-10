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
   CreateUserTH1D( "METSig_PAS"            ,    200 , 0       , 1000	 ); 
   CreateUserTH1D( "METPhi_PAS"		   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "METCharged_PAS"        ,    200 , 0       , 1000	 ); 
   CreateUserTH1D( "METChargedPhi_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "METType1_PAS"          ,    200 , 0       , 1000	 ); 
   CreateUserTH1D( "METType1Phi_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Pt1stJet_PAS"          ,    100 , 0       , 1000	 ); 
   CreateUserTH1D( "Pt2ndJet_PAS"          ,    100 , 0       , 1000	 ); 
   CreateUserTH1D( "Eta1stJet_PAS"         ,    100 , -5      , 5	 ); 
   CreateUserTH1D( "Eta2ndJet_PAS"         ,    100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Phi2ndJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "sTlep_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTjet_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PAS"                ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mjj_PAS"		   ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mee_PAS"		   ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "MTenu_PAS"             ,    200 , 0       , 1000	 ); 
   CreateUserTH1D( "Me1j1_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me1j2_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me2j1_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me2j2_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me1j_selected_PAS"     ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me2j_selected_PAS"     ,    200 , 0       , 2000     );
   CreateUserTH1D( "Mej_selected_avg_PAS"  ,    200 , 0       , 2000     );
   CreateUserTH1D( "Meejj_PAS"             ,    200 , 0       , 2000     );
		                           
   CreateUserTH1D( "Ptee_PAS"              ,    200 , 0       , 2000     );
   		                           
   CreateUserTH1D( "DCotTheta1stEle_PAS"   ,    100 , 0.0, 1.0);
   CreateUserTH1D( "Dist1stEle_PAS"        ,    100 , 0.0, 1.0);  
   CreateUserTH1D( "DCotTheta2ndEle_PAS"   ,    100 , 0.0, 1.0);
   CreateUserTH1D( "Dist2ndEle_PAS"        ,    100 , 0.0, 1.0);  
		                           
   CreateUserTH1D( "nVertex_PAS"           ,    31   , -0.5   , 30.5	 ) ; 
   CreateUserTH1D( "nVertex_good_PAS"      ,    31   , -0.5   , 30.5	 ) ; 
		                           
   CreateUserTH1D( "DR_Ele1Jet1_PAS"	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele1Jet2_PAS"	   , 	getHistoNBins("DR_Ele1Jet2"), getHistoMin("DR_Ele1Jet2"), getHistoMax("DR_Ele1Jet2")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet1_PAS"	   , 	getHistoNBins("DR_Ele2Jet1"), getHistoMin("DR_Ele2Jet1"), getHistoMax("DR_Ele2Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet2_PAS"	   , 	getHistoNBins("DR_Ele2Jet2"), getHistoMin("DR_Ele2Jet2"), getHistoMax("DR_Ele2Jet2")     ) ; 
   CreateUserTH1D( "DR_Jet1Jet2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
   CreateUserTH1D( "DR_Ele1Ele2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
   CreateUserTH1D( "minDR_EleJet_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 

   CreateUserTH2D( "Me1jVsMe2j_selected",     200, 0, 2000, 200, 0, 2000) ;
   CreateUserTH2D( "Me1jVsMe2j_rejected",     200, 0, 2000, 200, 0, 2000) ;

   CreateUserTH1D( "Mee_80_100_Preselection", 200, 60, 120 );
   CreateUserTH1D( "Mee_70_110_Preselection", 200, 60, 120 );
   
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
   CreateUserTH1D( "GeneratorWeight", 100, -2.0 , 2.0 );

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

     if ( isData && run > 175771 ) { // This trigger only available in 2011B
       if ( H_Photon30_CIdVL_8 > 0 && H_Photon30_CIdVL_8 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_8; }
     }
     
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

     if ( isData && run > 175771 ) { // This trigger only available in 2011B
       if ( H_Photon400_2      > 0 && H_Photon400_2      != 999 ) { N_Photon400     ++; PS_Photon400      = H_Photon400_2     ; }
     }
     
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
     
     int passHLT = 0;
     if ( min_prescale != 999999 ) { // if I found some suitable trigger that fired
       passHLT = 1;     
     } else {                        // if I was not able to find a suitable trigger that fired
       min_prescale = 0;            
     }

     if ( !isData ) {
       min_prescale = 1;
       passHLT    = 1;
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

     bool ele1_isEndcap = ( ele1_isEndcap1 || ele1_isEndcap2 ) ;
     bool ele2_isEndcap = ( ele2_isEndcap1 || ele2_isEndcap2 ) ;

     bool isEBEB = ( ele1_isBarrel && ele2_isBarrel ) ;
     bool isEBEE = ( ( ele1_isBarrel && ele2_isEndcap ) ||
		     ( ele2_isBarrel && ele1_isEndcap ) );
     bool isEEEE = ( ele1_isEndcap && ele2_isEndcap ) ;
     bool isEB   = ( isEBEB || isEBEE ) ;

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

     double fakeRateEffective  = fakeRate1 * fakeRate2;
     double eFakeRateEffective = fakeRateEffective * sqrt (  ( eFakeRate1 / fakeRate1 ) * ( eFakeRate1 / fakeRate1 ) +
							     ( eFakeRate2 / fakeRate2 ) * ( eFakeRate2 / fakeRate2 ) );

     //--------------------------------------------------------------------------
     // Calculate some variables:
     //--------------------------------------------------------------------------
     
     TLorentzVector loose_ele1, loose_ele2, jet1, jet2, met, muon;
     loose_ele1.SetPtEtaPhiM ( QCDFakeEle1_Pt , QCDFakeEle1_Eta , QCDFakeEle1_Phi , 0.0 );
     loose_ele2.SetPtEtaPhiM ( QCDFakeEle2_Pt , QCDFakeEle2_Eta , QCDFakeEle2_Phi , 0.0 );
     jet1.SetPtEtaPhiM       ( JetLooseEle1_Pt, JetLooseEle1_Eta, JetLooseEle1_Phi, 0.0 );
     jet2.SetPtEtaPhiM       ( JetLooseEle2_Pt, JetLooseEle2_Eta, JetLooseEle2_Phi, 0.0 );
     met.SetPtEtaPhiM        ( MET_Pt         , 0.0             , MET_Phi         , 0.0 );
     muon.SetPtEtaPhiM       ( Muon1_Pt       , Muon1_Eta       , Muon1_Phi       , 0.0 );

     TLorentzVector e1e2 = loose_ele1 + loose_ele2;
     TLorentzVector j1j2 = jet1 + jet2;
     TLorentzVector e1j1 = loose_ele1 + jet1;
     TLorentzVector e1j2 = loose_ele1 + jet2;
     TLorentzVector e2j1 = loose_ele2 + jet1;
     TLorentzVector e2j2 = loose_ele2 + jet2;
     TLorentzVector eejj = loose_ele1 + loose_ele2 + jet1 + jet2;

     double M_eejj = eejj.M();

     double DeltaR_Ele1Ele2 = loose_ele1.DeltaR( loose_ele2 ) ;

     M_e1j1 = e1j1.M();
     M_e1j2 = e1j2.M();
     M_e2j1 = e2j1.M();
     M_e2j2 = e2j2.M();
     M_e1e2 = e1e2.M();
     M_j1j2 = j1j2.M();

     Pt_e1e2 = e1e2.M();

     MT_Ele1MET = sqrt(2 * QCDFakeEle1_Pt * MET_Pt  * (1 - cos(loose_ele1.DeltaPhi ( met) ) ) );

     mDPhi_METEle1= fabs(loose_ele1.DeltaPhi ( met ));
     mDPhi_METJet1= fabs(jet1.DeltaPhi ( met ));
     mDPhi_METJet2= fabs(jet2.DeltaPhi ( met ));

     sT_eejj = QCDFakeEle1_Pt + QCDFakeEle2_Pt + JetLooseEle1_Pt + JetLooseEle2_Pt ;
     
     int nEle_QCDFake_Ana = 0;
     if      ( QCDFakeEle3_Pt > 35. ) nEle_QCDFake_Ana = 3;
     else if ( QCDFakeEle2_Pt > 35. ) nEle_QCDFake_Ana = 2;
     else if ( QCDFakeEle1_Pt > 35. ) nEle_QCDFake_Ana = 1;
     else                             nEle_QCDFake_Ana = 0;
     
     double DR_Ele1Ele2 = loose_ele1.DeltaR ( loose_ele2 ) ;
     DR_Ele1Jet1 = loose_ele1.DeltaR ( jet1 ) ;
     DR_Ele1Jet2 = loose_ele1.DeltaR ( jet2 ) ;
     DR_Ele2Jet1 = loose_ele2.DeltaR ( jet1 ) ;
     DR_Ele2Jet2 = loose_ele2.DeltaR ( jet2 ) ;
     DR_Jet1Jet2 = jet1.DeltaR( jet2 );

     // bug: don't have QCD btag values
     double JetLooseEle1_btagTCHE = 0.0;
     double JetLooseEle2_btagTCHE = 0.0;
     
     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // reweighting
     fillVariableWithValue ( "Reweighting", 1, fakeRateEffective * min_prescale ) ; 
     
     // JSON variable
     fillVariableWithValue(   "PassJSON"                      , passedJSON              , fakeRateEffective * min_prescale ) ; 

     // Dataset variable 
     fillVariableWithValue(   "PassDataset"                   , 1                       , fakeRateEffective * min_prescale ) ; 

     // Noise filters
     fillVariableWithValue(   "PassHBHENoiseFilter"           , PassHBHENoiseFilter     , fakeRateEffective * min_prescale ) ; 
     fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight , fakeRateEffective * min_prescale ) ; 
     // fillVariableWithValue(   "PassBPTX0"                     , PassBPTX0                  ,fakeRateEffective * min_prescale  ) ; 
     // fillVariableWithValue(   "PassPhysDecl"                  , PassPhysDecl               ,fakeRateEffective * min_prescale  ) ; 
     // fillVariableWithValue(   "PassBeamScraping"              , PassBeamScraping           ,fakeRateEffective * min_prescale  ) ; 
     // fillVariableWithValue(   "PassPrimaryVertex"             , PassPrimaryVertex          ,fakeRateEffective * min_prescale  ) ; 
     // fillVariableWithValue(   "PassBeamHaloFilterLoose"	      , PassBeamHaloFilterLoose	   ,fakeRateEffective * min_prescale  ) ; 
     // fillVariableWithValue(   "PassTrackingFailure"           , PassTrackingFailure        ,fakeRateEffective * min_prescale  ) ; 
     // fillVariableWithValue(   "PassCaloBoundaryDRFilter"      , PassCaloBoundaryDRFilter   ,fakeRateEffective * min_prescale  ) ; 
     // fillVariableWithValue(   "PassEcalMaskedCellDRFilter"    , PassEcalMaskedCellDRFilter ,fakeRateEffective * min_prescale  ) ; 



     // Electrons
     int PassNEle = 0;
     if ( nEle_QCDFake == 2 ) PassNEle = 1;
     
     // Muons
     int PassNMuon = 0;
     if ( nMuon_Ana == 0 ) PassNMuon = 1;

     fillVariableWithValue ( "PassHLT"                        , passHLT                 , fakeRateEffective * min_prescale ) ;
     
     // Electrons
     fillVariableWithValue(   "PassNEle"                      , PassNEle                , fakeRateEffective * min_prescale ) ;
     if ( nEle_QCDFake >= 1 ) { 
       fillVariableWithValue( "Ele1_Pt"                       , QCDFakeEle1_Pt          , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "Ele1_Eta"                      , QCDFakeEle1_Eta         , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "abs_Ele1_Eta"                  , fabs(QCDFakeEle1_Eta)   , fakeRateEffective * min_prescale ) ;
     }
     if ( nEle_QCDFake >= 2 ) { 
       fillVariableWithValue( "Ele2_Pt"                       , QCDFakeEle2_Pt          , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "Ele2_Eta"                      , QCDFakeEle2_Eta         , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "abs_Ele2_Eta"                  , fabs(QCDFakeEle2_Eta)   , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "M_e1e2"                        , M_e1e2                  , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "Pt_e1e2"                       , Pt_e1e2                 , fakeRateEffective * min_prescale ) ;
     }

     // Jets
     fillVariableWithValue(   "nJet"                          , nJetLooseEle_Ana        , fakeRateEffective * min_prescale ) ;
     if ( nJetLooseEle_Stored >= 1 ) { 						                
       fillVariableWithValue( "Jet1_Pt"                       , JetLooseEle1_Pt         , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "Jet1_Eta"                      , JetLooseEle1_Eta        , fakeRateEffective * min_prescale ) ;
     }
     if ( nJetLooseEle_Stored >= 2 ) { 
       fillVariableWithValue( "Jet2_Pt"                       , JetLooseEle2_Pt         , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "Jet2_Eta"                      , JetLooseEle2_Eta        , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "Pt_j1j2"                       , Pt_j1j2                 , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "M_j1j2"                        , M_j1j2                  , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "DR_Jet1Jet2"                   , DR_Jet1Jet2             , fakeRateEffective * min_prescale ) ;
     }

     // Muons
     fillVariableWithValue(   "PassNMuon"                     , PassNMuon               , fakeRateEffective * min_prescale ) ;

     // DeltaR
     if ( nEle_QCDFake >= 2 && nJetLooseEle_Stored >= 1) {
       fillVariableWithValue( "DR_Ele1Jet1"                   , DR_Ele1Jet1             , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "DR_Ele2Jet1"                   , DR_Ele2Jet1             , fakeRateEffective * min_prescale ) ;
       if(nJetLooseEle_Stored >= 2) {
	 fillVariableWithValue( "DR_Ele1Jet2"                 , DR_Ele1Jet2             , fakeRateEffective * min_prescale ) ;
	 fillVariableWithValue( "DR_Ele2Jet2"                 , DR_Ele2Jet2             , fakeRateEffective * min_prescale ) ;
       }
     }

     // sT
     if ( nEle_QCDFake >= 2 && nJetLooseEle_Stored >= 2) {
       fillVariableWithValue( "sT_eejj"                       , sT_eejj                 , fakeRateEffective  * min_prescale ) ;
     }      

     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------
     
     FillUserTH1D( "PileupWeight"   , -1.0 );
     FillUserTH1D( "GeneratorWeight", -1.0 );

     bool passed_preselection = ( passedAllPreviousCuts("M_e1e2") && passedCut ("M_e1e2") );
     
     if ( passed_preselection ) {


       
       double DR_Ele1Jet3 = 999.;
       double DR_Ele2Jet3 = 999.;
       double min_DR_EleJet = 999.;
       TLorentzVector eejj, e1e2mu, e1, j1, e2, j2,j3, mu, met;
       e1.SetPtEtaPhiM ( QCDFakeEle1_Pt, QCDFakeEle1_Eta, QCDFakeEle1_Phi, 0.0 );
       e2.SetPtEtaPhiM ( QCDFakeEle2_Pt, QCDFakeEle2_Eta, QCDFakeEle2_Phi, 0.0 );
       j1.SetPtEtaPhiM ( JetLooseEle1_Pt, JetLooseEle1_Eta, JetLooseEle1_Phi, 0.0 );
       j2.SetPtEtaPhiM ( JetLooseEle2_Pt, JetLooseEle2_Eta, JetLooseEle2_Phi, 0.0 );
       if ( nJetLooseEle_Ana > 2 ) {
	 j3.SetPtEtaPhiM ( JetLooseEle3_Pt, JetLooseEle3_Eta, JetLooseEle3_Phi, 0.0 );
	 DR_Ele1Jet3 = e1.DeltaR( j3 );
	 DR_Ele2Jet3 = e2.DeltaR( j3 );
       }
       
       
       if ( DR_Ele1Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet1;
       if ( DR_Ele1Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet2;
       if ( DR_Ele2Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet1;
       if ( DR_Ele2Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet2;
       if ( nJetLooseEle_Ana > 2 ) {
	 if ( DR_Ele1Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet3;
	 if ( DR_Ele2Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet3;
       }
       
       FillUserTH1D( "minDR_EleJet_PAS", min_DR_EleJet, pileup_weight * min_prescale * fakeRateEffective ) ;
       

       // TLorentzVector e1e2mu = loose_ele1 + loose_ele2 + muon ;
       // double MT_eemuMET = sqrt(2 * e1e2mu.Pt()    * MET_Pt  * (1 - cos(e1e2mu.DeltaPhi (met))));
       
       FillUserTH1D("nElectron_PAS"        , nEle_QCDFake                      , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("nMuon_PAS"            , nMuon_Stored                      , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("nJet_PAS"             , nJetLooseEle_Ana                  , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Pt1stEle_PAS"	   , QCDFakeEle1_Pt                    , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Eta1stEle_PAS"	   , QCDFakeEle1_Eta                   , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Phi1stEle_PAS"	   , QCDFakeEle1_Phi                   , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Pt2ndEle_PAS"	   , QCDFakeEle2_Pt                    , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Eta2ndEle_PAS"	   , QCDFakeEle2_Eta                   , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Phi2ndEle_PAS"	   , QCDFakeEle2_Phi                   , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Charge1stEle_PAS"	   , QCDFakeEle1_Charge                , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Charge2ndEle_PAS"	   , QCDFakeEle2_Charge                , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("MET_PAS"              , MET_Pt                            , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("METSig_PAS"           , PFMETSig                          , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("METPhi_PAS"	   , MET_Phi                           , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D( "METCharged_PAS"      , PFMETCharged                      , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D( "METChargedPhi_PAS"   , PFMETChargedPhi                   , pileup_weight * min_prescale * fakeRateEffective );   
       FillUserTH1D( "METType1_PAS"        , PFMETType1Cor                     , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D( "METType1Phi_PAS"     , PFMETPhiType1Cor                  , pileup_weight * min_prescale * fakeRateEffective );   
       FillUserTH1D("Pt1stJet_PAS"         , JetLooseEle1_Pt                   , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Pt2ndJet_PAS"         , JetLooseEle2_Pt                   , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Eta1stJet_PAS"        , JetLooseEle1_Eta                  , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Eta2ndJet_PAS"        , JetLooseEle2_Eta                  , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Phi1stJet_PAS"	   , JetLooseEle1_Phi                  , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Phi2ndJet_PAS"	   , JetLooseEle2_Phi                  , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("sTlep_PAS"            , QCDFakeEle1_Pt + QCDFakeEle2_Pt   , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("sTjet_PAS"            , JetLooseEle1_Pt + JetLooseEle2_Pt , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("sT_PAS"               , sT_eejj                           , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Mjj_PAS"		   , M_j1j2                            , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Mee_PAS"		   , M_e1e2                            , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D( "MTenu_PAS"           , MT_Ele1MET                        , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Me1j1_PAS"            , M_e1j1                            , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Me1j2_PAS"            , M_e1j2                            , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Me2j1_PAS"            , M_e2j1                            , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Me2j2_PAS"            , M_e2j2                            , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Ptee_PAS"             , Pt_e1e2                           , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DCotTheta1stEle_PAS"  , QCDFakeEle1_DCotTheta             , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Dist1stEle_PAS"       , QCDFakeEle1_Dist                  , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DCotTheta2ndEle_PAS"  , QCDFakeEle2_DCotTheta             , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Dist2ndEle_PAS"       , QCDFakeEle2_Dist                  , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("nVertex_PAS"          , nVertex                           , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("nVertex_good_PAS"     , nVertex_good                      , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DR_Ele1Jet1_PAS"	   , DR_Ele1Jet1                       , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DR_Ele1Jet2_PAS"	   , DR_Ele1Jet2                       , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DR_Ele2Jet1_PAS"	   , DR_Ele2Jet1                       , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DR_Ele2Jet2_PAS"	   , DR_Ele2Jet2                       , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DR_Jet1Jet2_PAS"	   , DR_Jet1Jet2                       , pileup_weight * min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DR_Ele1Ele2_PAS"	   , DR_Ele1Ele2                       , pileup_weight * min_prescale * fakeRateEffective ) ;

       if ( M_e1e2 > 80.0 && M_e1e2 < 100.0 ){
	 FillUserTH1D("Mee_80_100_Preselection", M_e1e2, pileup_weight * min_prescale * fakeRateEffective ) ;
	 if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_80_100_PAS", M_e1e2, pileup_weight * min_prescale * fakeRateEffective  ); 
	 else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_80_100_PAS", M_e1e2, pileup_weight * min_prescale * fakeRateEffective  ); 
	 else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_80_100_PAS", M_e1e2, pileup_weight * min_prescale * fakeRateEffective  ); 
	 if      ( isEB   ) FillUserTH1D( "Mee_EB_80_100_PAS"  , M_e1e2, pileup_weight * min_prescale * fakeRateEffective  ); 
       }

       if ( M_e1e2 > 70.0 && M_e1e2 < 110.0 ){
	 FillUserTH1D("Mee_70_110_Preselection", M_e1e2, pileup_weight * min_prescale * fakeRateEffective ) ;
	 if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_70_110_PAS", M_e1e2, pileup_weight * min_prescale * fakeRateEffective  ); 
	 else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_70_110_PAS", M_e1e2, pileup_weight * min_prescale * fakeRateEffective  ); 
	 else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_70_110_PAS", M_e1e2, pileup_weight * min_prescale * fakeRateEffective  ); 
	 if      ( isEB   ) FillUserTH1D( "Mee_EB_70_110_PAS"  , M_e1e2, pileup_weight * min_prescale * fakeRateEffective  ); 
       }
       
       if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_PAS", M_e1e2, pileup_weight * min_prescale * fakeRateEffective  ); 
       else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_PAS", M_e1e2, pileup_weight * min_prescale * fakeRateEffective  ); 
       else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_PAS", M_e1e2, pileup_weight * min_prescale * fakeRateEffective  ); 
       if      ( isEB   ) FillUserTH1D( "Mee_EB_PAS"  , M_e1e2, pileup_weight * min_prescale * fakeRateEffective  ); 
       
       FillUserTH1D("Meejj_PAS", M_eejj , pileup_weight * min_prescale * fakeRateEffective  );

       if ( fabs(M_e1j1-M_e2j2) < fabs(M_e1j2-M_e2j1) )  {

	 double M_ej_avg = (M_e1j1 + M_e2j2) / 2.0;

	 FillUserTH1D("Mej_selected_avg_PAS", M_ej_avg, pileup_weight * min_prescale * fakeRateEffective ) ;	   
	 FillUserTH1D("Me1j_selected_PAS"   , M_e1j1, pileup_weight * min_prescale * fakeRateEffective ) ;	   
	 FillUserTH1D("Me2j_selected_PAS"   , M_e2j2, pileup_weight * min_prescale * fakeRateEffective ) ;	   
	 FillUserTH2D( "Me1jVsMe2j_selected", M_e1j1, M_e2j2, pileup_weight * min_prescale * fakeRateEffective  ) ;
	 FillUserTH2D( "Me1jVsMe2j_rejected", M_e1j2, M_e2j1, pileup_weight * min_prescale * fakeRateEffective  ) ;
       }
       else {

	 double M_ej_avg = (M_e1j2 + M_e2j1) / 2.0;

	 FillUserTH1D("Mej_selected_avg_PAS", M_ej_avg, pileup_weight * min_prescale * fakeRateEffective ) ;	   
	 FillUserTH1D("Me1j_selected_PAS" , M_e1j2, pileup_weight * min_prescale * fakeRateEffective ) ;	   
	 FillUserTH1D("Me2j_selected_PAS" , M_e2j1, pileup_weight * min_prescale * fakeRateEffective ) ;	   
	 FillUserTH2D( "Me1jVsMe2j_selected", M_e1j2, M_e2j1, pileup_weight * min_prescale * fakeRateEffective  ) ;
	 FillUserTH2D( "Me1jVsMe2j_rejected", M_e1j1, M_e2j2, pileup_weight * min_prescale * fakeRateEffective  ) ;
       }
     }
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

