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
   
   // New movable break point

   double fakeRate_break_bar     = getPreCutValue4 ( "fakeRate_bar"  );
   double fakeRate_break_end1    = getPreCutValue4 ( "fakeRate_end1" );
   double fakeRate_break_end2    = getPreCutValue4 ( "fakeRate_end2" );
   
   double eFakeRate_break_bar    = getPreCutValue4 ( "eFakeRate_bar"  );
   double eFakeRate_break_end1   = getPreCutValue4 ( "eFakeRate_end1" );
   double eFakeRate_break_end2   = getPreCutValue4 ( "eFakeRate_end2" );
   
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
     
     double pileup_weight = getPileupWeight ( nPileUpInt_True, isData ) ;
     
     //--------------------------------------------------------------------------
     // Trigger
     //--------------------------------------------------------------------------

     short min_prescale = 0;
     int passTrigger = 0;

     if ( isData ) {

       if ( LooseEle1_hltPhotonPt > 0.0 ) { 
	 if ( H_Photon30_CIdVL > 0.1 && LooseEle1_hltPhotonPt >= 30.  && LooseEle1_hltPhotonPt < 50. ) { passTrigger = 1; min_prescale = H_Photon30_CIdVL; } 
	 if ( H_Photon50_CIdVL > 0.1 && LooseEle1_hltPhotonPt >= 50.  && LooseEle1_hltPhotonPt < 75. ) { passTrigger = 1; min_prescale = H_Photon50_CIdVL; } 
	 if ( H_Photon75_CIdVL > 0.1 && LooseEle1_hltPhotonPt >= 75.  && LooseEle1_hltPhotonPt < 90. ) { passTrigger = 1; min_prescale = H_Photon75_CIdVL; } 
	 if ( H_Photon90_CIdVL > 0.1 && LooseEle1_hltPhotonPt >= 90.  && LooseEle1_hltPhotonPt < 135.) { passTrigger = 1; min_prescale = H_Photon90_CIdVL; } 
	 if ( H_Photon135      > 0.1 && LooseEle1_hltPhotonPt >= 135. && LooseEle1_hltPhotonPt < 150.) { passTrigger = 1; min_prescale = H_Photon135     ; } 
	 if ( H_Photon150      > 0.1 && LooseEle1_hltPhotonPt >= 150.                                ) { passTrigger = 1; min_prescale = H_Photon150     ; } 
       }
       
     }  // end if (isData) 
					       
     else { 
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

     double break_point_ele1 = -1;
     if      ( ele1_isBarrel  ) break_point_ele1 = fakeRate_break_bar ;
     else if ( ele1_isEndcap1 ) break_point_ele1 = fakeRate_break_end1;
     else if ( ele1_isEndcap2 ) break_point_ele1 = fakeRate_break_end2;

     double break_point_ele2 = -1;
     if      ( ele2_isBarrel  ) break_point_ele2 = fakeRate_break_bar ;
     else if ( ele2_isEndcap1 ) break_point_ele2 = fakeRate_break_end1;
     else if ( ele2_isEndcap2 ) break_point_ele2 = fakeRate_break_end2;
     
     if ( LooseEle1_Pt < break_point_ele1 ){
       if ( ele1_isBarrel  ) {
	 fakeRate1  = fakeRate_low_bar_p0  + fakeRate_low_bar_p1  * LooseEle1_Pt;
	 eFakeRate1 = sqrt ((( eFakeRate_low_bar_p1  * LooseEle1_Pt       ) * 
			     ( eFakeRate_low_bar_p1  * LooseEle1_Pt       )) + 
			    (( eFakeRate_low_bar_p0  * eFakeRate_low_bar_p0 ) * 
			     ( eFakeRate_low_bar_p0  * eFakeRate_low_bar_p0 )));
       }
       if ( ele1_isEndcap1 ) {
	 fakeRate1  = fakeRate_low_end1_p0 + fakeRate_low_end1_p1 * LooseEle1_Pt;
	 eFakeRate1 = sqrt ((( eFakeRate_low_end1_p1  * LooseEle1_Pt       ) * 
			     ( eFakeRate_low_end1_p1  * LooseEle1_Pt       )) + 
			    (( eFakeRate_low_end1_p0  * eFakeRate_low_end1_p0 ) * 
			     ( eFakeRate_low_end1_p0  * eFakeRate_low_end1_p0 )));

       }
       if ( ele1_isEndcap2 ) {
	 fakeRate1 = fakeRate_low_end2_p0 + fakeRate_low_end2_p1 * LooseEle1_Pt;
	 eFakeRate1 = sqrt ((( eFakeRate_low_end2_p1  * LooseEle1_Pt       ) * 
			     ( eFakeRate_low_end2_p1  * LooseEle1_Pt       )) + 
			    (( eFakeRate_low_end2_p0  * eFakeRate_low_end2_p0 ) * 
			     ( eFakeRate_low_end2_p0  * eFakeRate_low_end2_p0 )));
       } 
     }
     else if  ( LooseEle1_Pt >= break_point_ele1 ){
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

     if ( LooseEle2_Pt < break_point_ele2 ){
       if ( ele1_isBarrel  ) {
	 fakeRate2  = fakeRate_low_bar_p0  + fakeRate_low_bar_p1  * LooseEle2_Pt;
	 eFakeRate2 = sqrt ((( eFakeRate_low_bar_p1  * LooseEle2_Pt       ) * 
			     ( eFakeRate_low_bar_p1  * LooseEle2_Pt       )) + 
			    (( eFakeRate_low_bar_p0  * eFakeRate_low_bar_p0 ) * 
			     ( eFakeRate_low_bar_p0  * eFakeRate_low_bar_p0 )));
       }
       if ( ele1_isEndcap1 ) {
	 fakeRate2 = fakeRate_low_end1_p0 + fakeRate_low_end1_p1 * LooseEle2_Pt;
	 eFakeRate2 = sqrt ((( eFakeRate_low_end1_p1  * LooseEle2_Pt       ) * 
			     ( eFakeRate_low_end1_p1  * LooseEle2_Pt       )) + 
			    (( eFakeRate_low_end1_p0  * eFakeRate_low_end1_p0 ) * 
			     ( eFakeRate_low_end1_p0  * eFakeRate_low_end1_p0 )));
       }
       if ( ele1_isEndcap2 ) {
	 fakeRate2 = fakeRate_low_end2_p0 + fakeRate_low_end2_p1 * LooseEle2_Pt;
	 eFakeRate2 = sqrt ((( eFakeRate_low_end2_p1  * LooseEle2_Pt       ) * 
			     ( eFakeRate_low_end2_p1  * LooseEle2_Pt       )) + 
			    (( eFakeRate_low_end2_p0  * eFakeRate_low_end2_p0 ) * 
			     ( eFakeRate_low_end2_p0  * eFakeRate_low_end2_p0 )));
       }
     }

     else if  ( LooseEle2_Pt >= break_point_ele2 ){
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
     // How many loose electrons have HEEP ID?
     //--------------------------------------------------------------------------

     int nPass = 0;
     if ( LooseEle1_PassID == 1 ) nPass ++;
     if ( LooseEle2_PassID == 1 ) nPass ++;

     //--------------------------------------------------------------------------
     // Calculate a few missing variables
     //--------------------------------------------------------------------------

     TLorentzVector loose_ele1, loose_ele2 , jet1, jet2;
     loose_ele1.SetPtEtaPhiM ( LooseEle1_Pt , LooseEle1_Eta , LooseEle1_Phi , 0.0 );
     loose_ele2.SetPtEtaPhiM ( LooseEle2_Pt , LooseEle2_Eta , LooseEle2_Phi , 0.0 );
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


     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "PassJSON"                      , passedJSON ) ; 

     // Noise filters
     fillVariableWithValue(   "PassHBHENoiseFilter"           , PassHBHENoiseFilter ) ; 
     fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight ) ; 

     //--------------------------------------------------------------------------							            
     // MET filters
     //--------------------------------------------------------------------------
           
     fillVariableWithValue(   "PassHBHENoiseFilter"	      , PassHBHENoiseFilter                              , pileup_weight * min_prescale  * fakeRateEffective);
     fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight                          , pileup_weight * min_prescale  * fakeRateEffective);
     fillVariableWithValue(   "PassBadEESupercrystalFilter"   , ( isData == 1 ) ? PassBadEESupercrystalFilter : 1, pileup_weight * min_prescale  * fakeRateEffective);
     fillVariableWithValue(   "PassBeamScraping"	      , ( isData == 1 ) ? PassBeamScraping	      : 1, pileup_weight * min_prescale  * fakeRateEffective);
     fillVariableWithValue(   "PassEcalDeadCellBoundEnergy"   , PassEcalDeadCellBoundEnergy                      , pileup_weight * min_prescale  * fakeRateEffective);
     fillVariableWithValue(   "PassEcalDeadCellTrigPrim"      , PassEcalDeadCellTrigPrim                         , pileup_weight * min_prescale  * fakeRateEffective);
     fillVariableWithValue(   "PassEcalLaserCorrFilter"       , ( isData == 1 ) ? PassEcalLaserCorrFilter     : 1, pileup_weight * min_prescale  * fakeRateEffective);
     fillVariableWithValue(   "PassPhysDecl"		      , ( isData == 1 ) ? PassPhysDecl		      : 1, pileup_weight * min_prescale  * fakeRateEffective);
     fillVariableWithValue(   "PassPrimaryVertex"	      , PassPrimaryVertex                                , pileup_weight * min_prescale  * fakeRateEffective);
     fillVariableWithValue(   "PassTrackingFailure"	      , ( isData == 1 ) ? PassTrackingFailure	      : 1, pileup_weight * min_prescale  * fakeRateEffective);
     
     
     fillVariableWithValue ( "PassHLT", passTrigger, pileup_weight * min_prescale * fakeRateEffective ) ;
     fillVariableWithValue ( "PFMET"  , PFMET_Type01XY_Pt, pileup_weight * min_prescale * fakeRateEffective ) ;
     
     // Electrons
     fillVariableWithValue(   "nEleLoose"                     , nLooseEle_ptCut     , pileup_weight * min_prescale  * fakeRateEffective );
     fillVariableWithValue(   "nEleTight"                     , nPass            , pileup_weight * min_prescale  * fakeRateEffective );
     if ( nLooseEle_ptCut >= 1 ) { 
       fillVariableWithValue( "Ele1_Pt"                       , LooseEle1_Pt   , pileup_weight * min_prescale  * fakeRateEffective ) ;
       fillVariableWithValue( "Ele1_Eta"                      , LooseEle1_Eta  , pileup_weight * min_prescale  * fakeRateEffective ) ;
     }
     if ( nLooseEle_ptCut >= 2 ) { 
       fillVariableWithValue( "Ele2_Pt"                       , LooseEle2_Pt   , pileup_weight * min_prescale  * fakeRateEffective ) ;
       fillVariableWithValue( "Ele2_Eta"                      , LooseEle2_Eta  , pileup_weight * min_prescale  * fakeRateEffective ) ;
       fillVariableWithValue( "M_e1e2"                        , loose_e1e2.M()   , pileup_weight * min_prescale  * fakeRateEffective );
       fillVariableWithValue( "Pt_e1e2"                       , loose_e1e2.Pt()  , pileup_weight * min_prescale  * fakeRateEffective );
     }

     // Jets
     fillVariableWithValue(   "nJet"                          , nJetLooseEle_ptCut , pileup_weight * min_prescale  * fakeRateEffective );
     if ( nJetLooseEle_store >= 1 ) {
       fillVariableWithValue( "Jet1_Pt"                       , JetLooseEle1_Pt   , pileup_weight * min_prescale  * fakeRateEffective ) ;
       fillVariableWithValue( "Jet1_Eta"                      , JetLooseEle1_Eta , pileup_weight * min_prescale  * fakeRateEffective ) ;
     }

     // Muons
     fillVariableWithValue(   "nMuon"                         , nMuon_ptCut           , pileup_weight * min_prescale  * fakeRateEffective ); 
										      			      
     // DeltaR
     if ( nLooseEle_ptCut >= 2 && nJetLooseEle_store >= 1) {
       double st = LooseEle1_Pt + LooseEle2_Pt + JetLooseEle1_Pt ;
       fillVariableWithValue( "DR_Ele1Jet1"                   , loose_ele1.DeltaR ( jet1 ), pileup_weight * min_prescale  * fakeRateEffective );
       fillVariableWithValue( "DR_Ele2Jet1"                   , loose_ele2.DeltaR ( jet1 ), pileup_weight * min_prescale  * fakeRateEffective );
       fillVariableWithValue( "sT_eej_200"                    , st                        , pileup_weight * min_prescale  * fakeRateEffective );
       fillVariableWithValue( "sT_eej_450"                    , st                        , pileup_weight * min_prescale  * fakeRateEffective );
       fillVariableWithValue( "sT_eej_850"                    , st                        , pileup_weight * min_prescale  * fakeRateEffective );
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

     bool passed_preselection = ( passedAllPreviousCuts("sT_eej_200") && passedCut ("sT_eej_200") );
     
     if ( passed_preselection ) {

       double sT_eej = LooseEle1_Pt + LooseEle2_Pt + JetLooseEle1_Pt ;

       FillUserTH1D("nElectron_PAS"        , nLooseEle_ptCut           , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("nMuon_PAS"            , nMuon_ptCut               , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("nJet_PAS"             , nJetLooseEle_ptCut        , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Pt1stEle_PAS"	   , LooseEle1_Pt              , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Eta1stEle_PAS"	   , LooseEle1_Eta             , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Phi1stEle_PAS"	   , LooseEle1_Phi             , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Pt2ndEle_PAS"	   , LooseEle2_Pt              , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Eta2ndEle_PAS"	   , LooseEle2_Eta             , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Phi2ndEle_PAS"	   , LooseEle2_Phi             , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Charge1stEle_PAS"	   , LooseEle1_Charge          , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Charge2ndEle_PAS"	   , LooseEle2_Charge          , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("MET_PAS"              , PFMET_Type01XY_Pt         , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("METPhi_PAS"	   , PFMET_Type01XY_Phi        , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Pt1stJet_PAS"         , JetLooseEle1_Pt           , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Eta1stJet_PAS"        , JetLooseEle1_Eta          , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Phi1stJet_PAS"	   , JetLooseEle1_Phi          , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("sT_PAS"               , sT_eej                    , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Mee_PAS"		   , loose_e1e2.M()            , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Ptee_PAS"             , loose_e1e2.Pt()           , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("nVertex_PAS"          , nVertex                   , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("DR_Ele1Jet1_PAS"	   , loose_ele1.DeltaR ( jet1 ), pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("DR_Ele2Jet1_PAS"	   , loose_ele2.DeltaR ( jet1 ), pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Me1j1_PAS"            , e1j1.M()                  , pileup_weight * min_prescale * fakeRateEffective );
       FillUserTH1D("Me2j1_PAS"            , e2j1.M()                  , pileup_weight * min_prescale * fakeRateEffective );
     }
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

