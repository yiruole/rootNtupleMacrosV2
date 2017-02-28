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
   
   fillSkim                         ( true  ) ;
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
   CreateUserTH1D( "Pt2ndJet_PAS"          ,    100 , 0       , 1000	 ); 
   CreateUserTH1D( "Eta1stJet_PAS"         ,    100 , -5      , 5	 ); 
   CreateUserTH1D( "Eta2ndJet_PAS"         ,    100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Phi2ndJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "sTlep_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTjet_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PAS"                ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PASandMee100"       ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mjj_PAS"		   ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mee_PAS"		   ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mee_PASandST445"       ,    200 , 0       , 2000	 ); 
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
		                           
   CreateUserTH1D( "nVertex_PAS"           ,    101   , -0.5   , 100.5	 ) ; 
		                           
   CreateUserTH1D( "DR_Ele1Jet1_PAS"	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele1Jet2_PAS"	   , 	getHistoNBins("DR_Ele1Jet2"), getHistoMin("DR_Ele1Jet2"), getHistoMax("DR_Ele1Jet2")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet1_PAS"	   , 	getHistoNBins("DR_Ele2Jet1"), getHistoMin("DR_Ele2Jet1"), getHistoMax("DR_Ele2Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet2_PAS"	   , 	getHistoNBins("DR_Ele2Jet2"), getHistoMin("DR_Ele2Jet2"), getHistoMax("DR_Ele2Jet2")     ) ; 
   CreateUserTH1D( "DR_Jet1Jet2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
   CreateUserTH1D( "DR_Ele1Ele2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
   CreateUserTH1D( "minDR_EleJet_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 

   CreateUserTH2D( "Me1jVsMe2j_selected",     200, 0, 2000, 200, 0, 2000) ;
   CreateUserTH2D( "Me1jVsMe2j_rejected",     200, 0, 2000, 200, 0, 2000) ;
   CreateUserTH2D( "sT_vs_Mee"          ,     200, 250, 750, 200, 50, 450) ;

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
     
     int passedJSON = passJSON ( run, ls , isData ) ;
     
     //--------------------------------------------------------------------------
     // Find the right prescale for this event
     //--------------------------------------------------------------------------
     
     int min_prescale = 0;
     int passTrigger  = 0;

     if ( LooseEle1_hltPhotonPt > 0.0 ) { 
       if ( H_Photon30_CIdVL > 0.1 && LooseEle1_hltPhotonPt >= 30.  && LooseEle1_hltPhotonPt < 50. ) { passTrigger = 1; min_prescale = H_Photon30_CIdVL; } 
       if ( H_Photon50_CIdVL > 0.1 && LooseEle1_hltPhotonPt >= 50.  && LooseEle1_hltPhotonPt < 75. ) { passTrigger = 1; min_prescale = H_Photon50_CIdVL; } 
       if ( H_Photon75_CIdVL > 0.1 && LooseEle1_hltPhotonPt >= 75.  && LooseEle1_hltPhotonPt < 90. ) { passTrigger = 1; min_prescale = H_Photon75_CIdVL; } 
       if ( H_Photon90_CIdVL > 0.1 && LooseEle1_hltPhotonPt >= 90.  && LooseEle1_hltPhotonPt < 135.) { passTrigger = 1; min_prescale = H_Photon90_CIdVL; } 
       if ( H_Photon135      > 0.1 && LooseEle1_hltPhotonPt >= 135. && LooseEle1_hltPhotonPt < 150.) { passTrigger = 1; min_prescale = H_Photon135     ; } 
       if ( H_Photon150      > 0.1 && LooseEle1_hltPhotonPt >= 150.                                ) { passTrigger = 1; min_prescale = H_Photon150     ; } 
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
     
     if ( LooseEle1_Pt < 100 ){
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
     else if  ( LooseEle1_Pt >= 100 ){
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

     if ( LooseEle2_Pt < 100 ){
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

     else if  ( LooseEle2_Pt >= 100 ){
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
     loose_ele1.SetPtEtaPhiM ( LooseEle1_Pt , LooseEle1_Eta , LooseEle1_Phi , 0.0 );
     loose_ele2.SetPtEtaPhiM ( LooseEle2_Pt , LooseEle2_Eta , LooseEle2_Phi , 0.0 );
     jet1.SetPtEtaPhiM       ( JetLooseEle1_Pt, JetLooseEle1_Eta, JetLooseEle1_Phi, 0.0 );
     jet2.SetPtEtaPhiM       ( JetLooseEle2_Pt, JetLooseEle2_Eta, JetLooseEle2_Phi, 0.0 );
     met.SetPtEtaPhiM        ( PFMET_Type01XY_Pt, 0.0             , PFMET_Type01XY_Phi , 0.0 );
     muon.SetPtEtaPhiM       ( Muon1_Pt         , Muon1_Eta       , Muon1_Phi          , 0.0 );

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

     MT_Ele1MET = sqrt(2 * LooseEle1_Pt * met.Pt()  * (1 - cos(loose_ele1.DeltaPhi ( met) ) ) );

     mDPhi_METEle1= fabs(loose_ele1.DeltaPhi ( met ));
     mDPhi_METJet1= fabs(jet1.DeltaPhi ( met ));
     mDPhi_METJet2= fabs(jet2.DeltaPhi ( met ));

     sT_eejj = LooseEle1_Pt + LooseEle2_Pt + JetLooseEle1_Pt + JetLooseEle2_Pt ;
     
     int nEle_Loose_Ana = 0;
     if      ( LooseEle3_Pt > 35. ) nEle_Loose_Ana = 3;
     else if ( LooseEle2_Pt > 35. ) nEle_Loose_Ana = 2;
     else if ( LooseEle1_Pt > 35. ) nEle_Loose_Ana = 1;
     else                             nEle_Loose_Ana = 0;
     
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

     // Noise filters
     fillVariableWithValue(   "PassHBHENoiseFilter"	      , PassHBHENoiseFilter                              , fakeRateEffective * min_prescale );
     fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight                          , fakeRateEffective * min_prescale );
     fillVariableWithValue(   "PassBadEESupercrystalFilter"   , ( isData == 1 ) ? PassBadEESupercrystalFilter : 1, fakeRateEffective * min_prescale );
     fillVariableWithValue(   "PassBeamScraping"	      , ( isData == 1 ) ? PassBeamScraping	      : 1, fakeRateEffective * min_prescale );
     fillVariableWithValue(   "PassEcalDeadCellBoundEnergy"   , PassEcalDeadCellBoundEnergy                      , fakeRateEffective * min_prescale );
     fillVariableWithValue(   "PassEcalDeadCellTrigPrim"      , PassEcalDeadCellTrigPrim                         , fakeRateEffective * min_prescale );
     fillVariableWithValue(   "PassEcalLaserCorrFilter"       , ( isData == 1 ) ? PassEcalLaserCorrFilter     : 1, fakeRateEffective * min_prescale );
     fillVariableWithValue(   "PassHcalLaserEventFilter"      , ( isData == 1 ) ? PassHcalLaserEventFilter    : 1, fakeRateEffective * min_prescale );
     fillVariableWithValue(   "PassPhysDecl"		      , ( isData == 1 ) ? PassPhysDecl		      : 1, fakeRateEffective * min_prescale );
     fillVariableWithValue(   "PassPrimaryVertex"	      , PassPrimaryVertex                                , fakeRateEffective * min_prescale );
     fillVariableWithValue(   "PassTrackingFailure"	      , ( isData == 1 ) ? PassTrackingFailure	      : 1, fakeRateEffective * min_prescale );
     
     // Electrons
     int PassNEle = 0;
     if ( nLooseEle_ptCut == 2 ) PassNEle = 1;
     
     double M_ej_avg;
     double M_ej_min;
     double M_ej_max;

     // Muons
     int PassNMuon = 0;
     if ( nMuon_ptCut == 0 ) PassNMuon = 1;

     fillVariableWithValue ( "PassHLT"                        , passTrigger             , fakeRateEffective * min_prescale ) ;
     										        
     // Electrons								        
     fillVariableWithValue(   "PassNEle"                      , PassNEle                , fakeRateEffective * min_prescale ) ;
     if ( nLooseEle_store >= 1 ) { 							        
       fillVariableWithValue( "Ele1_Pt"                       , LooseEle1_Pt            , fakeRateEffective * min_prescale ) ;
     }										        
     if ( nLooseEle_store >= 2 ) { 							        
       fillVariableWithValue( "Ele2_Pt"                       , LooseEle2_Pt            , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "M_e1e2"                        , M_e1e2                  , fakeRateEffective * min_prescale ) ;
     }

     // Jets
     fillVariableWithValue(   "nJet"                          , nJetLooseEle_ptCut      , fakeRateEffective * min_prescale ) ;
     if ( nJetLooseEle_store >= 1 ) { 						                
       fillVariableWithValue( "Jet1_Pt"                       , JetLooseEle1_Pt         , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "Jet1_Eta"                      , JetLooseEle1_Eta        , fakeRateEffective * min_prescale ) ;
     }
     
     if ( nJetLooseEle_store >= 2 ) { 
       fillVariableWithValue( "Jet2_Pt"                       , JetLooseEle2_Pt         , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "Jet2_Eta"                      , JetLooseEle2_Eta        , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "DR_Jet1Jet2"                   , DR_Jet1Jet2             , fakeRateEffective * min_prescale ) ;



       if ( nLooseEle_store >= 2 && nJetLooseEle_store >= 2) {
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
     }

     // Muons
     fillVariableWithValue(   "PassNMuon"                     , PassNMuon               , fakeRateEffective * min_prescale ) ;

     // DeltaR
     if ( nLooseEle_store >= 2 && nJetLooseEle_store >= 1) {
       fillVariableWithValue( "DR_Ele1Jet1"                   , DR_Ele1Jet1             , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "DR_Ele2Jet1"                   , DR_Ele2Jet1             , fakeRateEffective * min_prescale ) ;
       if(nJetLooseEle_store >= 2) {
	 fillVariableWithValue( "DR_Ele1Jet2"                 , DR_Ele1Jet2             , fakeRateEffective * min_prescale ) ;
	 fillVariableWithValue( "DR_Ele2Jet2"                 , DR_Ele2Jet2             , fakeRateEffective * min_prescale ) ;
       }
     }

     // sT
     if ( nLooseEle_store >= 2 && nJetLooseEle_store >= 2) {
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

       //--------------------------------------------------------------------------
       // Fill skim tree, if necessary
       //--------------------------------------------------------------------------
       
       fillSkimTree();
       
       double DR_Ele1Jet3 = 999.;
       double DR_Ele2Jet3 = 999.;
       double min_DR_EleJet = 999.;
       TLorentzVector eejj, e1e2mu, e1, j1, e2, j2,j3, mu, met;
       e1.SetPtEtaPhiM ( LooseEle1_Pt, LooseEle1_Eta, LooseEle1_Phi, 0.0 );
       e2.SetPtEtaPhiM ( LooseEle2_Pt, LooseEle2_Eta, LooseEle2_Phi, 0.0 );
       j1.SetPtEtaPhiM ( JetLooseEle1_Pt, JetLooseEle1_Eta, JetLooseEle1_Phi, 0.0 );
       j2.SetPtEtaPhiM ( JetLooseEle2_Pt, JetLooseEle2_Eta, JetLooseEle2_Phi, 0.0 );
       if ( nJetLooseEle_store > 2 ) {
	 j3.SetPtEtaPhiM ( JetLooseEle3_Pt, JetLooseEle3_Eta, JetLooseEle3_Phi, 0.0 );
	 DR_Ele1Jet3 = e1.DeltaR( j3 );
	 DR_Ele2Jet3 = e2.DeltaR( j3 );
       }
       
       
       if ( DR_Ele1Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet1;
       if ( DR_Ele1Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet2;
       if ( DR_Ele2Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet1;
       if ( DR_Ele2Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet2;
       if ( nJetLooseEle_store > 2 ) {
	 if ( DR_Ele1Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet3;
	 if ( DR_Ele2Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet3;
       }
       
       FillUserTH1D( "minDR_EleJet_PAS", min_DR_EleJet, min_prescale * fakeRateEffective ) ;
       
       FillUserTH1D("nElectron_PAS"        , nLooseEle_ptCut                   , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("nMuon_PAS"            , nMuon_ptCut                       , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("nJet_PAS"             , nJetLooseEle_ptCut                , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Pt1stEle_PAS"	   , LooseEle1_Pt                      , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Eta1stEle_PAS"	   , LooseEle1_Eta                     , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Phi1stEle_PAS"	   , LooseEle1_Phi                     , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Pt2ndEle_PAS"	   , LooseEle2_Pt                      , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Eta2ndEle_PAS"	   , LooseEle2_Eta                     , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Phi2ndEle_PAS"	   , LooseEle2_Phi                     , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Charge1stEle_PAS"	   , LooseEle1_Charge                  , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Charge2ndEle_PAS"	   , LooseEle2_Charge                  , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("MET_PAS"              , PFMET_Type01XY_Pt                 , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("METPhi_PAS"	   , PFMET_Type01XY_Phi                , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Pt1stJet_PAS"         , JetLooseEle1_Pt                   , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Pt2ndJet_PAS"         , JetLooseEle2_Pt                   , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Eta1stJet_PAS"        , JetLooseEle1_Eta                  , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Eta2ndJet_PAS"        , JetLooseEle2_Eta                  , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Phi1stJet_PAS"	   , JetLooseEle1_Phi                  , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Phi2ndJet_PAS"	   , JetLooseEle2_Phi                  , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("sTlep_PAS"            , LooseEle1_Pt + LooseEle2_Pt       , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("sTjet_PAS"            , JetLooseEle1_Pt + JetLooseEle2_Pt , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("sT_PAS"               , sT_eejj                           , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Mjj_PAS"		   , M_j1j2                            , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Mee_PAS"		   , M_e1e2                            , min_prescale * fakeRateEffective ) ;
       FillUserTH1D( "MTenu_PAS"           , MT_Ele1MET                        , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Me1j1_PAS"            , M_e1j1                            , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Me1j2_PAS"            , M_e1j2                            , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Me2j1_PAS"            , M_e2j1                            , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Me2j2_PAS"            , M_e2j2                            , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Ptee_PAS"             , Pt_e1e2                           , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta               , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Dist1stEle_PAS"       , LooseEle1_Dist                    , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DCotTheta2ndEle_PAS"  , LooseEle2_DCotTheta               , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Dist2ndEle_PAS"       , LooseEle2_Dist                    , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("nVertex_PAS"          , nVertex                           , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DR_Ele1Jet1_PAS"	   , DR_Ele1Jet1                       , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DR_Ele1Jet2_PAS"	   , DR_Ele1Jet2                       , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DR_Ele2Jet1_PAS"	   , DR_Ele2Jet1                       , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DR_Ele2Jet2_PAS"	   , DR_Ele2Jet2                       , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DR_Jet1Jet2_PAS"	   , DR_Jet1Jet2                       , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DR_Ele1Ele2_PAS"	   , DR_Ele1Ele2                       , min_prescale * fakeRateEffective ) ;
       FillUserTH2D("sT_vs_Mee"            , sT_eejj, M_e1e2                   , min_prescale * fakeRateEffective ) ;

       if ( sT_eejj > 445. ) { 
	 FillUserTH1D( "Mee_PASandST445", M_e1e2, min_prescale * fakeRateEffective ) ;
       }
       
       if ( M_e1e2 > 100. ) {
	 FillUserTH1D("sT_PASandMee100", sT_eejj, min_prescale * fakeRateEffective ) ;
       }

       if ( M_e1e2 > 80.0 && M_e1e2 < 100.0 ){
	 FillUserTH1D("Mee_80_100_Preselection", M_e1e2, min_prescale * fakeRateEffective ) ;
	 if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_80_100_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
	 else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_80_100_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
	 else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_80_100_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
	 if      ( isEB   ) FillUserTH1D( "Mee_EB_80_100_PAS"  , M_e1e2, min_prescale * fakeRateEffective  ); 
       }

       if ( M_e1e2 > 70.0 && M_e1e2 < 110.0 ){
	 FillUserTH1D("Mee_70_110_Preselection", M_e1e2, min_prescale * fakeRateEffective ) ;
	 if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_70_110_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
	 else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_70_110_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
	 else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_70_110_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
	 if      ( isEB   ) FillUserTH1D( "Mee_EB_70_110_PAS"  , M_e1e2, min_prescale * fakeRateEffective  ); 
       }
       
       if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
       else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
       else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
       if      ( isEB   ) FillUserTH1D( "Mee_EB_PAS"  , M_e1e2, min_prescale * fakeRateEffective  ); 
       
       FillUserTH1D("Meejj_PAS", M_eejj , min_prescale * fakeRateEffective  );

       if ( fabs(M_e1j1-M_e2j2) < fabs(M_e1j2-M_e2j1) )  {
	 FillUserTH1D("Mej_selected_avg_PAS", M_ej_avg, min_prescale * fakeRateEffective ) ;	   
	 FillUserTH1D("Me1j_selected_PAS"   , M_e1j1, min_prescale * fakeRateEffective ) ;	   
	 FillUserTH1D("Me2j_selected_PAS"   , M_e2j2, min_prescale * fakeRateEffective ) ;	   
	 FillUserTH2D( "Me1jVsMe2j_selected", M_e1j1, M_e2j2, min_prescale * fakeRateEffective  ) ;
	 FillUserTH2D( "Me1jVsMe2j_rejected", M_e1j2, M_e2j1, min_prescale * fakeRateEffective  ) ;
       }
       else {
	 FillUserTH1D("Mej_selected_avg_PAS", M_ej_avg, min_prescale * fakeRateEffective ) ;	   
	 FillUserTH1D("Me1j_selected_PAS" , M_e1j2, min_prescale * fakeRateEffective ) ;	   
	 FillUserTH1D("Me2j_selected_PAS" , M_e2j1, min_prescale * fakeRateEffective ) ;	   
	 FillUserTH2D( "Me1jVsMe2j_selected", M_e1j2, M_e2j1, min_prescale * fakeRateEffective  ) ;
	 FillUserTH2D( "Me1jVsMe2j_rejected", M_e1j1, M_e2j2, min_prescale * fakeRateEffective  ) ;
       }
     }
   } // End loop over events
   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

