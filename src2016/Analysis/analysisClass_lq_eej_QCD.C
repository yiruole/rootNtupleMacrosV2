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

   CreateUserTH1D( "sTfrac_Jet1_PAS"       ,   100  ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele1_PAS"       ,   100  ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele2_PAS"       ,   100  ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Jet_PAS"        ,   100  ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele_PAS"        ,   100  ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Jet1_PASandMee100"       ,   100  ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele1_PASandMee100"       ,   100  ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele2_PASandMee100"       ,   100  ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Jet_PASandMee100"        ,   100  ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele_PASandMee100"        ,   100  ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Jet1_ROI"       ,   100  ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele1_ROI"       ,   100  ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele2_ROI"       ,   100  ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Jet_ROI"        ,   100  ,  0.0    , 1.0      );
   CreateUserTH1D( "sTfrac_Ele_ROI"        ,   100  ,  0.0    , 1.0      );
   CreateUserTH1D( "nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "nMuon_PAS"             ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "nJet_PAS"              ,    10  , -0.5    , 9.5      );
   CreateUserTH1D( "nJet_PASandMee100"        ,    10  , -0.5    , 9.5      );
   CreateUserTH1D( "nJet_ROI"              ,    10  , -0.5    , 9.5      );
   CreateUserTH1D( "Pt1stEle_PAS"	   , 	100 , 0       , 1000     ); 
   CreateUserTH1D( "Pt1stEle_PASandMee100" , 	100 , 0       , 1000     ); 
   CreateUserTH1D( "Pt1stEle_ROI"	   , 	100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta1stEle_PAS"	   , 	100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stEle_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Pt2ndEle_PAS"	   , 	100 , 0       , 1000     ); 
   CreateUserTH1D( "Pt2ndEle_PASandMee100" , 	100 , 0       , 1000     ); 
   CreateUserTH1D( "Pt2ndEle_ROI"	   , 	100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta2ndEle_PAS"	   , 	100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi2ndEle_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Charge1stEle_PAS"	   , 	2   , -1.0001 , 1.0001	 ); 
   CreateUserTH1D( "Charge2ndEle_PAS"	   , 	2   , -1.0001 , 1.0001	 ); 
   CreateUserTH1D( "EleChargeSum_PAS"         ,    3   , -2.5    , 2.5  );
   CreateUserTH1D( "EleChargeSum_PASandMee100",    3   , -2.5    , 2.5  );
   CreateUserTH1D( "EleChargeSum_ROI"         ,    3   , -2.5    , 2.5  );
   CreateUserTH1D( "MET_PAS"               ,    200 , 0       , 1000	 ); 
   CreateUserTH1D( "MET_ROI"               ,    200 , 0       , 1000	 ); 
   CreateUserTH1D( "METPhi_PAS"		   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Pt1stJet_PAS"          ,    100 , 0       , 1000	 ); 
   CreateUserTH1D( "Pt1stJet_PASandMee100" ,    100 , 0       , 1000	 ); 
   CreateUserTH1D( "Pt1stJet_ROI"          ,    100 , 0       , 1000	 ); 
   CreateUserTH1D( "Eta1stJet_PAS"         ,    100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "sTlep_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTlep_PASandMee100"    ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTlep_ROI"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PAS"                ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PASandMee100"       ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PASandMee110"       ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PASandMee120"       ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PASandMee130"       ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PASandMee140"       ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PASandMee150"       ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PASandMee160"       ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PASandMee170"       ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PASandMee180"       ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PASandMee190"       ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PASandMee200"       ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_ROI"                ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mee_PAS"		   ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mee_ROI"		   ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mee_PASandST445"       ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "MTenu_PAS"             ,    200 , 0       , 1000	 ); 
   CreateUserTH1D( "Me1j1_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me1j1_PASandMee100"    ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me1j1_ROI"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me2j1_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me1j_selected_PAS"     ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Meej_PAS"              ,    400 , 0       , 4000     );
   CreateUserTH1D( "Meej_ROI"              ,    400 , 0       , 4000     );

   CreateUserTH1D( "Eta1stJet_ROI"                   ,    100   , -5      , 5	  ); 
   CreateUserTH1D( "Eta1stEle_ROI"	             , 	100    , -5      , 5	  ); 
   CreateUserTH1D( "Eta2ndEle_ROI"	             , 	100    , -5      , 5	  ); 

   CreateUserTH1D( "Phi1stJet_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
   CreateUserTH1D( "Phi1stEle_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
   CreateUserTH1D( "Phi2ndEle_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
		            			     
   CreateUserTH1D( "Ptee_PAS"              ,    200 , 0       , 2000     );
   CreateUserTH1D( "Ptee_PASandMee100"     ,    200 , 0       , 2000     );
   CreateUserTH1D( "Ptee_ROI"              ,    200 , 0       , 2000     );
   		                           
   CreateUserTH1D( "DCotTheta1stEle_PAS"   ,    100 , 0.0, 1.0);
   CreateUserTH1D( "Dist1stEle_PAS"        ,    100 , 0.0, 1.0);  
   CreateUserTH1D( "DCotTheta2ndEle_PAS"   ,    100 , 0.0, 1.0);
   CreateUserTH1D( "Dist2ndEle_PAS"        ,    100 , 0.0, 1.0);  
		                           
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


   CreateUserTH1D( "Mee_80_100_Preselection", 200, 60, 120 );
   CreateUserTH1D( "Mee_70_110_Preselection", 200, 60, 120 );
   CreateUserTH1D( "Mee_70_110_ST600_Preselection", 200, 60, 120 );
   
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
     fillVariableWithValue("nEle_hltMatched",-1, fakeRateEffective * min_prescale ) ;
     fillVariableWithValue("nJet_hltMatched",-1, fakeRateEffective * min_prescale ) ;
     			
							        
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

     double sT_eej = LooseEle1_Pt + LooseEle2_Pt + JetLooseEle1_Pt;

     // Muons
     fillVariableWithValue(   "PassNMuon"                     , PassNMuon               , fakeRateEffective * min_prescale ) ;

     // DeltaR
     if ( nLooseEle_store >= 2 && nJetLooseEle_store >= 1) {
       fillVariableWithValue( "DR_Ele1Jet1"                   , DR_Ele1Jet1             , fakeRateEffective * min_prescale ) ;
       fillVariableWithValue( "DR_Ele2Jet1"                   , DR_Ele2Jet1             , fakeRateEffective * min_prescale ) ;
     }

     // sT
     if ( nLooseEle_store >= 2 && nJetLooseEle_store >= 2) {
       fillVariableWithValue( "sT_eej"                        , sT_eej                  , fakeRateEffective  * min_prescale ) ;
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

     FillUserTH1D( "PileupWeight"   , -1.0 );
     FillUserTH1D( "GeneratorWeight", -1.0 );

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
       
       e1.SetPtEtaPhiM ( LooseEle1_Pt, LooseEle1_Eta, LooseEle1_Phi, 0.0 );
       e2.SetPtEtaPhiM ( LooseEle2_Pt, LooseEle2_Eta, LooseEle2_Phi, 0.0 );
       j1.SetPtEtaPhiM ( JetLooseEle1_Pt, JetLooseEle1_Eta, JetLooseEle1_Phi, 0.0 );
       mu.SetPtEtaPhiM ( Muon1_Pt, Muon1_Eta, Muon1_Phi, 0.0 );
       met.SetPtEtaPhiM ( PFMET_Type01XY_Pt, 0.0, PFMET_Type01XY_Phi, 0.0 );

       eej  = e1 + e2 + j1;
       ee   = e1 + e2;
       
       double DR_Ele1Ele2 = e1.DeltaR( e2 );       
       double M_eej  = eej.M();
       double DR_ZJ1 = ee.DeltaR ( j1 );
       
       //--------------------------------------------------------------------------
       // Preselection histograms
       //--------------------------------------------------------------------------

       FillUserTH1D("DR_Ele1Ele2_PAS"	   , DR_Ele1Ele2                        , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("EleChargeSum_PAS"     , LooseEle1_Charge + LooseEle2_Charge, min_prescale * fakeRateEffective ) ;
       FillUserTH1D("sT_PAS"               , sT_eejj                            , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("nElectron_PAS"        , nLooseEle_ptCut                    , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("nMuon_PAS"            , nMuon_ptCut                        , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("nJet_PAS"             , nJetLooseEle_ptCut                 , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Pt1stEle_PAS"	   , LooseEle1_Pt                       , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Eta1stEle_PAS"	   , LooseEle1_Eta                      , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Phi1stEle_PAS"	   , LooseEle1_Phi                      , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Pt2ndEle_PAS"	   , LooseEle2_Pt                       , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Eta2ndEle_PAS"	   , LooseEle2_Eta                      , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Phi2ndEle_PAS"	   , LooseEle2_Phi                      , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Charge1stEle_PAS"	   , LooseEle1_Charge                   , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Charge2ndEle_PAS"	   , LooseEle2_Charge                   , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("MET_PAS"              , PFMET_Type01XY_Pt                  , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("METPhi_PAS"	   , PFMET_Type01XY_Phi                 , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Pt1stJet_PAS"         , JetLooseEle1_Pt                    , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Eta1stJet_PAS"        , JetLooseEle1_Eta                   , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Phi1stJet_PAS"	   , JetLooseEle1_Phi                   , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("sTlep_PAS"            , LooseEle1_Pt + LooseEle2_Pt        , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Mee_PAS"		   , M_e1e2                             , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("MTenu_PAS"            , MT_Ele1MET                         , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Me1j1_PAS"            , M_e1j1                             , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Me2j1_PAS"            , M_e2j1                             , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Ptee_PAS"             , Pt_e1e2                            , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta                , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Dist1stEle_PAS"       , LooseEle1_Dist                     , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DCotTheta2ndEle_PAS"  , LooseEle2_DCotTheta                , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Dist2ndEle_PAS"       , LooseEle2_Dist                     , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("nVertex_PAS"          , nVertex                            , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DR_Ele1Jet1_PAS"	   , DR_Ele1Jet1                        , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DR_Ele2Jet1_PAS"	   , DR_Ele2Jet1                        , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("Meej_PAS"             , M_eej                              , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("DR_ZJet1_PAS"         , DR_ZJ1                             , min_prescale * fakeRateEffective ) ;
       FillUserTH1D("sTfrac_Jet1_PAS"      , JetLooseEle1_Pt / sT_eejj          , min_prescale * fakeRateEffective );
       FillUserTH1D("sTfrac_Ele1_PAS"      , LooseEle1_Pt / sT_eejj             , min_prescale * fakeRateEffective );
       FillUserTH1D("sTfrac_Ele2_PAS"      , LooseEle2_Pt / sT_eejj             , min_prescale * fakeRateEffective );
       FillUserTH1D("sTfrac_Jet_PAS"       , ( JetLooseEle1_Pt  ) / sT_eejj     , min_prescale * fakeRateEffective );
       FillUserTH1D("sTfrac_Ele_PAS"       , ( LooseEle1_Pt + LooseEle2_Pt ) / sT_eejj, min_prescale * fakeRateEffective );

       FillUserTH2D("MeeVsST_PAS" , M_e1e2, sT_eej, min_prescale * fakeRateEffective ) ;	   
       
       //--------------------------------------------------------------------------
       // Preselection + event type (EBEB, EEEB, EEEE, etc)
       //--------------------------------------------------------------------------

       if      ( isEB   ) FillUserTH1D( "Mee_EB_PAS"  , M_e1e2, min_prescale * fakeRateEffective  ); 
       if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
       else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
       else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
       
       //--------------------------------------------------------------------------
       // Preselection + high ST plot
       //--------------------------------------------------------------------------

       if ( sT_eej > 445. ) FillUserTH1D( "Mee_PASandST445", M_e1e2, min_prescale * fakeRateEffective ) ;

       //--------------------------------------------------------------------------
       // High M(ee) plots
       //--------------------------------------------------------------------------

       if ( M_e1e2 > 100. ) FillUserTH1D("sT_PASandMee100"   , sT_eej , min_prescale * fakeRateEffective  ); 
       if ( M_e1e2 > 110. ) FillUserTH1D("sT_PASandMee110"   , sT_eej , min_prescale * fakeRateEffective  ); 
       if ( M_e1e2 > 120. ) FillUserTH1D("sT_PASandMee120"   , sT_eej , min_prescale * fakeRateEffective  ); 
       if ( M_e1e2 > 130. ) FillUserTH1D("sT_PASandMee130"   , sT_eej , min_prescale * fakeRateEffective  ); 
       if ( M_e1e2 > 140. ) FillUserTH1D("sT_PASandMee140"   , sT_eej , min_prescale * fakeRateEffective  ); 
       if ( M_e1e2 > 150. ) FillUserTH1D("sT_PASandMee150"   , sT_eej , min_prescale * fakeRateEffective  ); 
       if ( M_e1e2 > 160. ) FillUserTH1D("sT_PASandMee160"   , sT_eej , min_prescale * fakeRateEffective  ); 
       if ( M_e1e2 > 170. ) FillUserTH1D("sT_PASandMee170"   , sT_eej , min_prescale * fakeRateEffective  ); 
       if ( M_e1e2 > 180. ) FillUserTH1D("sT_PASandMee180"   , sT_eej , min_prescale * fakeRateEffective  ); 
       if ( M_e1e2 > 190. ) FillUserTH1D("sT_PASandMee190"   , sT_eej , min_prescale * fakeRateEffective  ); 
       if ( M_e1e2 > 200. ) FillUserTH1D("sT_PASandMee200"   , sT_eej , min_prescale * fakeRateEffective  ); 
       
       
       if ( M_e1e2 > 100. ) { 

	 FillUserTH1D("Me1j1_PASandMee100"           , M_e1j1                              , min_prescale * fakeRateEffective ) ;
	 FillUserTH1D("Ptee_PASandMee100"            , Pt_e1e2                             , min_prescale * fakeRateEffective ) ;
	 FillUserTH2D("MeeVsST_PASandMee100"         , M_e1e2, sT_eej                      , min_prescale * fakeRateEffective ) ;	   
	 FillUserTH1D("nVertex_PASandMee100"         , nVertex                             , min_prescale * fakeRateEffective ) ;
	 FillUserTH1D("sT_PASandMee100"              , sT_eej                              , min_prescale * fakeRateEffective ) ;
	 FillUserTH1D("EleChargeSum_PASandMee100"    , LooseEle1_Charge + LooseEle2_Charge , min_prescale * fakeRateEffective ) ;
	 FillUserTH1D("nJet_PASandMee100"            , nJetLooseEle_ptCut                  , min_prescale * fakeRateEffective ) ;
	 FillUserTH1D("sTlep_PASandMee100"           , LooseEle1_Pt    + LooseEle2_Pt      , min_prescale * fakeRateEffective ) ;
	 FillUserTH1D("Pt1stEle_PASandMee100"        , LooseEle1_Pt                        , min_prescale * fakeRateEffective ) ;
	 FillUserTH1D("Pt2ndEle_PASandMee100"        , LooseEle2_Pt                        , min_prescale * fakeRateEffective ) ;
	 FillUserTH1D("Pt1stJet_PASandMee100"        , JetLooseEle1_Pt                     , min_prescale * fakeRateEffective ) ;

	 FillUserTH1D("sTfrac_Jet1_PASandMee100"     , JetLooseEle1_Pt / sT_eej                        , min_prescale * fakeRateEffective ) ;
	 FillUserTH1D("sTfrac_Ele1_PASandMee100"     , LooseEle1_Pt / sT_eej                           , min_prescale * fakeRateEffective ) ;
	 FillUserTH1D("sTfrac_Ele2_PASandMee100"     , LooseEle2_Pt / sT_eej                           , min_prescale * fakeRateEffective ) ;
	 FillUserTH1D("sTfrac_Jet_PASandMee100"      , ( JetLooseEle1_Pt )                   / sT_eej  , min_prescale * fakeRateEffective ) ;
	 FillUserTH1D("sTfrac_Ele_PASandMee100"      , ( LooseEle1_Pt + LooseEle2_Pt ) / sT_eej        , min_prescale * fakeRateEffective ) ;

       }

       //--------------------------------------------------------------------------
       // Preselection + M(ee) normalization region plots
       //--------------------------------------------------------------------------
       
       if ( M_e1e2 > 80.0 && M_e1e2 < 100.0 ){
	 FillUserTH1D("Mee_80_100_Preselection", M_e1e2, min_prescale * fakeRateEffective ) ;
	 if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_80_100_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
	 else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_80_100_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
	 else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_80_100_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
	 if      ( isEB   ) FillUserTH1D( "Mee_EB_80_100_PAS"  , M_e1e2, min_prescale * fakeRateEffective  ); 
       }

       if ( M_e1e2 > 70.0 && M_e1e2 < 110.0 ){
	 FillUserTH1D("Mee_70_110_Preselection", M_e1e2, min_prescale * fakeRateEffective ) ;
	 if ( sT_eej > 600 ) 	 FillUserTH1D("Mee_70_110_ST600_Preselection", M_e1e2, min_prescale * fakeRateEffective ) ;
	 if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_70_110_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
	 else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_70_110_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
	 else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_70_110_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
	 if      ( isEB   ) FillUserTH1D( "Mee_EB_70_110_PAS"  , M_e1e2, min_prescale * fakeRateEffective  ); 
       }

       //--------------------------------------------------------------------------
       // Region of interest plots
       //-------------------------------------------------------------------------- 
       
       if ( passed_region_of_interest ) { 
	 FillUserTH1D("Me1j1_ROI"           , M_e1j1                                         , min_prescale * fakeRateEffective );
	 FillUserTH1D("Ptee_ROI"            , Pt_e1e2                                        , min_prescale * fakeRateEffective );
	 FillUserTH1D("Eta1stJet_ROI"       , JetLooseEle1_Eta                               , min_prescale * fakeRateEffective );
	 FillUserTH1D("Eta1stEle_ROI"	    , LooseEle1_Eta                                  , min_prescale * fakeRateEffective );
	 FillUserTH1D("Eta2ndEle_ROI"	    , LooseEle2_Eta                                  , min_prescale * fakeRateEffective );
	 FillUserTH1D("Phi1stJet_ROI"       , JetLooseEle1_Phi                               , min_prescale * fakeRateEffective );
	 FillUserTH1D("Phi1stEle_ROI"	    , LooseEle1_Phi                                  , min_prescale * fakeRateEffective );
	 FillUserTH1D("Phi2ndEle_ROI"	    , LooseEle2_Phi                                  , min_prescale * fakeRateEffective );
	 FillUserTH2D("MeeVsST_ROI"         , M_e1e2                                , sT_eej , min_prescale * fakeRateEffective );	   
	 FillUserTH1D("Mee_ROI"		    , M_e1e2                                         , min_prescale * fakeRateEffective );
	 FillUserTH1D("nVertex_ROI"         , nVertex                                        , min_prescale * fakeRateEffective );
	 FillUserTH1D("EleChargeSum_ROI"    , LooseEle1_Charge + LooseEle2_Charge            , min_prescale * fakeRateEffective );
	 FillUserTH1D("nJet_ROI"            , nJetLooseEle_ptCut                             , min_prescale * fakeRateEffective );
	 FillUserTH1D("Meej_ROI"            , M_eej                                          , min_prescale * fakeRateEffective );
	 FillUserTH1D("DR_ZJet1_ROI"        , DR_ZJ1                                         , min_prescale * fakeRateEffective );
	 FillUserTH1D("MET_ROI"             , PFMET_Type01XY_Pt                              , min_prescale * fakeRateEffective );
	 FillUserTH1D("sT_ROI"              , sT_eej                                         , min_prescale * fakeRateEffective );
	 FillUserTH1D("sTlep_ROI"           , LooseEle1_Pt    + LooseEle2_Pt                 , min_prescale * fakeRateEffective );
	 FillUserTH1D("Pt1stEle_ROI"        , LooseEle1_Pt                                   , min_prescale * fakeRateEffective );
	 FillUserTH1D("Pt2ndEle_ROI"        , LooseEle2_Pt                                   , min_prescale * fakeRateEffective );
	 FillUserTH1D("Pt1stJet_ROI"        , JetLooseEle1_Pt                                , min_prescale * fakeRateEffective );
	 FillUserTH1D("sTfrac_Jet1_ROI"     , JetLooseEle1_Pt / sT_eej                       , min_prescale * fakeRateEffective );
	 FillUserTH1D("sTfrac_Ele1_ROI"     , LooseEle1_Pt / sT_eej                          , min_prescale * fakeRateEffective );
	 FillUserTH1D("sTfrac_Ele2_ROI"     , LooseEle2_Pt / sT_eej                          , min_prescale * fakeRateEffective );
	 FillUserTH1D("sTfrac_Jet_ROI"      , ( JetLooseEle1_Pt  )                  / sT_eej , min_prescale * fakeRateEffective );
	 FillUserTH1D("sTfrac_Ele_ROI"      , ( LooseEle1_Pt + LooseEle2_Pt )       / sT_eej , min_prescale * fakeRateEffective );
       }
     } // End preselection 
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

