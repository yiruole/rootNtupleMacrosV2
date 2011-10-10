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

   CreateUserTH1D( "Ele1_Pt"	           , 	getHistoNBins("Ele1_Pt"), getHistoMin("Ele1_Pt"), getHistoMax("Ele1_Pt")     ) ; 
   CreateUserTH1D( "Ele1_Eta"	           , 	getHistoNBins("Ele1_Eta"), getHistoMin("Ele1_Eta"), getHistoMax("Ele1_Eta")     ) ; 
   CreateUserTH1D( "Ele1_Phi"	           , 	60  , -3.1416 , +3.1416	 ) ; 
   CreateUserTH1D( "Ele1_Charge"     	   , 	2   , -1.0001 , 1.0001	 ) ; 
   CreateUserTH1D( "Ele2_Pt"	           , 	getHistoNBins("Ele2_Pt"), getHistoMin("Ele2_Pt"), getHistoMax("Ele2_Pt")     ) ; 
   CreateUserTH1D( "Ele2_Eta"	           , 	getHistoNBins("Ele2_Eta"), getHistoMin("Ele2_Eta"), getHistoMax("Ele2_Eta")     ) ; 
   CreateUserTH1D( "Ele2_Phi"	           , 	60  , -3.1416 , +3.1416	 ) ; 
   CreateUserTH1D( "Ele2_Charge"     	   , 	2   , -1.0001 , 1.0001	 ) ; 
		                     
   CreateUserTH1D( "MET_Pt"                ,    200, 0, 2000     ) ; 
   CreateUserTH1D( "MET_Phi"	     	   , 	60  , -3.1416 , +3.1416	 ) ; 
		                     
   CreateUserTH1D( "Jet1_Pt"               ,    getHistoNBins("Jet1_Pt"), getHistoMin("Jet1_Pt"), getHistoMax("Jet1_Pt")     ) ; 
   CreateUserTH1D( "Jet1_Eta"	           , 	getHistoNBins("Jet1_Eta"), getHistoMin("Jet1_Eta"), getHistoMax("Jet1_Eta")     ) ; 
   CreateUserTH1D( "Jet1_Phi"	           , 	60  , -3.1416 , +3.1416	 ) ;  
   CreateUserTH1D( "Jet1_btagTCHE"         ,    250 , 0       , 50	 ) ; 
   CreateUserTH1D( "Jet2_Pt"               ,    getHistoNBins("Jet2_Pt"), getHistoMin("Jet2_Pt"), getHistoMax("Jet2_Pt")     ) ; 
   CreateUserTH1D( "Jet2_Eta"	           , 	getHistoNBins("Jet2_Eta"), getHistoMin("Jet2_Eta"), getHistoMax("Jet2_Eta")     ) ; 
   CreateUserTH1D( "Jet2_Phi"	           , 	60  , -3.1416 , +3.1416	 ) ;  
   CreateUserTH1D( "Jet2_btagTCHE"         ,    250 , 0       , 50	 ) ; 
		                     
   CreateUserTH1D( "nEle"                  ,    getHistoNBins("nEle"), getHistoMin("nEle"), getHistoMax("nEle")     ) ; 
   CreateUserTH1D( "nMuon"                 ,    getHistoNBins("nMuon"), getHistoMin("nMuon"), getHistoMax("nMuon")     ) ; 
   CreateUserTH1D( "nJet"                  ,    getHistoNBins("nJet"), getHistoMin("nJet"), getHistoMax("nJet")     ) ; 
   CreateUserTH1D( "nJet_btagTCHE"         ,    getHistoNBins("nJet"), getHistoMin("nJet"), getHistoMax("nJet")     ) ; 
   CreateUserTH1D( "nVertex"               ,    31   , -0.5   , 30.5	 ) ; 
   CreateUserTH1D( "nVertex_good"          ,    31   , -0.5   , 30.5	 ) ; 
		                     
   CreateUserTH1D( "DR_Ele1Jet1"     	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele1Jet2"     	   , 	getHistoNBins("DR_Ele1Jet2"), getHistoMin("DR_Ele1Jet2"), getHistoMax("DR_Ele1Jet2")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet1"     	   , 	getHistoNBins("DR_Ele2Jet1"), getHistoMin("DR_Ele2Jet1"), getHistoMax("DR_Ele2Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet2"     	   , 	getHistoNBins("DR_Ele2Jet2"), getHistoMin("DR_Ele2Jet2"), getHistoMax("DR_Ele2Jet2")     ) ; 
   CreateUserTH1D( "DR_Jet1Jet2"     	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
		                     
   CreateUserTH1D( "M_e1e2"	           , 	getHistoNBins("M_e1e2"), getHistoMin("M_e1e2"), getHistoMax("M_e1e2")     ) ; 
   CreateUserTH1D( "sT_eejj"	           , 	getHistoNBins("sT_eejj"), getHistoMin("sT_eejj"), getHistoMax("sT_eejj")     ) ; 
   CreateUserTH1D( "Pt_e1e2"	           , 	getHistoNBins("Pt_e1e2"), getHistoMin("Pt_e1e2"), getHistoMax("Pt_e1e2")     ) ; 
   CreateUserTH1D( "Pt_j1j2"	           , 	getHistoNBins("Pt_j1j2"), getHistoMin("Pt_j1j2"), getHistoMax("Pt_j1j2")     ) ; 
   CreateUserTH1D( "M_j1j2"	           , 	getHistoNBins("M_j1j2"), getHistoMin("M_j1j2"), getHistoMax("M_j1j2")     ) ; 

   // CreateUserTH2D( "MejVsMej_selected",     200, 0, 2000, 200, 0, 2000) ;
   // CreateUserTH2D( "MejVsMej_rejected",     200, 0, 2000, 200, 0, 2000) ;

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
     
     int NPILEUP_AVE = int( (nPileUpInt_BXminus1 + nPileUpInt_BX0 + nPileUpInt_BXplus1)/3 );
     int NPILEUP_FINAL = min( NPILEUP_AVE , 25 );
     double weight = getPileupWeight ( NPILEUP_AVE, isData ) ;
     //double weight     = getPileupWeight ( nPileUpInteractions, isData ) ;

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "PassJSON"                      , passedJSON ) ; 

     // Noise filters
     fillVariableWithValue(   "PassHBHENoiseFilter"           , PassHBHENoiseFilter ) ; 
     fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight ) ; 

     // Fill HLT 
     int passHLT = 1;
     if ( isData ) { 
       passHLT = 0;
       if ( H_17_8_CIdT_CIsVL_TIdVL_TIsVL_5 == 1 || 
	    H_17_8_CIdT_CIsVL_TIdVL_TIsVL_6 == 1 || 
	    H_17_8_CIdT_CIsVL_TIdVL_TIsVL_7 == 1 || 
	    H_17_8_CIdT_CIsVL_TIdVL_TIsVL_8 == 1 || 
	    H_17_8_CIdT_TIdVL_CIsVL_TIsVL_1 == 1 || 
	    H_17_8_CIdT_TIdVL_CIsVL_TIsVL_2 == 1 || 
	    H_17_8_CIdT_TIdVL_CIsVL_TIsVL_3 == 1 || 
	    H_17_8_CIdT_TIdVL_CIsVL_TIsVL_4 == 1 || 
	    H_17_8_CIdT_TIdVL_CIsVL_TIsVL_5 == 1    ) {
	 passHLT = 1;
       }
     }
     
     fillVariableWithValue ( "PassHLT", passHLT ) ;

     // Electrons
     fillVariableWithValue(   "nEle"                          , nEle_Ana ) ;
     if ( nEle_Ana >= 1 ) { 
       fillVariableWithValue( "Ele1_Pt"                       , Ele1_Pt ) ;
       fillVariableWithValue( "Ele1_Eta"                      , Ele1_Eta ) ;
     }
     if ( nEle_Ana >= 2 ) { 
       fillVariableWithValue( "Ele2_Pt"                       , Ele2_Pt ) ;
       fillVariableWithValue( "Ele2_Eta"                      , Ele2_Eta ) ;
       fillVariableWithValue( "M_e1e2"                        , M_e1e2 ) ;
       fillVariableWithValue( "Pt_e1e2"                        , Pt_e1e2 ) ;
     }

     // Jets
     fillVariableWithValue(   "nJet"                          , nJet_Ana ) ;
     if ( nJet_Ana >= 1 ) { 
       fillVariableWithValue( "Jet1_Pt"                       , Jet1_Pt ) ;
       fillVariableWithValue( "Jet1_Eta"                      , Jet1_Eta ) ;
     }
     if ( nJet_Ana >= 2 ) { 
       fillVariableWithValue( "Jet2_Pt"                       , Jet2_Pt ) ;
       fillVariableWithValue( "Jet2_Eta"                      , Jet2_Eta ) ;
       fillVariableWithValue( "Pt_j1j2"                       , Pt_j1j2 ) ;
       fillVariableWithValue( "M_j1j2"                        , M_j1j2 ) ;
       fillVariableWithValue( "DR_Jet1Jet2"                        , DR_Jet1Jet2 ) ;
     }

     // Muons
     fillVariableWithValue(   "nMuon"                         , nMuon_Ana ) ;

     // DeltaR
     if ( nEle_Ana >= 2 && nJet_Ana >= 1) {
       fillVariableWithValue( "DR_Ele1Jet1"                   , DR_Ele1Jet1 ) ;
       fillVariableWithValue( "DR_Ele2Jet1"                   , DR_Ele2Jet1 ) ;
       if(nJet_Ana >= 2) {
	 fillVariableWithValue( "DR_Ele1Jet2"                 , DR_Ele1Jet2 ) ;
	 fillVariableWithValue( "DR_Ele2Jet2"                 , DR_Ele2Jet2 ) ;
       }
     }

     // sT
     if ( nEle_Ana >= 2 && nJet_Ana >= 2) {
       fillVariableWithValue( "sT_eejj"                      , sT_eejj ) ;
     }      

     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------
     
     bool passed_preselection = passedAllPreviousCuts("M_e1e2");
     bool passed_preselection_without_DR_cuts = passedAllPreviousCuts("DR_Ele1Jet1");

     if ( passed_preselection ) { 

       FillUserTH1D( "Ele1_Pt"	           , 	Ele1_Pt       , weight);
       FillUserTH1D( "Ele1_Eta"	           , 	Ele1_Eta      , weight);
       FillUserTH1D( "Ele1_Phi"	           , 	Ele1_Phi      , weight);
       FillUserTH1D( "Ele1_Charge"	   , 	Ele1_Charge   , weight);
       FillUserTH1D( "Ele2_Pt"	           , 	Ele2_Pt       , weight);
       FillUserTH1D( "Ele2_Eta"	           , 	Ele2_Eta      , weight);
       FillUserTH1D( "Ele2_Phi"	           , 	Ele2_Phi      , weight);
       FillUserTH1D( "Ele2_Charge"	   , 	Ele2_Charge   , weight);
       
       FillUserTH1D( "Jet1_Pt"             ,    Jet1_Pt       , weight);
       FillUserTH1D( "Jet1_Eta"	           , 	Jet1_Eta      , weight);
       FillUserTH1D( "Jet1_Phi"	           , 	Jet1_Phi      , weight);
       FillUserTH1D( "Jet1_btagTCHE"       ,    Jet1_btagTCHE , weight);
       FillUserTH1D( "Jet2_Pt"             ,    Jet2_Pt       , weight);
       FillUserTH1D( "Jet2_Eta"	           , 	Jet2_Eta      , weight);
       FillUserTH1D( "Jet2_Phi"	           , 	Jet2_Phi      , weight);
       FillUserTH1D( "Jet2_btagTCHE"       ,    Jet2_btagTCHE , weight);

       FillUserTH1D( "MET_Pt"	           , 	MET_Pt        , weight);
       FillUserTH1D( "MET_Phi"	           , 	MET_Phi       , weight);
       
       FillUserTH1D( "nEle"                ,    nEle_Ana      , weight);
       FillUserTH1D( "nJet"                ,    nJet_Ana      , weight);
       FillUserTH1D( "nJet_btagTCHE"       ,    nJet_btagTCHE_Ana , weight);
       FillUserTH1D( "nMuon"               ,    nMuon_Ana     , weight);
       FillUserTH1D( "nVertex"             ,    nVertex       , weight);
       FillUserTH1D( "nVertex_good"        ,    nVertex_good  , weight);
              
       FillUserTH1D( "Pt_e1e2"	           , 	Pt_e1e2       , weight);
       FillUserTH1D( "M_e1e2"	           , 	M_e1e2        , weight);
       FillUserTH1D( "sT_eejj"	           ,    sT_eejj       , weight);
       FillUserTH1D( "Pt_j1j2"	           , 	Pt_j1j2       , weight);
       FillUserTH1D( "M_j1j2"	           , 	M_j1j2        , weight);

       // Mej histos
       /*
       if ( fabs(M_e1j1-M_e2j2) < fabs(M_e1j2-M_e2j1) ) 
	 {
	   FillUserTH2D( "MejVsMej_selected", M_e1j1, M_e2j2, weight ) ;
	   FillUserTH2D( "MejVsMej_rejected", M_e1j2, M_e2j1, weight ) ;
	 }
       else
	 {
	   FillUserTH2D( "MejVsMej_selected", M_e1j2, M_e2j1, weight ) ;
	   FillUserTH2D( "MejVsMej_rejected", M_e1j1, M_e2j2, weight ) ;
	 }
       */
     }

     if ( passed_preselection_without_DR_cuts ) { 

       FillUserTH1D( "DR_Ele1Jet1"	   , 	DR_Ele1Jet1   , weight);
       FillUserTH1D( "DR_Ele1Jet2"	   , 	DR_Ele1Jet2   , weight);
       FillUserTH1D( "DR_Ele2Jet1"	   , 	DR_Ele2Jet1   , weight);
       FillUserTH1D( "DR_Ele2Jet2"	   , 	DR_Ele2Jet2   , weight);
       FillUserTH1D( "DR_Jet1Jet2"	   , 	DR_Jet1Jet2   , weight);

     }

   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

