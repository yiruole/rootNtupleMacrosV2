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
   fillAllPreviousCuts              ( !true  ) ;
   fillAllOtherCuts                 ( !true  ) ;
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------

   CreateUserTH1D( "Ele1_Pt"	           , 	getHistoNBins("Ele1_Pt"), getHistoMin("Ele1_Pt"), getHistoMax("Ele1_Pt")     ) ; 
   CreateUserTH1D( "Ele1_Eta"	           , 	getHistoNBins("Ele1_Eta"), getHistoMin("Ele1_Eta"), getHistoMax("Ele1_Eta")     ) ; 
   CreateUserTH1D( "Ele1_Phi"	           , 	60  , -3.1416 , +3.1416	 ) ; 
   CreateUserTH1D( "Ele1_Charge"	   , 	2   , -1.0001 , 1.0001	 ) ; 

   CreateUserTH1D( "MET_Pt"                ,    getHistoNBins("MET_Pt"), getHistoMin("MET_Pt"), getHistoMax("MET_Pt")     ) ; 
   CreateUserTH1D( "MET_Phi"		   , 	60  , -3.1416 , +3.1416	 ) ; 

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

   CreateUserTH1D( "DR_Ele1Jet1"	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele1Jet2"	   , 	getHistoNBins("DR_Ele1Jet2"), getHistoMin("DR_Ele1Jet2"), getHistoMax("DR_Ele1Jet2")     ) ; 
   CreateUserTH1D( "DR_Jet1Jet2"	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
   CreateUserTH1D( "mDPhi_METEle1"	   , 	getHistoNBins("mDPhi_METEle1"), getHistoMin("mDPhi_METEle1"), getHistoMax("mDPhi_METEle1")     ) ; 
   CreateUserTH1D( "mDPhi_METJet1"	   , 	getHistoNBins("mDPhi_METJet1"), getHistoMin("mDPhi_METJet1"), getHistoMax("mDPhi_METJet1")     ) ; 
   CreateUserTH1D( "mDPhi_METJet2"	   , 	getHistoNBins("mDPhi_METJet2"), getHistoMin("mDPhi_METJet2"), getHistoMax("mDPhi_METJet2")     ) ; 
   CreateUserTH1D( "mDEta_Jet1Jet2"	   , 	getHistoNBins("mDEta_Jet1Jet2"), getHistoMin("mDEta_Jet1Jet2"), getHistoMax("mDEta_Jet1Jet2")   ) ; 

   CreateUserTH1D( "MT_Ele1MET"	           , 	getHistoNBins("MT_Ele1MET"), getHistoMin("MT_Ele1MET"), getHistoMax("MT_Ele1MET")     ) ; 
   CreateUserTH1D( "Pt_Ele1MET"	           , 	getHistoNBins("Pt_Ele1MET"), getHistoMin("Pt_Ele1MET"), getHistoMax("Pt_Ele1MET")     ) ; 
   CreateUserTH1D( "Pt_j1j2"	           , 	getHistoNBins("Pt_j1j2"), getHistoMin("Pt_j1j2"), getHistoMax("Pt_j1j2")     ) ; 
   CreateUserTH1D( "sT_enujj"	           , 	getHistoNBins("sT_enujj"), getHistoMin("sT_enujj"), getHistoMax("sT_enujj")     ) ; 
   CreateUserTH1D( "M_j1j2"	           , 	getHistoNBins("M_j1j2"), getHistoMin("M_j1j2"), getHistoMax("M_j1j2")     ) ; 

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
     
     double weight     = getPileupWeight ( nPileUpInteractions, isData ) ;

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "PassJSON"                      , passedJSON ) ; 

     // Noise filters
     fillVariableWithValue(   "PassHBHENoiseFilter"           , PassHBHENoiseFilter ) ; 
     fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight ) ; 

     // Electrons
     fillVariableWithValue(   "nEle"                          , nEle ) ;
     if ( nEle >= 1 ) { 
       fillVariableWithValue( "Ele1_Pt"                       , Ele1_Pt ) ;
       fillVariableWithValue( "Ele1_Eta"                      , Ele1_Eta ) ;
     }

     // MET variables
     fillVariableWithValue(   "MET_Pt"                        , MET_Pt ) ;
     if ( nEle >= 1 ) { 
       fillVariableWithValue( "MT_Ele1MET"                    , MT_Ele1MET ) ;
       fillVariableWithValue( "Pt_Ele1MET"                    , Pt_Ele1MET ) ;
     }
     
     // Jets
     fillVariableWithValue(   "nJet"                          , nJet ) ;
     if ( nJet >= 1 ) { 
       fillVariableWithValue( "Jet1_Pt"                       , Jet1_Pt ) ;
       fillVariableWithValue( "Jet1_Eta"                      , Jet1_Eta ) ;
     }
     if ( nJet >= 2 ) { 
       fillVariableWithValue( "Jet2_Pt"                       , Jet2_Pt ) ;
       fillVariableWithValue( "Jet2_Eta"                      , Jet2_Eta ) ;
       fillVariableWithValue( "Pt_j1j2"                       , Pt_j1j2 ) ;
       fillVariableWithValue( "M_j1j2"                        , M_j1j2 ) ;
       fillVariableWithValue( "mDEta_Jet1Jet2"                , fabs(Jet1_Eta - Jet2_Eta) ) ;
     }

     // Muons
     fillVariableWithValue(   "nMuon"                         , nMuon ) ;

     // DeltaR
     if ( nEle >= 1 && nJet >= 1) {
       fillVariableWithValue( "DR_Ele1Jet1"                   , DR_Ele1Jet1 ) ;
       if(nJet >= 2) {
	 fillVariableWithValue( "DR_Ele1Jet2"                 , DR_Ele1Jet2 ) ;
       }
     }

     // DeltaPhi
     if ( nEle >= 1 ) {      
       fillVariableWithValue( "mDPhi_METEle1"                 , mDPhi_METEle1 ) ;
     }
     if ( nJet >= 1 ) { 
       fillVariableWithValue( "mDPhi_METJet1"                 , mDPhi_METJet1 ) ;
     }
     if ( nJet >= 2 ) { 
       fillVariableWithValue( "mDPhi_METJet2"                 , mDPhi_METJet2 ) ;
     }

     // sT
     if ( nEle >= 1 && nJet >= 2) {
       fillVariableWithValue( "sT_enujj"                      , sT_enujj ) ;
     }      

     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------
     
     bool passed_preselection = passedAllPreviousCuts("Pt_Ele1MET");
     bool passed_preselection_without_DR_DPhi_cuts = passedAllPreviousCuts("nMuon");

     if ( passed_preselection ) { 

       FillUserTH1D( "Ele1_Pt"	           , 	Ele1_Pt       , weight);
       FillUserTH1D( "Ele1_Eta"	           , 	Ele1_Eta      , weight);
       FillUserTH1D( "Ele1_Phi"	           , 	Ele1_Phi      , weight);
       FillUserTH1D( "Ele1_Charge"	   , 	Ele1_Charge   , weight);
       
       FillUserTH1D( "MET_Pt"              ,    MET_Pt        , weight);
       FillUserTH1D( "MET_Phi"		   , 	MET_Phi       , weight);
       
       FillUserTH1D( "Jet1_Pt"             ,    Jet1_Pt       , weight);
       FillUserTH1D( "Jet1_Eta"	           , 	Jet1_Eta      , weight);
       FillUserTH1D( "Jet1_Phi"	           , 	Jet1_Phi      , weight);
       FillUserTH1D( "Jet1_btagTCHE"       ,    Jet1_btagTCHE , weight);
       FillUserTH1D( "Jet2_Pt"             ,    Jet2_Pt       , weight);
       FillUserTH1D( "Jet2_Eta"	           , 	Jet2_Eta      , weight);
       FillUserTH1D( "Jet2_Phi"	           , 	Jet2_Phi      , weight);
       FillUserTH1D( "Jet2_btagTCHE"       ,    Jet2_btagTCHE , weight);
       
       FillUserTH1D( "nEle"                ,    nEle          , weight);
       FillUserTH1D( "nJet"                ,    nJet          , weight);
       FillUserTH1D( "nJet_btagTCHE"       ,    nJet_btagTCHE , weight);
       FillUserTH1D( "nVertex"             ,    nVertex       , weight);
       FillUserTH1D( "nVertex_good"        ,    nVertex_good  , weight);
              
       FillUserTH1D( "MT_Ele1MET"	   , 	MT_Ele1MET    , weight);
       FillUserTH1D( "Pt_Ele1MET"	   , 	Pt_Ele1MET    , weight);
       FillUserTH1D( "Pt_j1j2"	           , 	Pt_j1j2       , weight);
       FillUserTH1D( "M_j1j2"	           , 	M_j1j2        , weight);
       FillUserTH1D( "sT_enujj"	           ,    sT_enujj      , weight);
       
     }

     if ( passed_preselection_without_DR_DPhi_cuts ) { 

       FillUserTH1D( "nMuon"               ,    nMuon         , weight);

       FillUserTH1D( "DR_Ele1Jet1"	   , 	DR_Ele1Jet1   , weight);
       FillUserTH1D( "DR_Ele1Jet2"	   , 	DR_Ele1Jet2   , weight);
       FillUserTH1D( "DR_Jet1Jet2"	   , 	DR_Jet1Jet2   , weight);
       FillUserTH1D( "mDPhi_METEle1"	   , 	mDPhi_METEle1 , weight);
       FillUserTH1D( "mDPhi_METJet1"	   , 	mDPhi_METJet1 , weight);
       FillUserTH1D( "mDPhi_METJet2"	   , 	mDPhi_METJet2 , weight);
       FillUserTH1D( "mDEta_Jet1Jet2"	   , 	getVariableValue("mDEta_Jet1Jet2")   , weight);

     }

   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
