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
   fillAllPreviousCuts              ( true  ) ;
   fillAllOtherCuts                 ( true  ) ;
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------

//    CreateUserTH1D( "M_FatCaloCJet1FatCaloCJet2"	  , 
// 		   getHistoNBins("M_FatCaloCJet1FatCaloCJet2"), 
// 		   getHistoMin("M_FatCaloCJet1FatCaloCJet2"), 
// 		   getHistoMax("M_FatCaloCJet1FatCaloCJet2")     ) ; 

//    CreateUserTH1D( "M_FatPFJet1FatPFJet2"	  , 
// 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
// 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
// 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 

//    CreateUserTH1D( "DEta_FatPFJet1FatPFJet2"	  , 
// 		   getHistoNBins("DEta_FatPFJet1FatPFJet2"), 
// 		   getHistoMin("DEta_FatPFJet1FatPFJet2"), 
// 		   getHistoMax("DEta_FatPFJet1FatPFJet2")     ) ; 

//    CreateUserTH1D( "DPhi_FatPFJet1FatPFJet2"	  , 
// 		   getHistoNBins("DPhi_FatPFJet1FatPFJet2"), 
// 		   getHistoMin("DPhi_FatPFJet1FatPFJet2"), 
// 		   getHistoMax("DPhi_FatPFJet1FatPFJet2")     ) ; 

   //--------------------------------------------------------------------------
   // Precuts
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
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << "/" << nentries << std::endl;   

     //--------------------------------------------------------------------------
     // Reset the cuts
     //--------------------------------------------------------------------------

     resetCuts();

     //--------------------------------------------------------------------------
     // Check good run list
     //--------------------------------------------------------------------------
     
     int    passedJSON = passJSON ( run , ls , isData ) ;

     //--------------------------------------------------------------------------
     // Event-weight
     //--------------------------------------------------------------------------
     
     //     double event_weight = 1;

     //--------------------------------------------------------------------------
     // Do pileup re-weighting
     //--------------------------------------------------------------------------
     
     //      int NPILEUP_AVE = int( (nPileUpInt_BXminus1 + nPileUpInt_BX0 + nPileUpInt_BXplus1)/3 );
     //      int NPILEUP_FINAL = min( NPILEUP_AVE , 25 );
     //      event_weight = getPileupWeight ( NPILEUP_FINAL, isData ) ;
     //      //event_weight = getPileupWeight ( min(nPileUpInt_BX0,25), isData ) ;

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "PassJSON"                      , passedJSON ) ; 

     // Event filters
     fillVariableWithValue(   "PassPrimaryVertex"             , PassPrimaryVertex ) ; 

     //HLT     
     //fillVariableWithValue(   "HLT_PhysicsDST"                , (HLT_HT350 || HLT_FatJetMass300 || HLT_FatJetMass400) ) ; 
     fillVariableWithValue(   "HLT_HT350orM400"               , (HLT_HT350 || HLT_FatJetMass400) ) ; 
     fillVariableWithValue(   "HLT_M400"                      , (HLT_FatJetMass400) ) ; 

//      //CaloFatJets corrected
//      fillVariableWithValue(   "FatCaloCJet1_Pt"               , FatCaloCJet1_Pt ) ; 
//      fillVariableWithValue(   "FatCaloCJet1_Eta"              , FatCaloCJet1_Eta ) ; 
//      fillVariableWithValue(   "FatCaloCJet1_Phi"              , FatCaloCJet1_Phi ) ; 

//      fillVariableWithValue(   "FatCaloCJet2_Pt"               , FatCaloCJet2_Pt ) ; 
//      fillVariableWithValue(   "FatCaloCJet2_Eta"              , FatCaloCJet2_Eta ) ; 
//      fillVariableWithValue(   "FatCaloCJet2_Phi"              , FatCaloCJet2_Phi ) ; 

//      fillVariableWithValue(   "DR_FatCaloCJet1FatCaloCJet2"   , DR_FatCaloCJet1FatCaloCJet2 ) ; 
//      fillVariableWithValue(   "DPhi_FatCaloCJet1FatCaloCJet2" , DPhi_FatCaloCJet1FatCaloCJet2 ) ; 
//      fillVariableWithValue(   "DEta_FatCaloCJet1FatCaloCJet2" , DEta_FatCaloCJet1FatCaloCJet2 ) ; 
//      fillVariableWithValue(   "M_FatCaloCJet1FatCaloCJet2"    , M_FatCaloCJet1FatCaloCJet2 ) ; 

//      //CaloJets corrected
//      fillVariableWithValue(   "CaloCJet1_Pt"                  , CaloCJet1_Pt ) ; 
//      fillVariableWithValue(   "CaloCJet1_Eta"                 , CaloCJet1_Eta ) ; 
//      fillVariableWithValue(   "CaloCJet1_Phi"                 , CaloCJet1_Phi ) ; 
//      fillVariableWithValue(   "CaloCJet2_Pt"                  , CaloCJet2_Pt) ; 
//      fillVariableWithValue(   "CaloCJet2_Eta"                 , CaloCJet2_Eta ) ; 
//      fillVariableWithValue(   "CaloCJet2_Phi"                 , CaloCJet2_Phi ) ; 
//      fillVariableWithValue(   "DR_CaloCJet1CaloCJet2"         , DR_CaloCJet1CaloCJet2 ) ; 
//      fillVariableWithValue(   "DPhi_CaloCJet1CaloCJet2"       , DPhi_CaloCJet1CaloCJet2 ) ; 
//      fillVariableWithValue(   "DEta_CaloCJet1CaloCJet2"       , DEta_CaloCJet1CaloCJet2 ) ; 
//      fillVariableWithValue(   "M_CaloCJet1CaloCJet2"          , M_CaloCJet1CaloCJet2 ) ; 
//      fillVariableWithValue(   "HT_CaloCJets"                  , HT_CaloCJets ) ; 

     //PFFatJets
     fillVariableWithValue(   "FatPFJet1_Pt"    , FatPFJet1_Pt ) ; 
     fillVariableWithValue(   "FatPFJet1_Eta"   , FatPFJet1_Eta ) ; 
     fillVariableWithValue(   "FatPFJet1_Phi"   , FatPFJet1_Phi ) ; 

     fillVariableWithValue(   "FatPFJet2_Pt"    , FatPFJet2_Pt ) ; 
     fillVariableWithValue(   "FatPFJet2_Eta"   , FatPFJet2_Eta ) ; 
     fillVariableWithValue(   "FatPFJet2_Phi"   , FatPFJet2_Phi ) ; 

     fillVariableWithValue(   "DR_FatPFJet1FatPFJet2"   , DR_FatPFJet1FatPFJet2 ) ; 
     fillVariableWithValue(   "DPhi_FatPFJet1FatPFJet2"   , DPhi_FatPFJet1FatPFJet2 ) ; 
     fillVariableWithValue(   "DEta_FatPFJet1FatPFJet2"   , DEta_FatPFJet1FatPFJet2 ) ; 
     fillVariableWithValue(   "M_FatPFJet1FatPFJet2"   , M_FatPFJet1FatPFJet2 ) ; 
     fillVariableWithValue(   "M_PFJet1PFJet2"   , M_PFJet1PFJet2 ) ; 

     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------

//      if( passedCut("PassJSON") && passedCut("HLT_PhysicsDST") 
// 	 && passedCut("FatCaloCJet1_Pt") && passedCut("FatCaloCJet2_Pt") )
//        {
// 	 FillUserTH1D( "M_FatCaloCJet1FatCaloCJet2"	           , M_FatCaloCJet1FatCaloCJet2	 );	 
//        }
     
//      if( passedCut("PassJSON") && passedCut("HLT_HT350orM400") && passedCut("PassPrimaryVertex") 
// 	 && passedCut("FatPFJet1_Pt") && passedCut("FatPFJet2_Pt") )
//        {
// 	 FillUserTH1D( "M_FatPFJet1FatPFJet2"	           , M_FatPFJet1FatPFJet2	 );	 
// 	 FillUserTH1D( "DEta_FatPFJet1FatPFJet2"	   , DEta_FatPFJet1FatPFJet2	 );	 
// 	 FillUserTH1D( "DPhi_FatPFJet1FatPFJet2"	   , DPhi_FatPFJet1FatPFJet2	 );	 
//        }
     
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
