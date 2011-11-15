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

   //wrt HT150
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_HT150"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_HT150_M400"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_HT150_HT350"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_HT150_M400ORHT350"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 

   //wrt HT200
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_HT200"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_HT200_M400"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_HT200_HT350"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_HT200_M400ORHT350"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 

   //wrt HT250
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_HT250"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_HT250_M400"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_HT250_HT350"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_HT250_M400ORHT350"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 


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
     fillVariableWithValue(   "HLT_HT150"                     , HLT_HT150 ) ; 
     fillVariableWithValue(   "HLT_HT200"                     , HLT_HT200 ) ; 
     fillVariableWithValue(   "HLT_HT250"                     , HLT_HT250 ) ; 
     fillVariableWithValue(   "HLT_HT350"                     , HLT_HT350 ) ; 
     fillVariableWithValue(   "HLT_M300"                      , HLT_FatJetMass300 ) ; 
     fillVariableWithValue(   "HLT_M400"                      , HLT_FatJetMass400 ) ; 

     //PFFatJets
     fillVariableWithValue(   "FatPFJet1_Pt"                  , FatPFJet1_Pt ) ; 
     fillVariableWithValue(   "FatPFJet1_Eta"                 , FatPFJet1_Eta ) ; 
     fillVariableWithValue(   "FatPFJet2_Pt"                  , FatPFJet2_Pt ) ; 
     fillVariableWithValue(   "FatPFJet2_Eta"                 , FatPFJet2_Eta ) ; 
     fillVariableWithValue(   "DEta_FatPFJet1FatPFJet2"       , DEta_FatPFJet1FatPFJet2 ) ; 
     fillVariableWithValue(   "M_FatPFJet1FatPFJet2"          , M_FatPFJet1FatPFJet2 ) ; 

     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------

     if( passedCut("PassJSON") && passedCut("PassPrimaryVertex") 
	 && variableIsFilled("M_FatPFJet1FatPFJet2") 
	 && variableIsFilled("DEta_FatPFJet1FatPFJet2") 
	 && variableIsFilled("FatPFJet1_Pt") 
	 && variableIsFilled("FatPFJet1_Eta") 
	 && variableIsFilled("FatPFJet2_Pt") 
	 && variableIsFilled("FatPFJet2_Eta") 
	 )
       {
	 
	 bool passed_dijet_cuts = false;
	 if( passedCut("DEta_FatPFJet1FatPFJet2") 
	     && passedCut("FatPFJet1_Pt") 
	     && passedCut("FatPFJet1_Eta")
	     && passedCut("FatPFJet2_Pt") 
	     && passedCut("FatPFJet2_Eta")
	     )
	   {	    
	     passed_dijet_cuts = true;
	   }

	 if( passedCut("HLT_HT150") && passed_dijet_cuts )
	   {
	     FillUserTH1D( "M_FatPFJet1FatPFJet2_HT150" , getVariableValue("M_FatPFJet1FatPFJet2") ) ; 
	     if( passedCut("HLT_M400") )
	       FillUserTH1D( "M_FatPFJet1FatPFJet2_HT150_M400" , getVariableValue("M_FatPFJet1FatPFJet2") ) ; 
	     if( passedCut("HLT_HT350") )
	       FillUserTH1D( "M_FatPFJet1FatPFJet2_HT150_HT350" , getVariableValue("M_FatPFJet1FatPFJet2") ) ; 
	     if( passedCut("HLT_M400") || passedCut("HLT_HT350") )
	       FillUserTH1D( "M_FatPFJet1FatPFJet2_HT150_M400ORHT350" , getVariableValue("M_FatPFJet1FatPFJet2") ) ; 
	   }    
	 
	 if( passedCut("HLT_HT200") && passed_dijet_cuts )
	   {
	     FillUserTH1D( "M_FatPFJet1FatPFJet2_HT200" , getVariableValue("M_FatPFJet1FatPFJet2") ) ; 
	     if( passedCut("HLT_M400") )
	       FillUserTH1D( "M_FatPFJet1FatPFJet2_HT200_M400" , getVariableValue("M_FatPFJet1FatPFJet2") ) ; 
	     if( passedCut("HLT_HT350") )
	       FillUserTH1D( "M_FatPFJet1FatPFJet2_HT200_HT350" , getVariableValue("M_FatPFJet1FatPFJet2") ) ; 
	     if( passedCut("HLT_M400") || passedCut("HLT_HT350") )
	       FillUserTH1D( "M_FatPFJet1FatPFJet2_HT200_M400ORHT350" , getVariableValue("M_FatPFJet1FatPFJet2") ) ; 
	   }    
	 
	 if( passedCut("HLT_HT250") && passed_dijet_cuts )
	   {
	     FillUserTH1D( "M_FatPFJet1FatPFJet2_HT250" , getVariableValue("M_FatPFJet1FatPFJet2") ) ; 
	     if( passedCut("HLT_M400") )
	       FillUserTH1D( "M_FatPFJet1FatPFJet2_HT250_M400" , getVariableValue("M_FatPFJet1FatPFJet2") ) ; 
	     if( passedCut("HLT_HT350") )
	       FillUserTH1D( "M_FatPFJet1FatPFJet2_HT250_HT350" , getVariableValue("M_FatPFJet1FatPFJet2") ) ; 
	     if( passedCut("HLT_M400") || passedCut("HLT_HT350") )
	       FillUserTH1D( "M_FatPFJet1FatPFJet2_HT250_M400ORHT350" , getVariableValue("M_FatPFJet1FatPFJet2") ) ; 
	   }    
	 	 
       } 

   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
