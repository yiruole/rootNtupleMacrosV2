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
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;

  std::cout << "analysisClass::analysisClass(): ends " << std::endl;
}

analysisClass::~analysisClass()
{
}

void analysisClass::Loop() {
  
  //--------------------------------------------------------------------------
  // Decide which plots to save (default is to save everything)
  //--------------------------------------------------------------------------
  
  fillSkim                         (  true  ) ;
  fillAllPreviousCuts              ( !true  ) ;
  fillAllOtherCuts                 ( !true  ) ;
  fillAllSameLevelAndLowerLevelCuts( !true  ) ;
  fillAllCuts                      ( !true  ) ;
  
  //------------------------------------------------------------------
  // How many events to skim over?
  //------------------------------------------------------------------
  
  Long64_t nentries = fChain->GetEntries();
  std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

  //------------------------------------------------------------------
  // Loop over events
  //------------------------------------------------------------------

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    //------------------------------------------------------------------
    // ROOT loop preamble
    //------------------------------------------------------------------

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //------------------------------------------------------------------
    // Tell user how many events we've looped over
    //------------------------------------------------------------------

    if(jentry < 10 || jentry%5000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
    
    //-----------------------------------------------------------------
    // Get ready to fill variables 
    //-----------------------------------------------------------------

    resetCuts();
    
    //--------------------------------------------------------------------------
    // Check good run list
    //--------------------------------------------------------------------------
    
    int    passedJSON = passJSON ( run, ls , isData ) ;

    //--------------------------------------------------------------------------
    // Fill JSON variable
    //--------------------------------------------------------------------------
    fillVariableWithValue ("PassJSON", passedJSON); 

    //--------------------------------------------------------------------------
    // Fill noise filters
    //--------------------------------------------------------------------------
    // see: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    // we filled these at skim time
    fillVariableWithValue( "PassGlobalTightHalo2016Filter" , PassGlobalTightHalo2016Filter );
    fillVariableWithValue( "PassGoodVertices"	             , PassGoodVertices              );
    fillVariableWithValue( "PassHBHENoiseFilter"	         , PassHBHENoiseFilter           );
    fillVariableWithValue( "PassHBHENoiseIsoFilter"	       , PassHBHENoiseIsoFilter        );
    fillVariableWithValue( "PassBadEESupercrystalFilter"   , PassBadEESupercrystalFilter   );
    fillVariableWithValue( "PassEcalDeadCellTrigPrim"      , PassEcalDeadCellTrigPrim      );
    fillVariableWithValue( "PassChargedCandidateFilter"    , PassChargedCandidateFilter    );
    fillVariableWithValue( "PassBadPFMuonFilter"           , PassBadPFMuonFilter           );

    //------------------------------------------------------------------
    // Fill variables 
    //------------------------------------------------------------------
    fillVariableWithValue("nEle_ptCut"	      , nEle_ptCut	  );
    fillVariableWithValue("Ele1_PtHeep"  	    , Ele1_PtHeep   );	  
    fillVariableWithValue("PFMET_Type1XY_Pt" , PFMET_Type1XY_Pt );	  	  
    //fillVariableWithValue("PFMET_Type1_Pt" , PFMET_Type1_Pt );	  	  
    fillVariableWithValue("Jet1_Pt"   	      , Jet1_Pt   	  );
    fillVariableWithValue("Jet2_Pt"   	      , Jet2_Pt   	  );
    fillVariableWithValue("sT_enujj"   	      , sT_enujj   	  );
    fillVariableWithValue("MT_Ele1MET"	      , MT_Ele1MET	  );	  
    fillVariableWithValue("PassFilter"        , 1                 );
    
    //-----------------------------------------------------------------
    // Evaluate the cuts
    //-----------------------------------------------------------------    

    evaluateCuts();
    
    //------------------------------------------------------------------
    // If event passes, fill the tree
    //------------------------------------------------------------------

    if ( passedCut            ("PassFilter") &&
	 passedAllPreviousCuts("PassFilter") ){
      fillSkimTree();
    }

   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
