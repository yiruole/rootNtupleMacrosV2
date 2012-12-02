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
  std::cout << "analysisClass::~analysisClass(): begins " << std::endl;
  std::cout << "analysisClass::~analysisClass(): ends " << std::endl;
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

    if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
    
    //------------------------------------------------------------------
    // Fill variables
    //------------------------------------------------------------------

    fillVariableWithValue( "nLooseEle_ptCut"      , nLooseEle_ptCut       );
    fillVariableWithValue( "LooseEle1_hltPhotonPt", LooseEle1_hltPhotonPt );
    fillVariableWithValue( "LooseEle1_Pt"         , LooseEle1_Pt          );
    fillVariableWithValue( "LooseEle1_PassID"     , LooseEle1_PassID      );
    fillVariableWithValue( "JetLooseEle1_Pt"      , JetLooseEle1_Pt       );
    fillVariableWithValue( "JetLooseEle2_Pt"      , JetLooseEle2_Pt       );
    fillVariableWithValue( "MET_Pt"               , PFMET_Type01XY_Pt     );
    fillVariableWithValue( "sT_enujj"             , sT_enujj 	          );
    fillVariableWithValue( "MT_Ele1MET"           , M_e1e2                );
    
    //------------------------------------------------------------------
    // Evaluate
    //------------------------------------------------------------------
    
    evaluateCuts();

    //------------------------------------------------------------------
    // If event passes, fill the tree
    //------------------------------------------------------------------

    if ( passedCut("nLooseEle_ptCut"       ) && 
	 passedCut("LooseEle1_hltPhotonPt" ) && 
	 passedCut("LooseEle1_Pt"          ) && 
	 passedCut("LooseEle1_PassID"      ) && 
	 passedCut("JetLooseEle1_Pt"       ) && 
	 passedCut("JetLooseEle2_Pt"       ) && 
	 passedCut("MET_Pt"                ) && 
	 passedCut("sT_enujj"              ) && 
    	 passedCut("MT_Ele1MET"            ) ){
      fillSkimTree();
    }
	  
  } // End loop over events

  std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
