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

    fillVariableWithValue("nLooseEle_ptCut", nLooseEle_ptCut );

    //------------------------------------------------------------------
    // Evaluate
    //------------------------------------------------------------------
    
    evaluateCuts();

    //------------------------------------------------------------------
    // If event passes, fill the tree
    //------------------------------------------------------------------

    if(passedCut("nLooseEle_ptCut")) fillSkimTree();
    
  } // End loop over events

  std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
