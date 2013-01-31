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
  fillAllPreviousCuts              (  true  ) ;
  fillAllOtherCuts                 (  true  ) ;
  fillAllSameLevelAndLowerLevelCuts( !true  ) ;
  fillAllCuts                      (  true  ) ;

  double lj_filter_lead_ele_min_pt = getPreCutValue1("LJFilterLeadEleMinPt");
  
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
    // What HLT paths to apply?
    //------------------------------------------------------------------
    
    char trigger_name[200];
    
    if ( isData ){ 
      if      ( run >= 190456 && run <= 190738 ) sprintf(trigger_name, "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v3");
      else if ( run >= 190782 && run <= 191419 ) sprintf(trigger_name, "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v4");
      else if ( run >= 191691 && run <= 196531 ) sprintf(trigger_name, "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v5");
      else if ( run >= 198022                  ) sprintf(trigger_name, "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v6");
    } 
    else                                         sprintf(trigger_name, "HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v5");
    
    //------------------------------------------------------------------
    // Get trigger decision
    //------------------------------------------------------------------
    
    getTriggers ( HLTKey, HLTInsideDatasetTriggerNames, HLTInsideDatasetTriggerDecisions,  HLTInsideDatasetTriggerPrescales ) ; 
    bool trigger_fired = triggerFired ( trigger_name );

    //------------------------------------------------------------------
    // Also skim on electron filter
    //------------------------------------------------------------------
    
    bool pass_lj_filter = false;
    if (  ElectronPt -> size() > 0 ) 
      pass_lj_filter = ((*ElectronPt)[0] > lj_filter_lead_ele_min_pt );

    //------------------------------------------------------------------
    // Fill cut values
    //------------------------------------------------------------------

    fillVariableWithValue ( "PassTrigger" , int ( trigger_fired  ) );
    fillVariableWithValue ( "PassLJFilter", int ( pass_lj_filter ) );

    //------------------------------------------------------------------
    // Evaluate cuts
    //------------------------------------------------------------------

    evaluateCuts();

    //------------------------------------------------------------------
    // If event passes, fill the tree
    //------------------------------------------------------------------

    if ( trigger_fired && pass_lj_filter) fillSkimTree();

   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
