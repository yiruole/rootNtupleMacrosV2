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
    if (ientry < 0)
    {
      std::cout << "ERROR: Could not read from TTree; exiting." << std::endl;
      exit(-1);
    }
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (nb < 0)
    {
      std::cout << "ERROR: Could not read entry from TTree: read " << nb << "bytes; exiting." << std::endl;
      exit(-2);
    }

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
    // Do pileup re-weighting
    //--------------------------------------------------------------------------
    
    double pileup_weight = getPileupWeight ( nPileUpInt_True , isData ) ;
    
    //--------------------------------------------------------------------------
    // Get information about gen-level reweighting (should be for Sherpa only)
    //--------------------------------------------------------------------------

    double gen_weight = Weight;
    if ( isData ) gen_weight = 1.0;

    //--------------------------------------------------------------------------
    // Is this a barrel electron?
    //--------------------------------------------------------------------------

    bool Ele1_IsBarrel = bool ( fabs (Ele1_Eta) < 1.442 );

    //--------------------------------------------------------------------------
    // First variable to fill just shows the "reweighting".  Always passes.
    //--------------------------------------------------------------------------

    fillVariableWithValue ( "Reweighting", 1, gen_weight * pileup_weight  );

    //--------------------------------------------------------------------------
    // Fill JSON variable
    //--------------------------------------------------------------------------
    fillVariableWithValue ("PassJSON", passedJSON, gen_weight * pileup_weight  ); 

    //--------------------------------------------------------------------------
    // Fill noise filters
    //--------------------------------------------------------------------------
    // see: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    // we filled these at skim time
    fillVariableWithValue( "PassGlobalTightHalo2016Filter" , PassGlobalTightHalo2016Filter , gen_weight * pileup_weight  );
    fillVariableWithValue( "PassGoodVertices"	             , PassGoodVertices              , gen_weight * pileup_weight  );
    fillVariableWithValue( "PassHBHENoiseFilter"	         , PassHBHENoiseFilter           , gen_weight * pileup_weight  );
    fillVariableWithValue( "PassHBHENoiseIsoFilter"	       , PassHBHENoiseIsoFilter        , gen_weight * pileup_weight  );
    fillVariableWithValue( "PassBadEESupercrystalFilter"   , PassBadEESupercrystalFilter   , gen_weight * pileup_weight  );
    fillVariableWithValue( "PassEcalDeadCellTrigPrim"      , PassEcalDeadCellTrigPrim      , gen_weight * pileup_weight  );
    fillVariableWithValue( "PassChargedCandidateFilter"    , PassChargedCandidateFilter    , gen_weight * pileup_weight  );
    fillVariableWithValue( "PassBadPFMuonFilter"           , PassBadPFMuonFilter           , gen_weight * pileup_weight  );

    //------------------------------------------------------------------
    // Fill variables 
    //------------------------------------------------------------------
    //fillVariableWithValue("nMuon_ptCut"       , nMuon_ptCut             , gen_weight * pileup_weight );
    //fillVariableWithValue("nEle_ptCut"	      , nEle_ptCut	  , gen_weight * pileup_weight  );
    //fillVariableWithValue("Ele1_PtHeep"  	    , Ele1_PtHeep   , gen_weight * pileup_weight  );	  
    fillVariableWithValue("Ele1_Eta"          , Ele1_Eta                , gen_weight * pileup_weight );
    fillVariableWithValue("Ele1_IsBarrel"     , Ele1_IsBarrel           , gen_weight * pileup_weight );
    fillVariableWithValue("PFMET_Type1XY_Pt" , PFMET_Type1XY_Pt, gen_weight * pileup_weight   );	  	  
    //fillVariableWithValue("PFMET_Type1_Pt" , PFMET_Type1_Pt , gen_weight * pileup_weight  );	  	  
    fillVariableWithValue(   "mDeltaPhiMETEle"          , mDPhi_METEle1           , gen_weight * pileup_weight );
    fillVariableWithValue( "mDeltaPhiMET1stJet"       , mDPhi_METJet1           , gen_weight * pileup_weight );
    fillVariableWithValue(   "nJet"                     , nJet_ptCut              , gen_weight * pileup_weight );
    fillVariableWithValue("Jet1_Pt"   	      , Jet1_Pt   	  , gen_weight * pileup_weight  );
    fillVariableWithValue( "Jet1_Eta"                 , Jet1_Eta                , gen_weight * pileup_weight );
    fillVariableWithValue("Jet2_Pt"   	      , Jet2_Pt   	  , gen_weight * pileup_weight  );
    fillVariableWithValue("Jet2_Eta"                 , Jet2_Eta                , gen_weight * pileup_weight );
    fillVariableWithValue("sT_enujj"   	      , sT_enujj   	  , gen_weight * pileup_weight  );
    fillVariableWithValue ( "DR_Ele1Jet1"             , DR_Ele1Jet1             , gen_weight * pileup_weight );
    fillVariableWithValue ( "DR_Ele1Jet2"             , DR_Ele1Jet2             , gen_weight * pileup_weight );
    fillVariableWithValue("MT_Ele1MET"	      , MT_Ele1MET	  , gen_weight * pileup_weight  );	  
    fillVariableWithValue("PassFilter"        , 1             , gen_weight * pileup_weight      );
    
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
