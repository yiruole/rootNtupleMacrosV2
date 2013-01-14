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
  fillAllSameLevelAndLowerLevelCuts(  true  ) ;
  fillAllCuts                      (  true  ) ;
  
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
    
    //-----------------------------------------------------------------
    // Get ready to fill variables 
    //-----------------------------------------------------------------

    resetCuts();

    //------------------------------------------------------------------
    // Is the muon going to be Ele1 or Ele2?
    //------------------------------------------------------------------

    bool muon_is_Ele1 = false;
    bool muon_is_Ele2 = false;
    if   ( Muon1_Pt > 0.1 ) { 
      if   ( Muon1_Pt >  Ele1_Pt ) muon_is_Ele1 = true;
      else                         muon_is_Ele2 = true;
    }
    
    //------------------------------------------------------------------
    // If the muon is Ele2, then recalculate the Ele2 variables
    //------------------------------------------------------------------
    
    if ( muon_is_Ele2 ) {

      // Get object 4-vectors

      TLorentzVector muon, electron, jet1, jet2, met, e1e2, e2j1, e2j2;
      muon    .SetPtEtaPhiM ( Muon1_Pt, Muon1_Eta, Muon1_Phi, 0.0 );
      electron.SetPtEtaPhiM ( Ele1_Pt , Ele1_Eta , Ele1_Phi , 0.0 );
      jet1    .SetPtEtaPhiM ( Jet1_Pt , Jet1_Eta , Jet1_Phi , 0.0 );
      jet2    .SetPtEtaPhiM ( Jet2_Pt , Jet2_Eta , Jet2_Phi , 0.0 );
      met     .SetPtEtaPhiM ( PFMET_Type01XY_Pt, 0.0, PFMET_Type01XY_Phi, 0.0 );
      e1e2    = muon + electron;
      e2j1    = muon + jet1;
      e2j2    = muon + jet2;

      // These values can be re-calculated

      Ele2_Pt       = muon.Pt();
      Ele2_Eta      = muon.Eta();
      Ele2_Phi      = muon.Phi();
      Ele2_Charge   = Muon1_Charge;

      DR_Ele2Jet1   = muon.DeltaR ( jet1 );
      DR_Ele2Jet2   = muon.DeltaR ( jet2 );
      mDPhi_METEle2 = fabs ( muon.DeltaPhi ( met ));
      
      M_e1e2        = e1e2.M();
      M_e2j1        = e2j1.M();
      M_e2j2        = e2j2.M();
      sT_eejj       = electron.Pt() + muon.Pt() + jet1.Pt() + jet2.Pt();
      Pt_e1e2       = e1e2.Pt();

      if ( nEle_store < 2) nEle_store = 2;
      if ( nEle_ptCut < 2) nEle_ptCut = 2;

      // These values cannot be re-calculated with the information stored

      Ele2_MissingHits    = -999.;
      Ele2_DCotTheta      = -999.;
      Ele2_Dist           = -999.;
      Ele2_Energy         = -999.;
      Ele2_VtxD0          = -999.;
      Ele2_hltDoubleElePt = -999.;
      Ele2_hltEleSignalPt = -999.;
      Ele2_hltEleTTbarPt  = -999.;
      
    }

    else if ( muon_is_Ele1 ) { 
      
      // Get object 4-vectors

      TLorentzVector muon, electron, jet1, jet2, met, e1e2, e1j1, e1j2;
      muon    .SetPtEtaPhiM ( Muon1_Pt, Muon1_Eta, Muon1_Phi, 0.0 );
      electron.SetPtEtaPhiM ( Ele1_Pt , Ele1_Eta , Ele1_Phi , 0.0 );
      jet1    .SetPtEtaPhiM ( Jet1_Pt , Jet1_Eta , Jet1_Phi , 0.0 );
      jet2    .SetPtEtaPhiM ( Jet2_Pt , Jet2_Eta , Jet2_Phi , 0.0 );
      met     .SetPtEtaPhiM ( PFMET_Type01XY_Pt, 0.0, PFMET_Type01XY_Phi, 0.0 );
      e1e2    = muon + electron;
      e1j1    = muon + jet1;
      e1j2    = muon + jet2;

      // Demote the first electron (electron 1 -> electron 2)
      
      DR_Ele2Jet1          = DR_Ele1Jet1;         
      DR_Ele2Jet2	   = DR_Ele1Jet2;	   
      Ele2_Charge	   = Ele1_Charge;	   
      Ele2_DCotTheta	   = Ele1_DCotTheta;	   
      Ele2_Dist  	   = Ele1_Dist;	   
      Ele2_Energy	   = Ele1_Energy;	   
      Ele2_Eta		   = Ele1_Eta;		   
      Ele2_MissingHits	   = Ele1_MissingHits;	   
      Ele2_Phi		   = Ele1_Phi;		   
      Ele2_Pt		   = Ele1_Pt;		   
      Ele2_VtxD0	   = Ele1_VtxD0;	   
      Ele2_hltDoubleElePt  = Ele1_hltDoubleElePt; 
      Ele2_hltEleSignalPt  = Ele1_hltEleSignalPt; 
      Ele2_hltEleTTbarPt   = Ele1_hltEleTTbarPt;  
      mDPhi_METEle1	   = mDPhi_METEle1;	
      M_e2j1               = M_e1j1;
      M_e2j2               = M_e1j2;
      
      // These values can be re-calculated (muon 1 -> electron 1)

      Ele1_Pt       = muon.Pt();
      Ele1_Eta      = muon.Eta();
      Ele1_Phi      = muon.Phi();
      Ele1_Charge   = Muon1_Charge;
      MT_Ele1MET    = sqrt ( 2.0 * muon.Pt() * met.Pt() * ( 1.0 - cos ( met.DeltaPhi(muon))));

      DR_Ele1Jet1   = muon.DeltaR ( jet1 );
      DR_Ele1Jet2   = muon.DeltaR ( jet2 );
      mDPhi_METEle1 = fabs ( muon.DeltaPhi ( met ));
      
      M_e1e2        = e1e2.M();
      M_e1j1        = e1j1.M();
      M_e1j2        = e1j2.M();
      sT_eejj       = electron.Pt() + muon.Pt() + jet1.Pt() + jet2.Pt();
      Pt_e1e2       = e1e2.Pt();

      if ( nEle_store < 2) nEle_store = 2;
      if ( nEle_ptCut < 2) nEle_ptCut = 2;
      
      // These values cannot be re-calculated with the information stored

      Ele1_MissingHits    = -999.;
      Ele1_DCotTheta      = -999.;
      Ele1_Dist           = -999.;
      Ele1_Energy         = -999.;
      Ele1_ValidFrac      = -999.;
      Ele1_VtxD0          = -999.;
      Ele1_hltDoubleElePt = -999.;
      Ele1_hltEleSignalPt = -999.;
      Ele1_hltEleTTbarPt  = -999.;
      
    }

    // std::cout << "N(muon, pt cut) = " << nMuon_ptCut << std::endl;

    Ele2_ValidFrac = 999.0; // this is a tag

    //------------------------------------------------------------------
    // Fill variables 
    //------------------------------------------------------------------
    
    fillVariableWithValue("nEle_ptCut" , nEle_ptCut );
    fillVariableWithValue("nMuon_ptCut", nMuon_ptCut);
    fillVariableWithValue("Muon1_Pt"   , Muon1_Pt   );
    fillVariableWithValue("Muon2_Pt"   , Muon2_Pt   );
    fillVariableWithValue("Ele1_Pt"    , Ele1_Pt    );	  
    fillVariableWithValue("Ele2_Pt"    , Ele2_Pt    );	  	  
    fillVariableWithValue("Jet1_Pt"    , Jet1_Pt    );
    fillVariableWithValue("Jet2_Pt"    , Jet2_Pt    );
    fillVariableWithValue("sT_eejj"    , sT_eejj    );
    fillVariableWithValue("M_e1e2"     , M_e1e2	    );	  
    fillVariableWithValue("PassFilter" , 1          );

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
