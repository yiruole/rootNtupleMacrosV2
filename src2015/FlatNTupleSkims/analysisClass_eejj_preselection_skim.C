#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include "src/Ele27WPLooseTrigTurnOn.C"

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
  Long64_t nentries = GetTreeEntries();
  std::cout << "analysisClass::analysisClass(): nentries = " << nentries << std::endl;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    readerTools_->LoadEntry(jentry);
    //------------------------------------------------------------------
    // Tell user how many events we've looped over
    //------------------------------------------------------------------
    if(jentry < 10 || jentry%5000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << "/" << nentries << std::endl;   
    
    //-----------------------------------------------------------------
    // Get ready to fill variables 
    //-----------------------------------------------------------------

    resetCuts();

    //--------------------------------------------------------------------------
    // Check good run list
    //--------------------------------------------------------------------------
    
    int passedJSON = passJSON ( readerTools_->ReadValueBranch<Double_t>("run"),
        readerTools_->ReadValueBranch<Double_t>("ls"),
        isData() ) ;

    //--------------------------------------------------------------------------
    // Do pileup re-weighting
    //--------------------------------------------------------------------------
    
    //double pileup_weight = readerTools_->ReadValueBranch<Double_t>("puWeight");

    //double gen_weight = readerTools_->ReadValueBranch<Double_t>("Weight");
    //if ( isData() ) gen_weight = 1.0;
    // don't do pileup/genWeighting here
    float pileup_weight = 1.0;
    float gen_weight = 1.0;

    //// Ele2_ValidFrac==999 --> ttbar-type sample
    //if ( isData && Ele2_ValidFrac > 998. && getPreCutString1("TrigCorrForSingleLeptonFinalState")=="true"){
    //  // efficiency of electron firing the trigger is stored in hltEleTTbarPt of the muon
    //  // weight the event by 2-eff.
    //  gen_weight *= Ele1_Energy < -998 ? 2-Ele1_hltEleTTbarPt : 2-Ele2_hltEleTTbarPt;
    //}

    //// For TopPt reweighting
    //if (TopPtWeight != -999)
    //  gen_weight *= TopPtWeight;

    //------------------------------------------------------------------
    // Fill variables 
    //------------------------------------------------------------------

    //--------------------------------------------------------------------------
    // First variable to fill just shows the "reweighting".  Always passes.
    //--------------------------------------------------------------------------

    fillVariableWithValue ( "Reweighting", 1, gen_weight * pileup_weight);
    //continue;

    //--------------------------------------------------------------------------
    // Fill JSON variable
    //--------------------------------------------------------------------------

    // JSON variable
    fillVariableWithValue ("PassJSON", passedJSON, gen_weight * pileup_weight ); 

    //--------------------------------------------------------------------------
    // Fill noise filters
    //--------------------------------------------------------------------------

    // Noise/MET filters
    // see: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    // we filled these at skim time
    fillVariableWithValue("PassGlobalTightHalo2016Filter" , int(readerTools_->ReadValueBranch<Double_t>("PassGlobalTightHalo2016Filter")          == 1));
    fillVariableWithValue("PassGoodVertices"              , int(readerTools_->ReadValueBranch<Double_t>("PassGoodVertices")                       == 1));
    fillVariableWithValue("PassHBHENoiseFilter"           , int(readerTools_->ReadValueBranch<Double_t>("PassHBHENoiseFilter")                    == 1));
    fillVariableWithValue("PassHBHENoiseIsoFilter"        , int(readerTools_->ReadValueBranch<Double_t>("PassHBHENoiseIsoFilter")                 == 1));
    fillVariableWithValue("PassBadEESupercrystalFilter"   , int(readerTools_->ReadValueBranch<Double_t>("PassBadEESupercrystalFilter")            == 1));
    fillVariableWithValue("PassEcalDeadCellTrigPrim"      , int(readerTools_->ReadValueBranch<Double_t>("PassEcalDeadCellTrigPrim")               == 1));
    fillVariableWithValue("PassChargedCandidateFilter"    , int(readerTools_->ReadValueBranch<Double_t>("PassChargedCandidateFilter")             == 1));
    fillVariableWithValue("PassBadPFMuonFilter"           , int(readerTools_->ReadValueBranch<Double_t>("PassBadPFMuonFilter")                    == 1));


    //--------------------------------------------------------------------------
    // Calculate variables for trigger matching 
    //--------------------------------------------------------------------------

    int nEle_hltMatched = 0.0;
    //FIXME in reduced skim
    //if ( Ele1_hltEleSignalPt > 0.0 ) nEle_hltMatched++;
    //if ( Ele2_hltEleSignalPt > 0.0 ) nEle_hltMatched++;
    
    int nJet_hltMatched = 0.0;
    if ( readerTools_->ReadValueBranch<Double_t>("Jet1_hltJetPt") > 0.0 ) nJet_hltMatched++;
    if ( readerTools_->ReadValueBranch<Double_t>("Jet2_hltJetPt") > 0.0 ) nJet_hltMatched++;

    fillVariableWithValue("nEle_hltMatched",nEle_hltMatched, gen_weight * pileup_weight);
    fillVariableWithValue("nJet_hltMatched",nJet_hltMatched, gen_weight * pileup_weight);

    //--------------------------------------------------------------------------
    // Pass number of muons & electrons 
    // --> Special consideration if ttbar is derived from data
    //--------------------------------------------------------------------------

    //// Muons and electrons
    bool is_ttbar_from_data = false;
    if ( readerTools_->ReadValueBranch<Double_t>("Ele2_ValidFrac") > 998. ) is_ttbar_from_data = true;
    //
    int PassNEle = 0;
    //// nEle_ptCut are HEEP ID'ed electrons passing the Pt cut in the skim (which has been 10 GeV)
    //if ( !is_ttbar_from_data && nEle_ptCut == 2 ) PassNEle = 1;
    //if (  is_ttbar_from_data && nEle_ptCut == 2 ) PassNEle = 1;
    if ( readerTools_->ReadValueBranch<Double_t>("nEle_ptCut") == 2 ) PassNEle = 1;
    // SIC test
    //int PassNEle = 1;

    // remove muon veto for fake rate calc
    int PassNMuon = 0;
    //if ( !is_ttbar_from_data && nMuon_ptCut == 0 ) PassNMuon = 1;
    //if (  is_ttbar_from_data && nMuon_ptCut >  0 ) PassNMuon = 1;
    if (  readerTools_->ReadValueBranch<Double_t>("nMuon_ptCut") == 0 ) PassNMuon = 1;

    //fillVariableWithValue("PassNEle" , PassNEle , gen_weight * pileup_weight);
    //fillVariableWithValue("PassNMuon", PassNMuon, gen_weight * pileup_weight);
    
    //if(is_ttbar_from_data)
    //{
    //  fillVariableWithValue("nEle_ptCut" , nEle_ptCut , gen_weight * pileup_weight);
    //  fillVariableWithValue("nMuon_ptCut", nMuon_ptCut, gen_weight * pileup_weight);
    //}

    //--------------------------------------------------------------------------
    // Calculate electron-jet pair mass values
    //--------------------------------------------------------------------------

    double M_ej_avg, M_ej_min, M_ej_max;

    if ( readerTools_->ReadValueBranch<Double_t>("nEle_store") >= 2 &&
        readerTools_->ReadValueBranch<Double_t>("nJet_store") >= 2) {
      double M_e1j1 = readerTools_->ReadValueBranch<Double_t>("M_e1j1");
      double M_e1j2 = readerTools_->ReadValueBranch<Double_t>("M_e1j2");
      double M_e2j1 = readerTools_->ReadValueBranch<Double_t>("M_e2j1");
      double M_e2j2 = readerTools_->ReadValueBranch<Double_t>("M_e2j2");
      if ( fabs(M_e1j1-M_e2j2) < fabs(M_e1j2-M_e2j1) )  {
        M_ej_avg = (M_e1j1 + M_e2j2) / 2.0;
        if    ( M_e1j1 < M_e2j2 ) { M_ej_min = M_e1j1; M_ej_max = M_e2j2; }
        else                      { M_ej_min = M_e2j2; M_ej_max = M_e1j1; }
      }
      else { 
        M_ej_avg = (M_e1j2 + M_e2j1) / 2.0;
        if    ( M_e1j2 < M_e2j1 ) { M_ej_min = M_e1j2; M_ej_max = M_e2j1; }
        else                      { M_ej_min = M_e2j1; M_ej_max = M_e1j2; }
      }
    }

    //--------------------------------------------------------------------------
    // Fill electron variables 
    //--------------------------------------------------------------------------
    
    //if ( nEle_store >= 1 ) fillVariableWithValue( "Ele1_PtHeep", Ele1_PtHeep, gen_weight * pileup_weight  ) ;
    //if ( nEle_store >= 2 ) fillVariableWithValue( "Ele2_PtHeep", Ele2_PtHeep, gen_weight * pileup_weight  ) ;
		 
    //--------------------------------------------------------------------------
    // Fill jet variables 
    //--------------------------------------------------------------------------
		 		    
    // Jets								    
    fillVariableWithValue("nJet", readerTools_->ReadValueBranch<Double_t>("nJet_ptCut"), gen_weight * pileup_weight );
    if ( readerTools_->ReadValueBranch<Double_t>("nJet_store") >= 1 ) { 						    
      fillVariableWithValue( "Jet1_Pt"    , readerTools_->ReadValueBranch<Double_t>("Jet1_Pt")     , gen_weight * pileup_weight  ) ;
      fillVariableWithValue( "Jet1_Eta"   , readerTools_->ReadValueBranch<Double_t>("Jet1_Eta")    , gen_weight * pileup_weight  ) ;
    }
    if ( readerTools_->ReadValueBranch<Double_t>("nJet_store") >= 2 ) { 
      fillVariableWithValue( "Jet2_Pt"    , readerTools_->ReadValueBranch<Double_t>("Jet2_Pt")     , gen_weight * pileup_weight  ) ;
      fillVariableWithValue( "Jet2_Eta"   , readerTools_->ReadValueBranch<Double_t>("Jet2_Eta")    , gen_weight * pileup_weight  ) ;
      fillVariableWithValue( "DR_Jet1Jet2", readerTools_->ReadValueBranch<Double_t>("DR_Jet1Jet2") , gen_weight * pileup_weight  ) ;
    }

    //--------------------------------------------------------------------------
    // Fill DeltaR variables
    //--------------------------------------------------------------------------

    if ( readerTools_->ReadValueBranch<Double_t>("nEle_store") >= 2 &&
        readerTools_->ReadValueBranch<Double_t>("nJet_store") >= 1) {
      fillVariableWithValue( "DR_Ele1Jet1"  , readerTools_->ReadValueBranch<Double_t>("DR_Ele1Jet1") , gen_weight * pileup_weight  ) ;
      fillVariableWithValue( "DR_Ele2Jet1"  , readerTools_->ReadValueBranch<Double_t>("DR_Ele2Jet1") , gen_weight * pileup_weight  ) ;
      if(readerTools_->ReadValueBranch<Double_t>("nJet_store") >= 2) {
        fillVariableWithValue( "DR_Ele1Jet2", readerTools_->ReadValueBranch<Double_t>("DR_Ele1Jet2") , gen_weight * pileup_weight  ) ;
        fillVariableWithValue( "DR_Ele2Jet2", readerTools_->ReadValueBranch<Double_t>("DR_Ele2Jet2") , gen_weight * pileup_weight  ) ;
       }
     }


    //--------------------------------------------------------------------------
    // Multi-object variables
    //--------------------------------------------------------------------------

    if ( readerTools_->ReadValueBranch<Double_t>("nEle_store") >= 2 ) { 						    
      fillVariableWithValue( "M_e1e2"     , readerTools_->ReadValueBranch<Double_t>("M_e1e2") , gen_weight * pileup_weight  ) ;
      fillVariableWithValue( "M_e1e2_opt" , readerTools_->ReadValueBranch<Double_t>("M_e1e2") , gen_weight * pileup_weight  ) ;

      if ( readerTools_->ReadValueBranch<Double_t>("nJet_store") >= 2 ) { 
        fillVariableWithValue( "sT_eejj"    , readerTools_->ReadValueBranch<Double_t>("sT_eejj") , gen_weight * pileup_weight  ) ;
        fillVariableWithValue( "sT_eejj_opt", readerTools_->ReadValueBranch<Double_t>("sT_eejj") , gen_weight * pileup_weight  ) ;
        fillVariableWithValue( "Mej_min_opt", M_ej_min, gen_weight * pileup_weight  ) ;
      }      
    }

    //-----------------------------------------------------------------
    // Evaluate the cuts
    //-----------------------------------------------------------------    

    evaluateCuts();
    
    //------------------------------------------------------------------
    // If event passes, fill the tree
    //------------------------------------------------------------------

    //if ( passedCut            ("PassFilter") &&
    //    passedAllPreviousCuts("PassFilter") ){
    //  fillSkimTree();
    //}
    // XXX TEST
    //fillSkimTree();
    if ( passedCut            ("M_e1e2") &&
        passedAllPreviousCuts("M_e1e2") ){
      fillSkimTree();
    }
    // XXX TEST

   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
