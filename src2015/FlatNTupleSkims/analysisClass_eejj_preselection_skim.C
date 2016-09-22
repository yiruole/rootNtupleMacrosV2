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

    if(jentry < 10 || jentry%5000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << "/" << nentries << std::endl;   
    
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
    
    double pileup_weight = getPileupWeight ( nPileUpInt_True, isData ) ;

    double gen_weight = Weight;
    if ( isData ) gen_weight = 1.0;
    // Ele2_ValidFrac==999 --> ttbar-type sample
    if ( isData && Ele2_ValidFrac > 998. && getPreCutString1("TrigCorrForSingleLeptonFinalState")=="true"){
      // efficiency of electron firing the trigger is stored in hltEleTTbarPt of the muon
      // weight the event by 2-eff.
      gen_weight *= Ele1_Energy < -998 ? 2-Ele1_hltEleTTbarPt : 2-Ele2_hltEleTTbarPt;
    }
    //XXX For TopPt reweighting
    if (TopPtWeight != -999)
      gen_weight *= TopPtWeight;

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
    fillVariableWithValue(   "PassHBHENoiseFilter"	          , PassHBHENoiseFilter              , gen_weight * pileup_weight );
    fillVariableWithValue(   "PassHBHENoiseIsoFilter"	      , PassHBHENoiseIsoFilter             , gen_weight * pileup_weight);
    fillVariableWithValue(   "PassCSCBeamHaloFilterTight"    , PassCSCBeamHaloFilterTight        , gen_weight * pileup_weight );
    fillVariableWithValue(   "PassBadEESupercrystalFilter"   , PassBadEESupercrystalFilter       , gen_weight * pileup_weight );
    fillVariableWithValue(   "PassMuonTrackFilter"           , PassMuonTrackFilter               , gen_weight * pileup_weight );
    fillVariableWithValue(   "PassBadResolutionTrackFilter"  , PassBadResolutionTrackFilter      , gen_weight * pileup_weight );
    // no longer in 2016
    //fillVariableWithValue(   "PassPhysDecl"		      , ( isData == 1 ) ? PassPhysDecl		      : 1, gen_weight * pileup_weight);
    //fillVariableWithValue(   "PassPrimaryVertex"	    , PassPrimaryVertex                        , gen_weight * pileup_weight );

    // no longer in 2016
    ////--------------------------------------------------------------------------
    //// Exclude runs with bad beam spot ?
    ////--------------------------------------------------------------------------
    //bool passBadBeamspot = true;
    //if(run==259626 ||
    //    run==259636 ||
    //    run==259637 ||
    //    run==259681 ||
    //    run==259682 ||
    //    run==259683 ||
    //    run==259685)
    //  //continue;
    //  passBadBeamspot = false;

    //fillVariableWithValue(   "PassBadBeamspot"	    , passBadBeamspot                          , gen_weight * pileup_weight );

    //--------------------------------------------------------------------------
    // Fill HLT
    //--------------------------------------------------------------------------

    //int passHLT = 1;
    //if ( isData ) { 
    //  passHLT = 0;
    //  if ( H_Ele30_PFJet100_25 == 1 || H_Ele30_PFNoPUJet100_25  == 1 ){
    ////if ( H_DoubleEle33_CIdL_GsfIdVL == 1 ) { 
    //  	 passHLT = 1;
    //  }
    //}
    // no longer in 2016
    //// ignore Run2015C stuff
    //if(isData)
    //{
    //  if(run >= 254227 && run <= 254914) // in Run2015C 25 ns, there is no un-eta-restricted WPLoose path
    //    // continue;
    //    // no, treat this case like failing the HLT
    //    fillVariableWithValue("PassHLT",0,gen_weight*pileup_weight);
    //}

    int passHLT = 1;
    if ( isData ) { 
      passHLT = 0;
      if ( H_Ele27_WPTight == 1 || H_Photon175 == 1)
        passHLT = 1;
    }
    // FIXME: Update to new 2016 curve?
    else // using the turn-on in the MC?
    {
      // a la Z', throw a random number and if it's below the efficiency at this pt/eta, pass the event
      //   we get two chances to pass since we may have two electrons in the event
      passHLT = trigEle27::passTrig(Ele1_PtHeep,Ele1_SCEta) ? 1 : 0;
      //passHLT = trigEle27::passTrig(Ele1_Pt,Ele1_Eta) ? 1 : 0;
      if(!passHLT) // if the first one doesn't pass, try the second one
        passHLT = trigEle27::passTrig(Ele2_PtHeep,Ele2_SCEta) ? 1 : 0;
      //passHLT = trigEle27::passTrig(Ele2_Pt,Ele2_Eta) ? 1 : 0;
    }
    fillVariableWithValue ( "PassHLT", passHLT, gen_weight * pileup_weight  ) ;     

    //--------------------------------------------------------------------------
    // Calculate variables for trigger matching 
    //--------------------------------------------------------------------------

    int nEle_hltMatched = 0.0;
    if ( Ele1_hltEleSignalPt > 0.0 ) nEle_hltMatched++;
    if ( Ele2_hltEleSignalPt > 0.0 ) nEle_hltMatched++;
    
    int nJet_hltMatched = 0.0;
    if ( Jet1_hltNoPUJetPt > 0.0 || Jet1_hltJetPt > 0.0 ) nJet_hltMatched++;
    if ( Jet2_hltNoPUJetPt > 0.0 || Jet2_hltJetPt > 0.0 ) nJet_hltMatched++;

    fillVariableWithValue("nEle_hltMatched",nEle_hltMatched, gen_weight * pileup_weight);
    fillVariableWithValue("nJet_hltMatched",nJet_hltMatched, gen_weight * pileup_weight);

    //--------------------------------------------------------------------------
    // Pass number of muons & electrons 
    // --> Special consideration if ttbar is derived from data
    //--------------------------------------------------------------------------

    //// Muons and electrons
    bool is_ttbar_from_data = false;
    if ( Ele2_ValidFrac > 998. ) is_ttbar_from_data = true;
    //
    int PassNEle = 0;
    //// nEle_ptCut are HEEP ID'ed electrons passing the Pt cut in the skim (which has been 10 GeV)
    //if ( !is_ttbar_from_data && nEle_ptCut == 2 ) PassNEle = 1;
    //if (  is_ttbar_from_data && nEle_ptCut == 2 ) PassNEle = 1;
    if ( nEle_ptCut == 2 ) PassNEle = 1;
    // SIC test
    //int PassNEle = 1;

    int PassNMuon = 0;
    //if ( !is_ttbar_from_data && nMuon_ptCut == 0 ) PassNMuon = 1;
    //if (  is_ttbar_from_data && nMuon_ptCut >  0 ) PassNMuon = 1;
    if (  nMuon_ptCut == 0 ) PassNMuon = 1;

    fillVariableWithValue("PassNEle" , PassNEle , gen_weight * pileup_weight);
    fillVariableWithValue("PassNMuon", PassNMuon, gen_weight * pileup_weight);

    if(is_ttbar_from_data)
    {
      fillVariableWithValue("nEle_ptCut" , nEle_ptCut , gen_weight * pileup_weight);
      fillVariableWithValue("nMuon_ptCut", nMuon_ptCut, gen_weight * pileup_weight);
    }

    //--------------------------------------------------------------------------
    // Calculate electron-jet pair mass values
    //--------------------------------------------------------------------------

    double M_ej_avg, M_ej_min, M_ej_max;

    if ( nEle_store >= 2 && nJet_store >= 2) {
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
    
    if ( nEle_store >= 1 ) fillVariableWithValue( "Ele1_PtHeep", Ele1_PtHeep, gen_weight * pileup_weight  ) ;
    if ( nEle_store >= 2 ) fillVariableWithValue( "Ele2_PtHeep", Ele2_PtHeep, gen_weight * pileup_weight  ) ;
		 
    //--------------------------------------------------------------------------
    // Fill jet variables 
    //--------------------------------------------------------------------------
		 		    
    // Jets								    
    fillVariableWithValue("nJet", nJet_ptCut, gen_weight * pileup_weight );
    if ( nJet_store >= 1 ) { 						    
      fillVariableWithValue( "Jet1_Pt"    , Jet1_Pt     , gen_weight * pileup_weight  ) ;
      fillVariableWithValue( "Jet1_Eta"   , Jet1_Eta    , gen_weight * pileup_weight  ) ;
    }
    if ( nJet_store >= 2 ) { 
      fillVariableWithValue( "Jet2_Pt"    , Jet2_Pt     , gen_weight * pileup_weight  ) ;
      fillVariableWithValue( "Jet2_Eta"   , Jet2_Eta    , gen_weight * pileup_weight  ) ;
      fillVariableWithValue( "DR_Jet1Jet2", DR_Jet1Jet2 , gen_weight * pileup_weight  ) ;
    }

    //--------------------------------------------------------------------------
    // Fill DeltaR variables
    //--------------------------------------------------------------------------

    if ( nEle_store >= 2 && nJet_store >= 1) {
      fillVariableWithValue( "DR_Ele1Jet1"  , DR_Ele1Jet1 , gen_weight * pileup_weight  ) ;
      fillVariableWithValue( "DR_Ele2Jet1"  , DR_Ele2Jet1 , gen_weight * pileup_weight  ) ;
      if(nJet_store >= 2) {
        fillVariableWithValue( "DR_Ele1Jet2", DR_Ele1Jet2 , gen_weight * pileup_weight  ) ;
        fillVariableWithValue( "DR_Ele2Jet2", DR_Ele2Jet2 , gen_weight * pileup_weight  ) ;
       }
     }


    //--------------------------------------------------------------------------
    // Multi-object variables
    //--------------------------------------------------------------------------

    if ( nEle_store >= 2 ) { 						    
      fillVariableWithValue( "M_e1e2"     , M_e1e2 , gen_weight * pileup_weight  ) ;
      fillVariableWithValue( "M_e1e2_opt" , M_e1e2 , gen_weight * pileup_weight  ) ;

      if ( nJet_store >= 2 ) { 
        // SIC recompute sT using PtHeep. FIXME: this is now being done in skims
        //sT_eejj = Ele1_PtHeep+Ele2_PtHeep+Jet1_Pt+Jet2_Pt;
        fillVariableWithValue( "sT_eejj"    , sT_eejj , gen_weight * pileup_weight  ) ;
        fillVariableWithValue( "sT_eejj_opt", sT_eejj , gen_weight * pileup_weight  ) ;
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
    if ( passedCut            ("M_e1e2") &&
        passedAllPreviousCuts("M_e1e2") ){
      fillSkimTree();
    }

   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
