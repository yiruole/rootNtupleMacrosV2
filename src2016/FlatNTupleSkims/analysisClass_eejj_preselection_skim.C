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
  
  fillAllPreviousCuts              ( !true  ) ;
  fillAllOtherCuts                 ( !true  ) ;
  fillAllSameLevelAndLowerLevelCuts( !true  ) ;
  fillAllCuts                      ( !true  ) ;
  
  //--------------------------------------------------------------------------
  // Analysis year
  //--------------------------------------------------------------------------
  int analysisYear = getPreCutValue1("AnalysisYear");

  //------------------------------------------------------------------
  // How many events to skim over?
  //------------------------------------------------------------------
  Long64_t nentries = GetTreeEntries();
  std::cout << "analysisClass::analysisClass(): nentries = " << nentries << std::endl;
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    readerTools_->LoadEntry(jentry);
    //--------------------------------------------------------------------------
    // Tricky part: refine Weight branch for MC only
    //--------------------------------------------------------------------------
    float weight = 1.0;

    if(!isData())
      resetSkimTreeBranchAddress("Weight", &weight);
    //------------------------------------------------------------------
    // Tell user how many events we've looped over
    //------------------------------------------------------------------
    if(jentry < 10 || jentry%10000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << "/" << nentries << std::endl;   
    
    //-----------------------------------------------------------------
    // Get ready to fill variables 
    //-----------------------------------------------------------------

    resetCuts();

    float gen_weight = readerTools_->ReadValueBranch<Float_t>("Weight");
    float pileup_weight = 1.0;
    // add trigger scale factor
    // have to modify for trigger path selection below
    // TODO: better way of handling this
    bool ele1PassedHLTWPTight = readerTools_->ReadValueBranch<Bool_t>("Ele1_PassedHLTriggerWPTightFilter");
    bool ele1PassedHLTCaloIdVTGsfTrkIdT = readerTools_->ReadValueBranch<Bool_t>("Ele1_PassedHLTriggerCaloIdVTGsfTrkIdTFilter");
    bool ele1PassedHLTPhoton = readerTools_->ReadValueBranch<Bool_t>("Ele1_PassedHLTriggerPhotonFilter");
    bool ele2PassedHLTWPTight = readerTools_->ReadValueBranch<Bool_t>("Ele2_PassedHLTriggerWPTightFilter");
    bool ele2PassedHLTCaloIdVTGsfTrkIdT = readerTools_->ReadValueBranch<Bool_t>("Ele2_PassedHLTriggerCaloIdVTGsfTrkIdTFilter");
    bool ele2PassedHLTPhoton = readerTools_->ReadValueBranch<Bool_t>("Ele2_PassedHLTriggerPhotonFilter");
    if(ele1PassedHLTWPTight || ele1PassedHLTPhoton) {
      float trigSFEle1 = readerTools_->ReadValueBranch<Float_t>("Ele1_TrigSF");
      float trigSFEle1Err = readerTools_->ReadValueBranch<Float_t>("Ele1_TrigSF_Err");
      gen_weight*=trigSFEle1;
    }
    else if(ele2PassedHLTWPTight || ele2PassedHLTPhoton) {
      float trigSFEle2 = readerTools_->ReadValueBranch<Float_t>("Ele2_TrigSF");
      float trigSFEle2Err = readerTools_->ReadValueBranch<Float_t>("Ele2_TrigSF_Err");
      gen_weight*=trigSFEle2;
    }
    ////float totalScaleFactor = recoHeepSF*trigSFEle1*trigSFEle2;
    //////std::cout << "trigSFEle1*trigSFEle2 = " << trigSFEle1 << "*" << trigSFEle2 << " = " << trigSFEle1*trigSFEle2 << std::endl;
    //////std::cout << "totalScaleFactor = " << totalScaleFactor << std::endl;
    weight = gen_weight;
    fillVariableWithValue ( "Reweighting", 1, gen_weight * pileup_weight  );

    //--------------------------------------------------------------------------
    // Fill HLT
    //--------------------------------------------------------------------------
    std::string current_file_name ( readerTools_->GetTree()->GetCurrentFile()->GetName());
    int passHLT = 0;
    bool verboseTrigEff = false;
    if ( isData() ) {
      if(current_file_name.find("SinglePhoton") != std::string::npos) {
        if(analysisYear==2016) {
          if (readerTools_->ReadValueBranch<Float_t>("H_Photon175") == 1) // take events triggered by Photon175 only plus those triggered by Photon175 AND Ele27/Ele115
            passHLT = 1;
        }
        else {
          if (readerTools_->ReadValueBranch<Float_t>("H_Photon200") == 1) // take events triggered by Photon200 only plus those triggered by Photon200 AND Ele35
            passHLT = 1;
        }
      }
      else if(current_file_name.find("SingleElectron") != std::string::npos) {
        if(analysisYear==2016) {
          if (readerTools_->ReadValueBranch<Float_t>("H_Photon175") != 1 && 
              //(readerTools_->ReadValueBranch<Float_t>("H_Ele27_WPTight") == 1 || readerTools_->ReadValueBranch<Float_t>("H_Ele115_CIdVT_GsfIdT") == 1) ) // take events triggered only by Ele27 OR Ele115
              //(readerTools_->ReadValueBranch<Float_t>("H_Ele115_CIdVT_GsfIdT") == 1) ) // take events triggered only by Ele115
            (readerTools_->ReadValueBranch<Float_t>("H_Ele27_WPTight") == 1) )
              passHLT = 1;
        }
        else {
          if (readerTools_->ReadValueBranch<Float_t>("H_Photon200") != 1 && 
              readerTools_->ReadValueBranch<Float_t>("H_Ele35_WPTight") == 1 ) // take events triggered only by Ele35
            passHLT = 1;
        }
      }
      else if(analysisYear==2018) {
        if (readerTools_->ReadValueBranch<Float_t>("H_Photon200") == 1 ||
            readerTools_->ReadValueBranch<Float_t>("H_Ele32_WPTight") == 1 ||
            readerTools_->ReadValueBranch<Float_t>("H_Ele115_CIdVT_GsfIdT") == 1) // take events triggered by Photon200 OR Ele32 OR Ele115
          passHLT = 1;
      }
      //if(current_file_name.find("SingleElectron") != std::string::npos) {
      //  if(analysisYear==2016) {
      //    if (readerTools_->ReadValueBranch<Float_t>("H_Photon175") != 1 &&
      //        (readerTools_->ReadValueBranch<Float_t>("H_Ele27_WPTight") == 1 || readerTools_->ReadValueBranch<Float_t>("H_Ele115_CIdVT_GsfIdT") == 1) ) // take events triggered only by Ele27 OR Ele115
      //        //(readerTools_->ReadValueBranch<Float_t>("H_Ele115_CIdVT_GsfIdT") == 1) ) // take events triggered only by Ele115
      //      //(readerTools_->ReadValueBranch<Float_t>("H_Ele27_WPTight") == 1) ) // take events triggered only by Ele27
      //        passHLT = 1;
      //  }
      //  else {
      //    if (readerTools_->ReadValueBranch<Float_t>("H_Photon200") != 1 && 
      //        readerTools_->ReadValueBranch<Float_t>("H_Ele35_WPTight") == 1 ) // take events triggered only by Ele35
      //      passHLT = 1;
      //  }
      //}
      //else if(analysisYear==2018) {
      //  if (readerTools_->ReadValueBranch<Float_t>("H_Photon200") == 1 ||
      //      readerTools_->ReadValueBranch<Float_t>("H_Ele32_WPTight") == 1 ||
      //      readerTools_->ReadValueBranch<Float_t>("H_Ele115_CIdVT_GsfIdT") == 1) // take events triggered by Photon200 OR Ele32 OR Ele115
      //    passHLT = 1;
      //}
    }
    else
    {
      //// a la Z', throw a random number and if it's below the efficiency at this pt/eta, pass the event
      ////   we get two chances to pass since we may have two electrons in the event
      //// trigger efficiency is binned in SCEta and SCEt (uncorrected)
      //// FIXME if used with stock nano...
      //// uncorrect the Pt, if there is a nonzero ecorr stored
      //float ele1ECorr = readerTools_->ReadValueBranch<Float_t>("Ele1_ECorr");
      //float ele2ECorr = readerTools_->ReadValueBranch<Float_t>("Ele2_ECorr");
      //float ele1PtUncorr = ele1ECorr != 0 ? readerTools_->ReadValueBranch<Float_t>("Ele1_Pt")/ele1ECorr : readerTools_->ReadValueBranch<Float_t>("Ele1_Pt");
      //float ele2PtUncorr = ele2ECorr != 0 ? readerTools_->ReadValueBranch<Float_t>("Ele2_Pt")/ele1ECorr : readerTools_->ReadValueBranch<Float_t>("Ele2_Pt");
      ////std::cout << "INFO: ele1PtUncorr = Ele1_Pt / Ele1_ECorr = " << readerTools_->ReadValueBranch<Float_t>("Ele1_Pt") << " / " << readerTools_->ReadValueBranch<Float_t>("Ele1_ECorr") << " = " << ele1PtUncorr << std::endl;
      //passHLT = triggerEfficiency.PassTrigger(readerTools_->ReadValueBranch<Float_t>("Ele1_SCEta"),ele1PtUncorr,verboseTrigEff) ? 1 : 0;
      //if(!passHLT) { // if the first one doesn't pass, try the second one
      //  //std::cout << "INFO: ele2PtUncorr = Ele2_Pt / Ele2_ECorr = " << readerTools_->ReadValueBranch<Float_t>("Ele2_Pt") << " / " << readerTools_->ReadValueBranch<Float_t>("Ele2_ECorr") << " = " << ele2PtUncorr << std::endl;
      //  passHLT = triggerEfficiency.PassTrigger(readerTools_->ReadValueBranch<Float_t>("Ele2_SCEta"),ele2PtUncorr,verboseTrigEff) ? 1 : 0;
      //}
      // using trigger scale factors
      if(analysisYear==2016) {
        if (//readerTools_->ReadValueBranch<Float_t>("H_Photon175") == 1 ||
            // XXX SIC: exclude Ele115 for now
            //readerTools_->ReadValueBranch<Float_t>("H_Ele115_CIdVT_GsfIdT") == 1 ||
            readerTools_->ReadValueBranch<Float_t>("H_Ele27_WPTight") == 1 )
          passHLT = 1;
      }
      else if(analysisYear==2017) {
        if (readerTools_->ReadValueBranch<Float_t>("H_Photon200") == 1 ||
            readerTools_->ReadValueBranch<Float_t>("H_Ele35_WPTight") == 1 )
          passHLT = 1;
      }
      else if(analysisYear==2018) {
        if (readerTools_->ReadValueBranch<Float_t>("H_Photon200") == 1 ||
            readerTools_->ReadValueBranch<Float_t>("H_Ele115_CIdVT_GsfIdT") == 1 ||
            readerTools_->ReadValueBranch<Float_t>("H_Ele32_WPTight") == 1 )
          passHLT = 1;
      }
    }
    fillVariableWithValue ( "PassHLT", passHLT, gen_weight * pileup_weight  ) ;     

    //--------------------------------------------------------------------------
    // Fill JSON variable
    //--------------------------------------------------------------------------
    fillVariableWithValue ("PassJSON", readerTools_->ReadValueBranch<Bool_t>("PassJSON"), gen_weight * pileup_weight); 

    //--------------------------------------------------------------------------
    // Fill noise filters
    //--------------------------------------------------------------------------
    // see: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    // we filled these at skim time
    fillVariableWithValue("PassGlobalSuperTightHalo2016Filter" , readerTools_->ReadValueBranch<Bool_t>("PassGlobalSuperTightHalo2016Filter")     , gen_weight * pileup_weight);
    fillVariableWithValue("PassGoodVertices"                   , readerTools_->ReadValueBranch<Bool_t>("PassGoodVertices")                       , gen_weight * pileup_weight);
    fillVariableWithValue("PassHBHENoiseFilter"                , readerTools_->ReadValueBranch<Bool_t>("PassHBHENoiseFilter")                    , gen_weight * pileup_weight);
    fillVariableWithValue("PassHBHENoiseIsoFilter"             , readerTools_->ReadValueBranch<Bool_t>("PassHBHENoiseIsoFilter")                 , gen_weight * pileup_weight);
    // eBadScFilter not suggested for MC
    if(isData())
      fillVariableWithValue("PassBadEESupercrystalFilter"      , readerTools_->ReadValueBranch<Bool_t>("PassBadEESupercrystalFilter")            , gen_weight * pileup_weight);
    else
      fillVariableWithValue("PassBadEESupercrystalFilter"      , 1                                                                                , gen_weight * pileup_weight);
    fillVariableWithValue("PassEcalDeadCellTrigPrim"           , readerTools_->ReadValueBranch<Bool_t>("PassEcalDeadCellTrigPrim")               , gen_weight * pileup_weight);
    // not recommended
    //fillVariableWithValue("PassChargedCandidateFilter"         , int(readerTools_->ReadValueB<Float_t>("PassChargedCandidateFilter")            == 1), gen_weight * pileup_weight);
    fillVariableWithValue("PassBadPFMuonFilter"                , readerTools_->ReadValueBranch<Bool_t>("PassBadPFMuonFilter")                    , gen_weight * pileup_weight);
    // EcalBadCalibV2 for 2017, 2018
    if(analysisYear > 2016)
      fillVariableWithValue("PassEcalBadCalibV2Filter"         , readerTools_->ReadValueBranch<Bool_t>("PassEcalBadCalibV2Filter")               , gen_weight * pileup_weight);
    else
      fillVariableWithValue("PassEcalBadCalibV2Filter"         , 1                                                                                , gen_weight * pileup_weight);

    //--------------------------------------------------------------------------
    // Calculate variables for trigger matching 
    //--------------------------------------------------------------------------

    int nEle_hltMatched = 0.0;
    //FIXME in reduced skim
    //if ( Ele1_hltEleSignalPt > 0.0 ) nEle_hltMatched++;
    //if ( Ele2_hltEleSignalPt > 0.0 ) nEle_hltMatched++;
    
    int nJet_hltMatched = 0.0;
    if ( readerTools_->ReadValueBranch<Float_t>("Jet1_hltJetPt") > 0.0 ) nJet_hltMatched++;
    if ( readerTools_->ReadValueBranch<Float_t>("Jet2_hltJetPt") > 0.0 ) nJet_hltMatched++;

    fillVariableWithValue("nEle_hltMatched",nEle_hltMatched, gen_weight * pileup_weight);
    fillVariableWithValue("nJet_hltMatched",nJet_hltMatched, gen_weight * pileup_weight);

    //--------------------------------------------------------------------------
    // Pass number of muons & electrons 
    // --> Special consideration if ttbar is derived from data
    //--------------------------------------------------------------------------

    //// Muons and electrons
    bool is_ttbar_from_data = false;
    //FIXME
    //if ( readerTools_->ReadValueBranch<Float_t>("Ele2_ValidFrac") > 998. ) is_ttbar_from_data = true;
    //
    int PassNEle = 0;
    //// nEle_ptCut are HEEP ID'ed electrons passing the Pt cut in the skim (which has been 10 GeV)
    //if ( !is_ttbar_from_data && nEle_ptCut == 2 ) PassNEle = 1;
    //if (  is_ttbar_from_data && nEle_ptCut == 2 ) PassNEle = 1;
    //if ( readerTools_->ReadValueBranch<Float_t>("nEle_ptCut") == 2 ) PassNEle = 1;
    if ( readerTools_->ReadValueBranch<Int_t>("nEle_store") >= 2 ) PassNEle = 1;
    fillVariableWithValue("PassNEle" , PassNEle , gen_weight * pileup_weight);
    int PassNJet = 0;
    if( readerTools_->ReadValueBranch<Int_t>("nJet_store") >= 2) PassNJet = 1;
    fillVariableWithValue("PassNJet" , PassNJet , gen_weight * pileup_weight);

    // remove muon veto for fake rate calc
    int PassNMuon = 0;
    //if ( !is_ttbar_from_data && nMuon_ptCut == 0 ) PassNMuon = 1;
    //if (  is_ttbar_from_data && nMuon_ptCut >  0 ) PassNMuon = 1;
    if (  readerTools_->ReadValueBranch<Int_t>("nMuon_ptCut") == 0 ) PassNMuon = 1;

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

    if ( readerTools_->ReadValueBranch<Int_t>("nEle_store") >= 2 &&
        readerTools_->ReadValueBranch<Int_t>("nJet_store") >= 2) {
      double M_e1j1 = readerTools_->ReadValueBranch<Float_t>("M_e1j1");
      double M_e1j2 = readerTools_->ReadValueBranch<Float_t>("M_e1j2");
      double M_e2j1 = readerTools_->ReadValueBranch<Float_t>("M_e2j1");
      double M_e2j2 = readerTools_->ReadValueBranch<Float_t>("M_e2j2");
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
    
    if ( readerTools_->ReadValueBranch<Int_t>("nEle_store") >= 1 ) fillVariableWithValue( "Ele1_Pt", readerTools_->ReadValueBranch<Float_t>("Ele1_Pt"), gen_weight * pileup_weight  ) ;
    if ( readerTools_->ReadValueBranch<Int_t>("nEle_store") >= 2 ) fillVariableWithValue( "Ele2_Pt", readerTools_->ReadValueBranch<Float_t>("Ele2_Pt"), gen_weight * pileup_weight  ) ;
		 
    //--------------------------------------------------------------------------
    // Fill jet variables 
    //--------------------------------------------------------------------------
		 		    
    fillVariableWithValue("nJet", readerTools_->ReadValueBranch<Int_t>("nJet_ptCut"), gen_weight * pileup_weight );
    if ( readerTools_->ReadValueBranch<Int_t>("nJet_store") >= 1 ) { 						    
      fillVariableWithValue( "Jet1_Pt"    , readerTools_->ReadValueBranch<Float_t>("Jet1_Pt")     , gen_weight * pileup_weight  ) ;
      fillVariableWithValue( "Jet1_Eta"   , readerTools_->ReadValueBranch<Float_t>("Jet1_Eta")    , gen_weight * pileup_weight  ) ;
    }
    if ( readerTools_->ReadValueBranch<Int_t>("nJet_store") >= 2 ) { 
      fillVariableWithValue( "Jet2_Pt"    , readerTools_->ReadValueBranch<Float_t>("Jet2_Pt")     , gen_weight * pileup_weight  ) ;
      fillVariableWithValue( "Jet2_Eta"   , readerTools_->ReadValueBranch<Float_t>("Jet2_Eta")    , gen_weight * pileup_weight  ) ;
      fillVariableWithValue( "DR_Jet1Jet2", readerTools_->ReadValueBranch<Float_t>("DR_Jet1Jet2") , gen_weight * pileup_weight  ) ;
    }

    //--------------------------------------------------------------------------
    // Fill DeltaR variables
    //--------------------------------------------------------------------------

    //std::cout << "SethLog: nEle_store = " << readerTools_->ReadValueBranch<Float_t>("nEle_store") << "; nJet_store = " << readerTools_->ReadValueBranch<Float_t>("nJet_store") << std::endl;
    if ( readerTools_->ReadValueBranch<Int_t>("nEle_store") >= 2 &&
        readerTools_->ReadValueBranch<Int_t>("nJet_store") >= 1) {
      //std::cout << "SethLog: Fill DR_Ele1Jet1 = " << readerTools_->ReadValueBranch<Float_t>("DR_Ele1Jet1") << std::endl;
      fillVariableWithValue( "DR_Ele1Jet1"  , readerTools_->ReadValueBranch<Float_t>("DR_Ele1Jet1") , gen_weight * pileup_weight  ) ;
      fillVariableWithValue( "DR_Ele2Jet1"  , readerTools_->ReadValueBranch<Float_t>("DR_Ele2Jet1") , gen_weight * pileup_weight  ) ;
      if(readerTools_->ReadValueBranch<Int_t>("nJet_store") >= 2) {
        fillVariableWithValue( "DR_Ele1Jet2", readerTools_->ReadValueBranch<Float_t>("DR_Ele1Jet2") , gen_weight * pileup_weight  ) ;
        fillVariableWithValue( "DR_Ele2Jet2", readerTools_->ReadValueBranch<Float_t>("DR_Ele2Jet2") , gen_weight * pileup_weight  ) ;
       }
     }


    //--------------------------------------------------------------------------
    // Multi-object variables
    //--------------------------------------------------------------------------

    if ( readerTools_->ReadValueBranch<Int_t>("nEle_store") >= 2 ) { 						    
      fillVariableWithValue( "M_e1e2"     , readerTools_->ReadValueBranch<Float_t>("M_e1e2") , gen_weight * pileup_weight  ) ;

      if ( readerTools_->ReadValueBranch<Int_t>("nJet_store") >= 2 ) { 
        fillVariableWithValue( "sT_eejj"    , readerTools_->ReadValueBranch<Float_t>("sT_eejj") , gen_weight * pileup_weight  ) ;
      }      
    }

    //-----------------------------------------------------------------
    // Evaluate the cuts
    //-----------------------------------------------------------------    

    evaluateCuts();
    
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
