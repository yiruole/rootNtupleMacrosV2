#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TRandom3.h>
#include <map>
#include <stdio.h>

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  //STDOUT("analysisClass::analysisClass() was called");
}

analysisClass::~analysisClass()
{
  //STDOUT("analysisClass::~analysisClass() was called");
}

void analysisClass::Loop()
{
  //STDOUT("analysisClass::Loop() begins");

  if (fChain == 0) return;

  /*//------------------------------------------------------------------
   *
   *
   *      
   *      Get all Pre-cut values!
   *
   *
   *
   *///-----------------------------------------------------------------
  
   //--------------------------------------------------------------------------
   // Decide which plots to save (default is to save everything)
   //--------------------------------------------------------------------------
   
   fillSkim                         ( !true  ) ;
   fillAllPreviousCuts              ( !true  ) ;
   fillAllOtherCuts                 ( !true  ) ;
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;

  //-----------------------------------------------------------------
  // Electron cut values
  //-----------------------------------------------------------------
  
  double ele_PtCut_STORE = getPreCutValue2("ele_PtCut");
  double ele_PtCut_ANA   = getPreCutValue1("ele_PtCut");
  
  double eleEta_bar      = getPreCutValue1("eleEta_bar");
  double eleEta_end_min  = getPreCutValue1("eleEta_end");
  double eleEta_end_max  = getPreCutValue2("eleEta_end");
  
  if ( ele_PtCut_STORE > ele_PtCut_ANA ) {
    STDOUT("ERROR in Electron cut values: all storage cuts must be looser or equal to analysis cuts.");
    exit(0) ;
  }

  // For WP80

  double eleMissingHitsWP             = getPreCutValue1("eleMissingHitsWP"        );
  double eleDistWP                    = getPreCutValue1("eleDistWP"               );
  double eleDCotThetaWP               = getPreCutValue1("eleDCotThetaWP"          );
  double eleCombRelIsoWP_bar          = getPreCutValue1("eleCombRelIsoWP"         );
  double eleCombRelIsoWP_end          = getPreCutValue2("eleCombRelIsoWP"         );
  double eleSigmaIetaIetaWP_bar       = getPreCutValue1("eleSigmaIetaIetaWP"      );
  double eleSigmaIetaIetaWP_end       = getPreCutValue2("eleSigmaIetaIetaWP"      );
  double eleDeltaPhiTrkSCWP_bar       = getPreCutValue1("eleDeltaPhiTrkSCWP"      );
  double eleDeltaPhiTrkSCWP_end       = getPreCutValue2("eleDeltaPhiTrkSCWP"      );
  double eleDeltaEtaTrkSCWP_bar       = getPreCutValue1("eleDeltaEtaTrkSCWP"      );
  double eleDeltaEtaTrkSCWP_end       = getPreCutValue2("eleDeltaEtaTrkSCWP"      );
  double eleUseEcalDrivenWP           = getPreCutValue1("eleUseEcalDrivenWP"      );
  double eleUseHasMatchConvWP         = getPreCutValue1("eleUseHasMatchConvWP"    );

  // For HEEP 3.1

  double eleDeltaEtaTrkSCHeep_bar     = getPreCutValue1("eleDeltaEtaTrkSCHeep"    );
  double eleDeltaEtaTrkSCHeep_end     = getPreCutValue2("eleDeltaEtaTrkSCHeep"    );
  double eleDeltaPhiTrkSCHeep_bar     = getPreCutValue1("eleDeltaPhiTrkSCHeep"    );
  double eleDeltaPhiTrkSCHeep_end     = getPreCutValue2("eleDeltaPhiTrkSCHeep"    );
  double eleHoEHeep_bar               = getPreCutValue1("eleHoEHeep"              );
  double eleHoEHeep_end               = getPreCutValue2("eleHoEHeep"              );
  double eleE2x5OverE5x5Heep_bar      = getPreCutValue1("eleE2x5OverE5x5Heep"     );
  double eleE1x5OverE5x5Heep_bar      = getPreCutValue1("eleE1x5OverE5x5Heep"     );
  double eleSigmaIetaIetaHeep_end     = getPreCutValue2("eleSigmaIetaIetaHeep"    );
  double eleEcalHcalIsoHeep_1_bar     = getPreCutValue1("eleEcalHcalIsoHeep"      );
  double eleEcalHcalIsoHeep_2_bar     = getPreCutValue2("eleEcalHcalIsoHeep"      );
  double eleEcalHcalIsoHeep_1_end     = getPreCutValue3("eleEcalHcalIsoHeep"      );
  double eleEcalHcalIsoHeep_2_end     = getPreCutValue4("eleEcalHcalIsoHeep"      );
  double eleEcalHcalIsoHeep_PTthr_end = getPreCutValue2("eleEcalHcalIsoHeep_PTthr");
  double eleHcalIsoD2Heep_end         = getPreCutValue2("eleHcalIsoD2Heep"        );
  double eleTrkIsoHeep_bar            = getPreCutValue1("eleTrkIsoHeep"           );
  double eleTrkIsoHeep_end            = getPreCutValue2("eleTrkIsoHeep"           );
  double eleMissingHitsHeep           = getPreCutValue1("eleMissingHitsHeep"      );
  double eleUseEcalDrivenHeep         = getPreCutValue1("eleUseEcalDrivenHeep"    );

  // For HEEP 3.2

  double eleDeltaEtaTrkSCHeep32_bar     = getPreCutValue1("eleDeltaEtaTrkSCHeep32"    );
  double eleDeltaEtaTrkSCHeep32_end     = getPreCutValue2("eleDeltaEtaTrkSCHeep32"    );
  double eleDeltaPhiTrkSCHeep32_bar     = getPreCutValue1("eleDeltaPhiTrkSCHeep32"    );
  double eleDeltaPhiTrkSCHeep32_end     = getPreCutValue2("eleDeltaPhiTrkSCHeep32"    );
  double eleHoEHeep32_bar               = getPreCutValue1("eleHoEHeep32"              );
  double eleHoEHeep32_end               = getPreCutValue2("eleHoEHeep32"              );
  double eleE2x5OverE5x5Heep32_bar      = getPreCutValue1("eleE2x5OverE5x5Heep32"     );
  double eleE1x5OverE5x5Heep32_bar      = getPreCutValue1("eleE1x5OverE5x5Heep32"     );
  double eleSigmaIetaIetaHeep32_end     = getPreCutValue2("eleSigmaIetaIetaHeep32"    );
  double eleEcalHcalIsoHeep32_1_bar     = getPreCutValue1("eleEcalHcalIsoHeep32"      );
  double eleEcalHcalIsoHeep32_2_bar     = getPreCutValue2("eleEcalHcalIsoHeep32"      );
  double eleEcalHcalIsoHeep32_1_end     = getPreCutValue3("eleEcalHcalIsoHeep32"      );
  double eleEcalHcalIsoHeep32_2_end     = getPreCutValue4("eleEcalHcalIsoHeep32"      );
  double eleEcalHcalIsoHeep32_PTthr_end = getPreCutValue2("eleEcalHcalIsoHeep32_PTthr");
  //double eleHcalIsoD2Heep32_end         = getPreCutValue2("eleHcalIsoD2Heep32"        );
  double eleTrkIsoHeep32_bar            = getPreCutValue1("eleTrkIsoHeep32"           );
  double eleTrkIsoHeep32_end            = getPreCutValue2("eleTrkIsoHeep32"           );
  double eleMissingHitsHeep32           = getPreCutValue1("eleMissingHitsHeep32"      );
  double eleUseEcalDrivenHeep32         = getPreCutValue1("eleUseEcalDrivenHeep32"    );

  CreateUserTH1D("dphi_met_ele", 100, 0, 3.14159 );

  //-----------------------------------------------------------------
  // Which algorithms to use?
  //-----------------------------------------------------------------

  int    eleAlgorithm = (int) getPreCutValue1("eleAlgorithm");

  //-----------------------------------------------------------------
  // Counters
  //-----------------------------------------------------------------

  int n_photon30__160404_163869 = 0;
  int n_photon30_HEEP__160404_163869 = 0;
  int n_photon30_HEEP_ele27__160404_163869 = 0;
  double eff_eleIDIso__160404_163869 = 0.;

  //## Maps 
  // first  = key   = RunNumber 
  // second = value = NumberOfEvents 
  map<int, int> MapTrg; 
  map<int, int> MapDen; 
  map<int, int> MapNum;
  map<int, int> MapNumErrSquare;
    
  /*//------------------------------------------------------------------
   *
   *
   *      
   *      Start analysis loop!
   *
   *
   *
   *///-----------------------------------------------------------------


  Long64_t nentries = fChain->GetEntries();
  //Long64_t nentries = 200000;
  STDOUT("analysisClass::Loop(): nentries = " << nentries);

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) { // Begin of loop over events
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry < 10 || jentry%1000 == 0) STDOUT("analysisClass::Loop(): jentry = " << jentry << "/" << nentries );
    
    //-----------------------------------------------------------------
    // Do pileup re-weighting, if necessary
    //  --> To be done after the skim, so commented out for now
    //-----------------------------------------------------------------
    
    // double event_weight = getPileupWeight ( PileUpInteractions, isData ) ;
    
    //-----------------------------------------------------------------
    // Get trigger information, if necessary
    //-----------------------------------------------------------------

    if ( isData ) { 
      getTriggers ( HLTKey, HLTInsideDatasetTriggerNames, HLTInsideDatasetTriggerDecisions,  HLTInsideDatasetTriggerPrescales ) ;
    }
    
    //-----------------------------------------------------------------
    // Selection: Electrons
    //-----------------------------------------------------------------    

    vector<int> v_idx_ele_PtCut_IDISO_STORE;
    vector<int> v_idx_ele_PtCut_IDISO_ANA;
    vector<int> v_idx_ele_IDISO;

    //Loop over electrons
    for(int iele=0; iele<ElectronPt->size(); iele++){
      
      int passEleSel = 0;
      int isBarrel = 0;
      int isEndcap = 0;
      
      if( fabs( ElectronSCEta->at(iele) ) < eleEta_bar )       isBarrel = 1;
      if( fabs( ElectronSCEta->at(iele) ) > eleEta_end_min &&
	  fabs( ElectronSCEta->at(iele) ) < eleEta_end_max )   isEndcap = 1;

      //-----------------------------------------------------------------    
      // HEEP ID application 3.1
      //-----------------------------------------------------------------    

      if ( eleAlgorithm == 1 ) { 
	
	if(isBarrel) {		

	  if(   fabs(ElectronDeltaEtaTrkSC->at(iele)) < eleDeltaEtaTrkSCHeep_bar 
		&& fabs(ElectronDeltaPhiTrkSC->at(iele)) < eleDeltaPhiTrkSCHeep_bar 
		&& ElectronHoE->at(iele) < eleHoEHeep_bar 
		&& (ElectronE2x5OverE5x5->at(iele) >eleE2x5OverE5x5Heep_bar || ElectronE1x5OverE5x5->at(iele) > eleE1x5OverE5x5Heep_bar ) 
		&& ( ElectronEcalIsoDR03->at(iele)+ElectronHcalIsoD1DR03->at(iele) ) < eleEcalHcalIsoHeep_1_bar + eleEcalHcalIsoHeep_2_bar*ElectronPt->at(iele)
		&& ElectronTrkIsoDR03->at(iele) <eleTrkIsoHeep_bar 
		&& ElectronMissingHits->at(iele) == 0 
		&& ElectronHasEcalDrivenSeed->at(iele)
	     )
	    passEleSel = 1;		

	}//end barrel
	
	if(isEndcap) {		
	  
	  int passEcalHcalIsoCut=0;
	  if(ElectronPt->at(iele) < eleEcalHcalIsoHeep_PTthr_end && 
	     (ElectronEcalIsoDR03->at(iele)+ElectronHcalIsoD1DR03->at(iele)) < eleEcalHcalIsoHeep_1_end) 
	    passEcalHcalIsoCut=1;
	  if(ElectronPt->at(iele) > eleEcalHcalIsoHeep_PTthr_end && 
	     (ElectronEcalIsoDR03->at(iele)+ElectronHcalIsoD1DR03->at(iele)) < eleEcalHcalIsoHeep_1_end+eleEcalHcalIsoHeep_2_end*(ElectronPt->at(iele)-eleEcalHcalIsoHeep_PTthr_end) ) 
	    passEcalHcalIsoCut=1;
	  
	  if(fabs(ElectronDeltaEtaTrkSC->at(iele)) < eleDeltaEtaTrkSCHeep_end 
	     && fabs(ElectronDeltaPhiTrkSC->at(iele)) < eleDeltaPhiTrkSCHeep_end 
	     && ElectronHoE->at(iele) < eleHoEHeep_end 
	     && ElectronSigmaIEtaIEta->at(iele) < eleSigmaIetaIetaHeep_end 
	     && passEcalHcalIsoCut == 1
	     && ElectronHcalIsoD2DR03->at(iele) < eleHcalIsoD2Heep_end 
	     && ElectronTrkIsoDR03->at(iele) < eleTrkIsoHeep_end 
	     && ElectronMissingHits->at(iele) == 0 
	     && ElectronHasEcalDrivenSeed->at(iele)
	     )
	    passEleSel = 1;
	  
	}//end endcap
      }

      //-----------------------------------------------------------------    
      // WP80 ID application
      //-----------------------------------------------------------------    

      else if ( eleAlgorithm == 2 ) { 
	// ecal driven	    
	if( eleUseEcalDrivenWP && !ElectronHasEcalDrivenSeed->at(iele) ) continue;
	
	// isolation
	double ElectronCombRelIsoWP_bar  =  ( ElectronTrkIsoDR03->at(iele) 
					    + max( 0., ElectronEcalIsoDR03->at(iele) - 1. ) 
					    + ElectronHcalIsoDR03FullCone->at(iele) 
					    - rhoIso*TMath::Pi()*0.3*0.3 
					    ) / ElectronPt->at(iele) ;
	
	double ElectronCombRelIsoWP_end  =  ( ElectronTrkIsoDR03->at(iele) 
					      + ElectronEcalIsoDR03->at(iele) 
					      + ElectronHcalIsoDR03FullCone->at(iele) 
					      - rhoIso*TMath::Pi()*0.3*0.3 
					      ) / ElectronPt->at(iele) ;
	
	// conversions
	int isPhotConv = 0;
	if(eleUseHasMatchConvWP) {
	  if( ElectronHasMatchedConvPhot->at(iele) ) 
	    isPhotConv = 1;
	}
	else {
	  if( ElectronDist->at(iele) < eleDistWP && ElectronDCotTheta->at(iele) < eleDCotThetaWP )	
	    isPhotConv = 1;
	}

	if(isBarrel) {
	  
	  if( ElectronMissingHits->at(iele) <= eleMissingHitsWP              && 
	      isPhotConv == 0					             && 
	      ElectronCombRelIsoWP_bar < eleCombRelIsoWP_bar		     && 
	      ElectronSigmaIEtaIEta->at(iele) < eleSigmaIetaIetaWP_bar       && 
	      fabs(ElectronDeltaPhiTrkSC->at(iele)) < eleDeltaPhiTrkSCWP_bar && 
	      fabs(ElectronDeltaEtaTrkSC->at(iele)) < eleDeltaEtaTrkSCWP_bar  )
	    passEleSel = 1;		
	  
	}//end barrel
	
	if(isEndcap) {		
	    
	  if( ElectronMissingHits->at(iele) == eleMissingHitsWP              && 
	      isPhotConv == 0						     && 
	      ElectronCombRelIsoWP_end < eleCombRelIsoWP_end		     && 
	      ElectronSigmaIEtaIEta->at(iele) < eleSigmaIetaIetaWP_end 	     && 
	      fabs(ElectronDeltaPhiTrkSC->at(iele)) < eleDeltaPhiTrkSCWP_end && 
	      fabs(ElectronDeltaEtaTrkSC->at(iele)) < eleDeltaEtaTrkSCWP_end  )
	    passEleSel = 1;		
	  
	}//end endcap	
      }

      //-----------------------------------------------------------------    
      // HEEP ID application 3.2
      //-----------------------------------------------------------------    

      else if ( eleAlgorithm == 3 ) { 
	
	if(isBarrel) {
		
	  if(   fabs(ElectronDeltaEtaTrkSC->at(iele)) < eleDeltaEtaTrkSCHeep32_bar 
		&& fabs(ElectronDeltaPhiTrkSC->at(iele)) < eleDeltaPhiTrkSCHeep32_bar 
		&& ElectronHoE->at(iele) < eleHoEHeep32_bar 
		&& (ElectronE2x5OverE5x5->at(iele) >eleE2x5OverE5x5Heep32_bar || ElectronE1x5OverE5x5->at(iele) > eleE1x5OverE5x5Heep32_bar ) 
		&& ( ElectronEcalIsoDR03->at(iele)+ElectronHcalIsoD1DR03->at(iele) ) < eleEcalHcalIsoHeep32_1_bar + eleEcalHcalIsoHeep32_2_bar*ElectronPt->at(iele)
		&& ElectronTrkIsoDR03->at(iele) <eleTrkIsoHeep32_bar 
		&& ElectronMissingHits->at(iele) == 0 
		&& ElectronHasEcalDrivenSeed->at(iele)
	     )
	    passEleSel = 1;		

	}//end barrel
	
	if(isEndcap) {		
	  
	  int passEcalHcalIsoCut=0;
	  if(ElectronPt->at(iele) < eleEcalHcalIsoHeep32_PTthr_end && 
	     (ElectronEcalIsoDR03->at(iele)+ElectronHcalIsoD1DR03->at(iele)) < eleEcalHcalIsoHeep32_1_end) 
	    passEcalHcalIsoCut=1;
	  if(ElectronPt->at(iele) > eleEcalHcalIsoHeep32_PTthr_end && 
	     (ElectronEcalIsoDR03->at(iele)+ElectronHcalIsoD1DR03->at(iele)) < eleEcalHcalIsoHeep32_1_end+eleEcalHcalIsoHeep32_2_end*(ElectronPt->at(iele)-eleEcalHcalIsoHeep32_PTthr_end) ) 
	    passEcalHcalIsoCut=1;
	  
	  if(fabs(ElectronDeltaEtaTrkSC->at(iele)) < eleDeltaEtaTrkSCHeep32_end 
	     && fabs(ElectronDeltaPhiTrkSC->at(iele)) < eleDeltaPhiTrkSCHeep32_end 
	     && ElectronHoE->at(iele) < eleHoEHeep32_end 
	     && ElectronSigmaIEtaIEta->at(iele) < eleSigmaIetaIetaHeep32_end 
	     && passEcalHcalIsoCut == 1
	     //&& ElectronHcalIsoD2DR03->at(iele) < eleHcalIsoD2Heep32_end 
	     && ElectronTrkIsoDR03->at(iele) < eleTrkIsoHeep32_end 
	     && ElectronMissingHits->at(iele) == 0 
	     && ElectronHasEcalDrivenSeed->at(iele)
	     )
	    passEleSel = 1;
	  
	}//end endcap
      }
      
      if ( passEleSel ) { 
	v_idx_ele_IDISO.push_back ( iele ) ;
	if ( ElectronPt -> at (iele) >= ele_PtCut_STORE ) v_idx_ele_PtCut_IDISO_STORE.push_back ( iele ) ;
	if ( ElectronPt -> at (iele) >= ele_PtCut_ANA   ) v_idx_ele_PtCut_IDISO_ANA  .push_back ( iele ) ;
      }

    }
    
    //-----------------------------------------------------------------
    // Fill your single-object variables with values
    //-----------------------------------------------------------------

    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();
    
    fillVariableWithValue( "PassJSON", passJSON(run, ls, isData) );    
    // Set the value of the variableNames listed in the cutFile to their current value

    //event info
    fillVariableWithValue( "isData"   , isData     ) ;
    fillVariableWithValue( "bunch"    , bunch      ) ;
    fillVariableWithValue( "event"    , event      ) ;
    fillVariableWithValue( "ls"       , ls         ) ;
    fillVariableWithValue( "orbit"    , orbit      ) ;
    fillVariableWithValue( "run"      , run        ) ;

    // Trigger (L1 and HLT)
    if(isData==true)
      {
	fillVariableWithValue( "PassBPTX0", isBPTX0 ) ;
	fillVariableWithValue( "PassPhysDecl", isPhysDeclared ) ;
      }
    else
      {
	fillVariableWithValue( "PassBPTX0", true ) ;
	fillVariableWithValue( "PassPhysDecl", true ) ;
      }
    //triggerFired ("HLT_Photon30_CaloIdVL_v1")

    //Event filters at RECO level
    fillVariableWithValue( "PassBeamScraping", !isBeamScraping ) ;
    fillVariableWithValue( "PassPrimaryVertex", isPrimaryVertex ) ;
    fillVariableWithValue( "PassHBHENoiseFilter", passHBHENoiseFilter ) ;
    fillVariableWithValue( "PassBeamHaloFilterLoose", passBeamHaloFilterLoose ) ;
    fillVariableWithValue( "PassBeamHaloFilterTight", passBeamHaloFilterTight ) ;
    fillVariableWithValue( "PassTrackingFailure", !isTrackingFailure ) ;
    fillVariableWithValue( "PassCaloBoundaryDRFilter", passCaloBoundaryDRFilter ) ;
    fillVariableWithValue( "PassEcalMaskedCellDRFilter", passEcalMaskedCellDRFilter ) ;
    
    // Evaluate cuts (but do not apply them)
    evaluateCuts();

    double met = PFMET -> at(0);
    double met_phi = PFMETPhi -> at(0);

    //Basic Event Selection
    if( passedCut("PassJSON") 
	&& passedCut("PassBPTX0") 
	&& passedCut("PassBeamScraping") 
	&& passedCut("PassPrimaryVertex")
	&& passedCut("PassHBHENoiseFilter")
	&& passedCut("PassBeamHaloFilterTight") 
	)
      {      

	//### Fill Maps (full run range)
	
	if( triggerFired ("HLT_Photon30_CaloIdVL_v1") 
	    || triggerFired ("HLT_Photon30_CaloIdVL_v2") 
	    || triggerFired ("HLT_Photon30_CaloIdVL_v3") 
	    || triggerFired ("HLT_Photon30_CaloIdVL_v4") 
	    || triggerFired ("HLT_Photon30_CaloIdVL_v5") 
	    || triggerFired ("HLT_Photon30_CaloIdVL_v6") 
	    || triggerFired ("HLT_Photon30_CaloIdVL_v7") 
	    || triggerFired ("HLT_Photon30_CaloIdVL_v8") 
	    )
	  {
	
	    map<int,int>::iterator MapTrgIt = MapTrg.find(run);
	    map<int,int>::iterator MapTrgItEnd = MapTrg.end();
	    if( MapTrgIt == MapTrgItEnd )
	      {
		MapTrg.insert(pair<int,int>(run,1));
		//MapDen[run]=1;
	      }
	    else
	      MapTrgIt->second++;

	    TLorentzVector my_v_met, my_ele;
	    double my_delta_phi;
	    my_v_met.SetPtEtaPhiM ( met, 0.0, met_phi, 0.0);
	    
	    if( v_idx_ele_PtCut_IDISO_ANA.size()==1){
	      my_ele.SetPtEtaPhiM ( ElectronPt  -> at (v_idx_ele_PtCut_IDISO_ANA[0]),
				    ElectronEta -> at (v_idx_ele_PtCut_IDISO_ANA[0]),
				    ElectronPhi -> at (v_idx_ele_PtCut_IDISO_ANA[0]),
				    0.0 );
	      my_delta_phi = fabs(my_v_met.DeltaPhi ( my_ele ) );
	    }

	    
	    if( v_idx_ele_PtCut_IDISO_ANA.size()==1 && met > 30.0 && my_delta_phi > 2.5 )
	      {
		//denominator
		map<int,int>::iterator MapDenIt = MapDen.find(run);
		map<int,int>::iterator MapDenItEnd = MapDen.end();

		FillUserTH1D("dphi_met_ele", my_delta_phi );
		

		if( MapDenIt == MapDenItEnd )
		  {
		    MapDen.insert(pair<int,int>(run,1));
		    //MapDen[run]=1;
		  }
		else
		  MapDenIt->second++;

		CreateAndFillUserTH1D("ElePt_AfterPhoton30", 100, 0, 1000, ElectronPt -> at (v_idx_ele_PtCut_IDISO_ANA[0]) );
		CreateAndFillUserTH1D("MET_AfterPhoton30", 100, 0, 1000, PFMET->at(0) );

		//numerator
		map<int,int>::iterator MapNumIt = MapNum.find(run);
		map<int,int>::iterator MapNumErrSquareIt = MapNumErrSquare.find(run);
		map<int,int>::iterator MapNumItEnd = MapNum.end();

		//-- 160404-161176
		if( triggerFired ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1") )
		  {
		    if( MapNumIt == MapNumItEnd )
		      {
			MapNum.insert(pair<int,int>( run , triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1") ));
			//MapNum[run]=1;
			MapNumErrSquare.insert(pair<int,int>( run , pow(triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1"),2) ));
		      }
		    else
		      {
			MapNumIt->second += triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1");
			MapNumErrSquareIt->second += pow(triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1"),2);
		      }
		  }

		//-- 161216-163261
		if( triggerFired ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2") )
		  {
		    if( MapNumIt == MapNumItEnd )
		      {
			MapNum.insert(pair<int,int>( run , triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2") ));
			//MapNum[run]=1;
			MapNumErrSquare.insert(pair<int,int>( run , pow(triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2"),2) ));
		      }
		    else
		      {
			MapNumIt->second += triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2");
			MapNumErrSquareIt->second += pow(triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2"),2);
		      }
		  }

		//-- 163269-163869
		if( triggerFired ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3") )
		  {
		    if( MapNumIt == MapNumItEnd )
		      {
			MapNum.insert(pair<int,int>( run , triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3") ));
			//MapNum[run]=1;
			MapNumErrSquare.insert(pair<int,int>( run , pow(triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3"),2) ));
		      }
		    else
		      {
			MapNumIt->second += triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3");
			MapNumErrSquareIt->second += pow(triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3"),2);
		      }
		  }

		//-- 165088-165633
		if( triggerFired ("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3") )
		  {
		    if( MapNumIt == MapNumItEnd )
		      {
			MapNum.insert(pair<int,int>( run , triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3") ));
			//MapNum[run]=1;
			MapNumErrSquare.insert(pair<int,int>( run , pow(triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3"),2) ));
		      }
		    else
		      {
			MapNumIt->second += triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3");
			MapNumErrSquareIt->second += pow(triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3"),2);
		      }
		  }

		//-- 165970-166967
		if( triggerFired ("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4") )
		  {
		    if( MapNumIt == MapNumItEnd )
		      {
			MapNum.insert(pair<int,int>( run , triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4") ));
			//MapNum[run]=1;
			MapNumErrSquare.insert(pair<int,int>( run , pow(triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4"),2) ));
		      }
		    else
		      {
			MapNumIt->second += triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4");
			MapNumErrSquareIt->second += pow(triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4"),2);
		      }
		  }

		//-- 167039-167913
		if( triggerFired ("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v5") )
		  {
		    if( MapNumIt == MapNumItEnd )
		      {
			MapNum.insert(pair<int,int>( run , triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v5") ));
			//MapNum[run]=1;
			MapNumErrSquare.insert(pair<int,int>( run , pow(triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v5"),2) ));
		      }
		    else
		      {
			MapNumIt->second += triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v5");
			MapNumErrSquareIt->second += pow(triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v5"),2);
		      }
		  }

		//-- 170249-173198
		if( triggerFired ("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v6") )
		  {
		    if( MapNumIt == MapNumItEnd )
		      {
			MapNum.insert(pair<int,int>( run , triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v6") ));
			//MapNum[run]=1;
			MapNumErrSquare.insert(pair<int,int>( run , pow(triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v6"),2) ));
		      }
		    else
		      {
			MapNumIt->second += triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v6");
			MapNumErrSquareIt->second += pow(triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v6"),2);
		      }
		  }

		//-- 173236-178380
		if( triggerFired ("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v7") )
		  {
		    if( MapNumIt == MapNumItEnd )
		      {
			MapNum.insert(pair<int,int>( run , triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v7") ));
			//MapNum[run]=1;
			MapNumErrSquare.insert(pair<int,int>( run , pow(triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v7"),2) ));
		      }
		    else
		      {
			MapNumIt->second += triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v7");
			MapNumErrSquareIt->second += pow(triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v7"),2);
		      }
		  }

		//-- 178420-179889
		if( triggerFired ("HLT_Ele27_WP80_v2") )
		  {
		    if( MapNumIt == MapNumItEnd )
		      {
			MapNum.insert(pair<int,int>( run , triggerPrescale("HLT_Ele27_WP80_v2") ));
			//MapNum[run]=1;
			MapNumErrSquare.insert(pair<int,int>( run , pow(triggerPrescale("HLT_Ele27_WP80_v2"),2) ));
		      }
		    else
		      {
			MapNumIt->second += triggerPrescale("HLT_Ele27_WP80_v2");
			MapNumErrSquareIt->second += pow(triggerPrescale("HLT_Ele27_WP80_v2"),2);
		      }
		  }

		//-- 179959-180252
		if( triggerFired ("HLT_Ele27_WP80_v3") )
		  {
		    if( MapNumIt == MapNumItEnd )
		      {
			MapNum.insert(pair<int,int>( run , triggerPrescale("HLT_Ele27_WP80_v3") ));
			//MapNum[run]=1;
			MapNumErrSquare.insert(pair<int,int>( run , pow(triggerPrescale("HLT_Ele27_WP80_v3"),2) ));
		      }
		    else
		      {
			MapNumIt->second += triggerPrescale("HLT_Ele27_WP80_v3");
			MapNumErrSquareIt->second += pow(triggerPrescale("HLT_Ele27_WP80_v3"),2);
		      }
		  }

		// HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1    160404-161176
		// HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2    161216-163261
		// (HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1)  161216-163261
		// HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3    163269-163869
		// (HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2)  163269-163869
		// HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3    165088-165633
		// HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4    165970-166967
		// HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v5    167039-167913
		// HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v6    170249-173198
		// HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v7    173236-178380
		// HLT_Ele27_WP80_v2				 178420-179889
		// HLT_Ele27_WP80_v3				 179959-180252

	      }//end 1 HEEP requirement

	  }//end photon trigger fired
	//###

	//### run range: 160404-163869 (TEST)
	if( triggerFired ("HLT_Photon30_CaloIdVL_v1") 
	    || triggerFired ("HLT_Photon30_CaloIdVL_v2") 
	    || triggerFired ("HLT_Photon30_CaloIdVL_v3") 
	    )
	  {
	    n_photon30__160404_163869++;
	    
	    if( v_idx_ele_PtCut_IDISO_ANA.size()==1 )
	      {
		//denominator
		n_photon30_HEEP__160404_163869++; 	
		CreateAndFillUserTH1D("ElePt_AfterPhoton30__160404_163869", 100, 0, 1000, ElectronPt -> at (v_idx_ele_PtCut_IDISO_ANA[0]) );
		
		//numerator
		if( triggerFired ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1") )
		  n_photon30_HEEP_ele27__160404_163869 += triggerPrescale ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1" ); 
		
		if( triggerFired ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2") )
		  n_photon30_HEEP_ele27__160404_163869 += triggerPrescale ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2" ); 
		
		if( triggerFired ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3") )
		  n_photon30_HEEP_ele27__160404_163869 += triggerPrescale ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3" ); 
	      }//end 1 HEEP requirement
	  }//end photon trigger fired
	//###

      }//end Basic Event Selection


  } // End of loop over events

  //##################################################

  //## Printout (test)
  cout << "---- Test : run range = 160404-163869 ----" << endl;

  eff_eleIDIso__160404_163869 = double(n_photon30_HEEP_ele27__160404_163869) / double(n_photon30_HEEP__160404_163869);

  cout << "n_photon30__160404_163869: " << n_photon30__160404_163869 << endl;
  cout << "n_photon30_HEEP__160404_163869: " << n_photon30_HEEP__160404_163869 << endl;
  cout << "n_photon30_HEEP_ele27__160404_163869: " << n_photon30_HEEP_ele27__160404_163869 << endl;
  cout << "eff_eleIDIso__160404_163869: " << eff_eleIDIso__160404_163869 << endl;
  cout << endl;

  //## Printout (final results)
  cout << "---- Final Results (full run range) ----" << endl;

  map<int,int>::iterator MapTrgIt = MapTrg.begin();
  map<int,int>::iterator MapTrgItEnd = MapTrg.end();
  int sumTrg = 0 ;
  //   FILE * pFileTrg;
  //   pFileTrg = fopen ("MapTrg.txt","w");
  for(;MapTrgIt!=MapTrgItEnd;++MapTrgIt)
    {
      //cout << "run , N : " << MapTrgIt->first << " " << MapTrgIt->second << endl;
      sumTrg += MapTrgIt->second;
      cout << "MAPTRG: " << MapTrgIt->first << " " << MapTrgIt->second << " " << sqrt(MapTrgIt->second) << endl;
      //fprintf (pFileTrg, "%d %d %f\n",MapTrgIt->first,MapTrgIt->second,sqrt(MapTrgIt->second));
    }
  //fclose (pFileTrg);

  map<int,int>::iterator MapDenIt = MapDen.begin();
  map<int,int>::iterator MapDenItEnd = MapDen.end();
  int sumDen = 0 ;
  //   FILE * pFileDen;
  //   pFileDen = fopen ("MapDen.txt","w");
  for(;MapDenIt!=MapDenItEnd;++MapDenIt)
    {
      //cout << "run , N : " << MapDenIt->first << " " << MapDenIt->second << endl;
      sumDen += MapDenIt->second;
      cout << "MAPDEN: " << MapDenIt->first << " " << MapDenIt->second << " " << sqrt(MapDenIt->second) << endl;
      //fprintf (pFileDen, "%d %d %f\n",MapDenIt->first,MapDenIt->second,sqrt(MapDenIt->second));
    }
  //fclose (pFileDen);

  map<int,int>::iterator MapNumIt = MapNum.begin();
  map<int,int>::iterator MapNumErrSquareIt = MapNumErrSquare.begin();
  map<int,int>::iterator MapNumItEnd = MapNum.end();
  int sumNum = 0 ;
  float sumNumErrSquare = 0. ;
  //   FILE * pFileNum;
  //   pFileNum = fopen ("MapNum.txt","w");
  for(;MapNumIt!=MapNumItEnd;++MapNumIt)
    {            
      //cout << "run , N : " << MapNumIt->first << " " << MapNumIt->second << endl;
      sumNum += MapNumIt->second;
      sumNumErrSquare += MapNumErrSquareIt->second;
      cout << "MAPNUM: " << MapNumIt->first << " " << MapNumIt->second << " " << sqrt(MapNumErrSquareIt->second) << endl;
      //fprintf (pFileNum, "%d %d %f\n",MapNumIt->first,MapNumIt->second,sqrt(MapNumErrSquareIt->second));
      ++MapNumErrSquareIt;
    }
  //fclose (pFileNum);

  double eff = double(sumNum) / double(sumDen);

  cout << "sumTrg : " << sumTrg << endl;
  cout << "sumDen : " << sumDen << endl;
  cout << "sumNum +/- err: " << sumNum << " +/- " << sqrt(sumNumErrSquare) << endl;
  cout << "eff : "    << eff    << endl; 

  /*//------------------------------------------------------------------
   *
   *
   *      
   *      End analysis loop!
   *
   *
   *
   *///-----------------------------------------------------------------

  STDOUT("analysisClass::Loop() ends");
}
