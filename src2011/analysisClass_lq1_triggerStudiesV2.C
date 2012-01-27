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

  //-----------------------------------------------------------------
  // Which algorithms to use?
  //-----------------------------------------------------------------

  int    eleAlgorithm = (int) getPreCutValue1("eleAlgorithm");

  //-----------------------------------------------------------------
  // Counters
  //-----------------------------------------------------------------

  int    N_probe_PassEleOffline = 0;
  int    N_probe_PassEleOfflineAndWP80 = 0;
  int    N_probe_PassEleOfflineAndTag = 0;  

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
  //Long64_t nentries = 1000;
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

	bool PassHLT = false;
	if( triggerFired ("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v1") 
	    || triggerFired ("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v2") 
	    || triggerFired ("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v3") 
	    || triggerFired ("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v4") 
	    || triggerFired ("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v5") 
	    || triggerFired ("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v6") 
	    || triggerFired ("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v7") 
	    || triggerFired ("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v8") 
	    )
	  {
	    PassHLT = true;
	  }
	fillVariableWithValue( "PassHLT", PassHLT ) ;
      }
    else
      {
	fillVariableWithValue( "PassBPTX0", true ) ;
	fillVariableWithValue( "PassPhysDecl", true ) ;
	fillVariableWithValue( "PassHLT", true ) ;
      }


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


    //Basic Event Selection
    if( passedCut("PassJSON") 
	&& passedCut("PassBPTX0") 
	&& passedCut("PassBeamScraping") 
	&& passedCut("PassPrimaryVertex")
	&& passedCut("PassHBHENoiseFilter")
	&& passedCut("PassBeamHaloFilterTight") 
	&& passedCut("PassHLT")
	)
      {      
	//Loop over probes
	for(int iprobe=0; iprobe<HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilterPt->size(); iprobe++)
	  {

	    TLorentzVector probe;
	    probe.SetPtEtaPhiE(HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilterPt->at(iprobe),
			       HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilterEta->at(iprobe),
			       HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilterPhi->at(iprobe),
			       HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilterEnergy->at(iprobe)
			       );

	    CreateAndFillUserTH1D("Pt_Probe", 200, 0, 200, probe.Pt() );

	    //Loop over tags
 	    for(int itag=0; itag<HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilterPt->size(); itag++)
	      {

		TLorentzVector tag;
		tag.SetPtEtaPhiE(HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilterPt->at(itag),
				 HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilterEta->at(itag),
				 HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilterPhi->at(itag),
				 HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilterEnergy->at(itag)
				 );	    

		CreateAndFillUserTH1D("DR_ProbeVsTag", 100, 0, 10, probe.DeltaR(tag) );

		//-----------------
		//if the tag matches in DR with the probe --> move to the next tag candidate
		if( probe.DeltaR(tag) < 0.5) 
		  continue;
		//-----------------


		//Now we should have a good (tag-probe) pair


		bool IsProbeMatchedWithOfflineEle = false;
		bool IsProbeMatchedWithTriggerWP80 = false;
		bool IsProbeMatchedWithTriggerTag = false;
		bool IsTagMatchedWithOfflineEle = false;

		//Loop over offline electrons
		for(int iele=0; iele<v_idx_ele_PtCut_IDISO_ANA.size(); iele++)
		  {
		    TLorentzVector ele;
		    ele.SetPtEtaPhiE( ElectronPt->at(v_idx_ele_PtCut_IDISO_ANA[iele]),
				      ElectronEta->at(v_idx_ele_PtCut_IDISO_ANA[iele]),
				      ElectronPhi->at(v_idx_ele_PtCut_IDISO_ANA[iele]),
				      ElectronEnergy->at(v_idx_ele_PtCut_IDISO_ANA[iele])
				      );	    
		    
		    CreateAndFillUserTH1D("DR_ProbeVsEle", 100, 0, 10, probe.DeltaR(ele) );
		    CreateAndFillUserTH1D("DR_TagVsEle", 100, 0, 10, tag.DeltaR(ele) );

		    if( probe.DeltaR(ele) < 0.2 ) 
		      IsProbeMatchedWithOfflineEle = true;

		    if( tag.DeltaR(ele) < 0.2 ) 
		      IsTagMatchedWithOfflineEle = true;		    

		  }

		//Loop over trigger WP80 electrons
		for(int iwp80=0; iwp80<HLTEle27WP80TrackIsoFilterPt->size(); iwp80++)
		  {
		    TLorentzVector wp80;
		    wp80.SetPtEtaPhiE(HLTEle27WP80TrackIsoFilterPt->at(iwp80),
				      HLTEle27WP80TrackIsoFilterEta->at(iwp80),
				      HLTEle27WP80TrackIsoFilterPhi->at(iwp80),
				      HLTEle27WP80TrackIsoFilterEnergy->at(iwp80)
				      );	    
		    
		    CreateAndFillUserTH1D("DR_ProbeVsWP80", 100, 0, 10, probe.DeltaR(wp80) );

		    if( probe.DeltaR(wp80) < 0.2 ) 
		      IsProbeMatchedWithTriggerWP80 = true;
		  }

		//Loop over trigger tag
		for(int itag2=0; itag2<HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilterPt->size(); itag2++)
		  {
		    TLorentzVector tag2;
		    tag2.SetPtEtaPhiE(HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilterPt->at(itag2),
				      HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilterEta->at(itag2),
				      HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilterPhi->at(itag2),
				      HLTEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilterEnergy->at(itag2)
				      );	    

		    CreateAndFillUserTH1D("DR_ProbeVsTag2", 100, 0, 10, probe.DeltaR(tag2) );

		    if( probe.DeltaR(tag2) < 0.2 ) 
		      IsProbeMatchedWithTriggerTag = true;
		  }				

		TLorentzVector tag_probe_system;
		tag_probe_system = tag + probe;
		double mass = tag_probe_system.M();

		CreateAndFillUserTH1D("Mass_TagProbeSystem_NoCutOnProbe", 200, 0, 200, mass );

		if( IsProbeMatchedWithOfflineEle && IsTagMatchedWithOfflineEle )
		  {
		    CreateAndFillUserTH1D("Mass_TagProbeSystem_ProbeMatchedOfflineEle", 200, 0, 200, mass );

		    if( mass > 75 && mass < 95)
		      {
			CreateAndFillUserTH1D("Mass_TagProbeSystem_ForEfficiencyCalculation", 200, 0, 200, mass );

			N_probe_PassEleOffline++;
			
			if(IsProbeMatchedWithTriggerWP80)
			  N_probe_PassEleOfflineAndWP80++;
			
			if(IsProbeMatchedWithTriggerTag)
			  N_probe_PassEleOfflineAndTag++;
			
		      }//mass cut

		  }//IsProbeMatchedWithOfflineEle
		
	      }//end loop over tag


	  }//end loop over probe
	  


	//	fillVariableWithValue( "PassEcalMaskedCellDRFilter", passEcalMaskedCellDRFilter ) ;

	//triggerFired ("HLT_Photon30_CaloIdVL_v1") 


	//v_idx_ele_PtCut_IDISO_ANA.size()

	//CreateAndFillUserTH1D("ElePt_AfterPhoton30", 100, 0, 1000, ElectronPt -> at (v_idx_ele_PtCut_IDISO_ANA[0]) );

      }//end Basic Event Selection


  } // End of loop over events

  //##################################################

  /*//------------------------------------------------------------------
   *
   *
   *      
   *      End analysis loop!
   *
   *
   *
   *///-----------------------------------------------------------------

  //Printout
  double eff_WP80 = double(N_probe_PassEleOfflineAndWP80)/double(N_probe_PassEleOffline);
  double eff_Tag = double(N_probe_PassEleOfflineAndTag)/double(N_probe_PassEleOffline);

  cout << "N_probe_PassEleOffline: " << N_probe_PassEleOffline << endl;
  cout << "N_probe_PassEleOfflineAndWP80: " << N_probe_PassEleOfflineAndWP80 << endl;
  cout << "N_probe_PassEleOfflineAndTag: " << N_probe_PassEleOfflineAndTag << endl;
  cout << endl;
  cout << "eff_WP80: " << eff_WP80 << " +/- " << sqrt(eff_WP80 * (1- eff_WP80) / N_probe_PassEleOffline ) << endl;  
  cout << "eff_Tag: " << eff_Tag << " +/- " << sqrt(eff_Tag * (1- eff_Tag) / N_probe_PassEleOffline ) << endl;  

  STDOUT("analysisClass::Loop() ends");
}
