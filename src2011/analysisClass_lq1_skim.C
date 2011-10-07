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


//-----------------------------
//-----------------------------

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

  //-----------------------------------------------------------------
  // What kind of skim is this?
  //-----------------------------------------------------------------
    
  int    reducedSkimType        = (int) getPreCutValue1("reducedSkimType");

  //-----------------------------------------------------------------
  // CaloJet cut values
  //-----------------------------------------------------------------

  double caloJet_PtCut_STORE    = getPreCutValue2("calojet_PtCut"   );
  double caloJet_PtCut_ANA      = getPreCutValue1("calojet_PtCut" );
  double caloJet_EtaCut         = getPreCutValue1("calojet_EtaCut");
  double caloJet_eleDRCut       = getPreCutValue1("calojet_eleDRCut");
  double caloJet_elePtCut       = getPreCutValue1("calojet_elePtCut");

  if ( caloJet_PtCut_STORE    > caloJet_PtCut_ANA   ){
    STDOUT("ERROR in CaloJet cut values: all storage cuts must be looser or equal to analysis cuts.");
    exit(0) ;
  }

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

  //-----------------------------------------------------------------
  // Jet cut values
  //-----------------------------------------------------------------
  
  double jet_PtCut_STORE     = getPreCutValue2("jet_PtCut");
  double jet_PtCut_ANA       = getPreCutValue1("jet_PtCut");
  double jet_EtaCut          = getPreCutValue1("jet_EtaCut");
  double jet_TCHELCut        = getPreCutValue1("jet_TCHELCut");
  double jet_ele_DeltaRCut   = getPreCutValue1("jet_ele_DeltaRCut");

  double jet_PtCut_forMetScale     = getPreCutValue1("jet_PtCut_forMetScale");

  if ( jet_PtCut_STORE  > jet_PtCut_ANA    ){
    STDOUT("ERROR in Jet cut values: all storage cuts must be looser or equal to analysis cuts.");
    exit(0) ;
  }
  
  //-----------------------------------------------------------------
  // Muon cut values
  //-----------------------------------------------------------------

  double muon_PtCut_STORE   = getPreCutValue1("muon_PtCut");
  double muon_PtCut_ANA     = getPreCutValue2("muon_PtCut");
  double muNHits_minThresh  = getPreCutValue1("muNHits_minThresh");
  double muTrkD0Maximum     = getPreCutValue1("muTrkD0Maximum");
  double jet_muon_DeltaRCut = getPreCutValue1("jet_muon_DeltaRCut");

  if ( muon_PtCut_STORE > muon_PtCut_STORE ){ 
    STDOUT("ERROR in Muon cut values: all storage cuts must be looser or equal to analysis cuts.");
    exit(0) ;
  }

  //-----------------------------------------------------------------
  // Vertex cut values
  //-----------------------------------------------------------------

  double vertexMinimumNDOF = getPreCutValue1("vertexMinimumNDOF");
  double vertexMaxAbsZ     = getPreCutValue1("vertexMaxAbsZ");
  double vertexMaxd0       = getPreCutValue1("vertexMaxd0");
  
  //-----------------------------------------------------------------
  // Scaling values
  //-----------------------------------------------------------------
  
  double EleEnergyScale_EB    = getPreCutValue1("EleEnergyScale_EB");
  double EleEnergyScale_EE    = getPreCutValue1("EleEnergyScale_EE");
  double JetEnergyScale       = getPreCutValue1("JetEnergyScale");
  int    doJetOversmearing    = (int) getPreCutValue1("doJetOversmearing");
  double JetOversmearingSigma = getPreCutValue1("JetOversmearingSigma");

  //-----------------------------------------------------------------
  // Which algorithms to use?
  //-----------------------------------------------------------------

  int    jetAlgorithm = (int) getPreCutValue1("jetAlgorithm");
  int    metAlgorithm = (int) getPreCutValue1("metAlgorithm");
  int    eleAlgorithm = (int) getPreCutValue1("eleAlgorithm");
  
  // Random number generator for random scaling
  TRandom3 *randomNumGen = new TRandom3;
  randomNumGen->SetSeed();

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
  // Long64_t nentries = 5000;
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
    // Store variables: Jets
    //-----------------------------------------------------------------

    //## Define new jet collection
    std::auto_ptr<std::vector<double> >  JetPt          ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetPtRaw       ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetEnergy      ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetEta         ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetPhi         ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >     JetPassID      ( new std::vector<int>   ()  );
    std::auto_ptr<std::vector<double> >  JetTCHE        ( new std::vector<double>()  );

    std::auto_ptr<std::vector<double> >  caloJetPt      ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  caloJetPtRaw   ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  caloJetEnergy  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  caloJetEta     ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  caloJetPhi     ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >     caloJetPassID  ( new std::vector<int>   ()  );
    std::auto_ptr<std::vector<double> >  caloJetTCHE    ( new std::vector<double>()  );

    //std::auto_ptr<std::vector<double> >  JetChargedMuEnergyFraction  ( new std::vector<double>()  );

    if(jetAlgorithm==1) //PF jets
      {
	for (int ijet=0 ; ijet< PFJetPt->size() ; ijet++)
	  {
	    JetPt->push_back( PFJetPt->at(ijet) );
	    JetPtRaw->push_back( PFJetPtRaw->at(ijet) );
	    JetEnergy->push_back( PFJetEnergy->at(ijet) );
	    JetEta->push_back( PFJetEta->at(ijet) );
	    JetPhi->push_back( PFJetPhi->at(ijet) );
            JetPassID->push_back( PFJetPassLooseID->at(ijet) );
            JetTCHE->push_back( PFJetTrackCountingHighEffBTag->at(ijet) );
	    //JetChargedMuEnergyFraction->push_back( PFJetChargedMuEnergyFraction->at(ijet) );	    
	  }//end loop over pf jets
      }//end if "pf jets"

    for (int ijet=0 ; ijet < CaloJetPt->size() ; ijet++)
      {

  
	caloJetPt->push_back( CaloJetPt->at(ijet) );
	caloJetPtRaw->push_back( CaloJetPtRaw->at(ijet) );
	caloJetEnergy->push_back( CaloJetEnergy->at(ijet) );
	caloJetEta->push_back( CaloJetEta->at(ijet) );
	caloJetPhi->push_back( CaloJetPhi->at(ijet) );
	caloJetPassID->push_back( CaloJetPassLooseID->at(ijet) );
	caloJetTCHE->push_back( CaloJetTrackCountingHighEffBTag->at(ijet) );
	
	if(jetAlgorithm==2) //Calo jets
	  {
	    
	    JetPt->push_back( CaloJetPt->at(ijet) );
	    JetPtRaw->push_back( CaloJetPtRaw->at(ijet) );
	    JetEnergy->push_back( CaloJetEnergy->at(ijet) );
	    JetEta->push_back( CaloJetEta->at(ijet) );
	    JetPhi->push_back( CaloJetPhi->at(ijet) );
	    JetPassID->push_back( CaloJetPassLooseID->at(ijet) );
	    JetTCHE->push_back( CaloJetTrackCountingHighEffBTag->at(ijet) );
	  }

      }//end loop over calo jets

    //-----------------------------------------------------------------
    // Store variables: MET
    //-----------------------------------------------------------------

    //## Define new met collection
    double thisMET;
    double thisMETPhi;
    double thisSumET;
    double thisGenMET;
    double thisGenMETPhi;
    double thisGenSumET;

    if(isData==0)
      {
	thisGenMET    = GenMETTrue->at(0);
	thisGenMETPhi = GenMETPhiTrue->at(0);
	thisGenSumET  = GenSumETTrue->at(0);
      }
    else
      {
	thisGenMET    = -1;
	thisGenMETPhi = -1;
	thisGenSumET  = -1;
      }

    if(metAlgorithm==1) 	// --> PFMET
      {
	thisMET = PFMET->at(0);
	thisMETPhi = PFMETPhi->at(0);
	thisSumET = PFSumET->at(0);
      }
    if(metAlgorithm==2) 	// --> CaloMET
      {
	thisMET = CaloMET->at(0);
	thisMETPhi = CaloMETPhi->at(0);
	thisSumET = CaloSumET->at(0);
      }
    if(metAlgorithm==3) 	// --> PFMET (with type-1 corrections)
      {
	thisMET = PFMETType1Cor->at(0);
	thisMETPhi = PFMETPhiType1Cor->at(0);
	thisSumET = PFSumETType1Cor->at(0);
      }

    //-----------------------------------------------------------------
    // HEEP Electron Pt 
    //-----------------------------------------------------------------

    if(eleAlgorithm==1) {
	for(int iele=0; iele<ElectronPt->size(); iele++)
	  ElectronPt->at(iele) = ElectronPtHeep->at(iele);
    }
    
    //-----------------------------------------------------------------
    // Energy scaling: Electrons
    //-----------------------------------------------------------------

    //## EES and JES
    if( EleEnergyScale_EB != 1 || EleEnergyScale_EE != 1 )
      {
	for(int iele=0; iele<ElectronPt->size(); iele++)
	  {
	    if( fabs(ElectronEta->at(iele)) < eleEta_bar )
	      ElectronPt->at(iele) *= EleEnergyScale_EB;
	    if( fabs(ElectronEta->at(iele)) > eleEta_end_min && fabs(ElectronEta->at(iele)) < eleEta_end_max )
	      ElectronPt->at(iele) *= EleEnergyScale_EE;
	  }
      }

    //-----------------------------------------------------------------
    // Energy scaling: Jets
    //-----------------------------------------------------------------

    if( JetEnergyScale != 1 )
      { //use fix JES scaling passed from cut file

	for(int ijet=0; ijet<JetPt->size(); ijet++)
	  {
	    JetPt->at(ijet) *= JetEnergyScale;
	    JetEnergy->at(ijet) *= JetEnergyScale;
	  }
      }

    // Jet oversmearing
    std::auto_ptr<std::vector<double> >  JetOversmearingFactor ( new std::vector<double>(JetPt->size(), 1)  );
    if( !isData && doJetOversmearing )
      {
        for(int ijet=0; ijet<JetPt->size(); ijet++)
          {
            JetOversmearingFactor->at(ijet) = (1 + JetOversmearingSigma*randomNumGen->Gaus());
            JetPt->at(ijet) *= JetOversmearingFactor->at(ijet);
            JetEnergy->at(ijet) *= JetOversmearingFactor->at(ijet);
          }
      }

    //-----------------------------------------------------------------
    // Selection: Muons
    //-----------------------------------------------------------------

    vector<int> v_idx_muon_PtCut_IDISO_ANA;
    vector<int> v_idx_muon_PtCut_IDISO_STORE;

    // Loop over muons
    for(int imuon=0; imuon<MuonPt->size(); imuon++){

      bool ana_muon = true;

      if ( MuonPt         ->at(imuon)   <  muon_PtCut_STORE  ) continue;      
      if ( MuonTrkHits    ->at(imuon)   <  muNHits_minThresh ) continue;
      if ( fabs( MuonTrkD0->at(imuon) ) >= muTrkD0Maximum    ) continue;
      if ( MuonPassIso    ->at(imuon)   == 0                 ) continue;
      if ( MuonPassID     ->at(imuon)   == 0                 ) continue;
      if ( MuonPt         ->at(imuon)   <  muon_PtCut_ANA    ) ana_muon = false;

      v_idx_muon_PtCut_IDISO_STORE.push_back(imuon);
      if ( ana_muon ) v_idx_muon_PtCut_IDISO_ANA.push_back(imuon);
      
    }// end loop over muons

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
      // HEEP ID application
      //-----------------------------------------------------------------    

      if ( eleAlgorithm == 1 ) { 
	
	if(isBarrel) {		
	  if(   fabs(ElectronDeltaEtaTrkSC->at(iele)) < eleDeltaEtaTrkSCHeep_bar 
	     && fabs(ElectronDeltaPhiTrkSC->at(iele)) < eleDeltaPhiTrkSCHeep_bar 
	     && ElectronHoE->at(iele) < eleHoEHeep_bar 
	     && (ElectronE2x5OverE5x5->at(iele) >eleE2x5OverE5x5Heep_bar || ElectronE1x5OverE5x5->at(iele) > eleE1x5OverE5x5Heep_bar ) 
	     && ( ElectronEcalIsoDR03->at(iele)+ElectronHcalIsoD1DR03->at(iele) ) < eleEcalHcalIsoHeep_1_bar + eleEcalHcalIsoHeep_2_bar*ElectronPt->at(iele)
	     && ElectronTrkIsoDR03->at(iele) <eleTrkIsoHeep_bar 
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
      
      if ( passEleSel ) { 
	v_idx_ele_IDISO.push_back ( iele ) ;
	if ( ElectronPt -> at (iele) >= ele_PtCut_STORE ) v_idx_ele_PtCut_IDISO_STORE.push_back ( iele ) ;
	if ( ElectronPt -> at (iele) >= ele_PtCut_ANA   ) v_idx_ele_PtCut_IDISO_ANA  .push_back ( iele ) ;
      }

    }
    
    //-----------------------------------------------------------------
    // Selection: Analysis Jets
    //-----------------------------------------------------------------

    //## Jets
    vector<int> v_idx_jet_PtCut_STORE;
    vector<int> v_idx_jet_PtCut_ANA;

    vector<int> v_idx_jet_PtCut_noOverlap_ID_STORE;
    vector<int> v_idx_jet_PtCut_noOverlap_ID_TCHEL_STORE;

    vector<int> v_idx_jet_PtCut_noOverlap_ID_ANA;
    vector<int> v_idx_jet_PtCut_noOverlap_ID_TCHEL_ANA;

    vector<int> v_idx_calojet_PtCut_noOverlap_ID_STORE;
    vector<int> v_idx_calojet_PtCut_noOverlap_ID_ANA;

    // Loop over jets
    for(int ijet=0; ijet<JetPt->size(); ijet++)
      {
	bool ana_jet = true;

	//pT pre-cut on reco jets
	if ( JetPt->at(ijet) < jet_PtCut_STORE ) continue;
	if ( JetPt->at(ijet) < jet_PtCut_ANA   ) ana_jet = false;

	v_idx_jet_PtCut_STORE.push_back(ijet);
	if ( ana_jet ) v_idx_jet_PtCut_ANA.push_back(ijet);
      }

    //ele-jet overlap
    vector <int> jetFlagsEle(v_idx_jet_PtCut_STORE.size(), 0);
    int NjetflaggedEle = 0;
    for (int iele=0; iele<v_idx_ele_PtCut_IDISO_ANA.size(); iele++)
      {
	TLorentzVector ele;
        ele.SetPtEtaPhiM(ElectronPt ->at(v_idx_ele_PtCut_IDISO_ANA[iele]),
			 ElectronEta->at(v_idx_ele_PtCut_IDISO_ANA[iele]),
			 ElectronPhi->at(v_idx_ele_PtCut_IDISO_ANA[iele]),0);
	TLorentzVector jet;
	double minDR=9999.;
	int ijet_minDR = -1;
        for(int ijet=0; ijet<v_idx_jet_PtCut_STORE.size(); ijet++)
          {
	    if ( jetFlagsEle[ijet] == 1 )
	      continue;
            jet.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut_STORE    [ijet]),
			     JetEta->at(v_idx_jet_PtCut_STORE   [ijet]),
			     JetPhi->at(v_idx_jet_PtCut_STORE   [ijet]),
			     JetEnergy->at(v_idx_jet_PtCut_STORE[ijet]) );
	    double DR = jet.DeltaR(ele);
	    if (DR<minDR)
	      {
		minDR = DR;
		ijet_minDR = ijet;
	      }
	  }
	if ( minDR < jet_ele_DeltaRCut 
	     && ijet_minDR > -1 )
	  {
	    jetFlagsEle[ijet_minDR] = 1;
	    NjetflaggedEle++;
	  }
      }

    //muon-jet overlap
    vector <int> jetFlagsMuon(v_idx_jet_PtCut_STORE.size(), 0);
    int NjetflaggedMuon = 0;
    for (int imu=0; imu<v_idx_muon_PtCut_IDISO_ANA.size(); imu++)
      {
	TLorentzVector mu;
        mu.SetPtEtaPhiM(MuonPt ->at(v_idx_muon_PtCut_IDISO_ANA[imu]),
			MuonEta->at(v_idx_muon_PtCut_IDISO_ANA[imu]),
			MuonPhi->at(v_idx_muon_PtCut_IDISO_ANA[imu]),0);
	TLorentzVector jet;
	double minDR=9999.;
	int ijet_minDR = -1;
        for(int ijet=0; ijet<v_idx_jet_PtCut_STORE.size(); ijet++)
          {
	    if ( jetFlagsMuon[ijet] == 1 )
	      continue;
            jet.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut_STORE[ijet]),
			     JetEta->at(v_idx_jet_PtCut_STORE[ijet]),
			     JetPhi->at(v_idx_jet_PtCut_STORE[ijet]),
			     JetEnergy->at(v_idx_jet_PtCut_STORE[ijet]) );
	    double DR = jet.DeltaR(mu);
	    if (DR<minDR)
	      {
		minDR = DR;
		ijet_minDR = ijet;
	      }
	  }
	if ( minDR < jet_muon_DeltaRCut 
	     && ijet_minDR > -1 )
	  {
	    jetFlagsMuon[ijet_minDR] = 1;
	    NjetflaggedMuon++;
	  }
      }
    
    //Select final jets
    for(int ijet=0; ijet<v_idx_jet_PtCut_STORE.size(); ijet++) //pT pre-cut + no overlaps with electrons + jetID
      {
	bool passjetID = JetPassID->at(v_idx_jet_PtCut_STORE[ijet]);

	if( jetFlagsEle[ijet] == 0                                      /* NO overlap with electrons */ 
	    && jetFlagsMuon[ijet] == 0                                  /* NO overlap with muons */	    
	    && passjetID == true                                        /* pass JetID */
	    && fabs( JetEta->at(v_idx_jet_PtCut_STORE[ijet]) ) < jet_EtaCut 	/* pass Jet Eta cut */    
	    ) {                        
	  v_idx_jet_PtCut_noOverlap_ID_STORE.push_back(v_idx_jet_PtCut_STORE[ijet]);
	  if ( JetPt->at(v_idx_jet_PtCut_STORE[ijet]) > jet_PtCut_ANA )
	    v_idx_jet_PtCut_noOverlap_ID_ANA.push_back(v_idx_jet_PtCut_STORE[ijet]);
	}
	
	if( jetFlagsEle[ijet] == 0                                       /* NO overlap with electrons */ 
	    && jetFlagsMuon[ijet] == 0                                   /* NO overlap with muons */	    
	    && passjetID == true                                         /* pass JetID */
	    && fabs( JetEta ->at(v_idx_jet_PtCut_STORE[ijet]) ) < jet_EtaCut 	 /* pass Jet Eta cut */    
            && fabs( JetTCHE->at(v_idx_jet_PtCut_STORE[ijet]) ) > jet_TCHELCut /* TrackCountingHighEfficiency loose b-tag (see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagPerformance) */
	    ) {                        
	  v_idx_jet_PtCut_noOverlap_ID_TCHEL_STORE.push_back(v_idx_jet_PtCut_STORE[ijet]);	
	  if ( JetPt->at(v_idx_jet_PtCut_STORE[ijet]) > jet_PtCut_ANA )
	    v_idx_jet_PtCut_noOverlap_ID_TCHEL_ANA.push_back(v_idx_jet_PtCut_STORE[ijet]);	
	}
      } // End loop over jets

    //-----------------------------------------------------------------
    // Selection: CaloJets
    //-----------------------------------------------------------------
	
    for (int ijet = 0; ijet < caloJetPt->size(); ++ijet){
      bool passJetID   = caloJetPassID->at(ijet);
      bool passEtaCut  = fabs ( caloJetEta->at(ijet) ) < caloJet_EtaCut;
      bool passPtCut   = caloJetPt->at(ijet) > caloJet_PtCut_STORE;

      if ( !passJetID  ) continue;
      if ( !passEtaCut ) continue;
      if ( !passPtCut  ) continue;
      if ( v_idx_ele_IDISO.size() == 0 ) continue;

      TLorentzVector ele1, jet;

      ele1.SetPtEtaPhiM ( ElectronPt  -> at ( v_idx_ele_IDISO[0] ), 
			  ElectronEta -> at ( v_idx_ele_IDISO[0] ), 
			  ElectronPhi -> at ( v_idx_ele_IDISO[0] ), 0.0 ) ;

      if ( ele1.Pt() < caloJet_elePtCut ) continue;

      jet.SetPtEtaPhiM ( caloJetPt  -> at(ijet),
			 caloJetEta -> at(ijet),
			 caloJetPhi -> at(ijet), 0.0 ) ;
      
      double deltaR = jet.DeltaR ( ele1 ) ;
      
      if ( deltaR < caloJet_eleDRCut ) continue;

      v_idx_calojet_PtCut_noOverlap_ID_STORE.push_back ( ijet ) ;
      if ( caloJetPt->at(ijet) > caloJet_PtCut_ANA ) 
	v_idx_calojet_PtCut_noOverlap_ID_ANA.push_back ( ijet ) ;
    }

    //-----------------------------------------------------------------
    // MET energy scaling
    //-----------------------------------------------------------------

    if( JetEnergyScale != 1 || (!isData && doJetOversmearing) )
      { //use fix JES scaling passed from cut file

	TVector2 v_MET_old;
	TVector2 v_MET_new;

	//use only good jets (after electron-jet and muon-jet overlap) for re-doing MET
	for(int ijet=0; ijet<v_idx_jet_PtCut_noOverlap_ID_ANA.size(); ijet++)
	  {
	    TVector2 v_jet_pt_old;
	    TVector2 v_jet_pt_new;
	    v_jet_pt_old.SetMagPhi( JetPt->at(v_idx_jet_PtCut_noOverlap_ID_ANA[ijet])/JetEnergyScale/JetOversmearingFactor->at(v_idx_jet_PtCut_noOverlap_ID_ANA[ijet]) , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID_ANA[ijet]) );
	    v_jet_pt_new.SetMagPhi( JetPt->at(v_idx_jet_PtCut_noOverlap_ID_ANA[ijet]) , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID_ANA[ijet]) );
	    if ( v_jet_pt_old.Mod() < jet_PtCut_forMetScale ) continue;
	    v_MET_new += v_jet_pt_old - v_jet_pt_new;
	  }

	v_MET_old.SetMagPhi( thisMET , thisMETPhi );
	v_MET_new += v_MET_old;
 	thisMET = v_MET_new.Mod();
 	thisMETPhi = v_MET_new.Phi();

      }

    //-----------------------------------------------------------------
    // Selection: vertices
    //-----------------------------------------------------------------

    vector<int> v_idx_vertex_good;
    // loop over vertices
    for(int ivertex = 0; ivertex<VertexChi2->size(); ivertex++){
      if ( !(VertexIsFake->at(ivertex))
    	   && VertexNDF->at(ivertex) > vertexMinimumNDOF
    	   && fabs( VertexZ->at(ivertex) ) <= vertexMaxAbsZ
    	   && fabs( VertexRho->at(ivertex) ) <= vertexMaxd0 )
    	{
    	  v_idx_vertex_good.push_back(ivertex);
    	  //STDOUT("v_idx_vertex_good.size = "<< v_idx_vertex_good.size() );
    	}
    }

    //-----------------------------------------------------------------
    // Fill your single-object variables with values
    //-----------------------------------------------------------------

    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();
    
    // fillVariableWithValue( "PassJSON", passJSON(run, ls, isData) );    
    // Set the value of the variableNames listed in the cutFile to their current value

    //event info
    fillVariableWithValue( "isData"   , isData     ) ;
    fillVariableWithValue( "bunch"    , bunch      ) ;
    fillVariableWithValue( "event"    , event      ) ;
    fillVariableWithValue( "ls"       , ls         ) ;
    fillVariableWithValue( "orbit"    , orbit      ) ;
    fillVariableWithValue( "run"      , run        ) ;
    fillVariableWithValue( "ProcessID", ProcessID  ) ;
    fillVariableWithValue( "PtHat"    , PtHat      ) ;

    fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_1"     , -999 );
						              	   
    fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_1"     , -999 );
    fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_2"     , -999 );
    fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_3"     , -999 );
    fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_4"     , -999 );
    fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_5"     , -999 );
    fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_6"     , -999 );
    fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_7"     , -999 );
						              	   
    fillVariableWithValue("H_27_CIdVT_CIsT_TIdT_TIsT_1"     , -999 );
    fillVariableWithValue("H_27_CIdVT_CIsT_TIdT_TIsT_2"     , -999 );
    fillVariableWithValue("H_27_CIdVT_CIsT_TIdT_TIsT_3"     , -999 );
						              	   
    fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_1"     , -999 );
    fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_2"     , -999 );
    fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_3"     , -999 );
    fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_4"     , -999 );
    fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_5"     , -999 );
    fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_6"     , -999 );
    fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_7"     , -999 );
						              	   
    fillVariableWithValue("H_42_CIdVT_CIsT_TIdT_TIsT_1"     , -999 );
						              	   
    fillVariableWithValue("H_25_WP80_PFMT40_1"              , -999 );
						              	   
    fillVariableWithValue("H_27_WP80_PFMT50_1"              , -999 );
    fillVariableWithValue("H_27_WP80_PFMT50_2"              , -999 );
    fillVariableWithValue("H_27_WP80_PFMT50_3"              , -999 );
    fillVariableWithValue("H_27_WP80_PFMT50_4"              , -999 );

    fillVariableWithValue ("H_17_8_CIdT_TIdVL_CIsVL_TIsVL_1", -999 );
    fillVariableWithValue ("H_17_8_CIdT_TIdVL_CIsVL_TIsVL_2", -999 );
    fillVariableWithValue ("H_17_8_CIdT_TIdVL_CIsVL_TIsVL_3", -999 );
    fillVariableWithValue ("H_17_8_CIdT_TIdVL_CIsVL_TIsVL_4", -999 );
    fillVariableWithValue ("H_17_8_CIdT_TIdVL_CIsVL_TIsVL_5", -999 );

    fillVariableWithValue ("H_17_8_CIdT_CIsVL_TIdVL_TIsVL_5", -999 );
    fillVariableWithValue ("H_17_8_CIdT_CIsVL_TIdVL_TIsVL_6", -999 );
    fillVariableWithValue ("H_17_8_CIdT_CIsVL_TIdVL_TIsVL_7", -999 );
    fillVariableWithValue ("H_17_8_CIdT_CIsVL_TIdVL_TIsVL_8", -999 );

    fillVariableWithValue ("H_17_8_CIdL_CIsVL_1"            , -999 );
    fillVariableWithValue ("H_17_8_CIdL_CIsVL_2"            , -999 );
    fillVariableWithValue ("H_17_8_CIdL_CIsVL_3"            , -999 );
    fillVariableWithValue ("H_17_8_CIdL_CIsVL_4"            , -999 );
    fillVariableWithValue ("H_17_8_CIdL_CIsVL_5"            , -999 );
    fillVariableWithValue ("H_17_8_CIdL_CIsVL_6"            , -999 );
    fillVariableWithValue ("H_17_8_CIdL_CIsVL_7"            , -999 );
    fillVariableWithValue ("H_17_8_CIdL_CIsVL_8"            , -999 );
    
    fillVariableWithValue ("H_32_17_CIdT_CIsT_TIdT_TIsT_1"  , -999 );

    
    if ( isData ){
      if ( triggerFired ("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1") ) { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_1",triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1")); } else { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_1",-1 * triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1")); }
      			 									                                 												                               			     						
      if ( triggerFired ("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1") ) { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_1",triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1")); } else { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_1",-1 * triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1")); }
      if ( triggerFired ("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2") ) { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_2",triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2")); } else { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_2",-1 * triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2")); }
      if ( triggerFired ("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3") ) { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_3",triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3")); } else { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_3",-1 * triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3")); }
      if ( triggerFired ("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4") ) { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_4",triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4")); } else { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_4",-1 * triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4")); }
      if ( triggerFired ("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v5") ) { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_5",triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v5")); } else { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_5",-1 * triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v5")); }
      if ( triggerFired ("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v6") ) { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_6",triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v6")); } else { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_6",-1 * triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v6")); }
      if ( triggerFired ("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v7") ) { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_7",triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v7")); } else { fillVariableWithValue("H_15_CIdVT_CIsT_TIdT_TIsT_7",-1 * triggerPrescale("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v7")); }
      			 									                                 												                               			     						
      if ( triggerFired ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1") ) { fillVariableWithValue("H_27_CIdVT_CIsT_TIdT_TIsT_1",triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1")); } else { fillVariableWithValue("H_27_CIdVT_CIsT_TIdT_TIsT_1",-1 * triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1")); }
      if ( triggerFired ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2") ) { fillVariableWithValue("H_27_CIdVT_CIsT_TIdT_TIsT_2",triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2")); } else { fillVariableWithValue("H_27_CIdVT_CIsT_TIdT_TIsT_2",-1 * triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2")); }
      if ( triggerFired ("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3") ) { fillVariableWithValue("H_27_CIdVT_CIsT_TIdT_TIsT_3",triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3")); } else { fillVariableWithValue("H_27_CIdVT_CIsT_TIdT_TIsT_3",-1 * triggerPrescale("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3")); }
      			 									                                 												                               			     						
      if ( triggerFired ("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1") ) { fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_1",triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1")); } else { fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_1",-1 * triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1")); }
      if ( triggerFired ("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2") ) { fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_2",triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2")); } else { fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_2",-1 * triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2")); }
      if ( triggerFired ("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3") ) { fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_3",triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3")); } else { fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_3",-1 * triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3")); }
      if ( triggerFired ("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4") ) { fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_4",triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4")); } else { fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_4",-1 * triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4")); }
      if ( triggerFired ("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v5") ) { fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_5",triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v5")); } else { fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_5",-1 * triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v5")); }
      if ( triggerFired ("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v6") ) { fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_6",triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v6")); } else { fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_6",-1 * triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v6")); }
      if ( triggerFired ("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v7") ) { fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_7",triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v7")); } else { fillVariableWithValue("H_32_CIdVT_CIsT_TIdT_TIsT_7",-1 * triggerPrescale("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v7")); }
      			 									                                 												                              			     						
      if ( triggerFired ("HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1") ) { fillVariableWithValue("H_42_CIdVT_CIsT_TIdT_TIsT_1",triggerPrescale("HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1")); } else { fillVariableWithValue("H_42_CIdVT_CIsT_TIdT_TIsT_1",-1 * triggerPrescale("HLT_Ele42_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1")); }
      			                                                
      if ( triggerFired ("HLT_Ele25_WP80_PFMT40_v1"                     ) ) { fillVariableWithValue("H_25_WP80_PFMT40_1"         ,triggerPrescale("HLT_Ele25_WP80_PFMT40_v1"                     )); } else { fillVariableWithValue("H_25_WP80_PFMT40_1",-1 * triggerPrescale("HLT_Ele25_WP80_PFMT40_v1")); }
      			                                                			    																  						                          
      if ( triggerFired ("HLT_Ele27_WP80_PFMT50_v1"                     ) ) { fillVariableWithValue("H_27_WP80_PFMT50_1"         ,triggerPrescale("HLT_Ele27_WP80_PFMT50_v1"                     )); } else { fillVariableWithValue("H_27_WP80_PFMT50_1",-1 * triggerPrescale("HLT_Ele27_WP80_PFMT50_v1")); }
      if ( triggerFired ("HLT_Ele27_WP80_PFMT50_v2"                     ) ) { fillVariableWithValue("H_27_WP80_PFMT50_2"         ,triggerPrescale("HLT_Ele27_WP80_PFMT50_v2"                     )); } else { fillVariableWithValue("H_27_WP80_PFMT50_2",-1 * triggerPrescale("HLT_Ele27_WP80_PFMT50_v2")); }
      if ( triggerFired ("HLT_Ele27_WP80_PFMT50_v3"                     ) ) { fillVariableWithValue("H_27_WP80_PFMT50_3"         ,triggerPrescale("HLT_Ele27_WP80_PFMT50_v3"                     )); } else { fillVariableWithValue("H_27_WP80_PFMT50_3",-1 * triggerPrescale("HLT_Ele27_WP80_PFMT50_v3")); }
      if ( triggerFired ("HLT_Ele27_WP80_PFMT50_v4"                     ) ) { fillVariableWithValue("H_27_WP80_PFMT50_4"         ,triggerPrescale("HLT_Ele27_WP80_PFMT50_v4"                     )); } else { fillVariableWithValue("H_27_WP80_PFMT50_4",-1 * triggerPrescale("HLT_Ele27_WP80_PFMT50_v4")); }

      if ( triggerFired ("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v1") ) { fillVariableWithValue ("H_17_8_CIdT_TIdVL_CIsVL_TIsVL_1",triggerPrescale("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v1")); } else { fillVariableWithValue("H_17_8_CIdT_TIdVL_CIsVL_TIsVL_1",-1 * triggerPrescale("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v1")); }	
      if ( triggerFired ("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2") ) { fillVariableWithValue ("H_17_8_CIdT_TIdVL_CIsVL_TIsVL_2",triggerPrescale("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2")); } else { fillVariableWithValue("H_17_8_CIdT_TIdVL_CIsVL_TIsVL_2",-1 * triggerPrescale("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2")); }	
      if ( triggerFired ("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3") ) { fillVariableWithValue ("H_17_8_CIdT_TIdVL_CIsVL_TIsVL_3",triggerPrescale("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3")); } else { fillVariableWithValue("H_17_8_CIdT_TIdVL_CIsVL_TIsVL_3",-1 * triggerPrescale("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3")); }	
      if ( triggerFired ("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4") ) { fillVariableWithValue ("H_17_8_CIdT_TIdVL_CIsVL_TIsVL_4",triggerPrescale("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4")); } else { fillVariableWithValue("H_17_8_CIdT_TIdVL_CIsVL_TIsVL_4",-1 * triggerPrescale("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4")); }	
      if ( triggerFired ("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v5") ) { fillVariableWithValue ("H_17_8_CIdT_TIdVL_CIsVL_TIsVL_5",triggerPrescale("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v5")); } else { fillVariableWithValue("H_17_8_CIdT_TIdVL_CIsVL_TIsVL_5",-1 * triggerPrescale("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v5")); }	
      			 											      			                                      																		                                  				 											  
      if ( triggerFired ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5") ) { fillVariableWithValue ("H_17_8_CIdT_CIsVL_TIdVL_TIsVL_5",triggerPrescale("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5")); } else { fillVariableWithValue("H_17_8_CIdT_CIsVL_TIdVL_TIsVL_5",-1 * triggerPrescale("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5")); }		
      if ( triggerFired ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6") ) { fillVariableWithValue ("H_17_8_CIdT_CIsVL_TIdVL_TIsVL_6",triggerPrescale("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6")); } else { fillVariableWithValue("H_17_8_CIdT_CIsVL_TIdVL_TIsVL_6",-1 * triggerPrescale("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6")); }		
      if ( triggerFired ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7") ) { fillVariableWithValue ("H_17_8_CIdT_CIsVL_TIdVL_TIsVL_7",triggerPrescale("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7")); } else { fillVariableWithValue("H_17_8_CIdT_CIsVL_TIdVL_TIsVL_7",-1 * triggerPrescale("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7")); }		
      if ( triggerFired ("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8") ) { fillVariableWithValue ("H_17_8_CIdT_CIsVL_TIdVL_TIsVL_8",triggerPrescale("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8")); } else { fillVariableWithValue("H_17_8_CIdT_CIsVL_TIdVL_TIsVL_8",-1 * triggerPrescale("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8")); }		
      
      if ( triggerFired ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1") ) { fillVariableWithValue ("H_17_8_CIdL_CIsVL_1",triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1")); } else { fillVariableWithValue("H_17_8_CIdL_CIsVL_1",-1 * triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1")); }		
      if ( triggerFired ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2") ) { fillVariableWithValue ("H_17_8_CIdL_CIsVL_2",triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2")); } else { fillVariableWithValue("H_17_8_CIdL_CIsVL_2",-1 * triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2")); }		
      if ( triggerFired ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3") ) { fillVariableWithValue ("H_17_8_CIdL_CIsVL_3",triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3")); } else { fillVariableWithValue("H_17_8_CIdL_CIsVL_3",-1 * triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3")); }		
      if ( triggerFired ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4") ) { fillVariableWithValue ("H_17_8_CIdL_CIsVL_4",triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4")); } else { fillVariableWithValue("H_17_8_CIdL_CIsVL_4",-1 * triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4")); }		
      if ( triggerFired ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5") ) { fillVariableWithValue ("H_17_8_CIdL_CIsVL_5",triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5")); } else { fillVariableWithValue("H_17_8_CIdL_CIsVL_5",-1 * triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5")); }		
      if ( triggerFired ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6") ) { fillVariableWithValue ("H_17_8_CIdL_CIsVL_6",triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6")); } else { fillVariableWithValue("H_17_8_CIdL_CIsVL_6",-1 * triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6")); }			
      if ( triggerFired ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v7") ) { fillVariableWithValue ("H_17_8_CIdL_CIsVL_7",triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v7")); } else { fillVariableWithValue("H_17_8_CIdL_CIsVL_7",-1 * triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v7")); }			
      if ( triggerFired ("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v8") ) { fillVariableWithValue ("H_17_8_CIdL_CIsVL_8",triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v8")); } else { fillVariableWithValue("H_17_8_CIdL_CIsVL_8",-1 * triggerPrescale("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v8")); }			
      
      if ( triggerFired ("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_v1") ) { fillVariableWithValue ("H_32_17_CIdT_CIsT_TIdT_TIsT_1",triggerPrescale("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_v1")); } else { fillVariableWithValue("H_32_17_CIdT_CIsT_TIdT_TIsT_1",-1 * triggerPrescale("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_v1")); }
    }
      
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

    //Event filters at RECO level
    fillVariableWithValue( "PassBeamScraping", !isBeamScraping ) ;
    fillVariableWithValue( "PassPrimaryVertex", isPrimaryVertex ) ;
    fillVariableWithValue( "PassHBHENoiseFilter", passHBHENoiseFilter ) ;
    fillVariableWithValue( "PassBeamHaloFilterLoose", passBeamHaloFilterLoose ) ;
    fillVariableWithValue( "PassBeamHaloFilterTight", passBeamHaloFilterTight ) ;
    fillVariableWithValue( "PassTrackingFailure", !isTrackingFailure ) ;
    fillVariableWithValue( "PassCaloBoundaryDRFilter", passCaloBoundaryDRFilter ) ;
    fillVariableWithValue( "PassEcalMaskedCellDRFilter", passEcalMaskedCellDRFilter ) ;
    
    fillVariableWithValue( "nPileUpInt_BXminus1", -1 );
    fillVariableWithValue( "nPileUpInt_BX0"     , -1 );
    fillVariableWithValue( "nPileUpInt_BXplus1" , -1  );

    if ( isData == 0 ){
      for(int pu=0; pu<PileUpInteractions->size(); pu++) {
	if(PileUpOriginBX->at(pu) == -1 ) fillVariableWithValue( "nPileUpInt_BXminus1", PileUpInteractions->at(pu) ) ;	      
	if(PileUpOriginBX->at(pu) == 0  ) fillVariableWithValue( "nPileUpInt_BX0"     , PileUpInteractions->at(pu) ) ;	      
	if(PileUpOriginBX->at(pu) == 1  ) fillVariableWithValue( "nPileUpInt_BXplus1", PileUpInteractions->at(pu) ) ;	      	    
      }
    }
    
    // nVertex and pile-up
    fillVariableWithValue( "nVertex", VertexChi2->size() ) ;
    fillVariableWithValue( "nVertex_good", v_idx_vertex_good.size() ) ;
    
    // nEle
    fillVariableWithValue( "nEle_Ana"  , v_idx_ele_PtCut_IDISO_ANA.size() ) ;
    fillVariableWithValue( "nEle_Stored", v_idx_ele_PtCut_IDISO_STORE.size() ) ;

    // nJet
    fillVariableWithValue( "nJet_Ana", v_idx_jet_PtCut_noOverlap_ID_ANA.size() ) ;
    fillVariableWithValue( "nJet_Stored", v_idx_jet_PtCut_noOverlap_ID_STORE.size() ) ;
    fillVariableWithValue( "nJet_btagTCHE_Ana", v_idx_jet_PtCut_noOverlap_ID_TCHEL_ANA.size() ) ;

    // nMuon
    fillVariableWithValue( "nMuon_Ana", v_idx_muon_PtCut_IDISO_ANA.size() ) ;
    fillVariableWithValue( "nMuon_Stored", v_idx_muon_PtCut_IDISO_STORE.size() ) ;

    // MET
    fillVariableWithValue("MET_Pt", thisMET);
    fillVariableWithValue("MET_Phi", thisMETPhi);
    fillVariableWithValue("GenMET_Pt", thisGenMET);
    fillVariableWithValue("GenMET_Phi", thisGenMETPhi);  
    fillVariableWithValue("PFMETSig"         , (*PFMETSig         )[0] );
    fillVariableWithValue("PFMETSigMatrixDXX", (*PFMETSigMatrixDXX)[0] ) ;
    fillVariableWithValue("PFMETSigMatrixDXY", (*PFMETSigMatrixDXY)[0] ) ;
    fillVariableWithValue("PFMETSigMatrixDYX", (*PFMETSigMatrixDYX)[0] ) ;
    fillVariableWithValue("PFMETSigMatrixDYY", (*PFMETSigMatrixDYY)[0] ) ;
 
    // SUMET
    fillVariableWithValue("SumET", thisSumET);
    fillVariableWithValue("GenSumET", thisGenSumET);

     // 1st ele 
    double MT_METEle1, mDeltaPhiMETEle1 = -999;
    if( v_idx_ele_PtCut_IDISO_STORE.size() >= 1 )
      {
	fillVariableWithValue( "Ele1_Pt"       , ElectronPt                       -> at( v_idx_ele_PtCut_IDISO_STORE[0]) );
	fillVariableWithValue( "Ele1_Energy"   , ElectronEnergy                   -> at( v_idx_ele_PtCut_IDISO_STORE[0]) );
	fillVariableWithValue( "Ele1_Eta"      , ElectronEta                      -> at( v_idx_ele_PtCut_IDISO_STORE[0]) );
	fillVariableWithValue( "Ele1_Phi"      , ElectronPhi                      -> at( v_idx_ele_PtCut_IDISO_STORE[0]) );
	fillVariableWithValue( "Ele1_Charge"   , ElectronCharge                   -> at( v_idx_ele_PtCut_IDISO_STORE[0]) );
	fillVariableWithValue( "Ele1_Dist"     , ElectronDist                     -> at( v_idx_ele_PtCut_IDISO_STORE[0]) );
	fillVariableWithValue( "Ele1_DCotTheta", ElectronDCotTheta                -> at( v_idx_ele_PtCut_IDISO_STORE[0]) );
	fillVariableWithValue( "Ele1_VtxD0"    , ElectronVtxDistXY                -> at( v_idx_ele_PtCut_IDISO_STORE[0]) );
	fillVariableWithValue( "Ele1_ValidFrac", ElectronTrackValidFractionOfHits -> at( v_idx_ele_PtCut_IDISO_STORE[0]) );
	// DeltaPhi - MET vs 1st ele
	TVector2 v_MET;
	TVector2 v_ele;
	v_MET.SetMagPhi( thisMET , thisMETPhi);
	v_ele.SetMagPhi( ElectronPt->at(v_idx_ele_PtCut_IDISO_STORE[0]) , ElectronPhi->at(v_idx_ele_PtCut_IDISO_STORE[0]) );
	float deltaphi = v_MET.DeltaPhi(v_ele);
	fillVariableWithValue( "mDPhi_METEle1", fabs(deltaphi) );
        mDeltaPhiMETEle1 = fabs(deltaphi);

	// transverse mass enu
	MT_METEle1 = sqrt(2 * ElectronPt->at(v_idx_ele_PtCut_IDISO_STORE[0]) * thisMET * (1 - cos(deltaphi)) );
	fillVariableWithValue("MT_Ele1MET", MT_METEle1);

	//PT(e,nu)
	TVector2 v_ele_MET;
	v_ele_MET = v_ele + v_MET;
	fillVariableWithValue("Pt_Ele1MET", v_ele_MET.Mod());
      }

    // 2nd ele 
    double MT_METEle2, mDeltaPhiMETEle2 = -999;
    if( v_idx_ele_PtCut_IDISO_STORE.size() >= 2 )
      {
	fillVariableWithValue( "Ele2_Pt"       , ElectronPt                       -> at( v_idx_ele_PtCut_IDISO_STORE[1]) );
	fillVariableWithValue( "Ele2_Energy"   , ElectronEnergy                   -> at( v_idx_ele_PtCut_IDISO_STORE[1]) );
	fillVariableWithValue( "Ele2_Eta"      , ElectronEta                      -> at( v_idx_ele_PtCut_IDISO_STORE[1]) );
	fillVariableWithValue( "Ele2_Phi"      , ElectronPhi                      -> at( v_idx_ele_PtCut_IDISO_STORE[1]) );
	fillVariableWithValue( "Ele2_Charge"   , ElectronCharge                   -> at( v_idx_ele_PtCut_IDISO_STORE[1]) );
	fillVariableWithValue( "Ele2_Dist"     , ElectronDist                     -> at( v_idx_ele_PtCut_IDISO_STORE[1]) );
	fillVariableWithValue( "Ele2_DCotTheta", ElectronDCotTheta                -> at( v_idx_ele_PtCut_IDISO_STORE[1]) );
	fillVariableWithValue( "Ele2_VtxD0"    , ElectronVtxDistXY                -> at( v_idx_ele_PtCut_IDISO_STORE[1]) );
	fillVariableWithValue( "Ele2_ValidFrac", ElectronTrackValidFractionOfHits -> at( v_idx_ele_PtCut_IDISO_STORE[1]) );

	// DeltaPhi - MET vs 2nd ele
	TVector2 v_MET;
	TVector2 v_ele;
	v_MET.SetMagPhi( thisMET , thisMETPhi);
	v_ele.SetMagPhi( ElectronPt->at(v_idx_ele_PtCut_IDISO_STORE[1]) , ElectronPhi->at(v_idx_ele_PtCut_IDISO_STORE[1]) );
	float deltaphi = v_MET.DeltaPhi(v_ele);
	fillVariableWithValue( "mDPhi_METEle2", fabs(deltaphi) );
        mDeltaPhiMETEle2 = fabs(deltaphi);

	// transverse mass enu
	//MT_METEle2 = sqrt(2 * ElectronPt->at(v_idx_ele_PtCut_IDISO[1]) * thisMET * (1 - cos(deltaphi)) );
	//fillVariableWithValue("MT_Ele2MET", MT_METEle2);

	//PT(e,nu)
	//TVector2 v_ele_MET;
	//v_ele_MET = v_ele + v_MET;
	//fillVariableWithValue("Pt_Ele2MET", v_ele_MET.Mod());
      }

    // 1st jet 
    double mDeltaPhiMET1stJet = -999;
    if( v_idx_jet_PtCut_noOverlap_ID_STORE.size() >= 1 )        
      {
	fillVariableWithValue( "Jet1_Pt"      , JetPt    ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]) );
	fillVariableWithValue( "Jet1_Energy"  , JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]) );
	fillVariableWithValue( "Jet1_Eta"     , JetEta   ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]) );
	fillVariableWithValue( "Jet1_Phi"     , JetPhi   ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]) );
        fillVariableWithValue( "Jet1_btagTCHE", JetTCHE  ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]) );

	//DeltaPhi - MET vs 1st jet
	TVector2 v_MET;
	TVector2 v_jet;
	v_MET.SetMagPhi( 1 , thisMETPhi );
	v_jet.SetMagPhi( 1 , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]) );
	float deltaphi = v_MET.DeltaPhi(v_jet);
	fillVariableWithValue( "mDPhi_METJet1", fabs(deltaphi) );
        mDeltaPhiMET1stJet = fabs(deltaphi);
      }
    
    // 2nd jet 
    double mDeltaPhiMET2ndJet = -999;
    if( v_idx_jet_PtCut_noOverlap_ID_STORE.size() >= 2 )        
      {
	fillVariableWithValue( "Jet2_Pt"      , JetPt    ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]) );
	fillVariableWithValue( "Jet2_Energy"  , JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]) );
	fillVariableWithValue( "Jet2_Eta"     , JetEta   ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]) );
	fillVariableWithValue( "Jet2_Phi"     , JetPhi   ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]) );
        fillVariableWithValue( "Jet2_btagTCHE", JetTCHE  ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]) );

	//DeltaPhi - MET vs 2nd jet
	TVector2 v_MET;
	TVector2 v_jet;
	v_MET.SetMagPhi( 1 , thisMETPhi );
	v_jet.SetMagPhi( 1 , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]) );
	float deltaphi = v_MET.DeltaPhi(v_jet);
	fillVariableWithValue( "mDPhi_METJet2", fabs(deltaphi) );
        mDeltaPhiMET2ndJet = fabs(deltaphi);
      }

    // 3rd jet 
    double mDeltaPhiMET3rdJet = -999;
    if( v_idx_jet_PtCut_noOverlap_ID_STORE.size() >= 3 )        
      {
	fillVariableWithValue( "Jet3_Pt"      , JetPt     -> at(v_idx_jet_PtCut_noOverlap_ID_STORE[2]) );
	fillVariableWithValue( "Jet3_Energy"  , JetEnergy -> at(v_idx_jet_PtCut_noOverlap_ID_STORE[2]) );
	fillVariableWithValue( "Jet3_Eta"     , JetEta    -> at(v_idx_jet_PtCut_noOverlap_ID_STORE[2]) );
	fillVariableWithValue( "Jet3_Phi"     , JetPhi    -> at(v_idx_jet_PtCut_noOverlap_ID_STORE[2]) );
        fillVariableWithValue( "Jet3_btagTCHE", JetTCHE   -> at(v_idx_jet_PtCut_noOverlap_ID_STORE[2]) );

	//DeltaPhi - MET vs 3rd jet
	TVector2 v_MET;
	TVector2 v_jet;
	v_MET.SetMagPhi( 1 , thisMETPhi );
	v_jet.SetMagPhi( 1 , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID_STORE[2]) );
	float deltaphi = v_MET.DeltaPhi(v_jet);
	fillVariableWithValue( "mDPhi_METJet3", fabs(deltaphi) );
        mDeltaPhiMET3rdJet = fabs(deltaphi);
      }

    // CaloJets for HLT
    
    fillVariableWithValue ( "nCaloJet_Stored",  v_idx_calojet_PtCut_noOverlap_ID_STORE.size() );
    fillVariableWithValue ( "nCaloJet_Ana"   ,  v_idx_calojet_PtCut_noOverlap_ID_ANA.size() );
    
    if ( v_idx_calojet_PtCut_noOverlap_ID_STORE.size() >= 1 ) {
      fillVariableWithValue ( "CaloJet1_Pt" , caloJetPt  -> at ( v_idx_calojet_PtCut_noOverlap_ID_STORE[0] ));
      fillVariableWithValue ( "CaloJet1_Eta", caloJetEta -> at ( v_idx_calojet_PtCut_noOverlap_ID_STORE[0] ));
      fillVariableWithValue ( "CaloJet1_Phi", caloJetPhi -> at ( v_idx_calojet_PtCut_noOverlap_ID_STORE[0] ));
    }

    if ( v_idx_calojet_PtCut_noOverlap_ID_STORE.size() >= 2 ) {
      fillVariableWithValue ( "CaloJet2_Pt" , caloJetPt  -> at ( v_idx_calojet_PtCut_noOverlap_ID_STORE[1] ));
      fillVariableWithValue ( "CaloJet2_Eta", caloJetEta -> at ( v_idx_calojet_PtCut_noOverlap_ID_STORE[1] ));
      fillVariableWithValue ( "CaloJet2_Phi", caloJetPhi -> at ( v_idx_calojet_PtCut_noOverlap_ID_STORE[1] ));
    }

    if ( v_idx_calojet_PtCut_noOverlap_ID_STORE.size() >= 3 ) {
      fillVariableWithValue ( "CaloJet3_Pt" , caloJetPt  -> at ( v_idx_calojet_PtCut_noOverlap_ID_STORE[2] ));
      fillVariableWithValue ( "CaloJet3_Eta", caloJetEta -> at ( v_idx_calojet_PtCut_noOverlap_ID_STORE[2] ));
      fillVariableWithValue ( "CaloJet3_Phi", caloJetPhi -> at ( v_idx_calojet_PtCut_noOverlap_ID_STORE[2] ));
    }

    // 1st muon 
    double MT_METMuon1, mDeltaPhiMETMuon1 = -999;
    if( v_idx_muon_PtCut_IDISO_STORE.size() >= 1 )
      {
	fillVariableWithValue( "Muon1_Pt"    , MuonPt     ->at(v_idx_muon_PtCut_IDISO_STORE[0]) );
	fillVariableWithValue( "Muon1_Energy", MuonEnergy ->at(v_idx_muon_PtCut_IDISO_STORE[0]) );
	fillVariableWithValue( "Muon1_Eta"   , MuonEta    ->at(v_idx_muon_PtCut_IDISO_STORE[0]) );
	fillVariableWithValue( "Muon1_Phi"   , MuonPhi    ->at(v_idx_muon_PtCut_IDISO_STORE[0]) );
	fillVariableWithValue( "Muon1_Charge", MuonCharge ->at(v_idx_muon_PtCut_IDISO_STORE[0]) );
      }

    // define booleans

    //-----------------------------------------------------------------
    // Fill your multi-object variables with values
    //-----------------------------------------------------------------

    bool OneEle=false;
    bool TwoEles=false;
    bool OneJet=false;
    bool TwoJets=false;
    bool ThreeJets=false;
    if( v_idx_ele_PtCut_IDISO_STORE.size()        >= 1 ) OneEle = true;
    if( v_idx_ele_PtCut_IDISO_STORE.size()        >= 2 ) TwoEles = true;
    if( v_idx_jet_PtCut_noOverlap_ID_STORE.size() >= 1 ) OneJet = true;
    if( v_idx_jet_PtCut_noOverlap_ID_STORE.size() >= 2 ) TwoJets = true;
    if( v_idx_jet_PtCut_noOverlap_ID_STORE.size() >= 3 ) ThreeJets = true;

    if ( TwoEles ) {
      TLorentzVector ele1;
      ele1.SetPtEtaPhiM(ElectronPt ->at(v_idx_ele_PtCut_IDISO_STORE[0]),
			ElectronEta->at(v_idx_ele_PtCut_IDISO_STORE[0]),
			ElectronPhi->at(v_idx_ele_PtCut_IDISO_STORE[0]),0);

      TLorentzVector ele2;
      ele2.SetPtEtaPhiM(ElectronPt ->at(v_idx_ele_PtCut_IDISO_STORE[1]),
			ElectronEta->at(v_idx_ele_PtCut_IDISO_STORE[1]),
			ElectronPhi->at(v_idx_ele_PtCut_IDISO_STORE[1]),0);

      TLorentzVector ee = ele1 + ele2;
      
      fillVariableWithValue ("M_e1e2" , ee.M ()  );
      fillVariableWithValue ("Pt_e1e2", ee.Pt () );
    }

    //Delta R

    if ( (OneEle) && (OneJet) )
      {	

	TLorentzVector ele1;
	ele1.SetPtEtaPhiM(ElectronPt ->at(v_idx_ele_PtCut_IDISO_STORE[0]),
			  ElectronEta->at(v_idx_ele_PtCut_IDISO_STORE[0]),
			  ElectronPhi->at(v_idx_ele_PtCut_IDISO_STORE[0]),0);
	TLorentzVector jet1;
	jet1.SetPtEtaPhiE(JetPt ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]),
			  JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]) );

	TLorentzVector e1j1 = ele1 + jet1;

	fillVariableWithValue("M_e1j1", e1j1.M() );
	fillVariableWithValue("DR_Ele1Jet1", ele1.DeltaR(jet1));

	if ( (TwoJets) )
	  {
	    TLorentzVector jet2;
	    jet2.SetPtEtaPhiE(JetPt    ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]),
			      JetEta   ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]),
			      JetPhi   ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]),
			      JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]) );
	    
	    TLorentzVector ej2 = ele1 + jet2;

	    fillVariableWithValue("M_e1j2", ej2.M() );
	    fillVariableWithValue("DR_Ele1Jet2", ele1.DeltaR(jet2));	   
	  }

	if( (TwoEles) )
	  {
	    TLorentzVector ele2;
	    ele2.SetPtEtaPhiM(ElectronPt ->at(v_idx_ele_PtCut_IDISO_STORE[1]),
			      ElectronEta->at(v_idx_ele_PtCut_IDISO_STORE[1]),
			      ElectronPhi->at(v_idx_ele_PtCut_IDISO_STORE[1]),0);	    

	    TLorentzVector e2j1 = ele2 + jet1;

	    fillVariableWithValue("DR_Ele2Jet1", ele2.DeltaR(jet1));	
	    fillVariableWithValue("M_e2j1", e2j1.M() );   	    
	  }

	if( (TwoEles) && (TwoJets) )
	  {
	    TLorentzVector ele2;
	    ele2.SetPtEtaPhiM(ElectronPt ->at(v_idx_ele_PtCut_IDISO_STORE[1]),
			      ElectronEta->at(v_idx_ele_PtCut_IDISO_STORE[1]),
			      ElectronPhi->at(v_idx_ele_PtCut_IDISO_STORE[1]),0);	    
	    TLorentzVector jet2;
	    jet2.SetPtEtaPhiE(JetPt    ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]),
			      JetEta   ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]),
			      JetPhi   ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]),
			      JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]) );

	    TLorentzVector e2j2 = ele2 + jet2;

	    fillVariableWithValue("DR_Ele2Jet2", ele2.DeltaR(jet2));	   	    
	    fillVariableWithValue("M_e2j2", e2j2.M() );   	    
	    
	  }       	
      }
    
    // Mjj
    if (TwoJets)
      {
	TLorentzVector jet1, jet2, jj;
	jet1.SetPtEtaPhiE(JetPt     ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]),
			  JetEta    ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]),
			  JetPhi    ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]),
			  JetEnergy ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]) );
	jet2.SetPtEtaPhiE(JetPt     ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]),
			  JetEta    ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]),
			  JetPhi    ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]),
			  JetEnergy ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]) );
	jj = jet1+jet2;

	fillVariableWithValue("M_j1j2", jj.M());
	fillVariableWithValue("DR_Jet1Jet2", jet1.DeltaR(jet2));
	fillVariableWithValue("Pt_j1j2", jj.Pt());
      }

    if (ThreeJets)
      {
	TLorentzVector jet1, jet2, jet3, j1j3, j2j3;
	jet1.SetPtEtaPhiE(JetPt    ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]),
			  JetEta   ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]),
			  JetPhi   ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]),
			  JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]) );
	jet2.SetPtEtaPhiE(JetPt    ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]),
			  JetEta   ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]),
			  JetPhi   ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]),
			  JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]) );
	jet3.SetPtEtaPhiE(JetPt    ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[2]),
			  JetEta   ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[2]),
			  JetPhi   ->at(v_idx_jet_PtCut_noOverlap_ID_STORE[2]),
			  JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID_STORE[2]) );
	j1j3 = jet1+jet3;
	j2j3 = jet2+jet3;

	fillVariableWithValue("M_j1j3", j1j3.M());
	fillVariableWithValue("M_j2j3", j2j3.M());
      }

    // ST enujj
    if ( (OneEle) && (TwoJets) )
      {
	double sT_enujj =
	  ElectronPt->at(v_idx_ele_PtCut_IDISO_STORE[0]) +
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]) +
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]) +
	  thisMET;
	fillVariableWithValue("sT_enujj", sT_enujj);
      }

    if ( (TwoEles) && (TwoJets) )
      {
	double sT_eejj =
	  ElectronPt->at(v_idx_ele_PtCut_IDISO_STORE[0]) +
	  ElectronPt->at(v_idx_ele_PtCut_IDISO_STORE[1]) +
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID_STORE[0]) +
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID_STORE[1]) ;
	fillVariableWithValue("sT_eejj", sT_eejj);
      }

    // Evaluate cuts (but do not apply them)
    evaluateCuts();
    

    if ( reducedSkimType == 1 ) { 
      // Produce skim and reduced skim
      if( passedAllPreviousCuts("PassHBHENoiseFilter") 
	  && ( getVariableValue("nEle_Ana") == 1 ) 
	  && passedCut("Ele1_Pt")
	  && passedCut("MET_Pt")
	  && passedCut("Jet1_Pt")
	  && passedCut("Jet2_Pt")	
	  ) {
	fillSkimTree();
	fillReducedSkimTree();
      }
    }

    if ( reducedSkimType == 2 ) { 
      if( passedAllPreviousCuts("PassHBHENoiseFilter") 
	  && ( getVariableValue("nEle_Ana") == 2 ) 
	  && passedCut("Ele1_Pt")
	  && passedCut("Ele2_Pt")
	  && passedCut("Jet1_Pt")
	  && passedCut("Jet2_Pt")	
	  ) {
	fillSkimTree();
	fillReducedSkimTree();
      }
    }
  } // End of loop over events

  /*//------------------------------------------------------------------
   *
   *
   *      
   *      End analysis loop!
   *
   *
   *
   *///-----------------------------------------------------------------

  delete randomNumGen;

  STDOUT("analysisClass::Loop() ends");
}
