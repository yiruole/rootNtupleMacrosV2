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
  // Make histograms
  //-----------------------------------------------------------------    

  CreateUserTH1D("Ele1_Pt_Pass" , 100,0,1000);
  CreateUserTH1D("Ele1_Pt_Total", 100,0,1000);
  CreateUserTH1D("RecoOverGenPt", 100 , 0, 2.0 );
  CreateUserTH1D("RecoGenDR"    , 100 , 0, 0.25);
  CreateUserTH2D("RecoGenDR_vs_RecoOverGenPt", 100 , 0, 0.25, 100 , 0, 2.0);
  
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
  
  double ele_PtCut       = getPreCutValue1("ele_PtCut");
  
  double eleEta_bar      = getPreCutValue1("eleEta_bar");
  double eleEta_end_min  = getPreCutValue1("eleEta_end");
  double eleEta_end_max  = getPreCutValue2("eleEta_end");
  
  double GenRecoDeltaRMax   = getPreCutValue1("GenRecoDeltaRMax");
  double GenRecoPtTolerance = getPreCutValue1("GenRecoPtTolerance");

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
  //Long64_t nentries = 10;
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

    if(jetAlgorithm==2) //Calo jets
      {
	for (int ijet=0 ; ijet < CaloJetPt->size() ; ijet++)
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
    // Selection: GEN Electrons
    //-----------------------------------------------------------------
    
    std::vector<int> v_idx_gen_ele; 

    for (int igen = 0; igen < GenParticlePt->size(); ++igen ) {
      int status = (*GenParticleStatus)[igen];
      int pdgid  = (*GenParticlePdgId) [igen];
      
      if ( status     != 3 ) continue;
      if ( abs(pdgid) != 11) continue;

      v_idx_gen_ele.push_back ( igen );
      
    }
    
    if ( v_idx_gen_ele.size() != 1 ) {
      std::cout << "ERROR: There should only be one GEN electron!!!" << std::endl;
      exit(0);
    }

    //-----------------------------------------------------------------
    // Selection: RECO Electrons
    //-----------------------------------------------------------------    

    vector<int> v_idx_ele_PtCut_IDISO;
    vector<int> v_idx_ele_IDISO;	       
    vector<int> v_idx_ele;
    
    //Loop over electrons
    for(int iele=0; iele<ElectronPt->size(); iele++){
      
      int passEleSel = 0;
      int isBarrel = 0;
      int isEndcap = 0;

      if( fabs( ElectronSCEta->at(iele) ) < eleEta_bar )       isBarrel = 1;
      if( fabs( ElectronSCEta->at(iele) ) > eleEta_end_min &&
	  fabs( ElectronSCEta->at(iele) ) < eleEta_end_max )   isEndcap = 1;

      if (isBarrel || isEndcap) v_idx_ele.push_back ( iele ) ;

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
	if ( ElectronPt -> at (iele) >= ele_PtCut ) v_idx_ele_PtCut_IDISO.push_back ( iele ) ;
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
    for (int iele=0; iele<v_idx_ele_PtCut_IDISO.size(); iele++)
      {
	TLorentzVector ele;
        ele.SetPtEtaPhiM(ElectronPt ->at(v_idx_ele_PtCut_IDISO[iele]),
			 ElectronEta->at(v_idx_ele_PtCut_IDISO[iele]),
			 ElectronPhi->at(v_idx_ele_PtCut_IDISO[iele]),0);
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
    // Fill your event-object variables with values
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

    for(int pu=0; pu<PileUpInteractions->size(); pu++) {
      if(PileUpOriginBX->at(pu) == -1 ) fillVariableWithValue( "nPileUpInt_BXminus1", PileUpInteractions->at(pu) ) ;	      
      if(PileUpOriginBX->at(pu) == 0  ) fillVariableWithValue( "nPileUpInt_BX0"     , PileUpInteractions->at(pu) ) ;	      
      if(PileUpOriginBX->at(pu) == 1  ) fillVariableWithValue( "nPileUpInt_BXplus1", PileUpInteractions->at(pu) ) ;	      	    
    }
    
    // nVertex and pile-up
    fillVariableWithValue( "nVertex", VertexChi2->size() ) ;
    fillVariableWithValue( "nVertex_good", v_idx_vertex_good.size() ) ;
    
    //-----------------------------------------------------------------
    // Is the GEN electron even pointing at the ECAL?
    //-----------------------------------------------------------------

    int matched = 0;
    double deltaR = -999;
    double ptRatio = -999;

    TLorentzVector gen_ele, reco_ele;
    gen_ele.SetPtEtaPhiM ( GenParticlePt  ->at(v_idx_gen_ele[0]),
			   GenParticleEta ->at(v_idx_gen_ele[0]),
			   GenParticlePhi ->at(v_idx_gen_ele[0]), 0.0 );

    int gen_isBarrel = 0;
    int gen_isEndcap = 0;
    int gen_is_fiducial  = 0;

    if( fabs( gen_ele.Eta()  ) < eleEta_bar )       gen_isBarrel = 1;
    if( fabs( gen_ele.Eta()  ) > eleEta_end_min &&
	fabs( gen_ele.Eta()  ) < eleEta_end_max )   gen_isEndcap = 1;

    if ( gen_isBarrel + gen_isEndcap > 0 ) gen_is_fiducial = 1;
    fillVariableWithValue ("GenEleIsFiducial", gen_is_fiducial );
    
    //-----------------------------------------------------------------
    // Do GEN-RECO matching
    //-----------------------------------------------------------------
    
    if ( v_idx_ele.size() > 0 ) {
      reco_ele.SetPtEtaPhiM ( ElectronPt ->at(v_idx_ele[0]),
			      ElectronEta->at(v_idx_ele[0]),
			      ElectronPhi->at(v_idx_ele[0]), 0.0 );
      deltaR  = gen_ele.DeltaR ( reco_ele ) ;
      ptRatio = reco_ele.Pt() / gen_ele.Pt();
      if ( deltaR < GenRecoDeltaRMax && fabs ( ptRatio - 1.0) < GenRecoPtTolerance ) 
	matched = 1;
    }

    fillVariableWithValue("MatchedGenToReco", matched ) ;

    //-----------------------------------------------------------------
    // Does the matched RECO electron pass the ID?
    //-----------------------------------------------------------------

    int pass_id = 0;

    if ( matched == 1 && v_idx_ele_PtCut_IDISO.size() >= 1 ) { 
      if ( v_idx_ele[0] == v_idx_ele_PtCut_IDISO[0] ) pass_id = 1;
    }
    
    fillVariableWithValue("RecoEleHasID",pass_id ) ;

    //-----------------------------------------------------------------
    // Evaluate cuts
    //-----------------------------------------------------------------
    
    // Evaluate cuts (but do not apply them)
    evaluateCuts();
    
    if ( passedCut("MatchedGenToReco") )
      FillUserTH1D("Ele1_Pt_Total", gen_ele.Pt() ) ;
    
    if ( passedCut("RecoEleHasID") ) {
      FillUserTH1D("Ele1_Pt_Pass" , gen_ele.Pt() ) ;
      FillUserTH1D("RecoOverGenPt", ptRatio      ) ; 
      FillUserTH1D("RecoGenDR"    , deltaR       ) ;	
      FillUserTH2D("RecoGenDR_vs_RecoOverGenPt", deltaR, ptRatio);
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
