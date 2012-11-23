#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

#include "Collection.h"
#include "GenParticle.h"
#include "Electron.h"
#include "Muon.h"
#include "PFJet.h"
#include "GenJet.h"
#include "HLTFilterObject.h"
#include "HLTFilterObjectCollectionHelper.h"

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

void analysisClass::Loop(){

  std::cout << "analysisClass::Loop() begins" <<std::endl;   
  if (fChain == 0) return;
  
  //--------------------------------------------------------------------------
  // Verbose? Or not?
  //--------------------------------------------------------------------------
  
  bool verbose = !true;
  
  //--------------------------------------------------------------------------
  // Decide which plots to save (default is to save everything)
  //--------------------------------------------------------------------------
  
  fillSkim                         ( !true  ) ;
  fillAllPreviousCuts              ( !true  ) ;
  fillAllOtherCuts                 ( !true  ) ;
  fillAllSameLevelAndLowerLevelCuts( !true  ) ;
  fillAllCuts                      ( !true  ) ;

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
  // What reduced skim type?
  // - 0: QCD   (loose electron)
  // - 1: enujj (HEEP electron) 
  // - 2: eejj  (HEEP electron)
  //--------------------------------------------------------------------------

  int reducedSkimType = getPreCutValue1("reducedSkimType");

  //--------------------------------------------------------------------------
  // Cuts for physics objects selection
  //--------------------------------------------------------------------------

  // jet cuts
  double jet_PtCut               = getPreCutValue1("jet_PtCut");
  double jet_EtaCut              = getPreCutValue1("jet_EtaCut");
  double jet_ele_DeltaRCut       = getPreCutValue1("jet_ele_DeltaRCut") ;
  double jet_muon_DeltaRCut      = getPreCutValue1("jet_muon_DeltaRCut");
  double jet_hltMatch_DeltaRCut  = getPreCutValue1("jet_hltMatch_DeltaRCut");

  // muon cuts
  double muon_EtaCut             = getPreCutValue1("muon_EtaCut");
  double muon_PtCut              = getPreCutValue1("muon_PtCut");
  double muon_hltMatch_DeltaRCut = getPreCutValue1("muon_hltMatch_DeltaRCut");

  // electron cuts
  double ele_PtCut    	         = getPreCutValue1("ele_PtCut");
  double ele_hltMatch_DeltaRCut  = getPreCutValue1("ele_hltMatch_DeltaRCut");

  //--------------------------------------------------------------------------
  // Tell the user how many entries we'll look at
  //--------------------------------------------------------------------------

  Long64_t nentries = fChain->GetEntries();
  // Long64_t nentries = 20000;
  std::cout << "analysisClass::Loop(): nentries = " << fChain -> GetEntries() << std::endl;   
  
  //--------------------------------------------------------------------------
  // Create HLT collections in advance (won't need all of them)
  //--------------------------------------------------------------------------
  
  // QCD photon filters
  CollectionPtr c_hltPhoton30_CaloIdVL_QCD_all;
  CollectionPtr c_hltPhoton50_CaloIdVL_QCD_all;
  CollectionPtr c_hltPhoton75_CaloIdVL_QCD_all;
  CollectionPtr c_hltPhoton90_CaloIdVL_QCD_all;
  CollectionPtr c_hltPhoton135_QCD_all;	      
  CollectionPtr c_hltPhoton150_QCD_all;	      
  CollectionPtr c_hltPhoton160_QCD_all;	      

  // TTbar muon filter
  CollectionPtr c_hltMuon22_TTbar_all;
  CollectionPtr c_hltPhoton22_TTbar_all;

  // Signal
  CollectionPtr c_hltEle30_Signal_all;
  CollectionPtr c_hltPFJet25_Signal_all;
  CollectionPtr c_hltPFJetNoPU25_Signal_all;
  CollectionPtr c_hltPFJet100_Signal_all;
  CollectionPtr c_hltPFJetNoPU100_Signal_all;
  
  /*//------------------------------------------------------------------
   *
   *
   *      
   *      Start analysis loop!
   *
   *
   *
   *///-----------------------------------------------------------------
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    //-----------------------------------------------------------------
    // Print progress
    //-----------------------------------------------------------------
    
    if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
    
    //-----------------------------------------------------------------
    // Get access to HLT decisions
    //-----------------------------------------------------------------

    getTriggers ( HLTKey, HLTInsideDatasetTriggerNames, HLTInsideDatasetTriggerDecisions,  HLTInsideDatasetTriggerPrescales ) ; 
        
    //-----------------------------------------------------------------
    // Get access to HLT filter objects
    //-----------------------------------------------------------------
    
    HLTFilterObjectCollectionHelper helper (*this);

    if ( reducedSkimType == 0 ){ 
      
      // QCD photon triggers
      c_hltPhoton30_CaloIdVL_QCD_all = helper.GetHLTFilterObjects("hltEG30CaloIdVLHEFilter");
      c_hltPhoton50_CaloIdVL_QCD_all = helper.GetHLTFilterObjects("hltPhoton50CaloIdVLHEFilter");
      c_hltPhoton75_CaloIdVL_QCD_all = helper.GetHLTFilterObjects("hltPhoton75CaloIdVLHEFilter");
      c_hltPhoton90_CaloIdVL_QCD_all = helper.GetHLTFilterObjects("hltPhoton90CaloIdVLHEFilter");
      c_hltPhoton135_QCD_all	     = helper.GetHLTFilterObjects("hltPhoton135HEFilter");  
      c_hltPhoton150_QCD_all	     = helper.GetHLTFilterObjects("hltPhoton150HEFilter");  
      c_hltPhoton160_QCD_all	     = helper.GetHLTFilterObjects("hltPhoton160HEFilter");  
      
    }

    else if ( reducedSkimType == 1 || reducedSkimType == 2 ){
      
      // TTbar trigger
      c_hltMuon22_TTbar_all        = helper.GetHLTFilterObjects("hltEG22EtFilterL1Mu3p5EG12");
      c_hltPhoton22_TTbar_all      = helper.GetHLTFilterObjects("hltMu22Photon22CaloIdLHEFilter");
      
      // Ele+jets signal triggers
      c_hltEle30_Signal_all        = helper.GetHLTFilterObjects("hltEle30CaloIdVTTrkIdTDphiFilter");
      c_hltPFJet25_Signal_all      = helper.GetHLTFilterObjects("hltEle30CaloIdVTTrkIdTDiCentralPFJet25EleCleaned");
      c_hltPFJetNoPU25_Signal_all  = helper.GetHLTFilterObjects("hltEle30CaloIdVTTrkIdTDiCentralPFNoPUJet25EleCleaned");
      c_hltPFJet100_Signal_all     = helper.GetHLTFilterObjects("hltEle30CaloIdVTTrkIdTDiCentralPFJet100EleCleaned");
      c_hltPFJetNoPU100_Signal_all = helper.GetHLTFilterObjects("hltEle30CaloIdVTTrkIdTDiCentralPFNoPUJet100EleCleaned");
      
    }

    //-----------------------------------------------------------------
    // Define initial, inclusive collections for physics objects
    //-----------------------------------------------------------------

    CollectionPtr c_gen_all   ( new Collection(*this, GenParticlePt -> size()));
    CollectionPtr c_genJet_all( new Collection(*this, GenJetPt      -> size()));
    CollectionPtr c_ele_all   ( new Collection(*this, ElectronPt    -> size()));
    CollectionPtr c_muon_all  ( new Collection(*this, MuonPt        -> size()));
    CollectionPtr c_pfjet_all ( new Collection(*this, PFJetPt       -> size()));
    
    //-----------------------------------------------------------------
    // All skims need GEN particles/jets
    //-----------------------------------------------------------------

    CollectionPtr c_genEle_final = c_gen_all    -> SkimByID<GenParticle>(GEN_ELE_HARD_SCATTER);
    CollectionPtr c_genJet_final = c_genJet_all;
    
    //-----------------------------------------------------------------
    // QCD skims    (reducedSkimType = 0     ) have loose electrons
    // Signal skims (reducedSkimType = 1 or 2) have HEEP  electrons
    //-----------------------------------------------------------------

    CollectionPtr c_ele_final;
    CollectionPtr c_ele_final_ptCut;
    
    if ( reducedSkimType == 0 ){ 
      CollectionPtr c_ele_loose = c_ele_all   -> SkimByID  <Electron> ( HEEP_LOOSE );
      c_ele_final               = c_ele_loose;
      c_ele_final_ptCut         = c_ele_final -> SkimByMinPt<Electron>( ele_PtCut  );
    }

    else if ( reducedSkimType == 1 || reducedSkimType == 2 ){
      CollectionPtr c_ele_HEEP  = c_ele_all -> SkimByID <Electron> ( HEEP );
      c_ele_final               = c_ele_HEEP;
      c_ele_final_ptCut         = c_ele_final -> SkimByMinPt<Electron>( ele_PtCut  );
    }

    //-----------------------------------------------------------------
    // All skims need muons
    //-----------------------------------------------------------------
    
    CollectionPtr c_muon_eta2p1            = c_muon_all       -> SkimByEtaRange<Muon> ( -muon_EtaCut, muon_EtaCut );
    CollectionPtr c_muon_eta2p1_ID         = c_muon_eta2p1    -> SkimByID      <Muon> ( MUON_TIGHT_PFISO04 );
    CollectionPtr c_muon_final             = c_muon_eta2p1_ID;
    CollectionPtr c_muon_final_ptCut       = c_muon_final     -> SkimByMinPt   <Muon> ( muon_PtCut );
    
    //-----------------------------------------------------------------
    // All skims need PFJets
    //-----------------------------------------------------------------
    
    CollectionPtr c_pfjet_central                     = c_pfjet_all                          -> SkimByEtaRange   <PFJet>          ( -jet_EtaCut, jet_EtaCut   );
    CollectionPtr c_pfjet_central_ID                  = c_pfjet_central                      -> SkimByID         <PFJet>          ( PFJET_LOOSE );    
    CollectionPtr c_pfjet_central_ID_noMuonOverlap    = c_pfjet_central_ID                   -> SkimByVetoDRMatch<PFJet, Muon>    ( c_muon_final_ptCut   , jet_ele_DeltaRCut  );
    CollectionPtr c_pfjet_central_ID_noLeptonOverlap  = c_pfjet_central_ID_noMuonOverlap     -> SkimByVetoDRMatch<PFJet, Electron>( c_ele_final_ptCut    , jet_muon_DeltaRCut );
    CollectionPtr c_pfjet_final                       = c_pfjet_central_ID_noLeptonOverlap;
    CollectionPtr c_pfjet_final_ptCut                 = c_pfjet_final                        -> SkimByMinPt      <PFJet>          ( jet_PtCut );
    
    //-----------------------------------------------------------------
    // Get ready to fill variables 
    //-----------------------------------------------------------------

    resetCuts();

    //-----------------------------------------------------------------
    // Fill your single-object variables with values
    //-----------------------------------------------------------------
    
    fillVariableWithValue( "isData"   , isData     );
    fillVariableWithValue( "bunch"    , bunch      );
    fillVariableWithValue( "event"    , event      );
    fillVariableWithValue( "ls"       , ls         );
    fillVariableWithValue( "orbit"    , orbit      );
    fillVariableWithValue( "run"      , run        );
    fillVariableWithValue( "ProcessID", ProcessID  );
    fillVariableWithValue( "PtHat"    , PtHat      );
    fillVariableWithValue( "Weight"   , Weight     );

    //-----------------------------------------------------------------
    // Fill MET filter values
    //-----------------------------------------------------------------
    
    fillVariableWithValue("PassJSON"                   , passJSON(run, ls, isData)                       );    
    fillVariableWithValue("PassBPTX0"                  , isData == 1 ? int(isBPTX0        == 1) : 1      );
    fillVariableWithValue("PassPhysDecl"               , isData == 1 ? int(isPhysDeclared == 1) : 1      );
    fillVariableWithValue("PassBeamScraping"           , int(isBeamScraping                         == 0));
    fillVariableWithValue("PassPrimaryVertex"          , int(isPrimaryVertex                        == 1));
    fillVariableWithValue("PassBeamHaloFilterLoose"    , int(passBeamHaloFilterLoose                == 1));
    fillVariableWithValue("PassBeamHaloFilterTight"    , int(passBeamHaloFilterTight                == 1));
    fillVariableWithValue("PassHBHENoiseFilter"        , int(passHBHENoiseFilter                    == 1));
    fillVariableWithValue("PassHcalLaserEventFilter"   , int(passHcalLaserEventFilter               == 0));
    fillVariableWithValue("PassEcalDeadCellBoundEnergy", int(passEcalDeadCellBoundaryEnergyFilter   == 0));
    fillVariableWithValue("PassEcalDeadCellTrigPrim"   , int(passEcalDeadCellTriggerPrimitiveFilter == 0));
    fillVariableWithValue("PassTrackingFailure"        , int(passTrackingFailureFilter              == 0));
    fillVariableWithValue("PassBadEESupercrystalFilter", int(passBadEESupercrystalFilter            == 0));
    fillVariableWithValue("PassEcalLaserCorrFilter"    , int(passEcalLaserCorrFilter                == 0));

    //-----------------------------------------------------------------
    // Fill MET values
    //-----------------------------------------------------------------
    
    fillVariableWithValue("PFMET_Pt"            , PFMET               -> at (0));      
    fillVariableWithValue("PFMET_Phi"		, PFMETPhi	      -> at (0));
    fillVariableWithValue("PFMET_Type1_Pt"      , PFMETType1Cor       -> at (0));      
    fillVariableWithValue("PFMET_Type1_Phi" 	, PFMETPhiType1Cor    -> at (0));
    fillVariableWithValue("PFMET_Type01_Pt"     , PFMETType01Cor      -> at (0));      
    fillVariableWithValue("PFMET_Type01_Phi"	, PFMETPhiType01Cor   -> at (0));
    fillVariableWithValue("PFMET_Type01XY_Pt"   , PFMETType01XYCor    -> at (0));      
    fillVariableWithValue("PFMET_Type01XY_Phi"	, PFMETPhiType01XYCor -> at (0));
    
    if ( isData == 0 ) { 
      if ( reducedSkimType != 0 ){ 
	fillVariableWithValue("GenMET_Pt"		, GenMETTrue	      -> at (0));
	fillVariableWithValue("GenMET_Phi"       	, GenMETPhiTrue	      -> at (0));
      }
    }

    //-----------------------------------------------------------------
    // Fill pileup variables
    //-----------------------------------------------------------------
        
    fillVariableWithValue( "nPileUpInt_BXminus1", -1 );
    fillVariableWithValue( "nPileUpInt_BX0"     , -1 );
    fillVariableWithValue( "nPileUpInt_BXplus1" , -1 );
    fillVariableWithValue( "nVertex", VertexChi2->size() ) ;
    
    if ( isData == 0 ){
      for(int pu=0; pu<PileUpInteractions->size(); pu++) {
	if(PileUpOriginBX->at(pu) == 0  ) { 
	  fillVariableWithValue( "nPileUpInt_BX0" , PileUpInteractions    ->at(pu));
	  fillVariableWithValue( "nPileUpInt_True", PileUpInteractionsTrue->at(pu));
	}
	if(PileUpOriginBX->at(pu) == -1 ) fillVariableWithValue( "nPileUpInt_BXminus1", PileUpInteractions->at(pu));
	if(PileUpOriginBX->at(pu) == 1  ) fillVariableWithValue( "nPileUpInt_BXplus1" , PileUpInteractions->at(pu));
      }
    }

    //-----------------------------------------------------------------
    // How many ID'd objects are there?
    //-----------------------------------------------------------------

    int n_muon_store   = c_muon_final          -> GetSize();
    int n_ele_store    = c_ele_final           -> GetSize();
    int n_jet_store    = c_pfjet_final         -> GetSize();
    int n_genEle_store = c_genEle_final        -> GetSize();
    int n_genJet_store = c_genJet_final        -> GetSize();
    
    int n_muon_ptCut   = c_muon_final_ptCut    -> GetSize();
    int n_ele_ptCut    = c_ele_final_ptCut     -> GetSize();
    int n_jet_ptCut    = c_pfjet_final_ptCut   -> GetSize();

    //-----------------------------------------------------------------
    // All skims need GEN particles/jets
    //-----------------------------------------------------------------

    if ( reducedSkimType != 0 ) { 

      fillVariableWithValue("nGenJet_ptCut", n_genJet_store);
      fillVariableWithValue("nGenEle_ptCut", n_genEle_store);
      
      fillVariableWithValue("nGenJet_store", min(n_genJet_store,5));
      fillVariableWithValue("nGenEle_store", min(n_genEle_store,2));

      if ( n_genJet_store >= 1 ) { 
	GenJet genJet1 = c_genJet_final -> GetConstituent<GenJet>(0);
	fillVariableWithValue ( "GenJet1_Pt" , genJet1.Pt () );
	fillVariableWithValue ( "GenJet1_Eta", genJet1.Eta() );
	fillVariableWithValue ( "GenJet1_Phi", genJet1.Phi() );
	
	if ( n_genJet_store >= 2 ) { 
	  GenJet genJet2 = c_genJet_final -> GetConstituent<GenJet>(1);
	  fillVariableWithValue ( "GenJet2_Pt" , genJet2.Pt () );
	  fillVariableWithValue ( "GenJet2_Eta", genJet2.Eta() );
	  fillVariableWithValue ( "GenJet2_Phi", genJet2.Phi() );
	  
	  if ( n_genJet_store >= 3 ) { 
	    GenJet genJet3 = c_genJet_final -> GetConstituent<GenJet>(2);
	    fillVariableWithValue ( "GenJet3_Pt" , genJet3.Pt () );
	    fillVariableWithValue ( "GenJet3_Eta", genJet3.Eta() );
	    fillVariableWithValue ( "GenJet3_Phi", genJet3.Phi() );
	    
	    if ( n_genJet_store >= 4 ) { 
	      GenJet genJet4 = c_genJet_final -> GetConstituent<GenJet>(3);
	      fillVariableWithValue ( "GenJet4_Pt" , genJet4.Pt () );
	      fillVariableWithValue ( "GenJet4_Eta", genJet4.Eta() );
	      fillVariableWithValue ( "GenJet4_Phi", genJet4.Phi() );
	      
	      if ( n_genJet_store >= 5 ) { 
		GenJet genJet5 = c_genJet_final -> GetConstituent<GenJet>(4);
		fillVariableWithValue ( "GenJet5_Pt" , genJet5.Pt () );
		fillVariableWithValue ( "GenJet5_Eta", genJet5.Eta() );
		fillVariableWithValue ( "GenJet5_Phi", genJet5.Phi() );
	      }
	    }
	  }
	}
      }
      
      if ( n_genEle_store >= 1 ){ 
	GenParticle genEle1 = c_genEle_final -> GetConstituent<GenParticle>(0);
	fillVariableWithValue ( "GenEle1_Pt" , genEle1.Pt () );
	fillVariableWithValue ( "GenEle1_Eta", genEle1.Eta() );
	fillVariableWithValue ( "GenEle1_Phi", genEle1.Phi() );
	
	if ( n_genEle_store >= 2 ){ 
	  GenParticle genEle2 = c_genEle_final -> GetConstituent<GenParticle>(1);
	  fillVariableWithValue ( "GenEle2_Pt" , genEle2.Pt () );
	  fillVariableWithValue ( "GenEle2_Eta", genEle2.Eta() );
	  fillVariableWithValue ( "GenEle2_Phi", genEle2.Phi() );
	}
      }
    }

    //-----------------------------------------------------------------
    // All skims need muons
    //-----------------------------------------------------------------
    
    fillVariableWithValue ("nMuon_ptCut", n_muon_ptCut);
    fillVariableWithValue ("nMuon_store", min(n_muon_store,2));

    if ( n_muon_store >= 1 ){ 

      Muon muon1 = c_muon_final -> GetConstituent<Muon>(0);
      double hltMuon1Pt = triggerMatchPt<HLTFilterObject, Muon>(c_hltMuon22_TTbar_all, muon1, muon_hltMatch_DeltaRCut);
      fillVariableWithValue ("Muon1_Pt"       , muon1.Pt    ());
      fillVariableWithValue ("Muon1_Eta"      , muon1.Eta   ());
      fillVariableWithValue ("Muon1_Phi"      , muon1.Phi   ());
      fillVariableWithValue ("Muon1_Charge"   , muon1.Charge());
      fillVariableWithValue ("Muon1_hltMuonPt", hltMuon1Pt    );
      
      if ( n_muon_store >= 2 ){ 

	Muon muon2 = c_muon_final -> GetConstituent<Muon>(1);
	double hltMuon2Pt = triggerMatchPt<HLTFilterObject, Muon>(c_hltMuon22_TTbar_all, muon2, muon_hltMatch_DeltaRCut);
	fillVariableWithValue ("Muon2_Pt"       , muon2.Pt    ());
	fillVariableWithValue ("Muon2_Eta"      , muon2.Eta   ());
	fillVariableWithValue ("Muon2_Phi"      , muon2.Phi   ());
	fillVariableWithValue ("Muon2_Charge"   , muon2.Charge());
	fillVariableWithValue ("Muon2_hltMuonPt", hltMuon2Pt    );

      }
    }

    //-----------------------------------------------------------------
    // Fill variables for QCD skim (reducedSkimType == 0 )
    //-----------------------------------------------------------------

    if ( reducedSkimType == 0 ) { 

      fillVariableWithValue ("nLooseEle_store"   , min(n_ele_store,3));
      fillVariableWithValue ("nJetLooseEle_store", min(n_jet_store,5));
      fillVariableWithValue ("nLooseEle_ptCut"   , n_ele_ptCut );
      fillVariableWithValue ("nJetLooseEle_ptCut", n_jet_ptCut );

      int n_filters = ( c_hltPhoton30_CaloIdVL_QCD_all -> GetSize() + 
			c_hltPhoton50_CaloIdVL_QCD_all -> GetSize() + 
			c_hltPhoton75_CaloIdVL_QCD_all -> GetSize() + 
			c_hltPhoton90_CaloIdVL_QCD_all -> GetSize() + 
			c_hltPhoton135_QCD_all         -> GetSize() + 
			c_hltPhoton150_QCD_all         -> GetSize() + 
			c_hltPhoton160_QCD_all         -> GetSize() );
      
      if ( n_ele_store >= 1 ){
	Electron loose_ele1 = c_ele_final -> GetConstituent<Electron>(0);

	double hltPhotonPt = -999.;
	if ( n_filters != 0 ) {
	  double hltPhotonPt_array [7] =  { triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton30_CaloIdVL_QCD_all, loose_ele1, ele_hltMatch_DeltaRCut),
					    triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton50_CaloIdVL_QCD_all, loose_ele1, ele_hltMatch_DeltaRCut), 
					    triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton75_CaloIdVL_QCD_all, loose_ele1, ele_hltMatch_DeltaRCut), 
					    triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton90_CaloIdVL_QCD_all, loose_ele1, ele_hltMatch_DeltaRCut), 
					    triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton135_QCD_all        , loose_ele1, ele_hltMatch_DeltaRCut), 
					    triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton150_QCD_all        , loose_ele1, ele_hltMatch_DeltaRCut), 
					    triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton160_QCD_all        , loose_ele1, ele_hltMatch_DeltaRCut) };
	  for (int iFilter = 0; iFilter < 7; ++iFilter){
	    if ( hltPhotonPt_array[iFilter] > 0.0 ) hltPhotonPt = hltPhotonPt_array[iFilter];
	  }
	}
	
	fillVariableWithValue( "LooseEle1_PassID"       , loose_ele1.PassUserID ( HEEP ));
	fillVariableWithValue( "LooseEle1_Pt"           , loose_ele1.Pt()               );
	fillVariableWithValue( "LooseEle1_Energy"       , loose_ele1.CaloEnergy()       );
	fillVariableWithValue( "LooseEle1_Eta"          , loose_ele1.Eta()              );
	fillVariableWithValue( "LooseEle1_Phi"          , loose_ele1.Phi()              );
	fillVariableWithValue( "LooseEle1_Charge"       , loose_ele1.Charge()           );
	fillVariableWithValue( "LooseEle1_Dist"         , loose_ele1.Dist()             );
	fillVariableWithValue( "LooseEle1_DCotTheta"    , loose_ele1.DCotTheta()        );
	fillVariableWithValue( "LooseEle1_VtxD0"        , loose_ele1.LeadVtxDistXY()    );
	fillVariableWithValue( "LooseEle1_ValidFrac"    , loose_ele1.ValidFrac()        );
	fillVariableWithValue( "LooseEle1_MissingHits"  , loose_ele1.MissingHits()      );
	fillVariableWithValue( "LooseEle1_hltPhotonPt"  , hltPhotonPt );
	
	if ( n_ele_store >= 2 ){
	  
	  Electron loose_ele2 = c_ele_final -> GetConstituent<Electron>(1);
	  hltPhotonPt = -999.;
	  if ( n_filters != 0 ) {
	    double hltPhotonPt_array [7] =  { triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton30_CaloIdVL_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut),
					      triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton50_CaloIdVL_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut), 
					      triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton75_CaloIdVL_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut), 
					      triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton90_CaloIdVL_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut), 
					      triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton135_QCD_all        , loose_ele2, ele_hltMatch_DeltaRCut), 
					      triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton150_QCD_all        , loose_ele2, ele_hltMatch_DeltaRCut), 
					      triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton160_QCD_all        , loose_ele2, ele_hltMatch_DeltaRCut) };
	    for (int iFilter = 0; iFilter < 7; ++iFilter){
	      if ( hltPhotonPt_array[iFilter] > 0.0 ) hltPhotonPt = hltPhotonPt_array[iFilter];
	    }
	  }
	  
	  fillVariableWithValue( "LooseEle2_PassID"       , loose_ele2.PassUserID ( HEEP ));
	  fillVariableWithValue( "LooseEle2_Pt"           , loose_ele2.Pt()               );
	  fillVariableWithValue( "LooseEle2_Energy"       , loose_ele2.CaloEnergy()       );
	  fillVariableWithValue( "LooseEle2_Eta"          , loose_ele2.Eta()              );
	  fillVariableWithValue( "LooseEle2_Phi"          , loose_ele2.Phi()              );
	  fillVariableWithValue( "LooseEle2_Charge"       , loose_ele2.Charge()           );
	  fillVariableWithValue( "LooseEle2_Dist"         , loose_ele2.Dist()             );
	  fillVariableWithValue( "LooseEle2_DCotTheta"    , loose_ele2.DCotTheta()        );
	  fillVariableWithValue( "LooseEle2_VtxD0"        , loose_ele2.LeadVtxDistXY()    );
	  fillVariableWithValue( "LooseEle2_ValidFrac"    , loose_ele2.ValidFrac()        );
	  fillVariableWithValue( "LooseEle2_MissingHits"  , loose_ele2.MissingHits()      );
	  fillVariableWithValue( "LooseEle2_hltPhotonPt"  , hltPhotonPt );
	  
	  if ( n_ele_store >= 3 ){
	    Electron loose_ele3 = c_ele_final -> GetConstituent<Electron>(2);
	    hltPhotonPt = -999.;
	    if ( n_filters != 0 ) {
	      double hltPhotonPt_array [7] =  { triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton30_CaloIdVL_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut),
						triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton50_CaloIdVL_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut), 
						triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton75_CaloIdVL_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut), 
						triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton90_CaloIdVL_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut), 
						triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton135_QCD_all        , loose_ele2, ele_hltMatch_DeltaRCut), 
						triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton150_QCD_all        , loose_ele2, ele_hltMatch_DeltaRCut), 
						triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton160_QCD_all        , loose_ele2, ele_hltMatch_DeltaRCut) };
	      for (int iFilter = 0; iFilter < 7; ++iFilter){
		if ( hltPhotonPt_array[iFilter] > 0.0 ) hltPhotonPt = hltPhotonPt_array[iFilter];
	      }
	    }
	    fillVariableWithValue( "LooseEle3_PassID"       , loose_ele3.PassUserID ( HEEP ));
	    fillVariableWithValue( "LooseEle3_Pt"           , loose_ele3.Pt()               );
	    fillVariableWithValue( "LooseEle3_Energy"       , loose_ele3.CaloEnergy()       );
	    fillVariableWithValue( "LooseEle3_Eta"          , loose_ele3.Eta()              );
	    fillVariableWithValue( "LooseEle3_Phi"          , loose_ele3.Phi()              );
	    fillVariableWithValue( "LooseEle3_Charge"       , loose_ele3.Charge()           );
	    fillVariableWithValue( "LooseEle3_Dist"         , loose_ele3.Dist()             );
	    fillVariableWithValue( "LooseEle3_DCotTheta"    , loose_ele3.DCotTheta()        );
	    fillVariableWithValue( "LooseEle3_VtxD0"        , loose_ele3.LeadVtxDistXY()    );
	    fillVariableWithValue( "LooseEle3_ValidFrac"    , loose_ele3.ValidFrac()        );
	    fillVariableWithValue( "LooseEle3_MissingHits"  , loose_ele3.MissingHits()      );
	    fillVariableWithValue( "LooseEle3_hltPhotonPt"  , hltPhotonPt );
	  }
	}
      }

      if ( n_jet_store >= 1 ){
	PFJet loose_jet1 = c_pfjet_final -> GetConstituent<PFJet>(0);
	fillVariableWithValue( "JetLooseEle1_Pt"      , loose_jet1.Pt()                         );
	fillVariableWithValue( "JetLooseEle1_JECUnc"  , loose_jet1.JECUnc()                     );
	fillVariableWithValue( "JetLooseEle1_Energy"  , loose_jet1.Energy()                     );
	fillVariableWithValue( "JetLooseEle1_Eta"     , loose_jet1.Eta()                        );
	fillVariableWithValue( "JetLooseEle1_Phi"     , loose_jet1.Phi()                        );
        fillVariableWithValue( "JetLooseEle1_btagCSV" , loose_jet1.CombinedSecondaryVertexBTag());
	
	if ( n_jet_store >= 2 ){
	  PFJet loose_jet2 = c_pfjet_final -> GetConstituent<PFJet>(1);
	  fillVariableWithValue( "JetLooseEle2_Pt"      , loose_jet2.Pt()                         );
	  fillVariableWithValue( "JetLooseEle2_JECUnc"  , loose_jet2.JECUnc()                     );
	  fillVariableWithValue( "JetLooseEle2_Energy"  , loose_jet2.Energy()                     );
	  fillVariableWithValue( "JetLooseEle2_Eta"     , loose_jet2.Eta()                        );
	  fillVariableWithValue( "JetLooseEle2_Phi"     , loose_jet2.Phi()                        );
	  fillVariableWithValue( "JetLooseEle2_btagCSV" , loose_jet2.CombinedSecondaryVertexBTag());

	  if ( n_jet_store >= 3 ){
	    PFJet loose_jet3 = c_pfjet_final -> GetConstituent<PFJet>(2);
	    fillVariableWithValue( "JetLooseEle3_Pt"      , loose_jet3.Pt()                         );
	    fillVariableWithValue( "JetLooseEle3_JECUnc"  , loose_jet3.JECUnc()                     );
	    fillVariableWithValue( "JetLooseEle3_Energy"  , loose_jet3.Energy()                     );
	    fillVariableWithValue( "JetLooseEle3_Eta"     , loose_jet3.Eta()                        );
	    fillVariableWithValue( "JetLooseEle3_Phi"     , loose_jet3.Phi()                        );
	    fillVariableWithValue( "JetLooseEle3_btagCSV" , loose_jet3.CombinedSecondaryVertexBTag());

	    if ( n_jet_store >= 4 ){
	      PFJet loose_jet4 = c_pfjet_final -> GetConstituent<PFJet>(3);
	      fillVariableWithValue( "JetLooseEle4_Pt"      , loose_jet4.Pt()                         );
	      fillVariableWithValue( "JetLooseEle4_JECUnc"  , loose_jet4.JECUnc()                     );
	      fillVariableWithValue( "JetLooseEle4_Energy"  , loose_jet4.Energy()                     );
	      fillVariableWithValue( "JetLooseEle4_Eta"     , loose_jet4.Eta()                        );
	      fillVariableWithValue( "JetLooseEle4_Phi"     , loose_jet4.Phi()                        );
	      fillVariableWithValue( "JetLooseEle4_btagCSV" , loose_jet4.CombinedSecondaryVertexBTag());

	      if ( n_jet_store >= 5 ){
		PFJet loose_jet5 = c_pfjet_final -> GetConstituent<PFJet>(4);
		fillVariableWithValue( "JetLooseEle5_Pt"      , loose_jet5.Pt()                         );
		fillVariableWithValue( "JetLooseEle5_JECUnc"  , loose_jet5.JECUnc()                     );
		fillVariableWithValue( "JetLooseEle5_Energy"  , loose_jet5.Energy()                     );
		fillVariableWithValue( "JetLooseEle5_Eta"     , loose_jet5.Eta()                        );
		fillVariableWithValue( "JetLooseEle5_Phi"     , loose_jet5.Phi()                        );
		fillVariableWithValue( "JetLooseEle5_btagCSV" , loose_jet5.CombinedSecondaryVertexBTag());
	      }
	    }
	  }
	}
      }
    }

    //-----------------------------------------------------------------
    // Fill variables for signal-like skim (reducedSkimType == 1 or 2 )
    //-----------------------------------------------------------------
    
    // Electrons

    else if ( reducedSkimType == 1 || reducedSkimType == 2 ) { 
      
      fillVariableWithValue ("nEle_store" , min(n_ele_store,2) );
      fillVariableWithValue ("nJet_store" , min(n_jet_store,5) );
      fillVariableWithValue ("nEle_ptCut" , n_ele_ptCut );
      fillVariableWithValue ("nJet_ptCut" , n_jet_ptCut );

      if ( n_ele_store >= 1 ){
	Electron ele1 = c_ele_final -> GetConstituent<Electron>(0);
	double hltEle1Pt_signal = triggerMatchPt<HLTFilterObject, Electron>(c_hltEle30_Signal_all  , ele1, ele_hltMatch_DeltaRCut);
	double hltEle1Pt_ttbar  = triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton22_TTbar_all, ele1, ele_hltMatch_DeltaRCut);
	fillVariableWithValue( "Ele1_Pt"            , ele1.Pt()               );
	fillVariableWithValue( "Ele1_Energy"        , ele1.CaloEnergy()       );
	fillVariableWithValue( "Ele1_Eta"           , ele1.Eta()              );
	fillVariableWithValue( "Ele1_Phi"           , ele1.Phi()              );
	fillVariableWithValue( "Ele1_Charge"        , ele1.Charge()           );
	fillVariableWithValue( "Ele1_Dist"          , ele1.Dist()             );
	fillVariableWithValue( "Ele1_DCotTheta"     , ele1.DCotTheta()        );
	fillVariableWithValue( "Ele1_VtxD0"         , ele1.LeadVtxDistXY()    );
	fillVariableWithValue( "Ele1_ValidFrac"     , ele1.ValidFrac()        );
	fillVariableWithValue( "Ele1_MissingHits"   , ele1.MissingHits()      );
	fillVariableWithValue( "Ele1_hltEleSignalPt", hltEle1Pt_signal );
	fillVariableWithValue( "Ele1_hltEleTTbarPt" , hltEle1Pt_ttbar  );

	if ( n_ele_store >= 2 ){
	  Electron ele2 = c_ele_final -> GetConstituent<Electron>(1);
	  double hltEle2Pt_signal = triggerMatchPt<HLTFilterObject, Electron>(c_hltEle30_Signal_all  , ele2, ele_hltMatch_DeltaRCut);
	  double hltEle2Pt_ttbar  = triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton22_TTbar_all, ele2, ele_hltMatch_DeltaRCut);
	  fillVariableWithValue( "Ele2_Pt"           , ele2.Pt()               );
	  fillVariableWithValue( "Ele2_Energy"       , ele2.CaloEnergy()       );
	  fillVariableWithValue( "Ele2_Eta"          , ele2.Eta()              );
	  fillVariableWithValue( "Ele2_Phi"          , ele2.Phi()              );
	  fillVariableWithValue( "Ele2_Charge"       , ele2.Charge()           );
	  fillVariableWithValue( "Ele2_Dist"         , ele2.Dist()             );
	  fillVariableWithValue( "Ele2_DCotTheta"    , ele2.DCotTheta()        );
	  fillVariableWithValue( "Ele2_VtxD0"        , ele2.LeadVtxDistXY()    );
	  fillVariableWithValue( "Ele2_ValidFrac"    , ele2.ValidFrac()        );
	  fillVariableWithValue( "Ele2_MissingHits"  , ele2.MissingHits()      );
	  fillVariableWithValue( "Ele2_hltEleSignalPt", hltEle2Pt_signal );
	  fillVariableWithValue( "Ele2_hltEleTTbarPt" , hltEle2Pt_ttbar  );
	
	}
      }

      // Jets
      
      if ( n_jet_store >= 1 ){

	PFJet jet1 = c_pfjet_final -> GetConstituent<PFJet>(0);
	double hltJet1Pt     = triggerMatchPt<HLTFilterObject, PFJet >( c_hltPFJet25_Signal_all     , jet1, jet_hltMatch_DeltaRCut);
	double hltNoPUJet1Pt = triggerMatchPt<HLTFilterObject, PFJet >( c_hltPFJetNoPU25_Signal_all , jet1, jet_hltMatch_DeltaRCut);
	fillVariableWithValue( "Jet1_Pt"      , jet1.Pt()                         );
	fillVariableWithValue( "Jet1_JECUnc"  , jet1.JECUnc()                     );
	fillVariableWithValue( "Jet1_Energy"  , jet1.Energy()                     );
	fillVariableWithValue( "Jet1_Eta"     , jet1.Eta()                        );
	fillVariableWithValue( "Jet1_Phi"     , jet1.Phi()                        );
        fillVariableWithValue( "Jet1_btagCSV" , jet1.CombinedSecondaryVertexBTag());
	fillVariableWithValue( "Jet1_hltJetPt"	  , hltJet1Pt     );
	fillVariableWithValue( "Jet1_hltNoPUJetPt", hltNoPUJet1Pt );
	
	if ( n_jet_store >= 2 ){
	  PFJet jet2 = c_pfjet_final -> GetConstituent<PFJet>(1);
	  double hltJet2Pt     = triggerMatchPt<HLTFilterObject, PFJet >( c_hltPFJet25_Signal_all     , jet2, jet_hltMatch_DeltaRCut);
	  double hltNoPUJet2Pt = triggerMatchPt<HLTFilterObject, PFJet >( c_hltPFJetNoPU25_Signal_all , jet2, jet_hltMatch_DeltaRCut);
	  fillVariableWithValue( "Jet2_Pt"      , jet2.Pt()                         );
	  fillVariableWithValue( "Jet2_JECUnc"  , jet2.JECUnc()                     );
	  fillVariableWithValue( "Jet2_Energy"  , jet2.Energy()                     );
	  fillVariableWithValue( "Jet2_Eta"     , jet2.Eta()                        );
	  fillVariableWithValue( "Jet2_Phi"     , jet2.Phi()                        );
	  fillVariableWithValue( "Jet2_btagCSV" , jet2.CombinedSecondaryVertexBTag());
	  fillVariableWithValue( "Jet2_hltJetPt"    , hltJet2Pt     );
	  fillVariableWithValue( "Jet2_hltNoPUJetPt", hltNoPUJet2Pt );

	  if ( n_jet_store >= 3 ){
	    PFJet jet3 = c_pfjet_final -> GetConstituent<PFJet>(2);
	    double hltJet3Pt     = triggerMatchPt<HLTFilterObject, PFJet >( c_hltPFJet25_Signal_all     , jet3, jet_hltMatch_DeltaRCut);
	    double hltNoPUJet3Pt = triggerMatchPt<HLTFilterObject, PFJet >( c_hltPFJetNoPU25_Signal_all , jet3, jet_hltMatch_DeltaRCut);
	    fillVariableWithValue( "Jet3_Pt"      , jet3.Pt()                         );
	    fillVariableWithValue( "Jet3_JECUnc"  , jet3.JECUnc()                     );
	    fillVariableWithValue( "Jet3_Energy"  , jet3.Energy()                     );
	    fillVariableWithValue( "Jet3_Eta"     , jet3.Eta()                        );
	    fillVariableWithValue( "Jet3_Phi"     , jet3.Phi()                        );
	    fillVariableWithValue( "Jet3_btagCSV" , jet3.CombinedSecondaryVertexBTag());
	    fillVariableWithValue( "Jet3_hltJetPt"    , hltJet3Pt     );
	    fillVariableWithValue( "Jet3_hltNoPUJetPt", hltNoPUJet3Pt );

	    if ( n_jet_store >= 4 ){
	      PFJet jet4 = c_pfjet_final -> GetConstituent<PFJet>(3);
	      double hltJet4Pt     = triggerMatchPt<HLTFilterObject, PFJet >( c_hltPFJet25_Signal_all     , jet4, jet_hltMatch_DeltaRCut);
	      double hltNoPUJet4Pt = triggerMatchPt<HLTFilterObject, PFJet >( c_hltPFJetNoPU25_Signal_all , jet4, jet_hltMatch_DeltaRCut);
	      fillVariableWithValue( "Jet4_Pt"      , jet4.Pt()                         );
	      fillVariableWithValue( "Jet4_JECUnc"  , jet4.JECUnc()                     );
	      fillVariableWithValue( "Jet4_Energy"  , jet4.Energy()                     );
	      fillVariableWithValue( "Jet4_Eta"     , jet4.Eta()                        );
	      fillVariableWithValue( "Jet4_Phi"     , jet4.Phi()                        );
	      fillVariableWithValue( "Jet4_btagCSV" , jet4.CombinedSecondaryVertexBTag());
	      fillVariableWithValue( "Jet4_hltJetPt"    , hltJet4Pt     );
	      fillVariableWithValue( "Jet4_hltNoPUJetPt", hltNoPUJet4Pt );

	      if ( n_jet_store >= 5 ){
		PFJet jet5 = c_pfjet_final -> GetConstituent<PFJet>(4);
		double hltJet5Pt     = triggerMatchPt<HLTFilterObject, PFJet >( c_hltPFJet25_Signal_all     , jet5, jet_hltMatch_DeltaRCut);
		double hltNoPUJet5Pt = triggerMatchPt<HLTFilterObject, PFJet >( c_hltPFJetNoPU25_Signal_all , jet5, jet_hltMatch_DeltaRCut);
		fillVariableWithValue( "Jet5_Pt"      , jet5.Pt()                         );
		fillVariableWithValue( "Jet5_JECUnc"  , jet5.JECUnc()                     );
		fillVariableWithValue( "Jet5_Energy"  , jet5.Energy()                     );
		fillVariableWithValue( "Jet5_Eta"     , jet5.Eta()                        );
		fillVariableWithValue( "Jet5_Phi"     , jet5.Phi()                        );
		fillVariableWithValue( "Jet5_btagCSV" , jet5.CombinedSecondaryVertexBTag());
		fillVariableWithValue( "Jet5_hltJetPt"    , hltJet5Pt     );
		fillVariableWithValue( "Jet5_hltNoPUJetPt", hltNoPUJet5Pt );
	      }
	    }
	  }
	}
      }
    }

    //-----------------------------------------------------------------
    // Fill variables that depend on more than one object
    // All skims need this
    //-----------------------------------------------------------------

    TLorentzVector t_ele1, t_ele2, t_jet1, t_jet2, t_jet3;
    TLorentzVector t_MET;
    
    t_MET.SetPtEtaPhiM( PFMETType01XYCor -> at (0), 0.0, PFMETPhiType01XYCor -> at (0), 0.0 );

    if ( n_jet_store >= 1 ){

      PFJet jet1 = c_pfjet_final -> GetConstituent<PFJet>(0);
      t_jet1.SetPtEtaPhiM ( jet1.Pt(), jet1.Eta(), jet1.Phi(), 0.0 );

      fillVariableWithValue ("mDPhi_METJet1", fabs( t_MET.DeltaPhi ( t_jet1 )));
      
      if ( n_jet_store >= 2 ){
	
	PFJet jet2 = c_pfjet_final -> GetConstituent<PFJet>(1);
	t_jet2.SetPtEtaPhiM ( jet2.Pt(), jet2.Eta(), jet2.Phi(), 0.0 );
	TLorentzVector t_jet1jet2 = t_jet1 + t_jet2;

	fillVariableWithValue ("M_j1j2" , t_jet1jet2.M ());
	fillVariableWithValue ("Pt_j1j2", t_jet1jet2.Pt());
	fillVariableWithValue ("mDPhi_METJet2", fabs( t_MET.DeltaPhi ( t_jet2 )));
	fillVariableWithValue ("DR_Jet1Jet2"  , t_jet1.DeltaR( t_jet2 ));
	
	if ( n_jet_store >= 3 ){
	  
	  PFJet jet3 = c_pfjet_final -> GetConstituent<PFJet>(2);
	  t_jet3.SetPtEtaPhiM ( jet3.Pt(), jet3.Eta(), jet3.Phi(), 0.0 );
	  TLorentzVector t_jet1jet3 = t_jet1 + t_jet3;
	  TLorentzVector t_jet2jet3 = t_jet2 + t_jet3;

	  fillVariableWithValue ("M_j1j3", t_jet1jet3.M());
	  fillVariableWithValue ("M_j2j3", t_jet2jet3.M());
	  fillVariableWithValue ("mDPhi_METJet3", fabs( t_MET.DeltaPhi ( t_jet3 )));
	  
	}
      }
    }
    
    if ( n_ele_store >= 1 ) { 
      Electron ele1 = c_ele_final -> GetConstituent<Electron>(0);
      t_ele1.SetPtEtaPhiM ( ele1.Pt(), ele1.Eta(), ele1.Phi(), 0.0 );
      if ( n_ele_store >= 2 ) {
	Electron ele2 = c_ele_final -> GetConstituent<Electron>(1);
	t_ele2.SetPtEtaPhiM ( ele2.Pt(), ele2.Eta(), ele2.Phi(), 0.0 );

	TLorentzVector t_ele1ele2 = t_ele1 + t_ele2;
	fillVariableWithValue ("M_e1e2" , t_ele1ele2.M ());
	fillVariableWithValue ("Pt_e1e2", t_ele1ele2.Pt());
      }
    } 

    if ( n_ele_store >= 1 ){
      
      TLorentzVector t_ele1MET = t_ele1 + t_MET;
      fillVariableWithValue("mDPhi_METEle1", fabs ( t_MET.DeltaPhi(t_ele1)));
      fillVariableWithValue("MT_Ele1MET"   , t_ele1MET.Mt());
      fillVariableWithValue("Pt_Ele1MET"   , t_ele1MET.Pt());

      if ( n_jet_store >= 1 ){ 
	
	TLorentzVector t_ele1jet1 = t_ele1 + t_jet1;
	fillVariableWithValue("DR_Ele1Jet1", t_ele1.DeltaR ( t_jet1 ));
	fillVariableWithValue("M_e1j1"     , t_ele1jet1.M());

	if ( n_jet_store >= 2 ){ 
	  
	  TLorentzVector t_ele1jet2 = t_ele1 + t_jet2;
	  fillVariableWithValue("DR_Ele1Jet2", t_ele1.DeltaR ( t_jet2 ));
	  fillVariableWithValue("M_e1j2"     , t_ele1jet2.M());
	  fillVariableWithValue("sT_enujj"   , t_ele1.Pt() + t_MET.Pt() + t_jet1.Pt() + t_jet2.Pt());
	}
      }
    }


    if ( n_ele_store >= 2 ){
      fillVariableWithValue("mDPhi_METEle2", fabs ( t_MET.DeltaPhi(t_ele2)));

      if ( n_jet_store >= 1 ){ 
	
	TLorentzVector t_ele2jet1 = t_ele2 + t_jet1;
	fillVariableWithValue("DR_Ele2Jet1", t_ele2.DeltaR ( t_jet1 ));
	fillVariableWithValue("M_e2j1"     , t_ele2jet1.M());

	if ( n_jet_store >= 2 ){ 
	  
	  TLorentzVector t_ele2jet2 = t_ele2 + t_jet2;
	  fillVariableWithValue("DR_Ele2Jet2", t_ele2.DeltaR ( t_jet2 ));
	  fillVariableWithValue("M_e2j2"     , t_ele2jet2.M());
	  fillVariableWithValue("sT_eejj"    , t_ele1.Pt() + t_ele2.Pt() + t_jet1.Pt() + t_jet2.Pt());
	  
	}
      }
    }
    
    //-----------------------------------------------------------------
    // QCD triggers
    // - H_Photon30_CIdVL
    // - H_Photon50_CIdVL
    // - H_Photon75_CIdVL
    // - H_Photon90_CIdVL
    // - H_Photon135    
    // - H_Photon150    
    // - H_Photon160    
    //-----------------------------------------------------------------

    if ( reducedSkimType == 0 ) { 
    
      if ( ! isData ){
	fillTriggerVariable ( "HLT_Photon30_CaloIdVL_v13", "H_Photon30_CIdVL" );
	fillTriggerVariable ( "HLT_Photon50_CaloIdVL_v9" , "H_Photon50_CIdVL" );
	fillTriggerVariable ( "HLT_Photon75_CaloIdVL_v12", "H_Photon75_CIdVL" );
	fillTriggerVariable ( "HLT_Photon90_CaloIdVL_v9" , "H_Photon90_CIdVL" );
	fillTriggerVariable ( "HLT_Photon135_v6"         , "H_Photon135"      );
	fillTriggerVariable ( "HLT_Photon150_v3"         , "H_Photon150"      );
	fillTriggerVariable ( "HLT_Photon160_v3"         , "H_Photon160"      );
      }
      
      else if ( run >= 190456 && run <= 190738 ) {
	fillTriggerVariable ( "HLT_Photon30_CaloIdVL_v11", "H_Photon30_CIdVL" );
	fillTriggerVariable ( "HLT_Photon50_CaloIdVL_v7" , "H_Photon50_CIdVL" );
	fillTriggerVariable ( "HLT_Photon75_CaloIdVL_v10", "H_Photon75_CIdVL" );
	fillTriggerVariable ( "HLT_Photon90_CaloIdVL_v7" , "H_Photon90_CIdVL" );
	fillTriggerVariable ( "HLT_Photon135_v4"         , "H_Photon135"      );
	fillTriggerVariable ( "HLT_Photon150_v1"         , "H_Photon150"      );
	fillTriggerVariable ( "HLT_Photon160_v1"         , "H_Photon160"      );
      }
      
      else if ( run >= 190782 && run <= 191419 ) {
	fillTriggerVariable ( "HLT_Photon30_CaloIdVL_v12", "H_Photon30_CIdVL" );
	fillTriggerVariable ( "HLT_Photon50_CaloIdVL_v8" , "H_Photon50_CIdVL" );
	fillTriggerVariable ( "HLT_Photon75_CaloIdVL_v11", "H_Photon75_CIdVL" );
	fillTriggerVariable ( "HLT_Photon90_CaloIdVL_v8" , "H_Photon90_CIdVL" );
	fillTriggerVariable ( "HLT_Photon135_v5"         , "H_Photon135"      ); 
	fillTriggerVariable ( "HLT_Photon150_v2"         , "H_Photon150"      );
	fillTriggerVariable ( "HLT_Photon160_v2"         , "H_Photon160"      );
      }
      
      else if ( run >= 191691 && run <= 196531 ) {
	fillTriggerVariable ( "HLT_Photon30_CaloIdVL_v13", "H_Photon30_CIdVL" );
	fillTriggerVariable ( "HLT_Photon50_CaloIdVL_v9" , "H_Photon50_CIdVL" );
	fillTriggerVariable ( "HLT_Photon75_CaloIdVL_v12", "H_Photon75_CIdVL" );
	fillTriggerVariable ( "HLT_Photon90_CaloIdVL_v9" , "H_Photon90_CIdVL" );
	fillTriggerVariable ( "HLT_Photon135_v6"         , "H_Photon135"      );
	fillTriggerVariable ( "HLT_Photon150_v3"         , "H_Photon150"      );
	fillTriggerVariable ( "HLT_Photon160_v3"         , "H_Photon160"      );
      }
      
      else if ( run >= 198022 ) {
	fillTriggerVariable ( "HLT_Photon30_CaloIdVL_v14", "H_Photon30_CIdVL" );
	fillTriggerVariable ( "HLT_Photon50_CaloIdVL_v10", "H_Photon50_CIdVL" );
	fillTriggerVariable ( "HLT_Photon75_CaloIdVL_v13", "H_Photon75_CIdVL" );
	fillTriggerVariable ( "HLT_Photon90_CaloIdVL_v10", "H_Photon90_CIdVL" );
	fillTriggerVariable ( "HLT_Photon135_v7"         , "H_Photon135"      );
	fillTriggerVariable ( "HLT_Photon150_v4"         , "H_Photon150"      );
	fillTriggerVariable ( "HLT_Photon160_v4"         , "H_Photon160"      );
      }
      
      bool pass_trigger = ( getVariableValue("H_Photon30_CIdVL") > 0 || 
			    getVariableValue("H_Photon50_CIdVL") > 0 || 
			    getVariableValue("H_Photon75_CIdVL") > 0 || 
			    getVariableValue("H_Photon90_CIdVL") > 0 || 
			    getVariableValue("H_Photon135"     ) > 0 || 
			    getVariableValue("H_Photon150"     ) > 0 || 
			    getVariableValue("H_Photon160"     ) > 0 );
      
      fillVariableWithValue ("PassTrigger", pass_trigger ? 1 : 0 );
      
    }

    //-----------------------------------------------------------------
    // Signal triggers:
    // - HLT_Ele30_CaloIdVT_TrkIdT_PFJet100_PFJet25
    // - HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25
    // 
    // MuEG triggers
    // - HLT_Mu22_Photon22_CaloIdL
    //-----------------------------------------------------------------
    
    else if ( reducedSkimType == 1 || reducedSkimType == 2 ) { 
      
      if      ( ! isData )                       fillTriggerVariable( "HLT_Ele30_CaloIdVT_TrkIdT_PFJet100_PFJet25_v5", "H_Ele30_PFJet100_25");
      else if ( run >= 190456 && run <= 190738 ) fillTriggerVariable( "HLT_Ele30_CaloIdVT_TrkIdT_PFJet100_PFJet25_v3", "H_Ele30_PFJet100_25");
      else if ( run >= 190782 && run <= 191419 ) fillTriggerVariable( "HLT_Ele30_CaloIdVT_TrkIdT_PFJet100_PFJet25_v4", "H_Ele30_PFJet100_25");
      else if ( run >= 191691 && run <= 196027 ) fillTriggerVariable( "HLT_Ele30_CaloIdVT_TrkIdT_PFJet100_PFJet25_v5", "H_Ele30_PFJet100_25");
      
      if      ( ! isData )                       fillTriggerVariable( "HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v5", "H_Ele30_PFNoPUJet100_25");
      else if ( run >= 191691 && run <= 194225 ) fillTriggerVariable( "HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v4", "H_Ele30_PFNoPUJet100_25");
      else if ( run >= 194270 && run <= 196531 ) fillTriggerVariable( "HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v5", "H_Ele30_PFNoPUJet100_25");
      else if ( run >= 198022 && run <= 199608 ) fillTriggerVariable( "HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v6", "H_Ele30_PFNoPUJet100_25");
      else if ( run >= 199698 && run <= 202504 ) fillTriggerVariable( "HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v7", "H_Ele30_PFNoPUJet100_25");
      else if ( run >= 202970 )                  fillTriggerVariable( "HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v8", "H_Ele30_PFNoPUJet100_25");
      
      if      ( ! isData )                       fillTriggerVariable( "HLT_Mu22_Photon22_CaloIdL_v5", "H_Mu22_Photon22_CIdL");
      else if ( run >= 190456 && run <= 190738 ) fillTriggerVariable( "HLT_Mu22_Photon22_CaloIdL_v3", "H_Mu22_Photon22_CIdL");
      else if ( run >= 190782 && run <= 191419 ) fillTriggerVariable( "HLT_Mu22_Photon22_CaloIdL_v4", "H_Mu22_Photon22_CIdL");
      else if ( run >= 191691 && run <= 196531 ) fillTriggerVariable( "HLT_Mu22_Photon22_CaloIdL_v5", "H_Mu22_Photon22_CIdL");
      else if ( run >= 198022 && run <= 199608 ) fillTriggerVariable( "HLT_Mu22_Photon22_CaloIdL_v6", "H_Mu22_Photon22_CIdL");
      else if ( run >= 199698 )                  fillTriggerVariable( "HLT_Mu22_Photon22_CaloIdL_v7", "H_Mu22_Photon22_CIdL");
    }

    //-----------------------------------------------------------------
    // Evaluate the cuts
    //-----------------------------------------------------------------    

    evaluateCuts();

    //-----------------------------------------------------------------    
    // Fill the trees
    //-----------------------------------------------------------------    

    // QCD fake rate calculation skim
    
    if ( reducedSkimType == 0 ) { 
      if(passedCut("PassTrigger"      ) &&  
	 passedCut("nLooseEle_ptCut"  ) && 
	 passedCut("LooseEle1_Pt"     ) ){
	fillSkimTree();
	fillReducedSkimTree();
      }
    }    
 
    // enujj analysis skim

    else if ( reducedSkimType == 1 ) { 
      if( passedCut("nEle_store"       ) && 
	  passedCut("Ele1_Pt"          ) && 
	  passedCut("PFMET_Type01XY_Pt") && 
	  passedCut("Jet1_Pt"          ) && 
	  passedCut("Jet2_Pt"          ) ) {
	fillSkimTree();
	fillReducedSkimTree();
      }
    }

    // eejj analysis skim

    else if ( reducedSkimType == 2 ) { 
      if( passedCut("nEle_store"       ) && 
	  passedCut("Ele1_Pt"          ) &&
	  passedCut("Ele2_Pt"          ) && 
	  passedCut("Jet1_Pt"          ) &&
	  passedCut("Jet2_Pt"          ) ) {
	fillSkimTree();
	fillReducedSkimTree();
      }
    }
  } 
  
  std::cout << "analysisClass::Loop() ends" <<std::endl;   
  
}
