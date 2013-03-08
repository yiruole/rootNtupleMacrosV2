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

//--------------------------------------------------------------------------
// Function for trigger matching 
//--------------------------------------------------------------------------

template < class Object1, class Object2 > 
double triggerMatchPt ( const CollectionPtr & collection, Object2 & target_object, double delta_r_cut ){
  double matched_pt = -999.0;
  if ( collection ) { 
    int size = collection -> GetSize();
    if ( size > 0 ){ 
      Object1 matched_object = collection -> GetClosestInDR <Object1, Object2> ( target_object );
      double dr = matched_object.DeltaR ( & target_object );
      if ( dr < delta_r_cut ) { 
	matched_pt = matched_object.Pt();
      }
    }
  }
  return matched_pt;
}

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
  // Long64_t nentries = 1000;
  std::cout << "analysisClass::Loop(): nentries = " << fChain -> GetEntries() << std::endl;   
  
  //--------------------------------------------------------------------------
  // Create HLT collections in advance (won't need all of them)
  //--------------------------------------------------------------------------
  
  // Signal
  CollectionPtr c_hltEle30_Signal_all;
  CollectionPtr c_hltPFJetNoPU25_Signal_all;
  CollectionPtr c_hltPFJetNoPU100_Signal_all;
  CollectionPtr c_hltDoubleEle_Signal_all;

  // Tag and probe
  CollectionPtr c_hltEle27WP80_all;

  //--------------------------------------------------------------------------
  // Create historgams
  //--------------------------------------------------------------------------

  CreateUserTH1D( "nEle_all"              , 6, -0.5, 5.5 );
  CreateUserTH1D( "nEle_pass_ejj"         , 6, -0.5, 5.5 );
  CreateUserTH1D( "nEle_pass_eejj_veto"   , 6, -0.5, 5.5 );
  CreateUserTH1D( "nEle_pass_eejj_noveto" , 6, -0.5, 5.5 );
  CreateUserTH1D( "nEle_pass_enujj_veto"  , 6, -0.5, 5.5 );
  CreateUserTH1D( "nEle_pass_enujj_noveto", 6, -0.5, 5.5 );
  
  CreateUserTH1D( "nEle_ele1Pass"         , 6, -0.5, 5.5 );
  CreateUserTH1D( "nEle_ele2Pass_eejj"    , 6, -0.5, 5.5 );
  CreateUserTH1D( "nEle_passVeto_enujj"   , 6, -0.5, 5.5 );
  CreateUserTH1D( "nEle_passVeto_eejj"    , 6, -0.5, 5.5 );
  CreateUserTH1D( "nJet_jet1Pass"         , 17, -0.5, 16.5 );
  CreateUserTH1D( "nJet_jet2Pass"         , 17, -0.5, 16.5 );  
  CreateUserTH1D( "MET_metPass_enujj"     , 200, 0, 2000 );

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
    
    c_hltEle30_Signal_all        = helper.GetHLTFilterObjects("hltEle30CaloIdVTTrkIdTDphiFilter");
    c_hltPFJetNoPU25_Signal_all  = helper.GetHLTFilterObjects("hltEle30CaloIdVTTrkIdTDiCentralPFNoPUJet25EleCleaned");
    c_hltPFJetNoPU100_Signal_all = helper.GetHLTFilterObjects("hltEle30CaloIdVTTrkIdTDiCentralPFNoPUJet100EleCleaned");
    c_hltDoubleEle_Signal_all    = helper.GetHLTFilterObjects("hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter");
    
    //-----------------------------------------------------------------
    // Define initial, inclusive collections for physics objects
    //-----------------------------------------------------------------
    
    CollectionPtr c_ele_all   ( new Collection(*this, ElectronPt -> size()));
    CollectionPtr c_muon_all  ( new Collection(*this, MuonPt     -> size()));
    CollectionPtr c_pfjet_all ( new Collection(*this, PFJetPt    -> size()));
    
    //-----------------------------------------------------------------
    // ID electrons
    //-----------------------------------------------------------------

    CollectionPtr c_ele_final;
    CollectionPtr c_ele_final_ptCut;
    
    CollectionPtr c_ele_HEEP  = c_ele_all -> SkimByID <Electron> ( HEEP );
    c_ele_final               = c_ele_HEEP;
    c_ele_final_ptCut         = c_ele_final -> SkimByMinPt<Electron>( ele_PtCut  );
    
    //-----------------------------------------------------------------
    // ID muons
    //-----------------------------------------------------------------
    
    CollectionPtr c_muon_eta2p1      = c_muon_all       -> SkimByEtaRange<Muon> ( -muon_EtaCut, muon_EtaCut );
    CollectionPtr c_muon_eta2p1_ID   = c_muon_eta2p1    -> SkimByID      <Muon> ( MUON_TIGHT_PFISO04 );
    CollectionPtr c_muon_final       = c_muon_eta2p1_ID;
    CollectionPtr c_muon_final_ptCut = c_muon_final     -> SkimByMinPt   <Muon> ( muon_PtCut );
    
    //-----------------------------------------------------------------
    // ID jets
    //-----------------------------------------------------------------
    
    CollectionPtr c_pfjet_central                     = c_pfjet_all                          -> SkimByEtaRange   <PFJet>          ( -jet_EtaCut, jet_EtaCut   );
    CollectionPtr c_pfjet_central_ID                  = c_pfjet_central                      -> SkimByID         <PFJet>          ( PFJET_LOOSE );    
    CollectionPtr c_pfjet_central_ID_noMuonOverlap    = c_pfjet_central_ID                   -> SkimByVetoDRMatch<PFJet, Muon>    ( c_muon_final_ptCut   , jet_ele_DeltaRCut  );
    CollectionPtr c_pfjet_central_ID_noLeptonOverlap  = c_pfjet_central_ID_noMuonOverlap     -> SkimByVetoDRMatch<PFJet, Electron>( c_ele_final_ptCut    , jet_muon_DeltaRCut );
    CollectionPtr c_pfjet_final                       = c_pfjet_central_ID_noLeptonOverlap;
    CollectionPtr c_pfjet_final_ptCut                 = c_pfjet_final                        -> SkimByMinPt      <PFJet>          ( jet_PtCut );
    
    //------------------------------------------------------------------
    // Map triggers to decisions and prescales
    //------------------------------------------------------------------

    getTriggers ( HLTKey, HLTInsideDatasetTriggerNames, HLTInsideDatasetTriggerDecisions,  HLTInsideDatasetTriggerPrescales ) ; 

    //------------------------------------------------------------------
    // Get JSON and trigger decisions
    //------------------------------------------------------------------
    
    bool passed_json = passJSON ( run, ls , isData );
    bool trigger_fired = triggerFired ( "HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v8" );
    
    //-----------------------------------------------------------------
    // How many ID'd objects are there?
    //-----------------------------------------------------------------
    
    int n_ele   = c_ele_final   -> GetSize();
    int n_jet   = c_pfjet_final -> GetSize();

    //-----------------------------------------------------------------
    // Get ready to fill variables 
    //-----------------------------------------------------------------

    resetCuts();
    
    //------------------------------------------------------------------
    // Fill cut values
    //------------------------------------------------------------------

    fillVariableWithValue ( "PassJSON"    , int ( passed_json    ) );
    fillVariableWithValue ( "PassTrigger" , int ( trigger_fired  ) );
        
    //-----------------------------------------------------------------
    // Evaluate the cuts
    //-----------------------------------------------------------------    
    
    evaluateCuts();

    //-----------------------------------------------------------------    
    // Did we pass the minimum?
    //-----------------------------------------------------------------    
    
    bool pass_minimum = passedCut( "PassJSON" ) && passedCut ( "PassTrigger" );

    //-----------------------------------------------------------------    
    // Declare purity requirements
    //-----------------------------------------------------------------    

    bool jet1_pass              = false;
    bool jet2_pass              = false;
    bool ele1_pass              = false;
    bool ele2_pass_eejj         = false;
    bool met_pass_enujj         = false;
    bool ele_pass_veto_enujj    = false;
    bool ele_pass_veto_eejj     = false;

    bool passAll_ejj            = false;
    bool passAll_eejj_veto      = false;
    bool passAll_eejj_noveto    = false;
    bool passAll_enujj_veto     = false;
    bool passAll_enujj_noveto   = false;
    
    //-----------------------------------------------------------------    
    // Evaluate purity requirements
    //-----------------------------------------------------------------    
    
    // electrons

    if ( n_ele >= 1 ) { 
      Electron ele1 = c_ele_final -> GetConstituent<Electron>(0);
      if ( ele1.Pt() > 45. ) { 
	ele1_pass  = true;
      }
    }

    if ( n_ele >= 2 ) { 
      Electron ele2 = c_ele_final -> GetConstituent<Electron>(1);
      if ( ele2.Pt() > 45. ) ele2_pass_eejj = true;
    }

    if ( n_ele == 1 ) ele_pass_veto_enujj = true;
    if ( n_ele == 2 ) ele_pass_veto_eejj  = true;
    
    // jets ( same for both analyses ) 
    

    if ( n_jet >= 1 ){
      PFJet jet1 = c_pfjet_final -> GetConstituent<PFJet>(0);
      if ( jet1.Pt () >= 125. ) jet1_pass = true;
    }
    
    if ( n_jet >= 2 ){
      PFJet jet2 = c_pfjet_final -> GetConstituent<PFJet>(1);
      if ( jet2.Pt () >= 45. ) jet2_pass = true;
    }

    if ( PFMETType01XYCor -> at (0) > 45. ) met_pass_enujj = true;

    //-----------------------------------------------------------------    
    // Combine eejj purity requirements
    //-----------------------------------------------------------------    
    
    if ( ele1_pass           &&
	 jet1_pass           && 
	 jet2_pass           ){
      passAll_ejj          = true;
    }

    if ( ele1_pass           &&
	 ele2_pass_eejj      && 
	 jet1_pass           && 
	 jet2_pass           ){
      passAll_eejj_noveto  = true;
    }

    if ( ele1_pass           &&
	 ele2_pass_eejj      && 
	 ele_pass_veto_eejj  && 
	 jet1_pass           && 
	 jet2_pass           ){
      passAll_eejj_veto    = true;
    }

    if ( ele1_pass           &&
	 met_pass_enujj      && 
	 jet1_pass           && 
	 jet2_pass           ){
      passAll_enujj_noveto = true;
    }

    if ( ele1_pass           &&
	 ele_pass_veto_enujj && 
	 met_pass_enujj      && 
	 jet1_pass           && 
	 jet2_pass           ){
      passAll_enujj_veto   = true;
    }
    
    //-----------------------------------------------------------------    
    // Act on outcome of purity requirements
    //-----------------------------------------------------------------    
    
    FillUserTH1D ( "nEle_all" , n_ele );
    
    if ( passAll_ejj          ) FillUserTH1D ( "nEle_pass_ejj"          , n_ele );
    if ( passAll_eejj_veto    ) FillUserTH1D ( "nEle_pass_eejj_veto"    , n_ele );
    if ( passAll_eejj_noveto  ) FillUserTH1D ( "nEle_pass_eejj_noveto"  , n_ele );
    if ( passAll_enujj_veto   ) FillUserTH1D ( "nEle_pass_enujj_veto"   , n_ele );
    if ( passAll_enujj_noveto ) FillUserTH1D ( "nEle_pass_enujj_noveto" , n_ele );
    					                               
    if ( met_pass_enujj       ) FillUserTH1D ( "MET_metPass_enujj"      , PFMETType01XYCor -> at (0) );
    if ( jet1_pass            ) FillUserTH1D ( "nJet_jet1Pass"          , n_jet );
    if ( jet2_pass            ) FillUserTH1D ( "nJet_jet2Pass"          , n_jet );
    if ( ele1_pass            ) FillUserTH1D ( "nEle_ele1Pass"          , n_ele );
    if ( ele2_pass_eejj       ) FillUserTH1D ( "nEle_ele2Pass_eejj"     , n_ele );
    if ( ele_pass_veto_enujj  ) FillUserTH1D ( "nEle_passVeto_enujj"    , n_ele ); 
    if ( ele_pass_veto_eejj   ) FillUserTH1D ( "nEle_passVeto_eejj"     , n_ele ); 
    
  }

  std::cout << "analysisClass::Loop() ends" <<std::endl;   
  
}
