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
#include "HLTriggerObject.h"


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

  //--------------------------------------------------------------------------
  // Declare plots
  //--------------------------------------------------------------------------


  CreateUserTH1D ( "HLTPFJet1_Pt_pass_v1"   , 200, 0, 200.0);
  CreateUserTH1D ( "HLTPFJet1_Pt_pass_v2"   , 200, 0, 200.0);
  CreateUserTH1D ( "HLTPFJet1_Pt_pass_v3"   , 200, 0, 200.0);
  CreateUserTH1D ( "HLTPFJet1_Pt_pass_v4"   , 200, 0, 200.0);
  CreateUserTH1D ( "HLTPFJet1_Pt_pass_v5"   , 200, 0, 200.0);
  CreateUserTH1D ( "HLTPFJet1_Pt_pass_v6"   , 200, 0, 200.0);
  CreateUserTH1D ( "HLTPFJet1_Pt_pass_v7"   , 200, 0, 200.0);
  CreateUserTH1D ( "HLTPFJet1_Pt_pass_v8"   , 200, 0, 200.0);

  CreateUserTH1D( "PFJetMatched_Pt_pass_v1" , 200, 0, 200.0);
  CreateUserTH1D( "PFJetMatched_Pt_pass_v2" , 200, 0, 200.0);
  CreateUserTH1D( "PFJetMatched_Pt_pass_v3" , 200, 0, 200.0);
  CreateUserTH1D( "PFJetMatched_Pt_pass_v4" , 200, 0, 200.0);
  CreateUserTH1D( "PFJetMatched_Pt_pass_v5" , 200, 0, 200.0);
  CreateUserTH1D( "PFJetMatched_Pt_pass_v6" , 200, 0, 200.0);
  CreateUserTH1D( "PFJetMatched_Pt_pass_v7" , 200, 0, 200.0);
  CreateUserTH1D( "PFJetMatched_Pt_pass_v8" , 200, 0, 200.0);

  CreateUserTH1D( "PFJet1_Pt_pass_v1" , 200, 0, 200.0);
  CreateUserTH1D( "PFJet1_Pt_pass_v2" , 200, 0, 200.0);
  CreateUserTH1D( "PFJet1_Pt_pass_v3" , 200, 0, 200.0);
  CreateUserTH1D( "PFJet1_Pt_pass_v4" , 200, 0, 200.0);
  CreateUserTH1D( "PFJet1_Pt_pass_v5" , 200, 0, 200.0);
  CreateUserTH1D( "PFJet1_Pt_pass_v6" , 200, 0, 200.0);
  CreateUserTH1D( "PFJet1_Pt_pass_v7" , 200, 0, 200.0);
  CreateUserTH1D( "PFJet1_Pt_pass_v8" , 200, 0, 200.0);


  CreateUserTH1D ( "HLTPFJet1_Pt_total_v1"   , 200, 0, 200.0);
  CreateUserTH1D ( "HLTPFJet1_Pt_total_v2"   , 200, 0, 200.0);
  CreateUserTH1D ( "HLTPFJet1_Pt_total_v3"   , 200, 0, 200.0);
  CreateUserTH1D ( "HLTPFJet1_Pt_total_v4"   , 200, 0, 200.0);
  CreateUserTH1D ( "HLTPFJet1_Pt_total_v5"   , 200, 0, 200.0);
  CreateUserTH1D ( "HLTPFJet1_Pt_total_v6"   , 200, 0, 200.0);
  CreateUserTH1D ( "HLTPFJet1_Pt_total_v7"   , 200, 0, 200.0);
  CreateUserTH1D ( "HLTPFJet1_Pt_total_v8"   , 200, 0, 200.0);

  CreateUserTH1D( "PFJetMatched_Pt_total_v1" , 200, 0, 200.0);
  CreateUserTH1D( "PFJetMatched_Pt_total_v2" , 200, 0, 200.0);
  CreateUserTH1D( "PFJetMatched_Pt_total_v3" , 200, 0, 200.0);
  CreateUserTH1D( "PFJetMatched_Pt_total_v4" , 200, 0, 200.0);
  CreateUserTH1D( "PFJetMatched_Pt_total_v5" , 200, 0, 200.0);
  CreateUserTH1D( "PFJetMatched_Pt_total_v6" , 200, 0, 200.0);
  CreateUserTH1D( "PFJetMatched_Pt_total_v7" , 200, 0, 200.0);
  CreateUserTH1D( "PFJetMatched_Pt_total_v8" , 200, 0, 200.0);

  CreateUserTH1D( "PFJet1_Pt_total_v1" , 200, 0, 200.0);
  CreateUserTH1D( "PFJet1_Pt_total_v2" , 200, 0, 200.0);
  CreateUserTH1D( "PFJet1_Pt_total_v3" , 200, 0, 200.0);
  CreateUserTH1D( "PFJet1_Pt_total_v4" , 200, 0, 200.0);
  CreateUserTH1D( "PFJet1_Pt_total_v5" , 200, 0, 200.0);
  CreateUserTH1D( "PFJet1_Pt_total_v6" , 200, 0, 200.0);
  CreateUserTH1D( "PFJet1_Pt_total_v7" , 200, 0, 200.0);
  CreateUserTH1D( "PFJet1_Pt_total_v8" , 200, 0, 200.0);

  //--------------------------------------------------------------------------
  // Tell the user how many entries we'll look at
  //--------------------------------------------------------------------------

  // Long64_t nentries = fChain->GetEntries();
  Long64_t nentries = 10000;
  std::cout << "analysisClass::Loop(): nentries = " << fChain -> GetEntries() << std::endl;   
  
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
    //for (Long64_t jentry = 0; jentry<5;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //-----------------------------------------------------------------
    // Print progress
    //-----------------------------------------------------------------
    
    if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
    
    //-----------------------------------------------------------------
    // Get JSON information
    //-----------------------------------------------------------------

    int passedJSON = passJSON ( run, ls , isData ) ;
    if ( passedJSON == 0 ) continue;
    
    //-----------------------------------------------------------------
    // Get trigger information
    //-----------------------------------------------------------------

    getTriggers ( HLTKey, HLTInsideDatasetTriggerNames, HLTInsideDatasetTriggerDecisions,  HLTInsideDatasetTriggerPrescales ) ; 

    char ele_jet_trigger_name[200];
    char ele_jet_trigger_nopu_name[200];
    char ele_trigger_name    [200];
    
    int ele_jet_trigger_version = 1;
    int ele_jet_trigger_nopu_version = 1;
    int ele_trigger_version     = 1;

    //if      ( run >= 190456 && run <= 190738) ele_jet_trigger_version = 3;
    //else if ( run >= 190782 && run <= 191419) ele_jet_trigger_version = 4;
    //else if ( run >= 191691 && run <= 196027) ele_jet_trigger_version = 5;

    //if      ( run >= 191691 && run <= 194225 ) ele_jet_trigger_nopu_version = 4;
    //else if ( run >= 194270 && run <= 196531 ) ele_jet_trigger_nopu_version = 5;	
    //else if ( run >= 198022 && run <= 199608 ) ele_jet_trigger_nopu_version = 6;	
    //else if ( run >= 199698 && run <= 202504 ) ele_jet_trigger_nopu_version = 7;	
    //else if ( run >= 202970                  ) ele_jet_trigger_nopu_version = 8;

    //if      ( run >= 190456 && run <= 190738 ) ele_trigger_version = 3;
    //else if ( run >= 190782 && run <= 191419 ) ele_trigger_version = 4;
    //else if ( run >= 191691 && run <= 196531 ) ele_trigger_version = 5;
    //else if ( run >= 198022                  ) ele_trigger_version = 6;

    //sprintf(ele_jet_trigger_nopu_name, "HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v%d", ele_jet_trigger_nopu_version);	
    sprintf(ele_jet_trigger_name     , "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v%d"        , ele_jet_trigger_version     );	
    sprintf(ele_trigger_name         , "HLT_Ele27_eta2p1_WP85_Gsf_v%d"                         , ele_trigger_version         );
    // was HLT_Ele30_CaloIdVT_TrkIdT_v*
    
    bool ele_jet_trigger_fired    = triggerFired    ( ele_jet_trigger_name );
    bool ele_jet_trigger_enabled  = ( triggerPrescale ( ele_jet_trigger_name ) != -999 );
    //bool ele_jet_trigger_nopu_fired    = triggerFired    ( ele_jet_trigger_nopu_name );
    //bool ele_jet_trigger_nopu_enabled  = ( triggerPrescale ( ele_jet_trigger_nopu_name ) != -999 );
    bool ele_trigger_fired        = triggerFired    ( ele_trigger_name );
    bool ele_trigger_enabled      = ( triggerPrescale ( ele_trigger_name ) != -999 );
    
    assert ( ele_jet_trigger_version > 0 );
    assert ( ele_jet_trigger_enabled );
    assert ( ele_trigger_version > 0 );
    assert ( ele_trigger_enabled );

    //-----------------------------------------------------------------
    // Get HLT Filter objects
    //-----------------------------------------------------------------
    
    CollectionPtr triggerObjs_all (new Collection(*this, HLTriggerObjPt->size() ));
    CollectionPtr triggerObjs_eles_all  = triggerObjs_all->SkimByID<HLTriggerObject>(TRIGGER_PHOTON);
    CollectionPtr triggerObjs_jets_all  = triggerObjs_all->SkimByID<HLTriggerObject>(TRIGGER_JET);
    CollectionPtr triggerObjs_200jet_all (new Collection(*this, 0 ));
    triggerObjs_200jet_all->Clear();
    for(int i=0; i < triggerObjs_jets_all->GetSize(); ++i)
    {
      HLTriggerObject jet = triggerObjs_jets_all->GetConstituent<HLTriggerObject>(i);
      if(jet.PassedPathLastFilter(std::string(ele_jet_trigger_name)))
        triggerObjs_200jet_all->Append(jet.GetRawIndex());
    }
    //std::cout << "------------------------------------------------------------------------" << std::endl;
    //std::cout << "TriggerObjects size: " << triggerObjs_all->GetSize() << std::endl;
    //std::cout << "TriggerObjectsJet size: " << triggerObjs_jets_all->GetSize() << std::endl;
    //std::cout << "TriggerObjectsJet200 size: " << triggerObjs_200jet_all->GetSize() << std::endl;

    HLTriggerObject leadHLTJet = triggerObjs_jets_all -> GetConstituent<HLTriggerObject>(0);
    CollectionPtr triggerObjs_jets_no200Jet;
    if(triggerObjs_200jet_all->GetSize() > 0)
      triggerObjs_jets_no200Jet = triggerObjs_jets_all->SkimByVetoDRMatch<HLTriggerObject,HLTriggerObject>(triggerObjs_200jet_all,0.3);
    else
      triggerObjs_jets_no200Jet = CollectionPtr(new Collection(*this, 0));


    int n_trigger_no200jet_all  = triggerObjs_jets_no200Jet  -> GetSize();
    int n_trigger_200jet_all = triggerObjs_200jet_all -> GetSize();
    int n_trigger_jets_all = triggerObjs_jets_all -> GetSize();

    //-----------------------------------------------------------------
    // Matching to RECO objects
    //-----------------------------------------------------------------
    
    CollectionPtr ele_all   ( new Collection(*this, ElectronPt -> size() ));
    CollectionPtr ele_HEEP = ele_all -> SkimByID<Electron>( HEEP51 );


    CollectionPtr pfjet_all ( new Collection(*this, PFJetPt -> size() ));
    CollectionPtr pfjet_central = pfjet_all -> SkimByEtaRange<PFJet> ( -2.6,2.6);
    CollectionPtr pfjet_hltJetsMatched = pfjet_all -> SkimByRequireDRMatch<PFJet,HLTriggerObject>(triggerObjs_jets_all, 0.5);
    PFJet pfjet_matchedToLeadHLTJet = pfjet_hltJetsMatched -> GetClosestInDR<PFJet> ( leadHLTJet );

    //-----------------------------------------------------------------
    // Make turn-on curve
    //-----------------------------------------------------------------
    
    if ( ele_trigger_fired && n_trigger_jets_all >= 2 && ele_jet_trigger_enabled)
    {
      char total_name[200];

      sprintf(total_name, "PFJet1_Pt_total_v%d"      , ele_jet_trigger_version ); FillUserTH1D(total_name, pfjet_central -> GetConstituent<PFJet>(0).Pt());
      sprintf(total_name, "PFJetMatched_Pt_total_v%d", ele_jet_trigger_version ); FillUserTH1D(total_name, pfjet_matchedToLeadHLTJet.Pt());
      sprintf(total_name, "HLTPFJet1_Pt_total_v%d"   , ele_jet_trigger_version ); FillUserTH1D(total_name, leadHLTJet.Pt() );

      if ( ele_jet_trigger_fired ) { 
        char PFJet1_Pt_pass_name        [100]; sprintf(PFJet1_Pt_pass_name      , "PFJet1_Pt_pass_v%d"       , ele_jet_trigger_version);
        char PFJetMatched_Pt_pass_name  [100]; sprintf(PFJetMatched_Pt_pass_name, "PFJetMatched_Pt_pass_v%d" , ele_jet_trigger_version);
        char HLTPFJet1_Pt_pass_name     [100]; sprintf(HLTPFJet1_Pt_pass_name   , "HLTPFJet1_Pt_pass_v%d"    , ele_jet_trigger_version);

        FillUserTH1D(PFJet1_Pt_pass_name       , pfjet_central -> GetConstituent<PFJet>(0).Pt());
        FillUserTH1D(PFJetMatched_Pt_pass_name , pfjet_matchedToLeadHLTJet.Pt());
        FillUserTH1D(HLTPFJet1_Pt_pass_name    , leadHLTJet.Pt());
      }
    }

    //-----------------------------------------------------------------
    // Examine
    //-----------------------------------------------------------------

    if ( ele_jet_trigger_fired )
    {
      std::cout << "------------------------------------------------------------------------" << std::endl;
      std::cout << "Run = " << run << ", event = " << event << ", LS = " << ls << std::endl;
      triggerObjs_jets_all      -> examine<HLTriggerObject>("HLT Jet 50");
      pfjet_hltJetsMatched  -> examine<PFJet>          ("PFJets, matched to 50 GeV HLT jets");
      std::cout << "PFJet matched to lead HLT jet = " << pfjet_matchedToLeadHLTJet << std::endl;
      ele_all                -> examine<Electron>       ("Electrons");
    }

  } // End loop over events
  
  std::cout << "analysisClass::Loop() ends" <<std::endl;   
  
}
