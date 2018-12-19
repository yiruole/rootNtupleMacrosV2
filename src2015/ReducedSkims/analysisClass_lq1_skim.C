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
#include "LooseElectron.h"
#include "Muon.h"
#include "PFJet.h"
#include "GenJet.h"
#include "HLTriggerObject.h"
#include "HLTriggerObjectCollectionHelper.h"
// for 2016 L1 prescales
#include "src/EGPrescales2016.C"

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
  //setupTriggerReaders();
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

  //--------------------------------------------------------------------------
  // Plots
  //--------------------------------------------------------------------------
  CreateUserTH1D( "nEleNTuple"                       ,    10    , 0   , 10     );
  CreateUserTH1D( "nEleNrsk"                       ,    10    , 0   , 10     );
  CreateUserTH2D( "nEleNTupleVsNeleRsk"                 ,     10, 0, 10, 10, 0, 10) ;

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
  // - 3: single electron (HEEP)
  // - 4: single muon (tight)
  //--------------------------------------------------------------------------

  int reducedSkimType = getPreCutValue1("reducedSkimType");

  //--------------------------------------------------------------------------
  // Should we do any scaling / smearing for systematics?
  //--------------------------------------------------------------------------

  int electron_energy_scale_sign = int(getPreCutValue1("electron_energy_scale_sign" ));
  int pfjet_energy_scale_sign    = int(getPreCutValue1("pfjet_energy_scale_sign"    ));

  bool do_eer = bool ( int(getPreCutValue1("do_electron_energy_smear"   )) == 1 );
  bool do_jer = bool ( int(getPreCutValue1("do_pfjet_energy_smear"      )) == 1 );
  bool do_ees = bool ( electron_energy_scale_sign != 0 );
  bool do_jes = bool ( pfjet_energy_scale_sign    != 0 );
  
  //--------------------------------------------------------------------------
  // Cuts for physics objects selection
  //--------------------------------------------------------------------------

  // jet cuts
  double jet_PtCut               = getPreCutValue1("jet_PtCut");
  double jet_EtaCut              = getPreCutValue1("jet_EtaCut");
  double jet_HighEtaCut          = getPreCutValue1("jet_HighEtaCut");
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
  std::cout << "analysisClass::Loop(): nentries = " << fChain -> GetEntries() << std::endl;   
  
  //--------------------------------------------------------------------------
  // Create HLT collections in advance (won't need all of them)
  //--------------------------------------------------------------------------
  
  // QCD photon filters
  CollectionPtr c_hltPhoton22_QCD_all;
  CollectionPtr c_hltPhoton30_QCD_all;
  CollectionPtr c_hltPhoton36_QCD_all;
  CollectionPtr c_hltPhoton50_QCD_all;
  CollectionPtr c_hltPhoton75_QCD_all;
  CollectionPtr c_hltPhoton90_QCD_all;
  CollectionPtr c_hltPhoton120_QCD_all;	      
  CollectionPtr c_hltPhoton175_QCD_all;	      

  // TTbar muon filter
  //XXX SIC remove this for now
  //CollectionPtr c_hltMuon22_TTbar_all;
  //CollectionPtr c_hltPhoton22_TTbar_all;

  // muon filter
  CollectionPtr c_hltMuon_SingleMu_all;

  // Signal
  CollectionPtr c_hltEle45_Signal_all;
  CollectionPtr c_hltPFJet50_Signal_all;
  CollectionPtr c_hltPFJet200_Signal_all;
  CollectionPtr c_hltDoubleEle_Signal_all;

  // trigger
  CollectionPtr c_trigger_l3jets_all;

  // Tag and probe
  CollectionPtr c_hltEle27WP85Gsf_all;

  //--------------------------------------------------------------------------
  // For smearing systematics samples, you'll need a random number generator
  //--------------------------------------------------------------------------

  unsigned int seed = 987654321;
  TRandom3 * rootEngine = new TRandom3 ( seed ) ;

  //--------------------------------------------------------------------------
  // For smearing/scaling systematics samples, you'll need to calculate the effect on MET
  //--------------------------------------------------------------------------
  
  TLorentzVector v_delta_met;

  TLorentzVector v_PFMETRaw;
  TLorentzVector v_PFMETType1Cor;
  TLorentzVector v_PFMETType1XYCor;
  

  // prepare 2016 EG prescales
  EGPrescales2016 egPrescales2016;

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
    //if (ientry < 0) break;
    if (ientry < 0)
    {
      std::cout << "ERROR: Could not read from TTree; exiting." << std::endl;
      exit(-1);
    }
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (nb < 0)
    {
      std::cout << "ERROR: Could not read entry from TTree: read " << nb << "bytes; exiting." << std::endl;
      exit(-2);
    }

    //-----------------------------------------------------------------
    // Print progress
    //-----------------------------------------------------------------
    if(jentry < 10 || jentry%5000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << "/" << nentries << std::endl;
    //// run ls event
    //std::cout << static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(ls) << " " << static_cast<unsigned int>(event) << std::endl;

    //-----------------------------------------------------------------
    // Get access to HLT decisions
    //-----------------------------------------------------------------

    getTriggers(jentry); 
    //printTriggers();
    //printFiredTriggers();
        
    //-----------------------------------------------------------------
    // Get access to HLT filter objects
    //-----------------------------------------------------------------
    
    HLTriggerObjectCollectionHelper helper (*this,"");

    if ( reducedSkimType == 0 ){ 
      
      // QCD photon triggers
      // taken from singlePhoton triggers from menu here:
      //    https://cmsweb-testbed.cern.ch/confdb/#config=/cdaq/physics/Run2015/25ns14e33/v4.4.5/HLT/V1
      // was menu used in run 260627 (last of 2015 pp run, taken Nov. 3)

      c_hltPhoton22_QCD_all = helper.GetLastFilterObjectsByPath("HLT_Photon22_v");
      c_hltPhoton30_QCD_all = helper.GetLastFilterObjectsByPath("HLT_Photon30_v");
      c_hltPhoton36_QCD_all = helper.GetLastFilterObjectsByPath("HLT_Photon36_v");
      c_hltPhoton50_QCD_all = helper.GetLastFilterObjectsByPath("HLT_Photon50_v");
      c_hltPhoton75_QCD_all = helper.GetLastFilterObjectsByPath("HLT_Photon75_v");
      c_hltPhoton90_QCD_all = helper.GetLastFilterObjectsByPath("HLT_Photon90_v");
      c_hltPhoton120_QCD_all = helper.GetLastFilterObjectsByPath("HLT_Photon120_v");
      c_hltPhoton175_QCD_all = helper.GetLastFilterObjectsByPath("HLT_Photon175_v");
      
    }

    else if ( reducedSkimType == 1 || reducedSkimType == 2 || reducedSkimType == 3 || reducedSkimType == 4){
      
      // TTbar trigger: path was HLT_Mu22_Photon22_CaloIdL_v3+
      // saveTags is false for this path? how did that work?
      // also this is not selecting muons either...
      
      // SingleMu trigger: path was HLT_Mu40_eta2p1_v9+
      c_hltMuon_SingleMu_all       = helper.GetL3FilterObjectsByPath("HLT_Mu45_eta2p1_v"); // will do prefix matching

      // Ele+jets signal triggers
      CollectionPtr trigger_l3objects_all = helper.GetL3FilterObjectsByPath("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v");
      c_trigger_l3jets_all = trigger_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_JET);
      //c_trigger_l3jets_all->examine<HLTriggerObject>("c_trigger_l3jets_all");
      // which one passed the last filter?
      CollectionPtr trigger_lastObjects_all = helper.GetLastFilterObjectsByPath("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v");

      // electrons seem to come as TRIGGER_PHOTON most of the time
      c_hltEle45_Signal_all = trigger_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_PHOTON);
      // if not, try TRIGGER_ELECTRON
      if(c_hltEle45_Signal_all->GetSize() == 0)
       c_hltEle45_Signal_all = trigger_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_ELECTRON);
      // Note: could also be TRIGGER_CLUSTER?

      // jets
      c_hltPFJet200_Signal_all     =  trigger_lastObjects_all->SkimByID<HLTriggerObject>(TRIGGER_JET);
      // get rid of overlaps
      // XXX no, keep all L3 jets
      //c_hltPFJet50_Signal_all      =  c_trigger_l3jets_all -> SkimByVetoDRMatch<HLTriggerObject,HLTriggerObject>( c_hltPFJet200_Signal_all, 0.3 );

      // DoubleEle signal trigger
      CollectionPtr double_ele_l3objects_all    = helper.GetL3FilterObjectsByPath("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v");
      c_hltDoubleEle_Signal_all    = double_ele_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_PHOTON);
      // if not, try TRIGGER_ELECTRON
      if(c_hltDoubleEle_Signal_all->GetSize() == 0)
       c_hltDoubleEle_Signal_all = double_ele_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_ELECTRON);
      
      // Tag and probe trigger
      CollectionPtr tagProbe_l3objects_all;
      if(!isData())
        tagProbe_l3objects_all  = helper.GetL3FilterObjectsByPath("HLT_Ele27_WPLoose_Gsf_v");
      else
        tagProbe_l3objects_all  = helper.GetL3FilterObjectsByPath("HLT_Ele27_WPLoose_Gsf_v");
      c_hltEle27WP85Gsf_all = tagProbe_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_PHOTON);
      // if not, try TRIGGER_ELECTRON
      if(c_hltEle27WP85Gsf_all->GetSize() == 0)
       c_hltEle27WP85Gsf_all = tagProbe_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_ELECTRON);
      
    }


    //-----------------------------------------------------------------
    // Define initial, inclusive collections for physics objects
    //-----------------------------------------------------------------

    CollectionPtr c_gen_all   ( new Collection(*this, nGenPart));
    CollectionPtr c_ele_all   ( new Collection(*this, nElectron));
    //c_ele_all->examine<Electron>("c_ele_all = All reco electrons");
    //Electron ele1 = c_ele_all -> GetConstituent<Electron>(0);
    //for(unsigned int i=0; i<10; ++i) { 
    //  std::cout << "cut = " << i << " idLevel = " << ele1.GetNbitFromBitMap(i, 3) << std::endl;
    //}

    CollectionPtr c_muon_all  ( new Collection(*this, nMuon));
    CollectionPtr c_genJet_all( new Collection(*this, nGenJet));
    CollectionPtr c_pfjet_all ( new Collection(*this, nJet));

    //-----------------------------------------------------------------
    // All skims need GEN particles/jets
    //-----------------------------------------------------------------

    //c_gen_all->examine<GenParticle>("c_gen_all");
    CollectionPtr c_genEle_final = c_gen_all    -> SkimByID<GenParticle>(GEN_ELE_HARD_SCATTER);
    //c_genEle_final->examine<GenParticle>("c_genEle_final = c_gen_all after SkimByID<GenParticle>(GEN_ELE_HARD_SCATTER)");

    CollectionPtr c_genMu_final = c_gen_all    -> SkimByID<GenParticle>(GEN_MU_HARD_SCATTER);

    CollectionPtr c_genNu_final = c_gen_all    -> SkimByID<GenParticle>(GEN_NU_HARD_SCATTER);
    //c_genNu_final->examine<GenParticle>("c_genNu_final = c_gen_all after SkimByID<GenParticle>(GEN_NU_HARD_SCATTER)");
    CollectionPtr c_genJet_final = c_genJet_all;

    CollectionPtr c_genZgamma_final    = c_gen_all -> SkimByID<GenParticle>(GEN_ZGAMMA_HARD_SCATTER);
    //c_genZgamma_final->examine<GenParticle>("c_genZgamma_final = c_gen_all after SkimByID<GenParticle>(GEN_ZGAMMA_HARD_SCATTER)");
    CollectionPtr c_genW_final         = c_gen_all -> SkimByID<GenParticle>(GEN_W_HARD_SCATTER);
    //c_genW_final->examine<GenParticle>("c_genW_final = c_gen_all after SkimByID<GenParticle>(GEN_W_HARD_SCATTER)");
    //if(c_genW_final->GetSize()<1) {
    //  cout << "DID NOT FIND A W THIS EVENT! HERE ARE ALL THE GENPARTICLES" << endl;
    //  c_gen_all->examine<GenParticle>("c_gen_all");
    //}
    CollectionPtr c_genNuFromW_final   = c_gen_all -> SkimByID<GenParticle>(GEN_NU_FROM_W);
    CollectionPtr c_genEleFromW_final  = c_gen_all -> SkimByID<GenParticle>(GEN_ELE_FROM_W);
    CollectionPtr c_genEleFromDY_final = c_gen_all -> SkimByID<GenParticle>(GEN_ELE_FROM_DY);
    //c_genEleFromDY_final->examine<GenParticle>("c_genEleFromDY_final = c_gen_all after SkimByID<GenParticle>(GEN_ELE_FROM_DY)");
    //c_genEleFromW_final->examine<GenParticle>("c_genEleFromW_final = c_gen_all after SkimByID<GenParticle>(GEN_ELE_FROM_W)");
    //c_genNuFromW_final->examine<GenParticle>("c_genNuFromW_final = c_gen_all after SkimByID<GenParticle>(GEN_NU_FROM_W)");
    CollectionPtr c_genTop             = c_gen_all -> SkimByID<GenParticle>(GEN_TOP);
    CollectionPtr c_genTop_final       = c_genTop  -> SkimByID<GenParticle>(GEN_STATUS62);

    //-----------------------------------------------------------------
    // If this is MC, we're always going to smear the jets (do_jer = true)
    // Don't do it for data
    //-----------------------------------------------------------------

    if ( !isData() ) do_jer = true;
    else do_jer = false; // not for data
      
    //-----------------------------------------------------------------
    // Energy scaling and resolution smearing here
    //-----------------------------------------------------------------
    
    //SIC FIXME TODO: update this
    if ( do_eer || do_jer || do_ees || do_jes ) { 
      
      // If  you're scaling/smearing PFJets, recall that only jets with pt > 10 GeV affect the PFMET
      //  New: This is now 15 GeV
      // Also, only scale/smear the jets in our eta range (jets in the calorimeter crack are suspect)
      c_pfjet_all = c_pfjet_all -> SkimByMinPt   <PFJet>( 15.0 );
      c_pfjet_all = c_pfjet_all -> SkimByEtaRange<PFJet>( -jet_EtaCut, jet_EtaCut );
      
      // Set the PFMET difference to zero
      
      v_delta_met.SetPtEtaPhiM(0.,0.,0.,0.);

      // Do the energy scale / energy resolution operations
      // dR for matching = Rcone/2
      //   see: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
      
      if ( do_eer ) c_ele_all      -> MatchAndSmearEnergy <Electron, GenParticle> ( c_genEle_final, 0.4/2.0, rootEngine, v_delta_met );
      if ( do_jer ) c_pfjet_all    -> MatchAndSmearEnergy <PFJet   , GenJet     > ( c_genJet_final, 0.4/2.0, rootEngine, v_delta_met );
      if ( do_ees ) c_ele_all      -> ScaleEnergy <Electron> ( electron_energy_scale_sign, v_delta_met );
      if ( do_jes ) c_pfjet_all    -> ScaleEnergy <PFJet   > ( pfjet_energy_scale_sign   , v_delta_met );
      
      // Propagate the results to the PFMET

      //v_PFMETRaw        .SetPtEtaPhiM( (*PFMETRaw        )[0] , 0., (*PFMETPhiRaw        )[0] , 0. );
      v_PFMETType1Cor   .SetPtEtaPhiM( MET_pt, 0., MET_phi , 0. );
      //v_PFMETType1XYCor.SetPtEtaPhiM( (*PFMETType1XYCor)[0] , 0., (*PFMETPhiType1XYCor)[0] , 0. );
      
      //v_PFMETRaw         = v_PFMETRaw         + v_delta_met;
      v_PFMETType1Cor    = v_PFMETType1Cor    + v_delta_met;
      //v_PFMETType1XYCor = v_PFMETType1XYCor + v_delta_met;
      
      //(*PFMETRaw           )[0] = v_PFMETRaw        .Pt();
      //FIXME
      //(*PFMETType1Cor      )[0] = v_PFMETType1Cor   .Pt();
      //(*PFMETType1XYCor   )[0] = v_PFMETType1XYCor.Pt();
      //
      //(*PFMETPhiRaw        )[0] = v_PFMETRaw        .Phi();
      //FIXME
      //(*PFMETPhiType1Cor   )[0] = v_PFMETType1Cor   .Phi();
      //(*PFMETPhiType1XYCor)[0] = v_PFMETType1XYCor.Phi();
    }

    //-----------------------------------------------------------------
    // QCD skims    (reducedSkimType = 0     ) have loose electrons
    // Signal skims (reducedSkimType = 1 - 4 ) have HEEP  electrons
    //-----------------------------------------------------------------

    CollectionPtr c_ele_final;
    CollectionPtr c_ele_final_ptCut;
    
    if ( reducedSkimType == 0 ){ 
      CollectionPtr c_ele_loose = c_ele_all   -> SkimByID  <LooseElectron> ( FAKE_RATE_HEEP_LOOSE);
      c_ele_final               = c_ele_loose;
      c_ele_final_ptCut         = c_ele_final -> SkimByMinPt<LooseElectron>( ele_PtCut  );
    }

    else if ( reducedSkimType == 1 || reducedSkimType == 2 || reducedSkimType == 3 || reducedSkimType == 4 ){
      CollectionPtr c_ele_HEEP  = c_ele_all -> SkimByID <Electron> ( HEEP70 );
      //CollectionPtr c_ele_HEEP  = c_ele_all -> SkimByID <Electron> ( HEEP70_MANUAL , true );
      c_ele_final               = c_ele_HEEP;
      c_ele_final_ptCut         = c_ele_final -> SkimByMinPt<Electron>( ele_PtCut  );
    }
    // look at final electrons
    //c_ele_final->examine<Electron>("c_ele_final");
    //c_ele_final_ptCut->examine<Electron>("c_ele_final_ptCut");
    FillUserTH1D("nEleNTuple",c_ele_all->GetSize());
    FillUserTH1D("nEleNrsk",c_ele_final->GetSize());
    FillUserTH2D("nEleNTupleVsNeleRsk",c_ele_final->GetSize(),c_ele_all->GetSize());

    //-----------------------------------------------------------------
    // All skims need muons
    //-----------------------------------------------------------------

    CollectionPtr c_muon_eta               = c_muon_all       -> SkimByEtaRange<Muon> ( -muon_EtaCut, muon_EtaCut );
    //CollectionPtr c_muon_eta_IDTight       = c_muon_eta       -> SkimByID      <Muon> ( MUON_TIGHT_PFISO04TIGHT );
    CollectionPtr c_muon_eta_IDHighPt      = c_muon_eta       -> SkimByID      <Muon> ( MUON_HIGH_PT_TRKRELISO03 );
    //CollectionPtr c_muon_eta_IDLoose       = c_muon_eta       -> SkimByID      <Muon> ( MUON_LOOSE );
    CollectionPtr c_muon_final             = c_muon_eta_IDHighPt;
    //CollectionPtr c_muon_final             = c_muon_eta_IDLoose;
    CollectionPtr c_muon_final_ptCut       = c_muon_final     -> SkimByMinPt   <Muon> ( muon_PtCut );
    //c_muon_final->examine<Muon>("c_muon_final");
    
    //-----------------------------------------------------------------
    // All skims need PFJets
    //-----------------------------------------------------------------

    CollectionPtr c_pfjet_central                     = c_pfjet_all                          -> SkimByEtaRange   <PFJet>          ( -jet_EtaCut, jet_EtaCut   );
    CollectionPtr c_pfjet_central_ID                  = c_pfjet_central                      -> SkimByID         <PFJet>          ( PFJET_LOOSE );    
    CollectionPtr c_pfjet_central_ID_noMuonOverlap    = c_pfjet_central_ID                   -> SkimByVetoDRMatch<PFJet, Muon>    ( c_muon_final_ptCut   , jet_muon_DeltaRCut  );
    CollectionPtr c_pfjet_central_ID_noLeptonOverlap  = c_pfjet_central_ID_noMuonOverlap     -> SkimByVetoDRMatch<PFJet, Electron>( c_ele_final_ptCut    , jet_ele_DeltaRCut );
    CollectionPtr c_pfjet_final                       = c_pfjet_central_ID_noLeptonOverlap;
    CollectionPtr c_pfjet_final_ptCut                 = c_pfjet_final                        -> SkimByMinPt      <PFJet>          ( jet_PtCut );
    //c_pfjet_final->examine<PFJet>("c_pfjet_final");
    
    //-----------------------------------------------------------------
    // We need high-eta jets in order to look at boson recoil
    //-----------------------------------------------------------------
    
    CollectionPtr c_pfjet_highEta                    = c_pfjet_all                          -> SkimByEtaRange   <PFJet>          ( -jet_HighEtaCut, jet_HighEtaCut   );
    CollectionPtr c_pfjet_highEta_ID                 = c_pfjet_highEta                      -> SkimByID         <PFJet>          ( PFJET_LOOSE );    
    CollectionPtr c_pfjet_highEta_ID_noMuonOverlap   = c_pfjet_highEta_ID                   -> SkimByVetoDRMatch<PFJet, Muon>    ( c_muon_final_ptCut   , jet_muon_DeltaRCut  );
    CollectionPtr c_pfjet_highEta_ID_noLeptonOverlap = c_pfjet_highEta_ID_noMuonOverlap     -> SkimByVetoDRMatch<PFJet, Electron>( c_ele_final_ptCut    , jet_ele_DeltaRCut );
    CollectionPtr c_pfjet_highEta_final              = c_pfjet_highEta_ID_noLeptonOverlap;
    CollectionPtr c_pfjet_highEta_final_ptCut        = c_pfjet_highEta_final                -> SkimByMinPt      <PFJet>          ( jet_PtCut );
    
    //-----------------------------------------------------------------
    // Get ready to fill variables 
    //-----------------------------------------------------------------

    resetCuts();

    //-----------------------------------------------------------------
    // Fill your single-object variables with values
    //-----------------------------------------------------------------
    
    fillVariableWithValue( "isData"   , isData()   );
    //fillVariableWithValue( "bunch"    , bunch      );
    fillVariableWithValue( "event"    , event      );
    fillVariableWithValue( "ls"       , luminosityBlock         );
    //fillVariableWithValue( "orbit"    , orbit      );
    fillVariableWithValue( "run"      , run        );
    //fillVariableWithValue( "ProcessID", ProcessID  );
    //fillVariableWithValue( "PtHat"    , PtHat      );
    // if amcNLOWeight filled, use it _instead_ of the nominal weight
    //FIXME
    //fillVariableWithValue( "Weight"   , fabs(amcNLOWeight)==1 ? amcNLOWeight : Weight   );
    fillVariableWithValue( "Weight"   , genWeight   );
    //FIXME
    //fillVariableWithValue( "TopPtWeight",GenParticleTopPtWeight);

    //-----------------------------------------------------------------
    // Pass JSON
    //-----------------------------------------------------------------
    fillVariableWithValue("PassJSON"                   , passJSON(run, luminosityBlock, isData())                       );    

    //-----------------------------------------------------------------
    // Fill MET filter values
    // https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    //-----------------------------------------------------------------
    fillVariableWithValue("PassGlobalTightHalo2016Filter" , int(Flag_globalTightHalo2016Filter)         == 1);
    fillVariableWithValue("PassGoodVertices"              , int(Flag_goodVertices)                      == 1);
    fillVariableWithValue("PassHBHENoiseFilter"           , int(Flag_HBHENoiseFilter                    == 1));
    fillVariableWithValue("PassHBHENoiseIsoFilter"        , int(Flag_HBHENoiseIsoFilter                 == 1));
    fillVariableWithValue("PassBadEESupercrystalFilter"   , int(Flag_eeBadScFilter                      == 1));
    fillVariableWithValue("PassEcalDeadCellTrigPrim"      , int(Flag_EcalDeadCellTriggerPrimitiveFilter == 1));
    fillVariableWithValue("PassChargedCandidateFilter"    , int(Flag_BadChargedCandidateFilter)         == 1);
    fillVariableWithValue("PassBadPFMuonFilter"           , int(Flag_BadPFMuonFilter)                   == 1);

    //-----------------------------------------------------------------
    // Fill MET values
    //-----------------------------------------------------------------
    
    //fillVariableWithValue("PFMET_Raw_Pt"       , PFMETRaw            -> at (0));      
    //fillVariableWithValue("PFMET_Raw_Phi"	     , PFMETPhiRaw	       -> at (0));
    fillVariableWithValue("PFMET_Type1_Pt"     , MET_pt);      
    fillVariableWithValue("PFMET_Type1_Phi"    , MET_phi);
    //fillVariableWithValue("PFMET_Type1XY_Pt"   , PFMETType1XYCor    -> at (0));      
    //fillVariableWithValue("PFMET_Type1XY_Phi"  , PFMETPhiType1XYCor -> at (0));
    
    if ( !isData() ) { 
      if ( reducedSkimType != 0 ){ 
        fillVariableWithValue("GenMET_Pt"		, GenMET_pt);
        fillVariableWithValue("GenMET_Phi"       	, GenMET_phi);
      }
    }

    //-----------------------------------------------------------------
    // Fill pileup variables
    //-----------------------------------------------------------------
        
    fillVariableWithValue( "nPileUpInt_BXminus1", -1 );
    fillVariableWithValue( "nPileUpInt_BX0"     , -1 );
    fillVariableWithValue( "nPileUpInt_BXplus1" , -1 );
    fillVariableWithValue( "nVertex", PV_npvs ) ;
    
    //FIXME
    //if ( !isData() ){
    //  for(int pu=0; pu<PileUpInteractions->size(); pu++) {
    //    if(PileUpOriginBX->at(pu) == 0  ) { 
    //      fillVariableWithValue( "nPileUpInt_BX0" , PileUpInteractions    ->at(pu));
    //      fillVariableWithValue( "nPileUpInt_True", PileUpInteractionsTrue->at(pu));
    //    }
    //    if(PileUpOriginBX->at(pu) == -1 ) fillVariableWithValue( "nPileUpInt_BXminus1", PileUpInteractions->at(pu));
    //    if(PileUpOriginBX->at(pu) == 1  ) fillVariableWithValue( "nPileUpInt_BXplus1" , PileUpInteractions->at(pu));
    //  }
    //}

    //-----------------------------------------------------------------
    // How many ID'd objects are there?
    //-----------------------------------------------------------------

    //int n_muonLoose          = c_muon_eta_IDLoose            -> GetSize();
    int n_muonHighPt         = c_muon_eta_IDHighPt           -> GetSize();
    int n_muon_store         = c_muon_final                  -> GetSize();
    int n_ele_store          = c_ele_final                   -> GetSize();
    int n_jet_store          = c_pfjet_final                 -> GetSize();
    int n_jet_highEta_store  = c_pfjet_highEta_final         -> GetSize();
    int n_genEle_store       = c_genEle_final                -> GetSize();
    int n_genNu_store        = c_genNu_final                 -> GetSize();
    int n_genMu_store        = c_genMu_final                 -> GetSize();
    int n_genJet_store       = c_genJet_final                -> GetSize();
						             
    int n_muon_ptCut         = c_muon_final_ptCut            -> GetSize();
    int n_ele_ptCut          = c_ele_final_ptCut             -> GetSize();
    int n_jet_ptCut          = c_pfjet_final_ptCut           -> GetSize();
    int n_jet_highEta_ptCut  = c_pfjet_highEta_final_ptCut   -> GetSize();

    int n_genZgamma_store    = c_genZgamma_final             -> GetSize();
    int n_genW_store         = c_genW_final                  -> GetSize();
    int n_genNuFromW_store   = c_genNuFromW_final            -> GetSize();
    int n_genEleFromW_store  = c_genEleFromW_final           -> GetSize();
    int n_genEleFromDY_store = c_genEleFromDY_final          -> GetSize();
    
    //if(n_ele_ptCut < 1) {
    //  std::cout << "NO GOOD ELECTRON FOUND! " << 
    //    static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(luminosityBlock) << " " << static_cast<unsigned int>(event) << std::endl;
    //  c_ele_all->examine<Electron>("c_ele_all");
    //}
    //-----------------------------------------------------------------
    // All skims need GEN particles/jets
    //-----------------------------------------------------------------

    if ( reducedSkimType != 0 ) {

      fillVariableWithValue("nGenJet_ptCut", n_genJet_store);
      fillVariableWithValue("nGenEle_ptCut", n_genEle_store);
      fillVariableWithValue("nGenNu_ptCut",  n_genNu_store);
      fillVariableWithValue("nGenMu_ptCut",  n_genMu_store);

      fillVariableWithValue("nGenJet_store", min(n_genJet_store,5));
      fillVariableWithValue("nGenEle_store", min(n_genEle_store,2));
      fillVariableWithValue("nGenNu_store" , min(n_genNu_store,2));
      fillVariableWithValue("nGenMu_store" , min(n_genMu_store,3));

      fillVariableWithValue("nGenW_ptCut"	 , n_genW_store         );
      fillVariableWithValue("nGenDY_ptCut"	 , n_genZgamma_store    );  	     
      fillVariableWithValue("nGenNuFromW_ptCut"	 , n_genNuFromW_store   );
      fillVariableWithValue("nGenEleFromW_ptCut" , n_genEleFromW_store  );
      fillVariableWithValue("nGenEleFromDY_ptCut", n_genEleFromDY_store );

      fillVariableWithValue("nGenW_store"	 , min(n_genW_store        ,2));
      fillVariableWithValue("nGenDY_store"	 , min(n_genZgamma_store   ,2));  	     
      fillVariableWithValue("nGenNuFromW_store"	 , min(n_genNuFromW_store  ,2));
      fillVariableWithValue("nGenEleFromW_store" , min(n_genEleFromW_store ,2));
      fillVariableWithValue("nGenEleFromDY_store", min(n_genEleFromDY_store,4));

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
              else {
                std::cout << "This event: " <<
                  static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(luminosityBlock) << " " << static_cast<unsigned int>(event) <<
                  " had no 5th gen Jet; examine other GenJets." << std::endl;
                c_genJet_final->examine<GenJet>("finalGenJets");
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

      // neutrinos
      if ( n_genNu_store >= 1 ){ 
        GenParticle genNu1 = c_genNu_final -> GetConstituent<GenParticle>(0);
        fillVariableWithValue ( "GenNu1_Pt" , genNu1.Pt () );
        fillVariableWithValue ( "GenNu1_Eta", genNu1.Eta() );
        fillVariableWithValue ( "GenNu1_Phi", genNu1.Phi() );

        if ( n_genNu_store >= 2 ){ 
          GenParticle genNu2 = c_genNu_final -> GetConstituent<GenParticle>(1);
          fillVariableWithValue ( "GenNu2_Pt" , genNu2.Pt () );
          fillVariableWithValue ( "GenNu2_Eta", genNu2.Eta() );
          fillVariableWithValue ( "GenNu2_Phi", genNu2.Phi() );
        }
      }

      // muons
      if ( n_genMu_store >= 1 ){ 
        GenParticle genMu1 = c_genMu_final -> GetConstituent<GenParticle>(0);
        fillVariableWithValue ( "GenMu1_Pt" , genMu1.Pt () );
        fillVariableWithValue ( "GenMu1_Eta", genMu1.Eta() );
        fillVariableWithValue ( "GenMu1_Phi", genMu1.Phi() );

        if ( n_genMu_store >= 2 ){ 
          GenParticle genMu2 = c_genMu_final -> GetConstituent<GenParticle>(1);
          fillVariableWithValue ( "GenMu2_Pt" , genMu2.Pt () );
          fillVariableWithValue ( "GenMu2_Eta", genMu2.Eta() );
          fillVariableWithValue ( "GenMu2_Phi", genMu2.Phi() );
        }

        if ( n_genMu_store >= 3 ){ 
          GenParticle genMu3 = c_genMu_final -> GetConstituent<GenParticle>(2);
          fillVariableWithValue ( "GenMu3_Pt" , genMu3.Pt () );
          fillVariableWithValue ( "GenMu3_Eta", genMu3.Eta() );
          fillVariableWithValue ( "GenMu3_Phi", genMu3.Phi() );
        }
      }

      if ( n_genZgamma_store >= 1 ){ 
        GenParticle genZgamma1 = c_genZgamma_final -> GetConstituent<GenParticle>(0);
        fillVariableWithValue ( "GenZGamma1_Pt" , genZgamma1.Pt () );
        fillVariableWithValue ( "GenZGamma1_Eta", genZgamma1.Eta() );
        fillVariableWithValue ( "GenZGamma1_Phi", genZgamma1.Phi() );
        fillVariableWithValue ( "GenZGamma1_ID" , genZgamma1.PdgId());

        if ( n_genZgamma_store >= 2 ){ 
          GenParticle genZgamma2 = c_genZgamma_final -> GetConstituent<GenParticle>(1);
          fillVariableWithValue ( "GenZGamma2_Pt" , genZgamma2.Pt () );
          fillVariableWithValue ( "GenZGamma2_Eta", genZgamma2.Eta() );
          fillVariableWithValue ( "GenZGamma2_Phi", genZgamma2.Phi() );
          fillVariableWithValue ( "GenZGamma2_ID" , genZgamma2.PdgId());
        }
      }

      if ( n_genW_store >= 1 ){ 
        GenParticle genW1 = c_genW_final -> GetConstituent<GenParticle>(0);
        fillVariableWithValue ( "GenW1_Pt" , genW1.Pt () );
        fillVariableWithValue ( "GenW1_Eta", genW1.Eta() );
        fillVariableWithValue ( "GenW1_Phi", genW1.Phi() );
        fillVariableWithValue ( "GenW1_ID" , genW1.PdgId());

        if ( n_genW_store >= 2 ){ 
          GenParticle genW2 = c_genW_final -> GetConstituent<GenParticle>(1);
          fillVariableWithValue ( "GenW2_Pt" , genW2.Pt () );
          fillVariableWithValue ( "GenW2_Eta", genW2.Eta() );
          fillVariableWithValue ( "GenW2_Phi", genW2.Phi() );
          fillVariableWithValue ( "GenW2_ID" , genW2.PdgId());
        }
      }

      if ( n_genNuFromW_store >= 1 ){ 
        GenParticle genNuFromW1 = c_genNuFromW_final -> GetConstituent<GenParticle>(0);
        fillVariableWithValue ( "GenNuFromW1_Pt" , genNuFromW1.Pt () );
        fillVariableWithValue ( "GenNuFromW1_Eta", genNuFromW1.Eta() );
        fillVariableWithValue ( "GenNuFromW1_Phi", genNuFromW1.Phi() );
        fillVariableWithValue ( "GenNuFromW1_ID" , genNuFromW1.PdgId());

        if ( n_genNuFromW_store >= 2 ){ 
          GenParticle genNuFromW2 = c_genNuFromW_final -> GetConstituent<GenParticle>(1);
          fillVariableWithValue ( "GenNuFromW2_Pt" , genNuFromW2.Pt () );
          fillVariableWithValue ( "GenNuFromW2_Eta", genNuFromW2.Eta() );
          fillVariableWithValue ( "GenNuFromW2_Phi", genNuFromW2.Phi() );
          fillVariableWithValue ( "GenNuFromW2_ID" , genNuFromW2.PdgId());
        }
      }

      if ( n_genEleFromW_store >= 1 ){ 
        GenParticle genEleFromW1 = c_genEleFromW_final -> GetConstituent<GenParticle>(0);
        fillVariableWithValue ( "GenEleFromW1_Pt" , genEleFromW1.Pt () );
        fillVariableWithValue ( "GenEleFromW1_Eta", genEleFromW1.Eta() );
        fillVariableWithValue ( "GenEleFromW1_Phi", genEleFromW1.Phi() );
        fillVariableWithValue ( "GenEleFromW1_ID" , genEleFromW1.PdgId());

        if ( n_genEleFromW_store >= 2 ){ 
          GenParticle genEleFromW2 = c_genEleFromW_final -> GetConstituent<GenParticle>(1);
          fillVariableWithValue ( "GenEleFromW2_Pt" , genEleFromW2.Pt () );
          fillVariableWithValue ( "GenEleFromW2_Eta", genEleFromW2.Eta() );
          fillVariableWithValue ( "GenEleFromW2_Phi", genEleFromW2.Phi() );
          fillVariableWithValue ( "GenEleFromW2_ID" , genEleFromW2.PdgId());
        }
      }

      if ( n_genEleFromDY_store >= 1 ){ 
        GenParticle genEleFromDY1 = c_genEleFromDY_final -> GetConstituent<GenParticle>(0);
        fillVariableWithValue ( "GenEleFromDY1_Pt" , genEleFromDY1.Pt () );
        fillVariableWithValue ( "GenEleFromDY1_Eta", genEleFromDY1.Eta() );
        fillVariableWithValue ( "GenEleFromDY1_Phi", genEleFromDY1.Phi() );
        fillVariableWithValue ( "GenEleFromDY1_ID" , genEleFromDY1.PdgId());

        if ( n_genEleFromDY_store >= 2 ){ 
          GenParticle genEleFromDY2 = c_genEleFromDY_final -> GetConstituent<GenParticle>(1);
          fillVariableWithValue ( "GenEleFromDY2_Pt" , genEleFromDY2.Pt () );
          fillVariableWithValue ( "GenEleFromDY2_Eta", genEleFromDY2.Eta() );
          fillVariableWithValue ( "GenEleFromDY2_Phi", genEleFromDY2.Phi() );
          fillVariableWithValue ( "GenEleFromDY2_ID" , genEleFromDY2.PdgId());

          if ( n_genEleFromDY_store >= 3 ){ 
            GenParticle genEleFromDY3 = c_genEleFromDY_final -> GetConstituent<GenParticle>(2);
            fillVariableWithValue ( "GenEleFromDY3_Pt" , genEleFromDY3.Pt () );
            fillVariableWithValue ( "GenEleFromDY3_Eta", genEleFromDY3.Eta() );
            fillVariableWithValue ( "GenEleFromDY3_Phi", genEleFromDY3.Phi() );
            fillVariableWithValue ( "GenEleFromDY3_ID" , genEleFromDY3.PdgId());

            if ( n_genEleFromDY_store >= 4 ){ 
              GenParticle genEleFromDY4 = c_genEleFromDY_final -> GetConstituent<GenParticle>(3);
              fillVariableWithValue ( "GenEleFromDY4_Pt" , genEleFromDY4.Pt () );
              fillVariableWithValue ( "GenEleFromDY4_Eta", genEleFromDY4.Eta() );
              fillVariableWithValue ( "GenEleFromDY4_Phi", genEleFromDY4.Phi() );
              fillVariableWithValue ( "GenEleFromDY4_ID" , genEleFromDY4.PdgId());
            }
          }
        }
      }
    }
    
    //-----------------------------------------------------------------
    // All skims need muons
    //-----------------------------------------------------------------

    fillVariableWithValue ("nMuon_ptCut", n_muon_ptCut);
    //fillVariableWithValue ("nMuon_LooseId", n_muonLoose);
    fillVariableWithValue ("nMuon_HighPtId", n_muonHighPt);
    fillVariableWithValue ("nMuon_store", min(n_muon_store,3));

    if ( n_muon_store >= 1 ){ 

      Muon muon1 = c_muon_final -> GetConstituent<Muon>(0);
      //double hltTTMuon1Pt     = triggerMatchPt<HLTFilterObject, Muon>(c_hltMuon22_TTbar_all , muon1, muon_hltMatch_DeltaRCut);
      double hltSingleMuon1Pt = triggerMatchPt<HLTriggerObject, Muon>(c_hltMuon_SingleMu_all, muon1, muon_hltMatch_DeltaRCut);
      fillVariableWithValue ("Muon1_Pt"             , muon1.Pt      ());
      fillVariableWithValue ("Muon1_Eta"            , muon1.Eta     ());
      fillVariableWithValue ("Muon1_Phi"            , muon1.Phi     ());
      fillVariableWithValue ("Muon1_PtError"        , muon1.PtError ());
      fillVariableWithValue ("Muon1_Charge"         , muon1.Charge  ());
      fillVariableWithValue ("Muon1_hltSingleMuonPt", hltSingleMuon1Pt);
      
      if ( n_muon_store >= 2 ){ 

        Muon muon2 = c_muon_final -> GetConstituent<Muon>(1);
        //double hltMuon2Pt       = triggerMatchPt<HLTFilterObject, Muon>(c_hltMuon22_TTbar_all , muon2, muon_hltMatch_DeltaRCut);
        double hltSingleMuon2Pt = triggerMatchPt<HLTriggerObject, Muon>(c_hltMuon_SingleMu_all, muon2, muon_hltMatch_DeltaRCut);
        fillVariableWithValue ("Muon2_Pt"             , muon2.Pt      ());
        fillVariableWithValue ("Muon2_Eta"            , muon2.Eta     ());
        fillVariableWithValue ("Muon2_Phi"            , muon2.Phi     ());
        fillVariableWithValue ("Muon2_PtError"        , muon2.PtError ());
        fillVariableWithValue ("Muon2_Charge"         , muon2.Charge  ());
        fillVariableWithValue ("Muon2_hltSingleMuonPt", hltSingleMuon2Pt);

        if ( n_muon_store >= 3 ){ 

          Muon muon3 = c_muon_final -> GetConstituent<Muon>(2);
          //double hltMuon3Pt       = triggerMatchPt<HLTFilterObject, Muon>(c_hltMuon22_TTbar_all , muon3, muon_hltMatch_DeltaRCut);
          double hltSingleMuon3Pt = triggerMatchPt<HLTriggerObject, Muon>(c_hltMuon_SingleMu_all, muon3, muon_hltMatch_DeltaRCut);
          fillVariableWithValue ("Muon3_Pt"             , muon3.Pt      ());
          fillVariableWithValue ("Muon3_Eta"            , muon3.Eta     ());
          fillVariableWithValue ("Muon3_Phi"            , muon3.Phi     ());
          fillVariableWithValue ("Muon3_PtError"        , muon3.PtError ());
          fillVariableWithValue ("Muon3_Charge"         , muon3.Charge  ());
          fillVariableWithValue ("Muon3_hltSingleMuonPt", hltSingleMuon3Pt);
        }
      }
    }

    //-----------------------------------------------------------------
    // Fill variables for QCD skim (reducedSkimType == 0 )
    //-----------------------------------------------------------------

    if ( reducedSkimType == 0 ) { 

      fillVariableWithValue ("nLooseEle_store"   , min(n_ele_store,3));
      fillVariableWithValue ("nLooseEle_ID"      , n_ele_store);
      fillVariableWithValue ("nJetLooseEle_store", min(n_jet_store,5));
      fillVariableWithValue ("nJetLooseEle_etaIdLepCleaned", n_jet_store);
      fillVariableWithValue ("nLooseEle_ptCut"   , n_ele_ptCut );
      fillVariableWithValue ("nJetLooseEle_ptCut", n_jet_ptCut );

      int n_filters = (
          c_hltPhoton22_QCD_all -> GetSize() + 
          c_hltPhoton30_QCD_all -> GetSize() + 
          c_hltPhoton36_QCD_all -> GetSize() + 
          c_hltPhoton50_QCD_all -> GetSize() + 
          c_hltPhoton75_QCD_all -> GetSize() + 
          c_hltPhoton90_QCD_all -> GetSize() + 
          c_hltPhoton120_QCD_all         -> GetSize() + 
          c_hltPhoton175_QCD_all         -> GetSize()
          );

      if ( n_ele_store >= 1 ){
        Electron loose_ele1 = c_ele_final -> GetConstituent<Electron>(0);

        double hltPhotonPt = -999.;
        if ( n_filters != 0 ) {
          double hltPhotonPt_array [8] =  {
            triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton22_QCD_all, loose_ele1, ele_hltMatch_DeltaRCut),
            triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton30_QCD_all, loose_ele1, ele_hltMatch_DeltaRCut),
            triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton36_QCD_all, loose_ele1, ele_hltMatch_DeltaRCut),
            triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton50_QCD_all, loose_ele1, ele_hltMatch_DeltaRCut), 
            triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton75_QCD_all, loose_ele1, ele_hltMatch_DeltaRCut), 
            triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton90_QCD_all, loose_ele1, ele_hltMatch_DeltaRCut), 
            triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton120_QCD_all        , loose_ele1, ele_hltMatch_DeltaRCut), 
            triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton175_QCD_all        , loose_ele1, ele_hltMatch_DeltaRCut)
          };
          for (int iFilter = 0; iFilter < 8; ++iFilter){
            if ( hltPhotonPt_array[iFilter] > 0.0 ) hltPhotonPt = hltPhotonPt_array[iFilter];
          }
        }

        fillVariableWithValue( "LooseEle1_PassHEEPID"           , loose_ele1.PassUserID ( HEEP70 )  );
        fillVariableWithValue( "LooseEle1_Pt"                   , loose_ele1.Pt()                 );
        fillVariableWithValue( "LooseEle1_PtHeep"               , loose_ele1.PtHeep()                 );
        fillVariableWithValue( "LooseEle1_EcalDriven"           , loose_ele1.EcalSeed()           );
        fillVariableWithValue( "LooseEle1_DeltaEtaSeed"         , loose_ele1.DeltaEtaSeed()       );
        fillVariableWithValue( "LooseEle1_SCEnergy"             , loose_ele1.SCEnergy()           );
        fillVariableWithValue( "LooseEle1_SCEt"                 , loose_ele1.SCEnergy()/cosh(loose_ele1.SCEta()) );
        fillVariableWithValue( "LooseEle1_Full5x5E1x5OverE5x5"  , loose_ele1.Full5x5E1x5OverE5x5());
        fillVariableWithValue( "LooseEle1_Full5x5E2x5OverE5x5"  , loose_ele1.Full5x5E2x5OverE5x5());
        fillVariableWithValue( "LooseEle1_RhoForHeep"           , loose_ele1.RhoForHEEP());
        fillVariableWithValue( "LooseEle1_Energy"               , loose_ele1.CaloEnergy()         );
        fillVariableWithValue( "LooseEle1_Eta"                  , loose_ele1.Eta()                );
        fillVariableWithValue( "LooseEle1_Phi"                  , loose_ele1.Phi()                );
        fillVariableWithValue( "LooseEle1_SCEta"                , loose_ele1.SCEta()              );
        fillVariableWithValue( "LooseEle1_SCPhi"                , loose_ele1.SCPhi()              );
        fillVariableWithValue( "LooseEle1_Charge"               , loose_ele1.Charge()             );
        fillVariableWithValue( "LooseEle1_Dist"                 , loose_ele1.Dist()               );
        fillVariableWithValue( "LooseEle1_DCotTheta"            , loose_ele1.DCotTheta()          );
        fillVariableWithValue( "LooseEle1_MissingHits"          , loose_ele1.MissingHits()        );
        fillVariableWithValue( "LooseEle1_TrkPt"                , loose_ele1.TrackPt()            );
        fillVariableWithValue( "LooseEle1_TrkEta"               , loose_ele1.TrackEta()           );
        fillVariableWithValue( "LooseEle1_Full5x5SigmaIEtaIEta" , loose_ele1.Full5x5SigmaIEtaIEta());
        fillVariableWithValue( "LooseEle1_SigmaEtaEta"          , loose_ele1.SigmaEtaEta()        );

        fillVariableWithValue( "LooseEle1_DeltaPhiTrkSC" , loose_ele1.DeltaPhi()           );
        fillVariableWithValue( "LooseEle1_DeltaEtaTrkSC" , loose_ele1.DeltaEta()           );
        fillVariableWithValue( "LooseEle1_RawEnergy"     , loose_ele1.RawEnergy()          );
        fillVariableWithValue( "LooseEle1_NBrems"        , loose_ele1.NBrems()             );
        fillVariableWithValue( "LooseEle1_HoE"           , loose_ele1.HoE()                );
        fillVariableWithValue( "LooseEle1_HasMatchedPhot", loose_ele1.HasMatchedConvPhot() );
        fillVariableWithValue( "LooseEle1_FBrem"         , loose_ele1.FBrem()              );
        fillVariableWithValue( "LooseEle1_LeadVtxDistXY" , loose_ele1.LeadVtxDistXY()      );
        fillVariableWithValue( "LooseEle1_LeadVtxDistZ"  , loose_ele1.LeadVtxDistZ ()      );
        fillVariableWithValue( "LooseEle1_BeamSpotDXY"   , loose_ele1.BeamSpotDXY()        );
        fillVariableWithValue( "LooseEle1_BeamSpotDXYErr", loose_ele1.BeamSpotDXYErr()     );
        fillVariableWithValue( "LooseEle1_ValidFrac"     , loose_ele1.ValidFrac()          );
        fillVariableWithValue( "LooseEle1_Classif"       , loose_ele1.Classif()            );
        fillVariableWithValue( "LooseEle1_EOverP"        , loose_ele1.ESuperClusterOverP() );

        fillVariableWithValue( "LooseEle1_TrkIsolation"  , loose_ele1.TrkIsoDR03()          );
        fillVariableWithValue( "LooseEle1_TrkIsoHEEP7"   , loose_ele1.HEEP70TrackIsolation());
        fillVariableWithValue( "LooseEle1_EcalIsolation" , loose_ele1.EcalIsoDR03()         );
        fillVariableWithValue( "LooseEle1_HcalIsolation" , loose_ele1.HcalIsoD1DR03()       );
        fillVariableWithValue( "LooseEle1_CorrIsolation" , loose_ele1.HEEPCorrIsolation()   );
        fillVariableWithValue( "LooseEle1_PFCHIso03"     , loose_ele1.PFChargedHadronIso03());
        fillVariableWithValue( "LooseEle1_PFPhoIso03"    , loose_ele1.PFPhotonIso03       ());
        fillVariableWithValue( "LooseEle1_PFNHIso03"     , loose_ele1.PFNeutralHadronIso03());

        fillVariableWithValue("LooseEle1_GsfCtfScPixCharge", loose_ele1.GsfCtfScPixCharge()  );
        fillVariableWithValue("LooseEle1_GsfScPixCharge"   , loose_ele1.GsfScPixCharge()     );
        fillVariableWithValue("LooseEle1_GsfCtfCharge"     , loose_ele1.GsfCtfCharge()       );

        fillVariableWithValue( "LooseEle1_hltPhotonPt"  , hltPhotonPt );

        if ( n_ele_store >= 2 ){

          Electron loose_ele2 = c_ele_final -> GetConstituent<Electron>(1);
          hltPhotonPt = -999.;
          if ( n_filters != 0 ) {
            double hltPhotonPt_array [8] =  {
              triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton22_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut),
              triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton30_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut),
              triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton36_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut),
              triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton50_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut), 
              triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton75_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut), 
              triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton90_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut), 
              triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton120_QCD_all        , loose_ele2, ele_hltMatch_DeltaRCut), 
              triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton175_QCD_all        , loose_ele2, ele_hltMatch_DeltaRCut)
            };
            for (int iFilter = 0; iFilter < 8; ++iFilter){
              if ( hltPhotonPt_array[iFilter] > 0.0 ) hltPhotonPt = hltPhotonPt_array[iFilter];
            }
          }

          fillVariableWithValue( "LooseEle2_PassHEEPID"    , loose_ele2.PassUserID ( HEEP70 ));
          fillVariableWithValue( "LooseEle2_Pt"            , loose_ele2.Pt()                 );
          fillVariableWithValue( "LooseEle2_PtHeep"        , loose_ele2.PtHeep()             );
          fillVariableWithValue( "LooseEle2_EcalDriven"    , loose_ele2.EcalSeed()           );
          fillVariableWithValue( "LooseEle2_DeltaEtaSeed"  , loose_ele2.DeltaEtaSeed()       );
          fillVariableWithValue( "LooseEle2_SCEnergy"      , loose_ele2.SCEnergy()           );
          fillVariableWithValue( "LooseEle2_SCEt"                 , loose_ele2.SCEnergy()/cosh(loose_ele2.SCEta()) );
          fillVariableWithValue( "LooseEle2_Full5x5E1x5OverE5x5"  , loose_ele2.Full5x5E1x5OverE5x5());
          fillVariableWithValue( "LooseEle2_Full5x5E2x5OverE5x5"  , loose_ele2.Full5x5E2x5OverE5x5());
          fillVariableWithValue( "LooseEle2_RhoForHeep"  , loose_ele2.RhoForHEEP());
          fillVariableWithValue( "LooseEle2_Energy"        , loose_ele2.CaloEnergy()         );
          fillVariableWithValue( "LooseEle2_Eta"           , loose_ele2.Eta()                );
          fillVariableWithValue( "LooseEle2_Phi"           , loose_ele2.Phi()                );
          fillVariableWithValue( "LooseEle2_SCEta"         , loose_ele2.SCEta()              );
          fillVariableWithValue( "LooseEle2_SCPhi"         , loose_ele2.SCPhi()              );
          fillVariableWithValue( "LooseEle2_Charge"        , loose_ele2.Charge()             );
          fillVariableWithValue( "LooseEle2_Dist"          , loose_ele2.Dist()               );
          fillVariableWithValue( "LooseEle2_DCotTheta"     , loose_ele2.DCotTheta()          );
          fillVariableWithValue( "LooseEle2_MissingHits"   , loose_ele2.MissingHits()        );
          fillVariableWithValue( "LooseEle2_TrkPt"         , loose_ele2.TrackPt()            );
          fillVariableWithValue( "LooseEle2_TrkEta"        , loose_ele2.TrackEta()           );
          fillVariableWithValue( "LooseEle2_Full5x5SigmaIEtaIEta" , loose_ele2.Full5x5SigmaIEtaIEta());
          fillVariableWithValue( "LooseEle2_SigmaEtaEta"   , loose_ele2.SigmaEtaEta()        );

          fillVariableWithValue( "LooseEle2_DeltaPhiTrkSC" , loose_ele2.DeltaPhi()           );
          fillVariableWithValue( "LooseEle2_DeltaEtaTrkSC" , loose_ele2.DeltaEta()           );
          fillVariableWithValue( "LooseEle2_RawEnergy"     , loose_ele2.RawEnergy()          );
          fillVariableWithValue( "LooseEle2_NBrems"        , loose_ele2.NBrems()             );
          fillVariableWithValue( "LooseEle2_HoE"           , loose_ele2.HoE()                );
          fillVariableWithValue( "LooseEle2_HasMatchedPhot", loose_ele2.HasMatchedConvPhot() );
          fillVariableWithValue( "LooseEle2_FBrem"         , loose_ele2.FBrem()              );
          fillVariableWithValue( "LooseEle2_LeadVtxDistXY" , loose_ele2.LeadVtxDistXY()      );
          fillVariableWithValue( "LooseEle2_LeadVtxDistZ"  , loose_ele2.LeadVtxDistZ ()      );
          fillVariableWithValue( "LooseEle2_BeamSpotDXY"   , loose_ele2.BeamSpotDXY()        );
          fillVariableWithValue( "LooseEle2_BeamSpotDXYErr", loose_ele2.BeamSpotDXYErr()     );
          fillVariableWithValue( "LooseEle2_ValidFrac"     , loose_ele2.ValidFrac()          );
          fillVariableWithValue( "LooseEle2_Classif"       , loose_ele2.Classif()            );
          fillVariableWithValue( "LooseEle2_EOverP"        , loose_ele2.ESuperClusterOverP() );

          fillVariableWithValue( "LooseEle2_TrkIsolation"  , loose_ele2.TrkIsoDR03()          );
          fillVariableWithValue( "LooseEle2_TrkIsoHEEP7"   , loose_ele2.HEEP70TrackIsolation());
          fillVariableWithValue( "LooseEle2_EcalIsolation" , loose_ele2.EcalIsoDR03()         );
          fillVariableWithValue( "LooseEle2_HcalIsolation" , loose_ele2.HcalIsoD1DR03()       );
          fillVariableWithValue( "LooseEle2_CorrIsolation" , loose_ele2.HEEPCorrIsolation()   );
          fillVariableWithValue( "LooseEle2_PFCHIso03"     , loose_ele2.PFChargedHadronIso03());
          fillVariableWithValue( "LooseEle2_PFPhoIso03"    , loose_ele2.PFPhotonIso03       ());
          fillVariableWithValue( "LooseEle2_PFNHIso03"     , loose_ele2.PFNeutralHadronIso03());

          fillVariableWithValue("LooseEle2_GsfCtfScPixCharge", loose_ele2.GsfCtfScPixCharge()  );
          fillVariableWithValue("LooseEle2_GsfScPixCharge"   , loose_ele2.GsfScPixCharge()     );
          fillVariableWithValue("LooseEle2_GsfCtfCharge"     , loose_ele2.GsfCtfCharge()       );

          fillVariableWithValue( "LooseEle2_hltPhotonPt"  , hltPhotonPt );

          if ( n_ele_store >= 3 ){
            Electron loose_ele3 = c_ele_final -> GetConstituent<Electron>(2);
            hltPhotonPt = -999.;
            if ( n_filters != 0 ) {
              double hltPhotonPt_array [8] =  {
                triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton22_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut),
                triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton30_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut),
                triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton36_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut),
                triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton50_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut), 
                triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton75_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut), 
                triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton90_QCD_all, loose_ele2, ele_hltMatch_DeltaRCut), 
                triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton120_QCD_all        , loose_ele2, ele_hltMatch_DeltaRCut), 
                triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton175_QCD_all        , loose_ele2, ele_hltMatch_DeltaRCut)
              };
              for (int iFilter = 0; iFilter < 8; ++iFilter){
                if ( hltPhotonPt_array[iFilter] > 0.0 ) hltPhotonPt = hltPhotonPt_array[iFilter];
              }
            }

            fillVariableWithValue( "LooseEle3_PassHEEPID"    , loose_ele3.PassUserID ( HEEP70 )  );
            fillVariableWithValue( "LooseEle3_Pt"            , loose_ele3.Pt()                 );
            fillVariableWithValue( "LooseEle3_PtHeep"        , loose_ele3.PtHeep()                 );
            fillVariableWithValue( "LooseEle3_EcalDriven"           , loose_ele3.EcalSeed()           );
            fillVariableWithValue( "LooseEle3_DeltaEtaSeed"         , loose_ele3.DeltaEtaSeed()       );
            fillVariableWithValue( "LooseEle3_SCEnergy"             , loose_ele3.SCEnergy()           );
            fillVariableWithValue( "LooseEle3_SCEt"                 , loose_ele3.SCEnergy()/cosh(loose_ele3.SCEta()) );
            fillVariableWithValue( "LooseEle3_Full5x5E1x5OverE5x5"  , loose_ele3.Full5x5E1x5OverE5x5());
            fillVariableWithValue( "LooseEle3_Full5x5E2x5OverE5x5"  , loose_ele3.Full5x5E2x5OverE5x5());
            fillVariableWithValue( "LooseEle3_RhoForHeep"           , loose_ele3.RhoForHEEP());
            fillVariableWithValue( "LooseEle3_Energy"        , loose_ele3.CaloEnergy()         );
            fillVariableWithValue( "LooseEle3_Eta"           , loose_ele3.Eta()                );
            fillVariableWithValue( "LooseEle3_Phi"           , loose_ele3.Phi()                );
            fillVariableWithValue( "LooseEle3_SCEta"         , loose_ele3.SCEta()              );
            fillVariableWithValue( "LooseEle3_SCPhi"         , loose_ele3.SCPhi()              );
            fillVariableWithValue( "LooseEle3_Charge"        , loose_ele3.Charge()             );
            fillVariableWithValue( "LooseEle3_Dist"          , loose_ele3.Dist()               );
            fillVariableWithValue( "LooseEle3_DCotTheta"     , loose_ele3.DCotTheta()          );
            fillVariableWithValue( "LooseEle3_MissingHits"   , loose_ele3.MissingHits()        );
            fillVariableWithValue( "LooseEle3_TrkPt"         , loose_ele3.TrackPt()            );
            fillVariableWithValue( "LooseEle3_TrkEta"        , loose_ele3.TrackEta()           );
            fillVariableWithValue( "LooseEle3_Full5x5SigmaIEtaIEta" , loose_ele3.Full5x5SigmaIEtaIEta());
            fillVariableWithValue( "LooseEle3_SigmaEtaEta"   , loose_ele3.SigmaEtaEta()        );

            fillVariableWithValue( "LooseEle3_DeltaPhiTrkSC" , loose_ele3.DeltaPhi()           );
            fillVariableWithValue( "LooseEle3_DeltaEtaTrkSC" , loose_ele3.DeltaEta()           );
            fillVariableWithValue( "LooseEle3_RawEnergy"     , loose_ele3.RawEnergy()          );
            fillVariableWithValue( "LooseEle3_NBrems"        , loose_ele3.NBrems()             );
            fillVariableWithValue( "LooseEle3_HoE"           , loose_ele3.HoE()                );
            fillVariableWithValue( "LooseEle3_HasMatchedPhot", loose_ele3.HasMatchedConvPhot() );
            fillVariableWithValue( "LooseEle3_FBrem"         , loose_ele3.FBrem()              );
            fillVariableWithValue( "LooseEle3_LeadVtxDistXY" , loose_ele3.LeadVtxDistXY()      );
            fillVariableWithValue( "LooseEle3_LeadVtxDistZ"  , loose_ele3.LeadVtxDistZ ()      );
            fillVariableWithValue( "LooseEle3_BeamSpotDXY"   , loose_ele3.BeamSpotDXY()        );
            fillVariableWithValue( "LooseEle3_BeamSpotDXYErr", loose_ele3.BeamSpotDXYErr()     );
            fillVariableWithValue( "LooseEle3_ValidFrac"     , loose_ele3.ValidFrac()          );
            fillVariableWithValue( "LooseEle3_Classif"       , loose_ele3.Classif()            );
            fillVariableWithValue( "LooseEle3_EOverP"        , loose_ele3.ESuperClusterOverP() );

            fillVariableWithValue( "LooseEle3_TrkIsolation"  , loose_ele3.TrkIsoDR03()          );
            fillVariableWithValue( "LooseEle3_TrkIsoHEEP7"   , loose_ele3.HEEP70TrackIsolation());
            fillVariableWithValue( "LooseEle3_EcalIsolation" , loose_ele3.EcalIsoDR03()         );
            fillVariableWithValue( "LooseEle3_HcalIsolation" , loose_ele3.HcalIsoD1DR03()       );
            fillVariableWithValue( "LooseEle3_CorrIsolation" , loose_ele3.HEEPCorrIsolation()   );
            fillVariableWithValue( "LooseEle3_PFCHIso03"     , loose_ele3.PFChargedHadronIso03());
            fillVariableWithValue( "LooseEle3_PFPhoIso03"    , loose_ele3.PFPhotonIso03       ());
            fillVariableWithValue( "LooseEle3_PFNHIso03"     , loose_ele3.PFNeutralHadronIso03());

            fillVariableWithValue("LooseEle3_GsfCtfScPixCharge", loose_ele3.GsfCtfScPixCharge()  );
            fillVariableWithValue("LooseEle3_GsfScPixCharge"   , loose_ele3.GsfScPixCharge()     );
            fillVariableWithValue("LooseEle3_GsfCtfCharge"     , loose_ele3.GsfCtfCharge()       );

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
        fillVariableWithValue( "JetLooseEle1_btagCISV" , loose_jet1.CombinedInclusiveSecondaryVertexBTag());
        fillVariableWithValue( "JetLooseEle1_btagCMVA" , loose_jet1.CombinedMVABTag());

        if ( n_jet_store >= 2 ){
          PFJet loose_jet2 = c_pfjet_final -> GetConstituent<PFJet>(1);
          fillVariableWithValue( "JetLooseEle2_Pt"      , loose_jet2.Pt()                         );
          fillVariableWithValue( "JetLooseEle2_JECUnc"  , loose_jet2.JECUnc()                     );
          fillVariableWithValue( "JetLooseEle2_Energy"  , loose_jet2.Energy()                     );
          fillVariableWithValue( "JetLooseEle2_Eta"     , loose_jet2.Eta()                        );
          fillVariableWithValue( "JetLooseEle2_Phi"     , loose_jet2.Phi()                        );
          fillVariableWithValue( "JetLooseEle2_btagCISV" , loose_jet2.CombinedInclusiveSecondaryVertexBTag());
          fillVariableWithValue( "JetLooseEle2_btagCMVA" , loose_jet2.CombinedMVABTag());

          if ( n_jet_store >= 3 ){
            PFJet loose_jet3 = c_pfjet_final -> GetConstituent<PFJet>(2);
            fillVariableWithValue( "JetLooseEle3_Pt"      , loose_jet3.Pt()                         );
            fillVariableWithValue( "JetLooseEle3_JECUnc"  , loose_jet3.JECUnc()                     );
            fillVariableWithValue( "JetLooseEle3_Energy"  , loose_jet3.Energy()                     );
            fillVariableWithValue( "JetLooseEle3_Eta"     , loose_jet3.Eta()                        );
            fillVariableWithValue( "JetLooseEle3_Phi"     , loose_jet3.Phi()                        );
            fillVariableWithValue( "JetLooseEle3_btagCISV" , loose_jet3.CombinedInclusiveSecondaryVertexBTag());
            fillVariableWithValue( "JetLooseEle3_btagCMVA" , loose_jet3.CombinedMVABTag());

            if ( n_jet_store >= 4 ){
              PFJet loose_jet4 = c_pfjet_final -> GetConstituent<PFJet>(3);
              fillVariableWithValue( "JetLooseEle4_Pt"      , loose_jet4.Pt()                         );
              fillVariableWithValue( "JetLooseEle4_JECUnc"  , loose_jet4.JECUnc()                     );
              fillVariableWithValue( "JetLooseEle4_Energy"  , loose_jet4.Energy()                     );
              fillVariableWithValue( "JetLooseEle4_Eta"     , loose_jet4.Eta()                        );
              fillVariableWithValue( "JetLooseEle4_Phi"     , loose_jet4.Phi()                        );
              fillVariableWithValue( "JetLooseEle4_btagCISV" , loose_jet4.CombinedInclusiveSecondaryVertexBTag());
              fillVariableWithValue( "JetLooseEle4_btagCMVA" , loose_jet4.CombinedMVABTag());

              if ( n_jet_store >= 5 ){
                PFJet loose_jet5 = c_pfjet_final -> GetConstituent<PFJet>(4);
                fillVariableWithValue( "JetLooseEle5_Pt"      , loose_jet5.Pt()                         );
                fillVariableWithValue( "JetLooseEle5_JECUnc"  , loose_jet5.JECUnc()                     );
                fillVariableWithValue( "JetLooseEle5_Energy"  , loose_jet5.Energy()                     );
                fillVariableWithValue( "JetLooseEle5_Eta"     , loose_jet5.Eta()                        );
                fillVariableWithValue( "JetLooseEle5_Phi"     , loose_jet5.Phi()                        );
                fillVariableWithValue( "JetLooseEle5_btagCISV" , loose_jet5.CombinedInclusiveSecondaryVertexBTag());
                fillVariableWithValue( "JetLooseEle5_btagCMVA" , loose_jet5.CombinedMVABTag());
              }
            }
          }
        }
      }
    }

    //-----------------------------------------------------------------
    // Fill variables for signal-like skim (reducedSkimType == 1 - 4 )
    //-----------------------------------------------------------------
    
    // Electrons

    else if ( reducedSkimType == 1 || reducedSkimType == 2 || reducedSkimType == 3 || reducedSkimType == 4) { 
      
      fillVariableWithValue ("nEle_store"       , min(n_ele_store,2) );
      fillVariableWithValue ("nEle_ID"          , n_ele_store);
      fillVariableWithValue ("nJet_store"       , min(n_jet_store,5) );
      fillVariableWithValue ("nJet_etaIdLepCleaned"       , n_jet_store);
      fillVariableWithValue ("nHighEtaJet_store", min(n_jet_highEta_store, 1));
      fillVariableWithValue ("nEle_ptCut"       , n_ele_ptCut );
      fillVariableWithValue ("nJet_ptCut"       , n_jet_ptCut );
      fillVariableWithValue ("nHighEtaJet_ptCut", n_jet_highEta_ptCut );

      if ( n_ele_store >= 1 ){
        Electron ele1 = c_ele_final -> GetConstituent<Electron>(0);
        double hltEle1Pt_signal          = triggerMatchPt<HLTriggerObject, Electron>(c_hltEle45_Signal_all    , ele1, ele_hltMatch_DeltaRCut);
        //double hltEle1Pt_ttbar           = triggerMatchPt<HLTFilterObject, Electron>(c_hltPhoton22_TTbar_all  , ele1, ele_hltMatch_DeltaRCut);
        double hltEle1Pt_doubleEleSignal = triggerMatchPt<HLTriggerObject, Electron>(c_hltDoubleEle_Signal_all, ele1, ele_hltMatch_DeltaRCut);
        double hltEle1Pt_WP80            = triggerMatchPt<HLTriggerObject, Electron>(c_hltEle27WP85Gsf_all       , ele1, ele_hltMatch_DeltaRCut);

        fillVariableWithValue( "Ele1_Pt"            , ele1.Pt()                 );
        fillVariableWithValue( "Ele1_Energy"        , ele1.CaloEnergy()         );
        fillVariableWithValue( "Ele1_Eta"           , ele1.Eta()                );
        fillVariableWithValue( "Ele1_Phi"           , ele1.Phi()                );
        fillVariableWithValue( "Ele1_PtHeep"        , ele1.PtHeep()             );
        fillVariableWithValue( "Ele1_SCEt"          , ele1.SCEnergy()/cosh(ele1.SCEta()) );
        fillVariableWithValue( "Ele1_SCEta"         , ele1.SCEta()              );
        fillVariableWithValue( "Ele1_SCPhi"         , ele1.SCPhi()              );
        fillVariableWithValue( "Ele1_Charge"        , ele1.Charge()             );
        fillVariableWithValue( "Ele1_Dist"          , ele1.Dist()               );
        fillVariableWithValue( "Ele1_DCotTheta"     , ele1.DCotTheta()          );
        fillVariableWithValue( "Ele1_MissingHits"   , ele1.MissingHits()        );
        fillVariableWithValue( "Ele1_TrkPt"         , ele1.TrackPt()            );
        fillVariableWithValue( "Ele1_TrkEta"        , ele1.TrackEta()            );
        fillVariableWithValue( "Ele1_SigmaEtaEta"   , ele1.SigmaEtaEta()        );
        fillVariableWithValue( "Ele1_EcalDriven"    , ele1.EcalDriven()         );
        fillVariableWithValue( "Ele1_DeltaEtaSeed"  , ele1.DeltaEtaSeed()       );
        fillVariableWithValue( "Ele1_SCEnergy"      , ele1.SCEnergy()           );
        fillVariableWithValue( "Ele1_Full5x5SigmaIEtaIEta" , ele1.Full5x5SigmaIEtaIEta() );
        fillVariableWithValue( "Ele1_Full5x5E1x5OverE5x5"  , ele1.Full5x5E1x5OverE5x5()  );
        fillVariableWithValue( "Ele1_Full5x5E2x5OverE5x5"  , ele1.Full5x5E2x5OverE5x5()  );
        fillVariableWithValue( "Ele1_RhoForHEEP"    , ele1.RhoForHEEP()         );

        fillVariableWithValue( "Ele1_DeltaPhiTrkSC" , ele1.DeltaPhi()           );
        fillVariableWithValue( "Ele1_DeltaEtaTrkSC" , ele1.DeltaEta()           );
        fillVariableWithValue( "Ele1_RawEnergy"     , ele1.RawEnergy()          );
        fillVariableWithValue( "Ele1_NBrems"        , ele1.NBrems()             );
        fillVariableWithValue( "Ele1_HoE"           , ele1.HoE()                );
        fillVariableWithValue( "Ele1_HasMatchedPhot", ele1.HasMatchedConvPhot() );
        fillVariableWithValue( "Ele1_FBrem"         , ele1.FBrem()              );
        fillVariableWithValue( "Ele1_LeadVtxDistXY" , ele1.LeadVtxDistXY()      );
        fillVariableWithValue( "Ele1_LeadVtxDistZ"  , ele1.LeadVtxDistZ ()      );
        fillVariableWithValue( "Ele1_BeamSpotDXY"   , ele1.BeamSpotDXY()        );
        fillVariableWithValue( "Ele1_BeamSpotDXYErr", ele1.BeamSpotDXYErr()     );
        fillVariableWithValue( "Ele1_ValidFrac"     , ele1.ValidFrac()          );
        fillVariableWithValue( "Ele1_Classif"       , ele1.Classif()            );
        fillVariableWithValue( "Ele1_EOverP"        , ele1.ESuperClusterOverP() );

        fillVariableWithValue( "Ele1_TrkIsolation"  , ele1.TrkIsoDR03()         );
        fillVariableWithValue( "Ele1_TrkIsoHEEP7"   , ele1.HEEP70TrackIsolation());
        fillVariableWithValue( "Ele1_EcalIsolation" , ele1.EcalIsoDR03()        );
        fillVariableWithValue( "Ele1_HcalIsolation" , ele1.HcalIsoD1DR03()      );
        fillVariableWithValue( "Ele1_CorrIsolation" , ele1.HEEPCorrIsolation()  );
        fillVariableWithValue( "Ele1_PFCHIso03"     , ele1.PFChargedHadronIso03());
        fillVariableWithValue( "Ele1_PFPhoIso03"    , ele1.PFPhotonIso03       ());
        fillVariableWithValue( "Ele1_PFNHIso03"     , ele1.PFNeutralHadronIso03());

        fillVariableWithValue("Ele1_GsfCtfScPixCharge", ele1.GsfCtfScPixCharge()  );
        fillVariableWithValue("Ele1_GsfScPixCharge"   , ele1.GsfScPixCharge()     );
        fillVariableWithValue("Ele1_GsfCtfCharge"     , ele1.GsfCtfCharge()       );

        //fillVariableWithValue( "Ele1_hltEleSignalPt", hltEle1Pt_signal          );
        ////fillVariableWithValue( "Ele1_hltEleTTbarPt" , hltEle1Pt_ttbar           );
        //fillVariableWithValue( "Ele1_hltDoubleElePt", hltEle1Pt_doubleEleSignal ); 
        //fillVariableWithValue( "Ele1_hltEleWP80Pt"  , hltEle1Pt_WP80            );

        if ( n_ele_store >= 2 ){
          Electron ele2 = c_ele_final -> GetConstituent<Electron>(1);
          double hltEle2Pt_signal          = triggerMatchPt<HLTriggerObject, Electron>(c_hltEle45_Signal_all     , ele2, ele_hltMatch_DeltaRCut);
          //double hltEle2Pt_ttbar           = triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton22_TTbar_all   , ele2, ele_hltMatch_DeltaRCut);
          double hltEle2Pt_doubleEleSignal = triggerMatchPt<HLTriggerObject, Electron>(c_hltDoubleEle_Signal_all , ele2, ele_hltMatch_DeltaRCut);
          double hltEle2Pt_WP80            = triggerMatchPt<HLTriggerObject, Electron>(c_hltEle27WP85Gsf_all        , ele2, ele_hltMatch_DeltaRCut);

          fillVariableWithValue( "Ele2_Pt"            , ele2.Pt()                 );
          fillVariableWithValue( "Ele2_Energy"        , ele2.CaloEnergy()         );
          fillVariableWithValue( "Ele2_Eta"           , ele2.Eta()                );
          fillVariableWithValue( "Ele2_Phi"           , ele2.Phi()                );
          fillVariableWithValue( "Ele2_PtHeep"        , ele2.PtHeep()             );
          fillVariableWithValue( "Ele2_SCEt"          , ele2.SCEnergy()/cosh(ele2.SCEta()) );
          fillVariableWithValue( "Ele2_SCEta"         , ele2.SCEta()              );
          fillVariableWithValue( "Ele2_SCPhi"         , ele2.SCPhi()              );
          fillVariableWithValue( "Ele2_Charge"        , ele2.Charge()             );
          fillVariableWithValue( "Ele2_Dist"          , ele2.Dist()               );
          fillVariableWithValue( "Ele2_DCotTheta"     , ele2.DCotTheta()          );
          fillVariableWithValue( "Ele2_MissingHits"   , ele2.MissingHits()        );
          fillVariableWithValue( "Ele2_TrkPt"         , ele2.TrackPt()            );
          fillVariableWithValue( "Ele2_TrkEta"        , ele2.TrackEta()           );
          fillVariableWithValue( "Ele2_SigmaEtaEta"   , ele2.SigmaEtaEta()        );
          fillVariableWithValue( "Ele2_EcalDriven"    , ele2.EcalDriven()         );
          fillVariableWithValue( "Ele2_DeltaEtaSeed"  , ele2.DeltaEtaSeed()       );
          fillVariableWithValue( "Ele2_SCEnergy"      , ele2.SCEnergy()           );
          fillVariableWithValue( "Ele2_Full5x5SigmaIEtaIEta" , ele2.Full5x5SigmaIEtaIEta() );
          fillVariableWithValue( "Ele2_Full5x5E1x5OverE5x5"  , ele2.Full5x5E1x5OverE5x5()  );
          fillVariableWithValue( "Ele2_Full5x5E2x5OverE5x5"  , ele2.Full5x5E2x5OverE5x5()  );
          fillVariableWithValue( "Ele2_RhoForHEEP"    , ele2.RhoForHEEP()         );

          fillVariableWithValue( "Ele2_DeltaPhiTrkSC" , ele2.DeltaPhi()           );
          fillVariableWithValue( "Ele2_DeltaEtaTrkSC" , ele2.DeltaEta()           );
          fillVariableWithValue( "Ele2_RawEnergy"     , ele2.RawEnergy()          );
          fillVariableWithValue( "Ele2_NBrems"        , ele2.NBrems()             );
          fillVariableWithValue( "Ele2_HoE"           , ele2.HoE()                );
          fillVariableWithValue( "Ele2_HasMatchedPhot", ele2.HasMatchedConvPhot() );
          fillVariableWithValue( "Ele2_FBrem"         , ele2.FBrem()              );
          fillVariableWithValue( "Ele2_LeadVtxDistXY" , ele2.LeadVtxDistXY()      );
          fillVariableWithValue( "Ele2_LeadVtxDistZ"  , ele2.LeadVtxDistZ ()      );
          fillVariableWithValue( "Ele2_BeamSpotDXY"   , ele2.BeamSpotDXY()        );
          fillVariableWithValue( "Ele2_BeamSpotDXYErr", ele2.BeamSpotDXYErr()     );
          fillVariableWithValue( "Ele2_ValidFrac"     , ele2.ValidFrac()          );
          fillVariableWithValue( "Ele2_Classif"       , ele2.Classif()            );
          fillVariableWithValue( "Ele2_EOverP"        , ele2.ESuperClusterOverP() );

          fillVariableWithValue( "Ele2_TrkIsolation"  , ele2.TrkIsoDR03()         );
          fillVariableWithValue( "Ele2_TrkIsoHEEP7"   , ele2.HEEP70TrackIsolation());
          fillVariableWithValue( "Ele2_EcalIsolation" , ele2.EcalIsoDR03()        );
          fillVariableWithValue( "Ele2_HcalIsolation" , ele2.HcalIsoD1DR03()      );
          fillVariableWithValue( "Ele2_CorrIsolation" , ele2.HEEPCorrIsolation()  );
          fillVariableWithValue( "Ele2_PFCHIso03"     , ele2.PFChargedHadronIso03());
          fillVariableWithValue( "Ele2_PFPhoIso03"    , ele2.PFPhotonIso03       ());
          fillVariableWithValue( "Ele2_PFNHIso03"     , ele2.PFNeutralHadronIso03());

          fillVariableWithValue("Ele2_GsfCtfScPixCharge", ele2.GsfCtfScPixCharge()  );
          fillVariableWithValue("Ele2_GsfScPixCharge"   , ele2.GsfScPixCharge()     );
          fillVariableWithValue("Ele2_GsfCtfCharge"     , ele2.GsfCtfCharge()       );

          //fillVariableWithValue( "Ele2_hltEleSignalPt", hltEle2Pt_signal          );
          ////fillVariableWithValue( "Ele2_hltEleTTbarPt" , hltEle2Pt_ttbar           );
          //fillVariableWithValue( "Ele2_hltDoubleElePt", hltEle2Pt_doubleEleSignal ); 
          //fillVariableWithValue( "Ele2_hltEleWP80Pt"  , hltEle2Pt_WP80            );

          if ( n_ele_store >= 3 ){
            Electron ele3 = c_ele_final -> GetConstituent<Electron>(2);
            double hltEle3Pt_signal          = triggerMatchPt<HLTriggerObject, Electron>(c_hltEle45_Signal_all     , ele3, ele_hltMatch_DeltaRCut);
            //double hltEle3Pt_ttbar           = triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton22_TTbar_all   , ele3, ele_hltMatch_DeltaRCut);
            double hltEle3Pt_doubleEleSignal = triggerMatchPt<HLTriggerObject, Electron>(c_hltDoubleEle_Signal_all , ele3, ele_hltMatch_DeltaRCut);
            double hltEle3Pt_WP80            = triggerMatchPt<HLTriggerObject, Electron>(c_hltEle27WP85Gsf_all        , ele3, ele_hltMatch_DeltaRCut);

            fillVariableWithValue( "Ele3_Pt"            , ele3.Pt()                 );
            fillVariableWithValue( "Ele3_Energy"        , ele3.CaloEnergy()         );
            fillVariableWithValue( "Ele3_Eta"           , ele3.Eta()                );
            fillVariableWithValue( "Ele3_Phi"           , ele3.Phi()                );
            fillVariableWithValue( "Ele3_PtHeep"        , ele3.PtHeep()             );
            fillVariableWithValue( "Ele3_SCEt"          , ele3.SCEnergy()/cosh(ele3.SCEta()) );
            fillVariableWithValue( "Ele3_SCEta"         , ele3.SCEta()              );
            fillVariableWithValue( "Ele3_SCPhi"         , ele3.SCPhi()              );
            fillVariableWithValue( "Ele3_Charge"        , ele3.Charge()             );
            fillVariableWithValue( "Ele3_Dist"          , ele3.Dist()               );
            fillVariableWithValue( "Ele3_DCotTheta"     , ele3.DCotTheta()          );
            fillVariableWithValue( "Ele3_MissingHits"   , ele3.MissingHits()        );
            fillVariableWithValue( "Ele3_TrkPt"         , ele3.TrackPt()            );
            fillVariableWithValue( "Ele3_TrkEta"        , ele3.TrackEta()           );
            fillVariableWithValue( "Ele3_SigmaEtaEta"   , ele3.SigmaEtaEta()        );
            fillVariableWithValue( "Ele3_EcalDriven"    , ele3.EcalDriven()         );
            fillVariableWithValue( "Ele3_DeltaEtaSeed"  , ele3.DeltaEtaSeed()       );
            fillVariableWithValue( "Ele3_SCEnergy"      , ele3.SCEnergy()           );
            fillVariableWithValue( "Ele3_Full5x5SigmaIEtaIEta" , ele3.Full5x5SigmaIEtaIEta() );
            fillVariableWithValue( "Ele3_Full5x5E1x5OverE5x5"  , ele3.Full5x5E1x5OverE5x5()  );
            fillVariableWithValue( "Ele3_Full5x5E2x5OverE5x5"  , ele3.Full5x5E2x5OverE5x5()  );
            fillVariableWithValue( "Ele3_RhoForHEEP"    , ele3.RhoForHEEP()         );

            fillVariableWithValue( "Ele3_DeltaPhiTrkSC" , ele3.DeltaPhi()           );
            fillVariableWithValue( "Ele3_DeltaEtaTrkSC" , ele3.DeltaEta()           );
            fillVariableWithValue( "Ele3_RawEnergy"     , ele3.RawEnergy()          );
            fillVariableWithValue( "Ele3_NBrems"        , ele3.NBrems()             );
            fillVariableWithValue( "Ele3_HoE"           , ele3.HoE()                );
            fillVariableWithValue( "Ele3_HasMatchedPhot", ele3.HasMatchedConvPhot() );
            fillVariableWithValue( "Ele3_FBrem"         , ele3.FBrem()              );
            fillVariableWithValue( "Ele3_LeadVtxDistXY" , ele3.LeadVtxDistXY()      );
            fillVariableWithValue( "Ele3_LeadVtxDistZ"  , ele3.LeadVtxDistZ ()      );
            fillVariableWithValue( "Ele3_BeamSpotDXY"   , ele3.BeamSpotDXY()        );
            fillVariableWithValue( "Ele3_BeamSpotDXYErr", ele3.BeamSpotDXYErr()     );
            fillVariableWithValue( "Ele3_ValidFrac"     , ele3.ValidFrac()          );
            fillVariableWithValue( "Ele3_Classif"       , ele3.Classif()            );
            fillVariableWithValue( "Ele3_EOverP"        , ele3.ESuperClusterOverP() );

            fillVariableWithValue( "Ele3_TrkIsolation"  , ele3.TrkIsoDR03()         );
            fillVariableWithValue( "Ele3_TrkIsoHEEP7"   , ele3.HEEP70TrackIsolation());
            fillVariableWithValue( "Ele3_EcalIsolation" , ele3.EcalIsoDR03()        );
            fillVariableWithValue( "Ele3_HcalIsolation" , ele3.HcalIsoD1DR03()      );
            fillVariableWithValue( "Ele3_CorrIsolation" , ele3.HEEPCorrIsolation()  );
            fillVariableWithValue( "Ele3_PFCHIso03"     , ele3.PFChargedHadronIso03());
            fillVariableWithValue( "Ele3_PFPhoIso03"    , ele3.PFPhotonIso03       ());
            fillVariableWithValue( "Ele3_PFNHIso03"     , ele3.PFNeutralHadronIso03());

            fillVariableWithValue("Ele3_GsfCtfScPixCharge", ele3.GsfCtfScPixCharge()  );
            fillVariableWithValue("Ele3_GsfScPixCharge"   , ele3.GsfScPixCharge()     );
            fillVariableWithValue("Ele3_GsfCtfCharge"     , ele3.GsfCtfCharge()       );

            //fillVariableWithValue( "Ele3_hltEleSignalPt", hltEle3Pt_signal          );
            ////fillVariableWithValue( "Ele3_hltEleTTbarPt" , hltEle3Pt_ttbar           );
            //fillVariableWithValue( "Ele3_hltDoubleElePt", hltEle3Pt_doubleEleSignal ); 
            //fillVariableWithValue( "Ele3_hltEleWP80Pt"  , hltEle3Pt_WP80            );

          }
        }
      }

      // Jets

      if ( n_jet_highEta_store >= 1 ) { 
	PFJet jet1 = c_pfjet_highEta_final -> GetConstituent<PFJet>(0);
	fillVariableWithValue( "HighEtaJet1_Pt" , jet1.Pt () );
	fillVariableWithValue( "HighEtaJet1_Eta", jet1.Eta() );
	fillVariableWithValue( "HighEtaJet1_Phi", jet1.Phi() );
      }
      
      if ( n_jet_store >= 1 ){

	PFJet jet1 = c_pfjet_final -> GetConstituent<PFJet>(0);
  // leading HLT jet from 200 GeV collection
	double hltJet1Pt     = triggerMatchPt<HLTriggerObject, PFJet >( c_trigger_l3jets_all     , jet1, jet_hltMatch_DeltaRCut);
	fillVariableWithValue( "Jet1_Pt"      , jet1.Pt()                         );
	fillVariableWithValue( "Jet1_JECUnc"  , jet1.JECUnc()                     );
	fillVariableWithValue( "Jet1_Energy"  , jet1.Energy()                     );
	fillVariableWithValue( "Jet1_Eta"     , jet1.Eta()                        );
	fillVariableWithValue( "Jet1_Phi"     , jet1.Phi()                        );
  fillVariableWithValue( "Jet1_btagCISV" , jet1.CombinedInclusiveSecondaryVertexBTag());
  fillVariableWithValue( "Jet1_btagCMVA" , jet1.CombinedMVABTag());
	fillVariableWithValue( "Jet1_hltJetPt"	  , hltJet1Pt     );
	
	if ( n_jet_store >= 2 ){
	  PFJet jet2 = c_pfjet_final -> GetConstituent<PFJet>(1);
	  double hltJet2Pt     = triggerMatchPt<HLTriggerObject, PFJet >( c_trigger_l3jets_all     , jet2, jet_hltMatch_DeltaRCut);
	  fillVariableWithValue( "Jet2_Pt"      , jet2.Pt()                         );
	  fillVariableWithValue( "Jet2_JECUnc"  , jet2.JECUnc()                     );
	  fillVariableWithValue( "Jet2_Energy"  , jet2.Energy()                     );
	  fillVariableWithValue( "Jet2_Eta"     , jet2.Eta()                        );
	  fillVariableWithValue( "Jet2_Phi"     , jet2.Phi()                        );
	  fillVariableWithValue( "Jet2_btagCISV" , jet2.CombinedInclusiveSecondaryVertexBTag());
    fillVariableWithValue( "Jet2_btagCMVA" , jet2.CombinedMVABTag());
	  fillVariableWithValue( "Jet2_hltJetPt"    , hltJet2Pt     );

	  if ( n_jet_store >= 3 ){
	    PFJet jet3 = c_pfjet_final -> GetConstituent<PFJet>(2);
	    double hltJet3Pt     = triggerMatchPt<HLTriggerObject, PFJet >( c_trigger_l3jets_all     , jet3, jet_hltMatch_DeltaRCut);
	    fillVariableWithValue( "Jet3_Pt"      , jet3.Pt()                         );
	    fillVariableWithValue( "Jet3_JECUnc"  , jet3.JECUnc()                     );
	    fillVariableWithValue( "Jet3_Energy"  , jet3.Energy()                     );
	    fillVariableWithValue( "Jet3_Eta"     , jet3.Eta()                        );
	    fillVariableWithValue( "Jet3_Phi"     , jet3.Phi()                        );
	    fillVariableWithValue( "Jet3_btagCISV" , jet3.CombinedInclusiveSecondaryVertexBTag());
      fillVariableWithValue( "Jet3_btagCMVA" , jet3.CombinedMVABTag());
	    fillVariableWithValue( "Jet3_hltJetPt"    , hltJet3Pt     );

	    if ( n_jet_store >= 4 ){
	      PFJet jet4 = c_pfjet_final -> GetConstituent<PFJet>(3);
	      double hltJet4Pt     = triggerMatchPt<HLTriggerObject, PFJet >( c_trigger_l3jets_all     , jet4, jet_hltMatch_DeltaRCut);
	      fillVariableWithValue( "Jet4_Pt"      , jet4.Pt()                         );
	      fillVariableWithValue( "Jet4_JECUnc"  , jet4.JECUnc()                     );
	      fillVariableWithValue( "Jet4_Energy"  , jet4.Energy()                     );
	      fillVariableWithValue( "Jet4_Eta"     , jet4.Eta()                        );
	      fillVariableWithValue( "Jet4_Phi"     , jet4.Phi()                        );
	      fillVariableWithValue( "Jet4_btagCISV" , jet4.CombinedInclusiveSecondaryVertexBTag());
        fillVariableWithValue( "Jet4_btagCMVA" , jet4.CombinedMVABTag());
	      fillVariableWithValue( "Jet4_hltJetPt"    , hltJet4Pt     );

	      if ( n_jet_store >= 5 ){
		PFJet jet5 = c_pfjet_final -> GetConstituent<PFJet>(4);
		double hltJet5Pt     = triggerMatchPt<HLTriggerObject, PFJet >( c_trigger_l3jets_all     , jet5, jet_hltMatch_DeltaRCut);
		fillVariableWithValue( "Jet5_Pt"      , jet5.Pt()                         );
		fillVariableWithValue( "Jet5_JECUnc"  , jet5.JECUnc()                     );
		fillVariableWithValue( "Jet5_Energy"  , jet5.Energy()                     );
		fillVariableWithValue( "Jet5_Eta"     , jet5.Eta()                        );
		fillVariableWithValue( "Jet5_Phi"     , jet5.Phi()                        );
		fillVariableWithValue( "Jet5_btagCISV" , jet5.CombinedInclusiveSecondaryVertexBTag());
		fillVariableWithValue( "Jet5_btagCMVA" , jet5.CombinedMVABTag());
		fillVariableWithValue( "Jet5_hltJetPt"    , hltJet5Pt     );
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
    
    t_MET.SetPtEtaPhiM( MET_pt, 0.0, MET_phi, 0.0 );

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
      if ( reducedSkimType == 0 ) // for QCD skims, use the uncorrected Pt
        t_ele1.SetPtEtaPhiM ( ele1.SCEnergy()/cosh(ele1.SCEta()), ele1.Eta(), ele1.Phi(), 0.0 );
      if ( n_ele_store >= 2 ) {
        Electron ele2 = c_ele_final -> GetConstituent<Electron>(1);
        t_ele2.SetPtEtaPhiM ( ele2.Pt(), ele2.Eta(), ele2.Phi(), 0.0 );
        if ( reducedSkimType == 0 ) // for QCD skims, use the uncorrected Pt
          t_ele2.SetPtEtaPhiM ( ele2.SCEnergy()/cosh(ele2.SCEta()), ele2.Eta(), ele2.Phi(), 0.0 );

        TLorentzVector t_ele1ele2 = t_ele1 + t_ele2;
        fillVariableWithValue ("M_e1e2" , t_ele1ele2.M ());
        fillVariableWithValue ("Pt_e1e2", t_ele1ele2.Pt());
      }
    } 

    if ( n_ele_store >= 1 ){
      
      double MT_Ele1MET = sqrt ( 2.0 * t_ele1.Pt() * t_MET.Pt() * ( 1.0 - cos ( t_MET.DeltaPhi(t_ele1))));

      TLorentzVector t_ele1MET = t_ele1 + t_MET;
      fillVariableWithValue("mDPhi_METEle1", fabs ( t_MET.DeltaPhi(t_ele1)));
      fillVariableWithValue("MT_Ele1MET"   , MT_Ele1MET); 
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
    // - H_Photon22
    // - H_Photon30
    // - H_Photon36
    // - H_Photon50
    // - H_Photon75
    // - H_Photon90
    // - H_Photon120    
    // - H_Photon175    
    //-----------------------------------------------------------------

    if ( reducedSkimType == 0 ) { 
      if ( !isData() ) {
        fillTriggerVariable ( "HLT_Photon22"  , "H_Photon22"  );
        fillTriggerVariable ( "HLT_Photon30"  , "H_Photon30"  );
        fillTriggerVariable ( "HLT_Photon36"  , "H_Photon36"  );
        fillTriggerVariable ( "HLT_Photon50"  , "H_Photon50"  );
        fillTriggerVariable ( "HLT_Photon75"  , "H_Photon75"  );
        fillTriggerVariable ( "HLT_Photon90"  , "H_Photon90"  );
        fillTriggerVariable ( "HLT_Photon120" , "H_Photon120" );
        fillTriggerVariable ( "HLT_Photon175" , "H_Photon175" );
      }
      else {
        // same as MC
        // for 2016, need to feed in extra L1 prescales
        //FIXME
        int l1prescale = 1;
        //int l1prescale = egPrescales2016.LookupPrescale("SingleEG18", run, L1PrescaleColumn);
        fillTriggerVariable ( "HLT_Photon22"  , "H_Photon22"  , l1prescale);
        //l1prescale = egPrescales2016.LookupPrescale("SingleEG26", run, L1PrescaleColumn);
        fillTriggerVariable ( "HLT_Photon30"  , "H_Photon30"  , l1prescale);
        fillTriggerVariable ( "HLT_Photon36"  , "H_Photon36"  , l1prescale);
        //
        fillTriggerVariable ( "HLT_Photon50"  , "H_Photon50"  );
        fillTriggerVariable ( "HLT_Photon75"  , "H_Photon75"  );
        fillTriggerVariable ( "HLT_Photon90"  , "H_Photon90"  );
        fillTriggerVariable ( "HLT_Photon120" , "H_Photon120" );
        fillTriggerVariable ( "HLT_Photon175" , "H_Photon175" );
      }

      bool pass_trigger = ( getVariableValue("H_Photon22") > 0 || 
          getVariableValue("H_Photon30") > 0 || 
          getVariableValue("H_Photon36") > 0 || 
          getVariableValue("H_Photon50") > 0 || 
          getVariableValue("H_Photon75") > 0 || 
          getVariableValue("H_Photon90") > 0 || 
          getVariableValue("H_Photon120"     ) > 0 || 
          getVariableValue("H_Photon175"     ) > 0 );

      fillVariableWithValue ("PassTrigger", pass_trigger ? 1 : 0 );

    }

    //-----------------------------------------------------------------
    // Signal triggers:
    // - HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*
    // 
    // MuEG triggers
    // - HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL_v*
    //-----------------------------------------------------------------
    
    else if ( reducedSkimType == 1 || reducedSkimType == 2 || reducedSkimType == 3 || reducedSkimType == 4 ) { 
      
      // search for HLT path by prefix
      // in 2015 data, trigger is different
      // this exists in special from-RAW MC and data only
      //fillTriggerVariable( "HLT_Ele27_WPLoose_Gsf_v" , "H_Ele27_WPLoose" );
      //if      ( isData )
      //{
      //  if(run >= 254227 && run <= 254914) // in Run2015C 25 ns, there is no un-eta-restricted WPLoose path
      //    fillTriggerVariable( "HLT_Ele27_eta2p1_WPLoose_Gsf_v" , "H_Ele27_WPLoose_eta2p1" );
      //  else
      //  {
      //    fillTriggerVariable( "HLT_Ele27_WPLoose_Gsf_v" , "H_Ele27_WPLoose" );
      //    fillTriggerVariable( "HLT_Ele27_eta2p1_WPLoose_Gsf_v" , "H_Ele27_WPLoose_eta2p1" );
      //    fillTriggerVariable( "HLT_Ele27_WPTight_Gsf_v" , "H_Ele27_WPTight" );
      //  }
      //}

      // XXX Don't use the above for reHLT MC

      // just search by prefix
      fillTriggerVariable( "HLT_Ele27_eta2p1_WPLoose_Gsf" , "H_Ele27_WPLoose_eta2p1" );
      fillTriggerVariable( "HLT_Ele27_eta2p1_WPTight_Gsf" , "H_Ele27_WPTight_eta2p1" );
      if(triggerExists("HLT_Ele27_WPLoose_Gsf"))
          fillTriggerVariable( "HLT_Ele27_WPLoose_Gsf" , "H_Ele27_WPLoose" );
      if(triggerExists("HLT_Ele27_WPTight_Gsf"))
        fillTriggerVariable( "HLT_Ele27_WPTight_Gsf" , "H_Ele27_WPTight" );
      fillTriggerVariable( "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50", "H_Ele45_PFJet200_PFJet50");
      fillTriggerVariable( "HLT_Mu23NoFiltersNoVtx_Photon23_CaloIdL", "H_Mu23NoFiltNoVtx_Photon23_CIdL");
      if(triggerExists("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL"))
        fillTriggerVariable( "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", "H_DoubleEle33_CIdL_GsfIdVL" ); 
      fillTriggerVariable( "HLT_Mu45_eta2p1"  , "H_Mu45_eta2p1" );
      fillTriggerVariable( "HLT_Photon175" , "H_Photon175" );
      fillTriggerVariable( "HLT_Ele105_CaloIdVT_GsfTrkIdT" , "H_Ele105_CIdVT_GsfIdT");
      fillTriggerVariable( "HLT_Ele115_CaloIdVT_GsfTrkIdT" , "H_Ele115_CIdVT_GsfIdT");
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
          //passedCut("LooseEle1_Pt"     ) ){
        passedCut("LooseEle1_SCEt"     ) ){
        //c_ele_final->examine<Electron>("final electrons");
        //c_ele_final_ptCut->examine<Electron>("final electrons with Pt cut");
          fillSkimTree();
          fillReducedSkimTree();
      }
    }

    // enujj analysis skim
    else if ( reducedSkimType == 1 ) { 
      if( passedCut("nEle_ptCut"       ) && 
          //passedCut("Ele1_Pt"          ) && 
          passedCut("Ele1_SCEt"          ) && 
          passedCut("PFMET_Type1XY_Pt") && 
          passedCut("Jet1_Pt"          ) && 
          passedCut("Jet2_Pt"          ) && 
          passedCut("sT_enujj"         ) && 
          passedCut("MT_Ele1MET"       )) {
        fillSkimTree();
        fillReducedSkimTree();
      }
    }
    
    // eejj analysis skim
    else if ( reducedSkimType == 2 ) { 
      if( passedCut("nEle_ptCut"       ) && 
          //passedCut("Ele1_Pt"          ) &&
          //passedCut("Ele2_Pt"          ) && 
          passedCut("Ele1_SCEt"          ) && 
          passedCut("Ele2_SCEt"          ) && 
          passedCut("Jet1_Pt"          ) &&
          passedCut("Jet2_Pt"          ) &&
          passedCut("sT_eejj"          ) && 
          passedCut("M_e1e2"           )) {
        fillSkimTree();
        fillReducedSkimTree();
      }
    }
    
    // Single electron skim
    else if ( reducedSkimType == 3 ) { 
      if( passedCut("nEle_ptCut"       ) && 
          //passedCut("Ele1_Pt"          ) ){
          passedCut("Ele1_SCEt"          ) ){
        fillSkimTree();
        fillReducedSkimTree();
      }
    }
    
    // Single muon skim
    else if ( reducedSkimType == 4 ) { 
      if( passedCut("nMuon_ptCut"       ) && 
          passedCut("Muon1_Pt"          ) ){
        fillSkimTree();
        fillReducedSkimTree();
      }
    }


  }  // loop over entries
  
  std::cout << "analysisClass::Loop() ends" <<std::endl;   
  
}
