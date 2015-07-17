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
#include "HLTriggerObjectCollectionHelper.h"

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
  // Get electron ID
  //--------------------------------------------------------------------------

  int ElectronID_FromFile = (int) getPreCutValue1("electron_id");
  char electron_id_name[200];
  ID ElectronID = NULL_ID;
  
  if      ( ElectronID_FromFile == 1 ) {
    ElectronID = HEEP;
    sprintf(electron_id_name, "HEEP" );
  }
  
  else if ( ElectronID_FromFile == 2 ) {
    ElectronID = EGAMMA_TIGHT;
    sprintf(electron_id_name, "EGamma Tight" );
  }
  
  else if ( ElectronID_FromFile == 3 ) {
    ElectronID = EGAMMA_MEDIUM;
    sprintf(electron_id_name, "EGamma Medium" );
  }

  //--------------------------------------------------------------------------    
  // Get the trigger
  //--------------------------------------------------------------------------

  int Trigger_FromFile = (int) getPreCutValue1("trigger");
  int n_50jets_in_trigger = 0;
  int n_200jets_in_trigger = 0;
  int n_eles_in_trigger    = 0;
  char trigger_name         [200];
  char electron_filter_name [200];
  char jet50_filter_name    [200];
  char jet200_filter_name   [200];

  //--------------------------------------------------------------------------
  // Jet ID and muon ID, which should not be changed
  //--------------------------------------------------------------------------

  ID MuonID = MUON_TIGHT_PFISO04;
  ID JetID  = PFJET_LOOSE;

  //--------------------------------------------------------------------------  
  // Define electron pt cuts
  //--------------------------------------------------------------------------

  double ele_PtCut_ANA      = getPreCutValue1("ele_PtCut");

  //--------------------------------------------------------------------------
  // Define muon pt cuts
  //--------------------------------------------------------------------------

  double muon_PtCut_ANA     = getPreCutValue1("muon_PtCut");
  double muon_EtaCut_ANA    = getPreCutValue1("muon_EtaCut");

  //--------------------------------------------------------------------------
  // Define jet pt and eta cuts
  //--------------------------------------------------------------------------

  double jet_PtCut_ANA      = getPreCutValue1("jet_PtCut");
  double jet_EtaCut         = getPreCutValue1("jet_EtaCut");
  
  //--------------------------------------------------------------------------
  // Define jet overlap cuts
  //--------------------------------------------------------------------------

  double jet_ele_DeltaRCut  = getPreCutValue1("jet_ele_DeltaRCut" );
  double jet_muon_DeltaRCut = getPreCutValue1("jet_muon_DeltaRCut");
  
  //--------------------------------------------------------------------------
  // Define trigger matching cut
  //--------------------------------------------------------------------------

  double ele_triggerMatch_DeltaRMax = getPreCutValue1("ele_triggerMatch_DeltaRMax");
  double jet_triggerMatch_DeltaRMax = getPreCutValue1("jet_triggerMatch_DeltaRMax");

  //--------------------------------------------------------------------------
  // Tell the user!
  //--------------------------------------------------------------------------

  //std::cout << "----------------------------------------------"                << std::endl;
  //std::cout << "Trigger                      = " << trigger_name               << std::endl;
  //std::cout << "Electron trigger filter name = " << electron_filter_name       << std::endl;
  //std::cout << "Jet filter name, pt > 25     = " << jet25_filter_name          << std::endl;
  //std::cout << "Jet filter name, pt > 100    = " << jet100_filter_name         << std::endl;
  //std::cout << "N(electrons) in trigger      = " << n_eles_in_trigger          << std::endl;
  //std::cout << "N(jets) in trigger, pt > 25  = " << n_25jets_in_trigger        << std::endl;
  //std::cout << "N(jets) in trigger, pt > 100 = " << n_100jets_in_trigger       << std::endl;
  //std::cout << "----------------------------------------------"                << std::endl;
  //std::cout << "Electron ID algorithm        = " << electron_id_name           << std::endl;
  //std::cout << "Pt cut on electrons          = " << ele_PtCut_ANA              << std::endl;
  //std::cout << "----------------------------------------------"                << std::endl;
  //std::cout << "Muon ID algorithm            = " << MuonID                     << std::endl;
  //std::cout << "Pt cut on muons              = " << muon_PtCut_ANA             << std::endl;
  //std::cout << "----------------------------------------------"                << std::endl;
  //std::cout << "Jet  ID algorithm            = " << JetID                      << std::endl;
  //std::cout << "|eta| cut on jets            = " << jet_EtaCut                 << std::endl;
  //std::cout << "Pt cut on jets               = " << jet_PtCut_ANA              << std::endl;
  //std::cout << "DR(jet,electron) cut on jets = " << jet_ele_DeltaRCut          << std::endl;
  //std::cout << "DR(jet,muon    ) cut on jets = " << jet_ele_DeltaRCut          << std::endl;
  //std::cout << "----------------------------------------------"                << std::endl;
  //std::cout << "Ele trig match DeltaR max    = " << ele_triggerMatch_DeltaRMax << std::endl;
  //std::cout << "Jet trig match DeltaR max    = " << jet_triggerMatch_DeltaRMax << std::endl;
  //std::cout << "----------------------------------------------"                << std::endl;
  
  //--------------------------------------------------------------------------
  // Decide which plots to save (default is to save everything)
  //--------------------------------------------------------------------------
  
  fillSkim                         (  true  ) ;
  fillAllPreviousCuts              (  true  ) ;
  fillAllOtherCuts                 ( !true  ) ;
  fillAllSameLevelAndLowerLevelCuts( !true  ) ;
  fillAllCuts                      ( !true  ) ;

  //--------------------------------------------------------------------------
  // Create TH1D's: none yet
  //--------------------------------------------------------------------------

  CreateUserTH1D("nVertex1stGenMatchedEle_pass" ,51, -0.5, 50.5);
  CreateUserTH1D("nVertex1stGenMatchedEle_total",51, -0.5, 50.5);

  CreateUserTH1D("nPileup1stGenMatchedEle_pass" ,51, -0.5, 50.5);
  CreateUserTH1D("nPileup1stGenMatchedEle_total",51, -0.5, 50.5);

  CreateUserTH1D("nVertex2ndGenMatchedEle_pass" ,51, -0.5, 50.5);
  CreateUserTH1D("nVertex2ndGenMatchedEle_total",51, -0.5, 50.5);

  CreateUserTH1D("nPileup2ndGenMatchedEle_pass" ,51, -0.5, 50.5);
  CreateUserTH1D("nPileup2ndGenMatchedEle_total",51, -0.5, 50.5);

  CreateUserTH1D("Pt1stGenMatchedEle_pass" ,200,0,1000);
  CreateUserTH1D("Eta1stGenMatchedEle_pass",200,-5.0, 5.0);
  CreateUserTH1D("Phi1stGenMatchedEle_pass",200,-3.1416, 3.1416 );

  CreateUserTH1D("Pt1stGenMatchedEle_total" ,200,0,1000);
  CreateUserTH1D("Eta1stGenMatchedEle_total",200,-5.0, 5.0);
  CreateUserTH1D("Phi1stGenMatchedEle_total",200,-3.1416, 3.1416 );
  
  CreateUserTH1D("Pt2ndGenMatchedEle_pass" ,200,0,1000);
  CreateUserTH1D("Eta2ndGenMatchedEle_pass",200,-5.0, 5.0);
  CreateUserTH1D("Phi2ndGenMatchedEle_pass",200,-3.1416, 3.1416 );

  CreateUserTH1D("Pt2ndGenMatchedEle_total" ,200,0,1000);
  CreateUserTH1D("Eta2ndGenMatchedEle_total",200,-5.0, 5.0);
  CreateUserTH1D("Phi2ndGenMatchedEle_total",200,-3.1416, 3.1416 );

  CreateUserTH1D("Pt1stGenMatchedEle_barrel_pass" ,200,0,1000);
  CreateUserTH1D("Pt1stGenMatchedEle_endcap_pass" ,200,0,1000);
  CreateUserTH1D("Pt1stGenMatchedEle_endcap1_pass",200,0,1000);
  CreateUserTH1D("Pt1stGenMatchedEle_endcap2_pass",200,0,1000);

  CreateUserTH1D("Pt1stGenMatchedEle_barrel_total" ,200,0,1000);
  CreateUserTH1D("Pt1stGenMatchedEle_endcap_total" ,200,0,1000);
  CreateUserTH1D("Pt1stGenMatchedEle_endcap1_total",200,0,1000);
  CreateUserTH1D("Pt1stGenMatchedEle_endcap2_total",200,0,1000);

  CreateUserTH1D("Pt2ndGenMatchedEle_barrel_pass" ,200,0,1000);
  CreateUserTH1D("Pt2ndGenMatchedEle_endcap_pass" ,200,0,1000);
  CreateUserTH1D("Pt2ndGenMatchedEle_endcap1_pass",200,0,1000);
  CreateUserTH1D("Pt2ndGenMatchedEle_endcap2_pass",200,0,1000);

  CreateUserTH1D("Pt2ndGenMatchedEle_barrel_total" ,200,0,1000);
  CreateUserTH1D("Pt2ndGenMatchedEle_endcap_total" ,200,0,1000);
  CreateUserTH1D("Pt2ndGenMatchedEle_endcap1_total",200,0,1000);
  CreateUserTH1D("Pt2ndGenMatchedEle_endcap2_total",200,0,1000);

  CreateUserTH1D("Pt1stJet_total", 200, 0, 200 );
  CreateUserTH1D("Pt2ndJet_total", 200, 0, 200 );
  CreateUserTH1D("Pt1stJet_pass" , 200, 0, 200 );
  CreateUserTH1D("Pt2ndJet_pass" , 200, 0, 200 );

  CreateUserTH1D("Pt1stHLTJet_total", 200, 0, 200 );
  CreateUserTH1D("Pt2ndHLTJet_total", 200, 0, 200 );
  CreateUserTH1D("Pt1stHLTJet_pass" , 200, 0, 200 );
  CreateUserTH1D("Pt2ndHLTJet_pass" , 200, 0, 200 );
  
  CreateUserTH1D("NGenParticlesKept",500,0,500);
  CreateUserTH1D("NGenLQElectrons",10,0,10);
  CreateUserTH1D("NGenLQElectronsFiducial",10,0,10);

  CreateUserTH1D("NEventsTwoGenEleECAL",38,175,2075);
  CreateUserTH1D("NEventsTwoGenEleECALPlusTrig",38,175,2075);
  //--------------------------------------------------------------------------
  // Counting variables
  //--------------------------------------------------------------------------
  std::map<int,int> numEvents_withTwoGenElectronsECALFiducial_byMass;
  std::map<int,int> numEvents_withTwoGenElectronsECALFiducial_passesTrigger_byMass;
  // mass
  int sampleMass = 0;

  //--------------------------------------------------------------------------
  // Tell the user how many entries we'll look at
  //--------------------------------------------------------------------------

  Long64_t nentries = fChain->GetEntriesFast();
  //Long64_t nentries = 10;
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
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    //-----------------------------------------------------------------
    // Get trigger names
    //-----------------------------------------------------------------

    bool trigger_enabled = true;

    if      ( Trigger_FromFile == 1 ) {
      //if ( isData ) { 
      //  if      ( run >= 191691 && run <= 194225 ) sprintf(trigger_name, "HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v4");
      //  else trigger_enabled = false;
      //} else
      //sprintf(trigger_name, "HLT_Ele30_CaloIdVT_TrkIdT_PFNoPUJet100_PFNoPUJet25_v5");
      //sprintf(electron_filter_name,"hltEle30CaloIdVTTrkIdTDphiFilter"                       );
      sprintf(jet50_filter_name   ,"hltEle45CaloIdVTGsfTrkIdTDiCentralPFJet50EleCleaned"   );
      sprintf(jet200_filter_name  ,"hltEle45CaloIdVTGsfTrkIdTCentralPFJet200EleCleaned"  );
      sprintf(trigger_name, "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v1");
      n_50jets_in_trigger  = 2;
      n_200jets_in_trigger = 1;
      n_eles_in_trigger    = 1;
    }

    //else if ( Trigger_FromFile == 2 ) {
    //  //if ( isData ){ 
    //  //  if      ( run >= 190456 && run <= 190738 ) sprintf(trigger_name, "HLT_Ele30_CaloIdVT_TrkIdT_PFJet100_PFJet25_v3");
    //  //  else trigger_enabled = false;
    //  //} else
    //  sprintf(trigger_name, "HLT_Ele30_CaloIdVT_TrkIdT_PFJet100_PFJet25_v5");
    //  //sprintf(electron_filter_name,"hltEle30CaloIdVTTrkIdTDphiFilter"                    );
    //  //sprintf(jet25_filter_name   ,"hltEle30CaloIdVTTrkIdTDiCentralPFJet25EleCleaned"    );
    //  //sprintf(jet100_filter_name  ,"hltEle30CaloIdVTTrkIdTDiCentralPFJet100EleCleaned"   );
    //  n_25jets_in_trigger  = 2;
    //  n_100jets_in_trigger = 1;
    //  n_eles_in_trigger    = 1;
    //}

    //else if ( Trigger_FromFile == 3 ) {
    //  sprintf(trigger_name        ,"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v6" );
    //  //sprintf(electron_filter_name,"hltEle30CaloIdVTTrkIdTDphiFilter"      );
    //  //sprintf(jet25_filter_name   ,""                                      );
    //  //sprintf(jet100_filter_name  ,""                                      );
    //  n_25jets_in_trigger  = 0;
    //  n_100jets_in_trigger = 0;
    //  n_eles_in_trigger    = 2;
    //}

    //else if ( Trigger_FromFile == 4 ){
    //  sprintf(trigger_name        ,"HLT_Ele27_WP80_v10"         );
    //  //sprintf(electron_filter_name,"hltEle27WP80TrackIsoFilter" );
    //  //sprintf(jet25_filter_name   ,""                           );
    //  //sprintf(jet100_filter_name  ,""                           );
    //  n_25jets_in_trigger  = 0;
    //  n_100jets_in_trigger = 0;
    //  n_eles_in_trigger    = 1;
    //}
    
    //-----------------------------------------------------------------
    // Print progress
    //-----------------------------------------------------------------
    
    bool printEvery = false;
    if( (jentry < 10 || jentry%1000 == 0) || printEvery)
      std::cout << "--------------------> analysisClass::Loop(): jentry = " << jentry << std::endl;   
    
    //-----------------------------------------------------------------
    // Did the trigger fire?
    //-----------------------------------------------------------------

    getTriggers ( HLTKey, HLTInsideDatasetTriggerNames, HLTInsideDatasetTriggerDecisions,  HLTInsideDatasetTriggerPrescales ) ; 
    bool trigger_fired = triggerFired ( trigger_name );
    int trigger_prescale = triggerPrescale ( trigger_name );
    // check
    //std::cout << "did trigger fire? " << (trigger_fired==true ? "yes" : "no") << std::endl;
        
    //-----------------------------------------------------------------
    // Define initial, inclusive collections for RECO objects
    //-----------------------------------------------------------------

    CollectionPtr gen_all   ( new Collection(*this, GenParticlePt -> size() ));
    CollectionPtr ele_all   ( new Collection(*this, ElectronPt    -> size() ));
    CollectionPtr muon_all  ( new Collection(*this, MuonPt        -> size() ));
    CollectionPtr pfjet_all ( new Collection(*this, PFJetPt       -> size() ));

    //-----------------------------------------------------------------
    // Similarly, for HLT Filter objects
    //-----------------------------------------------------------------
    
    CollectionPtr trigger_50jet_all;
    CollectionPtr trigger_200jet_all;
    CollectionPtr trigger_50jet_noEle;
    CollectionPtr trigger_200jet_noEle;
    
    ////-----------------------------------------------------------------
    //// LQ mass (avg, round to nearest 25 GeV)
    ////-----------------------------------------------------------------
    //CollectionPtr lq_gen_all = gen_all->SkimByID<GenParticle>(GEN_LQ);
    //lq_gen_all->examine<GenParticle>("LQ gen particles");
    //float massF = 0;
    //for(int i=0;i<lq_gen_all->GetSize(); ++i)
    //{
    //  massF+=lq_gen_all->GetConstituent<GenParticle>(i).Mass();
    //}
    //massF/=lq_gen_all->GetSize();
    //int thisMass = round(massF);
    //// round to nearest 50 GeV
    //if(thisMass % 50 >= 26)
    //  thisMass = thisMass + (50 - (thisMass % 50));
    //else
    //  thisMass = thisMass - (thisMass % 50);
    //std::cout << "LQ mass= " << thisMass << std::endl;
    //if(thisMass != sampleMass) // sample changed
    //  sampleMass = thisMass;

    // get mass from filename
    std::string filename = fChain->GetCurrentFile()->GetName();
    int lastSlash = filename.rfind("/");
    filename = filename.substr(lastSlash);
    std::string mass = filename.substr(filename.find("M-")+2,filename.find("_Beta")-filename.find("M-")-2);
    // test
    //std::cout << "Current file: " << filename << std::endl;
    //std::cout << "CURRENT MASS: " << mass << std::endl;
    sampleMass = stoi(mass);

    //-----------------------------------------------------------------
    // Trigger object stuff
    //-----------------------------------------------------------------
    CollectionPtr trigger_ele_eleJetJet (new Collection(*this, ElectronHLTEleJetJetMatchPt -> size() ));
    //CollectionPtr trigger_ele_all = helper.GetHLTFilterObjects ( electron_filter_name );
  //  
  //  if ( n_25jets_in_trigger  > 0 ) { 
  //    trigger_25jet_all   = helper.GetHLTFilterObjects (jet25_filter_name );
  //    trigger_25jet_noEle = trigger_25jet_all -> SkimByVetoDRMatch<HLTFilterObject,HLTFilterObject>( trigger_ele_all, 0.3 );
  //  }
  //  else {
  //    trigger_25jet_all   = CollectionPtr(new Collection (*this, 0));
  //    trigger_25jet_noEle = CollectionPtr(new Collection (*this, 0));
  //  }

  //  if ( n_100jets_in_trigger > 0 ) {
  //    trigger_100jet_all   = helper.GetHLTFilterObjects (jet100_filter_name);
  //    trigger_100jet_noEle = trigger_100jet_all -> SkimByVetoDRMatch<HLTFilterObject,HLTFilterObject>( trigger_ele_all, 0.3 );
  //  }
  //  else {
  //    trigger_100jet_all   = CollectionPtr(new Collection (*this, 0));
  //    trigger_100jet_noEle = CollectionPtr(new Collection (*this, 0));
  //  }
    CollectionPtr trigger_all (new Collection(*this, HLTriggerObjPt->size() ));
    //HLTriggerObjectCollectionHelper helper (*this,"New741"); // use 741 HLT branches
    HLTriggerObjectCollectionHelper helper (*this,"");
    //CollectionPtr trigger_l3objects_all = helper.GetL3FilterObjectsByPath(trigger_name,true); // true->verbose for testing
    CollectionPtr trigger_l3objects_all = helper.GetL3FilterObjectsByPath(trigger_name);
    //trigger_l3objects_all->examine<HLTriggerObject>("HLT L3Filter-passing objects");
    CollectionPtr trigger_l3jets_all = trigger_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_JET);
    //trigger_l3jets_all->examine<HLTriggerObject>("HLT L3Filter-passing jets");
    // which one passed the last filter?
    CollectionPtr trigger_lastObjects = helper.GetLastFilterObjectsByPath(trigger_name);
    //trigger_lastObjects->examine<HLTriggerObject>("HLT LastFilter-passing objects");
    CollectionPtr trigger_lastJets = trigger_lastObjects->SkimByID<HLTriggerObject>(TRIGGER_JET);
    //trigger_lastJets->examine<HLTriggerObject>("HLT LastFilter-passing jets");
    //// get rid of overlaps
    //CollectionPtr trigger_l3jets_noLastJet = trigger_l3jets_all -> SkimByVetoDRMatch<HLTriggerObject,HLTriggerObject>( trigger_lastJets, 0.3 );
    //trigger_l3jets_noLastJet->examine<HLTriggerObject>("HLT L3-passing jets (no last jet)");

    // electrons seem to come as TRIGGER_PHOTON most of the time
    CollectionPtr trigger_ele_all = trigger_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_PHOTON);
    //trigger_ele_all->examine<HLTriggerObject>("HLT L3-passing electrons (with ID TRIGGER_PHOTON)");
    // if not, try TRIGGER_ELECTRON
    if(trigger_ele_all->GetSize() == 0)
      trigger_ele_all = trigger_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_ELECTRON);
    // Note: could also be TRIGGER_CLUSTER?

    //-----------------------------------------------------------------
    // Skim and examine GEN electrons:
    // - Electrons that come from LQs 
    // - Electrons that are fiducial in the detector
    //-----------------------------------------------------------------
    CollectionPtr gen_lq_electrons          = gen_all          -> SkimByID      <GenParticle> ( GEN_ELE_FROM_LQ  );
    CollectionPtr gen_lq_electrons_fiducial = gen_lq_electrons -> SkimByID      <GenParticle> ( GEN_ELE_FIDUCIAL );
    
    FillUserTH1D("NGenParticlesKept"   , gen_all->GetSize());
    FillUserTH1D("NGenLQElectrons"   , gen_lq_electrons->GetSize());
    FillUserTH1D("NGenLQElectronsFiducial"   , gen_lq_electrons_fiducial->GetSize());
    ////if(gen_all->GetSize() > 50)
    ////{
    //  std::cout << "Size of gen_all: " << gen_all->GetSize() << std::endl;
    //  std::cout << "Size of gen_lq_electrons: " << gen_lq_electrons->GetSize() << std::endl;
    //  std::cout << "Size of gen_lq_electrons_fiducial: " << gen_lq_electrons_fiducial->GetSize() << std::endl;
    //  // XXX SIC TEST
    //  for (int i_gen_all = 0; i_gen_all < gen_all->GetSize(); ++i_gen_all)
    //  {
    //    GenParticle genP = gen_all -> GetConstituent<GenParticle>(i_gen_all);
    //    //if(fabs(genP.PdgId())==42 && genP.Status()==62)
    //    std::cout << genP << std::endl;
    //  }
    //std::cout << "Size of gen_all: " << gen_all->GetSize() << std::endl;
    //std::cout << "Size of gen_lq_electrons: " << gen_lq_electrons->GetSize() << std::endl;
    //std::cout << "Size of gen_lq_electrons_fiducial: " << gen_lq_electrons_fiducial->GetSize() << std::endl;
    //  // XXX SIC TEST
    ////}

    // count the events for firing the trigger
    if(gen_lq_electrons_fiducial->GetSize() >=2)
    {
      numEvents_withTwoGenElectronsECALFiducial_byMass[sampleMass]++;
      FillUserTH1D("NEventsTwoGenEleECAL",sampleMass);
      if(trigger_fired)
      {
        numEvents_withTwoGenElectronsECALFiducial_passesTrigger_byMass[sampleMass]++;
        FillUserTH1D("NEventsTwoGenEleECALPlusTrig",sampleMass);
      }
    }


    //ele_all->examine<Electron>("reco electrons");
    //-----------------------------------------------------------------
    // RECO electrons that are matched to GEN electrons
    //-----------------------------------------------------------------

    CollectionPtr ele_genMatched = ele_all -> SkimByRequireDRMatch<Electron,GenParticle>( gen_lq_electrons, 0.3 );
    //ele_genMatched->examine<Electron>("Reco electrons matched to gen electrons");

    //-----------------------------------------------------------------
    // RECO electrons that pass the ID
    //-----------------------------------------------------------------
    
    //CollectionPtr ele_ID = ele_all -> SkimByID<Electron>( ElectronID , true);// verbose for testing
    CollectionPtr ele_ID = ele_all -> SkimByID<Electron>( ElectronID );
    //ele_ID->examine<Electron>("Reco electrons passing ID");

    //-----------------------------------------------------------------
    // You cut on the PT of electrons that pass the ID
    //-----------------------------------------------------------------

    double recoEleWithID1_Pt = 0.0;
    double recoEleWithID2_Pt = 0.0;

    if ( ele_ID -> GetSize() >= 1 ) recoEleWithID1_Pt = ele_ID -> GetConstituent<Electron>(0).Pt();
    if ( ele_ID -> GetSize() >= 2 ) recoEleWithID2_Pt = ele_ID -> GetConstituent<Electron>(1).Pt();
    
  //  //-----------------------------------------------------------------
  //  // How many of those RECO electrons with ID are matched to the trigger?
  //  //-----------------------------------------------------------------

  //  CollectionPtr ele_ID_triggerMatched = ele_ID -> SkimByRequireDRMatch<Electron, HLTFilterObject> ( trigger_ele_all, ele_triggerMatch_DeltaRMax);
    //CollectionPtr ele_ID_triggerMatched = ele_ID -> SkimByRequireDRMatch<Electron, Electron> ( trigger_ele_all, ele_triggerMatch_DeltaRMax);
    CollectionPtr ele_ID_triggerMatched = ele_ID -> SkimByRequireDRMatch<Electron, HLTriggerObject> ( trigger_ele_all, ele_triggerMatch_DeltaRMax);
    //ele_ID_triggerMatched->examine<Electron>("reco+ID electrons matching L3 filter electrons");
    int nRecoEleMatchedToTrig = ele_ID_triggerMatched->GetSize();
    //std::cout << "Number of reco ele with ID matched to EleJetJet HLT path: " << nRecoEleMatchedToTrig << std::endl;

    //double hltEleMatchedToRecoEleWithID1_Pt = 0.0;
    //double hltEleMatchedToRecoEleWithID2_Pt = 0.0;
    //if ( ele_ID -> GetSize() >= 1 ) hltEleMatchedToRecoEleWithID1_Pt = ele_ID -> GetConstituent<Electron>(0).HLTEleJetJetMatchPt();
    //if ( ele_ID -> GetSize() >= 2 ) hltEleMatchedToRecoEleWithID2_Pt = ele_ID -> GetConstituent<Electron>(1).HLTEleJetJetMatchPt();
    // this looks at the embedded match, but is apparently broken
    // to be fixed FIXME
    //// also count them
    //int nRecoEleMatchedToTrig = 0;
    //for (int i = 0; i < ele_ID->GetSize(); ++i)
    //{
    //  if(Trigger_FromFile == 1 ) //TODO: better way to make sure the embedded match is done with the same path as specified?
    //  {
    //    if(ele_ID -> GetConstituent<Electron>(i).IsHLTEleJetJetMatched())
    //      ++nRecoEleMatchedToTrig;
    //  }
    //}
    //std::cout << "Number of reco ele with ID matched to EleJetJet HLT path: " << nRecoEleMatchedToTrig << std::endl;

  //  
    //-----------------------------------------------------------------
    // RECO jets that pass the ID
    //-----------------------------------------------------------------

    CollectionPtr pfjet_ID         = pfjet_all -> SkimByID<PFJet> ( JetID ) ;
    CollectionPtr pfjet_ID_etaSkim = pfjet_ID -> SkimByEtaRange<PFJet> ( -2.6, 2.6 );

    //-----------------------------------------------------------------
    // You cut on the PT of jets that pass the ID
    //-----------------------------------------------------------------
    
    double recoJetWithID1_Pt = 0.0;
    double recoJetWithID2_Pt = 0.0;

    if ( pfjet_ID -> GetSize() >= 1 ) recoJetWithID1_Pt = pfjet_ID -> GetConstituent<PFJet>(0).Pt();
    if ( pfjet_ID -> GetSize() >= 2 ) recoJetWithID2_Pt = pfjet_ID -> GetConstituent<PFJet>(1).Pt();
    //std::cout << "Size of pfjet_all: " << pfjet_all->GetSize() << std::endl;
    //std::cout << "Size of pfjet_ID: " << pfjet_ID->GetSize() << std::endl;
    //std::cout << "Size of pfjet_ID_etaSkim: " << pfjet_ID_etaSkim->GetSize() << std::endl;

    //-----------------------------------------------------------------
    // How many of those RECO jets with ID are matched to the trigger?
    //-----------------------------------------------------------------

    //if(trigger_50jet_all->GetSize() <= 0)
    //{
    //  std::cout << "trigger_50jet_all has size: " << trigger_50jet_all->GetSize() << ", so skip event!" << std::endl;
    //  continue;
    //}
    CollectionPtr pfjet_ID_trigger50Matched  = pfjet_ID -> SkimByRequireDRMatch<PFJet, HLTriggerObject> (trigger_l3jets_all, jet_triggerMatch_DeltaRMax);
    CollectionPtr pfjet_ID_trigger200Matched = pfjet_ID -> SkimByRequireDRMatch<PFJet, HLTriggerObject> (trigger_lastJets, jet_triggerMatch_DeltaRMax);
    //std::cout << "Size of pfjet_ID_trigger50Matched: " << pfjet_ID_trigger50Matched->GetSize() << std::endl;
    //std::cout << "Size of pfjet_ID_trigger200Matched: " << pfjet_ID_trigger200Matched->GetSize() << std::endl;
    //pfjet_ID->examine<PFJet>("PFJETs with ID");

    //-----------------------------------------------------------------
    // Fill variables
    //-----------------------------------------------------------------

    fillVariableWithValue("NGenEleFiducial"	     , gen_lq_electrons_fiducial  -> GetSize());
    fillVariableWithValue("PassTrigger"		     , int(trigger_fired)                     );
    fillVariableWithValue("NRecoEleWithID"           , ele_ID                     -> GetSize());	
    fillVariableWithValue("RecoEleWithID1_Pt"	     , recoEleWithID1_Pt                      );
    fillVariableWithValue("RecoEleWithID2_Pt"	     , recoEleWithID2_Pt                      );
    fillVariableWithValue("NRecoEleMatchedToTrig"    , nRecoEleMatchedToTrig                );  
    fillVariableWithValue("NRecoJetWithID"	     , pfjet_ID                   -> GetSize());	
    fillVariableWithValue("RecoJetWithID1_Pt"	     , recoJetWithID1_Pt                      );
    fillVariableWithValue("RecoJetWithID2_Pt"	     , recoJetWithID2_Pt                      );   
    fillVariableWithValue("NRecoJetMatchedToTrig50"  , pfjet_ID_trigger50Matched  -> GetSize());	
    fillVariableWithValue("NRecoJetMatchedToTrig200" , pfjet_ID_trigger200Matched -> GetSize());	
    
    //-----------------------------------------------------------------
    // Evaluate the cuts
    //-----------------------------------------------------------------

    evaluateCuts();

    //-----------------------------------------------------------------
    // Turn-on curves for jets
    //-----------------------------------------------------------------

    
    /*
    int n_trigger_25jet_all = trigger_25jet_all -> GetSize();
    int n_trigger_ele_all   = trigger_ele_all   -> GetSize();
    for (int i_trigger_25jet_all = 0; i_trigger_25jet_all < n_trigger_25jet_all; ++i_trigger_25jet_all){
      HLTFilterObject trigger_25jet = trigger_25jet_all -> GetConstituent<HLTFilterObject>(i_trigger_25jet_all);
      for ( int i_trigger_ele_all = 0; i_trigger_ele_all < n_trigger_ele_all; ++i_trigger_ele_all ){
	HLTFilterObject trigger_ele = trigger_ele_all -> GetConstituent<HLTFilterObject>(i_trigger_ele_all);
	double dr = trigger_ele.DeltaR ( & trigger_25jet );
	if ( dr < 0.3 ) { 
	  std::cout << "ERROR:" << std::endl;
	  std::cout << "\t Electron passed : " << electron_filter_name << std::endl;
	  std::cout << "\t PFJet passed    : " << jet25_filter_name << std::endl;
	  std::cout << "\t Electron        = " << trigger_ele << std::endl;
	  std::cout << "\t PFJet           = " << trigger_25jet << std::endl;
	  std::cout << "DeltaR = " << dr << std::endl;
	  exit(0);
	}
      }
    }
    */
    
    //if ( trigger_enabled ) {
    //  if ( nRecoEleMatchedToTrig >= 1 && 
    //      trigger_l3jets_all -> GetSize() >= 2 ) {
    //    FillUserTH1D("Pt1stJet_total"   , pfjet_ID_etaSkim   -> GetConstituent<PFJet>          (0).Pt());
    //    std::cout << "Fill Pt1stHLTJet_total with: " << trigger_l3jets_all  -> GetConstituent<HLTriggerObject>(0).Pt() << std::endl;
    //    FillUserTH1D("Pt1stHLTJet_total", trigger_l3jets_all  -> GetConstituent<HLTriggerObject>(0).Pt());
    //    if ( trigger_fired ) { 
    //      FillUserTH1D("Pt1stJet_pass"   , pfjet_ID_etaSkim  -> GetConstituent<PFJet>          (0).Pt());
    //      FillUserTH1D("Pt1stHLTJet_pass", trigger_l3jets_all -> GetConstituent<HLTriggerObject>(0).Pt());
    //    }
    //  }

    //  if ( nRecoEleMatchedToTrig >= 1 &&
    //      trigger_lastJets -> GetSize() >  0 ) {

    //    HLTriggerObject lead_trigger_200jet = trigger_lastJets -> GetConstituent<HLTriggerObject>(0);
    //    CollectionPtr pfjet_ID_etaSkim_noLeadHLT = pfjet_ID_etaSkim -> SkimByVetoDRMatch<PFJet,HLTriggerObject>( lead_trigger_200jet, 0.5 );

    //XXX SIC FIXME: do we need the eta skim? the collection here can be empty and cause a crash
    //    FillUserTH1D("Pt2ndJet_total" , pfjet_ID_etaSkim_noLeadHLT    -> GetConstituent<PFJet>          (0).Pt());
    //    if ( trigger_fired ) { 
    //      FillUserTH1D("Pt2ndJet_pass", pfjet_ID_etaSkim_noLeadHLT    -> GetConstituent<PFJet>          (0).Pt());
    //    }
    //  }
    //}
    
    //-----------------------------------------------------------------
    // Fill electron 1 histograms
    //-----------------------------------------------------------------
    
    int nGenMatchedEle = ele_genMatched -> GetSize();

    if ( nGenMatchedEle > 0 ) {

      Electron ele_genMatched_1 = ele_genMatched -> GetConstituent<Electron>(0);
      
      FillUserTH1D("nVertex1stGenMatchedEle_total"   , VertexChi2->size());
      FillUserTH1D("nPileup1stGenMatchedEle_total"   , PileUpInteractions->at(0));
      FillUserTH1D("Pt1stGenMatchedEle_total"        , ele_genMatched_1.Pt()  );
      FillUserTH1D("Eta1stGenMatchedEle_total"       , ele_genMatched_1.Eta() );
      FillUserTH1D("Phi1stGenMatchedEle_total"       , ele_genMatched_1.Phi() );

      if      ( ele_genMatched_1.IsEB() )
	FillUserTH1D("Pt1stGenMatchedEle_barrel_total" , ele_genMatched_1.Pt() );
      
      else if ( ele_genMatched_1.IsEE() ){
	FillUserTH1D("Pt1stGenMatchedEle_endcap_total" , ele_genMatched_1.Pt() );
	if ( fabs(ele_genMatched_1.Eta()) <= 2.00 )  
	  FillUserTH1D("Pt1stGenMatchedEle_endcap1_total", ele_genMatched_1.Pt() );
	else                                            
	  FillUserTH1D("Pt1stGenMatchedEle_endcap2_total", ele_genMatched_1.Pt() );
      }
      
      if ( ele_genMatched_1.PassUserID ( ElectronID ) ){
	
	FillUserTH1D("nVertex1stGenMatchedEle_pass"   , VertexChi2->size());
	FillUserTH1D("nPileup1stGenMatchedEle_pass"   , PileUpInteractions->at(0));
	FillUserTH1D("Pt1stGenMatchedEle_pass"        , ele_genMatched_1.Pt()  );
	FillUserTH1D("Eta1stGenMatchedEle_pass"       , ele_genMatched_1.Eta() );
	FillUserTH1D("Phi1stGenMatchedEle_pass"       , ele_genMatched_1.Phi() );
	
	if      ( ele_genMatched_1.IsEB() )
	  FillUserTH1D("Pt1stGenMatchedEle_barrel_pass" , ele_genMatched_1.Pt() );
	else if ( ele_genMatched_1.IsEE() ){
	  FillUserTH1D("Pt1stGenMatchedEle_endcap_pass" , ele_genMatched_1.Pt() );
	  if ( fabs(ele_genMatched_1.Eta()) <= 2.00 ) 
	    FillUserTH1D("Pt1stGenMatchedEle_endcap1_pass", ele_genMatched_1.Pt() ); 
	  else                             
	    FillUserTH1D("Pt1stGenMatchedEle_endcap2_pass", ele_genMatched_1.Pt() );               
	}
      }
    }

    //-----------------------------------------------------------------
    // Fill electron 2 histograms
    //-----------------------------------------------------------------

    if ( nGenMatchedEle > 1 ) {

      Electron ele_genMatched_2 = ele_genMatched -> GetConstituent<Electron>(1);
      
      FillUserTH1D("nVertex2ndGenMatchedEle_total"   , VertexChi2->size());
      FillUserTH1D("nPileup2ndGenMatchedEle_total"   , PileUpInteractions->at(0));
      FillUserTH1D("Pt2ndGenMatchedEle_total"        , ele_genMatched_2.Pt()  );
      FillUserTH1D("Eta2ndGenMatchedEle_total"       , ele_genMatched_2.Eta() );
      FillUserTH1D("Phi2ndGenMatchedEle_total"       , ele_genMatched_2.Phi() );

      if      ( ele_genMatched_2.IsEB() )
	FillUserTH1D("Pt2ndGenMatchedEle_barrel_total" , ele_genMatched_2.Pt() );
      
      else if ( ele_genMatched_2.IsEE() ){
	FillUserTH1D("Pt2ndGenMatchedEle_endcap_total" , ele_genMatched_2.Pt() );
	if ( fabs(ele_genMatched_2.Eta()) <= 2.00 )  
	  FillUserTH1D("Pt2ndGenMatchedEle_endcap1_total", ele_genMatched_2.Pt() );
	else                                            
	  FillUserTH1D("Pt2ndGenMatchedEle_endcap2_total", ele_genMatched_2.Pt() );
      }
      
      if ( ele_genMatched_2.PassUserID ( ElectronID ) ){
	
	FillUserTH1D("nVertex2ndGenMatchedEle_pass"   , VertexChi2->size());
	FillUserTH1D("nPileup2ndGenMatchedEle_pass"   , PileUpInteractions->at(0));
	FillUserTH1D("Pt2ndGenMatchedEle_pass"        , ele_genMatched_2.Pt()  );
	FillUserTH1D("Eta2ndGenMatchedEle_pass"       , ele_genMatched_2.Eta() );
	FillUserTH1D("Phi2ndGenMatchedEle_pass"       , ele_genMatched_2.Phi() );
	
	if      ( ele_genMatched_2.IsEB() )
	  FillUserTH1D("Pt2ndGenMatchedEle_barrel_pass" , ele_genMatched_2.Pt() );
	else if ( ele_genMatched_2.IsEE() ){
	  FillUserTH1D("Pt2ndGenMatchedEle_endcap_pass" , ele_genMatched_2.Pt() );
	  if ( fabs(ele_genMatched_2.Eta()) <= 2.00 ) 
	    FillUserTH1D("Pt2ndGenMatchedEle_endcap1_pass", ele_genMatched_2.Pt() ); 
	  else                             
	    FillUserTH1D("Pt2ndGenMatchedEle_endcap2_pass", ele_genMatched_2.Pt() );               
	}
      }
    }
    
    //-----------------------------------------------------------------
    // (If verbose == true) Tell me about the electrons and jets 
    //-----------------------------------------------------------------

    if ( verbose ) { 
      std::cout << "------------------------------------------------------" << std::endl;
      gen_lq_electrons_fiducial  -> examine<GenParticle    >("GEN electrons from LQs");
      //trigger_ele_all            -> examine<HLTFilterObject>("HLT electrons");
      ele_ID                     -> examine<Electron       >("RECO electrons");
      ele_ID_triggerMatched      -> examine<Electron       >("RECO electrons, trigger matched");
      trigger_l3jets_all          -> examine<HLTriggerObject>("HLT PFJets");
      pfjet_ID                   -> examine<PFJet          >("RECO PFJets");
      pfjet_ID_trigger50Matched -> examine<PFJet          >("RECO PFJets, trigger 50 matched");
      pfjet_ID_trigger200Matched -> examine<PFJet          >("RECO PFJets, trigger 200 matched");
    }

  } // End loop over events
  
  std::cout << "analysisClass::Loop() ends" <<std::endl;   

  for(std::map<int,int>::iterator itr = numEvents_withTwoGenElectronsECALFiducial_byMass.begin();
      itr != numEvents_withTwoGenElectronsECALFiducial_byMass.end(); ++itr)
  {
    //int idx = std::distance(numEvents_withTwoGenElectronsECALFiducial_byMass.begin(),itr);
    int numPassTrig = numEvents_withTwoGenElectronsECALFiducial_passesTrigger_byMass[itr->first];
    float triggerEff = numPassTrig/(float)itr->second;
    std::cout << "Mass = " << itr->first << ", Trigger eff = " << numPassTrig << " / " << itr->second <<
      " = " << triggerEff << std::endl;
  }
}

