#define analysisClass_cxx
#include "analysisClass.h"
#include <typeinfo>
#include <variant>
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TLeaf.h>

#include "Collection.h"
#include "GenParticle.h"
#include "Electron.h"
#include "LooseElectron.h"
#include "Muon.h"
#include "PFJet.h"
#include "GenJet.h"
#include "HLTriggerObject.h"
#include "HLTriggerObjectCollectionHelper.h"
#include "HistoReader.h"
#include "ElectronScaleFactors.C"
#include "JMEUncertainties.C"

#include "CorrectionHandler.h"

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

void analysisClass::Loop()
{
  std::cout << "analysisClass::Loop(): begins " << std::endl;
  
  //--------------------------------------------------------------------------
  // Verbose? Or not?
  //--------------------------------------------------------------------------
  
  bool verbose = !true;
  
  //--------------------------------------------------------------------------
  // Decide which plots to save (default is to save everything)
  //--------------------------------------------------------------------------
  
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
  // Trigger scale factors
  //--------------------------------------------------------------------------
  std::string trigSFFileName = getPreCutString1("TriggerSFFileName");
  std::string trigSFBarrelHistName = hasPreCut("TriggerSFBarrelHistName") ? getPreCutString1("TriggerSFBarrelHistName") : "EGamma_SF2D";;
  std::string trigSFEndcapHistName = hasPreCut("TriggerSFEndcapHistName") ? getPreCutString1("TriggerSFEndcapHistName") : "EGamma_SF2D";;
  HistoReader triggerScaleFactorReader(trigSFFileName,trigSFBarrelHistName, trigSFEndcapHistName, false, false);

  //--------------------------------------------------------------------------
  // EGM scale factors
  //--------------------------------------------------------------------------
  std::string egmIDJSONFileName = getPreCutString1("EGMIDJSONFileName");
  auto cset_ulEGMID = CorrectionHandler::GetCorrectionFromFile(egmIDJSONFileName, "UL-Electron-ID-SF");

  //--------------------------------------------------------------------------
  // EGM scale factors
  //--------------------------------------------------------------------------
  std::string egmScaleJSONFileName = getPreCutString1("EGMScaleJSONFileName");
  auto cset_ulEGMScale = CorrectionHandler::GetCorrectionFromFile(egmScaleJSONFileName, "UL-EGM_ScaleUnc");

  //--------------------------------------------------------------------------
  // correctionlib B-tagging
  //--------------------------------------------------------------------------
  std::string btagJSONFileName = getPreCutString1("BTVJSONFileName");
  auto cset_btagDeepJetIncl = CorrectionHandler::GetCorrectionFromFile(btagJSONFileName, "deepJet_incl");
  auto cset_btagDeepJetComb = CorrectionHandler::GetCorrectionFromFile(btagJSONFileName, "deepJet_comb");
  auto cset_btagDeepJetMuJets = CorrectionHandler::GetCorrectionFromFile(btagJSONFileName, "deepJet_mujets");

  //--------------------------------------------------------------------------
  // Analysis year
  //--------------------------------------------------------------------------
  std::string analysisYear = getPreCutString1("AnalysisYear");
  int analysisYearInt;
  if(analysisYear.find("2016") != std::string::npos)
    analysisYearInt = std::stoi(analysisYear.substr(analysisYear.find("2016"), 4));
  else if(analysisYear.size() != 0)
    analysisYearInt = std::stoi(analysisYear);
  else {
    analysisYearInt = getPreCutValue1("AnalysisYear");
    analysisYear = to_string(analysisYearInt);
  }

  //--------------------------------------------------------------------------
  // correctionlib PU
  //--------------------------------------------------------------------------
  std::string puJSONFileName = getPreCutString1("PUJSONFileName");
  // need to use last 2 digits of analysisYear
  std::string puEntry = "Collisions"+to_string(analysisYearInt).replace(0, 2, "")+"_UltraLegacy_goldenJSON";
  auto cset_pu = CorrectionHandler::GetCorrectionFromFile(puJSONFileName, puEntry);
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
  // IDs
  //--------------------------------------------------------------------------
  std::string electronIDType     = getPreCutString1("electronIDType");
  std::string muonIDType         = getPreCutString1("muonIDType");
  std::string jetIDType          = getPreCutString1("jetIDType");
  if(electronIDType != "HEEP" && electronIDType != "EGMLoose") {
    STDOUT("electronIDType=" << electronIDType << " is unknown! Please implement it in the analysisClass code. Exiting.");
    exit(-5);
  }
  if(muonIDType != "HighPtTrkRelIso03" && muonIDType != "Tight" && muonIDType != "Loose") {
    STDOUT("muonIDType=" << muonIDType << " is unknown! Please implement it in the analysisClass code. Exiting.");
    exit(-5);
  }
  if(jetIDType != "PFJetTight") {
    STDOUT("jetIDType=" << jetIDType << " is unknown! Please implement it in the analysisClass code. Exiting.");
    exit(-5);
  }

  //--------------------------------------------------------------------------
  // Tell the user how many entries we'll look at
  //--------------------------------------------------------------------------

  Long64_t nentries = GetTreeEntries();
  std::cout << "analysisClass::analysisClass(): nentries = " << nentries << std::endl;
  
  //--------------------------------------------------------------------------
  // Create HLT collections in advance (won't need all of them)
  //--------------------------------------------------------------------------
  
  // trigger objects matching ID
  CollectionPtr c_hlTriggerObjects_id;

  //--------------------------------------------------------------------------
  // For smearing systematics samples, you'll need a random number generator
  //--------------------------------------------------------------------------

  unsigned int seed = 987654321;
  TRandom3 * rootEngine = new TRandom3 ( seed ) ;

  //--------------------------------------------------------------------------
  // correctionlib JME
  //--------------------------------------------------------------------------
  //std::string jmeJSONFileName = getPreCutString1("JMEJSONFileName");
  //CorrectionHandler cset_jmeAll(jmeJSONFileName);
  //--------------------------------------------------------------------------
  // For JES/JER/MET uncertainties
  //--------------------------------------------------------------------------
  std::string jecTextFilePath = getPreCutString1("JECTextFilePath");
  std::string jecTag = getPreCutString1("JECTag");
  std::string jerTextFilePath = getPreCutString1("JERTextFilePath");
  std::string jerTag = getPreCutString1("JERTag");
  
  static const std::vector<std::string> jesUncertainties {"AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias", "Fragmentation",
    "SinglePionECAL", "SinglePionHCAL", "FlavorQCD", "TimePtEta", "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF", "RelativePtBB",
    "RelativePtEC1", "RelativePtEC2", "RelativePtHF", "RelativeBal", "RelativeSample", "RelativeFSR", "RelativeStatFSR", "RelativeStatEC",
    "RelativeStatHF", "PileUpDataMC", "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF",
    "Total", "FlavorZJet", "FlavorPureGluon", "FlavorPureQuark", "FlavorPureBottom"};

  bool splitJER = false;
  JMEUncertainties jetUncertainties(JetVariationsCalculator(), this, jesUncertainties, jecTextFilePath, jerTextFilePath, jecTag, jerTag, splitJER);
  JMEUncertainties type1METUncertainties(Type1METVariationsCalculator(), this, jesUncertainties, jecTextFilePath, jerTextFilePath, jecTag, jerTag, splitJER);
  //-----------------------------------------------------------------
  // If this is MC, smear jets if requested
  // Don't do it for data
  //-----------------------------------------------------------------
  //FIXME TODO: implement the do_jer variable here?
  if (!isData()) {
    jetUncertainties.setSmearing();
    type1METUncertainties.setSmearing();
  }
  if(analysisYearInt == 2018) {
    jetUncertainties.setAddHEM2018Issue(true);
    type1METUncertainties.setAddHEM2018Issue(true);
  }

  /*//------------------------------------------------------------------
   *
   *
   *      
   *      Start analysis loop!
   *
   *
   *
   *///-----------------------------------------------------------------
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    readerTools_->LoadEntry(jentry);
    // run ls event
    //unsigned long long int event = readerTools_->ReadValueBranch<ULong64_t>("event");
    //unsigned int ls = readerTools_->ReadValueBranch<UInt_t>("luminosityBlock");
    //unsigned int run = readerTools_->ReadValueBranch<UInt_t>("run");
    //std::cout << run << " " << ls << " " << event << std::endl;
    ////std::string current_file_name ( readerTools_->GetTree()->GetCurrentFile()->GetName());
    ////cout << "Found the event! in file:" << current_file_name << endl;
    //if(jentry > 10000) continue;

    //-----------------------------------------------------------------
    // Print progress
    //-----------------------------------------------------------------
    if(jentry < 10 || jentry%5000 == 0) std::cout << "analysisClass:Loop(): jentry = " << jentry << "/" << nentries << std::endl;
    
    //-----------------------------------------------------------------
    // Get access to HLT decisions
    //-----------------------------------------------------------------

    getTriggers(jentry); 
    //printTriggers();
    //printFiredTriggers();

    //-----------------------------------------------------------------
    // Get access to HLT filter objects
    //-----------------------------------------------------------------
    HLTriggerObjectCollectionHelper helper(*this);

    std::vector<int> typeIds {11, 22};
    c_hlTriggerObjects_id = helper.GetTriggerObjectsByType(typeIds);
    //c_hlTriggerObjects_id->examine<HLTriggerObject>("c_hlTriggerObjectsid");

    //-----------------------------------------------------------------
    // Define initial, inclusive collections for physics objects
    //-----------------------------------------------------------------

    CollectionPtr c_gen_all(new Collection(readerTools_));
    if(!isData()) {
      c_gen_all.reset(new Collection(readerTools_, readerTools_->ReadValueBranch<UInt_t>("nGenPart")));
    }
    CollectionPtr c_ele_all   ( new Collection(readerTools_, readerTools_->ReadValueBranch<UInt_t>("nElectron")));
    //c_ele_all->examine<Electron>("c_ele_all = All reco electrons");
    //Electron ele1 = c_ele_all -> GetConstituent<Electron>(0);
    //for(unsigned int i=0; i<10; ++i) { 
    //  std::cout << "cut = " << i << " idLevel = " << ele1.GetNbitFromBitMap(i, 3) << std::endl;
    //}

    CollectionPtr c_muon_all  ( new Collection(readerTools_, readerTools_->ReadValueBranch<UInt_t>("nMuon")));
    CollectionPtr c_genJet_all(new Collection(readerTools_));
    if(!isData()) {
      c_genJet_all.reset(new Collection(readerTools_, readerTools_->ReadValueBranch<UInt_t>("nGenJet")));
    }
    //c_genJet_all->examine<GenJet>("c_genJet_all");
    CollectionPtr c_pfjet_all ( new Collection(readerTools_, readerTools_->ReadValueBranch<UInt_t>("nJet")));
    //  New: Cut all jet collections at 17 GeV as per 2016 custom skim
    //c_pfjet_all = c_pfjet_all -> SkimByMinPt   <PFJet>( 17.0 );
    //c_pfjet_all->examine<PFJet>("c_pfjet_all with Pt>17");

    // recall that only jets with pt > 10 GeV affect the PFMET
    // Also, only scale/smear the jets in our eta range (jets in the calorimeter crack are suspect)
    c_pfjet_all = c_pfjet_all -> SkimByEtaRange<PFJet>( -jet_EtaCut, jet_EtaCut );

    //-----------------------------------------------------------------
    // All skims need GEN particles/jets
    //-----------------------------------------------------------------

    //c_gen_all->examine<GenParticle>("c_gen_all");
    CollectionPtr c_genEle_final = c_gen_all    -> SkimByID<GenParticle>(GEN_ELE_HARDPROCESS_FINALSTATE);
    //c_genEle_final->examine<GenParticle>("c_genEle_final = c_gen_all after SkimByID<GenParticle>(GEN_ELE_HARDPROCESS_FINALSTATE)");

    CollectionPtr c_genEle_fromLQ_finalState = c_gen_all -> SkimByID<GenParticle>(GEN_ELE_FROM_LQ);
    c_genEle_fromLQ_finalState = c_genEle_fromLQ_finalState -> SkimByID<GenParticle>(GEN_ELE_HARDPROCESS_FINALSTATE);
    //c_genEle_fromLQ_finalState->examine<GenParticle>("c_genEle_fromLQ_finalState = (GEN_ELE_FROM_LQ && GEN_ELE_HARDPROCESS_FINALSTATE)");

    CollectionPtr c_genMu_final = c_gen_all    -> SkimByID<GenParticle>(GEN_MU_HARD_SCATTER);

    CollectionPtr c_genNu_final = c_gen_all    -> SkimByID<GenParticle>(GEN_NU_HARD_SCATTER);
    //c_genNu_final->examine<GenParticle>("c_genNu_final = c_gen_all after SkimByID<GenParticle>(GEN_NU_HARD_SCATTER)");
    //c_genJet_all->examine<GenJet>("c_genJet_all");
    CollectionPtr c_genJet_final = c_genJet_all;

    //c_genZgamma_final->examine<GenParticle>("c_genZgamma_final = c_gen_all after SkimByID<GenParticle>(GEN_ZGAMMA_HARD_SCATTER)");
    CollectionPtr c_genNuFromW_final   = c_gen_all -> SkimByID<GenParticle>(GEN_NU_FROM_W);
    CollectionPtr c_genTop             = c_gen_all -> SkimByID<GenParticle>(GEN_TOP);
    CollectionPtr c_genTop_final       = c_genTop  -> SkimByID<GenParticle>(GEN_STATUS62);

    CollectionPtr c_genLQ              = c_gen_all ->SkimByID<GenParticle>(GEN_LQ);
    CollectionPtr c_genLQ_final        = c_genLQ  -> SkimByID<GenParticle>(GEN_STATUS62);
    //c_genLQ_final->examine<GenParticle>("c_genLQ = c_gen_all after SkimByID GEN_LQ and GEN_STATUS62");

    CollectionPtr c_genQuark_hardScatter = c_gen_all -> SkimByID<GenParticle>(GEN_QUARK_HARD_SCATTER);
    //c_genQuark_hardScatter->examine<GenParticle>("c_genQuark_hardScatter");
    const CollectionPtr c_genQuark_hardScatter_const(c_genQuark_hardScatter);
    CollectionPtr c_genJetMatchedLQ = c_genJet_all -> SkimByRequireDRMatch<GenJet, GenParticle>(c_genQuark_hardScatter_const, 0.1);
   // c_genJetMatchedLQ->examine<GenJet>("c_genJetMatchedLQ");

    // old systematics handling for electron energy resolution and scale
    //std::vector<Electron> smearedEles;
    //std::vector<Electron> scaledUpEles;
    //std::vector<Electron> scaledDownEles;
    //if(!isData()) {
    //  smearedEles = c_ele_all->MatchAndSmearEnergy <Electron, GenParticle> ( c_genEle_final, 0.4/2.0, rootEngine, v_delta_met );
    //  scaledUpEles = c_ele_all->ScaleEnergy <Electron> (1 , v_delta_met );
    //  scaledDownEles = c_ele_all->ScaleEnergy <Electron> (-1, v_delta_met );
    //}

    //-----------------------------------------------------------------
    // QCD skims    (reducedSkimType = 0     ) have loose electrons
    // Signal skims (reducedSkimType = 1 - 4 ) have HEEP  electrons
    //-----------------------------------------------------------------

    CollectionPtr c_ele_HEEP;
    CollectionPtr c_ele_egammaLoose;
    CollectionPtr c_ele_final;
    CollectionPtr c_ele_final_ptCut;
    CollectionPtr c_ele_vLoose_ptCut;
    CollectionPtr c_ele_loose;
    CollectionPtr c_ele_vLoose;
    ID heepIdType = HEEP70;
    if(analysisYearInt == 2018)
      heepIdType = HEEP70_2018;

    if ( reducedSkimType == 0 ){ 
      if(electronIDType == "HEEP") {
        c_ele_loose = c_ele_all   -> SkimByID<LooseElectron>( FAKE_RATE_HEEP_LOOSE);
        c_ele_vLoose = c_ele_all   -> SkimByID<LooseElectron>(FAKE_RATE_VERY_LOOSE);
      }
      else if(electronIDType == "EGMLoose") {
        c_ele_loose = c_ele_all   -> SkimByID<LooseElectron>( FAKE_RATE_EGMLOOSE );
        c_ele_vLoose = c_ele_all   -> SkimByID<LooseElectron>(FAKE_RATE_VERY_LOOSE_EGMLOOSE);
      }
      c_ele_final               = c_ele_loose;
      c_ele_final_ptCut         = c_ele_final -> SkimByMinPt<LooseElectron>( ele_PtCut  );
      c_ele_vLoose_ptCut         = c_ele_vLoose -> SkimByMinPt<LooseElectron>( 10.0 );
      //c_ele_all->examine<Electron>("c_ele_all");
      //c_ele_loose->examine<Electron>("c_ele_loose");
      //c_ele_final_ptCut->examine<Electron>("c_ele_final_ptCut");
    }

    else if ( reducedSkimType == 1 || reducedSkimType == 2 || reducedSkimType == 3 || reducedSkimType == 4 ){
      if(electronIDType == "HEEP") {
        c_ele_HEEP  = c_ele_all -> SkimByID <Electron> ( heepIdType );
        //CollectionPtr c_ele_HEEP  = c_ele_all -> SkimByID <Electron> ( HEEP70_MANUAL , true );
        c_ele_final               = c_ele_HEEP;
        c_ele_vLoose = c_ele_all   -> SkimByID<LooseElectron>(FAKE_RATE_VERY_LOOSE);
      }
      else if(electronIDType == "EGMLoose") {
        c_ele_egammaLoose = c_ele_all -> SkimByID<Electron>(EGAMMA_LOOSE_HEEPETACUT);
        c_ele_final               = c_ele_egammaLoose;
        c_ele_vLoose = c_ele_all   -> SkimByID<LooseElectron>(FAKE_RATE_VERY_LOOSE_EGMLOOSE);
      }
      c_ele_final_ptCut         = c_ele_final -> SkimByMinPt<Electron>( ele_PtCut  );
      c_ele_vLoose_ptCut         = c_ele_vLoose -> SkimByMinPt<LooseElectron>( 10.0 );
    }
    // look at final electrons
    //Electron ele1_tmp;
    //if(c_ele_final->GetSize() > 0)
    //{
    //  ele1_tmp = c_ele_final -> GetConstituent<Electron>(0);
    //  if(ele1_tmp.Pt() > 2300)
    //  {
    //    std::cout << static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(luminosityBlock) << " " << static_cast<ULong64_t>(event) << std::endl;
    //    c_ele_all->examine<Electron>("c_ele_all");
    //    CollectionPtr c_ele_HEEP  = c_ele_all -> SkimByID <Electron> ( HEEP70_MANUAL , true );
    //    c_ele_HEEP->examine<Electron>("c_ele_HEEP_manual");
    //    c_ele_final->examine<Electron>("c_ele_final");
    //  }
    //}
    ////c_ele_final_ptCut->examine<Electron>("c_ele_final_ptCut");
    //c_ele_final->examine<Electron>("c_ele_final");
    FillUserTH1D("nEleNTuple",c_ele_all->GetSize());
    FillUserTH1D("nEleNrsk",c_ele_final->GetSize());
    FillUserTH2D("nEleNTupleVsNeleRsk",c_ele_final->GetSize(),c_ele_all->GetSize());

    //-----------------------------------------------------------------
    // All skims need muons
    //-----------------------------------------------------------------
    CollectionPtr c_muon_final;

    CollectionPtr c_muon_eta               = c_muon_all       -> SkimByEtaRange<Muon> ( -muon_EtaCut, muon_EtaCut );
    CollectionPtr c_muon_eta_IDLoose       = c_muon_eta       -> SkimByID      <Muon> ( MUON_LOOSE);
    CollectionPtr c_muon_eta_IDHighPt      = c_muon_eta       -> SkimByID      <Muon> ( MUON_HIGH_PT_TRKRELISO03);
    if(muonIDType == "HighPtTrkRelIso03") {
      c_muon_final             = c_muon_eta_IDHighPt;
    }
    else if(muonIDType == "Tight") {
      CollectionPtr c_muon_eta_IDTight       = c_muon_eta       -> SkimByID      <Muon> ( MUON_TIGHT_PFISO04TIGHT );
      c_muon_final             = c_muon_eta_IDTight;
    }
    else if(muonIDType == "Loose") {
      c_muon_final             = c_muon_eta_IDLoose;
    }
    //CollectionPtr c_muon_final             = c_muon_eta_IDLoose;
    CollectionPtr c_muon_final_ptCut       = c_muon_final     -> SkimByMinPt   <Muon> ( muon_PtCut );
    //c_muon_all->examine<Muon>("c_muon_all");
    //c_muon_final->examine<Muon>("c_muon_final");
    //c_muon_eta_IDLoose->examine<Muon>("c_muon_eta_IDLoose");

    //-----------------------------------------------------------------
    // All skims need PFJets
    //-----------------------------------------------------------------
    CollectionPtr c_pfjet_central_ID;

    //c_pfjet_all->examine<PFJet>("c_pfjet_all");
    CollectionPtr c_pfjet_central                     = c_pfjet_all                          -> SkimByEtaRange   <PFJet>          ( -jet_EtaCut, jet_EtaCut   );
    //c_pfjet_central->examine<PFJet>("c_pfjet_central");
    if(jetIDType == "PFJetTight")
      c_pfjet_central_ID                  = c_pfjet_central                      -> SkimByID         <PFJet>          ( PFJET_TIGHT );    
    //c_pfjet_central_ID->examine<PFJet>("c_pfjet_central_ID");
    CollectionPtr c_pfjet_central_ID_noMuonOverlap    = c_pfjet_central_ID                   -> SkimByVetoDRMatch<PFJet, Muon>    ( c_muon_final_ptCut   , jet_muon_DeltaRCut  );
    //c_pfjet_central_ID_noMuonOverlap->examine<PFJet>("c_pfjet_central_ID_noMuonOverlap [after muon cleaning]");
    CollectionPtr c_pfjet_central_ID_noLeptonOverlap  = c_pfjet_central_ID_noMuonOverlap     -> SkimByVetoDRMatch<PFJet, Electron>( c_ele_final_ptCut    , jet_ele_DeltaRCut );
    //c_pfjet_central_ID_noLeptonOverlap->examine<PFJet>("c_pfjet_central_ID_noLeptonOverlap [after electron cleaning]");
    CollectionPtr c_pfjet_final                       = c_pfjet_central_ID_noLeptonOverlap;
    //c_pfjet_final->examine<PFJet>("c_pfjet_final");
    CollectionPtr c_pfjet_final_ptCut                 = c_pfjet_final                        -> SkimByMinPt      <PFJet>          ( jet_PtCut );

    //c_pfjet_final_ptCut->examine<PFJet>("c_pfjet_final_ptCut");
    //if(c_pfjet_final->GetSize() > 4) {
    //  // run ls event
    //  //double run = readerTools_->ReadValueBranch<UInt_t>("run");
    //  double event = readerTools_->ReadValueBranch<ULong64_t>("event");
    //  //double ls = readerTools_->ReadValueBranch<UInt_t>("luminosityBlock");
    //  std::cout << static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(ls) << " " << static_cast<unsigned int>(event) << std::endl;
    //  c_pfjet_all->examine<PFJet>("c_pfjet_all");
    //  c_pfjet_final->examine<PFJet>("c_pfjet_final");
    //  c_ele_all->examine<Electron>("c_ele_all");
    //  c_ele_final->examine<Electron>("c_ele_final");
    //  c_muon_all->examine<Muon>("c_muon_all");
    //  c_muon_final->examine<Muon>("c_muon_final");
    //}
    //PFJet jet2 = c_pfjet_final -> GetConstituent<PFJet>(1);
    //if(jet2.Pt() < 50) {
    //  std::cout << "----> Interesting event!" << std::endl;
    //  // run ls event
    //  double run = readerTools_->ReadValueBranch<UInt_t>("run");
    //  double event = readerTools_->ReadValueBranch<ULong64_t>("event");
    //  double ls = readerTools_->ReadValueBranch<UInt_t>("luminosityBlock");
    //  std::cout << static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(ls) << " " << static_cast<unsigned int>(event) << std::endl;
    //}
    const CollectionPtr c_pfjet_all_const(c_pfjet_all);
    CollectionPtr c_pfJetMatchedLQ = c_genJetMatchedLQ -> SkimByRequireDRMatch<GenJet, PFJet>(c_pfjet_all_const, 0.1);
    //c_pfJetMatchedLQ->examine<PFJet>("c_pfJetMatchedLQ");

    //-----------------------------------------------------------------
    // Energy scaling and resolution smearing here
    //-----------------------------------------------------------------

    TLorentzVector v_delta_met;
    // Set the PFMET difference to zero
    v_delta_met.SetPtEtaPhiM(0.,0.,0.,0.);

    // Do the energy scale / energy resolution operations
    // dR for matching = Rcone/2
    //   see: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
    //std::string jerPtRes = jerTag+"_PtResolution_AK4PFchs";
    //std::string jerSF = jerTag+"_ScaleFactor_AK4PFchs";

    //actually, jet smearing/scaling won't work this way; have to do it like for electrons below
    //if ( do_jer ) c_pfjet_all    -> MatchAndSmearEnergy <PFJet   , GenJet     > ( c_genJet_final, 0.4/2.0, rootEngine, v_delta_met, cset_jmeAll.GetCorrection(jerPtRes).get(), cset_jmeAll.GetCorrection(jerSF).get() );
    //if ( do_jes ) c_pfjet_all    -> ScaleEnergy <PFJet   > ( pfjet_energy_scale_sign   , v_delta_met );
    jetUncertainties.ComputeJetVariations(c_pfjet_final_ptCut, c_genJet_all, !isData());
    //c_pfjet_final_ptCut->SetSystematics(c_pfjet_final->GetSystematicsNames(), c_pfjet_final->GetSystematics());

    map<string, double> metSystematics = type1METUncertainties.ComputeType1METVariations(c_pfjet_all, c_genJet_all, !isData());
    // get rid of "nominal" MET variations
    auto nominalItem = metSystematics.extract("Pt_nominal");
    nominalItem.key() = "Pt";
    metSystematics.insert(std::move(nominalItem));
    nominalItem = metSystematics.extract("Phi_nominal");
    nominalItem.key() = "Phi";
    metSystematics.insert(std::move(nominalItem));


    //-----------------------------------------------------------------
    // We need high-eta jets in order to look at boson recoil
    //-----------------------------------------------------------------
    CollectionPtr c_pfjet_highEta_ID;

    CollectionPtr c_pfjet_highEta                    = c_pfjet_all                          -> SkimByEtaRange   <PFJet>          ( -jet_HighEtaCut, jet_HighEtaCut   );
    if(jetIDType == "PFJetTight")
      c_pfjet_highEta_ID                 = c_pfjet_highEta                      -> SkimByID         <PFJet>          ( PFJET_TIGHT );    
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
    fillVariableWithValue( "event"    , readerTools_->ReadValueBranch<ULong64_t>("event")      );
    fillVariableWithValue( "ls"       , readerTools_->ReadValueBranch<UInt_t>("luminosityBlock")         );
    fillVariableWithValue( "run"      , readerTools_->ReadValueBranch<UInt_t>("run")        );
    //fillVariableWithValue( "ProcessID", ProcessID  );
    //fillVariableWithValue( "PtHat"    , PtHat      );
    float genWeight = -1.0;
    if(!isData()) {
      genWeight = readerTools_->ReadValueBranch<Float_t>("genWeight");
    }
    fillVariableWithValue( "Weight"   , genWeight   );
    // topPtWeights -- perhaps not needed since unused for 2016 analysis
    //fillVariableWithValue( "TopPtWeight",GenParticleTopPtWeight);
    // pileup
    float puWeight = 1.0;
    float puWeightDn = 1.0;
    float puWeightUp = 1.0;
    if(!isData()) {
      float nTrueInteractions = readerTools_->ReadValueBranch<Float_t>("Pileup_nTrueInt");
      puWeight =   cset_pu->evaluate({nTrueInteractions, "nominal"});
      puWeightUp = cset_pu->evaluate({nTrueInteractions, "up"});
      puWeightDn = cset_pu->evaluate({nTrueInteractions, "down"});  
      if(puWeight==0)
        std::cout << "Got puWeight = " << puWeight << "; run: " << getVariableValue<unsigned int>("run") << " ls: " << getVariableValue<unsigned int>("ls") << " event: " << getVariableValue<long long unsigned int>("event") << std::endl;
    }
    fillVariableWithValue( "puWeight"   , puWeight   );
    fillVariableWithValue( "puWeight_Up"   , puWeightUp   );
    fillVariableWithValue( "puWeight_Dn"   , puWeightDn   );
    // L1 prefiring
    float prefireDefault = 1.0;
    if(!isData()) {
      fillVariableWithValue("PrefireWeight",readerTools_->ReadValueBranch<Float_t>("L1PreFiringWeight_Nom"));
      fillVariableWithValue("PrefireWeight_Dn",readerTools_->ReadValueBranch<Float_t>("L1PreFiringWeight_Dn"));
      fillVariableWithValue("PrefireWeight_Up",readerTools_->ReadValueBranch<Float_t>("L1PreFiringWeight_Up"));
    }
    else {
      fillVariableWithValue("PrefireWeight", prefireDefault);
      fillVariableWithValue("PrefireWeight_Dn", prefireDefault);
      fillVariableWithValue("PrefireWeight_Up", prefireDefault);
    }

    //-----------------------------------------------------------------
    // ZVtx 2017
    //-----------------------------------------------------------------
    float zVtxSF = 1.0;
    if(analysisYearInt==2017)
      zVtxSF = ElectronScaleFactors2017::zVtxSF;
    fillVariableWithValue("ZVtxSF", zVtxSF);

    //-----------------------------------------------------------------
    // Pass JSON
    //-----------------------------------------------------------------
    fillVariableWithValue("PassJSON"                   , passJSON(getVariableValue<unsigned int>("run"), getVariableValue<unsigned int>("ls"), isData())                       );    

    //-----------------------------------------------------------------
    // Fill MET filter values
    // https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    //-----------------------------------------------------------------
    fillVariableWithValue("PassGlobalSuperTightHalo2016Filter" , readerTools_->ReadValueBranch<Bool_t>("Flag_globalSuperTightHalo2016Filter"));
    fillVariableWithValue("PassGoodVertices"              , readerTools_->ReadValueBranch<Bool_t>("Flag_goodVertices"));
    fillVariableWithValue("PassHBHENoiseFilter"           , readerTools_->ReadValueBranch<Bool_t>("Flag_HBHENoiseFilter"));
    fillVariableWithValue("PassHBHENoiseIsoFilter"        , readerTools_->ReadValueBranch<Bool_t>("Flag_HBHENoiseIsoFilter"));
    fillVariableWithValue("PassBadEESupercrystalFilter"   , readerTools_->ReadValueBranch<Bool_t>("Flag_eeBadScFilter"));
    fillVariableWithValue("PassEcalDeadCellTrigPrim"      , readerTools_->ReadValueBranch<Bool_t>("Flag_EcalDeadCellTriggerPrimitiveFilter"));
    std::string branchName = "Flag_BadChargedCandidateFilter";
    fillVariableWithValue("PassChargedCandidateFilter"    , readerTools_->ReadValueBranch<Bool_t>(branchName));
    fillVariableWithValue("PassBadPFMuonFilter"           , readerTools_->ReadValueBranch<Bool_t>("Flag_BadPFMuonFilter"));
    // for 2017 and 2018 only
    if(hasBranch("Flag_ecalBadCalibFilterV2"))
      fillVariableWithValue("PassEcalBadCalibV2Filter"    , readerTools_->ReadValueBranch<Bool_t>(branchName));
    else
      fillVariableWithValue("PassEcalBadCalibV2Filter"    , true);

    //-----------------------------------------------------------------
    // Fill MET values
    //-----------------------------------------------------------------
    if(reducedSkimType != 0) {
      // includes nominal
      for(const auto& syst : metSystematics)
        fillVariableWithValue("PFMET_Type1_"+syst.first, syst.second);
    }
    else {
      fillVariableWithValue("PFMET_Type1_Pt", metSystematics["Pt"]);
      fillVariableWithValue("PFMET_Type1_Phi", metSystematics["Phi"]);
    }

    if ( !isData() ) { 
      if ( reducedSkimType != 0 ){ 
        fillVariableWithValue("GenMET_Pt"		, readerTools_->ReadValueBranch<Float_t>("GenMET_pt"));
        fillVariableWithValue("GenMET_Phi"	, readerTools_->ReadValueBranch<Float_t>("GenMET_phi"));
        // add LHE variables if needed
        if(hasBranch("LHE_Vpt")) {
          fillVariableWithValue("LHE_Vpt"	  , readerTools_->ReadValueBranch<Float_t>("LHE_Vpt"));
          fillVariableWithValue("LHE_NpLO"	, static_cast<unsigned int>(readerTools_->ReadValueBranch<UChar_t>("LHE_NpLO")));
          fillVariableWithValue("LHE_NpNLO"	, static_cast<unsigned int>(readerTools_->ReadValueBranch<UChar_t>("LHE_NpNLO")));
          fillVariableWithValue("LHE_Njets"	, static_cast<unsigned int>(readerTools_->ReadValueBranch<UChar_t>("LHE_Njets")));
          fillVariableWithValue("LHE_Nglu"	, static_cast<unsigned int>(readerTools_->ReadValueBranch<UChar_t>("LHE_Nglu")));
        }
        if(hasBranch("LHEPdfWeight"))
          fillArrayVariableWithValue("LHEPdfWeight"	, readerTools_->ReadArrayBranch<Float_t>("LHEPdfWeight"));
        if(hasBranch("LHEScaleWeight"))
          fillArrayVariableWithValue("LHEScaleWeight"	, readerTools_->ReadArrayBranch<Float_t>("LHEScaleWeight"));
      }
    }

    //-----------------------------------------------------------------
    // Fill pileup variables
    //-----------------------------------------------------------------

    fillVariableWithValue( "nVertex", readerTools_->ReadValueBranch<Int_t>("PV_npvs"));
    float puNTrueInt = -1.0;
    if(!isData()) {
      puNTrueInt = readerTools_->ReadValueBranch<Float_t>("Pileup_nTrueInt");
    }
    fillVariableWithValue( "nPileUpInt_True", puNTrueInt);

    //-----------------------------------------------------------------
    // How many ID'd objects are there?
    //-----------------------------------------------------------------

    int n_muonLoose          = c_muon_eta_IDLoose            -> GetSize();
    int n_muonHighPt         = c_muon_eta_IDHighPt           -> GetSize();
    int n_muon_store         = c_muon_final                  -> GetSize();
    int n_ele_store          = c_ele_final                   -> GetSize();
    int n_jet_store          = c_pfjet_final_ptCut           -> GetSize();
    int n_jet_highEta_store  = c_pfjet_highEta_final         -> GetSize();
    int n_genEle_store       = c_genEle_final                -> GetSize();
    int n_genNu_store        = c_genNu_final                 -> GetSize();
    int n_genMu_store        = c_genMu_final                 -> GetSize();
    int n_genJet_store       = c_genJet_final                -> GetSize();

    int n_muon_ptCut         = c_muon_final_ptCut            -> GetSize();
    int n_ele_ptCut          = c_ele_final_ptCut             -> GetSize();
    int n_ele_vloose_ptCut   = c_ele_vLoose_ptCut            -> GetSize();
    int n_jet_ptCut          = c_pfjet_final_ptCut           -> GetSize();
    int n_jet_highEta_ptCut  = c_pfjet_highEta_final_ptCut   -> GetSize();

    int n_genNuFromW_store   = c_genNuFromW_final            -> GetSize();

    int n_genLQ_store = c_genLQ_final ->GetSize();

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

      fillVariableWithValue("nGenNuFromW_ptCut"	 , n_genNuFromW_store   );

      fillVariableWithValue("nGenNuFromW_store"	 , min(n_genNuFromW_store  ,2));

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
              //else {
              //  std::cout << "This event: " <<
              //    static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(luminosityBlock) << " " << static_cast<unsigned int>(event) <<
              //    " had no 5th gen Jet; examine other GenJets." << std::endl;
              //  c_genJet_final->examine<GenJet>("finalGenJets");
              //}
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

      if ( n_genLQ_store >= 1 ){ 
        GenParticle genLQ1 = c_genLQ_final -> GetConstituent<GenParticle>(0);
        fillVariableWithValue ( "GenLQ1_Pt" , genLQ1.Pt () );
        fillVariableWithValue ( "GenLQ1_Eta", genLQ1.Eta() );
        fillVariableWithValue ( "GenLQ1_Phi", genLQ1.Phi() );
        fillVariableWithValue ( "GenLQ1_Mass", genLQ1.Mass() );
        fillVariableWithValue ( "GenLQ1_ID" , genLQ1.PdgId());

        if ( n_genLQ_store >= 2 ){ 
          GenParticle genLQ2 = c_genLQ_final -> GetConstituent<GenParticle>(1);
          fillVariableWithValue ( "GenLQ2_Pt" , genLQ2.Pt () );
          fillVariableWithValue ( "GenLQ2_Eta", genLQ2.Eta() );
          fillVariableWithValue ( "GenLQ2_Phi", genLQ2.Phi() );
          fillVariableWithValue ( "GenLQ2_Mass", genLQ2.Mass() );
          fillVariableWithValue ( "GenLQ2_ID" , genLQ2.PdgId());
        }
      }

    }

    //-----------------------------------------------------------------
    // All skims need muons
    //-----------------------------------------------------------------

    fillVariableWithValue ("nMuon_ptCut", n_muon_ptCut);
    fillVariableWithValue ("nMuon_LooseId", n_muonLoose);
    fillVariableWithValue ("nMuon_HighPtId", n_muonHighPt);
    fillVariableWithValue ("nMuon_store", min(n_muon_store,3));

    if ( n_muon_store >= 1 ){ 

      Muon muon1 = c_muon_final -> GetConstituent<Muon>(0);
      //double hltSingleMuon1Pt = triggerMatchPt<HLTriggerObject, Muon>(c_hltMuon_SingleMu_all, muon1, muon_hltMatch_DeltaRCut);
      fillVariableWithValue ("Muon1_Pt"             , muon1.Pt      ());
      fillVariableWithValue ("Muon1_Eta"            , muon1.Eta     ());
      fillVariableWithValue ("Muon1_Phi"            , muon1.Phi     ());
      fillVariableWithValue ("Muon1_PtError"        , muon1.PtError ());
      fillVariableWithValue ("Muon1_Charge"         , muon1.Charge  ());
      //fillVariableWithValue ("Muon1_hltSingleMuonPt", hltSingleMuon1Pt);

      if ( n_muon_store >= 2 ){ 

        Muon muon2 = c_muon_final -> GetConstituent<Muon>(1);
        //double hltSingleMuon2Pt = triggerMatchPt<HLTriggerObject, Muon>(c_hltMuon_SingleMu_all, muon2, muon_hltMatch_DeltaRCut);
        fillVariableWithValue ("Muon2_Pt"             , muon2.Pt      ());
        fillVariableWithValue ("Muon2_Eta"            , muon2.Eta     ());
        fillVariableWithValue ("Muon2_Phi"            , muon2.Phi     ());
        fillVariableWithValue ("Muon2_PtError"        , muon2.PtError ());
        fillVariableWithValue ("Muon2_Charge"         , muon2.Charge  ());
        //fillVariableWithValue ("Muon2_hltSingleMuonPt", hltSingleMuon2Pt);

        if ( n_muon_store >= 3 ){ 

          Muon muon3 = c_muon_final -> GetConstituent<Muon>(2);
          //double hltSingleMuon3Pt = triggerMatchPt<HLTriggerObject, Muon>(c_hltMuon_SingleMu_all, muon3, muon_hltMatch_DeltaRCut);
          fillVariableWithValue ("Muon3_Pt"             , muon3.Pt      ());
          fillVariableWithValue ("Muon3_Eta"            , muon3.Eta     ());
          fillVariableWithValue ("Muon3_Phi"            , muon3.Phi     ());
          fillVariableWithValue ("Muon3_PtError"        , muon3.PtError ());
          fillVariableWithValue ("Muon3_Charge"         , muon3.Charge  ());
          //fillVariableWithValue ("Muon3_hltSingleMuonPt", hltSingleMuon3Pt);
        }
      }
    }

    //-----------------------------------------------------------------
    // Fill variables
    //-----------------------------------------------------------------

    // Electrons
    fillVariableWithValue ("nEle_store"       , min(n_ele_store,3) );
    fillVariableWithValue ("nEle_ID"          , n_ele_store);
    fillVariableWithValue ("nJet_store"       , min(n_jet_store,5) );
    fillVariableWithValue ("nJet_etaIdLepCleaned"       , n_jet_store);
    fillVariableWithValue ("nHighEtaJet_store", min(n_jet_highEta_store, 1));
    fillVariableWithValue ("nEle_ptCut"       , n_ele_ptCut );
    fillVariableWithValue ("nVLooseEle_ptCut"  , n_ele_vloose_ptCut );
    fillVariableWithValue ("nJet_ptCut"       , n_jet_ptCut );
    fillVariableWithValue ("nHighEtaJet_ptCut", n_jet_highEta_ptCut );

    //c_ele_all->examine<Electron>("c_ele_all");
    //c_ele_final->examine<Electron>("c_ele_final");
    if(getVariableValue<int>("nEle_store") < 1)
      fillVariableWithValue( "Ele1_Pt", -1.0);
    for(int iEle = 0; iEle < getVariableValue<int>("nEle_store"); ++iEle) {
      Electron ele = c_ele_final -> GetConstituent<Electron>(iEle);
      std::string prefix = "Ele"+std::to_string(iEle+1);
      unsigned int n_trigObjs = c_hlTriggerObjects_id -> GetSize();
      float matchedHLTriggerObjectPt = -999.;
      bool passedHLTriggerWPTightFilter = false;
      bool passedHLTriggerCaloIdVTGsfTrkIdTFilter = false;
      bool passedHLTriggerHighPtPhotonFilter = false;
      bool passedHLTriggerSinglePhotonFilter = false;
      if ( n_trigObjs > 0 ) {
        //c_hlTriggerObjects_id->examine<HLTriggerObject>("c_hlTriggerObjects_id");
        HLTriggerObject matchedObject;
        bool matched = ele.template MatchByDRAndDPt<HLTriggerObject>(c_hlTriggerObjects_id, matchedObject, ele_hltMatch_DeltaRCut, 0.5);
        //bool matched = ele.template MatchByDR<HLTriggerObject>(c_hlTriggerObjects_id, matchedObject, ele_hltMatch_DeltaRCut);
        if(matched) {
          matchedHLTriggerObjectPt = matchedObject.Pt();
          float dR = matchedObject.DeltaR(&ele);
          float dPt = matchedObject.DeltaPt(&ele);
          passedHLTriggerWPTightFilter = matchedObject.ObjectID()==11 && matchedObject.PassedFilterBit(1);
          passedHLTriggerCaloIdVTGsfTrkIdTFilter = matchedObject.ObjectID()==11 && matchedObject.PassedFilterBit(11);
          passedHLTriggerHighPtPhotonFilter = matchedObject.ObjectID()==11 && matchedObject.PassedFilterBit(13);
          bool passedSinglePhotonFilter = matchedObject.PassedFilterBit(1) || matchedObject.PassedFilterBit(2) || matchedObject.PassedFilterBit(3) || matchedObject.PassedFilterBit(4) || matchedObject.PassedFilterBit(5) || matchedObject.PassedFilterBit(6) || matchedObject.PassedFilterBit(7);
          passedHLTriggerSinglePhotonFilter = matchedObject.ObjectID()==22 && passedSinglePhotonFilter;
          //std::cout << "\t--> Matched a trigger object with pT=" << matchedHLTriggerObjectPt << " to electron; dR=" << dR << ", dPt=" << dPt << std::endl <<
          //  "\t\t" << matchedObject << "; passedWPTight=" << passedHLTriggerWPTightFilter << ", passedCaloId=" << passedHLTriggerCaloIdVTGsfTrkIdTFilter << 
          //  ", passedHighPtPhoton=" << passedHLTriggerHighPtPhotonFilter << ", passedHLTriggerSinglePhotonFilter=" << passedHLTriggerSinglePhotonFilter << 
          //  ", filterBits=" << matchedObject.FilterBits() << std::endl;
          c_hlTriggerObjects_id->RemoveConstituent(matchedObject); // remove so that we can't match it with another electron
        }
      }
      //std::cout << "\t" << ele << std::endl;
      if(reducedSkimType == 0)
        fillVariableWithValue( prefix+"_Pt"            , ele.PtUncorr()           );
      else
        fillVariableWithValue( prefix+"_Pt"            , ele.Pt()                 );
      fillVariableWithValue( prefix+"_SCEt"          , ele.SCEnergy()/cosh(ele.SCEta()) );
      fillVariableWithValue( prefix+"_PassHEEPID"    , ele.PassUserID ( heepIdType )  );
      fillVariableWithValue( prefix+"_PassEGMLooseID", ele.PassEGammaIDLoose()  );
      fillVariableWithValue( prefix+"_PassEGMMediumID", ele.PassEGammaIDMedium()  );
      fillVariableWithValue( prefix+"_ECorr"         , ele.ECorr()              );
      fillVariableWithValue( prefix+"_Eta"           , ele.Eta()                );
      fillVariableWithValue( prefix+"_Phi"           , ele.Phi()                );
      fillVariableWithValue( prefix+"_SCEta"         , ele.SCEta()              );
      fillVariableWithValue( prefix+"_Charge"        , ele.Charge()             );
      fillVariableWithValue( prefix+"_R9"            , ele.R9()                 );
      fillVariableWithValue( prefix+"_MissingHits"   , ele.MissingHits()        );
      fillVariableWithValue( prefix+"_Full5x5SigmaIEtaIEta" , ele.Full5x5SigmaIEtaIEta() );
      fillVariableWithValue( prefix+"_RhoForHEEP"    , ele.RhoForHEEP()         );

      fillVariableWithValue( prefix+"_DeltaEtaTrkSC" , ele.DeltaEta()           );
      fillVariableWithValue( prefix+"_HoE"           , ele.HoE()                );
      fillVariableWithValue( prefix+"_HasMatchedPhot", ele.HasMatchedConvPhot() );
      fillVariableWithValue( prefix+"_LeadVtxDistXY" , ele.LeadVtxDistXY()      );
      fillVariableWithValue( prefix+"_LeadVtxDistZ"  , ele.LeadVtxDistZ ()      );
      fillVariableWithValue( prefix+"_TrkIsolation"  , ele.TrkIsoDR03()         );
      fillVariableWithValue( prefix+"_TrkIsoHEEP7"   , ele.HEEP70TrackIsolation());
      fillVariableWithValue( prefix+"_EcalIsolation" , ele.EcalIsoDR03()        );
      fillVariableWithValue( prefix+"_HcalIsolation" , ele.HcalIsoD1DR03()      );
      fillVariableWithValue( prefix+"_CorrIsolation" , ele.HEEPCorrIsolation()  );
      fillVariableWithValue( prefix+"_PFRelIsoCh03"  , ele.PFRelIso03Charged());
      fillVariableWithValue( prefix+"_PFRelIsoAll03" , ele.PFRelIso03All());
      // HEEP ID cuts
      fillVariableWithValue( prefix+"_PassHEEPMinPtCut"                            ,ele.PassHEEPMinPtCut                            () );
      fillVariableWithValue( prefix+"_PassHEEPGsfEleSCEtaMultiRangeCut"            ,ele.PassHEEPGsfEleSCEtaMultiRangeCut            () ); 
      fillVariableWithValue( prefix+"_PassHEEPGsfEleDEtaInSeedCut"                 ,ele.PassHEEPGsfEleDEtaInSeedCut                 () ); 
      fillVariableWithValue( prefix+"_PassHEEPGsfEleDPhiInCut"                     ,ele.PassHEEPGsfEleDPhiInCut                     () ); 
      fillVariableWithValue( prefix+"_PassHEEPGsfEleFull5x5SigmaIEtaIEtaWithSatCut",ele.PassHEEPGsfEleFull5x5SigmaIEtaIEtaWithSatCut() ); 
      fillVariableWithValue( prefix+"_PassHEEPGsfEleFull5x5E2x5OverE5x5WithSatCut" ,ele.PassHEEPGsfEleFull5x5E2x5OverE5x5WithSatCut () ); 
      fillVariableWithValue( prefix+"_PassHEEPGsfEleTrkPtIsoCut"                   ,ele.PassHEEPGsfEleTrkPtIsoCut                   () ); 
      if(analysisYearInt == 2018) {
        fillVariableWithValue( prefix+"_PassHEEPGsfEleEmHadD1IsoRhoCut"              ,ele.PassHEEPGsfEleEmHadD1IsoRhoCut2018        () ); 
        fillVariableWithValue( prefix+"_PassHEEPGsfEleHadronicOverEMLinearCut"       ,ele.PassHEEPGsfEleHadronicOverEMLinearCut2018 () ); 
      }
      else {
        fillVariableWithValue( prefix+"_PassHEEPGsfEleEmHadD1IsoRhoCut"              ,ele.PassHEEPGsfEleEmHadD1IsoRhoCut              () ); 
        fillVariableWithValue( prefix+"_PassHEEPGsfEleHadronicOverEMLinearCut"       ,ele.PassHEEPGsfEleHadronicOverEMLinearCut       () ); 
      }
      fillVariableWithValue( prefix+"_PassHEEPGsfEleDxyCut"                        ,ele.PassHEEPGsfEleDxyCut                        () ); 
      fillVariableWithValue( prefix+"_PassHEEPGsfEleMissingHitsCut"                ,ele.PassHEEPGsfEleMissingHitsCut                () ); 
      fillVariableWithValue( prefix+"_PassHEEPEcalDrivenCut"                       ,ele.PassHEEPEcalDrivenCut                       () );
      // EGM Loose ID cuts
      fillVariableWithValue( prefix+"_PassEGammaLooseGsfEleFull5x5SigmaIEtaIEtaCut" ,ele.PassEGammaIDLooseGsfEleFull5x5SigmaIEtaIEtaCut () );
      fillVariableWithValue( prefix+"_PassEGammaLooseGsfEleDEtaInSeedCut"           ,ele.PassEGammaIDLooseGsfEleDEtaInSeedCut           () );
      fillVariableWithValue( prefix+"_PassEGammaLooseGsfEleDPhiInCut"               ,ele.PassEGammaIDLooseGsfEleDPhiInCut               () );
      fillVariableWithValue( prefix+"_PassEGammaLooseGsfEleHadronicOverEMScaledCut" ,ele.PassEGammaIDLooseGsfEleHadronicOverEMScaledCut () );
      fillVariableWithValue( prefix+"_PassEGammaLooseGsfEleEInverseMinusPInverseCut",ele.PassEGammaIDLooseGsfEleEInverseMinusPInverseCut() );
      fillVariableWithValue( prefix+"_PassEGammaLooseGsfEleRelPFIsoScaledCut"       ,ele.PassEGammaIDLooseGsfEleRelPFIsoScaledCut       () );
      fillVariableWithValue( prefix+"_PassEGammaLooseGsfEleConversionVetoCut"       ,ele.PassEGammaIDLooseGsfEleConversionVetoCut       () );
      fillVariableWithValue( prefix+"_PassEGammaLooseGsfEleMissingHitsCut"          ,ele.PassEGammaIDLooseGsfEleMissingHitsCut          () );

      fillVariableWithValue( prefix+"_MatchedHLTriggerObjectPt"  , matchedHLTriggerObjectPt );
      fillVariableWithValue( prefix+"_PassedHLTriggerWPTightFilter"  , passedHLTriggerWPTightFilter );
      fillVariableWithValue( prefix+"_PassedHLTriggerCaloIdVTGsfTrkIdTFilter"  , passedHLTriggerCaloIdVTGsfTrkIdTFilter );
      fillVariableWithValue( prefix+"_PassedHLTriggerHighPtPhotonFilter"  , passedHLTriggerHighPtPhotonFilter );
      fillVariableWithValue( prefix+"_PassedHLTriggerSinglePhotonFilter"  , passedHLTriggerSinglePhotonFilter );

      if(!isData() && reducedSkimType != 0) {
        fillVariableWithValue( prefix+"_Pt_EER_Up"  , ele.PtDESigmaUp()          );
        fillVariableWithValue( prefix+"_Pt_EER_Dn"  , ele.PtDESigmaDown()        );
        fillVariableWithValue( prefix+"_Pt_EES_Up"  , ele.Pt() * cset_ulEGMScale->evaluate({analysisYear, "scaleup",   ele.SCEta(), ele.SeedGain()}) );
        fillVariableWithValue( prefix+"_Pt_EES_Dn"  , ele.Pt() * cset_ulEGMScale->evaluate({analysisYear, "scaledown", ele.SCEta(), ele.SeedGain()}) );
      }
      if(!isData()) {
        fillVariableWithValue( prefix+"_TrigSF" , triggerScaleFactorReader.LookupValue(ele.SCEta(), ele.Pt()) );
        fillVariableWithValue( prefix+"_TrigSF_Err" , triggerScaleFactorReader.LookupValueError(ele.SCEta(), ele.Pt()) );
        //if(ele.PtUncorr() >= 10) {
        if(ele.Pt() >= 10) {
          float recoSF = ele.Pt() >= 20 ? cset_ulEGMID->evaluate({analysisYear, "sf", "RecoAbove20", ele.SCEta(), ele.Pt()}) : cset_ulEGMID->evaluate({analysisYear, "sf", "RecoBelow20", ele.SCEta(), ele.Pt()});
          float recoSFUp = ele.Pt() >= 20 ? cset_ulEGMID->evaluate({analysisYear, "sfup", "RecoAbove20", ele.SCEta(), ele.Pt()}) : cset_ulEGMID->evaluate({analysisYear, "sfup", "RecoBelow20", ele.SCEta(), ele.Pt()});
          float recoSFDn = ele.Pt() >= 20 ? cset_ulEGMID->evaluate({analysisYear, "sfdown", "RecoAbove20", ele.SCEta(), ele.Pt()}) : cset_ulEGMID->evaluate({analysisYear, "sfdown", "RecoBelow20", ele.SCEta(), ele.Pt()});
          fillVariableWithValue( prefix+"_RecoSF",    recoSF );
          fillVariableWithValue( prefix+"_RecoSF_Up", recoSFUp );
          fillVariableWithValue( prefix+"_RecoSF_Dn", recoSFDn );
          fillVariableWithValue( prefix+"_EGMLooseIDSF",    cset_ulEGMID->evaluate({analysisYear, "sf", "Loose", ele.SCEta(), ele.Pt()}));
          fillVariableWithValue( prefix+"_EGMLooseIDSF_Up", cset_ulEGMID->evaluate({analysisYear, "sfup", "Loose", ele.SCEta(), ele.Pt()}));
          fillVariableWithValue( prefix+"_EGMLooseIDSF_Dn", cset_ulEGMID->evaluate({analysisYear, "sfdown", "Loose", ele.SCEta(), ele.Pt()}));
        }
        fillVariableWithValue( prefix+"_HEEPSF",    ElectronScaleFactorsRunII::LookupHeepSF(ele.SCEta(), analysisYear) );
        fillVariableWithValue( prefix+"_HEEPSF_Err",ElectronScaleFactorsRunII::LookupHeepSFSyst(ele.SCEta(), ele.Pt(), analysisYear) );
      }
      else {
        fillVariableWithValue( prefix+"_TrigSF"          , float(1.0));
        fillVariableWithValue( prefix+"_TrigSF_Err"      , float(0.0));
        fillVariableWithValue( prefix+"_RecoSF"          , float(1.0));
        fillVariableWithValue( prefix+"_RecoSF_Up"       , float(0.0));
        fillVariableWithValue( prefix+"_RecoSF_Dn"       , float(0.0));
        fillVariableWithValue( prefix+"_EGMLooseIDSF"    , float(1.0));
        fillVariableWithValue( prefix+"_EGMLooseIDSF_Up" , float(0.0));
        fillVariableWithValue( prefix+"_EGMLooseIDSF_Dn" , float(0.0));
        fillVariableWithValue( prefix+"_HEEPSF"          , float(1.0));
        fillVariableWithValue( prefix+"_HEEPSF_Err"      , float(0.0));
      }
    }
    // trigger scale factor part
    // NB: need to change this if the trigger selection changes below
    float eventTriggerScaleFactor = 1.0;
    float eventTriggerScaleFactorErr = 0.0;
    bool ele1PassedHLTWPTight = variableIsFilled("Ele1_PassedHLTriggerWPTightFilter") ? getVariableValue<bool>("Ele1_PassedHLTriggerWPTightFilter") : false;
    bool ele1PassedHLTPhoton = variableIsFilled("Ele1_PassedHLTriggerHighPtPhotonFilter") ? getVariableValue<bool>("Ele1_PassedHLTriggerHighPtPhotonFilter") : false;
    bool ele2PassedHLTWPTight = variableIsFilled("Ele2_PassedHLTriggerWPTightFilter") ? getVariableValue<bool>("Ele2_PassedHLTriggerWPTightFilter") : false;
    bool ele2PassedHLTPhoton = variableIsFilled("Ele2_PassedHLTriggerHighPtPhotonFilter") ? getVariableValue<bool>("Ele2_PassedHLTriggerHighPtPhotonFilter") : false;
    if(ele1PassedHLTWPTight || ele1PassedHLTPhoton) {
      float trigSFEle1 = getVariableValue("Ele1_TrigSF");
      float trigSFEle1Err = getVariableValue("Ele1_TrigSF_Err");
      eventTriggerScaleFactor = trigSFEle1;
      eventTriggerScaleFactorErr = trigSFEle1Err;
    }
    else if(ele2PassedHLTWPTight || ele2PassedHLTPhoton) {
      float trigSFEle2 = getVariableValue("Ele2_TrigSF");
      float trigSFEle2Err = getVariableValue("Ele2_TrigSF_Err");
      eventTriggerScaleFactor = trigSFEle2;
      eventTriggerScaleFactorErr = trigSFEle2Err;
    }
    fillVariableWithValue("EventTriggerScaleFactor", eventTriggerScaleFactor);
    fillVariableWithValue("EventTriggerScaleFactorErr", eventTriggerScaleFactorErr);

    // Jets
    if ( n_jet_highEta_store >= 1 ) { 
      PFJet jet1 = c_pfjet_highEta_final -> GetConstituent<PFJet>(0);
      fillVariableWithValue( "HighEtaJet1_Pt" , jet1.Pt () );
      fillVariableWithValue( "HighEtaJet1_Eta", jet1.Eta() );
      fillVariableWithValue( "HighEtaJet1_Phi", jet1.Phi() );
    }

    std::vector<std::string> jetSystNames = c_pfjet_final_ptCut->GetSystematicsNames();
    for(int iJet = 0; iJet < getVariableValue<int>("nJet_store"); ++iJet) {
      PFJet jet = c_pfjet_final_ptCut -> GetConstituent<PFJet>(iJet);
      std::string prefix = "Jet"+std::to_string(iJet+1);
      int flavor = isData() ? 0 : jet.HadronFlavor();
      float pt = jet.GetSystematicVariation("Pt_nominal");
      fillVariableWithValue( prefix+"_Pt"          , pt                         );
      fillVariableWithValue( prefix+"_PtOriginal"  , jet.PtOriginal()            );
      fillVariableWithValue( prefix+"_Eta"         , jet.Eta()                  );
      fillVariableWithValue( prefix+"_Phi"         , jet.Phi()                  );
      fillVariableWithValue( prefix+"_btagDeepJet" , jet.DeepJetBTag()          );
      fillVariableWithValue( prefix+"_HadronFlavor", flavor                     );
      //fillVariableWithValue( prefix+"_qgl"	  , jet1.QuarkGluonLikelihood()     );
      float absEta = fabs(jet.Eta());
      // b-tagging scale factors and uncertainties
      if(!isData() && reducedSkimType != 0) {
        if(flavor == 0) {
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetIncl"  , cset_btagDeepJetIncl->evaluate({"central", "L", flavor, absEta, pt} ));
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetIncl_Up"  , cset_btagDeepJetIncl->evaluate({"up", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetIncl_UpCorrelated"  , cset_btagDeepJetIncl->evaluate({"up_correlated", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetIncl_UpUncorrelated"  , cset_btagDeepJetIncl->evaluate({"up_uncorrelated", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetIncl_Dn"  , cset_btagDeepJetIncl->evaluate({"down", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetIncl_DnCorrelated"  , cset_btagDeepJetIncl->evaluate({"down_correlated", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetIncl_DnUncorrelated"  , cset_btagDeepJetIncl->evaluate({"down_uncorrelated", "L", flavor, absEta, pt}) );
          // medium
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetIncl"  , cset_btagDeepJetIncl->evaluate({"central", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetIncl_Up"  , cset_btagDeepJetIncl->evaluate({"up", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetIncl_UpCorrelated"  , cset_btagDeepJetIncl->evaluate({"up_correlated", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetIncl_UpUncorrelated"  , cset_btagDeepJetIncl->evaluate({"up_uncorrelated", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetIncl_Dn"  , cset_btagDeepJetIncl->evaluate({"down", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetIncl_DnCorrelated"  , cset_btagDeepJetIncl->evaluate({"down_correlated", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetIncl_DnUncorrelated"  , cset_btagDeepJetIncl->evaluate({"down_uncorrelated", "M", flavor, absEta, pt}) );
        }
        else if(flavor==4 || flavor==5) {
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetComb"  , cset_btagDeepJetComb->evaluate({"central", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetComb_Up"  , cset_btagDeepJetComb->evaluate({"up", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetComb_UpCorrelated"  , cset_btagDeepJetComb->evaluate({"up_correlated", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetComb_UpUncorrelated"  , cset_btagDeepJetComb->evaluate({"up_uncorrelated", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetComb_Dn"  , cset_btagDeepJetComb->evaluate({"down", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetComb_DnCorrelated"  , cset_btagDeepJetComb->evaluate({"down_correlated", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetComb_DnUncorrelated"  , cset_btagDeepJetComb->evaluate({"down_uncorrelated", "L", flavor, absEta, pt}) );
          // medium
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetComb"  , cset_btagDeepJetComb->evaluate({"central", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetComb_Up"  , cset_btagDeepJetComb->evaluate({"up", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetComb_UpCorrelated"  , cset_btagDeepJetComb->evaluate({"up_correlated", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetComb_UpUncorrelated"  , cset_btagDeepJetComb->evaluate({"up_uncorrelated", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetComb_Dn"  , cset_btagDeepJetComb->evaluate({"down", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetComb_DnCorrelated"  , cset_btagDeepJetComb->evaluate({"down_correlated", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetComb_DnUncorrelated"  , cset_btagDeepJetComb->evaluate({"down_uncorrelated", "M", flavor, absEta, pt}) );
        }
        if(!isData() && iJet < 3) { // don't fill the systs for jets beyond the third
          for(int varIndex = 1; varIndex < jetSystNames.size(); ++varIndex) {
            string systName = jetSystNames[varIndex];
            if(systName.find("mass") != string::npos)
              continue;
            fillVariableWithValue( prefix+"_"+systName  , jet.GetSystematicVariation(systName));
          }
        }
      }
      if(c_pfJetMatchedLQ->Has<PFJet>(jet))
        fillVariableWithValue( prefix+"_LQMatched", true);
    }

    //-----------------------------------------------------------------
    // Fill variables that depend on more than one object
    // All skims need this
    //-----------------------------------------------------------------
    // strip off the beginning part for the 
    vector<string> compSystNames;
    for(const auto& systName : jetSystNames) {
      if(systName.find("mass") != string::npos)
        continue;
      compSystNames.push_back(systName.substr(3));
    }

    TLorentzVector t_ele1, t_ele2, t_jet1, t_jet2, t_jet3;
    TLorentzVector t_MET;
    TLorentzVector t_ele1SigmaUp, t_ele1SigmaDown, t_ele1ScaledUp, t_ele1ScaledDown;
    TLorentzVector t_ele2SigmaUp, t_ele2SigmaDown, t_ele2ScaledUp, t_ele2ScaledDown;
    std::vector<TLorentzVector> jet1FourVectors_ptVariations, jet2FourVectors_ptVariations, jet3FourVectors_ptVariations;

    t_MET.SetPtEtaPhiM( readerTools_->ReadValueBranch<Float_t>("MET_pt"), 0.0, readerTools_->ReadValueBranch<Float_t>("MET_phi"), 0.0 );

    if ( n_jet_store >= 1 ){

      PFJet jet1 = c_pfjet_final_ptCut -> GetConstituent<PFJet>(0);
      float jet1Pt = jet1.GetSystematicVariation("Pt_nominal");
      t_jet1.SetPtEtaPhiM ( jet1Pt, jet1.Eta(), jet1.Phi(), 0.0 );
      if(!isData()) {
        for(int varIndex = 1; varIndex < compSystNames.size(); ++varIndex) {
          string systName = jetSystNames[varIndex];
          TLorentzVector jet4Vec;
          jet4Vec.SetPtEtaPhiM(jet1.GetSystematicVariation(systName), jet1.Eta(), jet1.Phi(), 0.0);
          jet1FourVectors_ptVariations.push_back(jet4Vec);
        }
      }
      fillVariableWithValue ("mDPhi_METJet1", fabs( t_MET.DeltaPhi ( t_jet1 )));

      if ( n_jet_store >= 2 ){
        PFJet jet2 = c_pfjet_final_ptCut -> GetConstituent<PFJet>(1);
        float jet2Pt = jet2.GetSystematicVariation("Pt_nominal");
        t_jet2.SetPtEtaPhiM ( jet2Pt, jet2.Eta(), jet2.Phi(), 0.0 );
        if(!isData()) {
          for(int varIndex = 1; varIndex < compSystNames.size(); ++varIndex) {
            string systName = jetSystNames[varIndex];
            TLorentzVector jet4Vec;
            jet4Vec.SetPtEtaPhiM(jet2.GetSystematicVariation(systName), jet2.Eta(), jet2.Phi(), 0.0);
            jet2FourVectors_ptVariations.push_back(jet4Vec);
          }
        }
        TLorentzVector t_jet1jet2 = t_jet1 + t_jet2;

        fillVariableWithValue ("M_j1j2" , t_jet1jet2.M ());
        fillVariableWithValue ("Pt_j1j2", t_jet1jet2.Pt());
        fillVariableWithValue ("mDPhi_METJet2", fabs( t_MET.DeltaPhi ( t_jet2 )));
        fillVariableWithValue ("DR_Jet1Jet2"  , t_jet1.DeltaR( t_jet2 ));

        if ( n_jet_store >= 3 ){
          PFJet jet3 = c_pfjet_final_ptCut -> GetConstituent<PFJet>(2);
          float jet3Pt = jet3.GetSystematicVariation("Pt_nominal");
          t_jet3.SetPtEtaPhiM ( jet3Pt, jet3.Eta(), jet3.Phi(), 0.0 );
          if(!isData()) {
            for(int varIndex = 1; varIndex < compSystNames.size(); ++varIndex) {
              string systName = jetSystNames[varIndex];  // to lookup the syst, use the original name
              TLorentzVector jet4Vec;
              jet4Vec.SetPtEtaPhiM(jet3.GetSystematicVariation(systName), jet3.Eta(), jet3.Phi(), 0.0);
              jet3FourVectors_ptVariations.push_back(jet4Vec);
            }
          }
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
        t_ele1.SetPtEtaPhiM ( ele1.PtUncorr(), ele1.Eta(), ele1.Phi(), 0.0 );
      if(!isData() && reducedSkimType != 0) {
        t_ele1SigmaUp.SetPtEtaPhiM(   getVariableValue("Ele1_Pt_EER_Up"),    ele1.Eta(), ele1.Phi(), 0.0);
        t_ele1SigmaDown.SetPtEtaPhiM( getVariableValue("Ele1_Pt_EER_Dn"),  ele1.Eta(), ele1.Phi(), 0.0);
        t_ele1ScaledUp.SetPtEtaPhiM(  getVariableValue("Ele1_Pt_EES_Up"),    ele1.Eta(), ele1.Phi(), 0.0);
        t_ele1ScaledDown.SetPtEtaPhiM( getVariableValue("Ele1_Pt_EES_Dn"), ele1.Eta(), ele1.Phi(), 0.0);
      }
      if ( n_ele_store >= 2 ) {
        Electron ele2 = c_ele_final -> GetConstituent<Electron>(1);
        t_ele2.SetPtEtaPhiM ( ele2.Pt(), ele2.Eta(), ele2.Phi(), 0.0 );
        if ( reducedSkimType == 0 ) // for QCD skims, use the uncorrected Pt
          t_ele2.SetPtEtaPhiM ( ele2.PtUncorr(), ele2.Eta(), ele2.Phi(), 0.0 );
        if(!isData() && reducedSkimType != 0) {
          t_ele2SigmaUp.SetPtEtaPhiM(   getVariableValue("Ele2_Pt_EER_Up"),    ele2.Eta(), ele2.Phi(), 0.0);
          t_ele2SigmaDown.SetPtEtaPhiM( getVariableValue("Ele2_Pt_EER_Dn"),  ele2.Eta(), ele2.Phi(), 0.0);
          t_ele2ScaledUp.SetPtEtaPhiM(  getVariableValue("Ele2_Pt_EES_Up"),    ele2.Eta(), ele2.Phi(), 0.0);
          t_ele2ScaledDown.SetPtEtaPhiM( getVariableValue("Ele2_Pt_EES_Dn"), ele2.Eta(), ele2.Phi(), 0.0);
        }

        TLorentzVector t_ele1ele2 = t_ele1 + t_ele2;
        TLorentzVector t_ele1ele2_eerUp = t_ele1SigmaUp + t_ele2SigmaUp;
        TLorentzVector t_ele1ele2_eerDown = t_ele1SigmaDown + t_ele2SigmaDown;
        TLorentzVector t_ele1ele2_eesUp = t_ele1ScaledUp + t_ele2ScaledUp;
        TLorentzVector t_ele1ele2_eesDown = t_ele1ScaledDown + t_ele2ScaledDown;
        fillVariableWithValue ("M_e1e2" , t_ele1ele2.M ());
        fillVariableWithValue ("Pt_e1e2", t_ele1ele2.Pt());
        //std::cout << "Pt_e1e2 Pt=" << t_ele1ele2.Pt() << std::endl;
        // systs
        if(!isData() && reducedSkimType != 0) {
          fillVariableWithValue ("M_e1e2_EER_Up" , t_ele1ele2_eerUp.M ());
          fillVariableWithValue ("M_e1e2_EER_Dn" , t_ele1ele2_eerDown.M ());
          fillVariableWithValue ("M_e1e2_EES_Up" , t_ele1ele2_eesUp.M ());
          fillVariableWithValue ("M_e1e2_EES_Dn" , t_ele1ele2_eesDown.M ());
          fillVariableWithValue ("Pt_e1e2_EER_Up" , t_ele1ele2_eerUp.Pt ());
          fillVariableWithValue ("Pt_e1e2_EER_Dn" , t_ele1ele2_eerDown.Pt ());
          fillVariableWithValue ("Pt_e1e2_EES_Up" , t_ele1ele2_eesUp.Pt ());
          fillVariableWithValue ("Pt_e1e2_EES_Dn" , t_ele1ele2_eesDown.Pt ());
          //std::cout << "Pt_e1e2_EER Pt=" << t_ele1ele2_eer.Pt() << std::endl;
          //std::cout << "Pt_e1e2_EES_Up Pt=" << t_ele1ele2_eesUp.Pt() << std::endl;
          //std::cout << "Pt_e1e2_EES_Down Pt=" << t_ele1ele2_eesDown.Pt() << std::endl;
        }
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
        // systs
        if(!isData() && reducedSkimType != 0) {
          fillVariableWithValue("M_e1j1_EER_Up"        , (t_ele1SigmaUp+t_jet1).M());
          fillVariableWithValue("M_e1j1_EER_Dn"        , (t_ele1SigmaDown+t_jet1).M());
          fillVariableWithValue("M_e1j1_EES_Up"     , (t_ele1ScaledUp+t_jet1).M());
          fillVariableWithValue("M_e1j1_EES_Dn"     , (t_ele1ScaledDown+t_jet1).M());
          for(int idx = 1; idx < compSystNames.size(); ++idx)
            fillVariableWithValue("M_e1j1_"+compSystNames[idx], (t_ele1+jet1FourVectors_ptVariations[idx]).M());
        }

        if ( n_jet_store >= 2 ){ 

          TLorentzVector t_ele1jet2 = t_ele1 + t_jet2;
          fillVariableWithValue("DR_Ele1Jet2", t_ele1.DeltaR ( t_jet2 ));
          fillVariableWithValue("M_e1j2"     , t_ele1jet2.M());
          // systs
          if(!isData() && reducedSkimType != 0) {
            fillVariableWithValue("M_e1j2_EER_Up"        , (t_ele1SigmaUp+t_jet2).M());
            fillVariableWithValue("M_e1j2_EER_Dn"        , (t_ele1SigmaDown+t_jet2).M());
            fillVariableWithValue("M_e1j2_EES_Up"     , (t_ele1ScaledUp+t_jet2).M());
            fillVariableWithValue("M_e1j2_EES_Dn"     , (t_ele1ScaledDown+t_jet2).M());
            for(int idx = 1; idx < compSystNames.size(); ++idx)
              fillVariableWithValue("M_e1j2_"+compSystNames[idx], (t_ele1+jet2FourVectors_ptVariations[idx]).M());
            fillVariableWithValue("sT_enujj"   , t_ele1.Pt() + t_MET.Pt() + t_jet1.Pt() + t_jet2.Pt());
          }
        }
      }
    }


    if ( n_ele_store >= 2 ){
      fillVariableWithValue("mDPhi_METEle2", fabs ( t_MET.DeltaPhi(t_ele2)));

      if ( n_jet_store >= 1 ){ 

        TLorentzVector t_ele2jet1 = t_ele2 + t_jet1;
        fillVariableWithValue("DR_Ele2Jet1", t_ele2.DeltaR ( t_jet1 ));
        fillVariableWithValue("M_e2j1"     , t_ele2jet1.M());
        // systs
        if(!isData() && reducedSkimType != 0) {
          fillVariableWithValue("M_e2j1_EER_Up"        , (t_ele2SigmaUp+t_jet1).M());
          fillVariableWithValue("M_e2j1_EER_Dn"        , (t_ele2SigmaDown+t_jet1).M());
          fillVariableWithValue("M_e2j1_EES_Up"     , (t_ele2ScaledUp+t_jet1).M());
          fillVariableWithValue("M_e2j1_EES_Dn"     , (t_ele2ScaledDown+t_jet1).M());
          for(int idx = 1; idx < compSystNames.size(); ++idx)
            fillVariableWithValue("M_e2j1_"+compSystNames[idx], (t_ele2+jet1FourVectors_ptVariations[idx]).M());
        }

        if ( n_jet_store >= 2 ){ 

          TLorentzVector t_ele2jet2 = t_ele2 + t_jet2;
          fillVariableWithValue("DR_Ele2Jet2", t_ele2.DeltaR ( t_jet2 ));
          fillVariableWithValue("M_e2j2"     , t_ele2jet2.M());
          fillVariableWithValue("sT_eejj"    , t_ele1.Pt() + t_ele2.Pt() + t_jet1.Pt() + t_jet2.Pt());
          // systs
          if(!isData() && reducedSkimType != 0) {
            fillVariableWithValue("M_e2j2_EER_Up"        , (t_ele2SigmaUp+t_jet2).M());
            fillVariableWithValue("M_e2j2_EER_Dn"        , (t_ele2SigmaDown+t_jet2).M());
            fillVariableWithValue("M_e2j2_EES_Up"     , (t_ele2ScaledUp+t_jet2).M());
            fillVariableWithValue("M_e2j2_EES_Dn"     , (t_ele2ScaledDown+t_jet2).M());
            for(int idx = 1; idx < compSystNames.size(); ++idx) {
              fillVariableWithValue("M_e2j2_"+compSystNames[idx], (t_ele2+jet2FourVectors_ptVariations[idx]).M());
              fillVariableWithValue("sT_eejj_"+compSystNames[idx], t_ele1.Pt() + t_ele2.Pt() + jet1FourVectors_ptVariations[idx].Pt() + jet2FourVectors_ptVariations[idx].Pt());
            }
            fillVariableWithValue("sT_eejj_EER_Up"    , t_ele1SigmaUp.Pt() + t_ele2SigmaUp.Pt() + t_jet1.Pt() + t_jet2.Pt());
            fillVariableWithValue("sT_eejj_EER_Dn"    , t_ele1SigmaDown.Pt() + t_ele2SigmaDown.Pt() + t_jet1.Pt() + t_jet2.Pt());
            fillVariableWithValue("sT_eejj_EES_Up"    , t_ele1ScaledUp.Pt() + t_ele2ScaledUp.Pt() + t_jet1.Pt() + t_jet2.Pt());
            fillVariableWithValue("sT_eejj_EES_Dn"    , t_ele1ScaledDown.Pt() + t_ele2ScaledDown.Pt() + t_jet1.Pt() + t_jet2.Pt());
          }
        }
      }
    }
    //-----------------------------------------------------------------
    // QCD triggers 2016    2017        2018
    // - Photon22         
    //                      Photon25
    // - Photon30           
    //                      Photon33    Photon33
    // - Photon36           
    // - Photon50           Photon50    Photon50
    // - Photon75           Photon75    Photon75
    // - Photon90           Photon90    Photon90
    // - Photon120          Photon120   Photon120
    //                      Photon150   Photon150
    // - Photon175          Photon175   Photon175
    // - Photon200          Photon200   Photon200
    //-----------------------------------------------------------------

    if ( reducedSkimType == 0 ) { 
      // NB: for data, all prescales applied as averages in downstream analysis code
      if(triggerExists("HLT_Photon22"))
        fillTriggerVariable ( "HLT_Photon22" , "H_Photon22"  );
      else
        fillVariableWithValue( "H_Photon22", -1.0); 
      if(triggerExists("HLT_Photon25"))
        fillTriggerVariable ( "HLT_Photon25" , "H_Photon25"  );
      else
        fillVariableWithValue( "H_Photon25", -1.0); 
      if(triggerExists("HLT_Photon30"))
        fillTriggerVariable ( "HLT_Photon30" , "H_Photon30"  );
      else
        fillVariableWithValue( "H_Photon30", -1.0); 
      if(triggerExists("HLT_Photon33"))
        fillTriggerVariable ( "HLT_Photon33" , "H_Photon33"  );
      else
        fillVariableWithValue( "H_Photon33", -1.0); 
      if(triggerExists("HLT_Photon36"))
        fillTriggerVariable ( "HLT_Photon36" , "H_Photon36"  );
      else
        fillVariableWithValue( "H_Photon36", -1.0); 
      fillTriggerVariable ( "HLT_Photon50"   , "H_Photon50"  );
      fillTriggerVariable ( "HLT_Photon75"   , "H_Photon75"  );
      fillTriggerVariable ( "HLT_Photon90"   , "H_Photon90"  );
      fillTriggerVariable ( "HLT_Photon120"  , "H_Photon120" );
      if(triggerExists("HLT_Photon150"))
        fillTriggerVariable ( "HLT_Photon150", "H_Photon150" );
      else
        fillVariableWithValue( "H_Photon150", -1.0); 
      fillTriggerVariable ( "HLT_Photon175" , "H_Photon175" );
      if(triggerExists("HLT_Photon200"))
        fillTriggerVariable ( "HLT_Photon200" , "H_Photon200" );
      else
        fillVariableWithValue ( "H_Photon200" , -1.0 );

      bool pass_trigger = (
          getVariableValue("H_Photon22") > 0 || 
          getVariableValue("H_Photon25") > 0 || 
          getVariableValue("H_Photon30") > 0 || 
          getVariableValue("H_Photon33") > 0 || 
          getVariableValue("H_Photon36") > 0 || 
          getVariableValue("H_Photon50") > 0 || 
          getVariableValue("H_Photon75") > 0 || 
          getVariableValue("H_Photon90") > 0 || 
          getVariableValue("H_Photon120"     ) > 0 || 
          getVariableValue("H_Photon150"     ) > 0 || 
          getVariableValue("H_Photon175"     ) > 0 || 
          getVariableValue("H_Photon200"     ) > 0 );

      fillVariableWithValue ("PassTrigger", pass_trigger ? true : false );

    }

    else if ( reducedSkimType == 1 || reducedSkimType == 2 || reducedSkimType == 3 || reducedSkimType == 4 ) { 

      // search by prefix
      if(triggerExists("HLT_Ele27_WPTight_Gsf"))
        fillTriggerVariable( "HLT_Ele27_WPTight_Gsf" , "H_Ele27_WPTight" );
      if(triggerExists("HLT_Ele32_WPTight_Gsf"))
        fillTriggerVariable( "HLT_Ele32_WPTight_Gsf" , "H_Ele32_WPTight" );
      if(triggerExists("HLT_Ele35_WPTight_Gsf"))
        fillTriggerVariable( "HLT_Ele35_WPTight_Gsf" , "H_Ele35_WPTight" );
      // check that we have at least one WPTight trigger
      if(!triggerExists("HLT_Ele27_WPTight_Gsf") && !triggerExists("HLT_Ele32_WPTight_Gsf") && !triggerExists("HLT_Ele35_WPTight_Gsf"))
        exit(-5);
      // Ele115 is absent from first 5/fb of 2017
      if(triggerExists("HLT_Ele115_CaloIdVT_GsfTrkIdT"))
        fillTriggerVariable( "HLT_Ele115_CaloIdVT_GsfTrkIdT" , "H_Ele115_CIdVT_GsfIdT");
      if(triggerExists("HLT_Photon175"))
        fillTriggerVariable( "HLT_Photon175" , "H_Photon175" );
      if(triggerExists("HLT_Photon200"))
        fillTriggerVariable( "HLT_Photon200" , "H_Photon200" );
      // check that we have at least one photon trigger
      if(!triggerExists("HLT_Photon175") && !triggerExists("HLT_Photon200"))
        exit(-5);
      // other triggers
      if(triggerExists("HLT_Ele105_CaloIdVT_GsfTrkIdT"))
        fillTriggerVariable( "HLT_Ele105_CaloIdVT_GsfTrkIdT" , "H_Ele105_CIdVT_GsfIdT");
      if(triggerExists("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL"))
        fillTriggerVariable( "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", "H_DoubleEle33_CIdL_GsfIdVL" ); 

      //bool pass_lowPtEle = (triggerExists("HLT_Ele27_WPTight_Gsf") && getVariableValue("H_Ele27_WPTight") > 0) ||
      //  (triggerExists("HLT_Ele32_WPTight_Gsf") && getVariableValue("H_Ele32_WPTight") > 0) ||
      //  (triggerExists("HLT_Ele35_WPTight_Gsf") && getVariableValue("H_Ele35_WPTight") > 0);
      //bool pass_photon = (triggerExists("HLT_Photon175") && getVariableValue("H_Photon175") > 0) ||
      //  (triggerExists("HLT_Photon200") && getVariableValue("H_Photon200") > 0);
      // NB: need to change the EventTriggerScaleFactorErr above if the trigger selection changes below
      //bool pass_trigger = (
      //    pass_lowPtEle || 
      //    //getVariableValue("H_Ele115_CIdVT_GsfIdT") > 0 ||
      //    pass_photon);
      //fillVariableWithValue ("PassTrigger", pass_trigger ? true : false );

      // NB: need to change the EventTriggerScaleFactorErr above if the trigger selection changes below
      // disambiguate events with multiple triggers at the analysis stage
      bool passHLT = false;
      if(analysisYearInt==2016) {
        if (getVariableValue("H_Photon175") == 1 ||
            getVariableValue("H_Ele27_WPTight") == 1 )
          passHLT = true;
      }
      else if(analysisYearInt==2017) {
        if (getVariableValue("H_Photon200") == 1 ||
            getVariableValue("H_Ele35_WPTight") == 1 )
          passHLT = true;
      }
      else if(analysisYearInt==2018) {
        if (getVariableValue("H_Photon200") == 1 ||
            getVariableValue("H_Ele32_WPTight") == 1 )
          passHLT = true;
      }
      fillVariableWithValue ("PassTrigger", passHLT ? true : false );
    }

    //-----------------------------------------------------------------
    // Evaluate the cuts
    //-----------------------------------------------------------------    

    evaluateCuts();

  } // event loop
  std::cout << "analysisClass::Loop(): ends " << std::endl;
}
