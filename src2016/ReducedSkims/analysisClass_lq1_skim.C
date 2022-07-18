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

#include "correction.h"

//--------------------------------------------------------------------------
// Function for trigger matching 
//--------------------------------------------------------------------------

template < class Object1, class Object2 > 
double triggerMatchPt ( const CollectionPtr & collection, Object2 & target_object, double delta_r_cut, bool verbose=false ){
  double matched_pt = -999.0;
  if ( collection ) { 
    int size = collection -> GetSize();
    if ( size > 0 ){ 
      if(verbose) {
        std::cout << "triggerMatchPt(): try to find closest object in DR to object: " << target_object << "from collection: " << std::endl;
        collection->examine<HLTriggerObject>("trigObjs");
      }
      Object1 matched_object = collection -> GetClosestInDR <Object1, Object2> ( target_object );
      double dr = matched_object.DeltaR ( & target_object );
      if(verbose) {
        std::cout << "found matched_object: " << matched_object << " with dR=" << dr << std::endl;
      }
      if ( dr < delta_r_cut ) { 
        matched_pt = matched_object.Pt();
        if(verbose) {
          std::cout << "dr=" << dr << " < delta_r_cut=" << delta_r_cut << ", so matched_pt set to matched_object pt: " << matched_object << std::endl;
        }
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
  //std::string graphName = getPreCutString1("TriggerEfficiencyGraphName");
  HistoReader triggerScaleFactorReader(trigSFFileName,"SF_TH2F_Barrel","SF_TH2F_EndCap",false,false);

  //--------------------------------------------------------------------------
  // reco scale factors
  //--------------------------------------------------------------------------
  std::string recoSFFileName = getPreCutString1("RecoSFFileName");
  std::unique_ptr<HistoReader> recoScaleFactorReader = std::unique_ptr<HistoReader>(new HistoReader(recoSFFileName,"EGamma_SF2D","EGamma_SF2D",true,false));

  //--------------------------------------------------------------------------
  // ID scale factors
  //--------------------------------------------------------------------------
  std::string egmIDSFFileName = getPreCutString1("EGMIDSFFileName");
  std::unique_ptr<HistoReader> idScaleFactorReader = std::unique_ptr<HistoReader>(new HistoReader(egmIDSFFileName,"EGamma_SF2D","EGamma_SF2D",false,false));
  //FIXME TODO move to correctionlib

  //--------------------------------------------------------------------------
  // correctionlib PU
  //--------------------------------------------------------------------------
  std::string puJSONFileName = getPreCutString1("PUJSONFileName");
  auto cset_puAll = correction::CorrectionSet::from_file(puJSONFileName);
  auto cset_pu = cset_puAll->at("Collisions16_UltraLegacy_goldenJSON");

  //--------------------------------------------------------------------------
  // correctionlib B-tagging
  //--------------------------------------------------------------------------
  std::string btagJSONFileName = getPreCutString1("BTVJSONFileName");
  auto cset_btagAll = correction::CorrectionSet::from_file(btagJSONFileName);
  auto cset_btagDeepJetIncl = cset_btagAll->at("deepJet_incl");
  auto cset_btagDeepJetComb = cset_btagAll->at("deepJet_comb");
  auto cset_btagDeepJetMuJets = cset_btagAll->at("deepJet_mujets");

  //--------------------------------------------------------------------------
  // correctionlib JME
  //--------------------------------------------------------------------------
  std::string jmeJSONFileName = getPreCutString1("JMEJSONFileName");
  auto cset_jmeAll = correction::CorrectionSet::from_file(jmeJSONFileName);
  std::string jerTag = getPreCutString1("JERTag");

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
  // Analysis year
  //--------------------------------------------------------------------------
  int analysisYear = getPreCutValue1("AnalysisYear");

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
  
  // QCD photon filters
  CollectionPtr c_hltPhoton_QCD_all;

  // muon filter
  CollectionPtr c_hltMuon_SingleMu_all;

  // Signal
  //CollectionPtr c_hltEle45_Signal_all;
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
  //TLorentzVector v_PFMETType1XYCor;
  

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
    //// test
    //double event = readerTools_->ReadValueBranch<ULong64_t>("event");
    //double ls = readerTools_->ReadValueBranch<UInt_t>("luminosityBlock");
    //double run = readerTools_->ReadValueBranch<UInt_t>("run");
    //if(event!=21273 && event!=21288) continue;
    //if(ls!=227) continue;
    //// run ls event
    ////std::cout << static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(ls) << " " << static_cast<unsigned int>(event) << std::endl;
    ////std::string current_file_name ( readerTools_->GetTree()->GetCurrentFile()->GetName());
    ////cout << "Found the event! in file:" << current_file_name << endl;
    //if(jentry==0)
    //  std::cout << "WARNING WARNING WARNING -- ONLY RUNNING OVER FIRST 100 ENTRIES!" << std::endl;
    //if(jentry > 100) continue;
    ////test

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

    // QCD photon triggers
    std::vector<int> typeIds {11, 22};
    c_hltPhoton_QCD_all = helper.GetFilterObjectsByType(typeIds);
    //c_hltPhoton_QCD_all->examine<HLTriggerObject>("c_hltPhoton_QCD_all");

    //else if ( reducedSkimType == 1 || reducedSkimType == 2 || reducedSkimType == 3 || reducedSkimType == 4){

    //  // SingleMu trigger: path was HLT_Mu40_eta2p1_v9+
    //  //c_hltMuon_SingleMu_all       = helper.GetL3FilterObjectsByPath("HLT_Mu45_eta2p1_v"); // will do prefix matching

    //  //// Ele+jets signal triggers
    //  //CollectionPtr trigger_l3objects_all = helper.GetL3FilterObjectsByPath("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v");
    //  //c_trigger_l3jets_all = trigger_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_JET);
    //  ////c_trigger_l3jets_all->examine<HLTriggerObject>("c_trigger_l3jets_all");
    //  //// which one passed the last filter?
    //  ////CollectionPtr trigger_lastObjects_all = helper.GetLastFilterObjectsByPath("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v");

    //  //// electrons seem to come as TRIGGER_PHOTON most of the time
    //  //c_hltEle45_Signal_all = trigger_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_PHOTON);
    //  //// if not, try TRIGGER_ELECTRON
    //  //if(c_hltEle45_Signal_all->GetSize() == 0)
    //  //  c_hltEle45_Signal_all = trigger_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_ELECTRON);
    //  //// Note: could also be TRIGGER_CLUSTER?

    //  // jets
    //  //c_hltPFJet200_Signal_all     =  trigger_lastObjects_all->SkimByID<HLTriggerObject>(TRIGGER_JET);
    //  // get rid of overlaps
    //  // XXX no, keep all L3 jets
    //  //c_hltPFJet50_Signal_all      =  c_trigger_l3jets_all -> SkimByVetoDRMatch<HLTriggerObject,HLTriggerObject>( c_hltPFJet200_Signal_all, 0.3 );

    //  //// DoubleEle signal trigger
    //  //CollectionPtr double_ele_l3objects_all    = helper.GetL3FilterObjectsByPath("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v");
    //  //c_hltDoubleEle_Signal_all    = double_ele_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_PHOTON);
    //  //// if not, try TRIGGER_ELECTRON
    //  //if(c_hltDoubleEle_Signal_all->GetSize() == 0)
    //  //  c_hltDoubleEle_Signal_all = double_ele_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_ELECTRON);

    //  //// Tag and probe trigger
    //  //CollectionPtr tagProbe_l3objects_all;
    //  //if(!isData())
    //  //  tagProbe_l3objects_all  = helper.GetL3FilterObjectsByPath("HLT_Ele27_WPLoose_Gsf_v");
    //  //else
    //  //  tagProbe_l3objects_all  = helper.GetL3FilterObjectsByPath("HLT_Ele27_WPLoose_Gsf_v");
    //  //c_hltEle27WP85Gsf_all = tagProbe_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_PHOTON);
    //  //// if not, try TRIGGER_ELECTRON
    //  //if(c_hltEle27WP85Gsf_all->GetSize() == 0)
    //  //  c_hltEle27WP85Gsf_all = tagProbe_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_ELECTRON);

    //}


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

    //-----------------------------------------------------------------
    // If this is MC, smear jets if requested
    // Don't do it for data
    //-----------------------------------------------------------------
    if ( !isData() && do_jer) do_jer = true;
    else do_jer = false; // never for data

    //-----------------------------------------------------------------
    // Energy scaling and resolution smearing here
    //-----------------------------------------------------------------
    if ( do_eer || do_jer || do_ees || do_jes ) { 

      // If  you're scaling/smearing PFJets, recall that only jets with pt > 10 GeV affect the PFMET
      // Also, only scale/smear the jets in our eta range (jets in the calorimeter crack are suspect)
      c_pfjet_all = c_pfjet_all -> SkimByEtaRange<PFJet>( -jet_EtaCut, jet_EtaCut );

      // Set the PFMET difference to zero

      v_delta_met.SetPtEtaPhiM(0.,0.,0.,0.);

      // Do the energy scale / energy resolution operations
      // dR for matching = Rcone/2
      //   see: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
      std::string jerPtRes = jerTag+"_MC_PtResolution_AK4PFchs";
      std::string jerSF = jerTag+"_MC_ScaleFactor_AK4PFchs";

      if ( do_eer ) c_ele_all      -> MatchAndSmearEnergy <Electron, GenParticle> ( c_genEle_final, 0.4/2.0, rootEngine, v_delta_met );
      if ( do_jer ) c_pfjet_all    -> MatchAndSmearEnergy <PFJet   , GenJet     > ( c_genJet_final, 0.4/2.0, rootEngine, v_delta_met, cset_jmeAll->at(jerPtRes).get(), cset_jmeAll->at(jerSF).get() );
      if ( do_ees ) c_ele_all      -> ScaleEnergy <Electron> ( electron_energy_scale_sign, v_delta_met );
      if ( do_jes ) c_pfjet_all    -> ScaleEnergy <PFJet   > ( pfjet_energy_scale_sign   , v_delta_met );

      // Propagate the results to the PFMET

      v_PFMETType1Cor   .SetPtEtaPhiM( readerTools_->ReadValueBranch<Float_t>("MET_pt"), 0., readerTools_->ReadValueBranch<Float_t>("MET_phi"), 0. );
      //v_PFMETType1XYCor.SetPtEtaPhiM( (*PFMETType1XYCor)[0] , 0., (*PFMETPhiType1XYCor)[0] , 0. );

      v_PFMETType1Cor    = v_PFMETType1Cor    + v_delta_met;
      //v_PFMETType1XYCor = v_PFMETType1XYCor + v_delta_met;

      //FIXME
      //(*PFMETType1Cor      )[0] = v_PFMETType1Cor   .Pt();
      //(*PFMETType1XYCor   )[0] = v_PFMETType1XYCor.Pt();
      //
      //FIXME
      //(*PFMETPhiType1Cor   )[0] = v_PFMETType1Cor   .Phi();
      //(*PFMETPhiType1XYCor)[0] = v_PFMETType1XYCor.Phi();
    }

    // new systematics handling
    std::vector<Electron> smearedEles;
    std::vector<Electron> scaledUpEles;
    std::vector<Electron> scaledDownEles;
    if(!isData()) {
      smearedEles = c_ele_all->MatchAndSmearEnergy <Electron, GenParticle> ( c_genEle_final, 0.4/2.0, rootEngine, v_delta_met );
      scaledUpEles = c_ele_all->ScaleEnergy <Electron> (1 , v_delta_met );
      scaledDownEles = c_ele_all->ScaleEnergy <Electron> (-1, v_delta_met );
    }

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
    if(analysisYear == 2018)
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
        c_ele_egammaLoose = c_ele_all -> SkimByID<Electron>(EGAMMA_BUILTIN_LOOSE);
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
    //fillVariableWithValue( "bunch"    , bunch      );
    fillVariableWithValue( "event"    , readerTools_->ReadValueBranch<ULong64_t>("event")      );
    fillVariableWithValue( "ls"       , readerTools_->ReadValueBranch<UInt_t>("luminosityBlock")         );
    //fillVariableWithValue( "orbit"    , orbit      );
    fillVariableWithValue( "run"      , readerTools_->ReadValueBranch<UInt_t>("run")        );
    //fillVariableWithValue( "ProcessID", ProcessID  );
    //fillVariableWithValue( "PtHat"    , PtHat      );
    float genWeight = -1.0;
    if(!isData()) {
      genWeight = readerTools_->ReadValueBranch<Float_t>("genWeight");
    }
    fillVariableWithValue( "Weight"   , genWeight   );
    //FIXME -- topPtWeights -- perhaps not needed since unused for 2016 analysis
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
        std::cout << "Got puWeight = " << puWeight << "; run: " << getVariableValue("run") << " ls: " << getVariableValue("ls") << " event: " << getVariableValue("event") << std::endl;
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
    // Pass JSON
    //-----------------------------------------------------------------
    fillVariableWithValue("PassJSON"                   , passJSON(getVariableValue("run"), getVariableValue("ls"), isData())                       );    

    //-----------------------------------------------------------------
    // Fill MET filter values
    // https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    //-----------------------------------------------------------------
    fillVariableWithValue("PassGlobalSuperTightHalo2016Filter" , int(readerTools_->ReadValueBranch<Bool_t>("Flag_globalSuperTightHalo2016Filter")          == 1));
    fillVariableWithValue("PassGoodVertices"              , int(readerTools_->ReadValueBranch<Bool_t>("Flag_goodVertices")                       == 1));
    fillVariableWithValue("PassHBHENoiseFilter"           , int(readerTools_->ReadValueBranch<Bool_t>("Flag_HBHENoiseFilter")                    == 1));
    fillVariableWithValue("PassHBHENoiseIsoFilter"        , int(readerTools_->ReadValueBranch<Bool_t>("Flag_HBHENoiseIsoFilter")                 == 1));
    fillVariableWithValue("PassBadEESupercrystalFilter"   , int(readerTools_->ReadValueBranch<Bool_t>("Flag_eeBadScFilter")                      == 1));
    fillVariableWithValue("PassEcalDeadCellTrigPrim"      , int(readerTools_->ReadValueBranch<Bool_t>("Flag_EcalDeadCellTriggerPrimitiveFilter") == 1));
    std::string branchName = "Flag_BadChargedCandidateFilter";
    //std::string branchType = std::string(readerTools_->GetTree()->GetBranch(branchName.c_str())->GetLeaf(branchName.c_str())->GetTypeName());
    ////std::cout << "Found branchType=" << branchType << std::endl;
    //if(branchType=="Bool_t") {
    //  fillVariableWithValue("PassChargedCandidateFilter"    , int(readerTools_->ReadValueBranch<Bool_t>(branchName)          == 1));
    //  fillVariableWithValue("PassBadPFMuonFilter"           , int(readerTools_->ReadValueBranch<Bool_t>("Flag_BadPFMuonFilter")                    == 1));
    //}
    //else {
    //  fillVariableWithValue("PassChargedCandidateFilter"    , int(readerTools_->ReadValueBranch<UChar_t>(branchName)          == 1));
    //  fillVariableWithValue("PassBadPFMuonFilter"           , int(readerTools_->ReadValueBranch<UChar_t>("Flag_BadPFMuonFilter")                    == 1));
    //}
    fillVariableWithValue("PassChargedCandidateFilter"    , int(readerTools_->ReadValueBranch<Bool_t>(branchName)          == 1));
    fillVariableWithValue("PassBadPFMuonFilter"           , int(readerTools_->ReadValueBranch<Bool_t>("Flag_BadPFMuonFilter")                    == 1));
    // for 2017 and 2018 only
    if(hasBranch("Flag_ecalBadCalibFilterV2"))
      fillVariableWithValue("PassEcalBadCalibV2Filter"    , int(readerTools_->ReadValueBranch<Bool_t>(branchName)          == 1));
    else
      fillVariableWithValue("PassEcalBadCalibV2Filter"    , 1);

    //-----------------------------------------------------------------
    // Fill MET values
    //-----------------------------------------------------------------

    fillVariableWithValue("PFMET_Type1_Pt"     , readerTools_->ReadValueBranch<Float_t>("MET_pt"));      
    fillVariableWithValue("PFMET_Type1_Phi"    , readerTools_->ReadValueBranch<Float_t>("MET_phi"));
    //fillVariableWithValue("PFMET_Type1XY_Pt"   , PFMETType1XYCor    -> at (0));      
    //fillVariableWithValue("PFMET_Type1XY_Phi"  , PFMETPhiType1XYCor -> at (0));

    if ( !isData() ) { 
      if ( reducedSkimType != 0 ){ 
        fillVariableWithValue("GenMET_Pt"		, readerTools_->ReadValueBranch<Float_t>("GenMET_pt"));
        fillVariableWithValue("GenMET_Phi"	, readerTools_->ReadValueBranch<Float_t>("GenMET_phi"));
        // add LHE variables if needed
        if(hasBranch("LHE_Vpt")) {
          fillVariableWithValue("LHE_Vpt"	  , readerTools_->ReadValueBranch<Float_t>("LHE_Vpt"));
          fillVariableWithValue("LHE_NpLO"	, readerTools_->ReadValueBranch<UChar_t>("LHE_NpLO"));
          fillVariableWithValue("LHE_NpNLO"	, readerTools_->ReadValueBranch<UChar_t>("LHE_NpNLO"));
          fillVariableWithValue("LHE_Njets"	, readerTools_->ReadValueBranch<UChar_t>("LHE_Njets"));
          fillVariableWithValue("LHE_Nglu"	, readerTools_->ReadValueBranch<UChar_t>("LHE_Nglu"));
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

    //fillVariableWithValue( "nPileUpInt_BXminus1", -1 );
    //fillVariableWithValue( "nPileUpInt_BX0"     , -1 );
    //fillVariableWithValue( "nPileUpInt_BXplus1" , -1 );
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

    int n_muonLoose          = c_muon_eta_IDLoose            -> GetSize();
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

    int n_filters = c_hltPhoton_QCD_all -> GetSize();

    for(int iEle = 0; iEle < getVariableValue("nEle_store"); ++iEle) {
      Electron ele = c_ele_final -> GetConstituent<Electron>(iEle);
      std::string prefix = "Ele"+std::to_string(iEle+1);
      //if(fabs(ele1.Eta()) >= 1.4442 && fabs(ele1.Eta()) <= 1.566)
      //  c_ele_final->examine<Electron>("c_ele_final");
      //double hltEle1Pt_signal          = triggerMatchPt<HLTriggerObject, Electron>(c_hltEle45_Signal_all    , ele1, ele_hltMatch_DeltaRCut);
      //double hltEle1Pt_doubleEleSignal = triggerMatchPt<HLTriggerObject, Electron>(c_hltDoubleEle_Signal_all, ele1, ele_hltMatch_DeltaRCut);
      //double hltEle1Pt_WP80            = triggerMatchPt<HLTriggerObject, Electron>(c_hltEle27WP85Gsf_all       , ele1, ele_hltMatch_DeltaRCut);
      double hltPhotonPt = -999.;
      if ( n_filters != 0 ) {
        //std::cout << "we have at least one photon trigger object. try to triggerMatchPt" << std::endl;
        hltPhotonPt = triggerMatchPt<HLTriggerObject, Electron>(c_hltPhoton_QCD_all, ele, ele_hltMatch_DeltaRCut);
      }
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
      if(analysisYear == 2018) {
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

      fillVariableWithValue( prefix+"_hltPhotonPt"  , hltPhotonPt );

      //fillVariableWithValue( prefix+"_hltEleSignalPt", hltEle1Pt_signal          );
      //fillVariableWithValue( prefix+"_hltDoubleElePt", hltEle1Pt_doubleEleSignal ); 
      //fillVariableWithValue( prefix+"_hltEleWP80Pt"  , hltEle1Pt_WP80            );
      if(!isData() && reducedSkimType != 0) {
        Electron ele1smeared = *find(smearedEles.begin(), smearedEles.end(), ele);
        Electron ele1scaledUp = *find(scaledUpEles.begin(), scaledUpEles.end(), ele);
        Electron ele1scaledDown = *find(scaledDownEles.begin(), scaledDownEles.end(), ele);
        fillVariableWithValue( prefix+"_Pt_EER"      , ele1smeared.Pt()          );
        fillVariableWithValue( prefix+"_Pt_EES_Up"    , ele1scaledUp.Pt()         );
        fillVariableWithValue( prefix+"_Pt_EES_Dn"  , ele1scaledDown.Pt()       );
      }
      if(!isData()) {
        float elePtUncorr = ele.PtUncorr();
        fillVariableWithValue( prefix+"_TrigSF" , triggerScaleFactorReader.LookupValue(ele.SCEta(),elePtUncorr) );
        fillVariableWithValue( prefix+"_TrigSF_Err" , triggerScaleFactorReader.LookupValueError(ele.SCEta(),elePtUncorr) );
        fillVariableWithValue( prefix+"_RecoSF" , recoScaleFactorReader->LookupValue(ele.SCEta(),elePtUncorr) );
        fillVariableWithValue( prefix+"_RecoSF_Err" , recoScaleFactorReader->LookupValueError(ele.SCEta(),elePtUncorr) );
        fillVariableWithValue( prefix+"_HEEPSF"                                      ,ElectronScaleFactorsRunII::LookupHeepSF(ele.SCEta(), analysisYear) );
        fillVariableWithValue( prefix+"_HEEPSF_Err"                                  ,ElectronScaleFactorsRunII::LookupHeepSFSyst(ele.SCEta(), elePtUncorr, analysisYear) );
        fillVariableWithValue( prefix+"_EGMLooseIDSF", idScaleFactorReader->LookupValue(ele.SCEta(),elePtUncorr) );
        fillVariableWithValue( prefix+"_EGMLooseIDSF_Err", idScaleFactorReader->LookupValueError(ele.SCEta(),elePtUncorr) );
      }
      else {
        fillVariableWithValue( prefix+"_TrigSF"          , 1.0);
        fillVariableWithValue( prefix+"_TrigSF_Err"      , 0.0);
        fillVariableWithValue( prefix+"_RecoSF"          , 1.0);
        fillVariableWithValue( prefix+"_RecoSF_Err"      , 0.0);
        fillVariableWithValue( prefix+"_HEEPSF"          , 1.0);
        fillVariableWithValue( prefix+"_HEEPSF_Err"      , 0.0);
        fillVariableWithValue( prefix+"_EGMLooseIDSF"    , 1.0);
        fillVariableWithValue( prefix+"_EGMLooseIDSF_Err", 0.0);
      }
    }

    // Jets
    if ( n_jet_highEta_store >= 1 ) { 
      PFJet jet1 = c_pfjet_highEta_final -> GetConstituent<PFJet>(0);
      fillVariableWithValue( "HighEtaJet1_Pt" , jet1.Pt () );
      fillVariableWithValue( "HighEtaJet1_Eta", jet1.Eta() );
      fillVariableWithValue( "HighEtaJet1_Phi", jet1.Phi() );
    }

    for(int iJet = 0; iJet < getVariableValue("nJet_store"); ++iJet) {
      PFJet jet = c_pfjet_final -> GetConstituent<PFJet>(iJet);
      std::string prefix = "Jet"+std::to_string(iJet+1);
      // leading HLT jet from 200 GeV collection
      double hltJet1Pt     = triggerMatchPt<HLTriggerObject, PFJet >( c_trigger_l3jets_all     , jet, jet_hltMatch_DeltaRCut);
      int flavor = isData() ? 0 : jet.HadronFlavor();
      float pt = jet.Pt();
      fillVariableWithValue( prefix+"_Pt"          , pt                         );
      fillVariableWithValue( prefix+"_PtOriginal"  , jet.PtOrignal()            );
      fillVariableWithValue( prefix+"_Eta"         , jet.Eta()                  );
      fillVariableWithValue( prefix+"_Phi"         , jet.Phi()                  );
      fillVariableWithValue( prefix+"_btagDeepJet" , jet.DeepJetBTag()          );
      fillVariableWithValue( prefix+"_hltJetPt"	   , hltJet1Pt                  );
      fillVariableWithValue( prefix+"_HadronFlavor", flavor                     );
      //fillVariableWithValue( prefix+"_qgl"	  , jet1.QuarkGluonLikelihood()     );
      float absEta = fabs(jet.Eta());
      // b-tagging scale factors and uncertainties
      if(!isData() && reducedSkimType != 0) {
        if(flavor == 0) {
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetIncl"  , cset_btagDeepJetIncl->evaluate({"central", "L", flavor, absEta, pt}) );
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
          // mujets
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetMuJets"  , cset_btagDeepJetMuJets->evaluate({"central", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetMuJets_Up"  , cset_btagDeepJetMuJets->evaluate({"up", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetMuJets_UpCorrelated"  , cset_btagDeepJetMuJets->evaluate({"up_correlated", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetMuJets_UpUncorrelated"  , cset_btagDeepJetMuJets->evaluate({"up_uncorrelated", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetMuJets_Dn"  , cset_btagDeepJetMuJets->evaluate({"down", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetMuJets_DnCorrelated"  , cset_btagDeepJetMuJets->evaluate({"down_correlated", "L", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFLooseDeepJetMuJets_DnUncorrelated"  , cset_btagDeepJetMuJets->evaluate({"down_uncorrelated", "L", flavor, absEta, pt}) );
          // medium
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetMuJets"  , cset_btagDeepJetMuJets->evaluate({"central", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetMuJets_Up"  , cset_btagDeepJetMuJets->evaluate({"up", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetMuJets_UpCorrelated"  , cset_btagDeepJetMuJets->evaluate({"up_correlated", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetMuJets_UpUncorrelated"  , cset_btagDeepJetMuJets->evaluate({"up_uncorrelated", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetMuJets_Dn"  , cset_btagDeepJetMuJets->evaluate({"down", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetMuJets_DnCorrelated"  , cset_btagDeepJetMuJets->evaluate({"down_correlated", "M", flavor, absEta, pt}) );
          fillVariableWithValue( prefix+"_btagSFMediumDeepJetMuJets_DnUncorrelated"  , cset_btagDeepJetMuJets->evaluate({"down_uncorrelated", "M", flavor, absEta, pt}) );
        }
        //XXX SIC JUNE 2022 FIXME
        //fillVariableWithValue( prefix+"_Pt_JESTotal_Up"  , jet.PtJESTotalUp()        );
        //fillVariableWithValue( prefix+"_Pt_JESTotal_Dn", jet.PtJESTotalDown()        );
        //fillVariableWithValue( prefix+"_Pt_JER_Up"  , jet.PtJERUp()                  );
        //fillVariableWithValue( prefix+"_Pt_JER_Dn", jet.PtJERDown()                  );
        if(c_pfJetMatchedLQ->Has<PFJet>(jet))
          fillVariableWithValue( prefix+"_LQMatched", 1);
      }
    }

    //-----------------------------------------------------------------
    // Fill variables that depend on more than one object
    // All skims need this
    //-----------------------------------------------------------------

    TLorentzVector t_ele1, t_ele2, t_jet1, t_jet2, t_jet3;
    TLorentzVector t_MET;
    TLorentzVector t_ele1Smeared, t_ele1ScaledUp, t_ele1ScaledDown;
    TLorentzVector t_ele2Smeared, t_ele2ScaledUp, t_ele2ScaledDown;
    TLorentzVector t_jet1JESTotalUp, t_jet1JESTotalDown, t_jet1JERUp, t_jet1JERDown;
    TLorentzVector t_jet2JESTotalUp, t_jet2JESTotalDown, t_jet2JERUp, t_jet2JERDown;

    t_MET.SetPtEtaPhiM( readerTools_->ReadValueBranch<Float_t>("MET_pt"), 0.0, readerTools_->ReadValueBranch<Float_t>("MET_phi"), 0.0 );

    if ( n_jet_store >= 1 ){

      PFJet jet1 = c_pfjet_final -> GetConstituent<PFJet>(0);
      t_jet1.SetPtEtaPhiM ( jet1.Pt(), jet1.Eta(), jet1.Phi(), 0.0 );
      if(!isData()) {
        //XXX SIC JUNE 2022 FIXME
        //t_jet1JESTotalUp.SetPtEtaPhiM ( jet1.PtJESTotalUp(), jet1.Eta(), jet1.Phi(), 0.0 );
        //t_jet1JESTotalDown.SetPtEtaPhiM ( jet1.PtJESTotalDown(), jet1.Eta(), jet1.Phi(), 0.0 );
        //t_jet1JERUp.SetPtEtaPhiM ( jet1.PtJERUp(), jet1.Eta(), jet1.Phi(), 0.0 );
        //t_jet1JERDown.SetPtEtaPhiM ( jet1.PtJERDown(), jet1.Eta(), jet1.Phi(), 0.0 );
      }

      fillVariableWithValue ("mDPhi_METJet1", fabs( t_MET.DeltaPhi ( t_jet1 )));

      if ( n_jet_store >= 2 ){

        PFJet jet2 = c_pfjet_final -> GetConstituent<PFJet>(1);
        t_jet2.SetPtEtaPhiM ( jet2.Pt(), jet2.Eta(), jet2.Phi(), 0.0 );
        if(!isData()) {
        //XXX SIC JUNE 2022 FIXME
          //t_jet2JESTotalUp.SetPtEtaPhiM ( jet2.PtJESTotalUp(), jet2.Eta(), jet2.Phi(), 0.0 );
          //t_jet2JESTotalDown.SetPtEtaPhiM ( jet2.PtJESTotalDown(), jet2.Eta(), jet2.Phi(), 0.0 );
          //t_jet2JERUp.SetPtEtaPhiM ( jet2.PtJERUp(), jet2.Eta(), jet2.Phi(), 0.0 );
          //t_jet2JERDown.SetPtEtaPhiM ( jet2.PtJERDown(), jet2.Eta(), jet2.Phi(), 0.0 );
        }

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
        t_ele1.SetPtEtaPhiM ( ele1.PtUncorr(), ele1.Eta(), ele1.Phi(), 0.0 );
      if(!isData()) {
        Electron ele1smeared = *find(smearedEles.begin(), smearedEles.end(), ele1);
        Electron ele1scaledUp = *find(scaledUpEles.begin(), scaledUpEles.end(), ele1);
        Electron ele1scaledDown = *find(scaledDownEles.begin(), scaledDownEles.end(), ele1);
        t_ele1Smeared.SetPtEtaPhiM( ele1smeared.Pt(), ele1smeared.Eta(), ele1smeared.Phi(), 0.0);
        t_ele1ScaledUp.SetPtEtaPhiM( ele1scaledUp.Pt(), ele1scaledUp.Eta(), ele1scaledUp.Phi(), 0.0);
        t_ele1ScaledDown.SetPtEtaPhiM( ele1scaledDown.Pt(), ele1scaledDown.Eta(), ele1scaledDown.Phi(), 0.0);
      }
      if ( n_ele_store >= 2 ) {
        Electron ele2 = c_ele_final -> GetConstituent<Electron>(1);
        t_ele2.SetPtEtaPhiM ( ele2.Pt(), ele2.Eta(), ele2.Phi(), 0.0 );
        if ( reducedSkimType == 0 ) // for QCD skims, use the uncorrected Pt
          t_ele2.SetPtEtaPhiM ( ele2.PtUncorr(), ele2.Eta(), ele2.Phi(), 0.0 );
        if(!isData()) {
          Electron ele2smeared = *find(smearedEles.begin(), smearedEles.end(), ele2);
          Electron ele2scaledUp = *find(scaledUpEles.begin(), scaledUpEles.end(), ele2);
          Electron ele2scaledDown = *find(scaledDownEles.begin(), scaledDownEles.end(), ele2);
          t_ele2Smeared.SetPtEtaPhiM( ele2smeared.Pt(), ele2smeared.Eta(), ele2smeared.Phi(), 0.0);
          t_ele2ScaledUp.SetPtEtaPhiM( ele2scaledUp.Pt(), ele2scaledUp.Eta(), ele2scaledUp.Phi(), 0.0);
          t_ele2ScaledDown.SetPtEtaPhiM( ele2scaledDown.Pt(), ele2scaledDown.Eta(), ele2scaledDown.Phi(), 0.0);
        }

        TLorentzVector t_ele1ele2 = t_ele1 + t_ele2;
        TLorentzVector t_ele1ele2_eer = t_ele1Smeared + t_ele2Smeared;
        TLorentzVector t_ele1ele2_eesUp = t_ele1ScaledUp + t_ele2ScaledUp;
        TLorentzVector t_ele1ele2_eesDown = t_ele1ScaledDown + t_ele2ScaledDown;
        fillVariableWithValue ("M_e1e2" , t_ele1ele2.M ());
        fillVariableWithValue ("Pt_e1e2", t_ele1ele2.Pt());
        //std::cout << "Pt_e1e2 Pt=" << t_ele1ele2.Pt() << std::endl;
        // systs
        if(!isData() && reducedSkimType != 0) {
          fillVariableWithValue ("M_e1e2_EER" , t_ele1ele2_eer.M ());
          fillVariableWithValue ("M_e1e2_EES_Up" , t_ele1ele2_eesUp.M ());
          fillVariableWithValue ("M_e1e2_EES_Dn" , t_ele1ele2_eesDown.M ());
          fillVariableWithValue ("Pt_e1e2_EER" , t_ele1ele2_eer.Pt ());
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
          fillVariableWithValue("M_e1j1_EER"        , (t_ele1Smeared+t_jet1).M());
          fillVariableWithValue("M_e1j1_EES_Up"     , (t_ele1ScaledUp+t_jet1).M());
          fillVariableWithValue("M_e1j1_EES_Dn"     , (t_ele1ScaledDown+t_jet1).M());
          fillVariableWithValue("M_e1j1_JESTotal_Up", (t_ele1+t_jet1JESTotalUp).M());
          fillVariableWithValue("M_e1j1_JESTotal_Dn", (t_ele1+t_jet1JESTotalDown).M());
          fillVariableWithValue("M_e1j1_JER_Up"     , (t_ele1+t_jet1JERUp).M());
          fillVariableWithValue("M_e1j1_JER_Dn"     , (t_ele1+t_jet1JERDown).M());
        }

        if ( n_jet_store >= 2 ){ 

          TLorentzVector t_ele1jet2 = t_ele1 + t_jet2;
          fillVariableWithValue("DR_Ele1Jet2", t_ele1.DeltaR ( t_jet2 ));
          fillVariableWithValue("M_e1j2"     , t_ele1jet2.M());
          // systs
          if(!isData() && reducedSkimType != 0) {
            fillVariableWithValue("M_e1j2_EER"        , (t_ele1Smeared+t_jet2).M());
            fillVariableWithValue("M_e1j2_EES_Up"     , (t_ele1ScaledUp+t_jet2).M());
            fillVariableWithValue("M_e1j2_EES_Dn"     , (t_ele1ScaledDown+t_jet2).M());
            fillVariableWithValue("M_e1j2_JESTotal_Up", (t_ele1+t_jet2JESTotalUp).M());
            fillVariableWithValue("M_e1j2_JESTotal_Dn", (t_ele1+t_jet2JESTotalDown).M());
            fillVariableWithValue("M_e1j2_JER_Up"     , (t_ele1+t_jet2JERUp).M());
            fillVariableWithValue("M_e1j2_JER_Dn"     , (t_ele1+t_jet2JERDown).M());
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
          fillVariableWithValue("M_e2j1_EER"        , (t_ele2Smeared+t_jet1).M());
          fillVariableWithValue("M_e2j1_EES_Up"     , (t_ele2ScaledUp+t_jet1).M());
          fillVariableWithValue("M_e2j1_EES_Dn"     , (t_ele2ScaledDown+t_jet1).M());
          fillVariableWithValue("M_e2j1_JESTotal_Up", (t_ele2+t_jet1JESTotalUp).M());
          fillVariableWithValue("M_e2j1_JESTotal_Dn", (t_ele2+t_jet1JESTotalDown).M());
          fillVariableWithValue("M_e2j1_JER_Up"     , (t_ele2+t_jet1JERUp).M());
          fillVariableWithValue("M_e2j1_JER_Dn"     , (t_ele2+t_jet1JERDown).M());
        }

        if ( n_jet_store >= 2 ){ 

          TLorentzVector t_ele2jet2 = t_ele2 + t_jet2;
          fillVariableWithValue("DR_Ele2Jet2", t_ele2.DeltaR ( t_jet2 ));
          fillVariableWithValue("M_e2j2"     , t_ele2jet2.M());
          fillVariableWithValue("sT_eejj"    , t_ele1.Pt() + t_ele2.Pt() + t_jet1.Pt() + t_jet2.Pt());
          // systs
          if(!isData() && reducedSkimType != 0) {
            fillVariableWithValue("M_e2j2_EER"        , (t_ele2Smeared+t_jet2).M());
            fillVariableWithValue("M_e2j2_EES_Up"     , (t_ele2ScaledUp+t_jet2).M());
            fillVariableWithValue("M_e2j2_EES_Dn"     , (t_ele2ScaledDown+t_jet2).M());
            fillVariableWithValue("M_e2j2_JESTotal_Up", (t_ele2+t_jet2JESTotalUp).M());
            fillVariableWithValue("M_e2j2_JESTotal_Dn", (t_ele2+t_jet2JESTotalDown).M());
            fillVariableWithValue("M_e2j2_JER_Up"     , (t_ele2+t_jet2JERUp).M());
            fillVariableWithValue("M_e2j2_JER_Dn"     , (t_ele2+t_jet2JERDown).M());
            fillVariableWithValue("sT_eejj_EER"    , t_ele1Smeared.Pt() + t_ele2Smeared.Pt() + t_jet1.Pt() + t_jet2.Pt());
            fillVariableWithValue("sT_eejj_EES_Up"    , t_ele1ScaledUp.Pt() + t_ele2ScaledUp.Pt() + t_jet1.Pt() + t_jet2.Pt());
            fillVariableWithValue("sT_eejj_EES_Dn"    , t_ele1ScaledDown.Pt() + t_ele2ScaledDown.Pt() + t_jet1.Pt() + t_jet2.Pt());
            fillVariableWithValue("sT_eejj_JESTotal_Up"    , t_ele1.Pt() + t_ele2.Pt() + t_jet1JESTotalUp.Pt() + t_jet2JESTotalUp.Pt());
            fillVariableWithValue("sT_eejj_JESTotal_Dn"    , t_ele1.Pt() + t_ele2.Pt() + t_jet1JESTotalDown.Pt() + t_jet2JESTotalDown.Pt());
            fillVariableWithValue("sT_eejj_JER_Up"    , t_ele1.Pt() + t_ele2.Pt() + t_jet1JERUp.Pt() + t_jet2JERUp.Pt());
            fillVariableWithValue("sT_eejj_JER_Dn"    , t_ele1.Pt() + t_ele2.Pt() + t_jet1JERDown.Pt() + t_jet2JERDown.Pt());
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
        fillVariableWithValue( "H_Photon22", -1); 
      if(triggerExists("HLT_Photon25"))
        fillTriggerVariable ( "HLT_Photon25" , "H_Photon25"  );
      else
        fillVariableWithValue( "H_Photon25", -1); 
      if(triggerExists("HLT_Photon30"))
        fillTriggerVariable ( "HLT_Photon30" , "H_Photon30"  );
      else
        fillVariableWithValue( "H_Photon30", -1); 
      if(triggerExists("HLT_Photon33"))
        fillTriggerVariable ( "HLT_Photon33" , "H_Photon33"  );
      else
        fillVariableWithValue( "H_Photon33", -1); 
      if(triggerExists("HLT_Photon36"))
        fillTriggerVariable ( "HLT_Photon36" , "H_Photon36"  );
      else
        fillVariableWithValue( "H_Photon36", -1); 
      fillTriggerVariable ( "HLT_Photon50"   , "H_Photon50"  );
      fillTriggerVariable ( "HLT_Photon75"   , "H_Photon75"  );
      fillTriggerVariable ( "HLT_Photon90"   , "H_Photon90"  );
      fillTriggerVariable ( "HLT_Photon120"  , "H_Photon120" );
      if(triggerExists("HLT_Photon150"))
        fillTriggerVariable ( "HLT_Photon150", "H_Photon150" );
      else
        fillVariableWithValue( "H_Photon150", -1); 
      fillTriggerVariable ( "HLT_Photon175" , "H_Photon175" );
      if(triggerExists("HLT_Photon200"))
        fillTriggerVariable ( "HLT_Photon200" , "H_Photon200" );
      else
        fillVariableWithValue ( "H_Photon200" , -1 );

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

      fillVariableWithValue ("PassTrigger", pass_trigger ? 1 : 0 );

    }

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

      // just search by prefix
      //if(triggerExists("HLT_Ele27_WPLoose_Gsf"))
      //  fillTriggerVariable( "HLT_Ele27_WPLoose_Gsf" , "H_Ele27_WPLoose" );
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
      //fillTriggerVariable( "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50", "H_Ele45_PFJet200_PFJet50");
      if(triggerExists("HLT_Ele105_CaloIdVT_GsfTrkIdT"))
        fillTriggerVariable( "HLT_Ele105_CaloIdVT_GsfTrkIdT" , "H_Ele105_CIdVT_GsfIdT");
      if(triggerExists("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL"))
        fillTriggerVariable( "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", "H_DoubleEle33_CIdL_GsfIdVL" ); 
      //fillTriggerVariable( "HLT_Mu45_eta2p1"  , "H_Mu45_eta2p1" );

      //bool pass_lowPtEle = (triggerExists("HLT_Ele27_WPTight_Gsf") && getVariableValue("H_Ele27_WPTight") > 0) ||
      //  (triggerExists("HLT_Ele32_WPTight_Gsf") && getVariableValue("H_Ele32_WPTight") > 0) ||
      //  (triggerExists("HLT_Ele35_WPTight_Gsf") && getVariableValue("H_Ele35_WPTight") > 0);
      //bool pass_photon = (triggerExists("HLT_Photon175") && getVariableValue("H_Photon175") > 0) ||
      //  (triggerExists("HLT_Photon200") && getVariableValue("H_Photon200") > 0);
      //bool pass_trigger = (
      //    pass_lowPtEle || 
      //    getVariableValue("H_Ele115_CIdVT_GsfIdT") > 0 ||
      //    pass_photon);
      //fillVariableWithValue ("PassTrigger", pass_trigger ? 1 : 0 );
    }

    //-----------------------------------------------------------------
    // Evaluate the cuts
    //-----------------------------------------------------------------    

    evaluateCuts();

  } // event loop
  std::cout << "analysisClass::Loop(): ends " << std::endl;
}
