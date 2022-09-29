#define analysisClass_cxx
#include "analysisClass.h"
#include <typeinfo>
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
  // Analysis year
  //--------------------------------------------------------------------------
  int analysisYear = getPreCutValue1("AnalysisYear");

  //--------------------------------------------------------------------------
  // Cuts for physics objects selection
  //--------------------------------------------------------------------------
  // electron cuts
  double ele_PtCut    	         = getPreCutValue1("ele_PtCut");
  double ele_hltMatch_DeltaRCut  = getPreCutValue1("ele_hltMatch_DeltaRCut");

  //--------------------------------------------------------------------------
  // IDs
  //--------------------------------------------------------------------------
  std::string electronIDType     = getPreCutString1("electronIDType");
  if(electronIDType != "HEEP" && electronIDType != "EGMLoose") {
    STDOUT("electronIDType=" << electronIDType << " is unknown! Please implement it in the analysisClass code. Exiting.");
    exit(-5);
  }

  //--------------------------------------------------------------------------
  // Tell the user how many entries we'll look at
  //--------------------------------------------------------------------------
  Long64_t nentries = GetTreeEntries();
  std::cout << "analysisClass::analysisClass(): nentries = " << nentries << std::endl;
  
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
    // Define initial, inclusive collections for electrons
    //-----------------------------------------------------------------
    CollectionPtr c_ele_all   ( new Collection(readerTools_, readerTools_->ReadValueBranch<UInt_t>("nElectron")));
    //c_ele_all->examine<Electron>("c_ele_all = All reco electrons");
    //Electron ele1 = c_ele_all -> GetConstituent<Electron>(0);
    //for(unsigned int i=0; i<10; ++i) { 
    //  std::cout << "cut = " << i << " idLevel = " << ele1.GetNbitFromBitMap(i, 3) << std::endl;
    //}

    CollectionPtr c_ele_final;
    CollectionPtr c_ele_final_ptCut;
    CollectionPtr c_ele_vLoose_ptCut;
    CollectionPtr c_ele_loose;
    CollectionPtr c_ele_vLoose;

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
    FillUserTH1D("nEleNTuple",c_ele_all->GetSize());
    FillUserTH1D("nEleNrsk",c_ele_final->GetSize());
    FillUserTH2D("nEleNTupleVsNeleRsk",c_ele_final->GetSize(),c_ele_all->GetSize());

    //-----------------------------------------------------------------
    // Get ready to fill variables 
    //-----------------------------------------------------------------

    resetCuts();

    //-----------------------------------------------------------------
    // Pass JSON
    //-----------------------------------------------------------------
    fillVariableWithValue("PassJSON", passJSON(readerTools_->ReadValueBranch<UInt_t>("run"), readerTools_->ReadValueBranch<UInt_t>("luminosityBlock"), isData()));

    //-----------------------------------------------------------------
    // How many ID'd objects are there?
    //-----------------------------------------------------------------
    int n_ele_store          = c_ele_final                   -> GetSize();
    int n_ele_ptCut          = c_ele_final_ptCut             -> GetSize();
    int n_ele_vloose_ptCut   = c_ele_vLoose_ptCut            -> GetSize();
    int n_ele_vloose   = c_ele_vLoose            -> GetSize();

    //-----------------------------------------------------------------
    // Fill variables
    //-----------------------------------------------------------------
    fillVariableWithValue ("nVLooseEle"  , n_ele_vloose );
    fillVariableWithValue ("nVLooseEle_ptCut"  , n_ele_vloose_ptCut );
    fillVariableWithValue ("nLooseEle"  , n_ele_store );
    if ( n_ele_store >= 1 ) { 
      Electron ele = c_ele_final -> GetConstituent<Electron>(0);
      std::string prefix = "Ele1";
      fillVariableWithValue( prefix+"_Pt"            , ele.PtUncorr()           );
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

    bool pass_any_photon_trigger = (
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

    // just search by prefix
    //if(triggerExists("HLT_Ele27_WPLoose_Gsf"))
    //  fillTriggerVariable( "HLT_Ele27_WPLoose_Gsf" , "H_Ele27_WPLoose" );
    if(triggerExists("HLT_Ele27_WPTight_Gsf"))
      fillTriggerVariable( "HLT_Ele27_WPTight_Gsf" , "H_Ele27_WPTight" );
    else
      fillVariableWithValue( "H_Ele27_WPTight" , -1.0 );
    if(triggerExists("HLT_Ele32_WPTight_Gsf"))
      fillTriggerVariable( "HLT_Ele32_WPTight_Gsf" , "H_Ele32_WPTight" );
    else
      fillVariableWithValue( "H_Ele32_WPTight" , -1.0 );
    if(triggerExists("HLT_Ele35_WPTight_Gsf"))
      fillTriggerVariable( "HLT_Ele35_WPTight_Gsf" , "H_Ele35_WPTight" );
    else
      fillVariableWithValue( "H_Ele35_WPTight" , -1.0 );
    // check that we have at least one WPTight trigger
    if(!triggerExists("HLT_Ele27_WPTight_Gsf") && !triggerExists("HLT_Ele32_WPTight_Gsf") && !triggerExists("HLT_Ele35_WPTight_Gsf")) {
      STDOUT("Could not find any Ele WPTight trigger. Exiting.");
      exit(-5);
    }
    // Ele115 is absent from first 5/fb of 2017
    if(triggerExists("HLT_Ele115_CaloIdVT_GsfTrkIdT"))
      fillTriggerVariable( "HLT_Ele115_CaloIdVT_GsfTrkIdT" , "H_Ele115_CIdVT_GsfIdT");
    else
      fillVariableWithValue( "H_Ele115_CIdVT_GsfIdT" , -1.0 );
    if(triggerExists("HLT_Photon175"))
      fillTriggerVariable( "HLT_Photon175" , "H_Photon175" );
    else
      fillVariableWithValue( "H_Photon175" , -1.0 );
    if(triggerExists("HLT_Photon200"))
      fillTriggerVariable( "HLT_Photon200" , "H_Photon200" );
    else
      fillVariableWithValue( "H_Photon200" , -1.0 );
    // check that we have at least one photon trigger
    if(!triggerExists("HLT_Photon175") && !triggerExists("HLT_Photon200")) {
      STDOUT("Could not find a Photon175 or Photon200 trigger. Exiting.");
      exit(-5);
    }
    // other triggers
    //fillTriggerVariable( "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50", "H_Ele45_PFJet200_PFJet50");
    if(triggerExists("HLT_Ele105_CaloIdVT_GsfTrkIdT"))
      fillTriggerVariable( "HLT_Ele105_CaloIdVT_GsfTrkIdT" , "H_Ele105_CIdVT_GsfIdT");
    else
      fillVariableWithValue( "H_Ele105_CIdVT_GsfIdT" , -1.0 );
    if(triggerExists("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL"))
      fillTriggerVariable( "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", "H_DoubleEle33_CIdL_GsfIdVL" ); 
    else
      fillVariableWithValue( "H_DoubleEle33_CIdL_GsfIdVL", -1.0 );
    //fillTriggerVariable( "HLT_Mu45_eta2p1"  , "H_Mu45_eta2p1" );

    bool pass_lowPtEle = getVariableValue("H_Ele27_WPTight") > 0 ||
      getVariableValue("H_Ele32_WPTight") > 0 ||
      getVariableValue("H_Ele35_WPTight") > 0;
    bool pass_photon = getVariableValue("H_Photon175") > 0 || getVariableValue("H_Photon200") > 0;

    bool pass_trigger =
        pass_lowPtEle || 
        getVariableValue("H_Ele115_CIdVT_GsfIdT") > 0 ||
        pass_photon ||
        pass_any_photon_trigger;
    fillVariableWithValue ("PassTrigger", pass_trigger);

    //-----------------------------------------------------------------
    // Evaluate the cuts
    //-----------------------------------------------------------------    
    evaluateCuts();

  } // event loop
  std::cout << "analysisClass::Loop(): ends " << std::endl;
}
