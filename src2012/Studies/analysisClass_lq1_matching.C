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
    // Print progress
    //-----------------------------------------------------------------
    
    if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << "/" << nentries << std::endl;   
    
    //-----------------------------------------------------------------
    // Define initial, inclusive collections for all physics objects
    //-----------------------------------------------------------------

    CollectionPtr c_gen_all   ( new Collection(*this, GenParticlePt -> size()));
    CollectionPtr c_genJet_all( new Collection(*this, GenJetPt      -> size()));
    CollectionPtr c_ele_all   ( new Collection(*this, ElectronPt    -> size()));
    CollectionPtr c_muon_all  ( new Collection(*this, MuonPt        -> size()));
    CollectionPtr c_pfjet_all ( new Collection(*this, PFJetPt       -> size()));

    //-----------------------------------------------------------------
    // Store GEN electrons in c_genEle_final
    // You can find more GenParticle ID types, here: src/GenParticleIDs.C
    //-----------------------------------------------------------------

    CollectionPtr c_genEle_final = c_gen_all -> SkimByID<GenParticle>(GEN_ELE_HARD_SCATTER);

    //-----------------------------------------------------------------
    // Store GEN jets in c_genJet_final
    // There are no GenJet ID types coded for now
    //-----------------------------------------------------------------
    
    CollectionPtr c_genJet_final = c_genJet_all;

    //-----------------------------------------------------------------
    // Get HEEP RECO electrons that pass a pt cut: c_ele_final
    // You can find more Electron ID types, here: src/ElectronIDs.C
    //-----------------------------------------------------------------

    CollectionPtr c_ele_HEEP        = c_ele_all -> SkimByID <Electron> ( HEEP );
    CollectionPtr c_ele_final       = c_ele_HEEP -> SkimByMinPt<Electron>( ele_PtCut  );
    
    //-----------------------------------------------------------------
    // Get good muons that pass a pt cut: c_muon_final
    // You can find more Muon ID types, here: src/MuonIDs.C
    //-----------------------------------------------------------------

    CollectionPtr c_muon_eta2p1      = c_muon_all       -> SkimByEtaRange<Muon> ( -muon_EtaCut, muon_EtaCut );
    CollectionPtr c_muon_eta2p1_ID   = c_muon_eta2p1    -> SkimByID      <Muon> ( MUON_TIGHT_PFISO04 );
    CollectionPtr c_muon_final       = c_muon_eta2p1_ID -> SkimByMinPt   <Muon> ( muon_PtCut );
    
    //-----------------------------------------------------------------
    // Get PF jets that:
    // - Pass an ID
    // - Don't overlap with muons
    // - Don't overlap with electrons
    // - Pass a pt cut 
    // Store in c_pfjet_final 
    //-----------------------------------------------------------------

    CollectionPtr c_pfjet_central                     = c_pfjet_all                          -> SkimByEtaRange   <PFJet>          ( -jet_EtaCut, jet_EtaCut   );
    CollectionPtr c_pfjet_central_ID                  = c_pfjet_central                      -> SkimByID         <PFJet>          ( PFJET_LOOSE );    
    CollectionPtr c_pfjet_central_ID_noMuonOverlap    = c_pfjet_central_ID                   -> SkimByVetoDRMatch<PFJet, Muon>    ( c_muon_final , jet_ele_DeltaRCut  );
    CollectionPtr c_pfjet_central_ID_noLeptonOverlap  = c_pfjet_central_ID_noMuonOverlap     -> SkimByVetoDRMatch<PFJet, Electron>( c_ele_final  , jet_muon_DeltaRCut );
    CollectionPtr c_pfjet_final                       = c_pfjet_central_ID_noLeptonOverlap   -> SkimByMinPt      <PFJet>          ( jet_PtCut );
    
    //-----------------------------------------------------------------
    // Get ready to fill variables 
    //-----------------------------------------------------------------

    resetCuts();

    //-----------------------------------------------------------------
    // Fill the variables (blank for now)
    //-----------------------------------------------------------------

    //-----------------------------------------------------------------
    // Evaluate the cuts
    //-----------------------------------------------------------------    

    evaluateCuts();

    //-----------------------------------------------------------------
    // How many ID'd objects are there?
    //-----------------------------------------------------------------
    
    int n_ele    = c_ele_final    -> GetSize();
    int n_jet    = c_pfjet_final  -> GetSize();
    int n_genEle = c_genEle_final -> GetSize();
    int n_genJet = c_genJet_final -> GetSize();

    //-----------------------------------------------------------------
    // Example of how to loop over RECO jets
    //-----------------------------------------------------------------

    for ( int i_jet = 0; i_jet < n_jet ; ++i_jet ){
      PFJet jet = c_pfjet_final -> GetConstituent<PFJet>(i_jet);
      double jet_pt  = jet.Pt ();
      double jet_eta = jet.Eta();
      double jet_phi = jet.Phi();
    }
    
    //-----------------------------------------------------------------
    // Example of how to loop over GEN jets
    //-----------------------------------------------------------------

    for ( int i_genJet = 0; i_genJet < n_genJet ; ++i_genJet ){
      GenJet genJet  = c_genJet_final -> GetConstituent<GenJet>(i_genJet);
      double genJet_pt  = genJet.Pt ();
      double genJet_eta = genJet.Eta();
      double genJet_phi = genJet.Phi();
    }
    
    //-----------------------------------------------------------------
    // Example of how to loop over RECO electrons
    //-----------------------------------------------------------------

    for ( int i_ele = 0; i_ele < n_ele ; ++i_ele ){
      Electron ele      = c_ele_final -> GetConstituent<Electron>(i_ele);
      double ele_pt     = ele.Pt ();
      double ele_eta    = ele.Eta();
      double ele_phi    = ele.Phi();
      double ele_charge = ele.Charge();
    }
    
    //-----------------------------------------------------------------
    // Example of how to loop over GEN electrons
    //-----------------------------------------------------------------

    for ( int i_genEle = 0; i_genEle < n_genEle ; ++i_genEle ){
      GenParticle genEle = c_genEle_final -> GetConstituent<GenParticle>(i_genEle);
      double genEle_pt     = genEle.Pt ();
      double genEle_eta    = genEle.Eta();
      double genEle_phi    = genEle.Phi();
      int    genEle_pdg    = genEle.PdgId();
      int    genEle_charge = genEle_pdg < 0 ? 1 : -1;
    }
    
    //-----------------------------------------------------------------
    // Example of how to match each RECO jet to a GEN jet
    // (1) Loop over the RECO jets
    // (2) For each RECO jet, find the GEN jet that is closest to it
    //-----------------------------------------------------------------
    
    for ( int i_jet = 0; i_jet < n_jet ; ++i_jet ){
      PFJet jet = c_pfjet_final -> GetConstituent<PFJet>(i_jet);
      if ( n_genJet > 0 ) {
	GenJet closest_genJet = c_genJet_final -> GetClosestInDR<GenJet, PFJet> ( jet );
	double delta_r = jet.DeltaR ( & closest_genJet );
      
	// if you like, you can print out information about the matching
	
	if ( false ) { 
	  std::cout << "RECO jet #"        << i_jet + 1 << "    : " << jet << std::endl;
	  std::cout << "Matched GenJet: "                           << closest_genJet << ", DR = " << delta_r << std::endl;
	}
      }
    }

    //-----------------------------------------------------------------
    // Example of how to match each GEN jet to a RECO jet
    // (1) Loop over the GEN jets
    // (2) For each GEN jet, find the RECO jet that is closest to it
    // 
    // NOTE: You probably don't want to do this.
    // GenJets include lots of lower-pt jets from random gluons.
    // Most GEN jets will not have a corresponding RECO jet.
    //-----------------------------------------------------------------    

    for ( int i_genJet = 0; i_genJet < n_genJet ; ++i_genJet ){
      GenJet genJet = c_genJet_final -> GetConstituent<GenJet>(i_genJet);
      if ( n_jet > 0 ) {
	PFJet closest_jet = c_pfjet_final -> GetClosestInDR<PFJet,GenJet> ( genJet );
	double delta_r = genJet.DeltaR ( & closest_jet );
      }
    }

    //-----------------------------------------------------------------
    // Example of how to match each RECO electron to a GEN electron
    // (1) Loop over the RECO electrons
    // (2) For each RECO electron, find the GEN electron that is closest to it
    //-----------------------------------------------------------------
    
    for ( int i_ele = 0; i_ele < n_ele ; ++i_ele ){
      Electron ele = c_ele_final -> GetConstituent<Electron>(i_ele);
      if ( n_genEle > 0 ){
	GenParticle closest_genEle = c_genEle_final -> GetClosestInDR<GenParticle,Electron>( ele );
	double delta_r = ele.DeltaR ( & closest_genEle );
      }
    }
    
    //-----------------------------------------------------------------
    // Example of how to match each GEN electron to a RECO electron
    // (1) Loop over the GEN electrons
    // (2) For each GEN electron, find the RECO electron that is closest to it
    //-----------------------------------------------------------------
    
    for ( int i_genEle = 0; i_genEle < n_genEle ; ++i_genEle ){
      GenParticle genEle = c_genEle_final -> GetConstituent<GenParticle>(i_genEle);
      if ( n_ele > 0 ){
	Electron closest_ele = c_ele_final -> GetClosestInDR<Electron,GenParticle>( genEle );
	double delta_r = genEle.DeltaR ( & closest_ele );
      }
    }

    //-----------------------------------------------------------------
    // Example of how to print out information about the objects, if desired 
    //-----------------------------------------------------------------

    if ( false ) { 

      // Information about the run, event, ls
      std::cout << "Examining: " << run << ":" << ls << ":" << event << std::endl;   
	    
      // Information about the electrons

      c_genEle_final -> examine<GenParticle> ("GEN electrons" );
      c_ele_final    -> examine<Electron>    ("RECO electrons");
      
      // Information about the jets
      
      c_genJet_final -> examine<GenJet>      ("GEN jets");
      c_pfjet_final  -> examine<PFJet>       ("RECO jets");
    }
    
  } 
  
  std::cout << "analysisClass::Loop() ends" <<std::endl;   
  
}
