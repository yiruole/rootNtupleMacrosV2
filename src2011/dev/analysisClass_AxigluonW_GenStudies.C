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

  ////////////////////// User's code to get preCut values - BEGIN ///////////////

  int PDGID_AG =  getPreCutValue1("PDGID_axigluon");
  int PDGID_Wplus =  getPreCutValue1("PDGID_W");
  int PDGID_Wminus =  getPreCutValue2("PDGID_W");

  int PDGID_ELE =  getPreCutValue1("PDGID_ELE_MU_TAU");
  int PDGID_MU =  getPreCutValue2("PDGID_ELE_MU_TAU");
  int PDGID_TAU =  getPreCutValue3("PDGID_ELE_MU_TAU");
  int PDGID_NUELE =  getPreCutValue1("PDGID_NU_ELE_MU_TAU");
  int PDGID_NUMU =  getPreCutValue2("PDGID_NU_ELE_MU_TAU");
  int PDGID_NUTAU =  getPreCutValue3("PDGID_NU_ELE_MU_TAU");

  int PDGID_D =  getPreCutValue1("PDGID_D");
  int PDGID_DBAR =  getPreCutValue2("PDGID_D");
  int PDGID_U =  getPreCutValue1("PDGID_U");
  int PDGID_UBAR =  getPreCutValue2("PDGID_U");
  int PDGID_S =  getPreCutValue1("PDGID_S");
  int PDGID_SBAR =  getPreCutValue2("PDGID_S");
  int PDGID_C =  getPreCutValue1("PDGID_C");
  int PDGID_CBAR =  getPreCutValue2("PDGID_C");
  int PDGID_B =  getPreCutValue1("PDGID_B");
  int PDGID_BBAR =  getPreCutValue2("PDGID_B");
  int PDGID_T =  getPreCutValue1("PDGID_T");
  int PDGID_TBAR =  getPreCutValue2("PDGID_T");

  // Others
  int jetAlgorithm = getPreCutValue1("jetAlgorithm");
  int metAlgorithm = getPreCutValue1("metAlgorithm");
  int eleAlgorithm = getPreCutValue1("eleAlgorithm");
  int phoAlgorithm = getPreCutValue1("phoAlgorithm");
  double jet_ele_DeltaRcut = getPreCutValue1("jet_ele_DeltaRcut");
  double jet_muon_DeltaRcut = getPreCutValue1("jet_muon_DeltaRcut");
  double jet_pho_DeltaRcut = getPreCutValue1("jet_pho_DeltaRcut");

  // Electrons
  double elePtCut =  getPreCutValue1("elePtCut");
  double eleEta_bar = getPreCutValue1("eleEta_bar");
  double eleEta_end_min = getPreCutValue1("eleEta_end");
  double eleEta_end_max = getPreCutValue2("eleEta_end");
  double eleMissingHits = getPreCutValue1("eleMissingHits");
  double eleDist        = getPreCutValue1("eleDist");
  double eleDCotTheta   = getPreCutValue1("eleDCotTheta");
  double eleUseHasMatchConv = getPreCutValue1("eleUseHasMatchConv");
  double eleCombRelIso_bar = getPreCutValue1("eleCombRelIso");
  double eleCombRelIso_end = getPreCutValue2("eleCombRelIso");
  double eleSigmaIetaIeta_bar = getPreCutValue1("eleSigmaIetaIeta");
  double eleSigmaIetaIeta_end = getPreCutValue2("eleSigmaIetaIeta");
  double eleDeltaPhiTrkSC_bar = getPreCutValue1("eleDeltaPhiTrkSC");
  double eleDeltaPhiTrkSC_end = getPreCutValue2("eleDeltaPhiTrkSC");
  double eleDeltaEtaTrkSC_bar = getPreCutValue1("eleDeltaEtaTrkSC");
  double eleDeltaEtaTrkSC_end = getPreCutValue2("eleDeltaEtaTrkSC");
  double eleUseEcalDriven = getPreCutValue1("eleUseEcalDriven");

  double eleDeltaEtaTrkSCHeep_bar = getPreCutValue1("eleDeltaEtaTrkSCHeep");
  double eleDeltaEtaTrkSCHeep_end = getPreCutValue2("eleDeltaEtaTrkSCHeep");
  double eleDeltaPhiTrkSCHeep_bar = getPreCutValue1("eleDeltaPhiTrkSCHeep");
  double eleDeltaPhiTrkSCHeep_end = getPreCutValue2("eleDeltaPhiTrkSCHeep");
  double eleHoEHeep_bar = getPreCutValue1("eleHoEHeep");
  double eleHoEHeep_end = getPreCutValue2("eleHoEHeep");
  double eleE2x5OverE5x5Heep_bar = getPreCutValue1("eleE2x5OverE5x5Heep");
  //double eleE2x5OverE5x5Heep_end = getPreCutValue2("eleE2x5OverE5x5Heep");
  double eleE1x5OverE5x5Heep_bar = getPreCutValue1("eleE1x5OverE5x5Heep");
  //double eleE1x5OverE5x5Heep_end = getPreCutValue2("eleE1x5OverE5x5Heep");
  //double eleSigmaIetaIetaHeep_bar = getPreCutValue1("eleSigmaIetaIetaHeep");
  double eleSigmaIetaIetaHeep_end = getPreCutValue2("eleSigmaIetaIetaHeep");
  double eleEcalHcalIsoHeep_1_bar = getPreCutValue1("eleEcalHcalIsoHeep");
  double eleEcalHcalIsoHeep_2_bar = getPreCutValue2("eleEcalHcalIsoHeep");
  double eleEcalHcalIsoHeep_1_end = getPreCutValue3("eleEcalHcalIsoHeep");
  double eleEcalHcalIsoHeep_2_end = getPreCutValue4("eleEcalHcalIsoHeep");
  //double eleEcalHcalIso_PTthrHeep_bar = getPreCutValue1("eleEcalHcalIso_PTthrHeep");
  double eleEcalHcalIso_PTthrHeep_end = getPreCutValue2("eleEcalHcalIso_PTthrHeep");
  //double eleHcalIsoD2Heep_bar = getPreCutValue1("eleHcalIsoD2Heep");
  double eleHcalIsoD2Heep_end = getPreCutValue2("eleHcalIsoD2Heep");
  double eleTrkIsoHeep_bar = getPreCutValue1("eleTrkIsoHeep");
  double eleTrkIsoHeep_end = getPreCutValue2("eleTrkIsoHeep");
  double eleMissingHitsHeep = getPreCutValue1("eleMissingHitsHeep");
  double eleUseEcalDrivenHeep = getPreCutValue1("eleUseEcalDrivenHeep");

  //   int heepBitMask_EB  =  getPreCutValue1("heepBitMask_EBGapEE") ;
  //   int heepBitMask_GAP =  getPreCutValue2("heepBitMask_EBGapEE") ;
  //   int heepBitMask_EE  =  getPreCutValue3("heepBitMask_EBGapEE") ;

  // Muons 
  double muPtCut =  getPreCutValue1("muPtCut");
  double muEta_bar = getPreCutValue1("muEta_bar");
  double muEta_end_min = getPreCutValue1("muEta_end");
  double muEta_end_max = getPreCutValue2("muEta_end");
  double muNormChi2 = getPreCutValue1("muNormChi2");
  double muNValidHitsGlobal = getPreCutValue1("muNValidHitsGlobal");
  double muNMatchedStation = getPreCutValue1("muNMatchedStation");
  double muDxyTrack = getPreCutValue1("muDxyTrack");
  double muNValidHitsPixel = getPreCutValue1("muNValidHitsPixel");
  double muNValidHitsTrack = getPreCutValue1("muNValidHitsTrack");
  double muCombRelIso = getPreCutValue1("muCombRelIso");

  // Jets
  double jetPtCut = getPreCutValue1("jetPtCut");
  double jetEtaCut = getPreCutValue1("jetEtaCut");
  double jetTCHELCut = getPreCutValue1("jetTCHELCut");
  double jetIDloose = getPreCutValue1("jetIDloose");
  double jetIDtight = getPreCutValue1("jetIDtight");

  // Photons
  double phoPtCut = getPreCutValue1("phoPtCut");
  double phoEta_bar = getPreCutValue1("phoEta_bar");
  double phoEta_end_min = getPreCutValue1("phoEta_end");
  double phoEta_end_max = getPreCutValue2("phoEta_end");
  double phoEcalIso_bar_const = getPreCutValue1("phoEcalIso");
  double phoEcalIso_bar_ptCoeff = getPreCutValue2("phoEcalIso");
  double phoEcalIso_end_const = getPreCutValue3("phoEcalIso");
  double phoEcalIso_end_ptCoeff = getPreCutValue4("phoEcalIso");
  double phoHcalIso_bar_const = getPreCutValue1("phoHcalIso");
  double phoHcalIso_bar_ptCoeff = getPreCutValue2("phoHcalIso");
  double phoHcalIso_end_const = getPreCutValue3("phoHcalIso");
  double phoHcalIso_end_ptCoeff = getPreCutValue4("phoHcalIso");
  double phoTrkIso_bar_const = getPreCutValue1("phoTrkIso");
  double phoTrkIso_bar_ptCoeff = getPreCutValue2("phoTrkIso");
  double phoTrkIso_end_const = getPreCutValue3("phoTrkIso");
  double phoTrkIso_end_ptCoeff = getPreCutValue4("phoTrkIso");
  double phoHoE_bar = getPreCutValue1("phoHoE");
  double phoHoE_end = getPreCutValue1("phoHoE");
  double phoSigmaIetaIeta_bar = getPreCutValue1("phoSigmaIetaIeta");
  double phoSigmaIetaIeta_end = getPreCutValue2("phoSigmaIetaIeta");
  double phoUseMatchPromptEle = getPreCutValue1("phoUseMatchPromptEle");
  double phoUsePixelSeed = getPreCutValue1("phoUsePixelSeed");

  double vertexMinimumNDOF = getPreCutValue1("vertexMinimumNDOF");
  double vertexMaxAbsZ = getPreCutValue1("vertexMaxAbsZ");
  double vertexMaxd0 = getPreCutValue1("vertexMaxd0");

  double jet_PtCut_forMetScale = getPreCutValue1("jet_PtCut_forMetScale");
  double EleEnergyScale_EB=getPreCutValue1("EleEnergyScale_EB");
  double EleEnergyScale_EE=getPreCutValue1("EleEnergyScale_EE");
  double JetEnergyScale=getPreCutValue1("JetEnergyScale");
  int doJetOversmearing = getPreCutValue1("doJetOversmearing");
  double JetOversmearingSigma = getPreCutValue1("JetOversmearingSigma");

  ////////////////////// User's code to get preCut values - END /////////////////

  ////////////////////// User's code to book histos - BEGIN ///////////////////////

  //CreateUserTH1D("h1_LQGenEle_Pt", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
  //CreateUserTH2D("h2_pfMET_vs_neutrinoPt__sT", 100,0,1000, 100, 0, 1000);

  CreateUserTH1D("h1_Num_AG", 3, -0.5, 2.5);
  CreateUserTH1D("h1_Num_W", 3, -0.5, 2.5 );
  CreateUserTH1D("h1_Charge_W", 2, -1.01, 1.01);

  CreateUserTH1D("h1_Mass_AG", 1000, 0, 2000);
  CreateUserTH1D("h1_Mass_AG_fromqq", 1000, 0, 2000);
  CreateUserTH1D("h1_Mass_W", 100, 0, 200);
  CreateUserTH1D("h1_Mass_W_fromlnu", 100, 0, 200);
  CreateUserTH1D("h1_MT_W", 100, 0, 200);

  CreateUserTH1D("h1_Pt_AG", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_W", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_Wplus", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_Wminus", 100, 0, 1000);

  CreateUserTH1D("h1_Eta_AG", 100, -10, 10);
  CreateUserTH1D("h1_Eta_W", 100, -10, 10);
  CreateUserTH1D("h1_Eta_Wplus", 100, -10, 10);
  CreateUserTH1D("h1_Eta_Wminus", 100, -10, 10);

  CreateUserTH1D("h1_Pt_AG_q", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_AG_qbar", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_W_l", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_W_nu", 100, 0, 1000);

  CreateUserTH1D("h1_Eta_AG_q", 100, -10, 10);
  CreateUserTH1D("h1_Eta_AG_qbar", 100, -10, 10);
  CreateUserTH1D("h1_Eta_W_l", 100, -10, 10);
  CreateUserTH1D("h1_Eta_W_nu", 100, -10, 10);

  CreateUserTH1D("h1_Num_AG_Daughters", 20, 0, 20);		

  CreateUserTH1D("h1_GenMET", 100, 0, 1000);

  TH1F* h1_DecayProducts_AG = new TH1F("h1_DecayProducts_AG", "h1_DecayProducts_AG", 7, 0.5, 7.5);
  h1_DecayProducts_AG->Sumw2();
  h1_DecayProducts_AG->GetXaxis()->SetBinLabel(1,"ddbar");
  h1_DecayProducts_AG->GetXaxis()->SetBinLabel(2,"uubar");
  h1_DecayProducts_AG->GetXaxis()->SetBinLabel(3,"ssbar");
  h1_DecayProducts_AG->GetXaxis()->SetBinLabel(4,"ccbar");
  h1_DecayProducts_AG->GetXaxis()->SetBinLabel(5,"bbar");
  h1_DecayProducts_AG->GetXaxis()->SetBinLabel(6,"ttbar");
  h1_DecayProducts_AG->GetXaxis()->SetBinLabel(7,"total");

  TH1F* h1_DecayProducts_W = new TH1F("h1_DecayProducts_W", "h1_DecayProducts_W", 4, 0.5, 4.5);
  h1_DecayProducts_W->Sumw2();
  h1_DecayProducts_W->GetXaxis()->SetBinLabel(1,"e#nu");
  h1_DecayProducts_W->GetXaxis()->SetBinLabel(2,"#mu#nu");
  h1_DecayProducts_W->GetXaxis()->SetBinLabel(3,"#tau#nu");
  h1_DecayProducts_W->GetXaxis()->SetBinLabel(4,"total");


  // Random number generator
  TRandom3 *randomNumGen = new TRandom3;
  randomNumGen->SetSeed();

  CreateUserTH1D( "MuonAllPt"                  ,    200,  0,  2000  ); 
  CreateUserTH1D( "MuonAllEta"                 ,    100, -6, 6      ); 
  CreateUserTH1D( "EleAllPt"                   ,    200,  0, 2000   );
  CreateUserTH1D( "EleAllEta"                  ,    100, -6, 6      );
  CreateUserTH1D( "PhoAllPt"                   ,    200,  0, 2000   );
  CreateUserTH1D( "PhoAllEta"                  ,    100, -6, 6      );
  CreateUserTH1D( "JetAllPt"                   ,    200,  0, 2000   );
  CreateUserTH1D( "JetAllEta"                  ,    100, -6, 6      );
  CreateUserTH1D( "BJetAllPt"                  ,    200,  0, 2000   );
  CreateUserTH1D( "BJetAllEta"                 ,    100, -6, 6      );

  //Efficiency and acceptance

  CreateUserTH1D("h1_Pt_W_e", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_W_e_gencuts", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_W_e_recomatched", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_W_e_recomatched_IDISO", 100, 0, 1000);
  CreateUserTH1D("h1_Eta_W_e", 100, -10, 10);
  CreateUserTH1D("h1_Eta_W_e_gencuts", 100, -10, 10);
  CreateUserTH1D("h1_Eta_W_e_recomatched", 100, -10, 10);
  CreateUserTH1D("h1_Eta_W_e_recomatched_IDISO", 100, -10, 10);
  CreateUserTH1D("h1_minDR_elequark", 100, 0, 10);
  CreateUserTH1D("h1_minDR_elequark_gencuts", 100, 0, 10);
  CreateUserTH1D("h1_minDR_elequark_recomatched", 100, 0, 10);
  CreateUserTH1D("h1_minDR_elequark_recomatched_IDISO", 100, 0, 10);
  CreateUserTH1D("h1_nVtx_e", 50, 0, 50 );	 	    
  CreateUserTH1D("h1_nVtx_e_gencuts", 50 , 0, 50 );	 	    
  CreateUserTH1D("h1_nVtx_e_recomatched", 50 , 0, 50 );	 	    
  CreateUserTH1D("h1_nVtx_e_recomatched_IDISO", 50 , 0, 50 );	 	    

  CreateUserTH1D("h1_Pt_W_mu", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_W_mu_gencuts", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_W_mu_recomatched", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_W_mu_recomatched_IDISO", 100, 0, 1000);
  CreateUserTH1D("h1_Eta_W_mu", 100, -10, 10);
  CreateUserTH1D("h1_Eta_W_mu_gencuts", 100, -10, 10);
  CreateUserTH1D("h1_Eta_W_mu_recomatched", 100, -10, 10);
  CreateUserTH1D("h1_Eta_W_mu_recomatched_IDISO", 100, -10, 10);
  CreateUserTH1D("h1_minDR_muonquark", 100, 0, 10);
  CreateUserTH1D("h1_minDR_muonquark_gencuts", 100, 0, 10);
  CreateUserTH1D("h1_minDR_muonquark_recomatched", 100, 0, 10);
  CreateUserTH1D("h1_minDR_muonquark_recomatched_IDISO", 100, 0, 10);
  CreateUserTH1D("h1_nVtx_mu", 50, 0, 50 );	 	    
  CreateUserTH1D("h1_nVtx_mu_gencuts", 50 , 0, 50 );	 	    
  CreateUserTH1D("h1_nVtx_mu_recomatched", 50 , 0, 50 );	 	    
  CreateUserTH1D("h1_nVtx_mu_recomatched_IDISO", 50 , 0, 50 );	 	    

  ////////////////////// User's code to book histos - END ///////////////////////

  Long64_t nentries = fChain->GetEntriesFast();
  STDOUT("analysisClass::Loop(): nentries = " << nentries);

  ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
  ////// If the root version is updated and rootNtupleClass regenerated,     /////
  ////// these lines may need to be updated.                                 /////
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) { // Begin of loop over events
  // for (Long64_t jentry=0; jentry<1000;jentry++) { // Begin of loop over events
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry < 10 || jentry%1000 == 0) STDOUT("analysisClass::Loop(): jentry = " << jentry);
    // if (Cut(ientry) < 0) continue;

    ////////////////////// User's code to be done for every event - BEGIN ///////////////////////

    //## Define new jet collection
    std::auto_ptr<std::vector<double> >  JetPt  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetPtRaw  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetEnergy  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetEta  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetPhi  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >     JetPassLooseID  ( new std::vector<int>()  );
    std::auto_ptr<std::vector<int> >     JetPassTightID  ( new std::vector<int>()  );
    std::auto_ptr<std::vector<double> >  JetTCHE  ( new std::vector<double>()  );
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
            JetPassLooseID->push_back( PFJetPassLooseID->at(ijet) );
            JetPassTightID->push_back( PFJetPassTightID->at(ijet) );
            JetTCHE->push_back( PFJetTrackCountingHighEffBTag->at(ijet) );
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
            JetPassLooseID->push_back( CaloJetPassLooseID->at(ijet) );
            JetPassTightID->push_back( CaloJetPassTightID->at(ijet) );
            JetTCHE->push_back( CaloJetTrackCountingHighEffBTag->at(ijet) );
	  }//end loop over calo jets
      }//end if "calo jets"


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

    //## Define new electron Pt
    if(eleAlgorithm==2)
      {
	for(int iele=0; iele<ElectronPt->size(); iele++)
	  {
	    ElectronPt->at(iele) = ElectronPtHeep->at(iele);
	  }
      }

    //## EES and JES
    if( EleEnergyScale_EB != 1 || EleEnergyScale_EE != 1 )
      {
	for(int iele=0; iele<ElectronPt->size(); iele++)
	  {
	    if( fabs(ElectronEta->at(iele)) < eleEta_bar )
	      ElectronPt->at(iele) *= EleEnergyScale_EB;
	    if( fabs(ElectronEta->at(iele)) > eleEta_end_min && fabs(ElectronEta->at(iele)) < eleEta_end_max )
	      ElectronPt->at(iele) *= EleEnergyScale_EE;
	  }
      }
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

    ////////////////////// Reco Object Collections ///////////////////////

    //## Muons
    vector<int> v_idx_muon_PtCut;
    vector<int> v_idx_muon_PtCut_IDISO;

    // Loop over muons
    for(int imuon=0; imuon<MuonPt->size(); imuon++){
      
      // pT pre-cut on muon
      if ( MuonPt->at(imuon) < muPtCut ) continue;

      v_idx_muon_PtCut.push_back(imuon);

      // isolation
      double MuonCombRelIso  =  ( MuonTrackerIsoSumPT->at(imuon) 
				  + MuonEcalIso->at(imuon) 
				  + MuonHcalIso->at(imuon) 
				  - rhoIso*TMath::Pi()*0.3*0.3 
				  ) / MuonPt->at(imuon) ;

      // Pass Muon Selection
      if (   
	  //Is GlobalMuonPromptTight
	  MuonIsGlobal->at(imuon) 
	  && MuonGlobalChi2->at(imuon) < muNormChi2
	  && MuonGlobalTrkValidHits->at(imuon) > muNValidHitsGlobal
	  //Other requirements
	  && MuonStationMatches->at(imuon) > muNMatchedStation
	  && MuonPrimaryVertexDXY->at(imuon) < muDxyTrack
	  && MuonTrkPixelHitCount->at(imuon) > muNValidHitsPixel
	  && MuonTrkHitsTrackerOnly->at(imuon) > muNValidHitsTrack
	  //Isolation
	  && MuonCombRelIso < muCombRelIso
	  )
	{
	  v_idx_muon_PtCut_IDISO.push_back(imuon);

	  //histograms
	  FillUserTH1D( "MuonAllPt"                  ,    MuonPt->at(imuon)  ); 
	  FillUserTH1D( "MuonAllEta"                 ,    MuonEta->at(imuon) ); 
	}

    }// end loop over muons

    //## Electrons
    vector<int> v_idx_ele_PtCut;
    vector<int> v_idx_ele_PtCut_IDISO;

    //Loop over electrons
    for(int iele=0; iele<ElectronPt->size(); iele++)
      {
	
	if( eleAlgorithm == 1) //------> use EWK(WP80) SimpleCutBasedEleID2011
	  {

	    // pT pre-cut on ele
	    if( ElectronPt->at(iele) < elePtCut ) continue;

	    // ecal driven	    
	    if( eleUseEcalDriven && !ElectronHasEcalDrivenSeed->at(iele) )
	      continue;

	    v_idx_ele_PtCut.push_back(iele);

	    int passEleSel = 0;
	    int isBarrel = 0;
	    int isEndcap = 0;
	    int isPhotConv = 0;

	    // isolation
	    double ElectronCombRelIso_bar  =  ( ElectronTrkIsoDR03->at(iele) 
						+ max( 0., ElectronEcalIsoDR03->at(iele) - 1. ) 
						+ ElectronHcalIsoDR03FullCone->at(iele) 
						- rhoIso*TMath::Pi()*0.3*0.3 
						) / ElectronPt->at(iele) ;
	    
	    double ElectronCombRelIso_end  =  ( ElectronTrkIsoDR03->at(iele) 
						+ ElectronEcalIsoDR03->at(iele) 
						+ ElectronHcalIsoDR03FullCone->at(iele) 
						- rhoIso*TMath::Pi()*0.3*0.3 
						) / ElectronPt->at(iele) ;

	    // conversions
	    if(eleUseHasMatchConv)
	      {
		if( ElectronHasMatchedConvPhot->at(iele) )
		  isPhotConv = 1;
	      }
	    else 
	      {
		if( ElectronDist->at(iele) < eleDist && ElectronDCotTheta->at(iele) < eleDCotTheta )	      
		  isPhotConv = 1;
	      }

	    //barrel/endcap
	    if( fabs( ElectronSCEta->at(iele) ) < eleEta_bar ) 
	      isBarrel = 1;	    
	    if( fabs( ElectronSCEta->at(iele) ) > eleEta_end_min && fabs( ElectronSCEta->at(iele) ) < eleEta_end_max ) 
	      isEndcap = 1;
	    
	    if(isBarrel)
	      {		
		
		if( ElectronMissingHits->at(iele) == eleMissingHits 
		    && isPhotConv == 0
		    && ElectronCombRelIso_bar < eleCombRelIso_bar
		    && ElectronSigmaIEtaIEta->at(iele) < eleSigmaIetaIeta_bar 
		    && fabs(ElectronDeltaPhiTrkSC->at(iele)) < eleDeltaPhiTrkSC_bar 
		    && fabs(ElectronDeltaEtaTrkSC->at(iele)) < eleDeltaEtaTrkSC_bar 
		    )
		  passEleSel = 1;		
		
	      }//end barrel

	    if(isEndcap)
	      {		
		
		if( ElectronMissingHits->at(iele) == eleMissingHits 
		    && isPhotConv == 0
		    && ElectronCombRelIso_end < eleCombRelIso_end
		    && ElectronSigmaIEtaIEta->at(iele) < eleSigmaIetaIeta_end 
		    && fabs(ElectronDeltaPhiTrkSC->at(iele)) < eleDeltaPhiTrkSC_end 
		    && fabs(ElectronDeltaEtaTrkSC->at(iele)) < eleDeltaEtaTrkSC_end 
		    )
		  passEleSel = 1;		
		
	      }//end endcap
	    
	    // Pass Electron Selection
	    if ( passEleSel )
	      {
		v_idx_ele_PtCut_IDISO.push_back(iele);

		//histograms
		FillUserTH1D( "EleAllPt"                  ,    ElectronPt->at(iele)  ); 
		FillUserTH1D( "EleAllEta"                 ,    ElectronEta->at(iele) ); 
	      }
	    
	  }//------> use EWK(WP80) SimpleCutBasedEleID2011
	
	else if (eleAlgorithm == 2) //------> use HEEP variables
	  
	  {
	    
	    // pT pre-cut on ele
	    if( ElectronPt->at(iele) < elePtCut ) continue;
	    
	    // ecal driven	    
	    if( eleUseEcalDrivenHeep && !ElectronHasEcalDrivenSeed->at(iele) )
	      continue;

	    v_idx_ele_PtCut.push_back(iele);

	    int passEleSel = 0;
	    int isBarrel = 0;
	    int isEndcap = 0;
	    
	    if( fabs( ElectronSCEta->at(iele) ) < eleEta_bar ) 
	      isBarrel = 1;
	    
	    if( fabs( ElectronSCEta->at(iele) ) > eleEta_end_min && fabs( ElectronSCEta->at(iele) ) < eleEta_end_max ) 
	      isEndcap = 1;
	    
	    if(isBarrel)
	      {		

		if(fabs(ElectronDeltaEtaTrkSC->at(iele)) < eleDeltaEtaTrkSCHeep_bar 
		   && fabs(ElectronDeltaPhiTrkSC->at(iele)) < eleDeltaPhiTrkSCHeep_bar 
		   && ElectronHoE->at(iele) < eleHoEHeep_bar 
		   && (ElectronE2x5OverE5x5->at(iele) >eleE2x5OverE5x5Heep_bar || ElectronE1x5OverE5x5->at(iele) > eleE1x5OverE5x5Heep_bar ) 
		   && ElectronEcalIsoDR03->at(iele)+ElectronHcalIsoD1DR03->at(iele) 
		   < eleEcalHcalIsoHeep_1_bar + eleEcalHcalIsoHeep_2_bar*ElectronPt->at(iele) 
		   && ElectronTrkIsoDR03->at(iele) <eleTrkIsoHeep_bar
		   )
		  passEleSel = 1;		
		
	      }//end barrel
	
	    if(isEndcap)
	      {		
		
		int passEcalHcalIsoCut=0;
		if(ElectronPt->at(iele) < eleEcalHcalIso_PTthrHeep_end
		   && (ElectronEcalIsoDR03->at(iele)+ElectronHcalIsoD1DR03->at(iele)) < eleEcalHcalIsoHeep_1_end)
		  passEcalHcalIsoCut=1;
		if(ElectronPt->at(iele) > eleEcalHcalIso_PTthrHeep_end
		   && (ElectronEcalIsoDR03->at(iele)+ElectronHcalIsoD1DR03->at(iele)) 
		   < eleEcalHcalIsoHeep_1_end+eleEcalHcalIsoHeep_2_end*(ElectronPt->at(iele)-eleEcalHcalIso_PTthrHeep_end) )
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
	    
	    //Pass User Defined Electron Selection
	    if ( passEleSel )
	      {
		v_idx_ele_PtCut_IDISO.push_back(iele);

		//histograms
		FillUserTH1D( "EleAllPt"                  ,    ElectronPt->at(iele)  ); 
		FillUserTH1D( "EleAllEta"                 ,    ElectronEta->at(iele) ); 
	      }
	    
	  }//------> use HEEP variables

      } // End loop over electrons

    //## Photons
    vector<int> v_idx_pho_PtCut_IDISO;
    
    //Loop over photons
    for(int ipho=0; ipho<PhotonPt->size(); ipho++)
      {      
	// pT pre-cut on photon
	if( PhotonPt->at(ipho) < phoPtCut ) continue;
	
	if( phoAlgorithm == 1) //------> use TightPhoton (Twiki PhotonIDAnalysis) + electron veto
	  {
	    int passPhoSel = 0;
	    int isBarrel = 0;
	    int isEndcap = 0;
	    int isMatchEle = 0;

	    //matched prompt ele
	    if( phoUseMatchPromptEle==1 && phoUsePixelSeed==0 )
	      {
		if( PhotonHasMatchedPromptEle->at(ipho) )
		  isMatchEle = 1;
	      }
	    if( phoUseMatchPromptEle==0 && phoUsePixelSeed==1 )
	      {
		if( PhotonHasPixelSeed->at(ipho) )
		  isMatchEle = 1;
	      }

	    //barrel/endcap
	    if( fabs( PhotonSCeta->at(ipho) ) < phoEta_bar ) 
	      isBarrel = 1;	    
	    if( fabs( PhotonSCeta->at(ipho) ) > phoEta_end_min && fabs( PhotonSCeta->at(ipho) ) < phoEta_end_max ) 
	      isEndcap = 1;
	    
	    if(isBarrel)
	      {		
		
		if( isMatchEle == 0
		    && PhotonHoE->at(ipho) < phoHoE_bar
		    && PhotonSigmaIEtaIEta->at(ipho) < phoSigmaIetaIeta_bar 		    
		    && PhotonEcalIsoDR04->at(ipho) < phoEcalIso_bar_const + phoEcalIso_bar_ptCoeff * PhotonPt->at(ipho)  
		    && PhotonHcalIsoDR04->at(ipho) < phoHcalIso_bar_const + phoHcalIso_bar_ptCoeff * PhotonPt->at(ipho) 
		    && PhotonTrkIsoHollowDR04->at(ipho) < phoTrkIso_bar_const + phoTrkIso_bar_ptCoeff * PhotonPt->at(ipho)
		    )
		  passPhoSel = 1;		
		
	      }//end barrel

	    if(isEndcap)
	      {		

		if( isMatchEle == 0
		    && PhotonHoE->at(ipho) < phoHoE_end
		    && PhotonSigmaIEtaIEta->at(ipho) < phoSigmaIetaIeta_end 		    
		    && PhotonEcalIsoDR04->at(ipho) < phoEcalIso_end_const + phoEcalIso_end_ptCoeff * PhotonPt->at(ipho)  
		    && PhotonHcalIsoDR04->at(ipho) < phoHcalIso_end_const + phoHcalIso_end_ptCoeff * PhotonPt->at(ipho) 
		    && PhotonTrkIsoHollowDR04->at(ipho) < phoTrkIso_end_const + phoTrkIso_end_ptCoeff * PhotonPt->at(ipho)
		    )
		  passPhoSel = 1;		

	      }//end endcap
	    
	    // Pass Photon Selection
	    if ( passPhoSel )
	      {
		v_idx_pho_PtCut_IDISO.push_back(ipho);

		//histograms
		FillUserTH1D( "PhoAllPt"                  ,    PhotonPt->at(ipho)  ); 
		FillUserTH1D( "PhoAllEta"                 ,    PhotonEta->at(ipho) ); 
	      }
	    
	  }//------> use TightPhoton (Twiki PhotonIDAnalysis) + electron veto
      }

    //## Jets
    vector<int> v_idx_jet_PtCut;
    vector<int> v_idx_jet_PtCut_noOverlap_ID;
    vector<int> v_idx_jet_PtCut_noOverlap_ID_TCHEL;

    // Loop over jets
    for(int ijet=0; ijet<JetPt->size(); ijet++)
      {
	//pT pre-cut on reco jets
	if ( JetPt->at(ijet) < jetPtCut ) continue;
	v_idx_jet_PtCut.push_back(ijet);
      }

    //ele-jet overlap
    vector <int> jetFlagsEle(v_idx_jet_PtCut.size(), 0);
    int NjetflaggedEle = 0;
    for (int iele=0; iele<v_idx_ele_PtCut_IDISO.size(); iele++)
      {
	TLorentzVector ele;
        ele.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO[iele]),
			 ElectronEta->at(v_idx_ele_PtCut_IDISO[iele]),
			 ElectronPhi->at(v_idx_ele_PtCut_IDISO[iele]),0);
	TLorentzVector jet;
	double minDR=9999.;
	int ijet_minDR = -1;
        for(int ijet=0; ijet<v_idx_jet_PtCut.size(); ijet++)
          {
	    if ( jetFlagsEle[ijet] == 1 )
	      continue;
            jet.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut[ijet]),
			     JetEta->at(v_idx_jet_PtCut[ijet]),
			     JetPhi->at(v_idx_jet_PtCut[ijet]),
			     JetEnergy->at(v_idx_jet_PtCut[ijet]) );
	    double DR = jet.DeltaR(ele);
	    if (DR<minDR)
	      {
		minDR = DR;
		ijet_minDR = ijet;
	      }
	  }
	if ( minDR < jet_ele_DeltaRcut 
	     && ijet_minDR > -1 )
	  {
	    jetFlagsEle[ijet_minDR] = 1;
	    NjetflaggedEle++;
	  }
      }

    //muon-jet overlap
    vector <int> jetFlagsMuon(v_idx_jet_PtCut.size(), 0);
    int NjetflaggedMuon = 0;
    for (int imu=0; imu<v_idx_muon_PtCut_IDISO.size(); imu++)
      {
	TLorentzVector mu;
        mu.SetPtEtaPhiM(MuonPt->at(v_idx_muon_PtCut_IDISO[imu]),
			MuonEta->at(v_idx_muon_PtCut_IDISO[imu]),
			MuonPhi->at(v_idx_muon_PtCut_IDISO[imu]),0);
	TLorentzVector jet;
	double minDR=9999.;
	int ijet_minDR = -1;
        for(int ijet=0; ijet<v_idx_jet_PtCut.size(); ijet++)
          {
	    if ( jetFlagsMuon[ijet] == 1 )
	      continue;
            jet.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut[ijet]),
			     JetEta->at(v_idx_jet_PtCut[ijet]),
			     JetPhi->at(v_idx_jet_PtCut[ijet]),
			     JetEnergy->at(v_idx_jet_PtCut[ijet]) );
	    double DR = jet.DeltaR(mu);
	    if (DR<minDR)
	      {
		minDR = DR;
		ijet_minDR = ijet;
	      }
	  }
	if ( minDR < jet_muon_DeltaRcut 
	     && ijet_minDR > -1 )
	  {
	    jetFlagsMuon[ijet_minDR] = 1;
	    NjetflaggedMuon++;
	  }
      }

    //photon-jet overlap
    vector <int> jetFlagsPho(v_idx_jet_PtCut.size(), 0);
    int NjetflaggedPho = 0;
    for (int ipho=0; ipho<v_idx_pho_PtCut_IDISO.size(); ipho++)
      {                   
	TLorentzVector pho;
        pho.SetPtEtaPhiM(PhotonPt->at(v_idx_pho_PtCut_IDISO[ipho]),
			 PhotonEta->at(v_idx_pho_PtCut_IDISO[ipho]),
			 PhotonPhi->at(v_idx_pho_PtCut_IDISO[ipho]),0);
	TLorentzVector jet;
	double minDR=9999.;
	int ijet_minDR = -1;
        for(int ijet=0; ijet<v_idx_jet_PtCut.size(); ijet++)
          {
	    if ( jetFlagsPho[ijet] == 1 )
	      continue;
            jet.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut[ijet]),
			     JetEta->at(v_idx_jet_PtCut[ijet]),
			     JetPhi->at(v_idx_jet_PtCut[ijet]),
			     JetEnergy->at(v_idx_jet_PtCut[ijet]) );
	    double DR = jet.DeltaR(pho);
	    if (DR<minDR)
	      {
		minDR = DR;
		ijet_minDR = ijet;
	      }
	  }
	if ( minDR < jet_pho_DeltaRcut 
	     && ijet_minDR > -1 )
	  {
	    jetFlagsPho[ijet_minDR] = 1;
	    NjetflaggedPho++;
	  }
      }

    // Pass Jet Selection
    for(int ijet=0; ijet<v_idx_jet_PtCut.size(); ijet++) //pT pre-cut + no overlaps with electrons + jetID
      {
	bool passjetID = true;
	
	if( jetIDloose && !jetIDtight )
	  {
	    passjetID = JetPassLooseID->at(v_idx_jet_PtCut[ijet]);
	  }
	if( jetIDtight )
	  {
	    passjetID = JetPassTightID->at(v_idx_jet_PtCut[ijet]);
	  }

	if( jetFlagsEle[ijet] == 0                                      /* NO overlap with electrons */ 
	    && jetFlagsMuon[ijet] == 0                                  /* NO overlap with muons */	    
	    //&& jetFlagsPho[ijet] == 0                                 /* NO overlap with photons */	    
	    && passjetID == true                                        /* pass JetID */
	    && fabs( JetEta->at(v_idx_jet_PtCut[ijet]) ) < jetEtaCut 	/* pass Jet Eta cut */    
	    )                         
	  {
	    v_idx_jet_PtCut_noOverlap_ID.push_back(v_idx_jet_PtCut[ijet]);

	    //histograms
	    FillUserTH1D( "JetAllPt"                  ,    JetPt->at(ijet)  ); 
	    FillUserTH1D( "JetAllEta"                 ,    JetEta->at(ijet) ); 
	  }

	if( jetFlagsEle[ijet] == 0                                       /* NO overlap with electrons */ 
	    && jetFlagsMuon[ijet] == 0                                   /* NO overlap with muons */	    
	    //&& jetFlagsPho[ijet] == 0                                  /* NO overlap with photons */	    
	    && passjetID == true                                         /* pass JetID */
	    && fabs( JetEta->at(v_idx_jet_PtCut[ijet]) ) < jetEtaCut 	 /* pass Jet Eta cut */    
            && fabs( JetTCHE->at(v_idx_jet_PtCut[ijet]) ) > jetTCHELCut  /* TrackCountingHighEfficiency loose b-tag (see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagPerformance) */
	    )
	  {                         
	    v_idx_jet_PtCut_noOverlap_ID_TCHEL.push_back(v_idx_jet_PtCut[ijet]);	

	    //histograms
	    FillUserTH1D( "BJetAllPt"                  ,    JetPt->at(ijet)  ); 
	    FillUserTH1D( "BJetAllEta"                 ,    JetEta->at(ijet) ); 
	  }
      } // End loop over jets

    //## MET scale uncert.
    if( JetEnergyScale != 1 || (!isData && doJetOversmearing) )
      { //use fix JES scaling passed from cut file

	TVector2 v_MET_old;
	TVector2 v_MET_new;

	//use only good jets (after electron-jet and muon-jet overlap) for re-doing MET
	for(int ijet=0; ijet<v_idx_jet_PtCut_noOverlap_ID.size(); ijet++)
	  {
	    TVector2 v_jet_pt_old;
	    TVector2 v_jet_pt_new;
	    v_jet_pt_old.SetMagPhi( JetPt->at(v_idx_jet_PtCut_noOverlap_ID[ijet])/JetEnergyScale/JetOversmearingFactor->at(v_idx_jet_PtCut_noOverlap_ID[ijet]) , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[ijet]) );
	    v_jet_pt_new.SetMagPhi( JetPt->at(v_idx_jet_PtCut_noOverlap_ID[ijet]) , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[ijet]) );
	    //pT pre-cut on reco jets
	    if ( v_jet_pt_old.Mod() < jet_PtCut_forMetScale ) continue;
	    v_MET_new += v_jet_pt_old - v_jet_pt_new;
	  }

	//for MET energy scale
	v_MET_old.SetMagPhi( thisMET , thisMETPhi );
	v_MET_new += v_MET_old;
 	thisMET = v_MET_new.Mod();
 	thisMETPhi = v_MET_new.Phi();

	//## for debug
	//double METscale_diff = thisMET - v_MET_old.Mod() ;
	// 	cout << "old MET = " << v_MET_old.Mod()
	// 	     << " " << "new MET = " << thisMET
	// 	     << " " << "new-old = " << METscale_diff
	// 	     << endl;
	//CreateAndFillUserTH1D("h1_METscale_diff", 400, -100, 100, METscale_diff);
	//##
      }

    // vertices
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

    //Printout Gen level
//     if(isData==0 && (event==26937 || event==6534 || event==39495 || event==16006 || event==34295 || event==31772 ) )
//       {
// 	cout << endl;
// 	cout << "----------------------" << endl;
// 	cout << "event: " << event<< " , GenMET: " << GenMETTrue->at(0) << endl;


// 	for (int genp=0; genp<GenParticlePdgId->size(); genp++)
// 	  {	   
// 	    double mass = 0;
// 	    TLorentzVector particle;
// 	    particle.SetPtEtaPhiE(GenParticlePt->at(genp) , GenParticleEta->at(genp) , GenParticlePhi->at(genp) , GenParticleEnergy->at(genp) );

// 	    //--- printout
// 	    cout << "idx: " << genp
// 		 << " GenParticlePdgId: " << GenParticlePdgId->at(genp)
// 		 << " GenParticleMotherIndex: " << GenParticleMotherIndex->at(genp)
// 		 << " GenParticleStatus: " << GenParticleStatus->at(genp)
// 		 << " GenParticlePt: " << GenParticlePt->at(genp)
// 		 << " GenParticleEta: " << GenParticleEta->at(genp)
// 		 << " GenParticlePhi: " << GenParticlePhi->at(genp)
// 		 << " GenPArticleMass: " << particle.M()
// 		 << endl;	 
// 	    //---
// 	  }
//       }

    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();

    // Set the value of the variableNames listed in the cutFile to their current value
    //fillVariableWithValue( "PassJSON", passJSON(run, ls, isData) );    

    // Evaluate cuts (but do not apply them)
    evaluateCuts();

    int Num_AG = 0;
    int Num_AG_Daughters = 0;	      
    TLorentzVector v_AG;
    TLorentzVector v_AG_boosted;
    TLorentzVector v_AG_q;
    TLorentzVector v_AG_qbar;
    int            flavour_AG_q = 0;
    int            flavour_AG_qbar = 0;
    int            flavour_W_lepton = 0;
    int            flavour_W_neutrino = 0;

    int Num_W = 0; 
    int Wcharge = 0;
    int index_W = -1;
    TLorentzVector v_W;
    TLorentzVector v_W_l;
    TLorentzVector v_W_nu;
    double MT_W = 0;

    int T_EXISTS = 0;
    int TBAR_EXISTS = 0;

    if(isData==0)
      {
	for (int genp=0; genp<GenParticlePdgId->size(); genp++)
	  {

	    //--- printout
	    //  	    	     	    cout << "idx: " << genp
	    //  	    	     		 << " GenParticlePdgId: " << GenParticlePdgId->at(genp)
	    //  	    	     		 << " GenParticleMotherIndex: " << GenParticleMotherIndex->at(genp)
	    //  	    	     		 << " GenParticleStatus: " << GenParticleStatus->at(genp)
	    //  	    	     		 << " GenParticlePt: " << GenParticlePt->at(genp)
	    //  	    	     		 << endl;	    
	    //---

	    TLorentzVector current_p4;
	    current_p4.SetPtEtaPhiE(GenParticlePt->at(genp),
				    GenParticleEta->at(genp),
				    GenParticlePhi->at(genp),
				    GenParticleEnergy->at(genp)
				    );
	    
	    //AG
	    if(GenParticlePdgId->at(genp) == PDGID_AG && GenParticleStatus->at(genp) == 3)
	      {
		Num_AG++;
		Num_AG_Daughters = GenParticleNumDaught->at(genp);
		v_AG = current_p4;
	      }//

	    //qq from AG
	    if( GenParticleMotherIndex->at(genp) != -1 )
	      {  
		if( abs( GenParticlePdgId->at( GenParticleMotherIndex->at(genp) ) ) == PDGID_AG && GenParticleStatus->at(genp) == 3)
		  {
		    //First TOP = assuming that top comes always before than bottom in the gen particle list (should be always true)
		    if( GenParticlePdgId->at(genp) == PDGID_T )
		      {//t
			T_EXISTS = 1;
			v_AG_q = current_p4;
			flavour_AG_q = GenParticlePdgId->at(genp);			
		      }
		    if( GenParticlePdgId->at(genp) == PDGID_TBAR )
		      {//tbar
			TBAR_EXISTS = 1;
			v_AG_qbar = current_p4;
			flavour_AG_qbar = GenParticlePdgId->at(genp);
		      }

		    //Then D, U, S, C
		    if( ( GenParticlePdgId->at(genp) == PDGID_D 
			|| GenParticlePdgId->at(genp) == PDGID_U  
			|| GenParticlePdgId->at(genp) == PDGID_S  
			  || GenParticlePdgId->at(genp) == PDGID_C ) 
			&& T_EXISTS == 0
			) 
		      {//quark
			v_AG_q = current_p4;
			flavour_AG_q = GenParticlePdgId->at(genp);
		      }		    
		    if( ( GenParticlePdgId->at(genp) == PDGID_DBAR 
			|| GenParticlePdgId->at(genp) == PDGID_UBAR  
			|| GenParticlePdgId->at(genp) == PDGID_SBAR  
			  || GenParticlePdgId->at(genp) == PDGID_CBAR )
			&& TBAR_EXISTS == 0
			) 
		      {//antiquark
			v_AG_qbar = current_p4;
			flavour_AG_qbar = GenParticlePdgId->at(genp);
		      }
		    		   
		    //Finally BOTTOM
		    if( GenParticlePdgId->at(genp) == PDGID_B && T_EXISTS==0 )
		      {//b
			v_AG_q = current_p4;
			flavour_AG_q = GenParticlePdgId->at(genp);						
		      }
		    if( GenParticlePdgId->at(genp) == PDGID_BBAR && TBAR_EXISTS==0 )
		      {//bbar
			v_AG_qbar = current_p4;
			flavour_AG_qbar = GenParticlePdgId->at(genp);
		      }		    

		  }
	      }//

	    //W
	    if( abs(GenParticlePdgId->at(genp)) == abs(PDGID_Wplus) 
		&& GenParticleStatus->at(genp) == 3 
		&& abs( GenParticlePdgId->at( GenParticleMotherIndex->at(genp) ) ) != PDGID_AG 
		&& abs( GenParticlePdgId->at( GenParticleMotherIndex->at(genp) ) ) != abs(PDGID_T) 
		)
	      {
		Num_W++;
		v_W = current_p4;
		index_W = genp;

		//W+
		if(GenParticlePdgId->at(genp) == PDGID_Wplus)
		  Wcharge = 1;

		//W-
		if(GenParticlePdgId->at(genp) == PDGID_Wminus)
		  Wcharge = -1;
	      }

	    //leptons from W
	    if( GenParticleMotherIndex->at(genp) != -1 )
	      {  
		if( GenParticleMotherIndex->at(genp) == index_W )
		  {
		    if( abs(GenParticlePdgId->at(genp)) == abs(PDGID_ELE)
			|| abs(GenParticlePdgId->at(genp)) == abs(PDGID_MU)  
			|| abs(GenParticlePdgId->at(genp)) == abs(PDGID_TAU)  
			) 
		      {//lepton
			v_W_l = current_p4;
			flavour_W_lepton = GenParticlePdgId->at(genp);
		      }

		    if( abs(GenParticlePdgId->at(genp)) == abs(PDGID_NUELE)
			|| abs(GenParticlePdgId->at(genp)) == abs(PDGID_NUMU)  
			|| abs(GenParticlePdgId->at(genp)) == abs(PDGID_NUTAU)  
			) 
		      {//neutrino
			v_W_nu = current_p4;
			flavour_W_neutrino = GenParticlePdgId->at(genp);
		      }
		  }
	      }//

	    //MT
	    MT_W = sqrt( 2 * v_W_l.Pt() * v_W_nu.Pt() * (1 - cos(v_W_l.DeltaPhi(v_W_nu) ) ) );	    
	    
	  }//end loop over gen particles


// 	v_AG_boosted = v_AG;
// 	TVector3 v_beta_AG;
// 	v_beta_AG = v_AG.Vect();
// 	v_beta_AG.SetMag( v_AG.Beta() );		
// 	v_AG_boosted.Boost(v_beta_AG);

// 	cout << "v_AG.E(): " << v_AG.E() << endl;
// 	cout << "v_AG_boosted.E(): " << v_AG_boosted.E() << endl;


      }

    // Fill histograms and do analysis based on cut evaluation

    if(flavour_AG_q != -(flavour_AG_qbar))
      {
	cout << "---------" << endl;
	for (int genp1=0; genp1<GenParticlePdgId->size(); genp1++)
	  {			
	    cout << "idx: " << genp1
		 << " GenParticlePdgId: " << GenParticlePdgId->at(genp1)
		 << " GenParticleMotherIndex: " << GenParticleMotherIndex->at(genp1)
		 << " GenParticleStatus: " << GenParticleStatus->at(genp1)
		 << " GenParticlePt: " << GenParticlePt->at(genp1)
		 << endl;	    
	  }	
      }

    FillUserTH1D("h1_Num_AG", Num_AG);
    FillUserTH1D("h1_Num_W", Num_W);
    FillUserTH1D("h1_Charge_W", Wcharge);

    FillUserTH1D("h1_Mass_AG", v_AG.M() );
    FillUserTH1D("h1_Mass_AG_fromqq", (v_AG_q+v_AG_qbar).M() );
    FillUserTH1D("h1_Mass_W", v_W.M() );
    FillUserTH1D("h1_Mass_W_fromlnu", (v_W_l+v_W_nu).M() );
    FillUserTH1D("h1_MT_W", MT_W );

    FillUserTH1D("h1_Pt_AG", v_AG.Pt() );
    FillUserTH1D("h1_Pt_W", v_W.Pt() );
    if(	  Wcharge == 1  )
      FillUserTH1D("h1_Pt_Wplus", v_W.Pt() );
    if(	  Wcharge == -1  )
      FillUserTH1D("h1_Pt_Wminus", v_W.Pt() );
    
    FillUserTH1D("h1_Eta_AG", v_AG.Eta() );
    FillUserTH1D("h1_Eta_W", v_W.Eta() );
    if(	  Wcharge == 1  )
      FillUserTH1D("h1_Eta_Wplus", v_W.Eta() );
    if(	  Wcharge == -1  )
      FillUserTH1D("h1_Eta_Wminus", v_W.Eta() );

    FillUserTH1D("h1_Pt_AG_q", v_AG_q.Pt() );
    FillUserTH1D("h1_Pt_AG_qbar", v_AG_qbar.Pt() );
    FillUserTH1D("h1_Pt_W_l", v_W_l.Pt() );
    FillUserTH1D("h1_Pt_W_nu", v_W_nu.Pt() );
    
    FillUserTH1D("h1_Eta_AG_q", v_AG_q.Eta() );
    FillUserTH1D("h1_Eta_AG_qbar", v_AG_qbar.Eta() );
    FillUserTH1D("h1_Eta_W_l", v_W_l.Eta() );
    FillUserTH1D("h1_Eta_W_nu", v_W_nu.Eta() );

    if( flavour_AG_q == -(flavour_AG_qbar) )
      {
	h1_DecayProducts_AG->AddBinContent( abs(flavour_AG_q) );
	h1_DecayProducts_AG->AddBinContent( 7 );
      }
    else
      {
	cout << flavour_AG_q << endl;
	cout << flavour_AG_qbar << endl;
	cout << "THERE IS A PROBLEM!!!!!!!!!!!!!!!!!!!!" << endl;
      }

    if( fabs(flavour_W_lepton) == abs(PDGID_ELE) && fabs(flavour_W_neutrino) == abs(PDGID_NUELE) )
      h1_DecayProducts_W->AddBinContent( 1 );
    if( fabs(flavour_W_lepton) == abs(PDGID_MU) && fabs(flavour_W_neutrino) == abs(PDGID_NUMU) )
      h1_DecayProducts_W->AddBinContent( 2 );
    if( fabs(flavour_W_lepton) == abs(PDGID_TAU) && fabs(flavour_W_neutrino) == abs(PDGID_NUTAU) )
      h1_DecayProducts_W->AddBinContent( 3 );
    if( (fabs(flavour_W_lepton) == abs(PDGID_ELE) && fabs(flavour_W_neutrino) == abs(PDGID_NUELE)) 
	|| (fabs(flavour_W_lepton) == abs(PDGID_MU) && fabs(flavour_W_neutrino) == abs(PDGID_NUMU)) 
	|| (fabs(flavour_W_lepton) == abs(PDGID_TAU) && fabs(flavour_W_neutrino) == abs(PDGID_NUTAU)) )
      h1_DecayProducts_W->AddBinContent( 4 );

    FillUserTH1D("h1_Num_AG_Daughters", Num_AG_Daughters);		

    FillUserTH1D("h1_GenMET", GenMETTrue->at(0) );


    //### Calculate acceptance and efficiency
    
    
    //Electrons
    if( fabs(flavour_W_lepton) == abs(PDGID_ELE) && fabs(flavour_W_neutrino) == abs(PDGID_NUELE) )
      {	
	bool isInAcceptance             = false;
	bool isMatchedToRecoObject      = false;
	bool isMatchedToRecoIDISOObject = false;

	//acceptance cuts
	if( v_W_l.Pt() > elePtCut 
	    && ( 
		fabs(v_W_l.Eta()) < eleEta_bar 
		|| 
		( fabs(v_W_l.Eta()) > eleEta_end_min && fabs(v_W_l.Eta()) < eleEta_end_max ) 
		) 
	    )
	  {
	    isInAcceptance = true;
	  }

	//matched with reco object	
	int iele_minDR=-1;
	double minDR=9999.;	    
	for (int iele=0; iele<v_idx_ele_PtCut.size(); iele++)
	  {
	    TLorentzVector ele;
	    ele.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut[iele]),
			     ElectronEta->at(v_idx_ele_PtCut[iele]),
			     ElectronPhi->at(v_idx_ele_PtCut[iele]),0);
	    double DR = v_W_l.DeltaR(ele);
	    if (DR<minDR)
	      {
		minDR = DR;
		iele_minDR = iele;
	      }
	  }	
	double eleGen_DeltaRcut =0.1;
	CreateAndFillUserTH1D("h1_minDR_e", 1000, 0, 10, minDR);	    	
	if ( minDR < eleGen_DeltaRcut && iele_minDR > -1 )
	  {
	    CreateAndFillUserTH1D("h1_ErecoOverEgen_e", 100, 0, 3, ElectronEnergy->at(v_idx_ele_PtCut[iele_minDR]) / v_W_l.E() );	    	
	    isMatchedToRecoObject = true;
	  }      			     

	//matched with reco+ID+ISO object	
	int iele_minDR_IDISO=-1;
	double minDR_IDISO=9999.;	    
	for (int iele=0; iele<v_idx_ele_PtCut_IDISO.size(); iele++)
	  {
	    TLorentzVector ele;
	    ele.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO[iele]),
			     ElectronEta->at(v_idx_ele_PtCut_IDISO[iele]),
			     ElectronPhi->at(v_idx_ele_PtCut_IDISO[iele]),0);
	    double DR = v_W_l.DeltaR(ele);
	    if (DR<minDR_IDISO)
	      {
		minDR_IDISO = DR;
		iele_minDR_IDISO = iele;
	      }
	  }	
	double eleGen_DeltaRcut_IDISO =0.1;
	CreateAndFillUserTH1D("h1_minDR_e_IDISO", 1000, 0, 10, minDR_IDISO);	    	
	if ( minDR_IDISO < eleGen_DeltaRcut_IDISO && iele_minDR_IDISO > -1 )
	  {
	    CreateAndFillUserTH1D("h1_ErecoOverEgen_e_IDISO", 100, 0, 3, ElectronEnergy->at(v_idx_ele_PtCut[iele_minDR_IDISO]) / v_W_l.E() );	    	
	    isMatchedToRecoIDISOObject = true;
	  }      			     

	//Distance ele-quark
	double minDR_elequark = min( v_W_l.DeltaR(v_AG_q) , v_W_l.DeltaR(v_AG_qbar) );	    

	//Plots
 	FillUserTH1D("h1_Pt_W_e", v_W_l.Pt() );
 	FillUserTH1D("h1_Eta_W_e", v_W_l.Eta() );
	FillUserTH1D("h1_minDR_elequark", minDR_elequark );
 	FillUserTH1D("h1_nVtx_e", VertexChi2->size()  );	

	if( isInAcceptance )
	  {
 	    FillUserTH1D("h1_Pt_W_e_gencuts", v_W_l.Pt() );
 	    FillUserTH1D("h1_Eta_W_e_gencuts", v_W_l.Eta() );	    
	    FillUserTH1D("h1_minDR_elequark_gencuts", minDR_elequark );
	    FillUserTH1D("h1_nVtx_e_gencuts", VertexChi2->size()  );	
	  }
	if( isInAcceptance && isMatchedToRecoObject )
	  {
 	    FillUserTH1D("h1_Pt_W_e_recomatched", v_W_l.Pt() );
 	    FillUserTH1D("h1_Eta_W_e_recomatched", v_W_l.Eta() );	    	    
	    FillUserTH1D("h1_minDR_elequark_recomatched", minDR_elequark );
	    FillUserTH1D("h1_nVtx_e_recomatched", VertexChi2->size()  );	
	  }
	if( isInAcceptance && isMatchedToRecoObject && isMatchedToRecoIDISOObject )
	  {
 	    FillUserTH1D("h1_Pt_W_e_recomatched_IDISO", v_W_l.Pt() );
 	    FillUserTH1D("h1_Eta_W_e_recomatched_IDISO", v_W_l.Eta() );	   
	    FillUserTH1D("h1_minDR_elequark_recomatched_IDISO", minDR_elequark );
	    FillUserTH1D("h1_nVtx_e_recomatched_IDISO", VertexChi2->size()  );	 	    
	  }

      }


    //Muons
    if( fabs(flavour_W_lepton) == abs(PDGID_MU) && fabs(flavour_W_neutrino) == abs(PDGID_NUMU) )
      {	
	bool isInAcceptance             = false;
	bool isMatchedToRecoObject      = false;
	bool isMatchedToRecoIDISOObject = false;

	//acceptance cuts
	if( v_W_l.Pt() > muPtCut 
	    && ( 
		fabs(v_W_l.Eta()) < muEta_bar 
		|| 
		( fabs(v_W_l.Eta()) > muEta_end_min && fabs(v_W_l.Eta()) < muEta_end_max ) 
		) 
	    )
	  {
	    isInAcceptance = true;
	  }

	//matched with reco object	
	int imuon_minDR=-1;
	double minDR=9999.;	    
	for (int imuon=0; imuon<v_idx_muon_PtCut.size(); imuon++)
	  {
	    TLorentzVector muon;
	    muon.SetPtEtaPhiM(MuonPt->at(v_idx_muon_PtCut[imuon]),
			     MuonEta->at(v_idx_muon_PtCut[imuon]),
			     MuonPhi->at(v_idx_muon_PtCut[imuon]),0);
	    double DR = v_W_l.DeltaR(muon);
	    if (DR<minDR)
	      {
		minDR = DR;
		imuon_minDR = imuon;
	      }
	  }	
	double muonGen_DeltaRcut =0.1;
	CreateAndFillUserTH1D("h1_minDR_mu", 1000, 0, 10, minDR);	    	
	if ( minDR < muonGen_DeltaRcut && imuon_minDR > -1 )
	  {
	    CreateAndFillUserTH1D("h1_ErecoOverEgen_mu", 100, 0, 3, MuonP->at(v_idx_muon_PtCut[imuon_minDR]) / v_W_l.P() );	    	
	    isMatchedToRecoObject = true;
	  }      			     

	//matched with reco+ID+ISO object	
	int imuon_minDR_IDISO=-1;
	double minDR_IDISO=9999.;	    
	for (int imuon=0; imuon<v_idx_muon_PtCut_IDISO.size(); imuon++)
	  {
	    TLorentzVector muon;
	    muon.SetPtEtaPhiM(MuonPt->at(v_idx_muon_PtCut_IDISO[imuon]),
			     MuonEta->at(v_idx_muon_PtCut_IDISO[imuon]),
			     MuonPhi->at(v_idx_muon_PtCut_IDISO[imuon]),0);
	    double DR = v_W_l.DeltaR(muon);
	    if (DR<minDR_IDISO)
	      {
		minDR_IDISO = DR;
		imuon_minDR_IDISO = imuon;
	      }
	  }	
	double muonGen_DeltaRcut_IDISO =0.1;
	CreateAndFillUserTH1D("h1_minDR_mu_IDISO", 1000, 0, 10, minDR_IDISO);	    	
	if ( minDR_IDISO < muonGen_DeltaRcut_IDISO && imuon_minDR_IDISO > -1 )
	  {
	    CreateAndFillUserTH1D("h1_ErecoOverEgen_mu_IDISO", 100, 0, 3, MuonP->at(v_idx_muon_PtCut[imuon_minDR_IDISO]) / v_W_l.P() );	    	
	    isMatchedToRecoIDISOObject = true;
	  }      			     

	//Distance muon-quark
	double minDR_muonquark = min( v_W_l.DeltaR(v_AG_q) , v_W_l.DeltaR(v_AG_qbar) );	    

	//Plots
 	FillUserTH1D("h1_Pt_W_mu", v_W_l.Pt() );
 	FillUserTH1D("h1_Eta_W_mu", v_W_l.Eta() );
	FillUserTH1D("h1_minDR_muonquark", minDR_muonquark );
	FillUserTH1D("h1_nVtx_mu", VertexChi2->size()  );	 	    
	
	if( isInAcceptance )
	  {
 	    FillUserTH1D("h1_Pt_W_mu_gencuts", v_W_l.Pt() );
 	    FillUserTH1D("h1_Eta_W_mu_gencuts", v_W_l.Eta() );	    
	    FillUserTH1D("h1_minDR_muonquark_gencuts", minDR_muonquark );
	    FillUserTH1D("h1_nVtx_mu_gencuts", VertexChi2->size()  );	 	    
	  }
	if( isInAcceptance && isMatchedToRecoObject )
	  {
 	    FillUserTH1D("h1_Pt_W_mu_recomatched", v_W_l.Pt() );
 	    FillUserTH1D("h1_Eta_W_mu_recomatched", v_W_l.Eta() );	    	    
	    FillUserTH1D("h1_minDR_muonquark_recomatched", minDR_muonquark );
	    FillUserTH1D("h1_nVtx_mu_recomatched", VertexChi2->size()  );	 	    
	  }
	if( isInAcceptance && isMatchedToRecoObject && isMatchedToRecoIDISOObject )
	  {
 	    FillUserTH1D("h1_Pt_W_mu_recomatched_IDISO", v_W_l.Pt() );
 	    FillUserTH1D("h1_Eta_W_mu_recomatched_IDISO", v_W_l.Eta() );	   
	    FillUserTH1D("h1_minDR_muonquark_recomatched_IDISO", minDR_muonquark ); 	    
	    FillUserTH1D("h1_nVtx_mu_recomatched_IDISO", VertexChi2->size()  );	 	    
	  }

      }
    
    //### DiJet mass resolution studies

    //---------------------

    //1) matching quark with closest genjet
    //2) using 2 leading in pT genjets or reco jets
    //3) using fat genjets or fat recojets

    //---------------------

    //1) matching quark with closest genjet
    bool q_isMatchedToGenJet = false;
    bool qbar_isMatchedToGenJet = false;
    double q_genjet_DeltaRcut = 0.1;

    int igenJet_minDR_q_genJet=-1;
    double minDR_q_genJet=9999.;	    
    int igenJet_minDR_qbar_genJet=-1;
    double minDR_qbar_genJet=9999.;	    
	
    for(int igenjet=0; igenjet < GenJetPt->size(); igenjet++)
      {

	if( GenJetEMF->at(igenjet)>0.99 )
	  continue;
	
	TLorentzVector genjet;
	genjet.SetPtEtaPhiE(GenJetPt->at(igenjet),
			    GenJetEta->at(igenjet),
			    GenJetPhi->at(igenjet),
			    GenJetEnergy->at(igenjet)
			    );
	double DR_q = v_AG_q.DeltaR(genjet);
	if (DR_q<minDR_q_genJet)
	  {
	    minDR_q_genJet = DR_q;
	    igenJet_minDR_q_genJet = igenjet;
	  }
	
	double DR_qbar = v_AG_qbar.DeltaR(genjet);
	if (DR_qbar<minDR_qbar_genJet)
	  {
	    minDR_qbar_genJet = DR_qbar;
	    igenJet_minDR_qbar_genJet = igenjet;
	  }

      }
    
    if( igenJet_minDR_qbar_genJet!=-1 && igenJet_minDR_q_genJet !=-1 )
      {

	//genjet closest to gen quark

	TLorentzVector genjet_q_match;
	TLorentzVector genjet_qbar_match;
	TLorentzVector axigluon_genjets_match;
	genjet_q_match.SetPtEtaPhiE(GenJetPt->at(igenJet_minDR_q_genJet),
				    GenJetEta->at(igenJet_minDR_q_genJet),
				    GenJetPhi->at(igenJet_minDR_q_genJet),
				    GenJetEnergy->at(igenJet_minDR_q_genJet)
				    );
	genjet_qbar_match.SetPtEtaPhiE(GenJetPt->at(igenJet_minDR_qbar_genJet),
				       GenJetEta->at(igenJet_minDR_qbar_genJet),
				       GenJetPhi->at(igenJet_minDR_qbar_genJet),
				       GenJetEnergy->at(igenJet_minDR_qbar_genJet)
				       );	
	axigluon_genjets_match = genjet_q_match + genjet_qbar_match;

	CreateAndFillUserTH1D("h1_minDR_q_genjet", 1000, 0, 10, minDR_q_genJet);	    	
	CreateAndFillUserTH1D("h1_minDR_q_genjet", 1000, 0, 10, minDR_qbar_genJet);	    	
	CreateAndFillUserTH1D("h1_EgenjetOverEgen_q", 100, 0, 3, genjet_q_match.E() / v_AG_q.E() );	    	
	CreateAndFillUserTH1D("h1_EgenjetOverEgen_q", 100, 0, 3, genjet_qbar_match.E() / v_AG_qbar.E() );	    		    
	CreateAndFillUserTH1D("h1_Mjj_genJet_matched_quarks", 1000, 0, 2000, axigluon_genjets_match.M() );	    	

	if( flavour_AG_q == -(flavour_AG_qbar) )
	  {
	    
	    if( flavour_AG_q != fabs(PDGID_T)  )
	      {
		CreateAndFillUserTH1D("h1_minDR_notTop_genjet", 1000, 0, 10, minDR_q_genJet);	    	
		CreateAndFillUserTH1D("h1_minDR_notTop_genjet", 1000, 0, 10, minDR_qbar_genJet);	    	
		CreateAndFillUserTH1D("h1_EgenjetOverEgen_notTop", 100, 0, 3, GenJetEnergy->at(igenJet_minDR_q_genJet) / v_AG_q.E() );	    	
		CreateAndFillUserTH1D("h1_EgenjetOverEgen_notTop", 100, 0, 3, GenJetEnergy->at(igenJet_minDR_qbar_genJet) / v_AG_qbar.E() );	    		    
		CreateAndFillUserTH2D("h2_minDR_notTop_genjet_vs_EgenjetOverEgen_notTop", 
				      1000, 0, 3, 1000, 0, 10, GenJetEnergy->at(igenJet_minDR_q_genJet) / v_AG_q.E() , minDR_q_genJet );	 
		CreateAndFillUserTH2D("h2_minDR_notTop_genjet_vs_EgenjetOverEgen_notTop", 
				      1000, 0, 3, 1000, 0, 10, GenJetEnergy->at(igenJet_minDR_qbar_genJet) / v_AG_qbar.E() , minDR_qbar_genJet );	 

		CreateAndFillUserTH1D("h1_Mjj_genJet_matched_notTop", 1000, 0, 2000, axigluon_genjets_match.M() );	    	
		CreateAndFillUserTH2D("h2_DeltaRquarks_vs_Mjj_genJet_matched_notTop", 1000, 0, 2000, 100, 0, 10, axigluon_genjets_match.M(),v_AG_q.DeltaR(v_AG_qbar) );	 
	      }
	    
	    if( flavour_AG_q == fabs(PDGID_T) )
	      {
		CreateAndFillUserTH1D("h1_minDR_Top_genjet", 1000, 0, 10, minDR_q_genJet);	    	
		CreateAndFillUserTH1D("h1_minDR_Top_genjet", 1000, 0, 10, minDR_qbar_genJet);	    	
		CreateAndFillUserTH1D("h1_EgenjetOverEgen_Top", 100, 0, 3, GenJetEnergy->at(igenJet_minDR_q_genJet) / v_AG_q.E() );	    	
		CreateAndFillUserTH1D("h1_EgenjetOverEgen_Top", 100, 0, 3, GenJetEnergy->at(igenJet_minDR_qbar_genJet) / v_AG_qbar.E() );	    		    
		CreateAndFillUserTH1D("h1_Mjj_genJet_matched_Top", 1000, 0, 2000, axigluon_genjets_match.M() );	    	
	      }

	  }
	
      }
    

    //2) using 2 leading in pT genjets or reco jets
    TLorentzVector genjet_1stPt;
    TLorentzVector genjet_2ndPt;
    TLorentzVector axigluon_1st2nd_genjets;
    TLorentzVector recojet_1stPt;
    TLorentzVector recojet_2ndPt;
    TLorentzVector axigluon_1st2nd_recojets;

    vector<int> v_idx_good_genjets;
    for(int igenjet=0; igenjet < GenJetPt->size(); igenjet++)
      {
	
	TLorentzVector genjet;
	genjet.SetPtEtaPhiE(GenJetPt->at(igenjet),
			    GenJetEta->at(igenjet),
			    GenJetPhi->at(igenjet),
			    GenJetEnergy->at(igenjet)
			    );
	
	if( genjet.DeltaR(v_W_l) < 0.3 
	    || genjet.Pt() < 30 || fabs(genjet.Eta()) > 2.4 )
	  continue;
	
	v_idx_good_genjets.push_back(igenjet);
	
      }

    if( v_idx_good_genjets.size() >= 2)
      {
	genjet_1stPt.SetPtEtaPhiE(GenJetPt->at(v_idx_good_genjets[0]),
				  GenJetEta->at(v_idx_good_genjets[0]),
				  GenJetPhi->at(v_idx_good_genjets[0]),
				  GenJetEnergy->at(v_idx_good_genjets[0])
				  );
	genjet_2ndPt.SetPtEtaPhiE(GenJetPt->at(v_idx_good_genjets[1]),
				  GenJetEta->at(v_idx_good_genjets[1]),
				  GenJetPhi->at(v_idx_good_genjets[1]),
				  GenJetEnergy->at(v_idx_good_genjets[1])
				  );
	axigluon_1st2nd_genjets = genjet_1stPt + genjet_2ndPt ;

	//------------

	CreateAndFillUserTH1D("h1_Mjj_genJet_1st2nd", 1000, 0, 2000, axigluon_1st2nd_genjets.M() );	    	
	CreateAndFillUserTH1D("h1_Mjj_genJet_1st2nd_over_massAG", 100, 0, 5,  axigluon_1st2nd_genjets.M() /  v_AG.M()  );	    	

	if( flavour_AG_q != fabs(PDGID_T)  )
	  {
	    CreateAndFillUserTH1D("h1_Mjj_genJet_1st2nd_noTop", 1000, 0, 2000, axigluon_1st2nd_genjets.M() );	    	
	    CreateAndFillUserTH1D("h1_Mjj_genJet_1st2nd_over_massAG_noTop", 100, 0, 5,  axigluon_1st2nd_genjets.M() /  v_AG.M()  );	    	
	  }

	if( flavour_AG_q == fabs(PDGID_T)  )
	  {
	    CreateAndFillUserTH1D("h1_Mjj_genJet_1st2nd_Top", 1000, 0, 2000, axigluon_1st2nd_genjets.M() );	    	
	    CreateAndFillUserTH1D("h1_Mjj_genJet_1st2nd_over_massAG_Top", 100, 0, 5,  axigluon_1st2nd_genjets.M() /  v_AG.M()  );	    	
	  }

      }

    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2)
      {
	recojet_1stPt.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
				   JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
				   JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),
				   JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	recojet_2ndPt.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]),
				   JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]),
				   JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]),
				   JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	axigluon_1st2nd_recojets = recojet_1stPt + recojet_2ndPt;   

	//------------

	CreateAndFillUserTH1D("h1_Mjj_recoJet_1st2nd", 1000, 0, 2000, axigluon_1st2nd_recojets.M() );	    	
	CreateAndFillUserTH1D("h1_Mjj_recoJet_1st2nd_over_massAG", 100, 0, 5,  axigluon_1st2nd_recojets.M() /  v_AG.M()  );	    	


	if( flavour_AG_q != fabs(PDGID_T)  )
	  {
	    CreateAndFillUserTH1D("h1_Mjj_recoJet_1st2nd_noTop", 1000, 0, 2000, axigluon_1st2nd_recojets.M() );	    	
	    CreateAndFillUserTH1D("h1_Mjj_recoJet_1st2nd_over_massAG_noTop", 100, 0, 5,  axigluon_1st2nd_recojets.M() /  v_AG.M()  );	    	
	  }

	if( flavour_AG_q == fabs(PDGID_T)  )
	  {
	    CreateAndFillUserTH1D("h1_Mjj_recoJet_1st2nd_Top", 1000, 0, 2000, axigluon_1st2nd_recojets.M() );	    	
	    CreateAndFillUserTH1D("h1_Mjj_recoJet_1st2nd_over_massAG_Top", 100, 0, 5,  axigluon_1st2nd_recojets.M() /  v_AG.M()  );	    	
	  }

      }


    //3) using fat genjets or fat recojets
    TLorentzVector fatgenjet_1stPt_tmp;
    TLorentzVector fatgenjet_2ndPt_tmp;
    TLorentzVector fatgenjet_1stPt;
    TLorentzVector fatgenjet_2ndPt;
    TLorentzVector axigluon_1st2nd_fatgenjets;

    TLorentzVector fatrecojet_1stPt_tmp;
    TLorentzVector fatrecojet_2ndPt_tmp;
    TLorentzVector fatrecojet_1stPt;
    TLorentzVector fatrecojet_2ndPt;
    TLorentzVector axigluon_1st2nd_fatrecojets;

    double R_FAT = 1.1;

    //genjets
    if( v_idx_good_genjets.size() >= 2)
      {

	fatgenjet_1stPt_tmp = genjet_1stPt;      
	fatgenjet_2ndPt_tmp = genjet_2ndPt;      

	for(int igenjet=0; igenjet < GenJetPt->size(); igenjet++)
	  {
	    
	    TLorentzVector genjet;
	    genjet.SetPtEtaPhiE(GenJetPt->at(igenjet),
				GenJetEta->at(igenjet),
				GenJetPhi->at(igenjet),
				GenJetEnergy->at(igenjet)
				);
	    
	    if( genjet.DeltaR(v_W_l) < 0.3 
		|| genjet.Pt() < 20 
		|| fabs(genjet.Eta()) > 3 
		|| igenjet == v_idx_good_genjets[0] 
		|| igenjet == v_idx_good_genjets[1] 
		)
	      continue;
	    
	    if( genjet.DeltaR(genjet_1stPt) < R_FAT
		&& genjet.DeltaR(genjet_1stPt) < genjet.DeltaR(genjet_2ndPt) )
	      {
		fatgenjet_1stPt_tmp += genjet; 
	      }

	    if( genjet.DeltaR(genjet_2ndPt) < R_FAT
		&& genjet.DeltaR(genjet_2ndPt) < genjet.DeltaR(genjet_1stPt) )
	      {
		fatgenjet_2ndPt_tmp += genjet; 
	      }
	   	    
	  }
	
	//pt ordering
	if( fatgenjet_1stPt_tmp.Pt() > fatgenjet_2ndPt_tmp.Pt() )
	  {
	    fatgenjet_1stPt = fatgenjet_1stPt_tmp;
	    fatgenjet_2ndPt = fatgenjet_2ndPt_tmp;
	  }
	else
	  {
	    fatgenjet_2ndPt = fatgenjet_1stPt_tmp;
	    fatgenjet_1stPt = fatgenjet_2ndPt_tmp;
	  }

	axigluon_1st2nd_fatgenjets = fatgenjet_1stPt + fatgenjet_2ndPt ;

	//------------

	CreateAndFillUserTH1D("h1_Mjj_fatgenJet_1st2nd", 1000, 0, 2000, axigluon_1st2nd_fatgenjets.M() );	    	
	CreateAndFillUserTH1D("h1_Mjj_fatgenJet_1st2nd_over_massAG", 100, 0, 5,  axigluon_1st2nd_fatgenjets.M() /  v_AG.M()  );	    	

	if( flavour_AG_q != fabs(PDGID_T)  )
	  {
	    CreateAndFillUserTH1D("h1_Mjj_fatgenJet_1st2nd_noTop", 1000, 0, 2000, axigluon_1st2nd_fatgenjets.M() );	    	
	    CreateAndFillUserTH1D("h1_Mjj_fatgenJet_1st2nd_over_massAG_noTop", 100, 0, 5,  axigluon_1st2nd_fatgenjets.M() /  v_AG.M()  );	    	
	  }

	if( flavour_AG_q == fabs(PDGID_T)  )
	  {
	    CreateAndFillUserTH1D("h1_Mjj_fatgenJet_1st2nd_Top", 1000, 0, 2000, axigluon_1st2nd_fatgenjets.M() );	    	
	    CreateAndFillUserTH1D("h1_Mjj_fatgenJet_1st2nd_over_massAG_Top", 100, 0, 5,  axigluon_1st2nd_fatgenjets.M() /  v_AG.M()  );	    	
	  }

      }

    //recojets
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2)
      {

	fatrecojet_1stPt_tmp = recojet_1stPt;      
	fatrecojet_2ndPt_tmp = recojet_2ndPt;      

	for(int ijet=0; ijet<v_idx_jet_PtCut_noOverlap_ID.size(); ijet++)
	  {

	    TLorentzVector recojet;
	    recojet.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
			     JetEta->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
			     JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
			     JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[ijet]) );

	    if( ijet == v_idx_jet_PtCut_noOverlap_ID[0] 
		|| ijet == v_idx_jet_PtCut_noOverlap_ID[1] 
		)
	      continue;
	
	    if( recojet.DeltaR(recojet_1stPt) < R_FAT
		&& recojet.DeltaR(recojet_1stPt) < recojet.DeltaR(recojet_2ndPt) )
	      {
		fatrecojet_1stPt_tmp += recojet; 
	      }
	
	    if( recojet.DeltaR(recojet_2ndPt) < R_FAT
		&& recojet.DeltaR(recojet_2ndPt) < recojet.DeltaR(recojet_1stPt) )
	      {
		fatrecojet_2ndPt_tmp += recojet; 
	      }

	  }

	//pt ordering
	if( fatrecojet_1stPt_tmp.Pt() > fatrecojet_2ndPt_tmp.Pt() )
	  {
	    fatrecojet_1stPt = fatrecojet_1stPt_tmp;
	    fatrecojet_2ndPt = fatrecojet_2ndPt_tmp;
	  }
	else
	  {
	    fatrecojet_2ndPt = fatrecojet_1stPt_tmp;
	    fatrecojet_1stPt = fatrecojet_2ndPt_tmp;
	  }
    
	axigluon_1st2nd_fatrecojets = fatrecojet_1stPt + fatrecojet_2ndPt ;
    
	//------------

	CreateAndFillUserTH1D("h1_Mjj_fatrecoJet_1st2nd", 1000, 0, 2000, axigluon_1st2nd_fatrecojets.M() );	    	
	CreateAndFillUserTH1D("h1_Mjj_fatrecoJet_1st2nd_over_massAG", 100, 0, 5,  axigluon_1st2nd_fatrecojets.M() /  v_AG.M()  );	    	

	if( flavour_AG_q != fabs(PDGID_T)  )
	  {
	    CreateAndFillUserTH1D("h1_Mjj_fatrecoJet_1st2nd_noTop", 1000, 0, 2000, axigluon_1st2nd_fatrecojets.M() );	    	
	    CreateAndFillUserTH1D("h1_Mjj_fatrecoJet_1st2nd_over_massAG_noTop", 100, 0, 5,  axigluon_1st2nd_fatrecojets.M() /  v_AG.M()  );	    	
	  }

	if( flavour_AG_q == fabs(PDGID_T)  )
	  {
	    CreateAndFillUserTH1D("h1_Mjj_fatrecoJet_1st2nd_Top", 1000, 0, 2000, axigluon_1st2nd_fatrecojets.M() );	    	
	    CreateAndFillUserTH1D("h1_Mjj_fatrecoJet_1st2nd_over_massAG_Top", 100, 0, 5,  axigluon_1st2nd_fatrecojets.M() /  v_AG.M()  );	    	
	  }
      
      }

    //INFO
    //      // retrieve value of previously filled variables (after making sure that they were filled)
    //      double totpTEle;
    //      if ( variableIsFilled("pT1stEle") && variableIsFilled("pT2ndEle") )
    //        totpTEle = getVariableValue("pT1stEle")+getVariableValue("pT2ndEle");
    //      // reject events that did not pass level 0 cuts
    //      if( !passedCut("0") ) continue;
    //      // ......

    //after pre-selection
    //    if( passedAllPreviousCuts("Pt1stEle_PAS")
    //	&& variableIsFilled("MTenu_PAS") && variableIsFilled("sT_PAS")
    //FillUserTH1D("h1_MTenu_PAS_plus", getVariableValue("MTenu_PAS"));

    // Produce skim
    //if( passedAllPreviousCuts("minDRej") ) fillSkimTree();

    ////////////////////// User's code to be done for every event - END ///////////////////////

  } // End of loop over events


  ////////////////////// User's code to write histos - BEGIN ///////////////////////

  h1_DecayProducts_AG->Write();
  h1_DecayProducts_W->Write();

  ////////////////////// User's code to write histos - END ///////////////////////

  delete h1_DecayProducts_AG;
  
  //STDOUT("analysisClass::Loop() ends");
}
