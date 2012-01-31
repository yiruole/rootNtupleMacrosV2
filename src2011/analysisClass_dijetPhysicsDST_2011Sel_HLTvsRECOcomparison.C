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
#include "TH2.h"
#include "TMath.h"
#include "TProfile.h"

//For JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

//-----------------------------
//-----------------------------

Double_t median1(TH1 *h1) { 
   //compute the median for 1-d histogram h1 
   Int_t nbins = h1->GetXaxis()->GetNbins(); 
   Double_t *x = new Double_t[nbins]; 
   Double_t *y = new Double_t[nbins]; 
   for (Int_t i=0;i<nbins;i++) {
      x[i] = h1->GetXaxis()->GetBinCenter(i+1); 
      y[i] = h1->GetBinContent(i+1); 
   } 
   Double_t median = TMath::Median(nbins,x,y); 
   delete [] x; 
   delete [] y; 
   return median; 
} 

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

  //--------------------------------------------------------------------------
  //For JEC
  //--------------------------------------------------------------------------

  string L2Tag     = "/afs/cern.ch/user/s/santanas/public/JEC_txtfiles_GR_R_42_V19/GR_R_42_V19_AK5PF_L2Relative.txt"; 
  string L3Tag     = "/afs/cern.ch/user/s/santanas/public/JEC_txtfiles_GR_R_42_V19/GR_R_42_V19_AK5PF_L3Absolute.txt";
  string L23ResTag = "/afs/cern.ch/user/s/santanas/public/JEC_txtfiles_GR_R_42_V19/GR_R_42_V19_AK5PF_L2L3Residual.txt";
  
  JetCorrectorParameters *L2Par = new JetCorrectorParameters(L2Tag);
  JetCorrectorParameters *L3Par = new JetCorrectorParameters(L3Tag);
  JetCorrectorParameters *L23ResPar = new JetCorrectorParameters(L23ResTag);

  vector<JetCorrectorParameters> vParL2;
  vParL2.push_back(*L2Par);
  vector<JetCorrectorParameters> vParL3;
  vParL3.push_back(*L3Par);
  vector<JetCorrectorParameters> vParL23Res;
  vParL23Res.push_back(*L23ResPar);
  vector<JetCorrectorParameters> vParL2L3;
  vParL2L3.push_back(*L2Par);
  vParL2L3.push_back(*L3Par);
  vector<JetCorrectorParameters> vParAll;
  vParAll.push_back(*L2Par);
  vParAll.push_back(*L3Par);
  vParAll.push_back(*L23ResPar);

  FactorizedJetCorrector *JetCorrectorL2 = new FactorizedJetCorrector(vParL2);
  FactorizedJetCorrector *JetCorrectorL3 = new FactorizedJetCorrector(vParL3);
  FactorizedJetCorrector *JetCorrectorL23Res = new FactorizedJetCorrector(vParL23Res);
  FactorizedJetCorrector *JetCorrectorL2L3 = new FactorizedJetCorrector(vParL2L3);
  FactorizedJetCorrector *JetCorrectorAll = new FactorizedJetCorrector(vParAll);

  //   JetCorrectorL2->setJetEta(0.5);
  //   JetCorrectorL3->setJetEta(0.5);
  //   JetCorrectorL23Res->setJetEta(0.5);
  //   JetCorrectorAll->setJetEta(0.5);
  //   JetCorrectorL2->setJetPt(100);
  //   JetCorrectorL3->setJetPt(100);
  //   JetCorrectorL23Res->setJetPt(100);
  //   JetCorrectorAll->setJetPt(100);

  //   double corrL2 = JetCorrectorL2->getCorrection();
  //   double corrL3 = JetCorrectorL3->getCorrection();
  //   double corrL23Res = JetCorrectorL23Res->getCorrection();
  //   double corrAll = JetCorrectorAll->getCorrection();
  //   cout << "L2, L3, L23Res, Product :  " << corrL2 << ", " <<corrL3 << ", " << corrL23Res << ", " << corrL2*corrL3*corrL23Res << endl;
  //   cout << "All : " << corrAll << endl;
  
  //--------------------------------------------------------------------------

  //STDOUT("analysisClass::Loop() begins");

  if (fChain == 0) return;

   //--------------------------------------------------------------------------
   // Decide which plots to save (default is to save everything)
   //--------------------------------------------------------------------------
   
   fillSkim                         (  true  ) ;
   fillAllPreviousCuts              ( !true  ) ;
   fillAllOtherCuts                 ( !true  ) ;
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;

  ////////////////////// User's code to get preCut values - BEGIN ///////////////

  // Jets
  double pfjetPtCut = getPreCutValue1("pfjetPtCut");
  double pfjetEtaCut = getPreCutValue1("pfjetEtaCut");
  double pfjetIDloose = getPreCutValue1("pfjetIDloose");
  double pfjetIDtight = getPreCutValue1("pfjetIDtight");
  double pfjetPtCutForHT = getPreCutValue2("pfjetPtCut");
  double pfjetEtaCutForHT = getPreCutValue2("pfjetEtaCut");
  double pfjetDEtaForFatJet = getPreCutValue1("pfjetDEtaForFatJet");
  double pfjetRForFatJet = getPreCutValue1("pfjetRForFatJet");

  // Vertices
  double vertexMaxAbsZ = getPreCutValue1("vertexMaxAbsZ");
  double vertexMaxd0 = getPreCutValue1("vertexMaxd0");  

  //Others
  double applyPFJEC = getPreCutValue1("applyPFJEC");  
  double ptUE_DATA = getPreCutValue1("ptUE");  
  double ptUE_MC   = getPreCutValue2("ptUE");  
  double debug     = getPreCutValue1("debug");  
  double debugEv     = getPreCutValue2("debug");  

  ////////////////////// User's code to get preCut values - END /////////////////

  ////////////////////// User's code to book histos - BEGIN ///////////////////////

  CreateUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt10GeV_EtaL2p4"             ,    200, 0, 10      );
  CreateUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4"             ,    100, -1, 1      );
  CreateUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4"        ,    100, -1, 1      );
  CreateUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt30GeV_EtaL2p4"             ,    200, 0, 10      );
  CreateUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4"             ,    100, -1, 1      );
  CreateUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4"        ,    100, -1, 1      );
  CreateUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt50GeV_EtaL2p4"             ,    200, 0, 10      );
  CreateUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4"             ,    100, -1, 1      );
  CreateUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4"        ,    100, -1, 1      );
  CreateUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt100GeV_EtaL2p4"            ,    200, 0, 10      );
  CreateUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4"            ,    100, -1, 1      );
  CreateUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4"       ,    100, -1, 1      );

  CreateUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3"             ,    200, 0, 10      );
  CreateUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3"             ,    100, -1, 1      );
  CreateUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3"        ,    100, -1, 1      );
  CreateUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3"             ,    200, 0, 10      );
  CreateUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3"             ,    100, -1, 1      );
  CreateUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3"        ,    100, -1, 1      );
  CreateUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3"             ,    200, 0, 10      );
  CreateUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3"             ,    100, -1, 1      );
  CreateUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3"        ,    100, -1, 1      );
  CreateUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3"            ,    200, 0, 10      );
  CreateUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3"            ,    100, -1, 1      );
  CreateUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3"       ,    100, -1, 1      );

  CreateUserTH1D( "h_HLTPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"    , 100,  0, 1   );
  CreateUserTH1D( "h_HLTPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"        , 100,  0, 1   );
  CreateUserTH1D( "h_HLTPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4"                  , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_HLTPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"    , 100,  0, 1   );
  CreateUserTH1D( "h_HLTPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4"            , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_HLTPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"        , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"    , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"        , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4"                  , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_RECOPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"    , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4"            , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_RECOPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"        , 100,  0, 1   );

  CreateUserTH1D( "h_HLTPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"    , 100,  0, 1   );
  CreateUserTH1D( "h_HLTPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"        , 100,  0, 1   );
  CreateUserTH1D( "h_HLTPFJetNConstituents_smallestDR_Pt30GeV_Eta2p4_3"                  , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_HLTPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"    , 100,  0, 1   );
  CreateUserTH1D( "h_HLTPFJetChargedMultiplicity_smallestDR_Pt30GeV_Eta2p4_3"            , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_HLTPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"        , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"    , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"        , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetNConstituents_smallestDR_Pt30GeV_Eta2p4_3"                  , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_RECOPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"    , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetChargedMultiplicity_smallestDR_Pt30GeV_Eta2p4_3"            , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_RECOPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"        , 100,  0, 1   );

  CreateUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt10GeV_EtaL2p4"                    , 2, -0.5, 1.5, 2, -0.5, 1.5 );
  CreateUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt30GeV_EtaL2p4"                    , 2, -0.5, 1.5, 2, -0.5, 1.5 );
  CreateUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt50GeV_EtaL2p4"                    , 2, -0.5, 1.5, 2, -0.5, 1.5 );
  CreateUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt100GeV_EtaL2p4"                   , 2, -0.5, 1.5, 2, -0.5, 1.5 );
  CreateUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt10GeV_EtaL2p4"                    , 2, -0.5, 1.5, 2, -0.5, 1.5 );
  CreateUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt30GeV_EtaL2p4"                    , 2, -0.5, 1.5, 2, -0.5, 1.5 );
  CreateUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt50GeV_EtaL2p4"                    , 2, -0.5, 1.5, 2, -0.5, 1.5 );
  CreateUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt100GeV_EtaL2p4"                   , 2, -0.5, 1.5, 2, -0.5, 1.5 );

  CreateUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3"                   , 2, -0.5, 1.5, 2, -0.5, 1.5 );
  CreateUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3"                   , 2, -0.5, 1.5, 2, -0.5, 1.5 );
  CreateUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3"                   , 2, -0.5, 1.5, 2, -0.5, 1.5 );
  CreateUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3"                  , 2, -0.5, 1.5, 2, -0.5, 1.5 );
  CreateUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3"                   , 2, -0.5, 1.5, 2, -0.5, 1.5 );
  CreateUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3"                   , 2, -0.5, 1.5, 2, -0.5, 1.5 );
  CreateUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3"                   , 2, -0.5, 1.5, 2, -0.5, 1.5 );
  CreateUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3"                  , 2, -0.5, 1.5, 2, -0.5, 1.5 );

  CreateUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"             ,    200, 0, 10      );
  CreateUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"             ,    100, -1, 1      );
  CreateUserTH1D( "h_HLTPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"    , 100,  0, 1   );
  CreateUserTH1D( "h_HLTPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"        , 100,  0, 1   );
  CreateUserTH1D( "h_HLTPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"                  , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_HLTPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"    , 100,  0, 1   );
  CreateUserTH1D( "h_HLTPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"            , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_HLTPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"        , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"    , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"        , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"                  , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_RECOPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"    , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"            , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_RECOPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"        , 100,  0, 1   );

  CreateUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"             ,    200, 0, 10      );
  CreateUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"             ,    100, -1, 1      );
  CreateUserTH1D( "h_HLTPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"    , 100,  0, 1   );
  CreateUserTH1D( "h_HLTPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"        , 100,  0, 1   );
  CreateUserTH1D( "h_HLTPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"                  , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_HLTPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"    , 100,  0, 1   );
  CreateUserTH1D( "h_HLTPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"            , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_HLTPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"        , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"    , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"        , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"                  , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_RECOPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"    , 100,  0, 1   );
  CreateUserTH1D( "h_RECOPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"            , 201,  -0.5, 200.5   );
  CreateUserTH1D( "h_RECOPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"        , 100,  0, 1   );

  CreateUserTH1D( "h_HLT_PFJet1_Pt", getHistoNBins("HLT_PFJet1_Pt"), getHistoMin("HLT_PFJet1_Pt"), getHistoMax("HLT_PFJet1_Pt")  ) ; 
  CreateUserTH1D( "h_HLT_PFJet1_Eta", getHistoNBins("HLT_PFJet1_Eta"), getHistoMin("HLT_PFJet1_Eta"), getHistoMax("HLT_PFJet1_Eta")  ) ; 
  CreateUserTH1D( "h_HLT_PFJet2_Pt", getHistoNBins("HLT_PFJet2_Pt"), getHistoMin("HLT_PFJet2_Pt"), getHistoMax("HLT_PFJet2_Pt")  ) ; 
  CreateUserTH1D( "h_HLT_PFJet2_Eta", getHistoNBins("HLT_PFJet2_Eta"), getHistoMin("HLT_PFJet2_Eta"), getHistoMax("HLT_PFJet2_Eta")  ) ; 
  CreateUserTH1D( "h_HLT_DR_PFJet1PFJet2", getHistoNBins("HLT_DR_PFJet1PFJet2"), getHistoMin("HLT_DR_PFJet1PFJet2"), getHistoMax("HLT_DR_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_HLT_DPhi_PFJet1PFJet2", getHistoNBins("HLT_DPhi_PFJet1PFJet2"), getHistoMin("HLT_DPhi_PFJet1PFJet2"), getHistoMax("HLT_DPhi_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_HLT_DEta_PFJet1PFJet2", getHistoNBins("HLT_DEta_PFJet1PFJet2"), getHistoMin("HLT_DEta_PFJet1PFJet2"), getHistoMax("HLT_DEta_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_HLT_Pt_PFJet1PFJet2", getHistoNBins("HLT_Pt_PFJet1PFJet2"), getHistoMin("HLT_Pt_PFJet1PFJet2"), getHistoMax("HLT_Pt_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_HLT_M_PFJet1PFJet2", getHistoNBins("HLT_M_PFJet1PFJet2"), getHistoMin("HLT_M_PFJet1PFJet2"), getHistoMax("HLT_M_PFJet1PFJet2")  ) ; 

  CreateUserTH1D( "h_RECO_PFJet1_Pt", getHistoNBins("RECO_PFJet1_Pt"), getHistoMin("RECO_PFJet1_Pt"), getHistoMax("RECO_PFJet1_Pt")  ) ; 
  CreateUserTH1D( "h_RECO_PFJet1_Eta", getHistoNBins("RECO_PFJet1_Eta"), getHistoMin("RECO_PFJet1_Eta"), getHistoMax("RECO_PFJet1_Eta")  ) ; 
  CreateUserTH1D( "h_RECO_PFJet2_Pt", getHistoNBins("RECO_PFJet2_Pt"), getHistoMin("RECO_PFJet2_Pt"), getHistoMax("RECO_PFJet2_Pt")  ) ; 
  CreateUserTH1D( "h_RECO_PFJet2_Eta", getHistoNBins("RECO_PFJet2_Eta"), getHistoMin("RECO_PFJet2_Eta"), getHistoMax("RECO_PFJet2_Eta")  ) ; 
  CreateUserTH1D( "h_RECO_DR_PFJet1PFJet2", getHistoNBins("RECO_DR_PFJet1PFJet2"), getHistoMin("RECO_DR_PFJet1PFJet2"), getHistoMax("RECO_DR_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_RECO_DPhi_PFJet1PFJet2", getHistoNBins("RECO_DPhi_PFJet1PFJet2"), getHistoMin("RECO_DPhi_PFJet1PFJet2"), getHistoMax("RECO_DPhi_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_RECO_DEta_PFJet1PFJet2", getHistoNBins("RECO_DEta_PFJet1PFJet2"), getHistoMin("RECO_DEta_PFJet1PFJet2"), getHistoMax("RECO_DEta_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_RECO_Pt_PFJet1PFJet2", getHistoNBins("RECO_Pt_PFJet1PFJet2"), getHistoMin("RECO_Pt_PFJet1PFJet2"), getHistoMax("RECO_Pt_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_RECO_M_PFJet1PFJet2", getHistoNBins("RECO_M_PFJet1PFJet2"), getHistoMin("RECO_M_PFJet1PFJet2"), getHistoMax("RECO_M_PFJet1PFJet2")  ) ; 

  CreateUserTH1D( "h_HLT_PFJet1_Pt_MATCH", getHistoNBins("HLT_PFJet1_Pt"), getHistoMin("HLT_PFJet1_Pt"), getHistoMax("HLT_PFJet1_Pt")  ) ; 
  CreateUserTH1D( "h_HLT_PFJet1_Eta_MATCH", getHistoNBins("HLT_PFJet1_Eta"), getHistoMin("HLT_PFJet1_Eta"), getHistoMax("HLT_PFJet1_Eta")  ) ; 
  CreateUserTH1D( "h_HLT_PFJet2_Pt_MATCH", getHistoNBins("HLT_PFJet2_Pt"), getHistoMin("HLT_PFJet2_Pt"), getHistoMax("HLT_PFJet2_Pt")  ) ; 
  CreateUserTH1D( "h_HLT_PFJet2_Eta_MATCH", getHistoNBins("HLT_PFJet2_Eta"), getHistoMin("HLT_PFJet2_Eta"), getHistoMax("HLT_PFJet2_Eta")  ) ; 
  CreateUserTH1D( "h_HLT_DR_PFJet1PFJet2_MATCH", getHistoNBins("HLT_DR_PFJet1PFJet2"), getHistoMin("HLT_DR_PFJet1PFJet2"), getHistoMax("HLT_DR_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_HLT_DPhi_PFJet1PFJet2_MATCH", getHistoNBins("HLT_DPhi_PFJet1PFJet2"), getHistoMin("HLT_DPhi_PFJet1PFJet2"), getHistoMax("HLT_DPhi_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_HLT_DEta_PFJet1PFJet2_MATCH", getHistoNBins("HLT_DEta_PFJet1PFJet2"), getHistoMin("HLT_DEta_PFJet1PFJet2"), getHistoMax("HLT_DEta_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_HLT_Pt_PFJet1PFJet2_MATCH", getHistoNBins("HLT_Pt_PFJet1PFJet2"), getHistoMin("HLT_Pt_PFJet1PFJet2"), getHistoMax("HLT_Pt_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_HLT_M_PFJet1PFJet2_MATCH", getHistoNBins("HLT_M_PFJet1PFJet2"), getHistoMin("HLT_M_PFJet1PFJet2"), getHistoMax("HLT_M_PFJet1PFJet2")  ) ; 

  CreateUserTH1D( "h_RECO_PFJet1_Pt_MATCH", getHistoNBins("RECO_PFJet1_Pt"), getHistoMin("RECO_PFJet1_Pt"), getHistoMax("RECO_PFJet1_Pt")  ) ; 
  CreateUserTH1D( "h_RECO_PFJet1_Eta_MATCH", getHistoNBins("RECO_PFJet1_Eta"), getHistoMin("RECO_PFJet1_Eta"), getHistoMax("RECO_PFJet1_Eta")  ) ; 
  CreateUserTH1D( "h_RECO_PFJet2_Pt_MATCH", getHistoNBins("RECO_PFJet2_Pt"), getHistoMin("RECO_PFJet2_Pt"), getHistoMax("RECO_PFJet2_Pt")  ) ; 
  CreateUserTH1D( "h_RECO_PFJet2_Eta_MATCH", getHistoNBins("RECO_PFJet2_Eta"), getHistoMin("RECO_PFJet2_Eta"), getHistoMax("RECO_PFJet2_Eta")  ) ; 
  CreateUserTH1D( "h_RECO_DR_PFJet1PFJet2_MATCH", getHistoNBins("RECO_DR_PFJet1PFJet2"), getHistoMin("RECO_DR_PFJet1PFJet2"), getHistoMax("RECO_DR_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_RECO_DPhi_PFJet1PFJet2_MATCH", getHistoNBins("RECO_DPhi_PFJet1PFJet2"), getHistoMin("RECO_DPhi_PFJet1PFJet2"), getHistoMax("RECO_DPhi_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_RECO_DEta_PFJet1PFJet2_MATCH", getHistoNBins("RECO_DEta_PFJet1PFJet2"), getHistoMin("RECO_DEta_PFJet1PFJet2"), getHistoMax("RECO_DEta_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_RECO_Pt_PFJet1PFJet2_MATCH", getHistoNBins("RECO_Pt_PFJet1PFJet2"), getHistoMin("RECO_Pt_PFJet1PFJet2"), getHistoMax("RECO_Pt_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_RECO_M_PFJet1PFJet2_MATCH", getHistoNBins("RECO_M_PFJet1PFJet2"), getHistoMin("RECO_M_PFJet1PFJet2"), getHistoMax("RECO_M_PFJet1PFJet2")  ) ; 

  CreateUserTH1D( "h_MassBias_MATCH"  ,  100, -1, 1      );
  CreateUserTH1D( "h_DEtaBias_MATCH"  ,  100, -1, 1      );
  CreateUserTH1D( "h_HLT_M_PFJet1PFJet2_MATCH_MassBiasMore0p2", getHistoNBins("HLT_M_PFJet1PFJet2"), getHistoMin("HLT_M_PFJet1PFJet2"), getHistoMax("HLT_M_PFJet1PFJet2")  ) ; 
  CreateUserTH1D( "h_RECO_M_PFJet1PFJet2_MATCH_MassBiasMore0p2", getHistoNBins("RECO_M_PFJet1PFJet2"), getHistoMin("RECO_M_PFJet1PFJet2"), getHistoMax("RECO_M_PFJet1PFJet2")  ) ; 
  TProfile* p_MassBias_MATCH_vs_Mreco = new TProfile( "p_MassBias_MATCH_vs_Mreco" , "p_MassBias_MATCH_vs_Mreco", 
						      40, 0, 2000, 
						      -50, 50);

  CreateUserTH2D( "h2_FinalSelectionAgreement_recoY_vs_hltX"                                                , 2, -0.5, 1.5, 2, -0.5, 1.5 );

  TH1D *h_FirstMismatch_HLT1_RECO0 = new TH1D ("h_FirstMismatch_HLT1_RECO0","h_FirstMismatch_HLT1_RECO0",9,-0.5,8.5);
  h_FirstMismatch_HLT1_RECO0->Sumw2();
  h_FirstMismatch_HLT1_RECO0->GetXaxis()->SetBinLabel(1,"Ev. Filter");
  h_FirstMismatch_HLT1_RECO0->GetXaxis()->SetBinLabel(2,"Jet1Pt");
  h_FirstMismatch_HLT1_RECO0->GetXaxis()->SetBinLabel(3,"Jet1Eta");
  h_FirstMismatch_HLT1_RECO0->GetXaxis()->SetBinLabel(4,"Jet1ID");
  h_FirstMismatch_HLT1_RECO0->GetXaxis()->SetBinLabel(5,"Jet2Pt");
  h_FirstMismatch_HLT1_RECO0->GetXaxis()->SetBinLabel(6,"Jet2Eta");
  h_FirstMismatch_HLT1_RECO0->GetXaxis()->SetBinLabel(7,"Jet2ID");
  h_FirstMismatch_HLT1_RECO0->GetXaxis()->SetBinLabel(8,"DPhiJJ");
  h_FirstMismatch_HLT1_RECO0->GetXaxis()->SetBinLabel(9,"DEtaJJ");

  CreateUserTH1D("h_HLT_M_PFJet1PFJet2_HLT1_RECO0", getHistoNBins("HLT_M_PFJet1PFJet2"), getHistoMin("HLT_M_PFJet1PFJet2"), getHistoMax("HLT_M_PFJet1PFJet2")  ) ; 
  CreateUserTH1D("h_RECO_M_PFJet1PFJet2_HLT1_RECO0",  getHistoNBins("RECO_M_PFJet1PFJet2"), getHistoMin("RECO_M_PFJet1PFJet2"), getHistoMax("RECO_M_PFJet1PFJet2")  ) ; 
  CreateUserTH1D("h_HLT_DPhi_PFJet1PFJet2_HLT1_RECO0", getHistoNBins("HLT_DPhi_PFJet1PFJet2"), getHistoMin("HLT_DPhi_PFJet1PFJet2"), getHistoMax("HLT_DPhi_PFJet1PFJet2")  ) ; 
  CreateUserTH1D("h_RECO_DPhi_PFJet1PFJet2_HLT1_RECO0",  getHistoNBins("RECO_DPhi_PFJet1PFJet2"), getHistoMin("RECO_DPhi_PFJet1PFJet2"), getHistoMax("RECO_DPhi_PFJet1PFJet2")  ) ; 

  TH1D *h_FirstMismatch_HLT0_RECO1 = new TH1D ("h_FirstMismatch_HLT0_RECO1","h_FirstMismatch_HLT0_RECO1",9,-0.5,8.5);
  h_FirstMismatch_HLT0_RECO1->Sumw2();
  h_FirstMismatch_HLT0_RECO1->GetXaxis()->SetBinLabel(1,"Ev. Filter");
  h_FirstMismatch_HLT0_RECO1->GetXaxis()->SetBinLabel(2,"Jet1Pt");
  h_FirstMismatch_HLT0_RECO1->GetXaxis()->SetBinLabel(3,"Jet1Eta");
  h_FirstMismatch_HLT0_RECO1->GetXaxis()->SetBinLabel(4,"Jet1ID");
  h_FirstMismatch_HLT0_RECO1->GetXaxis()->SetBinLabel(5,"Jet2Pt");
  h_FirstMismatch_HLT0_RECO1->GetXaxis()->SetBinLabel(6,"Jet2Eta");
  h_FirstMismatch_HLT0_RECO1->GetXaxis()->SetBinLabel(7,"Jet2ID");
  h_FirstMismatch_HLT0_RECO1->GetXaxis()->SetBinLabel(8,"DPhiJJ");
  h_FirstMismatch_HLT0_RECO1->GetXaxis()->SetBinLabel(9,"DEtaJJ");

  CreateUserTH1D("h_HLT_M_PFJet1PFJet2_HLT0_RECO1", getHistoNBins("HLT_M_PFJet1PFJet2"), getHistoMin("HLT_M_PFJet1PFJet2"), getHistoMax("HLT_M_PFJet1PFJet2")  ) ; 
  CreateUserTH1D("h_RECO_M_PFJet1PFJet2_HLT0_RECO1",  getHistoNBins("RECO_M_PFJet1PFJet2"), getHistoMin("RECO_M_PFJet1PFJet2"), getHistoMax("RECO_M_PFJet1PFJet2")  ) ; 
  CreateUserTH1D("h_HLT_DPhi_PFJet1PFJet2_HLT0_RECO1", getHistoNBins("HLT_DPhi_PFJet1PFJet2"), getHistoMin("HLT_DPhi_PFJet1PFJet2"), getHistoMax("HLT_DPhi_PFJet1PFJet2")  ) ; 
  CreateUserTH1D("h_RECO_DPhi_PFJet1PFJet2_HLT0_RECO1",  getHistoNBins("RECO_DPhi_PFJet1PFJet2"), getHistoMin("RECO_DPhi_PFJet1PFJet2"), getHistoMax("RECO_DPhi_PFJet1PFJet2")  ) ; 

  //For Pileup studies
  TProfile* p_PTJetsMedian_vs_NVtx = new TProfile( "p_PTJetsMedian_vs_NVtx" , "p_PTJetsMedian_vs_NVtx", 51, -0.5, 50.5, 0, 200);
  TH1D *h_PtPFJetsInEvent1Ev = new TH1D ("h_PtPFJetsInEvent1Ev","h_PtPFJetsInEvent1Ev",10000,0,1000);
  h_PtPFJetsInEvent1Ev->Sumw2();
  float tobedone = 1;
  
  ////////////////////// User's code to book histos - END ///////////////////////

  Long64_t nentries = fChain->GetEntriesFast();
  STDOUT("analysisClass::Loop(): nentries = " << nentries);

  ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
  ////// If the root version is updated and rootNtupleClass regenerated,     /////
  ////// these lines may need to be updated.                                 /////
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) { // Begin of loop over events
  //for (Long64_t jentry=0; jentry<2000;jentry++) { // Begin of loop over events
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry < 10 || jentry%1000 == 0) STDOUT("analysisClass::Loop(): jentry = " << jentry);
    // if (Cut(ientry) < 0) continue;

    //int NPILEUP_AVE = int( (PileUpInteractions->at(0) + PileUpInteractions->at(1) + PileUpInteractions->at(2))/3 );
    //int NPILEUP_FINAL = min( NPILEUP_AVE , 25 );
    //double event_weight = getPileupWeight ( NPILEUP_FINAL, isData ) ;
    //double event_weight = getPileupWeight ( min(PileUpInteractions->at(1),25), isData ) ;

    ////////////////////// User's code to be done for every event - BEGIN ///////////////////////

    //## HLT PF jets with JEC
    std::auto_ptr<std::vector<double> >  HLTPFJetCorrPt  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  HLTPFJetCorrEnergy  ( new std::vector<double>()  );
    for(int ijet=0; ijet<HLTPFJetPt->size(); ijet++)
      {
	HLTPFJetCorrPt->push_back( HLTPFJetPt->at(ijet) );		
	HLTPFJetCorrEnergy->push_back( HLTPFJetEnergy->at(ijet) );		
      }

    //## Get pileup contribution for each event
    double pt_PU_UE = 0.;
    TH1D *h_PtPFJetsInEventTmp = new TH1D ("h_PtPFJetsInEventTmp","h_PtPFJetsInEventTmp",10000,0,1000);
    h_PtPFJetsInEventTmp->Sumw2();
    //Loop over PF jets
    for(int ijet=0; ijet<HLTPFJetPt->size(); ijet++)
      {
	//eta pre-cuts on jets
	if ( fabs( HLTPFJetEta->at(ijet) ) > pfjetEtaCut ) continue;
	
	h_PtPFJetsInEventTmp->Fill( HLTPFJetPt->at(ijet) );	    
      }
    pt_PU_UE = median1(h_PtPFJetsInEventTmp);
    delete h_PtPFJetsInEventTmp;
             
    //## Apply JEC to PF jets    
    if(applyPFJEC)
      {
	for(int ijet=0; ijet<HLTPFJetPt->size(); ijet++)
	  {	    
	    if(isData)  //Data: L1,L2,L3,L23Res
	      {
		//L1 FastJet-like correction
		double L1 = 1 - ( (pt_PU_UE - ptUE_DATA) / HLTPFJetPt->at(ijet) ) ;
		if(L1>1 || L1<=0) //pathological cases
		  L1 = 1;//don't do anything

		JetCorrectorAll->setJetEta( HLTPFJetEta->at(ijet) );
		JetCorrectorAll->setJetPt( L1 * HLTPFJetPt->at(ijet) );
		double thiscorrection = JetCorrectorAll->getCorrection();

		HLTPFJetCorrPt->at(ijet) *= L1 * thiscorrection;
		HLTPFJetCorrEnergy->at(ijet) *= L1 * thiscorrection;
	      }
	    else        //MC: L1,L2,L3 only (no residual)  
	      {
		//L1 FastJet-like correction
		double L1 = 1 - ( (pt_PU_UE - ptUE_MC) / HLTPFJetPt->at(ijet) ) ;
		if(L1>1 || L1<=0) //pathological cases
		  L1 = 1;//don't do anything

		JetCorrectorL2L3->setJetEta( HLTPFJetEta->at(ijet) );
		JetCorrectorL2L3->setJetPt( L1 * HLTPFJetPt->at(ijet) );
		double thiscorrection = JetCorrectorL2L3->getCorrection();

		HLTPFJetCorrPt->at(ijet) *= L1 * thiscorrection;
		HLTPFJetCorrEnergy->at(ijet) *= L1 * thiscorrection;
	      }
	  }
      }

    //#Get trigger information, if necessary                     
    //    if ( isData ) {
    vector<int>    *fakeHLTTriggerPrescales;
    fakeHLTTriggerPrescales = new std::vector<int> ( int (HLTTriggerNames -> size()), 1 );
    for (int ihlt=0 ; ihlt< HLTTriggerNames->size() ; ihlt++)
      {
	fakeHLTTriggerPrescales->push_back( 1 );
      }
    getTriggers ( HLTKey, HLTTriggerNames, HLTTriggerDecisions, fakeHLTTriggerPrescales ) ;
    delete fakeHLTTriggerPrescales;
    //    }
    //printTriggers();
    

    //#RECO vs HLT comparison jet-by-jet
    // Loop over RECO jets 
    for(int irecojet=0; irecojet<PFJetPtRaw->size(); irecojet++)
      {	
	TLorentzVector recopfjetraw;
	recopfjetraw.SetPtEtaPhiE(PFJetPtRaw->at(irecojet),
				  PFJetEta->at(irecojet),
				  PFJetPhi->at(irecojet),
				  PFJetEnergyRaw->at(irecojet) );
	
	double min_DeltaR_reco_hlt = 99999;
	double idx_hlt_min_DeltaR_fromreco = -1;

	//Reco JetID

	// PFJet ID ( https://twiki.cern.ch/twiki/bin/view/CMS/JetID#Documentation ) 
	// --> only defined for jet |eta| < 3
	
	bool passLooseID_recojet = true;
	bool passTightID_recojet = true;
	
	if( fabs( PFJetEta->at(irecojet) )<3 
	    && 
	    (
	     PFJetNeutralHadronEnergyFraction->at(irecojet) >= 0.99
	     || PFJetNeutralEmEnergyFraction->at(irecojet) >= 0.99
	     || PFJetNConstituents->at(irecojet) <= 1
	     )
	    )
	  {
	    passLooseID_recojet = false;
	  }
	
	if( fabs( PFJetEta->at(irecojet) )<3 
	    && 
	    (
	     PFJetNeutralHadronEnergyFraction->at(irecojet) >= 0.90
	     || PFJetNeutralEmEnergyFraction->at(irecojet) >= 0.90
	     || PFJetNConstituents->at(irecojet) <= 1
	     )
	    )
	  {
	    passTightID_recojet = false;
	  }
	
	if( fabs( PFJetEta->at(irecojet) )<2.4
	    &&
	    (
	     PFJetChargedHadronEnergyFraction->at(irecojet) <=0
	     || PFJetChargedMultiplicity->at(irecojet) <=0
	     || PFJetChargedEmEnergyFraction->at(irecojet) >=0.99
	     )
	    )
	  {
	    passLooseID_recojet = false;
	    passTightID_recojet = false;
	  }

	// Loop over HLT jets 
	for(int ihltjet=0; ihltjet<HLTPFJetPt->size(); ihltjet++)
	  {	  
	    TLorentzVector hltpfjetraw;
	    hltpfjetraw.SetPtEtaPhiE(HLTPFJetPt->at(ihltjet),
				     HLTPFJetEta->at(ihltjet),
				     HLTPFJetPhi->at(ihltjet),
				     HLTPFJetEnergy->at(ihltjet) );  	    

	    //-----------------------------------------------------

	    double DeltaR_reco_hlt = hltpfjetraw.DeltaR(recopfjetraw);

	    //hlt jet with smallest DR wrt reco jet
	    if( DeltaR_reco_hlt < min_DeltaR_reco_hlt )
	      {
		min_DeltaR_reco_hlt = DeltaR_reco_hlt;
		idx_hlt_min_DeltaR_fromreco = ihltjet; 
	      }
	    
	    // 	    if(event==113321535 && irecojet==0)
	    // 	      {
	    // 		cout << "HLTPFJetChargedMultiplicity->at(ihltjet)" << HLTPFJetChargedMultiplicity->at(ihltjet) << endl;
	    // 	      }

	  }//end loop over hlt pf jets

	double PtBias_reco_hlt = ( HLTPFJetPt->at(idx_hlt_min_DeltaR_fromreco) - recopfjetraw.Pt() ) / recopfjetraw.Pt();   

	TLorentzVector theHLTjet; 
	theHLTjet.SetPtEtaPhiE(HLTPFJetPt->at(idx_hlt_min_DeltaR_fromreco),
			       HLTPFJetEta->at(idx_hlt_min_DeltaR_fromreco),
			       HLTPFJetPhi->at(idx_hlt_min_DeltaR_fromreco),
			       HLTPFJetEnergy->at(idx_hlt_min_DeltaR_fromreco) );		
	double JetMassBias_reco_hlt = ( theHLTjet.M() - recopfjetraw.M() ) / recopfjetraw.M();	

	if ( fabs( PFJetEta->at(irecojet) ) < 2.4 )
	  {	    
	    if ( PFJetPtRaw->at(irecojet) > 10 )
	      {
		FillUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt10GeV_EtaL2p4"               ,  min_DeltaR_reco_hlt  );
		FillUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4"               ,  PtBias_reco_hlt  );
		FillUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4"          ,  JetMassBias_reco_hlt );
		FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt10GeV_EtaL2p4"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassLooseID->at(irecojet) );
		//FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt10GeV_EtaL2p4"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , passLooseID_recojet );
		FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt10GeV_EtaL2p4"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassTightID->at(irecojet) );
		//FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt10GeV_EtaL2p4"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , passTightID_recojet );
	      }			
	    
	    if ( PFJetPtRaw->at(irecojet) > 30 )
	      {
		FillUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt30GeV_EtaL2p4"             ,  min_DeltaR_reco_hlt  );
		FillUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4"             ,  PtBias_reco_hlt  );
		FillUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4"          ,  JetMassBias_reco_hlt );
		FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt30GeV_EtaL2p4"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassLooseID->at(irecojet) );
		//FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt30GeV_EtaL2p4"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , passLooseID_recojet );
		FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt30GeV_EtaL2p4"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassTightID->at(irecojet) );
		//FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt30GeV_EtaL2p4"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , passTightID_recojet );

		FillUserTH1D( "h_HLTPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"    , HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)      );
		FillUserTH1D( "h_HLTPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"        , HLTPFJetNeutralEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)       );
		FillUserTH1D( "h_HLTPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4"                  , HLTPFJetNConstituents->at(idx_hlt_min_DeltaR_fromreco)      );
		FillUserTH1D( "h_HLTPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"    , HLTPFJetChargedHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)      );
		FillUserTH1D( "h_HLTPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4"            , HLTPFJetChargedMultiplicity->at(idx_hlt_min_DeltaR_fromreco)      );
		FillUserTH1D( "h_HLTPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"        , HLTPFJetChargedEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)      );

		FillUserTH1D( "h_RECOPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"    , PFJetNeutralHadronEnergyFraction->at(irecojet)      );
		FillUserTH1D( "h_RECOPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"        , PFJetNeutralEmEnergyFraction->at(irecojet)      );
		FillUserTH1D( "h_RECOPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4"                  , PFJetNConstituents->at(irecojet)    );
		FillUserTH1D( "h_RECOPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"    , PFJetChargedHadronEnergyFraction->at(irecojet)  );
		FillUserTH1D( "h_RECOPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4"            , PFJetChargedMultiplicity->at(irecojet)     );
		FillUserTH1D( "h_RECOPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4"        , PFJetChargedEmEnergyFraction->at(irecojet)    );

		//if( passLooseID_recojet == 0 && HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) ==1 )
		if( PFJetPassLooseID->at(irecojet) == 0 && HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) ==1 )
		  {
		    if(debug)
		      {
			cout << "LooseRECO0HLT1 - ***** RECO Loose PFJetID = 0  &&  HLT Loose PFJetID = 1 *****" << endl;
			cout << "LooseRECO0HLT1 - run: " << run << " , ls: " << ls << " , event: " << event << endl;
			cout << "LooseRECO0HLT1 - RECOpT: " << recopfjetraw.Pt() << " , HLTpT: " << HLTPFJetPt->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "LooseRECO0HLT1 - RECOeta: " << recopfjetraw.Eta() << " , HLTeta: " << HLTPFJetEta->at(idx_hlt_min_DeltaR_fromreco) << endl;
			// 		    cout << "HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			// 		    cout << "HLTPFJetNeutralEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNeutralEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			// 		    cout << "HLTPFJetNConstituents->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNConstituents->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "LooseRECO0HLT1 - HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "LooseRECO0HLT1 - HLTPFJetChargedHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetChargedHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "LooseRECO1HLT1 - HLTPFJetChargedMultiplicity->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetChargedMultiplicity->at(idx_hlt_min_DeltaR_fromreco) << endl;
			//		    cout << "HLTPFJetChargedEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetChargedEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "LooseRECO0HLT1 - --" << endl;
			cout << "LooseRECO0HLT1 - PFJetNeutralHadronEnergyFraction->at(irecojet): " << PFJetNeutralHadronEnergyFraction->at(irecojet) << endl; 
			cout << "LooseRECO0HLT1 - PFJetChargedHadronEnergyFraction->at(irecojet): " << PFJetChargedHadronEnergyFraction->at(irecojet) << endl;
			cout << "LooseRECO0HLT1 - PFJetChargedMultiplicity->at(irecojet): " << PFJetChargedMultiplicity->at(irecojet) << endl;
		      }

		    FillUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"             ,  min_DeltaR_reco_hlt  );
		    FillUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"             ,  PtBias_reco_hlt  );
		    
		    FillUserTH1D( "h_HLTPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"    , HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)      );
		    FillUserTH1D( "h_HLTPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"        , HLTPFJetNeutralEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)       );
		    FillUserTH1D( "h_HLTPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"                  , HLTPFJetNConstituents->at(idx_hlt_min_DeltaR_fromreco)      );
		    FillUserTH1D( "h_HLTPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"    , HLTPFJetChargedHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)      );
		    FillUserTH1D( "h_HLTPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"            , HLTPFJetChargedMultiplicity->at(idx_hlt_min_DeltaR_fromreco)      );
		    FillUserTH1D( "h_HLTPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"        , HLTPFJetChargedEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)      );
		    
		    FillUserTH1D( "h_RECOPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"    , PFJetNeutralHadronEnergyFraction->at(irecojet)      );
		    FillUserTH1D( "h_RECOPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"        , PFJetNeutralEmEnergyFraction->at(irecojet)      );
		    FillUserTH1D( "h_RECOPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"                  , PFJetNConstituents->at(irecojet)    );
		    FillUserTH1D( "h_RECOPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"    , PFJetChargedHadronEnergyFraction->at(irecojet)  );
		    FillUserTH1D( "h_RECOPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"            , PFJetChargedMultiplicity->at(irecojet)     );
		    FillUserTH1D( "h_RECOPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID0_hltLID1"        , PFJetChargedEmEnergyFraction->at(irecojet)    );
		  }

		//if( passTightID_recojet == 0 && HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) ==1 )
		if( PFJetPassTightID->at(irecojet) == 0 && HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) ==1 )
		  {
		    if(debug)
		      {
			cout << "TightRECO0HLT1 - ***** RECO Tight PFJetID = 0  &&  HLT Tight PFJetID = 1 *****" << endl;
			cout << "TightRECO0HLT1 - run: " << run << " , ls: " << ls << " , event: " << event << endl;
			cout << "TightRECO0HLT1 - RECOpT: " << recopfjetraw.Pt() << " , HLTpT: " << HLTPFJetPt->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "TightRECO0HLT1 - RECOeta: " << recopfjetraw.Eta() << " , HLTeta: " << HLTPFJetEta->at(idx_hlt_min_DeltaR_fromreco) << endl;
			// 		    cout << "HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			// 		    cout << "HLTPFJetNeutralEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNeutralEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			// 		    cout << "HLTPFJetNConstituents->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNConstituents->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "TightRECO0HLT1 - HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "TightRECO0HLT1 - HLTPFJetChargedHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetChargedHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "TightRECO1HLT1 - HLTPFJetChargedMultiplicity->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetChargedMultiplicity->at(idx_hlt_min_DeltaR_fromreco) << endl;
			//		    cout << "HLTPFJetChargedEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetChargedEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "TightRECO0HLT1 - --" << endl;
			cout << "TightRECO0HLT1 - PFJetNeutralHadronEnergyFraction->at(irecojet): " << PFJetNeutralHadronEnergyFraction->at(irecojet) << endl; 
			cout << "TightRECO0HLT1 - PFJetChargedHadronEnergyFraction->at(irecojet): " << PFJetChargedHadronEnergyFraction->at(irecojet) << endl;
			cout << "TightRECO0HLT1 - PFJetChargedMultiplicity->at(irecojet): " << PFJetChargedMultiplicity->at(irecojet) << endl;
		      }
		  }

		//if( passLooseID_recojet == 1 && HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) ==0 )
		if( PFJetPassLooseID->at(irecojet) == 1 && HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) ==0 )
		  {	
		    if(debug)
		      {
			cout << "LooseRECO1HLT0 - ***** RECO Loose PFJetID = 1  &&  HLT Loose PFJetID = 0 *****" << endl;
			cout << "LooseRECO1HLT0 - run: " << run << " , ls: " << ls << " , event: " << event << endl;
			cout << "LooseRECO1HLT0 - RECOpT: " << recopfjetraw.Pt() << " , HLTpT: " << HLTPFJetPt->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "LooseRECO1HLT0 - RECOeta: " << recopfjetraw.Eta() << " , HLTeta: " << HLTPFJetEta->at(idx_hlt_min_DeltaR_fromreco) << endl;
			// 		    cout << "HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			// 		    cout << "HLTPFJetNeutralEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNeutralEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			// 		    cout << "HLTPFJetNConstituents->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNConstituents->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "LooseRECO1HLT0 - HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "LooseRECO1HLT0 - HLTPFJetChargedHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetChargedHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "LooseRECO1HLT0 - HLTPFJetChargedMultiplicity->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetChargedMultiplicity->at(idx_hlt_min_DeltaR_fromreco) << endl;
			//		    cout << "HLTPFJetChargedEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetChargedEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "LooseRECO1HLT0 - --" << endl;
			cout << "LooseRECO1HLT0 - PFJetNeutralHadronEnergyFraction->at(irecojet): " << PFJetNeutralHadronEnergyFraction->at(irecojet) << endl; 
			cout << "LooseRECO1HLT0 - PFJetChargedHadronEnergyFraction->at(irecojet): " << PFJetChargedHadronEnergyFraction->at(irecojet) << endl;
			cout << "LooseRECO1HLT0 - PFJetChargedMultiplicity->at(irecojet): " << PFJetChargedMultiplicity->at(irecojet) << endl;
		      }
		    
		    FillUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"             ,  min_DeltaR_reco_hlt  );
		    FillUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"             ,  PtBias_reco_hlt  );
		    
		    FillUserTH1D( "h_HLTPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"    , HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)      );
		    FillUserTH1D( "h_HLTPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"        , HLTPFJetNeutralEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)       );
		    FillUserTH1D( "h_HLTPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"                  , HLTPFJetNConstituents->at(idx_hlt_min_DeltaR_fromreco)      );
		    FillUserTH1D( "h_HLTPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"    , HLTPFJetChargedHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)      );
		    FillUserTH1D( "h_HLTPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"            , HLTPFJetChargedMultiplicity->at(idx_hlt_min_DeltaR_fromreco)      );
		    FillUserTH1D( "h_HLTPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"        , HLTPFJetChargedEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)      );
		    
		    FillUserTH1D( "h_RECOPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"    , PFJetNeutralHadronEnergyFraction->at(irecojet)      );
		    FillUserTH1D( "h_RECOPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"        , PFJetNeutralEmEnergyFraction->at(irecojet)      );
		    FillUserTH1D( "h_RECOPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"                  , PFJetNConstituents->at(irecojet)    );
		    FillUserTH1D( "h_RECOPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"    , PFJetChargedHadronEnergyFraction->at(irecojet)  );
		    FillUserTH1D( "h_RECOPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"            , PFJetChargedMultiplicity->at(irecojet)     );
		    FillUserTH1D( "h_RECOPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4_recoLID1_hltLID0"        , PFJetChargedEmEnergyFraction->at(irecojet)    );
		  }

		//if( passTightID_recojet == 1 && HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) ==0 )
		if( PFJetPassTightID->at(irecojet) == 1 && HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) ==0 )
		  {
		    if(debug)
		      {
			cout << "TightRECO1HLT0 - ***** RECO Tight PFJetID = 1  &&  HLT Tight PFJetID = 0 *****" << endl;
			cout << "TightRECO1HLT0 - run: " << run << " , ls: " << ls << " , event: " << event << endl;
			cout << "TightRECO1HLT0 - RECOpT: " << recopfjetraw.Pt() << " , HLTpT: " << HLTPFJetPt->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "TightRECO1HLT0 - RECOeta: " << recopfjetraw.Eta() << " , HLTeta: " << HLTPFJetEta->at(idx_hlt_min_DeltaR_fromreco) << endl;
			// 		    cout << "HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			// 		    cout << "HLTPFJetNeutralEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNeutralEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			// 		    cout << "HLTPFJetNConstituents->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNConstituents->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "TightRECO1HLT0 - HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "TightRECO1HLT0 - HLTPFJetChargedHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetChargedHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "TightRECO1HLT0 - HLTPFJetChargedMultiplicity->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetChargedMultiplicity->at(idx_hlt_min_DeltaR_fromreco) << endl;
			//		    cout << "HLTPFJetChargedEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco): " << HLTPFJetChargedEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco) << endl;
			cout << "TightRECO1HLT0 - --" << endl;
			cout << "TightRECO1HLT0 - PFJetNeutralHadronEnergyFraction->at(irecojet): " << PFJetNeutralHadronEnergyFraction->at(irecojet) << endl; 
			cout << "TightRECO1HLT0 - PFJetChargedHadronEnergyFraction->at(irecojet): " << PFJetChargedHadronEnergyFraction->at(irecojet) << endl;
			cout << "TightRECO1HLT0 - PFJetChargedMultiplicity->at(irecojet): " << PFJetChargedMultiplicity->at(irecojet) << endl;
		      }
		  }

	      }			
	    
	    if ( PFJetPtRaw->at(irecojet) > 50 )
	      {
		FillUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt50GeV_EtaL2p4"             ,  min_DeltaR_reco_hlt  );
		FillUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4"             ,  PtBias_reco_hlt  );
		FillUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4"          ,  JetMassBias_reco_hlt );
		FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt50GeV_EtaL2p4"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassLooseID->at(irecojet) );
		//FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt50GeV_EtaL2p4"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , passLooseID_recojet );
		FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt50GeV_EtaL2p4"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassTightID->at(irecojet) );
		//FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt50GeV_EtaL2p4"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , passTightID_recojet );
	      }			
	    
	    if ( PFJetPtRaw->at(irecojet) > 100 )
	      {
		FillUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt100GeV_EtaL2p4"             ,  min_DeltaR_reco_hlt  );
		FillUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4"             ,  PtBias_reco_hlt  );
		FillUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4"          ,  JetMassBias_reco_hlt );
		FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt100GeV_EtaL2p4"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassLooseID->at(irecojet) );
		//FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt100GeV_EtaL2p4"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , passLooseID_recojet );
		FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt100GeV_EtaL2p4"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassTightID->at(irecojet) );
		//FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt100GeV_EtaL2p4"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , passTightID_recojet );
	      }
	  }//end if |recojet eta| < 2.4

	if ( fabs( PFJetEta->at(irecojet) ) >= 2.4 && fabs( PFJetEta->at(irecojet) ) < 3 )
	  {	    
	    if ( PFJetPtRaw->at(irecojet) > 10 )
	      {
		FillUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3"             ,  min_DeltaR_reco_hlt  );
		FillUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3"             ,  PtBias_reco_hlt  );
		FillUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3"          ,  JetMassBias_reco_hlt );
		FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassLooseID->at(irecojet) );
		//FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , passLooseID_recojet );
		FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassTightID->at(irecojet) );
		//FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , passTightID_recojet );
	      }			
	    
	    if ( PFJetPtRaw->at(irecojet) > 30 )
	      {
		FillUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3"             ,  min_DeltaR_reco_hlt  );
		FillUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3"             ,  PtBias_reco_hlt  );
		FillUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3"          ,  JetMassBias_reco_hlt );
		FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassLooseID->at(irecojet) );
		//FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , passLooseID_recojet );
		FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassTightID->at(irecojet) );
		//FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , passTightID_recojet );

		FillUserTH1D( "h_HLTPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"    , HLTPFJetNeutralHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)      );
		FillUserTH1D( "h_HLTPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"        , HLTPFJetNeutralEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)       );
		FillUserTH1D( "h_HLTPFJetNConstituents_smallestDR_Pt30GeV_Eta2p4_3"                  , HLTPFJetNConstituents->at(idx_hlt_min_DeltaR_fromreco)      );
		FillUserTH1D( "h_HLTPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"    , HLTPFJetChargedHadronEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)      );
		FillUserTH1D( "h_HLTPFJetChargedMultiplicity_smallestDR_Pt30GeV_Eta2p4_3"            , HLTPFJetChargedMultiplicity->at(idx_hlt_min_DeltaR_fromreco)      );
		FillUserTH1D( "h_HLTPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"        , HLTPFJetChargedEmEnergyFraction->at(idx_hlt_min_DeltaR_fromreco)      );

		FillUserTH1D( "h_RECOPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"    , PFJetNeutralHadronEnergyFraction->at(irecojet)      );
		FillUserTH1D( "h_RECOPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"        , PFJetNeutralEmEnergyFraction->at(irecojet)      );
		FillUserTH1D( "h_RECOPFJetNConstituents_smallestDR_Pt30GeV_Eta2p4_3"                  , PFJetNConstituents->at(irecojet)    );
		FillUserTH1D( "h_RECOPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"    , PFJetChargedHadronEnergyFraction->at(irecojet)  );
		FillUserTH1D( "h_RECOPFJetChargedMultiplicity_smallestDR_Pt30GeV_Eta2p4_3"            , PFJetChargedMultiplicity->at(irecojet)     );
		FillUserTH1D( "h_RECOPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_Eta2p4_3"        , PFJetChargedEmEnergyFraction->at(irecojet)    );
	      }			
	    
	    if ( PFJetPtRaw->at(irecojet) > 50 )
	      {
		FillUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3"             ,  min_DeltaR_reco_hlt  );
		FillUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3"             ,  PtBias_reco_hlt  );
		FillUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3"          ,  JetMassBias_reco_hlt );
		FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassLooseID->at(irecojet) );
		//FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , passLooseID_recojet );
		FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassTightID->at(irecojet) );
		//FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , passTightID_recojet );
	      }			
	    
	    if ( PFJetPtRaw->at(irecojet) > 100 )
	      {
		FillUserTH1D( "h_DeltaR_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3"             ,  min_DeltaR_reco_hlt  );
		FillUserTH1D( "h_PtBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3"             ,  PtBias_reco_hlt  );
		FillUserTH1D( "h_JetMassBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3"          ,  JetMassBias_reco_hlt );
		FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassLooseID->at(irecojet) );
		//FillUserTH2D( "h2_LooseJetIDAgreement_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3"  ,  HLTPFJetPassLooseID->at(idx_hlt_min_DeltaR_fromreco) , passLooseID_recojet );
		FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , PFJetPassTightID->at(irecojet) );
		//FillUserTH2D( "h2_TightJetIDAgreement_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3"  ,  HLTPFJetPassTightID->at(idx_hlt_min_DeltaR_fromreco) , passTightID_recojet );
	      }
	  }//end if 2.4 < |recojet eta| < 3

      }//end loop over reco pf jets


    ////////////////////// Object Collections ///////////////////////

    
    //## HLT PF Jets
    vector<int> v_idx_hltpfjet_PtEtaCut;
    vector<int> v_idx_hltpfjet_PtEtaCut_ID;
    
    // Precuts
    for(int ijet=0; ijet<HLTPFJetCorrPt->size(); ijet++)
      {
	//pT/eta pre-cuts on jets
	if ( HLTPFJetCorrPt->at(ijet) < pfjetPtCut ) continue;
	if ( fabs( HLTPFJetEta->at(ijet) ) > pfjetEtaCut ) continue;
	v_idx_hltpfjet_PtEtaCut.push_back(ijet);
      }

    // Full selection
    double HT_HLTPFJets = 0.; 
    double MHT_HLTPFJets = 0.; 
    double MHTPhi_HLTPFJets = 0.; 
    TVector2 v_MHT_HLTPFJets; 
    v_MHT_HLTPFJets.SetMagPhi(0.,0.);
    for(int ijet=0; ijet<v_idx_hltpfjet_PtEtaCut.size(); ijet++) 
      {	
	bool passjetID = true;
	
	if( pfjetIDloose && !pfjetIDtight )
	  {
	    passjetID = HLTPFJetPassLooseID->at(v_idx_hltpfjet_PtEtaCut[ijet]);
	  }

	if( pfjetIDtight )
	  {
	    passjetID = HLTPFJetPassTightID->at(v_idx_hltpfjet_PtEtaCut[ijet]);
	  }

	if( passjetID == true )
	  {
	    v_idx_hltpfjet_PtEtaCut_ID.push_back(v_idx_hltpfjet_PtEtaCut[ijet]);

	    //HT and MHT calculation
	    if( HLTPFJetCorrPt->at(v_idx_hltpfjet_PtEtaCut[ijet]) > pfjetPtCutForHT  
		&& fabs( HLTPFJetEta->at(v_idx_hltpfjet_PtEtaCut[ijet]) ) < pfjetEtaCutForHT )
	      {
		HT_HLTPFJets += HLTPFJetCorrPt->at(v_idx_hltpfjet_PtEtaCut[ijet]);		

		TVector2 currentJet;                                                      
		currentJet.SetMagPhi( HLTPFJetCorrPt->at(v_idx_hltpfjet_PtEtaCut[ijet]) , HLTPFJetPhi->at(v_idx_hltpfjet_PtEtaCut[ijet]) );
		v_MHT_HLTPFJets += currentJet;//add to SumET vector
	      }
	  }
      } // End loop over hlt pf jets

    //(continue) MHT
    v_MHT_HLTPFJets.Rotate( TMath::Pi() );   
    MHT_HLTPFJets = v_MHT_HLTPFJets.Mod();
    MHTPhi_HLTPFJets = v_MHT_HLTPFJets.Phi_mpi_pi( v_MHT_HLTPFJets.Phi() );

    //## RECO PF Jets
    vector<int> v_idx_recopfjet_PtEtaCut;
    vector<int> v_idx_recopfjet_PtEtaCut_ID;
    
    // Precuts
    for(int ijet=0; ijet<PFJetPt->size(); ijet++)
      {
	//pT/eta pre-cuts on jets
	if ( PFJetPt->at(ijet) < pfjetPtCut ) continue;
	if ( fabs( PFJetEta->at(ijet) ) > pfjetEtaCut ) continue;
	v_idx_recopfjet_PtEtaCut.push_back(ijet);

	// 	//Test JEC
	//
	//      -- Validated!
	//
	// 	JetCorrectorL2->setJetEta( PFJetEta->at(ijet) );
	// 	JetCorrectorL2->setJetPt( PFJetL1FastJetJEC->at(ijet) * PFJetPtRaw->at(ijet) );
	// 	double thisL2correction = JetCorrectorL2->getCorrection();
	
	// 	JetCorrectorL3->setJetEta( PFJetEta->at(ijet) );
	// 	JetCorrectorL3->setJetPt( PFJetL1FastJetJEC->at(ijet) * thisL2correction * PFJetPtRaw->at(ijet) );
	// 	double thisL3correction = JetCorrectorL3->getCorrection();
	
	// 	JetCorrectorL23Res->setJetEta( PFJetEta->at(ijet) );
	// 	JetCorrectorL23Res->setJetPt( PFJetL1FastJetJEC->at(ijet) * thisL2correction * thisL3correction * PFJetPtRaw->at(ijet) );
	// 	double thisL23Rescorrection = JetCorrectorL23Res->getCorrection();
	
	// 	cout << "---------------" << endl;
	// 	cout << "From ntuple" << endl;
	// 	cout << "PT: " << PFJetPt->at(ijet) 
	// 	     << " PTraw: " << PFJetPtRaw->at(ijet)
	// 	     << " L1: " << PFJetL1FastJetJEC->at(ijet)
	// 	     << " L2: " << PFJetL2RelJEC->at(ijet)	  
	// 	     << " L3: " << PFJetL3AbsJEC->at(ijet)	  
	// 	     << " L2L3Res:" << PFJetL2L3ResJEC->at(ijet)
	// 	  //<< " PTraw * corr: " << PFJetPtRaw->at(ijet) * PFJetL1FastJetJEC->at(ijet) * PFJetL2RelJEC->at(ijet) * PFJetL3AbsJEC->at(ijet) * PFJetL2L3ResJEC->at(ijet)
	// 	     << endl;
	// 	cout << "From file" << endl;
	// 	cout << "PT: " << PFJetPt->at(ijet) 
	// 	     << " PTraw: " << PFJetPtRaw->at(ijet)
	// 	     << " L1: " << PFJetL1FastJetJEC->at(ijet)
	// 	     << " L2: " << thisL2correction	  
	// 	     << " L3: " << thisL3correction	  
	// 	     << " L2L3Res:" << thisL23Rescorrection	  
	// 	     << endl;
	// 	cout << "---------------" << endl;	
	//
      }
    
    // Full selection
    double HT_RECOPFJets = 0.; 
    double MHT_RECOPFJets = 0.; 
    double MHTPhi_RECOPFJets = 0.; 
    TVector2 v_MHT_RECOPFJets; 
    v_MHT_RECOPFJets.SetMagPhi(0.,0.);
    for(int ijet=0; ijet<v_idx_recopfjet_PtEtaCut.size(); ijet++) 
      {	
	bool passjetID = true;
	
	if( pfjetIDloose && !pfjetIDtight )
	  {
	    passjetID = PFJetPassLooseID->at(v_idx_recopfjet_PtEtaCut[ijet]);
	  }
	
	if( pfjetIDtight )
	  {
	    passjetID = PFJetPassTightID->at(v_idx_recopfjet_PtEtaCut[ijet]);
	  }
	
	if( passjetID == true )
	  {
	    v_idx_recopfjet_PtEtaCut_ID.push_back(v_idx_recopfjet_PtEtaCut[ijet]);
	    
	    //HT and MHT calculation
	    if( PFJetPt->at(v_idx_recopfjet_PtEtaCut[ijet]) > pfjetPtCutForHT  
		&& fabs( PFJetEta->at(v_idx_recopfjet_PtEtaCut[ijet]) ) < pfjetEtaCutForHT )
	      {
		HT_RECOPFJets += PFJetPt->at(v_idx_recopfjet_PtEtaCut[ijet]);		

		TVector2 currentJet;                                                      
		currentJet.SetMagPhi( PFJetPt->at(v_idx_recopfjet_PtEtaCut[ijet]) , PFJetPhi->at(v_idx_recopfjet_PtEtaCut[ijet]) );
		v_MHT_RECOPFJets += currentJet;//add to SumET vector
	      }
	  }
      } // End loop over reco pf jets

    //(continue) MHT
    v_MHT_RECOPFJets.Rotate( TMath::Pi() );   
    MHT_RECOPFJets = v_MHT_RECOPFJets.Mod();
    MHTPhi_RECOPFJets = v_MHT_RECOPFJets.Phi_mpi_pi( v_MHT_RECOPFJets.Phi() );
    
    //## Pixel HLT Vertices
    vector<int> v_idx_hltvertex_good;
      
    // Full selection
    for(int ivertex = 0; ivertex<HLTPixelVertexZCoord->size(); ivertex++)
      {
	
	double vertex_rho = sqrt( 
				 HLTPixelVertexXCoord->at(ivertex)*HLTPixelVertexXCoord->at(ivertex) +
				 HLTPixelVertexYCoord->at(ivertex)*HLTPixelVertexYCoord->at(ivertex) 
				 );
	
	if ( HLTPixelVertexIsValid->at(ivertex) 
	     && fabs( HLTPixelVertexZCoord->at(ivertex) ) <= vertexMaxAbsZ
	     && vertex_rho <= vertexMaxd0 
	     )
	  {
	    v_idx_hltvertex_good.push_back(ivertex);
	  }
      }

    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();

    // Set the value of the variableNames listed in the cutFile to their current value

    // Event filters
    fillVariableWithValue( "PassJSON", passJSON(run, ls, isData) );    

    int isHLTPrimaryVertex = 0;
    if ( v_idx_hltvertex_good.size() > 0 )
      isHLTPrimaryVertex = 1;
    fillVariableWithValue( "PassHLTPrimaryVertex", isHLTPrimaryVertex ) ;
    
    // Event info
    fillVariableWithValue( "isData", isData ) ;
    fillVariableWithValue( "bunch", bunch ) ;
    fillVariableWithValue( "event", event ) ;
    fillVariableWithValue( "ls", ls ) ;
    fillVariableWithValue( "orbit", orbit ) ;
    fillVariableWithValue( "run", run ) ;

    // Trigger (L1 and HLT)
    if(isData==true)
      {
	fillVariableWithValue( "PassBPTX0", isBPTX0 ) ;
      }
    else
      {
	fillVariableWithValue( "PassBPTX0", true ) ;
      }

    //Event filters at RECO level
    fillVariableWithValue( "PassBeamScraping", !isBeamScraping ) ;
    fillVariableWithValue( "PassRECOPrimaryVertex", isPrimaryVertex ) ;
    fillVariableWithValue( "PassHBHENoiseFilter", passHBHENoiseFilter ) ;
    fillVariableWithValue( "PassBeamHaloFilterTight", passBeamHaloFilterTight ) ;
    fillVariableWithValue( "PassTrackingFailure", !isTrackingFailure ) ;

    //Njets
    fillVariableWithValue( "HLT_N_PFJet", v_idx_hltpfjet_PtEtaCut_ID.size() ) ;
    fillVariableWithValue( "RECO_N_PFJet", v_idx_recopfjet_PtEtaCut_ID.size() ) ;

    //SumET and MET
    fillVariableWithValue( "HT_HLTPFJets", HT_HLTPFJets ) ;
    fillVariableWithValue( "MHT_HLTPFJets", MHT_HLTPFJets ) ;
    fillVariableWithValue( "MHTPhi_HLTPFJets", MHTPhi_HLTPFJets ) ;
    if( HT_HLTPFJets > 0 )
      fillVariableWithValue( "MHTSig_HLTPFJets", MHT_HLTPFJets/sqrt(HT_HLTPFJets) ) ;

    fillVariableWithValue( "HT_RECOPFJets", HT_RECOPFJets ) ;
    fillVariableWithValue( "MHT_RECOPFJets", MHT_RECOPFJets ) ;
    fillVariableWithValue( "MHTPhi_RECOPFJets", MHTPhi_RECOPFJets ) ;
    if( HT_RECOPFJets > 0 )
      fillVariableWithValue( "MHTSig_RECOPFJets", MHT_RECOPFJets/sqrt(HT_RECOPFJets) ) ;

    // 1st jet 
    if( HLTPFJetCorrPt->size() >=1 )
      {
	bool passjetID = true;
	
	if( pfjetIDloose && !pfjetIDtight )
	  {
	    passjetID = HLTPFJetPassLooseID->at(0);
	  }

	if( pfjetIDtight )
	  {
	    passjetID = HLTPFJetPassTightID->at(0);
	  }

	fillVariableWithValue( "HLT_PFJet1_Pt", HLTPFJetCorrPt->at(0) );
	fillVariableWithValue( "HLT_PFJet1_Energy", HLTPFJetCorrEnergy->at(0) );
	fillVariableWithValue( "HLT_PFJet1_Eta", HLTPFJetEta->at(0) );
	fillVariableWithValue( "HLT_PFJet1_Phi", HLTPFJetPhi->at(0) );
	fillVariableWithValue( "HLT_PFJet1_PassJetID", passjetID );
      }

    if( PFJetPt->size() >= 1 )        
      {
	bool passjetID = true;
	
	if( pfjetIDloose && !pfjetIDtight )
	  {
	    passjetID = PFJetPassLooseID->at(0);
	  }
	
	if( pfjetIDtight )
	  {
	    passjetID = PFJetPassTightID->at(0);
	  }

	fillVariableWithValue( "RECO_PFJet1_Pt", PFJetPt->at(0) );
	fillVariableWithValue( "RECO_PFJet1_Energy", PFJetEnergy->at(0) );
	fillVariableWithValue( "RECO_PFJet1_Eta", PFJetEta->at(0) );
	fillVariableWithValue( "RECO_PFJet1_Phi", PFJetPhi->at(0) );
	fillVariableWithValue( "RECO_PFJet1_PassJetID", passjetID );
      }

    // 2nd jet 
    if( HLTPFJetCorrPt->size() >= 2 )        
      {
	bool passjetID = true;
	
	if( pfjetIDloose && !pfjetIDtight )
	  {
	    passjetID = HLTPFJetPassLooseID->at(1);
	  }

	if( pfjetIDtight )
	  {
	    passjetID = HLTPFJetPassTightID->at(1);
	  }

	fillVariableWithValue( "HLT_PFJet2_Pt", HLTPFJetCorrPt->at(1) );
	fillVariableWithValue( "HLT_PFJet2_Energy", HLTPFJetCorrEnergy->at(1) );
	fillVariableWithValue( "HLT_PFJet2_Eta", HLTPFJetEta->at(1) );
	fillVariableWithValue( "HLT_PFJet2_Phi", HLTPFJetPhi->at(1) );
	fillVariableWithValue( "HLT_PFJet2_PassJetID", passjetID );

	if(debug)
	  {
	    if( HLTPFJetCorrPt->at(1) > HLTPFJetCorrPt->at(0) )
	      {
		cout << "diff , rel diff: " 
		     << HLTPFJetCorrPt->at(1)  - HLTPFJetCorrPt->at(0) 
		     << " , " 
		     << ( HLTPFJetCorrPt->at(1)  - HLTPFJetCorrPt->at(0) ) / HLTPFJetCorrPt->at(0)
		     << endl;
	      }
	  }
      }

    if( PFJetPt->size() >= 2 )        
      {
	bool passjetID = true;
	
	if( pfjetIDloose && !pfjetIDtight )
	  {
	    passjetID = PFJetPassLooseID->at(1);
	  }
	
	if( pfjetIDtight )
	  {
	    passjetID = PFJetPassTightID->at(1);
	  }

	fillVariableWithValue( "RECO_PFJet2_Pt", PFJetPt->at(1) );
	fillVariableWithValue( "RECO_PFJet2_Energy", PFJetEnergy->at(1) );
	fillVariableWithValue( "RECO_PFJet2_Eta", PFJetEta->at(1) );
	fillVariableWithValue( "RECO_PFJet2_Phi", PFJetPhi->at(1) );
	fillVariableWithValue( "RECO_PFJet2_PassJetID", passjetID );

	if(debug)
	  {
	    if( PFJetPt->at(1) > PFJetPt->at(0) )
	      cout << "cannot happen" << endl;
	  }
      }
       
    // define booleans
    bool TwoHLTPFJets=false;
    bool TwoRECOPFJets=false;
    if( HLTPFJetCorrPt->size() >= 2 )  TwoHLTPFJets = true;
    if( PFJetPt->size() >= 2 ) TwoRECOPFJets = true;

    if (TwoHLTPFJets)
      {
	TLorentzVector jet1, jet2, jj;
	jet1.SetPtEtaPhiE(HLTPFJetCorrPt->at(0),
			  HLTPFJetEta->at(0),
			  HLTPFJetPhi->at(0),
			  HLTPFJetCorrEnergy->at(0) );
	jet2.SetPtEtaPhiE(HLTPFJetCorrPt->at(1),
			  HLTPFJetEta->at(1),
			  HLTPFJetPhi->at(1),
			  HLTPFJetCorrEnergy->at(1) );
	jj = jet1+jet2;
	
	fillVariableWithValue("HLT_DR_PFJet1PFJet2", jet1.DeltaR(jet2));
	fillVariableWithValue("HLT_DPhi_PFJet1PFJet2", jet1.DeltaPhi(jet2));	
	fillVariableWithValue("HLT_DEta_PFJet1PFJet2", fabs(jet1.Eta()-jet2.Eta()) );
	fillVariableWithValue("HLT_Pt_PFJet1PFJet2", jj.Pt());
	fillVariableWithValue("HLT_M_PFJet1PFJet2", jj.M());
      }

    if (TwoRECOPFJets)
      {
	TLorentzVector jet1, jet2, jj;
	jet1.SetPtEtaPhiE(PFJetPt->at(0),
			  PFJetEta->at(0),
			  PFJetPhi->at(0),
			  PFJetEnergy->at(0) );
	jet2.SetPtEtaPhiE(PFJetPt->at(1),
			  PFJetEta->at(1),
			  PFJetPhi->at(1),
			  PFJetEnergy->at(1) );
	jj = jet1+jet2;

	fillVariableWithValue("RECO_DR_PFJet1PFJet2", jet1.DeltaR(jet2));
	fillVariableWithValue("RECO_DPhi_PFJet1PFJet2", jet1.DeltaPhi(jet2));	
	fillVariableWithValue("RECO_DEta_PFJet1PFJet2", fabs(jet1.Eta()-jet2.Eta()) );
	fillVariableWithValue("RECO_Pt_PFJet1PFJet2", jj.Pt());
	fillVariableWithValue("RECO_M_PFJet1PFJet2", jj.M());
      }

    // Evaluate cuts (but do not apply them)
    evaluateCuts();

    
    // Fill histograms and do analysis based on cut evaluation
    
    int pass_HLT_filters = 0;
    int pass_RECO_filters = 0;
    int pass_HLT_selection = 0;
    int pass_RECO_selection = 0;
    
    //HLT selection
    if( passedCut("PassJSON") 
	&& passedCut("PassHLTPrimaryVertex") 
	)
      {
	
	pass_HLT_filters = 1;
	
	if(  passedCut("HLT_PFJet1_Pt") 
	     && passedCut("HLT_PFJet2_Pt") 
	     && passedCut("HLT_PFJet1_Eta")
	     && passedCut("HLT_PFJet2_Eta")
	     && passedCut("HLT_PFJet1_PassJetID")
	     && passedCut("HLT_PFJet2_PassJetID")
	     && passedCut("HLT_DPhi_PFJet1PFJet2")
	     && passedCut("HLT_DEta_PFJet1PFJet2")
	     )
	  {
	    pass_HLT_selection = 1;
	  }
      }

    //RECO selection
    if( passedCut("PassJSON") 
	&& passedCut("PassBPTX0") 
	&& passedCut("PassBeamScraping") 
	&& passedCut("PassHBHENoiseFilter") 
	&& passedCut("PassBeamHaloFilterTight") 
	&& passedCut("PassTrackingFailure") 
	&& passedCut("PassRECOPrimaryVertex") 
	)
      {  
	
	pass_RECO_filters = 1;

	if( passedCut("RECO_PFJet1_Pt") 
	    && passedCut("RECO_PFJet2_Pt") 
	    && passedCut("RECO_PFJet1_Eta")
	    && passedCut("RECO_PFJet2_Eta")
	    && passedCut("RECO_PFJet1_PassJetID")
	    && passedCut("RECO_PFJet2_PassJetID")
	    && passedCut("RECO_DPhi_PFJet1PFJet2")
	    && passedCut("RECO_DEta_PFJet1PFJet2")
	    )
	  {
	    pass_RECO_selection = 1;
	  }
      }

    FillUserTH2D( "h2_FinalSelectionAgreement_recoY_vs_hltX" , pass_HLT_selection , pass_RECO_selection );    

    if( pass_HLT_selection == 1 )
      {
	FillUserTH1D("h_HLT_PFJet1_Pt", getVariableValue("HLT_PFJet1_Pt") );
	FillUserTH1D("h_HLT_PFJet1_Eta", getVariableValue("HLT_PFJet1_Eta") );
	FillUserTH1D("h_HLT_PFJet2_Pt", getVariableValue("HLT_PFJet2_Pt") );
	FillUserTH1D("h_HLT_PFJet2_Eta", getVariableValue("HLT_PFJet2_Eta") );
	FillUserTH1D("h_HLT_DR_PFJet1PFJet2", getVariableValue("HLT_DR_PFJet1PFJet2") );
	FillUserTH1D("h_HLT_DPhi_PFJet1PFJet2", getVariableValue("HLT_DPhi_PFJet1PFJet2") );
	FillUserTH1D("h_HLT_DEta_PFJet1PFJet2", getVariableValue("HLT_DEta_PFJet1PFJet2") );
	FillUserTH1D("h_HLT_Pt_PFJet1PFJet2", getVariableValue("HLT_Pt_PFJet1PFJet2") );
	FillUserTH1D("h_HLT_M_PFJet1PFJet2", getVariableValue("HLT_M_PFJet1PFJet2") );

	CreateAndFillUserTH1D("h_HLT_N_PFJet", 51, -0.5, 50.5, getVariableValue("HLT_N_PFJet") );
      }

    if( pass_RECO_selection == 1 )
      {
	FillUserTH1D("h_RECO_PFJet1_Pt", getVariableValue("RECO_PFJet1_Pt") );
	FillUserTH1D("h_RECO_PFJet1_Eta", getVariableValue("RECO_PFJet1_Eta") );
	FillUserTH1D("h_RECO_PFJet2_Pt", getVariableValue("RECO_PFJet2_Pt") );
	FillUserTH1D("h_RECO_PFJet2_Eta", getVariableValue("RECO_PFJet2_Eta") );
	FillUserTH1D("h_RECO_DR_PFJet1PFJet2", getVariableValue("RECO_DR_PFJet1PFJet2") );
	FillUserTH1D("h_RECO_DPhi_PFJet1PFJet2", getVariableValue("RECO_DPhi_PFJet1PFJet2") );
	FillUserTH1D("h_RECO_DEta_PFJet1PFJet2", getVariableValue("RECO_DEta_PFJet1PFJet2") );
	FillUserTH1D("h_RECO_Pt_PFJet1PFJet2", getVariableValue("RECO_Pt_PFJet1PFJet2") );
	FillUserTH1D("h_RECO_M_PFJet1PFJet2", getVariableValue("RECO_M_PFJet1PFJet2") );

	CreateAndFillUserTH1D("h_RECO_N_PFJet", 51, -0.5, 50.5, getVariableValue("RECO_N_PFJet") );
      }

    if( pass_HLT_selection == 1 && pass_RECO_selection == 1 )
      {       
	FillUserTH1D("h_HLT_PFJet1_Pt_MATCH", getVariableValue("HLT_PFJet1_Pt") );
	FillUserTH1D("h_HLT_PFJet1_Eta_MATCH", getVariableValue("HLT_PFJet1_Eta") );
	FillUserTH1D("h_HLT_PFJet2_Pt_MATCH", getVariableValue("HLT_PFJet2_Pt") );
	FillUserTH1D("h_HLT_PFJet2_Eta_MATCH", getVariableValue("HLT_PFJet2_Eta") );
	FillUserTH1D("h_HLT_DR_PFJet1PFJet2_MATCH", getVariableValue("HLT_DR_PFJet1PFJet2") );
	FillUserTH1D("h_HLT_DPhi_PFJet1PFJet2_MATCH", getVariableValue("HLT_DPhi_PFJet1PFJet2") );
	FillUserTH1D("h_HLT_DEta_PFJet1PFJet2_MATCH", getVariableValue("HLT_DEta_PFJet1PFJet2") );
	FillUserTH1D("h_HLT_Pt_PFJet1PFJet2_MATCH", getVariableValue("HLT_Pt_PFJet1PFJet2") );
	FillUserTH1D("h_HLT_M_PFJet1PFJet2_MATCH", getVariableValue("HLT_M_PFJet1PFJet2") );

	FillUserTH1D("h_RECO_PFJet1_Pt_MATCH", getVariableValue("RECO_PFJet1_Pt") );
	FillUserTH1D("h_RECO_PFJet1_Eta_MATCH", getVariableValue("RECO_PFJet1_Eta") );
	FillUserTH1D("h_RECO_PFJet2_Pt_MATCH", getVariableValue("RECO_PFJet2_Pt") );
	FillUserTH1D("h_RECO_PFJet2_Eta_MATCH", getVariableValue("RECO_PFJet2_Eta") );
	FillUserTH1D("h_RECO_DR_PFJet1PFJet2_MATCH", getVariableValue("RECO_DR_PFJet1PFJet2") );
	FillUserTH1D("h_RECO_DPhi_PFJet1PFJet2_MATCH", getVariableValue("RECO_DPhi_PFJet1PFJet2") );
	FillUserTH1D("h_RECO_DEta_PFJet1PFJet2_MATCH", getVariableValue("RECO_DEta_PFJet1PFJet2") );
	FillUserTH1D("h_RECO_Pt_PFJet1PFJet2_MATCH", getVariableValue("RECO_Pt_PFJet1PFJet2") );
	FillUserTH1D("h_RECO_M_PFJet1PFJet2_MATCH", getVariableValue("RECO_M_PFJet1PFJet2") );

	double mass_bias = ( getVariableValue("HLT_M_PFJet1PFJet2") - getVariableValue("RECO_M_PFJet1PFJet2") ) / getVariableValue("RECO_M_PFJet1PFJet2"); 
	double deta_bias = ( getVariableValue("HLT_DEta_PFJet1PFJet2") - getVariableValue("RECO_DEta_PFJet1PFJet2") ) / getVariableValue("RECO_DEta_PFJet1PFJet2"); 

	FillUserTH1D( "h_MassBias_MATCH"  ,  mass_bias  );
	FillUserTH1D( "h_DEtaBias_MATCH"  ,  deta_bias  );

	if( fabs(mass_bias) > 0.2 )
	  {
	    FillUserTH1D( "h_HLT_M_PFJet1PFJet2_MATCH_MassBiasMore0p2", getVariableValue("HLT_M_PFJet1PFJet2") );
	    FillUserTH1D( "h_RECO_M_PFJet1PFJet2_MATCH_MassBiasMore0p2", getVariableValue("RECO_M_PFJet1PFJet2") );
	  }

	p_MassBias_MATCH_vs_Mreco->Fill( getVariableValue("RECO_M_PFJet1PFJet2") , mass_bias );

      }

    if( pass_HLT_selection == 1 && pass_RECO_selection == 0 )
      {	
	int mismatch_filters = 0;
	int mismatch_Jet1Pt = 0;
	int mismatch_Jet1Eta = 0;
	int mismatch_Jet1ID = 0;
	int mismatch_Jet2Pt = 0;
	int mismatch_Jet2Eta = 0;
	int mismatch_Jet2ID = 0;
	int mismatch_DPhi = 0;
	int mismatch_DEta = 0;
	
	if( pass_RECO_filters != pass_HLT_filters )
	  mismatch_filters = 1;
	
	if( passedCut("HLT_PFJet1_Pt") != passedCut("RECO_PFJet1_Pt") )
	  mismatch_Jet1Pt = 1;
	
	if( passedCut("HLT_PFJet1_Eta") != passedCut("RECO_PFJet1_Eta") )
	  mismatch_Jet1Eta = 1;
	
	if( passedCut("HLT_PFJet1_PassJetID") != passedCut("RECO_PFJet1_PassJetID") )
	  mismatch_Jet1ID = 1;

	if( passedCut("HLT_PFJet2_Pt") != passedCut("RECO_PFJet2_Pt") )
	  mismatch_Jet2Pt = 1;
	
	if( passedCut("HLT_PFJet2_Eta") != passedCut("RECO_PFJet2_Eta") )
	  mismatch_Jet2Eta = 1;

	if( passedCut("HLT_PFJet2_PassJetID") != passedCut("RECO_PFJet2_PassJetID") )
	  mismatch_Jet2ID = 1;

	if( passedCut("HLT_DPhi_PFJet1PFJet2") != passedCut("RECO_DPhi_PFJet1PFJet2") )
	  mismatch_DPhi = 1;
	
	if( passedCut("HLT_DEta_PFJet1PFJet2") != passedCut("RECO_DEta_PFJet1PFJet2") )
	  mismatch_DEta = 1;

	if(debugEv)
	  {
	    cout << "SelectionHLT1RECO0 ----------" << endl;
	    cout << "SelectionHLT1RECO0 - run: " << run << " , ls: " << ls << " , event: " << event << endl;
	  }

	if(mismatch_filters==1)
	  {
	    h_FirstMismatch_HLT1_RECO0->Fill(0);

	    if(debugEv)
	      {
		cout << "SelectionHLT1RECO0 - PassHLTPrimaryVertex: " << getVariableValue("PassHLTPrimaryVertex") << endl;
		cout << "SelectionHLT1RECO0 - PassRECOPrimaryVertex: " << getVariableValue("PassRECOPrimaryVertex") << endl;
		cout << "SelectionHLT1RECO0 - PassBPTX0: " << getVariableValue("PassBPTX0") << endl;
		cout << "SelectionHLT1RECO0 - PassBeamScraping: " << getVariableValue("PassBeamScraping") << endl;
		cout << "SelectionHLT1RECO0 - PassHBHENoiseFilter: " << getVariableValue("PassHBHENoiseFilter") << endl;
		cout << "SelectionHLT1RECO0 - PassBeamHaloFilterTight: " << getVariableValue("PassBeamHaloFilterTight") << endl;
		cout << "SelectionHLT1RECO0 - PassTrackingFailure: " << getVariableValue("PassTrackingFailure") << endl;
	      }

	    CreateAndFillUserTH1D("h_HLT_DPhi_PFJet1PFJet2_HLT1_RECO0_mismatch_filters", getHistoNBins("HLT_DPhi_PFJet1PFJet2"), getHistoMin("HLT_DPhi_PFJet1PFJet2"), getHistoMax("HLT_DPhi_PFJet1PFJet2"), getVariableValue("HLT_DPhi_PFJet1PFJet2") );
	    CreateAndFillUserTH1D("h_RECO_DPhi_PFJet1PFJet2_HLT1_RECO0_mismatch_filters", getHistoNBins("RECO_DPhi_PFJet1PFJet2"), getHistoMin("RECO_DPhi_PFJet1PFJet2"), getHistoMax("RECO_DPhi_PFJet1PFJet2"), getVariableValue("RECO_DPhi_PFJet1PFJet2") );	    

	    CreateAndFillUserTH1D("h_HLT_N_PFJet_HLT1_RECO0_mismatch_filters", 51, -0.5, 50.5, v_idx_hltpfjet_PtEtaCut_ID.size() );
	    CreateAndFillUserTH1D("h_RECO_N_PFJet_HLT1_RECO0_mismatch_filters", 51, -0.5, 50.5, v_idx_recopfjet_PtEtaCut_ID.size() );
	  }
	if(mismatch_Jet1Pt==1 && mismatch_filters==0)
	  h_FirstMismatch_HLT1_RECO0->Fill(1);
	if(mismatch_Jet1Eta==1 && mismatch_filters==0 && mismatch_Jet1Pt==0)
	  h_FirstMismatch_HLT1_RECO0->Fill(2);
	if(mismatch_Jet1ID==1 && mismatch_Jet1Eta==0 && mismatch_filters==0 && mismatch_Jet1Pt==0)
	  h_FirstMismatch_HLT1_RECO0->Fill(3);
	if(mismatch_Jet2Pt==1 && mismatch_Jet1ID==0 && mismatch_Jet1Eta==0 && mismatch_filters==0 && mismatch_Jet1Pt==0)
	  h_FirstMismatch_HLT1_RECO0->Fill(4);
	if(mismatch_Jet2Eta==1 && mismatch_Jet2Pt==0 && mismatch_Jet1ID==0 && mismatch_Jet1Eta==0 && mismatch_filters==0 && mismatch_Jet1Pt==0)
	  h_FirstMismatch_HLT1_RECO0->Fill(5);
	if(mismatch_Jet2ID==1 && mismatch_Jet2Eta==0 && mismatch_Jet2Pt==0 && mismatch_Jet1ID==0 && mismatch_Jet1Eta==0 && mismatch_filters==0 && mismatch_Jet1Pt==0)
	  h_FirstMismatch_HLT1_RECO0->Fill(6);
	if(mismatch_DPhi==1 && mismatch_Jet2ID==0 && mismatch_Jet2Eta==0 && mismatch_Jet2Pt==0 && mismatch_Jet1ID==0 && mismatch_Jet1Eta==0 && mismatch_filters==0 && mismatch_Jet1Pt==0)
	  h_FirstMismatch_HLT1_RECO0->Fill(7);
	if(mismatch_DEta==1 && mismatch_DPhi==0 && mismatch_Jet2ID==0 && mismatch_Jet2Eta==0 && mismatch_Jet2Pt==0 && mismatch_Jet1ID==0 && mismatch_Jet1Eta==0 && mismatch_filters==0 && mismatch_Jet1Pt==0)
	  {
	    h_FirstMismatch_HLT1_RECO0->Fill(8);
	    if(debugEv)
	      {
		cout << "SelectionHLT1RECO0 -  Deta : RECO , HLT = " 
		     << getVariableValue("RECO_DEta_PFJet1PFJet2") << " , " << getVariableValue("HLT_DEta_PFJet1PFJet2") 
		     << endl;
		
		CreateAndFillUserTH1D("h_DEtaBias_PFJet1PFJet2_HLT1_RECO0_mismatch_deta", 100, 0, 10, fabs(getVariableValue("HLT_DEta_PFJet1PFJet2") - getVariableValue("RECO_DEta_PFJet1PFJet2")) );
		CreateAndFillUserTH1D("h_HLT_N_PFJet_HLT1_RECO0_mismatch_deta", 51, -0.5, 50.5, v_idx_hltpfjet_PtEtaCut_ID.size() );

// 		cout << "*** HLT ***" << endl;
// 		for(int ijet=0; ijet<min(5,int(HLTPFJetCorrPt->size())); ijet++)
// 		  {
// 		    cout << ijet << " - PT,ETA,PHI: " <<  HLTPFJetCorrPt->at(ijet) << " , " <<  HLTPFJetEta->at(ijet) << " , " << HLTPFJetPhi->at(ijet) << endl; 
// 		  }
// 		cout << "*** RECO ***" << endl;
// 		for(int ijet=0; ijet<min(5,int(PFJetPt->size())); ijet++)
// 		  {
// 		    cout << ijet << " - PT,ETA,PHI: " <<  PFJetPt->at(ijet) << " , " <<  PFJetEta->at(ijet) << " , " << PFJetPhi->at(ijet) << endl; 
// 		  }

	      }
	  }

	FillUserTH1D("h_HLT_M_PFJet1PFJet2_HLT1_RECO0", getVariableValue("HLT_M_PFJet1PFJet2") );
	FillUserTH1D("h_RECO_M_PFJet1PFJet2_HLT1_RECO0", getVariableValue("RECO_M_PFJet1PFJet2") );
	FillUserTH1D("h_HLT_DPhi_PFJet1PFJet2_HLT1_RECO0", getVariableValue("HLT_DPhi_PFJet1PFJet2") );
	FillUserTH1D("h_RECO_DPhi_PFJet1PFJet2_HLT1_RECO0", getVariableValue("RECO_DPhi_PFJet1PFJet2") );
      }

    if( pass_HLT_selection == 0 && pass_RECO_selection == 1 )
      {	
	int mismatch_filters = 0;
	int mismatch_Jet1Pt = 0;
	int mismatch_Jet1Eta = 0;
	int mismatch_Jet1ID = 0;
	int mismatch_Jet2Pt = 0;
	int mismatch_Jet2Eta = 0;
	int mismatch_Jet2ID = 0;
	int mismatch_DPhi = 0;
	int mismatch_DEta = 0;
	
	if( pass_RECO_filters != pass_HLT_filters )
	  mismatch_filters = 1;
	
	if( passedCut("HLT_PFJet1_Pt") != passedCut("RECO_PFJet1_Pt") )
	  mismatch_Jet1Pt = 1;
	
	if( passedCut("HLT_PFJet1_Eta") != passedCut("RECO_PFJet1_Eta") )
	  mismatch_Jet1Eta = 1;

	if( passedCut("HLT_PFJet1_PassJetID") != passedCut("RECO_PFJet1_PassJetID") )
	  mismatch_Jet1ID = 1;
	
	if( passedCut("HLT_PFJet2_Pt") != passedCut("RECO_PFJet2_Pt") )
	  mismatch_Jet2Pt = 1;
	
	if( passedCut("HLT_PFJet2_Eta") != passedCut("RECO_PFJet2_Eta") )
	  mismatch_Jet2Eta = 1;

	if( passedCut("HLT_PFJet2_PassJetID") != passedCut("RECO_PFJet2_PassJetID") )
	  mismatch_Jet2ID = 1;

	if( passedCut("HLT_DPhi_PFJet1PFJet2") != passedCut("RECO_DPhi_PFJet1PFJet2") )
	  mismatch_DPhi = 1;
	
	if( passedCut("HLT_DEta_PFJet1PFJet2") != passedCut("RECO_DEta_PFJet1PFJet2") )
	  mismatch_DEta = 1;

	if(debugEv)
	  {
	    cout << "SelectionHLT0RECO1 ----------" << endl;
	    cout << "SelectionHLT0RECO1 - run: " << run << " , ls: " << ls << " , event: " << event << endl;
	  }

	if(mismatch_filters==1)
	  {
	    h_FirstMismatch_HLT0_RECO1->Fill(0);

	    if(debugEv)
	      {
		cout << "SelectionHLT0RECO1 - PassHLTPrimaryVertex: " << getVariableValue("PassHLTPrimaryVertex") << endl;
		cout << "SelectionHLT0RECO1 - PassRECOPrimaryVertex: " << getVariableValue("PassRECOPrimaryVertex") << endl;
		cout << "SelectionHLT0RECO1 - PassBPTX0: " << getVariableValue("PassBPTX0") << endl;
		cout << "SelectionHLT0RECO1 - PassBeamScraping: " << getVariableValue("PassBeamScraping") << endl;
		cout << "SelectionHLT0RECO1 - PassHBHENoiseFilter: " << getVariableValue("PassHBHENoiseFilter") << endl;
		cout << "SelectionHLT0RECO1 - PassBeamHaloFilterTight: " << getVariableValue("PassBeamHaloFilterTight") << endl;
		cout << "SelectionHLT0RECO1 - PassTrackingFailure: " << getVariableValue("PassTrackingFailure") << endl;
	      }
	  }
	if(mismatch_Jet1Pt==1 && mismatch_filters==0)
	  h_FirstMismatch_HLT0_RECO1->Fill(1);
	if(mismatch_Jet1Eta==1 && mismatch_filters==0 && mismatch_Jet1Pt==0)
	  h_FirstMismatch_HLT0_RECO1->Fill(2);
	if(mismatch_Jet1ID==1 && mismatch_Jet1Eta==0 && mismatch_filters==0 && mismatch_Jet1Pt==0)
	  h_FirstMismatch_HLT0_RECO1->Fill(3);
	if(mismatch_Jet2Pt==1 && mismatch_Jet1ID==0 && mismatch_Jet1Eta==0 && mismatch_filters==0 && mismatch_Jet1Pt==0)
	  h_FirstMismatch_HLT0_RECO1->Fill(4);
	if(mismatch_Jet2Eta==1 && mismatch_Jet2Pt==0 && mismatch_Jet1ID==0 && mismatch_Jet1Eta==0 && mismatch_filters==0 && mismatch_Jet1Pt==0)
	  h_FirstMismatch_HLT0_RECO1->Fill(5);
	if(mismatch_Jet2ID==1 && mismatch_Jet2Eta==0 && mismatch_Jet2Pt==0 && mismatch_Jet1ID==0 && mismatch_Jet1Eta==0 && mismatch_filters==0 && mismatch_Jet1Pt==0)
	  h_FirstMismatch_HLT0_RECO1->Fill(6);
	if(mismatch_DPhi==1 && mismatch_Jet2ID==0 && mismatch_Jet2Eta==0 && mismatch_Jet2Pt==0 && mismatch_Jet1ID==0 && mismatch_Jet1Eta==0 && mismatch_filters==0 && mismatch_Jet1Pt==0)
	  h_FirstMismatch_HLT0_RECO1->Fill(7);
	if(mismatch_DEta==1 && mismatch_DPhi==0 && mismatch_Jet2ID==0 && mismatch_Jet2Eta==0 && mismatch_Jet2Pt==0 && mismatch_Jet1ID==0 && mismatch_Jet1Eta==0 && mismatch_filters==0 && mismatch_Jet1Pt==0)
	  {
	    h_FirstMismatch_HLT0_RECO1->Fill(8);
	    if(debugEv)
	      {
		cout << "SelectionHLT0RECO1 -  Deta : RECO , HLT = " 
		     << getVariableValue("RECO_DEta_PFJet1PFJet2") << " , " << getVariableValue("HLT_DEta_PFJet1PFJet2") 
		     << endl;
	      }
	  }

	FillUserTH1D("h_HLT_M_PFJet1PFJet2_HLT0_RECO1", getVariableValue("HLT_M_PFJet1PFJet2") );
	FillUserTH1D("h_RECO_M_PFJet1PFJet2_HLT0_RECO1", getVariableValue("RECO_M_PFJet1PFJet2") );
	FillUserTH1D("h_HLT_DPhi_PFJet1PFJet2_HLT0_RECO1", getVariableValue("HLT_DPhi_PFJet1PFJet2") );
	FillUserTH1D("h_RECO_DPhi_PFJet1PFJet2_HLT0_RECO1", getVariableValue("RECO_DPhi_PFJet1PFJet2") );
      }



    //INFO
    //      // retrieve value of previously filled variables (after making sure that they were filled)
    //      double totpTEle;
    //      if ( variableIsFilled("pT1stEle") && variableIsFilled("pT2ndEle") )
    //        totpTEle = getVariableValue("pT1stEle")+getVariableValue("pT2ndEle");
    //      // reject events that did not pass level 0 cuts
    //      if( !passedCut("0") ) continue;
    //      // ......

    //FillUserTH1D("h1_MTenu_PAS_minus", getVariableValue("MTenu_PAS"));
    //FillUserTH2D("h2_phi_VS_eta_1stEle", getVariableValue("Eta1stEle_PAS"), getVariableValue("Phi1stEle_PAS") );
    // 	    CreateAndFillUserTH1D("h1_ElectronRelIso_highMT", 1000, 0, 1, ElectronRelIso->at(myEle) );
    //if( passedAllPreviousCuts("d1_DPhi_METe_METj")
    //    && passedCut("nMuon_PtCut_IDISO")
    
    // Produce skim 
    //     if( //passedAllPreviousCuts("PassPrimaryVertex") 
    //  	passedCut("nPFJet")
    // 	&& passedCut("PFJet1_Pt")
    // 	&& passedCut("PFJet2_Pt")
    // 	) 
    fillSkimTree();
    
    // Produce reduced skim
    //     if( //passedAllPreviousCuts("PassPrimaryVertex") 
    //  	passedCut("nPFJet")
    // 	&& passedCut("PFJet1_Pt")
    // 	&& passedCut("PFJet2_Pt")
    // 	) 
    fillReducedSkimTree();
    
    ////////////////////// User's code to be done for every event - END ///////////////////////

  } // End of loop over events


  ////////////////////// User's code to write histos - BEGIN ///////////////////////

  h_FirstMismatch_HLT1_RECO0->Write();
  h_FirstMismatch_HLT0_RECO1->Write();

  p_MassBias_MATCH_vs_Mreco->Write();

  if(applyPFJEC==0)
    {
      p_PTJetsMedian_vs_NVtx->Write();
      h_PtPFJetsInEvent1Ev->Write();
    }

  delete p_PTJetsMedian_vs_NVtx;
  delete h_PtPFJetsInEvent1Ev;
  delete h_FirstMismatch_HLT1_RECO0;
  delete h_FirstMismatch_HLT0_RECO1;
  delete p_MassBias_MATCH_vs_Mreco;

  ////////////////////// User's code to write histos - END ///////////////////////

  //STDOUT("analysisClass::Loop() ends");
}
