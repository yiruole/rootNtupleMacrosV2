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


//-----------------------------
//### JetID ### --> see https://twiki.cern.ch/twiki/bin/view/CMS/ExoticaHighPtJets#JetId

// bool JetIdloose(double ak5JetJIDresEMF,double ak5JetJIDfHPD,int ak5JetJIDn90Hits, double ak5JetEta){
//   bool jetidloose=false;
//   bool jetidresEMF=true;
//
//   double fhpdmax = 0.98;
//   double n90hitsmin =1;
//   double emf_min = 0.01;
//
//   if(fabs(ak5JetEta)<2.6 && ak5JetJIDresEMF<=emf_min) jetidresEMF=false;
//
//   if(jetidresEMF && ak5JetJIDfHPD<fhpdmax && ak5JetJIDn90Hits>n90hitsmin) {
//     jetidloose=true;
//   }
//   return jetidloose;
// }
//
// bool JetIdtight(double ak5JetJIDresEMF,double ak5JetJIDfHPD,int ak5JetJIDn90Hits, double ak5JetEta, double ak5JetPt){
//   bool jetidtight=false;
//   bool jetidresEMF=true;
//   bool jetidfHPD_highPt=true;
//
//   double fhpdmax = 0.98;
//   double n90hitsmin =1;
//   double emf_min = 0.01;
//
//   if(fabs(ak5JetEta)<2.6 && ak5JetJIDresEMF<=emf_min) jetidresEMF=false;
//   if(fabs(ak5JetEta)<2.6 && ak5JetPt>80 && ak5JetJIDresEMF>=1) jetidresEMF=false;
//   if(ak5JetPt>25 && ak5JetJIDfHPD>=0.95) jetidfHPD_highPt=false;
//
//   if(jetidresEMF && jetidfHPD_highPt && ak5JetJIDfHPD<fhpdmax && ak5JetJIDn90Hits>n90hitsmin)
//     {
//       jetidtight=true;
//     }
//   return jetidtight;
// }
//
// bool PFJetIdloose(const double ak5ChargedHadronFraction, const double ak5ChargedEmFraction, const double ak5NeutralHadronFraction, const double ak5NeutralEmFraction, const double ak5JetEta){
//   bool jetidloose=false;
//   bool jetidChFrac=true;
//
//   double chHadFrac = 0.;
//   double chEmFrac = 0.99;
//   double neutHadFrac = 0.99;
//   double neutEmFrac = 0.99;
//
//   if(fabs(ak5JetEta)<2.4 && ak5ChargedHadronFraction<=chHadFrac && ak5ChargedEmFraction>=chEmFrac) jetidChFrac=false;
//
//   if(jetidChFrac && ak5NeutralHadronFraction<neutHadFrac && ak5NeutralEmFraction<neutEmFrac) {
//     jetidloose=true;
//   }
//   return jetidloose;
// }

//-----------------------------

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

  double ele_PtCut =  getPreCutValue1("ele_PtCut");
  double eleEta_bar = getPreCutValue1("eleEta_bar");
  double eleEta_end_min = getPreCutValue1("eleEta_end");
  double eleEta_end_max = getPreCutValue2("eleEta_end");

  double jet_PtCut =    getPreCutValue1("jet_PtCut");
  double jet_EtaCut = getPreCutValue1("jet_EtaCut");
  double jet_TCHELCut = getPreCutValue1("jet_TCHELCut");
  double jet_ele_DeltaRcut =   getPreCutValue1("jet_ele_DeltaRcut");
  double jet_PtCut_forMetScale =    getPreCutValue1("jet_PtCut_forMetScale");

  double EleEnergyScale_EB=getPreCutValue1("EleEnergyScale_EB");
  double EleEnergyScale_EE=getPreCutValue1("EleEnergyScale_EE");
  double JetEnergyScale=getPreCutValue1("JetEnergyScale");

  // Not used when using ElectronHeepID and heepBitMask // int eleIDType = (int) getPreCutValue1("eleIDType");
  int heepBitMask_EB  =  getPreCutValue1("heepBitMask_EBGapEE") ;
  int heepBitMask_GAP =  getPreCutValue2("heepBitMask_EBGapEE") ;
  int heepBitMask_EE  =  getPreCutValue3("heepBitMask_EBGapEE") ;

  double muon_PtCut = getPreCutValue1("muon_PtCut");
  double muFidRegion = getPreCutValue1("muFidRegion"); // currently unused !!!
  double muNHits_minThresh = getPreCutValue1("muNHits_minThresh");
  double muTrkD0Maximum = getPreCutValue1("muTrkD0Maximum");

  int doGenLevelStudies = getPreCutValue1("doGenLevelStudies");

  int jetAlgorithm = getPreCutValue1("jetAlgorithm");
  int metAlgorithm = getPreCutValue1("metAlgorithm");

  double vertexMinimumNDOF = getPreCutValue1("vertexMinimumNDOF");
  double vertexMaxAbsZ = getPreCutValue1("vertexMaxAbsZ");
  double vertexMaxd0 = getPreCutValue1("vertexMaxd0");

  double met_Thresh = getPreCutValue1("met_Thresh");
  double minMETPt1stEle_Thresh = getPreCutValue1("minMETPt1stEle_Thresh");
  double Pt1stEle_PAS_Thresh = getPreCutValue1("Pt1stEle_PAS_Thresh");
  double Pt1stEle_PAS_Thresh2 = getPreCutValue2("Pt1stEle_PAS_Thresh");
  double Pt1stJet_PAS_Thresh = getPreCutValue1("Pt1stJet_PAS_Thresh");
  double Pt2ndJet_PAS_Thresh = getPreCutValue1("Pt2ndJet_PAS_Thresh");
  double MTenu_Thresh = getPreCutValue1("MTenu_Thresh");
  double MTenu_Thresh2 = getPreCutValue2("MTenu_Thresh");
  double sT_Thresh = getPreCutValue1("sT_Thresh");
  double Mej_Thresh = getPreCutValue1("Mej_Thresh");

  int plotEleIDIsoVar = getPreCutValue1("plotEleIDIsoVar");

  int doPlot_Wmore0jet = getPreCutValue1("doPlot_Wmore0jet");
  int doPlot_Wmore1jet = getPreCutValue1("doPlot_Wmore1jet");

  int doExtraChecks = getPreCutValue1("doExtraChecks");

  int doPUMETSmearing = getPreCutValue1("doPUMETSmearing");
  double METxySigmaPerPU = getPreCutValue1("METxySigmaPerPU");

  ////////////////////// User's code to get preCut values - END /////////////////

  ////////////////////// User's code to book histos - BEGIN ///////////////////////

  if(doGenLevelStudies)
    {
      CreateUserTH1D("h_num_Neutrinos", 10,0,10);
      CreateUserTH1D("h_num_Ws", 10,0,10);
      CreateUserTH2D("h2_pfMET_vs_neutrinoPt", 100,0,1000, 100, 0, 1000);
      CreateUserTH1D("h_DeltaPhi_pfMET_neutrino", 100,0,3.1416);
      CreateUserTH1D("h_WsPt", 100,0,1000);
      CreateUserTH2D("h2_pfMET_vs_neutrinoPt__sT", 100,0,1000, 100, 0, 1000);
      CreateUserTH2D("h2_pfMET_vs_neutrinoPt_DphiMETeSmall__sT", 100,0,1000, 100, 0, 1000);
      CreateUserTH2D("h2_pfMET_vs_neutrinoPt_d2Small__sT", 100,0,1000, 100, 0, 1000);
      CreateUserTH2D("h2_pfMET_vs_neutrinoPt_d2DphiMETeLarge__sT", 100,0,1000, 100, 0, 1000);
      CreateUserTH1D("h_DeltaPhi_pfMET_neutrino__sT", 100,0,3.1416);
      CreateUserTH1D("h_DeltaPhi_pfMET_neutrino_DphiMETeSmall__sT", 100,0,3.1416);
      CreateUserTH1D("h_DeltaPhi_pfMET_neutrino_d2Small__sT", 100,0,3.1416);
      CreateUserTH1D("h_DeltaPhi_pfMET_neutrino_d2DphiMETeLarge__sT", 100,0,3.1416);
      CreateUserTH1D("h_WsPt__sT", 100,0,1000);
      CreateUserTH1D("h_WsPt_DphiMETeSmall__sT", 100,0,1000);
      CreateUserTH1D("h_WsPt_d2Small__sT", 100,0,1000);
      CreateUserTH1D("h_WsPt_d2DphiMETeLarge__sT", 100,0,1000);
      CreateUserTH1D("h1_LQGenEle_Pt", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
      CreateUserTH1D("h1_LQGenEle_Eta", getHistoNBins("Eta1stEle_PAS"), getHistoMin("Eta1stEle_PAS"), getHistoMax("Eta1stEle_PAS") );
      CreateUserTH1D("h1_LQGenEle_Phi", getHistoNBins("Phi1stEle_PAS"), getHistoMin("Phi1stEle_PAS"), getHistoMax("Phi1stEle_PAS") );
      CreateUserTH1D("h1_GenEle_Pt", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
      CreateUserTH1D("h1_GenEle_Eta", getHistoNBins("Eta1stEle_PAS"), getHistoMin("Eta1stEle_PAS"), getHistoMax("Eta1stEle_PAS") );
      CreateUserTH1D("h1_GenEle_Phi", getHistoNBins("Phi1stEle_PAS"), getHistoMin("Phi1stEle_PAS"), getHistoMax("Phi1stEle_PAS") );
      CreateUserTH1D("h1_PtReco_over_PtGen_Ele", 200, 0, 2 );
      CreateUserTH2D("h2_PtReco_over_PtGen_vs_GenEle_Pt"
		     , getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS")
		     , 200, 0, 2
		     );
      CreateUserTH1D("h1_deltaR_RecoEle_GenEle", 100, 0, 5  );
    }

  CreateUserTH2D("h2_MTnuj_vs_MET", 200,0,1000,200,0,1000);
  CreateUserTH2D("h2_ST_vs_MET", 200,0,1000,200,0,2000);
  CreateUserTH2D("h2_ST_vs_MTnuj", 200,0,1000,200,0,2000);
  CreateUserTH2D("h2_DeltaPhiMETEle_vs_MET", 200,0,1000,30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMET1stJet_vs_MET", 200,0,1000,30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMET2ndJet_vs_MET", 200,0,1000,30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet", 30, 0, 3.1416,30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet_plus", 30, 0, 3.1416,30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet_minus", 30, 0, 3.1416,30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMETEle_vs_MET2ndJet", 30, 0, 3.1416, 30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMET2ndJet_vs_MET1stJet", 30, 0, 3.1416, 30, 0, 3.1416);
  CreateUserTH2D("h2_MT_vs_etaEle", 100,-5,5,200,0,1000);
  CreateUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__minMETpTe", 30, 0, 3.1416, 30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__minMETpTe_plus", 30, 0, 3.1416,30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__minMETpTe_minus", 30, 0, 3.1416,30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__MTenu", 30, 0, 3.1416, 30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__MTenu_plus", 30, 0, 3.1416,30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__MTenu_minus", 30, 0, 3.1416,30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__sT", 30, 0, 3.1416, 30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__sT_plus", 30, 0, 3.1416,30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__sT_minus", 30, 0, 3.1416,30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMETEle_vs_MET2ndJet__sT", 30, 0, 3.1416, 30, 0, 3.1416);
  CreateUserTH2D("h2_DeltaPhiMET2ndJet_vs_MET1stJet__sT", 30, 0, 3.1416, 30, 0, 3.1416);

  CreateUserTH1D("h1_MTenu_PAS_plus", getHistoNBins("MTenu_PAS"), getHistoMin("MTenu_PAS"), getHistoMax("MTenu_PAS"));
  CreateUserTH1D("h1_MTenu_PAS_minus", getHistoNBins("MTenu_PAS"), getHistoMin("MTenu_PAS"), getHistoMax("MTenu_PAS"));
  CreateUserTH1D("h1_sT_PAS_plus", getHistoNBins("sT_PAS"), getHistoMin("sT_PAS"), getHistoMax("sT_PAS"));
  CreateUserTH1D("h1_sT_PAS_minus", getHistoNBins("sT_PAS"), getHistoMin("sT_PAS"), getHistoMax("sT_PAS"));
  CreateUserTH1D("h1_sT_AllNonSTCuts", getHistoNBins("sT_PAS"), getHistoMin("sT_PAS"), getHistoMax("sT_PAS"));

  CreateUserTH1D("h1_MTenu_PAS_DeltaPhiMET2ndJet_0_1", getHistoNBins("MTenu_PAS"), getHistoMin("MTenu_PAS"), getHistoMax("MTenu_PAS"));
  CreateUserTH1D("h1_MTenu_PAS_DeltaPhiMET2ndJet_1_2", getHistoNBins("MTenu_PAS"), getHistoMin("MTenu_PAS"), getHistoMax("MTenu_PAS"));
  CreateUserTH1D("h1_MTenu_PAS_DeltaPhiMET2ndJet_2_pi", getHistoNBins("MTenu_PAS"), getHistoMin("MTenu_PAS"), getHistoMax("MTenu_PAS"));

  CreateUserTH2D("h2_phi_VS_eta_1stEle", 100, -5, 5, 60, -3.1416, 3.1416);
  CreateUserTH2D("h2_phi_VS_eta_1stJet", 100, -5, 5, 60, -3.1416, 3.1416);
  CreateUserTH2D("h2_phi_VS_eta_2ndJet", 100, -5, 5, 60, -3.1416, 3.1416);

  CreateUserTH1D("h1_Phi1stEle_PAS_EleBarrel", getHistoNBins("Phi1stEle_PAS"), getHistoMin("Phi1stEle_PAS"), getHistoMax("Phi1stEle_PAS"));
  CreateUserTH1D("h1_Phi1stEle_PAS_EleEndcap", getHistoNBins("Phi1stEle_PAS"), getHistoMin("Phi1stEle_PAS"), getHistoMax("Phi1stEle_PAS"));
  CreateUserTH1D("h1_METPhi_PAS_EleBarrel", getHistoNBins("METPhi_PAS"), getHistoMin("METPhi_PAS"), getHistoMax("METPhi_PAS"));
  CreateUserTH1D("h1_METPhi_PAS_EleEndcap", getHistoNBins("METPhi_PAS"), getHistoMin("METPhi_PAS"), getHistoMax("METPhi_PAS"));

  CreateUserTH2D("h2_Phi1stEle_vs_METPhi", 60, -3.1416, 3.1416 , 60, -3.1416, 3.1416 );
  CreateUserTH2D("h2_Phi1stEle_vs_PtEleOverST", 10, 0, 1 , 60, -3.1416, 3.1416 );
  CreateUserTH2D("h2_METPhi_vs_PtEleOverST", 10, 0, 1 , 60, -3.1416, 3.1416 );

  CreateUserTH1D("h1_Charge1stEle_PAS__sT", getHistoNBins("Charge1stEle_PAS"), getHistoMin("Charge1stEle_PAS"), getHistoMax("Charge1stEle_PAS"));
  CreateUserTH1D("h1_Eta1stEle_PAS__sT", getHistoNBins("Eta1stEle_PAS"), getHistoMin("Eta1stEle_PAS"), getHistoMax("Eta1stEle_PAS"));

  CreateUserTH1D("h1_ElectronPt_all", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
  CreateUserTH1D("h1_ElectronEta_all", getHistoNBins("Eta1stEle_PAS"), getHistoMin("Eta1stEle_PAS"), getHistoMax("Eta1stEle_PAS") );
  CreateUserTH1D("h1_ElectronPhi_all", getHistoNBins("Phi1stEle_PAS"), getHistoMin("Phi1stEle_PAS"), getHistoMax("Phi1stEle_PAS") );

  if( plotEleIDIsoVar )
    {
      CreateUserTH1D("h1_ElectronDeltaEtaTrkSC_all", 200,-0.05,0.05 );
      CreateUserTH1D("h1_ElectronDeltaPhiTrkSC_all", 200,-0.5,0.5 );
      CreateUserTH1D("h1_ElectronHoE_all", 75,0,0.15 );
      CreateUserTH1D("h1_ElectronSigmaIEtaIEta_all", 100,0,0.1 );
      CreateUserTH1D("h1_ElectronEcalHcalIsoHeep_all", 500,0,100 );
      CreateUserTH1D("h1_ElectronEcalIsoHeep_all", 500,0,100 );
      CreateUserTH1D("h1_ElectronHcalIsoD1Heep_all", 500,0,100 );
      CreateUserTH1D("h1_ElectronHcalIsoD2Heep_all", 200,0,100 );
      CreateUserTH1D("h1_ElectronTrkIsoHeep_all", 200,0,100 );
      CreateUserTH1D("h1_ElectronE2x5OverE5x5_all", 100,0,1 );
      CreateUserTH1D("h1_ElectronE1x5OverE5x5_all", 100,0,1 );

      CreateUserTH1D("h1_ElectronDeltaEtaTrkSC_highMT", 200,-0.05,0.05 );
      CreateUserTH1D("h1_ElectronDeltaPhiTrkSC_highMT", 200,-0.5,0.5 );
      CreateUserTH1D("h1_ElectronHoE_highMT", 75,0,0.15 );
      CreateUserTH1D("h1_ElectronSigmaIEtaIEta_highMT", 100,0,0.1 );
      CreateUserTH1D("h1_ElectronEcalHcalIsoHeep_highMT", 500,0,100 );
      CreateUserTH1D("h1_ElectronHcalIsoD2Heep_highMT", 200,0,100 );
      CreateUserTH1D("h1_ElectronTrkIsoHeep_highMT", 200,0,100 );
      CreateUserTH1D("h1_ElectronE2x5OverE5x5_highMT", 100,0,1 );
      CreateUserTH1D("h1_ElectronE1x5OverE5x5_highMT", 100,0,1 );

      CreateUserTH1D("h1_ElectronDeltaEtaTrkSC_barrel_highMT", 200,-0.05,0.05 );
      CreateUserTH1D("h1_ElectronDeltaPhiTrkSC_barrel_highMT", 200,-0.5,0.5 );
      CreateUserTH1D("h1_ElectronHoE_barrel_highMT", 75,0,0.15 );
      CreateUserTH1D("h1_ElectronSigmaIEtaIEta_barrel_highMT", 100,0,0.1 );
      CreateUserTH1D("h1_ElectronEcalHcalIsoHeep_barrel_highMT", 500,0,100 );
      CreateUserTH1D("h1_ElectronHcalIsoD2Heep_barrel_highMT", 200,0,100 );
      CreateUserTH1D("h1_ElectronTrkIsoHeep_barrel_highMT", 200,0,100 );
      CreateUserTH1D("h1_ElectronE2x5OverE5x5_barrel_highMT", 100,0,1 );
      CreateUserTH1D("h1_ElectronE1x5OverE5x5_barrel_highMT", 100,0,1 );

      CreateUserTH1D("h1_ElectronDeltaEtaTrkSC_endcap_highMT", 200,-0.05,0.05 );
      CreateUserTH1D("h1_ElectronDeltaPhiTrkSC_endcap_highMT", 200,-0.5,0.5 );
      CreateUserTH1D("h1_ElectronHoE_endcap_highMT", 75,0,0.15 );
      CreateUserTH1D("h1_ElectronSigmaIEtaIEta_endcap_highMT", 100,0,0.1 );
      CreateUserTH1D("h1_ElectronEcalHcalIsoHeep_endcap_highMT", 500,0,100 );
      CreateUserTH1D("h1_ElectronHcalIsoD2Heep_endcap_highMT", 200,0,100 );
      CreateUserTH1D("h1_ElectronTrkIsoHeep_endcap_highMT", 200,0,100 );
      CreateUserTH1D("h1_ElectronE2x5OverE5x5_endcap_highMT", 100,0,1 );
      CreateUserTH1D("h1_ElectronE1x5OverE5x5_endcap_highMT", 100,0,1 );
    }

  if(doPlot_Wmore0jet)
    {
      CreateUserTH1D("h1_Pt1stEle_W0jet", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"));
      CreateUserTH1D("h1_Eta1stEle_W0jet", getHistoNBins("Eta1stEle_PAS"), getHistoMin("Eta1stEle_PAS"), getHistoMax("Eta1stEle_PAS"));
      CreateUserTH1D("h1_Phi1stEle_W0jet", getHistoNBins("Phi1stEle_PAS"), getHistoMin("Phi1stEle_PAS"), getHistoMax("Phi1stEle_PAS"));
      CreateUserTH1D("h1_MET_W0jet", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS"));
      CreateUserTH1D("h1_METPhi_W0jet", getHistoNBins("METPhi_PAS"), getHistoMin("METPhi_PAS"), getHistoMax("METPhi_PAS"));
      CreateUserTH1D("h1_MTenu_W0jet", getHistoNBins("MTenu_PAS"), getHistoMin("MTenu_PAS"), getHistoMax("MTenu_PAS"));
      CreateUserTH1D("h1_mDeltaPhiMETEle_W0jet", getHistoNBins("mDeltaPhiMETEle"), getHistoMin("mDeltaPhiMETEle"), getHistoMax("mDeltaPhiMETEle"));
      CreateUserTH2D("h2_EToverPT_vs_ET_1stEle_W0jet"
		     , getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS")
		     , 200, 0, 20
		     );
      CreateUserTH1D("h1_Pt1stEle_W0jet_barrel", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"));
      CreateUserTH1D("h1_Pt1stEle_W0jet_endcap", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"));
      CreateUserTH1D("h1_Energy1stEle_W0jet_barrel", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"));
      CreateUserTH1D("h1_Energy1stEle_W0jet_endcap", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"));
      CreateUserTH1D("h1_TrackPt1stEle_W0jet_barrel", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"));
      CreateUserTH1D("h1_TrackPt1stEle_W0jet_endcap", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"));
    }

  if(doPlot_Wmore1jet)
    {
      CreateUserTH1D("h1_Pt1stEle_W1jet", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"));
      CreateUserTH1D("h1_Eta1stEle_W1jet", getHistoNBins("Eta1stEle_PAS"), getHistoMin("Eta1stEle_PAS"), getHistoMax("Eta1stEle_PAS"));
      CreateUserTH1D("h1_Phi1stEle_W1jet", getHistoNBins("Phi1stEle_PAS"), getHistoMin("Phi1stEle_PAS"), getHistoMax("Phi1stEle_PAS"));
      CreateUserTH1D("h1_MET_W1jet", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS"));
      CreateUserTH1D("h1_METPhi_W1jet", getHistoNBins("METPhi_PAS"), getHistoMin("METPhi_PAS"), getHistoMax("METPhi_PAS"));
      CreateUserTH1D("h1_MTenu_W1jet", getHistoNBins("MTenu_PAS"), getHistoMin("MTenu_PAS"), getHistoMax("MTenu_PAS"));
      CreateUserTH1D("h1_mDeltaPhiMETEle_W1jet", getHistoNBins("mDeltaPhiMETEle"), getHistoMin("mDeltaPhiMETEle"), getHistoMax("mDeltaPhiMETEle"));
      CreateUserTH2D("h2_EToverPT_vs_ET_1stEle_W1jet"
		     , getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS")
		     , 200, 0, 20
		     );
      CreateUserTH1D("h1_Pt1stEle_W1jet_barrel", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"));
      CreateUserTH1D("h1_Pt1stEle_W1jet_endcap", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"));
      CreateUserTH1D("h1_Energy1stEle_W1jet_barrel", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"));
      CreateUserTH1D("h1_Energy1stEle_W1jet_endcap", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"));
      CreateUserTH1D("h1_TrackPt1stEle_W1jet_barrel", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"));
      CreateUserTH1D("h1_TrackPt1stEle_W1jet_endcap", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"));
    }


  if(doExtraChecks)
    {
      CreateUserTH2D("h2_Ngenlept_vs_Njets_highMT",
		     getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID"),
		     getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID")
		     );
    }
  TH1F* h1_Ngenlept__ee_emu_etau_eonly_noe_others = new TH1F("h1_Ngenlept__ee_emu_etau_eonly_noe_others", "h1_Ngenlept__ee_emu_etau_eonly_noe_others", 6, 0, 6);
  h1_Ngenlept__ee_emu_etau_eonly_noe_others->Sumw2();
  h1_Ngenlept__ee_emu_etau_eonly_noe_others->GetXaxis()->SetBinLabel(1,"ee");
  h1_Ngenlept__ee_emu_etau_eonly_noe_others->GetXaxis()->SetBinLabel(2,"e#mu");
  h1_Ngenlept__ee_emu_etau_eonly_noe_others->GetXaxis()->SetBinLabel(3,"e#tau");
  h1_Ngenlept__ee_emu_etau_eonly_noe_others->GetXaxis()->SetBinLabel(4,"e");
  h1_Ngenlept__ee_emu_etau_eonly_noe_others->GetXaxis()->SetBinLabel(5,"no e");
  h1_Ngenlept__ee_emu_etau_eonly_noe_others->GetXaxis()->SetBinLabel(6,"others");

  if(doExtraChecks)
    {
      CreateUserTH2D("h2_minDeltaPhiMETj_vs_minDeltaRej_highMT",
		     getHistoNBins("minDRej"), getHistoMin("minDRej"), getHistoMax("minDRej"),
		     getHistoNBins("mDeltaPhiMET1stJet_PAS"), getHistoMin("mDeltaPhiMET1stJet_PAS"), getHistoMax("mDeltaPhiMET1stJet_PAS")
		     );
      CreateUserTH1D("h1_Njet_highPt1stEle", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_NjetTCHELBTag_highPt1stEle", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_minDRej_highPt1stEle", getHistoNBins("minDRej"), getHistoMin("minDRej"), getHistoMax("minDRej") );
      CreateUserTH1D("h1_mDeltaPhiMETEle_highPt1stEle", getHistoNBins("mDeltaPhiMETEle"), getHistoMin("mDeltaPhiMETEle"), getHistoMax("mDeltaPhiMETEle") );
      CreateUserTH1D("h1_mDeltaPhiMET1stJet_highPt1stEle", getHistoNBins("mDeltaPhiMET1stJet"), getHistoMin("mDeltaPhiMET1stJet"), getHistoMax("mDeltaPhiMET1stJet") );
      CreateUserTH1D("h1_mDeltaPhiMET2ndJet_highPt1stEle", getHistoNBins("mDeltaPhiMET2ndJet"), getHistoMin("mDeltaPhiMET2ndJet"), getHistoMax("mDeltaPhiMET2ndJet") );
      CreateUserTH1D("h1_DeltaRjets_PAS_highPt1stEle", getHistoNBins("DeltaRjets_PAS"), getHistoMin("DeltaRjets_PAS"), getHistoMax("DeltaRjets_PAS") );
      CreateUserTH1D("h1_Vtxd01stEle_PAS_highPt1stEle", getHistoNBins("Vtxd01stEle_PAS"), getHistoMin("Vtxd01stEle_PAS"), getHistoMax("Vtxd01stEle_PAS") );
      CreateUserTH1D("h1_MissingHits1stEle_PAS_highPt1stEle", getHistoNBins("MissingHits1stEle_PAS"), getHistoMin("MissingHits1stEle_PAS"), getHistoMax("MissingHits1stEle_PAS") );
      CreateUserTH1D("h1_Dist1stEle_PAS_highPt1stEle", getHistoNBins("Dist1stEle_PAS"), getHistoMin("Dist1stEle_PAS"), getHistoMax("Dist1stEle_PAS") );
      CreateUserTH1D("h1_DCotTheta1stEle_PAS_highPt1stEle", getHistoNBins("DCotTheta1stEle_PAS"), getHistoMin("DCotTheta1stEle_PAS"), getHistoMax("DCotTheta1stEle_PAS") );
      CreateUserTH1D("h1_Conversion1stEle_highPt1stEle", 2, -0.5, 1.5 );

      CreateUserTH1D("h1_Njet_highPt1stJet", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_NjetTCHELBTag_highPt1stJet", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_minDRej_highPt1stJet", getHistoNBins("minDRej"), getHistoMin("minDRej"), getHistoMax("minDRej") );
      CreateUserTH1D("h1_mDeltaPhiMETEle_highPt1stJet", getHistoNBins("mDeltaPhiMETEle"), getHistoMin("mDeltaPhiMETEle"), getHistoMax("mDeltaPhiMETEle") );
      CreateUserTH1D("h1_mDeltaPhiMET1stJet_highPt1stJet", getHistoNBins("mDeltaPhiMET1stJet"), getHistoMin("mDeltaPhiMET1stJet"), getHistoMax("mDeltaPhiMET1stJet") );
      CreateUserTH1D("h1_mDeltaPhiMET2ndJet_highPt1stJet", getHistoNBins("mDeltaPhiMET2ndJet"), getHistoMin("mDeltaPhiMET2ndJet"), getHistoMax("mDeltaPhiMET2ndJet") );
      CreateUserTH1D("h1_DeltaRjets_PAS_highPt1stJet", getHistoNBins("DeltaRjets_PAS"), getHistoMin("DeltaRjets_PAS"), getHistoMax("DeltaRjets_PAS") );
      CreateUserTH1D("h1_Vtxd01stEle_PAS_highPt1stJet", getHistoNBins("Vtxd01stEle_PAS"), getHistoMin("Vtxd01stEle_PAS"), getHistoMax("Vtxd01stEle_PAS") );
      CreateUserTH1D("h1_MissingHits1stEle_PAS_highPt1stJet", getHistoNBins("MissingHits1stEle_PAS"), getHistoMin("MissingHits1stEle_PAS"), getHistoMax("MissingHits1stEle_PAS") );
      CreateUserTH1D("h1_Dist1stEle_PAS_highPt1stJet", getHistoNBins("Dist1stEle_PAS"), getHistoMin("Dist1stEle_PAS"), getHistoMax("Dist1stEle_PAS") );
      CreateUserTH1D("h1_DCotTheta1stEle_PAS_highPt1stJet", getHistoNBins("DCotTheta1stEle_PAS"), getHistoMin("DCotTheta1stEle_PAS"), getHistoMax("DCotTheta1stEle_PAS") );
      CreateUserTH1D("h1_Conversion1stEle_highPt1stJet", 2, -0.5, 1.5 );

      CreateUserTH1D("h1_Njet_highPt2ndJet", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_NjetTCHELBTag_highPt2ndJet", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_minDRej_highPt2ndJet", getHistoNBins("minDRej"), getHistoMin("minDRej"), getHistoMax("minDRej") );
      CreateUserTH1D("h1_mDeltaPhiMETEle_highPt2ndJet", getHistoNBins("mDeltaPhiMETEle"), getHistoMin("mDeltaPhiMETEle"), getHistoMax("mDeltaPhiMETEle") );
      CreateUserTH1D("h1_mDeltaPhiMET1stJet_highPt2ndJet", getHistoNBins("mDeltaPhiMET1stJet"), getHistoMin("mDeltaPhiMET1stJet"), getHistoMax("mDeltaPhiMET1stJet") );
      CreateUserTH1D("h1_mDeltaPhiMET2ndJet_highPt2ndJet", getHistoNBins("mDeltaPhiMET2ndJet"), getHistoMin("mDeltaPhiMET2ndJet"), getHistoMax("mDeltaPhiMET2ndJet") );
      CreateUserTH1D("h1_DeltaRjets_PAS_highPt2ndJet", getHistoNBins("DeltaRjets_PAS"), getHistoMin("DeltaRjets_PAS"), getHistoMax("DeltaRjets_PAS") );
      CreateUserTH1D("h1_Vtxd01stEle_PAS_highPt2ndJet", getHistoNBins("Vtxd01stEle_PAS"), getHistoMin("Vtxd01stEle_PAS"), getHistoMax("Vtxd01stEle_PAS") );
      CreateUserTH1D("h1_MissingHits1stEle_PAS_highPt2ndJet", getHistoNBins("MissingHits1stEle_PAS"), getHistoMin("MissingHits1stEle_PAS"), getHistoMax("MissingHits1stEle_PAS") );
      CreateUserTH1D("h1_Dist1stEle_PAS_highPt2ndJet", getHistoNBins("Dist1stEle_PAS"), getHistoMin("Dist1stEle_PAS"), getHistoMax("Dist1stEle_PAS") );
      CreateUserTH1D("h1_DCotTheta1stEle_PAS_highPt2ndJet", getHistoNBins("DCotTheta1stEle_PAS"), getHistoMin("DCotTheta1stEle_PAS"), getHistoMax("DCotTheta1stEle_PAS") );
      CreateUserTH1D("h1_Conversion1stEle_highPt2ndJet", 2, -0.5, 1.5 );

      CreateUserTH1D("h1_Njet_highMT", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_NjetTCHELBTag_highMT", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_minDRej_highMT", getHistoNBins("minDRej"), getHistoMin("minDRej"), getHistoMax("minDRej") );
      CreateUserTH1D("h1_mDeltaPhiMETEle_highMT", getHistoNBins("mDeltaPhiMETEle"), getHistoMin("mDeltaPhiMETEle"), getHistoMax("mDeltaPhiMETEle") );
      CreateUserTH1D("h1_mDeltaPhiMET1stJet_highMT", getHistoNBins("mDeltaPhiMET1stJet"), getHistoMin("mDeltaPhiMET1stJet"), getHistoMax("mDeltaPhiMET1stJet") );
      CreateUserTH1D("h1_mDeltaPhiMET2ndJet_highMT", getHistoNBins("mDeltaPhiMET2ndJet"), getHistoMin("mDeltaPhiMET2ndJet"), getHistoMax("mDeltaPhiMET2ndJet") );
      CreateUserTH1D("h1_DeltaRjets_PAS_highMT", getHistoNBins("DeltaRjets_PAS"), getHistoMin("DeltaRjets_PAS"), getHistoMax("DeltaRjets_PAS") );
      CreateUserTH1D("h1_Vtxd01stEle_PAS_highMT", getHistoNBins("Vtxd01stEle_PAS"), getHistoMin("Vtxd01stEle_PAS"), getHistoMax("Vtxd01stEle_PAS") );
      CreateUserTH1D("h1_MissingHits1stEle_PAS_highMT", getHistoNBins("MissingHits1stEle_PAS"), getHistoMin("MissingHits1stEle_PAS"), getHistoMax("MissingHits1stEle_PAS") );
      CreateUserTH1D("h1_Dist1stEle_PAS_highMT", getHistoNBins("Dist1stEle_PAS"), getHistoMin("Dist1stEle_PAS"), getHistoMax("Dist1stEle_PAS") );
      CreateUserTH1D("h1_DCotTheta1stEle_PAS_highMT", getHistoNBins("DCotTheta1stEle_PAS"), getHistoMin("DCotTheta1stEle_PAS"), getHistoMax("DCotTheta1stEle_PAS") );
      CreateUserTH1D("h1_Conversion1stEle_highMT", 2, -0.5, 1.5 );

      CreateUserTH1D("h1_Njet_highMej", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_NjetTCHELBTag_highMej", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_minDRej_highMej", getHistoNBins("minDRej"), getHistoMin("minDRej"), getHistoMax("minDRej") );
      CreateUserTH1D("h1_mDeltaPhiMETEle_highMej", getHistoNBins("mDeltaPhiMETEle"), getHistoMin("mDeltaPhiMETEle"), getHistoMax("mDeltaPhiMETEle") );
      CreateUserTH1D("h1_mDeltaPhiMET1stJet_highMej", getHistoNBins("mDeltaPhiMET1stJet"), getHistoMin("mDeltaPhiMET1stJet"), getHistoMax("mDeltaPhiMET1stJet") );
      CreateUserTH1D("h1_mDeltaPhiMET2ndJet_highMej", getHistoNBins("mDeltaPhiMET2ndJet"), getHistoMin("mDeltaPhiMET2ndJet"), getHistoMax("mDeltaPhiMET2ndJet") );
      CreateUserTH1D("h1_DeltaRjets_PAS_highMej", getHistoNBins("DeltaRjets_PAS"), getHistoMin("DeltaRjets_PAS"), getHistoMax("DeltaRjets_PAS") );
      CreateUserTH1D("h1_Vtxd01stEle_PAS_highMej", getHistoNBins("Vtxd01stEle_PAS"), getHistoMin("Vtxd01stEle_PAS"), getHistoMax("Vtxd01stEle_PAS") );
      CreateUserTH1D("h1_MissingHits1stEle_PAS_highMej", getHistoNBins("MissingHits1stEle_PAS"), getHistoMin("MissingHits1stEle_PAS"), getHistoMax("MissingHits1stEle_PAS") );
      CreateUserTH1D("h1_Dist1stEle_PAS_highMej", getHistoNBins("Dist1stEle_PAS"), getHistoMin("Dist1stEle_PAS"), getHistoMax("Dist1stEle_PAS") );
      CreateUserTH1D("h1_DCotTheta1stEle_PAS_highMej", getHistoNBins("DCotTheta1stEle_PAS"), getHistoMin("DCotTheta1stEle_PAS"), getHistoMax("DCotTheta1stEle_PAS") );
      CreateUserTH1D("h1_Conversion1stEle_highMej", 2, -0.5, 1.5 );

      CreateUserTH1D("h1_Pt1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Pt1stJet_PAS"), getHistoMin("Pt1stJet_PAS"), getHistoMax("Pt1stJet_PAS") );
      CreateUserTH1D("h1_Eta1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Eta1stJet_PAS"), getHistoMin("Eta1stJet_PAS"), getHistoMax("Eta1stJet_PAS") );
      CreateUserTH1D("h1_Phi1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Phi1stJet_PAS"), getHistoMin("Phi1stJet_PAS"), getHistoMax("Phi1stJet_PAS") );
      CreateUserTH1D("h1_CHF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NHF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_CEF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NEF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NCH1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NN1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NC1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_Mass1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 100 );
      CreateUserTH1D("h1_minDR_1stPFjet_CaloJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 2 );
      CreateUserTH1D("h1_PT1stJet_over_PTCaloJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 200, 0, 2 );
      CreateUserTH1D("h1_NC1stJet_over_n90HitsCaloJet_1stPFjet_CaloJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 5 );
      CreateUserTH1D("h1_Pt2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Pt2ndJet_PAS"), getHistoMin("Pt2ndJet_PAS"), getHistoMax("Pt2ndJet_PAS") );
      CreateUserTH1D("h1_Eta2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Eta2ndJet_PAS"), getHistoMin("Eta2ndJet_PAS"), getHistoMax("Eta2ndJet_PAS") );
      CreateUserTH1D("h1_Phi2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Phi2ndJet_PAS"), getHistoMin("Phi2ndJet_PAS"), getHistoMax("Phi2ndJet_PAS") );
      CreateUserTH1D("h1_CHF2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NHF2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_CEF2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NEF2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NCH2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NN2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NC2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_Mass2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 100 );
      CreateUserTH1D("h1_E1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
      CreateUserTH1D("h1_Pt1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
      CreateUserTH1D("h1_Eta1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Eta1stEle_PAS"), getHistoMin("Eta1stEle_PAS"), getHistoMax("Eta1stEle_PAS") );
      CreateUserTH1D("h1_Phi1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Phi1stEle_PAS"), getHistoMin("Phi1stEle_PAS"), getHistoMax("Phi1stEle_PAS") );
      CreateUserTH1D("h1_Charge1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Charge1stEle_PAS"), getHistoMin("Charge1stEle_PAS"), getHistoMax("Charge1stEle_PAS") );
      CreateUserTH1D("h1_mDeltaPhiMETEle_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("mDeltaPhiMETEle"), getHistoMin("mDeltaPhiMETEle"), getHistoMax("mDeltaPhiMETEle") );
      CreateUserTH1D("h1_mDeltaPhiEle1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("mDeltaPhiEle1stJet_PAS"), getHistoMin("mDeltaPhiEle1stJet_PAS"), getHistoMax("mDeltaPhiEle1stJet_PAS") );
      CreateUserTH1D("h1_mDeltaPhiEle2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("mDeltaPhiEle2ndJet_PAS"), getHistoMin("mDeltaPhiEle2ndJet_PAS"), getHistoMax("mDeltaPhiEle2ndJet_PAS") );
      CreateUserTH1D("h1_mDeltaPhiMET2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("mDeltaPhiMET2ndJet"), getHistoMin("mDeltaPhiMET2ndJet"), getHistoMax("mDeltaPhiMET2ndJet") );
      CreateUserTH1D("h1_Ptenu_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Ptenu_PAS"), getHistoMin("Ptenu_PAS"), getHistoMax("Ptenu_PAS") );
      CreateUserTH1D("h1_1stJet_PTOverPTPlusMET_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1.0 );
      CreateUserTH1D("h1_Conversion1stEle_highMej_mDeltaPhiMET1stJet_gt_2.5", 2, -0.5, 1.5 );
      CreateUserTH1D("h1_MET_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_CaloMET_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_TCMET_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_MET_over_CaloMET_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 5 );
      CreateUserTH1D("h1_Njet_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_NjetTCHELBTag_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_minDRej_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("minDRej"), getHistoMin("minDRej"), getHistoMax("minDRej") );
      CreateUserTH1D("h1_maxDRej_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("maxDRej"), getHistoMin("maxDRej"), getHistoMax("maxDRej") );
      CreateUserTH1D("h1_DRjets_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("DeltaRjets_PAS"), getHistoMin("DeltaRjets_PAS"), getHistoMax("DeltaRjets_PAS") );
      CreateUserTH1D("h1_MTenu_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("MTenu_PAS"), getHistoMin("MTenu_PAS"), getHistoMax("MTenu_PAS") );
      CreateUserTH2D("h2_EtaPhi1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", 100, -5, 5, 60, -3.1416, 3.1416 );
      CreateUserTH2D("h2_MinMaxMej_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Mej_1stPair_PAS"), getHistoMin("Mej_1stPair_PAS"), getHistoMax("Mej_1stPair_PAS"), getHistoNBins("Mej_2ndPair_PAS"), getHistoMin("Mej_2ndPair_PAS"), getHistoMax("Mej_2ndPair_PAS") );
      CreateUserTH2D("h2_MTnuj2_vs_Mej1_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Mej_1stPair_PAS"), getHistoMin("Mej_1stPair_PAS"), getHistoMax("Mej_1stPair_PAS"), getHistoNBins("MTnuj_1stPair_PAS"), getHistoMin("MTnuj_1stPair_PAS"), getHistoMax("MTnuj_1stPair_PAS") );
      CreateUserTH2D("h2_MTnuj1_vs_Mej2_highMej_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Mej_1stPair_PAS"), getHistoMin("Mej_1stPair_PAS"), getHistoMax("Mej_1stPair_PAS"), getHistoNBins("MTnuj_1stPair_PAS"), getHistoMin("MTnuj_1stPair_PAS"), getHistoMax("MTnuj_1stPair_PAS") );

      CreateUserTH1D("h1_Pt1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Pt1stJet_PAS"), getHistoMin("Pt1stJet_PAS"), getHistoMax("Pt1stJet_PAS") );
      CreateUserTH1D("h1_Eta1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Eta1stJet_PAS"), getHistoMin("Eta1stJet_PAS"), getHistoMax("Eta1stJet_PAS") );
      CreateUserTH1D("h1_Phi1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Phi1stJet_PAS"), getHistoMin("Phi1stJet_PAS"), getHistoMax("Phi1stJet_PAS") );
      CreateUserTH1D("h1_CHF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NHF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_CEF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NEF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NCH1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NN1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NC1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_Mass1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 100 );
      CreateUserTH1D("h1_Pt2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Pt2ndJet_PAS"), getHistoMin("Pt2ndJet_PAS"), getHistoMax("Pt2ndJet_PAS") );
      CreateUserTH1D("h1_Eta2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Eta2ndJet_PAS"), getHistoMin("Eta2ndJet_PAS"), getHistoMax("Eta2ndJet_PAS") );
      CreateUserTH1D("h1_Phi2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Phi2ndJet_PAS"), getHistoMin("Phi2ndJet_PAS"), getHistoMax("Phi2ndJet_PAS") );
      CreateUserTH1D("h1_CHF2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NHF2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_CEF2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NEF2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NCH2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NN2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NC2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_Mass2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 100 );
      CreateUserTH1D("h1_E1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
      CreateUserTH1D("h1_Pt1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
      CreateUserTH1D("h1_Eta1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Eta1stEle_PAS"), getHistoMin("Eta1stEle_PAS"), getHistoMax("Eta1stEle_PAS") );
      CreateUserTH1D("h1_Phi1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Phi1stEle_PAS"), getHistoMin("Phi1stEle_PAS"), getHistoMax("Phi1stEle_PAS") );
      CreateUserTH1D("h1_Charge1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Charge1stEle_PAS"), getHistoMin("Charge1stEle_PAS"), getHistoMax("Charge1stEle_PAS") );
      CreateUserTH1D("h1_mDeltaPhiMETEle_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("mDeltaPhiMETEle"), getHistoMin("mDeltaPhiMETEle"), getHistoMax("mDeltaPhiMETEle") );
      CreateUserTH1D("h1_mDeltaPhiEle1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("mDeltaPhiEle1stJet_PAS"), getHistoMin("mDeltaPhiEle1stJet_PAS"), getHistoMax("mDeltaPhiEle1stJet_PAS") );
      CreateUserTH1D("h1_mDeltaPhiEle2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("mDeltaPhiEle2ndJet_PAS"), getHistoMin("mDeltaPhiEle2ndJet_PAS"), getHistoMax("mDeltaPhiEle2ndJet_PAS") );
      CreateUserTH1D("h1_mDeltaPhiMET2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("mDeltaPhiMET2ndJet"), getHistoMin("mDeltaPhiMET2ndJet"), getHistoMax("mDeltaPhiMET2ndJet") );
      CreateUserTH1D("h1_Ptenu_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Ptenu_PAS"), getHistoMin("Ptenu_PAS"), getHistoMax("Ptenu_PAS") );
      CreateUserTH1D("h1_1stJet_PTOverPTPlusMET_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1.0 );
      CreateUserTH1D("h1_Conversion1stEle_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 2, -0.5, 1.5 );
      CreateUserTH1D("h1_MET_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_CaloMET_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_TCMET_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_MET_over_CaloMET_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 5 );
      CreateUserTH1D("h1_Njet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_NjetTCHELBTag_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_minDRej_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("minDRej"), getHistoMin("minDRej"), getHistoMax("minDRej") );
      CreateUserTH1D("h1_maxDRej_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("maxDRej"), getHistoMin("maxDRej"), getHistoMax("maxDRej") );
      CreateUserTH1D("h1_DRjets_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("DeltaRjets_PAS"), getHistoMin("DeltaRjets_PAS"), getHistoMax("DeltaRjets_PAS") );
      CreateUserTH1D("h1_MTenu_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("MTenu_PAS"), getHistoMin("MTenu_PAS"), getHistoMax("MTenu_PAS") );
      CreateUserTH2D("h2_EtaPhi1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", 100, -5, 5, 60, -3.1416, 3.1416 );
      CreateUserTH2D("h2_MinMaxMej_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Mej_1stPair_PAS"), getHistoMin("Mej_1stPair_PAS"), getHistoMax("Mej_1stPair_PAS"), getHistoNBins("Mej_2ndPair_PAS"), getHistoMin("Mej_2ndPair_PAS"), getHistoMax("Mej_2ndPair_PAS") );

      CreateUserTH1D("h1_Pt1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Pt1stJet_PAS"), getHistoMin("Pt1stJet_PAS"), getHistoMax("Pt1stJet_PAS") );
      CreateUserTH1D("h1_Eta1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Eta1stJet_PAS"), getHistoMin("Eta1stJet_PAS"), getHistoMax("Eta1stJet_PAS") );
      CreateUserTH1D("h1_Phi1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Phi1stJet_PAS"), getHistoMin("Phi1stJet_PAS"), getHistoMax("Phi1stJet_PAS") );
      CreateUserTH1D("h1_CHF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NHF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_CEF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NEF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NCH1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NN1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NC1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_Mass1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 100 );
      CreateUserTH1D("h1_Pt2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Pt2ndJet_PAS"), getHistoMin("Pt2ndJet_PAS"), getHistoMax("Pt2ndJet_PAS") );
      CreateUserTH1D("h1_Eta2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Eta2ndJet_PAS"), getHistoMin("Eta2ndJet_PAS"), getHistoMax("Eta2ndJet_PAS") );
      CreateUserTH1D("h1_Phi2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Phi2ndJet_PAS"), getHistoMin("Phi2ndJet_PAS"), getHistoMax("Phi2ndJet_PAS") );
      CreateUserTH1D("h1_CHF2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NHF2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_CEF2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NEF2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NCH2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NN2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NC2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_Mass2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 100 );
      CreateUserTH1D("h1_E1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
      CreateUserTH1D("h1_Pt1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
      CreateUserTH1D("h1_Eta1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Eta1stEle_PAS"), getHistoMin("Eta1stEle_PAS"), getHistoMax("Eta1stEle_PAS") );
      CreateUserTH1D("h1_Phi1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Phi1stEle_PAS"), getHistoMin("Phi1stEle_PAS"), getHistoMax("Phi1stEle_PAS") );
      CreateUserTH1D("h1_Charge1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Charge1stEle_PAS"), getHistoMin("Charge1stEle_PAS"), getHistoMax("Charge1stEle_PAS") );
      CreateUserTH1D("h1_mDeltaPhiMETEle_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("mDeltaPhiMETEle"), getHistoMin("mDeltaPhiMETEle"), getHistoMax("mDeltaPhiMETEle") );
      CreateUserTH1D("h1_mDeltaPhiEle1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("mDeltaPhiEle1stJet_PAS"), getHistoMin("mDeltaPhiEle1stJet_PAS"), getHistoMax("mDeltaPhiEle1stJet_PAS") );
      CreateUserTH1D("h1_mDeltaPhiEle2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("mDeltaPhiEle2ndJet_PAS"), getHistoMin("mDeltaPhiEle2ndJet_PAS"), getHistoMax("mDeltaPhiEle2ndJet_PAS") );
      CreateUserTH1D("h1_mDeltaPhiMET2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("mDeltaPhiMET2ndJet"), getHistoMin("mDeltaPhiMET2ndJet"), getHistoMax("mDeltaPhiMET2ndJet") );
      CreateUserTH1D("h1_Ptenu_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Ptenu_PAS"), getHistoMin("Ptenu_PAS"), getHistoMax("Ptenu_PAS") );
      CreateUserTH1D("h1_1stJet_PTOverPTPlusMET_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 1.0 );
      CreateUserTH1D("h1_Conversion1stEle_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 2, -0.5, 1.5 );
      CreateUserTH1D("h1_MET_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_CaloMET_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_TCMET_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_MET_over_CaloMET_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, 0, 5 );
      CreateUserTH1D("h1_Njet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_NjetTCHELBTag_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_minDRej_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("minDRej"), getHistoMin("minDRej"), getHistoMax("minDRej") );
      CreateUserTH1D("h1_maxDRej_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("maxDRej"), getHistoMin("maxDRej"), getHistoMax("maxDRej") );
      CreateUserTH1D("h1_DRjets_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("DeltaRjets_PAS"), getHistoMin("DeltaRjets_PAS"), getHistoMax("DeltaRjets_PAS") );
      CreateUserTH1D("h1_MTenu_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("MTenu_PAS"), getHistoMin("MTenu_PAS"), getHistoMax("MTenu_PAS") );
      CreateUserTH2D("h2_EtaPhi1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", 100, -5, 5, 60, -3.1416, 3.1416 );
      CreateUserTH2D("h2_MinMaxMej_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getHistoNBins("Mej_1stPair_PAS"), getHistoMin("Mej_1stPair_PAS"), getHistoMax("Mej_1stPair_PAS"), getHistoNBins("Mej_2ndPair_PAS"), getHistoMin("Mej_2ndPair_PAS"), getHistoMax("Mej_2ndPair_PAS") );

      CreateUserTH1D("h1_Pt1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("Pt1stJet_PAS"), getHistoMin("Pt1stJet_PAS"), getHistoMax("Pt1stJet_PAS") );
      CreateUserTH1D("h1_Eta1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("Eta1stJet_PAS"), getHistoMin("Eta1stJet_PAS"), getHistoMax("Eta1stJet_PAS") );
      CreateUserTH1D("h1_Phi1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("Phi1stJet_PAS"), getHistoMin("Phi1stJet_PAS"), getHistoMax("Phi1stJet_PAS") );
      CreateUserTH1D("h1_CHF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NHF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_CEF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NEF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NCH1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NN1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NC1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_Mass1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, 0, 100 );
      CreateUserTH1D("h1_Pt2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("Pt2ndJet_PAS"), getHistoMin("Pt2ndJet_PAS"), getHistoMax("Pt2ndJet_PAS") );
      CreateUserTH1D("h1_Eta2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("Eta2ndJet_PAS"), getHistoMin("Eta2ndJet_PAS"), getHistoMax("Eta2ndJet_PAS") );
      CreateUserTH1D("h1_Phi2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("Phi2ndJet_PAS"), getHistoMin("Phi2ndJet_PAS"), getHistoMax("Phi2ndJet_PAS") );
      CreateUserTH1D("h1_CHF2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NHF2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_CEF2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NEF2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, 0, 1 );
      CreateUserTH1D("h1_NCH2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NN2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NC2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_Mass2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", 100, 0, 100 );
      CreateUserTH1D("h1_E1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
      CreateUserTH1D("h1_Pt1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
      CreateUserTH1D("h1_Eta1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("Eta1stEle_PAS"), getHistoMin("Eta1stEle_PAS"), getHistoMax("Eta1stEle_PAS") );
      CreateUserTH1D("h1_Phi1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("Phi1stEle_PAS"), getHistoMin("Phi1stEle_PAS"), getHistoMax("Phi1stEle_PAS") );
      CreateUserTH1D("h1_Charge1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("Charge1stEle_PAS"), getHistoMin("Charge1stEle_PAS"), getHistoMax("Charge1stEle_PAS") );
      CreateUserTH1D("h1_mDeltaPhiMETEle_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("mDeltaPhiMETEle"), getHistoMin("mDeltaPhiMETEle"), getHistoMax("mDeltaPhiMETEle") );
      CreateUserTH1D("h1_mDeltaPhiEle1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("mDeltaPhiEle1stJet_PAS"), getHistoMin("mDeltaPhiEle1stJet_PAS"), getHistoMax("mDeltaPhiEle1stJet_PAS") );
      CreateUserTH1D("h1_mDeltaPhiEle2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("mDeltaPhiEle2ndJet_PAS"), getHistoMin("mDeltaPhiEle2ndJet_PAS"), getHistoMax("mDeltaPhiEle2ndJet_PAS") );
      CreateUserTH1D("h1_mDeltaPhiMET2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("mDeltaPhiMET2ndJet"), getHistoMin("mDeltaPhiMET2ndJet"), getHistoMax("mDeltaPhiMET2ndJet") );
      CreateUserTH1D("h1_Ptenu_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("Ptenu_PAS"), getHistoMin("Ptenu_PAS"), getHistoMax("Ptenu_PAS") );
      CreateUserTH1D("h1_1stJet_PTOverPTPlusMET_highMej_mDeltaPhiMET1stJet_le_2.5", 100, 0, 1.0 );
      CreateUserTH1D("h1_Conversion1stEle_highMej_mDeltaPhiMET1stJet_le_2.5", 2, -0.5, 1.5 );
      CreateUserTH1D("h1_MET_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_CaloMET_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_TCMET_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_MET_over_CaloMET_highMej_mDeltaPhiMET1stJet_le_2.5", 100, 0, 5 );
      CreateUserTH1D("h1_Njet_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_NjetTCHELBTag_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_minDRej_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("minDRej"), getHistoMin("minDRej"), getHistoMax("minDRej") );
      CreateUserTH1D("h1_maxDRej_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("maxDRej"), getHistoMin("maxDRej"), getHistoMax("maxDRej") );
      CreateUserTH1D("h1_DRjets_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("DeltaRjets_PAS"), getHistoMin("DeltaRjets_PAS"), getHistoMax("DeltaRjets_PAS") );
      CreateUserTH1D("h1_MTenu_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("MTenu_PAS"), getHistoMin("MTenu_PAS"), getHistoMax("MTenu_PAS") );
      CreateUserTH2D("h2_EtaPhi1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", 100, -5, 5, 60, -3.1416, 3.1416 );
      CreateUserTH2D("h2_MinMaxMej_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getHistoNBins("Mej_1stPair_PAS"), getHistoMin("Mej_1stPair_PAS"), getHistoMax("Mej_1stPair_PAS"), getHistoNBins("Mej_2ndPair_PAS"), getHistoMin("Mej_2ndPair_PAS"), getHistoMax("Mej_2ndPair_PAS") );

      CreateUserTH1D("h1_Njet_fullSel", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_NjetTCHELBTag_fullSel", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_minDRej_fullSel", getHistoNBins("minDRej"), getHistoMin("minDRej"), getHistoMax("minDRej") );
      CreateUserTH1D("h1_mDeltaPhiMETEle_fullSel", getHistoNBins("mDeltaPhiMETEle"), getHistoMin("mDeltaPhiMETEle"), getHistoMax("mDeltaPhiMETEle") );
      CreateUserTH1D("h1_mDeltaPhiMET1stJet_fullSel", getHistoNBins("mDeltaPhiMET1stJet"), getHistoMin("mDeltaPhiMET1stJet"), getHistoMax("mDeltaPhiMET1stJet") );
      CreateUserTH1D("h1_mDeltaPhiMET2ndJet_fullSel", getHistoNBins("mDeltaPhiMET2ndJet"), getHistoMin("mDeltaPhiMET2ndJet"), getHistoMax("mDeltaPhiMET2ndJet") );
      CreateUserTH1D("h1_DeltaRjets_PAS_fullSel", getHistoNBins("DeltaRjets_PAS"), getHistoMin("DeltaRjets_PAS"), getHistoMax("DeltaRjets_PAS") );
      CreateUserTH1D("h1_Vtxd01stEle_PAS_fullSel", getHistoNBins("Vtxd01stEle_PAS"), getHistoMin("Vtxd01stEle_PAS"), getHistoMax("Vtxd01stEle_PAS") );
      CreateUserTH1D("h1_MissingHits1stEle_PAS_fullSel", getHistoNBins("MissingHits1stEle_PAS"), getHistoMin("MissingHits1stEle_PAS"), getHistoMax("MissingHits1stEle_PAS") );
      CreateUserTH1D("h1_Dist1stEle_PAS_fullSel", getHistoNBins("Dist1stEle_PAS"), getHistoMin("Dist1stEle_PAS"), getHistoMax("Dist1stEle_PAS") );
      CreateUserTH1D("h1_DCotTheta1stEle_PAS_fullSel", getHistoNBins("DCotTheta1stEle_PAS"), getHistoMin("DCotTheta1stEle_PAS"), getHistoMax("DCotTheta1stEle_PAS") );
      CreateUserTH1D("h1_Conversion1stEle_fullSel", 2, -0.5, 1.5 );

      CreateUserTH1D("h1_Pt1stJet_PAS_Eta1stJetBump", getHistoNBins("Pt1stJet_PAS"), getHistoMin("Pt1stJet_PAS"), getHistoMax("Pt1stJet_PAS") );
      CreateUserTH1D("h1_Phi1stJet_PAS_Eta1stJetBump", getHistoNBins("Phi1stJet_PAS"), getHistoMin("Phi1stJet_PAS"), getHistoMax("Phi1stJet_PAS") );
      CreateUserTH1D("h1_CHF1stJet_Eta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_NHF1stJet_Eta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_CEF1stJet_Eta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_NEF1stJet_Eta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_NCH1stJet_Eta1stJetBump", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NN1stJet_Eta1stJetBump", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NC1stJet_Eta1stJetBump", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_Mass1stJet_Eta1stJetBump", 100, 0, 100 );
      CreateUserTH1D("h1_Pt2ndJet_PAS_Eta1stJetBump", getHistoNBins("Pt2ndJet_PAS"), getHistoMin("Pt2ndJet_PAS"), getHistoMax("Pt2ndJet_PAS") );
      CreateUserTH1D("h1_Eta2ndJet_PAS_Eta1stJetBump", getHistoNBins("Eta2ndJet_PAS"), getHistoMin("Eta2ndJet_PAS"), getHistoMax("Eta2ndJet_PAS") );
      CreateUserTH1D("h1_Phi2ndJet_PAS_Eta1stJetBump", getHistoNBins("Phi2ndJet_PAS"), getHistoMin("Phi2ndJet_PAS"), getHistoMax("Phi2ndJet_PAS") );
      CreateUserTH1D("h1_CHF2ndJet_Eta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_NHF2ndJet_Eta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_CEF2ndJet_Eta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_NEF2ndJet_Eta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_NCH2ndJet_Eta1stJetBump", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NN2ndJet_Eta1stJetBump", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NC2ndJet_Eta1stJetBump", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_Mass2ndJet_Eta1stJetBump", 100, 0, 100 );
      CreateUserTH1D("h1_E1stEle_PAS_Eta1stJetBump", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
      CreateUserTH1D("h1_Pt1stEle_PAS_Eta1stJetBump", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
      CreateUserTH1D("h1_Eta1stEle_PAS_Eta1stJetBump", getHistoNBins("Eta1stEle_PAS"), getHistoMin("Eta1stEle_PAS"), getHistoMax("Eta1stEle_PAS") );
      CreateUserTH1D("h1_Phi1stEle_PAS_Eta1stJetBump", getHistoNBins("Phi1stEle_PAS"), getHistoMin("Phi1stEle_PAS"), getHistoMax("Phi1stEle_PAS") );
      CreateUserTH1D("h1_Charge1stEle_PAS_Eta1stJetBump", getHistoNBins("Charge1stEle_PAS"), getHistoMin("Charge1stEle_PAS"), getHistoMax("Charge1stEle_PAS") );
      CreateUserTH1D("h1_mDeltaPhiMETEle_Eta1stJetBump", getHistoNBins("mDeltaPhiMETEle"), getHistoMin("mDeltaPhiMETEle"), getHistoMax("mDeltaPhiMETEle") );
      CreateUserTH1D("h1_mDeltaPhiEle1stJet_PAS_Eta1stJetBump", getHistoNBins("mDeltaPhiEle1stJet_PAS"), getHistoMin("mDeltaPhiEle1stJet_PAS"), getHistoMax("mDeltaPhiEle1stJet_PAS") );
      CreateUserTH1D("h1_mDeltaPhiEle2ndJet_PAS_Eta1stJetBump", getHistoNBins("mDeltaPhiEle2ndJet_PAS"), getHistoMin("mDeltaPhiEle2ndJet_PAS"), getHistoMax("mDeltaPhiEle2ndJet_PAS") );
      CreateUserTH1D("h1_mDeltaPhiMET2ndJet_Eta1stJetBump", getHistoNBins("mDeltaPhiMET2ndJet"), getHistoMin("mDeltaPhiMET2ndJet"), getHistoMax("mDeltaPhiMET2ndJet") );
      CreateUserTH1D("h1_Ptenu_PAS_Eta1stJetBump", getHistoNBins("Ptenu_PAS"), getHistoMin("Ptenu_PAS"), getHistoMax("Ptenu_PAS") );
      CreateUserTH1D("h1_1stJet_PTOverPTPlusMET_Eta1stJetBump", 100, 0, 1.0 );
      CreateUserTH1D("h1_Conversion1stEle_Eta1stJetBump", 2, -0.5, 1.5 );
      CreateUserTH1D("h1_MET_PAS_Eta1stJetBump", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_CaloMET_PAS_Eta1stJetBump", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_TCMET_PAS_Eta1stJetBump", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_MET_over_CaloMET_Eta1stJetBump", 100, 0, 5 );
      CreateUserTH1D("h1_Njet_Eta1stJetBump", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_NjetTCHELBTag_Eta1stJetBump", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_minDRej_Eta1stJetBump", getHistoNBins("minDRej"), getHistoMin("minDRej"), getHistoMax("minDRej") );
      CreateUserTH1D("h1_maxDRej_Eta1stJetBump", getHistoNBins("maxDRej"), getHistoMin("maxDRej"), getHistoMax("maxDRej") );
      CreateUserTH1D("h1_DRjets_Eta1stJetBump", getHistoNBins("DeltaRjets_PAS"), getHistoMin("DeltaRjets_PAS"), getHistoMax("DeltaRjets_PAS") );
      CreateUserTH1D("h1_MTenu_PAS_Eta1stJetBump", getHistoNBins("MTenu_PAS"), getHistoMin("MTenu_PAS"), getHistoMax("MTenu_PAS") );

      CreateUserTH1D("h1_Pt1stJet_PAS_OutsideEta1stJetBump", getHistoNBins("Pt1stJet_PAS"), getHistoMin("Pt1stJet_PAS"), getHistoMax("Pt1stJet_PAS") );
      CreateUserTH1D("h1_Phi1stJet_PAS_OutsideEta1stJetBump", getHistoNBins("Phi1stJet_PAS"), getHistoMin("Phi1stJet_PAS"), getHistoMax("Phi1stJet_PAS") );
      CreateUserTH1D("h1_CHF1stJet_OutsideEta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_NHF1stJet_OutsideEta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_CEF1stJet_OutsideEta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_NEF1stJet_OutsideEta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_NCH1stJet_OutsideEta1stJetBump", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NN1stJet_OutsideEta1stJetBump", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NC1stJet_OutsideEta1stJetBump", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_Mass1stJet_OutsideEta1stJetBump", 100, 0, 100 );
      CreateUserTH1D("h1_Pt2ndJet_PAS_OutsideEta1stJetBump", getHistoNBins("Pt2ndJet_PAS"), getHistoMin("Pt2ndJet_PAS"), getHistoMax("Pt2ndJet_PAS") );
      CreateUserTH1D("h1_Eta2ndJet_PAS_OutsideEta1stJetBump", getHistoNBins("Eta2ndJet_PAS"), getHistoMin("Eta2ndJet_PAS"), getHistoMax("Eta2ndJet_PAS") );
      CreateUserTH1D("h1_Phi2ndJet_PAS_OutsideEta1stJetBump", getHistoNBins("Phi2ndJet_PAS"), getHistoMin("Phi2ndJet_PAS"), getHistoMax("Phi2ndJet_PAS") );
      CreateUserTH1D("h1_CHF2ndJet_OutsideEta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_NHF2ndJet_OutsideEta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_CEF2ndJet_OutsideEta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_NEF2ndJet_OutsideEta1stJetBump", 100, 0, 1 );
      CreateUserTH1D("h1_NCH2ndJet_OutsideEta1stJetBump", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NN2ndJet_OutsideEta1stJetBump", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_NC2ndJet_OutsideEta1stJetBump", 100, -0.5, 99.5 );
      CreateUserTH1D("h1_Mass2ndJet_OutsideEta1stJetBump", 100, 0, 100 );
      CreateUserTH1D("h1_E1stEle_PAS_OutsideEta1stJetBump", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
      CreateUserTH1D("h1_Pt1stEle_PAS_OutsideEta1stJetBump", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
      CreateUserTH1D("h1_Eta1stEle_PAS_OutsideEta1stJetBump", getHistoNBins("Eta1stEle_PAS"), getHistoMin("Eta1stEle_PAS"), getHistoMax("Eta1stEle_PAS") );
      CreateUserTH1D("h1_Phi1stEle_PAS_OutsideEta1stJetBump", getHistoNBins("Phi1stEle_PAS"), getHistoMin("Phi1stEle_PAS"), getHistoMax("Phi1stEle_PAS") );
      CreateUserTH1D("h1_Charge1stEle_PAS_OutsideEta1stJetBump", getHistoNBins("Charge1stEle_PAS"), getHistoMin("Charge1stEle_PAS"), getHistoMax("Charge1stEle_PAS") );
      CreateUserTH1D("h1_mDeltaPhiMETEle_OutsideEta1stJetBump", getHistoNBins("mDeltaPhiMETEle"), getHistoMin("mDeltaPhiMETEle"), getHistoMax("mDeltaPhiMETEle") );
      CreateUserTH1D("h1_mDeltaPhiEle1stJet_PAS_OutsideEta1stJetBump", getHistoNBins("mDeltaPhiEle1stJet_PAS"), getHistoMin("mDeltaPhiEle1stJet_PAS"), getHistoMax("mDeltaPhiEle1stJet_PAS") );
      CreateUserTH1D("h1_mDeltaPhiEle2ndJet_PAS_OutsideEta1stJetBump", getHistoNBins("mDeltaPhiEle2ndJet_PAS"), getHistoMin("mDeltaPhiEle2ndJet_PAS"), getHistoMax("mDeltaPhiEle2ndJet_PAS") );
      CreateUserTH1D("h1_mDeltaPhiMET2ndJet_OutsideEta1stJetBump", getHistoNBins("mDeltaPhiMET2ndJet"), getHistoMin("mDeltaPhiMET2ndJet"), getHistoMax("mDeltaPhiMET2ndJet") );
      CreateUserTH1D("h1_Ptenu_PAS_OutsideEta1stJetBump", getHistoNBins("Ptenu_PAS"), getHistoMin("Ptenu_PAS"), getHistoMax("Ptenu_PAS") );
      CreateUserTH1D("h1_1stJet_PTOverPTPlusMET_OutsideEta1stJetBump", 100, 0, 1.0 );
      CreateUserTH1D("h1_Conversion1stEle_OutsideEta1stJetBump", 2, -0.5, 1.5 );
      CreateUserTH1D("h1_MET_PAS_OutsideEta1stJetBump", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_CaloMET_PAS_OutsideEta1stJetBump", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_TCMET_PAS_OutsideEta1stJetBump", getHistoNBins("MET_PAS"), getHistoMin("MET_PAS"), getHistoMax("MET_PAS") );
      CreateUserTH1D("h1_MET_over_CaloMET_OutsideEta1stJetBump", 100, 0, 5 );
      CreateUserTH1D("h1_Njet_OutsideEta1stJetBump", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_NjetTCHELBTag_OutsideEta1stJetBump", getHistoNBins("nJet_PtCut_noOvrlp_ID"), getHistoMin("nJet_PtCut_noOvrlp_ID"), getHistoMax("nJet_PtCut_noOvrlp_ID") );
      CreateUserTH1D("h1_minDRej_OutsideEta1stJetBump", getHistoNBins("minDRej"), getHistoMin("minDRej"), getHistoMax("minDRej") );
      CreateUserTH1D("h1_maxDRej_OutsideEta1stJetBump", getHistoNBins("maxDRej"), getHistoMin("maxDRej"), getHistoMax("maxDRej") );
      CreateUserTH1D("h1_DRjets_OutsideEta1stJetBump", getHistoNBins("DeltaRjets_PAS"), getHistoMin("DeltaRjets_PAS"), getHistoMax("DeltaRjets_PAS") );
      CreateUserTH1D("h1_MTenu_PAS_OutsideEta1stJetBump", getHistoNBins("MTenu_PAS"), getHistoMin("MTenu_PAS"), getHistoMax("MTenu_PAS") );
    }

  CreateUserTH1D("h1_MTenu_PAS_EleBarrel", getHistoNBins("MTenu_PAS"), getHistoMin("MTenu_PAS"), getHistoMax("MTenu_PAS"));
  CreateUserTH1D("h1_minDRej_EleBarrel", getHistoNBins("minDRej"), getHistoMin("minDRej"), getHistoMax("minDRej"));
  CreateUserTH1D("h1_MTenu_PAS_EleEndcap", getHistoNBins("MTenu_PAS"), getHistoMin("MTenu_PAS"), getHistoMax("MTenu_PAS"));
  CreateUserTH1D("h1_minDRej_EleEndcap", getHistoNBins("minDRej"), getHistoMin("minDRej"), getHistoMax("minDRej"));

  CreateUserTH2D("h2_EToverPT_vs_ET_1stEle"
		 , getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS")
		 , 200, 0, 20
		 );

  // Pile-up vertex multiplicity distribution
  TH1D *h1_nPUVertex = new TH1D("h1_nPUVertex","h1_nPUVertex",11,-0.5,10.5);
  h1_nPUVertex->SetBinContent(1,850);
  h1_nPUVertex->SetBinContent(2,1123);
  h1_nPUVertex->SetBinContent(3,778);
  h1_nPUVertex->SetBinContent(4,386);
  h1_nPUVertex->SetBinContent(5,135);
  h1_nPUVertex->SetBinContent(6,36);
  h1_nPUVertex->SetBinContent(7,4);
  h1_nPUVertex->SetBinContent(8,1);
  h1_nPUVertex->SetBinError(1,29.15476);
  h1_nPUVertex->SetBinError(2,33.51119);
  h1_nPUVertex->SetBinError(3,27.89265);
  h1_nPUVertex->SetBinError(4,19.64688);
  h1_nPUVertex->SetBinError(5,11.61895);
  h1_nPUVertex->SetBinError(6,6);
  h1_nPUVertex->SetBinError(7,2);
  h1_nPUVertex->SetBinError(8,1);
  h1_nPUVertex->SetEntries(3313);

  // Random number generator
  TRandom3 *randomNumGen = new TRandom3;
  randomNumGen->SetSeed();

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
    std::auto_ptr<std::vector<int> >     JetPassID  ( new std::vector<int>()  );
    std::auto_ptr<std::vector<double> >  JetTCHE  ( new std::vector<double>()  );

    if(jetAlgorithm==1) //PF jets
      {
	for (int ijet=0 ; ijet< PFJetPt->size() ; ijet++)
	  {
	    JetPt->push_back( PFJetPt->at(ijet) );
	    JetPtRaw->push_back( PFJetPtRaw->at(ijet) );
	    JetEnergy->push_back( PFJetEnergy->at(ijet) );
	    JetEta->push_back( PFJetEta->at(ijet) );
	    JetPhi->push_back( PFJetPhi->at(ijet) );
            JetPassID->push_back( PFJetPassLooseID->at(ijet) );
// 	    JetPassID->push_back(
// 				 PFJetIdloose(PFJetChargedHadronEnergyFraction->at(ijet),
// 					      PFJetChargedEmEnergyFraction->at(ijet),
// 					      PFJetNeutralHadronEnergyFraction->at(ijet),
// 					      PFJetNeutralEmEnergyFraction->at(ijet),
// 					      PFJetEta->at(ijet) )
// 				 );
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
            JetPassID->push_back( CaloJetPassLooseID->at(ijet) );
// 	    JetPassID->push_back(
// 				 JetIdloose(CaloJetresEMF->at(ijet),
// 					    CaloJetfHPD->at(ijet),
// 					    CaloJetn90Hits->at(ijet),
// 					    CaloJetEta->at(ijet) )
// 				 );
            JetTCHE->push_back( CaloJetTrackCountingHighEffBTag->at(ijet) );
	  }//end loop over calo jets
      }//end if "calo jets"

    //## Define new met collection
    double thisMET;
    double thisMETPhi;
    double thisSumET;
    double thisGenSumET;

    if(isData==0)
      thisGenSumET = GenSumETTrue->at(0);
    else
      thisGenSumET = -1;

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
    // --> TCMET
    //     thisMET = TCMET->at(0);
    //     thisMETPhi = TCMETPhi->at(0);

    // MET smearing due to pile-up
    if( !isData && doPUMETSmearing ) {

      double nPU = (int)(h1_nPUVertex->GetRandom()+0.5);
      double thisMETx = thisMET*cos(thisMETPhi) + sqrt(nPU)*METxySigmaPerPU*randomNumGen->Gaus();
      double thisMETy = thisMET*sin(thisMETPhi) + sqrt(nPU)*METxySigmaPerPU*randomNumGen->Gaus();
      thisMET = sqrt( thisMETx*thisMETx + thisMETy*thisMETy );
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


    // The following data were processed but were declared as bad in
    // https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions10/7TeV/StreamExpress/Cert_132440-149442_7TeV_StreamExpress_Collisions10_JSON_v3.txt
    int passGoodRunList = 1;
    //     // Skip runs declared bad for HCAL reasons
    //     if ( run == 146511 ) passGoodRunList = 0;
    //     if ( run == 146513 ) passGoodRunList = 0;
    //     if ( run == 146514 ) passGoodRunList = 0;
    //     if ( run == 146644 ) passGoodRunList = 0;
    //     // Skip runs declared bad for other reasons
    //     if ( run == 141874 ) passGoodRunList = 0;
    //     if ( run == 141876 ) passGoodRunList = 0;
    //     if ( run == 142414 ) passGoodRunList = 0;
    //     if ( run == 147929 and ls == 619 ) passGoodRunList = 0;

    //## HLT
    int PassTrig = 0;
    int HLTFromRun[4] = {getPreCutValue1("HLTFromRun"),
			 getPreCutValue2("HLTFromRun"),
			 getPreCutValue3("HLTFromRun"),
			 getPreCutValue4("HLTFromRun")};
    int HLTTrigger[4] = {getPreCutValue1("HLTTrigger"),
			 getPreCutValue2("HLTTrigger"),
			 getPreCutValue3("HLTTrigger"),
			 getPreCutValue4("HLTTrigger")};
    int HLTTrgUsed;
    for (int i=0; i<4; i++) {
      if ( !isData && i != 0) continue; // For MC use HLTPhoton15 as the cleaned trigger is not in MC yet as of July 20, 2010
      if ( HLTFromRun[i] <= run ) {
 	//if(jentry == 0 ) STDOUT("run, i, HLTTrigger[i], HLTFromRun[i] = "<<run<<"\t"<<i<<"\t"<<"\t"<<HLTTrigger[i]<<"\t"<<HLTFromRun[i]);
	if (HLTTrigger[i] > 0 && HLTTrigger[i] < HLTResults->size() ) {
	  PassTrig=HLTResults->at(HLTTrigger[i]);
	  HLTTrgUsed=HLTTrigger[i];
	} else {
	  STDOUT("ERROR: HLTTrigger out of range of HLTResults: HLTTrigger = "<<HLTTrigger[i] <<"and HLTResults size = "<< HLTResults->size());
	}
      }
    }
    if(jentry == 0 ) STDOUT("Run = "<<run <<", HLTTrgUsed is number = "<<HLTTrgUsed<<" of the list HLTPathsOfInterest");


    //## Electrons
    vector<int> v_idx_ele_all;
    vector<int> v_idx_ele_PtCut;
    vector<int> v_idx_ele_PtCut_IDISO_noOverlap;
    int heepBitMask;

    //Loop over electrons
    for(int iele=0; iele<ElectronPt->size(); iele++)
      {

	// Reject ECAL spikes
	if ( 1 - ElectronSCS4S1->at(iele) > 0.95 ) continue;

	//no cut on reco electrons
	v_idx_ele_all.push_back(iele);

	//pT pre-cut on ele
	if( ElectronPt->at(iele) < ele_PtCut ) continue;
	v_idx_ele_PtCut.push_back(iele);

	// get heepBitMask for EB, GAP, EE
	if( fabs(ElectronEta->at(iele)) < eleEta_bar )
	  {
	    heepBitMask = heepBitMask_EB;
	  }
	else if ( fabs(ElectronEta->at(iele)) > eleEta_end_min && fabs(ElectronEta->at(iele)) < eleEta_end_max )
	  {
	    heepBitMask = heepBitMask_EE;
	  }
	else {
	  heepBitMask = heepBitMask_GAP;
	}

	//ID + ISO + NO overlap with good muons
	// int eleID = ElectronPassID->at(iele);
	// if ( (eleID & 1<<eleIDType) > 0  && ElectronOverlaps->at(iele)==0 )
	if ( (ElectronHeepID->at(iele) & ~heepBitMask)==0x0
	     // && ElectronOverlaps->at(iele)==0 //## + NO overlap with good muons (removed by default) ##
	     )
	  {
	    //STDOUT("ElectronHeepID = " << hex << ElectronHeepID->at(iele) << " ; ElectronPassID = " << ElectronPassID->at(iele) )
	    v_idx_ele_PtCut_IDISO_noOverlap.push_back(iele);
	  }

        if( plotEleIDIsoVar ) {

          FillUserTH1D("h1_ElectronPt_all", ElectronPt->at(iele) );
          FillUserTH1D("h1_ElectronEta_all", ElectronEta->at(iele) );
          FillUserTH1D("h1_ElectronPhi_all", ElectronPhi->at(iele) );
          FillUserTH1D("h1_ElectronDeltaEtaTrkSC_all", ElectronDeltaEtaTrkSC->at(iele) );
          FillUserTH1D("h1_ElectronDeltaPhiTrkSC_all", ElectronDeltaPhiTrkSC->at(iele) );
          FillUserTH1D("h1_ElectronHoE_all", ElectronHoE->at(iele) );
          FillUserTH1D("h1_ElectronSigmaIEtaIEta_all", ElectronSigmaIEtaIEta->at(iele) );
          FillUserTH1D("h1_ElectronEcalHcalIsoHeep_all",  ElectronEcalIsoHeep->at(iele)+ElectronHcalIsoD1Heep->at(iele) );
          FillUserTH1D("h1_ElectronEcalIsoHeep_all",  ElectronEcalIsoHeep->at(iele) );
          FillUserTH1D("h1_ElectronHcalIsoD1Heep_all",  ElectronHcalIsoD1Heep->at(iele) );
          FillUserTH1D("h1_ElectronHcalIsoD2Heep_all",  ElectronHcalIsoD2Heep->at(iele) );
          FillUserTH1D("h1_ElectronTrkIsoHeep_all",  ElectronTrkIsoHeep->at(iele) );
          FillUserTH1D("h1_ElectronE2x5OverE5x5_all",  ElectronE2x5OverE5x5->at(iele) );
          FillUserTH1D("h1_ElectronE1x5OverE5x5_all",  ElectronE1x5OverE5x5->at(iele) );
        }

      } // End loop over electrons


    //## Jets
    vector<int> v_idx_jet_all;
    vector<int> v_idx_jet_PtCut;
    vector<int> v_idx_jet_PtCut_noOverlap;
    vector<int> v_idx_jet_PtCut_noOverlap_ID;
    vector<int> v_idx_jet_PtCut_noOverlap_ID_EtaCut;
    vector<int> v_idx_jet_PtCut_noOverlap_ID_EtaCut_TCHEL;

    // Loop over jets
    for(int ijet=0; ijet<JetPt->size(); ijet++)
      {
	//no cut on reco jets
	v_idx_jet_all.push_back(ijet);

	//pT pre-cut on reco jets
	if ( JetPt->at(ijet) < jet_PtCut ) continue;
	v_idx_jet_PtCut.push_back(ijet);
      }

    vector <int> jetFlags(v_idx_jet_PtCut.size(), 0);
    int Njetflagged = 0;
    for (int iele=0; iele<v_idx_ele_PtCut_IDISO_noOverlap.size(); iele++)
      {
	TLorentzVector ele;
        ele.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
			 ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
			 ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),0);
	TLorentzVector jet;
	double minDR=9999.;
	int ijet_minDR = -1;
        for(int ijet=0; ijet<v_idx_jet_PtCut.size(); ijet++)
          {
	    if ( jetFlags[ijet] == 1 )
	      continue;
            jet.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut[ijet]),
			     JetEta->at(v_idx_jet_PtCut[ijet]),
			     JetPhi->at(v_idx_jet_PtCut[ijet]),0);
	    double DR = jet.DeltaR(ele);
	    if (DR<minDR)
	      {
		minDR = DR;
		ijet_minDR = ijet;
	      }
	  }
	if ( minDR < jet_ele_DeltaRcut && ijet_minDR > -1)
	  {
	    jetFlags[ijet_minDR] = 1;
	    Njetflagged++;
	  }
      }

//     // printouts for jet cleaning
//     STDOUT("CLEANING ----------- v_idx_ele_PtCut_IDISO_noOverlap.size = "<< v_idx_ele_PtCut_IDISO_noOverlap.size() <<", Njetflagged = "<< Njetflagged<<", diff="<< v_idx_ele_PtCut_IDISO_noOverlap.size()-Njetflagged );
//     if( (v_idx_ele_PtCut_IDISO_noOverlap.size()-Njetflagged) == 1 )
//       {
// 	TLorentzVector thisele;
// 	for(int iele=0; iele<v_idx_ele_PtCut_IDISO_noOverlap.size(); iele++)
// 	  {
// 	    thisele.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
// 				 ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
// 				 ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),0);
// 	    STDOUT("CLEANING: e"<<iele+1<<" Pt, eta, phi = "  << ", "<<thisele.Pt()<<", "<< thisele.Eta() <<", "<< thisele.Phi());
// 	  }
// 	TLorentzVector thisjet;
// 	for(int ijet=0; ijet<v_idx_jet_PtCut.size(); ijet++)
// 	  {
// 	    thisjet.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut[ijet]),
// 				 JetEta->at(v_idx_jet_PtCut[ijet]),
// 				 JetPhi->at(v_idx_jet_PtCut[ijet]),0);
// 	    STDOUT("CLEANING: j"<<ijet+1<<" Pt, eta, phi = " << ", "<<thisjet.Pt()<<", "<< thisjet.Eta() <<", "<< thisjet.Phi()<<" jetFlags="<<jetFlags[ijet] );
// 	  }
//       } // printouts for jet cleaning

    for(int ijet=0; ijet<v_idx_jet_PtCut.size(); ijet++) //pT pre-cut + no overlaps with electrons + jetID
      {
	bool passjetID = JetPassID->at(v_idx_jet_PtCut[ijet]);

	// ---- use the flag stored in rootTuples
	//if( (JetOverlaps->at(v_idx_jet_PtCut[ijet]) & 1 << eleIDType) == 0  /* NO overlap with electrons */
	// ----

	if( jetFlags[ijet] == 0  )                         /* NO overlap with electrons */
	  //  && passjetID == true )                            /* pass JetID */
	  // && (caloJetOverlaps[ijet] & 1 << 5)==0 )         /* NO overlap with muons */
	  v_idx_jet_PtCut_noOverlap.push_back(v_idx_jet_PtCut[ijet]);

	if( jetFlags[ijet] == 0                           /* NO overlap with electrons */
	    && passjetID == true )                            /* pass JetID */
	  // && (caloJetOverlaps[ijet] & 1 << 5)==0 )         /* NO overlap with muons */
	  v_idx_jet_PtCut_noOverlap_ID.push_back(v_idx_jet_PtCut[ijet]);

	if( jetFlags[ijet] == 0                           /* NO overlap with electrons */
	    && passjetID == true                             /* pass JetID */
	    && fabs( JetEta->at(v_idx_jet_PtCut[ijet]) ) < jet_EtaCut )
	  // && (caloJetOverlaps[ijet] & 1 << 5)==0 )         /* NO overlap with muons */
	  v_idx_jet_PtCut_noOverlap_ID_EtaCut.push_back(v_idx_jet_PtCut[ijet]);

        if( jetFlags[ijet] == 0                           /* NO overlap with electrons */
            && passjetID == true                             /* pass JetID */
            && fabs( JetEta->at(v_idx_jet_PtCut[ijet]) ) < jet_EtaCut
            && fabs( JetTCHE->at(v_idx_jet_PtCut[ijet]) ) > jet_TCHELCut ) /* TrackCountingHighEfficiency loose b-tag (see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagPerformance) */
          // && (caloJetOverlaps[ijet] & 1 << 5)==0 )         /* NO overlap with muons */
          v_idx_jet_PtCut_noOverlap_ID_EtaCut_TCHEL.push_back(v_idx_jet_PtCut[ijet]);

	//NOTE: We should verify that caloJetOverlaps match with the code above
      } // End loop over jets


    //## MET scale uncert.
    if( JetEnergyScale != 1 )
      { //use fix JES scaling passed from cut file

	TVector2 v_MET_old;
	TVector2 v_MET_new;

	//use only good jets (after electron-jet overlap) for re-doing MET
	for(int ijet=0; ijet<v_idx_jet_PtCut_noOverlap_ID.size(); ijet++)
	  {
	    TVector2 v_jet_pt_old;
	    TVector2 v_jet_pt_new;
	    v_jet_pt_old.SetMagPhi( JetPt->at(v_idx_jet_PtCut_noOverlap_ID[ijet])/JetEnergyScale , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[ijet]) );
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



    //## Muons
    vector<int> v_idx_muon_all;
    vector<int> v_idx_muon_PtCut;
    vector<int> v_idx_muon_PtCut_IDISO;

    // Loop over muons
    for(int imuon=0; imuon<MuonPt->size(); imuon++){

      // no cut on reco muons
      v_idx_muon_all.push_back(imuon);

      if ( (*MuonPt)[imuon] < muon_PtCut) continue;

      // pT pre-cut on muons
      v_idx_muon_PtCut.push_back(imuon);

      if ( ((*MuonTrkHits)[imuon]  >= muNHits_minThresh  )
	   &&( fabs((*MuonTrkD0)[imuon]) < muTrkD0Maximum )
	   &&((*MuonPassIso)[imuon]==1 )
	   &&((*MuonPassID)[imuon]==1) )
	{
	  v_idx_muon_PtCut_IDISO.push_back(imuon);
	}

    }// end loop over muons


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


    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();


    // Set the value of the variableNames listed in the cutFile to their current value

    // Trigger (L1 and HLT)
    if(isData==true)
      {
        fillVariableWithValue( "PassGoodRunList", passGoodRunList );
	fillVariableWithValue( "PassBPTX0", isBPTX0 ) ;
	fillVariableWithValue( "PassPhysDecl", isPhysDeclared ) ;
      }
    else
      {
        fillVariableWithValue( "PassGoodRunList", 1 );
	fillVariableWithValue( "PassBPTX0", true ) ;
	fillVariableWithValue( "PassPhysDecl", true ) ;
      }

    fillVariableWithValue( "PassHLT", PassTrig ) ;
    fillVariableWithValue( "nVertex", VertexChi2->size() ) ;
    fillVariableWithValue( "nVertex_PAS", VertexChi2->size() ) ;
    fillVariableWithValue( "nVertex_good", v_idx_vertex_good.size() ) ;
    fillVariableWithValue( "nVertex_good_PAS", v_idx_vertex_good.size() ) ;

    //Event filters at RECO level
    fillVariableWithValue( "PassBeamScraping", !isBeamScraping ) ;
    fillVariableWithValue( "PassPrimaryVertex", isPrimaryVertex ) ;

    //### FIXME ###
    //Spring10 ntuple production
    //fillVariableWithValue( "PassHBHENoiseFilter", passLooseNoiseFilter ) ;
    //38x ntuple production
    //fillVariableWithValue( "PassHBHENoiseFilter", passHBHENoiseFilter ) ;

    // nEle
    fillVariableWithValue( "nEle_all", v_idx_ele_all.size() ) ;
    fillVariableWithValue( "nEle_PtCut", v_idx_ele_PtCut.size() ) ;
    fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp", v_idx_ele_PtCut_IDISO_noOverlap.size() ) ;

    // nJet
    fillVariableWithValue( "nJet_all", v_idx_jet_all.size() ) ;
    fillVariableWithValue( "nJet_PtCut", v_idx_jet_PtCut.size() ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlp", v_idx_jet_PtCut_noOverlap.size() ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlp_ID", v_idx_jet_PtCut_noOverlap_ID.size() ) ;
    fillVariableWithValue( "nJet_WithJetEtaCut", v_idx_jet_PtCut_noOverlap_ID_EtaCut.size() ) ;
    fillVariableWithValue( "nJet_TCHELBTag", v_idx_jet_PtCut_noOverlap_ID_EtaCut_TCHEL.size() ) ;

    // nMuon
    fillVariableWithValue( "nMuon_PtCut_IDISO", v_idx_muon_PtCut_IDISO.size() ) ;
    //PAS Sept 2010
    fillVariableWithValue( "nMuon_PtCut_IDISO_PAS", v_idx_muon_PtCut_IDISO.size() ) ;

    // MET
    fillVariableWithValue("MET", thisMET);
    //PAS Sept 2010
    fillVariableWithValue("MET_PAS", thisMET);
    fillVariableWithValue("METPhi_PAS", thisMETPhi);

    //SUMET
    fillVariableWithValue("SumET", thisSumET);
    fillVariableWithValue("GenSumET", thisGenSumET);

    // Loop over GenParticles to calculate the GenMET
    //     TLorentzVector nu_p4 = 0.;
    //     for(int ipart=0; ipart<GenParticlePt->size(); ipart++)
    //       {
    //         //if the particle is not a neutrino, skip it
    //         if ( abs(GenParticlePdgId->at(ipart))!=12 &&
    //              abs(GenParticlePdgId->at(ipart))!=14 &&
    //              abs(GenParticlePdgId->at(ipart))!=16 ) continue;
    //         TLorentzVector temp_p4;
    //         temp_p4.SetPtEtaPhiM(GenParticlePt->at(ipart),
    //                              GenParticleEta->at(ipart),
    //                              GenParticlePhi->at(ipart),0.);
    //         nu_p4 += temp_p4;
    //       } // End loop over GenParticles
    //
    //     if( nu_p4.Perp()>0. ) fillVariableWithValue("deltaMET", (thisMET-nu_p4.Perp())/nu_p4.Perp());

    // 1st ele and transverse mass enu
    double MT, DeltaPhiMETEle = -999;
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 1 )
      {
	fillVariableWithValue( "Pt1stEle_IDISO_NoOvrlp", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	fillVariableWithValue( "Eta1stEle_IDISO_NoOvrlp", ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	fillVariableWithValue( "mEta1stEle_IDISO_NoOvrlp", fabs(ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0])) );
        fillVariableWithValue( "minMETPt1stEle", min(thisMET, ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0])) );
       	//PAS Sept 2010
	fillVariableWithValue( "Pt1stEle_PAS", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	fillVariableWithValue( "Eta1stEle_PAS", ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
        fillVariableWithValue( "minMETPt1stEle_PAS", min(thisMET, ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0])) );
	fillVariableWithValue( "Phi1stEle_PAS", ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
        fillVariableWithValue( "Vtxd01stEle_PAS", ElectronVtxDistXY->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
        fillVariableWithValue( "MissingHits1stEle_PAS", ElectronMissingHits->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
        fillVariableWithValue( "Dist1stEle_PAS", fabs(ElectronDist->at(v_idx_ele_PtCut_IDISO_noOverlap[0])) );
        fillVariableWithValue( "DCotTheta1stEle_PAS", fabs(ElectronDCotTheta->at(v_idx_ele_PtCut_IDISO_noOverlap[0])) );

	// DeltaPhi - MET vs 1st ele
	TVector2 v_MET;
	TVector2 v_ele;
	v_MET.SetMagPhi( thisMET , thisMETPhi);
	v_ele.SetMagPhi( ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) , ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	float deltaphi = v_MET.DeltaPhi(v_ele);
	fillVariableWithValue( "mDeltaPhiMETEle", fabs(deltaphi) );
	//PAS Sept 2010
	fillVariableWithValue( "mDeltaPhiMETEle_PAS", fabs(deltaphi) );
        DeltaPhiMETEle = fabs(deltaphi);

	// transverse mass enu
	MT = sqrt(2 * ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) * thisMET * (1 - cos(deltaphi)) );
	fillVariableWithValue("MTenu", MT);
	//PAS Sept 2010
	fillVariableWithValue("MTenu_PAS", MT);

	//PT(e,nu)
	TVector2 v_ele_MET;
	v_ele_MET = v_ele + v_MET;
	fillVariableWithValue("Ptenu_PAS", v_ele_MET.Mod());
	fillVariableWithValue("Ptenu", v_ele_MET.Mod());

	//Electron Charge
	fillVariableWithValue( "Charge1stEle_PAS", ElectronCharge->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
      }


    // 1st jet and deltaphi jet-MET
    double DeltaPhiMET1stJet = -999;
    double Mass1stJet = -999;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 1 )
      {
	fillVariableWithValue( "Pt1stJet_noOvrlp_ID", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "Eta1stJet_noOvrlp_ID", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "mEta1stJet_noOvrlp_ID", fabs(JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0])) );
	//PAS Sept 2010
	fillVariableWithValue( "Pt1stJet_PAS", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "Eta1stJet_PAS", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "Phi1stJet_PAS", JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
        fillVariableWithValue( "TCHE1stJet_PAS", JetTCHE->at(v_idx_jet_PtCut_noOverlap_ID[0]) );

	//DeltaPhi - MET vs 1st jet
	TVector2 v_MET;
	TVector2 v_jet;
	v_MET.SetMagPhi( 1 , thisMETPhi );
	v_jet.SetMagPhi( 1 , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	float deltaphi = v_MET.DeltaPhi(v_jet);
	fillVariableWithValue( "mDeltaPhiMET1stJet", fabs(deltaphi) );
	//PAS Sept 2010
	fillVariableWithValue( "mDeltaPhiMET1stJet_PAS", fabs(deltaphi) );
        DeltaPhiMET1stJet = fabs(deltaphi);

	if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 1 )
	  {
	    //distance from (pi,0) in (DeltaPhiMETj, DeltaPhiMETe) plane
	    double d1_DPhi_METe_METj = sqrt( pow(TMath::Pi() - DeltaPhiMET1stJet , 2)
					     + pow(DeltaPhiMETEle , 2) );
	    //distance from (0,pi) in (DeltaPhiMETj, DeltaPhiMETe) plane
	    double d2_DPhi_METe_METj = sqrt( pow(DeltaPhiMET1stJet , 2)
					     + pow( TMath::Pi() - DeltaPhiMETEle , 2) );

	    fillVariableWithValue( "d1_DPhi_METe_METj", d1_DPhi_METe_METj );
	    fillVariableWithValue( "d2_DPhi_METe_METj", d2_DPhi_METe_METj );
	  }

	//jet mass
	TLorentzVector thisjet;
	thisjet.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			     JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			     JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			     JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[0]));
	Mass1stJet = thisjet.M();
      }


    // 2nd jet and deltaphi jet-MET
    double Mass2ndJet = -999;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2 )
      {
	fillVariableWithValue( "Pt2ndJet_noOvrlp_ID", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "Eta2ndJet_noOvrlp_ID", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "mEta2ndJet_noOvrlp_ID", fabs(JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1])) );
	fillVariableWithValue( "maxMEtaJets_noOvrlp_ID", max( getVariableValue("mEta1stJet_noOvrlp_ID"), getVariableValue("mEta2ndJet_noOvrlp_ID") ) );
	//PAS Sept 2010
	fillVariableWithValue( "Pt2ndJet_PAS", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "Eta2ndJet_PAS", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "Phi2ndJet_PAS", JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
        fillVariableWithValue( "TCHE2ndJet_PAS", JetTCHE->at(v_idx_jet_PtCut_noOverlap_ID[1]) );

	//DeltaPhi - MET vs 2nd jet
	TVector2 v_MET;
	TVector2 v_jet;
	v_MET.SetMagPhi( 1 , thisMETPhi );
	v_jet.SetMagPhi( 1 , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	float deltaphi = v_MET.DeltaPhi(v_jet);
	fillVariableWithValue( "mDeltaPhiMET2ndJet", fabs(deltaphi) );
	//PAS Sept 2010
	fillVariableWithValue( "mDeltaPhiMET2ndJet_PAS", fabs(deltaphi) );

	//jet mass
	TLorentzVector thisjet;
	thisjet.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			     JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			     JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			     JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[1]));
	Mass2ndJet = thisjet.M();
      }

    // define "1ele" and "2jets" booleans
    bool OneEle=false;
    bool TwoJets=false;
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() == 1 ) OneEle = true;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2 ) TwoJets = true;

    // Mjj
    if (TwoJets)
      {
	TLorentzVector jet1, jet2, jj;
	jet1.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),0);
	jet2.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]),0);
	jj = jet1+jet2;


	//PAS June 2010
	fillVariableWithValue("Mjj_PAS", jj.M());
	fillVariableWithValue("DeltaRjets_PAS", jet1.DeltaR(jet2));
      }

    // ST
    if ( (OneEle) && (TwoJets) )
      {
        TVector2 v_ele;
        TVector2 v_jet1;
        TVector2 v_jet2;
        v_ele.SetMagPhi( 1 , ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
        v_jet1.SetMagPhi( 1 , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
        v_jet2.SetMagPhi( 1 , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
        //DeltaPhi - 1st ele vs 1st jet
        float deltaphi = v_ele.DeltaPhi(v_jet1);
        fillVariableWithValue( "mDeltaPhiEle1stJet_PAS", fabs(deltaphi) );
        //DeltaPhi - 1st ele vs 2nd jet
        deltaphi = v_ele.DeltaPhi(v_jet2);
        fillVariableWithValue( "mDeltaPhiEle2ndJet_PAS", fabs(deltaphi) );

	double calc_sT =
	  ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) +
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) +
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) +
	  thisMET;
        fillVariableWithValue("sT_presel", calc_sT);
	fillVariableWithValue("sT", calc_sT);
	fillVariableWithValue("sT_MLQ200", calc_sT);
	fillVariableWithValue("sT_MLQ250", calc_sT);
	fillVariableWithValue("sT_MLQ280", calc_sT);
	fillVariableWithValue("sT_MLQ300", calc_sT);
	fillVariableWithValue("sT_MLQ320", calc_sT);
        fillVariableWithValue("sT_MLQ340", calc_sT);
        fillVariableWithValue("sT_MLQ370", calc_sT);
        fillVariableWithValue("sT_MLQ400", calc_sT);
        fillVariableWithValue("sT_MLQ450", calc_sT);
        fillVariableWithValue("sT_MLQ500", calc_sT);
	//PAS Sept 2010
	fillVariableWithValue("sT_PAS", calc_sT);
      }

    // ST leptons (electron and MET)
    if (OneEle)
      {
	double calc_sTlep =
	  ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) +
	  thisMET;
	fillVariableWithValue("sTlep_PAS", calc_sTlep);
      }

    // ST jets
    if (TwoJets)
      {
	double calc_sTjet =
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) +
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]);
	fillVariableWithValue("sTjet_PAS", calc_sTjet);
      }

    // Mej , MTnuj
    double Me1j1, Me1j2, MTn1j1, MTn1j2 = -999;
    if ( (OneEle) && (TwoJets) )
      {
	//invariant mass electron-jet
	TLorentzVector jet1, jet2, ele1;
	ele1.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),0);
	jet1.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),0);
	jet2.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]),0);
	TLorentzVector jet1ele1, jet2ele1;
	jet1ele1 = jet1 + ele1;
	jet2ele1 = jet2 + ele1;
	Me1j1 = jet1ele1.M();
	Me1j2 = jet2ele1.M();

	double deltaR_e1j1 = ele1.DeltaR(jet1);
	double deltaR_e1j2 = ele1.DeltaR(jet2);

	fillVariableWithValue("minDRej", min(deltaR_e1j1,deltaR_e1j2) );
        fillVariableWithValue("maxDRej", max(deltaR_e1j1,deltaR_e1j2) );

	if( Me1j1 > Me1j2 )
	  {
	    fillVariableWithValue("Mej_1stPair", Me1j1);
	    fillVariableWithValue("Mej_2ndPair", Me1j2);
	    //PAS June 2010
	    fillVariableWithValue("Mej_1stPair_PAS", Me1j1);
	    fillVariableWithValue("Mej_2ndPair_PAS", Me1j2);
	  }
	else
	  {
	    fillVariableWithValue("Mej_1stPair", Me1j2);
	    fillVariableWithValue("Mej_2ndPair", Me1j1);
	    //PAS June 2010
	    fillVariableWithValue("Mej_1stPair_PAS", Me1j2);
	    fillVariableWithValue("Mej_2ndPair_PAS", Me1j1);
	  }

	//transverse mass neutrino-jet
	TVector2 v_MET;
	TVector2 v_jet1;
	TVector2 v_jet2;
	v_MET.SetMagPhi( 1 , thisMETPhi);
	v_jet1.SetMagPhi(1 , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]));
	v_jet2.SetMagPhi(1 , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]));
	float deltaphi1 = v_MET.DeltaPhi(v_jet1);
	float deltaphi2 = v_MET.DeltaPhi(v_jet2);
	MTn1j1 = sqrt(2 * thisMET * JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) * (1 - cos(deltaphi1)) );
	MTn1j2 = sqrt(2 * thisMET * JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) * (1 - cos(deltaphi2)) );


	if( MTn1j1 > MTn1j2 )
	  {
	    fillVariableWithValue("MTnuj_1stPair", MTn1j1);
	    fillVariableWithValue("MTnuj_2ndPair", MTn1j2);
	    //PAS June 2010
	    fillVariableWithValue("MTnuj_1stPair_PAS", MTn1j1);
	    fillVariableWithValue("MTnuj_2ndPair_PAS", MTn1j2);
	  }
	else
	  {
	    fillVariableWithValue("MTnuj_1stPair", MTn1j2);
	    fillVariableWithValue("MTnuj_2ndPair", MTn1j1);
	    //PAS June 2010
	    fillVariableWithValue("MTnuj_1stPair_PAS", MTn1j2);
	    fillVariableWithValue("MTnuj_2ndPair_PAS", MTn1j1);
	  }

      }


    // Evaluate cuts (but do not apply them)
    evaluateCuts();

    // Fill histograms and do analysis based on cut evaluation
    //h_nEleFinal->Fill( ElectronPt->size() );


    //INFO
    //      // retrieve value of previously filled variables (after making sure that they were filled)
    //      double totpTEle;
    //      if ( variableIsFilled("pT1stEle") && variableIsFilled("pT2ndEle") )
    //        totpTEle = getVariableValue("pT1stEle")+getVariableValue("pT2ndEle");
    //      // reject events that did not pass level 0 cuts
    //      if( !passedCut("0") ) continue;
    //      // ......

    //---------------------------------------------
    //------------ Gen level studies --------------
    //---------------------------------------------
    int Ngenleptons = 0;
    int Ngenelectrons = 0;
    int Ngenmuons = 0;
    int Ngentaus = 0;

    if(doGenLevelStudies)
      {
	vector<TLorentzVector> Neutrinos;
	vector<TLorentzVector> Ws;
	vector<TLorentzVector> Electrons;
	vector<TLorentzVector> Muons;
	vector<TLorentzVector> Taus;

	if(isData==0)
	  {
	    for (int genp=0; genp<GenParticlePdgId->size(); genp++)
	      {
		TLorentzVector tmp;

		//printout
// 			    cout << "idx: " << genp
// 				 << " GenParticlePdgId: " << GenParticlePdgId->at(genp)
// 				 << " GenParticleMotherIndex: " << GenParticleMotherIndex->at(genp)
// 				 << " GenParticleStatus: " << GenParticleStatus->at(genp)
//                                  << " GenParticlePt: " << GenParticlePt->at(genp)
// 				 << endl;

		//Neutrinos and Ws
		if( GenParticleStatus->at(genp)==1 &&
		    ( fabs(GenParticlePdgId->at(genp))==12 || fabs(GenParticlePdgId->at(genp))==14 || fabs(GenParticlePdgId->at(genp))==16 )
		    )
		  {
		    tmp.SetPtEtaPhiE( GenParticlePt->at(genp), GenParticleEta->at(genp), GenParticlePhi->at(genp), GenParticleEnergy->at(genp) );
		    Neutrinos.push_back(tmp);

		    //Ws
		    int idx_Wcandidate  =  GenParticleMotherIndex->at( GenParticleMotherIndex->at(genp) );
		    if( fabs( GenParticlePdgId->at(idx_Wcandidate) ) == 24 )
		      {
			tmp.SetPtEtaPhiE( GenParticlePt->at(idx_Wcandidate), GenParticleEta->at(idx_Wcandidate),
					  GenParticlePhi->at(idx_Wcandidate), GenParticleEnergy->at(idx_Wcandidate) );
			Ws.push_back(tmp);
		      }
		  }

		//Electrons
		if( GenParticleStatus->at(genp)==3 &&
		    ( fabs(GenParticlePdgId->at(genp))==11 )
		    )
		  {
		    tmp.SetPtEtaPhiE( GenParticlePt->at(genp), GenParticleEta->at(genp), GenParticlePhi->at(genp), GenParticleEnergy->at(genp) );
		    Electrons.push_back(tmp);

                    FillUserTH1D("h1_GenEle_Pt", GenParticlePt->at(genp) );
                    FillUserTH1D("h1_GenEle_Eta", GenParticleEta->at(genp) );
                    FillUserTH1D("h1_GenEle_Phi", GenParticlePhi->at(genp) );

		    //Loop over electrons
		    double minDeltaR = 999;
		    int index_minDeltaR = -1;
		    for(int iele=0; iele<ElectronPt->size(); iele++)
		      {
			TLorentzVector electron;
			electron.SetPtEtaPhiM(ElectronPt->at(iele),
					      ElectronEta->at(iele),
					      ElectronPhi->at(iele),0);
			double deltaR = electron.DeltaR(tmp);
			if(deltaR < minDeltaR )
			  {
			    minDeltaR = deltaR;
			    index_minDeltaR = iele;
			  }
		      }
		    if(index_minDeltaR!=-1)
		      {
			FillUserTH1D("h1_deltaR_RecoEle_GenEle", minDeltaR  );
			FillUserTH1D("h1_PtReco_over_PtGen_Ele",  ElectronPt->at(index_minDeltaR) / tmp.Pt() );
			FillUserTH2D("h2_PtReco_over_PtGen_vs_GenEle_Pt", tmp.Pt() , ElectronPt->at(index_minDeltaR) / tmp.Pt()  );
		      }
		  }

		//Muons
		if( GenParticleStatus->at(genp)==3 &&
		    ( fabs(GenParticlePdgId->at(genp))==13 )
		    )
		  {
		    tmp.SetPtEtaPhiE( GenParticlePt->at(genp), GenParticleEta->at(genp), GenParticlePhi->at(genp), GenParticleEnergy->at(genp) );
		    Muons.push_back(tmp);
		  }

		//Taus
		if( GenParticleStatus->at(genp)==3 &&
		    ( fabs(GenParticlePdgId->at(genp))==15 )
		    )
		  {
		    tmp.SetPtEtaPhiE( GenParticlePt->at(genp), GenParticleEta->at(genp), GenParticlePhi->at(genp), GenParticleEnergy->at(genp) );
		    Taus.push_back(tmp);
		  }

                //Electrons from LQ
                if( GenParticleStatus->at(genp)==3 &&
                    ( fabs(GenParticlePdgId->at(genp))==11 ) &&
                    ( fabs(GenParticlePdgId->at( GenParticleMotherIndex->at(genp) ))==42 )
                    )
                  {
                    FillUserTH1D("h1_LQGenEle_Pt", GenParticlePt->at(genp) );
                    FillUserTH1D("h1_LQGenEle_Eta", GenParticleEta->at(genp) );
                    FillUserTH1D("h1_LQGenEle_Phi", GenParticlePhi->at(genp) );
                  }

	      }//end loop over gen particles

	    //printout
	    //	cout << "---------------------" << endl;

	  }//end if statement: "is MC"

	//calculate "approximate" genMET
	TLorentzVector NeutrinoSum;
	for( int nu=0; nu<Neutrinos.size(); nu++)
	  {
	    NeutrinoSum += Neutrinos[nu];
	  }

	//number of gen particles
	FillUserTH1D("h_num_Neutrinos", Neutrinos.size() );
	FillUserTH1D("h_num_Ws", Ws.size() );
	Ngenleptons = Neutrinos.size();
	Ngenelectrons = Electrons.size();
	Ngenmuons = Muons.size();
	Ngentaus = Taus.size();

	//pfMET vs NeutrinosPt
	FillUserTH2D("h2_pfMET_vs_neutrinoPt", NeutrinoSum.Pt(), thisMET );

	//Delta Phi(MET,Neutrinos)
	TVector2 v_METreco;
	TVector2 v_neutrino;
	v_METreco.SetMagPhi( 1 , thisMETPhi);
	v_neutrino.SetMagPhi( 1 , NeutrinoSum.Phi() );
	double DeltaPhiMETNeutrino = v_METreco.DeltaPhi(v_neutrino);
	FillUserTH1D("h_DeltaPhi_pfMET_neutrino", fabs(DeltaPhiMETNeutrino) );

	//Ws Pt
	for(int wboson=0; wboson < Ws.size(); wboson++)
	  FillUserTH1D("h_WsPt", Ws[wboson].Pt() );

	if( passedAllPreviousCuts("sT_MLQ300") && passedCut("sT_MLQ300")
	    //&& variableIsFilled("d1_DPhi_METe_METj")
	    && variableIsFilled("mDeltaPhiMETEle_PAS")
	    && variableIsFilled("d2_DPhi_METe_METj")
	    )
	  {
	    //pfMET vs NeutrinosPt
	    FillUserTH2D("h2_pfMET_vs_neutrinoPt__sT", NeutrinoSum.Pt(), thisMET );
	    if( getVariableValue("mDeltaPhiMETEle_PAS") < 1.2 )
	      FillUserTH2D("h2_pfMET_vs_neutrinoPt_DphiMETeSmall__sT", NeutrinoSum.Pt(), thisMET );
	    else if( getVariableValue("d2_DPhi_METe_METj") < 0.6 )
	      FillUserTH2D("h2_pfMET_vs_neutrinoPt_d2Small__sT", NeutrinoSum.Pt(), thisMET );
	    else
	      FillUserTH2D("h2_pfMET_vs_neutrinoPt_d2DphiMETeLarge__sT", NeutrinoSum.Pt(), thisMET );

	    //Delta Phi(MET,Neutrinos)
	    FillUserTH1D("h_DeltaPhi_pfMET_neutrino__sT", fabs(DeltaPhiMETNeutrino) );
	    if( getVariableValue("mDeltaPhiMETEle_PAS") < 1.2 )
	      FillUserTH1D("h_DeltaPhi_pfMET_neutrino_DphiMETeSmall__sT", fabs(DeltaPhiMETNeutrino) );
	    else if( getVariableValue("d2_DPhi_METe_METj") < 0.6 )
	      FillUserTH1D("h_DeltaPhi_pfMET_neutrino_d2Small__sT", fabs(DeltaPhiMETNeutrino) );
	    else
	      FillUserTH1D("h_DeltaPhi_pfMET_neutrino_d2DphiMETeLarge__sT", fabs(DeltaPhiMETNeutrino) );

	    //Ws Pt
	    for(int wboson=0; wboson < Ws.size(); wboson++)
	      {
		FillUserTH1D("h_WsPt__sT", Ws[wboson].Pt() );
		if( getVariableValue("mDeltaPhiMETEle_PAS") < 1.2 )
		  FillUserTH1D("h_WsPt_DphiMETeSmall__sT", Ws[wboson].Pt() );
		else if( getVariableValue("d2_DPhi_METe_METj") < 0.6 )
		  FillUserTH1D("h_WsPt_d2Small__sT", Ws[wboson].Pt() );
		else
		  FillUserTH1D("h_WsPt_d2DphiMETeLarge__sT", Ws[wboson].Pt() );
	      }
	  }

      }//end if do gen level studies

    //---------------------------------------------
    //---------------------------------------------
    //---------------------------------------------

    //after pre-selection
    if( passedAllPreviousCuts("Pt1stEle_PAS")
	&& variableIsFilled("MTenu_PAS") && variableIsFilled("sT_PAS")
	&& variableIsFilled("mDeltaPhiMET2ndJet_PAS")
	&& variableIsFilled("Eta1stEle_PAS")
	&& variableIsFilled("Phi1stEle_PAS")
	&& variableIsFilled("Eta1stJet_PAS")
	&& variableIsFilled("Phi1stJet_PAS")
	&& variableIsFilled("Eta2ndJet_PAS")
	&& variableIsFilled("Phi2ndJet_PAS")
	&& variableIsFilled("METPhi_PAS")
	&& variableIsFilled("Pt1stEle_PAS")
	&& variableIsFilled("minDRej")
	)
      {
	if( ElectronCharge->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) > 0 )
	  {
	    FillUserTH1D("h1_MTenu_PAS_plus", getVariableValue("MTenu_PAS"));
	    FillUserTH1D("h1_sT_PAS_plus", getVariableValue("sT_PAS"));
	  }
	else
	  {
	    FillUserTH1D("h1_MTenu_PAS_minus", getVariableValue("MTenu_PAS"));
	    FillUserTH1D("h1_sT_PAS_minus", getVariableValue("sT_PAS"));
	  }

	if( fabs( getVariableValue("mDeltaPhiMET2ndJet_PAS"))<=1 )
	  {
	    FillUserTH1D("h1_MTenu_PAS_DeltaPhiMET2ndJet_0_1", getVariableValue("MTenu_PAS"));
	  }
	if( fabs( getVariableValue("mDeltaPhiMET2ndJet_PAS"))>1 &&
	    fabs( getVariableValue("mDeltaPhiMET2ndJet_PAS"))<2 )
	  {
	    FillUserTH1D("h1_MTenu_PAS_DeltaPhiMET2ndJet_1_2", getVariableValue("MTenu_PAS"));
	  }
	if( fabs( getVariableValue("mDeltaPhiMET2ndJet_PAS"))>=2 )
	  {
	    FillUserTH1D("h1_MTenu_PAS_DeltaPhiMET2ndJet_2_pi", getVariableValue("MTenu_PAS"));
	  }

	FillUserTH2D("h2_phi_VS_eta_1stEle", getVariableValue("Eta1stEle_PAS"), getVariableValue("Phi1stEle_PAS") );
	FillUserTH2D("h2_phi_VS_eta_1stJet", getVariableValue("Eta1stJet_PAS"), getVariableValue("Phi1stJet_PAS") );
	FillUserTH2D("h2_phi_VS_eta_2ndJet", getVariableValue("Eta2ndJet_PAS"), getVariableValue("Phi2ndJet_PAS") );

	//ele and MET phi for electrons in barrel/endcap
	if( fabs(getVariableValue("Eta1stEle_PAS")) <= eleEta_bar )
	  {//barrel
	    FillUserTH1D("h1_Phi1stEle_PAS_EleBarrel", getVariableValue("Phi1stEle_PAS"));
	    FillUserTH1D("h1_METPhi_PAS_EleBarrel", getVariableValue("METPhi_PAS"));
	    FillUserTH1D("h1_MTenu_PAS_EleBarrel", getVariableValue("MTenu_PAS") );
	  }
	else
	  {//endcap
	    FillUserTH1D("h1_Phi1stEle_PAS_EleEndcap", getVariableValue("Phi1stEle_PAS"));
	    FillUserTH1D("h1_METPhi_PAS_EleEndcap", getVariableValue("METPhi_PAS"));
	    FillUserTH1D("h1_MTenu_PAS_EleEndcap", getVariableValue("MTenu_PAS") );
	  }

	FillUserTH2D("h2_Phi1stEle_vs_METPhi", getVariableValue("METPhi_PAS") , getVariableValue("Phi1stEle_PAS") );
	FillUserTH2D("h2_Phi1stEle_vs_PtEleOverST", getVariableValue("Pt1stEle_PAS")/getVariableValue("sT_PAS") , getVariableValue("Phi1stEle_PAS") );
	FillUserTH2D("h2_METPhi_vs_PtEleOverST", getVariableValue("Pt1stEle_PAS")/getVariableValue("sT_PAS") , getVariableValue("METPhi_PAS") );

	//HEEP electron id/isolation variables
	int myEle = v_idx_ele_PtCut_IDISO_noOverlap[0];

	//high MT tails
	if( getVariableValue("MTenu")>MTenu_Thresh2 && plotEleIDIsoVar )
	  {
	    FillUserTH1D("h1_ElectronDeltaEtaTrkSC_highMT", ElectronDeltaEtaTrkSC->at(myEle) );
	    FillUserTH1D("h1_ElectronDeltaPhiTrkSC_highMT", ElectronDeltaPhiTrkSC->at(myEle) );
	    FillUserTH1D("h1_ElectronHoE_highMT", ElectronHoE->at(myEle) );
	    FillUserTH1D("h1_ElectronSigmaIEtaIEta_highMT", ElectronSigmaIEtaIEta->at(myEle) );
	    FillUserTH1D("h1_ElectronEcalHcalIsoHeep_highMT",  ElectronEcalIsoHeep->at(myEle)+ElectronHcalIsoD1Heep->at(myEle) );
	    FillUserTH1D("h1_ElectronHcalIsoD2Heep_highMT",  ElectronHcalIsoD2Heep->at(myEle) );
	    FillUserTH1D("h1_ElectronTrkIsoHeep_highMT",  ElectronTrkIsoHeep->at(myEle) );
	    FillUserTH1D("h1_ElectronE2x5OverE5x5_highMT",  ElectronE2x5OverE5x5->at(myEle) );
	    FillUserTH1D("h1_ElectronE1x5OverE5x5_highMT",  ElectronE1x5OverE5x5->at(myEle) );

	    if( fabs( ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) ) <= eleEta_bar )
	      {
		FillUserTH1D("h1_ElectronDeltaEtaTrkSC_barrel_highMT", ElectronDeltaEtaTrkSC->at(myEle) );
		FillUserTH1D("h1_ElectronDeltaPhiTrkSC_barrel_highMT", ElectronDeltaPhiTrkSC->at(myEle) );
		FillUserTH1D("h1_ElectronHoE_barrel_highMT", ElectronHoE->at(myEle) );
		FillUserTH1D("h1_ElectronSigmaIEtaIEta_barrel_highMT", ElectronSigmaIEtaIEta->at(myEle) );
		FillUserTH1D("h1_ElectronEcalHcalIsoHeep_barrel_highMT",  ElectronEcalIsoHeep->at(myEle)+ElectronHcalIsoD1Heep->at(myEle) );
		FillUserTH1D("h1_ElectronHcalIsoD2Heep_barrel_highMT",  ElectronHcalIsoD2Heep->at(myEle) );
		FillUserTH1D("h1_ElectronTrkIsoHeep_barrel_highMT",  ElectronTrkIsoHeep->at(myEle) );
		FillUserTH1D("h1_ElectronE2x5OverE5x5_barrel_highMT",  ElectronE2x5OverE5x5->at(myEle) );
		FillUserTH1D("h1_ElectronE1x5OverE5x5_barrel_highMT",  ElectronE1x5OverE5x5->at(myEle) );
	      }
	    else
	      {
		FillUserTH1D("h1_ElectronDeltaEtaTrkSC_endcap_highMT", ElectronDeltaEtaTrkSC->at(myEle) );
		FillUserTH1D("h1_ElectronDeltaPhiTrkSC_endcap_highMT", ElectronDeltaPhiTrkSC->at(myEle) );
		FillUserTH1D("h1_ElectronHoE_endcap_highMT", ElectronHoE->at(myEle) );
		FillUserTH1D("h1_ElectronSigmaIEtaIEta_endcap_highMT", ElectronSigmaIEtaIEta->at(myEle) );
		FillUserTH1D("h1_ElectronEcalHcalIsoHeep_endcap_highMT",  ElectronEcalIsoHeep->at(myEle)+ElectronHcalIsoD1Heep->at(myEle) );
		FillUserTH1D("h1_ElectronHcalIsoD2Heep_endcap_highMT",  ElectronHcalIsoD2Heep->at(myEle) );
		FillUserTH1D("h1_ElectronTrkIsoHeep_endcap_highMT",  ElectronTrkIsoHeep->at(myEle) );
		FillUserTH1D("h1_ElectronE2x5OverE5x5_endcap_highMT",  ElectronE2x5OverE5x5->at(myEle) );
		FillUserTH1D("h1_ElectronE1x5OverE5x5_endcap_highMT",  ElectronE1x5OverE5x5->at(myEle) );
	      }

// 	    //XXXXXXXX DEBUG XXXXXXXXXX
// 	    CreateAndFillUserTH1D("h1_ElectronRelIso_highMT", 1000, 0, 1, ElectronRelIso->at(myEle) );
// 	    CreateAndFillUserTH1D("h1_ElectronVtxDistXY_highMT", 200, -0.01, 0.01, ElectronVtxDistXY->at(myEle) );
// 	    CreateAndFillUserTH1D("h1_ElectronVtxDistZ_highMT", 200, -0.1, 0.1, ElectronVtxDistZ->at(myEle) );
// 	    CreateAndFillUserTH1D("h1_ElectronVtxIndex_highMT", 10, 0, 10, ElectronVtxIndex->at(myEle) );
// 	    CreateAndFillUserTH1D("h1_ElectronCharge_highMT", 2, -1.001, 1.001, ElectronCharge->at(myEle) );
// 	    CreateAndFillUserTH1D("h1_ElectronMissingHits_highMT", 10, 0, 10, ElectronMissingHits->at(myEle) );
// 	    CreateAndFillUserTH1D("h1_ElectronSCkOutOfTime_highMT", 2, 0, 2, ElectronSCkOutOfTime->at(myEle) );

// 	    CreateAndFillUserTH1D("h1_ElectronET_highMT", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"), getVariableValue("Pt1stEle_PAS") );
// 	    CreateAndFillUserTH1D("h1_ElectronEta_highMT", getHistoNBins("Eta1stEle_PAS"), getHistoMin("Eta1stEle_PAS"), getHistoMax("Eta1stEle_PAS"), getVariableValue("Eta1stEle_PAS") );
// 	    CreateAndFillUserTH1D("h1_ElectronEToverPT_highMT", 200, 0, 20, double ( getVariableValue("Pt1stEle_PAS") / ElectronTrackPt->at(myEle) ) );
// 	    CreateAndFillUserTH1D("h1_nJet_PtCut_noOvrlp_ID_highMT", 20, 0, 20, getVariableValue("nJet_PtCut_noOvrlp_ID") );
// 	    CreateAndFillUserTH1D("h1_mDeltaPhiMET2ndJet_PAS_highMT", getHistoNBins("mDeltaPhiMET2ndJet_PAS"), getHistoMin("mDeltaPhiMET2ndJet_PAS"), getHistoMax("mDeltaPhiMET2ndJet_PAS"),
// 				  getVariableValue("mDeltaPhiMET2ndJet_PAS") );
// 	    CreateAndFillUserTH1D("h1_sT_highMT", getHistoNBins("sT"), getHistoMin("sT"), getHistoMax("sT"), getVariableValue("sT") );
// 	    CreateAndFillUserTH1D("h1_MET_highMT", getHistoNBins("MET"), getHistoMin("MET"), getHistoMax("MET"), getVariableValue("MET") );

// 	    if( getVariableValue("mDeltaPhiMET2ndJet_PAS") < 0.5 )
// 	      {
// 		CreateAndFillUserTH1D("h1_2ndJet_PTOverPTPlusMET_highMT", 100, 0, 1, double (getVariableValue("Pt2ndJet_PAS") / (getVariableValue("Pt2ndJet_PAS") + getVariableValue("MET")) ) );
// 	      }
// 	    //XXXXXXXXXXXXXXXXXXXXXXXXX
	  }

// 	//XXXXXXXX DEBUG XXXXXXXXXX
// 	//high PT tails
// 	if( getVariableValue("Pt1stEle_PAS")>Pt1stEle_PAS_Thresh2 && plotEleIDIsoVar )
// 	  {
// 	    CreateAndFillUserTH1D("h1_ElectronEcalHcalIsoHeep_highPTe", 500,0,100, ElectronEcalIsoHeep->at(myEle)+ElectronHcalIsoD1Heep->at(myEle) );
// 	    CreateAndFillUserTH1D("h1_ElectronHcalIsoD2Heep_highPTe", 200,0,100,  ElectronHcalIsoD2Heep->at(myEle) );
// 	    CreateAndFillUserTH1D("h1_ElectronTrkIsoHeep_highPTe",200,0,100,  ElectronTrkIsoHeep->at(myEle) );

// 	    CreateAndFillUserTH1D("h1_ElectronRelIso_highPTe", 1000, 0, 1, ElectronRelIso->at(myEle) );
// 	    CreateAndFillUserTH1D("h1_ElectronVtxDistXY_highPTe", 200, -0.01, 0.01, ElectronVtxDistXY->at(myEle) );
// 	    CreateAndFillUserTH1D("h1_ElectronVtxDistZ_highPTe", 200, -0.1, 0.1, ElectronVtxDistZ->at(myEle) );
// 	    CreateAndFillUserTH1D("h1_ElectronVtxIndex_highPTe", 10, 0, 10, ElectronVtxIndex->at(myEle) );
// 	    CreateAndFillUserTH1D("h1_ElectronCharge_highPTe", 2, -1.001, 1.001, ElectronCharge->at(myEle) );
// 	    CreateAndFillUserTH1D("h1_ElectronMissingHits_highPTe", 10, 0, 10, ElectronMissingHits->at(myEle) );
// 	    CreateAndFillUserTH1D("h1_ElectronSCkOutOfTime_highPTe", 2, 0, 2, ElectronSCkOutOfTime->at(myEle) );

// 	    CreateAndFillUserTH1D("h1_ElectronET_highPTe", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS"), getVariableValue("Pt1stEle_PAS") );
// 	    CreateAndFillUserTH1D("h1_ElectronEta_highPTe", getHistoNBins("Eta1stEle_PAS"), getHistoMin("Eta1stEle_PAS"), getHistoMax("Eta1stEle_PAS"), getVariableValue("Eta1stEle_PAS") );
// 	    CreateAndFillUserTH1D("h1_ElectronEToverPT_highPTe", 200, 0, 20, double ( getVariableValue("Pt1stEle_PAS") / ElectronTrackPt->at(myEle) ) );
// 	    CreateAndFillUserTH1D("h1_nJet_PtCut_noOvrlp_ID_highPTe", 20, 0, 20, getVariableValue("nJet_PtCut_noOvrlp_ID") );
// 	    CreateAndFillUserTH1D("h1_mDeltaPhiMET2ndJet_PAS_highPTe", getHistoNBins("mDeltaPhiMET2ndJet_PAS"), getHistoMin("mDeltaPhiMET2ndJet_PAS"), getHistoMax("mDeltaPhiMET2ndJet_PAS"),
// 				  getVariableValue("mDeltaPhiMET2ndJet_PAS") );
// 	    CreateAndFillUserTH1D("h1_sT_highPTe", getHistoNBins("sT"), getHistoMin("sT"), getHistoMax("sT"), getVariableValue("sT") );
// 	    CreateAndFillUserTH1D("h1_MET_highPTe", getHistoNBins("MET"), getHistoMin("MET"), getHistoMax("MET"), getVariableValue("MET") );

// 	    if( getVariableValue("mDeltaPhiMET2ndJet_PAS") < 0.5 )
// 	      {
// 		CreateAndFillUserTH1D("h1_2ndJet_PTOverPTPlusMET_highPTe", 100, 0, 1, double (getVariableValue("Pt2ndJet_PAS") / (getVariableValue("Pt2ndJet_PAS") + getVariableValue("MET")) ) );
// 	      }

// 	  }
// 	//XXXXXXXXXXXXXXXXXXXXXXXXX

	//ET/PT vs ET
	double EToverPT = ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) / ElectronTrackPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) ;
	FillUserTH2D("h2_EToverPT_vs_ET_1stEle", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]), EToverPT );

      }

    if( passedAllPreviousCuts("d1_DPhi_METe_METj")
	&& variableIsFilled("MTenu_PAS") && variableIsFilled("Eta1stEle_PAS")
	&& variableIsFilled("mDeltaPhiMETEle_PAS") && variableIsFilled("MET_PAS")
	&& variableIsFilled("mDeltaPhiMET1stJet_PAS") && variableIsFilled("mDeltaPhiMET2ndJet_PAS")
	&& variableIsFilled("sT_PAS")
	)
      {
	FillUserTH2D("h2_MTnuj_vs_MET", getVariableValue("MET_PAS"), getVariableValue("MTenu_PAS") );
	FillUserTH2D("h2_ST_vs_MET", getVariableValue("MET_PAS"), getVariableValue("sT_PAS") );
	FillUserTH2D("h2_ST_vs_MTnuj", getVariableValue("MTenu_PAS"), getVariableValue("sT_PAS") );
	FillUserTH2D("h2_DeltaPhiMETEle_vs_MET", getVariableValue("MET_PAS"),
		     fabs( getVariableValue("mDeltaPhiMETEle_PAS")) );
	FillUserTH2D("h2_DeltaPhiMET1stJet_vs_MET", getVariableValue("MET_PAS"),
		     fabs( getVariableValue("mDeltaPhiMET1stJet_PAS")) );
	FillUserTH2D("h2_DeltaPhiMET2ndJet_vs_MET", getVariableValue("MET_PAS") ,
		     fabs( getVariableValue("mDeltaPhiMET2ndJet_PAS")) );

	FillUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet",
		     fabs( getVariableValue("mDeltaPhiMET1stJet_PAS") ),
		     fabs( getVariableValue("mDeltaPhiMETEle_PAS")) );
	if( ElectronCharge->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) > 0)
	  {
	    FillUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet_plus",
			 fabs( getVariableValue("mDeltaPhiMET1stJet_PAS") ),
			 fabs( getVariableValue("mDeltaPhiMETEle_PAS")) );
	  }
	else
	  {
	    FillUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet_minus",
			 fabs( getVariableValue("mDeltaPhiMET1stJet_PAS") ),
			 fabs( getVariableValue("mDeltaPhiMETEle_PAS")) );
	  }

	FillUserTH2D("h2_DeltaPhiMETEle_vs_MET2ndJet",
		     fabs( getVariableValue("mDeltaPhiMET2ndJet_PAS")), fabs( getVariableValue("mDeltaPhiMETEle_PAS")) );
	FillUserTH2D("h2_DeltaPhiMET2ndJet_vs_MET1stJet",
		     fabs( getVariableValue("mDeltaPhiMET1stJet_PAS") ), fabs( getVariableValue("mDeltaPhiMET2ndJet_PAS")) );
	FillUserTH2D("h2_MT_vs_etaEle", getVariableValue("Eta1stEle_PAS") , getVariableValue("MTenu_PAS") );
      }

    if( passedAllPreviousCuts("minMETPt1stEle") && passedCut("minMETPt1stEle")
	&& variableIsFilled("mDeltaPhiMETEle_PAS")
	&& variableIsFilled("mDeltaPhiMET1stJet_PAS")
	)
      {
	FillUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__minMETpTe",
		     fabs( getVariableValue("mDeltaPhiMET1stJet_PAS") ),
		     fabs( getVariableValue("mDeltaPhiMETEle_PAS")) );
	if( ElectronCharge->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) > 0)
	  {
	    FillUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__minMETpTe_plus",
			 fabs( getVariableValue("mDeltaPhiMET1stJet_PAS") ),
			 fabs( getVariableValue("mDeltaPhiMETEle_PAS")) );
	  }
	else
	  {
	    FillUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__minMETpTe_minus",
			 fabs( getVariableValue("mDeltaPhiMET1stJet_PAS") ),
			 fabs( getVariableValue("mDeltaPhiMETEle_PAS")) );
	  }
      }

    if( passedAllPreviousCuts("MTenu") && passedCut("MTenu")
	&& variableIsFilled("mDeltaPhiMETEle_PAS")
	&& variableIsFilled("mDeltaPhiMET1stJet_PAS")
	)
      {
	FillUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__MTenu",
		     fabs( getVariableValue("mDeltaPhiMET1stJet_PAS") ),
		     fabs( getVariableValue("mDeltaPhiMETEle_PAS")) );
	if( ElectronCharge->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) > 0)
	  {
	    FillUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__MTenu_plus",
			 fabs( getVariableValue("mDeltaPhiMET1stJet_PAS") ),
			 fabs( getVariableValue("mDeltaPhiMETEle_PAS")) );
	  }
	else
	  {
	    FillUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__MTenu_minus",
			 fabs( getVariableValue("mDeltaPhiMET1stJet_PAS") ),
			 fabs( getVariableValue("mDeltaPhiMETEle_PAS")) );
	  }
      }

    if( passedAllPreviousCuts("sT_MLQ300") && passedCut("sT_MLQ300")
	&& variableIsFilled("mDeltaPhiMETEle_PAS")
	&& variableIsFilled("mDeltaPhiMET1stJet_PAS")
	&& variableIsFilled("mDeltaPhiMET2ndJet_PAS")
	&& variableIsFilled("Charge1stEle_PAS")
	&& variableIsFilled("Eta1stEle_PAS")
	)
      {
	FillUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__sT",
		     fabs( getVariableValue("mDeltaPhiMET1stJet_PAS") ),
		     fabs( getVariableValue("mDeltaPhiMETEle_PAS")) );

	FillUserTH2D("h2_DeltaPhiMETEle_vs_MET2ndJet__sT",
		     fabs( getVariableValue("mDeltaPhiMET2ndJet_PAS") ),
		     fabs( getVariableValue("mDeltaPhiMETEle_PAS")) );

	FillUserTH2D("h2_DeltaPhiMET2ndJet_vs_MET1stJet__sT",
		     fabs( getVariableValue("mDeltaPhiMET1stJet_PAS") ),
		     fabs( getVariableValue("mDeltaPhiMET2ndJet_PAS")) );

	if( ElectronCharge->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) > 0)
	  {
	    FillUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__sT_plus",
			 fabs( getVariableValue("mDeltaPhiMET1stJet_PAS") ),
			 fabs( getVariableValue("mDeltaPhiMETEle_PAS")) );
	  }
	else
	  {
	    FillUserTH2D("h2_DeltaPhiMETEle_vs_MET1stJet__sT_minus",
			 fabs( getVariableValue("mDeltaPhiMET1stJet_PAS") ),
			 fabs( getVariableValue("mDeltaPhiMETEle_PAS")) );
	  }

	FillUserTH1D("h1_Charge1stEle_PAS__sT", getVariableValue("Charge1stEle_PAS"));
	FillUserTH1D("h1_Eta1stEle_PAS__sT", getVariableValue("Eta1stEle_PAS"));

      }

    //minDeltaR for barrel and endcaps
    if( passedAllPreviousCuts("minDRej")
	&& variableIsFilled("minDRej")
	&& variableIsFilled("Eta1stEle_PAS")
	)
      {
	if( fabs(getVariableValue("Eta1stEle_PAS")) <= eleEta_bar )
	  {//barrel
	    FillUserTH1D("h1_minDRej_EleBarrel", getVariableValue("minDRej") );
	  }
	else
	  {//endcap
	    FillUserTH1D("h1_minDRej_EleEndcap", getVariableValue("minDRej") );
	  }
      }

    // Events passing all non-sT cuts
    if( passedAllPreviousCuts("sT_presel")
        && passedCut("nMuon_PtCut_IDISO")
        && passedCut("MTenu")
        && passedCut("minMETPt1stEle")
        )
      {
        FillUserTH1D("h1_sT_AllNonSTCuts", getVariableValue("sT_PAS"));

        if( isData ) {
          STDOUT("PassFullSelection: ----------- START ------------");

          STDOUT("PassFullSelection: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
          if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") )
            STDOUT("PassFullSelection: Pt1stEle_PAS,Eta1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS"));
          if( variableIsFilled("MET_PAS") )
            STDOUT("PassFullSelection: MET_PAS = "<<getVariableValue("MET_PAS"));
          if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") )
            STDOUT("PassFullSelection: Pt1stJet_PAS,Eta1stJet_PAS = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS"));
          if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") )
            STDOUT("PassFullSelection: Pt2ndJet_PAS,Eta2ndJet_PAS = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS"));
          if( variableIsFilled("MTenu_PAS") && variableIsFilled("Mjj_PAS") )
            STDOUT("PassFullSelection: MTenu_PAS,Mjj_PAS = "<<getVariableValue("MTenu_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
          if( variableIsFilled("sT_PAS") )
            STDOUT("PassFullSelection: sT_PAS = "<<getVariableValue("sT_PAS"));
          if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )
            STDOUT("PassFullSelection: Mej_1stPair_PAS,Mej_2ndPair_PAS = "
                   <<getVariableValue("Mej_1stPair_PAS")
                   <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
          if( variableIsFilled("MTnuj_1stPair_PAS") && variableIsFilled("MTnuj_2ndPair_PAS") )
            STDOUT("PassFullSelection: MTnuj_1stPair_PAS,MTnuj_2ndPair_PAS = "
                   <<getVariableValue("MTnuj_1stPair_PAS")
                   <<",\t"<<getVariableValue("MTnuj_2ndPair_PAS"));
          if( variableIsFilled("mDeltaPhiMETEle_PAS") && variableIsFilled("mDeltaPhiMET1stJet_PAS")
              && variableIsFilled("mDeltaPhiMET2ndJet_PAS") )
            STDOUT("PassFullSelection: mDeltaPhiMETEle_PAS,mDeltaPhiMET1stJet_PAS,mDeltaPhiMET2ndJet_PAS = "
                   <<getVariableValue("mDeltaPhiMETEle_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET1stJet_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET2ndJet_PAS") );
          if( variableIsFilled("nMuon_PtCut_IDISO") )
            STDOUT("PassFullSelection: nMuon_PtCut_IDISO = "<<getVariableValue("nMuon_PtCut_IDISO"));
          if( variableIsFilled("minDRej") )
            STDOUT("PassFullSelection: minDRej = "<<getVariableValue("minDRej"));


          STDOUT("PassFullSelection: ------------ END -------------");
        }
      }

    // Events in the tails of the MET distribution that pass event pre-selection
    if( isData && variableIsFilled("MET")
        && passedAllPreviousCuts("sT_presel") && passedCut("sT_presel")
      )
      {
        if( getVariableValue("MET")>met_Thresh ) {
          STDOUT("PassMETThreshold: ----------- START ------------");

          STDOUT("PassMETThreshold: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
          if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") )
            STDOUT("PassMETThreshold: Pt1stEle_PAS,Eta1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS"));
          if( variableIsFilled("MET_PAS") )
            STDOUT("PassMETThreshold: MET_PAS = "<<getVariableValue("MET_PAS"));
          if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") )
            STDOUT("PassMETThreshold: Pt1stJet_PAS,Eta1stJet_PAS = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS"));
          if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") )
            STDOUT("PassMETThreshold: Pt2ndJet_PAS,Eta2ndJet_PAS = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS"));
          if( variableIsFilled("MTenu_PAS") && variableIsFilled("Mjj_PAS") )
            STDOUT("PassMETThreshold: MTenu_PAS,Mjj_PAS = "<<getVariableValue("MTenu_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
          if( variableIsFilled("sT_PAS") )
            STDOUT("PassMETThreshold: sT_PAS = "<<getVariableValue("sT_PAS"));
          if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )
            STDOUT("PassMETThreshold: Mej_1stPair_PAS,Mej_2ndPair_PAS = "
                   <<getVariableValue("Mej_1stPair_PAS")
                   <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
          if( variableIsFilled("MTnuj_1stPair_PAS") && variableIsFilled("MTnuj_2ndPair_PAS") )
            STDOUT("PassMETThreshold: MTnuj_1stPair_PAS,MTnuj_2ndPair_PAS = "
                   <<getVariableValue("MTnuj_1stPair_PAS")
                   <<",\t"<<getVariableValue("MTnuj_2ndPair_PAS"));
          if( variableIsFilled("mDeltaPhiMETEle_PAS") && variableIsFilled("mDeltaPhiMET1stJet_PAS")
              && variableIsFilled("mDeltaPhiMET2ndJet_PAS") )
            STDOUT("PassMETThreshold: mDeltaPhiMETEle_PAS,mDeltaPhiMET1stJet_PAS,mDeltaPhiMET2ndJet_PAS = "
                   <<getVariableValue("mDeltaPhiMETEle_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET1stJet_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET2ndJet_PAS") );
          if( variableIsFilled("nMuon_PtCut_IDISO") )
            STDOUT("PassMETThreshold: nMuon_PtCut_IDISO = "<<getVariableValue("nMuon_PtCut_IDISO"));
          if( variableIsFilled("minDRej") )
            STDOUT("PassMETThreshold: minDRej = "<<getVariableValue("minDRej"));


          STDOUT("PassMETThreshold: ------------ END -------------");
        }
      }

    // Events in the tails of the minMETPt1stEle distribution that pass event pre-selection
    if( isData && variableIsFilled("minMETPt1stEle")
        && passedAllPreviousCuts("sT_presel") && passedCut("sT_presel")
      )
      {
        if( getVariableValue("minMETPt1stEle")>minMETPt1stEle_Thresh ) {
          STDOUT("PassMinMETPt1stEleThreshold: ----------- START ------------");

          STDOUT("PassMinMETPt1stEleThreshold: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
          if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") )
            STDOUT("PassMinMETPt1stEleThreshold: Pt1stEle_PAS,Eta1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS"));
          if( variableIsFilled("MET_PAS") )
            STDOUT("PassMinMETPt1stEleThreshold: MET_PAS = "<<getVariableValue("MET_PAS"));
          if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") )
            STDOUT("PassMinMETPt1stEleThreshold: Pt1stJet_PAS,Eta1stJet_PAS = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS"));
          if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") )
            STDOUT("PassMinMETPt1stEleThreshold: Pt2ndJet_PAS,Eta2ndJet_PAS = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS"));
          if( variableIsFilled("MTenu_PAS") && variableIsFilled("Mjj_PAS") )
            STDOUT("PassMinMETPt1stEleThreshold: MTenu_PAS,Mjj_PAS = "<<getVariableValue("MTenu_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
          if( variableIsFilled("sT_PAS") )
            STDOUT("PassMinMETPt1stEleThreshold: sT_PAS = "<<getVariableValue("sT_PAS"));
          if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )
            STDOUT("PassMinMETPt1stEleThreshold: Mej_1stPair_PAS,Mej_2ndPair_PAS = "
                   <<getVariableValue("Mej_1stPair_PAS")
                   <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
          if( variableIsFilled("MTnuj_1stPair_PAS") && variableIsFilled("MTnuj_2ndPair_PAS") )
            STDOUT("PassMinMETPt1stEleThreshold: MTnuj_1stPair_PAS,MTnuj_2ndPair_PAS = "
                   <<getVariableValue("MTnuj_1stPair_PAS")
                   <<",\t"<<getVariableValue("MTnuj_2ndPair_PAS"));
          if( variableIsFilled("mDeltaPhiMETEle_PAS") && variableIsFilled("mDeltaPhiMET1stJet_PAS")
              && variableIsFilled("mDeltaPhiMET2ndJet_PAS") )
            STDOUT("PassMinMETPt1stEleThreshold: mDeltaPhiMETEle_PAS,mDeltaPhiMET1stJet_PAS,mDeltaPhiMET2ndJet_PAS = "
                   <<getVariableValue("mDeltaPhiMETEle_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET1stJet_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET2ndJet_PAS") );
          if( variableIsFilled("nMuon_PtCut_IDISO") )
            STDOUT("PassMinMETPt1stEleThreshold: nMuon_PtCut_IDISO = "<<getVariableValue("nMuon_PtCut_IDISO"));
          if( variableIsFilled("minDRej") )
            STDOUT("PassMinMETPt1stEleThreshold: minDRej = "<<getVariableValue("minDRej"));


          STDOUT("PassMinMETPt1stEleThreshold: ------------ END -------------");
        }
      }

    // Events in the tails of the Pt1stEle_PAS distribution that pass event pre-selection
    if( isData && variableIsFilled("Pt1stEle_PAS")
        && passedAllPreviousCuts("sT_presel") && passedCut("sT_presel")
      )
      {
        if( getVariableValue("Pt1stEle_PAS")>Pt1stEle_PAS_Thresh ) {
          STDOUT("PassPt1stEleThreshold: ----------- START ------------");

          STDOUT("PassPt1stEleThreshold: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
          if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") )
            STDOUT("PassPt1stEleThreshold: Pt1stEle_PAS,Eta1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS"));
          if( variableIsFilled("MET_PAS") )
            STDOUT("PassPt1stEleThreshold: MET_PAS = "<<getVariableValue("MET_PAS"));
          if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") )
            STDOUT("PassPt1stEleThreshold: Pt1stJet_PAS,Eta1stJet_PAS = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS"));
          if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") )
            STDOUT("PassPt1stEleThreshold: Pt2ndJet_PAS,Eta2ndJet_PAS = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS"));
          if( variableIsFilled("MTenu_PAS") && variableIsFilled("Mjj_PAS") )
            STDOUT("PassPt1stEleThreshold: MTenu_PAS,Mjj_PAS = "<<getVariableValue("MTenu_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
          if( variableIsFilled("sT_PAS") )
            STDOUT("PassPt1stEleThreshold: sT_PAS = "<<getVariableValue("sT_PAS"));
          if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )
            STDOUT("PassPt1stEleThreshold: Mej_1stPair_PAS,Mej_2ndPair_PAS = "
                   <<getVariableValue("Mej_1stPair_PAS")
                   <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
          if( variableIsFilled("MTnuj_1stPair_PAS") && variableIsFilled("MTnuj_2ndPair_PAS") )
            STDOUT("PassPt1stEleThreshold: MTnuj_1stPair_PAS,MTnuj_2ndPair_PAS = "
                   <<getVariableValue("MTnuj_1stPair_PAS")
                   <<",\t"<<getVariableValue("MTnuj_2ndPair_PAS"));
          if( variableIsFilled("mDeltaPhiMETEle_PAS") && variableIsFilled("mDeltaPhiMET1stJet_PAS")
              && variableIsFilled("mDeltaPhiMET2ndJet_PAS") )
            STDOUT("PassPt1stEleThreshold: mDeltaPhiMETEle_PAS,mDeltaPhiMET1stJet_PAS,mDeltaPhiMET2ndJet_PAS = "
                   <<getVariableValue("mDeltaPhiMETEle_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET1stJet_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET2ndJet_PAS") );
          if( variableIsFilled("nMuon_PtCut_IDISO") )
            STDOUT("PassPt1stEleThreshold: nMuon_PtCut_IDISO = "<<getVariableValue("nMuon_PtCut_IDISO"));
          if( variableIsFilled("minDRej") )
            STDOUT("PassPt1stEleThreshold: minDRej = "<<getVariableValue("minDRej"));


          STDOUT("PassPt1stEleThreshold: ------------ END -------------");
        }
      }


    // Events in the tails of the Pt1stJet_PAS distribution that pass event pre-selection
    if( isData && variableIsFilled("Pt1stJet_PAS")
        && passedAllPreviousCuts("sT_presel") && passedCut("sT_presel")
      )
      {
        if( getVariableValue("Pt1stJet_PAS")>Pt1stJet_PAS_Thresh ) {
          STDOUT("PassPt1stJetThreshold: ----------- START ------------");

          STDOUT("PassPt1stJetThreshold: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
          if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") )
            STDOUT("PassPt1stJetThreshold: Pt1stEle_PAS,Eta1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS"));
          if( variableIsFilled("MET_PAS") )
            STDOUT("PassPt1stJetThreshold: MET_PAS = "<<getVariableValue("MET_PAS"));
          if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") )
            STDOUT("PassPt1stJetThreshold: Pt1stJet_PAS,Eta1stJet_PAS = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS"));
          if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") )
            STDOUT("PassPt1stJetThreshold: Pt2ndJet_PAS,Eta2ndJet_PAS = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS"));
          if( variableIsFilled("MTenu_PAS") && variableIsFilled("Mjj_PAS") )
            STDOUT("PassPt1stJetThreshold: MTenu_PAS,Mjj_PAS = "<<getVariableValue("MTenu_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
          if( variableIsFilled("sT_PAS") )
            STDOUT("PassPt1stJetThreshold: sT_PAS = "<<getVariableValue("sT_PAS"));
          if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )
            STDOUT("PassPt1stJetThreshold: Mej_1stPair_PAS,Mej_2ndPair_PAS = "
                   <<getVariableValue("Mej_1stPair_PAS")
                   <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
          if( variableIsFilled("MTnuj_1stPair_PAS") && variableIsFilled("MTnuj_2ndPair_PAS") )
            STDOUT("PassPt1stJetThreshold: MTnuj_1stPair_PAS,MTnuj_2ndPair_PAS = "
                   <<getVariableValue("MTnuj_1stPair_PAS")
                   <<",\t"<<getVariableValue("MTnuj_2ndPair_PAS"));
          if( variableIsFilled("mDeltaPhiMETEle_PAS") && variableIsFilled("mDeltaPhiMET1stJet_PAS")
              && variableIsFilled("mDeltaPhiMET2ndJet_PAS") )
            STDOUT("PassPt1stJetThreshold: mDeltaPhiMETEle_PAS,mDeltaPhiMET1stJet_PAS,mDeltaPhiMET2ndJet_PAS = "
                   <<getVariableValue("mDeltaPhiMETEle_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET1stJet_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET2ndJet_PAS") );
          if( variableIsFilled("nMuon_PtCut_IDISO") )
            STDOUT("PassPt1stJetThreshold: nMuon_PtCut_IDISO = "<<getVariableValue("nMuon_PtCut_IDISO"));
          if( variableIsFilled("minDRej") )
            STDOUT("PassPt1stJetThreshold: minDRej = "<<getVariableValue("minDRej"));

          STDOUT("PassPt1stJetThreshold: ------------ END -------------");
        }
      }


    // Events in the tails of the Pt2ndJet_PAS distribution that pass event pre-selection
    if( isData && variableIsFilled("Pt2ndJet_PAS")
        && passedAllPreviousCuts("sT_presel") && passedCut("sT_presel")
      )
      {
        if( getVariableValue("Pt2ndJet_PAS")>Pt2ndJet_PAS_Thresh ) {
          STDOUT("PassPt2ndJetThreshold: ----------- START ------------");

          STDOUT("PassPt2ndJetThreshold: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
          if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") )
            STDOUT("PassPt2ndJetThreshold: Pt1stEle_PAS,Eta1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS"));
          if( variableIsFilled("MET_PAS") )
            STDOUT("PassPt2ndJetThreshold: MET_PAS = "<<getVariableValue("MET_PAS"));
          if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") )
            STDOUT("PassPt2ndJetThreshold: Pt1stJet_PAS,Eta1stJet_PAS = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS"));
          if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") )
            STDOUT("PassPt2ndJetThreshold: Pt2ndJet_PAS,Eta2ndJet_PAS = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS"));
          if( variableIsFilled("MTenu_PAS") && variableIsFilled("Mjj_PAS") )
            STDOUT("PassPt2ndJetThreshold: MTenu_PAS,Mjj_PAS = "<<getVariableValue("MTenu_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
          if( variableIsFilled("sT_PAS") )
            STDOUT("PassPt2ndJetThreshold: sT_PAS = "<<getVariableValue("sT_PAS"));
          if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )
            STDOUT("PassPt2ndJetThreshold: Mej_1stPair_PAS,Mej_2ndPair_PAS = "
                   <<getVariableValue("Mej_1stPair_PAS")
                   <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
          if( variableIsFilled("MTnuj_1stPair_PAS") && variableIsFilled("MTnuj_2ndPair_PAS") )
            STDOUT("PassPt2ndJetThreshold: MTnuj_1stPair_PAS,MTnuj_2ndPair_PAS = "
                   <<getVariableValue("MTnuj_1stPair_PAS")
                   <<",\t"<<getVariableValue("MTnuj_2ndPair_PAS"));
          if( variableIsFilled("mDeltaPhiMETEle_PAS") && variableIsFilled("mDeltaPhiMET1stJet_PAS")
              && variableIsFilled("mDeltaPhiMET2ndJet_PAS") )
            STDOUT("PassPt2ndJetThreshold: mDeltaPhiMETEle_PAS,mDeltaPhiMET1stJet_PAS,mDeltaPhiMET2ndJet_PAS = "
                   <<getVariableValue("mDeltaPhiMETEle_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET1stJet_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET2ndJet_PAS") );
          if( variableIsFilled("nMuon_PtCut_IDISO") )
            STDOUT("PassPt2ndJetThreshold: nMuon_PtCut_IDISO = "<<getVariableValue("nMuon_PtCut_IDISO"));
          if( variableIsFilled("minDRej") )
            STDOUT("PassPt2ndJetThreshold: minDRej = "<<getVariableValue("minDRej"));

          STDOUT("PassPt2ndJetThreshold: ------------ END -------------");
        }
      }


    // Events in the tails of the MTenu distribution that pass event pre-selection
    if( isData && variableIsFilled("MTenu")
        && passedAllPreviousCuts("sT_presel") && passedCut("sT_presel")
      )
      {
        if( getVariableValue("MTenu")>MTenu_Thresh ) {
          STDOUT("PassMTenuThreshold: ----------- START ------------");

          STDOUT("PassMTenuThreshold: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
          if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") )
            STDOUT("PassMTenuThreshold: Pt1stEle_PAS,Eta1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS"));
          if( variableIsFilled("MET_PAS") )
            STDOUT("PassMTenuThreshold: MET_PAS = "<<getVariableValue("MET_PAS"));
          if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") )
            STDOUT("PassMTenuThreshold: Pt1stJet_PAS,Eta1stJet_PAS = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS"));
          if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") )
            STDOUT("PassMTenuThreshold: Pt2ndJet_PAS,Eta2ndJet_PAS = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS"));
          if( variableIsFilled("MTenu_PAS") && variableIsFilled("Mjj_PAS") )
            STDOUT("PassMTenuThreshold: MTenu_PAS,Mjj_PAS = "<<getVariableValue("MTenu_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
          if( variableIsFilled("sT_PAS") )
            STDOUT("PassMTenuThreshold: sT_PAS = "<<getVariableValue("sT_PAS"));
          if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )
            STDOUT("PassMTenuThreshold: Mej_1stPair_PAS,Mej_2ndPair_PAS = "
                   <<getVariableValue("Mej_1stPair_PAS")
                   <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
          if( variableIsFilled("MTnuj_1stPair_PAS") && variableIsFilled("MTnuj_2ndPair_PAS") )
            STDOUT("PassMTenuThreshold: MTnuj_1stPair_PAS,MTnuj_2ndPair_PAS = "
                   <<getVariableValue("MTnuj_1stPair_PAS")
                   <<",\t"<<getVariableValue("MTnuj_2ndPair_PAS"));
          if( variableIsFilled("mDeltaPhiMETEle_PAS") && variableIsFilled("mDeltaPhiMET1stJet_PAS")
              && variableIsFilled("mDeltaPhiMET2ndJet_PAS") )
            STDOUT("PassMTenuThreshold: mDeltaPhiMETEle_PAS,mDeltaPhiMET1stJet_PAS,mDeltaPhiMET2ndJet_PAS = "
                   <<getVariableValue("mDeltaPhiMETEle_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET1stJet_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET2ndJet_PAS") );
          if( variableIsFilled("nMuon_PtCut_IDISO") )
            STDOUT("PassMTenuThreshold: nMuon_PtCut_IDISO = "<<getVariableValue("nMuon_PtCut_IDISO"));
          if( variableIsFilled("minDRej") )
            STDOUT("PassMTenuThreshold: minDRej = "<<getVariableValue("minDRej"));


          STDOUT("PassMTenuThreshold: ------------ END -------------");
        }
      }

    // Events in the tails of the sT distribution that pass event pre-selection
    if( isData && variableIsFilled("sT")
        && passedAllPreviousCuts("sT_presel") && passedCut("sT_presel")
      )
      {
        if( getVariableValue("sT")>sT_Thresh ) {
          STDOUT("PasssTThreshold: ----------- START ------------");

          STDOUT("PasssTThreshold: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
          if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") )
            STDOUT("PasssTThreshold: Pt1stEle_PAS,Eta1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS"));
          if( variableIsFilled("MET_PAS") )
            STDOUT("PasssTThreshold: MET_PAS = "<<getVariableValue("MET_PAS"));
          if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") )
            STDOUT("PasssTThreshold: Pt1stJet_PAS,Eta1stJet_PAS = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS"));
          if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") )
            STDOUT("PasssTThreshold: Pt2ndJet_PAS,Eta2ndJet_PAS = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS"));
          if( variableIsFilled("MTenu_PAS") && variableIsFilled("Mjj_PAS") )
            STDOUT("PasssTThreshold: MTenu_PAS,Mjj_PAS = "<<getVariableValue("MTenu_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
          if( variableIsFilled("sT_PAS") )
            STDOUT("PasssTThreshold: sT_PAS = "<<getVariableValue("sT_PAS"));
          if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )
            STDOUT("PasssTThreshold: Mej_1stPair_PAS,Mej_2ndPair_PAS = "
                   <<getVariableValue("Mej_1stPair_PAS")
                   <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
          if( variableIsFilled("MTnuj_1stPair_PAS") && variableIsFilled("MTnuj_2ndPair_PAS") )
            STDOUT("PasssTThreshold: MTnuj_1stPair_PAS,MTnuj_2ndPair_PAS = "
                   <<getVariableValue("MTnuj_1stPair_PAS")
                   <<",\t"<<getVariableValue("MTnuj_2ndPair_PAS"));
          if( variableIsFilled("mDeltaPhiMETEle_PAS") && variableIsFilled("mDeltaPhiMET1stJet_PAS")
              && variableIsFilled("mDeltaPhiMET2ndJet_PAS") )
            STDOUT("PasssTThreshold: mDeltaPhiMETEle_PAS,mDeltaPhiMET1stJet_PAS,mDeltaPhiMET2ndJet_PAS = "
                   <<getVariableValue("mDeltaPhiMETEle_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET1stJet_PAS")
                   <<",\t"<<getVariableValue("mDeltaPhiMET2ndJet_PAS") );
          if( variableIsFilled("nMuon_PtCut_IDISO") )
            STDOUT("PasssTThreshold: nMuon_PtCut_IDISO = "<<getVariableValue("nMuon_PtCut_IDISO"));
          if( variableIsFilled("minDRej") )
            STDOUT("PasssTThreshold: minDRej = "<<getVariableValue("minDRej"));

          STDOUT("PasssTThreshold: ------------ END -------------");
        }
      }

    //W + >=0 jets
    if(doPlot_Wmore0jet)
      {
	if( passedAllPreviousCuts("nJet_all")
	    && variableIsFilled("MTenu_PAS") && variableIsFilled("Pt1stEle_PAS")
	    && variableIsFilled("Eta1stEle_PAS") && variableIsFilled("Phi1stEle_PAS")
	    && variableIsFilled("MET_PAS") && variableIsFilled("mDeltaPhiMETEle")
	    && variableIsFilled("METPhi_PAS")
	    )
	  {
	    FillUserTH1D("h1_Pt1stEle_W0jet", getVariableValue("Pt1stEle_PAS") );
	    FillUserTH1D("h1_Eta1stEle_W0jet", getVariableValue("Eta1stEle_PAS") );
	    FillUserTH1D("h1_Phi1stEle_W0jet", getVariableValue("Phi1stEle_PAS") );
	    FillUserTH1D("h1_MET_W0jet", getVariableValue("MET_PAS") );
	    FillUserTH1D("h1_METPhi_W0jet", getVariableValue("METPhi_PAS") );
	    FillUserTH1D("h1_MTenu_W0jet", getVariableValue("MTenu_PAS") );
	    FillUserTH1D("h1_mDeltaPhiMETEle_W0jet", getVariableValue("mDeltaPhiMETEle") );

	    if( fabs(getVariableValue("Eta1stEle_PAS")) <= eleEta_bar )
	      {
		FillUserTH1D("h1_Pt1stEle_W0jet_barrel", getVariableValue("Pt1stEle_PAS") );
		FillUserTH1D("h1_Energy1stEle_W0jet_barrel", ElectronCaloEnergy->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
		FillUserTH1D("h1_TrackPt1stEle_W0jet_barrel", ElectronTrackPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	      }
	    else
	      {
		FillUserTH1D("h1_Pt1stEle_W0jet_endcap", getVariableValue("Pt1stEle_PAS") );
		FillUserTH1D("h1_Energy1stEle_W0jet_endcap", ElectronCaloEnergy->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
		FillUserTH1D("h1_TrackPt1stEle_W0jet_endcap", ElectronTrackPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	      }

	    //ET/PT vs ET
	    double EToverPT = ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) / ElectronTrackPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) ;
	    FillUserTH2D("h2_EToverPT_vs_ET_1stEle_W0jet", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]), EToverPT );

	    if( isData && getVariableValue("Pt1stEle_PAS")>150 && getVariableValue("Pt1stEle_PAS")<175 )
	      {
		STDOUT("PassPtEleW0jetThreshold: ----------- START ------------");

		STDOUT("PassPtEleW0jetThreshold: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
		if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") )
		  STDOUT("PassPtEleW0jetThreshold: Pt1stEle_PAS,Eta1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS"));
		if( variableIsFilled("MET_PAS") )
		  STDOUT("PassPtEleW0jetThreshold: MET_PAS = "<<getVariableValue("MET_PAS"));
		if( variableIsFilled("MTenu_PAS") )
		  STDOUT("PassPtEleW0jetThreshold: MTenu_PAS = "<<getVariableValue("MTenu_PAS"));
		if( variableIsFilled("mDeltaPhiMETEle_PAS") )
		  STDOUT("PassPtEleW0jetThreshold: mDeltaPhiMETEle_PAS = "
			 <<getVariableValue("mDeltaPhiMETEle_PAS"));

		STDOUT("PassPtEleW0jetThreshold: ------------ END -------------");
	      }

	  }
      }

    //W + >=1 jets
    if(doPlot_Wmore1jet)
      {
	if( passedAllPreviousCuts("nJet_all") && passedCut("Pt1stJet_noOvrlp_ID") && passedCut("mEta1stJet_noOvrlp_ID")
	    && variableIsFilled("MTenu_PAS") && variableIsFilled("Pt1stEle_PAS")
	    && variableIsFilled("Eta1stEle_PAS") && variableIsFilled("Phi1stEle_PAS")
	    && variableIsFilled("MET_PAS") && variableIsFilled("mDeltaPhiMETEle")
	    && variableIsFilled("METPhi_PAS")
	    )
	  {
	    FillUserTH1D("h1_Pt1stEle_W1jet", getVariableValue("Pt1stEle_PAS") );
	    FillUserTH1D("h1_Eta1stEle_W1jet", getVariableValue("Eta1stEle_PAS") );
	    FillUserTH1D("h1_Phi1stEle_W1jet", getVariableValue("Phi1stEle_PAS") );
	    FillUserTH1D("h1_MET_W1jet", getVariableValue("MET_PAS") );
	    FillUserTH1D("h1_METPhi_W1jet", getVariableValue("METPhi_PAS") );
	    FillUserTH1D("h1_MTenu_W1jet", getVariableValue("MTenu_PAS") );
	    FillUserTH1D("h1_mDeltaPhiMETEle_W1jet", getVariableValue("mDeltaPhiMETEle") );

	    if( fabs(getVariableValue("Eta1stEle_PAS")) <= eleEta_bar )
	      {
		FillUserTH1D("h1_Pt1stEle_W1jet_barrel", getVariableValue("Pt1stEle_PAS") );
		FillUserTH1D("h1_Energy1stEle_W1jet_barrel", ElectronCaloEnergy->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
		FillUserTH1D("h1_TrackPt1stEle_W1jet_barrel", ElectronTrackPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	      }
	    else
	      {
		FillUserTH1D("h1_Pt1stEle_W1jet_endcap", getVariableValue("Pt1stEle_PAS") );
		FillUserTH1D("h1_Energy1stEle_W1jet_endcap", ElectronCaloEnergy->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
		FillUserTH1D("h1_TrackPt1stEle_W1jet_endcap", ElectronTrackPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	      }

	    //ET/PT vs ET
	    double EToverPT = ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) / ElectronTrackPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) ;
	    FillUserTH2D("h2_EToverPT_vs_ET_1stEle_W1jet", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]), EToverPT );

	    if( isData && getVariableValue("Pt1stEle_PAS")>150 && getVariableValue("Pt1stEle_PAS")<175 )
	      {
		STDOUT("PassPtEleW1jetThreshold: ----------- START ------------");

		STDOUT("PassPtEleW1jetThreshold: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
		if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") )
		  STDOUT("PassPtEleW1jetThreshold: Pt1stEle_PAS,Eta1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS"));
		if( variableIsFilled("MET_PAS") )
		  STDOUT("PassPtEleW1jetThreshold: MET_PAS = "<<getVariableValue("MET_PAS"));
		if( variableIsFilled("MTenu_PAS") )
		  STDOUT("PassPtEleW1jetThreshold: MTenu_PAS = "<<getVariableValue("MTenu_PAS"));
		if( variableIsFilled("mDeltaPhiMETEle_PAS") )
		  STDOUT("PassPtEleW1jetThreshold: mDeltaPhiMETEle_PAS = "
			 <<getVariableValue("mDeltaPhiMETEle_PAS"));

		STDOUT("PassPtEleW1jetThreshold: ------------ END -------------");
	      }

	  }
      }

    // Events with Pt1stEle>200 GeV passing enujj pre-selection
    if( isData && (
        (run==144089 && ls==257 && event==310079631) ||    // only in Nov4 re-reco
        (run==146944 && ls==207 && event==130150965) ||    // only in Sep17+PromptReco
        (run==147754 && ls==216 && event==258606954) ||    // only in Nov4 re-reco
        (run==147926 && ls==199 && event==186563265) ||    // only in Nov4 re-reco
        (run==148862 && ls==403 && event==596914505) ||    // both in Sep17+PromptReco and Nov4 re-reco
        (run==149011 && ls==518 && event==741762875) ||    // only in Nov4 re-reco
        (run==149011 && ls==521 && event==745453692) ||    // both in Sep17+PromptReco and Nov4 re-reco
        (run==149181 && ls==1075 && event==1060203426) ||  // only in Nov4 re-reco
        (run==149181 && ls==1611 && event==1548701490) ||  // only in Nov4 re-reco
        (run==149181 && ls==1644 && event==1576395911) )   // both in Sep17+PromptReco and Nov4 re-reco
      )
      {
        STDOUT("UserPt1stEle: ----------- START ------------");

        STDOUT("UserPt1stEle: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
        if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") && variableIsFilled("Phi1stEle_PAS") )
//           STDOUT("UserPt1stEle: Pt1stEle_PAS,Eta1stEle_PAS,Phi1stEle,kOutOfTime = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS")<<",\t"<<getVariableValue("Phi1stEle_PAS")<<",\t"<<ElectronSCkOutOfTime->at(v_idx_ele_PtCut_IDISO_noOverlap[0]));
          STDOUT("UserPt1stEle: Pt1stEle_PAS,Eta1stEle_PAS,Phi1stEle = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS")<<",\t"<<getVariableValue("Phi1stEle_PAS"));
        if( variableIsFilled("MET_PAS") && variableIsFilled("METPhi_PAS") )
          STDOUT("UserPt1stEle: MET_PAS, METPhi_PAS = "<<getVariableValue("MET_PAS")<<",\t"<<getVariableValue("METPhi_PAS"));
        if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") && variableIsFilled("Phi1stJet_PAS") )
          STDOUT("UserPt1stEle: Pt1stJet_PAS,Eta1stJet_PAS,Phi1stJet = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS")<<",\t"<<getVariableValue("Phi1stJet_PAS"));
        if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") && variableIsFilled("Phi2ndJet_PAS") )
          STDOUT("UserPt1stEle: Pt2ndJet_PAS,Eta2ndJet_PAS,Phi2ndJet = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS")<<",\t"<<getVariableValue("Phi2ndJet_PAS"));
        if( variableIsFilled("MTenu_PAS") && variableIsFilled("Mjj_PAS") )
          STDOUT("UserPt1stEle: MTenu_PAS,Mjj_PAS = "<<getVariableValue("MTenu_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
        if( variableIsFilled("sT_PAS") )
          STDOUT("UserPt1stEle: sT_PAS = "<<getVariableValue("sT_PAS"));
        if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )
          STDOUT("UserPt1stEle: Mej_1stPair_PAS,Mej_2ndPair_PAS = "
                 <<getVariableValue("Mej_1stPair_PAS")
                 <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
        if( variableIsFilled("MTnuj_1stPair_PAS") && variableIsFilled("MTnuj_2ndPair_PAS") )
          STDOUT("UserPt1stEle: MTnuj_1stPair_PAS,MTnuj_2ndPair_PAS = "
                 <<getVariableValue("MTnuj_1stPair_PAS")
                 <<",\t"<<getVariableValue("MTnuj_2ndPair_PAS"));
        if( variableIsFilled("mDeltaPhiMETEle_PAS") && variableIsFilled("mDeltaPhiMET1stJet_PAS")
            && variableIsFilled("mDeltaPhiMET2ndJet_PAS") )
          STDOUT("UserPt1stEle: mDeltaPhiMETEle_PAS,mDeltaPhiMET1stJet_PAS,mDeltaPhiMET2ndJet_PAS = "
                 <<getVariableValue("mDeltaPhiMETEle_PAS")
                 <<",\t"<<getVariableValue("mDeltaPhiMET1stJet_PAS")
                 <<",\t"<<getVariableValue("mDeltaPhiMET2ndJet_PAS") );
        if( variableIsFilled("nMuon_PtCut_IDISO") )
          STDOUT("UserPt1stEle: nMuon_PtCut_IDISO = "<<getVariableValue("nMuon_PtCut_IDISO"));
        if( variableIsFilled("minDRej") )
          STDOUT("UserPt1stEle: minDRej = "<<getVariableValue("minDRej"));


        STDOUT("UserPt1stEle: ------------ END -------------");
      }

    // Events with MTenu>200 GeV passing enujj pre-selection
    if( isData && (
        (run==142422 && ls==177 && event==116725694) ||   // both in Sep17+PromptReco and Nov4 re-reco
        (run==147115 && ls==275 && event==327457974) ||   // both in Sep17+PromptReco and Nov4 re-reco
        (run==147115 && ls==307 && event==363025814) ||   // only in Nov4 re-reco
        (run==147926 && ls==199 && event==186563265) ||   // only in Nov4 re-reco
        (run==147927 && ls==93 && event==98222777) ||     // both in Sep17+PromptReco and Nov4 re-reco
        (run==148029 && ls==53 && event==9851751) ||      // both in Sep17+PromptReco and Nov4 re-reco
        (run==148029 && ls==233 && event==174059244) ||   // both in Sep17+PromptReco and Nov4 re-reco
        (run==148822 && ls==213 && event==227962384) ||   // both in Sep17+PromptReco and Nov4 re-reco
        (run==148860 && ls==28 && event==38009489) ||     // both in Sep17+PromptReco and Nov4 re-reco
        (run==148864 && ls==572 && event==646881145) ||   // both in Sep17+PromptReco and Nov4 re-reco
        (run==149011 && ls==518 && event==741762875) ||   // only in Nov4 re-reco
        (run==149011 && ls==521 && event==745453692) ||   // both in Sep17+PromptReco and Nov4 re-reco
        (run==149181 && ls==292 && event==138356895) ||   // both in Sep17+PromptReco and Nov4 re-reco
        (run==149181 && ls==1611 && event==1548701490) || // only in Nov4 re-reco
        (run==149291 && ls==613 && event==639394420) ||   // both in Sep17+PromptReco and Nov4 re-reco
        (run==149291 && ls==767 && event==772369563) )    // both in Sep17+PromptReco and Nov4 re-reco
      )
      {
        STDOUT("UserMTenu: ----------- START ------------");

        STDOUT("UserMTenu: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
        if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") && variableIsFilled("Phi1stEle_PAS") )
//           STDOUT("UserMTenu: Pt1stEle_PAS,Eta1stEle_PAS,Phi1stEle,kOutOfTime = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS")<<",\t"<<getVariableValue("Phi1stEle_PAS")<<",\t"<<ElectronSCkOutOfTime->at(v_idx_ele_PtCut_IDISO_noOverlap[0]));
          STDOUT("UserMTenu: Pt1stEle_PAS,Eta1stEle_PAS,Phi1stEle = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS")<<",\t"<<getVariableValue("Phi1stEle_PAS"));
        if( variableIsFilled("MET_PAS") && variableIsFilled("METPhi_PAS") )
          STDOUT("UserMTenu: MET_PAS, METPhi_PAS = "<<getVariableValue("MET_PAS")<<",\t"<<getVariableValue("METPhi_PAS"));
        if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") && variableIsFilled("Phi1stJet_PAS") )
          STDOUT("UserMTenu: Pt1stJet_PAS,Eta1stJet_PAS,Phi1stJet = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS")<<",\t"<<getVariableValue("Phi1stJet_PAS"));
        if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") && variableIsFilled("Phi2ndJet_PAS") )
          STDOUT("UserMTenu: Pt2ndJet_PAS,Eta2ndJet_PAS,Phi2ndJet = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS")<<",\t"<<getVariableValue("Phi2ndJet_PAS"));
        if( variableIsFilled("MTenu_PAS") && variableIsFilled("Mjj_PAS") )
          STDOUT("UserMTenu: MTenu_PAS,Mjj_PAS = "<<getVariableValue("MTenu_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
        if( variableIsFilled("sT_PAS") )
          STDOUT("UserMTenu: sT_PAS = "<<getVariableValue("sT_PAS"));
        if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )
          STDOUT("UserMTenu: Mej_1stPair_PAS,Mej_2ndPair_PAS = "
                 <<getVariableValue("Mej_1stPair_PAS")
                 <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
        if( variableIsFilled("MTnuj_1stPair_PAS") && variableIsFilled("MTnuj_2ndPair_PAS") )
          STDOUT("UserMTenu: MTnuj_1stPair_PAS,MTnuj_2ndPair_PAS = "
                 <<getVariableValue("MTnuj_1stPair_PAS")
                 <<",\t"<<getVariableValue("MTnuj_2ndPair_PAS"));
        if( variableIsFilled("mDeltaPhiMETEle_PAS") && variableIsFilled("mDeltaPhiMET1stJet_PAS")
            && variableIsFilled("mDeltaPhiMET2ndJet_PAS") )
          STDOUT("UserMTenu: mDeltaPhiMETEle_PAS,mDeltaPhiMET1stJet_PAS,mDeltaPhiMET2ndJet_PAS = "
                 <<getVariableValue("mDeltaPhiMETEle_PAS")
                 <<",\t"<<getVariableValue("mDeltaPhiMET1stJet_PAS")
                 <<",\t"<<getVariableValue("mDeltaPhiMET2ndJet_PAS") );
        if( variableIsFilled("nMuon_PtCut_IDISO") )
          STDOUT("UserMTenu: nMuon_PtCut_IDISO = "<<getVariableValue("nMuon_PtCut_IDISO"));
        if( variableIsFilled("minDRej") )
          STDOUT("UserMTenu: minDRej = "<<getVariableValue("minDRej"));


        STDOUT("UserMTenu: ------------ END -------------");
      }

    //EXTRA CHECKS
    if( doExtraChecks
        && passedAllPreviousCuts("sT_presel") && passedCut("sT_presel")
	)
      {

	//FULL selection
	if( passedAllPreviousCuts("sT_presel")
            && passedCut("sT_presel")
            && passedCut("nMuon_PtCut_IDISO")
            && passedCut("MTenu")
	    && passedCut("minMETPt1stEle")
	    )
	  {
	    FillUserTH1D("h1_Njet_fullSel", getVariableValue("nJet_PtCut_noOvrlp_ID") );
            FillUserTH1D("h1_NjetTCHELBTag_fullSel", getVariableValue("nJet_TCHELBTag") );
	    FillUserTH1D("h1_minDRej_fullSel", getVariableValue("minDRej") );
	    FillUserTH1D("h1_mDeltaPhiMETEle_fullSel", getVariableValue("mDeltaPhiMETEle") );
	    FillUserTH1D("h1_mDeltaPhiMET1stJet_fullSel", getVariableValue("mDeltaPhiMET1stJet") );
	    FillUserTH1D("h1_mDeltaPhiMET2ndJet_fullSel", getVariableValue("mDeltaPhiMET2ndJet") );
	    FillUserTH1D("h1_DeltaRjets_PAS_fullSel", getVariableValue("DeltaRjets_PAS") );
            FillUserTH1D("h1_Vtxd01stEle_PAS_fullSel", getVariableValue("Vtxd01stEle_PAS") );
            FillUserTH1D("h1_MissingHits1stEle_PAS_fullSel", getVariableValue("MissingHits1stEle_PAS") );
            FillUserTH1D("h1_Dist1stEle_PAS_fullSel", getVariableValue("Dist1stEle_PAS") );
            FillUserTH1D("h1_DCotTheta1stEle_PAS_fullSel", getVariableValue("DCotTheta1stEle_PAS") );
            FillUserTH1D("h1_Conversion1stEle_fullSel", (getVariableValue("MissingHits1stEle_PAS")>=1 && getVariableValue("DCotTheta1stEle_PAS")<0.02 && getVariableValue("DCotTheta1stEle_PAS")<0.02) ? 1 : 0 );
	  }

	//High MT after pre-selection
        if ( variableIsFilled("MTenu") ) {
          if( getVariableValue("MTenu")>MTenu_Thresh )
            {
              //---------------------------------------------------------------
              //1D distributions
              FillUserTH1D("h1_Njet_highMT", getVariableValue("nJet_PtCut_noOvrlp_ID") );
              FillUserTH1D("h1_NjetTCHELBTag_highMT", getVariableValue("nJet_TCHELBTag") );
              FillUserTH1D("h1_minDRej_highMT", getVariableValue("minDRej") );
              FillUserTH1D("h1_mDeltaPhiMETEle_highMT", getVariableValue("mDeltaPhiMETEle") );
              FillUserTH1D("h1_mDeltaPhiMET1stJet_highMT", getVariableValue("mDeltaPhiMET1stJet") );
              FillUserTH1D("h1_mDeltaPhiMET2ndJet_highMT", getVariableValue("mDeltaPhiMET2ndJet") );
              FillUserTH1D("h1_DeltaRjets_PAS_highMT", getVariableValue("DeltaRjets_PAS") );
              FillUserTH1D("h1_Vtxd01stEle_PAS_highMT", getVariableValue("Vtxd01stEle_PAS") );
              FillUserTH1D("h1_MissingHits1stEle_PAS_highMT", getVariableValue("MissingHits1stEle_PAS") );
              FillUserTH1D("h1_Dist1stEle_PAS_highMT", getVariableValue("Dist1stEle_PAS") );
              FillUserTH1D("h1_DCotTheta1stEle_PAS_highMT", getVariableValue("DCotTheta1stEle_PAS") );
              FillUserTH1D("h1_Conversion1stEle_highMT", (getVariableValue("MissingHits1stEle_PAS")>=1 && getVariableValue("DCotTheta1stEle_PAS")<0.02 && getVariableValue("DCotTheta1stEle_PAS")<0.02) ? 1 : 0 );
              //---------------------------------------------------------------

              //---------------------------------------------------------------
              //2D distributions
              double minDeltaPhi = min( getVariableValue("mDeltaPhiMET1stJet_PAS"), getVariableValue("mDeltaPhiMET2ndJet_PAS") );
              FillUserTH2D("h2_minDeltaPhiMETj_vs_minDeltaRej_highMT", getVariableValue("minDRej"), minDeltaPhi );
              if(doGenLevelStudies && isData==0)
                {
                  FillUserTH2D("h2_Ngenlept_vs_Njets_highMT" , getVariableValue("nJet_PtCut_noOvrlp_ID") , Ngenleptons);

                  if( Ngenelectrons==2 )
                    h1_Ngenlept__ee_emu_etau_eonly_noe_others->AddBinContent(1);
                  else if( Ngenelectrons==1 && Ngenmuons==1 )
                    h1_Ngenlept__ee_emu_etau_eonly_noe_others->AddBinContent(2);
                  else if( Ngenelectrons==1 && Ngentaus==1 )
                    h1_Ngenlept__ee_emu_etau_eonly_noe_others->AddBinContent(3);
                  else if( Ngenelectrons==1 )
                    h1_Ngenlept__ee_emu_etau_eonly_noe_others->AddBinContent(4);
                  else if( Ngenelectrons==0 )
                    h1_Ngenlept__ee_emu_etau_eonly_noe_others->AddBinContent(5);
                  else
                    h1_Ngenlept__ee_emu_etau_eonly_noe_others->AddBinContent(6);
                  // 		cout << "leptons, electrons, muons, taus: "
                  // 		     << Ngenleptons << ", " << Ngenelectrons << ", " << Ngenmuons << ", " << Ngentaus << endl;
                  //---------------------------------------------------------------
                }
            }
        }
        //High Pt1stEle after pre-selection
        if ( variableIsFilled("Pt1stEle_PAS") ) {
          if( getVariableValue("Pt1stEle_PAS")>Pt1stEle_PAS_Thresh )
            {
              //---------------------------------------------------------------
              //1D distributions
              FillUserTH1D("h1_Njet_highPt1stEle", getVariableValue("nJet_PtCut_noOvrlp_ID") );
              FillUserTH1D("h1_NjetTCHELBTag_highPt1stEle", getVariableValue("nJet_TCHELBTag") );
              FillUserTH1D("h1_minDRej_highPt1stEle", getVariableValue("minDRej") );
              FillUserTH1D("h1_mDeltaPhiMETEle_highPt1stEle", getVariableValue("mDeltaPhiMETEle") );
              FillUserTH1D("h1_mDeltaPhiMET1stJet_highPt1stEle", getVariableValue("mDeltaPhiMET1stJet") );
              FillUserTH1D("h1_mDeltaPhiMET2ndJet_highPt1stEle", getVariableValue("mDeltaPhiMET2ndJet") );
              FillUserTH1D("h1_DeltaRjets_PAS_highPt1stEle", getVariableValue("DeltaRjets_PAS") );
              FillUserTH1D("h1_Vtxd01stEle_PAS_highPt1stEle", getVariableValue("Vtxd01stEle_PAS") );
              FillUserTH1D("h1_MissingHits1stEle_PAS_highPt1stEle", getVariableValue("MissingHits1stEle_PAS") );
              FillUserTH1D("h1_Dist1stEle_PAS_highPt1stEle", getVariableValue("Dist1stEle_PAS") );
              FillUserTH1D("h1_DCotTheta1stEle_PAS_highPt1stEle", getVariableValue("DCotTheta1stEle_PAS") );
              FillUserTH1D("h1_Conversion1stEle_highPt1stEle", (getVariableValue("MissingHits1stEle_PAS")>=1 && getVariableValue("DCotTheta1stEle_PAS")<0.02 && getVariableValue("DCotTheta1stEle_PAS")<0.02) ? 1 : 0 );
              //---------------------------------------------------------------
            }
        }
        //High Pt1stJet after pre-selection
        if ( variableIsFilled("Pt1stJet_PAS") ) {
          if( getVariableValue("Pt1stJet_PAS")>Pt1stJet_PAS_Thresh )
            {
              //---------------------------------------------------------------
              //1D distributions
              FillUserTH1D("h1_Njet_highPt1stJet", getVariableValue("nJet_PtCut_noOvrlp_ID") );
              FillUserTH1D("h1_NjetTCHELBTag_highPt1stJet", getVariableValue("nJet_TCHELBTag") );
              FillUserTH1D("h1_minDRej_highPt1stJet", getVariableValue("minDRej") );
              FillUserTH1D("h1_mDeltaPhiMETEle_highPt1stJet", getVariableValue("mDeltaPhiMETEle") );
              FillUserTH1D("h1_mDeltaPhiMET1stJet_highPt1stJet", getVariableValue("mDeltaPhiMET1stJet") );
              FillUserTH1D("h1_mDeltaPhiMET2ndJet_highPt1stJet", getVariableValue("mDeltaPhiMET2ndJet") );
              FillUserTH1D("h1_DeltaRjets_PAS_highPt1stJet", getVariableValue("DeltaRjets_PAS") );
              FillUserTH1D("h1_Vtxd01stEle_PAS_highPt1stJet", getVariableValue("Vtxd01stEle_PAS") );
              FillUserTH1D("h1_MissingHits1stEle_PAS_highPt1stJet", getVariableValue("MissingHits1stEle_PAS") );
              FillUserTH1D("h1_Dist1stEle_PAS_highPt1stJet", getVariableValue("Dist1stEle_PAS") );
              FillUserTH1D("h1_DCotTheta1stEle_PAS_highPt1stJet", getVariableValue("DCotTheta1stEle_PAS") );
              FillUserTH1D("h1_Conversion1stEle_highPt1stJet", (getVariableValue("MissingHits1stEle_PAS")>=1 && getVariableValue("DCotTheta1stEle_PAS")<0.02 && getVariableValue("DCotTheta1stEle_PAS")<0.02) ? 1 : 0 );
              //---------------------------------------------------------------
            }
        }
        //High Pt2ndJet after pre-selection
        if ( variableIsFilled("Pt2ndJet_PAS") ) {
          if( getVariableValue("Pt2ndJet_PAS")>Pt2ndJet_PAS_Thresh )
            {
              //---------------------------------------------------------------
              //1D distributions
              FillUserTH1D("h1_Njet_highPt2ndJet", getVariableValue("nJet_PtCut_noOvrlp_ID") );
              FillUserTH1D("h1_NjetTCHELBTag_highPt2ndJet", getVariableValue("nJet_TCHELBTag") );
              FillUserTH1D("h1_minDRej_highPt2ndJet", getVariableValue("minDRej") );
              FillUserTH1D("h1_mDeltaPhiMETEle_highPt2ndJet", getVariableValue("mDeltaPhiMETEle") );
              FillUserTH1D("h1_mDeltaPhiMET1stJet_highPt2ndJet", getVariableValue("mDeltaPhiMET1stJet") );
              FillUserTH1D("h1_mDeltaPhiMET2ndJet_highPt2ndJet", getVariableValue("mDeltaPhiMET2ndJet") );
              FillUserTH1D("h1_DeltaRjets_PAS_highPt2ndJet", getVariableValue("DeltaRjets_PAS") );
              FillUserTH1D("h1_Vtxd01stEle_PAS_highPt2ndJet", getVariableValue("Vtxd01stEle_PAS") );
              FillUserTH1D("h1_MissingHits1stEle_PAS_highPt2ndJet", getVariableValue("MissingHits1stEle_PAS") );
              FillUserTH1D("h1_Dist1stEle_PAS_highPt2ndJet", getVariableValue("Dist1stEle_PAS") );
              FillUserTH1D("h1_DCotTheta1stEle_PAS_highPt2ndJet", getVariableValue("DCotTheta1stEle_PAS") );
              FillUserTH1D("h1_Conversion1stEle_highPt2ndJet", (getVariableValue("MissingHits1stEle_PAS")>=1 && getVariableValue("DCotTheta1stEle_PAS")<0.02 && getVariableValue("DCotTheta1stEle_PAS")<0.02) ? 1 : 0 );
              //---------------------------------------------------------------
            }
        }
        //High Mej after pre-selection
        if ( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") ) {
          if( max(getVariableValue("Mej_1stPair_PAS"),getVariableValue("Mej_2ndPair_PAS"))>Mej_Thresh )
            {
              //---------------------------------------------------------------
              //1D distributions
              FillUserTH1D("h1_Njet_highMej", getVariableValue("nJet_PtCut_noOvrlp_ID") );
              FillUserTH1D("h1_NjetTCHELBTag_highMej", getVariableValue("nJet_TCHELBTag") );
              FillUserTH1D("h1_minDRej_highMej", getVariableValue("minDRej") );
              FillUserTH1D("h1_mDeltaPhiMETEle_highMej", getVariableValue("mDeltaPhiMETEle") );
              FillUserTH1D("h1_mDeltaPhiMET1stJet_highMej", getVariableValue("mDeltaPhiMET1stJet") );
              FillUserTH1D("h1_mDeltaPhiMET2ndJet_highMej", getVariableValue("mDeltaPhiMET2ndJet") );
              FillUserTH1D("h1_DeltaRjets_PAS_highMej", getVariableValue("DeltaRjets_PAS") );
              FillUserTH1D("h1_Vtxd01stEle_PAS_highMej", getVariableValue("Vtxd01stEle_PAS") );
              FillUserTH1D("h1_MissingHits1stEle_PAS_highMej", getVariableValue("MissingHits1stEle_PAS") );
              FillUserTH1D("h1_Dist1stEle_PAS_highMej", getVariableValue("Dist1stEle_PAS") );
              FillUserTH1D("h1_DCotTheta1stEle_PAS_highMej", getVariableValue("DCotTheta1stEle_PAS") );
              FillUserTH1D("h1_Conversion1stEle_highMej", (getVariableValue("MissingHits1stEle_PAS")>=1 && getVariableValue("DCotTheta1stEle_PAS")<0.02 && getVariableValue("DCotTheta1stEle_PAS")<0.02) ? 1 : 0 );
              //---------------------------------------------------------------
              // mDeltaPhiMET1stJet>2.5 && 20 <= PFJetNConstituents <= 30
              if( getVariableValue("mDeltaPhiMET1stJet")>2.5
		  // 		  && PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0])>=20
		  // 		  && PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0])<=30
		  ) {
                //---------------------------------------------------------------
                //1D distributions
                FillUserTH1D("h1_Pt1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Pt1stJet_PAS") );
                FillUserTH1D("h1_Eta1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Eta1stJet_PAS") );
                FillUserTH1D("h1_Phi1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Phi1stJet_PAS") );
                FillUserTH1D("h1_CHF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                FillUserTH1D("h1_NHF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", PFJetNeutralHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                FillUserTH1D("h1_CEF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                FillUserTH1D("h1_NEF1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", PFJetNeutralEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                FillUserTH1D("h1_NCH1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                FillUserTH1D("h1_NN1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0]) - PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                FillUserTH1D("h1_NC1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                FillUserTH1D("h1_Mass1stJet_highMej_mDeltaPhiMET1stJet_gt_2.5", Mass1stJet);
                FillUserTH1D("h1_Pt2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Pt2ndJet_PAS") );
                FillUserTH1D("h1_Eta2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Eta2ndJet_PAS") );
                FillUserTH1D("h1_Phi2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Phi2ndJet_PAS") );
                FillUserTH1D("h1_CHF2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                FillUserTH1D("h1_NHF2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", PFJetNeutralHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                FillUserTH1D("h1_CEF2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                FillUserTH1D("h1_NEF2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", PFJetNeutralEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                FillUserTH1D("h1_NCH2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                FillUserTH1D("h1_NN2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[1]) - PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                FillUserTH1D("h1_NC2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                FillUserTH1D("h1_Mass2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", Mass2ndJet);
                FillUserTH1D("h1_E1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", ElectronEnergy->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
                FillUserTH1D("h1_Pt1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Pt1stEle_PAS") );
                FillUserTH1D("h1_Eta1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Eta1stEle_PAS") );
                FillUserTH1D("h1_Phi1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Phi1stEle_PAS") );
                FillUserTH1D("h1_Charge1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Charge1stEle_PAS") );
                FillUserTH1D("h1_mDeltaPhiMETEle_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("mDeltaPhiMETEle") );
                FillUserTH1D("h1_mDeltaPhiEle1stJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("mDeltaPhiEle1stJet_PAS") );
                FillUserTH1D("h1_mDeltaPhiEle2ndJet_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("mDeltaPhiEle2ndJet_PAS") );
                FillUserTH1D("h1_mDeltaPhiMET2ndJet_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("mDeltaPhiMET2ndJet") );
                FillUserTH1D("h1_Ptenu_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Ptenu_PAS") );
                FillUserTH1D("h1_1stJet_PTOverPTPlusMET_highMej_mDeltaPhiMET1stJet_gt_2.5", (getVariableValue("Pt1stJet_PAS") / (getVariableValue("Pt1stJet_PAS") + getVariableValue("MET"))) );
                FillUserTH1D("h1_Conversion1stEle_highMej_mDeltaPhiMET1stJet_gt_2.5", (getVariableValue("MissingHits1stEle_PAS")>=1 && getVariableValue("DCotTheta1stEle_PAS")<0.02 && getVariableValue("DCotTheta1stEle_PAS")<0.02) ? 1 : 0 );
                FillUserTH1D("h1_MET_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("MET_PAS") );
                FillUserTH1D("h1_CaloMET_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", CaloMET->at(0) );
                FillUserTH1D("h1_TCMET_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", TCMET->at(0) );
                FillUserTH1D("h1_MET_over_CaloMET_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("MET_PAS") / CaloMET->at(0) );
                FillUserTH1D("h1_Njet_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("nJet_PtCut_noOvrlp_ID") );
                FillUserTH1D("h1_NjetTCHELBTag_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("nJet_TCHELBTag") );
                FillUserTH1D("h1_minDRej_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("minDRej") );
                FillUserTH1D("h1_maxDRej_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("maxDRej") );
                FillUserTH1D("h1_DRjets_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("DeltaRjets_PAS") );
		FillUserTH1D("h1_MTenu_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("MTenu_PAS") );
                //---------------------------------------------------------------
                //2D distributions
                FillUserTH2D("h2_EtaPhi1stEle_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Eta1stEle_PAS"), getVariableValue("Phi1stEle_PAS") );
                FillUserTH2D("h2_MinMaxMej_PAS_highMej_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Mej_1stPair_PAS"), getVariableValue("Mej_2ndPair_PAS") );
		FillUserTH2D("h2_MTnuj2_vs_Mej1_highMej_mDeltaPhiMET1stJet_gt_2.5", Me1j1 , MTn1j2 );
		FillUserTH2D("h2_MTnuj1_vs_Mej2_highMej_mDeltaPhiMET1stJet_gt_2.5", Me1j2 , MTn1j1 );
                //---------------------------------------------------------------
		//Comparison - CaloJet vs PFJet - for the 1st PF jet
		TLorentzVector the1stPFjet;
		the1stPFjet.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
					 JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
					 JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),0);
		double myMinDeltaR = 9999;
		double myIdxJet_MinDeltaR = -1;
		for (int myijet=0 ; myijet < CaloJetPt->size() ; myijet++)
		  {
		    TLorentzVector theCaloJet;
		    theCaloJet.SetPtEtaPhiM(CaloJetPt->at(myijet),
					    CaloJetEta->at(myijet),
					    CaloJetPhi->at(myijet),0);

		    double mydeltaR = the1stPFjet.DeltaR(theCaloJet);
		    if(mydeltaR < myMinDeltaR)
		      {
			myMinDeltaR = mydeltaR;
			myIdxJet_MinDeltaR = myijet;
		      }
		  }
		if(myIdxJet_MinDeltaR != -1)
		  {
		    double myPTratio = double(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0])/CaloJetPt->at(myIdxJet_MinDeltaR));
		    double myNCratio = double(PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0]))/double(CaloJetn90Hits->at(myIdxJet_MinDeltaR)) ;
		    // 		    cout << "---------" << endl;
		    // 		    cout << "myMinDeltaR: " << myMinDeltaR << endl;
		    // 		    cout << "PFJetPT/CaloJetPT: " << myPTratio << endl;
		    // 		    cout << "PFJetNC/CaloJetN90Hits: " << myNCratio << endl;
		    // 		    cout << "PFJetNConstituents 1st jet: " << PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0]) << endl;
		    // 		    cout << "CaloJetn90Hits matched CaloJet: " << CaloJetn90Hits->at(myIdxJet_MinDeltaR) << endl;
		    FillUserTH1D("h1_minDR_1stPFjet_CaloJet_highMej_mDeltaPhiMET1stJet_gt_2.5", myMinDeltaR );
		    FillUserTH1D("h1_PT1stJet_over_PTCaloJet_highMej_mDeltaPhiMET1stJet_gt_2.5", myPTratio );
		    FillUserTH1D("h1_NC1stJet_over_n90HitsCaloJet_1stPFjet_CaloJet_highMej_mDeltaPhiMET1stJet_gt_2.5", myNCratio );
		  }
		//---------------------------------------------------------------
                if( getVariableValue("Charge1stEle_PAS")>0 ) {
                  //---------------------------------------------------------------
                  //1D distributions
                  FillUserTH1D("h1_Pt1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Pt1stJet_PAS") );
                  FillUserTH1D("h1_Eta1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Eta1stJet_PAS") );
                  FillUserTH1D("h1_Phi1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Phi1stJet_PAS") );
                  FillUserTH1D("h1_CHF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                  FillUserTH1D("h1_NHF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNeutralHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                  FillUserTH1D("h1_CEF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                  FillUserTH1D("h1_NEF1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNeutralEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                  FillUserTH1D("h1_NCH1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                  FillUserTH1D("h1_NN1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0]) - PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                  FillUserTH1D("h1_NC1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                  FillUserTH1D("h1_Mass1stJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", Mass1stJet );
                  FillUserTH1D("h1_Pt2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Pt2ndJet_PAS") );
                  FillUserTH1D("h1_Eta2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Eta2ndJet_PAS") );
                  FillUserTH1D("h1_Phi2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Phi2ndJet_PAS") );
                  FillUserTH1D("h1_CHF2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                  FillUserTH1D("h1_NHF2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNeutralHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                  FillUserTH1D("h1_CEF2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                  FillUserTH1D("h1_NEF2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNeutralEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                  FillUserTH1D("h1_NCH2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                  FillUserTH1D("h1_NN2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[1]) - PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                  FillUserTH1D("h1_NC2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                  FillUserTH1D("h1_Mass2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", Mass2ndJet );
                  FillUserTH1D("h1_E1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", ElectronEnergy->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
                  FillUserTH1D("h1_Pt1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Pt1stEle_PAS") );
                  FillUserTH1D("h1_Eta1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Eta1stEle_PAS") );
                  FillUserTH1D("h1_Phi1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Phi1stEle_PAS") );
                  FillUserTH1D("h1_Charge1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Charge1stEle_PAS") );
                  FillUserTH1D("h1_mDeltaPhiMETEle_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("mDeltaPhiMETEle") );
                  FillUserTH1D("h1_mDeltaPhiEle1stJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("mDeltaPhiEle1stJet_PAS") );
                  FillUserTH1D("h1_mDeltaPhiEle2ndJet_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("mDeltaPhiEle2ndJet_PAS") );
                  FillUserTH1D("h1_mDeltaPhiMET2ndJet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("mDeltaPhiMET2ndJet") );
                  FillUserTH1D("h1_Ptenu_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Ptenu_PAS") );
                  FillUserTH1D("h1_1stJet_PTOverPTPlusMET_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", (getVariableValue("Pt1stJet_PAS") / (getVariableValue("Pt1stJet_PAS") + getVariableValue("MET"))) );
                  FillUserTH1D("h1_Conversion1stEle_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", (getVariableValue("MissingHits1stEle_PAS")>=1 && getVariableValue("DCotTheta1stEle_PAS")<0.02 && getVariableValue("DCotTheta1stEle_PAS")<0.02) ? 1 : 0 );
                  FillUserTH1D("h1_MET_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("MET_PAS") );
                  FillUserTH1D("h1_CaloMET_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", CaloMET->at(0) );
                  FillUserTH1D("h1_TCMET_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", TCMET->at(0) );
		  FillUserTH1D("h1_MET_over_CaloMET_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("MET_PAS") / CaloMET->at(0) );
                  FillUserTH1D("h1_Njet_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("nJet_PtCut_noOvrlp_ID") );
                  FillUserTH1D("h1_NjetTCHELBTag_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("nJet_TCHELBTag") );
                  FillUserTH1D("h1_minDRej_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("minDRej") );
                  FillUserTH1D("h1_maxDRej_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("maxDRej") );
                  FillUserTH1D("h1_DRjets_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("DeltaRjets_PAS") );
		  FillUserTH1D("h1_MTenu_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("MTenu_PAS") );
                  //---------------------------------------------------------------
                  //2D distributions
                  FillUserTH2D("h2_EtaPhi1stEle_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Eta1stEle_PAS"), getVariableValue("Phi1stEle_PAS") );
                  FillUserTH2D("h2_MinMaxMej_PAS_highMePlusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Mej_1stPair_PAS"), getVariableValue("Mej_2ndPair_PAS") );
                  //---------------------------------------------------------------
                } else {
                  //---------------------------------------------------------------
                  //1D distributions
                  FillUserTH1D("h1_Pt1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Pt1stJet_PAS") );
                  FillUserTH1D("h1_Eta1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Eta1stJet_PAS") );
                  FillUserTH1D("h1_Phi1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Phi1stJet_PAS") );
                  FillUserTH1D("h1_CHF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                  FillUserTH1D("h1_NHF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNeutralHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                  FillUserTH1D("h1_CEF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                  FillUserTH1D("h1_NEF1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNeutralEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                  FillUserTH1D("h1_NCH1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                  FillUserTH1D("h1_NN1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0]) - PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                  FillUserTH1D("h1_NC1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                  FillUserTH1D("h1_Mass1stJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", Mass1stJet );
                  FillUserTH1D("h1_Pt2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Pt2ndJet_PAS") );
                  FillUserTH1D("h1_Eta2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Eta2ndJet_PAS") );
                  FillUserTH1D("h1_Phi2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Phi2ndJet_PAS") );
                  FillUserTH1D("h1_CHF2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                  FillUserTH1D("h1_NHF2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNeutralHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                  FillUserTH1D("h1_CEF2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                  FillUserTH1D("h1_NEF2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNeutralEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                  FillUserTH1D("h1_NCH2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                  FillUserTH1D("h1_NN2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[1]) - PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                  FillUserTH1D("h1_NC2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                  FillUserTH1D("h1_Mass2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", Mass2ndJet );
                  FillUserTH1D("h1_E1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", ElectronEnergy->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
                  FillUserTH1D("h1_Pt1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Pt1stEle_PAS") );
                  FillUserTH1D("h1_Eta1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Eta1stEle_PAS") );
                  FillUserTH1D("h1_Phi1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Phi1stEle_PAS") );
                  FillUserTH1D("h1_Charge1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Charge1stEle_PAS") );
                  FillUserTH1D("h1_mDeltaPhiMETEle_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("mDeltaPhiMETEle") );
                  FillUserTH1D("h1_mDeltaPhiEle1stJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("mDeltaPhiEle1stJet_PAS") );
                  FillUserTH1D("h1_mDeltaPhiEle2ndJet_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("mDeltaPhiEle2ndJet_PAS") );
                  FillUserTH1D("h1_mDeltaPhiMET2ndJet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("mDeltaPhiMET2ndJet") );
                  FillUserTH1D("h1_Ptenu_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Ptenu_PAS") );
                  FillUserTH1D("h1_1stJet_PTOverPTPlusMET_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", (getVariableValue("Pt1stJet_PAS") / (getVariableValue("Pt1stJet_PAS") + getVariableValue("MET"))) );
                  FillUserTH1D("h1_Conversion1stEle_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", (getVariableValue("MissingHits1stEle_PAS")>=1 && getVariableValue("DCotTheta1stEle_PAS")<0.02 && getVariableValue("DCotTheta1stEle_PAS")<0.02) ? 1 : 0 );
                  FillUserTH1D("h1_MET_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("MET_PAS") );
                  FillUserTH1D("h1_CaloMET_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", CaloMET->at(0) );
                  FillUserTH1D("h1_TCMET_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", TCMET->at(0) );
		  FillUserTH1D("h1_MET_over_CaloMET_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("MET_PAS") / CaloMET->at(0) );
                  FillUserTH1D("h1_Njet_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("nJet_PtCut_noOvrlp_ID") );
                  FillUserTH1D("h1_NjetTCHELBTag_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("nJet_TCHELBTag") );
                  FillUserTH1D("h1_minDRej_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("minDRej") );
                  FillUserTH1D("h1_maxDRej_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("maxDRej") );
                  FillUserTH1D("h1_DRjets_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("DeltaRjets_PAS") );
		  FillUserTH1D("h1_MTenu_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("MTenu_PAS") );
                  //---------------------------------------------------------------
                  //2D distributions
                  FillUserTH2D("h2_EtaPhi1stEle_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Eta1stEle_PAS"), getVariableValue("Phi1stEle_PAS") );
                  FillUserTH2D("h2_MinMaxMej_PAS_highMeMinusj_mDeltaPhiMET1stJet_gt_2.5", getVariableValue("Mej_1stPair_PAS"), getVariableValue("Mej_2ndPair_PAS") );
                  //---------------------------------------------------------------
                }
                //---------------------------------------------------------------
                // Event printout
                if( isData ) {

                  STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5: ----------- START ------------");

                  STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
                  if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") && variableIsFilled("Phi1stEle_PAS") )
                    STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5: Pt1stEle_PAS,Eta1stEle_PAS,Phi1stEle_PAS,Charge1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS")<<",\t"<<getVariableValue("Phi1stEle_PAS")<<",\t"<<getVariableValue("Charge1stEle_PAS"));
                  if( variableIsFilled("MET_PAS") && variableIsFilled("METPhi_PAS") )
                    STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5: MET_PAS, METPhi_PAS = "<<getVariableValue("MET_PAS")<<",\t"<<getVariableValue("METPhi_PAS"));
                  if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") && variableIsFilled("Phi1stJet_PAS") )
                    STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5: Pt1stJet_PAS,Eta1stJet_PAS,Phi1stJet = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS")<<",\t"<<getVariableValue("Phi1stJet_PAS"));
                  if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") && variableIsFilled("Phi2ndJet_PAS") )
                    STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5: Pt2ndJet_PAS,Eta2ndJet_PAS,Phi2ndJet = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS")<<",\t"<<getVariableValue("Phi2ndJet_PAS"));
                  if( variableIsFilled("MTenu_PAS") && variableIsFilled("Mjj_PAS") )
                    STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5: MTenu_PAS,Mjj_PAS = "<<getVariableValue("MTenu_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
                  if( variableIsFilled("sT_PAS") )
                    STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5: sT_PAS = "<<getVariableValue("sT_PAS"));
                  if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )
                    STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5: Mej_1stPair_PAS,Mej_2ndPair_PAS = "
                           <<getVariableValue("Mej_1stPair_PAS")
                           <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
                  if( variableIsFilled("MTnuj_1stPair_PAS") && variableIsFilled("MTnuj_2ndPair_PAS") )
                    STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5: MTnuj_1stPair_PAS,MTnuj_2ndPair_PAS = "
                           <<getVariableValue("MTnuj_1stPair_PAS")
                           <<",\t"<<getVariableValue("MTnuj_2ndPair_PAS"));
                  if( variableIsFilled("mDeltaPhiMETEle_PAS") && variableIsFilled("mDeltaPhiMET1stJet_PAS")
                      && variableIsFilled("mDeltaPhiMET2ndJet_PAS") )
                    STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5: mDeltaPhiMETEle_PAS,mDeltaPhiMET1stJet_PAS,mDeltaPhiMET2ndJet_PAS = "
                           <<getVariableValue("mDeltaPhiMETEle_PAS")
                           <<",\t"<<getVariableValue("mDeltaPhiMET1stJet_PAS")
                           <<",\t"<<getVariableValue("mDeltaPhiMET2ndJet_PAS") );
                  if( variableIsFilled("nMuon_PtCut_IDISO") )
                    STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5: nMuon_PtCut_IDISO = "<<getVariableValue("nMuon_PtCut_IDISO"));
                  if( variableIsFilled("minDRej") )
                    STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5: minDRej = "<<getVariableValue("minDRej"));

                  STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5: ------------ END -------------");

                  int NC1stJet = PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0]);
                  if( NC1stJet>=20 && NC1stJet<=30 ) {

                    STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5NC1stJet20To30: ----------- START ------------");

                    STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5NC1stJet20To30: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
                    if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") && variableIsFilled("Phi1stEle_PAS") )
                      STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5NC1stJet20To30: Pt1stEle_PAS,Eta1stEle_PAS,Phi1stEle_PAS,Charge1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS")<<",\t"<<getVariableValue("Phi1stEle_PAS")<<",\t"<<getVariableValue("Charge1stEle_PAS"));
                    if( variableIsFilled("MET_PAS") && variableIsFilled("METPhi_PAS") )
                      STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5NC1stJet20To30: MET_PAS, METPhi_PAS = "<<getVariableValue("MET_PAS")<<",\t"<<getVariableValue("METPhi_PAS"));
                    if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") && variableIsFilled("Phi1stJet_PAS") )
                      STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5NC1stJet20To30: Pt1stJet_PAS,Eta1stJet_PAS,Phi1stJet = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS")<<",\t"<<getVariableValue("Phi1stJet_PAS"));
                    if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") && variableIsFilled("Phi2ndJet_PAS") )
                      STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5NC1stJet20To30: Pt2ndJet_PAS,Eta2ndJet_PAS,Phi2ndJet = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS")<<",\t"<<getVariableValue("Phi2ndJet_PAS"));
                    if( variableIsFilled("MTenu_PAS") && variableIsFilled("Mjj_PAS") )
                      STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5NC1stJet20To30: MTenu_PAS,Mjj_PAS = "<<getVariableValue("MTenu_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
                    if( variableIsFilled("sT_PAS") )
                      STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5NC1stJet20To30: sT_PAS = "<<getVariableValue("sT_PAS"));
                    if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )
                      STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5NC1stJet20To30: Mej_1stPair_PAS,Mej_2ndPair_PAS = "
                             <<getVariableValue("Mej_1stPair_PAS")
                             <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
                    if( variableIsFilled("MTnuj_1stPair_PAS") && variableIsFilled("MTnuj_2ndPair_PAS") )
                      STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5NC1stJet20To30: MTnuj_1stPair_PAS,MTnuj_2ndPair_PAS = "
                             <<getVariableValue("MTnuj_1stPair_PAS")
                             <<",\t"<<getVariableValue("MTnuj_2ndPair_PAS"));
                    if( variableIsFilled("mDeltaPhiMETEle_PAS") && variableIsFilled("mDeltaPhiMET1stJet_PAS")
                        && variableIsFilled("mDeltaPhiMET2ndJet_PAS") )
                      STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5NC1stJet20To30: mDeltaPhiMETEle_PAS,mDeltaPhiMET1stJet_PAS,mDeltaPhiMET2ndJet_PAS = "
                             <<getVariableValue("mDeltaPhiMETEle_PAS")
                             <<",\t"<<getVariableValue("mDeltaPhiMET1stJet_PAS")
                             <<",\t"<<getVariableValue("mDeltaPhiMET2ndJet_PAS") );
                    if( variableIsFilled("nMuon_PtCut_IDISO") )
                      STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5NC1stJet20To30: nMuon_PtCut_IDISO = "<<getVariableValue("nMuon_PtCut_IDISO"));
                    if( variableIsFilled("minDRej") )
                      STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5NC1stJet20To30: minDRej = "<<getVariableValue("minDRej"));

                    STDOUT("UserHighMejDeltaPhiMET1stJetGt2.5NC1stJet20To30: ------------ END -------------");
                  }
                }
              } else {
                //---------------------------------------------------------------
                //1D distributions
                FillUserTH1D("h1_Pt1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("Pt1stJet_PAS") );
                FillUserTH1D("h1_Eta1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("Eta1stJet_PAS") );
                FillUserTH1D("h1_Phi1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("Phi1stJet_PAS") );
                FillUserTH1D("h1_CHF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", PFJetChargedHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                FillUserTH1D("h1_NHF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", PFJetNeutralHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                FillUserTH1D("h1_CEF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", PFJetChargedEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                FillUserTH1D("h1_NEF1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", PFJetNeutralEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                FillUserTH1D("h1_NCH1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                FillUserTH1D("h1_NN1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0]) - PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                FillUserTH1D("h1_NC1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
                FillUserTH1D("h1_Mass1stJet_highMej_mDeltaPhiMET1stJet_le_2.5", Mass1stJet );
                FillUserTH1D("h1_Pt2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("Pt2ndJet_PAS") );
                FillUserTH1D("h1_Eta2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("Eta2ndJet_PAS") );
                FillUserTH1D("h1_Phi2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("Phi2ndJet_PAS") );
                FillUserTH1D("h1_CHF2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", PFJetChargedHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                FillUserTH1D("h1_NHF2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", PFJetNeutralHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                FillUserTH1D("h1_CEF2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", PFJetChargedEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                FillUserTH1D("h1_NEF2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", PFJetNeutralEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                FillUserTH1D("h1_NCH2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                FillUserTH1D("h1_NN2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[1]) - PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                FillUserTH1D("h1_NC2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
                FillUserTH1D("h1_Mass2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", Mass2ndJet );
                FillUserTH1D("h1_E1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", ElectronEnergy->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
                FillUserTH1D("h1_Pt1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("Pt1stEle_PAS") );
                FillUserTH1D("h1_Eta1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("Eta1stEle_PAS") );
                FillUserTH1D("h1_Phi1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("Phi1stEle_PAS") );
                FillUserTH1D("h1_Charge1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("Charge1stEle_PAS") );
                FillUserTH1D("h1_mDeltaPhiMETEle_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("mDeltaPhiMETEle") );
                FillUserTH1D("h1_mDeltaPhiEle1stJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("mDeltaPhiEle1stJet_PAS") );
                FillUserTH1D("h1_mDeltaPhiEle2ndJet_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("mDeltaPhiEle2ndJet_PAS") );
                FillUserTH1D("h1_mDeltaPhiMET2ndJet_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("mDeltaPhiMET2ndJet") );
                FillUserTH1D("h1_Ptenu_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("Ptenu_PAS") );
                FillUserTH1D("h1_1stJet_PTOverPTPlusMET_highMej_mDeltaPhiMET1stJet_le_2.5", (getVariableValue("Pt1stJet_PAS") / (getVariableValue("Pt1stJet_PAS") + getVariableValue("MET"))) );
                FillUserTH1D("h1_Conversion1stEle_highMej_mDeltaPhiMET1stJet_le_2.5", (getVariableValue("MissingHits1stEle_PAS")>=1 && getVariableValue("DCotTheta1stEle_PAS")<0.02 && getVariableValue("DCotTheta1stEle_PAS")<0.02) ? 1 : 0 );
                FillUserTH1D("h1_MET_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("MET_PAS") );
                FillUserTH1D("h1_CaloMET_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", CaloMET->at(0) );
                FillUserTH1D("h1_TCMET_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", TCMET->at(0) );
		FillUserTH1D("h1_MET_over_CaloMET_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("MET_PAS") / CaloMET->at(0) );
                FillUserTH1D("h1_Njet_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("nJet_PtCut_noOvrlp_ID") );
                FillUserTH1D("h1_NjetTCHELBTag_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("nJet_TCHELBTag") );
                FillUserTH1D("h1_minDRej_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("minDRej") );
                FillUserTH1D("h1_maxDRej_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("maxDRej") );
                FillUserTH1D("h1_DRjets_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("DeltaRjets_PAS") );
		FillUserTH1D("h1_MTenu_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("MTenu_PAS") );
                //---------------------------------------------------------------
                //2D distributions
                FillUserTH2D("h2_EtaPhi1stEle_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("Eta1stEle_PAS"), getVariableValue("Phi1stEle_PAS") );
                FillUserTH2D("h2_MinMaxMej_PAS_highMej_mDeltaPhiMET1stJet_le_2.5", getVariableValue("Mej_1stPair_PAS"), getVariableValue("Mej_2ndPair_PAS") );
                //---------------------------------------------------------------
              }
            }
        }
        // Eta1stJet peak at eta~0 after pre-selection
        if( getVariableValue("Eta1stJet_PAS")>0 && getVariableValue("Eta1stJet_PAS")<0.2 ) {
          //---------------------------------------------------------------
          //1D distributions
          FillUserTH1D("h1_Pt1stJet_PAS_Eta1stJetBump", getVariableValue("Pt1stJet_PAS") );
          FillUserTH1D("h1_Phi1stJet_PAS_Eta1stJetBump", getVariableValue("Phi1stJet_PAS") );
          FillUserTH1D("h1_CHF1stJet_Eta1stJetBump", PFJetChargedHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
          FillUserTH1D("h1_NHF1stJet_Eta1stJetBump", PFJetNeutralHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
          FillUserTH1D("h1_CEF1stJet_Eta1stJetBump", PFJetChargedEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
          FillUserTH1D("h1_NEF1stJet_Eta1stJetBump", PFJetNeutralEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
          FillUserTH1D("h1_NCH1stJet_Eta1stJetBump", PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
          FillUserTH1D("h1_NN1stJet_Eta1stJetBump", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0]) - PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
          FillUserTH1D("h1_NC1stJet_Eta1stJetBump", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
          FillUserTH1D("h1_Mass1stJet_Eta1stJetBump", Mass1stJet);
          FillUserTH1D("h1_Pt2ndJet_PAS_Eta1stJetBump", getVariableValue("Pt2ndJet_PAS") );
          FillUserTH1D("h1_Eta2ndJet_PAS_Eta1stJetBump", getVariableValue("Eta2ndJet_PAS") );
          FillUserTH1D("h1_Phi2ndJet_PAS_Eta1stJetBump", getVariableValue("Phi2ndJet_PAS") );
          FillUserTH1D("h1_CHF2ndJet_Eta1stJetBump", PFJetChargedHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
          FillUserTH1D("h1_NHF2ndJet_Eta1stJetBump", PFJetNeutralHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
          FillUserTH1D("h1_CEF2ndJet_Eta1stJetBump", PFJetChargedEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
          FillUserTH1D("h1_NEF2ndJet_Eta1stJetBump", PFJetNeutralEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
          FillUserTH1D("h1_NCH2ndJet_Eta1stJetBump", PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
          FillUserTH1D("h1_NN2ndJet_Eta1stJetBump", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[1]) - PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
          FillUserTH1D("h1_NC2ndJet_Eta1stJetBump", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
          FillUserTH1D("h1_Mass2ndJet_Eta1stJetBump", Mass2ndJet);
          FillUserTH1D("h1_E1stEle_PAS_Eta1stJetBump", ElectronEnergy->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
          FillUserTH1D("h1_Pt1stEle_PAS_Eta1stJetBump", getVariableValue("Pt1stEle_PAS") );
          FillUserTH1D("h1_Eta1stEle_PAS_Eta1stJetBump", getVariableValue("Eta1stEle_PAS") );
          FillUserTH1D("h1_Phi1stEle_PAS_Eta1stJetBump", getVariableValue("Phi1stEle_PAS") );
          FillUserTH1D("h1_Charge1stEle_PAS_Eta1stJetBump", getVariableValue("Charge1stEle_PAS") );
          FillUserTH1D("h1_mDeltaPhiMETEle_Eta1stJetBump", getVariableValue("mDeltaPhiMETEle") );
          FillUserTH1D("h1_mDeltaPhiEle1stJet_PAS_Eta1stJetBump", getVariableValue("mDeltaPhiEle1stJet_PAS") );
          FillUserTH1D("h1_mDeltaPhiEle2ndJet_PAS_Eta1stJetBump", getVariableValue("mDeltaPhiEle2ndJet_PAS") );
          FillUserTH1D("h1_mDeltaPhiMET2ndJet_Eta1stJetBump", getVariableValue("mDeltaPhiMET2ndJet") );
          FillUserTH1D("h1_Ptenu_PAS_Eta1stJetBump", getVariableValue("Ptenu_PAS") );
          FillUserTH1D("h1_1stJet_PTOverPTPlusMET_Eta1stJetBump", (getVariableValue("Pt1stJet_PAS") / (getVariableValue("Pt1stJet_PAS") + getVariableValue("MET"))) );
          FillUserTH1D("h1_Conversion1stEle_Eta1stJetBump", (getVariableValue("MissingHits1stEle_PAS")>=1 && getVariableValue("DCotTheta1stEle_PAS")<0.02 && getVariableValue("DCotTheta1stEle_PAS")<0.02) ? 1 : 0 );
          FillUserTH1D("h1_MET_PAS_Eta1stJetBump", getVariableValue("MET_PAS") );
          FillUserTH1D("h1_CaloMET_PAS_Eta1stJetBump", CaloMET->at(0) );
          FillUserTH1D("h1_TCMET_PAS_Eta1stJetBump", TCMET->at(0) );
	  FillUserTH1D("h1_MET_over_CaloMET_Eta1stJetBump", getVariableValue("MET_PAS") / CaloMET->at(0) );
          FillUserTH1D("h1_Njet_Eta1stJetBump", getVariableValue("nJet_PtCut_noOvrlp_ID") );
          FillUserTH1D("h1_NjetTCHELBTag_Eta1stJetBump", getVariableValue("nJet_TCHELBTag") );
          FillUserTH1D("h1_minDRej_Eta1stJetBump", getVariableValue("minDRej") );
          FillUserTH1D("h1_maxDRej_Eta1stJetBump", getVariableValue("maxDRej") );
          FillUserTH1D("h1_DRjets_Eta1stJetBump", getVariableValue("DeltaRjets_PAS") );
	  FillUserTH1D("h1_MTenu_PAS_Eta1stJetBump", getVariableValue("MTenu_PAS") );
          //---------------------------------------------------------------
          // Event printout
          if( isData ) {

            STDOUT("UserEta1stJetBump: ----------- START ------------");

            STDOUT("UserEta1stJetBump: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
            if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") && variableIsFilled("Phi1stEle_PAS") )
              STDOUT("UserEta1stJetBump: Pt1stEle_PAS,Eta1stEle_PAS,Phi1stEle_PAS,Charge1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS")<<",\t"<<getVariableValue("Phi1stEle_PAS")<<",\t"<<getVariableValue("Charge1stEle_PAS"));
            if( variableIsFilled("MET_PAS") && variableIsFilled("METPhi_PAS") )
              STDOUT("UserEta1stJetBump: MET_PAS, METPhi_PAS = "<<getVariableValue("MET_PAS")<<",\t"<<getVariableValue("METPhi_PAS"));
            if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") && variableIsFilled("Phi1stJet_PAS") )
              STDOUT("UserEta1stJetBump: Pt1stJet_PAS,Eta1stJet_PAS,Phi1stJet = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS")<<",\t"<<getVariableValue("Phi1stJet_PAS"));
            if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") && variableIsFilled("Phi2ndJet_PAS") )
              STDOUT("UserEta1stJetBump: Pt2ndJet_PAS,Eta2ndJet_PAS,Phi2ndJet = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS")<<",\t"<<getVariableValue("Phi2ndJet_PAS"));
            if( variableIsFilled("MTenu_PAS") && variableIsFilled("Mjj_PAS") )
              STDOUT("UserEta1stJetBump: MTenu_PAS,Mjj_PAS = "<<getVariableValue("MTenu_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
            if( variableIsFilled("sT_PAS") )
              STDOUT("UserEta1stJetBump: sT_PAS = "<<getVariableValue("sT_PAS"));
            if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )
              STDOUT("UserEta1stJetBump: Mej_1stPair_PAS,Mej_2ndPair_PAS = "
                     <<getVariableValue("Mej_1stPair_PAS")
                     <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
            if( variableIsFilled("MTnuj_1stPair_PAS") && variableIsFilled("MTnuj_2ndPair_PAS") )
              STDOUT("UserEta1stJetBump: MTnuj_1stPair_PAS,MTnuj_2ndPair_PAS = "
                     <<getVariableValue("MTnuj_1stPair_PAS")
                     <<",\t"<<getVariableValue("MTnuj_2ndPair_PAS"));
            if( variableIsFilled("mDeltaPhiMETEle_PAS") && variableIsFilled("mDeltaPhiMET1stJet_PAS")
                && variableIsFilled("mDeltaPhiMET2ndJet_PAS") )
              STDOUT("UserEta1stJetBump: mDeltaPhiMETEle_PAS,mDeltaPhiMET1stJet_PAS,mDeltaPhiMET2ndJet_PAS = "
                     <<getVariableValue("mDeltaPhiMETEle_PAS")
                     <<",\t"<<getVariableValue("mDeltaPhiMET1stJet_PAS")
                     <<",\t"<<getVariableValue("mDeltaPhiMET2ndJet_PAS") );
            if( variableIsFilled("nMuon_PtCut_IDISO") )
              STDOUT("UserEta1stJetBump: nMuon_PtCut_IDISO = "<<getVariableValue("nMuon_PtCut_IDISO"));
            if( variableIsFilled("minDRej") )
              STDOUT("UserEta1stJetBump: minDRej = "<<getVariableValue("minDRej"));

            STDOUT("UserEta1stJetBump: ------------ END -------------");
          }
        } else if ( (getVariableValue("Eta1stJet_PAS")>-0.6 && getVariableValue("Eta1stJet_PAS")<0) ||
                    (getVariableValue("Eta1stJet_PAS")>0.2 && getVariableValue("Eta1stJet_PAS")<0.6)   )
        {
          //---------------------------------------------------------------
          //1D distributions
          FillUserTH1D("h1_Pt1stJet_PAS_OutsideEta1stJetBump", getVariableValue("Pt1stJet_PAS") );
          FillUserTH1D("h1_Phi1stJet_PAS_OutsideEta1stJetBump", getVariableValue("Phi1stJet_PAS") );
          FillUserTH1D("h1_CHF1stJet_OutsideEta1stJetBump", PFJetChargedHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
          FillUserTH1D("h1_NHF1stJet_OutsideEta1stJetBump", PFJetNeutralHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
          FillUserTH1D("h1_CEF1stJet_OutsideEta1stJetBump", PFJetChargedEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
          FillUserTH1D("h1_NEF1stJet_OutsideEta1stJetBump", PFJetNeutralEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
          FillUserTH1D("h1_NCH1stJet_OutsideEta1stJetBump", PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
          FillUserTH1D("h1_NN1stJet_OutsideEta1stJetBump", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0]) - PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
          FillUserTH1D("h1_NC1stJet_OutsideEta1stJetBump", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
          FillUserTH1D("h1_Mass1stJet_OutsideEta1stJetBump", Mass1stJet );
          FillUserTH1D("h1_Pt2ndJet_PAS_OutsideEta1stJetBump", getVariableValue("Pt2ndJet_PAS") );
          FillUserTH1D("h1_Eta2ndJet_PAS_OutsideEta1stJetBump", getVariableValue("Eta2ndJet_PAS") );
          FillUserTH1D("h1_Phi2ndJet_PAS_OutsideEta1stJetBump", getVariableValue("Phi2ndJet_PAS") );
          FillUserTH1D("h1_CHF2ndJet_OutsideEta1stJetBump", PFJetChargedHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
          FillUserTH1D("h1_NHF2ndJet_OutsideEta1stJetBump", PFJetNeutralHadronEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
          FillUserTH1D("h1_CEF2ndJet_OutsideEta1stJetBump", PFJetChargedEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
          FillUserTH1D("h1_NEF2ndJet_OutsideEta1stJetBump", PFJetNeutralEmEnergyFraction->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
          FillUserTH1D("h1_NCH2ndJet_OutsideEta1stJetBump", PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
          FillUserTH1D("h1_NN2ndJet_OutsideEta1stJetBump", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[1]) - PFJetChargedMultiplicity->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
          FillUserTH1D("h1_NC2ndJet_OutsideEta1stJetBump", PFJetNConstituents->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
          FillUserTH1D("h1_Mass2ndJet_OutsideEta1stJetBump", Mass2ndJet );
          FillUserTH1D("h1_E1stEle_PAS_OutsideEta1stJetBump", ElectronEnergy->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
          FillUserTH1D("h1_Pt1stEle_PAS_OutsideEta1stJetBump", getVariableValue("Pt1stEle_PAS") );
          FillUserTH1D("h1_Eta1stEle_PAS_OutsideEta1stJetBump", getVariableValue("Eta1stEle_PAS") );
          FillUserTH1D("h1_Phi1stEle_PAS_OutsideEta1stJetBump", getVariableValue("Phi1stEle_PAS") );
          FillUserTH1D("h1_Charge1stEle_PAS_OutsideEta1stJetBump", getVariableValue("Charge1stEle_PAS") );
          FillUserTH1D("h1_mDeltaPhiMETEle_OutsideEta1stJetBump", getVariableValue("mDeltaPhiMETEle") );
          FillUserTH1D("h1_mDeltaPhiEle1stJet_PAS_OutsideEta1stJetBump", getVariableValue("mDeltaPhiEle1stJet_PAS") );
          FillUserTH1D("h1_mDeltaPhiEle2ndJet_PAS_OutsideEta1stJetBump", getVariableValue("mDeltaPhiEle2ndJet_PAS") );
          FillUserTH1D("h1_mDeltaPhiMET2ndJet_OutsideEta1stJetBump", getVariableValue("mDeltaPhiMET2ndJet") );
          FillUserTH1D("h1_Ptenu_PAS_OutsideEta1stJetBump", getVariableValue("Ptenu_PAS") );
          FillUserTH1D("h1_1stJet_PTOverPTPlusMET_OutsideEta1stJetBump", (getVariableValue("Pt1stJet_PAS") / (getVariableValue("Pt1stJet_PAS") + getVariableValue("MET"))) );
          FillUserTH1D("h1_Conversion1stEle_OutsideEta1stJetBump", (getVariableValue("MissingHits1stEle_PAS")>=1 && getVariableValue("DCotTheta1stEle_PAS")<0.02 && getVariableValue("DCotTheta1stEle_PAS")<0.02) ? 1 : 0 );
          FillUserTH1D("h1_MET_PAS_OutsideEta1stJetBump", getVariableValue("MET_PAS") );
          FillUserTH1D("h1_CaloMET_PAS_OutsideEta1stJetBump", CaloMET->at(0) );
          FillUserTH1D("h1_TCMET_PAS_OutsideEta1stJetBump", TCMET->at(0) );
	  FillUserTH1D("h1_MET_over_CaloMET_OutsideEta1stJetBump", getVariableValue("MET_PAS") / CaloMET->at(0) );
          FillUserTH1D("h1_Njet_OutsideEta1stJetBump", getVariableValue("nJet_PtCut_noOvrlp_ID") );
          FillUserTH1D("h1_NjetTCHELBTag_OutsideEta1stJetBump", getVariableValue("nJet_TCHELBTag") );
          FillUserTH1D("h1_minDRej_OutsideEta1stJetBump", getVariableValue("minDRej") );
          FillUserTH1D("h1_maxDRej_OutsideEta1stJetBump", getVariableValue("maxDRej") );
          FillUserTH1D("h1_DRjets_OutsideEta1stJetBump", getVariableValue("DeltaRjets_PAS") );
	  FillUserTH1D("h1_MTenu_PAS_OutsideEta1stJetBump", getVariableValue("MTenu_PAS") );
          //---------------------------------------------------------------
        }
      }//end do extra checks


    // Produce skim
    if( passedAllPreviousCuts("minDRej") ) fillSkimTree();

    ////////////////////// User's code to be done for every event - END ///////////////////////

  } // End of loop over events


  ////////////////////// User's code to write histos - BEGIN ///////////////////////

  h1_Ngenlept__ee_emu_etau_eonly_noe_others->Write();

  ////////////////////// User's code to write histos - END ///////////////////////

  delete h1_nPUVertex;
  delete randomNumGen;

  //STDOUT("analysisClass::Loop() ends");
}
