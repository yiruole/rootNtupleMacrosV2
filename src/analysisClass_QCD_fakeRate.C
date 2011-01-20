#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TRandom.h>


//-----------------------------

typedef vector<double>::const_iterator myiter;

struct ordering {
    bool operator ()(pair<size_t, myiter> const& a, pair<size_t, myiter> const& b) {
        return *a.second > *b.second;
    }
};

//-----------------------------
//### JetID ### --> see https://twiki.cern.ch/twiki/bin/view/CMS/ExoticaHighPtJets#JetId

bool JetIdloose(double ak5JetJIDresEMF,double ak5JetJIDfHPD,int ak5JetJIDn90Hits, double ak5JetEta){
  bool jetidloose=false;
  bool jetidresEMF=true;

  double fhpdmax = 0.98;
  double n90hitsmin =1;
  double emf_min = 0.01;

  if(fabs(ak5JetEta)<2.6 && ak5JetJIDresEMF<=emf_min) jetidresEMF=false;

  if(jetidresEMF && ak5JetJIDfHPD<fhpdmax && ak5JetJIDn90Hits>n90hitsmin) {
    jetidloose=true;
  }
  return jetidloose;
}

bool JetIdtight(double ak5JetJIDresEMF,double ak5JetJIDfHPD,int ak5JetJIDn90Hits, double ak5JetEta, double ak5JetPtRaw){
  bool jetidtight=false;
  bool jetidresEMF=true;
  bool jetidfHPD_highPt=true;

  double fhpdmax = 0.98;
  double n90hitsmin =1;
  double emf_min = 0.01;

  if(fabs(ak5JetEta)<2.6 && ak5JetJIDresEMF<=emf_min) jetidresEMF=false;
  if(fabs(ak5JetEta)<2.6 && ak5JetPtRaw>80 && ak5JetJIDresEMF>=1) jetidresEMF=false;
  if(ak5JetPtRaw>25 && ak5JetJIDfHPD>=0.95) jetidfHPD_highPt=false;

  if(jetidresEMF && jetidfHPD_highPt && ak5JetJIDfHPD<fhpdmax && ak5JetJIDn90Hits>n90hitsmin)
    {
      jetidtight=true;
    }
  return jetidtight;
}

bool PFJetIdloose(const double ak5ChargedHadronFraction, const double ak5ChargedEmFraction, const double ak5NeutralHadronFraction, const double ak5NeutralEmFraction, const double ak5JetEta){
  bool jetidloose=false;
  bool jetidChFrac=true;

  double chHadFrac = 0.;
  double chEmFrac = 0.99;
  double neutHadFrac = 0.99;
  double neutEmFrac = 0.99;

  if(fabs(ak5JetEta)<2.4 && ak5ChargedHadronFraction<=chHadFrac && ak5ChargedEmFraction>=chEmFrac) jetidChFrac=false;

  if(jetidChFrac && ak5NeutralHadronFraction<neutHadFrac && ak5NeutralEmFraction<neutEmFrac) {
    jetidloose=true;
  }
  return jetidloose;
}

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
  double eleDeltaEtaTrkSC_bar = getPreCutValue1("eleDeltaEtaTrkSC");
  double eleDeltaEtaTrkSC_end = getPreCutValue2("eleDeltaEtaTrkSC");
  double eleDeltaPhiTrkSC_bar = getPreCutValue1("eleDeltaPhiTrkSC");
  double eleDeltaPhiTrkSC_end = getPreCutValue2("eleDeltaPhiTrkSC");
  double eleHoE_bar = getPreCutValue1("eleHoE");
  double eleHoE_end = getPreCutValue2("eleHoE");
  double eleE2x5OverE5x5_bar = getPreCutValue1("eleE2x5OverE5x5");
  //double eleE2x5OverE5x5_end = getPreCutValue2("eleE2x5OverE5x5");
  double eleE1x5OverE5x5_bar = getPreCutValue1("eleE1x5OverE5x5");
  //double eleE1x5OverE5x5_end = getPreCutValue2("eleE1x5OverE5x5");
  //double eleSigmaIetaIeta_bar = getPreCutValue1("eleSigmaIetaIeta");
  double eleSigmaIetaIeta_end = getPreCutValue2("eleSigmaIetaIeta");
  double eleEcalHcalIso_1_bar = getPreCutValue1("eleEcalHcalIso");
  double eleEcalHcalIso_2_bar = getPreCutValue2("eleEcalHcalIso");
  double eleEcalHcalIso_1_end = getPreCutValue3("eleEcalHcalIso");
  double eleEcalHcalIso_2_end = getPreCutValue4("eleEcalHcalIso");
  //double eleEcalHcalIso_PTthr_bar = getPreCutValue1("eleEcalHcalIso_PTthr");
  double eleEcalHcalIso_PTthr_end = getPreCutValue2("eleEcalHcalIso_PTthr");
  //double eleHcalIsoD2_bar = getPreCutValue1("eleHcalIsoD2");
  double eleHcalIsoD2_end = getPreCutValue2("eleHcalIsoD2");
  double eleTrkIso_bar = getPreCutValue1("eleTrkIso");
  double eleTrkIso_end = getPreCutValue2("eleTrkIso");

  double eleEta_bar_forFR_1 = getPreCutValue1("eleEta_bar_forFR");
  double eleEta_bar_forFR_2 = getPreCutValue2("eleEta_bar_forFR");
  double eleEta_end_forFR_1 = getPreCutValue1("eleEta_end_forFR");
  double eleEta_end_forFR_2 = getPreCutValue2("eleEta_end_forFR");

  double jet_PtCut =    getPreCutValue1("jet_PtCut");
  double jet_PtCut_forMetScale =    getPreCutValue1("jet_PtCut_forMetScale");
  double jet_EtaCut = getPreCutValue1("jet_EtaCut");
  double jet_ele_DeltaRcut =   getPreCutValue1("jet_ele_DeltaRcut");

  double EleEnergyScale_EB=getPreCutValue1("EleEnergyScale_EB");
  double EleEnergyScale_EE=getPreCutValue1("EleEnergyScale_EE");
  double JetEnergyScale=getPreCutValue1("JetEnergyScale");

  int HLTPhoton30          = getPreCutValue1("HLTTrigger");
  int HLTPhoton50          = getPreCutValue2("HLTTrigger");
  int HLTPhoton70          = getPreCutValue3("HLTTrigger");
  int MaxRun_HLTPho30_Only = getPreCutValue1("HLTFromRun");
  int FirstRun_HLTPho70    = getPreCutValue3("HLTFromRun");

  // Not used when using ElectronHeepID and heepBitMask // int eleIDType = (int) getPreCutValue1("eleIDType");
  int heepBitMask_EB  =  getPreCutValue1("heepBitMask_EBGapEE") ;
  int heepBitMask_GAP =  getPreCutValue2("heepBitMask_EBGapEE") ;
  int heepBitMask_EE  =  getPreCutValue3("heepBitMask_EBGapEE") ;

  double muon_PtCut = getPreCutValue1("muon_PtCut");
  double muFidRegion = getPreCutValue1("muFidRegion"); // currently unused !!!
  double muNHits_minThresh = getPreCutValue1("muNHits_minThresh");
  double muTrkD0Maximum = getPreCutValue1("muTrkD0Maximum");

  int jetAlgorithm = getPreCutValue1("jetAlgorithm");
  int metAlgorithm = getPreCutValue1("metAlgorithm");
  int eleAlgorithm = getPreCutValue1("eleAlgorithm"); 

  double vertexMinimumNDOF = getPreCutValue1("vertexMinimumNDOF");
  double vertexMaxAbsZ = getPreCutValue1("vertexMaxAbsZ");
  double vertexMaxd0 = getPreCutValue1("vertexMaxd0");

  double met_Thresh = getPreCutValue1("met_Thresh");
  double minMETPt1stEle_Thresh = getPreCutValue1("minMETPt1stEle_Thresh");
  double Pt1stEle_PAS_Thresh = getPreCutValue1("Pt1stEle_PAS_Thresh");
  double Pt1stJet_PAS_Thresh = getPreCutValue1("Pt1stJet_PAS_Thresh");
  double Pt2ndJet_PAS_Thresh = getPreCutValue1("Pt2ndJet_PAS_Thresh");
  double MTenu_Thresh = getPreCutValue1("MTenu_Thresh");
  double sT_Thresh = getPreCutValue1("sT_Thresh");

  int doExtraChecks = getPreCutValue1("doExtraChecks");

  ////////////////////// User's code to get preCut values - END /////////////////

  ////////////////////// User's code to book histos - BEGIN ///////////////////////

  CreateUserTH1D("h1_SuperClusterPt", getHistoNBins("Pt1stSC"), getHistoMin("Pt1stSC"), getHistoMax("Pt1stSC"));
  CreateUserTH1D("h1_SuperClusterEta", getHistoNBins("Eta1stSC"), getHistoMin("Eta1stSC"), getHistoMax("Eta1stSC"));
  CreateUserTH2D("h2_SuperCluster_Pt_vs_Eta", 
		 getHistoNBins("Eta1stSC"), getHistoMin("Eta1stSC"), getHistoMax("Eta1stSC"),
		 getHistoNBins("Pt1stSC"), getHistoMin("Pt1stSC"), getHistoMax("Pt1stSC"));
  CreateUserTH1D("h1_SuperClusterPt_barrel1", getHistoNBins("Pt1stSC"), getHistoMin("Pt1stSC"), getHistoMax("Pt1stSC"));
  CreateUserTH1D("h1_SuperClusterPt_barrel2", getHistoNBins("Pt1stSC"), getHistoMin("Pt1stSC"), getHistoMax("Pt1stSC"));
  CreateUserTH1D("h1_SuperClusterPt_endcap1", getHistoNBins("Pt1stSC"), getHistoMin("Pt1stSC"), getHistoMax("Pt1stSC"));
  CreateUserTH1D("h1_SuperClusterPt_endcap2", getHistoNBins("Pt1stSC"), getHistoMin("Pt1stSC"), getHistoMax("Pt1stSC"));
  CreateUserTH2D("h2_SuperCluster_Phi_vs_Eta", 
		 getHistoNBins("Eta1stSC"), getHistoMin("Eta1stSC"), getHistoMax("Eta1stSC"),
		 getHistoNBins("Phi1stSC"), getHistoMin("Phi1stSC"), getHistoMax("Phi1stSC"));
  
  CreateUserTH1D("h1_ElePt", getHistoNBins("Pt1stSC"), getHistoMin("Pt1stSC"), getHistoMax("Pt1stSC"));
  CreateUserTH1D("h1_EleEta", getHistoNBins("Eta1stSC"), getHistoMin("Eta1stSC"), getHistoMax("Eta1stSC"));
  CreateUserTH2D("h2_Ele_Pt_vs_Eta", 
		 getHistoNBins("Eta1stSC"), getHistoMin("Eta1stSC"), getHistoMax("Eta1stSC"),
		 getHistoNBins("Pt1stSC"), getHistoMin("Pt1stSC"), getHistoMax("Pt1stSC"));
  CreateUserTH1D("h1_ElePt_barrel1", getHistoNBins("Pt1stSC"), getHistoMin("Pt1stSC"), getHistoMax("Pt1stSC"));
  CreateUserTH1D("h1_ElePt_barrel2", getHistoNBins("Pt1stSC"), getHistoMin("Pt1stSC"), getHistoMax("Pt1stSC"));
  CreateUserTH1D("h1_ElePt_endcap1", getHistoNBins("Pt1stSC"), getHistoMin("Pt1stSC"), getHistoMax("Pt1stSC"));
  CreateUserTH1D("h1_ElePt_endcap2", getHistoNBins("Pt1stSC"), getHistoMin("Pt1stSC"), getHistoMax("Pt1stSC"));
  CreateUserTH1D("h1_Ele_EtOverPt", 2000, 0, 20);
  CreateUserTH1D("h1_Ele_EtOverPt_barrel1", 2000, 0, 20);
  CreateUserTH1D("h1_Ele_EtOverPt_barrel2", 2000, 0, 20);
  CreateUserTH1D("h1_Ele_EtOverPt_endcap1", 2000, 0, 20);
  CreateUserTH1D("h1_Ele_EtOverPt_endcap2", 2000, 0, 20);
  CreateUserTH2D("h2_Ele_Phi_vs_Eta", 
		 getHistoNBins("Eta1stSC"), getHistoMin("Eta1stSC"), getHistoMax("Eta1stSC"),
		 getHistoNBins("Phi1stSC"), getHistoMin("Phi1stSC"), getHistoMax("Phi1stSC"));

  CreateUserTH2D("h2_SCPt_vs_ElePt", 
		 getHistoNBins("Pt1stSC"), getHistoMin("Pt1stSC"), getHistoMax("Pt1stSC"), 
		 getHistoNBins("Pt1stSC"), getHistoMin("Pt1stSC"), getHistoMax("Pt1stSC")); 
  CreateUserTH2D("h2_MET_vs_ElePt", 
		 getHistoNBins("Pt1stSC"), getHistoMin("Pt1stSC"), getHistoMax("Pt1stSC"), 
		 getHistoNBins("MET"), getHistoMin("MET"), getHistoMax("MET"));
  CreateUserTH2D("h2_MET_vs_SumET", 
		 getHistoNBins("SumET"), getHistoMin("SumET"), getHistoMax("SumET"),
		 getHistoNBins("MET"), getHistoMin("MET"), getHistoMax("MET"));

  ////////////////////// User's code to book histos - END ///////////////////////


  Long64_t nentries = fChain->GetEntriesFast();
  STDOUT("analysisClass::Loop(): nentries = " << nentries);

  ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
  ////// If the root version is updated and rootNtupleClass regenerated,     /////
  ////// these lines may need to be updated.                                 /////
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) { // Begin of loop over events
    //for (Long64_t jentry=0; jentry<1000;jentry++) { // Begin of loop over events
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry < 10 || jentry%1000 == 0) STDOUT("analysisClass::Loop(): jentry = " << jentry);
    // if (Cut(ientry) < 0) continue;

    ////////////////////// User's code to be done for every event - BEGIN ///////////////////////


    //## Define new jet collection
    std::auto_ptr<std::vector<double> >  JetPt  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetPtRaw  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetEta  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetPhi  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >     JetPassID  ( new std::vector<int>()  );

    if(jetAlgorithm==1) //PF jets
      {
	for (int ijet=0 ; ijet< PFJetPt->size() ; ijet++)
	  {
	    JetPt->push_back( PFJetPt->at(ijet) );
	    JetPtRaw->push_back( PFJetPtRaw->at(ijet) );
	    JetEta->push_back( PFJetEta->at(ijet) );
	    JetPhi->push_back( PFJetPhi->at(ijet) );
	    JetPassID->push_back(
				 PFJetIdloose(PFJetChargedHadronEnergyFraction->at(ijet),
					      PFJetChargedEmEnergyFraction->at(ijet),
					      PFJetNeutralHadronEnergyFraction->at(ijet),
					      PFJetNeutralEmEnergyFraction->at(ijet),
					      PFJetEta->at(ijet) )
				 );
	  }//end loop over pf jets
      }//end if "pf jets"

    if(jetAlgorithm==2) //Calo jets
      {
	for (int ijet=0 ; ijet < CaloJetPt->size() ; ijet++)
	  {
	    JetPt->push_back( CaloJetPt->at(ijet) );
	    JetPtRaw->push_back( CaloJetPtRaw->at(ijet) );
	    JetEta->push_back( CaloJetEta->at(ijet) );
	    JetPhi->push_back( CaloJetPhi->at(ijet) );
	    JetPassID->push_back(
				 JetIdloose(CaloJetresEMF->at(ijet),
					    CaloJetfHPD->at(ijet),
					    CaloJetn90Hits->at(ijet),
					    CaloJetEta->at(ijet) )
				 );
	  }//end loop over calo jets
      }//end if "calo jets"

    //## Define new met collection
    double thisMET;
    double thisSumET;
    double thisMETPhi;
    if(metAlgorithm==1) 	// --> PFMET
      {
	thisMET = PFMET->at(0);
	thisSumET = PFSumET->at(0);
	thisMETPhi = PFMETPhi->at(0);
      }
    if(metAlgorithm==2) 	// --> CaloMET
      {
	thisMET = CaloMET->at(0);
	thisSumET = CaloSumET->at(0);
	thisMETPhi = CaloMETPhi->at(0);
      }
    if(metAlgorithm==3) 	// --> PFMET (with type-1 corrections)
      {
	thisMET = PFMETType1Cor->at(0);
	thisMETPhi = PFMETPhiType1Cor->at(0);
      }
    // --> TCMET
    //     thisMET = TCMET->at(0);
    //     thisMETPhi = TCMETPhi->at(0);


     //## EES and JES
    if( EleEnergyScale_EB != 1 || EleEnergyScale_EE != 1 )
      {
	for(int iele=0; iele<SuperClusterPt->size(); iele++)
	  {
	    if( fabs(SuperClusterEta->at(iele)) < eleEta_bar )
	      SuperClusterPt->at(iele) *= EleEnergyScale_EB;
	    if( fabs(SuperClusterEta->at(iele)) > eleEta_end_min && fabs(SuperClusterEta->at(iele)) < eleEta_end_max )
	      SuperClusterPt->at(iele) *= EleEnergyScale_EE;
	  }
      }
    if( JetEnergyScale != 1 )
      { //use fix JES scaling passed from cut file

	for(int ijet=0; ijet<JetPt->size(); ijet++)
	  {
	    JetPt->at(ijet) *= JetEnergyScale;
	  }
      }

    //## Run selection
    int passGoodRunList = 1;

    //     if(isData)
    //       {       
    // 	//select first period of run 2010 (low luminosity)
    // 	if( run > 144114 ) passGoodRunList = 0;
    
    // 	//select second period of run 2010 (high luminosity) 
    // 	//if( run <= 144114 ) passGoodRunList = 0;		
    //       } 

    //## HLT
    //--> OLD HLT SELECTION BY RUN NUMBER
    //     int PassTrig = 0;
    //     int HLTFromRun[4] = {getPreCutValue1("HLTFromRun"),
    // 			 getPreCutValue2("HLTFromRun"),
    // 			 getPreCutValue3("HLTFromRun"),
    // 			 getPreCutValue4("HLTFromRun")};
    //     int HLTTrigger[4] = {getPreCutValue1("HLTTrigger"),
    // 			 getPreCutValue2("HLTTrigger"),
    // 			 getPreCutValue3("HLTTrigger"),
    // 			 getPreCutValue4("HLTTrigger")};
    //     int HLTTrgUsed;
    //     for (int i=0; i<4; i++) {
    //       if ( !isData && i != 0) continue; // For MC use HLTPhoton15 as the cleaned trigger is not in MC yet as of July 20, 2010
    //       if ( HLTFromRun[i] <= run ) {
    //  	//if(jentry == 0 ) STDOUT("run, i, HLTTrigger[i], HLTFromRun[i] = "<<run<<"\t"<<i<<"\t"<<"\t"<<HLTTrigger[i]<<"\t"<<HLTFromRun[i]);
    // 	if (HLTTrigger[i] > 0 && HLTTrigger[i] < HLTResults->size() ) {
    // 	  PassTrig=HLTResults->at(HLTTrigger[i]);
    // 	  HLTTrgUsed=HLTTrigger[i];
    // 	} else {
    // 	  STDOUT("ERROR: HLTTrigger out of range of HLTResults: HLTTrigger = "<<HLTTrigger[i] <<"and HLTResults size = "<< HLTResults->size());
    // 	}
    //       }
    //     }
    //     if(jentry == 0 ) STDOUT("Run = "<<run <<", HLTTrgUsed is number = "<<HLTTrgUsed<<" of the list HLTPathsOfInterest");

    //## Superclusters
    vector<int> v_idx_sc_all;
    vector<int> v_idx_sc_PtCut;
    vector<int> v_idx_sc_Iso;

    //Create vector with indices of supercluster ordered by pT
    vector<pair<size_t, myiter> > order(SuperClusterPt->size());
    size_t n = 0;
    for (myiter it = SuperClusterPt->begin(); it != SuperClusterPt->end(); ++it, ++n)
      order[n] = make_pair(n, it);
    sort(order.begin(), order.end(), ordering());

    for(int isc=0; isc<order.size(); isc++)
      {
	// 	cout << "index , pT: "
	// 	     << order[isc].first
	// 	     << " , "
	// 	     << *order[isc].second
	// 	     << endl;
	v_idx_sc_all.push_back(order[isc].first); //### All superclusters ordered by pT
      }

    for(int isc=0;isc<v_idx_sc_all.size();isc++){

      //pT cut + ECAL acceptance cut + remove spikes (all together)
      if ( 1 - SuperClusterS4S1->at(v_idx_sc_all[isc]) > 0.95 ) continue;

      bool Barrel = false;
      bool Endcap = false;
      if (fabs(SuperClusterEta->at(v_idx_sc_all[isc]))<eleEta_bar) Barrel = true;
      if ((fabs(SuperClusterEta->at(v_idx_sc_all[isc]))<eleEta_end_max)
	  &&(fabs(SuperClusterEta->at(v_idx_sc_all[isc]))>eleEta_end_min)) Endcap = true;
      if ( !Barrel && !Endcap) continue;

      bool PassPt = false;
      if ( SuperClusterPt->at(v_idx_sc_all[isc]) > getPreCutValue1("ele_PtCut") ) PassPt=true;
      if ( !PassPt ) continue;

      v_idx_sc_PtCut.push_back(v_idx_sc_all[isc]); //### Pt cut (+ no gaps + no spikes)

      //ID+Isolation together
      bool PassHoE = false;
      if ( SuperClusterHoE->at(v_idx_sc_all[isc])<0.05) PassHoE=true;
      if ( !PassHoE ) continue;

      bool PassEcalIso = false;
      if (Barrel && SuperClusterHEEPEcalIso->at(v_idx_sc_all[isc]) <(6+(0.01*SuperClusterPt->at(v_idx_sc_all[isc]))))
	PassEcalIso=true;
      if (Endcap && SuperClusterPt->at(v_idx_sc_all[isc])<50
	  && SuperClusterHEEPEcalIso->at(v_idx_sc_all[isc])<(6+(0.01*SuperClusterPt->at(v_idx_sc_all[isc]))))
	PassEcalIso=true;
      if (Endcap && SuperClusterPt->at(v_idx_sc_all[isc])>=50
	  && SuperClusterHEEPEcalIso->at(v_idx_sc_all[isc])<(6+(0.01*(SuperClusterPt->at(v_idx_sc_all[isc])-50))))
	PassEcalIso=true;
      if ( !PassEcalIso ) continue;

      v_idx_sc_Iso.push_back(v_idx_sc_all[isc]); //### Pt cut + ID+ISO

    }

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
	
	if(eleAlgorithm==1) 	//------> use HEEP mask
	  {
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
	    
	  }//end use HEEP mask

	else if (eleAlgorithm==2) //------> use HEEP variables
	  {
	    int passEleSel = 0;
	    int isBarrel = 0;
	    int isEndcap = 0;

	    if( fabs( ElectronSCEta->at(iele) ) < eleEta_bar ) 
	      isBarrel = 1;
	    
	    if( fabs( ElectronSCEta->at(iele) ) > eleEta_end_min && fabs( ElectronSCEta->at(iele) ) < eleEta_end_max ) 
	      isEndcap = 1;
	    
	    if(isBarrel)
	      {
		if(fabs(ElectronDeltaEtaTrkSC->at(iele)) < eleDeltaEtaTrkSC_bar //0.005
		   && fabs(ElectronDeltaPhiTrkSC->at(iele)) < eleDeltaPhiTrkSC_bar //0.09
		   && ElectronHoE->at(iele) < eleHoE_bar //0.05 
		   && (ElectronE2x5OverE5x5->at(iele) >eleE2x5OverE5x5_bar || ElectronE1x5OverE5x5->at(iele) > eleE1x5OverE5x5_bar ) //0.94  , 0.83
		   && ElectronEcalIsoHeep->at(iele)+ElectronHcalIsoD1Heep->at(iele) < eleEcalHcalIso_1_bar + eleEcalHcalIso_2_bar*ElectronPt->at(iele)  // 2+0.03*ElectronPt->at(iele)
		   && ElectronTrkIsoHeep->at(iele) <eleTrkIso_bar //7.5
		   )
		  passEleSel = 1;		
	      }//end barrel
	    
	    if(isEndcap)
	      {		
		int passEcalHcalIsoCut=0;
		if(ElectronPt->at(iele) < eleEcalHcalIso_PTthr_end // thr=50 
		   && (ElectronEcalIsoHeep->at(iele)+ElectronHcalIsoD1Heep->at(iele)) < eleEcalHcalIso_1_end) // value=2.5 
		  passEcalHcalIsoCut=1;
		if(ElectronPt->at(iele) > eleEcalHcalIso_PTthr_end // thr=50 
		   && (ElectronEcalIsoHeep->at(iele)+ElectronHcalIsoD1Heep->at(iele)) < eleEcalHcalIso_1_end+eleEcalHcalIso_2_end*(ElectronPt->at(iele)-eleEcalHcalIso_PTthr_end) ) //values=2.5, 0.03
		  passEcalHcalIsoCut=1;
		
		if(fabs(ElectronDeltaEtaTrkSC->at(iele)) < eleDeltaEtaTrkSC_end //0.007
		   && fabs(ElectronDeltaPhiTrkSC->at(iele)) < eleDeltaPhiTrkSC_end //0.09
		   && ElectronHoE->at(iele) < eleHoE_end //0.05 
		   && ElectronSigmaIEtaIEta->at(iele) < eleSigmaIetaIeta_end //0.03
		   && passEcalHcalIsoCut == 1
		   && ElectronHcalIsoD2Heep->at(iele) < eleHcalIsoD2_end //0.5
		   && ElectronTrkIsoHeep->at(iele) < eleTrkIso_end //15
		   )
		  passEleSel = 1;		
	      }//end endcap

	    //Pass User Defined Electron Selection
	    if ( passEleSel )
	      {
		v_idx_ele_PtCut_IDISO_noOverlap.push_back(iele);
	      }

	  }//end use HEEP variables

      } // End loop over electrons


    //## Jets
    vector<int> v_idx_jet_all;
    vector<int> v_idx_jet_PtCut;
    vector<int> v_idx_jet_PtCut_noOverlap;
    vector<int> v_idx_jet_PtCut_noOverlap_ID;
    vector<int> v_idx_jet_PtCut_noOverlap_ID_EtaCut;

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
    for (int isc=0; isc<v_idx_sc_Iso.size(); isc++)
      {
	TLorentzVector sc;
        sc.SetPtEtaPhiM(SuperClusterPt->at(v_idx_sc_Iso[isc]),
			SuperClusterEta->at(v_idx_sc_Iso[isc]),
			SuperClusterPhi->at(v_idx_sc_Iso[isc]),0);
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
	    double DR = jet.DeltaR(sc);
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


    //## HLT
    //--> NEW HLT SELECTION USING DYNAMIC PRESCALE
    int PassTrig = 0;
    double weight_HLT = 1;

    //HLT for MC samples
    if(isData==0)
      {
	PassTrig = 1;
      }

    
    //### TO BE COMMENTED OUT IF RUNNING ON SPRING10 MC NTUPLES !!! ###
    if(isData)
      {
	//cout << run << ", " << ls << ", " << event << endl;
	if( v_idx_sc_Iso.size()>=1 )
	  {
	    if( run <= MaxRun_HLTPho30_Only ) //---------------------------------> photon 30 only
	      {
		if( HLTPrescales->at( HLTPhoton30 ) == -1 )
		  HLTPrescales->at( HLTPhoton30 ) = 1 ;
		
		if( HLTResults->at( HLTPhoton30 ) ) //this trigger is un-prescaled in this run range
		  {
		    PassTrig = 1;
		    weight_HLT = HLTPrescales->at( HLTPhoton30 );
		    //cout << "Pho30 1st period" << endl;
		  }
	      }
	    else if(run > MaxRun_HLTPho30_Only && run < FirstRun_HLTPho70) //--> photon 30 and photon50
	      {
		if( HLTPrescales->at( HLTPhoton50 ) == -1 )
		  HLTPrescales->at( HLTPhoton50 ) = 1 ;
		if( HLTPrescales->at( HLTPhoton30 ) == -1 )
		  HLTPrescales->at( HLTPhoton30 ) = 1 ;
		
		if( HLTResults->at( HLTPhoton50 ) ) //this trigger is un-prescaled in this run range
		  {
		    PassTrig = 1;
		    weight_HLT = HLTPrescales->at( HLTPhoton50 );
		  }
		else if( HLTResults->at( HLTPhoton30 ) && SuperClusterPt->at(v_idx_sc_Iso[0])>30 && SuperClusterPt->at(v_idx_sc_Iso[0])<50 )
		  {
		    PassTrig = 1;
		    weight_HLT = HLTPrescales->at( HLTPhoton30 );
		  }
	      }
	    else //-------------------------------------------------------------> photon 70, photon50, and photon 30
	      {
		// 	    cout << "Pho70" << endl;
		// 	    cout << HLTResults->at( HLTPhoton70 ) << endl;
		// 	    cout << "Pho50" << endl;
		// 	    cout << HLTResults->at( HLTPhoton50 ) << endl;
		// 	    cout << "Pho30" << endl;
		// 	    cout << HLTResults->at( HLTPhoton30 ) << endl;
		if( HLTPrescales->at( HLTPhoton70 ) == -1 )
		  HLTPrescales->at( HLTPhoton70 ) = 1 ;
		if( HLTPrescales->at( HLTPhoton50 ) == -1 )
		  HLTPrescales->at( HLTPhoton50 ) = 1 ;
		if( HLTPrescales->at( HLTPhoton30 ) == -1 )
		  HLTPrescales->at( HLTPhoton30 ) = 1 ;
		
		if( HLTResults->at( HLTPhoton70 ) ) //this trigger is un-prescaled in this run range
		  {
		    PassTrig = 1;
		    weight_HLT = HLTPrescales->at( HLTPhoton70 );
		  }
		else if( HLTResults->at( HLTPhoton50 ) && SuperClusterPt->at(v_idx_sc_Iso[0])>50 && SuperClusterPt->at(v_idx_sc_Iso[0])<70 )
		  {
		    PassTrig = 1;
		    weight_HLT = HLTPrescales->at( HLTPhoton50 );
		  }
		else if( HLTResults->at( HLTPhoton30 ) && SuperClusterPt->at(v_idx_sc_Iso[0])>30 && SuperClusterPt->at(v_idx_sc_Iso[0])<50 )
		  {
		    PassTrig = 1;
		    weight_HLT = HLTPrescales->at( HLTPhoton30 );
		  }
	      }
	  }
      }//end use of HLT prescales
    

    //     //## Vertexes
    //     vector<int> v_idx_vertex_good;
    //     // loop over vertexes
    //     for(int ivertex = 0; ivertex<VertexChi2->size(); ivertex++){
    //       if ( !(VertexIsFake->at(ivertex))
    //     	   && VertexNDF->at(ivertex) > vertexMinimumNDOF
    //     	   && fabs( VertexZ->at(ivertex) ) <= vertexMaxAbsZ
    //     	   && fabs( VertexRho->at(ivertex) ) <= vertexMaxd0 )
    //     	{
    //     	  v_idx_vertex_good.push_back(ivertex);
    //     	  //STDOUT("v_idx_vertex_good.size = "<< v_idx_vertex_good.size() );
    //     	}
    //     }


    //FINAL EVENT WEIGHT: event weight from HLT prescale
    double p1 = 0;
    p1 = weight_HLT;
    //don't apply any fake rate probability
    //p1 = 1;


    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();

    //-----

    // Set the value of the variableNames listed in the cutFile to their current value


    /// Now fill all variables that are filled just once
    // Trigger (L1 and HLT)
    if(isData==true)
      {
        fillVariableWithValue( "PassGoodRunList", passGoodRunList );
	fillVariableWithValue( "PassBPTX0", isBPTX0 ) ;
	fillVariableWithValue( "PassPhysDecl", isPhysDeclared ) ;
	fillVariableWithValue( "PassHLT", PassTrig, p1 ) ;
      }
    else
      {
	fillVariableWithValue( "PassHLT", PassTrig, 1 ) ;
        fillVariableWithValue( "PassGoodRunList", 1 );
	fillVariableWithValue( "PassBPTX0", true ) ;
	fillVariableWithValue( "PassPhysDecl", true ) ;
      }

    //fillVariableWithValue( "nVertex_good", v_idx_vertex_good.size() ) ;

    //Event filters at RECO level
    fillVariableWithValue( "PassBeamScraping", !isBeamScraping ) ;
    fillVariableWithValue( "PassPrimaryVertex", isPrimaryVertex ) ;

    //Fill number of isolated superclusters per event
    fillVariableWithValue( "nIsoSC" , v_idx_sc_Iso.size(), p1 ) ;

    // nJet
    fillVariableWithValue( "nJet_PtCut_noOvrlp_ID", v_idx_jet_PtCut_noOverlap_ID.size(), p1 ) ;

    //Fill number of isolated electrons per event
    fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp" , v_idx_ele_PtCut_IDISO_noOverlap.size(), p1 ) ;

    // MET
    fillVariableWithValue("MET", thisMET, p1);
    fillVariableWithValue("METsel", thisMET, p1);
    fillVariableWithValue("SumET", thisSumET, p1);
    fillVariableWithValue("METPhi", thisMETPhi, p1);

    // 1st sc
    if( v_idx_sc_Iso.size() >= 1 )
      {
	fillVariableWithValue( "Pt1stSC", SuperClusterPt->at(v_idx_sc_Iso[0]), p1 );
	fillVariableWithValue( "Eta1stSC", SuperClusterEta->at(v_idx_sc_Iso[0]), p1 );
	fillVariableWithValue( "Phi1stSC", SuperClusterPhi->at(v_idx_sc_Iso[0]), p1 );

	//minDRscjets
	TLorentzVector sc;
	sc.SetPtEtaPhiM(SuperClusterPt->at(v_idx_sc_Iso[0]),
			SuperClusterEta->at(v_idx_sc_Iso[0]),
			SuperClusterPhi->at(v_idx_sc_Iso[0]),0);
	TLorentzVector jet;
	double minDR=9999.;
	int ijet_minDR = -1;
	for(int ijet=0; ijet<v_idx_jet_PtCut_noOverlap_ID.size(); ijet++)
	  {
	    jet.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
			     JetEta->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
			     JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),0);
	    double DR = jet.DeltaR(sc);
	    if (DR<minDR)
	      {
		minDR = DR;	       
		ijet_minDR = ijet;
	      }
	  }//end loop over jets	
	fillVariableWithValue( "minDRscjets", minDR, p1 );

	//deltaPhi sc and 1st jet
	double deltaphi = 999;
	if( v_idx_jet_PtCut_noOverlap_ID.size() >= 1 )
	  {
	    TLorentzVector jet1;
	    jet1.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			      JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			      JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),0);		
	    deltaphi = sc.DeltaPhi(jet1);
	  }
	fillVariableWithValue( "mDeltaPhiSC1stJet", fabs(deltaphi), p1);	    

	//scEToverMET
	fillVariableWithValue( "scEToverMET", double(SuperClusterPt->at(v_idx_sc_Iso[0])/thisMET), p1);	    

      }

    // 1st ele and transverse mass enu
    double MT, DeltaPhiMETEle = -999;
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 1 && v_idx_sc_Iso.size() >= 1 )
      {
	fillVariableWithValue( "Pt1stEle", ElectronSCPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]), p1 );
	fillVariableWithValue( "Eta1stEle", ElectronSCEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]), p1 );
	fillVariableWithValue( "Phi1stEle", ElectronSCPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]), p1 );

	// DeltaPhi - MET vs 1st ele
	TVector2 v_MET;
	TVector2 v_ele;
	v_MET.SetMagPhi( thisMET , thisMETPhi);
	v_ele.SetMagPhi( ElectronSCPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) , ElectronSCPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	float deltaphi = v_MET.DeltaPhi(v_ele);
	fillVariableWithValue( "mDeltaPhiMETEle", fabs(deltaphi), p1);
        DeltaPhiMETEle = fabs(deltaphi);

	// transverse mass enu
	MT = sqrt(2 * ElectronSCPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) * thisMET * (1 - cos(deltaphi)) );
	fillVariableWithValue("MTenu", MT, p1);

	// electron Charge
	fillVariableWithValue( "Charge1stEle", ElectronCharge->at(v_idx_ele_PtCut_IDISO_noOverlap[0]), p1);

	// deltaR ele, sc
	TLorentzVector ele, sc;
	ele.SetPtEtaPhiM(ElectronSCPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronSCEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronSCPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),0);
	sc.SetPtEtaPhiM(SuperClusterPt->at(v_idx_sc_Iso[0]),
			SuperClusterEta->at(v_idx_sc_Iso[0]),
			SuperClusterPhi->at(v_idx_sc_Iso[0]),0);
	double deltaR_ele_sc = ele.DeltaR(sc);
	fillVariableWithValue( "DeltaR_ele_sc", deltaR_ele_sc , p1);
	
	//deltaPhi ele and 1st jet
	double deltaphi2 = 999;
	if( v_idx_jet_PtCut_noOverlap_ID.size() >= 1 )
	  {
	    TLorentzVector jet1;
	    jet1.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			      JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			      JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),0);		
	    deltaphi2 = ele.DeltaPhi(jet1);
	  }
	fillVariableWithValue( "mDeltaPhiEle1stJet", fabs(deltaphi2), p1);	    

      }

    // 1st jet and deltaphi jet-MET
    double DeltaPhiMET1stJet = -999;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 1 )
      {
	fillVariableWithValue( "Pt1stJet", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]), p1 );
	fillVariableWithValue( "Eta1stJet", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]), p1 );

	//DeltaPhi - MET vs 1st jet
	TVector2 v_MET;
	TVector2 v_jet;
	v_MET.SetMagPhi( 1 , thisMETPhi );
	v_jet.SetMagPhi( 1 , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	float deltaphi = v_MET.DeltaPhi(v_jet);
	fillVariableWithValue( "mDeltaPhiMET1stJet", fabs(deltaphi), p1 );
        DeltaPhiMET1stJet = fabs(deltaphi);
      }


    // Evaluate cuts (but do not apply them)
    evaluateCuts();

    // Fill histograms and do analysis based on cut evaluation

    //pass supercluster selection
    if( passedAllPreviousCuts("nEle_PtCut_IDISO_noOvrlp") 
	&& variableIsFilled("Pt1stSC")  
	&& variableIsFilled("Eta1stSC")  
	&& variableIsFilled("Phi1stSC")  
	)
      {
	FillUserTH1D("h1_SuperClusterPt", getVariableValue("Pt1stSC"), p1);		     
	FillUserTH1D("h1_SuperClusterEta", getVariableValue("Eta1stSC"), p1);		     
	FillUserTH2D("h2_SuperCluster_Pt_vs_Eta", getVariableValue("Eta1stSC"), getVariableValue("Pt1stSC"), p1 );
	FillUserTH2D("h2_SuperCluster_Phi_vs_Eta", getVariableValue("Eta1stSC") , getVariableValue("Phi1stSC") , p1);

	if( fabs(getVariableValue("Eta1stSC")) < eleEta_bar_forFR_1  )//barrel 1
	  {
	    FillUserTH1D("h1_SuperClusterPt_barrel1", getVariableValue("Pt1stSC"), p1);		     
	  }
	if( fabs(getVariableValue("Eta1stSC")) >= eleEta_bar_forFR_1 && fabs(getVariableValue("Eta1stSC")) < eleEta_bar_forFR_2 )//barrel 2
	  {
	    FillUserTH1D("h1_SuperClusterPt_barrel2", getVariableValue("Pt1stSC"), p1);		     
	  }
	if( fabs(getVariableValue("Eta1stSC")) >= eleEta_end_forFR_1 && fabs(getVariableValue("Eta1stSC")) < eleEta_end_forFR_2 )//endcap 1
	  {
	    FillUserTH1D("h1_SuperClusterPt_endcap1", getVariableValue("Pt1stSC"), p1);		     
	  }
	if( fabs(getVariableValue("Eta1stSC")) >= eleEta_end_forFR_2 )//endcap 2
	  {
	    FillUserTH1D("h1_SuperClusterPt_endcap2", getVariableValue("Pt1stSC"), p1);		     
	  }
      }

    //pass electron selection
    if( passedCut("all") 
	&& variableIsFilled("Pt1stEle")  
	&& variableIsFilled("Eta1stEle") 
	&& variableIsFilled("Phi1stEle") 
	&& variableIsFilled("Pt1stSC")  
	&& variableIsFilled("Eta1stSC")  
	&& variableIsFilled("Phi1stSC")  
	&& variableIsFilled("MET")
	&& variableIsFilled("SumET")
	)
      {
	double EToverPT = double( getVariableValue("Pt1stEle") / ElectronTrackPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );

	FillUserTH1D("h1_ElePt", getVariableValue("Pt1stEle"), p1);		     
	FillUserTH1D("h1_EleEta", getVariableValue("Eta1stEle"), p1);		      
	FillUserTH2D("h2_Ele_Pt_vs_Eta", getVariableValue("Eta1stEle"), getVariableValue("Pt1stEle"), p1 );
	FillUserTH1D("h1_Ele_EtOverPt", EToverPT, p1 );
	FillUserTH2D("h2_Ele_Phi_vs_Eta", getVariableValue("Eta1stEle") , getVariableValue("Phi1stEle") , p1);

	if( fabs(getVariableValue("Eta1stSC")) < eleEta_bar_forFR_1  )//barrel 1
	  {
	    FillUserTH1D("h1_ElePt_barrel1", getVariableValue("Pt1stEle"), p1);		     
	    FillUserTH1D("h1_Ele_EtOverPt_barrel1", EToverPT, p1 );
	  }
	if( fabs(getVariableValue("Eta1stSC")) >= eleEta_bar_forFR_1 && fabs(getVariableValue("Eta1stSC")) < eleEta_bar_forFR_2 )//barrel 2
	  {
	    FillUserTH1D("h1_ElePt_barrel2", getVariableValue("Pt1stEle"), p1);		     
	    FillUserTH1D("h1_Ele_EtOverPt_barrel2", EToverPT, p1 );
	  }
	if( fabs(getVariableValue("Eta1stSC")) >= eleEta_end_forFR_1 && fabs(getVariableValue("Eta1stSC")) < eleEta_end_forFR_2 )//endcap 1
	  {
	    FillUserTH1D("h1_ElePt_endcap1", getVariableValue("Pt1stEle"), p1);		     
	    FillUserTH1D("h1_Ele_EtOverPt_endcap1", EToverPT, p1 );
	  }
	if( fabs(getVariableValue("Eta1stSC")) >= eleEta_end_forFR_2 )//endcap 2
	  {
	    FillUserTH1D("h1_Ele_EtOverPt_endcap2", EToverPT, p1 );
	    FillUserTH1D("h1_ElePt_endcap2", getVariableValue("Pt1stEle"), p1);		     
	  }

	FillUserTH2D("h2_SCPt_vs_ElePt", getVariableValue("Pt1stEle") , getVariableValue("Pt1stSC") , p1 );
	FillUserTH2D("h2_MET_vs_ElePt", getVariableValue("Pt1stEle") , getVariableValue("MET") , p1 );
	FillUserTH2D("h2_MET_vs_SumET", getVariableValue("SumET") , getVariableValue("MET") , p1);
      }


    //     // Produce skim
    //     if( passedAllPreviousCuts("minDRej") ) fillSkimTree();

    //INFO
    //      // retrieve value of previously filled variables (after making sure that they were filled)
    //      double totpTEle;
    //      if ( variableIsFilled("pT1stEle") && variableIsFilled("pT2ndEle") )
    //        totpTEle = getVariableValue("pT1stEle")+getVariableValue("pT2ndEle");
    //      // reject events that did not pass level 0 cuts
    //      if( !passedCut("0") ) continue;
    //      // ......

    ////////////////////// User's code to be done for every event - END ///////////////////////

  } // End of loop over events

  ////////////////////// User's code to write histos - BEGIN ///////////////////////

  ////////////////////// User's code to write histos - END ///////////////////////


  //STDOUT("analysisClass::Loop() ends");
}
