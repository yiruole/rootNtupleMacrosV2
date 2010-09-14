#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>




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
   
  ////////////////////// User's code to book histos - BEGIN ///////////////////////

  TH1F *h_Mej_PAS = new TH1F ("h_Mej_PAS","h_Mej_PAS",200,0,2000);  h_Mej_PAS->Sumw2();

  TH1F *h_ElectronDeltaEtaTrkSC_barrel_heep = new TH1F ("h_ElectronDeltaEtaTrkSC_barrel_heep","h_ElectronDeltaEtaTrkSC_barrel_heep",200,-0.05,0.05);
  TH1F *h_ElectronDeltaPhiTrkSC_barrel_heep = new TH1F ("h_ElectronDeltaPhiTrkSC_barrel_heep","h_ElectronDeltaPhiTrkSC_barrel_heep",200,-0.5,0.5);
  TH1F *h_ElectronHoE_barrel_heep = new TH1F ("h_ElectronHoE_barrel_heep","h_ElectronHoE_barrel_heep",75,0,0.15);
  TH1F *h_ElectronSigmaIEtaIEta_barrel_heep = new TH1F ("h_ElectronSigmaIEtaIEta_barrel_heep","h_ElectronSigmaIEtaIEta_barrel_heep",100,0,0.1);
  TH1F *h_ElectronE2x5OverE5x5_barrel_heep = new TH1F ("h_ElectronE2x5OverE5x5_barrel_heep","h_ElectronE2x5OverE5x5_barrel_heep",100,0,1);
  TH1F *h_ElectronE1x5OverE5x5_barrel_heep = new TH1F ("h_ElectronE1x5OverE5x5_barrel_heep","h_ElectronE1x5OverE5x5_barrel_heep",100,0,1);
  TH1F *h_ElectronEcalIsoHeep_barrel_heep = new TH1F ("h_ElectronEcalIsoHeep_barrel_heep","h2_ElectronEcalIsoHeep_barrel_heep",100,0,50);
  TH1F *h_ElectronHcalIsoD1Heep_barrel_heep = new TH1F ("h_ElectronHcalIsoD1Heep_barrel_heep","h_ElectronHcalIsoD1Heep_barrel_heep",100,0,50);
  TH1F *h_ElectronHcalIsoD2Heep_barrel_heep = new TH1F ("h_ElectronHcalIsoD2Heep_barrel_heep","h_ElectronHcalIsoD2Heep_barrel_heep",100,0,50);
  TH1F *h_ElectronTrkIsoHeep_barrel_heep = new TH1F ("h_ElectronTrkIsoHeep_barrel_heep","h_ElectronTrkIsoHeep_barrel_heep",100,0,50);
  
  h_ElectronDeltaEtaTrkSC_barrel_heep->Sumw2();
  h_ElectronDeltaPhiTrkSC_barrel_heep->Sumw2();
  h_ElectronHoE_barrel_heep->Sumw2();
  h_ElectronSigmaIEtaIEta_barrel_heep->Sumw2();
  h_ElectronE2x5OverE5x5_barrel_heep->Sumw2();
  h_ElectronE1x5OverE5x5_barrel_heep->Sumw2();
  h_ElectronEcalIsoHeep_barrel_heep->Sumw2();
  h_ElectronHcalIsoD1Heep_barrel_heep->Sumw2();
  h_ElectronHcalIsoD2Heep_barrel_heep->Sumw2();
  h_ElectronTrkIsoHeep_barrel_heep->Sumw2();

  TH1F *h_ElectronDeltaEtaTrkSC_endcap_heep = new TH1F ("h_ElectronDeltaEtaTrkSC_endcap_heep","h_ElectronDeltaEtaTrkSC_endcap_heep",200,-0.05,0.05);
  TH1F *h_ElectronDeltaPhiTrkSC_endcap_heep = new TH1F ("h_ElectronDeltaPhiTrkSC_endcap_heep","h_ElectronDeltaPhiTrkSC_endcap_heep",200,-0.5,0.5);
  TH1F *h_ElectronHoE_endcap_heep = new TH1F ("h_ElectronHoE_endcap_heep","h_ElectronHoE_endcap_heep",75,0,0.15);
  TH1F *h_ElectronSigmaIEtaIEta_endcap_heep = new TH1F ("h_ElectronSigmaIEtaIEta_endcap_heep","h_ElectronSigmaIEtaIEta_endcap_heep",100,0,0.1);
  TH1F *h_ElectronE2x5OverE5x5_endcap_heep = new TH1F ("h_ElectronE2x5OverE5x5_endcap_heep","h_ElectronE2x5OverE5x5_endcap_heep",100,0,1);
  TH1F *h_ElectronE1x5OverE5x5_endcap_heep = new TH1F ("h_ElectronE1x5OverE5x5_endcap_heep","h_ElectronE1x5OverE5x5_endcap_heep",100,0,1);
  TH1F *h_ElectronEcalIsoHeep_endcap_heep = new TH1F ("h_ElectronEcalIsoHeep_endcap_heep","h2_ElectronEcalIsoHeep_endcap_heep",500,0,100);
  TH1F *h_ElectronHcalIsoD1Heep_endcap_heep = new TH1F ("h_ElectronHcalIsoD1Heep_endcap_heep","h_ElectronHcalIsoD1Heep_endcap_heep",500,0,100);
  TH1F *h_ElectronHcalIsoD2Heep_endcap_heep = new TH1F ("h_ElectronHcalIsoD2Heep_endcap_heep","h_ElectronHcalIsoD2Heep_endcap_heep",200,0,100);
  TH1F *h_ElectronTrkIsoHeep_endcap_heep = new TH1F ("h_ElectronTrkIsoHeep_endcap_heep","h_ElectronTrkIsoHeep_endcap_heep",200,0,100);
  
  h_ElectronDeltaEtaTrkSC_endcap_heep->Sumw2();
  h_ElectronDeltaPhiTrkSC_endcap_heep->Sumw2();
  h_ElectronHoE_endcap_heep->Sumw2();
  h_ElectronSigmaIEtaIEta_endcap_heep->Sumw2();
  h_ElectronE2x5OverE5x5_endcap_heep->Sumw2();
  h_ElectronE1x5OverE5x5_endcap_heep->Sumw2();
  h_ElectronEcalIsoHeep_endcap_heep->Sumw2();
  h_ElectronHcalIsoD1Heep_endcap_heep->Sumw2();
  h_ElectronHcalIsoD2Heep_endcap_heep->Sumw2();
  h_ElectronTrkIsoHeep_endcap_heep->Sumw2();

  
  TH1F *h_CaloJetResEMF_looseID = new TH1F ("h_CaloJetResEMF_looseID","h_CaloJetResEMF_looseID",100,-2,2);
  TH1F *h_CaloJetFHPD_looseID = new TH1F ("h_CaloJetFHPD_looseID","h_CaloJetFHPD_looseID",100,-2,2);
  TH1F *h_CaloJetN90Hits_looseID = new TH1F ("h_CaloJetN90Hits_looseID","h_CaloJetN90Hits_looseID",50,0,100);

  h_CaloJetResEMF_looseID->Sumw2();
  h_CaloJetFHPD_looseID->Sumw2();
  h_CaloJetN90Hits_looseID->Sumw2();

  ////////////////////// User's code to book histos - END ///////////////////////

  ////////////////////// User's code to get preCut values - BEGIN ///////////////

  double eleEta_bar = getPreCutValue1("eleEta_bar");
  double eleEta_end_min = getPreCutValue1("eleEta_end");
  double eleEta_end_max = getPreCutValue2("eleEta_end");

  double EleEnergyScale_EB=getPreCutValue1("EleEnergyScale_EB");
  double EleEnergyScale_EE=getPreCutValue1("EleEnergyScale_EE");
  double JetEnergyScale=getPreCutValue1("JetEnergyScale");

  // Not used when using ElectronHeepID and heepBitMask // int eleIDType = (int) getPreCutValue1("eleIDType");
  int heepBitMask_EB  =  getPreCutValue1("heepBitMask_EBGapEE") ;
  int heepBitMask_GAP =  getPreCutValue2("heepBitMask_EBGapEE") ;

  int heepBitMask_EE  =  getPreCutValue3("heepBitMask_EBGapEE") ;

  int looseBitMask_EB       =  getPreCutValue1("looseBitMask_EBGapEE") ;
  int looseBitMask_GAP      =  getPreCutValue2("looseBitMask_EBGapEE") ;
  int looseBitMask_EE       =  getPreCutValue3("looseBitMask_EBGapEE") ;
  int looseBitMask_enabled  =  getPreCutValue4("looseBitMask_EBGapEE") ;

  ////////////////////// User's code to get preCut values - END /////////////////
    
  Long64_t nentries = fChain->GetEntriesFast();
  STDOUT("analysisClass::Loop(): nentries = " << nentries);   
  
  ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
  ////// If the root version is updated and rootNtupleClass regenerated,     /////
  ////// these lines may need to be updated.                                 /////    
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) { // Begin of loop over events
  //for (Long64_t jentry=0; jentry<100;jentry++) { // Begin of loop over events
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry < 10 || jentry%1000 == 0) STDOUT("analysisClass::Loop(): jentry = " << jentry);   
    // if (Cut(ientry) < 0) continue;
    
    ////////////////////// User's code to be done for every event - BEGIN ///////////////////////

    // EES and JES
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
      {
	for(int ijet=0; ijet<CaloJetPt->size(); ijet++)
	  {
	    CaloJetPt->at(ijet) *= JetEnergyScale;
	  }
      }

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

    // Superclusters
     vector<int> v_idx_sc_all;
     vector<int> v_idx_sc_PtCut;
     vector<int> v_idx_sc_Iso;

     ///Find two highest Pt SC, since they are not ordered 
     float scHighestPt = -99;
     float scNextPt = -99;
     int idx_scHighestPt = -1;
     int idx_scNextPt = -1;

     float scHighestPt_PtCut = -99;
     float scNextPt_PtCut = -99;
     int idx_scHighestPt_PtCut = -1;
     int idx_scNextPt_PtCut = -1;

     float scHighestPt_Iso = -99;
     float scNextPt_Iso = -99;
     int idx_scHighestPt_Iso = -1;
     int idx_scNextPt_Iso = -1;


     for(int isc=0;isc<SuperClusterPt->size();isc++){
      if ( 1 - SuperClusterS4S1->at(isc) > 0.95 ) continue; 
       bool PassPt = false;
       bool PassHoE = false;
      if ( SuperClusterPt->at(isc) > getPreCutValue1("ele_PtCut") ) PassPt=true;
      if (SuperClusterHoE->at(isc)<0.05) PassHoE=true;


      if (SuperClusterPt->at(isc)>scHighestPt){
	scNextPt = scHighestPt;
	idx_scNextPt = idx_scHighestPt;
	scHighestPt = SuperClusterPt->at(isc);
	idx_scHighestPt = isc;
      }
      else if (SuperClusterPt->at(isc)>scNextPt){
	scNextPt = SuperClusterPt->at(isc);
	idx_scNextPt = isc;
      }

      if (PassPt){
	if (SuperClusterPt->at(isc)>scHighestPt_PtCut){
	  scNextPt_PtCut = scHighestPt_PtCut;
	  idx_scNextPt_PtCut = idx_scHighestPt_PtCut;
	  scHighestPt_PtCut = SuperClusterPt->at(isc);
	  idx_scHighestPt_PtCut = isc;
	}
	else if (SuperClusterPt->at(isc)>scNextPt_PtCut){
	  scNextPt_PtCut = SuperClusterPt->at(isc);
	  idx_scNextPt_PtCut = isc;
	}
      }

      if (PassPt && PassHoE){
	if (SuperClusterPt->at(isc)>scHighestPt_Iso){
	  scNextPt_Iso = scHighestPt_Iso;
	  idx_scNextPt_Iso = idx_scHighestPt_Iso;
	  scHighestPt_Iso = SuperClusterPt->at(isc);
	  idx_scHighestPt_Iso = isc;
	}
	else if (SuperClusterPt->at(isc)>scNextPt_Iso){
	  scNextPt_Iso = SuperClusterPt->at(isc);
	  idx_scNextPt_Iso = isc;
	}
      }

     }
    if (idx_scHighestPt != -1) v_idx_sc_all.push_back(idx_scHighestPt);
    if (idx_scNextPt != -1)v_idx_sc_all.push_back(idx_scNextPt);
    if (idx_scHighestPt_PtCut != -1) v_idx_sc_PtCut.push_back(idx_scHighestPt_PtCut);
    if (idx_scNextPt_PtCut != -1)v_idx_sc_PtCut.push_back(idx_scNextPt_PtCut);
    if (idx_scHighestPt_Iso != -1) v_idx_sc_Iso.push_back(idx_scHighestPt_Iso);
    if (idx_scNextPt_Iso != -1)v_idx_sc_Iso.push_back(idx_scNextPt_Iso);

    //////now fill in the rest of the sc in whatever order they are
    for(int isc=0;isc<SuperClusterPt->size();isc++){
       if ( 1 - SuperClusterS4S1->at(isc) > 0.95 ) continue; 
       bool PassPt = false;
       bool PassHoE = false;
       if (isc!=idx_scHighestPt && isc!=idx_scNextPt) v_idx_sc_all.push_back(isc);
       if ( SuperClusterPt->at(isc) > getPreCutValue1("ele_PtCut") ) PassPt=true;
       if (SuperClusterHoE->at(isc)<0.05) PassHoE=true;

       if (isc!=idx_scHighestPt_PtCut && isc!=idx_scNextPt_PtCut && PassPt) v_idx_sc_PtCut.push_back(isc);
       if (isc!=idx_scHighestPt_Iso && isc!=idx_scNextPt_Iso && PassPt && PassHoE) v_idx_sc_Iso.push_back(isc);
    }

    // Electrons
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
	if( ElectronPt->at(iele) < getPreCutValue1("ele_PtCut") ) continue; 
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
	if ( (ElectronHeepID->at(iele) & ~heepBitMask)==0x0  && ElectronOverlaps->at(iele)==0 )
	  {
	    //STDOUT("ElectronHeepID = " << hex << ElectronHeepID->at(iele) << " ; ElectronPassID = " << ElectronPassID->at(iele) )
	    v_idx_ele_PtCut_IDISO_noOverlap.push_back(iele);
	  }

      } // End loop over electrons

    // tight-loose electrons, if enabled from cut file
    if ( looseBitMask_enabled == 1 && v_idx_ele_PtCut_IDISO_noOverlap.size() == 1 )
      {
	//STDOUT("v_idx_ele_PtCut_IDISO_noOverlap[0] = "<<v_idx_ele_PtCut_IDISO_noOverlap[0] << " - Pt = "<<ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]));
	bool loosePtLargerThanTightPt = true;
	for(int iele=0; iele<v_idx_ele_PtCut.size(); iele++)
	  {
	    
	    if (v_idx_ele_PtCut[iele] == v_idx_ele_PtCut_IDISO_noOverlap[0])
	      {
		loosePtLargerThanTightPt = false;
		continue;
	      }
	    // get looseBitMask for EB, GAP, EE 

	    int looseBitMask;
	    if( fabs(ElectronEta->at(v_idx_ele_PtCut[iele])) < eleEta_bar ) 
	      {
		looseBitMask = looseBitMask_EB;
	      }
	    else if ( fabs(ElectronEta->at(v_idx_ele_PtCut[iele])) > eleEta_end_min && fabs(ElectronEta->at(v_idx_ele_PtCut[iele])) < eleEta_end_max ) 
	      {
		looseBitMask = looseBitMask_EE;
	      }
	    else {
	      looseBitMask = looseBitMask_GAP;
	    }
	    if ( (ElectronHeepID->at(v_idx_ele_PtCut[iele]) & ~looseBitMask)==0x0  && ElectronOverlaps->at(v_idx_ele_PtCut[iele])==0 )
	      {
		if ( loosePtLargerThanTightPt )
		  {
		    v_idx_ele_PtCut_IDISO_noOverlap.insert(v_idx_ele_PtCut_IDISO_noOverlap.begin(),1,v_idx_ele_PtCut[iele]);
		  }
		else 
		  {
		    v_idx_ele_PtCut_IDISO_noOverlap.push_back(v_idx_ele_PtCut[iele]);
		  }
		break; // happy with one loose electron - Note: if you want more than 1 loose, pt sorting will not be OK with the code as is. 
	      }
	  }	
// 	for ( int i=0; i<v_idx_ele_PtCut_IDISO_noOverlap.size(); i++)
// 	  {
// 	    STDOUT("i="<<i<<", v_idx_ele_PtCut_IDISO_noOverlap[i] = "<<v_idx_ele_PtCut_IDISO_noOverlap[i] << ", Pt = "<<ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[i]));
// 	  }
      } // tight-loose electrons, if enabled from cut file

    // Jets
    vector<int> v_idx_jet_all;
    vector<int> v_idx_jet_PtCut;
    vector<int> v_idx_jet_PtCut_noOverlap;
    vector<int> v_idx_jet_PtCut_noOverlap_ID;
    vector<int> v_idx_jet_PtCut_noOverlap_ID_EtaCut;

    // Loop over jets
    for(int ijet=0; ijet<CaloJetPt->size(); ijet++)
      {
	//no cut on reco jets
	v_idx_jet_all.push_back(ijet);

	//pT pre-cut on reco jets
	if ( CaloJetPt->at(ijet) < getPreCutValue1("jet_PtCut") ) continue;
	v_idx_jet_PtCut.push_back(ijet);
      }

// 	//Checking overlap between electrons and jets
// 	int JetOverlapsWithEle = 0; //don't change the default (0) 
// 	float minDeltaR=9999.;
// 	TVector3 jet_vec;
// 	jet_vec.SetPtEtaPhi(CaloJetPt->at(ijet),CaloJetEta->at(ijet),CaloJetPhi->at(ijet));
// 	for (int i=0; i < v_idx_ele_PtCut_IDISO_noOverlap.size(); i++){
// 	  TVector3 ele_vec;	  
// 	  ele_vec.SetPtEtaPhi(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[i])
// 			      ,ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[i])
// 			      ,ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[i]));
// 	  double distance = jet_vec.DeltaR(ele_vec);
// 	  if (distance<minDeltaR) minDeltaR=distance;
// 	}
// 	if ( minDeltaR < getPreCutValue1("jet_ele_DeltaRcut") )
// 	  JetOverlapsWithEle = 1; //this jet overlaps with a good electron --> remove it from the analysis

// 	//pT pre-cut + no overlaps with electrons
// 	// ---- use the flag stored in rootTuples
// 	//if( ( CaloJetOverlaps->at(ijet) & 1 << eleIDType) == 0)/* NO overlap with electrons */  
// 	// && (caloJetOverlaps[ijet] & 1 << 5)==0 )/* NO overlap with muons */   
// 	// ----
// 	if( JetOverlapsWithEle == 0 )  /* NO overlap with electrons */  
// 	  v_idx_jet_PtCut_noOverlap.push_back(ijet);

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
            jet.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut[ijet]),
			     CaloJetEta->at(v_idx_jet_PtCut[ijet]),
			     CaloJetPhi->at(v_idx_jet_PtCut[ijet]),0);
	    double DR = jet.DeltaR(ele);
	    if (DR<minDR) 
	      {
		minDR = DR;
		ijet_minDR = ijet;
	      }
	  }
	if ( minDR < getPreCutValue1("jet_ele_DeltaRcut") && ijet_minDR > -1)
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
// 	    thisjet.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut[ijet]),
// 				 CaloJetEta->at(v_idx_jet_PtCut[ijet]),
// 				 CaloJetPhi->at(v_idx_jet_PtCut[ijet]),0);
// 	    STDOUT("CLEANING: j"<<ijet+1<<" Pt, eta, phi = " << ", "<<thisjet.Pt()<<", "<< thisjet.Eta() <<", "<< thisjet.Phi()<<" jetFlags="<<jetFlags[ijet] );
// 	  }
//       } // printouts for jet cleaning
    
    
    float JetEtaCutValue = getPreCutValue1("jet_EtaCut");
    for(int ijet=0; ijet<v_idx_jet_PtCut.size(); ijet++) //pT pre-cut + no overlaps with electrons + jetID
      {	
	bool passjetID = JetIdloose(CaloJetresEMF->at(v_idx_jet_PtCut[ijet]),CaloJetfHPD->at(v_idx_jet_PtCut[ijet]),CaloJetn90Hits->at(v_idx_jet_PtCut[ijet]), CaloJetEta->at(v_idx_jet_PtCut[ijet]));
	// ---- use the flag stored in rootTuples
	//if( (CaloJetOverlaps->at(v_idx_jet_PtCut[ijet]) & 1 << eleIDType) == 0  /* NO overlap with electrons */  
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
	    && fabs( CaloJetEta->at(v_idx_jet_PtCut[ijet]) ) < JetEtaCutValue )
	  // && (caloJetOverlaps[ijet] & 1 << 5)==0 )         /* NO overlap with muons */      
	  v_idx_jet_PtCut_noOverlap_ID_EtaCut.push_back(v_idx_jet_PtCut[ijet]);

	
	//NOTE: We should verify that caloJetOverlaps match with the code above
      } // End loop over jets
    


    //## Muons
    vector<int> v_idx_muon_all;
    vector<int> v_idx_muon_PtCut;
    vector<int> v_idx_muon_PtCut_IDISO;

    // Loop over muons
    for(int imuon=0; imuon<MuonPt->size(); imuon++)
      {
	//no cut on reco muons
	v_idx_muon_all.push_back(imuon);

	//pT pre-cut on ele
	if ( MuonPt->at(imuon) < getPreCutValue1("muon_PtCut") ) continue;
	 
	v_idx_muon_PtCut.push_back(imuon);

	//ID + ISO
        if ( MuonPassIso->at(imuon)==1 && MuonPassID->at(imuon)==1
             && fabs(MuonEta->at(imuon)) < getPreCutValue1("muFidRegion")
             && MuonTrkHits->at(imuon) >= getPreCutValue1("muNHits_minThresh")
             && fabs(MuonTrkD0->at(imuon)) < getPreCutValue1("muTrkD0Maximum") )
          v_idx_muon_PtCut_IDISO.push_back(imuon);

      } // End loop over muons


    /////  Define the fake rate for QCD and calculate prob. for each sc to be an ele
    double BarrelCross = 0.1;
    double BarrelSlope = 0;
    double EndcapCross = 0.1;
    double EndcapSlope = 0;

    double p1_all = 0, p2_all = 0, p3_all = 0;
    if (v_idx_sc_all.size()>=1){
      if (fabs(SuperClusterEta->at(v_idx_sc_all[0]))<eleEta_bar) p1_all = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_all[0]);
      if (fabs(SuperClusterEta->at(v_idx_sc_all[0]))>eleEta_end_min) p1_all = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_all[0]) ;
    }
    if (v_idx_sc_all.size()>=2){
      if (fabs(SuperClusterEta->at(v_idx_sc_all[1]))<eleEta_bar) p2_all = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_all[1]);
      if (fabs(SuperClusterEta->at(v_idx_sc_all[1]))>eleEta_end_min) p2_all = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_all[1]);
    }
    if (v_idx_sc_all.size()>=3){
      if (fabs(SuperClusterEta->at(v_idx_sc_all[2]))<eleEta_bar) p3_all = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_all[2]);
      if (fabs(SuperClusterEta->at(v_idx_sc_all[2]))>eleEta_end_min) p3_all = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_all[2]);
    }

    double p1_PtCut = 0, p2_PtCut = 0, p3_PtCut = 0;
    if (v_idx_sc_PtCut.size()>=1){
      if (fabs(SuperClusterEta->at(v_idx_sc_PtCut[0]))<eleEta_bar) p1_PtCut = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_PtCut[0]);
      if (fabs(SuperClusterEta->at(v_idx_sc_PtCut[0]))>eleEta_end_min) p1_PtCut = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_PtCut[0]) ;
    }
    if (v_idx_sc_PtCut.size()>=2){
      if (fabs(SuperClusterEta->at(v_idx_sc_PtCut[1]))<eleEta_bar) p2_PtCut = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_PtCut[1]);
      if (fabs(SuperClusterEta->at(v_idx_sc_PtCut[1]))>eleEta_end_min) p2_PtCut = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_PtCut[1]);
    }
    if (v_idx_sc_PtCut.size()>=3){
      if (fabs(SuperClusterEta->at(v_idx_sc_PtCut[2]))<eleEta_bar) p3_PtCut = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_PtCut[2]);
      if (fabs(SuperClusterEta->at(v_idx_sc_PtCut[2]))>eleEta_end_min) p3_PtCut = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_PtCut[2]);
    }

    double p1_Iso = 0, p2_Iso = 0, p3_Iso = 0;
    if (v_idx_sc_Iso.size()>=1){
      if (fabs(SuperClusterEta->at(v_idx_sc_Iso[0]))<eleEta_bar) p1_Iso = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_Iso[0]);
      if (fabs(SuperClusterEta->at(v_idx_sc_Iso[0]))>eleEta_end_min) p1_Iso = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_Iso[0]) ;
    }
    if (v_idx_sc_Iso.size()>=2){
      if (fabs(SuperClusterEta->at(v_idx_sc_Iso[1]))<eleEta_bar) p2_Iso = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_Iso[1]);
      if (fabs(SuperClusterEta->at(v_idx_sc_Iso[1]))>eleEta_end_min) p2_Iso = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_Iso[1]);
    }
    if (v_idx_sc_Iso.size()>=3){
      if (fabs(SuperClusterEta->at(v_idx_sc_Iso[2]))<eleEta_bar) p3_Iso = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_Iso[2]);
      if (fabs(SuperClusterEta->at(v_idx_sc_Iso[2]))>eleEta_end_min) p3_Iso = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_Iso[2]);
    }

    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();
    

    // Set the value of the variableNames listed in the cutFile to their current value
    
    // Trigger (L1 and HLT)
    if(isData==true)
      {
	fillVariableWithValue( "PassBPTX0", isBPTX0 ) ;
	fillVariableWithValue( "PassPhysDecl", isPhysDeclared ) ;       
      }
    else
      {
	fillVariableWithValue( "PassBPTX0", true ) ;
	fillVariableWithValue( "PassPhysDecl", true ) ;       
      }

    fillVariableWithValue( "PassHLT", PassTrig ) ;

    //Event filters at RECO level
    fillVariableWithValue( "PassBeamScraping", !isBeamScraping ) ;
    fillVariableWithValue( "PassPrimaryVertex", isPrimaryVertex ) ;
    fillVariableWithValue( "PassHBHENoiseFilter", passLooseNoiseFilter ) ;

    // nEle
    fillVariableWithValue( "nEle_all", v_idx_sc_all.size() ) ;
    fillVariableWithValue( "nEle_PtCut", v_idx_sc_PtCut.size() ) ;
    fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp", v_idx_sc_Iso.size()) ;

//     fillVariableWithValue( "nEle_all", v_idx_ele_all.size() ) ;
//     fillVariableWithValue( "nEle_PtCut", v_idx_ele_PtCut.size() ) ;
//     fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp", v_idx_ele_PtCut_IDISO_noOverlap.size()) ;
     
    // nJet
    fillVariableWithValue( "nJet_all", v_idx_jet_all.size() ) ;
    fillVariableWithValue( "nJet_PtCut", v_idx_jet_PtCut.size() ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlp", v_idx_jet_PtCut_noOverlap.size() ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlp_ID", v_idx_jet_PtCut_noOverlap_ID.size() ) ;
    //TwoEleOnly
    fillVariableWithValue( "nJet_TwoEleOnly_All", v_idx_jet_PtCut_noOverlap_ID.size() ) ;
    fillVariableWithValue( "nJet_TwoEleOnly_EtaCut", v_idx_jet_PtCut_noOverlap_ID_EtaCut.size() ) ;
    //PAS June 2010
    fillVariableWithValue( "nJet_PAS_All", v_idx_jet_PtCut_noOverlap_ID.size() ) ;
    fillVariableWithValue( "nJet_PAS_EtaCut", v_idx_jet_PtCut_noOverlap_ID_EtaCut.size() ) ;

    // MET
    //PAS June 2010
    fillVariableWithValue( "pfMET_PAS", PFMET->at(0) ) ;
    fillVariableWithValue( "tcMET_PAS", TCMET->at(0) ) ;
    fillVariableWithValue( "caloMET_PAS", CaloMET->at(0) ) ;

    // 1st ele
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stEle_IDISO_NoOvrlp", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	fillVariableWithValue( "Eta1stEle_IDISO_NoOvrlp", ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	fillVariableWithValue( "mEta1stEle_IDISO_NoOvrlp", fabs(ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0])) );
	//PAS June 2010
	fillVariableWithValue( "Pt1stEle_PAS", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	fillVariableWithValue( "Eta1stEle_PAS", ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
      }

    // 2nd ele
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 2 ) 
      {
	fillVariableWithValue( "Pt2ndEle_IDISO_NoOvrlp", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) );
	fillVariableWithValue( "Eta2ndEle_IDISO_NoOvrlp", ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) );
	fillVariableWithValue( "mEta2ndEle_IDISO_NoOvrlp", fabs(ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[1])) );
	fillVariableWithValue( "maxMEtaEles_IDISO_NoOvrl", max( getVariableValue("mEta1stEle_IDISO_NoOvrlp"), getVariableValue("mEta2ndEle_IDISO_NoOvrlp") ) );
	//PAS June 2010
	fillVariableWithValue( "Pt2ndEle_PAS", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) );
	fillVariableWithValue( "Eta2ndEle_PAS", ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) );
      }

    // 1st jet
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stJet_noOvrlp_ID", CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "Eta1stJet_noOvrlp_ID", CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "mEta1stJet_noOvrlp_ID", fabs(CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[0])) );
	//PAS June 2010
	fillVariableWithValue( "Pt1stJet_PAS", CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "Eta1stJet_PAS", CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
      }


    //cout << "2nd Jet" << endl;
    //## 2nd jet
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2 ) 
      {
	fillVariableWithValue( "Pt2ndJet_noOvrlp_ID", CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "Eta2ndJet_noOvrlp_ID", CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "mEta2ndJet_noOvrlp_ID", fabs(CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[1])) );
	fillVariableWithValue( "maxMEtaJets_noOvrlp_ID", max( getVariableValue("mEta1stJet_noOvrlp_ID"), getVariableValue("mEta2ndJet_noOvrlp_ID") ) );
	//PAS June 2010
	fillVariableWithValue( "Pt2ndJet_PAS", CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "Eta2ndJet_PAS", CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
      }

    //## define "2ele" and "2jets" booleans
    bool TwoEle=false;
    bool TwoJets=false;
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 2 ) TwoEle = true;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2 ) TwoJets = true;

    // ST
    double calc_sT=-999.; 
    if ( (TwoEle) && (TwoJets) ) 
      {
	calc_sT = 
	  ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) +
	  ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) +
	  CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) +
	  CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]);
	fillVariableWithValue("sT", calc_sT);
	fillVariableWithValue("sT_MLQ200", calc_sT);
	fillVariableWithValue("sT_MLQ250", calc_sT);
	fillVariableWithValue("sT_MLQ300", calc_sT);
	fillVariableWithValue("sT_MLQ400", calc_sT);       
	//PAS June 2010
	fillVariableWithValue("sT_PAS", calc_sT);
      }

    // ST jets
    if (TwoJets)
      {
	double calc_sTjet=-999.;
	calc_sTjet =
	  CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) + 
	  CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) ;
	fillVariableWithValue("sTjet_PAS", calc_sTjet);
      }

    // Mjj
    if (TwoJets)
      {
	TLorentzVector jet1, jet2, jj;
	jet1.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),0);
	jet2.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]),0);
	jj = jet1+jet2;
	//PAS June 2010
	fillVariableWithValue("Mjj_PAS", jj.M());
      }


    // Mee
    if (TwoEle)
      {
	TLorentzVector ele1, ele2, ee;
	ele1.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),0);
	ele2.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[1]),
			  ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[1]),
			  ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[1]),0);
	ee = ele1+ele2;
	fillVariableWithValue("Mee", ee.M());
	//TwoEleOnly
	fillVariableWithValue("Mee_TwoEleOnly", ee.M());
	fillVariableWithValue("Mee_presel", ee.M());
	//
	//PAS June 2010
	fillVariableWithValue("Mee_PAS", ee.M());

	double calc_sTele=-999.;
        calc_sTele =
	   ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) +
	   ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) ;
	fillVariableWithValue("sTele_PAS", calc_sTele);

	if(isData==true) 
	  {
	    STDOUT("Two electrons: Run, LS, Event = "<<run<<", "<<ls<<", "<<event);
	    STDOUT("Two electrons: M_ee, Pt_ee, Eta_ee, Phi_ee = "<<ee.M() <<", "<< ee.Pt() <<", "<< ee.Eta() <<", "<< ee.Phi());
	    STDOUT("Two electrons: 1st ele Pt, eta, phi = "<< ele1.Pt() <<", "<< ele1.Eta() <<", "<< ele1.Phi() );
	    STDOUT("Two electrons: 2nd ele Pt, eta, phi = "<< ele2.Pt() <<", "<< ele2.Eta() <<", "<< ele2.Phi() );
	  }
      }

    // Mej 
    double Me1j1, Me1j2, Me2j1, Me2j2 = -999;
    double deltaM_e1j1_e2j2 = 9999;
    double deltaM_e1j2_e2j1 = 9999;
    double Mej_1stPair = 0;
    double Mej_2ndPair = 0;
    double deltaR_e1j1 ;
    double deltaR_e2j2 ;
    double deltaR_e1j2 ;
    double deltaR_e2j1 ;
    if ( (TwoEle) && (TwoJets) ) // TwoEle and TwoJets
      {
	TLorentzVector jet1, jet2, ele1, ele2;
	ele1.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),0);
	ele2.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[1]),
			  ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[1]),
			  ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[1]),0);
	jet1.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),0);
	jet2.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]),0);
	TLorentzVector e1j1, e1j2, e2j1, e2j2;
	e1j1 = ele1 + jet1;
	e2j2 = ele2 + jet2;
	e2j1 = ele2 + jet1;
	e1j2 = ele1 + jet2;
	Me1j1 = e1j1.M();
	Me2j2 = e2j2.M();
	Me1j2 = e1j2.M();
	Me2j1 = e2j1.M();

	deltaM_e1j1_e2j2 = Me1j1 - Me2j2;
	deltaM_e1j2_e2j1 = Me1j2 - Me2j1;


	double deltaR_e1j1 = ele1.DeltaR(jet1);
	double deltaR_e2j2 = ele2.DeltaR(jet2);
	double deltaR_e1j2 = ele1.DeltaR(jet2);
	double deltaR_e2j1 = ele2.DeltaR(jet1);

// 	// Fill min DR between any of the 2 selected eles and any of the 2 selected jets
// 	double minDR_2ele_2jet = min ( min(deltaR_e1j1,deltaR_e2j2) , min(deltaR_e1j2,deltaR_e2j1) );
// 	fillVariableWithValue("minDR_2ele_2jet", minDR_2ele_2jet);

	if(fabs(deltaM_e1j1_e2j2) > fabs(deltaM_e1j2_e2j1))
	  {
	    Mej_1stPair = Me1j2;
	    Mej_2ndPair = Me2j1;
	    fillVariableWithValue("minDRej_selecPairs", min(deltaR_e1j2,deltaR_e2j1) );
	    fillVariableWithValue("minDRej_unselPairs", min(deltaR_e1j1,deltaR_e2j2) );
	  }
	else
	  {
	    Mej_1stPair = Me1j1;
	    Mej_2ndPair = Me2j2;
	    fillVariableWithValue("minDRej_selecPairs", min(deltaR_e1j1,deltaR_e2j2) );
	    fillVariableWithValue("minDRej_unselPairs", min(deltaR_e1j2,deltaR_e2j1) );
	  } 
	fillVariableWithValue("Mej_1stPair", Mej_1stPair);       
	fillVariableWithValue("Mej_2ndPair", Mej_2ndPair);
	//PAS June 2010
	h_Mej_PAS->Fill(Mej_1stPair);
	h_Mej_PAS->Fill(Mej_2ndPair);
	fillVariableWithValue("Mej_1stPair_PAS", Mej_1stPair);       
	fillVariableWithValue("Mej_2ndPair_PAS", Mej_2ndPair);

	// min and max DeltaR between electrons and any jet
	double minDeltaR_ej = 999999;
	double maxDeltaR_ej = -1;
	double thisMinDR, thisMaxDR, DR_thisjet_e1, DR_thisjet_e2;
	TLorentzVector thisjet;
	for(int ijet=0; ijet<v_idx_jet_PtCut_noOverlap_ID.size(); ijet++)
	  {
	    thisjet.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
				 CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
				 CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),0);
	    DR_thisjet_e1 = thisjet.DeltaR(ele1);
	    DR_thisjet_e2 = thisjet.DeltaR(ele2);
	    thisMinDR = min(DR_thisjet_e1, DR_thisjet_e2);
	    thisMaxDR = max(DR_thisjet_e1, DR_thisjet_e2);
	    if(thisMinDR < minDeltaR_ej)
	      minDeltaR_ej = thisMinDR;
	    if(thisMaxDR > maxDeltaR_ej)
	      maxDeltaR_ej = thisMaxDR;
	  } 
	fillVariableWithValue("minDeltaR_ej", minDeltaR_ej);
	fillVariableWithValue("maxDeltaR_ej", maxDeltaR_ej);

	// printouts for small Mej
	if(isData==true && ( Mej_1stPair<20 || Mej_2ndPair<20 ) ) // printouts for low Mej 
	  {
	    STDOUT("Mej<20GeV: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
	    STDOUT("Mej<20GeV: Mej_1stPair = "<<Mej_1stPair <<", Mej_2ndPair = "<< Mej_2ndPair );
	    STDOUT("Mej<20GeV: e1j1.M = "<<e1j1.M() <<", e2j2.M = "<<e2j2.M() <<", e1j2.M = "<<e1j2.M()  <<", e2j1.M = "<<e2j1.M()  );
	    STDOUT("Mej<20GeV: deltaM_e1j1_e2j2 = "<<deltaM_e1j1_e2j2 <<", deltaM_e1j2_e2j1 = "<<deltaM_e1j2_e2j1  );
	    STDOUT("Mej<20GeV: deltaR_e1j1 = "<<deltaR_e1j1 <<", deltaR_e2j2 = "<<deltaR_e2j2 <<", deltaR_e1j2 = "<<deltaR_e1j2  <<", deltaR_e2j1 = "<<deltaR_e2j1  );
// 	    STDOUT("Mej<20GeV: 1st ele Pt, eta, phi = "<< ele1.Pt() <<",\t"<< ele1.Eta() <<",\t"<< ele1.Phi() );
// 	    STDOUT("Mej<20GeV: 2nd ele Pt, eta, phi = "<< ele2.Pt() <<",\t"<< ele2.Eta() <<",\t"<< ele2.Phi() );
// 	    STDOUT("Mej<20GeV: 1st jet Pt, eta, phi = "<< jet1.Pt() <<",\t"<< jet1.Eta() <<",\t"<< jet1.Phi() );
// 	    STDOUT("Mej<20GeV: 2nd jet Pt, eta, phi = "<< jet2.Pt() <<",\t"<< jet2.Eta() <<",\t"<< jet2.Phi() );
	    TLorentzVector thisele;
	    for(int iele=0; iele<v_idx_ele_PtCut_IDISO_noOverlap.size(); iele++)
	      {
		thisele.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
				     ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
				     ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),0);
		STDOUT("Mej<20GeV: e"<<iele+1<<" Pt, eta, phi = " 
		       << ", "<<thisele.Pt()<<", "<< thisele.Eta() <<", "<< thisele.Phi()<<"; DR_j1, DR_j2 = "<< thisele.DeltaR(jet1)<<", "<<thisele.DeltaR(jet2));
	      }
	    TLorentzVector thisjet;
	    TLorentzVector thisjet_e1, thisjet_e2;
	    for(int ijet=0; ijet<v_idx_jet_PtCut_noOverlap_ID.size(); ijet++)
	      {
		thisjet.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
				     CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
				     CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),0);
		thisjet_e1 = thisjet + ele1;
		thisjet_e2 = thisjet + ele2;
		STDOUT("Mej<20GeV: j"<<ijet+1<<" Pt, eta, phi = " << ", "<<thisjet.Pt()<<", "<< thisjet.Eta() <<", "<< thisjet.Phi()<<"; DR_e1, DR_e2 = "<< thisjet.DeltaR(ele1)<<", "<<thisjet.DeltaR(ele2) << "; M_e1, M_e2 = " <<thisjet_e1.M() <<", "<<thisjet_e2.M() );
	      }
	  } // printouts for low Mej 

      } // TwoEle and TwoJets



    // Evaluate cuts (but do not apply them)
    evaluateCuts();

    // Fill histograms and do analysis based on cut evaluation
    //h_nEleFinal->Fill( ElectronPt->size() );
     

    //Large MET events passing pre-selection
    float METcut = 60;
    if( variableIsFilled("pfMET_PAS") && passedAllPreviousCuts("pfMET_PAS") && isData==true) 
      {

	if( getVariableValue("pfMET_PAS") > METcut )
	  {

	    STDOUT("pfMET>60GeV: ----------- START ------------");

	    STDOUT("pfMET>60GeV: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);	    
	    if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") )	      
	      STDOUT("pfMET>60GeV: Pt1stEle_PAS,Eta1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS"));
	    if( variableIsFilled("Pt2ndEle_PAS") && variableIsFilled("Eta2ndEle_PAS") )	      
	      STDOUT("pfMET>60GeV: Pt2ndEle_PAS,Eta2ndEle_PAS = "<<getVariableValue("Pt2ndEle_PAS")<<",\t"<<getVariableValue("Eta2ndEle_PAS"));
	    if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") )	      
	      STDOUT("pfMET>60GeV: Pt1stJet_PAS,Eta1stJet_PAS = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS"));
	    if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") )	      
	      STDOUT("pfMET>60GeV: Pt2ndJet_PAS,Eta2ndJet_PAS = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS"));
	    if( variableIsFilled("Mee_PAS") && variableIsFilled("Mjj_PAS") )	      
	      STDOUT("pfMET>60GeV: Mee_PAS,Mjj_PAS = "<<getVariableValue("Mee_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
	    if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )	      
	      STDOUT("pfMET>60GeV: Mej_1stPair_PAS,Mej_2ndPair_PAS = "<<getVariableValue("Mej_1stPair_PAS")
		     <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
	    if( variableIsFilled("pfMET_PAS") && variableIsFilled("caloMET_PAS") )	      
	      STDOUT("pfMET>60GeV: pfMET_PAS,caloMET_PAS = "<<getVariableValue("pfMET_PAS")<<",\t"<<getVariableValue("caloMET_PAS"));

	    STDOUT("pfMET>60GeV: ------------ END -------------");

	  }
      }

    //plot of electron ID/Isolation variables
    if( passedAllPreviousCuts("Pt1stEle_PAS") && (TwoEle) && (TwoJets) )
      {

	for(int iele=0 ; iele<v_idx_ele_PtCut_IDISO_noOverlap.size() ; iele++ )
	  {

	    int isBarrel = 0;
	    int isEndcap = 0;
	    
	    if( fabs( ElectronEta->at(iele) ) < 1.442 ) 
	      isBarrel = 1;
	    
	    if( fabs( ElectronEta->at(iele) ) > 1.560 && fabs( ElectronEta->at(iele) ) < 2.5 ) 
	      isEndcap = 1;
	    
	    if( isBarrel )
	      {
		h_ElectronDeltaEtaTrkSC_barrel_heep->Fill( ElectronDeltaEtaTrkSC->at(iele) ) ;
		h_ElectronDeltaPhiTrkSC_barrel_heep->Fill( ElectronDeltaPhiTrkSC->at(iele) ) ;
		h_ElectronHoE_barrel_heep->Fill( ElectronHoE->at(iele) );
		h_ElectronSigmaIEtaIEta_barrel_heep->Fill( ElectronSigmaIEtaIEta->at(iele) );
		h_ElectronE2x5OverE5x5_barrel_heep->Fill( ElectronE2x5OverE5x5->at(iele) );
		h_ElectronE1x5OverE5x5_barrel_heep->Fill( ElectronE1x5OverE5x5->at(iele) );
		h_ElectronEcalIsoHeep_barrel_heep->Fill( ElectronEcalIsoHeep->at(iele) );
		h_ElectronHcalIsoD1Heep_barrel_heep->Fill( ElectronHcalIsoD1Heep->at(iele) );
		h_ElectronHcalIsoD2Heep_barrel_heep->Fill( ElectronHcalIsoD2Heep->at(iele) );
		h_ElectronTrkIsoHeep_barrel_heep->Fill( ElectronTrkIsoHeep->at(iele) );
	      }

	    if( isEndcap )
	      {
		h_ElectronDeltaEtaTrkSC_endcap_heep->Fill( ElectronDeltaEtaTrkSC->at(iele) ) ;
		h_ElectronDeltaPhiTrkSC_endcap_heep->Fill( ElectronDeltaPhiTrkSC->at(iele) ) ;
		h_ElectronHoE_endcap_heep->Fill( ElectronHoE->at(iele) );
		h_ElectronSigmaIEtaIEta_endcap_heep->Fill( ElectronSigmaIEtaIEta->at(iele) );
		h_ElectronE2x5OverE5x5_endcap_heep->Fill( ElectronE2x5OverE5x5->at(iele) );
		h_ElectronE1x5OverE5x5_endcap_heep->Fill( ElectronE1x5OverE5x5->at(iele) );
		h_ElectronEcalIsoHeep_endcap_heep->Fill( ElectronEcalIsoHeep->at(iele) );
		h_ElectronHcalIsoD1Heep_endcap_heep->Fill( ElectronHcalIsoD1Heep->at(iele) );
		h_ElectronHcalIsoD2Heep_endcap_heep->Fill( ElectronHcalIsoD2Heep->at(iele) );
		h_ElectronTrkIsoHeep_endcap_heep->Fill( ElectronTrkIsoHeep->at(iele) );
	      }

	  }
      }
    
    if( passedAllPreviousCuts("Pt1stJet_PAS") && (TwoEle) && (TwoJets) )
      {
	//jet1
	h_CaloJetResEMF_looseID->Fill( CaloJetresEMF->at(v_idx_jet_PtCut_noOverlap_ID[0]) ) ;
	h_CaloJetFHPD_looseID->Fill( CaloJetfHPD->at(v_idx_jet_PtCut_noOverlap_ID[0]) ) ;
	h_CaloJetN90Hits_looseID->Fill( CaloJetn90Hits->at(v_idx_jet_PtCut_noOverlap_ID[0]) ) ;
    
	//jet2
	h_CaloJetResEMF_looseID->Fill( CaloJetresEMF->at(v_idx_jet_PtCut_noOverlap_ID[1]) ) ;
	h_CaloJetFHPD_looseID->Fill( CaloJetfHPD->at(v_idx_jet_PtCut_noOverlap_ID[1]) ) ;
	h_CaloJetN90Hits_looseID->Fill( CaloJetn90Hits->at(v_idx_jet_PtCut_noOverlap_ID[1]) ) ;
      }

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

  h_Mej_PAS->Write();

  h_ElectronDeltaEtaTrkSC_barrel_heep->Write();
  h_ElectronDeltaPhiTrkSC_barrel_heep->Write();
  h_ElectronHoE_barrel_heep->Write();
  h_ElectronSigmaIEtaIEta_barrel_heep->Write();
  h_ElectronE2x5OverE5x5_barrel_heep->Write();
  h_ElectronE1x5OverE5x5_barrel_heep->Write();
  h_ElectronEcalIsoHeep_barrel_heep->Write();
  h_ElectronHcalIsoD1Heep_barrel_heep->Write();
  h_ElectronHcalIsoD2Heep_barrel_heep->Write();
  h_ElectronTrkIsoHeep_barrel_heep->Write();

  h_ElectronDeltaEtaTrkSC_endcap_heep->Write();
  h_ElectronDeltaPhiTrkSC_endcap_heep->Write();
  h_ElectronHoE_endcap_heep->Write();
  h_ElectronSigmaIEtaIEta_endcap_heep->Write();
  h_ElectronE2x5OverE5x5_endcap_heep->Write();
  h_ElectronE1x5OverE5x5_endcap_heep->Write();
  h_ElectronEcalIsoHeep_endcap_heep->Write();
  h_ElectronHcalIsoD1Heep_endcap_heep->Write();
  h_ElectronHcalIsoD2Heep_endcap_heep->Write();
  h_ElectronTrkIsoHeep_endcap_heep->Write();

  h_CaloJetResEMF_looseID->Write();
  h_CaloJetFHPD_looseID->Write();
  h_CaloJetN90Hits_looseID->Write();


  ////////////////////// User's code to write histos - END ///////////////////////
  
  
  //STDOUT("analysisClass::Loop() ends");   
}
