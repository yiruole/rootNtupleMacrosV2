/*
Description and instructions:
This code creates vectors of HEEP electrons, isolated superclusters and cleaned jets (no overlaps) based on their index in the original collection.
It makes plots necessary for the calculation of the QCD fake rate, i.e. the pt/eta distributions of the HEEP electrons for the numerator and the pt/eta distributions of the isolated sc for the denominator.
A separate pyroot macro divides the two output histograms to get the fake rate.
There are plots for the pT of the leading electron and the ST in events passing the full selection criteria (ele criteria->sc criteria).
These plots have two versions, "actual" are events with esjj/eejj and "prob" are the weighted histograms for ssjj events based on the probability of one/both sc to be reco'ed as electrons.
The fake rate used for these histograms is taken from the cut file: cutTable_SCFakeRate.
Typically I run the code once to establish the fake rate, the re-run it with the correct rate in the cut file to accurately fill the "prob" histograms.
*/


#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

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

   TH1F *h_TrigDiff = new TH1F("TrigDiff","TrigDiff",3.0,-1.5,1.5);
   TH2F *h2_DebugTrig = new TH2F("DebugTrig","DebugTrig;HLTResults;HLTBits",2,0,2,2,0,2);

   TH1F *h_PtHat = new TH1F("","",400,0,2000);h_PtHat->Sumw2();
   
   TH1F *h_dPhi_JetSC = new TH1F("dPhi_JetSC","dPhi_JetSC",320,0,3.2); h_dPhi_JetSC->Sumw2();
   TH1F *h_dR_JetSC = new TH1F("dR_JetSC","dR_JetSC",600,0,3.0); h_dR_JetSC->Sumw2();
   TH1F *h_dR_SCSC = new TH1F("dR_SCSC","dR_SCSC",600,0,3.0); h_dR_SCSC->Sumw2();
   TH1F *h_dR_JetSC_2Jets = new TH1F("dR_JetSC_2Jets","dR_JetSC_2Jets",600,0,3.0); h_dR_JetSC_2Jets->Sumw2();
   TH2F *h_dR_JetSC_EcalIso = new TH2F("dR_JetSC_EcalIso","dR_JetSC_EcalIso",600,0,3.0,500,0,10); h_dR_JetSC->Sumw2();
   TH1F *h_NisoSC = new TH1F ("NisoSC","NisoSC",6,-0.5,5.5);  h_NisoSC->Sumw2();

   TH1F *h_goodEleSCPt = new TH1F ("goodEleSCPt","goodEleSCPt",1000,0,1000); h_goodEleSCPt->Sumw2();
   TH1F *h_goodEleSCEta = new TH1F ("goodEleSCEta","goodEleSCEta",100,-3.,3.); h_goodEleSCEta->Sumw2();
   TH1F *h_goodEleSCPt_2SC2Jets = new TH1F ("goodEleSCPt_2SC2Jets","goodEleSCPt_2SC2Jets",1000,0,1000); h_goodEleSCPt_2SC2Jets->Sumw2();
   TH1F *h_goodEleSCEta_2SC2Jets = new TH1F ("goodEleSCEta_2SC2Jets","goodEleSCEta_2SC2Jets",100,-3.,3.); h_goodEleSCEta_2SC2Jets->Sumw2();
   TH1F *h_goodEleSCPt_Barrel = new TH1F ("goodEleSCPt_Barrel","goodEleSCPt_Barrel",1000,0,1000); h_goodEleSCPt_Barrel->Sumw2();
   TH1F *h_goodEleSCPt_Barrel_2SC2Jets = new TH1F ("goodEleSCPt_Barrel_2SC2Jets","goodEleSCPt_Barrel_2SC2Jets",1000,0,1000); h_goodEleSCPt_Barrel_2SC2Jets->Sumw2();
   TH1F *h_goodEleSCPt_Endcap = new TH1F ("goodEleSCPt_Endcap","goodEleSCPt_Endcap",1000,0,1000); h_goodEleSCPt_Endcap->Sumw2();
   TH1F *h_goodEleSCPt_Endcap_2SC2Jets = new TH1F ("goodEleSCPt_Endcap_2SC2Jets","goodEleSCPt_Endcap_2SC2Jets",1000,0,1000); h_goodEleSCPt_Endcap_2SC2Jets->Sumw2();
   TH1F *h_goodElePt = new TH1F ("goodElePt","goodElePt",1000,0,1000); h_goodElePt->Sumw2();

   TH1F *h_goodSCPt = new TH1F ("goodSCPt","goodSCPt",1000,0,1000); h_goodSCPt->Sumw2();
   TH1F *h_goodSCEta = new TH1F ("goodSCEta","goodSCEta",100,-3.,3.); h_goodSCEta->Sumw2();
   TH1F *h_goodSCPt_2SC2Jets = new TH1F ("goodSCPt_2SC2Jets","goodSCPt_2SC2Jets",1000,0,1000); h_goodSCPt_2SC2Jets->Sumw2();
   TH1F *h_goodSCEta_2SC2Jets = new TH1F ("goodSCEta_2SC2Jets","goodSCEta_2SC2Jets",100,-3.,3.); h_goodSCEta_2SC2Jets->Sumw2();
   TH1F *h_goodSCPt_Barrel = new TH1F ("goodSCPt_Barrel","goodSCPt_Barrel",1000,0,1000); h_goodSCPt_Barrel->Sumw2();
   TH1F *h_goodSCPt_Barrel_2SC2Jets = new TH1F ("goodSCPt_Barrel_2SC2Jets","goodSCPt_Barrel_2SC2Jets",1000,0,1000); h_goodSCPt_Barrel_2SC2Jets->Sumw2();
   TH1F *h_goodSCPt_Endcap = new TH1F ("goodSCPt_Endcap","goodSCPt_Endcap",1000,0,1000); h_goodSCPt_Endcap->Sumw2();
   TH1F *h_goodSCPt_Endcap_2SC2Jets = new TH1F ("goodSCPt_Endcap_2SC2Jets","goodSCPt_Endcap_2SC2Jets",1000,0,1000); h_goodSCPt_Endcap_2SC2Jets->Sumw2();

   TH2F *h_scEcalIso_Et = new TH2F ("scEcalIso_Et","scEcalIso_Et",100,0,1000,500,0,20); h_scEcalIso_Et->Sumw2();
   TH1F *h_scHoE = new TH1F ("scHoE","scHoE",100,0,1.0); h_scHoE->Sumw2();
   TH1F *h_scTrkIso = new TH1F ("scTrkIso","scTrkIso",100,0,50); h_scTrkIso->Sumw2();
   TH1F *h_scIetaIeta = new TH1F ("scIetaIeta","scIetaIeta",100,-0.1,0.1); h_scIetaIeta->Sumw2();

   TH1F *h_eta_failHLT = new TH1F("eta_failHLT","eta_failHLT",500,-3.0,3.0);
   TH1F *h_phi_failHLT = new TH1F("phi_failHLT","phi_failHLT",100,-3.5,3.5);

   TH1F *h_probPt1stSc_esjj = new TH1F ("probPt1stSc_esjj","probPt1stSc_esjj",1000,0,1000);  h_probPt1stSc_esjj->Sumw2();  //N events based on fake rate
   TH1F *h_actualPt1stSc_esjj = new TH1F ("actualPt1stSc_esjj","actualPt1stSc_esjj",1000,0,1000);  h_actualPt1stSc_esjj->Sumw2();  //N events with at least 1 HEEP ele
   TH1F *h_probSt_esjj = new TH1F ("probSt_esjj","probSt_esjj",200,0,2000);  h_probSt_esjj->Sumw2();  //N events based on fake rate
   TH1F *h_actualSt_esjj = new TH1F ("actualSt_esjj","actualSt_esjj",200,0,2000);  h_actualSt_esjj->Sumw2();  //N events with at least 1 HEEP ele

   TH1F *h_probPt1stSc_eejj = new TH1F ("probPt1stSc_eejj","probPt1stSc_eejj",1000,0,1000);  h_probPt1stSc_eejj->Sumw2();  //N events based on fake rate
   TH1F *h_actualPt1stSc_eejj = new TH1F ("actualPt1stSc_eejj","actualPt1stSc_eejj",1000,0,1000);  h_actualPt1stSc_eejj->Sumw2();  //N events with at least 1 HEEP ele
   TH1F *h_probSt_eejj = new TH1F ("probSt_eejj","probSt_eejj",200,0,2000);  h_probSt_eejj->Sumw2();  //N events based on fake rate
   TH1F *h_actualSt_eejj = new TH1F ("actualSt_eejj","actualSt_eejj",200,0,2000);  h_actualSt_eejj->Sumw2();  //N events with at least 1 HEEP ele

   TH1F *h_SCpT_ssjj = new TH1F ("SCpT_ssjj","SCpT_ssjj",1000,0,1000); h_SCpT_ssjj->Sumw2();
   TH1F *h_SCeta_ssjj = new TH1F ("SCeta_ssjj","SCeta_ssjj",100,-3,3); h_SCeta_ssjj->Sumw2();
   TH1F *h_STss_ssjj = new TH1F ("STss_ssjj","STss_ssjj",100,-3,3); h_STss_ssjj->Sumw2();
   TH1F *h_STjj_ssjj = new TH1F ("STjj_ssjj","STjj_ssjj",100,-3,3); h_STjj_ssjj->Sumw2();

   /////////initialize variables
   double FailRate = 0;
   double NFailHLT = 0;
   double HasSC = 0;

  
  ////////////////////// User's code to book histos - END ///////////////////////
    
  Long64_t nentries = fChain->GetEntriesFast();
  STDOUT("analysisClass::Loop(): nentries = " << nentries);   
  
  ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
  ////// If the root version is updated and rootNtupleClass regenerated,     /////
  ////// these lines may need to be updated.                                 /////    
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) { // Begin of loop over events
    //for (Long64_t jentry=0; jentry<10;jentry++) { // Begin of loop over events
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry < 10 || jentry%10000 == 0) STDOUT("analysisClass::Loop(): jentry = " << jentry);   
    // if (Cut(ientry) < 0) continue;
    
    ////////////////////// User's code to be done for every event - BEGIN ///////////////////////

    //if (PtHat<30) continue;
    h_PtHat->Fill(PtHat);

    //## HLT
    bool PassTrig=HLTResults->at(1); // results of HLTPhoton15 
    //bool PassTrig=HLTBits->at(71); // results of HLTPhoton15 
    int TrigDiff = HLTBits->at(71) - HLTResults->at(1);

    // Electrons
    vector<int> v_idx_ele_all;
    vector<int> v_idx_ele_PtCut;
    vector<int> v_idx_ele_PtCut_IDISO_noOverlap;
    int eleIDType = (int) getPreCutValue1("eleIDType");
     
    //Loop over electrons
    for(int iele=0; iele<ElectronPt->size(); iele++)
      {
	//no cut on reco electrons
	v_idx_ele_all.push_back(iele); 

	//pT pre-cut on ele
	if( ElectronPt->at(iele) < getPreCutValue1("ele_PtCut") ) continue; 
	v_idx_ele_PtCut.push_back(iele);

// 	//ID + ISO + NO overlap with good muons 
// 	int eleID = ElectronPassID->at(iele);
// 	if ( (eleID & 1<<eleIDType) > 0  && ElectronOverlaps->at(iele)==0 )
// 	  {
// 	    v_idx_ele_PtCut_IDISO_noOverlap.push_back(iele);
// 	  }

	bool HEEP = false;
	if (fabs(ElectronEta->at(iele))<1.442){
	if (fabs(ElectronDeltaEtaTrkSC->at(iele))<0.005){
	  if (fabs(ElectronDeltaPhiTrkSC->at(iele))<0.09){
	    if (ElectronHoE->at(iele)<0.05){
		if ((ElectronE2x5OverE5x5->at(iele)>0.94)||(ElectronE1x5OverE5x5->at(iele)>0.83)){
		  if ((ElectronEcalIsoHeep->at(iele) + ElectronHcalIsoD1Heep->at(iele))<(2+0.03*ElectronPt->at(iele))){
		    if (ElectronTrkIsoHeep->at(iele)<7.5){
		      HEEP = true;
		    }
		  }
		}
	      }
	    }
	  }
	}

	double eleEt = ElectronPt->at(iele);
	if ((fabs(ElectronEta->at(iele))>1.56) && (fabs(ElectronEta->at(iele))<2.5)){
	if (fabs(ElectronDeltaEtaTrkSC->at(iele))<0.007){
	  if (fabs(ElectronDeltaPhiTrkSC->at(iele))<0.09){
	    if (ElectronHoE->at(iele)<0.05){
	      if (ElectronSigmaIEtaIEta->at(iele)<0.03){
		if (((ElectronEcalIsoHeep->at(iele) + ElectronHcalIsoD1Heep->at(iele))<2.5) || 
		    ((eleEt>50)&&(ElectronEcalIsoHeep->at(iele) + ElectronHcalIsoD1Heep->at(iele))< ((0.03*(eleEt-50))+2.5))){
		    if (ElectronTrkIsoHeep->at(iele)<15){
		      if (ElectronHcalIsoD2Heep->at(iele)<0.5){
		      HEEP = true;
		      }
		    }
		}
	      }
	    }
	  }
	}
	}

	if ( HEEP  && ElectronOverlaps->at(iele)==0 )
	  {
	    v_idx_ele_PtCut_IDISO_noOverlap.push_back(iele);
	  }

      } // End loop over electrons

    
    //////// Fill SC histos

     vector<int> v_idx_sc_iso;
     ///Find two highest Pt SC, since they are not ordered 
     float scHighestPt = -99;
     float scNextPt = -99;
     int idx_scHighestPt = -1;
     int idx_scNextPt = -1;
     for(int isc=0;isc<SuperClusterPt->size();isc++){
      if ( SuperClusterPt->at(isc) < getPreCutValue1("ele_PtCut") ) continue;
      if (SuperClusterHoE->at(isc)>0.05) continue;
      bool Barrel = false;
      bool Endcap = false;
      if (fabs(SuperClusterEta->at(isc))<1.442) Barrel = true;
      if ((fabs(SuperClusterEta->at(isc))<2.5)&&(fabs(SuperClusterEta->at(isc))>1.560)) Endcap = true;
      if (!Barrel && !Endcap) continue;
      if ((1-SuperClusterS4S1->at(isc))>0.95) continue; // spike cleaning
      if (Endcap && SuperClusterSigmaIEtaIEta->at(isc)>0.03) continue;
      if (Barrel && SuperClusterHEEPEcalIso->at(isc)>(6+(0.01*SuperClusterPt->at(isc)))) continue;
      if (Endcap && SuperClusterPt->at(isc)<50 && SuperClusterHEEPEcalIso->at(isc)>(6+(0.01*SuperClusterPt->at(isc)))) continue;
      if (Endcap && SuperClusterPt->at(isc)>=50 && SuperClusterHEEPEcalIso->at(isc)>(6+(0.01*(SuperClusterPt->at(isc)-50)))) continue;
      if (Barrel && SuperClusterHEEPTrkIso->at(isc)>7.5) continue;;
      if (Endcap && SuperClusterHEEPTrkIso->at(isc)>15) continue;
      if (SuperClusterPt->at(isc)>scHighestPt){
	scNextPt = scHighestPt;
	idx_scNextPt = idx_scNextPt;
	scHighestPt = SuperClusterPt->at(isc);
	idx_scHighestPt = isc;
      }
      else if (SuperClusterPt->at(isc)>scNextPt){
	scNextPt = SuperClusterPt->at(isc);
	idx_scNextPt = isc;
      }
     }
     if (idx_scHighestPt != -1) v_idx_sc_iso.push_back(idx_scHighestPt);
     if (idx_scNextPt != -1)v_idx_sc_iso.push_back(idx_scNextPt);

    //////now fill in the rest of the sc in whatever order they are
    for(int isc=0;isc<SuperClusterPt->size();isc++){
      if (isc==idx_scHighestPt || isc==idx_scNextPt) continue;
      if ( SuperClusterPt->at(isc) < getPreCutValue1("ele_PtCut") ) continue;
      bool scPassHoE= false;
      bool scPassSigmaEE= false;
      bool scPassEcalIso= false;
      bool scPassTrkIso= false;
      bool Barrel = false;
      bool Endcap = false;
      if (fabs(SuperClusterEta->at(isc))<1.442) Barrel = true;
      if ((fabs(SuperClusterEta->at(isc))<2.5)&&(fabs(SuperClusterEta->at(isc))>1.560)) Endcap = true;
      if (!Barrel && !Endcap) continue;
      if ((1-SuperClusterS4S1->at(isc))>0.95) continue; // spike cleaning
      if (SuperClusterHoE->at(isc)<0.05) scPassHoE=true;
      if (Endcap && SuperClusterSigmaIEtaIEta->at(isc)<0.03) scPassSigmaEE=true;
      if (Barrel) scPassSigmaEE=true;
      if (Barrel && SuperClusterHEEPEcalIso->at(isc)<(6+(0.01*SuperClusterPt->at(isc)))) scPassEcalIso = true;
      if (Endcap && SuperClusterPt->at(isc)<50 && SuperClusterHEEPEcalIso->at(isc)< (6+(0.01*SuperClusterPt->at(isc)))) scPassEcalIso = true;
      if (Endcap && SuperClusterPt->at(isc)>=50 && SuperClusterHEEPEcalIso->at(isc)< (6+(0.01*(SuperClusterPt->at(isc)-50)))) scPassEcalIso = true;
      if (Barrel && SuperClusterHEEPTrkIso->at(isc)<7.5) scPassTrkIso = true;
      if (Endcap && SuperClusterHEEPTrkIso->at(isc)<15) scPassTrkIso = true;
      if (scPassHoE && scPassSigmaEE && scPassEcalIso && scPassTrkIso ){
      //if (scPassHoE){
	v_idx_sc_iso.push_back(isc);
      }
    }

    vector<int> v_idx_sc_iso_barrel;
    for(int isc=0; isc<v_idx_sc_iso.size(); isc++)
      {
	if (fabs(SuperClusterEta->at(v_idx_sc_iso[isc]))>1.45) continue;
      v_idx_sc_iso_barrel.push_back(v_idx_sc_iso[isc]);
      }
    h_NisoSC->Fill(v_idx_sc_iso_barrel.size());

    // Jets
    vector<int> v_idx_jet_all;
    vector<int> v_idx_jet_PtCut;
    vector<int> v_idx_jet_PtCut_noOverlapSC;

    // Loop over jets
    for(int ijet=0; ijet<CaloJetPt->size(); ijet++)
      {

	//no cut on reco jets
	v_idx_jet_all.push_back(ijet);

	//pT pre-cut on reco jets
	if ( CaloJetPt->at(ijet) < getPreCutValue1("jet_PtCut") ) continue;

	v_idx_jet_PtCut.push_back(ijet);

	 ///// take the iso SC out of the jet collection
	 TVector3 jet_vec;
	 jet_vec.SetPtEtaPhi(CaloJetPt->at(ijet),CaloJetEta->at(ijet),CaloJetPhi->at(ijet));
	 float minDRsc = 99;
	 int idx_nearest_sc = -1;
	 for(int isc=0;isc<v_idx_sc_iso.size();isc++)
	     {
	       TVector3 sc_vec;
	       sc_vec.SetPtEtaPhi(SuperClusterPt->at(v_idx_sc_iso[isc]),
				  SuperClusterEta->at(v_idx_sc_iso[isc]),
				  SuperClusterPhi->at(v_idx_sc_iso[isc]));
	       float DR = jet_vec.DeltaR(sc_vec);
	       if (DR<minDRsc) {
		 minDRsc = DR;
		 idx_nearest_sc = isc;
	       }
	     }

	//pT pre-cut + no overlaps with electrons or superclusters
	if( ( CaloJetOverlaps->at(ijet) & 1 << eleIDType) == 0 && minDRsc>0.5 )/* NO overlap with electrons */  
	  v_idx_jet_PtCut_noOverlapSC.push_back(ijet);

      } // End loop over jets

     //find dR between isoSC and jets

	 for (int ijet=0; ijet<v_idx_jet_PtCut_noOverlapSC.size(); ijet++){
	   if (ijet>1) break;  // only care about two leading jets
	   TVector3 jet_vec;
	   jet_vec.SetPtEtaPhi(CaloJetPt->at(v_idx_jet_PtCut_noOverlapSC[ijet]),
			       CaloJetEta->at(v_idx_jet_PtCut_noOverlapSC[ijet]),
			       CaloJetPhi->at(v_idx_jet_PtCut_noOverlapSC[ijet]));
	   float minDR = 99;
	   int idx_nearest_sc = -1;
	   for(int isc=0;isc<v_idx_sc_iso.size();isc++)
	     {
	       if (isc>1) break; //only care about two leading SC
	       TVector3 sc_vec;
	       sc_vec.SetPtEtaPhi(SuperClusterPt->at(v_idx_sc_iso[isc]),
				  SuperClusterEta->at(v_idx_sc_iso[isc]),
				  SuperClusterPhi->at(v_idx_sc_iso[isc]));
	       float DR = jet_vec.DeltaR(sc_vec);
	       if (DR<minDR) {
		 minDR = DR;
		 idx_nearest_sc = isc;
	       }
	     }
	   if (minDR!=99) h_dR_JetSC->Fill(minDR);
	   if ((minDR!=99)&&(v_idx_jet_PtCut_noOverlapSC.size()>1)) h_dR_JetSC_2Jets->Fill(minDR);
	   if (idx_nearest_sc !=-1) h_dR_JetSC_EcalIso->Fill(minDR,SuperClusterHEEPEcalIso->at(v_idx_sc_iso[idx_nearest_sc]));
	 }


    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();
    

    // Set the value of the variableNames listed in the cutFile to their current value
    
    // HLT
    fillVariableWithValue( "HLT", PassTrig ) ;

    //## SC
    fillVariableWithValue( "nIsoSC", v_idx_sc_iso.size() );
    //fillVariableWithValue( "nIsoSC",v_idx_sc_iso_barrel.size() );

    // nEle
    fillVariableWithValue( "nEle_all", v_idx_ele_all.size() ) ;
    fillVariableWithValue( "nEle_PtCut", v_idx_ele_PtCut.size() ) ;
    fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp", v_idx_ele_PtCut_IDISO_noOverlap.size() ) ;
     
    // nJet
    fillVariableWithValue( "nJet_all", v_idx_jet_all.size() ) ;
    fillVariableWithValue( "nJet_PtCut", v_idx_jet_PtCut.size() ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlpSC", v_idx_jet_PtCut_noOverlapSC.size() ) ;

    // 1st SC
    if( v_idx_sc_iso.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stSC_ISO", SuperClusterPt->at(v_idx_sc_iso[0]) );
	fillVariableWithValue( "Eta1stSC_ISO", SuperClusterEta->at(v_idx_sc_iso[0]) );
	fillVariableWithValue( "mEta1stSC_ISO", fabs(SuperClusterEta->at(v_idx_sc_iso[0])) );
      }

    // 2nd SC
    if( v_idx_sc_iso.size() >= 2 ) 
      {
	fillVariableWithValue( "Pt2ndSC_ISO", SuperClusterPt->at(v_idx_sc_iso[1]) );
	fillVariableWithValue( "Eta2ndSC_ISO", SuperClusterEta->at(v_idx_sc_iso[1]) );
	fillVariableWithValue( "mEta2ndSC_ISO", fabs(SuperClusterEta->at(v_idx_sc_iso[1])) );
	fillVariableWithValue( "maxMEtaSC_ISO", max( getVariableValue("mEta1stSC_ISO"), getVariableValue("mEta2ndSC_ISO") ) );
	fillVariableWithValue( "STee", SuperClusterPt->at(v_idx_sc_iso[0])+SuperClusterPt->at(v_idx_sc_iso[1]));
      }

    // 1st jet
    if( v_idx_jet_PtCut_noOverlapSC.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stJet_noOvrlpSC", CaloJetPt->at(v_idx_jet_PtCut_noOverlapSC[0]) );
	fillVariableWithValue( "Eta1stJet_noOvrlpSC", CaloJetEta->at(v_idx_jet_PtCut_noOverlapSC[0]) );
	fillVariableWithValue( "mEta1stJet_noOvrlpSC", fabs(CaloJetEta->at(v_idx_jet_PtCut_noOverlapSC[0])) );
      }


    //cout << "2nd Jet" << endl;
    //## 2nd jet
    if( v_idx_jet_PtCut_noOverlapSC.size() >= 2 ) 
      {
	fillVariableWithValue( "Pt2ndJet_noOvrlpSC", CaloJetPt->at(v_idx_jet_PtCut_noOverlapSC[1]) );
	fillVariableWithValue( "Eta2ndJet_noOvrlpSC", CaloJetEta->at(v_idx_jet_PtCut_noOverlapSC[1]) );
	fillVariableWithValue( "mEta2ndJet_noOvrlpSC", fabs(CaloJetEta->at(v_idx_jet_PtCut_noOverlapSC[1])) );
	fillVariableWithValue( "maxMEtaJets_noOvrlpSC", max( getVariableValue("mEta1stJet_noOvrlpSC"), getVariableValue("mEta2ndJet_noOvrlpSC") ) );
	fillVariableWithValue( "STjj", CaloJetPt->at(v_idx_jet_PtCut_noOverlapSC[0])+CaloJetPt->at(v_idx_jet_PtCut_noOverlapSC[1]));
      }

    //## define "2ele" and "2jets" booleans
    bool TwoEle=false;
    bool TwoSC=false;
    bool TwoJets=false;
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 2 ) TwoEle = true;
    if( v_idx_jet_PtCut_noOverlapSC.size() >= 2 ) TwoJets = true;
    if( v_idx_sc_iso.size()>= 2) TwoSC = true;

    // ST
    double calc_sT=-999.; 
    if ( (TwoSC) && (TwoJets) ) 
      {
	calc_sT = 
	  SuperClusterPt->at(v_idx_sc_iso[0]) +
	  SuperClusterPt->at(v_idx_sc_iso[1]) +
	  CaloJetPt->at(v_idx_jet_PtCut_noOverlapSC[0]) +
	  CaloJetPt->at(v_idx_jet_PtCut_noOverlapSC[1]);
	fillVariableWithValue("sT", calc_sT);
      }

    // Mee
    double MassEE=0;
    if (TwoSC)
      {
	TLorentzVector ele1, ele2, ee;
	ele1.SetPtEtaPhiM(SuperClusterPt->at(v_idx_sc_iso[0]),
			  SuperClusterEta->at(v_idx_sc_iso[0]),
			  SuperClusterPhi->at(v_idx_sc_iso[0]),0);
	ele2.SetPtEtaPhiM(SuperClusterPt->at(v_idx_sc_iso[1]),
			  SuperClusterEta->at(v_idx_sc_iso[1]),
			  SuperClusterPhi->at(v_idx_sc_iso[1]),0);
	ee = ele1+ele2;
	MassEE=ee.M();
	fillVariableWithValue("Mee", ee.M());
	h_dR_SCSC->Fill(ele1.DeltaR(ele2));
      }

    // Mej 
    double M11, M12, M21, M22 = -999;
    double diff_11_22 = 9999;
    double diff_12_21 = 9999;
    if ( (TwoSC) && (TwoJets) ) 
      {
	TLorentzVector jet1, jet2, ele1, ele2;
	ele1.SetPtEtaPhiM(SuperClusterPt->at(v_idx_sc_iso[0]),
			  SuperClusterEta->at(v_idx_sc_iso[0]),
			  SuperClusterPhi->at(v_idx_sc_iso[0]),0);
	ele2.SetPtEtaPhiM(SuperClusterPt->at(v_idx_sc_iso[1]),
			  SuperClusterEta->at(v_idx_sc_iso[1]),
			  SuperClusterPhi->at(v_idx_sc_iso[1]),0);
	jet1.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlapSC[0]),
			  CaloJetEta->at(v_idx_jet_PtCut_noOverlapSC[0]),
			  CaloJetPhi->at(v_idx_jet_PtCut_noOverlapSC[0]),0);
	jet2.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlapSC[1]),
			  CaloJetEta->at(v_idx_jet_PtCut_noOverlapSC[1]),
			  CaloJetPhi->at(v_idx_jet_PtCut_noOverlapSC[1]),0);
	TLorentzVector jet1ele1, jet2ele1, jet1ele2, jet2ele2;
	jet1ele1 = jet1 + ele1;
	jet2ele1 = jet2 + ele1;
	jet1ele2 = jet1 + ele2;
	jet2ele2 = jet2 + ele2;
	M11 = jet1ele1.M();
	M21 = jet2ele1.M();
	M12 = jet1ele2.M();
	M22 = jet2ele2.M();

	diff_11_22 = M11 - M22;
	diff_12_21 = M12 - M21;

	if(diff_11_22 > diff_12_21)
	  {
	    fillVariableWithValue("Mej", M12);       
	    fillVariableWithValue("Mej", M21);
	  }
	else
	  {
	    fillVariableWithValue("Mej", M22);       
	    fillVariableWithValue("Mej", M11);
	  } 
      }


    // Evaluate cuts (but do not apply them)
    evaluateCuts();

     if (passedCut("nIsoSC")){
       HasSC++;
     }

     // Fill histograms and do analysis based on cut evaluation

     //calculate predicted and actual N events with HEEP ele
     if( passedCut("0")&&passedCut("1")&&passedCut("2")&&passedCut("3") ) 
       {
	 for (int iele=0; iele<v_idx_ele_PtCut_IDISO_noOverlap.size(); iele++){
	   h_actualPt1stSc_esjj->Fill(ElectronSCPt->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]));
	   h_actualSt_esjj->Fill(calc_sT);
	   if (v_idx_ele_PtCut_IDISO_noOverlap.size()>=2){
	     h_actualPt1stSc_eejj->Fill(ElectronSCPt->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]));
	     h_actualSt_eejj->Fill(calc_sT);
	   }
	 }

	 double probSC1 = 0, probSC2 = 0;
	 double BarrelCross = getPreCutValue1("FakeRate_BarrelCross");
	 double BarrelSlope = getPreCutValue1("FakeRate_BarrelSlope");
	 double EndcapCross = getPreCutValue1("FakeRate_EndcapCross");
	 double EndcapSlope = getPreCutValue1("FakeRate_EndcapSlope");

	 if (fabs(SuperClusterEta->at(v_idx_sc_iso[0]))<1.442) probSC1 = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_iso[0]);
	 if (fabs(SuperClusterEta->at(v_idx_sc_iso[1]))<1.442) probSC2 = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_iso[1]);
	 if (fabs(SuperClusterEta->at(v_idx_sc_iso[0]))>1.56 && fabs(SuperClusterEta->at(v_idx_sc_iso[0]))<2.5 && (SuperClusterPt->at(v_idx_sc_iso[0])<70)) probSC1 = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_iso[0]) ;
	 if (fabs(SuperClusterEta->at(v_idx_sc_iso[1]))>1.56 && fabs(SuperClusterEta->at(v_idx_sc_iso[1]))<2.5 && (SuperClusterPt->at(v_idx_sc_iso[0])<70) ) probSC2 = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_iso[1]);
	 if (fabs(SuperClusterEta->at(v_idx_sc_iso[0]))>1.56 && fabs(SuperClusterEta->at(v_idx_sc_iso[0]))<2.5 && (SuperClusterPt->at(v_idx_sc_iso[0])>=70)) probSC1 = EndcapCross + EndcapSlope*70 ;
	 if (fabs(SuperClusterEta->at(v_idx_sc_iso[1]))>1.56 && fabs(SuperClusterEta->at(v_idx_sc_iso[1]))<2.5 && (SuperClusterPt->at(v_idx_sc_iso[0])>=70) ) probSC2 = EndcapCross + EndcapSlope*70;
      
	 h_probPt1stSc_esjj->Fill(SuperClusterPt->at(v_idx_sc_iso[0]),probSC1+probSC2);
	 h_probSt_esjj->Fill(calc_sT,probSC1+probSC2);
	 h_probPt1stSc_eejj->Fill(SuperClusterPt->at(v_idx_sc_iso[0]),probSC1*probSC2);
	 h_probSt_eejj->Fill(calc_sT,probSC1*probSC2);

	 h_STss_ssjj->Fill(SuperClusterPt->at(v_idx_sc_iso[0])+SuperClusterPt->at(v_idx_sc_iso[1]));
	 h_STjj_ssjj->Fill(CaloJetPt->at(v_idx_jet_PtCut_noOverlapSC[0])+CaloJetPt->at(v_idx_jet_PtCut_noOverlapSC[1]));
	 for(int isc=0;isc<v_idx_sc_iso.size();isc++)
	   {
	     h_SCpT_ssjj->Fill(SuperClusterPt->at(v_idx_sc_iso[isc]));
	     h_SCeta_ssjj->Fill(SuperClusterEta->at(v_idx_sc_iso[isc]));
	   }

       }

     //// Fill fake rate plots
     if (MassEE>60 && MassEE<120) continue; //get rid of Zs
     if (PFMET->at(0) > 10) continue; // get rid of Ws
	 for(int isc=0;isc<v_idx_sc_iso.size();isc++)
	   {
	     //Require dR>3 between SC and JET
	     TVector3 sc_vec;
	     sc_vec.SetPtEtaPhi(SuperClusterPt->at(v_idx_sc_iso[isc]),
			   SuperClusterEta->at(v_idx_sc_iso[isc]),
			   SuperClusterPhi->at(v_idx_sc_iso[isc]));
	     double dPhi_SC_Jet=0;
	     for (int ijet=0;ijet<v_idx_jet_PtCut.size();ijet++){
	       TVector3 jet_vec;
	       jet_vec.SetPtEtaPhi(CaloJetPt->at(v_idx_jet_PtCut[ijet]),
			    CaloJetEta->at(v_idx_jet_PtCut[ijet]),
			    CaloJetPhi->at(v_idx_jet_PtCut[ijet]));
	       double deltaPhi=fabs(jet_vec.DeltaPhi(sc_vec));
	       if (deltaPhi>dPhi_SC_Jet)dPhi_SC_Jet=deltaPhi;
	     }
	     h_dPhi_JetSC->Fill(dPhi_SC_Jet);	
	     if (dPhi_SC_Jet<3.0) continue;

	     h_goodSCPt->Fill(SuperClusterPt->at(v_idx_sc_iso[isc]));
	     h_goodSCEta->Fill(SuperClusterEta->at(v_idx_sc_iso[isc]));

	     bool Barrel = false;
	     bool Endcap = false;
	     if (fabs(SuperClusterEta->at(v_idx_sc_iso[isc]))<1.45) Barrel = true;
	     if (fabs(SuperClusterEta->at(v_idx_sc_iso[isc]))>1.56 && fabs(SuperClusterEta->at(v_idx_sc_iso[isc]))<2.5) Endcap = true;

	     if (Barrel) h_goodSCPt_Barrel->Fill(SuperClusterPt->at(v_idx_sc_iso[isc]));
	     if (Barrel && TwoSC && TwoJets) h_goodSCPt_Barrel_2SC2Jets->Fill(SuperClusterPt->at(v_idx_sc_iso[isc]));

	     if (Endcap) h_goodSCPt_Endcap->Fill(SuperClusterPt->at(v_idx_sc_iso[isc]));
	     if (Endcap && TwoSC && TwoJets) h_goodSCPt_Endcap_2SC2Jets->Fill(SuperClusterPt->at(v_idx_sc_iso[isc]));

	     if (TwoSC && TwoJets) {
	       h_goodSCPt_2SC2Jets->Fill(SuperClusterPt->at(v_idx_sc_iso[isc]));
	       h_goodSCEta_2SC2Jets->Fill(SuperClusterEta->at(v_idx_sc_iso[isc]));
	     }

	     // check if there is a HEEP ele associated with this sc
	     double deltaR_ele_sc = 99;
	     int idx_HEEP = -1;
	     for(int iele=0;iele<v_idx_ele_PtCut_IDISO_noOverlap.size();iele++){
	       TVector3 ele_vec;
	       ele_vec.SetPtEtaPhi(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
			   ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
			   ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]));
	       double tmp_deltaR = ele_vec.DeltaR(sc_vec);
	       if (tmp_deltaR<deltaR_ele_sc){
		 deltaR_ele_sc = tmp_deltaR;
		 idx_HEEP = iele;
	       }
	     }
	     if (deltaR_ele_sc<0.5){
		h_goodEleSCPt->Fill(ElectronSCPt->at(v_idx_ele_PtCut_IDISO_noOverlap[idx_HEEP]));
		h_goodEleSCEta->Fill(ElectronSCEta->at(v_idx_ele_PtCut_IDISO_noOverlap[idx_HEEP]));
		
		bool eleBarrel = false;
		bool eleEndcap = false;
		if (fabs(ElectronSCEta->at(v_idx_ele_PtCut_IDISO_noOverlap[idx_HEEP]))<1.445) eleBarrel = true;
		if (fabs(ElectronSCEta->at(v_idx_ele_PtCut_IDISO_noOverlap[idx_HEEP]))>1.56 && fabs(ElectronSCEta->at(v_idx_ele_PtCut_IDISO_noOverlap[idx_HEEP]))<2.5) eleEndcap = true;
		
		if (eleBarrel) h_goodEleSCPt_Barrel->Fill(ElectronSCPt->at(v_idx_ele_PtCut_IDISO_noOverlap[idx_HEEP]));
		if (eleBarrel && TwoSC && TwoJets) h_goodEleSCPt_Barrel_2SC2Jets->Fill(ElectronSCPt->at(v_idx_ele_PtCut_IDISO_noOverlap[idx_HEEP]));
		
		if (eleEndcap) h_goodEleSCPt_Endcap->Fill(ElectronSCPt->at(v_idx_ele_PtCut_IDISO_noOverlap[idx_HEEP]));
		if (eleEndcap && TwoSC && TwoJets) h_goodEleSCPt_Endcap_2SC2Jets->Fill(ElectronSCPt->at(v_idx_ele_PtCut_IDISO_noOverlap[idx_HEEP]));
		
		if (TwoSC && TwoJets) {
		  h_goodEleSCPt_2SC2Jets->Fill(ElectronSCPt->at(v_idx_ele_PtCut_IDISO_noOverlap[idx_HEEP]));
		  h_goodEleSCEta_2SC2Jets->Fill(ElectronSCEta->at(v_idx_ele_PtCut_IDISO_noOverlap[idx_HEEP]));
		}
	     }
	   } // for sc


    ////////////////////// User's code to be done for every event - END ///////////////////////
    
  } // End of loop over events
  

  ////////////////////// User's code to write histos - BEGIN ///////////////////////

   h_TrigDiff->Write();
   h2_DebugTrig->Write();

   h_dPhi_JetSC->Write();
   h_dR_JetSC->Write();
   h_dR_SCSC->Write();
   h_dR_JetSC_2Jets->Write();
   h_dR_JetSC_EcalIso->Write();
   h_NisoSC->Write();
   h_goodEleSCPt->Write();
   h_goodEleSCEta->Write();
   h_goodEleSCPt_2SC2Jets->Write();
   h_goodEleSCEta_2SC2Jets->Write();
   h_goodEleSCPt_Barrel->Write();
   h_goodEleSCPt_Barrel_2SC2Jets->Write();
   h_goodEleSCPt_Endcap->Write();
   h_goodEleSCPt_Endcap_2SC2Jets->Write();
   h_goodSCPt->Write();
   h_goodSCEta->Write();
   h_goodSCPt_2SC2Jets->Write();
   h_goodSCEta_2SC2Jets->Write();
   h_goodSCPt_Barrel->Write();
   h_goodSCPt_Barrel_2SC2Jets->Write();
   h_goodSCPt_Endcap->Write();
   h_goodSCPt_Endcap_2SC2Jets->Write();

   h_scEcalIso_Et->Write();
   h_scHoE->Write();
   h_scTrkIso->Write();
   h_scIetaIeta->Write();

   h_eta_failHLT->Write();
   h_phi_failHLT->Write();

   h_probPt1stSc_esjj->Write();
   h_actualPt1stSc_esjj->Write();
   h_probSt_esjj->Write();
   h_actualSt_esjj->Write();
   h_probPt1stSc_eejj->Write();
   h_actualPt1stSc_eejj->Write();
   h_probSt_eejj->Write();
   h_actualSt_eejj->Write();

   FailRate = 100 * NFailHLT/HasSC;
   //cout << "NFail: " << NFailHLT << "\t" << "FailRate: " << FailRate << " %" << endl;
  
  //STDOUT("analysisClass::Loop() ends");   
}
