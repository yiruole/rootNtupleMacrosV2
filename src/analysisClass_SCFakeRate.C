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
   
   TH1F *h_dPhi_JetSC = new TH1F("dPhi_JetSC","dPhi_JetSC",650,0,6.5); h_dPhi_JetSC->Sumw2();
   TH1F *h_dR_JetSC = new TH1F("dR_JetSC","dR_JetSC",600,0,3.0); h_dR_JetSC->Sumw2();
   TH1F *h_NisoSC = new TH1F ("NisoSC","NisoSC",6,-0.5,5.5);  h_NisoSC->Sumw2();

   TH1F *h_goodEleSCPt = new TH1F ("goodEleSCPt","goodEleSCPt",100,0,100); h_goodEleSCPt->Sumw2();
   TH1F *h_goodEleSCEta = new TH1F ("goodEleSCEta","goodEleSCEta",100,-3.,3.); h_goodEleSCEta->Sumw2();
   TH1F *h_goodEleSCPt_Barrel = new TH1F ("goodEleSCPt_Barrel","goodEleSCPt_Barrel",100,0,100); h_goodEleSCPt_Barrel->Sumw2();
   TH1F *h_goodEleSCPt_Endcap = new TH1F ("goodEleSCPt_Endcap","goodEleSCPt_Endcap",100,0,100); h_goodEleSCPt_Endcap->Sumw2();

   TH1F *h_goodSCPt = new TH1F ("goodSCPt","goodSCPt",100,0,100); h_goodSCPt->Sumw2();
   TH1F *h_goodSCEta = new TH1F ("goodSCEta","goodSCEta",100,-3.,3.); h_goodSCEta->Sumw2();
   TH1F *h_goodSCPt_Barrel = new TH1F ("goodSCPt_Barrel","goodSCPt_Barrel",100,0,100); h_goodSCPt_Barrel->Sumw2();
   TH1F *h_goodSCPt_Endcap = new TH1F ("goodSCPt_Endcap","goodSCPt_Endcap",100,0,100); h_goodSCPt_Endcap->Sumw2();

   TH1F *h_loose_goodEleSCPt = new TH1F ("loose_goodEleSCPt","loose_goodEleSCPt",100,0,100); h_loose_goodEleSCPt->Sumw2();
   TH1F *h_loose_goodEleSCEta = new TH1F ("loose_goodEleSCEta","loose_goodEleSCEta",100,-3.,3.); h_loose_goodEleSCEta->Sumw2();
   TH1F *h_loose_goodEleSCPt_Barrel = new TH1F ("loose_goodEleSCPt_Barrel","loose_goodEleSCPt_Barrel",100,0,100); h_loose_goodEleSCPt_Barrel->Sumw2();
   TH1F *h_loose_goodEleSCPt_Endcap = new TH1F ("loose_goodEleSCPt_Endcap","loose_goodEleSCPt_Endcap",100,0,100); h_loose_goodEleSCPt_Endcap->Sumw2();

   TH1F *h_loose_goodSCPt = new TH1F ("loose_goodSCPt","loose_goodSCPt",100,0,100); h_loose_goodSCPt->Sumw2();
   TH1F *h_loose_goodSCEta = new TH1F ("loose_goodSCEta","loose_goodSCEta",100,-3.,3.); h_loose_goodSCEta->Sumw2();
   TH1F *h_loose_goodSCPt_Barrel = new TH1F ("loose_goodSCPt_Barrel","loose_goodSCPt_Barrel",100,0,100); h_loose_goodSCPt_Barrel->Sumw2();
   TH1F *h_loose_goodSCPt_Endcap = new TH1F ("loose_goodSCPt_Endcap","loose_goodSCPt_Endcap",100,0,100); h_loose_goodSCPt_Endcap->Sumw2();

   TH1F *h_eta_failHLT = new TH1F("eta_failHLT","eta_failHLT",500,-3.0,3.0);
   TH1F *h_phi_failHLT = new TH1F("phi_failHLT","phi_failHLT",100,-3.5,3.5);

   TH1F *h_probPt1stSc = new TH1F ("probPt1stSc","probPt1stSc",200,0,1000);  h_probPt1stSc->Sumw2();  //N events based on fake rate
   TH1F *h_actualPt1stSc = new TH1F ("actualPt1stSc","actualPt1stSc",200,0,1000);  h_actualPt1stSc->Sumw2();  //N events with at least 1 HEEP ele
   TH1F *h_probSt = new TH1F ("probSt","probSt",200,0,1000);  h_probSt->Sumw2();  //N events based on fake rate
   TH1F *h_actualSt = new TH1F ("actualSt","actualSt",200,0,1000);  h_actualSt->Sumw2();  //N events with at least 1 HEEP ele

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
  //for (Long64_t jentry=0; jentry<10000;jentry++) { // Begin of loop over events
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry < 10 || jentry%10000 == 0) STDOUT("analysisClass::Loop(): jentry = " << jentry);   
    // if (Cut(ientry) < 0) continue;
    
    ////////////////////// User's code to be done for every event - BEGIN ///////////////////////

    //## HLT

    bool PassTrig=HLTResults->at(1); // results of HLTPhoton15 
    //bool PassTrig=HLTBits->at(71); // results of HLTPhoton15 
    int TrigDiff = HLTBits->at(71) - HLTResults->at(1);
    h_TrigDiff->Fill(TrigDiff);
    h2_DebugTrig->Fill(HLTResults->at(1),HLTBits->at(71));

    // Electrons
    vector<int> v_idx_ele_all;
    vector<int> v_idx_ele_PtCut;
    vector<int> v_idx_ele_HEEP;
    vector<int> v_idx_ele_HEEP_loose;
    int eleIDType = (int) getPreCutValue1("eleIDType");
     
    //Loop over electrons
    for(int iele=0; iele<ElectronPt->size(); iele++)
      {
	//no cut on reco electrons
	v_idx_ele_all.push_back(iele); 

	//pT pre-cut on ele
	if( ElectronPt->at(iele) < getPreCutValue1("ele_PtCut") ) continue; 
	v_idx_ele_PtCut.push_back(iele);

	//ID + ISO + NO overlap with good muons 
	int eleID = ElectronPassID->at(iele);
	if ( (eleID & 1<<eleIDType) > 0  && ElectronOverlaps->at(iele)==0 )
	  {
	    v_idx_ele_HEEP.push_back(iele);
	  }

// 	bool HEEP = false;
// 	if (fabs(ElectronEta->at(iele))<1.442){
// 	if (fabs(ElectronDeltaEtaTrkSC->at(iele))<0.005){
// 	  if (fabs(ElectronDeltaPhiTrkSC->at(iele))<0.09){
// 	    if (ElectronHoE->at(iele)<0.05){
// 		if ((ElectronE2x5OverE5x5->at(iele)>0.94)||(ElectronE1x5OverE5x5->at(iele)>0.83)){
// 		  if ((ElectronEcalIsoHeep->at(iele) + ElectronHcalIsoD1Heep->at(iele))<(2+0.03*ElectronPt->at(iele))){
// 		    if (ElectronTrkIsoHeep->at(iele)<7.5){
// 		      HEEP = true;
// 		    }
// 		  }
// 		}
// 	      }
// 	    }
// 	  }
// 	}

// 	double eleEt = ElectronPt->at(iele);
// 	if ((fabs(ElectronEta->at(iele))>1.56) && (fabs(ElectronEta->at(iele))<2.5)){
// 	if (fabs(ElectronDeltaEtaTrkSC->at(iele))<0.007){
// 	  if (fabs(ElectronDeltaPhiTrkSC->at(iele))<0.09){
// 	    if (ElectronHoE->at(iele)<0.05){
// 	      if (ElectronSigmaIEtaIEta->at(iele)<0.03){
// 		if (((ElectronEcalIsoHeep->at(iele) + ElectronHcalIsoD1Heep->at(iele))<2.5) || 
// 		    ((eleEt>50)&&(ElectronEcalIsoHeep->at(iele) + ElectronHcalIsoD1Heep->at(iele))< ((0.03*(eleEt-50))+2.5))){
// 		    if (ElectronTrkIsoHeep->at(iele)<15){
// 		      if (ElectronHcalIsoD2Heep->at(iele)<0.5){
// 		      HEEP = true;
// 		      }
// 		    }
// 		}
// 	      }
// 	    }
// 	  }
// 	}
// 	}

// 	if ( HEEP  && ElectronOverlaps->at(iele)==0 )
// 	  {
// 	    v_idx_ele_HEEP.push_back(iele);
// 	  }

	bool HEEP_loose = false;
	if (fabs(ElectronEta->at(iele))<1.442){
	if (fabs(ElectronDeltaEtaTrkSC->at(iele))<0.005){
	  if (fabs(ElectronDeltaPhiTrkSC->at(iele))<0.09){
	    if (ElectronHoE->at(iele)<0.05){
		      HEEP_loose = true;
	      }
	    }
	  }
	}

	if ((fabs(ElectronEta->at(iele))>1.56) && (fabs(ElectronEta->at(iele))<2.5)){
	if (fabs(ElectronDeltaEtaTrkSC->at(iele))<0.007){
	  if (fabs(ElectronDeltaPhiTrkSC->at(iele))<0.09){
	    if (ElectronHoE->at(iele)<0.05){
	      if (ElectronSigmaIEtaIEta->at(iele)<0.03){
		      HEEP_loose = true;
	      }
	    }
	  }
	}
	}

	if ( HEEP_loose  && ElectronOverlaps->at(iele)==0 )
	  {
	    v_idx_ele_HEEP_loose.push_back(iele);
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
      if (SuperClusterHoE->at(isc)<0.05) {
	v_idx_sc_iso.push_back(isc);
      }
    }

    // Jets
    vector<int> v_idx_jet_all;
    vector<int> v_idx_jet_PtCut;
    vector<int> v_idx_jet_PtCut_noOvrlap;
 
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
	if( ( CaloJetOverlaps->at(ijet) & 1 << eleIDType) == 0 && minDRsc>0.3 )/* NO overlap with sc */  
	//if( ( CaloJetOverlaps->at(ijet) & 1 << eleIDType) == 0 )  
	  v_idx_jet_PtCut_noOvrlap.push_back(ijet);

      } // End loop over jets

     //find dR between isoSC and jets

    float smallest_ScJet_dR = 99;
    for (int ijet=0; ijet<v_idx_jet_PtCut_noOvrlap.size(); ijet++){
      TVector3 jet_vec;
      jet_vec.SetPtEtaPhi(CaloJetPt->at(v_idx_jet_PtCut_noOvrlap[ijet]),
			       CaloJetEta->at(v_idx_jet_PtCut_noOvrlap[ijet]),
			       CaloJetPhi->at(v_idx_jet_PtCut_noOvrlap[ijet]));
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
       if (minDR<smallest_ScJet_dR) smallest_ScJet_dR= minDR;
    }


    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();
    

    // Set the value of the variableNames listed in the cutFile to their current value
    
    // HLT
    fillVariableWithValue( "HLT", PassTrig ) ;

    // dR SC-Jet
    fillVariableWithValue( "dR_SCJet", smallest_ScJet_dR );

    //## SC
    fillVariableWithValue( "nIsoSC", v_idx_sc_iso.size() );
    //fillVariableWithValue( "nIsoSC",v_idx_sc_iso_barrel.size() );

    // nEle
    fillVariableWithValue( "nEle_all", v_idx_ele_all.size() ) ;
    fillVariableWithValue( "nEle_PtCut", v_idx_ele_PtCut.size() ) ;
    fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp", v_idx_ele_HEEP.size() ) ;
     
    // nJet
    fillVariableWithValue( "nJet_all", v_idx_jet_all.size() ) ;
    fillVariableWithValue( "nJet_PtCut", v_idx_jet_PtCut.size() ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlpSC", v_idx_jet_PtCut_noOvrlap.size() ) ;

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
      }

    // 1st jet
    if( v_idx_jet_PtCut_noOvrlap.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stJet_noOvrlpSC", CaloJetPt->at(v_idx_jet_PtCut_noOvrlap[0]) );
	fillVariableWithValue( "Eta1stJet_noOvrlpSC", CaloJetEta->at(v_idx_jet_PtCut_noOvrlap[0]) );
	fillVariableWithValue( "mEta1stJet_noOvrlpSC", fabs(CaloJetEta->at(v_idx_jet_PtCut_noOvrlap[0])) );
      }


    //cout << "2nd Jet" << endl;
    //## 2nd jet
    if( v_idx_jet_PtCut_noOvrlap.size() >= 2 ) 
      {
	fillVariableWithValue( "Pt2ndJet_noOvrlpSC", CaloJetPt->at(v_idx_jet_PtCut_noOvrlap[1]) );
	fillVariableWithValue( "Eta2ndJet_noOvrlpSC", CaloJetEta->at(v_idx_jet_PtCut_noOvrlap[1]) );
	fillVariableWithValue( "mEta2ndJet_noOvrlpSC", fabs(CaloJetEta->at(v_idx_jet_PtCut_noOvrlap[1])) );
	fillVariableWithValue( "maxMEtaJets_noOvrlpSC", max( getVariableValue("mEta1stJet_noOvrlpSC"), getVariableValue("mEta2ndJet_noOvrlpSC") ) );
      }

    //## define "2ele" and "2jets" booleans
    bool TwoEle=false;
    bool TwoSC=false;
    bool TwoJets=false;
    if( v_idx_ele_HEEP.size() >= 2 ) TwoEle = true;
    if( v_idx_jet_PtCut_noOvrlap.size() >= 2 ) TwoJets = true;
    if( v_idx_sc_iso.size()>= 2) TwoSC = true;

    // ST
    double calc_sT=-999.; 
    if ( (TwoSC) && (TwoJets) ) 
      {
	calc_sT = 
	  SuperClusterPt->at(v_idx_sc_iso[0]) +
	  SuperClusterPt->at(v_idx_sc_iso[1]) +
	  CaloJetPt->at(v_idx_jet_PtCut_noOvrlap[0]) +
	  CaloJetPt->at(v_idx_jet_PtCut_noOvrlap[1]);
	fillVariableWithValue("sT", calc_sT);
      }

    // Mee
    double MEE = 0;
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
	MEE = ee.M();
	fillVariableWithValue("Mee", ee.M());
      }

    // Evaluate cuts (but do not apply them)
    evaluateCuts();

     if (passedCut("nIsoSC")){
       HasSC++;
     }

     // Fill histograms and do analysis based on cut evaluation
     if (v_idx_sc_iso.size()>0 && !passedCut("HLT")){
       NFailHLT++;
       h_eta_failHLT->Fill(SuperClusterEta->at(v_idx_sc_iso[0]));
       h_phi_failHLT->Fill(SuperClusterPhi->at(v_idx_sc_iso[0]));
     }

     if (!passedCut("dR_SCJet")) continue;  //leading 2 sc not close to jets

     if( passedCut("0")&&passedCut("1")&&passedCut("2")&&passedCut("3") ) 
       {
	 if (v_idx_ele_HEEP_loose.size()==1){
	   h_actualPt1stSc->Fill(ElectronSCPt->at(v_idx_ele_HEEP_loose[0]));
	   h_actualSt->Fill(calc_sT);
	 }

	 double probSC1 = 0, probSC2 = 0;
	 double BarrelCross = getPreCutValue1("FakeRate_BarrelCross");
	 double BarrelSlope = getPreCutValue1("FakeRate_BarrelSlope");
	 double EndcapCross = getPreCutValue1("FakeRate_EndcapCross");
	 double EndcapSlope = getPreCutValue1("FakeRate_EndcapSlope");

	 if (fabs(SuperClusterEta->at(v_idx_sc_iso[0]))<1.442) probSC1 = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_iso[0]);
	 if (fabs(SuperClusterEta->at(v_idx_sc_iso[0]))>1.56) probSC1 = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_iso[0]) ;
	 if (fabs(SuperClusterEta->at(v_idx_sc_iso[1]))<1.442) probSC2 = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_iso[1]);
	 if (fabs(SuperClusterEta->at(v_idx_sc_iso[1]))>1.56) probSC2 = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_iso[1]);
      
	 h_probPt1stSc->Fill(SuperClusterPt->at(v_idx_sc_iso[0]),probSC1+probSC2);
	 h_probSt->Fill(calc_sT,probSC1+probSC2);
       }

     //// Fill fake rate pltos
     //

     if ( (MEE>60)&&(MEE<120)) continue; // require at most one loose ele to reject Zs
     if (PFMET->at(0) > 10) continue; // get rid of Ws
	 for(int isc=0;isc<v_idx_sc_iso.size();isc++)
	   {
	     //Require dR>3.0 between SC and JET
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
	     if (dPhi_SC_Jet<3.0) continue;
	     h_dPhi_JetSC->Fill(dPhi_SC_Jet);	

	     h_goodSCPt->Fill(SuperClusterPt->at(v_idx_sc_iso[isc]));
	     h_goodSCEta->Fill(SuperClusterEta->at(v_idx_sc_iso[isc]));
	     h_loose_goodSCPt->Fill(SuperClusterPt->at(v_idx_sc_iso[isc]));
	     h_loose_goodSCEta->Fill(SuperClusterEta->at(v_idx_sc_iso[isc]));

	     bool Barrel = false;
	     bool Endcap = false;
	     if (fabs(SuperClusterEta->at(v_idx_sc_iso[isc]))<1.45) Barrel = true;
	     if (fabs(SuperClusterEta->at(v_idx_sc_iso[isc]))>1.56 && fabs(SuperClusterEta->at(v_idx_sc_iso[isc]))<2.5) Endcap = true;

	     if (Barrel) h_goodSCPt_Barrel->Fill(SuperClusterPt->at(v_idx_sc_iso[isc]));
	     if (Endcap) h_goodSCPt_Endcap->Fill(SuperClusterPt->at(v_idx_sc_iso[isc]));
	     if (Barrel) h_loose_goodSCPt_Barrel->Fill(SuperClusterPt->at(v_idx_sc_iso[isc]));
	     if (Endcap) h_loose_goodSCPt_Endcap->Fill(SuperClusterPt->at(v_idx_sc_iso[isc]));

	     ///  see if there is a HEEP ele to match this
	     double deltaR_ele_sc = 99;
	     int idx_HEEP = -1;
	     for(int iele=0;iele<v_idx_ele_HEEP.size();iele++){
	       TVector3 ele_vec;
	       ele_vec.SetPtEtaPhi(ElectronPt->at(v_idx_ele_HEEP[iele]),
			   ElectronEta->at(v_idx_ele_HEEP[iele]),
			   ElectronPhi->at(v_idx_ele_HEEP[iele]));
	       double tmp_deltaR = ele_vec.DeltaR(sc_vec);
	       if (tmp_deltaR<deltaR_ele_sc){
		 deltaR_ele_sc = tmp_deltaR;
		 idx_HEEP = iele;
	       }
	     }
	     if (deltaR_ele_sc<0.3){
		h_goodEleSCPt->Fill(ElectronSCPt->at(v_idx_ele_HEEP[idx_HEEP]));
		h_goodEleSCEta->Fill(ElectronSCEta->at(v_idx_ele_HEEP[idx_HEEP]));
		
		bool Barrel = false;
		bool Endcap = false;
		if (fabs(ElectronSCEta->at(v_idx_ele_HEEP[idx_HEEP]))<1.45) Barrel = true;
		if (fabs(ElectronSCEta->at(v_idx_ele_HEEP[idx_HEEP]))>1.56 && fabs(ElectronSCEta->at(v_idx_ele_HEEP[idx_HEEP]))<2.5) Endcap = true;
		
		if (Barrel) h_goodEleSCPt_Barrel->Fill(ElectronSCPt->at(v_idx_ele_HEEP[idx_HEEP]));
		if (Endcap) h_goodEleSCPt_Endcap->Fill(ElectronSCPt->at(v_idx_ele_HEEP[idx_HEEP]));
		}


	     ///  see if there is a HEEP_loose ele to match this
	     double deltaR_loose_ele_sc = 99;
	     int idx_HEEP_loose = -1;
	     for(int iele=0;iele<v_idx_ele_HEEP_loose.size();iele++){
	       TVector3 loose_ele_vec;
	       loose_ele_vec.SetPtEtaPhi(ElectronPt->at(v_idx_ele_HEEP_loose[iele]),
			   ElectronEta->at(v_idx_ele_HEEP_loose[iele]),
			   ElectronPhi->at(v_idx_ele_HEEP_loose[iele]));
	       double tmp_deltaR = loose_ele_vec.DeltaR(sc_vec);
	       if (tmp_deltaR<deltaR_loose_ele_sc){
		 deltaR_loose_ele_sc = tmp_deltaR;
		 idx_HEEP_loose = iele;
	       }
	     }
	     if (deltaR_loose_ele_sc<0.3){
		h_loose_goodEleSCPt->Fill(ElectronSCPt->at(v_idx_ele_HEEP_loose[idx_HEEP_loose]));
		h_loose_goodEleSCEta->Fill(ElectronSCEta->at(v_idx_ele_HEEP_loose[idx_HEEP_loose]));
		
		bool Barrel = false;
		bool Endcap = false;
		if (fabs(ElectronSCEta->at(v_idx_ele_HEEP_loose[idx_HEEP_loose]))<1.45) Barrel = true;
		if (fabs(ElectronSCEta->at(v_idx_ele_HEEP_loose[idx_HEEP_loose]))>1.56 && fabs(ElectronSCEta->at(v_idx_ele_HEEP_loose[idx_HEEP_loose]))<2.5) Endcap = true;
		
		if (Barrel) h_loose_goodEleSCPt_Barrel->Fill(ElectronSCPt->at(v_idx_ele_HEEP_loose[idx_HEEP_loose]));
		if (Endcap) h_loose_goodEleSCPt_Endcap->Fill(ElectronSCPt->at(v_idx_ele_HEEP_loose[idx_HEEP_loose]));
		}
	   } // for sc


    ////////////////////// User's code to be done for every event - END ///////////////////////
    
  } // End of loop over events
  

  ////////////////////// User's code to write histos - BEGIN ///////////////////////

   h_TrigDiff->Write();
   h2_DebugTrig->Write();

   h_dPhi_JetSC->Write();
   h_dR_JetSC->Write();
   h_NisoSC->Write();
   h_goodEleSCPt->Write();
   h_goodEleSCEta->Write();
   h_goodEleSCPt_Barrel->Write();
   h_goodEleSCPt_Endcap->Write();
   h_loose_goodEleSCPt->Write();
   h_loose_goodEleSCEta->Write();
   h_loose_goodEleSCPt_Barrel->Write();
   h_loose_goodEleSCPt_Endcap->Write();
   h_goodSCPt->Write();
   h_goodSCEta->Write();
   h_goodSCPt_Barrel->Write();
   h_goodSCPt_Endcap->Write();
   h_loose_goodSCPt->Write();
   h_loose_goodSCEta->Write();
   h_loose_goodSCPt_Barrel->Write();
   h_loose_goodSCPt_Endcap->Write();

   h_eta_failHLT->Write();
   h_phi_failHLT->Write();

   h_probPt1stSc->Write();
   h_actualPt1stSc->Write();
   h_probSt->Write();
   h_actualSt->Write();

   FailRate = 100 * NFailHLT/HasSC;
   //cout << "NFail: " << NFailHLT << "\t" << "FailRate: " << FailRate << " %" << endl;
  
  //STDOUT("analysisClass::Loop() ends");   
}
