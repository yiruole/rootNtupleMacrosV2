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

  TH1F *h_Mmj_PAS = new TH1F ("h_Mmj_PAS","h_Mmj_PAS",200,0,2000);  h_Mmj_PAS->Sumw2();
  
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

  double muon_PtCut = getPreCutValue1("muon_PtCut");
  double muFidRegion = getPreCutValue1("muFidRegion"); // currently unused !!!
  double muNHits_minThresh = getPreCutValue1("muNHits_minThresh");
  double muTrkD0Maximum = getPreCutValue1("muTrkD0Maximum");

  ////////////////////// User's code to get preCut values - END /////////////////
    
  Long64_t nentries = fChain->GetEntriesFast();
  STDOUT("analysisClass::Loop(): nentries = " << nentries);   
  
  ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
  ////// If the root version is updated and rootNtupleClass regenerated,     /////
  ////// these lines may need to be updated.                                 /////    
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) { // Begin of loop over events
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

    // Electrons
    vector<int> v_idx_ele_all;
    vector<int> v_idx_ele_PtCut;
    vector<int> v_idx_ele_PtCut_IDISO_noOverlap;
    int heepBitMask;

    //Loop over electrons
    for(int iele=0; iele<ElectronPt->size(); iele++)
      {

	// Reject ECAL spikes
	// TEMPORARILY REMOVED TO RUN ON LQ MUON MC ROOTTUPLES	if ( 1 - ElectronSCS4S1->at(iele) > 0.95 ) continue; 

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

      } // End loop over jets
    

    // Muons
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
      
      if ( ((*MuonTrkHits)[imuon]  >= muNHits_minThresh  )&&( fabs((*MuonTrkD0)[imuon]) < muTrkD0Maximum ) &&((*MuonPassIso)[imuon]==1 ) &&((*MuonPassID)[imuon]==1) ) 
	{
	  v_idx_muon_PtCut_IDISO.push_back(imuon);
	}

    }// end loop over muons


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
//     fillVariableWithValue( "nEle_all", v_idx_ele_all.size() ) ;
//     fillVariableWithValue( "nEle_PtCut", v_idx_ele_PtCut.size() ) ;
//     fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp", v_idx_ele_PtCut_IDISO_noOverlap.size() ) ;

    // nMu
    fillVariableWithValue( "nMu_all", v_idx_muon_all.size() ) ;
    fillVariableWithValue( "nMu_PtCut", v_idx_muon_PtCut.size() ) ;
    fillVariableWithValue( "nMu_PtCut_IDISO", v_idx_muon_PtCut_IDISO.size() ) ;
     
    // nJet
    fillVariableWithValue( "nJet_all", v_idx_jet_all.size() ) ;
    fillVariableWithValue( "nJet_PtCut", v_idx_jet_PtCut.size() ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlp", v_idx_jet_PtCut_noOverlap.size() ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlp_ID", v_idx_jet_PtCut_noOverlap_ID.size() ) ;
    //TwoMuOnly
    fillVariableWithValue( "nJet_TwoMuOnly_All", v_idx_jet_PtCut_noOverlap_ID.size() ) ;
    fillVariableWithValue( "nJet_TwoMuOnly_EtaCut", v_idx_jet_PtCut_noOverlap_ID_EtaCut.size() ) ;
    //PAS June 2010
    fillVariableWithValue( "nJet_PAS_All", v_idx_jet_PtCut_noOverlap_ID.size() ) ;
    fillVariableWithValue( "nJet_PAS_EtaCut", v_idx_jet_PtCut_noOverlap_ID_EtaCut.size() ) ;

    // MET
    //PAS June 2010
    fillVariableWithValue( "pfMET_PAS", PFMET->at(0) ) ;
    fillVariableWithValue( "tcMET_PAS", TCMET->at(0) ) ;
    fillVariableWithValue( "caloMET_PAS", CaloMET->at(0) ) ;

    // 1st mu
    if( v_idx_muon_PtCut_IDISO.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stMu_IDISO", MuonPt->at(v_idx_muon_PtCut_IDISO[0]) );
	fillVariableWithValue( "Eta1stMu_IDISO", MuonEta->at(v_idx_muon_PtCut_IDISO[0]) );
	fillVariableWithValue( "mEta1stMu_IDISO", fabs(MuonEta->at(v_idx_muon_PtCut_IDISO[0])) );
	//PAS June 2010
	fillVariableWithValue( "Pt1stMu_PAS", MuonPt->at(v_idx_muon_PtCut_IDISO[0]) );
	fillVariableWithValue( "Eta1stMu_PAS", MuonEta->at(v_idx_muon_PtCut_IDISO[0]) );
      }

    // 2nd mu
    if( v_idx_muon_PtCut_IDISO.size() >= 2 ) 
      {
	fillVariableWithValue( "Pt2ndMu_IDISO", MuonPt->at(v_idx_muon_PtCut_IDISO[1]) );
	fillVariableWithValue( "Eta2ndMu_IDISO", MuonEta->at(v_idx_muon_PtCut_IDISO[1]) );
	fillVariableWithValue( "mEta2ndMu_IDISO", fabs(MuonEta->at(v_idx_muon_PtCut_IDISO[1])) );
	fillVariableWithValue( "maxMEtaMus_IDISO", max( getVariableValue("mEta1stMu_IDISO"), getVariableValue("mEta2ndMu_IDISO") ) );
	//PAS June 2010
	fillVariableWithValue( "Pt2ndMu_PAS", MuonPt->at(v_idx_muon_PtCut_IDISO[1]) );
	fillVariableWithValue( "Eta2ndMu_PAS", MuonEta->at(v_idx_muon_PtCut_IDISO[1]) );
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

    //## define "2mu" and "2jets" booleans
    bool TwoMu=false;
    bool TwoJets=false;
    if( v_idx_muon_PtCut_IDISO.size() >= 2 ) TwoMu = true;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2 ) TwoJets = true;

    // ST
    double calc_sT=-999.; 
    if ( (TwoMu) && (TwoJets) ) 
      {
	calc_sT = 
	  MuonPt->at(v_idx_muon_PtCut_IDISO[0]) +
	  MuonPt->at(v_idx_muon_PtCut_IDISO[1]) +
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


    // Mmm
    if (TwoMu)
      {
	TLorentzVector mu1, mu2, mm;
	mu1.SetPtEtaPhiM(MuonPt->at(v_idx_muon_PtCut_IDISO[0]),
			  MuonEta->at(v_idx_muon_PtCut_IDISO[0]),
			  MuonPhi->at(v_idx_muon_PtCut_IDISO[0]),0);
	mu2.SetPtEtaPhiM(MuonPt->at(v_idx_muon_PtCut_IDISO[1]),
			  MuonEta->at(v_idx_muon_PtCut_IDISO[1]),
			  MuonPhi->at(v_idx_muon_PtCut_IDISO[1]),0);
	mm = mu1+mu2;
	fillVariableWithValue("Mmm", mm.M());
	//TwoMuOnly
	fillVariableWithValue("Mmm_TwoMuOnly", mm.M());
	fillVariableWithValue("Mmm_presel", mm.M());
	//
	//PAS June 2010
	fillVariableWithValue("Mmm_PAS", mm.M());

	double calc_sTmuon=-999.;
        calc_sTmuon =
	   MuonPt->at(v_idx_muon_PtCut_IDISO[0]) +
	   MuonPt->at(v_idx_muon_PtCut_IDISO[1]) ;
	fillVariableWithValue("sTmuon_PAS", calc_sTmuon);

	if(isData==true) 
	  {
	    STDOUT("Two muons: Run, LS, Event = "<<run<<", "<<ls<<", "<<event);
	    STDOUT("Two muons: M_mm, Pt_mm, Eta_mm, Phi_mm = "<<mm.M() <<", "<< mm.Pt() <<", "<< mm.Eta() <<", "<< mm.Phi());
	    STDOUT("Two muons: 1st mu Pt, eta, phi = "<< mu1.Pt() <<", "<< mu1.Eta() <<", "<< mu1.Phi() );
	    STDOUT("Two muons: 2nd mu Pt, eta, phi = "<< mu2.Pt() <<", "<< mu2.Eta() <<", "<< mu2.Phi() );
	  }
      }

    // Mmj 
    double Mm1j1, Mm1j2, Mm2j1, Mm2j2 = -999;
    double deltaM_m1j1_m2j2 = 9999;
    double deltaM_m1j2_m2j1 = 9999;
    double Mmj_1stPair = 0;
    double Mmj_2ndPair = 0;
    double deltaR_m1j1 ;
    double deltaR_m2j2 ;
    double deltaR_m1j2 ;
    double deltaR_m2j1 ;
    if ( (TwoMu) && (TwoJets) ) // TwoMu and TwoJets
      {
	TLorentzVector jet1, jet2, mu1, mu2;
	mu1.SetPtEtaPhiM(MuonPt->at(v_idx_muon_PtCut_IDISO[0]),
			  MuonEta->at(v_idx_muon_PtCut_IDISO[0]),
			  MuonPhi->at(v_idx_muon_PtCut_IDISO[0]),0);
	mu2.SetPtEtaPhiM(MuonPt->at(v_idx_muon_PtCut_IDISO[1]),
			  MuonEta->at(v_idx_muon_PtCut_IDISO[1]),
			  MuonPhi->at(v_idx_muon_PtCut_IDISO[1]),0);
	jet1.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),0);
	jet2.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]),0);
	TLorentzVector m1j1, m1j2, m2j1, m2j2;
	m1j1 = mu1 + jet1;
	m2j2 = mu2 + jet2;
	m2j1 = mu2 + jet1;
	m1j2 = mu1 + jet2;
	Mm1j1 = m1j1.M();
	Mm2j2 = m2j2.M();
	Mm1j2 = m1j2.M();
	Mm2j1 = m2j1.M();

	deltaM_m1j1_m2j2 = Mm1j1 - Mm2j2;
	deltaM_m1j2_m2j1 = Mm1j2 - Mm2j1;


	double deltaR_m1j1 = mu1.DeltaR(jet1);
	double deltaR_m2j2 = mu2.DeltaR(jet2);
	double deltaR_m1j2 = mu1.DeltaR(jet2);
	double deltaR_m2j1 = mu2.DeltaR(jet1);

// 	// Fill min DR between any of the 2 selected eles and any of the 2 selected jets
// 	double minDR_2ele_2jet = min ( min(deltaR_m1j1,deltaR_m2j2) , min(deltaR_m1j2,deltaR_m2j1) );
// 	fillVariableWithValue("minDR_2ele_2jet", minDR_2ele_2jet);

	if(fabs(deltaM_m1j1_m2j2) > fabs(deltaM_m1j2_m2j1))
	  {
	    Mmj_1stPair = Mm1j2;
	    Mmj_2ndPair = Mm2j1;
	    fillVariableWithValue("minDRmj_selecPairs", min(deltaR_m1j2,deltaR_m2j1) );
	    fillVariableWithValue("minDRmj_unselPairs", min(deltaR_m1j1,deltaR_m2j2) );
	  }
	else
	  {
	    Mmj_1stPair = Mm1j1;
	    Mmj_2ndPair = Mm2j2;
	    fillVariableWithValue("minDRmj_selecPairs", min(deltaR_m1j1,deltaR_m2j2) );
	    fillVariableWithValue("minDRmj_unselPairs", min(deltaR_m1j2,deltaR_m2j1) );
	  } 
	fillVariableWithValue("Mmj_1stPair", Mmj_1stPair);       
	fillVariableWithValue("Mmj_2ndPair", Mmj_2ndPair);
	//PAS June 2010
	h_Mmj_PAS->Fill(Mmj_1stPair);
	h_Mmj_PAS->Fill(Mmj_2ndPair);
	fillVariableWithValue("Mmj_1stPair_PAS", Mmj_1stPair);       
	fillVariableWithValue("Mmj_2ndPair_PAS", Mmj_2ndPair);

	// min and max DeltaR between muons and any jet
	double minDeltaR_mj = 999999;
	double maxDeltaR_mj = -1;
	double thisMinDR, thisMaxDR, DR_thisjet_m1, DR_thisjet_m2;
	TLorentzVector thisjet;
	for(int ijet=0; ijet<v_idx_jet_PtCut_noOverlap_ID.size(); ijet++)
	  {
	    thisjet.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
				 CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
				 CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),0);
	    DR_thisjet_m1 = thisjet.DeltaR(mu1);
	    DR_thisjet_m2 = thisjet.DeltaR(mu2);
	    thisMinDR = min(DR_thisjet_m1, DR_thisjet_m2);
	    thisMaxDR = max(DR_thisjet_m1, DR_thisjet_m2);
	    if(thisMinDR < minDeltaR_mj)
	      minDeltaR_mj = thisMinDR;
	    if(thisMaxDR > maxDeltaR_mj)
	      maxDeltaR_mj = thisMaxDR;
	  } 
	fillVariableWithValue("minDeltaR_mj", minDeltaR_mj);
	fillVariableWithValue("maxDeltaR_mj", maxDeltaR_mj);

	// printouts for small Mmj
	if(isData==true && ( Mmj_1stPair<20 || Mmj_2ndPair<20 ) ) // printouts for low Mmj 
	  {
	    STDOUT("Mmj<20GeV: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
	    STDOUT("Mmj<20GeV: Mmj_1stPair = "<<Mmj_1stPair <<", Mmj_2ndPair = "<< Mmj_2ndPair );
	    STDOUT("Mmj<20GeV: m1j1.M = "<<m1j1.M() <<", m2j2.M = "<<m2j2.M() <<", m1j2.M = "<<m1j2.M()  <<", m2j1.M = "<<m2j1.M()  );
	    STDOUT("Mmj<20GeV: deltaM_m1j1_m2j2 = "<<deltaM_m1j1_m2j2 <<", deltaM_m1j2_m2j1 = "<<deltaM_m1j2_m2j1  );
	    STDOUT("Mmj<20GeV: deltaR_m1j1 = "<<deltaR_m1j1 <<", deltaR_m2j2 = "<<deltaR_m2j2 <<", deltaR_m1j2 = "<<deltaR_m1j2  <<", deltaR_m2j1 = "<<deltaR_m2j1  );
// 	    STDOUT("Mmj<20GeV: 1st mu Pt, eta, phi = "<< mu1.Pt() <<",\t"<< mu1.Eta() <<",\t"<< mu1.Phi() );
// 	    STDOUT("Mmj<20GeV: 2nd mu Pt, eta, phi = "<< mu2.Pt() <<",\t"<< mu2.Eta() <<",\t"<< mu2.Phi() );
// 	    STDOUT("Mmj<20GeV: 1st jet Pt, eta, phi = "<< jet1.Pt() <<",\t"<< jet1.Eta() <<",\t"<< jet1.Phi() );
// 	    STDOUT("Mmj<20GeV: 2nd jet Pt, eta, phi = "<< jet2.Pt() <<",\t"<< jet2.Eta() <<",\t"<< jet2.Phi() );
	    TLorentzVector thismuon;
	    for(int imuon=0; imuon<v_idx_muon_PtCut_IDISO.size(); imuon++)
	      {
		thismuon.SetPtEtaPhiM(MuonPt->at(v_idx_muon_PtCut_IDISO[imuon]),
				     MuonEta->at(v_idx_muon_PtCut_IDISO[imuon]),
				     MuonPhi->at(v_idx_muon_PtCut_IDISO[imuon]),0);
		STDOUT("Mmj<20GeV: e"<<imuon+1<<" Pt, eta, phi = " 
		       << ", "<<thismuon.Pt()<<", "<< thismuon.Eta() <<", "<< thismuon.Phi()<<"; DR_j1, DR_j2 = "<< thismuon.DeltaR(jet1)<<", "<<thismuon.DeltaR(jet2));
	      }
	    TLorentzVector thisjet;
	    TLorentzVector thisjet_m1, thisjet_m2;
	    for(int ijet=0; ijet<v_idx_jet_PtCut_noOverlap_ID.size(); ijet++)
	      {
		thisjet.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
				     CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
				     CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),0);
		thisjet_m1 = thisjet + mu1;
		thisjet_m2 = thisjet + mu2;
		STDOUT("Mmj<20GeV: j"<<ijet+1<<" Pt, eta, phi = " << ", "<<thisjet.Pt()<<", "<< thisjet.Eta() <<", "<< thisjet.Phi()<<"; DR_m1, DR_m2 = "<< thisjet.DeltaR(mu1)<<", "<<thisjet.DeltaR(mu2) << "; M_m1, M_m2 = " <<thisjet_m1.M() <<", "<<thisjet_m2.M() );
	      }
	  } // printouts for low Mmj 

      } // TwoMu and TwoJets



    // Evaluate cuts (but do not apply them)
    evaluateCuts();

    // Fill histograms and do analysis based on cut evaluation
    //h_nMuFinal->Fill( MuonPt->size() );
     

    //Large MET events passing pre-selection
    float METcut = 60;
    if( variableIsFilled("pfMET_PAS") && passedAllPreviousCuts("pfMET_PAS") && isData==true) 
      {

	if( getVariableValue("pfMET_PAS") > METcut )
	  {

	    STDOUT("pfMET>60GeV: ----------- START ------------");

	    STDOUT("pfMET>60GeV: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);	    
	    if( variableIsFilled("Pt1stMu_PAS") && variableIsFilled("Eta1stMu_PAS") )	      
	      STDOUT("pfMET>60GeV: Pt1stMu_PAS,Eta1stMu_PAS = "<<getVariableValue("Pt1stMu_PAS")<<",\t"<<getVariableValue("Eta1stMu_PAS"));
	    if( variableIsFilled("Pt2ndMu_PAS") && variableIsFilled("Eta2ndMu_PAS") )	      
	      STDOUT("pfMET>60GeV: Pt2ndMu_PAS,Eta2ndMu_PAS = "<<getVariableValue("Pt2ndMu_PAS")<<",\t"<<getVariableValue("Eta2ndMu_PAS"));
	    if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") )	      
	      STDOUT("pfMET>60GeV: Pt1stJet_PAS,Eta1stJet_PAS = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS"));
	    if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") )	      
	      STDOUT("pfMET>60GeV: Pt2ndJet_PAS,Eta2ndJet_PAS = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS"));
	    if( variableIsFilled("Mmm_PAS") && variableIsFilled("Mjj_PAS") )	      
	      STDOUT("pfMET>60GeV: Mmm_PAS,Mjj_PAS = "<<getVariableValue("Mmm_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
	    if( variableIsFilled("Mmj_1stPair_PAS") && variableIsFilled("Mmj_2ndPair_PAS") )	      
	      STDOUT("pfMET>60GeV: Mmj_1stPair_PAS,Mmj_2ndPair_PAS = "<<getVariableValue("Mmj_1stPair_PAS")
		     <<",\t"<<getVariableValue("Mmj_2ndPair_PAS"));
	    if( variableIsFilled("pfMET_PAS") && variableIsFilled("caloMET_PAS") )	      
	      STDOUT("pfMET>60GeV: pfMET_PAS,caloMET_PAS = "<<getVariableValue("pfMET_PAS")<<",\t"<<getVariableValue("caloMET_PAS"));

	    STDOUT("pfMET>60GeV: ------------ END -------------");

	  }
      }
    
    if( passedAllPreviousCuts("Pt1stJet_PAS") && (TwoMu) && (TwoJets) )
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

  h_Mmj_PAS->Write();

  h_CaloJetResEMF_looseID->Write();
  h_CaloJetFHPD_looseID->Write();
  h_CaloJetN90Hits_looseID->Write();


  ////////////////////// User's code to write histos - END ///////////////////////
  
  
  //STDOUT("analysisClass::Loop() ends");   
}
