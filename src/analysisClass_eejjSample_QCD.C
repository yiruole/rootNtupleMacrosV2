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

  double muon_PtCut = getPreCutValue1("muon_PtCut");
  double muFidRegion = getPreCutValue1("muFidRegion"); // currently unused !!!
  double muNHits_minThresh = getPreCutValue1("muNHits_minThresh");
  double muTrkD0Maximum = getPreCutValue1("muTrkD0Maximum");

  double BarrelCross = getPreCutValue1("fakeRate_Barrel");
  double BarrelSlope = getPreCutValue2("fakeRate_Barrel");
  double EndcapCross = getPreCutValue1("fakeRate_Endcap");
  double EndcapSlope = getPreCutValue2("fakeRate_Endcap");

  ////////////////////// User's code to get preCut values - END /////////////////
    
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
    if(jentry < 10 || jentry%1000 == 0) STDOUT("analysisClass::Loop(): jentry = " << jentry);   
    // if (Cut(ientry) < 0) continue;
    
    ////////////////////// User's code to be done for every event - BEGIN ///////////////////////

    //if (PtHat>=30) continue;

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
            jet.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut[ijet]),
			     CaloJetEta->at(v_idx_jet_PtCut[ijet]),
			     CaloJetPhi->at(v_idx_jet_PtCut[ijet]),0);
	    double DR = jet.DeltaR(sc);
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
      
      if ( ((*MuonTrkHits)[imuon]  >= muNHits_minThresh  )
	   &&( fabs((*MuonTrkD0)[imuon]) < muTrkD0Maximum ) 
	   &&((*MuonPassIso)[imuon]==1 ) 
	   &&((*MuonPassID)[imuon]==1) ) 
	{
	  v_idx_muon_PtCut_IDISO.push_back(imuon);
	}

    }// end loop over muons


    /////  Define the fake rate for QCD and calculate prob. for each sc to be an ele

    double p1 = 0, p2 = 0, p3 = 0;
    if (v_idx_sc_Iso.size()>=1){
      if (fabs(SuperClusterEta->at(v_idx_sc_Iso[0]))<eleEta_bar) p1 = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_Iso[0]);
      if (fabs(SuperClusterEta->at(v_idx_sc_Iso[0]))>eleEta_end_min) p1 = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_Iso[0]) ;
    }
    if (v_idx_sc_Iso.size()>=2){
      if (fabs(SuperClusterEta->at(v_idx_sc_Iso[1]))<eleEta_bar) p2 = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_Iso[1]);
      if (fabs(SuperClusterEta->at(v_idx_sc_Iso[1]))>eleEta_end_min) p2 = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_Iso[1]);
    }
    if (v_idx_sc_Iso.size()>=3){
      if (fabs(SuperClusterEta->at(v_idx_sc_Iso[2]))<eleEta_bar) p3 = BarrelCross + BarrelSlope*SuperClusterPt->at(v_idx_sc_Iso[2]);
      if (fabs(SuperClusterEta->at(v_idx_sc_Iso[2]))>eleEta_end_min) p3 = EndcapCross + EndcapSlope*SuperClusterPt->at(v_idx_sc_Iso[2]);
    }

    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();
    
    //-----

    // Set the value of the variableNames listed in the cutFile to their current value

    // original code from Ellie available at 
    // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/Leptoquarks/rootNtupleMacrosV2/src/analysisClass_eejjSample_QCD.C?revision=1.8&view=markup
    //here we report only the code to fill the "nEle_PtCut_IDISO_noOvrlp" variable
    
    //## fill bin Nele=0
    
    //     if (v_idx_sc_Iso.size()==0) fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp", v_idx_sc_Iso.size()) ;
    //     else if (v_idx_sc_Iso.size()==1) {
    //       fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp",0, 1-p1 ) ;
    //     }
    //     else if (v_idx_sc_Iso.size()==2) {
    //       fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp",0, (1-p1)*(1-p2) ) ;
    //     }
    //     else if (v_idx_sc_Iso.size()>2) {
    //       fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp",0, (1-p1)*(1-p2)*(1-p3) ) ;
    //     }

    //     evaluateCuts();
    //     resetCuts();

    //## fill bin Nele=1

    //     if (v_idx_sc_Iso.size()==1) fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp",1, p1 ) ;
    //     else if (v_idx_sc_Iso.size()==2) {
    //       fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp",1, p1 + p2 - 2*p1*p2 ) ;
    //     }
    //     else if (v_idx_sc_Iso.size()>2) {
    //       fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp",1, p1 + p2 + p3 - 2*p1*p2 - 2*p2*p3 - 2*p1*p3   ) ;
    //     }

    //     evaluateCuts();
    //     resetCuts();

    //## fill bin Nele=2

    //     if (v_idx_sc_Iso.size()==2) fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp",2, p1*p2) ;
    //     else if (v_idx_sc_Iso.size()>2) {
    //       fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp",2, p1*p3 + p1*p2 + p3*p2 ) ;
    //     }

    //     evaluateCuts();
    //     resetCuts();

    //## fill bin Nele=3

    //     if (v_idx_sc_Iso.size()>2) fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp",3, p1*p2*p3) ;

    //-----

    ///  Just fill all histograms assuming exactly 2 superclusters.
    //     if (v_idx_sc_all.size()>=2) fillVariableWithValue( "nEle_all",v_idx_sc_all.size() , p1*p2 ) ;
    //     if (v_idx_sc_PtCut.size()>=2) fillVariableWithValue( "nEle_PtCut",v_idx_sc_PtCut.size(), p1*p2) ;
    if (v_idx_sc_Iso.size()>=2) fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp",v_idx_sc_Iso.size(), p1*p2) ;


    /// Now fill all variables that are filled just once
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
    //fillVariableWithValue( "PassHBHENoiseFilter", passLooseNoiseFilter ) ;

    // nJet
    fillVariableWithValue( "nJet_all", v_idx_jet_all.size(), p1*p2 ) ;
    fillVariableWithValue( "nJet_PtCut", v_idx_jet_PtCut.size(), p1*p2 ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlp", v_idx_jet_PtCut_noOverlap.size(), p1*p2 ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlp_ID", v_idx_jet_PtCut_noOverlap_ID.size(), p1*p2 ) ;
    //TwoEleOnly
    fillVariableWithValue( "nJet_TwoEleOnly_All", v_idx_jet_PtCut_noOverlap_ID.size() , p1*p2) ;
    fillVariableWithValue( "nJet_TwoEleOnly_EtaCut", v_idx_jet_PtCut_noOverlap_ID_EtaCut.size(), p1*p2 ) ;
    //PAS June 2010
    fillVariableWithValue( "nJet_PAS_All", v_idx_jet_PtCut_noOverlap_ID.size(), p1*p2 ) ;
    fillVariableWithValue( "nJet_PAS_EtaCut", v_idx_jet_PtCut_noOverlap_ID_EtaCut.size(), p1*p2 ) ;

    // MET
    //PAS June 2010
    fillVariableWithValue( "pfMET_PAS", PFMET->at(0), p1*p2 ) ;
    fillVariableWithValue( "tcMET_PAS", TCMET->at(0), p1*p2 ) ;
    fillVariableWithValue( "caloMET_PAS", CaloMET->at(0), p1*p2 ) ;

    //ele Pt
    if( v_idx_sc_Iso.size() >= 2 ) 
      {
	fillVariableWithValue( "Pt1stEle_IDISO_NoOvrlp", SuperClusterPt->at(v_idx_sc_Iso[0]), p1*p2 );
	fillVariableWithValue( "Eta1stEle_IDISO_NoOvrlp", SuperClusterEta->at(v_idx_sc_Iso[0]), p1*p2 );
	fillVariableWithValue( "mEta1stEle_IDISO_NoOvrlp", fabs(SuperClusterEta->at(v_idx_sc_Iso[0])), p1*p2 );
	//PAS June 2010
	fillVariableWithValue( "Pt1stEle_PAS", SuperClusterPt->at(v_idx_sc_Iso[0]), p1*p2 );
	fillVariableWithValue( "Eta1stEle_PAS", SuperClusterEta->at(v_idx_sc_Iso[0]), p1*p2 );

	fillVariableWithValue( "Pt2ndEle_IDISO_NoOvrlp", SuperClusterPt->at(v_idx_sc_Iso[1]), p1*p2 );
	fillVariableWithValue( "Eta2ndEle_IDISO_NoOvrlp", SuperClusterEta->at(v_idx_sc_Iso[1]), p1*p2 );
	fillVariableWithValue( "mEta2ndEle_IDISO_NoOvrlp", fabs(SuperClusterEta->at(v_idx_sc_Iso[1])), p1*p2 );
	fillVariableWithValue( "maxMEtaEles_IDISO_NoOvrl", max( getVariableValue("mEta1stEle_IDISO_NoOvrlp"), getVariableValue("mEta2ndEle_IDISO_NoOvrlp") ) , p1*p2);
	//PAS June 2010
	fillVariableWithValue( "Pt2ndEle_PAS", SuperClusterPt->at(v_idx_sc_Iso[1]), p1*p2 );
	fillVariableWithValue( "Eta2ndEle_PAS", SuperClusterEta->at(v_idx_sc_Iso[1]), p1*p2 );
      }

    // 1st jet
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stJet_noOvrlp_ID", CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]), p1*p2 );
	fillVariableWithValue( "Eta1stJet_noOvrlp_ID", CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]) , p1*p2);
	fillVariableWithValue( "mEta1stJet_noOvrlp_ID", fabs(CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[0])), p1*p2 );
	//PAS June 2010
	fillVariableWithValue( "Pt1stJet_PAS", CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]), p1*p2 );
	fillVariableWithValue( "Eta1stJet_PAS", CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]), p1*p2 );
      }


    //cout << "2nd Jet" << endl;
    //## 2nd jet
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2 ) 
      {
	fillVariableWithValue( "Pt2ndJet_noOvrlp_ID", CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]), p1*p2 );
	fillVariableWithValue( "Eta2ndJet_noOvrlp_ID", CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]), p1*p2 );
	fillVariableWithValue( "mEta2ndJet_noOvrlp_ID", fabs(CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[1])), p1*p2 );
	fillVariableWithValue( "maxMEtaJets_noOvrlp_ID", max( getVariableValue("mEta1stJet_noOvrlp_ID"), getVariableValue("mEta2ndJet_noOvrlp_ID") ), p1*p2 );
	//PAS June 2010
	fillVariableWithValue( "Pt2ndJet_PAS", CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]), p1*p2 );
	fillVariableWithValue( "Eta2ndJet_PAS", CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]), p1*p2 );
      }

    //## define "2ele" and "2jets" booleans
    bool TwoEle=false;
    bool TwoJets=false;
    if( v_idx_sc_Iso.size() >= 2 ) TwoEle = true;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2 ) TwoJets = true;

    // ST
    double calc_sT=-999.; 
    if ( (TwoEle) && (TwoJets) ) 
      {
	calc_sT = 
	  SuperClusterPt->at(v_idx_sc_Iso[0]) +
	  SuperClusterPt->at(v_idx_sc_Iso[1]) +
	  CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) +
	  CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]);
	fillVariableWithValue("sT", calc_sT, p1*p2);
	fillVariableWithValue("sT_MLQ200", calc_sT, p1*p2);
	fillVariableWithValue("sT_MLQ250", calc_sT, p1*p2);
	fillVariableWithValue("sT_MLQ300", calc_sT, p1*p2);
	fillVariableWithValue("sT_MLQ400", calc_sT, p1*p2);       
	//PAS June 2010
	fillVariableWithValue("sT_PAS", calc_sT, p1*p2);
      }

    // ST jets
    if (TwoJets)
      {
	double calc_sTjet=-999.;
	calc_sTjet =
	  CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) + 
	  CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) ;
	fillVariableWithValue("sTjet_PAS", calc_sTjet, p1*p2);
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
	fillVariableWithValue("Mjj_PAS", jj.M(), p1*p2);
      }


    // Mee
    if (TwoEle)
      {
	TLorentzVector ele1, ele2, ee;
	ele1.SetPtEtaPhiM(SuperClusterPt->at(v_idx_sc_Iso[0]),
			  SuperClusterEta->at(v_idx_sc_Iso[0]),
			  SuperClusterPhi->at(v_idx_sc_Iso[0]),0);
	ele2.SetPtEtaPhiM(SuperClusterPt->at(v_idx_sc_Iso[1]),
			  SuperClusterEta->at(v_idx_sc_Iso[1]),
			  SuperClusterPhi->at(v_idx_sc_Iso[1]),0);
	ee = ele1+ele2;
	fillVariableWithValue("Mee", ee.M(), p1*p2);
	//TwoEleOnly
	fillVariableWithValue("Mee_TwoEleOnly", ee.M(), p1*p2);
	fillVariableWithValue("Mee_presel", ee.M(), p1*p2);
	//
	//PAS June 2010
	fillVariableWithValue("Mee_PAS", ee.M(), p1*p2);

	double calc_sTele=-999.;
        calc_sTele =
	  SuperClusterPt->at(v_idx_sc_Iso[0]) +
	  SuperClusterPt->at(v_idx_sc_Iso[1]) ;
	fillVariableWithValue("sTele_PAS", calc_sTele, p1*p2);

	// 	if(isData==true) 
	// 	  {
	// 	    STDOUT("Two electrons: Run, LS, Event = "<<run<<", "<<ls<<", "<<event);
	// 	    STDOUT("Two electrons: M_ee, Pt_ee, Eta_ee, Phi_ee = "<<ee.M() <<", "<< ee.Pt() <<", "<< ee.Eta() <<", "<< ee.Phi());
	// 	    STDOUT("Two electrons: 1st ele Pt, eta, phi = "<< ele1.Pt() <<", "<< ele1.Eta() <<", "<< ele1.Phi() );
	// 	    STDOUT("Two electrons: 2nd ele Pt, eta, phi = "<< ele2.Pt() <<", "<< ele2.Eta() <<", "<< ele2.Phi() );
	// 	  }
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
	ele1.SetPtEtaPhiM(SuperClusterPt->at(v_idx_sc_Iso[0]),
			  SuperClusterEta->at(v_idx_sc_Iso[0]),
			  SuperClusterPhi->at(v_idx_sc_Iso[0]),0);
	ele2.SetPtEtaPhiM(SuperClusterPt->at(v_idx_sc_Iso[1]),
			  SuperClusterEta->at(v_idx_sc_Iso[1]),
			  SuperClusterPhi->at(v_idx_sc_Iso[1]),0);
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
	    fillVariableWithValue("minDRej_selecPairs", min(deltaR_e1j2,deltaR_e2j1), p1*p2 );
	    fillVariableWithValue("minDRej_unselPairs", min(deltaR_e1j1,deltaR_e2j2), p1*p2 );
	  }
	else
	  {
	    Mej_1stPair = Me1j1;
	    Mej_2ndPair = Me2j2;
	    fillVariableWithValue("minDRej_selecPairs", min(deltaR_e1j1,deltaR_e2j2), p1*p2 );
	    fillVariableWithValue("minDRej_unselPairs", min(deltaR_e1j2,deltaR_e2j1), p1*p2 );
	  } 
	fillVariableWithValue("Mej_1stPair", Mej_1stPair, p1*p2);       
	fillVariableWithValue("Mej_2ndPair", Mej_2ndPair, p1*p2);
	//PAS June 2010
	h_Mej_PAS->Fill(Mej_1stPair);
	h_Mej_PAS->Fill(Mej_2ndPair);
	fillVariableWithValue("Mej_1stPair_PAS", Mej_1stPair, p1*p2);       
	fillVariableWithValue("Mej_2ndPair_PAS", Mej_2ndPair, p1*p2);

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
	fillVariableWithValue("minDeltaR_ej", minDeltaR_ej, p1*p2);
	fillVariableWithValue("maxDeltaR_ej", maxDeltaR_ej, p1*p2);

	// 	// printouts for small Mej
	// 	if(isData==true && ( Mej_1stPair<20 || Mej_2ndPair<20 ) ) // printouts for low Mej 
	// 	  {
	// 	    STDOUT("Mej<20GeV: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
	// 	    STDOUT("Mej<20GeV: Mej_1stPair = "<<Mej_1stPair <<", Mej_2ndPair = "<< Mej_2ndPair );
	// 	    STDOUT("Mej<20GeV: e1j1.M = "<<e1j1.M() <<", e2j2.M = "<<e2j2.M() <<", e1j2.M = "<<e1j2.M()  <<", e2j1.M = "<<e2j1.M()  );
	// 	    STDOUT("Mej<20GeV: deltaM_e1j1_e2j2 = "<<deltaM_e1j1_e2j2 <<", deltaM_e1j2_e2j1 = "<<deltaM_e1j2_e2j1  );
	// 	    STDOUT("Mej<20GeV: deltaR_e1j1 = "<<deltaR_e1j1 <<", deltaR_e2j2 = "<<deltaR_e2j2 <<", deltaR_e1j2 = "<<deltaR_e1j2  <<", deltaR_e2j1 = "<<deltaR_e2j1  );
	// // 	    STDOUT("Mej<20GeV: 1st ele Pt, eta, phi = "<< ele1.Pt() <<",\t"<< ele1.Eta() <<",\t"<< ele1.Phi() );
	// // 	    STDOUT("Mej<20GeV: 2nd ele Pt, eta, phi = "<< ele2.Pt() <<",\t"<< ele2.Eta() <<",\t"<< ele2.Phi() );
	// // 	    STDOUT("Mej<20GeV: 1st jet Pt, eta, phi = "<< jet1.Pt() <<",\t"<< jet1.Eta() <<",\t"<< jet1.Phi() );
	// // 	    STDOUT("Mej<20GeV: 2nd jet Pt, eta, phi = "<< jet2.Pt() <<",\t"<< jet2.Eta() <<",\t"<< jet2.Phi() );
	// 	    TLorentzVector thisele;
	// 	    for(int iele=0; iele<v_idx_sc_Iso.size(); iele++)
	// 	      {
	// 		thisele.SetPtEtaPhiM(SuperClusterPt->at(v_idx_sc_Iso[iele]),
	// 				     SuperClusterEta->at(v_idx_sc_Iso[iele]),
	// 				     SuperClusterPhi->at(v_idx_sc_Iso[iele]),0);
	// 		STDOUT("Mej<20GeV: e"<<iele+1<<" Pt, eta, phi = " 
	// 		       << ", "<<thisele.Pt()<<", "<< thisele.Eta() <<", "<< thisele.Phi()<<"; DR_j1, DR_j2 = "<< thisele.DeltaR(jet1)<<", "<<thisele.DeltaR(jet2));
	// 	      }
	// 	    TLorentzVector thisjet;
	// 	    TLorentzVector thisjet_e1, thisjet_e2;
	// 	    for(int ijet=0; ijet<v_idx_jet_PtCut_noOverlap_ID.size(); ijet++)
	// 	      {
	// 		thisjet.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
	// 				     CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
	// 				     CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),0);
	// 		thisjet_e1 = thisjet + ele1;
	// 		thisjet_e2 = thisjet + ele2;
	// 		STDOUT("Mej<20GeV: j"<<ijet+1<<" Pt, eta, phi = " << ", "<<thisjet.Pt()<<", "<< thisjet.Eta() <<", "<< thisjet.Phi()<<"; DR_e1, DR_e2 = "<< thisjet.DeltaR(ele1)<<", "<<thisjet.DeltaR(ele2) << "; M_e1, M_e2 = " <<thisjet_e1.M() <<", "<<thisjet_e2.M() );
	// 	      }
	// 	  } // printouts for low Mej 

      } // TwoEle and TwoJets



    // Evaluate cuts (but do not apply them)
    evaluateCuts();

    // Fill histograms and do analysis based on cut evaluation
    //h_nEleFinal->Fill( ElectronPt->size() );
     

    //     //Large MET events passing pre-selection
    //     float METcut = 60;
    //     if( variableIsFilled("pfMET_PAS") && passedAllPreviousCuts("pfMET_PAS") && isData==true) 
    //       {

    // 	if( getVariableValue("pfMET_PAS") > METcut )
    // 	  {

    // 	    STDOUT("pfMET>60GeV: ----------- START ------------");

    // 	    STDOUT("pfMET>60GeV: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);	    
    // 	    if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") )	      
    // 	      STDOUT("pfMET>60GeV: Pt1stEle_PAS,Eta1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS"));
    // 	    if( variableIsFilled("Pt2ndEle_PAS") && variableIsFilled("Eta2ndEle_PAS") )	      
    // 	      STDOUT("pfMET>60GeV: Pt2ndEle_PAS,Eta2ndEle_PAS = "<<getVariableValue("Pt2ndEle_PAS")<<",\t"<<getVariableValue("Eta2ndEle_PAS"));
    // 	    if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") )	      
    // 	      STDOUT("pfMET>60GeV: Pt1stJet_PAS,Eta1stJet_PAS = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS"));
    // 	    if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") )	      
    // 	      STDOUT("pfMET>60GeV: Pt2ndJet_PAS,Eta2ndJet_PAS = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS"));
    // 	    if( variableIsFilled("Mee_PAS") && variableIsFilled("Mjj_PAS") )	      
    // 	      STDOUT("pfMET>60GeV: Mee_PAS,Mjj_PAS = "<<getVariableValue("Mee_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
    // 	    if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )	      
    // 	      STDOUT("pfMET>60GeV: Mej_1stPair_PAS,Mej_2ndPair_PAS = "<<getVariableValue("Mej_1stPair_PAS")
    // 		     <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
    // 	    if( variableIsFilled("pfMET_PAS") && variableIsFilled("caloMET_PAS") )	      
    // 	      STDOUT("pfMET>60GeV: pfMET_PAS,caloMET_PAS = "<<getVariableValue("pfMET_PAS")<<",\t"<<getVariableValue("caloMET_PAS"));

    // 	    STDOUT("pfMET>60GeV: ------------ END -------------");

    // 	  }
    //       }


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

  ////////////////////// User's code to write histos - END ///////////////////////
  
  
  //STDOUT("analysisClass::Loop() ends");   
}
