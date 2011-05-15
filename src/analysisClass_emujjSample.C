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

// bool JetIdloose(double ak5JetJIDresEMF,double ak5JetJIDfHPD,int ak5JetJIDn90Hits, double ak5JetEta){
//   bool jetidloose=false;
//   bool jetidresEMF=true;

//   double fhpdmax = 0.98;
//   double n90hitsmin =1;
//   double emf_min = 0.01;

//   if(fabs(ak5JetEta)<2.6 && ak5JetJIDresEMF<=emf_min) jetidresEMF=false;

//   if(jetidresEMF && ak5JetJIDfHPD<fhpdmax && ak5JetJIDn90Hits>n90hitsmin) {
//     jetidloose=true;
//   }
//   return jetidloose;
// }

// bool JetIdtight(double ak5JetJIDresEMF,double ak5JetJIDfHPD,int ak5JetJIDn90Hits, double ak5JetEta, double ak5JetPtRaw){
//   bool jetidtight=false;
//   bool jetidresEMF=true;
//   bool jetidfHPD_highPt=true;

//   double fhpdmax = 0.98;
//   double n90hitsmin =1;
//   double emf_min = 0.01;

//   if(fabs(ak5JetEta)<2.6 && ak5JetJIDresEMF<=emf_min) jetidresEMF=false;
//   if(fabs(ak5JetEta)<2.6 && ak5JetPtRaw>80 && ak5JetJIDresEMF>=1) jetidresEMF=false;
//   if(ak5JetPtRaw>25 && ak5JetJIDfHPD>=0.95) jetidfHPD_highPt=false;

//   if(jetidresEMF && jetidfHPD_highPt && ak5JetJIDfHPD<fhpdmax && ak5JetJIDn90Hits>n90hitsmin) 
//     {
//       jetidtight=true;
//     }
//   return jetidtight;
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
   
  ////////////////////// User's code to book histos - BEGIN ///////////////////////

  TH1F *h_Mlj_PAS = new TH1F ("h_Mlj_PAS","h_Mlj_PAS",200,0,2000);  h_Mlj_PAS->Sumw2();

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

  int jetAlgorithm = getPreCutValue1("jetAlgorithm");
  int metAlgorithm = getPreCutValue1("metAlgorithm");

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
	    //     JetPassID->push_back(
	    //  PFJetIdloose(PFJetChargedHadronEnergyFraction->at(ijet),
	    //       PFJetChargedEmEnergyFraction->at(ijet),
	    //       PFJetNeutralHadronEnergyFraction->at(ijet),
	    //       PFJetNeutralEmEnergyFraction->at(ijet),
	    //       PFJetEta->at(ijet) )
	    //  );
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
	    //     JetPassID->push_back(
	    //  JetIdloose(CaloJetresEMF->at(ijet),
	    //     CaloJetfHPD->at(ijet),
	    //     CaloJetn90Hits->at(ijet),
	    //     CaloJetEta->at(ijet) )
	    //  );
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

    if(metAlgorithm==1) // --> PFMET
      {
	thisMET = PFMET->at(0);
	thisMETPhi = PFMETPhi->at(0);
	thisSumET = PFSumET->at(0);
      }
    if(metAlgorithm==2) // --> CaloMET
      {
	thisMET = CaloMET->at(0);
	thisMETPhi = CaloMETPhi->at(0);
	thisSumET = CaloSumET->at(0);
      }
    if(metAlgorithm==3) // --> PFMET (with type-1 corrections)
      {
	thisMET = PFMETType1Cor->at(0);
	thisMETPhi = PFMETPhiType1Cor->at(0);
	thisSumET = PFSumETType1Cor->at(0);
      }


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
	for(int ijet=0; ijet<JetPt->size(); ijet++)
	  {
	    JetPt->at(ijet) *= JetEnergyScale;
	    JetEnergy->at(ijet) *= JetEnergyScale;
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


    // Jets
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
	if ( JetPt->at(ijet) < getPreCutValue1("jet_PtCut") ) continue;
	v_idx_jet_PtCut.push_back(ijet);
      }

// 	//Checking overlap between electrons and jets
// 	int JetOverlapsWithEle = 0; //don't change the default (0) 
// 	float minDeltaR=9999.;
// 	TVector3 jet_vec;
// 	jet_vec.SetPtEtaPhi(JetPt->at(ijet),JetEta->at(ijet),JetPhi->at(ijet));
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
// 	//if( ( JetOverlaps->at(ijet) & 1 << eleIDType) == 0)/* NO overlap with electrons */  
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
// 	    thisjet.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut[ijet]),
// 				 JetEta->at(v_idx_jet_PtCut[ijet]),
// 				 JetPhi->at(v_idx_jet_PtCut[ijet]),0);
// 	    STDOUT("CLEANING: j"<<ijet+1<<" Pt, eta, phi = " << ", "<<thisjet.Pt()<<", "<< thisjet.Eta() <<", "<< thisjet.Phi()<<" jetFlags="<<jetFlags[ijet] );
// 	  }
//       } // printouts for jet cleaning
    
    // Flagging jets if they overlap with muons
    vector <int> jetFlags_muon(v_idx_jet_PtCut.size(), 0);
    int Njetflagged_muon = 0;
    for (int imuon=0; imuon<v_idx_muon_PtCut_IDISO.size(); imuon++)
      {
	TLorentzVector muon;
        muon.SetPtEtaPhiM(MuonPt->at(v_idx_muon_PtCut_IDISO[imuon]),
			  MuonEta->at(v_idx_muon_PtCut_IDISO[imuon]),
			  MuonPhi->at(v_idx_muon_PtCut_IDISO[imuon]),0);
	TLorentzVector jet;
	double minDR=9999.;
	int ijet_minDR = -1;    
        for(int ijet=0; ijet<v_idx_jet_PtCut.size(); ijet++)
          {
	    if ( jetFlags_muon[ijet] == 1 ) 
	      continue;
            jet.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut[ijet]),
			     JetEta->at(v_idx_jet_PtCut[ijet]),
			     JetPhi->at(v_idx_jet_PtCut[ijet]),0);
	    double DR = jet.DeltaR(muon);
	    if (DR<minDR) 
	      {
		minDR = DR;
		ijet_minDR = ijet;
	      }
	  }
	if ( minDR < getPreCutValue1("jet_muon_DeltaRcut") && ijet_minDR > -1)
	  {
	    jetFlags_muon[ijet_minDR] = 1;
	    Njetflagged_muon++;
	  }
      }

    
    float JetEtaCutValue = getPreCutValue1("jet_EtaCut");
    for(int ijet=0; ijet<v_idx_jet_PtCut.size(); ijet++) //pT pre-cut + no overlaps with electrons + jetID
      {	
	bool passjetID = JetPassID->at(v_idx_jet_PtCut[ijet]);
	//bool passjetID = JetIdloose(CaloJetresEMF->at(v_idx_jet_PtCut[ijet]),CaloJetfHPD->at(v_idx_jet_PtCut[ijet]),CaloJetn90Hits->at(v_idx_jet_PtCut[ijet]), CaloJetEta->at(v_idx_jet_PtCut[ijet]));
	// ---- use the flag stored in rootTuples
	//if( (JetOverlaps->at(v_idx_jet_PtCut[ijet]) & 1 << eleIDType) == 0  /* NO overlap with electrons */  
	// ----

	if( jetFlags[ijet] == 0                            /* NO overlap with electrons */  
	    && jetFlags_muon[ijet] == 0 )                  /* NO overlap with muons */
	  //  && passjetID == true )                            /* pass JetID */
	  // && (caloJetOverlaps[ijet] & 1 << 5)==0 )         /* NO overlap with muons */      
	  v_idx_jet_PtCut_noOverlap.push_back(v_idx_jet_PtCut[ijet]);

	if( jetFlags[ijet] == 0                            /* NO overlap with electrons */  
	    && jetFlags_muon[ijet] == 0                    /* NO overlap with muons */
	    && passjetID == true )                            /* pass JetID */
	  // && (caloJetOverlaps[ijet] & 1 << 5)==0 )         /* NO overlap with muons */      
	  v_idx_jet_PtCut_noOverlap_ID.push_back(v_idx_jet_PtCut[ijet]);

	if( jetFlags[ijet] == 0                           /* NO overlap with electrons */  
	    && jetFlags_muon[ijet] == 0                   /* NO overlap with muons */
	    && passjetID == true                             /* pass JetID */
	    && fabs( JetEta->at(v_idx_jet_PtCut[ijet]) ) < JetEtaCutValue )
	  // && (caloJetOverlaps[ijet] & 1 << 5)==0 )         /* NO overlap with muons */      
	  v_idx_jet_PtCut_noOverlap_ID_EtaCut.push_back(v_idx_jet_PtCut[ijet]);

	
	//NOTE: We should verify that caloJetOverlaps match with the code above
      } // End loop over jets
    


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
    //    fillVariableWithValue( "PassHBHENoiseFilter", passLooseNoiseFilter ) ;

    // nMu
    fillVariableWithValue( "nMu_all", MuonPt->size() ) ;
    fillVariableWithValue( "nMu_PtCut", v_idx_muon_PtCut.size() ) ;
    fillVariableWithValue( "nMu_PtCut_IDISO", v_idx_muon_PtCut_IDISO.size() ) ;

    // nEle
    fillVariableWithValue( "nEle_all", v_idx_ele_all.size() ) ;
    fillVariableWithValue( "nEle_PtCut", v_idx_ele_PtCut.size() ) ;
    fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp", v_idx_ele_PtCut_IDISO_noOverlap.size() ) ;

    // nLeptons
    fillVariableWithValue( "nLept_PtCut_IDISO", v_idx_ele_PtCut_IDISO_noOverlap.size()+v_idx_muon_PtCut_IDISO.size() ) ;
     
    // nJet
    fillVariableWithValue( "nJet_all", v_idx_jet_all.size() ) ;
    fillVariableWithValue( "nJet_PtCut", v_idx_jet_PtCut.size() ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlp", v_idx_jet_PtCut_noOverlap.size() ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlp_ID", v_idx_jet_PtCut_noOverlap_ID.size() ) ;
    //OneEleOneMu
    fillVariableWithValue( "nJet_OneEleOneMu_All", v_idx_jet_PtCut_noOverlap_ID.size() ) ;
    fillVariableWithValue( "nJet_OneEleOneMu_EtaCut", v_idx_jet_PtCut_noOverlap_ID_EtaCut.size() ) ;
    //PAS June 2010
    fillVariableWithValue( "nJet_PAS_All", v_idx_jet_PtCut_noOverlap_ID.size() ) ;
    fillVariableWithValue( "nJet_PAS_EtaCut", v_idx_jet_PtCut_noOverlap_ID_EtaCut.size() ) ;

    // MET
    //PAS June 2010
    fillVariableWithValue( "MET_PAS", thisMET ) ;
    fillVariableWithValue( "METPhi_PAS", thisMETPhi);
    fillVariableWithValue( "pfMET_PAS", PFMET->at(0) ) ;
    fillVariableWithValue( "tcMET_PAS", TCMET->at(0) ) ;
    fillVariableWithValue( "caloMET_PAS", CaloMET->at(0) ) ;
    TVector2 v_MET;
    v_MET.SetMagPhi( thisMET , thisMETPhi);

    // 1st ele
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stEle_IDISO_NoOvrlp", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	fillVariableWithValue( "Eta1stEle_IDISO_NoOvrlp", ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	fillVariableWithValue( "mEta1stEle_IDISO_NoOvrlp", fabs(ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0])) );
	//PAS June 2010
	fillVariableWithValue( "Pt1stEle_PAS", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	fillVariableWithValue( "Eta1stEle_PAS", ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
        fillVariableWithValue( "Phi1stEle_PAS", ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	// transverse mass ele+MET
	TVector2 v_ele;
	v_ele.SetMagPhi( ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) , ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	double deltaphi = v_MET.DeltaPhi(v_ele);
	double MT_eleMET = sqrt(2 * v_ele.Mod() * thisMET * (1 - cos(deltaphi)) );
	fillVariableWithValue("MT_eleMET", MT_eleMET);
      }

    // 1st muon
    if( v_idx_muon_PtCut_IDISO.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stMu_IDISO", MuonPt->at(v_idx_muon_PtCut_IDISO[0]) );
	fillVariableWithValue( "Eta1stMu_IDISO", MuonEta->at(v_idx_muon_PtCut_IDISO[0]) );
	fillVariableWithValue( "mEta1stMu_IDISO", fabs(MuonEta->at(v_idx_muon_PtCut_IDISO[0])) );
	//PAS June 2010
	fillVariableWithValue( "Pt1stMu_PAS", MuonPt->at(v_idx_muon_PtCut_IDISO[0]) );
	fillVariableWithValue( "Eta1stMu_PAS", MuonEta->at(v_idx_muon_PtCut_IDISO[0]) );
        fillVariableWithValue( "Phi1stMu_PAS", MuonPhi->at(v_idx_muon_PtCut_IDISO[0]) );
	// transverse mass ele+MET
	TVector2 v_muon;
	v_muon.SetMagPhi( MuonPt->at(v_idx_muon_PtCut_IDISO[0]) , MuonPhi->at(v_idx_muon_PtCut_IDISO[0]) );
	double deltaphi = v_MET.DeltaPhi(v_muon);
	double MT_muonMET = sqrt(2 * v_muon.Mod() * thisMET * (1 - cos(deltaphi)) );
	fillVariableWithValue("MT_muonMET", MT_muonMET);
      }

    // 1st jet
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stJet_noOvrlp_ID", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "Eta1stJet_noOvrlp_ID", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "mEta1stJet_noOvrlp_ID", fabs(JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0])) );
	//PAS June 2010
	fillVariableWithValue( "Pt1stJet_PAS", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "Eta1stJet_PAS", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
        fillVariableWithValue( "Phi1stJet_PAS", JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
      }


    //cout << "2nd Jet" << endl;
    //## 2nd jet
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2 ) 
      {
	fillVariableWithValue( "Pt2ndJet_noOvrlp_ID", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "Eta2ndJet_noOvrlp_ID", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "mEta2ndJet_noOvrlp_ID", fabs(JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1])) );
	fillVariableWithValue( "maxMEtaJets_noOvrlp_ID", max( getVariableValue("mEta1stJet_noOvrlp_ID"), getVariableValue("mEta2ndJet_noOvrlp_ID") ) );
	//PAS June 2010
	fillVariableWithValue( "Pt2ndJet_PAS", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "Eta2ndJet_PAS", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
        fillVariableWithValue( "Phi2ndJet_PAS", JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
      }

    //## define "1ele+1muon" and "2jets" booleans
    bool OneEleOneMu=false;
    bool TwoJets=false;
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 1 && v_idx_muon_PtCut_IDISO.size() >= 1 ) OneEleOneMu = true;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2 ) TwoJets = true;

    // ST
    double calc_sT=-999.; 
    if ( (OneEleOneMu) && (TwoJets) ) 
      {
	calc_sT = 
	  ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) +
	  MuonPt->at(v_idx_muon_PtCut_IDISO[0]) +
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) +
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]);
	fillVariableWithValue("sT", calc_sT);
	fillVariableWithValue("sT_MLQ250", calc_sT);
	fillVariableWithValue("sT_MLQ280", calc_sT);
	fillVariableWithValue("sT_MLQ300", calc_sT);
	fillVariableWithValue("sT_MLQ320", calc_sT);
	fillVariableWithValue("sT_MLQ340", calc_sT);
	fillVariableWithValue("sT_MLQ400", calc_sT);       
	//PAS June 2010
	fillVariableWithValue("sT_PAS", calc_sT);
      }

    // ST jets
    if (TwoJets)
      {
	double calc_sTjet=-999.;
	calc_sTjet =
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) + 
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) ;
	fillVariableWithValue("sTjet_PAS", calc_sTjet);
      }

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
      }


    // Memu
    if (OneEleOneMu)
      {
	fillVariableWithValue( "maxMEtaEleMu_IDISO", max( getVariableValue("mEta1stEle_IDISO_NoOvrlp"), getVariableValue("mEta1stMu_IDISO") ) );
	TLorentzVector ele1, mu1, emu;
	ele1.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),0);
	mu1.SetPtEtaPhiM(MuonPt->at(v_idx_muon_PtCut_IDISO[0]),
			 MuonEta->at(v_idx_muon_PtCut_IDISO[0]),
			 MuonPhi->at(v_idx_muon_PtCut_IDISO[0]),0);
	emu = ele1+mu1;
	fillVariableWithValue("Memu", emu.M());
	//OneEleOneMu
	fillVariableWithValue("Memu_OneEleOneMu", emu.M());
	fillVariableWithValue("Memu_presel", emu.M());
	//
	//PAS June 2010
	fillVariableWithValue("Memu_PAS", emu.M());

	double calc_sTemu=-999.;
        calc_sTemu =
	   ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) +
	   MuonPt->at(v_idx_muon_PtCut_IDISO[0]) ;
	fillVariableWithValue("sTemu_PAS", calc_sTemu);

	if(isData==true) 
	  {
	    STDOUT("OneEleOneMu: Run, LS, Event = "<<run<<", "<<ls<<", "<<event);
	    STDOUT("OneEleOneMu: M_emu, Pt_emu, Eta_emu, Phi_emu = "<<emu.M() <<", "<< emu.Pt() <<", "<< emu.Eta() <<", "<< emu.Phi());
	    STDOUT("OneEleOneMu: 1st ele  Pt, eta, phi = "<< ele1.Pt() <<", "<< ele1.Eta() <<", "<< ele1.Phi() );
	    STDOUT("OneEleOneMu: 1st muon Pt, eta, phi = "<< mu1.Pt() <<", "<< mu1.Eta() <<", "<< mu1.Phi() );
	  }
      }

    // Mlj 
    double Me1j1, Me1j2, Mm1j1, Mm1j2 = -999;
    double deltaM_e1j1_m1j2 = 9999;
    double deltaM_e1j2_m1j1 = 9999;
    double Mlj_1stPair = 0;
    double Mlj_2ndPair = 0;
    double deltaR_e1j1 ;
    double deltaR_m1j2 ;
    double deltaR_e1j2 ;
    double deltaR_m1j1 ;
    if ( (OneEleOneMu) && (TwoJets) ) // OneEleOneMu and TwoJets
      {
	TLorentzVector jet1, jet2, ele1, mu1;
	ele1.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),0);
	mu1.SetPtEtaPhiM(MuonPt->at(v_idx_muon_PtCut_IDISO[0]),
			 MuonEta->at(v_idx_muon_PtCut_IDISO[0]),
			 MuonPhi->at(v_idx_muon_PtCut_IDISO[0]),0);
	jet1.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),0);
	jet2.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]),0);
	TLorentzVector e1j1, e1j2, m1j1, m1j2;
	e1j1 = ele1 + jet1;
	m1j2 =  mu1 + jet2;
	m1j1 =  mu1 + jet1;
	e1j2 = ele1 + jet2;
	Me1j1 = e1j1.M();
	Mm1j2 = m1j2.M();
	Me1j2 = e1j2.M();
	Mm1j1 = m1j1.M();

	deltaM_e1j1_m1j2 = Me1j1 - Mm1j2;
	deltaM_e1j2_m1j1 = Me1j2 - Mm1j1;


	double deltaR_e1j1 = ele1.DeltaR(jet1);
	double deltaR_m1j2 =  mu1.DeltaR(jet2);
	double deltaR_e1j2 = ele1.DeltaR(jet2);
	double deltaR_m1j1 =  mu1.DeltaR(jet1);

// 	// Fill min DR between any of the 2 selected eles and any of the 2 selected jets
// 	double minDR_2ele_2jet = min ( min(deltaR_e1j1,deltaR_m1j2) , min(deltaR_e1j2,deltaR_m1j1) );
// 	fillVariableWithValue("minDR_2ele_2jet", minDR_2ele_2jet);

	if(fabs(deltaM_e1j1_m1j2) > fabs(deltaM_e1j2_m1j1))
	  {
	    Mlj_1stPair = Me1j2;
	    Mlj_2ndPair = Mm1j1;
	    fillVariableWithValue("minDRlj_selecPairs", min(deltaR_e1j2,deltaR_m1j1) );
	    fillVariableWithValue("minDRlj_unselPairs", min(deltaR_e1j1,deltaR_m1j2) );
	  }
	else
	  {
	    Mlj_1stPair = Me1j1;
	    Mlj_2ndPair = Mm1j2;
	    fillVariableWithValue("minDRlj_selecPairs", min(deltaR_e1j1,deltaR_m1j2) );
	    fillVariableWithValue("minDRlj_unselPairs", min(deltaR_e1j2,deltaR_m1j1) );
	  } 
	fillVariableWithValue("Mlj_1stPair", Mlj_1stPair);       
	fillVariableWithValue("Mlj_2ndPair", Mlj_2ndPair);
	//PAS June 2010
	h_Mlj_PAS->Fill(Mlj_1stPair);
	h_Mlj_PAS->Fill(Mlj_2ndPair);
	fillVariableWithValue("Mlj_1stPair_PAS", Mlj_1stPair);       
	fillVariableWithValue("Mlj_2ndPair_PAS", Mlj_2ndPair);

	// min and max DeltaR between electrons and any jet
	double minDeltaR_lj = 999999;
	double maxDeltaR_lj = -1;
	double thisMinDR, thisMaxDR, DR_thisjet_e1, DR_thisjet_m1;
	TLorentzVector thisjet;
	for(int ijet=0; ijet<v_idx_jet_PtCut_noOverlap_ID.size(); ijet++)
	  {
	    thisjet.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
				 JetEta->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
				 JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),0);
	    DR_thisjet_e1 = thisjet.DeltaR(ele1);
	    DR_thisjet_m1 = thisjet.DeltaR(mu1);
	    thisMinDR = min(DR_thisjet_e1, DR_thisjet_m1);
	    thisMaxDR = max(DR_thisjet_e1, DR_thisjet_m1);
	    if(thisMinDR < minDeltaR_lj)
	      minDeltaR_lj = thisMinDR;
	    if(thisMaxDR > maxDeltaR_lj)
	      maxDeltaR_lj = thisMaxDR;
	  } 
	fillVariableWithValue("minDeltaR_lj", minDeltaR_lj);
	fillVariableWithValue("maxDeltaR_lj", maxDeltaR_lj);

	// printouts for small Mlj
	if(isData==true && ( Mlj_1stPair<20 || Mlj_2ndPair<20 ) ) // printouts for low Mlj 
	  {
	    STDOUT("Mlj<20GeV: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
	    STDOUT("Mlj<20GeV: Mlj_1stPair = "<<Mlj_1stPair <<", Mlj_2ndPair = "<< Mlj_2ndPair );
	    STDOUT("Mlj<20GeV: e1j1.M = "<<e1j1.M() <<", m1j2.M = "<<m1j2.M() <<", e1j2.M = "<<e1j2.M()  <<", m1j1.M = "<<m1j1.M()  );
	    STDOUT("Mlj<20GeV: deltaM_e1j1_m1j2 = "<<deltaM_e1j1_m1j2 <<", deltaM_e1j2_m1j1 = "<<deltaM_e1j2_m1j1  );
	    STDOUT("Mlj<20GeV: deltaR_e1j1 = "<<deltaR_e1j1 <<", deltaR_m1j2 = "<<deltaR_m1j2 <<", deltaR_e1j2 = "<<deltaR_e1j2  <<", deltaR_m1j1 = "<<deltaR_m1j1  );
	    TLorentzVector thisele;
	    // Add printout for muons?
	    for(int iele=0; iele<v_idx_ele_PtCut_IDISO_noOverlap.size(); iele++)
	      {
		thisele.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
				     ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
				     ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),0);
		STDOUT("Mlj<20GeV: e"<<iele+1<<" Pt, eta, phi = " 
		       << ", "<<thisele.Pt()<<", "<< thisele.Eta() <<", "<< thisele.Phi()<<"; DR_j1, DR_j2 = "<< thisele.DeltaR(jet1)<<", "<<thisele.DeltaR(jet2));
	      }
	    TLorentzVector thisjet;
	    TLorentzVector thisjet_e1, thisjet_m1;
	    for(int ijet=0; ijet<v_idx_jet_PtCut_noOverlap_ID.size(); ijet++)
	      {
		thisjet.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
				     JetEta->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
				     JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),0);
		thisjet_e1 = thisjet + ele1;
		thisjet_m1 = thisjet +  mu1;
		STDOUT("Mlj<20GeV: j"<<ijet+1<<" Pt, eta, phi = " << ", "<<thisjet.Pt()<<", "<< thisjet.Eta() <<", "<< thisjet.Phi()<<"; DR_e1, DR_m1 = "<< thisjet.DeltaR(ele1)<<", "<<thisjet.DeltaR(mu1) << "; M_e1, M_m1 = " <<thisjet_e1.M() <<", "<<thisjet_m1.M() );
	      }
	  } // printouts for low Mlj 

      } // OneEleOneMu and TwoJets



    // Evaluate cuts (but do not apply them)
    evaluateCuts();

    // Fill histograms and do analysis based on cut evaluation
    //h_nEleFinal->Fill( ElectronPt->size() );
     

    // printouts for OneEleOneMuTwoJets
    if( isData==true && passedAllPreviousCuts("Memu") ) 
      {
	STDOUT("OneEleOneMuTwoJets: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
	STDOUT("OneEleOneMuTwoJets: sT, Memu, Mjj_PAS = "<< getVariableValue("sT") <<", "<< getVariableValue("Memu")<<", "<< getVariableValue("Mjj_PAS") );
	STDOUT("OneEleOneMuTwoJets: Mlj_1stPair, Mlj_2ndPair = "<< getVariableValue("Mlj_1stPair")<<", "<< getVariableValue("Mlj_2ndPair") );
	STDOUT("OneEleOneMuTwoJets: 1stEle Pt, eta, phi = "<<getVariableValue("Pt1stEle_PAS") <<", "<< getVariableValue("Eta1stEle_PAS") <<", "<< getVariableValue("Phi1stEle_PAS") );
	STDOUT("OneEleOneMuTwoJets: 1stMu  Pt, eta, phi = "<<getVariableValue("Pt1stMu_PAS") <<", "<< getVariableValue("Eta1stMu_PAS") <<", "<< getVariableValue("Phi1stMu_PAS") );
	STDOUT("OneEleOneMuTwoJets: 1stJet Pt, eta, phi = "<<getVariableValue("Pt1stJet_PAS") <<", "<< getVariableValue("Eta1stJet_PAS") <<", "<< getVariableValue("Phi1stJet_PAS") );
	STDOUT("OneEleOneMuTwoJets: 2ndJet Pt, eta, phi = "<<getVariableValue("Pt2ndJet_PAS") <<", "<< getVariableValue("Eta2ndJet_PAS") <<", "<< getVariableValue("Phi2ndJet_PAS") );
	STDOUT("OneEleOneMuTwoJets: minDRlj_selecPairs, minDRlj_unselPairs = "<<getVariableValue("minDRlj_selecPairs") <<", "<< getVariableValue("minDRlj_unselPairs") );
// 	if ( passedCut("Mee") )
// 	  {
// 	    STDOUT("PassedMeeAndAllPrevious: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event<<", sT = "<< getVariableValue("sT"));
// 	  }
      } // printouts for OneEleOneMuTwoJets:

    
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

  h_Mlj_PAS->Write();

  ////////////////////// User's code to write histos - END ///////////////////////
  
  
  //STDOUT("analysisClass::Loop() ends");   
}
