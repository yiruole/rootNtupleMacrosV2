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
  double jet_ele_DeltaRcut = getPreCutValue1("jet_ele_DeltaRcut");
  double jet_PtCut_forMetScale = getPreCutValue1("jet_PtCut_forMetScale");

  double EleEnergyScale_EB=getPreCutValue1("EleEnergyScale_EB");
  double EleEnergyScale_EE=getPreCutValue1("EleEnergyScale_EE");
  double JetEnergyScale=getPreCutValue1("JetEnergyScale");
  int doJetOversmearing = getPreCutValue1("doJetOversmearing");
  double JetOversmearingSigma = getPreCutValue1("JetOversmearingSigma");

  // Not used when using ElectronHeepID and heepBitMask // int eleIDType = (int) getPreCutValue1("eleIDType");
  int heepBitMask_EB  =  getPreCutValue1("heepBitMask_EBGapEE") ;
  int heepBitMask_GAP =  getPreCutValue2("heepBitMask_EBGapEE") ;
  int heepBitMask_EE  =  getPreCutValue3("heepBitMask_EBGapEE") ;

  double muon_PtCut = getPreCutValue1("muon_PtCut");
  double muNHits_minThresh = getPreCutValue1("muNHits_minThresh");
  double muTrkD0Maximum = getPreCutValue1("muTrkD0Maximum");
  double jet_muon_DeltaRcut = getPreCutValue1("jet_muon_DeltaRcut");

  int jetAlgorithm = getPreCutValue1("jetAlgorithm");
  int metAlgorithm = getPreCutValue1("metAlgorithm");
  int eleAlgorithm = getPreCutValue1("eleAlgorithm");

  double vertexMinimumNDOF = getPreCutValue1("vertexMinimumNDOF");
  double vertexMaxAbsZ = getPreCutValue1("vertexMaxAbsZ");
  double vertexMaxd0 = getPreCutValue1("vertexMaxd0");
  
  ////////////////////// User's code to get preCut values - END /////////////////

  ////////////////////// User's code to book histos - BEGIN ///////////////////////

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

    //double event_weight = getPileupWeight ( PileUpInteractions, isData ) ;

    ////////////////////// User's code to be done for every event - BEGIN ///////////////////////

    //## Define new jet collection
    std::auto_ptr<std::vector<double> >  JetPt  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetPtRaw  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetEnergy  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetEta  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetPhi  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >     JetPassID  ( new std::vector<int>()  );
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
            JetPassID->push_back( PFJetPassLooseID->at(ijet) );
            JetTCHE->push_back( PFJetTrackCountingHighEffBTag->at(ijet) );
	    //JetChargedMuEnergyFraction->push_back( PFJetChargedMuEnergyFraction->at(ijet) );	    
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

    //## Muons
    vector<int> v_idx_muon_PtCut_IDISO;

    // Loop over muons
    for(int imuon=0; imuon<MuonPt->size(); imuon++){
      
      // pT pre-cut on muon
      if ( MuonPt->at(imuon) < muon_PtCut ) continue;

      if (   
	     ( MuonTrkHits->at(imuon)  >= muNHits_minThresh  )
	  && ( fabs( MuonTrkD0->at(imuon) ) < muTrkD0Maximum )
	  && ( MuonPassIso->at(imuon)==1 )
	  && ( MuonPassID->at(imuon)==1 ) 
	  )
	{
	  v_idx_muon_PtCut_IDISO.push_back(imuon);
	}

    }// end loop over muons

    //## Electrons
    vector<int> v_idx_ele_PtCut_IDISO;
    int heepBitMask;

    //Loop over electrons
    for(int iele=0; iele<ElectronPt->size(); iele++)
      {

	// pT pre-cut on ele
	if( ElectronPt->at(iele) < ele_PtCut ) continue;

	// get heepBitMask for EB, GAP, EE
	if( fabs(ElectronEta->at(iele)) < eleEta_bar )
	  {
	    heepBitMask = heepBitMask_EB;
	  }
	else if ( fabs(ElectronEta->at(iele)) > eleEta_end_min && fabs(ElectronEta->at(iele)) < eleEta_end_max )
	  {
	    heepBitMask = heepBitMask_EE;
	  }
	else 
	  {
	    heepBitMask = heepBitMask_GAP;
	  }

	//ID + ISO
	if ( (ElectronHeepID->at(iele) & ~heepBitMask)==0x0 )
	  {
	    v_idx_ele_PtCut_IDISO.push_back(iele);
	  }

      } // End loop over electrons

    //## Jets
    vector<int> v_idx_jet_PtCut;
    vector<int> v_idx_jet_PtCut_noOverlap_ID;
    vector<int> v_idx_jet_PtCut_noOverlap_ID_TCHEL;

    // Loop over jets
    for(int ijet=0; ijet<JetPt->size(); ijet++)
      {
	//pT pre-cut on reco jets
	if ( JetPt->at(ijet) < jet_PtCut ) continue;
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

    //Select final jets
    for(int ijet=0; ijet<v_idx_jet_PtCut.size(); ijet++) //pT pre-cut + no overlaps with electrons + jetID
      {
	bool passjetID = JetPassID->at(v_idx_jet_PtCut[ijet]);

	if( jetFlagsEle[ijet] == 0                                      /* NO overlap with electrons */ 
	    && jetFlagsMuon[ijet] == 0                                  /* NO overlap with muons */	    
	    && passjetID == true                                        /* pass JetID */
	    && fabs( JetEta->at(v_idx_jet_PtCut[ijet]) ) < jet_EtaCut 	/* pass Jet Eta cut */    
	    )                         
	  v_idx_jet_PtCut_noOverlap_ID.push_back(v_idx_jet_PtCut[ijet]);

	if( jetFlagsEle[ijet] == 0                                       /* NO overlap with electrons */ 
	    && jetFlagsMuon[ijet] == 0                                   /* NO overlap with muons */	    
	    && passjetID == true                                         /* pass JetID */
	    && fabs( JetEta->at(v_idx_jet_PtCut[ijet]) ) < jet_EtaCut 	 /* pass Jet Eta cut */    
            && fabs( JetTCHE->at(v_idx_jet_PtCut[ijet]) ) > jet_TCHELCut /* TrackCountingHighEfficiency loose b-tag (see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagPerformance) */
	    )                         
	  v_idx_jet_PtCut_noOverlap_ID_TCHEL.push_back(v_idx_jet_PtCut[ijet]);	
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

    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();
    //fillVariableWithValue( "PassJSON", passJSON(run, ls, isData) );    

    // Set the value of the variableNames listed in the cutFile to their current value

    //event info
    fillVariableWithValue( "isData", isData ) ;
    fillVariableWithValue( "bunch", bunch ) ;
    fillVariableWithValue( "event", event ) ;
    fillVariableWithValue( "ls", ls ) ;
    fillVariableWithValue( "orbit", orbit ) ;
    fillVariableWithValue( "run", run ) ;
    fillVariableWithValue( "ProcessID", ProcessID ) ;
    fillVariableWithValue( "PtHat", PtHat ) ;

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

    //Event filters at RECO level
    fillVariableWithValue( "PassBeamScraping", !isBeamScraping ) ;
    fillVariableWithValue( "PassPrimaryVertex", isPrimaryVertex ) ;
    fillVariableWithValue( "PassHBHENoiseFilter", passHBHENoiseFilter ) ;
    fillVariableWithValue( "PassBeamHaloFilterLoose", passBeamHaloFilterLoose ) ;
    fillVariableWithValue( "PassBeamHaloFilterTight", passBeamHaloFilterTight ) ;

    // nVertex and pile-up
    fillVariableWithValue( "nVertex", VertexChi2->size() ) ;
    fillVariableWithValue( "nVertex_good", v_idx_vertex_good.size() ) ;
    fillVariableWithValue( "nPileUpInteractions", PileUpInteractions ) ;
    
    // nEle
    fillVariableWithValue( "nEle", v_idx_ele_PtCut_IDISO.size() ) ;

    // nJet
    fillVariableWithValue( "nJet", v_idx_jet_PtCut_noOverlap_ID.size() ) ;
    fillVariableWithValue( "nJet_btagTCHE", v_idx_jet_PtCut_noOverlap_ID_TCHEL.size() ) ;

    // nMuon
    fillVariableWithValue( "nMuon", v_idx_muon_PtCut_IDISO.size() ) ;

    // MET
    fillVariableWithValue("MET_Pt", thisMET);
    fillVariableWithValue("MET_Phi", thisMETPhi);
    fillVariableWithValue("GenMET_Pt", thisGenMET);
    fillVariableWithValue("GenMET_Phi", thisGenMETPhi);          

    // SUMET
    fillVariableWithValue("SumET", thisSumET);
    fillVariableWithValue("GenSumET", thisGenSumET);

     // 1st ele 
    double MT_METEle1, mDeltaPhiMETEle1 = -999;
    if( v_idx_ele_PtCut_IDISO.size() >= 1 )
      {
	fillVariableWithValue( "Ele1_Pt", ElectronPt->at(v_idx_ele_PtCut_IDISO[0]) );
	fillVariableWithValue( "Ele1_Energy", ElectronEnergy->at(v_idx_ele_PtCut_IDISO[0]) );
	fillVariableWithValue( "Ele1_Eta", ElectronEta->at(v_idx_ele_PtCut_IDISO[0]) );
	fillVariableWithValue( "Ele1_Phi", ElectronPhi->at(v_idx_ele_PtCut_IDISO[0]) );
	fillVariableWithValue( "Ele1_Charge", ElectronCharge->at(v_idx_ele_PtCut_IDISO[0]) );

	// DeltaPhi - MET vs 1st ele
	TVector2 v_MET;
	TVector2 v_ele;
	v_MET.SetMagPhi( thisMET , thisMETPhi);
	v_ele.SetMagPhi( ElectronPt->at(v_idx_ele_PtCut_IDISO[0]) , ElectronPhi->at(v_idx_ele_PtCut_IDISO[0]) );
	float deltaphi = v_MET.DeltaPhi(v_ele);
	fillVariableWithValue( "mDPhi_METEle1", fabs(deltaphi) );
        mDeltaPhiMETEle1 = fabs(deltaphi);

	// transverse mass enu
	MT_METEle1 = sqrt(2 * ElectronPt->at(v_idx_ele_PtCut_IDISO[0]) * thisMET * (1 - cos(deltaphi)) );
	fillVariableWithValue("MT_Ele1MET", MT_METEle1);

	//PT(e,nu)
	TVector2 v_ele_MET;
	v_ele_MET = v_ele + v_MET;
	fillVariableWithValue("Pt_Ele1MET", v_ele_MET.Mod());
      }

    // 2nd ele 
    double MT_METEle2, mDeltaPhiMETEle2 = -999;
    if( v_idx_ele_PtCut_IDISO.size() >= 2 )
      {
	fillVariableWithValue( "Ele2_Pt", ElectronPt->at(v_idx_ele_PtCut_IDISO[1]) );
	fillVariableWithValue( "Ele2_Energy", ElectronEnergy->at(v_idx_ele_PtCut_IDISO[1]) );
	fillVariableWithValue( "Ele2_Eta", ElectronEta->at(v_idx_ele_PtCut_IDISO[1]) );
	fillVariableWithValue( "Ele2_Phi", ElectronPhi->at(v_idx_ele_PtCut_IDISO[1]) );
	fillVariableWithValue( "Ele2_Charge", ElectronCharge->at(v_idx_ele_PtCut_IDISO[1]) );

	// DeltaPhi - MET vs 2nd ele
	TVector2 v_MET;
	TVector2 v_ele;
	v_MET.SetMagPhi( thisMET , thisMETPhi);
	v_ele.SetMagPhi( ElectronPt->at(v_idx_ele_PtCut_IDISO[1]) , ElectronPhi->at(v_idx_ele_PtCut_IDISO[1]) );
	float deltaphi = v_MET.DeltaPhi(v_ele);
	fillVariableWithValue( "mDPhi_METEle2", fabs(deltaphi) );
        mDeltaPhiMETEle2 = fabs(deltaphi);

	// transverse mass enu
	//MT_METEle2 = sqrt(2 * ElectronPt->at(v_idx_ele_PtCut_IDISO[1]) * thisMET * (1 - cos(deltaphi)) );
	//fillVariableWithValue("MT_Ele2MET", MT_METEle2);

	//PT(e,nu)
	//TVector2 v_ele_MET;
	//v_ele_MET = v_ele + v_MET;
	//fillVariableWithValue("Pt_Ele2MET", v_ele_MET.Mod());
      }

    // 1st jet 
    double mDeltaPhiMET1stJet = -999;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 1 )        
      {
	fillVariableWithValue( "Jet1_Pt", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "Jet1_Energy", JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "Jet1_Eta", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "Jet1_Phi", JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
        fillVariableWithValue( "Jet1_btagTCHE", JetTCHE->at(v_idx_jet_PtCut_noOverlap_ID[0]) );

	//DeltaPhi - MET vs 1st jet
	TVector2 v_MET;
	TVector2 v_jet;
	v_MET.SetMagPhi( 1 , thisMETPhi );
	v_jet.SetMagPhi( 1 , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	float deltaphi = v_MET.DeltaPhi(v_jet);
	fillVariableWithValue( "mDPhi_METJet1", fabs(deltaphi) );
        mDeltaPhiMET1stJet = fabs(deltaphi);
      }

    // 2nd jet 
    double mDeltaPhiMET2ndJet = -999;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2 )        
      {
	fillVariableWithValue( "Jet2_Pt", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "Jet2_Energy", JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "Jet2_Eta", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "Jet2_Phi", JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
        fillVariableWithValue( "Jet2_btagTCHE", JetTCHE->at(v_idx_jet_PtCut_noOverlap_ID[1]) );

	//DeltaPhi - MET vs 2nd jet
	TVector2 v_MET;
	TVector2 v_jet;
	v_MET.SetMagPhi( 1 , thisMETPhi );
	v_jet.SetMagPhi( 1 , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	float deltaphi = v_MET.DeltaPhi(v_jet);
	fillVariableWithValue( "mDPhi_METJet2", fabs(deltaphi) );
        mDeltaPhiMET2ndJet = fabs(deltaphi);
      }

    // 3rd jet 
    double mDeltaPhiMET3rdJet = -999;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 3 )        
      {
	fillVariableWithValue( "Jet3_Pt", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[2]) );
	fillVariableWithValue( "Jet3_Energy", JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[2]) );
	fillVariableWithValue( "Jet3_Eta", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[2]) );
	fillVariableWithValue( "Jet3_Phi", JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[2]) );
        fillVariableWithValue( "Jet3_btagTCHE", JetTCHE->at(v_idx_jet_PtCut_noOverlap_ID[2]) );

	//DeltaPhi - MET vs 3rd jet
	TVector2 v_MET;
	TVector2 v_jet;
	v_MET.SetMagPhi( 1 , thisMETPhi );
	v_jet.SetMagPhi( 1 , JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[2]) );
	float deltaphi = v_MET.DeltaPhi(v_jet);
	fillVariableWithValue( "mDPhi_METJet3", fabs(deltaphi) );
        mDeltaPhiMET3rdJet = fabs(deltaphi);
      }

    // 1st muon 
    double MT_METMuon1, mDeltaPhiMETMuon1 = -999;
    if( v_idx_muon_PtCut_IDISO.size() >= 1 )
      {
	fillVariableWithValue( "Muon1_Pt", MuonPt->at(v_idx_muon_PtCut_IDISO[0]) );
	fillVariableWithValue( "Muon1_Energy", MuonEnergy->at(v_idx_muon_PtCut_IDISO[0]) );
	fillVariableWithValue( "Muon1_Eta", MuonEta->at(v_idx_muon_PtCut_IDISO[0]) );
	fillVariableWithValue( "Muon1_Phi", MuonPhi->at(v_idx_muon_PtCut_IDISO[0]) );
	fillVariableWithValue( "Muon1_Charge", MuonCharge->at(v_idx_muon_PtCut_IDISO[0]) );
      }

    // define booleans
    bool OneEle=false;
    bool TwoEles=false;
    bool OneJet=false;
    bool TwoJets=false;
    bool ThreeJets=false;
    if( v_idx_ele_PtCut_IDISO.size() >= 1 ) OneEle = true;
    if( v_idx_ele_PtCut_IDISO.size() >= 2 ) TwoEles = true;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 1 ) OneJet = true;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2 ) TwoJets = true;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 3 ) ThreeJets = true;

    //Delta R
    if ( (OneEle) && (OneJet) )
      {	
	TLorentzVector ele1;
	ele1.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO[0]),
			  ElectronEta->at(v_idx_ele_PtCut_IDISO[0]),
			  ElectronPhi->at(v_idx_ele_PtCut_IDISO[0]),0);
	TLorentzVector jet1;
	jet1.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
       
	fillVariableWithValue("DR_Ele1Jet1", ele1.DeltaR(jet1));
	
	if ( (TwoJets) )
	  {
	    TLorentzVector jet2;
	    jet2.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			      JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			      JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			      JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	    
	    fillVariableWithValue("DR_Ele1Jet2", ele1.DeltaR(jet2));	   
	  }

	if( (TwoEles) )
	  {
	    TLorentzVector ele2;
	    ele2.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO[1]),
			      ElectronEta->at(v_idx_ele_PtCut_IDISO[1]),
			      ElectronPhi->at(v_idx_ele_PtCut_IDISO[1]),0);	    

	    fillVariableWithValue("DR_Ele2Jet1", ele2.DeltaR(jet1));	   	    
	  }

	if( (TwoEles) && (TwoJets) )
	  {
	    TLorentzVector ele2;
	    ele2.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO[1]),
			      ElectronEta->at(v_idx_ele_PtCut_IDISO[1]),
			      ElectronPhi->at(v_idx_ele_PtCut_IDISO[1]),0);	    
	    TLorentzVector jet2;
	    jet2.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			      JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			      JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			      JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[1]) );

	    fillVariableWithValue("DR_Ele2Jet2", ele2.DeltaR(jet2));	   	    
	  }       	
      }

    // Mjj
    if (TwoJets)
      {
	TLorentzVector jet1, jet2, jj;
	jet1.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	jet2.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	jj = jet1+jet2;

	fillVariableWithValue("M_j1j2", jj.M());
	fillVariableWithValue("DR_Jet1Jet2", jet1.DeltaR(jet2));
	fillVariableWithValue("Pt_j1j2", jj.Pt());
      }

    if (ThreeJets)
      {
	TLorentzVector jet1, jet2, jet3, j1j3, j2j3;
	jet1.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	jet2.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	jet3.SetPtEtaPhiE(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[2]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[2]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[2]),
			  JetEnergy->at(v_idx_jet_PtCut_noOverlap_ID[2]) );
	j1j3 = jet1+jet3;
	j2j3 = jet2+jet3;

	fillVariableWithValue("M_j1j3", j1j3.M());
	fillVariableWithValue("M_j2j3", j2j3.M());
      }

    // ST enujj
    if ( (OneEle) && (TwoJets) )
      {
	double sT_enujj =
	  ElectronPt->at(v_idx_ele_PtCut_IDISO[0]) +
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) +
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) +
	  thisMET;
	fillVariableWithValue("sT_enujj", sT_enujj);
      }
       
    // Evaluate cuts (but do not apply them)
    evaluateCuts();

    // Fill histograms and do analysis based on cut evaluation

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
    if( passedAllPreviousCuts("PassHBHENoiseFilter") 
	&& passedCut("nEle")
	&& passedCut("Ele1_Pt")
	&& passedCut("MET_Pt")
	&& passedCut("Jet1_Pt")
	&& passedCut("Jet2_Pt")	
	) 
      fillSkimTree();

    // Produce reduced skim
    if( passedAllPreviousCuts("PassHBHENoiseFilter") 
	&& passedCut("nEle")
	&& passedCut("Ele1_Pt")
	&& passedCut("MET_Pt")
	//&& passedCut("Jet1_Pt")
	//&& passedCut("Jet2_Pt")	
	) 
      fillReducedSkimTree();
    
    ////////////////////// User's code to be done for every event - END ///////////////////////

  } // End of loop over events


  ////////////////////// User's code to write histos - BEGIN ///////////////////////

  ////////////////////// User's code to write histos - END ///////////////////////

  delete randomNumGen;

  //STDOUT("analysisClass::Loop() ends");
}
