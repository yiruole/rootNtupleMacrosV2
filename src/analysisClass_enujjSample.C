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

bool JetIdtight(double ak5JetJIDresEMF,double ak5JetJIDfHPD,int ak5JetJIDn90Hits, double ak5JetEta, double ak5JetPt){
  bool jetidtight=false;
  bool jetidresEMF=true;
  bool jetidfHPD_highPt=true;

  double fhpdmax = 0.98;
  double n90hitsmin =1;
  double emf_min = 0.01;

  if(fabs(ak5JetEta)<2.6 && ak5JetJIDresEMF<=emf_min) jetidresEMF=false;
  if(fabs(ak5JetEta)<2.6 && ak5JetPt>80 && ak5JetJIDresEMF>=1) jetidresEMF=false;
  if(ak5JetPt>25 && ak5JetJIDfHPD>=0.95) jetidfHPD_highPt=false;

  if(jetidresEMF && jetidfHPD_highPt && ak5JetJIDfHPD<fhpdmax && ak5JetJIDn90Hits>n90hitsmin) 
    {
      jetidtight=true;
    }
  return jetidtight;
}

bool PFJetIdloose(const double ak5ChargedHadronFraction, const double ak5ChargedEmFraction, const double ak5NeutralHadronFraction, const double ak5NeutralEmFraction, const double ak5JetEta){
  bool jetidloose=false;
  bool jetidChHadFrac=true;

  double chHadFrac = 0.;
  double chEmFrac = 1.;
  double neutHadFrac = 1.;
  double neutEmFrac = 1.;

  if(fabs(ak5JetEta)<2.4 && ak5ChargedHadronFraction<=chHadFrac) jetidChHadFrac=false;

  if(jetidChHadFrac && ak5ChargedEmFraction<chEmFrac && ak5NeutralHadronFraction<neutHadFrac && ak5NeutralEmFraction<neutEmFrac) {
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
   
  ////////////////////// User's code to book histos - BEGIN ///////////////////////

  TH1F *h_Mej = new TH1F ("h_Mej","h_Mej",200,0,2000);
  h_Mej->Sumw2();

  TH1F *h_MTnuj = new TH1F ("h_MTnuj","h_MTnuj",200,0,2000);
  h_MTnuj->Sumw2();
  
  TH2D *h2_MTnuj_vs_MET = new TH2D ("h2_MTnuj_vs_MET","h2_MTnuj_vs_MET;#slash{E}_{T} [GeV];M_{T}(#nu,j) [GeV]",200,0,1000,200,0,1000);
  h2_MTnuj_vs_MET->Sumw2();
  
  TH2D *h2_ST_vs_MET = new TH2D ("h2_ST_vs_MET","h2_ST_vs_MET;#slash{E}_{T} [GeV];S_{T} [GeV]",200,0,1000,200,0,2000);
  h2_ST_vs_MET->Sumw2();
  
  TH2D *h2_ST_vs_MTnuj = new TH2D ("h2_ST_vs_MTnuj","h2_ST_vs_MTnuj;M_{T}(#nu,j) [GeV];S_{T} [GeV]",200,0,1000,200,0,2000);
  h2_ST_vs_MTnuj->Sumw2();
  
  TH2D *h2_DeltaPhiMETEle_vs_MET = new TH2D ("h2_DeltaPhiMETEle_vs_MET","h2_DeltaPhiMETEle_vs_MET;#slash{E}_{T} [GeV];#Delta#phi(#slash{E}_{T},e)",200,0,1000,50,0,4);
  h2_DeltaPhiMETEle_vs_MET->Sumw2();
  
  TH2D *h2_DeltaPhiMET1stJet_vs_MET = new TH2D ("h2_DeltaPhiMET1stJet_vs_MET","h2_DeltaPhiMET1stJet_vs_MET;#slash{E}_{T} [GeV];#Delta#phi(#slash{E}_{T},j1)",200,0,1000,50,0,4);
  h2_DeltaPhiMET1stJet_vs_MET->Sumw2();
  
  TH2D *h2_DeltaPhiMET2ndJet_vs_MET = new TH2D ("h2_DeltaPhiMET2ndJet_vs_MET","h2_DeltaPhiMET2ndJet_vs_MET;#slash{E}_{T} [GeV];#Delta#phi(#slash{E}_{T},j2)",200,0,1000,50,0,4);
  h2_DeltaPhiMET2ndJet_vs_MET->Sumw2();
  
  TH2D *h2_DeltaPhiMETEle_vs_MET1stJet = new TH2D ("h2_DeltaPhiMETEle_vs_MET1stJet","h2_DeltaPhiMETEle_vs_MET1stJet;#Delta#phi(#slash{E}_{T},j1);#Delta#phi(#slash{E}_{T},e)",50,0,4,50,0,4);
  h2_DeltaPhiMETEle_vs_MET1stJet->Sumw2();
  
  TH2D *h2_DeltaPhiMETEle_vs_MET2ndJet = new TH2D ("h2_DeltaPhiMETEle_vs_MET2ndJet","h2_DeltaPhiMETEle_vs_MET2ndJet;#Delta#phi(#slash{E}_{T},j2);#Delta#phi(#slash{E}_{T},e)",50,0,4,50,0,4);
  h2_DeltaPhiMETEle_vs_MET2ndJet->Sumw2();
  
  TH2D *h2_DeltaPhiMET2ndJet_vs_MET1stJet = new TH2D ("h2_DeltaPhiMET2ndJet_vs_MET1stJet","h2_DeltaPhiMET2ndJet_vs_MET1stJet;#Delta#phi(#slash{E}_{T},j1);#Delta#phi(#slash{E}_{T},j2)",50,0,4,50,0,4);
  h2_DeltaPhiMET2ndJet_vs_MET1stJet->Sumw2();
  ////////////////////// User's code to book histos - END ///////////////////////
    
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

    //## HLT
    int HLTTrigger = (int) getPreCutValue1("HLTTrigger");
    int PassTrig=HLTResults->at(HLTTrigger); // results of HLTPhoton15 

    // Electrons
    vector<int> v_idx_ele_all;
    vector<int> v_idx_ele_PtCut;
    vector<int> v_idx_ele_PtCut_IDISO_noOverlap;
    int eleIDType = (int) getPreCutValue1("eleIDType");
     
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

	//ID + ISO + NO overlap with good muons 
	int eleID = ElectronPassID->at(iele);
	if ( (eleID & 1<<eleIDType) > 0  && ElectronOverlaps->at(iele)==0 )
	  {
	    v_idx_ele_PtCut_IDISO_noOverlap.push_back(iele);
	  }

      } // End loop over electrons

    
    // Jets
    vector<int> v_idx_jet_all;
    vector<int> v_idx_jet_PtCut;
    vector<int> v_idx_jet_PtCut_noOverlap;
    vector<int> v_idx_jet_PtCut_noOverlap_ID;

    // Loop over jets
    for(int ijet=0; ijet<CaloJetPt->size(); ijet++)
      {

	//no cut on reco jets
	v_idx_jet_all.push_back(ijet);

	//pT pre-cut on reco jets
	if ( CaloJetPt->at(ijet) < getPreCutValue1("jet_PtCut") ) continue;
	v_idx_jet_PtCut.push_back(ijet);

        //Checking overlap between electrons and jets
        int JetOverlapsWithEle = 0; //don't change the default (0) 
        float minDeltaR=9999.;
        TVector3 jet_vec;
        jet_vec.SetPtEtaPhi(CaloJetPt->at(ijet),CaloJetEta->at(ijet),CaloJetPhi->at(ijet));
        for (int i=0; i < v_idx_ele_PtCut_IDISO_noOverlap.size(); i++){
          TVector3 ele_vec;       
          ele_vec.SetPtEtaPhi(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[i])
                              ,ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[i])
                              ,ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[i]));
          double distance = jet_vec.DeltaR(ele_vec);
          if (distance<minDeltaR) minDeltaR=distance;
        }
        if ( minDeltaR < getPreCutValue1("jet_ele_DeltaRcut") )
          JetOverlapsWithEle = 1; //this jet overlaps with a good electron --> remove it from the analysis

        //pT pre-cut + no overlaps with electrons
        // ---- use the flag stored in rootTuples
        //if( ( CaloJetOverlaps->at(ijet) & 1 << eleIDType) == 0)/* NO overlap with electrons */  
        // && (caloJetOverlaps[ijet] & 1 << 5)==0 )/* NO overlap with muons */   
        // ----
        if( JetOverlapsWithEle == 0 )  /* NO overlap with electrons */  
          v_idx_jet_PtCut_noOverlap.push_back(ijet);

        //pT pre-cut + no overlaps with electrons + jetID
        bool passjetID = JetIdloose(CaloJetresEMF->at(ijet),CaloJetfHPD->at(ijet),CaloJetn90Hits->at(ijet), CaloJetEta->at(ijet));
        //bool passjetID = PFJetIdloose(PFJetChargedHadronEnergyFraction->at(ijet),PFJetChargedEmEnergyFraction->at(ijet),PFJetNeutralHadronEnergyFraction->at(ijet),PFJetNeutralEmEnergyFraction->at(ijet),PFJetEta->at(ijet));
        // ---- use the flag stored in rootTuples
        //if( (CaloJetOverlaps->at(ijet) & 1 << eleIDType) == 0  /* NO overlap with electrons */  
        // ----
        if( JetOverlapsWithEle == 0                           /* NO overlap with electrons */  
            && passjetID == true )                            /* pass JetID */
          // && (caloJetOverlaps[ijet] & 1 << 5)==0 )         /* NO overlap with muons */      
          v_idx_jet_PtCut_noOverlap_ID.push_back(ijet);

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
    fillVariableWithValue( "nEle_all", v_idx_ele_all.size() ) ;
    fillVariableWithValue( "nEle_PtCut", v_idx_ele_PtCut.size() ) ;
    fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp", v_idx_ele_PtCut_IDISO_noOverlap.size() ) ;
     
    // nJet
    fillVariableWithValue( "nJet_all", v_idx_jet_all.size() ) ;
    fillVariableWithValue( "nJet_PtCut", v_idx_jet_PtCut.size() ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlp", v_idx_jet_PtCut_noOverlap.size() ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlp_ID", v_idx_jet_PtCut_noOverlap_ID.size() ) ;
    //PAS June 2010
    fillVariableWithValue( "nJet_PAS", v_idx_jet_PtCut_noOverlap_ID.size() ) ;

    // nMuon
    fillVariableWithValue( "nMuon_PtCut_IDISO", v_idx_muon_PtCut_IDISO.size() ) ;
    //PAS June 2010
    fillVariableWithValue( "nMuon_PtCut_IDISO_PAS", v_idx_muon_PtCut_IDISO.size() ) ;

    // MET
    // --> CaloMET
//     double thisMET = CaloMET->at(0);
//     double thisMETPhi = CaloMETPhi->at(0);
    //#
    // --> TCMET
//     double thisMET = TCMET->at(0);
//     double thisMETPhi = TCMETPhi->at(0);
    //#
    // --> PFMET
    double thisMET = PFMET->at(0);
    double thisMETPhi = PFMETPhi->at(0);
    //#
    fillVariableWithValue("MET", thisMET);

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

	// DeltaPhi - MET vs 1st ele
	TVector2 v_MET;
	TVector2 v_ele;
	v_MET.SetMagPhi( 1 , thisMETPhi);
	v_ele.SetMagPhi( 1 , ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) ); 
	float deltaphi = v_MET.DeltaPhi(v_ele);
	fillVariableWithValue( "mDeltaPhiMETEle", fabs(deltaphi) );
	//PAS June 2010
	fillVariableWithValue( "mDeltaPhiMETEle_PAS", fabs(deltaphi) );
        DeltaPhiMETEle = fabs(deltaphi);
        h2_DeltaPhiMETEle_vs_MET->Fill(thisMET, fabs(deltaphi) );

	// transverse mass enu
	MT = sqrt(2 * ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) * thisMET * (1 - cos(deltaphi)) );
	fillVariableWithValue("MTenu", MT);
        h2_MTnuj_vs_MET->Fill(thisMET,MT);
	//PAS June 2010
	fillVariableWithValue("MTenu_PAS", MT);
      }


    // 1st jet and deltaphi jet-MET
    double DeltaPhiMET1stJet = -999;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stJet_noOvrlp_ID", CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "Eta1stJet_noOvrlp_ID", CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "mEta1stJet_noOvrlp_ID", fabs(CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[0])) );
	//PAS June 2010
	fillVariableWithValue( "Pt1stJet_PAS", CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "Eta1stJet_PAS", CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]) );

	//DeltaPhi - MET vs 1st jet
	TVector2 v_MET;
	TVector2 v_jet;
	v_MET.SetMagPhi( 1 , thisMETPhi );
	v_jet.SetMagPhi( 1 , CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	float deltaphi = v_MET.DeltaPhi(v_jet);
	fillVariableWithValue( "mDeltaPhiMET1stJet", fabs(deltaphi) );
	//PAS June 2010
	fillVariableWithValue( "mDeltaPhiMET1stJet_PAS", fabs(deltaphi) );
        DeltaPhiMET1stJet = fabs(deltaphi);
        h2_DeltaPhiMET1stJet_vs_MET->Fill(thisMET, fabs(deltaphi) );
        if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 1 ) h2_DeltaPhiMETEle_vs_MET1stJet->Fill( fabs(deltaphi), DeltaPhiMETEle);
      }


    // 2nd jet and deltaphi jet-MET
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2 ) 
      {
	fillVariableWithValue( "Pt2ndJet_noOvrlp_ID", CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "Eta2ndJet_noOvrlp_ID", CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "mEta2ndJet_noOvrlp_ID", fabs(CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[1])) );
	fillVariableWithValue( "maxMEtaJets_noOvrlp_ID", max( getVariableValue("mEta1stJet_noOvrlp_ID"), getVariableValue("mEta2ndJet_noOvrlp_ID") ) );
	//PAS June 2010
	fillVariableWithValue( "Pt2ndJet_PAS", CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "Eta2ndJet_PAS", CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]) );

	//DeltaPhi - MET vs 2nd jet
	TVector2 v_MET;
	TVector2 v_jet;
	v_MET.SetMagPhi( 1 , thisMETPhi );
	v_jet.SetMagPhi( 1 , CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	float deltaphi = v_MET.DeltaPhi(v_jet);
	fillVariableWithValue( "mDeltaPhiMET2ndJet", fabs(deltaphi) );
	//PAS June 2010
	fillVariableWithValue( "mDeltaPhiMET2ndJet_PAS", fabs(deltaphi) );
        h2_DeltaPhiMET2ndJet_vs_MET->Fill(thisMET, fabs(deltaphi) );
        h2_DeltaPhiMET2ndJet_vs_MET1stJet->Fill( DeltaPhiMET1stJet, fabs(deltaphi) );
        if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 1 ) h2_DeltaPhiMETEle_vs_MET2ndJet->Fill( fabs(deltaphi), DeltaPhiMETEle);
      }

    // define "1ele" and "2jets" booleans
    bool OneEle=false;
    bool TwoJets=false;
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() == 1 ) OneEle = true;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2 ) TwoJets = true;


    // ST
    if ( (OneEle) && (TwoJets) ) 
      {
	double calc_sT = 
	  ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) +
	  CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) +
	  CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) + 
	  thisMET;
	fillVariableWithValue("sT", calc_sT);
        h2_ST_vs_MET->Fill(thisMET,calc_sT);
        h2_ST_vs_MTnuj->Fill(MT,calc_sT);
	//PAS June 2010
	fillVariableWithValue("sT_PAS", calc_sT);
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
	jet1.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),0);
	jet2.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  CaloJetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]),0);
	TLorentzVector jet1ele1, jet2ele1;
	jet1ele1 = jet1 + ele1;
	jet2ele1 = jet2 + ele1;
	Me1j1 = jet1ele1.M();
	Me1j2 = jet2ele1.M();

	//transverse mass neutrino-jet
	TVector2 v_MET;
	TVector2 v_jet1;
	TVector2 v_jet2;
	v_MET.SetMagPhi( 1 , thisMETPhi);
	v_jet1.SetMagPhi(1 , CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]));
	v_jet2.SetMagPhi(1 , CaloJetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]));
	float deltaphi1 = v_MET.DeltaPhi(v_jet1);
	float deltaphi2 = v_MET.DeltaPhi(v_jet2);
	MTn1j1 = sqrt(2 * thisMET * CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) * (1 - cos(deltaphi1)) );
	MTn1j2 = sqrt(2 * thisMET * CaloJetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) * (1 - cos(deltaphi2)) );
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

    if( passedCut("all") )
      {
	//Mej
	h_Mej->Fill(Me1j1);
	h_Mej->Fill(Me1j2);
	//MTnuj
	h_MTnuj->Fill(MTn1j1);
	h_MTnuj->Fill(MTn1j2);
      }
     
    ////////////////////// User's code to be done for every event - END ///////////////////////
    
  } // End of loop over events
  

  ////////////////////// User's code to write histos - BEGIN ///////////////////////

  h_Mej->Write();
  h_MTnuj->Write();
  h2_MTnuj_vs_MET->Write();
  h2_ST_vs_MET->Write();
  h2_ST_vs_MTnuj->Write();
  h2_DeltaPhiMETEle_vs_MET->Write();
  h2_DeltaPhiMET1stJet_vs_MET->Write();
  h2_DeltaPhiMET2ndJet_vs_MET->Write();
  h2_DeltaPhiMETEle_vs_MET1stJet->Write();
  h2_DeltaPhiMETEle_vs_MET2ndJet->Write();
  h2_DeltaPhiMET2ndJet_vs_MET1stJet->Write();

  ////////////////////// User's code to write histos - END ///////////////////////
  
  
  //STDOUT("analysisClass::Loop() ends");   
}
