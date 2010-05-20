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
   
  ////////////////////// User's code to book histos - BEGIN ///////////////////////

  //barrel - all   
  TH1F *h_ElectronPt_barrel_all = new TH1F ("h_ElectronPt_barrel_all","h_ElectronPt_barrel_all",100,0,500); 
  h_ElectronPt_barrel_all->Sumw2();
  //barrel - heep   
  TH1F *h_ElectronPt_barrel_heep = new TH1F ("h_ElectronPt_barrel_heep","h_ElectronPt_barrel_heep",100,0,500);
  h_ElectronPt_barrel_heep->Sumw2();
  
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
    vector<int> v_idx_jet_PtCut_noOverlapEle;

    // Loop over jets
    for(int ijet=0; ijet<CaloJetPt->size(); ijet++)
      {

	//no cut on reco jets
	v_idx_jet_all.push_back(ijet);

	//pT pre-cut on reco jets
	if ( CaloJetPt->at(ijet) < getPreCutValue1("jet_PtCut") ) continue;

	v_idx_jet_PtCut.push_back(ijet);

	//  //Disambiguation of electrons from jets
	//  float minDeltaR=9999;
	//  TVector3 jet_vec;
	//  jet_vec.SetPtEtaPhi(caloJetPt[ijet],caloJetEta[ijet],caloJetPhi[ijet]);
	//  for (int i=0; i < v_idx_ele_final.size(); i++){
	//    TVector3 ele_vec;
	//    ele_vec.SetPtEtaPhi(elePt[v_idx_ele_final[i]],eleEta[v_idx_ele_final[i]],elePhi[v_idx_ele_final[i]]);
	//    double distance = jet_vec.DeltaR(ele_vec);
	//    if (distance<minDeltaR) minDeltaR=distance;
	//   }
	//  if ( minDeltaR > deltaR_minCut )  v_idx_jet_final.push_back(ijet);

	//pT pre-cut + no overlaps with electrons
	if( ( CaloJetOverlaps->at(ijet) & 1 << eleIDType) == 0)/* NO overlap with electrons */  
	  // && (caloJetOverlaps[ijet] & 1 << 5)==0 )/* NO overlap with muons */   
	  v_idx_jet_PtCut_noOverlapEle.push_back(ijet);

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
	if ( MuonPassIso->at(imuon)==1 && MuonPassID->at(imuon)==1)
	  v_idx_muon_PtCut_IDISO.push_back(imuon);

	//  //barrel-endcap definition
	//  bool in_Barrel=false;
	//          bool in_Endcap=false;
	//  if( fabs(muonEta[imuon]) < getPreCutValue1("muonEta_bar") )  in_Barrel=true;
	//  if( ( fabs(muonEta[imuon]) > getPreCutValue1("muonEta_end") )
	//      && 
	//      ( fabs(muonEta[imuon]) < getPreCutValue2("muonEta_end") ) 
	//      )  in_Endcap=true;
	 
      } // End loop over muons


    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();
    

    // Set the value of the variableNames listed in the cutFile to their current value
    
    // HLT
    //fillVariableWithValue( "HLT", PassTrig ) ;

    // nEle
    fillVariableWithValue( "nEle_all", v_idx_ele_all.size() ) ;
    fillVariableWithValue( "nEle_PtCut", v_idx_ele_PtCut.size() ) ;
    fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp", v_idx_ele_PtCut_IDISO_noOverlap.size() ) ;
     
    // nJet
    fillVariableWithValue( "nJet_all", v_idx_jet_all.size() ) ;
    fillVariableWithValue( "nJet_PtCut", v_idx_jet_PtCut.size() ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlpEle", v_idx_jet_PtCut_noOverlapEle.size() ) ;

    // 1st ele
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stEle_IDISO_NoOvrlp", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	fillVariableWithValue( "Eta1stEle_IDISO_NoOvrlp", ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	fillVariableWithValue( "mEta1stEle_IDISO_NoOvrlp", fabs(ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0])) );
      }

    // 2nd ele
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 2 ) 
      {
	fillVariableWithValue( "Pt2ndEle_IDISO_NoOvrlp", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) );
	fillVariableWithValue( "Eta2ndEle_IDISO_NoOvrlp", ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) );
	fillVariableWithValue( "mEta2ndEle_IDISO_NoOvrlp", fabs(ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[1])) );
	fillVariableWithValue( "maxMEtaEles_IDISO_NoOvrl", max( getVariableValue("mEta1stEle_IDISO_NoOvrlp"), getVariableValue("mEta2ndEle_IDISO_NoOvrlp") ) );
      }

    // 1st jet
    if( v_idx_jet_PtCut_noOverlapEle.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stJet_noOvrlpEle", CaloJetPt->at(v_idx_jet_PtCut_noOverlapEle[0]) );
	fillVariableWithValue( "Eta1stJet_noOvrlpEle", CaloJetEta->at(v_idx_jet_PtCut_noOverlapEle[0]) );
	fillVariableWithValue( "mEta1stJet_noOvrlpEle", fabs(CaloJetEta->at(v_idx_jet_PtCut_noOverlapEle[0])) );
      }


    //cout << "2nd Jet" << endl;
    //## 2nd jet
    if( v_idx_jet_PtCut_noOverlapEle.size() >= 2 ) 
      {
	fillVariableWithValue( "Pt2ndJet_noOvrlpEle", CaloJetPt->at(v_idx_jet_PtCut_noOverlapEle[1]) );
	fillVariableWithValue( "Eta2ndJet_noOvrlpEle", CaloJetEta->at(v_idx_jet_PtCut_noOverlapEle[1]) );
	fillVariableWithValue( "mEta2ndJet_noOvrlpEle", fabs(CaloJetEta->at(v_idx_jet_PtCut_noOverlapEle[1])) );
	fillVariableWithValue( "maxMEtaJets_noOvrlpEle", max( getVariableValue("mEta1stJet_noOvrlpEle"), getVariableValue("mEta2ndJet_noOvrlpEle") ) );
      }

    //## define "2ele" and "2jets" booleans
    bool TwoEle=false;
    bool TwoJets=false;
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 2 ) TwoEle = true;
    if( v_idx_jet_PtCut_noOverlapEle.size() >= 2 ) TwoJets = true;

    // ST
    double calc_sT=-999.; 
    if ( (TwoEle) && (TwoJets) ) 
      {
	calc_sT = 
	  ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) +
	  ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) +
	  CaloJetPt->at(v_idx_jet_PtCut_noOverlapEle[0]) +
	  CaloJetPt->at(v_idx_jet_PtCut_noOverlapEle[1]);
	fillVariableWithValue("sT", calc_sT);
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
      }

    // Mej 
    double M11, M12, M21, M22 = -999;
    double diff_11_22 = 9999;
    double diff_12_21 = 9999;
    if ( (TwoEle) && (TwoJets) ) 
      {
	TLorentzVector jet1, jet2, ele1, ele2;
	ele1.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),0);
	ele2.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[1]),
			  ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[1]),
			  ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[1]),0);
	jet1.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlapEle[0]),
			  CaloJetEta->at(v_idx_jet_PtCut_noOverlapEle[0]),
			  CaloJetPhi->at(v_idx_jet_PtCut_noOverlapEle[0]),0);
	jet2.SetPtEtaPhiM(CaloJetPt->at(v_idx_jet_PtCut_noOverlapEle[1]),
			  CaloJetEta->at(v_idx_jet_PtCut_noOverlapEle[1]),
			  CaloJetPhi->at(v_idx_jet_PtCut_noOverlapEle[1]),0);
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
    
     
    ////////////////////// User's code to be done for every event - END ///////////////////////
    
  } // End of loop over events
  
  //////////write histos 
  
  
  //h_nEleFinal->Write();
  
  //barrel - all
  h_ElectronPt_barrel_all->Write();
  //barrel - heep
  h_ElectronPt_barrel_heep->Write();
  //INFO
  //    //pT of both electrons, to be built using the histograms produced automatically by baseClass
  //    TH1F * h_pTElectrons = new TH1F ("h_pTElectrons","", getHistoNBins("pT1stEle"), getHistoMin("pT1stEle"), getHistoMax("pT1stEle"));
  //    h_pTElectrons->Add( & getHisto_noCuts_or_skim("pT1stEle") ); // all histos can be retrieved, see other getHisto_xxxx methods in baseClass.h
  //    h_pTElectrons->Add( & getHisto_noCuts_or_skim("pT2ndEle") );
  //    //one could also do:  *h_pTElectrons = getHisto_noCuts_or_skim("pT1stEle") + getHisto_noCuts_or_skim("pT2ndEle");
  //    h_pTElectrons->Write();
  //    //one could also do:   const TH1F& h = getHisto_noCuts_or_skim// and use h
  
  //STDOUT("analysisClass::Loop() ends");   
}
