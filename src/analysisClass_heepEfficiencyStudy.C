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

  TH1F *h_N_ele_Gen = new TH1F ("N_ele_Gen","N_ele_Gen",5,-0.5,4.5);   h_N_ele_Gen->Sumw2();
  TH1F *h_N_ele_Gen_etaPtCut = new TH1F ("N_ele_Gen_etaPtCut","N_ele_Gen_etaPtCut",11,-0.5,10.5);   h_N_ele_Gen_etaPtCut->Sumw2();
  TH1F *h_N_ele_Gen_matched = new TH1F ("N_ele_Gen_matched","N_ele_Gen_matched",5,-0.5,4.5);  h_N_ele_Gen_matched->Sumw2();
  //TH1F *h_N_ele_Pt_Gen_matched_ID = new TH1F ("N_ele_Pt_Gen_matched_ID","N_ele_Pt_Gen_matched_ID",11,-0.5,10.5);   h_N_ele_Pt_Gen_matched_ID->Sumw2();
  TH1F *h_N_ele_Pt_Gen_matched_ID_ISO = new TH1F ("N_ele_Pt_Gen_matched_ID_ISO","N_ele_Pt_Gen_matched_ID_ISO",11,-0.5,10.5);   h_N_ele_Pt_Gen_matched_ID_ISO->Sumw2();

  TH1F *h_ele_Pt_Gen  = new TH1F ("ele_Pt_Gen","ele_Pt_Gen",400,0,4000);   h_ele_Pt_Gen->Sumw2();
  TH1F *h_ele_Pt_Gen_etaPtCut  = new TH1F ("ele_Pt_Gen_etaPtCut","ele_Pt_Gen_etaPtCut",400,0,4000);   h_ele_Pt_Gen_etaPtCut->Sumw2();
  TH1F *h_ele_Pt_Gen_matched  = new TH1F ("ele_Pt_Gen_matched","ele_Pt_Gen_matched",400,0,4000);   h_ele_Pt_Gen_matched->Sumw2();
  //TH1F *h_ele_Pt_Gen_matched_ID  = new TH1F ("ele_Pt_Gen_matched_ID","ele_Pt_Gen_matched_ID",400,0,4000);   h_ele_Pt_Gen_matched_ID->Sumw2();
  TH1F *h_ele_Pt_Gen_matched_ID_ISO = new TH1F ("ele_Pt_Gen_matched_ID_ISO","ele_Pt_Gen_matched_ID_ISO",400,0,4000);   h_ele_Pt_Gen_matched_ID_ISO->Sumw2();

  TH1F *h_ele_Eta_Gen  = new TH1F ("ele_Eta_Gen","ele_Eta_Gen",201,-10.05,10.05);   h_ele_Eta_Gen->Sumw2();
  TH1F *h_ele_Eta_Gen_etaPtCut  = new TH1F ("ele_Eta_Gen_etaPtCut","ele_Eta_Gen_etaPtCut",201,-10.05,10.05);   h_ele_Eta_Gen_etaPtCut->Sumw2();
  TH1F *h_ele_Eta_Gen_matched  = new TH1F ("ele_Eta_Gen_matched","ele_Eta_Gen_matched",201,-10.05,10.05);   h_ele_Eta_Gen_matched->Sumw2();
  //TH1F *h_ele_Eta_Gen_matched_ID  = new TH1F ("ele_Eta_Gen_matched_ID","ele_Eta_Gen_matched_ID",201,-10.05,10.05);   h_ele_Eta_Gen_matched_ID->Sumw2();
  TH1F *h_ele_Eta_Gen_matched_ID_ISO  = new TH1F ("ele_Eta_Gen_matched_ID_ISO","ele_Eta_Gen_matched_ID_ISO",201,-10.05,10.05);   h_ele_Eta_Gen_matched_ID_ISO->Sumw2();

  TH2F *h_ele_Eta_Pt_Gen = new TH2F ("ele_Eta_Pt_Gen","ele_Eta_Pt_Gen",400,0,4000,201,-10.05,10.05); h_ele_Eta_Pt_Gen->Sumw2();

  TH1F *h_DeltaR_Gen_Reco = new TH1F("DeltaR_Gen_Reco","DeltaR_Gen_Reco",500,0,10.);   h_DeltaR_Gen_Reco->Sumw2();
  TH1F *h_DeltaR_Gen_2ndReco = new TH1F("DeltaR_Gen_2ndReco","DeltaR_Gen_2ndReco",500,0,10.);   h_DeltaR_Gen_2ndReco->Sumw2();

  TH1F *h_ele_E_matched = new TH1F ("ele_E_matched","ele_E_matched",400,0,4000);   h_ele_E_matched->Sumw2();
  TH1F *h_ele_Pt_matched  = new TH1F ("ele_Pt_matched","ele_Pt_matched",400,0,4000);   h_ele_Pt_matched->Sumw2();
  TH1F *h_ele_Phi_matched  = new TH1F ("ele_Phi_matched","ele_Phi_matched",71,-3.55,3.55);   h_ele_Phi_matched->Sumw2();
  TH1F *h_ele_Eta_matched  = new TH1F ("ele_Eta_matched","ele_Eta_matched",201,-10.05,10.05);   h_ele_Eta_matched->Sumw2();
  //TH1F *h_ele_CaloEnergy_matched  = new TH1F ("ele_CaloEnergy_matched","ele_CaloEnergy_matched",400,0,4000);   h_ele_CaloEnergy_matched->Sumw2();
  TH1F *h_Energy_Res = new TH1F ("Energy_Res","Energy_Res",150,0,1.5);   h_Energy_Res->Sumw2();

  TH1F *h_N_recoEle = new TH1F ("N_recoEle","N_recoEle",5,-0.5,4.5);  h_N_recoEle->Sumw2();
  TH1F *h_N_recoEle_Heep = new TH1F ("N_recoEle_Heep","N_recoEle_Heep",5,-0.5,4.5);  h_N_recoEle_Heep->Sumw2();



  
  ////////////////////// User's code to book histos - END ///////////////////////
    

  ////////////////////// User's code get precut values - BEGIN ///////////////////////

  int electron_PID=int(getPreCutValue1("electronPID"));
  int MotherPID=int(getPreCutValue1("motherPID"));
  float ConeSizeMCmatch_cut=getPreCutValue1("coneSizeMCmatchCut");

  float elePt_cut=getPreCutValue1("elePtCut");
  float eleEta_cut1=getPreCutValue1("eleEtaCut");
  float eleEta_cut2=getPreCutValue2("eleEtaCut");
  float eleEta_cut3=getPreCutValue3("eleEtaCut");
  float genPartPt_cut=getPreCutValue1("genPartPtCut");

  int Ntot_recoEle_Heep = 0;
  int Ntot_recoEle_Heep_barrel = 0;
  int Ntot_recoEle_Heep_endcap = 0;

  ////////////////////// User's code get precut values - END ///////////////////////


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

    int nEleMatched = 0;
    int nGenEle = 0;
    int N_ele_Gen_etaPtCut = 0;
    int nEle_post_ID = 0;
    int nEle_post_ID_ISO = 0;

    vector<int> v_idx_reco_matched;

    for(int igen=0;igen<GenParticleEnergy->size();igen++) // Loop over gen particles
      {

	// select gen particles from their mother
	if( abs(GenParticlePdgId->at(igen))==electron_PID && abs(GenParticlePdgId->at(GenParticleMotherIndex->at(igen)))==MotherPID )
	  {
	    nGenEle++;
	    h_ele_Pt_Gen->Fill(GenParticlePt->at(igen));
	    h_ele_Eta_Gen->Fill(GenParticleEta->at(igen));
	    h_ele_Eta_Pt_Gen->Fill(GenParticlePt->at(igen),GenParticleEta->at(igen));


	    //skip events with gen particle with low pT or outside barrel or endcap acceptance
	    bool genEleInBarrel = false;
	    bool genEleInEndcap = false;
	    bool genEleInEcal = false;
	    if ( fabs(GenParticleEta->at(igen))<eleEta_cut1 )
	      genEleInBarrel = true;
	    if ( fabs(GenParticleEta->at(igen))>eleEta_cut2 && fabs(GenParticleEta->at(igen))<eleEta_cut3)
	      genEleInEndcap = true;
	    if ( genEleInBarrel || genEleInEndcap)
	      genEleInEcal = true;
	    if ( GenParticlePt->at(igen)<genPartPt_cut || !genEleInEcal )
	      continue;

	    N_ele_Gen_etaPtCut++;
	    h_ele_Pt_Gen_etaPtCut->Fill(GenParticlePt->at(igen));
	    h_ele_Eta_Gen_etaPtCut->Fill(GenParticleEta->at(igen));

	    TVector3 elegen;
	    elegen.SetPtEtaPhi(GenParticlePt->at(igen),
			       GenParticleEta->at(igen),
			       GenParticlePhi->at(igen));

	    ///////Calculating DeltaR and matching
	    float minDeltaR = 999;
	    float min2ndDeltaR = 999;
	    int idx_minDeltaR = -1;

	    for (int iele=0; iele<ElectronPt->size(); iele++) { // loop over reco particles for each gen particle
	      if ( ElectronPt->at(iele)<elePt_cut ) continue; 
	            
	      TVector3 ele;
	      ele.SetPtEtaPhi(ElectronPt->at(iele),
			      ElectronEta->at(iele),
			      ElectronPhi->at(iele));
	            
	      float DeltaR_ele_elegen = elegen.DeltaR(ele); // DeltaR between ele and elegen: TMath::Sqrt( deta*deta+dphi*dphi )
	            
	      if (DeltaR_ele_elegen < minDeltaR) {
		min2ndDeltaR=minDeltaR;
		minDeltaR=DeltaR_ele_elegen;
		idx_minDeltaR = iele;
	      }
	      else if (DeltaR_ele_elegen < min2ndDeltaR){
		min2ndDeltaR=DeltaR_ele_elegen;
	      }

	    } // end loop over reco particles for each gen particle

	    h_DeltaR_Gen_Reco->Fill(minDeltaR);
	    h_DeltaR_Gen_2ndReco->Fill(min2ndDeltaR);

	    //gen particle matched with reco candidate
	    if ( minDeltaR < ConeSizeMCmatch_cut )
	      {
		nEleMatched++;
		// the index of the matched reco particle is idx_minDeltaR;
		v_idx_reco_matched.push_back(idx_minDeltaR); //to be used later

		h_ele_Pt_Gen_matched->Fill(GenParticlePt->at(igen));
		h_ele_Eta_Gen_matched->Fill(GenParticleEta->at(igen));

		h_ele_E_matched->Fill(ElectronEnergy->at(idx_minDeltaR));
		h_ele_Pt_matched->Fill(ElectronPt->at(idx_minDeltaR));
		h_ele_Phi_matched->Fill(ElectronPhi->at(idx_minDeltaR));
		h_ele_Eta_matched->Fill(ElectronEta->at(idx_minDeltaR));
		//h_ele_CaloEnergy_matched->Fill(ElectronCaloEnergy->at(idx_minDeltaR));
		h_Energy_Res->Fill(ElectronEnergy->at(idx_minDeltaR)/GenParticleEnergy->at(igen));

		if (ElectronHeepID->at(idx_minDeltaR) == 0) {
		  nEle_post_ID_ISO++;
		  h_ele_Pt_Gen_matched_ID_ISO->Fill(GenParticlePt->at(igen));
		  h_ele_Eta_Gen_matched_ID_ISO->Fill(GenParticleEta->at(igen));
		}

	      } // end gen particle matched with reco candidate

	  } // select gen particles from their mother - end

      } // Loop over gen particles - end 

    h_N_ele_Gen->Fill(nGenEle);
    h_N_ele_Gen_etaPtCut->Fill(N_ele_Gen_etaPtCut);
    h_N_ele_Gen_matched->Fill(nEleMatched);
    //h_N_ele_Pt_Gen_matched_ID->Fill(nEle_post_ID);
    h_N_ele_Pt_Gen_matched_ID_ISO->Fill(nEle_post_ID_ISO);


    int N_recoEle = 0;
    int N_recoEle_Heep = 0;
    for (int iele=0; iele<ElectronPt->size(); iele++) { // loop over reco particles
      bool recoEleInBarrel = false;
      bool recoEleInEndcap = false;
      bool recoEleInEcal = false;
      if ( fabs(ElectronSCEta->at(iele))<eleEta_cut1 )
	recoEleInBarrel = true;
      if ( fabs(ElectronSCEta->at(iele))>eleEta_cut2 && fabs(ElectronSCEta->at(iele))<eleEta_cut3)
	recoEleInEndcap = true;
      if ( recoEleInBarrel || recoEleInEndcap)
	recoEleInEcal = true;
      if ( !recoEleInEcal ) // skip reco ele in gaps
	continue;
      N_recoEle++;
      if (ElectronHeepID->at(iele) == 0) {
	N_recoEle_Heep++;
	Ntot_recoEle_Heep++;
	if (recoEleInBarrel )
	  Ntot_recoEle_Heep_barrel++;
	if (recoEleInEndcap )
	  Ntot_recoEle_Heep_endcap++;
      }
    } // end loop over reco particles
    h_N_recoEle->Fill(N_recoEle);
    h_N_recoEle_Heep->Fill(N_recoEle_Heep);

    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();
    

    // Set the value of the variableNames listed in the cutFile to their current value
    
    // HLT
    //fillVariableWithValue( "HLT", PassTrig ) ;

//     // 1st ele
//     if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 1 ) 
//       {
// 	fillVariableWithValue( "Pt1stEle_IDISO_NoOvrlp", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
//       }



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

  STDOUT("Ntot_recoEle_Heep in ECAL, Barrel, Endcap = "<<Ntot_recoEle_Heep<<", "<< Ntot_recoEle_Heep_barrel<<", "<< Ntot_recoEle_Heep_endcap);  

  ////////////////////// User's code to write histos - BEGIN ///////////////////////

  h_N_ele_Gen->Write();
  h_N_ele_Gen_etaPtCut->Write();
  h_N_ele_Gen_matched->Write();
  //h_N_ele_Pt_Gen_matched_ID->Write();
  h_N_ele_Pt_Gen_matched_ID_ISO->Write();

  h_ele_Pt_Gen->Write();
  h_ele_Pt_Gen_etaPtCut->Write();
  h_ele_Pt_Gen_matched->Write();
  //h_ele_Pt_Gen_matched_ID->Write();
  h_ele_Pt_Gen_matched_ID_ISO->Write();

  h_ele_Eta_Gen->Write();
  h_ele_Eta_Gen_etaPtCut->Write();
  h_ele_Eta_Gen_matched->Write();
  //h_ele_Eta_Gen_matched_ID->Write();
  h_ele_Eta_Gen_matched_ID_ISO->Write();

  h_ele_Eta_Pt_Gen->Write();

  h_DeltaR_Gen_Reco->Write();
  h_DeltaR_Gen_2ndReco->Write();

  h_ele_E_matched->Write();
  h_ele_Pt_matched->Write();
  h_ele_Phi_matched->Write();
  h_ele_Eta_matched->Write();
  //h_ele_CaloEnergy_matched->Write();
  h_Energy_Res->Write();

  h_N_recoEle->Write();
  h_N_recoEle_Heep->Write();

  ////////////////////// User's code to write histos - END ///////////////////////
  
  
  //STDOUT("analysisClass::Loop() ends");   
}
