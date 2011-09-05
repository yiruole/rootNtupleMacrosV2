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

  int PDGID_AG =  getPreCutValue1("PDGID_axigluon");
  int PDGID_Wplus =  getPreCutValue1("PDGID_W");
  int PDGID_Wminus =  getPreCutValue2("PDGID_W");

  int PDGID_ELE =  getPreCutValue1("PDGID_ELE_MU_TAU");
  int PDGID_MU =  getPreCutValue2("PDGID_ELE_MU_TAU");
  int PDGID_TAU =  getPreCutValue3("PDGID_ELE_MU_TAU");
  int PDGID_NUELE =  getPreCutValue1("PDGID_NU_ELE_MU_TAU");
  int PDGID_NUMU =  getPreCutValue2("PDGID_NU_ELE_MU_TAU");
  int PDGID_NUTAU =  getPreCutValue3("PDGID_NU_ELE_MU_TAU");

  int PDGID_D =  getPreCutValue1("PDGID_D");
  int PDGID_DBAR =  getPreCutValue2("PDGID_D");
  int PDGID_U =  getPreCutValue1("PDGID_U");
  int PDGID_UBAR =  getPreCutValue2("PDGID_U");
  int PDGID_S =  getPreCutValue1("PDGID_S");
  int PDGID_SBAR =  getPreCutValue2("PDGID_S");
  int PDGID_C =  getPreCutValue1("PDGID_C");
  int PDGID_CBAR =  getPreCutValue2("PDGID_C");
  int PDGID_B =  getPreCutValue1("PDGID_B");
  int PDGID_BBAR =  getPreCutValue2("PDGID_B");
  int PDGID_T =  getPreCutValue1("PDGID_T");
  int PDGID_TBAR =  getPreCutValue2("PDGID_T");

  ////////////////////// User's code to get preCut values - END /////////////////

  ////////////////////// User's code to book histos - BEGIN ///////////////////////

  //CreateUserTH1D("h1_LQGenEle_Pt", getHistoNBins("Pt1stEle_PAS"), getHistoMin("Pt1stEle_PAS"), getHistoMax("Pt1stEle_PAS") );
  //CreateUserTH2D("h2_pfMET_vs_neutrinoPt__sT", 100,0,1000, 100, 0, 1000);

  CreateUserTH1D("h1_Num_AG", 3, -0.5, 2.5);
  CreateUserTH1D("h1_Num_W", 3, -0.5, 2.5 );
  CreateUserTH1D("h1_Charge_W", 2, -1.01, 1.01);

  CreateUserTH1D("h1_Mass_AG", 1000, 0, 2000);
  CreateUserTH1D("h1_Mass_AG_fromqq", 1000, 0, 2000);
  CreateUserTH1D("h1_Mass_W", 100, 0, 200);
  CreateUserTH1D("h1_Mass_W_fromlnu", 100, 0, 200);
  CreateUserTH1D("h1_MT_W", 100, 0, 200);

  CreateUserTH1D("h1_Pt_AG", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_W", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_Wplus", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_Wminus", 100, 0, 1000);

  CreateUserTH1D("h1_Eta_AG", 100, -10, 10);
  CreateUserTH1D("h1_Eta_W", 100, -10, 10);
  CreateUserTH1D("h1_Eta_Wplus", 100, -10, 10);
  CreateUserTH1D("h1_Eta_Wminus", 100, -10, 10);

  CreateUserTH1D("h1_Pt_AG_q", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_AG_qbar", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_W_l", 100, 0, 1000);
  CreateUserTH1D("h1_Pt_W_nu", 100, 0, 1000);

  CreateUserTH1D("h1_Eta_AG_q", 100, -10, 10);
  CreateUserTH1D("h1_Eta_AG_qbar", 100, -10, 10);
  CreateUserTH1D("h1_Eta_W_l", 100, -10, 10);
  CreateUserTH1D("h1_Eta_W_nu", 100, -10, 10);

  CreateUserTH1D("h1_Num_AG_Daughters", 20, 0, 20);		

  CreateUserTH1D("h1_GenMET", 100, 0, 1000);

  TH1F* h1_DecayProducts_AG = new TH1F("h1_DecayProducts_AG", "h1_DecayProducts_AG", 7, 0.5, 7.5);
  h1_DecayProducts_AG->Sumw2();
  h1_DecayProducts_AG->GetXaxis()->SetBinLabel(1,"ddbar");
  h1_DecayProducts_AG->GetXaxis()->SetBinLabel(2,"uubar");
  h1_DecayProducts_AG->GetXaxis()->SetBinLabel(3,"ssbar");
  h1_DecayProducts_AG->GetXaxis()->SetBinLabel(4,"ccbar");
  h1_DecayProducts_AG->GetXaxis()->SetBinLabel(5,"bbar");
  h1_DecayProducts_AG->GetXaxis()->SetBinLabel(6,"ttbar");
  h1_DecayProducts_AG->GetXaxis()->SetBinLabel(7,"total");

  TH1F* h1_DecayProducts_W = new TH1F("h1_DecayProducts_W", "h1_DecayProducts_W", 4, 0.5, 4.5);
  h1_DecayProducts_W->Sumw2();
  h1_DecayProducts_W->GetXaxis()->SetBinLabel(1,"e#nu");
  h1_DecayProducts_W->GetXaxis()->SetBinLabel(2,"#mu#nu");
  h1_DecayProducts_W->GetXaxis()->SetBinLabel(3,"#tau#nu");
  h1_DecayProducts_W->GetXaxis()->SetBinLabel(4,"total");

 
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

    ////////////////////// User's code to be done for every event - BEGIN ///////////////////////

    //## Define new jet collection
    //std::auto_ptr<std::vector<double> >  JetPt  ( new std::vector<double>()  );
    //std::auto_ptr<std::vector<int> >     JetPassID  ( new std::vector<int>()  );

    //FillUserTH1D("h1_ElectronPt_all", ElectronPt->at(iele) );

    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();

    // Set the value of the variableNames listed in the cutFile to their current value
    fillVariableWithValue( "PassJSON", passJSON(run, ls, isData) );    

    // Evaluate cuts (but do not apply them)
    evaluateCuts();

    int Num_AG = 0;
    int Num_AG_Daughters = 0;	      
    TLorentzVector v_AG;
    TLorentzVector v_AG_boosted;
    TLorentzVector v_AG_q;
    TLorentzVector v_AG_qbar;
    int            flavour_AG_q = 0;
    int            flavour_AG_qbar = 0;
    int            flavour_W_lepton = 0;
    int            flavour_W_neutrino = 0;

    int Num_W = 0; 
    int Wcharge = 0;
    int index_W = -1;
    TLorentzVector v_W;
    TLorentzVector v_W_l;
    TLorentzVector v_W_nu;
    double MT_W = 0;

    int T_EXISTS = 0;
    int TBAR_EXISTS = 0;

    if(isData==0)
      {
	for (int genp=0; genp<GenParticlePdgId->size(); genp++)
	  {

	    //--- printout
	    //  	    	     	    cout << "idx: " << genp
	    //  	    	     		 << " GenParticlePdgId: " << GenParticlePdgId->at(genp)
	    //  	    	     		 << " GenParticleMotherIndex: " << GenParticleMotherIndex->at(genp)
	    //  	    	     		 << " GenParticleStatus: " << GenParticleStatus->at(genp)
	    //  	    	     		 << " GenParticlePt: " << GenParticlePt->at(genp)
	    //  	    	     		 << endl;	    
	    //---

	    TLorentzVector current_p4;
	    current_p4.SetPtEtaPhiE(GenParticlePt->at(genp),
				    GenParticleEta->at(genp),
				    GenParticlePhi->at(genp),
				    GenParticleEnergy->at(genp)
				    );
	    
	    //AG
	    if(GenParticlePdgId->at(genp) == PDGID_AG && GenParticleStatus->at(genp) == 3)
	      {
		Num_AG++;
		Num_AG_Daughters = GenParticleNumDaught->at(genp);
		v_AG = current_p4;
	      }//

	    //qq from AG
	    if( GenParticleMotherIndex->at(genp) != -1 )
	      {  
		if( abs( GenParticlePdgId->at( GenParticleMotherIndex->at(genp) ) ) == PDGID_AG && GenParticleStatus->at(genp) == 3)
		  {
		    //First TOP = assuming that top comes always before than bottom in the gen particle list (should be always true)
		    if( GenParticlePdgId->at(genp) == PDGID_T )
		      {//t
			T_EXISTS = 1;
			v_AG_q = current_p4;
			flavour_AG_q = GenParticlePdgId->at(genp);			
		      }
		    if( GenParticlePdgId->at(genp) == PDGID_TBAR )
		      {//tbar
			TBAR_EXISTS = 1;
			v_AG_qbar = current_p4;
			flavour_AG_qbar = GenParticlePdgId->at(genp);
		      }

		    //Then D, U, S, C
		    if( ( GenParticlePdgId->at(genp) == PDGID_D 
			|| GenParticlePdgId->at(genp) == PDGID_U  
			|| GenParticlePdgId->at(genp) == PDGID_S  
			  || GenParticlePdgId->at(genp) == PDGID_C ) 
			&& T_EXISTS == 0
			) 
		      {//quark
			v_AG_q = current_p4;
			flavour_AG_q = GenParticlePdgId->at(genp);
		      }		    
		    if( ( GenParticlePdgId->at(genp) == PDGID_DBAR 
			|| GenParticlePdgId->at(genp) == PDGID_UBAR  
			|| GenParticlePdgId->at(genp) == PDGID_SBAR  
			  || GenParticlePdgId->at(genp) == PDGID_CBAR )
			&& TBAR_EXISTS == 0
			) 
		      {//antiquark
			v_AG_qbar = current_p4;
			flavour_AG_qbar = GenParticlePdgId->at(genp);
		      }
		    		   
		    //Finally BOTTOM
		    if( GenParticlePdgId->at(genp) == PDGID_B && T_EXISTS==0 )
		      {//b
			v_AG_q = current_p4;
			flavour_AG_q = GenParticlePdgId->at(genp);						
		      }
		    if( GenParticlePdgId->at(genp) == PDGID_BBAR && TBAR_EXISTS==0 )
		      {//bbar
			v_AG_qbar = current_p4;
			flavour_AG_qbar = GenParticlePdgId->at(genp);
		      }		    

		  }
	      }//

	    //W
	    if( abs(GenParticlePdgId->at(genp)) == abs(PDGID_Wplus) 
		&& GenParticleStatus->at(genp) == 3 
		&& abs( GenParticlePdgId->at( GenParticleMotherIndex->at(genp) ) ) != PDGID_AG 
		&& abs( GenParticlePdgId->at( GenParticleMotherIndex->at(genp) ) ) != abs(PDGID_T) 
		)
	      {
		Num_W++;
		v_W = current_p4;
		index_W = genp;

		//W+
		if(GenParticlePdgId->at(genp) == PDGID_Wplus)
		  Wcharge = 1;

		//W-
		if(GenParticlePdgId->at(genp) == PDGID_Wminus)
		  Wcharge = -1;
	      }

	    //leptons from W
	    if( GenParticleMotherIndex->at(genp) != -1 )
	      {  
		if( GenParticleMotherIndex->at(genp) == index_W )
		  {
		    if( abs(GenParticlePdgId->at(genp)) == abs(PDGID_ELE)
			|| abs(GenParticlePdgId->at(genp)) == abs(PDGID_MU)  
			|| abs(GenParticlePdgId->at(genp)) == abs(PDGID_TAU)  
			) 
		      {//lepton
			v_W_l = current_p4;
			flavour_W_lepton = GenParticlePdgId->at(genp);
		      }

		    if( abs(GenParticlePdgId->at(genp)) == abs(PDGID_NUELE)
			|| abs(GenParticlePdgId->at(genp)) == abs(PDGID_NUMU)  
			|| abs(GenParticlePdgId->at(genp)) == abs(PDGID_NUTAU)  
			) 
		      {//neutrino
			v_W_nu = current_p4;
			flavour_W_neutrino = GenParticlePdgId->at(genp);
		      }
		  }
	      }//

	    //MT
	    MT_W = sqrt( 2 * v_W_l.Pt() * v_W_nu.Pt() * (1 - cos(v_W_l.DeltaPhi(v_W_nu) ) ) );	    
	    
	  }//end loop over gen particles


// 	v_AG_boosted = v_AG;
// 	TVector3 v_beta_AG;
// 	v_beta_AG = v_AG.Vect();
// 	v_beta_AG.SetMag( v_AG.Beta() );		
// 	v_AG_boosted.Boost(v_beta_AG);

// 	cout << "v_AG.E(): " << v_AG.E() << endl;
// 	cout << "v_AG_boosted.E(): " << v_AG_boosted.E() << endl;


      }

    // Fill histograms and do analysis based on cut evaluation

    if(flavour_AG_q != -(flavour_AG_qbar))
      {
	cout << "---------" << endl;
	for (int genp1=0; genp1<GenParticlePdgId->size(); genp1++)
	  {			
	    cout << "idx: " << genp1
		 << " GenParticlePdgId: " << GenParticlePdgId->at(genp1)
		 << " GenParticleMotherIndex: " << GenParticleMotherIndex->at(genp1)
		 << " GenParticleStatus: " << GenParticleStatus->at(genp1)
		 << " GenParticlePt: " << GenParticlePt->at(genp1)
		 << endl;	    
	  }	
      }

    FillUserTH1D("h1_Num_AG", Num_AG);
    FillUserTH1D("h1_Num_W", Num_W);
    FillUserTH1D("h1_Charge_W", Wcharge);

    FillUserTH1D("h1_Mass_AG", v_AG.M() );
    FillUserTH1D("h1_Mass_AG_fromqq", (v_AG_q+v_AG_qbar).M() );
    FillUserTH1D("h1_Mass_W", v_W.M() );
    FillUserTH1D("h1_Mass_W_fromlnu", (v_W_l+v_W_nu).M() );
    FillUserTH1D("h1_MT_W", MT_W );

    FillUserTH1D("h1_Pt_AG", v_AG.Pt() );
    FillUserTH1D("h1_Pt_W", v_W.Pt() );
    if(	  Wcharge == 1  )
      FillUserTH1D("h1_Pt_Wplus", v_W.Pt() );
    if(	  Wcharge == -1  )
      FillUserTH1D("h1_Pt_Wminus", v_W.Pt() );
    
    FillUserTH1D("h1_Eta_AG", v_AG.Eta() );
    FillUserTH1D("h1_Eta_W", v_W.Eta() );
    if(	  Wcharge == 1  )
      FillUserTH1D("h1_Eta_Wplus", v_W.Eta() );
    if(	  Wcharge == -1  )
      FillUserTH1D("h1_Eta_Wminus", v_W.Eta() );

    FillUserTH1D("h1_Pt_AG_q", v_AG_q.Pt() );
    FillUserTH1D("h1_Pt_AG_qbar", v_AG_qbar.Pt() );
    FillUserTH1D("h1_Pt_W_l", v_W_l.Pt() );
    FillUserTH1D("h1_Pt_W_nu", v_W_nu.Pt() );
    
    FillUserTH1D("h1_Eta_AG_q", v_AG_q.Eta() );
    FillUserTH1D("h1_Eta_AG_qbar", v_AG_qbar.Eta() );
    FillUserTH1D("h1_Eta_W_l", v_W_l.Eta() );
    FillUserTH1D("h1_Eta_W_nu", v_W_nu.Eta() );

    if( flavour_AG_q == -(flavour_AG_qbar) )
      {
	h1_DecayProducts_AG->AddBinContent( abs(flavour_AG_q) );
	h1_DecayProducts_AG->AddBinContent( 7 );
      }
    else
      {
	cout << flavour_AG_q << endl;
	cout << flavour_AG_qbar << endl;
	cout << "stonzo!!!!!!!!!!!!!!!!!!!!" << endl;
      }

    if( fabs(flavour_W_lepton) == abs(PDGID_ELE) && fabs(flavour_W_neutrino) == abs(PDGID_NUELE) )
      h1_DecayProducts_W->AddBinContent( 1 );
    if( fabs(flavour_W_lepton) == abs(PDGID_MU) && fabs(flavour_W_neutrino) == abs(PDGID_NUMU) )
      h1_DecayProducts_W->AddBinContent( 2 );
    if( fabs(flavour_W_lepton) == abs(PDGID_TAU) && fabs(flavour_W_neutrino) == abs(PDGID_NUTAU) )
      h1_DecayProducts_W->AddBinContent( 3 );
    if( (fabs(flavour_W_lepton) == abs(PDGID_ELE) && fabs(flavour_W_neutrino) == abs(PDGID_NUELE)) 
	|| (fabs(flavour_W_lepton) == abs(PDGID_MU) && fabs(flavour_W_neutrino) == abs(PDGID_NUMU)) 
	|| (fabs(flavour_W_lepton) == abs(PDGID_TAU) && fabs(flavour_W_neutrino) == abs(PDGID_NUTAU)) )
      h1_DecayProducts_W->AddBinContent( 4 );

    FillUserTH1D("h1_Num_AG_Daughters", Num_AG_Daughters);		

    FillUserTH1D("h1_GenMET", GenMETTrue->at(0) );
    
    //INFO
    //      // retrieve value of previously filled variables (after making sure that they were filled)
    //      double totpTEle;
    //      if ( variableIsFilled("pT1stEle") && variableIsFilled("pT2ndEle") )
    //        totpTEle = getVariableValue("pT1stEle")+getVariableValue("pT2ndEle");
    //      // reject events that did not pass level 0 cuts
    //      if( !passedCut("0") ) continue;
    //      // ......

    //after pre-selection
    //    if( passedAllPreviousCuts("Pt1stEle_PAS")
    //	&& variableIsFilled("MTenu_PAS") && variableIsFilled("sT_PAS")
    //FillUserTH1D("h1_MTenu_PAS_plus", getVariableValue("MTenu_PAS"));

    // Produce skim
    //if( passedAllPreviousCuts("minDRej") ) fillSkimTree();

    ////////////////////// User's code to be done for every event - END ///////////////////////

  } // End of loop over events


  ////////////////////// User's code to write histos - BEGIN ///////////////////////

  h1_DecayProducts_AG->Write();
  h1_DecayProducts_W->Write();

  ////////////////////// User's code to write histos - END ///////////////////////

  delete h1_DecayProducts_AG;
  
  //STDOUT("analysisClass::Loop() ends");
}
