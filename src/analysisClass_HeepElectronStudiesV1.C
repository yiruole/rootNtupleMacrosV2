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
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;

  std::cout << "analysisClass::analysisClass(): ends " << std::endl;
}

analysisClass::~analysisClass()
{
  std::cout << "analysisClass::~analysisClass(): begins " << std::endl;

  std::cout << "analysisClass::~analysisClass(): ends " << std::endl;
}

void analysisClass::Loop()
{
   std::cout << "analysisClass::Loop() begins" <<std::endl;   
    
   if (fChain == 0) return;
   
   //////////book histos here

   // number of electrons
   //    TH1F *h_nEleFinal = new TH1F ("h_nEleFinal","",11,-0.5,10.5);
   //    h_nEleFinal->Sumw2();


   //barrel - all   
   TH1F *h_ElectronPt_barrel_all = new TH1F ("h_ElectronPt_barrel_all","h_ElectronPt_barrel_all",100,0,500);
   TH1F *h_ElectronSCEta_fabs_barrel_all = new TH1F ("h_ElectronSCEta_fabs_barrel_all","h_ElectronSCEta_fabs_barrel_all",50,0,5);
   TH1F *h_ElectronDeltaEtaTrkSC_barrel_all = new TH1F ("h_ElectronDeltaEtaTrkSC_barrel_all","h_ElectronDeltaEtaTrkSC_barrel_all",200,-0.05,0.05);
   TH1F *h_ElectronDeltaPhiTrkSC_barrel_all = new TH1F ("h_ElectronDeltaPhiTrkSC_barrel_all","h_ElectronDeltaPhiTrkSC_barrel_all",200,-0.5,0.5);
   TH1F *h_ElectronHoE_barrel_all = new TH1F ("h_ElectronHoE_barrel_all","h_ElectronHoE_barrel_all",75,0,0.15);
   TH1F *h_ElectronSigmaIEtaIEta_barrel_all = new TH1F ("h_ElectronSigmaIEtaIEta_barrel_all","h_ElectronSigmaIEtaIEta_barrel_all",100,0,0.1);
   TH2F *h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_barrel_all = new TH2F ("h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_barrel_all"
									,"h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_barrel_all"
									,100,0,1,100,0,1);
   TH2F *h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_barrel_all 
     = new TH2F ("h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_barrel_all"
		 ,"h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_barrel_all"
		 ,100,0,500,200,0,100);
   TH1F *h_ElectronHcalIsoD2Heep_barrel_all = new TH1F ("h_ElectronHcalIsoD2Heep_barrel_all","h_ElectronHcalIsoD2Heep_barrel_all",200,0,100);
   TH1F *h_ElectronTrkIsoHeep_barrel_all = new TH1F ("h_ElectronTrkIsoHeep_barrel_all","h_ElectronTrkIsoHeep_barrel_all",200,0,100);


   h_ElectronPt_barrel_all->Sumw2();
   h_ElectronSCEta_fabs_barrel_all->Sumw2();
   h_ElectronDeltaEtaTrkSC_barrel_all->Sumw2();
   h_ElectronDeltaPhiTrkSC_barrel_all->Sumw2();
   h_ElectronHoE_barrel_all->Sumw2();
   h_ElectronSigmaIEtaIEta_barrel_all->Sumw2();
   h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_barrel_all->Sumw2();
   h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_barrel_all->Sumw2(); 
   h_ElectronHcalIsoD2Heep_barrel_all->Sumw2();
   h_ElectronTrkIsoHeep_barrel_all->Sumw2();


   //barrel - heep   
   TH1F *h_ElectronPt_barrel_heep = new TH1F ("h_ElectronPt_barrel_heep","h_ElectronPt_barrel_heep",100,0,500);
   TH1F *h_ElectronSCEta_fabs_barrel_heep = new TH1F ("h_ElectronSCEta_fabs_barrel_heep","h_ElectronSCEta_fabs_barrel_heep",50,0,5);
   TH1F *h_ElectronDeltaEtaTrkSC_barrel_heep = new TH1F ("h_ElectronDeltaEtaTrkSC_barrel_heep","h_ElectronDeltaEtaTrkSC_barrel_heep",200,-0.05,0.05);
   TH1F *h_ElectronDeltaPhiTrkSC_barrel_heep = new TH1F ("h_ElectronDeltaPhiTrkSC_barrel_heep","h_ElectronDeltaPhiTrkSC_barrel_heep",200,-0.5,0.5);
   TH1F *h_ElectronHoE_barrel_heep = new TH1F ("h_ElectronHoE_barrel_heep","h_ElectronHoE_barrel_heep",75,0,0.15);
   TH1F *h_ElectronSigmaIEtaIEta_barrel_heep = new TH1F ("h_ElectronSigmaIEtaIEta_barrel_heep","h_ElectronSigmaIEtaIEta_barrel_heep",100,0,0.1);
   TH2F *h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_barrel_heep = new TH2F ("h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_barrel_heep"
									,"h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_barrel_heep"
									,100,0,1,100,0,1);
   TH2F *h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_barrel_heep 
     = new TH2F ("h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_barrel_heep"
		 ,"h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_barrel_heep"
		 ,100,0,500,200,0,100);
   TH1F *h_ElectronHcalIsoD2Heep_barrel_heep = new TH1F ("h_ElectronHcalIsoD2Heep_barrel_heep","h_ElectronHcalIsoD2Heep_barrel_heep",200,0,100);
   TH1F *h_ElectronTrkIsoHeep_barrel_heep = new TH1F ("h_ElectronTrkIsoHeep_barrel_heep","h_ElectronTrkIsoHeep_barrel_heep",200,0,100);

   h_ElectronPt_barrel_heep->Sumw2();
   h_ElectronSCEta_fabs_barrel_heep->Sumw2();
   h_ElectronDeltaEtaTrkSC_barrel_heep->Sumw2();
   h_ElectronDeltaPhiTrkSC_barrel_heep->Sumw2();
   h_ElectronHoE_barrel_heep->Sumw2();
   h_ElectronSigmaIEtaIEta_barrel_heep->Sumw2();
   h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_barrel_heep->Sumw2();
   h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_barrel_heep->Sumw2(); 
   h_ElectronHcalIsoD2Heep_barrel_heep->Sumw2();
   h_ElectronTrkIsoHeep_barrel_heep->Sumw2();


   //endcap - all   
   TH1F *h_ElectronPt_endcap_all = new TH1F ("h_ElectronPt_endcap_all","h_ElectronPt_endcap_all",100,0,500);
   TH1F *h_ElectronSCEta_fabs_endcap_all = new TH1F ("h_ElectronSCEta_fabs_endcap_all","h_ElectronSCEta_fabs_endcap_all",50,0,5);
   TH1F *h_ElectronDeltaEtaTrkSC_endcap_all = new TH1F ("h_ElectronDeltaEtaTrkSC_endcap_all","h_ElectronDeltaEtaTrkSC_endcap_all",200,-0.05,0.05);
   TH1F *h_ElectronDeltaPhiTrkSC_endcap_all = new TH1F ("h_ElectronDeltaPhiTrkSC_endcap_all","h_ElectronDeltaPhiTrkSC_endcap_all",200,-0.5,0.5);
   TH1F *h_ElectronHoE_endcap_all = new TH1F ("h_ElectronHoE_endcap_all","h_ElectronHoE_endcap_all",75,0,0.15);
   TH1F *h_ElectronSigmaIEtaIEta_endcap_all = new TH1F ("h_ElectronSigmaIEtaIEta_endcap_all","h_ElectronSigmaIEtaIEta_endcap_all",100,0,0.1);
   TH2F *h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_endcap_all = new TH2F ("h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_endcap_all"
									,"h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_endcap_all"
									,100,0,1,100,0,1);
   TH2F *h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_endcap_all 
     = new TH2F ("h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_endcap_all"
		 ,"h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_endcap_all"
		 ,100,0,500,200,0,100);
   TH1F *h_ElectronHcalIsoD2Heep_endcap_all = new TH1F ("h_ElectronHcalIsoD2Heep_endcap_all","h_ElectronHcalIsoD2Heep_endcap_all",200,0,100);
   TH1F *h_ElectronTrkIsoHeep_endcap_all = new TH1F ("h_ElectronTrkIsoHeep_endcap_all","h_ElectronTrkIsoHeep_endcap_all",200,0,100);

   h_ElectronPt_endcap_all->Sumw2();
   h_ElectronSCEta_fabs_endcap_all->Sumw2();
   h_ElectronDeltaEtaTrkSC_endcap_all->Sumw2();
   h_ElectronDeltaPhiTrkSC_endcap_all->Sumw2();
   h_ElectronHoE_endcap_all->Sumw2();
   h_ElectronSigmaIEtaIEta_endcap_all->Sumw2();
   h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_endcap_all->Sumw2();
   h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_endcap_all->Sumw2(); 
   h_ElectronHcalIsoD2Heep_endcap_all->Sumw2();
   h_ElectronTrkIsoHeep_endcap_all->Sumw2();


   //endcap - heep   
   TH1F *h_ElectronPt_endcap_heep = new TH1F ("h_ElectronPt_endcap_heep","h_ElectronPt_endcap_heep",100,0,500);
   TH1F *h_ElectronSCEta_fabs_endcap_heep = new TH1F ("h_ElectronSCEta_fabs_endcap_heep","h_ElectronSCEta_fabs_endcap_heep",50,0,5);
   TH1F *h_ElectronDeltaEtaTrkSC_endcap_heep = new TH1F ("h_ElectronDeltaEtaTrkSC_endcap_heep","h_ElectronDeltaEtaTrkSC_endcap_heep",200,-0.05,0.05);
   TH1F *h_ElectronDeltaPhiTrkSC_endcap_heep = new TH1F ("h_ElectronDeltaPhiTrkSC_endcap_heep","h_ElectronDeltaPhiTrkSC_endcap_heep",200,-0.5,0.5);
   TH1F *h_ElectronHoE_endcap_heep = new TH1F ("h_ElectronHoE_endcap_heep","h_ElectronHoE_endcap_heep",75,0,0.15);
   TH1F *h_ElectronSigmaIEtaIEta_endcap_heep = new TH1F ("h_ElectronSigmaIEtaIEta_endcap_heep","h_ElectronSigmaIEtaIEta_endcap_heep",100,0,0.1);
   TH2F *h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_endcap_heep = new TH2F ("h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_endcap_heep"
									,"h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_endcap_heep"
									,100,0,1,100,0,1);
   TH2F *h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_endcap_heep 
     = new TH2F ("h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_endcap_heep"
		 ,"h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_endcap_heep"
		 ,100,0,500,200,0,100);
   TH1F *h_ElectronHcalIsoD2Heep_endcap_heep = new TH1F ("h_ElectronHcalIsoD2Heep_endcap_heep","h_ElectronHcalIsoD2Heep_endcap_heep",200,0,100);
   TH1F *h_ElectronTrkIsoHeep_endcap_heep = new TH1F ("h_ElectronTrkIsoHeep_endcap_heep","h_ElectronTrkIsoHeep_endcap_heep",200,0,100);

   h_ElectronPt_endcap_heep->Sumw2();
   h_ElectronSCEta_fabs_endcap_heep->Sumw2();
   h_ElectronDeltaEtaTrkSC_endcap_heep->Sumw2();
   h_ElectronDeltaPhiTrkSC_endcap_heep->Sumw2();
   h_ElectronHoE_endcap_heep->Sumw2();
   h_ElectronSigmaIEtaIEta_endcap_heep->Sumw2();
   h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_endcap_heep->Sumw2();
   h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_endcap_heep->Sumw2(); 
   h_ElectronHcalIsoD2Heep_endcap_heep->Sumw2();
   h_ElectronTrkIsoHeep_endcap_heep->Sumw2();


   /////////initialize variables

   Long64_t nentries = fChain->GetEntriesFast();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
   ////// If the root version is updated and rootNtupleClass regenerated,     /////
   ////// these lines may need to be updated.                                 /////    
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << std::endl;   
     // if (Cut(ientry) < 0) continue;

     ////////////////////// User's code starts here ///////////////////////

     ///Stuff to be done every event
     
     //Loop over electrons
     for(int ele=0; ele<ElectronPt->size(); ele++)
       {

	 int isBarrel = 0;
	 int isEndcap = 0;

	 if( fabs( ElectronSCEta->at(ele) ) < 1.442 ) 
	   isBarrel = 1;

	 if( fabs( ElectronSCEta->at(ele) ) > 1.560 && fabs( ElectronSCEta->at(ele) ) < 2.5 ) 
	   isEndcap = 1;
	     
	 if(isBarrel)
	   {
	     h_ElectronPt_barrel_all->Fill( ElectronPt->at(ele) );
	     h_ElectronSCEta_fabs_barrel_all->Fill( fabs( ElectronSCEta->at(ele) ) );
	     h_ElectronDeltaEtaTrkSC_barrel_all->Fill( ElectronDeltaEtaTrkSC->at(ele) ) ;
	     h_ElectronDeltaPhiTrkSC_barrel_all->Fill( ElectronDeltaPhiTrkSC->at(ele) ) ;
	     h_ElectronHoE_barrel_all->Fill( ElectronHoE->at(ele) );
	     h_ElectronSigmaIEtaIEta_barrel_all->Fill( ElectronSigmaIEtaIEta->at(ele) );
	     h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_barrel_all->Fill( ElectronE1x5OverE5x5->at(ele) , ElectronE2x5OverE5x5->at(ele) );
	     h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_barrel_all
	       ->Fill( ElectronPt->at(ele) , ElectronEcalIsoHeep->at(ele)+ElectronHcalIsoD1Heep->at(ele) );
	     h_ElectronHcalIsoD2Heep_barrel_all->Fill( ElectronHcalIsoD2Heep->at(ele) );
	     h_ElectronTrkIsoHeep_barrel_all->Fill( ElectronTrkIsoHeep->at(ele) );
	     
	     if( ElectronHeepID->at(ele) == 0 )
	       {
		 h_ElectronPt_barrel_heep->Fill( ElectronPt->at(ele) );
		 h_ElectronSCEta_fabs_barrel_heep->Fill( fabs(ElectronSCEta->at(ele)) );
		 h_ElectronDeltaEtaTrkSC_barrel_heep->Fill( ElectronDeltaEtaTrkSC->at(ele) ) ;
		 h_ElectronDeltaPhiTrkSC_barrel_heep->Fill( ElectronDeltaPhiTrkSC->at(ele) ) ;
		 h_ElectronHoE_barrel_heep->Fill( ElectronHoE->at(ele) );
		 h_ElectronSigmaIEtaIEta_barrel_heep->Fill( ElectronSigmaIEtaIEta->at(ele) );
		 h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_barrel_heep->Fill( ElectronE1x5OverE5x5->at(ele) , ElectronE2x5OverE5x5->at(ele) );
		 h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_barrel_heep
		   ->Fill( ElectronPt->at(ele) , ElectronEcalIsoHeep->at(ele)+ElectronHcalIsoD1Heep->at(ele) );
		 h_ElectronHcalIsoD2Heep_barrel_heep->Fill( ElectronHcalIsoD2Heep->at(ele) );
		 h_ElectronTrkIsoHeep_barrel_heep->Fill( ElectronTrkIsoHeep->at(ele) );	       
	       }	     
	   }


	 if(isEndcap)
	   {
	     h_ElectronPt_endcap_all->Fill( ElectronPt->at(ele) );
	     h_ElectronSCEta_fabs_endcap_all->Fill( fabs(ElectronSCEta->at(ele)) );
	     h_ElectronDeltaEtaTrkSC_endcap_all->Fill( ElectronDeltaEtaTrkSC->at(ele) ) ;
	     h_ElectronDeltaPhiTrkSC_endcap_all->Fill( ElectronDeltaPhiTrkSC->at(ele) ) ;
	     h_ElectronHoE_endcap_all->Fill( ElectronHoE->at(ele) );
	     h_ElectronSigmaIEtaIEta_endcap_all->Fill( ElectronSigmaIEtaIEta->at(ele) );
	     h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_endcap_all->Fill( ElectronE1x5OverE5x5->at(ele) , ElectronE2x5OverE5x5->at(ele) );
	     h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_endcap_all
	       ->Fill( ElectronPt->at(ele) , ElectronEcalIsoHeep->at(ele)+ElectronHcalIsoD1Heep->at(ele) );
	     h_ElectronHcalIsoD2Heep_endcap_all->Fill( ElectronHcalIsoD2Heep->at(ele) );
	     h_ElectronTrkIsoHeep_endcap_all->Fill( ElectronTrkIsoHeep->at(ele) );
	     
	     if( ElectronHeepID->at(ele) == 0 )
	       {
		 h_ElectronPt_endcap_heep->Fill( ElectronPt->at(ele) );
		 h_ElectronSCEta_fabs_endcap_heep->Fill( fabs(ElectronSCEta->at(ele)) );
		 h_ElectronDeltaEtaTrkSC_endcap_heep->Fill( ElectronDeltaEtaTrkSC->at(ele) ) ;
		 h_ElectronDeltaPhiTrkSC_endcap_heep->Fill( ElectronDeltaPhiTrkSC->at(ele) ) ;
		 h_ElectronHoE_endcap_heep->Fill( ElectronHoE->at(ele) );
		 h_ElectronSigmaIEtaIEta_endcap_heep->Fill( ElectronSigmaIEtaIEta->at(ele) );
		 h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_endcap_heep->Fill( ElectronE1x5OverE5x5->at(ele) , ElectronE2x5OverE5x5->at(ele) );
		 h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_endcap_heep
		   ->Fill( ElectronPt->at(ele) , ElectronEcalIsoHeep->at(ele)+ElectronHcalIsoD1Heep->at(ele) );
		 h_ElectronHcalIsoD2Heep_endcap_heep->Fill( ElectronHcalIsoD2Heep->at(ele) );
		 h_ElectronTrkIsoHeep_endcap_heep->Fill( ElectronTrkIsoHeep->at(ele) );	       
	       }	     
	   }

	 
       }

     // Set the evaluation of the cuts to false and clear the variable values and filled status
     resetCuts();
     
     // Set the value of the variableNames listed in the cutFile to their current value
     //fillVariableWithValue("nEleFinal", ElectronPt->size() ) ;

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

     
     ////////////////////// User's code ends here ///////////////////////

   } // End loop over events

   //////////write histos 


   //h_nEleFinal->Write();

   //barrel - all
   h_ElectronPt_barrel_all->Write();
   h_ElectronSCEta_fabs_barrel_all->Write();
   h_ElectronDeltaEtaTrkSC_barrel_all->Write();
   h_ElectronDeltaPhiTrkSC_barrel_all->Write();
   h_ElectronHoE_barrel_all->Write();
   h_ElectronSigmaIEtaIEta_barrel_all->Write();
   h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_barrel_all->Write();
   h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_barrel_all->Write(); 
   h_ElectronHcalIsoD2Heep_barrel_all->Write();
   h_ElectronTrkIsoHeep_barrel_all->Write();

   //barrel - heep
   h_ElectronPt_barrel_heep->Write();
   h_ElectronSCEta_fabs_barrel_heep->Write();
   h_ElectronDeltaEtaTrkSC_barrel_heep->Write();
   h_ElectronDeltaPhiTrkSC_barrel_heep->Write();
   h_ElectronHoE_barrel_heep->Write();
   h_ElectronSigmaIEtaIEta_barrel_heep->Write();
   h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_barrel_heep->Write();
   h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_barrel_heep->Write(); 
   h_ElectronHcalIsoD2Heep_barrel_heep->Write();
   h_ElectronTrkIsoHeep_barrel_heep->Write();

   //endcap - all
   h_ElectronPt_endcap_all->Write();
   h_ElectronSCEta_fabs_endcap_all->Write();
   h_ElectronDeltaEtaTrkSC_endcap_all->Write();
   h_ElectronDeltaPhiTrkSC_endcap_all->Write();
   h_ElectronHoE_endcap_all->Write();
   h_ElectronSigmaIEtaIEta_endcap_all->Write();
   h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_endcap_all->Write();
   h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_endcap_all->Write(); 
   h_ElectronHcalIsoD2Heep_endcap_all->Write();
   h_ElectronTrkIsoHeep_endcap_all->Write();

   //endcap - heep
   h_ElectronPt_endcap_heep->Write();
   h_ElectronSCEta_fabs_endcap_heep->Write();
   h_ElectronDeltaEtaTrkSC_endcap_heep->Write();
   h_ElectronDeltaPhiTrkSC_endcap_heep->Write();
   h_ElectronHoE_endcap_heep->Write();
   h_ElectronSigmaIEtaIEta_endcap_heep->Write();
   h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_endcap_heep->Write();
   h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_endcap_heep->Write(); 
   h_ElectronHcalIsoD2Heep_endcap_heep->Write();
   h_ElectronTrkIsoHeep_endcap_heep->Write();


   //INFO
   //    //pT of both electrons, to be built using the histograms produced automatically by baseClass
   //    TH1F * h_pTElectrons = new TH1F ("h_pTElectrons","", getHistoNBins("pT1stEle"), getHistoMin("pT1stEle"), getHistoMax("pT1stEle"));
   //    h_pTElectrons->Add( & getHisto_noCuts_or_skim("pT1stEle") ); // all histos can be retrieved, see other getHisto_xxxx methods in baseClass.h
   //    h_pTElectrons->Add( & getHisto_noCuts_or_skim("pT2ndEle") );
   //    //one could also do:  *h_pTElectrons = getHisto_noCuts_or_skim("pT1stEle") + getHisto_noCuts_or_skim("pT2ndEle");
   //    h_pTElectrons->Write();
   //    //one could also do:   const TH1F& h = getHisto_noCuts_or_skim// and use h
   
   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
