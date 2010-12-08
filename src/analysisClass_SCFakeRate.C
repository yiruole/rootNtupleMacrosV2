/* This code chooses events that pass the criteria for the QCD fake rate calculation
It fills histograms to be used as the numerator (N_electrons) and denominator (N_superclusters)
First the electron collection is made, then the superclusters passing isolation, then the jet collection
with the isolated superclusters removed, then a NEW collection of superclusters is made which is the subset 
the isolated collection that are back to back with a jet (uncleaned of sc).
Finally the dR is calculated between the leading superclusters (up to 2) and the leading jets (cleaned of sc, 
up to 2) so that events can be rejected if the leadnig objects are too close.
*/


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

  CreateUserTH1D("SCPt_minus_ElePt", 100,0,100 );    
  CreateUserTH2D("SCPt_vs_ElePt", 100,0,100, 100,0,100 );    


   TH1F *h_TrigDiff = new TH1F("TrigDiff","TrigDiff",3.0,-1.5,1.5);
   TH2F *h2_DebugTrig = new TH2F("DebugTrig","DebugTrig;HLTResults;HLTBits",2,0,2,2,0,2);
   
   TH1F *h_dPhi_JetSC = new TH1F("dPhi_JetSC","dPhi_JetSC",650,0,6.5); h_dPhi_JetSC->Sumw2();
   TH1F *h_dR_EleSC = new TH1F("dR_EleSC","dR_EleSC",100,0,1.0); h_dR_EleSC->Sumw2();
   TH1F *h_dR_JetSC = new TH1F("dR_JetSC","dR_JetSC",600,0,3.0); h_dR_JetSC->Sumw2();
   TH1F *h_NisoSC = new TH1F ("NisoSC","NisoSC",6,-0.5,5.5);  h_NisoSC->Sumw2();
   TH1F *h_NisoSC_nearEle = new TH1F ("NisoSC_nearEle","NisoSC_nearEle",6,-0.5,5.5);  h_NisoSC_nearEle->Sumw2();

   TH1F *h_goodEleSCPt = new TH1F ("goodEleSCPt","goodEleSCPt",100,0,100); h_goodEleSCPt->Sumw2();
   TH1F *h_goodEleSCEta = new TH1F ("goodEleSCEta","goodEleSCEta",100,-3.,3.); h_goodEleSCEta->Sumw2();
   TH1F *h_goodEleSCPt_Barrel = new TH1F ("goodEleSCPt_Barrel","goodEleSCPt_Barrel",100,0,100); h_goodEleSCPt_Barrel->Sumw2();
   TH1F *h_goodEleSCPt_Endcap = new TH1F ("goodEleSCPt_Endcap","goodEleSCPt_Endcap",100,0,100); h_goodEleSCPt_Endcap->Sumw2();
   TH1F *h_goodEleSCPt_EcalIso_Barrel = new TH1F ("goodEleSCPt_EcalIso_Barrel","goodEleSCPt_EcalIso_Barrel",100,0,100); h_goodEleSCPt_EcalIso_Barrel->Sumw2();
   TH1F *h_goodEleSCPt_EcalIso_Endcap = new TH1F ("goodEleSCPt_EcalIso_Endcap","goodEleSCPt_EcalIso_Endcap",100,0,100); h_goodEleSCPt_EcalIso_Endcap->Sumw2();


   TH1F *h_goodEle_EoP = new TH1F ("goodEle_EoP","goodEle_EoP",500,0,5); h_goodEle_EoP->Sumw2();
   TH1F *h_goodEle_EoP_Barrel = new TH1F ("goodEle_EoP_Barrel","goodEle_EoP_Barrel",500,0,5); h_goodEle_EoP_Barrel->Sumw2();
   TH1F *h_goodEle_EoP_Endcap = new TH1F ("goodEle_EoP_Endcap","goodEle_EoP_Endcap",500,0,5); h_goodEle_EoP_Endcap->Sumw2();

   TH1F *h_goodSCPt = new TH1F ("goodSCPt","goodSCPt",100,0,100); h_goodSCPt->Sumw2();
   TH1F *h_goodSCEta = new TH1F ("goodSCEta","goodSCEta",100,-3.,3.); h_goodSCEta->Sumw2();
   TH1F *h_goodSCPt_Barrel = new TH1F ("goodSCPt_Barrel","goodSCPt_Barrel",100,0,100); h_goodSCPt_Barrel->Sumw2();
   TH1F *h_goodSCPt_Endcap = new TH1F ("goodSCPt_Endcap","goodSCPt_Endcap",100,0,100); h_goodSCPt_Endcap->Sumw2();
   TH1F *h_goodSCPt_EcalIso_Barrel = new TH1F ("goodSCPt_EcalIso_Barrel","goodSCPt_EcalIso_Barrel",100,0,100); h_goodSCPt_EcalIso_Barrel->Sumw2();
   TH1F *h_goodSCPt_EcalIso_Endcap = new TH1F ("goodSCPt_EcalIso_Endcap","goodSCPt_EcalIso_Endcap",100,0,100); h_goodSCPt_EcalIso_Endcap->Sumw2();

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

   TH1F *h_IsoSuperClusterHEEPEcalIso_barrel = new TH1F("IsoSuperClusterHEEPEcalIso_barrel","IsoSuperClusterHEEPEcalIso_barrel",100,0,10);
   TH1F *h_IsoSuperClusterEta_barrel = new TH1F("IsoSuperClusterEta_barrel","IsoSuperClusterEta_barrel",100,-3,3);
   TH1F *h_IsoSuperClusterHEEPTrkIso_barrel = new TH1F("IsoSuperClusterHEEPTrkIso_barrel","IsoSuperClusterHEEPTrkIso_barrel",100,0,15);
   TH1F *h_IsoSuperClusterHoE_barrel = new TH1F("IsoSuperClusterHoE_barrel","IsoSuperClusterHoE_barrel",200,0,0.2);
   TH1F *h_IsoSuperClusterPhi_barrel = new TH1F("IsoSuperClusterPhi_barrel","IsoSuperClusterPhi_barrel",100,-3.5,3.5);
   TH1F *h_IsoSuperClusterPt_barrel = new TH1F("IsoSuperClusterPt_barrel","IsoSuperClusterPt_barrel",100,0,400);
   TH1F *h_IsoSuperClusterS4S1_barrel = new TH1F("IsoSuperClusterS4S1_barrel","IsoSuperClusterS4S1_barrel",100,0,1.0);
   TH1F *h_IsoSuperClusterSigmaIEtaIEta_barrel = new TH1F("IsoSuperClusterSigmaIEtaIEta_barrel","IsoSuperClusterSigmaIEtaIEta_barrel",100,0,0.1);
   TH1F *h_IsoSuperClusterDrTrack1_barrel = new TH1F("IsoSuperClusterDrTrack1_barrel","IsoSuperClusterDrTrack1_barrel",300,0,3.0);

   TH1F *h_IsoSuperClusterHEEPEcalIso_endcap = new TH1F("IsoSuperClusterHEEPEcalIso_endcap","IsoSuperClusterHEEPEcalIso_endcap",100,0,10);
   TH1F *h_IsoSuperClusterEta_endcap = new TH1F("IsoSuperClusterEta_endcap","IsoSuperClusterEta_endcap",100,-3,3);
   TH1F *h_IsoSuperClusterHEEPTrkIso_endcap = new TH1F("IsoSuperClusterHEEPTrkIso_endcap","IsoSuperClusterHEEPTrkIso_endcap",100,0,15);
   TH1F *h_IsoSuperClusterHoE_endcap = new TH1F("IsoSuperClusterHoE_endcap","IsoSuperClusterHoE_endcap",200,0,0.2);
   TH1F *h_IsoSuperClusterPhi_endcap = new TH1F("IsoSuperClusterPhi_endcap","IsoSuperClusterPhi_endcap",100,-3.5,3.5);
   TH1F *h_IsoSuperClusterPt_endcap = new TH1F("IsoSuperClusterPt_endcap","IsoSuperClusterPt_endcap",100,0,400);
   TH1F *h_IsoSuperClusterS4S1_endcap = new TH1F("IsoSuperClusterS4S1_endcap","IsoSuperClusterS4S1_endcap",100,0,1.0);
   TH1F *h_IsoSuperClusterSigmaIEtaIEta_endcap = new TH1F("IsoSuperClusterSigmaIEtaIEta_endcap","IsoSuperClusterSigmaIEtaIEta_endcap",100,0,0.1);
   TH1F *h_IsoSuperClusterDrTrack1_endcap = new TH1F("IsoSuperClusterDrTrack1_endcap","IsoSuperClusterDrTrack1_endcap",300,0,3.0);

   TH2D *h_IsoSuperClusterPt_S4S1_barrel = new TH2D("IsoSuperClusterPt_S4S1_barrel","IsoSuperClusterPt_S4S1_barrel",100,0,400,110,0,1.1);
   TH2D *h_IsoSuperClusterPt_S4S1_endcap = new TH2D("IsoSuperClusterPt_S4S1_endcap","IsoSuperClusterPt_S4S1_endcap",100,0,400,110,0,1.1);
   TH2D *h_IsoSuperClusterPt_Phi_barrel = new TH2D("IsoSuperClusterPt_Phi_barrel","IsoSuperClusterPt_Phi_barrel",100,0,400,100,-3.5,3.5);
   TH2D *h_IsoSuperClusterPt_Phi_endcap = new TH2D("IsoSuperClusterPt_Phi_endcap","IsoSuperClusterPt_Phi_endcap",100,0,400,100,-3.5,3.5);
   TH2D *h_IsoSuperClusterEta_Phi = new TH2D("IsoSuperClusterEta_Phi","IsoSuperClusterEta_Phi",100,-3.5,3.5,100,-3.5,3.5);
  

   /////////initialize variables
   double FailRate = 0;
   double NFailHLT = 0;
   double HasSC = 0;

  
  ////////////////////// User's code to get preCut values - BEGIN ///////////////

  double eleEta_bar = getPreCutValue1("eleEta_bar");
  double eleEta_end_min = getPreCutValue1("eleEta_end");
  double eleEta_end_max = getPreCutValue2("eleEta_end");

  // Not used when using ElectronHeepID and heepBitMask // int eleIDType = (int) getPreCutValue1("eleIDType");
  int heepBitMask_EB  =  getPreCutValue1("heepBitMask_EBGapEE") ;
  int heepBitMask_GAP =  getPreCutValue2("heepBitMask_EBGapEE") ;

  int heepBitMask_EE  =  getPreCutValue3("heepBitMask_EBGapEE") ;

  int looseBitMask_EB       =  getPreCutValue1("looseBitMask_EBGapEE") ;
  int looseBitMask_GAP      =  getPreCutValue2("looseBitMask_EBGapEE") ;
  int looseBitMask_EE       =  getPreCutValue3("looseBitMask_EBGapEE") ;
  int looseBitMask_enabled  =  getPreCutValue4("looseBitMask_EBGapEE") ;

  double vertexMinimumNDOF = getPreCutValue1("vertexMinimumNDOF");
  double vertexMaxAbsZ = getPreCutValue1("vertexMaxAbsZ");
  double vertexMaxd0 = getPreCutValue1("vertexMaxd0");

  ////////////////////// User's code to get preCut values - END ///////////////

  Long64_t nentries = fChain->GetEntriesFast();
  STDOUT("analysisClass::Loop(): nentries = " << nentries);   
  
  ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
  ////// If the root version is updated and rootNtupleClass regenerated,     /////
  ////// these lines may need to be updated.                                 /////    
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) { // Begin of loop over events
  //for (Long64_t jentry=0; jentry<50;jentry++) { // Begin of loop over events
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry < 10 || jentry%10000 == 0) STDOUT("analysisClass::Loop(): jentry = " << jentry);   
    //if (Cut(ientry) < 0) continue;
    
    ////////////////////// User's code to be done for every event - BEGIN ///////////////////////

    //if (PtHat>=30) continue;   // to be used with samples that have lower PtHat cut only

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
      if ( !isData ) 
	{
	  PassTrig = 1; // for MC pass by default since most triggers are not in the currently used MC roottuples
	  if(jentry == 0 && i == 0) STDOUT("PassTrig set to 1 for MC for all events.");
	  continue;
	}
      if ( HLTFromRun[i] <= run ) {
 	if(jentry == 0 ) STDOUT("run, i, HLTTrigger[i], HLTFromRun[i] = "<<run<<"\t"<<i<<"\t"<<"\t"<<HLTTrigger[i]<<"\t"<<HLTFromRun[i]);
	if (HLTTrigger[i] > 0 && HLTTrigger[i] < HLTResults->size() ) {
	  PassTrig=HLTResults->at(HLTTrigger[i]);
	  HLTTrgUsed=HLTTrigger[i];
	} else {
	  STDOUT("ERROR: HLTTrigger out of range of HLTResults: HLTTrigger = "<<HLTTrigger[i] <<"and HLTResults size = "<< HLTResults->size());
	}
      }
    }
    if(jentry == 0 && isData) STDOUT("Run = "<<run <<", HLTTrgUsed is number = "<<HLTTrgUsed<<" of the list HLTPathsOfInterest");

    // Electrons
    // Choose the electrons for the numerator here
    vector<int> v_idx_ele_all;
    vector<int> v_idx_ele_PtCut;
    vector<int> v_idx_ele_HEEP;
    vector<int> v_idx_ele_HEEP_loose;
    int eleIDType = (int) getPreCutValue1("eleIDType");
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

	if ( (ElectronHeepID->at(iele) & ~heepBitMask)==0x0  && ElectronOverlaps->at(iele)==0 )
	  {
	    //STDOUT("ElectronHeepID = " << hex << ElectronHeepID->at(iele) << " ; ElectronPassID = " << ElectronPassID->at(iele) )
	    v_idx_ele_HEEP.push_back(iele);
	  }

	///// now look for loose eles
	    int looseBitMask;
	    if( fabs(ElectronEta->at(iele)) < eleEta_bar ) 
	      {
		looseBitMask = looseBitMask_EB;
	      }
	    else if ( fabs(ElectronEta->at(iele)) > eleEta_end_min && fabs(ElectronEta->at(iele)) < eleEta_end_max ) 
	      {
		looseBitMask = looseBitMask_EE;
	      }
	    else {
	      looseBitMask = looseBitMask_GAP;
	    }
	    if ( (ElectronHeepID->at(iele) & ~looseBitMask)==0x0  && ElectronOverlaps->at(iele)==0 )
	      {
		v_idx_ele_HEEP_loose.push_back(iele);
	      }

      } // End loop over electrons


    // Use the following code to apply HEEP cuts one by one, for example if you want to relax the pT cuts
// 	int eleID = ElectronPassID->at(iele);
// 	if ( (eleID & 1<<eleIDType) > 0  && ElectronOverlaps->at(iele)==0 )
// 	  {
// 	    v_idx_ele_HEEP.push_back(iele);
// 	  }

// 	bool HEEP = false;
// 	if (fabs(ElectronEta->at(iele))<eleEta_bar){
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
// 	if ((fabs(ElectronEta->at(iele))>eleEta_end_min) && (fabs(ElectronEta->at(iele))<eleEta_end_max)){
// 	if (fabs(ElectronDeltaEtaTrkSC->at(iele))<0.007){
// 	  if (fabs(ElectronDeltaPhiTrkSC->at(iele))<0.09){
// 	    if (ElectronHoE->at(iele)<0.05){
// 	      if (ElectronSigmaIEtaIEta->at(iele)<0.03){
// 		if (((ElectronEcalIsoHeep->at(iele) + ElectronHcalIsoD1Heep->at(iele))<eleEta_end_max) || 
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

// 	bool HEEP_loose = false;
// 	if (fabs(ElectronEta->at(iele))<eleEta_bar){
// 	if (fabs(ElectronDeltaEtaTrkSC->at(iele))<0.005){
// 	  if (fabs(ElectronDeltaPhiTrkSC->at(iele))<0.09){
// 	    if (ElectronHoE->at(iele)<0.05){
// 		      HEEP_loose = true;
// 	      }
// 	    }
// 	  }
// 	}

// 	if ((fabs(ElectronEta->at(iele))>eleEta_end_min) && (fabs(ElectronEta->at(iele))<eleEta_end_max)){
// 	if (fabs(ElectronDeltaEtaTrkSC->at(iele))<0.007){
// 	  if (fabs(ElectronDeltaPhiTrkSC->at(iele))<0.09){
// 	    if (ElectronHoE->at(iele)<0.05){
// 	      if (ElectronSigmaIEtaIEta->at(iele)<0.03){
// 		      HEEP_loose = true;
// 	      }
// 	    }
// 	  }
// 	}
// 	}

// 	if ( HEEP_loose  && ElectronOverlaps->at(iele)==0 )
// 	  {
// 	    v_idx_ele_HEEP_loose.push_back(iele);
// 	  }

//       } // End loop over electrons

    
    //////// Fill SC histos
    // Choose the superclusters for the denominator
    //  The superclusters are saved in the root tuples first from the barrel, 
    // then from the endcap.  So they are not ordered in pT.
    //  To deal with this, loop through all the superclusters and find the 2 highest pT
    //  then loop through the rest to get any other possible sc passing cuts.

     vector<int> v_idx_sc;

     ///Find two highest Pt SC, since they are not ordered 
     float scHighestPt = -99;
     float scNextPt = -99;
     int idx_scHighestPt = -1;
     int idx_scNextPt = -1;
     for(int isc=0;isc<SuperClusterPt->size();isc++){
      if ( SuperClusterPt->at(isc) < getPreCutValue1("ele_PtCut") ) continue;
      if ( 1 - SuperClusterS4S1->at(isc) > 0.95 ) continue; 
      if (SuperClusterHoE->at(isc)>0.05) continue;

      bool Barrel = false;
      bool Endcap = false;
      if (fabs(SuperClusterEta->at(isc))<1.442) Barrel = true;
      if ((fabs(SuperClusterEta->at(isc))<2.5)&&(fabs(SuperClusterEta->at(isc))>1.560)) Endcap = true;
      if (!Barrel && !Endcap) continue;

      if (Barrel && SuperClusterHEEPEcalIso->at(isc)>(6+(0.01*SuperClusterPt->at(isc)))) continue;
      if (Endcap && SuperClusterPt->at(isc)<50 && SuperClusterHEEPEcalIso->at(isc)>(6+(0.01*SuperClusterPt->at(isc)))) continue;
      if (Endcap && SuperClusterPt->at(isc)>=50 && SuperClusterHEEPEcalIso->at(isc)>(6+(0.01*(SuperClusterPt->at(isc)-50)))) continue;

      if (Barrel && SuperClusterSigmaIEtaIEta->at(isc)>0.014 ) continue;
      if (Endcap && SuperClusterSigmaIEtaIEta->at(isc)>0.035 ) continue;

      if (SuperClusterPt->at(isc)>scHighestPt){
	scNextPt = scHighestPt;
	idx_scNextPt = idx_scHighestPt;
	scHighestPt = SuperClusterPt->at(isc);
	idx_scHighestPt = isc;
      }
      else if (SuperClusterPt->at(isc)>scNextPt){
	scNextPt = SuperClusterPt->at(isc);
	idx_scNextPt = isc;
      }
     }
    if (idx_scHighestPt != -1) v_idx_sc.push_back(idx_scHighestPt);
    if (idx_scNextPt != -1)v_idx_sc.push_back(idx_scNextPt);

    //////now fill in the rest of the sc in whatever order they are
    for(int isc=0;isc<SuperClusterPt->size();isc++){
      if (isc==idx_scHighestPt || isc==idx_scNextPt) continue;
      if ( SuperClusterPt->at(isc) < getPreCutValue1("ele_PtCut") ) continue;
      if ( 1 - SuperClusterS4S1->at(isc) > 0.95 ) continue; 

      bool Barrel = false;
      bool Endcap = false;
      if (fabs(SuperClusterEta->at(isc))<1.442) Barrel = true;
      if ((fabs(SuperClusterEta->at(isc))<2.5)&&(fabs(SuperClusterEta->at(isc))>1.560)) Endcap = true;
      if (!Barrel && !Endcap) continue;

      if (Barrel && SuperClusterHEEPEcalIso->at(isc)>(6+(0.01*SuperClusterPt->at(isc)))) continue;
      if (Endcap && SuperClusterPt->at(isc)<50 && SuperClusterHEEPEcalIso->at(isc)>(6+(0.01*SuperClusterPt->at(isc)))) continue;
      if (Endcap && SuperClusterPt->at(isc)>=50 && SuperClusterHEEPEcalIso->at(isc)>(6+(0.01*(SuperClusterPt->at(isc)-50)))) continue;

      if (SuperClusterHoE->at(isc)<0.05) {
	v_idx_sc.push_back(isc);
      }
    }

    // Jets
    // Determine the jet collection that still has iso superclusters ("_PtCut") and cleaned of iso superclusters ("_PtCut_noOvrlap")
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
      } // jet loop


    vector <int> jetFlags(v_idx_jet_PtCut.size(), 0);
    int Njetflagged = 0;
    for (int isc=0; isc<v_idx_sc.size(); isc++)
      {
	TLorentzVector sc;
        sc.SetPtEtaPhiM(  SuperClusterPt->at(v_idx_sc[isc]),
			 SuperClusterEta->at(v_idx_sc[isc]),
			 SuperClusterPhi->at(v_idx_sc[isc]),0);
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
	if ( minDR < getPreCutValue1("jet_sc_DeltaRcut") && ijet_minDR > -1)
	  {
	    jetFlags[ijet_minDR] = 1;
	    Njetflagged++;
	  }
      } // sc loop

    for(int ijet=0; ijet<v_idx_jet_PtCut.size(); ijet++) //pT pre-cut + no overlaps with electrons + jetID
      {	
	bool passjetID = JetIdloose(CaloJetresEMF->at(ijet),CaloJetfHPD->at(ijet),CaloJetn90Hits->at(ijet), CaloJetEta->at(ijet));
	if( jetFlags[ijet] == 0                           /* NO overlap with sc */  
	    && passjetID == true )                            /* pass JetID */
	  // && (caloJetOverlaps[ijet] & 1 << 5)==0 )         /* NO overlap with muons */      
	  v_idx_jet_PtCut_noOvrlap.push_back(ijet);
	
      } // End loop over jets


    // make new collection of sc that have a jet (including iso sc) back to back.
     vector<int> v_idx_sc_dPhi;
    for(int isc=0;isc<v_idx_sc.size();isc++)
	   {
	     //Require dPhi>3.0 between SC and JET
	     TVector3 sc_vec;
	     sc_vec.SetPtEtaPhi(SuperClusterPt->at(v_idx_sc[isc]),
			   SuperClusterEta->at(v_idx_sc[isc]),
			   SuperClusterPhi->at(v_idx_sc[isc]));
	     double dPhi_SC_Jet=0;
	     for (int ijet=0;ijet<v_idx_jet_PtCut.size();ijet++){
	       TVector3 jet_vec;
	       jet_vec.SetPtEtaPhi(CaloJetPt->at(v_idx_jet_PtCut[ijet]),
			    CaloJetEta->at(v_idx_jet_PtCut[ijet]),
			    CaloJetPhi->at(v_idx_jet_PtCut[ijet]));
	       double deltaPhi=fabs(jet_vec.DeltaPhi(sc_vec));
	       if (deltaPhi>dPhi_SC_Jet)dPhi_SC_Jet=deltaPhi;
	     }
	     if (dPhi_SC_Jet>=3.0) v_idx_sc_dPhi.push_back(v_idx_sc[isc]);
	     h_dPhi_JetSC->Fill(dPhi_SC_Jet);	
	   }

    /// Make plots to make sure each sc matches only one electron
    for (int iele=0; iele<v_idx_ele_HEEP.size(); iele++){
    double minDr_Ele_Sc = 99;
    vector<int> closeSC;
	TLorentzVector ele_vec;
        ele_vec.SetPtEtaPhiM(  ElectronPt->at(v_idx_ele_HEEP[iele]),
			       ElectronEta->at(v_idx_ele_HEEP[iele]),
			       ElectronPhi->at(v_idx_ele_HEEP[iele]),0);
	for (int isc=0; isc<v_idx_sc_dPhi.size(); isc++){
	  TLorentzVector sc;
	  sc.SetPtEtaPhiM(  SuperClusterPt->at(v_idx_sc_dPhi[isc]),
			 SuperClusterEta->at(v_idx_sc_dPhi[isc]),
			 SuperClusterPhi->at(v_idx_sc_dPhi[isc]),0);
	  double dR = ele_vec.DeltaR(sc);
	  if (dR<minDr_Ele_Sc) minDr_Ele_Sc=dR;
	  if (dR<0.3) closeSC.push_back(iele);
	}
	h_dR_EleSC->Fill(minDr_Ele_Sc);
	h_NisoSC_nearEle->Fill(closeSC.size());
    }


     //find dR between leading isoSC and leading jets to reject events where they are too close
    float smallest_ScJet_dR = 99;
    for (int ijet=0; ijet<v_idx_jet_PtCut_noOvrlap.size(); ijet++){
      if (ijet>1) break; // only care about two leading jets
      TVector3 jet_vec;
      jet_vec.SetPtEtaPhi(CaloJetPt->at(v_idx_jet_PtCut_noOvrlap[ijet]),
			       CaloJetEta->at(v_idx_jet_PtCut_noOvrlap[ijet]),
			       CaloJetPhi->at(v_idx_jet_PtCut_noOvrlap[ijet]));
      float minDR = 99;
      int idx_nearest_sc = -1;
      for(int isc=0;isc<v_idx_sc_dPhi.size();isc++)
	  {
	    if (isc>1) break; //only care about two leading SC
	    TVector3 sc_vec;
	    sc_vec.SetPtEtaPhi(SuperClusterPt->at(v_idx_sc_dPhi[isc]),
				  SuperClusterEta->at(v_idx_sc_dPhi[isc]),
				  SuperClusterPhi->at(v_idx_sc_dPhi[isc]));
	    float DR = jet_vec.DeltaR(sc_vec);
	    if (DR<minDR) {
		 minDR = DR;
		 idx_nearest_sc = isc;
	      }
	  }
       if (minDR!=99) h_dR_JetSC->Fill(minDR);
       if (minDR<smallest_ScJet_dR) smallest_ScJet_dR= minDR;
    }

    // vertexes
    vector<int> v_idx_vertex_good;
    // loop over vertexes
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
    

    // Set the value of the variableNames listed in the cutFile to their current value
    
    // HLT
    fillVariableWithValue( "PassHLT", PassTrig ) ;

    // N_goodVertex
    fillVariableWithValue( "nVertex_good", v_idx_vertex_good.size() ) ;

    // dR SC-Jet
    fillVariableWithValue( "dR_SCJet", smallest_ScJet_dR );

    //## SC
    fillVariableWithValue( "nIsoSC", v_idx_sc.size() );

    // nEle
    fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp", v_idx_ele_HEEP.size() ) ;
     
    // nJet
    fillVariableWithValue( "nJet_PtCut", v_idx_jet_PtCut_noOvrlap.size() ) ;

    // dPhi SC-Jet
    fillVariableWithValue("nSC_BackToBack", v_idx_sc_dPhi.size());

    // 1st SC
    if( v_idx_sc_dPhi.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stSC_ISO", SuperClusterPt->at(v_idx_sc_dPhi[0]) );
	fillVariableWithValue( "Eta1stSC_ISO", SuperClusterEta->at(v_idx_sc_dPhi[0]) );
      }

    // 1st jet
    if( v_idx_jet_PtCut_noOvrlap.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stJet_noOvrlpSC", CaloJetPt->at(v_idx_jet_PtCut[0]) );
	fillVariableWithValue( "Eta1stJet_noOvrlpSC", CaloJetEta->at(v_idx_jet_PtCut[0]) );
      }

    // Mee
    float MEE = 0.05;  // set slightly larger than 0 so that events with 1 sc still pass the 0<MEE<60 cut
    if (v_idx_sc_dPhi.size()>= 2)
      {
	TLorentzVector ele1, ele2, ee;
	ele1.SetPtEtaPhiM(SuperClusterPt->at(v_idx_sc_dPhi[0]),
			  SuperClusterEta->at(v_idx_sc_dPhi[0]),
			  SuperClusterPhi->at(v_idx_sc_dPhi[0]),0);
	ele2.SetPtEtaPhiM(SuperClusterPt->at(v_idx_sc_dPhi[1]),
			  SuperClusterEta->at(v_idx_sc_dPhi[1]),
			  SuperClusterPhi->at(v_idx_sc_dPhi[1]),0);
	ee = ele1+ele2;
	MEE=ee.M();
      }
	fillVariableWithValue("Mee",MEE);

    // MET
    fillVariableWithValue("MET",PFMET->at(0));

    // Evaluate cuts (but do not apply them)
    evaluateCuts();

     if (passedCut("nIsoSC")){
       HasSC++;
     }

     // Fill histograms and do analysis based on cut evaluation
     if (v_idx_sc_dPhi.size()>0 && !passedCut("PassHLT")){
       NFailHLT++;
       h_eta_failHLT->Fill(SuperClusterEta->at(v_idx_sc_dPhi[0]));
       h_phi_failHLT->Fill(SuperClusterPhi->at(v_idx_sc_dPhi[0]));
     }


     /// Fill fake rate plots
     ///////////////////////////////////////
     if (passedCut("0")){  // use cut file to select cuts used for fake rate events
	 for(int isc=0;isc<v_idx_sc_dPhi.size();isc++)
	   {
	     h_goodSCPt->Fill(SuperClusterPt->at(v_idx_sc_dPhi[isc]));
	     h_goodSCEta->Fill(SuperClusterEta->at(v_idx_sc_dPhi[isc]));
	     h_loose_goodSCPt->Fill(SuperClusterPt->at(v_idx_sc_dPhi[isc]));
	     h_loose_goodSCEta->Fill(SuperClusterEta->at(v_idx_sc_dPhi[isc]));

	     bool Barrel = false;
	     bool Endcap = false;
	     if (fabs(SuperClusterEta->at(v_idx_sc_dPhi[isc]))<1.45) Barrel = true;
	     if (fabs(SuperClusterEta->at(v_idx_sc_dPhi[isc]))>eleEta_end_min && fabs(SuperClusterEta->at(v_idx_sc_dPhi[isc]))<eleEta_end_max) Endcap = true;
	     
	     //Fill the denominator histograms and all variable histograms
	     h_IsoSuperClusterEta_Phi->Fill(SuperClusterEta->at(v_idx_sc_dPhi[isc]),SuperClusterPhi->at(v_idx_sc_dPhi[isc]));
	     if (Barrel) {
	       h_goodSCPt_Barrel->Fill(SuperClusterPt->at(v_idx_sc_dPhi[isc]));
	       h_loose_goodSCPt_Barrel->Fill(SuperClusterPt->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterHEEPEcalIso_barrel->Fill(SuperClusterHEEPEcalIso->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterEta_barrel->Fill(SuperClusterEta->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterHEEPTrkIso_barrel->Fill(SuperClusterHEEPTrkIso->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterHoE_barrel->Fill(SuperClusterHoE->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterPhi_barrel->Fill(SuperClusterPhi->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterPt_barrel->Fill(SuperClusterPt->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterS4S1_barrel->Fill(SuperClusterS4S1->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterSigmaIEtaIEta_barrel->Fill(SuperClusterSigmaIEtaIEta->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterDrTrack1_barrel->Fill(SuperClusterDrTrack1->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterPt_S4S1_barrel->Fill(SuperClusterPt->at(v_idx_sc_dPhi[isc]),SuperClusterS4S1->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterPt_Phi_barrel->Fill(SuperClusterPt->at(v_idx_sc_dPhi[isc]),SuperClusterPhi->at(v_idx_sc_dPhi[isc]));
	     }

	     if (Endcap){
	       h_goodSCPt_Endcap->Fill(SuperClusterPt->at(v_idx_sc_dPhi[isc]));
	       h_loose_goodSCPt_Endcap->Fill(SuperClusterPt->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterHEEPEcalIso_endcap->Fill(SuperClusterHEEPEcalIso->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterEta_endcap->Fill(SuperClusterEta->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterHEEPTrkIso_endcap->Fill(SuperClusterHEEPTrkIso->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterHoE_endcap->Fill(SuperClusterHoE->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterPhi_endcap->Fill(SuperClusterPhi->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterPt_endcap->Fill(SuperClusterPt->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterS4S1_endcap->Fill(SuperClusterS4S1->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterSigmaIEtaIEta_endcap->Fill(SuperClusterSigmaIEtaIEta->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterDrTrack1_endcap->Fill(SuperClusterDrTrack1->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterPt_S4S1_endcap->Fill(SuperClusterPt->at(v_idx_sc_dPhi[isc]),SuperClusterS4S1->at(v_idx_sc_dPhi[isc]));
	       h_IsoSuperClusterPt_Phi_endcap->Fill(SuperClusterPt->at(v_idx_sc_dPhi[isc]),SuperClusterPhi->at(v_idx_sc_dPhi[isc]));
	     }

	       TVector3 sc_vec;
	       sc_vec.SetPtEtaPhi(SuperClusterPt->at(v_idx_sc_dPhi[isc]),
			   SuperClusterEta->at(v_idx_sc_dPhi[isc]),
			   SuperClusterPhi->at(v_idx_sc_dPhi[isc]));

	     ///  see if there is a HEEP ele to match this superclusters
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
	     if (deltaR_ele_sc<0.3){  // Fill the numerator histograms
// 	       STDOUT("FAKERATE: ------------------- Run, LS, Event = "<<run<<", "<<ls<<", "<<event);
// 	       STDOUT("FAKERATE: Barrel/Endcap, isc, v_idx_sc_dPhi[isc], SuperClusterPt, Eta, Phi = "<< Barrel<<"/"<<Endcap<<", "<<isc <<", "<< v_idx_sc_dPhi[isc] <<", "<< SuperClusterPt->at(v_idx_sc_dPhi[isc]) <<", "<< SuperClusterEta->at(v_idx_sc_dPhi[isc]) <<", "<< SuperClusterPhi->at(v_idx_sc_dPhi[isc]));
		h_goodEleSCPt->Fill(ElectronSCPt->at(v_idx_ele_HEEP[idx_HEEP]));
		h_goodEleSCEta->Fill(ElectronSCEta->at(v_idx_ele_HEEP[idx_HEEP]));
		
		bool Barrel = false;
		bool Endcap = false;
		if (fabs(ElectronSCEta->at(v_idx_ele_HEEP[idx_HEEP]))<eleEta_bar) Barrel = true;
		if (fabs(ElectronSCEta->at(v_idx_ele_HEEP[idx_HEEP]))>eleEta_end_min && fabs(ElectronSCEta->at(v_idx_ele_HEEP[idx_HEEP]))<eleEta_end_max) Endcap = true;
		//Fill EoP histograms to look for real electron contamination
		double E = ElectronPt->at(v_idx_ele_HEEP[idx_HEEP]);
		double P = ElectronTrackPt->at(v_idx_ele_HEEP[idx_HEEP]);
		double EoP = -99;
		if (P!=0) {
		  EoP = E / P;
		  h_goodEle_EoP->Fill(EoP);
		  if (Barrel) h_goodEle_EoP_Barrel->Fill(EoP);
		  if (Endcap) h_goodEle_EoP_Endcap->Fill(EoP);
		}
// 	       STDOUT("FAKERATE: Barrel/Endcap, idx_HEEP, v_idx_ele_HEEP[idx_HEEP], ElectronSCPt, Eta, Phi = "<< Barrel<<"/"<<Endcap <<", "<<idx_HEEP <<", "<< v_idx_ele_HEEP[idx_HEEP] <<", "<< ElectronSCPt->at(v_idx_ele_HEEP[idx_HEEP]) <<", "<< ElectronSCEta->at(v_idx_ele_HEEP[idx_HEEP]) <<", "<< ElectronSCPhi->at(v_idx_ele_HEEP[idx_HEEP]) );
	       FillUserTH1D("SCPt_minus_ElePt", ElectronSCPt->at(v_idx_ele_HEEP[idx_HEEP])-SuperClusterPt->at(v_idx_sc_dPhi[isc]) );    
	       FillUserTH2D("SCPt_vs_ElePt", ElectronSCPt->at(v_idx_ele_HEEP[idx_HEEP]), SuperClusterPt->at(v_idx_sc_dPhi[isc]) );    

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
		if (fabs(ElectronSCEta->at(v_idx_ele_HEEP_loose[idx_HEEP_loose]))<eleEta_bar) Barrel = true;
		if (fabs(ElectronSCEta->at(v_idx_ele_HEEP_loose[idx_HEEP_loose]))>eleEta_end_min && fabs(ElectronSCEta->at(v_idx_ele_HEEP_loose[idx_HEEP_loose]))<eleEta_end_max) Endcap = true;
		
		if (Barrel) h_loose_goodEleSCPt_Barrel->Fill(ElectronSCPt->at(v_idx_ele_HEEP_loose[idx_HEEP_loose]));
		if (Endcap) h_loose_goodEleSCPt_Endcap->Fill(ElectronSCPt->at(v_idx_ele_HEEP_loose[idx_HEEP_loose]));
		}
	   } // for sc

	 } // if passed cuts "ALL"
    ////////////////////// User's code to be done for every event - END ///////////////////////
  } // End of loop over events
  

  ////////////////////// User's code to write histos - BEGIN ///////////////////////

   h_TrigDiff->Write();
   h2_DebugTrig->Write();

   h_dPhi_JetSC->Write();
   h_dR_JetSC->Write();
   h_dR_EleSC->Write();
   h_NisoSC->Write();
   h_NisoSC_nearEle->Write();

   h_goodEleSCPt->Write();
   h_goodEleSCEta->Write();
   h_goodEleSCPt_Barrel->Write();
   h_goodEleSCPt_Endcap->Write();
   h_goodEle_EoP->Write();
   h_goodEle_EoP_Barrel->Write();
   h_goodEle_EoP_Endcap->Write();
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

  h_IsoSuperClusterHEEPEcalIso_barrel->Write();
  h_IsoSuperClusterEta_barrel->Write(); 
  h_IsoSuperClusterHEEPTrkIso_barrel->Write();
  h_IsoSuperClusterHoE_barrel->Write();
  h_IsoSuperClusterPhi_barrel->Write();
  h_IsoSuperClusterPt_barrel->Write();
  h_IsoSuperClusterS4S1_barrel->Write();
  h_IsoSuperClusterSigmaIEtaIEta_barrel->Write();
  h_IsoSuperClusterDrTrack1_barrel->Write();

  h_IsoSuperClusterHEEPEcalIso_endcap->Write();
  h_IsoSuperClusterEta_endcap->Write(); 
  h_IsoSuperClusterHEEPTrkIso_endcap->Write();
  h_IsoSuperClusterHoE_endcap->Write();
  h_IsoSuperClusterPhi_endcap->Write();
  h_IsoSuperClusterPt_endcap->Write();
  h_IsoSuperClusterS4S1_endcap->Write();
  h_IsoSuperClusterSigmaIEtaIEta_endcap->Write();
  h_IsoSuperClusterDrTrack1_endcap->Write();

  h_IsoSuperClusterPt_S4S1_barrel->Write();
  h_IsoSuperClusterPt_S4S1_endcap->Write();
  h_IsoSuperClusterPt_Phi_barrel->Write();
  h_IsoSuperClusterPt_Phi_endcap->Write();
  h_IsoSuperClusterEta_Phi->Write();

   FailRate = 100 * NFailHLT/HasSC;
   //cout << "NFail: " << NFailHLT << "\t" << "FailRate: " << FailRate << " %" << endl;
  
  //STDOUT("analysisClass::Loop() ends");   
}
