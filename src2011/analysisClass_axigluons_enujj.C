#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

double getPz( TLorentzVector lepton, double pxNeutrino, double pyNeutrino, double massW) {
 double pn=0., app=0., pznp=0., pznm=0.;
 double a=0., b=0., c=0.;

 app = pow(lepton.E(),2)+pow(pxNeutrino,2)+pow(pyNeutrino,2)-pow(lepton.Px()+pxNeutrino,2)-pow(lepton.Py()+pyNeutrino,2)-pow(lepton.Pz(),2)-pow(massW,2);
 a= pow(lepton.E(),2)-pow(lepton.Pz(),2);
 b= lepton.Pz()*app;
 c= ( pow(pxNeutrino,2)+pow(pyNeutrino,2) )*pow(lepton.E(),2) - pow(app,2)/4.;

 pznp = ( -b + sqrt( pow(b,2)-4.*a*c ) )/(2.*a);
 pznm = ( -b - sqrt( pow(b,2)-4.*a*c ) )/(2.*a);
 if ( pow(b,2)-4.*a*c < 0. ){ pznp=-b/(2.*a); pznm =-b/(2.*a); }

 if( fabs(pznp) < fabs(pznm) ){
   pn=pznp; }
 else{ pn=pznm;  }
 return pn;
}

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile){}

analysisClass::~analysisClass(){}

void analysisClass::Loop()
{
   std::cout << "analysisClass::Loop() begins" <<std::endl;   

   //--------------------------------------------------------------------------
   // Decide which plots to save (default is to save everything)
   //--------------------------------------------------------------------------
   
   fillSkim                         (  true  ) ;
   fillAllPreviousCuts              ( !true  ) ;
   fillAllOtherCuts                 (  true  ) ;
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------

   CreateUserTH1D( "Ele1_Pt"	           , 	getHistoNBins("Ele1_Pt"), getHistoMin("Ele1_Pt"), getHistoMax("Ele1_Pt")     ) ; 
   CreateUserTH1D( "Ele1_Eta"	           , 	getHistoNBins("Ele1_Eta"), getHistoMin("Ele1_Eta"), getHistoMax("Ele1_Eta")     ) ; 
   CreateUserTH1D( "Ele1_Phi"	           , 	60  , -3.1416 , +3.1416	 ) ; 
   CreateUserTH1D( "Ele1_Charge"	   , 	2   , -1.0001 , 1.0001	 ) ; 

   CreateUserTH1D( "MET_Pt"                ,    getHistoNBins("MET_Pt"), getHistoMin("MET_Pt"), getHistoMax("MET_Pt")     ) ; 
   CreateUserTH1D( "MET_Phi"		   , 	60  , -3.1416 , +3.1416	 ) ; 

   CreateUserTH1D( "Jet1_Pt"               ,    getHistoNBins("Jet1_Pt"), getHistoMin("Jet1_Pt"), getHistoMax("Jet1_Pt")     ) ; 
   CreateUserTH1D( "Jet1_Eta"	           , 	getHistoNBins("Jet1_Eta"), getHistoMin("Jet1_Eta"), getHistoMax("Jet1_Eta")     ) ; 
   CreateUserTH1D( "Jet1_Phi"	           , 	60  , -3.1416 , +3.1416	 ) ;  
   CreateUserTH1D( "Jet1_btagTCHE"         ,    250 , 0       , 50	 ) ; 
   CreateUserTH1D( "Jet2_Pt"               ,    getHistoNBins("Jet2_Pt"), getHistoMin("Jet2_Pt"), getHistoMax("Jet2_Pt")     ) ; 
   CreateUserTH1D( "Jet2_Eta"	           , 	getHistoNBins("Jet2_Eta"), getHistoMin("Jet2_Eta"), getHistoMax("Jet2_Eta")     ) ; 
   CreateUserTH1D( "Jet2_Phi"	           , 	60  , -3.1416 , +3.1416	 ) ;  
   CreateUserTH1D( "Jet2_btagTCHE"         ,    250 , 0       , 50	 ) ; 

   CreateUserTH1D( "nEle"                  ,    getHistoNBins("nEle"), getHistoMin("nEle"), getHistoMax("nEle")     ) ; 
   CreateUserTH1D( "nMuon"                 ,    getHistoNBins("nMuon"), getHistoMin("nMuon"), getHistoMax("nMuon")     ) ; 
   CreateUserTH1D( "nJet"                  ,    getHistoNBins("nJet"), getHistoMin("nJet"), getHistoMax("nJet")     ) ; 
   CreateUserTH1D( "nJet_btagTCHE"         ,    getHistoNBins("nJet"), getHistoMin("nJet"), getHistoMax("nJet")     ) ; 
   CreateUserTH1D( "nVertex"               ,    31   , -0.5   , 30.5	 ) ; 
   CreateUserTH1D( "nVertex_good"          ,    31   , -0.5   , 30.5	 ) ; 

   CreateUserTH1D( "DR_Ele1Jet1"	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele1Jet2"	   , 	getHistoNBins("DR_Ele1Jet2"), getHistoMin("DR_Ele1Jet2"), getHistoMax("DR_Ele1Jet2")     ) ; 
   CreateUserTH1D( "DR_Jet1Jet2"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
   CreateUserTH1D( "mDPhi_METEle1"	   , 	getHistoNBins("mDPhi_METEle1"), getHistoMin("mDPhi_METEle1"), getHistoMax("mDPhi_METEle1")     ) ; 
   CreateUserTH1D( "mDPhi_METJet1"	   , 	getHistoNBins("mDPhi_METJet1"), getHistoMin("mDPhi_METJet1"), getHistoMax("mDPhi_METJet1")     ) ; 
   CreateUserTH1D( "mDPhi_METJet2"	   , 	getHistoNBins("mDPhi_METJet2"), getHistoMin("mDPhi_METJet2"), getHistoMax("mDPhi_METJet2")     ) ; 
   CreateUserTH1D( "mDEta_Jet1Jet2"	   , 	getHistoNBins("mDEta_Jet1Jet2"), getHistoMin("mDEta_Jet1Jet2"), getHistoMax("mDEta_Jet1Jet2")   ) ; 
   CreateUserTH1D( "mDPhi_Jet1Jet2"	   , 	getHistoNBins("mDPhi_Jet1Jet2"), getHistoMin("mDPhi_Jet1Jet2"), getHistoMax("mDPhi_Jet1Jet2")     ) ; 
   CreateUserTH2D( "DR_Ele1Jet2_vs_DR_Ele1Jet1"	   
		   , getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     
		   , getHistoNBins("DR_Ele1Jet2"), getHistoMin("DR_Ele1Jet2"), getHistoMax("DR_Ele1Jet2")     
		   ) ;   
   CreateUserTH1D( "mDphi_BosonJet1"	   , 	getHistoNBins("mDphi_BosonJet1"), getHistoMin("mDphi_BosonJet1"), getHistoMax("mDphi_BosonJet1")     ) ; 
   CreateUserTH1D( "mDphi_BosonJet2"	   , 	getHistoNBins("mDphi_BosonJet2"), getHistoMin("mDphi_BosonJet2"), getHistoMax("mDphi_BosonJet2")     ) ; 
   CreateUserTH1D( "DR_BosonJet1"	   , 	getHistoNBins("DR_BosonJet1"), getHistoMin("DR_BosonJet1"), getHistoMax("DR_BosonJet1")     ) ; 
   CreateUserTH1D( "DR_BosonJet2"	   , 	getHistoNBins("DR_BosonJet2"), getHistoMin("DR_BosonJet2"), getHistoMax("DR_BosonJet2")     ) ; 
   CreateUserTH2D( "DR_BosonJet2_vs_DR_BosonJet1"	   
		   , getHistoNBins("DR_BosonJet1"), getHistoMin("DR_BosonJet1"), getHistoMax("DR_BosonJet1")     
		   , getHistoNBins("DR_BosonJet2"), getHistoMin("DR_BosonJet2"), getHistoMax("DR_BosonJet2")     
		   ) ;   
   CreateUserTH1D( "MT_Ele1MET"	           , 	getHistoNBins("MT_Ele1MET"), getHistoMin("MT_Ele1MET"), getHistoMax("MT_Ele1MET")     ) ; 
   CreateUserTH1D( "Pt_Ele1MET"	           , 	getHistoNBins("Pt_Ele1MET"), getHistoMin("Pt_Ele1MET"), getHistoMax("Pt_Ele1MET")     ) ; 
   CreateUserTH1D( "Pt_j1j2"	           , 	getHistoNBins("Pt_j1j2"), getHistoMin("Pt_j1j2"), getHistoMax("Pt_j1j2")     ) ; 
   CreateUserTH1D( "sT_enujj"	           , 	getHistoNBins("sT_enujj"), getHistoMin("sT_enujj"), getHistoMax("sT_enujj")     ) ; 
   CreateUserTH1D( "M_j1j2"	           , 	getHistoNBins("M_j1j2"), getHistoMin("M_j1j2"), getHistoMax("M_j1j2")     ) ; 
   CreateUserTH1D( "Jet1_Pt_over_Mj1j2"	   , 	getHistoNBins("Jet1_Pt_over_Mj1j2"), getHistoMin("Jet1_Pt_over_Mj1j2"), getHistoMax("Jet1_Pt_over_Mj1j2") ) ; 
   CreateUserTH1D( "Jet2_Pt_over_Mj1j2"	   , 	getHistoNBins("Jet2_Pt_over_Mj1j2"), getHistoMin("Jet2_Pt_over_Mj1j2"), getHistoMax("Jet2_Pt_over_Mj1j2") ) ; 
   CreateUserTH1D( "Jet1_E_over_Mj1j2"	   , 	getHistoNBins("Jet1_E_over_Mj1j2"), getHistoMin("Jet1_E_over_Mj1j2"), getHistoMax("Jet1_E_over_Mj1j2") ) ; 
   CreateUserTH1D( "Jet2_E_over_Mj1j2"	   , 	getHistoNBins("Jet2_E_over_Mj1j2"), getHistoMin("Jet2_E_over_Mj1j2"), getHistoMax("Jet2_E_over_Mj1j2") ) ; 
   CreateUserTH1D( "Jet1_E_over_Jet1_Pt"   , 	getHistoNBins("Jet1_E_over_Jet1_Pt"), getHistoMin("Jet1_E_over_Jet1_Pt"), getHistoMax("Jet1_E_over_Jet1_Pt") ) ; 
   CreateUserTH1D( "Jet2_E_over_Jet2_Pt"   , 	getHistoNBins("Jet2_E_over_Jet2_Pt"), getHistoMin("Jet2_E_over_Jet2_Pt"), getHistoMax("Jet2_E_over_Jet2_Pt") ) ; 


   //--------------------------------------------------------------------------
   // Precuts
   //--------------------------------------------------------------------------

   double Wmass = getPreCutValue1("Wmass");
   
   //--------------------------------------------------------------------------
   // Loop over the chain
   //--------------------------------------------------------------------------

   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntries();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << "/" << nentries << std::endl;   

     //--------------------------------------------------------------------------
     // Reset the cuts
     //--------------------------------------------------------------------------

     resetCuts();

     //--------------------------------------------------------------------------
     // Check good run list
     //--------------------------------------------------------------------------
     
     int    passedJSON = passJSON ( run, ls , isData ) ;


     //--------------------------------------------------------------------------
     // Event-weight
     //--------------------------------------------------------------------------
     
     double event_weight = 1;

     //--------------------------------------------------------------------------
     // Do pileup re-weighting
     //--------------------------------------------------------------------------
     
     int NPILEUP_AVE = int( (nPileUpInt_BXminus1 + nPileUpInt_BX0 + nPileUpInt_BXplus1)/3 );
     int NPILEUP_FINAL = min( NPILEUP_AVE , 25 );
     event_weight = getPileupWeight ( NPILEUP_FINAL, isData ) ;
     //event_weight = getPileupWeight ( min(nPileUpInt_BX0,25), isData ) ;

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "PassJSON"                      , passedJSON ) ; 

     // Noise filters
     fillVariableWithValue(   "PassHBHENoiseFilter"           , PassHBHENoiseFilter ) ; 
     fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight ) ; 
     fillVariableWithValue(   "PassTrackingFailure"           , PassTrackingFailure ) ; 


     // Electrons
     fillVariableWithValue(   "nEle"                          , nEle ) ;

     TLorentzVector Wboson;       

     if ( nEle >= 1 ) { 
       fillVariableWithValue( "Ele1_Pt"                       , Ele1_Pt ) ;
       fillVariableWithValue( "Ele1_Eta"                      , Ele1_Eta ) ;

       // Get Pz of W boson
       TVector2 Wboson2D;       
       TLorentzVector lepton;
       TVector2 lepton2D;       
       TVector2 neutrino2D;       

       lepton.SetPtEtaPhiE(Ele1_Pt, Ele1_Eta, Ele1_Phi, Ele1_Energy);       
       lepton2D.SetMagPhi(Ele1_Pt, Ele1_Phi);
       neutrino2D.SetMagPhi(MET_Pt, MET_Phi);
       Wboson2D = lepton2D + neutrino2D;  
       double PzW = getPz( lepton, neutrino2D.Px(), neutrino2D.Py(), Wmass);
       double EnergyW = sqrt( Wboson2D.Px()*Wboson2D.Px() + Wboson2D.Py()*Wboson2D.Py() + PzW*PzW + Wmass*Wmass );
       Wboson.SetPxPyPzE( Wboson2D.Px() , Wboson2D.Py() , PzW , EnergyW );       
     }

     // MET variables
     fillVariableWithValue(   "MET_Pt"                        , MET_Pt ) ;
     if ( nEle >= 1 ) { 
       fillVariableWithValue( "MT_Ele1MET"                    , MT_Ele1MET ) ;
       fillVariableWithValue( "Pt_Ele1MET"                    , Pt_Ele1MET ) ;
     }
     
     // Jets
     fillVariableWithValue(   "nJet"                          , nJet ) ;
     if ( nJet >= 1 ) { 
       fillVariableWithValue( "Jet1_Pt"                       , Jet1_Pt ) ;
       fillVariableWithValue( "Jet1_Eta"                      , Jet1_Eta ) ;
       fillVariableWithValue( "Jet1_btagTCHE"                 , Jet1_btagTCHE ) ;
       fillVariableWithValue( "Jet1_E_over_Jet1_Pt"           , Jet1_Energy / Jet1_Pt ) ;       
     }
     if ( nJet >= 2 ) { 
       fillVariableWithValue( "Jet2_Pt"                       , Jet2_Pt ) ;
       fillVariableWithValue( "Jet2_Eta"                      , Jet2_Eta ) ;
       fillVariableWithValue( "Jet2_btagTCHE"                 , Jet2_btagTCHE ) ;
       fillVariableWithValue( "Pt_j1j2"                       , Pt_j1j2 ) ;
       fillVariableWithValue( "M_j1j2"                        , M_j1j2 ) ;
       fillVariableWithValue( "M_j1j2_cut"                    , M_j1j2 ) ;
       fillVariableWithValue( "mDEta_Jet1Jet2"                , fabs(Jet1_Eta - Jet2_Eta) ) ;
       fillVariableWithValue( "Jet1_Pt_over_Mj1j2"            , Jet1_Pt / M_j1j2 ) ;
       fillVariableWithValue( "Jet2_Pt_over_Mj1j2"            , Jet2_Pt / M_j1j2 ) ;
       fillVariableWithValue( "Jet1_E_over_Mj1j2"             , Jet1_Energy / M_j1j2 ) ;
       fillVariableWithValue( "Jet2_E_over_Mj1j2"             , Jet2_Energy / M_j1j2 ) ;
       fillVariableWithValue( "Jet2_E_over_Jet2_Pt"           , Jet2_Energy / Jet2_Pt ) ;

       TVector2 jet1, jet2;
       jet1.SetMagPhi( Jet1_Pt , Jet1_Phi);
       jet2.SetMagPhi( Jet2_Pt , Jet2_Phi);
       fillVariableWithValue( "mDPhi_Jet1Jet2"                , fabs( jet1.DeltaPhi(jet2) ) ) ;
     }

     // Muons
     fillVariableWithValue(   "nMuon"                         , nMuon ) ;

     // DeltaR and DPhi(boson,jets)
     if ( nEle >= 1 && nJet >= 1) {
       fillVariableWithValue( "DR_Ele1Jet1"                   , DR_Ele1Jet1 ) ;
       TVector2 jet1 , ele1 , met, boson;
       jet1.SetMagPhi( Jet1_Pt , Jet1_Phi);
       ele1.SetMagPhi( Ele1_Pt , Ele1_Phi);
       met.SetMagPhi( MET_Pt , MET_Phi);
       boson = ele1 + met;
       fillVariableWithValue( "mDphi_BosonJet1"               , fabs( boson.DeltaPhi(jet1) ) ) ;

       TLorentzVector jet1_lorentz;
       jet1_lorentz.SetPtEtaPhiE(Jet1_Pt, Jet1_Eta, Jet1_Phi, Jet1_Energy);
       fillVariableWithValue( "DR_BosonJet1"                  , Wboson.DeltaR(jet1_lorentz) ) ;

       if(nJet >= 2) {
	 TVector2 jet2;
	 jet2.SetMagPhi( Jet2_Pt , Jet2_Phi);
	 fillVariableWithValue( "DR_Ele1Jet2"                 , DR_Ele1Jet2 ) ;
	 fillVariableWithValue( "mDphi_BosonJet2"             , fabs( boson.DeltaPhi(jet2) ) ) ;
	 fillVariableWithValue( "DR_Jet1Jet2"                 , DR_Jet1Jet2 ) ;

	 TLorentzVector jet2_lorentz;
	 jet2_lorentz.SetPtEtaPhiE(Jet2_Pt, Jet2_Eta, Jet2_Phi, Jet2_Energy);
	 fillVariableWithValue( "DR_BosonJet2"                , Wboson.DeltaR(jet2_lorentz) ) ;
       }
     }

     // DeltaPhi
     if ( nEle >= 1 ) {      
       fillVariableWithValue( "mDPhi_METEle1"                 , mDPhi_METEle1 ) ;
     }
     if ( nJet >= 1 ) { 
       fillVariableWithValue( "mDPhi_METJet1"                 , mDPhi_METJet1 ) ;
     }
     if ( nJet >= 2 ) { 
       fillVariableWithValue( "mDPhi_METJet2"                 , mDPhi_METJet2 ) ;
     }

     // sT
     if ( nEle >= 1 && nJet >= 2) {
       fillVariableWithValue( "sT_enujj"                      , sT_enujj ) ;
     }      
 
     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------
     
     bool passed_preselection = passedAllPreviousCuts("mDphi_BosonJet1");
     bool passed_preselection_without_DR_DPhi_cuts = passedAllPreviousCuts("nMuon");

     if ( passed_preselection ) { 

       FillUserTH1D( "Ele1_Pt"	           , 	Ele1_Pt       , event_weight);
       FillUserTH1D( "Ele1_Eta"	           , 	Ele1_Eta      , event_weight);
       FillUserTH1D( "Ele1_Phi"	           , 	Ele1_Phi      , event_weight);
       FillUserTH1D( "Ele1_Charge"	   , 	Ele1_Charge   , event_weight);
       
       FillUserTH1D( "MET_Pt"              ,    MET_Pt        , event_weight);
       FillUserTH1D( "MET_Phi"		   , 	MET_Phi       , event_weight);
       
       FillUserTH1D( "Jet1_Pt"             ,    Jet1_Pt       , event_weight);
       FillUserTH1D( "Jet1_Eta"	           , 	Jet1_Eta      , event_weight);
       FillUserTH1D( "Jet1_Phi"	           , 	Jet1_Phi      , event_weight);
       FillUserTH1D( "Jet1_btagTCHE"       ,    Jet1_btagTCHE , event_weight);
       FillUserTH1D( "Jet2_Pt"             ,    Jet2_Pt       , event_weight);
       FillUserTH1D( "Jet2_Eta"	           , 	Jet2_Eta      , event_weight);
       FillUserTH1D( "Jet2_Phi"	           , 	Jet2_Phi      , event_weight);
       FillUserTH1D( "Jet2_btagTCHE"       ,    Jet2_btagTCHE , event_weight);
       
       FillUserTH1D( "nEle"                ,    nEle          , event_weight);
       FillUserTH1D( "nJet"                ,    nJet          , event_weight);
       FillUserTH1D( "nJet_btagTCHE"       ,    nJet_btagTCHE , event_weight);
       FillUserTH1D( "nVertex"             ,    nVertex       , event_weight);
       FillUserTH1D( "nVertex_good"        ,    nVertex_good  , event_weight);
              
       FillUserTH1D( "MT_Ele1MET"	   , 	MT_Ele1MET    , event_weight);
       FillUserTH1D( "Pt_Ele1MET"	   , 	Pt_Ele1MET    , event_weight);
       FillUserTH1D( "Pt_j1j2"	           , 	Pt_j1j2       , event_weight);
       FillUserTH1D( "M_j1j2"	           , 	M_j1j2        , event_weight);
       FillUserTH1D( "sT_enujj"	           ,    sT_enujj      , event_weight);

       FillUserTH1D( "Jet1_Pt_over_Mj1j2"  , 	getVariableValue("Jet1_Pt_over_Mj1j2")   , event_weight);
       FillUserTH1D( "Jet2_Pt_over_Mj1j2"  , 	getVariableValue("Jet2_Pt_over_Mj1j2")   , event_weight);
       FillUserTH1D( "Jet1_E_over_Mj1j2"   , 	getVariableValue("Jet1_E_over_Mj1j2")   , event_weight);
       FillUserTH1D( "Jet2_E_over_Mj1j2"   , 	getVariableValue("Jet2_E_over_Mj1j2")   , event_weight);
       FillUserTH1D( "Jet1_E_over_Jet1_Pt" , 	getVariableValue("Jet1_E_over_Jet1_Pt")   , event_weight);
       FillUserTH1D( "Jet2_E_over_Jet2_Pt" , 	getVariableValue("Jet2_E_over_Jet2_Pt")   , event_weight);

       FillUserTH1D( "mDphi_BosonJet1"     , 	getVariableValue("mDphi_BosonJet1")   , event_weight);
       FillUserTH1D( "mDphi_BosonJet2"     , 	getVariableValue("mDphi_BosonJet2")   , event_weight);
       FillUserTH1D( "DR_BosonJet1"	   , 	getVariableValue("DR_BosonJet1")   , event_weight);
       FillUserTH1D( "DR_BosonJet2"	   , 	getVariableValue("DR_BosonJet2")   , event_weight);
       FillUserTH2D( "DR_BosonJet2_vs_DR_BosonJet1" , getVariableValue("DR_BosonJet1") , getVariableValue("DR_BosonJet2") , event_weight ); 	          

       FillUserTH1D( "nMuon"               ,    nMuon         , event_weight);

       FillUserTH1D( "DR_Ele1Jet1"	   , 	DR_Ele1Jet1   , event_weight);
       FillUserTH1D( "DR_Ele1Jet2"	   , 	DR_Ele1Jet2   , event_weight);
       FillUserTH2D( "DR_Ele1Jet2_vs_DR_Ele1Jet1" , DR_Ele1Jet1 , DR_Ele1Jet2 , event_weight ); 	   
       FillUserTH1D( "DR_Jet1Jet2"	   , 	DR_Jet1Jet2   , event_weight);
       FillUserTH1D( "mDPhi_METEle1"	   , 	mDPhi_METEle1 , event_weight);
       FillUserTH1D( "mDPhi_METJet1"	   , 	mDPhi_METJet1 , event_weight);
       FillUserTH1D( "mDPhi_METJet2"	   , 	mDPhi_METJet2 , event_weight);
       FillUserTH1D( "mDEta_Jet1Jet2"	   , 	getVariableValue("mDEta_Jet1Jet2")   , event_weight);
       FillUserTH1D( "mDPhi_Jet1Jet2"	   , 	getVariableValue("mDPhi_Jet1Jet2")   , event_weight);
     }

//      if ( passed_preselection_without_DR_DPhi_cuts ) { 
//      }

   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
