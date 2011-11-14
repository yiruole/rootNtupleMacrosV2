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
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile){}

analysisClass::~analysisClass(){}

void analysisClass::Loop()
{
   std::cout << "analysisClass::Loop() begins" <<std::endl;   

   //--------------------------------------------------------------------------
   // Decide which plots to save (default is to save everything)
   //--------------------------------------------------------------------------
   
   fillSkim                         ( !true  ) ;
   fillAllPreviousCuts              ( !true  ) ;
   fillAllOtherCuts                 (  true  ) ;
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------
   
   CreateUserTH1D( "nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "nMuon_PAS"             ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "nJet_PAS"              ,    10  , -0.5    , 9.5      );
   CreateUserTH1D( "Pt1stEle_PAS"	   , 	100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta1stEle_PAS"	   , 	100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stEle_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Pt2ndEle_PAS"	   , 	100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta2ndEle_PAS"	   , 	100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi2ndEle_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Charge1stEle_PAS"	   , 	2   , -1.0001 , 1.0001	 ); 
   CreateUserTH1D( "Charge2ndEle_PAS"	   , 	2   , -1.0001 , 1.0001	 ); 
   CreateUserTH1D( "MET_PAS"               ,    200 , 0       , 1000	 ); 
   CreateUserTH1D( "METPhi_PAS"		   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "METCharged_PAS"        ,    200 , 0       , 1000	 ); 
   CreateUserTH1D( "METChargedPhi_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "METType1_PAS"          ,    200 , 0       , 1000	 ); 
   CreateUserTH1D( "METType1Phi_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Pt1stJet_PAS"          ,    100 , 0       , 1000	 ); 
   CreateUserTH1D( "Pt2ndJet_PAS"          ,    100 , 0       , 1000	 ); 
   CreateUserTH1D( "Eta1stJet_PAS"         ,    100 , -5      , 5	 ); 
   CreateUserTH1D( "Eta2ndJet_PAS"         ,    100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Phi2ndJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "sTlep_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTjet_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PAS"                ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mjj_PAS"		   ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mee_PAS"		   ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "MTenu_PAS"                , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "Me1j1_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me1j2_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me2j1_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me2j2_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me1j_selected_PAS"     ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me2j_selected_PAS"     ,    200 , 0       , 2000     );
   CreateUserTH1D( "Meejj_PAS"             ,    200 , 0       , 2000     );
		                           
   CreateUserTH1D( "Ptee_PAS"              ,    200 , 0       , 2000     );
   		                           
   CreateUserTH1D( "DCotTheta1stEle_PAS"   ,    100 , 0.0, 1.0);
   CreateUserTH1D( "Dist1stEle_PAS"        ,    100 , 0.0, 1.0);  
   CreateUserTH1D( "DCotTheta2ndEle_PAS"   ,    100 , 0.0, 1.0);
   CreateUserTH1D( "Dist2ndEle_PAS"        ,    100 , 0.0, 1.0);  
		                           
   CreateUserTH1D( "nVertex_PAS"           ,    31   , -0.5   , 30.5	 ) ; 
   CreateUserTH1D( "nVertex_good_PAS"      ,    31   , -0.5   , 30.5	 ) ; 
		                           
   CreateUserTH1D( "DR_Ele1Jet1_PAS"	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele1Jet2_PAS"	   , 	getHistoNBins("DR_Ele1Jet2"), getHistoMin("DR_Ele1Jet2"), getHistoMax("DR_Ele1Jet2")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet1_PAS"	   , 	getHistoNBins("DR_Ele2Jet1"), getHistoMin("DR_Ele2Jet1"), getHistoMax("DR_Ele2Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet2_PAS"	   , 	getHistoNBins("DR_Ele2Jet2"), getHistoMin("DR_Ele2Jet2"), getHistoMax("DR_Ele2Jet2")     ) ; 
   CreateUserTH1D( "DR_Jet1Jet2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 

   CreateUserTH2D( "Me1jVsMe2j_selected",     200, 0, 2000, 200, 0, 2000) ;
   CreateUserTH2D( "Me1jVsMe2j_rejected",     200, 0, 2000, 200, 0, 2000) ;

   CreateUserTH1D( "MTeemunu_PAS"          ,    200 , 0       , 1000	 ); 
   
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
     // Do pileup re-weighting
     //--------------------------------------------------------------------------
     
     int NPILEUP_AVE = int( nPileUpInt_BX0 );
     int NPILEUP_FINAL = min( NPILEUP_AVE , 25 );
     double weight = getPileupWeight ( NPILEUP_FINAL, isData ) ;

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "PassJSON"                      , passedJSON ) ; 

     // Noise filters
     fillVariableWithValue(   "PassHBHENoiseFilter"           , PassHBHENoiseFilter ) ; 
     fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight ) ; 

     // Fill HLT 
     int passHLT = 1;
     if ( isData ) { 
       passHLT = 0;
       // if ( H_17_8_CIdT_CIsVL_TIdVL_TIsVL_5 == 1 || 
       // 	    H_17_8_CIdT_CIsVL_TIdVL_TIsVL_6 == 1 || 
       // 	    H_17_8_CIdT_CIsVL_TIdVL_TIsVL_7 == 1 || 
       // 	    H_17_8_CIdT_CIsVL_TIdVL_TIsVL_8 == 1 || 
       // 	    H_17_8_CIdT_TIdVL_CIsVL_TIsVL_1 == 1 || 
       // 	    H_17_8_CIdT_TIdVL_CIsVL_TIsVL_2 == 1 || 
       // 	    H_17_8_CIdT_TIdVL_CIsVL_TIsVL_3 == 1 || 
       // 	    H_17_8_CIdT_TIdVL_CIsVL_TIsVL_4 == 1 || 
       // 	    H_17_8_CIdT_TIdVL_CIsVL_TIsVL_5 == 1    ) {
       // 	 passHLT = 1;
       // }
       
       if ( H_DoubleEle33_CIdL_1  == 1 || 
       	    H_DoubleEle33_CIdL_2  == 1 || 
       	    H_DoubleEle33_CIdL_3  == 1 || 
       	    H_DoubleEle33_CIdL_4  == 1 || 
       	    H_DoubleEle33_CIdL_5  == 1 || 
       	    H_DoublePhoton33_1    == 1 || 
       	    H_DoublePhoton33_2    == 1 || 
       	    H_DoublePhoton33_3    == 1 || 
       	    H_DoublePhoton33_4    == 1 || 
       	    H_DoublePhoton33_5    == 1  ) {
       	 passHLT = 1;
       }
     }
     
     fillVariableWithValue ( "PassHLT", passHLT ) ;
     
     // Electrons
     fillVariableWithValue(   "nEle"                          , nEle_Ana ) ;
     if ( nEle_Ana >= 1 ) { 
       fillVariableWithValue( "Ele1_Pt"                       , Ele1_Pt ) ;
       fillVariableWithValue( "Ele1_Eta"                      , Ele1_Eta ) ;
     }
     if ( nEle_Ana >= 2 ) { 
       fillVariableWithValue( "Ele2_Pt"                       , Ele2_Pt ) ;
       fillVariableWithValue( "Ele2_Eta"                      , Ele2_Eta ) ;
       fillVariableWithValue( "M_e1e2"                        , M_e1e2 ) ;
       fillVariableWithValue( "Pt_e1e2"                       , Pt_e1e2 ) ;
     }

     // Jets
     fillVariableWithValue(   "nJet"                          , nJet_Ana ) ;
     if ( nJet_Ana >= 1 ) { 
       fillVariableWithValue( "Jet1_Pt"                       , Jet1_Pt ) ;
       fillVariableWithValue( "Jet1_Eta"                      , Jet1_Eta ) ;
     }
     if ( nJet_Ana >= 2 ) { 
       fillVariableWithValue( "Jet2_Pt"                       , Jet2_Pt ) ;
       fillVariableWithValue( "Jet2_Eta"                      , Jet2_Eta ) ;
       fillVariableWithValue( "Pt_j1j2"                       , Pt_j1j2 ) ;
       fillVariableWithValue( "M_j1j2"                        , M_j1j2 ) ;
       fillVariableWithValue( "DR_Jet1Jet2"                   , DR_Jet1Jet2 ) ;
     }

     // Muons
     fillVariableWithValue(   "nMuon"                         , nMuon_Ana ) ;

     // DeltaR
     if ( nEle_Ana >= 2 && nJet_Ana >= 1) {
       fillVariableWithValue( "DR_Ele1Jet1"                   , DR_Ele1Jet1 ) ;
       fillVariableWithValue( "DR_Ele2Jet1"                   , DR_Ele2Jet1 ) ;
       if(nJet_Ana >= 2) {
	 fillVariableWithValue( "DR_Ele1Jet2"                 , DR_Ele1Jet2 ) ;
	 fillVariableWithValue( "DR_Ele2Jet2"                 , DR_Ele2Jet2 ) ;
       }
     }

     // sT
     if ( nEle_Ana >= 2 && nJet_Ana >= 2) {
       fillVariableWithValue( "sT_eejj"                      , sT_eejj ) ;
     }      

     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------
     

     bool passed_preselection = ( passedAllPreviousCuts("M_e1e2") && passedCut ("M_e1e2") );
     
     if ( passed_preselection ) {

       FillUserTH1D("nElectron_PAS"        , nEle_Ana         , weight) ;
       FillUserTH1D("nMuon_PAS"            , nMuon_Stored     , weight) ;
       FillUserTH1D("nJet_PAS"             , nJet_Ana         , weight) ;
       FillUserTH1D("Pt1stEle_PAS"	   , Ele1_Pt          , weight) ;
       FillUserTH1D("Eta1stEle_PAS"	   , Ele1_Eta         , weight) ;
       FillUserTH1D("Phi1stEle_PAS"	   , Ele1_Phi         , weight) ;
       FillUserTH1D("Pt2ndEle_PAS"	   , Ele2_Pt          , weight) ;
       FillUserTH1D("Eta2ndEle_PAS"	   , Ele2_Eta         , weight) ;
       FillUserTH1D("Phi2ndEle_PAS"	   , Ele2_Phi         , weight) ;
       FillUserTH1D("Charge1stEle_PAS"	   , Ele1_Charge      , weight) ;
       FillUserTH1D("Charge2ndEle_PAS"	   , Ele2_Charge      , weight) ;
       FillUserTH1D("MET_PAS"              , MET_Pt           , weight) ;
       FillUserTH1D("METPhi_PAS"	   , MET_Phi          , weight) ;
       FillUserTH1D( "METCharged_PAS"      , PFMETCharged     , weight);
       FillUserTH1D( "METChargedPhi_PAS"   , PFMETChargedPhi  , weight);   
       FillUserTH1D( "METType1_PAS"        , PFMETType1Cor    , weight);
       FillUserTH1D( "METType1Phi_PAS"     , PFMETPhiType1Cor , weight);   
       FillUserTH1D("Pt1stJet_PAS"         , Jet1_Pt          , weight) ;
       FillUserTH1D("Pt2ndJet_PAS"         , Jet2_Pt          , weight) ;
       FillUserTH1D("Eta1stJet_PAS"        , Jet1_Eta         , weight) ;
       FillUserTH1D("Eta2ndJet_PAS"        , Jet2_Eta         , weight) ;
       FillUserTH1D("Phi1stJet_PAS"	   , Jet1_Phi         , weight) ;
       FillUserTH1D("Phi2ndJet_PAS"	   , Jet2_Phi         , weight) ;
       FillUserTH1D("sTlep_PAS"            , Ele1_Pt + Ele2_Pt, weight) ;
       FillUserTH1D("sTjet_PAS"            , Jet1_Pt + Jet2_Pt, weight) ;
       FillUserTH1D("sT_PAS"               , sT_eejj          , weight) ;
       FillUserTH1D("Mjj_PAS"		   , M_j1j2           , weight) ;
       FillUserTH1D("Mee_PAS"		   , M_e1e2           , weight) ;
       FillUserTH1D( "MTenu_PAS"           , MT_Ele1MET       , weight);
       FillUserTH1D("Me1j1_PAS"            , M_e1j1           , weight) ;
       FillUserTH1D("Me1j2_PAS"            , M_e1j2           , weight) ;
       FillUserTH1D("Me2j1_PAS"            , M_e2j1           , weight) ;
       FillUserTH1D("Me2j2_PAS"            , M_e2j2           , weight) ;
       FillUserTH1D("Ptee_PAS"             , Pt_e1e2          , weight) ;
       FillUserTH1D("DCotTheta1stEle_PAS"  , Ele1_DCotTheta   , weight) ;
       FillUserTH1D("Dist1stEle_PAS"       , Ele1_Dist        , weight) ;
       FillUserTH1D("DCotTheta2ndEle_PAS"  , Ele2_DCotTheta   , weight) ;
       FillUserTH1D("Dist2ndEle_PAS"       , Ele2_Dist        , weight) ;
       FillUserTH1D("nVertex_PAS"          , nVertex          , weight) ;
       FillUserTH1D("nVertex_good_PAS"     , nVertex_good     , weight) ;
       FillUserTH1D("DR_Ele1Jet1_PAS"	   , DR_Ele1Jet1      , weight) ;
       FillUserTH1D("DR_Ele1Jet2_PAS"	   , DR_Ele1Jet2      , weight) ;
       FillUserTH1D("DR_Ele2Jet1_PAS"	   , DR_Ele2Jet1      , weight) ;
       FillUserTH1D("DR_Ele2Jet2_PAS"	   , DR_Ele2Jet2      , weight) ;
       FillUserTH1D("DR_Jet1Jet2_PAS"	   , DR_Jet1Jet2      , weight) ;
       
       TLorentzVector eejj, e1e2mu, e1, j1, e2, j2, mu, met;
       e1.SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
       e2.SetPtEtaPhiM ( Ele2_Pt, Ele2_Eta, Ele2_Phi, 0.0 );
       j1.SetPtEtaPhiM ( Jet1_Pt, Jet1_Eta, Jet1_Phi, 0.0 );
       j2.SetPtEtaPhiM ( Jet2_Pt, Jet2_Eta, Jet2_Phi, 0.0 );
       mu.SetPtEtaPhiM ( Muon1_Pt, Muon1_Eta, Muon1_Phi, 0.0 );
       met.SetPtEtaPhiM ( MET_Pt, 0.0, MET_Phi, 0.0 );
       
       eejj = e1 + e2 + j1 + j2 ; 
       double M_eejj = eejj.M();

       // e1e2mu = e1 + e2 + mu;
       // double MT_eemuMET = sqrt(2 * e1e2mu.Pt()    * MET_Pt  * (1 - cos(e1e2mu.DeltaPhi (met))));
       // FillUserTH1D("MTeemunu_PAS"         , MT_eemuMET                        , weight );
	      
       FillUserTH1D("Meejj_PAS", M_eejj , weight );

       if ( fabs(M_e1j1-M_e2j2) < fabs(M_e1j2-M_e2j1) )  {
	 FillUserTH1D("Me1j_selected_PAS" , M_e1j1, weight) ;	   
	 FillUserTH1D("Me2j_selected_PAS" , M_e2j2, weight) ;	   
	 FillUserTH2D( "Me1jVsMe2j_selected", M_e1j1, M_e2j2, weight ) ;
	 FillUserTH2D( "Me1jVsMe2j_rejected", M_e1j2, M_e2j1, weight ) ;
       }
       else {
	 FillUserTH1D("Me1j_selected_PAS" , M_e1j2, weight) ;	   
	 FillUserTH1D("Me2j_selected_PAS" , M_e2j1, weight) ;	   
	 FillUserTH2D( "Me1jVsMe2j_selected", M_e1j2, M_e2j1, weight ) ;
	 FillUserTH2D( "Me1jVsMe2j_rejected", M_e1j1, M_e2j2, weight ) ;
       }
     }
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

