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
   
   fillSkim                         (  true  ) ;
   fillAllPreviousCuts              (  true  ) ;
   fillAllOtherCuts                 ( !true  ) ;
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------

   CreateUserTH1D( "Pt1stEle_PAS"	   , 	100 , 0       , 1000     ) ; 
   CreateUserTH1D( "Eta1stEle_PAS"	   , 	100 , -5      , 5	 ) ; 
   CreateUserTH1D( "Phi1stEle_PAS"	   , 	60  , -3.1416 , +3.1416	 ) ; 
   CreateUserTH1D( "Charge1stEle_PAS"	   , 	2   , -1.0001 , 1.0001	 ) ; 
   CreateUserTH1D( "MET_PAS"               ,    200 , 0       , 1000	 ) ; 
   CreateUserTH1D( "METPhi_PAS"		   , 	60  , -3.1416 , +3.1416	 ) ; 
   CreateUserTH1D( "minMETPt1stEle_PAS"    ,    200 , 0       , 1000	 ) ; 
   CreateUserTH1D( "Pt1stJet_PAS"          ,    100 , 0       , 1000	 ) ; 
   CreateUserTH1D( "Pt2ndJet_PAS"          ,    100 , 0       , 1000	 ) ; 
   CreateUserTH1D( "Eta1stJet_PAS"         ,    100 , -5      , 5	 ) ; 
   CreateUserTH1D( "Eta2ndJet_PAS"         ,    100 , -5      , 5	 ) ; 
   CreateUserTH1D( "Phi1stJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ) ; 
   CreateUserTH1D( "Phi2ndJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ) ; 
   CreateUserTH1D( "TCHE1stJet_PAS"        ,    100 , 0       , 20	 ) ; 
   CreateUserTH1D( "TCHE2ndJet_PAS"        ,    100 , 0       , 20	 ) ; 
   CreateUserTH1D( "nMuon_PtCut_IDISO_PAS" ,    16  , -0.5    , 15.5	 ) ; 
   CreateUserTH1D( "MTenu_PAS"             ,    200 , 0       , 1000	 ) ; 
   CreateUserTH1D( "Ptenu_PAS"		   , 	200 , 0       , 2000	 ) ; 
   CreateUserTH1D( "sTlep_PAS"             ,    200 , 0       , 2000	 ) ; 
   CreateUserTH1D( "sTjet_PAS"             ,    200 , 0       , 2000	 ) ; 
   CreateUserTH1D( "sT_PAS"                ,    200 , 0       , 2000	 ) ; 
   CreateUserTH1D( "Mjj_PAS"		   ,    200 , 0       , 2000	 ) ; 
   CreateUserTH1D( "Mej_1stPair_PAS"       ,    200 , 0       , 2000	 ) ; 
   CreateUserTH1D( "Mej_2ndPair_PAS"       ,    200 , 0       , 2000	 ) ; 
   CreateUserTH1D( "HcalIso1stEle_PAS"     ,    200 , 0       , 20       ) ;
   CreateUserTH1D( "EcalIso1stEle_PAS"     ,    200 , 0       , 20       ) ;
   CreateUserTH1D( "RelIso1stEle_PAS"      ,    200 , 0       , 1.0      ) ;
   CreateUserTH1D( "Mee_PAS"               ,    200 , 0       , 2000     ) ;
   CreateUserTH1D( "Ptee_PAS"              ,    200 , 0       , 2000     ) ;
   CreateUserTH1D( "DCotTheta1stEle_PAS"   ,    100 , 0.0, 1.0);
   CreateUserTH1D( "Dist1stEle_PAS"        ,    100 , 0.0, 1.0);  
   CreateUserTH1D( "mDPhi1stEleMET", 100, 0.,  3.14159 ) ;
   CreateUserTH1D( "mDPhi1stJetMET", 100, 0.,  3.14159 ) ;
   CreateUserTH1D( "mDPhi2ndJetMET", 100, 0.,  3.14159 ) ;
   CreateUserTH1D( "MT_GoodVtxLTE5", 200, 0.,  1000 ) ;
   CreateUserTH1D( "MT_GoodVtxGT5" , 200, 0.,  1000 ) ;

   
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
     // check if pass trigger
     //--------------------------------------------------------------------------

     bool passTrigger = true;
     if ( isData ) passTrigger = Decision_HLT_Ele22_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_CentralJet30_CentralJet25_PFMHT20_v2;

     //--------------------------------------------------------------------------
     // Trigger emulating variables
     //--------------------------------------------------------------------------

     fillVariableWithValue("PassTrigger"                     , (int) passTrigger ) ;
     fillVariableWithValue("Pt1stCentralCaloJet"             , CaloJet1_Pt  );
     fillVariableWithValue("Pt2ndCentralCaloJet"             , CaloJet2_Pt  );
     fillVariableWithValue("MET"                             , MET_Pt       );
     fillVariableWithValue("nEle_PtCut_IDISO_noOvrlp"        , nEle         ); 
     fillVariableWithValue("Pt1stEle_PtCut_IDISO_noOvrlp"    , Ele1_Pt      );
     
     //--------------------------------------------------------------------------
     // Do pileup re-weighting
     //--------------------------------------------------------------------------
     
     double weight = getPileupWeight ( nPileUpInteractions, isData ) ;

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "PassJSON"                      , passedJSON   ) ; 

     // Vertex variables
     fillVariableWithValue(   "nVertex_presel"                , nVertex      ) ;
     fillVariableWithValue(   "nVertex_good_presel"           , nVertex_good ) ;

     // Muon variables ( for veto ) 
     fillVariableWithValue(   "nMuon_PtCut_ID_ISO"            , nMuon        ) ;
			                                      
     // 1st Electron variables				      
     fillVariableWithValue(   "Pt1stEle_PtCut_IDISO_noOvrlp"  , Ele1_Pt      ) ;
     fillVariableWithValue(   "Eta1stEle_PtCut_IDISO_noOvrlp" , Ele1_Eta     ) ;
     fillVariableWithValue(   "Phi1stEle_PtCut_IDISO_noOvrlp" , Ele1_Phi     ) ;

     // MET variables	                                      
     fillVariableWithValue(   "MET_presel"                    , MET_Pt       ) ;
     fillVariableWithValue(   "METPhi"                        , MET_Phi      ) ;
     fillVariableWithValue(   "MT_presel"                     , MT_Ele1MET   ) ;
     
     // 1st JET variables                                     
     fillVariableWithValue(   "nJet_PtCut_ID_noOvrlp"         , nJet         );
     fillVariableWithValue(   "Pt1stJet_PtCut_ID_noOvrlp"     , Jet1_Pt      );
     fillVariableWithValue(   "Eta1stJet_PtCut_ID_noOvrlp"    , Jet1_Eta     );
     fillVariableWithValue(   "Phi1stJet_PtCut_ID_noOvrlp"    , Jet1_Phi     );
			                                      
     // 2nd JET variables                                     
     if ( nJet > 1 ) { 	                                      
       fillVariableWithValue( "Pt2ndJet_PtCut_ID_noOvrlp"     , Jet2_Pt      );
       fillVariableWithValue( "Eta2ndJet_PtCut_ID_noOvrlp"    , Jet2_Eta     );
       fillVariableWithValue( "Phi2ndJet_PtCut_ID_noOvrlp"    , Jet2_Phi     );

       fillVariableWithValue( "ST"                , sT_enujj      ) ;
       fillVariableWithValue( "ST_presel"         , sT_enujj      ) ;
       fillVariableWithValue( "mDeltaPhiMETEle"   , mDPhi_METEle1 ) ;
       fillVariableWithValue( "mDeltaPhiMET1stJet", mDPhi_METJet1 ) ;
     }

     if ( nEle > 2 ) { 
       FillUserTH1D( "Mee_PAS"   , M_e1e2  ) ;
       FillUserTH1D( "Ptee_PAS"  , Pt_e1e2 ) ;
     }

     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------
     
     bool passed_preselection = passedAllPreviousCuts("nVertex_presel");

     if ( passed_preselection ) { 
       
       FillUserTH1D( "Pt1stEle_PAS"	     , Ele1_Pt     ) ; 
       FillUserTH1D( "Eta1stEle_PAS"	     , Ele1_Eta    ) ;
       FillUserTH1D( "Phi1stEle_PAS"	     , Ele1_Phi    ) ;
       FillUserTH1D( "Charge1stEle_PAS"      , Ele1_Charge ) ;   
       FillUserTH1D( "MET_PAS"               , MET_Pt        ) ;
       FillUserTH1D( "METPhi_PAS"	     , MET_Phi       ) ;   
       FillUserTH1D( "minMETPt1stEle_PAS"    , TMath::Min ( Ele1_Pt, MET_Pt  ));
       FillUserTH1D( "Pt1stJet_PAS"          , Jet1_Pt           );
       FillUserTH1D( "Pt2ndJet_PAS"          , Jet2_Pt           );
       FillUserTH1D( "Eta1stJet_PAS"         , Jet1_Eta          );
       FillUserTH1D( "Eta2ndJet_PAS"         , Jet2_Eta          );
       FillUserTH1D( "Phi1stJet_PAS"         , Jet1_Phi          );
       FillUserTH1D( "Phi2ndJet_PAS"	     , Jet2_Phi          );
       FillUserTH1D( "TCHE1stJet_PAS"        , Jet1_btagTCHE     );
       FillUserTH1D( "TCHE2ndJet_PAS"        , Jet2_btagTCHE     );
       FillUserTH1D( "nMuon_PtCut_IDISO_PAS" , nMuon             ); 
       FillUserTH1D( "MTenu_PAS"             , MT_Ele1MET        );
       FillUserTH1D( "Ptenu_PAS"	     , Pt_Ele1MET        );
       FillUserTH1D( "sTlep_PAS"             , Ele1_Pt + MET_Pt  );
       FillUserTH1D( "sTjet_PAS"             , Jet1_Pt + Jet2_Pt );
       FillUserTH1D( "sT_PAS"                , sT_enujj      ) ;
       FillUserTH1D( "Mjj_PAS"	             , M_j1j2        ) ;   
       FillUserTH1D( "Mej_1stPair_PAS"       , M_ej1         ) ;
       FillUserTH1D( "Mej_2ndPair_PAS"       , M_ej2         ) ;
       FillUserTH1D( "HcalIso1stEle_PAS"   , Ele1_HcalIso ) ;
       FillUserTH1D( "EcalIso1stEle_PAS"   , Ele1_EcalIso ) ;
       FillUserTH1D( "RelIso1stEle_PAS"    , Ele1_RelIso ) ;
       FillUserTH1D( "DCotTheta1stEle_PAS" , Ele1_DCotTheta ) ;
       FillUserTH1D( "Dist1stEle_PAS"      , Ele1_Dist ) ;

       FillUserTH1D( "mDPhi1stEleMET", mDPhi_METEle1 ) ;
       FillUserTH1D( "mDPhi1stJetMET", mDPhi_METJet1 ) ;
       FillUserTH1D( "mDPhi2ndJetMET", mDPhi_METJet2 ) ;
       
       if ( nVertex_good <= 5 ) FillUserTH1D ( "MT_GoodVtxLTE5", MT_Ele1MET );
       if ( nVertex_good >  5 ) FillUserTH1D ( "MT_GoodVtxGT5" , MT_Ele1MET );
       
       // FillUserTH1D( "Vtxd01stEle_PAS"       , ) ;
       // FillUserTH1D( "MissingHits1stEle_PAS" , ) ;
       // FillUserTH1D( "Dist1stEle_PAS"        , ) ;
       // FillUserTH1D( "DCotTheta1stEle_PAS"   , ) ;
       // FillUserTH1D( "MTnuj_1stPair_PAS"     , ) ;
       // FillUserTH1D( "MTnuj_2ndPair_PAS"     , ) ;
       // FillUserTH1D( "nVertex_PAS"           , ) ;
       // FillUserTH1D( "nVertex_good_PAS"      , ) ;
     }
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
