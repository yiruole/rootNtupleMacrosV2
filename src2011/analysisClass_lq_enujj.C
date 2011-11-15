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
   // Set global variables
   //--------------------------------------------------------------------------
   
   TVector2 v_METCharged, v_METType1, v_ele;

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------

   CreateUserTH1D( "nElectron_PAS"            , 5   , -0.5    , 4.5      );
   CreateUserTH1D( "nMuon_PAS"                , 5   , -0.5    , 4.5      );
   CreateUserTH1D( "Pt1stEle_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta1stEle_PAS"	      , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stEle_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Charge1stEle_PAS"	      , 2   , -1.0001 , 1.0001	 ); 
   CreateUserTH1D( "MET_PAS"                  , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "METPhi_PAS"		      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "METCharged_PAS"           , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "METChargedPhi_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "METType1_PAS"             , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "METType1Phi_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "METSig_PAS"               , 100 , 0       , 200      );
   CreateUserTH1D( "minMETPt1stEle_PAS"       , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "Pt1stJet_PAS"             , 100 , 0       , 1000	 ); 
   CreateUserTH1D( "Pt2ndJet_PAS"             , 100 , 0       , 1000	 ); 
   CreateUserTH1D( "Eta1stJet_PAS"            , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Eta2ndJet_PAS"            , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stJet_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Phi2ndJet_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "TCHE1stJet_PAS"           , 100 , 0       , 20	 ); 
   CreateUserTH1D( "TCHE2ndJet_PAS"           , 100 , 0       , 20	 ); 
   CreateUserTH1D( "nMuon_PtCut_IDISO_PAS"    , 16  , -0.5    , 15.5	 ); 
   CreateUserTH1D( "MTenu_PAS"                , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "MT_charged_enu_PAS"       , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "MT_type1_enu_PAS"         , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "Ptenu_PAS"		      , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTlep_PAS"                , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTjet_PAS"                , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PAS"                   , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mjj_PAS"		      , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mej1_PAS"                 , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mej2_PAS"                 , 200 , 0       , 2000	 ); 
   CreateUserTH1D( "DCotTheta1stEle_PAS"      , 100 , 0.0     , 1.0      );
   CreateUserTH1D( "Dist1stEle_PAS"           , 100 , 0.0     , 1.0      );  
   CreateUserTH1D( "DR_Ele1Jet1_PAS"	      , 100 , 0       , 10       ); 
   CreateUserTH1D( "DR_Ele1Jet2_PAS"	      , 100 , 0       , 10       ); 
   CreateUserTH1D( "DR_Jet1Jet2_PAS"	      , 100 , 0       , 10       ); 
   CreateUserTH1D( "mDPhi1stEleMET_PAS"       , 100 , 0.      ,  3.14159 );
   CreateUserTH1D( "mDPhi1stJetMET_PAS"       , 100 , 0.      ,  3.14159 );
   CreateUserTH1D( "mDPhi2ndJetMET_PAS"       , 100 , 0.      ,  3.14159 );

   CreateUserTH1D( "MT_GoodVtxLTE3_PAS"       , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MT_GoodVtxGTE4_LTE8_PAS"  , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MT_GoodVtxGTE9_LTE15_PAS" , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MT_GoodVtxGTE16_PAS"      , 200 , 0.      ,  1000    );

   CreateUserTH1D( "MTCharged_GoodVtxLTE3_PAS"       , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTCharged_GoodVtxGTE4_LTE8_PAS"  , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTCharged_GoodVtxGTE9_LTE15_PAS" , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTCharged_GoodVtxGTE16_PAS"      , 200 , 0.      ,  1000    );

   CreateUserTH1D( "MTType1_GoodVtxLTE3_PAS"       , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTType1_GoodVtxGTE4_LTE8_PAS"  , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTType1_GoodVtxGTE9_LTE15_PAS" , 200 , 0.      ,  1000    );
   CreateUserTH1D( "MTType1_GoodVtxGTE16_PAS"      , 200 , 0.      ,  1000    );
   
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
     if(jentry < 10 || jentry%10000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << "/" << nentries << std::endl;   

     //--------------------------------------------------------------------------
     // Reset the cuts
     //--------------------------------------------------------------------------

     resetCuts();

     //--------------------------------------------------------------------------
     // Do pileup re-weighting
     //--------------------------------------------------------------------------
     
     // int NPILEUP_AVE = int( (nPileUpInt_BXminus1 + nPileUpInt_BX0 + nPileUpInt_BXplus1)/3 );
     int NPILEUP_AVE = int( nPileUpInt_BX0 );
     int NPILEUP_FINAL = min ( NPILEUP_AVE , 25 );
     double weight = getPileupWeight ( NPILEUP_FINAL, isData ) ;
     //double weight     = getPileupWeight ( nPileUpInteractions, isData ) ;

     //--------------------------------------------------------------------------
     // Check good run list
     //--------------------------------------------------------------------------
     
     int passedJSON = passJSON ( run, ls , isData ) ;
     if ( !isData ) passedJSON = 1;
     
     //--------------------------------------------------------------------------
     // Check HLT
     //--------------------------------------------------------------------------

     int passedHLT = 1;
     if ( isData ) {
       passedHLT = 0;
       if ( H_27_CIdVT_CIsT_TIdT_TIsT_1 == 1 || // 160405 - 161119, HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1
	    H_27_CIdVT_CIsT_TIdT_TIsT_2 == 1 || // 161217 - 163261, HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2
	    H_27_CIdVT_CIsT_TIdT_TIsT_3 == 1 || // 163270 - 163817, HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3
	    H_32_CIdVT_CIsT_TIdT_TIsT_3 == 1 || // 165088 - 165633, HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3
	    H_32_CIdVT_CIsT_TIdT_TIsT_4 == 1 || // 165970 - 166967, HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4
	    H_52_CIdVT_TIdT_3           == 1 || // 167039 - 167913, HLT_Ele52_CaloIdVT_TrkIdT_v3		     
	    H_52_CIdVT_TIdT_4           == 1 || // 170249 - 172952, HLT_Ele52_CaloIdVT_TrkIdT_v4		     
	    H_65_CIdVT_TIdT_3           == 1 || // 172953 - 173198, HLT_Ele65_CaloIdVT_TrkIdT_v3		     
	    H_65_CIdVT_TIdT_4           == 1 )  // 173236 - 173692, HLT_Ele65_CaloIdVT_TrkIdT_v4		     
	 passedHLT = 1; 
     }
     
     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "Reweighting"              , 1                       );
     fillVariableWithValue(   "PassJSON"                 , passedJSON              ); 
     									          
     // HLT variable							           
     fillVariableWithValue(   "PassHLT"                  , passedHLT               );
     
     // Filters
     fillVariableWithValue(   "PassHBHENoiseFilter"      , PassHBHENoiseFilter     );
     fillVariableWithValue(   "PassBeamHaloFilterTight"  , PassBeamHaloFilterTight );
									      
     // Muon variables ( for veto ) 					      
     fillVariableWithValue(   "nMuon"                    , nMuon_Ana               );
			                                      		                
     // 1st Electron variables				      		                
     fillVariableWithValue(   "nEle"                     , nEle_Ana                ); 
     fillVariableWithValue(   "Ele1_Pt"                  , Ele1_Pt                 );
     fillVariableWithValue(   "Ele1_Eta"                 , Ele1_Eta                );
									           
     // MET variables	                                      		           
     fillVariableWithValue(   "MET"                      , MET_Pt                  );
     fillVariableWithValue(   "mDeltaPhiMETEle"          , mDPhi_METEle1           );
     									           
     // 1st JET variables                                     		           
     fillVariableWithValue(   "nJet"                     , nJet_Ana                );
									           
     // 1st JET variables                                     		           
     if ( nJet_Stored > 0 ) { 						           
       fillVariableWithValue( "Jet1_Pt"                  , Jet1_Pt                 );
       fillVariableWithValue( "Jet1_Eta"                 , Jet1_Eta                );
       fillVariableWithValue( "mDeltaPhiMET1stJet"       , mDPhi_METJet1           );
     }									           
     									           
     // 2nd JET variables                                     		           
     if ( nJet_Stored > 1 ) { 	                                      	           
       fillVariableWithValue( "Jet2_Pt"                  , Jet2_Pt                 );
       fillVariableWithValue( "Jet2_Eta"                 , Jet2_Eta                );
       fillVariableWithValue( "ST"                       , sT_enujj                );
     }

     // 1 electron, 1 jet variables 
     if ( nEle_Ana > 0 && nJet_Ana > 0 ) { 
       fillVariableWithValue ( "DR_Ele1Jet1"             , DR_Ele1Jet1             ) ;
     }

     // 1 electron, 2 jet variables 
     if ( nEle_Ana > 0 && nJet_Ana > 1 ) { 
       fillVariableWithValue ( "DR_Ele1Jet2"             , DR_Ele1Jet2            ) ;
     }
     
     // Dummy variables
     fillVariableWithValue ("preselection",1);

     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     if (!isData && !passedCut ("PassJSON")){
       std::cout << "ERROR: This event did not pass the JSON file!" << std::endl;
       std::cout << "  isData = " << isData << std::endl;
       std::cout << "  passedJSON = " << passedJSON << std::endl;
     }
     
     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------
     
     bool passed_preselection = passedAllPreviousCuts("preselection");
       
     if ( passed_preselection ) { 

       bool use_charged_met = (PFMETCharged < MET_Pt);

       if ( use_charged_met ) v_METCharged.SetMagPhi(PFMETCharged , PFMETChargedPhi );
       else                   v_METCharged.SetMagPhi(MET_Pt       , MET_Phi         );
       
       v_METType1.SetMagPhi  (PFMETType1Cor, PFMETPhiType1Cor);
       v_ele.SetMagPhi       (Ele1_Pt      , Ele1_Phi        );
       
       double deltaphi_charged = v_METCharged.DeltaPhi(v_ele);
       double deltaphi_type1   = v_METType1  .DeltaPhi(v_ele);
       double MTCharged = sqrt(2 * Ele1_Pt * PFMETCharged  * (1 - cos(deltaphi_charged)) );
       double MTType1   = sqrt(2 * Ele1_Pt * PFMETType1Cor * (1 - cos(deltaphi_type1  )) );
              
       FillUserTH1D( "MT_charged_enu_PAS"    , MTCharged                      , weight);
       FillUserTH1D( "MT_type1_enu_PAS"      , MTType1                        , weight);
       FillUserTH1D( "nElectron_PAS"         , nEle_Ana                       , weight); 
       FillUserTH1D( "nMuon_PAS"             , nMuon_Ana                      , weight); 
       FillUserTH1D( "Pt1stEle_PAS"	     , Ele1_Pt                        , weight); 
       FillUserTH1D( "Eta1stEle_PAS"	     , Ele1_Eta                       , weight);
       FillUserTH1D( "Phi1stEle_PAS"	     , Ele1_Phi                       , weight);
       FillUserTH1D( "Charge1stEle_PAS"      , Ele1_Charge                    , weight);   
       FillUserTH1D( "METSig_PAS"	     , PFMETSig                       , weight);   
       FillUserTH1D( "MET_PAS"               , MET_Pt                         , weight);
       FillUserTH1D( "METPhi_PAS"	     , MET_Phi                        , weight);   
       FillUserTH1D( "METCharged_PAS"        , PFMETCharged                   , weight);
       FillUserTH1D( "METChargedPhi_PAS"     , PFMETChargedPhi                , weight);   
       FillUserTH1D( "METType1_PAS"          , PFMETType1Cor                  , weight);
       FillUserTH1D( "METType1Phi_PAS"       , PFMETPhiType1Cor               , weight);   
       FillUserTH1D( "minMETPt1stEle_PAS"    , TMath::Min ( Ele1_Pt, MET_Pt  ), weight);
       FillUserTH1D( "Pt1stJet_PAS"          , Jet1_Pt                        , weight);
       FillUserTH1D( "Pt2ndJet_PAS"          , Jet2_Pt                        , weight);
       FillUserTH1D( "Eta1stJet_PAS"         , Jet1_Eta                       , weight);
       FillUserTH1D( "Eta2ndJet_PAS"         , Jet2_Eta                       , weight);
       FillUserTH1D( "Phi1stJet_PAS"         , Jet1_Phi                       , weight);
       FillUserTH1D( "Phi2ndJet_PAS"	     , Jet2_Phi                       , weight);
       FillUserTH1D( "TCHE1stJet_PAS"        , Jet1_btagTCHE                  , weight);
       FillUserTH1D( "TCHE2ndJet_PAS"        , Jet2_btagTCHE                  , weight);
       FillUserTH1D( "nMuon_PtCut_IDISO_PAS" , nMuon_Ana                      , weight); 
       FillUserTH1D( "MTenu_PAS"             , MT_Ele1MET                     , weight);
       FillUserTH1D( "Ptenu_PAS"	     , Pt_Ele1MET                     , weight);
       FillUserTH1D( "sTlep_PAS"             , Ele1_Pt + MET_Pt               , weight);
       FillUserTH1D( "sTjet_PAS"             , Jet1_Pt + Jet2_Pt              , weight);
       FillUserTH1D( "sT_PAS"                , sT_enujj                       , weight);
       FillUserTH1D( "Mjj_PAS"	             , M_j1j2                         , weight);   
       FillUserTH1D( "DCotTheta1stEle_PAS"   , Ele1_DCotTheta                 , weight);
       FillUserTH1D( "Dist1stEle_PAS"        , Ele1_Dist                      , weight);
       FillUserTH1D( "mDPhi1stEleMET_PAS"    , mDPhi_METEle1                  , weight);
       FillUserTH1D( "mDPhi1stJetMET_PAS"    , mDPhi_METJet1                  , weight);
       FillUserTH1D( "mDPhi2ndJetMET_PAS"    , mDPhi_METJet2                  , weight); 
       FillUserTH1D( "Mej1_PAS"              , M_e1j1                         , weight);
       FillUserTH1D( "Mej2_PAS"              , M_e1j2                         , weight);
       FillUserTH1D( "DR_Ele1Jet1_PAS"	     , DR_Ele1Jet1                    , weight);
       FillUserTH1D( "DR_Ele1Jet2_PAS"	     , DR_Ele1Jet2                    , weight);
       FillUserTH1D( "DR_Jet1Jet2_PAS"	     , DR_Jet1Jet2                    , weight);

       if ( nVertex_good >= 0 && nVertex_good <= 3 ) {
	 FillUserTH1D( "MT_GoodVtxLTE3_PAS"              , MT_Ele1MET, weight ) ;
	 FillUserTH1D( "MTCharged_GoodVtxLTE3_PAS"       , MTCharged , weight ) ;
	 FillUserTH1D( "MTType1_GoodVtxLTE3_PAS"         , MTType1   , weight ) ;
       }						 
       							 
       if ( nVertex_good >= 4 && nVertex_good <= 8 ) {	 
	 FillUserTH1D( "MT_GoodVtxGTE4_LTE8_PAS"         , MT_Ele1MET, weight ) ;
	 FillUserTH1D( "MTCharged_GoodVtxGTE4_LTE8_PAS"  , MTCharged , weight ) ;
	 FillUserTH1D( "MTType1_GoodVtxGTE4_LTE8_PAS"    , MTType1   , weight ) ;
       }
       
       if ( nVertex_good >= 9 && nVertex_good <= 15) {
	 FillUserTH1D( "MT_GoodVtxGTE9_LTE15_PAS"        , MT_Ele1MET, weight ) ;
	 FillUserTH1D( "MTCharged_GoodVtxGTE9_LTE15_PAS" , MTCharged , weight ) ;
	 FillUserTH1D( "MTType1_GoodVtxGTE9_LTE15_PAS"   , MTType1   , weight ) ;
       }
       
       if ( nVertex_good >= 16                     ) {
	 FillUserTH1D( "MT_GoodVtxGTE16_PAS"             , MT_Ele1MET, weight ) ;
	 FillUserTH1D( "MTCharged_GoodVtxGTE16_PAS"      , MTCharged , weight ) ;
	 FillUserTH1D( "MTType1_GoodVtxGTE16_PAS"        , MTType1   , weight ) ;
       }
              
     }
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
