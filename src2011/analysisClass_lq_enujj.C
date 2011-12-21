#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

#include "Lumi3DReWeighting.h"

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
   // Get pre-cut values
   //--------------------------------------------------------------------------

   // dataset
   
   int dataset  = getPreCutValue1("dataset") ;
   bool select2011A = ( dataset == 0 );
   bool select2011B = ( dataset == 1 );
   bool select2011  = ( dataset == 2 );

   if ( ! select2011A &&
	! select2011B &&
	! select2011 ) {
     std::cout << "Error: Must choose dataset to be 0 (2011A), 1 (2011B), or 2 (all 2011)" << std::endl;
   }

   //--------------------------------------------------------------------------
   // Set global variables
   //--------------------------------------------------------------------------
   
   TVector2 v_METCharged, v_METType1, v_ele;
   
   //--------------------------------------------------------------------------
   // Pileup reweighting initialization
   //--------------------------------------------------------------------------

   // Lumi3DReWeighting lumiWeights = Lumi3DReWeighting("/afs/cern.ch/user/e/eberry/public/LQ_PILEUP/PUMC_dist.root", "/afs/cern.ch/user/e/eberry/public/LQ_PILEUP/PUData_dist.root", "pileup", "pileup");
   // lumiWeights.weight3D_init("/afs/cern.ch/user/e/eberry/public/LQ_PILEUP/Weight3D.root");
     
   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------

   CreateUserTH1D( "nElectron_PAS"            , 5   , -0.5    , 4.5      );
   CreateUserTH1D( "nMuon_PAS"                , 5   , -0.5    , 4.5      );
   CreateUserTH1D( "Pt1stEle_PAS"	      , 100 , 0       , 1000     ); 
   CreateUserTH1D( "Eta1stEle_PAS"	      , 100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stEle_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Charge1stEle_PAS"	      , 2   , -1.0001 , 1.0001	 ); 
   CreateUserTH1D( "MatchPhotonConv1stEle_PAS",	2   , -0.5    , 1.5      );
   CreateUserTH1D( "MatchPhotonConv2ndEle_PAS",	2   , -0.5    , 1.5      );
   CreateUserTH1D( "MET_PAS"                  , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "METPhi_PAS"		      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "METCharged_PAS"           , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "METChargedPhi_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "METType1_PAS"             , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "METType1Phi_PAS"	      , 60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "METSig_PAS"               , 100 , 0       , 800      );
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
   CreateUserTH1D( "minDR_EleJet_PAS"         , 100 , 0       , 10       );  
   CreateUserTH1D( "mDPhi1stEleMET_PAS"       , 100 , 0.      , 3.14159  );
   CreateUserTH1D( "mDPhi1stJetMET_PAS"       , 100 , 0.      , 3.14159  );
   CreateUserTH1D( "mDPhi2ndJetMET_PAS"       , 100 , 0.      , 3.14159  );
   CreateUserTH1D( "GeneratorWeight"          , 200 , -2.0    , 2.0      );
   CreateUserTH1D( "PileupWeight"             , 200 , -2.0    , 2.0      );

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

   CreateUserTH1D( "nVertex_PAS"           ,    31   , -0.5   , 30.5	 ) ; 
   CreateUserTH1D( "nVertex_good_PAS"      ,    31   , -0.5   , 30.5	 ) ; 
   
   CreateUserTH1D( "MTenu_50_110", 200, 40, 140 );

   CreateUserTH1D( "run_PAS"               ,    20000 , 160300  , 180300 );
   CreateUserTH1D( "run_HLT"               ,    20000 , 160300  , 180300 );
   
   //--------------------------------------------------------------------------
   // Loop over the chain
   //--------------------------------------------------------------------------

   if (fChain == 0) return;
   std::cout << "fChain = " << fChain << std::endl;

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
     
     // double pileup_weight = lumiWeights.weight3D (nPileUpInt_BXminus1, nPileUpInt_BX0, nPileUpInt_BXplus1 );
     // if ( isData ) pileup_weight = 1.0;
     // setPileupWeight ( pileup_weight ) ;

     int NPILEUP_AVE = int( nPileUpInt_BX0 );
     int NPILEUP_FINAL = min( NPILEUP_AVE , 25 );
     double pileup_weight = getPileupWeight ( NPILEUP_FINAL, isData ) ;
     
     // std::cout << "nPileUpInt_BX0 = " << nPileUpInt_BX0 << ", pileup_weight = " << pileup_weight << std::endl;
     // std::cout << "  Data file: " <<       pileupReweighter_.getPileupDataFile () << std::endl;
     // std::cout << "  MC   file: " <<       pileupReweighter_.getPileupMCFile   () << std::endl;
     // pileupReweighter_.printPileupWeights();
     // std::cout << "MC PDF at npileup = 11: " << pileupReweighter_.getMCPDF ( 11 )  << std::endl;
     
     //--------------------------------------------------------------------------
     // Get information about gen-level reweighting (should be for Sherpa only)
     //--------------------------------------------------------------------------

     double gen_weight = Weight;
     if ( isData ) gen_weight = 1.0;

     //--------------------------------------------------------------------------
     // Check good run list
     //--------------------------------------------------------------------------
     
     int passedJSON = passJSON ( run, ls , isData ) ;
     if ( !isData ) passedJSON = 1;
     
     //--------------------------------------------------------------------------
     // Check HLT
     //--------------------------------------------------------------------------
     
     int passedHLT = 1;
     if ( isData ){ 
       
       passedHLT = 0;
       
       if ( run < 175771 ){ // 2011A 
	 if ( H_27_CIdVT_CIsT_TIdT_TIsT_1  == 1|| // SingleElectron, 2011A
	      H_27_CIdVT_CIsT_TIdT_TIsT_2  == 1|| // SingleElectron, 2011A
	      H_27_CIdVT_CIsT_TIdT_TIsT_3  == 1|| // SingleElectron, 2011A
	      H_32_CIdVT_CIsT_TIdT_TIsT_3  == 1|| // SingleElectron, 2011A
	      H_17_HEEP_CJ30_25_MHT15_2    == 1|| // ElectronHad, 2011A
	      H_25_HEEP_CJ30_25_MHT20_4    == 1|| // ElectronHad, 2011A
	      H_22_HEEP_CJ30_25_MHT20_2    == 1|| // ElectronHad, 2011A
	      H_22_HEEP_CJ30_25_MHT20_4    == 1|| // ElectronHad, 2011A
	      H_27_HEEP_CJ30_25_MHT20_2    == 1){ // ElectronHad, 2011A
	   passedHLT = 1;
	 }
       }
       
       else if ( run > 175771 ){ // 2011 B 
	 if ( H_Ele30_HEEP_2CJ30_MHT25_1 == 1 ||  // ElectronHad, 2011B
	      H_Ele27_WP80_2CJ25_MHT15_4 == 1 ||  // ElectronHad, 2011B
	      H_Ele27_WP80_2CJ25_MHT15_5 == 1 ) { // ElectronHad, 2011B
	   passedHLT = 1;
	 }
       }
     }


     //--------------------------------------------------------------------------
     // Is this the dataset we want?
     //--------------------------------------------------------------------------
     
     int PassDataset = 1;
     if ( isData ) { 
       PassDataset = 0;
       if ( select2011A ){ 
	 if ( run >= 160329 && run <= 175770 ) PassDataset = 1;
       }
       if ( select2011B ){
	 if ( run >= 175832 && run <= 180296 ) PassDataset = 1;
       }
       if ( select2011 ) {
	 if ( run >= 160329 && run <= 180296 ) PassDataset = 1;
       }
     }

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "Reweighting"              , 1                       , gen_weight );
     fillVariableWithValue(   "PassJSON"                 , passedJSON              , gen_weight ); 
     				
     // Pass dataset variable     
     fillVariableWithValue(   "PassDataset"              , PassDataset             , gen_weight );

     // HLT variable							           
     fillVariableWithValue(   "PassHLT"                  , passedHLT               , gen_weight );
     
     // Filters
     fillVariableWithValue(   "PassHBHENoiseFilter"      , PassHBHENoiseFilter     , gen_weight );
     fillVariableWithValue(   "PassBeamHaloFilterTight"  , PassBeamHaloFilterTight , gen_weight );
									      
     // Muon variables ( for veto ) 					      
     fillVariableWithValue(   "nMuon"                    , nMuon_Ana               , gen_weight );
			                                      		                
     // 1st Electron variables				      		                
     fillVariableWithValue(   "nEle"                     , nEle_Ana                , gen_weight ); 
     fillVariableWithValue(   "Ele1_Pt"                  , Ele1_Pt                 , gen_weight );
     fillVariableWithValue(   "Ele1_Eta"                 , Ele1_Eta                , gen_weight );
									           
     // MET variables	                                      		           
     fillVariableWithValue(   "MET"                      , MET_Pt                  , gen_weight );
     fillVariableWithValue(   "mDeltaPhiMETEle"          , mDPhi_METEle1           , gen_weight );
     									           
     // 1st JET variables                                     		           
     fillVariableWithValue(   "nJet"                     , nJet_Ana                , gen_weight );
									           
     // 1st JET variables                                     		           
     if ( nJet_Stored > 0 ) { 						           
       fillVariableWithValue( "Jet1_Pt"                  , Jet1_Pt                 , gen_weight );
       fillVariableWithValue( "Jet1_Eta"                 , Jet1_Eta                , gen_weight );
       fillVariableWithValue( "mDeltaPhiMET1stJet"       , mDPhi_METJet1           , gen_weight );
     }									           
     									           
     // 2nd JET variables                                     		           
     if ( nJet_Stored > 1 ) { 	                                      	           
       fillVariableWithValue( "Jet2_Pt"                  , Jet2_Pt                 , gen_weight );
       fillVariableWithValue( "Jet2_Eta"                 , Jet2_Eta                , gen_weight );
       fillVariableWithValue( "ST"                       , sT_enujj                , gen_weight );
     }

     // 1 electron, 1 jet variables 
     if ( nEle_Ana > 0 && nJet_Ana > 0 ) { 
       fillVariableWithValue ( "DR_Ele1Jet1"             , DR_Ele1Jet1             , gen_weight );
     }

     // 1 electron, 2 jet variables 
     if ( nEle_Ana > 0 && nJet_Ana > 1 ) { 
       fillVariableWithValue ( "DR_Ele1Jet2"             , DR_Ele1Jet2             , gen_weight );
     }
     
     // Dummy variables
     fillVariableWithValue ("preselection",1, gen_weight ); 

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
     bool passed_minimum      = ( passedAllPreviousCuts("PassBeamHaloFilterTight") && passedCut ("PassBeamHaloFilterTight"));

     if ( passed_minimum && isData ){ 
       FillUserTH1D ("run_HLT", run );
     }

       
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

       double min_DR_EleJet = 999.0;
       double DR_Ele1Jet3 = 999.0;
       if ( nJet_Ana > 2 ) {
	 TLorentzVector ele1,  jet3;
	 ele1.SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
	 jet3.SetPtEtaPhiM ( Jet3_Pt, Jet3_Eta, Jet3_Phi, 0.0 );	 
	 DR_Ele1Jet3 = ele1.DeltaR ( jet3 ) ;
       }

       if ( DR_Ele1Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet1;
       if ( DR_Ele1Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet2;
       if ( nJet_Ana > 2 ) {
	 if ( DR_Ele1Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet3;
       }
              
       if ( isData )        FillUserTH1D("run_PAS", run ) ;
       FillUserTH1D( "MT_charged_enu_PAS"         , MTCharged                      , pileup_weight * gen_weight);
       FillUserTH1D( "MT_type1_enu_PAS"           , MTType1                        , pileup_weight * gen_weight);
       FillUserTH1D( "nElectron_PAS"              , nEle_Ana                       , pileup_weight * gen_weight); 
       FillUserTH1D( "nMuon_PAS"                  , nMuon_Ana                      , pileup_weight * gen_weight); 
       FillUserTH1D( "Pt1stEle_PAS"	          , Ele1_Pt                        , pileup_weight * gen_weight); 
       FillUserTH1D( "Eta1stEle_PAS"	          , Ele1_Eta                       , pileup_weight * gen_weight);
       FillUserTH1D( "Phi1stEle_PAS"	          , Ele1_Phi                       , pileup_weight * gen_weight);
       FillUserTH1D( "Charge1stEle_PAS"           , Ele1_Charge                    , pileup_weight * gen_weight);   
       FillUserTH1D( "METSig_PAS"	          , PFMETSig                       , pileup_weight * gen_weight);   
       FillUserTH1D( "MET_PAS"                    , MET_Pt                         , pileup_weight * gen_weight);
       FillUserTH1D( "METPhi_PAS"	          , MET_Phi                        , pileup_weight * gen_weight);   
       FillUserTH1D( "METCharged_PAS"             , PFMETCharged                   , pileup_weight * gen_weight);
       FillUserTH1D( "METChargedPhi_PAS"          , PFMETChargedPhi                , pileup_weight * gen_weight);   
       FillUserTH1D( "METType1_PAS"               , PFMETType1Cor                  , pileup_weight * gen_weight);
       FillUserTH1D( "METType1Phi_PAS"            , PFMETPhiType1Cor               , pileup_weight * gen_weight);   
       FillUserTH1D( "minMETPt1stEle_PAS"         , TMath::Min ( Ele1_Pt, MET_Pt  ), pileup_weight * gen_weight);
       FillUserTH1D( "Pt1stJet_PAS"               , Jet1_Pt                        , pileup_weight * gen_weight);
       FillUserTH1D( "Pt2ndJet_PAS"               , Jet2_Pt                        , pileup_weight * gen_weight);
       FillUserTH1D( "Eta1stJet_PAS"              , Jet1_Eta                       , pileup_weight * gen_weight);
       FillUserTH1D( "Eta2ndJet_PAS"              , Jet2_Eta                       , pileup_weight * gen_weight);
       FillUserTH1D( "Phi1stJet_PAS"              , Jet1_Phi                       , pileup_weight * gen_weight);
       FillUserTH1D( "Phi2ndJet_PAS"	          , Jet2_Phi                       , pileup_weight * gen_weight);
       FillUserTH1D( "TCHE1stJet_PAS"             , Jet1_btagTCHE                  , pileup_weight * gen_weight);
       FillUserTH1D( "TCHE2ndJet_PAS"             , Jet2_btagTCHE                  , pileup_weight * gen_weight);
       FillUserTH1D( "nMuon_PtCut_IDISO_PAS"      , nMuon_Ana                      , pileup_weight * gen_weight); 
       FillUserTH1D( "MTenu_PAS"                  , MT_Ele1MET                     , pileup_weight * gen_weight);
       FillUserTH1D( "Ptenu_PAS"	          , Pt_Ele1MET                     , pileup_weight * gen_weight);
       FillUserTH1D( "sTlep_PAS"                  , Ele1_Pt + MET_Pt               , pileup_weight * gen_weight);
       FillUserTH1D( "sTjet_PAS"                  , Jet1_Pt + Jet2_Pt              , pileup_weight * gen_weight);
       FillUserTH1D( "sT_PAS"                     , sT_enujj                       , pileup_weight * gen_weight);
       FillUserTH1D( "Mjj_PAS"	                  , M_j1j2                         , pileup_weight * gen_weight);   
       FillUserTH1D( "DCotTheta1stEle_PAS"        , Ele1_DCotTheta                 , pileup_weight * gen_weight);
       FillUserTH1D( "Dist1stEle_PAS"             , Ele1_Dist                      , pileup_weight * gen_weight);
       FillUserTH1D( "mDPhi1stEleMET_PAS"         , mDPhi_METEle1                  , pileup_weight * gen_weight);
       FillUserTH1D( "mDPhi1stJetMET_PAS"         , mDPhi_METJet1                  , pileup_weight * gen_weight);
       FillUserTH1D( "mDPhi2ndJetMET_PAS"         , mDPhi_METJet2                  , pileup_weight * gen_weight); 
       FillUserTH1D( "Mej1_PAS"                   , M_e1j1                         , pileup_weight * gen_weight);
       FillUserTH1D( "Mej2_PAS"                   , M_e1j2                         , pileup_weight * gen_weight);
       FillUserTH1D( "DR_Ele1Jet1_PAS"	          , DR_Ele1Jet1                    , pileup_weight * gen_weight);
       FillUserTH1D( "DR_Ele1Jet2_PAS"	          , DR_Ele1Jet2                    , pileup_weight * gen_weight);
       FillUserTH1D( "DR_Jet1Jet2_PAS"	          , DR_Jet1Jet2                    , pileup_weight * gen_weight);
       FillUserTH1D( "minDR_EleJet_PAS"           , min_DR_EleJet                  , pileup_weight * gen_weight);
       FillUserTH1D( "nVertex_PAS"                , nVertex                        , pileup_weight * gen_weight);
       FillUserTH1D( "nVertex_good_PAS"           , nVertex_good                   , pileup_weight * gen_weight);
       FillUserTH1D( "MatchPhotonConv1stEle_PAS"  , Ele1_MatchPhotConv             , pileup_weight * gen_weight);
       FillUserTH1D( "MatchPhotonConv2ndEle_PAS"  , Ele2_MatchPhotConv             , pileup_weight * gen_weight);
       FillUserTH1D( "GeneratorWeight"       , gen_weight                     );
       FillUserTH1D( "PileupWeight"          , pileup_weight                  ); 
       
       if ( MT_Ele1MET > 50 && MT_Ele1MET < 110 ){
	 FillUserTH1D( "MTenu_50_110", MT_Ele1MET,  pileup_weight * gen_weight ) ;
       }
       
       if ( nVertex_good >= 0 && nVertex_good <= 3 ) {
	 FillUserTH1D( "MT_GoodVtxLTE3_PAS"              , MT_Ele1MET, pileup_weight * gen_weight ) ;
	 FillUserTH1D( "MTCharged_GoodVtxLTE3_PAS"       , MTCharged , pileup_weight * gen_weight ) ;
	 FillUserTH1D( "MTType1_GoodVtxLTE3_PAS"         , MTType1   , pileup_weight * gen_weight ) ;
       }						 
       							 
       if ( nVertex_good >= 4 && nVertex_good <= 8 ) {	 
	 FillUserTH1D( "MT_GoodVtxGTE4_LTE8_PAS"         , MT_Ele1MET, pileup_weight * gen_weight ) ;
	 FillUserTH1D( "MTCharged_GoodVtxGTE4_LTE8_PAS"  , MTCharged , pileup_weight * gen_weight ) ;
	 FillUserTH1D( "MTType1_GoodVtxGTE4_LTE8_PAS"    , MTType1   , pileup_weight * gen_weight ) ;
       }
       
       if ( nVertex_good >= 9 && nVertex_good <= 15) {
	 FillUserTH1D( "MT_GoodVtxGTE9_LTE15_PAS"        , MT_Ele1MET, pileup_weight * gen_weight ) ;
	 FillUserTH1D( "MTCharged_GoodVtxGTE9_LTE15_PAS" , MTCharged , pileup_weight * gen_weight ) ;
	 FillUserTH1D( "MTType1_GoodVtxGTE9_LTE15_PAS"   , MTType1   , pileup_weight * gen_weight ) ;
       }
       
       if ( nVertex_good >= 16                     ) {
	 FillUserTH1D( "MT_GoodVtxGTE16_PAS"             , MT_Ele1MET, pileup_weight * gen_weight ) ;
	 FillUserTH1D( "MTCharged_GoodVtxGTE16_PAS"      , MTCharged , pileup_weight * gen_weight ) ;
	 FillUserTH1D( "MTType1_GoodVtxGTE16_PAS"        , MTType1   , pileup_weight * gen_weight ) ;
       }
              
     }
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
