#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TProfile.h>

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
   fillAllOtherCuts                 ( true  ) ;
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;

   //--------------------------------------------------------------------------
   // Any extra features
   //--------------------------------------------------------------------------
   
   TProfile * profile_run_vs_nvtx_HLT = new TProfile("run_vs_nvtx_HLT", "", 20000 , 160300  , 180300 );
   TProfile * profile_run_vs_nvtx_PAS = new TProfile("run_vs_nvtx_PAS", "", 20000 , 160300  , 180300 );
   
   //--------------------------------------------------------------------------
   // Get pre-cut values
   //--------------------------------------------------------------------------

   // eta boundaries

   double eleEta_bar_max = getPreCutValue1("eleEta_bar");
   double eleEta_end_min = getPreCutValue1("eleEta_end1");
   double eleEta_end_max = getPreCutValue2("eleEta_end2");

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------
   
   CreateUserTH1D( "ProcessID"             ,    21  , -0.5    , 20.5     );
   CreateUserTH1D( "ProcessID_PAS"         ,    21  , -0.5    , 20.5     );
   CreateUserTH1D( "ProcessID_ZWindow"     ,    21  , -0.5    , 20.5     );
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
   CreateUserTH1D( "Pt1stJet_PAS"          ,    100 , 0       , 1000	 ); 
   CreateUserTH1D( "Pt2ndJet_PAS"          ,    100 , 0       , 1000	 ); 
   CreateUserTH1D( "Eta1stJet_PAS"         ,    100 , -5      , 5	 ); 
   CreateUserTH1D( "Eta2ndJet_PAS"         ,    100 , -5      , 5	 ); 
   CreateUserTH1D( "Phi1stJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "Phi2ndJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
   CreateUserTH1D( "sTlep_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sTjet_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PAS"                ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "sT_PASandMee100"       ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mjj_PAS"		   ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mee_PAS"		   ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mee_PASandST445"       ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "MTenu_PAS"                , 200 , 0       , 1000	 ); 
   CreateUserTH1D( "Me1j1_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me1j2_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me2j1_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me2j2_PAS"             ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me1j_selected_PAS"     ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Me2j_selected_PAS"     ,    200 , 0       , 2000     );
   CreateUserTH1D( "Mej_selected_avg_PAS"  ,    200 , 0       , 2000     ); 
   CreateUserTH1D( "Mej_selected_min_PAS"  ,    200 , 0       , 2000     ); 
   CreateUserTH1D( "Mej_selected_max_PAS"  ,    200 , 0       , 2000     ); 
   CreateUserTH1D( "Mej_minmax_PAS"        ,    200 , 0       , 2000     ); 
   CreateUserTH1D( "Meejj_PAS"             ,    200 , 0       , 2000     );
   CreateUserTH1D( "run_PAS"               ,    20000 , 160300  , 180300 );
   CreateUserTH1D( "run_HLT"               ,    20000 , 160300  , 180300 );
		                           
   CreateUserTH1D( "Ptee_PAS"              ,    200 , 0       , 2000     );
   		                           
   CreateUserTH1D( "DCotTheta1stEle_PAS"   ,    100 , 0.0, 1.0);
   CreateUserTH1D( "Dist1stEle_PAS"        ,    100 , 0.0, 1.0);  
   CreateUserTH1D( "DCotTheta2ndEle_PAS"   ,    100 , 0.0, 1.0);
   CreateUserTH1D( "Dist2ndEle_PAS"        ,    100 , 0.0, 1.0);  
		                           
   CreateUserTH1D( "nVertex_PAS"           ,    101   , -0.5   , 100.5	 ) ; 
   
   CreateUserTH1D( "DR_Ele1Jet1_PAS"	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele1Jet2_PAS"	   , 	getHistoNBins("DR_Ele1Jet2"), getHistoMin("DR_Ele1Jet2"), getHistoMax("DR_Ele1Jet2")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet1_PAS"	   , 	getHistoNBins("DR_Ele2Jet1"), getHistoMin("DR_Ele2Jet1"), getHistoMax("DR_Ele2Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet2_PAS"	   , 	getHistoNBins("DR_Ele2Jet2"), getHistoMin("DR_Ele2Jet2"), getHistoMax("DR_Ele2Jet2")     ) ; 
   CreateUserTH1D( "DR_Jet1Jet2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
   CreateUserTH1D( "DR_Ele1Ele2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
   CreateUserTH1D( "minDR_EleJet_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 

   CreateUserTH2D( "Me1jVsMe2j_selected",     200, 0, 2000, 200, 0, 2000) ;
   CreateUserTH2D( "Me1jVsMe2j_rejected",     200, 0, 2000, 200, 0, 2000) ;
   CreateUserTH2D( "sT_vs_Mee"          ,     200, 250, 750, 200, 50, 450) ;

   CreateUserTH1D( "MTeemunu_PAS"          ,    200 , 0       , 1000	 ); 

   CreateUserTH1D( "Mee_80_100_Preselection", 200, 60, 120 );
   CreateUserTH1D( "Mee_70_110_Preselection", 200, 60, 120 );

   CreateUserTH1D( "Mee_70_110_Preselection_Process0", 200, 60, 120 );
   CreateUserTH1D( "Mee_70_110_Preselection_Process1", 200, 60, 120 );
   CreateUserTH1D( "Mee_70_110_Preselection_Process2", 200, 60, 120 );
   CreateUserTH1D( "Mee_70_110_Preselection_Process3", 200, 60, 120 );
   CreateUserTH1D( "Mee_70_110_Preselection_Process4", 200, 60, 120 );

   CreateUserTH1D( "Mee_EBEB_PAS"		   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EBEE_PAS"		   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EEEE_PAS"		   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EB_PAS" 		   ,    60 , 60       , 120	 ); 

   CreateUserTH1D( "Mee_EBEB_80_100_PAS"	   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EBEE_80_100_PAS"	   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EEEE_80_100_PAS"	   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EB_80_100_PAS" 	     	   ,    60 , 60       , 120	 ); 
   
   CreateUserTH1D( "Mee_EBEB_70_110_PAS"	   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EBEE_70_110_PAS"	   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EEEE_70_110_PAS"	   ,    60 , 60       , 120	 ); 
   CreateUserTH1D( "Mee_EB_70_110_PAS" 	     	   ,    60 , 60       , 120	 ); 

   CreateUserTH1D( "PileupWeight"   , 100, -10, 10 );
   CreateUserTH1D( "GeneratorWeight", 100, -2.0, 2.0);

   //--------------------------------------------------------------------------
   // Loop over the chain
   //--------------------------------------------------------------------------

   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntries();
   //Long64_t nentries = 5;
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
     
     double pileup_weight = getPileupWeight ( nPileUpInt_True, isData ) ;
     
     //--------------------------------------------------------------------------
     // Get information about gen-level reweighting (should be for Sherpa only)
     //--------------------------------------------------------------------------

     double gen_weight = Weight;
     if ( isData ) gen_weight = 1.0;

     // std::cout << "Gen weight = " << int ( 1.0 / gen_weight ) << std::endl;

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "PassJSON"                      , passedJSON, gen_weight * pileup_weight   ) ; 

     // Noise filters
     fillVariableWithValue(   "PassHBHENoiseFilter"	      , PassHBHENoiseFilter                              , gen_weight * pileup_weight );
     fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight                          , gen_weight * pileup_weight );
     fillVariableWithValue(   "PassBadEESupercrystalFilter"   , ( isData == 1 ) ? PassBadEESupercrystalFilter : 1, gen_weight * pileup_weight );
     fillVariableWithValue(   "PassBeamScraping"	      , ( isData == 1 ) ? PassBeamScraping	      : 1, gen_weight * pileup_weight );
     fillVariableWithValue(   "PassEcalDeadCellBoundEnergy"   , PassEcalDeadCellBoundEnergy                      , gen_weight * pileup_weight );
     fillVariableWithValue(   "PassEcalDeadCellTrigPrim"      , PassEcalDeadCellTrigPrim                         , gen_weight * pileup_weight );
     fillVariableWithValue(   "PassEcalLaserCorrFilter"       , ( isData == 1 ) ? PassEcalLaserCorrFilter     : 1, gen_weight * pileup_weight );
     fillVariableWithValue(   "PassHcalLaserEventFilter"      , ( isData == 1 ) ? PassHcalLaserEventFilter    : 1, gen_weight * pileup_weight );
     fillVariableWithValue(   "PassPhysDecl"		      , ( isData == 1 ) ? PassPhysDecl		      : 1, gen_weight * pileup_weight );
     fillVariableWithValue(   "PassPrimaryVertex"	      , PassPrimaryVertex                                , gen_weight * pileup_weight );
     fillVariableWithValue(   "PassTrackingFailure"	      , ( isData == 1 ) ? PassTrackingFailure	      : 1, gen_weight * pileup_weight );
     
     //fillVariableWithValue(   "PassBPTX0"                     , PassBPTX0                  , gen_weight  ) ; 
     //fillVariableWithValue(   "PassPhysDecl"                  , PassPhysDecl               , gen_weight  ) ; 
     //fillVariableWithValue(   "PassBeamScraping"              , PassBeamScraping           , gen_weight  ) ; 
     //fillVariableWithValue(   "PassPrimaryVertex"             , PassPrimaryVertex          , gen_weight  ) ; 
     //fillVariableWithValue(   "PassBeamHaloFilterLoose"	      , PassBeamHaloFilterLoose	   , gen_weight  ) ; 
     //fillVariableWithValue(   "PassTrackingFailure"           , PassTrackingFailure        , gen_weight  ) ; 
     //fillVariableWithValue(   "PassCaloBoundaryDRFilter"      , PassCaloBoundaryDRFilter   , gen_weight  ) ; 
     //fillVariableWithValue(   "PassEcalMaskedCellDRFilter"    , PassEcalMaskedCellDRFilter , gen_weight  ) ; 

     // Fill HLT

     int passHLT = 1;
     if ( isData ) { 
       passHLT = 0;
       if ( H_Ele30_PFJet100_25      == 1 ||
	    H_Ele30_PFNoPUJet100_25  == 1 ){
       	 passHLT = 1;
       }
     }

     int nEle_hltMatched = 0.0;
     if ( Ele1_hltEleSignalPt > 0.0 ) nEle_hltMatched++;
     if ( Ele2_hltEleSignalPt > 0.0 ) nEle_hltMatched++;
     
     int nJet_hltMatched = 0.0;
     if ( Jet1_hltNoPUJetPt > 0.0 || Jet1_hltJetPt > 0.0 ) nJet_hltMatched++;
     if ( Jet2_hltNoPUJetPt > 0.0 || Jet2_hltJetPt > 0.0 ) nJet_hltMatched++;
     
     fillVariableWithValue("nEle_hltMatched",nEle_hltMatched, gen_weight * pileup_weight  );
     fillVariableWithValue("nJet_hltMatched",nJet_hltMatched, gen_weight * pileup_weight  );

     // Electrons
     int PassNEle = 0;
     if ( nEle_ptCut == 2 ) PassNEle = 1;
     
     // Muons
     int PassNMuon = 0;
     if ( nMuon_ptCut == 0 ) PassNMuon = 1;
     
     fillVariableWithValue ( "Reweighting", 1, gen_weight * pileup_weight  );
     fillVariableWithValue ( "PassHLT", passHLT, gen_weight * pileup_weight  ) ;
     
     // Electrons
     fillVariableWithValue(   "PassNEle"                      , PassNEle    , gen_weight * pileup_weight  ) ;
     if ( nEle_store >= 1 ) { 
       fillVariableWithValue( "Ele1_Pt"                       , Ele1_Pt     , gen_weight * pileup_weight  ) ;
     }									    
     if ( nEle_store >= 2 ) { 						    
       fillVariableWithValue( "Ele2_Pt"                       , Ele2_Pt     , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "M_e1e2"                        , M_e1e2      , gen_weight * pileup_weight  ) ;
     }									    
									    
     // Jets								    
     fillVariableWithValue(   "nJet"                          , nJet_ptCut  , gen_weight * pileup_weight  ) ;
     if ( nJet_store >= 1 ) { 						    
       fillVariableWithValue( "Jet1_Pt"                       , Jet1_Pt     , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "Jet1_Eta"                      , Jet1_Eta    , gen_weight * pileup_weight  ) ;
     }
     if ( nJet_store >= 2 ) { 
       fillVariableWithValue( "Jet2_Pt"                       , Jet2_Pt     , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "Jet2_Eta"                      , Jet2_Eta    , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "DR_Jet1Jet2"                   , DR_Jet1Jet2 , gen_weight * pileup_weight  ) ;
     }

     // Muons
     fillVariableWithValue(   "PassNMuon"                     , PassNMuon   , gen_weight * pileup_weight  ) ;

     // DeltaR
     if ( nEle_store >= 2 && nJet_store >= 1) {
       fillVariableWithValue( "DR_Ele1Jet1"                   , DR_Ele1Jet1 , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "DR_Ele2Jet1"                   , DR_Ele2Jet1 , gen_weight * pileup_weight  ) ;
       if(nJet_store >= 2) {
	 fillVariableWithValue( "DR_Ele1Jet2"                 , DR_Ele1Jet2 , gen_weight * pileup_weight  ) ;
	 fillVariableWithValue( "DR_Ele2Jet2"                 , DR_Ele2Jet2 , gen_weight * pileup_weight  ) ;
       }
     }

     // sT
     double M_ej_avg;
     double M_ej_min;
     double M_ej_max;

     if ( nEle_store >= 2 && nJet_store >= 2) {
       if ( fabs(M_e1j1-M_e2j2) < fabs(M_e1j2-M_e2j1) )  {
	 M_ej_avg = (M_e1j1 + M_e2j2) / 2.0;
	 if    ( M_e1j1 < M_e2j2 ) { M_ej_min = M_e1j1; M_ej_max = M_e2j2; }
	 else                      { M_ej_min = M_e2j2; M_ej_max = M_e1j1; }
       }
       else { 
	 M_ej_avg = (M_e1j2 + M_e2j1) / 2.0;
	 if    ( M_e1j2 < M_e2j1 ) { M_ej_min = M_e1j2; M_ej_max = M_e2j1; }
	 else                      { M_ej_min = M_e2j1; M_ej_max = M_e1j2; }
       }
       
       fillVariableWithValue( "sT_eejj"                      , sT_eejj , gen_weight * pileup_weight  ) ;
     }      

     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();
     
     //--------------------------------------------------------------------------
     // Fill skim-level plots
     //--------------------------------------------------------------------------
     
     bool passed_minimum = ( passedAllPreviousCuts("PassTrackingFailure") && passedCut ("PassTrackingFailure"));

     if ( passed_minimum && isData ){ 
       FillUserTH1D ("run_HLT", run );
       profile_run_vs_nvtx_HLT -> Fill ( run, nVertex, 1 ) ;
     }

     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------

     FillUserTH1D( "PileupWeight"   , pileup_weight );
     FillUserTH1D( "GeneratorWeight", gen_weight ) ;

     bool passed_preselection = ( passedAllPreviousCuts("M_e1e2") && passedCut ("M_e1e2") );

     /*
     if ( isData ) {
       std::cout.precision(0);
       std::cout << fixed <<  "Run = " << run << ", event = " << event << ", ls = " << ls << std::endl;
       std::cout.precision(3);
       std::cout << fixed <<  "  Mej      = " << M_ej_avg << std::endl;
       std::cout << fixed <<  "  Mee      = " << M_e1e2 << std::endl;
       std::cout << fixed <<  "  sT       = " << sT_enujj << std::endl;
       std::cout << fixed <<  "  Ele1 Pt  = " << Ele1_Pt << "\t, Eta = " << Ele1_Eta << "\t, Phi = " << Ele1_Phi << std::endl;
       std::cout << fixed <<  "  Ele2 Pt  = " << Ele2_Pt << "\t, Eta = " << Ele2_Eta << "\t, Phi = " << Ele2_Phi << std::endl;
       std::cout << fixed <<  "  Jet1 Pt  = " << Jet1_Pt << "\t, Eta = " << Jet1_Eta << "\t, Phi = " << Jet1_Phi << std::endl;
       std::cout << fixed <<  "  Jet2 Pt  = " << Jet2_Pt << "\t, Eta = " << Jet2_Eta << "\t, Phi = " << Jet2_Phi << std::endl;
     }
     */
     
     if ( passed_preselection ) {

       FillUserTH1D( "ProcessID_PAS"      , ProcessID, pileup_weight * gen_weight ) ;
       
       //--------------------------------------------------------------------------
       // Fill skim tree, if necessary
       //--------------------------------------------------------------------------
       
       bool isEB1 = ( fabs(Ele1_Eta) < eleEta_bar_max  ) ;
       bool isEE1 = ( fabs(Ele1_Eta) > eleEta_end_min &&
		      fabs(Ele1_Eta) < eleEta_end_max ) ;
       
       bool isEB2 = ( fabs(Ele2_Eta) < eleEta_bar_max  ) ;
       bool isEE2 = ( fabs(Ele2_Eta) > eleEta_end_min &&
		      fabs(Ele2_Eta) < eleEta_end_max ) ;

       bool isEBEB = ( isEB1 && isEB2 ) ;
       bool isEBEE = ( ( isEB1 && isEE2 ) ||
		       ( isEE1 && isEB2 ) );
       bool isEEEE = ( isEE1  && isEE2  );
       bool isEB   = ( isEBEB || isEBEE );
       
       if ( isData ) { 
	 FillUserTH1D("run_PAS", run ) ;
	 profile_run_vs_nvtx_PAS -> Fill ( run, nVertex, 1 ) ;
       }
       
       FillUserTH1D("nElectron_PAS"        , nEle_ptCut         , pileup_weight * gen_weight ) ;
       FillUserTH1D("nMuon_PAS"            , nMuon_ptCut        , pileup_weight * gen_weight ) ;
       FillUserTH1D("nJet_PAS"             , nJet_ptCut         , pileup_weight * gen_weight ) ;
       FillUserTH1D("Pt1stEle_PAS"	   , Ele1_Pt            , pileup_weight * gen_weight ) ;
       FillUserTH1D("Eta1stEle_PAS"	   , Ele1_Eta           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Phi1stEle_PAS"	   , Ele1_Phi           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Pt2ndEle_PAS"	   , Ele2_Pt            , pileup_weight * gen_weight ) ;
       FillUserTH1D("Eta2ndEle_PAS"	   , Ele2_Eta           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Phi2ndEle_PAS"	   , Ele2_Phi           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Charge1stEle_PAS"	   , Ele1_Charge        , pileup_weight * gen_weight ) ;
       FillUserTH1D("Charge2ndEle_PAS"	   , Ele2_Charge        , pileup_weight * gen_weight ) ;
       FillUserTH1D("MET_PAS"              , PFMET_Type01XY_Pt  , pileup_weight * gen_weight ) ;
       FillUserTH1D("METPhi_PAS"	   , PFMET_Type01XY_Phi , pileup_weight * gen_weight ) ;
       FillUserTH1D("Pt1stJet_PAS"         , Jet1_Pt            , pileup_weight * gen_weight ) ;
       FillUserTH1D("Pt2ndJet_PAS"         , Jet2_Pt            , pileup_weight * gen_weight ) ;
       FillUserTH1D("Eta1stJet_PAS"        , Jet1_Eta           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Eta2ndJet_PAS"        , Jet2_Eta           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Phi1stJet_PAS"	   , Jet1_Phi           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Phi2ndJet_PAS"	   , Jet2_Phi           , pileup_weight * gen_weight ) ;
       FillUserTH1D("sTlep_PAS"            , Ele1_Pt + Ele2_Pt  , pileup_weight * gen_weight ) ;
       FillUserTH1D("sTjet_PAS"            , Jet1_Pt + Jet2_Pt  , pileup_weight * gen_weight ) ;
       FillUserTH1D("sT_PAS"               , sT_eejj            , pileup_weight * gen_weight ) ;
       FillUserTH1D("Mjj_PAS"		   , M_j1j2             , pileup_weight * gen_weight ) ;
       FillUserTH1D("Mee_PAS"		   , M_e1e2             , pileup_weight * gen_weight ) ;
       FillUserTH1D( "MTenu_PAS"           , MT_Ele1MET         , pileup_weight * gen_weight ) ;
       FillUserTH1D("Me1j1_PAS"            , M_e1j1             , pileup_weight * gen_weight ) ;
       FillUserTH1D("Me1j2_PAS"            , M_e1j2             , pileup_weight * gen_weight ) ;
       FillUserTH1D("Me2j1_PAS"            , M_e2j1             , pileup_weight * gen_weight ) ;
       FillUserTH1D("Me2j2_PAS"            , M_e2j2             , pileup_weight * gen_weight ) ;
       FillUserTH1D("Ptee_PAS"             , Pt_e1e2            , pileup_weight * gen_weight ) ;
       FillUserTH1D("DCotTheta1stEle_PAS"  , Ele1_DCotTheta     , pileup_weight * gen_weight ) ;
       FillUserTH1D("Dist1stEle_PAS"       , Ele1_Dist          , pileup_weight * gen_weight ) ;
       FillUserTH1D("DCotTheta2ndEle_PAS"  , Ele2_DCotTheta     , pileup_weight * gen_weight ) ;
       FillUserTH1D("Dist2ndEle_PAS"       , Ele2_Dist          , pileup_weight * gen_weight ) ;
       FillUserTH1D("nVertex_PAS"          , nVertex            , pileup_weight * gen_weight ) ;
       FillUserTH1D("DR_Ele1Jet1_PAS"	   , DR_Ele1Jet1        , pileup_weight * gen_weight ) ;
       FillUserTH1D("DR_Ele1Jet2_PAS"	   , DR_Ele1Jet2        , pileup_weight * gen_weight ) ;
       FillUserTH1D("DR_Ele2Jet1_PAS"	   , DR_Ele2Jet1        , pileup_weight * gen_weight ) ;
       FillUserTH1D("DR_Ele2Jet2_PAS"	   , DR_Ele2Jet2        , pileup_weight * gen_weight ) ;
       FillUserTH1D("DR_Jet1Jet2_PAS"	   , DR_Jet1Jet2        , pileup_weight * gen_weight ) ;
       FillUserTH2D("sT_vs_Mee"            , sT_eejj, M_e1e2    , pileup_weight * gen_weight ) ;
       
       if ( sT_eejj > 445. ) { 
	 FillUserTH1D( "Mee_PASandST445", M_e1e2, pileup_weight * gen_weight ) ;
       }
       
       if ( M_e1e2 > 100. ) { 
	 FillUserTH1D("sT_PASandMee100", sT_eejj, pileup_weight * gen_weight ) ;
       }

       if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_PAS", M_e1e2, pileup_weight * gen_weight ); 
       else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_PAS", M_e1e2, pileup_weight * gen_weight ); 
       else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_PAS", M_e1e2, pileup_weight * gen_weight ); 
       if      ( isEB   ) FillUserTH1D( "Mee_EB_PAS"  , M_e1e2, pileup_weight * gen_weight ); 

       if ( M_e1e2 > 80.0 && M_e1e2 < 100.0 ){
	 FillUserTH1D("Mee_80_100_Preselection", M_e1e2, pileup_weight * gen_weight ) ;
	 if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_80_100_PAS", M_e1e2, pileup_weight * gen_weight ); 
	 else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_80_100_PAS", M_e1e2, pileup_weight * gen_weight ); 
	 else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_80_100_PAS", M_e1e2, pileup_weight * gen_weight ); 
	 if      ( isEB   ) FillUserTH1D( "Mee_EB_80_100_PAS"  , M_e1e2, pileup_weight * gen_weight ); 
       } 

       if ( M_e1e2 > 70.0 && M_e1e2 < 110.0 ){
	 FillUserTH1D("Mee_70_110_Preselection", M_e1e2, pileup_weight * gen_weight ) ;
	 if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_70_110_PAS", M_e1e2, pileup_weight * gen_weight ); 
	 else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_70_110_PAS", M_e1e2, pileup_weight * gen_weight ); 
	 else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_70_110_PAS", M_e1e2, pileup_weight * gen_weight ); 
	 if      ( isEB   ) FillUserTH1D( "Mee_EB_70_110_PAS"  , M_e1e2, pileup_weight * gen_weight ); 

	 FillUserTH1D( "ProcessID_ZWindow", ProcessID, pileup_weight * gen_weight );
	 if ( ProcessID == 0 ) FillUserTH1D ( "Mee_70_110_Preselection_Process0", M_e1e2, pileup_weight * gen_weight );
	 if ( ProcessID == 1 ) FillUserTH1D ( "Mee_70_110_Preselection_Process1", M_e1e2, pileup_weight * gen_weight );
	 if ( ProcessID == 2 ) FillUserTH1D ( "Mee_70_110_Preselection_Process2", M_e1e2, pileup_weight * gen_weight );
	 if ( ProcessID == 3 ) FillUserTH1D ( "Mee_70_110_Preselection_Process3", M_e1e2, pileup_weight * gen_weight );
	 if ( ProcessID == 4 ) FillUserTH1D ( "Mee_70_110_Preselection_Process4", M_e1e2, pileup_weight * gen_weight );
       }
       
       double DR_Ele1Jet3 = 999.;
       double DR_Ele2Jet3 = 999.;
       double min_DR_EleJet = 999.;
       TLorentzVector eejj, e1e2mu, e1, j1, e2, j2,j3, mu, met;
       e1.SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
       e2.SetPtEtaPhiM ( Ele2_Pt, Ele2_Eta, Ele2_Phi, 0.0 );
       j1.SetPtEtaPhiM ( Jet1_Pt, Jet1_Eta, Jet1_Phi, 0.0 );
       j2.SetPtEtaPhiM ( Jet2_Pt, Jet2_Eta, Jet2_Phi, 0.0 );
       if ( nJet_store > 2 ) { 
	 j3.SetPtEtaPhiM ( Jet3_Pt, Jet3_Eta, Jet3_Phi, 0.0 );
	 DR_Ele1Jet3 = e1.DeltaR( j3 );
	 DR_Ele2Jet3 = e2.DeltaR( j3 );
       }
       
       if ( DR_Ele1Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet1;
       if ( DR_Ele1Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet2;
       if ( DR_Ele2Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet1;
       if ( DR_Ele2Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet2;
       if ( nJet_store > 2 ) {
	 if ( DR_Ele1Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet3;
	 if ( DR_Ele2Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet3;
       }
       
       FillUserTH1D( "minDR_EleJet_PAS", min_DR_EleJet, pileup_weight * gen_weight);
       
       
       mu.SetPtEtaPhiM ( Muon1_Pt, Muon1_Eta, Muon1_Phi, 0.0 );
       met.SetPtEtaPhiM ( PFMET_Type01XY_Pt, 0.0, PFMET_Type01XY_Phi, 0.0 );
       
       double DR_Ele1Ele2 = e1.DeltaR( e2 ) ;
       FillUserTH1D("DR_Ele1Ele2_PAS"	   , DR_Ele1Ele2      , pileup_weight * gen_weight ) ;
       
       eejj = e1 + e2 + j1 + j2 ; 
       double M_eejj = eejj.M();

       FillUserTH1D("Meejj_PAS", M_eejj , pileup_weight * gen_weight );

       
       if ( fabs(M_e1j1-M_e2j2) < fabs(M_e1j2-M_e2j1) )  {
	 FillUserTH1D("Mej_selected_avg_PAS", M_ej_avg, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Mej_selected_min_PAS", M_ej_min, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Mej_selected_max_PAS", M_ej_max, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Mej_minmax_PAS"      , M_ej_min, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Mej_minmax_PAS"      , M_ej_max, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Me1j_selected_PAS"   , M_e1j1  , pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Me2j_selected_PAS"   , M_e2j2  , pileup_weight * gen_weight) ;	   
	 FillUserTH2D( "Me1jVsMe2j_selected", M_e1j1  , M_e2j2, pileup_weight * gen_weight ) ;
	 FillUserTH2D( "Me1jVsMe2j_rejected", M_e1j2  , M_e2j1, pileup_weight * gen_weight ) ;
       }

       else {
	 FillUserTH1D("Mej_selected_avg_PAS", M_ej_avg, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Mej_selected_min_PAS", M_ej_min, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Mej_selected_max_PAS", M_ej_max, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Me1j_selected_PAS"   , M_e1j2, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Me2j_selected_PAS"   , M_e2j1, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Mej_minmax_PAS"      , M_ej_min, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Mej_minmax_PAS"      , M_ej_max, pileup_weight * gen_weight) ;	   
	 FillUserTH2D( "Me1jVsMe2j_selected", M_e1j2, M_e2j1, pileup_weight * gen_weight ) ;
	 FillUserTH2D( "Me1jVsMe2j_rejected", M_e1j1, M_e2j2, pileup_weight * gen_weight ) ;
       }
     }
   } // End loop over events
   
   output_root_ -> cd();
   profile_run_vs_nvtx_HLT -> Write();
   profile_run_vs_nvtx_PAS -> Write();
   
   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

