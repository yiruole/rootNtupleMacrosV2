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
   // Pileup reweighting initialization
   //--------------------------------------------------------------------------

   // Lumi3DReWeighting lumiWeights = Lumi3DReWeighting("/afs/cern.ch/user/e/eberry/public/LQ_PILEUP/pileup_truth_MC_Summer11_PU_S4_3DReweighting.root",
   // 						     "/afs/cern.ch/user/e/eberry/public/LQ_PILEUP/pileup_truth_finebin_2011_finebin.root",
   // 						     "pileup", "pileup");
   // lumiWeights.weight3D_init(1.0);

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
   // Create TH1D's
   //--------------------------------------------------------------------------
   
   // gap: 1.442 - 1.560

   // eleEta_bar            	1.442        -                    -               -               -1
   // eleEta_end1            	1.560        2.0                  -               -               -1
   // eleEta_end2            	2.000        2.5                  -               -               -1
   
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
   CreateUserTH1D( "METSig_PAS"            ,    200 , 0       , 1000	 ); 
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
   CreateUserTH1D( "Mej_selected_avg_PAS"  ,    200 , 0       , 2000     );
   CreateUserTH1D( "Meejj_PAS"             ,    200 , 0       , 2000     );
   CreateUserTH1D( "run_PAS"               ,    20000 , 160300  , 180300 );
   CreateUserTH1D( "run_HLT"               ,    20000 , 160300  , 180300 );
		                           
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
   CreateUserTH1D( "DR_Ele1Ele2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
   CreateUserTH1D( "minDR_EleJet_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 

   CreateUserTH2D( "Me1jVsMe2j_selected",     200, 0, 2000, 200, 0, 2000) ;
   CreateUserTH2D( "Me1jVsMe2j_rejected",     200, 0, 2000, 200, 0, 2000) ;

   CreateUserTH1D( "MTeemunu_PAS"          ,    200 , 0       , 1000	 ); 

   CreateUserTH1D( "Mee_80_100_Preselection", 200, 60, 120 );
   CreateUserTH1D( "Mee_70_110_Preselection", 200, 60, 120 );

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
     
     // double pileup_weight = lumiWeights.weight3D (nPileUpInt_BXminus1, nPileUpInt_BX0, nPileUpInt_BXplus1 );
     // if ( isData ) pileup_weight = 1.0;
     // setPileupWeight ( pileup_weight ) ;

     int NPILEUP_AVE = int( nPileUpInt_BX0 );
     int NPILEUP_FINAL = min( NPILEUP_AVE , 25 );
     double pileup_weight = getPileupWeight ( NPILEUP_FINAL, isData ) ;
      
     //--------------------------------------------------------------------------
     // Get information about gen-level reweighting (should be for Sherpa only)
     //--------------------------------------------------------------------------

     double gen_weight = Weight;
     if ( isData ) gen_weight = 1.0;

     // std::cout << "Gen weight = " << int ( 1.0 / gen_weight ) << std::endl;

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
	 PassDataset = 1;
       }
     }

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "PassJSON"                      , passedJSON, gen_weight   ) ; 

     // Dataset variable
     fillVariableWithValue(   "PassDataset"                      , PassDataset, gen_weight  ) ;

     // Noise filters
     fillVariableWithValue(   "PassHBHENoiseFilter"           , PassHBHENoiseFilter    , gen_weight   ) ; 
     fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight, gen_weight   ) ; 

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
       
       if ( H_DoubleEle33_CIdL_1  == 1 || 
       	    H_DoubleEle33_CIdL_2  == 1 || 
       	    H_DoubleEle33_CIdL_3  == 1 || 
       	    H_DoubleEle33_CIdL_4  == 1 || 
       	    H_DoubleEle33_CIdL_5  == 1 || 
	    H_DoubleEle33_CIdL_6  == 1 ||
	    H_DoubleEle33_CIdL_7  == 1 ||
	    H_DoubleEle33_CIdT_1  == 1 ||
	    H_DoubleEle33_CIdT_2  == 1 ||
	    H_DoubleEle33_CIdT_3  == 1 ||
       	    H_DoublePhoton33_1    == 1 || 
       	    H_DoublePhoton33_2    == 1 || 
       	    H_DoublePhoton33_3    == 1 || 
       	    H_DoublePhoton33_4    == 1 || 
       	    H_DoublePhoton33_5    == 1 ) {
       	 passHLT = 1;
       }
     }

     // What dataset is this?

     // Electrons
     int PassNEle = 0;
     if ( nEle_Ana == 2 ) PassNEle = 1;
     
     // Muons
     int PassNMuon = 0;
     if ( nMuon_Ana == 0 ) PassNMuon = 1;
     
     fillVariableWithValue ( "Reweighting", 1, gen_weight  );
     fillVariableWithValue ( "PassHLT", passHLT, gen_weight  ) ;
     
     // Electrons
     fillVariableWithValue(   "PassNEle"                      , PassNEle    , gen_weight  ) ;
     if ( nEle_Ana >= 1 ) { 
       fillVariableWithValue( "Ele1_Pt"                       , Ele1_Pt     , gen_weight  ) ;
       fillVariableWithValue( "Ele1_Eta"                      , Ele1_Eta    , gen_weight  ) ;
       fillVariableWithValue( "abs_Ele1_Eta"                  , fabs(Ele1_Eta), gen_weight  );
     }									    
     if ( nEle_Ana >= 2 ) { 						    
       fillVariableWithValue( "Ele2_Pt"                       , Ele2_Pt     , gen_weight  ) ;
       fillVariableWithValue( "Ele2_Eta"                      , Ele2_Eta    , gen_weight  ) ;
       fillVariableWithValue( "abs_Ele2_Eta"                  , fabs(Ele2_Eta), gen_weight  );
       fillVariableWithValue( "M_e1e2"                        , M_e1e2      , gen_weight  ) ;
       fillVariableWithValue( "Pt_e1e2"                       , Pt_e1e2     , gen_weight  ) ;
     }									    
									    
     // Jets								    
     fillVariableWithValue(   "nJet"                          , nJet_Ana    , gen_weight  ) ;
     if ( nJet_Ana >= 1 ) { 						    
       fillVariableWithValue( "Jet1_Pt"                       , Jet1_Pt     , gen_weight  ) ;
       fillVariableWithValue( "Jet1_Eta"                      , Jet1_Eta    , gen_weight  ) ;
     }
     if ( nJet_Ana >= 2 ) { 
       fillVariableWithValue( "Jet2_Pt"                       , Jet2_Pt     , gen_weight  ) ;
       fillVariableWithValue( "Jet2_Eta"                      , Jet2_Eta    , gen_weight  ) ;
       fillVariableWithValue( "Pt_j1j2"                       , Pt_j1j2     , gen_weight  ) ;
       fillVariableWithValue( "M_j1j2"                        , M_j1j2      , gen_weight  ) ;
       fillVariableWithValue( "DR_Jet1Jet2"                   , DR_Jet1Jet2 , gen_weight  ) ;
     }

     // Muons
     fillVariableWithValue(   "PassNMuon"                     , PassNMuon   , gen_weight  ) ;

     // DeltaR
     if ( nEle_Ana >= 2 && nJet_Ana >= 1) {
       fillVariableWithValue( "DR_Ele1Jet1"                   , DR_Ele1Jet1 , gen_weight  ) ;
       fillVariableWithValue( "DR_Ele2Jet1"                   , DR_Ele2Jet1 , gen_weight  ) ;
       if(nJet_Ana >= 2) {
	 fillVariableWithValue( "DR_Ele1Jet2"                 , DR_Ele1Jet2 , gen_weight  ) ;
	 fillVariableWithValue( "DR_Ele2Jet2"                 , DR_Ele2Jet2 , gen_weight  ) ;
       }
     }

     // sT
     if ( nEle_Ana >= 2 && nJet_Ana >= 2) {
       fillVariableWithValue( "sT_eejj"                      , sT_eejj, gen_weight  ) ;
     }      

     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Fill preselection plots
     //--------------------------------------------------------------------------

     FillUserTH1D( "PileupWeight"   , pileup_weight );
     FillUserTH1D( "GeneratorWeight", gen_weight ) ;
     
     bool passed_minimum      = ( passedAllPreviousCuts("PassBeamHaloFilterTight") && passedCut ("PassBeamHaloFilterTight"));

     if ( passed_minimum && isData ){ 
       FillUserTH1D ("run_HLT", run );
       profile_run_vs_nvtx_HLT -> Fill ( run, nVertex, 1 ) ;
     }

     bool passed_preselection = ( passedAllPreviousCuts("M_e1e2") && passedCut ("M_e1e2") );
     
     if ( passed_preselection ) {

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
       
       FillUserTH1D("nElectron_PAS"        , nEle_Ana         , pileup_weight * gen_weight ) ;
       FillUserTH1D("nMuon_PAS"            , nMuon_Stored     , pileup_weight * gen_weight ) ;
       FillUserTH1D("nJet_PAS"             , nJet_Ana         , pileup_weight * gen_weight ) ;
       FillUserTH1D("Pt1stEle_PAS"	   , Ele1_Pt          , pileup_weight * gen_weight ) ;
       FillUserTH1D("Eta1stEle_PAS"	   , Ele1_Eta         , pileup_weight * gen_weight ) ;
       FillUserTH1D("Phi1stEle_PAS"	   , Ele1_Phi         , pileup_weight * gen_weight ) ;
       FillUserTH1D("Pt2ndEle_PAS"	   , Ele2_Pt          , pileup_weight * gen_weight ) ;
       FillUserTH1D("Eta2ndEle_PAS"	   , Ele2_Eta         , pileup_weight * gen_weight ) ;
       FillUserTH1D("Phi2ndEle_PAS"	   , Ele2_Phi         , pileup_weight * gen_weight ) ;
       FillUserTH1D("Charge1stEle_PAS"	   , Ele1_Charge      , pileup_weight * gen_weight ) ;
       FillUserTH1D("Charge2ndEle_PAS"	   , Ele2_Charge      , pileup_weight * gen_weight ) ;
       FillUserTH1D("MET_PAS"              , MET_Pt           , pileup_weight * gen_weight ) ;
       FillUserTH1D("METSig_PAS"           , PFMETSig         , pileup_weight * gen_weight ) ;
       FillUserTH1D("METPhi_PAS"	   , MET_Phi          , pileup_weight * gen_weight ) ;
       FillUserTH1D("METCharged_PAS"       , PFMETCharged     , pileup_weight * gen_weight ) ;
       FillUserTH1D("METChargedPhi_PAS"    , PFMETChargedPhi  , pileup_weight * gen_weight ) ;   
       FillUserTH1D("METType1_PAS"         , PFMETType1Cor    , pileup_weight * gen_weight ) ;
       FillUserTH1D("METType1Phi_PAS"      , PFMETPhiType1Cor , pileup_weight * gen_weight ) ;   
       FillUserTH1D("Pt1stJet_PAS"         , Jet1_Pt          , pileup_weight * gen_weight ) ;
       FillUserTH1D("Pt2ndJet_PAS"         , Jet2_Pt          , pileup_weight * gen_weight ) ;
       FillUserTH1D("Eta1stJet_PAS"        , Jet1_Eta         , pileup_weight * gen_weight ) ;
       FillUserTH1D("Eta2ndJet_PAS"        , Jet2_Eta         , pileup_weight * gen_weight ) ;
       FillUserTH1D("Phi1stJet_PAS"	   , Jet1_Phi         , pileup_weight * gen_weight ) ;
       FillUserTH1D("Phi2ndJet_PAS"	   , Jet2_Phi         , pileup_weight * gen_weight ) ;
       FillUserTH1D("sTlep_PAS"            , Ele1_Pt + Ele2_Pt, pileup_weight * gen_weight ) ;
       FillUserTH1D("sTjet_PAS"            , Jet1_Pt + Jet2_Pt, pileup_weight * gen_weight ) ;
       FillUserTH1D("sT_PAS"               , sT_eejj          , pileup_weight * gen_weight ) ;
       FillUserTH1D("Mjj_PAS"		   , M_j1j2           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Mee_PAS"		   , M_e1e2           , pileup_weight * gen_weight ) ;
       FillUserTH1D( "MTenu_PAS"           , MT_Ele1MET       , pileup_weight * gen_weight ) ;
       FillUserTH1D("Me1j1_PAS"            , M_e1j1           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Me1j2_PAS"            , M_e1j2           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Me2j1_PAS"            , M_e2j1           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Me2j2_PAS"            , M_e2j2           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Ptee_PAS"             , Pt_e1e2          , pileup_weight * gen_weight ) ;
       FillUserTH1D("DCotTheta1stEle_PAS"  , Ele1_DCotTheta   , pileup_weight * gen_weight ) ;
       FillUserTH1D("Dist1stEle_PAS"       , Ele1_Dist        , pileup_weight * gen_weight ) ;
       FillUserTH1D("DCotTheta2ndEle_PAS"  , Ele2_DCotTheta   , pileup_weight * gen_weight ) ;
       FillUserTH1D("Dist2ndEle_PAS"       , Ele2_Dist        , pileup_weight * gen_weight ) ;
       FillUserTH1D("nVertex_PAS"          , nVertex          , pileup_weight * gen_weight ) ;
       FillUserTH1D("nVertex_good_PAS"     , nVertex_good     , pileup_weight * gen_weight ) ;
       FillUserTH1D("DR_Ele1Jet1_PAS"	   , DR_Ele1Jet1      , pileup_weight * gen_weight ) ;
       FillUserTH1D("DR_Ele1Jet2_PAS"	   , DR_Ele1Jet2      , pileup_weight * gen_weight ) ;
       FillUserTH1D("DR_Ele2Jet1_PAS"	   , DR_Ele2Jet1      , pileup_weight * gen_weight ) ;
       FillUserTH1D("DR_Ele2Jet2_PAS"	   , DR_Ele2Jet2      , pileup_weight * gen_weight ) ;
       FillUserTH1D("DR_Jet1Jet2_PAS"	   , DR_Jet1Jet2      , pileup_weight * gen_weight ) ;
       
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
       }
       
       double DR_Ele1Jet3 = 999.;
       double DR_Ele2Jet3 = 999.;
       double min_DR_EleJet = 999.;
       TLorentzVector eejj, e1e2mu, e1, j1, e2, j2,j3, mu, met;
       e1.SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
       e2.SetPtEtaPhiM ( Ele2_Pt, Ele2_Eta, Ele2_Phi, 0.0 );
       j1.SetPtEtaPhiM ( Jet1_Pt, Jet1_Eta, Jet1_Phi, 0.0 );
       j2.SetPtEtaPhiM ( Jet2_Pt, Jet2_Eta, Jet2_Phi, 0.0 );
       if ( nJet_Ana > 2 ) { 
	 j3.SetPtEtaPhiM ( Jet3_Pt, Jet3_Eta, Jet3_Phi, 0.0 );
	 DR_Ele1Jet3 = e1.DeltaR( j3 );
	 DR_Ele2Jet3 = e2.DeltaR( j3 );
       }
       
       if ( DR_Ele1Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet1;
       if ( DR_Ele1Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet2;
       if ( DR_Ele2Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet1;
       if ( DR_Ele2Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet2;
       if ( nJet_Ana > 2 ) {
	 if ( DR_Ele1Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet3;
	 if ( DR_Ele2Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet3;
       }
       
       FillUserTH1D( "minDR_EleJet_PAS", min_DR_EleJet, pileup_weight * gen_weight);
       
       
       mu.SetPtEtaPhiM ( Muon1_Pt, Muon1_Eta, Muon1_Phi, 0.0 );
       met.SetPtEtaPhiM ( MET_Pt, 0.0, MET_Phi, 0.0 );
       
       double DR_Ele1Ele2 = e1.DeltaR( e2 ) ;
       FillUserTH1D("DR_Ele1Ele2_PAS"	   , DR_Ele1Ele2      , pileup_weight * gen_weight ) ;
       
       eejj = e1 + e2 + j1 + j2 ; 
       double M_eejj = eejj.M();

       // e1e2mu = e1 + e2 + mu;
       // double MT_eemuMET = sqrt(2 * e1e2mu.Pt()    * MET_Pt  * (1 - cos(e1e2mu.DeltaPhi (met))));
       // FillUserTH1D("MTeemunu_PAS"         , MT_eemuMET                        , pileup_weight * gen_weight );
	      
       FillUserTH1D("Meejj_PAS", M_eejj , pileup_weight * gen_weight );
       
       if ( fabs(M_e1j1-M_e2j2) < fabs(M_e1j2-M_e2j1) )  {
	 
	 double M_ej_avg = (M_e1j1 + M_e2j2) / 2.0;

	 FillUserTH1D("Mej_selected_avg_PAS", M_ej_avg, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Me1j_selected_PAS"   , M_e1j1  , pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Me2j_selected_PAS"   , M_e2j2  , pileup_weight * gen_weight) ;	   
	 FillUserTH2D( "Me1jVsMe2j_selected", M_e1j1  , M_e2j2, pileup_weight * gen_weight ) ;
	 FillUserTH2D( "Me1jVsMe2j_rejected", M_e1j2  , M_e2j1, pileup_weight * gen_weight ) ;
       }

       else {

	 double M_ej_avg = (M_e1j2 + M_e2j1) / 2.0;

	 FillUserTH1D("Mej_selected_avg_PAS", M_ej_avg, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Me1j_selected_PAS"   , M_e1j2, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Me2j_selected_PAS"   , M_e2j1, pileup_weight * gen_weight) ;	   
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

