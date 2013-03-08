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
#include <TRandom3.h>

#include "Lumi3DReWeighting.h"

#include <algorithm>


void getJetScaleFactor ( float eta, float & mean, float &error ) { 
  
  float fabs_eta = fabs ( eta ) ;
  
  if      ( fabs_eta > 0.0 && fabs_eta <= 0.5 ) {
    mean  = 1.052;
    error = 0.063;
  }
  else if ( fabs_eta > 0.5 && fabs_eta <= 1.1 ) {
    mean  = 1.057;
    error = 0.057;
  }
  else if ( fabs_eta > 1.1 && fabs_eta <= 1.7 ) {
    mean  = 1.096;
    error = 0.065;
  }
  else if ( fabs_eta > 1.7 && fabs_eta <= 2.3 ) {
    mean  = 1.134;
    error = 0.094;
  }
  else {
    mean  = 1.288;
    error = 0.200;
  }
}

void getEleScaleFactor ( float eta, float & mean, float &error ) { 
  
  float fabs_eta = fabs ( eta ) ;
  
  if      ( fabs_eta > 0.0 && fabs_eta <= 1.45 ) {
    mean  = 1.004;
    error = 0.000;
  }
  if      ( fabs_eta > 1.55 && fabs_eta <= 3.0 ) {
    mean  = 1.041;
    error = 0.000;
  }
}

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile){}

analysisClass::~analysisClass(){}

void analysisClass::Loop()
{
   std::cout << "analysisClass::Loop() begins" <<std::endl;   

   //--------------------------------------------------------------------------
   // Final selection mass points
   //--------------------------------------------------------------------------

   const int n_lq_mass = 19;

   int LQ_MASS[n_lq_mass] = { 
     300,  350,  400, 450, 500, 550,  600,  650,
     700,  750,  800, 850, 900, 950, 1000, 1050,
     1100, 1150, 1200
   };

   std::vector<bool> passed_vector;
   
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

   // scaling

   double jet_energy_scale_sign     = (double) getPreCutValue1("jes_scale_sign") ;
   double ele_bar_energy_scale      = (double) getPreCutValue1("ees_bar_scale") ;
   double ele_end_energy_scale      = (double) getPreCutValue1("ees_end_scale") ;
   
   // jer / eer studies
   bool do_eer = ( (int) getPreCutValue1("do_eer") == 1 );
   bool do_jer = ( (int) getPreCutValue1("do_jer") == 1 );
   
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
   CreateUserTH1D( "Mjj_PAS"		   ,    200 , 0       , 2000	 ); 
   CreateUserTH1D( "Mee_PAS"		   ,    200 , 0       , 2000	 ); 
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
		                           
   CreateUserTH1D( "DR_Ele1Jet1_PAS"	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele1Jet2_PAS"	   , 	getHistoNBins("DR_Ele1Jet2"), getHistoMin("DR_Ele1Jet2"), getHistoMax("DR_Ele1Jet2")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet1_PAS"	   , 	getHistoNBins("DR_Ele2Jet1"), getHistoMin("DR_Ele2Jet1"), getHistoMax("DR_Ele2Jet1")     ) ; 
   CreateUserTH1D( "DR_Ele2Jet2_PAS"	   , 	getHistoNBins("DR_Ele2Jet2"), getHistoMin("DR_Ele2Jet2"), getHistoMax("DR_Ele2Jet2")     ) ; 
   CreateUserTH1D( "DR_Jet1Jet2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
   CreateUserTH1D( "DR_Ele1Ele2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
   CreateUserTH1D( "minDR_EleJet_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 

   CreateUserTH2D( "Me1jVsMe2j_selected",     200, 0, 2000, 200, 0, 2000) ;
   CreateUserTH2D( "Me1jVsMe2j_rejected",     200, 0, 2000, 200, 0, 2000) ;

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

   //--------------------------------------------------------------------------
   // Final selection plots
   //--------------------------------------------------------------------------

   char plot_name[100];

   for (int i_lq_mass = 0; i_lq_mass < n_lq_mass ; ++i_lq_mass ) { 
     int lq_mass = LQ_MASS[i_lq_mass];
     sprintf(plot_name, "Mej_selected_avg_LQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 50  , 0 , 2500 );
     sprintf(plot_name, "Mej_selected_min_LQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 50  , 0 , 2500 );
     sprintf(plot_name, "Mej_selected_max_LQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 50  , 0 , 2500 );
     sprintf(plot_name, "Mej_minmax_LQ%d"             , lq_mass ); CreateUserTH1D ( plot_name, 50  , 0 , 2500 );
     sprintf(plot_name, "sT_eejj_LQ%d"                , lq_mass ); CreateUserTH1D ( plot_name, 25  , 0 , 2500 );
     sprintf(plot_name, "Mee_LQ%d"                    , lq_mass ); CreateUserTH1D ( plot_name, 40  , 0 , 2000 );
     sprintf(plot_name, "Mej_selected_min_vs_max_LQ%d", lq_mass ); CreateUserTH2D ( plot_name, 50  , 0 , 1000, 50  , 0 , 1000 );
   }

   CreateUserTH1D( "PileupWeight"   , 100, -10, 10 );
   CreateUserTH1D( "GeneratorWeight", 100, -2.0, 2.0);
   
   //--------------------------------------------------------------------------
   // Loop over the chain
   //--------------------------------------------------------------------------

   if (fChain == 0) return;

   unsigned int seed = 987654321;
   TRandom3 * rootEngine = new TRandom3 ( seed ) ;
   
   Long64_t nentries = fChain->GetEntries();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   double jer_recojetpt_threshold = 10.0;
   double jer_genjetpt_threshold = 2.0;
   double jet_match_dr_threshold = 0.4;

   double eer_recoelept_threshold = 10.0;
   double eer_genelept_threshold = 2.0;   
   double ele_match_dr_threshold = 10.0;


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

     double pileup_weight = getPileupWeight ( nPileUpInt_True, isData ) ;
     
     if ( !isData ) { 

       //--------------------------------------------------------------------------
       // Identify barrel / endcap electron
       //--------------------------------------------------------------------------
       
       bool ele1_isBarrel  = false;
       bool ele1_isEndcap  = false;
       bool ele2_isBarrel  = false;
       bool ele2_isEndcap  = false;

       if( fabs( Ele1_Eta  ) < eleEta_bar_max )    ele1_isBarrel = true;
       if( fabs( Ele1_Eta  ) > eleEta_end_min &&
	   fabs( Ele1_Eta  ) < eleEta_end_max )    ele1_isEndcap = true;


       if( fabs( Ele2_Eta  ) < eleEta_bar_max )    ele2_isBarrel = true;
       if( fabs( Ele2_Eta  ) > eleEta_end_min &&
	   fabs( Ele2_Eta  ) < eleEta_end_max )    ele2_isEndcap = true;

       double ele1_energy_scale;
       double ele2_energy_scale;
       

       if( ele1_isBarrel ) {
	 ele1_energy_scale      = ele_bar_energy_scale;
       }		    
       else { 
	 ele1_energy_scale      = ele_end_energy_scale;
       }
       
       if( ele2_isBarrel ) {
	 ele2_energy_scale      = ele_bar_energy_scale;
       }		    
       else { 
	 ele2_energy_scale      = ele_end_energy_scale;
       }
       
       //--------------------------------------------------------------------------
       // JER / EER matching
       //--------------------------------------------------------------------------

       std::vector < float > v_genJet_Pt ;
       std::vector < float > v_genJet_Eta;
       std::vector < float > v_genJet_Phi;
       
       std::vector < float > v_recoJet_Pt ;
       std::vector < float > v_recoJet_Eta;
       std::vector < float > v_recoJet_Phi;
       
       std::vector < int > v_recoJet_index;
       std::vector < int > v_genJet_index;
       std::vector < int > v_recoJet_matchedGenIndex ( 5, -1 );
       std::vector < float > v_recoJet_bestDR        ( 5, 9999.0);
       
       
       std::vector < float > v_genEle_Pt;
       std::vector < float > v_genEle_Eta;
       std::vector < float > v_genEle_Phi;
       
       std::vector < float > v_recoEle_Pt;
       std::vector < float > v_recoEle_Eta;
       std::vector < float > v_recoEle_Phi;
       
       std::vector < int > v_recoEle_index;
       std::vector < int > v_genEle_index;
       std::vector < int > v_recoEle_matchedGenIndex ( 2, -1 );
       std::vector < float > v_recoEle_bestDR        ( 2, 9999.0 );

       int n_reco_jets ;
       int n_gen_jets  ;
       int n_reco_eles ;
       int n_gen_eles  ;

       
       if ( do_jer ) {
	 
	 // reco jetsl
	 
	 if ( Jet1_Pt > jer_recojetpt_threshold ) {
	   v_recoJet_Pt   .push_back ( Jet1_Pt  ) ;
	   v_recoJet_Eta  .push_back ( Jet1_Eta ) ;
	   v_recoJet_Phi  .push_back ( Jet1_Phi ) ;
	   v_recoJet_index.push_back ( 1 ) ;
	 }
	 
	 if ( Jet2_Pt > jer_recojetpt_threshold ) {
	   v_recoJet_Pt   .push_back ( Jet2_Pt  ) ;
	   v_recoJet_Eta  .push_back ( Jet2_Eta ) ;
	   v_recoJet_Phi  .push_back ( Jet2_Phi ) ;
	   v_recoJet_index.push_back ( 2 ) ;
	 }
	 
	 if ( Jet3_Pt > jer_recojetpt_threshold ) {
	   v_recoJet_Pt   .push_back ( Jet3_Pt  ) ;
	   v_recoJet_Eta  .push_back ( Jet3_Eta ) ;
	   v_recoJet_Phi  .push_back ( Jet3_Phi ) ;
	   v_recoJet_index.push_back ( 3 ) ;
	 }
	 
	 if ( Jet4_Pt > jer_recojetpt_threshold ) {
	   v_recoJet_Pt   .push_back ( Jet4_Pt  ) ;
	   v_recoJet_Eta  .push_back ( Jet4_Eta ) ;
	   v_recoJet_Phi  .push_back ( Jet4_Phi ) ;
	   v_recoJet_index.push_back ( 4 ) ;
	 }
	 
	 if ( Jet5_Pt > jer_recojetpt_threshold ) {
	   v_recoJet_Pt   .push_back ( Jet5_Pt  ) ;
	   v_recoJet_Eta  .push_back ( Jet5_Eta ) ;
	   v_recoJet_Phi  .push_back ( Jet5_Phi ) ;
	   v_recoJet_index.push_back ( 5 ) ;
	 }
	 
	 // gen jets 
	 
	 if ( GenJet1_Pt > jer_genjetpt_threshold ) {
	   v_genJet_Pt   .push_back ( GenJet1_Pt  ) ;
	   v_genJet_Eta  .push_back ( GenJet1_Eta ) ;
	   v_genJet_Phi  .push_back ( GenJet1_Phi ) ;
	   v_genJet_index.push_back ( 1 ) ;
	 }
	 
	 if ( GenJet2_Pt > jer_genjetpt_threshold ) {
	   v_genJet_Pt   .push_back ( GenJet2_Pt  ) ;
	   v_genJet_Eta  .push_back ( GenJet2_Eta ) ;
	   v_genJet_Phi  .push_back ( GenJet2_Phi ) ;
	   v_genJet_index.push_back ( 2 ) ;
	 }
	 
	 if ( GenJet3_Pt > jer_genjetpt_threshold ) {
	   v_genJet_Pt   .push_back ( GenJet3_Pt  ) ;
	   v_genJet_Eta  .push_back ( GenJet3_Eta ) ;
	   v_genJet_Phi  .push_back ( GenJet3_Phi ) ;
	   v_genJet_index.push_back ( 3 ) ;
	 }
	 
	 if ( GenJet4_Pt > jer_genjetpt_threshold ) {
	   v_genJet_Pt   .push_back ( GenJet4_Pt  ) ;
	   v_genJet_Eta  .push_back ( GenJet4_Eta ) ;
	   v_genJet_Phi  .push_back ( GenJet4_Phi ) ;
	   v_genJet_index.push_back ( 4 ) ;
	 }
	 
	 if ( GenJet5_Pt > jer_genjetpt_threshold ) {
	   v_genJet_Pt   .push_back ( GenJet5_Pt  ) ;
	   v_genJet_Eta  .push_back ( GenJet5_Eta ) ;
	   v_genJet_Phi  .push_back ( GenJet5_Phi ) ;
	   v_genJet_index.push_back ( 5 ) ;
	 }

	 n_reco_jets = v_recoJet_index.size() ; 
	 n_gen_jets  = v_genJet_index.size() ; 
	 
       }

       if ( do_eer ) { 
	 
	 // reco ele
	 
	 if ( Ele1_Pt > eer_recoelept_threshold ) {
	   v_recoEle_Pt   .push_back ( Ele1_Pt  ) ;
	   v_recoEle_Eta  .push_back ( Ele1_Eta ) ;
	   v_recoEle_Phi  .push_back ( Ele1_Phi ) ;
	   v_recoEle_index.push_back ( 1 ) ;
	 }
	 
	 if ( Ele2_Pt > eer_recoelept_threshold ) {
	   v_recoEle_Pt   .push_back ( Ele2_Pt  ) ;
	   v_recoEle_Eta  .push_back ( Ele2_Eta ) ;
	   v_recoEle_Phi  .push_back ( Ele2_Phi ) ;
	   v_recoEle_index.push_back ( 2 ) ;
	 }
	 
	 // gen ele
	 
	 if ( GenEle1_Pt > eer_genelept_threshold ) {
	   v_genEle_Pt   .push_back ( GenEle1_Pt  ) ;
	   v_genEle_Eta  .push_back ( GenEle1_Eta ) ;
	   v_genEle_Phi  .push_back ( GenEle1_Phi ) ;
	   v_genEle_index.push_back ( 1 ) ;
	 }
	 
	 if ( GenEle2_Pt > eer_genelept_threshold ) {
	   v_genEle_Pt   .push_back ( GenEle2_Pt  ) ;
	   v_genEle_Eta  .push_back ( GenEle2_Eta ) ;
	   v_genEle_Phi  .push_back ( GenEle2_Phi ) ;
	   v_genEle_index.push_back ( 2 ) ;
	 }
	 
	 n_reco_eles = v_recoEle_index.size() ; 
	 n_gen_eles = v_genEle_index.size() ; 


       }

       if ( do_jer ) { 
	 // jet matching
	 
	 TLorentzVector gen_jet, reco_jet;
	 
	 for (int i_reco_jet = 0; i_reco_jet < n_reco_jets; ++i_reco_jet){
	   for (int i_gen_jet = 0; i_gen_jet < n_gen_jets; ++i_gen_jet){
	     
	     gen_jet .SetPtEtaPhiM ( v_genJet_Pt [i_gen_jet  ], v_genJet_Eta [i_gen_jet ], v_genJet_Phi [i_gen_jet ], 0.0);
	     reco_jet.SetPtEtaPhiM ( v_recoJet_Pt[i_reco_jet ], v_recoJet_Eta[i_reco_jet], v_recoJet_Phi[i_reco_jet], 0.0);
	     
	     float gen_reco_dR = gen_jet.DeltaR( reco_jet ) ;
	     
	     if ( gen_reco_dR > jet_match_dr_threshold ) continue;
	     if ( gen_reco_dR > v_recoJet_bestDR [ i_reco_jet ] ) continue;
	     
	     v_recoJet_bestDR [ i_reco_jet ] = gen_reco_dR ;
	     v_recoJet_matchedGenIndex [i_reco_jet] = i_gen_jet;
	     
	   }
	 }

	 int n_jet_matched = 0;
	 
	 std::vector < bool >  jet_matched ( 6, false );

	 for (int i_reco_jet = 0; i_reco_jet < n_reco_jets; ++i_reco_jet){
	   
	   if ( v_recoJet_matchedGenIndex [i_reco_jet] < 0 ) continue;
	   
	   int index = v_recoJet_index [ i_reco_jet ];
	   
	   if ( jet_matched [ index ] ){
	     std::cout << "Reco jet #" << index << " was matched twice!" << std::endl;
	     exit(0);
	   }
	   
	   jet_matched [ index ] = true;
	   
	   // h_jet_dr -> Fill ( v_recoJet_bestDR [ i_reco_jet ] ) ;
	   
	   // n_jet_matched++;
	 }
	 // h_jet_nMatched -> Fill ( n_jet_matched ) ;

       }
       

       if ( do_eer ) { 
       
	 // ele matching
	 
	 TLorentzVector gen_ele, reco_ele;
	 
	 for (int i_reco_ele = 0; i_reco_ele < n_reco_eles; ++i_reco_ele){
	   for (int i_gen_ele = 0; i_gen_ele < n_gen_eles; ++i_gen_ele){
	     
	     gen_ele .SetPtEtaPhiM ( v_genEle_Pt [i_gen_ele  ], v_genEle_Eta [i_gen_ele ], v_genEle_Phi [i_gen_ele ], 0.0);
	     reco_ele.SetPtEtaPhiM ( v_recoEle_Pt[i_reco_ele ], v_recoEle_Eta[i_reco_ele], v_recoEle_Phi[i_reco_ele], 0.0);
	     
	     float gen_reco_dR      = gen_ele.DeltaR( reco_ele ) ;
	     float gen_reco_ptratio = ( reco_ele.Pt() - gen_ele.Pt() ) / gen_ele.Pt() ;
	     
	     if ( gen_reco_dR > ele_match_dr_threshold ) continue;
	     if ( gen_reco_dR > v_recoEle_bestDR [ i_reco_ele ] ) continue;
	     
	     v_recoEle_bestDR [ i_reco_ele ] = gen_reco_dR ;
	     v_recoEle_matchedGenIndex [i_reco_ele] = i_gen_ele;
	     
	   }
	 }
       
	 int n_ele_matched = 0;

	 std::vector < bool >  ele_matched ( 4, false );
	 
	 for (int i_reco_ele = 0; i_reco_ele < n_reco_eles; ++i_reco_ele){
	   
	   if ( v_recoEle_matchedGenIndex [i_reco_ele] < 0 ) continue;
	   
	   int index = v_recoEle_index [ i_reco_ele ];
	   
	   if ( ele_matched [ index ] ){
	     std::cout << "Reco ele #" << index << " was matched twice!" << std::endl;
	     exit(0);
	   }
	   
	   ele_matched [ index ] = true;
	   
	   // h_ele_dr -> Fill ( v_recoEle_bestDR [ i_reco_ele ] ) ;
	   
	   // n_ele_matched++;
	 }
	 
	 // h_ele_nMatched -> Fill ( n_ele_matched ) ;
       }
       
       //--------------------------------------------------------------------------
       // JES / JER
       //--------------------------------------------------------------------------
       
       TLorentzVector v_Jet1_old, v_Jet2_old, v_Jet3_old, v_Jet4_old, v_Jet5_old, v_MET_old;
       v_Jet1_old.SetPtEtaPhiM ( Jet1_Pt, Jet1_Eta, Jet1_Phi, 0.0 );
       v_Jet2_old.SetPtEtaPhiM ( Jet2_Pt, Jet2_Eta, Jet2_Phi, 0.0 );
       v_Jet3_old.SetPtEtaPhiM ( Jet3_Pt, Jet3_Eta, Jet3_Phi, 0.0 );
       v_Jet4_old.SetPtEtaPhiM ( Jet4_Pt, Jet4_Eta, Jet4_Phi, 0.0 );
       v_Jet5_old.SetPtEtaPhiM ( Jet5_Pt, Jet5_Eta, Jet5_Phi, 0.0 );
       v_MET_old.SetPtEtaPhiM  ( PFMET_Type01XY_Pt,  0.0     , PFMET_Type01XY_Phi , 0.0 );
       
       TLorentzVector v_Jet1_new, v_Jet2_new, v_Jet3_new, v_Jet4_new, v_Jet5_new, v_MET_new;       
       v_Jet1_new.SetPtEtaPhiM ( Jet1_Pt * ( 1.0 + ( jet_energy_scale_sign * Jet1_JECUnc ) ), Jet1_Eta, Jet1_Phi, 0.0 );
       v_Jet2_new.SetPtEtaPhiM ( Jet2_Pt * ( 1.0 + ( jet_energy_scale_sign * Jet2_JECUnc ) ), Jet2_Eta, Jet2_Phi, 0.0 );
       v_Jet3_new.SetPtEtaPhiM ( Jet3_Pt * ( 1.0 + ( jet_energy_scale_sign * Jet3_JECUnc ) ), Jet3_Eta, Jet3_Phi, 0.0 );
       v_Jet4_new.SetPtEtaPhiM ( Jet4_Pt * ( 1.0 + ( jet_energy_scale_sign * Jet4_JECUnc ) ), Jet4_Eta, Jet4_Phi, 0.0 );
       v_Jet5_new.SetPtEtaPhiM ( Jet5_Pt * ( 1.0 + ( jet_energy_scale_sign * Jet5_JECUnc ) ), Jet5_Eta, Jet5_Phi, 0.0 );
       
       if ( do_jer ) { 
    
	 std::vector < float > v_recoJet_newPt;
	 v_recoJet_newPt.push_back ( v_Jet1_new.Pt() );
	 v_recoJet_newPt.push_back ( v_Jet2_new.Pt() );
	 v_recoJet_newPt.push_back ( v_Jet3_new.Pt() );
	 v_recoJet_newPt.push_back ( v_Jet4_new.Pt() );
	 v_recoJet_newPt.push_back ( v_Jet5_new.Pt() );
	 
	 for (int i_reco_jet = 0; i_reco_jet < n_reco_jets; ++i_reco_jet){
	   
	   v_recoJet_newPt [ i_reco_jet ] = v_recoJet_Pt [ i_reco_jet ];
	   if ( v_recoJet_matchedGenIndex [i_reco_jet] < 0 ) continue;
	   
	   float gen_jet_pt   = v_genJet_Pt   [ v_recoJet_matchedGenIndex [i_reco_jet] ];
	   float reco_jet_pt  = v_recoJet_Pt  [ i_reco_jet ];
	   float reco_jet_eta = v_recoJet_Eta [ i_reco_jet ];
	   
	   float scale_factor, scale_factor_error;
	   getJetScaleFactor ( reco_jet_eta, scale_factor, scale_factor_error ) ;
	   
	   float smeared_scale_factor = scale_factor * rootEngine -> Gaus ( 1.0, scale_factor_error ) ;
	   float delta_pt = smeared_scale_factor * ( reco_jet_pt - gen_jet_pt ) ;
	   float new_jet_pt = max ( float (0.0), gen_jet_pt + delta_pt ) ;
	   
	   v_recoJet_newPt [ i_reco_jet ] = new_jet_pt ;

	 }
	 
	 
	 v_Jet1_new.SetPtEtaPhiM ( v_recoJet_newPt[0] , Jet1_Eta, Jet1_Phi, 0.0 );
	 v_Jet2_new.SetPtEtaPhiM ( v_recoJet_newPt[1] , Jet2_Eta, Jet2_Phi, 0.0 );
	 v_Jet3_new.SetPtEtaPhiM ( v_recoJet_newPt[2] , Jet3_Eta, Jet3_Phi, 0.0 );
	 v_Jet4_new.SetPtEtaPhiM ( v_recoJet_newPt[3] , Jet4_Eta, Jet4_Phi, 0.0 );
	 v_Jet5_new.SetPtEtaPhiM ( v_recoJet_newPt[4] , Jet5_Eta, Jet5_Phi, 0.0 );
	 
       }
       
       
       TLorentzVector v_Jet1_delta, v_Jet2_delta, v_Jet3_delta, v_Jet4_delta, v_Jet5_delta;
       v_Jet1_delta = v_Jet1_old - v_Jet1_new;
       v_Jet2_delta = v_Jet2_old - v_Jet2_new;
       v_Jet3_delta = v_Jet3_old - v_Jet3_new;
       v_Jet4_delta = v_Jet4_old - v_Jet4_new;
       v_Jet5_delta = v_Jet5_old - v_Jet5_new;

       if ( do_jer ){ 
	 if ( v_recoJet_matchedGenIndex [0] < 0  && v_Jet1_delta.Pt() > 0.00001 ||
	      v_recoJet_matchedGenIndex [1] < 0  && v_Jet2_delta.Pt() > 0.00001 ||
	      v_recoJet_matchedGenIndex [2] < 0  && v_Jet3_delta.Pt() > 0.00001 ||
	      v_recoJet_matchedGenIndex [3] < 0  && v_Jet4_delta.Pt() > 0.00001 ||
	      v_recoJet_matchedGenIndex [4] < 0  && v_Jet5_delta.Pt() > 0.00001 ){
	   std::cout << "ERROR: unmatched jet has changed pt!" << std::endl;
	   std::cout << "  v_recoJet_matchedGenIndex [0] = " << v_recoJet_matchedGenIndex [0] << std::endl;
	   std::cout << "  v_recoJet_matchedGenIndex [1] = " << v_recoJet_matchedGenIndex [1] << std::endl;
	   std::cout << "  v_recoJet_matchedGenIndex [2] = " << v_recoJet_matchedGenIndex [2] << std::endl;
	   std::cout << "  v_recoJet_matchedGenIndex [3] = " << v_recoJet_matchedGenIndex [3] << std::endl;
	   std::cout << "  v_recoJet_matchedGenIndex [4] = " << v_recoJet_matchedGenIndex [4] << std::endl;
	   std::cout << "  v_Jet1_delta.Pt() = " << v_Jet1_delta.Pt() << std::endl;
	   std::cout << "  v_Jet2_delta.Pt() = " << v_Jet2_delta.Pt() << std::endl;
	   std::cout << "  v_Jet3_delta.Pt() = " << v_Jet3_delta.Pt() << std::endl;
	   std::cout << "  v_Jet4_delta.Pt() = " << v_Jet4_delta.Pt() << std::endl;
	   std::cout << "  v_Jet5_delta.Pt() = " << v_Jet5_delta.Pt() << std::endl;
	   exit(0);
	 }
       }
       
       Jet1_Pt = v_Jet1_new.Pt();
       Jet2_Pt = v_Jet2_new.Pt();
       Jet3_Pt = v_Jet3_new.Pt();
       Jet4_Pt = v_Jet4_new.Pt();
       Jet5_Pt = v_Jet5_new.Pt();
       
       //--------------------------------------------------------------------------
       // Electron ES / Electron ER
       //--------------------------------------------------------------------------
       
       TLorentzVector v_Ele1_old, v_Ele2_old;
       v_Ele1_old.SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
       v_Ele2_old.SetPtEtaPhiM ( Ele2_Pt, Ele2_Eta, Ele2_Phi, 0.0 );
       
       TLorentzVector v_Ele1_new, v_Ele2_new;
       v_Ele1_new.SetPtEtaPhiM ( Ele1_Pt * ele1_energy_scale , Ele1_Eta, Ele1_Phi, 0.0 );
       v_Ele2_new.SetPtEtaPhiM ( Ele2_Pt * ele2_energy_scale , Ele2_Eta, Ele2_Phi, 0.0 );
       

       if ( do_eer ) { 
    
	 std::vector < float > v_recoEle_newPt;
	 v_recoEle_newPt.push_back ( v_Ele1_new.Pt() );
	 v_recoEle_newPt.push_back ( v_Ele2_new.Pt() );
	 
	 for (int i_reco_ele = 0; i_reco_ele < n_reco_eles; ++i_reco_ele){
	   
	   v_recoEle_newPt [ i_reco_ele ] = v_recoEle_Pt [ i_reco_ele ];
	   if ( v_recoEle_matchedGenIndex [i_reco_ele] < 0 ) continue;

	   float gen_ele_pt   = v_genEle_Pt   [ v_recoEle_matchedGenIndex [i_reco_ele] ];
	   float reco_ele_pt  = v_recoEle_Pt  [ i_reco_ele ] ;
	   float reco_ele_eta = v_recoEle_Eta [ i_reco_ele ] ;
	   
	   float scale_factor, scale_factor_error;
	   getEleScaleFactor ( reco_ele_eta, scale_factor, scale_factor_error ) ;
	   
	   float smeared_scale_factor = scale_factor * rootEngine -> Gaus ( 1.0, scale_factor_error ) ;	   
	   float delta_pt = smeared_scale_factor * ( reco_ele_pt - gen_ele_pt ) ;
	   float new_ele_pt = max ( float (0.0), gen_ele_pt + delta_pt ) ;

	   v_recoEle_newPt [ i_reco_ele ] = new_ele_pt ;

	 }
	 
	 
	 v_Ele1_new.SetPtEtaPhiM ( v_recoEle_newPt[0] , Ele1_Eta, Ele1_Phi, 0.0 );
	 v_Ele2_new.SetPtEtaPhiM ( v_recoEle_newPt[1] , Ele2_Eta, Ele2_Phi, 0.0 );
	 
       }

       TLorentzVector v_Ele1_delta, v_Ele2_delta;
       v_Ele1_delta = v_Ele1_old - v_Ele1_new;
       v_Ele2_delta = v_Ele2_old - v_Ele2_new;
       
       if ( do_eer ) { 
	 if ( v_recoEle_matchedGenIndex [0] < 0  && v_Ele1_delta.Pt() > 0.0001 ||
	      v_recoEle_matchedGenIndex [1] < 0  && v_Ele2_delta.Pt() > 0.0001 ){
	   exit(0);
	 }
       }
       
       Ele1_Pt = v_Ele1_new.Pt();
       Ele2_Pt = v_Ele2_new.Pt();

       // if ( do_eer ) {
       // 	 std::cout << "-----------------------------------------------" << std::endl;
       // 	 std::cout << "Ele 1 matched = " << v_recoEle_matchedGenIndex [0] << ", delta pt = " << v_Ele1_delta.Pt() << ", Old pt = " << v_Ele1_old.Pt() << ", New pt = " << v_Ele1_new.Pt() << std::endl;
       // 	 std::cout << "Ele 2 matched = " << v_recoEle_matchedGenIndex [1] << ", delta pt = " << v_Ele2_delta.Pt() << ", Old pt = " << v_Ele2_old.Pt() << ", New pt = " << v_Ele2_new.Pt() << std::endl;
       // 
       // }
       // 
       // 
       // if ( do_jer ) { 
       // 	 std::cout << "-----------------------------------------------" << std::endl;
       // 	 std::cout << "Jet 1 matched = " << v_recoJet_matchedGenIndex [0] << ", delta pt = " << v_Jet1_delta.Pt() << ", Old pt = " << v_Jet1_old.Pt() << ", New pt = " << v_Jet1_new.Pt() << std::endl;
       // 	 std::cout << "Jet 2 matched = " << v_recoJet_matchedGenIndex [1] << ", delta pt = " << v_Jet2_delta.Pt() << ", Old pt = " << v_Jet2_old.Pt() << ", New pt = " << v_Jet2_new.Pt() << std::endl;
       // 	 std::cout << "Jet 3 matched = " << v_recoJet_matchedGenIndex [2] << ", delta pt = " << v_Jet3_delta.Pt() << ", Old pt = " << v_Jet3_old.Pt() << ", New pt = " << v_Jet3_new.Pt() << std::endl;
       // 	 std::cout << "Jet 4 matched = " << v_recoJet_matchedGenIndex [3] << ", delta pt = " << v_Jet4_delta.Pt() << ", Old pt = " << v_Jet4_old.Pt() << ", New pt = " << v_Jet4_new.Pt() << std::endl;
       // 	 std::cout << "Jet 5 matched = " << v_recoJet_matchedGenIndex [4] << ", delta pt = " << v_Jet5_delta.Pt() << ", Old pt = " << v_Jet5_old.Pt() << ", New pt = " << v_Jet5_new.Pt() << std::endl;
       // }
       
       if ( Ele1_Pt == 0.0 ) {
	 std::cout << "ERROR! Ele1_Pt shifted to 0.0!!" << std::endl;
	 exit(0);
	   
       }

       //std::cout << "-----------------------------------------------" << std::endl;
       //std::cout << "Ele 1 matched = " << v_recoEle_matchedGenIndex [0] << ", delta pt = " << v_Ele1_delta.Pt() << ", Old pt = " << v_Ele1_old.Pt() << ", New pt = " << v_Ele1_new.Pt() << std::endl;
       //std::cout << "Ele 2 matched = " << v_recoEle_matchedGenIndex [1] << ", delta pt = " << v_Ele2_delta.Pt() << ", Old pt = " << v_Ele2_old.Pt() << ", New pt = " << v_Ele2_new.Pt() << std::endl;


       //--------------------------------------------------------------------------
       // Carry JES / JER / EES / EER over to MET
       //--------------------------------------------------------------------------
       
       v_MET_new = v_MET_old + v_Jet1_delta + v_Jet2_delta + v_Jet3_delta + v_Jet4_delta + v_Jet5_delta;  // Jet correction
       v_MET_new = v_MET_new + v_Ele1_delta + v_Ele2_delta;                                               // Ele correction
       
       PFMET_Type01XY_Pt  = v_MET_new.Pt() ;
       PFMET_Type01XY_Phi = v_MET_new.Phi();
       
       //--------------------------------------------------------------------------
       // Re-calculate variabled dependent on Electron Pt, Jet Pt, MET, or MET phi
       //--------------------------------------------------------------------------
       
       TLorentzVector v_j1j2 = v_Jet1_new + v_Jet2_new;
       TLorentzVector v_e1e2 = v_Ele1_new + v_Ele2_new;
       TLorentzVector v_e1j1 = v_Ele1_new + v_Jet1_new;
       TLorentzVector v_e1j2 = v_Ele1_new + v_Jet2_new;
       TLorentzVector v_e2j1 = v_Ele2_new + v_Jet1_new;
       TLorentzVector v_e2j2 = v_Ele2_new + v_Jet2_new;
       
       M_e1j1 = v_e1j1.M();
       M_e1j2 = v_e1j2.M();
       M_e2j1 = v_e2j1.M();
       M_e2j2 = v_e2j2.M();
       M_e1e2 = v_e1e2.M();
       M_j1j2 = v_j1j2.M();
       
       mDPhi_METEle1 = fabs ( v_MET_new.DeltaPhi ( v_Ele1_new ) );
       mDPhi_METJet1 = fabs ( v_MET_new.DeltaPhi ( v_Jet1_new ) );
       
       sT_eejj = Ele1_Pt + Ele2_Pt + Jet1_Pt + Jet2_Pt;
       
       if ( nEle_store > 0 ){ 
	 TVector2 v_ele;
	 TVector2 v_MET;
	 v_ele.SetMagPhi ( Ele1_Pt, Ele1_Phi );
	 v_MET.SetMagPhi ( PFMET_Type01XY_Pt , PFMET_Type01XY_Phi  );
	 float deltaphi = v_ele.DeltaPhi ( v_MET );
	 MT_Ele1MET = sqrt ( 2 * PFMET_Type01XY_Pt * Ele1_Pt * ( 1 - cos ( deltaphi ) ) );
       }
     }

     //--------------------------------------------------------------------------
     // Get information about gen-level reweighting (should be for Sherpa only)
     //--------------------------------------------------------------------------

     double gen_weight = Weight;
     if ( isData ) gen_weight = 1.0;

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue(   "PassJSON"                      , passedJSON, gen_weight * pileup_weight  ) ; 
     
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
     
     // Pass HLT
     
     int passHLT = 1;
     if ( isData ) { 
       passHLT = 0;
       if ( H_Ele30_PFJet100_25      == 1 ||
	    H_Ele30_PFNoPUJet100_25  == 1 ){
       	 passHLT = 1;
       }
     }

     // What dataset is this?

     // Electrons
     int PassNEle = 0;
     if ( nEle_store == 2 ) PassNEle = 1;
     
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
     fillVariableWithValue(   "nJet"                          , nJet_store  , gen_weight * pileup_weight  ) ;
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
       
       fillVariableWithValue( "sT_eejj"                      , sT_eejj, gen_weight * pileup_weight  ) ;
       
       //--------------------------------------------------------------------------
       // Fill final selection cuts
       //--------------------------------------------------------------------------
       
       char cut_name[100];
       for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
	 int lq_mass = LQ_MASS[i_lq_mass];
	 sprintf(cut_name, "M_e1e2_LQ%d"  , lq_mass ); fillVariableWithValue ( cut_name, M_e1e2  , gen_weight * pileup_weight  ) ;
	 sprintf(cut_name, "sT_eejj_LQ%d" , lq_mass ); fillVariableWithValue ( cut_name, sT_eejj , gen_weight * pileup_weight  ) ;
	 sprintf(cut_name, "min_M_ej_LQ%d", lq_mass ); fillVariableWithValue ( cut_name, M_ej_min, gen_weight * pileup_weight  ) ;
       }
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
     
     //--------------------------------------------------------------------------
     // Did we pass any final selections?
     //--------------------------------------------------------------------------

     passed_vector.clear();
     for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
       int lq_mass = LQ_MASS[i_lq_mass];
       char cut_name[200]; sprintf(cut_name, "M_e1e2_LQ%d", lq_mass );
       bool decision = bool ( passedAllPreviousCuts(cut_name) && passedCut (cut_name));
       passed_vector.push_back (decision);
     }
     
     bool passed_minimum      = ( passedAllPreviousCuts("PassTrackingFailure") && passedCut ("PassTrackingFailure"));
     bool passed_preselection = ( passedAllPreviousCuts("M_e1e2") && passedCut ("M_e1e2") );

     if ( passed_minimum && isData ){ 
       FillUserTH1D ("run_HLT", run );
       profile_run_vs_nvtx_HLT -> Fill ( run, nVertex, 1 ) ;
     }
     
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
       
       FillUserTH1D("nElectron_PAS"        , nEle_store        , pileup_weight * gen_weight ) ;
       FillUserTH1D("nMuon_PAS"            , nMuon_store       , pileup_weight * gen_weight ) ;
       FillUserTH1D("nJet_PAS"             , nJet_store        , pileup_weight * gen_weight ) ;
       FillUserTH1D("Pt1stEle_PAS"	   , Ele1_Pt           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Eta1stEle_PAS"	   , Ele1_Eta          , pileup_weight * gen_weight ) ;
       FillUserTH1D("Phi1stEle_PAS"	   , Ele1_Phi          , pileup_weight * gen_weight ) ;
       FillUserTH1D("Pt2ndEle_PAS"	   , Ele2_Pt           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Eta2ndEle_PAS"	   , Ele2_Eta          , pileup_weight * gen_weight ) ;
       FillUserTH1D("Phi2ndEle_PAS"	   , Ele2_Phi          , pileup_weight * gen_weight ) ;
       FillUserTH1D("Charge1stEle_PAS"	   , Ele1_Charge       , pileup_weight * gen_weight ) ;
       FillUserTH1D("Charge2ndEle_PAS"	   , Ele2_Charge       , pileup_weight * gen_weight ) ;
       FillUserTH1D("MET_PAS"              , PFMET_Type01XY_Pt , pileup_weight * gen_weight ) ;
       FillUserTH1D("METPhi_PAS"	   , PFMET_Type01XY_Phi, pileup_weight * gen_weight ) ;
       FillUserTH1D("Pt1stJet_PAS"         , Jet1_Pt           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Pt2ndJet_PAS"         , Jet2_Pt           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Eta1stJet_PAS"        , Jet1_Eta          , pileup_weight * gen_weight ) ;
       FillUserTH1D("Eta2ndJet_PAS"        , Jet2_Eta          , pileup_weight * gen_weight ) ;
       FillUserTH1D("Phi1stJet_PAS"	   , Jet1_Phi          , pileup_weight * gen_weight ) ;
       FillUserTH1D("Phi2ndJet_PAS"	   , Jet2_Phi          , pileup_weight * gen_weight ) ;
       FillUserTH1D("sTlep_PAS"            , Ele1_Pt + Ele2_Pt , pileup_weight * gen_weight ) ;
       FillUserTH1D("sTjet_PAS"            , Jet1_Pt + Jet2_Pt , pileup_weight * gen_weight ) ;
       FillUserTH1D("sT_PAS"               , sT_eejj           , pileup_weight * gen_weight ) ;
       FillUserTH1D("Mjj_PAS"		   , M_j1j2            , pileup_weight * gen_weight ) ;
       FillUserTH1D("Mee_PAS"		   , M_e1e2            , pileup_weight * gen_weight ) ;
       FillUserTH1D("Me1j1_PAS"            , M_e1j1            , pileup_weight * gen_weight ) ;
       FillUserTH1D("Me1j2_PAS"            , M_e1j2            , pileup_weight * gen_weight ) ;
       FillUserTH1D("Me2j1_PAS"            , M_e2j1            , pileup_weight * gen_weight ) ;
       FillUserTH1D("Me2j2_PAS"            , M_e2j2            , pileup_weight * gen_weight ) ;
       FillUserTH1D("Ptee_PAS"             , Pt_e1e2           , pileup_weight * gen_weight ) ;
       FillUserTH1D("DCotTheta1stEle_PAS"  , Ele1_DCotTheta    , pileup_weight * gen_weight ) ;
       FillUserTH1D("Dist1stEle_PAS"       , Ele1_Dist         , pileup_weight * gen_weight ) ;
       FillUserTH1D("DCotTheta2ndEle_PAS"  , Ele2_DCotTheta    , pileup_weight * gen_weight ) ;
       FillUserTH1D("Dist2ndEle_PAS"       , Ele2_Dist         , pileup_weight * gen_weight ) ;
       FillUserTH1D("nVertex_PAS"          , nVertex           , pileup_weight * gen_weight ) ;
       FillUserTH1D("DR_Ele1Jet1_PAS"	   , DR_Ele1Jet1       , pileup_weight * gen_weight ) ;
       FillUserTH1D("DR_Ele1Jet2_PAS"	   , DR_Ele1Jet2       , pileup_weight * gen_weight ) ;
       FillUserTH1D("DR_Ele2Jet1_PAS"	   , DR_Ele2Jet1       , pileup_weight * gen_weight ) ;
       FillUserTH1D("DR_Ele2Jet2_PAS"	   , DR_Ele2Jet2       , pileup_weight * gen_weight ) ;
       FillUserTH1D("DR_Jet1Jet2_PAS"	   , DR_Jet1Jet2       , pileup_weight * gen_weight ) ;
       
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
       if ( nJet_store > 2 && Jet3_Pt > 0.0) { 
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
       met.SetPtEtaPhiM (PFMET_Type01XY_Pt, 0.0, PFMET_Type01XY_Phi , 0.0 );
       
       double DR_Ele1Ele2 = e1.DeltaR( e2 ) ;
       FillUserTH1D("DR_Ele1Ele2_PAS"	   , DR_Ele1Ele2      , pileup_weight * gen_weight ) ;
       
       eejj = e1 + e2 + j1 + j2 ; 
       double M_eejj = eejj.M();

       FillUserTH1D("Meejj_PAS", M_eejj , pileup_weight * gen_weight );

       
       if ( fabs(M_e1j1-M_e2j2) < fabs(M_e1j2-M_e2j1) )  {
	 FillUserTH1D("Mej_selected_avg_PAS", M_ej_avg, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Me1j_selected_PAS"   , M_e1j1  , pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Me2j_selected_PAS"   , M_e2j2  , pileup_weight * gen_weight) ;	   
	 FillUserTH2D( "Me1jVsMe2j_selected", M_e1j1  , M_e2j2, pileup_weight * gen_weight ) ;
	 FillUserTH2D( "Me1jVsMe2j_rejected", M_e1j2  , M_e2j1, pileup_weight * gen_weight ) ;
       }

       else {
	 FillUserTH1D("Mej_selected_avg_PAS", M_ej_avg, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Me1j_selected_PAS"   , M_e1j2, pileup_weight * gen_weight) ;	   
	 FillUserTH1D("Me2j_selected_PAS"   , M_e2j1, pileup_weight * gen_weight) ;	   
	 FillUserTH2D( "Me1jVsMe2j_selected", M_e1j2, M_e2j1, pileup_weight * gen_weight ) ;
	 FillUserTH2D( "Me1jVsMe2j_rejected", M_e1j1, M_e2j2, pileup_weight * gen_weight ) ;
       }


       for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
	 int  lq_mass = LQ_MASS      [i_lq_mass];
	 bool pass    = passed_vector[i_lq_mass];
	 if ( !pass ) continue;
	 
	 sprintf(plot_name, "Mej_selected_avg_LQ%d"       , lq_mass ); FillUserTH1D ( plot_name, M_ej_avg          , pileup_weight * gen_weight );
	 sprintf(plot_name, "Mej_selected_min_LQ%d"       , lq_mass ); FillUserTH1D ( plot_name, M_ej_min          , pileup_weight * gen_weight );
	 sprintf(plot_name, "Mej_selected_max_LQ%d"       , lq_mass ); FillUserTH1D ( plot_name, M_ej_max          , pileup_weight * gen_weight );
	 sprintf(plot_name, "Mej_minmax_LQ%d"             , lq_mass ); FillUserTH1D ( plot_name, M_ej_min          , pileup_weight * gen_weight );
	 sprintf(plot_name, "Mej_minmax_LQ%d"             , lq_mass ); FillUserTH1D ( plot_name, M_ej_max          , pileup_weight * gen_weight );
	 sprintf(plot_name, "sT_eejj_LQ%d"                , lq_mass ); FillUserTH1D ( plot_name, sT_eejj           , pileup_weight * gen_weight );
	 sprintf(plot_name, "Mee_LQ%d"                    , lq_mass ); FillUserTH1D ( plot_name, M_e1e2            , pileup_weight * gen_weight );
	 sprintf(plot_name, "Mej_selected_min_vs_max_LQ%d", lq_mass ); FillUserTH2D ( plot_name, M_ej_min, M_ej_max, pileup_weight * gen_weight );
	 
       }
     }
   } // End loop over events

   if ( rootEngine ) delete rootEngine ;
   
   output_root_ -> cd();
   profile_run_vs_nvtx_HLT -> Write();
   profile_run_vs_nvtx_PAS -> Write();
   
   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

