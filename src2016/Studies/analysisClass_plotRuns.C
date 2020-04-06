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
    
   //--------------------------------------------------------------------------
   // Final selection mass points
   //--------------------------------------------------------------------------

   const int n_lq_mass = 19;

   int LQ_MASS[n_lq_mass] = { 
     300,  350,  400, 450, 500, 550,  600,  650,
     700,  750,  800, 850, 900, 950, 1000, 1050,
     1100, 1150, 1200
   };

   //--------------------------------------------------------------------------
   // Decide which plots to save (default is to save everything)
   //--------------------------------------------------------------------------
   
   fillSkim                         ( !true  ) ;
   fillAllPreviousCuts              ( !true  ) ;
   fillAllOtherCuts                 ( !true  ) ;
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;

   
   //////////book histos here
   CreateUserTH1D( "run_HLT"                         ,    15000 , 246000  , 261000 );

   /////////initialize variables

   Long64_t nentries = fChain->GetEntries();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   //--------------------------------------------------------------------------
   // Make event lists for various filters
   //--------------------------------------------------------------------------
   // XXX NB: the hcal noise filter event lists are incomplete!
   //    Need to rentuplize everything and run the filters at that point
   EventListHelper HBHENoiseRun2LooseEventList;
   HBHENoiseRun2LooseEventList.addFileToList("eventlist_hbher2l.txt");
   EventListHelper HBHENoiseRun2IsoEventList;
   HBHENoiseRun2IsoEventList.addFileToList("eventlist_hbheiso.txt");
   EventListHelper CSCBeamHaloTight2015EventList;
   CSCBeamHaloTight2015EventList.addFileToList("csc2015_Dec01.txt");
   EventListHelper ECALFourthBadSCEventList;
   ECALFourthBadSCEventList.addFileToList("ecalscn1043093_Dec01.txt");
   
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
     if ( isData && Ele2_ValidFrac > 998. ){
       gen_weight = 0.0;
       if      (  60.0 < M_e1e2 < 120. ) gen_weight = 0.61;
       else if ( 120.0 < M_e1e2 < 200. ) gen_weight = 0.42;
       else if ( 200.0 < M_e1e2        ) gen_weight = 0.42;
     }

     // std::cout << "Gen weight = " << int ( 1.0 / gen_weight ) << std::endl;

     //--------------------------------------------------------------------------
     // First variable to fill just shows the "reweighting".  Always passes.
     //--------------------------------------------------------------------------

     fillVariableWithValue ( "Reweighting", 1, gen_weight * pileup_weight  );

     //--------------------------------------------------------------------------
     // Fill JSON variable
     //--------------------------------------------------------------------------

     // JSON variable
     fillVariableWithValue ("PassJSON", passedJSON, gen_weight * pileup_weight); 

     //--------------------------------------------------------------------------
     // Fill noise filters
     //--------------------------------------------------------------------------

     // Noise/MET filters
     // see: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
     // XXX Right now, we have to read the bad events from event lists
     // txt files sent to me privately
     // see: https://indico.cern.ch/event/458729/contribution/8/attachments/1179252/1706416/metscan.pdf
     // Later this will be included in a re-reco
     //fillVariableWithValue(   "PassHBHENoiseFilter"	      , PassHBHENoiseFilter                              , gen_weight * pileup_weight );
     //fillVariableWithValue(   "PassHBHENoiseIsoFilter"	      , PassHBHENoiseIsoFilter                              , gen_weight * pileup_weight );
     //fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight                          , gen_weight * pileup_weight );
     //TODO to be implemented in skim?
     //fillVariableWithValue(   "PassBadEESupercrystalFilter"   , ( isData == 1 ) ? PassBadEESupercrystalFilter : 1, gen_weight * pileup_weight );
     
     //fillVariableWithValue(   "PassBeamScraping"	      , ( isData == 1 ) ? PassBeamScraping	      : 1, gen_weight * pileup_weight );
     //fillVariableWithValue(   "PassEcalDeadCellBoundEnergy"   , PassEcalDeadCellBoundEnergy                      , gen_weight * pileup_weight );
     //fillVariableWithValue(   "PassEcalDeadCellTrigPrim"      , PassEcalDeadCellTrigPrim                         , gen_weight * pileup_weight );
     //fillVariableWithValue(   "PassEcalLaserCorrFilter"       , ( isData == 1 ) ? PassEcalLaserCorrFilter     : 1, gen_weight * pileup_weight );
     //fillVariableWithValue(   "PassHcalLaserEventFilter"      , ( isData == 1 ) ? PassHcalLaserEventFilter    : 1, gen_weight * pileup_weight );
     fillVariableWithValue(   "PassHBHENoiseFilter"	        , HBHENoiseRun2LooseEventList.eventInList(run,ls,event)    ?  0 : 1  , gen_weight * pileup_weight );
     fillVariableWithValue(   "PassHBHENoiseIsoFilter"	    , HBHENoiseRun2IsoEventList.eventInList(run,ls,event)      ?  0 : 1  , gen_weight * pileup_weight );
     fillVariableWithValue(   "PassCSCBeamHaloFilterTight"  , CSCBeamHaloTight2015EventList.eventInList(run,ls,event)  ?  0 : 1  , gen_weight * pileup_weight );
     // for 4th bad SC
     fillVariableWithValue(   "PassBadEESupercrystalFilter" , ECALFourthBadSCEventList.eventInList(run,ls,event)       ?  0 : 1  , gen_weight * pileup_weight );
     //
     fillVariableWithValue(   "PassPhysDecl"		      , ( isData == 1 ) ? PassPhysDecl		      : 1, gen_weight * pileup_weight );
     fillVariableWithValue(   "PassPrimaryVertex"	      , PassPrimaryVertex                                , gen_weight * pileup_weight );
     //fillVariableWithValue(   "PassTrackingFailure"	      , ( isData == 1 ) ? PassTrackingFailure	      : 1, gen_weight * pileup_weight );
     
     //--------------------------------------------------------------------------
     // Fill HLT
     //--------------------------------------------------------------------------

     //int passHLT = 1;
     //if ( isData ) { 
     //  passHLT = 0;
     //  if ( H_Ele30_PFJet100_25 == 1 || H_Ele30_PFNoPUJet100_25  == 1 ){
	 ////if ( H_DoubleEle33_CIdL_GsfIdVL == 1 ) { 
     //  	 passHLT = 1;
     //  }
     //}
     int passHLT = 1;
     if ( isData ) { 
       passHLT = 0;
       if ( H_Ele45_PFJet200_PFJet50 == 1)
         passHLT = 1;
     }

     fillVariableWithValue ( "PassHLT", passHLT, gen_weight * pileup_weight  ) ;     

     //--------------------------------------------------------------------------
     // Calculate variables for trigger matching 
     //--------------------------------------------------------------------------

     int nEle_hltMatched = 0.0;
     if ( Ele1_hltEleSignalPt > 0.0 ) nEle_hltMatched++;
     if ( Ele2_hltEleSignalPt > 0.0 ) nEle_hltMatched++;
     
     int nJet_hltMatched = 0.0;
     if ( Jet1_hltNoPUJetPt > 0.0 || Jet1_hltJetPt > 0.0 ) nJet_hltMatched++;
     if ( Jet2_hltNoPUJetPt > 0.0 || Jet2_hltJetPt > 0.0 ) nJet_hltMatched++;

     fillVariableWithValue("nEle_hltMatched",nEle_hltMatched, gen_weight * pileup_weight  );
     fillVariableWithValue("nJet_hltMatched",nJet_hltMatched, gen_weight * pileup_weight  );
     
     //--------------------------------------------------------------------------
     // Pass number of muons & electrons 
     // --> Special consideration if ttbar is derived from data
     //--------------------------------------------------------------------------

     // Muons and electrons
     bool is_ttbar_from_data = false;
     if ( Ele2_ValidFrac > 998. ) is_ttbar_from_data = true;
     
     int PassNEle = 0;
     if ( !is_ttbar_from_data && nEle_ptCut == 2 ) PassNEle = 1;
     if (  is_ttbar_from_data && nEle_ptCut == 2 ) PassNEle = 1;

     int PassNMuon = 0;
     if ( !is_ttbar_from_data && nMuon_ptCut == 0 ) PassNMuon = 1;
     if (  is_ttbar_from_data && nMuon_ptCut >  0 ) PassNMuon = 1;

     fillVariableWithValue("PassNEle" , PassNEle , gen_weight * pileup_weight);
     fillVariableWithValue("PassNMuon", PassNMuon, gen_weight * pileup_weight);

     //--------------------------------------------------------------------------
     // Calculate electron-jet pair mass values
     //--------------------------------------------------------------------------
     
     double M_ej_avg, M_ej_min, M_ej_max;
     
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
     }
     
     double sT_zjj = Pt_e1e2 + Jet1_Pt + Jet2_Pt;

     //--------------------------------------------------------------------------
     // Fill electron variables 
     //--------------------------------------------------------------------------
     
     if ( nEle_store >= 1 ) fillVariableWithValue( "Ele1_Pt", Ele1_Pt, gen_weight * pileup_weight  ) ;
     if ( nEle_store >= 2 ) fillVariableWithValue( "Ele2_Pt", Ele2_Pt, gen_weight * pileup_weight  ) ;
			
     //--------------------------------------------------------------------------
     // Fill jet variables 
     //--------------------------------------------------------------------------
					    
     // Jets								    
     fillVariableWithValue("nJet", nJet_ptCut, gen_weight * pileup_weight );
     if ( nJet_store >= 1 ) { 						    
       fillVariableWithValue( "Jet1_Pt"    , Jet1_Pt     , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "Jet1_Eta"   , Jet1_Eta    , gen_weight * pileup_weight  ) ;
     }
     if ( nJet_store >= 2 ) { 
       fillVariableWithValue( "Jet2_Pt"    , Jet2_Pt     , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "Jet2_Eta"   , Jet2_Eta    , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "DR_Jet1Jet2", DR_Jet1Jet2 , gen_weight * pileup_weight  ) ;
     }

     //--------------------------------------------------------------------------
     // Fill DeltaR variables
     //--------------------------------------------------------------------------

     if ( nEle_store >= 2 && nJet_store >= 1) {
       fillVariableWithValue( "DR_Ele1Jet1"  , DR_Ele1Jet1 , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "DR_Ele2Jet1"  , DR_Ele2Jet1 , gen_weight * pileup_weight  ) ;
       if(nJet_store >= 2) {
	 fillVariableWithValue( "DR_Ele1Jet2", DR_Ele1Jet2 , gen_weight * pileup_weight  ) ;
	 fillVariableWithValue( "DR_Ele2Jet2", DR_Ele2Jet2 , gen_weight * pileup_weight  ) ;
       }
     }


     //--------------------------------------------------------------------------
     // Multi-object variables
     //--------------------------------------------------------------------------

     if ( nEle_store >= 2 ) { 						    
       fillVariableWithValue( "M_e1e2"     , M_e1e2 , gen_weight * pileup_weight  ) ;
       fillVariableWithValue( "M_e1e2_opt" , M_e1e2 , gen_weight * pileup_weight  ) ;
       
       if ( nJet_store >= 2 ) { 
	 fillVariableWithValue( "sT_eejj"    , sT_eejj , gen_weight * pileup_weight  ) ;
	 fillVariableWithValue( "sT_eejj_opt", sT_eejj , gen_weight * pileup_weight  ) ;
	 fillVariableWithValue( "Mej_min_opt", M_ej_min, gen_weight * pileup_weight  ) ;
       }      
     }

     //--------------------------------------------------------------------------
     // Fill final selection cuts
     //--------------------------------------------------------------------------

     char cut_name[100];
     for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
       int lq_mass = LQ_MASS[i_lq_mass];
       //XXX Only look at specific selections for now
       if(lq_mass!=300 && lq_mass!=600 && lq_mass!=650 && lq_mass!=650 && lq_mass!=1200) continue;
       // end Only look at specific selections for now
       sprintf(cut_name, "M_e1e2_LQ%d"  , lq_mass ); fillVariableWithValue ( cut_name, M_e1e2  , gen_weight * pileup_weight  ) ;
       sprintf(cut_name, "sT_eejj_LQ%d" , lq_mass ); fillVariableWithValue ( cut_name, sT_eejj , gen_weight * pileup_weight  ) ;
       sprintf(cut_name, "min_M_ej_LQ%d", lq_mass ); fillVariableWithValue ( cut_name, M_ej_min, gen_weight * pileup_weight  ) ;
     }
     
     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Did we at least pass the noise filtes?
     //--------------------------------------------------------------------------
     
     //bool passed_minimum = ( passedAllPreviousCuts("PassTrackingFailure") && passedCut ("PassTrackingFailure"));
     bool passed_minimum = ( passedAllPreviousCuts("PassPrimaryVertex") && passedCut ("PassPrimaryVertex"));
     
     //--------------------------------------------------------------------------
     // Did we pass preselection?
     //--------------------------------------------------------------------------
     
     bool passed_preselection = ( passedAllPreviousCuts("M_e1e2") && passedCut ("M_e1e2") );

     if ( passed_preselection ) {
       FillUserTH1D ("run_HLT", run );
     }


     ////////////////////// User's code ends here ///////////////////////

   } // End loop over events


   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
