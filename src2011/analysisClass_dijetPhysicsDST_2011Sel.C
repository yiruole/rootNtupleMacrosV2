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
   fillAllPreviousCuts              ( true  ) ;
   fillAllOtherCuts                 ( !true  ) ;
   fillAllSameLevelAndLowerLevelCuts( true  ) ;
   fillAllCuts                      ( !true  ) ;

   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------

   CreateUserTH1D( "M_PFJet1PFJet2"	  , 
 		   getHistoNBins("M_PFJet1PFJet2"), 
 		   getHistoMin("M_PFJet1PFJet2"), 
 		   getHistoMax("M_PFJet1PFJet2")     ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 

   CreateUserTH1D( "M_PFJet1PFJet2_nVtxL5"	  , 
 		   getHistoNBins("M_PFJet1PFJet2"), 
 		   getHistoMin("M_PFJet1PFJet2"), 
 		   getHistoMax("M_PFJet1PFJet2")     ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_nVtxL5"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 

   CreateUserTH1D( "M_PFJet1PFJet2_nVtx5To10"	  , 
 		   getHistoNBins("M_PFJet1PFJet2"), 
 		   getHistoMin("M_PFJet1PFJet2"), 
 		   getHistoMax("M_PFJet1PFJet2")     ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_nVtx5To10"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 

   CreateUserTH1D( "M_PFJet1PFJet2_nVtx10To15"	  , 
 		   getHistoNBins("M_PFJet1PFJet2"), 
 		   getHistoMin("M_PFJet1PFJet2"), 
 		   getHistoMax("M_PFJet1PFJet2")     ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_nVtx10To15"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 

   CreateUserTH1D( "M_PFJet1PFJet2_nVtxM15"	  , 
 		   getHistoNBins("M_PFJet1PFJet2"), 
 		   getHistoMin("M_PFJet1PFJet2"), 
 		   getHistoMax("M_PFJet1PFJet2")     ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_nVtxM15"	  , 
 		   getHistoNBins("M_FatPFJet1FatPFJet2"), 
 		   getHistoMin("M_FatPFJet1FatPFJet2"), 
 		   getHistoMax("M_FatPFJet1FatPFJet2")     ) ; 

   //Shifted in 0.3 - 1.3 range (divided by nominal signal mass)
   CreateUserTH1D( "M_PFJet1PFJet2_shifted_1"	  , 50, 0.3, 1.3  ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_shifted_1"   , 50, 0.3, 1.3  ) ; 
   CreateUserTH1D( "M_PFJet1PFJet2_shifted_2"	  , 50, 0.3, 1.3  ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_shifted_2"   , 50, 0.3, 1.3  ) ; 
   CreateUserTH1D( "M_PFJet1PFJet2_shifted_3"	  , 50, 0.3, 1.3  ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_shifted_3"   , 50, 0.3, 1.3  ) ; 
   CreateUserTH1D( "M_PFJet1PFJet2_shifted_4"	  , 50, 0.3, 1.3  ) ; 
   CreateUserTH1D( "M_FatPFJet1FatPFJet2_shifted_4"   , 50, 0.3, 1.3  ) ; 
   
   Double_t massBins[84] = {1, 3, 6, 10, 16, 23, 31, 40, 50, 61, 74, 88, 103, 119, 137, 156, 176, 197, 220, 244, 270, 296, 325, 354, 386, 419, 453, 489, 526, 565, 606, 649, 693, 740, 788, 838, 890, 944, 1000, 1058, 1118, 1181, 1246, 1313, 1383, 1455, 1530, 1607, 1687, 1770, 1856, 1945, 2037, 2132, 2231, 2332, 2438, 2546, 2659, 2775, 2895, 3019, 3147, 3279, 3416, 3558, 3704, 3854, 4010, 4171, 4337, 4509, 4686, 4869, 5058, 5253, 5455, 5663, 5877, 6099, 6328, 6564, 6808, 7000};
   
   TH1D *M_PFJet1PFJet2_VarBin = new TH1D ("M_PFJet1PFJet2_VarBin","M_PFJet1PFJet2_VarBin",83,massBins);
   M_PFJet1PFJet2_VarBin->Sumw2();
   TH1D *M_FatPFJet1FatPFJet2_VarBin = new TH1D ("M_FatPFJet1FatPFJet2_VarBin","M_FatPFJet1FatPFJet2_VarBin",83,massBins);
   M_FatPFJet1FatPFJet2_VarBin->Sumw2();

   TH1D *M_PFJet1PFJet2_VarBin_nVtxL5 = new TH1D ("M_PFJet1PFJet2_VarBin_nVtxL5","M_PFJet1PFJet2_VarBin_nVtxL5",83,massBins);
   M_PFJet1PFJet2_VarBin_nVtxL5->Sumw2();
   TH1D *M_FatPFJet1FatPFJet2_VarBin_nVtxL5 = new TH1D ("M_FatPFJet1FatPFJet2_VarBin_nVtxL5","M_FatPFJet1FatPFJet2_VarBin_nVtxL5",83,massBins);
   M_FatPFJet1FatPFJet2_VarBin_nVtxL5->Sumw2();

   TH1D *M_PFJet1PFJet2_VarBin_nVtx5To10 = new TH1D ("M_PFJet1PFJet2_VarBin_nVtx5To10","M_PFJet1PFJet2_VarBin_nVtx5To10",83,massBins);
   M_PFJet1PFJet2_VarBin_nVtx5To10->Sumw2();
   TH1D *M_FatPFJet1FatPFJet2_VarBin_nVtx5To10 = new TH1D ("M_FatPFJet1FatPFJet2_VarBin_nVtx5To10","M_FatPFJet1FatPFJet2_VarBin_nVtx5To10",83,massBins);
   M_FatPFJet1FatPFJet2_VarBin_nVtx5To10->Sumw2();

   TH1D *M_PFJet1PFJet2_VarBin_nVtx10To15 = new TH1D ("M_PFJet1PFJet2_VarBin_nVtx10To15","M_PFJet1PFJet2_VarBin_nVtx10To15",83,massBins);
   M_PFJet1PFJet2_VarBin_nVtx10To15->Sumw2();
   TH1D *M_FatPFJet1FatPFJet2_VarBin_nVtx10To15 = new TH1D ("M_FatPFJet1FatPFJet2_VarBin_nVtx10To15","M_FatPFJet1FatPFJet2_VarBin_nVtx10To15",83,massBins);
   M_FatPFJet1FatPFJet2_VarBin_nVtx10To15->Sumw2();

   TH1D *M_PFJet1PFJet2_VarBin_nVtxM15 = new TH1D ("M_PFJet1PFJet2_VarBin_nVtxM15","M_PFJet1PFJet2_VarBin_nVtxM15",83,massBins);
   M_PFJet1PFJet2_VarBin_nVtxM15->Sumw2();
   TH1D *M_FatPFJet1FatPFJet2_VarBin_nVtxM15 = new TH1D ("M_FatPFJet1FatPFJet2_VarBin_nVtxM15","M_FatPFJet1FatPFJet2_VarBin_nVtxM15",83,massBins);
   M_FatPFJet1FatPFJet2_VarBin_nVtxM15->Sumw2();
 
   //--------------------------------------------------------------------------
   // Precuts
   //--------------------------------------------------------------------------

   double printMass = getPreCutValue1("printMass");
   double nominalSignalMass1 = getPreCutValue1("nominalSignalMass");
   double nominalSignalMass2 = getPreCutValue2("nominalSignalMass");
   double nominalSignalMass3 = getPreCutValue3("nominalSignalMass");
   double nominalSignalMass4 = getPreCutValue4("nominalSignalMass");

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
     
     int passedJSON = passJSON ( run , ls , isData ) ;

     //--------------------------------------------------------------------------
     // Event-weight
     //--------------------------------------------------------------------------
     
     //     double event_weight = 1;

     //--------------------------------------------------------------------------
     // Do pileup re-weighting
     //--------------------------------------------------------------------------
     
     //      int NPILEUP_AVE = int( (nPileUpInt_BXminus1 + nPileUpInt_BX0 + nPileUpInt_BXplus1)/3 );
     //      int NPILEUP_FINAL = min( NPILEUP_AVE , 25 );
     //      event_weight = getPileupWeight ( NPILEUP_FINAL, isData ) ;
     //      //event_weight = getPileupWeight ( min(nPileUpInt_BX0,25), isData ) ;

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------

     // Gen level info     

     //gg --> G --> gg
     if( ProcessID == 392 && is_gg == 1 )
       fillVariableWithValue(   "ggFinalState"                      , 1 ) ;
     else
       fillVariableWithValue(   "ggFinalState"                      , 0 ) ;

     //qg --> Qstar --> qg
     if( (ProcessID == 147 || ProcessID == 148) && is_qg == 1 )
       fillVariableWithValue(   "qgFinalState"                      , 1 ) ;
     else
       fillVariableWithValue(   "qgFinalState"                      , 0 ) ;

     //qq --> G --> qq
     if( ProcessID == 391 && is_qq == 1 )
       fillVariableWithValue(   "qqFinalState"                      , 1 ) ;
     else
       fillVariableWithValue(   "qqFinalState"                      , 0 ) ;
     
     // JSON variable
     fillVariableWithValue(   "PassJSON"                      , passedJSON ) ; 
     //fillVariableWithValue(   "PassJSON"                      , PassJSON ) ;      

     // Event filters
     fillVariableWithValue(   "PassPrimaryVertex"             , PassPrimaryVertex ) ; 

     //Vertices
     fillVariableWithValue(   "nVertex_good"                  , nVertex_good ) ; 

     //HLT     
     fillVariableWithValue(   "HLT_HT350orM400"               , (HLT_HT350 || HLT_FatJetMass400) ) ; 
     fillVariableWithValue(   "HLT_M300"                      , (HLT_FatJetMass300) ) ; 
     fillVariableWithValue(   "HLT_M400"                      , (HLT_FatJetMass400) ) ; 
     fillVariableWithValue(   "HLT_HT350"                     , (HLT_HT350) ) ; 

     //PFJets
     fillVariableWithValue(   "PFJet1_Pt"    , PFJet1_Pt ) ; 
     fillVariableWithValue(   "PFJet1_Eta"    , PFJet1_Eta ) ; 
     fillVariableWithValue(   "PFJet1_PassJetID"    , PFJet1_PassJetID ) ; 

     fillVariableWithValue(   "PFJet2_Pt"    , PFJet2_Pt ) ; 
     fillVariableWithValue(   "PFJet2_Eta"    , PFJet2_Eta ) ; 
     fillVariableWithValue(   "PFJet2_PassJetID"    , PFJet2_PassJetID ) ; 

     fillVariableWithValue(   "DPhi_PFJet1PFJet2"    , DPhi_PFJet1PFJet2 ) ; 
     fillVariableWithValue(   "DEta_PFJet1PFJet2"    , DEta_PFJet1PFJet2 ) ; 
     fillVariableWithValue(   "M_PFJet1PFJet2"       , M_PFJet1PFJet2 ) ; 

     //PFFatJets
     fillVariableWithValue(   "DEta_FatPFJet1FatPFJet2"   , DEta_FatPFJet1FatPFJet2 ) ; 
     fillVariableWithValue(   "DPhi_FatPFJet1FatPFJet2"   , DPhi_FatPFJet1FatPFJet2 ) ; 
     fillVariableWithValue(   "M_FatPFJet1FatPFJet2"      , M_FatPFJet1FatPFJet2 ) ; 
     
     //--------------------------------------------------------------------------
     // Evaluate the cuts
     //--------------------------------------------------------------------------
     
     evaluateCuts();

     //--------------------------------------------------------------------------
     // Event selection
     //--------------------------------------------------------------------------

     bool pass_PFJet_selection       = false;
     bool pass_PFFatJet_selection    = false;

     //PFJet selection
     if( passedCut("0") 
	 && passedCut("1") )
       pass_PFJet_selection = true;

     //PFFatJet selection
     if( passedCut("0") 
	 && passedCut("2") )
       pass_PFFatJet_selection = true;

     //--------------------------------------------------------------------------
     // Fill histograms
     //--------------------------------------------------------------------------

     //PFJet selection
     if( pass_PFJet_selection )
       {
	 if(printMass)
	   {
	     cout << "PrintPFDiJetMass " << getVariableValue("M_PFJet1PFJet2") << endl;
	   }

	 FillUserTH1D( "M_PFJet1PFJet2"	                         , getVariableValue("M_PFJet1PFJet2") );	 	 
	 M_PFJet1PFJet2_VarBin->Fill( getVariableValue("M_PFJet1PFJet2") );

	 if(!isData)
	   {
	     FillUserTH1D( "M_PFJet1PFJet2_shifted_1" , getVariableValue("M_PFJet1PFJet2") / nominalSignalMass1 ); 	     
	     FillUserTH1D( "M_PFJet1PFJet2_shifted_2" , getVariableValue("M_PFJet1PFJet2") / nominalSignalMass2 ); 	     
	     FillUserTH1D( "M_PFJet1PFJet2_shifted_3" , getVariableValue("M_PFJet1PFJet2") / nominalSignalMass3 ); 	     
	     FillUserTH1D( "M_PFJet1PFJet2_shifted_4" , getVariableValue("M_PFJet1PFJet2") / nominalSignalMass4 ); 	     
	   }

	 if( getVariableValue("nVertex_good") < 5 )
	   {
	     FillUserTH1D( "M_PFJet1PFJet2_nVtxL5"	         , getVariableValue("M_PFJet1PFJet2") );	 	 	     
	     M_PFJet1PFJet2_VarBin_nVtxL5->Fill( getVariableValue("M_PFJet1PFJet2") );
	   }

	 if( getVariableValue("nVertex_good") >= 5 && getVariableValue("nVertex_good") < 10 )
	   {
	     FillUserTH1D( "M_PFJet1PFJet2_nVtx5To10"	         , getVariableValue("M_PFJet1PFJet2") );	 	 	     
	     M_PFJet1PFJet2_VarBin_nVtx5To10->Fill( getVariableValue("M_PFJet1PFJet2") );
	   }

	 if( getVariableValue("nVertex_good") >= 10 && getVariableValue("nVertex_good") < 15 )
	   {
	     FillUserTH1D( "M_PFJet1PFJet2_nVtx10To15"	         , getVariableValue("M_PFJet1PFJet2") );	 	 
	     M_PFJet1PFJet2_VarBin_nVtx10To15->Fill( getVariableValue("M_PFJet1PFJet2") );
	   }

	 if( getVariableValue("nVertex_good") >= 15 )
	   {
	     FillUserTH1D( "M_PFJet1PFJet2_nVtxM15"	         , getVariableValue("M_PFJet1PFJet2") );	 	 
	     M_PFJet1PFJet2_VarBin_nVtxM15->Fill( getVariableValue("M_PFJet1PFJet2") );
	   }
       }

     //PFFatJet selection
     if( pass_PFFatJet_selection )
       {
	 if(printMass)
	   {
	     cout << "PrintPFDiFatJetMass " << getVariableValue("M_FatPFJet1FatPFJet2") << endl;
	   }

	 FillUserTH1D( "M_FatPFJet1FatPFJet2"	 , getVariableValue("M_FatPFJet1FatPFJet2") );	 	 
	 M_FatPFJet1FatPFJet2_VarBin->Fill( getVariableValue("M_FatPFJet1FatPFJet2") );

	 if(!isData)
	   {
	     FillUserTH1D( "M_FatPFJet1FatPFJet2_shifted_1"  , getVariableValue("M_FatPFJet1FatPFJet2") / nominalSignalMass1 ); 
	     FillUserTH1D( "M_FatPFJet1FatPFJet2_shifted_2"  , getVariableValue("M_FatPFJet1FatPFJet2") / nominalSignalMass2 ); 
	     FillUserTH1D( "M_FatPFJet1FatPFJet2_shifted_3"  , getVariableValue("M_FatPFJet1FatPFJet2") / nominalSignalMass3 ); 
	     FillUserTH1D( "M_FatPFJet1FatPFJet2_shifted_4"  , getVariableValue("M_FatPFJet1FatPFJet2") / nominalSignalMass4 ); 
	   }

	 if( getVariableValue("nVertex_good") < 5 )
	   {
	     FillUserTH1D( "M_FatPFJet1FatPFJet2_nVtxL5"	 , getVariableValue("M_FatPFJet1FatPFJet2") );	 	 
	     M_FatPFJet1FatPFJet2_VarBin_nVtxL5->Fill( getVariableValue("M_FatPFJet1FatPFJet2") );
	   }

	 if( getVariableValue("nVertex_good") >= 5 && getVariableValue("nVertex_good") < 10 )
	   {
	     FillUserTH1D( "M_FatPFJet1FatPFJet2_nVtx5To10"	 , getVariableValue("M_FatPFJet1FatPFJet2") );	 	 
	     M_FatPFJet1FatPFJet2_VarBin_nVtx5To10->Fill( getVariableValue("M_FatPFJet1FatPFJet2") );
	   }

	 if( getVariableValue("nVertex_good") >= 10 && getVariableValue("nVertex_good") < 15 )
	   {
	     FillUserTH1D( "M_FatPFJet1FatPFJet2_nVtx10To15"	 , getVariableValue("M_FatPFJet1FatPFJet2") );	 	 
	     M_FatPFJet1FatPFJet2_VarBin_nVtx10To15->Fill( getVariableValue("M_FatPFJet1FatPFJet2") );
	   }

	 if( getVariableValue("nVertex_good") >= 15 )
	   {
	     FillUserTH1D( "M_FatPFJet1FatPFJet2_nVtxM15"	 , getVariableValue("M_FatPFJet1FatPFJet2") );	 	 
	     M_FatPFJet1FatPFJet2_VarBin_nVtxM15->Fill( getVariableValue("M_FatPFJet1FatPFJet2") );	     
	   }
       }

   } // End loop over events

   //////////write histos 

   M_PFJet1PFJet2_VarBin->Write();
   M_FatPFJet1FatPFJet2_VarBin->Write();

   M_PFJet1PFJet2_VarBin_nVtxL5->Write();
   M_FatPFJet1FatPFJet2_VarBin_nVtxL5->Write();

   M_PFJet1PFJet2_VarBin_nVtx5To10->Write();
   M_FatPFJet1FatPFJet2_VarBin_nVtx5To10->Write();

   M_PFJet1PFJet2_VarBin_nVtx10To15->Write();
   M_FatPFJet1FatPFJet2_VarBin_nVtx10To15->Write();

   M_PFJet1PFJet2_VarBin_nVtxM15->Write();
   M_FatPFJet1FatPFJet2_VarBin_nVtxM15->Write();

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
