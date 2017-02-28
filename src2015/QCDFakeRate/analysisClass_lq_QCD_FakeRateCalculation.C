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
  // Print info if asked
  //--------------------------------------------------------------------------

  bool verbose = false;

  //--------------------------------------------------------------------------
  // Decide which plots to save (default is to save everything)
  //--------------------------------------------------------------------------

  fillSkim                         ( !true  ) ;
  fillAllPreviousCuts              ( !true  ) ;
  fillAllOtherCuts                 ( !true  ) ; 
  fillAllSameLevelAndLowerLevelCuts( !true  ) ;
  fillAllCuts                      ( !true  ) ;

  //--------------------------------------------------------------------------
  // Get pre-cut values
  //--------------------------------------------------------------------------

  double eleEta_bar         = getPreCutValue1("eleEta_bar");
  double eleEta_end1_min    = getPreCutValue1("eleEta_end1");
  double eleEta_end1_max    = getPreCutValue2("eleEta_end1");
  double eleEta_end2_min    = getPreCutValue1("eleEta_end2");
  double eleEta_end2_max    = getPreCutValue2("eleEta_end2");

  //--------------------------------------------------------------------------
  // Create TH1D's
  //--------------------------------------------------------------------------

  // Debugging

  CreateUserTH1D("Total_minPrescale"                     ,   200 , 0.0      , 2000.0   ); CreateUserTH1D("Electrons_minPrescale"                       ,   200 , 0.0      , 2000.0   ); CreateUserTH1D("Jets_minPrescale"                       ,   200 , 0.0      , 2000.0   );
  CreateUserTH1D("Total_pileupWeight"                    ,   100 , 0.0      , 2.0      ); CreateUserTH1D("Electrons_pileupWeight"                      ,   100 , 0.0      , 2.0      ); CreateUserTH1D("Jets_pileupWeight"                      ,   100 , 0.0      , 2.0      );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D("Total_Bar_nVertex_PAS"                 ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Electrons_Bar_nVertex_PAS"                   ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Jets_Bar_nVertex_PAS"                   ,   100 , 0.0      , 100.0    ); 
  CreateUserTH1D("Total_End1_nVertex_PAS"                ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Electrons_End1_nVertex_PAS"                  ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Jets_End1_nVertex_PAS"                  ,   100 , 0.0      , 100.0    ); 
  CreateUserTH1D("Total_End2_nVertex_PAS"                ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Electrons_End2_nVertex_PAS"                  ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Jets_End2_nVertex_PAS"                  ,   100 , 0.0      , 100.0    ); 
                                                                                                                                                                                                                                                                                     
  // inclusive                                                                                                                                                                                                                                                                       
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_Bar_nElectron_PAS"              ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_Bar_nElectron_PAS"               ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_Bar_nElectron_PAS"               ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_Bar_Pt1stEle_PAS"               ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Bar_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "Jets_Bar_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
  CreateUserTH1D( "Total_Bar_hltPt1stEle_PAS"            ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Bar_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );  CreateUserTH1D( "Jets_Bar_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );
  CreateUserTH1D( "Total_Bar_Eta1stEle_PAS"	          ,    100 , -5      , 5	         );  CreateUserTH1D( "Electrons_Bar_Eta1stEle_PAS"	           ,    100 , -5      , 5          );  CreateUserTH1D( "Jets_Bar_Eta1stEle_PAS"	           ,    100 , -5      , 5          );
  CreateUserTH1D( "Total_Bar_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416     );  CreateUserTH1D( "Electrons_Bar_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416    );  CreateUserTH1D( "Jets_Bar_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416    );
  CreateUserTH1D( "Total_Bar_Charge1stEle_PAS"	          ,    2   , -1.0001 , 1.0001  );  CreateUserTH1D( "Electrons_Bar_Charge1stEle_PAS"            ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_Bar_Charge1stEle_PAS"            ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_Bar_DCotTheta1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_Bar_DCotTheta1stEle_PAS"         ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_Bar_DCotTheta1stEle_PAS"         ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_Bar_Dist1stEle_PAS"             ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_Bar_Dist1stEle_PAS"              ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_Bar_Dist1stEle_PAS"              ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_Bar_MET_PAS"                    ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_Bar_MET_PAS"                     ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_Bar_MET_PAS"                     ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_Bar_TrkIsoHEEP7_PAS"            ,    100 , 0.0     , 100.     );  CreateUserTH1D( "Electrons_Bar_TrkIsoHEEP7_PAS"              ,    100 , 0.0     , 100.    );  CreateUserTH1D( "Jets_Bar_TrkIsoHEEP7_PAS"              ,    100 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 100 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End1_nElectron_PAS"             ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End1_nElectron_PAS"              ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End1_nElectron_PAS"              ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End1_Pt1stEle_PAS"              ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End1_Pt1stEle_PAS"	           ,    1000, 0       , 1000       );  CreateUserTH1D( "Jets_End1_Pt1stEle_PAS"	           ,    1000, 0       , 1000       );
  CreateUserTH1D( "Total_End1_hltPt1stEle_PAS"           ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End1_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );  CreateUserTH1D( "Jets_End1_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );
  CreateUserTH1D( "Total_End1_Eta1stEle_PAS"	          ,    100 , -5      , 5	       );  CreateUserTH1D( "Electrons_End1_Eta1stEle_PAS"	           ,    100 , -5      , 5          );  CreateUserTH1D( "Jets_End1_Eta1stEle_PAS"	           ,    100 , -5      , 5          );
  CreateUserTH1D( "Total_End1_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416   );  CreateUserTH1D( "Electrons_End1_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416    );  CreateUserTH1D( "Jets_End1_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416    );
  CreateUserTH1D( "Total_End1_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End1_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End1_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End1_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End1_DCotTheta1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End1_DCotTheta1stEle_PAS"        ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End1_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End1_Dist1stEle_PAS"             ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End1_Dist1stEle_PAS"             ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End1_MET_PAS"                   ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End1_MET_PAS"                    ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End1_MET_PAS"                    ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End1_TrkIsoHEEP7_PAS"            ,    100 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End1_TrkIsoHEEP7_PAS"              ,    100 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End1_TrkIsoHEEP7_PAS"              ,    100 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End1_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 100 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End2_nElectron_PAS"             ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End2_nElectron_PAS"              ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End2_nElectron_PAS"              ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End2_Pt1stEle_PAS"              ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End2_Pt1stEle_PAS"	           ,    1000, 0       , 1000       );  CreateUserTH1D( "Jets_End2_Pt1stEle_PAS"	           ,    1000, 0       , 1000       );
  CreateUserTH1D( "Total_End2_hltPt1stEle_PAS"           ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End2_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );  CreateUserTH1D( "Jets_End2_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );
  CreateUserTH1D( "Total_End2_Eta1stEle_PAS"	          ,    100 , -5      , 5	       );  CreateUserTH1D( "Electrons_End2_Eta1stEle_PAS"	           ,    100 , -5      , 5          );  CreateUserTH1D( "Jets_End2_Eta1stEle_PAS"	           ,    100 , -5      , 5          );
  CreateUserTH1D( "Total_End2_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416   );  CreateUserTH1D( "Electrons_End2_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416    );  CreateUserTH1D( "Jets_End2_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416    );
  CreateUserTH1D( "Total_End2_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End2_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End2_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End2_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End2_DCotTheta1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End2_DCotTheta1stEle_PAS"        ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End2_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End2_Dist1stEle_PAS"             ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End2_Dist1stEle_PAS"             ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End2_MET_PAS"                   ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End2_MET_PAS"                    ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End2_MET_PAS"                    ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End2_TrkIsoHEEP7_PAS"            ,    100 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End2_TrkIsoHEEP7_PAS"              ,    100 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End2_TrkIsoHEEP7_PAS"              ,    100 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End2_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 100 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );
                                                                                                                                                                                                                                                                                     
  // 1 jet                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_Bar_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_Bar_1Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_Bar_1Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_Bar_1Jet_Pt1stEle_PAS"          ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Bar_1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_Bar_1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_Bar_1Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5             );  CreateUserTH1D( "Electrons_Bar_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5              );  CreateUserTH1D( "Jets_Bar_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5              );
  CreateUserTH1D( "Total_Bar_1Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_Bar_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416        );  CreateUserTH1D( "Jets_Bar_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416        );
  CreateUserTH1D( "Total_Bar_1Jet_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001      );  CreateUserTH1D( "Electrons_Bar_1Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_Bar_1Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_Bar_1Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_Bar_1Jet_DCotTheta1stEle_PAS"    ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_Bar_1Jet_DCotTheta1stEle_PAS"    ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_Bar_1Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_Bar_1Jet_Dist1stEle_PAS"         ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_Bar_1Jet_Dist1stEle_PAS"         ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_Bar_1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_Bar_1Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_Bar_1Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_Bar_1Jet_TrkIsoHEEP7_PAS"       ,    100 , 0.0     , 100.     );  CreateUserTH1D( "Electrons_Bar_1Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.    );  CreateUserTH1D( "Jets_Bar_1Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 100 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End1_1Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End1_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End1_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End1_1Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End1_1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_End1_1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_End1_1Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_End1_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );  CreateUserTH1D( "Jets_End1_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );
  CreateUserTH1D( "Total_End1_1Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_End1_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End1_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End1_1Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End1_1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End1_1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End1_1Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End1_1Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End1_1Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End1_1Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End1_1Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End1_1Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End1_1Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End1_1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End1_1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End1_1Jet_TrkIsoHEEP7_PAS"       ,    100 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End1_1Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End1_1Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 100 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End2_1Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End2_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End2_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End2_1Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End2_1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_End2_1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_End2_1Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_End2_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );  CreateUserTH1D( "Jets_End2_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );
  CreateUserTH1D( "Total_End2_1Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_End2_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End2_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End2_1Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End2_1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End2_1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End2_1Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End2_1Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End2_1Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End2_1Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End2_1Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End2_1Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End2_1Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End2_1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End2_1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End2_1Jet_TrkIsoHEEP7_PAS"       ,    100 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End2_1Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End2_1Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End2_1Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 100 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );
                                                                                                                                                                                                                                                                                     
  // 2 jet                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_Bar_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_Bar_2Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_Bar_2Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_Bar_2Jet_Pt1stEle_PAS"          ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Bar_2Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_Bar_2Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_Bar_2Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_Bar_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5              );  CreateUserTH1D( "Jets_Bar_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5              );
  CreateUserTH1D( "Total_Bar_2Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_Bar_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416        );  CreateUserTH1D( "Jets_Bar_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416        );
  CreateUserTH1D( "Total_Bar_2Jet_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001      );  CreateUserTH1D( "Electrons_Bar_2Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_Bar_2Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_Bar_2Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_Bar_2Jet_DCotTheta1stEle_PAS"    ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_Bar_2Jet_DCotTheta1stEle_PAS"    ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_Bar_2Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_Bar_2Jet_Dist1stEle_PAS"         ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_Bar_2Jet_Dist1stEle_PAS"         ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_Bar_2Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_Bar_2Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_Bar_2Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_Bar_2Jet_TrkIsoHEEP7_PAS"       ,    100 , 0.0     , 100.     );  CreateUserTH1D( "Electrons_Bar_2Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.    );  CreateUserTH1D( "Jets_Bar_2Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 100 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_2Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_2Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End1_2Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End1_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End1_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End1_2Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End1_2Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_End1_2Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_End1_2Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_End1_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );  CreateUserTH1D( "Jets_End1_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );
  CreateUserTH1D( "Total_End1_2Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_End1_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End1_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End1_2Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End1_2Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End1_2Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End1_2Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End1_2Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End1_2Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End1_2Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End1_2Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End1_2Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End1_2Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End1_2Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End1_2Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End1_2Jet_TrkIsoHEEP7_PAS"       ,    100 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End1_2Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End1_2Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 100 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_2Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_2Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End2_2Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End2_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End2_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End2_2Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End2_2Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_End2_2Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_End2_2Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_End2_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );  CreateUserTH1D( "Jets_End2_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );
  CreateUserTH1D( "Total_End2_2Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_End2_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End2_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End2_2Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End2_2Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End2_2Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End2_2Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End2_2Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End2_2Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End2_2Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End2_2Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End2_2Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End2_2Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End2_2Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End2_2Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End2_2Jet_TrkIsoHEEP7_PAS"       ,    100 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End2_2Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End2_2Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End2_2Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 100 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_2Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_2Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );
                                                                                                                                                                                                                                                                                     
  // 3 jet                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_Bar_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_Bar_3Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_Bar_3Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_Bar_3Jet_Pt1stEle_PAS"          ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Bar_3Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_Bar_3Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_Bar_3Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_Bar_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5              );  CreateUserTH1D( "Jets_Bar_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5              );
  CreateUserTH1D( "Total_Bar_3Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_Bar_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416        );  CreateUserTH1D( "Jets_Bar_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416        );
  CreateUserTH1D( "Total_Bar_3Jet_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001      );  CreateUserTH1D( "Electrons_Bar_3Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_Bar_3Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_Bar_3Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_Bar_3Jet_DCotTheta1stEle_PAS"    ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_Bar_3Jet_DCotTheta1stEle_PAS"    ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_Bar_3Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_Bar_3Jet_Dist1stEle_PAS"         ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_Bar_3Jet_Dist1stEle_PAS"         ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_Bar_3Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_Bar_3Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_Bar_3Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_Bar_3Jet_TrkIsoHEEP7_PAS"       ,    100 , 0.0     , 100.     );  CreateUserTH1D( "Electrons_Bar_3Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.    );  CreateUserTH1D( "Jets_Bar_3Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 100 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_3Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_3Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End1_3Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End1_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End1_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End1_3Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End1_3Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_End1_3Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_End1_3Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_End1_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );  CreateUserTH1D( "Jets_End1_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );
  CreateUserTH1D( "Total_End1_3Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_End1_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End1_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End1_3Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End1_3Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End1_3Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End1_3Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End1_3Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End1_3Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End1_3Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End1_3Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End1_3Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End1_3Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End1_3Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End1_3Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End1_3Jet_TrkIsoHEEP7_PAS"       ,    100 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End1_3Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End1_3Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 100 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_3Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_3Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End2_3Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End2_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End2_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End2_3Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End2_3Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_End2_3Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_End2_3Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_End2_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );  CreateUserTH1D( "Jets_End2_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );
  CreateUserTH1D( "Total_End2_3Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_End2_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End2_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End2_3Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End2_3Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End2_3Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End2_3Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End2_3Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End2_3Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End2_3Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End2_3Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End2_3Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End2_3Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End2_3Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End2_3Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End2_3Jet_TrkIsoHEEP7_PAS"       ,    100 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End2_3Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End2_3Jet_TrkIsoHEEP7_PAS"         ,    100 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End2_3Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 100 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_3Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_3Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     100 , 0.0     , 100.    );
                                                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                     
  // Pile-up [0,5]                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_Bar_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_Bar_PU0-8_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_Bar_PU0-8_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_Bar_PU0-8_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Bar_PU0-8_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_Bar_PU0-8_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_Bar_PU0-8_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_Bar_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );  CreateUserTH1D( "Jets_Bar_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );
  CreateUserTH1D( "Total_Bar_PU0-8_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_Bar_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_Bar_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_Bar_PU0-8_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001    );  CreateUserTH1D( "Electrons_Bar_PU0-8_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_Bar_PU0-8_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_Bar_PU0-8_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_Bar_PU0-8_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_Bar_PU0-8_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_Bar_PU0-8_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_Bar_PU0-8_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_Bar_PU0-8_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_Bar_PU0-8_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_Bar_PU0-8_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_Bar_PU0-8_MET_PAS"               ,    100 , 0.0     , 1000.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End1_PU0-8_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End1_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End1_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End1_PU0-8_Pt1stEle_PAS"        ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End1_PU0-8_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );  CreateUserTH1D( "Jets_End1_PU0-8_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );
  CreateUserTH1D( "Total_End1_PU0-8_Eta1stEle_PAS"	  ,    100 , -5      , 5	         );  CreateUserTH1D( "Electrons_End1_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );  CreateUserTH1D( "Jets_End1_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );
  CreateUserTH1D( "Total_End1_PU0-8_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416     );  CreateUserTH1D( "Electrons_End1_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End1_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End1_PU0-8_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001    );  CreateUserTH1D( "Electrons_End1_PU0-8_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End1_PU0-8_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End1_PU0-8_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End1_PU0-8_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End1_PU0-8_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End1_PU0-8_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End1_PU0-8_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End1_PU0-8_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End1_PU0-8_MET_PAS"             ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End1_PU0-8_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End1_PU0-8_MET_PAS"              ,    100 , 0.0     , 1000.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End2_PU0-8_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End2_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End2_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End2_PU0-8_Pt1stEle_PAS"        ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End2_PU0-8_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );  CreateUserTH1D( "Jets_End2_PU0-8_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );
  CreateUserTH1D( "Total_End2_PU0-8_Eta1stEle_PAS"	  ,    100 , -5      , 5	         );  CreateUserTH1D( "Electrons_End2_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );  CreateUserTH1D( "Jets_End2_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );
  CreateUserTH1D( "Total_End2_PU0-8_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416     );  CreateUserTH1D( "Electrons_End2_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End2_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End2_PU0-8_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001    );  CreateUserTH1D( "Electrons_End2_PU0-8_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End2_PU0-8_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End2_PU0-8_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End2_PU0-8_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End2_PU0-8_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End2_PU0-8_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End2_PU0-8_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End2_PU0-8_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End2_PU0-8_MET_PAS"             ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End2_PU0-8_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End2_PU0-8_MET_PAS"              ,    100 , 0.0     , 1000.    );
                                                                                                                                                                                                                                                                                     
  // Pile-up [6-10]                                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_Bar_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_Bar_PU9-UP_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_Bar_PU9-UP_nElectron_PAS"        ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_Bar_PU9-UP_Pt1stEle_PAS"        ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Bar_PU9-UP_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );  CreateUserTH1D( "Jets_Bar_PU9-UP_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );
  CreateUserTH1D( "Total_Bar_PU9-UP_Eta1stEle_PAS"	  ,    100 , -5      , 5	         );  CreateUserTH1D( "Electrons_Bar_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );  CreateUserTH1D( "Jets_Bar_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );
  CreateUserTH1D( "Total_Bar_PU9-UP_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416     );  CreateUserTH1D( "Electrons_Bar_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_Bar_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_Bar_PU9-UP_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001    );  CreateUserTH1D( "Electrons_Bar_PU9-UP_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_Bar_PU9-UP_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_Bar_PU9-UP_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_Bar_PU9-UP_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_Bar_PU9-UP_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_Bar_PU9-UP_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_Bar_PU9-UP_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_Bar_PU9-UP_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_Bar_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_Bar_PU9-UP_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_Bar_PU9-UP_MET_PAS"              ,    100 , 0.0     , 1000.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End1_PU9-UP_nElectron_PAS"      ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End1_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End1_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End1_PU9-UP_Pt1stEle_PAS"       ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End1_PU9-UP_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );  CreateUserTH1D( "Jets_End1_PU9-UP_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );
  CreateUserTH1D( "Total_End1_PU9-UP_Eta1stEle_PAS"	  ,    100 , -5      , 5	         );  CreateUserTH1D( "Electrons_End1_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	         );  CreateUserTH1D( "Jets_End1_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	         );
  CreateUserTH1D( "Total_End1_PU9-UP_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416     );  CreateUserTH1D( "Electrons_End1_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416    );  CreateUserTH1D( "Jets_End1_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416    );
  CreateUserTH1D( "Total_End1_PU9-UP_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001  );  CreateUserTH1D( "Electrons_End1_PU9-UP_Charge1stEle_PAS"    ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End1_PU9-UP_Charge1stEle_PAS"    ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End1_PU9-UP_DCotTheta1stEle_PAS",    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End1_PU9-UP_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End1_PU9-UP_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End1_PU9-UP_Dist1stEle_PAS"     ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End1_PU9-UP_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End1_PU9-UP_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End1_PU9-UP_MET_PAS"            ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End1_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End1_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1000.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End2_PU9-UP_nElectron_PAS"      ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End2_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End2_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End2_PU9-UP_Pt1stEle_PAS"       ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End2_PU9-UP_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );  CreateUserTH1D( "Jets_End2_PU9-UP_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );
  CreateUserTH1D( "Total_End2_PU9-UP_Eta1stEle_PAS"	  ,    100 , -5      , 5	         );  CreateUserTH1D( "Electrons_End2_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	         );  CreateUserTH1D( "Jets_End2_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	         );
  CreateUserTH1D( "Total_End2_PU9-UP_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416     );  CreateUserTH1D( "Electrons_End2_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416    );  CreateUserTH1D( "Jets_End2_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416    );
  CreateUserTH1D( "Total_End2_PU9-UP_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001  );  CreateUserTH1D( "Electrons_End2_PU9-UP_Charge1stEle_PAS"    ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End2_PU9-UP_Charge1stEle_PAS"    ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End2_PU9-UP_DCotTheta1stEle_PAS",    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End2_PU9-UP_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End2_PU9-UP_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End2_PU9-UP_Dist1stEle_PAS"     ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Electrons_End2_PU9-UP_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Jets_End2_PU9-UP_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );
  CreateUserTH1D( "Total_End2_PU9-UP_MET_PAS"            ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End2_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End2_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1000.    );

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
    // Check good run list
    //--------------------------------------------------------------------------

    int    passedJSON = passJSON ( int(run), int(ls) , int(isData) ) ;

    //--------------------------------------------------------------------------
    // Do pileup re-weighting
    //--------------------------------------------------------------------------

    double pileup_weight = getPileupWeight ( nPileUpInt_True, isData ) ;

    //--------------------------------------------------------------------------
    // Fill variables
    //--------------------------------------------------------------------------

    short min_prescale = 0;
    int passTrigger = 0;

    if ( isData ) {
      if ( LooseEle1_hltPhotonPt > 0.0 ) {  // enforce trigger object matching
        if ( H_Photon22   > 0.1 && LooseEle1_hltPhotonPt >= 22.  && LooseEle1_hltPhotonPt < 30. ) { passTrigger = 1; min_prescale = H_Photon22  ; } 
        if ( H_Photon30   > 0.1 && LooseEle1_hltPhotonPt >= 30.  && LooseEle1_hltPhotonPt < 36. ) { passTrigger = 1; min_prescale = H_Photon30  ; } 
        if ( H_Photon36   > 0.1 && LooseEle1_hltPhotonPt >= 36.  && LooseEle1_hltPhotonPt < 50. ) { passTrigger = 1; min_prescale = H_Photon36  ; } 
        if ( H_Photon50   > 0.1 && LooseEle1_hltPhotonPt >= 50.  && LooseEle1_hltPhotonPt < 75. ) { passTrigger = 1; min_prescale = H_Photon50  ; } 
        if ( H_Photon75   > 0.1 && LooseEle1_hltPhotonPt >= 75.  && LooseEle1_hltPhotonPt < 90. ) { passTrigger = 1; min_prescale = H_Photon75  ; } 
        if ( H_Photon90   > 0.1 && LooseEle1_hltPhotonPt >= 90.  && LooseEle1_hltPhotonPt < 120.) { passTrigger = 1; min_prescale = H_Photon90  ; } 
        if ( H_Photon120  > 0.1 && LooseEle1_hltPhotonPt >= 120. && LooseEle1_hltPhotonPt < 175.) { passTrigger = 1; min_prescale = H_Photon120 ; } 
        if ( H_Photon175  > 0.1 && LooseEle1_hltPhotonPt >= 175.)                                 { passTrigger = 1; min_prescale = H_Photon175 ; } 
      }
    }  // end if (isData) 

    else { 
      min_prescale = 1;
      passTrigger = 1 ; //FIXME TODO: would need proper trigger handling for MC
    }

    //--------------------------------------------------------------------------
    // Fill cut values
    //--------------------------------------------------------------------------

    // No cut
    fillVariableWithValue(   "Weighting"                     , 1                           , min_prescale * pileup_weight ); 

    // Trigger 
    fillVariableWithValue(   "PassTrigger"                   , passTrigger                 , min_prescale * pileup_weight ); 

    // JSON variable								            
    fillVariableWithValue(   "PassJSON"                      , passedJSON                  , min_prescale * pileup_weight ); 

    // Fill noise filters
    // see: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    // we filled these at skim time
    fillVariableWithValue( "PassGlobalTightHalo2016Filter" , PassGlobalTightHalo2016Filter  , min_prescale * pileup_weight );
    fillVariableWithValue( "PassGoodVertices"	             , PassGoodVertices               , min_prescale * pileup_weight );
    fillVariableWithValue( "PassHBHENoiseFilter"	         , PassHBHENoiseFilter            , min_prescale * pileup_weight );
    fillVariableWithValue( "PassHBHENoiseIsoFilter"	       , PassHBHENoiseIsoFilter         , min_prescale * pileup_weight );
    fillVariableWithValue( "PassBadEESupercrystalFilter"   , PassBadEESupercrystalFilter    , min_prescale * pileup_weight );
    fillVariableWithValue( "PassEcalDeadCellTrigPrim"      , PassEcalDeadCellTrigPrim       , min_prescale * pileup_weight );
    fillVariableWithValue( "PassChargedCandidateFilter"    , PassChargedCandidateFilter     , min_prescale * pileup_weight );
    fillVariableWithValue( "PassBadPFMuonFilter"           , PassBadPFMuonFilter            , min_prescale * pileup_weight );

    // Muon variables ( for veto ) 					      	       
    // remove muon veto
    //fillVariableWithValue(   "nMuon"                         , nMuon_ptCut                 , min_prescale * pileup_weight );

    // 1st Electron variables				      		                    
    fillVariableWithValue(   "nEle"                          , nLooseEle_ptCut             , min_prescale * pileup_weight );

    fillVariableWithValue(   "Pt1stEle"                      , LooseEle1_PtHeep            , min_prescale * pileup_weight );
    fillVariableWithValue(   "Eta1stEle"                     , LooseEle1_Eta               , min_prescale * pileup_weight );
    fillVariableWithValue(   "Phi1stEle"                     , LooseEle1_Phi               , min_prescale * pileup_weight );
    fillVariableWithValue(   "HLTPt1stEle"                   , LooseEle1_hltPhotonPt       , min_prescale * pileup_weight );
    fillVariableWithValue(   "H_Photon22"              , H_Photon22            , min_prescale * pileup_weight );
    fillVariableWithValue(   "H_Photon30"              , H_Photon30            , min_prescale * pileup_weight );
    fillVariableWithValue(   "H_Photon36"              , H_Photon36            , min_prescale * pileup_weight );
    fillVariableWithValue(   "H_Photon50"              , H_Photon50            , min_prescale * pileup_weight );   
    fillVariableWithValue(   "H_Photon75"              , H_Photon75            , min_prescale * pileup_weight );   
    fillVariableWithValue(   "H_Photon90"              , H_Photon90            , min_prescale * pileup_weight );   
    fillVariableWithValue(   "H_Photon120"             , H_Photon120           , min_prescale * pileup_weight );   
    fillVariableWithValue(   "H_Photon175"             , H_Photon175           , min_prescale * pileup_weight );    

    // 1st JET variables                                     		                    
    fillVariableWithValue(   "nJet"                          , nJetLooseEle_ptCut          , min_prescale * pileup_weight );

    // MET									            
    //fillVariableWithValue(   "MET"                           , PFMET_Type1_Pt           , min_prescale * pileup_weight );

    // Dummy variables
    // NB: QCD reduced skims contain loose electrons which have passed QCD loose FR ID 
    fillVariableWithValue(   "denominator"                   , 1                           , min_prescale * pileup_weight );

    // Debugging variables 
    // fillVariableWithValue ( "HLTPt1stEle"        , LooseEle1_hltPhotonPt , min_prescale * pileup_weight );
    // fillVariableWithValue ( "H_Photon30_CIdVL"  , H_Photon30_CIdVL      , min_prescale * pileup_weight );	      
    // fillVariableWithValue ( "H_Photon50_CIdVL"  , H_Photon50_CIdVL      , min_prescale * pileup_weight );		
    // fillVariableWithValue ( "H_Photon75_CIdVL"  , H_Photon75_CIdVL      , min_prescale * pileup_weight );		
    // fillVariableWithValue ( "H_Photon90_CIdVL"  , H_Photon90_CIdVL      , min_prescale * pileup_weight );		
    // fillVariableWithValue ( "H_Photon135"	 , H_Photon135	         , min_prescale * pileup_weight );		
    // fillVariableWithValue ( "H_Photon150"	 , H_Photon150	         , min_prescale * pileup_weight );		
    // fillVariableWithValue ( "H_Photon160"	 , H_Photon160	         , min_prescale * pileup_weight );		
    // fillVariableWithValue ( "MinPrescale"	 , min_prescale          , min_prescale * pileup_weight );
    // fillVariableWithValue ( "run"               , run                   , min_prescale * pileup_weight );
    // fillVariableWithValue ( "event"             , event                 , min_prescale * pileup_weight );
    // fillVariableWithValue ( "ls"                , ls                    , min_prescale * pileup_weight );

    //--------------------------------------------------------------------------
    // Evaluate the cuts
    //--------------------------------------------------------------------------

    if(run==284043)
    {
      std::cout << (uint64_t)run << " " << (uint64_t)ls << " " << (uint64_t)event << std::endl;
      evaluateCuts(true);
    }
    else
      evaluateCuts(false);

    //--------------------------------------------------------------------------
    // Fill preselection plots
    //--------------------------------------------------------------------------

    bool passed_denominator = passedAllPreviousCuts("denominator");

    bool isBarrel  = false;
    bool isEndcap1 = false;
    bool isEndcap2 = false;

    if( fabs( LooseEle1_Eta  ) < eleEta_bar )        isBarrel = true;
    if( fabs( LooseEle1_Eta  ) > eleEta_end1_min &&
        fabs( LooseEle1_Eta  ) < eleEta_end1_max )   isEndcap1 = true;
    if( fabs( LooseEle1_Eta  ) > eleEta_end2_min &&
        fabs( LooseEle1_Eta  ) < eleEta_end2_max )   isEndcap2 = true;

    if ( passed_denominator ) { 

      // debugging 
      // fillReducedSkimTree();

      FillUserTH1D("Total_minPrescale"                     ,   min_prescale );
      FillUserTH1D("Total_pileupWeight"                    ,   pileup_weight);

      if ( isBarrel ) {

        FillUserTH1D( "Total_Bar_nVertex_PAS"           , nVertex                      , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_Pt1stEle_PAS"	         , LooseEle1_Pt                 , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_hltPt1stEle_PAS"       , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_Eta1stEle_PAS"	 , LooseEle1_Eta                , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_Phi1stEle_PAS"	 , LooseEle1_Phi                , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
        FillUserTH1D( "Total_Bar_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        FillUserTH2D( "Total_Bar_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);

        if ( nJetLooseEle_ptCut >= 1 ) {
          FillUserTH1D( "Total_Bar_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_1Jet_Pt1stEle_PAS"	      , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_1Jet_Eta1stEle_PAS"	      , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_1Jet_Phi1stEle_PAS"	      , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_Bar_1Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_1Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }

        if ( nJetLooseEle_ptCut >= 2 ) {
          FillUserTH1D( "Total_Bar_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_2Jet_Pt1stEle_PAS"	      , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_2Jet_Eta1stEle_PAS"	      , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_2Jet_Phi1stEle_PAS"	      , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_Bar_2Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_2Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }

        if ( nJetLooseEle_ptCut >= 3 ) {
          FillUserTH1D( "Total_Bar_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_3Jet_Pt1stEle_PAS"	      , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_3Jet_Eta1stEle_PAS"	      , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_3Jet_Phi1stEle_PAS"	      , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_Bar_3Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_3Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }

        if ( nVertex >= 0 && nVertex <= 5 ){
          FillUserTH1D( "Total_Bar_PU0-8_nElectron_PAS"          , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU0-8_Pt1stEle_PAS"           , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU0-8_Eta1stEle_PAS"          , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU0-8_Phi1stEle_PAS"          , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU0-8_Charge1stEle_PAS"       , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_Bar_PU0-8_DCotTheta1stEle_PAS"    , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU0-8_Dist1stEle_PAS"         , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU0-8_MET_PAS"                , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        }						          

        if ( nVertex >= 6 && nVertex <= 10 ){      
          FillUserTH1D( "Total_Bar_PU9-UP_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU9-UP_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU9-UP_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU9-UP_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU9-UP_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_Bar_PU9-UP_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU9-UP_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU9-UP_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        }
      }

      if ( isEndcap1 ) {
        FillUserTH1D( "Total_End1_nVertex_PAS"           , nVertex                      , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End1_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End1_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End1_hltPt1stEle_PAS"       , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End1_Eta1stEle_PAS"	  , LooseEle1_Eta                , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End1_Phi1stEle_PAS"	  , LooseEle1_Phi                , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End1_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
        FillUserTH1D( "Total_End1_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End1_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End1_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End1_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        FillUserTH2D( "Total_End1_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);


        if ( nJetLooseEle_ptCut >= 1 ) {
          FillUserTH1D( "Total_End1_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_1Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End1_1Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_1Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }

        if ( nJetLooseEle_ptCut >= 2 ) {
          FillUserTH1D( "Total_End1_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_2Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End1_2Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_2Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }

        if ( nJetLooseEle_ptCut >= 3 ) {
          FillUserTH1D( "Total_End1_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_3Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End1_3Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_3Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }

        if ( nVertex >= 0 && nVertex <= 5 ){
          FillUserTH1D( "Total_End1_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End1_PU0-8_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU0-8_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU0-8_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        }						          

        if ( nVertex >= 6 && nVertex <= 10 ){      
          FillUserTH1D( "Total_End1_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End1_PU9-UP_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU9-UP_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU9-UP_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        }
      }

      if ( isEndcap2 ) {
        FillUserTH1D( "Total_End2_nVertex_PAS"           , nVertex                      , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End2_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End2_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End2_hltPt1stEle_PAS"       , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End2_Eta1stEle_PAS"	  , LooseEle1_Eta                , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End2_Phi1stEle_PAS"	  , LooseEle1_Phi                , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End2_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
        FillUserTH1D( "Total_End2_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End2_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End2_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End2_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        FillUserTH2D( "Total_End2_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);

        if ( nJetLooseEle_ptCut >= 1 ) {
          FillUserTH1D( "Total_End2_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_1Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End2_1Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_1Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End2_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }

        if ( nJetLooseEle_ptCut >= 2 ) {
          FillUserTH1D( "Total_End2_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_2Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End2_2Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_2Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End2_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }

        if ( nJetLooseEle_ptCut >= 3 ) {
          FillUserTH1D( "Total_End2_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_3Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End2_3Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_3Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End2_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }

        if ( nVertex >= 0 && nVertex <= 5 ){
          FillUserTH1D( "Total_End2_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End2_PU0-8_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU0-8_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU0-8_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        }						          

        if ( nVertex >= 6 && nVertex <= 10 ){      
          FillUserTH1D( "Total_End2_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End2_PU9-UP_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU9-UP_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU9-UP_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        }
      }

      //----------------------------------------------------------------------
      // HEEP' 7.0 ID part
      //  same as HEEPv7.0, but no TrkIso cut
      //----------------------------------------------------------------------
      //  Bools that are the same whether barrel or endcap
      bool pass_et            = bool ( LooseEle1_PtHeep              >  35.0 );
      bool pass_ecalDriven    = bool ( LooseEle1_EcalDriven        == 1    );
      bool pass_deltaPhi      = bool ( fabs (LooseEle1_DeltaPhiTrkSC) <  0.06 ); // dPhiSCTrkAtVtx
      bool pass_missingHits   = bool ( LooseEle1_MissingHits     <= 1    );

      // Bools that depend on barrel vs. endcap
      bool pass_deltaEtaSeed  = false;
      bool pass_sigmaIEtaIEta = false;
      bool pass_shape         = false;
      bool pass_shape1        = false;
      bool pass_shape2        = false;
      bool pass_caloIsolation = false;
      bool pass_dxy           = false;
      bool pass_hoe           = false;

      double caloIsolation = LooseEle1_EcalIsolation + LooseEle1_HcalIsolation;

      // Barrel electrons
      if ( fabs(LooseEle1_SCEta) < 1.4442 ){
        pass_deltaEtaSeed      = bool ( fabs(LooseEle1_DeltaEtaSeed )     < 0.004 );
        pass_hoe               = bool ( LooseEle1_HoE            < 1/LooseEle1_SCEnergy + 0.05 );
        pass_sigmaIEtaIEta     = true;
        pass_shape1            = bool ( LooseEle1_Full5x5E1x5OverE5x5        > 0.83  );
        pass_shape2            = bool ( LooseEle1_Full5x5E2x5OverE5x5        > 0.94  );
        pass_shape             = bool ( pass_shape1 || pass_shape2    );
        pass_caloIsolation     = bool ( caloIsolation < ( 2.0 + ( 0.03 * LooseEle1_PtHeep ) + (0.28 * LooseEle1_RhoForHeep ) ) );
        pass_dxy               = bool ( fabs(LooseEle1_LeadVtxDistXY) < 0.02  );
      }
      // Endcap electrons
      else if ( fabs(LooseEle1_SCEta) > 1.566 && fabs(LooseEle1_SCEta) < 2.5 ){ 
        pass_deltaEtaSeed      = bool ( fabs (LooseEle1_DeltaEtaSeed)     < 0.006 );
        pass_hoe               = bool ( LooseEle1_HoE            < 5/LooseEle1_SCEnergy + 0.05 );
        pass_sigmaIEtaIEta     = bool ( LooseEle1_Full5x5SigmaIEtaIEta       < 0.03  );
        pass_shape             = true;
        if   ( LooseEle1_PtHeep  < 50 ) {
          pass_caloIsolation = bool ( caloIsolation < ( 2.5 + 
                ( 0.28 * LooseEle1_RhoForHeep ) ) );
        }
        else                { 
          pass_caloIsolation = bool ( caloIsolation < ( 2.5 + 
                ( 0.28 * LooseEle1_RhoForHeep ) + 
                ( 0.03 * (LooseEle1_PtHeep - 50.0 ) ) ) );
        }
        pass_dxy               = bool ( fabs(LooseEle1_LeadVtxDistXY) < 0.05  );
      }
      // final decision
      bool passHEEPprime = (
          pass_et            && 
          pass_ecalDriven    && 
          pass_deltaEtaSeed  && 
          pass_deltaPhi      && 
          pass_hoe           && 
          pass_sigmaIEtaIEta && 
          pass_shape         && 
          pass_dxy           && 
          pass_missingHits   && 
          pass_caloIsolation ); 

      //XXX SIC TEST
      if(run==284043 && (ls==116 || ls==146 || ls==172 || ls==199 || ls==45 || ls==8) )
      {
              //std::cout << (uint64_t)run << " " << (uint64_t)ls << " " << (uint64_t)event << std::endl;
              std::cout << "\tPassEt: " << pass_et << " --> " << LooseEle1_PtHeep << std::endl
                << "\tPassEcalDriven: " << pass_ecalDriven << " --> " << LooseEle1_EcalDriven << std::endl
                << "\tPassDeltaPhi: " << pass_deltaPhi << " --> " << LooseEle1_DeltaPhiTrkSC << std::endl
                << "\tPassMissingHits: " << pass_missingHits << " --> " << LooseEle1_MissingHits << std::endl
                << "\tPassDeltaEtaSeed: " << pass_deltaEtaSeed << " --> " << LooseEle1_DeltaEtaSeed << std::endl;
              if ( fabs(LooseEle1_SCEta) > 1.566 && fabs(LooseEle1_SCEta) < 2.5 )
                std::cout << "\tPassSigmaIetaIeta: " << pass_sigmaIEtaIEta << " --> " << LooseEle1_Full5x5SigmaIEtaIEta << " [ENDCAP]"<< std::endl;
              else
                std::cout << "\tPassSigmaIetaIeta: " << pass_sigmaIEtaIEta << " --> " << LooseEle1_Full5x5SigmaIEtaIEta << " [BARREL]"<< std::endl;
              if ( fabs(LooseEle1_SCEta) < 1.4442 )
              {
                std::cout << "\tPassShowerShape: " << pass_shape << " --> E1x5/E5x5=" << LooseEle1_Full5x5E1x5OverE5x5 << ", E2x5/E5x5=" << LooseEle1_Full5x5E2x5OverE5x5 << " [BARREL]" << std::endl;
                std::cout << "\tPassCaloIso: " << pass_caloIsolation << " --> " << caloIsolation << " < " << ( 2.0 + ( 0.03 * LooseEle1_PtHeep ) + (0.28 * LooseEle1_RhoForHeep ) ) << std::endl;
                std::cout << "\tPassDxy: " << pass_dxy << " --> " << LooseEle1_LeadVtxDistXY << " < 0.02 " << std::endl
                  << "\tPassH/E: " << pass_hoe << " --> " << LooseEle1_HoE << " < " << 1/LooseEle1_SCEnergy + 0.05 << std::endl;
              }
              else
              {
                std::cout << "\tPassCaloIso: " << pass_caloIsolation << " --> " << caloIsolation << " < " << (LooseEle1_PtHeep<50 ? (2.5+( 0.28 * LooseEle1_RhoForHeep ) ) : (2.5+( 0.28 * LooseEle1_RhoForHeep ) + ( 0.03 * (LooseEle1_PtHeep - 50.0 ) )) ) << std::endl;
                std::cout << "\tPassDxy: " << pass_dxy << " --> " << LooseEle1_LeadVtxDistXY << " < 0.05" << std::endl
                  << "\tPassH/E: " << pass_hoe << " --> " << LooseEle1_HoE << " < " << 5/LooseEle1_SCEnergy + 0.05 << std::endl;
              }
              std::cout << "\tPassTrkIso: " << (LooseEle1_TrkIsoHEEP7 < 5.0) << " --> " << LooseEle1_TrkIsoHEEP7 << " < 5.0" << std::endl;
      }
      //XXX SIC TEST

      if ( passHEEPprime ) { 

        if ( isBarrel ) { 
          FillUserTH1D( "Electrons_Bar_nVertex_PAS"           , nVertex                      , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_Pt1stEle_PAS"	  , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_hltPt1stEle_PAS"	  , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_Eta1stEle_PAS"	  , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_Phi1stEle_PAS"	  , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Electrons_Bar_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);

          if ( nJetLooseEle_ptCut >= 1 ) {
            FillUserTH1D( "Electrons_Bar_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_1Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_1Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_1Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_Bar_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut >= 2 ) {
            FillUserTH1D( "Electrons_Bar_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_2Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_2Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_2Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_Bar_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut >= 3 ) {
            FillUserTH1D( "Electrons_Bar_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_3Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_3Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_3Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_Bar_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nVertex >= 0 && nVertex <= 5 ){
            FillUserTH1D( "Electrons_Bar_PU0-8_nElectron_PAS"          , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU0-8_Pt1stEle_PAS"           , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU0-8_Eta1stEle_PAS"          , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU0-8_Phi1stEle_PAS"          , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU0-8_Charge1stEle_PAS"       , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Electrons_Bar_PU0-8_DCotTheta1stEle_PAS"    , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU0-8_Dist1stEle_PAS"         , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU0-8_MET_PAS"                , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }						          

          if ( nVertex >= 6 && nVertex <= 10 ){      
            FillUserTH1D( "Electrons_Bar_PU9-UP_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU9-UP_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU9-UP_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU9-UP_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU9-UP_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Electrons_Bar_PU9-UP_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU9-UP_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU9-UP_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }	  
        }

        if ( isEndcap1 ) { 
          FillUserTH1D( "Electrons_End1_nVertex_PAS"          , nVertex                      , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End1_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End1_hltPt1stEle_PAS"	  , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End1_Pt1stEle_PAS"	  , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End1_Eta1stEle_PAS"	  , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End1_Phi1stEle_PAS"	  , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End1_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End1_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End1_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End1_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End1_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Electrons_End1_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);

          if ( nJetLooseEle_ptCut >= 1 ) {
            FillUserTH1D( "Electrons_End1_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_1Jet_Pt1stEle_PAS"        , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_1Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_1Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End1_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);

          }

          if ( nJetLooseEle_ptCut >= 2 ) {
            FillUserTH1D( "Electrons_End1_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_2Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_2Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_2Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End1_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut >= 3 ) {
            FillUserTH1D( "Electrons_End1_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_3Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_3Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_3Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End1_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nVertex >= 0 && nVertex <= 5 ){
            FillUserTH1D( "Electrons_End1_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Electrons_End1_PU0-8_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU0-8_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU0-8_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }						          

          if ( nVertex >= 6 && nVertex <= 10 ){      
            FillUserTH1D( "Electrons_End1_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Electrons_End1_PU9-UP_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU9-UP_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU9-UP_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }
        }

        if ( isEndcap2 ) { 
          FillUserTH1D( "Electrons_End2_nVertex_PAS"          , nVertex                      , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End2_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End2_Pt1stEle_PAS"	  , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End2_hltPt1stEle_PAS"	  , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End2_Eta1stEle_PAS"	  , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End2_Phi1stEle_PAS"	  , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End2_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End2_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End2_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End2_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End2_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Electrons_End2_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);

          if ( nJetLooseEle_ptCut >= 1 ) {
            FillUserTH1D( "Electrons_End2_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_1Jet_Pt1stEle_PAS"        , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_1Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_1Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End2_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut >= 2 ) {
            FillUserTH1D( "Electrons_End2_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_2Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_2Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_2Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End2_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut >= 3 ) {
            FillUserTH1D( "Electrons_End2_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_3Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_3Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_3Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End2_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nVertex >= 0 && nVertex <= 5 ){
            FillUserTH1D( "Electrons_End2_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Electrons_End2_PU0-8_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU0-8_Dist1stEle_PAS"        , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }						          

          if ( nVertex >= 6 && nVertex <= 10 ){      
            FillUserTH1D( "Electrons_End2_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Electrons_End2_PU9-UP_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU9-UP_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU9-UP_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }
        }

        //// for event list
        //// see if we passed full HEEP
        //if(LooseEle1_TrkIsoHEEP7 < 5.0)
        ////if(LooseEle1_PassHEEPID)
        //{
        //  if(run>=284043 && run<=284044)
        //  {
        //    if ( fabs(LooseEle1_SCEta) < 1.4442 )
        //    {
        //      std::cout << (uint64_t)run << " " << (uint64_t)ls << " " << (uint64_t)event << std::endl;

        //      if ( H_Photon22   > 0.1 && LooseEle1_hltPhotonPt >= 22.  && LooseEle1_hltPhotonPt < 30. ) { std::cout << "\tTrigger: HLT_Photon22"<<std::endl; } 
        //      if ( H_Photon30   > 0.1 && LooseEle1_hltPhotonPt >= 30.  && LooseEle1_hltPhotonPt < 36. ) { std::cout << "\tTrigger: HLT_Photon30"<<std::endl; } 
        //      if ( H_Photon36   > 0.1 && LooseEle1_hltPhotonPt >= 36.  && LooseEle1_hltPhotonPt < 50. ) { std::cout << "\tTrigger: HLT_Photon36"<<std::endl; } 
        //      if ( H_Photon50   > 0.1 && LooseEle1_hltPhotonPt >= 50.  && LooseEle1_hltPhotonPt < 75. ) { std::cout << "\tTrigger: HLT_Photon50"<<std::endl; } 
        //      if ( H_Photon75   > 0.1 && LooseEle1_hltPhotonPt >= 75.  && LooseEle1_hltPhotonPt < 90. ) { std::cout << "\tTrigger: HLT_Photon75"<<std::endl; } 
        //      if ( H_Photon90   > 0.1 && LooseEle1_hltPhotonPt >= 90.  && LooseEle1_hltPhotonPt < 120.) { std::cout << "\tTrigger: HLT_Photon90"<<std::endl; } 
        //      if ( H_Photon120  > 0.1 && LooseEle1_hltPhotonPt >= 120. && LooseEle1_hltPhotonPt < 175.) { std::cout << "\tTrigger: HLT_Photon120"<<std::endl; } 
        //      if ( H_Photon175  > 0.1 && LooseEle1_hltPhotonPt >= 175.)                                 { std::cout << "\tTrigger: HLT_Photon175"<<std::endl; } 
        //      std::cout << "\tPassEt: " << pass_et << " --> " << LooseEle1_PtHeep << std::endl
        //        << "\tPassEcalDriven: " << pass_ecalDriven << " --> " << LooseEle1_EcalDriven << std::endl
        //        << "\tPassDeltaPhi: " << pass_deltaPhi << " --> " << LooseEle1_DeltaPhiTrkSC << std::endl
        //        << "\tPassMissingHits: " << pass_missingHits << " --> " << LooseEle1_MissingHits << std::endl
        //        << "\tPassDeltaEtaSeed: " << pass_deltaEtaSeed << " --> " << LooseEle1_DeltaEtaSeed << std::endl;
        //      if ( fabs(LooseEle1_SCEta) > 1.566 && fabs(LooseEle1_SCEta) < 2.5 )
        //        std::cout << "\tPassSigmaIetaIeta: " << pass_sigmaIEtaIEta << " --> " << LooseEle1_Full5x5SigmaIEtaIEta << " [ENDCAP]"<< std::endl;
        //      else
        //        std::cout << "\tPassSigmaIetaIeta: " << pass_sigmaIEtaIEta << " --> " << LooseEle1_Full5x5SigmaIEtaIEta << " [BARREL]"<< std::endl;
        //      if ( fabs(LooseEle1_SCEta) < 1.4442 )
        //      {
        //        std::cout << "\tPassShowerShape: " << pass_shape << " --> E1x5/E5x5=" << LooseEle1_Full5x5E1x5OverE5x5 << ", E2x5/E5x5=" << LooseEle1_Full5x5E2x5OverE5x5 << " [BARREL]" << std::endl;
        //        std::cout << "\tPassCaloIso: " << pass_caloIsolation << " --> " << caloIsolation << " < " << ( 2.0 + ( 0.03 * LooseEle1_PtHeep ) + (0.28 * LooseEle1_RhoForHeep ) ) << std::endl;
        //        std::cout << "\tPassDxy: " << pass_dxy << " --> " << LooseEle1_LeadVtxDistXY << " < 0.02 " << std::endl
        //          << "\tPassH/E: " << pass_hoe << " --> " << LooseEle1_HoE << " < " << 1/LooseEle1_SCEnergy + 0.05 << std::endl;
        //      }
        //      else
        //      {
        //        std::cout << "\tPassCaloIso: " << pass_caloIsolation << " --> " << caloIsolation << " < " << (LooseEle1_PtHeep<50 ? (2.5+( 0.28 * LooseEle1_RhoForHeep ) ) : (2.5+( 0.28 * LooseEle1_RhoForHeep ) + ( 0.03 * (LooseEle1_PtHeep - 50.0 ) )) ) << std::endl;
        //        std::cout << "\tPassDxy: " << pass_dxy << " --> " << LooseEle1_LeadVtxDistXY << " < 0.05" << std::endl
        //          << "\tPassH/E: " << pass_hoe << " --> " << LooseEle1_HoE << " < " << 5/LooseEle1_SCEnergy + 0.05 << std::endl;
        //      }
        //      std::cout << "\tPassTrkIso: " << (LooseEle1_TrkIsoHEEP7 < 5.0) << " --> " << LooseEle1_TrkIsoHEEP7 << " < 5.0" << std::endl;
        //    }
        //  }
        //}
        ////

      }
      else if (pass_hoe && pass_caloIsolation) {  // "jets" must at least pass HoE and CaloIsolation cuts

        if ( isBarrel ) { 
          FillUserTH1D( "Jets_Bar_nVertex_PAS"           , nVertex                      , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_Pt1stEle_PAS"	  , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_hltPt1stEle_PAS"	  , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_Eta1stEle_PAS"	  , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_Phi1stEle_PAS"	  , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Jets_Bar_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);

          if ( nJetLooseEle_ptCut >= 1 ) {
            FillUserTH1D( "Jets_Bar_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_1Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_1Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_1Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_Bar_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut >= 2 ) {
            FillUserTH1D( "Jets_Bar_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_2Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_2Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_2Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_Bar_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut >= 3 ) {
            FillUserTH1D( "Jets_Bar_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_3Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_3Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_3Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_Bar_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nVertex >= 0 && nVertex <= 5 ){
            FillUserTH1D( "Jets_Bar_PU0-8_nElectron_PAS"          , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU0-8_Pt1stEle_PAS"           , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU0-8_Eta1stEle_PAS"          , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU0-8_Phi1stEle_PAS"          , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU0-8_Charge1stEle_PAS"       , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Jets_Bar_PU0-8_DCotTheta1stEle_PAS"    , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU0-8_Dist1stEle_PAS"         , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU0-8_MET_PAS"                , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }						          

          if ( nVertex >= 6 && nVertex <= 10 ){      
            FillUserTH1D( "Jets_Bar_PU9-UP_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU9-UP_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU9-UP_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU9-UP_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU9-UP_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Jets_Bar_PU9-UP_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU9-UP_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU9-UP_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }	  
        }

        if ( isEndcap1 ) { 
          FillUserTH1D( "Jets_End1_nVertex_PAS"          , nVertex                      , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End1_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End1_hltPt1stEle_PAS"	  , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End1_Pt1stEle_PAS"	  , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End1_Eta1stEle_PAS"	  , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End1_Phi1stEle_PAS"	  , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End1_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End1_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End1_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End1_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End1_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Jets_End1_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);

          if ( nJetLooseEle_ptCut >= 1 ) {
            FillUserTH1D( "Jets_End1_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_1Jet_Pt1stEle_PAS"        , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_1Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_1Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End1_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);

          }

          if ( nJetLooseEle_ptCut >= 2 ) {
            FillUserTH1D( "Jets_End1_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_2Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_2Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_2Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End1_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut >= 3 ) {
            FillUserTH1D( "Jets_End1_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_3Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_3Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_3Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End1_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nVertex >= 0 && nVertex <= 5 ){
            FillUserTH1D( "Jets_End1_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Jets_End1_PU0-8_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU0-8_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU0-8_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }						          

          if ( nVertex >= 6 && nVertex <= 10 ){      
            FillUserTH1D( "Jets_End1_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Jets_End1_PU9-UP_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU9-UP_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU9-UP_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }
        }

        if ( isEndcap2 ) { 
          FillUserTH1D( "Jets_End2_nVertex_PAS"          , nVertex                      , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End2_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End2_Pt1stEle_PAS"	  , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End2_hltPt1stEle_PAS"	  , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End2_Eta1stEle_PAS"	  , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End2_Phi1stEle_PAS"	  , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End2_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End2_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End2_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End2_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End2_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Jets_End2_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);

          if ( nJetLooseEle_ptCut >= 1 ) {
            FillUserTH1D( "Jets_End2_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_1Jet_Pt1stEle_PAS"        , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_1Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_1Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End2_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut >= 2 ) {
            FillUserTH1D( "Jets_End2_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_2Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_2Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_2Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End2_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut >= 3 ) {
            FillUserTH1D( "Jets_End2_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_3Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_3Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_3Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End2_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }

          if ( nVertex >= 0 && nVertex <= 5 ){
            FillUserTH1D( "Jets_End2_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Jets_End2_PU0-8_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU0-8_Dist1stEle_PAS"        , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }						          

          if ( nVertex >= 6 && nVertex <= 10 ){      
            FillUserTH1D( "Jets_End2_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Jets_End2_PU9-UP_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU9-UP_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU9-UP_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }
        }
      }
    }
  } // End loop over events

  std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
