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

float weightTruePileupV10toHcp53X(float input){
  float w[60] = {
    0.409409,
    0.527276,
    0.39328,
    0.507892,
    0.48029,
    0.787701,
    0.632356,
    0.618033,
    0.806089,
    1.14018,
    1.5788,
    1.93507,
    1.957,
    1.73004,
    1.46737,
    1.28278,
    1.18189,
    1.13388,
    1.12578,
    1.14415,
    1.16048,
    1.1618,
    1.15318,
    1.13405,
    1.09239,
    1.01915,
    0.914837,
    0.786744,
    0.644879,
    0.502039,
    0.371688,
    0.263586,
    0.18067,
    0.120472,
    0.0780184,
    0.0486113,
    0.0289039,
    0.0163367,
    0.00879674,
    0.00456046,
    0.0023098,
    0.00115977,
    0.000583207,
    0.000294815,
    0.000149865,
    7.62892e-05,
    3.87537e-05,
    1.96105e-05,
    9.87744e-06,
    4.95418e-06,
    2.47913e-06,
    1.23919e-06,
    6.19751e-07,
    3.10125e-07,
    1.54934e-07,
    7.71425e-08,
    3.8182e-08,
    1.87455e-08,
    9.10765e-09,
    9.19802e-09};
  return w[(int)floor(input)];
}


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

   CreateUserTH1D("Total_minPrescale"                     ,   200 , 0.0      , 2000.0   ); CreateUserTH1D("Pass_minPrescale"                       ,   200 , 0.0      , 2000.0   ); 
   CreateUserTH1D("Total_pileupWeight"                    ,   100 , 0.0      , 2.0      ); CreateUserTH1D("Pass_pileupWeight"                      ,   100 , 0.0      , 2.0      ); 
   
   CreateUserTH1D("Total_Bar_nVertex_PAS"                 ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Pass_Bar_nVertex_PAS"                   ,   100 , 0.0      , 100.0    );  
   CreateUserTH1D("Total_End1_nVertex_PAS"                ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Pass_End1_nVertex_PAS"                  ,   100 , 0.0      , 100.0    );  
   CreateUserTH1D("Total_End2_nVertex_PAS"                ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Pass_End2_nVertex_PAS"                  ,   100 , 0.0      , 100.0    );  
   
   // inclusive

   CreateUserTH1D( "Total_Bar_nElectron_PAS"              ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_Bar_nElectron_PAS"               ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_Bar_Pt1stEle_PAS"               ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_Bar_Pt1stEle_PAS"	           ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_Bar_hltPt1stEle_PAS"            ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_Bar_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_Bar_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_Bar_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_Bar_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_Bar_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_Bar_Charge1stEle_PAS"	          ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_Bar_Charge1stEle_PAS"            ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_Bar_DCotTheta1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_DCotTheta1stEle_PAS"         ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_Dist1stEle_PAS"             ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_Dist1stEle_PAS"              ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_MET_PAS"                    ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_Bar_MET_PAS"                     ,    100 , 0.0     , 1000.    );
						          										           
   CreateUserTH1D( "Total_End1_nElectron_PAS"             ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End1_nElectron_PAS"              ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End1_Pt1stEle_PAS"              ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_End1_Pt1stEle_PAS"	           ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_End1_hltPt1stEle_PAS"           ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_End1_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_End1_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End1_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End1_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End1_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End1_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End1_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End1_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_DCotTheta1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_Dist1stEle_PAS"             ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_MET_PAS"                   ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_End1_MET_PAS"                    ,    100 , 0.0     , 1000.    );

   CreateUserTH1D( "Total_End2_nElectron_PAS"             ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End2_nElectron_PAS"              ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End2_Pt1stEle_PAS"              ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_End2_Pt1stEle_PAS"	           ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_End2_hltPt1stEle_PAS"           ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_End2_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_End2_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End2_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End2_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End2_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End2_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End2_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End2_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_DCotTheta1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_Dist1stEle_PAS"             ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_MET_PAS"                   ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_End2_MET_PAS"                    ,    100 , 0.0     , 1000.    );

   // 1 jet 

   CreateUserTH1D( "Total_Bar_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_Bar_1Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_Bar_1Jet_Pt1stEle_PAS"          ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_Bar_1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_Bar_1Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_Bar_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_Bar_1Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_Bar_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_Bar_1Jet_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_Bar_1Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_Bar_1Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_1Jet_DCotTheta1stEle_PAS"    ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_1Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_1Jet_Dist1stEle_PAS"         ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_Bar_1Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );
						          										           
   CreateUserTH1D( "Total_End1_1Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End1_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End1_1Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_End1_1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_End1_1Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End1_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End1_1Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End1_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End1_1Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End1_1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End1_1Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_1Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_1Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_1Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_1Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_End1_1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );

   CreateUserTH1D( "Total_End2_1Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End2_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End2_1Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_End2_1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_End2_1Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End2_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End2_1Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End2_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End2_1Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End2_1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End2_1Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_1Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_1Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_1Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_1Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_End2_1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );

   // 2 jet 

   CreateUserTH1D( "Total_Bar_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_Bar_2Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_Bar_2Jet_Pt1stEle_PAS"          ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_Bar_2Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_Bar_2Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_Bar_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_Bar_2Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_Bar_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_Bar_2Jet_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_Bar_2Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_Bar_2Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_2Jet_DCotTheta1stEle_PAS"    ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_2Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_2Jet_Dist1stEle_PAS"         ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_2Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_Bar_2Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );
						          										           
   CreateUserTH1D( "Total_End1_2Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End1_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End1_2Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_End1_2Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_End1_2Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End1_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End1_2Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End1_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End1_2Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End1_2Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End1_2Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_2Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_2Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_2Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_2Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_End1_2Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );

   CreateUserTH1D( "Total_End2_2Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End2_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End2_2Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_End2_2Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_End2_2Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End2_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End2_2Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End2_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End2_2Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End2_2Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End2_2Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_2Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_2Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_2Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_2Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_End2_2Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );

   // 3 jet 

   CreateUserTH1D( "Total_Bar_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_Bar_3Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_Bar_3Jet_Pt1stEle_PAS"          ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_Bar_3Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_Bar_3Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_Bar_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_Bar_3Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_Bar_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_Bar_3Jet_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_Bar_3Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_Bar_3Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_3Jet_DCotTheta1stEle_PAS"    ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_3Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_3Jet_Dist1stEle_PAS"         ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_3Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_Bar_3Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );
						          										           
   CreateUserTH1D( "Total_End1_3Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End1_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End1_3Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_End1_3Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_End1_3Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End1_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End1_3Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End1_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End1_3Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End1_3Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End1_3Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_3Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_3Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_3Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_3Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_End1_3Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );

   CreateUserTH1D( "Total_End2_3Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End2_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End2_3Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_End2_3Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_End2_3Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End2_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End2_3Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End2_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End2_3Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End2_3Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End2_3Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_3Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_3Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_3Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_3Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_End2_3Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );


   // Runs 160406 - 166502
   
   CreateUserTH1D( "Total1_Bar_nElectron_PAS"             ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass1_Bar_nElectron_PAS"              ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total1_Bar_Pt1stEle_PAS"              ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass1_Bar_Pt1stEle_PAS"	           ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total1_Bar_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass1_Bar_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total1_Bar_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass1_Bar_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total1_Bar_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass1_Bar_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total1_Bar_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass1_Bar_DCotTheta1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total1_Bar_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass1_Bar_Dist1stEle_PAS"             ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total1_Bar_MET_PAS"                   ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass1_Bar_MET_PAS"                    ,    100 , 0.0     , 1000.    );
						          										           
   CreateUserTH1D( "Total1_End1_nElectron_PAS"            ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass1_End1_nElectron_PAS"             ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total1_End1_Pt1stEle_PAS"             ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass1_End1_Pt1stEle_PAS"	           ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total1_End1_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass1_End1_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total1_End1_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass1_End1_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total1_End1_Charge1stEle_PAS"         ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass1_End1_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total1_End1_DCotTheta1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass1_End1_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total1_End1_Dist1stEle_PAS"           ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass1_End1_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total1_End1_MET_PAS"                  ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass1_End1_MET_PAS"                   ,    100 , 0.0     , 1000.    );

   CreateUserTH1D( "Total1_End2_nElectron_PAS"            ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass1_End2_nElectron_PAS"             ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total1_End2_Pt1stEle_PAS"             ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass1_End2_Pt1stEle_PAS"	           ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total1_End2_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass1_End2_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total1_End2_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass1_End2_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total1_End2_Charge1stEle_PAS"         ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass1_End2_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total1_End2_DCotTheta1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass1_End2_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total1_End2_Dist1stEle_PAS"           ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass1_End2_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total1_End2_MET_PAS"                  ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass1_End2_MET_PAS"                   ,    100 , 0.0     , 1000.    );

   // Runs 166503 - onward
   
   CreateUserTH1D( "Total2_Bar_nElectron_PAS"             ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass2_Bar_nElectron_PAS"              ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total2_Bar_Pt1stEle_PAS"              ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass2_Bar_Pt1stEle_PAS"	           ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total2_Bar_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass2_Bar_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total2_Bar_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass2_Bar_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total2_Bar_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass2_Bar_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total2_Bar_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass2_Bar_DCotTheta1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total2_Bar_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass2_Bar_Dist1stEle_PAS"             ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total2_Bar_MET_PAS"                   ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass2_Bar_MET_PAS"                    ,    100 , 0.0     , 1000.    );
						          										           
   CreateUserTH1D( "Total2_End1_nElectron_PAS"            ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass2_End1_nElectron_PAS"             ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total2_End1_Pt1stEle_PAS"             ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass2_End1_Pt1stEle_PAS"	           ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total2_End1_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass2_End1_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total2_End1_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass2_End1_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total2_End1_Charge1stEle_PAS"         ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass2_End1_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total2_End1_DCotTheta1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass2_End1_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total2_End1_Dist1stEle_PAS"           ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass2_End1_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total2_End1_MET_PAS"                  ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass2_End1_MET_PAS"                   ,    100 , 0.0     , 1000.    );
   
   CreateUserTH1D( "Total2_End2_nElectron_PAS"            ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass2_End2_nElectron_PAS"             ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total2_End2_Pt1stEle_PAS"             ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass2_End2_Pt1stEle_PAS"	           ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total2_End2_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass2_End2_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total2_End2_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass2_End2_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total2_End2_Charge1stEle_PAS"         ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass2_End2_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total2_End2_DCotTheta1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass2_End2_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total2_End2_Dist1stEle_PAS"           ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass2_End2_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total2_End2_MET_PAS"                  ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass2_End2_MET_PAS"                   ,    100 , 0.0     , 1000.    );

   // Pile-up [0,5]

   CreateUserTH1D( "Total_Bar_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_Bar_PU0-8_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_Bar_PU0-8_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_Bar_PU0-8_Pt1stEle_PAS"	   ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_Bar_PU0-8_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_Bar_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	 );
   CreateUserTH1D( "Total_Bar_PU0-8_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_Bar_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_Bar_PU0-8_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_Bar_PU0-8_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_Bar_PU0-8_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_PU0-8_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_PU0-8_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_PU0-8_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_PU0-8_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_Bar_PU0-8_MET_PAS"               ,    100 , 0.0     , 1000.    );

   CreateUserTH1D( "Total_End1_PU0-8_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End1_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End1_PU0-8_Pt1stEle_PAS"        ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_End1_PU0-8_Pt1stEle_PAS"	   ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_End1_PU0-8_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End1_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	 );
   CreateUserTH1D( "Total_End1_PU0-8_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End1_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End1_PU0-8_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End1_PU0-8_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End1_PU0-8_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_PU0-8_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_PU0-8_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_PU0-8_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_PU0-8_MET_PAS"             ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_End1_PU0-8_MET_PAS"              ,    100 , 0.0     , 1000.    );

   CreateUserTH1D( "Total_End2_PU0-8_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End2_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End2_PU0-8_Pt1stEle_PAS"        ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_End2_PU0-8_Pt1stEle_PAS"	   ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_End2_PU0-8_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End2_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	 );
   CreateUserTH1D( "Total_End2_PU0-8_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End2_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End2_PU0-8_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End2_PU0-8_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End2_PU0-8_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_PU0-8_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_PU0-8_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_PU0-8_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_PU0-8_MET_PAS"             ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_End2_PU0-8_MET_PAS"              ,    100 , 0.0     , 1000.    );

   // Pile-up [6-10]

   CreateUserTH1D( "Total_Bar_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_Bar_PU9-UP_nElectron_PAS"        ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_Bar_PU9-UP_Pt1stEle_PAS"        ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_Bar_PU9-UP_Pt1stEle_PAS"	   ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_Bar_PU9-UP_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_Bar_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	 );
   CreateUserTH1D( "Total_Bar_PU9-UP_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_Bar_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_Bar_PU9-UP_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_Bar_PU9-UP_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_Bar_PU9-UP_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_PU9-UP_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_PU9-UP_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_PU9-UP_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_Bar_PU9-UP_MET_PAS"              ,    100 , 0.0     , 1000.    );

   CreateUserTH1D( "Total_End1_PU9-UP_nElectron_PAS"      ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End1_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End1_PU9-UP_Pt1stEle_PAS"       ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_End1_PU9-UP_Pt1stEle_PAS"	   ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_End1_PU9-UP_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End1_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	 );
   CreateUserTH1D( "Total_End1_PU9-UP_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End1_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End1_PU9-UP_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End1_PU9-UP_Charge1stEle_PAS"    ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End1_PU9-UP_DCotTheta1stEle_PAS",    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_PU9-UP_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_PU9-UP_Dist1stEle_PAS"     ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_PU9-UP_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_PU9-UP_MET_PAS"            ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_End1_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1000.    );

   CreateUserTH1D( "Total_End2_PU9-UP_nElectron_PAS"      ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End2_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End2_PU9-UP_Pt1stEle_PAS"       ,    1000, 0       , 1000     );  CreateUserTH1D( "Pass_End2_PU9-UP_Pt1stEle_PAS"	   ,    1000, 0       , 1000     );
   CreateUserTH1D( "Total_End2_PU9-UP_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End2_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	 );
   CreateUserTH1D( "Total_End2_PU9-UP_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End2_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End2_PU9-UP_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End2_PU9-UP_Charge1stEle_PAS"    ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End2_PU9-UP_DCotTheta1stEle_PAS",    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_PU9-UP_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_PU9-UP_Dist1stEle_PAS"     ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_PU9-UP_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_PU9-UP_MET_PAS"            ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Pass_End2_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1000.    );

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
     
     double pileup_weight;
     if ( isData == 0 ) pileup_weight = weightTruePileupV10toHcp53X ( nPileUpInt_True );
     else pileup_weight = 1.0;
     
     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------
     
     short min_prescale = 0;
     int passTrigger = 0;

     if ( isData ) {
       
       if ( LooseEle1_hltPhotonPt > 0.0 ) { 
       	 if ( H_Photon30_CIdVL > 0.1 && LooseEle1_hltPhotonPt >= 30.  && LooseEle1_hltPhotonPt < 50. ) { passTrigger = 1; min_prescale = H_Photon30_CIdVL; } 
       	 if ( H_Photon50_CIdVL > 0.1 && LooseEle1_hltPhotonPt >= 50.  && LooseEle1_hltPhotonPt < 75. ) { passTrigger = 1; min_prescale = H_Photon50_CIdVL; } 
       	 if ( H_Photon75_CIdVL > 0.1 && LooseEle1_hltPhotonPt >= 75.  && LooseEle1_hltPhotonPt < 90. ) { passTrigger = 1; min_prescale = H_Photon75_CIdVL; } 
       	 if ( H_Photon90_CIdVL > 0.1 && LooseEle1_hltPhotonPt >= 90.  && LooseEle1_hltPhotonPt < 135.) { passTrigger = 1; min_prescale = H_Photon90_CIdVL; } 
       	 if ( H_Photon135      > 0.1 && LooseEle1_hltPhotonPt >= 135. && LooseEle1_hltPhotonPt < 150.) { passTrigger = 1; min_prescale = H_Photon135     ; } 
       	 if ( H_Photon150      > 0.1 && LooseEle1_hltPhotonPt >= 150.                                ) { passTrigger = 1; min_prescale = H_Photon150     ; } 
       }
     }  // end if (isData) 
					       
     else { 
       min_prescale = 1;
       passTrigger = 1 ;
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
     										            
     // Filters									            
     fillVariableWithValue(   "PassHBHENoiseFilter"	      , PassHBHENoiseFilter                              , min_prescale * pileup_weight );
     fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight                          , min_prescale * pileup_weight );
     fillVariableWithValue(   "PassBadEESupercrystalFilter"   , ( isData == 1 ) ? PassBadEESupercrystalFilter : 1, min_prescale * pileup_weight );
     fillVariableWithValue(   "PassBeamScraping"	      , ( isData == 1 ) ? PassBeamScraping	      : 1, min_prescale * pileup_weight );
     fillVariableWithValue(   "PassEcalDeadCellBoundEnergy"   , PassEcalDeadCellBoundEnergy                      , min_prescale * pileup_weight );
     fillVariableWithValue(   "PassEcalDeadCellTrigPrim"      , PassEcalDeadCellTrigPrim                         , min_prescale * pileup_weight );
     fillVariableWithValue(   "PassEcalLaserCorrFilter"       , ( isData == 1 ) ? PassEcalLaserCorrFilter     : 1, min_prescale * pileup_weight );
     fillVariableWithValue(   "PassHcalLaserEventFilter"      , ( isData == 1 ) ? PassHcalLaserEventFilter    : 1, min_prescale * pileup_weight );
     fillVariableWithValue(   "PassPhysDecl"		      , ( isData == 1 ) ? PassPhysDecl		      : 1, min_prescale * pileup_weight );
     fillVariableWithValue(   "PassPrimaryVertex"	      , PassPrimaryVertex                                , min_prescale * pileup_weight );
     fillVariableWithValue(   "PassTrackingFailure"	      , ( isData == 1 ) ? PassTrackingFailure	      : 1, min_prescale * pileup_weight );
     
     // Muon variables ( for veto ) 					      	       
     fillVariableWithValue(   "nMuon"                         , nMuon_ptCut                 , min_prescale * pileup_weight );
			                                      		                    
     // 1st Electron variables				      		                    
     fillVariableWithValue(   "nEle"                          , nLooseEle_ptCut             , min_prescale * pileup_weight );
     
     fillVariableWithValue(   "Pt1stEle"                      , LooseEle1_Pt                , min_prescale * pileup_weight );
     fillVariableWithValue(   "Eta1stEle"                     , LooseEle1_Eta               , min_prescale * pileup_weight );
     fillVariableWithValue(   "Phi1stEle"                     , LooseEle1_Phi               , min_prescale * pileup_weight );
     fillVariableWithValue(   "HLTPt1stEle"                   , LooseEle1_hltPhotonPt       , min_prescale * pileup_weight );
     fillVariableWithValue(   "H_Photon30_CIdVL"              , H_Photon30_CIdVL            , min_prescale * pileup_weight );
     fillVariableWithValue(   "H_Photon50_CIdVL"              , H_Photon50_CIdVL            , min_prescale * pileup_weight );   
     fillVariableWithValue(   "H_Photon75_CIdVL"              , H_Photon75_CIdVL            , min_prescale * pileup_weight );   
     fillVariableWithValue(   "H_Photon90_CIdVL"              , H_Photon90_CIdVL            , min_prescale * pileup_weight );   
     fillVariableWithValue(   "H_Photon135"                   , H_Photon135                 , min_prescale * pileup_weight );   
     fillVariableWithValue(   "H_Photon150"                   , H_Photon150                 , min_prescale * pileup_weight );    
     
     // 1st JET variables                                     		                    
     fillVariableWithValue(   "nJet"                          , nJetLooseEle_ptCut          , min_prescale * pileup_weight );
										            
     // MET									            
     fillVariableWithValue(   "MET"                           , PFMET_Type01XY_Pt           , min_prescale * pileup_weight );
										            
     // Dummy variables								            
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
     
     evaluateCuts();

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
	 FillUserTH1D( "Total_Bar_MET_PAS"               , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);

	 if ( nJetLooseEle_ptCut >= 1 ) {
	   FillUserTH1D( "Total_Bar_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_1Jet_Pt1stEle_PAS"	      , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_1Jet_Eta1stEle_PAS"	      , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_1Jet_Phi1stEle_PAS"	      , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total_Bar_1Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_1Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_1Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 }

	 if ( nJetLooseEle_ptCut >= 2 ) {
	   FillUserTH1D( "Total_Bar_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_2Jet_Pt1stEle_PAS"	      , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_2Jet_Eta1stEle_PAS"	      , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_2Jet_Phi1stEle_PAS"	      , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total_Bar_2Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_2Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_2Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 }

	 if ( nJetLooseEle_ptCut >= 3 ) {
	   FillUserTH1D( "Total_Bar_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_3Jet_Pt1stEle_PAS"	      , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_3Jet_Eta1stEle_PAS"	      , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_3Jet_Phi1stEle_PAS"	      , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total_Bar_3Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_3Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_3Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 }

	 if ( isData && run <= 166502 ) { 
	   FillUserTH1D( "Total1_Bar_nElectron_PAS"               , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_Bar_Pt1stEle_PAS"                , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_Bar_Eta1stEle_PAS"	          , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_Bar_Phi1stEle_PAS"	          , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_Bar_Charge1stEle_PAS"            , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total1_Bar_DCotTheta1stEle_PAS"         , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_Bar_Dist1stEle_PAS"              , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_Bar_MET_PAS"                     , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 } 

	 if ( isData && run >= 166503 ) { 
	   FillUserTH1D( "Total2_Bar_nElectron_PAS"               , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_Bar_Pt1stEle_PAS"                , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_Bar_Eta1stEle_PAS"	          , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_Bar_Phi1stEle_PAS"	          , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_Bar_Charge1stEle_PAS"            , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total2_Bar_DCotTheta1stEle_PAS"         , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_Bar_Dist1stEle_PAS"              , LooseEle1_Dist               , min_prescale * pileup_weight);
 	   FillUserTH1D( "Total2_Bar_MET_PAS"                     , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 } 

	 if ( nVertex >= 0 && nVertex <= 5 ){
	   FillUserTH1D( "Total_Bar_PU0-8_nElectron_PAS"          , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_PU0-8_Pt1stEle_PAS"           , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_PU0-8_Eta1stEle_PAS"          , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_PU0-8_Phi1stEle_PAS"          , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_PU0-8_Charge1stEle_PAS"       , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total_Bar_PU0-8_DCotTheta1stEle_PAS"    , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_PU0-8_Dist1stEle_PAS"         , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_PU0-8_MET_PAS"                , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 }						          
							          
	 if ( nVertex >= 6 && nVertex <= 10 ){      
	   FillUserTH1D( "Total_Bar_PU9-UP_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_PU9-UP_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_PU9-UP_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_PU9-UP_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_PU9-UP_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total_Bar_PU9-UP_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_PU9-UP_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_Bar_PU9-UP_MET_PAS"               , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
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
	 FillUserTH1D( "Total_End1_MET_PAS"               , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 

	 if ( nJetLooseEle_ptCut >= 1 ) {
	   FillUserTH1D( "Total_End1_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_1Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total_End1_1Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_1Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_1Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 }

	 if ( nJetLooseEle_ptCut >= 2 ) {
	   FillUserTH1D( "Total_End1_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_2Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total_End1_2Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_2Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_2Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 }

	 if ( nJetLooseEle_ptCut >= 3 ) {
	   FillUserTH1D( "Total_End1_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_3Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total_End1_3Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_3Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_3Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 }


	 if ( isData && run <= 166502 ) { 
	   FillUserTH1D( "Total1_End1_nElectron_PAS"              , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_End1_Pt1stEle_PAS"               , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_End1_Eta1stEle_PAS"	          , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_End1_Phi1stEle_PAS"	          , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_End1_Charge1stEle_PAS"           , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total1_End1_DCotTheta1stEle_PAS"        , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_End1_Dist1stEle_PAS"             , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_End1_MET_PAS"                    , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 } 

	 if ( isData && run >= 166503 ) { 
	   FillUserTH1D( "Total2_End1_nElectron_PAS"              , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_End1_Pt1stEle_PAS"               , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_End1_Eta1stEle_PAS"	          , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_End1_Phi1stEle_PAS"	          , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_End1_Charge1stEle_PAS"           , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total2_End1_DCotTheta1stEle_PAS"        , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_End1_Dist1stEle_PAS"             , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_End1_MET_PAS"                    , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 } 

	 if ( nVertex >= 0 && nVertex <= 5 ){
	   FillUserTH1D( "Total_End1_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total_End1_PU0-8_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_PU0-8_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_PU0-8_MET_PAS"               , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 }						          
							          
	 if ( nVertex >= 6 && nVertex <= 10 ){      
	   FillUserTH1D( "Total_End1_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total_End1_PU9-UP_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_PU9-UP_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End1_PU9-UP_MET_PAS"              , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
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
	 FillUserTH1D( "Total_End2_MET_PAS"               , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 
	 if ( nJetLooseEle_ptCut >= 1 ) {
	   FillUserTH1D( "Total_End2_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_1Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total_End2_1Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_1Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_1Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 }

	 if ( nJetLooseEle_ptCut >= 2 ) {
	   FillUserTH1D( "Total_End2_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_2Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total_End2_2Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_2Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_2Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 }

	 if ( nJetLooseEle_ptCut >= 3 ) {
	   FillUserTH1D( "Total_End2_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_3Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total_End2_3Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_3Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_3Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 }
	 
	 if ( isData && run <= 166502 ) { 
	   FillUserTH1D( "Total1_End2_nElectron_PAS"              , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_End2_Pt1stEle_PAS"               , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_End2_Eta1stEle_PAS"	          , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_End2_Phi1stEle_PAS"	          , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_End2_Charge1stEle_PAS"           , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total1_End2_DCotTheta1stEle_PAS"        , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_End2_Dist1stEle_PAS"             , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total1_End2_MET_PAS"                    , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 } 

	 if ( isData && run >= 166503 ) { 
	   FillUserTH1D( "Total2_End2_nElectron_PAS"              , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_End2_Pt1stEle_PAS"               , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_End2_Eta1stEle_PAS"	          , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_End2_Phi1stEle_PAS"	          , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_End2_Charge1stEle_PAS"           , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total2_End2_DCotTheta1stEle_PAS"        , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_End2_Dist1stEle_PAS"             , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total2_End2_MET_PAS"                    , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 } 

	 if ( nVertex >= 0 && nVertex <= 5 ){
	   FillUserTH1D( "Total_End2_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total_End2_PU0-8_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_PU0-8_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_PU0-8_MET_PAS"               , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 }						          
							          
	 if ( nVertex >= 6 && nVertex <= 10 ){      
	   FillUserTH1D( "Total_End2_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
	   FillUserTH1D( "Total_End2_PU9-UP_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_PU9-UP_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Total_End2_PU9-UP_MET_PAS"              , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	 }
       }
       
       if ( LooseEle1_PassID == 1 ) { 
	 
	 if ( isBarrel ) { 
	   FillUserTH1D( "Pass_Bar_nVertex_PAS"           , nVertex                      , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_Bar_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_Bar_Pt1stEle_PAS"	  , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_Bar_hltPt1stEle_PAS"	  , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_Bar_Eta1stEle_PAS"	  , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_Bar_Phi1stEle_PAS"	  , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_Bar_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_Bar_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_Bar_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_Bar_MET_PAS"               , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);

	   if ( nJetLooseEle_ptCut >= 1 ) {
	     FillUserTH1D( "Pass_Bar_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_1Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_1Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_1Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_1Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   }

	   if ( nJetLooseEle_ptCut >= 2 ) {
	     FillUserTH1D( "Pass_Bar_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_2Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_2Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_2Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_2Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   }

	   if ( nJetLooseEle_ptCut >= 3 ) {
	     FillUserTH1D( "Pass_Bar_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_3Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_3Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_3Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_3Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   }
	   
	   if ( isData && run <= 166502 ) { 
	     FillUserTH1D( "Pass1_Bar_nElectron_PAS"               , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_Bar_Pt1stEle_PAS"                , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_Bar_Eta1stEle_PAS"	           , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_Bar_Phi1stEle_PAS"	           , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_Bar_Charge1stEle_PAS"            , LooseEle1_Charge             , min_prescale * pileup_weight);  
	     FillUserTH1D( "Pass1_Bar_DCotTheta1stEle_PAS"         , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_Bar_Dist1stEle_PAS"              , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_Bar_MET_PAS"                     , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   } 
	   
	   if ( isData && run >= 166503 ) { 
	     FillUserTH1D( "Pass2_Bar_nElectron_PAS"               , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_Bar_Pt1stEle_PAS"                , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_Bar_Eta1stEle_PAS"	           , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_Bar_Phi1stEle_PAS"	           , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_Bar_Charge1stEle_PAS"            , LooseEle1_Charge             , min_prescale * pileup_weight);  
	     FillUserTH1D( "Pass2_Bar_DCotTheta1stEle_PAS"         , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_Bar_Dist1stEle_PAS"              , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_Bar_MET_PAS"                     , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   } 

	   if ( nVertex >= 0 && nVertex <= 5 ){
	     FillUserTH1D( "Pass_Bar_PU0-8_nElectron_PAS"          , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_PU0-8_Pt1stEle_PAS"           , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_PU0-8_Eta1stEle_PAS"          , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_PU0-8_Phi1stEle_PAS"          , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_PU0-8_Charge1stEle_PAS"       , LooseEle1_Charge             , min_prescale * pileup_weight);  
	     FillUserTH1D( "Pass_Bar_PU0-8_DCotTheta1stEle_PAS"    , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_PU0-8_Dist1stEle_PAS"         , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_PU0-8_MET_PAS"                , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   }						          
	   
	   if ( nVertex >= 6 && nVertex <= 10 ){      
	     FillUserTH1D( "Pass_Bar_PU9-UP_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_PU9-UP_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_PU9-UP_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_PU9-UP_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_PU9-UP_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
	     FillUserTH1D( "Pass_Bar_PU9-UP_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_PU9-UP_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_Bar_PU9-UP_MET_PAS"               , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   }	  
	 }
	 
	 if ( isEndcap1 ) { 
	   FillUserTH1D( "Pass_End1_nVertex_PAS"          , nVertex                      , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End1_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End1_hltPt1stEle_PAS"	  , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End1_Pt1stEle_PAS"	  , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End1_Eta1stEle_PAS"	  , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End1_Phi1stEle_PAS"	  , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End1_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End1_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End1_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End1_MET_PAS"              , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   
	   if ( nJetLooseEle_ptCut >= 1 ) {
	     FillUserTH1D( "Pass_End1_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_1Jet_Pt1stEle_PAS"        , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_1Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_1Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_1Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	     
	   }

	   if ( nJetLooseEle_ptCut >= 2 ) {
	     FillUserTH1D( "Pass_End1_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_2Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_2Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_2Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_2Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   }

	   if ( nJetLooseEle_ptCut >= 3 ) {
	     FillUserTH1D( "Pass_End1_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_3Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_3Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_3Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_3Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   }

	   if ( isData && run <= 166502 ) { 
	     FillUserTH1D( "Pass1_End1_nElectron_PAS"              , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_End1_Pt1stEle_PAS"               , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_End1_Eta1stEle_PAS"	           , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_End1_Phi1stEle_PAS"	           , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_End1_Charge1stEle_PAS"           , LooseEle1_Charge             , min_prescale * pileup_weight);  
	     FillUserTH1D( "Pass1_End1_DCotTheta1stEle_PAS"        , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_End1_Dist1stEle_PAS"             , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_End1_MET_PAS"                    , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   } 
	   
	   if ( isData && run >= 166503 ) { 
	     FillUserTH1D( "Pass2_End1_nElectron_PAS"               , nLooseEle_ptCut             , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_End1_Pt1stEle_PAS"                , LooseEle1_Pt                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_End1_Eta1stEle_PAS"	           , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_End1_Phi1stEle_PAS"	           , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_End1_Charge1stEle_PAS"           , LooseEle1_Charge             , min_prescale * pileup_weight);  
	     FillUserTH1D( "Pass2_End1_DCotTheta1stEle_PAS"        , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_End1_Dist1stEle_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   } 

	   if ( nVertex >= 0 && nVertex <= 5 ){
	     FillUserTH1D( "Pass_End1_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
	     FillUserTH1D( "Pass_End1_PU0-8_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_PU0-8_Dist1stEle_PAS"        , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_PU0-8_MET_PAS"               , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   }						          
	   
	   if ( nVertex >= 6 && nVertex <= 10 ){      
	     FillUserTH1D( "Pass_End1_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
	     FillUserTH1D( "Pass_End1_PU9-UP_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_PU9-UP_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End1_PU9-UP_MET_PAS"              , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   }
	 }

	 if ( isEndcap2 ) { 
	   FillUserTH1D( "Pass_End2_nVertex_PAS"          , nVertex                      , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End2_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End2_Pt1stEle_PAS"	  , LooseEle1_Pt                 , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End2_hltPt1stEle_PAS"	  , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End2_Eta1stEle_PAS"	  , LooseEle1_Eta                , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End2_Phi1stEle_PAS"	  , LooseEle1_Phi                , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End2_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End2_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End2_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
	   FillUserTH1D( "Pass_End2_MET_PAS"              , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   
	   if ( nJetLooseEle_ptCut >= 1 ) {
	     FillUserTH1D( "Pass_End2_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_1Jet_Pt1stEle_PAS"        , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_1Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_1Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_1Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   }

	   if ( nJetLooseEle_ptCut >= 2 ) {
	     FillUserTH1D( "Pass_End2_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_2Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_2Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_2Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_2Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   }

	   if ( nJetLooseEle_ptCut >= 3 ) {
	     FillUserTH1D( "Pass_End2_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_3Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_3Jet_DCotTheta1stEle_PAS" , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_3Jet_Dist1stEle_PAS"      , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_3Jet_MET_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   }
	   
	   if ( isData && run <= 166502 ) { 
	     FillUserTH1D( "Pass1_End2_nElectron_PAS"              , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_End2_Pt1stEle_PAS"               , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_End2_Eta1stEle_PAS"	           , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_End2_Phi1stEle_PAS"	           , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_End2_Charge1stEle_PAS"           , LooseEle1_Charge             , min_prescale * pileup_weight);  
	     FillUserTH1D( "Pass1_End2_DCotTheta1stEle_PAS"        , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_End2_Dist1stEle_PAS"             , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass1_End2_MET_PAS"                    , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   } 
	   
	   if ( isData && run >= 166503 ) { 
	     FillUserTH1D( "Pass2_End2_nElectron_PAS"              , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_End2_Pt1stEle_PAS"               , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_End2_Eta1stEle_PAS"	           , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_End2_Phi1stEle_PAS"	           , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_End2_Charge1stEle_PAS"           , LooseEle1_Charge             , min_prescale * pileup_weight);  
	     FillUserTH1D( "Pass2_End2_DCotTheta1stEle_PAS"        , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass2_End2_Dist1stEle_PAS"             , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   } 

	   if ( nVertex >= 0 && nVertex <= 5 ){
	     FillUserTH1D( "Pass_End2_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
	     FillUserTH1D( "Pass_End2_PU0-8_DCotTheta1stEle_PAS"   , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_PU0-8_Dist1stEle_PAS"        , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   }						          
	   
	   if ( nVertex >= 6 && nVertex <= 10 ){      
	     FillUserTH1D( "Pass_End2_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
	     FillUserTH1D( "Pass_End2_PU9-UP_DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta          , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_PU9-UP_Dist1stEle_PAS"       , LooseEle1_Dist               , min_prescale * pileup_weight);
	     FillUserTH1D( "Pass_End2_PU9-UP_MET_PAS"              , PFMET_Type01XY_Pt            , min_prescale * pileup_weight);
	   }
	 }
       }
     }
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
