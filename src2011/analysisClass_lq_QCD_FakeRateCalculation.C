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
    
   double trigger_tolerance = getPreCutValue1("trigger_tolerance"); 
   
   //--------------------------------------------------------------------------
   // Create TH1D's
   //--------------------------------------------------------------------------

   // inclusive

   CreateUserTH1D( "Total_Bar_nElectron_PAS"              ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_Bar_nElectron_PAS"               ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_Bar_Pt1stEle_PAS"               ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_Bar_Pt1stEle_PAS"	           ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_Bar_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_Bar_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_Bar_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_Bar_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_Bar_Charge1stEle_PAS"	          ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_Bar_Charge1stEle_PAS"            ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_Bar_DCotTheta1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_DCotTheta1stEle_PAS"         ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_Dist1stEle_PAS"             ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_Dist1stEle_PAS"              ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_MET_PAS"                    ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_MET_PAS"                     ,    100 , 0.0     , 1.0      );
						          										           
   CreateUserTH1D( "Total_End1_nElectron_PAS"             ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End1_nElectron_PAS"              ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End1_Pt1stEle_PAS"              ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_End1_Pt1stEle_PAS"	           ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_End1_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End1_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End1_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End1_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End1_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End1_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End1_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_DCotTheta1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_Dist1stEle_PAS"             ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_MET_PAS"                   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_MET_PAS"                    ,    100 , 0.0     , 1.0      );

   CreateUserTH1D( "Total_End2_nElectron_PAS"             ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End2_nElectron_PAS"              ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End2_Pt1stEle_PAS"              ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_End2_Pt1stEle_PAS"	           ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_End2_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End2_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End2_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End2_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End2_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End2_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End2_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_DCotTheta1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_Dist1stEle_PAS"             ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_MET_PAS"                   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_MET_PAS"                    ,    100 , 0.0     , 1.0      );

   // 1 jet 

   CreateUserTH1D( "Total_Bar_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_Bar_1Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_Bar_1Jet_Pt1stEle_PAS"          ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_Bar_1Jet_Pt1stEle_PAS"	   ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_Bar_1Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_Bar_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_Bar_1Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_Bar_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_Bar_1Jet_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_Bar_1Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_Bar_1Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_1Jet_DCotTheta1stEle_PAS"    ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_1Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_1Jet_Dist1stEle_PAS"         ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_1Jet_MET_PAS"               ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_1Jet_MET_PAS"                ,    100 , 0.0     , 1.0      );
						          										           
   CreateUserTH1D( "Total_End1_1Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End1_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End1_1Jet_Pt1stEle_PAS"         ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_End1_1Jet_Pt1stEle_PAS"	   ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_End1_1Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End1_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End1_1Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End1_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End1_1Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End1_1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End1_1Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_1Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_1Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_1Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_1Jet_MET_PAS"              ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_1Jet_MET_PAS"               ,    100 , 0.0     , 1.0      );

   CreateUserTH1D( "Total_End2_1Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End2_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End2_1Jet_Pt1stEle_PAS"         ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_End2_1Jet_Pt1stEle_PAS"	   ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_End2_1Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End2_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End2_1Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End2_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End2_1Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End2_1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End2_1Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_1Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_1Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_1Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_1Jet_MET_PAS"              ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_1Jet_MET_PAS"               ,    100 , 0.0     , 1.0      );

   // 2 jet 

   CreateUserTH1D( "Total_Bar_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_Bar_2Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_Bar_2Jet_Pt1stEle_PAS"          ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_Bar_2Jet_Pt1stEle_PAS"	   ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_Bar_2Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_Bar_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_Bar_2Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_Bar_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_Bar_2Jet_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_Bar_2Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_Bar_2Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_2Jet_DCotTheta1stEle_PAS"    ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_2Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_2Jet_Dist1stEle_PAS"         ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_2Jet_MET_PAS"               ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_2Jet_MET_PAS"                ,    100 , 0.0     , 1.0      );
						          										           
   CreateUserTH1D( "Total_End1_2Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End1_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End1_2Jet_Pt1stEle_PAS"         ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_End1_2Jet_Pt1stEle_PAS"	   ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_End1_2Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End1_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End1_2Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End1_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End1_2Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End1_2Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End1_2Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_2Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_2Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_2Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_2Jet_MET_PAS"              ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_2Jet_MET_PAS"               ,    100 , 0.0     , 1.0      );

   CreateUserTH1D( "Total_End2_2Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End2_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End2_2Jet_Pt1stEle_PAS"         ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_End2_2Jet_Pt1stEle_PAS"	   ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_End2_2Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End2_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End2_2Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End2_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End2_2Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End2_2Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End2_2Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_2Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_2Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_2Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_2Jet_MET_PAS"              ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_2Jet_MET_PAS"               ,    100 , 0.0     , 1.0      );

   // 3 jet 

   CreateUserTH1D( "Total_Bar_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_Bar_3Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_Bar_3Jet_Pt1stEle_PAS"          ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_Bar_3Jet_Pt1stEle_PAS"	   ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_Bar_3Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_Bar_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_Bar_3Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_Bar_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_Bar_3Jet_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_Bar_3Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_Bar_3Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_3Jet_DCotTheta1stEle_PAS"    ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_3Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_3Jet_Dist1stEle_PAS"         ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_3Jet_MET_PAS"               ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_3Jet_MET_PAS"                ,    100 , 0.0     , 1.0      );
						          										           
   CreateUserTH1D( "Total_End1_3Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End1_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End1_3Jet_Pt1stEle_PAS"         ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_End1_3Jet_Pt1stEle_PAS"	   ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_End1_3Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End1_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End1_3Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End1_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End1_3Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End1_3Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End1_3Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_3Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_3Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_3Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_3Jet_MET_PAS"              ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_3Jet_MET_PAS"               ,    100 , 0.0     , 1.0      );

   CreateUserTH1D( "Total_End2_3Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End2_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End2_3Jet_Pt1stEle_PAS"         ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_End2_3Jet_Pt1stEle_PAS"	   ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_End2_3Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End2_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5        );
   CreateUserTH1D( "Total_End2_3Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End2_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End2_3Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End2_3Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End2_3Jet_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_3Jet_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_3Jet_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_3Jet_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_3Jet_MET_PAS"              ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_3Jet_MET_PAS"               ,    100 , 0.0     , 1.0      );


   // Runs 160406 - 166502
   
   CreateUserTH1D( "Total1_Bar_nElectron_PAS"             ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass1_Bar_nElectron_PAS"              ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total1_Bar_Pt1stEle_PAS"              ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass1_Bar_Pt1stEle_PAS"	           ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total1_Bar_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass1_Bar_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total1_Bar_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass1_Bar_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total1_Bar_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass1_Bar_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total1_Bar_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass1_Bar_DCotTheta1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total1_Bar_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass1_Bar_Dist1stEle_PAS"             ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total1_Bar_MET_PAS"                   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass1_Bar_MET_PAS"                    ,    100 , 0.0     , 1.0      );
						          										           
   CreateUserTH1D( "Total1_End1_nElectron_PAS"            ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass1_End1_nElectron_PAS"             ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total1_End1_Pt1stEle_PAS"             ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass1_End1_Pt1stEle_PAS"	           ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total1_End1_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass1_End1_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total1_End1_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass1_End1_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total1_End1_Charge1stEle_PAS"         ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass1_End1_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total1_End1_DCotTheta1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass1_End1_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total1_End1_Dist1stEle_PAS"           ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass1_End1_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total1_End1_MET_PAS"                  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass1_End1_MET_PAS"                   ,    100 , 0.0     , 1.0      );

   CreateUserTH1D( "Total1_End2_nElectron_PAS"            ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass1_End2_nElectron_PAS"             ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total1_End2_Pt1stEle_PAS"             ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass1_End2_Pt1stEle_PAS"	           ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total1_End2_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass1_End2_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total1_End2_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass1_End2_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total1_End2_Charge1stEle_PAS"         ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass1_End2_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total1_End2_DCotTheta1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass1_End2_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total1_End2_Dist1stEle_PAS"           ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass1_End2_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total1_End2_MET_PAS"                  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass1_End2_MET_PAS"                   ,    100 , 0.0     , 1.0      );

   // Runs 166503 - onward
   
   CreateUserTH1D( "Total2_Bar_nElectron_PAS"             ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass2_Bar_nElectron_PAS"              ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total2_Bar_Pt1stEle_PAS"              ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass2_Bar_Pt1stEle_PAS"	           ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total2_Bar_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass2_Bar_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total2_Bar_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass2_Bar_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total2_Bar_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass2_Bar_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total2_Bar_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass2_Bar_DCotTheta1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total2_Bar_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass2_Bar_Dist1stEle_PAS"             ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total2_Bar_MET_PAS"                   ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass2_Bar_MET_PAS"                    ,    100 , 0.0     , 1.0      );
						          										           
   CreateUserTH1D( "Total2_End1_nElectron_PAS"            ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass2_End1_nElectron_PAS"             ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total2_End1_Pt1stEle_PAS"             ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass2_End1_Pt1stEle_PAS"	           ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total2_End1_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass2_End1_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total2_End1_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass2_End1_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total2_End1_Charge1stEle_PAS"         ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass2_End1_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total2_End1_DCotTheta1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass2_End1_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total2_End1_Dist1stEle_PAS"           ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass2_End1_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total2_End1_MET_PAS"                  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass2_End1_MET_PAS"                   ,    100 , 0.0     , 1.0      );
   
   CreateUserTH1D( "Total2_End2_nElectron_PAS"            ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass2_End2_nElectron_PAS"             ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total2_End2_Pt1stEle_PAS"             ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass2_End2_Pt1stEle_PAS"	           ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total2_End2_Eta1stEle_PAS"	          ,    100 , -5      , 5	);  CreateUserTH1D( "Pass2_End2_Eta1stEle_PAS"	           ,    100 , -5      , 5        );
   CreateUserTH1D( "Total2_End2_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass2_End2_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total2_End2_Charge1stEle_PAS"         ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass2_End2_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total2_End2_DCotTheta1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass2_End2_DCotTheta1stEle_PAS"       ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total2_End2_Dist1stEle_PAS"           ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass2_End2_Dist1stEle_PAS"            ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total2_End2_MET_PAS"                  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass2_End2_MET_PAS"                   ,    100 , 0.0     , 1.0      );

   // Pile-up [0,5]

   CreateUserTH1D( "Total_Bar_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_Bar_PU0-8_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_Bar_PU0-8_Pt1stEle_PAS"         ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_Bar_PU0-8_Pt1stEle_PAS"	   ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_Bar_PU0-8_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_Bar_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	 );
   CreateUserTH1D( "Total_Bar_PU0-8_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_Bar_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_Bar_PU0-8_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_Bar_PU0-8_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_Bar_PU0-8_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_PU0-8_DCotTheta1stEle_PAS"   ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_PU0-8_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_PU0-8_Dist1stEle_PAS"        ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_PU0-8_MET_PAS"              ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_PU0-8_MET_PAS"               ,    100 , 0.0     , 1.0      );

   CreateUserTH1D( "Total_End1_PU0-8_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End1_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End1_PU0-8_Pt1stEle_PAS"        ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_End1_PU0-8_Pt1stEle_PAS"	   ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_End1_PU0-8_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End1_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	 );
   CreateUserTH1D( "Total_End1_PU0-8_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End1_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End1_PU0-8_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End1_PU0-8_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End1_PU0-8_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_PU0-8_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_PU0-8_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_PU0-8_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_PU0-8_MET_PAS"             ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_PU0-8_MET_PAS"              ,    100 , 0.0     , 1.0      );

   CreateUserTH1D( "Total_End2_PU0-8_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End2_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End2_PU0-8_Pt1stEle_PAS"        ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_End2_PU0-8_Pt1stEle_PAS"	   ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_End2_PU0-8_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End2_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	 );
   CreateUserTH1D( "Total_End2_PU0-8_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End2_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End2_PU0-8_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End2_PU0-8_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End2_PU0-8_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_PU0-8_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_PU0-8_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_PU0-8_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_PU0-8_MET_PAS"             ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_PU0-8_MET_PAS"              ,    100 , 0.0     , 1.0      );

   // Pile-up [6-10]

   CreateUserTH1D( "Total_Bar_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_Bar_PU9-UP_nElectron_PAS"        ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_Bar_PU9-UP_Pt1stEle_PAS"        ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_Bar_PU9-UP_Pt1stEle_PAS"	   ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_Bar_PU9-UP_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_Bar_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	 );
   CreateUserTH1D( "Total_Bar_PU9-UP_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_Bar_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_Bar_PU9-UP_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_Bar_PU9-UP_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_Bar_PU9-UP_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_PU9-UP_DCotTheta1stEle_PAS"  ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_PU9-UP_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_PU9-UP_Dist1stEle_PAS"       ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_Bar_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_Bar_PU9-UP_MET_PAS"              ,    100 , 0.0     , 1.0      );

   CreateUserTH1D( "Total_End1_PU9-UP_nElectron_PAS"      ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End1_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End1_PU9-UP_Pt1stEle_PAS"       ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_End1_PU9-UP_Pt1stEle_PAS"	   ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_End1_PU9-UP_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End1_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	 );
   CreateUserTH1D( "Total_End1_PU9-UP_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End1_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End1_PU9-UP_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End1_PU9-UP_Charge1stEle_PAS"    ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End1_PU9-UP_DCotTheta1stEle_PAS",    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_PU9-UP_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_PU9-UP_Dist1stEle_PAS"     ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_PU9-UP_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End1_PU9-UP_MET_PAS"            ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End1_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1.0      );

   CreateUserTH1D( "Total_End2_PU9-UP_nElectron_PAS"      ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Pass_End2_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );
   CreateUserTH1D( "Total_End2_PU9-UP_Pt1stEle_PAS"       ,    100 , 0       , 1000     );  CreateUserTH1D( "Pass_End2_PU9-UP_Pt1stEle_PAS"	   ,    100 , 0       , 1000     );
   CreateUserTH1D( "Total_End2_PU9-UP_Eta1stEle_PAS"	  ,    100 , -5      , 5	);  CreateUserTH1D( "Pass_End2_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	 );
   CreateUserTH1D( "Total_End2_PU9-UP_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416  );  CreateUserTH1D( "Pass_End2_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416  );
   CreateUserTH1D( "Total_End2_PU9-UP_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Pass_End2_PU9-UP_Charge1stEle_PAS"    ,    2   , -1.0001 , 1.0001   );
   CreateUserTH1D( "Total_End2_PU9-UP_DCotTheta1stEle_PAS",    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_PU9-UP_DCotTheta1stEle_PAS" ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_PU9-UP_Dist1stEle_PAS"     ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_PU9-UP_Dist1stEle_PAS"      ,    100 , 0.0     , 1.0      );
   CreateUserTH1D( "Total_End2_PU9-UP_MET_PAS"            ,    100 , 0.0     , 1.0      );  CreateUserTH1D( "Pass_End2_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1.0      );

   //--------------------------------------------------------------------------
   // Loop over the chain
   //--------------------------------------------------------------------------

   if (fChain == 0) return;
   
   Long64_t nentries = fChain->GetEntries();
   // Long64_t nentries = 100;
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
     
     int    passedJSON = passJSON ( run, ls , isData ) ;

     //--------------------------------------------------------------------------
     // Do pileup re-weighting
     //--------------------------------------------------------------------------
     
     int NPILEUP_AVE = int( nPileUpInt_BX0 );
     int NPILEUP_FINAL = min( NPILEUP_AVE , 25 );
     double pileup_weight = getPileupWeight ( NPILEUP_FINAL, isData ) ;

     //--------------------------------------------------------------------------
     // Skim bug: number of jets with no loose electron overlaps is not stored
     // Work around: count the number of jets with non-zero pT ( pT > 10 is ok )
     //--------------------------------------------------------------------------

     nJetLooseEle_Stored = 0;
     if       ( JetLooseEle3_Pt > 10.0 ) nJetLooseEle_Stored = 3 ;
     else if  ( JetLooseEle2_Pt > 10.0 ) nJetLooseEle_Stored = 2 ;
     else if  ( JetLooseEle1_Pt > 10.0 ) nJetLooseEle_Stored = 1 ;
     else                                nJetLooseEle_Stored = 0 ;

     //--------------------------------------------------------------------------
     // Fill variables
     //--------------------------------------------------------------------------
     
     int min_prescale;
     int passTrigger;
     std::string min_prescale_name;

     if ( isData ) {

       //--------------------------------------------------------------------------
       // 7 trigger paths (some with multiple versions)
       //--------------------------------------------------------------------------

       // Number of times a path fired per event ( should be 0 or 1 )

       int N_Photon30_CIdVL  = 0;
       int N_Photon50_CIdVL  = 0;
       int N_Photon75_CIdVL  = 0;
       int N_Photon90_CIdVL  = 0;
       int N_Photon125       = 0;
       int N_Photon135       = 0;
       int N_Photon400       = 0;

       // Trigger prescale in an event
       
       int PS_Photon30_CIdVL = 0;
       int PS_Photon50_CIdVL = 0;
       int PS_Photon75_CIdVL = 0;
       int PS_Photon90_CIdVL = 0;
       int PS_Photon125      = 0;
       int PS_Photon135      = 0;
       int PS_Photon400      = 0;

       //--------------------------------------------------------------------------
       // Find the right prescale for this event
       //--------------------------------------------------------------------------
       
       // Did the HLT_Photon30_CaloIdVL trigger fire?

       if ( H_Photon30_CIdVL_1 > 0 && H_Photon30_CIdVL_1 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_1; } 
       if ( H_Photon30_CIdVL_2 > 0 && H_Photon30_CIdVL_2 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_2; } 
       if ( H_Photon30_CIdVL_3 > 0 && H_Photon30_CIdVL_3 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_3; } 
       if ( H_Photon30_CIdVL_4 > 0 && H_Photon30_CIdVL_4 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_4; } 
       if ( H_Photon30_CIdVL_5 > 0 && H_Photon30_CIdVL_5 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_5; } 
       if ( H_Photon30_CIdVL_6 > 0 && H_Photon30_CIdVL_6 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_6; } 
       if ( H_Photon30_CIdVL_7 > 0 && H_Photon30_CIdVL_7 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_7; } 

       if ( isData && run > 175771 ) {
	 if ( H_Photon30_CIdVL_8 > 0 && H_Photon30_CIdVL_8 != 999 ) { N_Photon30_CIdVL++; PS_Photon30_CIdVL = H_Photon30_CIdVL_8; }
       }

       // Did the HLT_Photon50_CaloIdVL trigger fire?
       
       if ( H_Photon50_CIdVL_1 > 0 && H_Photon50_CIdVL_1 != 999 ) { N_Photon50_CIdVL++; PS_Photon50_CIdVL = H_Photon50_CIdVL_1; } 
       if ( H_Photon50_CIdVL_2 > 0 && H_Photon50_CIdVL_2 != 999 ) { N_Photon50_CIdVL++; PS_Photon50_CIdVL = H_Photon50_CIdVL_2; } 
       if ( H_Photon50_CIdVL_3 > 0 && H_Photon50_CIdVL_3 != 999 ) { N_Photon50_CIdVL++; PS_Photon50_CIdVL = H_Photon50_CIdVL_3; } 
       if ( H_Photon50_CIdVL_4 > 0 && H_Photon50_CIdVL_4 != 999 ) { N_Photon50_CIdVL++; PS_Photon50_CIdVL = H_Photon50_CIdVL_4; } 

       // Did the HLT_Photon75_CaloIdVL trigger fire?
       
       if ( H_Photon75_CIdVL_1 > 0 && H_Photon75_CIdVL_1 != 999 ) { N_Photon75_CIdVL++; PS_Photon75_CIdVL = H_Photon75_CIdVL_1; } 
       if ( H_Photon75_CIdVL_2 > 0 && H_Photon75_CIdVL_2 != 999 ) { N_Photon75_CIdVL++; PS_Photon75_CIdVL = H_Photon75_CIdVL_2; } 
       if ( H_Photon75_CIdVL_3 > 0 && H_Photon75_CIdVL_3 != 999 ) { N_Photon75_CIdVL++; PS_Photon75_CIdVL = H_Photon75_CIdVL_3; } 
       if ( H_Photon75_CIdVL_4 > 0 && H_Photon75_CIdVL_4 != 999 ) { N_Photon75_CIdVL++; PS_Photon75_CIdVL = H_Photon75_CIdVL_4; } 
       if ( H_Photon75_CIdVL_5 > 0 && H_Photon75_CIdVL_5 != 999 ) { N_Photon75_CIdVL++; PS_Photon75_CIdVL = H_Photon75_CIdVL_5; } 
       if ( H_Photon75_CIdVL_6 > 0 && H_Photon75_CIdVL_6 != 999 ) { N_Photon75_CIdVL++; PS_Photon75_CIdVL = H_Photon75_CIdVL_6; } 
       if ( H_Photon75_CIdVL_7 > 0 && H_Photon75_CIdVL_7 != 999 ) { N_Photon75_CIdVL++; PS_Photon75_CIdVL = H_Photon75_CIdVL_7; }

       // Did the HLT_Photon90_CaloIdVL trigger fire?
       
       if ( H_Photon90_CIdVL_1 > 0 && H_Photon90_CIdVL_1 != 999 ) { N_Photon90_CIdVL++; PS_Photon90_CIdVL = H_Photon90_CIdVL_1; } 
       if ( H_Photon90_CIdVL_2 > 0 && H_Photon90_CIdVL_2 != 999 ) { N_Photon90_CIdVL++; PS_Photon90_CIdVL = H_Photon90_CIdVL_2; } 
       if ( H_Photon90_CIdVL_3 > 0 && H_Photon90_CIdVL_3 != 999 ) { N_Photon90_CIdVL++; PS_Photon90_CIdVL = H_Photon90_CIdVL_3; } 
       if ( H_Photon90_CIdVL_4 > 0 && H_Photon90_CIdVL_4 != 999 ) { N_Photon90_CIdVL++; PS_Photon90_CIdVL = H_Photon90_CIdVL_4; } 

       // Did the HLT_Photon125 trigger fire?
       
       if ( H_Photon125_1      > 0 && H_Photon125_1      != 999 ) { N_Photon125     ++; PS_Photon125      = H_Photon125_1     ; } 
       if ( H_Photon125_2      > 0 && H_Photon125_2      != 999 ) { N_Photon125     ++; PS_Photon125      = H_Photon125_2     ; } 
   
       // Did the HLT_Photon135 trigger fire?

       if ( H_Photon135_1      > 0 && H_Photon135_1      != 999 ) { N_Photon135     ++; PS_Photon135      = H_Photon135_1     ; } 
       if ( H_Photon135_2      > 0 && H_Photon135_2      != 999 ) { N_Photon135     ++; PS_Photon135      = H_Photon135_2     ; } 
       
       // Did the HLT_Photon400 trigger fire?

       if ( H_Photon400_1      > 0 && H_Photon400_1      != 999 ) { N_Photon400     ++; PS_Photon400      = H_Photon400_1     ; } 

       if ( isData && run > 175771 ) {
	 if ( H_Photon400_2      > 0 && H_Photon400_2      != 999 ) { N_Photon400     ++; PS_Photon400      = H_Photon400_2     ; }
       }
       
       // Sanity check: make sure two versions of the same trigger didn't fire in the same event (impossible)
              
       if ( N_Photon30_CIdVL > 1 ) { std::cout << "ERROR: trigger overlap in N_Photon30_CIdVL" << std::endl; exit (0); }
       if ( N_Photon50_CIdVL > 1 ) { std::cout << "ERROR: trigger overlap in N_Photon50_CIdVL" << std::endl; exit (0); }
       if ( N_Photon75_CIdVL > 1 ) { std::cout << "ERROR: trigger overlap in N_Photon75_CIdVL" << std::endl; exit (0); }
       if ( N_Photon90_CIdVL > 1 ) { std::cout << "ERROR: trigger overlap in N_Photon90_CIdVL" << std::endl; exit (0); }
       if ( N_Photon125      > 1 ) { std::cout << "ERROR: trigger overlap in N_Photon125"      << std::endl; exit (0); }
       if ( N_Photon135      > 1 ) { std::cout << "ERROR: trigger overlap in N_Photon135"      << std::endl; exit (0); }
       if ( N_Photon400      > 1 ) { std::cout << "ERROR: trigger overlap in N_Photon400"      << std::endl; exit (0); }
      
       // What is the lowest-prescale trigger that this electron could have fired?
 
       min_prescale      = 999999;
       min_prescale_name = std::string("");

       if ( N_Photon30_CIdVL != 0 && QCDFakeEle1_Pt > 30. * trigger_tolerance  && PS_Photon30_CIdVL <= min_prescale ) { min_prescale = PS_Photon30_CIdVL; min_prescale_name = std::string("PS_Photon30_CIdVL"); }
       if ( N_Photon50_CIdVL != 0 && QCDFakeEle1_Pt > 50. * trigger_tolerance  && PS_Photon50_CIdVL <= min_prescale ) { min_prescale = PS_Photon50_CIdVL; min_prescale_name = std::string("PS_Photon50_CIdVL"); }
       if ( N_Photon75_CIdVL != 0 && QCDFakeEle1_Pt > 75. * trigger_tolerance  && PS_Photon75_CIdVL <= min_prescale ) { min_prescale = PS_Photon75_CIdVL; min_prescale_name = std::string("PS_Photon75_CIdVL"); }
       if ( N_Photon90_CIdVL != 0 && QCDFakeEle1_Pt > 90. * trigger_tolerance  && PS_Photon90_CIdVL <= min_prescale ) { min_prescale = PS_Photon90_CIdVL; min_prescale_name = std::string("PS_Photon90_CIdVL"); }
       if ( N_Photon125      != 0 && QCDFakeEle1_Pt > 125.* trigger_tolerance  && PS_Photon125      <= min_prescale ) { min_prescale = PS_Photon125     ; min_prescale_name = std::string("PS_Photon125"     ); }
       if ( N_Photon135      != 0 && QCDFakeEle1_Pt > 135.* trigger_tolerance  && PS_Photon135      <= min_prescale ) { min_prescale = PS_Photon135     ; min_prescale_name = std::string("PS_Photon135"     ); }
       if ( N_Photon400      != 0 && QCDFakeEle1_Pt > 400.* trigger_tolerance  && PS_Photon400      <= min_prescale ) { min_prescale = PS_Photon400     ; min_prescale_name = std::string("PS_Photon400"     ); }

       // If we find a suitable trigger, scale this event by that trigger's prescale

       passTrigger = 0;
       if ( min_prescale != 999999 ) {
	 passTrigger = 1;     
       } else {
	 min_prescale = 0;
       }
     }  // end if (isData) 
					       
     else { 
       min_prescale = 1;
       passTrigger = 1 ;
     }


     /*
       std::cout << "------------------------------------" << std::endl;
       std::cout << "WARN : more than one trigger fire in jentry " << jentry << std::endl;
       std::cout << "----------------------" << std::endl;
       std::cout << "Lead GSF electron = " << QCDFakeEle1_Pt    << std::endl;
       std::cout << "----------------------" << std::endl;
       std::cout << "PS_Photon30_CIdVL = " << PS_Photon30_CIdVL << std::endl;
       std::cout << "PS_Photon50_CIdVL = " << PS_Photon50_CIdVL << std::endl;
       std::cout << "PS_Photon75_CIdVL = " << PS_Photon75_CIdVL << std::endl;
       std::cout << "PS_Photon90_CIdVL = " << PS_Photon90_CIdVL << std::endl;
       std::cout << "PS_Photon125      = " << PS_Photon125      << std::endl;
       std::cout << "PS_Photon135      = " << PS_Photon135      << std::endl;
       std::cout << "PS_Photon400      = " << PS_Photon400      << std::endl;
       std::cout << "----------------------" << std::endl;
       std::cout << "Selected prescale = " << min_prescale_name << ", " << min_prescale << std::endl;
     */

     // trigger

     //--------------------------------------------------------------------------
     // Fill cut values
     //--------------------------------------------------------------------------

     fillVariableWithValue(   "PassTrigger"                   , passTrigger            , min_prescale  ); 
     
     // JSON variable
     fillVariableWithValue(   "PassJSON"                      , passedJSON             , min_prescale  ); 
     										       
     // Filters									       
     fillVariableWithValue(   "PassHBHENoiseFilter"           , PassHBHENoiseFilter    , min_prescale  );
     fillVariableWithValue(   "PassBeamHaloFilterTight"       , PassBeamHaloFilterTight, min_prescale  );
									      	       
     // Muon variables ( for veto ) 					      	       
     fillVariableWithValue(   "nMuon"                         , nMuon_Ana              , min_prescale  );
			                                      		               
     // 1st Electron variables				      		               
     fillVariableWithValue(   "nEle"                          , nEle_QCDFake           , min_prescale  ); 
     fillVariableWithValue(   "Pt1stEle"                      , QCDFakeEle1_Pt         , min_prescale  );

     // 1st JET variables                                     		               
     fillVariableWithValue(   "nJet"                          , nJetLooseEle_Stored    , min_prescale  );

     // MET
     fillVariableWithValue(   "MET"                           , MET_Pt                 , min_prescale  );

     // Dummy variables
     fillVariableWithValue(   "denominator"                   , 1                      , min_prescale  );

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

     if( fabs( QCDFakeEle1_Eta  ) < eleEta_bar )        isBarrel = true;
     if( fabs( QCDFakeEle1_Eta  ) > eleEta_end1_min &&
	 fabs( QCDFakeEle1_Eta  ) < eleEta_end1_max )   isEndcap1 = true;
     if( fabs( QCDFakeEle1_Eta  ) > eleEta_end2_min &&
	 fabs( QCDFakeEle1_Eta  ) < eleEta_end2_max )   isEndcap2 = true;
     
     if ( passed_denominator ) { 
       
       if ( isBarrel ) {
	 FillUserTH1D( "Total_Bar_nElectron_PAS"         , nEle_QCDFake                   , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_Bar_Pt1stEle_PAS"	         , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_Bar_Eta1stEle_PAS"	 , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_Bar_Phi1stEle_PAS"	 , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_Bar_Charge1stEle_PAS"      , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	 FillUserTH1D( "Total_Bar_DCotTheta1stEle_PAS"   , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_Bar_Dist1stEle_PAS"        , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_Bar_MET_PAS"               , MET_Pt                         , pileup_weight * min_prescale );

	 if ( nJetLooseEle_Stored >= 1 ) {
	   FillUserTH1D( "Total_Bar_1Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_1Jet_Pt1stEle_PAS"	      , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_1Jet_Eta1stEle_PAS"	      , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_1Jet_Phi1stEle_PAS"	      , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_1Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total_Bar_1Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_1Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_1Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	 }

	 if ( nJetLooseEle_Stored >= 2 ) {
	   FillUserTH1D( "Total_Bar_2Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_2Jet_Pt1stEle_PAS"	      , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_2Jet_Eta1stEle_PAS"	      , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_2Jet_Phi1stEle_PAS"	      , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_2Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total_Bar_2Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_2Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_2Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	 }

	 if ( nJetLooseEle_Stored >= 3 ) {
	   FillUserTH1D( "Total_Bar_3Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_3Jet_Pt1stEle_PAS"	      , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_3Jet_Eta1stEle_PAS"	      , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_3Jet_Phi1stEle_PAS"	      , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_3Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total_Bar_3Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_3Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_3Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	 }

	 if ( isData && run <= 166502 ) { 
	   FillUserTH1D( "Total1_Bar_nElectron_PAS"               , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_Bar_Pt1stEle_PAS"                , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_Bar_Eta1stEle_PAS"	          , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_Bar_Phi1stEle_PAS"	          , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_Bar_Charge1stEle_PAS"            , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total1_Bar_DCotTheta1stEle_PAS"         , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_Bar_Dist1stEle_PAS"              , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_Bar_MET_PAS"                     , MET_Pt                         , pileup_weight * min_prescale );
	 } 

	 if ( isData && run >= 166503 ) { 
	   FillUserTH1D( "Total2_Bar_nElectron_PAS"               , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_Bar_Pt1stEle_PAS"                , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_Bar_Eta1stEle_PAS"	          , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_Bar_Phi1stEle_PAS"	          , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_Bar_Charge1stEle_PAS"            , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total2_Bar_DCotTheta1stEle_PAS"         , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_Bar_Dist1stEle_PAS"              , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
 	   FillUserTH1D( "Total2_Bar_MET_PAS"                     , MET_Pt                         , pileup_weight * min_prescale );
	 } 

	 if ( nVertex_good >= 0 && nVertex_good <= 5 ){
	   FillUserTH1D( "Total_Bar_PU0-8_nElectron_PAS"          , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_PU0-8_Pt1stEle_PAS"           , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_PU0-8_Eta1stEle_PAS"          , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_PU0-8_Phi1stEle_PAS"          , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_PU0-8_Charge1stEle_PAS"       , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total_Bar_PU0-8_DCotTheta1stEle_PAS"    , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_PU0-8_Dist1stEle_PAS"         , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_PU0-8_MET_PAS"                , MET_Pt                         , pileup_weight * min_prescale );
	 }						          
							          
	 if ( nVertex_good >= 6 && nVertex_good <= 10 ){      
	   FillUserTH1D( "Total_Bar_PU9-UP_nElectron_PAS"         , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_PU9-UP_Pt1stEle_PAS"          , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_PU9-UP_Eta1stEle_PAS"         , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_PU9-UP_Phi1stEle_PAS"         , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_PU9-UP_Charge1stEle_PAS"      , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total_Bar_PU9-UP_DCotTheta1stEle_PAS"   , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_PU9-UP_Dist1stEle_PAS"        , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_Bar_PU9-UP_MET_PAS"               , MET_Pt                         , pileup_weight * min_prescale );
	 }
       }
       
       if ( isEndcap1 ) {
	 FillUserTH1D( "Total_End1_nElectron_PAS"         , nEle_QCDFake                   , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_End1_Pt1stEle_PAS"          , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_End1_Eta1stEle_PAS"	  , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_End1_Phi1stEle_PAS"	  , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_End1_Charge1stEle_PAS"      , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	 FillUserTH1D( "Total_End1_DCotTheta1stEle_PAS"   , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_End1_Dist1stEle_PAS"        , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_End1_MET_PAS"               , MET_Pt                         , pileup_weight * min_prescale );
	 

	 if ( nJetLooseEle_Stored >= 1 ) {
	   FillUserTH1D( "Total_End1_1Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_1Jet_Pt1stEle_PAS"	       , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_1Jet_Eta1stEle_PAS"       , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_1Jet_Phi1stEle_PAS"       , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_1Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total_End1_1Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_1Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_1Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	 }

	 if ( nJetLooseEle_Stored >= 2 ) {
	   FillUserTH1D( "Total_End1_2Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_2Jet_Pt1stEle_PAS"	       , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_2Jet_Eta1stEle_PAS"       , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_2Jet_Phi1stEle_PAS"       , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_2Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total_End1_2Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_2Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_2Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	 }

	 if ( nJetLooseEle_Stored >= 3 ) {
	   FillUserTH1D( "Total_End1_3Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_3Jet_Pt1stEle_PAS"	       , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_3Jet_Eta1stEle_PAS"       , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_3Jet_Phi1stEle_PAS"       , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_3Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total_End1_3Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_3Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_3Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	 }


	 if ( isData && run <= 166502 ) { 
	   FillUserTH1D( "Total1_End1_nElectron_PAS"              , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_End1_Pt1stEle_PAS"               , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_End1_Eta1stEle_PAS"	          , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_End1_Phi1stEle_PAS"	          , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_End1_Charge1stEle_PAS"           , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total1_End1_DCotTheta1stEle_PAS"        , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_End1_Dist1stEle_PAS"             , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_End1_MET_PAS"                    , MET_Pt                         , pileup_weight * min_prescale );
	 } 

	 if ( isData && run >= 166503 ) { 
	   FillUserTH1D( "Total2_End1_nElectron_PAS"              , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_End1_Pt1stEle_PAS"               , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_End1_Eta1stEle_PAS"	          , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_End1_Phi1stEle_PAS"	          , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_End1_Charge1stEle_PAS"           , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total2_End1_DCotTheta1stEle_PAS"        , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_End1_Dist1stEle_PAS"             , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_End1_MET_PAS"                    , MET_Pt                         , pileup_weight * min_prescale );
	 } 

	 if ( nVertex_good >= 0 && nVertex_good <= 5 ){
	   FillUserTH1D( "Total_End1_PU0-8_nElectron_PAS"         , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_PU0-8_Pt1stEle_PAS"          , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_PU0-8_Eta1stEle_PAS"         , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_PU0-8_Phi1stEle_PAS"         , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_PU0-8_Charge1stEle_PAS"      , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total_End1_PU0-8_DCotTheta1stEle_PAS"   , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_PU0-8_Dist1stEle_PAS"        , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_PU0-8_MET_PAS"               , MET_Pt                         , pileup_weight * min_prescale );
	 }						          
							          
	 if ( nVertex_good >= 6 && nVertex_good <= 10 ){      
	   FillUserTH1D( "Total_End1_PU9-UP_nElectron_PAS"        , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_PU9-UP_Pt1stEle_PAS"         , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_PU9-UP_Eta1stEle_PAS"        , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_PU9-UP_Phi1stEle_PAS"        , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_PU9-UP_Charge1stEle_PAS"     , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total_End1_PU9-UP_DCotTheta1stEle_PAS"  , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_PU9-UP_Dist1stEle_PAS"       , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End1_PU9-UP_MET_PAS"              , MET_Pt                         , pileup_weight * min_prescale );
	 }
       }
       
       if ( isEndcap2 ) {
	 FillUserTH1D( "Total_End2_nElectron_PAS"         , nEle_QCDFake                   , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_End2_Pt1stEle_PAS"          , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_End2_Eta1stEle_PAS"	  , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_End2_Phi1stEle_PAS"	  , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_End2_Charge1stEle_PAS"      , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	 FillUserTH1D( "Total_End2_DCotTheta1stEle_PAS"   , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_End2_Dist1stEle_PAS"        , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	 FillUserTH1D( "Total_End2_MET_PAS"               , MET_Pt                         , pileup_weight * min_prescale );
	 
	 if ( nJetLooseEle_Stored >= 1 ) {
	   FillUserTH1D( "Total_End2_1Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_1Jet_Pt1stEle_PAS"	       , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_1Jet_Eta1stEle_PAS"       , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_1Jet_Phi1stEle_PAS"       , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_1Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total_End2_1Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_1Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_1Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	 }

	 if ( nJetLooseEle_Stored >= 2 ) {
	   FillUserTH1D( "Total_End2_2Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_2Jet_Pt1stEle_PAS"	       , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_2Jet_Eta1stEle_PAS"       , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_2Jet_Phi1stEle_PAS"       , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_2Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total_End2_2Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_2Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_2Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	 }

	 if ( nJetLooseEle_Stored >= 3 ) {
	   FillUserTH1D( "Total_End2_3Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_3Jet_Pt1stEle_PAS"	       , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_3Jet_Eta1stEle_PAS"       , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_3Jet_Phi1stEle_PAS"       , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_3Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total_End2_3Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_3Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_3Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	 }
	 
	 if ( isData && run <= 166502 ) { 
	   FillUserTH1D( "Total1_End2_nElectron_PAS"              , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_End2_Pt1stEle_PAS"               , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_End2_Eta1stEle_PAS"	          , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_End2_Phi1stEle_PAS"	          , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_End2_Charge1stEle_PAS"           , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total1_End2_DCotTheta1stEle_PAS"        , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_End2_Dist1stEle_PAS"             , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total1_End2_MET_PAS"                    , MET_Pt                         , pileup_weight * min_prescale );
	 } 

	 if ( isData && run >= 166503 ) { 
	   FillUserTH1D( "Total2_End2_nElectron_PAS"              , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_End2_Pt1stEle_PAS"               , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_End2_Eta1stEle_PAS"	          , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_End2_Phi1stEle_PAS"	          , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_End2_Charge1stEle_PAS"           , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total2_End2_DCotTheta1stEle_PAS"        , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_End2_Dist1stEle_PAS"             , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total2_End2_MET_PAS"                    , MET_Pt                         , pileup_weight * min_prescale );
	 } 

	 if ( nVertex_good >= 0 && nVertex_good <= 5 ){
	   FillUserTH1D( "Total_End2_PU0-8_nElectron_PAS"         , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_PU0-8_Pt1stEle_PAS"          , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_PU0-8_Eta1stEle_PAS"         , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_PU0-8_Phi1stEle_PAS"         , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_PU0-8_Charge1stEle_PAS"      , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total_End2_PU0-8_DCotTheta1stEle_PAS"   , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_PU0-8_Dist1stEle_PAS"        , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_PU0-8_MET_PAS"               , MET_Pt                         , pileup_weight * min_prescale );
	 }						          
							          
	 if ( nVertex_good >= 6 && nVertex_good <= 10 ){      
	   FillUserTH1D( "Total_End2_PU9-UP_nElectron_PAS"        , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_PU9-UP_Pt1stEle_PAS"         , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_PU9-UP_Eta1stEle_PAS"        , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_PU9-UP_Phi1stEle_PAS"        , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_PU9-UP_Charge1stEle_PAS"     , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	   FillUserTH1D( "Total_End2_PU9-UP_DCotTheta1stEle_PAS"  , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_PU9-UP_Dist1stEle_PAS"       , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Total_End2_PU9-UP_MET_PAS"              , MET_Pt                         , pileup_weight * min_prescale );
	 }
       }
       
       if ( QCDFakeEle1_PassID == 1 ) { 
	 
	 if ( isBarrel ) { 
	   FillUserTH1D( "Pass_Bar_nElectron_PAS"         , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_Bar_Pt1stEle_PAS"	  , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_Bar_Eta1stEle_PAS"	  , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_Bar_Phi1stEle_PAS"	  , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_Bar_Charge1stEle_PAS"      , QCDFakeEle1_Charge             , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_Bar_DCotTheta1stEle_PAS"   , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_Bar_Dist1stEle_PAS"        , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_Bar_MET_PAS"               , MET_Pt                         , pileup_weight * min_prescale );

	   if ( nJetLooseEle_Stored >= 1 ) {
	     FillUserTH1D( "Pass_Bar_1Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_1Jet_Pt1stEle_PAS"	       , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_1Jet_Eta1stEle_PAS"       , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_1Jet_Phi1stEle_PAS"       , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_1Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_1Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_1Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_1Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	   }

	   if ( nJetLooseEle_Stored >= 2 ) {
	     FillUserTH1D( "Pass_Bar_2Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_2Jet_Pt1stEle_PAS"	       , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_2Jet_Eta1stEle_PAS"       , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_2Jet_Phi1stEle_PAS"       , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_2Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_2Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_2Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_2Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	   }

	   if ( nJetLooseEle_Stored >= 3 ) {
	     FillUserTH1D( "Pass_Bar_3Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_3Jet_Pt1stEle_PAS"	       , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_3Jet_Eta1stEle_PAS"       , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_3Jet_Phi1stEle_PAS"       , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_3Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_3Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_3Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_3Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	   }
	   
	   if ( isData && run <= 166502 ) { 
	     FillUserTH1D( "Pass1_Bar_nElectron_PAS"               , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_Bar_Pt1stEle_PAS"                , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_Bar_Eta1stEle_PAS"	           , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_Bar_Phi1stEle_PAS"	           , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_Bar_Charge1stEle_PAS"            , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	     FillUserTH1D( "Pass1_Bar_DCotTheta1stEle_PAS"         , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_Bar_Dist1stEle_PAS"              , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_Bar_MET_PAS"                     , MET_Pt                         , pileup_weight * min_prescale );
	   } 
	   
	   if ( isData && run >= 166503 ) { 
	     FillUserTH1D( "Pass2_Bar_nElectron_PAS"               , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_Bar_Pt1stEle_PAS"                , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_Bar_Eta1stEle_PAS"	           , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_Bar_Phi1stEle_PAS"	           , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_Bar_Charge1stEle_PAS"            , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	     FillUserTH1D( "Pass2_Bar_DCotTheta1stEle_PAS"         , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_Bar_Dist1stEle_PAS"              , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_Bar_MET_PAS"                     , MET_Pt                         , pileup_weight * min_prescale );
	   } 

	   if ( nVertex_good >= 0 && nVertex_good <= 5 ){
	     FillUserTH1D( "Pass_Bar_PU0-8_nElectron_PAS"          , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_PU0-8_Pt1stEle_PAS"           , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_PU0-8_Eta1stEle_PAS"          , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_PU0-8_Phi1stEle_PAS"          , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_PU0-8_Charge1stEle_PAS"       , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	     FillUserTH1D( "Pass_Bar_PU0-8_DCotTheta1stEle_PAS"    , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_PU0-8_Dist1stEle_PAS"         , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_PU0-8_MET_PAS"                , MET_Pt                         , pileup_weight * min_prescale );
	   }						          
	   
	   if ( nVertex_good >= 6 && nVertex_good <= 10 ){      
	     FillUserTH1D( "Pass_Bar_PU9-UP_nElectron_PAS"         , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_PU9-UP_Pt1stEle_PAS"          , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_PU9-UP_Eta1stEle_PAS"         , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_PU9-UP_Phi1stEle_PAS"         , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_PU9-UP_Charge1stEle_PAS"      , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	     FillUserTH1D( "Pass_Bar_PU9-UP_DCotTheta1stEle_PAS"   , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_PU9-UP_Dist1stEle_PAS"        , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_Bar_PU9-UP_MET_PAS"               , MET_Pt                         , pileup_weight * min_prescale );
	   }	  
	 }
	 
	 if ( isEndcap1 ) { 
	   FillUserTH1D( "Pass_End1_nElectron_PAS"        , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_End1_Pt1stEle_PAS"	  , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_End1_Eta1stEle_PAS"	  , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_End1_Phi1stEle_PAS"	  , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_End1_Charge1stEle_PAS"     , QCDFakeEle1_Charge             , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_End1_DCotTheta1stEle_PAS"  , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_End1_Dist1stEle_PAS"       , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_End1_MET_PAS"              , MET_Pt                         , pileup_weight * min_prescale );
	   
	   if ( nJetLooseEle_Stored >= 1 ) {
	     FillUserTH1D( "Pass_End1_1Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_1Jet_Pt1stEle_PAS"        , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_1Jet_Eta1stEle_PAS"       , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_1Jet_Phi1stEle_PAS"       , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_1Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_1Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_1Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_1Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	   }

	   if ( nJetLooseEle_Stored >= 2 ) {
	     FillUserTH1D( "Pass_End1_2Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_2Jet_Pt1stEle_PAS"	, QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_2Jet_Eta1stEle_PAS"       , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_2Jet_Phi1stEle_PAS"       , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_2Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_2Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_2Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_2Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	   }

	   if ( nJetLooseEle_Stored >= 3 ) {
	     FillUserTH1D( "Pass_End1_3Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_3Jet_Pt1stEle_PAS"	, QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_3Jet_Eta1stEle_PAS"       , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_3Jet_Phi1stEle_PAS"       , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_3Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_3Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_3Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_3Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	   }

	   if ( isData && run <= 166502 ) { 
	     FillUserTH1D( "Pass1_End1_nElectron_PAS"              , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_End1_Pt1stEle_PAS"               , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_End1_Eta1stEle_PAS"	           , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_End1_Phi1stEle_PAS"	           , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_End1_Charge1stEle_PAS"           , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	     FillUserTH1D( "Pass1_End1_DCotTheta1stEle_PAS"        , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_End1_Dist1stEle_PAS"             , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_End1_MET_PAS"                    , MET_Pt                         , pileup_weight * min_prescale );
	   } 
	   
	   if ( isData && run >= 166503 ) { 
	     FillUserTH1D( "Pass2_End1_nElectron_PAS"               , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_End1_Pt1stEle_PAS"                , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_End1_Eta1stEle_PAS"	           , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_End1_Phi1stEle_PAS"	           , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_End1_Charge1stEle_PAS"           , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	     FillUserTH1D( "Pass2_End1_DCotTheta1stEle_PAS"        , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_End1_Dist1stEle_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	   } 

	   if ( nVertex_good >= 0 && nVertex_good <= 5 ){
	     FillUserTH1D( "Pass_End1_PU0-8_nElectron_PAS"         , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_PU0-8_Pt1stEle_PAS"          , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_PU0-8_Eta1stEle_PAS"         , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_PU0-8_Phi1stEle_PAS"         , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_PU0-8_Charge1stEle_PAS"      , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	     FillUserTH1D( "Pass_End1_PU0-8_DCotTheta1stEle_PAS"   , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_PU0-8_Dist1stEle_PAS"        , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_PU0-8_MET_PAS"               , MET_Pt                         , pileup_weight * min_prescale );
	   }						          
	   
	   if ( nVertex_good >= 6 && nVertex_good <= 10 ){      
	     FillUserTH1D( "Pass_End1_PU9-UP_nElectron_PAS"        , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_PU9-UP_Pt1stEle_PAS"         , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_PU9-UP_Eta1stEle_PAS"        , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_PU9-UP_Phi1stEle_PAS"        , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_PU9-UP_Charge1stEle_PAS"     , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	     FillUserTH1D( "Pass_End1_PU9-UP_DCotTheta1stEle_PAS"  , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_PU9-UP_Dist1stEle_PAS"       , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End1_PU9-UP_MET_PAS"              , MET_Pt                         , pileup_weight * min_prescale );
	   }
	 }

	 if ( isEndcap2 ) { 
	   FillUserTH1D( "Pass_End2_nElectron_PAS"        , nEle_QCDFake                   , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_End2_Pt1stEle_PAS"	  , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_End2_Eta1stEle_PAS"	  , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_End2_Phi1stEle_PAS"	  , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_End2_Charge1stEle_PAS"     , QCDFakeEle1_Charge             , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_End2_DCotTheta1stEle_PAS"  , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_End2_Dist1stEle_PAS"       , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	   FillUserTH1D( "Pass_End2_MET_PAS"              , MET_Pt                         , pileup_weight * min_prescale );
	   
	   if ( nJetLooseEle_Stored >= 1 ) {
	     FillUserTH1D( "Pass_End2_1Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_1Jet_Pt1stEle_PAS"        , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_1Jet_Eta1stEle_PAS"       , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_1Jet_Phi1stEle_PAS"       , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_1Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_1Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_1Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_1Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	   }

	   if ( nJetLooseEle_Stored >= 2 ) {
	     FillUserTH1D( "Pass_End2_2Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_2Jet_Pt1stEle_PAS"	, QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_2Jet_Eta1stEle_PAS"       , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_2Jet_Phi1stEle_PAS"       , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_2Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_2Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_2Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_2Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	   }

	   if ( nJetLooseEle_Stored >= 3 ) {
	     FillUserTH1D( "Pass_End2_3Jet_nElectron_PAS"       , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_3Jet_Pt1stEle_PAS"	, QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_3Jet_Eta1stEle_PAS"       , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_3Jet_Phi1stEle_PAS"       , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_3Jet_Charge1stEle_PAS"    , QCDFakeEle1_Charge             , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_3Jet_DCotTheta1stEle_PAS" , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_3Jet_Dist1stEle_PAS"      , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_3Jet_MET_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	   }
	   
	   if ( isData && run <= 166502 ) { 
	     FillUserTH1D( "Pass1_End2_nElectron_PAS"              , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_End2_Pt1stEle_PAS"               , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_End2_Eta1stEle_PAS"	           , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_End2_Phi1stEle_PAS"	           , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_End2_Charge1stEle_PAS"           , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	     FillUserTH1D( "Pass1_End2_DCotTheta1stEle_PAS"        , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_End2_Dist1stEle_PAS"             , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass1_End2_MET_PAS"                    , MET_Pt                         , pileup_weight * min_prescale );
	   } 
	   
	   if ( isData && run >= 166503 ) { 
	     FillUserTH1D( "Pass2_End2_nElectron_PAS"               , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_End2_Pt1stEle_PAS"                , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_End2_Eta1stEle_PAS"	           , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_End2_Phi1stEle_PAS"	           , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_End2_Charge1stEle_PAS"           , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	     FillUserTH1D( "Pass2_End2_DCotTheta1stEle_PAS"        , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass2_End2_Dist1stEle_PAS"             , MET_Pt                         , pileup_weight * min_prescale );
	   } 

	   if ( nVertex_good >= 0 && nVertex_good <= 5 ){
	     FillUserTH1D( "Pass_End2_PU0-8_nElectron_PAS"         , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_PU0-8_Pt1stEle_PAS"          , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_PU0-8_Eta1stEle_PAS"         , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_PU0-8_Phi1stEle_PAS"         , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_PU0-8_Charge1stEle_PAS"      , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	     FillUserTH1D( "Pass_End2_PU0-8_DCotTheta1stEle_PAS"   , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_PU0-8_Dist1stEle_PAS"        , MET_Pt                         , pileup_weight * min_prescale );
	   }						          
	   
	   if ( nVertex_good >= 6 && nVertex_good <= 10 ){      
	     FillUserTH1D( "Pass_End2_PU9-UP_nElectron_PAS"        , nEle_QCDFake                   , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_PU9-UP_Pt1stEle_PAS"         , QCDFakeEle1_Pt                 , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_PU9-UP_Eta1stEle_PAS"        , QCDFakeEle1_Eta                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_PU9-UP_Phi1stEle_PAS"        , QCDFakeEle1_Phi                , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_PU9-UP_Charge1stEle_PAS"     , QCDFakeEle1_Charge             , pileup_weight * min_prescale );  
	     FillUserTH1D( "Pass_End2_PU9-UP_DCotTheta1stEle_PAS"  , QCDFakeEle1_DCotTheta          , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_PU9-UP_Dist1stEle_PAS"       , QCDFakeEle1_Dist               , pileup_weight * min_prescale );
	     FillUserTH1D( "Pass_End2_PU9-UP_MET_PAS"              , MET_Pt                         , pileup_weight * min_prescale );
	   }
	 }
       }
     }
   } // End loop over events

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
