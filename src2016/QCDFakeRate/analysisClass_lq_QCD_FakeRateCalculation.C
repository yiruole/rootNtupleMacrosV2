#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include "include/HistoReader.h"
// for scale factors
#include "ElectronScaleFactors.C"
// for prescales
#include "include/Run2PhotonTriggerPrescales.h"

bool isHEMElectron(float eta, float phi) {
  if(eta <= -1.3 && eta >= -3.0)
    if(phi <= -0.87 && phi >= -1.57)
      return true;
  return false;
}

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
  // Analysis year
  //--------------------------------------------------------------------------
  int analysisYear = getPreCutValue1("AnalysisYear");

  //--------------------------------------------------------------------------
  // reco scale factors
  //--------------------------------------------------------------------------
  std::string recoSFFileName = getPreCutString1("RecoSFFileName");
  std::unique_ptr<HistoReader> recoScaleFactorReader = std::unique_ptr<HistoReader>(new HistoReader(recoSFFileName,"EGamma_SF2D","EGamma_SF2D",true,false));

  //--------------------------------------------------------------------------
  // Photon trigger average prescales
  //--------------------------------------------------------------------------
  Run2PhotonTriggerPrescales run2PhotonTriggerPrescales;

  //--------------------------------------------------------------------------
  // Create TH1D's
  //--------------------------------------------------------------------------

  // Debugging

  CreateUserTH1D("Total_minPrescale"                     ,   200 , 0.0      , 2000.0   ); CreateUserTH1D("Electrons_minPrescale"                       ,   200 , 0.0      , 2000.0   ); CreateUserTH1D("Jets_minPrescale"                       ,   200 , 0.0      , 2000.0   );
  CreateUserTH1D("Total_pileupWeight"                    ,   100 , 0.0      , 2.0      ); CreateUserTH1D("Electrons_pileupWeight"                      ,   100 , 0.0      , 2.0      ); CreateUserTH1D("Jets_pileupWeight"                      ,   100 , 0.0      , 2.0      );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D("Total_nVertex_PAS"                 ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Electrons_nVertex_PAS"                   ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Jets_nVertex_PAS"                   ,   100 , 0.0      , 100.0    ); 
  CreateUserTH1D("Total_Bar_nVertex_PAS"                 ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Electrons_Bar_nVertex_PAS"                   ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Jets_Bar_nVertex_PAS"                   ,   100 , 0.0      , 100.0    ); 
  CreateUserTH1D("Total_End1_nVertex_PAS"                ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Electrons_End1_nVertex_PAS"                  ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Jets_End1_nVertex_PAS"                  ,   100 , 0.0      , 100.0    ); 
  CreateUserTH1D("Total_End2_nVertex_PAS"                ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Electrons_End2_nVertex_PAS"                  ,   100 , 0.0      , 100.0    ); CreateUserTH1D("Jets_End2_nVertex_PAS"                  ,   100 , 0.0      , 100.0    ); 
                                                                                                                                                                                                                                                                                     
  // inclusive                                                                                                                                                                                                                                                                       

  CreateUserTH1D( "Total_nElectron_PAS"              ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_nElectron_PAS"               ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_nElectron_PAS"               ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_Pt1stEle_PAS"               ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "Jets_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
  CreateUserTH1D( "Total_hltPt1stEle_PAS"            ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );  CreateUserTH1D( "Jets_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );
  CreateUserTH1D( "Total_Eta1stEle_PAS"	          ,    100 , -5      , 5	         );  CreateUserTH1D( "Electrons_Eta1stEle_PAS"	           ,    100 , -5      , 5          );  CreateUserTH1D( "Jets_Eta1stEle_PAS"	           ,    100 , -5      , 5          );
  CreateUserTH1D( "Total_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416     );  CreateUserTH1D( "Electrons_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416    );  CreateUserTH1D( "Jets_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416    );
  CreateUserTH1D( "Total_Charge1stEle_PAS"	          ,    2   , -1.0001 , 1.0001  );  CreateUserTH1D( "Electrons_Charge1stEle_PAS"            ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_Charge1stEle_PAS"            ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_MET_PAS"                    ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_MET_PAS"                     ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_MET_PAS"                     ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_TrkIsoHEEP7_PAS"            ,    200 , 0.0     , 100.     );  CreateUserTH1D( "Electrons_TrkIsoHEEP7_PAS"              ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Jets_TrkIsoHEEP7_PAS"              ,    200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_TrkIsoHEEP7vsPt_PAS"        ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_TrkIsoHEEP7vsHLTPt_PAS"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_Bar_nElectron_PAS"              ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_Bar_nElectron_PAS"               ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_Bar_nElectron_PAS"               ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_Bar_Pt1stEle_PAS"               ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Bar_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "Jets_Bar_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
  CreateUserTH1D( "Total_Bar_hltPt1stEle_PAS"            ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Bar_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );  CreateUserTH1D( "Jets_Bar_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );
  CreateUserTH1D( "Total_Bar_Eta1stEle_PAS"	          ,    100 , -5      , 5	         );  CreateUserTH1D( "Electrons_Bar_Eta1stEle_PAS"	           ,    100 , -5      , 5          );  CreateUserTH1D( "Jets_Bar_Eta1stEle_PAS"	           ,    100 , -5      , 5          );
  CreateUserTH1D( "Total_Bar_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416     );  CreateUserTH1D( "Electrons_Bar_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416    );  CreateUserTH1D( "Jets_Bar_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416    );
  CreateUserTH1D( "Total_Bar_Charge1stEle_PAS"	          ,    2   , -1.0001 , 1.0001  );  CreateUserTH1D( "Electrons_Bar_Charge1stEle_PAS"            ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_Bar_Charge1stEle_PAS"            ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_Bar_MET_PAS"                    ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_Bar_MET_PAS"                     ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_Bar_MET_PAS"                     ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_Bar_TrkIsoHEEP7_PAS"            ,    200 , 0.0     , 100.     );  CreateUserTH1D( "Electrons_Bar_TrkIsoHEEP7_PAS"              ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Jets_Bar_TrkIsoHEEP7_PAS"              ,    200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_TrkIsoHEEP7vsPt_PAS"        ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_TrkIsoHEEP7vsHLTPt_PAS"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_Bar_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_Bar_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_Bar_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_Bar_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_Bar_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_Bar_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_Bar_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_Bar_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_Bar_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End1_nElectron_PAS"             ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End1_nElectron_PAS"              ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End1_nElectron_PAS"              ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End1_Pt1stEle_PAS"              ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End1_Pt1stEle_PAS"	           ,    1000, 0       , 1000       );  CreateUserTH1D( "Jets_End1_Pt1stEle_PAS"	           ,    1000, 0       , 1000       );
  CreateUserTH1D( "Total_End1_hltPt1stEle_PAS"           ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End1_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );  CreateUserTH1D( "Jets_End1_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );
  CreateUserTH1D( "Total_End1_Eta1stEle_PAS"	          ,    100 , -5      , 5	       );  CreateUserTH1D( "Electrons_End1_Eta1stEle_PAS"	           ,    100 , -5      , 5          );  CreateUserTH1D( "Jets_End1_Eta1stEle_PAS"	           ,    100 , -5      , 5          );
  CreateUserTH1D( "Total_End1_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416   );  CreateUserTH1D( "Electrons_End1_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416    );  CreateUserTH1D( "Jets_End1_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416    );
  CreateUserTH1D( "Total_End1_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End1_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End1_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End1_MET_PAS"                   ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End1_MET_PAS"                    ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End1_MET_PAS"                    ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End1_TrkIsoHEEP7_PAS"            ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End1_TrkIsoHEEP7_PAS"              ,    200 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End1_TrkIsoHEEP7_PAS"              ,    200 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End1_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_TrkIsoHEEP7vsHLTPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_End1_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_End1_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_End1_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_End1_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_End1_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_End1_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_End1_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_End1_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_End1_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End2_nElectron_PAS"             ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End2_nElectron_PAS"              ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End2_nElectron_PAS"              ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End2_Pt1stEle_PAS"              ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End2_Pt1stEle_PAS"	           ,    1000, 0       , 1000       );  CreateUserTH1D( "Jets_End2_Pt1stEle_PAS"	           ,    1000, 0       , 1000       );
  CreateUserTH1D( "Total_End2_hltPt1stEle_PAS"           ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End2_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );  CreateUserTH1D( "Jets_End2_hltPt1stEle_PAS"	           ,    1000, 0       , 1000     );
  CreateUserTH1D( "Total_End2_Eta1stEle_PAS"	          ,    100 , -5      , 5	       );  CreateUserTH1D( "Electrons_End2_Eta1stEle_PAS"	           ,    100 , -5      , 5          );  CreateUserTH1D( "Jets_End2_Eta1stEle_PAS"	           ,    100 , -5      , 5          );
  CreateUserTH1D( "Total_End2_Phi1stEle_PAS"	          ,    60  , -3.1416 , +3.1416   );  CreateUserTH1D( "Electrons_End2_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416    );  CreateUserTH1D( "Jets_End2_Phi1stEle_PAS"	           ,    60  , -3.1416 , +3.1416    );
  CreateUserTH1D( "Total_End2_Charge1stEle_PAS"          ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End2_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End2_Charge1stEle_PAS"           ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End2_MET_PAS"                   ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End2_MET_PAS"                    ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End2_MET_PAS"                    ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End2_TrkIsoHEEP7_PAS"            ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End2_TrkIsoHEEP7_PAS"              ,    200 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End2_TrkIsoHEEP7_PAS"              ,    200 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End2_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_TrkIsoHEEP7vsHLTPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_End2_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_End2_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_End2_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_End2_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_End2_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_End2_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_End2_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_End2_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_End2_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
  // <= 1 jet                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_Bar_lte1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_Bar_lte1Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_Bar_lte1Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_Bar_lte1Jet_Pt1stEle_PAS"          ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Bar_lte1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_Bar_lte1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_Bar_lte1Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5             );  CreateUserTH1D( "Electrons_Bar_lte1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5              );  CreateUserTH1D( "Jets_Bar_lte1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5              );
  CreateUserTH1D( "Total_Bar_lte1Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_Bar_lte1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416        );  CreateUserTH1D( "Jets_Bar_lte1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416        );
  CreateUserTH1D( "Total_Bar_lte1Jet_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001      );  CreateUserTH1D( "Electrons_Bar_lte1Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_Bar_lte1Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_Bar_lte1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_Bar_lte1Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_Bar_lte1Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_Bar_lte1Jet_TrkIsoHEEP7_PAS"       ,    200 , 0.0     , 100.     );  CreateUserTH1D( "Electrons_Bar_lte1Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Jets_Bar_lte1Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_lte1Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_lte1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_lte1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_lte1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_Bar_lte1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_Bar_lte1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_Bar_lte1Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_Bar_lte1Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_Bar_lte1Jet_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_Bar_lte1Jet_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_Bar_lte1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_Bar_lte1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_Bar_lte1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End1_lte1Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End1_lte1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End1_lte1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End1_lte1Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End1_lte1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_End1_lte1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_End1_lte1Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_End1_lte1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );  CreateUserTH1D( "Jets_End1_lte1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );
  CreateUserTH1D( "Total_End1_lte1Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_End1_lte1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End1_lte1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End1_lte1Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End1_lte1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End1_lte1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End1_lte1Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End1_lte1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End1_lte1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End1_lte1Jet_TrkIsoHEEP7_PAS"       ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End1_lte1Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End1_lte1Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End1_lte1Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_lte1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_lte1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_lte1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_End1_lte1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_End1_lte1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_End1_lte1Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_End1_lte1Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_End1_lte1Jet_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_End1_lte1Jet_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_End1_lte1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_End1_lte1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_End1_lte1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End2_lte1Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End2_lte1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End2_lte1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End2_lte1Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End2_lte1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_End2_lte1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_End2_lte1Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_End2_lte1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );  CreateUserTH1D( "Jets_End2_lte1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );
  CreateUserTH1D( "Total_End2_lte1Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_End2_lte1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End2_lte1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End2_lte1Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End2_lte1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End2_lte1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End2_lte1Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End2_lte1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End2_lte1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End2_lte1Jet_TrkIsoHEEP7_PAS"       ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End2_lte1Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End2_lte1Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End2_lte1Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_lte1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_lte1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_lte1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_End2_lte1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_End2_lte1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_End2_lte1Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_End2_lte1Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_End2_lte1Jet_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_End2_lte1Jet_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_End2_lte1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_End2_lte1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_End2_lte1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
  // 1 jet                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_Bar_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_Bar_1Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_Bar_1Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_Bar_1Jet_Pt1stEle_PAS"          ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Bar_1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_Bar_1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_Bar_1Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5             );  CreateUserTH1D( "Electrons_Bar_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5              );  CreateUserTH1D( "Jets_Bar_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5              );
  CreateUserTH1D( "Total_Bar_1Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_Bar_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416        );  CreateUserTH1D( "Jets_Bar_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416        );
  CreateUserTH1D( "Total_Bar_1Jet_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001      );  CreateUserTH1D( "Electrons_Bar_1Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_Bar_1Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_Bar_1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_Bar_1Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_Bar_1Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_Bar_1Jet_TrkIsoHEEP7_PAS"       ,    200 , 0.0     , 100.     );  CreateUserTH1D( "Electrons_Bar_1Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Jets_Bar_1Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsHLTPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_1Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_1Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_1Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_1Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_Bar_1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_Bar_1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_Bar_1Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_Bar_1Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_Bar_1Jet_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_Bar_1Jet_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_Bar_1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_Bar_1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_Bar_1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End1_1Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End1_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End1_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End1_1Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End1_1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_End1_1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_End1_1Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_End1_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );  CreateUserTH1D( "Jets_End1_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );
  CreateUserTH1D( "Total_End1_1Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_End1_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End1_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End1_1Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End1_1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End1_1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End1_1Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End1_1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End1_1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End1_1Jet_TrkIsoHEEP7_PAS"       ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End1_1Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End1_1Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsHLTPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_1Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_1Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_1Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_1Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_End1_1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_End1_1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_End1_1Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_End1_1Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_End1_1Jet_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_End1_1Jet_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_End1_1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_End1_1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_End1_1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End2_1Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End2_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End2_1Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End2_1Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End2_1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_End2_1Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_End2_1Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_End2_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );  CreateUserTH1D( "Jets_End2_1Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );
  CreateUserTH1D( "Total_End2_1Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_End2_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End2_1Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End2_1Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End2_1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End2_1Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End2_1Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End2_1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End2_1Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End2_1Jet_TrkIsoHEEP7_PAS"       ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End2_1Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End2_1Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End2_1Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_1Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_1Jet_TrkIsoHEEP7vsHLTPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_1Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_1Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_1Jet_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_1Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_1Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_End2_1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_End2_1Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_End2_1Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_End2_1Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_End2_1Jet_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_End2_1Jet_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_End2_1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_End2_1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_End2_1Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
  // 2 jet                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_Bar_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_Bar_2Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_Bar_2Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_Bar_2Jet_Pt1stEle_PAS"          ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Bar_2Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_Bar_2Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_Bar_2Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_Bar_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5              );  CreateUserTH1D( "Jets_Bar_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5              );
  CreateUserTH1D( "Total_Bar_2Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_Bar_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416        );  CreateUserTH1D( "Jets_Bar_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416        );
  CreateUserTH1D( "Total_Bar_2Jet_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001      );  CreateUserTH1D( "Electrons_Bar_2Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_Bar_2Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_Bar_2Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_Bar_2Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_Bar_2Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_Bar_2Jet_TrkIsoHEEP7_PAS"       ,    200 , 0.0     , 100.     );  CreateUserTH1D( "Electrons_Bar_2Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Jets_Bar_2Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_2Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_2Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsHLTPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_2Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_2Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_2Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_2Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_Bar_2Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_Bar_2Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_Bar_2Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_Bar_2Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_Bar_2Jet_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_Bar_2Jet_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_Bar_2Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_Bar_2Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_Bar_2Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End1_2Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End1_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End1_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End1_2Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End1_2Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_End1_2Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_End1_2Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_End1_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );  CreateUserTH1D( "Jets_End1_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );
  CreateUserTH1D( "Total_End1_2Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_End1_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End1_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End1_2Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End1_2Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End1_2Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End1_2Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End1_2Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End1_2Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End1_2Jet_TrkIsoHEEP7_PAS"       ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End1_2Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End1_2Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_2Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_2Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsHLTPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_2Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_2Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_2Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_2Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_End1_2Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_End1_2Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_End1_2Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_End1_2Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_End1_2Jet_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_End1_2Jet_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_End1_2Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_End1_2Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_End1_2Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End2_2Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End2_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End2_2Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End2_2Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End2_2Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_End2_2Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_End2_2Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_End2_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );  CreateUserTH1D( "Jets_End2_2Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );
  CreateUserTH1D( "Total_End2_2Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_End2_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End2_2Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End2_2Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End2_2Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End2_2Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End2_2Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End2_2Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End2_2Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End2_2Jet_TrkIsoHEEP7_PAS"       ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End2_2Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End2_2Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End2_2Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_2Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_2Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_2Jet_TrkIsoHEEP7vsHLTPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_2Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_2Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_2Jet_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_2Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_2Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_2Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_End2_2Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_End2_2Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_End2_2Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_End2_2Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_End2_2Jet_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_End2_2Jet_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_End2_2Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_End2_2Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_End2_2Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
  // 3 jet                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_Bar_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_Bar_3Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_Bar_3Jet_nElectron_PAS"          ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_Bar_3Jet_Pt1stEle_PAS"          ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Bar_3Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_Bar_3Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_Bar_3Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_Bar_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5              );  CreateUserTH1D( "Jets_Bar_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5              );
  CreateUserTH1D( "Total_Bar_3Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_Bar_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416        );  CreateUserTH1D( "Jets_Bar_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416        );
  CreateUserTH1D( "Total_Bar_3Jet_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001      );  CreateUserTH1D( "Electrons_Bar_3Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_Bar_3Jet_Charge1stEle_PAS"       ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_Bar_3Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_Bar_3Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_Bar_3Jet_MET_PAS"                ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_Bar_3Jet_TrkIsoHEEP7_PAS"       ,    200 , 0.0     , 100.     );  CreateUserTH1D( "Electrons_Bar_3Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Jets_Bar_3Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_3Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_3Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsHLTPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_3Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_3Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_3Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_3Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_Bar_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_Bar_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_Bar_3Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_Bar_3Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_Bar_3Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_Bar_3Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_Bar_3Jet_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_Bar_3Jet_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_Bar_3Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_Bar_3Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_Bar_3Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End1_3Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End1_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End1_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End1_3Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End1_3Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_End1_3Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_End1_3Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_End1_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );  CreateUserTH1D( "Jets_End1_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );
  CreateUserTH1D( "Total_End1_3Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_End1_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End1_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End1_3Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End1_3Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End1_3Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End1_3Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End1_3Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End1_3Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End1_3Jet_TrkIsoHEEP7_PAS"       ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End1_3Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End1_3Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_3Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_3Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsHLTPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_3Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_3Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_3Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_3Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End1_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End1_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_End1_3Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_End1_3Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_End1_3Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_End1_3Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_End1_3Jet_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_End1_3Jet_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_End1_3Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_End1_3Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_End1_3Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End2_3Jet_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End2_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End2_3Jet_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End2_3Jet_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End2_3Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_End2_3Jet_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_End2_3Jet_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_End2_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );  CreateUserTH1D( "Jets_End2_3Jet_Eta1stEle_PAS"	   ,    100 , -5      , 5            );
  CreateUserTH1D( "Total_End2_3Jet_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_End2_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End2_3Jet_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End2_3Jet_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Electrons_End2_3Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End2_3Jet_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End2_3Jet_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End2_3Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End2_3Jet_MET_PAS"               ,    100 , 0.0     , 1000.    );
  CreateUserTH1D( "Total_End2_3Jet_TrkIsoHEEP7_PAS"       ,    200 , 0.0     , 100.    );  CreateUserTH1D( "Electrons_End2_3Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );  CreateUserTH1D( "Jets_End2_3Jet_TrkIsoHEEP7_PAS"         ,    200 , 0.0     , 100.   );
  CreateUserTH2D( "Total_End2_3Jet_TrkIsoHEEP7vsPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_3Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_3Jet_TrkIsoHEEP7vsPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_3Jet_TrkIsoHEEP7vsHLTPt_PAS"            ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_3Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_3Jet_TrkIsoHEEP7vsHLTPt_PAS"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_3Jet_TrkIsoHEEP7vsHLTPt_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_3Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_3Jet_TrkIsoHEEP7vsHLTPt_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"     ,    1000,0,1000, 200 , 0.0     , 100.     );  CreateUserTH2D( "Electrons_End2_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );  CreateUserTH2D( "Jets_End2_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"              ,1000,0,1000,     200 , 0.0     , 100.    );
  CreateUserTH2D( "Total_End2_3Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Electrons_End2_3Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.);   CreateUserTH2D( "Jets_End2_3Jet_TrkIsoHEEP7vsMTenu_PAS"     , 400, 0, 2000, 200, 0.0, 100.); 
  CreateUserTH1D( "Total_End2_3Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "ElectronsHEEP_End2_3Jet_MTenu_PAS"     , 400, 0, 2000);   CreateUserTH1D( "JetsSR_End2_3Jet_MTenu_PAS"     , 400, 0, 2000); CreateUserTH1D( "JetsSB_End2_3Jet_MTenu_PAS"     , 400, 0, 2000); 
  CreateUserTH1D( "ElectronsHEEP_End2_3Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );  CreateUserTH1D( "JetsSR_End2_3Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );CreateUserTH1D( "JetsSB_End2_3Jet_Pt1stEle_PAS"	           ,    1000, 0       , 1000         );
                                                                                                                                                                                                                                                                                     
                                                                                                                                                                                                                                                                                     
  // Pile-up [0,5]                                                                                                                                                                                                                                                                   
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_Bar_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_Bar_PU0-8_nElectron_PAS"         ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_Bar_PU0-8_nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_Bar_PU0-8_Pt1stEle_PAS"         ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Bar_PU0-8_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );  CreateUserTH1D( "Jets_Bar_PU0-8_Pt1stEle_PAS"	   ,    1000, 0       , 1000           );
  CreateUserTH1D( "Total_Bar_PU0-8_Eta1stEle_PAS"	  ,    100 , -5      , 5	           );  CreateUserTH1D( "Electrons_Bar_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );  CreateUserTH1D( "Jets_Bar_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );
  CreateUserTH1D( "Total_Bar_PU0-8_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416       );  CreateUserTH1D( "Electrons_Bar_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_Bar_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_Bar_PU0-8_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001    );  CreateUserTH1D( "Electrons_Bar_PU0-8_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_Bar_PU0-8_Charge1stEle_PAS"      ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_Bar_PU0-8_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_Bar_PU0-8_MET_PAS"               ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_Bar_PU0-8_MET_PAS"               ,    100 , 0.0     , 1000.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End1_PU0-8_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End1_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End1_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End1_PU0-8_Pt1stEle_PAS"        ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End1_PU0-8_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );  CreateUserTH1D( "Jets_End1_PU0-8_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );
  CreateUserTH1D( "Total_End1_PU0-8_Eta1stEle_PAS"	  ,    100 , -5      , 5	         );  CreateUserTH1D( "Electrons_End1_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );  CreateUserTH1D( "Jets_End1_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );
  CreateUserTH1D( "Total_End1_PU0-8_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416     );  CreateUserTH1D( "Electrons_End1_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End1_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End1_PU0-8_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001    );  CreateUserTH1D( "Electrons_End1_PU0-8_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End1_PU0-8_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End1_PU0-8_MET_PAS"             ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End1_PU0-8_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End1_PU0-8_MET_PAS"              ,    100 , 0.0     , 1000.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End2_PU0-8_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End2_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End2_PU0-8_nElectron_PAS"        ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End2_PU0-8_Pt1stEle_PAS"        ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End2_PU0-8_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );  CreateUserTH1D( "Jets_End2_PU0-8_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );
  CreateUserTH1D( "Total_End2_PU0-8_Eta1stEle_PAS"	  ,    100 , -5      , 5	         );  CreateUserTH1D( "Electrons_End2_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );  CreateUserTH1D( "Jets_End2_PU0-8_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );
  CreateUserTH1D( "Total_End2_PU0-8_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416     );  CreateUserTH1D( "Electrons_End2_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_End2_PU0-8_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_End2_PU0-8_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001    );  CreateUserTH1D( "Electrons_End2_PU0-8_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End2_PU0-8_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End2_PU0-8_MET_PAS"             ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End2_PU0-8_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End2_PU0-8_MET_PAS"              ,    100 , 0.0     , 1000.    );
                                                                                                                                                                                                                                                                                     
  // Pile-up [6-10]                                                                                                                                                                                                                                                                  
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_Bar_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_Bar_PU9-UP_nElectron_PAS"        ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_Bar_PU9-UP_nElectron_PAS"        ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_Bar_PU9-UP_Pt1stEle_PAS"        ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_Bar_PU9-UP_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );  CreateUserTH1D( "Jets_Bar_PU9-UP_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );
  CreateUserTH1D( "Total_Bar_PU9-UP_Eta1stEle_PAS"	  ,    100 , -5      , 5	         );  CreateUserTH1D( "Electrons_Bar_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );  CreateUserTH1D( "Jets_Bar_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	           );
  CreateUserTH1D( "Total_Bar_PU9-UP_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416     );  CreateUserTH1D( "Electrons_Bar_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );  CreateUserTH1D( "Jets_Bar_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416      );
  CreateUserTH1D( "Total_Bar_PU9-UP_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001    );  CreateUserTH1D( "Electrons_Bar_PU9-UP_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_Bar_PU9-UP_Charge1stEle_PAS"     ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_Bar_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_Bar_PU9-UP_MET_PAS"              ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_Bar_PU9-UP_MET_PAS"              ,    100 , 0.0     , 1000.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End1_PU9-UP_nElectron_PAS"      ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End1_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End1_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End1_PU9-UP_Pt1stEle_PAS"       ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End1_PU9-UP_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );  CreateUserTH1D( "Jets_End1_PU9-UP_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );
  CreateUserTH1D( "Total_End1_PU9-UP_Eta1stEle_PAS"	  ,    100 , -5      , 5	         );  CreateUserTH1D( "Electrons_End1_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	         );  CreateUserTH1D( "Jets_End1_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	         );
  CreateUserTH1D( "Total_End1_PU9-UP_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416     );  CreateUserTH1D( "Electrons_End1_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416    );  CreateUserTH1D( "Jets_End1_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416    );
  CreateUserTH1D( "Total_End1_PU9-UP_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001  );  CreateUserTH1D( "Electrons_End1_PU9-UP_Charge1stEle_PAS"    ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End1_PU9-UP_Charge1stEle_PAS"    ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End1_PU9-UP_MET_PAS"            ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End1_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End1_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1000.    );
                                                                                                                                                                                                                                                                                     
  CreateUserTH1D( "Total_End2_PU9-UP_nElectron_PAS"      ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Electrons_End2_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );  CreateUserTH1D( "Jets_End2_PU9-UP_nElectron_PAS"       ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "Total_End2_PU9-UP_Pt1stEle_PAS"       ,    1000, 0       , 1000     );  CreateUserTH1D( "Electrons_End2_PU9-UP_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );  CreateUserTH1D( "Jets_End2_PU9-UP_Pt1stEle_PAS"	   ,    1000, 0       , 1000         );
  CreateUserTH1D( "Total_End2_PU9-UP_Eta1stEle_PAS"	  ,    100 , -5      , 5	         );  CreateUserTH1D( "Electrons_End2_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	         );  CreateUserTH1D( "Jets_End2_PU9-UP_Eta1stEle_PAS"	   ,    100 , -5      , 5	         );
  CreateUserTH1D( "Total_End2_PU9-UP_Phi1stEle_PAS"	  ,    60  , -3.1416 , +3.1416     );  CreateUserTH1D( "Electrons_End2_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416    );  CreateUserTH1D( "Jets_End2_PU9-UP_Phi1stEle_PAS"	   ,    60  , -3.1416 , +3.1416    );
  CreateUserTH1D( "Total_End2_PU9-UP_Charge1stEle_PAS"	  ,    2   , -1.0001 , 1.0001  );  CreateUserTH1D( "Electrons_End2_PU9-UP_Charge1stEle_PAS"    ,    2   , -1.0001 , 1.0001   );  CreateUserTH1D( "Jets_End2_PU9-UP_Charge1stEle_PAS"    ,    2   , -1.0001 , 1.0001   );
  CreateUserTH1D( "Total_End2_PU9-UP_MET_PAS"            ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Electrons_End2_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1000.    );  CreateUserTH1D( "Jets_End2_PU9-UP_MET_PAS"             ,    100 , 0.0     , 1000.    );

  //--------------------------------------------------------------------------
  // Tell the user how many entries we'll look at
  //--------------------------------------------------------------------------

  Long64_t nentries = GetTreeEntries();
  std::cout << "analysisClass::analysisClass(): nentries = " << nentries << std::endl;

  //--------------------------------------------------------------------------
  // Loop over the chain
  //--------------------------------------------------------------------------
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    readerTools_->LoadEntry(jentry);
    //-----------------------------------------------------------------
    // Print progress
    //-----------------------------------------------------------------
    if(jentry < 10 || jentry%10000 == 0) std::cout << "analysisClass:Loop(): jentry = " << jentry << "/" << nentries << std::endl;
    //// run ls event
    //std::cout << static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(ls) << " " << static_cast<unsigned int>(event) << std::endl;

    //--------------------------------------------------------------------------
    // Reset the cuts
    //--------------------------------------------------------------------------

    resetCuts();

    //--------------------------------------------------------------------------
    // Check good run list
    //--------------------------------------------------------------------------

    double run = readerTools_->ReadValueBranch<Double_t>("run");
    int passedJSON = passJSON ( run,
        readerTools_->ReadValueBranch<Double_t>("ls"),
        isData() ) ;

    //--------------------------------------------------------------------------
    // Find the right prescale for this event
    //--------------------------------------------------------------------------
    double min_prescale = 1;
    int passTrigger = 0;

    double LooseEle1_hltPhotonPt = readerTools_->ReadValueBranch<Double_t>("LooseEle1_hltPhotonPt");

    std::string triggerName = "";
    if ( LooseEle1_hltPhotonPt > 0.0 ) {
      if(analysisYear==2016) {
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon22")   > 0.1 && LooseEle1_hltPhotonPt >= 22.  && LooseEle1_hltPhotonPt < 30. ) { passTrigger = 1; triggerName = "Photon22"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon30")   > 0.1 && LooseEle1_hltPhotonPt >= 30.  && LooseEle1_hltPhotonPt < 36. ) { passTrigger = 1; triggerName = "Photon30"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon36")   > 0.1 && LooseEle1_hltPhotonPt >= 36.  && LooseEle1_hltPhotonPt < 50. ) { passTrigger = 1; triggerName = "Photon36"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon50")   > 0.1 && LooseEle1_hltPhotonPt >= 50.  && LooseEle1_hltPhotonPt < 75. ) { passTrigger = 1; triggerName = "Photon50"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon75")   > 0.1 && LooseEle1_hltPhotonPt >= 75.  && LooseEle1_hltPhotonPt < 90. ) { passTrigger = 1; triggerName = "Photon75"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon90")   > 0.1 && LooseEle1_hltPhotonPt >= 90.  && LooseEle1_hltPhotonPt < 120.) { passTrigger = 1; triggerName = "Photon90"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon120")  > 0.1 && LooseEle1_hltPhotonPt >= 120. && LooseEle1_hltPhotonPt < 175.) { passTrigger = 1; triggerName = "Photon120"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon175")  > 0.1 && LooseEle1_hltPhotonPt >= 175.) { passTrigger = 1; triggerName = "Photon175"; } 
      }
      else if(analysisYear==2017) {
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon25")   > 0.1 && LooseEle1_hltPhotonPt >= 25.  && LooseEle1_hltPhotonPt < 33. ) { passTrigger = 1; triggerName = "Photon25"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon33")   > 0.1 && LooseEle1_hltPhotonPt >= 33.  && LooseEle1_hltPhotonPt < 50. ) { passTrigger = 1; triggerName = "Photon33"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon50")   > 0.1 && LooseEle1_hltPhotonPt >= 50.  && LooseEle1_hltPhotonPt < 75. ) { passTrigger = 1; triggerName = "Photon50"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon75")   > 0.1 && LooseEle1_hltPhotonPt >= 75.  && LooseEle1_hltPhotonPt < 90. ) { passTrigger = 1; triggerName = "Photon75"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon90")   > 0.1 && LooseEle1_hltPhotonPt >= 90.  && LooseEle1_hltPhotonPt < 120.) { passTrigger = 1; triggerName = "Photon90"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon120")  > 0.1 && LooseEle1_hltPhotonPt >= 120. && LooseEle1_hltPhotonPt < 150.) { passTrigger = 1; triggerName = "Photon120"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon150")  > 0.1 && LooseEle1_hltPhotonPt >= 150. && LooseEle1_hltPhotonPt < 175.) { passTrigger = 1; triggerName = "Photon150"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon175")  > 0.1 && LooseEle1_hltPhotonPt >= 175. && LooseEle1_hltPhotonPt < 200.) { passTrigger = 1; triggerName = "Photon175"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon200")  > 0.1 && LooseEle1_hltPhotonPt >= 200.) { passTrigger = 1; triggerName = "Photon200"; } 
      }
      else if(analysisYear==2018) {
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon33")   > 0.1 && LooseEle1_hltPhotonPt >= 33.  && LooseEle1_hltPhotonPt < 50. ) { passTrigger = 1; triggerName = "Photon33"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon50")   > 0.1 && LooseEle1_hltPhotonPt >= 50.  && LooseEle1_hltPhotonPt < 75. ) { passTrigger = 1; triggerName = "Photon50"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon75")   > 0.1 && LooseEle1_hltPhotonPt >= 75.  && LooseEle1_hltPhotonPt < 90. ) { passTrigger = 1; triggerName = "Photon75"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon90")   > 0.1 && LooseEle1_hltPhotonPt >= 90.  && LooseEle1_hltPhotonPt < 120.) { passTrigger = 1; triggerName = "Photon90"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon120")  > 0.1 && LooseEle1_hltPhotonPt >= 120. && LooseEle1_hltPhotonPt < 150.) { passTrigger = 1; triggerName = "Photon120"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon150")  > 0.1 && LooseEle1_hltPhotonPt >= 150. && LooseEle1_hltPhotonPt < 175.) { passTrigger = 1; triggerName = "Photon150"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon175")  > 0.1 && LooseEle1_hltPhotonPt >= 175. && LooseEle1_hltPhotonPt < 200.) { passTrigger = 1; triggerName = "Photon175"; } 
        if ( readerTools_->ReadValueBranch<Double_t>("H_Photon200")  > 0.1 && LooseEle1_hltPhotonPt >= 200.) { passTrigger = 1; triggerName = "Photon200"; } 
      }
    }
    if(isData() && passTrigger) {
      //std::cout << "INFO: lookup trigger name " << triggerName << " for year: " << year << std::endl;
      min_prescale = run2PhotonTriggerPrescales.LookupPrescale(analysisYear,triggerName);
    }

    //--------------------------------------------------------------------------
    // Do pileup re-weighting
    //--------------------------------------------------------------------------

    double pileup_weight = readerTools_->ReadValueBranch<Double_t>("puWeight");
    if ( isData() ) pileup_weight = 1.0;

    //--------------------------------------------------------------------------
    // Get information about gen-level reweighting (should be for Sherpa only)
    //--------------------------------------------------------------------------
    double gen_weight = readerTools_->ReadValueBranch<Double_t>("Weight");
    if ( isData() ) gen_weight = 1.0;

    //// TopPt reweight
    //// only valid for powheg
    //std::string current_file_name ( readerTools_->GetTree()->GetCurrentFile()->GetName());
    //if(current_file_name.find("TT_") != std::string::npos) {
    //  gen_weight*=TopPtWeight;
    //}

    // Electron scale factors for MC only
    if(!isData()) {
      bool verbose = false;
      float ele1ECorr = readerTools_->ReadValueBranch<Double_t>("LooseEle1_ECorr");
      float ele1PtUncorr = ele1ECorr != 0 ? readerTools_->ReadValueBranch<Double_t>("LooseEle1_Pt")/ele1ECorr : readerTools_->ReadValueBranch<Double_t>("LooseEle1_Pt");
      float recoSFEle1 = recoScaleFactorReader->LookupValue(readerTools_->ReadValueBranch<Double_t>("LooseEle1_SCEta"),ele1PtUncorr,verbose);
      //float heepSFEle1 = ElectronScaleFactors2016::LookupHeepSF(readerTools_->ReadValueBranch<Double_t>("LooseEle1_SCEta"));
      //float totalScaleFactor = recoSFEle1*heepSFEle1;
      //gen_weight*=totalScaleFactor;
      gen_weight*=recoSFEle1;
    }
    // add these to pileup weight
    pileup_weight*=gen_weight;

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
    fillVariableWithValue("PassGlobalSuperTightHalo2016Filter" , int(readerTools_->ReadValueBranch<Double_t>("PassGlobalSuperTightHalo2016Filter")     == 1), min_prescale * pileup_weight);
    fillVariableWithValue("PassGoodVertices"                   , int(readerTools_->ReadValueBranch<Double_t>("PassGoodVertices")                       == 1), min_prescale * pileup_weight);
    fillVariableWithValue("PassHBHENoiseFilter"                , int(readerTools_->ReadValueBranch<Double_t>("PassHBHENoiseFilter")                    == 1), min_prescale * pileup_weight);
    fillVariableWithValue("PassHBHENoiseIsoFilter"             , int(readerTools_->ReadValueBranch<Double_t>("PassHBHENoiseIsoFilter")                 == 1), min_prescale * pileup_weight);
    // eBadScFilter not suggested for MC
    if(isData())
      fillVariableWithValue("PassBadEESupercrystalFilter"      , int(readerTools_->ReadValueBranch<Double_t>("PassBadEESupercrystalFilter")            == 1), min_prescale * pileup_weight);
    else
      fillVariableWithValue("PassBadEESupercrystalFilter"      , 1                                                                                          , min_prescale * pileup_weight);
    fillVariableWithValue("PassEcalDeadCellTrigPrim"           , int(readerTools_->ReadValueBranch<Double_t>("PassEcalDeadCellTrigPrim")               == 1), min_prescale * pileup_weight);
    // not recommended
    //fillVariableWithValue("PassChargedCandidateFilter"         , int(readerTools_->ReadValueBranch<Double_t>("PassChargedCandidateFilter")             == 1), min_prescale * pileup_weight);
    fillVariableWithValue("PassBadPFMuonFilter"                , int(readerTools_->ReadValueBranch<Double_t>("PassBadPFMuonFilter")                    == 1), min_prescale * pileup_weight);
    // EcalBadCalibV2 for 2017, 2018
    if(analysisYear > 2016)
      fillVariableWithValue("PassEcalBadCalibV2Filter"         , int(readerTools_->ReadValueBranch<Double_t>("PassEcalBadCalibV2Filter")               == 1), min_prescale * pileup_weight);
    else
      fillVariableWithValue("PassEcalBadCalibV2Filter"         , 1                                                                                          , min_prescale * pileup_weight);

    // Muon variables ( for veto ) 					      	       
    // remove muon veto
    //fillVariableWithValue(   "nMuon"                   , readerTools_->ReadValueBranch<Double_t>("nMuon_ptCut")     , min_prescale * pileup_weight );
    double nLooseEle_ptCut = readerTools_->ReadValueBranch<Double_t>("nLooseEle_ptCut");
    double nVLooseEle_ptCut = readerTools_->ReadValueBranch<Double_t>("nVLooseEle_ptCut");
    double nJetLooseEle_ptCut = readerTools_->ReadValueBranch<Double_t>("nJetLooseEle_ptCut");
    // 1st Electron variables				      		              
    fillVariableWithValue(   "nEle"                    , nVLooseEle_ptCut , min_prescale * pileup_weight );

    // use uncorrected energy
    double LooseEle1_Pt = readerTools_->ReadValueBranch<Double_t>("LooseEle1_Pt");// loose ele Pt is now uncorrected /readerTools_->ReadValueBranch<Double_t>("LooseEle1_ECorr");
    double LooseEle2_Pt = readerTools_->ReadValueBranch<Double_t>("LooseEle2_Pt");// loose ele Pt is now uncorrected /readerTools_->ReadValueBranch<Double_t>("LooseEle2_ECorr");
    double LooseEle3_Pt = readerTools_->ReadValueBranch<Double_t>("LooseEle3_Pt");// loose ele Pt is now uncorrected /readerTools_->ReadValueBranch<Double_t>("LooseEle3_ECorr");
    //fillVariableWithValue(   "Pt1stEle"               , LooseEle1_SCEnergy/cosh(LooseEle1_SCEta)Heep            , min_prescale * pileup_weight );
    fillVariableWithValue(   "Pt1stEle"                , LooseEle1_Pt                , min_prescale * pileup_weight );
    fillVariableWithValue(   "SCEta1stEle"             , readerTools_->ReadValueBranch<Double_t>("LooseEle1_SCEta")             , min_prescale * pileup_weight );
    fillVariableWithValue(   "Phi1stEle"               , readerTools_->ReadValueBranch<Double_t>("LooseEle1_Phi")               , min_prescale * pileup_weight );
    fillVariableWithValue(   "HLTPt1stEle"             , LooseEle1_hltPhotonPt       , min_prescale * pileup_weight );
    if(analysisYear==2016)
      fillVariableWithValue(   "H_Photon22"              , readerTools_->ReadValueBranch<Double_t>("H_Photon22")            , min_prescale * pileup_weight );
    else
      fillVariableWithValue(   "H_Photon25"              , readerTools_->ReadValueBranch<Double_t>("H_Photon25")            , min_prescale * pileup_weight );
    if(analysisYear==2016)
      fillVariableWithValue(   "H_Photon30"              , readerTools_->ReadValueBranch<Double_t>("H_Photon30")            , min_prescale * pileup_weight );
    else
      fillVariableWithValue(   "H_Photon33"              , readerTools_->ReadValueBranch<Double_t>("H_Photon33")            , min_prescale * pileup_weight );
    if(analysisYear==2016)
      fillVariableWithValue(   "H_Photon36"              , readerTools_->ReadValueBranch<Double_t>("H_Photon36")            , min_prescale * pileup_weight );
    fillVariableWithValue(   "H_Photon50"              , readerTools_->ReadValueBranch<Double_t>("H_Photon50")            , min_prescale * pileup_weight );   
    fillVariableWithValue(   "H_Photon75"              , readerTools_->ReadValueBranch<Double_t>("H_Photon75")            , min_prescale * pileup_weight );   
    fillVariableWithValue(   "H_Photon90"              , readerTools_->ReadValueBranch<Double_t>("H_Photon90")            , min_prescale * pileup_weight );   
    fillVariableWithValue(   "H_Photon120"             , readerTools_->ReadValueBranch<Double_t>("H_Photon120")           , min_prescale * pileup_weight );   
    if(analysisYear>2016)
      fillVariableWithValue(   "H_Photon150"             , readerTools_->ReadValueBranch<Double_t>("H_Photon150")           , min_prescale * pileup_weight );   
    fillVariableWithValue(   "H_Photon175"             , readerTools_->ReadValueBranch<Double_t>("H_Photon175")           , min_prescale * pileup_weight );    
    if(analysisYear>2016)
      fillVariableWithValue(   "H_Photon200"             , readerTools_->ReadValueBranch<Double_t>("H_Photon200")           , min_prescale * pileup_weight );   

    // 1st JET variables                                     		                    
    fillVariableWithValue(   "nJet"                          , nJetLooseEle_ptCut          , min_prescale * pileup_weight );

    double PFMET_Type1_Pt  = readerTools_->ReadValueBranch<Double_t>("PFMET_Type1_Pt");
    double PFMET_Type1_Phi  = readerTools_->ReadValueBranch<Double_t>("PFMET_Type1_Phi");
    // MET									            
    fillVariableWithValue(   "MET"                           , PFMET_Type1_Pt           , min_prescale * pileup_weight );

    double LooseEle1_SCEta = readerTools_->ReadValueBranch<Double_t>("LooseEle1_SCEta");
    double LooseEle2_SCEta = readerTools_->ReadValueBranch<Double_t>("LooseEle2_SCEta");
    double LooseEle1_Eta = readerTools_->ReadValueBranch<Double_t>("LooseEle1_Eta");
    double LooseEle2_Eta = readerTools_->ReadValueBranch<Double_t>("LooseEle2_Eta");
    double LooseEle3_Eta = readerTools_->ReadValueBranch<Double_t>("LooseEle3_Eta");
    double LooseEle1_Phi = readerTools_->ReadValueBranch<Double_t>("LooseEle1_Phi");
    double LooseEle2_Phi = readerTools_->ReadValueBranch<Double_t>("LooseEle2_Phi");
    double LooseEle3_Phi = readerTools_->ReadValueBranch<Double_t>("LooseEle3_Phi");
    double JetLooseEle1_Eta = readerTools_->ReadValueBranch<Double_t>("JetLooseEle1_Eta");
    double JetLooseEle2_Eta = readerTools_->ReadValueBranch<Double_t>("JetLooseEle2_Eta");
    double JetLooseEle3_Eta = readerTools_->ReadValueBranch<Double_t>("JetLooseEle3_Eta");
    double JetLooseEle4_Eta = readerTools_->ReadValueBranch<Double_t>("JetLooseEle4_Eta");
    double JetLooseEle5_Eta = readerTools_->ReadValueBranch<Double_t>("JetLooseEle5_Eta");
    double JetLooseEle1_Phi = readerTools_->ReadValueBranch<Double_t>("JetLooseEle1_Phi");
    double JetLooseEle2_Phi = readerTools_->ReadValueBranch<Double_t>("JetLooseEle2_Phi");
    double JetLooseEle3_Phi = readerTools_->ReadValueBranch<Double_t>("JetLooseEle3_Phi");
    double JetLooseEle4_Phi = readerTools_->ReadValueBranch<Double_t>("JetLooseEle4_Phi");
    double JetLooseEle5_Phi = readerTools_->ReadValueBranch<Double_t>("JetLooseEle5_Phi");
    double JetLooseEle1_Pt = readerTools_->ReadValueBranch<Double_t>("JetLooseEle1_Pt");
    double JetLooseEle2_Pt = readerTools_->ReadValueBranch<Double_t>("JetLooseEle2_Pt");
    double JetLooseEle3_Pt = readerTools_->ReadValueBranch<Double_t>("JetLooseEle3_Pt");
    double JetLooseEle4_Pt = readerTools_->ReadValueBranch<Double_t>("JetLooseEle4_Pt");
    double JetLooseEle5_Pt = readerTools_->ReadValueBranch<Double_t>("JetLooseEle5_Pt");
    double DR_Ele1Jet1 = readerTools_->ReadValueBranch<Double_t>("DR_Ele1Jet1");
    double DR_Ele1Jet2 = readerTools_->ReadValueBranch<Double_t>("DR_Ele1Jet2");
    // max deltaR(ele,jets)
    double min_DR_EleJet = 999.0;
    double DR_Ele1Jet3 = 999.0;
    double DR_Ele1Jet4 = 999.0;
    double DR_Ele1Jet5 = 999.0;
    if ( nJetLooseEle_ptCut > 2 ) {
      TLorentzVector ele1,  jet3;
      ele1.SetPtEtaPhiM ( LooseEle1_Pt, LooseEle1_Eta, LooseEle1_Phi, 0.0 );
      jet3.SetPtEtaPhiM ( JetLooseEle3_Pt, JetLooseEle3_Eta, JetLooseEle3_Phi, 0.0 );	 
      DR_Ele1Jet3 = ele1.DeltaR ( jet3 ) ;
      if ( nJetLooseEle_ptCut > 3 ) {
        TLorentzVector ele1,  jet4;
        ele1.SetPtEtaPhiM ( LooseEle1_Pt, LooseEle1_Eta, LooseEle1_Phi, 0.0 );
        jet4.SetPtEtaPhiM ( JetLooseEle4_Pt, JetLooseEle4_Eta, JetLooseEle4_Phi, 0.0 );	 
        DR_Ele1Jet4 = ele1.DeltaR ( jet4 ) ;
        if ( nJetLooseEle_ptCut > 4 ) {
          TLorentzVector ele1,  jet5;
          ele1.SetPtEtaPhiM ( LooseEle1_Pt, LooseEle1_Eta, LooseEle1_Phi, 0.0 );
          jet5.SetPtEtaPhiM ( JetLooseEle5_Pt, JetLooseEle5_Eta, JetLooseEle5_Phi, 0.0 );	 
          DR_Ele1Jet5 = ele1.DeltaR ( jet5 ) ;
        }
      }
    }
    if ( nJetLooseEle_ptCut > 0 ) {
      if ( DR_Ele1Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet1;
      if ( nJetLooseEle_ptCut > 1 ) {
        if ( DR_Ele1Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet2;
        if ( nJetLooseEle_ptCut > 2 ) {
          if ( DR_Ele1Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet3;
          if ( nJetLooseEle_ptCut > 3 ) {
            if ( DR_Ele1Jet4 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet4;
            if ( nJetLooseEle_ptCut > 4 ) {
              if ( DR_Ele1Jet5 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet5;
            }
          }
        }
      }
    }
    fillVariableWithValue(   "minDREleJets"                          , min_DR_EleJet          , min_prescale * pileup_weight );

    double mDPhi_METEle1 = readerTools_->ReadValueBranch<Double_t>("mDPhi_METEle1");
    double nJetLooseEle_store = readerTools_->ReadValueBranch<Double_t>("nJetLooseEle_store");
    double nLooseEle_store = readerTools_->ReadValueBranch<Double_t>("nLooseEle_store");
    TLorentzVector loose_ele1, met;
    // need to use uncorrected Pt
    loose_ele1.SetPtEtaPhiM ( LooseEle1_Pt , LooseEle1_Eta , LooseEle1_Phi , 0.0 );
    met.SetPtEtaPhiM        ( PFMET_Type1_Pt         , 0.0             , PFMET_Type1_Phi         , 0.0 );
    mDPhi_METEle1= fabs(loose_ele1.DeltaPhi ( met ));
    fillVariableWithValue(   "mDeltaPhiMETEle"          , mDPhi_METEle1           , min_prescale * pileup_weight );

    //XXX FIXME: these got replaced in the tree previously
    //sT_enujj = LooseEle1_Pt + PFMET_Type1_Pt + JetLooseEle1_Pt + JetLooseEle2_Pt ;
    //MT_Ele1MET = sqrt(2 * LooseEle1_Pt * PFMET_Type1_Pt  * (1 - cos(loose_ele1.DeltaPhi ( met) ) ) );

    double MT_Jet1MET, MT_Jet2MET, MT_Ele1Jet1, MT_Ele1Jet2;//, MT_Ele1MET_Type01;
    //double mDPhi_METType01_Ele1, mDPhi_METType01_Jet1, mDPhi_METType01_Jet2;

    // 1st JET variables                                     		           
    if ( nJetLooseEle_store > 0 ) { 						           
      //fillVariableWithValue( "mDeltaPhiMET1stJet"       , mDPhi_METJet1           , min_prescale * pileup_weight );

      TVector2 v_MET;
      TVector2 v_jet;
      //TVector2 v_MET_Type01;
      //v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
      v_MET.SetMagPhi( PFMET_Type1_Pt , PFMET_Type1_Phi  );
      v_jet.SetMagPhi( JetLooseEle1_Pt, JetLooseEle1_Phi );
      //mDPhi_METType01_Jet1 = fabs(v_MET_Type01.DeltaPhi ( v_jet ));
      float deltaphi = v_MET.DeltaPhi(v_jet);
      MT_Jet1MET = sqrt ( 2 * JetLooseEle1_Pt * PFMET_Type1_Pt * ( 1 - cos ( deltaphi ) ) );
    }									           

    // 2nd JET variables                                     		           
    if ( nJetLooseEle_store > 1 ) { 	                                      	           
      //fillVariableWithValue( "ST"                       , sT_enujj                , min_prescale * pileup_weight );

      TVector2 v_MET;
      TVector2 v_jet;
      //TVector2 v_MET_Type01;
      //v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
      v_MET.SetMagPhi( PFMET_Type1_Pt , PFMET_Type1_Phi  );
      v_jet.SetMagPhi( JetLooseEle2_Pt, JetLooseEle2_Phi );
      //mDPhi_METType01_Jet2 = fabs(v_MET_Type01.DeltaPhi ( v_jet ));
      float deltaphi = v_MET.DeltaPhi(v_jet);
      MT_Jet2MET = sqrt ( 2 * JetLooseEle2_Pt * PFMET_Type1_Pt * ( 1 - cos ( deltaphi ) ) );
    }

    // 3rd JET variables 
    // if ( nJetLooseEle_store > 2 ) {
    //   fillVariableWithValue( "Jet3_Pt"                  , JetLooseEle3_Pt         , min_prescale * pileup_weight );
    //   fillVariableWithValue( "Jet3_Eta"                 , JetLooseEle3_Eta        , min_prescale * pileup_weight );
    // }

    //// 1 electron, 1 jet variables 
    //if ( nLooseEle_store > 0 && nJetLooseEle_store > 0 ) { 
    //  fillVariableWithValue ( "DR_Ele1Jet1"             , DR_Ele1Jet1             , min_prescale * pileup_weight );


    //  TVector2 v_ele;
    //  TVector2 v_jet1;
    //  //TVector2 v_MET_Type01;
    //  //v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
    //  //v_ele .SetMagPhi ( LooseEle1_Pt, LooseEle1_Phi );
    //  // need to use uncorrected Pt
    //  v_ele .SetMagPhi ( LooseEle1_SCEnergy/cosh(LooseEle1_SCEta), LooseEle1_Phi );
    //  v_jet1.SetMagPhi ( JetLooseEle1_Pt, JetLooseEle1_Phi );
    //  float deltaphi = v_ele.DeltaPhi ( v_jet1 );
    //  //MT_Ele1Jet2 = sqrt ( 2 * JetLooseEle1_Pt * LooseEle1_Pt * ( 1 - cos ( deltaphi ) ) );
    //  MT_Ele1Jet2 = sqrt ( 2 * JetLooseEle1_Pt * LooseEle1_SCEnergy/cosh(LooseEle1_SCEta) * ( 1 - cos ( deltaphi ) ) );

    //}

    //// 1 electron, 2 jet variables 
    //if ( nLooseEle_store > 0 && nJetLooseEle_store > 1 ) { 
    //  fillVariableWithValue ( "DR_Ele1Jet2"             , DR_Ele1Jet2             , min_prescale * pileup_weight );

    //  TVector2 v_ele;
    //  TVector2 v_jet2;
    //  //TVector2 v_MET_Type01;
    //  //v_MET_Type01.SetMagPhi( PFMET_Type01_Pt , PFMET_Type01_Phi  );
    //  //v_ele .SetMagPhi ( LooseEle1_Pt, LooseEle1_Phi );
    //  // need to use uncorrected Pt
    //  v_ele .SetMagPhi ( LooseEle1_SCEnergy/cosh(LooseEle1_SCEta), LooseEle1_Phi );
    //  v_jet2.SetMagPhi ( JetLooseEle2_Pt, JetLooseEle2_Phi );
    //  float deltaphi = v_ele.DeltaPhi ( v_jet2 );
    //  //MT_Ele1Jet2 = sqrt ( 2 * JetLooseEle2_Pt * LooseEle1_Pt * ( 1 - cos ( deltaphi ) ) );
    //  MT_Ele1Jet2 = sqrt ( 2 * JetLooseEle2_Pt * LooseEle1_SCEnergy/cosh(LooseEle1_SCEta) * ( 1 - cos ( deltaphi ) ) );

    //}

    double MT_JetMET;
    double Mej;
    bool mejSelectedJet1 = true;
    double M_e1j1 = readerTools_->ReadValueBranch<Double_t>("M_e1j1");
    double M_e1j2 = readerTools_->ReadValueBranch<Double_t>("M_e1j2");

    if ( fabs ( MT_Jet1MET - MT_Ele1Jet2 ) < fabs( MT_Jet2MET - MT_Ele1Jet1 )){
      MT_JetMET = MT_Jet1MET;
      Mej = M_e1j2;
      mejSelectedJet1 = false;
    } else { 
      MT_JetMET = MT_Jet2MET;
      Mej = M_e1j1;
    }	 

    double M_e1e2 = 0.0;
    if ( nLooseEle_store > 1) {
      TLorentzVector v_ele1, v_ele2, v_sum;
      v_ele1.SetPtEtaPhiM ( LooseEle1_Pt, LooseEle1_Eta, LooseEle1_Phi, 0.0 );
      v_ele2.SetPtEtaPhiM ( LooseEle2_Pt, LooseEle2_Eta, LooseEle2_Phi, 0.0 );
      v_sum = v_ele1+v_ele2;
      M_e1e2 = v_sum.M();
    }
    double M_e1e3 = 0.0;
    double M_e2e3 = 0.0;
    if ( nLooseEle_store > 2) {
      TLorentzVector v_ele1, v_ele2, v_ele3, v_sum;
      v_ele1.SetPtEtaPhiM ( LooseEle1_Pt, LooseEle1_Eta, LooseEle1_Phi, 0.0 );
      v_ele2.SetPtEtaPhiM ( LooseEle2_Pt, LooseEle2_Eta, LooseEle2_Phi, 0.0 );
      v_ele3.SetPtEtaPhiM ( LooseEle3_Pt, LooseEle3_Eta, LooseEle3_Phi, 0.0 );
      v_sum = v_ele1+v_ele3;
      M_e1e3 = v_sum.M();
      v_sum = v_ele2+v_ele3;
      M_e2e3 = v_sum.M();
    }
     

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

    //if(run==284043)
    //{
    //  std::cout << (uint64_t)run << " " << (uint64_t)ls << " " << (uint64_t)event << std::endl;
    //  evaluateCuts(true);
    //}
    //else
    //  evaluateCuts(false);
    evaluateCuts();

    //--------------------------------------------------------------------------
    // Fill preselection plots
    //--------------------------------------------------------------------------

    bool passed_denominator = passedAllPreviousCuts("denominator");

    bool isBarrel  = false;
    bool isEndcap1 = false;
    bool isEndcap2 = false;

    //if( fabs( LooseEle1_Eta  ) < eleEta_bar )        isBarrel = true;
    //if( fabs( LooseEle1_Eta  ) > eleEta_end1_min &&
    //    fabs( LooseEle1_Eta  ) < eleEta_end1_max )   isEndcap1 = true;
    //if( fabs( LooseEle1_Eta  ) > eleEta_end2_min &&
    //    fabs( LooseEle1_Eta  ) < eleEta_end2_max )   isEndcap2 = true;

    if( fabs( LooseEle1_SCEta  ) < eleEta_bar )        isBarrel = true;
    if( fabs( LooseEle1_SCEta  ) > eleEta_end1_min &&
        fabs( LooseEle1_SCEta  ) < eleEta_end1_max )   isEndcap1 = true;
    if( fabs( LooseEle1_SCEta  ) > eleEta_end2_min &&
        fabs( LooseEle1_SCEta  ) < eleEta_end2_max )   isEndcap2 = true;

    if ( passed_denominator ) { 
      double nVertex = readerTools_->ReadValueBranch<Double_t>("nVertex");
      double LooseEle1_Charge = readerTools_->ReadValueBranch<Double_t>("LooseEle1_Charge");
      //double LooseEle1_EcalDriven = readerTools_->ReadValueBranch<Double_t>("LooseEle1_EcalDriven");
      double LooseEle1_TrkIsoHEEP7 = readerTools_->ReadValueBranch<Double_t>("LooseEle1_TrkIsoHEEP7");
      double MT_Ele1MET = readerTools_->ReadValueBranch<Double_t>("MT_Ele1MET");

      // debugging 
      // fillReducedSkimTree();

      // run ls event
      //std::cout << "[QCD FR] passing run/ls/event: " << static_cast<int>(run) << " " << static_cast<int>(ls) << " " << ((unsigned int)event) << std::endl;

      FillUserTH1D("Total_minPrescale"                     ,   min_prescale );
      FillUserTH1D("Total_pileupWeight"                    ,   pileup_weight);

      FillUserTH1D( "Total_nVertex_PAS"           , nVertex                      , min_prescale * pileup_weight);
      FillUserTH1D( "Total_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
      FillUserTH1D( "Total_Pt1stEle_PAS"	         , LooseEle1_Pt                 , min_prescale * pileup_weight);
      FillUserTH1D( "Total_hltPt1stEle_PAS"       , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
      FillUserTH1D( "Total_Eta1stEle_PAS"	 , LooseEle1_Eta                , min_prescale * pileup_weight);
      FillUserTH1D( "Total_Phi1stEle_PAS"	 , LooseEle1_Phi                , min_prescale * pileup_weight);
      FillUserTH1D( "Total_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
      FillUserTH1D( "Total_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
      FillUserTH1D( "Total_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
      FillUserTH2D( "Total_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
      FillUserTH2D( "Total_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
      if(analysisYear==2018) {
        if(run >= 319077 || !isData()) {
          FillUserTH2D( "Total_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
            FillUserTH2D( "Total_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          else
            FillUserTH2D( "Total_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }
        else
          FillUserTH2D( "Total_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
      }
      FillUserTH2D( "Total_TrkIsoHEEP7vsMTenu_PAS"   , MT_Ele1MET, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
      FillUserTH1D( "Total_MTenu_PAS", MT_Ele1MET, min_prescale * pileup_weight);

      if ( isBarrel ) {

        FillUserTH1D( "Total_Bar_nVertex_PAS"           , nVertex                      , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_Pt1stEle_PAS"	         , LooseEle1_Pt                 , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_hltPt1stEle_PAS"       , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_Eta1stEle_PAS"	 , LooseEle1_Eta                , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_Phi1stEle_PAS"	 , LooseEle1_Phi                , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
        FillUserTH1D( "Total_Bar_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        FillUserTH2D( "Total_Bar_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        FillUserTH2D( "Total_Bar_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        if(analysisYear==2018) {
          if(run >= 319077 || !isData()) {
            FillUserTH2D( "Total_Bar_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
              FillUserTH2D( "Total_Bar_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            else
              FillUserTH2D( "Total_Bar_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          else
            FillUserTH2D( "Total_Bar_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }
        FillUserTH2D( "Total_Bar_TrkIsoHEEP7vsMTenu_PAS"   , MT_Ele1MET, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        FillUserTH1D( "Total_Bar_MTenu_PAS", MT_Ele1MET, min_prescale * pileup_weight);

        if ( nJetLooseEle_ptCut <= 1 ) {
          FillUserTH1D( "Total_Bar_lte1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_lte1Jet_Pt1stEle_PAS"	      , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_lte1Jet_Eta1stEle_PAS"	      , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_lte1Jet_Phi1stEle_PAS"	      , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_lte1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_Bar_lte1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_lte1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_Bar_lte1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Total_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Total_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Total_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Total_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          FillUserTH2D( "Total_Bar_lte1Jet_TrkIsoHEEP7vsMTenu_PAS"   , MT_Ele1MET, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_lte1Jet_MTenu_PAS", MT_Ele1MET, min_prescale * pileup_weight);
        }

        if ( nJetLooseEle_ptCut >= 1 ) {
          FillUserTH1D( "Total_Bar_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_1Jet_Pt1stEle_PAS"	      , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_1Jet_Eta1stEle_PAS"	      , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_1Jet_Phi1stEle_PAS"	      , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_Bar_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          FillUserTH2D( "Total_Bar_1Jet_TrkIsoHEEP7vsMTenu_PAS"   , MT_Ele1MET, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_1Jet_MTenu_PAS", MT_Ele1MET, min_prescale * pileup_weight);
        }

        if ( nJetLooseEle_ptCut >= 2 ) {
          FillUserTH1D( "Total_Bar_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_2Jet_Pt1stEle_PAS"	      , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_2Jet_Eta1stEle_PAS"	      , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_2Jet_Phi1stEle_PAS"	      , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_Bar_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          FillUserTH2D( "Total_Bar_2Jet_TrkIsoHEEP7vsMTenu_PAS"   , MT_Ele1MET, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_2Jet_MTenu_PAS", MT_Ele1MET, min_prescale * pileup_weight);
        }

        if ( nJetLooseEle_ptCut >= 3 ) {
          FillUserTH1D( "Total_Bar_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_3Jet_Pt1stEle_PAS"	      , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_3Jet_Eta1stEle_PAS"	      , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_3Jet_Phi1stEle_PAS"	      , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_Bar_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          FillUserTH2D( "Total_Bar_3Jet_TrkIsoHEEP7vsMTenu_PAS"   , MT_Ele1MET, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_3Jet_MTenu_PAS", MT_Ele1MET, min_prescale * pileup_weight);
        }

        if ( nVertex >= 0 && nVertex <= 5 ){
          FillUserTH1D( "Total_Bar_PU0-8_nElectron_PAS"          , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU0-8_Pt1stEle_PAS"           , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU0-8_Eta1stEle_PAS"          , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU0-8_Phi1stEle_PAS"          , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU0-8_Charge1stEle_PAS"       , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_Bar_PU0-8_MET_PAS"                , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        }						          

        if ( nVertex >= 6 && nVertex <= 10 ){      
          FillUserTH1D( "Total_Bar_PU9-UP_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU9-UP_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU9-UP_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU9-UP_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_Bar_PU9-UP_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
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
        FillUserTH1D( "Total_End1_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End1_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        FillUserTH2D( "Total_End1_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        if(analysisYear==2018) {
          if(run >= 319077 || !isData()) {
            FillUserTH2D( "Total_End1_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
              FillUserTH2D( "Total_End1_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            else
              FillUserTH2D( "Total_End1_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          else
            FillUserTH2D( "Total_End1_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }
        FillUserTH2D( "Total_End1_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);

        if ( nJetLooseEle_ptCut <= 1 ) {
          FillUserTH1D( "Total_End1_lte1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_lte1Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_lte1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_lte1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_lte1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End1_lte1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_lte1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End1_lte1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Total_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Total_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Total_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Total_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          FillUserTH2D( "Total_End1_lte1Jet_TrkIsoHEEP7vsMTenu_PAS"   , MT_Ele1MET, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_lte1Jet_MTenu_PAS", MT_Ele1MET, min_prescale * pileup_weight);
        }

        if ( nJetLooseEle_ptCut >= 1 ) {
          FillUserTH1D( "Total_End1_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_1Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End1_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          FillUserTH2D( "Total_End1_1Jet_TrkIsoHEEP7vsMTenu_PAS"   , MT_Ele1MET, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_1Jet_MTenu_PAS", MT_Ele1MET, min_prescale * pileup_weight);
        }

        if ( nJetLooseEle_ptCut >= 2 ) {
          FillUserTH1D( "Total_End1_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_2Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End1_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          FillUserTH2D( "Total_End1_2Jet_TrkIsoHEEP7vsMTenu_PAS"   , MT_Ele1MET, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_2Jet_MTenu_PAS", MT_Ele1MET, min_prescale * pileup_weight);
        }

        if ( nJetLooseEle_ptCut >= 3 ) {
          FillUserTH1D( "Total_End1_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_3Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End1_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          FillUserTH2D( "Total_End1_3Jet_TrkIsoHEEP7vsMTenu_PAS"   , MT_Ele1MET, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_3Jet_MTenu_PAS", MT_Ele1MET, min_prescale * pileup_weight);
        }

        if ( nVertex >= 0 && nVertex <= 5 ){
          FillUserTH1D( "Total_End1_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End1_PU0-8_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        }						          

        if ( nVertex >= 6 && nVertex <= 10 ){      
          FillUserTH1D( "Total_End1_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End1_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
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
        FillUserTH1D( "Total_End2_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        FillUserTH1D( "Total_End2_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        FillUserTH2D( "Total_End2_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        if(analysisYear==2018) {
          if(run >= 319077 || !isData()) {
            FillUserTH2D( "Total_End2_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
              FillUserTH2D( "Total_End2_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            else
              FillUserTH2D( "Total_End2_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          else
            FillUserTH2D( "Total_End2_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }
        FillUserTH2D( "Total_End2_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);

        if ( nJetLooseEle_ptCut <= 1 ) {
          FillUserTH1D( "Total_End2_lte1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_lte1Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_lte1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_lte1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_lte1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End2_lte1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_lte1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End2_lte1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Total_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Total_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Total_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Total_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          FillUserTH2D( "Total_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }

        if ( nJetLooseEle_ptCut >= 1 ) {
          FillUserTH1D( "Total_End2_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_1Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End2_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End2_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Total_End2_1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Total_End2_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Total_End2_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Total_End2_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          FillUserTH2D( "Total_End2_1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }

        if ( nJetLooseEle_ptCut >= 2 ) {
          FillUserTH1D( "Total_End2_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_2Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End2_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End2_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Total_End2_2Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Total_End2_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Total_End2_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Total_End2_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          FillUserTH2D( "Total_End2_2Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }

        if ( nJetLooseEle_ptCut >= 3 ) {
          FillUserTH1D( "Total_End2_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_3Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End2_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Total_End2_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Total_End2_3Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Total_End2_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Total_End2_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Total_End2_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          FillUserTH2D( "Total_End2_3Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }

        if ( nVertex >= 0 && nVertex <= 5 ){
          FillUserTH1D( "Total_End2_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End2_PU0-8_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        }						          

        if ( nVertex >= 6 && nVertex <= 10 ){      
          FillUserTH1D( "Total_End2_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Total_End2_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
          FillUserTH1D( "Total_End2_PU9-UP_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        }
      }

      //----------------------------------------------------------------------
      // HEEP' 7.0 ID part
      //  same as HEEPv7.0, but no TrkIso cut
      //----------------------------------------------------------------------
      bool passHEEPprime = (
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPMinPtCut") == 1                             &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleSCEtaMultiRangeCut") == 1             &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleDEtaInSeedCut") == 1                  &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleDPhiInCut") == 1                      &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleFull5x5SigmaIEtaIEtaWithSatCut") == 1 &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleFull5x5E2x5OverE5x5WithSatCut") == 1  &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleHadronicOverEMLinearCut") == 1        &&
          //readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleTrkPtIsoCut") == 1                    &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleEmHadD1IsoRhoCut") == 1               &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleDxyCut") == 1                         &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleMissingHitsCut") == 1                 &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPEcalDrivenCut") == 1
          );

      // full HEEP has TrkIso included
      bool passHEEP = readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPID") == 1;

      // for jets
      // they must pass all except deltaEtaSeed, deltaPhi, sigmaIEtaIEta, shape
      bool passJet = (
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPMinPtCut") == 1                             &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleSCEtaMultiRangeCut") == 1             &&
          //readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleDEtaInSeedCut") == 1                  &&
          //readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleDPhiInCut") == 1                      &&
          //readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleFull5x5SigmaIEtaIEtaWithSatCut") == 1 &&
          //readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleFull5x5E2x5OverE5x5WithSatCut") == 1  &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleHadronicOverEMLinearCut") == 1        &&
          //readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleTrkPtIsoCut") == 1                    &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleEmHadD1IsoRhoCut") == 1               &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleDxyCut") == 1                         &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPGsfEleMissingHitsCut") == 1                 &&
          readerTools_->ReadValueBranch<Double_t>("LooseEle1_PassHEEPEcalDrivenCut") == 1
          );

      ////XXX SIC TEST
      //if(run==284043 && (ls==116 || ls==146 || ls==172 || ls==199 || ls==45 || ls==8) )
      //{
      //        //std::cout << (uint64_t)run << " " << (uint64_t)ls << " " << (uint64_t)event << std::endl;
      //        std::cout << "\tPassEt: " << pass_et << " --> " << LooseEle1_SCEnergy/cosh(LooseEle1_SCEta)Heep << std::endl
      //          << "\tPassEcalDriven: " << pass_ecalDriven << " --> " << LooseEle1_EcalDriven << std::endl
      //          << "\tPassDeltaPhi: " << pass_deltaPhi << " --> " << LooseEle1_DeltaPhiTrkSC << std::endl
      //          << "\tPassMissingHits: " << pass_missingHits << " --> " << LooseEle1_MissingHits << std::endl
      //          << "\tPassDeltaEtaSeed: " << pass_deltaEtaSeed << " --> " << LooseEle1_DeltaEtaSeed << std::endl;
      //        if ( fabs(LooseEle1_SCEta) > 1.566 && fabs(LooseEle1_SCEta) < 2.5 )
      //          std::cout << "\tPassSigmaIetaIeta: " << pass_sigmaIEtaIEta << " --> " << LooseEle1_Full5x5SigmaIEtaIEta << " [ENDCAP]"<< std::endl;
      //        else
      //          std::cout << "\tPassSigmaIetaIeta: " << pass_sigmaIEtaIEta << " --> " << LooseEle1_Full5x5SigmaIEtaIEta << " [BARREL]"<< std::endl;
      //        if ( fabs(LooseEle1_SCEta) < 1.4442 )
      //        {
      //          std::cout << "\tPassShowerShape: " << pass_shape << " --> E1x5/E5x5=" << LooseEle1_Full5x5E1x5OverE5x5 << ", E2x5/E5x5=" << LooseEle1_Full5x5E2x5OverE5x5 << " [BARREL]" << std::endl;
      //          std::cout << "\tPassCaloIso: " << pass_caloIsolation << " --> " << caloIsolation << " < " << ( 2.0 + ( 0.03 * LooseEle1_SCEnergy/cosh(LooseEle1_SCEta)Heep ) + (0.28 * LooseEle1_RhoForHeep ) ) << std::endl;
      //          std::cout << "\tPassDxy: " << pass_dxy << " --> " << LooseEle1_LeadVtxDistXY << " < 0.02 " << std::endl
      //            << "\tPassH/E: " << pass_hoe << " --> " << LooseEle1_HoE << " < " << 1/LooseEle1_SCEnergy + 0.05 << std::endl;
      //        }
      //        else
      //        {
      //          std::cout << "\tPassCaloIso: " << pass_caloIsolation << " --> " << caloIsolation << " < " << (LooseEle1_SCEnergy/cosh(LooseEle1_SCEta)Heep<50 ? (2.5+( 0.28 * LooseEle1_RhoForHeep ) ) : (2.5+( 0.28 * LooseEle1_RhoForHeep ) + ( 0.03 * (LooseEle1_SCEnergy/cosh(LooseEle1_SCEta)Heep - 50.0 ) )) ) << std::endl;
      //          std::cout << "\tPassDxy: " << pass_dxy << " --> " << LooseEle1_LeadVtxDistXY << " < 0.05" << std::endl
      //            << "\tPassH/E: " << pass_hoe << " --> " << LooseEle1_HoE << " < " << 5/LooseEle1_SCEnergy + 0.05 << std::endl;
      //        }
      //        std::cout << "\tPassTrkIso: " << (LooseEle1_TrkIsoHEEP7 < 5.0) << " --> " << LooseEle1_TrkIsoHEEP7 << " < 5.0" << std::endl;
      //}
      ////XXX SIC TEST

      if ( passHEEPprime ) { 

        FillUserTH1D( "Electrons_nVertex_PAS"           , nVertex                      , min_prescale * pileup_weight);
        FillUserTH1D( "Electrons_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
        FillUserTH1D( "Electrons_Pt1stEle_PAS"	  , LooseEle1_Pt                 , min_prescale * pileup_weight);
        FillUserTH1D( "Electrons_hltPt1stEle_PAS"	  , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
        FillUserTH1D( "Electrons_Eta1stEle_PAS"	  , LooseEle1_Eta                , min_prescale * pileup_weight);
        FillUserTH1D( "Electrons_Phi1stEle_PAS"	  , LooseEle1_Phi                , min_prescale * pileup_weight);
        FillUserTH1D( "Electrons_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);
        FillUserTH1D( "Electrons_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        FillUserTH1D( "Electrons_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        FillUserTH2D( "Electrons_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        FillUserTH2D( "Electrons_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        if(analysisYear==2018) {
          if(run >= 319077 || !isData()) {
            FillUserTH2D( "Electrons_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
              FillUserTH2D( "Electrons_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            else
              FillUserTH2D( "Electrons_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          else
            FillUserTH2D( "Electrons_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }
        if(passHEEP) {
          FillUserTH1D( "ElectronsHEEP_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
          FillUserTH1D( "ElectronsHEEP_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
        }

        if ( isBarrel ) { 
          FillUserTH1D( "Electrons_Bar_nVertex_PAS"           , nVertex                      , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_Pt1stEle_PAS"	  , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_hltPt1stEle_PAS"	  , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_Eta1stEle_PAS"	  , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_Phi1stEle_PAS"	  , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_Bar_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Electrons_Bar_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Electrons_Bar_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Electrons_Bar_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Electrons_Bar_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Electrons_Bar_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Electrons_Bar_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          if(passHEEP) {
            FillUserTH1D( "ElectronsHEEP_Bar_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
            FillUserTH1D( "ElectronsHEEP_Bar_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut <= 1 ) {
            FillUserTH1D( "Electrons_Bar_lte1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_lte1Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_lte1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_lte1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_lte1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_lte1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_lte1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_Bar_lte1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Electrons_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Electrons_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Electrons_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Electrons_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            if(passHEEP) {
              FillUserTH1D( "ElectronsHEEP_Bar_lte1Jet_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
              FillUserTH1D( "ElectronsHEEP_Bar_lte1Jet_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 1 ) {
            FillUserTH1D( "Electrons_Bar_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_1Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_Bar_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_Bar_1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Electrons_Bar_1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Electrons_Bar_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Electrons_Bar_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Electrons_Bar_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            if(passHEEP) {
              FillUserTH1D( "ElectronsHEEP_Bar_1Jet_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
              FillUserTH1D( "ElectronsHEEP_Bar_1Jet_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 2 ) {
            FillUserTH1D( "Electrons_Bar_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_2Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_Bar_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_Bar_2Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Electrons_Bar_2Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Electrons_Bar_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Electrons_Bar_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Electrons_Bar_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            if(passHEEP) {
              FillUserTH1D( "ElectronsHEEP_Bar_2Jet_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
              FillUserTH1D( "ElectronsHEEP_Bar_2Jet_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 3 ) {
            FillUserTH1D( "Electrons_Bar_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_3Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_Bar_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_Bar_3Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Electrons_Bar_3Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Electrons_Bar_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Electrons_Bar_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Electrons_Bar_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            if(passHEEP) {
              FillUserTH1D( "ElectronsHEEP_Bar_3Jet_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
              FillUserTH1D( "ElectronsHEEP_Bar_3Jet_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
            }
          }

          if ( nVertex >= 0 && nVertex <= 5 ){
            FillUserTH1D( "Electrons_Bar_PU0-8_nElectron_PAS"          , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU0-8_Pt1stEle_PAS"           , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU0-8_Eta1stEle_PAS"          , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU0-8_Phi1stEle_PAS"          , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU0-8_Charge1stEle_PAS"       , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Electrons_Bar_PU0-8_MET_PAS"                , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }						          

          if ( nVertex >= 6 && nVertex <= 10 ){      
            FillUserTH1D( "Electrons_Bar_PU9-UP_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU9-UP_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU9-UP_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU9-UP_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_Bar_PU9-UP_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
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
          FillUserTH1D( "Electrons_End1_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End1_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Electrons_End1_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Electrons_End1_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Electrons_End1_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Electrons_End1_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Electrons_End1_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Electrons_End1_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          if(passHEEP) {
            FillUserTH1D( "ElectronsHEEP_End1_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
            FillUserTH1D( "ElectronsHEEP_End1_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut <= 1 ) {
            FillUserTH1D( "Electrons_End1_lte1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_lte1Jet_Pt1stEle_PAS"        , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_lte1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_lte1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_lte1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_lte1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_lte1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End1_lte1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Electrons_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Electrons_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Electrons_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Electrons_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            if(passHEEP) {
              FillUserTH1D( "ElectronsHEEP_End1_lte1Jet_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
              FillUserTH1D( "ElectronsHEEP_End1_lte1Jet_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
            }

          }

          if ( nJetLooseEle_ptCut >= 1 ) {
            FillUserTH1D( "Electrons_End1_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_1Jet_Pt1stEle_PAS"        , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End1_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End1_1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Electrons_End1_1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Electrons_End1_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Electrons_End1_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Electrons_End1_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            if(passHEEP) {
              FillUserTH1D( "ElectronsHEEP_End1_1Jet_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
              FillUserTH1D( "ElectronsHEEP_End1_1Jet_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
            }

          }

          if ( nJetLooseEle_ptCut >= 2 ) {
            FillUserTH1D( "Electrons_End1_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_2Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End1_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End1_2Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Electrons_End1_2Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Electrons_End1_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Electrons_End1_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Electrons_End1_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            if(passHEEP) {
              FillUserTH1D( "ElectronsHEEP_End1_2Jet_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
              FillUserTH1D( "ElectronsHEEP_End1_2Jet_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 3 ) {
            FillUserTH1D( "Electrons_End1_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_3Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End1_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End1_3Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Electrons_End1_3Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Electrons_End1_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Electrons_End1_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Electrons_End1_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            if(passHEEP) {
              FillUserTH1D( "ElectronsHEEP_End1_3Jet_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
              FillUserTH1D( "ElectronsHEEP_End1_3Jet_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
            }
          }

          if ( nVertex >= 0 && nVertex <= 5 ){
            FillUserTH1D( "Electrons_End1_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Electrons_End1_PU0-8_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }						          

          if ( nVertex >= 6 && nVertex <= 10 ){      
            FillUserTH1D( "Electrons_End1_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End1_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
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
          FillUserTH1D( "Electrons_End2_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Electrons_End2_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Electrons_End2_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Electrons_End2_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Electrons_End2_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Electrons_End2_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Electrons_End2_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Electrons_End2_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          if(passHEEP) {
            FillUserTH1D( "ElectronsHEEP_End2_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
            FillUserTH1D( "ElectronsHEEP_End2_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut <= 1 ) {
            FillUserTH1D( "Electrons_End2_lte1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_lte1Jet_Pt1stEle_PAS"        , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_lte1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_lte1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_lte1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_lte1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_lte1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End2_lte1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Electrons_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Electrons_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Electrons_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Electrons_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            if(passHEEP) {
              FillUserTH1D( "ElectronsHEEP_End2_lte1Jet_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
              FillUserTH1D( "ElectronsHEEP_End2_lte1Jet_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 1 ) {
            FillUserTH1D( "Electrons_End2_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_1Jet_Pt1stEle_PAS"        , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End2_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End2_1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Electrons_End2_1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Electrons_End2_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Electrons_End2_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Electrons_End2_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            if(passHEEP) {
              FillUserTH1D( "ElectronsHEEP_End2_1Jet_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
              FillUserTH1D( "ElectronsHEEP_End2_1Jet_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 2 ) {
            FillUserTH1D( "Electrons_End2_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_2Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End2_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End2_2Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Electrons_End2_2Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Electrons_End2_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Electrons_End2_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Electrons_End2_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            if(passHEEP) {
              FillUserTH1D( "ElectronsHEEP_End2_2Jet_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
              FillUserTH1D( "ElectronsHEEP_End2_2Jet_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 3 ) {
            FillUserTH1D( "Electrons_End2_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_3Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End2_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Electrons_End2_3Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Electrons_End2_3Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Electrons_End2_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Electrons_End2_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Electrons_End2_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            if(passHEEP) {
              FillUserTH1D( "ElectronsHEEP_End2_3Jet_MTenu_PAS"               , MT_Ele1MET            , min_prescale * pileup_weight);
              FillUserTH1D( "ElectronsHEEP_End2_3Jet_Pt1stEle_PAS"            , LooseEle1_Pt            , min_prescale * pileup_weight);
            }
          }

          if ( nVertex >= 0 && nVertex <= 5 ){
            FillUserTH1D( "Electrons_End2_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Electrons_End2_PU0-8_MET_PAS"        , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }						          

          if ( nVertex >= 6 && nVertex <= 10 ){      
            FillUserTH1D( "Electrons_End2_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Electrons_End2_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
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
      else if (passJet) {

        FillUserTH1D( "Jets_nVertex_PAS"           , nVertex                      , min_prescale * pileup_weight);
        FillUserTH1D( "Jets_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
        FillUserTH1D( "Jets_Pt1stEle_PAS"	  , LooseEle1_Pt                 , min_prescale * pileup_weight);
        FillUserTH1D( "Jets_hltPt1stEle_PAS"	  , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
        FillUserTH1D( "Jets_Eta1stEle_PAS"	  , LooseEle1_Eta                , min_prescale * pileup_weight);
        FillUserTH1D( "Jets_Phi1stEle_PAS"	  , LooseEle1_Phi                , min_prescale * pileup_weight);
        FillUserTH1D( "Jets_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);
        FillUserTH1D( "Jets_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
        FillUserTH1D( "Jets_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        FillUserTH2D( "Jets_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        FillUserTH2D( "Jets_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        if(analysisYear==2018) {
          if(run >= 319077 || !isData()) {
            FillUserTH2D( "Jets_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
              FillUserTH2D( "Jets_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            else
              FillUserTH2D( "Jets_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          else
            FillUserTH2D( "Jets_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        }
        FillUserTH2D( "Jets_TrkIsoHEEP7vsMTenu_PAS"   , MT_Ele1MET, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
        // for SR/SR division
        if(LooseEle1_TrkIsoHEEP7 < 5) {
          FillUserTH1D( "JetsSR_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "JetsSR_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
        }
        else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
        {
          FillUserTH1D( "JetsSB_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "JetsSB_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
        }

        if ( isBarrel ) { 
          FillUserTH1D( "Jets_Bar_nVertex_PAS"           , nVertex                      , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_Pt1stEle_PAS"	  , LooseEle1_Pt                 , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_hltPt1stEle_PAS"	  , LooseEle1_hltPhotonPt        , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_Eta1stEle_PAS"	  , LooseEle1_Eta                , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_Phi1stEle_PAS"	  , LooseEle1_Phi                , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_Bar_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Jets_Bar_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Jets_Bar_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Jets_Bar_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Jets_Bar_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Jets_Bar_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Jets_Bar_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          FillUserTH2D( "Jets_Bar_TrkIsoHEEP7vsMTenu_PAS"   , MT_Ele1MET, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          // for SR/SR division
          if(LooseEle1_TrkIsoHEEP7 < 5) {
            FillUserTH1D( "JetsSR_Bar_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "JetsSR_Bar_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
          }
          else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
          {
            FillUserTH1D( "JetsSB_Bar_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "JetsSB_Bar_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut <= 1 ) {
            FillUserTH1D( "Jets_Bar_lte1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_lte1Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_lte1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_lte1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_lte1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_lte1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_lte1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_Bar_lte1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Jets_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Jets_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Jets_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Jets_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            FillUserTH2D( "Jets_Bar_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            // for SR/SR division
            if(LooseEle1_TrkIsoHEEP7 < 5) {
              FillUserTH1D( "JetsSR_Bar_lte1Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSR_Bar_lte1Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
            else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
            {
              FillUserTH1D( "JetsSB_Bar_lte1Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSB_Bar_lte1Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 1 ) {
            FillUserTH1D( "Jets_Bar_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_1Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_Bar_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_Bar_1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Jets_Bar_1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Jets_Bar_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Jets_Bar_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Jets_Bar_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            // for SR/SR division
            if(LooseEle1_TrkIsoHEEP7 < 5) {
              FillUserTH1D( "JetsSR_Bar_1Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSR_Bar_1Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
            else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
            {
              FillUserTH1D( "JetsSB_Bar_1Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSB_Bar_1Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 2 ) {
            FillUserTH1D( "Jets_Bar_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_2Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_Bar_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_Bar_2Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Jets_Bar_2Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Jets_Bar_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Jets_Bar_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Jets_Bar_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            // for SR/SR division
            if(LooseEle1_TrkIsoHEEP7 < 5) {
              FillUserTH1D( "JetsSR_Bar_2Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSR_Bar_2Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
            else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
            {
              FillUserTH1D( "JetsSB_Bar_2Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSB_Bar_2Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 3 ) {
            FillUserTH1D( "Jets_Bar_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_3Jet_Pt1stEle_PAS"	       , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_Bar_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_Bar_3Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Jets_Bar_3Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Jets_Bar_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Jets_Bar_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Jets_Bar_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            // for SR/SR division
            if(LooseEle1_TrkIsoHEEP7 < 5) {
              FillUserTH1D( "JetsSR_Bar_3Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSR_Bar_3Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
            else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
            {
              FillUserTH1D( "JetsSB_Bar_3Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSB_Bar_3Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
          }

          if ( nVertex >= 0 && nVertex <= 5 ){
            FillUserTH1D( "Jets_Bar_PU0-8_nElectron_PAS"          , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU0-8_Pt1stEle_PAS"           , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU0-8_Eta1stEle_PAS"          , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU0-8_Phi1stEle_PAS"          , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU0-8_Charge1stEle_PAS"       , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Jets_Bar_PU0-8_MET_PAS"                , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }						          

          if ( nVertex >= 6 && nVertex <= 10 ){      
            FillUserTH1D( "Jets_Bar_PU9-UP_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU9-UP_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU9-UP_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU9-UP_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_Bar_PU9-UP_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
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
          FillUserTH1D( "Jets_End1_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End1_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Jets_End1_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Jets_End1_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Jets_End1_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Jets_End1_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Jets_End1_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Jets_End1_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          // for SR/SR division
          if(LooseEle1_TrkIsoHEEP7 < 5) {
            FillUserTH1D( "JetsSR_End1_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "JetsSR_End1_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
          }
          else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
          {
            FillUserTH1D( "JetsSB_End1_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "JetsSB_End1_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut <= 1 ) {
            FillUserTH1D( "Jets_End1_lte1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_lte1Jet_Pt1stEle_PAS"        , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_lte1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_lte1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_lte1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_lte1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_lte1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End1_lte1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Jets_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Jets_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Jets_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Jets_End1_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            // for SR/SR division
            if(LooseEle1_TrkIsoHEEP7 < 5) {
              FillUserTH1D( "JetsSR_End1_lte1Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSR_End1_lte1Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
            else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
            {
              FillUserTH1D( "JetsSB_End1_lte1Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSB_End1_lte1Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 1 ) {
            FillUserTH1D( "Jets_End1_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_1Jet_Pt1stEle_PAS"        , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End1_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End1_1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Jets_End1_1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Jets_End1_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Jets_End1_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Jets_End1_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            // for SR/SR division
            if(LooseEle1_TrkIsoHEEP7 < 5) {
              FillUserTH1D( "JetsSR_End1_1Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSR_End1_1Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
            else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
            {
              FillUserTH1D( "JetsSB_End1_1Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSB_End1_1Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 2 ) {
            FillUserTH1D( "Jets_End1_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_2Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End1_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End1_2Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Jets_End1_2Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Jets_End1_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Jets_End1_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Jets_End1_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            // for SR/SR division
            if(LooseEle1_TrkIsoHEEP7 < 5) {
              FillUserTH1D( "JetsSR_End1_2Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSR_End1_2Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
            else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
            {
              FillUserTH1D( "JetsSB_End1_2Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSB_End1_2Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 3 ) {
            FillUserTH1D( "Jets_End1_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_3Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End1_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End1_3Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Jets_End1_3Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Jets_End1_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Jets_End1_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Jets_End1_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            // for SR/SR division
            if(LooseEle1_TrkIsoHEEP7 < 5) {
              FillUserTH1D( "JetsSR_End1_3Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSR_End1_3Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
            else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
            {
              FillUserTH1D( "JetsSB_End1_3Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSB_End1_3Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
          }

          if ( nVertex >= 0 && nVertex <= 5 ){
            FillUserTH1D( "Jets_End1_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Jets_End1_PU0-8_MET_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }						          

          if ( nVertex >= 6 && nVertex <= 10 ){      
            FillUserTH1D( "Jets_End1_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End1_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
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
          FillUserTH1D( "Jets_End2_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          FillUserTH1D( "Jets_End2_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Jets_End2_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          FillUserTH2D( "Jets_End2_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          if(analysisYear==2018) {
            if(run >= 319077 || !isData()) {
              FillUserTH2D( "Jets_End2_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                FillUserTH2D( "Jets_End2_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              else
                FillUserTH2D( "Jets_End2_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            else
              FillUserTH2D( "Jets_End2_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
          }
          // for SR/SR division
          if(LooseEle1_TrkIsoHEEP7 < 5) {
            FillUserTH1D( "JetsSR_End2_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "JetsSR_End2_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
          }
          else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
          {
            FillUserTH1D( "JetsSB_End2_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "JetsSB_End2_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
          }

          if ( nJetLooseEle_ptCut <= 1 ) {
            FillUserTH1D( "Jets_End2_lte1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_lte1Jet_Pt1stEle_PAS"        , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_lte1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_lte1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_lte1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_lte1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_lte1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End2_lte1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Jets_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Jets_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Jets_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Jets_End2_lte1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            // for SR/SR division
            if(LooseEle1_TrkIsoHEEP7 < 5) {
              FillUserTH1D( "JetsSR_End2_lte1Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSR_End2_lte1Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
            else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
            {
              FillUserTH1D( "JetsSB_End2_lte1Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSB_End2_lte1Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 1 ) {
            FillUserTH1D( "Jets_End2_1Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_1Jet_Pt1stEle_PAS"        , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_1Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_1Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_1Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_1Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_1Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End2_1Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End2_1Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Jets_End2_1Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Jets_End2_1Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Jets_End2_1Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Jets_End2_1Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            // for SR/SR division
            if(LooseEle1_TrkIsoHEEP7 < 5) {
              FillUserTH1D( "JetsSR_End2_1Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSR_End2_1Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
            else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
            {
              FillUserTH1D( "JetsSB_End2_1Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSB_End2_1Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 2 ) {
            FillUserTH1D( "Jets_End2_2Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_2Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_2Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_2Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_2Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_2Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_2Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End2_2Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End2_2Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Jets_End2_2Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Jets_End2_2Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Jets_End2_2Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Jets_End2_2Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            // for SR/SR division
            if(LooseEle1_TrkIsoHEEP7 < 5) {
              FillUserTH1D( "JetsSR_End2_2Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSR_End2_2Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
            else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
            {
              FillUserTH1D( "JetsSB_End2_2Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSB_End2_2Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
          }

          if ( nJetLooseEle_ptCut >= 3 ) {
            FillUserTH1D( "Jets_End2_3Jet_nElectron_PAS"       , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_3Jet_Pt1stEle_PAS"	, LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_3Jet_Eta1stEle_PAS"       , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_3Jet_Phi1stEle_PAS"       , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_3Jet_Charge1stEle_PAS"    , LooseEle1_Charge             , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_3Jet_MET_PAS"             , PFMET_Type1_Pt            , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_3Jet_TrkIsoHEEP7_PAS"       , LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End2_3Jet_TrkIsoHEEP7vsPt_PAS"   , LooseEle1_Pt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            FillUserTH2D( "Jets_End2_3Jet_TrkIsoHEEP7vsHLTPt_PAS"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            if(analysisYear==2018) {
              if(run >= 319077 || !isData()) {
                FillUserTH2D( "Jets_End2_3Jet_TrkIsoHEEP7vsHLTPt_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                if(isHEMElectron(LooseEle1_SCEta, LooseEle1_Phi))
                  FillUserTH2D( "Jets_End2_3Jet_TrkIsoHEEP7vsHLTPt_HEMonly_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
                else
                  FillUserTH2D( "Jets_End2_3Jet_TrkIsoHEEP7vsHLTPt_noHEM_post319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
              }
              else
                FillUserTH2D( "Jets_End2_3Jet_TrkIsoHEEP7vsHLTPt_pre319077"   , LooseEle1_hltPhotonPt, LooseEle1_TrkIsoHEEP7            , min_prescale * pileup_weight);
            }
            // for SR/SR division
            if(LooseEle1_TrkIsoHEEP7 < 5) {
              FillUserTH1D( "JetsSR_End2_3Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSR_End2_3Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
            else if(LooseEle1_TrkIsoHEEP7 >= 10 && LooseEle1_TrkIsoHEEP7 < 20)
            {
              FillUserTH1D( "JetsSB_End2_3Jet_MTenu_PAS"               , PFMET_Type1_Pt            , min_prescale * pileup_weight);
              FillUserTH1D( "JetsSB_End2_3Jet_Pt1stEle_PAS"            , LooseEle1_Pt           , min_prescale * pileup_weight);
            }
          }

          if ( nVertex >= 0 && nVertex <= 5 ){
            FillUserTH1D( "Jets_End2_PU0-8_nElectron_PAS"         , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU0-8_Pt1stEle_PAS"          , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU0-8_Eta1stEle_PAS"         , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU0-8_Phi1stEle_PAS"         , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU0-8_Charge1stEle_PAS"      , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Jets_End2_PU0-8_MET_PAS"        , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }						          

          if ( nVertex >= 6 && nVertex <= 10 ){      
            FillUserTH1D( "Jets_End2_PU9-UP_nElectron_PAS"        , nLooseEle_ptCut              , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU9-UP_Pt1stEle_PAS"         , LooseEle1_Pt                 , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU9-UP_Eta1stEle_PAS"        , LooseEle1_Eta                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU9-UP_Phi1stEle_PAS"        , LooseEle1_Phi                , min_prescale * pileup_weight);
            FillUserTH1D( "Jets_End2_PU9-UP_Charge1stEle_PAS"     , LooseEle1_Charge             , min_prescale * pileup_weight);  
            FillUserTH1D( "Jets_End2_PU9-UP_MET_PAS"              , PFMET_Type1_Pt            , min_prescale * pileup_weight);
          }
        }
      }
    }
  } // End loop over events

  std::cout << "analysisClass::Loop() ends" <<std::endl;   
}
