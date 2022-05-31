#define analysisClass_cxx
#define USE_QCD_REDUCED_NTUPLE
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
// for fake rate
#include "include/QCDFakeRate.h"
// for prescales
#include "include/Run2PhotonTriggerPrescales.h"

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile){}

  analysisClass::~analysisClass(){}

void analysisClass::Loop()
{
  std::cout << "analysisClass::Loop() begins" <<std::endl;   

  //--------------------------------------------------------------------------
  // Decide which plots to save (default is to save everything)
  //--------------------------------------------------------------------------

  fillAllPreviousCuts              ( !true  ) ;
  fillAllOtherCuts                 ( !true  ) ;
  fillAllSameLevelAndLowerLevelCuts( !true  ) ;
  fillAllCuts                      ( !true  ) ;

  //--------------------------------------------------------------------------
  // Get pre-cut values
  //--------------------------------------------------------------------------
  // eta boundaries

  double eleEta_bar            = getPreCutValue1("eleEta_bar");
  double eleEta_end1_min       = getPreCutValue1("eleEta_end1");
  double eleEta_end1_max       = getPreCutValue2("eleEta_end1");
  double eleEta_end2_min       = getPreCutValue1("eleEta_end2");
  double eleEta_end2_max       = getPreCutValue2("eleEta_end2");

  // prescales
  Run2PhotonTriggerPrescales run2PhotonTriggerPrescales;

  //--------------------------------------------------------------------------
  // Analysis year
  //--------------------------------------------------------------------------
  int analysisYear = getPreCutValue1("AnalysisYear");

  //--------------------------------------------------------------------------
  // QCD Fake Rate loading part
  //--------------------------------------------------------------------------
  std::string qcdFileName = getPreCutString1("QCDFakeRateFilename");
  //HistoReader qcdFakeRateReader(qcdFileName,"fr2D_Bar_2Jet","fr2D_End_2Jet",true,false);
  std::vector<std::string> regionVec;
  if(analysisYear < 2018)
    regionVec = {"Bar_2Jet", "End1_2Jet", "End2_2Jet"};
  else
    regionVec = {
      "Bar_pre319077_2Jet",
      "End1_pre319077_2Jet",
      "End2_pre319077_2Jet",
      "Bar_noHEM_post319077_2Jet",
      "End1_noHEM_post319077_2Jet",
      "End2_noHEM_post319077_2Jet",
      "Bar_HEMonly_post319077_2Jet",
      "End1_HEMonly_post319077_2Jet",
      "End2_HEMonly_post319077_2Jet"};
  QCDFakeRate qcdFR(qcdFileName, regionVec, analysisYear);

  //--------------------------------------------------------------------------
  // Create TH1D's
  //--------------------------------------------------------------------------

  CreateUserTH1D( "EventCount"            ,    1 , 0       , 1	 );    

  CreateUserTH1D( "M_j1j3_PAS"            ,    200 , 0       , 2000	 );    
  CreateUserTH1D( "M_j2j3_PAS"            ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "M_e1j3_PAS"            ,    200 , 0       , 2000	 );    
  CreateUserTH1D( "M_e2j3_PAS"            ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "M_eejjj_PAS"           ,    500 , 0       , 5000	 ); 

  CreateUserTH1D( "M_j1j3_PASandMee100"   ,    200 , 0       , 2000	 );    
  CreateUserTH1D( "M_j2j3_PASandMee100"   ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "M_e1j3_PASandMee100"   ,    200 , 0       , 2000	 );    
  CreateUserTH1D( "M_e2j3_PASandMee100"   ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "M_eejjj_PASandMee100"  ,    500 , 0       , 5000	 ); 

  CreateUserTH1D( "M_j1j3_ROI"            ,    200 , 0       , 2000	 );    
  CreateUserTH1D( "M_j2j3_ROI"            ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "M_e1j3_ROI"            ,    200 , 0       , 2000	 );    
  CreateUserTH1D( "M_e2j3_ROI"            ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "M_eejjj_ROI"           ,    500 , 0       , 5000	 ); 

  CreateUserTH1D( "sTfrac_Jet1_PAS"       ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Jet2_PAS"       ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Ele1_PAS"       ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Ele2_PAS"       ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Jet_PAS"        ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Ele_PAS"        ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Jet1_PASandMee100"       ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Jet2_PASandMee100"       ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Ele1_PASandMee100"       ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Ele2_PASandMee100"       ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Jet_PASandMee100"        ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Ele_PASandMee100"        ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Jet1_ROI"       ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Jet2_ROI"       ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Ele1_ROI"       ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Ele2_ROI"       ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Jet_ROI"        ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "sTfrac_Ele_ROI"        ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "nMuon_PAS"             ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "nJet_PAS"              ,    10  , -0.5    , 9.5      );
  CreateUserTH1D( "nJet_PASandMee100"        ,    10  , -0.5    , 9.5      );
  CreateUserTH1D( "nJet_ROI"              ,    10  , -0.5    , 9.5      );
  CreateUserTH1D( "Pt1stEle_PAS"	   , 	300 , 0       , 3000     ); 
  CreateUserTH1D( "Pt1stEle_PASandMee100" , 	300 , 0       , 3000     ); 
  CreateUserTH1D( "Pt1stEle_ROI"	   , 	300 , 0       , 3000     ); 
  CreateUserTH1D( "Eta1stEle_PAS"	   , 	100 , -5      , 5	 ); 
  CreateUserTH1D( "Phi1stEle_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserTH1D( "Pt2ndEle_PAS"	   , 	300 , 0       , 3000     ); 
  CreateUserTH1D( "Pt2ndEle_PASandMee100" , 	300 , 0       , 3000     ); 
  CreateUserTH1D( "Pt2ndEle_ROI"	   , 	300 , 0       , 3000     ); 
  CreateUserTH1D( "Eta2ndEle_PAS"	   , 	100 , -5      , 5	 ); 
  CreateUserTH1D( "Phi2ndEle_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserTH1D( "Charge1stEle_PAS"	   , 	2   , -1.0001 , 1.0001	 ); 
  CreateUserTH1D( "Charge2ndEle_PAS"	   , 	2   , -1.0001 , 1.0001	 ); 
  CreateUserTH1D( "EleChargeSum_PAS"         ,    3   , -2.5    , 2.5  );
  CreateUserTH1D( "EleChargeSum_PASandMee100",    3   , -2.5    , 2.5  );
  CreateUserTH1D( "EleChargeSum_ROI"         ,    3   , -2.5    , 2.5  );
  CreateUserTH1D( "MET_PAS"               ,    200 , 0       , 1000	 ); 
  CreateUserTH1D( "MET_ROI"               ,    200 , 0       , 1000	 ); 
  CreateUserTH1D( "METPhi_PAS"		   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserTH1D( "Pt1stJet_PAS"          ,    300 , 0       , 3000	 ); 
  CreateUserTH1D( "Pt2ndJet_PAS"          ,    300 , 0       , 3000	 ); 
  CreateUserTH1D( "Pt1stJet_PASandMee100" ,    300 , 0       , 3000	 ); 
  CreateUserTH1D( "Pt2ndJet_PASandMee100" ,    300 , 0       , 3000	 ); 
  CreateUserTH1D( "Pt1stJet_ROI"          ,    300 , 0       , 3000	 ); 
  CreateUserTH1D( "Pt2ndJet_ROI"          ,    300 , 0       , 3000	 ); 
  CreateUserTH1D( "Eta1stJet_PAS"         ,    100 , -5      , 5	 ); 
  CreateUserTH1D( "Eta2ndJet_PAS"         ,    100 , -5      , 5	 ); 
  CreateUserTH1D( "Phi1stJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserTH1D( "Phi2ndJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserTH1D( "sTlep_PAS"             ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sTjet_PAS"             ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sTlep_PASandMee100"    ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sTjet_PASandMee100"    ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sTlep_ROI"             ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sTjet_ROI"             ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sT_PAS"                ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sT_zjj_PAS"                      ,    400   , 0       , 4000	  ); 
  CreateUserTH1D( "sT_zjj_PASandMee100"             ,    400   , 0       , 4000	  ); 
  CreateUserTH1D( "sT_zjj_ROI"                      ,    400   , 0       , 4000	  ); 
  CreateUserTH1D( "sT_PASandMee100"       ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sT_PASandMee110"       ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sT_PASandMee120"       ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sT_PASandMee130"       ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sT_PASandMee140"       ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sT_PASandMee150"       ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sT_PASandMee160"       ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sT_PASandMee170"       ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sT_PASandMee180"       ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sT_PASandMee190"       ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sT_PASandMee200"       ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "sT_ROI"                ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "Mjj_PAS"		   ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "Mjj_PASandMee100"	   ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "Mjj_ROI"		   ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "Mee_PAS"		   ,    2000 , 0       , 2000	 ); 
  CreateUserTH1D( "Mee_ROI"		   ,    2000 , 0       , 2000	 ); 
  CreateUserTH1D( "Mee_PASandST445"       ,    4000 , 0       , 4000	 ); 
  CreateUserTH1D( "MTenu_PAS"             ,    400 , 0       , 2000	 ); 
  CreateUserTH1D( "Me1j1_PAS"             ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "Me1j1_PASandMee100"    ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "Me1j1_ROI"             ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "Me1j2_PAS"             ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "Me2j1_PAS"             ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "Me2j2_PAS"             ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "Me1j_selected_PAS"     ,    400 , 0       , 4000	 ); 
  CreateUserTH1D( "Me2j_selected_PAS"     ,    400 , 0       , 4000     );
  CreateUserTH1D( "Mej_selected_min_PAS"  ,    400 , 0       , 4000     ); 
  CreateUserTH1D( "Mej_selected_max_PAS"  ,    400 , 0       , 4000     ); 
  CreateUserTH1D( "Mej_minmax_PAS"        ,    400 , 0       , 4000     ); 
  CreateUserTH1D( "Mej_selected_avg_PAS"  ,    400 , 0       , 4000     );
  CreateUserTH1D( "Mej_selected_avg_PASandMee100"  ,    400 , 0       , 4000     );
  CreateUserTH1D( "Mej_selected_avg_ROI"  ,    400 , 0       , 4000     );
  CreateUserTH1D( "Mejj_PAS"              ,    600 , 0       , 6000     );
  CreateUserTH1D( "Meej_PAS"              ,    600 , 0       , 6000     );
  CreateUserTH1D( "Meejj_ROI"             ,    600 , 0       , 6000     );
  CreateUserTH1D( "Mejj_ROI"              ,    600 , 0       , 6000     );
  CreateUserTH1D( "Meej_ROI"              ,    600 , 0       , 6000     );
  CreateUserTH1D( "Meejj_PAS"             ,    600 , 0       , 6000     );

  CreateUserTH1D( "Eta1stJet_ROI"                   ,    100   , -5      , 5	  ); 
  CreateUserTH1D( "Eta2ndJet_ROI"                   ,    100   , -5      , 5	  ); 
  CreateUserTH1D( "Eta1stEle_ROI"	             , 	100    , -5      , 5	  ); 
  CreateUserTH1D( "Eta2ndEle_ROI"	             , 	100    , -5      , 5	  ); 

  CreateUserTH1D( "Phi1stJet_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
  CreateUserTH1D( "Phi2ndJet_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
  CreateUserTH1D( "Phi1stEle_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
  CreateUserTH1D( "Phi2ndEle_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 

  CreateUserTH1D( "Ptj1j2j3_PAS"                    ,    400 , 0       , 4000     );
  CreateUserTH1D( "Ptj1j2_PAS"                      ,    400 , 0       , 4000     );
  CreateUserTH1D( "Ptj2j3_PAS"                      ,    400 , 0       , 4000     );
  CreateUserTH1D( "Ptj1j3_PAS"                      ,    400 , 0       , 4000     );

  CreateUserTH1D( "Ptee_Minus_Ptj1j2_PAS"           ,    200 , -500    , 500      );
  CreateUserTH1D( "Ptee_Minus_Ptj1j2j3_PAS"         ,    200 , -500    , 500      );


  CreateUserTH1D( "Ptj1j2j3_PASandMee100"           ,    200 , 0       , 2000     );
  CreateUserTH1D( "Ptj1j2_PASandMee100"             ,    200 , 0       , 2000     );
  CreateUserTH1D( "Ptj2j3_PASandMee100"             ,    200 , 0       , 2000     );
  CreateUserTH1D( "Ptj1j3_PASandMee100"             ,    200 , 0       , 2000     );

  CreateUserTH1D( "Ptee_Minus_Ptj1j2_PASandMee100"  ,    200 , -500    , 500      );
  CreateUserTH1D( "Ptee_Minus_Ptj1j2j3_PASandMee100",    200 , -500    , 500      );


  CreateUserTH1D( "Ptj1j2j3_ROI"                    ,    200 , 0       , 2000     );
  CreateUserTH1D( "Ptj1j2_ROI"                      ,    200 , 0       , 2000     );
  CreateUserTH1D( "Ptj2j3_ROI"                      ,    200 , 0       , 2000     );
  CreateUserTH1D( "Ptj1j3_ROI"                      ,    200 , 0       , 2000     );

  CreateUserTH1D( "Ptee_Minus_Ptj1j2_ROI"           ,    200 , -500    , 500      );
  CreateUserTH1D( "Ptee_Minus_Ptj1j2j3_ROI"         ,    200 , -500    , 500      );


  CreateUserTH1D( "Ptee_PAS"              ,    200 , 0       , 2000     );
  CreateUserTH1D( "Ptee_PASandMee100"     ,    200 , 0       , 2000     );
  CreateUserTH1D( "Ptee_ROI"              ,    200 , 0       , 2000     );

  CreateUserTH1D( "nVertex_PAS"                     ,    101   , -0.5   , 100.5	 ) ; 
  CreateUserTH1D( "nVertex_PASandMee100"            ,    101   , -0.5   , 100.5	 ) ; 
  CreateUserTH1D( "nVertex_ROI"                     ,    101   , -0.5   , 100.5	 ) ; 

  CreateUserTH1D( "DR_Ele1Jet1_PAS"	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
  CreateUserTH1D( "DR_Ele1Jet2_PAS"	   , 	getHistoNBins("DR_Ele1Jet2"), getHistoMin("DR_Ele1Jet2"), getHistoMax("DR_Ele1Jet2")     ) ; 
  CreateUserTH1D( "DR_Ele2Jet1_PAS"	   , 	getHistoNBins("DR_Ele2Jet1"), getHistoMin("DR_Ele2Jet1"), getHistoMax("DR_Ele2Jet1")     ) ; 
  CreateUserTH1D( "DR_Ele2Jet2_PAS"	   , 	getHistoNBins("DR_Ele2Jet2"), getHistoMin("DR_Ele2Jet2"), getHistoMax("DR_Ele2Jet2")     ) ; 
  CreateUserTH1D( "DR_Jet1Jet2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  CreateUserTH1D( "DR_Ele1Ele2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  CreateUserTH1D( "minDR_EleJet_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  CreateUserTH1D( "minDR_ZJet_PAS"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  CreateUserTH1D( "minDR_ZJet_ROI"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 

  CreateUserTH1D( "DR_ZJet1_PAS"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  CreateUserTH1D( "DR_ZJet1_ROI"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  CreateUserTH1D( "DR_ZJet2_PAS"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  CreateUserTH1D( "DR_ZJet2_ROI"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 


  CreateUserTH2D( "Me1jVsMe2j_selected",     200, 0, 2000, 200, 0, 2000) ;
  CreateUserTH2D( "Me1jVsMe2j_rejected",     200, 0, 2000, 200, 0, 2000) ;


  CreateUserTH2D( "MeeVsST_PAS"                 ,     200, 0, 2000, 200, 0, 2000) ;
  CreateUserTH2D( "MeeVsST_PASandMee100"        ,     200, 0, 2000, 200, 0, 2000) ;
  CreateUserTH2D( "MeeVsST_ROI"                 ,     200, 0, 2000, 200, 0, 2000) ;


  CreateUserTH1D( "Mee_80_100_Preselection", 200, 60, 120 );
  CreateUserTH1D( "Mee_70_110_Preselection", 200, 60, 120 );
  CreateUserTH1D( "Mee_70_110_ST600_Preselection", 200, 60, 120 );

  CreateUserTH1D( "Mee_EBEB_PAS"		   ,    2000 , 0       , 2000	 ); 
  CreateUserTH1D( "Mee_EBEE_PAS"		   ,    2000 , 0       , 2000	 ); 
  CreateUserTH1D( "Mee_EEEE_PAS"		   ,    2000 , 0       , 2000	 ); 
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
  CreateUserTH1D( "GeneratorWeight", 100, -2.0 , 2.0 );

  CreateUserTH1D("CorrIsolation_1stEle_PAS"                 , 200,-25.0 ,  25.0  ); CreateUserTH1D("CorrIsolation_2ndEle_PAS"                 , 200,-25.0 ,  25.0  );
  CreateUserTH1D("DeltaEtaTrkSC_1stEle_PAS"                 , 200, -0.01,   0.01 ); CreateUserTH1D("DeltaEtaTrkSC_2ndEle_PAS"                 , 200, -0.01,   0.01 );
  CreateUserTH1D("E2x5OverE5x5_1stEle_PAS"                  , 200,  0.0 ,   2.0  ); CreateUserTH1D("E2x5OverE5x5_2ndEle_PAS"                  , 200,  0.0 ,   2.0  );
  CreateUserTH1D("EcalIsolation_1stEle_PAS"                 , 200,  0.0 ,  20.0  ); CreateUserTH1D("EcalIsolation_2ndEle_PAS"                 , 200,  0.0 ,  20.0  );
  CreateUserTH1D("HcalIsolation_1stEle_PAS"                 , 200,  0.0 ,  20.0  ); CreateUserTH1D("HcalIsolation_2ndEle_PAS"                 , 200,  0.0 ,  20.0  );
  CreateUserTH1D("TrkIsolation_1stEle_PAS"                  , 200,  0.0,    5.0  ); CreateUserTH1D("TrkIsolation_2ndEle_PAS"                  , 200,  0.0,    5.0  );
  CreateUserTH1D("GsfScPixCharge_1stEle_PAS"                , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfScPixCharge_2ndEle_PAS"                , 2,   -0.5 ,   1.5  );
  CreateUserTH1D("HasMatchedPhot_1stEle_PAS"                , 2,   -0.5 ,   1.5  ); CreateUserTH1D("HasMatchedPhot_2ndEle_PAS"                , 2,   -0.5 ,   1.5  );
  CreateUserTH1D("HoE_1stEle_PAS"                           , 200,  0.0 ,   0.05 ); CreateUserTH1D("HoE_2ndEle_PAS"                           , 200,  0.0 ,   0.05 );
  CreateUserTH1D("LeadVtxDistXY_1stEle_PAS"                 , 200, -0.05,   0.05 ); CreateUserTH1D("LeadVtxDistXY_2ndEle_PAS"                 , 200, -0.05,   0.05 );
  CreateUserTH1D("LeadVtxDistZ_1stEle_PAS"                  , 200, -0.2 ,   0.2  ); CreateUserTH1D("LeadVtxDistZ_2ndEle_PAS"                  , 200, -0.2 ,   0.2  );
  CreateUserTH1D("MissingHits_1stEle_PAS"                   , 2  , -0.5,    1.5  ); CreateUserTH1D("MissingHits_2ndEle_PAS"                   , 2  , -0.5,    1.5  );
  CreateUserTH1D("SigmaIEtaIEta_Barrel_1stEle_PAS"          , 200,  0.0,    0.04 ); CreateUserTH1D("SigmaIEtaIEta_Barrel_2ndEle_PAS"          , 200,  0.0,    0.04 );
  CreateUserTH1D("SigmaIEtaIEta_Endcap_1stEle_PAS"          , 200,  0.0,    0.1  ); CreateUserTH1D("SigmaIEtaIEta_Endcap_2ndEle_PAS"          , 200,  0.0,    0.1  );

  CreateUserTH1D("CorrIsolation_1stEle_PASandMee100"        , 200,-25.0 ,  25.0  ); CreateUserTH1D("CorrIsolation_2ndEle_PASandMee100"        , 200,-25.0 ,  25.0  );
  CreateUserTH1D("DeltaEtaTrkSC_1stEle_PASandMee100"        , 200, -0.01,   0.01 ); CreateUserTH1D("DeltaEtaTrkSC_2ndEle_PASandMee100"        , 200, -0.01,   0.01 );
  CreateUserTH1D("E2x5OverE5x5_1stEle_PASandMee100"         , 200,  0.0 ,   2.0  ); CreateUserTH1D("E2x5OverE5x5_2ndEle_PASandMee100"         , 200,  0.0 ,   2.0  );
  CreateUserTH1D("EcalIsolation_1stEle_PASandMee100"        , 200,  0.0 ,  20.0  ); CreateUserTH1D("EcalIsolation_2ndEle_PASandMee100"        , 200,  0.0 ,  20.0  );
  CreateUserTH1D("HcalIsolation_1stEle_PASandMee100"        , 200,  0.0 ,  20.0  ); CreateUserTH1D("HcalIsolation_2ndEle_PASandMee100"        , 200,  0.0 ,  20.0  );
  CreateUserTH1D("TrkIsolation_1stEle_PASandMee100"         , 200,  0.0,    5.0  ); CreateUserTH1D("TrkIsolation_2ndEle_PASandMee100"         , 200,  0.0,    5.0  );
  CreateUserTH1D("HasMatchedPhot_1stEle_PASandMee100"       , 2,   -0.5 ,   1.5  ); CreateUserTH1D("HasMatchedPhot_2ndEle_PASandMee100"       , 2,   -0.5 ,   1.5  );
  CreateUserTH1D("HoE_1stEle_PASandMee100"                  , 200,  0.0 ,   0.05 ); CreateUserTH1D("HoE_2ndEle_PASandMee100"                  , 200,  0.0 ,   0.05 );
  CreateUserTH1D("LeadVtxDistXY_1stEle_PASandMee100"        , 200, -0.05,   0.05 ); CreateUserTH1D("LeadVtxDistXY_2ndEle_PASandMee100"        , 200, -0.05,   0.05 );
  CreateUserTH1D("LeadVtxDistZ_1stEle_PASandMee100"         , 200, -0.2 ,   0.2  ); CreateUserTH1D("LeadVtxDistZ_2ndEle_PASandMee100"         , 200, -0.2 ,   0.2  );
  CreateUserTH1D("MissingHits_1stEle_PASandMee100"          , 2  , -0.5,    1.5  ); CreateUserTH1D("MissingHits_2ndEle_PASandMee100"          , 2  , -0.5,    1.5  );
  CreateUserTH1D("SigmaIEtaIEta_Barrel_1stEle_PASandMee100" , 200,  0.0,    0.02 ); CreateUserTH1D("SigmaIEtaIEta_Barrel_2ndEle_PASandMee100" , 200,  0.0,    0.02 );
  CreateUserTH1D("SigmaIEtaIEta_Endcap_1stEle_PASandMee100" , 200,  0.0,    0.1  ); CreateUserTH1D("SigmaIEtaIEta_Endcap_2ndEle_PASandMee100" , 200,  0.0,    0.1  );


  CreateUserTH1D("CorrIsolation_1stEle_ROI"                 , 200,-25.0 ,  25.0  ); CreateUserTH1D("CorrIsolation_2ndEle_ROI"                 , 200,-25.0 ,  25.0  );
  CreateUserTH1D("DeltaEtaTrkSC_1stEle_ROI"                 , 200, -0.01,   0.01 ); CreateUserTH1D("DeltaEtaTrkSC_2ndEle_ROI"                 , 200, -0.01,   0.01 );
  CreateUserTH1D("E2x5OverE5x5_1stEle_ROI"                  , 200,  0.0 ,   2.0  ); CreateUserTH1D("E2x5OverE5x5_2ndEle_ROI"                  , 200,  0.0 ,   2.0  );
  CreateUserTH1D("EcalIsolation_1stEle_ROI"                 , 200,  0.0 ,  20.0  ); CreateUserTH1D("EcalIsolation_2ndEle_ROI"                 , 200,  0.0 ,  20.0  );
  CreateUserTH1D("HcalIsolation_1stEle_ROI"                 , 200,  0.0 ,  20.0  ); CreateUserTH1D("HcalIsolation_2ndEle_ROI"                 , 200,  0.0 ,  20.0  );
  CreateUserTH1D("TrkIsolation_1stEle_ROI"                  , 200,  0.0,    5.0  ); CreateUserTH1D("TrkIsolation_2ndEle_ROI"                  , 200,  0.0,    5.0  );
  CreateUserTH1D("HasMatchedPhot_1stEle_ROI"                , 2,   -0.5 ,   1.5  ); CreateUserTH1D("HasMatchedPhot_2ndEle_ROI"                , 2,   -0.5 ,   1.5  );
  CreateUserTH1D("HoE_1stEle_ROI"                           , 200,  0.0 ,   0.05 ); CreateUserTH1D("HoE_2ndEle_ROI"                           , 200,  0.0 ,   0.05 );
  CreateUserTH1D("LeadVtxDistXY_1stEle_ROI"                 , 200, -0.05,   0.05 ); CreateUserTH1D("LeadVtxDistXY_2ndEle_ROI"                 , 200, -0.05,   0.05 );
  CreateUserTH1D("LeadVtxDistZ_1stEle_ROI"                  , 200, -0.2 ,   0.2  ); CreateUserTH1D("LeadVtxDistZ_2ndEle_ROI"                  , 200, -0.2 ,   0.2  );
  CreateUserTH1D("MissingHits_1stEle_ROI"                   , 2  , -0.5,    1.5  ); CreateUserTH1D("MissingHits_2ndEle_ROI"                   , 2  , -0.5,    1.5  );
  CreateUserTH1D("SigmaIEtaIEta_Barrel_1stEle_ROI"          , 200,  0.0,    0.02 ); CreateUserTH1D("SigmaIEtaIEta_Barrel_2ndEle_ROI"          , 200,  0.0,    0.02 );
  CreateUserTH1D("SigmaIEtaIEta_Endcap_1stEle_ROI"          , 200,  0.0,    0.1  ); CreateUserTH1D("SigmaIEtaIEta_Endcap_2ndEle_ROI"          , 200,  0.0,    0.1  );
  // for scale factor dependence studies
  CreateUserTH1D( "Mee_NJetEq2_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_NJetEq3_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_NJetEq4_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_NJetEq5_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_NJetEq6_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_NJetEq7_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_NJetGeq3_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_NJetGeq4_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_sT300To500_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_sT500To750_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_sT750To1250_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_sT1250ToInf_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_MejMin100To200_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_MejMin200To300_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_MejMin300To400_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_MejMin400To500_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_MejMin500To650_PAS"		             ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "Mee_MejMin650ToInf_PAS"		             ,    200   , 0       , 2000	  ); 

  //------------------------------------------------------------------
  // How many events to skim over?
  //------------------------------------------------------------------
  Long64_t nentries = GetTreeEntries();
  std::cout << "analysisClass::analysisClass(): nentries = " << nentries << std::endl;
  

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    readerTools_->LoadEntry(jentry);
    //--------------------------------------------------------------------------
    // Tricky part: refine Weight branch for data only
    //--------------------------------------------------------------------------
    float weight = 1.0;

    if(isData())
      resetSkimTreeBranchAddress("Weight", &weight);
    //------------------------------------------------------------------
    // Tell user how many events we've looped over
    //------------------------------------------------------------------
    if(jentry < 10 || jentry%5000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << "/" << nentries << std::endl;   
    
    //--------------------------------------------------------------------------
    // Reset the cuts
    //--------------------------------------------------------------------------

    resetCuts();

    //--------------------------------------------------------------------------
    // Check good run list
    //--------------------------------------------------------------------------

    float run = readerTools_->ReadValueBranch<Float_t>("run");
    float ls = readerTools_->ReadValueBranch<Float_t>("ls");
    float event = readerTools_->ReadValueBranch<Float_t>("event");
    int passedJSON = passJSON ( run, ls, isData() ) ;

    //--------------------------------------------------------------------------
    // Find the right prescale for this event
    //--------------------------------------------------------------------------

    double min_prescale = 1;
    int passTrigger  = 0;

    float Ele1_hltPhotonPt = readerTools_->ReadValueBranch<Float_t>("Ele1_hltPhotonPt");
    
    std::string triggerName = "";

    if ( Ele1_hltPhotonPt > 0.0 ) {
      if(analysisYear==2016) {
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon22")   > 0.1 && Ele1_hltPhotonPt >= 22.  && Ele1_hltPhotonPt < 30. ) { passTrigger = 1; triggerName = "Photon22"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon30")   > 0.1 && Ele1_hltPhotonPt >= 30.  && Ele1_hltPhotonPt < 36. ) { passTrigger = 1; triggerName = "Photon30"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon36")   > 0.1 && Ele1_hltPhotonPt >= 36.  && Ele1_hltPhotonPt < 50. ) { passTrigger = 1; triggerName = "Photon36"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon50")   > 0.1 && Ele1_hltPhotonPt >= 50.  && Ele1_hltPhotonPt < 75. ) { passTrigger = 1; triggerName = "Photon50"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon75")   > 0.1 && Ele1_hltPhotonPt >= 75.  && Ele1_hltPhotonPt < 90. ) { passTrigger = 1; triggerName = "Photon75"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon90")   > 0.1 && Ele1_hltPhotonPt >= 90.  && Ele1_hltPhotonPt < 120.) { passTrigger = 1; triggerName = "Photon90"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon120")  > 0.1 && Ele1_hltPhotonPt >= 120. && Ele1_hltPhotonPt < 175.) { passTrigger = 1; triggerName = "Photon120"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon175")  > 0.1 && Ele1_hltPhotonPt >= 175.) { passTrigger = 1; triggerName = "Photon175"; } 
      }
      else if(analysisYear==2017) {
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon25")   > 0.1 && Ele1_hltPhotonPt >= 25.  && Ele1_hltPhotonPt < 33. ) { passTrigger = 1; triggerName = "Photon25"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon33")   > 0.1 && Ele1_hltPhotonPt >= 33.  && Ele1_hltPhotonPt < 50. ) { passTrigger = 1; triggerName = "Photon33"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon50")   > 0.1 && Ele1_hltPhotonPt >= 50.  && Ele1_hltPhotonPt < 75. ) { passTrigger = 1; triggerName = "Photon50"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon75")   > 0.1 && Ele1_hltPhotonPt >= 75.  && Ele1_hltPhotonPt < 90. ) { passTrigger = 1; triggerName = "Photon75"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon90")   > 0.1 && Ele1_hltPhotonPt >= 90.  && Ele1_hltPhotonPt < 120.) { passTrigger = 1; triggerName = "Photon90"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon120")  > 0.1 && Ele1_hltPhotonPt >= 120. && Ele1_hltPhotonPt < 150.) { passTrigger = 1; triggerName = "Photon120"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon150")  > 0.1 && Ele1_hltPhotonPt >= 150. && Ele1_hltPhotonPt < 175.) { passTrigger = 1; triggerName = "Photon150"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon175")  > 0.1 && Ele1_hltPhotonPt >= 175. && Ele1_hltPhotonPt < 200.) { passTrigger = 1; triggerName = "Photon175"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon200")  > 0.1 && Ele1_hltPhotonPt >= 200.) { passTrigger = 1; triggerName = "Photon200"; } 
      }
      else if(analysisYear==2018) {
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon33")   > 0.1 && Ele1_hltPhotonPt >= 33.  && Ele1_hltPhotonPt < 50. ) { passTrigger = 1; triggerName = "Photon33"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon50")   > 0.1 && Ele1_hltPhotonPt >= 50.  && Ele1_hltPhotonPt < 75. ) { passTrigger = 1; triggerName = "Photon50"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon75")   > 0.1 && Ele1_hltPhotonPt >= 75.  && Ele1_hltPhotonPt < 90. ) { passTrigger = 1; triggerName = "Photon75"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon90")   > 0.1 && Ele1_hltPhotonPt >= 90.  && Ele1_hltPhotonPt < 120.) { passTrigger = 1; triggerName = "Photon90"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon120")  > 0.1 && Ele1_hltPhotonPt >= 120. && Ele1_hltPhotonPt < 150.) { passTrigger = 1; triggerName = "Photon120"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon150")  > 0.1 && Ele1_hltPhotonPt >= 150. && Ele1_hltPhotonPt < 175.) { passTrigger = 1; triggerName = "Photon150"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon175")  > 0.1 && Ele1_hltPhotonPt >= 175. && Ele1_hltPhotonPt < 200.) { passTrigger = 1; triggerName = "Photon175"; } 
        if ( readerTools_->ReadValueBranch<Float_t>("H_Photon200")  > 0.1 && Ele1_hltPhotonPt >= 200.) { passTrigger = 1; triggerName = "Photon200"; } 
      }
    }
    if(isData() && passTrigger) {
      //std::cout << "INFO: lookup trigger name " << triggerName << " for year: " << year << std::endl;
      min_prescale = run2PhotonTriggerPrescales.LookupPrescale(analysisYear,triggerName);
    }
    //else {
    //  std::cout << "Photon with hltPt = " << Ele1_hltPhotonPt << " did not pass any of the single photon triggers.." << std::endl;
    //}
    ///XXX SIC TEST FIXME to remove trigger and prescale requirements
    //passTrigger=1;
    //if ( H_Photon22   > 0.1 ) { passTrigger = 1; min_prescale = H_Photon22  ; } 
    //if ( H_Photon30   > 0.1 ) { passTrigger = 1; min_prescale = H_Photon30  ; } 
    //if ( H_Photon36   > 0.1 ) { passTrigger = 1; min_prescale = H_Photon36  ; } 
    //if ( H_Photon50   > 0.1 ) { passTrigger = 1; min_prescale = H_Photon50  ; } 
    //if ( H_Photon75   > 0.1 ) { passTrigger = 1; min_prescale = H_Photon75  ; } 
    //if ( H_Photon90   > 0.1 ) { passTrigger = 1; min_prescale = H_Photon90  ; } 
    //if ( H_Photon120  > 0.1 ) { passTrigger = 1; min_prescale = H_Photon120 ; } 
    //if ( H_Photon175  > 0.1 ) { passTrigger = 1; min_prescale = H_Photon175 ; } 
    //min_prescale=1.0;
    //// test output
    //std::cout << "Ele1_hltPhotonPt: " << Ele1_hltPhotonPt << std::endl;
    //std::cout << "Ele2_hltPhotonPt: " << readerTools_->ReadValueBranch<Float_t>("Ele2_hltPhotonPt") << std::endl;
    //std::cout << "PassTrigger? " << passTrigger << "; prescale of H_Photon22 = " << readerTools_->ReadValueBranch<Float_t>("H_Photon22")  << "; min_prescale=" << min_prescale << std::endl;
    //std::cout << "PassTrigger? " << passTrigger << "; prescale of H_Photon30 = " << readerTools_->ReadValueBranch<Float_t>("H_Photon30")  << "; min_prescale=" << min_prescale << std::endl;
    //std::cout << "PassTrigger? " << passTrigger << "; prescale of H_Photon36 = " << readerTools_->ReadValueBranch<Float_t>("H_Photon36")  << "; min_prescale=" << min_prescale << std::endl;
    //std::cout << "PassTrigger? " << passTrigger << "; prescale of H_Photon50 = " << readerTools_->ReadValueBranch<Float_t>("H_Photon50")  << "; min_prescale=" << min_prescale << std::endl;
    //std::cout << "PassTrigger? " << passTrigger << "; prescale of H_Photon75 = " << readerTools_->ReadValueBranch<Float_t>("H_Photon75")  << "; min_prescale=" << min_prescale << std::endl;
    //std::cout << "PassTrigger? " << passTrigger << "; prescale of H_Photon90 = " << readerTools_->ReadValueBranch<Float_t>("H_Photon90")  << "; min_prescale=" << min_prescale << std::endl;
    //std::cout << "PassTrigger? " << passTrigger << "; prescale of H_Photon120= " << readerTools_->ReadValueBranch<Float_t>("H_Photon120") << "; min_prescale=" << min_prescale << std::endl;
    //std::cout << "PassTrigger? " << passTrigger << "; prescale of H_Photon175= " << readerTools_->ReadValueBranch<Float_t>("H_Photon175") << "; min_prescale=" << min_prescale << std::endl;
    //// test output

    //--------------------------------------------------------------------------
    // What kind of event is this?
    //   - Barrel
    //   - Endcap 1 (eta < 2.0)
    //   - Endcap 2 (eta > 2.0) 
    //--------------------------------------------------------------------------

    bool ele1_isBarrel  = false;
    bool ele1_isEndcap1 = false;
    bool ele1_isEndcap2 = false;
    bool ele2_isBarrel  = false;
    bool ele2_isEndcap1 = false;
    bool ele2_isEndcap2 = false;

    float Ele1_SCEta = readerTools_->ReadValueBranch<Float_t>("Ele1_SCEta");
    float Ele2_SCEta = readerTools_->ReadValueBranch<Float_t>("Ele2_SCEta");
    //float Ele1_PtHeep = readerTools_->ReadValueBranch<Float_t>("Ele1_PtHeep");
    //float Ele2_PtHeep = readerTools_->ReadValueBranch<Float_t>("Ele2_PtHeep");
    float Ele1_Pt = readerTools_->ReadValueBranch<Float_t>("Ele1_Pt");
    float Ele2_Pt = readerTools_->ReadValueBranch<Float_t>("Ele2_Pt");
    float Ele1_Phi = readerTools_->ReadValueBranch<Float_t>("Ele1_Phi");
    float Ele2_Phi = readerTools_->ReadValueBranch<Float_t>("Ele2_Phi");

    if( fabs( Ele1_SCEta  ) < eleEta_bar )        ele1_isBarrel  = true;
    if( fabs( Ele1_SCEta  ) > eleEta_end1_min &&
        fabs( Ele1_SCEta  ) < eleEta_end1_max )   ele1_isEndcap1 = true;
    if( fabs( Ele1_SCEta  ) > eleEta_end2_min &&
        fabs( Ele1_SCEta  ) < eleEta_end2_max )   ele1_isEndcap2 = true;

    if( fabs( Ele2_SCEta  ) < eleEta_bar )        ele2_isBarrel  = true;
    if( fabs( Ele2_SCEta  ) > eleEta_end1_min &&
        fabs( Ele2_SCEta  ) < eleEta_end1_max )   ele2_isEndcap1 = true;
    if( fabs( Ele2_SCEta  ) > eleEta_end2_min &&
        fabs( Ele2_SCEta  ) < eleEta_end2_max )   ele2_isEndcap2 = true;

    bool ele1_isEndcap = ( ele1_isEndcap1 || ele1_isEndcap2 ) ;
    bool ele2_isEndcap = ( ele2_isEndcap1 || ele2_isEndcap2 ) ;

    bool isEBEB = ( ele1_isBarrel && ele2_isBarrel ) ;
    bool isEBEE = ( ( ele1_isBarrel && ele2_isEndcap ) ||
        ( ele2_isBarrel && ele1_isEndcap ) );
    bool isEEEE = ( ele1_isEndcap && ele2_isEndcap ) ;
    bool isEB   = ( isEBEB || isEBEE ) ;

    //--------------------------------------------------------------------------
    // Make this a QCD fake rate calculation
    //--------------------------------------------------------------------------
    float fakeRate1 = -1;
    float fakeRate2 = -1;
    if(analysisYear < 2018) {
      fakeRate1 = qcdFR.GetFakeRate(Ele1_Pt, QCDFakeRate::GetFakeRateRegion(ele1_isBarrel, ele1_isEndcap1, ele1_isEndcap2));
      fakeRate2 = qcdFR.GetFakeRate(Ele2_Pt, QCDFakeRate::GetFakeRateRegion(ele2_isBarrel, ele2_isEndcap1, ele2_isEndcap2));
    }
    else {
      fakeRate1 = qcdFR.GetFakeRate(Ele1_Pt, QCDFakeRate::GetFakeRateRegion(ele1_isBarrel, ele1_isEndcap1, ele1_isEndcap2,
            Ele1_SCEta, Ele1_Phi, run));
      fakeRate2 = qcdFR.GetFakeRate(Ele2_Pt, QCDFakeRate::GetFakeRateRegion(ele2_isBarrel, ele2_isEndcap1, ele2_isEndcap2,
            Ele2_SCEta, Ele2_Phi, run ));
    }

    //--------------------------------------------------------------------------
    // Finally have the effective fake rate
    //--------------------------------------------------------------------------

    //FIXME: add error on fake rate as well
    //double fakeRateEffective  = fakeRate1 * fakeRate2;
    double fakeRateEffective  = fakeRate1/(1-fakeRate1); // require loose electron to fail HEEP ID
    //if(1-fakeRate1 <= 0)
    //{
    //  cout << "ERROR: Found fakeRate1: " << fakeRate1 << " for SCEta=" << Ele1_SCEta << " SCEt="
    //    << Ele1_SCEnergy/cosh(Ele1_SCEta) << "=" << Ele1_SCEnergy << "/" << 
    //    cosh(Ele1_SCEta) << endl;
    //}
    float nEle_store = readerTools_->ReadValueBranch<Float_t>("nEle_store");
    if ( nEle_store >= 2 ) { 							        
      fakeRateEffective *= fakeRate2/(1-fakeRate2);
    }
    // double eFakeRateEffective = fakeRateEffective * sqrt (  ( eFakeRate1 / fakeRate1 ) * ( eFakeRate1 / fakeRate1 ) +
    //					     ( eFakeRate2 / fakeRate2 ) * ( eFakeRate2 / fakeRate2 ) );
    double eFakeRateEffective = 0.0;

    weight = fakeRateEffective * min_prescale;

    //--------------------------------------------------------------------------
    // Fill variables
    //--------------------------------------------------------------------------
    FillUserTH1D("EventCount"           , 1                   , 1   ); 
    //if(fakeRateEffective * min_prescale == 0.0)
    //  std::cout << "!!!!!THIS EVENT HAD fakeRateEffective * min_prescale == 0.0: " << fakeRateEffective * min_prescale << std::endl;

    // reweighting
    fillVariableWithValue ( "Reweighting", 1, fakeRateEffective * min_prescale ) ; 

    // JSON variable
    fillVariableWithValue(   "PassJSON"                      , passedJSON              , fakeRateEffective * min_prescale ) ; 

    // Fill noise filters
    // see: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    // we filled these at skim time
    fillVariableWithValue("PassGlobalSuperTightHalo2016Filter" , int(readerTools_->ReadValueBranch<Float_t>("PassGlobalSuperTightHalo2016Filter")     == 1), fakeRateEffective * min_prescale);
    fillVariableWithValue("PassGoodVertices"                   , int(readerTools_->ReadValueBranch<Float_t>("PassGoodVertices")                       == 1), fakeRateEffective * min_prescale);
    fillVariableWithValue("PassHBHENoiseFilter"                , int(readerTools_->ReadValueBranch<Float_t>("PassHBHENoiseFilter")                    == 1), fakeRateEffective * min_prescale);
    fillVariableWithValue("PassHBHENoiseIsoFilter"             , int(readerTools_->ReadValueBranch<Float_t>("PassHBHENoiseIsoFilter")                 == 1), fakeRateEffective * min_prescale);
    // eBadScFilter not suggested for MC
    if(isData())
      fillVariableWithValue("PassBadEESupercrystalFilter"      , int(readerTools_->ReadValueBranch<Float_t>("PassBadEESupercrystalFilter")            == 1), fakeRateEffective * min_prescale);
    else
      fillVariableWithValue("PassBadEESupercrystalFilter"      , 1                                                                                          , fakeRateEffective * min_prescale);
    fillVariableWithValue("PassEcalDeadCellTrigPrim"           , int(readerTools_->ReadValueBranch<Float_t>("PassEcalDeadCellTrigPrim")               == 1), fakeRateEffective * min_prescale);
    // not recommended
    //fillVariableWithValue("PassChargedCandidateFilter"         , int(readerTools_->ReadValueBranch<Float_t>("PassChargedCandidateFilter")             == 1), fakeRateEffective * min_prescale);
    fillVariableWithValue("PassBadPFMuonFilter"                , int(readerTools_->ReadValueBranch<Float_t>("PassBadPFMuonFilter")                    == 1), fakeRateEffective * min_prescale);
    // EcalBadCalibV2 for 2017, 2018
    if(analysisYear > 2016)
      fillVariableWithValue("PassEcalBadCalibV2Filter"         , int(readerTools_->ReadValueBranch<Float_t>("PassEcalBadCalibV2Filter")               == 1), fakeRateEffective * min_prescale);
    else
      fillVariableWithValue("PassEcalBadCalibV2Filter"         , 1                                                                                          , fakeRateEffective * min_prescale);

    // Electrons
    int PassNEle = 0;
    float nEle_ptCut = readerTools_->ReadValueBranch<Float_t>("nEle_ptCut");
    if ( nEle_ptCut >= 2 ) PassNEle = 1;
    // we only look at events that have exactly two loose electrons (passing Pt>10)

    //FIXME fix the counting to avoid the below
    //// nEle_store is the number of loose eles with pt>10
    //if(nEle_store==3)
    //{
    //  // if 3 stored, 2 leads must pass Pt cut and third must fail
    //  if(Ele3_Pt<50 && Ele1_Pt>=50 && Ele2_Pt>=50)
    //    PassNEle=1;
    //}
    //else if(nEle_store==2)
    //{
    //  // if 2 stored, OK
    //  if(Ele1_Pt>=50 && Ele2_Pt>=50)
    //      PassNEle=1;
    //}

    double M_ej_avg;
    double M_ej_min;
    double M_ej_max;

    // Muons
    int PassNMuon = 0;
    float nMuon_ptCut = readerTools_->ReadValueBranch<Float_t>("nMuon_ptCut");
    if ( nMuon_ptCut == 0 ) PassNMuon = 1;

    fillVariableWithValue ( "PassHLT"                        , passTrigger             , fakeRateEffective * min_prescale ) ;
    fillVariableWithValue("nEle_hltMatched",-1, fakeRateEffective * min_prescale ) ;
    fillVariableWithValue("nJet_hltMatched",-1, fakeRateEffective * min_prescale ) ;


    // Electrons								        
    //if(nEle_ptCut > 1) std::cout << "INFO: nEle_ptCut=" << nEle_ptCut << "; PassNEle=" << PassNEle << std::endl;
    fillVariableWithValue(   "PassNEle"                      , PassNEle                , fakeRateEffective * min_prescale ) ;
    float Ele1_Eta = readerTools_->ReadValueBranch<Float_t>("Ele1_Eta");
    float Ele1_PassHEEPID = readerTools_->ReadValueBranch<Float_t>("Ele1_PassHEEPID");
    float Ele2_Eta = readerTools_->ReadValueBranch<Float_t>("Ele2_Eta");
    float Ele2_PassHEEPID = readerTools_->ReadValueBranch<Float_t>("Ele2_PassHEEPID");
    float Pt_e1e2 = readerTools_->ReadValueBranch<Float_t>("Pt_e1e2");
    float M_e1e2 = readerTools_->ReadValueBranch<Float_t>("M_e1e2");
    if ( nEle_store >= 1 ) { 							        
      fillVariableWithValue( "Ele1_Eta"                      , Ele1_Eta            , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "Ele1_Pt"                       , Ele1_Pt            , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "Ele1_PassHEEPID"               , Ele1_PassHEEPID                , fakeRateEffective * min_prescale ) ;
    }										        
    if ( nEle_store >= 2 ) { 							        
      fillVariableWithValue( "Ele2_Eta"                        , Ele2_Eta            , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "Ele2_Pt"                       , Ele2_Pt                 , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "Ele2_PassHEEPID"                 , Ele2_PassHEEPID                , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "Pt_e1e2"                       , Pt_e1e2             , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "M_e1e2"                        , M_e1e2                  , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "M_e1e2_opt"                    , M_e1e2                  , fakeRateEffective * min_prescale ) ;
    }

    // Jets
    float nJet_ptCut = readerTools_->ReadValueBranch<Float_t>("nJet_ptCut");
    float nJet_store = readerTools_->ReadValueBranch<Float_t>("nJet_store");
    float Jet1_Pt = readerTools_->ReadValueBranch<Float_t>("Jet1_Pt");
    float Jet2_Pt = readerTools_->ReadValueBranch<Float_t>("Jet2_Pt");
    float Jet1_Eta = readerTools_->ReadValueBranch<Float_t>("Jet1_Eta");
    float Jet2_Eta = readerTools_->ReadValueBranch<Float_t>("Jet2_Eta");
    float Jet1_Phi = readerTools_->ReadValueBranch<Float_t>("Jet1_Phi");
    float Jet2_Phi = readerTools_->ReadValueBranch<Float_t>("Jet2_Phi");
    float DR_Jet1Jet2 = readerTools_->ReadValueBranch<Float_t>("DR_Jet1Jet2");
    fillVariableWithValue(   "nJet"                          , nJet_ptCut      , fakeRateEffective * min_prescale ) ;
    if ( nJet_store >= 1 ) { 						                
      fillVariableWithValue( "Jet1_Pt"                       , Jet1_Pt         , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "Jet1_Eta"                      , Jet1_Eta        , fakeRateEffective * min_prescale ) ;
    }

    float M_e1j1 = readerTools_->ReadValueBranch<Float_t>("M_e1j1");
    float M_e1j2 = readerTools_->ReadValueBranch<Float_t>("M_e1j2");
    float M_e2j1 = readerTools_->ReadValueBranch<Float_t>("M_e2j1");
    float M_e2j2 = readerTools_->ReadValueBranch<Float_t>("M_e2j2");
    if ( nJet_store >= 2 ) { 
      fillVariableWithValue( "Jet2_Pt"                       , Jet2_Pt         , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "Jet2_Eta"                      , Jet2_Eta        , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "DR_Jet1Jet2"                   , DR_Jet1Jet2             , fakeRateEffective * min_prescale ) ;

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
    }

    double sT_zjj = Pt_e1e2 + Jet1_Pt + Jet2_Pt;

    // Muons
    fillVariableWithValue(   "PassNMuon"                     , PassNMuon               , fakeRateEffective * min_prescale ) ;

    // DeltaR
    float DR_Ele1Jet1 = readerTools_->ReadValueBranch<Float_t>("DR_Ele1Jet1");
    float DR_Ele1Jet2 = readerTools_->ReadValueBranch<Float_t>("DR_Ele1Jet2");
    float DR_Ele2Jet1 = readerTools_->ReadValueBranch<Float_t>("DR_Ele2Jet1");
    float DR_Ele2Jet2 = readerTools_->ReadValueBranch<Float_t>("DR_Ele2Jet2");
    if ( nEle_store >= 2 && nJet_store >= 1) {
      fillVariableWithValue( "DR_Ele1Jet1"                   , DR_Ele1Jet1             , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "DR_Ele2Jet1"                   , DR_Ele2Jet1             , fakeRateEffective * min_prescale ) ;
      if(nJet_store >= 2) {
        fillVariableWithValue( "DR_Ele1Jet2"                 , DR_Ele1Jet2             , fakeRateEffective * min_prescale ) ;
        fillVariableWithValue( "DR_Ele2Jet2"                 , DR_Ele2Jet2             , fakeRateEffective * min_prescale ) ;
      }
    }

    // sT
    float sT_eejj = readerTools_->ReadValueBranch<Float_t>("sT_eejj");
    if ( nEle_store >= 2 && nJet_store >= 2) {
      // SIC recompute sT using PtHeep. FIXME: this is now being done in skims
      //sT_eejj = Ele1_PtHeep+Ele2_PtHeep+Jet1_Pt+Jet2_Pt;
      fillVariableWithValue( "sT_eejj"                       , sT_eejj                 , fakeRateEffective  * min_prescale ) ;
      fillVariableWithValue( "sT_eejj_opt"                   , sT_eejj                 , fakeRateEffective  * min_prescale ) ;
      fillVariableWithValue( "Mej_min_opt"                   , M_ej_min                , fakeRateEffective  * min_prescale ) ;
    }      

    // Dummy variables
    fillVariableWithValue ("preselection", 1, fakeRateEffective * min_prescale ); 

    //--------------------------------------------------------------------------
    // Evaluate the cuts
    //--------------------------------------------------------------------------

    evaluateCuts();
    //if(passedCut("PassNEle")) std::cout << "DID MANAGE TO PassNEle !" <<
    //  "The weight of this event: fakeRateEffective=" << fakeRateEffective << " x min_prescale=" <<
    //    min_prescale << " = " << fakeRateEffective*min_prescale << std::endl;

    //--------------------------------------------------------------------------
    // Did we at least pass the noise filtes?
    //--------------------------------------------------------------------------

    bool passed_minimum = ( passedAllPreviousCuts("PassBadPFMuonFilter") && passedCut ("PassBadPFMuonFilter"));

    //--------------------------------------------------------------------------
    // Did we pass preselection?
    //--------------------------------------------------------------------------

    bool passed_preselection = ( passedAllPreviousCuts("M_e1e2") && passedCut ("M_e1e2") );

    //--------------------------------------------------------------------------
    // Are we in the region of interest?
    //--------------------------------------------------------------------------

    bool passed_region_of_interest = bool ( passed_preselection && M_e1e2 > 170. && sT_eejj > 900.0 );

    //--------------------------------------------------------------------------
    // Fill plots with no selection applied
    //--------------------------------------------------------------------------

    FillUserTH1D( "PileupWeight"   , -1.0 );
    FillUserTH1D( "GeneratorWeight", -1.0 );

    //--------------------------------------------------------------------------
    // Fill preselection plots
    //--------------------------------------------------------------------------

    if ( passed_preselection ) {
      //FIXME DEBUG
      //std::cout << "The weight of this event: fakeRateEffective=" << fakeRateEffective << " x min_prescale=" << min_prescale << " = " << fakeRateEffective*min_prescale << std::endl;
      //
      //std::cout << "We have electrons which look like (ptHEEP,eta,phi): (" << Ele1_PtHeep << "," << Ele1_Eta << "," << Ele1_Phi << "); FR = " << qcdFakeRate.GetFakeRate(Ele1_Eta,Ele2_PtHeep) << std::endl;
      //std::cout << "We have electrons which look like (pt,eta,phi): (" << Ele1_Pt << "," << Ele1_Eta << "," << Ele1_Phi << "); FR = " << qcdFakeRate.GetFakeRate(Ele1_Eta,Ele2_Pt) << std::endl;
      //std::cout << "Used fake rate=" << fakeRate1 << std::endl;
      ////float fakeRate2 = qcdFakeRate.GetFakeRate(Ele2_Eta,Ele2_PtHeep);
      //std::cout << "We have electrons which look like (pt,eta,phi): (" << Ele2_Pt << "," << Ele2_Eta << "," << Ele2_Phi << ")" << std::endl;
      //if(nEle_store > 2)
      //  std::cout << "We have electrons which look like (pt,eta,phi): (" << Ele3_Pt << "," << Ele3_Eta << "," << Ele3_Phi << ")" << std::endl;

      //--------------------------------------------------------------------------
      // Recalculate some variables
      //--------------------------------------------------------------------------
      float Ele1_Pt = readerTools_->ReadValueBranch<Float_t>("Ele1_Pt");
      float Ele2_Pt = readerTools_->ReadValueBranch<Float_t>("Ele2_Pt");
      float Ele1_Phi = readerTools_->ReadValueBranch<Float_t>("Ele1_Phi");
      float Ele2_Phi = readerTools_->ReadValueBranch<Float_t>("Ele2_Phi");
      float Muon1_Pt = readerTools_->ReadValueBranch<Float_t>("Muon1_Pt");
      float Muon2_Pt = readerTools_->ReadValueBranch<Float_t>("Muon2_Pt");
      float Muon1_Eta = readerTools_->ReadValueBranch<Float_t>("Muon1_Eta");
      float Muon2_Eta = readerTools_->ReadValueBranch<Float_t>("Muon2_Eta");
      float Muon1_Phi = readerTools_->ReadValueBranch<Float_t>("Muon1_Phi");
      float Muon2_Phi = readerTools_->ReadValueBranch<Float_t>("Muon2_Phi");
      float PFMET_Type1_Pt = readerTools_->ReadValueBranch<Float_t>("PFMET_Type1_Pt");
      float PFMET_Type1_Phi = readerTools_->ReadValueBranch<Float_t>("PFMET_Type1_Phi");

      TLorentzVector e1, j1, e2, j2,j3, mu, met;
      TLorentzVector eejj, e1e2mu;
      TLorentzVector eej, ejj, ee;
      TLorentzVector e1j3, e2j3, j1j3, j2j3, j1j2, j1j2j3, eejjj;

      e1.SetPtEtaPhiM ( Ele1_Pt, Ele1_SCEta, Ele1_Phi, 0.0 );
      e2.SetPtEtaPhiM ( Ele2_Pt, Ele2_SCEta, Ele2_Phi, 0.0 );
      j1.SetPtEtaPhiM ( Jet1_Pt, Jet1_Eta, Jet1_Phi, 0.0 );
      j2.SetPtEtaPhiM ( Jet2_Pt, Jet2_Eta, Jet2_Phi, 0.0 );
      mu.SetPtEtaPhiM ( Muon1_Pt, Muon1_Eta, Muon1_Phi, 0.0 );
      met.SetPtEtaPhiM ( PFMET_Type1_Pt, 0.0, PFMET_Type1_Phi, 0.0 );

      eejj = e1 + e2 + j1 + j2 ; 
      eej  = e1 + e2 + j1;
      ejj  = e1 + j1 + j2;
      ee   = e1 + e2;
      j1j2 = j1 + j2;

      double min_DR_EleJet;
      double DR_Ele1Jet3;
      double DR_Ele2Jet3;
      double DR_Ele1Ele2 = e1.DeltaR( e2 );

      double M_eejj = eejj.M();
      double M_eej  = eej.M();
      double M_ejj  = ejj.M();

      double Pt_j1j2 = j1j2.Pt();
      double Pt_j1j3;
      double Pt_j2j3;
      double Pt_j1j2j3;

      double M_j1j3, M_j2j3;
      double M_e1j3, M_e2j3;
      double M_eejjj;

      double min_DeltaR_Zj = 999.;
      if ( ee.DeltaR ( j1 ) < min_DeltaR_Zj ) min_DeltaR_Zj = ee.DeltaR ( j1 );
      if ( ee.DeltaR ( j2 ) < min_DeltaR_Zj ) min_DeltaR_Zj = ee.DeltaR ( j2 );
      double DR_ZJ1 = ee.DeltaR ( j1 );
      double DR_ZJ2 = ee.DeltaR ( j2 );

      float Jet3_Pt = readerTools_->ReadValueBranch<Float_t>("Jet3_Pt");
      float Jet3_Eta = readerTools_->ReadValueBranch<Float_t>("Jet3_Eta");
      float Jet3_Phi = readerTools_->ReadValueBranch<Float_t>("Jet3_Phi");
      if ( nJet_ptCut > 2 ) { 
        j3.SetPtEtaPhiM ( Jet3_Pt, Jet3_Eta, Jet3_Phi, 0.0 );

        e1j3 = e1 + j3;
        e2j3 = e2 + j3;
        j1j3 = j1 + j3;
        j2j3 = j2 + j3;
        j1j2j3 = j1 + j2 + j3;
        eejjj = e1 + e2 + j1 + j2 + j3;

        DR_Ele1Jet3 = e1.DeltaR( j3 );
        DR_Ele2Jet3 = e2.DeltaR( j3 );
        M_e1j3 = e1j3.M();
        M_e2j3 = e2j3.M();
        M_j1j3 = j1j3.M();
        M_j2j3 = j2j3.M();
        Pt_j1j3 = j1j3.Pt();
        Pt_j2j3 = j2j3.Pt();
        Pt_j1j2j3 = j1j2j3.Pt();
        M_eejjj = eejjj.M();

      }

      if ( DR_Ele1Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet1;
      if ( DR_Ele1Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet2;
      if ( DR_Ele2Jet1 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet1;
      if ( DR_Ele2Jet2 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet2;
      if ( nJet_ptCut > 2 ) {
        if ( DR_Ele1Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet3;
        if ( DR_Ele2Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet3;
      }

      //--------------------------------------------------------------------------
      // Electron quality histograms (preselection)
      //--------------------------------------------------------------------------
       float Ele1_CorrIsolation        = readerTools_->ReadValueBranch<Float_t>("Ele1_CorrIsolation")      ; 
       float Ele1_DeltaEtaTrkSC        = readerTools_->ReadValueBranch<Float_t>("Ele1_DeltaEtaTrkSC")      ; 
       float Ele1_EcalIsolation        = readerTools_->ReadValueBranch<Float_t>("Ele1_EcalIsolation")      ; 
       float Ele1_HcalIsolation        = readerTools_->ReadValueBranch<Float_t>("Ele1_HcalIsolation")      ; 
       float Ele1_TrkIsolation         = readerTools_->ReadValueBranch<Float_t>("Ele1_TrkIsolation")       ; 
       float Ele1_HasMatchedPhot       = readerTools_->ReadValueBranch<Float_t>("Ele1_HasMatchedPhot")     ; 
       float Ele1_HoE                  = readerTools_->ReadValueBranch<Float_t>("Ele1_HoE")                ; 
       float Ele1_LeadVtxDistXY        = readerTools_->ReadValueBranch<Float_t>("Ele1_LeadVtxDistXY")      ; 
       float Ele1_LeadVtxDistZ         = readerTools_->ReadValueBranch<Float_t>("Ele1_LeadVtxDistZ")       ; 
       float Ele1_MissingHits          = readerTools_->ReadValueBranch<Float_t>("Ele1_MissingHits")        ; 
       float Ele1_SCEta                = readerTools_->ReadValueBranch<Float_t>("Ele1_SCEta");
       float Ele1_Full5x5SigmaIEtaIEta = readerTools_->ReadValueBranch<Float_t>("Ele1_Full5x5SigmaIEtaIEta");
       float Ele1_Charge               = readerTools_->ReadValueBranch<Float_t>("Ele1_Charge");
       //float Ele1_PtHeep               = readerTools_->ReadValueBranch<Float_t>("Ele1_PtHeep");

      FillUserTH1D("CorrIsolation_1stEle_PAS"         , Ele1_CorrIsolation                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("DeltaEtaTrkSC_1stEle_PAS"         , Ele1_DeltaEtaTrkSC                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("EcalIsolation_1stEle_PAS"         , Ele1_EcalIsolation                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("HcalIsolation_1stEle_PAS"         , Ele1_HcalIsolation                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("TrkIsolation_1stEle_PAS"          , Ele1_TrkIsolation                   , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("HasMatchedPhot_1stEle_PAS"        , Ele1_HasMatchedPhot                 , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("HoE_1stEle_PAS"                   , Ele1_HoE                            , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("LeadVtxDistXY_1stEle_PAS"         , Ele1_LeadVtxDistXY                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("LeadVtxDistZ_1stEle_PAS"          , Ele1_LeadVtxDistZ                   , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("MissingHits_1stEle_PAS"           , Ele1_MissingHits                    , min_prescale * fakeRateEffective   ); 
      if ( fabs(Ele1_SCEta) < eleEta_bar ) { 
        FillUserTH1D("SigmaIEtaIEta_Barrel_1stEle_PAS", Ele1_Full5x5SigmaIEtaIEta           , min_prescale * fakeRateEffective   ); 
      }
      else if ( fabs(Ele1_SCEta) > eleEta_end1_min && fabs(Ele2_SCEta) < eleEta_end2_max ){
        FillUserTH1D("SigmaIEtaIEta_Endcap_1stEle_PAS", Ele1_Full5x5SigmaIEtaIEta           , min_prescale * fakeRateEffective   ); 
      }

      float Ele2_CorrIsolation        = readerTools_->ReadValueBranch<Float_t>("Ele2_CorrIsolation")      ; 
      float Ele2_DeltaEtaTrkSC        = readerTools_->ReadValueBranch<Float_t>("Ele2_DeltaEtaTrkSC")      ; 
      float Ele2_EcalIsolation        = readerTools_->ReadValueBranch<Float_t>("Ele2_EcalIsolation")      ; 
      float Ele2_HcalIsolation        = readerTools_->ReadValueBranch<Float_t>("Ele2_HcalIsolation")      ; 
      float Ele2_TrkIsolation         = readerTools_->ReadValueBranch<Float_t>("Ele2_TrkIsolation")       ; 
      float Ele2_HasMatchedPhot       = readerTools_->ReadValueBranch<Float_t>("Ele2_HasMatchedPhot")     ; 
      float Ele2_HoE                  = readerTools_->ReadValueBranch<Float_t>("Ele2_HoE")                ; 
      float Ele2_LeadVtxDistXY        = readerTools_->ReadValueBranch<Float_t>("Ele2_LeadVtxDistXY")      ; 
      float Ele2_LeadVtxDistZ         = readerTools_->ReadValueBranch<Float_t>("Ele2_LeadVtxDistZ")       ; 
      float Ele2_MissingHits          = readerTools_->ReadValueBranch<Float_t>("Ele2_MissingHits")        ; 
      float Ele2_SCEta                = readerTools_->ReadValueBranch<Float_t>("Ele2_SCEta");
      float Ele2_Full5x5SigmaIEtaIEta = readerTools_->ReadValueBranch<Float_t>("Ele2_Full5x5SigmaIEtaIEta");
      float Ele2_Charge               = readerTools_->ReadValueBranch<Float_t>("Ele2_Charge");
      //float Ele2_PtHeep               = readerTools_->ReadValueBranch<Float_t>("Ele2_PtHeep");

      FillUserTH1D("CorrIsolation_2ndEle_PAS"         , Ele2_CorrIsolation                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("DeltaEtaTrkSC_2ndEle_PAS"         , Ele2_DeltaEtaTrkSC                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("EcalIsolation_2ndEle_PAS"         , Ele2_EcalIsolation                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("HcalIsolation_2ndEle_PAS"         , Ele2_HcalIsolation                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("TrkIsolation_2ndEle_PAS"          , Ele2_TrkIsolation                   , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("HasMatchedPhot_2ndEle_PAS"        , Ele2_HasMatchedPhot                 , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("HoE_2ndEle_PAS"                   , Ele2_HoE                            , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("LeadVtxDistXY_2ndEle_PAS"         , Ele2_LeadVtxDistXY                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("LeadVtxDistZ_2ndEle_PAS"          , Ele2_LeadVtxDistZ                   , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("MissingHits_2ndEle_PAS"           , Ele2_MissingHits                    , min_prescale * fakeRateEffective   ); 
      if ( fabs(Ele2_SCEta) < eleEta_bar ) { 
        FillUserTH1D("SigmaIEtaIEta_Barrel_2ndEle_PAS", Ele2_Full5x5SigmaIEtaIEta           , min_prescale * fakeRateEffective   ); 
      }
      else if ( fabs(Ele2_SCEta) > eleEta_end1_min && fabs(Ele2_SCEta) < eleEta_end2_max ){
        FillUserTH1D("SigmaIEtaIEta_Endcap_2ndEle_PAS", Ele2_Full5x5SigmaIEtaIEta           , min_prescale * fakeRateEffective   ); 
      }

      //--------------------------------------------------------------------------
      // Preselection histograms
      //--------------------------------------------------------------------------
      float M_j1j2 = readerTools_->ReadValueBranch<Float_t>("M_j1j2");
      float nVertex = readerTools_->ReadValueBranch<Float_t>("nVertex");

      FillUserTH1D( "Ptj1j2_PAS"           , Pt_j1j2                        , min_prescale * fakeRateEffective ) ;
      FillUserTH1D( "Ptee_Minus_Ptj1j2_PAS", Pt_e1e2 - Pt_j1j2              , min_prescale * fakeRateEffective ) ;

      FillUserTH1D("minDR_EleJet_PAS"     , min_DR_EleJet                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("DR_Ele1Ele2_PAS"	   , DR_Ele1Ele2                        , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("EleChargeSum_PAS"     , Ele1_Charge + Ele2_Charge, min_prescale * fakeRateEffective ) ;

      FillUserTH1D("nElectron_PAS"        , nEle_ptCut                    , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("nMuon_PAS"            , nMuon_ptCut                        , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("nJet_PAS"             , nJet_ptCut                 , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Pt1stEle_PAS"	   , Ele1_Pt                       , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Eta1stEle_PAS"	   , Ele1_SCEta                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Phi1stEle_PAS"	   , Ele1_Phi                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Pt2ndEle_PAS"	   , Ele2_Pt                       , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Eta2ndEle_PAS"	   , Ele2_SCEta                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Phi2ndEle_PAS"	   , Ele2_Phi                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Charge1stEle_PAS"	   , Ele1_Charge                   , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Charge2ndEle_PAS"	   , Ele2_Charge                   , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("MET_PAS"              , PFMET_Type1_Pt                  , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("METPhi_PAS"	   , PFMET_Type1_Phi                 , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Pt1stJet_PAS"         , Jet1_Pt                    , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Pt2ndJet_PAS"         , Jet2_Pt                    , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Eta1stJet_PAS"        , Jet1_Eta                   , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Eta2ndJet_PAS"        , Jet2_Eta                   , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Phi1stJet_PAS"	   , Jet1_Phi                   , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Phi2ndJet_PAS"	   , Jet2_Phi                   , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("sTlep_PAS"            , Ele1_Pt + Ele2_Pt        , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("sTjet_PAS"            , Jet1_Pt + Jet2_Pt  , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("sT_PAS"               , sT_eejj                            , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("sT_zjj_PAS"           , sT_zjj                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Mjj_PAS"		   , M_j1j2                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Mee_PAS"		   , M_e1e2                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("MTenu_PAS"            , readerTools_->ReadValueBranch<Float_t>("MT_Ele1MET")                   , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Me1j1_PAS"            , M_e1j1                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Me1j2_PAS"            , M_e1j2                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Me2j1_PAS"            , M_e2j1                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Me2j2_PAS"            , M_e2j2                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Ptee_PAS"             , Pt_e1e2                            , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("nVertex_PAS"           , nVertex                         , min_prescale * fakeRateEffective );
      FillUserTH1D("DR_Ele1Jet1_PAS"	   , DR_Ele1Jet1                        , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("DR_Ele1Jet2_PAS"	   , DR_Ele1Jet2                        , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("DR_Ele2Jet1_PAS"	   , DR_Ele2Jet1                        , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("DR_Ele2Jet2_PAS"	   , DR_Ele2Jet2                        , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("DR_Jet1Jet2_PAS"	   , DR_Jet1Jet2                        , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Meejj_PAS"            , M_eejj                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Meej_PAS"             , M_eej                              , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Mejj_PAS"             , M_ejj                              , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("minDR_ZJet_PAS"       , min_DeltaR_Zj                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("DR_ZJet1_PAS"         , DR_ZJ1                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("DR_ZJet2_PAS"         , DR_ZJ2                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Mej_selected_avg_PAS" , M_ej_avg                           , min_prescale * fakeRateEffective ) ;	   
      FillUserTH1D("Mej_selected_min_PAS" , M_ej_min                           , min_prescale * fakeRateEffective ) ;	   
      FillUserTH1D("Mej_selected_max_PAS" , M_ej_max                           , min_prescale * fakeRateEffective ) ;	   
      FillUserTH1D("Mej_minmax_PAS"       , M_ej_min                           , min_prescale * fakeRateEffective ) ;	   
      FillUserTH1D("Mej_minmax_PAS"       , M_ej_max                           , min_prescale * fakeRateEffective ) ;	   

      FillUserTH2D("MeeVsST_PAS" , M_e1e2, sT_eejj, min_prescale * fakeRateEffective ) ;	   
      // scale factor dependence histos
      if ( nJet_ptCut == 2 )
        FillUserTH1D("Mee_NJetEq2_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if( nJet_ptCut == 3 )
        FillUserTH1D("Mee_NJetEq3_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if( nJet_ptCut == 4 )
        FillUserTH1D("Mee_NJetEq4_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if( nJet_ptCut == 5 )
        FillUserTH1D("Mee_NJetEq5_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if( nJet_ptCut == 6 )
        FillUserTH1D("Mee_NJetEq6_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if( nJet_ptCut == 7 )
        FillUserTH1D("Mee_NJetEq7_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      //
      if ( nJet_ptCut >= 3 )
        FillUserTH1D("Mee_NJetGeq3_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      if ( nJet_ptCut >= 4 )
        FillUserTH1D("Mee_NJetGeq4_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      //
      if (sT_eejj >= 300 && sT_eejj < 500)
        FillUserTH1D("Mee_sT300To500_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if (sT_eejj >= 500 && sT_eejj < 750)
        FillUserTH1D("Mee_sT500To750_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if (sT_eejj >= 750 && sT_eejj < 1250)
        FillUserTH1D("Mee_sT750To1250_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if (sT_eejj >= 1250)
        FillUserTH1D("Mee_sT1250ToInf_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      //
      if (M_ej_min >= 100 && M_ej_min < 200)
        FillUserTH1D("Mee_MejMin100To200_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if (M_ej_min >= 200 && M_ej_min < 300)
        FillUserTH1D("Mee_MejMin200To300_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if (M_ej_min >= 300 && M_ej_min < 400)
        FillUserTH1D("Mee_MejMin300To400_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if (M_ej_min >= 400 && M_ej_min < 500)
        FillUserTH1D("Mee_MejMin400To500_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if (M_ej_min >= 500 && M_ej_min < 650)
        FillUserTH1D("Mee_MejMin500To650_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if (M_ej_min >= 650)
        FillUserTH1D("Mee_MejMin650ToInf_PAS", M_e1e2                         , min_prescale * fakeRateEffective );

      //--------------------------------------------------------------------------
      // Mass-pairing histograms at preselection
      //--------------------------------------------------------------------------

      if ( fabs(M_e1j1-M_e2j2) < fabs(M_e1j2-M_e2j1) )  {
        FillUserTH1D("Me1j_selected_PAS"   , M_e1j1,         min_prescale * fakeRateEffective );	   
        FillUserTH1D("Me2j_selected_PAS"   , M_e2j2,         min_prescale * fakeRateEffective );	   
        FillUserTH2D("Me1jVsMe2j_selected" , M_e1j1, M_e2j2, min_prescale * fakeRateEffective );
        FillUserTH2D("Me1jVsMe2j_rejected" , M_e1j2, M_e2j1, min_prescale * fakeRateEffective );
      }
      else {
        FillUserTH1D("Me1j_selected_PAS"   , M_e1j2,         min_prescale * fakeRateEffective );	   
        FillUserTH1D("Me2j_selected_PAS"   , M_e2j1,         min_prescale * fakeRateEffective );	   
        FillUserTH2D("Me1jVsMe2j_selected" , M_e1j2, M_e2j1, min_prescale * fakeRateEffective );
        FillUserTH2D("Me1jVsMe2j_rejected" , M_e1j1, M_e2j2, min_prescale * fakeRateEffective );
      }

      //--------------------------------------------------------------------------
      // Preselection + N(Jet) > 2 
      //--------------------------------------------------------------------------

      if ( nJet_ptCut > 2 ){ 
        FillUserTH1D( "M_e1j3_PAS"  , M_e1j3, min_prescale * fakeRateEffective  ); 
        FillUserTH1D( "M_e2j3_PAS"  , M_e2j3, min_prescale * fakeRateEffective  ); 
        FillUserTH1D( "M_j1j3_PAS"  , M_j1j3, min_prescale * fakeRateEffective  ); 
        FillUserTH1D( "M_j2j3_PAS"  , M_j2j3, min_prescale * fakeRateEffective  ); 
        FillUserTH1D( "M_eejjj_PAS" , M_eejjj,min_prescale * fakeRateEffective  ); 

        FillUserTH1D( "Ptj1j2j3_PAS"            , Pt_j1j2j3           , min_prescale * fakeRateEffective ) ;
        FillUserTH1D( "Ptj2j3_PAS"              , Pt_j2j3             , min_prescale * fakeRateEffective ) ;
        FillUserTH1D( "Ptj1j3_PAS"              , Pt_j1j3             , min_prescale * fakeRateEffective ) ;
        FillUserTH1D( "Ptee_Minus_Ptj1j2j3_PAS" , Pt_e1e2 - Pt_j1j2j3 , min_prescale * fakeRateEffective ) ;
      }

      //--------------------------------------------------------------------------
      // Preselection + event type (EBEB, EEEB, EEEE, etc)
      //--------------------------------------------------------------------------

      if      ( isEB   ) FillUserTH1D( "Mee_EB_PAS"  , M_e1e2, min_prescale * fakeRateEffective  ); 
      if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
      else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
      else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 

      //--------------------------------------------------------------------------
      // Preselection + high ST plot
      //--------------------------------------------------------------------------

      if ( sT_eejj > 445. ) FillUserTH1D( "Mee_PASandST445", M_e1e2, min_prescale * fakeRateEffective ) ;

      //--------------------------------------------------------------------------
      // High M(ee) plots
      //--------------------------------------------------------------------------

      if ( M_e1e2 > 100. ) FillUserTH1D("sT_PASandMee100"   , sT_eejj , min_prescale * fakeRateEffective  ); 
      if ( M_e1e2 > 110. ) FillUserTH1D("sT_PASandMee110"   , sT_eejj , min_prescale * fakeRateEffective  ); 
      if ( M_e1e2 > 120. ) FillUserTH1D("sT_PASandMee120"   , sT_eejj , min_prescale * fakeRateEffective  ); 
      if ( M_e1e2 > 130. ) FillUserTH1D("sT_PASandMee130"   , sT_eejj , min_prescale * fakeRateEffective  ); 
      if ( M_e1e2 > 140. ) FillUserTH1D("sT_PASandMee140"   , sT_eejj , min_prescale * fakeRateEffective  ); 
      if ( M_e1e2 > 150. ) FillUserTH1D("sT_PASandMee150"   , sT_eejj , min_prescale * fakeRateEffective  ); 
      if ( M_e1e2 > 160. ) FillUserTH1D("sT_PASandMee160"   , sT_eejj , min_prescale * fakeRateEffective  ); 
      if ( M_e1e2 > 170. ) FillUserTH1D("sT_PASandMee170"   , sT_eejj , min_prescale * fakeRateEffective  ); 
      if ( M_e1e2 > 180. ) FillUserTH1D("sT_PASandMee180"   , sT_eejj , min_prescale * fakeRateEffective  ); 
      if ( M_e1e2 > 190. ) FillUserTH1D("sT_PASandMee190"   , sT_eejj , min_prescale * fakeRateEffective  ); 
      if ( M_e1e2 > 200. ) FillUserTH1D("sT_PASandMee200"   , sT_eejj , min_prescale * fakeRateEffective  ); 


      if ( M_e1e2 > 100. ) { 

        FillUserTH1D("CorrIsolation_1stEle_PASandMee100"         , Ele1_CorrIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("DeltaEtaTrkSC_1stEle_PASandMee100"         , Ele1_DeltaEtaTrkSC                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("EcalIsolation_1stEle_PASandMee100"         , Ele1_EcalIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HcalIsolation_1stEle_PASandMee100"         , Ele1_HcalIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("TrkIsolation_1stEle_PASandMee100"          , Ele1_TrkIsolation                   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HasMatchedPhot_1stEle_PASandMee100"        , Ele1_HasMatchedPhot                 , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HoE_1stEle_PASandMee100"                   , Ele1_HoE                            , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("LeadVtxDistXY_1stEle_PASandMee100"         , Ele1_LeadVtxDistXY                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("LeadVtxDistZ_1stEle_PASandMee100"          , Ele1_LeadVtxDistZ                   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("MissingHits_1stEle_PASandMee100"           , Ele1_MissingHits                    , min_prescale * fakeRateEffective   ); 
        if ( fabs(Ele1_SCEta) < eleEta_bar ) { 
          FillUserTH1D("SigmaIEtaIEta_Barrel_1stEle_PASandMee100", Ele1_Full5x5SigmaIEtaIEta           , min_prescale * fakeRateEffective   ); 
        }
        else if ( fabs(Ele1_SCEta) > eleEta_end1_min && fabs(Ele2_SCEta) < eleEta_end2_max ){
          FillUserTH1D("SigmaIEtaIEta_Endcap_1stEle_PASandMee100", Ele1_Full5x5SigmaIEtaIEta           , min_prescale * fakeRateEffective   ); 
        }

        FillUserTH1D("CorrIsolation_2ndEle_PASandMee100"         , Ele2_CorrIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("DeltaEtaTrkSC_2ndEle_PASandMee100"         , Ele2_DeltaEtaTrkSC                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("EcalIsolation_2ndEle_PASandMee100"         , Ele2_EcalIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HcalIsolation_2ndEle_PASandMee100"         , Ele2_HcalIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("TrkIsolation_2ndEle_PASandMee100"          , Ele2_TrkIsolation                   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HasMatchedPhot_2ndEle_PASandMee100"        , Ele2_HasMatchedPhot                 , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HoE_2ndEle_PASandMee100"                   , Ele2_HoE                            , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("LeadVtxDistXY_2ndEle_PASandMee100"         , Ele2_LeadVtxDistXY                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("LeadVtxDistZ_2ndEle_PASandMee100"          , Ele2_LeadVtxDistZ                   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("MissingHits_2ndEle_PASandMee100"           , Ele2_MissingHits                    , min_prescale * fakeRateEffective   ); 
        if ( fabs(Ele2_SCEta) < eleEta_bar ) { 
          FillUserTH1D("SigmaIEtaIEta_Barrel_2ndEle_PASandMee100", Ele2_Full5x5SigmaIEtaIEta           , min_prescale * fakeRateEffective   ); 
        }
        else if ( fabs(Ele2_SCEta) > eleEta_end1_min && fabs(Ele2_SCEta) < eleEta_end2_max ){
          FillUserTH1D("SigmaIEtaIEta_Endcap_2ndEle_PASandMee100", Ele2_Full5x5SigmaIEtaIEta           , min_prescale * fakeRateEffective   ); 
        }

        FillUserTH1D("Me1j1_PASandMee100"           , M_e1j1                              , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Ptee_PASandMee100"            , Pt_e1e2                             , min_prescale * fakeRateEffective ) ;
        FillUserTH2D("MeeVsST_PASandMee100" , M_e1e2, sT_eejj, min_prescale * fakeRateEffective ) ;	   
        FillUserTH1D("sT_zjj_PASandMee100"          , sT_zjj                              , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("nVertex_PASandMee100"         , nVertex                             , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sT_PASandMee100"              , sT_eejj                             , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("EleChargeSum_PASandMee100"    , Ele1_Charge + Ele2_Charge , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("nJet_PASandMee100"            , nJet_ptCut                  , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sTlep_PASandMee100"           , Ele1_Pt    + Ele2_Pt      , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sTjet_PASandMee100"           , Jet1_Pt + Jet2_Pt   , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Mjj_PASandMee100"             , M_j1j2                              , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Pt1stEle_PASandMee100"        , Ele1_Pt                        , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Pt2ndEle_PASandMee100"        , Ele2_Pt                        , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Pt1stJet_PASandMee100"        , Jet1_Pt                     , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Pt2ndJet_PASandMee100"        , Jet2_Pt                     , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Mej_selected_avg_PASandMee100", M_ej_avg                            , min_prescale * fakeRateEffective ) ;

        FillUserTH1D("sTfrac_Jet1_PASandMee100"     , Jet1_Pt / sT_eejj                       , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sTfrac_Jet2_PASandMee100"     , Jet2_Pt / sT_eejj                       , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sTfrac_Ele1_PASandMee100"     , Ele1_Pt / sT_eejj                          , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sTfrac_Ele2_PASandMee100"     , Ele2_Pt / sT_eejj                          , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sTfrac_Jet_PASandMee100"      , ( Jet1_Pt + Jet2_Pt ) / sT_eejj , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sTfrac_Ele_PASandMee100"      , ( Ele1_Pt + Ele2_Pt ) / sT_eejj       , min_prescale * fakeRateEffective ) ;

        FillUserTH1D("Ptj1j2_PASandMee100"            , Pt_j1j2                        ,  min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Ptee_Minus_Ptj1j2_PASandMee100" , Pt_e1e2 - Pt_j1j2              ,  min_prescale * fakeRateEffective ) ;

        if ( nJet_ptCut > 2 ) { 	 
          FillUserTH1D( "M_j1j3_PASandMee100" , M_j1j3, min_prescale * fakeRateEffective ) ;
          FillUserTH1D( "M_j2j3_PASandMee100" , M_j2j3, min_prescale * fakeRateEffective ) ;
          FillUserTH1D( "M_e1j3_PASandMee100" , M_e1j3, min_prescale * fakeRateEffective ) ;
          FillUserTH1D( "M_e2j3_PASandMee100" , M_e2j3, min_prescale * fakeRateEffective ) ;
          FillUserTH1D( "M_eejjj_PASandMee100", M_eejjj,min_prescale * fakeRateEffective ) ;

          FillUserTH1D( "Ptj1j2j3_PASandMee100"            , Pt_j1j2j3           , min_prescale * fakeRateEffective ) ;
          FillUserTH1D( "Ptj2j3_PASandMee100"              , Pt_j2j3             , min_prescale * fakeRateEffective ) ;
          FillUserTH1D( "Ptj1j3_PASandMee100"              , Pt_j1j3             , min_prescale * fakeRateEffective ) ;
          FillUserTH1D( "Ptee_Minus_Ptj1j2j3_PASandMee100" , Pt_e1e2 - Pt_j1j2j3 , min_prescale * fakeRateEffective ) ;
        }
      }

      //--------------------------------------------------------------------------
      // Preselection + M(ee) normalization region plots
      //--------------------------------------------------------------------------

      if ( M_e1e2 > 80.0 && M_e1e2 < 100.0 ){
        FillUserTH1D("Mee_80_100_Preselection", M_e1e2, min_prescale * fakeRateEffective ) ;
        if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_80_100_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
        else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_80_100_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
        else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_80_100_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
        if      ( isEB   ) FillUserTH1D( "Mee_EB_80_100_PAS"  , M_e1e2, min_prescale * fakeRateEffective  ); 
      }

      if ( M_e1e2 > 70.0 && M_e1e2 < 110.0 ){
        FillUserTH1D("Mee_70_110_Preselection", M_e1e2, min_prescale * fakeRateEffective ) ;
        if ( sT_eejj > 600 ) 	 FillUserTH1D("Mee_70_110_ST600_Preselection", M_e1e2, min_prescale * fakeRateEffective ) ;
        if      ( isEBEB ) FillUserTH1D( "Mee_EBEB_70_110_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
        else if ( isEBEE ) FillUserTH1D( "Mee_EBEE_70_110_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
        else if ( isEEEE ) FillUserTH1D( "Mee_EEEE_70_110_PAS", M_e1e2, min_prescale * fakeRateEffective  ); 
        if      ( isEB   ) FillUserTH1D( "Mee_EB_70_110_PAS"  , M_e1e2, min_prescale * fakeRateEffective  ); 
      }

      //--------------------------------------------------------------------------
      // Region of interest plots
      //-------------------------------------------------------------------------- 

      if ( passed_region_of_interest ) { 


        FillUserTH1D("CorrIsolation_1stEle_ROI"         , Ele1_CorrIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("DeltaEtaTrkSC_1stEle_ROI"         , Ele1_DeltaEtaTrkSC                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("EcalIsolation_1stEle_ROI"         , Ele1_EcalIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HcalIsolation_1stEle_ROI"         , Ele1_HcalIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("TrkIsolation_1stEle_ROI"          , Ele1_TrkIsolation                   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HasMatchedPhot_1stEle_ROI"        , Ele1_HasMatchedPhot                 , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HoE_1stEle_ROI"                   , Ele1_HoE                            , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("LeadVtxDistXY_1stEle_ROI"         , Ele1_LeadVtxDistXY                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("LeadVtxDistZ_1stEle_ROI"          , Ele1_LeadVtxDistZ                   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("MissingHits_1stEle_ROI"           , Ele1_MissingHits                    , min_prescale * fakeRateEffective   ); 
        if ( fabs(Ele1_SCEta) < eleEta_bar ) { 
          FillUserTH1D("SigmaIEtaIEta_Barrel_1stEle_ROI", Ele1_Full5x5SigmaIEtaIEta           , min_prescale * fakeRateEffective   ); 
        }
        else if ( fabs(Ele1_SCEta) > eleEta_end1_min && fabs(Ele2_SCEta) > eleEta_end2_max ){
          FillUserTH1D("SigmaIEtaIEta_Endcap_1stEle_ROI", Ele1_Full5x5SigmaIEtaIEta           , min_prescale * fakeRateEffective   ); 
        }

        FillUserTH1D("CorrIsolation_2ndEle_ROI"         , Ele2_CorrIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("DeltaEtaTrkSC_2ndEle_ROI"         , Ele2_DeltaEtaTrkSC                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("EcalIsolation_2ndEle_ROI"         , Ele2_EcalIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HcalIsolation_2ndEle_ROI"         , Ele2_HcalIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("TrkIsolation_2ndEle_ROI"          , Ele2_TrkIsolation                   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HasMatchedPhot_2ndEle_ROI"        , Ele2_HasMatchedPhot                 , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HoE_2ndEle_ROI"                   , Ele2_HoE                            , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("LeadVtxDistXY_2ndEle_ROI"         , Ele2_LeadVtxDistXY                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("LeadVtxDistZ_2ndEle_ROI"          , Ele2_LeadVtxDistZ                   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("MissingHits_2ndEle_ROI"           , Ele2_MissingHits                    , min_prescale * fakeRateEffective   ); 
        if ( fabs(Ele2_SCEta) < eleEta_bar ) { 
          FillUserTH1D("SigmaIEtaIEta_Barrel_2ndEle_ROI", Ele2_Full5x5SigmaIEtaIEta           , min_prescale * fakeRateEffective   ); 
        }
        else if ( fabs(Ele2_SCEta) > eleEta_end1_min && fabs(Ele2_SCEta) < eleEta_end2_max ){
          FillUserTH1D("SigmaIEtaIEta_Endcap_2ndEle_ROI", Ele2_Full5x5SigmaIEtaIEta           , min_prescale * fakeRateEffective   ); 
        }

        FillUserTH1D("Me1j1_ROI"           , M_e1j1                                         , min_prescale * fakeRateEffective );
        FillUserTH1D("Ptee_ROI"            , Pt_e1e2                                        , min_prescale * fakeRateEffective );
        FillUserTH1D("Eta1stJet_ROI"       , Jet1_Eta                               , min_prescale * fakeRateEffective );
        FillUserTH1D("Eta2ndJet_ROI"       , Jet2_Eta                               , min_prescale * fakeRateEffective );
        FillUserTH1D("Eta1stEle_ROI"	    , Ele1_SCEta                                  , min_prescale * fakeRateEffective );
        FillUserTH1D("Eta2ndEle_ROI"	    , Ele2_SCEta                                  , min_prescale * fakeRateEffective );
        FillUserTH1D("Phi1stJet_ROI"       , Jet1_Phi                               , min_prescale * fakeRateEffective );
        FillUserTH1D("Phi2ndJet_ROI"       , Jet2_Phi                               , min_prescale * fakeRateEffective );
        FillUserTH1D("Phi1stEle_ROI"	    , Ele1_Phi                                  , min_prescale * fakeRateEffective );
        FillUserTH1D("Phi2ndEle_ROI"	    , Ele2_Phi                                  , min_prescale * fakeRateEffective );
        FillUserTH2D("MeeVsST_ROI"         , M_e1e2                                , sT_eejj, min_prescale * fakeRateEffective );	   
        FillUserTH1D("Mee_ROI"		    , M_e1e2                                         , min_prescale * fakeRateEffective );
        FillUserTH1D("sT_zjj_ROI"          , sT_zjj                                         , min_prescale * fakeRateEffective );
        FillUserTH1D("nVertex_ROI"         , nVertex                                        , min_prescale * fakeRateEffective );
        FillUserTH1D("EleChargeSum_ROI"    , Ele1_Charge + Ele2_Charge            , min_prescale * fakeRateEffective );
        FillUserTH1D("nJet_ROI"            , nJet_ptCut                             , min_prescale * fakeRateEffective );
        FillUserTH1D("Mej_selected_avg_ROI", M_ej_avg                                       , min_prescale * fakeRateEffective );
        FillUserTH1D("Meejj_ROI"           , M_eejj                                         , min_prescale * fakeRateEffective );
        FillUserTH1D("Meej_ROI"            , M_eej                                          , min_prescale * fakeRateEffective );
        FillUserTH1D("Mejj_ROI"            , M_ejj                                          , min_prescale * fakeRateEffective );
        FillUserTH1D("minDR_ZJet_ROI"      , min_DeltaR_Zj                                  , min_prescale * fakeRateEffective );
        FillUserTH1D("DR_ZJet1_ROI"        , DR_ZJ1                                         , min_prescale * fakeRateEffective );
        FillUserTH1D("DR_ZJet2_ROI"        , DR_ZJ2                                         , min_prescale * fakeRateEffective );
        FillUserTH1D("MET_ROI"             , PFMET_Type1_Pt                                 , min_prescale * fakeRateEffective );
        FillUserTH1D("Mjj_ROI"             , M_j1j2                                         , min_prescale * fakeRateEffective );
        FillUserTH1D("sT_ROI"              , sT_eejj                                        , min_prescale * fakeRateEffective );
        FillUserTH1D("sTlep_ROI"           , Ele1_Pt    + Ele2_Pt                 , min_prescale * fakeRateEffective );
        FillUserTH1D("sTjet_ROI"           , Jet1_Pt + Jet2_Pt              , min_prescale * fakeRateEffective );
        FillUserTH1D("Pt1stEle_ROI"        , Ele1_Pt                                   , min_prescale * fakeRateEffective );
        FillUserTH1D("Pt2ndEle_ROI"        , Ele2_Pt                                   , min_prescale * fakeRateEffective );
        FillUserTH1D("Pt1stJet_ROI"        , Jet1_Pt                                , min_prescale * fakeRateEffective );
        FillUserTH1D("Pt2ndJet_ROI"        , Jet2_Pt                                , min_prescale * fakeRateEffective );
        FillUserTH1D("sTfrac_Jet1_ROI"     , Jet1_Pt / sT_eejj                      , min_prescale * fakeRateEffective );
        FillUserTH1D("sTfrac_Jet2_ROI"     , Jet2_Pt / sT_eejj                      , min_prescale * fakeRateEffective );
        FillUserTH1D("sTfrac_Ele1_ROI"     , Ele1_Pt / sT_eejj                         , min_prescale * fakeRateEffective );
        FillUserTH1D("sTfrac_Ele2_ROI"     , Ele2_Pt / sT_eejj                         , min_prescale * fakeRateEffective );
        FillUserTH1D("sTfrac_Jet_ROI"      , ( Jet1_Pt + Jet2_Pt ) / sT_eejj, min_prescale * fakeRateEffective );
        FillUserTH1D("sTfrac_Ele_ROI"      , ( Ele1_Pt + Ele2_Pt )       / sT_eejj, min_prescale * fakeRateEffective );
        FillUserTH1D("Ptj1j2_ROI"            , Pt_j1j2                                      , min_prescale * fakeRateEffective );
        FillUserTH1D("Ptee_Minus_Ptj1j2_ROI" , Pt_e1e2 - Pt_j1j2                            , min_prescale * fakeRateEffective );

        if ( nJet_ptCut > 2 ) { 	 
          FillUserTH1D( "M_e1j3_ROI" , M_e1j3, min_prescale * fakeRateEffective ) ;
          FillUserTH1D( "M_e2j3_ROI" , M_e2j3, min_prescale * fakeRateEffective ) ;
          FillUserTH1D( "M_j1j3_ROI" , M_j1j3, min_prescale * fakeRateEffective ) ;
          FillUserTH1D( "M_j2j3_ROI" , M_j2j3, min_prescale * fakeRateEffective ) ;
          FillUserTH1D( "M_eejjj_ROI", M_eejjj,min_prescale * fakeRateEffective ) ;
          FillUserTH1D( "Ptj1j2j3_ROI"            , Pt_j1j2j3           , min_prescale * fakeRateEffective );
          FillUserTH1D( "Ptj2j3_ROI"              , Pt_j2j3             , min_prescale * fakeRateEffective );
          FillUserTH1D( "Ptj1j3_ROI"              , Pt_j1j3             , min_prescale * fakeRateEffective );
          FillUserTH1D( "Ptee_Minus_Ptj1j2j3_ROI" , Pt_e1e2 - Pt_j1j2j3 , min_prescale * fakeRateEffective );
        }
      }


    } // End preselection 
  } // End loop over events

  std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

