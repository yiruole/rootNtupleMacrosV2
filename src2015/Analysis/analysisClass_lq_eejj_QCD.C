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

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile){}

  analysisClass::~analysisClass(){}

void analysisClass::Loop()
{
  std::cout << "analysisClass::Loop() begins" <<std::endl;   

  //--------------------------------------------------------------------------
  // Final selection mass points
  //--------------------------------------------------------------------------

  //const int n_lq_mass = 19;

  //int LQ_MASS[n_lq_mass] = { 
  //  300,  350,  400, 450, 500, 550,  600,  650,
  //  700,  750,  800, 850, 900, 950, 1000, 1050,
  //  1100, 1150, 1200
  //};

  const int n_lq_mass = 37;
  int LQ_MASS[n_lq_mass] = { 
    200,  250,
    300,  350,  400, 450, 500, 550,  600,  650,
    700,  750,  800, 850, 900, 950, 1000, 1050,
    1100, 1150, 1200, 1250, 1300, 1350, 1400, 1450,
    1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850,
    1900, 1950, 2000
  };

  //const int n_lq_mass = 1;
  //int LQ_MASS[n_lq_mass] = { 6502012 };

  std::vector<bool> passed_vector;

  char cut_name[100];
  //--------------------------------------------------------------------------
  // Decide which plots to save (default is to save everything)
  //--------------------------------------------------------------------------

  fillSkim                         ( !true  ) ;
  fillAllPreviousCuts              ( true  ) ;
  fillAllOtherCuts                 ( true  ) ;
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

  //--------------------------------------------------------------------------
  // QCD Fake Rate loading part
  //--------------------------------------------------------------------------
  std::string qcdFileName = getPreCutString1("QCDFakeRateFileName");
  //FIXME TODO eventually remove hardcoding of graph name
  QCDFakeRate qcdFakeRate(qcdFileName);

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

  //CreateUserTH1D( "M_j1j3_ROI"            ,    200 , 0       , 2000	 );    
  //CreateUserTH1D( "M_j2j3_ROI"            ,    200 , 0       , 2000	 ); 
  //CreateUserTH1D( "M_e1j3_ROI"            ,    200 , 0       , 2000	 );    
  //CreateUserTH1D( "M_e2j3_ROI"            ,    200 , 0       , 2000	 ); 
  //CreateUserTH1D( "M_eejjj_ROI"           ,    500 , 0       , 5000	 ); 

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
  //CreateUserTH1D( "sTfrac_Jet1_ROI"       ,   100  ,  0.0    , 1.0      );
  //CreateUserTH1D( "sTfrac_Jet2_ROI"       ,   100  ,  0.0    , 1.0      );
  //CreateUserTH1D( "sTfrac_Ele1_ROI"       ,   100  ,  0.0    , 1.0      );
  //CreateUserTH1D( "sTfrac_Ele2_ROI"       ,   100  ,  0.0    , 1.0      );
  //CreateUserTH1D( "sTfrac_Jet_ROI"        ,   100  ,  0.0    , 1.0      );
  //CreateUserTH1D( "sTfrac_Ele_ROI"        ,   100  ,  0.0    , 1.0      );
  CreateUserTH1D( "nElectron_PAS"         ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "nMuon_PAS"             ,    5   , -0.5    , 4.5      );
  CreateUserTH1D( "nJet_PAS"              ,    10  , -0.5    , 9.5      );
  CreateUserTH1D( "nJet_PASandMee100"        ,    10  , -0.5    , 9.5      );
  //CreateUserTH1D( "nJet_ROI"              ,    10  , -0.5    , 9.5      );
  CreateUserTH1D( "Pt1stEle_PAS"	   , 	100 , 0       , 1000     ); 
  CreateUserTH1D( "PtHeep1stEle_PAS"	   , 	100 , 0       , 1000     ); 
  CreateUserTH1D( "Pt1stEle_PASandMee100" , 	100 , 0       , 1000     ); 
  //CreateUserTH1D( "Pt1stEle_ROI"	   , 	100 , 0       , 1000     ); 
  CreateUserTH1D( "Eta1stEle_PAS"	   , 	100 , -5      , 5	 ); 
  CreateUserTH1D( "SCEta1stEle_PAS"	   , 	100 , -5      , 5	 ); 
  CreateUserTH1D( "DeltaEtaEleTrk1stEle_Presel", 400, -0.5,   0.5 );
  CreateUserTH1D( "Phi1stEle_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserTH1D( "Pt2ndEle_PAS"	   , 	300 , 0       , 3000     ); 
  CreateUserTH1D( "PtHeep2ndEle_PAS"	   , 	300 , 0       , 3000     ); 
  CreateUserTH1D( "Pt2ndEle_PASandMee100" , 	300 , 0       , 3000     ); 
  //CreateUserTH1D( "Pt2ndEle_ROI"	   , 	300 , 0       , 3000     ); 
  CreateUserTH1D( "Eta2ndEle_PAS"	   , 	100 , -5      , 5	 ); 
  CreateUserTH1D( "SCEta2ndEle_PAS"	   , 	100 , -5      , 5	 ); 
  CreateUserTH1D( "DeltaEtaEleTrk2ndEle_Presel", 400, -0.5,   0.5 );
  CreateUserTH1D( "Phi2ndEle_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserTH1D( "Charge1stEle_PAS"	   , 	2   , -1.0001 , 1.0001	 ); 
  CreateUserTH1D( "Charge2ndEle_PAS"	   , 	2   , -1.0001 , 1.0001	 ); 
  CreateUserTH1D( "EleChargeSum_PAS"         ,    3   , -2.5    , 2.5  );
  CreateUserTH1D( "EleChargeSum_PASandMee100",    3   , -2.5    , 2.5  );
  //CreateUserTH1D( "EleChargeSum_ROI"         ,    3   , -2.5    , 2.5  );
  CreateUserTH1D( "MET_PAS"               ,    200 , 0       , 1000	 ); 
  //CreateUserTH1D( "MET_ROI"               ,    200 , 0       , 1000	 ); 
  CreateUserTH1D( "METPhi_PAS"		   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserTH1D( "Pt1stJet_PAS"          ,    100 , 0       , 1000	 ); 
  CreateUserTH1D( "Pt2ndJet_PAS"          ,    100 , 0       , 1000	 ); 
  CreateUserTH1D( "Pt1stJet_PASandMee100" ,    100 , 0       , 1000	 ); 
  CreateUserTH1D( "Pt2ndJet_PASandMee100" ,    100 , 0       , 1000	 ); 
  //CreateUserTH1D( "Pt1stJet_ROI"          ,    100 , 0       , 1000	 ); 
  //CreateUserTH1D( "Pt2ndJet_ROI"          ,    100 , 0       , 1000	 ); 
  CreateUserTH1D( "Eta1stJet_PAS"         ,    100 , -5      , 5	 ); 
  CreateUserTH1D( "Eta2ndJet_PAS"         ,    100 , -5      , 5	 ); 
  CreateUserTH1D( "Phi1stJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserTH1D( "Phi2ndJet_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserTH1D( "sTlep_PAS"             ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "sTjet_PAS"             ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "sTlep_PASandMee100"    ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "sTjet_PASandMee100"    ,    200 , 0       , 2000	 ); 
  //CreateUserTH1D( "sTlep_ROI"             ,    200 , 0       , 2000	 ); 
  //CreateUserTH1D( "sTjet_ROI"             ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "sT_PAS"                ,    300   , 0       , 3000	 ); 
  CreateUserTH1D( "sT_zjj_PAS"            ,    300   , 0       , 3000	  ); 
  CreateUserTH1D( "sT_zjj_PASandMee100"   ,    300   , 0       , 3000	  ); 
  //CreateUserTH1D( "sT_zjj_ROI"                      ,    200   , 0       , 2000	  ); 
  CreateUserTH1D( "sT_PASandMee100"       ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "sT_PASandMee110"       ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "sT_PASandMee120"       ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "sT_PASandMee130"       ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "sT_PASandMee140"       ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "sT_PASandMee150"       ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "sT_PASandMee160"       ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "sT_PASandMee170"       ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "sT_PASandMee180"       ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "sT_PASandMee190"       ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "sT_PASandMee200"       ,    200 , 0       , 2000	 ); 
  //CreateUserTH1D( "sT_ROI"                ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "Mjj_PAS"		   ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "Mjj_PASandMee100"	   ,    200 , 0       , 2000	 ); 
  //CreateUserTH1D( "Mjj_ROI"		   ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "Mee_PAS"		   ,    200 , 0       , 2000	 ); 
  //CreateUserTH1D( "Mee_ROI"		   ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "Mee_PASandST445"       ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "MTenu_PAS"             ,    200 , 0       , 1000	 ); 
  CreateUserTH1D( "Me1j1_PAS"             ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "Me1j1_PASandMee100"    ,    200 , 0       , 2000	 ); 
  //CreateUserTH1D( "Me1j1_ROI"             ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "Me1j2_PAS"             ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "Me2j1_PAS"             ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "Me2j2_PAS"             ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "Me1j_selected_PAS"     ,    200 , 0       , 2000	 ); 
  CreateUserTH1D( "Me2j_selected_PAS"     ,    200 , 0       , 2000     );
  CreateUserTH1D( "Mej_selected_min_PAS"  ,    200 , 0       , 2000     ); 
  CreateUserTH1D( "Mej_selected_max_PAS"  ,    200 , 0       , 2000     ); 
  CreateUserTH1D( "Mej_minmax_PAS"        ,    200 , 0       , 2000     ); 
  CreateUserTH1D( "Mej_selected_avg_PAS"  ,    200 , 0       , 2000     );
  CreateUserTH1D( "Mej_selected_avg_PASandMee100"  ,    200 , 0       , 2000     );
  //CreateUserTH1D( "Mej_selected_avg_ROI"  ,    200 , 0       , 2000     );
  //CreateUserTH1D( "Meejj_ROI"             ,    400 , 0       , 4000     );
  CreateUserTH1D( "Mejj_PAS"              ,    400 , 0       , 4000     );
  CreateUserTH1D( "Meej_PAS"              ,    400 , 0       , 4000     );
  //CreateUserTH1D( "Mejj_ROI"              ,    400 , 0       , 4000     );
  //CreateUserTH1D( "Meej_ROI"              ,    400 , 0       , 4000     );
  CreateUserTH1D( "Meejj_PAS"             ,    400 , 0       , 4000     );
  // muon kinematics
  CreateUserTH1D( "Pt1stMuon_PAS"	   , 	100 , 0       , 1000     ); 
  CreateUserTH1D( "Eta1stMuon_PAS"	   , 	100 , -5      , 5	 ); 
  CreateUserTH1D( "Phi1stMuon_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 
  CreateUserTH1D( "Pt2ndMuon_PAS"	   , 	300 , 0       , 3000     ); 
  CreateUserTH1D( "Eta2ndMuon_PAS"	   , 	100 , -5      , 5	 ); 
  CreateUserTH1D( "Phi2ndMuon_PAS"	   , 	60  , -3.1416 , +3.1416	 ); 

  //CreateUserTH1D( "Eta1stJet_ROI"                   ,    100   , -5      , 5	  ); 
  //CreateUserTH1D( "Eta2ndJet_ROI"                   ,    100   , -5      , 5	  ); 
  //CreateUserTH1D( "Eta1stEle_ROI"	             , 	100    , -5      , 5	  ); 
  //CreateUserTH1D( "Eta2ndEle_ROI"	             , 	100    , -5      , 5	  ); 

  //CreateUserTH1D( "Phi1stJet_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
  //CreateUserTH1D( "Phi2ndJet_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
  //CreateUserTH1D( "Phi1stEle_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 
  //CreateUserTH1D( "Phi2ndEle_ROI"	             , 	 60    , -3.1416 , +3.1416  ); 

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


  //CreateUserTH1D( "Ptj1j2j3_ROI"                    ,    200 , 0       , 2000     );
  //CreateUserTH1D( "Ptj1j2_ROI"                      ,    200 , 0       , 2000     );
  //CreateUserTH1D( "Ptj2j3_ROI"                      ,    200 , 0       , 2000     );
  //CreateUserTH1D( "Ptj1j3_ROI"                      ,    200 , 0       , 2000     );

  //CreateUserTH1D( "Ptee_Minus_Ptj1j2_ROI"           ,    200 , -500    , 500      );
  //CreateUserTH1D( "Ptee_Minus_Ptj1j2j3_ROI"         ,    200 , -500    , 500      );


  CreateUserTH1D( "Ptee_PAS"              ,    200 , 0       , 2000     );
  CreateUserTH1D( "Ptee_PASandMee100"     ,    200 , 0       , 2000     );
  //CreateUserTH1D( "Ptee_ROI"              ,    200 , 0       , 2000     );

  CreateUserTH1D( "DCotTheta1stEle_PAS"   ,    100 , 0.0, 1.0);
  CreateUserTH1D( "Dist1stEle_PAS"        ,    100 , 0.0, 1.0);  
  CreateUserTH1D( "DCotTheta2ndEle_PAS"   ,    100 , 0.0, 1.0);
  CreateUserTH1D( "Dist2ndEle_PAS"        ,    100 , 0.0, 1.0);  

  CreateUserTH1D( "nVertex_PAS"                     ,    101   , -0.5   , 100.5	 ) ; 
  CreateUserTH1D( "nVertex_PASandMee100"            ,    101   , -0.5   , 100.5	 ) ; 
  //CreateUserTH1D( "nVertex_ROI"                     ,    101   , -0.5   , 100.5	 ) ; 

  CreateUserTH1D( "DR_Ele1Jet1_PAS"	   , 	getHistoNBins("DR_Ele1Jet1"), getHistoMin("DR_Ele1Jet1"), getHistoMax("DR_Ele1Jet1")     ) ; 
  CreateUserTH1D( "DR_Ele1Jet2_PAS"	   , 	getHistoNBins("DR_Ele1Jet2"), getHistoMin("DR_Ele1Jet2"), getHistoMax("DR_Ele1Jet2")     ) ; 
  CreateUserTH1D( "DR_Ele2Jet1_PAS"	   , 	getHistoNBins("DR_Ele2Jet1"), getHistoMin("DR_Ele2Jet1"), getHistoMax("DR_Ele2Jet1")     ) ; 
  CreateUserTH1D( "DR_Ele2Jet2_PAS"	   , 	getHistoNBins("DR_Ele2Jet2"), getHistoMin("DR_Ele2Jet2"), getHistoMax("DR_Ele2Jet2")     ) ; 
  CreateUserTH1D( "DR_Jet1Jet2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  CreateUserTH1D( "DR_Ele1Ele2_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  CreateUserTH1D( "minDR_EleJet_PAS"	   , 	getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  CreateUserTH1D( "minDR_ZJet_PAS"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  //CreateUserTH1D( "minDR_ZJet_ROI"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 

  CreateUserTH1D( "DR_ZJet1_PAS"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  //CreateUserTH1D( "DR_ZJet1_ROI"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  CreateUserTH1D( "DR_ZJet2_PAS"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
  //CreateUserTH1D( "DR_ZJet2_ROI"        ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 


  CreateUserTH2D( "Me1jVsMe2j_selected",     200, 0, 2000, 200, 0, 2000) ;
  CreateUserTH2D( "Me1jVsMe2j_rejected",     200, 0, 2000, 200, 0, 2000) ;


  CreateUserTH2D( "MeeVsST_PAS"                 ,     200, 0, 2000, 200, 0, 2000) ;
  CreateUserTH2D( "MeeVsST_PASandMee100"        ,     200, 0, 2000, 200, 0, 2000) ;
  //CreateUserTH2D( "MeeVsST_ROI"                 ,     200, 0, 2000, 200, 0, 2000) ;


  CreateUserTH1D( "Mee_80_100_Preselection", 200, 60, 120 );
  CreateUserTH1D( "Mee_70_110_Preselection", 200, 60, 120 );
  CreateUserTH1D( "Mee_70_110_ST600_Preselection", 200, 60, 120 );

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
  CreateUserTH1D( "GeneratorWeight", 100, -2.0 , 2.0 );

  CreateUserTH1D("BeamSpotDXY_1stEle_PAS"                   , 200,  0.0 ,   0.5  ); CreateUserTH1D("BeamSpotDXY_2ndEle_PAS"                   , 200,  0.0 ,   0.5  );
  CreateUserTH1D("Classif_1stEle_PAS"                       , 5  , -0.5 ,   4.5  ); CreateUserTH1D("Classif_2ndEle_PAS"                       , 5  , -0.5 ,   4.5  );
  CreateUserTH1D("CorrIsolation_1stEle_PAS"                 , 200,-25.0 ,  25.0  ); CreateUserTH1D("CorrIsolation_2ndEle_PAS"                 , 200,-25.0 ,  25.0  );
  CreateUserTH1D("DeltaEtaTrkSC_1stEle_PAS"                 , 200, -0.01,   0.01 ); CreateUserTH1D("DeltaEtaTrkSC_2ndEle_PAS"                 , 200, -0.01,   0.01 );
  CreateUserTH1D("DeltaPhiTrkSC_1stEle_PAS"                 , 200, -0.1 ,   0.1  ); CreateUserTH1D("DeltaPhiTrkSC_2ndEle_PAS"                 , 200, -0.1 ,   0.1  );
  CreateUserTH1D("Full5x5E1x5OverE5x5_1stEle_PAS"                  , 200,  0.0 ,   2.0  ); CreateUserTH1D("Full5x5E1x5OverE5x5_2ndEle_PAS"                  , 200,  0.0 ,   2.0  );
  CreateUserTH1D("Full5x5E2x5OverE5x5_1stEle_PAS"                  , 200,  0.0 ,   2.0  ); CreateUserTH1D("Full5x5E2x5OverE5x5_2ndEle_PAS"                  , 200,  0.0 ,   2.0  );
  CreateUserTH1D("EcalIsolation_1stEle_PAS"                 , 200,  0.0 ,  20.0  ); CreateUserTH1D("EcalIsolation_2ndEle_PAS"                 , 200,  0.0 ,  20.0  );
  CreateUserTH1D("HcalIsolation_1stEle_PAS"                 , 200,  0.0 ,  20.0  ); CreateUserTH1D("HcalIsolation_2ndEle_PAS"                 , 200,  0.0 ,  20.0  );
  CreateUserTH1D("TrkIsolation_1stEle_PAS"                  , 200,  0.0,    5.0  ); CreateUserTH1D("TrkIsolation_2ndEle_PAS"                  , 200,  0.0,    5.0  );
  CreateUserTH1D("Energy_1stEle_PAS"                        , 200,  0.0 ,3000.0  ); CreateUserTH1D("Energy_2ndEle_PAS"                        , 200,  0.0 ,3000.0  );
  CreateUserTH1D("FBrem_1stEle_PAS"                         , 200,-10.0 ,  10.0  ); CreateUserTH1D("FBrem_2ndEle_PAS"                         , 200,-10.0 ,  10.0  );
  CreateUserTH1D("GsfCtfCharge_1stEle_PAS"                  , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfCtfCharge_2ndEle_PAS"                  , 2,   -0.5 ,   1.5  );
  CreateUserTH1D("GsfCtfScPixCharge_1stEle_PAS"             , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfCtfScPixCharge_2ndEle_PAS"             , 2,   -0.5 ,   1.5  );
  CreateUserTH1D("GsfScPixCharge_1stEle_PAS"                , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfScPixCharge_2ndEle_PAS"                , 2,   -0.5 ,   1.5  );
  CreateUserTH1D("HasMatchedPhot_1stEle_PAS"                , 2,   -0.5 ,   1.5  ); CreateUserTH1D("HasMatchedPhot_2ndEle_PAS"                , 2,   -0.5 ,   1.5  );
  CreateUserTH1D("HoE_1stEle_PAS"                           , 200,  0.0 ,   0.05 ); CreateUserTH1D("HoE_2ndEle_PAS"                           , 200,  0.0 ,   0.05 );
  CreateUserTH1D("LeadVtxDistXY_1stEle_PAS"                 , 200, -0.05,   0.05 ); CreateUserTH1D("LeadVtxDistXY_2ndEle_PAS"                 , 200, -0.05,   0.05 );
  CreateUserTH1D("LeadVtxDistZ_1stEle_PAS"                  , 200, -0.2 ,   0.2  ); CreateUserTH1D("LeadVtxDistZ_2ndEle_PAS"                  , 200, -0.2 ,   0.2  );
  CreateUserTH1D("MissingHits_1stEle_PAS"                   , 2  , -0.5,    1.5  ); CreateUserTH1D("MissingHits_2ndEle_PAS"                   , 2  , -0.5,    1.5  );
  CreateUserTH1D("NBrems_1stEle_PAS"                        , 11 , -0.5,   10.5  ); CreateUserTH1D("NBrems_2ndEle_PAS"                        , 11 , -0.5,   10.5  );
  CreateUserTH1D("EnergyORawEnergy_1stEle_PAS"              , 200,  0.9,    1.4  ); CreateUserTH1D("EnergyORawEnergy_2ndEle_PAS"              , 200,  0.9,    1.4  );
  CreateUserTH1D("SigmaEtaEta_Barrel_1stEle_PAS"            , 200,  0.0,    0.02 ); CreateUserTH1D("SigmaEtaEta_Barrel_2ndEle_PAS"            , 200,  0.0,    0.02 );
  CreateUserTH1D("SigmaEtaEta_Endcap_1stEle_PAS"            , 200,  0.0,    0.1  ); CreateUserTH1D("SigmaEtaEta_Endcap_2ndEle_PAS"            , 200,  0.0,    0.1  );
  CreateUserTH1D("Full5x5SigmaIEtaIEta_Barrel_1stEle_PAS"   , 200,  0.0,    0.04 ); CreateUserTH1D("Full5x5SigmaIEtaIEta_Barrel_2ndEle_PAS"   , 200,  0.0,    0.04 );
  CreateUserTH1D("Full5x5SigmaIEtaIEta_Endcap_1stEle_PAS"   , 200,  0.0,    0.1  ); CreateUserTH1D("Full5x5SigmaIEtaIEta_Endcap_2ndEle_PAS"   , 200,  0.0,    0.1  );
  CreateUserTH1D("TrkPtOPt_1stEle_PAS"                      , 200,  0.0,  100.0  ); CreateUserTH1D("TrkPtOPt_2ndEle_PAS"                      , 200,  0.0,  100.0  );
  CreateUserTH1D("ValidFrac_1stEle_PAS"                     , 200,  0.0 ,   2.0  ); CreateUserTH1D("ValidFrac_2ndEle_PAS"                     , 200,  0.0 ,   2.0  );

  CreateUserTH1D("BeamSpotDXY_1stEle_PASandMee100"          , 200,  0.0 ,   0.5  ); CreateUserTH1D("BeamSpotDXY_2ndEle_PASandMee100"          , 200,  0.0 ,   0.5  );
  CreateUserTH1D("Classif_1stEle_PASandMee100"              , 5  , -0.5 ,   4.5  ); CreateUserTH1D("Classif_2ndEle_PASandMee100"              , 5  , -0.5 ,   4.5  );
  CreateUserTH1D("CorrIsolation_1stEle_PASandMee100"        , 200,-25.0 ,  25.0  ); CreateUserTH1D("CorrIsolation_2ndEle_PASandMee100"        , 200,-25.0 ,  25.0  );
  CreateUserTH1D("DeltaEtaTrkSC_1stEle_PASandMee100"        , 200, -0.01,   0.01 ); CreateUserTH1D("DeltaEtaTrkSC_2ndEle_PASandMee100"        , 200, -0.01,   0.01 );
  CreateUserTH1D("DeltaPhiTrkSC_1stEle_PASandMee100"        , 200, -0.1 ,   0.1  ); CreateUserTH1D("DeltaPhiTrkSC_2ndEle_PASandMee100"        , 200, -0.1 ,   0.1  );
  CreateUserTH1D("Full5x5E1x5OverE5x5_1stEle_PASandMee100"         , 200,  0.0 ,   2.0  ); CreateUserTH1D("Full5x5E1x5OverE5x5_2ndEle_PASandMee100"         , 200,  0.0 ,   2.0  );
  CreateUserTH1D("Full5x5E2x5OverE5x5_1stEle_PASandMee100"         , 200,  0.0 ,   2.0  ); CreateUserTH1D("Full5x5E2x5OverE5x5_2ndEle_PASandMee100"         , 200,  0.0 ,   2.0  );
  CreateUserTH1D("EcalIsolation_1stEle_PASandMee100"        , 200,  0.0 ,  20.0  ); CreateUserTH1D("EcalIsolation_2ndEle_PASandMee100"        , 200,  0.0 ,  20.0  );
  CreateUserTH1D("HcalIsolation_1stEle_PASandMee100"        , 200,  0.0 ,  20.0  ); CreateUserTH1D("HcalIsolation_2ndEle_PASandMee100"        , 200,  0.0 ,  20.0  );
  CreateUserTH1D("TrkIsolation_1stEle_PASandMee100"         , 200,  0.0,    5.0  ); CreateUserTH1D("TrkIsolation_2ndEle_PASandMee100"         , 200,  0.0,    5.0  );
  CreateUserTH1D("Energy_1stEle_PASandMee100"               , 200,  0.0 ,3000.0  ); CreateUserTH1D("Energy_2ndEle_PASandMee100"               , 200,  0.0 ,3000.0  );
  CreateUserTH1D("FBrem_1stEle_PASandMee100"                , 200,-10.0 ,  10.0  ); CreateUserTH1D("FBrem_2ndEle_PASandMee100"                , 200,-10.0 ,  10.0  );
  CreateUserTH1D("GsfCtfCharge_1stEle_PASandMee100"         , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfCtfCharge_2ndEle_PASandMee100"         , 2,   -0.5 ,   1.5  );
  CreateUserTH1D("GsfCtfScPixCharge_1stEle_PASandMee100"    , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfCtfScPixCharge_2ndEle_PASandMee100"    , 2,   -0.5 ,   1.5  );
  CreateUserTH1D("GsfScPixCharge_1stEle_PASandMee100"       , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfScPixCharge_2ndEle_PASandMee100"       , 2,   -0.5 ,   1.5  );
  CreateUserTH1D("HasMatchedPhot_1stEle_PASandMee100"       , 2,   -0.5 ,   1.5  ); CreateUserTH1D("HasMatchedPhot_2ndEle_PASandMee100"       , 2,   -0.5 ,   1.5  );
  CreateUserTH1D("HoE_1stEle_PASandMee100"                  , 200,  0.0 ,   0.05 ); CreateUserTH1D("HoE_2ndEle_PASandMee100"                  , 200,  0.0 ,   0.05 );
  CreateUserTH1D("LeadVtxDistXY_1stEle_PASandMee100"        , 200, -0.05,   0.05 ); CreateUserTH1D("LeadVtxDistXY_2ndEle_PASandMee100"        , 200, -0.05,   0.05 );
  CreateUserTH1D("LeadVtxDistZ_1stEle_PASandMee100"         , 200, -0.2 ,   0.2  ); CreateUserTH1D("LeadVtxDistZ_2ndEle_PASandMee100"         , 200, -0.2 ,   0.2  );
  CreateUserTH1D("MissingHits_1stEle_PASandMee100"          , 2  , -0.5,    1.5  ); CreateUserTH1D("MissingHits_2ndEle_PASandMee100"          , 2  , -0.5,    1.5  );
  CreateUserTH1D("NBrems_1stEle_PASandMee100"               , 11 , -0.5,   10.5  ); CreateUserTH1D("NBrems_2ndEle_PASandMee100"               , 11 , -0.5,   10.5  );
  CreateUserTH1D("EnergyORawEnergy_1stEle_PASandMee100"     , 200,  0.9,    1.4  ); CreateUserTH1D("EnergyORawEnergy_2ndEle_PASandMee100"     , 200,  0.9,    1.4  );
  CreateUserTH1D("SigmaEtaEta_Barrel_1stEle_PASandMee100"   , 200,  0.0,    0.02 ); CreateUserTH1D("SigmaEtaEta_Barrel_2ndEle_PASandMee100"   , 200,  0.0,    0.02 );
  CreateUserTH1D("SigmaEtaEta_Endcap_1stEle_PASandMee100"   , 200,  0.0,    0.1  ); CreateUserTH1D("SigmaEtaEta_Endcap_2ndEle_PASandMee100"   , 200,  0.0,    0.1  );
  CreateUserTH1D("Full5x5SigmaIEtaIEta_Barrel_1stEle_PASandMee100" , 200,  0.0,    0.02 ); CreateUserTH1D("Full5x5SigmaIEtaIEta_Barrel_2ndEle_PASandMee100" , 200,  0.0,    0.02 );
  CreateUserTH1D("Full5x5SigmaIEtaIEta_Endcap_1stEle_PASandMee100" , 200,  0.0,    0.1  ); CreateUserTH1D("Full5x5SigmaIEtaIEta_Endcap_2ndEle_PASandMee100" , 200,  0.0,    0.1  );
  CreateUserTH1D("TrkPtOPt_1stEle_PASandMee100"             , 200,  0.0,  100.0  ); CreateUserTH1D("TrkPtOPt_2ndEle_PASandMee100"             , 200,  0.0,  100.0  );
  CreateUserTH1D("ValidFrac_1stEle_PASandMee100"            , 200,  0.0 ,   2.0  ); CreateUserTH1D("ValidFrac_2ndEle_PASandMee100"            , 200,  0.0 ,   2.0  );


  //CreateUserTH1D("BeamSpotDXY_1stEle_ROI"                   , 200,  0.0 ,   0.5  ); CreateUserTH1D("BeamSpotDXY_2ndEle_ROI"                   , 200,  0.0 ,   0.5  );
  //CreateUserTH1D("Classif_1stEle_ROI"                       , 5  , -0.5 ,   4.5  ); CreateUserTH1D("Classif_2ndEle_ROI"                       , 5  , -0.5 ,   4.5  );
  //CreateUserTH1D("CorrIsolation_1stEle_ROI"                 , 200,-25.0 ,  25.0  ); CreateUserTH1D("CorrIsolation_2ndEle_ROI"                 , 200,-25.0 ,  25.0  );
  //CreateUserTH1D("DeltaEtaTrkSC_1stEle_ROI"                 , 200, -0.01,   0.01 ); CreateUserTH1D("DeltaEtaTrkSC_2ndEle_ROI"                 , 200, -0.01,   0.01 );
  //CreateUserTH1D("DeltaPhiTrkSC_1stEle_ROI"                 , 200, -0.1 ,   0.1  ); CreateUserTH1D("DeltaPhiTrkSC_2ndEle_ROI"                 , 200, -0.1 ,   0.1  );
  //CreateUserTH1D("E1x5OverE5x5_1stEle_ROI"                  , 200,  0.0 ,   2.0  ); CreateUserTH1D("E1x5OverE5x5_2ndEle_ROI"                  , 200,  0.0 ,   2.0  );
  //CreateUserTH1D("E2x5OverE5x5_1stEle_ROI"                  , 200,  0.0 ,   2.0  ); CreateUserTH1D("E2x5OverE5x5_2ndEle_ROI"                  , 200,  0.0 ,   2.0  );
  //CreateUserTH1D("EcalIsolation_1stEle_ROI"                 , 200,  0.0 ,  20.0  ); CreateUserTH1D("EcalIsolation_2ndEle_ROI"                 , 200,  0.0 ,  20.0  );
  //CreateUserTH1D("HcalIsolation_1stEle_ROI"                 , 200,  0.0 ,  20.0  ); CreateUserTH1D("HcalIsolation_2ndEle_ROI"                 , 200,  0.0 ,  20.0  );
  //CreateUserTH1D("TrkIsolation_1stEle_ROI"                  , 200,  0.0,    5.0  ); CreateUserTH1D("TrkIsolation_2ndEle_ROI"                  , 200,  0.0,    5.0  );
  //CreateUserTH1D("Energy_1stEle_ROI"                        , 200,  0.0 ,3000.0  ); CreateUserTH1D("Energy_2ndEle_ROI"                        , 200,  0.0 ,3000.0  );
  //CreateUserTH1D("FBrem_1stEle_ROI"                         , 200,-10.0 ,  10.0  ); CreateUserTH1D("FBrem_2ndEle_ROI"                         , 200,-10.0 ,  10.0  );
  //CreateUserTH1D("GsfCtfCharge_1stEle_ROI"                  , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfCtfCharge_2ndEle_ROI"                  , 2,   -0.5 ,   1.5  );
  //CreateUserTH1D("GsfCtfScPixCharge_1stEle_ROI"             , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfCtfScPixCharge_2ndEle_ROI"             , 2,   -0.5 ,   1.5  );
  //CreateUserTH1D("GsfScPixCharge_1stEle_ROI"                , 2,   -0.5 ,   1.5  ); CreateUserTH1D("GsfScPixCharge_2ndEle_ROI"                , 2,   -0.5 ,   1.5  );
  //CreateUserTH1D("HasMatchedPhot_1stEle_ROI"                , 2,   -0.5 ,   1.5  ); CreateUserTH1D("HasMatchedPhot_2ndEle_ROI"                , 2,   -0.5 ,   1.5  );
  //CreateUserTH1D("HoE_1stEle_ROI"                           , 200,  0.0 ,   0.05 ); CreateUserTH1D("HoE_2ndEle_ROI"                           , 200,  0.0 ,   0.05 );
  //CreateUserTH1D("LeadVtxDistXY_1stEle_ROI"                 , 200, -0.05,   0.05 ); CreateUserTH1D("LeadVtxDistXY_2ndEle_ROI"                 , 200, -0.05,   0.05 );
  //CreateUserTH1D("LeadVtxDistZ_1stEle_ROI"                  , 200, -0.2 ,   0.2  ); CreateUserTH1D("LeadVtxDistZ_2ndEle_ROI"                  , 200, -0.2 ,   0.2  );
  //CreateUserTH1D("MissingHits_1stEle_ROI"                   , 2  , -0.5,    1.5  ); CreateUserTH1D("MissingHits_2ndEle_ROI"                   , 2  , -0.5,    1.5  );
  //CreateUserTH1D("NBrems_1stEle_ROI"                        , 11 , -0.5,   10.5  ); CreateUserTH1D("NBrems_2ndEle_ROI"                        , 11 , -0.5,   20.5  );
  //CreateUserTH1D("EnergyORawEnergy_1stEle_ROI"              , 200,  0.9,    1.4  ); CreateUserTH1D("EnergyORawEnergy_2ndEle_ROI"              , 200,  0.9,    1.4  );
  //CreateUserTH1D("SigmaEtaEta_Barrel_1stEle_ROI"            , 200,  0.0,    0.02 ); CreateUserTH1D("SigmaEtaEta_Barrel_2ndEle_ROI"            , 200,  0.0,    0.02 );
  //CreateUserTH1D("SigmaEtaEta_Endcap_1stEle_ROI"            , 200,  0.0,    0.1  ); CreateUserTH1D("SigmaEtaEta_Endcap_2ndEle_ROI"            , 200,  0.0,    0.1  );
  //CreateUserTH1D("SigmaIEtaIEta_Barrel_1stEle_ROI"          , 200,  0.0,    0.02 ); CreateUserTH1D("SigmaIEtaIEta_Barrel_2ndEle_ROI"          , 200,  0.0,    0.02 );
  //CreateUserTH1D("SigmaIEtaIEta_Endcap_1stEle_ROI"          , 200,  0.0,    0.1  ); CreateUserTH1D("SigmaIEtaIEta_Endcap_2ndEle_ROI"          , 200,  0.0,    0.1  );
  //CreateUserTH1D("TrkPtOPt_1stEle_ROI"                      , 200,  0.0,  100.0  ); CreateUserTH1D("TrkPtOPt_2ndEle_ROI"                      , 200,  0.0,  100.0  );
  //CreateUserTH1D("ValidFrac_1stEle_ROI"                     , 200,  0.0 ,   2.0  ); CreateUserTH1D("ValidFrac_2ndEle_ROI"                     , 200,  0.0 ,   2.0  );
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
  // 3D opt cut space
  CreateUserTH3D( "OptimizationCutSpace", 200, 0, 2000, 200, 0, 2000, 200, 0, 2000);


  //--------------------------------------------------------------------------
  // Final selection plots
  //--------------------------------------------------------------------------
  bool doFinalSelections = false;
  // check if there is a final Mej specific in cutfile for any LQ mass
  for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
    int lq_mass = LQ_MASS[i_lq_mass];
    sprintf(cut_name, "min_M_ej_LQ%d"   , lq_mass );
    if(hasCut(cut_name)) {
      doFinalSelections = true;
      break;
    }
  }

  char plot_name[100];

  for (int i_lq_mass = 0; i_lq_mass < n_lq_mass ; ++i_lq_mass ) { 
    int lq_mass = LQ_MASS[i_lq_mass];
    sprintf(plot_name, "Mej_selected_avg_LQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 60  , 0 , 3000 );
    sprintf(plot_name, "Mej_selected_min_LQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 60  , 0 , 3000 );
    sprintf(plot_name, "Mej_selected_max_LQ%d"       , lq_mass ); CreateUserTH1D ( plot_name, 60  , 0 , 3000 );
    sprintf(plot_name, "Mej_minmax_LQ%d"             , lq_mass ); CreateUserTH1D ( plot_name, 60  , 0 , 3000 );
    sprintf(plot_name, "sT_eejj_LQ%d"                , lq_mass ); CreateUserTH1D ( plot_name, 30  , 0 , 3000 );
    sprintf(plot_name, "Mee_LQ%d"                    , lq_mass ); CreateUserTH1D ( plot_name, 40  , 0 , 2000 );
    sprintf(plot_name, "Mej_selected_min_vs_max_LQ%d", lq_mass ); CreateUserTH2D ( plot_name, 150  , 0 , 3000, 150  , 0 , 3000 );
    sprintf(plot_name, "DR_Ele1Jet1_LQ%d"            , lq_mass ); CreateUserTH1D ( plot_name, 
        getHistoNBins("DR_Ele1Jet1"), 
        getHistoMin  ("DR_Ele1Jet1"), 
        getHistoMax  ("DR_Ele1Jet1"));

    sprintf(plot_name, "BeamSpotDXY_1stEle_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   0.5  );
    sprintf(plot_name, "Classif_1stEle_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name , 5  , -0.5 ,   4.5  );
    sprintf(plot_name, "CorrIsolation_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200,-25.0 ,  25.0  );
    sprintf(plot_name, "DeltaEtaTrkSC_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.01,   0.01 );
    sprintf(plot_name, "DeltaPhiTrkSC_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.1 ,   0.1  );
    sprintf(plot_name, "Full5x5E1x5OverE5x5_1stEle_LQ%d"   , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );
    sprintf(plot_name, "Full5x5E2x5OverE5x5_1stEle_LQ%d"   , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );
    sprintf(plot_name, "EcalIsolation_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,  20.0  );
    sprintf(plot_name, "HcalIsolation_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,  20.0  );
    sprintf(plot_name, "TrkIsolation_1stEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    5.0  );
    sprintf(plot_name, "Energy_1stEle_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,3000.0  );
    sprintf(plot_name, "FBrem_1stEle_LQ%d"                , lq_mass ); CreateUserTH1D( plot_name , 200,-10.0 ,  10.0  );
    sprintf(plot_name, "GsfCtfCharge_1stEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
    sprintf(plot_name, "GsfCtfScPixCharge_1stEle_LQ%d"    , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
    sprintf(plot_name, "GsfScPixCharge_1stEle_LQ%d"       , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
    sprintf(plot_name, "HasMatchedPhot_1stEle_LQ%d"       , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
    sprintf(plot_name, "HoE_1stEle_LQ%d"                  , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   0.05 );
    sprintf(plot_name, "LeadVtxDistXY_1stEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.05,   0.05 );
    sprintf(plot_name, "LeadVtxDistZ_1stEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 200, -0.2 ,   0.2  );
    sprintf(plot_name, "MissingHits_1stEle_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name , 2  , -0.5,    1.5  );
    sprintf(plot_name, "NBrems_1stEle_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name , 11 , -0.5,   10.5  );
    sprintf(plot_name, "EnergyORawEnergy_1stEle_LQ%d"     , lq_mass ); CreateUserTH1D( plot_name , 200,  0.9,    1.4  );
    sprintf(plot_name, "SigmaEtaEta_Barrel_1stEle_LQ%d"   , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.02 );
    sprintf(plot_name, "SigmaEtaEta_Endcap_1stEle_LQ%d"   , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.1  );
    sprintf(plot_name, "Full5x5SigmaIEtaIEta_Barrel_1stEle_LQ%d" , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.02 );
    sprintf(plot_name, "Full5x5SigmaIEtaIEta_Endcap_1stEle_LQ%d" , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.1  );
    sprintf(plot_name, "TrkPtOPt_1stEle_LQ%d"             , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,  100.0  );
    sprintf(plot_name, "ValidFrac_1stEle_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );

    sprintf(plot_name, "BeamSpotDXY_2ndEle_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   0.5  );
    sprintf(plot_name, "Classif_2ndEle_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name , 5  , -0.5 ,   4.5  );
    sprintf(plot_name, "CorrIsolation_2ndEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200,-25.0 ,  25.0  );
    sprintf(plot_name, "DeltaEtaTrkSC_2ndEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.01,   0.01 );
    sprintf(plot_name, "DeltaPhiTrkSC_2ndEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.1 ,   0.1  );
    sprintf(plot_name, "Full5x5E1x5OverE5x5_2ndEle_LQ%d"  , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );
    sprintf(plot_name, "Full5x5E2x5OverE5x5_2ndEle_LQ%d"  , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );
    sprintf(plot_name, "EcalIsolation_2ndEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,  20.0  );
    sprintf(plot_name, "HcalIsolation_2ndEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,  20.0  );
    sprintf(plot_name, "TrkIsolation_2ndEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    5.0  );
    sprintf(plot_name, "Energy_2ndEle_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,3000.0  );
    sprintf(plot_name, "FBrem_2ndEle_LQ%d"                , lq_mass ); CreateUserTH1D( plot_name , 200,-10.0 ,  10.0  );
    sprintf(plot_name, "GsfCtfCharge_2ndEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
    sprintf(plot_name, "GsfCtfScPixCharge_2ndEle_LQ%d"    , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
    sprintf(plot_name, "GsfScPixCharge_2ndEle_LQ%d"       , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
    sprintf(plot_name, "HasMatchedPhot_2ndEle_LQ%d"       , lq_mass ); CreateUserTH1D( plot_name , 2,   -0.5 ,   1.5  );
    sprintf(plot_name, "HoE_2ndEle_LQ%d"                  , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   0.05 );
    sprintf(plot_name, "LeadVtxDistXY_2ndEle_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name , 200, -0.05,   0.05 );
    sprintf(plot_name, "LeadVtxDistZ_2ndEle_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name , 200, -0.2 ,   0.2  );
    sprintf(plot_name, "MissingHits_2ndEle_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name , 2  , -0.5,    1.5  );
    sprintf(plot_name, "NBrems_2ndEle_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name , 11 , -0.5,   10.5  );
    sprintf(plot_name, "EnergyORawEnergy_2ndEle_LQ%d"     , lq_mass ); CreateUserTH1D( plot_name , 200,  0.9,    1.4  );
    sprintf(plot_name, "SigmaEtaEta_Barrel_2ndEle_LQ%d"   , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.02 );
    sprintf(plot_name, "SigmaEtaEta_Endcap_2ndEle_LQ%d"   , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.1  );
    sprintf(plot_name, "Full5x5SigmaIEtaIEta_Barrel_2ndEle_LQ%d" , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.02 );
    sprintf(plot_name, "Full5x5SigmaIEtaIEta_Endcap_2ndEle_LQ%d" , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,    0.1  );
    sprintf(plot_name, "TrkPtOPt_2ndEle_LQ%d"             , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0,  100.0  );
    sprintf(plot_name, "ValidFrac_2ndEle_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name , 200,  0.0 ,   2.0  );

    sprintf(plot_name, "EleChargeSum_LQ%d"        , lq_mass ); CreateUserTH1D( plot_name ,      3 , -2.5    , 2.5      );
    sprintf(plot_name, "sTfrac_Jet1_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
    sprintf(plot_name, "sTfrac_Jet2_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
    sprintf(plot_name, "sTfrac_Ele1_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
    sprintf(plot_name, "sTfrac_Ele2_LQ%d"         , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
    sprintf(plot_name, "sTfrac_Jet_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
    sprintf(plot_name, "sTfrac_Ele_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1.0      );
    sprintf(plot_name, "nJet_LQ%d"                , lq_mass ); CreateUserTH1D( plot_name ,     10 , -0.5    , 9.5      );
    sprintf(plot_name, "Pt1stEle_LQ%d"	           , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
    sprintf(plot_name, "Eta1stEle_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
    sprintf(plot_name, "Phi1stEle_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , +3.1416  ); 
    sprintf(plot_name, "Pt2ndEle_LQ%d"	           , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
    sprintf(plot_name, "Eta2ndEle_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
    sprintf(plot_name, "Phi2ndEle_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , 3.1416   ); 
    sprintf(plot_name, "MET_LQ%d"                 , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 1000     ); 
    sprintf(plot_name, "Pt1stJet_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
    sprintf(plot_name, "Pt2ndJet_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
    sprintf(plot_name, "Eta1stJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
    sprintf(plot_name, "Eta2ndJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
    sprintf(plot_name, "Phi1stJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , 3.1416   ); 
    sprintf(plot_name, "Phi2ndJet_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , 3.1416   ); 
    sprintf(plot_name, "sTlep_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
    sprintf(plot_name, "sTjet_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
    sprintf(plot_name, "sT_zjj_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
    sprintf(plot_name, "Mjj_LQ%d"		   , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
    sprintf(plot_name, "Meejj_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    400 ,  0.0    , 4000     );
    sprintf(plot_name, "Mejj_LQ%d"                , lq_mass ); CreateUserTH1D( plot_name ,    400 ,  0.0    , 4000     );
    sprintf(plot_name, "Meej_LQ%d"                , lq_mass ); CreateUserTH1D( plot_name ,    400 ,  0.0    , 4000     );
    sprintf(plot_name, "Ptj1j2j3_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Ptj1j2_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Ptj2j3_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Ptj1j3_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Ptee_Minus_Ptj1j2_LQ%d"   , lq_mass ); CreateUserTH1D( plot_name ,    200 , -500    ,  500     );
    sprintf(plot_name, "Ptee_Minus_Ptj1j2j3_LQ%d" , lq_mass ); CreateUserTH1D( plot_name ,    200 , -500    ,  500     );
    sprintf(plot_name, "Ptee_LQ%d"                , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Me1j1_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Me1j2_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Me2j1_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "Me2j2_LQ%d"               , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );
    sprintf(plot_name, "M_j1j3_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );    
    sprintf(plot_name, "M_j2j3_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
    sprintf(plot_name, "M_e1j3_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     );    
    sprintf(plot_name, "M_e2j3_LQ%d"              , lq_mass ); CreateUserTH1D( plot_name ,    200 ,  0.0    , 2000     ); 
    sprintf(plot_name, "M_eejjj_LQ%d"             , lq_mass ); CreateUserTH1D( plot_name ,    500 ,  0.0    , 5000     ); 
    sprintf(plot_name, "nVertex_LQ%d"             , lq_mass ); CreateUserTH1D( plot_name ,    101 , -0.5    ,  100.5   ); 
    sprintf(plot_name, "MeeVsST_LQ%d"             , lq_mass ); CreateUserTH2D( plot_name ,    200 ,  0.0, 2000, 200, 0, 2000) ;
    sprintf(plot_name, "minDR_ZJet_LQ%d"          , lq_mass ); CreateUserTH1D( plot_name ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
    sprintf(plot_name, "DR_ZJet1_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
    sprintf(plot_name, "DR_ZJet2_LQ%d"            , lq_mass ); CreateUserTH1D( plot_name ,    getHistoNBins("DR_Jet1Jet2"), getHistoMin("DR_Jet1Jet2"), getHistoMax("DR_Jet1Jet2")     ) ; 
    // muon kinematics
    sprintf(plot_name, "Pt1stMuon_LQ%d"	           , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
    sprintf(plot_name, "Eta1stMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
    sprintf(plot_name, "Phi1stMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , +3.1416  ); 
    sprintf(plot_name, "Pt2ndMuon_LQ%d"	           , lq_mass ); CreateUserTH1D( plot_name ,    100 ,  0.0    , 1000     ); 
    sprintf(plot_name, "Eta2ndMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,    100 , -5.0    , 5.0      ); 
    sprintf(plot_name, "Phi2ndMuon_LQ%d"           , lq_mass ); CreateUserTH1D( plot_name ,     60 , -3.1416 , 3.1416   ); 

  }
  // for SF at final selections
  CreateUserTH1D( "Mee_70_110_LQ300", 200, 60, 120 );
  CreateUserTH1D( "Mee_70_110_LQ600", 200, 60, 120 );
  CreateUserTH1D( "Mee_70_110_LQ800", 200, 60, 120 );
  CreateUserTH1D( "Mee_70_110_LQ900", 200, 60, 120 );
  CreateUserTH1D( "Mee_70_110_LQ1000", 200, 60, 120 );

  //--------------------------------------------------------------------------
  // Loop over the chain
  //--------------------------------------------------------------------------

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntries();
  std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
    {
      std::cout << "ERROR: Could not read from TTree; exiting." << std::endl;
      exit(-1);
    }
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (nb < 0)
    {
      std::cout << "ERROR: Could not read entry from TTree: read " << nb << "bytes; exiting." << std::endl;
      exit(-2);
    }

    if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << "/" << nentries << std::endl;   

    //--------------------------------------------------------------------------
    // Reset the cuts
    //--------------------------------------------------------------------------

    resetCuts();

    //--------------------------------------------------------------------------
    // Check good run list
    //--------------------------------------------------------------------------

    int passedJSON = passJSON ( run, ls , isData ) ;

    //--------------------------------------------------------------------------
    // Find the right prescale for this event
    //--------------------------------------------------------------------------

    int min_prescale = 0;
    int passTrigger  = 0;

    if ( LooseEle1_hltPhotonPt > 0.0 ) { 
      if ( H_Photon22   > 0.1 && LooseEle1_hltPhotonPt >= 22.  && LooseEle1_hltPhotonPt < 30. ) { passTrigger = 1; min_prescale = H_Photon22  ; } 
      if ( H_Photon30   > 0.1 && LooseEle1_hltPhotonPt >= 30.  && LooseEle1_hltPhotonPt < 36. ) { passTrigger = 1; min_prescale = H_Photon30  ; } 
      if ( H_Photon36   > 0.1 && LooseEle1_hltPhotonPt >= 36.  && LooseEle1_hltPhotonPt < 50. ) { passTrigger = 1; min_prescale = H_Photon36  ; } 
      if ( H_Photon50   > 0.1 && LooseEle1_hltPhotonPt >= 50.  && LooseEle1_hltPhotonPt < 75. ) { passTrigger = 1; min_prescale = H_Photon50  ; } 
      if ( H_Photon75   > 0.1 && LooseEle1_hltPhotonPt >= 75.  && LooseEle1_hltPhotonPt < 90. ) { passTrigger = 1; min_prescale = H_Photon75  ; } 
      if ( H_Photon90   > 0.1 && LooseEle1_hltPhotonPt >= 90.  && LooseEle1_hltPhotonPt < 120.) { passTrigger = 1; min_prescale = H_Photon90  ; } 
      if ( H_Photon120  > 0.1 && LooseEle1_hltPhotonPt >= 120. && LooseEle1_hltPhotonPt < 175.) { passTrigger = 1; min_prescale = H_Photon120 ; } 
      if ( H_Photon175  > 0.1 && LooseEle1_hltPhotonPt >= 175.) { passTrigger = 1; min_prescale = H_Photon175 ; } 
    }
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
    // test output
    //std::cout << " prescale of H_Photon22 = " << H_Photon22  << std::endl;
    //std::cout << " prescale of H_Photon30 = " << H_Photon30  << std::endl;
    //std::cout << " prescale of H_Photon36 = " << H_Photon36  << std::endl;
    //std::cout << " prescale of H_Photon50 = " << H_Photon50  << std::endl;
    //std::cout << " prescale of H_Photon75 = " << H_Photon75  << std::endl;
    //std::cout << " prescale of H_Photon90 = " << H_Photon90  << std::endl;
    //std::cout << " prescale of H_Photon120= " << H_Photon120 << std::endl;
    //std::cout << " prescale of H_Photon175= " << H_Photon175 << std::endl;
    // test output

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

    if( fabs( LooseEle1_SCEta  ) < eleEta_bar )        ele1_isBarrel  = true;
    if( fabs( LooseEle1_SCEta  ) > eleEta_end1_min &&
        fabs( LooseEle1_SCEta  ) < eleEta_end1_max )   ele1_isEndcap1 = true;
    if( fabs( LooseEle1_SCEta  ) > eleEta_end2_min &&
        fabs( LooseEle1_SCEta  ) < eleEta_end2_max )   ele1_isEndcap2 = true;

    if( fabs( LooseEle2_SCEta  ) < eleEta_bar )        ele2_isBarrel  = true;
    if( fabs( LooseEle2_SCEta  ) > eleEta_end1_min &&
        fabs( LooseEle2_SCEta  ) < eleEta_end1_max )   ele2_isEndcap1 = true;
    if( fabs( LooseEle2_SCEta  ) > eleEta_end2_min &&
        fabs( LooseEle2_SCEta  ) < eleEta_end2_max )   ele2_isEndcap2 = true;

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
    //float fakeRate1 = qcdFakeRate.GetFakeRate(LooseEle1_SCEta,LooseEle1_PtHeep);
    //float fakeRate2 = qcdFakeRate.GetFakeRate(LooseEle2_SCEta,LooseEle2_PtHeep);
    float fakeRate1 = qcdFakeRate.GetFakeRate(LooseEle1_SCEta,LooseEle1_SCEnergy/cosh(LooseEle1_SCEta));
    float fakeRate2 = qcdFakeRate.GetFakeRate(LooseEle2_SCEta,LooseEle2_SCEnergy/cosh(LooseEle2_SCEta));

    //--------------------------------------------------------------------------
    // Finally have the effective fake rate
    //--------------------------------------------------------------------------

    //FIXME: add error on fake rate as well
    //double fakeRateEffective  = fakeRate1 * fakeRate2;
    double fakeRateEffective  = fakeRate1/(1-fakeRate1); // require loose electron to fail HEEP ID
    //if(1-fakeRate1 <= 0)
    //{
    //  cout << "ERROR: Found fakeRate1: " << fakeRate1 << " for SCEta=" << LooseEle1_SCEta << " SCEt="
    //    << LooseEle1_SCEnergy/cosh(LooseEle1_SCEta) << "=" << LooseEle1_SCEnergy << "/" << 
    //    cosh(LooseEle1_SCEta) << endl;
    //}
    if ( nLooseEle_store >= 2 ) { 							        
      fakeRateEffective *= fakeRate2/(1-fakeRate2);
    }
    // double eFakeRateEffective = fakeRateEffective * sqrt (  ( eFakeRate1 / fakeRate1 ) * ( eFakeRate1 / fakeRate1 ) +
    //					     ( eFakeRate2 / fakeRate2 ) * ( eFakeRate2 / fakeRate2 ) );
    double eFakeRateEffective = 0.0;
    // FIXME XXX TESTING
    //fakeRateEffective  = 1.0;
    //--------------------------------------------------------------------------
    // Fill variables
    //--------------------------------------------------------------------------
    FillUserTH1D("EventCount"           , 1                   , 1   ); 
    //if(fakeRateEffective * min_prescale != 1.0)
    //  std::cout << "!!!!!THIS EVENT HAD fakeRateEffective * min_prescale != 1.0: " << fakeRateEffective * min_prescale << std::endl;

    // reweighting
    fillVariableWithValue ( "Reweighting", 1, fakeRateEffective * min_prescale ) ; 

    // JSON variable
    fillVariableWithValue(   "PassJSON"                      , passedJSON              , fakeRateEffective * min_prescale ) ; 

    //--------------------------------------------------------------------------
    // Fill noise filters
    //--------------------------------------------------------------------------
    // see: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    // we filled these at skim time
    fillVariableWithValue( "PassGlobalTightHalo2016Filter" , PassGlobalTightHalo2016Filter  , fakeRateEffective * min_prescale );
    fillVariableWithValue( "PassGoodVertices"	             , PassGoodVertices               , fakeRateEffective * min_prescale );
    fillVariableWithValue( "PassHBHENoiseFilter"	         , PassHBHENoiseFilter            , fakeRateEffective * min_prescale );
    fillVariableWithValue( "PassHBHENoiseIsoFilter"	       , PassHBHENoiseIsoFilter         , fakeRateEffective * min_prescale );
    fillVariableWithValue( "PassBadEESupercrystalFilter"   , PassBadEESupercrystalFilter    , fakeRateEffective * min_prescale );
    fillVariableWithValue( "PassEcalDeadCellTrigPrim"      , PassEcalDeadCellTrigPrim       , fakeRateEffective * min_prescale );
    fillVariableWithValue( "PassChargedCandidateFilter"    , PassChargedCandidateFilter     , fakeRateEffective * min_prescale );
    fillVariableWithValue( "PassBadPFMuonFilter"           , PassBadPFMuonFilter            , fakeRateEffective * min_prescale );
    //

    // Electrons
    int PassNEle = 0;
    //if ( nLooseEle_ptCut == 2 ) PassNEle = 1;
    if ( nLooseEle_ptCut >= 2 ) PassNEle = 1;
    // we only look at events that have exactly two loose electrons (passing Pt>10)

    //FIXME fix the counting to avoid the below
    //// nLooseEle_store is the number of loose eles with pt>10
    //if(nLooseEle_store==3)
    //{
    //  // if 3 stored, 2 leads must pass Pt cut and third must fail
    //  if(LooseEle3_Pt<50 && LooseEle1_Pt>=50 && LooseEle2_Pt>=50)
    //    PassNEle=1;
    //}
    //else if(nLooseEle_store==2)
    //{
    //  // if 2 stored, OK
    //  if(LooseEle1_Pt>=50 && LooseEle2_Pt>=50)
    //      PassNEle=1;
    //}

    double M_ej_avg;
    double M_ej_min;
    double M_ej_max;

    // Muons
    int PassNMuon = 0;
    if ( nMuon_ptCut == 0 ) PassNMuon = 1;
    //FIXME SIC TEST
    //PassNMuon = 1;

    fillVariableWithValue ( "PassHLT"                        , passTrigger             , fakeRateEffective * min_prescale ) ;
    fillVariableWithValue("nEle_hltMatched",-1, fakeRateEffective * min_prescale ) ;
    fillVariableWithValue("nJet_hltMatched",-1, fakeRateEffective * min_prescale ) ;


    // Electrons								        
    fillVariableWithValue(   "PassNEle"                      , PassNEle                , fakeRateEffective * min_prescale ) ;
    if ( nLooseEle_store >= 1 ) { 							        
      fillVariableWithValue( "Ele1_Eta"                        , LooseEle1_Eta            , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "Ele1_SCEt"                       , LooseEle1_SCEnergy/cosh(LooseEle1_SCEta)            , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "Ele1_PassHEEPID"                 , LooseEle1_PassHEEPID                , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "Ele1_AbsDeltaEtaEleTrk"          , fabs(LooseEle1_Eta-LooseEle1_TrkEta), fakeRateEffective * min_prescale );
    }										        
    if ( nLooseEle_store >= 2 ) { 							        
      fillVariableWithValue( "Ele2_Eta"                        , LooseEle2_Eta            , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "Ele2_SCEt"                       , LooseEle2_SCEnergy/cosh(LooseEle2_SCEta)            , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "Ele2_PassHEEPID"                 , LooseEle2_PassHEEPID                , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "Ele2_AbsDeltaEtaEleTrk"          , fabs(LooseEle2_Eta-LooseEle2_TrkEta), fakeRateEffective * min_prescale );
      fillVariableWithValue( "Pt_e1e2"                       , Pt_e1e2             , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "M_e1e2"                        , M_e1e2                  , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "M_e1e2_opt"                    , M_e1e2                  , fakeRateEffective * min_prescale ) ;
    }

    // Jets
    fillVariableWithValue(   "nJet"                          , nJetLooseEle_ptCut      , fakeRateEffective * min_prescale ) ;
    if ( nJetLooseEle_store >= 1 ) { 						                
      fillVariableWithValue( "Jet1_Pt"                       , JetLooseEle1_Pt         , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "Jet1_Eta"                      , JetLooseEle1_Eta        , fakeRateEffective * min_prescale ) ;
    }

    if ( nJetLooseEle_store >= 2 ) { 
      fillVariableWithValue( "Jet2_Pt"                       , JetLooseEle2_Pt         , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "Jet2_Eta"                      , JetLooseEle2_Eta        , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "DR_Jet1Jet2"                   , DR_Jet1Jet2             , fakeRateEffective * min_prescale ) ;



      if ( nLooseEle_store >= 2 && nJetLooseEle_store >= 2) {
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

    double sT_zjj = Pt_e1e2 + JetLooseEle1_Pt + JetLooseEle2_Pt;

    // Muons
    fillVariableWithValue(   "PassNMuon"                     , PassNMuon               , fakeRateEffective * min_prescale ) ;

    // DeltaR
    if ( nLooseEle_store >= 2 && nJetLooseEle_store >= 1) {
      fillVariableWithValue( "DR_Ele1Jet1"                   , DR_Ele1Jet1             , fakeRateEffective * min_prescale ) ;
      fillVariableWithValue( "DR_Ele2Jet1"                   , DR_Ele2Jet1             , fakeRateEffective * min_prescale ) ;
      if(nJetLooseEle_store >= 2) {
        fillVariableWithValue( "DR_Ele1Jet2"                 , DR_Ele1Jet2             , fakeRateEffective * min_prescale ) ;
        fillVariableWithValue( "DR_Ele2Jet2"                 , DR_Ele2Jet2             , fakeRateEffective * min_prescale ) ;
      }
    }

    // sT
    if ( nLooseEle_store >= 2 && nJetLooseEle_store >= 2) {
      // SIC recompute sT using PtHeep. FIXME: this is now being done in skims
      //sT_eejj = LooseEle1_PtHeep+LooseEle2_PtHeep+JetLooseEle1_Pt+JetLooseEle2_Pt;
      fillVariableWithValue( "sT_eejj"                       , sT_eejj                 , fakeRateEffective  * min_prescale ) ;
      fillVariableWithValue( "sT_eejj_opt"                   , sT_eejj                 , fakeRateEffective  * min_prescale ) ;
      fillVariableWithValue( "Mej_min_opt"                   , M_ej_min                , fakeRateEffective  * min_prescale ) ;
    }      


    //--------------------------------------------------------------------------
    // Fill final selection cuts
    //--------------------------------------------------------------------------

    if(doFinalSelections)
    {
      for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
        int lq_mass = LQ_MASS[i_lq_mass];
        sprintf(cut_name, "M_e1e2_LQ%d"  , lq_mass ); fillVariableWithValue ( cut_name, M_e1e2  , fakeRateEffective  * min_prescale ) ;
        sprintf(cut_name, "sT_eejj_LQ%d" , lq_mass ); fillVariableWithValue ( cut_name, sT_eejj , fakeRateEffective  * min_prescale ) ;
        sprintf(cut_name, "min_M_ej_LQ%d", lq_mass ); fillVariableWithValue ( cut_name, M_ej_min, fakeRateEffective  * min_prescale ) ;
      }
    }

    //--------------------------------------------------------------------------
    // Evaluate the cuts
    //--------------------------------------------------------------------------

    evaluateCuts();

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
    // Did we pass any final selections?
    //--------------------------------------------------------------------------

    passed_vector.clear();
    if(doFinalSelections) 
    {
      for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
        int lq_mass = LQ_MASS[i_lq_mass];
        //sprintf(cut_name, "M_e1e2_LQ%d", lq_mass );
        sprintf(cut_name, "min_M_ej_LQ%d", lq_mass ); // this is actually the last cut in the cut file...!
        bool decision = bool ( passedAllPreviousCuts(cut_name) && passedCut (cut_name));
        passed_vector.push_back (decision);
      }
    }

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
      //std::cout << "We have electrons which look like (ptHEEP,eta,phi): (" << LooseEle1_PtHeep << "," << LooseEle1_Eta << "," << LooseEle1_Phi << "); FR = " << qcdFakeRate.GetFakeRate(LooseEle1_Eta,LooseEle2_PtHeep) << std::endl;
      //std::cout << "We have electrons which look like (pt,eta,phi): (" << LooseEle1_Pt << "," << LooseEle1_Eta << "," << LooseEle1_Phi << "); FR = " << qcdFakeRate.GetFakeRate(LooseEle1_Eta,LooseEle2_Pt) << std::endl;
      //std::cout << "Used fake rate=" << fakeRate1 << std::endl;
      ////float fakeRate2 = qcdFakeRate.GetFakeRate(LooseEle2_Eta,LooseEle2_PtHeep);
      //std::cout << "We have electrons which look like (pt,eta,phi): (" << LooseEle2_Pt << "," << LooseEle2_Eta << "," << LooseEle2_Phi << ")" << std::endl;
      //if(nLooseEle_store > 2)
      //  std::cout << "We have electrons which look like (pt,eta,phi): (" << LooseEle3_Pt << "," << LooseEle3_Eta << "," << LooseEle3_Phi << ")" << std::endl;

      //--------------------------------------------------------------------------
      // Recalculate some variables
      //--------------------------------------------------------------------------

      TLorentzVector e1, j1, e2, j2,j3, mu, met;
      TLorentzVector eejj, e1e2mu;
      TLorentzVector eej, ejj, ee;
      TLorentzVector e1j3, e2j3, j1j3, j2j3, j1j2, j1j2j3, eejjj;

      e1.SetPtEtaPhiM ( LooseEle1_Pt, LooseEle1_SCEta, LooseEle1_Phi, 0.0 );
      e2.SetPtEtaPhiM ( LooseEle2_Pt, LooseEle2_SCEta, LooseEle2_Phi, 0.0 );
      j1.SetPtEtaPhiM ( JetLooseEle1_Pt, JetLooseEle1_Eta, JetLooseEle1_Phi, 0.0 );
      j2.SetPtEtaPhiM ( JetLooseEle2_Pt, JetLooseEle2_Eta, JetLooseEle2_Phi, 0.0 );
      mu.SetPtEtaPhiM ( Muon1_Pt, Muon1_Eta, Muon1_Phi, 0.0 );
      //met.SetPtEtaPhiM ( PFMET_Type01XY_Pt, 0.0, PFMET_Type01XY_Phi, 0.0 );

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

      if ( nJetLooseEle_ptCut > 2 ) { 
        j3.SetPtEtaPhiM ( JetLooseEle3_Pt, JetLooseEle3_Eta, JetLooseEle3_Phi, 0.0 );

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
      if ( nJetLooseEle_ptCut > 2 ) {
        if ( DR_Ele1Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele1Jet3;
        if ( DR_Ele2Jet3 < min_DR_EleJet ) min_DR_EleJet = DR_Ele2Jet3;
      }

      //--------------------------------------------------------------------------
      // Electron quality histograms (preselection)
      //--------------------------------------------------------------------------

      FillUserTH1D("BeamSpotDXY_1stEle_PAS"           , LooseEle1_BeamSpotDXY                    , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("Classif_1stEle_PAS"               , LooseEle1_Classif                        , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("CorrIsolation_1stEle_PAS"         , LooseEle1_CorrIsolation                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("DeltaEtaTrkSC_1stEle_PAS"         , LooseEle1_DeltaEtaTrkSC                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("DeltaPhiTrkSC_1stEle_PAS"         , LooseEle1_DeltaPhiTrkSC                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("Full5x5E1x5OverE5x5_1stEle_PAS"   , LooseEle1_Full5x5E1x5OverE5x5            , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("Full5x5E2x5OverE5x5_1stEle_PAS"   , LooseEle1_Full5x5E2x5OverE5x5            , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("EcalIsolation_1stEle_PAS"         , LooseEle1_EcalIsolation                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("HcalIsolation_1stEle_PAS"         , LooseEle1_HcalIsolation                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("TrkIsolation_1stEle_PAS"          , LooseEle1_TrkIsolation                   , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("Energy_1stEle_PAS"                , LooseEle1_Energy                         , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("FBrem_1stEle_PAS"                 , LooseEle1_FBrem                          , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("GsfCtfCharge_1stEle_PAS"          , LooseEle1_GsfCtfCharge                   , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("GsfCtfScPixCharge_1stEle_PAS"     , LooseEle1_GsfCtfScPixCharge              , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("GsfScPixCharge_1stEle_PAS"        , LooseEle1_GsfScPixCharge                 , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("HasMatchedPhot_1stEle_PAS"        , LooseEle1_HasMatchedPhot                 , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("HoE_1stEle_PAS"                   , LooseEle1_HoE                            , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("LeadVtxDistXY_1stEle_PAS"         , LooseEle1_LeadVtxDistXY                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("LeadVtxDistZ_1stEle_PAS"          , LooseEle1_LeadVtxDistZ                   , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("MissingHits_1stEle_PAS"           , LooseEle1_MissingHits                    , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("NBrems_1stEle_PAS"                , LooseEle1_NBrems                         , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("ValidFrac_1stEle_PAS"             , LooseEle1_ValidFrac                      , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("EnergyORawEnergy_1stEle_PAS"      , LooseEle1_Energy / LooseEle1_RawEnergy   , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("TrkPtOPt_1stEle_PAS"              , LooseEle1_TrkPt  / LooseEle1_Pt          , min_prescale * fakeRateEffective   ); 
      if ( fabs(LooseEle1_SCEta) < eleEta_bar ) { 
        FillUserTH1D("SigmaEtaEta_Barrel_1stEle_PAS"  , LooseEle1_SigmaEtaEta                    , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("Full5x5SigmaIEtaIEta_Barrel_1stEle_PAS", LooseEle1_Full5x5SigmaIEtaIEta                  , min_prescale * fakeRateEffective   ); 
      }
      else if ( fabs(LooseEle1_SCEta) > eleEta_end1_min && fabs(LooseEle2_SCEta) < eleEta_end2_max ){
        FillUserTH1D("SigmaEtaEta_Endcap_1stEle_PAS"  , LooseEle1_SigmaEtaEta                    , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("Full5x5SigmaIEtaIEta_Endcap_1stEle_PAS", LooseEle1_Full5x5SigmaIEtaIEta    , min_prescale * fakeRateEffective   ); 
      }

      FillUserTH1D("BeamSpotDXY_2ndEle_PAS"           , LooseEle2_BeamSpotDXY                    , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("Classif_2ndEle_PAS"               , LooseEle2_Classif                        , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("CorrIsolation_2ndEle_PAS"         , LooseEle2_CorrIsolation                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("DeltaEtaTrkSC_2ndEle_PAS"         , LooseEle2_DeltaEtaTrkSC                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("DeltaPhiTrkSC_2ndEle_PAS"         , LooseEle2_DeltaPhiTrkSC                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("Full5x5E1x5OverE5x5_2ndEle_PAS"   , LooseEle2_Full5x5E1x5OverE5x5            , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("Full5x5E2x5OverE5x5_2ndEle_PAS"   , LooseEle2_Full5x5E2x5OverE5x5            , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("EcalIsolation_2ndEle_PAS"         , LooseEle2_EcalIsolation                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("HcalIsolation_2ndEle_PAS"         , LooseEle2_HcalIsolation                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("TrkIsolation_2ndEle_PAS"          , LooseEle2_TrkIsolation                   , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("Energy_2ndEle_PAS"                , LooseEle2_Energy                         , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("FBrem_2ndEle_PAS"                 , LooseEle2_FBrem                          , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("GsfCtfCharge_2ndEle_PAS"          , LooseEle2_GsfCtfCharge                   , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("GsfCtfScPixCharge_2ndEle_PAS"     , LooseEle2_GsfCtfScPixCharge              , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("GsfScPixCharge_2ndEle_PAS"        , LooseEle2_GsfScPixCharge                 , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("HasMatchedPhot_2ndEle_PAS"        , LooseEle2_HasMatchedPhot                 , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("HoE_2ndEle_PAS"                   , LooseEle2_HoE                            , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("LeadVtxDistXY_2ndEle_PAS"         , LooseEle2_LeadVtxDistXY                  , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("LeadVtxDistZ_2ndEle_PAS"          , LooseEle2_LeadVtxDistZ                   , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("MissingHits_2ndEle_PAS"           , LooseEle2_MissingHits                    , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("NBrems_2ndEle_PAS"                , LooseEle2_NBrems                         , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("ValidFrac_2ndEle_PAS"             , LooseEle2_ValidFrac                      , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("EnergyORawEnergy_2ndEle_PAS"      , LooseEle2_Energy / LooseEle2_RawEnergy   , min_prescale * fakeRateEffective   ); 
      FillUserTH1D("TrkPtOPt_2ndEle_PAS"              , LooseEle2_TrkPt  / LooseEle2_Pt          , min_prescale * fakeRateEffective   ); 
      if ( fabs(LooseEle2_SCEta) < eleEta_bar ) { 
        FillUserTH1D("SigmaEtaEta_Barrel_2ndEle_PAS"  , LooseEle2_SigmaEtaEta                    , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("Full5x5SigmaIEtaIEta_Barrel_2ndEle_PAS", LooseEle2_Full5x5SigmaIEtaIEta    , min_prescale * fakeRateEffective   ); 
      }
      else if ( fabs(LooseEle2_SCEta) > eleEta_end1_min && fabs(LooseEle2_SCEta) < eleEta_end2_max ){
        FillUserTH1D("SigmaEtaEta_Endcap_2ndEle_PAS"  , LooseEle2_SigmaEtaEta                    , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("Full5x5SigmaIEtaIEta_Endcap_2ndEle_PAS", LooseEle2_Full5x5SigmaIEtaIEta    , min_prescale * fakeRateEffective   ); 
      }

      //--------------------------------------------------------------------------
      // Preselection histograms
      //--------------------------------------------------------------------------

      FillUserTH1D( "Ptj1j2_PAS"           , Pt_j1j2                        , min_prescale * fakeRateEffective ) ;
      FillUserTH1D( "Ptee_Minus_Ptj1j2_PAS", Pt_e1e2 - Pt_j1j2              , min_prescale * fakeRateEffective ) ;

      FillUserTH1D("minDR_EleJet_PAS"     , min_DR_EleJet                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("DR_Ele1Ele2_PAS"	   , DR_Ele1Ele2                        , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("EleChargeSum_PAS"     , LooseEle1_Charge + LooseEle2_Charge, min_prescale * fakeRateEffective ) ;

      FillUserTH1D("nElectron_PAS"        , nLooseEle_ptCut                    , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("nMuon_PAS"            , nMuon_ptCut                        , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("nJet_PAS"             , nJetLooseEle_ptCut                 , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Pt1stEle_PAS"	   , LooseEle1_Pt                       , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("PtHeep1stEle_PAS"	   , LooseEle1_PtHeep                       , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Eta1stEle_PAS"	   , LooseEle1_Eta                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("SCEta1stEle_PAS"	   , LooseEle1_SCEta                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("DeltaEtaEleTrk1stEle_Presel"       , fabs(LooseEle1_Eta-LooseEle1_TrkEta)                   , min_prescale * fakeRateEffective);
      FillUserTH1D("Phi1stEle_PAS"	   , LooseEle1_Phi                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Pt2ndEle_PAS"	   , LooseEle2_Pt                       , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("PtHeep2ndEle_PAS"	   , LooseEle2_PtHeep                       , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Eta2ndEle_PAS"	   , LooseEle2_Eta                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("SCEta2ndEle_PAS"	   , LooseEle2_SCEta                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("DeltaEtaEleTrk2ndEle_Presel"       , fabs(LooseEle2_Eta-LooseEle2_TrkEta)                   , min_prescale * fakeRateEffective);
      FillUserTH1D("Phi2ndEle_PAS"	   , LooseEle2_Phi                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Charge1stEle_PAS"	   , LooseEle1_Charge                   , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Charge2ndEle_PAS"	   , LooseEle2_Charge                   , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("MET_PAS"              , PFMET_Type1XY_Pt                  , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("METPhi_PAS"	   , PFMET_Type1XY_Phi                 , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Pt1stJet_PAS"         , JetLooseEle1_Pt                    , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Pt2ndJet_PAS"         , JetLooseEle2_Pt                    , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Eta1stJet_PAS"        , JetLooseEle1_Eta                   , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Eta2ndJet_PAS"        , JetLooseEle2_Eta                   , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Phi1stJet_PAS"	   , JetLooseEle1_Phi                   , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Phi2ndJet_PAS"	   , JetLooseEle2_Phi                   , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("sTlep_PAS"            , LooseEle1_Pt + LooseEle2_Pt        , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("sTjet_PAS"            , JetLooseEle1_Pt + JetLooseEle2_Pt  , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("sT_PAS"               , sT_eejj                            , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("sT_zjj_PAS"           , sT_zjj                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Mjj_PAS"		   , M_j1j2                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Mee_PAS"		   , M_e1e2                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("MTenu_PAS"            , MT_Ele1MET                         , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Me1j1_PAS"            , M_e1j1                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Me1j2_PAS"            , M_e1j2                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Me2j1_PAS"            , M_e2j1                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Me2j2_PAS"            , M_e2j2                             , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Ptee_PAS"             , Pt_e1e2                            , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("DCotTheta1stEle_PAS"  , LooseEle1_DCotTheta                , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Dist1stEle_PAS"       , LooseEle1_Dist                     , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("DCotTheta2ndEle_PAS"  , LooseEle2_DCotTheta                , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Dist2ndEle_PAS"       , LooseEle2_Dist                     , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("nVertex_PAS"          , nVertex                            , min_prescale * fakeRateEffective ) ;
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
      // muon kinematics
      FillUserTH1D("Pt1stMuon_PAS"	   , Muon1_Pt                       , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Eta1stMuon_PAS"	   , Muon1_Eta                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Phi1stMuon_PAS"	   , Muon1_Phi                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Pt2ndMuon_PAS"	   , Muon2_Pt                       , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Eta2ndMuon_PAS"	   , Muon2_Eta                      , min_prescale * fakeRateEffective ) ;
      FillUserTH1D("Phi2ndMuon_PAS"	   , Muon2_Phi                      , min_prescale * fakeRateEffective ) ;

      FillUserTH2D("MeeVsST_PAS" , M_e1e2, sT_eejj, min_prescale * fakeRateEffective ) ;	   
      // scale factor dependence histos
      if ( nJetLooseEle_ptCut == 2 )
        FillUserTH1D("Mee_NJetEq2_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if( nJetLooseEle_ptCut == 3 )
        FillUserTH1D("Mee_NJetEq3_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if( nJetLooseEle_ptCut == 4 )
        FillUserTH1D("Mee_NJetEq4_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if( nJetLooseEle_ptCut == 5 )
        FillUserTH1D("Mee_NJetEq5_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if( nJetLooseEle_ptCut == 6 )
        FillUserTH1D("Mee_NJetEq6_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      else if( nJetLooseEle_ptCut == 7 )
        FillUserTH1D("Mee_NJetEq7_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      //
      if ( nJetLooseEle_ptCut >= 3 )
        FillUserTH1D("Mee_NJetGeq3_PAS", M_e1e2                         , min_prescale * fakeRateEffective );
      if ( nJetLooseEle_ptCut >= 4 )
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
      //
      FillUserTH3D("OptimizationCutSpace",sT_eejj,M_ej_min,M_e1e2, min_prescale * fakeRateEffective );

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

      if ( nJetLooseEle_ptCut > 2 ){ 
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

        FillUserTH1D("BeamSpotDXY_1stEle_PASandMee100"           , LooseEle1_BeamSpotDXY                    , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("Classif_1stEle_PASandMee100"               , LooseEle1_Classif                        , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("CorrIsolation_1stEle_PASandMee100"         , LooseEle1_CorrIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("DeltaEtaTrkSC_1stEle_PASandMee100"         , LooseEle1_DeltaEtaTrkSC                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("DeltaPhiTrkSC_1stEle_PASandMee100"         , LooseEle1_DeltaPhiTrkSC                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("Full5x5E1x5OverE5x5_1stEle_PASandMee100"   , LooseEle1_Full5x5E1x5OverE5x5            , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("Full5x5E2x5OverE5x5_1stEle_PASandMee100"   , LooseEle1_Full5x5E2x5OverE5x5            , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("EcalIsolation_1stEle_PASandMee100"         , LooseEle1_EcalIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HcalIsolation_1stEle_PASandMee100"         , LooseEle1_HcalIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("TrkIsolation_1stEle_PASandMee100"          , LooseEle1_TrkIsolation                   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("Energy_1stEle_PASandMee100"                , LooseEle1_Energy                         , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("FBrem_1stEle_PASandMee100"                 , LooseEle1_FBrem                          , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("GsfCtfCharge_1stEle_PASandMee100"          , LooseEle1_GsfCtfCharge                   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("GsfCtfScPixCharge_1stEle_PASandMee100"     , LooseEle1_GsfCtfScPixCharge              , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("GsfScPixCharge_1stEle_PASandMee100"        , LooseEle1_GsfScPixCharge                 , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HasMatchedPhot_1stEle_PASandMee100"        , LooseEle1_HasMatchedPhot                 , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HoE_1stEle_PASandMee100"                   , LooseEle1_HoE                            , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("LeadVtxDistXY_1stEle_PASandMee100"         , LooseEle1_LeadVtxDistXY                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("LeadVtxDistZ_1stEle_PASandMee100"          , LooseEle1_LeadVtxDistZ                   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("MissingHits_1stEle_PASandMee100"           , LooseEle1_MissingHits                    , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("NBrems_1stEle_PASandMee100"                , LooseEle1_NBrems                         , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("ValidFrac_1stEle_PASandMee100"             , LooseEle1_ValidFrac                      , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("EnergyORawEnergy_1stEle_PASandMee100"      , LooseEle1_Energy / LooseEle1_RawEnergy   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("TrkPtOPt_1stEle_PASandMee100"              , LooseEle1_TrkPt  / LooseEle1_Pt          , min_prescale * fakeRateEffective   ); 
        if ( fabs(LooseEle1_SCEta) < eleEta_bar ) { 
          FillUserTH1D("SigmaEtaEta_Barrel_1stEle_PASandMee100"  , LooseEle1_SigmaEtaEta                    , min_prescale * fakeRateEffective   ); 
          FillUserTH1D("Full5x5SigmaIEtaIEta_Barrel_1stEle_PASandMee100", LooseEle1_Full5x5SigmaIEtaIEta    , min_prescale * fakeRateEffective   ); 
        }
        else if ( fabs(LooseEle1_SCEta) > eleEta_end1_min && fabs(LooseEle2_SCEta) < eleEta_end2_max ){
          FillUserTH1D("SigmaEtaEta_Endcap_1stEle_PASandMee100"  , LooseEle1_SigmaEtaEta                    , min_prescale * fakeRateEffective   ); 
          FillUserTH1D("Full5x5SigmaIEtaIEta_Endcap_1stEle_PASandMee100", LooseEle1_Full5x5SigmaIEtaIEta    , min_prescale * fakeRateEffective   ); 
        }

        FillUserTH1D("BeamSpotDXY_2ndEle_PASandMee100"           , LooseEle2_BeamSpotDXY                    , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("Classif_2ndEle_PASandMee100"               , LooseEle2_Classif                        , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("CorrIsolation_2ndEle_PASandMee100"         , LooseEle2_CorrIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("DeltaEtaTrkSC_2ndEle_PASandMee100"         , LooseEle2_DeltaEtaTrkSC                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("DeltaPhiTrkSC_2ndEle_PASandMee100"         , LooseEle2_DeltaPhiTrkSC                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("Full5x5E1x5OverE5x5_2ndEle_PASandMee100"   , LooseEle2_Full5x5E1x5OverE5x5            , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("Full5x5E2x5OverE5x5_2ndEle_PASandMee100"   , LooseEle2_Full5x5E2x5OverE5x5            , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("EcalIsolation_2ndEle_PASandMee100"         , LooseEle2_EcalIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HcalIsolation_2ndEle_PASandMee100"         , LooseEle2_HcalIsolation                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("TrkIsolation_2ndEle_PASandMee100"          , LooseEle2_TrkIsolation                   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("Energy_2ndEle_PASandMee100"                , LooseEle2_Energy                         , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("FBrem_2ndEle_PASandMee100"                 , LooseEle2_FBrem                          , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("GsfCtfCharge_2ndEle_PASandMee100"          , LooseEle2_GsfCtfCharge                   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("GsfCtfScPixCharge_2ndEle_PASandMee100"     , LooseEle2_GsfCtfScPixCharge              , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("GsfScPixCharge_2ndEle_PASandMee100"        , LooseEle2_GsfScPixCharge                 , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HasMatchedPhot_2ndEle_PASandMee100"        , LooseEle2_HasMatchedPhot                 , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("HoE_2ndEle_PASandMee100"                   , LooseEle2_HoE                            , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("LeadVtxDistXY_2ndEle_PASandMee100"         , LooseEle2_LeadVtxDistXY                  , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("LeadVtxDistZ_2ndEle_PASandMee100"          , LooseEle2_LeadVtxDistZ                   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("MissingHits_2ndEle_PASandMee100"           , LooseEle2_MissingHits                    , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("NBrems_2ndEle_PASandMee100"                , LooseEle2_NBrems                         , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("ValidFrac_2ndEle_PASandMee100"             , LooseEle2_ValidFrac                      , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("EnergyORawEnergy_2ndEle_PASandMee100"      , LooseEle2_Energy / LooseEle2_RawEnergy   , min_prescale * fakeRateEffective   ); 
        FillUserTH1D("TrkPtOPt_2ndEle_PASandMee100"              , LooseEle2_TrkPt  / LooseEle2_Pt          , min_prescale * fakeRateEffective   ); 
        if ( fabs(LooseEle2_SCEta) < eleEta_bar ) { 
          FillUserTH1D("SigmaEtaEta_Barrel_2ndEle_PASandMee100"  , LooseEle2_SigmaEtaEta                    , min_prescale * fakeRateEffective   ); 
          FillUserTH1D("Full5x5SigmaIEtaIEta_Barrel_2ndEle_PASandMee100", LooseEle2_Full5x5SigmaIEtaIEta    , min_prescale * fakeRateEffective   ); 
        }
        else if ( fabs(LooseEle2_SCEta) > eleEta_end1_min && fabs(LooseEle2_SCEta) < eleEta_end2_max ){
          FillUserTH1D("SigmaEtaEta_Endcap_2ndEle_PASandMee100"  , LooseEle2_SigmaEtaEta                    , min_prescale * fakeRateEffective   ); 
          FillUserTH1D("Full5x5SigmaIEtaIEta_Endcap_2ndEle_PASandMee100", LooseEle2_Full5x5SigmaIEtaIEta    , min_prescale * fakeRateEffective   ); 
        }

        FillUserTH1D("Me1j1_PASandMee100"           , M_e1j1                              , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Ptee_PASandMee100"            , Pt_e1e2                             , min_prescale * fakeRateEffective ) ;
        FillUserTH2D("MeeVsST_PASandMee100" , M_e1e2, sT_eejj, min_prescale * fakeRateEffective ) ;	   
        FillUserTH1D("sT_zjj_PASandMee100"          , sT_zjj                              , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("nVertex_PASandMee100"         , nVertex                             , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sT_PASandMee100"              , sT_eejj                             , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("EleChargeSum_PASandMee100"    , LooseEle1_Charge + LooseEle2_Charge , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("nJet_PASandMee100"            , nJetLooseEle_ptCut                  , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sTlep_PASandMee100"           , LooseEle1_Pt    + LooseEle2_Pt      , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sTjet_PASandMee100"           , JetLooseEle1_Pt + JetLooseEle2_Pt   , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Mjj_PASandMee100"             , M_j1j2                              , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Pt1stEle_PASandMee100"        , LooseEle1_Pt                        , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Pt2ndEle_PASandMee100"        , LooseEle2_Pt                        , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Pt1stJet_PASandMee100"        , JetLooseEle1_Pt                     , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Pt2ndJet_PASandMee100"        , JetLooseEle2_Pt                     , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Mej_selected_avg_PASandMee100", M_ej_avg                            , min_prescale * fakeRateEffective ) ;

        FillUserTH1D("sTfrac_Jet1_PASandMee100"     , JetLooseEle1_Pt / sT_eejj                       , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sTfrac_Jet2_PASandMee100"     , JetLooseEle2_Pt / sT_eejj                       , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sTfrac_Ele1_PASandMee100"     , LooseEle1_Pt / sT_eejj                          , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sTfrac_Ele2_PASandMee100"     , LooseEle2_Pt / sT_eejj                          , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sTfrac_Jet_PASandMee100"      , ( JetLooseEle1_Pt + JetLooseEle2_Pt ) / sT_eejj , min_prescale * fakeRateEffective ) ;
        FillUserTH1D("sTfrac_Ele_PASandMee100"      , ( LooseEle1_Pt + LooseEle2_Pt ) / sT_eejj       , min_prescale * fakeRateEffective ) ;

        FillUserTH1D("Ptj1j2_PASandMee100"            , Pt_j1j2                        ,  min_prescale * fakeRateEffective ) ;
        FillUserTH1D("Ptee_Minus_Ptj1j2_PASandMee100" , Pt_e1e2 - Pt_j1j2              ,  min_prescale * fakeRateEffective ) ;

        if ( nJetLooseEle_ptCut > 2 ) { 	 
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

      //if ( passed_region_of_interest ) { 

      //  FillUserTH1D("BeamSpotDXY_1stEle_ROI"           , LooseEle1_BeamSpotDXY                    , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("Classif_1stEle_ROI"               , LooseEle1_Classif                        , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("CorrIsolation_1stEle_ROI"         , LooseEle1_CorrIsolation                  , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("DeltaEtaTrkSC_1stEle_ROI"         , LooseEle1_DeltaEtaTrkSC                  , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("DeltaPhiTrkSC_1stEle_ROI"         , LooseEle1_DeltaPhiTrkSC                  , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("E1x5OverE5x5_1stEle_ROI"          , LooseEle1_E1x5OverE5x5                   , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("E2x5OverE5x5_1stEle_ROI"          , LooseEle1_E2x5OverE5x5                   , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("EcalIsolation_1stEle_ROI"         , LooseEle1_EcalIsolation                  , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("HcalIsolation_1stEle_ROI"         , LooseEle1_HcalIsolation                  , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("TrkIsolation_1stEle_ROI"          , LooseEle1_TrkIsolation                   , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("Energy_1stEle_ROI"                , LooseEle1_Energy                         , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("FBrem_1stEle_ROI"                 , LooseEle1_FBrem                          , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("GsfCtfCharge_1stEle_ROI"          , LooseEle1_GsfCtfCharge                   , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("GsfCtfScPixCharge_1stEle_ROI"     , LooseEle1_GsfCtfScPixCharge              , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("GsfScPixCharge_1stEle_ROI"        , LooseEle1_GsfScPixCharge                 , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("HasMatchedPhot_1stEle_ROI"        , LooseEle1_HasMatchedPhot                 , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("HoE_1stEle_ROI"                   , LooseEle1_HoE                            , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("LeadVtxDistXY_1stEle_ROI"         , LooseEle1_LeadVtxDistXY                  , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("LeadVtxDistZ_1stEle_ROI"          , LooseEle1_LeadVtxDistZ                   , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("MissingHits_1stEle_ROI"           , LooseEle1_MissingHits                    , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("NBrems_1stEle_ROI"                , LooseEle1_NBrems                         , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("ValidFrac_1stEle_ROI"             , LooseEle1_ValidFrac                      , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("EnergyORawEnergy_1stEle_ROI"      , LooseEle1_Energy / LooseEle1_RawEnergy   , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("TrkPtOPt_1stEle_ROI"              , LooseEle1_TrkPt  / LooseEle1_Pt          , min_prescale * fakeRateEffective   ); 
      //  if ( fabs(LooseEle1_SCEta) < eleEta_bar ) { 
      //    FillUserTH1D("SigmaEtaEta_Barrel_1stEle_ROI"  , LooseEle1_SigmaEtaEta                    , min_prescale * fakeRateEffective   ); 
      //    FillUserTH1D("SigmaIEtaIEta_Barrel_1stEle_ROI", LooseEle1_SigmaIEtaIEta                  , min_prescale * fakeRateEffective   ); 
      //  }
      //  else if ( fabs(LooseEle1_SCEta) > eleEta_end1_min && fabs(LooseEle2_SCEta) > eleEta_end2_max ){
      //    FillUserTH1D("SigmaEtaEta_Endcap_1stEle_ROI"  , LooseEle1_SigmaEtaEta                    , min_prescale * fakeRateEffective   ); 
      //    FillUserTH1D("SigmaIEtaIEta_Endcap_1stEle_ROI", LooseEle1_SigmaIEtaIEta                  , min_prescale * fakeRateEffective   ); 
      //  }

      //  FillUserTH1D("BeamSpotDXY_2ndEle_ROI"           , LooseEle2_BeamSpotDXY                    , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("Classif_2ndEle_ROI"               , LooseEle2_Classif                        , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("CorrIsolation_2ndEle_ROI"         , LooseEle2_CorrIsolation                  , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("DeltaEtaTrkSC_2ndEle_ROI"         , LooseEle2_DeltaEtaTrkSC                  , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("DeltaPhiTrkSC_2ndEle_ROI"         , LooseEle2_DeltaPhiTrkSC                  , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("E1x5OverE5x5_2ndEle_ROI"          , LooseEle2_E1x5OverE5x5                   , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("E2x5OverE5x5_2ndEle_ROI"          , LooseEle2_E2x5OverE5x5                   , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("EcalIsolation_2ndEle_ROI"         , LooseEle2_EcalIsolation                  , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("HcalIsolation_2ndEle_ROI"         , LooseEle2_HcalIsolation                  , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("TrkIsolation_2ndEle_ROI"          , LooseEle2_TrkIsolation                   , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("Energy_2ndEle_ROI"                , LooseEle2_Energy                         , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("FBrem_2ndEle_ROI"                 , LooseEle2_FBrem                          , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("GsfCtfCharge_2ndEle_ROI"          , LooseEle2_GsfCtfCharge                   , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("GsfCtfScPixCharge_2ndEle_ROI"     , LooseEle2_GsfCtfScPixCharge              , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("GsfScPixCharge_2ndEle_ROI"        , LooseEle2_GsfScPixCharge                 , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("HasMatchedPhot_2ndEle_ROI"        , LooseEle2_HasMatchedPhot                 , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("HoE_2ndEle_ROI"                   , LooseEle2_HoE                            , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("LeadVtxDistXY_2ndEle_ROI"         , LooseEle2_LeadVtxDistXY                  , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("LeadVtxDistZ_2ndEle_ROI"          , LooseEle2_LeadVtxDistZ                   , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("MissingHits_2ndEle_ROI"           , LooseEle2_MissingHits                    , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("NBrems_2ndEle_ROI"                , LooseEle2_NBrems                         , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("ValidFrac_2ndEle_ROI"             , LooseEle2_ValidFrac                      , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("EnergyORawEnergy_2ndEle_ROI"      , LooseEle2_Energy / LooseEle2_RawEnergy   , min_prescale * fakeRateEffective   ); 
      //  FillUserTH1D("TrkPtOPt_2ndEle_ROI"              , LooseEle2_TrkPt  / LooseEle2_Pt          , min_prescale * fakeRateEffective   ); 
      //  if ( fabs(LooseEle2_Eta) < eleEta_bar ) { 
      //    FillUserTH1D("SigmaEtaEta_Barrel_2ndEle_ROI"  , LooseEle2_SigmaEtaEta                    , min_prescale * fakeRateEffective   ); 
      //    FillUserTH1D("SigmaIEtaIEta_Barrel_2ndEle_ROI", LooseEle2_SigmaIEtaIEta                  , min_prescale * fakeRateEffective   ); 
      //  }
      //  else if ( fabs(LooseEle2_Eta) > eleEta_end1_min && fabs(LooseEle2_Eta) < eleEta_end2_max ){
      //    FillUserTH1D("SigmaEtaEta_Endcap_2ndEle_ROI"  , LooseEle2_SigmaEtaEta                    , min_prescale * fakeRateEffective   ); 
      //    FillUserTH1D("SigmaIEtaIEta_Endcap_2ndEle_ROI", LooseEle2_SigmaIEtaIEta                  , min_prescale * fakeRateEffective   ); 
      //  }

      //  FillUserTH1D("Me1j1_ROI"           , M_e1j1                                         , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Ptee_ROI"            , Pt_e1e2                                        , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Eta1stJet_ROI"       , JetLooseEle1_Eta                               , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Eta2ndJet_ROI"       , JetLooseEle2_Eta                               , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Eta1stEle_ROI"	    , LooseEle1_SCEta                                  , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Eta2ndEle_ROI"	    , LooseEle2_SCEta                                  , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Phi1stJet_ROI"       , JetLooseEle1_Phi                               , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Phi2ndJet_ROI"       , JetLooseEle2_Phi                               , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Phi1stEle_ROI"	    , LooseEle1_Phi                                  , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Phi2ndEle_ROI"	    , LooseEle2_Phi                                  , min_prescale * fakeRateEffective );
      //  FillUserTH2D("MeeVsST_ROI"         , M_e1e2                                , sT_eejj, min_prescale * fakeRateEffective );	   
      //  FillUserTH1D("Mee_ROI"		    , M_e1e2                                         , min_prescale * fakeRateEffective );
      //  FillUserTH1D("sT_zjj_ROI"          , sT_zjj                                         , min_prescale * fakeRateEffective );
      //  FillUserTH1D("nVertex_ROI"         , nVertex                                        , min_prescale * fakeRateEffective );
      //  FillUserTH1D("EleChargeSum_ROI"    , LooseEle1_Charge + LooseEle2_Charge            , min_prescale * fakeRateEffective );
      //  FillUserTH1D("nJet_ROI"            , nJetLooseEle_ptCut                             , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Mej_selected_avg_ROI", M_ej_avg                                       , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Meejj_ROI"           , M_eejj                                         , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Meej_ROI"            , M_eej                                          , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Mejj_ROI"            , M_ejj                                          , min_prescale * fakeRateEffective );
      //  FillUserTH1D("minDR_ZJet_ROI"      , min_DeltaR_Zj                                  , min_prescale * fakeRateEffective );
      //  FillUserTH1D("DR_ZJet1_ROI"        , DR_ZJ1                                         , min_prescale * fakeRateEffective );
      //  FillUserTH1D("DR_ZJet2_ROI"        , DR_ZJ2                                         , min_prescale * fakeRateEffective );
      //  FillUserTH1D("MET_ROI"             , PFMET_Type01XY_Pt                              , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Mjj_ROI"             , M_j1j2                                         , min_prescale * fakeRateEffective );
      //  FillUserTH1D("sT_ROI"              , sT_eejj                                        , min_prescale * fakeRateEffective );
      //  FillUserTH1D("sTlep_ROI"           , LooseEle1_Pt    + LooseEle2_Pt                 , min_prescale * fakeRateEffective );
      //  FillUserTH1D("sTjet_ROI"           , JetLooseEle1_Pt + JetLooseEle2_Pt              , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Pt1stEle_ROI"        , LooseEle1_Pt                                   , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Pt2ndEle_ROI"        , LooseEle2_Pt                                   , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Pt1stJet_ROI"        , JetLooseEle1_Pt                                , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Pt2ndJet_ROI"        , JetLooseEle2_Pt                                , min_prescale * fakeRateEffective );
      //  FillUserTH1D("sTfrac_Jet1_ROI"     , JetLooseEle1_Pt / sT_eejj                      , min_prescale * fakeRateEffective );
      //  FillUserTH1D("sTfrac_Jet2_ROI"     , JetLooseEle2_Pt / sT_eejj                      , min_prescale * fakeRateEffective );
      //  FillUserTH1D("sTfrac_Ele1_ROI"     , LooseEle1_Pt / sT_eejj                         , min_prescale * fakeRateEffective );
      //  FillUserTH1D("sTfrac_Ele2_ROI"     , LooseEle2_Pt / sT_eejj                         , min_prescale * fakeRateEffective );
      //  FillUserTH1D("sTfrac_Jet_ROI"      , ( JetLooseEle1_Pt + JetLooseEle2_Pt ) / sT_eejj, min_prescale * fakeRateEffective );
      //  FillUserTH1D("sTfrac_Ele_ROI"      , ( LooseEle1_Pt + LooseEle2_Pt )       / sT_eejj, min_prescale * fakeRateEffective );
      //  FillUserTH1D("Ptj1j2_ROI"            , Pt_j1j2                                      , min_prescale * fakeRateEffective );
      //  FillUserTH1D("Ptee_Minus_Ptj1j2_ROI" , Pt_e1e2 - Pt_j1j2                            , min_prescale * fakeRateEffective );

      //  if ( nJetLooseEle_ptCut > 2 ) { 	 
      //    FillUserTH1D( "M_e1j3_ROI" , M_e1j3, min_prescale * fakeRateEffective ) ;
      //    FillUserTH1D( "M_e2j3_ROI" , M_e2j3, min_prescale * fakeRateEffective ) ;
      //    FillUserTH1D( "M_j1j3_ROI" , M_j1j3, min_prescale * fakeRateEffective ) ;
      //    FillUserTH1D( "M_j2j3_ROI" , M_j2j3, min_prescale * fakeRateEffective ) ;
      //    FillUserTH1D( "M_eejjj_ROI", M_eejjj,min_prescale * fakeRateEffective ) ;
      //    FillUserTH1D( "Ptj1j2j3_ROI"            , Pt_j1j2j3           , min_prescale * fakeRateEffective );
      //    FillUserTH1D( "Ptj2j3_ROI"              , Pt_j2j3             , min_prescale * fakeRateEffective );
      //    FillUserTH1D( "Ptj1j3_ROI"              , Pt_j1j3             , min_prescale * fakeRateEffective );
      //    FillUserTH1D( "Ptee_Minus_Ptj1j2j3_ROI" , Pt_e1e2 - Pt_j1j2j3 , min_prescale * fakeRateEffective );
      //  }
      //}

      //-------------------------------------------------------------------------- 
      // Final selection plots
      //-------------------------------------------------------------------------- 
      // now, we must have an Mej cut and optimization must be off to have final selections enabled
      doFinalSelections = doFinalSelections && !isOptimizationEnabled();

      if(doFinalSelections)
      {
        for (int i_lq_mass = 0; i_lq_mass < n_lq_mass; ++i_lq_mass ){ 
          int lq_mass = LQ_MASS[i_lq_mass];
          bool pass = passed_vector[i_lq_mass];
          if ( !pass ) continue;

          sprintf(plot_name, "Mej_selected_avg_LQ%d"       , lq_mass ); FillUserTH1D ( plot_name, M_ej_avg          , min_prescale * fakeRateEffective);
          sprintf(plot_name, "Mej_selected_min_LQ%d"       , lq_mass ); FillUserTH1D ( plot_name, M_ej_min          , min_prescale * fakeRateEffective);
          sprintf(plot_name, "Mej_selected_max_LQ%d"       , lq_mass ); FillUserTH1D ( plot_name, M_ej_max          , min_prescale * fakeRateEffective);
          sprintf(plot_name, "Mej_minmax_LQ%d"             , lq_mass ); FillUserTH1D ( plot_name, M_ej_min          , min_prescale * fakeRateEffective);
          sprintf(plot_name, "Mej_minmax_LQ%d"             , lq_mass ); FillUserTH1D ( plot_name, M_ej_max          , min_prescale * fakeRateEffective);
          sprintf(plot_name, "sT_eejj_LQ%d"                , lq_mass ); FillUserTH1D ( plot_name, sT_eejj           , min_prescale * fakeRateEffective);
          sprintf(plot_name, "Mee_LQ%d"                    , lq_mass ); FillUserTH1D ( plot_name, M_e1e2            , min_prescale * fakeRateEffective);
          sprintf(plot_name, "DR_Ele1Jet1_LQ%d"            , lq_mass ); FillUserTH1D ( plot_name, DR_Ele1Jet1       , min_prescale * fakeRateEffective);
          sprintf(plot_name, "Mej_selected_min_vs_max_LQ%d", lq_mass ); FillUserTH2D ( plot_name, M_ej_min, M_ej_max, min_prescale * fakeRateEffective);

          sprintf(plot_name, "BeamSpotDXY_1stEle_LQ%d"        , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_BeamSpotDXY               , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "Classif_1stEle_LQ%d"            , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_Classif                   , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "CorrIsolation_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_CorrIsolation             , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "DeltaEtaTrkSC_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_DeltaEtaTrkSC             , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "DeltaPhiTrkSC_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_DeltaPhiTrkSC             , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "Full5x5E1x5OverE5x5_1stEle_LQ%d", lq_mass );   FillUserTH1D(plot_name,  LooseEle1_Full5x5E1x5OverE5x5       , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "Full5x5E2x5OverE5x5_1stEle_LQ%d", lq_mass );   FillUserTH1D(plot_name,  LooseEle1_Full5x5E2x5OverE5x5       , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "EcalIsolation_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_EcalIsolation             , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "HcalIsolation_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_HcalIsolation             , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "TrkIsolation_1stEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_TrkIsolation              , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "Energy_1stEle_LQ%d"             , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_Energy                    , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "FBrem_1stEle_LQ%d"              , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_FBrem                     , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "GsfCtfCharge_1stEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_GsfCtfCharge              , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "GsfCtfScPixCharge_1stEle_LQ%d"  , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_GsfCtfScPixCharge         , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "GsfScPixCharge_1stEle_LQ%d"     , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_GsfScPixCharge            , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "HasMatchedPhot_1stEle_LQ%d"     , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_HasMatchedPhot            , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "HoE_1stEle_LQ%d"                , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_HoE                       , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "LeadVtxDistXY_1stEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_LeadVtxDistXY             , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "LeadVtxDistZ_1stEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_LeadVtxDistZ              , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "MissingHits_1stEle_LQ%d"        , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_MissingHits               , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "NBrems_1stEle_LQ%d"             , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_NBrems                    , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "ValidFrac_1stEle_LQ%d"          , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_ValidFrac                 , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "EnergyORawEnergy_1stEle_LQ%d"   , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_Energy / LooseEle1_RawEnergy   , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "TrkPtOPt_1stEle_LQ%d"           , lq_mass );   FillUserTH1D(plot_name,  LooseEle1_TrkPt  / LooseEle1_Pt          , min_prescale * fakeRateEffective ); 

          if ( fabs(LooseEle1_Eta) < eleEta_bar ) { 
            sprintf(plot_name, "SigmaEtaEta_Barrel_1stEle_LQ%d"  , lq_mass ); FillUserTH1D( plot_name , LooseEle1_SigmaEtaEta   , min_prescale * fakeRateEffective    ); 
            sprintf(plot_name, "Full5x5SigmaIEtaIEta_Barrel_1stEle_LQ%d", lq_mass ); FillUserTH1D( plot_name , LooseEle1_Full5x5SigmaIEtaIEta , min_prescale * fakeRateEffective    ); 
          }
          else if ( fabs(LooseEle1_Eta) > eleEta_end1_min && fabs(LooseEle2_Eta) < eleEta_end2_max ){
            sprintf(plot_name, "SigmaEtaEta_Endcap_1stEle_LQ%d"  , lq_mass ); FillUserTH1D( plot_name , LooseEle1_SigmaEtaEta   , min_prescale * fakeRateEffective    ); 
            sprintf(plot_name, "Full5x5SigmaIEtaIEta_Endcap_1stEle_LQ%d", lq_mass ); FillUserTH1D( plot_name , LooseEle1_Full5x5SigmaIEtaIEta , min_prescale * fakeRateEffective    ); 
          }

          sprintf(plot_name, "BeamSpotDXY_2ndEle_LQ%d"        , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_BeamSpotDXY               , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "Classif_2ndEle_LQ%d"            , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_Classif                   , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "CorrIsolation_2ndEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_CorrIsolation             , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "DeltaEtaTrkSC_2ndEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_DeltaEtaTrkSC             , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "DeltaPhiTrkSC_2ndEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_DeltaPhiTrkSC             , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "Full5x5E1x5OverE5x5_2ndEle_LQ%d", lq_mass );   FillUserTH1D(plot_name,  LooseEle2_Full5x5E1x5OverE5x5       , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "Full5x5E2x5OverE5x5_2ndEle_LQ%d", lq_mass );   FillUserTH1D(plot_name,  LooseEle2_Full5x5E2x5OverE5x5       , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "EcalIsolation_2ndEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_EcalIsolation             , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "HcalIsolation_2ndEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_HcalIsolation             , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "TrkIsolation_2ndEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_TrkIsolation              , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "Energy_2ndEle_LQ%d"             , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_Energy                    , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "FBrem_2ndEle_LQ%d"              , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_FBrem                     , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "GsfCtfCharge_2ndEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_GsfCtfCharge              , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "GsfCtfScPixCharge_2ndEle_LQ%d"  , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_GsfCtfScPixCharge         , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "GsfScPixCharge_2ndEle_LQ%d"     , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_GsfScPixCharge            , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "HasMatchedPhot_2ndEle_LQ%d"     , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_HasMatchedPhot            , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "HoE_2ndEle_LQ%d"                , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_HoE                       , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "LeadVtxDistXY_2ndEle_LQ%d"      , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_LeadVtxDistXY             , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "LeadVtxDistZ_2ndEle_LQ%d"       , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_LeadVtxDistZ              , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "MissingHits_2ndEle_LQ%d"        , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_MissingHits               , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "NBrems_2ndEle_LQ%d"             , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_NBrems                    , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "ValidFrac_2ndEle_LQ%d"          , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_ValidFrac                 , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "EnergyORawEnergy_2ndEle_LQ%d"   , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_Energy / LooseEle2_RawEnergy   , min_prescale * fakeRateEffective ); 
          sprintf(plot_name, "TrkPtOPt_2ndEle_LQ%d"           , lq_mass );   FillUserTH1D(plot_name,  LooseEle2_TrkPt  / LooseEle2_Pt          , min_prescale * fakeRateEffective ); 

          if ( fabs(LooseEle2_Eta) < eleEta_bar ) { 
            sprintf(plot_name, "SigmaEtaEta_Barrel_2ndEle_LQ%d"  , lq_mass ); FillUserTH1D( plot_name , LooseEle2_SigmaEtaEta   , min_prescale * fakeRateEffective    ); 
            sprintf(plot_name, "Full5x5SigmaIEtaIEta_Barrel_2ndEle_LQ%d", lq_mass ); FillUserTH1D( plot_name , LooseEle2_Full5x5SigmaIEtaIEta , min_prescale * fakeRateEffective    ); 
          }
          else if ( fabs(LooseEle2_Eta) > eleEta_end2_min && fabs(LooseEle2_Eta) < eleEta_end2_max ){
            sprintf(plot_name, "SigmaEtaEta_Endcap_2ndEle_LQ%d"  , lq_mass ); FillUserTH1D( plot_name , LooseEle2_SigmaEtaEta   , min_prescale * fakeRateEffective    ); 
            sprintf(plot_name, "Full5x5SigmaIEtaIEta_Endcap_2ndEle_LQ%d", lq_mass ); FillUserTH1D( plot_name , LooseEle2_Full5x5SigmaIEtaIEta , min_prescale * fakeRateEffective    ); 
          }

          sprintf(plot_name, "Me1j1_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , M_e1j1                                         , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Me1j2_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , M_e1j2                                         , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Me2j1_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , M_e2j1                                         , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Me2j2_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , M_e2j2                                         , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Ptee_LQ%d"              , lq_mass ); FillUserTH1D( plot_name , Pt_e1e2                                        , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Eta1stJet_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , JetLooseEle1_Eta                               , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Eta2ndJet_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , JetLooseEle2_Eta                               , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Eta1stEle_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , LooseEle1_Eta                                  , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Eta2ndEle_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , LooseEle2_Eta                                  , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Phi1stJet_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , JetLooseEle1_Phi                               , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Phi2ndJet_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , JetLooseEle2_Phi                               , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Phi1stEle_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , LooseEle1_Phi                                  , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Phi2ndEle_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , LooseEle2_Phi                                  , min_prescale * fakeRateEffective );
          sprintf(plot_name, "MeeVsST_LQ%d"           , lq_mass ); FillUserTH2D( plot_name , M_e1e2, sT_eejj                                , min_prescale * fakeRateEffective );	   
          sprintf(plot_name, "sT_zjj_LQ%d"            , lq_mass ); FillUserTH1D( plot_name , sT_zjj                                         , min_prescale * fakeRateEffective );
          sprintf(plot_name, "nVertex_LQ%d"           , lq_mass ); FillUserTH1D( plot_name , nVertex                                        , min_prescale * fakeRateEffective );
          sprintf(plot_name, "nJet_LQ%d"              , lq_mass ); FillUserTH1D( plot_name , nJetLooseEle_ptCut                             , min_prescale * fakeRateEffective );
          sprintf(plot_name, "EleChargeSum_LQ%d"      , lq_mass ); FillUserTH1D( plot_name , LooseEle1_Charge + LooseEle2_Charge            , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Meejj_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , M_eejj                                         , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Meej_LQ%d"              , lq_mass ); FillUserTH1D( plot_name , M_eej                                          , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Mejj_LQ%d"              , lq_mass ); FillUserTH1D( plot_name , M_ejj                                          , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Mjj_LQ%d"               , lq_mass ); FillUserTH1D( plot_name , M_j1j2                                         , min_prescale * fakeRateEffective );
          sprintf(plot_name, "minDR_ZJet_LQ%d"        , lq_mass ); FillUserTH1D( plot_name , min_DeltaR_Zj                                  , min_prescale * fakeRateEffective );
          sprintf(plot_name, "DR_ZJet1_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , DR_ZJ1                                         , min_prescale * fakeRateEffective );
          sprintf(plot_name, "DR_ZJet2_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , DR_ZJ2                                         , min_prescale * fakeRateEffective );
          sprintf(plot_name, "MET_LQ%d"               , lq_mass ); FillUserTH1D( plot_name , PFMET_Type1XY_Pt                              , min_prescale * fakeRateEffective );
          sprintf(plot_name, "sTlep_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , LooseEle1_Pt + LooseEle2_Pt                    , min_prescale * fakeRateEffective );
          sprintf(plot_name, "sTjet_LQ%d"             , lq_mass ); FillUserTH1D( plot_name , JetLooseEle1_Pt + JetLooseEle2_Pt              , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Pt1stEle_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , LooseEle1_Pt                                   , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Pt2ndEle_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , LooseEle2_Pt                                   , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Pt1stJet_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , JetLooseEle1_Pt                                , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Pt2ndJet_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , JetLooseEle2_Pt                                , min_prescale * fakeRateEffective );
          sprintf(plot_name, "sTfrac_Jet1_LQ%d"       , lq_mass ); FillUserTH1D( plot_name , JetLooseEle1_Pt / sT_eejj                      , min_prescale * fakeRateEffective );
          sprintf(plot_name, "sTfrac_Jet2_LQ%d"       , lq_mass ); FillUserTH1D( plot_name , JetLooseEle2_Pt / sT_eejj                      , min_prescale * fakeRateEffective );
          sprintf(plot_name, "sTfrac_Ele1_LQ%d"       , lq_mass ); FillUserTH1D( plot_name , LooseEle1_Pt / sT_eejj                         , min_prescale * fakeRateEffective );
          sprintf(plot_name, "sTfrac_Ele2_LQ%d"       , lq_mass ); FillUserTH1D( plot_name , LooseEle2_Pt / sT_eejj                         , min_prescale * fakeRateEffective );
          sprintf(plot_name, "sTfrac_Jet_LQ%d"        , lq_mass ); FillUserTH1D( plot_name , ( JetLooseEle1_Pt + JetLooseEle2_Pt ) / sT_eejj, min_prescale * fakeRateEffective );
          sprintf(plot_name, "sTfrac_Ele_LQ%d"        , lq_mass ); FillUserTH1D( plot_name , ( LooseEle1_Pt + LooseEle2_Pt ) / sT_eejj      , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Ptj1j2_LQ%d"            , lq_mass ); FillUserTH1D( plot_name , Pt_j1j2                                        , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Ptee_Minus_Ptj1j2_LQ%d" , lq_mass ); FillUserTH1D( plot_name , Pt_e1e2 - Pt_j1j2                              , min_prescale * fakeRateEffective );
          // muon kinematics
          sprintf(plot_name, "Pt1stMuon_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Muon1_Pt                                   , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Pt2ndMuon_LQ%d"          , lq_mass ); FillUserTH1D( plot_name , Muon2_Pt                                   , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Eta1stMuon_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Muon1_Eta                                  , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Eta2ndMuon_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Muon2_Eta                                  , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Phi1stMuon_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Muon1_Phi                                  , min_prescale * fakeRateEffective );
          sprintf(plot_name, "Phi2ndMuon_LQ%d"         , lq_mass ); FillUserTH1D( plot_name , Muon2_Phi                                  , min_prescale * fakeRateEffective );

        } // End final selection

        if(passedCut("sT_eejj_LQ300") && passedCut("min_M_ej_LQ300"))
          FillUserTH1D("Mee_70_110_LQ300", M_e1e2 , min_prescale * fakeRateEffective );
        if(passedCut("sT_eejj_LQ600") && passedCut("min_M_ej_LQ600"))
          FillUserTH1D("Mee_70_110_LQ600", M_e1e2 , min_prescale * fakeRateEffective );
        if(passedCut("sT_eejj_LQ800") && passedCut("min_M_ej_LQ800"))
          FillUserTH1D("Mee_70_110_LQ800", M_e1e2 , min_prescale * fakeRateEffective );
        if(passedCut("sT_eejj_LQ900") && passedCut("min_M_ej_LQ900"))
          FillUserTH1D("Mee_70_110_LQ900", M_e1e2 , min_prescale * fakeRateEffective );
        if(passedCut("sT_eejj_LQ1000") && passedCut("min_M_ej_LQ1000"))
          FillUserTH1D("Mee_70_110_LQ1000", M_e1e2 , min_prescale * fakeRateEffective );

      } // End do final selections

    } // End preselection 
  } // End loop over events

  std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

