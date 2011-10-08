//////////////////////////////////////////////////////////////////////////
//
// 'MJJ ANALYSIS'
// 
// 09/2011 - Francesco Santanastasio 
//
// http://root.cern.ch/root/html/tutorials/roofit/index.html
//
//  using namespace RooFit
// .include /afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms8/include/
//
//  fit options at http://root.cern.ch/root/html/RooAbsPdf.html --> fitTo
//
/////////////////////////////////////////////////////////////////////////


// //#ifndef __CINT__
// #include "RooGlobalFunc.h"
// //#endif
// #include "RooRealVar.h"
// #include "RooDataSet.h"
// #include "RooDataHist.h"
// #include "RooGaussian.h"
// #include "TCanvas.h"
// #include "RooPlot.h"
// #include "TTree.h"
// #include "TH1D.h"
// #include "TRandom.h"

//using namespace RooFit ;

//void MjjAnalysis()
{

  gROOT->Reset();

  //TFile inputfile("/home/santanas/Axigluons/data/output_fromAFS/axigluons_enujj/5fb-1_Summer11MC_AxigluonW_enujj_27092011/output_cutTable_axigluons_enujj/analysisClass_axigluons_enujj_plots.root");
  TFile inputfile("/home/santanas/Axigluons/data/output_fromAFS/axigluons_enujj/5fb-1_Summer11MC_AxigluonW_enujj_noSHERPA_07102011/output_cutTable_axigluons_enujj/analysisClass_axigluons_enujj_plots.root");

  //bkg
  TH1D *h_mjj = (TH1D*)inputfile.Get( "histo1D__ALLBKG__M_j1j2" ); 
  //TH1D *h_mjj = (TH1D*)inputfile.Get( "histo1D__WJet_Madgraph__M_j1j2" ); 
  //TH1D *h_mjj = (TH1D*)inputfile.Get( "histo1D__TTbar_Madgraph__M_j1j2" ); 
  h_mjj->Rebin(1);
  //h_mjj->Draw("HISTE");
  //TODO Create histogram with different binning

  //signal
  TH1D *h_mjj_AG150 = (TH1D*)inputfile.Get( "histo1D__AxigluonW_M150__M_j1j2" ); 
  TH1D *h_mjj_AG500 = (TH1D*)inputfile.Get( "histo1D__AxigluonW_M500__M_j1j2" ); 
  TH1D *h_mjj_AG1000 = (TH1D*)inputfile.Get( "histo1D__AxigluonW_M1000__M_j1j2" ); 
  TH1D *h_mjj_AG1500 = (TH1D*)inputfile.Get( "histo1D__AxigluonW_M1500__M_j1j2" ); 

  // Declare observable mass
  RooRealVar mass("mass","dijet mass (GeV)",200,1000);

  // Make plot of binned dataset showing Poisson error bars (RooFit default)
  RooPlot* frame = mass.frame(Title("Dijet mass distribution background"));
 
  // Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
  RooDataHist RooHist_mjj("RooHist_mjj","RooHist_mjj",mass,Import(*h_mjj));
  RooDataHist RooHist_mjj_AG150("RooHist_mjj_AG150","RooHist_mjj_AG150",mass,Import(*h_mjj_AG150));
  RooDataHist RooHist_mjj_AG500("RooHist_mjj_AG500","RooHist_mjj_AG500",mass,Import(*h_mjj_AG500));
  RooDataHist RooHist_mjj_AG1000("RooHist_mjj_AG1000","RooHist_mjj_AG1000",mass,Import(*h_mjj_AG1000));
  RooDataHist RooHist_mjj_AG1500("RooHist_mjj_AG1500","RooHist_mjj_AG1500",mass,Import(*h_mjj_AG1500));

  // To construct a proper p.d.f, the formula expression is explicitly normalized internally by dividing 
  // it by a numeric integral of the expresssion over x in the range
  //RooRealVar CoM_energy("CoM_energy","CoM_energy",7000) ;

  //Default PDF (DEF)
  RooRealVar P0_DEF("P0_DEF","P0_DEF",100,0.,100000) ;
  RooRealVar P1_DEF("P1_DEF","P1_DEF",3.89,-1000,1000) ;
  RooRealVar P2_DEF("P2_DEF","P2_DEF",6.5,-100.,100) ;
  RooRealVar P3_DEF("P3_DEF","P3_DEF",0.54,-100.,100) ;
  // fix parameters
  //   P1_DEF.setConstant(kTRUE);
  //   P2_DEF.setConstant(kTRUE);
  //   P3_DEF.setConstant(kTRUE);
  //RooGenericPdf pdf_DEF("pdf_DEF","pdf_DEF","P0_DEF * pow( (1-mass/7000) , P1_DEF ) / pow( (mass/7000) , P2_DEF + P3_DEF * log(mass/7000) )",RooArgSet(mass,P0_DEF,P1_DEF,P2_DEF,P3_DEF)) ; 
  RooGenericPdf pdf_DEF("pdf_DEF","pdf_DEF","pow( (1-mass/7000) , P1_DEF ) / pow( (mass/7000) , P2_DEF + P3_DEF * log(mass/7000) )",RooArgSet(mass,P1_DEF,P2_DEF,P3_DEF)) ; 
  RooExtendPdf pdf_DEF_EXT("pdf_DEF_EXT","pdf_DEF_EXT",pdf_DEF,P0_DEF);

  //Alternative 1 (SYSA)
  RooRealVar P0_SYSA("P0_SYSA","P0_SYSA",0.1,0.,100) ;
  RooRealVar P1_SYSA("P1_SYSA","P1_SYSA",1,-100,100) ;
  RooRealVar P2_SYSA("P2_SYSA","P2_SYSA",1,-100,100) ;
  RooRealVar P3_SYSA("P3_SYSA","P3_SYSA",1,0.,100) ;
  RooGenericPdf pdf_SYSA("pdf_SYSA","pdf_SYSA","P0_SYSA * pow( 1 - mass/7000 + P3_SYSA * pow( mass/7000 , 2), P1_SYSA  ) / pow(mass,P2_SYSA)",RooArgSet(mass,P0_SYSA,P1_SYSA,P2_SYSA,P3_SYSA)) ;

  //Alternative 2 (SYSB)
  RooRealVar P0_SYSB("P0_SYSB","P0_SYSB",0.1,0.,100) ;
  RooRealVar P1_SYSB("P1_SYSB","P1_SYSB",1,-100,100) ;
  RooRealVar P2_SYSB("P2_SYSB","P2_SYSB",1,-100,100) ;
  RooGenericPdf pdf_SYSB("pdf_SYSB","pdf_SYSB","P0_SYSB * pow( 1- mass/7000 , P1_SYSB ) / pow(mass,P2_SYSB)",RooArgSet(mass,P0_SYSB,P1_SYSB,P2_SYSB)) ;

  //Signal 
  RooHistPdf pdf_AG150("pdf_AG150","pdf_AG150",RooArgSet(mass),RooHist_mjj_AG150);
  RooHistPdf pdf_AG500("pdf_AG500","pdf_AG500",RooArgSet(mass),RooHist_mjj_AG500);
  RooHistPdf pdf_AG1000("pdf_AG1000","pdf_AG1000",RooArgSet(mass),RooHist_mjj_AG1000);
  RooHistPdf pdf_AG1500("pdf_AG1500","pdf_AG1500",RooArgSet(mass),RooHist_mjj_AG1500);

  //Fits

  //test
  bool do_GenData = 0;

  //bkg only
  bool plot_data = 0;
  bool do_DEF = 0;
  bool do_SYSA = 0;
  bool do_SYSB = 0;

  //signal and background
  bool do_DEF_mod = 1;

  if(plot_data)
    {
      RooHist_mjj.plotOn(frame,Name("RooHist_mjj"),SumW2Error(kTRUE)); 
    }

  if(do_GenData)
    {
      RooDataSet* data = pdf_DEF_EXT.generate(mass,100000) ;
      data->plotOn(frame,Name("mydata")); 
      RooFitResult* fitres_DEF_data = pdf_DEF_EXT.fitTo(*data,SumW2Error(kTRUE),PrintLevel(1),Save(kTRUE)) ;
      fitres_DEF_data->Print("v");
      pdf_DEF_EXT.plotOn(frame,Name("pdf_DEF_EXT")) ;
      pdf_DEF_EXT.paramOn(frame) ;
      double chi2_DEF = frame->chiSquare("pdf_DEF_EXT","mydata",4) ;
      cout << "chi2 (DEF): " << chi2_DEF << endl;
      cout << endl;
    }

  if(do_DEF)
    {
      cout << endl;
      cout << "===================================================" << endl;
      cout << "===================== FIT DEF =====================" << endl;
      cout << "===================================================" << endl;
      cout << endl;
      RooFitResult* fitres_DEF = pdf_DEF_EXT.fitTo(RooHist_mjj,SumW2Error(kTRUE),PrintLevel(1),Save(kTRUE)) ;
      fitres_DEF->Print("v");
      pdf_DEF_EXT.plotOn(frame,Name("pdf_DEF_EXT")) ;
      pdf_DEF_EXT.paramOn(frame) ;
      double chi2_DEF = frame->chiSquare("pdf_DEF_EXT","RooHist_mjj",4) ;
      cout << "chi2 (DEF): " << chi2_DEF << endl;
      cout << endl;
    }
 
  if(do_DEF_mod)   
    {
      cout << endl;
      cout << "===================================================" << endl;
      cout << "================ FIT DEF (modified) ===============" << endl;
      cout << "===================================================" << endl;
      cout << endl;
      //my background model
      RooFitResult* fitres_DEF = pdf_DEF_EXT.fitTo(RooHist_mjj,SumW2Error(kTRUE),PrintLevel(1),Save(kTRUE)) ;
      //my signal+background model
      //double N_AG500 = RooHist_mjj_AG500.sumEntries();
      double N_AG500 = RooHist_mjj_AG500.sumEntries() * 1;    
      double N_BKG = P0_DEF.getVal();
      double frac_AG500_init = N_AG500 / ( N_AG500 + N_BKG);
      RooRealVar frac_AG500("frac_AG500","fraction of signal AG500", frac_AG500_init , 0., 1.) ;
      RooAddPdf  pdf_AG500_plus_DEF_EXT("pdf_AG500_plus_DEF_EXT", "pdf_AG500_plus_DEF_EXT", RooArgList(pdf_AG500,pdf_DEF_EXT), RooArgList(frac_AG500));

      // bkg only
      //       RooDataSet* data = pdf_DEF_EXT.generate(mass,P0_DEF.getVal()) ;
      //       data->plotOn(frame,Name("mydata")); 
      //       RooFitResult* fitres_DEF_data = pdf_DEF_EXT.fitTo(*data,SumW2Error(kFALSE),PrintLevel(1),Save(kTRUE)) ;
      //       fitres_DEF_data->Print("v");
      //       pdf_DEF_EXT.plotOn(frame,Name("pdf_DEF_EXT")) ;
      //       pdf_DEF_EXT.paramOn(frame) ;
      //       double chi2_DEF = frame->chiSquare("pdf_DEF_EXT","mydata",4) ;
      //       cout << "chi2 (DEF): " << chi2_DEF << endl;
      //       cout << endl;           

      // signal + bkg
      RooDataSet* data = pdf_AG500_plus_DEF_EXT.generate(mass, int( N_AG500 + N_BKG ) ) ;
      data->plotOn(frame,Name("mydata")); 
      RooFitResult* fitres_DEF_data = pdf_AG500_plus_DEF_EXT.fitTo(*data,SumW2Error(kFALSE),PrintLevel(1),Save(kTRUE)) ;
      fitres_DEF_data->Print("v");
      pdf_AG500_plus_DEF_EXT.plotOn(frame,Name("pdf_AG500_plus_DEF_EXT"),LineColor(kRed)) ;
      pdf_AG500_plus_DEF_EXT.plotOn(frame,Name("pdf_AG500_plus_DEF_EXT"),Components(RooArgSet(pdf_AG500)),DrawOption("F"),FillColor(kGreen)) ;
      pdf_AG500_plus_DEF_EXT.plotOn(frame,Name("pdf_AG500_plus_DEF_EXT"),Components(RooArgSet(pdf_DEF_EXT)),LineColor(kBlue),LineStyle(kDashed)) ;
      pdf_AG500_plus_DEF_EXT.paramOn(frame) ;
      double chi2_AG500_plus_DEF = frame->chiSquare("pdf_AG500_plus_DEF_EXT","mydata",4) ;
      cout << "chi2 (AG500 + DEF): " << chi2_AG500_plus_DEF << endl;
      cout << endl;           

      cout << "N_AG500 true : " << N_AG500 << endl;
      cout << "N_AG500 fit : " << data->sumEntries() * frac_AG500.getVal() << endl;
      cout << "N_BKG true : " << N_BKG << endl;
      cout << "N_BKG fit : " << P0_DEF.getVal() << endl;

    }

  if(do_SYSA)
    {
      cout << endl;
      cout << "====================================================" << endl;
      cout << "===================== FIT SYSA =====================" << endl;
      cout << "====================================================" << endl;
      cout << endl;
      RooFitResult* fitres_SYSA = pdf_SYSA.fitTo(RooHist_mjj,SumW2Error(kTRUE),PrintLevel(-1),Save(kTRUE)) ;      
      fitres_SYSA->Print("v");
      pdf_SYSA.plotOn(frame,Name("pdf_SYSA"),LineColor(kGreen)) ;
      pdf_SYSA.paramOn(frame) ;
      double chi2_SYSA = frame->chiSquare("pdf_SYSA","RooHist_mjj",4) ;
      cout << "chi2 (SYSA): " << chi2_SYSA << endl;
      cout << endl;
    }

  if(do_SYSB)
    {
      cout << endl;
      cout << "====================================================" << endl;
      cout << "===================== FIT SYSB =====================" << endl;
      cout << "====================================================" << endl;
      cout << endl;
      RooFitResult* fitres_SYSB = pdf_SYSB.fitTo(RooHist_mjj,SumW2Error(kTRUE),PrintLevel(-1),Save(kTRUE)) ;
      fitres_SYSB->Print("v");
      pdf_SYSB.plotOn(frame,Name("pdf_SYSB"),LineColor(kRed)) ;
      pdf_SYSB.paramOn(frame) ;
      double chi2_SYSB = frame->chiSquare("pdf_SYSB","RooHist_mjj",3) ;
      cout << "chi2 (SYSB): " << chi2_SYSB << endl; 
    }

  //Draw frame
  frame->GetYaxis()->SetRangeUser(0.01,30000) ;
  frame->Draw();
 
}
