{
  gStyle->Reset("Default");
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(10);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetStatColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  gStyle->SetPadTickY(1);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetStatFont(132);
  gStyle->SetTitleFont(132, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetTitleXOffset(0.9);
//   gStyle->SetTitleYOffset(1.05);
  gStyle->SetLabelFont(132, "XYZ");
  gStyle->SetLabelSize(0.05, "XYZ");
  
  TCanvas c1;
  c1.SetGridx();
  c1.SetGridy();

  //-----------------------------------------------
  //USER PART (modify here)  -  BEGIN

  //available optimized sT cuts vs MLQ
  const int Nmass = 10;
  float MLQ[Nmass]=  {200, 250, 280, 300, 320, 340, 370, 400, 450, 500};
  float STcut[Nmass]={350, 410, 460, 490, 520, 540, 570, 600, 640, 670};

  //MLQ values for which you need an optimal sT cut
  const int Nmass_extrapolate = 10;
  float MLQ_extrapolate[Nmass_extrapolate]={200, 250, 280, 300, 320, 340, 370, 400, 450, 500};

  //range for fit
  float minMLQfit = MLQ[0];
  float maxMLQfit = MLQ[Nmass-1];

  //USER PART                -  END
  //-----------------------------------------------

  TGraph g_STcut_vs_MLQ(Nmass, MLQ, STcut);

  g_STcut_vs_MLQ->SetTitle("");
  g_STcut_vs_MLQ->SetMarkerStyle(20);
  g_STcut_vs_MLQ->SetMarkerColor(kBlue);
  g_STcut_vs_MLQ->GetXaxis()->SetTitle("M_{LQ} [GeV]");
  g_STcut_vs_MLQ->GetYaxis()->SetTitle("Optimized S_{T} cut [GeV]");
  g_STcut_vs_MLQ->Draw("ap");

  TF1 *fitST = new TF1("fitST","pol2",minMLQfit,maxMLQfit);
  fitST->SetLineStyle(2);
  g_STcut_vs_MLQ->Fit(fitST);
//   fitST->Draw("same");


//   TLatex l;
//   l.SetTextAlign(12);
//   l.SetTextSize(0.04);
//   l.SetTextFont(132);
//   l.SetNDC();
//   l.DrawLatex(0.15,0.75,"Optimized S_{T} cut vs M_{LQ}");

  cout << endl;
  cout << "#########################" << endl;
  for(int mass=0; mass < Nmass_extrapolate; mass++)
    {


      cout << "The optimized sT cut for MLQ of "
	   << MLQ_extrapolate[mass] << " GeV is : "
	   << fitST.Eval(MLQ_extrapolate[mass]) << " GeV" << endl;

    }
  cout << "#########################" << endl;

  c1.SaveAs("STcut_vs_MLQ.eps");

}//end macro
