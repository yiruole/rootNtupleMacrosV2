{

  gROOT->Reset();
  gStyle->SetTitle("");

  TCanvas c1;
  c1.SetGridx();
  c1.SetGridy();

  //-----------------------------------------------
  //USER PART (modify here)  -  BEGIN

  //available optimized sT cuts vs MLQ
  const int Nmass = 5;
  float MLQ[Nmass]={200,250,280,300, 320}; 
  float STcut[Nmass]={350,410,460,490,520}; 

  //MLQ values for which you need an optimal sT cut
  const int Nmass_extrapolate = 1;
  float MLQ_extrapolate[Nmass_extrapolate]={340}; 

  //range for fit
  float minMLQfit = MLQ[0];
  float maxMLQfit = MLQ[Nmass-1];

  //USER PART                -  END  
  //-----------------------------------------------

  TGraph g_STcut_vs_MLQ(Nmass, MLQ, STcut);

  g_STcut_vs_MLQ->SetTitle("");
  g_STcut_vs_MLQ->SetMarkerStyle(20);
  g_STcut_vs_MLQ->GetXaxis()->SetTitle("M_{LQ} (GeV/c^{2})");
  g_STcut_vs_MLQ->GetYaxis()->SetTitle("sT cut (GeV/c)");
  g_STcut_vs_MLQ->Draw("ap");

  TF1 *fitST = new TF1("fitST","pol2",minMLQfit,maxMLQfit);
  g_STcut_vs_MLQ->Fit(fitST);
  fitST->Draw("same");


  TLatex l;
  l.SetTextAlign(12);
  l.SetTextSize(0.04);
  l.SetTextFont(62);
  l.SetNDC();
  l.DrawLatex(0.15,0.75,"Optimized sT cut vs M_{LQ}"); 

  cout << endl;
  cout << "#########################" << endl;
  for(int mass=0; mass < Nmass_extrapolate; mass++)
    {


      cout << "The optimized sT cut for MLQ of "
	   << MLQ_extrapolate[mass] << " GeV/c2 is : " 
	   << fitST.Eval(MLQ_extrapolate[mass]) << " GeV/c" << endl;

    }
  cout << "#########################" << endl;

}//end macro
