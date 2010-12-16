void style() {

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
     
//    gStyle->SetOptStat(kFALSE);
   gStyle->SetOptTitle(0); // for PAS/AN
//    gStyle->SetPadLeftMargin(0.13); // for Noise-only
//    gStyle->SetPadRightMargin(0.07); // for Noise-only
//    gStyle->SetStatFont(42);
//    gStyle->SetTitleFont(42);
//    gStyle->SetTitleFont(42, "XYZ");
//    gStyle->SetLabelFont(42, "XYZ");
}

void overlay_plots(const string& fFile0, const string& fFile1, const string& fPlot0, const string& fPlot1, const double fXmin, const double fXmax, const string& fXAxisLabel, const string& fName, const int logY, const int rescale, const double cutValue0, const double cutValue1, const string& label) {
  
  TH1F *h[2];
 
  TFile file0(fFile0.c_str());
  h[0] = (TH1F*)file0.Get(fPlot0.c_str()); 

  TFile file1(fFile1.c_str());
  h[1] = (TH1F*)file1.Get(fPlot1.c_str()); 
   

  h[0]->Rebin(5);
  h[1]->Rebin(5);

  //    double ymax = (h[0]->GetMaximum()>h[1]->GetMaximum()) ? h[0]->GetMaximum() : h[1]->GetMaximum();
   
  h[1]->GetXaxis()->SetTitle(fXAxisLabel.c_str());
  h[1]->GetXaxis()->SetRangeUser(fXmin,fXmax);
  //    h[0]->GetYaxis()->SetRangeUser(0.5,1.25*ymax);

  h[1]->SetTitleOffset(1.2,"X");
  h[1]->GetXaxis()->SetTitleSize(0.04);
  h[1]->GetYaxis()->SetTitleSize(0.04);

  if(rescale) 
    {  
      double scale = h[1]->GetEntries()/h[0]->GetEntries();
      h[0]->Scale(scale);
    }

  TCanvas *c = new TCanvas("c","",800,800);
  c->cd();

  h[1]->SetLineWidth(3);
  //    h[1]->SetLineStyle(3);
  h[1]->SetLineColor(kBlack);
  //    h[1]->SetMarkerSize(.8);
  h[1]->SetMarkerStyle(20);
  h[1]->SetMarkerColor(kBlack);
  h[1]->Draw();
  
  h[0]->SetLineWidth(3);
  //    h[0]->SetLineStyle(4);
  h[0]->SetLineColor(kRed);
  h[0]->SetFillColor(kRed);
  h[0]->SetFillStyle(3002);
  //    h[0]->SetMarkerSize(.6);
  //    h[0]->SetMarkerStyle(26);
  //    h[0]->SetMarkerColor(kRed);
  h[0]->Draw("sameshist");
   
  //update the current pad, needed to modify statboxes
  gPad->Update();
   
  // get the statboxes and set color
  TPaveStats *st1 = (TPaveStats*)h[0]->GetListOfFunctions()->FindObject("stats");
  st1->SetTextColor(kRed);
    st1->SetLineColor(kRed);
  st1->SetOptStat(1111111);
  TPaveStats *st2 = (TPaveStats*)h[1]->GetListOfFunctions()->FindObject("stats");
  st2->SetTextColor(kBlack);
  st2->SetLineColor(kBlack);
  st2->SetOptStat(1111111);

  // set the position of the statboxes
  double x1 = st1->GetX1NDC();
  double y1 = st1->GetY1NDC();
  double x2 = st1->GetX2NDC();
  double y2 = st1->GetY2NDC();
  //double xx = x2-x1;
  double yy = y2-y1;
  st2->SetX1NDC(x1);
  st2->SetY1NDC(y1-yy);
  st2->SetX2NDC(x2);
  st2->SetY2NDC(y1);
  gPad->Modified();
   
  TLegend *legend = new TLegend(.4,.91,.75,.99);
  legend->SetBorderSize(1);
  legend->SetFillColor(0);
  //    legend->SetFillStyle(0);
  legend->AddEntry(h[0],"ALL MC background","l");
  legend->AddEntry(h[1],"DATA","l");
  legend->Draw();

  //   TLine line0(cutValue0,0,cutValue0,h[1]->GetMaximum());
  //   line0.SetLineColor(kRed);
  //   line0.Draw("same");
  
  //   TLine line1(cutValue1,0,cutValue1,h[1]->GetMaximum());
  //   line1.SetLineColor(kRed);
  //   line1.Draw("same");
   
  TLatex l;
  l.SetTextAlign(12);
  l.SetTextSize(0.04);
  l.SetTextFont(62);
  l.SetNDC();
  l.DrawLatex(0.6,0.8,label.c_str());

  if(logY==1)
    c->SetLogy();
   
  string fileName = fName;
  c->SaveAs(fileName.c_str());
   
  delete legend;
  delete c;
}


void makePlots() {
   // turn on/off batch mode
   gROOT->SetBatch(kTRUE);
   // set ROOT style
   style();
   //********************************************
   // root files
   //********************************************
   // *** MC ***
   string file0 = "/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/33.2pb-1_sT_presel_250_Zrescale1.20_Wrescale1.06_WithEleIDIsoVar_pcuscms46_noLQ//analysisClass_enujjSample_plots.root";
   // *** data ***
   string file1 = "/home/santanas/Leptoquarks/data/output_fromAFS/enujj_analysis/33.2pb-1_sT_presel_250_Zrescale1.20_Wrescale1.06_WithEleIDIsoVar_pcuscms46_noLQ//analysisClass_enujjSample_plots.root";
   //********************************************
   // make plots
   //********************************************

   TFile file0_(file0.c_str());
   TFile file1_(file1.c_str());

   //## barrel ##
   overlay_plots(file0, file1, "histo1D__ALLBKG__h1_ElectronDeltaEtaTrkSC_highMT", "histo1D__DATA__h1_ElectronDeltaEtaTrkSC_highMT", -0.01, 0.01, "ElectronDeltaEtaTrkSC", "h_ElectronDeltaEtaTrkSC.png",0, 1, 0.005, -0.005,  "barrel+endcap");
   overlay_plots(file0, file1, "histo1D__ALLBKG__h1_ElectronDeltaPhiTrkSC_highMT", "histo1D__DATA__h1_ElectronDeltaPhiTrkSC_highMT", -0.1, 0.1, "ElectronDeltaPhiTrkSC", "h_ElectronDeltaPhiTrkSC.png",0, 1, 0.005, -0.005,  "barrel+endcap");
   overlay_plots(file0, file1, "histo1D__ALLBKG__h1_ElectronHoE_highMT", "histo1D__DATA__h1_ElectronHoE_highMT", 0, 0.1, "ElectronHoE", "h_ElectronHoE.png",0, 1, 0.005, -0.005,  "barrel+endcap");
   overlay_plots(file0, file1, "histo1D__ALLBKG__h1_ElectronSigmaIEtaIEta_highMT", "histo1D__DATA__h1_ElectronSigmaIEtaIEta_highMT", 0, 0.1, "ElectronSigmaIEtaIEta", "h_ElectronSigmaIEtaIEta.png",0, 1, 0.005, -0.005,  "barrel+endcap");
   overlay_plots(file0, file1, "histo1D__ALLBKG__h1_ElectronEcalHcalIsoHeep_highMT", "histo1D__DATA__h1_ElectronEcalHcalIsoHeep_highMT", 0, 10, "ElectronEcalHcalIsoHeep", "h_ElectronEcalHcalIsoHeep.png",0, 1, 0.005, -0.005,  "barrel+endcap");
   overlay_plots(file0, file1, "histo1D__ALLBKG__h1_ElectronHcalIsoD2Heep_highMT", "histo1D__DATA__h1_ElectronHcalIsoD2Heep_highMT", 0, 10, "ElectronHcalIsoD2Heep", "h_ElectronHcalIsoD2Heep.png",0, 1, 0.005, -0.005,  "barrel+endcap");
   overlay_plots(file0, file1, "histo1D__ALLBKG__h1_ElectronTrkIsoHeep_highMT", "histo1D__DATA__h1_ElectronTrkIsoHeep_highMT", 0, 10, "ElectronTrkIsoHeep", "h_ElectronTrkIsoHeep.png",0, 1, 0.005, -0.005,  "barrel+endcap");
   overlay_plots(file0, file1, "histo1D__ALLBKG__h1_ElectronE2x5OverE5x5_highMT", "histo1D__DATA__h1_ElectronE2x5OverE5x5_highMT", 0, 1.2, "ElectronE2x5OverE5x5", "h_ElectronE2x5OverE5x5.png",0, 1, 0.005, -0.005,  "barrel+endcap");
   overlay_plots(file0, file1, "histo1D__ALLBKG__h1_ElectronE1x5OverE5x5_highMT", "histo1D__DATA__h1_ElectronE1x5OverE5x5_highMT", 0, 1.2, "ElectronE1x5OverE5x5", "h_ElectronE1x5OverE5x5.png",0, 1, 0.005, -0.005,  "barrel+endcap");


}
