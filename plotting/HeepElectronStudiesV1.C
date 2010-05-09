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
  legend->AddEntry(h[0],"HEEP electrons","l");
  legend->AddEntry(h[1],"ALL electrons","l");
  legend->Draw();

  TLine line0(cutValue0,0,cutValue0,h[1]->GetMaximum());
  line0.SetLineColor(kRed);
  line0.Draw("same");

  TLine line1(cutValue1,0,cutValue1,h[1]->GetMaximum());
  line1.SetLineColor(kRed);
  line1.Draw("same");
   
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
   string file0 = "output.root";
   // *** data ***
   string file1 = "output.root";
   //********************************************
   // make plots
   //********************************************

   TFile file0_(file0.c_str());
   TFile file1_(file1.c_str());

   //## barrel ##
   overlay_plots(file0, file1, "h_ElectronPt_barrel_heep", "h_ElectronPt_barrel_all", 0, 500, "ElectronPt [GeV]", "plots/h_ElectronPt_barrel.png",1, 0, 25, -999, "barrel");
   overlay_plots(file0, file1, "h_ElectronSCEta_fabs_barrel_heep", "h_ElectronSCEta_fabs_barrel_all", 0, 5, "ElectronSCEta_fabs", "plots/h_ElectronSCEta_fabs_barrel.png",0, 0, 1.442, -999,  "barrel");
   overlay_plots(file0, file1, "h_ElectronDeltaEtaTrkSC_barrel_heep", "h_ElectronDeltaEtaTrkSC_barrel_all", -0.05, 0.05, "ElectronDeltaEtaTrkSC", "plots/h_ElectronDeltaEtaTrkSC_barrel.png",1, 0, 0.005, -0.005,  "barrel");
   overlay_plots(file0, file1, "h_ElectronDeltaPhiTrkSC_barrel_heep", "h_ElectronDeltaPhiTrkSC_barrel_all", -0.5, 0.5, "ElectronDeltaPhiTrkSC", "plots/h_ElectronDeltaPhiTrkSC_barrel.png",1, 0, 0.09, -0.09,  "barrel");
   overlay_plots(file0, file1, "h_ElectronHoE_barrel_heep", "h_ElectronHoE_barrel_all", 0, 0.15, "ElectronHoE", "plots/h_ElectronHoE_barrel.png",1, 0, 0.05, -999,  "barrel");
   overlay_plots(file0, file1, "h_ElectronSigmaIEtaIEta_barrel_heep", "h_ElectronSigmaIEtaIEta_barrel_all", 0, 0.05, "ElectronSigmaIEtaIEta", "plots/h_ElectronSigmaIEtaIEta_barrel.png",1, 0, -999, -999,  "barrel");
   overlay_plots(file0, file1, "h_ElectronHcalIsoD2Heep_barrel_heep", "h_ElectronHcalIsoD2Heep_barrel_all", 0, 10, "ElectronHcalIsoD2Heep [GeV]", "plots/h_ElectronHcalIsoD2Heep_barrel.png",1, 0, -999, -999,  "barrel");
   overlay_plots(file0, file1, "h_ElectronTrkIsoHeep_barrel_heep", "h_ElectronTrkIsoHeep_barrel_all", 0, 40, "ElectronTrkIsoHeep [GeV]", "plots/h_ElectronTrkIsoHeep_barrel.png",1, 0, 7.5, -999,  "barrel");

   // ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt
   TCanvas c1_barrel;
   TH2F *h1_barrel[2];
   h1_barrel[0] = (TH2F*)file0_.Get("h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_barrel_all"); 
   h1_barrel[1] = (TH2F*)file0_.Get("h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_barrel_heep"); 

   h1_barrel[0]->SetStats(0);
   h1_barrel[0]->SetMarkerStyle(20);
   h1_barrel[0]->SetMarkerSize(0.4);
   h1_barrel[0]->GetYaxis()->SetRangeUser(0,10);
   h1_barrel[0]->GetXaxis()->SetRangeUser(0,200);
   h1_barrel[0]->GetXaxis()->SetTitle("ElectronPt [GeV]");   
   h1_barrel[0]->GetYaxis()->SetTitle("ElectronEcalIsoHeep_plus_HcalIsoD1Heep [GeV]");   
   h1_barrel[0]->Draw();
   h1_barrel[1]->SetMarkerStyle(20);
   h1_barrel[1]->SetMarkerSize(0.4);
   h1_barrel[1]->SetMarkerColor(kRed);
   h1_barrel[1]->SetLineColor(kRed);
   h1_barrel[1]->Draw("same");   
   TLegend *legend1_barrel = new TLegend(.4,.91,.75,.99);
   legend1_barrel->SetBorderSize(1);
   legend1_barrel->SetFillColor(0);
   legend1_barrel->AddEntry(h1_barrel[1],"HEEP electrons","l");
   legend1_barrel->AddEntry(h1_barrel[0],"ALL electrons","l");
   legend1_barrel->Draw();
   TLine line1_barrel(25,2+0.03*25, 200, 2+0.03*200);
   line1_barrel.SetLineColor(kRed);
   line1_barrel.Draw("same");
   TLatex l1_barrel;
   l1_barrel.SetTextAlign(12);
   l1_barrel.SetTextSize(0.04);
   l1_barrel.SetTextFont(62);
   l1_barrel.SetNDC();
   l1_barrel.DrawLatex(0.6,0.8,"barrel");
   c1_barrel.SaveAs("plots/h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_barrel.png");


   // ElectronE2x5OverE5x5_vs_E1x5OverE5x5 
   TCanvas c2_barrel;
   TH2F *h2_barrel[2];
   h2_barrel[0] = (TH2F*)file0_.Get("h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_barrel_all"); 
   h2_barrel[1] = (TH2F*)file0_.Get("h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_barrel_heep"); 

   h2_barrel[0]->SetStats(0);
   h2_barrel[0]->SetMarkerStyle(20);
   h2_barrel[0]->SetMarkerSize(0.4);
   h2_barrel[0]->GetYaxis()->SetRangeUser(0.8,1);
   h2_barrel[0]->GetXaxis()->SetRangeUser(0.7,1);
   h2_barrel[0]->GetXaxis()->SetTitle("E1x5OverE5x5");   
   h2_barrel[0]->GetYaxis()->SetTitle("ElectronE2x5OverE5x5");   
   h2_barrel[0]->Draw();
   h2_barrel[1]->SetMarkerStyle(20);
   h2_barrel[1]->SetMarkerSize(0.4);
   h2_barrel[1]->SetMarkerColor(kRed);
   h2_barrel[1]->SetLineColor(kRed);
   h2_barrel[1]->Draw("same");   
   TLegend *legend2_barrel = new TLegend(.4,.91,.75,.99);
   legend2_barrel->SetBorderSize(1);
   legend2_barrel->SetFillColor(0);
   legend2_barrel->AddEntry(h2_barrel[1],"HEEP electrons","l");
   legend2_barrel->AddEntry(h2_barrel[0],"ALL electrons","l");
   legend2_barrel->Draw();
   TLine line2_barrel(0.7,0.94, 1, 0.94);
   line2_barrel.SetLineColor(kRed);
   line2_barrel.Draw("same");
   TLine line2_barrel_(0.83,0.8, 0.83, 1);
   line2_barrel_.SetLineColor(kRed);
   line2_barrel_.Draw("same");
   TLatex l2_barrel;
   l2_barrel.SetTextAlign(12);
   l2_barrel.SetTextSize(0.04);
   l2_barrel.SetTextFont(62);
   l2_barrel.SetNDC();
   l2_barrel.DrawLatex(0.6,0.2,"barrel");
   c2_barrel.SaveAs("plots/h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_barrel.png");

   TCanvas c3_barrel;
   TH2F *h3_barrel;
   h3_barrel = (TH2F*)file0_.Get("h2_HEEPcheck_barrel"); 
   h3_barrel->Draw("colztext");
   l2_barrel.DrawLatex(0.6,0.2,"barrel");
   c3_barrel.SaveAs("plots/h2_HEEPcheck_barrel.png");
   

   //## endcap ##
   overlay_plots(file0, file1, "h_ElectronPt_endcap_heep", "h_ElectronPt_endcap_all", 0, 500, "ElectronPt [GeV]", "plots/h_ElectronPt_endcap.png",1, 0, 25, -999, "endcap");
   overlay_plots(file0, file1, "h_ElectronSCEta_fabs_endcap_heep", "h_ElectronSCEta_fabs_endcap_all", 0, 5, "ElectronSCEta_fabs", "plots/h_ElectronSCEta_fabs_endcap.png",0, 0, 1.560, 2.5,  "endcap");
   overlay_plots(file0, file1, "h_ElectronDeltaEtaTrkSC_endcap_heep", "h_ElectronDeltaEtaTrkSC_endcap_all", -0.05, 0.05, "ElectronDeltaEtaTrkSC", "plots/h_ElectronDeltaEtaTrkSC_endcap.png",1, 0, 0.007, -0.007,  "endcap");
   overlay_plots(file0, file1, "h_ElectronDeltaPhiTrkSC_endcap_heep", "h_ElectronDeltaPhiTrkSC_endcap_all", -0.5, 0.5, "ElectronDeltaPhiTrkSC", "plots/h_ElectronDeltaPhiTrkSC_endcap.png",1, 0, 0.09, -0.09,  "endcap");
   overlay_plots(file0, file1, "h_ElectronHoE_endcap_heep", "h_ElectronHoE_endcap_all", 0, 0.15, "ElectronHoE", "plots/h_ElectronHoE_endcap.png",1, 0, 0.05, -999,  "endcap");
   overlay_plots(file0, file1, "h_ElectronSigmaIEtaIEta_endcap_heep", "h_ElectronSigmaIEtaIEta_endcap_all", 0, 0.1, "ElectronSigmaIEtaIEta", "plots/h_ElectronSigmaIEtaIEta_endcap.png",1, 0, 0.03, -999,  "endcap");
   overlay_plots(file0, file1, "h_ElectronHcalIsoD2Heep_endcap_heep", "h_ElectronHcalIsoD2Heep_endcap_all", 0, 10, "ElectronHcalIsoD2Heep [GeV]", "plots/h_ElectronHcalIsoD2Heep_endcap.png",1, 0, 0.5, -999,  "endcap");
   overlay_plots(file0, file1, "h_ElectronTrkIsoHeep_endcap_heep", "h_ElectronTrkIsoHeep_endcap_all", 0, 40, "ElectronTrkIsoHeep [GeV]", "plots/h_ElectronTrkIsoHeep_endcap.png",1, 0, 15, -999,  "endcap");

   // ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt
   TCanvas c1_endcap;
   TH2F *h1_endcap[2];
   h1_endcap[0] = (TH2F*)file0_.Get("h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_endcap_all"); 
   h1_endcap[1] = (TH2F*)file0_.Get("h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_endcap_heep"); 

   h1_endcap[0]->SetStats(0);
   h1_endcap[0]->SetMarkerStyle(20);
   h1_endcap[0]->SetMarkerSize(0.4);
   h1_endcap[0]->GetYaxis()->SetRangeUser(0,10);
   h1_endcap[0]->GetXaxis()->SetRangeUser(0,200);
   h1_endcap[0]->GetXaxis()->SetTitle("ElectronPt [GeV]");   
   h1_endcap[0]->GetYaxis()->SetTitle("ElectronEcalIsoHeep_plus_HcalIsoD1Heep [GeV]");   
   h1_endcap[0]->Draw();
   h1_endcap[1]->SetMarkerStyle(20);
   h1_endcap[1]->SetMarkerSize(0.4);
   h1_endcap[1]->SetMarkerColor(kRed);
   h1_endcap[1]->SetLineColor(kRed);
   h1_endcap[1]->Draw("same");   
   TLegend *legend1_endcap = new TLegend(.4,.91,.75,.99);
   legend1_endcap->SetBorderSize(1);
   legend1_endcap->SetFillColor(0);
   legend1_endcap->AddEntry(h1_endcap[1],"HEEP electrons","l");
   legend1_endcap->AddEntry(h1_endcap[0],"ALL electrons","l");
   legend1_endcap->Draw();
   TLine line1_endcap(0,2.5, 50, 2.5);
   line1_endcap.SetLineColor(kRed);
   line1_endcap.Draw("same");
   TLine line1_endcap_(50,2.5, 200, 2.5+0.03*(200-50));
   line1_endcap_.SetLineColor(kRed);
   line1_endcap_.Draw("same");
   TLatex l1_endcap;
   l1_endcap.SetTextAlign(12);
   l1_endcap.SetTextSize(0.04);
   l1_endcap.SetTextFont(62);
   l1_endcap.SetNDC();
   l1_endcap.DrawLatex(0.6,0.8,"endcap");
   c1_endcap.SaveAs("plots/h2_ElectronEcalIsoHeep_plus_HcalIsoD1Heep_vs_ElectronPt_endcap.png");


   // ElectronE2x5OverE5x5_vs_E1x5OverE5x5 
   TCanvas c2_endcap;
   TH2F *h2_endcap[2];
   h2_endcap[0] = (TH2F*)file0_.Get("h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_endcap_all"); 
   h2_endcap[1] = (TH2F*)file0_.Get("h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_endcap_heep"); 

   h2_endcap[0]->SetStats(0);
   h2_endcap[0]->SetMarkerStyle(20);
   h2_endcap[0]->SetMarkerSize(0.4);
   h2_endcap[0]->GetYaxis()->SetRangeUser(0.8,1);
   h2_endcap[0]->GetXaxis()->SetRangeUser(0.7,1);
   h2_endcap[0]->GetXaxis()->SetTitle("E1x5OverE5x5");   
   h2_endcap[0]->GetYaxis()->SetTitle("ElectronE2x5OverE5x5");   
   h2_endcap[0]->Draw();
   h2_endcap[1]->SetMarkerStyle(20);
   h2_endcap[1]->SetMarkerSize(0.4);
   h2_endcap[1]->SetMarkerColor(kRed);
   h2_endcap[1]->SetLineColor(kRed);
   h2_endcap[1]->Draw("same");   
   TLegend *legend2_endcap = new TLegend(.4,.91,.75,.99);
   legend2_endcap->SetBorderSize(1);
   legend2_endcap->SetFillColor(0);
   legend2_endcap->AddEntry(h2_endcap[1],"HEEP electrons","l");
   legend2_endcap->AddEntry(h2_endcap[0],"ALL electrons","l");
   legend2_endcap->Draw();
   //    TLine line2_endcap(0.7,0.94, 1, 0.94);
   //    line2_endcap.SetLineColor(kRed);
   //    line2_endcap.Draw("same");
   //    TLine line2_endcap_(0.83,0.8, 0.83, 1);
   //    line2_endcap_.SetLineColor(kRed);
   //    line2_endcap_.Draw("same");
   TLatex l2_endcap;
   l2_endcap.SetTextAlign(12);
   l2_endcap.SetTextSize(0.04);
   l2_endcap.SetTextFont(62);
   l2_endcap.SetNDC();
   l2_endcap.DrawLatex(0.6,0.2,"endcap");
   c2_endcap.SaveAs("plots/h2_ElectronE2x5OverE5x5_vs_E1x5OverE5x5_endcap.png");

   TCanvas c3_endcap;
   TH2F *h3_endcap;
   h3_endcap = (TH2F*)file0_.Get("h2_HEEPcheck_endcap"); 
   h3_endcap->Draw("colztext");
   l2_endcap.DrawLatex(0.6,0.2,"endcap");
   c3_endcap.SaveAs("plots/h2_HEEPcheck_endcap.png");

}
