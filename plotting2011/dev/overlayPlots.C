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

void overlay_plots(const string& fFile0, const string& fFile1, const string& fFile2, const string& fFile3, const string& fFile4, const string& fPlot0, const string& fPlot1, const string& fPlot2, const string& fPlot3, const string& fPlot4, const double fXmin, const double fXmax, const string& fXAxisLabel, const string& fName, const int logY, const int rescale, const double cutValue0, const double cutValue1, const string& label) {
  
  TH1F *h[5];
 
  TFile file0(fFile0.c_str());
  h[0] = (TH1F*)file0.Get(fPlot0.c_str()); 

  TFile file1(fFile1.c_str());
  h[1] = (TH1F*)file1.Get(fPlot1.c_str()); 

  TFile file2(fFile2.c_str());
  h[2] = (TH1F*)file2.Get(fPlot2.c_str()); 

  TFile file3(fFile3.c_str());
  h[3] = (TH1F*)file3.Get(fPlot3.c_str()); 
   
  TFile file4(fFile4.c_str());
  h[4] = (TH1F*)file4.Get(fPlot4.c_str()); 

  h[0]->Rebin(1);
  h[1]->Rebin(1);
  h[2]->Rebin(1);
  h[3]->Rebin(1);
  h[4]->Rebin(1);

  //    double ymax = (h[0]->GetMaximum()>h[1]->GetMaximum()) ? h[0]->GetMaximum() : h[1]->GetMaximum();
   
  h[0]->GetXaxis()->SetTitle(fXAxisLabel.c_str());
  h[0]->GetXaxis()->SetRangeUser(fXmin,fXmax);
  //    h[0]->GetYaxis()->SetRangeUser(0.5,1.25*ymax);

  h[0]->SetTitleOffset(1.2,"X");
  h[0]->GetXaxis()->SetTitleSize(0.04);
  h[0]->GetYaxis()->SetTitleSize(0.04);

  if(rescale) 
    {  

       if (h[0]->Integral()!=0)
	h[0]->Scale(1/h[0]->Integral());

       //----------------------------       
//        if (h[1]->Integral()!=0)
// 	 h[1]->Scale(1/h[1]->Integral());

//        if (h[2]->Integral()!=0)
// 	 h[2]->Scale(1/h[2]->Integral());
       //----------------------------       

       if (h[3]->Integral()!=0)
	 h[3]->Scale(1/h[3]->Integral());
 
       if (h[4]->Integral()!=0)
	 h[4]->Scale(1/h[4]->Integral());

    }

  TCanvas *c = new TCanvas("c","",800,800);
  c->cd();
  c->Range(-0.4638998,-0.01028103,4.210331,0.09111118);
  c->SetRightMargin(0.2286432);

  h[0]->SetLineWidth(3);
  //    h[0]->SetLineStyle(4);
  h[0]->SetLineColor(kRed);
  h[0]->SetFillColor(kRed);
  h[0]->SetFillStyle(3002);
  //    h[0]->SetMarkerSize(.6);
  //    h[0]->SetMarkerStyle(26);
  //    h[0]->SetMarkerColor(kRed);
  h[0]->Draw("hist");

  h[1]->SetLineWidth(3);
  //    h[1]->SetLineStyle(3);
  h[1]->SetLineColor(kBlack);
  //    h[1]->SetMarkerSize(.8);
  h[1]->SetMarkerStyle(20);
  h[1]->SetMarkerColor(kBlack);
  //----------------------------       
  //h[1]->Draw("sameshist");
  //----------------------------       

  h[2]->SetLineWidth(3);
  h[2]->SetLineStyle(3);
  h[2]->SetLineColor(kBlack);
  //    h[1]->SetMarkerSize(.8);
  h[2]->SetMarkerStyle(20);
  h[2]->SetMarkerColor(kBlack);
  //----------------------------       
  //h[2]->Draw("sameshist");
  //----------------------------       

  h[3]->SetLineWidth(3);
  h[3]->SetLineStyle(4);
  h[3]->SetLineColor(kBlack);
  //    h[1]->SetMarkerSize(.8);
  h[3]->SetMarkerStyle(20);
  h[3]->SetMarkerColor(kBlack);
  h[3]->Draw("sameshist");

  h[4]->SetLineWidth(3);
  h[4]->SetLineStyle(5);
  h[4]->SetLineColor(kBlack);
  //    h[1]->SetMarkerSize(.8);
  h[4]->SetMarkerStyle(20);
  h[4]->SetMarkerColor(kBlack);
  h[4]->Draw("sameshist");
    
  double maximum = 0;
  if( h[0]->GetMaximum() > maximum )
    maximum = h[0]->GetMaximum() ;
  //-------------------------------
  //   if( h[1]->GetMaximum() > maximum )
  //     maximum = h[1]->GetMaximum();
  //   if( h[2]->GetMaximum() > maximum )
  //     maximum = h[2]->GetMaximum();
  //-------------------------------
  if( h[3]->GetMaximum() > maximum )
    maximum = h[3]->GetMaximum();
  if( h[4]->GetMaximum() > maximum )
    maximum = h[4]->GetMaximum();
  
  h[0]->GetYaxis()->SetRangeUser(0,maximum+maximum*0.1);

  //update the current pad, needed to modify statboxes
  gPad->Update();
   
  // get the statboxes and set color
  TPaveStats *st1 = (TPaveStats*)h[0]->GetListOfFunctions()->FindObject("stats");
  st1->SetTextColor(kRed);
    st1->SetLineColor(kRed);
  st1->SetOptStat(1111111);
  //------------------------------
  //   TPaveStats *st2 = (TPaveStats*)h[1]->GetListOfFunctions()->FindObject("stats");
  //   st2->SetTextColor(kBlack);
  //   st2->SetLineColor(kBlack);
  //   st2->SetOptStat(1111111);
  //   TPaveStats *st3 = (TPaveStats*)h[2]->GetListOfFunctions()->FindObject("stats");
  //   st3->SetTextColor(kBlack);
  //   st3->SetLineColor(kBlack);
  //   st3->SetOptStat(1111111);
  //------------------------------
  TPaveStats *st4 = (TPaveStats*)h[3]->GetListOfFunctions()->FindObject("stats");
  st4->SetTextColor(kBlack);
  st4->SetLineColor(kBlack);
  st4->SetOptStat(1111111);
  TPaveStats *st5 = (TPaveStats*)h[4]->GetListOfFunctions()->FindObject("stats");
  st5->SetTextColor(kBlack);
  st5->SetLineColor(kBlack);
  st5->SetOptStat(1111111);


  // set the position of the statboxes
  double x1 = st1->GetX1NDC();
  double y1 = st1->GetY1NDC();
  double x2 = st1->GetX2NDC();
  double y2 = st1->GetY2NDC();
  //double xx = x2-x1;
  double yy = y2-y1;
  //--------------------------
  //   st2->SetX1NDC(x1);
  //   st2->SetY1NDC(y1-yy);
  //   st2->SetX2NDC(x2);
  //   st2->SetY2NDC(y1);
  //   st3->SetX1NDC(x1);
  //   st3->SetY1NDC(y1-yy-yy);
  //   st3->SetX2NDC(x2);
  //   st3->SetY2NDC(y1-yy);
  //--------------------------
  st4->SetX1NDC(x1);
  st4->SetY1NDC(y1-yy-yy-yy);
  st4->SetX2NDC(x2);
  st4->SetY2NDC(y1-yy-yy);
  st5->SetX1NDC(x1);
  st5->SetY1NDC(y1-yy-yy-yy-yy);
  st5->SetX2NDC(x2);
  st5->SetY2NDC(y1-yy-yy-yy);
  gPad->Modified();
   
  TLegend *legend = new TLegend(.4,.91,.75,.99);
  legend->SetBorderSize(1);
  legend->SetFillColor(0);
  //    legend->SetFillStyle(0);
  legend->AddEntry(h[0],"SM background","l");
  //--------------------------
  //   legend->AddEntry(h[1],"Axigluon + W , M = 150 GeV","l");
  //   legend->AddEntry(h[2],"Axigluon + W , M = 500 GeV","l");
  //--------------------------
  legend->AddEntry(h[3],"Axigluon + W , M = 1000 GeV","l");
  legend->AddEntry(h[4],"Axigluon + W , M = 1500 GeV","l");
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
   //string file0 = "$AXIGLUONDATA/axigluons_enujj/5fb-1_Summer11MC_AxigluonW_enujj_13092011/output_cutTable_axigluons_enujj/analysisClass_axigluons_enujj_plots.root";
   string file0 = "$AXIGLUONDATA/axigluons_enujj/5fb-1_Summer11MC_AxigluonW_enujj_Mjj500cut_13092011_2/output_cutTable_axigluons_enujj/analysisClass_axigluons_enujj_plots.root";
   // *** signal ***
   string file1 = file0;
   string file2 = file0;
   string file3 = file0;
   string file4 = file0;
//    string file1 = "$AXIGLUONDATA/axigluons_lnujj/AxigluonW_GenStudies_13092011/AG150WGen.root";
//    string file2 = "$AXIGLUONDATA/axigluons_lnujj/AxigluonW_GenStudies_13092011/AG500WGen.root";
//    string file3 = "$AXIGLUONDATA/axigluons_lnujj/AxigluonW_GenStudies_13092011/AG1000WGen.root";
//    string file4 = "$AXIGLUONDATA/axigluons_lnujj/AxigluonW_GenStudies_13092011/AG1500WGen.root";

   //********************************************
   // make plots
   //********************************************

   TFile file0_(file0.c_str());
   TFile file1_(file1.c_str());

   
   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__Ele1_Pt" , 
		 "histo1D__AxigluonW_M150__Ele1_Pt", "histo1D__AxigluonW_M500__Ele1_Pt", "histo1D__AxigluonW_M1000__Ele1_Pt", "histo1D__AxigluonW_M1500__Ele1_Pt" ,
		 0, 600, "Ele1 Pt (GeV)", "Ele1_Pt.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__Ele1_Eta" , 
		 "histo1D__AxigluonW_M150__Ele1_Eta", "histo1D__AxigluonW_M500__Ele1_Eta", "histo1D__AxigluonW_M1000__Ele1_Eta", "histo1D__AxigluonW_M1500__Ele1_Eta" ,
		 -3, 3, "Ele1 \\eta (GeV)", "Ele1_Eta.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__Jet1_Pt" , 
		 "histo1D__AxigluonW_M150__Jet1_Pt", "histo1D__AxigluonW_M500__Jet1_Pt", "histo1D__AxigluonW_M1000__Jet1_Pt", "histo1D__AxigluonW_M1500__Jet1_Pt" ,
		 0, 1000, "Jet1 Pt (GeV)", "Jet1_Pt.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__Jet2_Pt" , 
		 "histo1D__AxigluonW_M150__Jet2_Pt", "histo1D__AxigluonW_M500__Jet2_Pt", "histo1D__AxigluonW_M1000__Jet2_Pt", "histo1D__AxigluonW_M1500__Jet2_Pt" ,
		 0, 1000, "Jet2 Pt (GeV)", "Jet2_Pt.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__Jet1_Eta" , 
		 "histo1D__AxigluonW_M150__Jet1_Eta", "histo1D__AxigluonW_M500__Jet1_Eta", "histo1D__AxigluonW_M1000__Jet1_Eta", "histo1D__AxigluonW_M1500__Jet1_Eta" ,
		 -3, 3, "Jet1 \\eta (GeV)", "Jet1_Eta.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__Jet2_Eta" , 
		 "histo1D__AxigluonW_M150__Jet2_Eta", "histo1D__AxigluonW_M500__Jet2_Eta", "histo1D__AxigluonW_M1000__Jet2_Eta", "histo1D__AxigluonW_M1500__Jet2_Eta" ,
		 -3, 3, "Jet2 \\eta (GeV)", "Jet2_Eta.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__mDEta_Jet1Jet2" , 
		 "histo1D__AxigluonW_M150__mDEta_Jet1Jet2", "histo1D__AxigluonW_M500__mDEta_Jet1Jet2", "histo1D__AxigluonW_M1000__mDEta_Jet1Jet2", "histo1D__AxigluonW_M1500__mDEta_Jet1Jet2" ,
		 0, 5, "|\\Delta\\eta (j1,j2)|", "mDEta_Jet1Jet2.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__mDPhi_Jet1Jet2" , 
		 "histo1D__AxigluonW_M150__mDPhi_Jet1Jet2", "histo1D__AxigluonW_M500__mDPhi_Jet1Jet2", "histo1D__AxigluonW_M1000__mDPhi_Jet1Jet2", "histo1D__AxigluonW_M1500__mDPhi_Jet1Jet2" ,
		 0, 3.1416, "|\\Delta\\phi (j1,j2)|", "mDPhi_Jet1Jet2.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__mDphi_BosonJet1" , 
		 "histo1D__AxigluonW_M150__mDphi_BosonJet1", "histo1D__AxigluonW_M500__mDphi_BosonJet1", "histo1D__AxigluonW_M1000__mDphi_BosonJet1", "histo1D__AxigluonW_M1500__mDphi_BosonJet1" ,
		 0, 3.1416, "|\\Delta\\phi (boson,j1)|", "mDphi_BosonJet1.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__mDphi_BosonJet2" , 
		 "histo1D__AxigluonW_M150__mDphi_BosonJet2", "histo1D__AxigluonW_M500__mDphi_BosonJet2", "histo1D__AxigluonW_M1000__mDphi_BosonJet2", "histo1D__AxigluonW_M1500__mDphi_BosonJet2" ,
		 0, 3.1416, "|\\Delta\\phi (boson,j2)|", "mDphi_BosonJet2.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__Jet1_Pt_over_Mj1j2" , 
		 "histo1D__AxigluonW_M150__Jet1_Pt_over_Mj1j2", "histo1D__AxigluonW_M500__Jet1_Pt_over_Mj1j2", "histo1D__AxigluonW_M1000__Jet1_Pt_over_Mj1j2", "histo1D__AxigluonW_M1500__Jet1_Pt_over_Mj1j2" ,
		 0, 10, "Jet1 Pt / M_j1j2", "Jet1_Pt_over_Mj1j2.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__Jet2_Pt_over_Mj1j2" , 
		 "histo1D__AxigluonW_M150__Jet2_Pt_over_Mj1j2", "histo1D__AxigluonW_M500__Jet2_Pt_over_Mj1j2", "histo1D__AxigluonW_M1000__Jet2_Pt_over_Mj1j2", "histo1D__AxigluonW_M1500__Jet2_Pt_over_Mj1j2" ,
		 0, 10, "Jet2 Pt / M_j1j2", "Jet2_Pt_over_Mj1j2.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__M_j1j2" , 
		 "histo1D__AxigluonW_M150__M_j1j2", "histo1D__AxigluonW_M500__M_j1j2", "histo1D__AxigluonW_M1000__M_j1j2", "histo1D__AxigluonW_M1500__M_j1j2" ,
		 0, 2000, "M_j1j2 (GeV)", "M_j1j2.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__MT_Ele1MET" , 
		 "histo1D__AxigluonW_M150__MT_Ele1MET", "histo1D__AxigluonW_M500__MT_Ele1MET", "histo1D__AxigluonW_M1000__MT_Ele1MET", "histo1D__AxigluonW_M1500__MT_Ele1MET" ,
		 0, 500, "M_{T}(e,MET) (GeV)", "MT_Ele1MET.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__MET_Pt" , 
		 "histo1D__AxigluonW_M150__MET_Pt", "histo1D__AxigluonW_M500__MET_Pt", "histo1D__AxigluonW_M1000__MET_Pt", "histo1D__AxigluonW_M1500__MET_Pt" ,
		 0, 600, "PFMET (GeV)", "MET_Pt.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__Pt_j1j2" , 
		 "histo1D__AxigluonW_M150__Pt_j1j2", "histo1D__AxigluonW_M500__Pt_j1j2", "histo1D__AxigluonW_M1000__Pt_j1j2", "histo1D__AxigluonW_M1500__Pt_j1j2" ,
		 0, 600, "Pt(j1,j2) (GeV)", "Pt_j1j2.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__Jet1_E_over_Mj1j2" , 
		 "histo1D__AxigluonW_M150__Jet1_E_over_Mj1j2", "histo1D__AxigluonW_M500__Jet1_E_over_Mj1j2", "histo1D__AxigluonW_M1000__Jet1_E_over_Mj1j2", "histo1D__AxigluonW_M1500__Jet1_E_over_Mj1j2" ,
		 0, 10, "Jet1 Energy / M_j1j2", "Jet1_E_over_Mj1j2.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__Jet1_E_over_Jet1_Pt" , 
		 "histo1D__AxigluonW_M150__Jet1_E_over_Jet1_Pt", "histo1D__AxigluonW_M500__Jet1_E_over_Jet1_Pt", "histo1D__AxigluonW_M1000__Jet1_E_over_Jet1_Pt", "histo1D__AxigluonW_M1500__Jet1_E_over_Jet1_Pt" ,
		 0, 10, "Jet1 Energy / Jet1 Pt", "Jet1_E_over_Jet1_Pt.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__Jet2_E_over_Mj1j2" , 
		 "histo1D__AxigluonW_M150__Jet2_E_over_Mj1j2", "histo1D__AxigluonW_M500__Jet2_E_over_Mj1j2", "histo1D__AxigluonW_M1000__Jet2_E_over_Mj1j2", "histo1D__AxigluonW_M1500__Jet2_E_over_Mj1j2" ,
		 0, 10, "Jet2 Energy / M_j1j2", "Jet2_E_over_Mj1j2.png",0, 1, -99, -99,  "");

   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__Jet2_E_over_Jet2_Pt" , 
		 "histo1D__AxigluonW_M150__Jet2_E_over_Jet2_Pt", "histo1D__AxigluonW_M500__Jet2_E_over_Jet2_Pt", "histo1D__AxigluonW_M1000__Jet2_E_over_Jet2_Pt", "histo1D__AxigluonW_M1500__Jet2_E_over_Jet2_Pt" ,
		 0, 10, "Jet2 Energy / Jet2 Pt", "Jet2_E_over_Jet2_Pt.png",0, 1, -99, -99,  "");


  
   

   /* 
   overlay_plots(file0, file1, file2, file3, file4, "histo1D__ALLBKG__M_j1j2" , 
		 "h1_Mass_AG", "h1_Mass_AG", "h1_Mass_AG", "h1_Mass_AG" ,
		 0, 2000, "Axigluon Mass (GeV)", "M_j1j2.png",0, 1, -99, -99,  "");

   */



}
