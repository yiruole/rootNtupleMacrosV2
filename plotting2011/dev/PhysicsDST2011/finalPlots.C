{
  gROOT->Reset();
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit(1111);

  //DATA - HLT
  TFile data_hlt ("$DIJETDATA/dijets_PhysicsDST/117pb-1_JECL123Res__Fall11MC_JECL123__31_01_2012/finalResults_DATA.root");

  //MC - HLT
  TFile QstarToJJ_M_500_qg_hlt ("$DIJETDATA/dijets_PhysicsDST/117pb-1_JECL123Res__Fall11MC_JECL123__31_01_2012/finalResults_QstarToJJ_M-500_qg.root");
  TFile QstarToJJ_M_700_qg_hlt ("$DIJETDATA/dijets_PhysicsDST/117pb-1_JECL123Res__Fall11MC_JECL123__31_01_2012/finalResults_QstarToJJ_M-700_qg.root");
  TFile QstarToJJ_M_1200_qg_hlt ("$DIJETDATA/dijets_PhysicsDST/117pb-1_JECL123Res__Fall11MC_JECL123__31_01_2012/finalResults_QstarToJJ_M-1200_qg.root");
  TFile QstarToJJ_M_2000_qg_hlt ("$DIJETDATA/dijets_PhysicsDST/117pb-1_JECL123Res__Fall11MC_JECL123__31_01_2012/finalResults_QstarToJJ_M-2000_qg.root");

  TFile RSGravitonToJJ_M_700_gg_hlt ("$DIJETDATA/dijets_PhysicsDST/117pb-1_JECL123Res__Fall11MC_JECL123__31_01_2012/finalResults_RSGravitonToJJ_M-700_gg.root");
  TFile RSGravitonToJJ_M_1200_gg_hlt ("$DIJETDATA/dijets_PhysicsDST/117pb-1_JECL123Res__Fall11MC_JECL123__31_01_2012/finalResults_RSGravitonToJJ_M-1200_gg.root");
  TFile RSGravitonToJJ_M_2000_gg_hlt ("$DIJETDATA/dijets_PhysicsDST/117pb-1_JECL123Res__Fall11MC_JECL123__31_01_2012/finalResults_RSGravitonToJJ_M-2000_gg.root");

  TFile RSGravitonToJJ_M_700_qq_hlt ("$DIJETDATA/dijets_PhysicsDST/117pb-1_JECL123Res__Fall11MC_JECL123__31_01_2012/finalResults_RSGravitonToJJ_M-700_qq.root");
  TFile RSGravitonToJJ_M_1200_qq_hlt ("$DIJETDATA/dijets_PhysicsDST/117pb-1_JECL123Res__Fall11MC_JECL123__31_01_2012/finalResults_RSGravitonToJJ_M-1200_qq.root");
  TFile RSGravitonToJJ_M_2000_qq_hlt ("$DIJETDATA/dijets_PhysicsDST/117pb-1_JECL123Res__Fall11MC_JECL123__31_01_2012/finalResults_RSGravitonToJJ_M-2000_qq.root");

  //Ranges
  double XMAX = 4000;
  double YMAX = 0.3;
  double min_mjj_FatPF = 526;

  //DATA - FatJets
  TH1D *h_mjj_FatPF_VarBin_data_hlt = (TH1D*)data_hlt.Get( "M_FatPFJet1FatPFJet2_VarBin" ); 

  //DATA - PFJets
  TH1D *h_mjj_PF_VarBin_data_hlt = (TH1D*)data_hlt.Get( "M_PFJet1PFJet2_VarBin" ); 

  //MC - FatJets
  TH1D *h_mjj_FatPF_VarBin_QstarToJJ_M_500_qg_hlt = (TH1D*)QstarToJJ_M_500_qg_hlt.Get( "M_FatPFJet1FatPFJet2_VarBin" ); 
  TH1D *h_mjj_FatPF_VarBin_QstarToJJ_M_700_qg_hlt = (TH1D*)QstarToJJ_M_700_qg_hlt.Get( "M_FatPFJet1FatPFJet2_VarBin" ); 
  TH1D *h_mjj_FatPF_VarBin_QstarToJJ_M_1200_qg_hlt = (TH1D*)QstarToJJ_M_1200_qg_hlt.Get( "M_FatPFJet1FatPFJet2_VarBin" ); 
  TH1D *h_mjj_FatPF_VarBin_QstarToJJ_M_2000_qg_hlt = (TH1D*)QstarToJJ_M_2000_qg_hlt.Get( "M_FatPFJet1FatPFJet2_VarBin" ); 

  TH1D *h_mjj_FatPF_VarBin_RSGravitonToJJ_M_700_gg_hlt = (TH1D*)RSGravitonToJJ_M_700_gg_hlt.Get( "M_FatPFJet1FatPFJet2_VarBin" ); 
  TH1D *h_mjj_FatPF_VarBin_RSGravitonToJJ_M_1200_gg_hlt = (TH1D*)RSGravitonToJJ_M_1200_gg_hlt.Get( "M_FatPFJet1FatPFJet2_VarBin" ); 
  TH1D *h_mjj_FatPF_VarBin_RSGravitonToJJ_M_2000_gg_hlt = (TH1D*)RSGravitonToJJ_M_2000_gg_hlt.Get( "M_FatPFJet1FatPFJet2_VarBin" ); 

  TH1D *h_mjj_FatPF_VarBin_RSGravitonToJJ_M_700_qq_hlt = (TH1D*)RSGravitonToJJ_M_700_qq_hlt.Get( "M_FatPFJet1FatPFJet2_VarBin" ); 
  TH1D *h_mjj_FatPF_VarBin_RSGravitonToJJ_M_1200_qq_hlt = (TH1D*)RSGravitonToJJ_M_1200_qq_hlt.Get( "M_FatPFJet1FatPFJet2_VarBin" ); 
  TH1D *h_mjj_FatPF_VarBin_RSGravitonToJJ_M_2000_qq_hlt = (TH1D*)RSGravitonToJJ_M_2000_qq_hlt.Get( "M_FatPFJet1FatPFJet2_VarBin" ); 

  //MC - PFJets
  TH1D *h_mjj_PF_VarBin_QstarToJJ_M_500_qg_hlt = (TH1D*)QstarToJJ_M_500_qg_hlt.Get( "M_PFJet1PFJet2_VarBin" ); 
  TH1D *h_mjj_PF_VarBin_QstarToJJ_M_700_qg_hlt = (TH1D*)QstarToJJ_M_700_qg_hlt.Get( "M_PFJet1PFJet2_VarBin" ); 
  TH1D *h_mjj_PF_VarBin_QstarToJJ_M_1200_qg_hlt = (TH1D*)QstarToJJ_M_1200_qg_hlt.Get( "M_PFJet1PFJet2_VarBin" ); 
  TH1D *h_mjj_PF_VarBin_QstarToJJ_M_2000_qg_hlt = (TH1D*)QstarToJJ_M_2000_qg_hlt.Get( "M_PFJet1PFJet2_VarBin" ); 

  TH1D *h_mjj_PF_VarBin_RSGravitonToJJ_M_700_gg_hlt = (TH1D*)RSGravitonToJJ_M_700_gg_hlt.Get( "M_PFJet1PFJet2_VarBin" ); 
  TH1D *h_mjj_PF_VarBin_RSGravitonToJJ_M_1200_gg_hlt = (TH1D*)RSGravitonToJJ_M_1200_gg_hlt.Get( "M_PFJet1PFJet2_VarBin" ); 
  TH1D *h_mjj_PF_VarBin_RSGravitonToJJ_M_2000_gg_hlt = (TH1D*)RSGravitonToJJ_M_2000_gg_hlt.Get( "M_PFJet1PFJet2_VarBin" ); 

  TH1D *h_mjj_PF_VarBin_RSGravitonToJJ_M_700_qq_hlt = (TH1D*)RSGravitonToJJ_M_700_qq_hlt.Get( "M_PFJet1PFJet2_VarBin" ); 
  TH1D *h_mjj_PF_VarBin_RSGravitonToJJ_M_1200_qq_hlt = (TH1D*)RSGravitonToJJ_M_1200_qq_hlt.Get( "M_PFJet1PFJet2_VarBin" ); 
  TH1D *h_mjj_PF_VarBin_RSGravitonToJJ_M_2000_qq_hlt = (TH1D*)RSGravitonToJJ_M_2000_qq_hlt.Get( "M_PFJet1PFJet2_VarBin" ); 


  //data
  TCanvas c0;
  c0->SetLogy();
  h_mjj_FatPF_VarBin_data_hlt->GetXaxis()->SetRangeUser(0,XMAX);
  h_mjj_FatPF_VarBin_data_hlt->GetYaxis()->SetRangeUser(0.01,h_mjj_FatPF_VarBin_data_hlt->GetMaximum()*10);
  h_mjj_FatPF_VarBin_data_hlt->GetXaxis()->SetTitle("FatPFJetMass@HLT [GeV]");
  h_mjj_FatPF_VarBin_data_hlt->SetMarkerStyle(20);
  h_mjj_FatPF_VarBin_data_hlt->SetMarkerSize(0.5);
  h_mjj_FatPF_VarBin_data_hlt->Draw("HISTE");

  TLine line_FatPF_VarBin_data_hlt(min_mjj_FatPF,0,min_mjj_FatPF,h_mjj_FatPF_VarBin_data_hlt->GetMaximum());
  line_FatPF_VarBin_data_hlt->SetLineColor(2);
  line_FatPF_VarBin_data_hlt->Draw();

  c0.SaveAs("data_hlt.eps");

  //qg
  TCanvas c1;
  h_mjj_FatPF_VarBin_QstarToJJ_M_500_qg_hlt->GetXaxis()->SetRangeUser(0,XMAX);
  h_mjj_FatPF_VarBin_QstarToJJ_M_500_qg_hlt->GetXaxis()->SetTitle("FatPFJetMass@HLT [GeV]");
  h_mjj_FatPF_VarBin_QstarToJJ_M_500_qg_hlt->SetLineColor(1);
  h_mjj_FatPF_VarBin_QstarToJJ_M_700_qg_hlt->SetLineColor(2);
  h_mjj_FatPF_VarBin_QstarToJJ_M_1200_qg_hlt->SetLineColor(3);
  h_mjj_FatPF_VarBin_QstarToJJ_M_2000_qg_hlt->SetLineColor(4);

  h_mjj_FatPF_VarBin_QstarToJJ_M_500_qg_hlt->SetStats(0);
  h_mjj_FatPF_VarBin_QstarToJJ_M_500_qg_hlt->DrawNormalized("HISTE");
  h_mjj_FatPF_VarBin_QstarToJJ_M_700_qg_hlt->DrawNormalized("HISTEsame");
  h_mjj_FatPF_VarBin_QstarToJJ_M_1200_qg_hlt->DrawNormalized("HISTEsame");
  h_mjj_FatPF_VarBin_QstarToJJ_M_2000_qg_hlt->DrawNormalized("HISTEsame");

  TLegend legend1 (0.55,0.55,0.85,0.85);
  legend1.SetFillColor(kWhite);
  legend1.SetMargin(0.2);
  legend1.AddEntry(h_mjj_FatPF_VarBin_QstarToJJ_M_500_qg_hlt,"QstarToJJ_M_500_qg_hlt","l");
  legend1.AddEntry(h_mjj_FatPF_VarBin_QstarToJJ_M_700_qg_hlt,"QstarToJJ_M_700_qg_hlt","l");
  legend1.AddEntry(h_mjj_FatPF_VarBin_QstarToJJ_M_1200_qg_hlt,"QstarToJJ_M_1200_qg_hlt","l");
  legend1.AddEntry(h_mjj_FatPF_VarBin_QstarToJJ_M_2000_qg_hlt,"QstarToJJ_M_2000_qg_hlt","l");
  legend1.Draw();
  c1.Update();
  gPad->RedrawAxis();
  gPad->Modified();

  c1.SaveAs("Qstar_qg.eps");

  //gg
  TCanvas c2;
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_700_gg_hlt->GetXaxis()->SetRangeUser(0,XMAX);
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_700_gg_hlt->GetXaxis()->SetTitle("FatPFJetMass@HLT [GeV]");
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_700_gg_hlt->SetLineColor(2);
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_1200_gg_hlt->SetLineColor(3);
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_2000_gg_hlt->SetLineColor(4);

  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_700_gg_hlt->SetStats(0);
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_700_gg_hlt->DrawNormalized("HISTE");
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_1200_gg_hlt->DrawNormalized("HISTEsame");
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_2000_gg_hlt->DrawNormalized("HISTEsame");

  TLegend legend2 (0.55,0.55,0.85,0.85);
  legend2.SetFillColor(kWhite);
  legend2.SetMargin(0.2);
  legend2.AddEntry(h_mjj_FatPF_VarBin_RSGravitonToJJ_M_700_gg_hlt,"RSGravitonToJJ_M_700_gg_hlt","l");
  legend2.AddEntry(h_mjj_FatPF_VarBin_RSGravitonToJJ_M_1200_gg_hlt,"RSGravitonToJJ_M_1200_gg_hlt","l");
  legend2.AddEntry(h_mjj_FatPF_VarBin_RSGravitonToJJ_M_2000_gg_hlt,"RSGravitonToJJ_M_2000_gg_hlt","l");
  legend2.Draw();
  c2.Update();
  gPad->RedrawAxis();
  gPad->Modified();

  c2.SaveAs("RSGraviton_gg.eps");

  //qq
  TCanvas c3;
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_700_qq_hlt->GetXaxis()->SetRangeUser(0,XMAX);
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_700_qq_hlt->GetXaxis()->SetTitle("FatPFJetMass@HLT [GeV]");
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_700_qq_hlt->SetLineColor(2);
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_1200_qq_hlt->SetLineColor(3);
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_2000_qq_hlt->SetLineColor(4);

  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_700_qq_hlt->SetStats(0);
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_700_qq_hlt->DrawNormalized("HISTE");
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_1200_qq_hlt->DrawNormalized("HISTEsame");
  h_mjj_FatPF_VarBin_RSGravitonToJJ_M_2000_qq_hlt->DrawNormalized("HISTEsame");

  TLegend legend3 (0.55,0.55,0.85,0.85);
  legend3.SetFillColor(kWhite);
  legend3.SetMargin(0.2);
  legend3.AddEntry(h_mjj_FatPF_VarBin_RSGravitonToJJ_M_700_qq_hlt,"RSGravitonToJJ_M_700_qq_hlt","l");
  legend3.AddEntry(h_mjj_FatPF_VarBin_RSGravitonToJJ_M_1200_qq_hlt,"RSGravitonToJJ_M_1200_qq_hlt","l");
  legend3.AddEntry(h_mjj_FatPF_VarBin_RSGravitonToJJ_M_2000_qq_hlt,"RSGravitonToJJ_M_2000_qq_hlt","l");
  legend3.Draw();
  c3.Update();
  gPad->RedrawAxis();
  gPad->Modified();

  c3.SaveAs("RSGraviton_qq.eps");

}
