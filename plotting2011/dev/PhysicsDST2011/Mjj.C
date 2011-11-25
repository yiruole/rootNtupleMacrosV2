{
  gROOT->Reset();
  gStyle->SetOptStat(1111111);
  gStyle->SetOptFit(1111);

  //NoJEC
  //TFile datafile_hlt("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/117pb-1__15_11_2011/PhysicsDST_allEvents.root");
  //L23Res JEC
  //TFile datafile_hlt("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/117pb-1_JECL23Res__19_11_2011/PhysicsDST_allEvents.root");
  //L123Res JEC
  TFile datafile_hlt("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/117pb-1_JECL123Res__24_11_2011/PhysicsDST_allEvents.root");

  TFile datafile_reco("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/DiJetHighMass_ControlPlots_Chiyoung_3p4fb-1/histograms_data_HT_340fb_Fat_ak5.root");

  //NoJEC
  //   TFile Qstar500ForTest_hlt("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/QstarJJ_ForTest__16_11_2011/QstarJJ500.root");
  //   TFile Qstar700ForTest_hlt("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/QstarJJ_ForTest__16_11_2011/QstarJJ700.root");
  //   TFile Qstar1200ForTest_hlt("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/QstarJJ_ForTest__16_11_2011/QstarJJ1200.root");
  //L23 JEC
  //   TFile Qstar500ForTest_hlt("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/QstarJJ_ForTest_JECL23__18_11_2011/QstarJJ500.root");
  //   TFile Qstar700ForTest_hlt("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/QstarJJ_ForTest_JECL23__18_11_2011/QstarJJ700.root");
  //   TFile Qstar1200ForTest_hlt("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/QstarJJ_ForTest_JECL23__18_11_2011/QstarJJ1200.root");
  //L123 JEC
  TFile Qstar500ForTest_hlt("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/QstarJJ_ForTest_JECL123__25_11_2011/QstarJJ500.root");
  TFile Qstar700ForTest_hlt("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/QstarJJ_ForTest_JECL123__25_11_2011/QstarJJ700.root");
  TFile Qstar1200ForTest_hlt("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/QstarJJ_ForTest_JECL123__25_11_2011/QstarJJ1200.root");

  TFile Qstar500_reco("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/DiJetHighMass_ControlPlots_Chiyoung_3p4fb-1/histograms_mc_QstarToJJ_M-500_TuneD6T_ak5.root");
  TFile Qstar700_reco("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/DiJetHighMass_ControlPlots_Chiyoung_3p4fb-1/histograms_mc_QstarToJJ_M-700_TuneD6T_ak5.root");
  TFile Qstar1200_reco("/afs/cern.ch/user/s/santanas/scratch0/DiJets/data/output_fromAFS/dijets_PhysicsDST/DiJetHighMass_ControlPlots_Chiyoung_3p4fb-1/histograms_mc_QstarToJJ_M-1200_TuneD6T_ak5.root");

  //Mass Ranges
  double min_mjj_FatPF = 526;
  double min_mjj_FatCaloC = 325;

  //DATA - FatPFJets
  TH1D *h_mjj_FatPF_DATA = (TH1D*)datafile_hlt.Get( "M_FatPFJet1FatPFJet2" ); 
  TH1D *h_mjj_PF_DATA = (TH1D*)datafile_hlt.Get( "M_PFJet1PFJet2" ); 
  TH1D *h_mjj_FatPF_VarBin_DATA = (TH1D*)datafile_hlt.Get( "M_FatPFJet1FatPFJet2_VarBin" ); 
  TH1D *h_mjj_PF_VarBin_DATA = (TH1D*)datafile_hlt.Get( "M_PFJet1PFJet2_VarBin" ); 
  TH1D *h_Ratio_mjj_FatPF_VarBin_DATA_HLT_over_RECO = (TH1D*) h_mjj_FatPF_VarBin_DATA->Clone();

  //DATA - FatCaloCJets
  TH1D *h_mjj_FatCaloC_DATA = (TH1D*)datafile_hlt.Get( "M_FatCaloCJet1FatCaloCJet2" ); 
  TH1D *h_mjj_CaloC_DATA = (TH1D*)datafile_hlt.Get( "M_CaloCJet1CaloCJet2" ); 
  TH1D *h_mjj_FatCaloC_VarBin_DATA = (TH1D*)datafile_hlt.Get( "M_FatCaloCJet1FatCaloCJet2_VarBin" ); 
  TH1D *h_mjj_CaloC_VarBin_DATA = (TH1D*)datafile_hlt.Get( "M_CaloCJet1CaloCJet2_VarBin" ); 

  //DATA - FatPFJets - RECO
  TH1D *h_mjj_FatPF_VarBin_DATA_RECO = (TH1D*)datafile_reco.Get( "h_DijetMass_data_fat" ); 

  //MC - QstarToJJ - For Tests 
  TH1D *h_mjj_FatPF_VarBin_Qstar500 = (TH1D*)Qstar500ForTest_hlt.Get( "M_FatPFJet1FatPFJet2_VarBin" ); 
  TH1D *h_mjj_FatPF_VarBin_Qstar700 = (TH1D*)Qstar700ForTest_hlt.Get( "M_FatPFJet1FatPFJet2_VarBin" ); 
  TH1D *h_mjj_FatPF_VarBin_Qstar1200 = (TH1D*)Qstar1200ForTest_hlt.Get( "M_FatPFJet1FatPFJet2_VarBin" ); 

  //MC - QstarToJJ - RECO
  TH1D *h_mjj_FatPF_VarBin_Qstar500_RECO = (TH1D*)Qstar500_reco.Get( "h_DijetMass_data_fat" ); 
  TH1D *h_mjj_FatPF_VarBin_Qstar700_RECO = (TH1D*)Qstar700_reco.Get( "h_DijetMass_data_fat" ); 
  TH1D *h_mjj_FatPF_VarBin_Qstar1200_RECO = (TH1D*)Qstar1200_reco.Get( "h_DijetMass_data_fat" ); 


  //--- DATA - FatPFJets
  TCanvas c1;
  c1.Divide(2,1);

  c1.cd(1)->SetLogy();
  h_mjj_FatPF_DATA->SetMarkerStyle(20);
  h_mjj_FatPF_DATA->SetMarkerSize(0.5);
  h_mjj_FatPF_DATA->GetXaxis()->SetRangeUser(200,4000);
  h_mjj_FatPF_DATA->GetXaxis()->SetTitle("FatPFDiJetMass@HLT [GeV]");
  h_mjj_FatPF_DATA->Draw();

  TLine line_mjj_FatPF_DATA(min_mjj_FatPF,0,min_mjj_FatPF,h_mjj_FatPF_DATA->GetMaximum());
  line_mjj_FatPF_DATA->SetLineColor(2);
  line_mjj_FatPF_DATA->Draw();

  c1.cd(2)->SetLogy();
  h_mjj_FatPF_VarBin_DATA->SetMarkerStyle(20);
  h_mjj_FatPF_VarBin_DATA->SetMarkerSize(0.5);
  h_mjj_FatPF_VarBin_DATA->GetXaxis()->SetRangeUser(220,4000);
  h_mjj_FatPF_VarBin_DATA->GetXaxis()->SetTitle("FatPFDiJetMass@HLT [GeV]");
  h_mjj_FatPF_VarBin_DATA->Draw();

  TLine line_mjj_FatPF_VarBin_DATA(min_mjj_FatPF,0,min_mjj_FatPF,h_mjj_FatPF_VarBin_DATA->GetMaximum());
  line_mjj_FatPF_VarBin_DATA->SetLineColor(2);
  line_mjj_FatPF_VarBin_DATA->Draw();

  //--- DATA - FatCaloCJets
  TCanvas c2;
  c2.Divide(2,1);

  c2.cd(1)->SetLogy();
  h_mjj_FatCaloC_DATA->SetMarkerStyle(20);
  h_mjj_FatCaloC_DATA->SetMarkerSize(0.5);
  h_mjj_FatCaloC_DATA->GetXaxis()->SetRangeUser(200,4000);
  h_mjj_FatCaloC_DATA->GetXaxis()->SetTitle("FatCaloCDiJetMass@HLT [GeV]");
  h_mjj_FatCaloC_DATA->Draw();

  TLine line_mjj_FatCaloC_DATA(min_mjj_FatCaloC,0,min_mjj_FatCaloC,h_mjj_FatCaloC_DATA->GetMaximum());
  line_mjj_FatCaloC_DATA->SetLineColor(2);
  line_mjj_FatCaloC_DATA->Draw();

  c2.cd(2)->SetLogy();
  h_mjj_FatCaloC_VarBin_DATA->SetMarkerStyle(20);
  h_mjj_FatCaloC_VarBin_DATA->SetMarkerSize(0.5);
  h_mjj_FatCaloC_VarBin_DATA->GetXaxis()->SetRangeUser(220,4000);
  h_mjj_FatCaloC_VarBin_DATA->GetXaxis()->SetTitle("FatCaloCDiJetMass@HLT [GeV]");
  h_mjj_FatCaloC_VarBin_DATA->Draw();

  TLine line_mjj_FatCaloC_VarBin_DATA(min_mjj_FatCaloC,0,min_mjj_FatCaloC,h_mjj_FatCaloC_VarBin_DATA->GetMaximum());
  line_mjj_FatCaloC_VarBin_DATA->SetLineColor(2);
  line_mjj_FatCaloC_VarBin_DATA->Draw();

  //--- Comparison DATA - FatPFJets (RECO vs HLT)
  TCanvas c3;
  c3.SetLogy();
  h_mjj_FatPF_VarBin_DATA->GetXaxis()->SetRangeUser( 0 , 4000 );
  h_mjj_FatPF_VarBin_DATA->GetYaxis()->SetRangeUser( 0.01 , h_mjj_FatPF_VarBin_DATA->GetMaximum()*10 );
  h_mjj_FatPF_VarBin_DATA->Draw();
  int NeventsDATA = h_mjj_FatPF_VarBin_DATA->Integral(36,83);
  //double scale_factor = NeventsDATA / h_mjj_FatPF_VarBin_DATA_RECO->Integral();
  double scale_factor = 0.1176 / 3.5; //lumi
  cout << "NeventsDATA (PhysicsDST): " << NeventsDATA << endl;
  cout << "Scale Factor: " << scale_factor << endl;
  h_mjj_FatPF_VarBin_DATA_RECO->SetLineColor(2);
  h_mjj_FatPF_VarBin_DATA_RECO->Scale( scale_factor );
  h_mjj_FatPF_VarBin_DATA_RECO->Draw("samesHISTE");
  
  TLegend legend3 (0.35,0.65,0.75,0.85);
  legend3.SetFillColor(kWhite);
  legend3.SetMargin(0.2);
  legend3.AddEntry(h_mjj_FatPF_VarBin_DATA,"117 pb-1 - HLT (L123Res JEC)","pl");
  legend3.AddEntry(h_mjj_FatPF_VarBin_DATA_RECO,"3.5 fb-1 - RECO (L123Res JEC), scaled to Lumi","pl");
  legend3.Draw();
  c3.Update();
  gPad->RedrawAxis();
  gPad->Modified();

  //--- Ratio HLT / RECO
  TCanvas c4;
  c4.SetGridx();
  c4.SetGridy();
  h_Ratio_mjj_FatPF_VarBin_DATA_HLT_over_RECO->GetXaxis()->SetTitle("FatPFDiJetMass [GeV]");
  h_Ratio_mjj_FatPF_VarBin_DATA_HLT_over_RECO->Divide(h_mjj_FatPF_VarBin_DATA_RECO);
  h_Ratio_mjj_FatPF_VarBin_DATA_HLT_over_RECO->SetMarkerStyle(20);
  TF1* pol1 = new TF1("pol1","pol1",850,2500);
  h_Ratio_mjj_FatPF_VarBin_DATA_HLT_over_RECO->Fit(pol1,"R");
  //h_Ratio_mjj_FatPF_VarBin_DATA_HLT_over_RECO->Draw();
  h_Ratio_mjj_FatPF_VarBin_DATA_HLT_over_RECO->GetYaxis()->SetRangeUser(0,2);
  h_Ratio_mjj_FatPF_VarBin_DATA_HLT_over_RECO->GetXaxis()->SetRangeUser( 0 , 4000);

  TLegend legend4 (0.15,0.7,0.4,0.85);
  legend4.SetFillColor(kWhite);
  legend4.SetMargin(0.2);
  legend4.AddEntry(h_Ratio_mjj_FatPF_VarBin_DATA_HLT_over_RECO,"HLT / RECO","pl");
  legend4.Draw();
  c4.Update();
  gPad->RedrawAxis();
  gPad->Modified();
  

 
  //--- Comparison MC - FatPFJets (RECO vs HLT) (only for 1200 GeV makes sense)
  TCanvas c5;
  //c5.SetLogy();
  h_mjj_FatPF_VarBin_Qstar1200->GetXaxis()->SetTitle("FatPFDiJetMass [GeV]");
  h_mjj_FatPF_VarBin_Qstar1200->GetXaxis()->SetRangeUser( 0 , 1600 );
  h_mjj_FatPF_VarBin_Qstar1200->GetYaxis()->SetRangeUser( 0., h_mjj_FatPF_VarBin_Qstar1200->GetMaximum()*2 );
  h_mjj_FatPF_VarBin_Qstar1200->Draw();
  int NeventsMC = h_mjj_FatPF_VarBin_Qstar1200->Integral(36,83);
  double scale_factor_MC = NeventsMC / h_mjj_FatPF_VarBin_Qstar1200_RECO->Integral();
  cout << "NeventsMC (PhysicsDST): " << NeventsMC << endl;
  cout << "Scale Factor: " << scale_factor_MC << endl;
  h_mjj_FatPF_VarBin_Qstar1200_RECO->SetLineColor(2);
  h_mjj_FatPF_VarBin_Qstar1200_RECO->Scale( scale_factor_MC );
  h_mjj_FatPF_VarBin_Qstar1200_RECO->Draw("samesHISTE");
  
  TLegend legend5 (0.15,0.65,0.6,0.85);
  legend5.SetFillColor(kWhite);
  legend5.SetMargin(0.2);
  legend5.AddEntry(h_mjj_FatPF_VarBin_Qstar1200,"QstarJJ M1200 - HLT (L123 JEC)","pl");
  legend5.AddEntry(h_mjj_FatPF_VarBin_Qstar1200_RECO,"QstarJJ M1200 - RECO (L123 JEC), norm. to HLT","pl");
  legend5.Draw();
  c5.Update();
  gPad->RedrawAxis();
  gPad->Modified();





  
}
