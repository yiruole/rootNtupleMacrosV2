{
  gROOT->Reset();
  gStyle->SetOptFit(1111);

  //TFile file("HLT_RECO_HT_1kevents_TightJetID.root");
  //TFile file("HLT_vs_RECO/HTDATA_HLTvsRECO.root"); //DATA (22k events)
  //TFile file("HLT_vs_RECO/QstarToJJ500MC_HLTvsRECO.root"); //MC (86k events)

  //2011 Selection (only two leading jets, deltaphi cut)
  //TFile file("HLT_vs_RECO/HTDATA_HLTvsRECO_2011Sel.root"); //DATA (22k events)
  //TFile file("HLT_vs_RECO/HTDATA_HLTvsRECO_2011Sel_v2.root"); //DATA (22k events, added JetID in "first mismatch" plot)
  //TFile file("/afs/cern.ch/user/s/santanas/scratch0/DiJets/rootNtupleAnalyzerV2/analysisClass_dijetPhysicsDST_2011Sel_HLTvsRECOcomparison.root"); //DATA (22k events, added JetID in "first mismatch" plot, jetmass bias)
  //TFile file("HLT_vs_RECO/QstarToJJ500MC_HLTvsRECO_2011Sel.root"); //MC (86k events)
  TFile file("HLT_vs_RECO/QstarToJJ500MC_HLTvsRECO_2011Sel_v2.root"); //MC (86k events, added JetID in "first mismatch" plot)
  double Case=1;// DATA(0) , MC(1)

  //----------------------------------------------------------------------------------

  TH1D *hClone_HLT_M_PFJet1PFJet2 = (TH1D*)file.Get( "h_HLT_M_PFJet1PFJet2" )->Clone(); 
  TH1D *hClone_RECO_M_PFJet1PFJet2 = (TH1D*)file.Get( "h_RECO_M_PFJet1PFJet2" )->Clone(); 

  TH2D *h2Clone_FinalSelectionAgreement_recoY_vs_hltX = (TH2D*)file.Get( "h2_FinalSelectionAgreement_recoY_vs_hltX" )->Clone(); 

  TH1D *hClone_HLT_M_PFJet1PFJet2_MATCH = (TH1D*)file.Get( "h_HLT_M_PFJet1PFJet2_MATCH" )->Clone(); 
  TH1D *hClone_RECO_M_PFJet1PFJet2_MATCH = (TH1D*)file.Get( "h_RECO_M_PFJet1PFJet2_MATCH" )->Clone(); 

  TH1D *hClone_MassBias_MATCH = (TH1D*)file.Get( "h_MassBias_MATCH" )->Clone(); 

  TProfile *pClone_MassBias_MATCH_vs_Mreco = (TProfile*)file.Get( "p_MassBias_MATCH_vs_Mreco" )->Clone();

  TH1D *hClone_HLT_M_PFJet1PFJet2_MATCH_MassBiasMore0p2 = (TH1D*)file.Get( "h_HLT_M_PFJet1PFJet2_MATCH_MassBiasMore0p2" )->Clone(); 
  TH1D *hClone_RECO_M_PFJet1PFJet2_MATCH_MassBiasMore0p2 = (TH1D*)file.Get( "h_RECO_M_PFJet1PFJet2_MATCH_MassBiasMore0p2" )->Clone(); 

  TH1D *hClone_HLT_M_PFJet1PFJet2_HLT1_RECO0 = (TH1D*)file.Get( "h_HLT_M_PFJet1PFJet2_HLT1_RECO0" )->Clone(); 
  TH1D *hClone_RECO_M_PFJet1PFJet2_HLT1_RECO0 = (TH1D*)file.Get( "h_RECO_M_PFJet1PFJet2_HLT1_RECO0" )->Clone(); 

  TH1D *hClone_FirstMismatch_HLT1_RECO0 = (TH1D*)file.Get( "h_FirstMismatch_HLT1_RECO0" )->Clone(); 
  TH1D *hClone_FirstMismatch_HLT0_RECO1 = (TH1D*)file.Get( "h_FirstMismatch_HLT0_RECO1" )->Clone(); 

  TH1D *hClone_HLT_M_PFJet1PFJet2_HLT0_RECO1 = (TH1D*)file.Get( "h_HLT_M_PFJet1PFJet2_HLT0_RECO1" )->Clone(); 
  TH1D *hClone_RECO_M_PFJet1PFJet2_HLT0_RECO1 = (TH1D*)file.Get( "h_RECO_M_PFJet1PFJet2_HLT0_RECO1" )->Clone(); 
  
  TH1D *hClone_HLTPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4 = (TH1D*)file.Get( "h_HLTPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4" )->Clone();
  TH1D *hClone_HLTPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4 = (TH1D*)file.Get( "h_HLTPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4" )->Clone();
  TH1D *hClone_HLTPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4 = (TH1D*)file.Get( "h_HLTPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4" )->Clone();
  TH1D *hClone_HLTPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4 = (TH1D*)file.Get( "h_HLTPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4" )->Clone();
  TH1D *hClone_HLTPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4 = (TH1D*)file.Get( "h_HLTPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4" )->Clone();
  TH1D *hClone_HLTPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4 = (TH1D*)file.Get( "h_HLTPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4" )->Clone();

  TH1D *hClone_RECOPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4 = (TH1D*)file.Get( "h_RECOPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4" )->Clone();
  TH1D *hClone_RECOPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4 = (TH1D*)file.Get( "h_RECOPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4" )->Clone();
  TH1D *hClone_RECOPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4 = (TH1D*)file.Get( "h_RECOPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4" )->Clone();
  TH1D *hClone_RECOPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4 = (TH1D*)file.Get( "h_RECOPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4" )->Clone();
  TH1D *hClone_RECOPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4 = (TH1D*)file.Get( "h_RECOPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4" )->Clone();
  TH1D *hClone_RECOPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4 = (TH1D*)file.Get( "h_RECOPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4" )->Clone();

  TH1D *hClone_PtBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4 = (TH1D*)file.Get( "h_PtBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4" )->Clone();
  TH1D *hClone_PtBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4 = (TH1D*)file.Get( "h_PtBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4" )->Clone();
  TH1D *hClone_PtBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4 = (TH1D*)file.Get( "h_PtBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4" )->Clone();
  TH1D *hClone_PtBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4 = (TH1D*)file.Get( "h_PtBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4" )->Clone();
  TH1D *hClone_PtBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3 = (TH1D*)file.Get( "h_PtBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3" )->Clone();
  TH1D *hClone_PtBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3 = (TH1D*)file.Get( "h_PtBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3" )->Clone();
  TH1D *hClone_PtBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3 = (TH1D*)file.Get( "h_PtBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3" )->Clone();
  TH1D *hClone_PtBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3 = (TH1D*)file.Get( "h_PtBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3" )->Clone();

  /*
  TH1D *hClone_JetMassBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4 = (TH1D*)file.Get( "h_JetMassBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4" )->Clone();
  TH1D *hClone_JetMassBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4 = (TH1D*)file.Get( "h_JetMassBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4" )->Clone();
  TH1D *hClone_JetMassBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4 = (TH1D*)file.Get( "h_JetMassBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4" )->Clone();
  TH1D *hClone_JetMassBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4 = (TH1D*)file.Get( "h_JetMassBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4" )->Clone();
  TH1D *hClone_JetMassBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3 = (TH1D*)file.Get( "h_JetMassBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3" )->Clone();
  TH1D *hClone_JetMassBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3 = (TH1D*)file.Get( "h_JetMassBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3" )->Clone();
  TH1D *hClone_JetMassBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3 = (TH1D*)file.Get( "h_JetMassBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3" )->Clone();
  TH1D *hClone_JetMassBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3 = (TH1D*)file.Get( "h_JetMassBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3" )->Clone();
  */


  //--------------------------------------------------------------------

  TCanvas c1;
  hClone_HLT_M_PFJet1PFJet2->GetXaxis()->SetTitle("M(Jet1,Jet2) [GeV]");
  hClone_HLT_M_PFJet1PFJet2->SetLineColor(1);
  hClone_RECO_M_PFJet1PFJet2->SetLineColor(2);
  hClone_HLT_M_PFJet1PFJet2->Draw();
  hClone_RECO_M_PFJet1PFJet2->Draw("sames");
  hClone_HLT_M_PFJet1PFJet2->GetYaxis()->SetRangeUser(0,TMath::Max( hClone_RECO_M_PFJet1PFJet2->GetMaximum()*1.2,hClone_HLT_M_PFJet1PFJet2->GetMaximum()*1.2));

  TLegend legend1 (0.35,0.65,0.75,0.85);
  legend1.SetFillColor(kWhite);
  legend1.SetMargin(0.2);
  legend1.AddEntry(hClone_HLT_M_PFJet1PFJet2,"HLT","pl");
  legend1.AddEntry(hClone_RECO_M_PFJet1PFJet2,"RECO","pl");
  legend1.Draw();

  if(Case==0)
    c1.SaveAs("DATA_MJJ_HLT_vs_RECO.eps");
  else if(Case==1)
    c1.SaveAs("MC_MJJ_HLT_vs_RECO.eps");

  TCanvas c2;
  h2Clone_FinalSelectionAgreement_recoY_vs_hltX->GetXaxis()->SetTitle("Pass HLT selection");
  h2Clone_FinalSelectionAgreement_recoY_vs_hltX->GetYaxis()->SetTitle("Pass RECO selection");
  h2Clone_FinalSelectionAgreement_recoY_vs_hltX->SetMarkerSize(5);
  h2Clone_FinalSelectionAgreement_recoY_vs_hltX->Draw("textcolz");

  if(Case==0)
    c2.SaveAs("DATA_FinalSelectionAgreement.eps");
  else if(Case==1)
    c2.SaveAs("MC_FinalSelectionAgreement.eps");

  TCanvas c3;
  hClone_HLT_M_PFJet1PFJet2_MATCH->GetXaxis()->SetTitle("M(Jet1,Jet2) [GeV]");
  hClone_HLT_M_PFJet1PFJet2_MATCH->SetLineColor(1);
  hClone_RECO_M_PFJet1PFJet2_MATCH->SetLineColor(2);
  hClone_HLT_M_PFJet1PFJet2_MATCH->Draw();
  hClone_RECO_M_PFJet1PFJet2_MATCH->Draw("sames");
  hClone_HLT_M_PFJet1PFJet2_MATCH->GetYaxis()->SetRangeUser(0,TMath::Max( hClone_RECO_M_PFJet1PFJet2_MATCH->GetMaximum()*1.2,hClone_HLT_M_PFJet1PFJet2_MATCH->GetMaximum()*1.2));

  TLegend legend3 (0.35,0.65,0.75,0.85);
  legend3.SetFillColor(kWhite);
  legend3.SetMargin(0.2);
  legend3.AddEntry(hClone_HLT_M_PFJet1PFJet2_MATCH,"HLT","pl");
  legend3.AddEntry(hClone_RECO_M_PFJet1PFJet2_MATCH,"RECO","pl");
  legend3.Draw();

  if(Case==0)
    c3.SaveAs("DATA_MJJ_HLT_vs_RECO_MATCH.eps");
  else if(Case==1)
    c3.SaveAs("MC_MJJ_HLT_vs_RECO_MATCH.eps");

  TCanvas c4;
  TF1 *f_gaus = new TF1("f_gaus","gaus",-0.1,0.1);
  f_gaus->SetLineColor(2);
  f_gaus->SetLineWidth(2);
  hClone_MassBias_MATCH->GetXaxis()->SetTitle("(M_{HLT} - M_{RECO}) / M_{RECO}");
  hClone_MassBias_MATCH->SetMarkerStyle(20);
  hClone_MassBias_MATCH->Fit(f_gaus,"R");
  hClone_MassBias_MATCH->Draw();

  if(Case==0)
    c4.SaveAs("DATA_MJJBias_MATCH.eps");
  else if(Case==1)
    c4.SaveAs("MC_MJJBias_MATCH.eps");

  TCanvas c5;
  c5->SetGridx();
  c5->SetGridy();
  TF1 *f_pol0 = new TF1("f_pol0","pol0",400,1300);
  f_pol0->SetLineColor(2);
  f_pol0->SetLineWidth(2);
  pClone_MassBias_MATCH_vs_Mreco->GetYaxis()->SetRangeUser(-0.1,0.1);
  pClone_MassBias_MATCH_vs_Mreco->GetYaxis()->SetTitle("(M_{HLT} - M_{RECO}) / M_{RECO}");
  pClone_MassBias_MATCH_vs_Mreco->GetYaxis()->SetTitleOffset(1.2);
  pClone_MassBias_MATCH_vs_Mreco->GetXaxis()->SetTitle("M_{RECO}");
  pClone_MassBias_MATCH_vs_Mreco->SetMarkerStyle(20);
  pClone_MassBias_MATCH_vs_Mreco->Fit(f_pol0,"R");
  pClone_MassBias_MATCH_vs_Mreco->Draw();

  if(Case==0)
    c5.SaveAs("DATA_MJJBiasVsMJJ_MATCH.eps");
  else if(Case==1)
    c5.SaveAs("MC_MJJBiasVsMJJ_MATCH.eps");

  TCanvas c6;
  hClone_HLT_M_PFJet1PFJet2_MATCH_MassBiasMore0p2->GetXaxis()->SetTitle("M(Jet1,Jet2) [GeV]");
  hClone_HLT_M_PFJet1PFJet2_MATCH_MassBiasMore0p2->SetLineColor(1);
  hClone_RECO_M_PFJet1PFJet2_MATCH_MassBiasMore0p2->SetLineColor(2);
  hClone_HLT_M_PFJet1PFJet2_MATCH_MassBiasMore0p2->Draw();
  hClone_RECO_M_PFJet1PFJet2_MATCH_MassBiasMore0p2->Draw("sames");
  hClone_HLT_M_PFJet1PFJet2_MATCH_MassBiasMore0p2->GetYaxis()->SetRangeUser(0,TMath::Max( hClone_RECO_M_PFJet1PFJet2_MATCH_MassBiasMore0p2->GetMaximum()*1.2,hClone_HLT_M_PFJet1PFJet2_MATCH_MassBiasMore0p2->GetMaximum()*1.2));

  TLegend legend6 (0.35,0.65,0.75,0.85);
  legend6.SetFillColor(kWhite);
  legend6.SetMargin(0.2);
  legend6.AddEntry(hClone_HLT_M_PFJet1PFJet2_MATCH_MassBiasMore0p2,"HLT","pl");
  legend6.AddEntry(hClone_RECO_M_PFJet1PFJet2_MATCH_MassBiasMore0p2,"RECO","pl");
  legend6.Draw();

  if(Case==0)
    c6.SaveAs("DATA_MJJ_HLT_vs_RECO_MATCH_MassBiasMore0p2.eps");
  else if(Case==1)
    c6.SaveAs("MC_MJJ_HLT_vs_RECO_MATCH_MassBiasMore0p2.eps");

  TCanvas c7;
  hClone_HLT_M_PFJet1PFJet2_HLT1_RECO0->GetXaxis()->SetTitle("M(Jet1,Jet2) [GeV]");
  hClone_HLT_M_PFJet1PFJet2_HLT1_RECO0->SetLineColor(1);
  hClone_RECO_M_PFJet1PFJet2_HLT1_RECO0->SetLineColor(2);
  hClone_HLT_M_PFJet1PFJet2_HLT1_RECO0->Draw();
  //hClone_RECO_M_PFJet1PFJet2_HLT1_RECO0->Draw("sames");
  //hClone_HLT_M_PFJet1PFJet2_HLT1_RECO0->GetYaxis()->SetRangeUser(0,TMath::Max( hClone_RECO_M_PFJet1PFJet2_HLT1_RECO0->GetMaximum()*1.2,hClone_HLT_M_PFJet1PFJet2_HLT1_RECO0->GetMaximum()*1.2));

  TLegend legend7 (0.35,0.65,0.75,0.85);
  legend7.SetFillColor(kWhite);
  legend7.SetMargin(0.2);
  legend7.AddEntry(hClone_HLT_M_PFJet1PFJet2_HLT1_RECO0,"HLT","pl");
  //legend7.AddEntry(hClone_RECO_M_PFJet1PFJet2_HLT1_RECO0,"RECO","pl");
  legend7.Draw();

  if(Case==0)
    c7.SaveAs("DATA_MJJ_HLT_HLT1_RECO0.eps");
  else if(Case==1)
    c7.SaveAs("MC_MJJ_HLT_HLT1_RECO0.eps");

  TCanvas c8;
  c8->SetGridx();
  c8->SetGridy();  
  hClone_FirstMismatch_HLT1_RECO0->SetMarkerStyle(20);
  hClone_FirstMismatch_HLT1_RECO0->SetFillColor(2);
  hClone_FirstMismatch_HLT1_RECO0->SetFillStyle(3002);
  hClone_FirstMismatch_HLT1_RECO0->Draw("HISTE");

  if(Case==0)
    c8.SaveAs("DATA_FirstMismatch_HLT1_RECO0.eps");
  else if(Case==1)
    c8.SaveAs("MC_FirstMismatch_HLT1_RECO0.eps");

  TCanvas c9;
  c9->SetGridx();
  c9->SetGridy();  
  hClone_FirstMismatch_HLT0_RECO1->SetMarkerStyle(20);
  hClone_FirstMismatch_HLT0_RECO1->SetFillColor(2);
  hClone_FirstMismatch_HLT0_RECO1->SetFillStyle(3002);
  hClone_FirstMismatch_HLT0_RECO1->Draw("HISTE");

  if(Case==0)
    c9.SaveAs("DATA_FirstMismatch_HLT0_RECO1.eps");
  else if(Case==1)
    c9.SaveAs("MC_FirstMismatch_HLT0_RECO1.eps");

  TCanvas c10;
  hClone_RECO_M_PFJet1PFJet2_HLT0_RECO1->GetXaxis()->SetTitle("M(Jet1,Jet2) [GeV]");
  hClone_RECO_M_PFJet1PFJet2_HLT0_RECO1->SetLineColor(2);
  hClone_HLT_M_PFJet1PFJet2_HLT0_RECO1->SetLineColor(1);
  hClone_RECO_M_PFJet1PFJet2_HLT0_RECO1->Draw();
  //hClone_HLT_M_PFJet1PFJet2_HLT0_RECO1->Draw("sames");

  TLegend legend10 (0.35,0.65,0.75,0.85);
  legend10.SetFillColor(kWhite);
  legend10.SetMargin(0.2);
  legend10.AddEntry(hClone_RECO_M_PFJet1PFJet2_HLT0_RECO1,"RECO","pl");
  legend10.Draw();

  if(Case==0)
    c10.SaveAs("DATA_MJJ_RECO_HLT0_RECO1.eps");
  else if(Case==1)
    c10.SaveAs("MC_MJJ_RECO_HLT0_RECO1.eps");

  TCanvas c11;
  c11.Divide(3,2);

  c11.cd(1);
  hClone_HLTPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->SetLineColor(1);  
  hClone_RECOPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->SetLineColor(2);  
  hClone_HLTPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->Draw();
  hClone_RECOPFJetNeutralHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->Draw("sames");

  c11.cd(2);
  hClone_HLTPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->SetLineColor(1);  
  hClone_RECOPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->SetLineColor(2);  
  hClone_HLTPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->Draw();
  hClone_RECOPFJetNeutralEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->Draw("sames");

  c11.cd(3);
  hClone_HLTPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4->GetXaxis()->SetRangeUser(0,100);
  hClone_HLTPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4->SetLineColor(1);  
  hClone_RECOPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4->SetLineColor(2);  
  hClone_HLTPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4->Draw();
  hClone_RECOPFJetNConstituents_smallestDR_Pt30GeV_EtaL2p4->Draw("sames");

  c11.cd(4);
  hClone_HLTPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->GetYaxis()->SetRangeUser(0,hClone_RECOPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->GetMaximum()*1.2);
  hClone_HLTPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->SetLineColor(1);  
  hClone_RECOPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->SetLineColor(2);  
  hClone_HLTPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->Draw();
  hClone_RECOPFJetChargedHadronEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->Draw("sames");

  c11.cd(5);
  hClone_HLTPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4->GetXaxis()->SetRangeUser(0,100);
  hClone_HLTPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4->SetLineColor(1);  
  hClone_RECOPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4->SetLineColor(2);  
  hClone_HLTPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4->Draw();
  hClone_RECOPFJetChargedMultiplicity_smallestDR_Pt30GeV_EtaL2p4->Draw("sames");

  c11.cd(6);
  c11.cd(6)->SetLogy();
  c11.cd(6)->SetLogx();
  hClone_HLTPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->GetYaxis()->SetRangeUser(0.1,hClone_HLTPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->GetMaximum()*10);
  hClone_HLTPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->SetLineColor(1);  
  hClone_RECOPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->SetLineColor(2);  
  hClone_HLTPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->Draw();
  hClone_RECOPFJetChargedEmEnergyFraction_smallestDR_Pt30GeV_EtaL2p4->Draw("sames");

  if(Case==0)
    c11.SaveAs("DATA_JetID_HLT_vs_RECO.eps");
  else if(Case==1)
    c11.SaveAs("MC_JetID_HLT_vs_RECO.eps");

  TCanvas c12;
  c12.Divide(2,1);

  c12.cd(1);
  hClone_PtBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4->SetLineColor(1);
  hClone_PtBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4->SetMarkerColor(1);
  hClone_PtBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4->SetLineColor(2);
  hClone_PtBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4->SetMarkerColor(2);
  hClone_PtBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4->SetLineColor(3);
  hClone_PtBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4->SetMarkerColor(3);
  hClone_PtBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4->SetLineColor(4);
  hClone_PtBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4->SetMarkerColor(4);

  hClone_PtBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4->DrawNormalized();
  hClone_PtBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4->DrawNormalized("sames");
  hClone_PtBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4->DrawNormalized("sames");
  hClone_PtBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4->DrawNormalized("sames");

  TLegend legend12_a (0.35,0.65,0.75,0.85);
  legend12_a.SetFillColor(kWhite);
  legend12_a.SetMargin(0.2);
  legend12_a.AddEntry(hClone_PtBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4,"PT>10 GeV","pl");
  legend12_a.AddEntry(hClone_PtBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4,"PT>30 GeV","pl");
  legend12_a.AddEntry(hClone_PtBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4,"PT>50 GeV","pl");
  legend12_a.AddEntry(hClone_PtBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4,"PT>100 GeV","pl");
  legend12_a.Draw();

  c12.cd(2);
  hClone_PtBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3->SetLineColor(1);
  hClone_PtBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3->SetMarkerColor(1);
  hClone_PtBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3->SetLineColor(2);
  hClone_PtBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3->SetMarkerColor(2);
  hClone_PtBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3->SetLineColor(3);
  hClone_PtBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3->SetMarkerColor(3);
  hClone_PtBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3->SetLineColor(4);
  hClone_PtBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3->SetMarkerColor(4);

  hClone_PtBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3->DrawNormalized();
  hClone_PtBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3->DrawNormalized("sames");
  hClone_PtBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3->DrawNormalized("sames");
  hClone_PtBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3->DrawNormalized("sames");

  TLegend legend12_b (0.35,0.65,0.75,0.85);
  legend12_b.SetFillColor(kWhite);
  legend12_b.SetMargin(0.2);
  legend12_b.AddEntry(hClone_PtBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3,"PT>10 GeV","pl");
  legend12_b.AddEntry(hClone_PtBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3,"PT>30 GeV","pl");
  legend12_b.AddEntry(hClone_PtBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3,"PT>50 GeV","pl");
  legend12_b.AddEntry(hClone_PtBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3,"PT>100 GeV","pl");
  legend12_b.Draw();

  if(Case==0)
    c12.SaveAs("DATA_JetPtBias_HLT_vs_RECO.eps");
  else if(Case==1)
    c12.SaveAs("MC_JetPtBias_HLT_vs_RECO.eps");

  /*
  TCanvas c13;
  c13.Divide(2,1);

  c13.cd(1);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4->SetLineColor(1);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4->SetMarkerColor(1);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4->SetLineColor(2);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4->SetMarkerColor(2);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4->SetLineColor(3);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4->SetMarkerColor(3);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4->SetLineColor(4);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4->SetMarkerColor(4);

  hClone_JetMassBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4->DrawNormalized();
  hClone_JetMassBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4->DrawNormalized("sames");
  hClone_JetMassBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4->DrawNormalized("sames");
  hClone_JetMassBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4->DrawNormalized("sames");

  TLegend legend12_a (0.35,0.65,0.75,0.85);
  legend12_a.SetFillColor(kWhite);
  legend12_a.SetMargin(0.2);
  legend12_a.AddEntry(hClone_JetMassBias_reco_hlt_smallestDR_Pt10GeV_EtaL2p4,"PT>10 GeV","pl");
  legend12_a.AddEntry(hClone_JetMassBias_reco_hlt_smallestDR_Pt30GeV_EtaL2p4,"PT>30 GeV","pl");
  legend12_a.AddEntry(hClone_JetMassBias_reco_hlt_smallestDR_Pt50GeV_EtaL2p4,"PT>50 GeV","pl");
  legend12_a.AddEntry(hClone_JetMassBias_reco_hlt_smallestDR_Pt100GeV_EtaL2p4,"PT>100 GeV","pl");
  legend12_a.Draw();

  c13.cd(2);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3->SetLineColor(1);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3->SetMarkerColor(1);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3->SetLineColor(2);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3->SetMarkerColor(2);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3->SetLineColor(3);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3->SetMarkerColor(3);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3->SetLineColor(4);
  hClone_JetMassBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3->SetMarkerColor(4);

  hClone_JetMassBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3->DrawNormalized();
  hClone_JetMassBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3->DrawNormalized("sames");
  hClone_JetMassBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3->DrawNormalized("sames");
  hClone_JetMassBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3->DrawNormalized("sames");

  TLegend legend12_b (0.35,0.65,0.75,0.85);
  legend12_b.SetFillColor(kWhite);
  legend12_b.SetMargin(0.2);
  legend12_b.AddEntry(hClone_JetMassBias_reco_hlt_smallestDR_Pt10GeV_Eta2p4_3,"PT>10 GeV","pl");
  legend12_b.AddEntry(hClone_JetMassBias_reco_hlt_smallestDR_Pt30GeV_Eta2p4_3,"PT>30 GeV","pl");
  legend12_b.AddEntry(hClone_JetMassBias_reco_hlt_smallestDR_Pt50GeV_Eta2p4_3,"PT>50 GeV","pl");
  legend12_b.AddEntry(hClone_JetMassBias_reco_hlt_smallestDR_Pt100GeV_Eta2p4_3,"PT>100 GeV","pl");
  legend12_b.Draw();

  if(Case==0)
    c13.SaveAs("DATA_JetMassBias_HLT_vs_RECO.eps");
  else if(Case==1)
    c13.SaveAs("MC_JetMassBias_HLT_vs_RECO.eps");
  */

}
