{

  gROOT->Reset();

  TFile f("/data/santanas/Ntuples/2011_logs/RootNtuple-HEAD_Nov3-2011BData-PhysicsDST-179959-180282_dijets_allEvents_14112011/output_cutTable_dijetPhysicsDST_skim/analysisClass_dijetPhysicsDST_skim___PhysicsDST__Run2011B-v1__RAW.root");

  //--------------------------
  TCanvas c1;  
  c1.Divide(3,2);

  c1.cd(1);
  PFJetPrecutsNeutralHadronEnergyFraction->GetXaxis()->SetRangeUser(0,1);
  PFJetPrecutsNeutralHadronEnergyFraction->Draw("HISTE");

  c1.cd(2);
  PFJetPrecutsNeutralEmEnergyFraction->GetXaxis()->SetRangeUser(0,1);
  PFJetPrecutsNeutralEmEnergyFraction->Draw("HISTE");

  c1.cd(3);
  PFJetPrecutsNConstituents->GetXaxis()->SetRangeUser(0,100);
  PFJetPrecutsNConstituents->Draw("HISTE");

  c1.cd(4);
  PFJetPrecutsChargedHadronEnergyFraction->GetXaxis()->SetRangeUser(0,1);
  PFJetPrecutsChargedHadronEnergyFraction->Draw("HISTE");

  c1.cd(5);
  PFJetPrecutsChargedEmEnergyFraction->GetXaxis()->SetRangeUser(0,1);
  PFJetPrecutsChargedEmEnergyFraction->Draw("HISTE");

  c1.cd(6);
  PFJetPrecutsChargedMultiplicity->GetXaxis()->SetRangeUser(0,50);
  PFJetPrecutsChargedMultiplicity->Draw("HISTE");

  //--------------------------
  TCanvas c2;  
  c2.Divide(2,1);
  
  c2.cd(1)->SetLogy();
  PFJetPrecutsPt->Draw("");
  PFJetIDPt->SetLineColor(2);
  PFJetIDPt->Draw("sames");

  TLegend *legend2_1 = new TLegend(0.3,0.7,0.7,0.8);
  legend2_1->SetBorderSize(1);
  legend2_1->SetFillColor(0);
  legend2_1->AddEntry(PFJetPrecutsPt,"All HLT PFjets","l");
  legend2_1->AddEntry(PFJetIDPt,"With Loose JetID","l");
  legend2_1->Draw();

  c2.cd(2);
  PFJetPrecutsEta->GetXaxis()->SetRangeUser(-4.5,4.5);
  PFJetPrecutsEta->Draw("");
  PFJetIDEta->SetLineColor(2);
  PFJetIDEta->Draw("sames");

  TLegend *legend2_1 = new TLegend(0.3,0.7,0.7,0.8);
  legend2_1->SetBorderSize(1);
  legend2_1->SetFillColor(0);
  legend2_1->AddEntry(PFJetPrecutsEta,"All HLT PFjets","l");
  legend2_1->AddEntry(PFJetIDEta,"With Loose JetID","l");
  legend2_1->Draw();
   
  //--------------------------
  TCanvas c3;  
  c3.Divide(2,1);

  c3.cd(1);   
  CaloCorrJetPrecutsEnergyFractionEm->GetXaxis()->SetRangeUser(-0.5,1.5);
  CaloCorrJetPrecutsEnergyFractionEm->Draw("HISTE");

  c3.cd(2);   
  CaloCorrJetPrecutsEnergyFractionHadronic->GetXaxis()->SetRangeUser(-0.5,1.5);
  CaloCorrJetPrecutsEnergyFractionHadronic->Draw("HISTE");

  //--------------------------
  TCanvas c4;  
  c4.Divide(2,1);
  
  c4.cd(1)->SetLogy();
  CaloCorrJetPrecutsPt->Draw("");

  c4.cd(2);
  CaloCorrJetPrecutsEta->GetXaxis()->SetRangeUser(-4.5,4.5);
  CaloCorrJetPrecutsEta->Draw("");

}
