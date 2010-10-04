{

  gROOT->Reset();
  gStyle->SetPalette(1);
 

  TFile f_all("/Users/santanas/Desktop/QCDresults/2.9pb-1_v9_nod1d2.root");
  TFile f_QCD_data("/Users/santanas/Desktop/QCDresults/2.9pb-1_v9_nod1d2_QCD.root");

  TEllipse el_d1(TMath::Pi(),0.,.6,0,0,90,90);
  el_d1.SetLineColor(kRed);
  el_d1.SetLineWidth(3);
  el_d1.SetFillStyle(0);

  TEllipse el_d2(0,TMath::Pi(),.6,0,180,270,90);
  el_d2.SetLineColor(kRed);
  el_d2.SetLineWidth(3);
  el_d2.SetFillStyle(0);


  TCanvas c_presel;
  c_presel.Divide(2,2);

  f_QCD_data.cd();

  c_presel.cd(1);
  histo2D__DATA__h2_DeltaPhiMETEle_vs_MET1stJet->DrawClone("colz");
  el_d1.Draw();
  el_d2.Draw();

  f_all.cd();

  c_presel.cd(2);
  histo2D__LQenujj_M200__h2_DeltaPhiMETEle_vs_MET1stJet->DrawClone("colz");
  el_d1.Draw();
  el_d2.Draw();

  c_presel.cd(3);
  histo2D__TTbar_Madgraph__h2_DeltaPhiMETEle_vs_MET1stJet->DrawClone("colz");
  el_d1.Draw();
  el_d2.Draw();

  c_presel.cd(4);
  histo2D__WJetAlpgen__h2_DeltaPhiMETEle_vs_MET1stJet->DrawClone("colz");
  el_d1.Draw();
  el_d2.Draw();

  c_presel.SaveAs("deltaPhi2_presel.gif");


  TCanvas c_MTenu;
  c_MTenu.Divide(2,2);

  f_QCD_data.cd();

  c_MTenu.cd(1);
  histo2D__DATA__h2_DeltaPhiMETEle_vs_MET1stJet__MTenu->DrawClone("colz");
  el_d1.Draw();
  el_d2.Draw();

  f_all.cd();

  c_MTenu.cd(2);
  histo2D__LQenujj_M200__h2_DeltaPhiMETEle_vs_MET1stJet__MTenu->DrawClone("colz");
  el_d1.Draw();
  el_d2.Draw();

  c_MTenu.cd(3);
  histo2D__TTbar_Madgraph__h2_DeltaPhiMETEle_vs_MET1stJet__MTenu->DrawClone("colz");
  el_d1.Draw();
  el_d2.Draw();

  c_MTenu.cd(4);
  histo2D__WJetAlpgen__h2_DeltaPhiMETEle_vs_MET1stJet__MTenu->DrawClone("colz");
  el_d1.Draw();
  el_d2.Draw();

  c_MTenu.SaveAs("deltaPhi2_MTenu.gif");

  TCanvas c_sT;
  c_sT.Divide(2,2);

  f_QCD_data.cd();

  c_sT.cd(1);
  histo2D__DATA__h2_DeltaPhiMETEle_vs_MET1stJet__sT->DrawClone("colz");
  el_d1.Draw();
  el_d2.Draw();

  f_all.cd();

  c_sT.cd(2);
  histo2D__LQenujj_M200__h2_DeltaPhiMETEle_vs_MET1stJet__sT->DrawClone("colz");
  el_d1.Draw();
  el_d2.Draw();

  c_sT.cd(3);
  histo2D__TTbar_Madgraph__h2_DeltaPhiMETEle_vs_MET1stJet__sT->DrawClone("colz");
  el_d1.Draw();
  el_d2.Draw();

  c_sT.cd(4);
  histo2D__WJetAlpgen__h2_DeltaPhiMETEle_vs_MET1stJet__sT->DrawClone("colz");
  el_d1.Draw();
  el_d2.Draw();

  c_sT.SaveAs("deltaPhi2_sT.gif");

}
