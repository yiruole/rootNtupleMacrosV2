
void retrieve(){
  
  //TFile * file = new TFile ("Mej_minmax_PAS_eejj_WZSherpa_noNSigma.root");
  TFile * file = new TFile ("sT_PAS_eejj_WZSherpa_noNSigma.root");

  bool verbose = false;


  TCanvas * canvas = (TCanvas*) file -> Get ("c1_n2");
  TPad * pad = (TPad*) canvas -> GetPrimitive("pad1");
  TList * primitives = pad -> GetListOfPrimitives();
  TIter next(primitives);
  
  TObject * object = 0;

  bool found_zjets = false;
  bool found_other = false;
  bool found_qcd   = false;
  bool found_ttbar = false;
  bool found_data  = false;
  bool found_lq    = false;

  bool made_print = false;
  
  while ( object = next() ){
    
    if ( ! object -> InheritsFrom(TH1F::Class())) continue;
    
    TH1F * hist = (TH1F*) object;

    if ( !made_print ) { 
      canvas -> Update();
      std::cout << "// Number of bins = " << hist -> GetNbinsX() << std::endl;
      std::cout << "// x axis minimum = " << hist -> GetXaxis() -> GetXmin() << std::endl;
      std::cout << "// x axis maximum = " << hist -> GetXaxis() -> GetXmax() << std::endl;
      std::cout << std::endl;
      made_print = true;
    }
    
    char name[100];
    
    if      ( hist -> GetFillStyle () == 3004 ) {
      if ( found_zjets ) continue;
      found_zjets = true;
      sprintf(name, "Z/#gamma* + jets");
    }

    else if ( hist -> GetFillStyle () == 3008 ) {
      if ( found_other ) continue;
      found_other = true;
      sprintf(name, "Other backgrounds");
    }
    
    else if ( hist -> GetFillColor () == kCyan) { 
      if ( found_qcd ) continue;
      found_qcd = true;
      sprintf(name, "QCD multijets");
    }
    
    else if ( hist -> GetFillColor () == kBlue) {
      if ( found_ttbar ) continue;
      found_ttbar = true;
      sprintf(name, "t#bar{t} [Data driven e-#mu]");
    }
    
    else if ( hist -> GetLineColor () == 1    ) {
      if ( found_data ) continue;
      found_data = true;
      sprintf(name, "Data, 5.0 fb^{-1}");
    }
    
    else if ( hist -> GetLineColor () == 8    ) {
      if ( found_lq ) continue;
      found_lq = true;
      sprintf(name, "M_{LQ} = 400 GeV");
    }
    
    else continue;
    
    std::cout << "\"" << name << "\":" << std::endl;
    if ( verbose ) { 
      std::cout << "Line   = " << hist -> GetLineColor () << std::endl;
      std::cout << "Color  = " << hist -> GetFillColor () << std::endl;
      std::cout << "Style  = " << hist -> GetFillStyle () << std::endl;
      std::cout << "NBins  = " << hist -> GetNbinsX() << std::endl;
      std::cout << "Edge1  = " << hist -> GetBinLowEdge (1 ) << std::endl;
      std::cout << "Edge91 = " << hist -> GetBinLowEdge (91) << std::endl;
      std::cout << "Cont2  = " << hist -> GetBinContent (2 ) << std::endl;
      std::cout << "Cont60 = " << hist -> GetBinContent (60) << std::endl;
    }

    int last_bin = hist -> GetNbinsX();
    char tmp_bin_content[100]; sprintf(tmp_bin_content, "float bin_content[%d] = { ", last_bin);
    
    std::string bin_content (tmp_bin_content);
    


    for (int i = 1; i < (last_bin + 1); ++i){
      double content = hist -> GetBinContent (i);
      char content_char[100]; 
      if (i != last_bin ) sprintf(content_char,"%1.2f, ", content );
      else                sprintf(content_char,"%1.2f " , content );
      bin_content = bin_content + std::string(content_char) ;
    }
    bin_content = bin_content + std::string ("};");
    std::cout << bin_content << std::endl;
    std::cout << std::endl;

  }
  
  
  

}
