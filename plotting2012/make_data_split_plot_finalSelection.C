//------------------------------------------------------------------------
// User has to set this information
//------------------------------------------------------------------------

const int n_run_periods = 20;

int run_max[20] = { 193621, 194533, 195304, 195915, 196452, 196531, 199319, 199804, 200244, 201196, 
                    202014, 202305, 203002, 204577, 205666, 206246, 206598, 207273, 207905, 208686 };



int run_min[20] = { 190645, 193834, 194619, 195378, 195916, 196453, 198049, 199336, 199812, 200245, 
                    201197, 202016, 202314, 203894, 204599, 205667, 206257, 206605, 207279, 207920 };

float run_lumi[20] = { 890.608, 1060.486, 1085.701, 1028.249, 1065.344, 188.894, 1064.787, 1005.727, 1039.579, 1033.637, 
                       1012.554, 1014.667, 725.718, 1013.083, 1020.154, 1052.137, 1041.412, 1061.761, 1092.669, 992.480 };

bool rescale = true;
double y_axis_minimum = 0.;
double effective_lumi = 1000.;
const char* file_path  = "/mnt/lxplus/work/LQDATA//eejj_analysis/eejj//output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root";


TH1F* addZeroPoints(TH1F* hist){
  
  TH1F* clone_hist = hist -> Clone();

  int nbins = clone_hist -> GetNbinsX();
  
  for (int i = 0; i <= nbins; ++i){
    double old_content = clone_hist -> GetBinContent ( i );
    double old_error   = clone_hist -> GetBinError   ( i );
    double new_content = old_content;
    double new_error   = old_error;
    
    if ( old_content == 0.0 ) new_error = 1.0;
    
    clone_hist -> SetBinContent ( i, new_content );
    clone_hist -> SetBinError   ( i, new_error   );
  }
  return clone_hist;
}

void make_data_split_plot_finalSelection( int cut_level, double y_axis_maximum ) { 

  //------------------------------------------------------------------------
  // Find limits of run eras
  //------------------------------------------------------------------------
  
  int era2012A_min_bin, era2012B_min_bin, era2012C_min_bin, era2012D_min_bin;
  int era2012A_max_bin, era2012B_max_bin, era2012C_max_bin, era2012D_max_bin;
  
  for ( int i_run_period = 0; i_run_period < n_run_periods; ++i_run_period){
    int tmp_run_min = run_min[i_run_period];
    int tmp_run_max = run_max[i_run_period];

    if ( tmp_run_min == 190645 ) era2012A_min_bin = i_run_period;
    if ( tmp_run_max == 193621 ) era2012A_max_bin = i_run_period;
    if ( tmp_run_min == 193834 ) era2012B_min_bin = i_run_period;
    if ( tmp_run_max == 196531 ) era2012B_max_bin = i_run_period;
    if ( tmp_run_min == 198049 ) era2012C_min_bin = i_run_period;
    if ( tmp_run_max == 203002 ) era2012C_max_bin = i_run_period;
    if ( tmp_run_min == 203894 ) era2012D_min_bin = i_run_period;
    if ( tmp_run_max == 208686 ) era2012D_max_bin = i_run_period;

  }

  //------------------------------------------------------------------------
  // Get the file and histogram
  //------------------------------------------------------------------------

  char save_name [200]; sprintf(save_name , "nEventsPassing_LQM%d.png", cut_level);
  char hist_name [100]; sprintf(hist_name , "histo1D__DATA__split_1fb_LQ%d", cut_level);
  
  TFile * file = new TFile ( file_path );
  TH1F  * hist = (TH1F*) file -> Get(hist_name);

  hist = addZeroPoints ( hist );


  //------------------------------------------------------------------------
  // Count events per run period
  //------------------------------------------------------------------------

  double era_2012A_lumi = 0.89;
  double era_2012B_lumi = 4.43;
  double era_2012C_lumi = 6.88;
  double era_2012D_lumi = 7.23;
  double total_lumi = 19.5;

  double era2012A_integral = hist -> Integral ( era2012A_min_bin + 1, era2012A_max_bin + 1) ;
  double era2012B_integral = hist -> Integral ( era2012B_min_bin + 1, era2012B_max_bin + 1) ;
  double era2012C_integral = hist -> Integral ( era2012C_min_bin + 1, era2012C_max_bin + 1) ;
  double era2012D_integral = hist -> Integral ( era2012D_min_bin + 1, era2012D_max_bin + 1) ;
  double total_integral    = ( era2012A_integral +
			       era2012B_integral +
			       era2012C_integral +
			       era2012D_integral );


  double era2012A_error = TMath::Sqrt ( era2012A_integral  + 1) ; 
  double era2012B_error = TMath::Sqrt ( era2012B_integral  + 1) ; 
  double era2012C_error = TMath::Sqrt ( era2012C_integral  + 1) ; 
  double era2012D_error = TMath::Sqrt ( era2012D_integral  + 1) ; 
  double total_error    = TMath::Sqrt ( total_integral + 1 ) ;

  double era2012A_xsection = era2012A_integral / era_2012A_lumi;
  double era2012B_xsection = era2012B_integral / era_2012B_lumi;
  double era2012C_xsection = era2012C_integral / era_2012C_lumi;
  double era2012D_xsection = era2012D_integral / era_2012D_lumi;
  double total_xsection = total_integral / total_lumi;
  
  double era2012A_xsection_error = era2012A_error / era_2012A_lumi;
  double era2012B_xsection_error = era2012B_error / era_2012B_lumi;
  double era2012C_xsection_error = era2012C_error / era_2012C_lumi;
  double era2012D_xsection_error = era2012D_error / era_2012D_lumi;
  double total_xsection_error = total_error / total_lumi;


  std::cout << "\n\n" << std::endl;
  
  std::cout << "N(2012A) = " << era2012A_integral <<  std::endl;
  std::cout << "N(2012B) = " << era2012B_integral <<  std::endl;
  std::cout << "N(2012C) = " << era2012C_integral <<  std::endl;
  std::cout << "N(2012D) = " << era2012D_integral <<  std::endl;
  std::cout << "Total    = " << total_integral    <<  std::endl;

  std::cout << "\n\n" << std::endl;

  std::cout.precision(1);
  std::cout << std::fixed << "2012A xsection = " << era2012A_xsection  << " +/- " << era2012A_xsection_error << " fb" << std::endl;
  std::cout << std::fixed << "2012B xsection = " << era2012B_xsection  << " +/- " << era2012B_xsection_error << " fb" << std::endl;
  std::cout << std::fixed << "2012C xsection = " << era2012C_xsection  << " +/- " << era2012C_xsection_error << " fb" << std::endl;
  std::cout << std::fixed << "2012D xsection = " << era2012D_xsection  << " +/- " << era2012D_xsection_error << " fb" << std::endl;
  std::cout << std::fixed << "Total xsection = " << total_xsection << " +/- " << total_xsection_error << " fb " << std::endl;

  std::cout << "\n\n" << std::endl;
  
  bool verbose = false;
  if ( verbose ) { 
    std::cout << "2012A, bin-by-bin:" << std::endl;
    double total_2012A = 0;
    for (int ibin = era2012A_min_bin + 1; ibin <= era2012A_max_bin + 1; ++ibin ){
      double content  = hist -> Integral (ibin, ibin );
      double low_edge = hist -> GetBinLowEdge ( ibin ) ;
      double high_edge = hist -> GetBinLowEdge ( ibin + 1 ) ;
      total_2012A += content;
      std::cout << "\t" << ibin << ": " << content << ", low edge = " << low_edge << ", high edge = " << high_edge << std::endl;
    }
    std::cout << "2012A total = " << total_2012A << std::endl;
    
    std::cout << "2012B, bin-by-bin:" << std::endl;
    double total_2012B = 0;
    for (int ibin = era2012B_min_bin + 1; ibin <= era2012B_max_bin + 1; ++ibin ){
      double content  = hist -> Integral (ibin, ibin );
      double low_edge = hist -> GetBinLowEdge ( ibin ) ;
      double high_edge = hist -> GetBinLowEdge ( ibin + 1 ) ;
      total_2012B += content;
      std::cout << "\t" << ibin << ": " << content << ", low edge = " << low_edge << ", high edge = " << high_edge << std::endl;
    }
    std::cout << "2012B total = " << total_2012B << std::endl;
    
    std::cout << "2012C, bin-by-bin:" << std::endl;
    double total_2012C = 0;
    for (int ibin = era2012C_min_bin + 1; ibin <= era2012C_max_bin + 1; ++ibin ){
      double content  = hist -> Integral (ibin, ibin );
      double low_edge = hist -> GetBinLowEdge ( ibin ) ;
      double high_edge = hist -> GetBinLowEdge ( ibin + 1 ) ;
    total_2012C += content;
    std::cout << "\t" << ibin << ": " << content << ", low edge = " << low_edge << ", high edge = " << high_edge << std::endl;
    }
    std::cout << "2012C total = " << total_2012C << std::endl;
    
    std::cout << "2012D, bin-by-bin:" << std::endl;
    double total_2012D = 0;
    for (int ibin = era2012D_min_bin + 1; ibin <= era2012D_max_bin + 1; ++ibin ){
      double content  = hist -> Integral (ibin, ibin );
      double low_edge = hist -> GetBinLowEdge ( ibin ) ;
      double high_edge = hist -> GetBinLowEdge ( ibin + 1 ) ;
      total_2012D += content;
      std::cout << "\t" << ibin << ": " << content << ", low edge = " << low_edge << ", high edge = " << high_edge << std::endl;
    }
    std::cout << "2012D total = " << total_2012D << std::endl;
  }


  //------------------------------------------------------------------------
  // Find where era splits occur on the canvas
  //------------------------------------------------------------------------

  double era2012A_min = hist -> GetBinLowEdge ( era2012A_min_bin  + 1);
  double era2012A_max = hist -> GetBinLowEdge ( era2012B_min_bin  + 1);
  
  double era2012B_min = hist -> GetBinLowEdge ( era2012B_min_bin  + 1);
  double era2012B_max = hist -> GetBinLowEdge ( era2012C_min_bin  + 1);

  double era2012C_min = hist -> GetBinLowEdge ( era2012C_min_bin  + 1);
  double era2012C_max = hist -> GetBinLowEdge ( era2012D_min_bin  + 1);

  double era2012D_min = hist -> GetBinLowEdge ( era2012D_min_bin  + 1);
  double era2012D_max = hist -> GetBinLowEdge ( era2012D_max_bin  + 2);
  
  //------------------------------------------------------------------------
  // Global style settings
  //------------------------------------------------------------------------

  gStyle -> SetOptStat(""); // Do not show statistics info
  gStyle -> SetOptFit ()  ; // Do show fit info
  //gStyle -> SetOptFit (0); // Do not show fit info
  
  //------------------------------------------------------------------------
  // Style histogram
  //------------------------------------------------------------------------
  
  for (int i = 0; i < n_run_periods; ++i) {
    int bin = hist -> FindBin (i);
    char bin_label[200]; sprintf(bin_label,"%d - %d, L = %1.1f pb^{-1}", run_min[i], run_max[i], run_lumi[i] );
    hist -> GetXaxis() -> SetBinLabel (bin, bin_label);
  }

  hist -> SetMaximum(y_axis_maximum);
  hist -> SetMinimum(y_axis_minimum);

  hist -> GetXaxis() -> SetRangeUser(-0.5, 18.5);
  hist -> GetXaxis() -> LabelsOption("v");

  hist -> SetLineColor (kBlack);
  hist -> SetLineWidth (2.0);

  //------------------------------------------------------------------------
  // Label histogram axes
  //------------------------------------------------------------------------
  
  char y_axis_label[200]; 
  if ( rescale ) sprintf(y_axis_label, "N(Pass) #times %1.1f fb^{-1} / L_{int}", effective_lumi / 1000. );
  else           sprintf(y_axis_label, "N(Pass)");
  hist -> GetYaxis() -> SetTitleOffset ( 1.2 );
  hist -> GetYaxis() -> SetTitle( y_axis_label );

  char x_axis_label[200]; sprintf(x_axis_label, "Run period");
  hist -> GetXaxis() -> SetTitleOffset ( 6.25 );
  hist -> GetXaxis() -> SetTitle( x_axis_label );

  //------------------------------------------------------------------------
  // Style canvas
  //------------------------------------------------------------------------

  TCanvas * c = new TCanvas();
  c -> SetBottomMargin (0.38);  
  c -> SetGridx();
  c -> SetGridy();

  //------------------------------------------------------------------------
  // Style fit function
  //------------------------------------------------------------------------

  TF1 * func2012A = new TF1 ("func2012A","[0]", era2012A_min, era2012A_max );
  TF1 * func2012B = new TF1 ("func2012B","[0]", era2012B_min, era2012B_max );
  TF1 * func2012C = new TF1 ("func2012C","[0]", era2012C_min, era2012C_max );
  TF1 * func2012D = new TF1 ("func2012D","[0]", era2012D_min, era2012D_max );
  TF1 * funcTotal = new TF1 ("funcTotal","[0] + [1] * x", era2012A_min, era2012D_max );
  
  funcTotal -> SetLineColor(kBlue);
  funcTotal -> SetParName ( 0, "Offset");
  funcTotal -> SetParName ( 1, "Slope" );
  
  //------------------------------------------------------------------------
  // Scale each hist bin so that they all have the same "effective" lumi
  //------------------------------------------------------------------------

  if ( rescale ) { 
    for (int i = 0; i < n_run_periods; ++i) {
      int bin = hist -> FindBin (i);
      
      double lumi = run_lumi[i];
      double ratio = effective_lumi / lumi;
      
      double old_content = hist -> GetBinContent ( bin ) ;
      double old_error   = hist -> GetBinError   ( bin ) ;
      double new_content = old_content * ratio;
      double new_error   = old_error   * ratio;
      
      hist -> SetBinContent ( bin, new_content ) ;
      hist -> SetBinError   ( bin, new_error   ) ;
    }
  }

  //------------------------------------------------------------------------
  // Perform the fit (draw the histogram) and update the canvas/pad
  //------------------------------------------------------------------------
  
  hist -> Fit ( "funcTotal", "QR"  );
  hist -> Fit ( "func2012A", "QR+" );
  hist -> Fit ( "func2012B", "QR+" );
  hist -> Fit ( "func2012C", "QR+" );
  hist -> Fit ( "func2012D", "QR+" );
  
  //------------------------------------------------------------------------
  // Update the histogram and get the coordinates
  //------------------------------------------------------------------------

  gPad -> Update();

  double x_min = c -> GetUxmin();
  double x_max = c -> GetUxmax();
  double y_min = c -> GetUymin();
  double y_max = c -> GetUymax();
  
  //------------------------------------------------------------------------
  // Move fit info panel
  //------------------------------------------------------------------------
  
  if ( gStyle -> GetOptFit () != 0 ){
    TPaveStats *st = (TPaveStats*) hist -> FindObject("stats");
    
    st -> SetX1NDC ( 0.640805 );
    st -> SetY1NDC ( 0.737288 );
    st -> SetX2NDC ( 0.87069  );
    st -> SetY2NDC ( 0.885593 );
  }
  
  //------------------------------------------------------------------------
  // Draw the lines separating the run eras
  //------------------------------------------------------------------------

  TLine * line_2012A  = new TLine ( era2012A_min, y_min, era2012A_min, y_max ) ;
  TLine * line_2012AB = new TLine ( era2012B_min, y_min, era2012B_min, y_max ) ;
  TLine * line_2012BC = new TLine ( era2012C_min, y_min, era2012C_min, y_max ) ;
  TLine * line_2012CD = new TLine ( era2012D_min, y_min, era2012D_min, y_max ) ;
  TLine * line_2012D  = new TLine ( era2012D_max, y_min, era2012D_max, y_max ) ;

  line_2012A  -> SetLineColor(kBlue);
  line_2012AB -> SetLineColor(kBlue);
  line_2012BC -> SetLineColor(kBlue);
  line_2012CD -> SetLineColor(kBlue);
  line_2012D  -> SetLineColor(kBlue);

  line_2012A  -> SetLineWidth(3.0);
  line_2012AB -> SetLineWidth(3.0);
  line_2012BC -> SetLineWidth(3.0);
  line_2012CD -> SetLineWidth(3.0);
  line_2012D  -> SetLineWidth(3.0);

  line_2012A  -> Draw("SAME");
  line_2012AB -> Draw("SAME");
  line_2012BC -> Draw("SAME");
  line_2012CD -> Draw("SAME");
  line_2012D  -> Draw("SAME");

  //------------------------------------------------------------------------
  // Draw fit uncertainty boxes
  //------------------------------------------------------------------------

  TBox * box_fit2012A = new TBox ( era2012A_min, func2012A -> GetParameter(0) - func2012A -> GetParError (0), era2012A_max, func2012A -> GetParameter(0) + func2012A -> GetParError (0) );
  TBox * box_fit2012B = new TBox ( era2012B_min, func2012B -> GetParameter(0) - func2012B -> GetParError (0), era2012B_max, func2012B -> GetParameter(0) + func2012B -> GetParError (0) );
  TBox * box_fit2012C = new TBox ( era2012C_min, func2012C -> GetParameter(0) - func2012C -> GetParError (0), era2012C_max, func2012C -> GetParameter(0) + func2012C -> GetParError (0) );
  TBox * box_fit2012D = new TBox ( era2012D_min, func2012D -> GetParameter(0) - func2012D -> GetParError (0), era2012D_max, func2012D -> GetParameter(0) + func2012D -> GetParError (0) );
  
  box_fit2012A -> SetFillColor(kRed);
  box_fit2012A -> SetLineColor(kRed);
  box_fit2012A -> SetFillStyle(3004);

  box_fit2012B -> SetFillColor(kRed);
  box_fit2012B -> SetLineColor(kRed);
  box_fit2012B -> SetFillStyle(3004);

  box_fit2012C -> SetFillColor(kRed);
  box_fit2012C -> SetLineColor(kRed);
  box_fit2012C -> SetFillStyle(3004);

  box_fit2012D -> SetFillColor(kRed);
  box_fit2012D -> SetLineColor(kRed);
  box_fit2012D -> SetFillStyle(3004);
  

  box_fit2012A -> Draw("SAME");
  box_fit2012B -> Draw("SAME");
  box_fit2012C -> Draw("SAME");
  box_fit2012D -> Draw("SAME");
  
  //------------------------------------------------------------------------
  // Label the lines separating the run eras
  //------------------------------------------------------------------------
  
  TLatex * texA = new TLatex ();
  texA -> SetTextAngle(90.);
  texA -> SetTextSize( 0.03);
  texA -> DrawLatex( line_2012AB -> GetX1() - 0.1 , 0.05 * ( y_max - y_min) + y_min , "2012A" );
  texA -> DrawLatex( line_2012BC -> GetX1() - 0.1 , 0.05 * ( y_max - y_min) + y_min , "2012B" );
  texA -> DrawLatex( line_2012CD -> GetX1() - 0.1 , 0.05 * ( y_max - y_min) + y_min , "2012C" );
  texA -> DrawLatex( line_2012D  -> GetX1() - 0.1 , 0.05 * ( y_max - y_min) + y_min , "2012D" );

  char plot_title[100]; sprintf(plot_title, "Number of events passing M(LQ) = %d selection", cut_level );
  TLatex * texB = new TLatex ();
  texB -> SetTextSize( 0.05 );
  texB -> DrawLatex ( x_min, 1.05 * ( y_max - y_min) + y_min, plot_title );

  //------------------------------------------------------------------------
  // Save the canvas
  //------------------------------------------------------------------------
  
  c -> SaveAs ( save_name );

}
