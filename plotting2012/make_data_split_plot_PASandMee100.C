//------------------------------------------------------------------------
// User has to set this information
//------------------------------------------------------------------------

const int n_run_periods = 38;

float run_lumi[n_run_periods] = { 511.298, 379.311, 504.478, 556.008, 545.397, 540.304, 586.352, 510.617, 551.351, 558.733, 
				  75.434, 517.230, 547.557, 511.426, 565.942, 532.062, 527.294, 584.215, 556.217, 500.730,  
				  509.010, 512.555, 515.621, 516.810, 560.513, 593.284, 533.723, 524.672, 516.731, 507.937, 
				  574.634, 500.076, 536.876, 531.969, 564.997, 598.525, 553.151, 176.608 };                 

int run_max[n_run_periods] = { 191718, 193621, 194223, 194533, 194912, 195304, 195552, 195930, 196239, 196453, 
			       196531, 198941, 199319, 199571, 199812, 200042, 200369, 200991, 201278, 201705, 
			       202045, 202209, 202477, 203002, 204250, 204601, 205339, 205694, 206187, 206389, 
			       206512, 206745, 207214, 207372, 207515, 208307, 208487, 208686 };

int run_min[n_run_periods] = { 190645, 191720, 193834, 194224, 194619, 194914, 195378, 195633, 195937, 196249, 
			       196495, 198049, 198954, 199336, 199572, 199833, 200049, 200381, 200992, 201554, 
			       201706, 202054, 202237, 202478, 203894, 204511, 205086, 205344, 205718, 206188, 
			       206391, 206513, 206859, 207217, 207397, 207517, 208339, 208509 };

bool rescale = true;
double y_axis_minimum = 0.;
double y_axis_maximum = 150.;
double effective_lumi = 500.;
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


void make_data_split_plot_PASandMee100(){

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

  char save_name [200]; sprintf(save_name , "nEventsPassing_PASandMee100.png");
  char hist_name [100]; sprintf(hist_name , "histo1D__DATA__split_PASandMee100");
  
  TFile * file = new TFile ( file_path );
  TH1F  * hist = (TH1F*) file -> Get(hist_name);

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

  // hist -> GetXaxis() -> SetRangeUser(-0.5, 19.5);
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
  

  hist -> Fit ( "funcTotal", "R"  );
  hist -> Fit ( "func2012A", "R+" );
  hist -> Fit ( "func2012B", "R+" );
  hist -> Fit ( "func2012C", "R+" );
  hist -> Fit ( "func2012D", "R+" );
  
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
    
    st -> SetX1NDC ( 0.619253 );
    st -> SetY1NDC ( 0.737288 );
    st -> SetX2NDC ( 0.849138 );
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

  char plot_title[100]; sprintf(plot_title, "Number of events passing preselection and M(ee) > 100");
  TLatex * texB = new TLatex ();
  texB -> SetTextSize( 0.05 );
  texB -> DrawLatex ( x_min, 1.05 * ( y_max - y_min) + y_min, plot_title );

  //------------------------------------------------------------------------
  // Save the canvas
  //------------------------------------------------------------------------
  
  c -> SaveAs ( save_name );

}
