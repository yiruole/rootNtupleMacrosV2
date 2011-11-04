#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TH2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TRandom3.h>

//-----------------------------
//-----------------------------


analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  //STDOUT("analysisClass::analysisClass() was called");
}

analysisClass::~analysisClass()
{
  //STDOUT("analysisClass::~analysisClass() was called");
}

void analysisClass::Loop()
{
  //STDOUT("analysisClass::Loop() begins");

  if (fChain == 0) return;

   //--------------------------------------------------------------------------
   // Decide which plots to save (default is to save everything)
   //--------------------------------------------------------------------------
   
   fillSkim                         (  true  ) ;
   fillAllPreviousCuts              ( !true  ) ;
   fillAllOtherCuts                 ( !true  ) ;
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;

  ////////////////////// User's code to get preCut values - BEGIN ///////////////

  // Jets
  double pfjetPtCut = getPreCutValue1("pfjetPtCut");
  double pfjetEtaCut = getPreCutValue1("pfjetEtaCut");
  double pfjetIDloose = getPreCutValue1("pfjetIDloose");
  double pfjetIDtight = getPreCutValue1("pfjetIDtight");

  double caloRawjetPtCut = getPreCutValue1("caloRawjetPtCut");
  double caloRawjetEtaCut = getPreCutValue1("caloRawjetEtaCut");
  double caloRawjetIDloose = getPreCutValue1("caloRawjetIDloose");
  double caloRawjetIDtight = getPreCutValue1("caloRawjetIDtight");

  double caloCorrjetPtCut = getPreCutValue1("caloCorrjetPtCut");
  double caloCorrjetEtaCut = getPreCutValue1("caloCorrjetEtaCut");
  double caloCorrjetIDloose = getPreCutValue1("caloCorrjetIDloose");
  double caloCorrjetIDtight = getPreCutValue1("caloCorrjetIDtight");

  // Vertices
  double vertexMaxAbsZ = getPreCutValue1("vertexMaxAbsZ");
  double vertexMaxd0 = getPreCutValue1("vertexMaxd0");  
  
  ////////////////////// User's code to get preCut values - END /////////////////

  ////////////////////// User's code to book histos - BEGIN ///////////////////////

  // Random number generator
  CreateUserTH1D( "PFJetPrecutsPt"              ,    200,  0, 2000   );
  CreateUserTH1D( "PFJetPrecutsEta"             ,    100, -6, 6      );
  CreateUserTH1D( "PFJetPrecutsNeutralHadronEnergyFraction"             ,    200,  -2, 2   );
  CreateUserTH1D( "PFJetPrecutsNeutralEmEnergyFraction"                 ,    200,  -2, 2   );
  CreateUserTH1D( "PFJetPrecutsNConstituents"                           ,    201,  -0.5, 200.5   );
  CreateUserTH1D( "PFJetPrecutsChargedHadronEnergyFraction"             ,    200,  -2, 2   );
  CreateUserTH1D( "PFJetPrecutsChargedMultiplicity"                     ,    201,  -0.5, 200.5   );
  CreateUserTH1D( "PFJetPrecutsChargedEmEnergyFraction"                 ,    200,  -2, 2   );

  CreateUserTH1D( "PFJetIDPt"                   ,    200,  0, 2000   );
  CreateUserTH1D( "PFJetIDEta"                  ,    100, -6, 6      );

  CreateUserTH1D( "CaloRawJetAllPt"              ,    200,  0, 2000   );
  CreateUserTH1D( "CaloRawJetAllEta"             ,    100, -6, 6      );

  CreateUserTH1D( "CaloCorrJetAllPt"             ,    200,  0, 2000   );
  CreateUserTH1D( "CaloCorrJetAllEta"            ,    100, -6, 6      );

  ////////////////////// User's code to book histos - END ///////////////////////

  Long64_t nentries = fChain->GetEntriesFast();
  STDOUT("analysisClass::Loop(): nentries = " << nentries);

  ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
  ////// If the root version is updated and rootNtupleClass regenerated,     /////
  ////// these lines may need to be updated.                                 /////
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) { // Begin of loop over events
  // for (Long64_t jentry=0; jentry<1000;jentry++) { // Begin of loop over events
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry < 10 || jentry%1000 == 0) STDOUT("analysisClass::Loop(): jentry = " << jentry);
    // if (Cut(ientry) < 0) continue;

    //int NPILEUP_AVE = int( (PileUpInteractions->at(0) + PileUpInteractions->at(1) + PileUpInteractions->at(2))/3 );
    //int NPILEUP_FINAL = min( NPILEUP_AVE , 25 );
    //double event_weight = getPileupWeight ( NPILEUP_FINAL, isData ) ;
    //double event_weight = getPileupWeight ( min(PileUpInteractions->at(1),25), isData ) ;

    ////////////////////// User's code to be done for every event - BEGIN ///////////////////////

    ////////////////////// Reco Object Collections ///////////////////////

    //## PF Jets
    vector<int> v_idx_pfjet_PtEtaCut;
    vector<int> v_idx_pfjet_PtEtaCut_ID;

    // Precuts
    for(int ijet=0; ijet<HLTPFJetPt->size(); ijet++)
      {
	//pT/eta pre-cuts on jets
	if ( HLTPFJetPt->at(ijet) < pfjetPtCut ) continue;
	if ( fabs( HLTPFJetEta->at(ijet) ) > pfjetEtaCut ) continue;
	v_idx_pfjet_PtEtaCut.push_back(ijet);
	
	//histograms
	FillUserTH1D( "PFJetPrecutsPt"                                      ,    HLTPFJetPt->at(ijet)  ); 
	FillUserTH1D( "PFJetPrecutsEta"                                     ,    HLTPFJetEta->at(ijet) ); 	
	FillUserTH1D( "PFJetPrecutsNeutralHadronEnergyFraction"             ,    HLTPFJetNeutralHadronEnergyFraction->at(ijet)  ); 
	FillUserTH1D( "PFJetPrecutsNeutralEmEnergyFraction"                 ,    HLTPFJetNeutralEmEnergyFraction->at(ijet)     );  
	FillUserTH1D( "PFJetPrecutsNConstituents"                           ,    HLTPFJetNConstituents->at(ijet)  ); 
	FillUserTH1D( "PFJetPrecutsChargedHadronEnergyFraction"             ,    HLTPFJetChargedHadronEnergyFraction->at(ijet)  ); 
	FillUserTH1D( "PFJetPrecutsChargedMultiplicity"                     ,    HLTPFJetChargedMultiplicity->at(ijet)  ); 
	FillUserTH1D( "PFJetPrecutsChargedEmEnergyFraction"                 ,    HLTPFJetChargedEmEnergyFraction->at(ijet)  ); 
      }

    // Full selection
    for(int ijet=0; ijet<v_idx_pfjet_PtEtaCut.size(); ijet++) 
      {
	bool passjetID = true;
	
	if( pfjetIDloose && !pfjetIDtight )
	  {
	    passjetID = HLTPFJetPassLooseID->at(v_idx_pfjet_PtEtaCut[ijet]);
	  }

	if( pfjetIDtight )
	  {
	    passjetID = HLTPFJetPassTightID->at(v_idx_pfjet_PtEtaCut[ijet]);
	  }

	//passjetID = true;

	if( passjetID == true )
	  {
	    v_idx_pfjet_PtEtaCut_ID.push_back(v_idx_pfjet_PtEtaCut[ijet]);

	    //histograms
	    FillUserTH1D( "PFJetIDPt"                  ,    HLTPFJetPt->at(v_idx_pfjet_PtEtaCut[ijet]) ); 
	    FillUserTH1D( "PFJetIDEta"                 ,    HLTPFJetEta->at(v_idx_pfjet_PtEtaCut[ijet]) ); 
	  }

      } // End loop over jets


    //## Pixel Vertices
    vector<int> v_idx_vertex_good;
      
    // Full selection
    for(int ivertex = 0; ivertex<HLTPixelVertexZCoord->size(); ivertex++){
      
      double vertex_rho = sqrt( 
			       HLTPixelVertexXCoord->at(ivertex)*HLTPixelVertexXCoord->at(ivertex) +
			       HLTPixelVertexYCoord->at(ivertex)*HLTPixelVertexYCoord->at(ivertex) 
			       );
      
      if ( HLTPixelVertexIsValid->at(ivertex) 
	   && fabs( HLTPixelVertexZCoord->at(ivertex) ) <= vertexMaxAbsZ
	   && vertex_rho <= vertexMaxd0 
	   )
	{
	  v_idx_vertex_good.push_back(ivertex);
	}
    }    

    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();
    //fillVariableWithValue( "PassJSON", passJSON(run, ls, isData) );    

    // Set the value of the variableNames listed in the cutFile to their current value

    //event info
    fillVariableWithValue( "isData", isData ) ;
    fillVariableWithValue( "bunch", bunch ) ;
    fillVariableWithValue( "event", event ) ;
    fillVariableWithValue( "ls", ls ) ;
    fillVariableWithValue( "orbit", orbit ) ;
    fillVariableWithValue( "run", run ) ;

//     // Trigger (HLT)
//     HLT_HT350                       -inf         +inf               -               -               0       2 -0.5 1.5	      SAVE
//       HLT_FatJetMass300               -inf         +inf               -               -               0       2 -0.5 1.5	      SAVE
//       HLT_FatJetMass400               -inf         +inf               -               -               0       2 -0.5 1.5	      SAVE


//     if(isData==true)
//       {
//       }
//     else
//       {
//       }

    //Event filters at RECO level
    int isPrimaryVertex = 0;
    if ( v_idx_vertex_good.size() > 0 )
      isPrimaryVertex = 1;
    fillVariableWithValue( "PassPrimaryVertex", isPrimaryVertex ) ;

    // nVertex
    fillVariableWithValue( "nVertex", HLTPixelVertexZCoord->size() ) ;
    fillVariableWithValue( "nVertex_good", v_idx_vertex_good.size() ) ;

    // 1st Vertex
    if( v_idx_vertex_good.size() >= 1 )
      {
	fillVariableWithValue( "Vertex1_X", HLTPixelVertexXCoord->at(v_idx_vertex_good[0]) ) ;
	fillVariableWithValue( "Vertex1_Y", HLTPixelVertexYCoord->at(v_idx_vertex_good[0]) ) ;
	fillVariableWithValue( "Vertex1_Z", HLTPixelVertexZCoord->at(v_idx_vertex_good[0]) ) ;
      }
 
    // nJet
    fillVariableWithValue( "nPFJet", v_idx_pfjet_PtEtaCut_ID.size() ) ;

    // 1st jet 
    if( v_idx_pfjet_PtEtaCut_ID.size() >= 1 )        
      {
	fillVariableWithValue( "PFJet1_Pt", HLTPFJetPt->at(v_idx_pfjet_PtEtaCut_ID[0]) );
	fillVariableWithValue( "PFJet1_Energy", HLTPFJetEnergy->at(v_idx_pfjet_PtEtaCut_ID[0]) );
	fillVariableWithValue( "PFJet1_Eta", HLTPFJetEta->at(v_idx_pfjet_PtEtaCut_ID[0]) );
	fillVariableWithValue( "PFJet1_Phi", HLTPFJetPhi->at(v_idx_pfjet_PtEtaCut_ID[0]) );
      }

    // 2nd jet 
    if( v_idx_pfjet_PtEtaCut_ID.size() >= 2 )        
      {
	fillVariableWithValue( "PFJet2_Pt", HLTPFJetPt->at(v_idx_pfjet_PtEtaCut_ID[1]) );
	fillVariableWithValue( "PFJet2_Energy", HLTPFJetEnergy->at(v_idx_pfjet_PtEtaCut_ID[1]) );
	fillVariableWithValue( "PFJet2_Eta", HLTPFJetEta->at(v_idx_pfjet_PtEtaCut_ID[1]) );
	fillVariableWithValue( "PFJet2_Phi", HLTPFJetPhi->at(v_idx_pfjet_PtEtaCut_ID[1]) );
      }

    // define booleans
    bool TwoPFJets=false;
    if( v_idx_pfjet_PtEtaCut_ID.size() >= 2 ) TwoPFJets = true;

    // Mjj
    if (TwoPFJets)
      {
	TLorentzVector jet1, jet2, jj;
	jet1.SetPtEtaPhiE(HLTPFJetPt->at(v_idx_pfjet_PtEtaCut_ID[0]),
			  HLTPFJetEta->at(v_idx_pfjet_PtEtaCut_ID[0]),
			  HLTPFJetPhi->at(v_idx_pfjet_PtEtaCut_ID[0]),
			  HLTPFJetEnergy->at(v_idx_pfjet_PtEtaCut_ID[0]) );
	jet2.SetPtEtaPhiE(HLTPFJetPt->at(v_idx_pfjet_PtEtaCut_ID[1]),
			  HLTPFJetEta->at(v_idx_pfjet_PtEtaCut_ID[1]),
			  HLTPFJetPhi->at(v_idx_pfjet_PtEtaCut_ID[1]),
			  HLTPFJetEnergy->at(v_idx_pfjet_PtEtaCut_ID[1]) );
	jj = jet1+jet2;
	
	fillVariableWithValue("M_PFJet1PFJet2", jj.M());
	fillVariableWithValue("Pt_PFJet1PFJet2", jj.Pt());
	fillVariableWithValue("DR_PFJet1PFJet2", jet1.DeltaR(jet2));
	fillVariableWithValue("DEta_PFJet1PFJet2", fabs(jet1.Eta()-jet2.Eta()) );
	fillVariableWithValue("DPhi_PFJet1PFJet2", jet1.DeltaPhi(jet2));	
      }
    
    // Evaluate cuts (but do not apply them)
    evaluateCuts();
    
    // Fill histograms and do analysis based on cut evaluation
    
    //INFO
    //      // retrieve value of previously filled variables (after making sure that they were filled)
    //      double totpTEle;
    //      if ( variableIsFilled("pT1stEle") && variableIsFilled("pT2ndEle") )
    //        totpTEle = getVariableValue("pT1stEle")+getVariableValue("pT2ndEle");
    //      // reject events that did not pass level 0 cuts
    //      if( !passedCut("0") ) continue;
    //      // ......

    //FillUserTH1D("h1_MTenu_PAS_minus", getVariableValue("MTenu_PAS"));
    //FillUserTH2D("h2_phi_VS_eta_1stEle", getVariableValue("Eta1stEle_PAS"), getVariableValue("Phi1stEle_PAS") );
    // 	    CreateAndFillUserTH1D("h1_ElectronRelIso_highMT", 1000, 0, 1, ElectronRelIso->at(myEle) );
    //if( passedAllPreviousCuts("d1_DPhi_METe_METj")
    //    && passedCut("nMuon_PtCut_IDISO")
    
    // Produce skim 
    if( //passedAllPreviousCuts("PassPrimaryVertex") 
 	passedCut("nPFJet")
	&& passedCut("PFJet1_Pt")
	&& passedCut("PFJet2_Pt")
	) 
      fillSkimTree();
    
    // Produce reduced skim
    if( //passedAllPreviousCuts("PassPrimaryVertex") 
 	passedCut("nPFJet")
	&& passedCut("PFJet1_Pt")
	&& passedCut("PFJet2_Pt")
	) 
      fillReducedSkimTree();
    
    ////////////////////// User's code to be done for every event - END ///////////////////////

  } // End of loop over events


  ////////////////////// User's code to write histos - BEGIN ///////////////////////

  ////////////////////// User's code to write histos - END ///////////////////////

  //STDOUT("analysisClass::Loop() ends");
}
