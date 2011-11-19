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

//For JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

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

  //--------------------------------------------------------------------------
  //For JEC
  //--------------------------------------------------------------------------

  string L2Tag     = "/afs/cern.ch/user/s/santanas/public/JEC_txtfiles_GR_R_42_V19/GR_R_42_V19_AK5PF_L2Relative.txt"; 
  string L3Tag     = "/afs/cern.ch/user/s/santanas/public/JEC_txtfiles_GR_R_42_V19/GR_R_42_V19_AK5PF_L3Absolute.txt";
  string L23ResTag = "/afs/cern.ch/user/s/santanas/public/JEC_txtfiles_GR_R_42_V19/GR_R_42_V19_AK5PF_L2L3Residual.txt";
  
  JetCorrectorParameters *L2Par = new JetCorrectorParameters(L2Tag);
  JetCorrectorParameters *L3Par = new JetCorrectorParameters(L3Tag);
  JetCorrectorParameters *L23ResPar = new JetCorrectorParameters(L23ResTag);

  vector<JetCorrectorParameters> vParL2;
  vParL2.push_back(*L2Par);
  vector<JetCorrectorParameters> vParL3;
  vParL3.push_back(*L3Par);
  vector<JetCorrectorParameters> vParL23Res;
  vParL23Res.push_back(*L23ResPar);
  vector<JetCorrectorParameters> vParAll;
  vParAll.push_back(*L2Par);
  vParAll.push_back(*L3Par);
  vParAll.push_back(*L23ResPar);

  FactorizedJetCorrector *JetCorrectorL2 = new FactorizedJetCorrector(vParL2);
  FactorizedJetCorrector *JetCorrectorL3 = new FactorizedJetCorrector(vParL3);
  FactorizedJetCorrector *JetCorrectorL23Res = new FactorizedJetCorrector(vParL23Res);
  FactorizedJetCorrector *JetCorrectorAll = new FactorizedJetCorrector(vParAll);

  //   JetCorrectorL2->setJetEta(0.5);
  //   JetCorrectorL3->setJetEta(0.5);
  //   JetCorrectorL23Res->setJetEta(0.5);
  //   JetCorrectorAll->setJetEta(0.5);
  //   JetCorrectorL2->setJetPt(100);
  //   JetCorrectorL3->setJetPt(100);
  //   JetCorrectorL23Res->setJetPt(100);
  //   JetCorrectorAll->setJetPt(100);

  //   double corrL2 = JetCorrectorL2->getCorrection();
  //   double corrL3 = JetCorrectorL3->getCorrection();
  //   double corrL23Res = JetCorrectorL23Res->getCorrection();
  //   double corrAll = JetCorrectorAll->getCorrection();
  //   cout << "L2, L3, L23Res, Product :  " << corrL2 << ", " <<corrL3 << ", " << corrL23Res << ", " << corrL2*corrL3*corrL23Res << endl;
  //   cout << "All : " << corrAll << endl;
  
  //--------------------------------------------------------------------------

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
  double pfjetPtCutForHT = getPreCutValue2("pfjetPtCut");
  double pfjetEtaCutForHT = getPreCutValue2("pfjetEtaCut");
  double pfjetDEtaForFatJet = getPreCutValue1("pfjetDEtaForFatJet");
  double pfjetRForFatJet = getPreCutValue1("pfjetRForFatJet");

  double caloRawjetPtCut = getPreCutValue1("caloRawjetPtCut");
  double caloRawjetEtaCut = getPreCutValue1("caloRawjetEtaCut");
  double caloRawjetIDloose = getPreCutValue1("caloRawjetIDloose");
  double caloRawjetIDtight = getPreCutValue1("caloRawjetIDtight");
  double caloRawjetPtCutForHT = getPreCutValue2("caloRawjetPtCut");
  double caloRawjetEtaCutForHT = getPreCutValue2("caloRawjetEtaCut");

  double caloCorrjetPtCut = getPreCutValue1("caloCorrjetPtCut");
  double caloCorrjetEtaCut = getPreCutValue1("caloCorrjetEtaCut");
  double caloCorrjetIDloose = getPreCutValue1("caloCorrjetIDloose");
  double caloCorrjetIDtight = getPreCutValue1("caloCorrjetIDtight");
  double caloCorrjetPtCutForHT = getPreCutValue2("caloCorrjetPtCut");
  double caloCorrjetEtaCutForHT = getPreCutValue2("caloCorrjetEtaCut");
  double caloCorrjetDEtaForFatJet = getPreCutValue1("caloCorrjetDEtaForFatJet");
  double caloCorrjetRForFatJet = getPreCutValue1("caloCorrjetRForFatJet");

  // Vertices
  double vertexMaxAbsZ = getPreCutValue1("vertexMaxAbsZ");
  double vertexMaxd0 = getPreCutValue1("vertexMaxd0");  

  //Others
  double applyPFJEC = getPreCutValue1("applyPFJEC");  
  
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

  CreateUserTH1D( "CaloCorrJetPrecutsPt"              ,    200,  0, 2000   );
  CreateUserTH1D( "CaloCorrJetPrecutsEta"             ,    100, -6, 6      );
  CreateUserTH1D( "CaloCorrJetPrecutsEnergyFractionEm"         , 200,  -2, 2   );
  CreateUserTH1D( "CaloCorrJetPrecutsEnergyFractionHadronic"   , 200,  -2, 2   );
  CreateUserTH1D( "CaloCorrJetIDPt"                   ,    200,  0, 2000   );
  CreateUserTH1D( "CaloCorrJetIDEta"                  ,    100, -6, 6      );

  ////////////////////// User's code to book histos - END ///////////////////////

  Long64_t nentries = fChain->GetEntriesFast();
  STDOUT("analysisClass::Loop(): nentries = " << nentries);

  ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
  ////// If the root version is updated and rootNtupleClass regenerated,     /////
  ////// these lines may need to be updated.                                 /////
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) { // Begin of loop over events
  //for (Long64_t jentry=0; jentry<2000;jentry++) { // Begin of loop over events
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

    //## Apply JEC to PF jets    
    if(applyPFJEC)
      {
	for(int ijet=0; ijet<HLTPFJetPt->size(); ijet++)
	  {
	    if(isData)  //Data: L2,L3,L23Res
	      {
		JetCorrectorAll->setJetEta( HLTPFJetEta->at(ijet) );
		JetCorrectorAll->setJetPt( HLTPFJetPt->at(ijet) );
		double corrAll = JetCorrectorAll->getCorrection();
		
		HLTPFJetPt->at(ijet) *= corrAll;
		HLTPFJetEnergy->at(ijet) *= corrAll;
	      }
	    else        //MC: L2,L3 only (no residual)
	      {
		JetCorrectorL2->setJetEta( HLTPFJetEta->at(ijet) );
		JetCorrectorL2->setJetPt( HLTPFJetPt->at(ijet) );
		JetCorrectorL3->setJetEta( HLTPFJetEta->at(ijet) );
		JetCorrectorL3->setJetPt( HLTPFJetPt->at(ijet) );
		double corrAll = JetCorrectorL2->getCorrection() * JetCorrectorL3->getCorrection();
		
		HLTPFJetPt->at(ijet) *= corrAll;
		HLTPFJetEnergy->at(ijet) *= corrAll;
	      }
	  }
      }

    //-----------------------------------------------------------------    
    // Get trigger information, if necessary                     
    //-----------------------------------------------------------------   

    //    if ( isData ) {
    vector<int>    *fakeHLTTriggerPrescales;
    fakeHLTTriggerPrescales = new std::vector<int> ( int (HLTTriggerNames -> size()), 1 );
    for (int ihlt=0 ; ihlt< HLTTriggerNames->size() ; ihlt++)
      {
	fakeHLTTriggerPrescales->push_back( 1 );
      }
    getTriggers ( HLTKey, HLTTriggerNames, HLTTriggerDecisions, fakeHLTTriggerPrescales ) ;
    delete fakeHLTTriggerPrescales;
    //    }
    
    //printTriggers();
    
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
    double HT_PFJets = 0.; 
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

	if( passjetID == true )
	  {
	    v_idx_pfjet_PtEtaCut_ID.push_back(v_idx_pfjet_PtEtaCut[ijet]);

	    //Histograms
	    FillUserTH1D( "PFJetIDPt"                  ,    HLTPFJetPt->at(v_idx_pfjet_PtEtaCut[ijet]) ); 
	    FillUserTH1D( "PFJetIDEta"                 ,    HLTPFJetEta->at(v_idx_pfjet_PtEtaCut[ijet]) ); 

	    //HT calculation
	    if( HLTPFJetPt->at(v_idx_pfjet_PtEtaCut[ijet]) > pfjetPtCutForHT  
		&& fabs( HLTPFJetEta->at(v_idx_pfjet_PtEtaCut[ijet]) ) < pfjetEtaCutForHT )
	      {
		HT_PFJets += HLTPFJetPt->at(v_idx_pfjet_PtEtaCut[ijet]);		
	      }
	  }

      } // End loop over jets


    //## Calo Corr Jets
    vector<int> v_idx_caloCorrjet_PtEtaCut;
    vector<int> v_idx_caloCorrjet_PtEtaCut_ID;

    // Precuts
    for(int ijet=0; ijet<HLTCaloJetCorrPt->size(); ijet++)
      {
	//pT/eta pre-cuts on jets
	if ( HLTCaloJetCorrPt->at(ijet) < caloCorrjetPtCut ) continue;
	if ( fabs( HLTCaloJetCorrEta->at(ijet) ) > caloCorrjetEtaCut ) continue;
	v_idx_caloCorrjet_PtEtaCut.push_back(ijet);
	
	//histograms
	FillUserTH1D( "CaloCorrJetPrecutsPt"                                      ,    HLTCaloJetCorrPt->at(ijet)  ); 
	FillUserTH1D( "CaloCorrJetPrecutsEta"                                     ,    HLTCaloJetCorrEta->at(ijet) ); 	
	FillUserTH1D( "CaloCorrJetPrecutsEnergyFractionEm"                        ,    HLTCaloJetCorrEnergyFractionEm->at(ijet) );
	FillUserTH1D( "CaloCorrJetPrecutsEnergyFractionHadronic"                  ,    HLTCaloJetCorrEnergyFractionHadronic->at(ijet) );
      }

    // Full selection
    double HT_CaloCorrJets = 0.; 
    for(int ijet=0; ijet<v_idx_caloCorrjet_PtEtaCut.size(); ijet++) 
      {
	bool passjetID = true;
	
	if( caloCorrjetIDloose && !caloCorrjetIDtight )
	  {
	    passjetID = true;
	    //FIXME
	  }

	if( caloCorrjetIDtight )
	  {
	    passjetID = true;
	    //FIXME
	  }

	if( passjetID == true )
	  {
	    v_idx_caloCorrjet_PtEtaCut_ID.push_back(v_idx_caloCorrjet_PtEtaCut[ijet]);

	    //Histograms
	    FillUserTH1D( "CaloCorrJetIDPt"                  ,    HLTCaloJetCorrPt->at(v_idx_caloCorrjet_PtEtaCut[ijet]) ); 
	    FillUserTH1D( "CaloCorrJetIDEta"                 ,    HLTCaloJetCorrEta->at(v_idx_caloCorrjet_PtEtaCut[ijet]) ); 

	    //HT calculation
	    if( HLTCaloJetCorrPt->at(v_idx_caloCorrjet_PtEtaCut[ijet]) > caloCorrjetPtCutForHT  
		&& fabs( HLTCaloJetCorrEta->at(v_idx_caloCorrjet_PtEtaCut[ijet]) ) < caloCorrjetEtaCutForHT )
	      {
		HT_CaloCorrJets += HLTCaloJetCorrPt->at(v_idx_caloCorrjet_PtEtaCut[ijet]);		
	      }
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


    // Set the value of the variableNames listed in the cutFile to their current value

    // Event filters
    fillVariableWithValue( "PassJSON", passJSON(run, ls, isData) );    

    int isPrimaryVertex = 0;
    if ( v_idx_vertex_good.size() > 0 )
      isPrimaryVertex = 1;
    fillVariableWithValue( "PassPrimaryVertex", isPrimaryVertex ) ;

    // Event info
    fillVariableWithValue( "isData", isData ) ;
    fillVariableWithValue( "bunch", bunch ) ;
    fillVariableWithValue( "event", event ) ;
    fillVariableWithValue( "ls", ls ) ;
    fillVariableWithValue( "orbit", orbit ) ;
    fillVariableWithValue( "run", run ) ;
    
    // Trigger (HLT)
    //     if(isData)
    //       {
    fillVariableWithValue("HLT_FatJetMass300" , triggerFired("DST_FatJetMass300_DR1p1_Deta2p0_v1") );
    fillVariableWithValue("HLT_FatJetMass400" , triggerFired("DST_FatJetMass400_DR1p1_Deta2p0_RunPF_v1") );
    fillVariableWithValue("HLT_HT150" , triggerFired("HLT_HT150_v11") );
    fillVariableWithValue("HLT_HT200" , triggerFired("HLT_HT200_v11") );
    fillVariableWithValue("HLT_HT250" , triggerFired("HLT_HT250_v11") );
    fillVariableWithValue("HLT_HT350" , triggerFired("DST_HT350_RunPF_v1") );
    //       }
    //     else
    //       {
    // 	fillVariableWithValue("HLT_FatJetMass300" , 1 );
    // 	fillVariableWithValue("HLT_FatJetMass400" , 1 );
    // 	fillVariableWithValue("HLT_HT150" , 1 );
    // 	fillVariableWithValue("HLT_HT200" , 1 );
    // 	fillVariableWithValue("HLT_HT250" , 1 );
    // 	fillVariableWithValue("HLT_HT350" , 1 );
    //       }

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
    fillVariableWithValue( "HT_PFJets", HT_PFJets ) ;

    fillVariableWithValue( "nCaloCJet", v_idx_caloCorrjet_PtEtaCut_ID.size() ) ;
    fillVariableWithValue( "HT_CaloCJets", HT_CaloCorrJets ) ;

    // 1st jet 
    if( v_idx_pfjet_PtEtaCut_ID.size() >= 1 )        
      {
	fillVariableWithValue( "PFJet1_Pt", HLTPFJetPt->at(v_idx_pfjet_PtEtaCut_ID[0]) );
	fillVariableWithValue( "PFJet1_Energy", HLTPFJetEnergy->at(v_idx_pfjet_PtEtaCut_ID[0]) );
	fillVariableWithValue( "PFJet1_Eta", HLTPFJetEta->at(v_idx_pfjet_PtEtaCut_ID[0]) );
	fillVariableWithValue( "PFJet1_Phi", HLTPFJetPhi->at(v_idx_pfjet_PtEtaCut_ID[0]) );
      }

    if( v_idx_caloCorrjet_PtEtaCut_ID.size() >= 1 )        
      {
	fillVariableWithValue( "CaloCJet1_Pt", HLTCaloJetCorrPt->at(v_idx_caloCorrjet_PtEtaCut_ID[0]) );
	fillVariableWithValue( "CaloCJet1_Energy", HLTCaloJetCorrEnergy->at(v_idx_caloCorrjet_PtEtaCut_ID[0]) );
	fillVariableWithValue( "CaloCJet1_Eta", HLTCaloJetCorrEta->at(v_idx_caloCorrjet_PtEtaCut_ID[0]) );
	fillVariableWithValue( "CaloCJet1_Phi", HLTCaloJetCorrPhi->at(v_idx_caloCorrjet_PtEtaCut_ID[0]) );
      }

    // 2nd jet 
    if( v_idx_pfjet_PtEtaCut_ID.size() >= 2 )        
      {
	fillVariableWithValue( "PFJet2_Pt", HLTPFJetPt->at(v_idx_pfjet_PtEtaCut_ID[1]) );
	fillVariableWithValue( "PFJet2_Energy", HLTPFJetEnergy->at(v_idx_pfjet_PtEtaCut_ID[1]) );
	fillVariableWithValue( "PFJet2_Eta", HLTPFJetEta->at(v_idx_pfjet_PtEtaCut_ID[1]) );
	fillVariableWithValue( "PFJet2_Phi", HLTPFJetPhi->at(v_idx_pfjet_PtEtaCut_ID[1]) );
      }

    if( v_idx_caloCorrjet_PtEtaCut_ID.size() >= 2 )        
      {
	fillVariableWithValue( "CaloCJet2_Pt", HLTCaloJetCorrPt->at(v_idx_caloCorrjet_PtEtaCut_ID[1]) );
	fillVariableWithValue( "CaloCJet2_Energy", HLTCaloJetCorrEnergy->at(v_idx_caloCorrjet_PtEtaCut_ID[1]) );
	fillVariableWithValue( "CaloCJet2_Eta", HLTCaloJetCorrEta->at(v_idx_caloCorrjet_PtEtaCut_ID[1]) );
	fillVariableWithValue( "CaloCJet2_Phi", HLTCaloJetCorrPhi->at(v_idx_caloCorrjet_PtEtaCut_ID[1]) );
      }

    // define booleans
    bool TwoPFJets=false;
    if( v_idx_pfjet_PtEtaCut_ID.size() >= 2 ) TwoPFJets = true;
    bool TwoCaloCorrJets=false;
    if( v_idx_caloCorrjet_PtEtaCut_ID.size() >= 2 ) TwoCaloCorrJets = true;


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

	// FatJetMass Algorithm (taken from HLTrigger/JetMET/src/HLTFatJetMassFilter.cc)
	TLorentzVector fatjet1, fatjet2, fatjj;
	if( fabs(jet1.Eta()-jet2.Eta()) < pfjetDEtaForFatJet )
	  {
	    for(int ijet=0; ijet<v_idx_pfjet_PtEtaCut_ID.size(); ijet++) 
	      {
		TLorentzVector currentjet;
		currentjet.SetPtEtaPhiE(HLTPFJetPt->at(v_idx_pfjet_PtEtaCut_ID[ijet]),
					HLTPFJetEta->at(v_idx_pfjet_PtEtaCut_ID[ijet]),
					HLTPFJetPhi->at(v_idx_pfjet_PtEtaCut_ID[ijet]),
					HLTPFJetEnergy->at(v_idx_pfjet_PtEtaCut_ID[ijet]) );
		
		double DeltaR1 = currentjet.DeltaR(jet1);
		double DeltaR2 = currentjet.DeltaR(jet2);
		if(DeltaR1 < DeltaR2 && DeltaR1 < pfjetRForFatJet) {
		  fatjet1 += currentjet;
		} else if(DeltaR2 < pfjetRForFatJet) {
		  fatjet2 += currentjet;
		}			
	      } //end creation of 2 fat jets
	    
	    fatjj = fatjet1 + fatjet2;
	    
	    fillVariableWithValue( "FatPFJet1_Pt", fatjet1.Pt() );
	    fillVariableWithValue( "FatPFJet1_Energy", fatjet1.Energy() );
	    fillVariableWithValue( "FatPFJet1_Eta", fatjet1.Eta() );
	    fillVariableWithValue( "FatPFJet1_Phi", fatjet1.Phi() );
	    
	    fillVariableWithValue( "FatPFJet2_Pt", fatjet2.Pt() );
	    fillVariableWithValue( "FatPFJet2_Energy", fatjet2.Energy() );
	    fillVariableWithValue( "FatPFJet2_Eta", fatjet2.Eta() );
	    fillVariableWithValue( "FatPFJet2_Phi", fatjet2.Phi() );
	    
	    fillVariableWithValue("M_FatPFJet1FatPFJet2", fatjj.M());
	    fillVariableWithValue("Pt_FatPFJet1FatPFJet2", fatjj.Pt());
	    fillVariableWithValue("DR_FatPFJet1FatPFJet2", fatjet1.DeltaR(fatjet2));
	    fillVariableWithValue("DEta_FatPFJet1FatPFJet2", fabs(fatjet1.Eta()-fatjet2.Eta()) );
	    fillVariableWithValue("DPhi_FatPFJet1FatPFJet2", fatjet1.DeltaPhi(fatjet2));		   
	  }
      }

    if (TwoCaloCorrJets)
      {
	TLorentzVector jet1, jet2, jj;
	jet1.SetPtEtaPhiE(HLTCaloJetCorrPt->at(v_idx_caloCorrjet_PtEtaCut_ID[0]),
			  HLTCaloJetCorrEta->at(v_idx_caloCorrjet_PtEtaCut_ID[0]),
			  HLTCaloJetCorrPhi->at(v_idx_caloCorrjet_PtEtaCut_ID[0]),
			  HLTCaloJetCorrEnergy->at(v_idx_caloCorrjet_PtEtaCut_ID[0]) );
	jet2.SetPtEtaPhiE(HLTCaloJetCorrPt->at(v_idx_caloCorrjet_PtEtaCut_ID[1]),
			  HLTCaloJetCorrEta->at(v_idx_caloCorrjet_PtEtaCut_ID[1]),
			  HLTCaloJetCorrPhi->at(v_idx_caloCorrjet_PtEtaCut_ID[1]),
			  HLTCaloJetCorrEnergy->at(v_idx_caloCorrjet_PtEtaCut_ID[1]) );
	jj = jet1+jet2;
	
	fillVariableWithValue("M_CaloCJet1CaloCJet2", jj.M());
	fillVariableWithValue("Pt_CaloCJet1CaloCJet2", jj.Pt());
	fillVariableWithValue("DR_CaloCJet1CaloCJet2", jet1.DeltaR(jet2));
	fillVariableWithValue("DEta_CaloCJet1CaloCJet2", fabs(jet1.Eta()-jet2.Eta()) );
	fillVariableWithValue("DPhi_CaloCJet1CaloCJet2", jet1.DeltaPhi(jet2));	

	// FatJetMass Algorithm (taken from HLTrigger/JetMET/src/HLTFatJetMassFilter.cc)
	TLorentzVector fatjet1, fatjet2, fatjj;
	if( fabs(jet1.Eta()-jet2.Eta()) < caloCorrjetDEtaForFatJet 
	    && jet1.Pt()>30 && jet2.Pt()>30 && fabs(jet1.Eta())<3 && fabs(jet2.Eta())<3 // doublejetcentral30 requirement
	    )  
	  {
	    for(int ijet=0; ijet<v_idx_caloCorrjet_PtEtaCut_ID.size(); ijet++) 
	      {
		TLorentzVector currentjet;
		currentjet.SetPtEtaPhiE(HLTCaloJetCorrPt->at(v_idx_caloCorrjet_PtEtaCut_ID[ijet]),
					HLTCaloJetCorrEta->at(v_idx_caloCorrjet_PtEtaCut_ID[ijet]),
					HLTCaloJetCorrPhi->at(v_idx_caloCorrjet_PtEtaCut_ID[ijet]),
					HLTCaloJetCorrEnergy->at(v_idx_caloCorrjet_PtEtaCut_ID[ijet]) );
		
		double DeltaR1 = currentjet.DeltaR(jet1);
		double DeltaR2 = currentjet.DeltaR(jet2);
		if(DeltaR1 < DeltaR2 && DeltaR1 < caloCorrjetRForFatJet) {
		  fatjet1 += currentjet;
		} else if(DeltaR2 < caloCorrjetRForFatJet) {
		  fatjet2 += currentjet;
		}			
	      } //end creation of 2 fat jets
	    
	    fatjj = fatjet1 + fatjet2;
	    
	    fillVariableWithValue( "FatCaloCJet1_Pt", fatjet1.Pt() );
	    fillVariableWithValue( "FatCaloCJet1_Energy", fatjet1.Energy() );
	    fillVariableWithValue( "FatCaloCJet1_Eta", fatjet1.Eta() );
	    fillVariableWithValue( "FatCaloCJet1_Phi", fatjet1.Phi() );
	    
	    fillVariableWithValue( "FatCaloCJet2_Pt", fatjet2.Pt() );
	    fillVariableWithValue( "FatCaloCJet2_Energy", fatjet2.Energy() );
	    fillVariableWithValue( "FatCaloCJet2_Eta", fatjet2.Eta() );
	    fillVariableWithValue( "FatCaloCJet2_Phi", fatjet2.Phi() );
	    
	    fillVariableWithValue("M_FatCaloCJet1FatCaloCJet2", fatjj.M());
	    fillVariableWithValue("Pt_FatCaloCJet1FatCaloCJet2", fatjj.Pt());
	    fillVariableWithValue("DR_FatCaloCJet1FatCaloCJet2", fatjet1.DeltaR(fatjet2));
	    fillVariableWithValue("DEta_FatCaloCJet1FatCaloCJet2", fabs(fatjet1.Eta()-fatjet2.Eta()) );
	    fillVariableWithValue("DPhi_FatCaloCJet1FatCaloCJet2", fatjet1.DeltaPhi(fatjet2));		   
	  }
      }
    
    // Evaluate cuts (but do not apply them)
    evaluateCuts();
    
    // Fill histograms and do analysis based on cut evaluation
    
    //     if( variableIsFilled("HLT_HT350") &&  variableIsFilled("HLT_FatJetMass400") && variableIsFilled("PFJet1_Pt")    
    // 	&& max( getVariableValue("HLT_HT350"), getVariableValue("HLT_FatJetMass400") ) == 0 
    // 	&& getVariableValue("PFJet1_Pt")>0 )
    //       {
    // 	std::cout << "************** next event" << std::endl;
    // 	printTriggers();
    // 	std::cout << getVariableValue("PFJet1_Pt") << std::endl;
    // 	std::cout << "**************" << std::endl;
    //       }
    
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
    //     if( //passedAllPreviousCuts("PassPrimaryVertex") 
    //  	passedCut("nPFJet")
    // 	&& passedCut("PFJet1_Pt")
    // 	&& passedCut("PFJet2_Pt")
    // 	) 
    fillSkimTree();
    
    // Produce reduced skim
    //     if( //passedAllPreviousCuts("PassPrimaryVertex") 
    //  	passedCut("nPFJet")
    // 	&& passedCut("PFJet1_Pt")
    // 	&& passedCut("PFJet2_Pt")
    // 	) 
    fillReducedSkimTree();
    
    ////////////////////// User's code to be done for every event - END ///////////////////////

  } // End of loop over events


  ////////////////////// User's code to write histos - BEGIN ///////////////////////

  ////////////////////// User's code to write histos - END ///////////////////////

  //STDOUT("analysisClass::Loop() ends");
}
