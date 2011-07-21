#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TMatrixDSym.h>
#include <TMath.h>

//////////////////////////////////////////////////////////////////
// Solve system of non-linear n equations with n variables - BEGIN

// GSL headers
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlin.h>

struct rparams
{
  //  double a0;
  double pxj1;
  double pyj1;
  double pzj1;
  double pxj2;
  double pyj2;
  double pzj2;
  double metx;
  double mety;
};

int rosenbrock_f (const gsl_vector * x, void *params, gsl_vector * f)
{
  double pxj1 = ((struct rparams *) params)->pxj1;
  double pyj1 = ((struct rparams *) params)->pyj1;
  double pzj1 = ((struct rparams *) params)->pzj1;
  double pxj2 = ((struct rparams *) params)->pxj2;
  double pyj2 = ((struct rparams *) params)->pyj2;
  double pzj2 = ((struct rparams *) params)->pzj2;
  double metx = ((struct rparams *) params)->metx;
  double mety = ((struct rparams *) params)->mety;
  
  const double pxlq1 = gsl_vector_get (x, 0);
  const double pylq1 = gsl_vector_get (x, 1);
  const double pzlq1 = gsl_vector_get (x, 2);
  const double mlq   = gsl_vector_get (x, 3);

//   const double b2 = bx*bx + by*by + bz*bz;
//   const double gm = 1 / sqrt(1-b2);  
//   TMatrixDSym bst(4,4);
//   bst(0,0)=gm;     bst(0,1)=-bx*gm;                bst(0,2)=-by*gm;              bst(0,3)=-bz*gm; 
//                    bst(1,1)=1+(gm-1)*(bx*bx/b2);   bst(1,2)=(gm-1)*(bx*by/b2);   bst(1,3)=(gm-1)*(bx*bz/b2); 
//                                      		   bst(2,2)=1+(gm-1)*(by*by/b2); bst(2,3)=(gm-1)*(by*bz/b2); 
//                                       						 bst(3,3)=1+(gm-1)*(bz*bz/b2); 

  TLorentzVector LQ1;
  TLorentzVector LQ2;
  LQ1.SetPxPyPzE( pxlq1, pylq1, pzlq1, sqrt(pxlq1*pxlq1+pylq1*pylq1+pzlq1*pzlq1+mlq*mlq));
  LQ2.SetPxPyPzE(-pxlq1,-pylq1,-pzlq1, sqrt(pxlq1*pxlq1+pylq1*pylq1+pzlq1*pzlq1+mlq*mlq));
  cout<<"lq1     Px, Py, Pz, M = "<<pxlq1<<"\t"<<pylq1<<"\t"<<pzlq1<<"\t"<<LQ1.M()<<endl;
  double bx = LQ1.Px()/LQ1.E();
  double by = LQ1.Py()/LQ1.E();
  double bz = LQ1.Pz()/LQ1.E();

  TLorentzVector j1,j2;
  j1.SetPxPyPzE(pxj1,pyj1,pzj1,sqrt(pxj1*pxj1+pyj1*pyj1+pzj1*pzj1));
  j2.SetPxPyPzE(pxj2,pyj2,pzj2,sqrt(pxj2*pxj2+pyj2*pyj2+pzj2*pzj2));
  TLorentzVector j1_SR1,j2_SR2;
  j1_SR1=j1; j1_SR1.Boost(-bx,-by,-bz);
  j2_SR2=j2; j2_SR2.Boost( bx, by, bz);
  TLorentzVector nu1_SR1,nu2_SR2;
  nu1_SR1=-j1_SR1; nu1_SR1.SetE(j1_SR1.E());
  nu2_SR2=-j2_SR2; nu2_SR2.SetE(j2_SR2.E());
  TLorentzVector LQ1_SR1;   LQ1_SR1=j1_SR1+nu1_SR1;
  TLorentzVector LQ2_SR2;   LQ2_SR2=j2_SR2+nu2_SR2;
  TLorentzVector nu1,nu2;
  nu1=nu1_SR1; nu1.Boost( bx, by, bz);
  nu2=nu2_SR2; nu2.Boost(-bx,-by,-bz);
  cout<<"---------------------"<<endl;
  cout<<"bx by bz and b = "<<bx<<"\t"<<by<<"\t"<<bz<<"\t"<<sqrt(bx*bx+by*by+bz*bz)<<endl;
  cout<<"j1     Px, Py, Pz, P = "<<j1.Px()<<"\t"<<j1.Py()<<"\t"<<j1.Pz()<<"\t"<<j1.P()<<endl;
  cout<<"j2     Px, Py, Pz, P = "<<j2.Px()<<"\t"<<j2.Py()<<"\t"<<j2.Pz()<<"\t"<<j2.P()<<endl;
  cout<<"nu1     Px, Py, Pz, P = "<<nu1.Px()<<"\t"<<nu1.Py()<<"\t"<<nu1.Pz()<<"\t"<<nu1.P()<<endl;
  cout<<"nu2     Px, Py, Pz, P = "<<nu2.Px()<<"\t"<<nu2.Py()<<"\t"<<nu2.Pz()<<"\t"<<nu2.P()<<endl;
  cout<<"nu1+nu2 Px Py = "<<nu1.Px()+nu2.Px()<<"\t"<<nu1.Py()+nu2.Py()<<endl;
  cout<<"j1_SR1 Px, Py, Pz, P = "<<j1_SR1.Px()<<"\t"<<j1_SR1.Py()<<"\t"<<j1_SR1.Pz()<<"\t"<<j1_SR1.P()<<endl;
  cout<<"j2_SR2 Px, Py, Pz, P = "<<j2_SR2.Px()<<"\t"<<j2_SR2.Py()<<"\t"<<j2_SR2.Pz()<<"\t"<<j2_SR2.P()<<endl;
  cout<<"nu1_SR1 Px, Py, Pz, P = "<<nu1_SR1.Px()<<"\t"<<nu1_SR1.Py()<<"\t"<<nu1_SR1.Pz()<<"\t"<<nu1_SR1.P()<<endl;
  cout<<"nu2_SR2 Px, Py, Pz, P = "<<nu2_SR2.Px()<<"\t"<<nu2_SR2.Py()<<"\t"<<nu2_SR2.Pz()<<"\t"<<nu2_SR2.P()<<endl;
  cout <<"LQ1.M LQ2.M = "<<LQ1.M()<<"\t"<<LQ2.M()<<endl;
  cout <<"LQ1_SR1.M LQ2_SR2.M = "<<LQ1_SR1.M()<<"\t"<<LQ2_SR2.M()<<endl;
  cout<<"---------------------"<<endl;

  const double y0 = LQ1_SR1.M() - fabs(mlq);
  const double y1 = LQ2_SR2.M() - fabs(mlq);
  const double y2 = nu1.Px()+nu2.Px()-metx;
  const double y3 = nu1.Py()+nu2.Py()-mety;

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);
  gsl_vector_set (f, 2, y2);
  gsl_vector_set (f, 3, y3);
  
  return GSL_SUCCESS;
}

void print_state (size_t iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3u x = %.3f %.3f %.3f %.3f "
	          "f(x) = %.3e %.3e %.3e %.3e \n",
	  (uint) iter,
	  gsl_vector_get (s->x, 0),
	  gsl_vector_get (s->x, 1),
	  gsl_vector_get (s->x, 2),
	  gsl_vector_get (s->x, 3),
	  gsl_vector_get (s->f, 0),
	  gsl_vector_get (s->f, 1),
	  gsl_vector_get (s->f, 2),
	  gsl_vector_get (s->f, 3)
	  );
}

// Solve system of non-linear n equations with n variables - END
////////////////////////////////////////////////////////////////



//-----------------------------
//### JetID ### --> see https://twiki.cern.ch/twiki/bin/view/CMS/ExoticaHighPtJets#JetId

bool JetIdloose(double ak5JetJIDresEMF,double ak5JetJIDfHPD,int ak5JetJIDn90Hits, double ak5JetEta){
  bool jetidloose=false;
  bool jetidresEMF=true;

  double fhpdmax = 0.98;
  double n90hitsmin =1;
  double emf_min = 0.01;

  if(fabs(ak5JetEta)<2.6 && ak5JetJIDresEMF<=emf_min) jetidresEMF=false;

  if(jetidresEMF && ak5JetJIDfHPD<fhpdmax && ak5JetJIDn90Hits>n90hitsmin) {
    jetidloose=true;
  }
  return jetidloose;
}

bool JetIdtight(double ak5JetJIDresEMF,double ak5JetJIDfHPD,int ak5JetJIDn90Hits, double ak5JetEta, double ak5JetPtRaw){
  bool jetidtight=false;
  bool jetidresEMF=true;
  bool jetidfHPD_highPt=true;

  double fhpdmax = 0.98;
  double n90hitsmin =1;
  double emf_min = 0.01;

  if(fabs(ak5JetEta)<2.6 && ak5JetJIDresEMF<=emf_min) jetidresEMF=false;
  if(fabs(ak5JetEta)<2.6 && ak5JetPtRaw>80 && ak5JetJIDresEMF>=1) jetidresEMF=false;
  if(ak5JetPtRaw>25 && ak5JetJIDfHPD>=0.95) jetidfHPD_highPt=false;

  if(jetidresEMF && jetidfHPD_highPt && ak5JetJIDfHPD<fhpdmax && ak5JetJIDn90Hits>n90hitsmin) 
    {
      jetidtight=true;
    }
  return jetidtight;
}

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
   
  ////////////////////// User's code to book histos - BEGIN ///////////////////////

  TH1F *h_Mej_PAS = new TH1F ("h_Mej_PAS","h_Mej_PAS",200,0,2000);  h_Mej_PAS->Sumw2();

  CreateUserTH1D("sT_fullSeleNoPreStCut", 200, 0, 2000);

  // tests/studies
  CreateUserTH1D("LQ1_Eta", 100, -5, 5);
  CreateUserTH1D("LQ2_Eta", 100, -5, 5);
  CreateUserTH1D("LQ1Eta_plus_LQ2Eta", 100, -7, 7);
  CreateUserTH2D("LQ1Eta_vs_LQ2Eta", 100, -7, 7, 100, -7, 7);
  CreateUserTH1D("LQ12Pz", 200, -1000, 1000);
  CreateUserTH1D("LQ12Beta", 50, 0, 1);
  CreateUserTH1D("LQ12Betaz", 100, -1, 1);
  CreateUserTH1D("LQ12Betaz_abs", 50, 0, 1);
  CreateUserTH1D("LQ12Betat", 50, 0, 1);
  CreateUserTH1D("LQ1Betaz", 100, -1, 1);
  CreateUserTH1D("LQ1Betax", 100, -1, 1);
  CreateUserTH1D("LQ1Betay", 100, -1, 1);
  CreateUserTH1D("LQ1Betat", 50, 0, 1);
  CreateUserTH1D("MR", 100, 0, 2000);
  CreateUserTH1D("MTR", 100, 0, 2000);
  CreateUserTH1D("Razor", 100, 0, 3);
  CreateUserTH1D("Razor_LQ12Betaz06", 100, 0, 3);
  CreateUserTH1D("Razor_LQ1Betat01", 100, 0, 3);

  ////////////////////// User's code to book histos - END ///////////////////////

  ////////////////////// User's code to get preCut values - BEGIN ///////////////

  double eleEta_bar = getPreCutValue1("eleEta_bar");
  double eleEta_end_min = getPreCutValue1("eleEta_end");
  double eleEta_end_max = getPreCutValue2("eleEta_end");

  double EleEnergyScale_EB=getPreCutValue1("EleEnergyScale_EB");
  double EleEnergyScale_EE=getPreCutValue1("EleEnergyScale_EE");
  double JetEnergyScale=getPreCutValue1("JetEnergyScale");

  // Not used when using ElectronHeepID and heepBitMask // int eleIDType = (int) getPreCutValue1("eleIDType");
  int heepBitMask_EB  =  getPreCutValue1("heepBitMask_EBGapEE") ;
  int heepBitMask_GAP =  getPreCutValue2("heepBitMask_EBGapEE") ;

  int heepBitMask_EE  =  getPreCutValue3("heepBitMask_EBGapEE") ;

  int looseBitMask_EB       =  getPreCutValue1("looseBitMask_EBGapEE") ;
  int looseBitMask_GAP      =  getPreCutValue2("looseBitMask_EBGapEE") ;
  int looseBitMask_EE       =  getPreCutValue3("looseBitMask_EBGapEE") ;
  int looseBitMask_enabled  =  getPreCutValue4("looseBitMask_EBGapEE") ;

  double muon_PtCut = getPreCutValue1("muon_PtCut");
  double muFidRegion = getPreCutValue1("muFidRegion"); // currently unused !!!
  double muNHits_minThresh = getPreCutValue1("muNHits_minThresh");
  double muTrkD0Maximum = getPreCutValue1("muTrkD0Maximum");

  int jetAlgorithm = getPreCutValue1("jetAlgorithm");

  double vertexMinimumNDOF = getPreCutValue1("vertexMinimumNDOF");
  double vertexMaxAbsZ = getPreCutValue1("vertexMaxAbsZ");
  double vertexMaxd0 = getPreCutValue1("vertexMaxd0");

  ////////////////////// User's code to get preCut values - END /////////////////
    
  Long64_t nentries = fChain->GetEntriesFast();
  STDOUT("analysisClass::Loop(): nentries = " << nentries);   
  
  ////// The following ~7 lines have been taken from rootNtupleClass->Loop() /////
  ////// If the root version is updated and rootNtupleClass regenerated,     /////
  ////// these lines may need to be updated.                                 /////    
  Long64_t nbytes = 0, nb = 0;
  //  for (Long64_t jentry=0; jentry<nentries;jentry++) { // Begin of loop over events
  for (Long64_t jentry=0; jentry<2000;jentry++) { // Begin of loop over events
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(jentry < 10 || jentry%1000 == 0) STDOUT("analysisClass::Loop(): jentry = " << jentry);   
    // if (Cut(ientry) < 0) continue;
    
    ////////////////////// User's code to be done for every event - BEGIN ///////////////////////

    //## Define new jet collection
    std::auto_ptr<std::vector<double> >  JetPt  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetPtRaw  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetEnergy  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetEta  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<double> >  JetPhi  ( new std::vector<double>()  );
    std::auto_ptr<std::vector<int> >     JetPassID  ( new std::vector<int>()  );
    std::auto_ptr<std::vector<double> >  JetTCHE  ( new std::vector<double>()  );

    if(jetAlgorithm==1) //PF jets
      {
	for (int ijet=0 ; ijet< PFJetPt->size() ; ijet++)
	  {
	    JetPt->push_back( PFJetPt->at(ijet) );
	    JetPtRaw->push_back( PFJetPtRaw->at(ijet) );
	    JetEnergy->push_back( PFJetEnergy->at(ijet) );
	    JetEta->push_back( PFJetEta->at(ijet) );
	    JetPhi->push_back( PFJetPhi->at(ijet) );
            JetPassID->push_back( PFJetPassLooseID->at(ijet) );
	    //     JetPassID->push_back(
	    //  PFJetIdloose(PFJetChargedHadronEnergyFraction->at(ijet),
	    //       PFJetChargedEmEnergyFraction->at(ijet),
	    //       PFJetNeutralHadronEnergyFraction->at(ijet),
	    //       PFJetNeutralEmEnergyFraction->at(ijet),
	    //       PFJetEta->at(ijet) )
	    //  );
            JetTCHE->push_back( PFJetTrackCountingHighEffBTag->at(ijet) );
	  }//end loop over pf jets
      }//end if "pf jets"

    if(jetAlgorithm==2) //Calo jets
      {
	for (int ijet=0 ; ijet < CaloJetPt->size() ; ijet++)
	  {
	    JetPt->push_back( CaloJetPt->at(ijet) );
	    JetPtRaw->push_back( CaloJetPtRaw->at(ijet) );
	    JetEnergy->push_back( CaloJetEnergy->at(ijet) );
	    JetEta->push_back( CaloJetEta->at(ijet) );
	    JetPhi->push_back( CaloJetPhi->at(ijet) );
            JetPassID->push_back( CaloJetPassLooseID->at(ijet) );
	    //     JetPassID->push_back(
	    //  JetIdloose(CaloJetresEMF->at(ijet),
	    //     CaloJetfHPD->at(ijet),
	    //     CaloJetn90Hits->at(ijet),
	    //     CaloJetEta->at(ijet) )
	    //  );
            JetTCHE->push_back( CaloJetTrackCountingHighEffBTag->at(ijet) );
	  }//end loop over calo jets
      }//end if "calo jets"


    // EES and JES
    if( EleEnergyScale_EB != 1 || EleEnergyScale_EE != 1 )
      {
	for(int iele=0; iele<ElectronPt->size(); iele++)
	  {
	    if( fabs(ElectronEta->at(iele)) < eleEta_bar )
	      ElectronPt->at(iele) *= EleEnergyScale_EB;
	    if( fabs(ElectronEta->at(iele)) > eleEta_end_min && fabs(ElectronEta->at(iele)) < eleEta_end_max )
	      ElectronPt->at(iele) *= EleEnergyScale_EE;
	  }
      }
    if( JetEnergyScale != 1 )
      {
	for(int ijet=0; ijet<JetPt->size(); ijet++)
	  {
	    JetPt->at(ijet) *= JetEnergyScale;
	  }
      }

    //## HLT
    int PassTrig = 0;
    int HLTFromRun[4] = {getPreCutValue1("HLTFromRun"),
			 getPreCutValue2("HLTFromRun"),
			 getPreCutValue3("HLTFromRun"),
			 getPreCutValue4("HLTFromRun")};
    int HLTTrigger[4] = {getPreCutValue1("HLTTrigger"),
			 getPreCutValue2("HLTTrigger"),
			 getPreCutValue3("HLTTrigger"),
			 getPreCutValue4("HLTTrigger")};
    int HLTTrgUsed;
    for (int i=0; i<4; i++) {
      if ( !isData && i != 0) continue; // For MC use HLTPhoton15 as the cleaned trigger is not in MC yet as of July 20, 2010
      if ( HLTFromRun[i] <= run ) {
 	//if(jentry == 0 ) STDOUT("run, i, HLTTrigger[i], HLTFromRun[i] = "<<run<<"\t"<<i<<"\t"<<"\t"<<HLTTrigger[i]<<"\t"<<HLTFromRun[i]);
	if (HLTTrigger[i] > 0 && HLTTrigger[i] < HLTResults->size() ) {
	  PassTrig=HLTResults->at(HLTTrigger[i]);
	  HLTTrgUsed=HLTTrigger[i];
	} else {
	  STDOUT("ERROR: HLTTrigger out of range of HLTResults: HLTTrigger = "<<HLTTrigger[i] <<"and HLTResults size = "<< HLTResults->size());
	}
      }
    }
    if(jentry == 0 ) STDOUT("Run = "<<run <<", HLTTrgUsed is number = "<<HLTTrgUsed<<" of the list HLTPathsOfInterest");

//     if ( isData==true)
//       {
// 	// Skip runs declared bad for HCAL reasons
// 	if ( run == 146511 ) continue;
// 	if ( run == 146513 ) continue;
// 	if ( run == 146514 ) continue;
// 	if ( run == 146644 ) continue;
// 	// Skip runs declared bad for other reasons
// 	if ( run == 141874 ) continue;
// 	if ( run == 141876 ) continue;
// 	if ( run == 142414 ) continue;
// 	if ( run == 147929 and ls == 619 ) continue;
//       }


    int LQ1=-1, LQ2=-1;
    for(int igen=0;igen<GenParticleEnergy->size();igen++) // Loop over gen particles
      {
	int LQ_PID=42;
        if( abs(GenParticlePdgId->at(igen))==LQ_PID ) 
	  {
	    if (LQ1 == -1) 
	      {
		LQ1=igen;
	      }
	    else if (LQ2 == -1)
	      {
		LQ2=igen;
	      }
	  }
      }
    //STDOUT("LQ1, LQ2= "<<LQ1<<", "<<LQ2);

    //GenParticlePt->at(igen)

    // Electrons
    vector<int> v_idx_ele_all;
    vector<int> v_idx_ele_PtCut;
    vector<int> v_idx_ele_PtCut_IDISO_noOverlap;
    int heepBitMask;

    //Loop over electrons
    for(int iele=0; iele<ElectronPt->size(); iele++)
      {

	// Reject ECAL spikes
	if ( 1 - ElectronSCS4S1->at(iele) > 0.95 ) continue; 

	//no cut on reco electrons
	v_idx_ele_all.push_back(iele); 

	//pT pre-cut on ele
	if( ElectronPt->at(iele) < getPreCutValue1("ele_PtCut") ) continue; 
	v_idx_ele_PtCut.push_back(iele);

	// get heepBitMask for EB, GAP, EE 
	if( fabs(ElectronEta->at(iele)) < eleEta_bar ) 
	  {
	    heepBitMask = heepBitMask_EB;
	  }
	else if ( fabs(ElectronEta->at(iele)) > eleEta_end_min && fabs(ElectronEta->at(iele)) < eleEta_end_max ) 
	  {
	    heepBitMask = heepBitMask_EE;
	  }
	else {
	    heepBitMask = heepBitMask_GAP;
	}

	//ID + ISO + NO overlap with good muons 
	// int eleID = ElectronPassID->at(iele);
	// if ( (eleID & 1<<eleIDType) > 0  && ElectronOverlaps->at(iele)==0 )
	if ( (ElectronHeepID->at(iele) & ~heepBitMask)==0x0  && ElectronOverlaps->at(iele)==0 )
	  {
	    //STDOUT("ElectronHeepID = " << hex << ElectronHeepID->at(iele) << " ; ElectronPassID = " << ElectronPassID->at(iele) )
	    v_idx_ele_PtCut_IDISO_noOverlap.push_back(iele);
	  }

      } // End loop over electrons

    // tight-loose electrons, if enabled from cut file
    if ( looseBitMask_enabled == 1 && v_idx_ele_PtCut_IDISO_noOverlap.size() == 1 )
      {
	//STDOUT("v_idx_ele_PtCut_IDISO_noOverlap[0] = "<<v_idx_ele_PtCut_IDISO_noOverlap[0] << " - Pt = "<<ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]));
	bool loosePtLargerThanTightPt = true;
	for(int iele=0; iele<v_idx_ele_PtCut.size(); iele++)
	  {
	    
	    if (v_idx_ele_PtCut[iele] == v_idx_ele_PtCut_IDISO_noOverlap[0])
	      {
		loosePtLargerThanTightPt = false;
		continue;
	      }
	    // get looseBitMask for EB, GAP, EE 

	    int looseBitMask;
	    if( fabs(ElectronEta->at(v_idx_ele_PtCut[iele])) < eleEta_bar ) 
	      {
		looseBitMask = looseBitMask_EB;
	      }
	    else if ( fabs(ElectronEta->at(v_idx_ele_PtCut[iele])) > eleEta_end_min && fabs(ElectronEta->at(v_idx_ele_PtCut[iele])) < eleEta_end_max ) 
	      {
		looseBitMask = looseBitMask_EE;
	      }
	    else {
	      looseBitMask = looseBitMask_GAP;
	    }
	    if ( (ElectronHeepID->at(v_idx_ele_PtCut[iele]) & ~looseBitMask)==0x0  && ElectronOverlaps->at(v_idx_ele_PtCut[iele])==0 )
	      {
		if ( loosePtLargerThanTightPt )
		  {
		    v_idx_ele_PtCut_IDISO_noOverlap.insert(v_idx_ele_PtCut_IDISO_noOverlap.begin(),1,v_idx_ele_PtCut[iele]);
		  }
		else 
		  {
		    v_idx_ele_PtCut_IDISO_noOverlap.push_back(v_idx_ele_PtCut[iele]);
		  }
		break; // happy with one loose electron - Note: if you want more than 1 loose, pt sorting will not be OK with the code as is. 
	      }
	  }	
// 	for ( int i=0; i<v_idx_ele_PtCut_IDISO_noOverlap.size(); i++)
// 	  {
// 	    STDOUT("i="<<i<<", v_idx_ele_PtCut_IDISO_noOverlap[i] = "<<v_idx_ele_PtCut_IDISO_noOverlap[i] << ", Pt = "<<ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[i]));
// 	  }
      } // tight-loose electrons, if enabled from cut file

    // Jets
    vector<int> v_idx_jet_all;
    vector<int> v_idx_jet_PtCut;
    vector<int> v_idx_jet_PtCut_noOverlap;
    vector<int> v_idx_jet_PtCut_noOverlap_ID;
    vector<int> v_idx_jet_PtCut_noOverlap_ID_EtaCut;

    // Loop over jets
    for(int ijet=0; ijet<JetPt->size(); ijet++)
      {
	//no cut on reco jets
	v_idx_jet_all.push_back(ijet);

	//pT pre-cut on reco jets
	if ( JetPt->at(ijet) < getPreCutValue1("jet_PtCut") ) continue;
	v_idx_jet_PtCut.push_back(ijet);
      }

// 	//Checking overlap between electrons and jets
// 	int JetOverlapsWithEle = 0; //don't change the default (0) 
// 	float minDeltaR=9999.;
// 	TVector3 jet_vec;
// 	jet_vec.SetPtEtaPhi(JetPt->at(ijet),JetEta->at(ijet),JetPhi->at(ijet));
// 	for (int i=0; i < v_idx_ele_PtCut_IDISO_noOverlap.size(); i++){
// 	  TVector3 ele_vec;	  
// 	  ele_vec.SetPtEtaPhi(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[i])
// 			      ,ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[i])
// 			      ,ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[i]));
// 	  double distance = jet_vec.DeltaR(ele_vec);
// 	  if (distance<minDeltaR) minDeltaR=distance;
// 	}
// 	if ( minDeltaR < getPreCutValue1("jet_ele_DeltaRcut") )
// 	  JetOverlapsWithEle = 1; //this jet overlaps with a good electron --> remove it from the analysis

// 	//pT pre-cut + no overlaps with electrons
// 	// ---- use the flag stored in rootTuples
// 	//if( ( JetOverlaps->at(ijet) & 1 << eleIDType) == 0)/* NO overlap with electrons */  
// 	// && (caloJetOverlaps[ijet] & 1 << 5)==0 )/* NO overlap with muons */   
// 	// ----
// 	if( JetOverlapsWithEle == 0 )  /* NO overlap with electrons */  
// 	  v_idx_jet_PtCut_noOverlap.push_back(ijet);

    vector <int> jetFlags(v_idx_jet_PtCut.size(), 0);
    int Njetflagged = 0;
    for (int iele=0; iele<v_idx_ele_PtCut_IDISO_noOverlap.size(); iele++)
      {
	TLorentzVector ele;
        ele.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
			 ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
			 ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),0);
	TLorentzVector jet;
	double minDR=9999.;
	int ijet_minDR = -1;    
        for(int ijet=0; ijet<v_idx_jet_PtCut.size(); ijet++)
          {
	    if ( jetFlags[ijet] == 1 ) 
	      continue;
            jet.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut[ijet]),
			     JetEta->at(v_idx_jet_PtCut[ijet]),
			     JetPhi->at(v_idx_jet_PtCut[ijet]),0);
	    double DR = jet.DeltaR(ele);
	    if (DR<minDR) 
	      {
		minDR = DR;
		ijet_minDR = ijet;
	      }
	  }
	if ( minDR < getPreCutValue1("jet_ele_DeltaRcut") && ijet_minDR > -1)
	  {
	    jetFlags[ijet_minDR] = 1;
	    Njetflagged++;
	  }
      }

//     // printouts for jet cleaning
//     STDOUT("CLEANING ----------- v_idx_ele_PtCut_IDISO_noOverlap.size = "<< v_idx_ele_PtCut_IDISO_noOverlap.size() <<", Njetflagged = "<< Njetflagged<<", diff="<< v_idx_ele_PtCut_IDISO_noOverlap.size()-Njetflagged );
//     if( (v_idx_ele_PtCut_IDISO_noOverlap.size()-Njetflagged) == 1 ) 
//       {
// 	TLorentzVector thisele;
// 	for(int iele=0; iele<v_idx_ele_PtCut_IDISO_noOverlap.size(); iele++)
// 	  {
// 	    thisele.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
// 				 ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
// 				 ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),0);
// 	    STDOUT("CLEANING: e"<<iele+1<<" Pt, eta, phi = "  << ", "<<thisele.Pt()<<", "<< thisele.Eta() <<", "<< thisele.Phi());
// 	  }
// 	TLorentzVector thisjet;
// 	for(int ijet=0; ijet<v_idx_jet_PtCut.size(); ijet++)
// 	  {
// 	    thisjet.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut[ijet]),
// 				 JetEta->at(v_idx_jet_PtCut[ijet]),
// 				 JetPhi->at(v_idx_jet_PtCut[ijet]),0);
// 	    STDOUT("CLEANING: j"<<ijet+1<<" Pt, eta, phi = " << ", "<<thisjet.Pt()<<", "<< thisjet.Eta() <<", "<< thisjet.Phi()<<" jetFlags="<<jetFlags[ijet] );
// 	  }
//       } // printouts for jet cleaning
    
    
    float JetEtaCutValue = getPreCutValue1("jet_EtaCut");
    for(int ijet=0; ijet<v_idx_jet_PtCut.size(); ijet++) //pT pre-cut + no overlaps with electrons + jetID
      {	
	bool passjetID = JetPassID->at(v_idx_jet_PtCut[ijet]);
	//bool passjetID = JetIdloose(CaloJetresEMF->at(v_idx_jet_PtCut[ijet]),CaloJetfHPD->at(v_idx_jet_PtCut[ijet]),CaloJetn90Hits->at(v_idx_jet_PtCut[ijet]), CaloJetEta->at(v_idx_jet_PtCut[ijet]));
	// ---- use the flag stored in rootTuples
	//if( (JetOverlaps->at(v_idx_jet_PtCut[ijet]) & 1 << eleIDType) == 0  /* NO overlap with electrons */  
	// ----

	if( jetFlags[ijet] == 0  )                         /* NO overlap with electrons */  
	  //  && passjetID == true )                            /* pass JetID */
	  // && (caloJetOverlaps[ijet] & 1 << 5)==0 )         /* NO overlap with muons */      
	  v_idx_jet_PtCut_noOverlap.push_back(v_idx_jet_PtCut[ijet]);

	if( jetFlags[ijet] == 0                           /* NO overlap with electrons */  
	    && passjetID == true )                            /* pass JetID */
	  // && (caloJetOverlaps[ijet] & 1 << 5)==0 )         /* NO overlap with muons */      
	  v_idx_jet_PtCut_noOverlap_ID.push_back(v_idx_jet_PtCut[ijet]);

	if( jetFlags[ijet] == 0                           /* NO overlap with electrons */  
	    && passjetID == true                             /* pass JetID */
	    && fabs( JetEta->at(v_idx_jet_PtCut[ijet]) ) < JetEtaCutValue )
	  // && (caloJetOverlaps[ijet] & 1 << 5)==0 )         /* NO overlap with muons */      
	  v_idx_jet_PtCut_noOverlap_ID_EtaCut.push_back(v_idx_jet_PtCut[ijet]);

	
	//NOTE: We should verify that caloJetOverlaps match with the code above
      } // End loop over jets
    

    // Muons
    vector<int> v_idx_muon_all;
    vector<int> v_idx_muon_PtCut;
    vector<int> v_idx_muon_PtCut_IDISO;
    
    // Loop over muons  
    for(int imuon=0; imuon<MuonPt->size(); imuon++){

      // no cut on reco muons
      v_idx_muon_all.push_back(imuon);

      if ( (*MuonPt)[imuon] < muon_PtCut) continue;

      // pT pre-cut on muons
      v_idx_muon_PtCut.push_back(imuon);
      
      if ( ((*MuonTrkHits)[imuon]  >= muNHits_minThresh  )&&( fabs((*MuonTrkD0)[imuon]) < muTrkD0Maximum ) &&((*MuonPassIso)[imuon]==1 ) &&((*MuonPassID)[imuon]==1) ) 
	{
	  v_idx_muon_PtCut_IDISO.push_back(imuon);
	}

    }// end loop over muons

//     // vertexes
//     vector<int> v_idx_vertex_good;
//     // loop over vertexes
//     for(int ivertex = 0; ivertex<VertexChi2->size(); ivertex++){
//       if ( !(VertexIsFake->at(ivertex)) 
// 	   && VertexNDF->at(ivertex) > vertexMinimumNDOF 
// 	   && fabs( VertexZ->at(ivertex) ) <= vertexMaxAbsZ
// 	   && fabs( VertexRho->at(ivertex) ) <= vertexMaxd0 )
// 	{
// 	  v_idx_vertex_good.push_back(ivertex);
// 	  //STDOUT("v_idx_vertex_good.size = "<< v_idx_vertex_good.size() );
// 	}
//     }


    // Set the evaluation of the cuts to false and clear the variable values and filled status
    resetCuts();
    

    // Set the value of the variableNames listed in the cutFile to their current value
    
    // Trigger (L1 and HLT)
    if(isData==true)
      {
	fillVariableWithValue( "PassBPTX0", isBPTX0 ) ;
	fillVariableWithValue( "PassPhysDecl", isPhysDeclared ) ;       
      }
    else
      {
	fillVariableWithValue( "PassBPTX0", true ) ;
	fillVariableWithValue( "PassPhysDecl", true ) ;       
      }

    fillVariableWithValue( "PassHLT", PassTrig ) ;
//     fillVariableWithValue( "nVertex_good", v_idx_vertex_good.size() ) ;

    //Event filters at RECO level
    fillVariableWithValue( "PassBeamScraping", !isBeamScraping ) ;
    fillVariableWithValue( "PassPrimaryVertex", isPrimaryVertex ) ;
    //    fillVariableWithValue( "PassHBHENoiseFilter", passLooseNoiseFilter ) ;

    // nMu
//     fillVariableWithValue( "nMu_all", MuonPt->size() ) ;
//     fillVariableWithValue( "nMu_PtCut", v_idx_muon_PtCut.size() ) ;
//     fillVariableWithValue( "nMu_PtCut_IDISO", v_idx_muon_PtCut_IDISO.size() ) ;

    // nEle
    fillVariableWithValue( "nEle_all", v_idx_ele_all.size() ) ;
    fillVariableWithValue( "nEle_PtCut", v_idx_ele_PtCut.size() ) ;
    fillVariableWithValue( "nEle_PtCut_IDISO_noOvrlp", v_idx_ele_PtCut_IDISO_noOverlap.size() ) ;
     
    // nJet
    fillVariableWithValue( "nJet_all", v_idx_jet_all.size() ) ;
    fillVariableWithValue( "nJet_PtCut", v_idx_jet_PtCut.size() ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlp", v_idx_jet_PtCut_noOverlap.size() ) ;
    fillVariableWithValue( "nJet_PtCut_noOvrlp_ID", v_idx_jet_PtCut_noOverlap_ID.size() ) ;
    //TwoEleOnly
    fillVariableWithValue( "nJet_TwoEleOnly_All", v_idx_jet_PtCut_noOverlap_ID.size() ) ;
    fillVariableWithValue( "nJet_TwoEleOnly_EtaCut", v_idx_jet_PtCut_noOverlap_ID_EtaCut.size() ) ;
    //PAS June 2010
    fillVariableWithValue( "nJet_PAS_All", v_idx_jet_PtCut_noOverlap_ID.size() ) ;
    fillVariableWithValue( "nJet_PAS_EtaCut", v_idx_jet_PtCut_noOverlap_ID_EtaCut.size() ) ;

    // MET
    //PAS June 2010
    fillVariableWithValue( "pfMET_PAS", PFMET->at(0) ) ;
    fillVariableWithValue( "tcMET_PAS", TCMET->at(0) ) ;
    fillVariableWithValue( "caloMET_PAS", CaloMET->at(0) ) ;
    fillVariableWithValue( "pfMETPhi_PAS", PFMETPhi->at(0) ) ;

    // 1st ele
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stEle_IDISO_NoOvrlp", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	fillVariableWithValue( "Eta1stEle_IDISO_NoOvrlp", ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	fillVariableWithValue( "mEta1stEle_IDISO_NoOvrlp", fabs(ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0])) );
	//PAS June 2010
	fillVariableWithValue( "Pt1stEle_PAS", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	fillVariableWithValue( "Eta1stEle_PAS", ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
	fillVariableWithValue( "Phi1stEle_PAS", ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) );
      }

    // 2nd ele
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 2 ) 
      {
	fillVariableWithValue( "Pt2ndEle_IDISO_NoOvrlp", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) );
	fillVariableWithValue( "Eta2ndEle_IDISO_NoOvrlp", ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) );
	fillVariableWithValue( "mEta2ndEle_IDISO_NoOvrlp", fabs(ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[1])) );
	fillVariableWithValue( "maxMEtaEles_IDISO_NoOvrl", max( getVariableValue("mEta1stEle_IDISO_NoOvrlp"), getVariableValue("mEta2ndEle_IDISO_NoOvrlp") ) );
	//PAS June 2010
	fillVariableWithValue( "Pt2ndEle_PAS", ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) );
	fillVariableWithValue( "Eta2ndEle_PAS", ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) );
	fillVariableWithValue( "Phi2ndEle_PAS", ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) );
      }

    // 1st jet
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 1 ) 
      {
	fillVariableWithValue( "Pt1stJet_noOvrlp_ID", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "Eta1stJet_noOvrlp_ID", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "mEta1stJet_noOvrlp_ID", fabs(JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0])) );
	//PAS June 2010
	fillVariableWithValue( "Pt1stJet_PAS", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "Eta1stJet_PAS", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
	fillVariableWithValue( "Phi1stJet_PAS", JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]) );
      }


    //cout << "2nd Jet" << endl;
    //## 2nd jet
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2 ) 
      {
	fillVariableWithValue( "Pt2ndJet_noOvrlp_ID", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "Eta2ndJet_noOvrlp_ID", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "mEta2ndJet_noOvrlp_ID", fabs(JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1])) );
	fillVariableWithValue( "maxMEtaJets_noOvrlp_ID", max( getVariableValue("mEta1stJet_noOvrlp_ID"), getVariableValue("mEta2ndJet_noOvrlp_ID") ) );
	//PAS June 2010
	fillVariableWithValue( "Pt2ndJet_PAS", JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "Eta2ndJet_PAS", JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
	fillVariableWithValue( "Phi2ndJet_PAS", JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]) );
      }

    //## define "2ele" and "2jets" booleans
    bool TwoEle=false;
    bool TwoJets=false;
    if( v_idx_ele_PtCut_IDISO_noOverlap.size() >= 2 ) TwoEle = true;
    if( v_idx_jet_PtCut_noOverlap_ID.size() >= 2 ) TwoJets = true;

    // ST
    double calc_sT=-999.; 
    if ( (TwoEle) && (TwoJets) ) 
      {
	calc_sT = 
	  ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) +
	  ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) +
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) +
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]);
	fillVariableWithValue("sT_preliminary", calc_sT);
	fillVariableWithValue("sT", calc_sT);
	fillVariableWithValue("sT_MLQ150", calc_sT);
	fillVariableWithValue("sT_MLQ200", calc_sT);
	fillVariableWithValue("sT_MLQ250", calc_sT);
	fillVariableWithValue("sT_MLQ280", calc_sT);
	fillVariableWithValue("sT_MLQ300", calc_sT);
	fillVariableWithValue("sT_MLQ320", calc_sT);
	fillVariableWithValue("sT_MLQ340", calc_sT);
	fillVariableWithValue("sT_MLQ370", calc_sT);
	fillVariableWithValue("sT_MLQ400", calc_sT);       
	fillVariableWithValue("sT_MLQ450", calc_sT);       
	fillVariableWithValue("sT_MLQ500", calc_sT);       
	//PAS June 2010
	fillVariableWithValue("sT_PAS", calc_sT);
      }

    // ST jets
    if (TwoJets)
      {
	double calc_sTjet=-999.;
	calc_sTjet =
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]) + 
	  JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]) ;
	fillVariableWithValue("sTjet_PAS", calc_sTjet);
      }

    // Mjj
    if (TwoJets)
      {
	TLorentzVector jet1, jet2, jj;
	jet1.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),0);
	jet2.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]),0);
	jj = jet1+jet2;
	//PAS June 2010
	fillVariableWithValue("Mjj_PAS", jj.M());
      }


    // Mee
    if (TwoEle)
      {
	TLorentzVector ele1, ele2, ee;
	ele1.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),0);
	ele2.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[1]),
			  ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[1]),
			  ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[1]),0);
	ee = ele1+ele2;
	fillVariableWithValue("Mee", ee.M());
	//TwoEleOnly
	fillVariableWithValue("Mee_TwoEleOnly", ee.M());
	fillVariableWithValue("Mee_presel", ee.M());
	//
	//PAS June 2010
	fillVariableWithValue("Mee_PAS", ee.M());

	double calc_sTele=-999.;
        calc_sTele =
	   ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]) +
	   ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[1]) ;
	fillVariableWithValue("sTele_PAS", calc_sTele);

// 	if(isData==true) 
// 	  {
// 	    STDOUT("Two electrons: Run, LS, Event = "<<run<<", "<<ls<<", "<<event);
// 	    STDOUT("Two electrons: M_ee, Pt_ee, Eta_ee, Phi_ee = "<<ee.M() <<", "<< ee.Pt() <<", "<< ee.Eta() <<", "<< ee.Phi());
// 	    STDOUT("Two electrons: 1st ele Pt, eta, phi = "<< ele1.Pt() <<", "<< ele1.Eta() <<", "<< ele1.Phi() );
// 	    STDOUT("Two electrons: 2nd ele Pt, eta, phi = "<< ele2.Pt() <<", "<< ele2.Eta() <<", "<< ele2.Phi() );
// 	  }
      }

    // Mej 
    double Me1j1, Me1j2, Me2j1, Me2j2 = -999;
    double deltaM_e1j1_e2j2 = 9999;
    double deltaM_e1j2_e2j1 = 9999;
    double Mej_1stPair = 0;
    double Mej_2ndPair = 0;
    double deltaR_e1j1 ;
    double deltaR_e2j2 ;
    double deltaR_e1j2 ;
    double deltaR_e2j1 ;
    if ( (TwoEle) && (TwoJets) ) // TwoEle and TwoJets
      {
	TLorentzVector jet1, jet2, ele1, ele2;
	ele1.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),
			  ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[0]),0);
	ele2.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[1]),
			  ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[1]),
			  ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[1]),0);
	jet1.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[0]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[0]),0);
	jet2.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetEta->at(v_idx_jet_PtCut_noOverlap_ID[1]),
			  JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[1]),0);
 	cout<<"jet1 Px, Py, Pz, P = "<<jet1.Px()<<" "<<jet1.Py()<<" "<<jet1.Pz()<<" "<<jet1.P()<<endl;
 	cout<<"jet2 Px, Py, Pz, P = "<<jet2.Px()<<" "<<jet2.Py()<<" "<<jet2.Pz()<<" "<<jet2.P()<<endl;
	TLorentzVector e1j1, e1j2, e2j1, e2j2;
	e1j1 = ele1 + jet1;
	e2j2 = ele2 + jet2;
	e2j1 = ele2 + jet1;
	e1j2 = ele1 + jet2;
	Me1j1 = e1j1.M();
	Me2j2 = e2j2.M();
	Me1j2 = e1j2.M();
	Me2j1 = e2j1.M();

	deltaM_e1j1_e2j2 = Me1j1 - Me2j2;
	deltaM_e1j2_e2j1 = Me1j2 - Me2j1;


	double deltaR_e1j1 = ele1.DeltaR(jet1);
	double deltaR_e2j2 = ele2.DeltaR(jet2);
	double deltaR_e1j2 = ele1.DeltaR(jet2);
	double deltaR_e2j1 = ele2.DeltaR(jet1);

	// tests/studies
	TLorentzVector LQ12 = ele1+ele2+jet1+jet2;
	FillUserTH1D("LQ12Pz", LQ12.Pz());
	FillUserTH1D("LQ12Beta", LQ12.Beta());
	FillUserTH1D("LQ12Betaz", LQ12.Pz()/LQ12.E() );
	FillUserTH1D("LQ12Betaz_abs", sqrt(pow(LQ12.Pz()/LQ12.E(),2)) );
	FillUserTH1D("LQ12Betat", LQ12.Pt()/LQ12.E() );
	
	double MR = pow((jet1.E()*jet2.Pz()-jet2.E()*jet1.Pz()),2) / ( pow((jet1.Pz()-jet2.Pz()),2) - pow((jet1.E()-jet2.E()),2) );
	MR = 2*sqrt(MR);
	TVector2 v_MET, ele1T, ele2T, jet1T, jet2T;
	ele1T.Set( ele1.Px(), ele1.Py() );
	ele2T.Set( ele2.Px(), ele2.Py() );
	jet1T.Set( jet1.Px(), jet1.Py() );
	jet2T.Set( jet2.Px(), jet2.Py() );
	v_MET.SetMagPhi( PFMET->at(0), PFMETPhi->at(0) );
	v_MET = v_MET + ele1T + ele2T;
	double MTR = 0.5 * v_MET.Mod() * (jet1T.Mod()+jet2T.Mod()) - 0.5 * ( v_MET*(jet1T+jet2T) ) ;
	MTR = sqrt(MTR);
	FillUserTH1D("MR", MR );
	FillUserTH1D("MTR", MTR );
	FillUserTH1D("Razor", MTR/MR );
	if ( fabs(LQ12.Pz()/LQ12.E())>0.6 ) FillUserTH1D("Razor_LQ12Betaz06", MTR/MR);

// 	// Fill min DR between any of the 2 selected eles and any of the 2 selected jets
// 	double minDR_2ele_2jet = min ( min(deltaR_e1j1,deltaR_e2j2) , min(deltaR_e1j2,deltaR_e2j1) );
// 	fillVariableWithValue("minDR_2ele_2jet", minDR_2ele_2jet);

	if(fabs(deltaM_e1j1_e2j2) > fabs(deltaM_e1j2_e2j1))
	  {
	    Mej_1stPair = Me1j2;
	    Mej_2ndPair = Me2j1;
	    fillVariableWithValue("minDRej_selecPairs", min(deltaR_e1j2,deltaR_e2j1) );
	    fillVariableWithValue("minDRej_unselPairs", min(deltaR_e1j1,deltaR_e2j2) );
	    // tests/studies
	    FillUserTH1D( "LQ1_Eta", e1j2.Eta() );
	    FillUserTH1D( "LQ2_Eta", e2j1.Eta() );
	    FillUserTH1D( "LQ1Eta_plus_LQ2Eta", e1j2.Eta()+e2j1.Eta() );
	    FillUserTH2D( "LQ1Eta_vs_LQ2Eta", e1j2.Eta(), e2j1.Eta() );
	    FillUserTH1D( "LQ1Betaz", e1j2.Pz()/e1j2.E() );
	    FillUserTH1D( "LQ1Betax", e1j2.Px()/e1j2.E() );
	    FillUserTH1D( "LQ1Betay", e1j2.Py()/e1j2.E() );
	    FillUserTH1D( "LQ1Betat", e1j2.Pt()/e1j2.E() );
	    if ( fabs(e1j2.Pt()/e1j2.E())<0.1 ) FillUserTH1D("Razor_LQ1Betat01", MTR/MR);
	  }
	else
	  {
	    Mej_1stPair = Me1j1;
	    Mej_2ndPair = Me2j2;
	    fillVariableWithValue("minDRej_selecPairs", min(deltaR_e1j1,deltaR_e2j2) );
	    fillVariableWithValue("minDRej_unselPairs", min(deltaR_e1j2,deltaR_e2j1) );
	    // tests/studies
	    FillUserTH1D( "LQ1_Eta", e1j1.Eta() );
	    FillUserTH1D( "LQ2_Eta", e2j2.Eta() );
	    FillUserTH1D( "LQ1Eta_plus_LQ2Eta", e1j1.Eta()+e2j2.Eta() );
	  } 


	fillVariableWithValue("Mej_1stPair", Mej_1stPair);       
	fillVariableWithValue("Mej_2ndPair", Mej_2ndPair);
	fillVariableWithValue("DeltaMej", Mej_1stPair-Mej_2ndPair);
	fillVariableWithValue("DeltaMejRel", 2*(Mej_1stPair-Mej_2ndPair)/(Mej_1stPair+Mej_2ndPair));
	fillVariableWithValue("DeltaMejRel_1", 2*(Mej_1stPair-Mej_2ndPair)/(Mej_1stPair+Mej_2ndPair));
	fillVariableWithValue("DeltaMejRel_2", 2*(Mej_1stPair-Mej_2ndPair)/(Mej_1stPair+Mej_2ndPair));
	fillVariableWithValue("DeltaMejRel_3", 2*(Mej_1stPair-Mej_2ndPair)/(Mej_1stPair+Mej_2ndPair));
	//PAS June 2010
	h_Mej_PAS->Fill(Mej_1stPair);
	h_Mej_PAS->Fill(Mej_2ndPair);
	fillVariableWithValue("Mej_1stPair_PAS", Mej_1stPair);       
	fillVariableWithValue("Mej_2ndPair_PAS", Mej_2ndPair);

	// min and max DeltaR between electrons and any jet
	double minDeltaR_ej = 999999;
	double maxDeltaR_ej = -1;
	double thisMinDR, thisMaxDR, DR_thisjet_e1, DR_thisjet_e2;
	TLorentzVector thisjet;
	for(int ijet=0; ijet<v_idx_jet_PtCut_noOverlap_ID.size(); ijet++)
	  {
	    thisjet.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
				 JetEta->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
				 JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),0);
	    DR_thisjet_e1 = thisjet.DeltaR(ele1);
	    DR_thisjet_e2 = thisjet.DeltaR(ele2);
	    thisMinDR = min(DR_thisjet_e1, DR_thisjet_e2);
	    thisMaxDR = max(DR_thisjet_e1, DR_thisjet_e2);
	    if(thisMinDR < minDeltaR_ej)
	      minDeltaR_ej = thisMinDR;
	    if(thisMaxDR > maxDeltaR_ej)
	      maxDeltaR_ej = thisMaxDR;
	  } 
	fillVariableWithValue("minDeltaR_ej", minDeltaR_ej);
	fillVariableWithValue("maxDeltaR_ej", maxDeltaR_ej);

// 	// printouts for TwoElesTwoJets
// 	//	if(isData==true && ( Mej_1stPair<20 || Mej_2ndPair<20 ) ) // printouts for low Mej 
// 	if( isData==true ) 
// 	  {
// 	    STDOUT("TwoElesTwoJets: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
// 	    STDOUT("TwoElesTwoJets: sT, Mee, Mjj_PAS = "<< getVariableValue("sT") <<", "<< getVariableValue("Mee")<<", "<< getVariableValue("Mjj_PAS") );
// 	    STDOUT("TwoElesTwoJets: Mej_1stPair = "<<Mej_1stPair <<", Mej_2ndPair = "<< Mej_2ndPair );
// 	    STDOUT("TwoElesTwoJets: e1j1.M = "<<e1j1.M() <<", e2j2.M = "<<e2j2.M() <<", e1j2.M = "<<e1j2.M()  <<", e2j1.M = "<<e2j1.M()  );
// 	    STDOUT("TwoElesTwoJets: deltaM_e1j1_e2j2 = "<<deltaM_e1j1_e2j2 <<", deltaM_e1j2_e2j1 = "<<deltaM_e1j2_e2j1  );
// 	    STDOUT("TwoElesTwoJets: deltaR_e1j1 = "<<deltaR_e1j1 <<", deltaR_e2j2 = "<<deltaR_e2j2 <<", deltaR_e1j2 = "<<deltaR_e1j2  <<", deltaR_e2j1 = "<<deltaR_e2j1  );
// 	    TLorentzVector thisele;
// 	    for(int iele=0; iele<v_idx_ele_PtCut_IDISO_noOverlap.size(); iele++)
// 	      {
// 		thisele.SetPtEtaPhiM(ElectronPt->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
// 				     ElectronEta->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),
// 				     ElectronPhi->at(v_idx_ele_PtCut_IDISO_noOverlap[iele]),0);
// 		STDOUT("TwoElesTwoJets: e"<<iele+1<<" Pt, eta, phi = "<<thisele.Pt()<<", "<< thisele.Eta() <<", "<< thisele.Phi()<<"; DR_j1, DR_j2 = "<< thisele.DeltaR(jet1)<<", "<<thisele.DeltaR(jet2));
// 	      }
// 	    TLorentzVector thisjet;
// 	    TLorentzVector thisjet_e1, thisjet_e2;
// 	    for(int ijet=0; ijet<v_idx_jet_PtCut_noOverlap_ID.size(); ijet++)
// 	      {
// 		thisjet.SetPtEtaPhiM(JetPt->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
// 				     JetEta->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),
// 				     JetPhi->at(v_idx_jet_PtCut_noOverlap_ID[ijet]),0);
// 		thisjet_e1 = thisjet + ele1;
// 		thisjet_e2 = thisjet + ele2;
// 		STDOUT("TwoElesTwoJets: j"<<ijet+1<<" Pt, eta, phi = "<<thisjet.Pt()<<", "<< thisjet.Eta() <<", "<< thisjet.Phi()<<"; DR_e1, DR_e2 = "<< thisjet.DeltaR(ele1)<<", "<<thisjet.DeltaR(ele2) << "; M_e1, M_e2 = " <<thisjet_e1.M() <<", "<<thisjet_e2.M() );
// 	      }
// 	  } // printouts for TwoElesTwoJets:

      } // TwoEle and TwoJets



    // Evaluate cuts (but do not apply them)
    evaluateCuts();

    // Fill histograms and do analysis based on cut evaluation
    //h_nEleFinal->Fill( ElectronPt->size() );
     

    //Large MET events passing pre-selection
    float METcut = 60;
    if( variableIsFilled("pfMET_PAS") && passedAllPreviousCuts("pfMET_PAS") && isData==true) 
      {

	if( getVariableValue("pfMET_PAS") > METcut )
	  {

	    STDOUT("pfMET>60GeV: ----------- START ------------");

	    STDOUT("pfMET>60GeV: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);	    
	    if( variableIsFilled("Pt1stEle_PAS") && variableIsFilled("Eta1stEle_PAS") )	      
	      STDOUT("pfMET>60GeV: Pt1stEle_PAS,Eta1stEle_PAS = "<<getVariableValue("Pt1stEle_PAS")<<",\t"<<getVariableValue("Eta1stEle_PAS"));
	    if( variableIsFilled("Pt2ndEle_PAS") && variableIsFilled("Eta2ndEle_PAS") )	      
	      STDOUT("pfMET>60GeV: Pt2ndEle_PAS,Eta2ndEle_PAS = "<<getVariableValue("Pt2ndEle_PAS")<<",\t"<<getVariableValue("Eta2ndEle_PAS"));
	    if( variableIsFilled("Pt1stJet_PAS") && variableIsFilled("Eta1stJet_PAS") )	      
	      STDOUT("pfMET>60GeV: Pt1stJet_PAS,Eta1stJet_PAS = "<<getVariableValue("Pt1stJet_PAS")<<",\t"<<getVariableValue("Eta1stJet_PAS"));
	    if( variableIsFilled("Pt2ndJet_PAS") && variableIsFilled("Eta2ndJet_PAS") )	      
	      STDOUT("pfMET>60GeV: Pt2ndJet_PAS,Eta2ndJet_PAS = "<<getVariableValue("Pt2ndJet_PAS")<<",\t"<<getVariableValue("Eta2ndJet_PAS"));
	    if( variableIsFilled("Mee_PAS") && variableIsFilled("Mjj_PAS") )	      
	      STDOUT("pfMET>60GeV: Mee_PAS,Mjj_PAS = "<<getVariableValue("Mee_PAS")<<",\t"<<getVariableValue("Mjj_PAS"));
	    if( variableIsFilled("Mej_1stPair_PAS") && variableIsFilled("Mej_2ndPair_PAS") )	      
	      STDOUT("pfMET>60GeV: Mej_1stPair_PAS,Mej_2ndPair_PAS = "<<getVariableValue("Mej_1stPair_PAS")
		     <<",\t"<<getVariableValue("Mej_2ndPair_PAS"));
	    if( variableIsFilled("pfMET_PAS") && variableIsFilled("caloMET_PAS") )	      
	      STDOUT("pfMET>60GeV: pfMET_PAS,caloMET_PAS = "<<getVariableValue("pfMET_PAS")<<",\t"<<getVariableValue("caloMET_PAS"));

	    STDOUT("pfMET>60GeV: ------------ END -------------");

	  }
      }


    // printouts for TwoElesTwoJets
    if( isData==true && passedAllPreviousCuts("Mee") ) 
      {
	STDOUT("TwoElesTwoJets: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event);
	STDOUT("TwoElesTwoJets: sT, Mee, Mjj_PAS = "<< getVariableValue("sT") <<", "<< getVariableValue("Mee")<<", "<< getVariableValue("Mjj_PAS") );
	STDOUT("TwoElesTwoJets: Mej_1stPair, Mej_2ndPair = "<< getVariableValue("Mej_1stPair")<<", "<< getVariableValue("Mej_2ndPair") );
	STDOUT("TwoElesTwoJets: 1stEle Pt, eta, phi = "<<getVariableValue("Pt1stEle_PAS") <<", "<< getVariableValue("Eta1stEle_PAS") <<", "<< getVariableValue("Phi1stEle_PAS") );
	STDOUT("TwoElesTwoJets: 2ndEle Pt, eta, phi = "<<getVariableValue("Pt2ndEle_PAS") <<", "<< getVariableValue("Eta2ndEle_PAS") <<", "<< getVariableValue("Phi2ndEle_PAS") );
	STDOUT("TwoElesTwoJets: 1stJet Pt, eta, phi = "<<getVariableValue("Pt1stJet_PAS") <<", "<< getVariableValue("Eta1stJet_PAS") <<", "<< getVariableValue("Phi1stJet_PAS") );
	STDOUT("TwoElesTwoJets: 2ndJet Pt, eta, phi = "<<getVariableValue("Pt2ndJet_PAS") <<", "<< getVariableValue("Eta2ndJet_PAS") <<", "<< getVariableValue("Phi2ndJet_PAS") );
	if ( passedCut("Mee") )
	  {
	    STDOUT("PassedMeeAndAllPrevious: Run, LS, Event = "<<run<<",\t"<<ls<<",\t"<<event<<", sT = "<< getVariableValue("sT"));
	  }
      } // printouts for TwoElesTwoJets:

    if( passedAllPreviousCuts("sT_preliminary") && passedCut("Mee") )
      {
	FillUserTH1D( "sT_fullSeleNoPreStCut", getVariableValue("sT") );
      }

    //INFO
    //      // retrieve value of previously filled variables (after making sure that they were filled)
    //      double totpTEle;
    //      if ( variableIsFilled("pT1stEle") && variableIsFilled("pT2ndEle") ) 
    //        totpTEle = getVariableValue("pT1stEle")+getVariableValue("pT2ndEle");
    //      // reject events that did not pass level 0 cuts
    //      if( !passedCut("0") ) continue;
    //      // ......

    //////////////////////////////////////////////////////////////////
    // Solve system of non-linear n equations with n variables - BEGIN
    // See http://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Root_002dFinding.html
    // and http://www.gnu.org/software/gsl/manual/html_node/Example-programs-for-Multidimensional-Root-finding.html

    if(passedAllPreviousCuts("sT_preliminary"))
      {
	const gsl_multiroot_fsolver_type *T;
	gsl_multiroot_fsolver *s;
  
    
	const size_t n = 4;

	double pxj1 = getVariableValue("Pt1stJet_PAS")*cos(getVariableValue("Phi1stJet_PAS"));
	double pyj1 = getVariableValue("Pt1stJet_PAS")*sin(getVariableValue("Phi1stJet_PAS"));
	double thj1 = 2*atan(exp(-1*getVariableValue("Eta1stJet_PAS")));
	double pzj1 = getVariableValue("Pt1stJet_PAS")/tan(thj1);
	STDOUT("pxj1, pyj1, pzj1 = "<<pxj1<<" "<<pyj1<<" "<<pzj1<<" thj1 = "<<thj1);
	
	double pxj2 = getVariableValue("Pt2ndJet_PAS")*cos(getVariableValue("Phi2ndJet_PAS"));
	double pyj2 = getVariableValue("Pt2ndJet_PAS")*sin(getVariableValue("Phi2ndJet_PAS"));
	double thj2 = 2*atan(exp(-1*getVariableValue("Eta2ndJet_PAS")));
	double pzj2 = getVariableValue("Pt2ndJet_PAS")/tan(thj2);
	STDOUT("pxj1, pyj1, pzj1 = "<<pxj2<<" "<<pyj2<<" "<<pzj2<<" thj2 = "<<thj2);
	
// 	double metx =  getVariableValue("pfMET_PAS")*cos(getVariableValue("pfMETPhi_PAS"));
// 	double mety =  getVariableValue("pfMET_PAS")*sin(getVariableValue("pfMETPhi_PAS"));
// 	STDOUT("MET, METPhi, metx, mety = "<<getVariableValue("pfMET_PAS")<<" "<<getVariableValue("pfMETPhi_PAS")<<" "<<metx<<" "<<mety)	
	double metx =  getVariableValue("Pt1stEle_PAS")*cos(getVariableValue("Phi1stEle_PAS"))+getVariableValue("Pt2ndEle_PAS")*cos(getVariableValue("Phi2ndEle_PAS"));
	double mety =  getVariableValue("Pt1stEle_PAS")*sin(getVariableValue("Phi1stEle_PAS"))+getVariableValue("Pt2ndEle_PAS")*sin(getVariableValue("Phi2ndEle_PAS"));

	///////////////////////////////////////////////////////////////////////////
	//test of using a well-defined event with no z-boost of the LQ1+LQ2 system
	pxj1=430; pyj1=0; pzj1=0; pxj2=0; pyj2=430; pzj2=0; 
	TLorentzVector j1,j2;
	j1.SetPxPyPzE(pxj1,pyj1,pzj1,sqrt(pxj1*pxj1+pyj1*pyj1+pzj1*pzj1));
	j2.SetPxPyPzE(pxj2,pyj2,pzj2,sqrt(pxj2*pxj2+pyj2*pyj2+pzj2*pzj2));
	TLorentzVector nu1,nu2;
	nu1=-j1; nu1.SetE(j1.E());
	nu2=-j2; nu2.SetE(j2.E());
	double bx,by,bz; 
	bx=-0.3;
	by=0.;
	bz=0.;
	j1.Boost( bx, by, bz);
	j2.Boost(-bx,-by,-bz);
	nu1.Boost( bx, by, bz);
	nu2.Boost(-bx,-by,-bz);
	pxj1=j1.Px(); 
	pyj1=j1.Py(); 
	pzj1=j1.Pz(); 
	pxj2=j2.Px(); 
	pyj2=j2.Py(); 
	pzj2=j2.Pz(); 
	metx=nu1.Px()+nu2.Px();
	mety=nu1.Py()+nu2.Py();
        TLorentzVector LQ1,LQ2;
	LQ1=j1+nu1;
	LQ2=j2+nu2;
	TLorentzVector j1j2; j1j2 = j1+j2;
	TLorentzVector nu1nu2; nu1nu2 = nu1+nu2;
	cout<<"set values of bx, by, bz and b = "<<bx<<"\t"<<by<<"\t"<<bz<<"\t"<<sqrt(bx*bx+by*by+bz*bz)<<endl;
	cout<<"set values of j1     Px, Py, Pz, P, E = "<<j1.Px()<<"\t"<<j1.Py()<<"\t"<<j1.Pz()<<"\t"<<j1.P()<<"\t"<<j1.E()<<endl;
	cout<<"set values of j2     Px, Py, Pz, P, E = "<<j2.Px()<<"\t"<<j2.Py()<<"\t"<<j2.Pz()<<"\t"<<j2.P()<<"\t"<<j2.E()<<endl;
	cout<<"set values of nu1    Px, Py, Pz, P, E = "<<nu1.Px()<<"\t"<<nu1.Py()<<"\t"<<nu1.Pz()<<"\t"<<nu1.P()<<"\t"<<nu1.E()<<endl;
	cout<<"set values of nu2    Px, Py, Pz, P, E = "<<nu2.Px()<<"\t"<<nu2.Py()<<"\t"<<nu2.Pz()<<"\t"<<nu2.P()<<"\t"<<nu2.E()<<endl;
	cout<<"set values of metx mety = "<<metx<<"\t"<<mety<<endl;
	cout<<"set values of j1j2   Px, Py, Pz, P, E, M = "<<j1j2.Px()<<"\t"<<j1j2.Py()<<"\t"<<j1j2.Pz()<<"\t"<<j1j2.P()<<"\t"<<j1j2.E()<<"\t"<<j1j2.M()<<endl;
	cout<<"set values of nu1nu2   Px, Py, Pz, P, E, M = "<<nu1nu2.Px()<<"\t"<<nu1nu2.Py()<<"\t"<<nu1nu2.Pz()<<"\t"<<nu1nu2.P()<<"\t"<<nu1nu2.E()<<"\t"<<nu1nu2.M()<<endl;
	cout<<"set values of LQ1    Px, Py, Pz, P, E, M = "<<LQ1.Px()<<"\t"<<LQ1.Py()<<"\t"<<LQ1.Pz()<<"\t"<<LQ1.P()<<"\t"<<LQ1.E()<<"\t"<<LQ1.M()<<endl;
	cout<<"set values of LQ2    Px, Py, Pz, P, E, M = "<<LQ2.Px()<<"\t"<<LQ2.Py()<<"\t"<<LQ2.Pz()<<"\t"<<LQ2.P()<<"\t"<<LQ2.E()<<"\t"<<LQ2.M()<<endl;
	// end of using a well-defined event with no z-boost of the LQ1+LQ2 system
	//////////////////////////////////////////////////////////////////////////

	struct rparams p = {pxj1, pyj1, pzj1, pxj2, pyj2, pzj2, metx, mety};
	gsl_multiroot_function f = {&rosenbrock_f, n, &p};
	
	double x_init[4] = {150, 100, 150, 1000};
	gsl_vector *x = gsl_vector_alloc (n);

	vector<double> MLQs;
	// loop over a grid of starting points
	for (x_init[0] = -100; x_init[0] <= 100; x_init[0] += 100) {
	  for (x_init[1] = -100; x_init[1] <= 100; x_init[1] += 100) {
	    for (x_init[2] = -100; x_init[2] <= 100; x_init[2] += 100) {
	      for (x_init[3] =  100; x_init[3] <= 1000; x_init[3] += 300) {
		
		int status;
		size_t iter = 0;
		
		gsl_vector_set (x, 0, x_init[0]);
		gsl_vector_set (x, 1, x_init[1]);
		gsl_vector_set (x, 2, x_init[2]);
		gsl_vector_set (x, 3, x_init[3]);
		
		T = gsl_multiroot_fsolver_hybrids;
		s = gsl_multiroot_fsolver_alloc (T, n);
		gsl_multiroot_fsolver_set (s, &f, x);
		
		print_state (iter, s);
		
		do
		  {
		    iter++;
		    status = gsl_multiroot_fsolver_iterate (s);
		    
		    print_state (iter, s);
		    
		    if (status)   /* check if solver is stuck */
		      break;
		    
		    status =
		      gsl_multiroot_test_residual (s->f, 1e-7);
		  }
		while (status == GSL_CONTINUE && iter < 1000);
		
		printf ("status = %s\n", gsl_strerror (status));
		
		//	gsl_multiroot_fsolver_free (s);
		//	gsl_vector_free (x);
		
		if(status == GSL_SUCCESS)  MLQs.push_back(gsl_vector_get (s->x, 3));
		std::cout << "The End" << std::endl;
	      }
	    }
	  }
	} // end loop over a grid of starting points

	cout<<"MLQs = ";
	for (vector<double>::iterator it=MLQs.begin(); it<MLQs.end(); it++){
	  cout<<*it<<" ";
	}
	cout <<endl;
	
      }
    // Solve system of non-linear n equations with n variables - END
    ////////////////////////////////////////////////////////////////
    
    
    
    ////////////////////// User's code to be done for every event - END ///////////////////////
    
  } // End of loop over events
  

  ////////////////////// User's code to write histos - BEGIN ///////////////////////

  h_Mej_PAS->Write();

  ////////////////////// User's code to write histos - END ///////////////////////
  
  
  //STDOUT("analysisClass::Loop() ends");   
}
