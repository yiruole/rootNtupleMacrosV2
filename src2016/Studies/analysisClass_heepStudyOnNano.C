#define analysisClass_cxx
#include "analysisClass.h"
#include <typeinfo>
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TLeaf.h>

#include "Collection.h"
#include "GenParticle.h"
#include "Electron.h"
#include "LooseElectron.h"
#include "Muon.h"
#include "PFJet.h"
#include "GenJet.h"
#include "HLTriggerObject.h"
#include "HLTriggerObjectCollectionHelper.h"
#include "HistoReader.h"
#include "ElectronScaleFactors.C"

//--------------------------------------------------------------------------
// Function for trigger matching 
//--------------------------------------------------------------------------

template < class Object1, class Object2 > 
double triggerMatchPt ( const CollectionPtr & collection, Object2 & target_object, double delta_r_cut, bool verbose=false ){
  double matched_pt = -999.0;
  if ( collection ) { 
    int size = collection -> GetSize();
    if ( size > 0 ){ 
      if(verbose) {
        std::cout << "triggerMatchPt(): try to find closest object in DR to object: " << target_object << "from collection: " << std::endl;
        collection->examine<HLTriggerObject>("trigObjs");
      }
      Object1 matched_object = collection -> GetClosestInDR <Object1, Object2> ( target_object );
      double dr = matched_object.DeltaR ( & target_object );
      if(verbose) {
        std::cout << "found matched_object: " << matched_object << " with dR=" << dr << std::endl;
      }
      if ( dr < delta_r_cut ) { 
        matched_pt = matched_object.Pt();
        if(verbose) {
          std::cout << "dr=" << dr << " < delta_r_cut=" << delta_r_cut << ", so matched_pt set to matched_object pt: " << matched_object << std::endl;
        }
      }
    }
  }
  return matched_pt;
}

bool isElectronBarrel(float superclusterEta) {
  return fabs(superclusterEta) < 1.4442;
}

bool isElectronEndcap(float superclusterEta) {
  return fabs(superclusterEta) > 1.566 && fabs(superclusterEta) < 2.5;
}

int getJetLQMotherIdx(CollectionPtr c_genJet_all, CollectionPtr c_pfjet_final, int jetRawIndex, CollectionPtr c_genQuark_hardProcess) {
  GenJet matchedGenJet = c_genJet_all->GetConstituent<GenJet>(c_pfjet_final->GetConstituent<PFJet>(jetRawIndex).MatchedGenJetIndex());
  GenParticle matchedGenParticle = c_genQuark_hardProcess->GetClosestInDR<GenParticle, GenJet>(matchedGenJet);
  if(matchedGenParticle.Pt() > 0) {
    matchedGenParticle.PassUserID(GEN_FROM_LQ);
    return matchedGenParticle.MotherLQIndex();
  }
  else // consider this a failed match
    return -999;
}

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile) 
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile)
{
  std::cout << "analysisClass::analysisClass(): begins " << std::endl;
  std::cout << "analysisClass::analysisClass(): ends " << std::endl;
}

analysisClass::~analysisClass()
{ 
  std::cout << "analysisClass::~analysisClass(): begins " << std::endl;
  std::cout << "analysisClass::~analysisClass(): ends " << std::endl;
}

void analysisClass::Loop()
{
  std::cout << "analysisClass::Loop(): begins " << std::endl;
  
  //--------------------------------------------------------------------------
  // Verbose? Or not?
  //--------------------------------------------------------------------------
  
  bool verbose = !true;
  
  //--------------------------------------------------------------------------
  // Decide which plots to save (default is to save everything)
  //--------------------------------------------------------------------------
  
  fillSkim                         ( !true  ) ;
  fillAllPreviousCuts              ( !true  ) ;
  fillAllOtherCuts                 ( !true  ) ;
  fillAllSameLevelAndLowerLevelCuts( true  ) ;
  fillAllCuts                      ( !true  ) ;
  fillAllSameLevelCuts             ( true  ) ;

  //--------------------------------------------------------------------------
  // Plots
  //--------------------------------------------------------------------------
  CreateUserTH1D( "nJetFinal"                       ,    10    , 0   , 10     );
  CreateUserTH1D( "nJetFinalMatchedLQ"              ,    10    , 0   , 10     );
  CreateUserTH1D( "idxJetFinalMatchedLQ"            ,    10    , 0   , 10     );
  CreateUserTH1D( "genJetsNearEGLooseEle_matchedToLQ", 2, 0, 2);
  CreateUserTH1D( "genJetsNearEGLooseEle_matchedToLQ_elePassHeep", 2, 0, 2);
  CreateUserTH1D( "genJetsNearEGLooseEle_unmatchedToLQ_elePassHeep", 2, 0, 2);
  CreateUserTH1D( "nGenJetsNearEGLooseEle"          ,    10    , 0   , 10     );
  CreateUserTH1D( "DeltaPtMatchJetEle_elePassHeep"                   ,    200    , 0   , 1000    );
  CreateUserTH1D( "DeltaPtMatchJetEle_eleFailHeep"                   ,    200    , 0   , 1000    );
  CreateUserTH1D( "DeltaPtMatchJetEle_eleFailHeepCaloIso"                   ,    200    , 0   , 1000    );
  CreateUserTH1D( "DeltaRMatchJetEle_elePassHeep"                   ,    100    , 0   , 1    );
  CreateUserTH1D( "DeltaRMatchJetEle_eleFailHeep"                   ,    100    , 0   , 1    );
  CreateUserTH1D( "DeltaRMatchJetEle_eleFailHeepCaloIso"                   ,    100    , 0   , 1    );
  CreateUserTH2D( "DeltaRVsDeltaPtMatchJetEle_elePassHeep", 200, 0, 1000, 100, 0, 1);
  CreateUserTH2D( "DeltaRVsDeltaPtMatchJetEle_eleFailHeep", 200, 0, 1000, 100, 0, 1);
  CreateUserTH2D( "DeltaRVsDeltaPtMatchJetEle_eleFailHeepCaloIso", 200, 0, 1000, 100, 0, 1);
  CreateUserTH1D( "JetFinalPtJet1"                  ,    200    , 0   , 1000    );
  CreateUserTH1D( "JetFinalPtJet2"                  ,    200    , 0   , 1000    );
  CreateUserTH1D( "JetFinalPtJet3"                  ,    200    , 0   , 1000    );
  CreateUserTH1D( "JetFinalPtJet4"                  ,    200    , 0   , 1000    );
  CreateUserTH1D( "JetFinalPtJet5"                  ,    200    , 0   , 1000    );
  CreateUserTH1D( "JetFinalPtJet6"                  ,    200    , 0   , 1000    );
  CreateUserTH1D( "nEleReco"                       ,    10    , 0   , 10     );
  CreateUserTH1D( "nEleReco_ptCut"                       ,    10    , 0   , 10     );
  CreateUserTH1D( "nEleEGLoosePt35EtaCut"                       ,    10    , 0   , 10     );
  CreateUserTH1D( "nEleEGLoosePt35"                       ,    10    , 0   , 10     );
  CreateUserTH1D( "nEleEGLoose"                       ,    10    , 0   , 10     );
  CreateUserTH1D( "nEleHEEP"                       ,    10    , 0   , 10     );
  CreateUserTH1D( "nEleFinal"                       ,    10    , 0   , 10     );
  CreateUserTH1D( "ele1PassHeep"                   ,     2    , 0   , 2     );
  CreateUserTH1D( "ele2PassHeep"                   ,     2    , 0   , 2     );
  CreateUserTH1D( "ele3PassHeep"                   ,     2    , 0   , 2     );
  CreateUserTH2D( "nEleNTupleVsNeleRsk"            ,     10, 0, 10, 10, 0, 10) ;
  CreateUserTH1D( "heepElectronsCategory"          ,    10    , 0   , 10     );
  CreateUserTH1D( "deltaREleJetsPassHEEP"         ,    100    , 0   , 10     );
  CreateUserTH1D( "deltaREleJetsFailHEEP"         ,    100    , 0   , 10     );
  CreateUserTH1D( "deltaREleJetsBarrelPassHEEP"         ,    100    , 0   , 10     );
  CreateUserTH1D( "deltaREleJetsBarrelFailHEEP"         ,    100    , 0   , 10     );
  CreateUserTH1D( "deltaREleJetsEndcapPassHEEP"         ,    100    , 0   , 10     );
  CreateUserTH1D( "deltaREleJetsEndcapFailHEEP"         ,    100    , 0   , 10     );
  CreateUserTH1D( "deltaREle1Jets"         ,    100    , 0   , 10     );
  CreateUserTH1D( "deltaREle2Jets"         ,    100    , 0   , 10     );
  CreateUserTH1D( "deltaREle2JetsExcludeClosestJet"         ,    100    , 0   , 10     );
  CreateUserTH1D( "deltaREle2PassHEEPJetsExcludeClosestJet"         ,    100    , 0   , 10     );
  CreateUserTH1D( "deltaREle2FailHEEPCaloIsoJetsExcludeClosestJet"         ,    100    , 0   , 10     );
  CreateUserTH1D( "deltaREle2Jets_highET", 100, 0, 10);
  CreateUserTH1D( "deltaREle2Jets_highET_ExcludeClosestJet", 100, 0, 10);
  CreateUserTH1D( "deltaREle2PassHEEPJets_highET_ExcludeClosestJet", 100, 0, 10);
  CreateUserTH1D( "deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet", 100, 0, 10);
  CreateUserTH1D( "deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_sameLQMother", 100, 0, 10);
  CreateUserTH1D( "deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_diffLQMother", 100, 0, 10);
  CreateUserTH1D( "deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_noLQMother", 100, 0, 10);
  CreateUserTH1D( "deltaREle2Jets_highET_barrel", 100, 0, 10);
  CreateUserTH1D( "deltaREle2Jets_highET_endcap", 100, 0, 10);
  CreateUserTH1D( "deltaREle2Jets_lowET", 100, 0, 10);
  CreateUserTH1D( "deltaREle2Jets_lowET_ExcludeClosestJet", 100, 0, 10);
  CreateUserTH1D( "deltaREle2PassHEEPJets_lowET_ExcludeClosestJet", 100, 0, 10);
  CreateUserTH1D( "deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet", 100, 0, 10);
  CreateUserTH1D( "deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_sameLQMother", 100, 0, 10);
  CreateUserTH1D( "deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_diffLQMother", 100, 0, 10);
  CreateUserTH1D( "deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_noLQMother", 100, 0, 10);
  CreateUserTH1D( "deltaREle2Jets_lowET_barrel", 100, 0, 10);
  CreateUserTH1D( "deltaREle2Jets_lowET_endcap", 100, 0, 10);
  CreateUserTH1D( "deltaR_Ele1JetSameLQ", 100, 0, 10);
  CreateUserTH1D( "deltaR_Ele2JetSameLQ", 100, 0, 10);
  CreateUserTH1D( "deltaR_Ele2JetSameLQ_excludeClosestJet", 100, 0, 10);
  //
  CreateUserTH2D( "deltaREle2JetVsPt_sameLQmother", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_diffLQmother", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_noLQmother", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_passHEEP_sameLQmother", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_passHEEP_diffLQmother", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_passHEEP_noLQmother", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_failHEEP_sameLQmother", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_failHEEP_diffLQmother", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_failHEEP_noLQmother", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_sameLQmother_excludeClosestJet", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_diffLQmother_excludeClosestJet", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_noLQmother_excludeClosestJet", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_passHEEP_sameLQmother_excludeClosestJet", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_passHEEP_diffLQmother_excludeClosestJet", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_passHEEP_noLQmother_excludeClosestJet", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_failHEEP_sameLQmother_excludeClosestJet", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_failHEEP_diffLQmother_excludeClosestJet", 200, 0, 2000, 100, 0, 10);
  CreateUserTH2D( "deltaREle2JetVsPt_failHEEP_noLQmother_excludeClosestJet", 200, 0, 2000, 100, 0, 10);
  //
  CreateUserTH1D( "EMHadD1IsoElectron1"                 ,    1000    , 0   , 100     );
  CreateUserTH1D( "EMHadD1IsoElectron2"                 ,    1000    , 0   , 100     );
  CreateUserTH1D( "EMHadD1IsoElectron1PassHEEP"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "EMHadD1IsoElectron1FailHEEP"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "EMHadD1IsoElectron2PassHEEP"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "EMHadD1IsoElectron2FailHEEP"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "EMHadD1IsoBarrel"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "EMHadD1IsoBarrelPassHEEP"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "EMHadD1IsoBarrelFailHEEP"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "EMHadD1IsoEndcap"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "EMHadD1IsoEndcapPassHEEP"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "EMHadD1IsoEndcapFailHEEP"         ,    1000    , 0   , 100     );
  //
  CreateUserTH1D( "TrkIsoElectron1"                 ,    1000    , 0   , 100     );
  CreateUserTH1D( "TrkIsoElectron2"                 ,    1000    , 0   , 100     );
  CreateUserTH1D( "TrkIsoElectron1PassHEEP"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "TrkIsoElectron1FailHEEP"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "TrkIsoElectron2PassHEEP"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "TrkIsoElectron2FailHEEP"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "TrkIsoBarrel"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "TrkIsoBarrelPassHEEP"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "TrkIsoBarrelFailHEEP"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "TrkIsoEndcap"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "TrkIsoEndcapPassHEEP"         ,    1000    , 0   , 100     );
  CreateUserTH1D( "TrkIsoEndcapFailHEEP"         ,    1000    , 0   , 100     );
  //
  CreateUserTH1D( "ECALDrivenElectron1"                   ,     2    , 0   , 2     );
  CreateUserTH1D( "ECALDrivenElectron2"                   ,     2    , 0   , 2     );
  CreateUserTH1D( "ECALDrivenElectron1PassHEEP"           ,     2    , 0   , 2     );
  CreateUserTH1D( "ECALDrivenElectron1FailHEEP"           ,     2    , 0   , 2     );
  CreateUserTH1D( "ECALDrivenElectron2PassHEEP"           ,     2    , 0   , 2     );
  CreateUserTH1D( "ECALDrivenElectron2FailHEEP"           ,     2    , 0   , 2     );
  // SC eta
  CreateUserTH1D( "SCEtaElectron1"                           ,     100    , -5   , 5     );
  CreateUserTH1D( "SCEtaElectron2"                           ,     100    , -5   , 5     );
  CreateUserTH1D( "SCEtaElectron1PassHEEP"                   ,     100    , -5   , 5     );
  CreateUserTH1D( "SCEtaElectron1FailHEEP"                   ,     100    , -5   , 5     );
  CreateUserTH1D( "SCEtaElectron2PassHEEP"                   ,     100    , -5   , 5     );
  CreateUserTH1D( "SCEtaElectron2FailHEEP"                   ,     100    , -5   , 5     );
  CreateUserTH1D( "SCEtaElectron"                           ,     100    , -5   , 5     );
  CreateUserTH1D( "SCEtaElectronPassHEEP"                           ,     100    , -5   , 5     );
  CreateUserTH1D( "SCEtaElectronFailHEEP"                           ,     100    , -5   , 5     );
  // SC Et
  CreateUserTH1D( "SCEtElectron1"                 ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtElectron2"                 ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtElectron"                 ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtElectronPassHEEP"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtElectronPassLoose"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtElectronFailHEEP"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtElectronFailLoose"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtElectron1PassHEEP"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtElectron1PassLoose"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtElectron2PassHEEP"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtElectron2PassLoose"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtElectron1FailHEEP"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtElectron1FailLoose"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtElectron2FailHEEP"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtElectron2FailLoose"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtBarrel"                 ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtBarrelPassHEEP"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtBarrelPassLoose"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtBarrelFailHEEP"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtBarrelFailLoose"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtEndcap"                 ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtEndcapPassHEEP"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtEndcapFailHEEP"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtEndcapPassLoose"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEtEndcapFailLoose"         ,    200    , 0   , 2000    );
  // SC Energy
  CreateUserTH1D( "SCEnergyElectron1"                 ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEnergyElectron2"                 ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEnergyElectron"                 ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEnergyElectronPassHEEP"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEnergyElectronFailHEEP"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEnergyBarrel"                 ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEnergyBarrelPassHEEP"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEnergyBarrelFailHEEP"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEnergyEndcap"                 ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEnergyEndcapPassHEEP"         ,    200    , 0   , 2000    );
  CreateUserTH1D( "SCEnergyEndcapFailHEEP"         ,    200    , 0   , 2000    );

  /*//------------------------------------------------------------------
   *
   *
   *      
   *      Get all Pre-cut values!
   *
   *
   *
   *///-----------------------------------------------------------------

  //--------------------------------------------------------------------------
  // What reduced skim type?
  // - 0: QCD   (loose electron)
  // - 1: enujj (HEEP electron) 
  // - 2: eejj  (HEEP electron)
  // - 3: single electron (HEEP)
  // - 4: single muon (tight)
  //--------------------------------------------------------------------------

  int reducedSkimType = getPreCutValue1("reducedSkimType");

  //--------------------------------------------------------------------------
  // Trigger scale factors
  //--------------------------------------------------------------------------
  std::string trigSFFileName = getPreCutString1("TriggerSFFileName");
  //std::string graphName = getPreCutString1("TriggerEfficiencyGraphName");
  HistoReader triggerScaleFactorReader(trigSFFileName,"SF_TH2F_Barrel","SF_TH2F_EndCap",false,false);

  //--------------------------------------------------------------------------
  // reco scale factors
  //--------------------------------------------------------------------------
  std::string recoSFFileName = getPreCutString1("RecoSFFileName");
  std::unique_ptr<HistoReader> recoScaleFactorReader = std::unique_ptr<HistoReader>(new HistoReader(recoSFFileName,"EGamma_SF2D","EGamma_SF2D",true,false));

  //--------------------------------------------------------------------------
  // Should we do any scaling / smearing for systematics?
  //--------------------------------------------------------------------------

  int electron_energy_scale_sign = int(getPreCutValue1("electron_energy_scale_sign" ));
  int pfjet_energy_scale_sign    = int(getPreCutValue1("pfjet_energy_scale_sign"    ));

  bool do_eer = bool ( int(getPreCutValue1("do_electron_energy_smear"   )) == 1 );
  bool do_jer = bool ( int(getPreCutValue1("do_pfjet_energy_smear"      )) == 1 );
  bool do_ees = bool ( electron_energy_scale_sign != 0 );
  bool do_jes = bool ( pfjet_energy_scale_sign    != 0 );
  
  //--------------------------------------------------------------------------
  // Analysis year
  //--------------------------------------------------------------------------
  int analysisYear = getPreCutValue1("AnalysisYear");

  //--------------------------------------------------------------------------
  // Cuts for physics objects selection
  //--------------------------------------------------------------------------

  // jet cuts
  double jet_PtCut               = getPreCutValue1("jet_PtCut");
  double jet_EtaCut              = getPreCutValue1("jet_EtaCut");
  double jet_HighEtaCut          = getPreCutValue1("jet_HighEtaCut");
  double jet_ele_DeltaRCut       = getPreCutValue1("jet_ele_DeltaRCut") ;
  double jet_muon_DeltaRCut      = getPreCutValue1("jet_muon_DeltaRCut");
  double jet_hltMatch_DeltaRCut  = getPreCutValue1("jet_hltMatch_DeltaRCut");

  // muon cuts
  double muon_EtaCut             = getPreCutValue1("muon_EtaCut");
  double muon_PtCut              = getPreCutValue1("muon_PtCut");
  double muon_hltMatch_DeltaRCut = getPreCutValue1("muon_hltMatch_DeltaRCut");

  // electron cuts
  double ele_PtCut    	         = getPreCutValue1("ele_PtCut");
  double ele_hltMatch_DeltaRCut  = getPreCutValue1("ele_hltMatch_DeltaRCut");
  float ele2_PtCut               = getPreCutValue1("ele2_PtCut");

  //--------------------------------------------------------------------------
  // Tell the user how many entries we'll look at
  //--------------------------------------------------------------------------

  Long64_t nentries = GetTreeEntries();
  std::cout << "analysisClass::analysisClass(): nentries = " << nentries << std::endl;
  
  //--------------------------------------------------------------------------
  // Create HLT collections in advance (won't need all of them)
  //--------------------------------------------------------------------------
  
  // QCD photon filters
  CollectionPtr c_hltPhoton_QCD_all;

  // muon filter
  CollectionPtr c_hltMuon_SingleMu_all;

  // Signal
  //CollectionPtr c_hltEle45_Signal_all;
  CollectionPtr c_hltPFJet50_Signal_all;
  CollectionPtr c_hltPFJet200_Signal_all;
  CollectionPtr c_hltDoubleEle_Signal_all;

  // trigger
  CollectionPtr c_trigger_l3jets_all;

  // Tag and probe
  CollectionPtr c_hltEle27WP85Gsf_all;

  //--------------------------------------------------------------------------
  // For smearing systematics samples, you'll need a random number generator
  //--------------------------------------------------------------------------

  unsigned int seed = 987654321;
  TRandom3 * rootEngine = new TRandom3 ( seed ) ;

  //--------------------------------------------------------------------------
  // For smearing/scaling systematics samples, you'll need to calculate the effect on MET
  //--------------------------------------------------------------------------
  
  TLorentzVector v_delta_met;

  TLorentzVector v_PFMETRaw;
  TLorentzVector v_PFMETType1Cor;
  //TLorentzVector v_PFMETType1XYCor;
  

  /*//------------------------------------------------------------------
   *
   *
   *      
   *      Start analysis loop!
   *
   *
   *
   *///-----------------------------------------------------------------
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    readerTools_->LoadEntry(jentry);
    //// test
    double event = readerTools_->ReadValueBranch<ULong64_t>("event");
    double ls = readerTools_->ReadValueBranch<UInt_t>("luminosityBlock");
    double run = readerTools_->ReadValueBranch<UInt_t>("run");
    //if(event!=21273 && event!=21288) continue;
    //if(ls!=227) continue;
    // run ls event
    //if(jentry > 1000) continue;
    std::cout << "------" << std::endl << static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(ls) << " " << static_cast<unsigned int>(event) << std::endl << "------" << std::endl;;
    ////std::string current_file_name ( readerTools_->GetTree()->GetCurrentFile()->GetName());
    ////cout << "Found the event! in file:" << current_file_name << endl;
    //if(jentry==0)
    //  std::cout << "WARNING WARNING WARNING -- ONLY RUNNING OVER FIRST 100 ENTRIES!" << std::endl;
    ////test

    //-----------------------------------------------------------------
    // Print progress
    //-----------------------------------------------------------------
    if(jentry < 10 || jentry%5000 == 0) std::cout << "analysisClass:Loop(): jentry = " << jentry << "/" << nentries << std::endl;
    
    //-----------------------------------------------------------------
    // Get access to HLT decisions
    //-----------------------------------------------------------------

    getTriggers(jentry); 
    //printTriggers();
    //printFiredTriggers();

    //-----------------------------------------------------------------
    // Get access to HLT filter objects
    //-----------------------------------------------------------------
    HLTriggerObjectCollectionHelper helper(*this);

    if ( reducedSkimType == 0 ){ 

      // QCD photon triggers
      std::vector<int> typeIds {11, 22};
      c_hltPhoton_QCD_all = helper.GetFilterObjectsByType(typeIds);
      //c_hltPhoton_QCD_all->examine<HLTriggerObject>("c_hltPhoton_QCD_all");

    }

    //else if ( reducedSkimType == 1 || reducedSkimType == 2 || reducedSkimType == 3 || reducedSkimType == 4){

    //  // SingleMu trigger: path was HLT_Mu40_eta2p1_v9+
    //  //c_hltMuon_SingleMu_all       = helper.GetL3FilterObjectsByPath("HLT_Mu45_eta2p1_v"); // will do prefix matching

    //  //// Ele+jets signal triggers
    //  //CollectionPtr trigger_l3objects_all = helper.GetL3FilterObjectsByPath("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v");
    //  //c_trigger_l3jets_all = trigger_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_JET);
    //  ////c_trigger_l3jets_all->examine<HLTriggerObject>("c_trigger_l3jets_all");
    //  //// which one passed the last filter?
    //  ////CollectionPtr trigger_lastObjects_all = helper.GetLastFilterObjectsByPath("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v");

    //  //// electrons seem to come as TRIGGER_PHOTON most of the time
    //  //c_hltEle45_Signal_all = trigger_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_PHOTON);
    //  //// if not, try TRIGGER_ELECTRON
    //  //if(c_hltEle45_Signal_all->GetSize() == 0)
    //  //  c_hltEle45_Signal_all = trigger_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_ELECTRON);
    //  //// Note: could also be TRIGGER_CLUSTER?

    //  // jets
    //  //c_hltPFJet200_Signal_all     =  trigger_lastObjects_all->SkimByID<HLTriggerObject>(TRIGGER_JET);
    //  // get rid of overlaps
    //  // XXX no, keep all L3 jets
    //  //c_hltPFJet50_Signal_all      =  c_trigger_l3jets_all -> SkimByVetoDRMatch<HLTriggerObject,HLTriggerObject>( c_hltPFJet200_Signal_all, 0.3 );

    //  //// DoubleEle signal trigger
    //  //CollectionPtr double_ele_l3objects_all    = helper.GetL3FilterObjectsByPath("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v");
    //  //c_hltDoubleEle_Signal_all    = double_ele_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_PHOTON);
    //  //// if not, try TRIGGER_ELECTRON
    //  //if(c_hltDoubleEle_Signal_all->GetSize() == 0)
    //  //  c_hltDoubleEle_Signal_all = double_ele_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_ELECTRON);

    //  //// Tag and probe trigger
    //  //CollectionPtr tagProbe_l3objects_all;
    //  //if(!isData())
    //  //  tagProbe_l3objects_all  = helper.GetL3FilterObjectsByPath("HLT_Ele27_WPLoose_Gsf_v");
    //  //else
    //  //  tagProbe_l3objects_all  = helper.GetL3FilterObjectsByPath("HLT_Ele27_WPLoose_Gsf_v");
    //  //c_hltEle27WP85Gsf_all = tagProbe_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_PHOTON);
    //  //// if not, try TRIGGER_ELECTRON
    //  //if(c_hltEle27WP85Gsf_all->GetSize() == 0)
    //  //  c_hltEle27WP85Gsf_all = tagProbe_l3objects_all->SkimByID<HLTriggerObject>(TRIGGER_ELECTRON);

    //}


    //-----------------------------------------------------------------
    // Define initial, inclusive collections for physics objects
    //-----------------------------------------------------------------

    CollectionPtr c_gen_all(new Collection(readerTools_));
    if(!isData()) {
      c_gen_all.reset(new Collection(readerTools_, readerTools_->ReadValueBranch<UInt_t>("nGenPart")));
    }
    CollectionPtr c_ele_all   ( new Collection(readerTools_, readerTools_->ReadValueBranch<UInt_t>("nElectron")));
    //c_ele_all->examine<Electron>("c_ele_all = All reco electrons");
    //Electron ele1 = c_ele_all -> GetConstituent<Electron>(0);
    //for(unsigned int i=0; i<10; ++i) { 
    //  std::cout << "cut = " << i << " idLevel = " << ele1.GetNbitFromBitMap(i, 3) << std::endl;
    //}

    CollectionPtr c_muon_all  ( new Collection(readerTools_, readerTools_->ReadValueBranch<UInt_t>("nMuon")));
    CollectionPtr c_genJet_all(new Collection(readerTools_));
    if(!isData()) {
      c_genJet_all.reset(new Collection(readerTools_, readerTools_->ReadValueBranch<UInt_t>("nGenJet")));
    }
    c_genJet_all->examine<GenJet>("c_genJet_all");
    CollectionPtr c_pfjet_all ( new Collection(readerTools_, readerTools_->ReadValueBranch<UInt_t>("nJet")));
    //  New: Cut all jet collections at 17 GeV as per 2016 custom skim
    //c_pfjet_all = c_pfjet_all -> SkimByMinPt   <PFJet>( 17.0 );
    //c_pfjet_all->examine<PFJet>("c_pfjet_all with Pt>17");

    //-----------------------------------------------------------------
    // All skims need GEN particles/jets
    //-----------------------------------------------------------------

    //c_gen_all->examine<GenParticle>("c_gen_all");
    CollectionPtr c_genEle_final = c_gen_all    -> SkimByID<GenParticle>(GEN_ELE_HARDPROCESS_FINALSTATE);
    c_genEle_final->examine<GenParticle>("c_genEle_final = c_gen_all after SkimByID<GenParticle>(GEN_ELE_HARDPROCESS_FINALSTATE)");

    CollectionPtr c_genEle_fromLQ_finalState = c_gen_all -> SkimByID<GenParticle>(GEN_ELE_FROM_LQ);
    c_genEle_fromLQ_finalState = c_genEle_fromLQ_finalState -> SkimByID<GenParticle>(GEN_ELE_HARDPROCESS_FINALSTATE);
    c_genEle_fromLQ_finalState->examine<GenParticle>("c_genEle_fromLQ_finalState = (GEN_ELE_FROM_LQ && GEN_ELE_HARDPROCESS_FINALSTATE)");

    CollectionPtr c_genPartMatchEle = c_gen_all -> SkimByRequireDRMatch<GenParticle, GenParticle>(c_genEle_final, 0.4);
    c_genPartMatchEle->examine<GenParticle>("c_genPartMatchDR04GenEleFinal");

    CollectionPtr c_genMu_final = c_gen_all    -> SkimByID<GenParticle>(GEN_MU_HARD_SCATTER);

    CollectionPtr c_genNu_final = c_gen_all    -> SkimByID<GenParticle>(GEN_NU_HARD_SCATTER);
    //c_genNu_final->examine<GenParticle>("c_genNu_final = c_gen_all after SkimByID<GenParticle>(GEN_NU_HARD_SCATTER)");
    //c_genJet_all->examine<GenJet>("c_genJet_all");
    CollectionPtr c_genJet_final = c_genJet_all;

    //c_genZgamma_final->examine<GenParticle>("c_genZgamma_final = c_gen_all after SkimByID<GenParticle>(GEN_ZGAMMA_HARD_SCATTER)");
    CollectionPtr c_genNuFromW_final   = c_gen_all -> SkimByID<GenParticle>(GEN_NU_FROM_W);
    CollectionPtr c_genTop             = c_gen_all -> SkimByID<GenParticle>(GEN_TOP);
    CollectionPtr c_genTop_final       = c_genTop  -> SkimByID<GenParticle>(GEN_STATUS62);

    CollectionPtr c_genLQ              = c_gen_all ->SkimByID<GenParticle>(GEN_LQ);
    CollectionPtr c_genLQ_final        = c_genLQ  -> SkimByID<GenParticle>(GEN_STATUS62);
    //c_genLQ_final->examine<GenParticle>("c_genLQ = c_gen_all after SkimByID GEN_LQ and GEN_STATUS62");

    CollectionPtr c_genQuark_hardScatter = c_gen_all -> SkimByID<GenParticle>(GEN_QUARK_HARD_SCATTER);
    CollectionPtr c_genQuark_hardProcess = c_gen_all -> SkimByID<GenParticle>(GEN_QUARK_HARD_PROCESS);
    //c_genQuark_hardScatter->examine<GenParticle>("c_genQuark_hardScatter");
    //const CollectionPtr c_genQuark_hardScatter_const(c_genQuark_hardScatter);
    //CollectionPtr c_genJetMatchedLQ = c_genJet_all -> SkimByRequireDRMatch<GenJet, GenParticle>(c_genQuark_hardScatter_const, 0.3);
    //c_genJetMatchedLQ->examine<GenJet>("c_genJetMatchedLQ");
    //
    const CollectionPtr c_genQuark_hardProcess_const(c_genQuark_hardProcess);
    CollectionPtr c_genJetMatchedLQ = c_genJet_all -> SkimByRequireDRMatch<GenJet, GenParticle>(c_genQuark_hardProcess_const, 0.3);
    c_genJetMatchedLQ->examine<GenJet>("c_genJetMatchedLQ(quarkFromHardprocess)");
    CollectionPtr c_genQuark_hardProcess_matchedGenJet = c_genQuark_hardProcess->SkimByRequireDRMatch<GenParticle, GenJet>(c_genJet_all, 0.3);

    //-----------------------------------------------------------------
    // If this is MC, smear jets if requested (not needed for NanoAOD-based setup)
    // Don't do it for data
    //-----------------------------------------------------------------
    if ( !isData() && do_jer) do_jer = true;
    else do_jer = false; // never for data

    //-----------------------------------------------------------------
    // Energy scaling and resolution smearing here
    //-----------------------------------------------------------------

    //SIC: JER/JES already replaced by JER/JES variations from nanoAOD-tools; eventually EER/EES should be there as well
    if ( do_eer || do_jer || do_ees || do_jes ) { 

      // If  you're scaling/smearing PFJets, recall that only jets with pt > 10 GeV affect the PFMET
      // Also, only scale/smear the jets in our eta range (jets in the calorimeter crack are suspect)
      c_pfjet_all = c_pfjet_all -> SkimByEtaRange<PFJet>( -jet_EtaCut, jet_EtaCut );

      // Set the PFMET difference to zero

      v_delta_met.SetPtEtaPhiM(0.,0.,0.,0.);

      // Do the energy scale / energy resolution operations
      // dR for matching = Rcone/2
      //   see: https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution

      if ( do_eer ) c_ele_all      -> MatchAndSmearEnergy <Electron, GenParticle> ( c_genEle_final, 0.4/2.0, rootEngine, v_delta_met );
      if ( do_jer ) c_pfjet_all    -> MatchAndSmearEnergy <PFJet   , GenJet     > ( c_genJet_final, 0.4/2.0, rootEngine, v_delta_met );
      if ( do_ees ) c_ele_all      -> ScaleEnergy <Electron> ( electron_energy_scale_sign, v_delta_met );
      if ( do_jes ) c_pfjet_all    -> ScaleEnergy <PFJet   > ( pfjet_energy_scale_sign   , v_delta_met );

      // Propagate the results to the PFMET

      v_PFMETType1Cor   .SetPtEtaPhiM( readerTools_->ReadValueBranch<Float_t>("MET_pt"), 0., readerTools_->ReadValueBranch<Float_t>("MET_phi"), 0. );
      //v_PFMETType1XYCor.SetPtEtaPhiM( (*PFMETType1XYCor)[0] , 0., (*PFMETPhiType1XYCor)[0] , 0. );

      v_PFMETType1Cor    = v_PFMETType1Cor    + v_delta_met;
      //v_PFMETType1XYCor = v_PFMETType1XYCor + v_delta_met;

      //FIXME
      //(*PFMETType1Cor      )[0] = v_PFMETType1Cor   .Pt();
      //(*PFMETType1XYCor   )[0] = v_PFMETType1XYCor.Pt();
      //
      //FIXME
      //(*PFMETPhiType1Cor   )[0] = v_PFMETType1Cor   .Phi();
      //(*PFMETPhiType1XYCor)[0] = v_PFMETType1XYCor.Phi();
    }

    // new systematics handling
    std::vector<Electron> smearedEles;
    std::vector<Electron> scaledUpEles;
    std::vector<Electron> scaledDownEles;
    if(!isData()) {
      smearedEles = c_ele_all->MatchAndSmearEnergy <Electron, GenParticle> ( c_genEle_final, 0.4/2.0, rootEngine, v_delta_met );
      scaledUpEles = c_ele_all->ScaleEnergy <Electron> (1 , v_delta_met );
      scaledDownEles = c_ele_all->ScaleEnergy <Electron> (-1, v_delta_met );
    }

    //-----------------------------------------------------------------
    // QCD skims    (reducedSkimType = 0     ) have loose electrons
    // Signal skims (reducedSkimType = 1 - 4 ) have HEEP  electrons
    //-----------------------------------------------------------------

    CollectionPtr c_ele_HEEP;
    CollectionPtr c_ele_egammaLoose;
    CollectionPtr c_ele_egammaLoose_ptCut;
    CollectionPtr c_ele_final;
    CollectionPtr c_ele_final_ptCut;
    CollectionPtr c_ele_vLoose_ptCut;
    ID heepIdType = HEEP70;
    if(analysisYear == 2018)
      heepIdType = HEEP70_2018;

    if ( reducedSkimType == 0 ){ 
      CollectionPtr c_ele_loose = c_ele_all   -> SkimByID<LooseElectron>( FAKE_RATE_HEEP_LOOSE);
      c_ele_final               = c_ele_loose;
      c_ele_final_ptCut         = c_ele_final -> SkimByMinPt<LooseElectron>( ele_PtCut  );
      CollectionPtr c_ele_vLoose = c_ele_all   -> SkimByID<LooseElectron>(FAKE_RATE_VERY_LOOSE);
      c_ele_vLoose_ptCut         = c_ele_vLoose -> SkimByMinPt<LooseElectron>( 10.0 );
      //c_ele_all->examine<Electron>("c_ele_all");
      //c_ele_loose->examine<Electron>("c_ele_loose");
      //c_ele_final_ptCut->examine<Electron>("c_ele_final_ptCut");
    }

    else if ( reducedSkimType == 1 || reducedSkimType == 2 || reducedSkimType == 3 || reducedSkimType == 4 ){
      c_ele_HEEP  = c_ele_all -> SkimByID <Electron> ( heepIdType );
      //CollectionPtr c_ele_HEEP  = c_ele_all -> SkimByID <Electron> ( HEEP70_MANUAL , true );
      //c_ele_final               = c_ele_HEEP;
      c_ele_egammaLoose = c_ele_all -> SkimByID<Electron>(EGAMMA_BUILTIN_LOOSE);
      //c_ele_final               = c_ele_egammaLoose;
      c_ele_final = c_ele_all;
      c_ele_final_ptCut         = c_ele_final -> SkimByMinPt<Electron>( ele_PtCut  );
      CollectionPtr c_ele_vLoose = c_ele_all   -> SkimByID<LooseElectron>(FAKE_RATE_VERY_LOOSE);
      c_ele_vLoose_ptCut         = c_ele_vLoose -> SkimByMinPt<LooseElectron>( 10.0 );
      c_ele_egammaLoose_ptCut    = c_ele_egammaLoose -> SkimByMinPt<Electron>( 35  );
    }
    // look at final electrons
    //Electron ele1_tmp;
    //if(c_ele_final->GetSize() > 0)
    //{
    //  ele1_tmp = c_ele_final -> GetConstituent<Electron>(0);
    //  if(ele1_tmp.Pt() > 2300)
    //  {
    //    std::cout << static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(luminosityBlock) << " " << static_cast<ULong64_t>(event) << std::endl;
    //    c_ele_all->examine<Electron>("c_ele_all");
    //    CollectionPtr c_ele_HEEP  = c_ele_all -> SkimByID <Electron> ( HEEP70_MANUAL , true );
    //    c_ele_HEEP->examine<Electron>("c_ele_HEEP_manual");
    //    c_ele_final->examine<Electron>("c_ele_final");
    //  }
    //}
    ////c_ele_final_ptCut->examine<Electron>("c_ele_final_ptCut");
    CollectionPtr c_genJetVetoGenElectron = c_genJet_all -> SkimByVetoDRMatch<GenJet, GenParticle>(c_genEle_final, 0.1);
    //CollectionPtr c_genJetsNearElectron = c_genJetVetoGenElectron -> SkimByRequireDRMatch<GenJet, Electron>(c_ele_all, 0.3);
    //c_genJetsNearElectron->examine<GenJet>("c_genJetsNearElectron");
    CollectionPtr c_genJetsNearEGLElectron = c_genJetVetoGenElectron -> SkimByRequireDRMatch<GenJet, Electron>(c_ele_egammaLoose, 0.3);
    c_genJetsNearEGLElectron->examine<GenJet>("c_genJetsNearEGLooseElectron");
    for(int i=0; i<c_genJetsNearEGLElectron->GetSize(); ++i) {
      GenJet jet = c_genJetsNearEGLElectron -> GetConstituent<GenJet>(i);
      bool lqMatchedJet = c_genJetMatchedLQ->Has(jet);
      if(lqMatchedJet)
        FillUserTH1D("genJetsNearEGLooseEle_matchedToLQ", 1);
      else
        FillUserTH1D("genJetsNearEGLooseEle_matchedToLQ", 0);
      for(int j=0; j<c_ele_egammaLoose->GetSize(); ++j) {
        Electron ele = c_ele_egammaLoose->GetConstituent<Electron>(j);
        float dr = ele.DeltaR(&jet);
        if(dr < 0.3) {
          if(lqMatchedJet) {
            if(ele.PassHEEPID())
              FillUserTH1D("genJetsNearEGLooseEle_matchedToLQ_elePassHeep", 1);
            else
              FillUserTH1D("genJetsNearEGLooseEle_matchedToLQ_elePassHeep", 0);
          }
          else {
            if(ele.PassHEEPID())
              FillUserTH1D("genJetsNearEGLooseEle_unmatchedToLQ_elePassHeep", 1);
            else
              FillUserTH1D("genJetsNearEGLooseEle_unmatchedToLQ_elePassHeep", 0);
          }
        }
      }
    }
    FillUserTH1D("nGenJetsNearEGLooseEle",c_genJetsNearEGLElectron->GetSize());

    c_ele_all->examine<Electron>("c_ele_all");
    FillUserTH1D("nEleReco",c_ele_all->GetSize());
    FillUserTH1D("nEleReco_ptCut",c_ele_final_ptCut->GetSize());
    FillUserTH1D("nEleHEEP",c_ele_HEEP->GetSize());
    int nLooseElePassPtPassEta = 0;
    for(int i=0; i<c_ele_egammaLoose_ptCut->GetSize(); ++i) {
      Electron ele = c_ele_egammaLoose_ptCut -> GetConstituent<Electron>(i);
      if(ele.PassHEEPMinPtCut()  && ele.PassHEEPGsfEleSCEtaMultiRangeCut())
        nLooseElePassPtPassEta++;
    }
    FillUserTH1D("nEleEGLoosePt35EtaCut", nLooseElePassPtPassEta);
    FillUserTH1D("nEleEGLoosePt35",c_ele_egammaLoose_ptCut->GetSize());
    FillUserTH1D("nEleEGLoose",c_ele_egammaLoose->GetSize());
    FillUserTH1D("nEleFinal",c_ele_final->GetSize());
    FillUserTH2D("nEleNTupleVsNeleRsk",c_ele_final->GetSize(),c_ele_all->GetSize());

    for(int i=0; i<c_ele_final->GetSize(); ++i) {
      Electron ele = c_ele_final -> GetConstituent<Electron>(i);
      std::string electronNumber = std::to_string(i+1);
      if(!ele.PassHEEPMinPtCut() || !ele.PassHEEPGsfEleSCEtaMultiRangeCut())
        continue; // only consider electrons which pass pT/acceptance cuts
      //if(ele.SCEt() <= 150) continue;
      //if(ele.SCEt() > ele2_PtCut) continue;
      //if(fabs(ele.SCEta()) >= 1.4442) continue;
      bool passHeep = ele.PassHEEPID();
      bool passLoose = ele.PassEGammaIDLoose();
      FillUserTH1D("SCEtaElectron", ele.SCEta());
      FillUserTH1D("SCEtElectron", ele.SCEt());
      FillUserTH1D("SCEnergyElectron", ele.SCEnergy());
      if(isElectronBarrel(ele.SCEta())) {
          FillUserTH1D("SCEtBarrel", ele.SCEt());
          FillUserTH1D("SCEnergyBarrel", ele.SCEnergy());
          FillUserTH1D("TrkIsoBarrel", ele.HEEP70TrackIsolation());
          FillUserTH1D("EMHadD1IsoBarrel", ele.HEEPCaloIsolation());
      }
      else if(ele.PassHEEPGsfEleSCEtaMultiRangeCut()) {
          FillUserTH1D("SCEtEndcap", ele.SCEt());
          FillUserTH1D("SCEnergyEndcap", ele.SCEnergy());
          FillUserTH1D("TrkIsoEndcap", ele.HEEP70TrackIsolation());
          FillUserTH1D("EMHadD1IsoEndcap", ele.HEEPCaloIsolation());
      }
      if(passHeep) {
        FillUserTH1D("SCEtaElectronPassHEEP", ele.SCEta());
        FillUserTH1D("SCEtElectronPassHEEP", ele.SCEt());
        FillUserTH1D("SCEnergyElectronPassHEEP", ele.SCEnergy());
        if(isElectronBarrel(ele.SCEta())) {
          FillUserTH1D("EMHadD1IsoBarrelPassHEEP", ele.HEEPCaloIsolation());
          FillUserTH1D("TrkIsoBarrelPassHEEP", ele.HEEP70TrackIsolation());
          FillUserTH1D("SCEtBarrelPassHEEP", ele.SCEt());
          FillUserTH1D("SCEnergyBarrelPassHEEP", ele.SCEnergy());
        }
        else {
          FillUserTH1D("EMHadD1IsoEndcapPassHEEP", ele.HEEPCaloIsolation());
          FillUserTH1D("TrkIsoEndcapPassHEEP", ele.HEEP70TrackIsolation());
          FillUserTH1D("SCEtEndcapPassHEEP", ele.SCEt());
          FillUserTH1D("SCEnergyEndcapPassHEEP", ele.SCEnergy());
        }
      }
      else {
        FillUserTH1D("SCEtaElectronFailHEEP", ele.SCEta());
        FillUserTH1D("SCEtElectronFailHEEP", ele.SCEt());
        FillUserTH1D("SCEnergyElectronFailHEEP", ele.SCEnergy());
        if(isElectronBarrel(ele.SCEta())) {
          FillUserTH1D("EMHadD1IsoBarrelFailHEEP", ele.HEEPCaloIsolation());
          FillUserTH1D("TrkIsoBarrelFailHEEP", ele.HEEP70TrackIsolation());
          FillUserTH1D("SCEtBarrelFailHEEP", ele.SCEt());
          FillUserTH1D("SCEnergyBarrelFailHEEP", ele.SCEnergy());
        }
        else if(ele.PassHEEPGsfEleSCEtaMultiRangeCut()) {
          FillUserTH1D("EMHadD1IsoEndcapFailHEEP", ele.HEEPCaloIsolation());
          FillUserTH1D("TrkIsoEndcapFailHEEP", ele.HEEP70TrackIsolation());
          FillUserTH1D("SCEnergyEndcapFailHEEP", ele.SCEnergy());
          FillUserTH1D("SCEtEndcapFailHEEP", ele.SCEt());
        }
      }
      if(passLoose) {
        FillUserTH1D("SCEtElectronPassLoose", ele.SCEt());
        if(isElectronBarrel(ele.SCEta())) {
          FillUserTH1D("SCEtBarrelPassLoose", ele.SCEt());
        }
        else
          FillUserTH1D("SCEtEndcapPassLoose", ele.SCEt());
      }
      else {
        FillUserTH1D("SCEtElectronFailLoose", ele.SCEt());
        if(isElectronBarrel(ele.SCEta())) {
          FillUserTH1D("SCEtBarrelFailLoose", ele.SCEt());
        }
        else
          FillUserTH1D("SCEtEndcapFailLoose", ele.SCEt());
      }
      if(i > 1)
        break;
      FillUserTH1D(("EMHadD1IsoElectron"+electronNumber).c_str(), ele.HEEPCaloIsolation());
      FillUserTH1D(("TrkIsoElectron"+electronNumber).c_str(), ele.HEEP70TrackIsolation());
      FillUserTH1D(("ECALDrivenElectron"+electronNumber).c_str(), ele.PassHEEPEcalDrivenCut());
      FillUserTH1D(("SCEtaElectron"+electronNumber).c_str(), ele.SCEta());
      FillUserTH1D(("SCEtElectron"+electronNumber).c_str(), ele.SCEt());
      FillUserTH1D(("SCEnergyElectron"+electronNumber).c_str(), ele.SCEnergy());
      if(passHeep) {
        FillUserTH1D(("ele"+electronNumber+"PassHeep").c_str(), 1);
        FillUserTH1D(("EMHadD1IsoElectron"+electronNumber+"PassHEEP").c_str(), ele.HEEPCaloIsolation());
        FillUserTH1D(("TrkIsoElectron"+electronNumber+"PassHEEP").c_str(), ele.HEEP70TrackIsolation());
        FillUserTH1D(("ECALDrivenElectron"+electronNumber+"PassHEEP").c_str(), ele.PassHEEPEcalDrivenCut());
        FillUserTH1D(("SCEtaElectron"+electronNumber+"PassHEEP").c_str(), ele.SCEta());
        FillUserTH1D(("SCEtElectron"+electronNumber+"PassHEEP").c_str(), ele.SCEt());
      }
      else {
        FillUserTH1D(("ele"+electronNumber+"PassHeep").c_str(), 0);
        FillUserTH1D(("EMHadD1IsoElectron"+electronNumber+"FailHEEP").c_str(), ele.HEEPCaloIsolation());
        FillUserTH1D(("TrkIsoElectron"+electronNumber+"FailHEEP").c_str(), ele.HEEP70TrackIsolation());
        FillUserTH1D(("ECALDrivenElectron"+electronNumber+"FailHEEP").c_str(), ele.PassHEEPEcalDrivenCut());
        FillUserTH1D(("SCEtaElectron"+electronNumber+"FailHEEP").c_str(), ele.SCEta());
        FillUserTH1D(("SCEtElectron"+electronNumber+"FailHEEP").c_str(), ele.SCEt());
      }
      if(passLoose) {
        FillUserTH1D(("SCEtElectron"+electronNumber+"PassLoose").c_str(), ele.SCEt());
      }
      else {
        FillUserTH1D(("SCEtElectron"+electronNumber+"FailLoose").c_str(), ele.SCEt());
      }
    }

    if( c_ele_all->GetSize() >= 1) {
      Electron ele1 = c_ele_all -> GetConstituent<Electron>(0);
      bool ele1PassHeep = ele1.PassHEEPID();
      if ( c_ele_all->GetSize() >= 2 ) {
        Electron ele2 = c_ele_all -> GetConstituent<Electron>(1);
        bool ele2PassHeep = ele2.PassHEEPID();
        if(ele1PassHeep && ele2PassHeep)
          FillUserTH1D("heepElectronsCategory", 1);
        if ( c_ele_all->GetSize() >= 3 ) {
          Electron ele3 = c_ele_all -> GetConstituent<Electron>(2);
          bool ele3PassHeep = ele3.PassHEEPID();
          if(ele2PassHeep && ele3PassHeep)
            FillUserTH1D("heepElectronsCategory", 2);
          if(ele1PassHeep && ele3PassHeep)
            FillUserTH1D("heepElectronsCategory", 3);
        }
      }
    }

    //-----------------------------------------------------------------
    // All skims need muons
    //-----------------------------------------------------------------

    CollectionPtr c_muon_eta               = c_muon_all       -> SkimByEtaRange<Muon> ( -muon_EtaCut, muon_EtaCut );
    //CollectionPtr c_muon_eta_IDTight       = c_muon_eta       -> SkimByID      <Muon> ( MUON_TIGHT_PFISO04TIGHT );
    CollectionPtr c_muon_eta_IDHighPt      = c_muon_eta       -> SkimByID      <Muon> ( MUON_HIGH_PT_TRKRELISO03);
    CollectionPtr c_muon_eta_IDLoose       = c_muon_eta       -> SkimByID      <Muon> ( MUON_LOOSE);
    CollectionPtr c_muon_final             = c_muon_eta_IDHighPt;
    //CollectionPtr c_muon_final             = c_muon_eta_IDLoose;
    CollectionPtr c_muon_final_ptCut       = c_muon_final     -> SkimByMinPt   <Muon> ( muon_PtCut );
    //c_muon_all->examine<Muon>("c_muon_all");
    //c_muon_final->examine<Muon>("c_muon_final");
    //c_muon_eta_IDLoose->examine<Muon>("c_muon_eta_IDLoose");

    //-----------------------------------------------------------------
    // All skims need PFJets
    //-----------------------------------------------------------------

    //c_pfjet_all->examine<PFJet>("c_pfjet_all");
    CollectionPtr c_pfjet_central                     = c_pfjet_all                          -> SkimByEtaRange   <PFJet>          ( -jet_EtaCut, jet_EtaCut   );
    //c_pfjet_central->examine<PFJet>("c_pfjet_central");
    CollectionPtr c_pfjet_central_ID                  = c_pfjet_central                      -> SkimByID         <PFJet>          ( PFJET_TIGHT );    
    c_pfjet_central_ID->examine<PFJet>("c_pfjet_central_ID");
    //CollectionPtr c_pfjet_central_ID_noMuonOverlap    = c_pfjet_central_ID                   -> SkimByVetoDRMatch<PFJet, Muon>    ( c_muon_final_ptCut   , jet_muon_DeltaRCut  );
    CollectionPtr c_pfjet_central_ID_noMuonOverlap    = c_pfjet_central_ID                   -> SkimByVetoDRMatch<PFJet, Muon>    ( c_muon_final   , jet_muon_DeltaRCut  );
    //c_pfjet_central_ID_noMuonOverlap->examine<PFJet>("c_pfjet_central_ID_noMuonOverlap [after muon cleaning]");
    //CollectionPtr c_pfjet_central_ID_noLeptonOverlap  = c_pfjet_central_ID_noMuonOverlap     -> SkimByVetoDRMatch<PFJet, Electron>( c_ele_final_ptCut    , jet_ele_DeltaRCut );
    //CollectionPtr c_pfjet_central_ID_noLeptonOverlap  = c_pfjet_central_ID_noMuonOverlap     -> SkimByVetoDRMatch<PFJet, Electron>( c_ele_final    , jet_ele_DeltaRCut );
    CollectionPtr c_pfjet_central_ID_noLeptonOverlap  = c_pfjet_central_ID_noMuonOverlap     -> SkimByVetoDRMatch<PFJet, Electron>( c_ele_egammaLoose    , jet_ele_DeltaRCut );
    //c_pfjet_central_ID_noLeptonOverlap->examine<PFJet>("c_pfjet_central_ID_noLeptonOverlap [after electron cleaning]");
    //CollectionPtr c_pfjet_final                       = c_pfjet_central_ID_noLeptonOverlap;
    CollectionPtr c_pfjet_final = c_pfjet_central_ID;
    c_pfjet_final->examine<PFJet>("c_pfjet_final");
    CollectionPtr c_pfjet_final_ptCut                 = c_pfjet_final                        -> SkimByMinPt      <PFJet>          ( jet_PtCut );
    c_pfjet_final_ptCut->examine<PFJet>("c_pfjet_final_ptCut");
    //if(c_pfjet_final->GetSize() > 4) {
    //  // run ls event
    //  //double run = readerTools_->ReadValueBranch<UInt_t>("run");
    //  double event = readerTools_->ReadValueBranch<ULong64_t>("event");
    //  //double ls = readerTools_->ReadValueBranch<UInt_t>("luminosityBlock");
    //  std::cout << static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(ls) << " " << static_cast<unsigned int>(event) << std::endl;
    //  c_pfjet_all->examine<PFJet>("c_pfjet_all");
    //  c_pfjet_final->examine<PFJet>("c_pfjet_final");
    //  c_ele_all->examine<Electron>("c_ele_all");
    //  c_ele_final->examine<Electron>("c_ele_final");
    //  c_muon_all->examine<Muon>("c_muon_all");
    //  c_muon_final->examine<Muon>("c_muon_final");
    //}
    //PFJet jet2 = c_pfjet_final -> GetConstituent<PFJet>(1);
    //if(jet2.Pt() < 50) {
    //  std::cout << "----> Interesting event!" << std::endl;
    //  // run ls event
    //  double run = readerTools_->ReadValueBranch<UInt_t>("run");
    //  double event = readerTools_->ReadValueBranch<ULong64_t>("event");
    //  double ls = readerTools_->ReadValueBranch<UInt_t>("luminosityBlock");
    //  std::cout << static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(ls) << " " << static_cast<unsigned int>(event) << std::endl;
    //}
    //const CollectionPtr c_pfjet_final_const(c_pfjet_final);
    //CollectionPtr c_pfJetMatchedLQ = c_genJetMatchedLQ -> SkimByRequireDRMatch<GenJet, PFJet>(c_pfjet_final_const, 0.3);
    CollectionPtr c_pfjet_central_ID_eleOverlap  = c_pfjet_central_ID_noMuonOverlap     -> SkimByRequireDRMatch<PFJet, Electron>( c_ele_egammaLoose    , jet_ele_DeltaRCut );
    c_pfjet_central_ID_eleOverlap->examine<PFJet>("c_pfjet_central_ID_eleOverlap (PFJets near loose ele)");
    CollectionPtr c_pfjet_central_ID_genEleLQOverlap = c_pfjet_central_ID_noMuonOverlap->SkimByRequireDRMatch<PFJet, GenParticle>( c_genEle_fromLQ_finalState, jet_ele_DeltaRCut);
    c_pfjet_central_ID_genEleLQOverlap->examine<PFJet>("c_pfjet_central_ID_genEleLQOverlap (PFJets near gen ele from LQ)");
    for(int i=0; i<c_pfjet_central_ID_genEleLQOverlap->GetSize(); ++i) {
      PFJet jet = c_pfjet_central_ID_genEleLQOverlap -> GetConstituent<PFJet>(i);
      for(int j=0; j<c_ele_egammaLoose->GetSize(); ++j) {
        Electron ele = c_ele_egammaLoose->GetConstituent<Electron>(j);
        float dr = ele.DeltaR(&jet);
        float deltaPt = jet.Pt()-ele.Pt();
        if(dr < 0.1) {
          if(!ele.PassHEEPID()) {
            std::cout << "Interesting! PFJet matched with loose ele, but ele failed HEEP. deltaR = " << dr << ", deltaPt = " << deltaPt << std::endl;
            std::cout << jet << std::endl;
            std::cout << ele << std::endl;
            FillUserTH1D("DeltaRMatchJetEle_eleFailHeep", dr);
            FillUserTH1D("DeltaPtMatchJetEle_eleFailHeep", deltaPt);
            FillUserTH2D("DeltaRVsDeltaPtMatchJetEle_eleFailHeep", deltaPt, dr);
            if(!ele.PassHEEPGsfEleEmHadD1IsoRhoCut()) {
              FillUserTH1D("DeltaRMatchJetEle_eleFailHeepCaloIso", dr);
              FillUserTH1D("DeltaPtMatchJetEle_eleFailHeepCaloIso", deltaPt);
              FillUserTH2D("DeltaRVsDeltaPtMatchJetEle_eleFailHeepCaloIso", deltaPt, dr);
            }
          }
          else {
            FillUserTH1D("DeltaRMatchJetEle_elePassHeep", dr);
            FillUserTH1D("DeltaPtMatchJetEle_elePassHeep", deltaPt);
            FillUserTH2D("DeltaRVsDeltaPtMatchJetEle_elePassHeep", deltaPt, dr);
          }
        }
      }
    }

    CollectionPtr c_pfJetMatchedLQ = c_pfjet_final -> SkimByRequireDRMatch<PFJet, GenJet>(c_genJetMatchedLQ, 0.3);
    c_pfJetMatchedLQ->examine<PFJet>("c_pfJetMatchedLQ");

    FillUserTH1D("nJetFinalMatchedLQ",c_pfJetMatchedLQ->GetSize());
    FillUserTH1D("nJetFinal",c_pfjet_final->GetSize());
    for(int i=0; i<c_pfjet_final->GetSize(); ++i) {
      if(i > 5) continue;
      PFJet jet = c_pfjet_final -> GetConstituent<PFJet>(i);
      FillUserTH1D(("JetFinalPtJet"+std::to_string(i+1)).c_str(), jet.Pt());
    }
    for(int i=0; i<c_pfJetMatchedLQ->GetSize(); ++i) {
      PFJet jet = c_pfJetMatchedLQ -> GetConstituent<PFJet>(i);
      int idx = jet.GetRawIndex();
      if(idx >= 0)
        FillUserTH1D("idxJetFinalMatchedLQ", idx);
    }
    //-----------------------------------------------------------------
    // We need high-eta jets in order to look at boson recoil
    //-----------------------------------------------------------------

    CollectionPtr c_pfjet_highEta                    = c_pfjet_all                          -> SkimByEtaRange   <PFJet>          ( -jet_HighEtaCut, jet_HighEtaCut   );
    CollectionPtr c_pfjet_highEta_ID                 = c_pfjet_highEta                      -> SkimByID         <PFJet>          ( PFJET_TIGHT );    
    CollectionPtr c_pfjet_highEta_ID_noMuonOverlap   = c_pfjet_highEta_ID                   -> SkimByVetoDRMatch<PFJet, Muon>    ( c_muon_final_ptCut   , jet_muon_DeltaRCut  );
    CollectionPtr c_pfjet_highEta_ID_noLeptonOverlap = c_pfjet_highEta_ID_noMuonOverlap     -> SkimByVetoDRMatch<PFJet, Electron>( c_ele_final_ptCut    , jet_ele_DeltaRCut );
    CollectionPtr c_pfjet_highEta_final              = c_pfjet_highEta_ID_noLeptonOverlap;
    CollectionPtr c_pfjet_highEta_final_ptCut        = c_pfjet_highEta_final                -> SkimByMinPt      <PFJet>          ( jet_PtCut );

    //-----------------------------------------------------------------
    // Get ready to fill variables 
    //-----------------------------------------------------------------

    resetCuts();

    //-----------------------------------------------------------------
    // Fill your single-object variables with values
    //-----------------------------------------------------------------

    fillVariableWithValue( "isData"   , isData()   );
    //fillVariableWithValue( "bunch"    , bunch      );
    fillVariableWithValue( "event"    , readerTools_->ReadValueBranch<ULong64_t>("event")      );
    fillVariableWithValue( "ls"       , readerTools_->ReadValueBranch<UInt_t>("luminosityBlock")         );
    //fillVariableWithValue( "orbit"    , orbit      );
    fillVariableWithValue( "run"      , readerTools_->ReadValueBranch<UInt_t>("run")        );
    //fillVariableWithValue( "ProcessID", ProcessID  );
    //fillVariableWithValue( "PtHat"    , PtHat      );
    float genWeight = -1.0;
    if(!isData()) {
      genWeight = readerTools_->ReadValueBranch<Float_t>("genWeight");
    }
    fillVariableWithValue( "Weight"   , genWeight   );
    //FIXME -- topPtWeights -- perhaps not needed since unused for 2016 analysis
    //fillVariableWithValue( "TopPtWeight",GenParticleTopPtWeight);
    // pileup
    float puWeight = -1.0;
    float puWeightDn = -1.0;
    float puWeightUp = -1.0;
    if(!isData()) {
      puWeight = readerTools_->ReadValueBranch<Float_t>("puWeight");
      puWeightUp = readerTools_->ReadValueBranch<Float_t>("puWeightUp");
      puWeightDn = readerTools_->ReadValueBranch<Float_t>("puWeightDown");
      if(puWeight==0)
        std::cout << "Got puWeight = " << puWeight << "; run: " << getVariableValue("run") << " ls: " << getVariableValue("ls") << " event: " << getVariableValue("event") << std::endl;
    }
    fillVariableWithValue( "puWeight"   , puWeight   );
    fillVariableWithValue( "puWeight_Up"   , puWeightUp   );
    fillVariableWithValue( "puWeight_Dn"   , puWeightDn   );
    // L1 prefiring
    float prefireDefault = -1.0;
    if(hasBranch("L1PreFiringWeight_Nom"))
      fillVariableWithValue("PrefireWeight",readerTools_->ReadValueBranch<Float_t>("L1PreFiringWeight_Nom"));
    else
      fillVariableWithValue("PrefireWeight", prefireDefault);
    if(hasBranch("L1PreFiringWeight_Dn"))
      fillVariableWithValue("PrefireWeight_Dn",readerTools_->ReadValueBranch<Float_t>("L1PreFiringWeight_Dn"));
    else
      fillVariableWithValue("PrefireWeight_Dn", prefireDefault);
    if(hasBranch("L1PreFiringWeight_Up"))
      fillVariableWithValue("PrefireWeight_Up",readerTools_->ReadValueBranch<Float_t>("L1PreFiringWeight_Up"));
    else
      fillVariableWithValue("PrefireWeight_Up", prefireDefault);
    // back to old weights as test
    //fillVariableWithValue("PrefireWeight",readerTools_->ReadValueBranch<Float_t>("L1PreFiringWeight_Nom"));
    //fillVariableWithValue("PrefireWeight_Dn",readerTools_->ReadValueBranch<Float_t>("L1PreFiringWeight_Dn"));
    //fillVariableWithValue("PrefireWeight_Up",readerTools_->ReadValueBranch<Float_t>("L1PreFiringWeight_Up"));


    //-----------------------------------------------------------------
    // Pass JSON
    //-----------------------------------------------------------------
    fillVariableWithValue("PassJSON"                   , passJSON(getVariableValue("run"), getVariableValue("ls"), isData())                       );    

    //-----------------------------------------------------------------
    // Fill MET filter values
    // https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    //-----------------------------------------------------------------
    fillVariableWithValue("PassGlobalSuperTightHalo2016Filter" , int(readerTools_->ReadValueBranch<Bool_t>("Flag_globalSuperTightHalo2016Filter")          == 1));
    fillVariableWithValue("PassGoodVertices"              , int(readerTools_->ReadValueBranch<Bool_t>("Flag_goodVertices")                       == 1));
    fillVariableWithValue("PassHBHENoiseFilter"           , int(readerTools_->ReadValueBranch<Bool_t>("Flag_HBHENoiseFilter")                    == 1));
    fillVariableWithValue("PassHBHENoiseIsoFilter"        , int(readerTools_->ReadValueBranch<Bool_t>("Flag_HBHENoiseIsoFilter")                 == 1));
    fillVariableWithValue("PassBadEESupercrystalFilter"   , int(readerTools_->ReadValueBranch<Bool_t>("Flag_eeBadScFilter")                      == 1));
    fillVariableWithValue("PassEcalDeadCellTrigPrim"      , int(readerTools_->ReadValueBranch<Bool_t>("Flag_EcalDeadCellTriggerPrimitiveFilter") == 1));
    std::string branchName = "Flag_BadChargedCandidateFilter";
    //std::string branchType = std::string(readerTools_->GetTree()->GetBranch(branchName.c_str())->GetLeaf(branchName.c_str())->GetTypeName());
    ////std::cout << "Found branchType=" << branchType << std::endl;
    //if(branchType=="Bool_t") {
    //  fillVariableWithValue("PassChargedCandidateFilter"    , int(readerTools_->ReadValueBranch<Bool_t>(branchName)          == 1));
    //  fillVariableWithValue("PassBadPFMuonFilter"           , int(readerTools_->ReadValueBranch<Bool_t>("Flag_BadPFMuonFilter")                    == 1));
    //}
    //else {
    //  fillVariableWithValue("PassChargedCandidateFilter"    , int(readerTools_->ReadValueBranch<UChar_t>(branchName)          == 1));
    //  fillVariableWithValue("PassBadPFMuonFilter"           , int(readerTools_->ReadValueBranch<UChar_t>("Flag_BadPFMuonFilter")                    == 1));
    //}
    fillVariableWithValue("PassChargedCandidateFilter"    , int(readerTools_->ReadValueBranch<Bool_t>(branchName)          == 1));
    fillVariableWithValue("PassBadPFMuonFilter"           , int(readerTools_->ReadValueBranch<Bool_t>("Flag_BadPFMuonFilter")                    == 1));
    // for 2017 and 2018 only
    if(hasBranch("Flag_ecalBadCalibFilterV2"))
      fillVariableWithValue("PassEcalBadCalibV2Filter"    , int(readerTools_->ReadValueBranch<Bool_t>(branchName)          == 1));
    else
      fillVariableWithValue("PassEcalBadCalibV2Filter"    , 1);

    //-----------------------------------------------------------------
    // Fill MET values
    //-----------------------------------------------------------------

    fillVariableWithValue("PFMET_Type1_Pt"     , readerTools_->ReadValueBranch<Float_t>("MET_pt"));      
    fillVariableWithValue("PFMET_Type1_Phi"    , readerTools_->ReadValueBranch<Float_t>("MET_phi"));
    //fillVariableWithValue("PFMET_Type1XY_Pt"   , PFMETType1XYCor    -> at (0));      
    //fillVariableWithValue("PFMET_Type1XY_Phi"  , PFMETPhiType1XYCor -> at (0));

    if ( !isData() ) { 
      if ( reducedSkimType != 0 ){ 
        fillVariableWithValue("GenMET_Pt"		, readerTools_->ReadValueBranch<Float_t>("GenMET_pt"));
        fillVariableWithValue("GenMET_Phi"	, readerTools_->ReadValueBranch<Float_t>("GenMET_phi"));
        // add LHE variables if needed
        if(hasBranch("LHE_Vpt")) {
          fillVariableWithValue("LHE_Vpt"	  , readerTools_->ReadValueBranch<Float_t>("LHE_Vpt"));
          fillVariableWithValue("LHE_NpLO"	, readerTools_->ReadValueBranch<UChar_t>("LHE_NpLO"));
          fillVariableWithValue("LHE_NpNLO"	, readerTools_->ReadValueBranch<UChar_t>("LHE_NpNLO"));
          fillVariableWithValue("LHE_Njets"	, readerTools_->ReadValueBranch<UChar_t>("LHE_Njets"));
          fillVariableWithValue("LHE_Nglu"	, readerTools_->ReadValueBranch<UChar_t>("LHE_Nglu"));
        }
        if(hasBranch("LHEPdfWeight"))
          fillArrayVariableWithValue("LHEPdfWeight"	, readerTools_->ReadArrayBranch<Float_t>("LHEPdfWeight"));
        if(hasBranch("LHEScaleWeight"))
          fillArrayVariableWithValue("LHEScaleWeight"	, readerTools_->ReadArrayBranch<Float_t>("LHEScaleWeight"));
      }
    }

    //-----------------------------------------------------------------
    // Fill pileup variables
    //-----------------------------------------------------------------

    fillVariableWithValue( "nVertex", readerTools_->ReadValueBranch<Int_t>("PV_npvs"));
    float puNTrueInt = -1.0;
    if(!isData()) {
      puNTrueInt = readerTools_->ReadValueBranch<Float_t>("Pileup_nTrueInt");
    }
    fillVariableWithValue( "nPileUpInt_True", puNTrueInt);

    //fillVariableWithValue( "nPileUpInt_BXminus1", -1 );
    //fillVariableWithValue( "nPileUpInt_BX0"     , -1 );
    //fillVariableWithValue( "nPileUpInt_BXplus1" , -1 );
    //FIXME
    //if ( !isData() ){
    //  for(int pu=0; pu<PileUpInteractions->size(); pu++) {
    //    if(PileUpOriginBX->at(pu) == 0  ) { 
    //      fillVariableWithValue( "nPileUpInt_BX0" , PileUpInteractions    ->at(pu));
    //      fillVariableWithValue( "nPileUpInt_True", PileUpInteractionsTrue->at(pu));
    //    }
    //    if(PileUpOriginBX->at(pu) == -1 ) fillVariableWithValue( "nPileUpInt_BXminus1", PileUpInteractions->at(pu));
    //    if(PileUpOriginBX->at(pu) == 1  ) fillVariableWithValue( "nPileUpInt_BXplus1" , PileUpInteractions->at(pu));
    //  }
    //}

    //-----------------------------------------------------------------
    // How many ID'd objects are there?
    //-----------------------------------------------------------------

    int n_muonLoose          = c_muon_eta_IDLoose            -> GetSize();
    int n_muonHighPt         = c_muon_eta_IDHighPt           -> GetSize();
    int n_muon_store         = c_muon_final                  -> GetSize();
    int n_ele_store          = c_ele_final                   -> GetSize();
    int n_jet_store          = c_pfjet_final                 -> GetSize();
    int n_jet_highEta_store  = c_pfjet_highEta_final         -> GetSize();
    int n_genEle_store       = c_genEle_final                -> GetSize();
    int n_genNu_store        = c_genNu_final                 -> GetSize();
    int n_genMu_store        = c_genMu_final                 -> GetSize();
    int n_genJet_store       = c_genJet_final                -> GetSize();

    int n_muon_ptCut         = c_muon_final_ptCut            -> GetSize();
    int n_ele_ptCut          = c_ele_final_ptCut             -> GetSize();
    int n_ele_vloose_ptCut   = c_ele_vLoose_ptCut            -> GetSize();
    int n_jet_ptCut          = c_pfjet_final_ptCut           -> GetSize();
    int n_jet_highEta_ptCut  = c_pfjet_highEta_final_ptCut   -> GetSize();

    int n_genNuFromW_store   = c_genNuFromW_final            -> GetSize();

    int n_genLQ_store = c_genLQ_final ->GetSize();

    //if(n_ele_ptCut < 1) {
    //  std::cout << "NO GOOD ELECTRON FOUND! " << 
    //    static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(luminosityBlock) << " " << static_cast<unsigned int>(event) << std::endl;
    //  c_ele_all->examine<Electron>("c_ele_all");
    //}
    //-----------------------------------------------------------------
    // All skims need GEN particles/jets
    //-----------------------------------------------------------------

    if ( reducedSkimType != 0 ) {

      fillVariableWithValue("nGenJet_ptCut", n_genJet_store);
      fillVariableWithValue("nGenEle_ptCut", n_genEle_store);
      fillVariableWithValue("nGenNu_ptCut",  n_genNu_store);
      fillVariableWithValue("nGenMu_ptCut",  n_genMu_store);

      fillVariableWithValue("nGenJet_store", min(n_genJet_store,5));
      fillVariableWithValue("nGenEle_store", min(n_genEle_store,2));
      fillVariableWithValue("nGenNu_store" , min(n_genNu_store,2));
      fillVariableWithValue("nGenMu_store" , min(n_genMu_store,3));

      fillVariableWithValue("nGenNuFromW_ptCut"	 , n_genNuFromW_store   );

      fillVariableWithValue("nGenNuFromW_store"	 , min(n_genNuFromW_store  ,2));

      if ( n_genJet_store >= 1 ) { 
        GenJet genJet1 = c_genJet_final -> GetConstituent<GenJet>(0);
        fillVariableWithValue ( "GenJet1_Pt" , genJet1.Pt () );
        fillVariableWithValue ( "GenJet1_Eta", genJet1.Eta() );
        fillVariableWithValue ( "GenJet1_Phi", genJet1.Phi() );

        if ( n_genJet_store >= 2 ) { 
          GenJet genJet2 = c_genJet_final -> GetConstituent<GenJet>(1);
          fillVariableWithValue ( "GenJet2_Pt" , genJet2.Pt () );
          fillVariableWithValue ( "GenJet2_Eta", genJet2.Eta() );
          fillVariableWithValue ( "GenJet2_Phi", genJet2.Phi() );

          if ( n_genJet_store >= 3 ) { 
            GenJet genJet3 = c_genJet_final -> GetConstituent<GenJet>(2);
            fillVariableWithValue ( "GenJet3_Pt" , genJet3.Pt () );
            fillVariableWithValue ( "GenJet3_Eta", genJet3.Eta() );
            fillVariableWithValue ( "GenJet3_Phi", genJet3.Phi() );

            if ( n_genJet_store >= 4 ) { 
              GenJet genJet4 = c_genJet_final -> GetConstituent<GenJet>(3);
              fillVariableWithValue ( "GenJet4_Pt" , genJet4.Pt () );
              fillVariableWithValue ( "GenJet4_Eta", genJet4.Eta() );
              fillVariableWithValue ( "GenJet4_Phi", genJet4.Phi() );

              if ( n_genJet_store >= 5 ) { 
                GenJet genJet5 = c_genJet_final -> GetConstituent<GenJet>(4);
                fillVariableWithValue ( "GenJet5_Pt" , genJet5.Pt () );
                fillVariableWithValue ( "GenJet5_Eta", genJet5.Eta() );
                fillVariableWithValue ( "GenJet5_Phi", genJet5.Phi() );
              }
              //else {
              //  std::cout << "This event: " <<
              //    static_cast<unsigned int>(run) << " " << static_cast<unsigned int>(luminosityBlock) << " " << static_cast<unsigned int>(event) <<
              //    " had no 5th gen Jet; examine other GenJets." << std::endl;
              //  c_genJet_final->examine<GenJet>("finalGenJets");
              //}
            }
          }
        }
      }

      if ( n_genEle_store >= 1 ){ 
        GenParticle genEle1 = c_genEle_final -> GetConstituent<GenParticle>(0);
        fillVariableWithValue ( "GenEle1_Pt" , genEle1.Pt () );
        fillVariableWithValue ( "GenEle1_Eta", genEle1.Eta() );
        fillVariableWithValue ( "GenEle1_Phi", genEle1.Phi() );

        if ( n_genEle_store >= 2 ){ 
          GenParticle genEle2 = c_genEle_final -> GetConstituent<GenParticle>(1);
          fillVariableWithValue ( "GenEle2_Pt" , genEle2.Pt () );
          fillVariableWithValue ( "GenEle2_Eta", genEle2.Eta() );
          fillVariableWithValue ( "GenEle2_Phi", genEle2.Phi() );
        }
      }

      // neutrinos
      if ( n_genNu_store >= 1 ){ 
        GenParticle genNu1 = c_genNu_final -> GetConstituent<GenParticle>(0);
        fillVariableWithValue ( "GenNu1_Pt" , genNu1.Pt () );
        fillVariableWithValue ( "GenNu1_Eta", genNu1.Eta() );
        fillVariableWithValue ( "GenNu1_Phi", genNu1.Phi() );

        if ( n_genNu_store >= 2 ){ 
          GenParticle genNu2 = c_genNu_final -> GetConstituent<GenParticle>(1);
          fillVariableWithValue ( "GenNu2_Pt" , genNu2.Pt () );
          fillVariableWithValue ( "GenNu2_Eta", genNu2.Eta() );
          fillVariableWithValue ( "GenNu2_Phi", genNu2.Phi() );
        }
      }

      // muons
      if ( n_genMu_store >= 1 ){ 
        GenParticle genMu1 = c_genMu_final -> GetConstituent<GenParticle>(0);
        fillVariableWithValue ( "GenMu1_Pt" , genMu1.Pt () );
        fillVariableWithValue ( "GenMu1_Eta", genMu1.Eta() );
        fillVariableWithValue ( "GenMu1_Phi", genMu1.Phi() );

        if ( n_genMu_store >= 2 ){ 
          GenParticle genMu2 = c_genMu_final -> GetConstituent<GenParticle>(1);
          fillVariableWithValue ( "GenMu2_Pt" , genMu2.Pt () );
          fillVariableWithValue ( "GenMu2_Eta", genMu2.Eta() );
          fillVariableWithValue ( "GenMu2_Phi", genMu2.Phi() );
        }

        if ( n_genMu_store >= 3 ){ 
          GenParticle genMu3 = c_genMu_final -> GetConstituent<GenParticle>(2);
          fillVariableWithValue ( "GenMu3_Pt" , genMu3.Pt () );
          fillVariableWithValue ( "GenMu3_Eta", genMu3.Eta() );
          fillVariableWithValue ( "GenMu3_Phi", genMu3.Phi() );
        }
      }

      if ( n_genNuFromW_store >= 1 ){ 
        GenParticle genNuFromW1 = c_genNuFromW_final -> GetConstituent<GenParticle>(0);
        fillVariableWithValue ( "GenNuFromW1_Pt" , genNuFromW1.Pt () );
        fillVariableWithValue ( "GenNuFromW1_Eta", genNuFromW1.Eta() );
        fillVariableWithValue ( "GenNuFromW1_Phi", genNuFromW1.Phi() );
        fillVariableWithValue ( "GenNuFromW1_ID" , genNuFromW1.PdgId());

        if ( n_genNuFromW_store >= 2 ){ 
          GenParticle genNuFromW2 = c_genNuFromW_final -> GetConstituent<GenParticle>(1);
          fillVariableWithValue ( "GenNuFromW2_Pt" , genNuFromW2.Pt () );
          fillVariableWithValue ( "GenNuFromW2_Eta", genNuFromW2.Eta() );
          fillVariableWithValue ( "GenNuFromW2_Phi", genNuFromW2.Phi() );
          fillVariableWithValue ( "GenNuFromW2_ID" , genNuFromW2.PdgId());
        }
      }

      if ( n_genLQ_store >= 1 ){ 
        GenParticle genLQ1 = c_genLQ_final -> GetConstituent<GenParticle>(0);
        fillVariableWithValue ( "GenLQ1_Pt" , genLQ1.Pt () );
        fillVariableWithValue ( "GenLQ1_Eta", genLQ1.Eta() );
        fillVariableWithValue ( "GenLQ1_Phi", genLQ1.Phi() );
        fillVariableWithValue ( "GenLQ1_Mass", genLQ1.Mass() );
        fillVariableWithValue ( "GenLQ1_ID" , genLQ1.PdgId());

        if ( n_genLQ_store >= 2 ){ 
          GenParticle genLQ2 = c_genLQ_final -> GetConstituent<GenParticle>(1);
          fillVariableWithValue ( "GenLQ2_Pt" , genLQ2.Pt () );
          fillVariableWithValue ( "GenLQ2_Eta", genLQ2.Eta() );
          fillVariableWithValue ( "GenLQ2_Phi", genLQ2.Phi() );
          fillVariableWithValue ( "GenLQ2_Mass", genLQ2.Mass() );
          fillVariableWithValue ( "GenLQ2_ID" , genLQ2.PdgId());
        }
      }

    }

    //-----------------------------------------------------------------
    // All skims need muons
    //-----------------------------------------------------------------

    fillVariableWithValue ("nMuon_ptCut", n_muon_ptCut);
    fillVariableWithValue ("nMuon_LooseId", n_muonLoose);
    fillVariableWithValue ("nMuon_HighPtId", n_muonHighPt);
    fillVariableWithValue ("nMuon_store", min(n_muon_store,3));

    if ( n_muon_store >= 1 ){ 

      Muon muon1 = c_muon_final -> GetConstituent<Muon>(0);
      //double hltSingleMuon1Pt = triggerMatchPt<HLTriggerObject, Muon>(c_hltMuon_SingleMu_all, muon1, muon_hltMatch_DeltaRCut);
      fillVariableWithValue ("Muon1_Pt"             , muon1.Pt      ());
      fillVariableWithValue ("Muon1_Eta"            , muon1.Eta     ());
      fillVariableWithValue ("Muon1_Phi"            , muon1.Phi     ());
      fillVariableWithValue ("Muon1_PtError"        , muon1.PtError ());
      fillVariableWithValue ("Muon1_Charge"         , muon1.Charge  ());
      //fillVariableWithValue ("Muon1_hltSingleMuonPt", hltSingleMuon1Pt);

      if ( n_muon_store >= 2 ){ 

        Muon muon2 = c_muon_final -> GetConstituent<Muon>(1);
        //double hltSingleMuon2Pt = triggerMatchPt<HLTriggerObject, Muon>(c_hltMuon_SingleMu_all, muon2, muon_hltMatch_DeltaRCut);
        fillVariableWithValue ("Muon2_Pt"             , muon2.Pt      ());
        fillVariableWithValue ("Muon2_Eta"            , muon2.Eta     ());
        fillVariableWithValue ("Muon2_Phi"            , muon2.Phi     ());
        fillVariableWithValue ("Muon2_PtError"        , muon2.PtError ());
        fillVariableWithValue ("Muon2_Charge"         , muon2.Charge  ());
        //fillVariableWithValue ("Muon2_hltSingleMuonPt", hltSingleMuon2Pt);

        if ( n_muon_store >= 3 ){ 

          Muon muon3 = c_muon_final -> GetConstituent<Muon>(2);
          //double hltSingleMuon3Pt = triggerMatchPt<HLTriggerObject, Muon>(c_hltMuon_SingleMu_all, muon3, muon_hltMatch_DeltaRCut);
          fillVariableWithValue ("Muon3_Pt"             , muon3.Pt      ());
          fillVariableWithValue ("Muon3_Eta"            , muon3.Eta     ());
          fillVariableWithValue ("Muon3_Phi"            , muon3.Phi     ());
          fillVariableWithValue ("Muon3_PtError"        , muon3.PtError ());
          fillVariableWithValue ("Muon3_Charge"         , muon3.Charge  ());
          //fillVariableWithValue ("Muon3_hltSingleMuonPt", hltSingleMuon3Pt);
        }
      }
    }

    //-----------------------------------------------------------------
    // Fill variables for signal-like skim (reducedSkimType == 1 - 4 )
    //-----------------------------------------------------------------

    // Electrons
    if ( reducedSkimType == 1 || reducedSkimType == 2 || reducedSkimType == 3 || reducedSkimType == 4) { 

      fillVariableWithValue ("nEle_store"       , min(n_ele_store,3) );
      fillVariableWithValue ("nEle_ID"          , n_ele_store);
      fillVariableWithValue ("nJet_store"       , min(n_jet_store,5) );
      fillVariableWithValue ("nJet_etaIdLepCleaned"       , n_jet_store);
      fillVariableWithValue ("nHighEtaJet_store", min(n_jet_highEta_store, 1));
      fillVariableWithValue ("nEle_ptCut"       , n_ele_ptCut );
      fillVariableWithValue ("nJet_ptCut"       , n_jet_ptCut );
      fillVariableWithValue ("nHighEtaJet_ptCut", n_jet_highEta_ptCut );

      if ( n_ele_store >= 1 ){
        Electron ele1 = c_ele_final -> GetConstituent<Electron>(0);
        //if(fabs(ele1.Eta()) >= 1.4442 && fabs(ele1.Eta()) <= 1.566)
        //  c_ele_final->examine<Electron>("c_ele_final");
        //double hltEle1Pt_signal          = triggerMatchPt<HLTriggerObject, Electron>(c_hltEle45_Signal_all    , ele1, ele_hltMatch_DeltaRCut);
        //double hltEle1Pt_doubleEleSignal = triggerMatchPt<HLTriggerObject, Electron>(c_hltDoubleEle_Signal_all, ele1, ele_hltMatch_DeltaRCut);
        //double hltEle1Pt_WP80            = triggerMatchPt<HLTriggerObject, Electron>(c_hltEle27WP85Gsf_all       , ele1, ele_hltMatch_DeltaRCut);

        fillVariableWithValue( "Ele1_Pt"            , ele1.Pt()                 );
        fillVariableWithValue( "Ele1_SCEt"          , ele1.SCEnergy()/cosh(ele1.SCEta())               );
        fillVariableWithValue( "Ele1_ECorr"         , ele1.ECorr()              );
        fillVariableWithValue( "Ele1_Eta"           , ele1.Eta()                );
        fillVariableWithValue( "Ele1_Phi"           , ele1.Phi()                );
        fillVariableWithValue( "Ele1_SCEta"         , ele1.SCEta()              );
        fillVariableWithValue( "Ele1_Charge"        , ele1.Charge()             );
        fillVariableWithValue( "Ele1_R9"            , ele1.R9()                 );
        fillVariableWithValue( "Ele1_MissingHits"   , ele1.MissingHits()        );
        fillVariableWithValue( "Ele1_Full5x5SigmaIEtaIEta" , ele1.Full5x5SigmaIEtaIEta() );
        fillVariableWithValue( "Ele1_RhoForHEEP"    , ele1.RhoForHEEP()         );

        fillVariableWithValue( "Ele1_DeltaEtaTrkSC" , ele1.DeltaEta()           );
        fillVariableWithValue( "Ele1_HoE"           , ele1.HoE()                );
        fillVariableWithValue( "Ele1_HasMatchedPhot", ele1.HasMatchedConvPhot() );
        fillVariableWithValue( "Ele1_LeadVtxDistXY" , ele1.LeadVtxDistXY()      );
        fillVariableWithValue( "Ele1_LeadVtxDistZ"  , ele1.LeadVtxDistZ ()      );

        fillVariableWithValue( "Ele1_TrkIsolation"  , ele1.TrkIsoDR03()         );
        fillVariableWithValue( "Ele1_TrkIsoHEEP7"   , ele1.HEEP70TrackIsolation());
        fillVariableWithValue( "Ele1_EcalIsolation" , ele1.EcalIsoDR03()        );
        fillVariableWithValue( "Ele1_HcalIsolation" , ele1.HcalIsoD1DR03()      );
        fillVariableWithValue( "Ele1_CorrIsolation" , ele1.HEEPCorrIsolation()  );
        fillVariableWithValue( "Ele1_PFRelIso03Charged"     , ele1.PFRelIso03Charged());
        fillVariableWithValue( "Ele1_PFRelIso03All"     , ele1.PFRelIso03All());

        fillVariableWithValue( "Ele1_PassHEEPMinPtCut"                            ,ele1.PassHEEPMinPtCut                            () );
        fillVariableWithValue( "Ele1_PassHEEPGsfEleSCEtaMultiRangeCut"            ,ele1.PassHEEPGsfEleSCEtaMultiRangeCut            () ); 
        fillVariableWithValue( "Ele1_PassHEEPGsfEleDEtaInSeedCut"                 ,ele1.PassHEEPGsfEleDEtaInSeedCut                 () ); 
        fillVariableWithValue( "Ele1_PassHEEPGsfEleDPhiInCut"                     ,ele1.PassHEEPGsfEleDPhiInCut                     () ); 
        fillVariableWithValue( "Ele1_PassHEEPGsfEleFull5x5SigmaIEtaIEtaWithSatCut",ele1.PassHEEPGsfEleFull5x5SigmaIEtaIEtaWithSatCut() ); 
        fillVariableWithValue( "Ele1_PassHEEPGsfEleFull5x5E2x5OverE5x5WithSatCut" ,ele1.PassHEEPGsfEleFull5x5E2x5OverE5x5WithSatCut () ); 
        fillVariableWithValue( "Ele1_PassHEEPGsfEleTrkPtIsoCut"                   ,ele1.PassHEEPGsfEleTrkPtIsoCut                   () ); 
        if(analysisYear == 2018) {
          fillVariableWithValue( "Ele1_PassHEEPGsfEleEmHadD1IsoRhoCut"              ,ele1.PassHEEPGsfEleEmHadD1IsoRhoCut2018        () ); 
          fillVariableWithValue( "Ele1_PassHEEPGsfEleHadronicOverEMLinearCut"       ,ele1.PassHEEPGsfEleHadronicOverEMLinearCut2018 () ); 
        }
        else {
          fillVariableWithValue( "Ele1_PassHEEPGsfEleEmHadD1IsoRhoCut"              ,ele1.PassHEEPGsfEleEmHadD1IsoRhoCut              () ); 
          fillVariableWithValue( "Ele1_PassHEEPGsfEleHadronicOverEMLinearCut"       ,ele1.PassHEEPGsfEleHadronicOverEMLinearCut       () ); 
        }
        fillVariableWithValue( "Ele1_PassHEEPGsfEleDxyCut"                        ,ele1.PassHEEPGsfEleDxyCut                        () ); 
        fillVariableWithValue( "Ele1_PassHEEPGsfEleMissingHitsCut"                ,ele1.PassHEEPGsfEleMissingHitsCut                () ); 
        fillVariableWithValue( "Ele1_PassHEEPEcalDrivenCut"                       ,ele1.PassHEEPEcalDrivenCut                       () );

        //fillVariableWithValue( "Ele1_hltEleSignalPt", hltEle1Pt_signal          );
        //fillVariableWithValue( "Ele1_hltDoubleElePt", hltEle1Pt_doubleEleSignal ); 
        //fillVariableWithValue( "Ele1_hltEleWP80Pt"  , hltEle1Pt_WP80            );
        //if(!isData()) {
        //  Electron ele1smeared = *find(smearedEles.begin(), smearedEles.end(), ele1);
        //  Electron ele1scaledUp = *find(scaledUpEles.begin(), scaledUpEles.end(), ele1);
        //  Electron ele1scaledDown = *find(scaledDownEles.begin(), scaledDownEles.end(), ele1);
        //  fillVariableWithValue( "Ele1_Pt_EER"      , ele1smeared.Pt()          );
        //  fillVariableWithValue( "Ele1_Pt_EES_Up"    , ele1scaledUp.Pt()         );
        //  fillVariableWithValue( "Ele1_Pt_EES_Dn"  , ele1scaledDown.Pt()       );
        //  float elePtUncorr = ele1.ECorr() != 0 ? ele1.Pt()/ele1.ECorr() : ele1.Pt();
        //  fillVariableWithValue( "Ele1_TrigSF" , triggerScaleFactorReader.LookupValue(ele1.SCEta(),elePtUncorr) );
        //  fillVariableWithValue( "Ele1_TrigSF_Err" , triggerScaleFactorReader.LookupValueError(ele1.SCEta(),elePtUncorr) );
        //  fillVariableWithValue( "Ele1_RecoSF" , recoScaleFactorReader->LookupValue(ele1.SCEta(),elePtUncorr) );
        //  fillVariableWithValue( "Ele1_RecoSF_Err" , recoScaleFactorReader->LookupValueError(ele1.SCEta(),elePtUncorr) );
        //  fillVariableWithValue( "Ele1_HEEPSF"                                      ,ElectronScaleFactorsRunII::LookupHeepSF(ele1.SCEta(), analysisYear) );
        //  fillVariableWithValue( "Ele1_HEEPSF_Err"                                  ,ElectronScaleFactorsRunII::LookupHeepSFSyst(ele1.SCEta(), elePtUncorr, analysisYear) );
        //}

        if ( n_ele_store >= 2 ){
          Electron ele2 = c_ele_final -> GetConstituent<Electron>(1);
          //double hltEle2Pt_signal          = triggerMatchPt<HLTriggerObject, Electron>(c_hltEle45_Signal_all     , ele2, ele_hltMatch_DeltaRCut);
          //double hltEle2Pt_doubleEleSignal = triggerMatchPt<HLTriggerObject, Electron>(c_hltDoubleEle_Signal_all , ele2, ele_hltMatch_DeltaRCut);
          //double hltEle2Pt_WP80            = triggerMatchPt<HLTriggerObject, Electron>(c_hltEle27WP85Gsf_all        , ele2, ele_hltMatch_DeltaRCut);

          fillVariableWithValue( "Ele2_Pt"            , ele2.Pt()                 );
          fillVariableWithValue( "Ele2_SCEt"          , ele2.SCEnergy()/cosh(ele2.SCEta())               );
          fillVariableWithValue( "Ele2_ECorr"         , ele2.ECorr()              );
          fillVariableWithValue( "Ele2_Eta"           , ele2.Eta()                );
          fillVariableWithValue( "Ele2_Phi"           , ele2.Phi()                );
          fillVariableWithValue( "Ele2_SCEta"         , ele2.SCEta()              );
          fillVariableWithValue( "Ele2_Charge"        , ele2.Charge()             );
          fillVariableWithValue( "Ele2_R9"            , ele2.R9()                 );
          fillVariableWithValue( "Ele2_MissingHits"   , ele2.MissingHits()        );
          fillVariableWithValue( "Ele2_Full5x5SigmaIEtaIEta" , ele2.Full5x5SigmaIEtaIEta() );
          fillVariableWithValue( "Ele2_RhoForHEEP"    , ele2.RhoForHEEP()         );

          fillVariableWithValue( "Ele2_DeltaEtaTrkSC" , ele2.DeltaEta()           );
          fillVariableWithValue( "Ele2_HoE"           , ele2.HoE()                );
          fillVariableWithValue( "Ele2_HasMatchedPhot", ele2.HasMatchedConvPhot() );
          fillVariableWithValue( "Ele2_LeadVtxDistXY" , ele2.LeadVtxDistXY()      );
          fillVariableWithValue( "Ele2_LeadVtxDistZ"  , ele2.LeadVtxDistZ ()      );

          fillVariableWithValue( "Ele2_TrkIsolation"  , ele2.TrkIsoDR03()         );
          fillVariableWithValue( "Ele2_TrkIsoHEEP7"   , ele2.HEEP70TrackIsolation());
          fillVariableWithValue( "Ele2_EcalIsolation" , ele2.EcalIsoDR03()        );
          fillVariableWithValue( "Ele2_HcalIsolation" , ele2.HcalIsoD1DR03()      );
          fillVariableWithValue( "Ele2_CorrIsolation" , ele2.HEEPCorrIsolation()  );
        fillVariableWithValue( "Ele2_PFRelIso03Charged"     , ele2.PFRelIso03Charged());
        fillVariableWithValue( "Ele2_PFRelIso03All"     , ele2.PFRelIso03All());

          fillVariableWithValue( "Ele2_PassHEEPMinPtCut"                            ,ele2.PassHEEPMinPtCut                            () );
          fillVariableWithValue( "Ele2_PassHEEPGsfEleSCEtaMultiRangeCut"            ,ele2.PassHEEPGsfEleSCEtaMultiRangeCut            () ); 
          fillVariableWithValue( "Ele2_PassHEEPGsfEleDEtaInSeedCut"                 ,ele2.PassHEEPGsfEleDEtaInSeedCut                 () ); 
          fillVariableWithValue( "Ele2_PassHEEPGsfEleDPhiInCut"                     ,ele2.PassHEEPGsfEleDPhiInCut                     () ); 
          fillVariableWithValue( "Ele2_PassHEEPGsfEleFull5x5SigmaIEtaIEtaWithSatCut",ele2.PassHEEPGsfEleFull5x5SigmaIEtaIEtaWithSatCut() ); 
          fillVariableWithValue( "Ele2_PassHEEPGsfEleFull5x5E2x5OverE5x5WithSatCut" ,ele2.PassHEEPGsfEleFull5x5E2x5OverE5x5WithSatCut () ); 
          fillVariableWithValue( "Ele2_PassHEEPGsfEleTrkPtIsoCut"                   ,ele2.PassHEEPGsfEleTrkPtIsoCut                   () ); 
          if(analysisYear == 2018) {
            fillVariableWithValue( "Ele2_PassHEEPGsfEleEmHadD1IsoRhoCut"              ,ele2.PassHEEPGsfEleEmHadD1IsoRhoCut2018        () ); 
            fillVariableWithValue( "Ele2_PassHEEPGsfEleHadronicOverEMLinearCut"       ,ele2.PassHEEPGsfEleHadronicOverEMLinearCut2018 () ); 
          }
          else {
            fillVariableWithValue( "Ele2_PassHEEPGsfEleEmHadD1IsoRhoCut"              ,ele2.PassHEEPGsfEleEmHadD1IsoRhoCut              () ); 
            fillVariableWithValue( "Ele2_PassHEEPGsfEleHadronicOverEMLinearCut"       ,ele2.PassHEEPGsfEleHadronicOverEMLinearCut       () ); 
          }
          fillVariableWithValue( "Ele2_PassHEEPGsfEleDxyCut"                        ,ele2.PassHEEPGsfEleDxyCut                        () ); 
          fillVariableWithValue( "Ele2_PassHEEPGsfEleMissingHitsCut"                ,ele2.PassHEEPGsfEleMissingHitsCut                () ); 
          fillVariableWithValue( "Ele2_PassHEEPEcalDrivenCut"                       ,ele2.PassHEEPEcalDrivenCut                       () );

          //fillVariableWithValue( "Ele2_hltEleSignalPt", hltEle2Pt_signal          );
          //fillVariableWithValue( "Ele2_hltDoubleElePt", hltEle2Pt_doubleEleSignal ); 
          //fillVariableWithValue( "Ele2_hltEleWP80Pt"  , hltEle2Pt_WP80            );
          //if(!isData()) {
          //  Electron ele2smeared = *find(smearedEles.begin(), smearedEles.end(), ele2);
          //  Electron ele2scaledUp = *find(scaledUpEles.begin(), scaledUpEles.end(), ele2);
          //  Electron ele2scaledDown = *find(scaledDownEles.begin(), scaledDownEles.end(), ele2);
          //  fillVariableWithValue( "Ele2_Pt_EER"      , ele2smeared.Pt()          );
          //  fillVariableWithValue( "Ele2_Pt_EES_Up"    , ele2scaledUp.Pt()         );
          //  fillVariableWithValue( "Ele2_Pt_EES_Dn"  , ele2scaledDown.Pt()       );
          //  float elePtUncorr = ele2.ECorr() != 0 ? ele2.Pt()/ele2.ECorr() : ele2.Pt();
          //  fillVariableWithValue( "Ele2_TrigSF" , triggerScaleFactorReader.LookupValue(ele2.SCEta(),elePtUncorr) );
          //  fillVariableWithValue( "Ele2_TrigSF_Err" , triggerScaleFactorReader.LookupValueError(ele2.SCEta(),elePtUncorr) );
          //  fillVariableWithValue( "Ele2_RecoSF" , recoScaleFactorReader->LookupValue(ele2.SCEta(),elePtUncorr) );
          //  fillVariableWithValue( "Ele2_RecoSF_Err" , recoScaleFactorReader->LookupValueError(ele2.SCEta(),elePtUncorr) );
          //  fillVariableWithValue( "Ele2_HEEPSF"                                      ,ElectronScaleFactorsRunII::LookupHeepSF(ele2.SCEta(), analysisYear) );
          //  fillVariableWithValue( "Ele2_HEEPSF_Err"                                  ,ElectronScaleFactorsRunII::LookupHeepSFSyst(ele2.SCEta(), elePtUncorr, analysisYear) );
          //}

          if ( n_ele_store >= 3 ){
            Electron ele3 = c_ele_final -> GetConstituent<Electron>(2);
            //double hltEle3Pt_signal          = triggerMatchPt<HLTriggerObject, Electron>(c_hltEle45_Signal_all     , ele3, ele_hltMatch_DeltaRCut);
            //double hltEle3Pt_doubleEleSignal = triggerMatchPt<HLTriggerObject, Electron>(c_hltDoubleEle_Signal_all , ele3, ele_hltMatch_DeltaRCut);
            //double hltEle3Pt_WP80            = triggerMatchPt<HLTriggerObject, Electron>(c_hltEle27WP85Gsf_all        , ele3, ele_hltMatch_DeltaRCut);

            fillVariableWithValue( "Ele3_Pt"            , ele3.Pt()                 );
            //fillVariableWithValue( "Ele3_SCEt"          , ele3.SCEnergy()/cosh(ele3.SCEta())               );
            fillVariableWithValue( "Ele3_ECorr"         , ele3.ECorr()              );
            fillVariableWithValue( "Ele3_Eta"           , ele3.Eta()                );
            fillVariableWithValue( "Ele3_Phi"           , ele3.Phi()                );
            fillVariableWithValue( "Ele3_SCEta"         , ele3.SCEta()              );
            fillVariableWithValue( "Ele3_Charge"        , ele3.Charge()             );
            fillVariableWithValue( "Ele3_R9"            , ele3.R9()                 );
            fillVariableWithValue( "Ele3_MissingHits"   , ele3.MissingHits()        );
            fillVariableWithValue( "Ele3_Full5x5SigmaIEtaIEta" , ele3.Full5x5SigmaIEtaIEta() );
            fillVariableWithValue( "Ele3_RhoForHEEP"    , ele3.RhoForHEEP()         );

            fillVariableWithValue( "Ele3_DeltaEtaTrkSC" , ele3.DeltaEta()           );
            fillVariableWithValue( "Ele3_HoE"           , ele3.HoE()                );
            fillVariableWithValue( "Ele3_HasMatchedPhot", ele3.HasMatchedConvPhot() );
            fillVariableWithValue( "Ele3_LeadVtxDistXY" , ele3.LeadVtxDistXY()      );
            fillVariableWithValue( "Ele3_LeadVtxDistZ"  , ele3.LeadVtxDistZ ()      );

            fillVariableWithValue( "Ele3_TrkIsolation"  , ele3.TrkIsoDR03()         );
            fillVariableWithValue( "Ele3_TrkIsoHEEP7"   , ele3.HEEP70TrackIsolation());
            fillVariableWithValue( "Ele3_EcalIsolation" , ele3.EcalIsoDR03()        );
            fillVariableWithValue( "Ele3_HcalIsolation" , ele3.HcalIsoD1DR03()      );
            fillVariableWithValue( "Ele3_CorrIsolation" , ele3.HEEPCorrIsolation()  );
            fillVariableWithValue( "Ele3_PFRelIso03Charged"     , ele3.PFRelIso03Charged());
            fillVariableWithValue( "Ele3_PFRelIso03All"     , ele3.PFRelIso03All());

            fillVariableWithValue( "Ele3_PassHEEPMinPtCut"                            ,ele3.PassHEEPMinPtCut                            () );
            fillVariableWithValue( "Ele3_PassHEEPGsfEleSCEtaMultiRangeCut"            ,ele3.PassHEEPGsfEleSCEtaMultiRangeCut            () ); 
            fillVariableWithValue( "Ele3_PassHEEPGsfEleDEtaInSeedCut"                 ,ele3.PassHEEPGsfEleDEtaInSeedCut                 () ); 
            fillVariableWithValue( "Ele3_PassHEEPGsfEleDPhiInCut"                     ,ele3.PassHEEPGsfEleDPhiInCut                     () ); 
            fillVariableWithValue( "Ele3_PassHEEPGsfEleFull5x5SigmaIEtaIEtaWithSatCut",ele3.PassHEEPGsfEleFull5x5SigmaIEtaIEtaWithSatCut() ); 
            fillVariableWithValue( "Ele3_PassHEEPGsfEleFull5x5E2x5OverE5x5WithSatCut" ,ele3.PassHEEPGsfEleFull5x5E2x5OverE5x5WithSatCut () ); 
            fillVariableWithValue( "Ele3_PassHEEPGsfEleTrkPtIsoCut"                   ,ele3.PassHEEPGsfEleTrkPtIsoCut                   () ); 
            if(analysisYear == 2018) {
              fillVariableWithValue( "Ele3_PassHEEPGsfEleEmHadD1IsoRhoCut"              ,ele3.PassHEEPGsfEleEmHadD1IsoRhoCut2018        () ); 
              fillVariableWithValue( "Ele3_PassHEEPGsfEleHadronicOverEMLinearCut"       ,ele3.PassHEEPGsfEleHadronicOverEMLinearCut2018 () ); 
            }
            else {
              fillVariableWithValue( "Ele3_PassHEEPGsfEleEmHadD1IsoRhoCut"              ,ele3.PassHEEPGsfEleEmHadD1IsoRhoCut              () ); 
              fillVariableWithValue( "Ele3_PassHEEPGsfEleHadronicOverEMLinearCut"       ,ele3.PassHEEPGsfEleHadronicOverEMLinearCut       () ); 
            }
            fillVariableWithValue( "Ele3_PassHEEPGsfEleDxyCut"                        ,ele3.PassHEEPGsfEleDxyCut                        () ); 
            fillVariableWithValue( "Ele3_PassHEEPGsfEleMissingHitsCut"                ,ele3.PassHEEPGsfEleMissingHitsCut                () ); 
            fillVariableWithValue( "Ele3_PassHEEPEcalDrivenCut"                       ,ele3.PassHEEPEcalDrivenCut                       () );

            //fillVariableWithValue( "Ele3_hltEleSignalPt", hltEle3Pt_signal          );
            //fillVariableWithValue( "Ele3_hltDoubleElePt", hltEle3Pt_doubleEleSignal ); 
            //fillVariableWithValue( "Ele3_hltEleWP80Pt"  , hltEle3Pt_WP80            );
            //if(!isData()) {
            //  Electron ele3smeared = *find(smearedEles.begin(), smearedEles.end(), ele3);
            //  Electron ele3scaledUp = *find(scaledUpEles.begin(), scaledUpEles.end(), ele3);
            //  Electron ele3scaledDown = *find(scaledDownEles.begin(), scaledDownEles.end(), ele3);
            //  fillVariableWithValue( "Ele3_Pt_EER"      , ele3smeared.Pt()          );
            //  fillVariableWithValue( "Ele3_Pt_EES_Up"    , ele3scaledUp.Pt()         );
            //  fillVariableWithValue( "Ele3_Pt_EES_Dn"  , ele3scaledDown.Pt()       );
            //  float elePtUncorr = ele3.ECorr() != 0 ? ele3.Pt()/ele3.ECorr() : ele3.Pt();
            //  fillVariableWithValue( "Ele3_TrigSF" , triggerScaleFactorReader.LookupValue(ele3.SCEta(),elePtUncorr) );
            //  fillVariableWithValue( "Ele3_TrigSF_Err" , triggerScaleFactorReader.LookupValueError(ele3.SCEta(),elePtUncorr) );
            //  fillVariableWithValue( "Ele3_RecoSF" , recoScaleFactorReader->LookupValue(ele3.SCEta(),elePtUncorr) );
            //  fillVariableWithValue( "Ele3_RecoSF_Err" , recoScaleFactorReader->LookupValueError(ele3.SCEta(),elePtUncorr) );
            //  fillVariableWithValue( "Ele3_HEEPSF"                                      ,ElectronScaleFactorsRunII::LookupHeepSF(ele3.SCEta(), analysisYear) );
            //  fillVariableWithValue( "Ele3_HEEPSF_Err"                                  ,ElectronScaleFactorsRunII::LookupHeepSFSyst(ele3.SCEta(), elePtUncorr, analysisYear) );
            //}

          }
        }
      }

      // Jets

      if ( n_jet_highEta_store >= 1 ) { 
        PFJet jet1 = c_pfjet_highEta_final -> GetConstituent<PFJet>(0);
        fillVariableWithValue( "HighEtaJet1_Pt" , jet1.Pt () );
        fillVariableWithValue( "HighEtaJet1_Eta", jet1.Eta() );
        fillVariableWithValue( "HighEtaJet1_Phi", jet1.Phi() );
      }

      if ( n_jet_store >= 1 ){

        PFJet jet1 = c_pfjet_final -> GetConstituent<PFJet>(0);
        // leading HLT jet from 200 GeV collection
        double hltJet1Pt     = triggerMatchPt<HLTriggerObject, PFJet >( c_trigger_l3jets_all     , jet1, jet_hltMatch_DeltaRCut);
        fillVariableWithValue( "Jet1_Pt"      , jet1.Pt()                         );
        fillVariableWithValue( "Jet1_Eta"     , jet1.Eta()                        );
        fillVariableWithValue( "Jet1_Phi"     , jet1.Phi()                        );
        fillVariableWithValue( "Jet1_btagDeepCSV" , jet1.DeepCSVBTag()            );
        fillVariableWithValue( "Jet1_btagDeepJet" , jet1.DeepJetBTag()            );
        fillVariableWithValue( "Jet1_hltJetPt"	  , hltJet1Pt     );
        //fillVariableWithValue( "Jet1_qgl"	  , jet1.QuarkGluonLikelihood()     );
        if(!isData()) {
          fillVariableWithValue( "Jet1_btagSFLooseDeepCSV"  , jet1.DeepCSVBTagSFLoose()   );
          fillVariableWithValue( "Jet1_btagSFMediumDeepCSV"  , jet1.DeepCSVBTagSFMedium()   );
          fillVariableWithValue( "Jet1_btagSFLooseDeepJet"  , jet1.DeepJetBTagSFLoose()   );
          fillVariableWithValue( "Jet1_btagSFLooseDeepJet_Up"  , jet1.DeepJetBTagSFLooseUp()   );
          fillVariableWithValue( "Jet1_btagSFLooseDeepJet_Dn"  , jet1.DeepJetBTagSFLooseDown()   );
          fillVariableWithValue( "Jet1_btagSFMediumDeepJet"  , jet1.DeepJetBTagSFMedium()   );
          fillVariableWithValue( "Jet1_Pt_JESTotal_Up"  , jet1.PtJESTotalUp()        );
          fillVariableWithValue( "Jet1_Pt_JESTotal_Dn", jet1.PtJESTotalDown()        );
          fillVariableWithValue( "Jet1_Pt_JER_Up"  , jet1.PtJERUp()                  );
          fillVariableWithValue( "Jet1_Pt_JER_Dn", jet1.PtJERDown()                  );
          if(c_pfJetMatchedLQ->Has<PFJet>(jet1))
            fillVariableWithValue( "Jet1_LQMatched", 1);
        }

        if ( n_jet_store >= 2 ){
          PFJet jet2 = c_pfjet_final -> GetConstituent<PFJet>(1);
          double hltJet2Pt     = triggerMatchPt<HLTriggerObject, PFJet >( c_trigger_l3jets_all     , jet2, jet_hltMatch_DeltaRCut);
          fillVariableWithValue( "Jet2_Pt"      , jet2.Pt()                         );
          fillVariableWithValue( "Jet2_Eta"     , jet2.Eta()                        );
          fillVariableWithValue( "Jet2_Phi"     , jet2.Phi()                        );
          fillVariableWithValue( "Jet2_btagDeepCSV" , jet2.DeepCSVBTag()            );
          fillVariableWithValue( "Jet2_btagDeepJet" , jet2.DeepJetBTag()            );
          fillVariableWithValue( "Jet2_hltJetPt"    , hltJet2Pt     );
          //fillVariableWithValue( "Jet2_qgl"	  , jet2.QuarkGluonLikelihood()     );
          if(!isData()) {
            fillVariableWithValue( "Jet2_btagSFLooseDeepCSV"  , jet2.DeepCSVBTagSFLoose()   );
            fillVariableWithValue( "Jet2_btagSFMediumDeepCSV"  , jet2.DeepCSVBTagSFMedium()   );
            fillVariableWithValue( "Jet2_btagSFLooseDeepJet"  , jet2.DeepJetBTagSFLoose()   );
            fillVariableWithValue( "Jet2_btagSFLooseDeepJet_Up"  , jet2.DeepJetBTagSFLooseUp()   );
            fillVariableWithValue( "Jet2_btagSFLooseDeepJet_Dn"  , jet2.DeepJetBTagSFLooseDown()   );
            fillVariableWithValue( "Jet2_btagSFMediumDeepJet"  , jet2.DeepJetBTagSFMedium()   );
            fillVariableWithValue( "Jet2_Pt_JESTotal_Up"  , jet2.PtJESTotalUp()        );
            fillVariableWithValue( "Jet2_Pt_JESTotal_Dn", jet2.PtJESTotalDown()        );
            fillVariableWithValue( "Jet2_Pt_JER_Up"  , jet2.PtJERUp()                  );
            fillVariableWithValue( "Jet2_Pt_JER_Dn", jet2.PtJERDown()                  );
            if(c_pfJetMatchedLQ->Has<PFJet>(jet2))
              fillVariableWithValue( "Jet2_LQMatched", 1);
          }

          if ( n_jet_store >= 3 ){
            PFJet jet3 = c_pfjet_final -> GetConstituent<PFJet>(2);
            double hltJet3Pt     = triggerMatchPt<HLTriggerObject, PFJet >( c_trigger_l3jets_all     , jet3, jet_hltMatch_DeltaRCut);
            fillVariableWithValue( "Jet3_Pt"      , jet3.Pt()                         );
            fillVariableWithValue( "Jet3_Eta"     , jet3.Eta()                        );
            fillVariableWithValue( "Jet3_Phi"     , jet3.Phi()                        );
            fillVariableWithValue( "Jet3_btagDeepCSV" , jet3.DeepCSVBTag()            );
            fillVariableWithValue( "Jet3_btagDeepJet" , jet3.DeepJetBTag()            );
            fillVariableWithValue( "Jet3_hltJetPt"    , hltJet3Pt     );
            //fillVariableWithValue( "Jet3_qgl"	  , jet3.QuarkGluonLikelihood()     );
            if(!isData()) {
              fillVariableWithValue( "Jet3_btagSFLooseDeepCSV"  , jet3.DeepCSVBTagSFLoose()   );
              fillVariableWithValue( "Jet3_btagSFMediumDeepCSV"  , jet3.DeepCSVBTagSFMedium()   );
              fillVariableWithValue( "Jet3_btagSFLooseDeepJet"  , jet3.DeepJetBTagSFLoose()   );
              fillVariableWithValue( "Jet3_btagSFLooseDeepJet_Up"  , jet3.DeepJetBTagSFLooseUp()   );
              fillVariableWithValue( "Jet3_btagSFLooseDeepJet_Dn"  , jet3.DeepJetBTagSFLooseDown()   );
              fillVariableWithValue( "Jet3_btagSFMediumDeepJet"  , jet3.DeepJetBTagSFMedium()   );
              fillVariableWithValue( "Jet3_Pt_JESTotal_Up"  , jet3.PtJESTotalUp()        );
              fillVariableWithValue( "Jet3_Pt_JESTotal_Dn", jet3.PtJESTotalDown()        );
              fillVariableWithValue( "Jet3_Pt_JER_Up"  , jet3.PtJERUp()                  );
              fillVariableWithValue( "Jet3_Pt_JER_Dn", jet3.PtJERDown()                  );
              if(c_pfJetMatchedLQ->Has<PFJet>(jet3))
                fillVariableWithValue( "Jet3_LQMatched", 1);
            }

            if ( n_jet_store >= 4 ){
              PFJet jet4 = c_pfjet_final -> GetConstituent<PFJet>(3);
              double hltJet4Pt     = triggerMatchPt<HLTriggerObject, PFJet >( c_trigger_l3jets_all     , jet4, jet_hltMatch_DeltaRCut);
              fillVariableWithValue( "Jet4_Pt"      , jet4.Pt()                         );
              fillVariableWithValue( "Jet4_Eta"     , jet4.Eta()                        );
              fillVariableWithValue( "Jet4_Phi"     , jet4.Phi()                        );
              fillVariableWithValue( "Jet4_btagDeepCSV" , jet4.DeepCSVBTag()            );
              fillVariableWithValue( "Jet4_btagDeepJet" , jet4.DeepJetBTag()            );
              fillVariableWithValue( "Jet4_hltJetPt"    , hltJet4Pt     );
              //fillVariableWithValue( "Jet4_qgl"	  , jet4.QuarkGluonLikelihood()     );
              if(!isData()) {
                fillVariableWithValue( "Jet4_btagSFLooseDeepCSV"  , jet4.DeepCSVBTagSFLoose()   );
                fillVariableWithValue( "Jet4_btagSFMediumDeepCSV"  , jet4.DeepCSVBTagSFMedium()   );
                fillVariableWithValue( "Jet4_btagSFLooseDeepJet"  , jet4.DeepJetBTagSFLoose()   );
                fillVariableWithValue( "Jet4_btagSFLooseDeepJet_Up"  , jet4.DeepJetBTagSFLooseUp()   );
                fillVariableWithValue( "Jet4_btagSFLooseDeepJet_Dn"  , jet4.DeepJetBTagSFLooseDown()   );
                fillVariableWithValue( "Jet4_btagSFMediumDeepJet"  , jet4.DeepJetBTagSFMedium()   );
                fillVariableWithValue( "Jet4_Pt_JESTotal_Up"  , jet4.PtJESTotalUp()        );
                fillVariableWithValue( "Jet4_Pt_JESTotal_Dn", jet4.PtJESTotalDown()        );
                fillVariableWithValue( "Jet4_Pt_JER_Up"  , jet4.PtJERUp()                  );
                fillVariableWithValue( "Jet4_Pt_JER_Dn", jet4.PtJERDown()                  );
                if(c_pfJetMatchedLQ->Has<PFJet>(jet4))
                  fillVariableWithValue( "Jet4_LQMatched", 1);
              }

              if ( n_jet_store >= 5 ){
                PFJet jet5 = c_pfjet_final -> GetConstituent<PFJet>(4);
                double hltJet5Pt     = triggerMatchPt<HLTriggerObject, PFJet >( c_trigger_l3jets_all     , jet5, jet_hltMatch_DeltaRCut);
                fillVariableWithValue( "Jet5_Pt"      , jet5.Pt()                         );
                fillVariableWithValue( "Jet5_Eta"     , jet5.Eta()                        );
                fillVariableWithValue( "Jet5_Phi"     , jet5.Phi()                        );
                fillVariableWithValue( "Jet5_btagDeepCSV" , jet5.DeepCSVBTag()            );
                fillVariableWithValue( "Jet5_btagDeepJet" , jet5.DeepJetBTag()            );
                fillVariableWithValue( "Jet5_hltJetPt"    , hltJet5Pt     );
                //fillVariableWithValue( "Jet5_qgl"	  , jet5.QuarkGluonLikelihood()     );
                if(!isData()) {
                  fillVariableWithValue( "Jet5_btagSFLooseDeepCSV"  , jet5.DeepCSVBTagSFLoose()   );
                  fillVariableWithValue( "Jet5_btagSFMediumDeepCSV"  , jet5.DeepCSVBTagSFMedium()   );
                  fillVariableWithValue( "Jet5_btagSFLooseDeepJet"  , jet5.DeepJetBTagSFLoose()   );
                  fillVariableWithValue( "Jet5_btagSFLooseDeepJet_Up"  , jet5.DeepJetBTagSFLooseUp()   );
                  fillVariableWithValue( "Jet5_btagSFLooseDeepJet_Dn"  , jet5.DeepJetBTagSFLooseDown()   );
                  fillVariableWithValue( "Jet5_btagSFMediumDeepJet"  , jet5.DeepJetBTagSFMedium()   );
                  fillVariableWithValue( "Jet5_Pt_JESTotal_Up"  , jet5.PtJESTotalUp()        );
                  fillVariableWithValue( "Jet5_Pt_JESTotal_Dn", jet5.PtJESTotalDown()        );
                  fillVariableWithValue( "Jet5_Pt_JER_Up"  , jet5.PtJERUp()                  );
                  fillVariableWithValue( "Jet5_Pt_JER_Dn", jet5.PtJERDown()                  );
                  if(c_pfJetMatchedLQ->Has<PFJet>(jet5))
                    fillVariableWithValue( "Jet5_LQMatched", 1);
                }
              }
            }
          }
        }
      }
    }

    //-----------------------------------------------------------------
    // Fill variables that depend on more than one object
    // All skims need this
    //-----------------------------------------------------------------

    std::map<short, short> jetRawIndicesMap;
    TLorentzVector t_ele1, t_ele2, t_jet1, t_jet2, t_jet3, t_jet4, t_jet5, t_jet6;
    TLorentzVector t_MET;
    TLorentzVector t_ele1Smeared, t_ele1ScaledUp, t_ele1ScaledDown;
    TLorentzVector t_ele2Smeared, t_ele2ScaledUp, t_ele2ScaledDown;
    TLorentzVector t_jet1JESTotalUp, t_jet1JESTotalDown, t_jet1JERUp, t_jet1JERDown;
    TLorentzVector t_jet2JESTotalUp, t_jet2JESTotalDown, t_jet2JERUp, t_jet2JERDown;

    t_MET.SetPtEtaPhiM( readerTools_->ReadValueBranch<Float_t>("MET_pt"), 0.0, readerTools_->ReadValueBranch<Float_t>("MET_phi"), 0.0 );

    if ( n_jet_store >= 1 ){

      PFJet jet1 = c_pfjet_final -> GetConstituent<PFJet>(0);
      jetRawIndicesMap[0] = jet1.GetRawIndex();
      t_jet1.SetPtEtaPhiM ( jet1.Pt(), jet1.Eta(), jet1.Phi(), 0.0 );
      
      //if(!isData()) {
      //  t_jet1JESTotalUp.SetPtEtaPhiM ( jet1.PtJESTotalUp(), jet1.Eta(), jet1.Phi(), 0.0 );
      //  t_jet1JESTotalDown.SetPtEtaPhiM ( jet1.PtJESTotalDown(), jet1.Eta(), jet1.Phi(), 0.0 );
      //  t_jet1JERUp.SetPtEtaPhiM ( jet1.PtJERUp(), jet1.Eta(), jet1.Phi(), 0.0 );
      //  t_jet1JERDown.SetPtEtaPhiM ( jet1.PtJERDown(), jet1.Eta(), jet1.Phi(), 0.0 );
      //}

      fillVariableWithValue ("mDPhi_METJet1", fabs( t_MET.DeltaPhi ( t_jet1 )));

      if ( n_jet_store >= 2 ){

        PFJet jet2 = c_pfjet_final -> GetConstituent<PFJet>(1);
        jetRawIndicesMap[1] = jet2.GetRawIndex();
        t_jet2.SetPtEtaPhiM ( jet2.Pt(), jet2.Eta(), jet2.Phi(), 0.0 );
        //if(!isData()) {
        //  t_jet2JESTotalUp.SetPtEtaPhiM ( jet2.PtJESTotalUp(), jet2.Eta(), jet2.Phi(), 0.0 );
        //  t_jet2JESTotalDown.SetPtEtaPhiM ( jet2.PtJESTotalDown(), jet2.Eta(), jet2.Phi(), 0.0 );
        //  t_jet2JERUp.SetPtEtaPhiM ( jet2.PtJERUp(), jet2.Eta(), jet2.Phi(), 0.0 );
        //  t_jet2JERDown.SetPtEtaPhiM ( jet2.PtJERDown(), jet2.Eta(), jet2.Phi(), 0.0 );
        //}

        TLorentzVector t_jet1jet2 = t_jet1 + t_jet2;

        fillVariableWithValue ("M_j1j2" , t_jet1jet2.M ());
        fillVariableWithValue ("Pt_j1j2", t_jet1jet2.Pt());
        fillVariableWithValue ("mDPhi_METJet2", fabs( t_MET.DeltaPhi ( t_jet2 )));
        fillVariableWithValue ("DR_Jet1Jet2"  , t_jet1.DeltaR( t_jet2 ));

        if ( n_jet_store >= 3 ){

          PFJet jet3 = c_pfjet_final -> GetConstituent<PFJet>(2);
          jetRawIndicesMap[2] = jet3.GetRawIndex();
          t_jet3.SetPtEtaPhiM ( jet3.Pt(), jet3.Eta(), jet3.Phi(), 0.0 );
          TLorentzVector t_jet1jet3 = t_jet1 + t_jet3;
          TLorentzVector t_jet2jet3 = t_jet2 + t_jet3;

          fillVariableWithValue ("M_j1j3", t_jet1jet3.M());
          fillVariableWithValue ("M_j2j3", t_jet2jet3.M());
          fillVariableWithValue ("mDPhi_METJet3", fabs( t_MET.DeltaPhi ( t_jet3 )));

          if ( n_jet_store >= 4 ){
            PFJet jet4 = c_pfjet_final -> GetConstituent<PFJet>(3);
            jetRawIndicesMap[3] = jet4.GetRawIndex();
            t_jet4.SetPtEtaPhiM ( jet4.Pt(), jet4.Eta(), jet4.Phi(), 0.0 );
            if ( n_jet_store >= 5 ){
              PFJet jet5 = c_pfjet_final -> GetConstituent<PFJet>(4);
              jetRawIndicesMap[4] = jet5.GetRawIndex();
              t_jet5.SetPtEtaPhiM ( jet5.Pt(), jet5.Eta(), jet5.Phi(), 0.0 );
              if ( n_jet_store >= 6 ){
                PFJet jet6 = c_pfjet_final -> GetConstituent<PFJet>(5);
                jetRawIndicesMap[5] = jet6.GetRawIndex();
                t_jet6.SetPtEtaPhiM ( jet6.Pt(), jet6.Eta(), jet6.Phi(), 0.0 );
              }
            }
          }
        }
      }
    }

    if ( n_ele_store >= 1 ) { 
      Electron ele1 = c_ele_final -> GetConstituent<Electron>(0);
      t_ele1.SetPtEtaPhiM ( ele1.Pt(), ele1.Eta(), ele1.Phi(), 0.0 );
      if ( reducedSkimType == 0 ) // for QCD skims, use the uncorrected Pt
        t_ele1.SetPtEtaPhiM ( ele1.Pt(), ele1.Eta(), ele1.Phi(), 0.0 );
      //if(!isData()) {
      //  Electron ele1smeared = *find(smearedEles.begin(), smearedEles.end(), ele1);
      //  Electron ele1scaledUp = *find(scaledUpEles.begin(), scaledUpEles.end(), ele1);
      //  Electron ele1scaledDown = *find(scaledDownEles.begin(), scaledDownEles.end(), ele1);
      //  t_ele1Smeared.SetPtEtaPhiM( ele1smeared.Pt(), ele1smeared.Eta(), ele1smeared.Phi(), 0.0);
      //  t_ele1ScaledUp.SetPtEtaPhiM( ele1scaledUp.Pt(), ele1scaledUp.Eta(), ele1scaledUp.Phi(), 0.0);
      //  t_ele1ScaledDown.SetPtEtaPhiM( ele1scaledDown.Pt(), ele1scaledDown.Eta(), ele1scaledDown.Phi(), 0.0);
      //}
      if ( n_ele_store >= 2 ) {
        Electron ele2 = c_ele_final -> GetConstituent<Electron>(1);
        t_ele2.SetPtEtaPhiM ( ele2.Pt(), ele2.Eta(), ele2.Phi(), 0.0 );
        if ( reducedSkimType == 0 ) // for QCD skims, use the uncorrected Pt
          t_ele2.SetPtEtaPhiM ( ele2.Pt(), ele2.Eta(), ele2.Phi(), 0.0 );
        //if(!isData()) {
        //  Electron ele2smeared = *find(smearedEles.begin(), smearedEles.end(), ele2);
        //  Electron ele2scaledUp = *find(scaledUpEles.begin(), scaledUpEles.end(), ele2);
        //  Electron ele2scaledDown = *find(scaledDownEles.begin(), scaledDownEles.end(), ele2);
        //  t_ele2Smeared.SetPtEtaPhiM( ele2smeared.Pt(), ele2smeared.Eta(), ele2smeared.Phi(), 0.0);
        //  t_ele2ScaledUp.SetPtEtaPhiM( ele2scaledUp.Pt(), ele2scaledUp.Eta(), ele2scaledUp.Phi(), 0.0);
        //  t_ele2ScaledDown.SetPtEtaPhiM( ele2scaledDown.Pt(), ele2scaledDown.Eta(), ele2scaledDown.Phi(), 0.0);
        //}

        TLorentzVector t_ele1ele2 = t_ele1 + t_ele2;
        //TLorentzVector t_ele1ele2_eer = t_ele1Smeared + t_ele2Smeared;
        //TLorentzVector t_ele1ele2_eesUp = t_ele1ScaledUp + t_ele2ScaledUp;
        //TLorentzVector t_ele1ele2_eesDown = t_ele1ScaledDown + t_ele2ScaledDown;
        fillVariableWithValue ("M_e1e2" , t_ele1ele2.M ());
        fillVariableWithValue ("Pt_e1e2", t_ele1ele2.Pt());
        //std::cout << "Pt_e1e2 Pt=" << t_ele1ele2.Pt() << std::endl;
        // systs
        //if(!isData()) {
        //  fillVariableWithValue ("M_e1e2_EER" , t_ele1ele2_eer.M ());
        //  fillVariableWithValue ("M_e1e2_EES_Up" , t_ele1ele2_eesUp.M ());
        //  fillVariableWithValue ("M_e1e2_EES_Dn" , t_ele1ele2_eesDown.M ());
        //  fillVariableWithValue ("Pt_e1e2_EER" , t_ele1ele2_eer.Pt ());
        //  fillVariableWithValue ("Pt_e1e2_EES_Up" , t_ele1ele2_eesUp.Pt ());
        //  fillVariableWithValue ("Pt_e1e2_EES_Dn" , t_ele1ele2_eesDown.Pt ());
        //  //std::cout << "Pt_e1e2_EER Pt=" << t_ele1ele2_eer.Pt() << std::endl;
        //  //std::cout << "Pt_e1e2_EES_Up Pt=" << t_ele1ele2_eesUp.Pt() << std::endl;
        //  //std::cout << "Pt_e1e2_EES_Down Pt=" << t_ele1ele2_eesDown.Pt() << std::endl;
        //}
      }
    } 

    if ( n_ele_store >= 1 ){

      double MT_Ele1MET = sqrt ( 2.0 * t_ele1.Pt() * t_MET.Pt() * ( 1.0 - cos ( t_MET.DeltaPhi(t_ele1))));

      TLorentzVector t_ele1MET = t_ele1 + t_MET;
      fillVariableWithValue("mDPhi_METEle1", fabs ( t_MET.DeltaPhi(t_ele1)));
      fillVariableWithValue("MT_Ele1MET"   , MT_Ele1MET); 
      fillVariableWithValue("Pt_Ele1MET"   , t_ele1MET.Pt());

      Electron ele1 = c_ele_final -> GetConstituent<Electron>(0);
      bool ele1PassHeep = ele1.PassHEEPID();
      int jetIdx = 0;
      if ( n_jet_store >= 1 ){ 
        int ele1LQMotherIdx = -1, jetLQMotherIdx = -1;
        if(ele1.MatchedGenParticleIdx() >= 0) {
          GenParticle matchedGenParticleEle = c_gen_all->GetConstituent<GenParticle>(ele1.MatchedGenParticleIdx());
          matchedGenParticleEle.PassUserID(GEN_ELE_FROM_LQ);
          ele1LQMotherIdx =  matchedGenParticleEle.MotherLQIndex();
        }
        jetLQMotherIdx = getJetLQMotherIdx(c_genJet_all, c_pfjet_final, jetRawIndicesMap[jetIdx], c_genQuark_hardProcess);

        TLorentzVector t_ele1jet1 = t_ele1 + t_jet1;
        fillVariableWithValue("DR_Ele1Jet1", t_ele1.DeltaR ( t_jet1 ));
        fillVariableWithValue("M_e1j1"     , t_ele1jet1.M());
        FillUserTH1D("deltaREle1Jets", t_ele1.DeltaR ( t_jet1 ));
        if(ele1LQMotherIdx == jetLQMotherIdx)
          FillUserTH1D("deltaR_Ele1JetSameLQ", t_ele1.DeltaR ( t_jet1 ));
        if(ele1PassHeep) {
          FillUserTH1D("deltaREleJetsPassHEEP", t_ele1.DeltaR ( t_jet1 ));
          if(isElectronBarrel(ele1.SCEta()))
            FillUserTH1D("deltaREleJetsBarrelPassHEEP", t_ele1.DeltaR ( t_jet1 ));
          else
            FillUserTH1D("deltaREleJetsEndcapPassHEEP", t_ele1.DeltaR ( t_jet1 ));
        }
        else {
          FillUserTH1D("deltaREleJetsFailHEEP", t_ele1.DeltaR ( t_jet1 ));
          if(ele1.PassHEEPGsfEleSCEtaMultiRangeCut()) {
            if(isElectronBarrel(ele1.SCEta()))
              FillUserTH1D("deltaREleJetsBarrelFailHEEP", t_ele1.DeltaR ( t_jet1 ));
            else
              FillUserTH1D("deltaREleJetsEndcapFailHEEP", t_ele1.DeltaR ( t_jet1 ));
          }
        }
        // systs
        //if(!isData()) {
        //  fillVariableWithValue("M_e1j1_EER"        , (t_ele1Smeared+t_jet1).M());
        //  fillVariableWithValue("M_e1j1_EES_Up"     , (t_ele1ScaledUp+t_jet1).M());
        //  fillVariableWithValue("M_e1j1_EES_Dn"     , (t_ele1ScaledDown+t_jet1).M());
        //  fillVariableWithValue("M_e1j1_JESTotal_Up", (t_ele1+t_jet1JESTotalUp).M());
        //  fillVariableWithValue("M_e1j1_JESTotal_Dn", (t_ele1+t_jet1JESTotalDown).M());
        //  fillVariableWithValue("M_e1j1_JER_Up"     , (t_ele1+t_jet1JERUp).M());
        //  fillVariableWithValue("M_e1j1_JER_Dn"     , (t_ele1+t_jet1JERDown).M());
        //}

        if ( n_jet_store >= 2 ){ 
          ++jetIdx;
          jetLQMotherIdx = getJetLQMotherIdx(c_genJet_all, c_pfjet_final, jetRawIndicesMap[jetIdx], c_genQuark_hardProcess);

          TLorentzVector t_ele1jet2 = t_ele1 + t_jet2;
          fillVariableWithValue("DR_Ele1Jet2", t_ele1.DeltaR ( t_jet2 ));
          fillVariableWithValue("M_e1j2"     , t_ele1jet2.M());
          if(ele1LQMotherIdx == jetLQMotherIdx)
            FillUserTH1D("deltaR_Ele1JetSameLQ", t_ele1.DeltaR ( t_jet2 ));
          FillUserTH1D("deltaREle1Jets", t_ele1.DeltaR ( t_jet2 ));
          if(ele1PassHeep) {
            FillUserTH1D("deltaREleJetsPassHEEP", t_ele1.DeltaR ( t_jet2 ));
            if(isElectronBarrel(ele1.SCEta()))
              FillUserTH1D("deltaREleJetsBarrelPassHEEP", t_ele1.DeltaR ( t_jet2 ));
            else
              FillUserTH1D("deltaREleJetsEndcapPassHEEP", t_ele1.DeltaR ( t_jet2 ));
          }
          else {
            FillUserTH1D("deltaREleJetsFailHEEP", t_ele1.DeltaR ( t_jet2 ));
            if(ele1.PassHEEPGsfEleSCEtaMultiRangeCut()) {
              if(isElectronBarrel(ele1.SCEta()))
                FillUserTH1D("deltaREleJetsBarrelFailHEEP", t_ele1.DeltaR ( t_jet2 ));
              else
                FillUserTH1D("deltaREleJetsEndcapFailHEEP", t_ele1.DeltaR ( t_jet2 ));
            }
          }
          // systs
          //if(!isData()) {
          //  fillVariableWithValue("M_e1j2_EER"        , (t_ele1Smeared+t_jet2).M());
          //  fillVariableWithValue("M_e1j2_EES_Up"     , (t_ele1ScaledUp+t_jet2).M());
          //  fillVariableWithValue("M_e1j2_EES_Dn"     , (t_ele1ScaledDown+t_jet2).M());
          //  fillVariableWithValue("M_e1j2_JESTotal_Up", (t_ele1+t_jet2JESTotalUp).M());
          //  fillVariableWithValue("M_e1j2_JESTotal_Dn", (t_ele1+t_jet2JESTotalDown).M());
          //  fillVariableWithValue("M_e1j2_JER_Up"     , (t_ele1+t_jet2JERUp).M());
          //  fillVariableWithValue("M_e1j2_JER_Dn"     , (t_ele1+t_jet2JERDown).M());
          //  fillVariableWithValue("sT_enujj"   , t_ele1.Pt() + t_MET.Pt() + t_jet1.Pt() + t_jet2.Pt());
          //}

          if ( n_jet_store >= 3 ){ 
            ++jetIdx;
            jetLQMotherIdx = getJetLQMotherIdx(c_genJet_all, c_pfjet_final, jetRawIndicesMap[jetIdx], c_genQuark_hardProcess);
            TLorentzVector t_ele1jet3 = t_ele1 + t_jet3;
            FillUserTH1D("deltaREle1Jets", t_ele1.DeltaR ( t_jet3 ));
            if(ele1LQMotherIdx == jetLQMotherIdx)
              FillUserTH1D("deltaR_Ele1JetSameLQ", t_ele1.DeltaR ( t_jet3 ));
            if ( n_jet_store >= 4 ){ 
              ++jetIdx;
              jetLQMotherIdx = getJetLQMotherIdx(c_genJet_all, c_pfjet_final, jetRawIndicesMap[jetIdx], c_genQuark_hardProcess);
              TLorentzVector t_ele1jet4 = t_ele1 + t_jet4;
              FillUserTH1D("deltaREle1Jets", t_ele1.DeltaR ( t_jet4 ));
              if(ele1LQMotherIdx == jetLQMotherIdx)
                FillUserTH1D("deltaR_Ele1JetSameLQ", t_ele1.DeltaR ( t_jet4 ));
              if ( n_jet_store >= 5 ){ 
                ++jetIdx;
                jetLQMotherIdx = getJetLQMotherIdx(c_genJet_all, c_pfjet_final, jetRawIndicesMap[jetIdx], c_genQuark_hardProcess);
                TLorentzVector t_ele1jet5 = t_ele1 + t_jet5;
                FillUserTH1D("deltaREle1Jets", t_ele1.DeltaR ( t_jet5 ));
                if(ele1LQMotherIdx == jetLQMotherIdx)
                  FillUserTH1D("deltaR_Ele1JetSameLQ", t_ele1.DeltaR ( t_jet5 ));
                if ( n_jet_store >= 6 ){ 
                  ++jetIdx;
                  jetLQMotherIdx = getJetLQMotherIdx(c_genJet_all, c_pfjet_final, jetRawIndicesMap[jetIdx], c_genQuark_hardProcess);
                  TLorentzVector t_ele1jet6 = t_ele1 + t_jet6;
                  FillUserTH1D("deltaREle1Jets", t_ele1.DeltaR ( t_jet6 ));
                  if(ele1LQMotherIdx == jetLQMotherIdx)
                    FillUserTH1D("deltaR_Ele1JetSameLQ", t_ele1.DeltaR ( t_jet6 ));
                }
              }
            }
          }
        }
      }
    }


    if ( n_ele_store >= 2 ){
      fillVariableWithValue("mDPhi_METEle2", fabs ( t_MET.DeltaPhi(t_ele2)));

      Electron ele2 = c_ele_final -> GetConstituent<Electron>(1);
      bool ele2PassLoose = ele2.PassEGammaIDLoose();
      bool ele2PassHeep = ele2.PassHEEPID();
      bool ele2PassHeepCaloIso = ele2.PassHEEPGsfEleEmHadD1IsoRhoCut();
      bool ele2PassEt = ele2.SCEt() > ele2_PtCut;
      int jetIdx = 0;
      if ( n_jet_store >= 1 ){ 
        PFJet closestJet = c_pfjet_final -> GetClosestInDR<PFJet, Electron>(ele2);
        int closestJetIndex = -1;
        if(closestJet.DeltaR(&ele2) < 0.3)
          closestJetIndex = closestJet.GetRawIndex();
        STDOUT("For ele2, closest jet has rawIndex=" << closestJetIndex);
        int ele2LQMotherIdx = -1, jetLQMotherIdx = -2;
        if(ele2.MatchedGenParticleIdx() >= 0) {
          GenParticle matchedGenParticleEle = c_gen_all->GetConstituent<GenParticle>(ele2.MatchedGenParticleIdx());
          //STDOUT("Electron2: matched gen particle Idx=" << ele2.MatchedGenParticleIdx() <<
          //    " is from LQ? " << matchedGenParticleEle.PassUserID(GEN_ELE_FROM_LQ));
          //STDOUT("\tindex of ele LQ mother: " << matchedGenParticleEle.MotherLQIndex());
          matchedGenParticleEle.PassUserID(GEN_ELE_FROM_LQ);
          ele2LQMotherIdx =  matchedGenParticleEle.MotherLQIndex();
          //STDOUT("\tindex of ele LQ mother: " << ele2LQMotherIdx);
        }
        //GenJet matchedGenJet = c_genJet_all->GetConstituent<GenJet>(c_pfjet_final->GetConstituent<PFJet>(jetRawIndicesMap[jetIdx]).MatchedGenJetIndex());
        //GenParticle matchedGenParticle = c_genQuark_hardProcess->GetClosestInDR<GenParticle, GenJet>(matchedGenJet);
        ////STDOUT("Jet1: matched genJet Idx=" << c_pfjet_final->GetConstituent<PFJet>(jetRawIndicesMap[0]).MatchedGenJetIndex());
        ////STDOUT("\tmatched genParticle Idx=" << matchedGenParticle.GetRawIndex() << "; is from LQ? " << matchedGenParticle.PassUserID(GEN_FROM_LQ));
        ////STDOUT("\tindex of LQ mother: " << matchedGenParticle.MotherLQIndex());
        //jetLQMotherIdx = matchedGenParticle.MotherLQIndex();
        jetLQMotherIdx = getJetLQMotherIdx(c_genJet_all, c_pfjet_final, jetRawIndicesMap[jetIdx], c_genQuark_hardProcess);
        //STDOUT("\tindex of jet LQ mother: " << jetLQMotherIdx);

        TLorentzVector t_ele2jet1 = t_ele2 + t_jet1;
        fillVariableWithValue("DR_Ele2Jet1", t_ele2.DeltaR ( t_jet1 ));
        fillVariableWithValue("M_e2j1"     , t_ele2jet1.M());
        if(ele2LQMotherIdx == jetLQMotherIdx) {
          FillUserTH1D("deltaR_Ele2JetSameLQ", t_ele2.DeltaR ( t_jet1 ));
          FillUserTH2D("deltaREle2JetVsPt_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
          if(ele2PassHeep)
            FillUserTH2D("deltaREle2JetVsPt_passHEEP_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
          else
            FillUserTH2D("deltaREle2JetVsPt_failHEEP_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
        }
        else if(jetLQMotherIdx > -1) {
          FillUserTH2D("deltaREle2JetVsPt_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
          if(ele2PassHeep)
            FillUserTH2D("deltaREle2JetVsPt_passHEEP_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
          else
            FillUserTH2D("deltaREle2JetVsPt_failHEEP_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
        }
        else {
          FillUserTH2D("deltaREle2JetVsPt_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
          if(ele2PassHeep)
            FillUserTH2D("deltaREle2JetVsPt_passHEEP_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
          else
            FillUserTH2D("deltaREle2JetVsPt_failHEEP_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
        }
        if(closestJetIndex != jetRawIndicesMap[jetIdx]) {
          if(ele2LQMotherIdx == jetLQMotherIdx) {
            FillUserTH1D("deltaR_Ele2JetSameLQ_excludeClosestJet", t_ele2.DeltaR ( t_jet1 ));
            FillUserTH2D("deltaREle2JetVsPt_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
            if(ele2PassHeep)
              FillUserTH2D("deltaREle2JetVsPt_passHEEP_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
            else
              FillUserTH2D("deltaREle2JetVsPt_failHEEP_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
          }
          else if(jetLQMotherIdx > -1) {
            FillUserTH2D("deltaREle2JetVsPt_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
            if(ele2PassHeep)
              FillUserTH2D("deltaREle2JetVsPt_passHEEP_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
            else
              FillUserTH2D("deltaREle2JetVsPt_failHEEP_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
          }
          else {
            FillUserTH2D("deltaREle2JetVsPt_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
            if(ele2PassHeep)
              FillUserTH2D("deltaREle2JetVsPt_passHEEP_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
            else
              FillUserTH2D("deltaREle2JetVsPt_failHEEP_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet1 ));
          }
        }
        if(ele2PassEt) {
          FillUserTH1D("deltaREle2Jets_highET", t_ele2.DeltaR ( t_jet1 ));
          if(closestJetIndex != jetRawIndicesMap[jetIdx]) {
            FillUserTH1D("deltaREle2Jets_highET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet1 ));
            if(ele2PassHeep)
              FillUserTH1D("deltaREle2PassHEEPJets_highET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet1 ));
            else if(!ele2PassHeepCaloIso) {
              FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet", t_ele2.DeltaR ( t_jet1 ));
              if(jetLQMotherIdx == ele2LQMotherIdx)
                FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_sameLQMother", t_ele2.DeltaR ( t_jet1 ));
              else if(jetLQMotherIdx > -1)
                FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_diffLQMother", t_ele2.DeltaR ( t_jet1 ));
              else
                FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_noLQMother", t_ele2.DeltaR ( t_jet1 ));
            }
          }
          if(isElectronBarrel(ele2.SCEta()))
            FillUserTH1D("deltaREle2Jets_highET_barrel", t_ele2.DeltaR ( t_jet1 ));
          else if(isElectronEndcap(ele2.SCEta()))
            FillUserTH1D("deltaREle2Jets_highET_endcap", t_ele2.DeltaR ( t_jet1 ));
        }
        else {
          FillUserTH1D("deltaREle2Jets_lowET", t_ele2.DeltaR ( t_jet1 ));
          if(closestJetIndex != jetRawIndicesMap[jetIdx]) {
            FillUserTH1D("deltaREle2Jets_lowET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet1 ));
            if(ele2PassHeep)
              FillUserTH1D("deltaREle2PassHEEPJets_lowET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet1 ));
            else if(!ele2PassHeepCaloIso) {
              FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet", t_ele2.DeltaR ( t_jet1 ));
              if(jetLQMotherIdx == ele2LQMotherIdx)
                FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_sameLQMother", t_ele2.DeltaR ( t_jet1 ));
              else if(jetLQMotherIdx > -1)
                FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_diffLQMother", t_ele2.DeltaR ( t_jet1 ));
              else
                FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_noLQMother", t_ele2.DeltaR ( t_jet1 ));
            }
            if(!ele2PassHeep) {
              if(t_ele2.DeltaR ( t_jet1 ) < 0.4 && ele2.SCEt() > 35 && ele2PassLoose)
                STDOUT("--> found event with lowET ele2 failing heep; jet1 is the second-closest jet with deltaR = " << t_ele2.DeltaR ( t_jet1 ));
            }
          }
          if(isElectronBarrel(ele2.SCEta()))
            FillUserTH1D("deltaREle2Jets_lowET_barrel", t_ele2.DeltaR ( t_jet1 ));
          else if(isElectronEndcap(ele2.SCEta()))
            FillUserTH1D("deltaREle2Jets_lowET_endcap", t_ele2.DeltaR ( t_jet1 ));
        }
        if(closestJetIndex != jetRawIndicesMap[jetIdx]) {
          FillUserTH1D("deltaREle2JetsExcludeClosestJet", t_ele2.DeltaR ( t_jet1 ));
          if(ele2PassHeep)
            FillUserTH1D("deltaREle2PassHEEPJetsExcludeClosestJet", t_ele2.DeltaR ( t_jet1 ));
          else if(!ele2PassHeepCaloIso)
            FillUserTH1D("deltaREle2FailHEEPCaloIsoJetsExcludeClosestJet", t_ele2.DeltaR ( t_jet1 ));
        }
        FillUserTH1D("deltaREle2Jets", t_ele2.DeltaR ( t_jet1 ));
        if(ele2PassHeep) {
          FillUserTH1D("deltaREleJetsPassHEEP", t_ele1.DeltaR ( t_jet1 ));
          if(isElectronBarrel(ele2.SCEta()))
            FillUserTH1D("deltaREleJetsBarrelPassHEEP", t_ele1.DeltaR ( t_jet1 ));
          else
            FillUserTH1D("deltaREleJetsEndcapPassHEEP", t_ele1.DeltaR ( t_jet1 ));
        }
        else {
          FillUserTH1D("deltaREleJetsFailHEEP", t_ele1.DeltaR ( t_jet1 ));
          if(ele2.PassHEEPGsfEleSCEtaMultiRangeCut()) {
            if(isElectronBarrel(ele2.SCEta()))
              FillUserTH1D("deltaREleJetsBarrelFailHEEP", t_ele1.DeltaR ( t_jet1 ));
            else
              FillUserTH1D("deltaREleJetsEndcapFailHEEP", t_ele1.DeltaR ( t_jet1 ));
          }
        }
        // systs
        //if(!isData()) {
        //  fillVariableWithValue("M_e2j1_EER"        , (t_ele2Smeared+t_jet1).M());
        //  fillVariableWithValue("M_e2j1_EES_Up"     , (t_ele2ScaledUp+t_jet1).M());
        //  fillVariableWithValue("M_e2j1_EES_Dn"     , (t_ele2ScaledDown+t_jet1).M());
        //  fillVariableWithValue("M_e2j1_JESTotal_Up", (t_ele2+t_jet1JESTotalUp).M());
        //  fillVariableWithValue("M_e2j1_JESTotal_Dn", (t_ele2+t_jet1JESTotalDown).M());
        //  fillVariableWithValue("M_e2j1_JER_Up"     , (t_ele2+t_jet1JERUp).M());
        //  fillVariableWithValue("M_e2j1_JER_Dn"     , (t_ele2+t_jet1JERDown).M());
        //}

        if ( n_jet_store >= 2 ){ 
          ++jetIdx;
          jetLQMotherIdx = getJetLQMotherIdx(c_genJet_all, c_pfjet_final, jetRawIndicesMap[jetIdx], c_genQuark_hardProcess);

          TLorentzVector t_ele2jet2 = t_ele2 + t_jet2;
          if(ele2LQMotherIdx == jetLQMotherIdx) {
            FillUserTH1D("deltaR_Ele2JetSameLQ", t_ele2.DeltaR ( t_jet2 ));
            FillUserTH2D("deltaREle2JetVsPt_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
            if(ele2PassHeep)
              FillUserTH2D("deltaREle2JetVsPt_passHEEP_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
            else
              FillUserTH2D("deltaREle2JetVsPt_failHEEP_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
          }
          else if(jetLQMotherIdx > -1) {
            FillUserTH2D("deltaREle2JetVsPt_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
            if(ele2PassHeep)
              FillUserTH2D("deltaREle2JetVsPt_passHEEP_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
            else
              FillUserTH2D("deltaREle2JetVsPt_failHEEP_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
          }
          else {
            FillUserTH2D("deltaREle2JetVsPt_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
            if(ele2PassHeep)
              FillUserTH2D("deltaREle2JetVsPt_passHEEP_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
            else
              FillUserTH2D("deltaREle2JetVsPt_failHEEP_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
          }
          if(closestJetIndex != jetRawIndicesMap[jetIdx]) {
            if(ele2LQMotherIdx == jetLQMotherIdx) {
              FillUserTH1D("deltaR_Ele2JetSameLQ_excludeClosestJet", t_ele2.DeltaR ( t_jet2 ));
              FillUserTH2D("deltaREle2JetVsPt_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
              if(ele2PassHeep)
                FillUserTH2D("deltaREle2JetVsPt_passHEEP_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
              else
                FillUserTH2D("deltaREle2JetVsPt_failHEEP_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
            }
            else if(jetLQMotherIdx > -1) {
              FillUserTH2D("deltaREle2JetVsPt_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
              if(ele2PassHeep)
                FillUserTH2D("deltaREle2JetVsPt_passHEEP_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
              else
                FillUserTH2D("deltaREle2JetVsPt_failHEEP_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
            }
            else {
              FillUserTH2D("deltaREle2JetVsPt_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
              if(ele2PassHeep)
                FillUserTH2D("deltaREle2JetVsPt_passHEEP_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
              else
                FillUserTH2D("deltaREle2JetVsPt_failHEEP_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet2 ));
            }
          }
          if(ele2PassEt) {
            FillUserTH1D("deltaREle2Jets_highET", t_ele2.DeltaR ( t_jet2 ));
            if(closestJetIndex != jetRawIndicesMap[1]) {
              FillUserTH1D("deltaREle2Jets_highET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet2 ));
              if(ele2PassHeep)
                FillUserTH1D("deltaREle2PassHEEPJets_highET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet2 ));
              else if(!ele2PassHeepCaloIso) {
                FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet", t_ele2.DeltaR ( t_jet2 ));
                if(jetLQMotherIdx == ele2LQMotherIdx)
                  FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_sameLQMother", t_ele2.DeltaR ( t_jet2 ));
                else if(jetLQMotherIdx > -1)
                  FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_diffLQMother", t_ele2.DeltaR ( t_jet2 ));
                else
                  FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_noLQMother", t_ele2.DeltaR ( t_jet2 ));
              }
            }
            if(isElectronBarrel(ele2.SCEta()))
              FillUserTH1D("deltaREle2Jets_highET_barrel", t_ele2.DeltaR ( t_jet2 ));
            else if(isElectronEndcap(ele2.SCEta()))
              FillUserTH1D("deltaREle2Jets_highET_endcap", t_ele2.DeltaR ( t_jet2 ));
          }
          else {
            FillUserTH1D("deltaREle2Jets_lowET", t_ele2.DeltaR ( t_jet2 ));
            if(closestJetIndex != jetRawIndicesMap[1]) {
              FillUserTH1D("deltaREle2Jets_lowET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet2 ));
              if(ele2PassHeep)
                FillUserTH1D("deltaREle2PassHEEPJets_lowET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet2 ));
              else if(!ele2PassHeepCaloIso) {
                FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet", t_ele2.DeltaR ( t_jet2 ));
                if(jetLQMotherIdx == ele2LQMotherIdx)
                  FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_sameLQMother", t_ele2.DeltaR ( t_jet2 ));
                else if(jetLQMotherIdx > -1)
                  FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_diffLQMother", t_ele2.DeltaR ( t_jet2 ));
                else
                  FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_noLQMother", t_ele2.DeltaR ( t_jet2 ));
              }
              if(!ele2PassHeep) {
                if(t_ele2.DeltaR ( t_jet2 ) < 0.4 && ele2.SCEt() > 35 && ele2PassLoose)
                  STDOUT("--> found event with lowET ele2 failing heep; jet2 is the second-closest jet with deltaR = " << t_ele2.DeltaR ( t_jet2 ));
              }
            }
            if(isElectronBarrel(ele2.SCEta()))
              FillUserTH1D("deltaREle2Jets_lowET_barrel", t_ele2.DeltaR ( t_jet2 ));
            else if(isElectronEndcap(ele2.SCEta()))
              FillUserTH1D("deltaREle2Jets_lowET_endcap", t_ele2.DeltaR ( t_jet2 ));
          }
          fillVariableWithValue("DR_Ele2Jet2", t_ele2.DeltaR ( t_jet2 ));
          fillVariableWithValue("M_e2j2"     , t_ele2jet2.M());
          fillVariableWithValue("sT_eejj"    , t_ele1.Pt() + t_ele2.Pt() + t_jet1.Pt() + t_jet2.Pt());
          if(closestJetIndex != jetRawIndicesMap[1]) {
            FillUserTH1D("deltaREle2JetsExcludeClosestJet", t_ele2.DeltaR ( t_jet2 ));
            if(ele2PassHeep)
              FillUserTH1D("deltaREle2PassHEEPJetsExcludeClosestJet", t_ele2.DeltaR ( t_jet2 ));
            else if(!ele2PassHeepCaloIso)
              FillUserTH1D("deltaREle2FailHEEPCaloIsoJetsExcludeClosestJet", t_ele2.DeltaR ( t_jet2 ));
          }
          FillUserTH1D("deltaREle2Jets", t_ele2.DeltaR ( t_jet2 ));
          if(ele2PassHeep) {
            FillUserTH1D("deltaREleJetsPassHEEP", t_ele1.DeltaR ( t_jet2 ));
            if(isElectronBarrel(ele2.SCEta()))
              FillUserTH1D("deltaREleJetsBarrelPassHEEP", t_ele1.DeltaR ( t_jet2 ));
            else
              FillUserTH1D("deltaREleJetsEndcapPassHEEP", t_ele1.DeltaR ( t_jet2 ));
          }
          else {
            FillUserTH1D("deltaREleJetsFailHEEP", t_ele1.DeltaR ( t_jet2 ));
            if(ele2.PassHEEPGsfEleSCEtaMultiRangeCut()) {
              if(isElectronBarrel(ele2.SCEta()))
                FillUserTH1D("deltaREleJetsBarrelFailHEEP", t_ele1.DeltaR ( t_jet2 ));
              else
                FillUserTH1D("deltaREleJetsEndcapFailHEEP", t_ele1.DeltaR ( t_jet2 ));
            }
          }
          // systs
          //if(!isData()) {
          //  fillVariableWithValue("M_e2j2_EER"        , (t_ele2Smeared+t_jet2).M());
          //  fillVariableWithValue("M_e2j2_EES_Up"     , (t_ele2ScaledUp+t_jet2).M());
          //  fillVariableWithValue("M_e2j2_EES_Dn"     , (t_ele2ScaledDown+t_jet2).M());
          //  fillVariableWithValue("M_e2j2_JESTotal_Up", (t_ele2+t_jet2JESTotalUp).M());
          //  fillVariableWithValue("M_e2j2_JESTotal_Dn", (t_ele2+t_jet2JESTotalDown).M());
          //  fillVariableWithValue("M_e2j2_JER_Up"     , (t_ele2+t_jet2JERUp).M());
          //  fillVariableWithValue("M_e2j2_JER_Dn"     , (t_ele2+t_jet2JERDown).M());
          //  fillVariableWithValue("sT_eejj_EER"    , t_ele1Smeared.Pt() + t_ele2Smeared.Pt() + t_jet1.Pt() + t_jet2.Pt());
          //  fillVariableWithValue("sT_eejj_EES_Up"    , t_ele1ScaledUp.Pt() + t_ele2ScaledUp.Pt() + t_jet1.Pt() + t_jet2.Pt());
          //  fillVariableWithValue("sT_eejj_EES_Dn"    , t_ele1ScaledDown.Pt() + t_ele2ScaledDown.Pt() + t_jet1.Pt() + t_jet2.Pt());
          //  fillVariableWithValue("sT_eejj_JESTotal_Up"    , t_ele1.Pt() + t_ele2.Pt() + t_jet1JESTotalUp.Pt() + t_jet2JESTotalUp.Pt());
          //  fillVariableWithValue("sT_eejj_JESTotal_Dn"    , t_ele1.Pt() + t_ele2.Pt() + t_jet1JESTotalDown.Pt() + t_jet2JESTotalDown.Pt());
          //  fillVariableWithValue("sT_eejj_JER_Up"    , t_ele1.Pt() + t_ele2.Pt() + t_jet1JERUp.Pt() + t_jet2JERUp.Pt());
          //  fillVariableWithValue("sT_eejj_JER_Dn"    , t_ele1.Pt() + t_ele2.Pt() + t_jet1JERDown.Pt() + t_jet2JERDown.Pt());
          //}
          if ( n_jet_store >= 3 ){ 
            ++jetIdx;
            jetLQMotherIdx = getJetLQMotherIdx(c_genJet_all, c_pfjet_final, jetRawIndicesMap[jetIdx], c_genQuark_hardProcess);
            TLorentzVector t_ele2jet3 = t_ele2 + t_jet3;
            if(ele2LQMotherIdx == jetLQMotherIdx) {
              FillUserTH1D("deltaR_Ele2JetSameLQ", t_ele2.DeltaR ( t_jet3 ));
              FillUserTH2D("deltaREle2JetVsPt_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
              if(ele2PassHeep)
                FillUserTH2D("deltaREle2JetVsPt_passHEEP_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
              else
                FillUserTH2D("deltaREle2JetVsPt_failHEEP_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
            }
            else if(jetLQMotherIdx > -1) {
              FillUserTH2D("deltaREle2JetVsPt_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
              if(ele2PassHeep)
                FillUserTH2D("deltaREle2JetVsPt_passHEEP_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
              else
                FillUserTH2D("deltaREle2JetVsPt_failHEEP_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
            }
            else {
              FillUserTH2D("deltaREle2JetVsPt_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
              if(ele2PassHeep)
                FillUserTH2D("deltaREle2JetVsPt_passHEEP_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
              else
                FillUserTH2D("deltaREle2JetVsPt_failHEEP_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
            }
            if(closestJetIndex != jetRawIndicesMap[jetIdx]) {
              if(ele2LQMotherIdx == jetLQMotherIdx) {
                FillUserTH1D("deltaR_Ele2JetSameLQ_excludeClosestJet", t_ele2.DeltaR ( t_jet3 ));
                FillUserTH2D("deltaREle2JetVsPt_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
                if(ele2PassHeep)
                  FillUserTH2D("deltaREle2JetVsPt_passHEEP_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
                else
                  FillUserTH2D("deltaREle2JetVsPt_failHEEP_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
              }
              else if(jetLQMotherIdx > -1) {
                FillUserTH2D("deltaREle2JetVsPt_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
                if(ele2PassHeep)
                  FillUserTH2D("deltaREle2JetVsPt_passHEEP_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
                else
                  FillUserTH2D("deltaREle2JetVsPt_failHEEP_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
              }
              else {
                FillUserTH2D("deltaREle2JetVsPt_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
                if(ele2PassHeep)
                  FillUserTH2D("deltaREle2JetVsPt_passHEEP_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
                else
                  FillUserTH2D("deltaREle2JetVsPt_failHEEP_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet3 ));
              }
            }
            if(ele2PassEt) {
              FillUserTH1D("deltaREle2Jets_highET", t_ele2.DeltaR ( t_jet3 ));
              if(closestJetIndex != jetRawIndicesMap[2]) {
                FillUserTH1D("deltaREle2Jets_highET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet3 ));
                if(ele2PassHeep)
                  FillUserTH1D("deltaREle2PassHEEPJets_highET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet3 ));
                else if(!ele2PassHeepCaloIso) {
                  FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet", t_ele2.DeltaR ( t_jet3 ));
                  if(jetLQMotherIdx == ele2LQMotherIdx)
                    FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_sameLQMother", t_ele2.DeltaR ( t_jet3 ));
                  else if(jetLQMotherIdx > -1)
                    FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_diffLQMother", t_ele2.DeltaR ( t_jet3 ));
                  else
                    FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_noLQMother", t_ele2.DeltaR ( t_jet3 ));
                }
              }
              if(isElectronBarrel(ele2.SCEta()))
                FillUserTH1D("deltaREle2Jets_highET_barrel", t_ele2.DeltaR ( t_jet3 ));
              else if(isElectronEndcap(ele2.SCEta()))
                FillUserTH1D("deltaREle2Jets_highET_endcap", t_ele2.DeltaR ( t_jet3 ));
            }
            else {
              FillUserTH1D("deltaREle2Jets_lowET", t_ele2.DeltaR ( t_jet3 ));
              if(closestJetIndex != jetRawIndicesMap[2]) {
                FillUserTH1D("deltaREle2Jets_lowET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet3 ));
                if(ele2PassHeep)
                  FillUserTH1D("deltaREle2PassHEEPJets_lowET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet3 ));
                else if(!ele2PassHeepCaloIso) {
                  FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet", t_ele2.DeltaR ( t_jet3 ));
                  if(jetLQMotherIdx == ele2LQMotherIdx)
                    FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_sameLQMother", t_ele2.DeltaR ( t_jet3 ));
                  else if(jetLQMotherIdx > -1)
                    FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_diffLQMother", t_ele2.DeltaR ( t_jet3 ));
                  else
                    FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_noLQMother", t_ele2.DeltaR ( t_jet3 ));
                }
                if(!ele2PassHeep) {
                  if(t_ele2.DeltaR ( t_jet3 ) < 0.4 && ele2.SCEt() > 35 && ele2PassLoose)
                    STDOUT("--> found event with lowET ele2 failing heep; jet3 is the second-closest jet with deltaR = " << t_ele2.DeltaR ( t_jet3 ));
                }
              }
              if(isElectronBarrel(ele2.SCEta()))
                FillUserTH1D("deltaREle2Jets_lowET_barrel", t_ele2.DeltaR ( t_jet3 ));
              else if(isElectronEndcap(ele2.SCEta()))
                FillUserTH1D("deltaREle2Jets_lowET_endcap", t_ele2.DeltaR ( t_jet3 ));
            }
            if(closestJetIndex != jetRawIndicesMap[2]) {
              FillUserTH1D("deltaREle2JetsExcludeClosestJet", t_ele2.DeltaR ( t_jet3 ));
              if(ele2PassHeep)
                FillUserTH1D("deltaREle2PassHEEPJetsExcludeClosestJet", t_ele2.DeltaR ( t_jet3 ));
              else if(!ele2PassHeepCaloIso)
                FillUserTH1D("deltaREle2FailHEEPCaloIsoJetsExcludeClosestJet", t_ele2.DeltaR ( t_jet3 ));
            }
            FillUserTH1D("deltaREle2Jets", t_ele2.DeltaR ( t_jet3 ));
            if ( n_jet_store >= 4 ){ 
              ++jetIdx;
              jetLQMotherIdx = getJetLQMotherIdx(c_genJet_all, c_pfjet_final, jetRawIndicesMap[jetIdx], c_genQuark_hardProcess);
              if(ele2LQMotherIdx == jetLQMotherIdx) {
                FillUserTH1D("deltaR_Ele2JetSameLQ", t_ele2.DeltaR ( t_jet4 ));
                FillUserTH2D("deltaREle2JetVsPt_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
                if(ele2PassHeep)
                  FillUserTH2D("deltaREle2JetVsPt_passHEEP_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
                else
                  FillUserTH2D("deltaREle2JetVsPt_failHEEP_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
              }
              else if(jetLQMotherIdx > -1) {
                FillUserTH2D("deltaREle2JetVsPt_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
                if(ele2PassHeep)
                  FillUserTH2D("deltaREle2JetVsPt_passHEEP_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
                else
                  FillUserTH2D("deltaREle2JetVsPt_failHEEP_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
              }
              else {
                FillUserTH2D("deltaREle2JetVsPt_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
                if(ele2PassHeep)
                  FillUserTH2D("deltaREle2JetVsPt_passHEEP_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
                else
                  FillUserTH2D("deltaREle2JetVsPt_failHEEP_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
              }
              if(closestJetIndex != jetRawIndicesMap[jetIdx]) {
                if(ele2LQMotherIdx == jetLQMotherIdx) {
                  FillUserTH1D("deltaR_Ele2JetSameLQ_excludeClosestJet", t_ele2.DeltaR ( t_jet4 ));
                  FillUserTH2D("deltaREle2JetVsPt_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
                  if(ele2PassHeep)
                    FillUserTH2D("deltaREle2JetVsPt_passHEEP_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
                  else
                    FillUserTH2D("deltaREle2JetVsPt_failHEEP_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
                }
                else if(jetLQMotherIdx > -1) {
                  FillUserTH2D("deltaREle2JetVsPt_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
                  if(ele2PassHeep)
                    FillUserTH2D("deltaREle2JetVsPt_passHEEP_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
                  else
                    FillUserTH2D("deltaREle2JetVsPt_failHEEP_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
                }
                else {
                  FillUserTH2D("deltaREle2JetVsPt_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
                  if(ele2PassHeep)
                    FillUserTH2D("deltaREle2JetVsPt_passHEEP_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
                  else
                    FillUserTH2D("deltaREle2JetVsPt_failHEEP_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet4 ));
                }
              }
              if(ele2PassEt) {
                FillUserTH1D("deltaREle2Jets_highET", t_ele2.DeltaR ( t_jet4 ));
                if(closestJetIndex != jetRawIndicesMap[3]) {
                  FillUserTH1D("deltaREle2Jets_highET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet4 ));
                  if(ele2PassHeep)
                    FillUserTH1D("deltaREle2PassHEEPJets_highET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet4 ));
                  else if(!ele2PassHeepCaloIso) {
                    FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet", t_ele2.DeltaR ( t_jet4 ));
                    if(jetLQMotherIdx == ele2LQMotherIdx)
                      FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_sameLQMother", t_ele2.DeltaR ( t_jet4 ));
                    else if(jetLQMotherIdx > -1)
                      FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_diffLQMother", t_ele2.DeltaR ( t_jet4 ));
                    else
                      FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_noLQMother", t_ele2.DeltaR ( t_jet4 ));
                  }
                }
                if(isElectronBarrel(ele2.SCEta()))
                  FillUserTH1D("deltaREle2Jets_highET_barrel", t_ele2.DeltaR ( t_jet4 ));
                else if(isElectronEndcap(ele2.SCEta()))
                  FillUserTH1D("deltaREle2Jets_highET_endcap", t_ele2.DeltaR ( t_jet4 ));
              }
              else {
                FillUserTH1D("deltaREle2Jets_lowET", t_ele2.DeltaR ( t_jet4 ));
                if(closestJetIndex != jetRawIndicesMap[3]) {
                  FillUserTH1D("deltaREle2Jets_lowET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet4 ));
                  if(ele2PassHeep)
                    FillUserTH1D("deltaREle2PassHEEPJets_lowET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet4 ));
                  else if(!ele2PassHeepCaloIso) {
                    FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet", t_ele2.DeltaR ( t_jet4 ));
                    if(jetLQMotherIdx == ele2LQMotherIdx)
                      FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_sameLQMother", t_ele2.DeltaR ( t_jet4 ));
                    else if(jetLQMotherIdx > -1)
                      FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_diffLQMother", t_ele2.DeltaR ( t_jet4 ));
                    else
                      FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_noLQMother", t_ele2.DeltaR ( t_jet4 ));
                  }
                  if(!ele2PassHeep) {
                    if(t_ele2.DeltaR ( t_jet4 ) < 0.4 && ele2.SCEt() > 35 && ele2PassLoose)
                      STDOUT("--> found event with lowET ele2 failing heep; jet4 is the second-closest jet with deltaR = " << t_ele2.DeltaR ( t_jet4 ));
                  }
                }
                if(isElectronBarrel(ele2.SCEta()))
                  FillUserTH1D("deltaREle2Jets_lowET_barrel", t_ele2.DeltaR ( t_jet4 ));
                else if(isElectronEndcap(ele2.SCEta()))
                  FillUserTH1D("deltaREle2Jets_lowET_endcap", t_ele2.DeltaR ( t_jet4 ));
              }
              TLorentzVector t_ele2jet4 = t_ele2 + t_jet4;
              if(closestJetIndex != jetRawIndicesMap[3]) {
                FillUserTH1D("deltaREle2JetsExcludeClosestJet", t_ele2.DeltaR ( t_jet4 ));
                if(ele2PassHeep)
                  FillUserTH1D("deltaREle2PassHEEPJetsExcludeClosestJet", t_ele2.DeltaR ( t_jet4 ));
                else if(!ele2PassHeepCaloIso)
                  FillUserTH1D("deltaREle2FailHEEPCaloIsoJetsExcludeClosestJet", t_ele2.DeltaR ( t_jet4 ));
              }
              FillUserTH1D("deltaREle2Jets", t_ele2.DeltaR ( t_jet4 ));
              if ( n_jet_store >= 5 ){ 
                ++jetIdx;
                jetLQMotherIdx = getJetLQMotherIdx(c_genJet_all, c_pfjet_final, jetRawIndicesMap[jetIdx], c_genQuark_hardProcess);
                TLorentzVector t_ele2jet5 = t_ele2 + t_jet5;
                if(ele2LQMotherIdx == jetLQMotherIdx) {
                  FillUserTH1D("deltaR_Ele2JetSameLQ", t_ele2.DeltaR ( t_jet5 ));
                  FillUserTH2D("deltaREle2JetVsPt_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                  if(ele2PassHeep)
                    FillUserTH2D("deltaREle2JetVsPt_passHEEP_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                  else
                    FillUserTH2D("deltaREle2JetVsPt_failHEEP_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                }
                else if(jetLQMotherIdx > -1) {
                  FillUserTH2D("deltaREle2JetVsPt_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                  if(ele2PassHeep)
                    FillUserTH2D("deltaREle2JetVsPt_passHEEP_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                  else
                    FillUserTH2D("deltaREle2JetVsPt_failHEEP_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                }
                else {
                  FillUserTH2D("deltaREle2JetVsPt_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                  if(ele2PassHeep)
                    FillUserTH2D("deltaREle2JetVsPt_passHEEP_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                  else
                    FillUserTH2D("deltaREle2JetVsPt_failHEEP_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                }
                if(closestJetIndex != jetRawIndicesMap[jetIdx]) {
                  if(ele2LQMotherIdx == jetLQMotherIdx) {
                    FillUserTH1D("deltaR_Ele2JetSameLQ_excludeClosestJet", t_ele2.DeltaR ( t_jet5 ));
                    FillUserTH2D("deltaREle2JetVsPt_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                    if(ele2PassHeep)
                      FillUserTH2D("deltaREle2JetVsPt_passHEEP_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                    else
                      FillUserTH2D("deltaREle2JetVsPt_failHEEP_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                  }
                  else if(jetLQMotherIdx > -1) {
                    FillUserTH2D("deltaREle2JetVsPt_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                    if(ele2PassHeep)
                      FillUserTH2D("deltaREle2JetVsPt_passHEEP_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                    else
                      FillUserTH2D("deltaREle2JetVsPt_failHEEP_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                  }
                  else {
                    FillUserTH2D("deltaREle2JetVsPt_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                    if(ele2PassHeep)
                      FillUserTH2D("deltaREle2JetVsPt_passHEEP_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                    else
                      FillUserTH2D("deltaREle2JetVsPt_failHEEP_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet5 ));
                  }
                }
                if(ele2PassEt) {
                  FillUserTH1D("deltaREle2Jets_highET", t_ele2.DeltaR ( t_jet5 ));
                  if(closestJetIndex != jetRawIndicesMap[4]) {
                    FillUserTH1D("deltaREle2Jets_highET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet5 ));
                    if(ele2PassHeep)
                      FillUserTH1D("deltaREle2PassHEEPJets_highET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet5 ));
                    else if(!ele2PassHeepCaloIso) {
                      FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet", t_ele2.DeltaR ( t_jet5 ));
                      if(jetLQMotherIdx == ele2LQMotherIdx)
                        FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_sameLQMother", t_ele2.DeltaR ( t_jet5 ));
                      else if(jetLQMotherIdx > -1)
                        FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_diffLQMother", t_ele2.DeltaR ( t_jet5 ));
                      else
                        FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_noLQMother", t_ele2.DeltaR ( t_jet5 ));
                    }
                  }
                  if(isElectronBarrel(ele2.SCEta()))
                    FillUserTH1D("deltaREle2Jets_highET_barrel", t_ele2.DeltaR ( t_jet5 ));
                  else if(isElectronEndcap(ele2.SCEta()))
                    FillUserTH1D("deltaREle2Jets_highET_endcap", t_ele2.DeltaR ( t_jet5 ));
                }
                else {
                  FillUserTH1D("deltaREle2Jets_lowET", t_ele2.DeltaR ( t_jet5 ));
                  if(closestJetIndex != jetRawIndicesMap[4]) {
                    FillUserTH1D("deltaREle2Jets_lowET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet5 ));
                    if(ele2PassHeep)
                      FillUserTH1D("deltaREle2PassHEEPJets_lowET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet5 ));
                    else if(!ele2PassHeepCaloIso) {
                      FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet", t_ele2.DeltaR ( t_jet5 ));
                      if(jetLQMotherIdx == ele2LQMotherIdx)
                        FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_sameLQMother", t_ele2.DeltaR ( t_jet5 ));
                      else if(jetLQMotherIdx > -1)
                        FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_diffLQMother", t_ele2.DeltaR ( t_jet5 ));
                      else
                        FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_noLQMother", t_ele2.DeltaR ( t_jet5 ));
                    }
                    if(!ele2PassHeep) {
                      if(t_ele2.DeltaR ( t_jet5 ) < 0.4 && ele2.SCEt() > 35 && ele2PassLoose)
                        STDOUT("--> found event with lowET ele2 failing heep; jet5 is the second-closest jet with deltaR = " << t_ele2.DeltaR ( t_jet5 ));
                    }
                  }
                  if(isElectronBarrel(ele2.SCEta()))
                    FillUserTH1D("deltaREle2Jets_lowET_barrel", t_ele2.DeltaR ( t_jet5 ));
                  else if(isElectronEndcap(ele2.SCEta()))
                    FillUserTH1D("deltaREle2Jets_lowET_endcap", t_ele2.DeltaR ( t_jet5 ));
                }
                if(closestJetIndex != jetRawIndicesMap[4]) {
                  FillUserTH1D("deltaREle2JetsExcludeClosestJet", t_ele2.DeltaR ( t_jet5 ));
                  if(ele2PassHeep)
                    FillUserTH1D("deltaREle2PassHEEPJetsExcludeClosestJet", t_ele2.DeltaR ( t_jet5 ));
                  else if(!ele2PassHeepCaloIso)
                    FillUserTH1D("deltaREle2FailHEEPCaloIsoJetsExcludeClosestJet", t_ele2.DeltaR ( t_jet5 ));
                }
                FillUserTH1D("deltaREle2Jets", t_ele2.DeltaR ( t_jet5 ));
                if ( n_jet_store >= 6 ){ 
                  ++jetIdx;
                  jetLQMotherIdx = getJetLQMotherIdx(c_genJet_all, c_pfjet_final, jetRawIndicesMap[jetIdx], c_genQuark_hardProcess);
                  TLorentzVector t_ele2jet6 = t_ele2 + t_jet6;
                  if(ele2LQMotherIdx == jetLQMotherIdx) {
                    FillUserTH1D("deltaR_Ele2JetSameLQ", t_ele2.DeltaR ( t_jet6 ));
                    FillUserTH2D("deltaREle2JetVsPt_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                    if(ele2PassHeep)
                      FillUserTH2D("deltaREle2JetVsPt_passHEEP_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                    else
                      FillUserTH2D("deltaREle2JetVsPt_failHEEP_sameLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                  }
                  else if(jetLQMotherIdx > -1) {
                    FillUserTH2D("deltaREle2JetVsPt_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                    if(ele2PassHeep)
                      FillUserTH2D("deltaREle2JetVsPt_passHEEP_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                    else
                      FillUserTH2D("deltaREle2JetVsPt_failHEEP_diffLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                  }
                  else {
                    FillUserTH2D("deltaREle2JetVsPt_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                    if(ele2PassHeep)
                      FillUserTH2D("deltaREle2JetVsPt_passHEEP_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                    else
                      FillUserTH2D("deltaREle2JetVsPt_failHEEP_noLQmother", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                  }
                  if(closestJetIndex != jetRawIndicesMap[jetIdx]) {
                    if(ele2LQMotherIdx == jetLQMotherIdx) {
                      FillUserTH1D("deltaR_Ele2JetSameLQ_excludeClosestJet", t_ele2.DeltaR ( t_jet6 ));
                      FillUserTH2D("deltaREle2JetVsPt_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                      if(ele2PassHeep)
                        FillUserTH2D("deltaREle2JetVsPt_passHEEP_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                      else
                        FillUserTH2D("deltaREle2JetVsPt_failHEEP_sameLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                    }
                    else if(jetLQMotherIdx > -1) {
                      FillUserTH2D("deltaREle2JetVsPt_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                      if(ele2PassHeep)
                        FillUserTH2D("deltaREle2JetVsPt_passHEEP_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                      else
                        FillUserTH2D("deltaREle2JetVsPt_failHEEP_diffLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                    }
                    else {
                      FillUserTH2D("deltaREle2JetVsPt_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                      if(ele2PassHeep)
                        FillUserTH2D("deltaREle2JetVsPt_passHEEP_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                      else
                        FillUserTH2D("deltaREle2JetVsPt_failHEEP_noLQmother_excludeClosestJet", ele2.SCEt(), t_ele2.DeltaR ( t_jet6 ));
                    }
                  }
                  if(ele2PassEt) {
                    FillUserTH1D("deltaREle2Jets_highET", t_ele2.DeltaR ( t_jet6 ));
                    if(closestJetIndex != jetRawIndicesMap[5]) {
                      FillUserTH1D("deltaREle2Jets_highET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet6 ));
                      if(ele2PassHeep)
                        FillUserTH1D("deltaREle2PassHEEPJets_highET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet6 ));
                      else if(!ele2PassHeepCaloIso) {
                        FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet", t_ele2.DeltaR ( t_jet6 ));
                        if(jetLQMotherIdx == ele2LQMotherIdx)
                          FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_sameLQMother", t_ele2.DeltaR ( t_jet6 ));
                        else if(jetLQMotherIdx > -1)
                          FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_diffLQMother", t_ele2.DeltaR ( t_jet6 ));
                        else
                          FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_highETExcludeClosestJet_noLQMother", t_ele2.DeltaR ( t_jet6 ));
                      }
                    }
                    if(isElectronBarrel(ele2.SCEta()))
                      FillUserTH1D("deltaREle2Jets_highET_barrel", t_ele2.DeltaR ( t_jet6 ));
                    else if(isElectronEndcap(ele2.SCEta()))
                      FillUserTH1D("deltaREle2Jets_highET_endcap", t_ele2.DeltaR ( t_jet6 ));
                  }
                  else {
                    FillUserTH1D("deltaREle2Jets_lowET", t_ele2.DeltaR ( t_jet6 ));
                    if(closestJetIndex != jetRawIndicesMap[5]) {
                      FillUserTH1D("deltaREle2Jets_lowET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet6 ));
                      if(ele2PassHeep)
                        FillUserTH1D("deltaREle2PassHEEPJets_lowET_ExcludeClosestJet", t_ele2.DeltaR ( t_jet6 ));
                      else if(!ele2PassHeepCaloIso) {
                        FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet", t_ele2.DeltaR ( t_jet6 ));
                        if(jetLQMotherIdx == ele2LQMotherIdx)
                          FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_sameLQMother", t_ele2.DeltaR ( t_jet6 ));
                        else if(jetLQMotherIdx > -1)
                          FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_diffLQMother", t_ele2.DeltaR ( t_jet6 ));
                        else
                          FillUserTH1D("deltaREle2FailHEEPCaloIsoJets_lowETExcludeClosestJet_noLQMother", t_ele2.DeltaR ( t_jet6 ));
                      }
                      if(!ele2PassHeep) {
                        if(t_ele2.DeltaR ( t_jet6 ) < 0.4 && ele2.SCEt() > 35 && ele2PassLoose)
                          STDOUT("--> found event with lowET ele2 failing heep; jet6 is the second-closest jet with deltaR = " << t_ele2.DeltaR ( t_jet6 ));
                      }
                    }
                    if(isElectronBarrel(ele2.SCEta()))
                      FillUserTH1D("deltaREle2Jets_lowET_barrel", t_ele2.DeltaR ( t_jet6 ));
                    else if(isElectronEndcap(ele2.SCEta()))
                      FillUserTH1D("deltaREle2Jets_lowET_endcap", t_ele2.DeltaR ( t_jet6 ));
                  }
                  if(closestJetIndex != jetRawIndicesMap[5]) {
                    FillUserTH1D("deltaREle2JetsExcludeClosestJet", t_ele2.DeltaR ( t_jet6 ));
                    if(ele2PassHeep)
                      FillUserTH1D("deltaREle2PassHEEPJetsExcludeClosestJet", t_ele2.DeltaR ( t_jet6 ));
                    else if(!ele2PassHeepCaloIso)
                      FillUserTH1D("deltaREle2FailHEEPCaloIsoJetsExcludeClosestJet", t_ele2.DeltaR ( t_jet6 ));
                  }
                  FillUserTH1D("deltaREle2Jets", t_ele2.DeltaR ( t_jet6 ));
                }
              }
            }
          }
        }
      }
    }
    //-----------------------------------------------------------------
    // QCD triggers 2016    2017        2018
    // - Photon22         
    //                      Photon25
    // - Photon30           
    //                      Photon33    Photon33
    // - Photon36           
    // - Photon50           Photon50    Photon50
    // - Photon75           Photon75    Photon75
    // - Photon90           Photon90    Photon90
    // - Photon120          Photon120   Photon120
    //                      Photon150   Photon150
    // - Photon175          Photon175   Photon175
    // - Photon200          Photon200   Photon200
    //-----------------------------------------------------------------

    if ( reducedSkimType == 0 ) { 
      // NB: for data, all prescales applied as averages in downstream analysis code
      if(triggerExists("HLT_Photon22"))
        fillTriggerVariable ( "HLT_Photon22" , "H_Photon22"  );
      else
        fillVariableWithValue( "H_Photon22", -1); 
      if(triggerExists("HLT_Photon25"))
        fillTriggerVariable ( "HLT_Photon25" , "H_Photon25"  );
      else
        fillVariableWithValue( "H_Photon25", -1); 
      if(triggerExists("HLT_Photon30"))
        fillTriggerVariable ( "HLT_Photon30" , "H_Photon30"  );
      else
        fillVariableWithValue( "H_Photon30", -1); 
      if(triggerExists("HLT_Photon33"))
        fillTriggerVariable ( "HLT_Photon33" , "H_Photon33"  );
      else
        fillVariableWithValue( "H_Photon33", -1); 
      if(triggerExists("HLT_Photon36"))
        fillTriggerVariable ( "HLT_Photon36" , "H_Photon36"  );
      else
        fillVariableWithValue( "H_Photon36", -1); 
      fillTriggerVariable ( "HLT_Photon50"   , "H_Photon50"  );
      fillTriggerVariable ( "HLT_Photon75"   , "H_Photon75"  );
      fillTriggerVariable ( "HLT_Photon90"   , "H_Photon90"  );
      fillTriggerVariable ( "HLT_Photon120"  , "H_Photon120" );
      if(triggerExists("HLT_Photon150"))
        fillTriggerVariable ( "HLT_Photon150", "H_Photon150" );
      else
        fillVariableWithValue( "H_Photon150", -1); 
      fillTriggerVariable ( "HLT_Photon175" , "H_Photon175" );
      if(triggerExists("HLT_Photon200"))
        fillTriggerVariable ( "HLT_Photon200" , "H_Photon200" );
      else
        fillVariableWithValue ( "H_Photon200" , -1 );

      bool pass_trigger = (
          getVariableValue("H_Photon22") > 0 || 
          getVariableValue("H_Photon25") > 0 || 
          getVariableValue("H_Photon30") > 0 || 
          getVariableValue("H_Photon33") > 0 || 
          getVariableValue("H_Photon36") > 0 || 
          getVariableValue("H_Photon50") > 0 || 
          getVariableValue("H_Photon75") > 0 || 
          getVariableValue("H_Photon90") > 0 || 
          getVariableValue("H_Photon120"     ) > 0 || 
          getVariableValue("H_Photon150"     ) > 0 || 
          getVariableValue("H_Photon175"     ) > 0 || 
          getVariableValue("H_Photon200"     ) > 0 );

      fillVariableWithValue ("PassTrigger", pass_trigger ? 1 : 0 );

    }

    else if ( reducedSkimType == 1 || reducedSkimType == 2 || reducedSkimType == 3 || reducedSkimType == 4 ) { 

      // search for HLT path by prefix
      // in 2015 data, trigger is different
      // this exists in special from-RAW MC and data only
      //fillTriggerVariable( "HLT_Ele27_WPLoose_Gsf_v" , "H_Ele27_WPLoose" );
      //if      ( isData )
      //{
      //  if(run >= 254227 && run <= 254914) // in Run2015C 25 ns, there is no un-eta-restricted WPLoose path
      //    fillTriggerVariable( "HLT_Ele27_eta2p1_WPLoose_Gsf_v" , "H_Ele27_WPLoose_eta2p1" );
      //  else
      //  {
      //    fillTriggerVariable( "HLT_Ele27_WPLoose_Gsf_v" , "H_Ele27_WPLoose" );
      //    fillTriggerVariable( "HLT_Ele27_eta2p1_WPLoose_Gsf_v" , "H_Ele27_WPLoose_eta2p1" );
      //    fillTriggerVariable( "HLT_Ele27_WPTight_Gsf_v" , "H_Ele27_WPTight" );
      //  }
      //}

      // just search by prefix
      //if(triggerExists("HLT_Ele27_WPLoose_Gsf"))
      //  fillTriggerVariable( "HLT_Ele27_WPLoose_Gsf" , "H_Ele27_WPLoose" );
      if(triggerExists("HLT_Ele27_WPTight_Gsf"))
        fillTriggerVariable( "HLT_Ele27_WPTight_Gsf" , "H_Ele27_WPTight" );
      if(triggerExists("HLT_Ele32_WPTight_Gsf"))
        fillTriggerVariable( "HLT_Ele32_WPTight_Gsf" , "H_Ele32_WPTight" );
      if(triggerExists("HLT_Ele35_WPTight_Gsf"))
        fillTriggerVariable( "HLT_Ele35_WPTight_Gsf" , "H_Ele35_WPTight" );
      // check that we have at least one WPTight trigger
      if(!triggerExists("HLT_Ele27_WPTight_Gsf") && !triggerExists("HLT_Ele32_WPTight_Gsf") && !triggerExists("HLT_Ele35_WPTight_Gsf"))
        exit(-5);
      // Ele115 is absent from first 5/fb of 2017
      if(triggerExists("HLT_Ele115_CaloIdVT_GsfTrkIdT"))
        fillTriggerVariable( "HLT_Ele115_CaloIdVT_GsfTrkIdT" , "H_Ele115_CIdVT_GsfIdT");
      if(triggerExists("HLT_Photon175"))
        fillTriggerVariable( "HLT_Photon175" , "H_Photon175" );
      if(triggerExists("HLT_Photon200"))
        fillTriggerVariable( "HLT_Photon200" , "H_Photon200" );
      // check that we have at least one photon trigger
      if(!triggerExists("HLT_Photon175") && !triggerExists("HLT_Photon200"))
        exit(-5);
      // other triggers
      //fillTriggerVariable( "HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50", "H_Ele45_PFJet200_PFJet50");
      if(triggerExists("HLT_Ele105_CaloIdVT_GsfTrkIdT"))
        fillTriggerVariable( "HLT_Ele105_CaloIdVT_GsfTrkIdT" , "H_Ele105_CIdVT_GsfIdT");
      if(triggerExists("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL"))
        fillTriggerVariable( "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL", "H_DoubleEle33_CIdL_GsfIdVL" ); 
      //fillTriggerVariable( "HLT_Mu45_eta2p1"  , "H_Mu45_eta2p1" );
    }

    //-----------------------------------------------------------------
    // Evaluate the cuts
    //-----------------------------------------------------------------    

    evaluateCuts();


    //-----------------------------------------------------------------    
    // Fill the trees
    //-----------------------------------------------------------------    
    //// QCD fake rate calculation skim
    //if ( reducedSkimType == 0 ) { 
    //  if(passedCut("PassTrigger"      ) &&  
    //      passedCut("nVLooseEle_ptCut"  ) && 
    //      //passedCut("LooseEle1_Pt"     ) ){
    //    passedCut("LooseEle1_Pt"     ) ){
    //      //c_ele_final->examine<Electron>("final electrons");
    //      //c_ele_final_ptCut->examine<Electron>("final electrons with Pt cut");
    //      fillSkimTree();
    //      fillReducedSkimTree();
    //    }
    //}

    //// enujj analysis skim
    //else if ( reducedSkimType == 1 ) { 
    //  if( passedCut("nEle_ptCut"       ) && 
    //      //passedCut("Ele1_Pt"          ) && 
    //      passedCut("Ele1_Pt"          ) && 
    //      passedCut("PFMET_Type1XY_Pt") && 
    //      passedCut("Jet1_Pt"          ) && 
    //      passedCut("Jet2_Pt"          ) && 
    //      passedCut("sT_enujj"         ) && 
    //      passedCut("MT_Ele1MET"       )) {
    //    fillSkimTree();
    //    fillReducedSkimTree();
    //  }
    //}

    //// eejj analysis skim
    //else if ( reducedSkimType == 2 ) { 
    //  if( passedCut("nEle_ptCut"       ) && 
    //      //passedCut("Ele1_Pt"          ) &&
    //      //passedCut("Ele2_Pt"          ) && 
    //      passedCut("Ele1_Pt"          ) && 
    //      passedCut("Ele2_Pt"          ) && 
    //      passedCut("Jet1_Pt"          ) &&
    //      passedCut("Jet2_Pt"          ) &&
    //      passedCut("sT_eejj"          ) && 
    //      passedCut("M_e1e2"           )) {
    //    fillSkimTree();
    //    fillReducedSkimTree();
    //  }
    //}

    //// Single muon skim
    //else if ( reducedSkimType == 4 ) { 
    //  if( passedCut("nMuon_ptCut"       ) && 
    //      passedCut("Muon1_Pt"          ) ){
    //    fillSkimTree();
    //    fillReducedSkimTree();
    //  }
    //}

  } // event loop
  std::cout << "analysisClass::Loop(): ends " << std::endl;
}
