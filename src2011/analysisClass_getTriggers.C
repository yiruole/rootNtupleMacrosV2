#define analysisClass_cxx
#include "analysisClass.h"
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

analysisClass::analysisClass(string * inputList, string * cutFile, string * treeName, string * outputFileName, string * cutEfficFile)
  :baseClass(inputList, cutFile, treeName, outputFileName, cutEfficFile){}

analysisClass::~analysisClass(){}

void analysisClass::Loop()
{
   std::cout << "analysisClass::Loop() begins" <<std::endl;   

   //--------------------------------------------------------------------------
   // Decide which plots to save (default is to save everything)
   //--------------------------------------------------------------------------
   
   fillSkim                         ( !true  ) ;
   fillAllPreviousCuts              ( !true  ) ;
   fillAllOtherCuts                 ( !true  ) ;
   fillAllSameLevelAndLowerLevelCuts( !true  ) ;
   fillAllCuts                      ( !true  ) ;

   //--------------------------------------------------------------------------
   // Loop over the chain
   //--------------------------------------------------------------------------

   if (fChain == 0) return;
   
   Long64_t nentries = fChain -> GetEntries();
   std::cout << "analysisClass::Loop(): nentries = " << nentries << std::endl;   

   std::map<std::string, std::pair <int,int>  > trigger_name_map;
   std::map<std::string, std::pair <int,int>  >::iterator it;
   std::map<std::string, std::pair <int,int>  >::iterator it_end = trigger_name_map.end();

   std::string old_key("");

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     if(jentry < 10 || jentry%1000 == 0) std::cout << "analysisClass::Loop(): jentry = " << jentry << "/" << nentries << std::endl;

     bool is_the_old_key = ( HLTKey -> compare ( old_key ) == 0 );

     if ( is_the_old_key ) continue;

     std::vector<std::string>::iterator inside_name      = HLTInsideDatasetTriggerNames -> begin();
     std::vector<std::string>::iterator inside_name_end  = HLTInsideDatasetTriggerNames -> end();

     std::vector<std::string>::iterator outside_name     = HLTOutsideDatasetTriggerNames -> begin();
     std::vector<std::string>::iterator outside_name_end = HLTOutsideDatasetTriggerNames -> end();

     for (; inside_name != inside_name_end; ++inside_name ) {
       it = trigger_name_map.find ( *inside_name );
       if ( it == trigger_name_map.end() ) trigger_name_map[*inside_name] = std::pair<int,int> ( run, run ) ;
       else { 
	 if ( trigger_name_map[*inside_name].first  > run ) trigger_name_map[*inside_name].first  = run;
	 if ( trigger_name_map[*inside_name].second < run ) trigger_name_map[*inside_name].second = run;
       }
     }

     for (; outside_name != outside_name_end; ++outside_name ) {
       it = trigger_name_map.find ( *outside_name );
       if ( it == trigger_name_map.end() ) trigger_name_map[*outside_name] = std::pair<int,int> ( run, run ) ;
       else { 
	 if ( trigger_name_map[*outside_name].first  > run ) trigger_name_map[*outside_name].first  = run;
	 if ( trigger_name_map[*outside_name].second < run ) trigger_name_map[*outside_name].second = run;
       }
     }
     
     old_key = *HLTKey;

   } // End loop over events

   it = trigger_name_map.begin();
   it_end = trigger_name_map.end();

   std::cout << "I found these triggers" << std::endl;

   for (; it != it_end; ++it){
     std::cout << "  FOUND_A_TRIGGER " << it -> second.first << ", " << it -> second.second << ": \t" << it -> first << std::endl;
   }

   std::cout << "analysisClass::Loop() ends" <<std::endl;   
}

