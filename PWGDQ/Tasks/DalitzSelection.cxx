// This task selects electrons from Dalitz decay and writes bits on every track for selection/rejection in later analysis
// Input are: names of the electron cuts, names of the pair cuts
// For every cut, form every possible electron pair and tag both electrons if the pair passes the pair cut
// produces pair and single track qa plots for every (electron cut, pair cut) combination


#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"



using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::soa;
using std::array;


namespace o2::aod
{
namespace dqdalitzflags
{
DECLARE_SOA_COLUMN(IsEventSelectedDal, isEventSelectedDal, int);
DECLARE_SOA_COLUMN(IsBarrelSelectedDal, isBarrelSelectedDal, int);
DECLARE_SOA_COLUMN(IsReducedEventSelectedDal, isReducedEventSelectedDal, int);
DECLARE_SOA_COLUMN(IsReducedBarrelSelectedDal, isReducedBarrelSelectedDal, int);
}
DECLARE_SOA_TABLE(EventSelectedDal, "AOD", "DQEVENTSELDAL", dqdalitzflags::IsEventSelectedDal);
DECLARE_SOA_TABLE(BarrelSelectedDal, "AOD", "DQTRACKSELDAL", dqdalitzflags::IsBarrelSelectedDal);
DECLARE_SOA_TABLE(ReducedEventSelectedDal, "AOD", "DQREDEVSELDAL", dqdalitzflags::IsReducedEventSelectedDal);
DECLARE_SOA_TABLE(ReducedBarrelSelectedDal, "AOD", "DQREDTRSELDAL", dqdalitzflags::IsReducedBarrelSelectedDal);
}



using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyReducedEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;

using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA,
                                aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullMu,
                                aod::pidTPCFullKa, aod::pidTPCFullPr,
                                aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullMu,
                                aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;            
using MyBarrelSelectedTracks = soa::Join<MyBarrelTracks, o2::aod::BarrelSelectedDal>;
using BarrelReducedTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID>;
using MyBarrelSelectedReducedTracks =  soa::Join<BarrelReducedTracks, o2::aod::ReducedBarrelSelectedDal>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::Collision;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkReducedEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkReducedTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;


// Variables and functions common for full and reduced tracks

std::map<int, int> eventmap; //only select events with at least 2 dalitz candidates
std::map<int, int> trackmap;
std::map<int, uint8_t> dalitzmap; 
int nCuts;
bool fQA;

HistogramManager* fHistMan;
TList* statsList;

//Cuts for tracks selection
std::vector<AnalysisCompositeCut> fTrackSelCuts;
AnalysisCompositeCut* fEventCut;

// Cuts for pairing
std::vector<AnalysisCompositeCut> fTrackCuts; //! Barrel track cuts
std::vector<AnalysisCompositeCut> fPairCuts; //! Dalitz pair cuts

template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
void runTrackSelection(TEvent const& collisions, TTracks const& tracksBarrel);
template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
void runDalitzPairing(TEvent const& event, TTracks const& tracksP, TTracks const& tracksN);
void InitTrackSelection(o2::framework::InitContext& context, std::string fConfigTrackCuts, std::string fConfigEventCuts);
void InitPairing(o2::framework::InitContext& context, std::string fConfigTrackCuts, std::string fConfigPairCuts, bool cfgQA);
void DefineHistograms(TString histClasses, o2::framework::InitContext& context);


///////////////////////////////////////
// Track selection for full tracks ////
///////////////////////////////////////

struct dalitzTrackSelection {
  Produces<aod::BarrelSelectedDal> trackSel;
  Produces<aod::EventSelectedDal> eventSel;
  
  Configurable<std::string> fConfigTrackCuts{"cfgDalitzTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandardNoINT7", "Event selection"};
  Configurable<float> fConfigBarrelTrackPtLow{"cfgBarrelLowPt", 0.1f, "Low pt cut for tracks in the barrel"};
  
  Partition<MyBarrelTracks> barrelTracksFilter = o2::aod::track::pt >= fConfigBarrelTrackPtLow  && nabs(o2::aod::track::eta) <= 0.9f && nabs(o2::aod::pidtpc::tpcNSigmaEl) <= 3.0f && o2::aod::track::tpcChi2NCl < 4.0f && o2::aod::track::itsChi2NCl < 36.0f;
  
  
  void init(o2::framework::InitContext& context)
  {
    InitTrackSelection(context, fConfigTrackCuts, fConfigEventCuts);
  }   
  
  void processFullTracks(MyEvents const& collisions, MyBarrelTracks const& tracks)
  {
    eventmap.clear();
    trackmap.clear();
    
    runTrackSelection<gkEventFillMap,gkTrackFillMap>(collisions, barrelTracksFilter);
    
    for (auto& collision : collisions) {
      eventSel(eventmap[collision.globalIndex()]);
    }
    
    trackSel.reserve(tracks.size());
    for (auto& track : tracks) {
      trackSel(trackmap[track.globalIndex()]);
    }

  }
  
  void processDummy(MyEvents&)
  {
  }

  PROCESS_SWITCH(dalitzTrackSelection, processFullTracks, "Run dalitz track selection on full tracks", false);
  PROCESS_SWITCH(dalitzTrackSelection, processDummy, "process dummy", false);
  
};



/////////////////////////////////////////
// Track selection for reduced track  ///
/////////////////////////////////////////

struct dalitzReducedTrackSelection {
  Produces<aod::ReducedBarrelSelectedDal> trackSel;
  Produces<aod::ReducedEventSelectedDal> eventSel;
  
  Configurable<std::string> fConfigTrackCuts{"cfgDalitzTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandardNoINT7", "Event selection"};
  Configurable<float> fConfigBarrelTrackPtLow{"cfgBarrelLowPt", 0.1f, "Low pt cut for tracks in the barrel"};
  
  Partition<BarrelReducedTracks> barrelReducedTracksFilter = o2::aod::reducedtrack::pt >= fConfigBarrelTrackPtLow  && nabs(o2::aod::reducedtrack::eta) <= 0.9f && nabs(o2::aod::pidtpc::tpcNSigmaEl) <= 3.0f && o2::aod::track::tpcChi2NCl < 4.0f && o2::aod::track::itsChi2NCl < 36.0f;
    
  void init(o2::framework::InitContext& context)
  {
    InitTrackSelection(context, fConfigTrackCuts, fConfigEventCuts);
  }   
  
  void processReducedTracks(MyReducedEvents const& collisions, BarrelReducedTracks const& tracks)
  {
    eventmap.clear();
    trackmap.clear();
    
    runTrackSelection<gkReducedEventFillMap,gkReducedTrackFillMap>(collisions, barrelReducedTracksFilter);
    
    for (auto& collision : collisions) {
      eventSel(eventmap[collision.globalIndex()]);
    }
    for (auto& track : tracks) {
      trackSel(trackmap[track.globalIndex()]);
    } 
  }
  
  void processDummy(MyEvents&)
  {
  }

  PROCESS_SWITCH(dalitzReducedTrackSelection, processReducedTracks, "Run dalitz track selection on reduced tracks", false);
  PROCESS_SWITCH(dalitzReducedTrackSelection, processDummy, "process dummy", false);
  
};



////////////////////////////////////
///// Pairing for full tracks //////
////////////////////////////////////

struct dalitzPairing {  
  Produces<o2::aod::DalitzBits> dalitzbits;
   
  Filter eventSel = o2::aod::dqdalitzflags::isEventSelectedDal == 1;
  Partition<MyBarrelSelectedTracks> barrelTrackSelP = o2::aod::dqdalitzflags::isBarrelSelectedDal > 0 && o2::aod::track::signed1Pt > 0.f;  
  Partition<MyBarrelSelectedTracks> barrelTrackSelN = o2::aod::dqdalitzflags::isBarrelSelectedDal > 0 && o2::aod::track::signed1Pt < 0.f;  

  //Configurables
  Configurable<std::string> fConfigTrackCuts{"cfgDalitzTrackCuts", "", "Dalitz track selection cuts, separated by a comma"};
  Configurable<std::string> fConfigPairCuts{"cfgDalitzPairCuts", "", "Dalitz pair selection cuts"};
  Configurable<bool> cfgQA{"cfgQA", true, "QA histograms"};

  OutputObj<THashList> fOutputList{"output"}; //! the histogram manager output list
  OutputObj<TList> fStatsList{"Statistics"};  //! skimming statistics

  void init(o2::framework::InitContext& context)
  {
    InitPairing(context, fConfigTrackCuts, fConfigPairCuts, cfgQA); 
    if (context.mOptions.get<bool>("processDummy")) return;
    fStatsList.setObject(statsList); 
    fOutputList.setObject(fHistMan->GetMainHistogramList());  
  }

  void processFullTracks(soa::Filtered<soa::Join<MyEvents, o2::aod::EventSelectedDal>> const& collisions, MyBarrelSelectedTracks const& tracks)
  {
    dalitzmap.clear();

    for (auto& collision : collisions) {
      auto groupedTracksP = barrelTrackSelP->sliceByCached(aod::track::collisionId, collision.globalIndex());
      auto groupedTracksN = barrelTrackSelN->sliceByCached(aod::track::collisionId, collision.globalIndex());
      runDalitzPairing<gkEventFillMap,gkTrackFillMap>(collision, groupedTracksP, groupedTracksN);
    }
    
    //Fill dalitz bits
    for (auto& track : tracks) {
      dalitzbits(dalitzmap[track.globalIndex()]);
    } 
  }
 
  void processDummy(MyEvents&)
  {
  }

  PROCESS_SWITCH(dalitzPairing, processFullTracks, "Run Dalitz selection on AO2D tables", false);
  PROCESS_SWITCH(dalitzPairing, processDummy, "Do nothing", false);


};


////////////////////////////////////////
////// Pairing for reduced tracks //////
////////////////////////////////////////

struct dalitzReducedPairing {  
  Produces<o2::aod::DalitzBitsReduced> dalitzbits;
   
  Filter eventSel = o2::aod::dqdalitzflags::isReducedEventSelectedDal == 1;
  Partition<MyBarrelSelectedReducedTracks> barrelReducedTrackSelP = o2::aod::dqdalitzflags::isReducedBarrelSelectedDal > 0 && o2::aod::reducedtrack::sign > 0.f;
  Partition<MyBarrelSelectedReducedTracks> barrelReducedTrackSelN = o2::aod::dqdalitzflags::isReducedBarrelSelectedDal > 0 && o2::aod::reducedtrack::sign < 0.f; 

  //Configurables
  Configurable<std::string> fConfigTrackCuts{"cfgDalitzTrackCuts", "", "Dalitz track selection cuts, separated by a comma"};
  Configurable<std::string> fConfigPairCuts{"cfgDalitzPairCuts", "", "Dalitz pair selection cuts"};
  Configurable<bool> cfgQA{"cfgQA", true, "QA histograms"};
  
  OutputObj<THashList> fOutputList{"output"}; //! the histogram manager output list
  OutputObj<TList> fStatsList{"Statistics"};  //! skimming statistics


  void init(o2::framework::InitContext& context)
  {
    InitPairing(context, fConfigTrackCuts, fConfigPairCuts, cfgQA); 
    if (context.mOptions.get<bool>("processDummy")) return; 
    fStatsList.setObject(statsList);  
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  void processReducedTracks(soa::Filtered<soa::Join<MyReducedEvents, o2::aod::ReducedEventSelectedDal>> const& events, MyBarrelSelectedReducedTracks const& tracks)
  {
    dalitzmap.clear();
    
    for (auto& event : events) {
      auto groupedTracksP = barrelReducedTrackSelP->sliceByCached(aod::reducedtrack::reducedeventId, event.globalIndex());
      auto groupedTracksN = barrelReducedTrackSelN->sliceByCached(aod::reducedtrack::reducedeventId, event.globalIndex());
      runDalitzPairing<gkReducedEventFillMap,gkReducedTrackFillMap>(event, groupedTracksP, groupedTracksN);
    }
    
    for (auto& track : tracks) {
      dalitzbits(dalitzmap[track.globalIndex()]);
    }
  }
  
  void processDummy(MyEvents&)
  {
  }

  PROCESS_SWITCH(dalitzReducedPairing, processReducedTracks, "Run Dalitz selection on skimmed tables", false);
  PROCESS_SWITCH(dalitzReducedPairing, processDummy, "Do nothing", false);


};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dalitzTrackSelection>(cfgc),
    adaptAnalysisTask<dalitzReducedTrackSelection>(cfgc),
    adaptAnalysisTask<dalitzPairing>(cfgc),
    adaptAnalysisTask<dalitzReducedPairing>(cfgc)    
  };
}


/////////////////////////////////////////////////////
///// Templated function for track selection ////////
/////////////////////////////////////////////////////

template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
void runTrackSelection(TEvent const& collisions, TTracks const& tracksBarrel)
{
  uint8_t filterMap = uint8_t(0);
  int CollisionId = -1;
  bool isEventSelected = true;
  uint8_t eventOneTrackSelected = uint8_t(0);

  VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);

  for (auto& track : tracksBarrel) {
    filterMap = uint8_t(0);
    int trackCollisionId;
    if constexpr (static_cast<bool>(TEventFillMap & VarManager::ObjTypes::Collision)) {
      trackCollisionId = track.collisionId();
    }  
    if constexpr (static_cast<bool>(TEventFillMap & VarManager::ObjTypes::ReducedEvent)) {
      trackCollisionId = track.reducedeventId();
    }
      
    if (trackCollisionId >= 0) {
      // fill event information which might be needed in histograms or cuts that combine track and event properties
      if (trackCollisionId != CollisionId) { // check if the track belongs to a different event than the previous one
        CollisionId = trackCollisionId;
        if constexpr (static_cast<bool>(TEventFillMap & VarManager::ObjTypes::Collision)) {
          auto collision = track.template collision_as<TEvent>();
          VarManager::FillEvent<TEventFillMap>(collision);
        }
        if constexpr (static_cast<bool>(TEventFillMap & VarManager::ObjTypes::ReducedEvent)) {
          auto collision = track.template reducedevent_as<TEvent>();
          VarManager::FillEvent<TEventFillMap>(collision);
        }
          
        isEventSelected = fEventCut->IsSelected(VarManager::fgValues);
        eventOneTrackSelected = uint8_t(0);
      }
      if (isEventSelected) {
        VarManager::FillTrack<TTrackFillMap>(track);
        int i = 0;
        for (auto cut = fTrackSelCuts.begin(); cut != fTrackSelCuts.end(); ++cut, ++i) {
          if ((*cut).IsSelected(VarManager::fgValues)) {
            filterMap |= (uint8_t(1) << i);
          }
        }
        if (filterMap) {// Fill event selection
          if (eventOneTrackSelected & filterMap) {
            eventmap[CollisionId] = 1;
          }
          eventOneTrackSelected = eventOneTrackSelected | filterMap;            
        }
      }
      trackmap[track.globalIndex()] = static_cast<int>(filterMap);
    }
  } // end loop over tracks
}


/////////////////////////////////////////////////////
////// Templated function for track pairing /////////
/////////////////////////////////////////////////////

template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
void runDalitzPairing(TEvent const& event, TTracks const& tracksP, TTracks const& tracksN)
{
  const int TPairType = VarManager::kJpsiToEE;
  for (auto& [track1,track2] : o2::soa::combinations(CombinationsFullIndexPolicy(tracksP, tracksN))) {
      
    uint8_t twoTracksFilterMap;
    if constexpr (static_cast<bool>(TEventFillMap & VarManager::ObjTypes::Collision)) {
      twoTracksFilterMap = uint8_t(track1.isBarrelSelectedDal()) & uint8_t(track2.isBarrelSelectedDal());
    }
    if constexpr (static_cast<bool>(TEventFillMap & VarManager::ObjTypes::ReducedEvent)) {
      twoTracksFilterMap = uint8_t(track1.isReducedBarrelSelectedDal()) & uint8_t(track2.isReducedBarrelSelectedDal());
    }
    if (!twoTracksFilterMap) continue;
      
    // pairing
    VarManager::FillPair<TPairType, TTrackFillMap>(track1, track2);  
    uint8_t track1Untagged = uint8_t(0);
    uint8_t track2Untagged = uint8_t(0);
      
    // Fill pair selection map and fill pair histogram
    for (int icut = 0; icut < nCuts; icut++) {
      if (!(twoTracksFilterMap & (uint8_t(1) << icut))) continue;
      AnalysisCompositeCut pairCut = fPairCuts.at(icut);
      if (pairCut.IsSelected(VarManager::fgValues)) {
        if (fQA) {         
          AnalysisCompositeCut trackCut = fTrackCuts.at(icut);
          fHistMan->FillHistClass(Form("Pair_%s_%s", trackCut.GetName(), pairCut.GetName()), VarManager::fgValues);  
        }     
          
        // Check if tracks were already tagged       
        bool b1 = dalitzmap[track1.globalIndex()] & (uint8_t(1) << icut);
        bool b2 = dalitzmap[track2.globalIndex()] & (uint8_t(1) << icut);
        if(!b1) {
          track1Untagged |= (uint8_t(1) << icut);
          ((TH1I*) statsList->At(0))->Fill(icut);
        }
        if(!b2) {
          track2Untagged |= (uint8_t(1) << icut);
          ((TH1I*) statsList->At(0))->Fill(icut);
        }          
      }
    }
      
    // Tag tracks which are not already tagged           
    dalitzmap[track1.globalIndex()] |= track1Untagged;
    dalitzmap[track2.globalIndex()] |= track2Untagged;
      
    // Fill track histograms if not already tagged
    if (fQA) {
      VarManager::FillTrack<TTrackFillMap>(track1);
      for (int icut = 0; icut < nCuts; icut++) {
        if (track1Untagged & (uint8_t(1) << icut)) {
          AnalysisCompositeCut trackCut = fTrackCuts.at(icut);
          AnalysisCompositeCut pairCut = fPairCuts.at(icut);
          fHistMan->FillHistClass(Form("TrackBarrel_%s_%s", trackCut.GetName(), pairCut.GetName()), VarManager::fgValues);
        }
      }
      VarManager::FillTrack<TTrackFillMap>(track2);
      for (int icut = 0; icut < nCuts; icut++) {
        if (track2Untagged & (uint8_t(1) << icut)) {
          AnalysisCompositeCut trackCut = fTrackCuts.at(icut);
          AnalysisCompositeCut pairCut = fPairCuts.at(icut);
          fHistMan->FillHistClass(Form("TrackBarrel_%s_%s", trackCut.GetName(), pairCut.GetName()), VarManager::fgValues);
        }
      }
    } 
      
  } //end of tracksP,N loop
}



/////////////////////////////////////////////////
////// Initialization for track selection ///////
/////////////////////////////////////////////////

void InitTrackSelection(o2::framework::InitContext& context, std::string fConfigTrackCuts, std::string fConfigEventCuts) 
{

  if (context.mOptions.get<bool>("processDummy")) return;

  // Track cuts
  TString cutNamesStr = (TString) fConfigTrackCuts.data();
  if (!cutNamesStr.IsNull()) {
    std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
    for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
      fTrackSelCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
    }
  }

  // Event cuts
  fEventCut = new AnalysisCompositeCut(true);
  TString eventCutStr = (TString) fConfigEventCuts.data();
  fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data())); 
     
  VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

}


/////////////////////////////////////////////////
/////// Initialization for track pairing ////////
/////////////////////////////////////////////////

void InitPairing(o2::framework::InitContext& context, std::string fConfigTrackCuts, std::string fConfigPairCuts, bool cfgQA) 
{
  
  if (context.mOptions.get<bool>("processDummy")) return;
  
  // Barrel track cuts
  TString cutNamesStr = (TString) fConfigTrackCuts.data();
  if (!cutNamesStr.IsNull()) {
    std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
    for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
      fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
    }
  }
    
  // Pair cuts
  TString cutNamesPairStr = (TString) fConfigPairCuts.data();
  if (!cutNamesPairStr.IsNull()) {
    std::unique_ptr<TObjArray> objArray(cutNamesPairStr.Tokenize(","));
    for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
      fPairCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
    }
  }
       
  VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill

  VarManager::SetDefaultVarNames();
  fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
  fHistMan->SetUseDefaultVariableNames(kTRUE);
  fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

  // Create the histogram class names to be added to the histogram manager
  TString histClasses = "";

  if (cfgQA) { 
    for (int icut = 0; icut < nCuts; icut++) {
      AnalysisCompositeCut trackCut = fTrackCuts.at(icut);
      AnalysisCompositeCut pairCut = fPairCuts.at(icut);
      histClasses += Form("TrackBarrel_%s_%s;", trackCut.GetName(), pairCut.GetName());
      histClasses += Form("Pair_%s_%s;", trackCut.GetName(), pairCut.GetName());
    }
  }
    
  DefineHistograms(histClasses, context);                   // define all histograms
  VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill

  if(fTrackCuts.size() != fPairCuts.size() ) {
    std::cout<<"WARNING: YOU NEED THE SAME NUMBER OF TRACK AND PAIR CUTS"<<std::endl;
  }   
  
  nCuts = std::min(fTrackCuts.size(), fPairCuts.size());
  fQA = cfgQA;
}




////////////////////////////////////////////////
//////// Define histograms (pairing) ///////////
////////////////////////////////////////////////

void DefineHistograms(TString histClasses, o2::framework::InitContext& context)
{
  std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
  for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
    TString classStr = objArray->At(iclass)->GetName();
    fHistMan->AddHistClass(classStr.Data());

    if (classStr.Contains("Event")) {
      dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "event", "");
    }

    if (classStr.Contains("Track")) {
      dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track","kine,tpcpid,tofpid,dca");
    }
      
    if (classStr.Contains("Pair")) {
      dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "pair_barrel","");
    }
  }

  statsList = new TList();
  statsList->SetOwner(kTRUE);

  // Dalitz selection statistics: one bin for each (track,pair) selection
  TH1I* histTracks = new TH1I("TrackStats", "Dalitz selection statistics", nCuts, -0.5, nCuts - 0.5 );
  for (int icut = 0; icut < nCuts; icut++) {
    AnalysisCompositeCut trackCut = fTrackCuts.at(icut);
    AnalysisCompositeCut pairCut = fPairCuts.at(icut);
    histTracks->GetXaxis()->SetBinLabel(icut, Form("%s_%s", trackCut.GetName(), pairCut.GetName()));
  }
  statsList->Add(histTracks);    
}

