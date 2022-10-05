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
}
DECLARE_SOA_TABLE(EventSelectedDal, "AOD", "DQEVENTSELDAL", dqdalitzflags::IsEventSelectedDal);
DECLARE_SOA_TABLE(BarrelSelectedDal, "AOD", "DQTRACKSELDAL", dqdalitzflags::IsBarrelSelectedDal);
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
using MyBarrelSelectedReducedTracks =  soa::Join<BarrelReducedTracks, o2::aod::BarrelSelectedDal>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::Collision;
constexpr static uint32_t gkTrackFillMap = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackPID;
constexpr static uint32_t gkReducedEventFillMap = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended;
constexpr static uint32_t gkReducedTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;




struct dalitzTrackSelection {
  Produces<aod::BarrelSelectedDal> trackSel;
  Produces<aod::EventSelectedDal> eventSel;
  
  
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandardNoINT7", "Event selection"};
  Configurable<float> fConfigBarrelTrackPtLow{"cfgBarrelLowPt", 0.1f, "Low pt cut for tracks in the barrel"};
  
  Partition<MyBarrelTracks> barrelTracksFilter = o2::aod::track::pt >= fConfigBarrelTrackPtLow  && nabs(o2::aod::track::eta) <= 0.9f && nabs(o2::aod::pidtpc::tpcNSigmaEl) <= 3.0f && o2::aod::track::tpcChi2NCl < 4.0f && o2::aod::track::itsChi2NCl < 36.0f;
  
  Partition<BarrelReducedTracks> barrelReducedTracksFilter = o2::aod::reducedtrack::pt >= fConfigBarrelTrackPtLow  && nabs(o2::aod::reducedtrack::eta) <= 0.9f && nabs(o2::aod::pidtpc::tpcNSigmaEl) <= 3.0f && o2::aod::track::tpcChi2NCl < 4.0f && o2::aod::track::itsChi2NCl < 36.0f;
  
  std::vector<AnalysisCompositeCut> fTrackCuts;
  AnalysisCompositeCut* fEventCut;              //! Event selection cut
  std::map<int, int> eventmap; //only select events with at least 2 dalitz candidates
  std::map<int, int> trackmap;
  
  void init(o2::framework::InitContext&)
  {
    // Track cuts
    TString cutNamesStr = fConfigTrackCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }

    // Event cuts
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data())); 
    
    
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }   

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runTrackSelection(TEvent const& collisions, TTracks const& tracksBarrel)
  {
    uint16_t filterMap = uint16_t(0);
    int CollisionId = -1;
    bool isEventSelected = true;
    uint16_t eventOneTrackSelected = uint16_t(0);

    VarManager::ResetValues(0, VarManager::kNBarrelTrackVariables);

    for (auto& track : tracksBarrel) {
      filterMap = uint16_t(0);
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
          eventOneTrackSelected = uint16_t(0);
        }
        if (isEventSelected) {
          VarManager::FillTrack<TTrackFillMap>(track);
          int i = 0;
          for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); ++cut, ++i) {
            if ((*cut).IsSelected(VarManager::fgValues)) {
              filterMap |= (uint16_t(1) << i);
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
  /*void processReducedTracks(MyReducedEvents const& collisions, BarrelReducedTracks const& tracks)
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

  }*/

  PROCESS_SWITCH(dalitzTrackSelection, processFullTracks, "Run dalitz track selection on full tracks", false);
  //PROCESS_SWITCH(dalitzTrackSelection, processReducedTracks, "Run dalitz track selection on reduced tracks", false);
  //PROCESS_SWITCH(dalitzTrackSelection, processDummy, "process dummy", false);
  
};





struct dalitzPairing {  
  Produces<o2::aod::DalitzBits> dalitzbits;
   
  Filter eventSel = o2::aod::dqdalitzflags::isEventSelectedDal == 1;
  Partition<MyBarrelSelectedTracks> barrelTrackSelP = o2::aod::dqdalitzflags::isBarrelSelectedDal > 0 && o2::aod::track::signed1Pt > 0.f;  
  Partition<MyBarrelSelectedReducedTracks> barrelReducedTrackSelP = o2::aod::dqdalitzflags::isBarrelSelectedDal > 0 && o2::aod::reducedtrack::sign > 0.f;
  Partition<MyBarrelSelectedTracks> barrelTrackSelN = o2::aod::dqdalitzflags::isBarrelSelectedDal > 0 && o2::aod::track::signed1Pt < 0.f;  
  Partition<MyBarrelSelectedReducedTracks> barrelReducedTrackSelN = o2::aod::dqdalitzflags::isBarrelSelectedDal > 0 && o2::aod::reducedtrack::sign < 0.f; 

  //Configurables
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Dalitz track selection cuts, separated by a comma"};
  Configurable<std::string> fConfigPairCuts{"cfgPairCuts", "", "Dalitz pair selection cuts"};
  Configurable<bool> fQA{"cfgQA", true, "QA histograms"};

  std::vector<AnalysisCompositeCut> fTrackCuts; //! Barrel track cuts
  std::vector<AnalysisCompositeCut> fPairCuts; //! Dalitz pair cuts

  OutputObj<THashList> fOutputList{"output"}; //! the histogram manager output list
  OutputObj<TList> fStatsList{"Statistics"};  //! skimming statistics
  HistogramManager* fHistMan;

  std::map<int, uint32_t> dalitzmap; 

  void init(o2::framework::InitContext& context)
  {
    DefineCuts();

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Create the histogram class names to be added to the histogram manager
    TString histClasses = "";

    if (fQA) { 
      for (auto& cutT : fTrackCuts) {
        for (auto& cutP : fPairCuts) {// single track and pair histograms
          histClasses += Form("TrackBarrel_%s_%s;", cutT.GetName(), cutP.GetName());
          histClasses += Form("Pair_%s_%s;", cutT.GetName(), cutP.GetName());
        }
      }
    }
    
    DefineHistograms(histClasses, context);                   // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
    
  }

  void DefineCuts()
  {
    // Barrel track cuts
    TString cutNamesStr = fConfigTrackCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    
    // Pair cuts
    TString cutNamesPairStr = fConfigPairCuts.value;
    if (!cutNamesPairStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesPairStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fPairCuts.push_back(*dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
       
    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

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

    fStatsList.setObject(new TList());
    fStatsList->SetOwner(kTRUE);

    // Dalitz selection statistics: one bin for each (track,pair) selection
    TH1I* histTracks = new TH1I("TrackStats", "Dalitz selection statistics", fTrackCuts.size()*fPairCuts.size(), -0.5, fTrackCuts.size()*fPairCuts.size() - 0.5 );
    int ib = 1;      
    for (auto cutP = fPairCuts.begin(); cutP != fPairCuts.end(); cutP++, ib++) {
      for (auto cutT = fTrackCuts.begin(); cutT != fTrackCuts.end(); cutT++, ib++) {
        histTracks->GetXaxis()->SetBinLabel(ib, Form("%s_%s", (*cutT).GetName(), (*cutP).GetName()));
      }
    }
    fStatsList->Add(histTracks);    
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, typename TEvent, typename TTracks>
  void runDalitzSelection(TEvent const& event, TTracks const& tracksP, TTracks const& tracksN)
  {
    const int TPairType = VarManager::kJpsiToEE;
    for (auto& [track1,track2] : o2::soa::combinations(CombinationsFullIndexPolicy(tracksP, tracksN))) {
      
      uint8_t twoTracksFilterMap = uint8_t(track1.isBarrelSelectedDal()) & uint8_t(track2.isBarrelSelectedDal());
      if (!twoTracksFilterMap) continue;
      
      // pairing
      VarManager::FillPair<TPairType, TTrackFillMap>(track1, track2);    

      // pair cuts
      int iT = 0;  
      for (auto cutT = fTrackCuts.begin(); cutT != fTrackCuts.end(); cutT++, iT++) {
        if (!(twoTracksFilterMap & (uint8_t(1) << iT))) continue;
        int iP = 0;
        for (auto cutP = fPairCuts.begin(); cutP != fPairCuts.end(); cutP++, iP++) {
          if ((*cutP).IsSelected(VarManager::fgValues)) {

            // tag the tracks and fill hists
            int bitNo = fTrackCuts.size()*iP+iT;
            if (fQA) { 
              fHistMan->FillHistClass(Form("Pair_%s_%s", (*cutT).GetName(), (*cutP).GetName()), VarManager::fgValues);  
            } 

            bool b1 = dalitzmap[track1.globalIndex()] & (uint32_t(1) << bitNo);
            if (!b1) { //avoid double counting in histograms
              dalitzmap[track1.globalIndex()] |= (uint32_t(1) << bitNo);
              if (fQA) {
                VarManager::FillTrack<TTrackFillMap>(track1);
                fHistMan->FillHistClass(Form("TrackBarrel_%s_%s", (*cutT).GetName(), (*cutP).GetName()), VarManager::fgValues);
              }  
              ((TH1I*) fStatsList->At(0))->Fill(bitNo);
            } //endif not already filled
              	
            bool b2 = dalitzmap[track2.globalIndex()] & (uint32_t(1) << bitNo);
            if (!b2) { //avoid double counting in histograms
              dalitzmap[track2.globalIndex()] |= (uint32_t(1) << bitNo);
              if (fQA) {
                VarManager::FillTrack<TTrackFillMap>(track2);
                fHistMan->FillHistClass(Form("TrackBarrel_%s_%s", (*cutT).GetName(), (*cutP).GetName()), VarManager::fgValues);
              }
              ((TH1I*) fStatsList->At(0))->Fill(bitNo);
            } //endif not already filled
         	
          } // end if selected     
      	} //end of pair cuts
      } //end of track cuts
    } //end of tracksP,N loop
  }

  void processFullTracks(soa::Filtered<soa::Join<MyEvents, o2::aod::EventSelectedDal>> const& collisions, MyBarrelSelectedTracks const& tracks)
  {
    if(fTrackCuts.size()==0 || fPairCuts.size()==0 ) {
      std::cout<<"WARNING: YOU NEED THE SAME NUMBER OF TRACK AND PAIR CUTS"<<std::endl;
    }   

    dalitzmap.clear();

    for (auto& collision : collisions) {
      auto groupedTracksP = barrelTrackSelP->sliceByCached(aod::track::collisionId, collision.globalIndex());
      auto groupedTracksN = barrelTrackSelN->sliceByCached(aod::track::collisionId, collision.globalIndex());
      runDalitzSelection<gkEventFillMap,gkTrackFillMap>(collision, groupedTracksP, groupedTracksN);
    }
    
    //Fill dalitz bits
    for (auto& track : tracks) {
      dalitzbits(dalitzmap[track.globalIndex()]);
    } 
  }
 
  
 /*void processReducedTracks(soa::Filtered<soa::Join<MyReducedEvents, o2::aod::EventSelectedDal>> const& events, MyBarrelSelectedReducedTracks const& tracks)
  {
    if(fTrackCuts.size() != fPairCuts.size() ) {
      std::cout<<"WARNING: YOU NEED THE SAME NUMBER OF TRACK AND PAIR CUTS"<<std::endl;
    }
    dalitzmap.clear();
    
    for (auto& event : events) {
      auto groupedTracksP = barrelReducedTrackSelP->sliceByCached(aod::reducedtrack::reducedeventId, event.globalIndex());
      auto groupedTracksN = barrelReducedTrackSelN->sliceByCached(aod::reducedtrack::reducedeventId, event.globalIndex());
      runDalitzSelection<gkReducedEventFillMap,gkReducedTrackFillMap>(event, groupedTracksP, groupedTracksN);
    }
    
    for (auto& track : tracks) {
      dalitzbits(dalitzmap[track.globalIndex()]);
    }
  }*/
  
  void processDummy(soa::Filtered<soa::Join<MyEvents, o2::aod::EventSelectedDal>> const& collisions, soa::Join<MyBarrelTracks, o2::aod::BarrelSelectedDal> const& tracks)
  {

  }

  PROCESS_SWITCH(dalitzPairing, processFullTracks, "Run Dalitz selection on AO2D tables", false);
  //PROCESS_SWITCH(dalitzPairing, processReducedTracks, "Run Dalitz selection on skimmed tables", false);
  PROCESS_SWITCH(dalitzPairing, processDummy, "Do nothing", false);


};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dalitzTrackSelection>(cfgc),
    adaptAnalysisTask<dalitzPairing>(cfgc)
  };
}



