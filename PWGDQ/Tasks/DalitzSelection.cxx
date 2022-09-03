// This task selects electrons from Dalitz decay and writes bits on every track for selection/rejection in later analysis
// Input are: names of the electron cuts, names of the pair cuts
// For every cut, form every possible electron pair and tag both electrons if the pair passes the pair cut
// produces pair and single track qa plots for every (electron cut, pair cut) combination



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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;
using namespace o2::soa;
using std::array;

using FullTracksExt = soa::Join<aod::FullTracks, aod::TracksCov, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullMu,
                                aod::pidTPCFullKa, aod::pidTPCFullPr,
                                aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullMu,
                                aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyReducedEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using BarrelReducedTracksWithCov = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID>;
using BarrelReducedTracks = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelPID>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::Collision;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID;

constexpr static uint32_t gkReducedEventFillMapWithCov = VarManager::ObjTypes::ReducedEvent | VarManager::ObjTypes::ReducedEventExtended | VarManager::ObjTypes::ReducedEventVtxCov;
constexpr static uint32_t gkReducedTrackFillMap = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelPID;
constexpr static uint32_t gkReducedTrackFillMapWithCov = VarManager::ObjTypes::ReducedTrack | VarManager::ObjTypes::ReducedTrackBarrel | VarManager::ObjTypes::ReducedTrackBarrelCov | VarManager::ObjTypes::ReducedTrackBarrelPID;


struct dalitzSelection {
  Preslice<FullTracksExt> perCollision = aod::track::collisionId;
  Produces<o2::aod::DalitzBits> dalitzbits;

  //Configurables
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Dalitz track selection cuts, separated by a comma"};
  Configurable<std::string> fConfigPairCuts{"cfgPairCuts", "", "Dalitz pair selection cuts"};
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandardNoINT7", "Event selection"};

  AnalysisCompositeCut* fEventCut;              //! Event selection cut
  std::vector<AnalysisCompositeCut> fTrackCuts; //! Barrel track cuts
  std::vector<AnalysisCompositeCut> fPairCuts; //! Dalitz pair cuts

  OutputObj<THashList> fOutputList{"output"}; //! the histogram manager output list
  OutputObj<TList> fStatsList{"Statistics"};  //! skimming statistics
  HistogramManager* fHistMan;

  int NFilledTracks;
  int globalCounter;
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
    histClasses += "Event_AfterCuts;";

    for (auto& cutT : fTrackCuts) {
      for (auto& cutP : fPairCuts) {// single track and pair histograms
        histClasses += Form("TrackBarrel_%s_%s;", cutT.GetName(), cutP.GetName());
        histClasses += Form("Pair_%s_%s;", cutT.GetName(), cutP.GetName());
      }
    }

    DefineHistograms(histClasses);                   // define all histograms
    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
    
    VarManager::SetupTwoProngDCAFitter(5.0f, true, 200.0f, 4.0f, 1.0e-3f, 0.9f, true);

    globalCounter = 0;
  }

  void DefineCuts()
  {
    // Event cuts
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));

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

  void DefineHistograms(TString histClasses)
  {
    std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
    for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
      TString classStr = objArray->At(iclass)->GetName();
      fHistMan->AddHistClass(classStr.Data());

      if (classStr.Contains("Event")) {
        dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "event", "");
      }

      if (classStr.Contains("Track")) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track","kine,dca,tpcpid,tofpid");
      }
      
      if (classStr.Contains("Pair")) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "pair_barrel","vertexing-barrel");
      }

    }


    fStatsList.setObject(new TList());
    fStatsList->SetOwner(kTRUE);

    std::vector<TString> eventLabels{"BCs", "Collisions before filtering", "Before cuts", "After cuts"};
    TH1I* histEvents = new TH1I("EventStats", "Event statistics", eventLabels.size(), -0.5, eventLabels.size() - 0.5);
    int ib = 1;
    for (auto label = eventLabels.begin(); label != eventLabels.end(); label++, ib++) {
      histEvents->GetXaxis()->SetBinLabel(ib, (*label).Data());
    }
    fStatsList->Add(histEvents);

    // Dalitz selection statistics: one bin for each (track,pair) selection
    TH1I* histTracks = new TH1I("TrackStats", "Dalitz selection statistics", fTrackCuts.size()*fPairCuts.size(), -0.5, fTrackCuts.size()*fPairCuts.size() - 0.5 );
    ib = 1;      
    for (auto cutP = fPairCuts.begin(); cutP != fPairCuts.end(); cutP++, ib++) {
      for (auto cutT = fTrackCuts.begin(); cutT != fTrackCuts.end(); cutT++, ib++) {
        histTracks->GetXaxis()->SetBinLabel(ib, Form("%s_%s", (*cutT).GetName(), (*cutP).GetName()));
      }
    }

    fStatsList->Add(histTracks);
  }

  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, bool processReduced, typename TEvent, typename TTracks>
  void runDalitzSelection(TEvent const& event, TTracks const& tracks)
  {
    const int TPairType = VarManager::kJpsiToEE;
   
    VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    VarManager::FillEvent<TEventFillMap>(event);
    ((TH1I*)fStatsList->At(0))->Fill(2.0);

    if (fEventCut->IsSelected(VarManager::fgValues)) {
      ((TH1I*)fStatsList->At(0))->Fill(3.0);
      fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);

      uint8_t track1TempFilterMap = 0;
      uint8_t track2TempFilterMap = 0;
      uint8_t twoTracksFilterMap = 0;
      for (auto& [track1,track2] : o2::soa::combinations(CombinationsUpperIndexPolicy(tracks, tracks))) {
      
        if (track1.sign()*track2.sign() > 0) continue; // loops on combinations +-
        track1TempFilterMap = uint8_t(0);
        track2TempFilterMap = uint8_t(0);
        VarManager::FillTrack<TTrackFillMap>(track1);

        // apply track cuts  to positive track
        int i = 0;
        for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, i++) {
          if ((*cut).IsSelected(VarManager::fgValues)) {
            track1TempFilterMap |= (uint8_t(1) << i);
          }
        } 
        if (!track1TempFilterMap) {
          continue;
        }

        // apply track cuts  to negative track
        i = 0;
        VarManager::FillTrack<TTrackFillMap>(track2);     	
        for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, i++) {
      	  if ((*cut).IsSelected(VarManager::fgValues)) {
      	    track2TempFilterMap |= (uint8_t(1) << i);
      	  }
        }
         
      
        twoTracksFilterMap = track1TempFilterMap & track2TempFilterMap;
        if (!twoTracksFilterMap) { //do not pair if there is no common track filter
          continue;
        }
      	
        // pairing
        VarManager::FillPair<TPairType, TTrackFillMap>(track1, track2);
        if constexpr (static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::TrackCov) || static_cast<bool>(TTrackFillMap & VarManager::ObjTypes::ReducedTrackBarrelCov)) {
          VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, track1, track2);
        }      

        // pair cuts
        int iT = 0;  
        for (auto cutT = fTrackCuts.begin(); cutT != fTrackCuts.end(); cutT++, iT++) {
       	  if (!(twoTracksFilterMap & (uint8_t(1) << iT))) continue;
       	  int iP = 0;
      	  for (auto cutP = fPairCuts.begin(); cutP != fPairCuts.end(); cutP++, iP++) {
            if ((*cutP).IsSelected(VarManager::fgValues)) {

              // tag the tracks and fill hists
              int bitNo = fTrackCuts.size()*iP+iT;
              fHistMan->FillHistClass(Form("Pair_%s_%s", (*cutT).GetName(), (*cutP).GetName()), VarManager::fgValues);   

              bool b1 = dalitzmap[track1.globalIndex()] & (uint32_t(1) << bitNo);
              if (!b1) { //avoid double counting in histograms
                VarManager::FillTrack<TTrackFillMap>(track1);
                dalitzmap[track1.globalIndex()] |= (uint32_t(1) << bitNo);
                fHistMan->FillHistClass(Form("TrackBarrel_%s_%s", (*cutT).GetName(), (*cutP).GetName()), VarManager::fgValues);
                ((TH1I*) fStatsList->At(1))->Fill(bitNo);
              }
              	
              bool b2 = dalitzmap[track2.globalIndex()] & (uint32_t(1) << bitNo);
              if (!b2) { //avoid double counting in histograms
                VarManager::FillTrack<TTrackFillMap>(track2);
                dalitzmap[track2.globalIndex()] |= (uint32_t(1) << bitNo);
                fHistMan->FillHistClass(Form("TrackBarrel_%s_%s", (*cutT).GetName(), (*cutP).GetName()), VarManager::fgValues);
                ((TH1I*) fStatsList->At(1))->Fill(bitNo);
              } 
         	
            } // end if      
      	  } //end of pair cuts
        } //end of track cuts
      	        
      } //end of tracksP,N loop
      
    }
   
  }

  void processFullTracks(MyEvents const& collisions, FullTracksExt const& tracks)
  {
    if(fTrackCuts.size()==0 || fPairCuts.size()==0 ) {
      std::cout<<"WARNING: YOU NEED A TRACK AND A PAIR CUT"<<std::endl;
      return;
    }   

    dalitzmap.clear();

    for (auto& collision : collisions) {
      auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());
      runDalitzSelection<gkEventFillMap,gkTrackFillMapWithCov,false>(collision,groupedTracks);
    }
    
    //Fill dalitz bits
    for (auto& track : tracks) {
      dalitzbits(dalitzmap[track.globalIndex()]);
    } 
  }
  
  /*
 void processReducedTracksWithCov(MyReducedEvents::iterator const& event, BarrelReducedTracksWithCov const& tracks)
  {
    if(fTrackCuts.size()==0 || fPairCuts.size()==0 ) {
      std::cout<<"WARNING: YOU NEED A TRACK AND A PAIR CUT"<<std::endl;
      return;
    }
    dalitzmap.clear();
    
    runDalitzSelection<gkReducedEventFillMapWithCov,gkReducedTrackFillMapWithCov,true>(event,tracks);
    for (auto& track : tracks) {
      dalitzbits(dalitzmap[track.globalIndex()]);
    }
    
  }
 void processReducedTracks(MyReducedEvents::iterator const& event, BarrelReducedTracks const& tracks)
  {
    if(fTrackCuts.size()==0 || fPairCuts.size()==0 ) {
      std::cout<<"WARNING: YOU NEED A TRACK AND A PAIR CUT"<<std::endl;
      return;
    }
    dalitzmap.clear();
    runDalitzSelection<gkReducedEventFillMapWithCov,gkReducedTrackFillMap,true>(event,tracks);
    
    for (auto& track : tracks) {
      dalitzbits(dalitzmap[track.globalIndex()]);
    }
  }*/

  PROCESS_SWITCH(dalitzSelection, processFullTracks, "Run Dalitz selection on AO2D tables", false);
  //PROCESS_SWITCH(dalitzSelection, processReducedTracks, "Run Dalitz selection on skimmed tables", false);
  //PROCESS_SWITCH(dalitzSelection, processReducedTracksWithCov, "Run Dalitz selection on skimmed tables with track cov information stored", false);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dalitzSelection>(cfgc)};
}



