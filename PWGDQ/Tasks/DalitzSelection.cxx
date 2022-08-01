// This task selects electrons from Dalitz decay and writes bits on every track for selection/rejection in later analysis
// Input are: names of the electron cuts, names of the pair cuts, event selection
// For every cut, form every possible electron pair and tag both electrons if the pair passes the pair cut
// produces pair and single track qa plots for every (electron cut, pair cut) combinaison
// For example, we can choose max 4 pair cut and 8 electron cut, the associated bit would be 8*(NpairCut) + (NelectronCut) + 32 (if we are in the analysis stage)
// Could also allow for less cuts if we want more space for other tags
// Should not erase the former bits when running on skimmed table (since the first bits might have used looser electron cuts which are no more available in the reduced tracks)
// Still need a way to remember which bit corresponds to which cut in later analysis
// Should try not to repeat twice the track cuts (once here and once in the table maker)

// Add an option to have different cuts for pairing ? 


#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/Core/PID/PIDResponse.h"
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
using std::array;

using FullTracksExt = soa::Join<aod::FullTracks, aod::TracksCov, aod::TracksExtra, aod::TracksExtended,
                                aod::pidTPCFullEl, aod::pidTPCFullPi,
                                aod::pidTPCFullKa, aod::pidTPCFullPr,
                                aod::pidTOFFullEl, aod::pidTOFFullPi,
                                aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyEventsExt = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;
using BarrelTracksExt = soa::Join<aod::ReducedTracks, aod::ReducedTracksBarrel, aod::ReducedTracksBarrelCov, aod::ReducedTracksBarrelPID>;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::BC | VarManager::ObjTypes::Collision;
constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID;


struct dalitzSelection {

  Produces<o2::aod::DalitzBits> dalitzbits;

  //Configurables
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Dalitz track selection cuts, separated by a comma"};
  Configurable<std::string> fConfigPairCuts{"cfgPairCuts", "", "Dalitz pair selection cuts"};

  AnalysisCompositeCut* fEventCut;              //! Event selection cut
  std::vector<AnalysisCompositeCut> fTrackCuts; //! Barrel track cuts
  std::vector<AnalysisCompositeCut> fPairCuts; //! Dalitz pair cuts

	Partition<Tracks> tracksP = track::sign > 0;
	Partition<Tracks> tracksN = track::sign < 0;


  void init(o2::framework::InitContext& context)
  {
    DefineCuts();

    VarManager::SetDefaultVarNames();
    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    // Create the histogram class names to be added to the histogram manager
    TString histClasses = "";

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

			// No need for event qa since it will be done in the table maker

      if (classStr.Contains("Track")) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track","kine,dca,its,tpcpid,tofpid");
      }
      
      if (classStr.Contains("Pair")) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "pair_barrel","vertexing-barrel");
      }

    }

    // Dalitz selection statistics: one bin for each (track,pair) selection
    TH1I* histTracks = new TH1I("TrackStats", "Dalitz selection statistics", fTrackCuts.size()*fPairCuts.size(), -0.5, fTrackCuts.size()*fPairCuts.size() - 0.5 );
    ib = 1;
    for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, ib++) {
      histTracks->GetXaxis()->SetBinLabel(ib, (*cut).GetName());
    }

  }



  template <uint32_t TEventFillMap, uint32_t TTrackFillMap, bool processReduced, typename TEvent, typename TTracks>
  void runDalitzSelection(TEvent const& event, TTracks const& tracks)
  {
    TPairType = VarManager::kJpsiToEE;
		if(fTrackCuts.size()==0 || fPairCuts.size()==0 ) return;

		// Event selection (no need to select Dalitz on unselected events)
	  VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
    VarManager::FillEvent<TEventFillMap>(event);
    if (!fEventCut->IsSelected(VarManager::fgValues)) {
      return;
    }
    
    std::map<int, uint64_t> dalitzmap;

    uint8_t track1TempFilterMap = 0;
    uint8_t track2TempFilterMap = 0;
    uint8_t twoTracksFilterMap = 0;

		for (auto& [track1,track2] : combinations(tracksP(tracks),tracksN(tracks))) { // loops on combinations +-
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
      VarManager::FillPairVertexing<TPairType, TEventFillMap, TTrackFillMap>(event, track1, track2);
      
      // pair cuts
      int iT = 0;
      for (auto cutT = fTrackCuts.begin(); cutT != fTrackCuts.end(); cutT++, iT++) {
       	if (!(twoTracksFilterMap & (uint8_t(1) << iT))) continue;
       	int iP = 0;
      	for (auto cutP = fPairCuts.begin(); cutP != fPairCuts.end(); cutP++, iP++) {
        	if ((*cutP).IsSelected(VarManager::fgValues)) {// tag the tracks and fill hists
        		
        		int bitNo = 8*iP+iT+32*processReduced;
        		
						fHistMan->FillHistClass(Form("Pair_%s_%s;", cutT.GetName(), cutP.GetName()), VarManager::fgValues);   
						     		
        		bool b1 = track1.dalitzBits() & (uint64(1) << bitNo);
        		if (!b1) {//avoid double counting in histograms
        		  VarManager::FillTrack<TTrackFillMap>(track1);
        			dalitzmap[track1.globalIndex()] |= (uint64(1) << bitNo);
        			fHistMan->FillHistClass(Form("TrackBarrel_%s_%s;", cutT.GetName(), cutP.GetName()), VarManager::fgValues);
        			((TH1I*) histTracks)->Fill(bitNo);
        		}
        		
        		bool b2 = track2.dalitzBits() & (uint64(1) << bitNo);
        		if (!b2) {//avoid double counting in histograms
        		  VarManager::FillTrack<TTrackFillMap>(track2);
        			dalitzmap[track2.globalIndex()] |= (uint64(1) << bitNo);
        			fHistMan->FillHistClass(Form("TrackBarrel_%s_%s;", cutT.GetName(), cutP.GetName()), VarManager::fgValues);
        			((TH1I*) histTracks)->Fill(bitNo);
        		} 
        		
        	}      
      	} //end of pair cuts
      } //end of track cuts
      	        
		} //end of tracksP,N loop
		
		for (auto& track : tracks) {
      dalitzbits(dalitzmap[track.globalIndex()]);
    } 

  }

  void processFullTracks(aod::Collision::iterator const& event, FullTracksExt const& tracks)
  {
    runDalitzSelection<gkEventFillMap,gkTrackFillMapWithCov,false>(event,tracks);
  }
  void processReducedTracks(MyEventsExt::iterator const& event, BarrelTracksExt const& tracks)
  {
    runDalitzSelection<gkEventFillMap,gkTrackFillMapWithCov,true>(event,tracks);
  }

  PROCESS_SWITCH(AnalysisEventSelection, processFullTracks, "Run Dalitz selection on AO2D tables", false);
  PROCESS_SWITCH(AnalysisEventSelection, processReducedTracks, "Run Dalitz selection on skimmed tables", false);


};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<dalitzSelection>(cfgc, TaskName{"dalitz-selector"});
}



