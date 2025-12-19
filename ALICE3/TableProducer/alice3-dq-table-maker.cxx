// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
//
// Contact: Ionut Cristian Arsene iarsene@cern.ch, i.c.arsene@fys.uio.no
//          Alexander Tiekoetter (alexander.tiekoetter@cern.ch)
/// \file alice3-dq-table-maker.cxx
/// \brief DQ table maker for ALICE 3

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "PWGDQ/Core/MuonMatchingMlResponse.h"
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGDQ/DataModel/ReducedTablesAlice3.h"

#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFPIDTrk.h"
#include "ALICE3/DataModel/collisionAlice3.h"
#include "ALICE3/DataModel/tracksAlice3.h"

#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsGlobalTracking/RecoContainer.h"
#include "DataFormatsGlobalTracking/RecoContainerCreateTracksVariadic.h"
#include "DataFormatsITSMFT/ROFRecord.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsVertexing/PVertexerParams.h"
#include "DetectorsVertexing/VertexTrackMatcher.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/DataTypes.h"
#include "Framework/runDataProcessing.h"
#include "MathUtils/Primitive2D.h"
#include "ReconstructionDataFormats/PrimaryVertex.h"
#include "ReconstructionDataFormats/VtxTrackIndex.h"
#include "ReconstructionDataFormats/VtxTrackRef.h"

#include "TList.h"
#include "THashList.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, 
                                 aod::TracksCov, aod::TracksAlice3, aod::TracksExtraA3, 
                                 aod::UpgradeTofs, aod::UpgradeRichs, aod::UpgradeRichSignals,
                                 aod::UpgradeTrkPids, aod::UpgradeTrkPidSignals, 
                                 aod::McTrackLabels>;

using MyEvents = soa::Join<aod::Collisions, aod::CollisionsAlice3, aod::McCollisionLabels>;
using MyEventsMC = aod::McCollisions;

constexpr static uint32_t gkEventFillMap = VarManager::ObjTypes::Collision;
constexpr static uint32_t gkEventMcFillMap = VarManager::ObjTypes::CollisionMC;

constexpr static uint32_t gkTrackFillMapWithCov = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID;

struct Alice3DQTableMaker {
  
  Produces<ReducedA3MCEvents> eventMC;
  Produces<ReducedA3MCTracks> trackMC;

  Produces<ReducedA3Events> event;
  Produces<ReducedA3EventsVtxCov> eventVtxCov;
  Produces<ReducedA3EventsInfo> eventInfo;
  Produces<ReducedA3MCEventLabels> eventMClabels;

  Produces<ReducedA3TracksBarrelInfo> trackBarrelInfo;
  Produces<ReducedA3Tracks> trackBasic;
  Produces<ReducedA3TracksBarrel> trackBarrel;
  Produces<ReducedA3TracksBarrelCov> trackBarrelCov;
  Produces<ReducedA3TracksAssoc> trackBarrelAssoc;
  Produces<ReducedA3TracksBarrelLabels> trackBarrelLabels;
  
  Produces<ReducedA3PIDTOF> trackPIDTOF;
  Produces<ReducedA3PIDRich> trackPIDRich;
  Produces<ReducedA3PIDRichSignals> trackPIDRichSig;
  Produces<ReducedA3PIDOT> trackPIDOT;
  

  OutputObj<THashList> fOutputList{"output"};
  OutputObj<TList> fStatsList{"Statistics"}; //! skimming statistics

  HistogramManager* fHistMan;

  // Event and track AnalysisCut configurables
  struct : ConfigurableGroup {
    Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "", "Event selection"};
    Configurable<std::string> fConfigTrackCuts{"cfgBarrelTrackCuts", "", "barrel track cut"};
    Configurable<std::string> fConfigEventCutsJSON{"cfgEventCutsJSON", "", "Additional event selection in JSON format"};
    Configurable<std::string> fConfigTrackCutsJSON{"cfgBarrelTrackCutsJSON", "", "Additional list of barrel track cuts in JSON format"};
  } fConfigCuts;

  // MC signals to be skimmed
  Configurable<std::string> fConfigMCSignals{"cfgMCsignals", "", "Comma separated list of MC signals"};
  Configurable<std::string> fConfigMCSignalsJSON{"cfgMCsignalsJSON", "", "Additional list of MC signals via JSON"};

  // Steer QA output
  struct : ConfigurableGroup {
    Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
    Configurable<bool> fConfigDetailedQA{"cfgDetailedQA", false, "If true, include more QA histograms (BeforeCuts classes)"};
    Configurable<std::string> fConfigAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> fConfigAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> fConfigAddMuonHistogram{"cfgAddMuonHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> fConfigAddMCTruthHistogram{"cfgAddMCTruthHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> fConfigAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};
  } fConfigHistOutput;

  AnalysisCompositeCut* fEventCut;               //! Event selection cut
  std::vector<AnalysisCompositeCut*> fTrackCuts; //! Barrel track cuts
  std::vector<AnalysisCompositeCut*> fMuonCuts;  //! Muon track cuts

  bool fDoDetailedQA = false;

  std::vector<MCSignal*> fMCSignals;
  std::map<uint64_t, int> fLabelsMap;
  std::map<uint64_t, int> fLabelsMapReversed;
  std::map<uint64_t, uint16_t> fMCFlags;
  std::map<uint32_t, uint32_t> fCollIndexMap;             // key: old collision index, value: skimmed collision index
  std::map<uint32_t, uint32_t> fTrackIndexMap;            // key: old track global index, value: new track global index

  void init(InitContext& context)
  {
    bool isProcessSkimmingEnabled = context.mOptions.get<bool>("processSkimming");

    if(!isProcessSkimmingEnabled) 
      LOG(fatal) << "No process function was enabled ALICE 3 TableMaker";
    
    VarManager::SetDefaultVarNames(); // Important that this is called before DefineCuts() !!!

    DefineCuts();

    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(kTRUE);
    fHistMan->SetDefaultVarNames(VarManager::fgVariableNames, VarManager::fgVariableUnits);

    if (fConfigHistOutput.fConfigQA && fConfigHistOutput.fConfigDetailedQA) {
      fDoDetailedQA = true;
    }

    TString histClasses = "";

    if (fDoDetailedQA) {
      histClasses += "Event_BeforeCuts;";
    }
    
    if (fConfigHistOutput.fConfigQA) {
      histClasses += "Event_AfterCuts;";
      histClasses += "Event_MCTruth;";
    }

    if(isProcessSkimmingEnabled) {
      if (fDoDetailedQA) {
        histClasses += "TrackBarrel_BeforeCuts;";
      }
      
      if (fConfigHistOutput.fConfigQA) {
        for (auto& cut : fTrackCuts) {
          histClasses += Form("TrackBarrel_%s;", cut->GetName());
        }
      }
    }

    TString configNamesStr = fConfigMCSignals.value;
    std::unique_ptr<TObjArray> objArray(configNamesStr.Tokenize(","));

    if (objArray->GetEntries() > 0) {
      for (int isig = 0; isig < objArray->GetEntries(); ++isig) {
        MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objArray->At(isig)->GetName());
        if (sig) {
          fMCSignals.push_back(sig);
        }
      }
    }

    TString addMCSignalsStr = fConfigMCSignalsJSON.value;

    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());

      for (auto& mcIt : addMCSignals) {
        if (mcIt) {
          fMCSignals.push_back(mcIt);
        }
      }
    }

    for(auto& mcIt : fMCSignals) {
      if (fConfigHistOutput.fConfigQA) {
        histClasses += Form("MCTruth_%s;", mcIt->GetName());
      }
      if (fDoDetailedQA) {
        if (isProcessSkimmingEnabled) {
          for (auto& cut : fTrackCuts) {
            histClasses += Form("TrackBarrel_%s_%s;", cut->GetName(), mcIt->GetName());
          }
        }
      }
    }

    DefineHistograms(histClasses);
    
    TString addHistsStr = fConfigHistOutput.fConfigAddJSONHistograms.value;
    if (fConfigHistOutput.fConfigQA && addHistsStr != "") 
    {
      dqhistograms::AddHistogramsFromJSON(fHistMan, addHistsStr.Data());
    }

    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  void DefineCuts()
  {
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigCuts.fConfigEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));

    TString addEvCutsStr = fConfigCuts.fConfigEventCutsJSON.value;
    if (addEvCutsStr != "") {
      std::vector<AnalysisCut*> addEvCuts = dqcuts::GetCutsFromJSON(addEvCutsStr.Data());
      for (auto& cutIt : addEvCuts) {
        fEventCut->AddCut(cutIt);
      }
    }

    // Barrel track cuts
    TString cutNamesStr = fConfigCuts.fConfigTrackCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    // Additional Barrel track cuts via JSON
    TString addTrackCutsStr = fConfigCuts.fConfigTrackCutsJSON.value;
    if (addTrackCutsStr != "") {
      std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(addTrackCutsStr.Data());
      for (auto& t : addTrackCuts) {
        fTrackCuts.push_back(reinterpret_cast<AnalysisCompositeCut*>(t));
      }
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

  void DefineHistograms(TString histClasses)
  {
    std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
    for (Int_t iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
      TString classStr = objArray->At(iclass)->GetName();
      if (fConfigHistOutput.fConfigQA) {
        fHistMan->AddHistClass(classStr.Data());
      }

      TString histEventName = fConfigHistOutput.fConfigAddEventHistogram.value;
      if (classStr.Contains("Event")) {
        if (fConfigHistOutput.fConfigQA && !classStr.Contains("MCTruth")) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "event", histEventName);
        } else {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "event", "generator");
        }
      }

      TString histTrackName = fConfigHistOutput.fConfigAddTrackHistogram.value;
      if (classStr.Contains("Track")) {
        if (fConfigHistOutput.fConfigQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", histTrackName);
        }
      }

      TString histMCTruthName = fConfigHistOutput.fConfigAddMCTruthHistogram.value;
      if (classStr.Contains("MCTruth") && !classStr.Contains("Event")) {
        if (fConfigHistOutput.fConfigQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "mctruth_track", histMCTruthName);
        }
      }
    }

    // create statistics histograms (event, tracks, muons, MCsignals)
    fStatsList.setObject(new TList());
    fStatsList->SetOwner(kTRUE);
    std::vector<TString> eventLabels{"Collisions before filtering", "Before cuts", "After cuts"};
    TH2I* histEvents = new TH2I("EventStats", "Event statistics", eventLabels.size(), -0.5, eventLabels.size() - 0.5, o2::aod::evsel::kNsel + 1, -0.5, (float)o2::aod::evsel::kNsel + 0.5);
    int ib = 1;
    for (auto label = eventLabels.begin(); label != eventLabels.end(); label++, ib++) {
      histEvents->GetXaxis()->SetBinLabel(ib, (*label).Data());
    }
    for (int ib = 1; ib <= o2::aod::evsel::kNsel; ib++) {
      histEvents->GetYaxis()->SetBinLabel(ib, o2::aod::evsel::selectionLabels[ib - 1]);
    }
    histEvents->GetYaxis()->SetBinLabel(o2::aod::evsel::kNsel + 1, "Total");
    fStatsList->Add(histEvents);

    // Track statistics: one bin for each track selection and 5 bins for V0 tags (gamma, K0s, Lambda, anti-Lambda, Omega)
    TH1I* histTracks = new TH1I("TrackStats", "Track statistics", fTrackCuts.size() + 5.0, -0.5, fTrackCuts.size() - 0.5 + 5.0);
    ib = 1;
    for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, ib++) {
      histTracks->GetXaxis()->SetBinLabel(ib, (*cut)->GetName());
    }
    const char* v0TagNames[5] = {"Photon conversion", "K^{0}_{s}", "#Lambda", "#bar{#Lambda}", "#Omega"};
    for (int ib = 0; ib < 5; ib++) {
      histTracks->GetXaxis()->SetBinLabel(fTrackCuts.size() + 1 + ib, v0TagNames[ib]);
    }
    fStatsList->Add(histTracks);

    TH1I* histMCsignals = new TH1I("MCsignals", "MC signals", fMCSignals.size() + 1, -0.5, fMCSignals.size() - 0.5 + 1.0);
    ib = 1;
    for (auto signal = fMCSignals.begin(); signal != fMCSignals.end(); signal++, ib++) {
      histMCsignals->GetXaxis()->SetBinLabel(ib, (*signal)->GetName());
    }
    histMCsignals->GetXaxis()->SetBinLabel(fMCSignals.size() + 1, "Others (matched to reco tracks)");
    fStatsList->Add(histMCsignals);
  }

  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  void skimMCCollisions(MyEventsMC const& mcCollisions)
  {
    // skim MC collisions
    // NOTE: So far, all MC collisions are skimmed. In case there will be filtering based on MC collisions,
    //       one has to do a mapping of the old vs new indices so that the skimmed labels are properly updated.
    VarManager::ResetValues(0, VarManager::kNVars);

    for (auto& mcCollision : mcCollisions) {
      VarManager::FillEventAlice3<gkEventMcFillMap>(mcCollision);

      fHistMan->FillHistClass("Event_MCTruth", VarManager::fgValues);

      eventMC(mcCollision.generatorsID(), mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(),
              mcCollision.t(), mcCollision.weight(), mcCollision.impactParameter()); // TODO: Determine and fill multiplicity values
    }
  }

  void skimMCParticles(aod::McParticles const& mcTracks, MyEventsMC const&)
  {
    // Select MC particles which fulfill at least one of the user specified MC signals
    // In this function we just fill a map with the labels of selected particles, not creating the tables themselves.
    //  The reason is that in the skims we will additionally add any MC label connected to selected reconstructed tracks
    //      which were not selected already via the MC signals

    fLabelsMap.clear();
    fLabelsMapReversed.clear();
    fMCFlags.clear();

    uint16_t mcflags = static_cast<uint16_t>(0); // flags which will hold the decisions for each MC signal
    int trackCounter = 0;

    for (auto& mctrack : mcTracks) {
      // check all the requested MC signals and fill the decision bit map
      mcflags = 0;
      int i = 0;
      for (auto& sig : fMCSignals) {
        bool checked = false;
        if constexpr (soa::is_soa_filtered_v<aod::McParticles>) {
          auto mctrack_raw = mcTracks.rawIteratorAt(mctrack.globalIndex());
          checked = sig->CheckSignal(true, mctrack_raw);
        } else {
          checked = sig->CheckSignal(true, mctrack);
        }
        if (checked) {
          mcflags |= (static_cast<uint16_t>(1) << i);
        }
        i++;
      }

      // if no MC signals were matched, continue
      if (mcflags == 0) {
        continue;
      }

      // If this MC track was not already added to the map, add it now
      if (fLabelsMap.find(mctrack.globalIndex()) == fLabelsMap.end()) {
        fLabelsMap[mctrack.globalIndex()] = trackCounter;
        fLabelsMapReversed[trackCounter] = mctrack.globalIndex();
        fMCFlags[mctrack.globalIndex()] = mcflags;
        trackCounter++;

        // fill histograms for each of the signals, if found
        if (fConfigHistOutput.fConfigQA) {
          VarManager::FillTrackMC(mcTracks, mctrack);
          auto mcCollision = mctrack.template mcCollision_as<MyEventsMC>();
          VarManager::FillEvent<gkEventMcFillMap>(mcCollision);
          int j = 0;
          for (auto signal = fMCSignals.begin(); signal != fMCSignals.end(); signal++, j++) {
            if (mcflags & (static_cast<uint16_t>(1) << j)) {
              fHistMan->FillHistClass(Form("MCTruth_%s", (*signal)->GetName()), VarManager::fgValues);
            }
          }
        }
      }
    }
  }

  void skimCollisions(MyEvents const& collisions)
  {
    // Skim reconstructed collisions which are selected by the user specified cuts
    // Create a collision index map to relate between the "old" AO2D indices and the skimmed ones
    fCollIndexMap.clear();

    // Loop over collisions
    for (const auto& collision : collisions) {
      
      (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(o2::aod::evsel::kNsel));

      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEventAlice3<gkEventFillMap>(collision); // extract event information and place it in the fValues array

      if (collision.has_mcCollision()) {
        auto mcCollision = collision.template mcCollision_as<MyEventsMC>();
        VarManager::FillEventAlice3<gkEventMcFillMap>(mcCollision);
      }
      
      if (fDoDetailedQA) {
        fHistMan->FillHistClass("Event_BeforeCuts", VarManager::fgValues);
      }

      // Apply the user specified event selection
      if (!fEventCut->IsSelected(VarManager::fgValues)) {
        continue;
      }

      (reinterpret_cast<TH2I*>(fStatsList->At(0)))->Fill(3.0, static_cast<float>(o2::aod::evsel::kNsel));

      // Fill historams after event cuts
      fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);

      event(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(), 
            collision.collisionTime(), collision.collisionTimeRes(), collision.multDensity());
      
      eventVtxCov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());
      eventMClabels(collision.mcCollisionId(), collision.mcMask());
      eventInfo(collision.globalIndex());

      // add an element for this collision into the map
      fCollIndexMap[collision.globalIndex()] = event.lastIndex();
    }
  }

  void skimTracks(MyEvents::iterator const& collision, MyBarrelTracks const& /*tracks*/, TrackAssoc const& assocs, aod::McParticles const& mcTracks)
  {
    // Skim the barrel track associations
    // Apply track cuts for each collision association and if it passes the cuts, we skim it.
    // NOTE: If selection cuts include conditions on quantities dependent on the associated collision (e.g. DCA),
    //         one track may pass for some association and fail for others.
    //       Tracks are written only once in the skims, even if they contribute to more than one association
    //         so in case of multiple associations, the variables depending on the collision association (e.g. DCA, secondary vertexing, etc)
    //         have to be recomputed at analysis time for each association.

    uint64_t trackFilteringTag = static_cast<uint64_t>(0);
    uint32_t trackTempFilterMap = static_cast<uint32_t>(0);
    uint16_t mcflags = static_cast<uint16_t>(0);
    int trackCounter = fLabelsMap.size();

    for (const auto& assoc : assocs) {

      auto track = assoc.template track_as<MyBarrelTracks>();
      
      if (fCollIndexMap.find(track.collisionId()) == fCollIndexMap.end()) {
        continue;
      }

      trackFilteringTag = static_cast<uint64_t>(0);
      trackTempFilterMap = static_cast<uint32_t>(0);

      // Compute track quantities and fill histograms
      VarManager::FillTrackAlice3<gkTrackFillMapWithCov>(track);

      if (fDoDetailedQA) {
        fHistMan->FillHistClass("TrackBarrel_BeforeCuts", VarManager::fgValues);
      }

      int i = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, i++) {
        if ((*cut)->IsSelected(VarManager::fgValues)) {
          trackTempFilterMap |= (static_cast<uint32_t>(1) << i);
          if (fConfigHistOutput.fConfigQA) {
            fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut)->GetName()), VarManager::fgValues);
          }
          (reinterpret_cast<TH1I*>(fStatsList->At(1)))->Fill(static_cast<float>(i));
        }
      }
      if (!trackTempFilterMap) {
        continue;
      }

      // If this track is already present in the index map, it means it was already skimmed,
      // so we just store the association and we skip the track
      if (fTrackIndexMap.find(track.globalIndex()) != fTrackIndexMap.end()) {
        trackBarrelAssoc(fCollIndexMap[collision.globalIndex()], fTrackIndexMap[track.globalIndex()]);
        continue;
      }

      trackFilteringTag |= (static_cast<uint64_t>(trackTempFilterMap) << VarManager::kBarrelUserCutsBits); // BIT13-...:  user track filters

      // NOTE: The collision ID that is written in the table is the one originally assigned in the AOD.
      //       However, in data analysis one should loop over associations, so this one should not be used.
      //      In the case of Run2-like analysis, there will be no associations, so this ID will be the one originally assigned in the AO2Ds (updated for the skims)
      uint32_t reducedEventIdx = fCollIndexMap[track.collisionId()];

      // NOTE: trackBarrelInfo stores the index of the collision as in AO2D (for use in some cases where the analysis on skims is done
      //   in workflows where the original AO2Ds are also present)
      //trackBarrelInfo(track.collisionId(), collision.posX(), collision.posY(), collision.posZ(), track.globalIndex());
      trackBasic(reducedEventIdx, trackFilteringTag, track.pt(), track.eta(), track.phi(), track.sign(), 0);

      trackBarrel(track.x(), track.alpha(), track.y(), track.z(), track.snp(), track.tgl(), track.signed1Pt(),
                  track.flags(), track.itsClusterMap(), track.itsChi2NCl(),
                  track.isReconstructed(), track.nSiliconHits(), track.nTPCHits(), track.length(), track.dcaXY(), 
                  track.dcaZ());

      if constexpr (static_cast<bool>(gkTrackFillMapWithCov & VarManager::ObjTypes::TrackCov)) {
        trackBarrelCov(track.cYY(), track.cZY(), track.cZZ(), track.cSnpY(), track.cSnpZ(),
                       track.cSnpSnp(), track.cTglY(), track.cTglZ(), track.cTglSnp(), track.cTglTgl(),
                       track.c1PtY(), track.c1PtZ(), track.c1PtSnp(), track.c1PtTgl(), track.c1Pt21Pt2());
      }

      if constexpr (static_cast<bool>(gkTrackFillMapWithCov & VarManager::ObjTypes::TrackPID)) {
        
        trackPIDTOF(track.tofEventTime(), track.tofEventTimeErr(),
                    track.nSigmaElectronInnerTOF(), track.nSigmaMuonInnerTOF(), track.nSigmaPionInnerTOF(),
                    track.nSigmaKaonInnerTOF(), track.nSigmaProtonInnerTOF(), track.nSigmaDeuteronInnerTOF(), 
                    track.nSigmaTritonInnerTOF(), track.nSigmaHelium3InnerTOF(), track.nSigmaAlphaInnerTOF(),
                    track.innerTOFTrackTimeReco(), track.innerTOFTrackLengthReco(),
                    track.nSigmaElectronOuterTOF(), track.nSigmaMuonOuterTOF(), track.nSigmaPionOuterTOF(),
                    track.nSigmaKaonOuterTOF(), track.nSigmaProtonOuterTOF(), track.nSigmaDeuteronOuterTOF(), 
                    track.nSigmaTritonOuterTOF(), track.nSigmaHelium3OuterTOF(), track.nSigmaAlphaOuterTOF(),
                    track.outerTOFTrackTimeReco(), track.outerTOFTrackLengthReco());

        trackPIDRich(track.nSigmaElectronRich(), track.nSigmaMuonRich(), track.nSigmaPionRich(),
                     track.nSigmaKaonRich(), track.nSigmaProtonRich(), track.nSigmaDeuteronRich(),
                     track.nSigmaTritonRich(), track.nSigmaHelium3Rich(), track.nSigmaAlphaRich());
        
        trackPIDRichSig(track.hasSig(), track.hasSigInGas(), 
                        track.hasSigEl(), track.hasSigMu(), track.hasSigPi(),
                        track.hasSigKa(), track.hasSigPr(), track.hasSigDe(),
                        track.hasSigTr(), track.hasSigHe3(), track.hasSigAl());
        
        trackPIDOT(track.timeOverThresholdBarrel(),
                   track.nSigmaTrkEl(), track.nSigmaTrkMu(), track.nSigmaTrkPi(),
                   track.nSigmaTrkKa(), track.nSigmaTrkPr(), track.nSigmaTrkDe(),
                   track.nSigmaTrkTr(), track.nSigmaTrkHe(), track.nSigmaTrkAl());
      }

      fTrackIndexMap[track.globalIndex()] = trackBasic.lastIndex();

      // Check whether the MCParticle corresponding to this reconstructed track was already selected for skimming
      // If not, add it to the skimming map
      if (!track.has_mcParticle()) {
        trackBarrelLabels(-1, 0, 0); // this is the case when there is no matched MCParticle
      } else {
        auto mctrack = track.template mcParticle_as<aod::McParticles>();
        VarManager::FillTrackMC(mcTracks, mctrack);

        mcflags = 0;
        int i = 0; // runs over the MC signals
        int j = 0; // runs over the track cuts
        // check all the specified signals and fill histograms for MC truth matched tracks
        for (auto& sig : fMCSignals) {
          if (sig->CheckSignal(true, mctrack)) {
            mcflags |= (static_cast<uint16_t>(1) << i);
            // If detailed QA is on, fill histograms for each MC signal and track cut combination
            if (fDoDetailedQA) {
              j = 0;
              for (auto& cut : fTrackCuts) {
                if (trackTempFilterMap & (uint8_t(1) << j)) {
                  fHistMan->FillHistClass(Form("TrackBarrel_%s_%s", cut->GetName(), sig->GetName()), VarManager::fgValues); // fill the reconstructed truth
                }
                j++;
              }
            }
          }
          i++;
        }

        // if the MC truth particle corresponding to this reconstructed track is not already written,
        //   add it to the skimmed stack
        if (!(fLabelsMap.find(mctrack.globalIndex()) != fLabelsMap.end())) {
          fLabelsMap[mctrack.globalIndex()] = trackCounter;
          fLabelsMapReversed[trackCounter] = mctrack.globalIndex();
          fMCFlags[mctrack.globalIndex()] = mcflags;
          trackCounter++;
        }

        trackBarrelLabels(fLabelsMap.find(mctrack.globalIndex())->second, track.mcMask(), mcflags);
      }

      trackBarrelAssoc(fCollIndexMap[collision.globalIndex()], fTrackIndexMap[track.globalIndex()]);
    } // end loop over associations
  } // end skimTracks

  Preslice<aod::McParticles_001> perMcCollision = aod::mcparticle::mcCollisionId;
  Preslice<MyBarrelTracks> perCollisionTracks = aod::track::collisionId;

  void fullSkimming(MyEvents const& collisions,
                    MyBarrelTracks const& tracksBarrel, aod::TrackAssoc const& trackAssocs, 
                    MyEventsMC const& mcCollisions, aod::McParticles const& mcParticles)
  {
    eventMC.reserve(mcCollisions.size());
    skimMCCollisions(mcCollisions);

    // skim collisions
    event.reserve(collisions.size());
    eventVtxCov.reserve(collisions.size());
    eventMClabels.reserve(collisions.size());
    eventInfo.reserve(collisions.size());
    
    skimCollisions(collisions);

    if (fCollIndexMap.size() == 0) 
      return;

    skimMCParticles(mcParticles, mcCollisions);

    if constexpr (static_cast<bool>(gkTrackFillMapWithCov)) {
      fTrackIndexMap.clear();
      trackBasic.reserve(tracksBarrel.size());
      trackBarrel.reserve(tracksBarrel.size());
      trackBarrelCov.reserve(tracksBarrel.size());
      trackPIDTOF.reserve(tracksBarrel.size());
      trackPIDRich.reserve(tracksBarrel.size());
      trackPIDRichSig.reserve(tracksBarrel.size());
      trackPIDOT.reserve(tracksBarrel.size());
      trackBarrelAssoc.reserve(tracksBarrel.size());
      trackBarrelLabels.reserve(tracksBarrel.size());
    }

    if (fCollIndexMap.size() > 0) {

      for (auto const& [origIdx, skimIdx] : fCollIndexMap) 
      {
        auto collision = collisions.rawIteratorAt(origIdx);

        if constexpr (static_cast<bool>(gkTrackFillMapWithCov)) {
          auto groupedTrackIndices = trackAssocs.sliceBy(trackIndicesPerCollision, origIdx);

          skimTracks(collision, tracksBarrel, groupedTrackIndices, mcParticles);
        }
      }
    }

    // Loop over the label map, create the mother/daughter relationships if these exist and write the skimmed MC stack
    for (const auto& [newLabel, oldLabel] : fLabelsMapReversed) {
      auto mctrack = mcParticles.iteratorAt(oldLabel);
      uint16_t mcflags = fMCFlags.find(oldLabel)->second;

      std::vector<int> mothers;
      if (mctrack.has_mothers()) {
        for (auto& m : mctrack.mothersIds()) {
          if (m < mcParticles.size()) { // protect against bad mother indices
            if (fLabelsMap.find(m) != fLabelsMap.end()) {
              mothers.push_back(fLabelsMap.find(m)->second);
            }
          } else {
            cout << "Mother label (" << m << ") exceeds the McParticles size (" << mcParticles.size() << ")" << endl;
            cout << " Check the MC generator" << endl;
          }
        }
      }

      // TODO: Check that the daughter slice in the skimmed table works as expected
      //       Note that not all daughters from the original table are preserved in the skimmed MC stack
      std::vector<int> daughters;
      if (mctrack.has_daughters()) {
        for (int d = mctrack.daughtersIds()[0]; d <= mctrack.daughtersIds()[1]; ++d) {
          // TODO: remove this check as soon as issues with MC production are fixed
          if (d < mcParticles.size()) { // protect against bad daughter indices
            if (fLabelsMap.find(d) != fLabelsMap.end()) {
              daughters.push_back(fLabelsMap.find(d)->second);
            }
          } else {
            cout << "Daughter label (" << d << ") exceeds the McParticles size (" << mcParticles.size() << ")" << endl;
            cout << " Check the MC generator" << endl;
          }
        }
      }
      int daughterRange[2] = {-1, -1};
      if (daughters.size() > 0) {
        daughterRange[0] = daughters[0];
        daughterRange[1] = daughters[daughters.size() - 1];
      }

      // NOTE: Here we assume that MC collisions are not filtered, so there is no new vs old index map for translation
      auto mcCollision = mctrack.template mcCollision_as<MyEventsMC>();
      trackMC(mcCollision.globalIndex(), mctrack.pdgCode(), mctrack.statusCode(), mctrack.flags(),
              mothers, daughterRange,
              mctrack.weight(), mctrack.pt(), mctrack.eta(), mctrack.phi(), mctrack.e(),
              mctrack.vx(), mctrack.vy(), mctrack.vz(), mctrack.vt(), mcflags);

      for (unsigned int isig = 0; isig < fMCSignals.size(); isig++) {
        if (mcflags & (static_cast<uint16_t>(1) << isig)) {
          (reinterpret_cast<TH1I*>(fStatsList->At(2)))->Fill(static_cast<float>(isig));
        }
      }
      if (mcflags == 0) {
        (reinterpret_cast<TH1I*>(fStatsList->At(2)))->Fill(static_cast<float>(fMCSignals.size()));
      }
    }
  }

  void processSkimming(MyEvents const& collisions,
                       MyBarrelTracks const& tracksBarrel, aod::TrackAssoc const& trackAssocs, 
                       MyEventsMC const& mcCollisions, aod::McParticles const& mcParticles)
  {
    fullSkimming(collisions, tracksBarrel, trackAssocs, mcCollisions, mcParticles);
  }

  PROCESS_SWITCH(Alice3DQTableMaker, processSkimming, "Build DQ skimmed data model for ALICE3", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Alice3DQTableMaker>(cfgc)
  };
}