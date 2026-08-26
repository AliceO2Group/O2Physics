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
/// \author Ionut Cristian Arsene <iarsene@cern.ch>, Oslo
/// \author Alexander Tiekoetter <alexander.tiekoetter@cern.ch>, Muenster
/// \brief Skimming Task for DQ Table Maker
/// \file alice3DqTableMaker.cxx

#include "PWGDQ/Core/AnalysisCompositeCut.h"
#include "PWGDQ/Core/AnalysisCut.h"
#include "PWGDQ/Core/CutsLibrary.h"
#include "PWGDQ/Core/HistogramManager.h"
#include "PWGDQ/Core/HistogramsLibrary.h"
#include "PWGDQ/Core/MCSignal.h"
#include "PWGDQ/Core/MCSignalLibrary.h"
#include "PWGDQ/Core/VarManager.h"

#include "ALICE3/DataModel/OTFRICH.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/ReducedTablesAlice3.h"
#include "ALICE3/DataModel/collisionAlice3.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <THashList.h>
#include <TList.h>
#include <TObjArray.h>
#include <TString.h>

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using MyBarrelTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                 aod::TracksCov, aod::TracksAlice3, aod::TracksExtraA3,
                                 aod::UpgradeTofs, aod::UpgradeRichs, aod::UpgradeRichSignals,
                                 aod::McTrackLabels>;

using MyEvents = soa::Join<aod::Collisions, aod::CollisionsAlice3, aod::McCollisionLabels>;
using MyEventsMC = aod::McCollisions;

constexpr static uint32_t GkEventFillMap = VarManager::ObjTypes::Collision;
constexpr static uint32_t GkEventMcFillMap = VarManager::ObjTypes::CollisionMC;

constexpr static uint32_t GkTrackFillMapWithCov = VarManager::ObjTypes::Track | VarManager::ObjTypes::TrackExtra | VarManager::ObjTypes::TrackDCA | VarManager::ObjTypes::TrackSelection | VarManager::ObjTypes::TrackCov | VarManager::ObjTypes::TrackPID;

namespace dqefficiency_helpers
{
inline float* varValues() { return static_cast<float*>(VarManager::fgValues); }
inline TString* varNames() { return static_cast<TString*>(VarManager::fgVariableNames); }
inline TString* varUnits() { return static_cast<TString*>(VarManager::fgVariableUnits); }
} // namespace dqefficiency_helpers

struct Alice3DqTableMaker {

  Produces<ReA3MCEvents> eventMC;
  Produces<ReA3MCTracks> trackMC;

  Produces<ReA3Events> event;
  Produces<ReducedA3EventsVtxCov> eventVtxCov;
  Produces<ReducedA3MCEventLabels> eventMClabels;

  Produces<ReA3Tracks> trackBasic;
  Produces<ReducedA3TracksBarrel> trackBarrel;
  Produces<ReducedA3TracksBarrelCov> trackBarrelCov;
  Produces<ReducedA3TracksAssoc> trackBarrelAssoc;
  Produces<ReducedA3TracksBarrelLabels> trackBarrelLabels;

  Produces<ReducedA3PIDTOF> trackPIDTOF;
  Produces<ReducedA3PIDRich> trackPIDRich;
  Produces<ReducedA3PIDRichSignals> trackPIDRichSig;

  OutputObj<THashList> fOutputList{"output"};
  OutputObj<TList> fStatsList{"Statistics"}; //! skimming statistics

  HistogramManager* fHistMan = nullptr;

  // Event and track AnalysisCut configurables
  struct : ConfigurableGroup {
    Configurable<std::string> cfgEventCuts{"cfgEventCuts", "", "Event selection"};
    Configurable<std::string> cfgBarrelTrackCuts{"cfgBarrelTrackCuts", "", "barrel track cut"};
    Configurable<std::string> cfgEventCutsJSON{"cfgEventCutsJSON", "", "Additional event selection in JSON format"};
    Configurable<std::string> cfgBarrelTrackCutsJSON{"cfgBarrelTrackCutsJSON", "", "Additional list of barrel track cuts in JSON format"};
  } fConfigCuts;

  // MC signals to be skimmed
  Configurable<std::string> cfgMCsignals{"cfgMCsignals", "", "Comma separated list of MC signals"};
  Configurable<std::string> cfgMCsignalsJSON{"cfgMCsignalsJSON", "", "Additional list of MC signals via JSON"};

  // Steer QA output
  struct : ConfigurableGroup {
    Configurable<bool> cfgQA{"cfgQA", false, "If true, fill QA histograms"};
    Configurable<bool> cfgDetailedQA{"cfgDetailedQA", false, "If true, include more QA histograms (BeforeCuts classes)"};
    Configurable<std::string> cfgAddEventHistogram{"cfgAddEventHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> cfgAddTrackHistogram{"cfgAddTrackHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> cfgAddMCTruthHistogram{"cfgAddMCTruthHistogram", "", "Comma separated list of histograms"};
    Configurable<std::string> cfgAddJSONHistograms{"cfgAddJSONHistograms", "", "Histograms in JSON format"};
  } fConfigHistOutput;

  AnalysisCompositeCut* fEventCut = nullptr;     //! Event selection cut
  std::vector<AnalysisCompositeCut*> fTrackCuts; //! Barrel track cuts

  bool fDoDetailedQA = false;

  std::vector<MCSignal*> fMCSignals;
  std::map<uint64_t, int> fLabelsMap;
  std::map<uint64_t, int> fLabelsMapReversed;
  std::map<uint64_t, uint16_t> fMCFlags;
  std::map<uint32_t, uint32_t> fCollIndexMap;  // key: old collision index, value: skimmed collision index
  std::map<uint32_t, uint32_t> fTrackIndexMap; // key: old track global index, value: new track global index

  void init(InitContext& context)
  {
    bool isProcessSkimmingEnabled = context.mOptions.get<bool>("processSkimming");

    if (!isProcessSkimmingEnabled)
      LOG(fatal) << "No process function was enabled ALICE 3 TableMaker";

    VarManager::SetDefaultVarNames(); // Important that this is called before defineCuts() !!!

    defineCuts();

    fHistMan = new HistogramManager("analysisHistos", "aa", VarManager::kNVars);
    fHistMan->SetUseDefaultVariableNames(true);
    fHistMan->SetDefaultVarNames(dqefficiency_helpers::varNames(), dqefficiency_helpers::varUnits());

    if (fConfigHistOutput.cfgQA && fConfigHistOutput.cfgDetailedQA) {
      fDoDetailedQA = true;
    }

    TString histClasses = "";

    if (fDoDetailedQA) {
      histClasses += "Event_BeforeCuts;";
    }

    if (fConfigHistOutput.cfgQA) {
      histClasses += "Event_AfterCuts;";
      histClasses += "Event_MCTruth;";
    }

    if (isProcessSkimmingEnabled) {
      if (fDoDetailedQA) {
        histClasses += "TrackBarrel_BeforeCuts;";
      }

      if (fConfigHistOutput.cfgQA) {
        for (const auto& cut : fTrackCuts) {
          histClasses += Form("TrackBarrel_%s;", cut->GetName());
        }
      }
    }

    TString configNamesStr = cfgMCsignals.value;
    std::unique_ptr<TObjArray> objArray(configNamesStr.Tokenize(","));

    if (objArray->GetEntries() > 0) {
      for (int isig = 0; isig < objArray->GetEntries(); ++isig) {
        MCSignal* sig = o2::aod::dqmcsignals::GetMCSignal(objArray->At(isig)->GetName());
        if (sig) {
          fMCSignals.push_back(sig);
        }
      }
    }

    TString addMCSignalsStr = cfgMCsignalsJSON.value;

    if (addMCSignalsStr != "") {
      std::vector<MCSignal*> addMCSignals = dqmcsignals::GetMCSignalsFromJSON(addMCSignalsStr.Data());

      for (const auto& mcIt : addMCSignals) {
        if (mcIt) {
          fMCSignals.push_back(mcIt);
        }
      }
    }

    for (const auto& mcIt : fMCSignals) {
      if (fConfigHistOutput.cfgQA) {
        histClasses += Form("MCTruth_%s;", mcIt->GetName());
      }
      if (fDoDetailedQA) {
        if (isProcessSkimmingEnabled) {
          for (const auto& cut : fTrackCuts) {
            histClasses += Form("TrackBarrel_%s_%s;", cut->GetName(), mcIt->GetName());
          }
        }
      }
    }

    defineHistograms(histClasses);

    TString addHistsStr = fConfigHistOutput.cfgAddJSONHistograms.value;
    if (fConfigHistOutput.cfgQA && addHistsStr != "") {
      dqhistograms::AddHistogramsFromJSON(fHistMan, addHistsStr.Data());
    }

    VarManager::SetUseVars(fHistMan->GetUsedVars()); // provide the list of required variables so that VarManager knows what to fill
    fOutputList.setObject(fHistMan->GetMainHistogramList());
  }

  void defineCuts()
  {
    fEventCut = new AnalysisCompositeCut(true);
    TString eventCutStr = fConfigCuts.cfgEventCuts.value;
    fEventCut->AddCut(dqcuts::GetAnalysisCut(eventCutStr.Data()));

    TString addEvCutsStr = fConfigCuts.cfgEventCutsJSON.value;
    if (addEvCutsStr != "") {
      std::vector<AnalysisCut*> addEvCuts = dqcuts::GetCutsFromJSON(addEvCutsStr.Data());
      for (const auto& cutIt : addEvCuts) {
        fEventCut->AddCut(cutIt);
      }
    }

    // Barrel track cuts
    TString cutNamesStr = fConfigCuts.cfgBarrelTrackCuts.value;
    if (!cutNamesStr.IsNull()) {
      std::unique_ptr<TObjArray> objArray(cutNamesStr.Tokenize(","));
      for (int icut = 0; icut < objArray->GetEntries(); ++icut) {
        fTrackCuts.push_back(dqcuts::GetCompositeCut(objArray->At(icut)->GetName()));
      }
    }
    // Additional Barrel track cuts via JSON
    TString addTrackCutsStr = fConfigCuts.cfgBarrelTrackCutsJSON.value;
    if (addTrackCutsStr != "") {
      std::vector<AnalysisCut*> addTrackCuts = dqcuts::GetCutsFromJSON(addTrackCutsStr.Data());
      for (const auto& t : addTrackCuts) {
        fTrackCuts.push_back(dynamic_cast<AnalysisCompositeCut*>(t));
      }
    }

    VarManager::SetUseVars(AnalysisCut::fgUsedVars); // provide the list of required variables so that VarManager knows what to fill
  }

  void defineHistograms(const TString& histClasses)
  {
    std::unique_ptr<TObjArray> objArray(histClasses.Tokenize(";"));
    for (int iclass = 0; iclass < objArray->GetEntries(); ++iclass) {
      TString classStr = objArray->At(iclass)->GetName();
      if (fConfigHistOutput.cfgQA) {
        fHistMan->AddHistClass(classStr.Data());
      }

      TString histEventName = fConfigHistOutput.cfgAddEventHistogram.value;
      if (classStr.Contains("Event")) {
        if (fConfigHistOutput.cfgQA && !classStr.Contains("MCTruth")) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "event", histEventName);
        } else {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "event", "generator");
        }
      }

      TString histTrackName = fConfigHistOutput.cfgAddTrackHistogram.value;
      if (classStr.Contains("Track")) {
        if (fConfigHistOutput.cfgQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "track", histTrackName);
        }
      }

      TString histMCTruthName = fConfigHistOutput.cfgAddMCTruthHistogram.value;
      if (classStr.Contains("MCTruth") && !classStr.Contains("Event")) {
        if (fConfigHistOutput.cfgQA) {
          dqhistograms::DefineHistograms(fHistMan, objArray->At(iclass)->GetName(), "mctruth_track", histMCTruthName);
        }
      }
    }

    // create statistics histograms (event, tracks, MCsignals)
    fStatsList.setObject(new TList());
    fStatsList->SetOwner(true);
    std::vector<TString> eventLabels{"Collisions before filtering", "Before cuts", "After cuts"};
    TH2I* histEvents = new TH2I("EventStats", "Event statistics", eventLabels.size(), -0.5, eventLabels.size() - 0.5, o2::aod::evsel::kNsel + 1, -0.5, (float)o2::aod::evsel::kNsel + 0.5);
    int ibX = 1;
    for (auto label = eventLabels.begin(); label != eventLabels.end(); label++, ibX++) {
      histEvents->GetXaxis()->SetBinLabel(ibX, (*label).Data());
    }
    for (int ibY = 1; ibY <= o2::aod::evsel::kNsel; ibY++) {
      histEvents->GetYaxis()->SetBinLabel(ibY, o2::aod::evsel::selectionLabels[ibY - 1]);
    }
    histEvents->GetYaxis()->SetBinLabel(o2::aod::evsel::kNsel + 1, "Total");
    fStatsList->Add(histEvents);

    // Track statistics: one bin for each track selection and 5 bins for V0 tags (gamma, K0s, Lambda, anti-Lambda, Omega)
    TH1I* histTracks = new TH1I("TrackStats", "Track statistics", fTrackCuts.size() + 5.0, -0.5, fTrackCuts.size() - 0.5 + 5.0);
    ibX = 1;
    for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, ibX++) {
      histTracks->GetXaxis()->SetBinLabel(ibX, (*cut)->GetName());
    }
    constexpr int NV0Tags = 5;
    const std::array<std::string, NV0Tags> v0TagNames = {"Photon conversion", "K^{0}_{s}", "#Lambda", "#bar{#Lambda}", "#Omega"};
    for (int ibY = 0; ibY < NV0Tags; ibY++) {
      histTracks->GetXaxis()->SetBinLabel(fTrackCuts.size() + 1 + ibY, v0TagNames[ibY].c_str());
    }
    fStatsList->Add(histTracks);

    TH1I* histMCsignals = new TH1I("MCsignals", "MC signals", fMCSignals.size() + 1, -0.5, fMCSignals.size() - 0.5 + 1.0);
    ibX = 1;
    for (auto signal = fMCSignals.begin(); signal != fMCSignals.end(); signal++, ibX++) {
      histMCsignals->GetXaxis()->SetBinLabel(ibX, (*signal)->GetName());
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

    for (const auto& mcCollision : mcCollisions) {
      VarManager::FillEventAlice3<GkEventMcFillMap>(mcCollision);

      fHistMan->FillHistClass("Event_MCTruth", dqefficiency_helpers::varValues());

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

    uint16_t mcflags{0}; // flags which will hold the decisions for each MC signal
    int trackCounter = 0;

    for (const auto& mctrack : mcTracks) {
      // check all the requested MC signals and fill the decision bit map
      mcflags = 0;
      int i = 0;
      for (const auto& sig : fMCSignals) {
        bool checked = false;
        if constexpr (soa::is_soa_filtered_v<aod::McParticles>) {
          auto mcTrackRaw = mcTracks.rawIteratorAt(mctrack.globalIndex());
          checked = sig->CheckSignal(true, mcTrackRaw);
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
      const auto mcTrackIndex = mctrack.globalIndex();
      const bool inserted = fLabelsMap.try_emplace(mcTrackIndex, trackCounter).second;

      if (inserted) {
        fLabelsMapReversed.try_emplace(trackCounter, mcTrackIndex);
        fMCFlags.try_emplace(mcTrackIndex, mcflags);
        ++trackCounter;

        // fill histograms for each of the signals, if found
        if (fConfigHistOutput.cfgQA) {
          VarManager::FillTrackMC(mcTracks, mctrack);
          auto mcCollision = mctrack.template mcCollision_as<MyEventsMC>();
          VarManager::FillEvent<GkEventMcFillMap>(mcCollision);

          int j = 0;
          for (auto signal = fMCSignals.begin(); signal != fMCSignals.end(); ++signal, ++j) {
            if ((mcflags & (static_cast<uint16_t>(1) << j)) != 0u) {
              fHistMan->FillHistClass(Form("MCTruth_%s", (*signal)->GetName()), dqefficiency_helpers::varValues());
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

      (dynamic_cast<TH2I*>(fStatsList->At(0)))->Fill(1.0, static_cast<float>(o2::aod::evsel::kNsel));

      VarManager::ResetValues(0, VarManager::kNEventWiseVariables);
      VarManager::FillEventAlice3<GkEventFillMap>(collision); // extract event information and place it in the fValues array

      if (collision.has_mcCollision()) {
        auto mcCollision = collision.template mcCollision_as<MyEventsMC>();
        VarManager::FillEventAlice3<GkEventMcFillMap>(mcCollision);
      }

      if (fDoDetailedQA) {
        fHistMan->FillHistClass("Event_BeforeCuts", dqefficiency_helpers::varValues());
      }

      // Apply the user specified event selection
      if (!fEventCut->IsSelected(dqefficiency_helpers::varValues())) {
        continue;
      }

      (dynamic_cast<TH2I*>(fStatsList->At(0)))->Fill(3.0, static_cast<float>(o2::aod::evsel::kNsel));

      // Fill historams after event cuts
      fHistMan->FillHistClass("Event_AfterCuts", VarManager::fgValues);

      event(collision.posX(), collision.posY(), collision.posZ(), collision.numContrib(),
            collision.collisionTime(), collision.collisionTimeRes(), collision.multDensity());

      eventVtxCov(collision.covXX(), collision.covXY(), collision.covXZ(), collision.covYY(), collision.covYZ(), collision.covZZ(), collision.chi2());
      eventMClabels(collision.mcCollisionId(), collision.mcMask());

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

    uint64_t trackFilteringTag{0};
    uint32_t trackTempFilterMap{0};
    uint16_t mcflags{0};
    int trackCounter = fLabelsMap.size();

    for (const auto& assoc : assocs) {

      auto track = assoc.template track_as<MyBarrelTracks>();

      if (!fCollIndexMap.contains(track.collisionId())) {
        continue;
      }

      trackFilteringTag = static_cast<uint64_t>(0);
      trackTempFilterMap = static_cast<uint32_t>(0);

      // Compute track quantities and fill histograms
      VarManager::FillTrackAlice3<GkTrackFillMapWithCov>(track);

      if (fDoDetailedQA) {
        fHistMan->FillHistClass("TrackBarrel_BeforeCuts", dqefficiency_helpers::varValues());
      }

      int n = 0;
      for (auto cut = fTrackCuts.begin(); cut != fTrackCuts.end(); cut++, n++) {
        if ((*cut)->IsSelected(dqefficiency_helpers::varValues())) {
          trackTempFilterMap |= (static_cast<uint32_t>(1) << n);
          if (fConfigHistOutput.cfgQA) {
            fHistMan->FillHistClass(Form("TrackBarrel_%s", (*cut)->GetName()), dqefficiency_helpers::varValues());
          }
          (dynamic_cast<TH1I*>(fStatsList->At(1)))->Fill(static_cast<float>(n));
        }
      }
      if (trackTempFilterMap == 0u) {
        continue;
      }

      // If this track is already present in the index map, it means it was already skimmed,
      // so we just store the association and we skip the track
      if (fTrackIndexMap.contains(track.globalIndex())) {
        trackBarrelAssoc(fCollIndexMap[collision.globalIndex()], fTrackIndexMap[track.globalIndex()]);
        continue;
      }

      trackFilteringTag |= (static_cast<uint64_t>(trackTempFilterMap) << VarManager::kBarrelUserCutsBits); // BIT13-...:  user track filters

      // NOTE: The collision ID that is written in the table is the one originally assigned in the AOD.
      //       However, in data analysis one should loop over associations, so this one should not be used.
      //      In the case of Run2-like analysis, there will be no associations, so this ID will be the one originally assigned in the AO2Ds (updated for the skims)
      uint32_t reducedEventIdx = fCollIndexMap[track.collisionId()];

      trackBasic(reducedEventIdx, trackFilteringTag, track.pt(), track.eta(), track.phi(), track.sign(), 0);

      trackBarrel(track.x(), track.alpha(), track.y(), track.z(), track.snp(), track.tgl(), track.signed1Pt(),
                  track.flags(), track.itsClusterMap(), track.itsChi2NCl(),
                  track.isReconstructed(), track.nSiliconHits(), track.nTPCHits(), track.length(), track.dcaXY(),
                  track.dcaZ());

      if constexpr (static_cast<bool>(GkTrackFillMapWithCov & VarManager::ObjTypes::TrackCov)) {
        trackBarrelCov(track.cYY(), track.cZY(), track.cZZ(), track.cSnpY(), track.cSnpZ(),
                       track.cSnpSnp(), track.cTglY(), track.cTglZ(), track.cTglSnp(), track.cTglTgl(),
                       track.c1PtY(), track.c1PtZ(), track.c1PtSnp(), track.c1PtTgl(), track.c1Pt21Pt2());
      }

      if constexpr (static_cast<bool>(GkTrackFillMapWithCov & VarManager::ObjTypes::TrackPID)) {

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
      }

      fTrackIndexMap[track.globalIndex()] = trackBasic.lastIndex();

      // Check whether the MCParticle corresponding to this reconstructed track was already selected for skimming
      // If not, add it to the skimming map
      if (!track.has_mcParticle()) {
        trackBarrelLabels(-1, 0, 0); // this is the case when there is no matched MCParticle
      } else {
        auto mctrack = track.template mcParticle_as<aod::McParticles>();
        VarManager::FillTrackMC(mcTracks, mctrack);
        VarManager::FillResolutions(mctrack, track);

        mcflags = 0;
        int i = 0; // runs over the MC signals
        int j = 0; // runs over the track cuts
        // check all the specified signals and fill histograms for MC truth matched tracks
        for (const auto& sig : fMCSignals) {
          if (sig->CheckSignal(true, mctrack)) {
            mcflags |= (static_cast<uint16_t>(1) << i);
            // If detailed QA is on, fill histograms for each MC signal and track cut combination
            if (fDoDetailedQA) {
              j = 0;
              for (const auto& cut : fTrackCuts) {
                if ((trackTempFilterMap & (uint8_t(1) << j)) != 0u) {
                  fHistMan->FillHistClass(Form("TrackBarrel_%s_%s", cut->GetName(), sig->GetName()), dqefficiency_helpers::varValues()); // fill the reconstructed truth
                }
                j++;
              }
            }
          }
          i++;
        }

        // if the MC truth particle corresponding to this reconstructed track is not already written,
        //   add it to the skimmed stack
        if (!(fLabelsMap.contains(mctrack.globalIndex()))) {
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

    skimCollisions(collisions);

    if (fCollIndexMap.empty())
      return;

    skimMCParticles(mcParticles, mcCollisions);

    if constexpr (static_cast<bool>(GkTrackFillMapWithCov)) {
      fTrackIndexMap.clear();
      trackBasic.reserve(tracksBarrel.size());
      trackBarrel.reserve(tracksBarrel.size());
      trackBarrelCov.reserve(tracksBarrel.size());
      trackPIDTOF.reserve(tracksBarrel.size());
      trackPIDRich.reserve(tracksBarrel.size());
      trackPIDRichSig.reserve(tracksBarrel.size());
      trackBarrelAssoc.reserve(tracksBarrel.size());
      trackBarrelLabels.reserve(tracksBarrel.size());
    }

    if (!fCollIndexMap.empty()) {

      for (auto const& [origIdx, skimIdx] : fCollIndexMap) {
        auto collision = collisions.rawIteratorAt(origIdx);

        if constexpr (static_cast<bool>(GkTrackFillMapWithCov)) {
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
        for (const auto& m : mctrack.mothersIds()) {
          if (m < mcParticles.size()) { // protect against bad mother indices
            if (fLabelsMap.contains(m)) {
              mothers.push_back(fLabelsMap.find(m)->second);
            }
          } else {
            LOG(warn) << "Mother label (" << m << ") exceeds the McParticles size (" << mcParticles.size() << ")";
            LOG(warn) << " Check the MC generator";
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
            if (fLabelsMap.contains(d)) {
              daughters.push_back(fLabelsMap.find(d)->second);
            }
          } else {
            LOG(warn) << "Daughter label (" << d << ") exceeds the McParticles size (" << mcParticles.size() << ")";
            LOG(warn) << " Check the MC generator";
          }
        }
      }
      constexpr int NDaughters = 2;
      std::array<int, NDaughters> daughterRange = {-1, -1};
      if (!daughters.empty()) {
        daughterRange[0] = daughters[0];
        daughterRange[1] = daughters[daughters.size() - 1];
      }

      // NOTE: Here we assume that MC collisions are not filtered, so there is no new vs old index map for translation
      auto mcCollision = mctrack.template mcCollision_as<MyEventsMC>();
      trackMC(mcCollision.globalIndex(), mctrack.pdgCode(), mctrack.statusCode(), mctrack.flags(),
              mothers, daughterRange.data(),
              mctrack.weight(), mctrack.pt(), mctrack.eta(), mctrack.phi(), mctrack.e(),
              mctrack.vx(), mctrack.vy(), mctrack.vz(), mctrack.vt(), mcflags);

      for (unsigned int isig = 0; isig < fMCSignals.size(); isig++) {
        if ((mcflags & (static_cast<uint16_t>(1) << isig)) != 0u) {
          (dynamic_cast<TH1I*>(fStatsList->At(2)))->Fill(static_cast<float>(isig));
        }
      }
      if (mcflags == 0) {
        (dynamic_cast<TH1I*>(fStatsList->At(2)))->Fill(static_cast<float>(fMCSignals.size()));
      }
    }
  }

  void processSkimming(MyEvents const& collisions,
                       MyBarrelTracks const& tracksBarrel, aod::TrackAssoc const& trackAssocs,
                       MyEventsMC const& mcCollisions, aod::McParticles const& mcParticles)
  {
    fullSkimming(collisions, tracksBarrel, trackAssocs, mcCollisions, mcParticles);
  }

  PROCESS_SWITCH(Alice3DqTableMaker, processSkimming, "Build DQ skimmed data model for ALICE3", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Alice3DqTableMaker>(cfgc)};
}
