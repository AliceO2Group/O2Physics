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

#include <cmath>

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/Configurable.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"
#include "Index.h"
#include "TDatabasePDG.h"

#include "bestCollisionTable.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;

AxisSpec ZAxis = {301, -30.1, 30.1};
AxisSpec DeltaZAxis = {61, -6.1, 6.1};
AxisSpec DCAAxis = {601, -3.01, 3.01};
AxisSpec EtaAxis = {22, -2.2, 2.2};
// AxisSpec MultAxis = {301, -0.5, 300.5};
AxisSpec PhiAxis = {629, 0, 2 * M_PI};
AxisSpec PtAxis = {2401, -0.005, 24.005};
AxisSpec PtAxis_wide = {1041, -0.05, 104.05};

static constexpr TrackSelectionFlags::flagtype trackSelectionITS =
  TrackSelectionFlags::kITSNCls | TrackSelectionFlags::kITSChi2NDF |
  TrackSelectionFlags::kITSHits;

static constexpr TrackSelectionFlags::flagtype trackSelectionTPC =
  TrackSelectionFlags::kTPCNCls |
  TrackSelectionFlags::kTPCCrossedRowsOverNCls |
  TrackSelectionFlags::kTPCChi2NDF;

static constexpr TrackSelectionFlags::flagtype trackSelectionDCA =
  TrackSelectionFlags::kDCAz | TrackSelectionFlags::kDCAxy;

using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;

struct MultiplicityCounter {
  SliceCache cache;
  Preslice<aod::Tracks> perCol = aod::track::collisionId;
  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;

  Service<O2DatabasePDG> pdg;

  Configurable<float> estimatorEta{"estimatorEta", 1.0, "eta range for INEL>0 sample definition"};
  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  Configurable<bool> fillResponse{"fillResponse", false, "Fill response matrix"};
  ConfigurableAxis multBinning{"multBinning", {301, -0.5, 300.5}, ""};
  ConfigurableAxis centBinning{"centBinning", {VARIABLE_WIDTH, 0, 10, 20, 30, 40, 50, 60, 70, 80, 100}, ""};

  HistogramRegistry registry{
    "registry",
    {
      {"Events/BCSelection", ";status;count", {HistType::kTH1F, {{3, 0.5, 3.5}}}} //
    }                                                                             //
  };

  std::vector<int> usedTracksIds;
  std::vector<int> usedTracksIdsDF;

  void init(InitContext&)
  {
    AxisSpec MultAxis = {multBinning, "N_{trk}"};
    AxisSpec CentAxis = {centBinning, "centrality"};

    auto hstat = registry.get<TH1>(HIST("Events/BCSelection"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "Good BCs");
    x->SetBinLabel(2, "BCs with collisions");
    x->SetBinLabel(3, "BCs with pile-up/splitting");

    if (doprocessEventStat) {
      registry.add({"Events/Control/Chi2", " ; #chi^2", {HistType::kTH1F, {{101, -0.1, 10.1}}}});
      registry.add({"Events/Control/TimeResolution", " ; t (ms)", {HistType::kTH1F, {{1001, -0.1, 100.1}}}});
    }
    if (doprocessEventStatCentralityFT0C || doprocessEventStatCentralityFT0M) {
      registry.add({"Events/Centrality/Control/Chi2", " ; #chi^2; centrality", {HistType::kTH2F, {{101, -0.1, 10.1}, CentAxis}}});
      registry.add({"Events/Centrality/Control/TimeResolution", " ; t (ms); centrality", {HistType::kTH2F, {{1001, -0.1, 100.1}, CentAxis}}});
    }

    if (doprocessCounting) {
      registry.add({"Events/Selection", ";status;events", {HistType::kTH1F, {{4, 0.5, 4.5}}}});
      auto hstat = registry.get<TH1>(HIST("Events/Selection"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All");
      x->SetBinLabel(2, "Selected");
      x->SetBinLabel(3, "Selected INEL>0");
      x->SetBinLabel(4, "Rejected");

      registry.add({"Events/NtrkZvtx", "; N_{trk}; Z_{vtx} (cm); events", {HistType::kTH2F, {MultAxis, ZAxis}}});
      registry.add({"Tracks/EtaZvtx", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"Tracks/EtaZvtx_gt0", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"Tracks/PhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      registry.add({"Tracks/Control/PtEta", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, EtaAxis}}});
      registry.add({"Tracks/Control/DCAXYPt", " ; p_{T} (GeV/c) ; DCA_{XY} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
      registry.add({"Tracks/Control/DCAZPt", " ; p_{T} (GeV/c) ; DCA_{Z} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
      registry.add({"Tracks/Control/ReassignedDCAXYPt", " ; p_{T} (GeV/c) ; DCA_{XY} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
      registry.add({"Tracks/Control/ReassignedDCAZPt", " ; p_{T} (GeV/c) ; DCA_{Z} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
      registry.add({"Tracks/Control/ExtraDCAXYPt", " ; p_{T} (GeV/c) ; DCA_{XY} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
      registry.add({"Tracks/Control/ExtraDCAZPt", " ; p_{T} (GeV/c) ; DCA_{Z} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
      registry.add({"Tracks/Control/ExtraTracksEtaZvtx", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"Tracks/Control/ExtraTracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      registry.add({"Tracks/Control/ReassignedTracksEtaZvtx", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"Tracks/Control/ReassignedTracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      registry.add({"Tracks/Control/ReassignedVertexCorr", "; Z_{vtx}^{orig} (cm); Z_{vtx}^{re} (cm)", {HistType::kTH2F, {ZAxis, ZAxis}}});
    }

    if (doprocessCountingCentralityFT0C || doprocessCountingCentralityFT0M) {
      registry.add({"Events/Centrality/Selection", ";status;centrality;events", {HistType::kTH2F, {{3, 0.5, 3.5}, CentAxis}}});
      auto hstat = registry.get<TH2>(HIST("Events/Centrality/Selection"));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "All");
      x->SetBinLabel(2, "Selected");
      x->SetBinLabel(3, "Rejected");

      registry.add({"Events/Centrality/NtrkZvtx", "; N_{trk}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtx", "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/PhiEta", "; #varphi; #eta; centrality", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/PtEta", " ; p_{T} (GeV/c); #eta; centrality", {HistType::kTHnSparseF, {PtAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/DCAXYPt", " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/DCAZPt", " ; p_{T} (GeV/c) ; DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ReassignedDCAXYPt", " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ReassignedDCAZPt", " ; p_{T} (GeV/c) ; DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ExtraDCAXYPt", " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ExtraDCAZPt", " ; p_{T} (GeV/c) ; DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ExtraTracksEtaZvtx", "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ExtraTracksPhiEta", "; #varphi; #eta; centrality", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ReassignedTracksEtaZvtx", "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ReassignedTracksPhiEta", "; #varphi; #eta; centrality", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/ReassignedVertexCorr", "; Z_{vtx}^{orig} (cm); Z_{vtx}^{re} (cm); centrality", {HistType::kTHnSparseF, {ZAxis, ZAxis, CentAxis}}});
    }

    if (doprocessGen) {
      registry.add({"Events/NtrkZvtxGen", "; N_{trk}; Z_{vtx} (cm); events", {HistType::kTH2F, {MultAxis, ZAxis}}});
      registry.add({"Events/NtrkZvtxGen_t", "; N_{part}; Z_{vtx} (cm); events", {HistType::kTH2F, {MultAxis, ZAxis}}});
      registry.add({"Tracks/EtaZvtxGen", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"Tracks/EtaZvtxGen_t", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"Tracks/EtaZvtxGen_gt0", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"Tracks/EtaZvtxGen_gt0t", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      registry.add({"Tracks/Control/PtEtaGen", " ; p_{T} (GeV/c) ; #eta", {HistType::kTH2F, {PtAxis, EtaAxis}}});

      registry.add({"Tracks/PhiEtaGen", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      registry.add({"Tracks/Control/PhiEtaGenDuplicates", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      registry.add({"Tracks/Control/PhiEtaDuplicates", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      registry.add({"Events/Efficiency", "; status; events", {HistType::kTH1F, {{5, 0.5, 5.5}}}});
      registry.add({"Events/NotFoundEventZvtx", " ; Z_{vtx} (cm)", {HistType::kTH1F, {ZAxis}}});

      if (fillResponse) {
        registry.add({"Events/Response", " ; N_{rec}; N_{gen}; Z_{vtx} (cm)", {HistType::kTHnSparseF, {MultAxis, MultAxis, ZAxis}}});
        registry.add({"Events/EfficiencyMult", " ; N_{gen}; Z_{vtx} (cm)", {HistType::kTH2F, {MultAxis, ZAxis}}});
        registry.add({"Events/SplitMult", " ; N_{gen} ; Z_{vtx} (cm)", {HistType::kTH2F, {MultAxis, ZAxis}}});
      }

      auto heff = registry.get<TH1>(HIST("Events/Efficiency"));
      auto* x = heff->GetXaxis();
      x->SetBinLabel(1, "Generated");
      x->SetBinLabel(2, "Generated INEL>0");
      x->SetBinLabel(3, "Reconstructed");
      x->SetBinLabel(4, "Selected");
      x->SetBinLabel(5, "Selected INEL>0");
    }

    if (doprocessGenFT0C || doprocessGenFT0M) {
      registry.add({"Events/Centrality/NtrkZvtxGen", "; N_{trk}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});
      registry.add({"Events/Centrality/NtrkZvtxGen_t", "; N_{part}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen", "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen_t", "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen_gt0", "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/EtaZvtxGen_gt0t", "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/PtEtaGen", " ; p_{T} (GeV/c) ; #eta; centrality", {HistType::kTHnSparseF, {PtAxis, EtaAxis, CentAxis}}});

      registry.add({"Tracks/Centrality/PhiEtaGen", "; #varphi; #eta; tracks", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/PhiEtaGenDuplicates", "; #varphi; #eta; tracks", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
      registry.add({"Tracks/Centrality/Control/PhiEtaDuplicates", "; #varphi; #eta; tracks", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
      registry.add({"Events/Centrality/Efficiency", "; status; centrality; events", {HistType::kTH2F, {{5, 0.5, 5.5}, CentAxis}}});
      registry.add({"Events/Centrality/NotFoundEventZvtx", " ; Z_{vtx} (cm); centrality; events", {HistType::kTH2F, {ZAxis, CentAxis}}});

      if (fillResponse) {
        registry.add({"Events/Centrality/Response", " ; N_{rec}; N_{gen}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, MultAxis, ZAxis, CentAxis}}});
        registry.add({"Events/Centrality/EfficiencyMult", " ; N_{gen}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});
        registry.add({"Events/Centrality/SplitMult", " ; N_{gen} ; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});
      }

      auto heff = registry.get<TH2>(HIST("Events/Centrality/Efficiency"));
      auto* x = heff->GetXaxis();
      x->SetBinLabel(1, "Generated");
      x->SetBinLabel(2, "Generated INEL>0");
      x->SetBinLabel(3, "Reconstructed");
      x->SetBinLabel(4, "Selected");
      x->SetBinLabel(5, "Selected INEL>0");
    }

    if (doprocessTrackEfficiency) {
      registry.add({"Tracks/Control/PtGen", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
      registry.add({"Tracks/Control/PtEfficiency", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
      registry.add({"Tracks/Control/PtEfficiencySecondaries", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
    }
    if (doprocessTrackEfficiencyIndexed) {
      registry.add({"Tracks/Control/PtGenI", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
      registry.add({"Tracks/Control/PtEfficiencyI", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
      registry.add({"Tracks/Control/PtEfficiencyINoEtaCut", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
      registry.add({"Tracks/Control/PtEfficiencyISecondaries", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
      registry.add({"Tracks/Control/Mask", " ; bit", {HistType::kTH1F, {{17, -0.5, 16.5}}}});
      registry.add({"Tracks/Control/ITSClusters", " ; layer", {HistType::kTH1F, {{8, 0.5, 8.5}}}});
    }
  }

  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  template <typename C>
  void processEventStatGeneral(FullBCs const& bcs, C const& collisions)
  {
    constexpr bool hasCentrality = C::template contains<aod::CentFT0Cs>() || C::template contains<aod::CentFT0Ms>();
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto& bc : bcs) {
      if (!useEvSel || (bc.selection()[evsel::kIsBBT0A] &
                        bc.selection()[evsel::kIsBBT0C]) != 0) {
        registry.fill(HIST("Events/BCSelection"), 1.);
        cols.clear();
        for (auto& collision : collisions) {
          if (collision.has_foundBC()) {
            if (collision.foundBCId() == bc.globalIndex()) {
              cols.emplace_back(collision);
            }
          } else if (collision.bcId() == bc.globalIndex()) {
            cols.emplace_back(collision);
          }
        }
        LOGP(debug, "BC {} has {} collisions", bc.globalBC(), cols.size());
        if (!cols.empty()) {
          registry.fill(HIST("Events/BCSelection"), 2.);
          if (cols.size() > 1) {
            registry.fill(HIST("Events/BCSelection"), 3.);
          }
        }
        for (auto& col : cols) {
          float c = -1;
          if constexpr (hasCentrality) {
            if constexpr (C::template contains<aod::CentFT0Cs>()) {
              c = col.centFT0C();
            } else if constexpr (C::template contains<aod::CentFT0Ms>()) {
              c = col.centFT0M();
            }
            registry.fill(HIST("Events/Centrality/Control/Chi2"), col.chi2(), c);
            registry.fill(HIST("Events/Centrality/Control/TimeResolution"), col.collisionTimeRes(), c);
          } else {
            registry.fill(HIST("Events/Control/Chi2"), col.chi2());
            registry.fill(HIST("Events/Control/TimeResolution"), col.collisionTimeRes());
          }
        }
      }
    }
  }

  void processEventStat(
    FullBCs const& bcs,
    soa::Join<aod::Collisions, aod::EvSels> const& collisions)
  {
    processEventStatGeneral(bcs, collisions);
  }

  PROCESS_SWITCH(MultiplicityCounter, processEventStat, "Collect event sample stats", false);

  void processEventStatCentralityFT0C(
    FullBCs const& bcs,
    soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs> const& collisions)
  {
    processEventStatGeneral(bcs, collisions);
  }

  PROCESS_SWITCH(MultiplicityCounter, processEventStatCentralityFT0C, "Collect event sample stats (FT0C binned)", false);

  void processEventStatCentralityFT0M(
    FullBCs const& bcs,
    soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms> const& collisions)
  {
    processEventStatGeneral(bcs, collisions);
  }

  PROCESS_SWITCH(MultiplicityCounter, processEventStatCentralityFT0M, "Collect event sample stats (FT0M binned)", false);

  // clean up used Ids each dataframe (default process is always executed first)
  void process(aod::Collisions const&)
  {
    usedTracksIdsDF.clear();
  }

  expressions::Filter trackSelectionProper = ((aod::track::trackCutFlag & trackSelectionITS) == trackSelectionITS) &&
                                             ifnode((aod::track::detectorMap & (uint8_t)o2::aod::track::TPC) == (uint8_t)o2::aod::track::TPC,
                                                    (aod::track::trackCutFlag & trackSelectionTPC) == trackSelectionTPC,
                                                    true) &&
                                             ((aod::track::trackCutFlag & trackSelectionDCA) == trackSelectionDCA);

  using ExTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
  using FiTracks = soa::Filtered<ExTracks>;

  expressions::Filter atrackFilter = (aod::track::bestCollisionId >= 0) &&
                                     (nabs(aod::track::bestDCAZ) <= 2.f) &&
                                     (nabs(aod::track::bestDCAXY) <= ((0.0105f + 0.0350f / npow(aod::track::pts, 1.1f))));
  using ExCols = soa::Join<aod::Collisions, aod::EvSels>;

  template <typename C>
  void processCountingGeneral(
    typename C::iterator const& collision,
    FiTracks const& tracks,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    float c = -1;
    constexpr bool hasCentrality = C::template contains<aod::CentFT0Cs>() || C::template contains<aod::CentFT0Ms>();
    if constexpr (hasCentrality) {
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        c = collision.centFT0C();
      } else if (C::template contains<aod::CentFT0Ms>()) {
        c = collision.centFT0M();
      }
      registry.fill(HIST("Events/Centrality/Selection"), 1., c);
    } else {
      registry.fill(HIST("Events/Selection"), 1.);
    }

    if (!useEvSel || collision.sel8()) {
      if constexpr (hasCentrality) {
        registry.fill(HIST("Events/Centrality/Selection"), 2., c);
      } else {
        registry.fill(HIST("Events/Selection"), 2.);
      }
      auto z = collision.posZ();
      usedTracksIds.clear();

      auto Ntrks = 0;
      for (auto& track : atracks) {
        auto otrack = track.track_as<FiTracks>();
        usedTracksIds.emplace_back(track.trackId());
        if (std::abs(otrack.eta()) < estimatorEta) {
          ++Ntrks;
        }
        if constexpr (hasCentrality) {
          registry.fill(HIST("Tracks/Centrality/EtaZvtx"), otrack.eta(), z, c);
          registry.fill(HIST("Tracks/Centrality/PhiEta"), otrack.phi(), otrack.eta(), c);
          registry.fill(HIST("Tracks/Centrality/Control/PtEta"), otrack.pt(), otrack.eta(), c);
          registry.fill(HIST("Tracks/Centrality/Control/DCAXYPt"), otrack.pt(), track.bestDCAXY(), c);
          registry.fill(HIST("Tracks/Centrality/Control/DCAZPt"), otrack.pt(), track.bestDCAZ(), c);
        } else {
          registry.fill(HIST("Tracks/EtaZvtx"), otrack.eta(), z);
          registry.fill(HIST("Tracks/PhiEta"), otrack.phi(), otrack.eta());
          registry.fill(HIST("Tracks/Control/PtEta"), otrack.pt(), otrack.eta());
          registry.fill(HIST("Tracks/Control/DCAXYPt"), otrack.pt(), track.bestDCAXY());
          registry.fill(HIST("Tracks/Control/DCAZPt"), otrack.pt(), track.bestDCAZ());
        }
        if (!otrack.has_collision()) {
          if constexpr (hasCentrality) {
            registry.fill(HIST("Tracks/Centrality/Control/ExtraTracksEtaZvtx"), otrack.eta(), z, c);
            registry.fill(HIST("Tracks/Centrality/Control/ExtraTracksPhiEta"), otrack.phi(), otrack.eta(), c);
            registry.fill(HIST("Tracks/Centrality/Control/ExtraDCAXYPt"), otrack.pt(), track.bestDCAXY(), c);
            registry.fill(HIST("Tracks/Centrality/Control/ExtraDCAZPt"), otrack.pt(), track.bestDCAZ(), c);
          } else {
            registry.fill(HIST("Tracks/Control/ExtraTracksEtaZvtx"), otrack.eta(), z);
            registry.fill(HIST("Tracks/Control/ExtraTracksPhiEta"), otrack.phi(), otrack.eta());
            registry.fill(HIST("Tracks/Control/ExtraDCAXYPt"), otrack.pt(), track.bestDCAXY());
            registry.fill(HIST("Tracks/Control/ExtraDCAZPt"), otrack.pt(), track.bestDCAZ());
          }
        } else if (otrack.collisionId() != track.bestCollisionId()) {
          usedTracksIdsDF.emplace_back(track.trackId());
          if constexpr (hasCentrality) {
            registry.fill(HIST("Tracks/Centrality/Control/ReassignedTracksEtaZvtx"), otrack.eta(), z, c);
            registry.fill(HIST("Tracks/Centrality/Control/ReassignedTracksPhiEta"), otrack.phi(), otrack.eta(), c);
            registry.fill(HIST("Tracks/Centrality/Control/ReassignedVertexCorr"), otrack.collision_as<C>().posZ(), z, c);
            registry.fill(HIST("Tracks/Centrality/Control/ReassignedDCAXYPt"), otrack.pt(), track.bestDCAXY(), c);
            registry.fill(HIST("Tracks/Centrality/Control/ReassignedDCAZPt"), otrack.pt(), track.bestDCAZ(), c);
          } else {
            registry.fill(HIST("Tracks/Control/ReassignedTracksEtaZvtx"), otrack.eta(), z);
            registry.fill(HIST("Tracks/Control/ReassignedTracksPhiEta"), otrack.phi(), otrack.eta());
            registry.fill(HIST("Tracks/Control/ReassignedVertexCorr"), otrack.collision_as<C>().posZ(), z);
            registry.fill(HIST("Tracks/Control/ReassignedDCAXYPt"), otrack.pt(), track.bestDCAXY());
            registry.fill(HIST("Tracks/Control/ReassignedDCAZPt"), otrack.pt(), track.bestDCAZ());
          }
        }
      }
      for (auto& track : tracks) {
        if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
          continue;
        }
        if (std::find(usedTracksIdsDF.begin(), usedTracksIdsDF.end(), track.globalIndex()) != usedTracksIdsDF.end()) {
          continue;
        }
        if (std::abs(track.eta()) < estimatorEta) {
          ++Ntrks;
        }
        if constexpr (hasCentrality) {
          registry.fill(HIST("Tracks/Centrality/EtaZvtx"), track.eta(), z, c);
          registry.fill(HIST("Tracks/Centrality/PhiEta"), track.phi(), track.eta(), c);
          registry.fill(HIST("Tracks/Centrality/Control/PtEta"), track.pt(), track.eta(), c);
          registry.fill(HIST("Tracks/Centrality/Control/DCAXYPt"), track.pt(), track.dcaXY(), c);
          registry.fill(HIST("Tracks/Centrality/Control/DCAZPt"), track.pt(), track.dcaZ(), c);
        } else {
          registry.fill(HIST("Tracks/EtaZvtx"), track.eta(), z);
          registry.fill(HIST("Tracks/PhiEta"), track.phi(), track.eta());
          registry.fill(HIST("Tracks/Control/PtEta"), track.pt(), track.eta());
          registry.fill(HIST("Tracks/Control/DCAXYPt"), track.pt(), track.dcaXY());
          registry.fill(HIST("Tracks/Control/DCAZPt"), track.pt(), track.dcaZ());
        }
      }
      if constexpr (hasCentrality) {
        registry.fill(HIST("Events/Centrality/NtrkZvtx"), Ntrks, z, c);
      } else {
        if (Ntrks > 0) {
          registry.fill(HIST("Events/Selection"), 3.);
          for (auto& track : atracks) {
            registry.fill(HIST("Tracks/EtaZvtx_gt0"), track.track_as<FiTracks>().eta(), z);
          }
          for (auto& track : tracks) {
            if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
              continue;
            }
            if (std::find(usedTracksIdsDF.begin(), usedTracksIdsDF.end(), track.globalIndex()) != usedTracksIdsDF.end()) {
              continue;
            }
            registry.fill(HIST("Tracks/EtaZvtx_gt0"), track.eta(), z);
          }
        }
        registry.fill(HIST("Events/NtrkZvtx"), Ntrks, z);
      }
    } else {
      if constexpr (hasCentrality) {
        registry.fill(HIST("Events/Centrality/Selection"), 3., c);
      } else {
        registry.fill(HIST("Events/Selection"), 4.);
      }
    }
  }

  void processCounting(
    ExCols::iterator const& collision,
    FiTracks const& tracks,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    processCountingGeneral<ExCols>(collision, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCounting, "Count tracks", false);

  using ExColsCentFT0C = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>;
  void processCountingCentralityFT0C(
    ExColsCentFT0C::iterator const& collision,
    FiTracks const& tracks,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    processCountingGeneral<ExColsCentFT0C>(collision, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCountingCentralityFT0C, "Count tracks in centrality bins", false);

  using ExColsCentFT0M = soa::Join<aod::Collisions, aod::CentFT0Ms, aod::EvSels>;
  void processCountingCentralityFT0M(
    ExColsCentFT0M::iterator const& collision,
    FiTracks const& tracks,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    processCountingGeneral<ExColsCentFT0M>(collision, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCountingCentralityFT0M, "Count tracks in centrality bins", false);

  using Particles = soa::Filtered<aod::McParticles>;
  using LabeledTracksEx = soa::Join<LabeledTracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
  using Particle = Particles::iterator;
  using ParticlesI = soa::Join<aod::McParticles, aod::ParticlesToTracks>;
  expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;
  Partition<Particles> mcSample = nabs(aod::mcparticle::eta) < estimatorEta;
  Partition<ParticlesI> primariesI = ((aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) &&
                                     (nabs(aod::mcparticle::eta) < estimatorEta);

  void processTrackEfficiencyIndexed(
    soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const&, ParticlesI const&,
    soa::Filtered<LabeledTracksEx> const& tracks)
  {
    for (auto& collision : collisions) {
      if (useEvSel && !collision.sel8()) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto mcCollision = collision.mcCollision();
      auto particlesI = primariesI->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      particlesI.bindExternalIndices(&tracks);

      for (auto& particle : particlesI) {
        auto charge = 0.;
        auto p = pdg->GetParticle(particle.pdgCode());
        if (p != nullptr) {
          charge = p->Charge();
        }
        if (std::abs(charge) < 3.) {
          continue;
        }
        registry.fill(HIST("Tracks/Control/PtGenI"), particle.pt());
        if (particle.has_tracks()) {
          auto counted = false;
          auto countedNoEtaCut = false;
          auto counter = 0;
          auto relatedTracks = particle.tracks_as<soa::Filtered<LabeledTracksEx>>();
          for (auto& track : relatedTracks) {
            if ((track.trackCutFlag() & TrackSelectionFlags::kDCAz) != TrackSelectionFlags::kDCAz) {
              continue;
            }
            ++counter;
            if (!countedNoEtaCut) {
              registry.fill(HIST("Tracks/Control/PtEfficiencyINoEtaCut"), particle.pt());
              countedNoEtaCut = true;
            }
            if (std::abs(track.eta()) < estimatorEta) {
              if (!counted) {
                registry.fill(HIST("Tracks/Control/PtEfficiencyI"), particle.pt());
                counted = true;
              }
            }
            if (counter > 1) {
              registry.fill(HIST("Tracks/Control/PtEfficiencyISecondaries"), particle.pt());
            }
          }
          if (counter > 1) {
            for (auto& track : relatedTracks) {
              if ((track.trackCutFlag() & TrackSelectionFlags::kDCAz) != TrackSelectionFlags::kDCAz) {
                continue;
              }
              for (auto layer = 0; layer < 7; ++layer) {
                if (track.itsClusterMap() & (uint8_t(1) << layer)) {
                  registry.fill(HIST("Tracks/Control/ITSClusters"), layer + 1);
                }
              }
              auto hasbit = false;
              for (auto bit = 0; bit < 16; ++bit) {
                if (track.mcMask() & (uint8_t(1) << bit)) {
                  registry.fill(HIST("Tracks/Control/Mask"), bit);
                  hasbit = true;
                }
              }
              if (!hasbit) {
                registry.fill(HIST("Tracks/Control/Mask"), 16);
              }
            }
          }
          if (relatedTracks.size() > 1) {
            registry.fill(HIST("Tracks/Control/PhiEtaGenDuplicates"), particle.phi(), particle.eta());
            for (auto& track : relatedTracks) {
              if ((track.trackCutFlag() & TrackSelectionFlags::kDCAz) != TrackSelectionFlags::kDCAz) {
                continue;
              }
              registry.fill(HIST("Tracks/Control/PhiEtaDuplicates"), track.phi(), track.eta());
            }
          }
        }
      }
    }
  }

  PROCESS_SWITCH(MultiplicityCounter, processTrackEfficiencyIndexed, "Calculate tracking efficiency vs pt (indexed)", false);

  Partition<soa::Filtered<LabeledTracksEx>> lsample = nabs(aod::track::eta) < estimatorEta;
  void processTrackEfficiency(
    soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions,
    aod::McCollisions const&, Particles const& mcParticles,
    soa::Filtered<LabeledTracksEx> const&,
    soa::SmallGroups<aod::ReassignedTracksCore> const& atracks)
  {
    for (auto& collision : collisions) {
      if (useEvSel && !collision.sel8()) {
        continue;
      }
      if (!collision.has_mcCollision()) {
        continue;
      }
      auto mcCollision = collision.mcCollision();
      auto particles = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      auto tracks = lsample->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      tracks.bindExternalIndices(&mcParticles);

      usedTracksIds.clear();
      for (auto& track : atracks) {
        auto ttrack = track.track_as<soa::Filtered<LabeledTracksEx>>();
        usedTracksIds.emplace_back(ttrack.globalIndex());
        if (ttrack.has_mcParticle()) {
          registry.fill(HIST("Tracks/Control/PtEfficiency"), ttrack.mcParticle_as<Particles>().pt());
        } else {
          registry.fill(HIST("Tracks/Control/PtEfficiencySecondaries"), ttrack.pt());
        }
      }
      for (auto& track : tracks) {
        if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
          continue;
        }
        if (track.has_mcParticle()) {
          registry.fill(HIST("Tracks/Control/PtEfficiency"), track.mcParticle_as<Particles>().pt());
        } else {
          registry.fill(HIST("Tracks/Control/PtEfficiencySecondaries"), track.pt());
        }
      }

      for (auto& particle : particles) {
        if (!particle.producedByGenerator()) {
          continue;
        }
        auto charge = 0.;
        auto p = pdg->GetParticle(particle.pdgCode());
        if (p != nullptr) {
          charge = p->Charge();
        }
        if (std::abs(charge) < 3.) {
          continue;
        }
        registry.fill(HIST("Tracks/Control/PtGen"), particle.pt());
      }
    }
  }

  PROCESS_SWITCH(MultiplicityCounter, processTrackEfficiency, "Calculate tracking efficiency vs pt", false);

  template <typename C>
  void processGenGeneral(
    aod::McCollisions::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks)
  {
    constexpr bool hasCentrality = C::template contains<aod::CentFT0Cs>() || C::template contains<aod::CentFT0Ms>();
    float c_rec = -1;
    float c_gen = -1;
    // add generated centrality estimation

    auto perCollisionMCSample = mcSample->sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
    auto nCharged = 0;
    for (auto& particle : perCollisionMCSample) {
      auto charge = 0.;
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < 3.) {
        continue;
      }
      nCharged++;
    }
    if constexpr (hasCentrality) {
      registry.fill(HIST("Events/Centrality/NtrkZvtxGen_t"), nCharged, mcCollision.posZ(), c_gen);
      registry.fill(HIST("Events/Centrality/Efficiency"), 1., c_gen);
    } else {
      registry.fill(HIST("Events/NtrkZvtxGen_t"), nCharged, mcCollision.posZ());
      registry.fill(HIST("Events/Efficiency"), 1.);
    }

    if (nCharged > 0) {
      if constexpr (hasCentrality) {
        registry.fill(HIST("Events/Centrality/Efficiency"), 2., c_gen);
      } else {
        registry.fill(HIST("Events/Efficiency"), 2.);
      }
    }
    bool atLeastOne = false;
    bool atLeastOne_gt0 = false;
    auto moreThanOne = 0;
    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());

    auto Nrec = 0;

    for (auto& collision : collisions) {
      c_rec = -1;
      if constexpr (hasCentrality) {
        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          c_rec = collision.centFT0C();
        } else if (C::template contains<aod::CentFT0Ms>()) {
          c_rec = collision.centFT0M();
        }
      }
      if constexpr (hasCentrality) {
        registry.fill(HIST("Events/Centrality/Efficiency"), 3., c_gen);
      } else {
        registry.fill(HIST("Events/Efficiency"), 3.);
      }
      if (!useEvSel || collision.sel8()) {
        auto perCollisionSample = tracks.sliceBy(perCol, collision.globalIndex());
        if constexpr (hasCentrality) {
          registry.fill(HIST("Events/Centrality/Efficiency"), 4., c_gen);
        } else {
          registry.fill(HIST("Events/Efficiency"), 4.);
        }

        if (perCollisionSample.size() > 0) {
          if constexpr (hasCentrality) {
            registry.fill(HIST("Events/Centrality/Efficiency"), 5., c_gen);
          } else {
            registry.fill(HIST("Events/Efficiency"), 5.);
          }
        }
        ++moreThanOne;
        atLeastOne = true;
        for (auto& t : perCollisionSample) {
          if (std::abs(t.eta()) < estimatorEta) {
            ++Nrec;
          }
        }
        if (perCollisionSample.size() > 0) {
          atLeastOne_gt0 = true;
        }
        if constexpr (hasCentrality) {
          registry.fill(HIST("Events/Centrality/NtrkZvtxGen"), Nrec, collision.posZ(), c_rec);
        } else {
          registry.fill(HIST("Events/NtrkZvtxGen"), Nrec, collision.posZ());
        }
      }
    }

    if (fillResponse) {
      for (auto& collision : collisions) {
        c_rec = -1;
        if constexpr (hasCentrality) {
          if constexpr (C::template contains<aod::CentFT0Cs>()) {
            c_rec = collision.centFT0C();
          } else if (C::template contains<aod::CentFT0Ms>()) {
            c_rec = collision.centFT0M();
          }
        }
        auto perCollisionSample = tracks.sliceBy(perCol, collision.globalIndex());
        Nrec = 0;
        for (auto& t : perCollisionSample) {
          if (std::abs(t.eta()) < estimatorEta) {
            ++Nrec;
          }
        }
        if constexpr (hasCentrality) {
          registry.fill(HIST("Events/Centrality/Response"), Nrec, nCharged, mcCollision.posZ(), c_rec);
          registry.fill(HIST("Events/Centrality/EfficiencyMult"), nCharged, mcCollision.posZ(), c_rec);
        } else {
          registry.fill(HIST("Events/Response"), Nrec, nCharged, mcCollision.posZ());
          registry.fill(HIST("Events/EfficiencyMult"), nCharged, mcCollision.posZ());
        }
      }
      if (moreThanOne > 1) {
        if constexpr (hasCentrality) {
          registry.fill(HIST("Events/Centrality/SplitMult"), nCharged, mcCollision.posZ(), c_gen);
        } else {
          registry.fill(HIST("Events/SplitMult"), nCharged, mcCollision.posZ());
        }
      }
    }

    if (collisions.size() == 0) {
      if constexpr (hasCentrality) {
        registry.fill(HIST("Events/Centrality/NotFoundEventZvtx"), mcCollision.posZ(), c_gen);
      } else {
        registry.fill(HIST("Events/NotFoundEventZvtx"), mcCollision.posZ());
      }
    }

    for (auto& particle : particles) {
      auto p = pdg->GetParticle(particle.pdgCode());
      auto charge = 0.;
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < 3.) {
        continue;
      }
      if constexpr (hasCentrality) {
        registry.fill(HIST("Tracks/Centrality/EtaZvtxGen_t"), particle.eta(), mcCollision.posZ(), c_gen);
        registry.fill(HIST("Tracks/Centrality/Control/PtEtaGen"), particle.pt(), particle.eta(), c_gen);
      } else {
        registry.fill(HIST("Tracks/EtaZvtxGen_t"), particle.eta(), mcCollision.posZ());
        registry.fill(HIST("Tracks/Control/PtEtaGen"), particle.pt(), particle.eta());
      }
      if (perCollisionMCSample.size() > 0) {
        if constexpr (hasCentrality) {
          registry.fill(HIST("Tracks/Centrality/EtaZvtxGen_gt0t"), particle.eta(), mcCollision.posZ(), c_gen);
        } else {
          registry.fill(HIST("Tracks/EtaZvtxGen_gt0t"), particle.eta(), mcCollision.posZ());
        }
      }
      if (atLeastOne) {
        if constexpr (hasCentrality) {
          registry.fill(HIST("Tracks/Centrality/EtaZvtxGen"), particle.eta(), mcCollision.posZ(), c_gen);
          if (atLeastOne_gt0) {
            registry.fill(HIST("Tracks/Centrality/EtaZvtxGen_gt0"), particle.eta(), mcCollision.posZ(), c_gen);
          }
          registry.fill(HIST("Tracks/Centrality/PhiEtaGen"), particle.phi(), particle.eta(), c_gen);
        } else {
          registry.fill(HIST("Tracks/EtaZvtxGen"), particle.eta(), mcCollision.posZ());
          if (atLeastOne_gt0) {
            registry.fill(HIST("Tracks/EtaZvtxGen_gt0"), particle.eta(), mcCollision.posZ());
          }
          registry.fill(HIST("Tracks/PhiEtaGen"), particle.phi(), particle.eta());
        }
      }
    }
  }

  void processGen(
    aod::McCollisions::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExCols, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks)
  {
    processGenGeneral<ExCols>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGen, "Process generator-level info", false);

  void processGenFT0C(
    aod::McCollisions::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks)
  {
    processGenGeneral<ExColsCentFT0C>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0C, "Process generator-level info (FT0C centrality)", false);

  void processGenFT0M(
    aod::McCollisions::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks)
  {
    processGenGeneral<ExColsCentFT0M>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0M, "Process generator-level info (FT0M centrality)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityCounter>(cfgc)};
}
