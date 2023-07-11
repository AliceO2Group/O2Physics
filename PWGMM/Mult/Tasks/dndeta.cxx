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
using ReTracks = soa::Join<aod::ReassignedTracksCore, aod::ReassignedTracksExtra>;

namespace
{
template <typename T>
static constexpr bool hasCent()
{
  if constexpr (!soa::is_soa_join_v<T>) {
    return false;
  } else if (T::template contains<aod::HepMCHeavyIons>()) {
    return true;
  } else {
    return false;
  }
}
} // namespace

struct MultiplicityCounter {
  SliceCache cache;
  Preslice<aod::Tracks> perCol = aod::track::collisionId;
  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<ReTracks> perColU = aod::track::bestCollisionId;

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
  std::vector<int> usedTracksIdsDFMC;
  std::vector<int> usedTracksIdsDFMCEff;

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
      hstat = registry.get<TH1>(HIST("Events/Selection"));
      x = hstat->GetXaxis();
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
      hstat = registry.get<TH2>(HIST("Events/Centrality/Selection"));
      x = hstat->GetXaxis();
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
      x = heff->GetXaxis();
      x->SetBinLabel(1, "Generated");
      x->SetBinLabel(2, "Generated INEL>0");
      x->SetBinLabel(3, "Reconstructed");
      x->SetBinLabel(4, "Selected");
      x->SetBinLabel(5, "Selected INEL>0");
    }

    if (doprocessGenFT0C || doprocessGenFT0M || doprocessGenFT0Chi || doprocessGenFT0Mhi) {
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
      x = heff->GetXaxis();
      x->SetBinLabel(1, "Generated");
      x->SetBinLabel(2, "Generated INEL>0");
      x->SetBinLabel(3, "Reconstructed");
      x->SetBinLabel(4, "Selected");
      x->SetBinLabel(5, "Selected INEL>0");
    }

    if (doprocessTrackEfficiency) {
      registry.add({"Tracks/Control/PtGen", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
      registry.add({"Tracks/Control/PtGenNoEtaCut", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
      registry.add({"Tracks/Control/PtEfficiency", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
      registry.add({"Tracks/Control/PtEfficiencyNoEtaCut", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
      registry.add({"Tracks/Control/PtEfficiencyFakes", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
    }
    if (doprocessTrackEfficiencyIndexed) {
      registry.add({"Tracks/Control/PtGenI", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
      registry.add({"Tracks/Control/PtGenINoEtaCut", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
      registry.add({"Tracks/Control/PtEfficiencyI", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
      registry.add({"Tracks/Control/PtEfficiencyINoEtaCut", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
      registry.add({"Tracks/Control/PtEfficiencyISecondaries", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
      registry.add({"Tracks/Control/PtEfficiencyISecondariesNoEtaCut", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxis}}});
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
      if (!useEvSel || (bc.selection_bit(evsel::kIsBBT0A) &&
                        bc.selection_bit(evsel::kIsBBT0C)) != 0) {
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
    usedTracksIdsDFMC.clear();
    usedTracksIdsDFMCEff.clear();
  }

  expressions::Filter trackSelectionProper = ((aod::track::trackCutFlag & trackSelectionITS) == trackSelectionITS) &&
                                             ifnode(ncheckbit(aod::track::detectorMap, (uint8_t)o2::aod::track::TPC),
                                                    ncheckbit(aod::track::trackCutFlag, trackSelectionTPC),
                                                    true) &&
                                             ncheckbit(aod::track::trackCutFlag, trackSelectionDCA);

  expressions::Filter atrackFilter = (aod::track::bestCollisionId >= 0) &&
                                     (nabs(aod::track::bestDCAZ) <= 2.f) &&
                                     (nabs(aod::track::bestDCAXY) <= ((0.004f + 0.013f / npow(aod::track::pts, 1.1f))));

  using ExTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
  using FiTracks = soa::Filtered<ExTracks>;
  using FiReTracks = soa::Filtered<ReTracks>;

  using ExCols = soa::Join<aod::Collisions, aod::EvSels>;

  template <typename C>
  void processCountingGeneral(
    typename C::iterator const& collision,
    FiTracks const& tracks,
    soa::SmallGroups<ReTracks> const& atracks)
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
    soa::SmallGroups<ReTracks> const& atracks)
  {
    processCountingGeneral<ExCols>(collision, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCounting, "Count tracks", false);

  using ExColsCentFT0C = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>;
  void processCountingCentralityFT0C(
    ExColsCentFT0C::iterator const& collision,
    FiTracks const& tracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    processCountingGeneral<ExColsCentFT0C>(collision, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCountingCentralityFT0C, "Count tracks in centrality bins", false);

  using ExColsCentFT0M = soa::Join<aod::Collisions, aod::CentFT0Ms, aod::EvSels>;
  void processCountingCentralityFT0M(
    ExColsCentFT0M::iterator const& collision,
    FiTracks const& tracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    processCountingGeneral<ExColsCentFT0M>(collision, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCountingCentralityFT0M, "Count tracks in centrality bins", false);

  using Particles = soa::Filtered<aod::McParticles>;
  using LabeledTracksEx = soa::Join<LabeledTracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
  using FiLTracks = soa::Filtered<LabeledTracksEx>;
  using ParticlesI = soa::Filtered<soa::Join<aod::McParticles, aod::ParticlesToTracks>>;
  expressions::Filter primaries = (aod::mcparticle::flags & (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary) == (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary;

  template <typename C, typename MC>
  void processTrackEfficiencyIndexedGeneral(
    typename soa::Join<C, aod::McCollisionLabels>::iterator const& collision,
    MC const&, ParticlesI const& particles,
    FiLTracks const& tracks)
  {
    if (useEvSel && !collision.sel8()) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }
    auto mcCollision = collision.mcCollision();
    auto sample = particles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);

    for (auto& particle : sample) {
      auto charge = 0.;
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < 3.) {
        continue;
      }
      registry.fill(HIST("Tracks/Control/PtGenINoEtaCut"), particle.pt());
      if (std::abs(particle.eta()) < estimatorEta) {
        registry.fill(HIST("Tracks/Control/PtGenI"), particle.pt());
      }
      if (particle.has_tracks()) {
        auto counted = false;
        auto countedNoEtaCut = false;
        auto counter = 0;
        auto relatedTracks = particle.template filtered_tracks_as<FiLTracks>();
        for (auto const& track : relatedTracks) {
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
            registry.fill(HIST("Tracks/Control/PtEfficiencyISecondariesNoEtaCut"), particle.pt());
            if (std::abs(track.eta()) < estimatorEta) {
              registry.fill(HIST("Tracks/Control/PtEfficiencyISecondaries"), particle.pt());
            }
          }
        }
        if (counter > 1) {
          for (auto const& track : relatedTracks) {
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
          for (auto const& track : relatedTracks) {
            registry.fill(HIST("Tracks/Control/PhiEtaDuplicates"), track.phi(), track.eta());
          }
        }
      }
    }
  }

  void processTrackEfficiencyIndexed(
    soa::Join<ExCols, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, ParticlesI const& particles,
    FiLTracks const& tracks)
  {
    processTrackEfficiencyIndexedGeneral<ExCols, aod::McCollisions>(collision, mccollisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processTrackEfficiencyIndexed, "Calculate tracking efficiency vs pt (indexed)", false);

  template <typename C, typename MC>
  void processTrackEfficiencyGeneral(
    typename soa::Join<C, aod::McCollisionLabels>::iterator const& collision,
    MC const&, Particles const& particles,
    FiLTracks const& tracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    constexpr bool hasCentrality = C::template contains<aod::CentFT0Cs>() || C::template contains<aod::CentFT0Ms>() || hasCent<MC>();

    if (useEvSel && !collision.sel8()) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }
    float c_rec = -1;
    float c_gen = -1;
    if constexpr (hasCentrality) {
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        c_rec = collision.centFT0C();
      } else if (C::template contains<aod::CentFT0Ms>()) {
        c_rec = collision.centFT0M();
      }
    }
    auto mcCollision = collision.mcCollision();
    if constexpr (hasCent<MC>()) {
      c_gen = mcCollision.centrality();
    }

    auto particlesPerCol = particles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);

    usedTracksIds.clear();
    for (auto const& track : atracks) {
      auto otrack = track.template track_as<FiLTracks>();
      usedTracksIds.emplace_back(track.trackId());
      if (otrack.collisionId() != track.bestCollisionId()) {
        usedTracksIdsDFMCEff.emplace_back(track.trackId());
      }
      if (otrack.has_mcParticle()) {
        auto particle = otrack.mcParticle_as<Particles>();
        registry.fill(HIST("Tracks/Control/PtEfficiencyNoEtaCut"), particle.pt());
        if (std::abs(otrack.eta()) < estimatorEta) {
          registry.fill(HIST("Tracks/Control/PtEfficiency"), particle.pt());
        }
      } else {
        registry.fill(HIST("Tracks/Control/PtEfficiencyFakes"), otrack.pt());
      }
    }
    for (auto const& track : tracks) {
      if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
        continue;
      }
      if (std::find(usedTracksIdsDFMCEff.begin(), usedTracksIdsDFMCEff.end(), track.globalIndex()) != usedTracksIdsDFMCEff.end()) {
        continue;
      }
      if (track.has_mcParticle()) {
        auto particle = track.template mcParticle_as<Particles>();
        registry.fill(HIST("Tracks/Control/PtEfficiencyNoEtaCut"), particle.pt());
        if (std::abs(track.eta()) < estimatorEta) {
          registry.fill(HIST("Tracks/Control/PtEfficiency"), particle.pt());
        }
      } else {
        registry.fill(HIST("Tracks/Control/PtEfficiencyFakes"), track.pt());
      }
    }

    for (auto& particle : particlesPerCol) {
      auto charge = 0.;
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < 3.) {
        continue;
      }
      registry.fill(HIST("Tracks/Control/PtGenNoEtaCut"), particle.pt());
      if (std::abs(particle.eta()) < estimatorEta) {
        registry.fill(HIST("Tracks/Control/PtGen"), particle.pt());
      }
    }
  }

  void processTrackEfficiency(
    soa::Join<ExCols, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, Particles const& mcParticles,
    FiLTracks const& filtracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    processTrackEfficiencyGeneral<ExCols, aod::McCollisions>(collision, mccollisions, mcParticles, filtracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processTrackEfficiency, "Calculate tracking efficiency vs pt", false);

  template <typename MC, typename C>
  void processGenGeneral(
    typename MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks)
  {
    constexpr bool hasCentrality = C::template contains<aod::CentFT0Cs>() || C::template contains<aod::CentFT0Ms>() || hasCent<MC>();

    float c_rec = -1;
    float c_gen = -1;
    // add generated centrality estimation
    if constexpr (hasCent<MC>()) {
      c_gen = mcCollision.centrality();
    }

    auto nCharged = 0;
    for (auto& particle : particles) {
      auto charge = 0.;
      auto p = pdg->GetParticle(particle.pdgCode());
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < 3.) {
        continue;
      }
      if (std::abs(particle.eta()) >= estimatorEta) {
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
    std::vector<int> NrecPerCol;
    std::vector<float> c_recPerCol;

    for (auto& collision : collisions) {
      usedTracksIds.clear();
      c_rec = -1;
      if constexpr (hasCentrality) {
        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          c_rec = collision.centFT0C();
        } else if (C::template contains<aod::CentFT0Ms>()) {
          c_rec = collision.centFT0M();
        }
        c_recPerCol.emplace_back(c_rec);
        registry.fill(HIST("Events/Centrality/Efficiency"), 3., c_gen);
      } else {
        registry.fill(HIST("Events/Efficiency"), 3.);
      }
      if (!useEvSel || collision.sel8()) {
        Nrec = 0;
        ++moreThanOne;
        atLeastOne = true;

        auto perCollisionSample = tracks.sliceBy(perCol, collision.globalIndex());
        auto perCollisionASample = atracks.sliceBy(perColU, collision.globalIndex());
        for (auto const& track : perCollisionASample) {
          auto otrack = track.template track_as<FiTracks>();
          usedTracksIds.emplace_back(track.trackId());
          if (otrack.collisionId() != track.bestCollisionId()) {
            usedTracksIdsDFMC.emplace_back(track.trackId());
          }
          if (std::abs(otrack.eta()) < estimatorEta) {
            ++Nrec;
          }
        }
        for (auto const& track : perCollisionSample) {
          if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
            continue;
          }
          if (std::find(usedTracksIdsDFMC.begin(), usedTracksIdsDFMC.end(), track.globalIndex()) != usedTracksIdsDFMC.end()) {
            continue;
          }
          if (std::abs(track.eta()) < estimatorEta) {
            ++Nrec;
          }
        }
        NrecPerCol.emplace_back(Nrec);

        if constexpr (hasCentrality) {
          registry.fill(HIST("Events/Centrality/Efficiency"), 4., c_gen);
        } else {
          registry.fill(HIST("Events/Efficiency"), 4.);
        }

        if (Nrec > 0) {
          if constexpr (hasCentrality) {
            registry.fill(HIST("Events/Centrality/Efficiency"), 5., c_gen);
          } else {
            registry.fill(HIST("Events/Efficiency"), 5.);
          }
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
      for (auto i = 0U; i < NrecPerCol.size(); ++i) {
        if constexpr (hasCentrality) {
          registry.fill(HIST("Events/Centrality/Response"), NrecPerCol[i], nCharged, mcCollision.posZ(), c_recPerCol[i]);
          registry.fill(HIST("Events/Centrality/EfficiencyMult"), nCharged, mcCollision.posZ(), c_recPerCol[i]);
        } else {
          registry.fill(HIST("Events/Response"), NrecPerCol[i], nCharged, mcCollision.posZ());
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
      if (nCharged > 0) {
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

  using MC = aod::McCollisions; // soa::Join<aod::McCollisions, aod::HepMCXSections>;
  void processGen(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExCols, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks)
  {
    processGenGeneral<MC, ExCols>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGen, "Process generator-level info", false);

  void processGenFT0C(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks)
  {
    processGenGeneral<MC, ExColsCentFT0C>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0C, "Process generator-level info (FT0C centrality)", false);

  void processGenFT0M(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks)
  {
    processGenGeneral<MC, ExColsCentFT0M>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0M, "Process generator-level info (FT0M centrality)", false);

  using MChi = soa::Join<aod::McCollisions, aod::HepMCHeavyIons>;

  void processGenFT0Chi(
    MChi::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks)
  {
    processGenGeneral<MChi, ExColsCentFT0C>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0Chi, "Process generator-level info (FT0C centrality, HI)", false);

  void processGenFT0Mhi(
    MChi::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks)
  {
    processGenGeneral<MChi, ExColsCentFT0M>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0Mhi, "Process generator-level info (FT0M centrality, HI)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityCounter>(cfgc)};
}
