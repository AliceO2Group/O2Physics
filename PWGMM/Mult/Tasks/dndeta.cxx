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
#include <TDatabasePDG.h>
#include <TPDGCode.h>

#include "bestCollisionTable.h"
#include "Axes.h"
#include "Functions.h"
#include "Selections.h"
#include "Histograms.h"

#include <ranges>

using namespace o2;
using namespace o2::aod::track;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace pwgmm::mult;
using namespace pwgmm::mult::histograms;

using LabeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
using ReTracks = soa::Join<aod::ReassignedTracksCore, aod::ReassignedTracksExtra>;

struct MultiplicityCounter {
  SliceCache cache;
  Preslice<aod::Tracks> perCol = aod::track::collisionId;
  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;
  PresliceUnsorted<ReTracks> perColU = aod::track::bestCollisionId;

  Service<o2::framework::O2DatabasePDG> pdg;

  Configurable<float> estimatorEta{"estimatorEta", 1.0, "eta range for INEL>0 sample definition"};
  Configurable<float> dcaZ{"dcaZ", 0.2f, "Custom DCA Z cut (ignored if negative)"};
  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  Configurable<bool> fillResponse{"fillResponse", false, "Fill response matrix"};
  Configurable<bool> responseStudy{"responseStudy", false, "Fill multi-estimator response"};
  ConfigurableAxis multBinning{"multBinning", {301, -0.5, 300.5}, ""};
  ConfigurableAxis centBinning{"centBinning", {VARIABLE_WIDTH, 0, 10, 20, 30, 40, 50, 60, 70, 80, 100}, ""};

  HistogramRegistry commonRegistry{
    "Common",
    {
      {"Events/BCSelection", ";status;count", {HistType::kTH1F, {{3, 0.5, 3.5}}}} //
    } //
  };

  HistogramRegistry inclusiveRegistry{
    InclusivePrefix.data(),
    {}
  };
  HistogramRegistry binnedRegistry{
    BinnedPrefix.data(),
    {}
  };

  std::vector<int> usedTracksIds;
  std::vector<int> usedTracksIdsDF;
  std::vector<int> usedTracksIdsDFMC;
  std::vector<int> usedTracksIdsDFMCEff;

  void init(InitContext&)
  {
    AxisSpec MultAxis = {multBinning};
    AxisSpec CentAxis = {centBinning, "centrality"};

    auto hstat = inclusiveRegistry.get<TH1>(HIST("Events/BCSelection"));
    auto* x = hstat->GetXaxis();
    x->SetBinLabel(1, "Good BCs");
    x->SetBinLabel(2, "BCs with collisions");
    x->SetBinLabel(3, "BCs with pile-up/splitting");

    if (doprocessEventStat) {
      inclusiveRegistry.add({"Events/Control/Chi2", " ; #chi^2", {HistType::kTH1F, {{101, -0.1, 10.1}}}});
      inclusiveRegistry.add({"Events/Control/TimeResolution", " ; t (ms)", {HistType::kTH1F, {{1001, -0.1, 100.1}}}});
    }
    if (doprocessEventStatCentralityFT0C || doprocessEventStatCentralityFT0M) {
      binnedRegistry.add({"Events/Control/Chi2", " ; #chi^2; centrality", {HistType::kTH2F, {{101, -0.1, 10.1}, CentAxis}}});
      binnedRegistry.add({"Events/Control/TimeResolution", " ; t (ms); centrality", {HistType::kTH2F, {{1001, -0.1, 100.1}, CentAxis}}});
    }

    if (doprocessCounting || doprocessCountingNoAmb) {
      inclusiveRegistry.add({"Events/Selection", ";status;events", {HistType::kTH1F, {{static_cast<int>(EvSelBins::kRejected), 0.5, static_cast<float>(EvSelBins::kRejected) + 0.5}}}});
      hstat = inclusiveRegistry.get<TH1>(HIST("Events/Selection"));
      x = hstat->GetXaxis();
      x->SetBinLabel(static_cast<int>(EvSelBins::kAll), "All");
      x->SetBinLabel(static_cast<int>(EvSelBins::kSelected), "Selected");
      x->SetBinLabel(static_cast<int>(EvSelBins::kSelectedgt0), "Selected INEL>0");
      x->SetBinLabel(static_cast<int>(EvSelBins::kSelectedPVgt0), "Selected INEL>0 (PV)");
      x->SetBinLabel(static_cast<int>(EvSelBins::kRejected), "Rejected");

      inclusiveRegistry.add({"Events/NtrkZvtx", "; N_{trk}; Z_{vtx} (cm); events", {HistType::kTH2F, {MultAxis, ZAxis}}});
      inclusiveRegistry.add({"Tracks/EtaZvtx", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({"Tracks/EtaZvtx_gt0", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({"Tracks/EtaZvtx_PVgt0", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({"Tracks/PhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      inclusiveRegistry.add({"Tracks/Control/PtEta", " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, EtaAxis}}});
      inclusiveRegistry.add({"Tracks/Control/DCAXYPt", " ; p_{T} (GeV/c) ; DCA_{XY} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
      inclusiveRegistry.add({"Tracks/Control/DCAZPt", " ; p_{T} (GeV/c) ; DCA_{Z} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
      if (doprocessCounting) {
        inclusiveRegistry.add({"Tracks/Control/ReassignedDCAXYPt", " ; p_{T} (GeV/c) ; DCA_{XY} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
        inclusiveRegistry.add({"Tracks/Control/ReassignedDCAZPt", " ; p_{T} (GeV/c) ; DCA_{Z} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
        inclusiveRegistry.add({"Tracks/Control/ExtraDCAXYPt", " ; p_{T} (GeV/c) ; DCA_{XY} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
        inclusiveRegistry.add({"Tracks/Control/ExtraDCAZPt", " ; p_{T} (GeV/c) ; DCA_{Z} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
        inclusiveRegistry.add({"Tracks/Control/ExtraTracksEtaZvtx", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
        inclusiveRegistry.add({"Tracks/Control/ExtraTracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
        inclusiveRegistry.add({"Tracks/Control/ReassignedTracksEtaZvtx", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
        inclusiveRegistry.add({"Tracks/Control/ReassignedTracksPhiEta", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
        inclusiveRegistry.add({"Tracks/Control/ReassignedVertexCorr", "; Z_{vtx}^{orig} (cm); Z_{vtx}^{re} (cm)", {HistType::kTH2F, {ZAxis, ZAxis}}});
      }
    }

    if (doprocessCountingCentralityFT0C || doprocessCountingCentralityFT0M || doprocessCountingCentralityFT0CNoAmb || doprocessCountingCentralityFT0MNoAmb) {
      binnedRegistry.add({"Events/Selection", ";status;centrality;events", {HistType::kTH2F, {{3, 0.5, 3.5}, CentAxis}}});
      hstat = binnedRegistry.get<TH2>(HIST("Events/Selection"));
      x = hstat->GetXaxis();
      x->SetBinLabel(1, "All");
      x->SetBinLabel(2, "Selected");
      x->SetBinLabel(3, "Rejected");

      binnedRegistry.add({"Events/NtrkZvtx", "; N_{trk}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({"Tracks/EtaZvtx", "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({"Tracks/PhiEta", "; #varphi; #eta; centrality", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
      binnedRegistry.add({"Tracks/Control/PtEta", " ; p_{T} (GeV/c); #eta; centrality", {HistType::kTHnSparseF, {PtAxis, EtaAxis, CentAxis}}});
      binnedRegistry.add({"Tracks/Control/DCAXYPt", " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
      binnedRegistry.add({"Tracks/Control/DCAZPt", " ; p_{T} (GeV/c) ; DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
      if (doprocessCountingCentralityFT0C || doprocessCountingCentralityFT0M) {
        binnedRegistry.add({"Tracks/Control/ReassignedDCAXYPt", " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
        binnedRegistry.add({"Tracks/Control/ReassignedDCAZPt", " ; p_{T} (GeV/c) ; DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
        binnedRegistry.add({"Tracks/Control/ExtraDCAXYPt", " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
        binnedRegistry.add({"Tracks/Control/ExtraDCAZPt", " ; p_{T} (GeV/c) ; DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
        binnedRegistry.add({"Tracks/Control/ExtraTracksEtaZvtx", "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
        binnedRegistry.add({"Tracks/Control/ExtraTracksPhiEta", "; #varphi; #eta; centrality", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
        binnedRegistry.add({"Tracks/Control/ReassignedTracksEtaZvtx", "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
        binnedRegistry.add({"Tracks/Control/ReassignedTracksPhiEta", "; #varphi; #eta; centrality", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
        binnedRegistry.add({"Tracks/Control/ReassignedVertexCorr", "; Z_{vtx}^{orig} (cm); Z_{vtx}^{re} (cm); centrality", {HistType::kTHnSparseF, {ZAxis, ZAxis, CentAxis}}});
      }
    }

    if (doprocessGen || doprocessGenNoAmb) {
      inclusiveRegistry.add({"Events/NtrkZvtxGen", "; N_{trk}; Z_{vtx} (cm); events", {HistType::kTH2F, {MultAxis, ZAxis}}});
      inclusiveRegistry.add({"Events/NtrkZvtxGen_t", "; N_{part}; Z_{vtx} (cm); events", {HistType::kTH2F, {MultAxis, ZAxis}}});
      inclusiveRegistry.add({"Tracks/EtaZvtxGen", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({"Tracks/EtaZvtxGen_t", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({"Tracks/EtaZvtxGen_gt0", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({"Tracks/EtaZvtxGen_PVgt0", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({"Tracks/EtaZvtxGen_gt0t", "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({"Tracks/Control/PtEtaGen", " ; p_{T} (GeV/c) ; #eta", {HistType::kTH2F, {PtAxis, EtaAxis}}});

      inclusiveRegistry.add({"Tracks/PhiEtaGen", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});

      inclusiveRegistry.add({"Events/Efficiency", "; status; events", {HistType::kTH1F, {{static_cast<int>(EvEffBins::kSelectedPVgt0), 0.5, static_cast<float>(EvEffBins::kSelectedPVgt0) + 0.5}}}});
      inclusiveRegistry.add({"Events/NotFoundEventZvtx", " ; Z_{vtx} (cm)", {HistType::kTH1F, {ZAxis}}});

      if (fillResponse) {
        inclusiveRegistry.add({"Events/Response", " ; N_{rec}; N_{gen}; Z_{vtx} (cm)", {HistType::kTHnSparseF, {MultAxis, MultAxis, ZAxis}}});
        inclusiveRegistry.add({"Events/EfficiencyMult", " ; N_{gen}; Z_{vtx} (cm)", {HistType::kTH2F, {MultAxis, ZAxis}}});
        inclusiveRegistry.add({"Events/SplitMult", " ; N_{gen} ; Z_{vtx} (cm)", {HistType::kTH2F, {MultAxis, ZAxis}}});
        if (responseStudy) {
          inclusiveRegistry.add({"Events/Control/MultiResponse", " ; N_{gen}; N_{rec}; N_{PV cont}; N_{FT0A}; N_{FT0C}; N_{FDA}; N_{FDC}; Z_{vtx} (cm)", {HistType::kTHnSparseF, {MultAxis, MultAxis, MultAxis, FT0AAxis, FT0CAxis, FDDAxis, FDDAxis, ZAxis}}});
        }
      }

      auto heff = inclusiveRegistry.get<TH1>(HIST("Events/Efficiency"));
      x = heff->GetXaxis();
      x->SetBinLabel(static_cast<int>(EvEffBins::kGen), "Generated");
      x->SetBinLabel(static_cast<int>(EvEffBins::kGengt0), "Generated INEL>0");
      x->SetBinLabel(static_cast<int>(EvEffBins::kRec), "Reconstructed");
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelected), "Selected");
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelectedgt0), "Selected INEL>0");
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelectedPVgt0), "Selected INEL>0 (PV)");
    }

    if (doprocessGenFT0C || doprocessGenFT0M || doprocessGenFT0Chi || doprocessGenFT0Mhi ||
        doprocessGenFT0CNoAmb || doprocessGenFT0MNoAmb || doprocessGenFT0ChiNoAmb || doprocessGenFT0MhiNoAmb) {
      binnedRegistry.add({"Events/NtrkZvtxGen", "; N_{trk}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({"Events/NtrkZvtxGen_t", "; N_{part}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({"Tracks/EtaZvtxGen", "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({"Tracks/EtaZvtxGen_t", "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({"Tracks/EtaZvtxGen_gt0", "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({"Tracks/EtaZvtxGen_PVgt0", "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({"Tracks/EtaZvtxGen_gt0t", "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({"Tracks/Control/PtEtaGen", " ; p_{T} (GeV/c) ; #eta; centrality", {HistType::kTHnSparseF, {PtAxis, EtaAxis, CentAxis}}});

      binnedRegistry.add({"Tracks/PhiEtaGen", "; #varphi; #eta; tracks", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
      binnedRegistry.add({"Tracks/Control/PhiEtaGenDuplicates", "; #varphi; #eta; centrality", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
      binnedRegistry.add({"Tracks/Control/PhiEtaDuplicates", "; #varphi; #eta; centrality", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
      binnedRegistry.add({"Events/Efficiency", "; status; centrality; events", {HistType::kTH2F, {{static_cast<int>(EvEffBins::kSelectedPVgt0), 0.5, static_cast<float>(EvEffBins::kSelectedPVgt0) + 0.5}, CentAxis}}});
      binnedRegistry.add({"Events/NotFoundEventZvtx", " ; Z_{vtx} (cm); centrality; events", {HistType::kTH2F, {ZAxis, CentAxis}}});

      if (fillResponse) {
        binnedRegistry.add({"Events/Response", " ; N_{rec}; N_{gen}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, MultAxis, ZAxis, CentAxis}}});
        binnedRegistry.add({"Events/EfficiencyMult", " ; N_{gen}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});
        binnedRegistry.add({"Events/SplitMult", " ; N_{gen} ; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});
        if (responseStudy) {
          binnedRegistry.add({"Events/Control/MultiResponse", " ; N_{gen}; N_{rec}, N_{PV cont}; N_{FT0A}; N_{FT0C}; N_{FDA}; N_{FDC}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, MultAxis, MultAxis, FT0AAxis, FT0CAxis, FDDAxis, FDDAxis, ZAxis, CentAxis}}});
        }
      }

      auto heff = binnedRegistry.get<TH2>(HIST("Events/Efficiency"));
      x = heff->GetXaxis();
      x->SetBinLabel(static_cast<int>(EvEffBins::kGen), "Generated");
      x->SetBinLabel(static_cast<int>(EvEffBins::kGengt0), "Generated INEL>0");
      x->SetBinLabel(static_cast<int>(EvEffBins::kRec), "Reconstructed");
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelected), "Selected");
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelectedgt0), "Selected INEL>0");
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelectedPVgt0), "Selected INEL>0 (PV)");
    }

    if (doprocessTrackEfficiency || doprocessTrackEfficiencyNoAmb) {
      inclusiveRegistry.add({"Tracks/Control/PtGen", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({"Tracks/Control/PtGenNoEtaCut", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({"Tracks/Control/PtEfficiency", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({"Tracks/Control/PtEfficiencyNoEtaCut", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({"Tracks/Control/PtEfficiencyFakes", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      for (auto i = 0u; i < speciesIds.size(); ++i) {
        inclusiveRegistry.add({(std::string("Tracks/Control/") + std::string(species[i]) + "/PtGen").c_str(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
        inclusiveRegistry.add({(std::string("Tracks/Control/") + std::string(species[i]) + "/PtEfficiency").c_str(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      }
    }
    if (doprocessTrackEfficiencyIndexed) {
      inclusiveRegistry.add({"Tracks/Control/PhiEtaGenDuplicates", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      inclusiveRegistry.add({"Tracks/Control/PhiEtaDuplicates", "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      inclusiveRegistry.add({"Tracks/Control/PtGenI", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({"Tracks/Control/PtGenINoEtaCut", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({"Tracks/Control/PtEfficiencyI", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({"Tracks/Control/PtEfficiencyINoEtaCut", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({"Tracks/Control/PtEfficiencyISecondaries", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({"Tracks/Control/PtEfficiencyISecondariesNoEtaCut", " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({"Tracks/Control/Mask", " ; bit", {HistType::kTH1F, {{17, -0.5, 16.5}}}});
      inclusiveRegistry.add({"Tracks/Control/ITSClusters", " ; layer", {HistType::kTH1F, {{8, 0.5, 8.5}}}});
      for (auto i = 0u; i < speciesIds.size(); ++i) {
        inclusiveRegistry.add({(std::string("Tracks/Control/") + std::string(species[i]) + "/PtGenI").c_str(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
        inclusiveRegistry.add({(std::string("Tracks/Control/") + std::string(species[i]) + "/PtEfficiencyI").c_str(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      }
    }
  }

  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  template <typename C>
  void processEventStatGeneral(FullBCs const& bcs, C const& collisions)
  {
    constexpr bool hasCentrality = C::template contains<aod::CentFT0Cs>() || C::template contains<aod::CentFT0Ms>();
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto& bc : bcs) {
      if (!useEvSel || (bc.selection_bit(aod::evsel::kIsBBT0A) &&
                        bc.selection_bit(aod::evsel::kIsBBT0C)) != 0) {
        commonRegistry.fill(HIST("Events/BCSelection"), 1.);
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
          commonRegistry.fill(HIST("Events/BCSelection"), 2.);
          if (cols.size() > 1) {
            commonRegistry.fill(HIST("Events/BCSelection"), 3.);
          }
        }
        for (auto& col : cols) {
          if constexpr (hasCentrality) {
            float c = -1;
            if constexpr (C::template contains<aod::CentFT0Cs>()) {
              c = col.centFT0C();
            } else if constexpr (C::template contains<aod::CentFT0Ms>()) {
              c = col.centFT0M();
            }
            binnedRegistry.fill(HIST("Events/Control/Chi2"), col.chi2(), c);
            binnedRegistry.fill(HIST("Events/Control/TimeResolution"), col.collisionTimeRes(), c);
          } else {
            inclusiveRegistry.fill(HIST("Events/Control/Chi2"), col.chi2());
            inclusiveRegistry.fill(HIST("Events/Control/TimeResolution"), col.collisionTimeRes());
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

  // require a mix of ITS+TPC and ITS-only tracks
  expressions::Filter fTrackSelectionITS = ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::ITS) &&
                                           ncheckbit(aod::track::trackCutFlag, trackSelectionITS);
  expressions::Filter fTrackSelectionTPC = ifnode(ncheckbit(aod::track::v001::detectorMap, (uint8_t)o2::aod::track::TPC),
                                                  ncheckbit(aod::track::trackCutFlag, trackSelectionTPC), true);
  expressions::Filter fTrackSelectionDCA = ifnode(dcaZ.node() > 0.f, nabs(aod::track::dcaZ) <= dcaZ && ncheckbit(aod::track::trackCutFlag, trackSelectionDCAXYonly),
                                                  ncheckbit(aod::track::trackCutFlag, trackSelectionDCA));

  expressions::Filter fAtrackAssigned = (aod::track::bestCollisionId >= 0);
  expressions::Filter fAtrackDCA = ifnode(dcaZ.node() > 0.f, nabs(aod::track::bestDCAZ) <= dcaZ, nabs(aod::track::bestDCAZ) <= 2.0f) &&
                                   nabs(aod::track::bestDCAXY) <= ((0.004f + 0.013f / npow(aod::track::pts, 1.1f)));

  using ExTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
  using FiTracks = soa::Filtered<ExTracks>;
  using FiReTracks = soa::Filtered<ReTracks>;

  using ExCols = soa::Join<aod::Collisions, aod::EvSels>;

  // PV contributors for INEL>0 (PV) collision sample definition
  Partition<FiTracks> pvContribTracksIUEta1 = (nabs(aod::track::eta) < 1.0f) && ((aod::track::flags & (uint32_t)o2::aod::track::PVContributor) == (uint32_t)o2::aod::track::PVContributor);

  template <typename C>
  void processCountingGeneral(
    typename C::iterator const& collision,
    FiTracks const& tracks)
  {
    float c = -1;
    constexpr bool hasCentrality = C::template contains<aod::CentFT0Cs>() || C::template contains<aod::CentFT0Ms>();
    if constexpr (hasCentrality) {
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        c = collision.centFT0C();
      } else if (C::template contains<aod::CentFT0Ms>()) {
        c = collision.centFT0M();
      }
      binnedRegistry.fill(HIST("Events/Selection"), 1., c);
    } else {
      inclusiveRegistry.fill(HIST("Events/Selection"), static_cast<float>(EvSelBins::kAll));
    }

    if (!useEvSel || collision.sel8()) {
      if constexpr (hasCentrality) {
        binnedRegistry.fill(HIST("Events/Selection"), 2., c);
      } else {
        inclusiveRegistry.fill(HIST("Events/Selection"), static_cast<float>(EvSelBins::kSelected));
      }
      auto z = collision.posZ();
      usedTracksIds.clear();

      auto groupPVContrib = pvContribTracksIUEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      auto Ntrks = 0;
      for (auto& track : tracks) {
        if (std::abs(track.eta()) < estimatorEta) {
          ++Ntrks;
        }
        if constexpr (hasCentrality) {
          binnedRegistry.fill(HIST("Tracks/EtaZvtx"), track.eta(), z, c);
          binnedRegistry.fill(HIST("Tracks/PhiEta"), track.phi(), track.eta(), c);
          binnedRegistry.fill(HIST("Tracks/Control/PtEta"), track.pt(), track.eta(), c);
          binnedRegistry.fill(HIST("Tracks/Control/DCAXYPt"), track.pt(), track.dcaXY(), c);
          binnedRegistry.fill(HIST("Tracks/Control/DCAZPt"), track.pt(), track.dcaZ(), c);
        } else {
          inclusiveRegistry.fill(HIST("Tracks/EtaZvtx"), track.eta(), z);
          inclusiveRegistry.fill(HIST("Tracks/PhiEta"), track.phi(), track.eta());
          inclusiveRegistry.fill(HIST("Tracks/Control/PtEta"), track.pt(), track.eta());
          inclusiveRegistry.fill(HIST("Tracks/Control/DCAXYPt"), track.pt(), track.dcaXY());
          inclusiveRegistry.fill(HIST("Tracks/Control/DCAZPt"), track.pt(), track.dcaZ());
        }
      }
      if constexpr (hasCentrality) {
        binnedRegistry.fill(HIST("Events/NtrkZvtx"), Ntrks, z, c);
      } else {
        if (Ntrks > 0 || groupPVContrib.size() > 0) {
          if (groupPVContrib.size() > 0) {
            inclusiveRegistry.fill(HIST("Events/Selection"), static_cast<float>(EvSelBins::kSelectedPVgt0));
          }
          if (Ntrks > 0) {
            inclusiveRegistry.fill(HIST("Events/Selection"), static_cast<float>(EvSelBins::kSelectedgt0));
          }
          for (auto& track : tracks) {
            if (Ntrks > 0) {
              inclusiveRegistry.fill(HIST("Tracks/EtaZvtx_gt0"), track.eta(), z);
            }
            if (groupPVContrib.size() > 0) {
              inclusiveRegistry.fill(HIST("Tracks/EtaZvtx_PVgt0"), track.eta(), z);
            }
          }
        }
        inclusiveRegistry.fill(HIST("Events/NtrkZvtx"), Ntrks, z);
      }
    } else {
      if constexpr (hasCentrality) {
        binnedRegistry.fill(HIST("Events/Selection"), 3., c);
      } else {
        inclusiveRegistry.fill(HIST("Events/Selection"), static_cast<float>(EvSelBins::kRejected));
      }
    }
  }

  template <typename C>
  void processCountingGeneralwAmbiguous(
    typename C::iterator const& collision,
    FiTracks const& tracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    float c = -1;
    if constexpr (hasRecoCent<C>()) {
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        c = collision.centFT0C();
      } else if (C::template contains<aod::CentFT0Ms>()) {
        c = collision.centFT0M();
      }
      binnedRegistry.fill(HIST("Events/Selection"), 1., c);
    } else {
      inclusiveRegistry.fill(HIST("Events/Selection"), static_cast<float>(EvSelBins::kAll));
    }

    if (!useEvSel || collision.sel8()) {
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST("Events/Selection"), 2., c);
      } else {
        inclusiveRegistry.fill(HIST("Events/Selection"), static_cast<float>(EvSelBins::kSelected));
      }
      auto z = collision.posZ();
      usedTracksIds.clear();

      auto groupPVContrib = pvContribTracksIUEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      auto Ntrks = 0;

      for (auto& track : atracks) {
        auto otrack = track.track_as<FiTracks>();
        // same filtering for ambiguous as for general
        if (!otrack.hasITS()) {
          continue;
        }
        if ((otrack.trackCutFlag() & trackSelectionITS) == 0) {
          continue;
        }
        if (otrack.hasTPC()) {
          if ((otrack.trackCutFlag() & trackSelectionTPC) == 0) {
            continue;
          }
        }
        usedTracksIds.emplace_back(track.trackId());
        if (std::abs(otrack.eta()) < estimatorEta) {
          ++Ntrks;
        }
        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST("Tracks/EtaZvtx"), otrack.eta(), z, c);
          binnedRegistry.fill(HIST("Tracks/PhiEta"), otrack.phi(), otrack.eta(), c);
          binnedRegistry.fill(HIST("Tracks/Control/PtEta"), otrack.pt(), otrack.eta(), c);
          binnedRegistry.fill(HIST("Tracks/Control/DCAXYPt"), otrack.pt(), track.bestDCAXY(), c);
          binnedRegistry.fill(HIST("Tracks/Control/DCAZPt"), otrack.pt(), track.bestDCAZ(), c);
        } else {
          inclusiveRegistry.fill(HIST("Tracks/EtaZvtx"), otrack.eta(), z);
          inclusiveRegistry.fill(HIST("Tracks/PhiEta"), otrack.phi(), otrack.eta());
          inclusiveRegistry.fill(HIST("Tracks/Control/PtEta"), otrack.pt(), otrack.eta());
          inclusiveRegistry.fill(HIST("Tracks/Control/DCAXYPt"), otrack.pt(), track.bestDCAXY());
          inclusiveRegistry.fill(HIST("Tracks/Control/DCAZPt"), otrack.pt(), track.bestDCAZ());
        }
        if (!otrack.has_collision()) {
          if constexpr (hasRecoCent<C>()) {
            binnedRegistry.fill(HIST("Tracks/Control/ExtraTracksEtaZvtx"), otrack.eta(), z, c);
            binnedRegistry.fill(HIST("Tracks/Control/ExtraTracksPhiEta"), otrack.phi(), otrack.eta(), c);
            binnedRegistry.fill(HIST("Tracks/Control/ExtraDCAXYPt"), otrack.pt(), track.bestDCAXY(), c);
            binnedRegistry.fill(HIST("Tracks/Control/ExtraDCAZPt"), otrack.pt(), track.bestDCAZ(), c);
          } else {
            inclusiveRegistry.fill(HIST("Tracks/Control/ExtraTracksEtaZvtx"), otrack.eta(), z);
            inclusiveRegistry.fill(HIST("Tracks/Control/ExtraTracksPhiEta"), otrack.phi(), otrack.eta());
            inclusiveRegistry.fill(HIST("Tracks/Control/ExtraDCAXYPt"), otrack.pt(), track.bestDCAXY());
            inclusiveRegistry.fill(HIST("Tracks/Control/ExtraDCAZPt"), otrack.pt(), track.bestDCAZ());
          }
        } else if (otrack.collisionId() != track.bestCollisionId()) {
          usedTracksIdsDF.emplace_back(track.trackId());
          if constexpr (hasRecoCent<C>()) {
            binnedRegistry.fill(HIST("Tracks/Control/ReassignedTracksEtaZvtx"), otrack.eta(), z, c);
            binnedRegistry.fill(HIST("Tracks/Control/ReassignedTracksPhiEta"), otrack.phi(), otrack.eta(), c);
            binnedRegistry.fill(HIST("Tracks/Control/ReassignedVertexCorr"), otrack.collision_as<C>().posZ(), z, c);
            binnedRegistry.fill(HIST("Tracks/Control/ReassignedDCAXYPt"), otrack.pt(), track.bestDCAXY(), c);
            binnedRegistry.fill(HIST("Tracks/Control/ReassignedDCAZPt"), otrack.pt(), track.bestDCAZ(), c);
          } else {
            inclusiveRegistry.fill(HIST("Tracks/Control/ReassignedTracksEtaZvtx"), otrack.eta(), z);
            inclusiveRegistry.fill(HIST("Tracks/Control/ReassignedTracksPhiEta"), otrack.phi(), otrack.eta());
            inclusiveRegistry.fill(HIST("Tracks/Control/ReassignedVertexCorr"), otrack.collision_as<C>().posZ(), z);
            inclusiveRegistry.fill(HIST("Tracks/Control/ReassignedDCAXYPt"), otrack.pt(), track.bestDCAXY());
            inclusiveRegistry.fill(HIST("Tracks/Control/ReassignedDCAZPt"), otrack.pt(), track.bestDCAZ());
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
        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST("Tracks/EtaZvtx"), track.eta(), z, c);
          binnedRegistry.fill(HIST("Tracks/PhiEta"), track.phi(), track.eta(), c);
          binnedRegistry.fill(HIST("Tracks/Control/PtEta"), track.pt(), track.eta(), c);
          binnedRegistry.fill(HIST("Tracks/Control/DCAXYPt"), track.pt(), track.dcaXY(), c);
          binnedRegistry.fill(HIST("Tracks/Control/DCAZPt"), track.pt(), track.dcaZ(), c);
        } else {
          inclusiveRegistry.fill(HIST("Tracks/EtaZvtx"), track.eta(), z);
          inclusiveRegistry.fill(HIST("Tracks/PhiEta"), track.phi(), track.eta());
          inclusiveRegistry.fill(HIST("Tracks/Control/PtEta"), track.pt(), track.eta());
          inclusiveRegistry.fill(HIST("Tracks/Control/DCAXYPt"), track.pt(), track.dcaXY());
          inclusiveRegistry.fill(HIST("Tracks/Control/DCAZPt"), track.pt(), track.dcaZ());
        }
      }
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST("Events/NtrkZvtx"), Ntrks, z, c);
      } else {
        if (Ntrks > 0 || groupPVContrib.size() > 0) {
          if (groupPVContrib.size() > 0) {
            inclusiveRegistry.fill(HIST("Events/Selection"), static_cast<float>(EvSelBins::kSelectedPVgt0));
          }
          if (Ntrks > 0) {
            inclusiveRegistry.fill(HIST("Events/Selection"), static_cast<float>(EvSelBins::kSelectedgt0));
          }
          for (auto& track : atracks) {
            if (Ntrks > 0) {
              inclusiveRegistry.fill(HIST("Tracks/EtaZvtx_gt0"), track.track_as<FiTracks>().eta(), z);
            }
            if (groupPVContrib.size() > 0) {
              inclusiveRegistry.fill(HIST("Tracks/EtaZvtx_PVgt0"), track.track_as<FiTracks>().eta(), z);
            }
          }
          for (auto& track : tracks) {
            if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
              continue;
            }
            if (std::find(usedTracksIdsDF.begin(), usedTracksIdsDF.end(), track.globalIndex()) != usedTracksIdsDF.end()) {
              continue;
            }
            if (Ntrks > 0) {
              inclusiveRegistry.fill(HIST("Tracks/EtaZvtx_gt0"), track.eta(), z);
            }
            if (groupPVContrib.size() > 0) {
              inclusiveRegistry.fill(HIST("Tracks/EtaZvtx_PVgt0"), track.eta(), z);
            }
          }
        }
        inclusiveRegistry.fill(HIST("Events/NtrkZvtx"), Ntrks, z);
      }
    } else {
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST("Events/Selection"), 3., c);
      } else {
        inclusiveRegistry.fill(HIST("Events/Selection"), static_cast<float>(EvSelBins::kRejected));
      }
    }
  }

  void processCounting(
    ExCols::iterator const& collision,
    FiTracks const& tracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    processCountingGeneralwAmbiguous<ExCols>(collision, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCounting, "Count tracks", false);

  void processCountingNoAmb(
    ExCols::iterator const& collision,
    FiTracks const& tracks)
  {
    processCountingGeneral<ExCols>(collision, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCountingNoAmb, "Count tracks w/o ambiguous", false);

  using ExColsCentFT0C = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>;
  void processCountingCentralityFT0C(
    ExColsCentFT0C::iterator const& collision,
    FiTracks const& tracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    processCountingGeneralwAmbiguous<ExColsCentFT0C>(collision, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCountingCentralityFT0C, "Count tracks in FT0C centrality bins", false);

  void processCountingCentralityFT0CNoAmb(
    ExColsCentFT0C::iterator const& collision,
    FiTracks const& tracks)
  {
    processCountingGeneral<ExColsCentFT0C>(collision, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCountingCentralityFT0CNoAmb, "Count tracks in FT0C centrality bins w/o ambiguous", false);

  using ExColsCentFT0M = soa::Join<aod::Collisions, aod::CentFT0Ms, aod::EvSels>;
  void processCountingCentralityFT0M(
    ExColsCentFT0M::iterator const& collision,
    FiTracks const& tracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    processCountingGeneralwAmbiguous<ExColsCentFT0M>(collision, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCountingCentralityFT0M, "Count tracks in FT0M centrality bins", false);

  void processCountingCentralityFT0MNoAmb(
    ExColsCentFT0M::iterator const& collision,
    FiTracks const& tracks)
  {
    processCountingGeneral<ExColsCentFT0M>(collision, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCountingCentralityFT0MNoAmb, "Count tracks in FT0M centrality bins w/o ambiguous", false);

  using Particles = soa::Filtered<aod::McParticles>;
  using LabeledTracksEx = soa::Join<LabeledTracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
  using FiLTracks = soa::Filtered<LabeledTracksEx>;
  using ParticlesI = soa::Filtered<soa::Join<aod::McParticles, aod::ParticlesToTracks>>;
  expressions::Filter primaries = ncheckbit(aod::mcparticle::flags, (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary);

  template <typename C, typename MC>
  void processTrackEfficiencyIndexedGeneral(
    typename soa::Join<C, aod::McCollisionLabels>::iterator const& collision,
    MC const&, ParticlesI const& particles,
    FiLTracks const& /*tracks*/)
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
      inclusiveRegistry.fill(HIST("Tracks/Control/PtGenINoEtaCut"), particle.pt());

      if (std::abs(particle.eta()) < estimatorEta) {
        inclusiveRegistry.fill(HIST("Tracks/Control/PtGenI"), particle.pt());
        if (particle.pdgCode() == speciesIds[0]) {
          inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[0]) + HIST("/PtGenI"), particle.pt());
        } else if (particle.pdgCode() == speciesIds[1]) {
          inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[1]) + HIST("/PtGenI"), particle.pt());
        } else if (particle.pdgCode() == speciesIds[2]) {
          inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[2]) + HIST("/PtGenI"), particle.pt());
        } else if (particle.pdgCode() == speciesIds[3]) {
          inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[3]) + HIST("/PtGenI"), particle.pt());
        }
      }
      if (particle.has_tracks()) {
        auto counted = false;
        auto countedNoEtaCut = false;
        auto counter = 0;
        auto relatedTracks = particle.template filtered_tracks_as<FiLTracks>();
        for (auto const& track : relatedTracks) {
          ++counter;
          if (!countedNoEtaCut) {
            inclusiveRegistry.fill(HIST("Tracks/Control/PtEfficiencyINoEtaCut"), particle.pt());
            countedNoEtaCut = true;
          }
          if (std::abs(track.eta()) < estimatorEta) {
            if (!counted) {
              inclusiveRegistry.fill(HIST("Tracks/Control/PtEfficiencyI"), particle.pt());
              if (particle.pdgCode() == speciesIds[0]) {
                inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[0]) + HIST("/PtEfficiencyI"), particle.pt());
              } else if (particle.pdgCode() == speciesIds[1]) {
                inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[1]) + HIST("/PtEfficiencyI"), particle.pt());
              } else if (particle.pdgCode() == speciesIds[2]) {
                inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[2]) + HIST("/PtEfficiencyI"), particle.pt());
              } else if (particle.pdgCode() == speciesIds[3]) {
                inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[3]) + HIST("/PtEfficiencyI"), particle.pt());
              }
              counted = true;
            }
          }
          if (counter > 1) {
            inclusiveRegistry.fill(HIST("Tracks/Control/PtEfficiencyISecondariesNoEtaCut"), particle.pt());
            if (std::abs(track.eta()) < estimatorEta) {
              inclusiveRegistry.fill(HIST("Tracks/Control/PtEfficiencyISecondaries"), particle.pt());
            }
          }
        }
        if (counter > 1) {
          for (auto const& track : relatedTracks) {
            for (auto layer = 0; layer < 7; ++layer) {
              if (track.itsClusterMap() & (uint8_t(1) << layer)) {
                inclusiveRegistry.fill(HIST("Tracks/Control/ITSClusters"), layer + 1);
              }
            }
            auto hasbit = false;
            for (auto bit = 0; bit < 16; ++bit) {
              if (track.mcMask() & (uint8_t(1) << bit)) {
                inclusiveRegistry.fill(HIST("Tracks/Control/Mask"), bit);
                hasbit = true;
              }
            }
            if (!hasbit) {
              inclusiveRegistry.fill(HIST("Tracks/Control/Mask"), 16);
            }
          }
        }
        if (relatedTracks.size() > 1) {
          inclusiveRegistry.fill(HIST("Tracks/Control/PhiEtaGenDuplicates"), particle.phi(), particle.eta());
          for (auto const& track : relatedTracks) {
            inclusiveRegistry.fill(HIST("Tracks/Control/PhiEtaDuplicates"), track.phi(), track.eta());
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
  void processTrackEfficiencyGeneralAmbiguous(
    typename soa::Join<C, aod::McCollisionLabels>::iterator const& collision,
    MC const&, Particles const& particles,
    FiLTracks const& tracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    if (useEvSel && !collision.sel8()) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }
    float c_rec = -1;
    float c_gen = -1;
    if constexpr (hasRecoCent<C>()) {
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        c_rec = collision.centFT0C();
      } else if (C::template contains<aod::CentFT0Ms>()) {
        c_rec = collision.centFT0M();
      }
    }
    auto mcCollision = collision.mcCollision();
    if constexpr (hasSimCent<MC>()) {
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
        inclusiveRegistry.fill(HIST("Tracks/Control/PtEfficiencyNoEtaCut"), particle.pt());
        if (std::abs(otrack.eta()) < estimatorEta) {
          inclusiveRegistry.fill(HIST("Tracks/Control/PtEfficiency"), particle.pt());
          if (particle.pdgCode() == speciesIds[0]) {
            inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[0]) + HIST("/PtEfficiency"), particle.pt());
          } else if (particle.pdgCode() == speciesIds[1]) {
            inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[1]) + HIST("/PtEfficiency"), particle.pt());
          } else if (particle.pdgCode() == speciesIds[2]) {
            inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[2]) + HIST("/PtEfficiency"), particle.pt());
          } else if (particle.pdgCode() == speciesIds[3]) {
            inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[3]) + HIST("/PtEfficiency"), particle.pt());
          }
        }
      } else {
        inclusiveRegistry.fill(HIST("Tracks/Control/PtEfficiencyFakes"), otrack.pt());
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
        inclusiveRegistry.fill(HIST("Tracks/Control/PtEfficiencyNoEtaCut"), particle.pt());
        if (std::abs(track.eta()) < estimatorEta) {
          inclusiveRegistry.fill(HIST("Tracks/Control/PtEfficiency"), particle.pt());
          if (particle.pdgCode() == speciesIds[0]) {
            inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[0]) + HIST("/PtEfficiency"), particle.pt());
          } else if (particle.pdgCode() == speciesIds[1]) {
            inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[1]) + HIST("/PtEfficiency"), particle.pt());
          } else if (particle.pdgCode() == speciesIds[2]) {
            inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[2]) + HIST("/PtEfficiency"), particle.pt());
          } else if (particle.pdgCode() == speciesIds[3]) {
            inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[3]) + HIST("/PtEfficiency"), particle.pt());
          }
        }
      } else {
        inclusiveRegistry.fill(HIST("Tracks/Control/PtEfficiencyFakes"), track.pt());
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
      inclusiveRegistry.fill(HIST("Tracks/Control/PtGenNoEtaCut"), particle.pt());
      if (std::abs(particle.eta()) < estimatorEta) {
        inclusiveRegistry.fill(HIST("Tracks/Control/PtGen"), particle.pt());
        if (particle.pdgCode() == speciesIds[0]) {
          inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[0]) + HIST("/PtGen"), particle.pt());
        } else if (particle.pdgCode() == speciesIds[1]) {
          inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[1]) + HIST("/PtGen"), particle.pt());
        } else if (particle.pdgCode() == speciesIds[2]) {
          inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[2]) + HIST("/PtGen"), particle.pt());
        } else if (particle.pdgCode() == speciesIds[3]) {
          inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[3]) + HIST("/PtGen"), particle.pt());
        }
      }
    }
  }

  template <typename C, typename MC>
  void processTrackEfficiencyGeneral(
    typename soa::Join<C, aod::McCollisionLabels>::iterator const& collision,
    MC const&, Particles const& particles,
    FiLTracks const& tracks)
  {
    if (useEvSel && !collision.sel8()) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }
    float c_rec = -1;
    float c_gen = -1;
    if constexpr (hasRecoCent<C>()) {
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        c_rec = collision.centFT0C();
      } else if (C::template contains<aod::CentFT0Ms>()) {
        c_rec = collision.centFT0M();
      }
    }
    auto mcCollision = collision.mcCollision();
    if constexpr (hasSimCent<MC>()) {
      c_gen = mcCollision.centrality();
    }

    auto particlesPerCol = particles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);

    usedTracksIds.clear();
    for (auto const& track : tracks) {
      if (track.has_mcParticle()) {
        auto particle = track.template mcParticle_as<Particles>();
        inclusiveRegistry.fill(HIST("Tracks/Control/PtEfficiencyNoEtaCut"), particle.pt());
        if (std::abs(track.eta()) < estimatorEta) {
          inclusiveRegistry.fill(HIST("Tracks/Control/PtEfficiency"), particle.pt());
          if (particle.pdgCode() == speciesIds[0]) {
            inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[0]) + HIST("/PtEfficiency"), particle.pt());
          } else if (particle.pdgCode() == speciesIds[1]) {
            inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[1]) + HIST("/PtEfficiency"), particle.pt());
          } else if (particle.pdgCode() == speciesIds[2]) {
            inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[2]) + HIST("/PtEfficiency"), particle.pt());
          } else if (particle.pdgCode() == speciesIds[3]) {
            inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[3]) + HIST("/PtEfficiency"), particle.pt());
          }
        }
      } else {
        inclusiveRegistry.fill(HIST("Tracks/Control/PtEfficiencyFakes"), track.pt());
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
      inclusiveRegistry.fill(HIST("Tracks/Control/PtGenNoEtaCut"), particle.pt());
      if (std::abs(particle.eta()) < estimatorEta) {
        inclusiveRegistry.fill(HIST("Tracks/Control/PtGen"), particle.pt());
        if (particle.pdgCode() == speciesIds[0]) {
          inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[0]) + HIST("/PtGen"), particle.pt());
        } else if (particle.pdgCode() == speciesIds[1]) {
          inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[1]) + HIST("/PtGen"), particle.pt());
        } else if (particle.pdgCode() == speciesIds[2]) {
          inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[2]) + HIST("/PtGen"), particle.pt());
        } else if (particle.pdgCode() == speciesIds[3]) {
          inclusiveRegistry.fill(HIST("Tracks/Control/") + HIST(species[3]) + HIST("/PtGen"), particle.pt());
        }
      }
    }
  }

  void processTrackEfficiency(
    soa::Join<ExCols, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, Particles const& mcParticles,
    FiLTracks const& filtracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    processTrackEfficiencyGeneralAmbiguous<ExCols, aod::McCollisions>(collision, mccollisions, mcParticles, filtracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processTrackEfficiency, "Calculate tracking efficiency vs pt", false);

  void processTrackEfficiencyNoAmb(
    soa::Join<ExCols, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, Particles const& mcParticles,
    FiLTracks const& filtracks)
  {
    processTrackEfficiencyGeneral<ExCols, aod::McCollisions>(collision, mccollisions, mcParticles, filtracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processTrackEfficiencyNoAmb, "Calculate tracking efficiency vs pt w/o ambiguous", false);

  template <typename MC, typename C>
  void processGenGeneral(
    typename MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const* atracks)
  {
    float c_gen = -1;
    // add generated centrality estimation
    if constexpr (hasSimCent<MC>()) {
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
    if constexpr (hasRecoCent<C>()) {
      binnedRegistry.fill(HIST("Events/NtrkZvtxGen_t"), nCharged, mcCollision.posZ(), c_gen);
      binnedRegistry.fill(HIST("Events/Efficiency"), static_cast<float>(EvEffBins::kGen), c_gen);
    } else {
      inclusiveRegistry.fill(HIST("Events/NtrkZvtxGen_t"), nCharged, mcCollision.posZ());
      inclusiveRegistry.fill(HIST("Events/Efficiency"), static_cast<float>(EvEffBins::kGen));
    }

    if (nCharged > 0) {
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST("Events/Efficiency"), static_cast<float>(EvEffBins::kGengt0), c_gen);
      } else {
        inclusiveRegistry.fill(HIST("Events/Efficiency"), static_cast<float>(EvEffBins::kGengt0));
      }
    }
    bool atLeastOne = false;
    bool atLeastOne_gt0 = false;
    bool atLeastOne_PVgt0 = false;
    auto moreThanOne = 0;
    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());

    auto Nrec = 0;
    std::vector<int> NrecPerCol;
    std::vector<float> c_recPerCol;
    std::vector<int> NPVPerCol;
    std::vector<float> NFT0APerCol;
    std::vector<float> NFT0CPerCol;
    std::vector<float> NFDDAPerCol;
    std::vector<float> NFDDCPerCol;

    for (auto& collision : collisions) {
      usedTracksIds.clear();
      float c_rec = -1;
      if constexpr (hasRecoCent<C>()) {
        if constexpr (C::template contains<aod::CentFT0Cs>()) {
          c_rec = collision.centFT0C();
        } else if (C::template contains<aod::CentFT0Ms>()) {
          c_rec = collision.centFT0M();
        }
        c_recPerCol.emplace_back(c_rec);
        binnedRegistry.fill(HIST("Events/Efficiency"), static_cast<float>(EvEffBins::kRec), c_gen);
      } else {
        inclusiveRegistry.fill(HIST("Events/Efficiency"), static_cast<float>(EvEffBins::kRec));
      }
      if (!useEvSel || collision.sel8()) {
        Nrec = 0;
        ++moreThanOne;
        atLeastOne = true;

        auto groupPVcontrib = pvContribTracksIUEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
        if (groupPVcontrib.size() > 0) {
          if constexpr (hasRecoCent<C>()) {
            binnedRegistry.fill(HIST("Events/Efficiency"), static_cast<float>(EvEffBins::kSelectedPVgt0), c_gen);
          } else {
            inclusiveRegistry.fill(HIST("Events/Efficiency"), static_cast<float>(EvEffBins::kSelectedPVgt0));
          }
        }

        if (atracks != nullptr) {
          auto perCollisionASample = atracks->sliceBy(perColU, collision.globalIndex());
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
        }
        auto perCollisionSample = tracks.sliceBy(perCol, collision.globalIndex());
        for (auto const& track : perCollisionSample) {
          if (atracks != nullptr) {
            if (std::find(usedTracksIds.begin(), usedTracksIds.end(), track.globalIndex()) != usedTracksIds.end()) {
              continue;
            }
            if (std::find(usedTracksIdsDFMC.begin(), usedTracksIdsDFMC.end(), track.globalIndex()) != usedTracksIdsDFMC.end()) {
              continue;
            }
          }
          if (std::abs(track.eta()) < estimatorEta) {
            ++Nrec;
          }
        }
        NrecPerCol.emplace_back(Nrec);
        NPVPerCol.emplace_back(collision.numContrib());
        if (collision.has_foundFT0()) {
          auto ft0 = collision.foundFT0();
          float tA = 0;
          float tC = 0;
          for (auto i = 0u; i < ft0.amplitudeA().size(); ++i) {
            tA += ft0.amplitudeA()[i];
          }
          for (auto i = 0u; i < ft0.amplitudeC().size(); ++i) {
            tC += ft0.amplitudeC()[i];
          }
          NFT0APerCol.emplace_back(tA);
          NFT0CPerCol.emplace_back(tC);
        } else {
          NFT0APerCol.emplace_back(-1);
          NFT0CPerCol.emplace_back(-1);
        }
        if (collision.has_foundFDD()) {
          auto fdd = collision.foundFDD();
          float tA = 0;
          float tC = 0;
          for (auto i = 0u; i < 8; ++i) {
            tA += fdd.chargeA()[i];
          }
          for (auto i = 0u; i < 8; ++i) {
            tC += fdd.chargeC()[i];
          }
          NFDDAPerCol.emplace_back(tA);
          NFDDCPerCol.emplace_back(tC);
        } else {
          NFDDAPerCol.emplace_back(-1);
          NFDDCPerCol.emplace_back(-1);
        }

        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST("Events/Efficiency"), static_cast<float>(EvEffBins::kSelected), c_gen);
        } else {
          inclusiveRegistry.fill(HIST("Events/Efficiency"), static_cast<float>(EvEffBins::kSelected));
        }

        if (Nrec > 0) {
          if constexpr (hasRecoCent<C>()) {
            binnedRegistry.fill(HIST("Events/Efficiency"), static_cast<float>(EvEffBins::kSelectedgt0), c_gen);
          } else {
            inclusiveRegistry.fill(HIST("Events/Efficiency"), static_cast<float>(EvEffBins::kSelectedgt0));
          }
          atLeastOne_gt0 = true;
        }
        if (groupPVcontrib.size() > 0) {
          if constexpr (hasRecoCent<C>()) {
            binnedRegistry.fill(HIST("Events/Efficiency"), static_cast<float>(EvEffBins::kSelectedPVgt0), c_gen);
          } else {
            inclusiveRegistry.fill(HIST("Events/Efficiency"), static_cast<float>(EvEffBins::kSelectedPVgt0));
          }
          atLeastOne_PVgt0 = true;
        }

        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST("Events/NtrkZvtxGen"), Nrec, collision.posZ(), c_rec);
        } else {
          inclusiveRegistry.fill(HIST("Events/NtrkZvtxGen"), Nrec, collision.posZ());
        }
      }
    }

    if (fillResponse) {
      for (auto i = 0U; i < NrecPerCol.size(); ++i) {
        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST("Events/Response"), NrecPerCol[i], nCharged, mcCollision.posZ(), c_recPerCol[i]);
          binnedRegistry.fill(HIST("Events/EfficiencyMult"), nCharged, mcCollision.posZ(), c_recPerCol[i]);
          if (responseStudy) {
            binnedRegistry.fill(HIST("Events/Control/MultiResponse"), nCharged, NrecPerCol[i], NPVPerCol[i], NFT0APerCol[i], NFT0CPerCol[i], NFDDAPerCol[i], NFDDCPerCol[i], mcCollision.posZ(), c_recPerCol[i]);
          }
        } else {
          inclusiveRegistry.fill(HIST("Events/Response"), NrecPerCol[i], nCharged, mcCollision.posZ());
          inclusiveRegistry.fill(HIST("Events/EfficiencyMult"), nCharged, mcCollision.posZ());
          if (responseStudy) {
            inclusiveRegistry.fill(HIST("Events/Control/MultiResponse"), nCharged, NrecPerCol[i], NPVPerCol[i], NFT0APerCol[i], NFT0CPerCol[i], NFDDAPerCol[i], NFDDCPerCol[i], mcCollision.posZ());
          }
        }
      }
      if (moreThanOne > 1) {
        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST("Events/SplitMult"), nCharged, mcCollision.posZ(), c_gen);
        } else {
          inclusiveRegistry.fill(HIST("Events/SplitMult"), nCharged, mcCollision.posZ());
        }
      }
    }

    if (collisions.size() == 0) {
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST("Events/NotFoundEventZvtx"), mcCollision.posZ(), c_gen);
      } else {
        inclusiveRegistry.fill(HIST("Events/NotFoundEventZvtx"), mcCollision.posZ());
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
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST("Tracks/EtaZvtxGen_t"), particle.eta(), mcCollision.posZ(), c_gen);
        binnedRegistry.fill(HIST("Tracks/Control/PtEtaGen"), particle.pt(), particle.eta(), c_gen);
      } else {
        inclusiveRegistry.fill(HIST("Tracks/EtaZvtxGen_t"), particle.eta(), mcCollision.posZ());
        inclusiveRegistry.fill(HIST("Tracks/Control/PtEtaGen"), particle.pt(), particle.eta());
      }
      if (nCharged > 0) {
        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST("Tracks/EtaZvtxGen_gt0t"), particle.eta(), mcCollision.posZ(), c_gen);
        } else {
          inclusiveRegistry.fill(HIST("Tracks/EtaZvtxGen_gt0t"), particle.eta(), mcCollision.posZ());
        }
      }
      if (atLeastOne) {
        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST("Tracks/EtaZvtxGen"), particle.eta(), mcCollision.posZ(), c_gen);
          if (atLeastOne_gt0) {
            binnedRegistry.fill(HIST("Tracks/EtaZvtxGen_gt0"), particle.eta(), mcCollision.posZ(), c_gen);
          }
          if (atLeastOne_PVgt0) {
            binnedRegistry.fill(HIST("Tracks/EtaZvtxGen_PVgt0"), particle.eta(), mcCollision.posZ(), c_gen);
          }
          binnedRegistry.fill(HIST("Tracks/PhiEtaGen"), particle.phi(), particle.eta(), c_gen);
        } else {
          inclusiveRegistry.fill(HIST("Tracks/EtaZvtxGen"), particle.eta(), mcCollision.posZ());
          if (atLeastOne_gt0) {
            inclusiveRegistry.fill(HIST("Tracks/EtaZvtxGen_gt0"), particle.eta(), mcCollision.posZ());
          }
          if (atLeastOne_PVgt0) {
            inclusiveRegistry.fill(HIST("Tracks/EtaZvtxGen_PVgt0"), particle.eta(), mcCollision.posZ());
          }
          inclusiveRegistry.fill(HIST("Tracks/PhiEtaGen"), particle.phi(), particle.eta());
        }
      }
    }
  }

  using MC = aod::McCollisions; // soa::Join<aod::McCollisions, aod::HepMCXSections>;
  void processGen(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExCols, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MC, ExCols>(mcCollision, collisions, particles, tracks, &atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGen, "Process generator-level info", false);

  void processGenNoAmb(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExCols, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MC, ExCols>(mcCollision, collisions, particles, tracks, nullptr);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenNoAmb, "Process generator-level info w/o ambiguous", false);

  void processGenFT0C(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MC, ExColsCentFT0C>(mcCollision, collisions, particles, tracks, &atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0C, "Process generator-level info (FT0C centrality)", false);

  void processGenFT0CNoAmb(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MC, ExColsCentFT0C>(mcCollision, collisions, particles, tracks, nullptr);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0CNoAmb, "Process generator-level info (FT0C centrality) w/o ambiguous", false);

  void processGenFT0M(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MC, ExColsCentFT0M>(mcCollision, collisions, particles, tracks, &atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0M, "Process generator-level info (FT0M centrality)", false);

  void processGenFT0MNoAmb(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MC, ExColsCentFT0M>(mcCollision, collisions, particles, tracks, nullptr);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0MNoAmb, "Process generator-level info (FT0M centrality) w/o ambiguous", false);

  using MChi = soa::Join<aod::McCollisions, aod::HepMCHeavyIons>;

  void processGenFT0Chi(
    MChi::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MChi, ExColsCentFT0C>(mcCollision, collisions, particles, tracks, &atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0Chi, "Process generator-level info (FT0C centrality, HI)", false);

  void processGenFT0ChiNoAmb(
    MChi::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MChi, ExColsCentFT0C>(mcCollision, collisions, particles, tracks, nullptr);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0ChiNoAmb, "Process generator-level info (FT0C centrality, HI) w/o ambiguous", false);

  void processGenFT0Mhi(
    MChi::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MChi, ExColsCentFT0M>(mcCollision, collisions, particles, tracks, &atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0Mhi, "Process generator-level info (FT0M centrality, HI)", false);

  void processGenFT0MhiNoAmb(
    MChi::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MChi, ExColsCentFT0M>(mcCollision, collisions, particles, tracks, nullptr);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0MhiNoAmb, "Process generator-level info (FT0M centrality, HI) w/o ambiguous", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityCounter>(cfgc)};
}
