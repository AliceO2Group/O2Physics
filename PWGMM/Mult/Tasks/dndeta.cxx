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
#include <TDatabasePDG.h>
#include <TPDGCode.h>

#include "Gencentralities.h"
#include "Index.h"
#include "bestCollisionTable.h"

#include "Axes.h"
#include "Functions.h"
#include "Selections.h"
#include "Histograms.h"

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
  ConfigurableAxis multBinning{"multBinning", {301, -0.5, 300.5}, "Multiplicity axis binning"};
  ConfigurableAxis centBinning{"centBinning", {VARIABLE_WIDTH, 0, 10, 20, 30, 40, 50, 60, 70, 80, 100}, "Centrality axis binning"};
  ConfigurableAxis occuBinning{"occuBinning", {VARIABLE_WIDTH, 0, 500, 1000, 2000, 5000, 10000}, "Occupancy axis binning"}; // Pb-Pb default

  Configurable<bool> fillResponse{"fillResponse", true, "Fill response matrix"};
  Configurable<bool> useProcId{"use-process-id", true, "Use process ID from generator"};
  Configurable<bool> addFT0{"addFT0", false, "add FT0 estimators"};
  Configurable<bool> addFDD{"addFDD", false, "add FDD estimators"};

  Configurable<bool> useEvSel{"useEvSel", true, "Use event selection"};
  Configurable<bool> excludeTFborder{"excludeTFborder", true, "Exclude TF border"};
  Configurable<bool> excludeITSROFborder{"excludeITSROFborder", true, "Exclude ITS RO frame border"};
  Configurable<bool> requireFT0PVcoincidence{"requireFT0PVcoincidence", true, "Require coincidence between FT0 and PV"};
  Configurable<bool> rejectITSonly{"rejectITSonly", true, "Reject ITS-only vertex"};
  Configurable<bool> requireVtxTOFMatched{"checkTOFMatch", true, "Consider only vertex with TOF match"};

  template <typename C>
  inline bool isCollisionSelected(C const& collision)
  {
    return collision.selection_bit(aod::evsel::kIsTriggerTVX) &&
           (!excludeTFborder || collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) &&
           (!excludeITSROFborder || collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) &&
           (!requireFT0PVcoincidence || collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) &&
           (!rejectITSonly || collision.selection_bit(aod::evsel::kIsVertexITSTPC)) &&
           (!requireVtxTOFMatched || collision.selection_bit(aod::evsel::kIsVertexTOFmatched));
  }

  template <typename B>
  inline bool isBCSelected(B const& bc)
  {
    return bc.selection_bit(aod::evsel::kIsTriggerTVX) &&
           (!excludeTFborder || bc.selection_bit(aod::evsel::kNoTimeFrameBorder)) &&
           (!excludeITSROFborder || bc.selection_bit(aod::evsel::kNoITSROFrameBorder));
  }

  HistogramRegistry commonRegistry{
    "Common",
    {
      {BCSelection.data(), ";status;count", {HistType::kTH1F, {{3, 0.5, 3.5}}}} //
    },                                                                          //
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true //
  };

  HistogramRegistry inclusiveRegistry{
    InclusivePrefix.data(),
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true};
  HistogramRegistry binnedRegistry{
    BinnedPrefix.data(),
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    false,
    true};

  std::vector<int> usedTracksIds;
  std::vector<int> usedTracksIdsDF;
  std::vector<int> usedTracksIdsDFMC;
  std::vector<int> usedTracksIdsDFMCEff;

  void init(InitContext&)
  {
    AxisSpec MultAxis = {multBinning};
    AxisSpec CentAxis = {centBinning, "centrality"};
    AxisSpec OccuAxis = {occuBinning, "occupancy"};
    {
      auto hstat = commonRegistry.get<TH1>(HIST(BCSelection));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "Good BCs");
      x->SetBinLabel(2, "BCs with collisions");
      x->SetBinLabel(3, "BCs with pile-up/splitting");
    }

    if (doprocessEventStat) {
      inclusiveRegistry.add({EventChi2.data(), " ; #chi^2", {HistType::kTH2F, {{101, -0.1, 10.1}, OccuAxis}}});
      inclusiveRegistry.add({EventTimeRes.data(), " ; t (ms)", {HistType::kTH2F, {{1001, -0.1, 100.1}, OccuAxis}}});
    }
    if (doprocessEventStatCentralityFT0C || doprocessEventStatCentralityFT0M) {
      binnedRegistry.add({EventChi2.data(), " ; #chi^2; centrality", {HistType::kTHnSparseF, {{101, -0.1, 10.1}, CentAxis, OccuAxis}}});
      binnedRegistry.add({EventTimeRes.data(), " ; t (ms); centrality", {HistType::kTHnSparseF, {{1001, -0.1, 100.1}, CentAxis, OccuAxis}}});
    }

    if (doprocessCountingAmbiguous || doprocessCounting) {
      inclusiveRegistry.add({EventSelection.data(), ";status;occupancy;events", {HistType::kTH2F, {{static_cast<int>(EvSelBins::kRejected), 0.5, static_cast<float>(EvSelBins::kRejected) + 0.5}, OccuAxis}}});
      auto hstat = inclusiveRegistry.get<TH2>(HIST(EventSelection));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(static_cast<int>(EvSelBins::kAll), EvSelBinLabels[static_cast<int>(EvSelBins::kAll)].data());
      x->SetBinLabel(static_cast<int>(EvSelBins::kSelected), EvSelBinLabels[static_cast<int>(EvSelBins::kSelected)].data());
      x->SetBinLabel(static_cast<int>(EvSelBins::kSelectedgt0), EvSelBinLabels[static_cast<int>(EvSelBins::kSelectedgt0)].data());
      x->SetBinLabel(static_cast<int>(EvSelBins::kSelectedPVgt0), EvSelBinLabels[static_cast<int>(EvSelBins::kSelectedPVgt0)].data());
      x->SetBinLabel(static_cast<int>(EvSelBins::kRejected), EvSelBinLabels[static_cast<int>(EvSelBins::kRejected)].data());

      inclusiveRegistry.add({NtrkZvtx.data(), "; N_{trk}; Z_{vtx} (cm);occupancy; events", {HistType::kTHnSparseF, {MultAxis, ZAxis, OccuAxis}}});
      inclusiveRegistry.add({NpvcZvtx.data(), "; N_{PVc}; Z_{vtx} (cm);occupancy; events", {HistType::kTHnSparseF, {MultAxis, ZAxis, OccuAxis}}});
      inclusiveRegistry.add({EtaZvtx.data(), "; #eta; Z_{vtx} (cm);occupancy; tracks", {HistType::kTHnSparseF, {EtaAxis, ZAxis, OccuAxis}}});
      inclusiveRegistry.add({EtaZvtx_gt0.data(), "; #eta; Z_{vtx} (cm);occupancy; tracks", {HistType::kTHnSparseF, {EtaAxis, ZAxis, OccuAxis}}});
      inclusiveRegistry.add({EtaZvtx_PVgt0.data(), "; #eta; Z_{vtx} (cm);occupancy; tracks", {HistType::kTHnSparseF, {EtaAxis, ZAxis, OccuAxis}}});
      inclusiveRegistry.add({PhiEta.data(), "; #varphi; #eta;occupancy; tracks", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, OccuAxis}}});
      inclusiveRegistry.add({PtEta.data(), " ; p_{T} (GeV/c);occupancy; #eta", {HistType::kTHnSparseF, {PtAxis, EtaAxis, OccuAxis}}});
      inclusiveRegistry.add({DCAXYPt.data(), " ; p_{T} (GeV/c) ; DCA_{XY} (cm);occupancy", {HistType::kTHnSparseF, {PtAxis, DCAAxis, OccuAxis}}});
      inclusiveRegistry.add({DCAZPt.data(), " ; p_{T} (GeV/c) ; DCA_{Z} (cm);occupancy", {HistType::kTHnSparseF, {PtAxis, DCAAxis, OccuAxis}}});
      if (doprocessCountingAmbiguous) {
        inclusiveRegistry.add({ReassignedDCAXYPt.data(), " ; p_{T} (GeV/c) ; DCA_{XY} (cm);occupancy", {HistType::kTHnSparseF, {PtAxis, DCAAxis, OccuAxis}}});
        inclusiveRegistry.add({ReassignedDCAZPt.data(), " ; p_{T} (GeV/c) ; DCA_{Z} (cm);occupancy", {HistType::kTHnSparseF, {PtAxis, DCAAxis, OccuAxis}}});
        inclusiveRegistry.add({ExtraDCAXYPt.data(), " ; p_{T} (GeV/c) ; DCA_{XY} (cm);occupancy", {HistType::kTHnSparseF, {PtAxis, DCAAxis, OccuAxis}}});
        inclusiveRegistry.add({ExtraDCAZPt.data(), " ; p_{T} (GeV/c) ; DCA_{Z} (cm);occupancy", {HistType::kTHnSparseF, {PtAxis, DCAAxis, OccuAxis}}});
        inclusiveRegistry.add({ExtraEtaZvtx.data(), "; #eta; Z_{vtx} (cm);occupancy; tracks", {HistType::kTHnSparseF, {EtaAxis, ZAxis, OccuAxis}}});
        inclusiveRegistry.add({ExtraPhiEta.data(), "; #varphi; #eta;occupancy; tracks", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, OccuAxis}}});
        inclusiveRegistry.add({ReassignedEtaZvtx.data(), "; #eta; Z_{vtx} (cm);occupancy; tracks", {HistType::kTHnSparseF, {EtaAxis, ZAxis, OccuAxis}}});
        inclusiveRegistry.add({ReassignedPhiEta.data(), "; #varphi; #eta;occupancy; tracks", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, OccuAxis}}});
        inclusiveRegistry.add({ReassignedZvtxCorr.data(), "; Z_{vtx}^{orig} (cm); Z_{vtx}^{re} (cm);occupancy", {HistType::kTHnSparseF, {ZAxis, ZAxis, OccuAxis}}});
      }
    }

    if (doprocessCountingAmbiguousCentralityFT0C || doprocessCountingAmbiguousCentralityFT0M || doprocessCountingCentralityFT0C || doprocessCountingCentralityFT0M) {
      binnedRegistry.add({EventSelection.data(), ";status;centrality;occupancy;events", {HistType::kTHnSparseF, {{static_cast<int>(EvSelBins::kRejected), 0.5, static_cast<float>(EvSelBins::kRejected) + 0.5}, CentAxis, OccuAxis}}});
      auto hstat = binnedRegistry.get<THnSparse>(HIST(EventSelection));
      auto* x = hstat->GetAxis(0);
      x->SetBinLabel(static_cast<int>(EvSelBins::kAll), EvSelBinLabels[static_cast<int>(EvSelBins::kAll)].data());
      x->SetBinLabel(static_cast<int>(EvSelBins::kSelected), EvSelBinLabels[static_cast<int>(EvSelBins::kSelected)].data());
      x->SetBinLabel(static_cast<int>(EvSelBins::kSelectedgt0), EvSelBinLabels[static_cast<int>(EvSelBins::kSelectedgt0)].data());
      x->SetBinLabel(static_cast<int>(EvSelBins::kSelectedPVgt0), EvSelBinLabels[static_cast<int>(EvSelBins::kSelectedPVgt0)].data());
      x->SetBinLabel(static_cast<int>(EvSelBins::kRejected), EvSelBinLabels[static_cast<int>(EvSelBins::kRejected)].data());

      binnedRegistry.add({NtrkZvtx.data(), "; N_{trk}; Z_{vtx} (cm); centrality;occupancy", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis, OccuAxis}}});
      binnedRegistry.add({NpvcZvtx.data(), "; N_{PVc}; Z_{vtx} (cm); centrality;occupancy", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis, OccuAxis}}});
      binnedRegistry.add({EtaZvtx.data(), "; #eta; Z_{vtx} (cm); centrality;occupancy", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis, OccuAxis}}});
      binnedRegistry.add({EtaZvtx_gt0.data(), "; #eta; Z_{vtx} (cm); centrality;occupancy", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis, OccuAxis}}});
      binnedRegistry.add({EtaZvtx_PVgt0.data(), "; #eta; Z_{vtx} (cm); centrality;occupancy", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis, OccuAxis}}});
      binnedRegistry.add({PhiEta.data(), "; #varphi; #eta; centrality;occupancy", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis, OccuAxis}}});
      binnedRegistry.add({PtEta.data(), " ; p_{T} (GeV/c); #eta; centrality;occupancy", {HistType::kTHnSparseF, {PtAxis, EtaAxis, CentAxis, OccuAxis}}});
      binnedRegistry.add({DCAXYPt.data(), " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality;occupancy", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis, OccuAxis}}});
      binnedRegistry.add({DCAZPt.data(), " ; p_{T} (GeV/c) ; DCA_{Z} (cm); centrality;occupancy", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis, OccuAxis}}});
      if (doprocessCountingAmbiguousCentralityFT0C || doprocessCountingAmbiguousCentralityFT0M) {
        binnedRegistry.add({ReassignedDCAXYPt.data(), " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality;occupancy", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis, OccuAxis}}});
        binnedRegistry.add({ReassignedDCAZPt.data(), " ; p_{T} (GeV/c) ; DCA_{Z} (cm); centrality;occupancy", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis, OccuAxis}}});
        binnedRegistry.add({ExtraDCAXYPt.data(), " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality;occupancy", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis, OccuAxis}}});
        binnedRegistry.add({ExtraDCAZPt.data(), " ; p_{T} (GeV/c) ; DCA_{Z} (cm); centrality;occupancy", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis, OccuAxis}}});
        binnedRegistry.add({ExtraEtaZvtx.data(), "; #eta; Z_{vtx} (cm); centrality;occupancy", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis, OccuAxis}}});
        binnedRegistry.add({ExtraPhiEta.data(), "; #varphi; #eta; centrality;occupancy", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis, OccuAxis}}});
        binnedRegistry.add({ReassignedEtaZvtx.data(), "; #eta; Z_{vtx} (cm); centrality;occupancy", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis, OccuAxis}}});
        binnedRegistry.add({ReassignedPhiEta.data(), "; #varphi; #eta; centrality;occupancy", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis, OccuAxis}}});
        binnedRegistry.add({ReassignedZvtxCorr.data(), "; Z_{vtx}^{orig} (cm); Z_{vtx}^{re} (cm); centrality;occupancy", {HistType::kTHnSparseF, {ZAxis, ZAxis, CentAxis, OccuAxis}}});
      }
    }

    if (doprocessGenAmbiguous || doprocessGen || doprocessGenAmbiguousEx || doprocessGenEx) {
      std::string effLabels{" ; N_{gen}; Z_{vtx} (cm)"};
      std::vector<AxisSpec> effAxes{MultAxis, ZAxis};
      if ((doprocessGenAmbiguousEx || doprocessGenEx) && useProcId) {
        effLabels += " ; process ID";
        effAxes.push_back(ProcAxis);
      }
      inclusiveRegistry.add({NtrkZvtxGen.data(), "; N_{trk}; Z_{vtx} (cm); occupancy; events", {HistType::kTHnSparseF, {MultAxis, ZAxis, OccuAxis}}});
      inclusiveRegistry.add({NtrkZvtxGen_t.data(), effLabels.c_str(), {HistType::kTHnSparseF, effAxes}});
      inclusiveRegistry.add({NpvcZvxtGen.data(), "; N_{PVc}; Z_{vtx} (cm); occupancy; events", {HistType::kTHnSparseF, {MultAxis, ZAxis, OccuAxis}}});
      inclusiveRegistry.add({EtaZvtxGen.data(), "; #eta; Z_{vtx} (cm); occupancy; tracks", {HistType::kTHnSparseF, {EtaAxis, ZAxis, OccuAxis}}});
      inclusiveRegistry.add({EtaZvtxGen_t.data(), "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({EtaZvtxGen_gt0.data(), "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({EtaZvtxGen_PVgt0.data(), "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({EtaZvtxGen_gt0t.data(), "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({PtEtaGen.data(), " ; p_{T} (GeV/c) ; #eta", {HistType::kTH2F, {PtAxis, EtaAxis}}});
      inclusiveRegistry.add({PhiEtaGen.data(), "; #varphi; #eta; tracks", {HistType::kTHnSparseF, {PhiAxis, EtaAxis}}});

      inclusiveRegistry.add({Efficiency.data(), "; status; occupancy; events", {HistType::kTH2F, {{static_cast<int>(EvEffBins::kSelectedPVgt0), 0.5, static_cast<float>(EvEffBins::kSelectedPVgt0) + 0.5}, OccuAxis}}});
      inclusiveRegistry.add({NotFoundZvtx.data(), " ; Z_{vtx} (cm) ", {HistType::kTH1F, {ZAxis}}});

      if (fillResponse) {
        inclusiveRegistry.add({EfficiencyMult.data(), effLabels.c_str(), {HistType::kTHnSparseF, effAxes}});
        inclusiveRegistry.add({SplitMult.data(), " ; N_{gen} ; Z_{vtx} (cm)", {HistType::kTH2F, {MultAxis, ZAxis}}});

        std::string reLabels{" ; N_{rec}; N_{PV cont}; N_{gen}"};
        std::vector<AxisSpec> reAxes{MultAxis, MultAxis, MultAxis};
        if (addFT0) {
          reLabels += " ; N_{FT0A}; N_{FT0C}";
          reAxes.push_back(FT0AAxis);
          reAxes.push_back(FT0CAxis);
        }
        if (addFDD) {
          reLabels += " ; N_{FDA}; N_{FDC}";
          reAxes.push_back(FDAAxis);
          reAxes.push_back(FDCAxis);
        }
        reLabels += " ; Z_{vtx} (cm); occupancy";
        reAxes.push_back(ZAxis);
        reAxes.push_back(OccuAxis);
        inclusiveRegistry.add({Response.data(), reLabels.c_str(), {HistType::kTHnSparseF, reAxes}});
      }

      auto heff = inclusiveRegistry.get<TH2>(HIST(Efficiency));
      auto* x = heff->GetXaxis();
      x->SetBinLabel(static_cast<int>(EvEffBins::kGen), EvEffBinLabels[static_cast<int>(EvEffBins::kGen)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kGengt0), EvEffBinLabels[static_cast<int>(EvEffBins::kGengt0)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kRec), EvEffBinLabels[static_cast<int>(EvEffBins::kRec)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelected), EvEffBinLabels[static_cast<int>(EvEffBins::kSelected)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelectedgt0), EvEffBinLabels[static_cast<int>(EvEffBins::kSelectedgt0)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelectedPVgt0), EvEffBinLabels[static_cast<int>(EvEffBins::kSelectedPVgt0)].data());
    }

    if (doprocessGenAmbiguousFT0C || doprocessGenAmbiguousFT0M || doprocessGenAmbiguousFT0Cplus || doprocessGenAmbiguousFT0Mplus || doprocessGenAmbiguousFT0Chi || doprocessGenAmbiguousFT0Mhi ||
        doprocessGenFT0C || doprocessGenFT0M || doprocessGenFT0Cplus || doprocessGenFT0Mplus || doprocessGenFT0Chi || doprocessGenFT0Mhi || doprocessGenAmbiguousExFT0C || doprocessGenAmbiguousExFT0M ||
        doprocessGenExFT0C || doprocessGenExFT0M) {
      std::string effLabels{" ; N_{gen}; Z_{vtx} (cm); centrality"};
      std::vector<AxisSpec> effAxes{MultAxis, ZAxis, CentAxis};
      if ((doprocessGenAmbiguousExFT0C || doprocessGenAmbiguousExFT0M || doprocessGenExFT0C || doprocessGenExFT0M) && useProcId) {
        effLabels += " ; process ID";
        effAxes.push_back(ProcAxis);
      }
      binnedRegistry.add({NtrkZvtxGen.data(), "; N_{trk}; Z_{vtx} (cm); centrality; occupance", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis, OccuAxis}}});
      binnedRegistry.add({NtrkZvtxGen_t.data(), effLabels.c_str(), {HistType::kTHnSparseF, effAxes}});
      binnedRegistry.add({NpvcZvxtGen.data(), "; N_{PVc}; Z_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis, OccuAxis}}});
      binnedRegistry.add({EtaZvtxGen.data(), "; #eta; Z_{vtx} (cm); centrality; occupancy", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis, OccuAxis}}});
      binnedRegistry.add({EtaZvtxGen_t.data(), "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({EtaZvtxGen_gt0.data(), "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({EtaZvtxGen_PVgt0.data(), "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({EtaZvtxGen_gt0t.data(), "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({PtEtaGen.data(), " ; p_{T} (GeV/c) ; #eta; centrality", {HistType::kTHnSparseF, {PtAxis, EtaAxis, CentAxis}}});
      binnedRegistry.add({PhiEtaGen.data(), "; #varphi; #eta; tracks", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});

      binnedRegistry.add({Efficiency.data(), "; status; centrality; occupancy; events", {HistType::kTHnSparseF, {{static_cast<int>(EvEffBins::kSelectedPVgt0), 0.5, static_cast<float>(EvEffBins::kSelectedPVgt0) + 0.5}, CentAxis, OccuAxis}}});
      binnedRegistry.add({NotFoundZvtx.data(), " ; Z_{vtx} (cm); centrality; events", {HistType::kTH2F, {ZAxis, CentAxis}}});

      if (fillResponse) {
        binnedRegistry.add({EfficiencyMult.data(), effLabels.c_str(), {HistType::kTHnSparseF, effAxes}});
        binnedRegistry.add({SplitMult.data(), " ; N_{gen} ; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});

        std::string reLabels{" ; N_{rec}; N_{PV cont}; N_{gen}"};
        std::vector<AxisSpec> reAxes{MultAxis, MultAxis, MultAxis};
        if (addFT0) {
          reLabels += " ; N_{FT0A}; N_{FT0C}";
          reAxes.push_back(FT0AAxis);
          reAxes.push_back(FT0CAxis);
        }
        if (addFDD) {
          reLabels += " ; N_{FDA}; N_{FDC}";
          reAxes.push_back(FDAAxis);
          reAxes.push_back(FDCAxis);
        }
        reLabels += " ; Z_{vtx} (cm); centrality; occupancy";
        reAxes.push_back(ZAxis);
        reAxes.push_back(CentAxis);
        reAxes.push_back(OccuAxis);
        binnedRegistry.add({Response.data(), reLabels.c_str(), {HistType::kTHnSparseF, reAxes}});
      }

      auto heff = binnedRegistry.get<THnSparse>(HIST(Efficiency));
      auto* x = heff->GetAxis(0);
      x->SetBinLabel(static_cast<int>(EvEffBins::kGen), EvEffBinLabels[static_cast<int>(EvEffBins::kGen)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kGengt0), EvEffBinLabels[static_cast<int>(EvEffBins::kGengt0)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kRec), EvEffBinLabels[static_cast<int>(EvEffBins::kRec)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelected), EvEffBinLabels[static_cast<int>(EvEffBins::kSelected)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelectedgt0), EvEffBinLabels[static_cast<int>(EvEffBins::kSelectedgt0)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelectedPVgt0), EvEffBinLabels[static_cast<int>(EvEffBins::kSelectedPVgt0)].data());
    }

    if (doprocessTrackEfficiencyAmbiguous || doprocessTrackEfficiency) {
      inclusiveRegistry.add({PtGen.data(), " ; p_{T} (GeV/c); occupancy", {HistType::kTH2F, {PtAxisEff, OccuAxis}}});
      inclusiveRegistry.add({PtGenNoEtaCut.data(), " ; p_{T} (GeV/c); occupancy", {HistType::kTH2F, {PtAxisEff, OccuAxis}}});
      inclusiveRegistry.add({PtEfficiency.data(), " ; p_{T} (GeV/c); occupancy", {HistType::kTH2F, {PtAxisEff, OccuAxis}}});
      inclusiveRegistry.add({PtEfficiencyNoEtaCut.data(), " ; p_{T} (GeV/c); occupancy", {HistType::kTH2F, {PtAxisEff, OccuAxis}}});
      inclusiveRegistry.add({PtEfficiencyFakes.data(), " ; p_{T} (GeV/c); occupancy", {HistType::kTH2F, {PtAxisEff, OccuAxis}}});
      for (auto i = 0u; i < speciesIds.size(); ++i) {
        inclusiveRegistry.add({fmt::format(PtGenF.data(), species[i]).c_str(), " ; p_{T} (GeV/c); occupancy", {HistType::kTH2F, {PtAxisEff, OccuAxis}}});
        inclusiveRegistry.add({fmt::format(PtEfficiencyF.data(), species[i]).c_str(), " ; p_{T} (GeV/c); occupancy", {HistType::kTH2F, {PtAxisEff, OccuAxis}}});
      }
    }
    if (doprocessTrackEfficiencyAmbiguousCentralityFT0M || doprocessTrackEfficiencyCentralityFT0M || doprocessTrackEfficiencyAmbiguousCentralityFT0C || doprocessTrackEfficiencyCentralityFT0C) {
      binnedRegistry.add({PtGen.data(), " ; p_{T} (GeV/c); centrality; occupancy", {HistType::kTHnSparseF, {PtAxisEff, CentAxis, OccuAxis}}});
      binnedRegistry.add({PtGenNoEtaCut.data(), " ; p_{T} (GeV/c); centrality; occupancy", {HistType::kTHnSparseF, {PtAxisEff, CentAxis, OccuAxis}}});
      binnedRegistry.add({PtEfficiency.data(), " ; p_{T} (GeV/c); centrality; occupancy", {HistType::kTHnSparseF, {PtAxisEff, CentAxis, OccuAxis}}});
      binnedRegistry.add({PtEfficiencyNoEtaCut.data(), " ; p_{T} (GeV/c); centrality; occupancy", {HistType::kTHnSparseF, {PtAxisEff, CentAxis, OccuAxis}}});
      binnedRegistry.add({PtEfficiencyFakes.data(), " ; p_{T} (GeV/c); centrality; occupancy", {HistType::kTHnSparseF, {PtAxisEff, CentAxis, OccuAxis}}});
      for (auto i = 0u; i < speciesIds.size(); ++i) {
        binnedRegistry.add({fmt::format(PtGenF.data(), species[i]).c_str(), " ; p_{T} (GeV/c); centrality; occupancy", {HistType::kTHnSparseF, {PtAxisEff, CentAxis, OccuAxis}}});
        binnedRegistry.add({fmt::format(PtEfficiencyF.data(), species[i]).c_str(), " ; p_{T} (GeV/c); centrality; occupancy", {HistType::kTHnSparseF, {PtAxisEff, CentAxis, OccuAxis}}});
      }
    }
    if (doprocessTrackEfficiencyIndexed) {
      inclusiveRegistry.add({PhiEtaGenDuplicates.data(), "; #varphi; #eta; occupancy; tracks", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, OccuAxis}}});
      inclusiveRegistry.add({PhiEtaDuplicates.data(), "; #varphi; #eta; occupancy; tracks", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, OccuAxis}}});
      inclusiveRegistry.add({PtGenIdx.data(), " ; p_{T} (GeV/c); occupancy", {HistType::kTH2F, {PtAxisEff, OccuAxis}}});
      inclusiveRegistry.add({PtGenIdxNoEtaCut.data(), " ; p_{T} (GeV/c); occupancy", {HistType::kTH2F, {PtAxisEff, OccuAxis}}});
      inclusiveRegistry.add({PtEfficiencyIdx.data(), " ; p_{T} (GeV/c); occupancy", {HistType::kTH2F, {PtAxisEff, OccuAxis}}});
      inclusiveRegistry.add({PtEfficiencyIdxNoEtaCut.data(), " ; p_{T} (GeV/c); occupancy", {HistType::kTH2F, {PtAxisEff, OccuAxis}}});
      inclusiveRegistry.add({PtEfficiencySecondariesIdx.data(), " ; p_{T} (GeV/c); occupancy", {HistType::kTH2F, {PtAxisEff, OccuAxis}}});
      inclusiveRegistry.add({PtEfficiencySecondariesIdxNoEtaCut.data(), " ; p_{T} (GeV/c); occupancy", {HistType::kTH2F, {PtAxisEff, OccuAxis}}});
      inclusiveRegistry.add({Mask.data(), " ; bit; occupancy", {HistType::kTH2F, {{17, -0.5, 16.5}, OccuAxis}}});
      inclusiveRegistry.add({ITSlayers.data(), " ; layer; occupancy", {HistType::kTH2F, {{8, 0.5, 8.5}, OccuAxis}}});
      for (auto i = 0u; i < speciesIds.size(); ++i) {
        inclusiveRegistry.add({fmt::format(PtGenIdxF.data(), species[i]).c_str(), " ; p_{T} (GeV/c); occupancy", {HistType::kTH2F, {PtAxisEff, OccuAxis}}});
        inclusiveRegistry.add({fmt::format(PtEfficiencyIdxF.data(), species[i]).c_str(), " ; p_{T} (GeV/c); occupancy", {HistType::kTH2F, {PtAxisEff, OccuAxis}}});
      }
    }
    if (doprocessTrackEfficiencyIndexedCentralityFT0M || doprocessTrackEfficiencyIndexedCentralityFT0C) {
      binnedRegistry.add({PhiEtaGenDuplicates.data(), "; #varphi; #eta; centrality; occupancy; tracks", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis, OccuAxis}}});
      binnedRegistry.add({PhiEtaDuplicates.data(), "; #varphi; #eta; centrality; occupancy; tracks", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis, OccuAxis}}});
      binnedRegistry.add({PtGenIdx.data(), " ; p_{T} (GeV/c); centrality; occupancy", {HistType::kTHnSparseF, {PtAxisEff, CentAxis, OccuAxis}}});
      binnedRegistry.add({PtGenIdxNoEtaCut.data(), " ; p_{T} (GeV/c); centrality; occupancy", {HistType::kTHnSparseF, {PtAxisEff, CentAxis, OccuAxis}}});
      binnedRegistry.add({PtEfficiencyIdx.data(), " ; p_{T} (GeV/c); centrality; occupancy", {HistType::kTHnSparseF, {PtAxisEff, CentAxis, OccuAxis}}});
      binnedRegistry.add({PtEfficiencyIdxNoEtaCut.data(), " ; p_{T} (GeV/c); centrality; occupancy", {HistType::kTHnSparseF, {PtAxisEff, CentAxis, OccuAxis}}});
      binnedRegistry.add({PtEfficiencySecondariesIdx.data(), " ; p_{T} (GeV/c); centrality; occupancy", {HistType::kTHnSparseF, {PtAxisEff, CentAxis, OccuAxis}}});
      binnedRegistry.add({PtEfficiencySecondariesIdxNoEtaCut.data(), " ; p_{T} (GeV/c); centrality; occupancy", {HistType::kTHnSparseF, {PtAxisEff, CentAxis, OccuAxis}}});
      binnedRegistry.add({Mask.data(), " ; bit; centrality; occupancy", {HistType::kTHnSparseF, {{17, -0.5, 16.5}, CentAxis, OccuAxis}}});
      binnedRegistry.add({ITSlayers.data(), " ; layer; centrality; occupancy", {HistType::kTHnSparseF, {{8, 0.5, 8.5}, CentAxis, OccuAxis}}});
      for (auto i = 0u; i < speciesIds.size(); ++i) {
        binnedRegistry.add({fmt::format(PtGenIdxF.data(), species[i]).c_str(), " ; p_{T} (GeV/c); centrality; occupancy", {HistType::kTHnSparseF, {PtAxisEff, CentAxis, OccuAxis}}});
        binnedRegistry.add({fmt::format(PtEfficiencyIdxF.data(), species[i]).c_str(), " ; p_{T} (GeV/c); centrality; occupancy", {HistType::kTHnSparseF, {PtAxisEff, CentAxis, OccuAxis}}});
      }
    }
  }

  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  template <typename C>
  void processEventStatGeneral(FullBCs const& bcs, C const& collisions)
  {
    std::vector<int64_t> colids;
    for (auto& bc : bcs) {
      if (!useEvSel || isBCSelected(bc)) {
        commonRegistry.fill(HIST(BCSelection), 1.);
        colids.clear();
        for (auto& collision : collisions) {
          if (collision.has_foundBC()) {
            if (collision.foundBCId() == bc.globalIndex()) {
              colids.push_back(collision.globalIndex());
            }
          } else if (collision.bcId() == bc.globalIndex()) {
            colids.push_back(collision.globalIndex());
          }
        }
        LOGP(debug, "BC {} has {} collisions", bc.globalBC(), colids.size());
        if (!colids.empty()) {
          commonRegistry.fill(HIST(BCSelection), 2.);
          if (colids.size() > 1) {
            commonRegistry.fill(HIST(BCSelection), 3.);
          }
        }
        auto col = collisions.begin();
        for (auto& colid : colids) {
          col.moveByIndex(colid - col.globalIndex());
          float c = getRecoCent(c);
          if constexpr (has_reco_cent<C>) {
            binnedRegistry.fill(HIST(EventChi2), col.chi2(), c, col.trackOccupancyInTimeRange());
            binnedRegistry.fill(HIST(EventTimeRes), col.collisionTimeRes(), c, col.trackOccupancyInTimeRange());
          } else {
            inclusiveRegistry.fill(HIST(EventChi2), col.chi2(), col.trackOccupancyInTimeRange());
            inclusiveRegistry.fill(HIST(EventTimeRes), col.collisionTimeRes(), col.trackOccupancyInTimeRange());
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

  // require a mix of ITS+TPC and ITS-only tracks (filters on the same table are automatically combined with &&)
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

  template <typename C, bool fillHistos = true, typename T>
  int countTracks(T const& tracks, float z, float c, float o)
  {
    auto Ntrks = 0;
    for (auto& track : tracks) {
      if (std::abs(track.eta()) < estimatorEta) {
        ++Ntrks;
      }
      if constexpr (fillHistos) {
        if constexpr (has_reco_cent<C>) {
          binnedRegistry.fill(HIST(EtaZvtx), track.eta(), z, c, o);
          binnedRegistry.fill(HIST(PhiEta), track.phi(), track.eta(), c, o);
          binnedRegistry.fill(HIST(PtEta), track.pt(), track.eta(), c, o);
          binnedRegistry.fill(HIST(DCAXYPt), track.pt(), track.dcaXY(), c, o);
          binnedRegistry.fill(HIST(DCAZPt), track.pt(), track.dcaZ(), c, o);
        } else {
          inclusiveRegistry.fill(HIST(EtaZvtx), track.eta(), z, o);
          inclusiveRegistry.fill(HIST(PhiEta), track.phi(), track.eta(), o);
          inclusiveRegistry.fill(HIST(PtEta), track.pt(), track.eta(), o);
          inclusiveRegistry.fill(HIST(DCAXYPt), track.pt(), track.dcaXY(), o);
          inclusiveRegistry.fill(HIST(DCAZPt), track.pt(), track.dcaZ(), o);
        }
      }
    }
    return Ntrks;
  }

  template <typename C>
  void processCountingGeneral(
    typename C::iterator const& collision,
    FiTracks const& tracks)
  {
    float c = getRecoCent(collision);
    float o = collision.trackOccupancyInTimeRange();
    if constexpr (has_reco_cent<C>) {
      binnedRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kAll), c, o);
    } else {
      inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kAll), o);
    }

    if (!useEvSel || isCollisionSelected(collision)) {
      if constexpr (has_reco_cent<C>) {
        binnedRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelected), c, o);
      } else {
        inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelected), o);
      }
      auto z = collision.posZ();
      usedTracksIds.clear();

      auto groupPVContrib = pvContribTracksIUEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto INELgt0PV = groupPVContrib.size() > 0;

      auto Ntrks = countTracks<C>(tracks, z, c, o);
      if constexpr (has_reco_cent<C>) {
        if (Ntrks > 0 || INELgt0PV) {
          if (INELgt0PV) {
            binnedRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelectedPVgt0), c, o);
          }
          if (Ntrks > 0) {
            binnedRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelectedgt0), c, o);
          }
          for (auto& track : tracks) {
            if (Ntrks > 0) {
              binnedRegistry.fill(HIST(EtaZvtx_gt0), track.eta(), z, c, o);
            }
            if (INELgt0PV) {
              binnedRegistry.fill(HIST(EtaZvtx_PVgt0), track.eta(), z, c, o);
            }
          }
        }
        binnedRegistry.fill(HIST(NtrkZvtx), Ntrks, z, c, o);
        binnedRegistry.fill(HIST(NpvcZvtx), groupPVContrib.size(), z, c, o);
      } else {
        if (Ntrks > 0 || INELgt0PV) {
          if (INELgt0PV) {
            inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelectedPVgt0), o);
          }
          if (Ntrks > 0) {
            inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelectedgt0), o);
          }
          for (auto& track : tracks) {
            if (Ntrks > 0) {
              inclusiveRegistry.fill(HIST(EtaZvtx_gt0), track.eta(), z, o);
            }
            if (INELgt0PV) {
              inclusiveRegistry.fill(HIST(EtaZvtx_PVgt0), track.eta(), z, o);
            }
          }
        }
        inclusiveRegistry.fill(HIST(NtrkZvtx), Ntrks, z, o);
        inclusiveRegistry.fill(HIST(NpvcZvtx), groupPVContrib.size(), z, o);
      }
    } else {
      if constexpr (has_reco_cent<C>) {
        binnedRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kRejected), c, o);
      } else {
        inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kRejected), o);
      }
    }
  }

  template <typename C, bool fillHistos = true, typename T, typename AT>
  int countTracksAmbiguous(T const& tracks, AT const& atracks, float z, float c, float o)
  {
    auto Ntrks = 0;
    for (auto& track : atracks) {
      auto otrack = track.template track_as<T>();
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
      if (fillHistos) {
        if constexpr (has_reco_cent<C>) {
          binnedRegistry.fill(HIST(EtaZvtx), otrack.eta(), z, c, o);
          binnedRegistry.fill(HIST(PhiEta), otrack.phi(), otrack.eta(), c, o);
          binnedRegistry.fill(HIST(PtEta), otrack.pt(), otrack.eta(), c, o);
          binnedRegistry.fill(HIST(DCAXYPt), otrack.pt(), track.bestDCAXY(), c, o);
          binnedRegistry.fill(HIST(DCAZPt), otrack.pt(), track.bestDCAZ(), c, o);
        } else {
          inclusiveRegistry.fill(HIST(EtaZvtx), otrack.eta(), z, o);
          inclusiveRegistry.fill(HIST(PhiEta), otrack.phi(), otrack.eta(), o);
          inclusiveRegistry.fill(HIST(PtEta), otrack.pt(), otrack.eta(), o);
          inclusiveRegistry.fill(HIST(DCAXYPt), otrack.pt(), track.bestDCAXY(), o);
          inclusiveRegistry.fill(HIST(DCAZPt), otrack.pt(), track.bestDCAZ(), o);
        }
      }
      if (otrack.has_collision() && otrack.collisionId() != track.bestCollisionId()) {
        usedTracksIdsDF.emplace_back(track.trackId());
        if constexpr (fillHistos) {
          if constexpr (has_reco_cent<C>) {
            binnedRegistry.fill(HIST(ReassignedEtaZvtx), otrack.eta(), z, c, o);
            binnedRegistry.fill(HIST(ReassignedPhiEta), otrack.phi(), otrack.eta(), c, o);
            binnedRegistry.fill(HIST(ReassignedZvtxCorr), otrack.template collision_as<C>().posZ(), z, c, o);
            binnedRegistry.fill(HIST(ReassignedDCAXYPt), otrack.pt(), track.bestDCAXY(), c, o);
            binnedRegistry.fill(HIST(ReassignedDCAZPt), otrack.pt(), track.bestDCAZ(), c, o);
          } else {
            inclusiveRegistry.fill(HIST(ReassignedEtaZvtx), otrack.eta(), z, o);
            inclusiveRegistry.fill(HIST(ReassignedPhiEta), otrack.phi(), otrack.eta(), o);
            inclusiveRegistry.fill(HIST(ReassignedZvtxCorr), otrack.template collision_as<C>().posZ(), z, o);
            inclusiveRegistry.fill(HIST(ReassignedDCAXYPt), otrack.pt(), track.bestDCAXY(), o);
            inclusiveRegistry.fill(HIST(ReassignedDCAZPt), otrack.pt(), track.bestDCAZ(), o);
          }
        }
      } else if (!otrack.has_collision()) {
        if constexpr (fillHistos) {
          if constexpr (has_reco_cent<C>) {
            binnedRegistry.fill(HIST(ExtraEtaZvtx), otrack.eta(), z, c, o);
            binnedRegistry.fill(HIST(ExtraPhiEta), otrack.phi(), otrack.eta(), c, o);
            binnedRegistry.fill(HIST(ExtraDCAXYPt), otrack.pt(), track.bestDCAXY(), c, o);
            binnedRegistry.fill(HIST(ExtraDCAZPt), otrack.pt(), track.bestDCAZ(), c, o);
          } else {
            inclusiveRegistry.fill(HIST(ExtraEtaZvtx), otrack.eta(), z, o);
            inclusiveRegistry.fill(HIST(ExtraPhiEta), otrack.phi(), otrack.eta(), o);
            inclusiveRegistry.fill(HIST(ExtraDCAXYPt), otrack.pt(), track.bestDCAXY(), o);
            inclusiveRegistry.fill(HIST(ExtraDCAZPt), otrack.pt(), track.bestDCAZ(), o);
          }
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
      if constexpr (fillHistos) {
        if constexpr (has_reco_cent<C>) {
          binnedRegistry.fill(HIST(EtaZvtx), track.eta(), z, c, o);
          binnedRegistry.fill(HIST(PhiEta), track.phi(), track.eta(), c, o);
          binnedRegistry.fill(HIST(PtEta), track.pt(), track.eta(), c, o);
          binnedRegistry.fill(HIST(DCAXYPt), track.pt(), track.dcaXY(), c, o);
          binnedRegistry.fill(HIST(DCAZPt), track.pt(), track.dcaZ(), c, o);
        } else {
          inclusiveRegistry.fill(HIST(EtaZvtx), track.eta(), z, o);
          inclusiveRegistry.fill(HIST(PhiEta), track.phi(), track.eta(), o);
          inclusiveRegistry.fill(HIST(PtEta), track.pt(), track.eta(), o);
          inclusiveRegistry.fill(HIST(DCAXYPt), track.pt(), track.dcaXY(), o);
          inclusiveRegistry.fill(HIST(DCAZPt), track.pt(), track.dcaZ(), o);
        }
      }
    }
    return Ntrks;
  }

  template <typename C>
  void processCountingGeneralwAmbiguous(
    typename C::iterator const& collision,
    FiTracks const& tracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    float c = getRecoCent(collision);
    float o = collision.trackOccupancyInTimeRange();
    if constexpr (has_reco_cent<C>) {
      binnedRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kAll), c, o);
    } else {
      inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kAll), o);
    }

    if (!useEvSel || isCollisionSelected(collision)) {
      if constexpr (has_reco_cent<C>) {
        binnedRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelected), c, o);
      } else {
        inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelected), o);
      }
      auto z = collision.posZ();
      usedTracksIds.clear();

      auto groupPVContrib = pvContribTracksIUEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto INELgt0PV = groupPVContrib.size() > 0;

      auto Ntrks = countTracksAmbiguous<C>(tracks, atracks, z, c, o);
      if constexpr (has_reco_cent<C>) {
        if (Ntrks > 0 || INELgt0PV) {
          if (INELgt0PV) {
            binnedRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelectedPVgt0), c, o);
          }
          if (Ntrks > 0) {
            binnedRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelectedgt0), c, o);
          }
          for (auto& track : atracks) {
            if (Ntrks > 0) {
              binnedRegistry.fill(HIST(EtaZvtx_gt0), track.track_as<FiTracks>().eta(), z, c, o);
            }
            if (INELgt0PV) {
              binnedRegistry.fill(HIST(EtaZvtx_PVgt0), track.track_as<FiTracks>().eta(), z, c, o);
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
              binnedRegistry.fill(HIST(EtaZvtx_gt0), track.eta(), z, c, o);
            }
            if (INELgt0PV) {
              binnedRegistry.fill(HIST(EtaZvtx_PVgt0), track.eta(), z, c, o);
            }
          }
        }
        binnedRegistry.fill(HIST(NtrkZvtx), Ntrks, z, c, o);
        binnedRegistry.fill(HIST(NpvcZvtx), groupPVContrib.size(), z, c, o);
      } else {
        if (Ntrks > 0 || INELgt0PV) {
          if (INELgt0PV) {
            inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelectedPVgt0), o);
          }
          if (Ntrks > 0) {
            inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelectedgt0), o);
          }
          for (auto& track : atracks) {
            if (Ntrks > 0) {
              inclusiveRegistry.fill(HIST(EtaZvtx_gt0), track.track_as<FiTracks>().eta(), z, o);
            }
            if (INELgt0PV) {
              inclusiveRegistry.fill(HIST(EtaZvtx_PVgt0), track.track_as<FiTracks>().eta(), z, o);
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
              inclusiveRegistry.fill(HIST(EtaZvtx_gt0), track.eta(), z, o);
            }
            if (INELgt0PV) {
              inclusiveRegistry.fill(HIST(EtaZvtx_PVgt0), track.eta(), z, o);
            }
          }
        }
        inclusiveRegistry.fill(HIST(NtrkZvtx), Ntrks, z, o);
        inclusiveRegistry.fill(HIST(NpvcZvtx), groupPVContrib.size(), z, o);
      }
    } else {
      if constexpr (has_reco_cent<C>) {
        binnedRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kRejected), c, o);
      } else {
        inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kRejected), o);
      }
    }
  }

  void processCountingAmbiguous(
    ExCols::iterator const& collision,
    FiTracks const& tracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    processCountingGeneralwAmbiguous<ExCols>(collision, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCountingAmbiguous, "Count tracks", false);

  void processCounting(
    ExCols::iterator const& collision,
    FiTracks const& tracks)
  {
    processCountingGeneral<ExCols>(collision, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCounting, "Count tracks w/o ambiguous", false);

  using ExColsCentFT0C = soa::Join<aod::Collisions, aod::CentFT0Cs, aod::EvSels>;
  void processCountingAmbiguousCentralityFT0C(
    ExColsCentFT0C::iterator const& collision,
    FiTracks const& tracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    processCountingGeneralwAmbiguous<ExColsCentFT0C>(collision, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCountingAmbiguousCentralityFT0C, "Count tracks in FT0C centrality bins", false);

  void processCountingCentralityFT0C(
    ExColsCentFT0C::iterator const& collision,
    FiTracks const& tracks)
  {
    processCountingGeneral<ExColsCentFT0C>(collision, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCountingCentralityFT0C, "Count tracks in FT0C centrality bins w/o ambiguous", false);

  using ExColsCentFT0M = soa::Join<aod::Collisions, aod::CentFT0Ms, aod::EvSels>;
  void processCountingAmbiguousCentralityFT0M(
    ExColsCentFT0M::iterator const& collision,
    FiTracks const& tracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    processCountingGeneralwAmbiguous<ExColsCentFT0M>(collision, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCountingAmbiguousCentralityFT0M, "Count tracks in FT0M centrality bins", false);

  void processCountingCentralityFT0M(
    ExColsCentFT0M::iterator const& collision,
    FiTracks const& tracks)
  {
    processCountingGeneral<ExColsCentFT0M>(collision, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processCountingCentralityFT0M, "Count tracks in FT0M centrality bins w/o ambiguous", false);

  using Particles = soa::Filtered<aod::McParticles>;
  using LabeledTracksEx = soa::Join<LabeledTracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
  using FiLTracks = soa::Filtered<LabeledTracksEx>;
  using ParticlesI = soa::Filtered<soa::Join<aod::McParticles, aod::ParticlesToTracks>>;
  expressions::Filter primaries = ncheckbit(aod::mcparticle::flags, (uint8_t)o2::aod::mcparticle::enums::PhysicalPrimary);

  bool isChargedParticle(int code)
  {
    auto p = pdg->GetParticle(code);
    auto charge = 0.;
    if (p != nullptr) {
      charge = p->Charge();
    }
    return std::abs(charge) >= 3.;
  }

  template <typename C, typename MC>
  void processTrackEfficiencyIndexedGeneral(
    typename soa::Join<C, aod::McCollisionLabels>::iterator const& collision,
    MC const&, ParticlesI const& particles,
    FiLTracks const& /*tracks*/)
  {
    if (useEvSel && !isCollisionSelected(collision)) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }
    float c_rec = getRecoCent(collision);
    float o = collision.trackOccupancyInTimeRange();
    auto mcCollision = collision.mcCollision();
    float c_gen = getSimCent(mcCollision);
    if (c_gen < 0 && c_rec >= 0) {
      c_gen = c_rec;
    }

    auto sample = particles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);

    for (auto& particle : sample) {
      if (!isChargedParticle(particle.pdgCode())) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        binnedRegistry.fill(HIST(PtGenIdxNoEtaCut), particle.pt(), c_gen, o);
        if (std::abs(particle.eta()) < estimatorEta) {
          binnedRegistry.fill(HIST(PtGenIdx), particle.pt(), c_gen, o);
          if (particle.pdgCode() == speciesIds[0]) {
            binnedRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtGenIdxSuff), particle.pt(), c_gen, o);
          } else if (particle.pdgCode() == speciesIds[1]) {
            binnedRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtGenIdxSuff), particle.pt(), c_gen, o);
          } else if (particle.pdgCode() == speciesIds[2]) {
            binnedRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtGenIdxSuff), particle.pt(), c_gen, o);
          } else if (particle.pdgCode() == speciesIds[3]) {
            binnedRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtGenIdxSuff), particle.pt(), c_gen, o);
          }
        }
      } else {
        inclusiveRegistry.fill(HIST(PtGenIdxNoEtaCut), particle.pt(), o);
        if (std::abs(particle.eta()) < estimatorEta) {
          inclusiveRegistry.fill(HIST(PtGenIdx), particle.pt(), o);
          if (particle.pdgCode() == speciesIds[0]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtGenIdxSuff), particle.pt(), o);
          } else if (particle.pdgCode() == speciesIds[1]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtGenIdxSuff), particle.pt(), o);
          } else if (particle.pdgCode() == speciesIds[2]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtGenIdxSuff), particle.pt(), o);
          } else if (particle.pdgCode() == speciesIds[3]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtGenIdxSuff), particle.pt(), o);
          }
        }
      }

      if (particle.has_tracks()) {
        auto counted = false;
        auto countedNoEtaCut = false;
        auto counter = 0;
        auto relatedTracks = particle.template filtered_tracks_as<FiLTracks>();
        for (auto const& track : relatedTracks) {
          ++counter;
          if constexpr (has_reco_cent<C>) {
            if (!countedNoEtaCut) {
              binnedRegistry.fill(HIST(PtEfficiencyIdxNoEtaCut), particle.pt(), c_gen, o);
              countedNoEtaCut = true;
            }
            if (std::abs(track.eta()) < estimatorEta) {
              if (!counted) {
                binnedRegistry.fill(HIST(PtEfficiencyIdx), particle.pt(), c_gen, o);
                if (particle.pdgCode() == speciesIds[0]) {
                  binnedRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtEffIdxSuff), particle.pt(), c_gen, o);
                } else if (particle.pdgCode() == speciesIds[1]) {
                  binnedRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtEffIdxSuff), particle.pt(), c_gen, o);
                } else if (particle.pdgCode() == speciesIds[2]) {
                  binnedRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtEffIdxSuff), particle.pt(), c_gen, o);
                } else if (particle.pdgCode() == speciesIds[3]) {
                  binnedRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtEffIdxSuff), particle.pt(), c_gen, o);
                }
                counted = true;
              }
            }
            if (counter > 1) {
              binnedRegistry.fill(HIST(PtEfficiencySecondariesIdxNoEtaCut), particle.pt(), c_gen, o);
              if (std::abs(track.eta()) < estimatorEta) {
                binnedRegistry.fill(HIST(PtEfficiencySecondariesIdx), particle.pt(), c_gen, o);
              }
            }
          } else {
            if (!countedNoEtaCut) {
              inclusiveRegistry.fill(HIST(PtEfficiencyIdxNoEtaCut), particle.pt(), o);
              countedNoEtaCut = true;
            }
            if (std::abs(track.eta()) < estimatorEta) {
              if (!counted) {
                inclusiveRegistry.fill(HIST(PtEfficiencyIdx), particle.pt(), o);
                if (particle.pdgCode() == speciesIds[0]) {
                  inclusiveRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtEffIdxSuff), particle.pt(), o);
                } else if (particle.pdgCode() == speciesIds[1]) {
                  inclusiveRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtEffIdxSuff), particle.pt(), o);
                } else if (particle.pdgCode() == speciesIds[2]) {
                  inclusiveRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtEffIdxSuff), particle.pt(), o);
                } else if (particle.pdgCode() == speciesIds[3]) {
                  inclusiveRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtEffIdxSuff), particle.pt(), o);
                }
                counted = true;
              }
            }
            if (counter > 1) {
              inclusiveRegistry.fill(HIST(PtEfficiencySecondariesIdxNoEtaCut), particle.pt(), o);
              if (std::abs(track.eta()) < estimatorEta) {
                inclusiveRegistry.fill(HIST(PtEfficiencySecondariesIdx), particle.pt(), o);
              }
            }
          }
        }
        if constexpr (has_reco_cent<C>) {
          for (auto const& track : relatedTracks) {
            for (auto layer = 0; layer < 7; ++layer) {
              if (track.itsClusterMap() & (uint8_t(1) << layer)) {
                binnedRegistry.fill(HIST(ITSlayers), layer + 1, c_gen, o);
              }
            }
            auto hasbit = false;
            for (auto bit = 0; bit < 16; ++bit) {
              if (track.mcMask() & (uint8_t(1) << bit)) {
                binnedRegistry.fill(HIST(Mask), bit, c_gen, o);
                hasbit = true;
              }
            }
            if (!hasbit) {
              binnedRegistry.fill(HIST(Mask), 16, c_gen, o);
            }
          }
          if (relatedTracks.size() > 1) {
            binnedRegistry.fill(HIST(PhiEtaGenDuplicates), particle.phi(), particle.eta(), c_gen, o);
            for (auto const& track : relatedTracks) {
              binnedRegistry.fill(HIST(PhiEtaDuplicates), track.phi(), track.eta(), c_gen, o);
            }
          }
        } else {
          for (auto const& track : relatedTracks) {
            for (auto layer = 0; layer < 7; ++layer) {
              if (track.itsClusterMap() & (uint8_t(1) << layer)) {
                inclusiveRegistry.fill(HIST(ITSlayers), layer + 1, o);
              }
            }
            auto hasbit = false;
            for (auto bit = 0; bit < 16; ++bit) {
              if (track.mcMask() & (uint8_t(1) << bit)) {
                inclusiveRegistry.fill(HIST(Mask), bit, o);
                hasbit = true;
              }
            }
            if (!hasbit) {
              inclusiveRegistry.fill(HIST(Mask), 16, o);
            }
          }
          if (relatedTracks.size() > 1) {
            inclusiveRegistry.fill(HIST(PhiEtaGenDuplicates), particle.phi(), particle.eta(), o);
            for (auto const& track : relatedTracks) {
              inclusiveRegistry.fill(HIST(PhiEtaDuplicates), track.phi(), track.eta(), o);
            }
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

  void processTrackEfficiencyIndexedCentralityFT0M(
    soa::Join<ExColsCentFT0M, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, ParticlesI const& particles,
    FiLTracks const& tracks)
  {
    processTrackEfficiencyIndexedGeneral<ExColsCentFT0M, aod::McCollisions>(collision, mccollisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processTrackEfficiencyIndexedCentralityFT0M, "Calculate tracking efficiency vs pt (indexed, FT0M binned)", false);

  void processTrackEfficiencyIndexedCentralityFT0C(
    soa::Join<ExColsCentFT0C, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, ParticlesI const& particles,
    FiLTracks const& tracks)
  {
    processTrackEfficiencyIndexedGeneral<ExColsCentFT0C, aod::McCollisions>(collision, mccollisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processTrackEfficiencyIndexedCentralityFT0C, "Calculate tracking efficiency vs pt (indexed, FT0C binned)", false);

  template <typename C, typename MC>
  void processTrackEfficiencyGeneralAmbiguous(
    typename soa::Join<C, aod::McCollisionLabels>::iterator const& collision,
    MC const&, Particles const& particles,
    FiLTracks const& tracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    if (useEvSel && !isCollisionSelected(collision)) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }
    float c_rec = getRecoCent(collision);
    float o = collision.trackOccupancyInTimeRange();
    auto mcCollision = collision.mcCollision();
    float c_gen = getSimCent(mcCollision);
    if (c_gen < 0 && c_rec >= 0) {
      c_gen = c_rec;
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
        if constexpr (has_reco_cent<C>) {
          binnedRegistry.fill(HIST(PtEfficiencyNoEtaCut), particle.pt(), c_gen, o);
        } else {
          inclusiveRegistry.fill(HIST(PtEfficiencyNoEtaCut), particle.pt(), o);
        }
        if (std::abs(otrack.eta()) < estimatorEta) {
          if constexpr (has_reco_cent<C>) {
            binnedRegistry.fill(HIST(PtEfficiency), particle.pt(), c_gen, o);
            if (particle.pdgCode() == speciesIds[0]) {
              binnedRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtEffSuff), particle.pt(), c_gen, o);
            } else if (particle.pdgCode() == speciesIds[1]) {
              binnedRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtEffSuff), particle.pt(), c_gen, o);
            } else if (particle.pdgCode() == speciesIds[2]) {
              binnedRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtEffSuff), particle.pt(), c_gen, o);
            } else if (particle.pdgCode() == speciesIds[3]) {
              binnedRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtEffSuff), particle.pt(), c_gen, o);
            }
          } else {
            inclusiveRegistry.fill(HIST(PtEfficiency), particle.pt(), o);
            if (particle.pdgCode() == speciesIds[0]) {
              inclusiveRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtEffSuff), particle.pt(), o);
            } else if (particle.pdgCode() == speciesIds[1]) {
              inclusiveRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtEffSuff), particle.pt(), o);
            } else if (particle.pdgCode() == speciesIds[2]) {
              inclusiveRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtEffSuff), particle.pt(), o);
            } else if (particle.pdgCode() == speciesIds[3]) {
              inclusiveRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtEffSuff), particle.pt(), o);
            }
          }
        }
      } else {
        if constexpr (has_reco_cent<C>) {
          binnedRegistry.fill(HIST(PtEfficiencyFakes), otrack.pt(), c_gen, o);
        } else {
          inclusiveRegistry.fill(HIST(PtEfficiencyFakes), otrack.pt(), o);
        }
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
        if constexpr (has_reco_cent<C>) {
          binnedRegistry.fill(HIST(PtEfficiencyNoEtaCut), particle.pt(), c_gen, o);
        } else {
          inclusiveRegistry.fill(HIST(PtEfficiencyNoEtaCut), particle.pt(), o);
        }
        if (std::abs(track.eta()) < estimatorEta) {
          if constexpr (has_reco_cent<C>) {
            binnedRegistry.fill(HIST(PtEfficiency), particle.pt(), c_gen, o);
            if (particle.pdgCode() == speciesIds[0]) {
              binnedRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtEffSuff), particle.pt(), c_gen, o);
            } else if (particle.pdgCode() == speciesIds[1]) {
              binnedRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtEffSuff), particle.pt(), c_gen, o);
            } else if (particle.pdgCode() == speciesIds[2]) {
              binnedRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtEffSuff), particle.pt(), c_gen, o);
            } else if (particle.pdgCode() == speciesIds[3]) {
              binnedRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtEffSuff), particle.pt(), c_gen, o);
            }
          } else {
            inclusiveRegistry.fill(HIST(PtEfficiency), particle.pt(), o);
            if (particle.pdgCode() == speciesIds[0]) {
              inclusiveRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtEffSuff), particle.pt(), o);
            } else if (particle.pdgCode() == speciesIds[1]) {
              inclusiveRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtEffSuff), particle.pt(), o);
            } else if (particle.pdgCode() == speciesIds[2]) {
              inclusiveRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtEffSuff), particle.pt(), o);
            } else if (particle.pdgCode() == speciesIds[3]) {
              inclusiveRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtEffSuff), particle.pt(), o);
            }
          }
        }
      } else {
        if constexpr (has_reco_cent<C>) {
          binnedRegistry.fill(HIST(PtEfficiencyFakes), track.pt(), c_gen, o);
        } else {
          inclusiveRegistry.fill(HIST(PtEfficiencyFakes), track.pt(), o);
        }
      }
    }

    for (auto& particle : particlesPerCol) {
      if (!isChargedParticle(particle.pdgCode())) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        binnedRegistry.fill(HIST(PtGenNoEtaCut), particle.pt(), c_gen, o);
        if (std::abs(particle.eta()) < estimatorEta) {
          binnedRegistry.fill(HIST(PtGen), particle.pt(), c_gen, o);
          if (particle.pdgCode() == speciesIds[0]) {
            binnedRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtGenSuff), particle.pt(), c_gen, o);
          } else if (particle.pdgCode() == speciesIds[1]) {
            binnedRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtGenSuff), particle.pt(), c_gen, o);
          } else if (particle.pdgCode() == speciesIds[2]) {
            binnedRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtGenSuff), particle.pt(), c_gen, o);
          } else if (particle.pdgCode() == speciesIds[3]) {
            binnedRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtGenSuff), particle.pt(), c_gen, o);
          }
        }
      } else {
        inclusiveRegistry.fill(HIST(PtGenNoEtaCut), particle.pt());
        if (std::abs(particle.eta()) < estimatorEta) {
          inclusiveRegistry.fill(HIST(PtGen), particle.pt(), o);
          if (particle.pdgCode() == speciesIds[0]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtGenSuff), particle.pt(), o);
          } else if (particle.pdgCode() == speciesIds[1]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtGenSuff), particle.pt(), o);
          } else if (particle.pdgCode() == speciesIds[2]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtGenSuff), particle.pt(), o);
          } else if (particle.pdgCode() == speciesIds[3]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtGenSuff), particle.pt(), o);
          }
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
    if (useEvSel && !isCollisionSelected(collision)) {
      return;
    }
    if (!collision.has_mcCollision()) {
      return;
    }
    float c_rec = getRecoCent(collision);
    float o = collision.trackOccupancyInTimeRange();
    auto mcCollision = collision.mcCollision();
    float c_gen = getSimCent(mcCollision);
    if (c_gen < 0 && c_rec >= 0) {
      c_gen = c_rec;
    }

    auto particlesPerCol = particles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);

    usedTracksIds.clear();
    for (auto const& track : tracks) {
      if (track.has_mcParticle()) {
        auto particle = track.template mcParticle_as<Particles>();
        if constexpr (has_reco_cent<C>) {
          binnedRegistry.fill(HIST(PtEfficiencyNoEtaCut), particle.pt(), c_gen, o);
          if (std::abs(track.eta()) < estimatorEta) {
            binnedRegistry.fill(HIST(PtEfficiency), particle.pt(), c_gen, o);
            if (particle.pdgCode() == speciesIds[0]) {
              binnedRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtEffSuff), particle.pt(), c_gen, o);
            } else if (particle.pdgCode() == speciesIds[1]) {
              binnedRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtEffSuff), particle.pt(), c_gen, o);
            } else if (particle.pdgCode() == speciesIds[2]) {
              binnedRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtEffSuff), particle.pt(), c_gen, o);
            } else if (particle.pdgCode() == speciesIds[3]) {
              binnedRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtEffSuff), particle.pt(), c_gen, o);
            }
          }
        } else {
          inclusiveRegistry.fill(HIST(PtEfficiencyNoEtaCut), particle.pt(), o);
          if (std::abs(track.eta()) < estimatorEta) {
            inclusiveRegistry.fill(HIST(PtEfficiency), particle.pt(), o);
            if (particle.pdgCode() == speciesIds[0]) {
              inclusiveRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtEffSuff), particle.pt(), o);
            } else if (particle.pdgCode() == speciesIds[1]) {
              inclusiveRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtEffSuff), particle.pt(), o);
            } else if (particle.pdgCode() == speciesIds[2]) {
              inclusiveRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtEffSuff), particle.pt(), o);
            } else if (particle.pdgCode() == speciesIds[3]) {
              inclusiveRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtEffSuff), particle.pt(), o);
            }
          }
        }
      } else {
        if constexpr (has_reco_cent<C>) {
          binnedRegistry.fill(HIST(PtEfficiencyFakes), track.pt(), c_gen, o);
        } else {
          inclusiveRegistry.fill(HIST(PtEfficiencyFakes), track.pt(), o);
        }
      }
    }

    for (auto& particle : particlesPerCol) {
      if (!isChargedParticle(particle.pdgCode())) {
        continue;
      }
      if constexpr (has_reco_cent<C>) {
        binnedRegistry.fill(HIST(PtGenNoEtaCut), particle.pt(), c_gen, o);
        if (std::abs(particle.eta()) < estimatorEta) {
          binnedRegistry.fill(HIST(PtGen), particle.pt(), c_gen, o);
          if (particle.pdgCode() == speciesIds[0]) {
            binnedRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtGenSuff), particle.pt(), c_gen, o);
          } else if (particle.pdgCode() == speciesIds[1]) {
            binnedRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtGenSuff), particle.pt(), c_gen, o);
          } else if (particle.pdgCode() == speciesIds[2]) {
            binnedRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtGenSuff), particle.pt(), c_gen, o);
          } else if (particle.pdgCode() == speciesIds[3]) {
            binnedRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtGenSuff), particle.pt(), c_gen, o);
          }
        }
      } else {
        inclusiveRegistry.fill(HIST(PtGenNoEtaCut), particle.pt(), o);
        if (std::abs(particle.eta()) < estimatorEta) {
          inclusiveRegistry.fill(HIST(PtGen), particle.pt(), o);
          if (particle.pdgCode() == speciesIds[0]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtGenSuff), particle.pt(), o);
          } else if (particle.pdgCode() == speciesIds[1]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtGenSuff), particle.pt(), o);
          } else if (particle.pdgCode() == speciesIds[2]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtGenSuff), particle.pt(), o);
          } else if (particle.pdgCode() == speciesIds[3]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtGenSuff), particle.pt(), o);
          }
        }
      }
    }
  }

  void processTrackEfficiencyAmbiguous(
    soa::Join<ExCols, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, Particles const& mcParticles,
    FiLTracks const& filtracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    processTrackEfficiencyGeneralAmbiguous<ExCols, aod::McCollisions>(collision, mccollisions, mcParticles, filtracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processTrackEfficiencyAmbiguous, "Calculate tracking efficiency vs pt", false);

  void processTrackEfficiencyAmbiguousCentralityFT0M(
    soa::Join<ExColsCentFT0M, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, Particles const& mcParticles,
    FiLTracks const& filtracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    processTrackEfficiencyGeneralAmbiguous<ExColsCentFT0M, aod::McCollisions>(collision, mccollisions, mcParticles, filtracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processTrackEfficiencyAmbiguousCentralityFT0M, "Calculate tracking efficiency vs pt (FT0M binned)", false);

  void processTrackEfficiencyAmbiguousCentralityFT0C(
    soa::Join<ExColsCentFT0C, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, Particles const& mcParticles,
    FiLTracks const& filtracks,
    soa::SmallGroups<ReTracks> const& atracks)
  {
    processTrackEfficiencyGeneralAmbiguous<ExColsCentFT0C, aod::McCollisions>(collision, mccollisions, mcParticles, filtracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processTrackEfficiencyAmbiguousCentralityFT0C, "Calculate tracking efficiency vs pt (FT0C binned)", false);

  void processTrackEfficiency(
    soa::Join<ExCols, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, Particles const& mcParticles,
    FiLTracks const& filtracks)
  {
    processTrackEfficiencyGeneral<ExCols, aod::McCollisions>(collision, mccollisions, mcParticles, filtracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processTrackEfficiency, "Calculate tracking efficiency vs pt w/o ambiguous", false);

  void processTrackEfficiencyCentralityFT0M(
    soa::Join<ExColsCentFT0M, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, Particles const& mcParticles,
    FiLTracks const& filtracks)
  {
    processTrackEfficiencyGeneral<ExColsCentFT0M, aod::McCollisions>(collision, mccollisions, mcParticles, filtracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processTrackEfficiencyCentralityFT0M, "Calculate tracking efficiency vs pt w/o ambiguous (FT0M binned)", false);

  void processTrackEfficiencyCentralityFT0C(
    soa::Join<ExColsCentFT0C, aod::McCollisionLabels>::iterator const& collision,
    aod::McCollisions const& mccollisions, Particles const& mcParticles,
    FiLTracks const& filtracks)
  {
    processTrackEfficiencyGeneral<ExColsCentFT0C, aod::McCollisions>(collision, mccollisions, mcParticles, filtracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processTrackEfficiencyCentralityFT0C, "Calculate tracking efficiency vs pt w/o ambiguous (FT0C binned)", false);

  template <typename CIT>
  void fillFIT(CIT const& collision, std::vector<float>& ft0as, std::vector<float>& ft0cs, std::vector<float>& fddas, std::vector<float>& fddcs)
  {
    if (addFT0 && collision.has_foundFT0()) {
      auto ft0 = collision.foundFT0();
      float tA = 0;
      float tC = 0;
      for (auto i = 0u; i < ft0.amplitudeA().size(); ++i) {
        tA += ft0.amplitudeA()[i];
      }
      for (auto i = 0u; i < ft0.amplitudeC().size(); ++i) {
        tC += ft0.amplitudeC()[i];
      }
      ft0as.emplace_back(tA);
      ft0cs.emplace_back(tC);
    } else {
      ft0as.emplace_back(-1);
      ft0cs.emplace_back(-1);
    }
    if (addFDD && collision.has_foundFDD()) {
      auto fdd = collision.foundFDD();
      float tA = 0;
      float tC = 0;
      for (auto i = 0u; i < 8; ++i) {
        tA += fdd.chargeA()[i];
      }
      for (auto i = 0u; i < 8; ++i) {
        tC += fdd.chargeC()[i];
      }
      fddas.emplace_back(tA);
      fddcs.emplace_back(tC);
    } else {
      fddas.emplace_back(-1);
      fddcs.emplace_back(-1);
    }
  }

  template <bool hasCent, typename Ps>
  void fillParticleHistos(Ps const& particles, float z,
                          int const nCharged, float c,
                          bool const atLeastOne, bool const atLeastOne_gt0, bool const atLeastOne_PVgt0)
  {
    for (auto& particle : particles) {
      if (!isChargedParticle(particle.pdgCode())) {
        continue;
      }
      if constexpr (hasCent) {
        binnedRegistry.fill(HIST(EtaZvtxGen_t), particle.eta(), z, c);
        binnedRegistry.fill(HIST(PtEtaGen), particle.pt(), particle.eta(), c);
      } else {
        inclusiveRegistry.fill(HIST(EtaZvtxGen_t), particle.eta(), z);
        inclusiveRegistry.fill(HIST(PtEtaGen), particle.pt(), particle.eta());
      }
      if (nCharged > 0) {
        if constexpr (hasCent) {
          binnedRegistry.fill(HIST(EtaZvtxGen_gt0t), particle.eta(), z, c);
        } else {
          inclusiveRegistry.fill(HIST(EtaZvtxGen_gt0t), particle.eta(), z);
        }
      }
      if (atLeastOne) {
        if constexpr (hasCent) {
          binnedRegistry.fill(HIST(EtaZvtxGen), particle.eta(), z, c);
          if (atLeastOne_gt0) {
            binnedRegistry.fill(HIST(EtaZvtxGen_gt0), particle.eta(), z, c);
          }
          if (atLeastOne_PVgt0) {
            binnedRegistry.fill(HIST(EtaZvtxGen_PVgt0), particle.eta(), z, c);
          }
          binnedRegistry.fill(HIST(PhiEtaGen), particle.phi(), particle.eta(), c);
        } else {
          inclusiveRegistry.fill(HIST(EtaZvtxGen), particle.eta(), z);
          if (atLeastOne_gt0) {
            inclusiveRegistry.fill(HIST(EtaZvtxGen_gt0), particle.eta(), z);
          }
          if (atLeastOne_PVgt0) {
            inclusiveRegistry.fill(HIST(EtaZvtxGen_PVgt0), particle.eta(), z);
          }
          inclusiveRegistry.fill(HIST(PhiEtaGen), particle.phi(), particle.eta());
        }
      }
    }
  }

  template <typename Ps>
  int countParticles(Ps const& particles)
  {
    auto nCharged = 0;
    for (auto& particle : particles) {
      if (!isChargedParticle(particle.pdgCode())) {
        continue;
      }
      if (std::abs(particle.eta()) >= estimatorEta) {
        continue;
      }
      nCharged++;
    }
    return nCharged;
  }

  std::vector<int> NrecPerCol;
  std::vector<float> c_recPerCol;
  std::vector<float> OccuPerCol;
  std::vector<int> NPVPerCol;
  std::vector<float> NFT0APerCol;
  std::vector<float> NFT0CPerCol;
  std::vector<float> NFDDAPerCol;
  std::vector<float> NFDDCPerCol;

  template <typename MC, typename C>
  void processGenGeneralAmbiguous(
    typename MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks)
  {
    float c_gen = -1;
    if constexpr (has_Centrality<MC>) {
      c_gen = getSimCent(mcCollision);
    } else {
      if constexpr (has_FT0C<C>) {
        c_gen = getGenCentFT0C(mcCollision);
      } else if (has_FT0M<C>) {
        c_gen = getGenCentFT0M(mcCollision);
      }
    }

    NrecPerCol.clear();
    c_recPerCol.clear();
    OccuPerCol.clear();
    NPVPerCol.clear();
    NFT0APerCol.clear();
    NFT0CPerCol.clear();
    NFDDAPerCol.clear();
    NFDDCPerCol.clear();

    bool atLeastOne = false;
    bool atLeastOne_gt0 = false;
    bool atLeastOne_PVgt0 = false;
    auto moreThanOne = 0;
    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());

    if constexpr (has_reco_cent<C>) {
      float min_c_rec = 2e2;
      for (auto& collision : collisions) {
        if (!useEvSel || isCollisionSelected(collision)) {
          float c = getRecoCent(collision);
          if (c < min_c_rec) {
            min_c_rec = c;
          }
        }
      }
      if constexpr (!has_Centrality<MC>) {
        if (c_gen < 0) {
          c_gen = min_c_rec; // if there is no generator centrality info, fall back to reco (from the largest reco collision)
        }
      }
    }
    float o_max = -1;
    for (auto& collision : collisions) {
      if (!useEvSel || isCollisionSelected(collision)) {
        auto o = collision.trackOccupancyInTimeRange();
        if (o > o_max) {
          o_max = o;
        }
      }
    }

    for (auto& collision : collisions) {
      usedTracksIds.clear();
      float c_rec = getRecoCent(collision);
      float o = collision.trackOccupancyInTimeRange();
      if constexpr (has_reco_cent<C>) {
        binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kRec), c_gen, o);
      } else {
        inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kRec), o);
      }
      if (!useEvSel || isCollisionSelected(collision)) {
        c_recPerCol.emplace_back(c_rec);
        OccuPerCol.emplace_back(o);
        auto z = collision.posZ();
        ++moreThanOne;
        atLeastOne = true;

        auto groupPVcontrib = pvContribTracksIUEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
        if (groupPVcontrib.size() > 0) {
          if constexpr (has_reco_cent<C>) {
            binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedPVgt0), c_gen, o);
          } else {
            inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedPVgt0), o);
          }
          atLeastOne_PVgt0 = true;
        }

        auto perCollisionASample = atracks.sliceBy(perColU, collision.globalIndex());
        auto perCollisionSample = tracks.sliceBy(perCol, collision.globalIndex());
        auto Nrec = countTracksAmbiguous<C, false>(perCollisionSample, perCollisionASample, z, c_rec, o);
        NrecPerCol.emplace_back(Nrec);
        NPVPerCol.emplace_back(groupPVcontrib.size());
        fillFIT(collision, NFT0APerCol, NFT0CPerCol, NFDDAPerCol, NFDDCPerCol);

        if constexpr (has_reco_cent<C>) {
          binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelected), c_gen, o);
        } else {
          inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelected), o);
        }

        if (Nrec > 0) {
          if constexpr (has_reco_cent<C>) {
            binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedgt0), c_gen, o);
          } else {
            inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedgt0), o);
          }
          atLeastOne_gt0 = true;
        }

        if constexpr (has_reco_cent<C>) {
          binnedRegistry.fill(HIST(NtrkZvtxGen), Nrec, collision.posZ(), c_rec, o);
          binnedRegistry.fill(HIST(NpvcZvxtGen), groupPVcontrib.size(), collision.posZ(), c_rec, o);
        } else {
          inclusiveRegistry.fill(HIST(NtrkZvtxGen), Nrec, collision.posZ());
          inclusiveRegistry.fill(HIST(NpvcZvxtGen), groupPVcontrib.size(), collision.posZ(), o);
        }
      }
    }

    auto nCharged = countParticles(particles);
    if constexpr (has_reco_cent<C>) {
      binnedRegistry.fill(HIST(NtrkZvtxGen_t), nCharged, mcCollision.posZ(), c_gen);
      binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGen), c_gen, o_max);
    } else {
      inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGen), o_max);
      if (useProcId) {
        if constexpr (has_hepmc_pid<MC>) {
          inclusiveRegistry.fill(HIST(NtrkZvtxGen_t), nCharged, mcCollision.posZ(), mcCollision.processId());
        } else {
          LOGP(fatal, "useProcId = true, but MC collision does not have HepMC info");
        }
      } else {
        inclusiveRegistry.fill(HIST(NtrkZvtxGen_t), nCharged, mcCollision.posZ());
      }
    }

    if (nCharged > 0) {
      if constexpr (has_reco_cent<C>) {
        binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGengt0), c_gen, o_max);
      } else {
        inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGengt0), o_max);
      }
    }

    if (fillResponse) {
      for (auto i = 0U; i < NrecPerCol.size(); ++i) {
        if constexpr (has_reco_cent<C>) {
          if (useProcId) {
            if constexpr (has_hepmc_pid<MC>) {
              binnedRegistry.fill(HIST(EfficiencyMult), nCharged, mcCollision.posZ(), c_recPerCol[i], mcCollision.processId());
            } else {
              LOGP(fatal, "useProcId = true, but MC collision does not have HepMC info");
            }
          } else {
            binnedRegistry.fill(HIST(EfficiencyMult), nCharged, mcCollision.posZ(), c_recPerCol[i]);
          }
          if (addFT0 && !addFDD) {
            binnedRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, NFT0APerCol[i], NFT0CPerCol[i], mcCollision.posZ(), c_recPerCol[i], OccuPerCol[i]);
          } else if (addFDD && !addFT0) {
            binnedRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, NFDDAPerCol[i], NFDDCPerCol[i], mcCollision.posZ(), c_recPerCol[i], OccuPerCol[i]);
          } else if (addFT0 && addFDD) {
            binnedRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, NFT0APerCol[i], NFT0CPerCol[i], NFDDAPerCol[i], NFDDCPerCol[i], mcCollision.posZ(), c_recPerCol[i], OccuPerCol[i]);
          } else {
            binnedRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, mcCollision.posZ(), c_recPerCol[i], OccuPerCol[i]);
          }
        } else {
          if (useProcId) {
            if constexpr (has_hepmc_pid<MC>) {
              inclusiveRegistry.fill(HIST(EfficiencyMult), nCharged, mcCollision.posZ(), mcCollision.processId());
            } else {
              LOGP(fatal, "useProcId = true, but MC collision does not have HepMC info");
            }
          } else {
            inclusiveRegistry.fill(HIST(EfficiencyMult), nCharged, mcCollision.posZ());
          }
          if (addFT0 && !addFDD) {
            inclusiveRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, NFT0APerCol[i], NFT0CPerCol[i], mcCollision.posZ(), OccuPerCol[i]);
          } else if (addFDD && !addFT0) {
            inclusiveRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, NFDDAPerCol[i], NFDDCPerCol[i], mcCollision.posZ(), OccuPerCol[i]);
          } else if (addFT0 && addFDD) {
            inclusiveRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, NFT0APerCol[i], NFT0CPerCol[i], NFDDAPerCol[i], NFDDCPerCol[i], mcCollision.posZ(), OccuPerCol[i]);
          } else {
            inclusiveRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, mcCollision.posZ(), OccuPerCol[i]);
          }
        }
      }
      if (moreThanOne > 1) {
        if constexpr (has_reco_cent<C>) {
          binnedRegistry.fill(HIST(SplitMult), nCharged, mcCollision.posZ(), c_gen);
        } else {
          inclusiveRegistry.fill(HIST(SplitMult), nCharged, mcCollision.posZ());
        }
      }
    }

    if (collisions.size() == 0) {
      if constexpr (has_reco_cent<C>) {
        binnedRegistry.fill(HIST(NotFoundZvtx), mcCollision.posZ(), c_gen);
      } else {
        inclusiveRegistry.fill(HIST(NotFoundZvtx), mcCollision.posZ());
      }
    }

    auto zmc = mcCollision.posZ();
    fillParticleHistos<has_reco_cent<C>>(particles, zmc, nCharged, c_gen, atLeastOne, atLeastOne_gt0, atLeastOne_PVgt0);
  }

  template <typename MC, typename C>
  void processGenGeneral(
    typename MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks)
  {
    float c_gen = -1;
    if constexpr (has_Centrality<MC>) {
      c_gen = getSimCent(mcCollision);
    } else {
      if constexpr (has_FT0C<C>) {
        c_gen = getGenCentFT0C(mcCollision);
      } else if (has_FT0M<C>) {
        c_gen = getGenCentFT0M(mcCollision);
      }
    }

    bool atLeastOne = false;
    bool atLeastOne_gt0 = false;
    bool atLeastOne_PVgt0 = false;
    auto moreThanOne = 0;
    LOGP(debug, "MC col {} has {} reco cols", mcCollision.globalIndex(), collisions.size());

    NrecPerCol.clear();
    c_recPerCol.clear();
    OccuPerCol.clear();
    NPVPerCol.clear();
    NFT0APerCol.clear();
    NFT0CPerCol.clear();
    NFDDAPerCol.clear();
    NFDDCPerCol.clear();

    if constexpr (has_reco_cent<C>) {
      float min_c_rec = 2e2;
      for (auto& collision : collisions) {
        if (!useEvSel || isCollisionSelected(collision)) {
          float c = getRecoCent(collision);
          if (c < min_c_rec) {
            min_c_rec = c;
          }
        }
      }
      if constexpr (!has_Centrality<MC>) {
        if (c_gen < 0) {
          c_gen = min_c_rec; // if there is no generator centrality info, fall back to reco (from the largest reco collision)
        }
      }
    }
    float o_max = -1;
    for (auto& collision : collisions) {
      if (!useEvSel || isCollisionSelected(collision)) {
        auto o = collision.trackOccupancyInTimeRange();
        if (o > o_max) {
          o_max = o;
        }
      }
    }

    for (auto& collision : collisions) {
      usedTracksIds.clear();
      float c_rec = getRecoCent(collision);
      float o = collision.trackOccupancyInTimeRange();
      if constexpr (has_reco_cent<C>) {
        binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kRec), c_gen, o);
      } else {
        inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kRec), o);
      }
      if (!useEvSel || isCollisionSelected(collision)) {
        c_recPerCol.emplace_back(c_rec);
        auto z = collision.posZ();
        ++moreThanOne;
        atLeastOne = true;

        auto groupPVcontrib = pvContribTracksIUEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
        if (groupPVcontrib.size() > 0) {
          if constexpr (has_reco_cent<C>) {
            binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedPVgt0), c_gen, o);
          } else {
            inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedPVgt0), o);
          }
          atLeastOne_PVgt0 = true;
        }

        auto perCollisionSample = tracks.sliceBy(perCol, collision.globalIndex());
        auto Nrec = countTracks<C, false>(perCollisionSample, z, c_rec, o);
        NrecPerCol.emplace_back(Nrec);
        NPVPerCol.emplace_back(groupPVcontrib.size());
        fillFIT(collision, NFT0APerCol, NFT0CPerCol, NFDDAPerCol, NFDDCPerCol);

        if constexpr (has_reco_cent<C>) {
          binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelected), c_gen, o);
        } else {
          inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelected), o);
        }

        if (Nrec > 0) {
          if constexpr (has_reco_cent<C>) {
            binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedgt0), c_gen, o);
          } else {
            inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedgt0), o);
          }
          atLeastOne_gt0 = true;
        }

        if constexpr (has_reco_cent<C>) {
          binnedRegistry.fill(HIST(NtrkZvtxGen), Nrec, collision.posZ(), c_rec, o);
          binnedRegistry.fill(HIST(NpvcZvxtGen), groupPVcontrib.size(), collision.posZ(), c_rec, o);
        } else {
          inclusiveRegistry.fill(HIST(NtrkZvtxGen), Nrec, collision.posZ(), o);
          inclusiveRegistry.fill(HIST(NpvcZvxtGen), groupPVcontrib.size(), collision.posZ(), o);
        }
      }
    }

    auto nCharged = countParticles(particles);
    if constexpr (has_reco_cent<C>) {
      binnedRegistry.fill(HIST(NtrkZvtxGen_t), nCharged, mcCollision.posZ(), c_gen);
      binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGen), c_gen, o_max);
    } else {
      inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGen), o_max);
      if (useProcId) {
        if constexpr (has_hepmc_pid<MC>) {
          inclusiveRegistry.fill(HIST(NtrkZvtxGen_t), nCharged, mcCollision.posZ(), mcCollision.processId());
        } else {
          LOGP(fatal, "useProcId = true, but MC collision does not have HepMC info");
        }
      } else {
        inclusiveRegistry.fill(HIST(NtrkZvtxGen_t), nCharged, mcCollision.posZ());
      }
    }

    if (nCharged > 0) {
      if constexpr (has_reco_cent<C>) {
        binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGengt0), c_gen, o_max);
      } else {
        inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGengt0), o_max);
      }
    }

    if (fillResponse) {
      for (auto i = 0U; i < NrecPerCol.size(); ++i) {
        if constexpr (has_reco_cent<C>) {
          if (useProcId) {
            if constexpr (has_hepmc_pid<MC>) {
              binnedRegistry.fill(HIST(EfficiencyMult), nCharged, mcCollision.posZ(), c_recPerCol[i], mcCollision.processId());
            } else {
              LOGP(fatal, "useProcId = true, but MC collision does not have HepMC info");
            }
          } else {
            binnedRegistry.fill(HIST(EfficiencyMult), nCharged, mcCollision.posZ(), c_recPerCol[i]);
          }
          if (addFT0 && !addFDD) {
            binnedRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, NFT0APerCol[i], NFT0CPerCol[i], mcCollision.posZ(), c_recPerCol[i], OccuPerCol[i]);
          } else if (addFDD && !addFT0) {
            binnedRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, NFDDAPerCol[i], NFDDCPerCol[i], mcCollision.posZ(), c_recPerCol[i], OccuPerCol[i]);
          } else if (addFT0 && addFDD) {
            binnedRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, NFT0APerCol[i], NFT0CPerCol[i], NFDDAPerCol[i], NFDDCPerCol[i], mcCollision.posZ(), c_recPerCol[i], OccuPerCol[i]);
          } else {
            binnedRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, mcCollision.posZ(), c_recPerCol[i], OccuPerCol[i]);
          }
        } else {
          if (useProcId) {
            if constexpr (has_hepmc_pid<MC>) {
              inclusiveRegistry.fill(HIST(EfficiencyMult), nCharged, mcCollision.posZ(), mcCollision.processId());
            } else {
              LOGP(fatal, "useProcId = true, but MC collision does not have HepMC info");
            }
          } else {
            inclusiveRegistry.fill(HIST(EfficiencyMult), nCharged, mcCollision.posZ());
          }
          if (addFT0 && !addFDD) {
            inclusiveRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, NFT0APerCol[i], NFT0CPerCol[i], mcCollision.posZ(), OccuPerCol[i]);
          } else if (addFDD && !addFT0) {
            inclusiveRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, NFDDAPerCol[i], NFDDCPerCol[i], mcCollision.posZ(), OccuPerCol[i]);
          } else if (addFT0 && addFDD) {
            inclusiveRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, NFT0APerCol[i], NFT0CPerCol[i], NFDDAPerCol[i], NFDDCPerCol[i], mcCollision.posZ(), OccuPerCol[i]);
          } else {
            inclusiveRegistry.fill(HIST(Response), NrecPerCol[i], NPVPerCol[i], nCharged, mcCollision.posZ(), OccuPerCol[i]);
          }
        }
      }
      if (moreThanOne > 1) {
        if constexpr (has_reco_cent<C>) {
          binnedRegistry.fill(HIST(SplitMult), nCharged, mcCollision.posZ(), c_gen);
        } else {
          inclusiveRegistry.fill(HIST(SplitMult), nCharged, mcCollision.posZ());
        }
      }
    }

    if (collisions.size() == 0) {
      if constexpr (has_reco_cent<C>) {
        binnedRegistry.fill(HIST(NotFoundZvtx), mcCollision.posZ(), c_gen);
      } else {
        inclusiveRegistry.fill(HIST(NotFoundZvtx), mcCollision.posZ());
      }
    }
    auto zmc = mcCollision.posZ();
    fillParticleHistos<has_reco_cent<C>>(particles, zmc, nCharged, c_gen, atLeastOne, atLeastOne_gt0, atLeastOne_PVgt0);
  }

  using MC = aod::McCollisions; // soa::Join<aod::McCollisions, aod::HepMCXSections>;
  void processGenAmbiguous(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExCols, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneralAmbiguous<MC, ExCols>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenAmbiguous, "Process generator-level info", false);

  void processGen(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExCols, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MC, ExCols>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGen, "Process generator-level info w/o ambiguous", false);

  using MCex = soa::Join<aod::McCollisions, aod::HepMCXSections>;
  void processGenAmbiguousEx(
    MCex::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExCols, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneralAmbiguous<MCex, ExCols>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenAmbiguousEx, "Process generator-level info", false);

  void processGenEx(
    MCex::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExCols, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MCex, ExCols>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenEx, "Process generator-level info w/o ambiguous", false);

  void processGenAmbiguousFT0C(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneralAmbiguous<MC, ExColsCentFT0C>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenAmbiguousFT0C, "Process generator-level info (FT0C centrality)", false);

  void processGenFT0C(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MC, ExColsCentFT0C>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0C, "Process generator-level info (FT0C centrality) w/o ambiguous", false);

  void processGenAmbiguousFT0Cplus(
    soa::Join<MC, aod::GenCents>::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneralAmbiguous<soa::Join<MC, aod::GenCents>, ExColsCentFT0C>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenAmbiguousFT0Cplus, "Process generator-level info (FT0C centrality + gen level)", false);

  void processGenFT0Cplus(
    soa::Join<MC, aod::GenCents>::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<soa::Join<MC, aod::GenCents>, ExColsCentFT0C>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0Cplus, "Process generator-level info (FT0C centrality + gen level) w/o ambiguous", false);

  void processGenAmbiguousExFT0C(
    MCex::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneralAmbiguous<MCex, ExColsCentFT0C>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenAmbiguousExFT0C, "Process generator-level info (FT0C centrality)", false);

  void processGenExFT0C(
    MCex::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MCex, ExColsCentFT0C>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenExFT0C, "Process generator-level info (FT0C centrality) w/o ambiguous", false);

  void processGenAmbiguousFT0M(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneralAmbiguous<MC, ExColsCentFT0M>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenAmbiguousFT0M, "Process generator-level info (FT0M centrality)", false);

  void processGenFT0M(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MC, ExColsCentFT0M>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0M, "Process generator-level info (FT0M centrality) w/o ambiguous", false);

  void processGenAmbiguousFT0Mplus(
    soa::Join<MC, aod::GenCents>::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneralAmbiguous<soa::Join<MC, aod::GenCents>, ExColsCentFT0M>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenAmbiguousFT0Mplus, "Process generator-level info (FT0M centrality + gen level)", false);

  void processGenFT0Mplus(
    soa::Join<MC, aod::GenCents>::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<soa::Join<MC, aod::GenCents>, ExColsCentFT0M>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0Mplus, "Process generator-level info (FT0M centrality + gen level) w/o ambiguous", false);

  void processGenAmbiguousExFT0M(
    MCex::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneralAmbiguous<MCex, ExColsCentFT0M>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenAmbiguousExFT0M, "Process generator-level info (FT0M centrality)", false);

  void processGenExFT0M(
    MCex::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MCex, ExColsCentFT0M>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenExFT0M, "Process generator-level info (FT0M centrality) w/o ambiguous", false);

  using MChi = soa::Join<aod::McCollisions, aod::HepMCHeavyIons>;

  void processGenAmbiguousFT0Chi(
    MChi::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneralAmbiguous<MChi, ExColsCentFT0C>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenAmbiguousFT0Chi, "Process generator-level info (FT0C centrality, HI)", false);

  void processGenFT0Chi(
    MChi::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MChi, ExColsCentFT0C>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0Chi, "Process generator-level info (FT0C centrality, HI) w/o ambiguous", false);

  void processGenAmbiguousFT0Mhi(
    MChi::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneralAmbiguous<MChi, ExColsCentFT0M>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenAmbiguousFT0Mhi, "Process generator-level info (FT0M centrality, HI)", false);

  void processGenFT0Mhi(
    MChi::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MChi, ExColsCentFT0M>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0Mhi, "Process generator-level info (FT0M centrality, HI) w/o ambiguous", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityCounter>(cfgc)};
}
