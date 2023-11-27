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
  Configurable<bool> useEvSel{"useEvSel", true, "use event selection"};
  Configurable<bool> fillResponse{"fillResponse", false, "Fill response matrix"};
  Configurable<bool> responseStudy{"responseStudy", false, "Fill multi-estimator response"};
  ConfigurableAxis multBinning{"multBinning", {301, -0.5, 300.5}, ""};
  ConfigurableAxis centBinning{"centBinning", {VARIABLE_WIDTH, 0, 10, 20, 30, 40, 50, 60, 70, 80, 100}, ""};

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
    {
      auto hstat = commonRegistry.get<TH1>(HIST(BCSelection));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, "Good BCs");
      x->SetBinLabel(2, "BCs with collisions");
      x->SetBinLabel(3, "BCs with pile-up/splitting");
    }

    if (doprocessEventStat) {
      inclusiveRegistry.add({EventChi2.data(), " ; #chi^2", {HistType::kTH1F, {{101, -0.1, 10.1}}}});
      inclusiveRegistry.add({EventTimeRes.data(), " ; t (ms)", {HistType::kTH1F, {{1001, -0.1, 100.1}}}});
    }
    if (doprocessEventStatCentralityFT0C || doprocessEventStatCentralityFT0M) {
      binnedRegistry.add({EventChi2.data(), " ; #chi^2; centrality", {HistType::kTH2F, {{101, -0.1, 10.1}, CentAxis}}});
      binnedRegistry.add({EventTimeRes.data(), " ; t (ms); centrality", {HistType::kTH2F, {{1001, -0.1, 100.1}, CentAxis}}});
    }

    if (doprocessCounting || doprocessCountingNoAmb) {
      inclusiveRegistry.add({EventSelection.data(), ";status;events", {HistType::kTH1F, {{static_cast<int>(EvSelBins::kRejected), 0.5, static_cast<float>(EvSelBins::kRejected) + 0.5}}}});
      auto hstat = inclusiveRegistry.get<TH1>(HIST(EventSelection));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(static_cast<int>(EvSelBins::kAll), EvSelBinLabels[static_cast<int>(EvSelBins::kAll)].data());
      x->SetBinLabel(static_cast<int>(EvSelBins::kSelected), EvSelBinLabels[static_cast<int>(EvSelBins::kSelected)].data());
      x->SetBinLabel(static_cast<int>(EvSelBins::kSelectedgt0), EvSelBinLabels[static_cast<int>(EvSelBins::kSelectedgt0)].data());
      x->SetBinLabel(static_cast<int>(EvSelBins::kSelectedPVgt0), EvSelBinLabels[static_cast<int>(EvSelBins::kSelectedPVgt0)].data());
      x->SetBinLabel(static_cast<int>(EvSelBins::kRejected), EvSelBinLabels[static_cast<int>(EvSelBins::kRejected)].data());

      inclusiveRegistry.add({NtrkZvtx.data(), "; N_{trk}; Z_{vtx} (cm); events", {HistType::kTH2F, {MultAxis, ZAxis}}});
      inclusiveRegistry.add({EtaZvtx.data(), "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({EtaZvtx_gt0.data(), "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({EtaZvtx_PVgt0.data(), "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({PhiEta.data(), "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      inclusiveRegistry.add({PtEta.data(), " ; p_{T} (GeV/c); #eta", {HistType::kTH2F, {PtAxis, EtaAxis}}});
      inclusiveRegistry.add({DCAXYPt.data(), " ; p_{T} (GeV/c) ; DCA_{XY} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
      inclusiveRegistry.add({DCAZPt.data(), " ; p_{T} (GeV/c) ; DCA_{Z} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
      if (doprocessCounting) {
        inclusiveRegistry.add({ReassignedDCAXYPt.data(), " ; p_{T} (GeV/c) ; DCA_{XY} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
        inclusiveRegistry.add({ReassignedDCAZPt.data(), " ; p_{T} (GeV/c) ; DCA_{Z} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
        inclusiveRegistry.add({ExtraDCAXYPt.data(), " ; p_{T} (GeV/c) ; DCA_{XY} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
        inclusiveRegistry.add({ExtraDCAZPt.data(), " ; p_{T} (GeV/c) ; DCA_{Z} (cm)", {HistType::kTH2F, {PtAxis, DCAAxis}}});
        inclusiveRegistry.add({ExtraEtaZvtx.data(), "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
        inclusiveRegistry.add({ExtraPhiEta.data(), "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
        inclusiveRegistry.add({ReassignedEtaZvtx.data(), "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
        inclusiveRegistry.add({ReassignedPhiEta.data(), "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
        inclusiveRegistry.add({ReassignedZvtxCorr.data(), "; Z_{vtx}^{orig} (cm); Z_{vtx}^{re} (cm)", {HistType::kTH2F, {ZAxis, ZAxis}}});
      }
    }

    if (doprocessCountingCentralityFT0C || doprocessCountingCentralityFT0M || doprocessCountingCentralityFT0CNoAmb || doprocessCountingCentralityFT0MNoAmb) {
      binnedRegistry.add({EventSelection.data(), ";status;centrality;events", {HistType::kTH2F, {{3, 0.5, 3.5}, CentAxis}}});
      auto hstat = binnedRegistry.get<TH2>(HIST(EventSelection));
      auto* x = hstat->GetXaxis();
      x->SetBinLabel(1, EvSelBinLabels[static_cast<int>(EvSelBins::kAll)].data());
      x->SetBinLabel(2, EvSelBinLabels[static_cast<int>(EvSelBins::kSelected)].data());
      x->SetBinLabel(3, EvSelBinLabels[static_cast<int>(EvSelBins::kRejected)].data());

      binnedRegistry.add({NtrkZvtx.data(), "; N_{trk}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({EtaZvtx.data(), "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({PhiEta.data(), "; #varphi; #eta; centrality", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
      binnedRegistry.add({PtEta.data(), " ; p_{T} (GeV/c); #eta; centrality", {HistType::kTHnSparseF, {PtAxis, EtaAxis, CentAxis}}});
      binnedRegistry.add({DCAXYPt.data(), " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
      binnedRegistry.add({DCAZPt.data(), " ; p_{T} (GeV/c) ; DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
      if (doprocessCountingCentralityFT0C || doprocessCountingCentralityFT0M) {
        binnedRegistry.add({ReassignedDCAXYPt.data(), " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
        binnedRegistry.add({ReassignedDCAZPt.data(), " ; p_{T} (GeV/c) ; DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
        binnedRegistry.add({ExtraDCAXYPt.data(), " ; p_{T} (GeV/c) ; DCA_{XY} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
        binnedRegistry.add({ExtraDCAZPt.data(), " ; p_{T} (GeV/c) ; DCA_{Z} (cm); centrality", {HistType::kTHnSparseF, {PtAxis, DCAAxis, CentAxis}}});
        binnedRegistry.add({ExtraEtaZvtx.data(), "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
        binnedRegistry.add({ExtraPhiEta.data(), "; #varphi; #eta; centrality", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
        binnedRegistry.add({ReassignedEtaZvtx.data(), "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
        binnedRegistry.add({ReassignedPhiEta.data(), "; #varphi; #eta; centrality", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
        binnedRegistry.add({ReassignedZvtxCorr.data(), "; Z_{vtx}^{orig} (cm); Z_{vtx}^{re} (cm); centrality", {HistType::kTHnSparseF, {ZAxis, ZAxis, CentAxis}}});
      }
    }

    if (doprocessGen || doprocessGenNoAmb) {
      inclusiveRegistry.add({NtrkZvtxGen.data(), "; N_{trk}; Z_{vtx} (cm); events", {HistType::kTH2F, {MultAxis, ZAxis}}});
      inclusiveRegistry.add({NtrkZvtxGen_t.data(), "; N_{part}; Z_{vtx} (cm); events", {HistType::kTH2F, {MultAxis, ZAxis}}});
      inclusiveRegistry.add({EtaZvtxGen.data(), "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({EtaZvtxGen_t.data(), "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({EtaZvtxGen_gt0.data(), "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({EtaZvtxGen_PVgt0.data(), "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({EtaZvtxGen_gt0t.data(), "; #eta; Z_{vtx} (cm); tracks", {HistType::kTH2F, {EtaAxis, ZAxis}}});
      inclusiveRegistry.add({PtEtaGen.data(), " ; p_{T} (GeV/c) ; #eta", {HistType::kTH2F, {PtAxis, EtaAxis}}});

      inclusiveRegistry.add({PhiEtaGen.data(), "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});

      inclusiveRegistry.add({Efficiency.data(), "; status; events", {HistType::kTH1F, {{static_cast<int>(EvEffBins::kSelectedPVgt0), 0.5, static_cast<float>(EvEffBins::kSelectedPVgt0) + 0.5}}}});
      inclusiveRegistry.add({NotFoundZvtx.data(), " ; Z_{vtx} (cm)", {HistType::kTH1F, {ZAxis}}});

      if (fillResponse) {
        inclusiveRegistry.add({Response.data(), " ; N_{rec}; N_{gen}; Z_{vtx} (cm)", {HistType::kTHnSparseF, {MultAxis, MultAxis, ZAxis}}});
        inclusiveRegistry.add({EfficiencyMult.data(), " ; N_{gen}; Z_{vtx} (cm)", {HistType::kTH2F, {MultAxis, ZAxis}}});
        inclusiveRegistry.add({SplitMult.data(), " ; N_{gen} ; Z_{vtx} (cm)", {HistType::kTH2F, {MultAxis, ZAxis}}});
        if (responseStudy) {
          inclusiveRegistry.add({MultiResponse.data(), " ; N_{gen}; N_{rec}; N_{PV cont}; N_{FT0A}; N_{FT0C}; N_{FDA}; N_{FDC}; Z_{vtx} (cm)", {HistType::kTHnSparseF, {MultAxis, MultAxis, MultAxis, FT0AAxis, FT0CAxis, FDDAxis, FDDAxis, ZAxis}}});
        }
      }

      auto heff = inclusiveRegistry.get<TH1>(HIST(Efficiency));
      auto* x = heff->GetXaxis();
      x->SetBinLabel(static_cast<int>(EvEffBins::kGen), EvEffBinLabels[static_cast<int>(EvEffBins::kGen)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kGengt0), EvEffBinLabels[static_cast<int>(EvEffBins::kGengt0)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kRec), EvEffBinLabels[static_cast<int>(EvEffBins::kRec)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelected), EvEffBinLabels[static_cast<int>(EvEffBins::kSelected)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelectedgt0), EvEffBinLabels[static_cast<int>(EvEffBins::kSelectedgt0)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelectedPVgt0), EvEffBinLabels[static_cast<int>(EvEffBins::kSelectedPVgt0)].data());
    }

    if (doprocessGenFT0C || doprocessGenFT0M || doprocessGenFT0Chi || doprocessGenFT0Mhi ||
        doprocessGenFT0CNoAmb || doprocessGenFT0MNoAmb || doprocessGenFT0ChiNoAmb || doprocessGenFT0MhiNoAmb) {
      binnedRegistry.add({NtrkZvtxGen.data(), "; N_{trk}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({NtrkZvtxGen_t.data(), "; N_{part}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({EtaZvtxGen.data(), "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({EtaZvtxGen_t.data(), "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({EtaZvtxGen_gt0.data(), "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({EtaZvtxGen_PVgt0.data(), "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({EtaZvtxGen_gt0t.data(), "; #eta; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {EtaAxis, ZAxis, CentAxis}}});
      binnedRegistry.add({PtEtaGen.data(), " ; p_{T} (GeV/c) ; #eta; centrality", {HistType::kTHnSparseF, {PtAxis, EtaAxis, CentAxis}}});

      binnedRegistry.add({PhiEtaGen.data(), "; #varphi; #eta; tracks", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
      binnedRegistry.add({PhiEtaGenDuplicates.data(), "; #varphi; #eta; centrality", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
      binnedRegistry.add({PhiEtaDuplicates.data(), "; #varphi; #eta; centrality", {HistType::kTHnSparseF, {PhiAxis, EtaAxis, CentAxis}}});
      binnedRegistry.add({Efficiency.data(), "; status; centrality; events", {HistType::kTH2F, {{static_cast<int>(EvEffBins::kSelectedPVgt0), 0.5, static_cast<float>(EvEffBins::kSelectedPVgt0) + 0.5}, CentAxis}}});
      binnedRegistry.add({NotFoundZvtx.data(), " ; Z_{vtx} (cm); centrality; events", {HistType::kTH2F, {ZAxis, CentAxis}}});

      if (fillResponse) {
        binnedRegistry.add({Response.data(), " ; N_{rec}; N_{gen}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, MultAxis, ZAxis, CentAxis}}});
        binnedRegistry.add({EfficiencyMult.data(), " ; N_{gen}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});
        binnedRegistry.add({SplitMult.data(), " ; N_{gen} ; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, ZAxis, CentAxis}}});
        if (responseStudy) {
          binnedRegistry.add({MultiResponse.data(), " ; N_{gen}; N_{rec}, N_{PV cont}; N_{FT0A}; N_{FT0C}; N_{FDA}; N_{FDC}; Z_{vtx} (cm); centrality", {HistType::kTHnSparseF, {MultAxis, MultAxis, MultAxis, FT0AAxis, FT0CAxis, FDDAxis, FDDAxis, ZAxis, CentAxis}}});
        }
      }

      auto heff = binnedRegistry.get<TH2>(HIST(Efficiency));
      auto* x = heff->GetXaxis();
      x->SetBinLabel(static_cast<int>(EvEffBins::kGen), EvEffBinLabels[static_cast<int>(EvEffBins::kGen)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kGengt0), EvEffBinLabels[static_cast<int>(EvEffBins::kGengt0)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kRec), EvEffBinLabels[static_cast<int>(EvEffBins::kRec)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelected), EvEffBinLabels[static_cast<int>(EvEffBins::kSelected)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelectedgt0), EvEffBinLabels[static_cast<int>(EvEffBins::kSelectedgt0)].data());
      x->SetBinLabel(static_cast<int>(EvEffBins::kSelectedPVgt0), EvEffBinLabels[static_cast<int>(EvEffBins::kSelectedPVgt0)].data());
    }

    if (doprocessTrackEfficiency || doprocessTrackEfficiencyNoAmb) {
      inclusiveRegistry.add({PtGen.data(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({PtGenNoEtaCut.data(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({PtEfficiency.data(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({PtEfficiencyNoEtaCut.data(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({PtEfficiencyFakes.data(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      for (auto i = 0u; i < speciesIds.size(); ++i) {
        inclusiveRegistry.add({fmt::format(PtGenF.data(), species[i]).c_str(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
        inclusiveRegistry.add({fmt::format(PtEfficiencyF.data(), species[i]).c_str(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      }
    }
    if (doprocessTrackEfficiencyIndexed) {
      inclusiveRegistry.add({PhiEtaGenDuplicates.data(), "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      inclusiveRegistry.add({PhiEtaDuplicates.data(), "; #varphi; #eta; tracks", {HistType::kTH2F, {PhiAxis, EtaAxis}}});
      inclusiveRegistry.add({PtGenIdx.data(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({PtGenIdxNoEtaCut.data(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({PtEfficiencyIdx.data(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({PtEfficiencyIdxNoEtaCut.data(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({PtEfficiencySecondariesIdx.data(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({PtEfficiencySecondariesIdxNoEtaCut.data(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      inclusiveRegistry.add({Mask.data(), " ; bit", {HistType::kTH1F, {{17, -0.5, 16.5}}}});
      inclusiveRegistry.add({ITSlayers.data(), " ; layer", {HistType::kTH1F, {{8, 0.5, 8.5}}}});
      for (auto i = 0u; i < speciesIds.size(); ++i) {
        inclusiveRegistry.add({fmt::format(PtGenIdxF.data(), species[i]).c_str(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
        inclusiveRegistry.add({fmt::format(PtEfficiencyIdxF.data(), species[i]).c_str(), " ; p_{T} (GeV/c)", {HistType::kTH1F, {PtAxisEff}}});
      }
    }
  }

  using FullBCs = soa::Join<aod::BCsWithTimestamps, aod::BcSels>;
  template <typename C>
  void processEventStatGeneral(FullBCs const& bcs, C const& collisions)
  {
    std::vector<typename std::decay_t<decltype(collisions)>::iterator> cols;
    for (auto& bc : bcs) {
      if (!useEvSel || (bc.selection_bit(aod::evsel::kNoITSROFrameBorder) &&
                        bc.selection_bit(aod::evsel::kIsBBT0A) &&
                        bc.selection_bit(aod::evsel::kIsBBT0C)) != 0) {
        commonRegistry.fill(HIST(BCSelection), 1.);
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
          commonRegistry.fill(HIST(BCSelection), 2.);
          if (cols.size() > 1) {
            commonRegistry.fill(HIST(BCSelection), 3.);
          }
        }
        for (auto& col : cols) {
          if constexpr (hasRecoCent<C>()) {
            float c = -1;
            if constexpr (C::template contains<aod::CentFT0Cs>()) {
              c = col.centFT0C();
            } else if constexpr (C::template contains<aod::CentFT0Ms>()) {
              c = col.centFT0M();
            }
            binnedRegistry.fill(HIST(EventChi2), col.chi2(), c);
            binnedRegistry.fill(HIST(EventTimeRes), col.collisionTimeRes(), c);
          } else {
            inclusiveRegistry.fill(HIST(EventChi2), col.chi2());
            inclusiveRegistry.fill(HIST(EventTimeRes), col.collisionTimeRes());
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

  template <typename C>
  void processCountingGeneral(
    typename C::iterator const& collision,
    FiTracks const& tracks)
  {
    float c = -1;
    if constexpr (hasRecoCent<C>()) {
      if constexpr (C::template contains<aod::CentFT0Cs>()) {
        c = collision.centFT0C();
      } else if (C::template contains<aod::CentFT0Ms>()) {
        c = collision.centFT0M();
      }
      binnedRegistry.fill(HIST(EventSelection), 1., c);
    } else {
      inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kAll));
    }

    if (!useEvSel || (collision.sel8() && collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST(EventSelection), 2., c);
      } else {
        inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelected));
      }
      auto z = collision.posZ();
      usedTracksIds.clear();

      auto groupPVContrib = pvContribTracksIUEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      auto Ntrks = 0;
      for (auto& track : tracks) {
        if (std::abs(track.eta()) < estimatorEta) {
          ++Ntrks;
        }
        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST(EtaZvtx), track.eta(), z, c);
          binnedRegistry.fill(HIST(PhiEta), track.phi(), track.eta(), c);
          binnedRegistry.fill(HIST(PtEta), track.pt(), track.eta(), c);
          binnedRegistry.fill(HIST(DCAXYPt), track.pt(), track.dcaXY(), c);
          binnedRegistry.fill(HIST(DCAZPt), track.pt(), track.dcaZ(), c);
        } else {
          inclusiveRegistry.fill(HIST(EtaZvtx), track.eta(), z);
          inclusiveRegistry.fill(HIST(PhiEta), track.phi(), track.eta());
          inclusiveRegistry.fill(HIST(PtEta), track.pt(), track.eta());
          inclusiveRegistry.fill(HIST(DCAXYPt), track.pt(), track.dcaXY());
          inclusiveRegistry.fill(HIST(DCAZPt), track.pt(), track.dcaZ());
        }
      }
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST(NtrkZvtx), Ntrks, z, c);
      } else {
        if (Ntrks > 0 || groupPVContrib.size() > 0) {
          if (groupPVContrib.size() > 0) {
            inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelectedPVgt0));
          }
          if (Ntrks > 0) {
            inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelectedgt0));
          }
          for (auto& track : tracks) {
            if (Ntrks > 0) {
              inclusiveRegistry.fill(HIST(EtaZvtx_gt0), track.eta(), z);
            }
            if (groupPVContrib.size() > 0) {
              inclusiveRegistry.fill(HIST(EtaZvtx_PVgt0), track.eta(), z);
            }
          }
        }
        inclusiveRegistry.fill(HIST(NtrkZvtx), Ntrks, z);
      }
    } else {
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST(EventSelection), 3., c);
      } else {
        inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kRejected));
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
      binnedRegistry.fill(HIST(EventSelection), 1., c);
    } else {
      inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kAll));
    }

    if (!useEvSel || (collision.sel8() && collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST(EventSelection), 2., c);
      } else {
        inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelected));
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
          binnedRegistry.fill(HIST(EtaZvtx), otrack.eta(), z, c);
          binnedRegistry.fill(HIST(PhiEta), otrack.phi(), otrack.eta(), c);
          binnedRegistry.fill(HIST(PtEta), otrack.pt(), otrack.eta(), c);
          binnedRegistry.fill(HIST(DCAXYPt), otrack.pt(), track.bestDCAXY(), c);
          binnedRegistry.fill(HIST(DCAZPt), otrack.pt(), track.bestDCAZ(), c);
        } else {
          inclusiveRegistry.fill(HIST(EtaZvtx), otrack.eta(), z);
          inclusiveRegistry.fill(HIST(PhiEta), otrack.phi(), otrack.eta());
          inclusiveRegistry.fill(HIST(PtEta), otrack.pt(), otrack.eta());
          inclusiveRegistry.fill(HIST(DCAXYPt), otrack.pt(), track.bestDCAXY());
          inclusiveRegistry.fill(HIST(DCAZPt), otrack.pt(), track.bestDCAZ());
        }
        if (!otrack.has_collision()) {
          if constexpr (hasRecoCent<C>()) {
            binnedRegistry.fill(HIST(ExtraEtaZvtx), otrack.eta(), z, c);
            binnedRegistry.fill(HIST(ExtraPhiEta), otrack.phi(), otrack.eta(), c);
            binnedRegistry.fill(HIST(ExtraDCAXYPt), otrack.pt(), track.bestDCAXY(), c);
            binnedRegistry.fill(HIST(ExtraDCAZPt), otrack.pt(), track.bestDCAZ(), c);
          } else {
            inclusiveRegistry.fill(HIST(ExtraEtaZvtx), otrack.eta(), z);
            inclusiveRegistry.fill(HIST(ExtraPhiEta), otrack.phi(), otrack.eta());
            inclusiveRegistry.fill(HIST(ExtraDCAXYPt), otrack.pt(), track.bestDCAXY());
            inclusiveRegistry.fill(HIST(ExtraDCAZPt), otrack.pt(), track.bestDCAZ());
          }
        } else if (otrack.collisionId() != track.bestCollisionId()) {
          usedTracksIdsDF.emplace_back(track.trackId());
          if constexpr (hasRecoCent<C>()) {
            binnedRegistry.fill(HIST(ReassignedEtaZvtx), otrack.eta(), z, c);
            binnedRegistry.fill(HIST(ReassignedPhiEta), otrack.phi(), otrack.eta(), c);
            binnedRegistry.fill(HIST(ReassignedZvtxCorr), otrack.collision_as<C>().posZ(), z, c);
            binnedRegistry.fill(HIST(ReassignedDCAXYPt), otrack.pt(), track.bestDCAXY(), c);
            binnedRegistry.fill(HIST(ReassignedDCAZPt), otrack.pt(), track.bestDCAZ(), c);
          } else {
            inclusiveRegistry.fill(HIST(ReassignedEtaZvtx), otrack.eta(), z);
            inclusiveRegistry.fill(HIST(ReassignedPhiEta), otrack.phi(), otrack.eta());
            inclusiveRegistry.fill(HIST(ReassignedZvtxCorr), otrack.collision_as<C>().posZ(), z);
            inclusiveRegistry.fill(HIST(ReassignedDCAXYPt), otrack.pt(), track.bestDCAXY());
            inclusiveRegistry.fill(HIST(ReassignedDCAZPt), otrack.pt(), track.bestDCAZ());
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
          binnedRegistry.fill(HIST(EtaZvtx), track.eta(), z, c);
          binnedRegistry.fill(HIST(PhiEta), track.phi(), track.eta(), c);
          binnedRegistry.fill(HIST(PtEta), track.pt(), track.eta(), c);
          binnedRegistry.fill(HIST(DCAXYPt), track.pt(), track.dcaXY(), c);
          binnedRegistry.fill(HIST(DCAZPt), track.pt(), track.dcaZ(), c);
        } else {
          inclusiveRegistry.fill(HIST(EtaZvtx), track.eta(), z);
          inclusiveRegistry.fill(HIST(PhiEta), track.phi(), track.eta());
          inclusiveRegistry.fill(HIST(PtEta), track.pt(), track.eta());
          inclusiveRegistry.fill(HIST(DCAXYPt), track.pt(), track.dcaXY());
          inclusiveRegistry.fill(HIST(DCAZPt), track.pt(), track.dcaZ());
        }
      }
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST(NtrkZvtx), Ntrks, z, c);
      } else {
        if (Ntrks > 0 || groupPVContrib.size() > 0) {
          if (groupPVContrib.size() > 0) {
            inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelectedPVgt0));
          }
          if (Ntrks > 0) {
            inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kSelectedgt0));
          }
          for (auto& track : atracks) {
            if (Ntrks > 0) {
              inclusiveRegistry.fill(HIST(EtaZvtx_gt0), track.track_as<FiTracks>().eta(), z);
            }
            if (groupPVContrib.size() > 0) {
              inclusiveRegistry.fill(HIST(EtaZvtx_PVgt0), track.track_as<FiTracks>().eta(), z);
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
              inclusiveRegistry.fill(HIST(EtaZvtx_gt0), track.eta(), z);
            }
            if (groupPVContrib.size() > 0) {
              inclusiveRegistry.fill(HIST(EtaZvtx_PVgt0), track.eta(), z);
            }
          }
        }
        inclusiveRegistry.fill(HIST(NtrkZvtx), Ntrks, z);
      }
    } else {
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST(EventSelection), 3., c);
      } else {
        inclusiveRegistry.fill(HIST(EventSelection), static_cast<float>(EvSelBins::kRejected));
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
    if (useEvSel && !(collision.sel8() && collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
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
      inclusiveRegistry.fill(HIST(PtGenIdxNoEtaCut), particle.pt());

      if (std::abs(particle.eta()) < estimatorEta) {
        inclusiveRegistry.fill(HIST(PtGenIdx), particle.pt());
        if (particle.pdgCode() == speciesIds[0]) {
          inclusiveRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtGenIdxSuff), particle.pt());
        } else if (particle.pdgCode() == speciesIds[1]) {
          inclusiveRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtGenIdxSuff), particle.pt());
        } else if (particle.pdgCode() == speciesIds[2]) {
          inclusiveRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtGenIdxSuff), particle.pt());
        } else if (particle.pdgCode() == speciesIds[3]) {
          inclusiveRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtGenIdxSuff), particle.pt());
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
            inclusiveRegistry.fill(HIST(PtEfficiencyIdxNoEtaCut), particle.pt());
            countedNoEtaCut = true;
          }
          if (std::abs(track.eta()) < estimatorEta) {
            if (!counted) {
              inclusiveRegistry.fill(HIST(PtEfficiencyIdx), particle.pt());
              if (particle.pdgCode() == speciesIds[0]) {
                inclusiveRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtEffIdxSuff), particle.pt());
              } else if (particle.pdgCode() == speciesIds[1]) {
                inclusiveRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtEffIdxSuff), particle.pt());
              } else if (particle.pdgCode() == speciesIds[2]) {
                inclusiveRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtEffIdxSuff), particle.pt());
              } else if (particle.pdgCode() == speciesIds[3]) {
                inclusiveRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtEffIdxSuff), particle.pt());
              }
              counted = true;
            }
          }
          if (counter > 1) {
            inclusiveRegistry.fill(HIST(PtEfficiencySecondariesIdxNoEtaCut), particle.pt());
            if (std::abs(track.eta()) < estimatorEta) {
              inclusiveRegistry.fill(HIST(PtEfficiencySecondariesIdx), particle.pt());
            }
          }
        }
        if (counter > 1) {
          for (auto const& track : relatedTracks) {
            for (auto layer = 0; layer < 7; ++layer) {
              if (track.itsClusterMap() & (uint8_t(1) << layer)) {
                inclusiveRegistry.fill(HIST(ITSlayers), layer + 1);
              }
            }
            auto hasbit = false;
            for (auto bit = 0; bit < 16; ++bit) {
              if (track.mcMask() & (uint8_t(1) << bit)) {
                inclusiveRegistry.fill(HIST(Mask), bit);
                hasbit = true;
              }
            }
            if (!hasbit) {
              inclusiveRegistry.fill(HIST(Mask), 16);
            }
          }
        }
        if (relatedTracks.size() > 1) {
          inclusiveRegistry.fill(HIST(PhiEtaGenDuplicates), particle.phi(), particle.eta());
          for (auto const& track : relatedTracks) {
            inclusiveRegistry.fill(HIST(PhiEtaDuplicates), track.phi(), track.eta());
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
    if (useEvSel && !(collision.sel8() && collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
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
        inclusiveRegistry.fill(HIST(PtEfficiencyNoEtaCut), particle.pt());
        if (std::abs(otrack.eta()) < estimatorEta) {
          inclusiveRegistry.fill(HIST(PtEfficiency), particle.pt());
          if (particle.pdgCode() == speciesIds[0]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtEffSuff), particle.pt());
          } else if (particle.pdgCode() == speciesIds[1]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtEffSuff), particle.pt());
          } else if (particle.pdgCode() == speciesIds[2]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtEffSuff), particle.pt());
          } else if (particle.pdgCode() == speciesIds[3]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtEffSuff), particle.pt());
          }
        }
      } else {
        inclusiveRegistry.fill(HIST(PtEfficiencyFakes), otrack.pt());
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
        inclusiveRegistry.fill(HIST(PtEfficiencyNoEtaCut), particle.pt());
        if (std::abs(track.eta()) < estimatorEta) {
          inclusiveRegistry.fill(HIST(PtEfficiency), particle.pt());
          if (particle.pdgCode() == speciesIds[0]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtEffSuff), particle.pt());
          } else if (particle.pdgCode() == speciesIds[1]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtEffSuff), particle.pt());
          } else if (particle.pdgCode() == speciesIds[2]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtEffSuff), particle.pt());
          } else if (particle.pdgCode() == speciesIds[3]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtEffSuff), particle.pt());
          }
        }
      } else {
        inclusiveRegistry.fill(HIST(PtEfficiencyFakes), track.pt());
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
      inclusiveRegistry.fill(HIST(PtGenNoEtaCut), particle.pt());
      if (std::abs(particle.eta()) < estimatorEta) {
        inclusiveRegistry.fill(HIST(PtGen), particle.pt());
        if (particle.pdgCode() == speciesIds[0]) {
          inclusiveRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtGenSuff), particle.pt());
        } else if (particle.pdgCode() == speciesIds[1]) {
          inclusiveRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtGenSuff), particle.pt());
        } else if (particle.pdgCode() == speciesIds[2]) {
          inclusiveRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtGenSuff), particle.pt());
        } else if (particle.pdgCode() == speciesIds[3]) {
          inclusiveRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtGenSuff), particle.pt());
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
    if (useEvSel && !(collision.sel8() && collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
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
        inclusiveRegistry.fill(HIST(PtEfficiencyNoEtaCut), particle.pt());
        if (std::abs(track.eta()) < estimatorEta) {
          inclusiveRegistry.fill(HIST(PtEfficiency), particle.pt());
          if (particle.pdgCode() == speciesIds[0]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtEffSuff), particle.pt());
          } else if (particle.pdgCode() == speciesIds[1]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtEffSuff), particle.pt());
          } else if (particle.pdgCode() == speciesIds[2]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtEffSuff), particle.pt());
          } else if (particle.pdgCode() == speciesIds[3]) {
            inclusiveRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtEffSuff), particle.pt());
          }
        }
      } else {
        inclusiveRegistry.fill(HIST(PtEfficiencyFakes), track.pt());
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
      inclusiveRegistry.fill(HIST(PtGenNoEtaCut), particle.pt());
      if (std::abs(particle.eta()) < estimatorEta) {
        inclusiveRegistry.fill(HIST(PtGen), particle.pt());
        if (particle.pdgCode() == speciesIds[0]) {
          inclusiveRegistry.fill(HIST(prefix) + HIST(species[0]) + HIST(PtGenSuff), particle.pt());
        } else if (particle.pdgCode() == speciesIds[1]) {
          inclusiveRegistry.fill(HIST(prefix) + HIST(species[1]) + HIST(PtGenSuff), particle.pt());
        } else if (particle.pdgCode() == speciesIds[2]) {
          inclusiveRegistry.fill(HIST(prefix) + HIST(species[2]) + HIST(PtGenSuff), particle.pt());
        } else if (particle.pdgCode() == speciesIds[3]) {
          inclusiveRegistry.fill(HIST(prefix) + HIST(species[3]) + HIST(PtGenSuff), particle.pt());
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

  template <typename CIT>
  void fillFIT(CIT const& collision, std::vector<float>& ft0as, std::vector<float>& ft0cs, std::vector<float>& fddas, std::vector<float>& fddcs)
  {
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
      ft0as.emplace_back(tA);
      ft0cs.emplace_back(tC);
    } else {
      ft0as.emplace_back(-1);
      ft0cs.emplace_back(-1);
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
      fddas.emplace_back(tA);
      fddcs.emplace_back(tC);
    } else {
      fddas.emplace_back(-1);
      fddcs.emplace_back(-1);
    }
  }

  template <bool hasRecoCent, typename Ps, typename MCIT>
  void countParticles(Ps const& particles, MCIT const& mcCollision, int const nCharged, HistogramRegistry& binnedRegistry, HistogramRegistry& inclusiveRegistry, float c_gen,
                      bool const atLeastOne,
                      bool const atLeastOne_gt0,
                      bool const atLeastOne_PVgt0)
  {
    for (auto& particle : particles) {
      auto p = pdg->GetParticle(particle.pdgCode());
      auto charge = 0.;
      if (p != nullptr) {
        charge = p->Charge();
      }
      if (std::abs(charge) < 3.) {
        continue;
      }
      if constexpr (hasRecoCent) {
        binnedRegistry.fill(HIST(EtaZvtxGen_t), particle.eta(), mcCollision.posZ(), c_gen);
        binnedRegistry.fill(HIST(PtEtaGen), particle.pt(), particle.eta(), c_gen);
      } else {
        inclusiveRegistry.fill(HIST(EtaZvtxGen_t), particle.eta(), mcCollision.posZ());
        inclusiveRegistry.fill(HIST(PtEtaGen), particle.pt(), particle.eta());
      }
      if (nCharged > 0) {
        if constexpr (hasRecoCent) {
          binnedRegistry.fill(HIST(EtaZvtxGen_gt0t), particle.eta(), mcCollision.posZ(), c_gen);
        } else {
          inclusiveRegistry.fill(HIST(EtaZvtxGen_gt0t), particle.eta(), mcCollision.posZ());
        }
      }
      if (atLeastOne) {
        if constexpr (hasRecoCent) {
          binnedRegistry.fill(HIST(EtaZvtxGen), particle.eta(), mcCollision.posZ(), c_gen);
          if (atLeastOne_gt0) {
            binnedRegistry.fill(HIST(EtaZvtxGen_gt0), particle.eta(), mcCollision.posZ(), c_gen);
          }
          if (atLeastOne_PVgt0) {
            binnedRegistry.fill(HIST(EtaZvtxGen_PVgt0), particle.eta(), mcCollision.posZ(), c_gen);
          }
          binnedRegistry.fill(HIST(PhiEtaGen), particle.phi(), particle.eta(), c_gen);
        } else {
          inclusiveRegistry.fill(HIST(EtaZvtxGen), particle.eta(), mcCollision.posZ());
          if (atLeastOne_gt0) {
            inclusiveRegistry.fill(HIST(EtaZvtxGen_gt0), particle.eta(), mcCollision.posZ());
          }
          if (atLeastOne_PVgt0) {
            inclusiveRegistry.fill(HIST(EtaZvtxGen_PVgt0), particle.eta(), mcCollision.posZ());
          }
          inclusiveRegistry.fill(HIST(PhiEtaGen), particle.phi(), particle.eta());
        }
      }
    }
  }

  template <typename MC, typename C>
  void processGenGeneralAmbiguous(
    typename MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks)
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
      binnedRegistry.fill(HIST(NtrkZvtxGen_t), nCharged, mcCollision.posZ(), c_gen);
      binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGen), c_gen);
    } else {
      inclusiveRegistry.fill(HIST(NtrkZvtxGen_t), nCharged, mcCollision.posZ());
      inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGen));
    }

    if (nCharged > 0) {
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGengt0), c_gen);
      } else {
        inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGengt0));
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
        binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kRec), c_gen);
      } else {
        inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kRec));
      }
      if (!useEvSel || (collision.sel8() && collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        Nrec = 0;
        ++moreThanOne;
        atLeastOne = true;

        auto groupPVcontrib = pvContribTracksIUEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
        if (groupPVcontrib.size() > 0) {
          if constexpr (hasRecoCent<C>()) {
            binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedPVgt0), c_gen);
          } else {
            inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedPVgt0));
          }
          atLeastOne_PVgt0 = true;
        }

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
        auto perCollisionSample = tracks.sliceBy(perCol, collision.globalIndex());
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
        NPVPerCol.emplace_back(collision.numContrib());
        if (responseStudy) {
          fillFIT(collision, NFT0APerCol, NFT0CPerCol, NFDDAPerCol, NFDDCPerCol);
        }

        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelected), c_gen);
        } else {
          inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelected));
        }

        if (Nrec > 0) {
          if constexpr (hasRecoCent<C>()) {
            binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedgt0), c_gen);
          } else {
            inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedgt0));
          }
          atLeastOne_gt0 = true;
        }

        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST(NtrkZvtxGen), Nrec, collision.posZ(), c_rec);
        } else {
          inclusiveRegistry.fill(HIST(NtrkZvtxGen), Nrec, collision.posZ());
        }
      }
    }

    if (fillResponse) {
      for (auto i = 0U; i < NrecPerCol.size(); ++i) {
        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST(Response), NrecPerCol[i], nCharged, mcCollision.posZ(), c_recPerCol[i]);
          binnedRegistry.fill(HIST(EfficiencyMult), nCharged, mcCollision.posZ(), c_recPerCol[i]);
          if (responseStudy) {
            binnedRegistry.fill(HIST(MultiResponse), nCharged, NrecPerCol[i], NPVPerCol[i], NFT0APerCol[i], NFT0CPerCol[i], NFDDAPerCol[i], NFDDCPerCol[i], mcCollision.posZ(), c_recPerCol[i]);
          }
        } else {
          inclusiveRegistry.fill(HIST(Response), NrecPerCol[i], nCharged, mcCollision.posZ());
          inclusiveRegistry.fill(HIST(EfficiencyMult), nCharged, mcCollision.posZ());
          if (responseStudy) {
            inclusiveRegistry.fill(HIST(MultiResponse), nCharged, NrecPerCol[i], NPVPerCol[i], NFT0APerCol[i], NFT0CPerCol[i], NFDDAPerCol[i], NFDDCPerCol[i], mcCollision.posZ());
          }
        }
      }
      if (moreThanOne > 1) {
        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST(SplitMult), nCharged, mcCollision.posZ(), c_gen);
        } else {
          inclusiveRegistry.fill(HIST(SplitMult), nCharged, mcCollision.posZ());
        }
      }
    }

    if (collisions.size() == 0) {
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST(NotFoundZvtx), mcCollision.posZ(), c_gen);
      } else {
        inclusiveRegistry.fill(HIST(NotFoundZvtx), mcCollision.posZ());
      }
    }

    countParticles<hasRecoCent<C>()>(particles, mcCollision, nCharged, binnedRegistry, inclusiveRegistry, c_gen, atLeastOne, atLeastOne_gt0, atLeastOne_PVgt0);
  }

  template <typename MC, typename C>
  void processGenGeneral(
    typename MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks)
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
      binnedRegistry.fill(HIST(NtrkZvtxGen_t), nCharged, mcCollision.posZ(), c_gen);
      binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGen), c_gen);
    } else {
      inclusiveRegistry.fill(HIST(NtrkZvtxGen_t), nCharged, mcCollision.posZ());
      inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGen));
    }

    if (nCharged > 0) {
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGengt0), c_gen);
      } else {
        inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kGengt0));
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
        binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kRec), c_gen);
      } else {
        inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kRec));
      }
      if (!useEvSel || (collision.sel8() && collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
        Nrec = 0;
        ++moreThanOne;
        atLeastOne = true;

        auto groupPVcontrib = pvContribTracksIUEta1->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
        if (groupPVcontrib.size() > 0) {
          if constexpr (hasRecoCent<C>()) {
            binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedPVgt0), c_gen);
          } else {
            inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedPVgt0));
          }
          atLeastOne_PVgt0 = true;
        }

        auto perCollisionSample = tracks.sliceBy(perCol, collision.globalIndex());
        for (auto const& track : perCollisionSample) {
          if (std::abs(track.eta()) < estimatorEta) {
            ++Nrec;
          }
        }
        NrecPerCol.emplace_back(Nrec);
        NPVPerCol.emplace_back(collision.numContrib());
        if (responseStudy) {
          fillFIT(collision, NFT0APerCol, NFT0CPerCol, NFDDAPerCol, NFDDCPerCol);
        }

        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelected), c_gen);
        } else {
          inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelected));
        }

        if (Nrec > 0) {
          if constexpr (hasRecoCent<C>()) {
            binnedRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedgt0), c_gen);
          } else {
            inclusiveRegistry.fill(HIST(Efficiency), static_cast<float>(EvEffBins::kSelectedgt0));
          }
          atLeastOne_gt0 = true;
        }

        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST(NtrkZvtxGen), Nrec, collision.posZ(), c_rec);
        } else {
          inclusiveRegistry.fill(HIST(NtrkZvtxGen), Nrec, collision.posZ());
        }
      }
    }

    if (fillResponse) {
      for (auto i = 0U; i < NrecPerCol.size(); ++i) {
        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST(Response), NrecPerCol[i], nCharged, mcCollision.posZ(), c_recPerCol[i]);
          binnedRegistry.fill(HIST(EfficiencyMult), nCharged, mcCollision.posZ(), c_recPerCol[i]);
          if (responseStudy) {
            binnedRegistry.fill(HIST(MultiResponse), nCharged, NrecPerCol[i], NPVPerCol[i], NFT0APerCol[i], NFT0CPerCol[i], NFDDAPerCol[i], NFDDCPerCol[i], mcCollision.posZ(), c_recPerCol[i]);
          }
        } else {
          inclusiveRegistry.fill(HIST(Response), NrecPerCol[i], nCharged, mcCollision.posZ());
          inclusiveRegistry.fill(HIST(EfficiencyMult), nCharged, mcCollision.posZ());
          if (responseStudy) {
            inclusiveRegistry.fill(HIST(MultiResponse), nCharged, NrecPerCol[i], NPVPerCol[i], NFT0APerCol[i], NFT0CPerCol[i], NFDDAPerCol[i], NFDDCPerCol[i], mcCollision.posZ());
          }
        }
      }
      if (moreThanOne > 1) {
        if constexpr (hasRecoCent<C>()) {
          binnedRegistry.fill(HIST(SplitMult), nCharged, mcCollision.posZ(), c_gen);
        } else {
          inclusiveRegistry.fill(HIST(SplitMult), nCharged, mcCollision.posZ());
        }
      }
    }

    if (collisions.size() == 0) {
      if constexpr (hasRecoCent<C>()) {
        binnedRegistry.fill(HIST(NotFoundZvtx), mcCollision.posZ(), c_gen);
      } else {
        inclusiveRegistry.fill(HIST(NotFoundZvtx), mcCollision.posZ());
      }
    }

    countParticles<hasRecoCent<C>()>(particles, mcCollision, nCharged, binnedRegistry, inclusiveRegistry, c_gen, atLeastOne, atLeastOne_gt0, atLeastOne_PVgt0);
  }

  using MC = aod::McCollisions; // soa::Join<aod::McCollisions, aod::HepMCXSections>;
  void processGen(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExCols, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneralAmbiguous<MC, ExCols>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGen, "Process generator-level info", false);

  void processGenNoAmb(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExCols, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MC, ExCols>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenNoAmb, "Process generator-level info w/o ambiguous", false);

  void processGenFT0C(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneralAmbiguous<MC, ExColsCentFT0C>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0C, "Process generator-level info (FT0C centrality)", false);

  void processGenFT0CNoAmb(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MC, ExColsCentFT0C>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0CNoAmb, "Process generator-level info (FT0C centrality) w/o ambiguous", false);

  void processGenFT0M(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneralAmbiguous<MC, ExColsCentFT0M>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0M, "Process generator-level info (FT0M centrality)", false);

  void processGenFT0MNoAmb(
    MC::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MC, ExColsCentFT0M>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0MNoAmb, "Process generator-level info (FT0M centrality) w/o ambiguous", false);

  using MChi = soa::Join<aod::McCollisions, aod::HepMCHeavyIons>;

  void processGenFT0Chi(
    MChi::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneralAmbiguous<MChi, ExColsCentFT0C>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0Chi, "Process generator-level info (FT0C centrality, HI)", false);

  void processGenFT0ChiNoAmb(
    MChi::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0C, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MChi, ExColsCentFT0C>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0ChiNoAmb, "Process generator-level info (FT0C centrality, HI) w/o ambiguous", false);

  void processGenFT0Mhi(
    MChi::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, FiReTracks const& atracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneralAmbiguous<MChi, ExColsCentFT0M>(mcCollision, collisions, particles, tracks, atracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0Mhi, "Process generator-level info (FT0M centrality, HI)", false);

  void processGenFT0MhiNoAmb(
    MChi::iterator const& mcCollision,
    o2::soa::SmallGroups<soa::Join<ExColsCentFT0M, aod::McCollisionLabels>> const& collisions,
    Particles const& particles, FiTracks const& tracks, aod::FT0s const&, aod::FDDs const&)
  {
    processGenGeneral<MChi, ExColsCentFT0M>(mcCollision, collisions, particles, tracks);
  }

  PROCESS_SWITCH(MultiplicityCounter, processGenFT0MhiNoAmb, "Process generator-level info (FT0M centrality, HI) w/o ambiguous", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MultiplicityCounter>(cfgc)};
}
