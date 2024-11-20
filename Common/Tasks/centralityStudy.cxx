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
// This task does dedicated centrality studies for understanding the
// Run 3 Pb-Pb centrality selections in 2023 data. It is compatible with
// derived data.

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace o2;
using namespace o2::framework;

using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;

struct centralityStudy {
  // Raw multiplicities
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables
  Configurable<bool> do2DPlots{"do2DPlots", true, "0 - no, 1 - yes"};
  Configurable<bool> applySel8{"applySel8", true, "0 - no, 1 - yes"};
  Configurable<bool> applyVtxZ{"applyVtxZ", true, "0 - no, 1 - yes"};

  // Apply extra event selections
  Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border"};
  Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border"};
  Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", true, "require events with at least one ITS-TPC track"};
  Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", true, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference"};
  Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", true, "require events with at least one of vertex contributors matched to TOF"};
  Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", true, "require events with at least one of vertex contributors matched to TRD"};
  Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", true, "reject collisions in case of pileup with another collision in the same foundBC"};

  Configurable<bool> rejectITSinROFpileupStandard{"rejectITSinROFpileupStandard", false, "reject collisions in case of in-ROF ITS pileup (standard)"};
  Configurable<bool> rejectITSinROFpileupStrict{"rejectITSinROFpileupStrict", false, "reject collisions in case of in-ROF ITS pileup (strict)"};
  Configurable<bool> rejectCollInTimeRangeNarrow{"rejectCollInTimeRangeNarrow", false, "reject if extra colls in time range (narrow)"};

  Configurable<bool> selectUPCcollisions{"selectUPCcollisions", false, "select collisions tagged with UPC flag"};

  Configurable<bool> selectCollidingBCs{"selectCollidingBCs", true, "BC analysis: select colliding BCs"};
  Configurable<bool> selectTVX{"selectTVX", true, "BC analysis: select TVX"};
  Configurable<bool> selectFV0OrA{"selectFV0OrA", true, "BC analysis: select FV0OrA"};
  Configurable<float> vertexZwithT0{"vertexZwithT0", 1000.0f, "require a certain vertex-Z in BC analysis"};

  Configurable<float> minTimeDelta{"minTimeDelta", -1.0f, "reject collision if another collision is this close or less in time"};
  Configurable<float> minFT0CforVertexZ{"minFT0CforVertexZ", 250, "minimum FT0C for vertex-Z profile calculation"};

  Configurable<float> scaleSignalFT0C{"scaleSignalFT0C", 1.00f, "scale FT0C signal for convenience"};
  Configurable<float> scaleSignalFT0M{"scaleSignalFT0M", 1.00f, "scale FT0M signal for convenience"};
  Configurable<float> scaleSignalFV0A{"scaleSignalFV0A", 1.00f, "scale FV0A signal for convenience"};

  // Configurable Axes
  ConfigurableAxis axisMultFT0C{"axisMultFT0C", {2000, 0, 100000}, "FT0C amplitude"};
  ConfigurableAxis axisMultPVContributors{"axisMultPVContributors", {200, 0, 6000}, "Number of PV Contributors"};

  // For one-dimensional plots, where binning is no issue
  ConfigurableAxis axisMultUltraFineFV0A{"axisMultUltraFineFV0A", {60000, 0, 60000}, "FV0A amplitude"};
  ConfigurableAxis axisMultUltraFineFT0M{"axisMultUltraFineFT0M", {50000, 0, 200000}, "FT0M amplitude"};
  ConfigurableAxis axisMultUltraFineFT0C{"axisMultUltraFineFT0C", {60000, 0, 60000}, "FT0C amplitude"};
  ConfigurableAxis axisMultUltraFinePVContributors{"axisMultUltraFinePVContributors", {10000, 0, 10000}, "Number of PV Contributors"};

  ConfigurableAxis axisMultITSOnly{"axisMultITSOnly", {200, 0, 6000}, "Number of ITS only tracks"};
  ConfigurableAxis axisMultITSTPC{"axisMultITSTPC", {200, 0, 6000}, "Number of ITSTPC matched tracks"};

  // For centrality studies if requested
  ConfigurableAxis axisCentrality{"axisCentrality", {100, 0, 100}, "FT0C percentile"};
  ConfigurableAxis axisPVChi2{"axisPVChi2", {300, 0, 30}, "FT0C percentile"};
  ConfigurableAxis axisDeltaTime{"axisDeltaTime", {300, 0, 300}, "#Delta time"};

  // For profile Z
  ConfigurableAxis axisPVz{"axisPVz", {400, -20.0f, +20.0f}, "PVz (cm)"};

  void init(InitContext&)
  {
    if (doprocessCollisions || doprocessCollisionsWithCentrality) {
      const AxisSpec axisCollisions{100, -0.5f, 99.5f, "Number of collisions"};
      histos.add("hCollisionSelection", "hCollisionSelection", kTH1D, {{20, -0.5f, +19.5f}});
      histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
      histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
      histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(3, "posZ cut");
      histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
      histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
      histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(6, "kIsVertexITSTPC");
      histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(7, "kIsGoodZvtxFT0vsPV");
      histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(8, "kIsVertexTOFmatched");
      histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(9, "kIsVertexTRDmatched");
      histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(10, "kNoSameBunchPileup");
      histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(11, "Neighbour rejection");
      histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(12, "no ITS in-ROF pileup (standard)");
      histos.get<TH1>(HIST("hCollisionSelection"))->GetXaxis()->SetBinLabel(13, "no ITS in-ROF pileup (strict)");

      histos.add("hFT0C_Collisions", "hFT0C_Collisions", kTH1D, {axisMultUltraFineFT0C});
      histos.add("hFT0M_Collisions", "hFT0M_Collisions", kTH1D, {axisMultUltraFineFT0M});
      histos.add("hFV0A_Collisions", "hFV0A_Collisions", kTH1D, {axisMultUltraFineFV0A});
      histos.add("hNPVContributors", "hNPVContributors", kTH1D, {axisMultUltraFinePVContributors});

      histos.add("hFT0CvsPVz_Collisions_All", "hFT0CvsPVz_Collisions_All", kTProfile, {axisPVz});
      histos.add("hFT0CvsPVz_Collisions", "hFT0CvsPVz_Collisions", kTProfile, {axisPVz});
    }

    if (doprocessBCs) {
      histos.add("hBCSelection", "hBCSelection", kTH1D, {{20, -0.5, 19.5f}});
      histos.add("hFT0C_BCs", "hFT0C_BCs", kTH1D, {axisMultUltraFineFT0C});
      histos.add("hFT0M_BCs", "hFT0M_BCs", kTH1D, {axisMultUltraFineFT0M});
      histos.add("hFV0A_BCs", "hFV0A_BCs", kTH1D, {axisMultUltraFineFV0A});

      histos.add("hFT0CvsPVz_BCs_All", "hFT0CvsPVz_BCs_All", kTProfile, {axisPVz});
      histos.add("hFT0CvsPVz_BCs", "hFT0CvsPVz_BCs", kTProfile, {axisPVz});

      histos.add("hVertexZ_BCvsCO", "hVertexZ_BCvsCO", kTH2D, {axisPVz, axisPVz});
    }

    if (do2DPlots) {
      histos.add("hFT0CvsNContribs", "hFT0CvsNContribs", kTH2F, {axisMultPVContributors, axisMultFT0C});
      histos.add("hMatchedVsITSOnly", "hMatchedVsITSOnly", kTH2F, {axisMultITSOnly, axisMultITSTPC});
    }

    if (doprocessCollisionsWithCentrality) {
      // in case requested: do vs centrality debugging
      histos.add("hCentrality", "hCentrality", kTH1F, {axisCentrality});
      histos.add("hNContribsVsCentrality", "hNContribsVsCentrality", kTH2F, {axisCentrality, axisMultPVContributors});
      histos.add("hNITSTPCTracksVsCentrality", "hNITSTPCTracksVsCentrality", kTH2F, {axisCentrality, axisMultPVContributors});
      histos.add("hNITSOnlyTracksVsCentrality", "hNITSOnlyTracksVsCentrality", kTH2F, {axisCentrality, axisMultPVContributors});
      histos.add("hNGlobalTracksVsCentrality", "hNGlobalTracksVsCentrality", kTH2F, {axisCentrality, axisMultPVContributors});
      histos.add("hPVChi2VsCentrality", "hPVChi2VsCentrality", kTH2F, {axisCentrality, axisPVChi2});
      histos.add("hDeltaTimeVsCentrality", "hDeltaTimeVsCentrality", kTH2F, {axisCentrality, axisDeltaTime});
    }
  }

  template <typename TCollision>
  void genericProcessCollision(TCollision collision)
  // process this collisions
  {

    histos.fill(HIST("hCollisionSelection"), 0); // all collisions
    if (applySel8 && !collision.multSel8())
      return;
    histos.fill(HIST("hCollisionSelection"), 1);
    if (applyVtxZ && TMath::Abs(collision.multPVz()) > 10)
      return;
    histos.fill(HIST("hCollisionSelection"), 2);

    // _______________________________________________________
    // Extra event selections start here
    if (rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 3 /* Not at ITS ROF border */);

    if (rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 4 /* Not at TF border */);

    if (requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 5 /* Contains at least one ITS-TPC track */);

    if (requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 6 /* PV position consistency check */);

    if (requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 7 /* PV with at least one contributor matched with TOF */);

    if (requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 8 /* PV with at least one contributor matched with TRD */);

    if (rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 9 /* Not at same bunch pile-up */);

    // do this only if information is available
    if constexpr (requires { collision.timeToNext(); }) {
      float timeToNeighbour = TMath::Min(
        std::abs(collision.timeToNext()),
        std::abs(collision.timeToPrevious()));
      histos.fill(HIST("hDeltaTimeVsCentrality"), collision.centFT0C(), timeToNeighbour);
      if (timeToNeighbour < minTimeDelta) {
        return;
      }
      histos.fill(HIST("hCollisionSelection"), 10 /* has suspicious neighbour */);
    }

    if (rejectITSinROFpileupStandard && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 11 /* Not ITS ROF pileup (standard) */);

    if (rejectITSinROFpileupStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 12 /* Not ITS ROF pileup (strict) */);

    if (selectUPCcollisions && collision.flags() < 1) { // if zero then NOT upc, otherwise UPC
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 13 /* is UPC event */);

    if (rejectCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return;
    }
    histos.fill(HIST("hCollisionSelection"), 14 /* Not ITS ROF pileup (strict) */);

    // if we got here, we also finally fill the FT0C histogram, please
    histos.fill(HIST("hNPVContributors"), collision.multPVTotalContributors());
    histos.fill(HIST("hFT0C_Collisions"), collision.multFT0C() * scaleSignalFT0C);
    histos.fill(HIST("hFT0M_Collisions"), (collision.multFT0A() + collision.multFT0C()) * scaleSignalFT0M);
    histos.fill(HIST("hFV0A_Collisions"), collision.multFV0A() * scaleSignalFV0A);
    histos.fill(HIST("hFT0CvsPVz_Collisions_All"), collision.multPVz(), collision.multFT0C() * scaleSignalFT0C);
    if (collision.multFT0C() > minFT0CforVertexZ) {
      histos.fill(HIST("hFT0CvsPVz_Collisions"), collision.multPVz(), collision.multFT0C() * scaleSignalFT0C);
    }
    if (do2DPlots) {
      histos.fill(HIST("hFT0CvsNContribs"), collision.multNTracksPV(), collision.multFT0C() * scaleSignalFT0C);
      histos.fill(HIST("hMatchedVsITSOnly"), collision.multNTracksITSOnly(), collision.multNTracksITSTPC());
    }

    // if the table has centrality information
    if constexpr (requires { collision.centFT0C(); }) {
      // process FT0C centrality plots
      histos.fill(HIST("hCentrality"), collision.centFT0C());
      histos.fill(HIST("hNContribsVsCentrality"), collision.centFT0C(), collision.multPVTotalContributors());
      histos.fill(HIST("hNITSTPCTracksVsCentrality"), collision.centFT0C(), collision.multNTracksITSTPC());
      histos.fill(HIST("hNITSOnlyTracksVsCentrality"), collision.centFT0C(), collision.multNTracksITSOnly());
      histos.fill(HIST("hNGlobalTracksVsCentrality"), collision.centFT0C(), collision.multNTracksGlobal());
      histos.fill(HIST("hPVChi2VsCentrality"), collision.centFT0C(), collision.multPVChi2());
    }
  }

  void processCollisions(soa::Join<aod::Mults, aod::MultsExtra, aod::MultSelections>::iterator const& collision)
  {
    genericProcessCollision(collision);
  }

  void processCollisionsWithCentrality(soa::Join<aod::Mults, aod::MultsExtra, aod::MultSelections, aod::CentFT0Cs, aod::MultsGlobal>::iterator const& collision)
  {
    genericProcessCollision(collision);
  }

  void processCollisionsWithCentralityWithNeighbours(soa::Join<aod::Mults, aod::MultsExtra, aod::MultSelections, aod::CentFT0Cs, aod::MultsGlobal, aod::MultNeighs>::iterator const& collision)
  {
    genericProcessCollision(collision);
  }

  void processBCs(soa::Join<aod::BC2Mults, aod::MultBCs>::iterator const& multbc, soa::Join<aod::Mults, aod::MultsExtra, aod::MultSelections, aod::CentFT0Cs, aod::MultsGlobal> const&)
  {
    // process BCs, calculate FT0C distribution
    // conditionals suggested by FIT team (Jacek O. et al)
    histos.fill(HIST("hBCSelection"), 0); // all BCs
    if (selectCollidingBCs && !multbc.multBCColliding())
      return;
    histos.fill(HIST("hBCSelection"), 1); // colliding
    if (selectTVX && !multbc.multBCTVX())
      return;
    histos.fill(HIST("hBCSelection"), 2); // TVX
    if (selectFV0OrA && !multbc.multBCFV0OrA())
      return;
    histos.fill(HIST("hBCSelection"), 3); // FV0OrA
    if (vertexZwithT0 < 100.0f) {
      if (!multbc.multBCFT0PosZValid())
        return;
      if (TMath::Abs(multbc.multBCFT0PosZ()) > vertexZwithT0)
        return;
    }
    histos.fill(HIST("hBCSelection"), 4); // FV0OrA

    // if we got here, we also finally fill the FT0C histogram, please
    histos.fill(HIST("hFT0C_BCs"), multbc.multBCFT0C() * scaleSignalFT0C);
    histos.fill(HIST("hFT0M_BCs"), (multbc.multBCFT0A() + multbc.multBCFT0C()) * scaleSignalFT0M);
    histos.fill(HIST("hFV0A_BCs"), multbc.multBCFV0A() * scaleSignalFV0A);
    if (multbc.multBCFT0PosZValid()) {
      histos.fill(HIST("hFT0CvsPVz_BCs_All"), multbc.multBCFT0PosZ(), multbc.multBCFT0C() * scaleSignalFT0C);
      if (multbc.multBCFT0C() > minFT0CforVertexZ) {
        histos.fill(HIST("hFT0CvsPVz_BCs"), multbc.multBCFT0PosZ(), multbc.multBCFT0C() * scaleSignalFT0C);
      }
    }

    if (multbc.has_ft0Mult()) {
      auto multco = multbc.ft0Mult_as<soa::Join<aod::Mults, aod::MultsExtra, aod::MultSelections, aod::CentFT0Cs, aod::MultsGlobal>>();
      if (multbc.multBCFT0PosZValid()) {
        histos.fill(HIST("hVertexZ_BCvsCO"), multco.multPVz(), multbc.multBCFT0PosZ());
      }
    }
  }

  PROCESS_SWITCH(centralityStudy, processCollisions, "per-collision analysis", false);
  PROCESS_SWITCH(centralityStudy, processCollisionsWithCentrality, "per-collision analysis", true);
  PROCESS_SWITCH(centralityStudy, processCollisionsWithCentralityWithNeighbours, "per-collision analysis", false);
  PROCESS_SWITCH(centralityStudy, processBCs, "per-BC analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<centralityStudy>(cfgc)};
}
