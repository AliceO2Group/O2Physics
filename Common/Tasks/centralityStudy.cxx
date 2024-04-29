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

  // Configurable Axes
  ConfigurableAxis axisMultFT0C{"axisMultFT0C", {2000, 0, 100000}, "FT0C amplitude"};
  ConfigurableAxis axisMultPVContributors{"axisMultPVContributors", {200, 0, 6000}, "Number of PV Contributors"};

  // For one-dimensional plots, where binning is no issue
  ConfigurableAxis axisMultUltraFineFT0C{"axisMultUltraFineFT0C", {60000, 0, 60000}, "FT0C amplitude"};
  ConfigurableAxis axisMultUltraFinePVContributors{"axisMultUltraFinePVContributors", {10000, 0, 10000}, "Number of PV Contributors"};

  ConfigurableAxis axisMultITSOnly{"axisMultITSOnly", {200, 0, 6000}, "Number of ITS only tracks"};
  ConfigurableAxis axisMultITSTPC{"axisMultITSTPC", {200, 0, 6000}, "Number of ITSTPC matched tracks"};

  void init(InitContext&)
  {
    if (doprocessCollisions) {
      const AxisSpec axisCollisions{100, -0.5f, 99.5f, "Number of collisions"};
      histos.add("hCollisionSelection", "hCollisionSelection", kTH1D, {{10, -0.5f, +9.5f}});
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

      histos.add("hFT0C_Collisions", "hFT0C_Collisions", kTH1D, {axisMultUltraFineFT0C});
      histos.add("hNPVContributors", "hNPVContributors", kTH1D, {axisMultUltraFinePVContributors});
    }

    if (doprocessBCs) {
      histos.add("hBCSelection", "hBCSelection", kTH1D, {{10, -0.5, 9.5f}});
      histos.add("hFT0C_BCs", "hFT0C_BCs", kTH1D, {axisMultUltraFineFT0C});
    }

    if (do2DPlots) {
      histos.add("hFT0CvsNContribs", "hFT0CvsNContribs", kTH2F, {axisMultPVContributors, axisMultFT0C});
      histos.add("hMatchedVsITSOnly", "hMatchedVsITSOnly", kTH2F, {axisMultITSOnly, axisMultITSTPC});
    }
  }

  void processCollisions(soa::Join<aod::Mults, aod::MultsExtra, aod::MultSelections>::iterator const& collision)
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

    // if we got here, we also finally fill the FT0C histogram, please
    histos.fill(HIST("hNPVContributors"), collision.multPVTotalContributors());
    histos.fill(HIST("hFT0C_Collisions"), collision.multFT0C());

    if (do2DPlots) {
      histos.fill(HIST("hFT0CvsNContribs"), collision.multNTracksPV(), collision.multFT0C());
      histos.fill(HIST("hMatchedVsITSOnly"), collision.multNTracksITSOnly(), collision.multNTracksITSTPC());
    }
  }

  void processBCs(aod::MultsBC::iterator const& multbc)
  {
    // process BCs, calculate FT0C distribution
    // conditionals suggested by FIT team (Jacek O. et al)
    histos.fill(HIST("hBCSelection"), 0); // all BCs
    if (!multbc.multBCColliding())
      return;
    histos.fill(HIST("hBCSelection"), 1); // colliding
    if (!multbc.multBCTVX())
      return;
    histos.fill(HIST("hBCSelection"), 2); // TVX
    if (!multbc.multBCFV0OrA())
      return;
    histos.fill(HIST("hBCSelection"), 3); // FV0OrA

    // if we got here, we also finally fill the FT0C histogram, please
    histos.fill(HIST("hFT0C_BCs"), multbc.multBCFT0C());
  }

  PROCESS_SWITCH(centralityStudy, processCollisions, "per-collision analysis", true);
  PROCESS_SWITCH(centralityStudy, processBCs, "per-BC analysis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<centralityStudy>(cfgc)};
}
