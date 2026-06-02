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

/// \file deltaanalysis.cxx
/// \brief  Delta(1232) resonance analysis via proton-pion invariant mass reconstruction with advance PID and background rejection cuts
/// \author Durgesh Bhatt <durgesh.bhatt@cern.ch>

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <THnSparse.h>
#include <TMath.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace
{
constexpr float massProton = o2::constants::physics::MassProton;
constexpr float massPion = o2::constants::physics::MassPionCharged;
} // namespace
namespace delta_analysis
{
static constexpr int PdgProton{2212};
static constexpr int PdgPion{211};
static constexpr int PdgDeltaPlusPlus{2224};
static constexpr int PdgDeltaZero{2114};

enum CentEstimator : int {
  kFT0M = 0,
  kFT0A = 1,
  kFT0C = 2,
  kFV0A = 3,
  kNTPV = 4
};
} // namespace delta_analysis

struct DeltaAnalysis {

  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // FIX name/configurable: single space between type and name; member name == JSON key; lowerCamelCase
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted |z-vertex| range [cm]"};
  Configurable<bool> applyOccupancyInTimeRangeCut{"applyOccupancyInTimeRangeCut", false, "Apply occupancy-in-time-range cut"};
  Configurable<int> cfgOccupancyMin{"cfgOccupancyMin", 0, "Minimum track occupancy in time range"};
  Configurable<int> cfgOccupancyMax{"cfgOccupancyMax", 9999, "Maximum track occupancy in time range"};

  Configurable<int> cfgCentralityEstimator{"cfgCentralityEstimator", 0, "Centrality estimator: 0=FT0M 1=FT0A 2=FT0C 3=FV0A 4=NTPV"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0.f, "Minimum centrality percentile"};
  Configurable<float> cfgCentMax{"cfgCentMax", 100.f, "Maximum centrality percentile"};

  Configurable<float> cfgCutPt{"cfgCutPt", 0.2f, "Minimum pT of daughter track [GeV/c]"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Maximum |eta| of daughter track"};
  // CHANGED: replaced single symmetric cfgCutY with independent cfgMinY / cfgMaxY
  Configurable<float> cfgMinY{"cfgMinY", -0.5f, "Minimum rapidity of reconstructed Delta"};
  Configurable<float> cfgMaxY{"cfgMaxY", 0.5f, "Maximum rapidity of reconstructed Delta"};

  Configurable<int> cfgMinITSClusters{"cfgMinITSClusters", 5, "Minimum ITS clusters"};
  Configurable<int> cfgMinTPCClusters{"cfgMinTPCClusters", 70, "Minimum TPC clusters found"};
  Configurable<int> numberOfInvMassBins{"numberOfInvMassBins", 120, "Number of bins along the invariant mass axis"};
  Configurable<int> cfgMinTPCCrossedRows{"cfgMinTPCCrossedRows", 70, "Minimum TPC crossed pad rows"};
  Configurable<float> cfgMinCrossedRowsOverFindable{"cfgMinCrossedRowsOverFindable", 0.8f, "Minimum ratio crossedRows/findableClusters"};
  Configurable<int> cfgMaxTPCSharedClusters{"cfgMaxTPCSharedClusters", 0, "Maximum TPC shared clusters"};
  Configurable<float> cfgMaxTPCChi2NCl{"cfgMaxTPCChi2NCl", 4.0f, "Maximum TPC chi2/NCl"};
  Configurable<float> cfgMaxITSChi2NCl{"cfgMaxITSChi2NCl", 36.0f, "Maximum ITS chi2/NCl"};
  Configurable<bool> requirePrimaryTrack{"requirePrimaryTrack", true, "Require isPrimaryTrack flag"};
  Configurable<bool> requireGlobalTrackNoDCA{"requireGlobalTrackNoDCA", true, "Require isGlobalTrackWoDCA flag"};
  Configurable<bool> requirePVContributor{"requirePVContributor", true, "Require PV-contributor flag"};

  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 0.1f, "Maximum |DCAz| for all tracks [cm]"};
  Configurable<std::vector<float>> protonDCAPtBinEdges{"protonDCAPtBinEdges", {0.0f, 0.5f, 1.0f, 2.0f, 1000.f}, "Proton pT bin edges for DCAxy cut [GeV/c]"};
  Configurable<std::vector<float>> protonMaxDCAxyPerPtBin{"protonMaxDCAxyPerPtBin", {0.10f, 0.08f, 0.05f, 0.05f}, "Max |DCAxy| for proton per pT bin [cm]"};
  Configurable<std::vector<float>> pionDCAPtBinEdges{"pionDCAPtBinEdges", {0.0f, 0.5f, 1.0f, 2.0f, 1000.f}, "Pion pT bin edges for DCAxy cut [GeV/c]"};
  Configurable<std::vector<float>> pionMaxDCAxyPerPtBin{"pionMaxDCAxyPerPtBin", {0.20f, 0.15f, 0.10f, 0.08f}, "Max |DCAxy| for pion per pT bin [cm]"};

  Configurable<bool> useTPCOnlyPID{"useTPCOnlyPID", false, "Use TPC-only PID (ignore TOF even if present)"};
  Configurable<bool> requireTOFForProton{"requireTOFForProton", false, "Require TOF signal for proton candidates"};
  Configurable<bool> requireTOFForPion{"requireTOFForPion", false, "Require TOF signal for pion candidates"};

  Configurable<float> minProtonMomentum{"minProtonMomentum", 0.f, "Minimum proton momentum for TOF PID [GeV/c]"};
  Configurable<float> minTPCNSigmaProton{"minTPCNSigmaProton", -6.0f, "Minimum (lower bound) TPC nSigma for proton"};
  Configurable<float> minTOFNSigmaProton{"minTOFNSigmaProton", -6.0f, "Minimum (lower bound) TOF nSigma for proton"};
  Configurable<float> minCombinedNSigmaProton{"minCombinedNSigmaProton", -6.0f, "Minimum combined nSigma for proton (asymmetric mode)"};
  Configurable<double> maxTPCNSigmaProton{"maxTPCNSigmaProton", 3.0, "Maximum |TPC nSigma| for proton"};
  Configurable<double> combinedNSigmaCutProton{"combinedNSigmaCutProton", 3.0, "Circular TPC+TOF combined nSigma cut for proton. Negative = asymmetric mode."};
  Configurable<std::vector<float>> protonTPCPIDMomentumBins{"protonTPCPIDMomentumBins", {0.f, 0.5f, 0.7f, 0.8f}, "Proton TPC PID momentum bin edges [GeV/c]"};
  Configurable<std::vector<float>> protonTPCNSigmaCutPerBin{"protonTPCNSigmaCutPerBin", {5.f, 3.5f, 2.5f}, "Maximum TPC nSigma for proton per momentum bin"};
  Configurable<std::vector<float>> protonTOFPIDMomentumBins{"protonTOFPIDMomentumBins", {0.f, 999.f}, "Proton TOF PID momentum bin edges [GeV/c]"};
  Configurable<std::vector<float>> protonTOFNSigmaCutPerBin{"protonTOFNSigmaCutPerBin", {3.0f}, "Maximum TOF nSigma for proton per momentum bin"};

  Configurable<float> minPionMomentum{"minPionMomentum", 0.f, "Minimum pion momentum for TOF PID [GeV/c]"};
  Configurable<float> minTPCNSigmaPion{"minTPCNSigmaPion", -6.0f, "Minimum (lower bound) TPC nSigma for pion"};
  Configurable<float> minTOFNSigmaPion{"minTOFNSigmaPion", -6.0f, "Minimum (lower bound) TOF nSigma for pion"};
  Configurable<float> minCombinedNSigmaPion{"minCombinedNSigmaPion", -6.0f, "Minimum combined nSigma for pion (asymmetric mode)"};
  Configurable<double> maxTPCNSigmaPion{"maxTPCNSigmaPion", 3.0, "Maximum |TPC nSigma| for pion"};
  Configurable<double> combinedNSigmaCutPion{"combinedNSigmaCutPion", 3.0, "Circular TPC+TOF combined nSigma cut for pion. Negative = asymmetric mode."};
  Configurable<std::vector<float>> pionTPCPIDMomentumBins{"pionTPCPIDMomentumBins", {0.f, 0.25f, 0.4f, 0.5f}, "Pion TPC PID momentum bin edges [GeV/c]"};
  Configurable<std::vector<float>> pionTPCNSigmaCutPerBin{"pionTPCNSigmaCutPerBin", {5.f, 3.5f, 2.5f}, "Maximum TPC nSigma for pion per momentum bin"};
  Configurable<std::vector<float>> pionTOFPIDMomentumBins{"pionTOFPIDMomentumBins", {0.f, 999.f}, "Pion TOF PID momentum bin edges [GeV/c]"};
  Configurable<std::vector<float>> pionTOFNSigmaCutPerBin{"pionTOFNSigmaCutPerBin", {3.0f}, "Maximum TOF nSigma for pion per momentum bin"};

  Configurable<float> tpcNSigmaVetoThreshold{"tpcNSigmaVetoThreshold", 3.0f, "Reject track if TPC nSigma of a competing species is below this value"};
  Configurable<float> tofNSigmaVetoThreshold{"tofNSigmaVetoThreshold", 3.0f, "Reject track if TOF nSigma of a competing species is below this value"};

  Configurable<bool> applyDeepAngleCut{"applyDeepAngleCut", false, "Apply minimum opening-angle cut (removes split-track background)"};
  Configurable<double> deepAngleCutValue{"deepAngleCutValue", 0.04, "Minimum opening angle between proton and pion [rad]"};

  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per signal event"};

  Configurable<bool> enableRotationalBackground{"enableRotationalBackground", false, "Compute rotational background by rotating pion phi near PI"};
  Configurable<int> numberOfRotations{"numberOfRotations", 10, "Number of pion-phi rotations for background"};
  Configurable<float> rotationAngleWindow{"rotationAngleWindow", 6.f, "Pion rotated by angles within PI +/- PI/rotationAngleWindow"};

  ConfigurableAxis cfgPtAxis{"cfgPtAxis", {VARIABLE_WIDTH, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8, 3.2, 3.6, 4.0, 5.0, 7.0, 10.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis cfgCentAxis{"cfgCentAxis", {VARIABLE_WIDTH, 0.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f, 80.f, 90.f, 100.f}, "Centrality (%)"};
  ConfigurableAxis cfgVtxAxis{"cfgVtxAxis", {VARIABLE_WIDTH, -12.f, -10.f, -9.f, -8.f, -7.f, -6.f, -5.f, -4.f, -3.f, -2.f, -1.f, 0.f, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 7.f, 8.f, 9.f, 10.f, 12.f}, "Vertex z [cm]"};
  ConfigurableAxis cfgRapAxis{"cfgRapAxis", {20, -1.0, 1.0}, "Rapidity y"};

  std::vector<float> mProtonTPCMomBins{};
  std::vector<float> mProtonTPCNSigCuts{};
  std::vector<float> mProtonTOFMomBins{};
  std::vector<float> mProtonTOFNSigCuts{};
  std::vector<float> mPionTPCMomBins{};
  std::vector<float> mPionTPCNSigCuts{};
  std::vector<float> mPionTOFMomBins{};
  std::vector<float> mPionTOFNSigCuts{};
  std::vector<float> mProtonDCAPtEdges{};
  std::vector<float> mProtonMaxDCAxy{};
  std::vector<float> mPionDCAPtEdges{};
  std::vector<float> mPionMaxDCAxy{};

  void init(InitContext const&)
  {
    auto chkSz = [](auto const& cuts, auto const& bins, const char* label) {
      const auto& cutsVec = cuts.value;
      const auto& binsVec = bins.value;
      if (cutsVec.size() != binsVec.size() - 1) {
        LOG(fatal) << "PID/DCA bin-size mismatch for: " << label
                   << " (cuts=" << cutsVec.size()
                   << ", bins=" << binsVec.size() << ")";
      }
    };
    chkSz(protonTPCNSigmaCutPerBin, protonTPCPIDMomentumBins, "proton TPC PID");
    chkSz(protonTOFNSigmaCutPerBin, protonTOFPIDMomentumBins, "proton TOF PID");
    chkSz(pionTPCNSigmaCutPerBin, pionTPCPIDMomentumBins, "pion TPC PID");
    chkSz(pionTOFNSigmaCutPerBin, pionTOFPIDMomentumBins, "pion TOF PID");
    chkSz(protonMaxDCAxyPerPtBin, protonDCAPtBinEdges, "proton DCA");
    chkSz(pionMaxDCAxyPerPtBin, pionDCAPtBinEdges, "pion DCA");

    mProtonTPCMomBins = protonTPCPIDMomentumBins.value;
    mProtonTPCNSigCuts = protonTPCNSigmaCutPerBin.value;
    mProtonTOFMomBins = protonTOFPIDMomentumBins.value;
    mProtonTOFNSigCuts = protonTOFNSigmaCutPerBin.value;
    mPionTPCMomBins = pionTPCPIDMomentumBins.value;
    mPionTPCNSigCuts = pionTPCNSigmaCutPerBin.value;
    mPionTOFMomBins = pionTOFPIDMomentumBins.value;
    mPionTOFNSigCuts = pionTOFNSigmaCutPerBin.value;
    mProtonDCAPtEdges = protonDCAPtBinEdges.value;
    mProtonMaxDCAxy = protonMaxDCAxyPerPtBin.value;
    mPionDCAPtEdges = pionDCAPtBinEdges.value;
    mPionMaxDCAxy = pionMaxDCAxyPerPtBin.value;

    const AxisSpec ptAxis{200, 0., 10., "p_{T} (GeV/c)"};
    const AxisSpec massAxis{numberOfInvMassBins, 0.8, 1.8, "M_{inv} (GeV/#it{c}^{2})"};
    const AxisSpec centAxis{cfgCentAxis, "Centrality (%)"};
    const AxisSpec vtxAxis{cfgVtxAxis, "Vertex z [cm]"};
    const AxisSpec rapAxis{cfgRapAxis, "Rapidity y"};
    const AxisSpec nSigmaTPCaxis{100, -10., 10., "n#sigma^{TPC}"};
    const AxisSpec nSigmaTOFaxis{100, -10., 10., "n#sigma^{TOF}"};
    const AxisSpec momentumAxis{200, 0., 10., "p (GeV/#it{c})"};
    const AxisSpec ptForPIDAxis{200, 0., 10., "p_{T} (GeV/#it{c})"};
    const AxisSpec dcaXYaxis{200, -0.2f, 0.2f, "DCA_{xy} (cm)"};
    const AxisSpec dcaZaxis{200, -1.2f, 1.2f, "DCA_{z} (cm)"};
    const AxisSpec tpcClusAxis{200, 0, 200, "N_{clusters}^{TPC}"};
    const AxisSpec tpcRowsAxis{200, 0, 200, "Crossed rows^{TPC}"};
    const AxisSpec occupancyAxis{200, 0, 10000, "Track occupancy [-40, 100] ns"};
    const AxisSpec etaAxis{100, -1.0, 1.0, "#eta"};
    const AxisSpec openAngleAxis{100, 0., o2::constants::math::PI, "Opening angle (rad)"};

    histos.add("Event/hVtxZ", "Vertex z; z (cm)", kTH1F, {{400, -20., 20.}});
    histos.add("Event/hNcontributor", "PV contributors; N", kTH1F, {{2001, -0.5f, 2000.5f}});
    histos.add("Event/hCentrality", "Centrality", kTH1F, {centAxis});
    histos.add("Event/hOccupancy", "Occupancy in time range", kTH1F, {occupancyAxis});

    histos.add("CentQA/hCentralityVsVtxZ", "Centrality vs vertex z", kTH2F, {vtxAxis, centAxis});
    histos.add("CentQA/hCentralityVsOccupancy", "Centrality vs occupancy", kTH2F, {occupancyAxis, centAxis});
    histos.add("CentQA/hEventCountVsCentrality", "Event count vs centrality", kTH1F, {centAxis});

    histos.add("OccupancyQA/hOccupancyVsCentralityBefore", "Occupancy vs centrality (before cut)", kTH2F, {centAxis, occupancyAxis});
    histos.add("OccupancyQA/hOccupancyVsVtxZBefore", "Occupancy vs vertex z (before cut)", kTH2F, {vtxAxis, occupancyAxis});
    histos.add("OccupancyQA/hOccupancyVsCentralityAfter", "Occupancy vs centrality (after cut)", kTH2F, {centAxis, occupancyAxis});
    histos.add("OccupancyQA/hOccupancyVsVtxZAfter", "Occupancy vs vertex z (after cut)", kTH2F, {vtxAxis, occupancyAxis});

    histos.add("QAbefore/Proton/tpcNSigmaVsMomentum", "TPC n#sigma proton vs p (before cuts)", kTH2F, {momentumAxis, nSigmaTPCaxis});
    histos.add("QAbefore/Proton/tofNSigmaVsMomentum", "TOF n#sigma proton vs p (before cuts)", kTH2F, {momentumAxis, nSigmaTOFaxis});
    histos.add("QAbefore/Proton/tofNSigmaVsTPCNSigma", "TOF vs TPC n#sigma proton (before cuts)", kTH2F, {nSigmaTPCaxis, nSigmaTOFaxis});
    histos.add("QAbefore/Pion/tpcNSigmaVsMomentum", "TPC n#sigma pion vs p (before cuts)", kTH2F, {momentumAxis, nSigmaTPCaxis});
    histos.add("QAbefore/Pion/tofNSigmaVsMomentum", "TOF n#sigma pion vs p (before cuts)", kTH2F, {momentumAxis, nSigmaTOFaxis});
    histos.add("QAbefore/Pion/tofNSigmaVsTPCNSigma", "TOF vs TPC n#sigma pion (before cuts)", kTH2F, {nSigmaTPCaxis, nSigmaTOFaxis});

    histos.add("QAafter/Proton/dcaXYvsPt", "Proton DCA_{xy} vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, dcaXYaxis});
    histos.add("QAafter/Proton/dcaZvsPt", "Proton DCA_{z} vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, dcaZaxis});
    histos.add("QAafter/Proton/tpcNSigmaVsMomentum", "Proton TPC n#sigma vs p (after cuts)", kTH2F, {momentumAxis, nSigmaTPCaxis});
    histos.add("QAafter/Proton/tpcNSigmaVsPt", "Proton TPC n#sigma vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
    histos.add("QAafter/Proton/tpcNSigmaVsCentrality", "Proton TPC n#sigma vs centrality", kTH2F, {centAxis, nSigmaTPCaxis});
    histos.add("QAafter/Proton/tofNSigmaVsMomentum", "Proton TOF n#sigma vs p (after cuts)", kTH2F, {momentumAxis, nSigmaTOFaxis});
    histos.add("QAafter/Proton/tofNSigmaVsPt", "Proton TOF n#sigma vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, nSigmaTOFaxis});
    histos.add("QAafter/Proton/tofNSigmaVsCentrality", "Proton TOF n#sigma vs centrality", kTH2F, {centAxis, nSigmaTOFaxis});
    histos.add("QAafter/Proton/tofNSigmaVsTPCNSigma", "Proton TOF vs TPC n#sigma (after cuts)", kTH2F, {nSigmaTPCaxis, nSigmaTOFaxis});
    histos.add("QAafter/Proton/tpcNSigmaPionContamVsPt", "Proton track: TPC n#sigma pion contamination", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
    histos.add("QAafter/Proton/tpcNSigmaKaonContamVsPt", "Proton track: TPC n#sigma kaon contamination", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
    histos.add("QAafter/Proton/tofNSigmaPionContamVsMomentum", "Proton track: TOF n#sigma pion contamination", kTH2F, {momentumAxis, nSigmaTOFaxis});
    histos.add("QAafter/Proton/tpcCrossedRowsVsPt", "Proton TPC crossed rows vs p_{T}", kTH2F, {ptForPIDAxis, tpcRowsAxis});
    histos.add("QAafter/Proton/tpcClustersFoundVsPt", "Proton TPC clusters found vs p_{T}", kTH2F, {ptForPIDAxis, tpcClusAxis});
    histos.add("QAafter/Proton/dcaXYdist", "Proton DCA_{xy} distribution (fine bins)", kTH1F, {dcaXYaxis});
    histos.add("QAafter/Proton/dcaZdist", "Proton DCA_{z} distribution (fine bins)", kTH1F, {dcaZaxis});

    histos.add("QAafter/Pion/dcaXYvsPt", "Pion DCA_{xy} vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, dcaXYaxis});
    histos.add("QAafter/Pion/dcaZvsPt", "Pion DCA_{z} vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, dcaZaxis});
    histos.add("QAafter/Pion/tpcNSigmaVsMomentum", "Pion TPC n#sigma vs p (after cuts)", kTH2F, {momentumAxis, nSigmaTPCaxis});
    histos.add("QAafter/Pion/tpcNSigmaVsPt", "Pion TPC n#sigma vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
    histos.add("QAafter/Pion/tpcNSigmaVsCentrality", "Pion TPC n#sigma vs centrality", kTH2F, {centAxis, nSigmaTPCaxis});
    histos.add("QAafter/Pion/tofNSigmaVsMomentum", "Pion TOF n#sigma vs p (after cuts)", kTH2F, {momentumAxis, nSigmaTOFaxis});
    histos.add("QAafter/Pion/tofNSigmaVsPt", "Pion TOF n#sigma vs p_{T} (after cuts)", kTH2F, {ptForPIDAxis, nSigmaTOFaxis});
    histos.add("QAafter/Pion/tofNSigmaVsCentrality", "Pion TOF n#sigma vs centrality", kTH2F, {centAxis, nSigmaTOFaxis});
    histos.add("QAafter/Pion/tofNSigmaVsTPCNSigma", "Pion TOF vs TPC n#sigma (after cuts)", kTH2F, {nSigmaTPCaxis, nSigmaTOFaxis});
    histos.add("QAafter/Pion/tpcNSigmaProtonContamVsPt", "Pion track: TPC n#sigma proton contamination", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
    histos.add("QAafter/Pion/tpcNSigmaKaonContamVsPt", "Pion track: TPC n#sigma kaon contamination", kTH2F, {ptForPIDAxis, nSigmaTPCaxis});
    histos.add("QAafter/Pion/tofNSigmaProtonContamVsMomentum", "Pion track: TOF n#sigma proton contamination", kTH2F, {momentumAxis, nSigmaTOFaxis});
    histos.add("QAafter/Pion/tpcCrossedRowsVsPt", "Pion TPC crossed rows vs p_{T}", kTH2F, {ptForPIDAxis, tpcRowsAxis});
    histos.add("QAafter/Pion/tpcClustersFoundVsPt", "Pion TPC clusters found vs p_{T}", kTH2F, {ptForPIDAxis, tpcClusAxis});
    histos.add("QAafter/Pion/dcaXYdist", "Pion DCA_{xy} distribution (fine bins)", kTH1F, {dcaXYaxis});
    histos.add("QAafter/Pion/dcaZdist", "Pion DCA_{z} distribution (fine bins)", kTH1F, {dcaZaxis});

    histos.add("QAChecks/Pair/hOpeningAngleBefore", "Opening angle (before cuts)", kTH1F, {openAngleAxis});
    histos.add("QAChecks/Pair/hOpeningAngleAfter", "Opening angle (after cuts)", kTH1F, {openAngleAxis});

    histos.add("Analysis/hDeltaPlusPlusInvMass", "#Delta^{++} (#rightarrow p + #pi^{+}) invariant mass", kTH2F, {ptAxis, massAxis});
    histos.add("Analysis/hAntiDeltaPlusPlusInvMass", "#bar{#Delta}^{++} (#rightarrow #bar{p} + #pi^{-}) invariant mass", kTH2F, {ptAxis, massAxis});
    histos.add("Analysis/hDeltaZeroInvMass", "#Delta^{0} (#rightarrow p + #pi^{-}) invariant mass", kTH2F, {ptAxis, massAxis});
    histos.add("Analysis/hAntiDeltaZeroInvMass", "#bar{#Delta}^{0} (#rightarrow #bar{p} + #pi^{+}) invariant mass", kTH2F, {ptAxis, massAxis});

    if (enableRotationalBackground) {
      histos.add("Analysis/hDeltaPlusPlusInvMassRot", "#Delta^{++} invariant mass - rotational background", kTH2F, {ptAxis, massAxis});
      histos.add("Analysis/hAntiDeltaPlusPlusInvMassRot", "#bar{#Delta}^{++} invariant mass - rotational background", kTH2F, {ptAxis, massAxis});
      histos.add("Analysis/hDeltaZeroInvMassRot", "#Delta^{0} invariant mass - rotational background", kTH2F, {ptAxis, massAxis});
      histos.add("Analysis/hAntiDeltaZeroInvMassRot", "#bar{#Delta}^{0} invariant mass - rotational background", kTH2F, {ptAxis, massAxis});
    }

    if (doprocessMixedEvent) {
      histos.add("Analysis/hDeltaPlusPlusInvMassEM", "#Delta^{++} invariant mass - event mixing", kTH2F, {ptAxis, massAxis});
      histos.add("Analysis/hAntiDeltaPlusPlusInvMassEM", "#bar{#Delta}^{++} invariant mass - event mixing", kTH2F, {ptAxis, massAxis});
      histos.add("Analysis/hDeltaZeroInvMassEM", "#Delta^{0} invariant mass - event mixing", kTH2F, {ptAxis, massAxis});
      histos.add("Analysis/hAntiDeltaZeroInvMassEM", "#bar{#Delta}^{0} invariant mass - event mixing", kTH2F, {ptAxis, massAxis});
    }

    histos.add("THnSparse/hDeltaPlusPlus", "THnSparse #Delta^{++} same-event", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
    histos.add("THnSparse/hAntiDeltaPlusPlus", "THnSparse #bar{#Delta}^{++} same-event", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
    histos.add("THnSparse/hDeltaZero", "THnSparse #Delta^{0} same-event", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
    histos.add("THnSparse/hAntiDeltaZero", "THnSparse #bar{#Delta}^{0} same-event", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});

    if (enableRotationalBackground) {
      histos.add("THnSparse/hDeltaPlusPlusRot", "THnSparse #Delta^{++} rotational background", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("THnSparse/hAntiDeltaPlusPlusRot", "THnSparse #bar{#Delta}^{++} rotational background", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("THnSparse/hDeltaZeroRot", "THnSparse #Delta^{0} rotational background", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("THnSparse/hAntiDeltaZeroRot", "THnSparse #bar{#Delta}^{0} rotational background", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
    }

    if (doprocessMixedEvent) {
      histos.add("THnSparse/hDeltaPlusPlusEM", "THnSparse #Delta^{++} event mixing", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("THnSparse/hAntiDeltaPlusPlusEM", "THnSparse #bar{#Delta}^{++} event mixing", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("THnSparse/hDeltaZeroEM", "THnSparse #Delta^{0} event mixing", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("THnSparse/hAntiDeltaZeroEM", "THnSparse #bar{#Delta}^{0} event mixing", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
    }

    if (doprocessMC) {
      histos.add("THnSparse/hDeltaPlusPlusMC", "THnSparse #Delta^{++} MC reconstructed", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("THnSparse/hAntiDeltaPlusPlusMC", "THnSparse #bar{#Delta}^{++} MC reconstructed", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("THnSparse/hDeltaZeroMC", "THnSparse #Delta^{0} MC reconstructed", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("THnSparse/hAntiDeltaZeroMC", "THnSparse #bar{#Delta}^{0} MC reconstructed", kTHnSparseF, {massAxis, ptAxis, centAxis, rapAxis});
      histos.add("THnSparse/hDeltaPlusPlusGen", "THnSparse #Delta^{++} generated", kTHnSparseF, {ptAxis, centAxis, rapAxis});
      histos.add("THnSparse/hAntiDeltaPlusPlusGen", "THnSparse #bar{#Delta}^{++} generated", kTHnSparseF, {ptAxis, centAxis, rapAxis});
      histos.add("THnSparse/hDeltaZeroGen", "THnSparse #Delta^{0} generated", kTHnSparseF, {ptAxis, centAxis, rapAxis});
      histos.add("THnSparse/hAntiDeltaZeroGen", "THnSparse #bar{#Delta}^{0} generated", kTHnSparseF, {ptAxis, centAxis, rapAxis});

      histos.add("Analysis/hDeltaPlusPlusInvMassMC", "#Delta^{++} invariant mass (MC truth-matched)", kTH2F, {ptAxis, massAxis});
      histos.add("Analysis/hAntiDeltaPlusPlusInvMassMC", "#bar{#Delta}^{++} invariant mass (MC truth-matched)", kTH2F, {ptAxis, massAxis});
      histos.add("Analysis/hDeltaZeroInvMassMC", "#Delta^{0} invariant mass (MC truth-matched)", kTH2F, {ptAxis, massAxis});
      histos.add("Analysis/hAntiDeltaZeroInvMassMC", "#bar{#Delta}^{0} invariant mass (MC truth-matched)", kTH2F, {ptAxis, massAxis});
      histos.add("Analysis/hDeltaPlusPlusInvMassGen", "#Delta^{++} invariant mass - generated", kTH2F, {ptAxis, massAxis});
      histos.add("Analysis/hAntiDeltaPlusPlusInvMassGen", "#bar{#Delta}^{++} invariant mass - generated", kTH2F, {ptAxis, massAxis});
      histos.add("Analysis/hDeltaZeroInvMassGen", "#Delta^{0} invariant mass - generated", kTH2F, {ptAxis, massAxis});
      histos.add("Analysis/hAntiDeltaZeroInvMassGen", "#bar{#Delta}^{0} invariant mass - generated", kTH2F, {ptAxis, massAxis});

      histos.add("QAChecks/hRecProtonDeltaPlusPlus", "Rec proton from #Delta^{++}", kTH1F, {ptAxis});
      histos.add("QAChecks/hRecProtonAntiDeltaPlusPlus", "Rec proton from #bar{#Delta}^{++}", kTH1F, {ptAxis});
      histos.add("QAChecks/hRecProtonDeltaZero", "Rec proton from #Delta^{0}", kTH1F, {ptAxis});
      histos.add("QAChecks/hRecProtonAntiDeltaZero", "Rec proton from #bar{#Delta}^{0}", kTH1F, {ptAxis});
      histos.add("QAChecks/hRecPionDeltaPlusPlus", "Rec pion from #Delta^{++}", kTH1F, {ptAxis});
      histos.add("QAChecks/hRecPionAntiDeltaPlusPlus", "Rec pion from #bar{#Delta}^{++}", kTH1F, {ptAxis});
      histos.add("QAChecks/hRecPionDeltaZero", "Rec pion from #Delta^{0}", kTH1F, {ptAxis});
      histos.add("QAChecks/hRecPionAntiDeltaZero", "Rec pion from #bar{#Delta}^{0}", kTH1F, {ptAxis});
      histos.add("QAChecks/hGenProtonDeltaPlusPlus", "Gen proton from #Delta^{++}", kTH1F, {ptAxis});
      histos.add("QAChecks/hGenProtonAntiDeltaPlusPlus", "Gen proton from #bar{#Delta}^{++}", kTH1F, {ptAxis});
      histos.add("QAChecks/hGenProtonDeltaZero", "Gen proton from #Delta^{0}", kTH1F, {ptAxis});
      histos.add("QAChecks/hGenProtonAntiDeltaZero", "Gen proton from #bar{#Delta}^{0}", kTH1F, {ptAxis});
      histos.add("QAChecks/hGenPionDeltaPlusPlus", "Gen pion from #Delta^{++}", kTH1F, {ptAxis});
      histos.add("QAChecks/hGenPionAntiDeltaPlusPlus", "Gen pion from #bar{#Delta}^{++}", kTH1F, {ptAxis});
      histos.add("QAChecks/hGenPionDeltaZero", "Gen pion from #Delta^{0}", kTH1F, {ptAxis});
      histos.add("QAChecks/hGenPionAntiDeltaZero", "Gen pion from #bar{#Delta}^{0}", kTH1F, {ptAxis});
    }
  } // end init()

  template <typename CollisionType>
  float getCentrality(CollisionType const& collision)
  {
    switch (cfgCentralityEstimator) {
      case delta_analysis::kFT0M:
        return collision.centFT0M();
      case delta_analysis::kFT0A:
        return collision.centFT0A();
      case delta_analysis::kFT0C:
        return collision.centFT0C();
      case delta_analysis::kFV0A:
        return collision.centFV0A();
      case delta_analysis::kNTPV:
        return collision.centNTPV();
      default:
        return collision.centFT0M();
    }
  }

  template <typename CollisionType>
  bool passesEventSelection(CollisionType const& collision)
  {
    if (!collision.sel8())
      return false;
    if (std::abs(collision.posZ()) > cfgCutVertex)
      return false;
    if (applyOccupancyInTimeRangeCut) {
      const int occ = collision.trackOccupancyInTimeRange();
      if (occ < cfgOccupancyMin || occ > cfgOccupancyMax)
        return false;
    }
    const float cent = getCentrality(collision);
    if (cent < cfgCentMin || cent > cfgCentMax)
      return false;
    return true;
  }

  template <typename TrackType>
  bool passesBasicTrackSelection(TrackType const& track)
  {
    if (track.itsNCls() < cfgMinITSClusters)
      return false;
    if (track.tpcNClsShared() > cfgMaxTPCSharedClusters)
      return false;
    if (track.tpcNClsFound() < cfgMinTPCClusters)
      return false;
    if (track.tpcNClsCrossedRows() < cfgMinTPCCrossedRows)
      return false;
    if (track.tpcNClsCrossedRows() < cfgMinCrossedRowsOverFindable * track.tpcNClsFindable())
      return false;
    if (track.tpcChi2NCl() > cfgMaxTPCChi2NCl)
      return false;
    if (track.itsChi2NCl() > cfgMaxITSChi2NCl)
      return false;
    if (requirePrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (requireGlobalTrackNoDCA && !track.isGlobalTrackWoDCA())
      return false;
    if (requirePVContributor && !track.isPVContributor())
      return false;
    return true;
  }

  template <typename TrackType>
  bool passesProtonDCASelection(TrackType const& track)
  {
    const int nBins = static_cast<int>(mProtonDCAPtEdges.size()) - 1;
    const float pt = track.pt();
    bool passed = false;
    for (int i = 0; i < nBins; ++i) {
      if (pt >= mProtonDCAPtEdges[i] && pt < mProtonDCAPtEdges[i + 1] &&
          std::abs(track.dcaXY()) < mProtonMaxDCAxy[i]) {
        passed = true;
        break;
      }
    }
    return passed && (std::abs(track.dcaZ()) < cfgCutDCAz);
  }

  template <typename TrackType>
  bool passesPionDCASelection(TrackType const& track)
  {
    const int nBins = static_cast<int>(mPionDCAPtEdges.size()) - 1;
    const float pt = track.pt();
    bool passed = false;
    for (int i = 0; i < nBins; ++i) {
      if (pt >= mPionDCAPtEdges[i] && pt < mPionDCAPtEdges[i + 1] &&
          std::abs(track.dcaXY()) < mPionMaxDCAxy[i]) {
        passed = true;
        break;
      }
    }
    return passed && (std::abs(track.dcaZ()) < cfgCutDCAz);
  }

  template <typename TrackType>
  bool passesProtonPID(TrackType const& track, float totalMomentum)
  {
    bool tpcPassed{false}, tofPassed{false};
    const int nTPCBins = static_cast<int>(mProtonTPCMomBins.size());
    const int nTOFBins = static_cast<int>(mProtonTOFMomBins.size());
    const float tpcNSigPi = std::abs(track.tpcNSigmaPi());
    const float tpcNSigPr = std::abs(track.tpcNSigmaPr());
    const float tofNSigPi = std::abs(track.tofNSigmaPi());
    const float tofNSigPr = std::abs(track.tofNSigmaPr());
    const float combinedNSigPr = tpcNSigPr * tpcNSigPr + tofNSigPr * tofNSigPr;
    const float combinedNSigPi = tpcNSigPi * tpcNSigPi + tofNSigPi * tofNSigPi;
    const float circularCutSq = static_cast<float>(combinedNSigmaCutProton * combinedNSigmaCutProton);

    const float circularVetoCutSq = tpcNSigmaVetoThreshold * tpcNSigmaVetoThreshold + tofNSigmaVetoThreshold * tofNSigmaVetoThreshold;

    if (!useTPCOnlyPID && track.hasTOF()) {
      if (combinedNSigmaCutProton < 0 && totalMomentum >= minProtonMomentum) {
        if (track.tofNSigmaPr() < minTOFNSigmaProton)
          return false;
        for (int i = 0; i < nTOFBins - 1; ++i) {
          if (totalMomentum >= mProtonTOFMomBins[i] && totalMomentum < mProtonTOFMomBins[i + 1] &&
              tofNSigPr < mProtonTOFNSigCuts[i] && tofNSigPi > tofNSigmaVetoThreshold) {
            tofPassed = true;
            break;
          }
        }
        if (track.tpcNSigmaPr() < minCombinedNSigmaProton)
          return false;
        if (tpcNSigPr < static_cast<float>(maxTPCNSigmaProton) && tpcNSigPi > tpcNSigmaVetoThreshold)
          tpcPassed = true;
      } else if (combinedNSigmaCutProton > 0 && totalMomentum >= minProtonMomentum) {
        if (combinedNSigPr < circularCutSq && combinedNSigPi > circularVetoCutSq) {
          tpcPassed = true;
          tofPassed = true;
        }
      }
      if (totalMomentum < minProtonMomentum && tpcNSigPr < static_cast<float>(maxTPCNSigmaProton)) {
        tpcPassed = true;
        tofPassed = true;
      }
    } else {
      tofPassed = true;
      if (track.tpcNSigmaPr() < minTPCNSigmaProton)
        return false;
      for (int i = 0; i < nTPCBins - 1; ++i) {
        if (totalMomentum >= mProtonTPCMomBins[i] && totalMomentum < mProtonTPCMomBins[i + 1] &&
            tpcNSigPr < mProtonTPCNSigCuts[i] && tpcNSigPi > tpcNSigmaVetoThreshold) {
          tpcPassed = true;
          break;
        }
      }
    }
    return tpcPassed && tofPassed;
  }

  template <typename TrackType>
  bool passesPionPID(TrackType const& track, float totalMomentum)
  {
    bool tpcPassed{false}, tofPassed{false};
    const int nTPCBins = static_cast<int>(mPionTPCMomBins.size());
    const int nTOFBins = static_cast<int>(mPionTOFMomBins.size());
    const float tpcNSigPi = std::abs(track.tpcNSigmaPi());
    const float tpcNSigPr = std::abs(track.tpcNSigmaPr());
    const float tofNSigPi = std::abs(track.tofNSigmaPi());
    const float tofNSigPr = std::abs(track.tofNSigmaPr());
    const float combinedNSigPi = tpcNSigPi * tpcNSigPi + tofNSigPi * tofNSigPi;
    const float combinedNSigPr = tpcNSigPr * tpcNSigPr + tofNSigPr * tofNSigPr;
    const float circularCutSq = static_cast<float>(combinedNSigmaCutPion * combinedNSigmaCutPion);
    const float circularVetoCutSq = tpcNSigmaVetoThreshold * tpcNSigmaVetoThreshold + tofNSigmaVetoThreshold * tofNSigmaVetoThreshold;

    if (!useTPCOnlyPID && track.hasTOF()) {
      if (combinedNSigmaCutPion < 0 && totalMomentum >= minPionMomentum) {
        if (track.tofNSigmaPi() < minTOFNSigmaPion)
          return false;
        for (int i = 0; i < nTOFBins - 1; ++i) {
          if (totalMomentum >= mPionTOFMomBins[i] && totalMomentum < mPionTOFMomBins[i + 1] &&
              tofNSigPi < mPionTOFNSigCuts[i] && tofNSigPr > tofNSigmaVetoThreshold) {
            tofPassed = true;
            break;
          }
        }
        if (track.tpcNSigmaPi() < minCombinedNSigmaPion)
          return false;
        if (tpcNSigPi < static_cast<float>(maxTPCNSigmaPion) && tpcNSigPr > tpcNSigmaVetoThreshold)
          tpcPassed = true;
      } else if (combinedNSigmaCutPion > 0 && totalMomentum >= minPionMomentum) {
        if (combinedNSigPi < circularCutSq && combinedNSigPr > circularVetoCutSq) {
          tpcPassed = true;
          tofPassed = true;
        }
      }
      if (totalMomentum < minPionMomentum && tpcNSigPi < static_cast<float>(maxTPCNSigmaPion)) {
        tpcPassed = true;
        tofPassed = true;
      }
    } else {
      tofPassed = true;
      if (track.tpcNSigmaPi() < minTPCNSigmaPion)
        return false;
      for (int i = 0; i < nTPCBins - 1; ++i) {
        if (totalMomentum >= mPionTPCMomBins[i] && totalMomentum < mPionTPCMomBins[i + 1] &&
            tpcNSigPi < mPionTPCNSigCuts[i] && tpcNSigPr > tpcNSigmaVetoThreshold) {
          tpcPassed = true;
          break;
        }
      }
    }
    return tpcPassed && tofPassed;
  }

  template <typename TrackType>
  void fillQAProton(TrackType const& track, float totalMomentum, float centralityPercent)
  {
    const float pt = track.pt();
    const float tpcNSigPr = track.tpcNSigmaPr();
    histos.fill(HIST("QAafter/Proton/dcaXYvsPt"), pt, track.dcaXY());
    histos.fill(HIST("QAafter/Proton/dcaZvsPt"), pt, track.dcaZ());
    histos.fill(HIST("QAafter/Proton/dcaXYdist"), track.dcaXY());
    histos.fill(HIST("QAafter/Proton/dcaZdist"), track.dcaZ());
    histos.fill(HIST("QAafter/Proton/tpcNSigmaVsMomentum"), totalMomentum, tpcNSigPr);
    histos.fill(HIST("QAafter/Proton/tpcNSigmaVsPt"), pt, tpcNSigPr);
    histos.fill(HIST("QAafter/Proton/tpcNSigmaVsCentrality"), centralityPercent, tpcNSigPr);
    histos.fill(HIST("QAafter/Proton/tpcNSigmaPionContamVsPt"), pt, track.tpcNSigmaPi());
    histos.fill(HIST("QAafter/Proton/tpcNSigmaKaonContamVsPt"), pt, track.tpcNSigmaKa());
    histos.fill(HIST("QAafter/Proton/tpcCrossedRowsVsPt"), pt, track.tpcNClsCrossedRows());
    histos.fill(HIST("QAafter/Proton/tpcClustersFoundVsPt"), pt, track.tpcNClsFound());
    if (!useTPCOnlyPID && track.hasTOF()) {
      const float tofNSigPr = track.tofNSigmaPr();
      histos.fill(HIST("QAafter/Proton/tofNSigmaVsMomentum"), totalMomentum, tofNSigPr);
      histos.fill(HIST("QAafter/Proton/tofNSigmaVsPt"), pt, tofNSigPr);
      histos.fill(HIST("QAafter/Proton/tofNSigmaVsCentrality"), centralityPercent, tofNSigPr);
      histos.fill(HIST("QAafter/Proton/tofNSigmaVsTPCNSigma"), tpcNSigPr, tofNSigPr);
      histos.fill(HIST("QAafter/Proton/tofNSigmaPionContamVsMomentum"), totalMomentum, track.tofNSigmaPi());
    }
  }

  template <typename TrackType>
  void fillQAPion(TrackType const& track, float totalMomentum, float centralityPercent)
  {
    const float pt = track.pt();
    const float tpcNSigPi = track.tpcNSigmaPi();
    histos.fill(HIST("QAafter/Pion/dcaXYvsPt"), pt, track.dcaXY());
    histos.fill(HIST("QAafter/Pion/dcaZvsPt"), pt, track.dcaZ());
    histos.fill(HIST("QAafter/Pion/dcaXYdist"), track.dcaXY());
    histos.fill(HIST("QAafter/Pion/dcaZdist"), track.dcaZ());
    histos.fill(HIST("QAafter/Pion/tpcNSigmaVsMomentum"), totalMomentum, tpcNSigPi);
    histos.fill(HIST("QAafter/Pion/tpcNSigmaVsPt"), pt, tpcNSigPi);
    histos.fill(HIST("QAafter/Pion/tpcNSigmaVsCentrality"), centralityPercent, tpcNSigPi);
    histos.fill(HIST("QAafter/Pion/tpcNSigmaProtonContamVsPt"), pt, track.tpcNSigmaPr());
    histos.fill(HIST("QAafter/Pion/tpcNSigmaKaonContamVsPt"), pt, track.tpcNSigmaKa());
    histos.fill(HIST("QAafter/Pion/tpcCrossedRowsVsPt"), pt, track.tpcNClsCrossedRows());
    histos.fill(HIST("QAafter/Pion/tpcClustersFoundVsPt"), pt, track.tpcNClsFound());
    if (!useTPCOnlyPID && track.hasTOF()) {
      const float tofNSigPi = track.tofNSigmaPi();
      histos.fill(HIST("QAafter/Pion/tofNSigmaVsMomentum"), totalMomentum, tofNSigPi);
      histos.fill(HIST("QAafter/Pion/tofNSigmaVsPt"), pt, tofNSigPi);
      histos.fill(HIST("QAafter/Pion/tofNSigmaVsCentrality"), centralityPercent, tofNSigPi);
      histos.fill(HIST("QAafter/Pion/tofNSigmaVsTPCNSigma"), tpcNSigPi, tofNSigPi);
      histos.fill(HIST("QAafter/Pion/tofNSigmaProtonContamVsMomentum"), totalMomentum, track.tofNSigmaPr());
    }
  }

  void fillDeltaHistogramSameEvent(int protonSign, int pionSign, float pairPt, float pairMass, float centrality, float rapidity)
  {
    if (protonSign > 0) {
      if (pionSign > 0) {
        histos.fill(HIST("Analysis/hDeltaPlusPlusInvMass"), pairPt, pairMass);
        histos.fill(HIST("THnSparse/hDeltaPlusPlus"), pairMass, pairPt, centrality, rapidity);
      } else {
        histos.fill(HIST("Analysis/hDeltaZeroInvMass"), pairPt, pairMass);
        histos.fill(HIST("THnSparse/hDeltaZero"), pairMass, pairPt, centrality, rapidity);
      }
    } else {
      if (pionSign < 0) {
        histos.fill(HIST("Analysis/hAntiDeltaPlusPlusInvMass"), pairPt, pairMass);
        histos.fill(HIST("THnSparse/hAntiDeltaPlusPlus"), pairMass, pairPt, centrality, rapidity);
      } else {
        histos.fill(HIST("Analysis/hAntiDeltaZeroInvMass"), pairPt, pairMass);
        histos.fill(HIST("THnSparse/hAntiDeltaZero"), pairMass, pairPt, centrality, rapidity);
      }
    }
  }

  void fillDeltaHistogramMixedEvent(int protonSign, int pionSign, float pairPt, float pairMass, float centrality, float rapidity)
  {
    if (protonSign > 0) {
      if (pionSign > 0) {
        histos.fill(HIST("Analysis/hDeltaPlusPlusInvMassEM"), pairPt, pairMass);
        histos.fill(HIST("THnSparse/hDeltaPlusPlusEM"), pairMass, pairPt, centrality, rapidity);
      } else {
        histos.fill(HIST("Analysis/hDeltaZeroInvMassEM"), pairPt, pairMass);
        histos.fill(HIST("THnSparse/hDeltaZeroEM"), pairMass, pairPt, centrality, rapidity);
      }
    } else {
      if (pionSign < 0) {
        histos.fill(HIST("Analysis/hAntiDeltaPlusPlusInvMassEM"), pairPt, pairMass);
        histos.fill(HIST("THnSparse/hAntiDeltaPlusPlusEM"), pairMass, pairPt, centrality, rapidity);
      } else {
        histos.fill(HIST("Analysis/hAntiDeltaZeroInvMassEM"), pairPt, pairMass);
        histos.fill(HIST("THnSparse/hAntiDeltaZeroEM"), pairMass, pairPt, centrality, rapidity);
      }
    }
  }

  void fillDeltaHistogramMC(int protonSign, int pionSign, float pairPt, float pairMass, float centrality, float rapidity)
  {
    if (protonSign > 0) {
      if (pionSign > 0) {
        histos.fill(HIST("Analysis/hDeltaPlusPlusInvMassMC"), pairPt, pairMass);
        histos.fill(HIST("THnSparse/hDeltaPlusPlusMC"), pairMass, pairPt, centrality, rapidity);
      } else {
        histos.fill(HIST("Analysis/hDeltaZeroInvMassMC"), pairPt, pairMass);
        histos.fill(HIST("THnSparse/hDeltaZeroMC"), pairMass, pairPt, centrality, rapidity);
      }
    } else {
      if (pionSign < 0) {
        histos.fill(HIST("Analysis/hAntiDeltaPlusPlusInvMassMC"), pairPt, pairMass);
        histos.fill(HIST("THnSparse/hAntiDeltaPlusPlusMC"), pairMass, pairPt, centrality, rapidity);
      } else {
        histos.fill(HIST("Analysis/hAntiDeltaZeroInvMassMC"), pairPt, pairMass);
        histos.fill(HIST("THnSparse/hAntiDeltaZeroMC"), pairMass, pairPt, centrality, rapidity);
      }
    }
  }

  void fillRotationalBackground(int protonSign, int pionSign, float pxProton, float pyProton, float pzProton, float pionPhi, float pionPt, float pzPion, float centrality)
  {
    if (numberOfRotations <= 0)
      return;
    const float weight = 1.f / static_cast<float>(numberOfRotations);
    const float rotWindowHalf = o2::constants::math::PI / static_cast<float>(rotationAngleWindow);

    for (int iRot = 0; iRot < numberOfRotations; ++iRot) {
      float rotAngle;
      if (numberOfRotations == 1) {
        rotAngle = o2::constants::math::PI;
      } else {
        rotAngle = (o2::constants::math::PI - rotWindowHalf) +
                   static_cast<float>(iRot) * (2.f * rotWindowHalf / static_cast<float>(numberOfRotations - 1));
      }
      const float newPhi = RecoDecay::constrainAngle(pionPhi + rotAngle, 0.f);
      const float pxPionRot = pionPt * std::cos(newPhi);
      const float pyPionRot = pionPt * std::sin(newPhi);
      const std::array<std::array<float, 3>, 2> rotMomenta = {
        std::array<float, 3>{pxProton, pyProton, pzProton},
        std::array<float, 3>{pxPionRot, pyPionRot, pzPion}};
      const float rotMass = RecoDecay::m(rotMomenta, std::array{massProton, massPion});
      const float rotPt = RecoDecay::pt(std::array{pxProton + pxPionRot, pyProton + pyPionRot});
      const float rotY = RecoDecay::y(
        std::array{pxProton + pxPionRot, pyProton + pyPionRot, pzProton + pzPion}, rotMass);
      if (rotY < cfgMinY || rotY > cfgMaxY)
        continue;

      if (protonSign > 0) {
        if (pionSign > 0) {
          histos.fill(HIST("Analysis/hDeltaPlusPlusInvMassRot"), rotPt, rotMass, weight);
          histos.fill(HIST("THnSparse/hDeltaPlusPlusRot"), rotMass, rotPt, centrality, rotY, weight);
        } else {
          histos.fill(HIST("Analysis/hDeltaZeroInvMassRot"), rotPt, rotMass, weight);
          histos.fill(HIST("THnSparse/hDeltaZeroRot"), rotMass, rotPt, centrality, rotY, weight);
        }
      } else {
        if (pionSign < 0) {
          histos.fill(HIST("Analysis/hAntiDeltaPlusPlusInvMassRot"), rotPt, rotMass, weight);
          histos.fill(HIST("THnSparse/hAntiDeltaPlusPlusRot"), rotMass, rotPt, centrality, rotY, weight);
        } else {
          histos.fill(HIST("Analysis/hAntiDeltaZeroInvMassRot"), rotPt, rotMass, weight);
          histos.fill(HIST("THnSparse/hAntiDeltaZeroRot"), rotMass, rotPt, centrality, rotY, weight);
        }
      }
    }
  }

  void fillPairQABefore(float pxPr, float pyPr, float pzPr, float pxPi, float pyPi, float pzPi, float momPr, float momPi)
  {
    const float cosAngle = std::clamp((pxPr * pxPi + pyPr * pyPi + pzPr * pzPi) / (momPr * momPi), -1.f, 1.f);
    histos.fill(HIST("QAChecks/Pair/hOpeningAngleBefore"), std::acos(cosAngle));
  }

  void fillPairQAAfter(float pxPr, float pyPr, float pzPr, float pxPi, float pyPi, float pzPi, float momPr, float momPi)
  {
    const float cosAngle = std::clamp((pxPr * pxPi + pyPr * pyPi + pzPr * pzPi) / (momPr * momPi), -1.f, 1.f);
    histos.fill(HIST("QAChecks/Pair/hOpeningAngleAfter"), std::acos(cosAngle));
  }

  struct TrackCandidate {
    float px, py, pz, pt, eta, phi, dcaXY, dcaZ;
    float mom;
    int sign;
    bool hasTOF;
    float tpcNSigmaPr, tpcNSigmaPi;
    uint64_t globalIndex;
  };

  template <typename TrackCollection>
  void buildCandidatePools(TrackCollection const& tracks, std::vector<TrackCandidate>& protons, std::vector<TrackCandidate>& pions)
  {
    protons.clear();
    pions.clear();
    for (auto const& track : tracks) {
      if (!passesBasicTrackSelection(track))
        continue;
      const float mom = RecoDecay::p(track.px(), track.py(), track.pz());
      const bool isProton = passesProtonPID(track, mom) && passesProtonDCASelection(track) &&
                            !(requireTOFForProton && !track.hasTOF());
      const bool isPion = passesPionPID(track, mom) && passesPionDCASelection(track) &&
                          !(requireTOFForPion && !track.hasTOF());
      if (!isProton && !isPion)
        continue;
      TrackCandidate cand;
      cand.px = track.px();
      cand.py = track.py();
      cand.pz = track.pz();
      cand.pt = track.pt();
      cand.eta = track.eta();
      cand.phi = track.phi();
      cand.dcaXY = track.dcaXY();
      cand.dcaZ = track.dcaZ();
      cand.mom = mom;
      cand.sign = track.sign();
      cand.hasTOF = track.hasTOF();
      cand.tpcNSigmaPr = track.tpcNSigmaPr();
      cand.tpcNSigmaPi = track.tpcNSigmaPi();
      cand.globalIndex = track.globalIndex();
      if (isProton)
        protons.push_back(cand);
      if (isPion)
        pions.push_back(cand);
    }
  }

  template <bool isMixed>
  void fillInvariantMassHistogramsFromPools(
    std::vector<TrackCandidate> const& protonPool,
    std::vector<TrackCandidate> const& pionPool,
    float centrality,
    bool fillPairQA = false)
  {
    for (auto const& protonCand : protonPool) {
      for (auto const& pionCand : pionPool) {

        if constexpr (!isMixed) {
          if (protonCand.globalIndex == pionCand.globalIndex)
            continue;
        }

        const float pxPr = protonCand.px, pyPr = protonCand.py, pzPr = protonCand.pz;
        const float pxPi = pionCand.px, pyPi = pionCand.py, pzPi = pionCand.pz;

        if constexpr (!isMixed) {
          if (fillPairQA)
            fillPairQABefore(pxPr, pyPr, pzPr, pxPi, pyPi, pzPi, protonCand.mom, pionCand.mom);
        }

        if (applyDeepAngleCut) {
          const float cosAngle = std::clamp(
            (pxPr * pxPi + pyPr * pyPi + pzPr * pzPi) / (protonCand.mom * pionCand.mom), -1.f, 1.f);
          if (std::acos(cosAngle) < static_cast<float>(deepAngleCutValue))
            continue;
        }

        const std::array<std::array<float, 3>, 2> bothMomenta = {
          std::array<float, 3>{pxPr, pyPr, pzPr},
          std::array<float, 3>{pxPi, pyPi, pzPi}};
        const float pairMass = RecoDecay::m(bothMomenta, std::array{massProton, massPion});
        const float pairPt = RecoDecay::pt(std::array{pxPr + pxPi, pyPr + pyPi});
        const float pairY = RecoDecay::y(std::array{pxPr + pxPi, pyPr + pyPi, pzPr + pzPi}, pairMass);

        if (pairY < cfgMinY || pairY > cfgMaxY)
          continue;

        if constexpr (!isMixed) {
          if (fillPairQA)
            fillPairQAAfter(pxPr, pyPr, pzPr, pxPi, pyPi, pzPi, protonCand.mom, pionCand.mom);
        }

        if constexpr (isMixed) {
          fillDeltaHistogramMixedEvent(protonCand.sign, pionCand.sign, pairPt, pairMass, centrality, pairY);
        } else {
          fillDeltaHistogramSameEvent(protonCand.sign, pionCand.sign, pairPt, pairMass, centrality, pairY);
          if (enableRotationalBackground)
            fillRotationalBackground(protonCand.sign, pionCand.sign,
                                     pxPr, pyPr, pzPr, pionCand.phi, pionCand.pt, pzPi, centrality);
        }
      }
    }
  }

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPt);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::CentNTPVs>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr>>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::McTrackLabels>>;

  Preslice<TrackCandidates> perCol = aod::track::collisionId;
  Preslice<TrackCandidatesMC> perColMC = aod::track::collisionId;

  // ── Event-mixing binning policies (one per centrality estimator) ────────────
  using BinningTypeFT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  using BinningTypeFT0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0A>;
  using BinningTypeFT0C = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  using BinningTypeFV0A = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFV0A>;
  using BinningTypeNTPV = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentNTPV>;

  BinningTypeFT0M binningFT0M{{cfgVtxAxis, cfgCentAxis}, true};
  BinningTypeFT0A binningFT0A{{cfgVtxAxis, cfgCentAxis}, true};
  BinningTypeFT0C binningFT0C{{cfgVtxAxis, cfgCentAxis}, true};
  BinningTypeFV0A binningFV0A{{cfgVtxAxis, cfgCentAxis}, true};
  BinningTypeNTPV binningNTPV{{cfgVtxAxis, cfgCentAxis}, true};

  SameKindPair<EventCandidates, TrackCandidates, BinningTypeFT0M> pairFT0M{binningFT0M, cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidates, TrackCandidates, BinningTypeFT0A> pairFT0A{binningFT0A, cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidates, TrackCandidates, BinningTypeFT0C> pairFT0C{binningFT0C, cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidates, TrackCandidates, BinningTypeFV0A> pairFV0A{binningFV0A, cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidates, TrackCandidates, BinningTypeNTPV> pairNTPV{binningNTPV, cfgNoMixedEvents, -1, &cache};

  std::vector<TrackCandidate> mProtonPool{};
  std::vector<TrackCandidate> mPionPool{};

  void processSameEvent(EventCandidates const& collisions,
                        TrackCandidates const& tracks,
                        aod::BCs const&)
  {
    for (auto const& collision : collisions) {
      if (!passesEventSelection(collision))
        continue;
      const float centrality = getCentrality(collision);
      const int occupancy = collision.trackOccupancyInTimeRange();

      histos.fill(HIST("OccupancyQA/hOccupancyVsCentralityBefore"), centrality, occupancy);
      histos.fill(HIST("OccupancyQA/hOccupancyVsVtxZBefore"), collision.posZ(), occupancy);
      histos.fill(HIST("Event/hNcontributor"), collision.numContrib());
      histos.fill(HIST("Event/hVtxZ"), collision.posZ());
      histos.fill(HIST("Event/hCentrality"), centrality);
      histos.fill(HIST("Event/hOccupancy"), occupancy);
      histos.fill(HIST("CentQA/hCentralityVsVtxZ"), collision.posZ(), centrality);
      histos.fill(HIST("CentQA/hCentralityVsOccupancy"), occupancy, centrality);
      histos.fill(HIST("CentQA/hEventCountVsCentrality"), centrality);
      histos.fill(HIST("OccupancyQA/hOccupancyVsCentralityAfter"), centrality, occupancy);
      histos.fill(HIST("OccupancyQA/hOccupancyVsVtxZAfter"), collision.posZ(), occupancy);

      const uint64_t collIdx = collision.globalIndex();
      auto perColTracks = tracks.sliceBy(perCol, collIdx);
      perColTracks.bindExternalIndices(&tracks);

      for (auto const& track : perColTracks) {
        const float mom = RecoDecay::p(track.px(), track.py(), track.pz());
        histos.fill(HIST("QAbefore/Proton/tpcNSigmaVsMomentum"), mom, track.tpcNSigmaPr());
        if (track.hasTOF()) {
          histos.fill(HIST("QAbefore/Proton/tofNSigmaVsMomentum"), mom, track.tofNSigmaPr());
          histos.fill(HIST("QAbefore/Proton/tofNSigmaVsTPCNSigma"), track.tpcNSigmaPr(), track.tofNSigmaPr());
        }
        histos.fill(HIST("QAbefore/Pion/tpcNSigmaVsMomentum"), mom, track.tpcNSigmaPi());
        if (track.hasTOF()) {
          histos.fill(HIST("QAbefore/Pion/tofNSigmaVsMomentum"), mom, track.tofNSigmaPi());
          histos.fill(HIST("QAbefore/Pion/tofNSigmaVsTPCNSigma"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        }
        if (!passesBasicTrackSelection(track))
          continue;
        if (passesProtonPID(track, mom) && passesProtonDCASelection(track))
          fillQAProton(track, mom, centrality);
        if (passesPionPID(track, mom) && passesPionDCASelection(track))
          fillQAPion(track, mom, centrality);
      }
      buildCandidatePools(perColTracks, mProtonPool, mPionPool);
      fillInvariantMassHistogramsFromPools<false>(mProtonPool, mPionPool, centrality, true);
    }
  }
  PROCESS_SWITCH(DeltaAnalysis, processSameEvent, "Process same event", true);

  // ── Shared mixed-event logic ─────────────────────────────────────────────
  template <typename PairType>
  void runMixedEvent(PairType& mixingPair)
  {
    for (auto const& [c1, tracks1, c2, tracks2] : mixingPair) {
      if (!passesEventSelection(c1) || !passesEventSelection(c2))
        continue;
      const float centrality = getCentrality(c1);
      std::vector<TrackCandidate> protonPool1, pionPool1, protonPool2, pionPool2;
      buildCandidatePools(tracks1, protonPool1, pionPool1);
      buildCandidatePools(tracks2, protonPool2, pionPool2);
      fillInvariantMassHistogramsFromPools<true>(protonPool1, pionPool2, centrality);
      fillInvariantMassHistogramsFromPools<true>(protonPool2, pionPool1, centrality);
    }
  }

  void processMixedEvent(EventCandidates const&, TrackCandidates const&)
  {
    switch (cfgCentralityEstimator) {
      case delta_analysis::kFT0M:
        runMixedEvent(pairFT0M);
        break;
      case delta_analysis::kFT0A:
        runMixedEvent(pairFT0A);
        break;
      case delta_analysis::kFT0C:
        runMixedEvent(pairFT0C);
        break;
      case delta_analysis::kFV0A:
        runMixedEvent(pairFV0A);
        break;
      case delta_analysis::kNTPV:
        runMixedEvent(pairNTPV);
        break;
      default:
        runMixedEvent(pairFT0M);
        break;
    }
  }
  PROCESS_SWITCH(DeltaAnalysis, processMixedEvent, "Process mixed event", true);

  void processMC(soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As, aod::CentNTPVs> const& collisions, aod::BCs const&, TrackCandidatesMC const& tracks, aod::McParticles const& mcParticles)
  {
    constexpr float kGenCentrality = 1.f;

    for (auto const& collision : collisions) {
      if (!collision.sel8() || std::abs(collision.posZ()) > cfgCutVertex)
        continue;
      const float centrality = getCentrality(collision);
      histos.fill(HIST("Event/hNcontributor"), collision.numContrib());
      histos.fill(HIST("Event/hVtxZ"), collision.posZ());

      const uint64_t collIdx = collision.globalIndex();
      auto perColTracks = tracks.sliceBy(perColMC, collIdx);
      perColTracks.bindExternalIndices(&tracks);

      for (auto const& t0 : perColTracks) {
        if (!passesBasicTrackSelection(t0) || !t0.has_mcParticle())
          continue;
        const auto mcTrack = t0.mcParticle();
        const float mom = RecoDecay::p(t0.px(), t0.py(), t0.pz());
        if (std::abs(mcTrack.pdgCode()) == delta_analysis::PdgProton && passesProtonPID(t0, mom)) {
          for (const auto& mother : mcTrack.mothers_as<aod::McParticles>()) {
            if (std::abs(mother.pdgCode()) == delta_analysis::PdgDeltaPlusPlus) {
              if (t0.sign() < 0)
                histos.fill(HIST("QAChecks/hRecProtonAntiDeltaPlusPlus"), t0.pt());
              else
                histos.fill(HIST("QAChecks/hRecProtonDeltaPlusPlus"), t0.pt());
            } else if (std::abs(mother.pdgCode()) == delta_analysis::PdgDeltaZero) {
              if (t0.sign() < 0)
                histos.fill(HIST("QAChecks/hRecProtonAntiDeltaZero"), t0.pt());
              else
                histos.fill(HIST("QAChecks/hRecProtonDeltaZero"), t0.pt());
            }
          }
        }
        if (std::abs(mcTrack.pdgCode()) == delta_analysis::PdgPion && passesPionPID(t0, mom)) {
          for (const auto& mother : mcTrack.mothers_as<aod::McParticles>()) {
            if (std::abs(mother.pdgCode()) == delta_analysis::PdgDeltaPlusPlus) {
              if (t0.sign() < 0)
                histos.fill(HIST("QAChecks/hRecPionAntiDeltaPlusPlus"), t0.pt());
              else
                histos.fill(HIST("QAChecks/hRecPionDeltaPlusPlus"), t0.pt());
            } else if (std::abs(mother.pdgCode()) == delta_analysis::PdgDeltaZero) {
              if (t0.sign() < 0)
                histos.fill(HIST("QAChecks/hRecPionAntiDeltaZero"), t0.pt());
              else
                histos.fill(HIST("QAChecks/hRecPionDeltaZero"), t0.pt());
            }
          }
        }
      }

      for (auto const& [t0, t1] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(perColTracks, perColTracks))) {
        if (t0.globalIndex() == t1.globalIndex())
          continue;
        if (!passesBasicTrackSelection(t0) || !passesBasicTrackSelection(t1))
          continue;
        const float momT0 = RecoDecay::p(t0.px(), t0.py(), t0.pz());
        const float momT1 = RecoDecay::p(t1.px(), t1.py(), t1.pz());
        if (!passesProtonPID(t0, momT0) || !passesPionPID(t1, momT1))
          continue;
        if (!passesProtonDCASelection(t0) || !passesPionDCASelection(t1))
          continue;
        if (!t0.has_mcParticle() || !t1.has_mcParticle())
          continue;
        const auto mcProton = t0.mcParticle();
        const auto mcPion = t1.mcParticle();
        if (std::abs(mcProton.pdgCode()) != delta_analysis::PdgProton)
          continue;
        if (std::abs(mcPion.pdgCode()) != delta_analysis::PdgPion)
          continue;
        bool foundMother = false;
        for (const auto& motherPr : mcProton.mothers_as<aod::McParticles>()) {
          for (const auto& motherPi : mcPion.mothers_as<aod::McParticles>()) {
            if (motherPr != motherPi)
              continue;
            if (std::abs(motherPr.pdgCode()) != delta_analysis::PdgDeltaPlusPlus &&
                std::abs(motherPr.pdgCode()) != delta_analysis::PdgDeltaZero)
              continue;
            foundMother = true;
            break;
          }
          if (foundMother)
            break;
        }
        if (!foundMother)
          continue;
        const std::array<std::array<float, 3>, 2> momenta = {
          std::array<float, 3>{t0.px(), t0.py(), t0.pz()},
          std::array<float, 3>{t1.px(), t1.py(), t1.pz()}};
        const float pairMass = RecoDecay::m(momenta, std::array{massProton, massPion});
        const float pairPt = RecoDecay::pt(std::array{t0.px() + t1.px(), t0.py() + t1.py()});
        const float pairY = RecoDecay::y(std::array{t0.px() + t1.px(), t0.py() + t1.py(), t0.pz() + t1.pz()}, pairMass);
        if (pairY < cfgMinY || pairY > cfgMaxY)
          continue;
        fillDeltaHistogramMC(t0.sign(), t1.sign(), pairPt, pairMass, centrality, pairY);
      }
    }

    for (auto const& mcParticle : mcParticles) {
      const int pdg = mcParticle.pdgCode();
      if (std::abs(pdg) != delta_analysis::PdgDeltaPlusPlus &&
          std::abs(pdg) != delta_analysis::PdgDeltaZero)
        continue;
      if (mcParticle.y() < cfgMinY || mcParticle.y() > cfgMaxY)
        continue;
      const auto daughters = mcParticle.daughters_as<aod::McParticles>();
      bool hasPr = false, hasPi = false;
      float ptPr = -999.f, ptPi = -999.f;
      double eTotal = 0.;
      for (const auto& d : daughters) {
        if (std::abs(d.pdgCode()) == delta_analysis::PdgProton) {
          hasPr = true;
          ptPr = d.pt();
          eTotal += d.e();
        } else if (std::abs(d.pdgCode()) == delta_analysis::PdgPion) {
          hasPi = true;
          ptPi = d.pt();
          eTotal += d.e();
        }
      }
      if (!hasPr || !hasPi)
        continue;
      const float pSq = mcParticle.p() * mcParticle.p();
      const float genMass = std::sqrt(std::max(0., eTotal * eTotal - static_cast<double>(pSq)));
      const float genPt = mcParticle.pt();
      const float genY = mcParticle.y();
      if (pdg == delta_analysis::PdgDeltaPlusPlus) {
        histos.fill(HIST("Analysis/hDeltaPlusPlusInvMassGen"), genPt, genMass);
        histos.fill(HIST("THnSparse/hDeltaPlusPlusGen"), genPt, kGenCentrality, genY);
        histos.fill(HIST("QAChecks/hGenProtonDeltaPlusPlus"), ptPr);
        histos.fill(HIST("QAChecks/hGenPionDeltaPlusPlus"), ptPi);
      } else if (pdg == -delta_analysis::PdgDeltaPlusPlus) {
        histos.fill(HIST("Analysis/hAntiDeltaPlusPlusInvMassGen"), genPt, genMass);
        histos.fill(HIST("THnSparse/hAntiDeltaPlusPlusGen"), genPt, kGenCentrality, genY);
        histos.fill(HIST("QAChecks/hGenProtonAntiDeltaPlusPlus"), ptPr);
        histos.fill(HIST("QAChecks/hGenPionAntiDeltaPlusPlus"), ptPi);
      } else if (pdg == delta_analysis::PdgDeltaZero) {
        histos.fill(HIST("Analysis/hDeltaZeroInvMassGen"), genPt, genMass);
        histos.fill(HIST("THnSparse/hDeltaZeroGen"), genPt, kGenCentrality, genY);
        histos.fill(HIST("QAChecks/hGenProtonDeltaZero"), ptPr);
        histos.fill(HIST("QAChecks/hGenPionDeltaZero"), ptPi);
      } else if (pdg == -delta_analysis::PdgDeltaZero) {
        histos.fill(HIST("Analysis/hAntiDeltaZeroInvMassGen"), genPt, genMass);
        histos.fill(HIST("THnSparse/hAntiDeltaZeroGen"), genPt, kGenCentrality, genY);
        histos.fill(HIST("QAChecks/hGenProtonAntiDeltaZero"), ptPr);
        histos.fill(HIST("QAChecks/hGenPionAntiDeltaZero"), ptPi);
      }
    }
  }
  PROCESS_SWITCH(DeltaAnalysis, processMC, "Process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // FIX name/o2-task: TaskName is redundant when it equals the derived device name
  return WorkflowSpec{adaptAnalysisTask<DeltaAnalysis>(cfgc)};
}
