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
///
/// \file doubleResonanceScan.cxx
/// \brief Resonance Scanner with ResoTracks and ResoMicroTracks
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
/// \since 27/03/2025
///

#include <vector>
#include <TLorentzVector.h>

#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "CommonConstants/MathConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;
// Extract STEP
// Handle resomicrotracks
struct DoubleResonanceScan {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Configurables
  struct : ConfigurableGroup {
    ConfigurableAxis cfgBinsPt{"cfgBinsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
    ConfigurableAxis cfgBinsPtQA{"cfgBinsPtQA", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0}, "Binning of the pT axis"};
    ConfigurableAxis cfgBinsCent{"cfgBinsCent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};
    ConfigurableAxis cfgBinsVtxZ{"cfgBinsVtxZ", {VARIABLE_WIDTH, -10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0}, "Binning of the z-vertex axis"};
    Configurable<int> cNbinsDiv{"cNbinsDiv", 1, "Integer to divide the number of bins"};
    Configurable<int> cNbinsDivQA{"cNbinsDivQA", 1, "Integer to divide the number of bins for QA"};
    Configurable<int> cfgInvMassNBins1{"cfgInvMassNBins1", 400, "Number of bins for the invariant mass pair"};
    Configurable<float> cfgInvMassPairStart1{"cfgInvMassPairStart1", 0.9, "Start of the invariant mass pair"};
    Configurable<float> cfgInvMassPairEnd1{"cfgInvMassPairEnd1", 1.3, "Start of the invariant mass pair"};
    Configurable<int> cfgInvMassNBins2{"cfgInvMassNBins2", 400, "Number of bins for the invariant mass pair"};
    Configurable<float> cfgInvMassPairStart2{"cfgInvMassPairStart2", 0.9, "Start of the invariant mass pair"};
    Configurable<float> cfgInvMassPairEnd2{"cfgInvMassPairEnd2", 1.3, "Start of the invariant mass pair"};
    Configurable<int> cfgInvMassNBinsReso{"cfgInvMassNBinsReso", 800, "Number of bins for the invariant mass final pair"};
    Configurable<float> cfgInvMassPairStartReso{"cfgInvMassPairStartReso", 2.4, "Start of the invariant mass final pair"};
    Configurable<float> cfgInvMassPairEndReso{"cfgInvMassPairEndReso", 4.0, "Start of the invariant mass final pair"};
  } AxisConfig;

  struct : ConfigurableGroup {
    Configurable<bool> cfgFillQAPlots{"cfgFillQAPlots", true, "Fill QA plots"};
  } AnalysisConfig;

  // Configurable for min pT cut
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};

  // Track selection
  // primary track condition
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};

  // DCA Selections
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.15, "Track DCAr cut to PV Maximum"};
  Configurable<bool> cUsePtDependentDCArCut{"cUsePtDependentDCArCut", false, "Use Pt dependent DCAr cut"};
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 0.15, "Track DCAz cut to PV Maximum"};

  // PID selection
  Configurable<int> cfgFirstDaughter{"cfgFirstDaughter", 0, "code of the first daughter, 0: pion, 1: kaon, 2: proton"};
  Configurable<int> cfgSecondDaughter{"cfgSecondDaughter", 0, "code of the second daughter, 0: pion, 1: kaon, 2: proton"};
  Configurable<int> cfgThirdDaughter{"cfgThirdDaughter", 0, "code of the third daughter, 0: pion, 1: kaon, 2: proton"};
  Configurable<int> cfgFourthDaughter{"cfgFourthDaughter", 0, "code of the fourth daughter, 0: pion, 1: kaon, 2: proton"};
  // PID selection values
  Configurable<float> nSigmaCutTPC{"nSigmaCutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nSigmaCutTOF{"nSigmaCutTOF", 3.0, "Value of the TOF Nsigma cut, if negative, TOF is not used"};

  struct : ConfigurableGroup {
    Configurable<std::vector<float>> cfgPairMassesLow{"cfgPairMassesLow", {1.01, 1.01}, "Low mass cut for pair_1 pair_2"};
    Configurable<std::vector<float>> cfgPairMassesHigh{"cfgPairMassesHigh", {1.03, 1.03}, "High mass cut for pair_1 pair_2"};
    Configurable<std::vector<float>> cfgPairOACut{"cfgPairOACut", {0.04, 0.04}, "Opening angle cut for pair_1 pair_2"};
  } PairCuts;

  struct : ConfigurableGroup {
    Configurable<float> cfgPairOALow{"cfgPairOALow", -999, "Low opening angle cut for pair_1 pair_2"};
    Configurable<float> cfgPairOAHigh{"cfgPairOAHigh", 999, "High opening angle cut for pair_1 pair_2"};
  } ResoCuts;

  /// Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - z-vertex"};

  /// Rotation background
  struct : ConfigurableGroup {
    Configurable<bool> cfgFillRotBkg{"cfgFillRotBkg", true, "Fill rotated background"};
    Configurable<float> cfgMinRot{"cfgMinRot", 5.0 * constants::math::PI / 6.0, "Minimum of rotation"};
    Configurable<float> cfgMaxRot{"cfgMaxRot", 7.0 * constants::math::PI / 6.0, "Maximum of rotation"};
    Configurable<int> cfgNrotBkg{"cfgNrotBkg", 9, "Number of rotated copies (background) per each original candidate"};
  } BkgEstimationConfig;

  float mass1{MassPionCharged}, mass2{MassPionCharged}, mass3{MassPionCharged}, mass4{MassPionCharged};

  float getMassFromCode(int code)
  {
    switch (code) {
      case 0:
        return MassPionCharged;
      case 1:
        return MassKaonCharged;
      case 2:
        return MassProton;
      default:
        return 0.f;
    }
  }
  void init(o2::framework::InitContext&)
  {
    LOG(info) << "Initializing DoubleResonanceScan";
    LOG(info) << "First Daughter: " << cfgFirstDaughter << " Second Daughter: " << cfgSecondDaughter << " Third Daughter: " << cfgThirdDaughter << " Fourth Daughter: " << cfgFourthDaughter;
    mass1 = getMassFromCode(cfgFirstDaughter);
    mass2 = getMassFromCode(cfgSecondDaughter);
    mass3 = getMassFromCode(cfgThirdDaughter);
    mass4 = getMassFromCode(cfgFourthDaughter);
    LOG(info) << "Masses: " << mass1 << " " << mass2 << " " << mass3 << " " << mass4;

    AxisSpec centAxis = {AxisConfig.cfgBinsCent, "T0M (%)"};
    AxisSpec vtxzAxis = {AxisConfig.cfgBinsVtxZ, "Z Vertex (cm)"};
    AxisSpec ptAxis = {AxisConfig.cfgBinsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {AxisConfig.cfgBinsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec dcaAxis = {200 / AxisConfig.cNbinsDivQA, 0, 2, "DCA (cm)"};
    AxisSpec dcaxyAxis = {100 / AxisConfig.cNbinsDivQA, 0, 1, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {100 / AxisConfig.cNbinsDivQA, 0, 1, "DCA_{#it{z}} (cm)"};
    AxisSpec yAxis = {50, -1, 1, "Rapidity"};
    AxisSpec invMassAxisPair1 = {AxisConfig.cfgInvMassNBins1, AxisConfig.cfgInvMassPairStart1, AxisConfig.cfgInvMassPairEnd1, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec invMassAxisPair2 = {AxisConfig.cfgInvMassNBins2, AxisConfig.cfgInvMassPairStart2, AxisConfig.cfgInvMassPairEnd2, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec invMassAxisReso = {AxisConfig.cfgInvMassNBinsReso, AxisConfig.cfgInvMassPairStartReso, AxisConfig.cfgInvMassPairEndReso, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec pidQAAxis = {130 / AxisConfig.cNbinsDivQA, -6.5, 6.5};
    // Event QA
    histos.add("hVertexZ", "hVertexZ", HistType::kTH1F, {{300, -15., 15.}});
    histos.add("hMultiplicityPercent", "Multiplicity Percentile", kTH1F, {{120, 0.0f, 120.0f}});

    if (AnalysisConfig.cfgFillQAPlots) {
      // First Daughter
      histos.add("hEta_1", "Eta distribution", kTH1F, {yAxis});
      histos.add("hPhi_1", "Phi distribution", kTH1F, {{200 / AxisConfig.cNbinsDivQA, 0, TwoPI}});
      histos.add("hPt_1", "Pt distribution", kTH1F, {ptAxisQA});
      histos.add("hDCAr_1", "DCAr distribution", kTH1F, {dcaxyAxis});
      histos.add("hDCAz_1", "DCAz distribution", kTH1F, {dcazAxis});
      histos.add("hNsigmaTPC_1", "nSigmaTPC distribution", kTH1F, {pidQAAxis});
      histos.add("hNsigmaTOF_1", "nSigmaTOF distribution", kTH1F, {pidQAAxis});

      // Second Daughter
      histos.add("hEta_2", "Eta distribution", kTH1F, {yAxis});
      histos.add("hPhi_2", "Phi distribution", kTH1F, {{200 / AxisConfig.cNbinsDivQA, 0, TwoPI}});
      histos.add("hPt_2", "Pt distribution", kTH1F, {ptAxisQA});
      histos.add("hDCAr_2", "DCAr distribution", kTH1F, {dcaxyAxis});
      histos.add("hDCAz_2", "DCAz distribution", kTH1F, {dcazAxis});
      histos.add("hNsigmaTPC_2", "nSigmaTPC distribution", kTH1F, {pidQAAxis});
      histos.add("hNsigmaTOF_2", "nSigmaTOF distribution", kTH1F, {pidQAAxis});

      // Third Daughter
      histos.add("hEta_3", "Eta distribution", kTH1F, {yAxis});
      histos.add("hPhi_3", "Phi distribution", kTH1F, {{200 / AxisConfig.cNbinsDivQA, 0, TwoPI}});
      histos.add("hPt_3", "Pt distribution", kTH1F, {ptAxisQA});
      histos.add("hDCAr_3", "DCAr distribution", kTH1F, {dcaxyAxis});
      histos.add("hDCAz_3", "DCAz distribution", kTH1F, {dcazAxis});
      histos.add("hNsigmaTPC_3", "nSigmaTPC distribution", kTH1F, {pidQAAxis});
      histos.add("hNsigmaTOF_3", "nSigmaTOF distribution", kTH1F, {pidQAAxis});

      // Forth Daughter
      histos.add("hEta_4", "Eta distribution", kTH1F, {yAxis});
      histos.add("hPhi_4", "Phi distribution", kTH1F, {{200 / AxisConfig.cNbinsDivQA, 0, TwoPI}});
      histos.add("hPt_4", "Pt distribution", kTH1F, {ptAxisQA});
      histos.add("hDCAr_4", "DCAr distribution", kTH1F, {dcaxyAxis});
      histos.add("hDCAz_4", "DCAz distribution", kTH1F, {dcazAxis});
      histos.add("hNsigmaTPC_4", "nSigmaTPC distribution", kTH1F, {pidQAAxis});
      histos.add("hNsigmaTOF_4", "nSigmaTOF distribution", kTH1F, {pidQAAxis});

      // First Pair
      histos.add("hPairInvMass_1", "Invariant mass distribution", kTH1F, {invMassAxisPair1});
      histos.add("hPairPt_1", "Pt distribution", kTH1F, {ptAxis});
      histos.add("hPairOA_1", "Opening angle distribution", kTH1F, {AxisSpec{100, 0, PI, "Opening Angle (rad)"}});

      // Second Pair
      histos.add("hPairInvMass_2", "Invariant mass distribution", kTH1F, {invMassAxisPair2});
      histos.add("hPairPt_2", "Pt distribution", kTH1F, {ptAxis});
      histos.add("hPairOA_2", "Opening angle distribution", kTH1F, {AxisSpec{100, 0, PI, "Opening Angle (rad)"}});

      // Resonance
      histos.add("h2PairInvMass", "Invariant mass distribution", kTH2F, {invMassAxisPair1, invMassAxisPair2});
      histos.add("hResoInvMass", "Invariant mass distribution", kTH1F, {invMassAxisReso});
      histos.add("hResoOA", "Opening angle distribution", kTH1F, {AxisSpec{100, 0, PI, "Opening Angle (rad)"}});
      histos.add("THnResoInvMass", "Invariant mass distribution with other axes", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisReso});
      histos.add("THnResoInvMassBkg", "Invariant mass distribution with other axes", HistType::kTHnSparseD, {centAxis, ptAxis, invMassAxisReso});
    }

    LOG(info) << "Size of the histograms in double resonance scan with table combination:";
    histos.print();
  }

  template <bool IsResoMicrotrack, typename TrackType>
  bool trackCut(const TrackType track)
  {
    if constexpr (!IsResoMicrotrack) {
      if (std::abs(track.pt()) < cMinPtcut)
        return false;
      if (cUsePtDependentDCArCut) {
        // Tuned on the LHC22f anchored MC LHC23d1d on primary pions. 7 Sigmas of the resolution
        if (std::abs(track.dcaXY()) > (0.004 + (0.013 / track.pt())))
          return false;
      } else {
        if (std::abs(track.dcaXY()) > cMaxDCArToPVcut)
          return false;
      }
      if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut)
        return false;
      if (cfgPrimaryTrack && !track.isPrimaryTrack())
        return false;
      if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
        return false;
      if (cfgPVContributor && !track.isPVContributor())
        return false;
    } else {
      if (std::abs(track.pt()) < cMinPtcut)
        return false;
      if (cUsePtDependentDCArCut) {
        if (o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(track.trackSelectionFlags()) > -Epsilon)
          return false;
      } else {
        if (o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(track.trackSelectionFlags()) > cMaxDCArToPVcut - Epsilon)
          return false;
      }
      if (o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(track.trackSelectionFlags()) > cMaxDCAzToPVcut - Epsilon)
        return false;
      if (cfgPrimaryTrack && !track.isPrimaryTrack())
        return false;
      if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
        return false;
      if (cfgPVContributor && !track.isPVContributor())
        return false;
    }
    return true;
  }

  template <bool IsResoMicrotrack, typename T>
  float getPIDTPC(const T& candidate, int particleType)
  {
    if constexpr (!IsResoMicrotrack) {
      // switch based on particleType
      switch (particleType) {
        case 0: // pion
          return candidate.tpcNSigmaPi();
        case 1: // kaon
          return candidate.tpcNSigmaKa();
        case 2: // proton
          return candidate.tpcNSigmaPr();
        default:
          return -999;
      }
    } else {
      switch (particleType) {
        case 0: // pion
          return o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(candidate.pidNSigmaPiFlag());
        case 1: // kaon
          return o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(candidate.pidNSigmaKaFlag());
        case 2: // proton
          return o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(candidate.pidNSigmaPrFlag());
        default:
          return -999;
      }
    }
  }

  template <bool IsResoMicrotrack, typename T>
  float getPIDTOF(const T& candidate, int particleType)
  {
    if constexpr (!IsResoMicrotrack) {
      // switch based on particleType
      switch (particleType) {
        case 0: // pion
          return candidate.hasTOF() ? candidate.tofNSigmaPi() : -999;
        case 1: // kaon
          return candidate.hasTOF() ? candidate.tofNSigmaKa() : -999;
        case 2: // proton
          return candidate.hasTOF() ? candidate.tofNSigmaPr() : -999;
        default:
          return -999;
      }
    } else {
      switch (particleType) {
        case 0: // pion
          return candidate.hasTOF() ? o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(candidate.pidNSigmaPiFlag()) : -999;
        case 1: // kaon
          return candidate.hasTOF() ? o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(candidate.pidNSigmaKaFlag()) : -999;
        case 2: // proton
          return candidate.hasTOF() ? o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(candidate.pidNSigmaPrFlag()) : -999;
        default:
          return -999;
      }
    }
  }

  template <bool IsResoMicrotrack, typename T>
  bool selectionPID(const T& candidate, int particleType)
  {
    bool tpcPass = std::abs(getPIDTPC<IsResoMicrotrack>(candidate, particleType)) < nSigmaCutTPC;
    bool tofPass = (nSigmaCutTOF > 0) ? std::abs(getPIDTOF<IsResoMicrotrack>(candidate, particleType)) < nSigmaCutTOF : true;
    if (tpcPass && tofPass) {
      return true;
    }
    return false;
  }

  template <bool IsResoMicrotrack, typename TracksType>
  std::vector<int> selectTrackIndicesWithQA(const TracksType& dTracks, int daughterType, int daughterIndex)
  {
    std::vector<int> selectedIndices;
    for (auto const& track : dTracks) {
      if (!trackCut<IsResoMicrotrack>(track)) {
        continue;
      }
      if (!selectionPID<IsResoMicrotrack>(track, daughterType)) {
        continue;
      }

      if (AnalysisConfig.cfgFillQAPlots) {
        auto dcaXY = -999;
        auto dcaZ = -999;
        if constexpr (!IsResoMicrotrack) {
          dcaXY = track.dcaXY();
          dcaZ = track.dcaZ();
        } else {
          dcaXY = o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(track.trackSelectionFlags());
          dcaZ = o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(track.trackSelectionFlags());
        }
        // FIXME: Apply better method
        switch (daughterIndex) {
          case 0: // First
            histos.fill(HIST("hEta_1"), track.eta());
            histos.fill(HIST("hPhi_1"), track.phi());
            histos.fill(HIST("hPt_1"), track.pt());
            histos.fill(HIST("hDCAr_1"), dcaXY);
            histos.fill(HIST("hDCAz_1"), dcaZ);
            histos.fill(HIST("hNsigmaTPC_1"), getPIDTPC<IsResoMicrotrack>(track, daughterType));
            histos.fill(HIST("hNsigmaTOF_1"), getPIDTOF<IsResoMicrotrack>(track, daughterType));
            break;
          case 1: // Second
            histos.fill(HIST("hEta_2"), track.eta());
            histos.fill(HIST("hPhi_2"), track.phi());
            histos.fill(HIST("hPt_2"), track.pt());
            histos.fill(HIST("hDCAr_2"), dcaXY);
            histos.fill(HIST("hDCAz_2"), dcaZ);
            histos.fill(HIST("hNsigmaTPC_2"), getPIDTPC<IsResoMicrotrack>(track, daughterType));
            histos.fill(HIST("hNsigmaTOF_2"), getPIDTOF<IsResoMicrotrack>(track, daughterType));
            break;
          case 2: // Third
            histos.fill(HIST("hEta_3"), track.eta());
            histos.fill(HIST("hPhi_3"), track.phi());
            histos.fill(HIST("hPt_3"), track.pt());
            histos.fill(HIST("hDCAr_3"), dcaXY);
            histos.fill(HIST("hDCAz_3"), dcaZ);
            histos.fill(HIST("hNsigmaTPC_3"), getPIDTPC<IsResoMicrotrack>(track, daughterType));
            histos.fill(HIST("hNsigmaTOF_3"), getPIDTOF<IsResoMicrotrack>(track, daughterType));
            break;
          case 3: // Forth
            histos.fill(HIST("hEta_4"), track.eta());
            histos.fill(HIST("hPhi_4"), track.phi());
            histos.fill(HIST("hPt_4"), track.pt());
            histos.fill(HIST("hDCAr_4"), dcaXY);
            histos.fill(HIST("hDCAz_4"), dcaZ);
            histos.fill(HIST("hNsigmaTPC_4"), getPIDTPC<IsResoMicrotrack>(track, daughterType));
            histos.fill(HIST("hNsigmaTOF_4"), getPIDTOF<IsResoMicrotrack>(track, daughterType));
            break;
          default:
            break;
        }
      }
      selectedIndices.push_back(track.index());
    }
    return selectedIndices;
  }

  bool isPairSelected(const TLorentzVector& lv1, const TLorentzVector& lv2, int pairType = 0)
  {
    TLorentzVector lvSum = lv1 + lv2;
    // Mass window cut
    auto pairMass = lvSum.M();
    auto pairMassesLow = PairCuts.cfgPairMassesLow.value;
    auto pairMassesHigh = PairCuts.cfgPairMassesHigh.value;
    if ((pairMassesLow[pairType] > 0) || (pairMassesHigh[pairType] > 0)) {
      if (pairMass < pairMassesLow[pairType] || pairMass > pairMassesHigh[pairType]) {
        return false;
      }
    }
    // Opening angle cut
    double angle = lv1.Vect().Angle(lv2.Vect());
    auto angleCut = PairCuts.cfgPairOACut.value;
    if (angleCut[pairType] > 0) {
      if (angle < angleCut[pairType]) {
        return false;
      }
    }
    if (AnalysisConfig.cfgFillQAPlots) {
      // FIXME: Apply better method
      if (pairType > 0) {
        histos.fill(HIST("hPairInvMass_2"), pairMass);
        histos.fill(HIST("hPairPt_2"), lvSum.Pt());
        histos.fill(HIST("hPairOA_2"), angle);
      } else {
        histos.fill(HIST("hPairInvMass_1"), pairMass);
        histos.fill(HIST("hPairPt_1"), lvSum.Pt());
        histos.fill(HIST("hPairOA_1"), angle);
      }
    }
    return true;
  }

  template <typename TracksType>
  std::vector<std::pair<int, int>>
    getSelectedTrackPairs(const TracksType& dTracks,
                          const std::vector<int>& indicesA,
                          const std::vector<int>& indicesB,
                          float massA, float massB,
                          int pairType)
  {
    std::vector<std::pair<int, int>> selectedPairs;
    TLorentzVector lv1, lv2;
    for (const auto& indexA : indicesA) {
      for (const auto& indexB : indicesB) {
        if (indexA == indexB) {
          continue;
        }

        auto trackA = dTracks.rawIteratorAt(indexA);
        auto trackB = dTracks.rawIteratorAt(indexB);

        lv1.SetXYZM(trackA.px(), trackA.py(), trackA.pz(), massA);
        lv2.SetXYZM(trackB.px(), trackB.py(), trackB.pz(), massB);

        if (!isPairSelected(lv1, lv2, pairType)) {
          continue;
        }
        selectedPairs.emplace_back(indexA, indexB);
      }
    }
    return selectedPairs;
  }

  bool isResoSelected(const TLorentzVector& par1, const TLorentzVector& pair2)
  {
    // Opening angle (3D)
    double oa = par1.Vect().Angle(pair2.Vect());
    if (oa < ResoCuts.cfgPairOALow || oa > ResoCuts.cfgPairOAHigh) {
      return false;
    }
    // Rapidity cut
    TLorentzVector lvTotal = par1 + pair2;
    if (lvTotal.Rapidity() < -0.5 || lvTotal.Rapidity() > 0.5) {
      return false;
    }
    histos.fill(HIST("hResoOA"), oa);
    return true;
  }

  template <bool IsMC, bool IsMix, bool IsResoMicrotrack, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks)
  {
    auto multiplicity = collision.cent();
    histos.fill(HIST("hVertexZ"), collision.posZ());
    histos.fill(HIST("hMultiplicityPercent"), multiplicity);

    auto trackIndices1 = selectTrackIndicesWithQA<IsResoMicrotrack>(dTracks, cfgFirstDaughter, 0);
    auto trackIndices2 = selectTrackIndicesWithQA<IsResoMicrotrack>(dTracks, cfgSecondDaughter, 1);
    auto trackIndices3 = selectTrackIndicesWithQA<IsResoMicrotrack>(dTracks, cfgThirdDaughter, 2);
    auto trackIndices4 = selectTrackIndicesWithQA<IsResoMicrotrack>(dTracks, cfgFourthDaughter, 3);

    // First resconstructed pair
    auto selectedPairs1 = getSelectedTrackPairs(dTracks, trackIndices1, trackIndices2, mass1, mass2, 0);
    auto selectedPairs2 = getSelectedTrackPairs(dTracks, trackIndices3, trackIndices4, mass3, mass4, 1);

    // Resonance loop of selectedPairs1 and selectedPairs2
    for (const auto& pair1 : selectedPairs1) {
      const auto& i1 = pair1.first;
      const auto& i2 = pair1.second;
      for (const auto& pair2 : selectedPairs2) {
        const auto& j1 = pair2.first;
        const auto& j2 = pair2.second;
        // Remove the same track
        if (i1 == j1 || i1 == j2 || i2 == j1 || i2 == j2) {
          continue;
        }
        auto t1 = dTracks.rawIteratorAt(i1);
        auto t2 = dTracks.rawIteratorAt(i2);
        auto t3 = dTracks.rawIteratorAt(j1);
        auto t4 = dTracks.rawIteratorAt(j2);

        TLorentzVector lv1, lv2, lv3, lv4, lvPair1, lvPair2, lvTotal, lResonanceRot;
        lv1.SetXYZM(t1.px(), t1.py(), t1.pz(), mass1);
        lv2.SetXYZM(t2.px(), t2.py(), t2.pz(), mass2);
        lv3.SetXYZM(t3.px(), t3.py(), t3.pz(), mass3);
        lv4.SetXYZM(t4.px(), t4.py(), t4.pz(), mass4);
        lvPair1 = lv1 + lv2;
        lvPair2 = lv3 + lv4;
        if (!isResoSelected(lvPair1, lvPair2))
          continue;
        lvTotal = lv1 + lv2 + lv3 + lv4;
        histos.fill(HIST("h2PairInvMass"), lvPair1.M(), lvPair2.M());
        histos.fill(HIST("hResoInvMass"), lvTotal.M());
        histos.fill(HIST("THnResoInvMass"), multiplicity, lvTotal.Pt(), lvTotal.M());
        if (BkgEstimationConfig.cfgFillRotBkg) {
          for (int i = 0; i < BkgEstimationConfig.cfgNrotBkg; i++) {
            auto lRotAngle = BkgEstimationConfig.cfgMinRot + i * ((BkgEstimationConfig.cfgMaxRot - BkgEstimationConfig.cfgMinRot) / (BkgEstimationConfig.cfgNrotBkg - 1));
            lvPair2.RotateZ(lRotAngle);
            lResonanceRot = lvPair1 + lvPair2;
            histos.fill(HIST("THnResoInvMassBkg"), multiplicity, lResonanceRot.Pt(), lResonanceRot.M());
          }
        }
      }
    }
  }

  void processDummy(aod::ResoCollision const& /*collisions*/)
  {
  }
  PROCESS_SWITCH(DoubleResonanceScan, processDummy, "Process Dummy", true);

  void processResoTracks(aod::ResoCollision const& collision, aod::ResoTracks const& resotracks)
  {
    fillHistograms<false, false, false>(collision, resotracks);
  }
  PROCESS_SWITCH(DoubleResonanceScan, processResoTracks, "Process ResoTracks", false);

  void processResoMicroTracks(aod::ResoCollision const& collision, aod::ResoMicroTracks const& resomicrotracks)
  {
    fillHistograms<false, false, true>(collision, resomicrotracks);
  }
  PROCESS_SWITCH(DoubleResonanceScan, processResoMicroTracks, "Process ResoMicroTracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<DoubleResonanceScan>(cfgc)}; }
