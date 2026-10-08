// Copyright 2019-2020 CERN and copyright holders of the ALICE Collaboration.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

/// \file f0980FlowSP.cxx
/// \brief pi+pi- invariant-mass analysis and f0(980) v2 with the scalar-product method

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h>

#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
using P4 = ROOT::Math::PxPyPzMVector;
constexpr float kInvalidQ = -900.f;
} // namespace

struct F0980FlowSP {
  SliceCache cache;
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Pair and daughter selection. Defaults follow rho770analysis.cxx.
  Configurable<float> cfgMinPt{"cfgMinPt", 0.15f, "Minimum daughter transverse momentum"};
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8f, "Maximum daughter pseudorapidity"};
  Configurable<float> cfgMinRap{"cfgMinRap", -0.5f, "Minimum pair rapidity"};
  Configurable<float> cfgMaxRap{"cfgMaxRap", 0.5f, "Maximum pair rapidity"};
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Require primary track"};
  Configurable<bool> cfgGlobalWoDCA{"cfgGlobalWoDCA", true, "Require global track without DCA cuts"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "Require PV contributor"};
  Configurable<float> cfgMaxDCAxy{"cfgMaxDCAxy", 0.1f, "Maximum |DCAxy| (cm)"};
  Configurable<float> cfgMaxDCAz{"cfgMaxDCAz", 0.1f, "Maximum |DCAz| (cm)"};
  Configurable<float> cfgMaxTPCChi2NCl{"cfgMaxTPCChi2NCl", 4.f, "Maximum TPC chi2 per cluster"};
  Configurable<float> cfgMaxITSChi2NCl{"cfgMaxITSChi2NCl", 36.f, "Maximum ITS chi2 per cluster"};
  Configurable<int> cfgMinTPCRows{"cfgMinTPCRows", 70, "Minimum TPC crossed rows"};
  Configurable<int> cfgMinITSClusters{"cfgMinITSClusters", 5, "Minimum ITS clusters"};

  Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "Require sel8"};
  Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", true, "Require kNoSameBunchPileup"};
  Configurable<bool> cfgRequireNoTimeFrameBorder{"cfgRequireNoTimeFrameBorder", true, "Require kNoTimeFrameBorder"};
  Configurable<bool> cfgRequireNoITSROFrameBorder{"cfgRequireNoITSROFrameBorder", true, "Require kNoITSROFrameBorder"};
  Configurable<bool> cfgRequireITSTPCVertex{"cfgRequireITSTPCVertex", true, "Require kIsVertexITSTPC"};

  Configurable<int> cfgPIDMode{"cfgPIDMode", 3, "Pion PID: 0 TPC-or-TOF, 1 TPC only, 2 combined, 3 TPC with TOF veto"};
  Configurable<float> cfgTPCNSigmaPi{"cfgTPCNSigmaPi", 5.f, "TPC pion n-sigma with TOF"};
  Configurable<float> cfgTPCNSigmaPiNoTOF{"cfgTPCNSigmaPiNoTOF", 2.f, "TPC pion n-sigma without TOF"};
  Configurable<float> cfgTOFNSigmaPi{"cfgTOFNSigmaPi", 3.f, "TOF pion n-sigma"};
  Configurable<float> cfgCombinedNSigmaPi{"cfgCombinedNSigmaPi", 3.f, "Combined TPC-TOF pion n-sigma"};
  Configurable<uint64_t> cfgPairLogInterval{"cfgPairLogInterval", 100000, "Print one accepted pi-pi pair every N accepted pairs (0 disables logging)"};

  // The vector tables contain one entry per harmonic requested from qVectorsTable.
  Configurable<int> cfgHarmonicIndex{"cfgHarmonicIndex", 0, "Index of n=2 in Q-vector vector columns"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0.f, "Minimum FT0C centrality for OO"};
  Configurable<float> cfgCentMax{"cfgCentMax", 100.f, "Maximum FT0C centrality for OO"};

  ConfigurableAxis axisMass{"axisMass", {400, 0.2, 2.2}, "m_{pi pi} (GeV/c^{2})"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0., 0.5, 1., 2., 3., 4., 6., 10., 13., 20.}, "pair pT (GeV/c)"};
  ConfigurableAxis axisCent{"axisCent", {VARIABLE_WIDTH, 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 100.}, "FT0C centrality (%)"};
  ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0., 1., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 105.}, "multiplicity percentile"};

  uint64_t acceptedPairCount = 0;

  using PPCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::PVMults>;
  using OOCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs,
                                 aod::QvectorFT0CVecs, aod::QvectorFT0AVecs,
                                 aod::QvectorTPCposVecs>;
  using TracksTPC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                              aod::TrackSelection, aod::pidTPCPi>;
  using TracksTPCTOF = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                 aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi>;

  void init(InitContext const&)
  {
    // Counts and weighted SP numerators have identical axes. Their ratio is
    // <u_2 q_A*> in each (mass, pT, centrality) bin.
    registry.add("PP/hMassUS", "pp unlike-sign;mass;pT;multiplicity", HistType::kTHnSparseF, {axisMass, axisPt, axisMult});
    registry.add("PP/hMassLSpp", "pp ++;mass;pT;multiplicity", HistType::kTHnSparseF, {axisMass, axisPt, axisMult});
    registry.add("PP/hMassLSmm", "pp --;mass;pT;multiplicity", HistType::kTHnSparseF, {axisMass, axisPt, axisMult});

    registry.add("OO/hMassUS", "OO unlike-sign;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("OO/hMassLSpp", "OO ++;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("OO/hMassLSmm", "OO --;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("OO/hSpNumUS", "sum Re(u2 qFT0C*), unlike-sign;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("OO/hSpNumLSpp", "sum Re(u2 qFT0C*), ++;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("OO/hSpNumLSmm", "sum Re(u2 qFT0C*), --;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});

    registry.add("OO/pQAB", "<qFT0C qFT0A*>;centrality", HistType::kTProfile, {axisCent});
    registry.add("OO/pQAC", "<qFT0C qTPCpos*>;centrality", HistType::kTProfile, {axisCent});
    registry.add("OO/pQBC", "<qFT0A qTPCpos*>;centrality", HistType::kTProfile, {axisCent});
    registry.add("OO/pQABIm", "Im(qFT0C qFT0A*);centrality", HistType::kTProfile, {axisCent});
    registry.add("OO/pQACIm", "Im(qFT0C qTPCpos*);centrality", HistType::kTProfile, {axisCent});
    registry.add("OO/pQBCIm", "Im(qFT0A qTPCpos*);centrality", HistType::kTProfile, {axisCent});

    registry.add("QA/hEventCounter", "event counter", HistType::kTH1F, {{8, -0.5, 7.5}});
    registry.add("QA/hQCalibration", "Q-vector calibration flag", HistType::kTH1F, {{2, -0.5, 1.5}});
    registry.add("QA/hQFT0C", "corrected q_{2,FT0C};q_{x};q_{y}", HistType::kTH2F, {{200, -1., 1.}, {200, -1., 1.}});
    registry.add("QA/TrackSelected/hPt", "track-selected candidates;p_{T} (GeV/c);counts", HistType::kTH1F, {{200, 0., 20.}});
    registry.add("QA/TrackSelected/hEta", "track-selected candidates;#eta;counts", HistType::kTH1F, {{160, -0.8, 0.8}});
    registry.add("QA/TrackSelected/hDCAxy", "track-selected candidates;DCA_{xy} (cm);counts", HistType::kTH1F, {{200, -0.1, 0.1}});
    registry.add("QA/TrackSelected/hDCAz", "track-selected candidates;DCA_{z} (cm);counts", HistType::kTH1F, {{200, -0.1, 0.1}});
    registry.add("QA/TrackSelected/hTPCCrossedRows", "track-selected candidates;TPC crossed rows;counts", HistType::kTH1F, {{165, -0.5, 164.5}});
    registry.add("QA/TrackSelected/hTPCNSigmaPiVsP", "track-selected candidates before pion PID;p (GeV/c);n#sigma_{TPC}^{#pi}", HistType::kTH2F, {{200, 0., 20.}, {120, -6., 6.}});
    registry.add("QA/TrackSelected/hTPCSignalVsP", "track-selected candidates before pion PID;p (GeV/c);TPC dE/dx", HistType::kTH2F, {{200, 0., 20.}, {500, 0., 500.}});
    registry.add("QA/SelectedPions/hTPCNSigmaPiVsP", "selected pion candidates;p (GeV/c);n#sigma_{TPC}^{#pi}", HistType::kTH2F, {{200, 0., 20.}, {120, -6., 6.}});
    registry.add("QA/SelectedPions/hTPCSignalVsP", "selected pion candidates after pion PID;p (GeV/c);TPC dE/dx", HistType::kTH2F, {{200, 0., 20.}, {500, 0., 500.}});
  }

  template <typename T>
  bool selectTrack(T const& track) const
  {
    if (track.pt() < cfgMinPt.value || std::abs(track.eta()) > cfgMaxEta.value) {
      return false;
    }
    if (std::abs(track.dcaXY()) >= cfgMaxDCAxy.value || std::abs(track.dcaZ()) >= cfgMaxDCAz.value) {
      return false;
    }
    if (cfgPrimaryTrack.value && !track.isPrimaryTrack()) {
      return false;
    }
    if (cfgGlobalWoDCA.value && !track.isGlobalTrackWoDCA()) {
      return false;
    }
    if (cfgPVContributor.value && !track.isPVContributor()) {
      return false;
    }
    return track.tpcChi2NCl() < cfgMaxTPCChi2NCl.value &&
           track.itsChi2NCl() < cfgMaxITSChi2NCl.value &&
           track.tpcNClsCrossedRows() > cfgMinTPCRows.value &&
           track.itsNCls() > cfgMinITSClusters.value;
  }

  template <typename C>
  bool selectEvent(C const& collision) const
  {
    return (!cfgRequireSel8.value || collision.sel8()) &&
           (!cfgRequireNoSameBunchPileup.value || collision.selection_bit(aod::evsel::kNoSameBunchPileup)) &&
           (!cfgRequireNoTimeFrameBorder.value || collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) &&
           (!cfgRequireNoITSROFrameBorder.value || collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) &&
           (!cfgRequireITSTPCVertex.value || collision.selection_bit(aod::evsel::kIsVertexITSTPC));
  }

  template <typename T>
  bool selectPion(T const& track) const
  {
    const float tpc = std::abs(track.tpcNSigmaPi());
    const float tof = std::abs(track.tofNSigmaPi());
    switch (cfgPIDMode.value) {
      case 0:
        return tpc < cfgTPCNSigmaPi.value || (track.hasTOF() && tof < cfgTOFNSigmaPi.value);
      case 1:
        return tpc < cfgTPCNSigmaPi.value;
      case 2:
        return track.hasTOF() && std::hypot(tpc, tof) < cfgCombinedNSigmaPi.value;
      case 3:
        return track.hasTOF() ? (tpc < cfgTPCNSigmaPi.value && tof < cfgTOFNSigmaPi.value) : (tpc < cfgTPCNSigmaPiNoTOF.value);
      default:
        return false;
    }
  }

  template <typename T>
  bool selectPionTPCOnly(T const& track) const
  {
    return std::abs(track.tpcNSigmaPi()) < cfgTPCNSigmaPiNoTOF.value;
  }

  template <typename Tracks, typename TrackSelect, typename PionSelect, typename FillPair>
  void loopPairsSelected(Tracks const& tracks, TrackSelect&& trackSelect, PionSelect&& pionSelect, FillPair&& fillPair)
  {
    struct SelectedPion {
      float px;
      float py;
      float pz;
      int8_t sign;
    };

    // Apply daughter-quality and PID selections once per track before building
    // pairs.  In OO events this avoids constructing all unfiltered O(N^2)
    // combinations and repeatedly evaluating the same track cuts.
    std::vector<SelectedPion> pions;
    pions.reserve(tracks.size());
    for (auto const& track : tracks) {
      if (!trackSelect(track)) {
        continue;
      }
      registry.fill(HIST("QA/TrackSelected/hPt"), track.pt());
      registry.fill(HIST("QA/TrackSelected/hEta"), track.eta());
      registry.fill(HIST("QA/TrackSelected/hDCAxy"), track.dcaXY());
      registry.fill(HIST("QA/TrackSelected/hDCAz"), track.dcaZ());
      registry.fill(HIST("QA/TrackSelected/hTPCCrossedRows"), track.tpcNClsCrossedRows());
      registry.fill(HIST("QA/TrackSelected/hTPCNSigmaPiVsP"), track.p(), track.tpcNSigmaPi());
      registry.fill(HIST("QA/TrackSelected/hTPCSignalVsP"), track.p(), track.tpcSignal());
      if (!pionSelect(track)) {
        continue;
      }
      registry.fill(HIST("QA/SelectedPions/hTPCNSigmaPiVsP"), track.p(), track.tpcNSigmaPi());
      registry.fill(HIST("QA/SelectedPions/hTPCSignalVsP"), track.p(), track.tpcSignal());
      pions.push_back({track.px(), track.py(), track.pz(), static_cast<int8_t>(track.sign())});
    }

    const double mpi = o2::constants::physics::MassPionCharged;
    for (size_t i = 0; i < pions.size(); ++i) {
      for (size_t j = i + 1; j < pions.size(); ++j) {
      auto const& track1 = pions[i];
      auto const& track2 = pions[j];
      P4 pion1(track1.px, track1.py, track1.pz, mpi);
      P4 pion2(track2.px, track2.py, track2.pz, mpi);
      auto pair = pion1 + pion2;
      if (pair.Rapidity() < cfgMinRap.value || pair.Rapidity() > cfgMaxRap.value) {
        continue;
      }
      ++acceptedPairCount;
      if (cfgPairLogInterval.value > 0 && acceptedPairCount % cfgPairLogInterval.value == 0) {
        LOGF(info, "pi-pi pair progress: accepted=%llu, mass=%.6f GeV/c^2, pt=%.6f GeV/c, signs=(%d,%d)",
             static_cast<unsigned long long>(acceptedPairCount), pair.M(), pair.Pt(),
             static_cast<int>(track1.sign), static_cast<int>(track2.sign));
      }
      fillPair(track1.sign, track2.sign, pair);
      }
    }
  }

  template <typename Tracks, typename FillPair>
  void loopPairs(Tracks const& tracks, FillPair&& fillPair)
  {
    loopPairsSelected(tracks, [&](auto const& track) { return selectTrack(track); },
                      [&](auto const& track) { return selectPion(track); },
                      std::forward<FillPair>(fillPair));
  }

  void processPP(PPCollisions::iterator const& collision, TracksTPC const& tracks)
  {
    registry.fill(HIST("QA/hEventCounter"), 0.);
    if (!selectEvent(collision) || !collision.isInelGt0()) {
      return;
    }
    registry.fill(HIST("QA/hEventCounter"), 1.);
    const float multiplicity = collision.multNTracksPVeta1();
    loopPairsSelected(tracks, [&](auto const& track) { return selectTrack(track); },
                      [&](auto const& track) { return selectPionTPCOnly(track); },
                      [&](int8_t sign1, int8_t sign2, P4 const& pair) {
      if (sign1 * sign2 < 0) {
        registry.fill(HIST("PP/hMassUS"), pair.M(), pair.Pt(), multiplicity);
      } else if (sign1 > 0 && sign2 > 0) {
        registry.fill(HIST("PP/hMassLSpp"), pair.M(), pair.Pt(), multiplicity);
      } else if (sign1 < 0 && sign2 < 0) {
        registry.fill(HIST("PP/hMassLSmm"), pair.M(), pair.Pt(), multiplicity);
      }
                      });
  }
  PROCESS_SWITCH(F0980FlowSP, processPP, "Process pp invariant-mass analysis", false);

  void processOO(OOCollisions::iterator const& collision, TracksTPCTOF const& tracks)
  {
    registry.fill(HIST("QA/hEventCounter"), 2.);
    const float centrality = collision.centFT0C();
    if (!selectEvent(collision) || centrality < cfgCentMin || centrality >= cfgCentMax) {
      return;
    }

    const int ih = cfgHarmonicIndex;
    if (ih < 0 || ih >= static_cast<int>(collision.qvecFT0CReVec().size()) ||
        ih >= static_cast<int>(collision.qvecFT0AReVec().size()) ||
        ih >= static_cast<int>(collision.qvecTPCposReVec().size())) {
      return;
    }
    const bool calibrated = collision.isCalibrated();
    registry.fill(HIST("QA/hQCalibration"), calibrated);
    if (!calibrated) {
      return;
    }

    const double qAx = collision.qvecFT0CReVec()[ih];
    const double qAy = collision.qvecFT0CImVec()[ih];
    const double qBx = collision.qvecFT0AReVec()[ih];
    const double qBy = collision.qvecFT0AImVec()[ih];
    const double qCx = collision.qvecTPCposReVec()[ih];
    const double qCy = collision.qvecTPCposImVec()[ih];
    if (qAx < kInvalidQ || qAy < kInvalidQ || qBx < kInvalidQ || qBy < kInvalidQ || qCx < kInvalidQ || qCy < kInvalidQ) {
      return;
    }

    registry.fill(HIST("QA/hEventCounter"), 3.);
    registry.fill(HIST("QA/hQFT0C"), qAx, qAy);
    registry.fill(HIST("OO/pQAB"), centrality, qAx * qBx + qAy * qBy);
    registry.fill(HIST("OO/pQAC"), centrality, qAx * qCx + qAy * qCy);
    registry.fill(HIST("OO/pQBC"), centrality, qBx * qCx + qBy * qCy);
    registry.fill(HIST("OO/pQABIm"), centrality, qAy * qBx - qAx * qBy);
    registry.fill(HIST("OO/pQACIm"), centrality, qAy * qCx - qAx * qCy);
    registry.fill(HIST("OO/pQBCIm"), centrality, qBy * qCx - qBx * qCy);

    loopPairs(tracks, [&](int8_t sign1, int8_t sign2, P4 const& pair) {
      const double ux = std::cos(2. * pair.Phi());
      const double uy = std::sin(2. * pair.Phi());
      const double sp = ux * qAx + uy * qAy;
      if (sign1 * sign2 < 0) {
        registry.fill(HIST("OO/hMassUS"), pair.M(), pair.Pt(), centrality);
        registry.fill(HIST("OO/hSpNumUS"), pair.M(), pair.Pt(), centrality, sp);
      } else if (sign1 > 0 && sign2 > 0) {
        registry.fill(HIST("OO/hMassLSpp"), pair.M(), pair.Pt(), centrality);
        registry.fill(HIST("OO/hSpNumLSpp"), pair.M(), pair.Pt(), centrality, sp);
      } else if (sign1 < 0 && sign2 < 0) {
        registry.fill(HIST("OO/hMassLSmm"), pair.M(), pair.Pt(), centrality);
        registry.fill(HIST("OO/hSpNumLSmm"), pair.M(), pair.Pt(), centrality, sp);
      }
    });
  }
  PROCESS_SWITCH(F0980FlowSP, processOO, "Process OO invariant-mass scalar-product analysis", false);

  void processOOTPCOnly(OOCollisions::iterator const& collision, TracksTPC const& tracks)
  {
    registry.fill(HIST("QA/hEventCounter"), 2.);
    const float centrality = collision.centFT0C();
    if (!selectEvent(collision) || centrality < cfgCentMin || centrality >= cfgCentMax) {
      return;
    }

    const int ih = cfgHarmonicIndex;
    if (ih < 0 || ih >= static_cast<int>(collision.qvecFT0CReVec().size()) ||
        ih >= static_cast<int>(collision.qvecFT0AReVec().size()) ||
        ih >= static_cast<int>(collision.qvecTPCposReVec().size()) || !collision.isCalibrated()) {
      return;
    }
    const double qAx = collision.qvecFT0CReVec()[ih];
    const double qAy = collision.qvecFT0CImVec()[ih];
    const double qBx = collision.qvecFT0AReVec()[ih];
    const double qBy = collision.qvecFT0AImVec()[ih];
    const double qCx = collision.qvecTPCposReVec()[ih];
    const double qCy = collision.qvecTPCposImVec()[ih];
    if (qAx < kInvalidQ || qAy < kInvalidQ || qBx < kInvalidQ || qBy < kInvalidQ || qCx < kInvalidQ || qCy < kInvalidQ) {
      return;
    }

    registry.fill(HIST("QA/hEventCounter"), 3.);
    registry.fill(HIST("QA/hQCalibration"), true);
    registry.fill(HIST("QA/hQFT0C"), qAx, qAy);
    registry.fill(HIST("OO/pQAB"), centrality, qAx * qBx + qAy * qBy);
    registry.fill(HIST("OO/pQAC"), centrality, qAx * qCx + qAy * qCy);
    registry.fill(HIST("OO/pQBC"), centrality, qBx * qCx + qBy * qCy);
    registry.fill(HIST("OO/pQABIm"), centrality, qAy * qBx - qAx * qBy);
    registry.fill(HIST("OO/pQACIm"), centrality, qAy * qCx - qAx * qCy);
    registry.fill(HIST("OO/pQBCIm"), centrality, qBy * qCx - qBx * qCy);

    loopPairsSelected(tracks, [&](auto const& track) { return selectTrack(track); },
                      [&](auto const& track) { return selectPionTPCOnly(track); },
                      [&](int8_t sign1, int8_t sign2, P4 const& pair) {
                        const double sp = std::cos(2. * pair.Phi()) * qAx + std::sin(2. * pair.Phi()) * qAy;
                        if (sign1 * sign2 < 0) {
                          registry.fill(HIST("OO/hMassUS"), pair.M(), pair.Pt(), centrality);
                          registry.fill(HIST("OO/hSpNumUS"), pair.M(), pair.Pt(), centrality, sp);
                        } else if (sign1 > 0 && sign2 > 0) {
                          registry.fill(HIST("OO/hMassLSpp"), pair.M(), pair.Pt(), centrality);
                          registry.fill(HIST("OO/hSpNumLSpp"), pair.M(), pair.Pt(), centrality, sp);
                        } else if (sign1 < 0 && sign2 < 0) {
                          registry.fill(HIST("OO/hMassLSmm"), pair.M(), pair.Pt(), centrality);
                          registry.fill(HIST("OO/hSpNumLSmm"), pair.M(), pair.Pt(), centrality, sp);
                        }
                      });
  }
  PROCESS_SWITCH(F0980FlowSP, processOOTPCOnly, "Process OO with TPC-only pion PID", true);

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<F0980FlowSP>(cfgc, TaskName{"lf-f0980-flow-sp"})};
}
