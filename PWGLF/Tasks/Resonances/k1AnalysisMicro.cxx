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
/// \file k1AnalysisMicro.cxx
/// \brief Reconstruction of track-track decay resonance candidates
/// \author Su-Jeong Ji <su-jeong.ji@cern.ch>, Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///

#include "PWGLF/DataModel/LFResonanceTables.h"

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
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
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/VectorUtil.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <set>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
using namespace o2::constants::math;

struct K1AnalysisMicro {
  // Module-initializer v001 tables; full tracks keep their unversioned schema as a fallback.
  using ResoCollisions = aod::ResoCollisions_001;
  using ResoMCCols = soa::Join<ResoCollisions, aod::ResoMCCollisions_001>;
  using ResoTracks = aod::ResoTracks; // no v001 exists; K1 does not need ResoTrackTracks (trackId unused)
  using ResoMicroTracks = aod::ResoMicroTracks_001;
  using ResoMCTracks = soa::Join<ResoTracks, aod::ResoMCTracks>;
  using ResoMCMicroTracks = soa::Join<ResoMicroTracks, aod::ResoMCMicroTracks_001>;
  using ResoMCParents = aod::ResoMCParents_001;

  enum BinAnti : unsigned int {
    kNormal = 0,
    kAnti,
    kNAEnd
  };
  enum BinType : unsigned int {
    kK1P = 0,
    kK1N,
    kK1P_Mix,
    kK1N_Mix,
    kK1P_GenINEL10,
    kK1N_GenINEL10,
    kK1P_GenINELgt10,
    kK1N_GenINELgt10,
    kK1P_GenTrig10,
    kK1N_GenTrig10,
    kK1P_GenEvtSel,
    kK1N_GenEvtSel,
    kK1P_Rec,
    kK1N_Rec,
    kTYEnd
  };
  enum class Species : int {
    Pion = 0,
    Kaon = 1
  };
  // Last stage passed by a track; the cut-flow histogram is filled directly from this value.
  enum TrackStage : int {
    kTrkInput = 0,
    kTrkPt,
    kTrkEta,
    kTrkDCAxy,
    kTrkDCAz,
    kTrkFlags,
    kTrkClusters,
    kTrkTOFRequired,
    kTrkPID,
    kTrkNStages
  };
  enum class QAFolder {
    Before, // QA/*: before the candidate cuts
    After,  // QAcut/*: after the candidate cuts
    MC      // QAMC/*: matched K1 truth candidates
  };
  // Resolved PID cut of one species at a given pT.
  struct PIDCut {
    double tpcMax = 0.;
    double tofMax = 0.;
    double combined = 0.;
    bool tofRequired = false;
  };

  static constexpr float DisabledCut = -999.f; // an optional cut with this value is off and not evaluated
  static constexpr double MassRho770 = 0.77526; // PDG 2024, not available in o2::constants::physics
  static constexpr double DCAGridStep = 0.025;  // v001 micro DCA encoding, lower-inclusive bins up to DCAGridMax
  static constexpr double DCAGridMax = 0.15;
  static constexpr double PIDGridStart = 2.0; // v001 micro nSigma encoding: 0.25 bins in [2.0, 3.5]
  static constexpr double PIDGridStep = 0.25;
  static constexpr double PIDGridMax = 3.5;
  static constexpr double GridTolerance = 1e-4;
  static constexpr int NCandidateStages = 12;

  SliceCache cache;
  // Registered only to enable the slice cache that SameKindPair (event mixing) needs, as in Xi1820Analysis
  Preslice<ResoTracks> perResoCollisionTrack = aod::resodaughter::resoCollisionId;
  Preslice<ResoMicroTracks> perResoCollisionMicroTrack = aod::resodaughter::resoCollisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Configurable<int> cfgTruthDebug{"cfgTruthDebug", 0, "Maximum logged matched candidates per truth channel"};
  std::array<int, 3> truthDebugCounts{};

  //// Configurables
  Configurable<int> cNbinsDiv{"cNbinsDiv", 1, "Integer to divide the number of bins"};
  /// Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 5, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 20.0f, 40.0f, 60.0f, 80.0f, 100.0f, 200.0f, 99999.f}, "Mixing bins - multiplicity"};

  // Event selection (a group without prefix keeps the plain key names)
  struct : ConfigurableGroup {
    Configurable<bool> cRecoINELgt0{"cRecoINELgt0", false, "Apply reconstructed INEL>0 selection"};
    Configurable<bool> cMCINELgt0{"cMCINELgt0", false, "Require generator INEL>0 in MC processes"};
    Configurable<bool> cMCVtxIn10{"cMCVtxIn10", false, "Require generator |vz| < 10 cm in MC processes"};
  } eventCuts;

  /// Track selections (common for pion and kaon, -999 switches an optional cut off)
  struct : ConfigurableGroup {
    Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
    Configurable<float> cMaxEtacut{"cMaxEtacut", -999.f, "Track maximum |eta| cut (-999: off)"};
    // DCAr to PV
    Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.1, "Track DCAr cut to PV Maximum"};
    // DCAz to PV
    Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 0.1, "Track DCAz cut to PV Maximum"};
    Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};
    Configurable<bool> cfgUsePtDepDCA{"cfgUsePtDepDCA", false, "Use pT dependent DCA cut instead of the fixed maximum"};
    Configurable<float> cDCAToPVByPtP0{"cDCAToPVByPtP0", 0.004f, "pT dependent DCA cut = P0 + coefficient / pT^power (cm)"};
    Configurable<float> cDCAToPVByPtCoeff{"cDCAToPVByPtCoeff", 0.013f, "Coefficient in the pT dependent DCA cut"};
    Configurable<float> cDCAToPVByPtPower{"cDCAToPVByPtPower", 1.f, "Power in the pT dependent DCA cut"};
    Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
    Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
    Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                      // kGoldenChi2 | kDCAxy | kDCAz
    Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};          // PV Contriuibutor
    Configurable<bool> cfgUseTPCRefit{"cfgUseTPCRefit", false, "Require TPC Refit"};
    Configurable<bool> cfgUseITSRefit{"cfgUseITSRefit", false, "Require ITS Refit"};
    Configurable<int> cfgTPCcluster{"cfgTPCcluster", 0, "Number of TPC cluster (found clusters, ResoTracks only)"};
    Configurable<int> cfgTPCCrossedRowsMin{"cfgTPCCrossedRowsMin", 0, "Minimum number of TPC crossed rows"};
    Configurable<int> cfgITSNClsMin{"cfgITSNClsMin", 0, "Minimum number of ITS clusters (ResoMicroTracks only)"};
    Configurable<bool> cfgHasTOF{"cfgHasTOF", false, "Require TOF"};
  } trackCuts;

  /// PID Selections
  Configurable<bool> cByPassTOF{"cByPassTOF", false, "Bypass the TOF nSigma selection"};
  struct : ConfigurableGroup {
    Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 3.0, "TPC nSigma cut for Pion (-999: off)"};    // TPC
    Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion (-999: off)"};    // TOF
    Configurable<double> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", -999, "Combined nSigma cut for Pion"}; // Combined
    Configurable<bool> cUseOnlyTOFTrackPi{"cUseOnlyTOFTrackPi", false, "Use only TOF track for PID selection"}; // Use only TOF track for Pion PID selection
    Configurable<bool> cPionUsePtDepPID{"cPionUsePtDepPID", false, "Use pT-dependent PID cuts for pion"};
    Configurable<std::vector<float>> cPionPIDPtBins{"cPionPIDPtBins", {0.0f, 0.5f, 0.8f, 2.0f, 999.0f}, "pT bin edges for pion PID cuts"};
    Configurable<std::vector<float>> cPionTPCNSigmaCuts{"cPionTPCNSigmaCuts", {3.0f, 3.0f, 2.0f, 2.0f}, "TPC NSigma cuts per pT bin (pion)"};
    Configurable<std::vector<float>> cPionTOFNSigmaCuts{"cPionTOFNSigmaCuts", {3.0f, 3.0f, 3.0f, 3.0f}, "TOF NSigma cuts per pT bin (pion)"};
    Configurable<std::vector<int>> cPionTOFRequired{"cPionTOFRequired", {0, 0, 1, 1}, "Require TOF per pT bin (pion)"};
  } pionPID;
  struct : ConfigurableGroup {
    Configurable<double> cMaxTPCnSigmaKaon{"cMaxTPCnSigmaKaon", 3.0, "TPC nSigma cut for Kaon (-999: off)"};    // TPC
    Configurable<double> cMaxTOFnSigmaKaon{"cMaxTOFnSigmaKaon", 3.0, "TOF nSigma cut for Kaon (-999: off)"};    // TOF
    Configurable<double> nsigmaCutCombinedKaon{"nsigmaCutCombinedKaon", -999, "Combined nSigma cut for Kaon"}; // Combined
    Configurable<bool> cUseOnlyTOFTrackKa{"cUseOnlyTOFTrackKa", false, "Use only TOF track for PID selection"}; // Use only TOF track for Kaon PID selection
    Configurable<bool> cKaonUsePtDepPID{"cKaonUsePtDepPID", false, "Use pT-dependent PID cuts for kaon"};
    Configurable<std::vector<float>> cKaonPIDPtBins{"cKaonPIDPtBins", {0.0f, 0.5f, 0.8f, 2.0f, 999.0f}, "pT bin edges for kaon PID cuts"};
    Configurable<std::vector<float>> cKaonTPCNSigmaCuts{"cKaonTPCNSigmaCuts", {3.0f, 3.0f, 2.0f, 2.0f}, "TPC NSigma cuts per pT bin (kaon)"};
    Configurable<std::vector<float>> cKaonTOFNSigmaCuts{"cKaonTOFNSigmaCuts", {3.0f, 3.0f, 3.0f, 3.0f}, "TOF NSigma cuts per pT bin (kaon)"};
    Configurable<std::vector<int>> cKaonTOFRequired{"cKaonTOFRequired", {0, 0, 1, 1}, "Require TOF per pT bin (kaon)"};
  } kaonPID;

  Configurable<bool> additionalQAplots{"additionalQAplots", true, "Additional QA plots"};

  // Secondary selection (-999 switches a cut off; the values it needs are then not computed)
  struct : ConfigurableGroup {
    Configurable<double> cMinSecondaryPtCut{"cMinSecondaryPtCut", 0.5, "Min pT cut for secondary selection"};
    Configurable<bool> cfgModeK892orRho{"cfgModeK892orRho", false, "Secondary scenario for K892 (true) or Rho (false)"};
    Configurable<double> cSecondaryMasswindow{"cSecondaryMasswindow", -999, "Secondary inv mass selection window"};
    Configurable<double> cMinAnotherSecondaryMassCut{"cMinAnotherSecondaryMassCut", -999, "Min inv. mass selection of another secondary scenario"};
    Configurable<double> cMaxAnotherSecondaryMassCut{"cMaxAnotherSecondaryMassCut", -999, "MAx inv. mass selection of another secondary scenario"};
    Configurable<double> cMinPiKaMassCut{"cMinPiKaMassCut", -999, "bPion-Kaon pair inv mass selection minimum"};
    Configurable<double> cMaxPiKaMassCut{"cMaxPiKaMassCut", -999, "bPion-Kaon pair inv mass selection maximum"};
    Configurable<double> cMinAngle{"cMinAngle", -999, "Minimum angle between the secondary resonance and the bachelor"};
    Configurable<double> cMaxAngle{"cMaxAngle", -999, "Maximum angle between the secondary resonance and the bachelor"};
    Configurable<double> cMinPairAsym{"cMinPairAsym", -999, "Minimum pair asymmetry"};
    Configurable<double> cMaxPairAsym{"cMaxPairAsym", -999, "Maximum pair asymmetry"};
  } secondaryCuts;

  // K1 selection
  Configurable<double> cK1MaxRap{"cK1MaxRap", 0.5, "K1 maximum rapidity"};
  Configurable<double> cK1MinRap{"cK1MinRap", -0.5, "K1 minimum rapidity"};

  // A cut is on unless it carries the disabled value (tolerant to the float parsing of the JSON value).
  static bool isCutEnabled(float value)
  {
    return value > DisabledCut + 1.f;
  }

  // v001 micro values are lower-inclusive bin edges: a maximum cut on the grid keeps bins below it.
  static bool passesBinnedMax(double decoded, double cut)
  {
    return decoded < cut - Epsilon;
  }

  // Minimum cut on the grid keeps the bin starting at the cut.
  static bool passesBinnedMin(double decoded, double cut)
  {
    return decoded >= cut - Epsilon;
  }

  template <bool IsResoMicrotrack>
  static bool passesMax(double value, double cut)
  {
    if constexpr (IsResoMicrotrack) {
      return passesBinnedMax(value, cut);
    } else {
      return value < cut;
    }
  }

  static bool isInRange(double value, double minimum, double maximum)
  {
    if (isCutEnabled(minimum) && value < minimum) {
      return false;
    }
    if (isCutEnabled(maximum) && value > maximum) {
      return false;
    }
    return true;
  }

  static bool isInWindow(double value, double center, double width)
  {
    return std::abs(value - center) < width;
  }

  // Preserve pT-bin membership [low, high).
  static int getPtBinIndex(float pt, const std::vector<float>& ptBins)
  {
    for (std::size_t i = 1; i < ptBins.size(); ++i) {
      if (pt >= ptBins[i - 1] && pt < ptBins[i]) {
        return static_cast<int>(i - 1);
      }
    }
    return -1;
  }

  // Derived once in init(): which candidate cuts are switched on.
  bool secondaryWindowOn = false;
  bool anotherMassCutOn = false;
  bool piKaMassCutOn = false;
  bool angleCutOn = false;
  bool pairAsymCutOn = false;

  void init(o2::framework::InitContext&)
  {
    const int sameEventModes = static_cast<int>(doprocessResoTracks) + static_cast<int>(doprocessResoMicroTracks) +
                               static_cast<int>(doprocessMC) + static_cast<int>(doprocessMCMicro);
    const int mixedEventModes = static_cast<int>(doprocessME) + static_cast<int>(doprocessMEMicro);
    if (sameEventModes > 1 || mixedEventModes > 1) {
      LOG(fatal) << "Enable at most one same-event mode and one mixing mode";
    }

    secondaryWindowOn = isCutEnabled(secondaryCuts.cSecondaryMasswindow);
    anotherMassCutOn = isCutEnabled(secondaryCuts.cMinAnotherSecondaryMassCut) || isCutEnabled(secondaryCuts.cMaxAnotherSecondaryMassCut);
    piKaMassCutOn = isCutEnabled(secondaryCuts.cMinPiKaMassCut) || isCutEnabled(secondaryCuts.cMaxPiKaMassCut);
    angleCutOn = isCutEnabled(secondaryCuts.cMinAngle) || isCutEnabled(secondaryCuts.cMaxAngle);
    pairAsymCutOn = isCutEnabled(secondaryCuts.cMinPairAsym) || isCutEnabled(secondaryCuts.cMaxPairAsym);

    // Consistency of the pT dependent PID configuration
    if (pionPID.cPionUsePtDepPID) {
      const auto& bins = pionPID.cPionPIDPtBins.value;
      if (bins.size() < 2 || pionPID.cPionTPCNSigmaCuts.value.size() != bins.size() - 1 ||
          pionPID.cPionTOFNSigmaCuts.value.size() != bins.size() - 1 || pionPID.cPionTOFRequired.value.size() != bins.size() - 1) {
        LOG(fatal) << "Pion pT dependent PID vectors must have (number of pT bin edges - 1) entries";
      }
    }
    if (kaonPID.cKaonUsePtDepPID) {
      const auto& bins = kaonPID.cKaonPIDPtBins.value;
      if (bins.size() < 2 || kaonPID.cKaonTPCNSigmaCuts.value.size() != bins.size() - 1 ||
          kaonPID.cKaonTOFNSigmaCuts.value.size() != bins.size() - 1 || kaonPID.cKaonTOFRequired.value.size() != bins.size() - 1) {
        LOG(fatal) << "Kaon pT dependent PID vectors must have (number of pT bin edges - 1) entries";
      }
    }
    if (cByPassTOF && (pionPID.cUseOnlyTOFTrackPi || kaonPID.cUseOnlyTOFTrackKa)) {
      LOG(warning) << "cByPassTOF skips the TOF nSigma cut, but cUseOnlyTOFTrack* still requires a TOF signal";
    }

    // Micro tracks store quantised DCA and nSigma: a cut off the grid would silently act as a different cut.
    if (doprocessResoMicroTracks || doprocessMCMicro || doprocessMEMicro) {
      auto checkDCAGrid = [](const char* name, double cut) {
        const double nearest = std::min(std::max(std::round(cut / DCAGridStep) * DCAGridStep, 0.), DCAGridMax);
        if (std::abs(cut - nearest) > GridTolerance) {
          LOG(fatal) << name << " = " << cut << " is not on the quantised DCA grid (multiples of " << DCAGridStep << " up to " << DCAGridMax << "); nearest value: " << nearest;
        }
      };
      auto checkPIDGrid = [](const char* name, double cut) {
        const double nearest = std::min(std::max(PIDGridStart + std::round((cut - PIDGridStart) / PIDGridStep) * PIDGridStep, PIDGridStart), PIDGridMax);
        if (std::abs(cut - nearest) > GridTolerance) {
          LOG(fatal) << name << " = " << cut << " is not on the quantised nSigma grid ({2.0, 2.25, ..., 3.5}); nearest value: " << nearest;
        }
      };
      if (trackCuts.cfgUsePtDepDCA) {
        LOG(info) << "Micro tracks use the producer pT dependent DCA flags (0.004 + 0.013 / pT); cDCAToPVByPt* are ignored";
        if (std::abs(trackCuts.cDCAToPVByPtP0 - 0.004f) > 1e-6f || std::abs(trackCuts.cDCAToPVByPtCoeff - 0.013f) > 1e-6f || std::abs(trackCuts.cDCAToPVByPtPower - 1.f) > 1e-6f) {
          LOG(warning) << "cDCAToPVByPt* differ from the producer defaults, but micro tracks always use the producer formula";
        }
      } else {
        if (isCutEnabled(trackCuts.cMaxDCArToPVcut)) {
          checkDCAGrid("cMaxDCArToPVcut", trackCuts.cMaxDCArToPVcut);
        }
        if (isCutEnabled(trackCuts.cMaxDCAzToPVcut)) {
          checkDCAGrid("cMaxDCAzToPVcut", trackCuts.cMaxDCAzToPVcut);
        }
      }
      if (isCutEnabled(trackCuts.cMinDCAzToPVcut)) {
        checkDCAGrid("cMinDCAzToPVcut", trackCuts.cMinDCAzToPVcut);
      }
      if (isCutEnabled(pionPID.cMaxTPCnSigmaPion) && !pionPID.cPionUsePtDepPID) {
        checkPIDGrid("cMaxTPCnSigmaPion", pionPID.cMaxTPCnSigmaPion);
      }
      if (isCutEnabled(pionPID.cMaxTOFnSigmaPion) && !pionPID.cPionUsePtDepPID) {
        checkPIDGrid("cMaxTOFnSigmaPion", pionPID.cMaxTOFnSigmaPion);
      }
      if (isCutEnabled(kaonPID.cMaxTPCnSigmaKaon) && !kaonPID.cKaonUsePtDepPID) {
        checkPIDGrid("cMaxTPCnSigmaKaon", kaonPID.cMaxTPCnSigmaKaon);
      }
      if (isCutEnabled(kaonPID.cMaxTOFnSigmaKaon) && !kaonPID.cKaonUsePtDepPID) {
        checkPIDGrid("cMaxTOFnSigmaKaon", kaonPID.cMaxTOFnSigmaKaon);
      }
      if (pionPID.cPionUsePtDepPID) {
        for (const auto cut : pionPID.cPionTPCNSigmaCuts.value) {
          if (isCutEnabled(cut)) {
            checkPIDGrid("cPionTPCNSigmaCuts", cut);
          }
        }
        for (const auto cut : pionPID.cPionTOFNSigmaCuts.value) {
          if (isCutEnabled(cut)) {
            checkPIDGrid("cPionTOFNSigmaCuts", cut);
          }
        }
      }
      if (kaonPID.cKaonUsePtDepPID) {
        for (const auto cut : kaonPID.cKaonTPCNSigmaCuts.value) {
          if (isCutEnabled(cut)) {
            checkPIDGrid("cKaonTPCNSigmaCuts", cut);
          }
        }
        for (const auto cut : kaonPID.cKaonTOFNSigmaCuts.value) {
          if (isCutEnabled(cut)) {
            checkPIDGrid("cKaonTOFNSigmaCuts", cut);
          }
        }
      }
      if (pionPID.nsigmaCutCombinedPion > 0 || kaonPID.nsigmaCutCombinedKaon > 0) {
        LOG(warning) << "nsigmaCutCombined* on micro tracks uses quantised nSigma values (approximate)";
      }
    }

    std::vector<double> centBinning = {0., 1., 5., 10., 15., 20., 25., 30., 35., 40., 45., 50., 55., 60., 65., 70., 80., 90., 100., 200.};
    AxisSpec centAxis = {centBinning, "T0M (%)"};
    AxisSpec ptAxis = {150, 0, 15, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec dcaxyAxis = {300, 0, 3, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {500, 0, 5, "DCA_{#it{z}} (cm)"};
    AxisSpec invMassAxisK892 = {1400 / cNbinsDiv, 0.6, 2.0, "Invariant Mass (GeV/#it{c}^2)"};   // K(892)0
    AxisSpec invMassAxisRho = {2000 / cNbinsDiv, 0.0, 2.0, "Invariant Mass (GeV/#it{c}^2)"};    // rho
    AxisSpec invMassAxisReso = {1600 / cNbinsDiv, 0.9f, 2.5f, "Invariant Mass (GeV/#it{c}^2)"}; // K1
    AxisSpec invMassAxisScan = {250, 0, 2.5, "Invariant Mass (GeV/#it{c}^2)"};                  // For selection
    AxisSpec pidQAAxis = {130, -6.5, 6.5};
    AxisSpec dataTypeAxis = {9, 0, 9, "Histogram types"};
    AxisSpec mcTypeAxis = {4, 0, 4, "Histogram types"};

    // THnSparse
    AxisSpec axisAnti = {BinAnti::kNAEnd, 0, BinAnti::kNAEnd, "Type of bin: Normal or Anti"};
    AxisSpec axisType = {BinType::kTYEnd, 0, BinType::kTYEnd, "Type of bin with charge and mix"};
    AxisSpec mcLabelAxis = {5, -0.5, 4.5, "MC Label"};

    // Micro-only instrumentation: category 0 includes all combinations, not just unmatched.
    auto trackFlow = histos.add<TH2>("CutFlow/tracks", "Micro tracks, once per selected collision;stage;species", HistType::kTH2D, {{static_cast<int>(kTrkNStages), -0.5, static_cast<int>(kTrkNStages) - 0.5}, {2, -0.5, 1.5}});
    const std::array<const char*, TrackStage::kTrkNStages> trackLabels{"input", "pT", "eta", "DCAxy", "DCAz", "track flags", "clusters / crossed rows", "TOF required", "PID"};
    for (size_t i = 0; i < trackLabels.size(); ++i) {
      trackFlow->GetXaxis()->SetBinLabel(i + 1, trackLabels[i]);
    }
    trackFlow->GetYaxis()->SetBinLabel(1, "pion");
    trackFlow->GetYaxis()->SetBinLabel(2, "kaon");
    auto candidateFlow = histos.add<TH2>("CutFlow/candidates", "Unordered micro triplets;stage;category", HistType::kTH2D, {{NCandidateStages, -0.5, NCandidateStages - 0.5}, {3, -0.5, 2.5}});
    const std::array<const char*, NCandidateStages> candidateLabels{"input unordered triplets", "distinct pion IDs", "pion selection (quality+PID)", "pion pair constructed", "pair pT", "secondary mass window (rho mode)", "three distinct IDs", "kaon selection (quality+PID)", "K1 rapidity", "candidate cuts", "final US", "final LS"};
    for (size_t i = 0; i < candidateLabels.size(); ++i) {
      candidateFlow->GetXaxis()->SetBinLabel(i + 1, candidateLabels[i]);
    }
    candidateFlow->GetYaxis()->SetBinLabel(1, "all");
    candidateFlow->GetYaxis()->SetBinLabel(2, "rho K");
    candidateFlow->GetYaxis()->SetBinLabel(3, "K* pi");
    if (doprocessMCMicro) {
      auto mothers = histos.add<TH1>("CutFlow/uniqueMothersPerCollision", "Final unique K1 IDs summed over reconstructed collisions (not globally deduplicated)", HistType::kTH1D, {{2, 0.5, 2.5}});
      mothers->GetXaxis()->SetBinLabel(1, "rho K");
      mothers->GetXaxis()->SetBinLabel(2, "K* pi");
    }
    if (doprocessMCTrue) {
      auto generated = histos.add<TH2>("CutFlow/generated", "K1 parent rows conditional on selected reconstructed events;stage;immediate channel", HistType::kTH2D, {{2, -0.5, 1.5}, {3, -0.5, 2.5}});
      generated->GetXaxis()->SetBinLabel(1, "all K1 parent rows");
      generated->GetXaxis()->SetBinLabel(2, "K1 rapidity window");
      generated->GetYaxis()->SetBinLabel(1, "other / unresolved");
      generated->GetYaxis()->SetBinLabel(2, "rho K");
      generated->GetYaxis()->SetBinLabel(3, "K* pi");
    }

    // DCA QA
    // Primary pion
    histos.add("QA/trkppionDCAxy", "DCAxy disstribution of primary pion candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QA/trkppionDCAz", "DCAz disstribution of primary pion candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QA/trkppionpT", "pT distribution of primary pion candidates", HistType::kTH1F, {ptAxis});
    histos.add("QA/trkppionTPCPID", "TPC PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkppionTOFPID", "TOF PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkppionTPCTOFPID", "TPC-TOF PID map of primary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    histos.add("QAcut/trkppionDCAxy", "DCAxy distribution of primary pion candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAcut/trkppionDCAz", "DCAz distribution of primary pion candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAcut/trkppionpT", "pT distribution of primary pion candidates", HistType::kTH1F, {ptAxis});
    histos.add("QAcut/trkppionTPCPID", "TPC PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkppionTOFPID", "TOF PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkppionTPCTOFPID", "TPC-TOF PID map of primary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    // Secondary pion
    histos.add("QA/trkspionDCAxy", "DCAxy distribution of secondary pion candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QA/trkspionDCAz", "DCAz distribution of secondary pion candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QA/trkspionpT", "pT distribution of secondary pion candidates", HistType::kTH1F, {ptAxis});
    histos.add("QA/trkspionTPCPID", "TPC PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkspionTOFPID", "TOF PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkspionTPCTOFPID", "TPC-TOF PID map of secondary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    histos.add("QAcut/trkspionDCAxy", "DCAxy distribution of secondary pion candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAcut/trkspionDCAz", "DCAz distribution of secondary pion candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAcut/trkspionpT", "pT distribution of secondary pion candidates", HistType::kTH1F, {ptAxis});
    histos.add("QAcut/trkspionTPCPID", "TPC PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkspionTOFPID", "TOF PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkspionTPCTOFPID", "TPC-TOF PID map of secondary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    // Kaon
    histos.add("QA/trkkaonDCAxy", "DCAxy distribution of kaon candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QA/trkkaonDCAz", "DCAz distribution of kaon candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QA/trkkaonpT", "pT distribution of kaon candidates", HistType::kTH1F, {ptAxis});
    histos.add("QA/trkkaonTPCPID", "TPC PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkkaonTOFPID", "TOF PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QA/trkkaonTPCTOFPID", "TPC-TOF PID map of kaon candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    histos.add("QAcut/trkkaonDCAxy", "DCAxy distribution of kaon candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAcut/trkkaonDCAz", "DCAz distribution of kaon candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAcut/trkkaonpT", "pT distribution of kaon candidates", HistType::kTH1F, {ptAxis});
    histos.add("QAcut/trkkaonTPCPID", "TPC PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkkaonTOFPID", "TOF PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
    histos.add("QAcut/trkkaonTPCTOFPID", "TPC-TOF PID map of kaon candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

    // K1
    histos.add("QA/K1OA", "Opening angle of K1(1270)", HistType::kTH1F, {AxisSpec{100, 0, 3.14, "Opening angle"}});
    histos.add("QA/K1PairAsym", "Pair asymmetry of K1(1270)", HistType::kTH1F, {AxisSpec{100, -1, 1, "Pair asymmetry"}});
    histos.add("QA/hInvmassK892_Rho", "Invariant mass of K(892)0 vs Rho(770)", HistType::kTH2F, {invMassAxisK892, invMassAxisRho});
    histos.add("QA/hInvmassSecon_PiKa", "Invariant mass of secondary resonance vs pion-kaon", HistType::kTH2F, {invMassAxisRho, invMassAxisK892});
    histos.add("QA/hInvmassSecon", "Invariant mass of secondary resonance", HistType::kTH1F, {invMassAxisRho});
    histos.add("QA/hpT_Secondary", "pT distribution of secondary resonance", HistType::kTH1F, {ptAxis});

    histos.add("QAcut/K1OA", "Opening angle of K1(1270)", HistType::kTH1F, {AxisSpec{100, 0, 3.14, "Opening angle"}});
    histos.add("QAcut/K1PairAsym", "Pair asymmetry of K1(1270)", HistType::kTH1F, {AxisSpec{100, -1, 1, "Pair asymmetry"}});
    histos.add("QAcut/hInvmassK892_Rho", "Invariant mass of K(892)0 vs Rho(770)", HistType::kTH2F, {invMassAxisK892, invMassAxisRho});
    histos.add("QAcut/hInvmassSecon_PiKa", "Invariant mass of secondary resonance vs pion-kaon", HistType::kTH2F, {invMassAxisRho, invMassAxisK892});
    histos.add("QAcut/hInvmassSecon", "Invariant mass of secondary resonance", HistType::kTH1F, {invMassAxisRho});
    histos.add("QAcut/hpT_Secondary", "pT distribution of secondary resonance", HistType::kTH1F, {ptAxis});

    // Invariant mass
    histos.add("hInvmass_K1", "Invariant mass of K1(1270) (US)", HistType::kTHnSparseD, {axisAnti, axisType, centAxis, ptAxis, invMassAxisReso});
    histos.add("hInvmass_K1_LS", "Invariant mass of K1(1270) (LS)", HistType::kTHnSparseD, {axisAnti, axisType, centAxis, ptAxis, invMassAxisReso});
    histos.add("hInvmass_K1_Mix", "Invariant mass of K1(1270) (ME)", HistType::kTHnSparseD, {axisAnti, axisType, centAxis, ptAxis, invMassAxisReso});
    // Mass QA (quick check)
    histos.add("k1invmass", "Invariant mass of K1(1270) (US)", HistType::kTH1F, {invMassAxisReso});
    histos.add("k1invmass_LS", "Invariant mass of K1(1270) (LS)", HistType::kTH1F, {invMassAxisReso});
    histos.add("k1invmass_Mix", "Invariant mass of K1(1270) (ME)", HistType::kTH1F, {invMassAxisReso});

    // MC
    if (doprocessMC || doprocessMCMicro) {
      AxisSpec channelAxis = {3, -0.5, 2.5, "0: non-K1, 1: rho K, 2: K* pi"};
      histos.add("MCReco/collisions", "Selected reconstructed MC collisions", HistType::kTH1D, {{1, 0, 1}});
      histos.add("MCReco/microTracks", "Input micro tracks in selected MC collisions", HistType::kTH1D, {{1, 0, 1}});
      histos.add("MCReco/channel", "All selected pi-pi-K combinations by truth channel", HistType::kTH1D, {channelAxis});
      histos.add("MCReco/mass", "Reconstructed mass by truth channel", HistType::kTH2D, {channelAxis, invMassAxisReso});
      histos.add("MCReco/pt", "Reconstructed pT by truth channel", HistType::kTH2D, {channelAxis, ptAxis});
      histos.add("MCReco/piPiMass", "pi-pi mass by truth channel", HistType::kTH2D, {channelAxis, invMassAxisRho});
      histos.add("MCReco/pi1KMass", "First pion-kaon mass by truth channel", HistType::kTH2D, {channelAxis, invMassAxisK892});
      histos.add("MCReco/pi2KMass", "Second pion-kaon mass by truth channel", HistType::kTH2D, {channelAxis, invMassAxisK892});
      histos.add("k1invmass_MC", "Invariant mass of K1(1270)", HistType::kTH1F, {invMassAxisReso});
      histos.add("k1invmass_MC_noK1", "Invariant mass of K1(1270)", HistType::kTH1F, {invMassAxisReso});

      histos.add("QAMC/trkppionDCAxy", "DCAxy distribution of primary pion candidates", HistType::kTH1F, {dcaxyAxis});
      histos.add("QAMC/trkppionDCAz", "DCAz distribution of primary pion candidates", HistType::kTH1F, {dcazAxis});
      histos.add("QAMC/trkppionpT", "pT distribution of primary pion candidates", HistType::kTH1F, {ptAxis});
      histos.add("QAMC/trkppionTPCPID", "TPC PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkppionTOFPID", "TOF PID of primary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkppionTPCTOFPID", "TPC-TOF PID map of primary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

      histos.add("QAMC/trkspionDCAxy", "DCAxy distribution of secondary pion candidates", HistType::kTH1F, {dcaxyAxis});
      histos.add("QAMC/trkspionDCAz", "DCAz distribution of secondary pion candidates", HistType::kTH1F, {dcazAxis});
      histos.add("QAMC/trkspionpT", "pT distribution of secondary pion candidates", HistType::kTH1F, {ptAxis});
      histos.add("QAMC/trkspionTPCPID", "TPC PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkspionTOFPID", "TOF PID of secondary pion candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkspionTPCTOFPID", "TPC-TOF PID map of secondary pion candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

      histos.add("QAMC/trkkaonDCAxy", "DCAxy distribution of kaon candidates", HistType::kTH1F, {dcaxyAxis});
      histos.add("QAMC/trkkaonDCAz", "DCAz distribution of kaon candidates", HistType::kTH1F, {dcazAxis});
      histos.add("QAMC/trkkaonpT", "pT distribution of kaon candidates", HistType::kTH1F, {ptAxis});
      histos.add("QAMC/trkkaonTPCPID", "TPC PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkkaonTOFPID", "TOF PID of kaon candidates", HistType::kTH2F, {ptAxis, pidQAAxis});
      histos.add("QAMC/trkkaonTPCTOFPID", "TPC-TOF PID map of kaon candidates", HistType::kTH2F, {pidQAAxis, pidQAAxis});

      histos.add("QAMC/K1OA", "Opening angle of K1(1270)", HistType::kTH1F, {AxisSpec{100, 0, 3.14, "Opening angle"}});
      histos.add("QAMC/K1PairAsym", "Pair asymmetry of K1(1270)", HistType::kTH1F, {AxisSpec{100, -1, 1, "Pair asymmetry"}});
      histos.add("QAMC/hInvmassK892_Rho", "Invariant mass of K(892)0 vs Rho(770)", HistType::kTH2F, {invMassAxisK892, invMassAxisRho});
      histos.add("QAMC/hInvmassSecon_PiKa", "Invariant mass of secondary resonance vs pion-kaon", HistType::kTH2F, {invMassAxisRho, invMassAxisK892});
      histos.add("QAMC/hInvmassSecon", "Invariant mass of secondary resonance", HistType::kTH1F, {invMassAxisRho});
      histos.add("QAMC/hpT_Secondary", "pT distribution of secondary resonance", HistType::kTH1F, {ptAxis});
    } // doprocessMC
    if (doprocessMCTrue) {
      AxisSpec channelAxis = {3, -0.5, 2.5, "0: other/unresolved, 1: rho K, 2: K* pi"};
      histos.add("MCGen/chargeChannel", "K1 parents in selected reconstructed events, inside the K1 rapidity window", HistType::kTH2D, {{2, -1.5, 1.5, "K1 charge"}, channelAxis});
      histos.add("MCGen/ptChannel", "Generated K1 pT by immediate decay channel", HistType::kTH2D, {channelAxis, ptAxis});
    }
    // Print output histograms statistics
    LOG(info) << "Size of the histograms in K1 Analysis Task";
    histos.print();
  } // init

  // Resolve the PID cut of one species at a given pT; false if the pT is outside all pT-dependent bins.
  template <Species S>
  bool getPIDCut(float pt, PIDCut& cut)
  {
    if constexpr (S == Species::Pion) {
      cut.tpcMax = pionPID.cMaxTPCnSigmaPion;
      cut.tofMax = pionPID.cMaxTOFnSigmaPion;
      cut.combined = pionPID.nsigmaCutCombinedPion;
      cut.tofRequired = false;
      if (pionPID.cPionUsePtDepPID) {
        const int ptBin = getPtBinIndex(pt, pionPID.cPionPIDPtBins.value);
        if (ptBin < 0) {
          return false;
        }
        const auto bin = static_cast<std::size_t>(ptBin);
        cut.tpcMax = pionPID.cPionTPCNSigmaCuts.value[bin];
        cut.tofMax = pionPID.cPionTOFNSigmaCuts.value[bin];
        cut.tofRequired = pionPID.cPionTOFRequired.value[bin] != 0;
      }
    } else {
      cut.tpcMax = kaonPID.cMaxTPCnSigmaKaon;
      cut.tofMax = kaonPID.cMaxTOFnSigmaKaon;
      cut.combined = kaonPID.nsigmaCutCombinedKaon;
      cut.tofRequired = false;
      if (kaonPID.cKaonUsePtDepPID) {
        const int ptBin = getPtBinIndex(pt, kaonPID.cKaonPIDPtBins.value);
        if (ptBin < 0) {
          return false;
        }
        const auto bin = static_cast<std::size_t>(ptBin);
        cut.tpcMax = kaonPID.cKaonTPCNSigmaCuts.value[bin];
        cut.tofMax = kaonPID.cKaonTOFNSigmaCuts.value[bin];
        cut.tofRequired = kaonPID.cKaonTOFRequired.value[bin] != 0;
      }
    }
    return true;
  }

  // Track quality selection shared by pion and kaon. Returns the last stage that was passed.
  // Full tracks store exact values; micro tracks store quantised DCA (see LFResonanceTables.h).
  template <bool IsResoMicrotrack, typename TrackType>
  int trackQualityStage(const TrackType& track)
  {
    const double pt = track.pt();
    const double dcaXY = track.dcaXY();
    const double dcaZ = track.dcaZ();
    // Invalid micro DCA codes decode to NaN
    if (!std::isfinite(pt) || !std::isfinite(track.eta()) || !std::isfinite(dcaXY) || !std::isfinite(dcaZ)) {
      return kTrkInput;
    }
    if (std::abs(pt) < trackCuts.cMinPtcut) {
      return kTrkInput;
    }
    if (isCutEnabled(trackCuts.cMaxEtacut) && !(std::abs(track.eta()) < trackCuts.cMaxEtacut)) {
      return kTrkPt;
    }

    if (trackCuts.cfgUsePtDepDCA) {
      if constexpr (IsResoMicrotrack) {
        if (!track.passedPtDependentDCAxy()) {
          return kTrkEta;
        }
        if (!track.passedPtDependentDCAz()) {
          return kTrkDCAxy;
        }
      } else {
        const double dcaPtCut = trackCuts.cDCAToPVByPtP0 + trackCuts.cDCAToPVByPtCoeff * std::pow(pt, -static_cast<double>(trackCuts.cDCAToPVByPtPower));
        if (!(std::abs(dcaXY) < dcaPtCut)) {
          return kTrkEta;
        }
        if (!(std::abs(dcaZ) < dcaPtCut)) {
          return kTrkDCAxy;
        }
      }
    } else {
      if (isCutEnabled(trackCuts.cMaxDCArToPVcut)) {
        if constexpr (IsResoMicrotrack) {
          if (!passesBinnedMax(dcaXY, trackCuts.cMaxDCArToPVcut)) {
            return kTrkEta;
          }
        } else {
          if (!(std::abs(dcaXY) <= trackCuts.cMaxDCArToPVcut)) {
            return kTrkEta;
          }
        }
      }
      if (isCutEnabled(trackCuts.cMaxDCAzToPVcut)) {
        if constexpr (IsResoMicrotrack) {
          if (!passesBinnedMax(dcaZ, trackCuts.cMaxDCAzToPVcut)) {
            return kTrkDCAxy;
          }
        } else {
          if (!(std::abs(dcaZ) <= trackCuts.cMaxDCAzToPVcut)) {
            return kTrkDCAxy;
          }
        }
      }
    }
    if (isCutEnabled(trackCuts.cMinDCAzToPVcut)) {
      if constexpr (IsResoMicrotrack) {
        if (!passesBinnedMin(dcaZ, trackCuts.cMinDCAzToPVcut)) {
          return kTrkDCAxy;
        }
      } else {
        if (!(std::abs(dcaZ) >= trackCuts.cMinDCAzToPVcut)) {
          return kTrkDCAxy;
        }
      }
    }

    // Track flags
    if ((trackCuts.cfgPrimaryTrack && !track.isPrimaryTrack()) ||
        (trackCuts.cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA()) ||
        (trackCuts.cfgGlobalTrack && !track.isGlobalTrack()) ||
        (trackCuts.cfgPVContributor && !track.isPVContributor()) ||
        (trackCuts.cfgUseITSRefit && !track.passedITSRefit()) ||
        (trackCuts.cfgUseTPCRefit && !track.passedTPCRefit())) {
      return kTrkDCAz;
    }

    // Clusters: found clusters exist only in ResoTracks, ITS clusters only in ResoMicroTracks
    if constexpr (!IsResoMicrotrack) {
      if constexpr (requires { track.tpcNClsFound(); }) {
        if (track.tpcNClsFound() < trackCuts.cfgTPCcluster) {
          return kTrkFlags;
        }
      }
    }
    if constexpr (requires { track.tpcNClsCrossedRows(); }) {
      if (track.tpcNClsCrossedRows() < trackCuts.cfgTPCCrossedRowsMin) {
        return kTrkFlags;
      }
    }
    if constexpr (IsResoMicrotrack) {
      if constexpr (requires { track.itsNCls(); }) {
        if (track.itsNCls() < trackCuts.cfgITSNClsMin) {
          return kTrkFlags;
        }
      }
    }
    return kTrkClusters;
  }

  // TOF signal requirement of the track (global, per species, or per pT bin)
  template <Species S, typename TrackType>
  bool passesTOFRequired(const TrackType& track)
  {
    bool required = trackCuts.cfgHasTOF;
    if constexpr (S == Species::Pion) {
      required = required || pionPID.cUseOnlyTOFTrackPi;
    } else {
      required = required || kaonPID.cUseOnlyTOFTrackKa;
    }
    PIDCut cut;
    // A pT outside all bins is rejected by passesPID
    if (!cByPassTOF && getPIDCut<S>(track.pt(), cut) && cut.tofRequired) {
      required = true;
    }
    return !required || track.hasTOF();
  }

  // PID selection: the same code for full and micro tracks, only the comparison is quantisation aware
  template <bool IsResoMicrotrack, Species S, typename TrackType>
  bool passesPID(const TrackType& track)
  {
    PIDCut cut;
    if (!getPIDCut<S>(track.pt(), cut)) {
      return false;
    }
    const bool hasTOF = track.hasTOF();
    double tpcNSigma = std::numeric_limits<double>::quiet_NaN();
    double tofNSigma = std::numeric_limits<double>::quiet_NaN(); // TOF value is only valid with hasTOF
    if constexpr (S == Species::Pion) {
      tpcNSigma = track.tpcNSigmaPi();
      if (hasTOF) {
        tofNSigma = track.tofNSigmaPi();
      }
    } else {
      tpcNSigma = track.tpcNSigmaKa();
      if (hasTOF) {
        tofNSigma = track.tofNSigmaKa();
      }
    }
    if (isCutEnabled(cut.tpcMax) && !passesMax<IsResoMicrotrack>(std::abs(tpcNSigma), cut.tpcMax)) {
      return false;
    }
    // Missing TOF is handled by passesTOFRequired; here the TPC alone decides
    if (cByPassTOF || !hasTOF) {
      return true;
    }
    bool tofPassed = !isCutEnabled(cut.tofMax) || passesMax<IsResoMicrotrack>(std::abs(tofNSigma), cut.tofMax);
    if (!tofPassed && cut.combined > 0 && tpcNSigma * tpcNSigma + tofNSigma * tofNSigma < cut.combined * cut.combined) {
      tofPassed = true;
    }
    return tofPassed;
  }

  // Full selection stage of a track (quality, TOF requirement, PID)
  template <bool IsResoMicrotrack, Species S, typename TrackType>
  int trackSelectionStage(const TrackType& track)
  {
    const int qualityStage = trackQualityStage<IsResoMicrotrack>(track);
    if (qualityStage < kTrkClusters) {
      return qualityStage;
    }
    if (!passesTOFRequired<S>(track)) {
      return kTrkClusters;
    }
    if (!passesPID<IsResoMicrotrack, S>(track)) {
      return kTrkTOFRequired;
    }
    return kTrkPID;
  }

  template <bool IsResoMicrotrack, Species S, typename TrackType>
  bool selectTrack(const TrackType& track)
  {
    return trackSelectionStage<IsResoMicrotrack, S>(track) == kTrkPID;
  }

  // Selection cache of one track slice. The row number of a grouped slice is global, hence the offset.
  template <typename TrackType>
  std::size_t getCacheIndex(const TrackType& track, int64_t firstIndex, std::size_t size)
  {
    const int64_t index = static_cast<int64_t>(track.index()) - firstIndex;
    if (index < 0 || index >= static_cast<int64_t>(size)) {
      LOG(fatal) << "Track index " << track.index() << " is outside the selection cache [" << firstIndex << ", " << firstIndex + static_cast<int64_t>(size) << ")";
    }
    return static_cast<std::size_t>(index);
  }

  template <bool IsResoMicrotrack, Species S, bool FillCutFlow, typename TracksType>
  std::vector<uint8_t> buildSelectionCache(const TracksType& tracks, int64_t firstIndex)
  {
    std::vector<uint8_t> selected(tracks.size(), 0);
    for (const auto& track : tracks) {
      const int stage = trackSelectionStage<IsResoMicrotrack, S>(track);
      selected[getCacheIndex(track, firstIndex, selected.size())] = (stage == kTrkPID) ? 1 : 0;
      if constexpr (FillCutFlow) {
        for (int i = 0; i <= stage; ++i) {
          histos.fill(HIST("CutFlow/tracks"), i, static_cast<int>(S));
        }
      }
    }
    return selected;
  }

  enum class K1TruthChannel {
    None = 0,
    RhoK = 1,
    KStarPi = 2
  };

  template <typename Track>
  bool hasSibling(const Track& directDaughter, int resonanceId)
  {
    if (resonanceId < 0) {
      return false;
    }
    const auto siblings = directDaughter.siblingIds();
    return siblings[0] == resonanceId || siblings[1] == resonanceId;
  }

  template <typename Track, typename Kaon>
  bool matchesKStarPi(const Track& resonancePion, const Track& directPion, const Kaon& kaon)
  {
    const int charge = kaon.pdgCode() > 0 ? 1 : -1;
    if (resonancePion.motherId() != kaon.motherId() ||
        resonancePion.motherId() == directPion.motherId()) {
      return false;
    }
    if (resonancePion.motherPDG() != charge * kK0Star892 || kaon.motherPDG() != charge * kK0Star892) {
      return false;
    }
    if (resonancePion.pdgCode() != -charge * kPiPlus || directPion.pdgCode() != charge * kPiPlus ||
        directPion.motherPDG() != charge * Pdg::kK1_1270Plus) {
      return false;
    }
    return hasSibling(directPion, kaon.motherId());
  }

  template <typename Track, typename Kaon>
  K1TruthChannel classifyK1Truth(const Track& pion1, const Track& pion2, const Kaon& kaon)
  {
    if (std::abs(pion1.pdgCode()) != kPiPlus || std::abs(pion2.pdgCode()) != kPiPlus ||
        std::abs(kaon.pdgCode()) != kKPlus) {
      return K1TruthChannel::None;
    }
    if (pion1.motherId() < 0 || pion2.motherId() < 0 || kaon.motherId() < 0) {
      return K1TruthChannel::None;
    }
    const int charge = kaon.pdgCode() > 0 ? 1 : -1;
    const bool rhoPions = pion1.motherId() == pion2.motherId() &&
                          pion1.motherPDG() == kRho770_0 && pion2.motherPDG() == kRho770_0 &&
                          pion1.pdgCode() == -pion2.pdgCode();
    if (rhoPions && kaon.motherPDG() == charge * Pdg::kK1_1270Plus &&
        kaon.motherId() != pion1.motherId() && hasSibling(kaon, pion1.motherId())) {
      return K1TruthChannel::RhoK;
    }
    if (matchesKStarPi(pion1, pion2, kaon) || matchesKStarPi(pion2, pion1, kaon)) {
      return K1TruthChannel::KStarPi;
    }
    return K1TruthChannel::None;
  }

  // Track QA of a pion; isPrimary selects the trkppion (first) or trkspion (second) histograms
  template <QAFolder Folder, typename TrackType>
  void fillPionQA(const TrackType& track, bool isPrimary)
  {
    const bool hasTOF = track.hasTOF();
    if (isPrimary) {
      if constexpr (Folder == QAFolder::Before) {
        histos.fill(HIST("QA/trkppionTPCPID"), track.pt(), track.tpcNSigmaPi());
        if (hasTOF) {
          histos.fill(HIST("QA/trkppionTOFPID"), track.pt(), track.tofNSigmaPi());
          histos.fill(HIST("QA/trkppionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        }
        histos.fill(HIST("QA/trkppionpT"), track.pt());
        histos.fill(HIST("QA/trkppionDCAxy"), track.dcaXY());
        histos.fill(HIST("QA/trkppionDCAz"), track.dcaZ());
      } else if constexpr (Folder == QAFolder::After) {
        histos.fill(HIST("QAcut/trkppionTPCPID"), track.pt(), track.tpcNSigmaPi());
        if (hasTOF) {
          histos.fill(HIST("QAcut/trkppionTOFPID"), track.pt(), track.tofNSigmaPi());
          histos.fill(HIST("QAcut/trkppionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        }
        histos.fill(HIST("QAcut/trkppionpT"), track.pt());
        histos.fill(HIST("QAcut/trkppionDCAxy"), track.dcaXY());
        histos.fill(HIST("QAcut/trkppionDCAz"), track.dcaZ());
      } else {
        histos.fill(HIST("QAMC/trkppionTPCPID"), track.pt(), track.tpcNSigmaPi());
        if (hasTOF) {
          histos.fill(HIST("QAMC/trkppionTOFPID"), track.pt(), track.tofNSigmaPi());
          histos.fill(HIST("QAMC/trkppionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        }
        histos.fill(HIST("QAMC/trkppionpT"), track.pt());
        histos.fill(HIST("QAMC/trkppionDCAxy"), track.dcaXY());
        histos.fill(HIST("QAMC/trkppionDCAz"), track.dcaZ());
      }
    } else {
      if constexpr (Folder == QAFolder::Before) {
        histos.fill(HIST("QA/trkspionTPCPID"), track.pt(), track.tpcNSigmaPi());
        if (hasTOF) {
          histos.fill(HIST("QA/trkspionTOFPID"), track.pt(), track.tofNSigmaPi());
          histos.fill(HIST("QA/trkspionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        }
        histos.fill(HIST("QA/trkspionpT"), track.pt());
        histos.fill(HIST("QA/trkspionDCAxy"), track.dcaXY());
        histos.fill(HIST("QA/trkspionDCAz"), track.dcaZ());
      } else if constexpr (Folder == QAFolder::After) {
        histos.fill(HIST("QAcut/trkspionTPCPID"), track.pt(), track.tpcNSigmaPi());
        if (hasTOF) {
          histos.fill(HIST("QAcut/trkspionTOFPID"), track.pt(), track.tofNSigmaPi());
          histos.fill(HIST("QAcut/trkspionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        }
        histos.fill(HIST("QAcut/trkspionpT"), track.pt());
        histos.fill(HIST("QAcut/trkspionDCAxy"), track.dcaXY());
        histos.fill(HIST("QAcut/trkspionDCAz"), track.dcaZ());
      } else {
        histos.fill(HIST("QAMC/trkspionTPCPID"), track.pt(), track.tpcNSigmaPi());
        if (hasTOF) {
          histos.fill(HIST("QAMC/trkspionTOFPID"), track.pt(), track.tofNSigmaPi());
          histos.fill(HIST("QAMC/trkspionTPCTOFPID"), track.tpcNSigmaPi(), track.tofNSigmaPi());
        }
        histos.fill(HIST("QAMC/trkspionpT"), track.pt());
        histos.fill(HIST("QAMC/trkspionDCAxy"), track.dcaXY());
        histos.fill(HIST("QAMC/trkspionDCAz"), track.dcaZ());
      }
    }
  }

  // Track QA of the bachelor kaon
  template <QAFolder Folder, typename TrackType>
  void fillKaonQA(const TrackType& track)
  {
    const bool hasTOF = track.hasTOF();
    if constexpr (Folder == QAFolder::Before) {
      histos.fill(HIST("QA/trkkaonTPCPID"), track.pt(), track.tpcNSigmaKa());
      if (hasTOF) {
        histos.fill(HIST("QA/trkkaonTOFPID"), track.pt(), track.tofNSigmaKa());
        histos.fill(HIST("QA/trkkaonTPCTOFPID"), track.tpcNSigmaKa(), track.tofNSigmaKa());
      }
      histos.fill(HIST("QA/trkkaonpT"), track.pt());
      histos.fill(HIST("QA/trkkaonDCAxy"), track.dcaXY());
      histos.fill(HIST("QA/trkkaonDCAz"), track.dcaZ());
    } else if constexpr (Folder == QAFolder::After) {
      histos.fill(HIST("QAcut/trkkaonTPCPID"), track.pt(), track.tpcNSigmaKa());
      if (hasTOF) {
        histos.fill(HIST("QAcut/trkkaonTOFPID"), track.pt(), track.tofNSigmaKa());
        histos.fill(HIST("QAcut/trkkaonTPCTOFPID"), track.tpcNSigmaKa(), track.tofNSigmaKa());
      }
      histos.fill(HIST("QAcut/trkkaonpT"), track.pt());
      histos.fill(HIST("QAcut/trkkaonDCAxy"), track.dcaXY());
      histos.fill(HIST("QAcut/trkkaonDCAz"), track.dcaZ());
    } else {
      histos.fill(HIST("QAMC/trkkaonTPCPID"), track.pt(), track.tpcNSigmaKa());
      if (hasTOF) {
        histos.fill(HIST("QAMC/trkkaonTOFPID"), track.pt(), track.tofNSigmaKa());
        histos.fill(HIST("QAMC/trkkaonTPCTOFPID"), track.tpcNSigmaKa(), track.tofNSigmaKa());
      }
      histos.fill(HIST("QAMC/trkkaonpT"), track.pt());
      histos.fill(HIST("QAMC/trkkaonDCAxy"), track.dcaXY());
      histos.fill(HIST("QAMC/trkkaonDCAz"), track.dcaZ());
    }
  }

  template <bool IsMC, bool IsMix, bool IsResoMicrotrack, typename CollisionType, typename TracksType>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksType& dTracks2)
  {
    if (dTracks1.size() == 0 || dTracks2.size() == 0) {
      return;
    }
    // Sets are local to this reconstructed collision: IDs cannot leak across DFs.
    // Source-file/DF deduplication across split collisions belongs in the audit.
    std::array<std::set<int>, 3> matchedMothers;

    // Selection cache: every track is selected once, not once per pair x bachelor.
    // dTracks1: bachelor kaons, dTracks2: pions (different collisions in mixed events).
    constexpr bool FillCutFlow = IsResoMicrotrack && !IsMix;
    const int64_t firstKaonIndex = dTracks1.begin().index();
    const int64_t firstPionIndex = dTracks2.begin().index();
    const auto kaonSelected = buildSelectionCache<IsResoMicrotrack, Species::Kaon, FillCutFlow>(dTracks1, firstKaonIndex);
    const auto pionSelected = buildSelectionCache<IsResoMicrotrack, Species::Pion, FillCutFlow>(dTracks2, firstPionIndex);

    // Values needed only by switched-on cuts or QA are computed only then
    const bool isK892Mode = secondaryCuts.cfgModeK892orRho;
    const bool fillQA = !IsMix && additionalQAplots;
    const bool needAngle = IsMC || fillQA || angleCutOn;
    const bool needPairAsym = IsMC || fillQA || pairAsymCutOn;
    // K892 mode: the K* candidate is (trk1, K), rho mode: the rho is (trk1, trk2)
    const bool needMass13 = IsMC || fillQA || (isK892Mode ? secondaryWindowOn : anotherMassCutOn) || (isK892Mode && (needAngle || needPairAsym));
    const bool needMass23 = IsMC || fillQA || piKaMassCutOn;
    const bool rhoWindowOn = secondaryWindowOn && !isK892Mode;

    auto multiplicity = collision.cent();
    ROOT::Math::PxPyPzMVector lDecayDaughter1, lDecayDaughter2, lResonanceSecondary, lDecayDaughter_bach, lResonanceK1, lPair13, lPair23;
    // Unordered pion pairs: each (pion, pion, kaon) triplet is filled once.
    // Here trk1 is the pion with the lower index; the roles are assigned once the bachelor is known.
    for (const auto& [trk1, trk2] : combinations(CombinationsStrictlyUpperIndexPolicy(dTracks2, dTracks2))) {
      // trk1: pion, trk2: pion, bTrack: kaon
      const bool pionsSelected = pionSelected[getCacheIndex(trk1, firstPionIndex, pionSelected.size())] && pionSelected[getCacheIndex(trk2, firstPionIndex, pionSelected.size())];
      bool pairPt = false;
      bool rhoWindow = true;
      if (pionsSelected) {
        // Resonance reconstruction
        lDecayDaughter1.SetCoordinates(trk1.px(), trk1.py(), trk1.pz(), MassPionCharged);
        lDecayDaughter2.SetCoordinates(trk2.px(), trk2.py(), trk2.pz(), MassPionCharged);
        lResonanceSecondary = lDecayDaughter1 + lDecayDaughter2;
        pairPt = !(lResonanceSecondary.Pt() < secondaryCuts.cMinSecondaryPtCut);
        rhoWindow = !rhoWindowOn || isInWindow(lResonanceSecondary.M(), MassRho770, secondaryCuts.cSecondaryMasswindow);
      }
      if constexpr (FillCutFlow) {
        // Early stages count potential triplets: each pair carries N bachelor trials.
        // This preserves the pair-first reconstruction and avoids a new cubic data loop.
        // Distinct pion IDs are guaranteed by the strictly upper index policy (stage 1 is always passed).
        const int lastStage = !pionsSelected ? 1 : !pairPt ? 3 : !rhoWindow ? 4 : 5;
        for (int stage = 0; stage <= lastStage; ++stage) {
          histos.fill(HIST("CutFlow/candidates"), stage, 0, static_cast<double>(dTracks1.size()));
        }
        if constexpr (IsMC) {
          // Match before rejecting quality/pT so both channels have an upstream numerator.
          if (std::abs(trk1.pdgCode()) == kPiPlus && trk1.pdgCode() == -trk2.pdgCode()) {
            for (const auto& bachelor : dTracks1) {
              const auto channel = classifyK1Truth(trk1, trk2, bachelor);
              if (channel != K1TruthChannel::None) {
                for (int stage = 0; stage <= lastStage; ++stage) {
                  histos.fill(HIST("CutFlow/candidates"), stage, static_cast<int>(channel));
                }
              }
            }
          }
        }
      }
      if (!pionsSelected) {
        continue;
      }

      if (fillQA) {
        fillPionQA<QAFolder::Before>(trk1, true);
        fillPionQA<QAFolder::Before>(trk2, false);
      }

      if (!pairPt) {
        continue;
      }

      if (fillQA) {
        histos.fill(HIST("QA/hInvmassSecon"), lResonanceSecondary.M());
      }
      if constexpr (IsMC) {
        histos.fill(HIST("QAMC/hpT_Secondary"), lResonanceSecondary.Pt());
      }
      // Secondary mass window (rho mode): the bachelor loop is skipped for rejected pairs
      if (!rhoWindow) {
        continue;
      }

      for (const auto& bTrack : dTracks1) {
        if (bTrack.index() == trk1.index() || bTrack.index() == trk2.index()) {
          continue;
        }
        K1TruthChannel flowChannel = K1TruthChannel::None;
        if constexpr (IsMC && IsResoMicrotrack && !IsMix) {
          flowChannel = classifyK1Truth(trk1, trk2, bTrack);
        }
        auto countCandidate = [&](int stage) {
          if constexpr (FillCutFlow) {
            histos.fill(HIST("CutFlow/candidates"), stage, 0);
            if constexpr (IsMC) {
              if (flowChannel != K1TruthChannel::None) {
                histos.fill(HIST("CutFlow/candidates"), stage, static_cast<int>(flowChannel));
              }
            }
          }
        };
        countCandidate(6);
        if (!kaonSelected[getCacheIndex(bTrack, firstKaonIndex, kaonSelected.size())]) {
          continue;
        }
        countCandidate(7);

        if (fillQA) {
          fillKaonQA<QAFolder::Before>(bTrack);
        }

        // Canonical assignment of the pion roles, once the bachelor is known.
        // Unlike-sign pair: the pion with the sign opposite to the kaon is pion 1 (K*0 partner), the other is pion 2.
        // Like-sign pair (the rule is ambiguous): the pion with the lower index is pion 1.
        const bool isUnlikeSign = trk1.sign() * trk2.sign() < 0;
        const bool swapPions = isUnlikeSign && trk1.sign() == bTrack.sign();
        const auto& pion1 = swapPions ? trk2 : trk1;
        const auto& pion2 = swapPions ? trk1 : trk2;
        const auto& lPion1 = swapPions ? lDecayDaughter2 : lDecayDaughter1;
        const auto& lPion2 = swapPions ? lDecayDaughter1 : lDecayDaughter2;

        // K1 reconstruction
        lDecayDaughter_bach.SetCoordinates(bTrack.px(), bTrack.py(), bTrack.pz(), MassKaonCharged);
        lResonanceK1 = lResonanceSecondary + lDecayDaughter_bach;

        // Cuts
        if (lResonanceK1.Rapidity() > cK1MaxRap || lResonanceK1.Rapidity() < cK1MinRap) {
          continue;
        }
        countCandidate(8);

        double mass13 = 0.;
        double mass23 = 0.;
        double lK1Angle = 0.;
        double lPairAsym = 0.;
        if (needMass13) {
          lPair13 = lPion1 + lDecayDaughter_bach;
          mass13 = lPair13.M();
        }
        if (needMass23) {
          lPair23 = lPion2 + lDecayDaughter_bach;
          mass23 = lPair23.M();
        }
        // Rho mode: secondary = (trk1, trk2) against the bachelor. K892 mode: secondary = (trk1, K) against trk2.
        if (needAngle) {
          lK1Angle = isK892Mode ? ROOT::Math::VectorUtil::Angle(lPair13, lPion2) : ROOT::Math::VectorUtil::Angle(lResonanceSecondary, lDecayDaughter_bach);
        }
        if (needPairAsym) {
          lPairAsym = isK892Mode ? (lPair13.E() - lPion2.E()) / (lPair13.E() + lPion2.E())
                                 : (lResonanceSecondary.E() - lDecayDaughter_bach.E()) / (lResonanceSecondary.E() + lDecayDaughter_bach.E());
        }

        // QA histogram before the candidate cuts
        if (fillQA) {
          histos.fill(HIST("QA/K1OA"), lK1Angle);
          histos.fill(HIST("QA/K1PairAsym"), lPairAsym);
          histos.fill(HIST("QA/hInvmassK892_Rho"), mass13, lResonanceSecondary.M());
          histos.fill(HIST("QA/hInvmassSecon_PiKa"), lResonanceSecondary.M(), mass23);
          histos.fill(HIST("QA/hpT_Secondary"), lResonanceSecondary.Pt());
        }

        // Candidate cuts (each one is evaluated only if switched on)
        if (isK892Mode && secondaryWindowOn && (!isInWindow(mass13, MassK0Star892, secondaryCuts.cSecondaryMasswindow) || pion1.sign() == bTrack.sign())) {
          continue;
        }
        if (anotherMassCutOn && !isInRange(isK892Mode ? lResonanceSecondary.M() : mass13, secondaryCuts.cMinAnotherSecondaryMassCut, secondaryCuts.cMaxAnotherSecondaryMassCut)) {
          continue;
        }
        if (piKaMassCutOn && !isInRange(mass23, secondaryCuts.cMinPiKaMassCut, secondaryCuts.cMaxPiKaMassCut)) {
          continue;
        }
        if (angleCutOn && !isInRange(lK1Angle, secondaryCuts.cMinAngle, secondaryCuts.cMaxAngle)) {
          continue;
        }
        if (pairAsymCutOn && !isInRange(lPairAsym, secondaryCuts.cMinPairAsym, secondaryCuts.cMaxPairAsym)) {
          continue;
        }
        countCandidate(9);

        // QA histograms after the candidate cuts
        if (fillQA) {
          fillPionQA<QAFolder::After>(pion1, true);
          fillPionQA<QAFolder::After>(pion2, false);
          fillKaonQA<QAFolder::After>(bTrack);
          histos.fill(HIST("QAcut/K1OA"), lK1Angle);
          histos.fill(HIST("QAcut/K1PairAsym"), lPairAsym);
          histos.fill(HIST("QAcut/hInvmassK892_Rho"), mass13, lResonanceSecondary.M());
          histos.fill(HIST("QAcut/hInvmassSecon_PiKa"), lResonanceSecondary.M(), mass23);
          histos.fill(HIST("QAcut/hInvmassSecon"), lResonanceSecondary.M());
          histos.fill(HIST("QAcut/hpT_Secondary"), lResonanceSecondary.Pt());
        }

        countCandidate(isUnlikeSign ? 10 : 11);
        if constexpr (IsMC && IsResoMicrotrack && !IsMix) {
          if (flowChannel != K1TruthChannel::None) {
            const int mother = flowChannel == K1TruthChannel::RhoK ? bTrack.motherId() :
                               std::abs(pion1.motherPDG()) == Pdg::kK1_1270Plus ? pion1.motherId() : pion2.motherId();
            if (matchedMothers[static_cast<int>(flowChannel)].insert(mother).second) {
              histos.fill(HIST("CutFlow/uniqueMothersPerCollision"), static_cast<int>(flowChannel));
            }
          }
        }

        if constexpr (!IsMix) {
          unsigned int typeK1 = bTrack.sign() > 0 ? BinType::kK1P : BinType::kK1N;
          unsigned int typeNormal = BinAnti::kNormal;
          if (isUnlikeSign) {
            histos.fill(HIST("k1invmass"), lResonanceK1.M());
            histos.fill(HIST("hInvmass_K1"), typeNormal, typeK1, multiplicity, lResonanceK1.Pt(), lResonanceK1.M());
          } else {
            histos.fill(HIST("k1invmass_LS"), lResonanceK1.M());
            histos.fill(HIST("hInvmass_K1_LS"), typeNormal, typeK1, multiplicity, lResonanceK1.Pt(), lResonanceK1.M());
          }

          if constexpr (IsMC) {
            const auto channel = classifyK1Truth(pion1, pion2, bTrack);
            const int channelBin = static_cast<int>(channel);
            histos.fill(HIST("MCReco/channel"), channelBin);
            histos.fill(HIST("MCReco/mass"), channelBin, lResonanceK1.M());
            histos.fill(HIST("MCReco/pt"), channelBin, lResonanceK1.Pt());
            histos.fill(HIST("MCReco/piPiMass"), channelBin, lResonanceSecondary.M());
            histos.fill(HIST("MCReco/pi1KMass"), channelBin, mass13);
            histos.fill(HIST("MCReco/pi2KMass"), channelBin, mass23);
            if (channel != K1TruthChannel::None) {
              if (truthDebugCounts[channelBin] < cfgTruthDebug) {
                ++truthDebugCounts[channelBin];
                LOGF(info, "K1Truth channel=%d collision=%lld tracks=(%lld,%lld,%lld) pdg=(%d,%d,%d) mothers=(%d,%d,%d) motherPDG=(%d,%d,%d) siblings=((%d,%d),(%d,%d),(%d,%d))",
                     channelBin, static_cast<long long>(collision.globalIndex()),
                     static_cast<long long>(pion1.globalIndex()), static_cast<long long>(pion2.globalIndex()), static_cast<long long>(bTrack.globalIndex()),
                     pion1.pdgCode(), pion2.pdgCode(), bTrack.pdgCode(), pion1.motherId(), pion2.motherId(), bTrack.motherId(),
                     pion1.motherPDG(), pion2.motherPDG(), bTrack.motherPDG(),
                     pion1.siblingIds()[0], pion1.siblingIds()[1], pion2.siblingIds()[0], pion2.siblingIds()[1], bTrack.siblingIds()[0], bTrack.siblingIds()[1]);
              }
              typeK1 = bTrack.sign() > 0 ? BinType::kK1P_Rec : BinType::kK1N_Rec;
              histos.fill(HIST("hInvmass_K1"), typeNormal, typeK1, multiplicity, lResonanceK1.Pt(), lResonanceK1.M());
              histos.fill(HIST("k1invmass_MC"), lResonanceK1.M());
              histos.fill(HIST("QAMC/K1OA"), lK1Angle);
              histos.fill(HIST("QAMC/K1PairAsym"), lPairAsym);
              histos.fill(HIST("QAMC/hInvmassK892_Rho"), mass13, lResonanceSecondary.M());
              histos.fill(HIST("QAMC/hInvmassSecon_PiKa"), lResonanceSecondary.M(), mass23);
              histos.fill(HIST("QAMC/hInvmassSecon"), lResonanceSecondary.M());
              histos.fill(HIST("QAMC/hpT_Secondary"), lResonanceSecondary.Pt());

              // PID QA primary and secondary pion
              fillPionQA<QAFolder::MC>(pion1, true);
              fillPionQA<QAFolder::MC>(pion2, false);
              fillKaonQA<QAFolder::MC>(bTrack);
            } else {
              histos.fill(HIST("k1invmass_MC_noK1"), lResonanceK1.M());
            }
          } // IsMC
        } else {
          unsigned int typeK1 = bTrack.sign() > 0 ? BinType::kK1P_Mix : BinType::kK1N_Mix;
          unsigned int typeNormal = BinAnti::kNormal;
          histos.fill(HIST("hInvmass_K1_Mix"), typeNormal, typeK1, multiplicity, lResonanceK1.Pt(), lResonanceK1.M());
          histos.fill(HIST("k1invmass_Mix"), lResonanceK1.M());
        }
      } // bTrack
    }
  } // fillHistograms

  template <typename CollisionType>
  bool passesEventCuts(const CollisionType& collision)
  {
    return !(eventCuts.cRecoINELgt0 && !collision.isRecINELgt0());
  }

  template <typename CollisionType>
  bool passesMCEventCuts(const CollisionType& collision)
  {
    if (eventCuts.cMCINELgt0 && !collision.isINELgt0()) {
      return false;
    }
    if (eventCuts.cMCVtxIn10 && !collision.isVtxIn10()) {
      return false;
    }
    return true;
  }

  void processResoTracks(ResoCollisions::iterator const& collision,
                         ResoTracks const& resotracks)
  {
    if (!passesEventCuts(collision)) {
      return;
    }
    fillHistograms<false, false, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processResoTracks, "Process ResoTracks", false);

  void processResoMicroTracks(ResoCollisions::iterator const& collision,
                              ResoMicroTracks const& resomicrotracks)
  {
    if (!passesEventCuts(collision)) {
      return;
    }
    fillHistograms<false, false, true>(collision, resomicrotracks, resomicrotracks);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processResoMicroTracks, "Process ResoMicroTracks", true);

  void processMC(ResoMCCols::iterator const& collision,
                 ResoMCTracks const& resotracks)
  {
    if (!passesEventCuts(collision) || !passesMCEventCuts(collision)) {
      return;
    }
    histos.fill(HIST("MCReco/collisions"), 0.5);
    fillHistograms<true, false, false>(collision, resotracks, resotracks);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processMC, "Process Event for MC", false);

  void processMCMicro(ResoMCCols::iterator const& collision, ResoMCMicroTracks const& tracks)
  {
    // The modular producer already selected these reconstructed collisions.
    // Apply precisely the same reconstruction loop as the frozen data baseline.
    if (!passesEventCuts(collision) || !passesMCEventCuts(collision)) {
      return;
    }
    histos.fill(HIST("MCReco/collisions"), 0.5);
    histos.fill(HIST("MCReco/microTracks"), 0.5, tracks.size());
    fillHistograms<true, false, true>(collision, tracks, tracks);
  }
  PROCESS_SWITCH(K1AnalysisMicro, processMCMicro, "Process reconstructed MC with micro v001 tables", false);

  void processMCTrue(ResoMCCols::iterator const& collision, ResoMCParents const& resoParents)
  {
    if (!passesEventCuts(collision) || !passesMCEventCuts(collision)) {
      return;
    }
    // Parents belong to selected reconstructed events; split reco collisions
    // repeat parent sets. This is not an unconditional generated denominator.
    for (const auto& part : resoParents) {
      if (std::abs(part.pdgCode()) != Pdg::kK1_1270Plus) {
        continue;
      }
      const int charge = part.pdgCode() > 0 ? 1 : -1;
      const int daughter1 = part.daughterPDG1();
      const int daughter2 = part.daughterPDG2();
      K1TruthChannel channel = K1TruthChannel::None;
      if ((daughter1 == kRho770_0 && daughter2 == charge * kKPlus) ||
          (daughter2 == kRho770_0 && daughter1 == charge * kKPlus)) {
        channel = K1TruthChannel::RhoK;
      } else if ((daughter1 == charge * kK0Star892 && daughter2 == charge * kPiPlus) ||
                 (daughter2 == charge * kK0Star892 && daughter1 == charge * kPiPlus)) {
        channel = K1TruthChannel::KStarPi;
      }
      histos.fill(HIST("CutFlow/generated"), 0, static_cast<int>(channel));
      if (part.y() < cK1MinRap || part.y() > cK1MaxRap) {
        continue;
      }
      histos.fill(HIST("CutFlow/generated"), 1, static_cast<int>(channel));
      // Keep other/unresolved immediate decays too; never require both pairs.
      histos.fill(HIST("MCGen/chargeChannel"), charge, static_cast<int>(channel));
      histos.fill(HIST("MCGen/ptChannel"), static_cast<int>(channel), part.pt());
    }
  }
  PROCESS_SWITCH(K1AnalysisMicro, processMCTrue, "Process generated K1 in selected events with v001 parents", false);

  // Processing Event Mixing
  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processME(ResoCollisions const& collisions, ResoTracks const& resotracks)
  {
    auto tracksTuple = std::make_tuple(resotracks);
    BinningTypeVtxZT0M colBinning{{cfgVtxBins, cfgMultBins}, true};
    SameKindPair<ResoCollisions, ResoTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (!passesEventCuts(collision1) || !passesEventCuts(collision2)) {
        continue;
      }
      fillHistograms<false, true, false>(collision1, tracks1, tracks2);
    }
  };
  PROCESS_SWITCH(K1AnalysisMicro, processME, "Process EventMixing light without partition", false);

  // Processing Event Mixing -- Micro
  void processMEMicro(ResoCollisions const& collisions, ResoMicroTracks const& resomicrotracks)
  {
    auto tracksTuple = std::make_tuple(resomicrotracks);
    BinningTypeVtxZT0M colBinning{{cfgVtxBins, cfgMultBins}, true};
    SameKindPair<ResoCollisions, ResoMicroTracks, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, collisions, tracksTuple, &cache}; // -1 is the number of the bin to skip

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (!passesEventCuts(collision1) || !passesEventCuts(collision2)) {
        continue;
      }
      fillHistograms<false, true, true>(collision1, tracks1, tracks2);
    }
  };
  PROCESS_SWITCH(K1AnalysisMicro, processMEMicro, "Process EventMixing light without partition", false);
}; // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<K1AnalysisMicro>(cfgc)};
}
