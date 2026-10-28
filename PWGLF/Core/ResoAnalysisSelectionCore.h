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
/// \file ResoAnalysisSelectionCore.h
/// \brief Event, track-quality and PID selection of resonance daughters from the reduced v001 resonance tables
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///
/// The same selection code serves full ResoTracks (exact values) and ResoMicroTracks (quantised DCA and nSigma,
/// see LFResonanceTables.h); only the comparisons are quantisation aware. A cut set to DisabledCut is off.

#ifndef PWGLF_CORE_RESOANALYSISSELECTIONCORE_H_
#define PWGLF_CORE_RESOANALYSISSELECTIONCORE_H_

#include <CommonConstants/MathConstants.h>
#include <Framework/Configurable.h>
#include <Framework/Logger.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace o2::analysis::resonance
{

inline constexpr float DisabledCut = -999.f; // an optional cut with this value is off and not evaluated
inline constexpr double DCAGridStep = 0.025; // v001 micro DCA encoding, lower-inclusive bins up to DCAGridMax
inline constexpr double DCAGridMax = 0.15;
inline constexpr double PIDGridStart = 2.0; // v001 micro nSigma encoding: 0.25 bins in [2.0, 3.5]
inline constexpr double PIDGridStep = 0.25;
inline constexpr double PIDGridMax = 3.5;
inline constexpr double GridTolerance = 1e-4;
inline constexpr std::size_t MinPtBinEdges = 2;     // a pT dependent PID table needs at least one bin
inline constexpr float ProducerDCAPtP0 = 0.004f;    // resonanceModuleInitializer cfgTightDCAOffset default
inline constexpr float ProducerDCAPtCoeff = 0.013f; // resonanceModuleInitializer cfgTightDCAPtCoefficient default
inline constexpr float ProducerDCAPtPower = 1.f;    // resonanceModuleInitializer cfgTightDCAPtPower default
inline constexpr float ConfigTolerance = 1e-6f;

// Last stage passed by a track; a cut-flow histogram can be filled directly from this value.
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

// Resolved PID cut of one species at a given pT.
struct PIDCut {
  double tpcMax = 0.;
  double tofMax = 0.;
  double combined = 0.;
  bool tofRequired = false;
};

// A cut is on unless it carries the disabled value (tolerant to the float parsing of the JSON value).
inline bool isCutEnabled(float value)
{
  return value > DisabledCut + 1.f;
}

// v001 micro values are lower-inclusive bin edges: a maximum cut on the grid keeps bins below it.
inline bool passesBinnedMax(double decoded, double cut)
{
  return decoded < cut - o2::constants::math::Epsilon;
}

// Minimum cut on the grid keeps the bin starting at the cut.
inline bool passesBinnedMin(double decoded, double cut)
{
  return decoded >= cut - o2::constants::math::Epsilon;
}

template <bool IsResoMicrotrack>
bool passesMax(double value, double cut)
{
  if constexpr (IsResoMicrotrack) {
    return passesBinnedMax(value, cut);
  } else {
    return value < cut;
  }
}

inline bool isInRange(double value, double minimum, double maximum)
{
  if (isCutEnabled(minimum) && value < minimum) {
    return false;
  }
  if (isCutEnabled(maximum) && value > maximum) {
    return false;
  }
  return true;
}

inline bool isInWindow(double value, double center, double width)
{
  return std::abs(value - center) < width;
}

// Preserve pT-bin membership [low, high).
inline int getPtBinIndex(float pt, const std::vector<float>& ptBins)
{
  for (std::size_t i = 1; i < ptBins.size(); ++i) {
    if (pt >= ptBins[i - 1] && pt < ptBins[i]) {
      return static_cast<int>(i - 1);
    }
  }
  return -1;
}

// Configurable groups without prefix: the JSON keys are the plain configurable names.

/// Event selection
struct EventCuts : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<bool> cRecoINELgt0{"cRecoINELgt0", false, "Apply reconstructed INEL>0 selection"};
  o2::framework::Configurable<bool> cMCINELgt0{"cMCINELgt0", false, "Require generator INEL>0 in MC processes"};
  o2::framework::Configurable<bool> cMCVtxIn10{"cMCVtxIn10", false, "Require generator |vz| < 10 cm in MC processes"};
};

/// Track selections (common for all daughter species, -999 switches an optional cut off)
struct TrackCuts : o2::framework::ConfigurableGroup {
  o2::framework::Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
  o2::framework::Configurable<float> cMaxEtacut{"cMaxEtacut", -999.f, "Track maximum |eta| cut (-999: off)"};
  // DCAr to PV
  o2::framework::Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.1, "Track DCAr cut to PV Maximum"};
  // DCAz to PV
  o2::framework::Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 0.1, "Track DCAz cut to PV Maximum"};
  o2::framework::Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};
  o2::framework::Configurable<bool> cfgUsePtDepDCA{"cfgUsePtDepDCA", false, "Use pT dependent DCA cut instead of the fixed maximum"};
  o2::framework::Configurable<float> cDCAToPVByPtP0{"cDCAToPVByPtP0", 0.004f, "pT dependent DCA cut = P0 + coefficient / pT^power (cm)"};
  o2::framework::Configurable<float> cDCAToPVByPtCoeff{"cDCAToPVByPtCoeff", 0.013f, "Coefficient in the pT dependent DCA cut"};
  o2::framework::Configurable<float> cDCAToPVByPtPower{"cDCAToPVByPtPower", 1.f, "Power in the pT dependent DCA cut"};
  o2::framework::Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  o2::framework::Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  o2::framework::Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                      // kGoldenChi2 | kDCAxy | kDCAz
  o2::framework::Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};          // PV Contriuibutor
  o2::framework::Configurable<bool> cfgUseTPCRefit{"cfgUseTPCRefit", false, "Require TPC Refit"};
  o2::framework::Configurable<bool> cfgUseITSRefit{"cfgUseITSRefit", false, "Require ITS Refit"};
  o2::framework::Configurable<int> cfgTPCcluster{"cfgTPCcluster", 0, "Number of TPC cluster (found clusters, ResoTracks only)"};
  o2::framework::Configurable<int> cfgTPCCrossedRowsMin{"cfgTPCCrossedRowsMin", 0, "Minimum number of TPC crossed rows"};
  o2::framework::Configurable<int> cfgITSNClsMin{"cfgITSNClsMin", 0, "Minimum number of ITS clusters (ResoMicroTracks only)"};
  o2::framework::Configurable<bool> cfgHasTOF{"cfgHasTOF", false, "Require TOF"};
};

/// PID selection of one daughter species, filled by the task from its own (species-named) configurables.
/// The configurable names are used only in the configuration messages.
struct PIDCutConfig {
  std::string species; // e.g. "Pion", used in messages
  double maxTPCnSigma = DisabledCut;
  double maxTOFnSigma = DisabledCut;
  double combinedNSigma = DisabledCut; // combined TPC-TOF cut, on when > 0
  bool onlyTOFTracks = false;          // require a TOF signal
  bool usePtDependent = false;         // use the pT binned cuts below instead of the fixed maxima
  std::vector<float> ptBins;           // bin edges; the other vectors have one entry per bin
  std::vector<float> tpcNSigmaCuts;
  std::vector<float> tofNSigmaCuts;
  std::vector<int> tofRequired;
  std::string maxTPCName; // configurable names for the messages
  std::string maxTOFName;
  std::string tpcCutsName;
  std::string tofCutsName;
};

/// Event, track-quality, TOF-requirement and PID selection of resonance daughters.
/// The species index is the position of its PIDCutConfig in init().
class ResoAnalysisSelectionCore
{
 public:
  // byPassTOF skips the TOF nSigma cut and the pT binned TOF requirement; microTracks enables the grid checks.
  void init(EventCuts const& eventCuts, TrackCuts const& trackCuts, std::vector<PIDCutConfig> pidCuts, bool byPassTOF, bool microTracks)
  {
    mEventCuts = eventCuts;
    mTrackCuts = trackCuts;
    mPID = std::move(pidCuts);
    mByPassTOF = byPassTOF;
    checkConfiguration(microTracks);
  }

  template <typename CollisionType>
  bool passesEventCuts(const CollisionType& collision)
  {
    return !(mEventCuts.cRecoINELgt0 && !collision.isRecINELgt0());
  }

  template <typename CollisionType>
  bool passesMCEventCuts(const CollisionType& collision)
  {
    if (mEventCuts.cMCINELgt0 && !collision.isINELgt0()) {
      return false;
    }
    if (mEventCuts.cMCVtxIn10 && !collision.isVtxIn10()) {
      return false;
    }
    return true;
  }

  // Resolve the PID cut of one species at a given pT; false if the pT is outside all pT-dependent bins.
  bool getPIDCut(int species, float pt, PIDCut& cut)
  {
    const auto& config = mPID[species];
    cut.tpcMax = config.maxTPCnSigma;
    cut.tofMax = config.maxTOFnSigma;
    cut.combined = config.combinedNSigma;
    cut.tofRequired = false;
    if (config.usePtDependent) {
      const int ptBin = getPtBinIndex(pt, config.ptBins);
      if (ptBin < 0) {
        return false;
      }
      const auto bin = static_cast<std::size_t>(ptBin);
      cut.tpcMax = config.tpcNSigmaCuts[bin];
      cut.tofMax = config.tofNSigmaCuts[bin];
      cut.tofRequired = config.tofRequired[bin] != 0;
    }
    return true;
  }

  // Track quality selection shared by all species. Returns the last stage that was passed.
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
    if (std::abs(pt) < mTrackCuts.cMinPtcut) {
      return kTrkInput;
    }
    if (isCutEnabled(mTrackCuts.cMaxEtacut) && !(std::abs(track.eta()) < mTrackCuts.cMaxEtacut)) {
      return kTrkPt;
    }

    if (mTrackCuts.cfgUsePtDepDCA) {
      if constexpr (IsResoMicrotrack) {
        if (!track.passedPtDependentDCAxy()) {
          return kTrkEta;
        }
        if (!track.passedPtDependentDCAz()) {
          return kTrkDCAxy;
        }
      } else {
        const double dcaPtCut = mTrackCuts.cDCAToPVByPtP0 + mTrackCuts.cDCAToPVByPtCoeff * std::pow(pt, -static_cast<double>(mTrackCuts.cDCAToPVByPtPower));
        if (!(std::abs(dcaXY) < dcaPtCut)) {
          return kTrkEta;
        }
        if (!(std::abs(dcaZ) < dcaPtCut)) {
          return kTrkDCAxy;
        }
      }
    } else {
      if (isCutEnabled(mTrackCuts.cMaxDCArToPVcut)) {
        if constexpr (IsResoMicrotrack) {
          if (!passesBinnedMax(dcaXY, mTrackCuts.cMaxDCArToPVcut)) {
            return kTrkEta;
          }
        } else {
          if (!(std::abs(dcaXY) <= mTrackCuts.cMaxDCArToPVcut)) {
            return kTrkEta;
          }
        }
      }
      if (isCutEnabled(mTrackCuts.cMaxDCAzToPVcut)) {
        if constexpr (IsResoMicrotrack) {
          if (!passesBinnedMax(dcaZ, mTrackCuts.cMaxDCAzToPVcut)) {
            return kTrkDCAxy;
          }
        } else {
          if (!(std::abs(dcaZ) <= mTrackCuts.cMaxDCAzToPVcut)) {
            return kTrkDCAxy;
          }
        }
      }
    }
    if (isCutEnabled(mTrackCuts.cMinDCAzToPVcut)) {
      if constexpr (IsResoMicrotrack) {
        if (!passesBinnedMin(dcaZ, mTrackCuts.cMinDCAzToPVcut)) {
          return kTrkDCAxy;
        }
      } else {
        if (!(std::abs(dcaZ) >= mTrackCuts.cMinDCAzToPVcut)) {
          return kTrkDCAxy;
        }
      }
    }

    // Track flags
    if ((mTrackCuts.cfgPrimaryTrack && !track.isPrimaryTrack()) ||
        (mTrackCuts.cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA()) ||
        (mTrackCuts.cfgGlobalTrack && !track.isGlobalTrack()) ||
        (mTrackCuts.cfgPVContributor && !track.isPVContributor()) ||
        (mTrackCuts.cfgUseITSRefit && !track.passedITSRefit()) ||
        (mTrackCuts.cfgUseTPCRefit && !track.passedTPCRefit())) {
      return kTrkDCAz;
    }

    // Clusters: found clusters exist only in ResoTracks, ITS clusters only in ResoMicroTracks
    if constexpr (!IsResoMicrotrack) {
      if constexpr (requires { track.tpcNClsFound(); }) {
        if (track.tpcNClsFound() < mTrackCuts.cfgTPCcluster) {
          return kTrkFlags;
        }
      }
    }
    if constexpr (requires { track.tpcNClsCrossedRows(); }) {
      if (track.tpcNClsCrossedRows() < mTrackCuts.cfgTPCCrossedRowsMin) {
        return kTrkFlags;
      }
    }
    if constexpr (IsResoMicrotrack) {
      if constexpr (requires { track.itsNCls(); }) {
        if (track.itsNCls() < mTrackCuts.cfgITSNClsMin) {
          return kTrkFlags;
        }
      }
    }
    return kTrkClusters;
  }

  // TOF signal requirement of the track (global, per species, or per pT bin)
  template <typename TrackType>
  bool passesTOFRequired(int species, const TrackType& track)
  {
    bool required = mTrackCuts.cfgHasTOF || mPID[species].onlyTOFTracks;
    PIDCut cut;
    // A pT outside all bins is rejected by passesPID
    if (!mByPassTOF && getPIDCut(species, track.pt(), cut) && cut.tofRequired) {
      required = true;
    }
    return !required || track.hasTOF();
  }

  // PID selection from the nSigma values of the species; tofNSigma is used only with hasTOF.
  template <bool IsResoMicrotrack>
  bool passesPID(int species, float pt, bool hasTOF, double tpcNSigma, double tofNSigma)
  {
    PIDCut cut;
    if (!getPIDCut(species, pt, cut)) {
      return false;
    }
    if (isCutEnabled(cut.tpcMax) && !passesMax<IsResoMicrotrack>(std::abs(tpcNSigma), cut.tpcMax)) {
      return false;
    }
    // Missing TOF is handled by passesTOFRequired; here the TPC alone decides
    if (mByPassTOF || !hasTOF) {
      return true;
    }
    bool tofPassed = !isCutEnabled(cut.tofMax) || passesMax<IsResoMicrotrack>(std::abs(tofNSigma), cut.tofMax);
    if (!tofPassed && cut.combined > 0 && tpcNSigma * tpcNSigma + tofNSigma * tofNSigma < cut.combined * cut.combined) {
      tofPassed = true;
    }
    return tofPassed;
  }

 private:
  void checkConfiguration(bool microTracks)
  {
    // Consistency of the pT dependent PID configuration
    bool anyOnlyTOF = false;
    for (const auto& config : mPID) {
      if (config.usePtDependent) {
        const auto& bins = config.ptBins;
        if (bins.size() < MinPtBinEdges || config.tpcNSigmaCuts.size() != bins.size() - 1 ||
            config.tofNSigmaCuts.size() != bins.size() - 1 || config.tofRequired.size() != bins.size() - 1) {
          LOG(fatal) << config.species << " pT dependent PID vectors must have (number of pT bin edges - 1) entries";
        }
      }
      anyOnlyTOF = anyOnlyTOF || config.onlyTOFTracks;
    }
    if (mByPassTOF && anyOnlyTOF) {
      LOG(warning) << "cByPassTOF skips the TOF nSigma selection, but cUseOnlyTOFTrack* still requires a TOF signal";
    }

    // Micro tracks store quantised DCA and nSigma: a cut off the grid would silently act as a different cut.
    if (!microTracks) {
      return;
    }
    auto checkDCAGrid = [](const char* name, double cut) {
      const double nearest = std::min(std::max(std::round(cut / DCAGridStep) * DCAGridStep, 0.), DCAGridMax);
      if (std::abs(cut - nearest) > GridTolerance) {
        LOG(fatal) << name << " = " << cut << " is not on the quantised DCA grid (multiples of " << DCAGridStep << " up to " << DCAGridMax << "); nearest value: " << nearest;
      }
    };
    auto checkPIDGrid = [](const std::string& name, double cut) {
      const double nearest = std::min(std::max(PIDGridStart + std::round((cut - PIDGridStart) / PIDGridStep) * PIDGridStep, PIDGridStart), PIDGridMax);
      if (std::abs(cut - nearest) > GridTolerance) {
        LOG(fatal) << name << " = " << cut << " is not on the quantised nSigma grid ({2.0, 2.25, ..., 3.5}); nearest value: " << nearest;
      }
    };
    if (mTrackCuts.cfgUsePtDepDCA) {
      LOG(info) << "Micro tracks use the producer pT dependent DCA flags (0.004 + 0.013 / pT); cDCAToPVByPt* are ignored";
      if (std::abs(mTrackCuts.cDCAToPVByPtP0 - ProducerDCAPtP0) > ConfigTolerance || std::abs(mTrackCuts.cDCAToPVByPtCoeff - ProducerDCAPtCoeff) > ConfigTolerance || std::abs(mTrackCuts.cDCAToPVByPtPower - ProducerDCAPtPower) > ConfigTolerance) {
        LOG(warning) << "cDCAToPVByPt* differ from the producer defaults, but micro tracks always use the producer formula";
      }
    } else {
      if (isCutEnabled(mTrackCuts.cMaxDCArToPVcut)) {
        checkDCAGrid("cMaxDCArToPVcut", mTrackCuts.cMaxDCArToPVcut);
      }
      if (isCutEnabled(mTrackCuts.cMaxDCAzToPVcut)) {
        checkDCAGrid("cMaxDCAzToPVcut", mTrackCuts.cMaxDCAzToPVcut);
      }
    }
    if (isCutEnabled(mTrackCuts.cMinDCAzToPVcut)) {
      checkDCAGrid("cMinDCAzToPVcut", mTrackCuts.cMinDCAzToPVcut);
    }
    bool anyCombined = false;
    for (const auto& config : mPID) {
      if (isCutEnabled(config.maxTPCnSigma) && !config.usePtDependent) {
        checkPIDGrid(config.maxTPCName, config.maxTPCnSigma);
      }
      if (isCutEnabled(config.maxTOFnSigma) && !config.usePtDependent) {
        checkPIDGrid(config.maxTOFName, config.maxTOFnSigma);
      }
      anyCombined = anyCombined || config.combinedNSigma > 0;
    }
    for (const auto& config : mPID) {
      if (!config.usePtDependent) {
        continue;
      }
      for (const auto& cut : config.tpcNSigmaCuts) {
        if (isCutEnabled(cut)) {
          checkPIDGrid(config.tpcCutsName, cut);
        }
      }
      for (const auto& cut : config.tofNSigmaCuts) {
        if (isCutEnabled(cut)) {
          checkPIDGrid(config.tofCutsName, cut);
        }
      }
    }
    if (anyCombined) {
      LOG(warning) << "nsigmaCutCombined* on micro tracks uses quantised nSigma values (approximate)";
    }
  }

  EventCuts mEventCuts;
  TrackCuts mTrackCuts;
  std::vector<PIDCutConfig> mPID;
  bool mByPassTOF = false;
};

} // namespace o2::analysis::resonance

#endif // PWGLF_CORE_RESOANALYSISSELECTIONCORE_H_
