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

/// \file tpcSkimsTableCreator.cxx
/// \brief Task to produce table with clean selections for TPC PID calibration
///
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>
/// \author Christian Sonnabend <christian.sonnabend@cern.ch>
/// \author Jeremy Wilkinson <jeremy.wilkinson@cern.ch>
/// \author Oleksii Lubynets <oleksii.lubynets@cern.ch>

#include "tpcSkimsTableCreator.h"

#include "utilsTpcSkimsTableCreator.h"

#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/OccupancyTables.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTPCBase.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <TH2.h>
#include <TObject.h>
#include <TRandom3.h>

#include <array>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <numeric>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::dataformats;
using namespace o2::dpg_tpcskimstablecreator;

using V0sWithID = soa::Join<aod::V0Datas, aod::V0MapID, aod::V0TOFNSigmas>;
using CascsWithID = soa::Join<aod::CascDatas, aod::CascMapID, aod::CascTOFNSigmas>;

template <typename T>
concept V0OrCasc = std::same_as<T, V0sWithID::iterator> || std::same_as<T, CascsWithID::iterator>;

struct CfgCommon : ConfigurableGroup {
  /// Configurables general
  Configurable<float> nClNorm{"nClNorm", 152., "Number of cluster normalization. Run 2: 159, Run 3 152"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<int> trackSelection{"trackSelection", 0, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  /// Configurables adaptive downsampling
  Configurable<bool> applyAdaptiveDownsampling{"applyAdaptiveDownsampling", false, "Apply downsampling driven by expected momentum distribution"};
  Configurable<std::string> ccdbPathForAdaptDwnsplFile{"ccdbPathForAdaptDwnsplFile", "", "CCDB path for file with expected-momentum histograms"};
  Configurable<std::string> adaptDwnsplFileName{"adaptDwnsplFileName", "AnalysisResults.root", "Name of the file with expected-momentum histograms"};
  Configurable<float> adaptDwnsplHistoScaleFactor{"adaptDwnsplHistoScaleFactor", 1.f, "Scaling factor for the histograms with expected-momentum"};
  Configurable<float> adaptDwnsplTargetYield{"adaptDwnsplTargetYield", 5000000.f, "Target yield of each particle species per log10(p) unity"};
  Configurable<bool> saveReserveQaHisto{"saveReserveQaHisto", true, "Flag to save the DF-wise ratio of output table size to that of input table"};
  // Configurables for run condtion table
  Configurable<std::string> rctLabel{"rctLabel", "CBT_hadronPID", "select 1 [CBT, CBT_hadronPID, CBT_muon_glo] see O2Physics/Common/CCDB/RCTSelectionFlags.h"};
  Configurable<bool> checkZdc{"checkZdc", false, "set ZDC flag for PbPb"};
  Configurable<bool> treatLimitedAcceptanceAsBad{"treatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
  Configurable<bool> requireGoodRct{"requireGoodRct", false, "require good detector flag in run condtion table"};
  // Configurable for the path of CCDB General Run Parameters LHC Interface information
  Configurable<std::string> ccdbPathGrpLhcIf{"ccdbPathGrpLhcIf", "GLO/Config/GRPLHCIF", "Path on the CCDB for the GRPLHCIF object"};
  // Configurables flagging saving trees and/or HistogramSpec
  Configurable<bool> saveTrees{"saveTrees", true, "Flag to save trees"};
  Configurable<bool> saveMomentumHistograms{"saveMomentumHistograms", false, "Flag to save momentum histograms"};
  ConfigurableAxis axisMomentumHistograms{"axisMomentumHistograms", {200, std::log10(0.1), std::log10(20.)}, "log_{10}(#it{p} / (GeV/#it{c}))"};
};

struct TreeWriterTpcV0 {

  Produces<o2::aod::SkimmedTPCV0Tree> rowTPCTree;
  Produces<o2::aod::SkimmedTPCV0TreeWithTrkQA> rowTPCTreeWithTrkQA;

  CfgCommon cfgsCommon;

  /// Configurables general
  Configurable<float> nSigmaTofDauTrackPi{"nSigmaTofDauTrackPi", 999.f, "n-sigma TOF cut on the pion daughter tracks"};
  Configurable<float> nSigmaTofDauTrackPr{"nSigmaTofDauTrackPr", 999.f, "n-sigma TOF cut on the proton daughter tracks"};
  Configurable<float> nSigmaTofDauTrackKa{"nSigmaTofDauTrackKa", 999.f, "n-sigma TOF cut on the kaon daughter tracks"};
  Configurable<bool> rejectNoTofDauTrackPi{"rejectNoTofDauTrackPi", false, "reject not matched to TOF pion daughter tracks"};
  Configurable<bool> rejectNoTofDauTrackPr{"rejectNoTofDauTrackPr", false, "reject not matched to TOF proton daughter tracks"};
  Configurable<bool> rejectNoTofDauTrackKa{"rejectNoTofDauTrackKa", false, "reject not matched to TOF kaon daughter tracks"};

  /// Configurables adaptive downsampling
  Configurable<std::string> adaptDwnsplHistoName{"adaptDwnsplHistoName", "tree-writer-tpc-v0/hParticleVsMomentum", "Name of the histogram with expected-momentum of particles"};
  // Configurables for output tables reservation size
  Configurable<float> reserveV0Ratio{"reserveV0Ratio", 0.05, "Ratio of how many tracks from V0s are expected in the output table to the input V0 table size"};
  Configurable<float> reserveCascRatio{"reserveCascRatio", 0.0025, "Ratio of how many tracks from cascades are expected in the output table to the input Cascade table size"};

  HistogramRegistry registry{"registry", {}};

  // an arbitrary value of N sigma TOF assigned by TOF task to tracks which are not matched to TOF hits
  constexpr static float NSigmaTofUnmatched{o2::aod::v0data::kNoTOFValue};
  const float nSigmaTofUnmatchedEqualityTolerance{std::fabs(NSigmaTofUnmatched) / 1e4f}; // NOLINT(cppcoreguidelines-avoid-const-or-ref-data-members): should be constexpr, but constexpr-ness of std::fabs() is not guaranteed
  // an arbitrary value of "N sigma TOF" assigned to electorns (for uniformity reasons)
  constexpr static float NSigmaTofElectorn{1000.f};

  enum { // Reconstructed V0 and cascade
    MotherUndef = -1,
    MotherGamma,
    MotherK0S,
    MotherLambda,
    MotherAntiLambda,
    MotherOmega,
    MotherAntiOmega
  };

  enum {
    DaughterUndef = -1,
    DaughterElectron,
    DaughterPion,
    DaughterKaon,
    DaughterProton
  };

  struct V0Daughter {
    double tpcNSigma{UndefValueDouble};
    double tofNSigma{UndefValueDouble};
    double itsNSigma{UndefValueDouble};
    double tpcExpSignal{UndefValueDouble};
    o2::track::PID::ID id{0};
    double nSigmaTofDauTrack{UndefValueDouble};
    bool rejectNoTofDauTrack{false};
  };

  struct V0Mother {
    int posDaughterId{UndefValueInt};
    int negDaughterId{UndefValueInt};
  };

  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  ctpRateFetcher mRateFetcher;

  o2::aod::rctsel::RCTFlagsChecker rctChecker;

  TRandom3* fRndm = new TRandom3();

  TH2* mHistMomentum{};

  double mTargetYieldPerBin{};

  using Trks = soa::Join<aod::Tracks, aod::V0Bits, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TrackSelection>;
  using TrksWithDEdxCorrection = soa::Join<Trks, aod::DEdxsCorrected>;
  using Colls = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;
  using MyBCTable = soa::Join<aod::BCsWithTimestamps, aod::BCTFinfoTable>;
  using TrksTmo = soa::Join<Trks, aod::TrackToTmo>;
  using TrksTmoWithDEdxCorrection = soa::Join<Trks, aod::DEdxsCorrected, aod::TrackToTmo>;

  Preslice<V0sWithID> perCollisionV0s = aod::v0data::collisionId;
  Preslice<CascsWithID> perCollisionCascs = aod::cascdata::collisionId;

  void init(o2::framework::InitContext&)
  {
    const std::array<bool, 7> doprocess{doprocessStandard, doprocessStandardWithCorrecteddEdx, doprocessWithdEdxTrQA, doprocessWithdEdxTrQAWithCorrecteddEdx, doprocessWithTrQA, doprocessWithTrQAWithCorrecteddEdx, doprocessDummy};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function should be enabled");
    }

    if (!(cfgsCommon.saveTrees || cfgsCommon.saveMomentumHistograms)) {
      LOGP(fatal, "At least one out of trees or momentum histograms should be saved");
    }

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);

    rctChecker.init(cfgsCommon.rctLabel, cfgsCommon.checkZdc, cfgsCommon.treatLimitedAcceptanceAsBad);

    if (cfgsCommon.saveReserveQaHisto) {
      registry.add("hV0OutputRatio", "V0 out/in ratio;V0 out/in ratio;Entries", {HistType::kTH1F, {{100, 0, reserveV0Ratio}}});
      registry.add("hCascOutputRatio", "Casc out/in ratio;Casc out/in ratio;Entries", {HistType::kTH1F, {{100, 0, reserveCascRatio}}});
    }

    if (cfgsCommon.saveMomentumHistograms) {
      registry.add("hParticleVsMomentum", "hParticleVsMomentum;log_{10}(#it{p} / (GeV/#it{c}));PID id", {HistType::kTH2D, {{cfgsCommon.axisMomentumHistograms}, {2 * o2::track::PID::NIDs, -0.5, 2 * o2::track::PID::NIDs - 0.5}}});
      for (int iParticle = 0; iParticle < o2::track::PID::NIDs; ++iParticle) {
        registry.get<TH2>(HIST("hParticleVsMomentum"))->GetYaxis()->SetBinLabel(iParticle + 1, (static_cast<std::string>(ParticleNames.at(iParticle)) + "^{+}").c_str());
        registry.get<TH2>(HIST("hParticleVsMomentum"))->GetYaxis()->SetBinLabel(iParticle + 1 + o2::track::PID::NIDs, (static_cast<std::string>(ParticleNames.at(iParticle)) + "^{-}").c_str());
      }
    }

    if (cfgsCommon.applyAdaptiveDownsampling && cfgsCommon.saveTrees) {
      initAdaptDwnsplHistograms(cfgsCommon.ccdbPathForAdaptDwnsplFile, cfgsCommon.adaptDwnsplFileName, adaptDwnsplHistoName, mHistMomentum);
      if (cfgsCommon.adaptDwnsplHistoScaleFactor != 1.f) {
        mHistMomentum->Scale(cfgsCommon.adaptDwnsplHistoScaleFactor);
      }
      mTargetYieldPerBin = cfgsCommon.adaptDwnsplTargetYield * mHistMomentum->GetXaxis()->GetBinWidth(1);
    }
  }

  template <bool IsCorrectedDeDx, typename V0Casc, typename T>
  V0Daughter createV0Daughter(const V0Casc& v0Casc, const T& track, const int motherId, const int daughterId, const bool isPositive = true)
  {
    switch (daughterId) {
      case DaughterElectron:
        return V0Daughter{track.tpcNSigmaEl(), getStrangenessTofNSigma(v0Casc, motherId, daughterId, isPositive), track.itsNSigmaEl(), track.tpcExpSignalEl(tpcSignalGeneric<IsCorrectedDeDx>(track)), PidElectron, NSigmaTofElectorn + 1.f, false};
      case DaughterPion:
        return V0Daughter{track.tpcNSigmaPi(), getStrangenessTofNSigma(v0Casc, motherId, daughterId, isPositive), track.itsNSigmaPi(), track.tpcExpSignalPi(tpcSignalGeneric<IsCorrectedDeDx>(track)), PidPion, nSigmaTofDauTrackPi, rejectNoTofDauTrackPi};
      case DaughterProton:
        return V0Daughter{track.tpcNSigmaPr(), getStrangenessTofNSigma(v0Casc, motherId, daughterId, isPositive), track.itsNSigmaPr(), track.tpcExpSignalPr(tpcSignalGeneric<IsCorrectedDeDx>(track)), PidProton, nSigmaTofDauTrackPr, rejectNoTofDauTrackPr};
      case DaughterKaon:
        return V0Daughter{track.tpcNSigmaKa(), getStrangenessTofNSigma(v0Casc, motherId, daughterId, isPositive), track.itsNSigmaKa(), track.tpcExpSignalKa(tpcSignalGeneric<IsCorrectedDeDx>(track)), PidKaon, nSigmaTofDauTrackKa, rejectNoTofDauTrackKa};
      default: {
        LOGP(fatal, "createV0Daughter: unknown daughterId");
        return {};
      }
    }
  }

  V0Mother createV0Mother(const int motherId)
  {
    switch (motherId) {
      case MotherGamma:
        return V0Mother{.posDaughterId = DaughterElectron, .negDaughterId = DaughterElectron};
      case MotherK0S:
        return V0Mother{.posDaughterId = DaughterPion, .negDaughterId = DaughterPion};
      case MotherLambda:
        return V0Mother{.posDaughterId = DaughterProton, .negDaughterId = DaughterPion};
      case MotherAntiLambda:
        return V0Mother{.posDaughterId = DaughterPion, .negDaughterId = DaughterProton};
      case MotherOmega:
        return V0Mother{.posDaughterId = DaughterUndef, .negDaughterId = DaughterKaon};
      case MotherAntiOmega:
        return V0Mother{.posDaughterId = DaughterKaon, .negDaughterId = DaughterUndef};
      default: {
        LOGP(fatal, "createV0Mother: unknown motherId");
        return {};
      }
    }
  }

  float getStrangenessTofNSigma(V0sWithID::iterator const& v0, const int motherId, const int daughterId, const bool isPositive)
  {
    if (motherId == MotherGamma && daughterId == DaughterElectron) {
      return NSigmaTofElectorn;
    }
    if (motherId == MotherK0S && daughterId == DaughterPion) {
      return isPositive ? v0.tofNSigmaK0PiPlus() : v0.tofNSigmaK0PiMinus();
    }
    if (motherId == MotherLambda) {
      if (daughterId == DaughterProton && isPositive) {
        return v0.tofNSigmaLaPr();
      }
      if (daughterId == DaughterPion && !isPositive) {
        return v0.tofNSigmaLaPi();
      }
    }
    if (motherId == MotherAntiLambda) {
      if (daughterId == DaughterProton && !isPositive) {
        return v0.tofNSigmaALaPr();
      }
      if (daughterId == DaughterPion && isPositive) {
        return v0.tofNSigmaALaPi();
      }
    }

    LOGP(fatal, "getStrangenessTofNSigma for V0: wrong combination of motherId, daughterId and sign");
    return UndefValueFloat;
  }

  float getStrangenessTofNSigma(CascsWithID::iterator const& casc, const int motherId, const int daughterId, bool)
  {
    if ((motherId == MotherOmega || motherId == MotherAntiOmega) && daughterId == DaughterKaon) {
      return casc.tofNSigmaOmKa();
    }

    LOGP(fatal, "getStrangenessTofNSigma for cascade: wrong combination of motherId and daughterId");
    return UndefValueFloat;
  }

  template <bool DoUseCorrectedDeDx, int ModeId, typename T, typename C, typename V0Casc>
  void fillSkimmedV0Table(V0Casc const& v0casc, T const& track, aod::TracksQA const& trackQA, const bool existTrkQA, C const& collision, const float nSigmaTPC, const float nSigmaTOF, const float nSigmaITS, const float dEdxExp, const o2::track::PID::ID id, const int runnumber, const float hadronicRate, const int bcGlobalIndex, const int bcTimeFrameId, const int bcBcInTimeFrame, const OccupancyValues& occValues, const bool isGoodRctEvent)
  {
    const double ncl = track.tpcNClsFound();
    const double nclPID = track.tpcNClsPID();
    const double nclShared = track.tpcNClsShared();
    const double p = track.tpcInnerParam();
    const double mass = o2::track::pid_constants::sMasses[id];
    const double bg = p / mass;
    const int multTPC = collision.multTPC();
    const auto trackOcc = collision.trackOccupancyInTimeRange();
    const auto ft0Occ = collision.ft0cOccupancyInTimeRange();
    const auto occMedianTime = collision.occupancyMedianTime();

    const float alpha = v0casc.alpha();
    const float qt = v0casc.qtarm();
    const float cosPA = getCosPA(v0casc, collision);
    const float pT = v0casc.pt();
    const float v0radius = getRadius(v0casc);
    const float gammapsipair = v0casc.psipair();

    const float usedDedx = tpcSignalGeneric<DoUseCorrectedDeDx>(track);
    float tpcdEdxNorm{UndefValueFloat};
    if constexpr (ModeId != ModeStandard) {
      tpcdEdxNorm = existTrkQA ? trackQA.tpcdEdxNorm() : UndefValueFloat;
    }
    if constexpr (ModeId == ModeStandard || ModeId == ModeWithdEdxTrkQA) {
      rowTPCTree(isGoodRctEvent,
                 usedDedx,
                 1. / dEdxExp,
                 track.tpcInnerParam(),
                 track.tgl(),
                 track.signed1Pt(),
                 track.eta(),
                 track.phi(),
                 track.y(),
                 mass,
                 bg,
                 multTPC / MultiplicityNorm,
                 std::sqrt(cfgsCommon.nClNorm / ncl),
                 nclPID,
                 nclShared,
                 id,
                 nSigmaTPC,
                 nSigmaTOF,
                 nSigmaITS,
                 runnumber,
                 trackOcc,
                 ft0Occ,
                 occMedianTime,
                 hadronicRate,
                 tpcdEdxNorm,
                 alpha,
                 qt,
                 cosPA,
                 pT,
                 v0radius,
                 gammapsipair);
    } else {
      rowTPCTreeWithTrkQA(isGoodRctEvent,
                          usedDedx,
                          1. / dEdxExp,
                          track.tpcInnerParam(),
                          track.tgl(),
                          track.signed1Pt(),
                          track.eta(),
                          track.phi(),
                          track.y(),
                          mass,
                          bg,
                          multTPC / MultiplicityNorm,
                          std::sqrt(cfgsCommon.nClNorm / ncl),
                          nclPID,
                          nclShared,
                          id,
                          nSigmaTPC,
                          nSigmaTOF,
                          nSigmaITS,
                          runnumber,
                          trackOcc,
                          ft0Occ,
                          occMedianTime,
                          hadronicRate,
                          tpcdEdxNorm,
                          alpha,
                          qt,
                          cosPA,
                          pT,
                          v0radius,
                          gammapsipair,
                          bcGlobalIndex,
                          bcTimeFrameId,
                          bcBcInTimeFrame,
                          existTrkQA ? trackQA.tpcClusterByteMask() : UndefValueInt,
                          existTrkQA ? trackQA.tpcdEdxMax0R() : UndefValueInt,
                          existTrkQA ? trackQA.tpcdEdxMax1R() : UndefValueInt,
                          existTrkQA ? trackQA.tpcdEdxMax2R() : UndefValueInt,
                          existTrkQA ? trackQA.tpcdEdxMax3R() : UndefValueInt,
                          existTrkQA ? trackQA.tpcdEdxTot0R() : UndefValueInt,
                          existTrkQA ? trackQA.tpcdEdxTot1R() : UndefValueInt,
                          existTrkQA ? trackQA.tpcdEdxTot2R() : UndefValueInt,
                          existTrkQA ? trackQA.tpcdEdxTot3R() : UndefValueInt,
                          occValues.tmoPrimUnfm80,
                          occValues.tmoFV0AUnfm80,
                          occValues.tmoFT0AUnfm80,
                          occValues.tmoFT0CUnfm80,
                          occValues.tmoRT0V0PrimUnfm80,
                          occValues.twmoPrimUnfm80,
                          occValues.twmoFV0AUnfm80,
                          occValues.twmoFT0AUnfm80,
                          occValues.twmoFT0CUnfm80,
                          occValues.twmoRT0V0PrimUnfm80);
    }
  } /// fillSkimmedV0Table

  /// Evaluate cosPA of v0 or casc
  template <V0OrCasc V0Casc, typename CollisionType>
  double getCosPA(V0Casc const& v0casc, CollisionType const& collision)
  {
    if constexpr (std::is_same_v<V0Casc, V0sWithID::iterator>) {
      return v0casc.v0cosPA();
    } else {
      return v0casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
    }
  }

  /// Evaluate radius of v0 or casc
  template <V0OrCasc V0Casc>
  double getRadius(V0Casc const& v0casc)
  {
    if constexpr (std::is_same_v<V0Casc, V0sWithID::iterator>) {
      return v0casc.v0radius();
    } else {
      return v0casc.cascradius();
    }
  }

  /// Evaluate add id of v0 or casc
  template <V0OrCasc V0Casc>
  int getAddId(V0Casc const& v0casc)
  {
    if constexpr (std::is_same_v<V0Casc, V0sWithID::iterator>) {
      return v0casc.v0addid();
    } else {
      return v0casc.cascaddid();
    }
  }

  template <bool IsCorrectedDeDx, int ModeId, typename TrksType, typename BCType, typename TrkQAType>
  void runV0(Colls const& collisions, TrksType const& myTracks, V0sWithID const& myV0s, CascsWithID const& myCascs, TrkQAType const& tracksQA)
  {
    if (collisions.size() != 0) {
      const auto timestamp = collisions.begin().bc_as<BCType>().timestamp();
      LOG(info) << "fRndm is seeded with timestamp of the 0-th bunch crossing (" << timestamp << ")";
      fRndm->SetSeed(timestamp);
    }

    constexpr bool IsWithTrackQa = ModeId != ModeStandard;

    std::vector<int64_t> labelTrack2TrackQA;
    if constexpr (IsWithTrackQa) {
      labelTrack2TrackQA.resize(myTracks.size(), -1);
      for (const auto& trackQA : tracksQA) {
        const int64_t trackId = trackQA.trackId();
        labelTrack2TrackQA.at(trackId) = trackQA.globalIndex();
      }
    }

    const auto tracksWithITSPid = soa::Attach<TrksType,
                                              aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaPi,
                                              aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr>(myTracks);

    int nV0Entries{0};
    int nCascEntries{0};

    const auto expectedOutputTableSize = static_cast<int64_t>(reserveV0Ratio * myV0s.size() + reserveCascRatio * myCascs.size());
    if constexpr (ModeId == ModeWithdEdxTrkQA || ModeId == ModeStandard) {
      rowTPCTree.reserve(expectedOutputTableSize);
    } else {
      rowTPCTreeWithTrkQA.reserve(expectedOutputTableSize);
    }

    std::string irSource{};
    bool isFirstCollision{true};
    for (const auto& collision : collisions) {
      if (!isEventSelected(collision, cfgsCommon.applyEvSel)) {
        continue;
      }
      const bool isGoodRctEvent = rctChecker.checkTable(collision);
      if (cfgsCommon.requireGoodRct && !isGoodRctEvent) {
        continue;
      }

      const auto v0s = myV0s.sliceBy(perCollisionV0s, static_cast<int>(collision.globalIndex()));
      const auto cascs = myCascs.sliceBy(perCollisionCascs, static_cast<int>(collision.globalIndex()));
      const auto bc = collision.bc_as<BCType>();
      if (isFirstCollision) {
        irSource = evaluateIrSource(ccdb, cfgsCommon.ccdbPathGrpLhcIf, bc.timestamp());
      }
      isFirstCollision = false;
      const int runnumber = bc.runNumber();
      const auto hadronicRate = !irSource.empty() ? mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * OneToKilo : 0.;
      const int bcGlobalIndex = bc.globalIndex();
      int bcTimeFrameId{}, bcBcInTimeFrame{};
      if constexpr (ModeId == ModeWithdEdxTrkQA || ModeId == ModeStandard) {
        bcTimeFrameId = UndefValueInt;
        bcBcInTimeFrame = UndefValueInt;
      } else if constexpr (ModeId == ModeWithTrkQA) {
        bcTimeFrameId = bc.tfId();
        bcBcInTimeFrame = bc.bcInTF();
      }

      auto getTrackQA = [&](const TrksType::iterator& track) {
        if constexpr (!IsWithTrackQa) {
          return std::make_pair(aod::TracksQA{}, false);
        } else {
          const auto trackGlobalIndex = track.globalIndex();
          const auto label = labelTrack2TrackQA.at(trackGlobalIndex);
          const bool existTrkQA = (label != -1);
          const int64_t trkIndex = existTrkQA ? label : 0;
          const aod::TracksQA& trkQA = tracksQA.iteratorAt(trkIndex);

          return std::make_pair(trkQA, existTrkQA);
        }
      };

      auto fillDaughterTrack = [&](const auto& mother, const TrksType::iterator& dauTrack, const auto& v0, const bool isPositive) {
        const auto [trackQAInstance, existTrkQA] = getTrackQA(dauTrack);
        const auto trackId = dauTrack.globalIndex();
        const auto dauTrackWithITSPid = tracksWithITSPid.rawIteratorAt(trackId);
        const auto v0Id = getAddId(v0);
        const V0Mother v0Mother = createV0Mother(v0Id);
        const auto daughterId = isPositive ? v0Mother.posDaughterId : v0Mother.negDaughterId;
        const V0Daughter daughter = createV0Daughter<IsCorrectedDeDx>(v0, dauTrackWithITSPid, v0Id, daughterId, isPositive);

        const bool passTrackSelection = isTrackSelected(dauTrack, cfgsCommon.trackSelection);
        const bool passNSigmaTofCut = std::fabs(daughter.tofNSigma) < daughter.nSigmaTofDauTrack || std::fabs(daughter.tofNSigma - NSigmaTofUnmatched) < nSigmaTofUnmatchedEqualityTolerance;
        const bool passMatchTofRequirement = !daughter.rejectNoTofDauTrack || std::fabs(daughter.tofNSigma - NSigmaTofUnmatched) > nSigmaTofUnmatchedEqualityTolerance;
        if (passTrackSelection && passNSigmaTofCut && passMatchTofRequirement) {
          if (cfgsCommon.saveMomentumHistograms) {
            const double log10p = std::log10(dauTrack.tpcInnerParam());
            const int pidShift = isPositive ? 0 : o2::track::PID::NIDs;
            registry.fill(HIST("hParticleVsMomentum"), log10p, pidShift + daughter.id);
          }
          if (cfgsCommon.saveTrees) {
            if (cfgsCommon.applyAdaptiveDownsampling && !performAdaptiveDownsampling(fRndm, dauTrack.tpcInnerParam(), daughter.id, mTargetYieldPerBin, mHistMomentum, isPositive)) {
              return false;
            }
            OccupancyValues occValues{};
            if constexpr (ModeId == ModeWithTrkQA) {
              evaluateOccupancyVariables(dauTrack, occValues);
            }
            fillSkimmedV0Table<IsCorrectedDeDx, ModeId>(mother, dauTrack, trackQAInstance, existTrkQA, collision, daughter.tpcNSigma, daughter.tofNSigma, daughter.itsNSigma, daughter.tpcExpSignal, daughter.id, runnumber, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame, occValues, isGoodRctEvent);
            return true;
          }
        }
        return false;
      };

      /// Loop over v0 candidates
      for (const auto& v0 : v0s) {
        const auto v0Id = v0.v0addid();
        if (v0Id == MotherUndef) {
          continue;
        }
        const auto posTrack = v0.posTrack_as<TrksType>();
        const auto negTrack = v0.negTrack_as<TrksType>();

        if (fillDaughterTrack(v0, posTrack, v0, true)) {
          ++nV0Entries;
        }
        if (fillDaughterTrack(v0, negTrack, v0, false)) {
          ++nV0Entries;
        }
      }

      /// Loop over cascade candidates
      for (const auto& casc : cascs) {
        const auto cascId = casc.cascaddid();
        if (cascId == MotherUndef) {
          continue;
        }
        const auto bachTrack = casc.bachelor_as<TrksType>();
        // Omega and antiomega
        const auto isDaughterPositive = cascId == MotherAntiOmega;
        if (fillDaughterTrack(casc, bachTrack, casc, isDaughterPositive)) {
          ++nCascEntries;
        }
      }
    }
    LOG(info) << "runV0() summary:";
    LOG(info) << "V0 table size = " << myV0s.size();
    LOG(info) << "Cascade table size = " << myCascs.size();
    LOG(info) << "nV0Entries = " << nV0Entries;
    LOG(info) << "nCascEntries = " << nCascEntries;
    LOG(info) << "nV0Entries / V0 table size = " << static_cast<double>(nV0Entries) / myV0s.size();
    LOG(info) << "nCascEntries / Cascade table size = " << static_cast<double>(nCascEntries) / myCascs.size();

    if (cfgsCommon.saveReserveQaHisto) {
      registry.fill(HIST("hV0OutputRatio"), static_cast<double>(nV0Entries) / myV0s.size());
      registry.fill(HIST("hCascOutputRatio"), static_cast<double>(nCascEntries) / myCascs.size());
    }
  } /// runV0

  void processStandard(Colls const& collisions,
                       Trks const& myTracks,
                       V0sWithID const& myV0s,
                       CascsWithID const& myCascs,
                       aod::BCsWithTimestamps const&)
  {
    runV0<false, ModeStandard, Trks, aod::BCsWithTimestamps>(collisions, myTracks, myV0s, myCascs, static_cast<TObject*>(nullptr));
  }
  PROCESS_SWITCH(TreeWriterTpcV0, processStandard, "Standard V0 Samples for PID", true);

  void processStandardWithCorrecteddEdx(Colls const& collisions,
                                        TrksWithDEdxCorrection const& myTracks,
                                        V0sWithID const& myV0s,
                                        CascsWithID const& myCascs,
                                        aod::BCsWithTimestamps const&)
  {
    runV0<true, ModeStandard, TrksWithDEdxCorrection, aod::BCsWithTimestamps>(collisions, myTracks, myV0s, myCascs, static_cast<TObject*>(nullptr));
  }
  PROCESS_SWITCH(TreeWriterTpcV0, processStandardWithCorrecteddEdx, "Standard V0 Samples for PID with corrected dEdx", false);

  void processWithdEdxTrQA(Colls const& collisions,
                           Trks const& myTracks,
                           V0sWithID const& myV0s,
                           CascsWithID const& myCascs,
                           aod::BCsWithTimestamps const&,
                           aod::TracksQAVersion const& tracksQA)
  {
    runV0<false, ModeWithdEdxTrkQA, Trks, aod::BCsWithTimestamps>(collisions, myTracks, myV0s, myCascs, tracksQA);
  }
  PROCESS_SWITCH(TreeWriterTpcV0, processWithdEdxTrQA, "Standard V0 Samples with dEdx from Track QA for PID", false);

  void processWithdEdxTrQAWithCorrecteddEdx(Colls const& collisions,
                                            TrksWithDEdxCorrection const& myTracks,
                                            V0sWithID const& myV0s,
                                            CascsWithID const& myCascs,
                                            aod::BCsWithTimestamps const&,
                                            aod::TracksQAVersion const& tracksQA)
  {
    runV0<true, ModeWithdEdxTrkQA, TrksWithDEdxCorrection, aod::BCsWithTimestamps>(collisions, myTracks, myV0s, myCascs, tracksQA);
  }
  PROCESS_SWITCH(TreeWriterTpcV0, processWithdEdxTrQAWithCorrecteddEdx, "Standard V0 Samples with dEdx from Track QA for PID with corrected dEdx", false);

  void processWithTrQA(Colls const& collisions,
                       TrksTmo const& myTracks,
                       V0sWithID const& myV0s,
                       CascsWithID const& myCascs,
                       MyBCTable const&,
                       aod::TracksQAVersion const& tracksQA,
                       TrackMeanOccs const&)
  {
    runV0<false, ModeWithTrkQA, TrksTmo, MyBCTable>(collisions, myTracks, myV0s, myCascs, tracksQA);
  }
  PROCESS_SWITCH(TreeWriterTpcV0, processWithTrQA, "Standard V0 Samples with Track QA for PID", false);

  void processWithTrQAWithCorrecteddEdx(Colls const& collisions,
                                        TrksTmoWithDEdxCorrection const& myTracks,
                                        V0sWithID const& myV0s,
                                        CascsWithID const& myCascs,
                                        MyBCTable const&,
                                        aod::TracksQAVersion const& tracksQA,
                                        TrackMeanOccs const&)
  {
    runV0<true, ModeWithTrkQA, TrksTmoWithDEdxCorrection, MyBCTable>(collisions, myTracks, myV0s, myCascs, tracksQA);
  }
  PROCESS_SWITCH(TreeWriterTpcV0, processWithTrQAWithCorrecteddEdx, "Standard V0 Samples with Track QA for PID with corrected dEdx", false);

  void processDummy(Colls const&) {}
  PROCESS_SWITCH(TreeWriterTpcV0, processDummy, "Dummy function", false);

}; /// struct TreeWriterTpcV0

struct TreeWriterTpcTof {

  Produces<o2::aod::SkimmedTPCTOFTree> rowTPCTOFTree;
  Produces<o2::aod::SkimmedTPCTOFTreeWithTrkQA> rowTPCTOFTreeWithTrkQA;

  CfgCommon cfgsCommon;

  /// Triton
  Configurable<float> maxMomTPCOnlyTr{"maxMomTPCOnlyTr", 1.5, "Maximum momentum for TPC only cut triton"};
  Configurable<float> maxMomHardCutOnlyTr{"maxMomHardCutOnlyTr", 50, "Maximum TPC inner momentum for triton"};
  Configurable<float> nSigmaTPCOnlyTr{"nSigmaTPCOnlyTr", 4., "number of sigma for TPC only cut triton"};
  Configurable<float> nSigmaTpcTpctofTr{"nSigmaTpcTpctofTr", 4., "number of sigma for TPC cut for TPC and TOF combined triton"};
  Configurable<float> nSigmaTofTpctofTr{"nSigmaTofTpctofTr", 3., "number of sigma for TOF cut for TPC and TOF combined triton"};
  /// Deuteron
  Configurable<float> maxMomTPCOnlyDe{"maxMomTPCOnlyDe", 1.0, "Maximum momentum for TPC only cut deuteron"};
  Configurable<float> maxMomHardCutOnlyDe{"maxMomHardCutOnlyDe", 50, "Maximum TPC inner momentum for deuteron"};
  Configurable<float> nSigmaTPCOnlyDe{"nSigmaTPCOnlyDe", 4., "number of sigma for TPC only cut deuteron"};
  Configurable<float> nSigmaTpcTpctofDe{"nSigmaTpcTpctofDe", 4., "number of sigma for TPC cut for TPC and TOF combined deuteron"};
  Configurable<float> nSigmaTofTpctofDe{"nSigmaTofTpctofDe", 3., "number of sigma for TOF cut for TPC and TOF combined deuteron"};
  /// Proton
  Configurable<float> maxMomTPCOnlyPr{"maxMomTPCOnlyPr", 0.6, "Maximum momentum for TPC only cut proton"};
  Configurable<float> maxMomHardCutOnlyPr{"maxMomHardCutOnlyPr", 50, "Maximum TPC inner momentum for proton"};
  Configurable<float> nSigmaTPCOnlyPr{"nSigmaTPCOnlyPr", 4., "number of sigma for TPC only cut proton"};
  Configurable<float> nSigmaTpcTpctofPr{"nSigmaTpcTpctofPr", 4., "number of sigma for TPC cut for TPC and TOF combined proton"};
  Configurable<float> nSigmaTofTpctofPr{"nSigmaTofTpctofPr", 3., "number of sigma for TOF cut for TPC and TOF combined proton"};
  /// Kaon
  Configurable<float> maxMomTPCOnlyKa{"maxMomTPCOnlyKa", 0.3, "Maximum momentum for TPC only cut kaon"};
  Configurable<float> maxMomHardCutOnlyKa{"maxMomHardCutOnlyKa", 50, "Maximum TPC inner momentum for kaons"};
  Configurable<float> nSigmaTPCOnlyKa{"nSigmaTPCOnlyKa", 4., "number of sigma for TPC only cut kaon"};
  Configurable<float> nSigmaTpcTpctofKa{"nSigmaTpcTpctofKa", 4., "number of sigma for TPC cut for TPC and TOF combined kaon"};
  Configurable<float> nSigmaTofTpctofKa{"nSigmaTofTpctofKa", 3., "number of sigma for TOF cut for TPC and TOF combined kaon"};
  /// Pion
  Configurable<float> maxMomTPCOnlyPi{"maxMomTPCOnlyPi", 0.5, "Maximum momentum for TPC only cut pion"};
  Configurable<float> maxMomHardCutOnlyPi{"maxMomHardCutOnlyPi", 50, "Maximum TPC inner momentum for pion"};
  Configurable<float> nSigmaTPCOnlyPi{"nSigmaTPCOnlyPi", 4., "number of sigma for TPC only cut pion"};
  Configurable<float> nSigmaTpcTpctofPi{"nSigmaTpcTpctofPi", 4., "number of sigma for TPC cut for TPC and TOF combined pion"};
  Configurable<float> nSigmaTofTpctofPi{"nSigmaTofTpctofPi", 4., "number of sigma for TOF cut for TPC and TOF combined pion"};
  // Configurable for output table reservation size
  Configurable<float> reserveTrackRatio{"reserveTrackRatio", 0.003, "Ratio of how many tracks are expected in the output table to the input Tracks table size"};
  // Configurables for run condtion table
  Configurable<std::string> adaptDwnsplHistoName{"adaptDwnsplHistoName", "tree-writer-tpc-tof/hParticleVsMomentum", "Name of the histogram with expected-momentum of particles"};

  HistogramRegistry registry{"registry", {}};

  struct TofTrack {
    double maxMomHardCutOnly{UndefValueDouble};
    double maxMomTPCOnly{UndefValueDouble};
    double tpcNSigma{UndefValueDouble};
    double nSigmaTPCOnly{UndefValueDouble};
    double tofNSigma{UndefValueDouble};
    double itsNSigma{UndefValueDouble};
    double tpcExpSignal{UndefValueDouble};
    o2::track::PID::ID pid{0};
    double nSigmaTofTpctof{UndefValueDouble};
    double nSigmaTpcTpctof{UndefValueDouble};
  };

  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  ctpRateFetcher mRateFetcher;

  o2::aod::rctsel::RCTFlagsChecker rctChecker;

  TRandom3* fRndm = new TRandom3();

  TH2* mHistMomentum{};

  double mTargetYieldPerBin{};

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra,
                         aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa,
                         aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr,
                         aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa,
                         aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr,
                         aod::TrackSelection>;
  using TrksWithDEdxCorrection = soa::Join<Trks, aod::DEdxsCorrected>;
  using Colls = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;
  using MyBCTable = soa::Join<aod::BCsWithTimestamps, aod::BCTFinfoTable>;
  using TrksTmo = soa::Join<Trks, aod::TrackToTmo>;
  using TrksTmoWithDEdxCorrection = soa::Join<Trks, aod::DEdxsCorrected, aod::TrackToTmo>;

  Preslice<Trks> perCollisionTracks = aod::track::collisionId;
  Preslice<TrksWithDEdxCorrection> perCollisionTracksWithCorrecteddEdx = aod::track::collisionId;
  Preslice<TrksTmo> perCollisionTracksTmo = aod::track::collisionId;
  Preslice<TrksTmoWithDEdxCorrection> perCollisionTracksTmoWithCorrecteddEdx = aod::track::collisionId;

  void init(o2::framework::InitContext&)
  {
    const std::array<bool, 7> doprocess{doprocessStandard, doprocessStandardWithCorrecteddEdx, doprocessWithdEdxTrQA, doprocessWithdEdxTrQAWithCorrecteddEdx, doprocessWithTrQA, doprocessWithTrQAWithCorrecteddEdx, doprocessDummy};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function should be enabled");
    }

    if (!(cfgsCommon.saveTrees || cfgsCommon.saveMomentumHistograms)) {
      LOGP(fatal, "At least one out of trees or momentum histograms should be saved");
    }

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);

    rctChecker.init(cfgsCommon.rctLabel, cfgsCommon.checkZdc, cfgsCommon.treatLimitedAcceptanceAsBad);

    if (cfgsCommon.saveReserveQaHisto) {
      registry.add("hTrackOutputRatio", "Track out/in ratio;Track out/in ratio;Entries", {HistType::kTH1F, {{100, 0, reserveTrackRatio}}});
    }

    if (cfgsCommon.saveMomentumHistograms) {
      registry.add("hParticleVsMomentum", "hParticleVsMomentum;log_{10}(#it{p} / (GeV/#it{c}));PID id", {HistType::kTH2D, {{cfgsCommon.axisMomentumHistograms}, {2 * o2::track::PID::NIDs, -0.5, 2 * o2::track::PID::NIDs - 0.5}}});
      for (int iParticle = 0; iParticle < o2::track::PID::NIDs; ++iParticle) {
        registry.get<TH2>(HIST("hParticleVsMomentum"))->GetYaxis()->SetBinLabel(iParticle + 1, (static_cast<std::string>(ParticleNames.at(iParticle)) + "^{+}").c_str());
        registry.get<TH2>(HIST("hParticleVsMomentum"))->GetYaxis()->SetBinLabel(iParticle + 1 + o2::track::PID::NIDs, (static_cast<std::string>(ParticleNames.at(iParticle)) + "^{-}").c_str());
      }
    }

    if (cfgsCommon.applyAdaptiveDownsampling && cfgsCommon.saveTrees) {
      initAdaptDwnsplHistograms(cfgsCommon.ccdbPathForAdaptDwnsplFile, cfgsCommon.adaptDwnsplFileName, adaptDwnsplHistoName, mHistMomentum);
      if (cfgsCommon.adaptDwnsplHistoScaleFactor != 1.f) {
        mHistMomentum->Scale(cfgsCommon.adaptDwnsplHistoScaleFactor);
      }
      mTargetYieldPerBin = cfgsCommon.adaptDwnsplTargetYield * mHistMomentum->GetXaxis()->GetBinWidth(1);
    }
  }

  template <bool DoCorrectDeDx, int ModeId, typename T, typename C>
  void fillSkimmedTpcTofTable(T const& track, aod::TracksQA const& trackQA, const bool existTrkQA, C const& collision, const float nSigmaTPC, const float nSigmaTOF, const float nSigmaITS, const float dEdxExp, const o2::track::PID::ID id, const int runnumber, const double hadronicRate, const int bcGlobalIndex, const int bcTimeFrameId, const int bcBcInTimeFrame, const OccupancyValues& occValues, const bool isGoodRctEvent)
  {
    const double ncl = track.tpcNClsFound();
    const double nclPID = track.tpcNClsPID();
    const double nclShared = track.tpcNClsShared();
    const double p = track.tpcInnerParam();
    const double mass = o2::track::pid_constants::sMasses[id];
    const double bg = p / mass;
    const int multTPC = collision.multTPC();
    const auto trackOcc = collision.trackOccupancyInTimeRange();
    const auto ft0Occ = collision.ft0cOccupancyInTimeRange();
    const auto occMedianTime = collision.occupancyMedianTime();

    const float usedEdx = tpcSignalGeneric<DoCorrectDeDx>(track);
    float tpcdEdxNorm{UndefValueFloat};
    if constexpr (ModeId != ModeStandard) {
      tpcdEdxNorm = existTrkQA ? trackQA.tpcdEdxNorm() : UndefValueFloat;
    }
    if (ModeId == ModeStandard || ModeId == ModeWithdEdxTrkQA) {
      rowTPCTOFTree(isGoodRctEvent,
                    usedEdx,
                    1. / dEdxExp,
                    track.tpcInnerParam(),
                    track.tgl(),
                    track.signed1Pt(),
                    track.eta(),
                    track.phi(),
                    track.y(),
                    mass,
                    bg,
                    multTPC / MultiplicityNorm,
                    std::sqrt(cfgsCommon.nClNorm / ncl),
                    nclPID,
                    nclShared,
                    id,
                    nSigmaTPC,
                    nSigmaTOF,
                    nSigmaITS,
                    runnumber,
                    trackOcc,
                    ft0Occ,
                    occMedianTime,
                    hadronicRate,
                    tpcdEdxNorm);
    } else {
      rowTPCTOFTreeWithTrkQA(isGoodRctEvent,
                             usedEdx,
                             1. / dEdxExp,
                             track.tpcInnerParam(),
                             track.tgl(),
                             track.signed1Pt(),
                             track.eta(),
                             track.phi(),
                             track.y(),
                             mass,
                             bg,
                             multTPC / MultiplicityNorm,
                             std::sqrt(cfgsCommon.nClNorm / ncl),
                             nclPID,
                             nclShared,
                             id,
                             nSigmaTPC,
                             nSigmaTOF,
                             nSigmaITS,
                             runnumber,
                             trackOcc,
                             ft0Occ,
                             occMedianTime,
                             hadronicRate,
                             tpcdEdxNorm,
                             bcGlobalIndex,
                             bcTimeFrameId,
                             bcBcInTimeFrame,
                             existTrkQA ? trackQA.tpcClusterByteMask() : UndefValueInt,
                             existTrkQA ? trackQA.tpcdEdxMax0R() : UndefValueInt,
                             existTrkQA ? trackQA.tpcdEdxMax1R() : UndefValueInt,
                             existTrkQA ? trackQA.tpcdEdxMax2R() : UndefValueInt,
                             existTrkQA ? trackQA.tpcdEdxMax3R() : UndefValueInt,
                             existTrkQA ? trackQA.tpcdEdxTot0R() : UndefValueInt,
                             existTrkQA ? trackQA.tpcdEdxTot1R() : UndefValueInt,
                             existTrkQA ? trackQA.tpcdEdxTot2R() : UndefValueInt,
                             existTrkQA ? trackQA.tpcdEdxTot3R() : UndefValueInt,
                             occValues.tmoPrimUnfm80,
                             occValues.tmoFV0AUnfm80,
                             occValues.tmoFT0AUnfm80,
                             occValues.tmoFT0CUnfm80,
                             occValues.tmoRT0V0PrimUnfm80,
                             occValues.twmoPrimUnfm80,
                             occValues.twmoFV0AUnfm80,
                             occValues.twmoFT0AUnfm80,
                             occValues.twmoFT0CUnfm80,
                             occValues.twmoRT0V0PrimUnfm80);
    }
  } /// fillSkimmedTpcTofTable

  template <bool IsCorrectedDeDx, int ModeId, typename TrksType, typename BCType, typename TrkQAType, typename TMOType>
  void runTof(Colls const& collisions, TrksType const& myTracks, TrkQAType const& tracksQA, Preslice<TrksType> const& perCollisionTracksType, TMOType const& trackMeanOccs)
  {
    if (collisions.size() != 0) {
      const auto timestamp = collisions.begin().bc_as<BCType>().timestamp();
      LOG(info) << "fRndm is seeded with timestamp of the 0-th bunch crossing (" << timestamp << ")";
      fRndm->SetSeed(timestamp);
    }

    constexpr bool IsWithTrackQa = ModeId != ModeStandard;

    std::vector<int64_t> labelTrack2TrackQA;
    if constexpr (IsWithTrackQa) {
      labelTrack2TrackQA.resize(myTracks.size(), -1);
      for (const auto& trackQA : tracksQA) {
        const int64_t trackId = trackQA.trackId();
        labelTrack2TrackQA.at(trackId) = trackQA.globalIndex();
      }
    }

    const int64_t expectedOutputTableSize = static_cast<int64_t>(reserveTrackRatio * myTracks.size());
    if constexpr (ModeId == ModeWithdEdxTrkQA || ModeId == ModeStandard) {
      rowTPCTOFTree.reserve(expectedOutputTableSize);
    } else {
      rowTPCTOFTreeWithTrkQA.reserve(expectedOutputTableSize);
    }

    std::string irSource{};
    bool isFirstCollision{true};
    for (const auto& collision : collisions) {
      const auto tracks = myTracks.sliceBy(perCollisionTracksType, collision.globalIndex());
      if (!isEventSelected(collision, cfgsCommon.applyEvSel)) {
        continue;
      }
      const bool isGoodRctEvent = rctChecker.checkTable(collision);
      if (cfgsCommon.requireGoodRct && !isGoodRctEvent) {
        continue;
      }

      auto tracksWithITSPid = soa::Attach<TrksType,
                                          aod::pidits::ITSNSigmaPi, aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr,
                                          aod::pidits::ITSNSigmaDe, aod::pidits::ITSNSigmaTr>(tracks);

      if constexpr (ModeId == ModeWithTrkQA) {
        tracksWithITSPid.bindExternalIndices(&trackMeanOccs);
      }

      const auto bc = collision.bc_as<BCType>();
      if (isFirstCollision) {
        irSource = evaluateIrSource(ccdb, cfgsCommon.ccdbPathGrpLhcIf, bc.timestamp());
      }
      isFirstCollision = false;
      const int runnumber = bc.runNumber();
      const auto hadronicRate = !irSource.empty() ? mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * OneToKilo : 0.;
      const int bcGlobalIndex = bc.globalIndex();
      int bcTimeFrameId{}, bcBcInTimeFrame{};
      if constexpr (ModeId == ModeStandard || ModeId == ModeWithdEdxTrkQA) {
        bcTimeFrameId = UndefValueInt;
        bcBcInTimeFrame = UndefValueInt;
      } else {
        bcTimeFrameId = bc.tfId();
        bcBcInTimeFrame = bc.bcInTF();
      }
      for (auto const& trk : tracksWithITSPid) {
        if (!isTrackSelected(trk, cfgsCommon.trackSelection)) {
          continue;
        }
        // get the corresponding trackQA using labelTracks2TracKQA and get variables of interest
        aod::TracksQA trackQA{};
        bool existTrkQA{false};
        if constexpr (IsWithTrackQa) {
          const auto label = labelTrack2TrackQA.at(trk.globalIndex());
          existTrkQA = (label != -1);
          const int64_t trkIndex = existTrkQA ? label : 0;
          trackQA = tracksQA.iteratorAt(trkIndex);
        }

        TofTrack tofTriton(maxMomHardCutOnlyTr, maxMomTPCOnlyTr, trk.tpcNSigmaTr(), nSigmaTPCOnlyTr, trk.tofNSigmaTr(), trk.itsNSigmaTr(), trk.tpcExpSignalTr(tpcSignalGeneric<IsCorrectedDeDx>(trk)), PidTriton, nSigmaTofTpctofTr, nSigmaTpcTpctofTr);

        TofTrack tofDeuteron(maxMomHardCutOnlyDe, maxMomTPCOnlyDe, trk.tpcNSigmaDe(), nSigmaTPCOnlyDe, trk.tofNSigmaDe(), trk.itsNSigmaDe(), trk.tpcExpSignalDe(tpcSignalGeneric<IsCorrectedDeDx>(trk)), PidDeuteron, nSigmaTofTpctofDe, nSigmaTpcTpctofDe);

        TofTrack tofProton(maxMomHardCutOnlyPr, maxMomTPCOnlyPr, trk.tpcNSigmaPr(), nSigmaTPCOnlyPr, trk.tofNSigmaPr(), trk.itsNSigmaPr(), trk.tpcExpSignalPr(tpcSignalGeneric<IsCorrectedDeDx>(trk)), PidProton, nSigmaTofTpctofPr, nSigmaTpcTpctofPr);

        TofTrack tofKaon(maxMomHardCutOnlyKa, maxMomTPCOnlyKa, trk.tpcNSigmaKa(), nSigmaTPCOnlyKa, trk.tofNSigmaKa(), trk.itsNSigmaKa(), trk.tpcExpSignalKa(tpcSignalGeneric<IsCorrectedDeDx>(trk)), PidKaon, nSigmaTofTpctofKa, nSigmaTpcTpctofKa);

        TofTrack tofPion(maxMomHardCutOnlyPi, maxMomTPCOnlyPi, trk.tpcNSigmaPi(), nSigmaTPCOnlyPi, trk.tofNSigmaPi(), trk.itsNSigmaPi(), trk.tpcExpSignalPi(tpcSignalGeneric<IsCorrectedDeDx>(trk)), PidPion, nSigmaTofTpctofPi, nSigmaTpcTpctofPi);

        OccupancyValues occValues{};
        if constexpr (ModeId == ModeWithTrkQA) {
          evaluateOccupancyVariables(trk, occValues);
        }

        for (const auto& tofTrack : {&tofTriton, &tofDeuteron, &tofProton, &tofKaon, &tofPion}) {
          const bool passMomHardCut = trk.tpcInnerParam() < tofTrack->maxMomHardCutOnly;
          const bool passMomTpcOnly = trk.tpcInnerParam() <= tofTrack->maxMomTPCOnly && std::fabs(tofTrack->tpcNSigma) < tofTrack->nSigmaTPCOnly;
          const bool passMomTpcTof = trk.tpcInnerParam() > tofTrack->maxMomTPCOnly && std::fabs(tofTrack->tofNSigma) < tofTrack->nSigmaTofTpctof && std::fabs(tofTrack->tpcNSigma) < tofTrack->nSigmaTpcTpctof;
          const bool isPositive = (trk.signed1Pt() > 0);
          if (passMomHardCut && (passMomTpcOnly || passMomTpcTof)) {
            if (cfgsCommon.saveMomentumHistograms) {
              const double log10p = std::log10(trk.tpcInnerParam());
              const int pidShift = isPositive ? 0 : o2::track::PID::NIDs;
              registry.fill(HIST("hParticleVsMomentum"), log10p, pidShift + tofTrack->pid);
            }
            if (cfgsCommon.saveTrees) {
              if (cfgsCommon.applyAdaptiveDownsampling && !performAdaptiveDownsampling(fRndm, trk.tpcInnerParam(), tofTrack->pid, mTargetYieldPerBin, mHistMomentum, isPositive)) {
                continue;
              }
              fillSkimmedTpcTofTable<IsCorrectedDeDx, ModeId>(trk, trackQA, existTrkQA, collision, tofTrack->tpcNSigma, tofTrack->tofNSigma, tofTrack->itsNSigma, tofTrack->tpcExpSignal, tofTrack->pid, runnumber, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame, occValues, isGoodRctEvent);
            }
          }
        }
      } /// Loop tracks
    }
    LOG(info) << "runTof() summary:";
    LOG(info) << "Track table size = " << myTracks.size();
    LOG(info) << "nTrackEntries = " << rowTPCTOFTree.lastIndex() + 1;
    LOG(info) << "nTrackEntries / Track table size = " << static_cast<double>((rowTPCTOFTree.lastIndex() + 1)) / myTracks.size();

    if (cfgsCommon.saveReserveQaHisto) {
      registry.fill(HIST("hTrackOutputRatio"), static_cast<double>((rowTPCTOFTree.lastIndex() + 1)) / myTracks.size());
    }
  } /// runTof

  void processStandard(Colls const& collisions,
                       Trks const& myTracks,
                       aod::BCsWithTimestamps const&)
  {
    runTof<false, ModeStandard, Trks, aod::BCsWithTimestamps>(collisions, myTracks, static_cast<TObject*>(nullptr), perCollisionTracks, static_cast<TObject*>(nullptr));
  }
  PROCESS_SWITCH(TreeWriterTpcTof, processStandard, "Standard Samples for PID", true);

  void processStandardWithCorrecteddEdx(Colls const& collisions,
                                        TrksWithDEdxCorrection const& myTracks,
                                        aod::BCsWithTimestamps const&)
  {
    runTof<true, ModeStandard, TrksWithDEdxCorrection, aod::BCsWithTimestamps>(collisions, myTracks, static_cast<TObject*>(nullptr), perCollisionTracksWithCorrecteddEdx, static_cast<TObject*>(nullptr));
  }
  PROCESS_SWITCH(TreeWriterTpcTof, processStandardWithCorrecteddEdx, "Standard Samples for PID with corrected dEdx", false);

  void processWithdEdxTrQA(Colls const& collisions,
                           Trks const& myTracks,
                           aod::BCsWithTimestamps const&,
                           aod::TracksQAVersion const& tracksQA)
  {
    runTof<false, ModeWithdEdxTrkQA, Trks, aod::BCsWithTimestamps>(collisions, myTracks, tracksQA, perCollisionTracks, static_cast<TObject*>(nullptr));
  }
  PROCESS_SWITCH(TreeWriterTpcTof, processWithdEdxTrQA, "Samples for PID with TrackQA info", false);

  void processWithdEdxTrQAWithCorrecteddEdx(Colls const& collisions,
                                            TrksWithDEdxCorrection const& myTracks,
                                            aod::BCsWithTimestamps const&,
                                            aod::TracksQAVersion const& tracksQA)
  {
    runTof<true, ModeWithdEdxTrkQA, TrksWithDEdxCorrection, aod::BCsWithTimestamps>(collisions, myTracks, tracksQA, perCollisionTracksWithCorrecteddEdx, static_cast<TObject*>(nullptr));
  }
  PROCESS_SWITCH(TreeWriterTpcTof, processWithdEdxTrQAWithCorrecteddEdx, "Samples for PID with TrackQA info with corrected dEdx", false);

  void processWithTrQA(Colls const& collisions,
                       TrksTmo const& myTracks,
                       MyBCTable const&,
                       aod::TracksQAVersion const& tracksQA,
                       TrackMeanOccs const& trackMeanOccs)
  {
    runTof<false, ModeWithTrkQA, TrksTmo, MyBCTable>(collisions, myTracks, tracksQA, perCollisionTracksTmo, trackMeanOccs);
  }
  PROCESS_SWITCH(TreeWriterTpcTof, processWithTrQA, "Samples for PID with TrackQA info", false);

  void processWithTrQAWithCorrecteddEdx(Colls const& collisions,
                                        TrksTmoWithDEdxCorrection const& myTracks,
                                        MyBCTable const&,
                                        aod::TracksQAVersion const& tracksQA,
                                        TrackMeanOccs const& trackMeanOccs)
  {
    runTof<true, ModeWithTrkQA, TrksTmoWithDEdxCorrection, MyBCTable>(collisions, myTracks, tracksQA, perCollisionTracksTmoWithCorrecteddEdx, trackMeanOccs);
  }
  PROCESS_SWITCH(TreeWriterTpcTof, processWithTrQAWithCorrecteddEdx, "Samples for PID with TrackQA info with correced dEdx", false);

  void processDummy(Colls const&) {}
  PROCESS_SWITCH(TreeWriterTpcTof, processDummy, "Dummy function", false);

}; /// struct TreeWriterTpcTof
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<TreeWriterTpcTof>(cfgc)};
  workflow.push_back(adaptAnalysisTask<TreeWriterTpcV0>(cfgc));
  return workflow;
}
