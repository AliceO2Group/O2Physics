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

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/OccupancyTables.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTPCBase.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>

#include <TRandom3.h>

#include <cmath>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::constants::physics;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::dataformats;
using namespace o2::dpg_tpcskimstablecreator;

struct TreeWriterTpcV0 {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  using Trks = soa::Join<aod::Tracks, aod::V0Bits, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TrackSelection>;
  using TrksWithDEdxCorrection = soa::Join<aod::Tracks, aod::V0Bits, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TrackSelection, aod::DEdxsCorrected>;
  using Colls = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;
  using MyBCTable = soa::Join<aod::BCsWithTimestamps, aod::BCTFinfoTable>;
  using V0sWithID = soa::Join<aod::V0Datas, aod::V0MapID, aod::V0TOFNSigmas>;
  using CascsWithID = soa::Join<aod::CascDatas, aod::CascMapID, aod::CascTOFNSigmas>;

  /// Tables to be produced
  Produces<o2::aod::SkimmedTPCV0Tree> rowTPCTree;
  Produces<o2::aod::SkimmedTPCV0TreeWithdEdxTrkQA> rowTPCTreeWithdEdxTrkQA;
  Produces<o2::aod::SkimmedTPCV0TreeWithTrkQA> rowTPCTreeWithTrkQA;

  constexpr static o2::track::PID::ID PidElectron{o2::track::PID::Electron};
  constexpr static o2::track::PID::ID PidPion{o2::track::PID::Pion};
  constexpr static o2::track::PID::ID PidKaon{o2::track::PID::Kaon};
  constexpr static o2::track::PID::ID PidProton{o2::track::PID::Proton};

  // an arbitrary value of N sigma TOF assigned by TOF task to tracks which are not matched to TOF hits
  constexpr static float NSigmaTofUnmatched{o2::aod::v0data::kNoTOFValue};
  constexpr static float NSigmaTofUnmatchedEqualityTolerance{std::fabs(NSigmaTofUnmatched) / 1e4f};

  // an arbitrary value of "N sigma TOF" assigned to electorns (for uniformity reasons)
  constexpr static float NSigmaTofElectorn{1000.f};

  /// Configurables
  Configurable<float> nSigmaTofDauTrackPi{"nSigmaTofDauTrackPi", 999.f, "n-sigma TOF cut on the pion daughter tracks"};
  Configurable<float> nSigmaTofDauTrackPr{"nSigmaTofDauTrackPr", 999.f, "n-sigma TOF cut on the proton daughter tracks"};
  Configurable<float> nSigmaTofDauTrackKa{"nSigmaTofDauTrackKa", 999.f, "n-sigma TOF cut on the kaon daughter tracks"};
  Configurable<bool> rejectNoTofDauTrackPi{"rejectNoTofDauTrackPi", false, "reject not matched to TOF pion daughter tracks"};
  Configurable<bool> rejectNoTofDauTrackPr{"rejectNoTofDauTrackPr", false, "reject not matched to TOF proton daughter tracks"};
  Configurable<bool> rejectNoTofDauTrackKa{"rejectNoTofDauTrackKa", false, "reject not matched to TOF kaon daughter tracks"};
  Configurable<float> nClNorm{"nClNorm", 152., "Number of cluster normalization. Run 2: 159, Run 3 152"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<int> trackSelection{"trackSelection", 0, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};
  /// Configurables downsampling
  Configurable<double> dwnSmplFactorPi{"dwnSmplFactorPi", 1., "downsampling factor for pions, default fraction to keep is 1."};
  Configurable<double> dwnSmplFactorPr{"dwnSmplFactorPr", 1., "downsampling factor for protons, default fraction to keep is 1."};
  Configurable<double> dwnSmplFactorEl{"dwnSmplFactorEl", 1., "downsampling factor for electrons, default fraction to keep is 1."};
  Configurable<double> dwnSmplFactorKa{"dwnSmplFactorKa", 1., "downsampling factor for kaons, default fraction to keep is 1."};
  Configurable<float> sqrtSNN{"sqrtSNN", 5360., "sqrt(s_NN), used for downsampling with the Tsallis distribution"};
  Configurable<float> downsamplingTsalisPions{"downsamplingTsalisPions", -1., "Downsampling factor to reduce the number of pions"};
  Configurable<float> downsamplingTsalisProtons{"downsamplingTsalisProtons", -1., "Downsampling factor to reduce the number of protons"};
  Configurable<float> downsamplingTsalisElectrons{"downsamplingTsalisElectrons", -1., "Downsampling factor to reduce the number of electrons"};
  Configurable<float> downsamplingTsalisKaons{"downsamplingTsalisKaons", -1., "Downsampling factor to reduce the number of kaons"};
  Configurable<float> maxPt4dwnsmplTsalisPions{"maxPt4dwnsmplTsalisPions", 100., "Maximum Pt for applying downsampling factor of pions"};
  Configurable<float> maxPt4dwnsmplTsalisProtons{"maxPt4dwnsmplTsalisProtons", 100., "Maximum Pt for applying downsampling factor of protons"};
  Configurable<float> maxPt4dwnsmplTsalisElectrons{"maxPt4dwnsmplTsalisElectrons", 100., "Maximum Pt for applying  downsampling factor of electrons"};
  Configurable<float> maxPt4dwnsmplTsalisKaons{"maxPt4dwnsmplTsalisKaons", 100., "Maximum Pt for applying  downsampling factor of kaons"};

  enum { // Reconstructed V0 and cascade
    MotherUndef = -1,
    MotherGamma = 0,
    MotherK0S,
    MotherLambda,
    MotherAntiLambda,
    MotherOmega,
    MotherAntiOmega
  };

  enum {
    DaughterElectron = 0,
    DaughterPion,
    DaughterKaon,
    DaughterProton
  };

  ctpRateFetcher mRateFetcher;

  struct V0Daughter {
    double downsamplingTsalis{UndefValueDouble};
    double mass{UndefValueDouble};
    double maxPt4dwnsmplTsalis{UndefValueDouble};
    double tpcNSigma{UndefValueDouble};
    double tofNSigma{UndefValueDouble};
    double tpcExpSignal{UndefValueDouble};
    o2::track::PID::ID id{0};
    double dwnSmplFactor{UndefValueDouble};
    double nSigmaTofDauTrack{UndefValueDouble};
    bool rejectNoTofDauTrack{false};
  };

  template <bool IsCorrectedDeDx, typename V0Casc, typename T>
  V0Daughter createV0Daughter(const V0Casc& v0Casc, const T& track, const int motherId, const int daughterId, const bool isPositive = true)
  {
    switch (daughterId) {
      case DaughterElectron:
        return V0Daughter{downsamplingTsalisElectrons, MassElectron, maxPt4dwnsmplTsalisElectrons, track.tpcNSigmaEl(), getStrangenessTofNSigma(v0Casc, motherId, daughterId, isPositive), track.tpcExpSignalEl(tpcSignalGeneric<IsCorrectedDeDx>(track)), PidElectron, dwnSmplFactorEl, NSigmaTofElectorn + 1.f, false};
      case DaughterPion:
        return V0Daughter{downsamplingTsalisPions, MassPiPlus, maxPt4dwnsmplTsalisPions, track.tpcNSigmaPi(), getStrangenessTofNSigma(v0Casc, motherId, daughterId, isPositive), track.tpcExpSignalPi(tpcSignalGeneric<IsCorrectedDeDx>(track)), PidPion, dwnSmplFactorPi, nSigmaTofDauTrackPi, rejectNoTofDauTrackPi};
      case DaughterProton:
        return V0Daughter{downsamplingTsalisProtons, MassProton, maxPt4dwnsmplTsalisProtons, track.tpcNSigmaPr(), getStrangenessTofNSigma(v0Casc, motherId, daughterId, isPositive), track.tpcExpSignalPr(tpcSignalGeneric<IsCorrectedDeDx>(track)), PidProton, dwnSmplFactorPr, nSigmaTofDauTrackPr, rejectNoTofDauTrackPr};
      case DaughterKaon:
        return V0Daughter{downsamplingTsalisKaons, MassKPlus, maxPt4dwnsmplTsalisKaons, track.tpcNSigmaKa(), getStrangenessTofNSigma(v0Casc, motherId, daughterId, isPositive), track.tpcExpSignalKa(tpcSignalGeneric<IsCorrectedDeDx>(track)), PidKaon, dwnSmplFactorKa, nSigmaTofDauTrackKa, rejectNoTofDauTrackKa};
      default: {
        LOGP(fatal, "createV0Daughter: unknown daughterId");
        return V0Daughter();
      }
    }
  }

  struct V0Mother {
    int posDaughterId{UndefValueInt};
    int negDaughterId{UndefValueInt};
  };

  V0Mother createV0Mother(const int motherId)
  {
    switch (motherId) {
      case MotherGamma:
        return V0Mother{DaughterElectron, DaughterElectron};
      case MotherK0S:
        return V0Mother{DaughterPion, DaughterPion};
      case MotherLambda:
        return V0Mother{DaughterProton, DaughterPion};
      case MotherAntiLambda:
        return V0Mother{DaughterPion, DaughterProton};
      default: {
        LOGP(fatal, "createV0Mother: unknown motherId");
        return V0Mother();
      }
    }
  }

  float getStrangenessTofNSigma(V0sWithID::iterator const& v0, const int motherId, const int daughterId, const bool isPositive)
  {
    if (motherId == MotherGamma && daughterId == DaughterElectron) {
      return NSigmaTofElectorn;
    } else if (motherId == MotherK0S && daughterId == DaughterPion) {
      if (isPositive)
        return v0.tofNSigmaK0PiPlus();
      else
        return v0.tofNSigmaK0PiMinus();
    } else if (motherId == MotherLambda) {
      if (daughterId == DaughterProton && isPositive)
        return v0.tofNSigmaLaPr();
      else if (daughterId == DaughterPion && !isPositive)
        return v0.tofNSigmaLaPi();
    } else if (motherId == MotherAntiLambda) {
      if (daughterId == DaughterProton && !isPositive)
        return v0.tofNSigmaALaPr();
      else if (daughterId == DaughterPion && isPositive)
        return v0.tofNSigmaALaPi();
    }

    LOGP(fatal, "getStrangenessTofNSigma for V0: wrong combination of motherId, daughterId and sign");
    return UndefValueFloat;
  }

  float getStrangenessTofNSigma(CascsWithID::iterator const& casc, const int motherId, const int daughterId, bool)
  {
    if ((motherId == MotherOmega || motherId == MotherAntiOmega) && daughterId == DaughterKaon)
      return casc.tofNSigmaOmKa();

    LOGP(fatal, "getStrangenessTofNSigma for cascade: wrong combination of motherId and daughterId");
    return UndefValueFloat;
  }

  template <bool DoUseCorrectedDeDx, int ModeId, typename T, typename C, typename V0Casc>
  void fillSkimmedV0Table(V0Casc const& v0casc, T const& track, aod::TracksQA const& trackQA, const bool existTrkQA, C const& collision, const float nSigmaTPC, const float nSigmaTOF, const float dEdxExp, const o2::track::PID::ID id, const int runnumber, const double dwnSmplFactor, const float hadronicRate, const int bcGlobalIndex, const int bcTimeFrameId, const int bcBcInTimeFrame)
  {
    const double ncl = track.tpcNClsFound();
    const double nclPID = track.tpcNClsFindableMinusPID();
    const double p = track.tpcInnerParam();
    const double mass = o2::track::pid_constants::sMasses[id];
    const double bg = p / mass;
    const int multTPC = collision.multTPC();
    const auto trackOcc = collision.trackOccupancyInTimeRange();
    const auto ft0Occ = collision.ft0cOccupancyInTimeRange();

    const float alpha = v0casc.alpha();
    const float qt = v0casc.qtarm();
    const float cosPA = getCosPA(v0casc, collision);
    const float pT = v0casc.pt();
    const float v0radius = getRadius(v0casc);
    const float gammapsipair = v0casc.psipair();

    const double pseudoRndm = track.pt() * 1000. - static_cast<int64_t>(track.pt() * 1000);
    if (pseudoRndm < dwnSmplFactor) {
      float usedDedx;
      if constexpr (DoUseCorrectedDeDx) {
        usedDedx = track.tpcSignalCorrected();
      } else {
        usedDedx = track.tpcSignal();
      }
      if constexpr (ModeId == ModeStandard) {
        rowTPCTree(usedDedx,
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
                   std::sqrt(nClNorm / ncl),
                   nclPID,
                   id,
                   nSigmaTPC,
                   nSigmaTOF,
                   runnumber,
                   trackOcc,
                   ft0Occ,
                   hadronicRate,
                   alpha,
                   qt,
                   cosPA,
                   pT,
                   v0radius,
                   gammapsipair);
      } else if constexpr (ModeId == ModeWithdEdxTrkQA) {
        rowTPCTreeWithdEdxTrkQA(usedDedx,
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
                                std::sqrt(nClNorm / ncl),
                                nclPID,
                                id,
                                nSigmaTPC,
                                nSigmaTOF,
                                runnumber,
                                trackOcc,
                                ft0Occ,
                                hadronicRate,
                                alpha,
                                qt,
                                cosPA,
                                pT,
                                v0radius,
                                gammapsipair,
                                existTrkQA ? trackQA.tpcdEdxNorm() : UndefValueFloat);
      } else if constexpr (ModeId == ModeWithTrkQA) {
        rowTPCTreeWithTrkQA(usedDedx,
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
                            std::sqrt(nClNorm / ncl),
                            nclPID,
                            id,
                            nSigmaTPC,
                            nSigmaTOF,
                            runnumber,
                            trackOcc,
                            ft0Occ,
                            hadronicRate,
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
                            existTrkQA ? trackQA.tpcdEdxNorm() : UndefValueFloat);
      }
    }
  }

  TRandom3* fRndm = new TRandom3(0);

  void init(o2::framework::InitContext&)
  {
    const std::array<bool, 7> doprocess{doprocessStandard, doprocessStandardWithCorrecteddEdx, doprocessWithdEdxTrQA, doprocessWithdEdxTrQAWithCorrecteddEdx, doprocessWithTrQA, doprocessWithTrQAWithCorrecteddEdx, doprocessDummy};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function should be enabled");
    }

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);
  }

  /// Evaluate cosPA of the v0
  template <typename CollisionType>
  double getCosPA(V0sWithID::iterator const& v0, CollisionType const&)
  {
    return v0.v0cosPA();
  }

  /// Evaluate cosPA of the cascade
  template <typename CollisionType>
  double getCosPA(CascsWithID::iterator const& casc, CollisionType const& collision)
  {
    return casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
  }

  /// Evaluate radius of the v0
  double getRadius(V0sWithID::iterator const& v0)
  {
    return v0.v0radius();
  }

  /// Evaluate radius of the cascade
  double getRadius(CascsWithID::iterator const& casc)
  {
    return casc.cascradius();
  }

  Preslice<V0sWithID> perCollisionV0s = aod::v0data::collisionId;
  Preslice<CascsWithID> perCollisionCascs = aod::cascdata::collisionId;

  template <bool IsCorrectedDeDx, int ModeId, typename TrksType, typename BCType, typename TrkQAType>
  void runV0(Colls const& collisions, TrksType const& myTracks, V0sWithID const& myV0s, CascsWithID const& myCascs, TrkQAType const& tracksQA)
  {
    constexpr bool IsWithTrackQa = ModeId != ModeStandard;

    std::vector<int64_t> labelTrack2TrackQA;
    if constexpr (IsWithTrackQa) {
      labelTrack2TrackQA.resize(myTracks.size(), -1);
      for (const auto& trackQA : tracksQA) {
        const int64_t trackId = trackQA.trackId();
        labelTrack2TrackQA.at(trackId) = trackQA.globalIndex();
      }
    }
    for (const auto& collision : collisions) {
      /// Check event slection
      if (!isEventSelected(collision, applyEvSel)) {
        continue;
      }
      const auto& v0s = myV0s.sliceBy(perCollisionV0s, collision.globalIndex());
      const auto& cascs = myCascs.sliceBy(perCollisionCascs, collision.globalIndex());
      const auto& bc = collision.bc_as<BCType>();
      const int runnumber = bc.runNumber();
      const float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * OneToKilo;
      const int bcGlobalIndex = bc.globalIndex();
      int bcTimeFrameId, bcBcInTimeFrame;
      if constexpr (ModeId == ModeWithdEdxTrkQA || ModeId == ModeStandard) {
        bcTimeFrameId = UndefValueInt;
        bcBcInTimeFrame = UndefValueInt;
        if constexpr (ModeId == ModeWithdEdxTrkQA) {
          rowTPCTreeWithdEdxTrkQA.reserve(2 * v0s.size() + cascs.size());
        } else if (ModeId == ModeStandard) {
          rowTPCTree.reserve(2 * v0s.size() + cascs.size());
        }
      } else if constexpr (ModeId == ModeWithTrkQA) {
        bcTimeFrameId = bc.tfId();
        bcBcInTimeFrame = bc.bcInTF();
        rowTPCTreeWithTrkQA.reserve(2 * v0s.size() + cascs.size());
      }

      auto fillDaughterTrack = [&](const auto& mother, const TrksType::iterator& dauTrack, const V0Daughter& daughter, const aod::TracksQA& trackQAInstance, const bool existTrkQA) {
        const bool passTrackSelection = isTrackSelected(dauTrack, trackSelection);
        const bool passDownsamplig = downsampleTsalisCharged(fRndm, dauTrack.pt(), daughter.downsamplingTsalis, daughter.mass, sqrtSNN, daughter.maxPt4dwnsmplTsalis);
        const bool passNSigmaTofCut = std::fabs(daughter.tofNSigma) < daughter.nSigmaTofDauTrack || std::fabs(daughter.tofNSigma - NSigmaTofUnmatched) < NSigmaTofUnmatchedEqualityTolerance;
        const bool passMatchTofRequirement = !daughter.rejectNoTofDauTrack || std::fabs(daughter.tofNSigma - NSigmaTofUnmatched) > NSigmaTofUnmatchedEqualityTolerance;
        if (passTrackSelection && passDownsamplig && passNSigmaTofCut && passMatchTofRequirement) {
          fillSkimmedV0Table<IsCorrectedDeDx, ModeId>(mother, dauTrack, trackQAInstance, existTrkQA, collision, daughter.tpcNSigma, daughter.tofNSigma, daughter.tpcExpSignal, daughter.id, runnumber, daughter.dwnSmplFactor, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        }
      };

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

      /// Loop over v0 candidates
      for (const auto& v0 : v0s) {
        const auto v0Id = v0.v0addid();
        if (v0Id == MotherUndef) {
          continue;
        }
        const auto& posTrack = v0.posTrack_as<TrksType>();
        const auto& negTrack = v0.negTrack_as<TrksType>();

        const auto& [posTrackQA, existPosTrkQA] = getTrackQA(posTrack);
        const auto& [negTrackQA, existNegTrkQA] = getTrackQA(negTrack);

        const V0Mother v0Mother = createV0Mother(v0Id);
        const V0Daughter posDaughter = createV0Daughter<IsCorrectedDeDx>(v0, posTrack, v0Id, v0Mother.posDaughterId, true);
        const V0Daughter negDaughter = createV0Daughter<IsCorrectedDeDx>(v0, negTrack, v0Id, v0Mother.negDaughterId, false);

        fillDaughterTrack(v0, posTrack, posDaughter, posTrackQA, existPosTrkQA);
        fillDaughterTrack(v0, negTrack, negDaughter, negTrackQA, existNegTrkQA);
      }

      /// Loop over cascade candidates
      for (const auto& casc : cascs) {
        const auto cascId = casc.cascaddid();
        if (cascId == MotherUndef) {
          continue;
        }
        const auto& bachTrack = casc.bachelor_as<TrksType>();
        const V0Daughter bachDaughter = createV0Daughter<IsCorrectedDeDx>(casc, bachTrack, cascId, DaughterKaon);
        const auto& [bachTrackQA, existBachTrkQA] = getTrackQA(bachTrack);
        // Omega and antiomega
        fillDaughterTrack(casc, bachTrack, bachDaughter, bachTrackQA, existBachTrkQA);
      }
    }
  }

  /// Apply a track quality selection with a filter!
  void processStandard(Colls const& collisions, Trks const& myTracks, V0sWithID const& myV0s, CascsWithID const& myCascs, aod::BCsWithTimestamps const&)
  {
    runV0<false, ModeStandard, Trks, aod::BCsWithTimestamps>(collisions, myTracks, myV0s, myCascs, static_cast<TObject*>(nullptr));
  } /// process Standard
  PROCESS_SWITCH(TreeWriterTpcV0, processStandard, "Standard V0 Samples for PID", true);

  void processStandardWithCorrecteddEdx(Colls const& collisions, TrksWithDEdxCorrection const& myTracks, V0sWithID const& myV0s, CascsWithID const& myCascs, aod::BCsWithTimestamps const&)
  {
    runV0<true, ModeStandard, TrksWithDEdxCorrection, aod::BCsWithTimestamps>(collisions, myTracks, myV0s, myCascs, static_cast<TObject*>(nullptr));
  } /// process StandardWithCorrecteddEdx
  PROCESS_SWITCH(TreeWriterTpcV0, processStandardWithCorrecteddEdx, "Standard V0 Samples for PID with corrected dEdx", false);

  void processWithdEdxTrQA(Colls const& collisions, Trks const& myTracks, V0sWithID const& myV0s, CascsWithID const& myCascs, aod::BCsWithTimestamps const&, aod::TracksQAVersion const& tracksQA)
  {
    runV0<false, ModeWithdEdxTrkQA, Trks, aod::BCsWithTimestamps>(collisions, myTracks, myV0s, myCascs, tracksQA);
  } /// process with dEdx from TrackQA
  PROCESS_SWITCH(TreeWriterTpcV0, processWithdEdxTrQA, "Standard V0 Samples with dEdx from Track QA for PID", false);

  void processWithdEdxTrQAWithCorrecteddEdx(Colls const& collisions, TrksWithDEdxCorrection const& myTracks, V0sWithID const& myV0s, CascsWithID const& myCascs, aod::BCsWithTimestamps const&, aod::TracksQAVersion const& tracksQA)
  {
    runV0<true, ModeWithdEdxTrkQA, TrksWithDEdxCorrection, aod::BCsWithTimestamps>(collisions, myTracks, myV0s, myCascs, tracksQA);
  } /// process with dEdx from TrackQA
  PROCESS_SWITCH(TreeWriterTpcV0, processWithdEdxTrQAWithCorrecteddEdx, "Standard V0 Samples with dEdx from Track QA for PID with corrected dEdx", false);

  void processWithTrQA(Colls const& collisions, Trks const& myTracks, V0sWithID const& myV0s, CascsWithID const& myCascs, MyBCTable const&, aod::TracksQAVersion const& tracksQA)
  {
    runV0<false, ModeWithTrkQA, Trks, MyBCTable>(collisions, myTracks, myV0s, myCascs, tracksQA);
  } /// process with TrackQA
  PROCESS_SWITCH(TreeWriterTpcV0, processWithTrQA, "Standard V0 Samples with Track QA for PID", false);

  void processWithTrQAWithCorrecteddEdx(Colls const& collisions, TrksWithDEdxCorrection const& myTracks, V0sWithID const& myV0s, CascsWithID const& myCascs, MyBCTable const&, aod::TracksQAVersion const& tracksQA)
  {
    runV0<true, ModeWithTrkQA, TrksWithDEdxCorrection, MyBCTable>(collisions, myTracks, myV0s, myCascs, tracksQA);
  } /// process with TrackQA
  PROCESS_SWITCH(TreeWriterTpcV0, processWithTrQAWithCorrecteddEdx, "Standard V0 Samples with Track QA for PID with corrected dEdx", false);

  void processDummy(Colls const&) {}
  PROCESS_SWITCH(TreeWriterTpcV0, processDummy, "Dummy function", false);

}; /// struct TreeWriterTpcV0

struct TreeWriterTpcTof {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra,
                         aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa,
                         aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr,
                         aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa,
                         aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr,
                         aod::TrackSelection>;
  using TrksWithDEdxCorrection = soa::Join<aod::Tracks, aod::TracksExtra,
                                           aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa,
                                           aod::pidTPCFullPr, aod::pidTPCFullDe, aod::pidTPCFullTr,
                                           aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa,
                                           aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr,
                                           aod::TrackSelection, aod::DEdxsCorrected>;
  using Colls = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;
  using MyBCTable = soa::Join<aod::BCsWithTimestamps, aod::BCTFinfoTable>;

  /// Tables to be produced
  Produces<o2::aod::SkimmedTPCTOFTree> rowTPCTOFTree;
  Produces<o2::aod::SkimmedTPCTOFTreeWithdEdxTrkQA> rowTPCTOFTreeWithdEdxTrkQA;
  Produces<o2::aod::SkimmedTPCTOFTreeWithTrkQA> rowTPCTOFTreeWithTrkQA;

  constexpr static o2::track::PID::ID PidPion{o2::track::PID::Pion};
  constexpr static o2::track::PID::ID PidKaon{o2::track::PID::Kaon};
  constexpr static o2::track::PID::ID PidProton{o2::track::PID::Proton};
  constexpr static o2::track::PID::ID PidDeuteron{o2::track::PID::Deuteron};
  constexpr static o2::track::PID::ID PidTriton{o2::track::PID::Triton};

  /// Configurables
  Configurable<float> nClNorm{"nClNorm", 152., "Number of cluster normalization. Run 2: 159, Run 3 152"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};
  /// Triton
  Configurable<float> maxMomTPCOnlyTr{"maxMomTPCOnlyTr", 1.5, "Maximum momentum for TPC only cut triton"};
  Configurable<float> maxMomHardCutOnlyTr{"maxMomHardCutOnlyTr", 50, "Maximum TPC inner momentum for triton"};
  Configurable<float> nSigmaTPCOnlyTr{"nSigmaTPCOnlyTr", 4., "number of sigma for TPC only cut triton"};
  Configurable<float> nSigmaTpcTpctofTr{"nSigmaTpcTpctofTr", 4., "number of sigma for TPC cut for TPC and TOF combined triton"};
  Configurable<float> nSigmaTofTpctofTr{"nSigmaTofTpctofTr", 3., "number of sigma for TOF cut for TPC and TOF combined triton"};
  Configurable<double> dwnSmplFactorTr{"dwnSmplFactorTr", 1., "downsampling factor for triton, default fraction to keep is 1."};
  /// Deuteron
  Configurable<float> maxMomTPCOnlyDe{"maxMomTPCOnlyDe", 1.0, "Maximum momentum for TPC only cut deuteron"};
  Configurable<float> maxMomHardCutOnlyDe{"maxMomHardCutOnlyDe", 50, "Maximum TPC inner momentum for deuteron"};
  Configurable<float> nSigmaTPCOnlyDe{"nSigmaTPCOnlyDe", 4., "number of sigma for TPC only cut deuteron"};
  Configurable<float> nSigmaTpcTpctofDe{"nSigmaTpcTpctofDe", 4., "number of sigma for TPC cut for TPC and TOF combined deuteron"};
  Configurable<float> nSigmaTofTpctofDe{"nSigmaTofTpctofDe", 3., "number of sigma for TOF cut for TPC and TOF combined deuteron"};
  Configurable<double> dwnSmplFactorDe{"dwnSmplFactorDe", 1., "downsampling factor for deuteron, default fraction to keep is 1."};
  /// Proton
  Configurable<float> maxMomTPCOnlyPr{"maxMomTPCOnlyPr", 0.6, "Maximum momentum for TPC only cut proton"};
  Configurable<float> nSigmaTPCOnlyPr{"nSigmaTPCOnlyPr", 4., "number of sigma for TPC only cut proton"};
  Configurable<float> nSigmaTpcTpctofPr{"nSigmaTpcTpctofPr", 4., "number of sigma for TPC cut for TPC and TOF combined proton"};
  Configurable<float> nSigmaTofTpctofPr{"nSigmaTofTpctofPr", 3., "number of sigma for TOF cut for TPC and TOF combined proton"};
  Configurable<double> dwnSmplFactorPr{"dwnSmplFactorPr", 1., "downsampling factor for protons, default fraction to keep is 1."};
  /// Kaon
  Configurable<float> maxMomTPCOnlyKa{"maxMomTPCOnlyKa", 0.3, "Maximum momentum for TPC only cut kaon"};
  Configurable<float> maxMomHardCutOnlyKa{"maxMomHardCutOnlyKa", 50, "Maximum TPC inner momentum for kaons"};
  Configurable<float> nSigmaTPCOnlyKa{"nSigmaTPCOnlyKa", 4., "number of sigma for TPC only cut kaon"};
  Configurable<float> nSigmaTpcTpctofKa{"nSigmaTpcTpctofKa", 4., "number of sigma for TPC cut for TPC and TOF combined kaon"};
  Configurable<float> nSigmaTofTpctofKa{"nSigmaTofTpctofKa", 3., "number of sigma for TOF cut for TPC and TOF combined kaon"};
  Configurable<double> dwnSmplFactorKa{"dwnSmplFactorKa", 1., "downsampling factor for kaons, default fraction to keep is 1."};
  /// Pion
  Configurable<float> maxMomTPCOnlyPi{"maxMomTPCOnlyPi", 0.5, "Maximum momentum for TPC only cut pion"};
  Configurable<float> nSigmaTPCOnlyPi{"nSigmaTPCOnlyPi", 4., "number of sigma for TPC only cut pion"};
  Configurable<float> nSigmaTpcTpctofPi{"nSigmaTpcTpctofPi", 4., "number of sigma for TPC cut for TPC and TOF combined pion"};
  Configurable<float> nSigmaTofTpctofPi{"nSigmaTofTpctofPi", 4., "number of sigma for TOF cut for TPC and TOF combined pion"};
  Configurable<double> dwnSmplFactorPi{"dwnSmplFactorPi", 1., "downsampling factor for pions, default fraction to keep is 1."};
  /// pT dependent downsampling
  Configurable<float> sqrtSNN{"sqrtSNN", 5360., "sqrt(s_NN), used for downsampling with the Tsallis distribution"};
  Configurable<float> downsamplingTsalisTritons{"downsamplingTsalisTritons", -1., "Downsampling factor to reduce the number of tritons"};
  Configurable<float> downsamplingTsalisDeuterons{"downsamplingTsalisDeuterons", -1., "Downsampling factor to reduce the number of deuterons"};
  Configurable<float> downsamplingTsalisProtons{"downsamplingTsalisProtons", -1., "Downsampling factor to reduce the number of protons"};
  Configurable<float> downsamplingTsalisKaons{"downsamplingTsalisKaons", -1., "Downsampling factor to reduce the number of kaons"};
  Configurable<float> downsamplingTsalisPions{"downsamplingTsalisPions", -1., "Downsampling factor to reduce the number of pions"};

  Filter trackFilter = (trackSelection.node() == static_cast<int>(TrackSelectionNoCut)) ||
                       ((trackSelection.node() == static_cast<int>(TrackSelectionGlobalTrack)) && requireGlobalTrackInFilter()) ||
                       ((trackSelection.node() == static_cast<int>(TrackSelectionTrackWoPtEta)) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((trackSelection.node() == static_cast<int>(TrackSelectionGlobalTrackWoDCA)) && requireGlobalTrackWoDCAInFilter()) ||
                       ((trackSelection.node() == static_cast<int>(TrackSelectionQualityTracks)) && requireQualityTracksInFilter()) ||
                       ((trackSelection.node() == static_cast<int>(TrackSelectionInAcceptanceTracks)) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));

  ctpRateFetcher mRateFetcher;

  struct TofTrack {
    bool isApplyHardCutOnly;
    double maxMomHardCutOnly;
    double maxMomTPCOnly;
    double tpcNSigma;
    double nSigmaTPCOnly;
    double downsamplingTsalis;
    double mass;
    double tofNSigma;
    double itsNSigma;
    double tpcExpSignal;
    o2::track::PID::ID pid;
    double dwnSmplFactor;
    double nSigmaTofTpctof;
    double nSigmaTpcTpctof;
  };

  TRandom3* fRndm = new TRandom3(0);

  template <bool DoCorrectDeDx, int ModeId, typename T, typename C>
  void fillSkimmedTpcTofTable(T const& track, aod::TracksQA const& trackQA, const bool existTrkQA, C const& collision, const float nSigmaTPC, const float nSigmaTOF, const float nSigmaITS, const float dEdxExp, const o2::track::PID::ID id, const int runnumber, const double dwnSmplFactor, const double hadronicRate, const int bcGlobalIndex, const int bcTimeFrameId, const int bcBcInTimeFrame)
  {
    const double ncl = track.tpcNClsFound();
    const double nclPID = track.tpcNClsFindableMinusPID();
    const double p = track.tpcInnerParam();
    const double mass = o2::track::pid_constants::sMasses[id];
    const double bg = p / mass;
    const int multTPC = collision.multTPC();
    const auto trackOcc = collision.trackOccupancyInTimeRange();
    const auto ft0Occ = collision.ft0cOccupancyInTimeRange();

    const double pseudoRndm = track.pt() * 1000. - static_cast<int64_t>(track.pt() * 1000);
    if (pseudoRndm < dwnSmplFactor) {
      float usedEdx;
      if constexpr (DoCorrectDeDx) {
        usedEdx = track.tpcSignalCorrected();
      } else {
        usedEdx = track.tpcSignal();
      }
      if (ModeId == ModeStandard) {
        rowTPCTOFTree(usedEdx,
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
                      std::sqrt(nClNorm / ncl),
                      nclPID,
                      id,
                      nSigmaTPC,
                      nSigmaTOF,
                      runnumber,
                      trackOcc,
                      ft0Occ,
                      hadronicRate,
                      nSigmaITS);
      } else if constexpr (ModeId == ModeWithdEdxTrkQA) {
        rowTPCTOFTreeWithdEdxTrkQA(usedEdx,
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
                                   std::sqrt(nClNorm / ncl),
                                   nclPID,
                                   id,
                                   nSigmaTPC,
                                   nSigmaTOF,
                                   runnumber,
                                   trackOcc,
                                   ft0Occ,
                                   hadronicRate,
                                   nSigmaITS,
                                   existTrkQA ? trackQA.tpcdEdxNorm() : UndefValueFloat);
      } else if constexpr (ModeId == ModeWithTrkQA) {
        rowTPCTOFTreeWithTrkQA(usedEdx,
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
                               std::sqrt(nClNorm / ncl),
                               nclPID,
                               id,
                               nSigmaTPC,
                               nSigmaTOF,
                               runnumber,
                               trackOcc,
                               ft0Occ,
                               hadronicRate,
                               nSigmaITS,
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
                               existTrkQA ? trackQA.tpcdEdxNorm() : UndefValueFloat);
      }
    }
  }

  void init(o2::framework::InitContext&)
  {
    const std::array<bool, 7> doprocess{doprocessStandard, doprocessStandardWithCorrecteddEdx, doprocessWithdEdxTrQA, doprocessWithdEdxTrQAWithCorrecteddEdx, doprocessWithTrQA, doprocessWithTrQAWithCorrecteddEdx, doprocessDummy};
    if (std::accumulate(doprocess.begin(), doprocess.end(), 0) != 1) {
      LOGP(fatal, "One and only one process function should be enabled");
    }

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);
  }

  Preslice<Trks> perCollisionTracks = aod::track::collisionId;
  Preslice<TrksWithDEdxCorrection> perCollisionTracksWithCorrecteddEdx = aod::track::collisionId;

  template <bool IsCorrectedDeDx, int ModeId, typename TrksType, typename BCType, typename TrkQAType>
  void runTof(Colls const& collisions, TrksType const& myTracks, TrkQAType const& tracksQA, Preslice<TrksType> const& perCollisionTracksType)
  {
    constexpr bool IsWithTrackQa = ModeId != ModeStandard;

    std::vector<int64_t> labelTrack2TrackQA;
    if constexpr (IsWithTrackQa) {
      labelTrack2TrackQA.resize(myTracks.size(), -1);
      for (const auto& trackQA : tracksQA) {
        const int64_t trackId = trackQA.trackId();
        labelTrack2TrackQA.at(trackId) = trackQA.globalIndex();
      }
    }
    for (const auto& collision : collisions) {
      /// Check event selection
      const auto& tracks = myTracks.sliceBy(perCollisionTracksType, collision.globalIndex());
      if (!isEventSelected(collision, applyEvSel)) {
        continue;
      }
      const auto& tracksWithITSPid = soa::Attach<TrksType,
                                                 aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi,
                                                 aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe,
                                                 aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe, aod::pidits::ITSNSigmaAl>(tracks);

      const auto& bc = collision.bc_as<BCType>();
      const int runnumber = bc.runNumber();
      const float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * OneToKilo;
      const int bcGlobalIndex = bc.globalIndex();
      int bcTimeFrameId, bcBcInTimeFrame;
      if constexpr (ModeId == ModeStandard || ModeId == ModeWithdEdxTrkQA) {
        bcTimeFrameId = UndefValueInt;
        bcBcInTimeFrame = UndefValueInt;
        if constexpr (ModeId == ModeStandard) {
          rowTPCTOFTree.reserve(tracks.size());
        } else {
          rowTPCTOFTreeWithdEdxTrkQA.reserve(tracks.size());
        }
      } else {
        bcTimeFrameId = bc.tfId();
        bcBcInTimeFrame = bc.bcInTF();
        rowTPCTOFTreeWithTrkQA.reserve(tracks.size());
      }
      for (auto const& trk : tracksWithITSPid) {
        if (!isTrackSelected(trk, trackSelection)) {
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

        TofTrack tofTriton(true, maxMomHardCutOnlyTr, maxMomTPCOnlyTr, trk.tpcNSigmaTr(), nSigmaTPCOnlyTr, downsamplingTsalisTritons, MassTriton, trk.tofNSigmaTr(), trk.itsNSigmaTr(), trk.tpcExpSignalTr(tpcSignalGeneric<IsCorrectedDeDx>(trk)), PidTriton, dwnSmplFactorTr, nSigmaTofTpctofTr, nSigmaTpcTpctofTr);

        TofTrack tofDeuteron(true, maxMomHardCutOnlyDe, maxMomTPCOnlyDe, trk.tpcNSigmaDe(), nSigmaTPCOnlyDe, downsamplingTsalisDeuterons, MassDeuteron, trk.tofNSigmaDe(), trk.itsNSigmaDe(), trk.tpcExpSignalDe(tpcSignalGeneric<IsCorrectedDeDx>(trk)), PidDeuteron, dwnSmplFactorDe, nSigmaTofTpctofDe, nSigmaTpcTpctofDe);

        TofTrack tofProton(false, UndefValueDouble, maxMomTPCOnlyPr, trk.tpcNSigmaPr(), nSigmaTPCOnlyPr, downsamplingTsalisProtons, MassProton, trk.tofNSigmaPr(), trk.itsNSigmaPr(), trk.tpcExpSignalPr(tpcSignalGeneric<IsCorrectedDeDx>(trk)), PidProton, dwnSmplFactorPr, nSigmaTofTpctofPr, nSigmaTpcTpctofPr);

        TofTrack tofKaon(true, maxMomHardCutOnlyKa, maxMomTPCOnlyKa, trk.tpcNSigmaKa(), nSigmaTPCOnlyKa, downsamplingTsalisKaons, MassKPlus, trk.tofNSigmaKa(), trk.itsNSigmaKa(), trk.tpcExpSignalKa(tpcSignalGeneric<IsCorrectedDeDx>(trk)), PidKaon, dwnSmplFactorKa, nSigmaTofTpctofKa, nSigmaTpcTpctofKa);

        TofTrack tofPion(false, UndefValueDouble, maxMomTPCOnlyPi, trk.tpcNSigmaPi(), nSigmaTPCOnlyPi, downsamplingTsalisPions, MassPiPlus, trk.tofNSigmaPi(), trk.itsNSigmaPi(), trk.tpcExpSignalPi(tpcSignalGeneric<IsCorrectedDeDx>(trk)), PidPion, dwnSmplFactorPi, nSigmaTofTpctofPi, nSigmaTpcTpctofPi);

        for (const auto& tofTrack : {&tofTriton, &tofDeuteron, &tofProton, &tofKaon, &tofPion}) {
          if ((!tofTrack->isApplyHardCutOnly || trk.tpcInnerParam() < tofTrack->maxMomHardCutOnly) &&
              ((trk.tpcInnerParam() <= tofTrack->maxMomTPCOnly && std::fabs(tofTrack->tpcNSigma) < tofTrack->nSigmaTPCOnly) ||
               (trk.tpcInnerParam() > tofTrack->maxMomTPCOnly && std::fabs(tofTrack->tofNSigma) < tofTrack->nSigmaTofTpctof && std::fabs(tofTrack->tpcNSigma) < tofTrack->nSigmaTpcTpctof)) &&
              downsampleTsalisCharged(fRndm, trk.pt(), tofTrack->downsamplingTsalis, tofTrack->mass, sqrtSNN)) {
            fillSkimmedTpcTofTable<IsCorrectedDeDx, ModeId>(trk, trackQA, existTrkQA, collision, tofTrack->tpcNSigma, tofTrack->tofNSigma, tofTrack->itsNSigma, tofTrack->tpcExpSignal, tofTrack->pid, runnumber, tofTrack->dwnSmplFactor, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
          }
        }
      } /// Loop tracks
    }
  }

  void processStandard(Colls const& collisions, Trks const& myTracks, aod::BCsWithTimestamps const&)
  {
    runTof<false, ModeStandard, Trks, aod::BCsWithTimestamps>(collisions, myTracks, static_cast<TObject*>(nullptr), perCollisionTracks);
  } /// process
  PROCESS_SWITCH(TreeWriterTpcTof, processStandard, "Standard Samples for PID", true);

  void processStandardWithCorrecteddEdx(Colls const& collisions, TrksWithDEdxCorrection const& myTracks, aod::BCsWithTimestamps const&)
  {
    runTof<true, ModeStandard, TrksWithDEdxCorrection, aod::BCsWithTimestamps>(collisions, myTracks, static_cast<TObject*>(nullptr), perCollisionTracksWithCorrecteddEdx);
  } /// process
  PROCESS_SWITCH(TreeWriterTpcTof, processStandardWithCorrecteddEdx, "Standard Samples for PID with corrected dEdx", false);

  void processWithdEdxTrQA(Colls const& collisions, Trks const& myTracks, aod::BCsWithTimestamps const&, aod::TracksQAVersion const& tracksQA)
  {
    runTof<false, ModeWithdEdxTrkQA, Trks, aod::BCsWithTimestamps>(collisions, myTracks, tracksQA, perCollisionTracks);
  } /// process
  PROCESS_SWITCH(TreeWriterTpcTof, processWithdEdxTrQA, "Samples for PID with TrackQA info", false);

  void processWithdEdxTrQAWithCorrecteddEdx(Colls const& collisions, TrksWithDEdxCorrection const& myTracks, aod::BCsWithTimestamps const&, aod::TracksQAVersion const& tracksQA)
  {
    runTof<true, ModeWithdEdxTrkQA, TrksWithDEdxCorrection, aod::BCsWithTimestamps>(collisions, myTracks, tracksQA, perCollisionTracksWithCorrecteddEdx);
  } /// process
  PROCESS_SWITCH(TreeWriterTpcTof, processWithdEdxTrQAWithCorrecteddEdx, "Samples for PID with TrackQA info with corrected dEdx", false);

  void processWithTrQA(Colls const& collisions, Trks const& myTracks, MyBCTable const&, aod::TracksQAVersion const& tracksQA)
  {
    runTof<false, ModeWithTrkQA, Trks, MyBCTable>(collisions, myTracks, tracksQA, perCollisionTracks);
  } /// process
  PROCESS_SWITCH(TreeWriterTpcTof, processWithTrQA, "Samples for PID with TrackQA info", false);

  void processWithTrQAWithCorrecteddEdx(Colls const& collisions, TrksWithDEdxCorrection const& myTracks, MyBCTable const&, aod::TracksQAVersion const& tracksQA)
  {
    runTof<true, ModeWithTrkQA, TrksWithDEdxCorrection, MyBCTable>(collisions, myTracks, tracksQA, perCollisionTracksWithCorrecteddEdx);
  } /// process
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
