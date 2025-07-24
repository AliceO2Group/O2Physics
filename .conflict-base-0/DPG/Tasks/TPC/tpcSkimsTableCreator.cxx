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

#include "tpcSkimsTableCreator.h"

#include <CCDB/BasicCCDBManager.h>

#include <cmath>
#include <string>
#include <vector>
/// ROOT
#include "TRandom3.h"
/// O2
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
/// O2Physics
#include "PWGDQ/DataModel/ReducedInfoTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/OccupancyTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTPCBase.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::dataformats;

struct TreeWriterTpcV0 {

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  using Trks = soa::Join<aod::Tracks, aod::V0Bits, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TrackSelection>;
  using TrksWithDEdxCorrection = soa::Join<aod::Tracks, aod::V0Bits, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TrackSelection, aod::DEdxsCorrected>;
  using Colls = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;
  using MyBCTable = soa::Join<aod::BCsWithTimestamps, aod::BCTFinfoTable>;
  using V0sWithID = soa::Join<aod::V0Datas, aod::V0MapID>;

  /// Tables to be produced
  Produces<o2::aod::SkimmedTPCV0Tree> rowTPCTree;
  Produces<o2::aod::SkimmedTPCV0TreeWithdEdxTrkQA> rowTPCTreeWithdEdxTrkQA;
  Produces<o2::aod::SkimmedTPCV0TreeWithTrkQA> rowTPCTreeWithTrkQA;

  /// Configurables
  Configurable<float> nSigmaTOFdautrack{"nSigmaTOFdautrack", 999., "n-sigma TOF cut on the proton daughter tracks. Set 999 to switch it off."};
  Configurable<float> nClNorm{"nClNorm", 152., "Number of cluster normalization. Run 2: 159, Run 3 152"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};
  /// Configurables downsampling
  Configurable<double> dwnSmplFactor_Pi{"dwnSmplFactor_Pi", 1., "downsampling factor for pions, default fraction to keep is 1."};
  Configurable<double> dwnSmplFactor_Pr{"dwnSmplFactor_Pr", 1., "downsampling factor for protons, default fraction to keep is 1."};
  Configurable<double> dwnSmplFactor_El{"dwnSmplFactor_El", 1., "downsampling factor for electrons, default fraction to keep is 1."};
  Configurable<float> sqrtSNN{"sqrt_s_NN", 0., "sqrt(s_NN), used for downsampling with the Tsallis distribution"};
  Configurable<float> downsamplingTsalisPions{"downsamplingTsalisPions", -1., "Downsampling factor to reduce the number of pions"};
  Configurable<float> downsamplingTsalisProtons{"downsamplingTsalisProtons", -1., "Downsampling factor to reduce the number of protons"};
  Configurable<float> downsamplingTsalisElectrons{"downsamplingTsalisElectrons", -1., "Downsampling factor to reduce the number of electrons"};
  Configurable<float> maxPt4dwnsmplTsalisPions{"maxPt4dwnsmplTsalisPions", 100., "Maximum Pt for applying downsampling factor of pions"};
  Configurable<float> maxPt4dwnsmplTsalisProtons{"maxPt4dwnsmplTsalisProtons", 100., "Maximum Pt for applying downsampling factor of protons"};
  Configurable<float> maxPt4dwnsmplTsalisElectrons{"maxPt4dwnsmplTsalisElectrons", 100., "Maximum Pt for applying  downsampling factor of electrons"};

  Filter trackFilter = (trackSelection.node() == 0) ||
                       ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                       ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                       ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                       ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));

  ctpRateFetcher mRateFetcher;

  /// Funktion to fill skimmed tables
  template <bool doUseCorreceddEdx = false, typename T, typename C, typename V0>
  void fillSkimmedV0Table(V0 const& v0, T const& track, C const& collision, const float nSigmaTPC, const float nSigmaTOF, const float dEdxExp, const o2::track::PID::ID id, int runnumber, double dwnSmplFactor, float hadronicRate)
  {

    const double ncl = track.tpcNClsFound();
    const double nclPID = track.tpcNClsFindableMinusPID();
    const double p = track.tpcInnerParam();
    const double mass = o2::track::pid_constants::sMasses[id];
    const double bg = p / mass;
    const int multTPC = collision.multTPC();
    auto trackocc = collision.trackOccupancyInTimeRange();
    auto ft0occ = collision.ft0cOccupancyInTimeRange();

    const float alpha = v0.alpha();
    const float qt = v0.qtarm();
    const float cosPA = v0.v0cosPA();
    const float pT = v0.pt();
    const float v0radius = v0.v0radius();
    const float gammapsipair = v0.psipair();

    const double pseudoRndm = track.pt() * 1000. - static_cast<int64_t>(track.pt() * 1000);
    if (pseudoRndm < dwnSmplFactor) {
      float usedDedx;
      if constexpr (doUseCorreceddEdx) {
        usedDedx = track.tpcSignalCorrected();
      } else {
        usedDedx = track.tpcSignal();
      }
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
                 multTPC / 11000.,
                 std::sqrt(nClNorm / ncl),
                 nclPID,
                 id,
                 nSigmaTPC,
                 nSigmaTOF,
                 alpha,
                 qt,
                 cosPA,
                 pT,
                 v0radius,
                 gammapsipair,
                 runnumber,
                 trackocc,
                 ft0occ,
                 hadronicRate);
    }
  };

  template <bool doUseCorreceddEdx = false, typename T, typename TQA, typename C, typename V0>
  void fillSkimmedV0TableWithdEdxTrQA(V0 const& v0, T const& track, TQA const& trackQA, bool existTrkQA, C const& collision, const float nSigmaTPC, const float nSigmaTOF, const float dEdxExp, const o2::track::PID::ID id, int runnumber, double dwnSmplFactor, float hadronicRate)
  {

    const double ncl = track.tpcNClsFound();
    const double nclPID = track.tpcNClsFindableMinusPID();
    const double p = track.tpcInnerParam();
    const double mass = o2::track::pid_constants::sMasses[id];
    const double bg = p / mass;
    const int multTPC = collision.multTPC();
    auto trackocc = collision.trackOccupancyInTimeRange();
    auto ft0occ = collision.ft0cOccupancyInTimeRange();

    const float alpha = v0.alpha();
    const float qt = v0.qtarm();
    const float cosPA = v0.v0cosPA();
    const float pT = v0.pt();
    const float v0radius = v0.v0radius();
    const float gammapsipair = v0.psipair();

    const double pseudoRndm = track.pt() * 1000. - static_cast<int64_t>(track.pt() * 1000);
    if (pseudoRndm < dwnSmplFactor) {
      float usedDedx;
      if constexpr (doUseCorreceddEdx) {
        usedDedx = track.tpcSignalCorrected();
      } else {
        usedDedx = track.tpcSignal();
      }
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
                              multTPC / 11000.,
                              std::sqrt(nClNorm / ncl),
                              nclPID,
                              id,
                              nSigmaTPC,
                              nSigmaTOF,
                              alpha,
                              qt,
                              cosPA,
                              pT,
                              v0radius,
                              gammapsipair,
                              runnumber,
                              trackocc,
                              ft0occ,
                              hadronicRate,
                              existTrkQA ? trackQA.tpcdEdxNorm() : -999);
    }
  };

  /// Function to fill skimmed tables
  template <bool doUseCorreceddEdx = false, typename T, typename TQA, typename C, typename V0>
  void fillSkimmedV0TableWithTrQA(V0 const& v0, T const& track, TQA const& trackQA, bool existTrkQA, C const& collision, const float nSigmaTPC, const float nSigmaTOF, const float dEdxExp, const o2::track::PID::ID id, int runnumber, double dwnSmplFactor, float hadronicRate, int bcGlobalIndex, int bcTimeFrameId, int bcBcInTimeFrame)
  {

    const double ncl = track.tpcNClsFound();
    const double nclPID = track.tpcNClsFindableMinusPID();
    const double p = track.tpcInnerParam();
    const double mass = o2::track::pid_constants::sMasses[id];
    const double bg = p / mass;
    const int multTPC = collision.multTPC();
    auto trackocc = collision.trackOccupancyInTimeRange();
    auto ft0occ = collision.ft0cOccupancyInTimeRange();

    const float alpha = v0.alpha();
    const float qt = v0.qtarm();
    const float cosPA = v0.v0cosPA();
    const float pT = v0.pt();
    const float v0radius = v0.v0radius();
    const float gammapsipair = v0.psipair();

    const double pseudoRndm = track.pt() * 1000. - static_cast<int64_t>(track.pt() * 1000);
    if (pseudoRndm < dwnSmplFactor) {
      float usedDedx;
      if constexpr (doUseCorreceddEdx) {
        usedDedx = track.tpcSignalCorrected();
      } else {
        usedDedx = track.tpcSignal();
      }
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
                          multTPC / 11000.,
                          std::sqrt(nClNorm / ncl),
                          nclPID,
                          id,
                          nSigmaTPC,
                          nSigmaTOF,
                          alpha,
                          qt,
                          cosPA,
                          pT,
                          v0radius,
                          gammapsipair,
                          runnumber,
                          trackocc,
                          ft0occ,
                          hadronicRate,
                          bcGlobalIndex,
                          bcTimeFrameId,
                          bcBcInTimeFrame,
                          existTrkQA ? trackQA.tpcClusterByteMask() : -999,
                          existTrkQA ? trackQA.tpcdEdxMax0R() : -999,
                          existTrkQA ? trackQA.tpcdEdxMax1R() : -999,
                          existTrkQA ? trackQA.tpcdEdxMax2R() : -999,
                          existTrkQA ? trackQA.tpcdEdxMax3R() : -999,
                          existTrkQA ? trackQA.tpcdEdxTot0R() : -999,
                          existTrkQA ? trackQA.tpcdEdxTot1R() : -999,
                          existTrkQA ? trackQA.tpcdEdxTot2R() : -999,
                          existTrkQA ? trackQA.tpcdEdxTot3R() : -999,
                          existTrkQA ? trackQA.tpcdEdxNorm() : -999);
    }
  };

  double tsalisCharged(double pt, double mass, double sqrts)
  {
    const double a = 6.81, b = 59.24;
    const double c = 0.082, d = 0.151;
    const double mt = std::sqrt(mass * mass + pt * pt);
    const double n = a + b / sqrts;
    const double T = c + d / sqrts;
    const double p0 = n * T;
    const double result = std::pow((1. + mt / p0), -n);
    return result;
  };

  /// Random downsampling trigger function using Tsalis/Hagedorn spectra fit (sqrt(s) = 62.4 GeV to 13 TeV)
  /// as in https://iopscience.iop.org/article/10.1088/2399-6528/aab00f/pdf
  TRandom3* fRndm = new TRandom3(0);
  bool downsampleTsalisCharged(double pt, double factor1Pt, double sqrts, double mass, double maxPt)
  {
    if (factor1Pt < 0.) {
      return true;
    }
    if (pt > maxPt) {
      return true;
    }
    const double prob = tsalisCharged(pt, mass, sqrts) * pt;
    const double probNorm = tsalisCharged(1., mass, sqrts);
    if ((fRndm->Rndm() * ((prob / probNorm) * pt * pt)) > factor1Pt) {
      return false;
    } else {
      return true;
    }
  };

  /// Event selection
  template <typename CollisionType, typename TrackType>
  bool isEventSelected(const CollisionType& collision, const TrackType& /*tracks*/)
  {
    if (applyEvSel == 1) {
      if (!collision.sel7()) {
        return false;
      }
    } else if (applyEvSel == 2) {
      if (!collision.sel8()) {
        return false;
      }
    }
    return true;
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);
  }

  /// Apply a track quality selection with a filter!
  void processStandard(Colls::iterator const& collision, soa::Filtered<Trks> const& tracks, V0sWithID const& v0s, aod::BCsWithTimestamps const&)
  {
    /// Check event slection
    if (!isEventSelected(collision, tracks)) {
      return;
    }
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    const int runnumber = bc.runNumber();
    float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * 1.e-3;

    rowTPCTree.reserve(tracks.size());

    /// Loop over v0 candidates
    for (const auto& v0 : v0s) {
      auto posTrack = v0.posTrack_as<soa::Filtered<Trks>>();
      auto negTrack = v0.negTrack_as<soa::Filtered<Trks>>();
      if (v0.v0addid() == -1) {
        continue;
      }
      // gamma
      if (static_cast<bool>(posTrack.pidbit() & (1 << 0)) && static_cast<bool>(negTrack.pidbit() & (1 << 0))) {
        if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisElectrons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Electron], maxPt4dwnsmplTsalisElectrons)) {
          fillSkimmedV0Table(v0, posTrack, collision, posTrack.tpcNSigmaEl(), posTrack.tofNSigmaEl(), posTrack.tpcExpSignalEl(posTrack.tpcSignal()), o2::track::PID::Electron, runnumber, dwnSmplFactor_El, hadronicRate);
        }
        if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisElectrons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Electron], maxPt4dwnsmplTsalisElectrons)) {
          fillSkimmedV0Table(v0, negTrack, collision, negTrack.tpcNSigmaEl(), negTrack.tofNSigmaEl(), negTrack.tpcExpSignalEl(negTrack.tpcSignal()), o2::track::PID::Electron, runnumber, dwnSmplFactor_El, hadronicRate);
        }
      }
      // Ks0
      if (static_cast<bool>(posTrack.pidbit() & (1 << 1)) && static_cast<bool>(negTrack.pidbit() & (1 << 1))) {
        if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
          fillSkimmedV0Table(v0, posTrack, collision, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSignalPi(posTrack.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
        }
        if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
          fillSkimmedV0Table(v0, negTrack, collision, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSignalPi(negTrack.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
        }
      }
      // Lambda
      if (static_cast<bool>(posTrack.pidbit() & (1 << 2)) && static_cast<bool>(negTrack.pidbit() & (1 << 2))) {
        if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton], maxPt4dwnsmplTsalisProtons)) {
          if (std::abs(posTrack.tofNSigmaPr()) <= nSigmaTOFdautrack) {
            fillSkimmedV0Table(v0, posTrack, collision, posTrack.tpcNSigmaPr(), posTrack.tofNSigmaPr(), posTrack.tpcExpSignalPr(posTrack.tpcSignal()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
          }
        }
        if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
          fillSkimmedV0Table(v0, negTrack, collision, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSignalPi(negTrack.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
        }
      }
      // Antilambda
      if (static_cast<bool>(posTrack.pidbit() & (1 << 3)) && static_cast<bool>(negTrack.pidbit() & (1 << 3))) {
        if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
          fillSkimmedV0Table(v0, posTrack, collision, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSignalPi(posTrack.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
        }
        if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton], maxPt4dwnsmplTsalisProtons)) {
          if (std::abs(negTrack.tofNSigmaPr()) <= nSigmaTOFdautrack) {
            fillSkimmedV0Table(v0, negTrack, collision, negTrack.tpcNSigmaPr(), negTrack.tofNSigmaPr(), negTrack.tpcExpSignalPr(negTrack.tpcSignal()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
          }
        }
      }
    }
  } /// process Standard
  PROCESS_SWITCH(TreeWriterTpcV0, processStandard, "Standard V0 Samples for PID", true);

  void processStandardWithCorrecteddEdx(Colls::iterator const& collision, soa::Filtered<TrksWithDEdxCorrection> const& tracks, V0sWithID const& v0s, aod::BCsWithTimestamps const&)
  {
    /// Check event slection
    if (!isEventSelected(collision, tracks)) {
      return;
    }
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    const int runnumber = bc.runNumber();
    float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * 1.e-3;

    rowTPCTree.reserve(tracks.size());

    /// Loop over v0 candidates
    for (const auto& v0 : v0s) {
      auto posTrack = v0.posTrack_as<soa::Filtered<TrksWithDEdxCorrection>>();
      auto negTrack = v0.negTrack_as<soa::Filtered<TrksWithDEdxCorrection>>();
      if (v0.v0addid() == -1) {
        continue;
      }
      // gamma
      if (static_cast<bool>(posTrack.pidbit() & (1 << 0)) && static_cast<bool>(negTrack.pidbit() & (1 << 0))) {
        if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisElectrons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Electron], maxPt4dwnsmplTsalisElectrons)) {
          fillSkimmedV0Table<true>(v0, posTrack, collision, posTrack.tpcNSigmaEl(), posTrack.tofNSigmaEl(), posTrack.tpcExpSignalEl(posTrack.tpcSignalCorrected()), o2::track::PID::Electron, runnumber, dwnSmplFactor_El, hadronicRate);
        }
        if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisElectrons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Electron], maxPt4dwnsmplTsalisElectrons)) {
          fillSkimmedV0Table<true>(v0, negTrack, collision, negTrack.tpcNSigmaEl(), negTrack.tofNSigmaEl(), negTrack.tpcExpSignalEl(negTrack.tpcSignalCorrected()), o2::track::PID::Electron, runnumber, dwnSmplFactor_El, hadronicRate);
        }
      }
      // Ks0
      if (static_cast<bool>(posTrack.pidbit() & (1 << 1)) && static_cast<bool>(negTrack.pidbit() & (1 << 1))) {
        if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
          fillSkimmedV0Table<true>(v0, posTrack, collision, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSignalPi(posTrack.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
        }
        if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
          fillSkimmedV0Table<true>(v0, negTrack, collision, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSignalPi(negTrack.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
        }
      }
      // Lambda
      if (static_cast<bool>(posTrack.pidbit() & (1 << 2)) && static_cast<bool>(negTrack.pidbit() & (1 << 2))) {
        if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton], maxPt4dwnsmplTsalisProtons)) {
          if (std::abs(posTrack.tofNSigmaPr()) <= nSigmaTOFdautrack) {
            fillSkimmedV0Table<true>(v0, posTrack, collision, posTrack.tpcNSigmaPr(), posTrack.tofNSigmaPr(), posTrack.tpcExpSignalPr(posTrack.tpcSignalCorrected()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
          }
        }
        if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
          fillSkimmedV0Table<true>(v0, negTrack, collision, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSignalPi(negTrack.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
        }
      }
      // Antilambda
      if (static_cast<bool>(posTrack.pidbit() & (1 << 3)) && static_cast<bool>(negTrack.pidbit() & (1 << 3))) {
        if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
          fillSkimmedV0Table<true>(v0, posTrack, collision, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSignalPi(posTrack.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
        }
        if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton], maxPt4dwnsmplTsalisProtons)) {
          if (std::abs(negTrack.tofNSigmaPr()) <= nSigmaTOFdautrack) {
            fillSkimmedV0Table<true>(v0, negTrack, collision, negTrack.tpcNSigmaPr(), negTrack.tofNSigmaPr(), negTrack.tpcExpSignalPr(negTrack.tpcSignalCorrected()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
          }
        }
      }
    }
  } /// process Standard
  PROCESS_SWITCH(TreeWriterTpcV0, processStandardWithCorrecteddEdx, "Standard V0 Samples for PID with corrected dEdx", false);

  Preslice<Trks> perCollisionTracks = aod::track::collisionId;
  Preslice<V0sWithID> perCollisionV0s = aod::v0data::collisionId;
  void processWithdEdxTrQA(Colls const& collisions, Trks const& myTracks, V0sWithID const& myV0s, aod::BCsWithTimestamps const&, aod::TracksQAVersion const& tracksQA)
  {
    std::vector<int64_t> labelTrack2TrackQA;
    labelTrack2TrackQA.clear();
    labelTrack2TrackQA.resize(myTracks.size(), -1);
    for (const auto& trackQA : tracksQA) {
      int64_t trackId = trackQA.trackId();
      int64_t trackQAIndex = trackQA.globalIndex();
      labelTrack2TrackQA[trackId] = trackQAIndex;
    }
    for (const auto& collision : collisions) {
      auto tracks = myTracks.sliceBy(perCollisionTracks, collision.globalIndex());
      auto v0s = myV0s.sliceBy(perCollisionV0s, collision.globalIndex());
      /// Check event slection
      if (!isEventSelected(collision, tracks)) {
        continue;
      }
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      const int runnumber = bc.runNumber();
      float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * 1.e-3;
      rowTPCTreeWithTrkQA.reserve(tracks.size());
      /// Loop over v0 candidates
      for (const auto& v0 : v0s) {
        auto posTrack = v0.posTrack_as<Trks>();
        auto negTrack = v0.negTrack_as<Trks>();
        if (v0.v0addid() == -1) {
          continue;
        }
        aod::TracksQA posTrackQA;
        aod::TracksQA negTrackQA;
        bool existPosTrkQA;
        bool existNegTrkQA;
        if (labelTrack2TrackQA[posTrack.globalIndex()] != -1) {
          posTrackQA = tracksQA.iteratorAt(labelTrack2TrackQA[posTrack.globalIndex()]);
          existPosTrkQA = true;
        } else {
          posTrackQA = tracksQA.iteratorAt(0);
          existPosTrkQA = false;
        }
        if (labelTrack2TrackQA[negTrack.globalIndex()] != -1) {
          negTrackQA = tracksQA.iteratorAt(labelTrack2TrackQA[negTrack.globalIndex()]);
          existNegTrkQA = true;
        } else {
          negTrackQA = tracksQA.iteratorAt(0);
          existNegTrkQA = false;
        }

        // gamma
        if (static_cast<bool>(posTrack.pidbit() & (1 << 0)) && static_cast<bool>(negTrack.pidbit() & (1 << 0))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisElectrons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Electron], maxPt4dwnsmplTsalisElectrons)) {
            fillSkimmedV0TableWithdEdxTrQA(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaEl(), posTrack.tofNSigmaEl(), posTrack.tpcExpSignalEl(posTrack.tpcSignal()), o2::track::PID::Electron, runnumber, dwnSmplFactor_El, hadronicRate);
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisElectrons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Electron], maxPt4dwnsmplTsalisElectrons)) {
            fillSkimmedV0TableWithdEdxTrQA(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaEl(), negTrack.tofNSigmaEl(), negTrack.tpcExpSignalEl(negTrack.tpcSignal()), o2::track::PID::Electron, runnumber, dwnSmplFactor_El, hadronicRate);
          }
        }
        // Ks0
        if (static_cast<bool>(posTrack.pidbit() & (1 << 1)) && static_cast<bool>(negTrack.pidbit() & (1 << 1))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithdEdxTrQA(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSignalPi(posTrack.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithdEdxTrQA(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSignalPi(negTrack.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
          }
        }
        // Lambda
        if (static_cast<bool>(posTrack.pidbit() & (1 << 2)) && static_cast<bool>(negTrack.pidbit() & (1 << 2))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton], maxPt4dwnsmplTsalisProtons)) {
            if (std::abs(posTrack.tofNSigmaPr()) <= nSigmaTOFdautrack) {
              fillSkimmedV0TableWithdEdxTrQA(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaPr(), posTrack.tofNSigmaPr(), posTrack.tpcExpSignalPr(posTrack.tpcSignal()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
            }
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithdEdxTrQA(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSignalPi(negTrack.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
          }
        }
        // Antilambda
        if (static_cast<bool>(posTrack.pidbit() & (1 << 3)) && static_cast<bool>(negTrack.pidbit() & (1 << 3))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithdEdxTrQA(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSignalPi(posTrack.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton], maxPt4dwnsmplTsalisProtons)) {
            if (std::abs(negTrack.tofNSigmaPr()) <= nSigmaTOFdautrack) {
              fillSkimmedV0TableWithdEdxTrQA(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaPr(), negTrack.tofNSigmaPr(), negTrack.tpcExpSignalPr(negTrack.tpcSignal()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
            }
          }
        }
      }
    }
  } /// process with dEdx from TrackQA
  PROCESS_SWITCH(TreeWriterTpcV0, processWithdEdxTrQA, "Standard V0 Samples with dEdx from Track QA for PID", false);

  Preslice<TrksWithDEdxCorrection> perCollisionTracksWithNewDEdx = aod::track::collisionId;
  void processWithdEdxTrQAWithCorrecteddEdx(Colls const& collisions, TrksWithDEdxCorrection const& myTracks, V0sWithID const& myV0s, aod::BCsWithTimestamps const&, aod::TracksQAVersion const& tracksQA)
  {
    std::vector<int64_t> labelTrack2TrackQA;
    labelTrack2TrackQA.clear();
    labelTrack2TrackQA.resize(myTracks.size(), -1);
    for (const auto& trackQA : tracksQA) {
      int64_t trackId = trackQA.trackId();
      int64_t trackQAIndex = trackQA.globalIndex();
      labelTrack2TrackQA[trackId] = trackQAIndex;
    }
    for (const auto& collision : collisions) {
      auto tracks = myTracks.sliceBy(perCollisionTracksWithNewDEdx, collision.globalIndex());
      auto v0s = myV0s.sliceBy(perCollisionV0s, collision.globalIndex());
      /// Check event slection
      if (!isEventSelected(collision, tracks)) {
        continue;
      }
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      const int runnumber = bc.runNumber();
      float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * 1.e-3;
      rowTPCTreeWithTrkQA.reserve(tracks.size());
      /// Loop over v0 candidates
      for (const auto& v0 : v0s) {
        auto posTrack = v0.posTrack_as<TrksWithDEdxCorrection>();
        auto negTrack = v0.negTrack_as<TrksWithDEdxCorrection>();
        if (v0.v0addid() == -1) {
          continue;
        }
        aod::TracksQA posTrackQA;
        aod::TracksQA negTrackQA;
        bool existPosTrkQA;
        bool existNegTrkQA;
        if (labelTrack2TrackQA[posTrack.globalIndex()] != -1) {
          posTrackQA = tracksQA.iteratorAt(labelTrack2TrackQA[posTrack.globalIndex()]);
          existPosTrkQA = true;
        } else {
          posTrackQA = tracksQA.iteratorAt(0);
          existPosTrkQA = false;
        }
        if (labelTrack2TrackQA[negTrack.globalIndex()] != -1) {
          negTrackQA = tracksQA.iteratorAt(labelTrack2TrackQA[negTrack.globalIndex()]);
          existNegTrkQA = true;
        } else {
          negTrackQA = tracksQA.iteratorAt(0);
          existNegTrkQA = false;
        }

        // gamma
        if (static_cast<bool>(posTrack.pidbit() & (1 << 0)) && static_cast<bool>(negTrack.pidbit() & (1 << 0))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisElectrons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Electron], maxPt4dwnsmplTsalisElectrons)) {
            fillSkimmedV0TableWithdEdxTrQA<true>(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaEl(), posTrack.tofNSigmaEl(), posTrack.tpcExpSignalEl(posTrack.tpcSignalCorrected()), o2::track::PID::Electron, runnumber, dwnSmplFactor_El, hadronicRate);
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisElectrons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Electron], maxPt4dwnsmplTsalisElectrons)) {
            fillSkimmedV0TableWithdEdxTrQA<true>(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaEl(), negTrack.tofNSigmaEl(), negTrack.tpcExpSignalEl(negTrack.tpcSignalCorrected()), o2::track::PID::Electron, runnumber, dwnSmplFactor_El, hadronicRate);
          }
        }
        // Ks0
        if (static_cast<bool>(posTrack.pidbit() & (1 << 1)) && static_cast<bool>(negTrack.pidbit() & (1 << 1))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithdEdxTrQA<true>(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSignalPi(posTrack.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithdEdxTrQA<true>(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSignalPi(negTrack.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
          }
        }
        // Lambda
        if (static_cast<bool>(posTrack.pidbit() & (1 << 2)) && static_cast<bool>(negTrack.pidbit() & (1 << 2))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton], maxPt4dwnsmplTsalisProtons)) {
            if (std::abs(posTrack.tofNSigmaPr()) <= nSigmaTOFdautrack) {
              fillSkimmedV0TableWithdEdxTrQA<true>(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaPr(), posTrack.tofNSigmaPr(), posTrack.tpcExpSignalPr(posTrack.tpcSignal()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
            }
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithdEdxTrQA<true>(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSignalPi(negTrack.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
          }
        }
        // Antilambda
        if (static_cast<bool>(posTrack.pidbit() & (1 << 3)) && static_cast<bool>(negTrack.pidbit() & (1 << 3))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithdEdxTrQA<true>(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSignalPi(posTrack.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton], maxPt4dwnsmplTsalisProtons)) {
            if (std::abs(negTrack.tofNSigmaPr()) <= nSigmaTOFdautrack) {
              fillSkimmedV0TableWithdEdxTrQA<true>(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaPr(), negTrack.tofNSigmaPr(), negTrack.tpcExpSignalPr(negTrack.tpcSignalCorrected()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
            }
          }
        }
      }
    }
  } /// process with dEdx from TrackQA
  PROCESS_SWITCH(TreeWriterTpcV0, processWithdEdxTrQAWithCorrecteddEdx, "Standard V0 Samples with dEdx from Track QA for PID with corrected dEdx", false);

  void processWithTrQA(Colls const& collisions, Trks const& myTracks, V0sWithID const& myV0s, MyBCTable const&, aod::TracksQAVersion const& tracksQA)
  {
    std::vector<int64_t> labelTrack2TrackQA;
    labelTrack2TrackQA.clear();
    labelTrack2TrackQA.resize(myTracks.size(), -1);
    for (const auto& trackQA : tracksQA) {
      int64_t trackId = trackQA.trackId();
      int64_t trackQAIndex = trackQA.globalIndex();
      labelTrack2TrackQA[trackId] = trackQAIndex;
    }
    for (const auto& collision : collisions) {
      auto tracks = myTracks.sliceBy(perCollisionTracks, collision.globalIndex());
      auto v0s = myV0s.sliceBy(perCollisionV0s, collision.globalIndex());
      /// Check event slection
      if (!isEventSelected(collision, tracks)) {
        continue;
      }
      auto bc = collision.bc_as<MyBCTable>();
      const int runnumber = bc.runNumber();
      float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * 1.e-3;
      const int bcGlobalIndex = bc.globalIndex();
      const int bcTimeFrameId = bc.tfId();
      const int bcBcInTimeFrame = bc.bcInTF();
      rowTPCTreeWithTrkQA.reserve(tracks.size());
      /// Loop over v0 candidates
      for (const auto& v0 : v0s) {
        auto posTrack = v0.posTrack_as<Trks>();
        auto negTrack = v0.negTrack_as<Trks>();
        if (v0.v0addid() == -1) {
          continue;
        }
        aod::TracksQA posTrackQA;
        aod::TracksQA negTrackQA;
        bool existPosTrkQA;
        bool existNegTrkQA;
        if (labelTrack2TrackQA[posTrack.globalIndex()] != -1) {
          posTrackQA = tracksQA.iteratorAt(labelTrack2TrackQA[posTrack.globalIndex()]);
          existPosTrkQA = true;
        } else {
          posTrackQA = tracksQA.iteratorAt(0);
          existPosTrkQA = false;
        }
        if (labelTrack2TrackQA[negTrack.globalIndex()] != -1) {
          negTrackQA = tracksQA.iteratorAt(labelTrack2TrackQA[negTrack.globalIndex()]);
          existNegTrkQA = true;
        } else {
          negTrackQA = tracksQA.iteratorAt(0);
          existNegTrkQA = false;
        }

        // gamma
        if (static_cast<bool>(posTrack.pidbit() & (1 << 0)) && static_cast<bool>(negTrack.pidbit() & (1 << 0))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisElectrons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Electron], maxPt4dwnsmplTsalisElectrons)) {
            fillSkimmedV0TableWithTrQA(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaEl(), posTrack.tofNSigmaEl(), posTrack.tpcExpSignalEl(posTrack.tpcSignal()), o2::track::PID::Electron, runnumber, dwnSmplFactor_El, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisElectrons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Electron], maxPt4dwnsmplTsalisElectrons)) {
            fillSkimmedV0TableWithTrQA(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaEl(), negTrack.tofNSigmaEl(), negTrack.tpcExpSignalEl(negTrack.tpcSignal()), o2::track::PID::Electron, runnumber, dwnSmplFactor_El, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
          }
        }
        // Ks0
        if (static_cast<bool>(posTrack.pidbit() & (1 << 1)) && static_cast<bool>(negTrack.pidbit() & (1 << 1))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithTrQA(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSignalPi(posTrack.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithTrQA(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSignalPi(negTrack.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
          }
        }
        // Lambda
        if (static_cast<bool>(posTrack.pidbit() & (1 << 2)) && static_cast<bool>(negTrack.pidbit() & (1 << 2))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton], maxPt4dwnsmplTsalisProtons)) {
            if (std::abs(posTrack.tofNSigmaPr()) <= nSigmaTOFdautrack) {
              fillSkimmedV0TableWithTrQA(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaPr(), posTrack.tofNSigmaPr(), posTrack.tpcExpSignalPr(posTrack.tpcSignal()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
            }
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithTrQA(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSignalPi(negTrack.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
          }
        }
        // Antilambda
        if (static_cast<bool>(posTrack.pidbit() & (1 << 3)) && static_cast<bool>(negTrack.pidbit() & (1 << 3))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithTrQA(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSignalPi(posTrack.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton], maxPt4dwnsmplTsalisProtons)) {
            if (std::abs(negTrack.tofNSigmaPr()) <= nSigmaTOFdautrack) {
              fillSkimmedV0TableWithTrQA(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaPr(), negTrack.tofNSigmaPr(), negTrack.tpcExpSignalPr(negTrack.tpcSignal()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
            }
          }
        }
      }
    }
  } /// process with TrackQA
  PROCESS_SWITCH(TreeWriterTpcV0, processWithTrQA, "Standard V0 Samples with Track QA for PID", false);

  void processWithTrQAWithCorrecteddEdx(Colls const& collisions, TrksWithDEdxCorrection const& myTracks, V0sWithID const& myV0s, MyBCTable const&, aod::TracksQAVersion const& tracksQA)
  {
    std::vector<int64_t> labelTrack2TrackQA;
    labelTrack2TrackQA.clear();
    labelTrack2TrackQA.resize(myTracks.size(), -1);
    for (const auto& trackQA : tracksQA) {
      int64_t trackId = trackQA.trackId();
      int64_t trackQAIndex = trackQA.globalIndex();
      labelTrack2TrackQA[trackId] = trackQAIndex;
    }
    for (const auto& collision : collisions) {
      auto tracks = myTracks.sliceBy(perCollisionTracksWithNewDEdx, collision.globalIndex());
      auto v0s = myV0s.sliceBy(perCollisionV0s, collision.globalIndex());
      /// Check event slection
      if (!isEventSelected(collision, tracks)) {
        continue;
      }
      auto bc = collision.bc_as<MyBCTable>();
      const int runnumber = bc.runNumber();
      float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * 1.e-3;
      const int bcGlobalIndex = bc.globalIndex();
      const int bcTimeFrameId = bc.tfId();
      const int bcBcInTimeFrame = bc.bcInTF();
      rowTPCTreeWithTrkQA.reserve(tracks.size());
      /// Loop over v0 candidates
      for (const auto& v0 : v0s) {
        auto posTrack = v0.posTrack_as<TrksWithDEdxCorrection>();
        auto negTrack = v0.negTrack_as<TrksWithDEdxCorrection>();
        if (v0.v0addid() == -1) {
          continue;
        }
        aod::TracksQA posTrackQA;
        aod::TracksQA negTrackQA;
        bool existPosTrkQA;
        bool existNegTrkQA;
        if (labelTrack2TrackQA[posTrack.globalIndex()] != -1) {
          posTrackQA = tracksQA.iteratorAt(labelTrack2TrackQA[posTrack.globalIndex()]);
          existPosTrkQA = true;
        } else {
          posTrackQA = tracksQA.iteratorAt(0);
          existPosTrkQA = false;
        }
        if (labelTrack2TrackQA[negTrack.globalIndex()] != -1) {
          negTrackQA = tracksQA.iteratorAt(labelTrack2TrackQA[negTrack.globalIndex()]);
          existNegTrkQA = true;
        } else {
          negTrackQA = tracksQA.iteratorAt(0);
          existNegTrkQA = false;
        }

        // gamma
        if (static_cast<bool>(posTrack.pidbit() & (1 << 0)) && static_cast<bool>(negTrack.pidbit() & (1 << 0))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisElectrons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Electron], maxPt4dwnsmplTsalisElectrons)) {
            fillSkimmedV0TableWithTrQA<true>(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaEl(), posTrack.tofNSigmaEl(), posTrack.tpcExpSignalEl(posTrack.tpcSignalCorrected()), o2::track::PID::Electron, runnumber, dwnSmplFactor_El, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisElectrons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Electron], maxPt4dwnsmplTsalisElectrons)) {
            fillSkimmedV0TableWithTrQA<true>(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaEl(), negTrack.tofNSigmaEl(), negTrack.tpcExpSignalEl(negTrack.tpcSignalCorrected()), o2::track::PID::Electron, runnumber, dwnSmplFactor_El, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
          }
        }
        // Ks0
        if (static_cast<bool>(posTrack.pidbit() & (1 << 1)) && static_cast<bool>(negTrack.pidbit() & (1 << 1))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithTrQA<true>(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSignalPi(posTrack.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithTrQA<true>(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSignalPi(negTrack.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
          }
        }
        // Lambda
        if (static_cast<bool>(posTrack.pidbit() & (1 << 2)) && static_cast<bool>(negTrack.pidbit() & (1 << 2))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton], maxPt4dwnsmplTsalisProtons)) {
            if (std::abs(posTrack.tofNSigmaPr()) <= nSigmaTOFdautrack) {
              fillSkimmedV0TableWithTrQA<true>(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaPr(), posTrack.tofNSigmaPr(), posTrack.tpcExpSignalPr(posTrack.tpcSignalCorrected()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
            }
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithTrQA<true>(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaPi(), negTrack.tofNSigmaPi(), negTrack.tpcExpSignalPi(negTrack.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
          }
        }
        // Antilambda
        if (static_cast<bool>(posTrack.pidbit() & (1 << 3)) && static_cast<bool>(negTrack.pidbit() & (1 << 3))) {
          if (downsampleTsalisCharged(posTrack.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion], maxPt4dwnsmplTsalisPions)) {
            fillSkimmedV0TableWithTrQA<true>(v0, posTrack, posTrackQA, existPosTrkQA, collision, posTrack.tpcNSigmaPi(), posTrack.tofNSigmaPi(), posTrack.tpcExpSignalPi(posTrack.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
          }
          if (downsampleTsalisCharged(negTrack.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton], maxPt4dwnsmplTsalisProtons)) {
            if (std::abs(negTrack.tofNSigmaPr()) <= nSigmaTOFdautrack) {
              fillSkimmedV0TableWithTrQA<true>(v0, negTrack, negTrackQA, existNegTrkQA, collision, negTrack.tpcNSigmaPr(), negTrack.tofNSigmaPr(), negTrack.tpcExpSignalPr(negTrack.tpcSignalCorrected()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
            }
          }
        }
      }
    }
  } /// process with TrackQA
  PROCESS_SWITCH(TreeWriterTpcV0, processWithTrQAWithCorrecteddEdx, "Standard V0 Samples with Track QA for PID with corrected dEdx", false);

  void processDummy(Colls const&) {}
  PROCESS_SWITCH(TreeWriterTpcV0, processDummy, "Dummy function", false);

}; /// struct TreeWriterTpcV0

struct TreeWriterTPCTOF {

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

  /// Configurables
  Configurable<float> nClNorm{"nClNorm", 152., "Number of cluster normalization. Run 2: 159, Run 3 152"};
  Configurable<int> applyEvSel{"applyEvSel", 2, "Flag to apply rapidity cut: 0 -> no event selection, 1 -> Run 2 event selection, 2 -> Run 3 event selection"};
  Configurable<int> applyTrkSel{"applyTrkSel", 1, "Flag to apply track selection: 0 -> no track selection, 1 -> track selection"};
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<std::string> irSource{"irSource", "T0VTX", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};
  /// Triton
  Configurable<float> maxMomTPCOnlyTr{"maxMomTPCOnlyTr", 1.5, "Maximum momentum for TPC only cut triton"};
  Configurable<float> maxMomHardCutOnlyTr{"maxMomHardCutOnlyTr", 50, "Maximum TPC inner momentum for triton"};
  Configurable<float> nSigmaTPCOnlyTr{"nSigmaTPCOnlyTr", 4., "number of sigma for TPC only cut triton"};
  Configurable<float> nSigmaTPC_TPCTOF_Tr{"nSigmaTPC_TPCTOF_Tr", 4., "number of sigma for TPC cut for TPC and TOF combined triton"};
  Configurable<float> nSigmaTOF_TPCTOF_Tr{"nSigmaTOF_TPCTOF_Tr", 3., "number of sigma for TOF cut for TPC and TOF combined triton"};
  Configurable<double> dwnSmplFactor_Tr{"dwnSmplFactor_Tr", 1., "downsampling factor for triton, default fraction to keep is 1."};
  /// Deuteron
  Configurable<float> maxMomTPCOnlyDe{"maxMomTPCOnlyDe", 1.0, "Maximum momentum for TPC only cut deuteron"};
  Configurable<float> maxMomHardCutOnlyDe{"maxMomHardCutOnlyDe", 50, "Maximum TPC inner momentum for deuteron"};
  Configurable<float> nSigmaTPCOnlyDe{"nSigmaTPCOnlyDe", 4., "number of sigma for TPC only cut deuteron"};
  Configurable<float> nSigmaTPC_TPCTOF_De{"nSigmaTPC_TPCTOF_De", 4., "number of sigma for TPC cut for TPC and TOF combined deuteron"};
  Configurable<float> nSigmaTOF_TPCTOF_De{"nSigmaTOF_TPCTOF_De", 3., "number of sigma for TOF cut for TPC and TOF combined deuteron"};
  Configurable<double> dwnSmplFactor_De{"dwnSmplFactor_De", 1., "downsampling factor for deuteron, default fraction to keep is 1."};
  /// Proton
  Configurable<float> maxMomTPCOnlyPr{"maxMomTPCOnlyPr", 0.6, "Maximum momentum for TPC only cut proton"};
  Configurable<float> nSigmaTPCOnlyPr{"nSigmaTPCOnlyPr", 4., "number of sigma for TPC only cut proton"};
  Configurable<float> nSigmaTPC_TPCTOF_Pr{"nSigmaTPC_TPCTOF_Pr", 4., "number of sigma for TPC cut for TPC and TOF combined proton"};
  Configurable<float> nSigmaTOF_TPCTOF_Pr{"nSigmaTOF_TPCTOF_Pr", 3., "number of sigma for TOF cut for TPC and TOF combined proton"};
  Configurable<double> dwnSmplFactor_Pr{"dwnSmplFactor_Pr", 1., "downsampling factor for protons, default fraction to keep is 1."};
  /// Kaon
  Configurable<float> maxMomTPCOnlyKa{"maxMomTPCOnlyKa", 0.3, "Maximum momentum for TPC only cut kaon"};
  Configurable<float> maxMomHardCutOnlyKa{"maxMomHardCutOnlyKa", 50, "Maximum TPC inner momentum for kaons"};
  Configurable<float> nSigmaTPCOnlyKa{"nSigmaTPCOnlyKa", 4., "number of sigma for TPC only cut kaon"};
  Configurable<float> nSigmaTPC_TPCTOF_Ka{"nSigmaTPC_TPCTOF_Ka", 4., "number of sigma for TPC cut for TPC and TOF combined kaon"};
  Configurable<float> nSigmaTOF_TPCTOF_Ka{"nSigmaTOF_TPCTOF_Ka", 3., "number of sigma for TOF cut for TPC and TOF combined kaon"};
  Configurable<double> dwnSmplFactor_Ka{"dwnSmplFactor_Ka", 1., "downsampling factor for kaons, default fraction to keep is 1."};
  /// Pion
  Configurable<float> maxMomTPCOnlyPi{"maxMomTPCOnlyPi", 0.5, "Maximum momentum for TPC only cut pion"};
  Configurable<float> nSigmaTPCOnlyPi{"nSigmaTPCOnlyPi", 4., "number of sigma for TPC only cut pion"};
  Configurable<float> nSigmaTPC_TPCTOF_Pi{"nSigmaTPC_TPCTOF_Pi", 4., "number of sigma for TPC cut for TPC and TOF combined pion"};
  Configurable<float> nSigmaTOF_TPCTOF_Pi{"nSigmaTOF_TPCTOF_Pi", 4., "number of sigma for TOF cut for TPC and TOF combined pion"};
  Configurable<double> dwnSmplFactor_Pi{"dwnSmplFactor_Pi", 1., "downsampling factor for pions, default fraction to keep is 1."};
  /// pT dependent downsampling
  Configurable<float> sqrtSNN{"sqrt_s_NN", 0., "sqrt(s_NN), used for downsampling with the Tsallis distribution"};
  Configurable<float> downsamplingTsalisTritons{"downsamplingTsalisTritons", -1., "Downsampling factor to reduce the number of tritons"};
  Configurable<float> downsamplingTsalisDeuterons{"downsamplingTsalisDeuterons", -1., "Downsampling factor to reduce the number of deuterons"};
  Configurable<float> downsamplingTsalisProtons{"downsamplingTsalisProtons", -1., "Downsampling factor to reduce the number of protons"};
  Configurable<float> downsamplingTsalisKaons{"downsamplingTsalisKaons", -1., "Downsampling factor to reduce the number of kaons"};
  Configurable<float> downsamplingTsalisPions{"downsamplingTsalisPions", -1., "Downsampling factor to reduce the number of pions"};

  Filter trackFilter = (trackSelection.node() == 0) ||
                       ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                       ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                       ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                       ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));

  ctpRateFetcher mRateFetcher;

  double tsalisCharged(double pt, double mass, double sqrts)
  {
    const double a = 6.81, b = 59.24;
    const double c = 0.082, d = 0.151;
    double mt = std::sqrt(mass * mass + pt * pt);
    double n = a + b / sqrts;
    double T = c + d / sqrts;
    double p0 = n * T;
    double result = std::pow((1. + mt / p0), -n);
    return result;
  };

  /// Random downsampling trigger function using Tsalis/Hagedorn spectra fit (sqrt(s) = 62.4 GeV to 13 TeV)
  /// as in https://iopscience.iop.org/article/10.1088/2399-6528/aab00f/pdf
  TRandom3* fRndm = new TRandom3(0);
  bool downsampleTsalisCharged(double pt, float factor1Pt, double sqrts, double mass)
  {
    if (factor1Pt < 0.) {
      return true;
    }
    const double prob = tsalisCharged(pt, mass, sqrts) * pt;
    const double probNorm = tsalisCharged(1., mass, sqrts);
    if ((fRndm->Rndm() * ((prob / probNorm) * pt * pt)) > factor1Pt) {
      return false;
    } else {
      return true;
    }
  };

  /// Function to fill trees
  template <bool doCorrectdEdx = false, typename T, typename C>
  void fillSkimmedTPCTOFTable(T const& track, C const& collision, const float nSigmaTPC, const float nSigmaTOF, const float dEdxExp, const o2::track::PID::ID id, int runnumber, double dwnSmplFactor, double hadronicRate)
  {

    const double ncl = track.tpcNClsFound();
    const double nclPID = track.tpcNClsFindableMinusPID();
    const double p = track.tpcInnerParam();
    const double mass = o2::track::pid_constants::sMasses[id];
    const double bg = p / mass;
    const int multTPC = collision.multTPC();
    auto trackocc = collision.trackOccupancyInTimeRange();
    auto ft0occ = collision.ft0cOccupancyInTimeRange();

    const double pseudoRndm = track.pt() * 1000. - static_cast<int64_t>(track.pt() * 1000);
    if (pseudoRndm < dwnSmplFactor) {
      float usedEdx;
      if constexpr (doCorrectdEdx) {
        usedEdx = track.tpcSignalCorrected();
      } else {
        usedEdx = track.tpcSignal();
      }
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
                    multTPC / 11000.,
                    std::sqrt(nClNorm / ncl),
                    nclPID,
                    id,
                    nSigmaTPC,
                    nSigmaTOF,
                    runnumber,
                    trackocc,
                    ft0occ,
                    hadronicRate);
    }
  };
  template <bool doCorrectdEdx = false, typename T, typename TQA, typename C>
  void fillSkimmedTPCTOFTableWithdEdxTrkQA(T const& track, TQA const& trackQA, bool existTrkQA, C const& collision, const float nSigmaTPC, const float nSigmaTOF, const float nSigmaITS, const float dEdxExp, const o2::track::PID::ID id, int runnumber, double dwnSmplFactor, double hadronicRate)
  {

    const double ncl = track.tpcNClsFound();
    const double nclPID = track.tpcNClsFindableMinusPID();
    const double p = track.tpcInnerParam();
    const double mass = o2::track::pid_constants::sMasses[id];
    const double bg = p / mass;
    const int multTPC = collision.multTPC();
    auto trackocc = collision.trackOccupancyInTimeRange();
    auto ft0occ = collision.ft0cOccupancyInTimeRange();

    const double pseudoRndm = track.pt() * 1000. - static_cast<int64_t>(track.pt() * 1000);
    if (pseudoRndm < dwnSmplFactor) {
      float usedEdx;
      if constexpr (doCorrectdEdx) {
        usedEdx = track.tpcSignalCorrected();
      } else {
        usedEdx = track.tpcSignal();
      }
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
                                 multTPC / 11000.,
                                 std::sqrt(nClNorm / ncl),
                                 nclPID,
                                 id,
                                 nSigmaTPC,
                                 nSigmaTOF,
                                 nSigmaITS,
                                 runnumber,
                                 trackocc,
                                 ft0occ,
                                 hadronicRate,
                                 existTrkQA ? trackQA.tpcdEdxNorm() : -999);
    }
  };
  /// Function to fill trees
  template <bool doCorrectdEdx = false, typename T, typename TQA, typename C>
  void fillSkimmedTPCTOFTableWithTrkQA(T const& track, TQA const& trackQA, bool existTrkQA, C const& collision, const float nSigmaTPC, const float nSigmaTOF, const float nSigmaITS, const float dEdxExp, const o2::track::PID::ID id, int runnumber, double dwnSmplFactor, double hadronicRate, int bcGlobalIndex, int bcTimeFrameId, int bcBcInTimeFrame)
  {

    const double ncl = track.tpcNClsFound();
    const double nclPID = track.tpcNClsFindableMinusPID();
    const double p = track.tpcInnerParam();
    const double mass = o2::track::pid_constants::sMasses[id];
    const double bg = p / mass;
    const int multTPC = collision.multTPC();
    auto trackocc = collision.trackOccupancyInTimeRange();
    auto ft0occ = collision.ft0cOccupancyInTimeRange();

    const double pseudoRndm = track.pt() * 1000. - static_cast<int64_t>(track.pt() * 1000);
    if (pseudoRndm < dwnSmplFactor) {
      float usedEdx;
      if constexpr (doCorrectdEdx) {
        usedEdx = track.tpcSignalCorrected();
      } else {
        usedEdx = track.tpcSignal();
      }
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
                             multTPC / 11000.,
                             std::sqrt(nClNorm / ncl),
                             nclPID,
                             id,
                             nSigmaTPC,
                             nSigmaTOF,
                             nSigmaITS,
                             runnumber,
                             trackocc,
                             ft0occ,
                             hadronicRate,
                             bcGlobalIndex,
                             bcTimeFrameId,
                             bcBcInTimeFrame,
                             existTrkQA ? trackQA.tpcClusterByteMask() : -999,
                             existTrkQA ? trackQA.tpcdEdxMax0R() : -999,
                             existTrkQA ? trackQA.tpcdEdxMax1R() : -999,
                             existTrkQA ? trackQA.tpcdEdxMax2R() : -999,
                             existTrkQA ? trackQA.tpcdEdxMax3R() : -999,
                             existTrkQA ? trackQA.tpcdEdxTot0R() : -999,
                             existTrkQA ? trackQA.tpcdEdxTot1R() : -999,
                             existTrkQA ? trackQA.tpcdEdxTot2R() : -999,
                             existTrkQA ? trackQA.tpcdEdxTot3R() : -999,
                             existTrkQA ? trackQA.tpcdEdxNorm() : -999);
    }
  };

  /// Event selection
  template <typename CollisionType, typename TrackType>
  bool isEventSelected(const CollisionType& collision, const TrackType& /*tracks*/)
  {
    if (applyEvSel == 1) {
      if (!collision.sel7()) {
        return false;
      }
    } else if (applyEvSel == 2) {
      if (!collision.sel8()) {
        return false;
      }
    }
    return true;
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);
  }

  void processStandard(Colls::iterator const& collision, soa::Filtered<Trks> const& tracks, aod::BCsWithTimestamps const&)
  {
    /// Check event selection
    if (!isEventSelected(collision, tracks)) {
      return;
    }
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    const int runnumber = bc.runNumber();
    float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * 1.e-3;

    rowTPCTOFTree.reserve(tracks.size());
    for (auto const& trk : tracks) {
      /// Fill tree for tritons
      if (trk.tpcInnerParam() < maxMomHardCutOnlyTr && trk.tpcInnerParam() <= maxMomTPCOnlyTr && std::abs(trk.tpcNSigmaTr()) < nSigmaTPCOnlyTr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Triton])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaTr(), trk.tofNSigmaTr(), trk.tpcExpSignalTr(trk.tpcSignal()), o2::track::PID::Triton, runnumber, dwnSmplFactor_Tr, hadronicRate);
      } else if (trk.tpcInnerParam() < maxMomHardCutOnlyTr && trk.tpcInnerParam() > maxMomTPCOnlyTr && std::abs(trk.tofNSigmaTr()) < nSigmaTOF_TPCTOF_Tr && std::abs(trk.tpcNSigmaTr()) < nSigmaTPC_TPCTOF_Tr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Triton])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaTr(), trk.tofNSigmaTr(), trk.tpcExpSignalTr(trk.tpcSignal()), o2::track::PID::Triton, runnumber, dwnSmplFactor_Tr, hadronicRate);
      }
      /// Fill tree for deuterons
      if (trk.tpcInnerParam() < maxMomHardCutOnlyDe && trk.tpcInnerParam() <= maxMomTPCOnlyDe && std::abs(trk.tpcNSigmaDe()) < nSigmaTPCOnlyDe && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Deuteron])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaDe(), trk.tofNSigmaDe(), trk.tpcExpSignalDe(trk.tpcSignal()), o2::track::PID::Deuteron, runnumber, dwnSmplFactor_De, hadronicRate);
      } else if (trk.tpcInnerParam() < maxMomHardCutOnlyDe && trk.tpcInnerParam() > maxMomTPCOnlyDe && std::abs(trk.tofNSigmaDe()) < nSigmaTOF_TPCTOF_De && std::abs(trk.tpcNSigmaDe()) < nSigmaTPC_TPCTOF_De && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Deuteron])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaDe(), trk.tofNSigmaDe(), trk.tpcExpSignalDe(trk.tpcSignal()), o2::track::PID::Deuteron, runnumber, dwnSmplFactor_De, hadronicRate);
      }
      /// Fill tree for protons
      if (trk.tpcInnerParam() <= maxMomTPCOnlyPr && std::abs(trk.tpcNSigmaPr()) < nSigmaTPCOnlyPr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaPr(), trk.tofNSigmaPr(), trk.tpcExpSignalPr(trk.tpcSignal()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
      } else if (trk.tpcInnerParam() > maxMomTPCOnlyPr && std::abs(trk.tofNSigmaPr()) < nSigmaTOF_TPCTOF_Pr && std::abs(trk.tpcNSigmaPr()) < nSigmaTPC_TPCTOF_Pr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaPr(), trk.tofNSigmaPr(), trk.tpcExpSignalPr(trk.tpcSignal()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
      }
      /// Fill tree for kaons
      if (trk.tpcInnerParam() < maxMomHardCutOnlyKa && trk.tpcInnerParam() <= maxMomTPCOnlyKa && std::abs(trk.tpcNSigmaKa()) < nSigmaTPCOnlyKa && downsampleTsalisCharged(trk.pt(), downsamplingTsalisKaons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Kaon])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaKa(), trk.tofNSigmaKa(), trk.tpcExpSignalKa(trk.tpcSignal()), o2::track::PID::Kaon, runnumber, dwnSmplFactor_Ka, hadronicRate);
      } else if (trk.tpcInnerParam() < maxMomHardCutOnlyKa && trk.tpcInnerParam() > maxMomTPCOnlyKa && std::abs(trk.tofNSigmaKa()) < nSigmaTOF_TPCTOF_Ka && std::abs(trk.tpcNSigmaKa()) < nSigmaTPC_TPCTOF_Ka && downsampleTsalisCharged(trk.pt(), downsamplingTsalisKaons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Kaon])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaKa(), trk.tofNSigmaKa(), trk.tpcExpSignalKa(trk.tpcSignal()), o2::track::PID::Kaon, runnumber, dwnSmplFactor_Ka, hadronicRate);
      }
      /// Fill tree pions
      if (trk.tpcInnerParam() <= maxMomTPCOnlyPi && std::abs(trk.tpcNSigmaPi()) < nSigmaTPCOnlyPi && downsampleTsalisCharged(trk.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaPi(), trk.tofNSigmaPi(), trk.tpcExpSignalPi(trk.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
      } else if (trk.tpcInnerParam() > maxMomTPCOnlyPi && std::abs(trk.tofNSigmaPi()) < nSigmaTOF_TPCTOF_Pi && std::abs(trk.tpcNSigmaPi()) < nSigmaTPC_TPCTOF_Pi && downsampleTsalisCharged(trk.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion])) {
        fillSkimmedTPCTOFTable(trk, collision, trk.tpcNSigmaPi(), trk.tofNSigmaPi(), trk.tpcExpSignalPi(trk.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
      }
    } /// Loop tracks
  } /// process
  PROCESS_SWITCH(TreeWriterTPCTOF, processStandard, "Standard Samples for PID", true);

  void processStandard2(Colls::iterator const& collision, soa::Filtered<TrksWithDEdxCorrection> const& tracks, aod::BCsWithTimestamps const&)
  {
    /// Check event selection
    if (!isEventSelected(collision, tracks)) {
      return;
    }
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    const int runnumber = bc.runNumber();
    float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * 1.e-3;

    rowTPCTOFTree.reserve(tracks.size());
    for (auto const& trk : tracks) {
      /// Fill tree for tritons/*  */
      if (trk.tpcInnerParam() < maxMomHardCutOnlyTr && trk.tpcInnerParam() <= maxMomTPCOnlyTr && std::abs(trk.tpcNSigmaTr()) < nSigmaTPCOnlyTr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Triton])) {
        fillSkimmedTPCTOFTable<true>(trk, collision, trk.tpcNSigmaTr(), trk.tofNSigmaTr(), trk.tpcExpSignalTr(trk.tpcSignalCorrected()), o2::track::PID::Triton, runnumber, dwnSmplFactor_Tr, hadronicRate);
      } else if (trk.tpcInnerParam() < maxMomHardCutOnlyTr && trk.tpcInnerParam() > maxMomTPCOnlyTr && std::abs(trk.tofNSigmaTr()) < nSigmaTOF_TPCTOF_Tr && std::abs(trk.tpcNSigmaTr()) < nSigmaTPC_TPCTOF_Tr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Triton])) {
        fillSkimmedTPCTOFTable<true>(trk, collision, trk.tpcNSigmaTr(), trk.tofNSigmaTr(), trk.tpcExpSignalTr(trk.tpcSignalCorrected()), o2::track::PID::Triton, runnumber, dwnSmplFactor_Tr, hadronicRate);
      }
      /// Fill tree for deuterons
      if (trk.tpcInnerParam() < maxMomHardCutOnlyDe && trk.tpcInnerParam() <= maxMomTPCOnlyDe && std::abs(trk.tpcNSigmaDe()) < nSigmaTPCOnlyDe && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Deuteron])) {
        fillSkimmedTPCTOFTable<true>(trk, collision, trk.tpcNSigmaDe(), trk.tofNSigmaDe(), trk.tpcExpSignalDe(trk.tpcSignalCorrected()), o2::track::PID::Deuteron, runnumber, dwnSmplFactor_De, hadronicRate);
      } else if (trk.tpcInnerParam() < maxMomHardCutOnlyDe && trk.tpcInnerParam() > maxMomTPCOnlyDe && std::abs(trk.tofNSigmaDe()) < nSigmaTOF_TPCTOF_De && std::abs(trk.tpcNSigmaDe()) < nSigmaTPC_TPCTOF_De && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Deuteron])) {
        fillSkimmedTPCTOFTable<true>(trk, collision, trk.tpcNSigmaDe(), trk.tofNSigmaDe(), trk.tpcExpSignalDe(trk.tpcSignalCorrected()), o2::track::PID::Deuteron, runnumber, dwnSmplFactor_De, hadronicRate);
      }
      /// Fill tree for protons
      if (trk.tpcInnerParam() <= maxMomTPCOnlyPr && std::abs(trk.tpcNSigmaPr()) < nSigmaTPCOnlyPr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton])) {
        fillSkimmedTPCTOFTable<true>(trk, collision, trk.tpcNSigmaPr(), trk.tofNSigmaPr(), trk.tpcExpSignalPr(trk.tpcSignalCorrected()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
      } else if (trk.tpcInnerParam() > maxMomTPCOnlyPr && std::abs(trk.tofNSigmaPr()) < nSigmaTOF_TPCTOF_Pr && std::abs(trk.tpcNSigmaPr()) < nSigmaTPC_TPCTOF_Pr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton])) {
        fillSkimmedTPCTOFTable<true>(trk, collision, trk.tpcNSigmaPr(), trk.tofNSigmaPr(), trk.tpcExpSignalPr(trk.tpcSignalCorrected()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
      }
      /// Fill tree for kaons
      if (trk.tpcInnerParam() < maxMomHardCutOnlyKa && trk.tpcInnerParam() <= maxMomTPCOnlyKa && std::abs(trk.tpcNSigmaKa()) < nSigmaTPCOnlyKa && downsampleTsalisCharged(trk.pt(), downsamplingTsalisKaons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Kaon])) {
        fillSkimmedTPCTOFTable<true>(trk, collision, trk.tpcNSigmaKa(), trk.tofNSigmaKa(), trk.tpcExpSignalKa(trk.tpcSignalCorrected()), o2::track::PID::Kaon, runnumber, dwnSmplFactor_Ka, hadronicRate);
      } else if (trk.tpcInnerParam() < maxMomHardCutOnlyKa && trk.tpcInnerParam() > maxMomTPCOnlyKa && std::abs(trk.tofNSigmaKa()) < nSigmaTOF_TPCTOF_Ka && std::abs(trk.tpcNSigmaKa()) < nSigmaTPC_TPCTOF_Ka && downsampleTsalisCharged(trk.pt(), downsamplingTsalisKaons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Kaon])) {
        fillSkimmedTPCTOFTable<true>(trk, collision, trk.tpcNSigmaKa(), trk.tofNSigmaKa(), trk.tpcExpSignalKa(trk.tpcSignalCorrected()), o2::track::PID::Kaon, runnumber, dwnSmplFactor_Ka, hadronicRate);
      }
      /// Fill tree pions
      if (trk.tpcInnerParam() <= maxMomTPCOnlyPi && std::abs(trk.tpcNSigmaPi()) < nSigmaTPCOnlyPi && downsampleTsalisCharged(trk.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion])) {
        fillSkimmedTPCTOFTable<true>(trk, collision, trk.tpcNSigmaPi(), trk.tofNSigmaPi(), trk.tpcExpSignalPi(trk.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
      } else if (trk.tpcInnerParam() > maxMomTPCOnlyPi && std::abs(trk.tofNSigmaPi()) < nSigmaTOF_TPCTOF_Pi && std::abs(trk.tpcNSigmaPi()) < nSigmaTPC_TPCTOF_Pi && downsampleTsalisCharged(trk.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion])) {
        fillSkimmedTPCTOFTable<true>(trk, collision, trk.tpcNSigmaPi(), trk.tofNSigmaPi(), trk.tpcExpSignalPi(trk.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
      }
    } /// Loop tracks
  } /// process
  PROCESS_SWITCH(TreeWriterTPCTOF, processStandard2, "Standard Samples for PID with corrected dEdx", false);

  Preslice<Trks> perCollisionTracks = aod::track::collisionId;
  void processWithdEdxTrQA(Colls const& collisions, Trks const& myTracks, aod::BCsWithTimestamps const&, aod::TracksQAVersion const& tracksQA)
  {
    std::vector<int64_t> labelTrack2TrackQA;
    labelTrack2TrackQA.clear();
    labelTrack2TrackQA.resize(myTracks.size(), -1);
    for (const auto& trackQA : tracksQA) {
      int64_t trackId = trackQA.trackId();
      int64_t trackQAIndex = trackQA.globalIndex();
      labelTrack2TrackQA[trackId] = trackQAIndex;
    }
    for (const auto& collision : collisions) {
      auto tracks = myTracks.sliceBy(perCollisionTracks, collision.globalIndex());
      auto tracksWithITSPid = soa::Attach<Trks,
                                          aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi,
                                          aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe,
                                          aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe, aod::pidits::ITSNSigmaAl>(tracks);
      /// Check event selection
      if (!isEventSelected(collision, tracks)) {
        continue;
      }
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      const int runnumber = bc.runNumber();
      float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * 1.e-3;
      rowTPCTOFTreeWithTrkQA.reserve(tracks.size());
      for (auto const& trk : tracksWithITSPid) {
        if (!((trackSelection == 0) ||
              ((trackSelection == 1) && trk.isGlobalTrack()) ||
              ((trackSelection == 2) && trk.isGlobalTrackWoPtEta()) ||
              ((trackSelection == 3) && trk.isGlobalTrackWoDCA()) ||
              ((trackSelection == 4) && trk.isQualityTrack()) ||
              ((trackSelection == 5) && trk.isInAcceptanceTrack()))) {
          continue;
        }
        // get the corresponding trackQA using labelTracks2TracKQA and get variables of interest
        aod::TracksQA trackQA;
        bool existTrkQA;
        if (labelTrack2TrackQA[trk.globalIndex()] != -1) {
          trackQA = tracksQA.iteratorAt(labelTrack2TrackQA[trk.globalIndex()]);
          existTrkQA = true;
        } else {
          trackQA = tracksQA.iteratorAt(0);
          existTrkQA = false;
        }
        /// Fill tree for tritons
        if (trk.tpcInnerParam() < maxMomHardCutOnlyTr && trk.tpcInnerParam() <= maxMomTPCOnlyTr && std::abs(trk.tpcNSigmaTr()) < nSigmaTPCOnlyTr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Triton])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaTr(), trk.tofNSigmaTr(), trk.itsNSigmaTr(), trk.tpcExpSignalTr(trk.tpcSignal()), o2::track::PID::Triton, runnumber, dwnSmplFactor_Tr, hadronicRate);
        } else if (trk.tpcInnerParam() < maxMomHardCutOnlyTr && trk.tpcInnerParam() > maxMomTPCOnlyTr && std::abs(trk.tofNSigmaTr()) < nSigmaTOF_TPCTOF_Tr && std::abs(trk.tpcNSigmaTr()) < nSigmaTPC_TPCTOF_Tr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Triton])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaTr(), trk.tofNSigmaTr(), trk.itsNSigmaTr(), trk.tpcExpSignalTr(trk.tpcSignal()), o2::track::PID::Triton, runnumber, dwnSmplFactor_Tr, hadronicRate);
        }
        /// Fill tree for deuterons
        if (trk.tpcInnerParam() < maxMomHardCutOnlyDe && trk.tpcInnerParam() <= maxMomTPCOnlyDe && std::abs(trk.tpcNSigmaDe()) < nSigmaTPCOnlyDe && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Deuteron])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaDe(), trk.tofNSigmaDe(), trk.itsNSigmaDe(), trk.tpcExpSignalDe(trk.tpcSignal()), o2::track::PID::Deuteron, runnumber, dwnSmplFactor_De, hadronicRate);
        } else if (trk.tpcInnerParam() < maxMomHardCutOnlyDe && trk.tpcInnerParam() > maxMomTPCOnlyDe && std::abs(trk.tofNSigmaDe()) < nSigmaTOF_TPCTOF_De && std::abs(trk.tpcNSigmaDe()) < nSigmaTPC_TPCTOF_De && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Deuteron])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaDe(), trk.tofNSigmaDe(), trk.itsNSigmaDe(), trk.tpcExpSignalDe(trk.tpcSignal()), o2::track::PID::Deuteron, runnumber, dwnSmplFactor_De, hadronicRate);
        }
        /// Fill tree for protons
        if (trk.tpcInnerParam() <= maxMomTPCOnlyPr && std::abs(trk.tpcNSigmaPr()) < nSigmaTPCOnlyPr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPr(), trk.tofNSigmaPr(), trk.itsNSigmaPr(), trk.tpcExpSignalPr(trk.tpcSignal()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
        } else if (trk.tpcInnerParam() > maxMomTPCOnlyPr && std::abs(trk.tofNSigmaPr()) < nSigmaTOF_TPCTOF_Pr && std::abs(trk.tpcNSigmaPr()) < nSigmaTPC_TPCTOF_Pr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPr(), trk.tofNSigmaPr(), trk.itsNSigmaPr(), trk.tpcExpSignalPr(trk.tpcSignal()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
        }
        /// Fill tree for kaons
        if (trk.tpcInnerParam() < maxMomHardCutOnlyKa && trk.tpcInnerParam() <= maxMomTPCOnlyKa && std::abs(trk.tpcNSigmaKa()) < nSigmaTPCOnlyKa && downsampleTsalisCharged(trk.pt(), downsamplingTsalisKaons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Kaon])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaKa(), trk.tofNSigmaKa(), trk.itsNSigmaKa(), trk.tpcExpSignalKa(trk.tpcSignal()), o2::track::PID::Kaon, runnumber, dwnSmplFactor_Ka, hadronicRate);
        } else if (trk.tpcInnerParam() < maxMomHardCutOnlyKa && trk.tpcInnerParam() > maxMomTPCOnlyKa && std::abs(trk.tofNSigmaKa()) < nSigmaTOF_TPCTOF_Ka && std::abs(trk.tpcNSigmaKa()) < nSigmaTPC_TPCTOF_Ka && downsampleTsalisCharged(trk.pt(), downsamplingTsalisKaons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Kaon])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaKa(), trk.tofNSigmaKa(), trk.itsNSigmaKa(), trk.tpcExpSignalKa(trk.tpcSignal()), o2::track::PID::Kaon, runnumber, dwnSmplFactor_Ka, hadronicRate);
        }
        /// Fill tree pions
        if (trk.tpcInnerParam() <= maxMomTPCOnlyPi && std::abs(trk.tpcNSigmaPi()) < nSigmaTPCOnlyPi && downsampleTsalisCharged(trk.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPi(), trk.tofNSigmaPi(), trk.itsNSigmaPi(), trk.tpcExpSignalPi(trk.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
        } else if (trk.tpcInnerParam() > maxMomTPCOnlyPi && std::abs(trk.tofNSigmaPi()) < nSigmaTOF_TPCTOF_Pi && std::abs(trk.tpcNSigmaPi()) < nSigmaTPC_TPCTOF_Pi && downsampleTsalisCharged(trk.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPi(), trk.tofNSigmaPi(), trk.itsNSigmaPi(), trk.tpcExpSignalPi(trk.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
        }
      } /// Loop tracks
    }
  } /// process
  PROCESS_SWITCH(TreeWriterTPCTOF, processWithdEdxTrQA, "Samples for PID with TrackQA info", false);

  Preslice<TrksWithDEdxCorrection> perCollisionTracksWithCorrecteddEdx = aod::track::collisionId;
  void processWithdEdxTrQAWithCorrecteddEdx(Colls const& collisions, TrksWithDEdxCorrection const& myTracks, aod::BCsWithTimestamps const&, aod::TracksQAVersion const& tracksQA)
  {
    std::vector<int64_t> labelTrack2TrackQA;
    labelTrack2TrackQA.clear();
    labelTrack2TrackQA.resize(myTracks.size(), -1);
    for (const auto& trackQA : tracksQA) {
      int64_t trackId = trackQA.trackId();
      int64_t trackQAIndex = trackQA.globalIndex();
      labelTrack2TrackQA[trackId] = trackQAIndex;
    }
    for (const auto& collision : collisions) {
      auto tracks = myTracks.sliceBy(perCollisionTracksWithCorrecteddEdx, collision.globalIndex());
      auto tracksWithITSPid = soa::Attach<TrksWithDEdxCorrection,
                                          aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi,
                                          aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe,
                                          aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe, aod::pidits::ITSNSigmaAl>(tracks);
      /// Check event selection
      if (!isEventSelected(collision, tracks)) {
        continue;
      }
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      const int runnumber = bc.runNumber();
      float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * 1.e-3;
      rowTPCTOFTreeWithTrkQA.reserve(tracks.size());
      for (auto const& trk : tracksWithITSPid) {
        if (!((trackSelection == 0) ||
              ((trackSelection == 1) && trk.isGlobalTrack()) ||
              ((trackSelection == 2) && trk.isGlobalTrackWoPtEta()) ||
              ((trackSelection == 3) && trk.isGlobalTrackWoDCA()) ||
              ((trackSelection == 4) && trk.isQualityTrack()) ||
              ((trackSelection == 5) && trk.isInAcceptanceTrack()))) {
          continue;
        }
        // get the corresponding trackQA using labelTracks2TracKQA and get variables of interest
        aod::TracksQA trackQA;
        bool existTrkQA;
        if (labelTrack2TrackQA[trk.globalIndex()] != -1) {
          trackQA = tracksQA.iteratorAt(labelTrack2TrackQA[trk.globalIndex()]);
          existTrkQA = true;
        } else {
          trackQA = tracksQA.iteratorAt(0);
          existTrkQA = false;
        }
        /// Fill tree for tritons
        if (trk.tpcInnerParam() < maxMomHardCutOnlyTr && trk.tpcInnerParam() <= maxMomTPCOnlyTr && std::abs(trk.tpcNSigmaTr()) < nSigmaTPCOnlyTr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Triton])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaTr(), trk.tofNSigmaTr(), trk.itsNSigmaTr(), trk.tpcExpSignalTr(trk.tpcSignalCorrected()), o2::track::PID::Triton, runnumber, dwnSmplFactor_Tr, hadronicRate);
        } else if (trk.tpcInnerParam() < maxMomHardCutOnlyTr && trk.tpcInnerParam() > maxMomTPCOnlyTr && std::abs(trk.tofNSigmaTr()) < nSigmaTOF_TPCTOF_Tr && std::abs(trk.tpcNSigmaTr()) < nSigmaTPC_TPCTOF_Tr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Triton])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaTr(), trk.tofNSigmaTr(), trk.itsNSigmaTr(), trk.tpcExpSignalTr(trk.tpcSignalCorrected()), o2::track::PID::Triton, runnumber, dwnSmplFactor_Tr, hadronicRate);
        }
        /// Fill tree for deuterons
        if (trk.tpcInnerParam() < maxMomHardCutOnlyDe && trk.tpcInnerParam() <= maxMomTPCOnlyDe && std::abs(trk.tpcNSigmaDe()) < nSigmaTPCOnlyDe && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Deuteron])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaDe(), trk.tofNSigmaDe(), trk.itsNSigmaDe(), trk.tpcExpSignalDe(trk.tpcSignalCorrected()), o2::track::PID::Deuteron, runnumber, dwnSmplFactor_De, hadronicRate);
        } else if (trk.tpcInnerParam() < maxMomHardCutOnlyDe && trk.tpcInnerParam() > maxMomTPCOnlyDe && std::abs(trk.tofNSigmaDe()) < nSigmaTOF_TPCTOF_De && std::abs(trk.tpcNSigmaDe()) < nSigmaTPC_TPCTOF_De && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Deuteron])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaDe(), trk.tofNSigmaDe(), trk.itsNSigmaDe(), trk.tpcExpSignalDe(trk.tpcSignalCorrected()), o2::track::PID::Deuteron, runnumber, dwnSmplFactor_De, hadronicRate);
        }
        /// Fill tree for protons
        if (trk.tpcInnerParam() <= maxMomTPCOnlyPr && std::abs(trk.tpcNSigmaPr()) < nSigmaTPCOnlyPr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPr(), trk.tofNSigmaPr(), trk.itsNSigmaPr(), trk.tpcExpSignalPr(trk.tpcSignalCorrected()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
        } else if (trk.tpcInnerParam() > maxMomTPCOnlyPr && std::abs(trk.tofNSigmaPr()) < nSigmaTOF_TPCTOF_Pr && std::abs(trk.tpcNSigmaPr()) < nSigmaTPC_TPCTOF_Pr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPr(), trk.tofNSigmaPr(), trk.itsNSigmaPr(), trk.tpcExpSignalPr(trk.tpcSignalCorrected()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate);
        }
        /// Fill tree for kaons
        if (trk.tpcInnerParam() < maxMomHardCutOnlyKa && trk.tpcInnerParam() <= maxMomTPCOnlyKa && std::abs(trk.tpcNSigmaKa()) < nSigmaTPCOnlyKa && downsampleTsalisCharged(trk.pt(), downsamplingTsalisKaons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Kaon])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaKa(), trk.tofNSigmaKa(), trk.itsNSigmaKa(), trk.tpcExpSignalKa(trk.tpcSignalCorrected()), o2::track::PID::Kaon, runnumber, dwnSmplFactor_Ka, hadronicRate);
        } else if (trk.tpcInnerParam() < maxMomHardCutOnlyKa && trk.tpcInnerParam() > maxMomTPCOnlyKa && std::abs(trk.tofNSigmaKa()) < nSigmaTOF_TPCTOF_Ka && std::abs(trk.tpcNSigmaKa()) < nSigmaTPC_TPCTOF_Ka && downsampleTsalisCharged(trk.pt(), downsamplingTsalisKaons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Kaon])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaKa(), trk.tofNSigmaKa(), trk.itsNSigmaKa(), trk.tpcExpSignalKa(trk.tpcSignalCorrected()), o2::track::PID::Kaon, runnumber, dwnSmplFactor_Ka, hadronicRate);
        }
        /// Fill tree pions
        if (trk.tpcInnerParam() <= maxMomTPCOnlyPi && std::abs(trk.tpcNSigmaPi()) < nSigmaTPCOnlyPi && downsampleTsalisCharged(trk.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPi(), trk.tofNSigmaPi(), trk.itsNSigmaPi(), trk.tpcExpSignalPi(trk.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
        } else if (trk.tpcInnerParam() > maxMomTPCOnlyPi && std::abs(trk.tofNSigmaPi()) < nSigmaTOF_TPCTOF_Pi && std::abs(trk.tpcNSigmaPi()) < nSigmaTPC_TPCTOF_Pi && downsampleTsalisCharged(trk.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion])) {
          fillSkimmedTPCTOFTableWithdEdxTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPi(), trk.tofNSigmaPi(), trk.itsNSigmaPi(), trk.tpcExpSignalPi(trk.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate);
        }
      } /// Loop tracks
    }
  } /// process
  PROCESS_SWITCH(TreeWriterTPCTOF, processWithdEdxTrQAWithCorrecteddEdx, "Samples for PID with TrackQA info with corrected dEdx", false);

  void processWithTrQA(Colls const& collisions, Trks const& myTracks, MyBCTable const&, aod::TracksQAVersion const& tracksQA)
  {
    std::vector<int64_t> labelTrack2TrackQA;
    labelTrack2TrackQA.clear();
    labelTrack2TrackQA.resize(myTracks.size(), -1);
    for (const auto& trackQA : tracksQA) {
      int64_t trackId = trackQA.trackId();
      int64_t trackQAIndex = trackQA.globalIndex();
      labelTrack2TrackQA[trackId] = trackQAIndex;
    }
    for (const auto& collision : collisions) {
      auto tracks = myTracks.sliceBy(perCollisionTracks, collision.globalIndex());
      auto tracksWithITSPid = soa::Attach<Trks,
                                          aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi,
                                          aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe,
                                          aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe, aod::pidits::ITSNSigmaAl>(tracks);
      /// Check event selection
      if (!isEventSelected(collision, tracks)) {
        continue;
      }
      auto bc = collision.bc_as<MyBCTable>();
      const int runnumber = bc.runNumber();
      float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * 1.e-3;
      const int bcGlobalIndex = bc.globalIndex();
      const int bcTimeFrameId = bc.tfId();
      const int bcBcInTimeFrame = bc.bcInTF();
      rowTPCTOFTreeWithTrkQA.reserve(tracks.size());
      for (auto const& trk : tracksWithITSPid) {
        if (!((trackSelection == 0) ||
              ((trackSelection == 1) && trk.isGlobalTrack()) ||
              ((trackSelection == 2) && trk.isGlobalTrackWoPtEta()) ||
              ((trackSelection == 3) && trk.isGlobalTrackWoDCA()) ||
              ((trackSelection == 4) && trk.isQualityTrack()) ||
              ((trackSelection == 5) && trk.isInAcceptanceTrack()))) {
          continue;
        }
        // get the corresponding trackQA using labelTracks2TracKQA and get variables of interest
        aod::TracksQA trackQA;
        bool existTrkQA;
        if (labelTrack2TrackQA[trk.globalIndex()] != -1) {
          trackQA = tracksQA.iteratorAt(labelTrack2TrackQA[trk.globalIndex()]);
          existTrkQA = true;
        } else {
          trackQA = tracksQA.iteratorAt(0);
          existTrkQA = false;
        }
        /// Fill tree for tritons
        if (trk.tpcInnerParam() < maxMomHardCutOnlyTr && trk.tpcInnerParam() <= maxMomTPCOnlyTr && std::abs(trk.tpcNSigmaTr()) < nSigmaTPCOnlyTr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Triton])) {
          fillSkimmedTPCTOFTableWithTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaTr(), trk.tofNSigmaTr(), trk.itsNSigmaTr(), trk.tpcExpSignalTr(trk.tpcSignal()), o2::track::PID::Triton, runnumber, dwnSmplFactor_Tr, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        } else if (trk.tpcInnerParam() < maxMomHardCutOnlyTr && trk.tpcInnerParam() > maxMomTPCOnlyTr && std::abs(trk.tofNSigmaTr()) < nSigmaTOF_TPCTOF_Tr && std::abs(trk.tpcNSigmaTr()) < nSigmaTPC_TPCTOF_Tr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Triton])) {
          fillSkimmedTPCTOFTableWithTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaTr(), trk.tofNSigmaTr(), trk.itsNSigmaTr(), trk.tpcExpSignalTr(trk.tpcSignal()), o2::track::PID::Triton, runnumber, dwnSmplFactor_Tr, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        }
        /// Fill tree for deuterons
        if (trk.tpcInnerParam() < maxMomHardCutOnlyDe && trk.tpcInnerParam() <= maxMomTPCOnlyDe && std::abs(trk.tpcNSigmaDe()) < nSigmaTPCOnlyDe && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Deuteron])) {
          fillSkimmedTPCTOFTableWithTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaDe(), trk.tofNSigmaDe(), trk.itsNSigmaDe(), trk.tpcExpSignalDe(trk.tpcSignal()), o2::track::PID::Deuteron, runnumber, dwnSmplFactor_De, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        } else if (trk.tpcInnerParam() < maxMomHardCutOnlyDe && trk.tpcInnerParam() > maxMomTPCOnlyDe && std::abs(trk.tofNSigmaDe()) < nSigmaTOF_TPCTOF_De && std::abs(trk.tpcNSigmaDe()) < nSigmaTPC_TPCTOF_De && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Deuteron])) {
          fillSkimmedTPCTOFTableWithTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaDe(), trk.tofNSigmaDe(), trk.itsNSigmaDe(), trk.tpcExpSignalDe(trk.tpcSignal()), o2::track::PID::Deuteron, runnumber, dwnSmplFactor_De, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        }
        /// Fill tree for protons
        if (trk.tpcInnerParam() <= maxMomTPCOnlyPr && std::abs(trk.tpcNSigmaPr()) < nSigmaTPCOnlyPr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton])) {
          fillSkimmedTPCTOFTableWithTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPr(), trk.tofNSigmaPr(), trk.itsNSigmaPr(), trk.tpcExpSignalPr(trk.tpcSignal()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        } else if (trk.tpcInnerParam() > maxMomTPCOnlyPr && std::abs(trk.tofNSigmaPr()) < nSigmaTOF_TPCTOF_Pr && std::abs(trk.tpcNSigmaPr()) < nSigmaTPC_TPCTOF_Pr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton])) {
          fillSkimmedTPCTOFTableWithTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPr(), trk.tofNSigmaPr(), trk.itsNSigmaPr(), trk.tpcExpSignalPr(trk.tpcSignal()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        }
        /// Fill tree for kaons
        if (trk.tpcInnerParam() < maxMomHardCutOnlyKa && trk.tpcInnerParam() <= maxMomTPCOnlyKa && std::abs(trk.tpcNSigmaKa()) < nSigmaTPCOnlyKa && downsampleTsalisCharged(trk.pt(), downsamplingTsalisKaons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Kaon])) {
          fillSkimmedTPCTOFTableWithTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaKa(), trk.tofNSigmaKa(), trk.itsNSigmaKa(), trk.tpcExpSignalKa(trk.tpcSignal()), o2::track::PID::Kaon, runnumber, dwnSmplFactor_Ka, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        } else if (trk.tpcInnerParam() < maxMomHardCutOnlyKa && trk.tpcInnerParam() > maxMomTPCOnlyKa && std::abs(trk.tofNSigmaKa()) < nSigmaTOF_TPCTOF_Ka && std::abs(trk.tpcNSigmaKa()) < nSigmaTPC_TPCTOF_Ka && downsampleTsalisCharged(trk.pt(), downsamplingTsalisKaons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Kaon])) {
          fillSkimmedTPCTOFTableWithTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaKa(), trk.tofNSigmaKa(), trk.itsNSigmaKa(), trk.tpcExpSignalKa(trk.tpcSignal()), o2::track::PID::Kaon, runnumber, dwnSmplFactor_Ka, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        }
        /// Fill tree pions
        if (trk.tpcInnerParam() <= maxMomTPCOnlyPi && std::abs(trk.tpcNSigmaPi()) < nSigmaTPCOnlyPi && downsampleTsalisCharged(trk.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion])) {
          fillSkimmedTPCTOFTableWithTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPi(), trk.tofNSigmaPi(), trk.itsNSigmaPi(), trk.tpcExpSignalPi(trk.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        } else if (trk.tpcInnerParam() > maxMomTPCOnlyPi && std::abs(trk.tofNSigmaPi()) < nSigmaTOF_TPCTOF_Pi && std::abs(trk.tpcNSigmaPi()) < nSigmaTPC_TPCTOF_Pi && downsampleTsalisCharged(trk.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion])) {
          fillSkimmedTPCTOFTableWithTrkQA(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPi(), trk.tofNSigmaPi(), trk.itsNSigmaPi(), trk.tpcExpSignalPi(trk.tpcSignal()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        }
      } /// Loop tracks
    }
  } /// process
  PROCESS_SWITCH(TreeWriterTPCTOF, processWithTrQA, "Samples for PID with TrackQA info", false);

  void processWithTrQAWithCorrecteddEdx(Colls const& collisions, TrksWithDEdxCorrection const& myTracks, MyBCTable const&, aod::TracksQAVersion const& tracksQA)
  {
    std::vector<int64_t> labelTrack2TrackQA;
    labelTrack2TrackQA.clear();
    labelTrack2TrackQA.resize(myTracks.size(), -1);
    for (const auto& trackQA : tracksQA) {
      int64_t trackId = trackQA.trackId();
      int64_t trackQAIndex = trackQA.globalIndex();
      labelTrack2TrackQA[trackId] = trackQAIndex;
    }
    for (const auto& collision : collisions) {
      auto tracks = myTracks.sliceBy(perCollisionTracks, collision.globalIndex());
      auto tracksWithITSPid = soa::Attach<TrksWithDEdxCorrection,
                                          aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi,
                                          aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe,
                                          aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe, aod::pidits::ITSNSigmaAl>(tracks);
      /// Check event selection
      if (!isEventSelected(collision, tracks)) {
        continue;
      }
      auto bc = collision.bc_as<MyBCTable>();
      const int runnumber = bc.runNumber();
      float hadronicRate = mRateFetcher.fetch(ccdb.service, bc.timestamp(), runnumber, irSource) * 1.e-3;
      const int bcGlobalIndex = bc.globalIndex();
      const int bcTimeFrameId = bc.tfId();
      const int bcBcInTimeFrame = bc.bcInTF();
      rowTPCTOFTreeWithTrkQA.reserve(tracks.size());
      for (auto const& trk : tracksWithITSPid) {
        if (!((trackSelection == 0) ||
              ((trackSelection == 1) && trk.isGlobalTrack()) ||
              ((trackSelection == 2) && trk.isGlobalTrackWoPtEta()) ||
              ((trackSelection == 3) && trk.isGlobalTrackWoDCA()) ||
              ((trackSelection == 4) && trk.isQualityTrack()) ||
              ((trackSelection == 5) && trk.isInAcceptanceTrack()))) {
          continue;
        }
        // get the corresponding trackQA using labelTracks2TracKQA and get variables of interest
        aod::TracksQA trackQA;
        bool existTrkQA;
        if (labelTrack2TrackQA[trk.globalIndex()] != -1) {
          trackQA = tracksQA.iteratorAt(labelTrack2TrackQA[trk.globalIndex()]);
          existTrkQA = true;
        } else {
          trackQA = tracksQA.iteratorAt(0);
          existTrkQA = false;
        }
        /// Fill tree for tritons
        if (trk.tpcInnerParam() < maxMomHardCutOnlyTr && trk.tpcInnerParam() <= maxMomTPCOnlyTr && std::abs(trk.tpcNSigmaTr()) < nSigmaTPCOnlyTr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Triton])) {
          fillSkimmedTPCTOFTableWithTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaTr(), trk.tofNSigmaTr(), trk.itsNSigmaTr(), trk.tpcExpSignalTr(trk.tpcSignalCorrected()), o2::track::PID::Triton, runnumber, dwnSmplFactor_Tr, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        } else if (trk.tpcInnerParam() < maxMomHardCutOnlyTr && trk.tpcInnerParam() > maxMomTPCOnlyTr && std::abs(trk.tofNSigmaTr()) < nSigmaTOF_TPCTOF_Tr && std::abs(trk.tpcNSigmaTr()) < nSigmaTPC_TPCTOF_Tr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Triton])) {
          fillSkimmedTPCTOFTableWithTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaTr(), trk.tofNSigmaTr(), trk.itsNSigmaTr(), trk.tpcExpSignalTr(trk.tpcSignalCorrected()), o2::track::PID::Triton, runnumber, dwnSmplFactor_Tr, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        }
        /// Fill tree for deuterons
        if (trk.tpcInnerParam() < maxMomHardCutOnlyDe && trk.tpcInnerParam() <= maxMomTPCOnlyDe && std::abs(trk.tpcNSigmaDe()) < nSigmaTPCOnlyDe && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Deuteron])) {
          fillSkimmedTPCTOFTableWithTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaDe(), trk.tofNSigmaDe(), trk.itsNSigmaDe(), trk.tpcExpSignalDe(trk.tpcSignalCorrected()), o2::track::PID::Deuteron, runnumber, dwnSmplFactor_De, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        } else if (trk.tpcInnerParam() < maxMomHardCutOnlyDe && trk.tpcInnerParam() > maxMomTPCOnlyDe && std::abs(trk.tofNSigmaDe()) < nSigmaTOF_TPCTOF_De && std::abs(trk.tpcNSigmaDe()) < nSigmaTPC_TPCTOF_De && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Deuteron])) {
          fillSkimmedTPCTOFTableWithTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaDe(), trk.tofNSigmaDe(), trk.itsNSigmaDe(), trk.tpcExpSignalDe(trk.tpcSignalCorrected()), o2::track::PID::Deuteron, runnumber, dwnSmplFactor_De, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        }
        /// Fill tree for protons
        if (trk.tpcInnerParam() <= maxMomTPCOnlyPr && std::abs(trk.tpcNSigmaPr()) < nSigmaTPCOnlyPr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton])) {
          fillSkimmedTPCTOFTableWithTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPr(), trk.tofNSigmaPr(), trk.itsNSigmaPr(), trk.tpcExpSignalPr(trk.tpcSignalCorrected()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        } else if (trk.tpcInnerParam() > maxMomTPCOnlyPr && std::abs(trk.tofNSigmaPr()) < nSigmaTOF_TPCTOF_Pr && std::abs(trk.tpcNSigmaPr()) < nSigmaTPC_TPCTOF_Pr && downsampleTsalisCharged(trk.pt(), downsamplingTsalisProtons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Proton])) {
          fillSkimmedTPCTOFTableWithTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPr(), trk.tofNSigmaPr(), trk.itsNSigmaPr(), trk.tpcExpSignalPr(trk.tpcSignalCorrected()), o2::track::PID::Proton, runnumber, dwnSmplFactor_Pr, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        }
        /// Fill tree for kaons
        if (trk.tpcInnerParam() < maxMomHardCutOnlyKa && trk.tpcInnerParam() <= maxMomTPCOnlyKa && std::abs(trk.tpcNSigmaKa()) < nSigmaTPCOnlyKa && downsampleTsalisCharged(trk.pt(), downsamplingTsalisKaons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Kaon])) {
          fillSkimmedTPCTOFTableWithTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaKa(), trk.tofNSigmaKa(), trk.itsNSigmaKa(), trk.tpcExpSignalKa(trk.tpcSignalCorrected()), o2::track::PID::Kaon, runnumber, dwnSmplFactor_Ka, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        } else if (trk.tpcInnerParam() < maxMomHardCutOnlyKa && trk.tpcInnerParam() > maxMomTPCOnlyKa && std::abs(trk.tofNSigmaKa()) < nSigmaTOF_TPCTOF_Ka && std::abs(trk.tpcNSigmaKa()) < nSigmaTPC_TPCTOF_Ka && downsampleTsalisCharged(trk.pt(), downsamplingTsalisKaons, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Kaon])) {
          fillSkimmedTPCTOFTableWithTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaKa(), trk.tofNSigmaKa(), trk.itsNSigmaKa(), trk.tpcExpSignalKa(trk.tpcSignalCorrected()), o2::track::PID::Kaon, runnumber, dwnSmplFactor_Ka, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        }
        /// Fill tree pions
        if (trk.tpcInnerParam() <= maxMomTPCOnlyPi && std::abs(trk.tpcNSigmaPi()) < nSigmaTPCOnlyPi && downsampleTsalisCharged(trk.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion])) {
          fillSkimmedTPCTOFTableWithTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPi(), trk.tofNSigmaPi(), trk.itsNSigmaPi(), trk.tpcExpSignalPi(trk.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        } else if (trk.tpcInnerParam() > maxMomTPCOnlyPi && std::abs(trk.tofNSigmaPi()) < nSigmaTOF_TPCTOF_Pi && std::abs(trk.tpcNSigmaPi()) < nSigmaTPC_TPCTOF_Pi && downsampleTsalisCharged(trk.pt(), downsamplingTsalisPions, sqrtSNN, o2::track::pid_constants::sMasses[o2::track::PID::Pion])) {
          fillSkimmedTPCTOFTableWithTrkQA<true>(trk, trackQA, existTrkQA, collision, trk.tpcNSigmaPi(), trk.tofNSigmaPi(), trk.itsNSigmaPi(), trk.tpcExpSignalPi(trk.tpcSignalCorrected()), o2::track::PID::Pion, runnumber, dwnSmplFactor_Pi, hadronicRate, bcGlobalIndex, bcTimeFrameId, bcBcInTimeFrame);
        }
      } /// Loop tracks
    }
  } /// process
  PROCESS_SWITCH(TreeWriterTPCTOF, processWithTrQAWithCorrecteddEdx, "Samples for PID with TrackQA info with correced dEdx", false);

  void processDummy(Colls const&) {}
  PROCESS_SWITCH(TreeWriterTPCTOF, processDummy, "Dummy function", false);

}; /// struct TreeWriterTPCTOF
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<TreeWriterTPCTOF>(cfgc)};
  workflow.push_back(adaptAnalysisTask<TreeWriterTpcV0>(cfgc));
  return workflow;
}
