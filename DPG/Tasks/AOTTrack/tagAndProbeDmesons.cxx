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

/// \file tagAndProbeDmesons.cxx
/// \brief Task for tracking efficiency studies with tag-and-probe using 3-prong D-meson decays
///
/// \author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN

#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/Core/TrackSelection.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace tagandprobe
{
enum TagChannels : uint8_t {
  DplusToKPiPi = 0,
  DsOrDplusToKKPi,
  DstarPlusToDzeroPi,
  DstarMinusToDzeroBarPi,
  NTagChannels
};

enum TrackTypes : uint8_t {
  GlobalWoDca = 0,
  GlobalWoDcaWoIts,
  GlobalWoDcaWoTpc,
  NTrackTypes
};

static constexpr int nBinsPt = 7;
static constexpr int nCutVars = 6;
static constexpr int nCutVarsDzero = 9;
constexpr float binsPt[nBinsPt + 1] = {0., 1., 2., 4., 6., 10., 20., 1000.};
auto vecBinsPt = std::vector<float>{binsPt, binsPt + nBinsPt + 1};

// default values for the cuts
constexpr float cuts[nBinsPt][nCutVars] = {{0.1f, 1.5f, 0.01f, 0.01f, 2.f, 2.f},
                                           {0.1f, 1.5f, 0.01f, 0.01f, 2.f, 2.f},
                                           {0.1f, 1.5f, 0.02f, 0.02f, 2.f, 2.f},
                                           {0.1f, 1.5f, 0.02f, 0.02f, 2.f, 2.f},
                                           {0.1f, 1.5f, 0.04f, 0.04f, 2.f, 2.f},
                                           {0.1f, 1.5f, 0.04f, 0.04f, 2.f, 2.f},
                                           {0.1f, 1.5f, 0.06f, 0.06f, 2.f, 2.f}};

constexpr float cutsDzero[nBinsPt][nCutVarsDzero] = {{1.815f, 1.915f, 0.01f, 0.01f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                                     {1.815f, 1.915f, 0.01f, 0.01f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                                     {1.815f, 1.915f, 0.02f, 0.02f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                                     {1.815f, 1.915f, 0.02f, 0.02f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                                     {1.815f, 1.915f, 0.04f, 0.04f, 2.f, 2.f, 0.f, 0.95f, 0.95f},
                                                     {1.815f, 1.915f, 0.04f, 0.04f, 2.f, 2.f, 0.f, 0.95f, 0.95f},
                                                     {1.815f, 1.915f, 0.06f, 0.06f, 2.f, 2.f, 0.f, 0.95f, 0.95f}};

static const std::vector<std::string> labelsPt{};
static const std::vector<std::string> labelsCutVar = {"minMass", "maxMass", "decayLength", "decayLengthXY", "normDecayLength", "normDecayLengthXY"};
static const std::vector<std::string> labelsCutVarDzero = {"minMass", "maxMass", "decayLength", "decayLengthXY", "normDecayLength", "normDecayLengthXY", "impParProd", "cosPointing", "cosPointingXY"};

DECLARE_SOA_INDEX_COLUMN(Collision, collision);                   //! Collision index
DECLARE_SOA_INDEX_COLUMN_FULL(Track0, track0, int, Tracks, "_0"); //! Index to first track
DECLARE_SOA_INDEX_COLUMN_FULL(Track1, track1, int, Tracks, "_1"); //! Index to second track
} // namespace tagandprobe

DECLARE_SOA_TABLE(PiPiFromDpTags, "AOD", "PIPIFROMDPTAG", //! Table for same sign 2-pion vertices used as tags
                  soa::Index<>,
                  aod::tagandprobe::CollisionId,
                  aod::tagandprobe::Track0Id,
                  aod::tagandprobe::Track1Id);
DECLARE_SOA_TABLE(KaKaFromDspTags, "AOD", "KAKAFROMDSPTAG", //! Table for opposite sign 2-kaon vertices used as tags
                  soa::Index<>,
                  aod::tagandprobe::CollisionId,
                  aod::tagandprobe::Track0Id,
                  aod::tagandprobe::Track1Id,
                  soa::Marker<1>);
DECLARE_SOA_TABLE(PiKaFromDzTags, "AOD", "PIKAFROMDZTAG", //! Table for opposite sign pion(+)-kaon(-) vertices used as tags
                  soa::Index<>,
                  aod::tagandprobe::CollisionId,
                  aod::tagandprobe::Track0Id,
                  aod::tagandprobe::Track1Id,
                  soa::Marker<2>);
DECLARE_SOA_TABLE(KaPiFromDzTags, "AOD", "KAPIFROMDZTAG", //! Table for opposite sign kaon(+)-pion(-) vertices used as tags
                  soa::Index<>,
                  aod::tagandprobe::CollisionId,
                  aod::tagandprobe::Track0Id,
                  aod::tagandprobe::Track1Id,
                  soa::Marker<3>);
} // namespace o2::aod

/// Reconstruction of 2-prong displaced vertices (very good quality and purity)
/// 1) K∓K± for φ from Ds± or D± → φπ± decays
/// 2) π±π± for D± → K∓π±π± decays
/// 3) K∓π± for D0 from D±* → D0π± decays
struct TagTwoProngDisplacedVertices {

  Produces<aod::PiPiFromDpTags> tagPiPiTable;
  Produces<aod::KaKaFromDspTags> tagKaKaTable;
  Produces<aod::KaPiFromDzTags> tagKaPiTable;
  Produces<aod::PiKaFromDzTags> tagPiKaTable;
  SliceCache cache;

  Configurable<bool> applyTofPid{"applyTofPid", true, "flag to enable TOF PID selection"};
  Configurable<float> trackNumSigmaTof{"trackNumSigmaTof", 3.f, "number of sigma for TOF PID compatibility"};
  Configurable<float> trackNumSigmaTpc{"trackNumSigmaTpc", 3.f, "number of sigma for TOF PID compatibility"};
  Configurable<float> trackDcaXyMin{"trackDcaXyMin", 0.002f, "minimum DCAxy for tracks with pT < 2 GeV/c"};
  Configurable<float> trackPtMin{"trackPtMin", 0.4f, "minimum track pT"};

  Configurable<std::vector<float>> binsPtPiPiFromDplus{"binsPtPiPiFromDplus", std::vector<float>{aod::tagandprobe::vecBinsPt}, "pT bin limits for pipi pairs from D+ decays"};
  Configurable<std::vector<float>> binsKaKaFromDsOrDplus{"binsKaKaFromDsOrDplus", std::vector<float>{aod::tagandprobe::vecBinsPt}, "pT bin limits for KK pairs from Ds or D+ decays"};
  Configurable<std::vector<float>> binsPtDzeroFromDstar{"binsPtDzeroFromDstar", std::vector<float>{aod::tagandprobe::vecBinsPt}, "pT bin limits for Kpi pairs from D0 <- D*+ decays"};

  Configurable<LabeledArray<float>> cutsPiPiFromDplus{"cutsPiPiFromDplus", {aod::tagandprobe::cuts[0], aod::tagandprobe::nBinsPt, aod::tagandprobe::nCutVars, aod::tagandprobe::labelsPt, aod::tagandprobe::labelsCutVar}, "Selections for pipi pairs from D+ decays"};
  Configurable<LabeledArray<float>> cutsKaKaFromDsOrDplus{"cutsKaKaFromDsOrDplus", {aod::tagandprobe::cuts[0], aod::tagandprobe::nBinsPt, aod::tagandprobe::nCutVars, aod::tagandprobe::labelsPt, aod::tagandprobe::labelsCutVar}, "Selections for KK pairs from Ds or D+ decays"};
  Configurable<LabeledArray<float>> cutsDzeroFromDstar{"cutsDzeroFromDstar", {aod::tagandprobe::cutsDzero[0], aod::tagandprobe::nBinsPt, aod::tagandprobe::nCutVarsDzero, aod::tagandprobe::labelsPt, aod::tagandprobe::labelsCutVarDzero}, "Selections for Kpi pairs from D0 <- D*+ decays"};

  using TracksWithSelAndDca = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksDCA, aod::TracksExtra, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using CollisionsWithEvSel = soa::Join<aod::Collisions, aod::EvSels>;

  Filter evSelFilter = aod::evsel::sel8 == true;                                                                                                              // simple event selection
  Filter collisionFilter = nabs(aod::collision::posZ) < 10.f;                                                                                                 // simple event selection
  Filter trackFilter = requireGlobalTrackWoDCAInFilter() && aod::track::pt > trackPtMin && (nabs(aod::track::dcaXY) > trackDcaXyMin || aod::track::pt > 2.f); // for the tag, we only consider global tracks with large dcaXY (low pT only)
  using TracksWithSelAndDcaFiltered = soa::Filtered<TracksWithSelAndDca>;
  using CollisionsFiltered = soa::Filtered<CollisionsWithEvSel>;

  // in the partition we only apply TPC PID
  Preslice<TracksWithSelAndDcaFiltered> perCollision = aod::track::collisionId;
  Partition<TracksWithSelAndDcaFiltered> positivePions = aod::track::signed1Pt > 0.f && nabs(aod::pidtpc::tpcNSigmaPi) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaFiltered> negativePions = aod::track::signed1Pt < 0.f && nabs(aod::pidtpc::tpcNSigmaPi) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaFiltered> positiveKaons = aod::track::signed1Pt > 0.f && nabs(aod::pidtpc::tpcNSigmaKa) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaFiltered> negativeKaons = aod::track::signed1Pt < 0.f && nabs(aod::pidtpc::tpcNSigmaKa) < trackNumSigmaTpc;

  ccdb::CcdbApi ccdbApi;
  Service<ccdb::BasicCCDBManager> ccdb;
  vertexing::DCAFitterN<2> vertexer;
  int runNumber{0};

  std::array<std::array<double, 2>, aod::tagandprobe::TagChannels::NTagChannels> masses = {std::array{constants::physics::MassPionCharged, constants::physics::MassPionCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassKaonCharged},
                                                                                           std::array{constants::physics::MassPionCharged, constants::physics::MassKaonCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassPionCharged}};

  std::array<LabeledArray<float>, aod::tagandprobe::TagChannels::NTagChannels> topologicalCuts{};
  std::array<std::vector<float>, aod::tagandprobe::TagChannels::NTagChannels> ptBinsForTopologicalCuts{};

  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    std::string ccdbUrl = "http://alice-ccdb.cern.ch";
    ccdb->setURL(ccdbUrl.data());
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdbApi.init(ccdbUrl.data());

    vertexer.setPropagateToPCA(true);
    vertexer.setMaxR(200.f);
    vertexer.setMaxDZIni(4.f);
    vertexer.setMinParamChange(1.e-3);
    vertexer.setMinRelChi2Change(0.9f);
    vertexer.setUseAbsDCA(false);

    topologicalCuts = {cutsPiPiFromDplus, cutsKaKaFromDsOrDplus, cutsDzeroFromDstar, cutsDzeroFromDstar};
    ptBinsForTopologicalCuts = {binsPtPiPiFromDplus, binsKaKaFromDsOrDplus, binsPtDzeroFromDstar, binsPtDzeroFromDstar};

    const AxisSpec axisPt{250, 0.f, 50.f};
    const AxisSpec axisMassPiPi{280, 0.1f, 1.5f};
    const AxisSpec axisMassKaKa{200, constants::physics::MassPhi - 0.04f, constants::physics::MassPhi + 0.04f};
    const AxisSpec axisMassKaPi{200, constants::physics::MassD0 - 0.05f, constants::physics::MassD0 + 0.05f};

    if (doprocessPiPiFromDplus) {
      registry.add<TH2>("hMassPiPiVsPt", ";#it{p}_{T}(#pi#pi) (GeV/#it{c}); #it{M}(#pi#pi) (GeV/#it{c}^{2})", HistType::kTH2D, {axisPt, axisMassPiPi});
    }
    if (doprocessKaKaFromDsOrDplus) {
      registry.add<TH2>("hMassKaKaVsPt", ";#it{p}_{T}(KK) (GeV/#it{c}); #it{M}(KK) (GeV/#it{c}^{2})", HistType::kTH2D, {axisPt, axisMassKaKa});
    }
    if (doprocessKaPiFromDstar) {
      registry.add<TH2>("hMassKaPiVsPt", ";#it{p}_{T}(K#pi) (GeV/#it{c}); #it{M}(K#pi) (GeV/#it{c}^{2})", HistType::kTH2D, {axisPt, axisMassKaPi});
    }
  }

  /// Finds pT bin in an array.
  /// \param bins  array of pT bins
  /// \param value  pT
  /// \return index of the pT bin
  template <typename T1, typename T2>
  int findBin(T1 const& binsPt, T2 value)
  {
    if (value < binsPt->front()) {
      return -1;
    }
    if (value >= binsPt->back()) {
      return -1;
    }
    return std::distance(binsPt->begin(), std::upper_bound(binsPt->begin(), binsPt->end(), value)) - 1;
  }

  template <typename Pvec>
  bool isSelectedInvariantMass(const Pvec& pVecTrackFirst,
                               const Pvec& pVecTrackSecond,
                               const uint8_t channel,
                               float& invMass2)
  {
    auto arrMomentum = std::array{pVecTrackFirst, pVecTrackSecond};
    auto pVec = RecoDecay::pVec(pVecTrackFirst, pVecTrackSecond);
    auto ptBin = findBin(&ptBinsForTopologicalCuts[channel], RecoDecay::pt(pVec));
    if (ptBin == -1) {
      return false;
    }

    auto invMassMin = topologicalCuts[channel].get(ptBin, 0u);
    auto invMassMax = topologicalCuts[channel].get(ptBin, 1u);
    invMass2 = RecoDecay::m2(arrMomentum, masses[channel]);
    if (invMass2 > invMassMax * invMassMax || invMass2 < invMassMin * invMassMin) {
      return false;
    }
    return true;
  }

  template <typename TTrack>
  bool isSelectedPidTof(const TTrack& track,
                        const uint8_t channel)
  {
    if (!track.hasTOF()) { // TOF not forced anyway
      return true;
    }

    switch (channel) {
      case aod::tagandprobe::TagChannels::DplusToKPiPi: {
        if (std::abs(track.tofNSigmaPi()) < trackNumSigmaTof) {
          return true;
        }
      }
      case aod::tagandprobe::TagChannels::DsOrDplusToKKPi: {
        if (std::abs(track.tofNSigmaKa()) < trackNumSigmaTof) {
          return true;
        }
      }
      case aod::tagandprobe::TagChannels::DstarPlusToDzeroPi: {
        if ((track.signed1Pt() > 0 && std::abs(track.tofNSigmaPi()) < trackNumSigmaTof) || (track.signed1Pt() < 0 && std::abs(track.tofNSigmaKa()) < trackNumSigmaTof)) {
          return true;
        }
      }
      case aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi: {
        if ((track.signed1Pt() < 0 && std::abs(track.tofNSigmaPi()) < trackNumSigmaTof) || (track.signed1Pt() > 0 && std::abs(track.tofNSigmaKa()) < trackNumSigmaTof)) {
          return true;
        }
      }
    }
    return false;
  }

  template <typename PV, typename SV, typename CovMatSV, typename PVec>
  bool isSelectedTopology(const PV& primVtx,
                          const SV& secVtx,
                          const CovMatSV& covMatrixSecVtx,
                          const PVec& pVec,
                          std::array<float, 2>& trackDcaXy,
                          const uint8_t channel)
  {

    auto ptBin = findBin(&ptBinsForTopologicalCuts[channel], RecoDecay::pt(pVec));
    if (ptBin == -1) {
      return false;
    }

    std::array<float, 3> pvCoord = {primVtx.getX(), primVtx.getY(), primVtx.getZ()};
    auto decLen = RecoDecay::distance(pvCoord, secVtx);
    if (decLen < topologicalCuts[channel].get(ptBin, 2u)) {
      return false;
    }

    auto covMatrixPV = primVtx.getCov();

    auto decLenXy = RecoDecay::distanceXY(pvCoord, secVtx);
    if (decLenXy < topologicalCuts[channel].get(ptBin, 3u)) {
      return false;
    }

    float phi, theta;
    getPointDirection(pvCoord, secVtx, phi, theta);
    auto errorDecLen = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSecVtx, phi, theta));
    if (decLen / errorDecLen < topologicalCuts[channel].get(ptBin, 4u)) {
      return false;
    }

    auto errorDecLenXy = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.f) + getRotatedCovMatrixXX(covMatrixSecVtx, phi, 0.f));
    if (decLenXy / errorDecLenXy < topologicalCuts[channel].get(ptBin, 5u)) {
      return false;
    }

    // only for D0 meson
    if (channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi || channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi) {
      if (trackDcaXy[0] * trackDcaXy[1] > topologicalCuts[channel].get(ptBin, 6u)) {
        return false;
      }
      auto cpa = RecoDecay::cpa(pvCoord, secVtx, pVec);
      if (cpa < topologicalCuts[channel].get(ptBin, 7u)) {
        return false;
      }
      auto cpaXy = RecoDecay::cpaXY(pvCoord, secVtx, pVec);
      if (cpaXy < topologicalCuts[channel].get(ptBin, 8u)) {
        return false;
      }
    }

    return true;
  }

  template <typename CCollision, typename TTracks>
  void computeCombinatorialSameCharge(CCollision const& collision,
                                      TTracks const& tracks, // pool of tracks
                                      const uint8_t channel,
                                      float& bz)
  {
    for (auto trackFirst = tracks.begin(); trackFirst != tracks.end(); ++trackFirst) {

      if (applyTofPid && !isSelectedPidTof(trackFirst, channel)) {
        continue;
      }

      for (auto trackSecond = trackFirst + 1; trackSecond != tracks.end(); ++trackSecond) {

        if (applyTofPid && !isSelectedPidTof(trackSecond, channel)) {
          continue;
        }

        float invMass2{0.f};
        std::array<float, 3> pVecTrackFirst{trackFirst.px(), trackFirst.py(), trackFirst.pz()};
        std::array<float, 3> pVecTrackSecond{trackSecond.px(), trackSecond.py(), trackSecond.pz()};

        if (!isSelectedInvariantMass(pVecTrackFirst, pVecTrackSecond, channel, invMass2)) {
          continue;
        }

        auto trackParCovFirst = getTrackParCov(trackFirst);
        auto trackParCovSecond = getTrackParCov(trackSecond);

        int nVertices{0};
        try {
          nVertices = vertexer.process(trackParCovFirst, trackParCovSecond);
        } catch (...) {
          LOG(error) << "Exception caught in DCA fitter process call!";
          continue;
        }
        if (nVertices == 0) {
          continue;
        }

        auto pVec = RecoDecay::pVec(pVecTrackFirst, pVecTrackSecond);
        auto primVtx = getPrimaryVertex(collision);
        const auto& secVtx = vertexer.getPCACandidate();
        const auto& covMatrixPCA = vertexer.calcPCACovMatrixFlat();
        std::array<float, 2> trackDcaXy{trackFirst.dcaXY(), trackSecond.dcaXY()};
        if (!isSelectedTopology(primVtx, secVtx, covMatrixPCA, pVec, trackDcaXy, channel)) {
          continue;
        }

        registry.fill(HIST("hMassPiPiVsPt"), RecoDecay::pt(pVec), std::sqrt(invMass2)); // only channel with same sign tracks for the moment
        tagPiPiTable(trackFirst.collisionId(), trackFirst.globalIndex(), trackSecond.globalIndex());
      }
    }
  }

  template <typename CCollision, typename TTracks>
  void computeCombinatorialOppositeCharge(CCollision const& collision,
                                          TTracks const& tracksPos,
                                          TTracks const& tracksNeg,
                                          const uint8_t channel,
                                          float& bz)
  {
    for (const auto& trackPos : tracksPos) {

      if (applyTofPid && !isSelectedPidTof(trackPos, channel)) {
        continue;
      }

      for (const auto& trackNeg : tracksNeg) {

        if (applyTofPid && !isSelectedPidTof(trackNeg, channel)) {
          continue;
        }

        float invMass2{0.f};
        std::array<float, 3> pVecTrackPos{trackPos.px(), trackPos.py(), trackPos.pz()};
        std::array<float, 3> pVecTrackNeg{trackNeg.px(), trackNeg.py(), trackNeg.pz()};
        if (!isSelectedInvariantMass(pVecTrackPos, pVecTrackNeg, channel, invMass2)) {
          continue;
        }

        auto trackParCovPos = getTrackParCov(trackPos);
        auto trackParCovNeg = getTrackParCov(trackNeg);

        int nVertices{0};
        try {
          nVertices = vertexer.process(trackParCovPos, trackParCovNeg);
        } catch (...) {
          LOG(error) << "Exception caught in DCA fitter process call!";
          continue;
        }
        if (nVertices == 0) {
          continue;
        }

        auto pVec = RecoDecay::pVec(pVecTrackPos, pVecTrackNeg);
        auto primVtx = getPrimaryVertex(collision);
        const auto& secVtx = vertexer.getPCACandidate();
        const auto& covMatrixPCA = vertexer.calcPCACovMatrixFlat();
        std::array<float, 2> trackDcaXy{trackPos.dcaXY(), trackNeg.dcaXY()};
        if (!isSelectedTopology(primVtx, secVtx, covMatrixPCA, pVec, trackDcaXy, channel)) {
          continue;
        }

        if (channel == aod::tagandprobe::TagChannels::DsOrDplusToKKPi) {
          registry.fill(HIST("hMassKaKaVsPt"), RecoDecay::pt(pVec), std::sqrt(invMass2));
          tagKaKaTable(trackPos.collisionId(), trackPos.globalIndex(), trackNeg.globalIndex());
        } else if (channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi) {
          registry.fill(HIST("hMassKaPiVsPt"), RecoDecay::pt(pVec), std::sqrt(invMass2));
          tagPiKaTable(trackPos.collisionId(), trackPos.globalIndex(), trackNeg.globalIndex());
        } else if (channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi) {
          registry.fill(HIST("hMassKaPiVsPt"), RecoDecay::pt(pVec), std::sqrt(invMass2));
          tagKaPiTable(trackPos.collisionId(), trackPos.globalIndex(), trackNeg.globalIndex());
        }
      }
    }
  }

  void processPiPiFromDplus(CollisionsFiltered::iterator const& collision,
                            TracksWithSelAndDcaFiltered const& tracks,
                            aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    auto groupPositive = positivePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialSameCharge(collision, groupPositive, aod::tagandprobe::TagChannels::DplusToKPiPi, bz);

    auto groupNegative = negativePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialSameCharge(collision, groupNegative, aod::tagandprobe::TagChannels::DplusToKPiPi, bz);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processPiPiFromDplus, "Process pipi combinatorial to tag pion pairs from D+ decays", true);

  void processKaKaFromDsOrDplus(CollisionsFiltered::iterator const& collision,
                                TracksWithSelAndDcaFiltered const& tracks,
                                aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    auto groupPositive = positiveKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupNegative = negativeKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialOppositeCharge(collision, groupPositive, groupNegative, aod::tagandprobe::TagChannels::DsOrDplusToKKPi, bz);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processKaKaFromDsOrDplus, "Process KK combinatorial to tag kaon pairs from Ds+/D+ decays", false);

  void processKaPiFromDstar(CollisionsFiltered::iterator const& collision,
                            TracksWithSelAndDcaFiltered const& tracks,
                            aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    float bz{0};
    if (runNumber != bc.runNumber()) {
      parameters::GRPMagField* grpo = ccdb->getForTimeStamp<parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo != nullptr) {
        base::Propagator::initFieldFromGRP(grpo);
        bz = base::Propagator::Instance()->getNominalBz();
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      runNumber = bc.runNumber();
    }

    auto groupPionPositive = positivePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupPionNegative = negativePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupKaonPositive = positiveKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupKaonNegative = negativeKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialOppositeCharge(collision, groupPionPositive, groupKaonNegative, aod::tagandprobe::TagChannels::DstarPlusToDzeroPi, bz);
    computeCombinatorialOppositeCharge(collision, groupKaonPositive, groupPionNegative, aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi, bz);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processKaPiFromDstar, "Process Kpi combinatorial to tag D0 from D*+ decays", false);
};

/// Probe third track reconstruction efficiency with different selections
struct ProbeThirdTrack {

  using TracksWithSelAndDca = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksDCA, aod::TracksExtra>;

  Preslice<aod::PiPiFromDpTags> tagsPiPiPerCollision = aod::tagandprobe::collisionId;
  Preslice<aod::KaKaFromDspTags> tagsKaKaPerCollision = aod::tagandprobe::collisionId;
  Preslice<aod::PiKaFromDzTags> tagsPiKaPerCollision = aod::tagandprobe::collisionId;
  Preslice<aod::KaPiFromDzTags> tagsKaPiPerCollision = aod::tagandprobe::collisionId;
  Preslice<aod::TrackAssoc> trackIndicesPerCollision = aod::track_association::collisionId;

  std::array<std::array<double, 3>, aod::tagandprobe::TagChannels::NTagChannels> masses = {std::array{constants::physics::MassPionCharged, constants::physics::MassPionCharged, constants::physics::MassKaonCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassKaonCharged, constants::physics::MassPionCharged},
                                                                                           std::array{constants::physics::MassPionCharged, constants::physics::MassKaonCharged, constants::physics::MassPionCharged},
                                                                                           std::array{constants::physics::MassKaonCharged, constants::physics::MassPionCharged, constants::physics::MassPionCharged}};

  std::array<TrackSelection, aod::tagandprobe::TrackTypes::NTrackTypes> trackSelector{}; // define the track selectors

  std::array<std::array<std::shared_ptr<TH2>, aod::tagandprobe::TrackTypes::NTrackTypes>, aod::tagandprobe::TagChannels::NTagChannels> histos{};
  HistogramRegistry registry{"registry"};

  void init(InitContext&)
  {
    // ITS-TPC tracks (global tracks)
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetTrackType(o2::aod::track::TrackTypeEnum::Track);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetPtRange(0.1f, 1e10f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetEtaRange(-0.8f, 0.8f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetRequireITSRefit(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetRequireTPCRefit(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetRequireGoldenChi2(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetMinNCrossedRowsTPC(70);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetMinNCrossedRowsOverFindableClustersTPC(0.8f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetMaxChi2PerClusterTPC(4.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetRequireHitsInITSLayers(1, {0, 1, 2}); // one hit in any IB layer
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetMaxChi2PerClusterITS(36.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDca].SetMaxDcaZ(2.f);

    // TPC tracks (global tracks without ITS requirements)
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetTrackType(o2::aod::track::TrackTypeEnum::Track);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetPtRange(0.1f, 1e10f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetEtaRange(-0.8f, 0.8f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetRequireTPCRefit(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetRequireGoldenChi2(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetMinNCrossedRowsTPC(70);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetMinNCrossedRowsOverFindableClustersTPC(0.8f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoIts].SetMaxChi2PerClusterTPC(4.f);

    // ITS tracks (global tracks without TPC requirements)
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetTrackType(o2::aod::track::TrackTypeEnum::Track);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetPtRange(0.1f, 1e10f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetEtaRange(-0.8f, 0.8f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetRequireITSRefit(true);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetRequireHitsInITSLayers(1, {0, 1, 2}); // one hit in any SPD layer
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetMaxChi2PerClusterITS(36.f);
    trackSelector[aod::tagandprobe::TrackTypes::GlobalWoDcaWoTpc].SetMaxDcaZ(2.f);

    const AxisSpec axisPt{250, 0.f, 25.f};
    std::array<AxisSpec, aod::tagandprobe::TagChannels::NTagChannels> axisMass = {AxisSpec{450, 1.65f, 2.10f}, AxisSpec{450, 1.65f, 2.10f}, AxisSpec{360, 0.f, 0.18f}, AxisSpec{360, 0.f, 0.18f}};

    std::string trackTypes[aod::tagandprobe::TrackTypes::NTrackTypes] = {"ItsTpc", "Tpc", "Its"};
    std::string tagChannels[aod::tagandprobe::TagChannels::NTagChannels] = {"DplusToKPiPi", "DsOrDplusToKKPi", "DstarPlusToDzeroPi", "DstarMinusToDzeroBarPi"};

    for (int iChannel{0}; iChannel < aod::tagandprobe::TagChannels::NTagChannels; ++iChannel) {
      for (int iTrackType{0}; iTrackType < aod::tagandprobe::TrackTypes::NTrackTypes; ++iTrackType) {
        histos[iChannel][iTrackType] = registry.add<TH2>(Form("h%sVsPtProbe_%s", tagChannels[iChannel].data(), trackTypes[iTrackType].data()), ";#it{p}_{T}(probe track) (GeV/#it{c}); #it{M} (GeV/#it{c}^{2})", HistType::kTH2D, {axisPt, axisMass[iChannel]});
      }
    }
  }

  template <typename TTrack>
  float computeInvariantMass(TTrack const& trackFirst, TTrack const& trackSecond, TTrack const& trackThird, const uint8_t channel)
  {
    std::array<float, 3> pVecTrackFirst{trackFirst.px(), trackFirst.py(), trackFirst.pz()};
    std::array<float, 3> pVecTrackSecond{trackSecond.px(), trackSecond.py(), trackSecond.pz()};
    std::array<float, 3> pVecTrackThird{trackThird.px(), trackThird.py(), trackThird.pz()};
    auto arrMomentum = std::array{pVecTrackFirst, pVecTrackSecond, pVecTrackThird};
    float invMass = RecoDecay::m(arrMomentum, masses[channel]);
    if (channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi || channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi) {
      auto arrMomentumDzero = std::array{pVecTrackFirst, pVecTrackSecond};
      auto massesDzeroDau = std::array{masses[channel][0], masses[channel][1]};
      float invMassDzero = RecoDecay::m(arrMomentumDzero, massesDzeroDau);
      invMass -= invMassDzero;
    }

    return invMass;
  }

  template <typename TTrackIndices, typename TTrack, typename TTracks>
  void loopOverThirdTrack(TTrackIndices const& groupedTrackThirdIndices, TTracks const& tracks, TTrack const& trackFirst, TTrack const& trackSecond, const uint8_t channel)
  {
    for (const auto& trackIndex : groupedTrackThirdIndices) {
      auto trackThird = trackIndex.template track_as<TTracks>();
      if (trackThird.globalIndex() == trackFirst.globalIndex() || trackThird.globalIndex() == trackSecond.globalIndex()) {
        continue;
      }
      if (channel == aod::tagandprobe::TagChannels::DplusToKPiPi && trackThird.signed1Pt() * trackFirst.signed1Pt() > 0.) { // must be opposite sign
        continue;
      }
      if (channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi && trackThird.signed1Pt() < 0.) { // must be positive
        continue;
      }
      if (channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi && trackThird.signed1Pt() > 0.) { // must be negative
        continue;
      }
      auto ptTrackThird = trackThird.pt();
      auto invMass = computeInvariantMass(trackFirst, trackSecond, trackThird, channel);
      if ((channel == aod::tagandprobe::TagChannels::DstarPlusToDzeroPi || channel == aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi) && invMass > 0.18f) {
        continue;
      } else if (invMass < 1.65f || invMass > 2.10f) {
        continue;
      }
      for (int iTrackType{0}; iTrackType < aod::tagandprobe::TrackTypes::NTrackTypes; ++iTrackType) {
        if (trackSelector[iTrackType].IsSelected(trackThird)) {
          histos[channel][iTrackType]->Fill(ptTrackThird, invMass);
        }
      }
    }
  }

  void processCombinatorialDplusToKaPiPi(aod::Collisions const& collisions,
                                         aod::PiPiFromDpTags const& tagsPiPi,
                                         aod::TrackAssoc const& trackIndices,
                                         TracksWithSelAndDca const& tracks)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      // D+ -> pi+pi+K- and c.c.
      auto groupedTagsPiPi = tagsPiPi.sliceBy(tagsPiPiPerCollision, thisCollId);
      for (const auto& tagPiPi : groupedTagsPiPi) {
        auto trackFirst = tagPiPi.track0_as<TracksWithSelAndDca>();
        auto trackSecond = tagPiPi.track1_as<TracksWithSelAndDca>();
        loopOverThirdTrack(groupedTrackIndices, tracks, trackFirst, trackSecond, aod::tagandprobe::TagChannels::DplusToKPiPi);
      }
    }
  }
  PROCESS_SWITCH(ProbeThirdTrack, processCombinatorialDplusToKaPiPi, "Process combinatorial of tagged 2-pion vertices with additional track", true);

  void processCombinatorialDsToPhiPi(aod::Collisions const& collisions,
                                     aod::KaKaFromDspTags const& tagsKaKa,
                                     aod::TrackAssoc const& trackIndices,
                                     TracksWithSelAndDca const& tracks)
  {
    for (const auto& collision : collisions) {
      // Ds+/D+ -> phi(->K+K-)pi+ and c.c.
      auto thisCollId = collision.globalIndex();
      auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      auto groupedTagsKaKa = tagsKaKa.sliceBy(tagsKaKaPerCollision, thisCollId);
      for (const auto& tagKaKa : groupedTagsKaKa) {
        auto trackFirst = tagKaKa.track0_as<TracksWithSelAndDca>();
        auto trackSecond = tagKaKa.track1_as<TracksWithSelAndDca>();
        loopOverThirdTrack(groupedTrackIndices, tracks, trackFirst, trackSecond, aod::tagandprobe::TagChannels::DsOrDplusToKKPi);
      }
    }
  }
  PROCESS_SWITCH(ProbeThirdTrack, processCombinatorialDsToPhiPi, "Process combinatorial of tagged 2-kaon (phi) vertices with additional track", true);

  void processCombinatorialDstarToDzeroPi(aod::Collisions const& collisions,
                                          aod::PiKaFromDzTags const& tagsPiKa,
                                          aod::KaPiFromDzTags const& tagsKaPi,
                                          aod::TrackAssoc const& trackIndices,
                                          TracksWithSelAndDca const& tracks)
  {
    for (const auto& collision : collisions) {
      auto thisCollId = collision.globalIndex();
      auto groupedTrackIndices = trackIndices.sliceBy(trackIndicesPerCollision, thisCollId);
      // D*+ -> D0(->pi+K-)pi+
      auto groupedTagsPiKa = tagsPiKa.sliceBy(tagsPiKaPerCollision, thisCollId);
      for (const auto& tagPiKa : groupedTagsPiKa) {
        auto trackFirst = tagPiKa.track0_as<TracksWithSelAndDca>();
        auto trackSecond = tagPiKa.track1_as<TracksWithSelAndDca>();
        loopOverThirdTrack(groupedTrackIndices, tracks, trackFirst, trackSecond, aod::tagandprobe::TagChannels::DstarPlusToDzeroPi);
      }
      // D*- -> D0bar(->K+pi-)pi-
      auto groupedTagsKaPi = tagsKaPi.sliceBy(tagsKaPiPerCollision, thisCollId);
      for (const auto& tagKaPi : groupedTagsKaPi) {
        auto trackFirst = tagKaPi.track0_as<TracksWithSelAndDca>();
        auto trackSecond = tagKaPi.track1_as<TracksWithSelAndDca>();
        loopOverThirdTrack(groupedTrackIndices, tracks, trackFirst, trackSecond, aod::tagandprobe::TagChannels::DstarMinusToDzeroBarPi);
      }
    }
  }
  PROCESS_SWITCH(ProbeThirdTrack, processCombinatorialDstarToDzeroPi, "Process combinatorial of tagged pion-kaon (D0) vertices with additional track", true);

  void processDummy(aod::Collisions const&) {}
  PROCESS_SWITCH(ProbeThirdTrack, processDummy, "Dummy process function that does nothing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<TagTwoProngDisplacedVertices>(cfgc));
  workflow.push_back(adaptAnalysisTask<ProbeThirdTrack>(cfgc));
  return workflow;
}
