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

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod {
  namespace tagandprobe {
    enum TagChannel : uint8_t {
      DplusToKPiPi = 0,
      DsOrDplusToKKPi,
      DstarPlusToDzeroPi,
      DstarMinusToDzeroBarPi,
      NTagChannels
    };

    static constexpr int nBinsPt = 7;
    static constexpr int nCutVars = 4;
    static constexpr int nCutVarsDzero = 7;
    constexpr float binsPt[nBinsPt + 1] = {0., 1., 2., 4., 6., 10., 20., 1000.};
    auto vecBinsPt = std::vector<float>{binsPt, binsPt + nBinsPt + 1};

    // default values for the cuts
    constexpr float cuts[nBinsPt][nCutVars] = {{0.01f, 0.01f, 2.f, 2.f},
                                               {0.01f, 0.01f, 2.f, 2.f},
                                               {0.02f, 0.02f, 2.f, 2.f},
                                               {0.02f, 0.02f, 2.f, 2.f},
                                               {0.04f, 0.04f, 2.f, 2.f},
                                               {0.04f, 0.04f, 2.f, 2.f},
                                               {0.06f, 0.06f, 2.f, 2.f}};

    constexpr float cutsDzero[nBinsPt][nCutVarsDzero] = {{0.01f, 0.01f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                                         {0.01f, 0.01f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                                         {0.02f, 0.02f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                                         {0.02f, 0.02f, 2.f, 2.f, 0.f, 0.90f, 0.90f},
                                                         {0.04f, 0.04f, 2.f, 2.f, 0.f, 0.95f, 0.95f},
                                                         {0.04f, 0.04f, 2.f, 2.f, 0.f, 0.95f, 0.95f},
                                                         {0.06f, 0.06f, 2.f, 2.f, 0.f, 0.95f, 0.95f}};

    static const std::vector<std::string> labelsPt{};
    static const std::vector<std::string> labelsCutVar = {"decayLength", "decayLengthXY", "normDecayLength", "normDecayLengthXY"};
    static const std::vector<std::string> labelsCutVarDzero = {"decayLength", "decayLengthXY", "normDecayLength", "normDecayLengthXY", "impParProd", "cosPointing", "cosPointingXY"};

    DECLARE_SOA_INDEX_COLUMN(Collision, collision);                   //! Collision index
    DECLARE_SOA_INDEX_COLUMN_FULL(Track0, track0, int, Tracks, "_0"); //! Index to first track
    DECLARE_SOA_INDEX_COLUMN_FULL(Track1, track1, int, Tracks, "_1"); //! Index to second track
    DECLARE_SOA_COLUMN(Channel, channel, uint8_t);                    //! tag channel
  }

  DECLARE_SOA_TABLE(Dmeson2VtxTags, "AOD", "DMESON2VTXTAG", //! Table for 2-track vertices used as tag
                    soa::Index<>,
                    tagandprobe::CollisionId,
                    tagandprobe::Track0Id,
                    tagandprobe::Track1Id,
                    tagandprobe::Channel);
}

/// Reconstruction of 2-prong displaced vertices (very good quality and purity)
/// 1) K∓K± for φ from Ds± or D± → φπ± decays
/// 2) π±π± for D± → K∓π±π± decays
/// 3) K∓π± for D0 from D±* → D0π± decays
struct TagTwoProngDisplacedVertices {

  Produces<Dmeson2VtxTags> tagTable;

  Configurable<bool> applyTofPid{"applyTofPid", true, "flag to enable TOF PID selection"};
  Configurable<float> trackNumSigmaTof{"trackNumSigmaTof", 3.f, "number of sigma for TOF PID compatibility"};
  Configurable<float> trackNumSigmaTpc{"trackNumSigmaTpc", 3.f, "number of sigma for TOF PID compatibility"};
  Configurable<float> trackDcaXyMin{"trackDcaXyMin", 0.002f, "minimum DCAxy for tracks with pT < 2 GeV/c"};
  Configurable<float> trackPtMin{"trackPtMin", 0.4f, "minimum track pT"};

  Configurable<std::vector<float>> binsPtPiPiFromDplus{"binsPtPiPiFromDplus", std::vector<float>{tagandprobe::vecBinsPt}, "pT bin limits for pipi pairs from D+ decays"};
  Configurable<std::vector<float>> binsKaKaFromDsOrDplus{"binsKaKaFromDsOrDplus", std::vector<float>{tagandprobe::vecBinsPt}, "pT bin limits for KK pairs from Ds or D+ decays"};
  Configurable<std::vector<float>> binsPtDzeroFromDstar{"binsPtDzeroFromDstar", std::vector<float>{tagandprobe::vecBinsPt}, "pT bin limits for Kpi pairs from D0 <- D*+ decays"};

  Configurable<LabeledArray<float>> cutsPiPiFromDplus{"cutsPiPiFromDplus", {tagandprobe::cuts[0], tagandprobe::nBinsPt, tagandprobe::nCutVars, tagandprobe::labelsPt, tagandprobe::labelsCutVar}, "Selections for pipi pairs from D+ decays"};
  Configurable<LabeledArray<float>> cutsKaKaFromDsOrDplus{"cutsKaKaFromDsOrDplus", {tagandprobe::cuts[0], tagandprobe::nBinsPt, tagandprobe::nCutVars, tagandprobe::labelsPt, tagandprobe::labelsCutVar}, "Selections for KK pairs from Ds or D+ decays"};
  Configurable<LabeledArray<float>> cutsDzeroFromDstar{"cutsDzeroFromDstar", {tagandprobe::cutsDzero[0], tagandprobe::nBinsPt, tagandprobe::nCutVarsDzero, tagandprobe::labelsPt, tagandprobe::labelsCutVarDzero}, "Selections for Kpi pairs from D0 <- D*+ decays"};

  using TracksWithSelAndDca = soa::Join<Tracks, TracksCov, TracksDCA, TracksExtra, TrackSelection, pidTPCFullPi, pidTOFFullPi, pidTPCFullKa, pidTOFFullKa>;
  using CollisionsWithEvSel = soa::Join<Collisions, EvSels>;

  Filter evSelFilter = evsel::sel8 == true; // simple event selection
  Filter collisionFilter = nabs(collision::posZ) < 10.f; // simple event selection
  Filter trackFilter = requireGlobalTrackWoDCAInFilter() && aod::track::pt > trackPtMin && (nabs(aod::track::dcaXY) > trackDcaXyMin || aod::track::pt > 2.f); // for the tag, we only consider global tracks with lareg dcaXY (low pT only)
  using TracksWithSelAndDcaFiltered = soa::Filtered<TracksWithSelAndDca>;
  using CollisionsFiltered = soa::Filtered<CollisionsWithEvSel>;

  // in the partition we only apply TPC PID
  Partition<TracksWithSelAndDcaFiltered> positivePions = aod::track::signed1Pt > 0. && nabs(pidtpc::tpcNSigmaPi) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaFiltered> negativePions = aod::track::signed1Pt < 0. && nabs(pidtpc::tpcNSigmaPi) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaFiltered> positiveKaons = aod::track::signed1Pt > 0. && nabs(pidtpc::tpcNSigmaKa) < trackNumSigmaTpc;
  Partition<TracksWithSelAndDcaFiltered> negativeKaons = aod::track::signed1Pt < 0. && nabs(pidtpc::tpcNSigmaKa) < trackNumSigmaTpc;

  SliceCache cache;
  Preslice<TracksWithSelAndDcaFiltered> perCollision = aod::track::collisionId;

  ccdb::CcdbApi ccdbApi;
  Service<ccdb::BasicCCDBManager> ccdb;
  vertexing::DCAFitterN<2> vertexer;
  int runNumber{0};

  std::array<std::array<double, 2>, tagandprobe::TagChannel::NTagChannels> masses = {std::array{constants::physics::MassPionCharged, constants::physics::MassPionCharged},
                                                                               std::array{constants::physics::MassKaonCharged, constants::physics::MassKaonCharged},
                                                                               std::array{constants::physics::MassPionCharged, constants::physics::MassKaonCharged},
                                                                               std::array{constants::physics::MassKaonCharged, constants::physics::MassPionCharged}};

  std::array<LabeledArray<float>, tagandprobe::TagChannel::NTagChannels> topologicalCuts{};
  std::array<std::vector<float>, tagandprobe::TagChannel::NTagChannels> ptBinsForTopologicalCuts{};

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

  template<typename Pvec>
  bool isSelectedInvariantMass(const Pvec& pVecTrackFirst,
                               const Pvec& pVecTrackSecond,
                               const uint8_t channel,
                               float& invMassMin,
                               float& invMassMax,
                               float& invMass2)
  {
    auto arrMomentum = std::array{pVecTrackFirst, pVecTrackSecond};
    invMass2 = RecoDecay::m2(arrMomentum, masses[channel]);
    if (invMass2 > invMassMax*invMassMax || invMass2 < invMassMin*invMassMin) {
      return false;
    }
    return true;
  }

  template<typename TTrack>
  bool isSelectedPidTof(const TTrack& track,
                        const uint8_t channel)
  {
    if (!track.hasTOF()) { // TOF not forced anyway
      return true;
    }

    switch (channel) {
      case tagandprobe::TagChannel::DplusToKPiPi:
      {
        if (std::abs(track.tofNSigmaPi()) < trackNumSigmaTof) {
          return true;
        }
      }
      case tagandprobe::TagChannel::DsOrDplusToKKPi:
      {
        if (std::abs(track.tofNSigmaKa()) < trackNumSigmaTof) {
          return true;
        }
      }
      case tagandprobe::TagChannel::DstarPlusToDzeroPi:
      {
        if ((track.signed1Pt() > 0 && std::abs(track.tofNSigmaPi()) < trackNumSigmaTof) || (track.signed1Pt() < 0 && std::abs(track.tofNSigmaKa()) < trackNumSigmaTof)) {
          return true;
        }
      }
      case tagandprobe::TagChannel::DstarMinusToDzeroBarPi:
      {
        if ((track.signed1Pt() < 0 && std::abs(track.tofNSigmaPi()) < trackNumSigmaTof) || (track.signed1Pt() > 0 && std::abs(track.tofNSigmaKa()) < trackNumSigmaTof)) {
          return true;
        }
      }
    }
    return false;
  }

  template<typename PV, typename SV, typename CovMatSV, typename PVec>
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
    if (decLen < topologicalCuts[channel].get(ptBin, 0u)) {
      return false;
    }

    auto covMatrixPV = primVtx.getCov();

    auto decLenXy = RecoDecay::distanceXY(pvCoord, secVtx);
    if (decLenXy < topologicalCuts[channel].get(ptBin, 1u)) {
      return false;
    }

    float phi, theta;
    getPointDirection(pvCoord, secVtx, phi, theta);
    auto errorDecLen = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, theta) + getRotatedCovMatrixXX(covMatrixSecVtx, phi, theta));
    if (decLen/errorDecLen < topologicalCuts[channel].get(ptBin, 2u)) {
      return false;
    }

    auto errorDecLenXy = std::sqrt(getRotatedCovMatrixXX(covMatrixPV, phi, 0.f) + getRotatedCovMatrixXX(covMatrixSecVtx, phi, 0.f));
    if (decLenXy/errorDecLenXy < topologicalCuts[channel].get(ptBin, 3u)) {
      return false;
    }

    // only for D0 meson
    if (channel == tagandprobe::TagChannel::DstarPlusToDzeroPi || channel == tagandprobe::TagChannel::DstarMinusToDzeroBarPi) {
      if (trackDcaXy[0] * trackDcaXy[1] > topologicalCuts[channel].get(ptBin, 4u)) {
        return false;
      }
      auto cpa = RecoDecay::cpa(pvCoord, secVtx, pVec);
      if (cpa < topologicalCuts[channel].get(ptBin, 5u)) {
        return false;
      }
      auto cpaXy = RecoDecay::cpaXY(pvCoord, secVtx, pVec);
      if (cpaXy < topologicalCuts[channel].get(ptBin, 6u)) {
        return false;
      }
    }

    return true;
  }

  template<typename Collision, typename TTracks>
  void computeCombinatorialSameCharge(Collision const& collision,
                                      TTracks const& tracks, // pool of tracks
                                      const uint8_t channel,
                                      float& invMassMin,
                                      float& invMassMax,
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
        if (!isSelectedInvariantMass(pVecTrackFirst, pVecTrackSecond, channel, invMassMin, invMassMax, invMass2)) {
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
        tagTable(trackFirst.collisionId(), trackFirst.globalIndex(), trackSecond.globalIndex(), channel);
      }
    }
  }

  template<typename Collision, typename TTracks>
  void computeCombinatorialOppositeCharge(Collision const& collision,
                                          TTracks const& tracksPos,
                                          TTracks const& tracksNeg,
                                          const uint8_t channel,
                                          float& invMassMin,
                                          float& invMassMax,
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
        if (!isSelectedInvariantMass(pVecTrackPos, pVecTrackNeg, channel, invMassMin, invMassMax, invMass2)) {
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

        if (channel == tagandprobe::TagChannel::DsOrDplusToKKPi) {
          registry.fill(HIST("hMassKaKaVsPt"), RecoDecay::pt(pVec), std::sqrt(invMass2));
        } else if (channel == tagandprobe::TagChannel::DstarPlusToDzeroPi || channel == tagandprobe::TagChannel::DstarMinusToDzeroBarPi) {
          registry.fill(HIST("hMassKaPiVsPt"), RecoDecay::pt(pVec), std::sqrt(invMass2));
        }
        tagTable(trackPos.collisionId(), trackPos.globalIndex(), trackNeg.globalIndex(), channel);
      }
    }
  }

  void processPiPiFromDplus(CollisionsFiltered::iterator const& collision,
                            TracksWithSelAndDcaFiltered const& tracks,
                            BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<BCsWithTimestamps>();
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

    // kinematic limits for massPiPi in D+ decays
    float invMassMin{0.1f};
    float invMassMax{1.5f};

    auto groupPositive = positivePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialSameCharge(collision, groupPositive, tagandprobe::TagChannel::DplusToKPiPi, invMassMin, invMassMax, bz);

    auto groupNegative = negativePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialSameCharge(collision, groupNegative, tagandprobe::TagChannel::DplusToKPiPi, invMassMin, invMassMax, bz);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processPiPiFromDplus, "Process pipi combinatorial to tag pion pairs from D+ decays", true);

  void processKaKaFromDsOrDplus(CollisionsFiltered::iterator const& collision,
                                TracksWithSelAndDcaFiltered const& tracks,
                                BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<BCsWithTimestamps>();
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

    // We expect the phi resonance, so we cut around it with a good margin
    float invMassMin{constants::physics::MassPhi - 0.04f};
    float invMassMax{constants::physics::MassPhi + 0.04f};

    auto groupPositive = positiveKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupNegative = negativeKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialOppositeCharge(collision, groupPositive, groupNegative, tagandprobe::TagChannel::DsOrDplusToKKPi, invMassMin, invMassMax, bz);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processKaKaFromDsOrDplus, "Process KK combinatorial to tag kaon pairs from Ds+/D+ decays", true);

  void processKaPiFromDstar(CollisionsFiltered::iterator const& collision,
                            TracksWithSelAndDcaFiltered const& tracks,
                            BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<BCsWithTimestamps>();
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

    // We expect the phi resonance, so we cut around it with a good margin
    float invMassMin{constants::physics::MassD0 - 0.05f};
    float invMassMax{constants::physics::MassD0 + 0.05f};

    auto groupPionPositive = positivePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupPionNegative = negativePions->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupKaonPositive = positiveKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto groupKaonNegative = negativeKaons->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    computeCombinatorialOppositeCharge(collision, groupPionPositive, groupKaonNegative, tagandprobe::TagChannel::DstarPlusToDzeroPi, invMassMin, invMassMax, bz);
    computeCombinatorialOppositeCharge(collision, groupKaonPositive, groupPionNegative, tagandprobe::TagChannel::DstarMinusToDzeroBarPi, invMassMin, invMassMax, bz);
  }
  PROCESS_SWITCH(TagTwoProngDisplacedVertices, processKaPiFromDstar, "Process Kpi combinatorial to tag D0 from D*+ decays", true);
};

// /// Probe third track reconstruction efficiency with different selections
// struct ProbeThirdTrack {

// }

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{};
  workflow.push_back(adaptAnalysisTask<TagTwoProngDisplacedVertices>(cfgc));
  return workflow;
}