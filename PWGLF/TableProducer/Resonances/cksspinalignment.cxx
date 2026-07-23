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

/// \file cksspinalignment.cxx
/// \brief Reduced table producer for K0s and charged pions for later charged K* reconstruction
///
/// \author prottay.das@cern.ch

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFCKSSpinalignmentTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/RCTSelectionFlags.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>

#include <Math/Vector4D.h>

#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;

struct cksspinalignment {

  Produces<aod::KShortpionEvents> kshortpionEvent;
  Produces<aod::KShortTracks> kshortTrack;
  Produces<aod::PionTracks> pionTrack;

  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;

  Configurable<bool> useNoCollInTimeRangeStandard{"useNoCollInTimeRangeStandard", false, "Apply kNoCollInTimeRangeStandard selection bit"};
  Configurable<bool> useGoodITSLayersAll{"useGoodITSLayersAll", true, "Apply kIsGoodITSLayersAll selection bit"};
  Configurable<int> cfgCutOccupancy{"cfgCutOccupancy", 2000, "Occupancy cut"};

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 80.0f, "Accepted maximum centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 0.0f, "Accepted minimum centrality"};

  Configurable<float> cfgCutPt{"cfgCutPt", 0.2f, "Prefilter minimum pT for bachelor pion"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Prefilter maximum eta for bachelor pion"};

  struct : ConfigurableGroup {
    Configurable<bool> itsPIDSelection{"itsPIDSelection", true, "Apply ITS PID for bachelor pion"};
    Configurable<float> lowITSPIDNsigma{"lowITSPIDNsigma", -3.0f, "Lower ITS n-sigma pion cut"};
    Configurable<float> highITSPIDNsigma{"highITSPIDNsigma", 3.0f, "Upper ITS n-sigma pion cut"};
    Configurable<int> itsclusterPiMeson{"itsclusterPiMeson", 5, "Minimum ITS clusters for bachelor pion"};
    Configurable<int> tpcCrossedRowsPiMeson{"tpcCrossedRowsPiMeson", 80, "Minimum TPC crossed rows for bachelor pion"};
    Configurable<float> cutDCAxyPiMeson{"cutDCAxyPiMeson", 0.1f, "Maximum DCAxy for bachelor pion"};
    Configurable<float> cutDCAzPiMeson{"cutDCAzPiMeson", 0.1f, "Maximum DCAz for bachelor pion"};
    Configurable<float> cutEtaPiMeson{"cutEtaPiMeson", 0.8f, "Maximum eta for bachelor pion"};
    Configurable<float> cutPTPiMeson{"cutPTPiMeson", 0.2f, "Minimum pT for bachelor pion"};
    Configurable<bool> usePID{"usePID", true, "Apply pion PID"};
    Configurable<float> nsigmaCutTPCPiMeson{"nsigmaCutTPCPiMeson", 3.0f, "Maximum TPC n-sigma pion cut"};
    Configurable<float> nsigmaCutTOFPiMeson{"nsigmaCutTOFPiMeson", 3.0f, "Maximum TOF n-sigma pion cut"};
    Configurable<float> cutTOFBetaPiMeson{"cutTOFBetaPiMeson", 0.5f, "Minimum TOF beta for bachelor pion"};
    Configurable<float> pSwitchPID{"pSwitchPID", 0.5f, "pT switch for pT-dependent pion PID"};
  } grpPion;

  enum PionPidBits : uint8_t {
    kPID1 = 1u << 0,
    kPID2 = 1u << 1,
    kPID3 = 1u << 2,
    kPID4 = 1u << 3
  };

  Configurable<float> confV0PtMin{"confV0PtMin", 0.0f, "Minimum K0s pT"};
  Configurable<float> confV0PtMax{"confV0PtMax", 1000.0f, "Maximum K0s pT"};
  Configurable<float> confV0Rap{"confV0Rap", 0.8f, "Maximum K0s rapidity"};
  Configurable<float> confV0DCADaughMax{"confV0DCADaughMax", 1.0f, "Maximum DCA between V0 daughters"};
  Configurable<double> confV0CPAMin{"confV0CPAMin", 0.9998, "Minimum K0s CPA"};
  Configurable<float> confV0TranRadV0Min{"confV0TranRadV0Min", 1.5f, "Minimum K0s transverse radius"};
  Configurable<float> confV0TranRadV0Max{"confV0TranRadV0Max", 100.0f, "Maximum K0s transverse radius"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 1.2, "Maximum K0s DCA to PV"};
  Configurable<float> cMinV0DCAPi{"cMinV0DCAPi", 0.05f, "Minimum V0 daughter DCA to PV"};
  Configurable<float> cMaxV0LifeTime{"cMaxV0LifeTime", 50.0f, "Maximum K0s lifetime"};
  Configurable<float> qtArmenterosMin{"qtArmenterosMin", 0.2f, "Minimum Armenteros qT/|alpha| for K0s"};

  Configurable<float> confDaughEta{"confDaughEta", 0.8f, "Maximum V0 daughter eta"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.2f, "Minimum V0 daughter pion pT"};
  Configurable<float> confDaughTPCnclsMin{"confDaughTPCnclsMin", 50.0f, "Minimum V0 daughter TPC clusters"};
  Configurable<float> confDaughTPCncrwsMin{"confDaughTPCncrwsMin", 70.0f, "Minimum V0 daughter TPC crossed rows"};
  Configurable<float> confDaughPIDCuts{"confDaughPIDCuts", 3.0f, "TPC pion PID cut for V0 daughters"};

  Configurable<float> cfgK0sMassMin{"cfgK0sMassMin", 0.45f, "Minimum K0s invariant mass"};
  Configurable<float> cfgK0sMassMax{"cfgK0sMassMax", 0.55f, "Maximum K0s invariant mass"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  RCTFlagsChecker rctChecker;

  void init(o2::framework::InitContext&)
  {
    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);

    histos.add("hCent", "hCent", kTH1F, {{16, 0.0f, 80.0f}});
    histos.add("hEvtSelInfo", "hEvtSelInfo", kTH1F, {{6, 0.0f, 6.0f}});
    histos.add("hTrkSelInfo", "hTrkSelInfo", kTH1F, {{5, 0.0f, 5.0f}});
    histos.add("hV0Info", "hV0Info", kTH1F, {{5, 0.0f, 5.0f}});
    histos.add("hKShortMass", "hKShortMass;M_{#pi^{+}#pi^{-}} (GeV/#it{c}^{2});Counts", kTH1F, {{200, 0.4f, 0.6f}});
    histos.add("hNStoredK0s", "hNStoredK0s;N_{K^{0}_{S}};Events", kTH1F, {{100, 0.0f, 100.0f}});
    histos.add("hNStoredPions", "hNStoredPions;N_{#pi};Events", kTH1F, {{500, 0.0f, 500.0f}});
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    return candidate.isGlobalTrack() &&
           candidate.isPVContributor() &&
           candidate.itsNCls() >= grpPion.itsclusterPiMeson &&
           candidate.tpcNClsCrossedRows() >= grpPion.tpcCrossedRowsPiMeson &&
           std::abs(candidate.dcaXY()) <= grpPion.cutDCAxyPiMeson &&
           std::abs(candidate.dcaZ()) <= grpPion.cutDCAzPiMeson &&
           std::abs(candidate.eta()) <= grpPion.cutEtaPiMeson &&
           candidate.pt() >= grpPion.cutPTPiMeson;
  }

  template <typename T>
  bool selectionPID(const T& candidate)
  {
    const float nTPC = candidate.tpcNSigmaPi();

    if (!candidate.hasTOF()) {
      return std::abs(nTPC) < grpPion.nsigmaCutTPCPiMeson;
    }

    if (candidate.beta() <= grpPion.cutTOFBetaPiMeson) {
      return false;
    }

    return std::abs(nTPC) < grpPion.nsigmaCutTPCPiMeson &&
           std::abs(candidate.tofNSigmaPi()) < grpPion.nsigmaCutTOFPiMeson;
  }

  template <typename T>
  bool selectionPID2(const T& candidate)
  {
    const float nTPC = candidate.tpcNSigmaPi();

    if (!candidate.hasTOF()) {
      return std::abs(nTPC) < grpPion.nsigmaCutTPCPiMeson;
    }

    if (candidate.beta() <= grpPion.cutTOFBetaPiMeson) {
      return false;
    }

    const float nTOF = candidate.tofNSigmaPi();
    const float nCombined = std::sqrt(nTPC * nTPC + nTOF * nTOF);

    return nCombined < grpPion.nsigmaCutTOFPiMeson;
  }

  template <typename T>
  bool selectionPID3(const T& candidate)
  {
    const float pt = candidate.pt();
    const float nTPC = candidate.tpcNSigmaPi();

    if (pt < grpPion.pSwitchPID) {
      return std::abs(nTPC) < grpPion.nsigmaCutTPCPiMeson;
    }

    if (!candidate.hasTOF()) {
      return false;
    }

    if (candidate.beta() <= grpPion.cutTOFBetaPiMeson) {
      return false;
    }

    const float nTOF = candidate.tofNSigmaPi();
    const float nCombined = std::sqrt(nTPC * nTPC + nTOF * nTOF);

    return nCombined < grpPion.nsigmaCutTOFPiMeson;
  }

  template <typename T>
  bool selectionPID4(const T& candidate)
  {
    const float pt = candidate.pt();
    const float nTPC = candidate.tpcNSigmaPi();

    if (pt < grpPion.pSwitchPID) {
      return std::abs(nTPC) < grpPion.nsigmaCutTPCPiMeson;
    }

    if (!candidate.hasTOF()) {
      return false;
    }

    if (candidate.beta() <= grpPion.cutTOFBetaPiMeson) {
      return false;
    }

    const float nTOF = candidate.tofNSigmaPi();

    return std::abs(nTPC) < grpPion.nsigmaCutTPCPiMeson &&
           std::abs(nTOF) < grpPion.nsigmaCutTOFPiMeson;
  }

  template <typename T>
  uint8_t pionPidMask(const T& trk)
  {
    uint8_t mask = 0;

    if (selectionPID(trk)) {
      mask |= kPID1;
    }

    if (selectionPID2(trk)) {
      mask |= kPID2;
    }

    if (selectionPID3(trk)) {
      mask |= kPID3;
    }

    if (selectionPID4(trk)) {
      mask |= kPID4;
    }

    return mask;
  }

  struct StoredK0s {
    float px;
    float py;
    float pz;
    float mass;
    float cospa;
    float radius;
    float dcaPositive;
    float dcaNegative;
    float dcaBetweenDaughters;
    // float lifetime;
    int64_t positiveIndex;
    int64_t negativeIndex;
  };

  struct StoredPion {
    float px;
    float py;
    float pz;
    int8_t charge;
    int64_t index;
    uint8_t pidMask;
  };

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (aod::cent::centFT0C < cfgCutCentralityMax && aod::cent::centFT0C >= cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && aod::track::pt > cfgCutPt);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::EPCalibrationTables, aod::FT0Mults, aod::TPCMults, aod::CentFT0Ms, aod::CentFT0As>>;
  using AllTrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTOFbeta>>;
  using ResoV0s = aod::V0Datas;

  template <typename Collision, typename V0>
  bool selectionV0(Collision const& collision, V0 const& candidate)
  {
    const float pT = candidate.pt();
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = std::abs(candidate.dcaV0daughters());
    const float cpav0 = candidate.v0cosPA();
    const float ctauKShort = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0;

    if (std::abs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }

    if (pT < confV0PtMin || pT > confV0PtMax) {
      return false;
    }

    if (dcaDaughv0 > confV0DCADaughMax) {
      return false;
    }

    if (cpav0 < confV0CPAMin) {
      return false;
    }

    if (tranRad < confV0TranRadV0Min || tranRad > confV0TranRadV0Max) {
      return false;
    }

    if (std::abs(ctauKShort) > cMaxV0LifeTime) {
      return false;
    }

    if ((candidate.qtarm() / std::abs(candidate.alpha())) < qtArmenterosMin) {
      return false;
    }

    if (std::abs(candidate.yK0Short()) > confV0Rap) {
      return false;
    }

    return true;
  }

  template <typename V0>
  bool isSelectedV0Daughter(V0 const& candidate)
  {
    auto postrack = candidate.template posTrack_as<AllTrackCandidates>();
    auto negtrack = candidate.template negTrack_as<AllTrackCandidates>();

    constexpr float minCrossedRowsOverFindable = 0.8f;

    if (postrack.sign() < 0 || negtrack.sign() > 0) {
      return false;
    }

    if (postrack.tpcNClsCrossedRows() < confDaughTPCncrwsMin ||
        negtrack.tpcNClsCrossedRows() < confDaughTPCncrwsMin) {
      return false;
    }

    if (postrack.tpcNClsFound() < confDaughTPCnclsMin ||
        negtrack.tpcNClsFound() < confDaughTPCnclsMin) {
      return false;
    }

    if (postrack.tpcCrossedRowsOverFindableCls() < minCrossedRowsOverFindable ||
        negtrack.tpcCrossedRowsOverFindableCls() < minCrossedRowsOverFindable) {
      return false;
    }

    if (std::abs(postrack.tpcNSigmaPi()) > confDaughPIDCuts ||
        std::abs(negtrack.tpcNSigmaPi()) > confDaughPIDCuts) {
      return false;
    }

    if (candidate.positivept() < cfgDaughPiPt ||
        candidate.negativept() < cfgDaughPiPt) {
      return false;
    }

    if (std::abs(candidate.positiveeta()) > confDaughEta ||
        std::abs(candidate.negativeeta()) > confDaughEta) {
      return false;
    }

    if (std::abs(candidate.dcapostopv()) < cMinV0DCAPi ||
        std::abs(candidate.dcanegtopv()) < cMinV0DCAPi) {
      return false;
    }

    return true;
  }

  template <typename Collision, typename V0>
  bool selectionK0s(Collision const& collision, V0 const& v0)
  {
    if (!isSelectedV0Daughter(v0)) {
      return false;
    }

    if (!selectionV0(collision, v0)) {
      return false;
    }

    if (v0.mK0Short() < cfgK0sMassMin || v0.mK0Short() > cfgK0sMassMax) {
      return false;
    }

    return true;
  }

  void processData(EventCandidates::iterator const& collision,
                   AllTrackCandidates const& tracks,
                   ResoV0s const& v0s)
  {
    o2::aod::ITSResponse itsResponse;

    const float centrality = collision.centFT0C();
    const float vz = collision.posZ();
    const int occupancy = collision.trackOccupancyInTimeRange();

    const float psiFT0C = collision.psiFT0C();
    const float psiFT0A = collision.psiFT0A();
    const float psiTPC = collision.psiTPC();

    histos.fill(HIST("hEvtSelInfo"), 0.5);

    if ((rctCut.requireRCTFlagChecker && !rctChecker(collision)) ||
        !collision.selection_bit(aod::evsel::kNoSameBunchPileup) ||
        !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) ||
        !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) ||
        !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) ||
        (useNoCollInTimeRangeStandard &&
         !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) ||
        !collision.sel8() ||
        (useGoodITSLayersAll &&
         !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) ||
        occupancy >= cfgCutOccupancy) {
      return;
    }

    histos.fill(HIST("hEvtSelInfo"), 1.5);

    if (!collision.triggereventep()) {
      return;
    }

    histos.fill(HIST("hEvtSelInfo"), 2.5);
    histos.fill(HIST("hCent"), centrality);

    std::vector<StoredK0s> selectedK0s;
    std::vector<StoredPion> selectedPions;

    for (const auto& track : tracks) {
      histos.fill(HIST("hTrkSelInfo"), 0.5);

      if (!selectionTrack(track)) {
        continue;
      }

      histos.fill(HIST("hTrkSelInfo"), 1.5);

      const float nSigmaITS =
        itsResponse.nSigmaITS<o2::track::PID::Pion>(track);

      if (grpPion.itsPIDSelection &&
          (nSigmaITS <= grpPion.lowITSPIDNsigma ||
           nSigmaITS >= grpPion.highITSPIDNsigma)) {
        continue;
      }

      histos.fill(HIST("hTrkSelInfo"), 2.5);

      const uint8_t mask = pionPidMask(track);

      if (grpPion.usePID && mask == 0) {
        continue;
      }

      const int8_t charge = static_cast<int8_t>(track.sign());

      if (charge == 0) {
        continue;
      }

      histos.fill(HIST("hTrkSelInfo"), 3.5);

      selectedPions.push_back({track.px(),
                               track.py(),
                               track.pz(),
                               charge,
                               track.globalIndex(),
                               mask});
    }

    for (const auto& v0 : v0s) {
      histos.fill(HIST("hV0Info"), 0.5);

      if (!selectionK0s(collision, v0)) {
        continue;
      }

      auto postrack = v0.template posTrack_as<AllTrackCandidates>();
      auto negtrack = v0.template negTrack_as<AllTrackCandidates>();

      ROOT::Math::PxPyPzMVector pionPos(v0.pxpos(), v0.pypos(), v0.pzpos(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector pionNeg(v0.pxneg(), v0.pyneg(), v0.pzneg(), o2::constants::physics::MassPionCharged);
      ROOT::Math::PxPyPzMVector k0s = pionPos + pionNeg;

      // const float lifetime =
      // v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) *
      // o2::constants::physics::MassK0;

      selectedK0s.push_back({static_cast<float>(k0s.Px()),
                             static_cast<float>(k0s.Py()),
                             static_cast<float>(k0s.Pz()),
                             static_cast<float>(k0s.M()),
                             static_cast<float>(v0.v0cosPA()),
                             static_cast<float>(v0.v0radius()),
                             static_cast<float>(std::abs(v0.dcapostopv())),
                             static_cast<float>(std::abs(v0.dcanegtopv())),
                             static_cast<float>(std::abs(v0.dcaV0daughters())),
                             // lifetime,
                             postrack.globalIndex(),
                             negtrack.globalIndex()});

      histos.fill(HIST("hV0Info"), 1.5);
      histos.fill(HIST("hKShortMass"), k0s.M());
    }

    histos.fill(HIST("hNStoredK0s"), selectedK0s.size());
    histos.fill(HIST("hNStoredPions"), selectedPions.size());

    if (selectedK0s.empty() && selectedPions.empty()) {
      return;
    }

    histos.fill(HIST("hEvtSelInfo"), 3.5);

    kshortpionEvent(centrality,
                    vz,
                    psiFT0C,
                    psiFT0A,
                    psiTPC);

    const int64_t indexEvent = kshortpionEvent.lastIndex();

    for (const auto& k0s : selectedK0s) {
      kshortTrack(indexEvent,
                  k0s.cospa,
                  k0s.radius,
                  k0s.dcaPositive,
                  k0s.dcaNegative,
                  k0s.dcaBetweenDaughters,
                  // k0s.lifetime,
                  k0s.px,
                  k0s.py,
                  k0s.pz,
                  k0s.mass,
                  k0s.positiveIndex,
                  k0s.negativeIndex);
    }

    for (const auto& pion : selectedPions) {
      pionTrack(indexEvent,
                pion.px,
                pion.py,
                pion.pz,
                pion.charge,
                pion.index,
                pion.pidMask);
    }

    histos.fill(HIST("hEvtSelInfo"), 4.5);
  }

  PROCESS_SWITCH(cksspinalignment, processData, "Process data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<cksspinalignment>(cfgc)};
}
