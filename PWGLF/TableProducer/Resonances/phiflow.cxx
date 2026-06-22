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

/// \file phiflow.cxx
/// \brief Table producer for Phi FLow
///
/// \author prottay.das@cern.ch

#include "PWGLF/DataModel/LFPhiFlowTables.h"
#include "PWGLF/DataModel/SPCalibrationTables.h"

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

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::rctsel;

struct phiflow {

  Produces<aod::KaonkaonEvents> kaonkaonEvent;
  Produces<aod::KaonTracks> kaonTrack;

  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;

  Configurable<bool> useNoCollInTimeRangeStandard{"useNoCollInTimeRangeStandard", false, "Apply kNoCollInTimeRangeStandard selection bit"};
  Configurable<bool> useGoodITSLayersAll{"useGoodITSLayersAll", true, "Apply kIsGoodITSLayersAll selection bit"};
  Configurable<int> cfgCutOccupancy{"cfgCutOccupancy", 2000, "Occupancy cut"};

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 80.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 0.0f, "Accepted minimum Centrality"};

  // Configs for track filtering
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "Charge cut on daughter"};
  Configurable<float> cfgCutPt{"cfgCutPt", 0.2, "Pt cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 0.2, "DCAxy cut on daughter track"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 0.2, "DCAz cut on daughter track"};
  Configurable<float> nsigmaCutTPC{"nsigmaCutTPC", 3.0, "Maximum nsigma cut TPC for filtered kaon track"};

  // Configs for kaon
  struct : ConfigurableGroup {
    Configurable<bool> itsPIDSelection{"itsPIDSelection", true, "PID ITS"};
    Configurable<float> lowITSPIDNsigma{"lowITSPIDNsigma", -3.0, "lower cut on PID nsigma for ITS"};
    Configurable<float> highITSPIDNsigma{"highITSPIDNsigma", 3.0, "higher cut on PID nsigma for ITS"};
    Configurable<int> itsclusterKaMeson{"itsclusterKaMeson", 5, "Minimum number of ITS cluster for kaon meson track"};
    Configurable<int> tpcCrossedRowsKaMeson{"tpcCrossedRowsKaMeson", 80, "Minimum number of TPC Crossed Rows for kaon meson track"};
    Configurable<float> cutDCAxyKaMeson{"cutDCAxyKaMeson", 0.1, "Maximum DCAxy for kaon meson track"};
    Configurable<float> cutDCAzKaMeson{"cutDCAzKaMeson", 0.1, "Maximum DCAz for kaon meson track"};
    Configurable<float> cutEtaKaMeson{"cutEtaKaMeson", 0.8, "Maximum eta for kaon meson track"};
    Configurable<float> cutPTKaMeson{"cutPTKaMeson", 0.8, "Minimum pt for kaon meson track"};
    Configurable<bool> usePID{"usePID", true, "Flag for using PID selection for kaon meson track"};
    Configurable<float> nsigmaCutTPCKaMeson{"nsigmaCutTPCKaMeson", 3.0, "Maximum nsigma cut TPC for kaon meson track"};
    Configurable<float> nsigmaCutTOFKaMeson{"nsigmaCutTOFKaMeson", 3.0, "Maximum nsigma cut TOF for kaon meson track"};
    Configurable<float> cutTOFBetaKaMeson{"cutTOFBetaKaMeson", 1.0, "Maximum beta cut for kaon meson track"};
    Configurable<bool> isDeepAngle{"isDeepAngle", true, "Deep Angle cut"};
    Configurable<double> cutDeepAngle{"cutDeepAngle", 0.04, "Deep Angle cut value"};
  } grpKaon;

  enum KaonPidBits : uint8_t {
    kPID1 = 1u << 0, // selectionPID
    kPID2 = 1u << 1, // selectionPID2
    kPID3 = 1u << 2, // selectionPID3
    kPID4 = 1u << 3  // selectionPID4
  };

  // configurable for chargedkstar
  Configurable<float> cfgPhiMassMin{"cfgPhiMassMin", 0.9f, "Phi mass min"};
  Configurable<float> cfgPhiMassMax{"cfgPhiMassMax", 1.2f, "Phi mass max"};

  Configurable<int> iMNbins{"iMNbins", 120, "Number of bins in invariant mass"};
  Configurable<float> lbinIM{"lbinIM", 0.9, "lower bin value in IM histograms"};
  Configurable<float> hbinIM{"hbinIM", 1.2, "higher bin value in IM histograms"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  RCTFlagsChecker rctChecker;
  void init(o2::framework::InitContext&)
  {
    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);
    AxisSpec thnAxisInvMass{iMNbins, lbinIM, hbinIM, "#it{M} (GeV/#it{c}^{2})"};

    histos.add("hCent", "hCent", kTH1F, {{8, 0, 80.0}});
    histos.add("hEvtSelInfo", "hEvtSelInfo", kTH1F, {{5, 0, 5.0}});
    histos.add("hTrkSelInfo", "hTrkSelInfo", kTH1F, {{10, 0, 10.0}});
    histos.add("hPhiMass", "hPhiMass", kTH1F, {thnAxisInvMass});
    histos.add("hInfo", "hInfo", kTH1F, {{5, 0, 5.0}});
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    return candidate.isGlobalTrack() &&
           candidate.isPVContributor() &&
           candidate.itsNCls() >= grpKaon.itsclusterKaMeson &&
           candidate.tpcNClsCrossedRows() > grpKaon.tpcCrossedRowsKaMeson &&
           std::abs(candidate.dcaXY()) <= grpKaon.cutDCAxyKaMeson &&
           std::abs(candidate.dcaZ()) <= grpKaon.cutDCAzKaMeson &&
           std::abs(candidate.eta()) <= grpKaon.cutEtaKaMeson &&
           candidate.pt() >= grpKaon.cutPTKaMeson;
  }

  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < grpKaon.nsigmaCutTPCKaMeson) {
      return true;
    }
    if (candidate.hasTOF() && candidate.beta() > grpKaon.cutTOFBetaKaMeson && std::abs(candidate.tpcNSigmaKa()) < grpKaon.nsigmaCutTPCKaMeson && std::abs(candidate.tofNSigmaKa()) < grpKaon.nsigmaCutTOFKaMeson) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectionPID2(const T& candidate)
  {
    if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < grpKaon.nsigmaCutTPCKaMeson) {
      return true;
    }
    if (candidate.hasTOF() && candidate.beta() > grpKaon.cutTOFBetaKaMeson && std::sqrt(candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) < grpKaon.nsigmaCutTOFKaMeson) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectionPID3(const T& candidate)
  {
    constexpr float pSwitch = 0.5f;

    const float pt = candidate.pt();
    const float nTPC = candidate.tpcNSigmaKa();

    if (pt < pSwitch && !candidate.hasTOF()) {
      return std::abs(nTPC) < grpKaon.nsigmaCutTPCKaMeson;
    }

    if (!candidate.hasTOF()) {
      return false;
    }

    const float nTOF = candidate.tofNSigmaKa();
    const float nCombined = std::sqrt(nTPC * nTPC + nTOF * nTOF);

    return nCombined < grpKaon.nsigmaCutTOFKaMeson;
  }

  template <typename T>
  bool selectionPID4(const T& candidate)
  {
    const float px = candidate.px();
    const float py = candidate.py();
    // const float pz = candidate.pz();
    const float pt = std::sqrt(px * px + py * py);

    constexpr float pSwitch = 0.5f; // GeV/c

    const float nTPC = candidate.tpcNSigmaKa();

    // Low momentum: TPC-only, TOF not required and not used
    if (pt < pSwitch) {
      return std::abs(nTPC) < grpKaon.nsigmaCutTPCKaMeson; // e.g. 3
    }

    // High momentum: TOF hit mandatory + separate 3σ cuts
    if (!candidate.hasTOF()) {
      return false;
    }

    const float nTOF = candidate.tofNSigmaKa();
    return (std::abs(nTPC) < grpKaon.nsigmaCutTPCKaMeson) &&
           (std::abs(nTOF) < grpKaon.nsigmaCutTOFKaMeson);
  }

  template <typename T>
  uint8_t kaonPidMask(const T& trk)
  {
    uint8_t m = 0;

    if (selectionPID(trk)) {
      m |= kPID1;
    }
    if (selectionPID2(trk)) {
      m |= kPID2;
    }
    if (selectionPID3(trk)) {
      m |= kPID3;
    }
    if (selectionPID4(trk)) {
      m |= kPID4;
    }
    return m;
  }

  // deep angle cut on pair to remove photon conversion
  template <typename T1, typename T2>
  bool selectionPair(const T1& candidate1, const T2& candidate2)
  {
    const double pt1 = candidate1.pt();
    const double pt2 = candidate2.pt();
    const double pz1 = candidate1.pz();
    const double pz2 = candidate2.pz();
    const double p1 = candidate1.p();
    const double p2 = candidate2.p();

    if (p1 <= 0.0 || p2 <= 0.0) {
      return false;
    }

    const double cosAngle = std::clamp((pt1 * pt2 + pz1 * pz2) / (p1 * p2), -1.0, 1.0);
    const double angle = std::acos(cosAngle);

    return !(grpKaon.isDeepAngle && angle < grpKaon.cutDeepAngle);
  }

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (aod::cent::centFT0C < cfgCutCentralityMax && aod::cent::centFT0C > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && (aod::track::pt) > cfgCutPt);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  Filter PIDcutFilter = nabs(aod::pidtpc::tpcNSigmaKa) < nsigmaCutTPC;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::FT0Mults, aod::TPCMults, aod::SPCalibrationTables>>;
  using AllTrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTOFbeta>>;
  SliceCache cache;
  Partition<AllTrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<AllTrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  static constexpr double massKa = o2::constants::physics::MassKPlus;

  void processData(EventCandidates::iterator const& collision,
                   AllTrackCandidates const& /*tracks*/)
  {
    o2::aod::ITSResponse itsResponse;

    const auto centrality = collision.centFT0C();
    const auto vz = collision.posZ();
    const int occupancy = collision.trackOccupancyInTimeRange();

    const auto qxZDCA = collision.qxZDCA();
    const auto qxZDCC = collision.qxZDCC();
    const auto qyZDCA = collision.qyZDCA();
    const auto qyZDCC = collision.qyZDCC();

    histos.fill(HIST("hEvtSelInfo"), 0.5);

    if (!((!rctCut.requireRCTFlagChecker || rctChecker(collision)) &&
          collision.selection_bit(aod::evsel::kNoSameBunchPileup) &&
          collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) &&
          (!useNoCollInTimeRangeStandard || collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) &&
          collision.sel8() &&
          (!useGoodITSLayersAll || collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) &&
          occupancy < cfgCutOccupancy)) {
      return;
    }

    histos.fill(HIST("hEvtSelInfo"), 1.5);
    histos.fill(HIST("hCent"), centrality);

    std::vector<ROOT::Math::PtEtaPhiMVector> phiresonance;
    std::vector<ROOT::Math::PtEtaPhiMVector> phiresonanced1;
    std::vector<ROOT::Math::PtEtaPhiMVector> phiresonanced2;
    std::vector<int64_t> Phid1Index;
    std::vector<int64_t> Phid2Index;
    std::vector<uint16_t> PhiPairPidMask;

    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId,
                                                collision.globalIndex(),
                                                cache);

    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId,
                                                collision.globalIndex(),
                                                cache);

    for (const auto& track1 : posThisColl) {
      histos.fill(HIST("hTrkSelInfo"), 0.5);

      if (!selectionTrack(track1)) {
        continue;
      }

      histos.fill(HIST("hTrkSelInfo"), 1.5);

      const float nSigmaITS1 = itsResponse.nSigmaITS<o2::track::PID::Kaon>(track1);
      if (grpKaon.itsPIDSelection &&
          (nSigmaITS1 <= grpKaon.lowITSPIDNsigma ||
           nSigmaITS1 >= grpKaon.highITSPIDNsigma)) {
        continue;
      }

      histos.fill(HIST("hTrkSelInfo"), 2.5);

      const uint8_t mask = kaonPidMask(track1);
      if (grpKaon.usePID && mask == 0) {
        continue;
      }

      histos.fill(HIST("hTrkSelInfo"), 3.5);

      const auto track1ID = track1.globalIndex();

      for (const auto& track2 : negThisColl) {
        histos.fill(HIST("hTrkSelInfo"), 4.5);

        if (!selectionTrack(track2)) {
          continue;
        }

        histos.fill(HIST("hTrkSelInfo"), 5.5);

        const float nSigmaITS2 = itsResponse.nSigmaITS<o2::track::PID::Kaon>(track2);
        if (grpKaon.itsPIDSelection &&
            (nSigmaITS2 <= grpKaon.lowITSPIDNsigma ||
             nSigmaITS2 >= grpKaon.highITSPIDNsigma)) {
          continue;
        }

        histos.fill(HIST("hTrkSelInfo"), 6.5);

        const uint8_t mask2 = kaonPidMask(track2);
        if (grpKaon.usePID && mask2 == 0) {
          continue;
        }

        histos.fill(HIST("hTrkSelInfo"), 7.5);

        const auto track2ID = track2.globalIndex();

        if (track2ID == track1ID) {
          continue;
        }

        if (!selectionPair(track1, track2)) {
          continue;
        }

        ROOT::Math::PxPyPzMVector KaonPlus(track1.px(), track1.py(), track1.pz(), massKa);
        ROOT::Math::PxPyPzMVector KaonMinus(track2.px(), track2.py(), track2.pz(), massKa);
        auto PhiMesonMother = KaonPlus + KaonMinus;

        histos.fill(HIST("hPhiMass"), PhiMesonMother.M());

        ROOT::Math::PtEtaPhiMVector temp1(track1.pt(), track1.eta(), track1.phi(), massKa);
        ROOT::Math::PtEtaPhiMVector temp2(track2.pt(), track2.eta(), track2.phi(), massKa);
        ROOT::Math::PtEtaPhiMVector temp3(PhiMesonMother.pt(), PhiMesonMother.eta(), PhiMesonMother.phi(), PhiMesonMother.M());

        phiresonanced1.push_back(temp1);
        phiresonanced2.push_back(temp2);
        phiresonance.push_back(temp3);

        Phid1Index.push_back(track1ID);
        Phid2Index.push_back(track2ID);

        const uint16_t pairPidMask = (static_cast<uint16_t>(mask2) << 8) | mask;
        PhiPairPidMask.push_back(pairPidMask);
      }
    }

    // Optional: if no accepted K+K- pair, do not write this event
    if (phiresonance.empty()) {
      return;
    }

    kaonkaonEvent(centrality, vz, qxZDCA, qxZDCC, qyZDCA, qyZDCC);
    const auto indexEvent = kaonkaonEvent.lastIndex();

    for (size_t i = 0; i < phiresonance.size(); ++i) {
      const auto& phiCand = phiresonance[i];
      const auto& d1 = phiresonanced1[i];
      const auto& d2 = phiresonanced2[i];

      kaonTrack(indexEvent,
                d1.Px(), d1.Py(), d1.Pz(),
                d2.Px(), d2.Py(), d2.Pz(),
                phiCand.M(),
                Phid1Index[i],
                Phid2Index[i],
                PhiPairPidMask[i]);
    }
  }

  PROCESS_SWITCH(phiflow, processData, "Process data", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<phiflow>(cfgc)};
}
