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
/// \brief Table producer for Charged KStar spin alignment
///
/// \author prottay.das@cern.ch

#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "PWGLF/DataModel/LFCKSSpinalignmentTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGMM/Mult/DataModel/Index.h" // for Particles2Tracks table

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector2D.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"

#include <fairlogger/Logger.h>

#include <string>
#include <tuple>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace o2::aod::rctsel;

struct cksspinalignment {

  Produces<aod::KShortpionEvents> kshortpionEvent;
  Produces<aod::KShortTracks> kshortTrack;
  Produces<aod::PionTracks> pionTrack;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;

  Configurable<int> cfgCutOccupancy{"cfgCutOccupancy", 2000, "Occupancy cut"};

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 80.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 0.0f, "Accepted minimum Centrality"};

  // Configs for track
  Configurable<float> cfgCutPt{"cfgCutPt", 0.2, "Pt cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};

  // Configs for pion
  struct : ConfigurableGroup {
    Configurable<bool> itsPIDSelection{"itsPIDSelection", true, "PID ITS"};
    Configurable<float> lowITSPIDNsigma{"lowITSPIDNsigma", -3.0, "lower cut on PID nsigma for ITS"};
    Configurable<float> highITSPIDNsigma{"highITSPIDNsigma", 3.0, "higher cut on PID nsigma for ITS"};
    Configurable<int> itsclusterPiMeson{"itsclusterPiMeson", 5, "Minimum number of ITS cluster for pi meson track"};
    Configurable<int> tpcCrossedRowsPiMeson{"tpcCrossedRowsPiMeson", 80, "Minimum number of TPC Crossed Rows for pi meson track"};
    Configurable<float> cutDCAxyPiMeson{"cutDCAxyPiMeson", 0.1, "Maximum DCAxy for pi meson track"};
    Configurable<float> cutDCAzPiMeson{"cutDCAzPiMeson", 0.1, "Maximum DCAz for pi meson track"};
    Configurable<float> cutEtaPiMeson{"cutEtaPiMeson", 0.8, "Maximum eta for pi meson track"};
    Configurable<float> cutPTPiMeson{"cutPTPiMeson", 0.8, "Maximum pt for pi meson track"};
    Configurable<bool> usePID{"usePID", true, "Flag for using PID selection for pi meson track"};
    Configurable<float> nsigmaCutTPCPiMeson{"nsigmaCutTPCPiMeson", 3.0, "Maximum nsigma cut TPC for pi meson track"};
    Configurable<float> nsigmaCutTOFPiMeson{"nsigmaCutTOFPiMeson", 3.0, "Maximum nsigma cut TOF for pi meson track"};
    Configurable<float> cutTOFBetaPiMeson{"cutTOFBetaPiMeson", 3.0, "Maximum beta cut for pi meson track"};
  } grpPion;

  // Configs for V0
  Configurable<float> confV0PtMin{"confV0PtMin", 0.f, "Minimum transverse momentum of V0"};
  Configurable<float> confV0PtMax{"confV0PtMax", 1000.f, "Maximum transverse momentum of V0"};
  Configurable<float> confV0Rap{"confV0Rap", 0.8f, "Rapidity range of V0"};
  Configurable<float> confV0DCADaughMax{"confV0DCADaughMax", 1.0f, "Maximum DCA between the V0 daughters"};
  Configurable<double> confV0CPAMin{"confV0CPAMin", 0.9998f, "Minimum CPA of V0"};
  Configurable<float> confV0TranRadV0Min{"confV0TranRadV0Min", 1.5f, "Minimum transverse radius"};
  Configurable<float> confV0TranRadV0Max{"confV0TranRadV0Max", 100.f, "Maximum transverse radius"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 1.2, "Maximum V0 DCA to PV"};
  Configurable<float> cMinV0DCAPi{"cMinV0DCAPi", 0.05, "Minimum V0 daughters DCA to PV for Pi"};
  Configurable<float> cMaxV0LifeTime{"cMaxV0LifeTime", 50, "Maximum V0 life time"};
  Configurable<float> qtArmenterosMin{"qtArmenterosMin", 0.2, "Minimum armenteros cut for K0s"};
  // config for V0 daughters
  Configurable<float> confDaughEta{"confDaughEta", 0.8f, "V0 Daugh sel: max eta"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.2, "minimum daughter pion pt"};
  Configurable<float> confDaughTPCnclsMin{"confDaughTPCnclsMin", 50.f, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<float> confDaughTPCncrwsMin{"confDaughTPCncrwsMin", 70.f, "V0 Daugh sel: Min. nCrws TPC"};
  Configurable<float> confDaughPIDCuts{"confDaughPIDCuts", 3, "PID selections for Kshortpion daughters"};

  Configurable<int> iMNbins{"iMNbins", 50, "Number of bins in invariant mass"};
  Configurable<float> lbinIM{"lbinIM", 1.09, "lower bin value in IM histograms"};
  Configurable<float> hbinIM{"hbinIM", 1.14, "higher bin value in IM histograms"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  RCTFlagsChecker rctChecker;
  void init(o2::framework::InitContext&)
  {
    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerZDCCheck, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);
    AxisSpec thnAxisInvMass{iMNbins, lbinIM, hbinIM, "#it{M} (GeV/#it{c}^{2})"};
    histos.add("hEvtSelInfo", "hEvtSelInfo", kTH1F, {{5, 0, 5.0}});
    histos.add("hTrkSelInfo", "hTrkSelInfo", kTH1F, {{5, 0, 5.0}});
    histos.add("hKShortMass", "hKShortMass", kTH1F, {thnAxisInvMass});
    histos.add("hV0Info", "hV0Info", kTH1F, {{5, 0, 5.0}});
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() >= grpPion.itsclusterPiMeson && candidate.tpcNClsCrossedRows() > grpPion.tpcCrossedRowsPiMeson && std::abs(candidate.dcaXY()) <= grpPion.cutDCAxyPiMeson && std::abs(candidate.dcaZ()) <= grpPion.cutDCAzPiMeson && std::abs(candidate.eta()) <= grpPion.cutEtaPiMeson && candidate.pt() >= grpPion.cutPTPiMeson) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < grpPion.nsigmaCutTPCPiMeson) {
      return true;
    }
    if (candidate.hasTOF() && candidate.beta() > grpPion.cutTOFBetaPiMeson && std::abs(candidate.tpcNSigmaPi()) < grpPion.nsigmaCutTPCPiMeson && std::abs(candidate.tofNSigmaPi()) < grpPion.nsigmaCutTOFPiMeson) {
      return true;
    }
    return false;
  }

  template <typename Collision, typename V0>
  bool selectionV0(Collision const& collision, V0 const& candidate)
  {
    if (std::abs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }
    const float pT = candidate.pt();
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = std::abs(candidate.dcaV0daughters());
    const float cpav0 = candidate.v0cosPA();
    float ctauKShort = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * (o2::constants::physics::MassK0);

    if (pT < confV0PtMin) {
      return false;
    }
    if (pT > confV0PtMax) {
      return false;
    }
    if (dcaDaughv0 > confV0DCADaughMax) {
      return false;
    }
    if (cpav0 < confV0CPAMin) {
      return false;
    }
    if (tranRad < confV0TranRadV0Min) {
      return false;
    }
    if (tranRad > confV0TranRadV0Max) {
      return false;
    }
    if (std::abs(ctauKShort) > cMaxV0LifeTime) {
      return false;
    }
    if ((candidate.qtarm() / std::abs(candidate.alpha())) < qtArmenterosMin) {
      return false;
    }
    if (std::abs(candidate.yK0Short()) > confV0Rap) { // use full rapidity 0.8 for K0s
      return false;
    }
    return true;
  }

  template <typename V0>
  bool isSelectedV0Daughter(V0 const& candidate)
  {
    auto postrack = candidate.template posTrack_as<AllTrackCandidates>();
    auto negtrack = candidate.template negTrack_as<AllTrackCandidates>();

    const auto ncrfc = 0.8;

    if (postrack.sign() < 0 || negtrack.sign() > 0) {
      return false;
    }
    if (postrack.tpcNClsCrossedRows() < confDaughTPCncrwsMin || negtrack.tpcNClsCrossedRows() < confDaughTPCncrwsMin) {
      return false;
    }
    if (postrack.tpcNClsFound() < confDaughTPCnclsMin || negtrack.tpcNClsFound() < confDaughTPCnclsMin) {
      return false;
    }
    if (postrack.tpcCrossedRowsOverFindableCls() < ncrfc || negtrack.tpcCrossedRowsOverFindableCls() < ncrfc) {
      return false;
    }
    if (std::abs(postrack.tpcNSigmaPi()) > confDaughPIDCuts || std::abs(negtrack.tpcNSigmaPi()) > confDaughPIDCuts) {
      return false;
    }
    if (candidate.positivept() < cfgDaughPiPt || candidate.negativept() < cfgDaughPiPt) {
      return false;
    }
    if (std::abs(candidate.positiveeta()) > confDaughEta || std::abs(candidate.negativeeta()) > confDaughEta) {
      return false;
    }
    if (std::abs(candidate.dcapostopv()) < cMinV0DCAPi || std::abs(candidate.dcanegtopv()) < cMinV0DCAPi) {
      return false;
    }

    return true;
  }

  std::tuple<int, bool> getK0sTags(const auto& v0, const auto& collision)
  {
    // auto postrack = v0.template posTrack_as<AllTrackCandidates>();
    // auto negtrack = v0.template negTrack_as<AllTrackCandidates>();

    int kshortTag = 0;

    if (isSelectedV0Daughter(v0) && v0.mK0Short() > lbinIM && v0.mK0Short() < hbinIM) {
      kshortTag = 1;
    }

    if (!kshortTag) {
      return {0, false}; // No valid tags
    }

    if (!selectionV0(collision, v0)) {
      return {0, false}; // Fails selection
    }

    return {kshortTag, true}; // Valid candidate
  }

  ROOT::Math::PxPyPzMVector kshort, pion, pionbach, antiPion;
  ROOT::Math::PxPyPzMVector kshortDummy, pionDummy;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && (aod::track::pt) > cfgCutPt);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::EPCalibrationTables, aod::FT0Mults, aod::TPCMults, aod::CentFT0Ms, aod::CentFT0As>>;
  using AllTrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTOFbeta>>;
  using ResoV0s = aod::V0Datas;
  void processData(EventCandidates::iterator const& collision, AllTrackCandidates const& tracks, ResoV0s const& V0s)
  {
    o2::aod::ITSResponse itsResponse;
    std::vector<ROOT::Math::PxPyPzMVector> kshortMother, pionBachelor;
    std::vector<float> v0Cospa = {};
    std::vector<float> v0Radius = {};
    std::vector<float> dcaPositive = {};
    std::vector<float> dcaNegative = {};
    std::vector<int> positiveIndex = {};
    std::vector<int> negativeIndex = {};
    std::vector<float> dcaBetweenDaughter = {};
    std::vector<float> v0Lifetime = {};
    // std::vector<float> armenteros = {};
    std::vector<float> pionBachelorIndex = {};
    // std::vector<float> pionBachelorSign = {};
    std::vector<float> pionBachelorTPC = {};
    std::vector<float> pionBachelorTOF = {};
    std::vector<float> pionBachelorTOFHit = {};

    int numbV0 = 0;
    auto centrality = collision.centFT0C();
    auto vz = collision.posZ();
    int occupancy = collision.trackOccupancyInTimeRange();
    auto psiFT0C = collision.psiFT0C();
    auto psiFT0A = collision.psiFT0A();
    auto psiTPC = collision.psiTPC();
    // auto psiTPCR = collision.psiTPCR();
    // auto psiTPCL = collision.psiTPCL();
    histos.fill(HIST("hEvtSelInfo"), 0.5);
    if ((rctCut.requireRCTFlagChecker && rctChecker(collision)) && collision.selection_bit(aod::evsel::kNoSameBunchPileup) && collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoITSROFrameBorder) && collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard) && collision.sel8() && collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll) && occupancy < cfgCutOccupancy) {
      histos.fill(HIST("hEvtSelInfo"), 1.5);
      if (collision.triggereventep()) {
        histos.fill(HIST("hEvtSelInfo"), 2.5);

        for (const auto& track1 : tracks) {
          histos.fill(HIST("hTrkSelInfo"), 0.5);
          if (!selectionTrack(track1)) {
            continue;
          }
          histos.fill(HIST("hTrkSelInfo"), 1.5);

          if (grpPion.itsPIDSelection && !(itsResponse.nSigmaITS<o2::track::PID::Pion>(track1) > grpPion.lowITSPIDNsigma && itsResponse.nSigmaITS<o2::track::PID::Pion>(track1) < grpPion.highITSPIDNsigma)) {
            continue;
          }
          histos.fill(HIST("hTrkSelInfo"), 2.5);

          if (grpPion.usePID && !selectionPID(track1)) {
            continue;
          }
          histos.fill(HIST("hTrkSelInfo"), 3.5);

          auto track1ID = track1.globalIndex();
          auto track1sign = track1.sign();
          if (track1sign == 0)
            continue;
          auto track1nsigTPC = track1.tpcNSigmaPi();
          auto track1nsigTOF = -999.9;
          auto track1TOFHit = -1;
          if (track1.hasTOF()) {
            track1TOFHit = 1;
            track1nsigTOF = track1.tofNSigmaPi();
            histos.fill(HIST("hTrkSelInfo"), 4.5);
          }
          pionbach = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), o2::constants::physics::MassPionCharged);
          pionBachelor.push_back(pionbach);
          pionBachelorIndex.push_back(track1ID);
          // pionBachelorSign.push_back(track1sign);
          pionBachelorTPC.push_back(track1nsigTPC);
          pionBachelorTOF.push_back(track1nsigTOF);
          pionBachelorTOFHit.push_back(track1TOFHit);
        }
        for (const auto& v0 : V0s) {
          histos.fill(HIST("hV0Info"), 0.5);
          auto [kshortTag, isValid] = getK0sTags(v0, collision);
          if (kshortTag && isValid) {
            histos.fill(HIST("hV0Info"), 1.5);
            auto postrack1 = v0.template posTrack_as<AllTrackCandidates>();
            auto negtrack1 = v0.template negTrack_as<AllTrackCandidates>();
            positiveIndex.push_back(postrack1.globalIndex());
            negativeIndex.push_back(negtrack1.globalIndex());
            v0Cospa.push_back(v0.v0cosPA());
            v0Radius.push_back(v0.v0radius());
            dcaPositive.push_back(std::abs(v0.dcapostopv()));
            dcaNegative.push_back(std::abs(v0.dcanegtopv()));
            dcaBetweenDaughter.push_back(std::abs(v0.dcaV0daughters()));
            v0Lifetime.push_back(v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * (o2::constants::physics::MassK0));
            // armenteros.push_back((v0.qtarm() / std::abs(v0.alpha())));

            pion = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), o2::constants::physics::MassPionCharged);
            antiPion = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), o2::constants::physics::MassPionCharged);
            kshort = pion + antiPion;
            // chargedkstar = kshort + pionbach;
            kshortMother.push_back(kshort);
            // chargedkstarMother.push_back(chargedkstar);
            histos.fill(HIST("hKShortMass"), kshort.M());
          }
          numbV0 = numbV0 + 1;
        }
        if (numbV0 > 1 && v0Cospa.size() > 1) {
          histos.fill(HIST("hEvtSelInfo"), 3.5);
          kshortpionEvent(centrality, vz, collision.index(), psiFT0C, psiFT0A, psiTPC);
          auto indexEvent = kshortpionEvent.lastIndex();
          //// Fill track table for Charged KStar//////////////////
          for (auto icks = kshortMother.begin(); icks != kshortMother.end(); ++icks) {
            auto iter = std::distance(kshortMother.begin(), icks);
            kshortDummy = kshortMother.at(iter);
            // chargedkstarDummy = chargedkstarMother.at(iter);

            kshortTrack(indexEvent, v0Cospa.at(iter), v0Radius.at(iter), dcaPositive.at(iter), dcaNegative.at(iter), dcaBetweenDaughter.at(iter), v0Lifetime.at(iter), kshortDummy.Px(), kshortDummy.Py(), kshortDummy.Pz(), kshortDummy.M(), positiveIndex.at(iter), negativeIndex.at(iter));
          }
          for (auto ipi = pionBachelor.begin(); ipi != pionBachelor.end(); ++ipi) {
            auto iterpi = std::distance(pionBachelor.begin(), ipi);
            pionDummy = pionBachelor.at(iterpi);

            pionTrack(indexEvent, pionDummy.Px(), pionDummy.Py(), pionDummy.Pz(), pionBachelorTPC.at(iterpi), pionBachelorTOFHit.at(iterpi), pionBachelorTOF.at(iterpi), pionBachelorIndex.at(iterpi));
          }
        }
      }
    }
  }
  PROCESS_SWITCH(cksspinalignment, processData, "Process data", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<cksspinalignment>(cfgc)};
}
