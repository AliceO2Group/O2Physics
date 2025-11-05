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

/// \file taskLambdaSpinCorr.cxx
/// \brief Analysis task for Lambda spin spin correlation
///
/// \author prottay.das@cern.ch
/// \author sourav.kundu@cern.ch

#include "PWGLF/DataModel/LFSpincorrelationTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGMM/Mult/DataModel/Index.h" // for Particles2Tracks table

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
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

struct lambdaspincorrelation {

  Produces<aod::LambdaEvents> lambdaEvent;
  Produces<aod::LambdaPairs> lambdaPair;
  Produces<aod::LambdaEventmcs> lambdaEventmc;
  Produces<aod::LambdaPairmcs> lambdaPairmc;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerZDCCheck{"cfgEvtRCTFlagCheckerZDCCheck", false, "Evt sel: RCT flag checker ZDC check"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;
  Configurable<bool> useNoCollInTimeRangeStandard{"useNoCollInTimeRangeStandard", false, "Apply kNoCollInTimeRangeStandard selection bit"};
  Configurable<bool> useGoodITSLayersAll{"useGoodITSLayersAll", true, "Apply kIsGoodITSLayersAll selection bit"};
  // mixing
  // Produce derived tables
  Configurable<int> cfgCutOccupancy{"cfgCutOccupancy", 2000, "Occupancy cut"};
  ConfigurableAxis axisVertex{"axisVertex", {5, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {8, 0, 80}, "multiplicity percentile for bin"};

  // events
  Configurable<float> cfgEventTypepp{"cfgEventTypepp", false, "Type of collisions"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 80.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 0.0f, "Accepted minimum Centrality"};

  // Configs for track
  Configurable<float> cfgCutPt{"cfgCutPt", 0.2, "Pt cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};

  // Configs for V0
  Configurable<float> confV0PtMin{"confV0PtMin", 0.f, "Minimum transverse momentum of V0"};
  Configurable<float> confV0PtMax{"confV0PtMax", 1000.f, "Maximum transverse momentum of V0"};
  Configurable<float> confV0Rap{"confV0Rap", 0.8f, "Rapidity range of V0"};
  Configurable<float> confV0DCADaughMax{"confV0DCADaughMax", 1.0f, "Maximum DCA between the V0 daughters"};
  Configurable<double> confV0CPAMin{"confV0CPAMin", 0.9998f, "Minimum CPA of V0"};
  Configurable<float> confV0TranRadV0Min{"confV0TranRadV0Min", 1.5f, "Minimum transverse radius"};
  Configurable<float> confV0TranRadV0Max{"confV0TranRadV0Max", 100.f, "Maximum transverse radius"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 1.2, "Maximum V0 DCA to PV"};
  Configurable<float> cMinV0DCAPr{"cMinV0DCAPr", 0.05, "Minimum V0 daughters DCA to PV for Pr"};
  Configurable<float> cMinV0DCAPi{"cMinV0DCAPi", 0.05, "Minimum V0 daughters DCA to PV for Pi"};
  Configurable<float> cMaxV0LifeTime{"cMaxV0LifeTime", 50, "Maximum V0 life time"};

  // config for V0 daughters
  Configurable<float> confDaughEta{"confDaughEta", 0.8f, "V0 Daugh sel: max eta"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.2, "minimum daughter proton pt"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.2, "minimum daughter pion pt"};
  Configurable<float> confDaughTPCnclsMin{"confDaughTPCnclsMin", 50.f, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<float> confDaughPIDCuts{"confDaughPIDCuts", 3, "PID selections for Lambda daughters"};

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
    histos.add("hLambdaMass", "hLambdaMass", kTH1F, {thnAxisInvMass});
    histos.add("hV0Info", "hV0Info", kTH1F, {{5, 0, 5.0}});
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
    float ctauLambda = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * (o2::constants::physics::MassLambda);

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
    if (std::abs(ctauLambda) > cMaxV0LifeTime) {
      return false;
    }
    if (std::abs(candidate.yLambda()) > confV0Rap) {
      return false;
    }
    return true;
  }

  template <typename V0, typename T>
  bool isSelectedV0Daughter(V0 const& candidate, T const& track, int pid)
  {
    const auto tpcNClsF = track.tpcNClsFound();
    const auto ncr = 70;
    const auto ncrfc = 0.8;

    if (track.tpcNClsCrossedRows() < ncr) {
      return false;
    }
    if (tpcNClsF < confDaughTPCnclsMin) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < ncrfc) {
      return false;
    }

    if (pid == 0 && std::abs(track.tpcNSigmaPr()) > confDaughPIDCuts) {
      return false;
    }
    if (pid == 1 && std::abs(track.tpcNSigmaPi()) > confDaughPIDCuts) {
      return false;
    }
    if (pid == 0 && (candidate.positivept() < cfgDaughPrPt || candidate.negativept() < cfgDaughPiPt)) {
      return false;
    }
    if (pid == 1 && (candidate.positivept() < cfgDaughPiPt || candidate.negativept() < cfgDaughPrPt)) {
      return false;
    }
    if (std::abs(candidate.positiveeta()) > confDaughEta || std::abs(candidate.negativeeta()) > confDaughEta) {
      return false;
    }

    if (pid == 0 && (std::abs(candidate.dcapostopv()) < cMinV0DCAPr || std::abs(candidate.dcanegtopv()) < cMinV0DCAPi)) {
      return false;
    }
    if (pid == 1 && (std::abs(candidate.dcapostopv()) < cMinV0DCAPi || std::abs(candidate.dcanegtopv()) < cMinV0DCAPr)) {
      return false;
    }

    return true;
  }

  std::tuple<int, int, bool> getLambdaTags(const auto& v0, const auto& collision)
  {
    auto postrack = v0.template posTrack_as<AllTrackCandidates>();
    auto negtrack = v0.template negTrack_as<AllTrackCandidates>();

    int lambdaTag = 0;
    int aLambdaTag = 0;

    const auto signpos = postrack.sign();
    const auto signneg = negtrack.sign();

    if (signpos < 0 || signneg > 0) {
      return {0, 0, false}; // Invalid candidate
    }

    if (isSelectedV0Daughter(v0, postrack, 0) && isSelectedV0Daughter(v0, negtrack, 1) && v0.mLambda() > lbinIM && v0.mLambda() < hbinIM) {
      lambdaTag = 1;
    }
    if (isSelectedV0Daughter(v0, negtrack, 0) && isSelectedV0Daughter(v0, postrack, 1) && v0.mAntiLambda() > lbinIM && v0.mAntiLambda() < hbinIM) {
      aLambdaTag = 1;
    }

    if (!lambdaTag && !aLambdaTag) {
      return {0, 0, false}; // No valid tags
    }

    if (!selectionV0(collision, v0)) {
      return {0, 0, false}; // Fails selection
    }

    return {lambdaTag, aLambdaTag, true}; // Valid candidate
  }

  ROOT::Math::PxPyPzMVector lambda, antiLambda, proton, pion, antiProton, antiPion;
  ROOT::Math::PxPyPzMVector lambdaDummy, pionDummy, protonDummy;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs, aod::CentFT0Ms>>;
  using AllTrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr>;
  using ResoV0s = aod::V0Datas;

  void processData(EventCandidates::iterator const& collision, AllTrackCandidates const&, ResoV0s const& V0s)
  {
    std::vector<ROOT::Math::PxPyPzMVector> lambdaMother, protonDaughter, pionDaughter;
    std::vector<int> v0Status = {};
    std::vector<bool> doubleStatus = {};
    std::vector<float> v0Cospa = {};
    std::vector<float> v0Radius = {};
    std::vector<float> dcaPositive = {};
    std::vector<float> dcaNegative = {};
    std::vector<int> positiveIndex = {};
    std::vector<int> negativeIndex = {};
    std::vector<float> dcaBetweenDaughter = {};
    int numbV0 = 0;
    // LOGF(info, "event collisions: (%d)", collision.index());
    auto centrality = collision.centFT0C();
    if (cfgEventTypepp)
      centrality = collision.centFT0M();
    auto vz = collision.posZ();
    int occupancy = collision.trackOccupancyInTimeRange();
    histos.fill(HIST("hEvtSelInfo"), 0.5);
    // if ((!rctCut.requireRCTFlagChecker || rctChecker(collision)) && collision.selection_bit(aod::evsel::kNoSameBunchPileup) && collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoITSROFrameBorder) && collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard) && collision.sel8() && collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll) && occupancy < cfgCutOccupancy) {
    if ((!rctCut.requireRCTFlagChecker || rctChecker(collision)) && collision.selection_bit(aod::evsel::kNoSameBunchPileup) && collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) && collision.selection_bit(aod::evsel::kNoTimeFrameBorder) && collision.selection_bit(aod::evsel::kNoITSROFrameBorder) && (!useNoCollInTimeRangeStandard || collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) && collision.sel8() && (!useGoodITSLayersAll || collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) && occupancy < cfgCutOccupancy) {
      histos.fill(HIST("hEvtSelInfo"), 1.5);
      for (const auto& v0 : V0s) {
        // LOGF(info, "v0 index 0 : (%d)", v0.index());
        auto [lambdaTag, aLambdaTag, isValid] = getLambdaTags(v0, collision);
        if (isValid) {
          // LOGF(info, "v0 index 1 : (%d)", v0.index());
          if (lambdaTag) {
            histos.fill(HIST("hV0Info"), 0.5);
          }
          if (aLambdaTag) {
            histos.fill(HIST("hV0Info"), 1.5);
          }
          if (lambdaTag && aLambdaTag) {
            doubleStatus.push_back(true);
            if (std::abs(v0.mLambda() - 1.1154) < std::abs(v0.mAntiLambda() - 1.1154)) {
              lambdaTag = true;
              aLambdaTag = false;
            } else {
              lambdaTag = false;
              aLambdaTag = true;
            }
          } else {
            doubleStatus.push_back(false);
          }
          if (lambdaTag) {
            histos.fill(HIST("hV0Info"), 2.5);
          }
          if (aLambdaTag) {
            histos.fill(HIST("hV0Info"), 3.5);
          }
          // LOGF(info, "v0 index2: (%d)", v0.index());
          auto postrack1 = v0.template posTrack_as<AllTrackCandidates>();
          auto negtrack1 = v0.template negTrack_as<AllTrackCandidates>();
          positiveIndex.push_back(postrack1.globalIndex());
          negativeIndex.push_back(negtrack1.globalIndex());
          v0Cospa.push_back(v0.v0cosPA());
          v0Radius.push_back(v0.v0radius());
          dcaPositive.push_back(std::abs(v0.dcapostopv()));
          dcaNegative.push_back(std::abs(v0.dcanegtopv()));
          dcaBetweenDaughter.push_back(std::abs(v0.dcaV0daughters()));
          if (lambdaTag) {
            v0Status.push_back(0);
            proton = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), o2::constants::physics::MassProton);
            antiPion = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), o2::constants::physics::MassPionCharged);
            lambda = proton + antiPion;
            lambdaMother.push_back(lambda);
            protonDaughter.push_back(proton);
            pionDaughter.push_back(antiPion);
            histos.fill(HIST("hLambdaMass"), lambda.M());
          } else if (aLambdaTag) {
            v0Status.push_back(1);
            antiProton = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), o2::constants::physics::MassProton);
            pion = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), o2::constants::physics::MassPionCharged);
            antiLambda = antiProton + pion;
            lambdaMother.push_back(antiLambda);
            protonDaughter.push_back(antiProton);
            pionDaughter.push_back(pion);
            histos.fill(HIST("hLambdaMass"), lambda.M());
          }
          numbV0 = numbV0 + 1;
        }
      }
      if (numbV0 > 1 && v0Cospa.size() > 1) {
        histos.fill(HIST("hEvtSelInfo"), 2.5);
        lambdaEvent(centrality, vz);
        auto indexEvent = lambdaEvent.lastIndex();
        //// Fill track table for V0//////////////////
        for (auto if1 = lambdaMother.begin(); if1 != lambdaMother.end(); ++if1) {
          auto i5 = std::distance(lambdaMother.begin(), if1);
          lambdaDummy = lambdaMother.at(i5);
          protonDummy = protonDaughter.at(i5);
          pionDummy = pionDaughter.at(i5);
          lambdaPair(indexEvent, v0Status.at(i5), doubleStatus.at(i5), v0Cospa.at(i5), v0Radius.at(i5), dcaPositive.at(i5), dcaNegative.at(i5), dcaBetweenDaughter.at(i5), lambdaDummy.Pt(), lambdaDummy.Eta(), lambdaDummy.Phi(), lambdaDummy.M(), protonDummy.Pt(), protonDummy.Eta(), protonDummy.Phi(), positiveIndex.at(i5), negativeIndex.at(i5));
        }
      }
    }
  }
  PROCESS_SWITCH(lambdaspincorrelation, processData, "Process data", false);

  void processMc(EventCandidates::iterator const& collision, AllTrackCandidates const&, ResoV0s const& V0s)
  {
    std::vector<ROOT::Math::PxPyPzMVector> lambdaMother, protonDaughter, pionDaughter;
    std::vector<int> v0Status = {};
    std::vector<bool> doubleStatus = {};
    std::vector<float> v0Cospa = {};
    std::vector<float> v0Radius = {};
    std::vector<float> dcaPositive = {};
    std::vector<float> dcaNegative = {};
    std::vector<int> positiveIndex = {};
    std::vector<int> negativeIndex = {};
    std::vector<float> dcaBetweenDaughter = {};
    int numbV0 = 0;
    // LOGF(info, "event collisions: (%d)", collision.index());
    auto centrality = collision.centFT0C();
    if (cfgEventTypepp)
      centrality = collision.centFT0M();
    auto vz = collision.posZ();
    int occupancy = collision.trackOccupancyInTimeRange();
    histos.fill(HIST("hEvtSelInfo"), 0.5);
    if ((rctCut.requireRCTFlagChecker && rctChecker(collision)) && collision.selection_bit(aod::evsel::kNoSameBunchPileup) && collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) && collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard) && collision.sel8() && collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll) && occupancy < cfgCutOccupancy) {
      histos.fill(HIST("hEvtSelInfo"), 1.5);
      for (const auto& v0 : V0s) {
        // LOGF(info, "v0 index 0 : (%d)", v0.index());
        auto [lambdaTag, aLambdaTag, isValid] = getLambdaTags(v0, collision);
        if (isValid) {
          // LOGF(info, "v0 index 1 : (%d)", v0.index());
          if (lambdaTag) {
            histos.fill(HIST("hV0Info"), 0.5);
          }
          if (aLambdaTag) {
            histos.fill(HIST("hV0Info"), 1.5);
          }
          if (lambdaTag && aLambdaTag) {
            doubleStatus.push_back(true);
            if (std::abs(v0.mLambda() - 1.1154) < std::abs(v0.mAntiLambda() - 1.1154)) {
              lambdaTag = true;
              aLambdaTag = false;
            } else {
              lambdaTag = false;
              aLambdaTag = true;
            }
          } else {
            doubleStatus.push_back(false);
          }
          if (lambdaTag) {
            histos.fill(HIST("hV0Info"), 2.5);
          }
          if (aLambdaTag) {
            histos.fill(HIST("hV0Info"), 3.5);
          }
          // LOGF(info, "v0 index2: (%d)", v0.index());
          auto postrack1 = v0.template posTrack_as<AllTrackCandidates>();
          auto negtrack1 = v0.template negTrack_as<AllTrackCandidates>();
          positiveIndex.push_back(postrack1.globalIndex());
          negativeIndex.push_back(negtrack1.globalIndex());
          v0Cospa.push_back(v0.v0cosPA());
          v0Radius.push_back(v0.v0radius());
          dcaPositive.push_back(std::abs(v0.dcapostopv()));
          dcaNegative.push_back(std::abs(v0.dcanegtopv()));
          dcaBetweenDaughter.push_back(std::abs(v0.dcaV0daughters()));
          if (lambdaTag) {
            v0Status.push_back(0);
            proton = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), o2::constants::physics::MassProton);
            antiPion = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), o2::constants::physics::MassPionCharged);
            lambda = proton + antiPion;
            lambdaMother.push_back(lambda);
            protonDaughter.push_back(proton);
            pionDaughter.push_back(antiPion);
            histos.fill(HIST("hLambdaMass"), lambda.M());
          } else if (aLambdaTag) {
            v0Status.push_back(1);
            antiProton = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), o2::constants::physics::MassProton);
            pion = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), o2::constants::physics::MassPionCharged);
            antiLambda = antiProton + pion;
            lambdaMother.push_back(antiLambda);
            protonDaughter.push_back(antiProton);
            pionDaughter.push_back(pion);
            histos.fill(HIST("hLambdaMass"), lambda.M());
          }
          numbV0 = numbV0 + 1;
        }
      }
      if (numbV0 > 1 && v0Cospa.size() > 1) {
        histos.fill(HIST("hEvtSelInfo"), 2.5);
        lambdaEventmc(centrality, vz);
        auto indexEvent = lambdaEventmc.lastIndex();
        //// Fill track table for V0//////////////////
        for (auto if1 = lambdaMother.begin(); if1 != lambdaMother.end(); ++if1) {
          auto i5 = std::distance(lambdaMother.begin(), if1);
          lambdaDummy = lambdaMother.at(i5);
          protonDummy = protonDaughter.at(i5);
          pionDummy = pionDaughter.at(i5);
          lambdaPairmc(indexEvent, v0Status.at(i5), doubleStatus.at(i5), v0Cospa.at(i5), v0Radius.at(i5), dcaPositive.at(i5), dcaNegative.at(i5), dcaBetweenDaughter.at(i5), lambdaDummy.Pt(), lambdaDummy.Eta(), lambdaDummy.Phi(), lambdaDummy.M(), protonDummy.Pt(), protonDummy.Eta(), protonDummy.Phi(), positiveIndex.at(i5), negativeIndex.at(i5));
        }
      }
    }
  }
  PROCESS_SWITCH(lambdaspincorrelation, processMc, "Process montecarlo", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<lambdaspincorrelation>(cfgc)};
}
