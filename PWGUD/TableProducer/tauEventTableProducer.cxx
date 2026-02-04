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
//
/// \file tauEventTableProducer.cxx
/// \brief Produces derived table from UD tables
///
/// \author Roman Lavicka <roman.lavicka@cern.ch>, Austrian Academy of Sciences & SMI
/// \since  09.04.2025
//

// C++ headers
#include <algorithm>
#include <random>
#include <set>
#include <utility>
#include <vector>

// O2 headers
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

// O2Physics headers
#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"
#include "PWGUD/DataModel/TauEventTables.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct TauEventTableProducer {
  Produces<o2::aod::TauTwoTracks> tauTwoTracks;
  Produces<o2::aod::TrueTauTwoTracks> trueTauTwoTracks;

  // Global varialbes
  Service<o2::framework::O2DatabasePDG> pdg;
  SGSelector sgSelector;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // declare configurables
  Configurable<bool> verboseInfo{"verboseInfo", false, {"Print general info to terminal; default it false."}};

  struct : ConfigurableGroup {
    Configurable<int> whichGapSide{"whichGapSide", 2, {"0 for side A, 1 for side C, 2 for both sides"}};
    Configurable<bool> useTrueGap{"useTrueGap", true, {"Calculate gapSide for a given FV0/FT0/ZDC thresholds"}};
    Configurable<int> cutNumContribs{"cutNumContribs", 2, {"How many contributors event has"}};
    Configurable<bool> useNumContribs{"useNumContribs", false, {"Use coll.numContribs as event cut"}};
    Configurable<int> cutRecoFlag{"cutRecoFlag", 1, {"0 = std mode, 1 = upc mode"}};
    Configurable<bool> useRecoFlag{"useRecoFlag", false, {"Use coll.flags as event cut"}};
    Configurable<int> cutRCTflag{"cutRCTflag", 0, {"0 = off, 1 = CBT, 2 = CBT+ZDC, 3 = CBThadron, 4 = CBThadron+ZDC"}};
    Configurable<float> cutTrueGapSideFV0{"cutTrueGapSideFV0", 180000, "FV0A threshold for SG selector"};
    Configurable<float> cutTrueGapSideFT0A{"cutTrueGapSideFT0A", 150., "FT0A threshold for SG selector"};
    Configurable<float> cutTrueGapSideFT0C{"cutTrueGapSideFT0C", 50., "FT0C threshold for SG selector"};
    Configurable<float> cutTrueGapSideZDC{"cutTrueGapSideZDC", 10000., "ZDC threshold for SG selector. 0 is <1n, 4.2 is <2n, 6.7 is <3n, 9.5 is <4n, 12.5 is <5n"};
    Configurable<float> cutFITtime{"cutFITtime", 40., "Maximum FIT time allowed. Default is 40ns"};
    Configurable<bool> cutEvTFb{"cutEvTFb", true, {"Event selection bit kNoTimeFrameBorder"}};
    Configurable<bool> cutEvITSROFb{"cutEvITSROFb", true, {"Event selection bit kNoITSROFrameBorder"}};
    Configurable<bool> cutEvSbp{"cutEvSbp", true, {"Event selection bit kNoSameBunchPileup"}};
    Configurable<bool> cutEvZvtxFT0vPV{"cutEvZvtxFT0vPV", false, {"Event selection bit kIsGoodZvtxFT0vsPV"}};
    Configurable<bool> cutEvVtxITSTPC{"cutEvVtxITSTPC", true, {"Event selection bit kIsVertexITSTPC"}};
    Configurable<float> cutEvOccupancy{"cutEvOccupancy", 100000., "Maximum allowed occupancy"};
    Configurable<bool> cutEvTrs{"cutEvTrs", false, {"Event selection bit kNoCollInTimeRangeStandard"}};
    Configurable<bool> cutEvTrofs{"cutEvTrofs", false, {"Event selection bit kNoCollInRofStandard"}};
    Configurable<bool> cutEvHmpr{"cutEvHmpr", false, {"Event selection bit kNoHighMultCollInPrevRof"}};
  } cutSample;

  struct : ConfigurableGroup {
    Configurable<bool> applyGlobalTrackSelection{"applyGlobalTrackSelection", false, {"Applies cut on here defined global tracks"}};
    Configurable<float> cutMinPt{"cutMinPt", 0.1f, {"Global track cut"}};
    Configurable<float> cutMaxPt{"cutMaxPt", 1e10f, {"Global track cut"}};
    Configurable<float> cutMinEta{"cutMinEta", -0.8f, {"Global track cut"}};
    Configurable<float> cutMaxEta{"cutMaxEta", 0.8f, {"Global track cut"}};
    Configurable<float> cutMaxDCAz{"cutMaxDCAz", 2.f, {"Global track cut"}};
    Configurable<float> cutMaxDCAxy{"cutMaxDCAxy", 1e10f, {"Global track cut"}};
    Configurable<bool> applyPtDependentDCAxy{"applyPtDependentDCAxy", false, {"Global track cut"}};
    Configurable<bool> cutHasITS{"cutHasITS", true, {"Global track cut"}};
    Configurable<int> cutMinITSnCls{"cutMinITSnCls", 1, {"Global track cut"}};
    Configurable<float> cutMaxITSchi2{"cutMaxITSchi2", 36.f, {"Global track cut"}};
    Configurable<int> cutITShitsRule{"cutITShitsRule", 0, {"Global track cut"}};
    Configurable<bool> cutHasTPC{"cutHasTPC", true, {"Global track cut"}};
    Configurable<int> cutMinTPCnCls{"cutMinTPCnCls", 1, {"Global track cut"}};
    Configurable<int> cutMinTPCnClsXrows{"cutMinTPCnClsXrows", 70, {"Global track cut"}};
    Configurable<float> cutMinTPCnClsXrowsOverNcls{"cutMinTPCnClsXrowsOverNcls", 0.8f, {"Global track cut"}};
    Configurable<float> cutMaxTPCchi2{"cutMaxTPCchi2", 4.f, {"Global track cut"}};
    Configurable<bool> cutGoodITSTPCmatching{"cutGoodITSTPCmatching", true, {"Global track cut"}};
    Configurable<float> cutMaxTOFchi2{"cutMaxTOFchi2", 3.f, {"Global track cut"}};
  } cutGlobalTrack;

  struct : ConfigurableGroup {
    Configurable<bool> preselUseTrackPID{"preselUseTrackPID", true, {"Apply weak PID check on tracks."}};
    Configurable<int> preselNgoodPVtracs{"preselNgoodPVtracs", 2, {"How many good PV tracks to select."}};
    Configurable<float> preselMinElectronNsigmaEl{"preselMinElectronNsigmaEl", 4.0, {"Good el candidate hypo in. Upper n sigma cut on el hypo of selected electron. What is more goes away."}};
    Configurable<float> preselMaxElectronNsigmaEl{"preselMaxElectronNsigmaEl", -2.0, {"Good el candidate hypo in. Lower n sigma cut on el hypo of selected electron. What is less goes away."}};
    Configurable<bool> preselElectronHasTOF{"preselElectronHasTOF", true, {"Electron candidated is required to hit TOF."}};
    Configurable<float> preselMinPionNsigmaEl{"preselMinPionNsigmaEl", 5.0, {"Good pi candidate hypo in. Upper n sigma cut on pi hypo of selected electron. What is more goes away."}};
    Configurable<float> preselMaxPionNsigmaEl{"preselMaxPionNsigmaEl", -5.0, {"Good pi candidate hypo in. Lower n sigma cut on pi hypo of selected electron. What is less goes away."}};
    Configurable<float> preselMinMuonNsigmaEl{"preselMinMuonNsigmaEl", 5.0, {"Good pi candidate hypo in. Upper n sigma cut on pi hypo of selected electron. What is more goes away."}};
    Configurable<float> preselMaxMuonNsigmaEl{"preselMaxMuonNsigmaEl", -5.0, {"Good pi candidate hypo in. Lower n sigma cut on pi hypo of selected electron. What is less goes away."}};
    Configurable<bool> preselMupionHasTOF{"preselMupionHasTOF", true, {"Mupion candidate is required to hit TOF."}};
  } cutPreselect;

  using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags>;
  using FullSGUDCollisions = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDCollisionSelExtras, aod::SGCollisions, aod::UDZdcsReduced>;
  using FullSGUDCollision = FullSGUDCollisions::iterator;
  using FullMCUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags, aod::UDMcTrackLabels>;
  using FullMCSGUDCollisions = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDCollisionSelExtras, aod::SGCollisions, aod::UDMcCollsLabels>;
  using FullMCSGUDCollision = FullMCSGUDCollisions::iterator;

  // init
  void init(InitContext&)
  {
    if (verboseInfo)
      printMediumMessage("INIT METHOD");

    mySetITShitsRule(cutGlobalTrack.cutITShitsRule);

    histos.add("Truth/hTroubles", "Counter of unwanted issues;;Number of  troubles (-)", HistType::kTH1D, {{15, 0.5, 15.5}});

  } // end init

  template <typename C>
  bool isGoodFITtime(C const& coll, float maxFITtime)
  {

    // FTOA
    if ((std::abs(coll.timeFT0A()) > maxFITtime) && coll.timeFT0A() > -998.)
      return false;

    // FTOC
    if ((std::abs(coll.timeFT0C()) > maxFITtime) && coll.timeFT0C() > -998.)
      return false;

    return true;
  }

  template <typename C>
  bool isGoodRCTflag(C const& coll)
  {
    switch (cutSample.cutRCTflag) {
      case 1:
        return sgSelector.isCBTOk(coll);
      case 2:
        return sgSelector.isCBTZdcOk(coll);
      case 3:
        return sgSelector.isCBTHadronOk(coll);
      case 4:
        return sgSelector.isCBTHadronZdcOk(coll);
      default:
        return true;
    }
  }

  template <typename C>
  bool isGoodROFtime(C const& coll)
  {

    // kNoTimeFrameBorder
    if (cutSample.cutEvTFb && !coll.tfb())
      return false;

    // kNoITSROFrameBorder
    if (cutSample.cutEvITSROFb && !coll.itsROFb())
      return false;

    // kNoSameBunchPileup
    if (cutSample.cutEvSbp && !coll.sbp())
      return false;

    // kIsGoodZvtxFT0vsPV
    if (cutSample.cutEvZvtxFT0vPV && !coll.zVtxFT0vPV())
      return false;

    // kIsVertexITSTPC
    if (cutSample.cutEvVtxITSTPC && !coll.vtxITSTPC())
      return false;

    // Occupancy
    if (coll.occupancyInTime() > cutSample.cutEvOccupancy)
      return false;

    // kNoCollInTimeRangeStandard
    if (cutSample.cutEvTrs && !coll.trs())
      return false;

    // kNoCollInRofStandard
    if (cutSample.cutEvTrofs && !coll.trofs())
      return false;

    // kNoHighMultCollInPrevRof
    if (cutSample.cutEvHmpr && !coll.hmpr())
      return false;

    return true;
  }

  std::vector<std::pair<int8_t, std::set<uint8_t>>> cutMyRequiredITSHits{};

  void mySetRequireHitsInITSLayers(int8_t minNRequiredHits, std::set<uint8_t> requiredLayers)
  {
    // layer 0 corresponds to the the innermost ITS layer
    cutMyRequiredITSHits.push_back(std::make_pair(minNRequiredHits, requiredLayers));
  }

  void mySetITShitsRule(int matching)
  {
    switch (matching) {
      case 0: // Run3ITSibAny
        mySetRequireHitsInITSLayers(1, {0, 1, 2});
        break;
      case 1: // Run3ITSibTwo
        mySetRequireHitsInITSLayers(2, {0, 1, 2});
        break;
      case 2: // Run3ITSallAny
        mySetRequireHitsInITSLayers(1, {0, 1, 2, 3, 4, 5, 6});
        break;
      case 3: // Run3ITSall7Layers
        mySetRequireHitsInITSLayers(7, {0, 1, 2, 3, 4, 5, 6});
        break;
      default:
        LOG(fatal) << "You chose wrong ITS matching";
        break;
    }
  }

  bool isFulfillsITSHitRequirementsReinstatement(uint8_t itsClusterMap) const
  {
    constexpr uint8_t kBit = 1;
    for (const auto& kITSrequirement : cutMyRequiredITSHits) {
      auto hits = std::count_if(kITSrequirement.second.begin(), kITSrequirement.second.end(), [&](auto&& requiredLayer) { return itsClusterMap & (kBit << requiredLayer); });
      if ((kITSrequirement.first == -1) && (hits > 0)) {
        return false; // no hits were required in specified layers
      } else if (hits < kITSrequirement.first) {
        return false; // not enough hits found in specified layers
      }
    }
    return true;
  }

  template <typename T>
  bool isGlobalTrackReinstatement(T const& track)
  {
    // kInAcceptance copy
    if (track.pt() < cutGlobalTrack.cutMinPt || track.pt() > cutGlobalTrack.cutMaxPt)
      return false;
    if (eta(track.px(), track.py(), track.pz()) < cutGlobalTrack.cutMinEta || eta(track.px(), track.py(), track.pz()) > cutGlobalTrack.cutMaxEta)
      return false;
    // kPrimaryTracks
    // GoldenChi2 cut is only for Run 2
    if (std::abs(track.dcaZ()) > cutGlobalTrack.cutMaxDCAz)
      return false;
    if (cutGlobalTrack.applyPtDependentDCAxy) {
      float maxDCA = 0.0182f + 0.0350f / std::pow(track.pt(), 1.01f);
      if (std::abs(track.dcaXY()) > maxDCA)
        return false;
    } else {
      if (std::abs(track.dcaXY()) > cutGlobalTrack.cutMaxDCAxy)
        return false;
    }
    // kQualityTrack
    // TrackType is always 1 as per definition of processed Run3 AO2Ds
    // ITS
    if (cutGlobalTrack.cutHasITS && !track.hasITS())
      return false; // ITS refit
    if (track.itsNCls() < cutGlobalTrack.cutMinITSnCls)
      return false;
    if (track.itsChi2NCl() > cutGlobalTrack.cutMaxITSchi2)
      return false;
    if (!isFulfillsITSHitRequirementsReinstatement(track.itsClusterMap()))
      return false;
    //  TPC
    if (cutGlobalTrack.cutHasTPC && !track.hasTPC())
      return false; // TPC refit
    if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < cutGlobalTrack.cutMinTPCnCls)
      return false; // tpcNClsFound()
    if (track.tpcNClsCrossedRows() < cutGlobalTrack.cutMinTPCnClsXrows)
      return false;
    if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < cutGlobalTrack.cutMinTPCnClsXrowsOverNcls)
      return false;
    if (track.tpcChi2NCl() > cutGlobalTrack.cutMaxTPCchi2)
      return false; // TPC chi2
    if (cutGlobalTrack.cutGoodITSTPCmatching) {
      if (track.itsChi2NCl() < 0.)
        return false; // good ITS-TPC matching means ITS ch2 is not negative
    }
    //  TOF
    if (track.hasTOF()) {
      if (track.tpcChi2NCl() > cutGlobalTrack.cutMaxTOFchi2)
        return false; // TOF chi2
    }

    return true;
  }

  template <typename T>
  bool isElectronCandidate(T const& electronCandidate)
  // Loose criterium to find electron-like particle
  // Requiring TOF to avoid double-counting pions/electrons and for better timing
  {
    if (electronCandidate.tpcNSigmaEl() < cutPreselect.preselMaxElectronNsigmaEl || electronCandidate.tpcNSigmaEl() > cutPreselect.preselMinElectronNsigmaEl)
      return false;
    if (cutPreselect.preselElectronHasTOF && !electronCandidate.hasTOF())
      return false;
    return true;
  }

  template <typename T>
  bool isMuPionCandidate(T const& muPionCandidate)
  // Loose criterium to find muon/pion-like particle
  // Requiring TOF for better timing
  {
    if (muPionCandidate.tpcNSigmaMu() < cutPreselect.preselMaxMuonNsigmaEl || muPionCandidate.tpcNSigmaMu() > cutPreselect.preselMinMuonNsigmaEl)
      return false;
    if (muPionCandidate.tpcNSigmaPi() < cutPreselect.preselMaxPionNsigmaEl || muPionCandidate.tpcNSigmaPi() > cutPreselect.preselMinPionNsigmaEl)
      return false;
    if (cutPreselect.preselMupionHasTOF && !muPionCandidate.hasTOF())
      return false;
    return true;
  }

  void processDataSG(FullSGUDCollision const& collision,
                     FullUDTracks const& tracks)
  {

    if (!isGoodRCTflag(collision))
      return;

    if (!isGoodROFtime(collision))
      return;

    int gapSide = collision.gapSide();
    int trueGapSide = sgSelector.trueGap(collision, cutSample.cutTrueGapSideFV0, cutSample.cutTrueGapSideFT0A, cutSample.cutTrueGapSideFT0C, cutSample.cutTrueGapSideZDC);
    if (cutSample.useTrueGap)
      gapSide = trueGapSide;
    if (gapSide != cutSample.whichGapSide)
      return;

    if (!isGoodFITtime(collision, cutSample.cutFITtime))
      return;

    if (cutSample.useNumContribs && (collision.numContrib() != cutSample.cutNumContribs))
      return;

    if (cutSample.useRecoFlag && (collision.flags() != cutSample.cutRecoFlag))
      return;

    int countTracksPerCollision = 0;
    int countGoodNonPVtracks = 0;
    int countGoodPVtracks = 0;
    std::vector<int> vecTrkIdx;
    // Loop over tracks with selections
    for (const auto& track : tracks) {
      countTracksPerCollision++;
      if (!isGlobalTrackReinstatement(track))
        continue;
      if (!track.isPVContributor()) {
        countGoodNonPVtracks++;
        continue;
      }
      countGoodPVtracks++;
      vecTrkIdx.push_back(track.index());
    } // Loop over tracks with selections

    // Apply weak condition on track PID
    int countPVGTel = 0;
    int countPVGTmupi = 0;
    if (countGoodPVtracks == 2) {
      for (const auto& vecMember : vecTrkIdx) {
        const auto& thisTrk = tracks.iteratorAt(vecMember);
        if (isElectronCandidate(thisTrk)) {
          countPVGTel++;
          continue;
        }
        if (isMuPionCandidate(thisTrk)) {
          countPVGTmupi++;
        }
      }
    }

    if (cutPreselect.preselUseTrackPID ? ((countPVGTel == 2 && countPVGTmupi == 0) || (countPVGTel == 1 && countPVGTmupi == 1)) : countGoodPVtracks == cutPreselect.preselNgoodPVtracs) {
      const auto& trk1 = tracks.iteratorAt(vecTrkIdx[0]);
      const auto& trk2 = tracks.iteratorAt(vecTrkIdx[1]);

      float px[2] = {trk1.px(), trk2.px()};
      float py[2] = {trk1.py(), trk2.py()};
      float pz[2] = {trk1.pz(), trk2.pz()};
      int sign[2] = {trk1.sign(), trk2.sign()};
      float dcaxy[2] = {trk1.dcaXY(), trk2.dcaXY()};
      float dcaz[2] = {trk1.dcaZ(), trk2.dcaZ()};
      float trkTimeRes[2] = {trk1.trackTimeRes(), trk2.trackTimeRes()};
      uint32_t itsClusterSizesTrk1 = trk1.itsClusterSizes();
      uint32_t itsClusterSizesTrk2 = trk2.itsClusterSizes();
      float tpcSignal[2] = {trk1.tpcSignal(), trk2.tpcSignal()};
      float tpcEl[2] = {trk1.tpcNSigmaEl(), trk2.tpcNSigmaEl()};
      float tpcMu[2] = {trk1.tpcNSigmaMu(), trk2.tpcNSigmaMu()};
      float tpcPi[2] = {trk1.tpcNSigmaPi(), trk2.tpcNSigmaPi()};
      float tpcKa[2] = {trk1.tpcNSigmaKa(), trk2.tpcNSigmaKa()};
      float tpcPr[2] = {trk1.tpcNSigmaPr(), trk2.tpcNSigmaPr()};
      float tpcIP[2] = {trk1.tpcInnerParam(), trk2.tpcInnerParam()};
      float tofSignal[2] = {trk1.tofSignal(), trk2.tofSignal()};
      float tofEl[2] = {trk1.tofNSigmaEl(), trk2.tofNSigmaEl()};
      float tofMu[2] = {trk1.tofNSigmaMu(), trk2.tofNSigmaMu()};
      float tofPi[2] = {trk1.tofNSigmaPi(), trk2.tofNSigmaPi()};
      float tofKa[2] = {trk1.tofNSigmaKa(), trk2.tofNSigmaKa()};
      float tofPr[2] = {trk1.tofNSigmaPr(), trk2.tofNSigmaPr()};
      float tofEP[2] = {trk1.tofExpMom(), trk2.tofExpMom()};
      //      float infoZDC[4] = {-999., -999., -999., -999.};
      //      if constexpr (requires { collision.udZdcsReduced(); }) {
      //        infoZDC[0] = collision.energyCommonZNA();
      //        infoZDC[1] = collision.energyCommonZNC();
      //        infoZDC[2] = collision.timeZNA();
      //        infoZDC[3] = collision.timeZNC();
      //      }
      float infoZDC[4] = {collision.energyCommonZNA(), collision.energyCommonZNC(), collision.timeZNA(), collision.timeZNC()};

      tauTwoTracks(collision.runNumber(), collision.globalBC(), countTracksPerCollision, collision.numContrib(), countGoodNonPVtracks, collision.posX(), collision.posY(), collision.posZ(),
                   collision.flags(), collision.occupancyInTime(), collision.hadronicRate(), collision.trs(), collision.trofs(), collision.hmpr(),
                   collision.tfb(), collision.itsROFb(), collision.sbp(), collision.zVtxFT0vPV(), collision.vtxITSTPC(),
                   collision.totalFT0AmplitudeA(), collision.totalFT0AmplitudeC(), collision.totalFV0AmplitudeA(), infoZDC[0], infoZDC[1],
                   collision.timeFT0A(), collision.timeFT0C(), collision.timeFV0A(), infoZDC[2], infoZDC[3],
                   px, py, pz, sign, dcaxy, dcaz, trkTimeRes,
                   itsClusterSizesTrk1, itsClusterSizesTrk2,
                   tpcSignal, tpcEl, tpcMu, tpcPi, tpcKa, tpcPr, tpcIP,
                   tofSignal, tofEl, tofMu, tofPi, tofKa, tofPr, tofEP);
    }
  }
  PROCESS_SWITCH(TauEventTableProducer, processDataSG, "Iterate UD tables with measured data created by SG-Candidate-Producer.", false);

  PresliceUnsorted<aod::UDMcParticles> partPerMcCollision = aod::udmcparticle::udMcCollisionId;
  PresliceUnsorted<FullMCSGUDCollisions> colPerMcCollision = aod::udcollision::udMcCollisionId;
  PresliceUnsorted<FullMCUDTracks> trackPerMcParticle = aod::udmctracklabel::udMcParticleId;
  Preslice<FullMCUDTracks> trackPerCollision = aod::udtrack::udCollisionId; // sorted preslice used because the pair track-collision is already sorted in processDataSG function

  void processMonteCarlo(aod::UDMcCollisions const& mccollisions,
                         aod::UDMcParticles const& parts,
                         FullMCSGUDCollisions const& recolls,
                         FullMCUDTracks const& trks)
  {
    // start loop over generated collisions
    for (const auto& mccoll : mccollisions) {

      // prepare local variables for output table
      int32_t runNumber = -999;
      int bc = -999;
      int nTrks[3] = {-999, -999, -999}; // totalTracks, numContrib, globalNonPVtracks
      float vtxPos[3] = {-999., -999., -999.};
      int recoMode = -999;
      int occupancy = -999.;
      double hadronicRate = -999.;
      int bcSels[8] = {-999, -999, -999, -999, -999, -999, -999, -999};
      float amplitudesFIT[3] = {-999., -999., -999.}; // FT0A, FT0C, FV0
      float timesFIT[3] = {-999., -999., -999.};      // FT0A, FT0C, FV0

      float px[2] = {-999., -999.};
      float py[2] = {-999., -999.};
      float pz[2] = {-999., -999.};
      int sign[2] = {-999, -999};
      float dcaxy[2] = {-999., -999.};
      float dcaz[2] = {-999., -999.};
      float trkTimeRes[2] = {-999., -999.};
      uint32_t itsClusterSizesTrk1 = 4294967295;
      uint32_t itsClusterSizesTrk2 = 4294967295;
      float tpcSignal[2] = {-999, -999};
      float tpcEl[2] = {-999, -999};
      float tpcMu[2] = {-999, -999};
      float tpcPi[2] = {-999, -999};
      float tpcKa[2] = {-999, -999};
      float tpcPr[2] = {-999, -999};
      float tpcIP[2] = {-999, -999};
      float tofSignal[2] = {-999, -999};
      float tofEl[2] = {-999, -999};
      float tofMu[2] = {-999, -999};
      float tofPi[2] = {-999, -999};
      float tofKa[2] = {-999, -999};
      float tofPr[2] = {-999, -999};
      float tofEP[2] = {-999, -999};

      int trueChannel = -1;
      bool trueHasRecoColl = false;
      float trueTauX[2] = {-999., -999.};
      float trueTauY[2] = {-999., -999.};
      float trueTauZ[2] = {-999., -999.};
      float trueDaugX[2] = {-999., -999.};
      float trueDaugY[2] = {-999., -999.};
      float trueDaugZ[2] = {-999., -999.};
      int trueDaugPdgCode[2] = {-999, -999};
      bool problem = false;

      // find reconstructed collisions associated to the generated collision
      auto const& collFromMcColls = recolls.sliceBy(colPerMcCollision, mccoll.globalIndex());
      // check the generated collision was reconstructed
      if (collFromMcColls.size() > 0) { // get the truth and reco-level info
        trueHasRecoColl = true;
        // check there is exactly one reco-level collision associated to generated collision
        if (collFromMcColls.size() > 1) {
          if (verboseInfo)
            printLargeMessage("Truth collision has more than 1 reco collision. Skipping this event.");
          histos.get<TH1>(HIST("Truth/hTroubles"))->Fill(1);
          problem = true;
          continue;
        }
        // grap reco-level collision
        auto const& collFromMcColl = collFromMcColls.iteratorAt(0);
        // grab tracks from the reco-level collision to get info to match measured data tables (processDataSG function)
        auto const& trksFromColl = trks.sliceBy(trackPerCollision, collFromMcColl.globalIndex());
        int countTracksPerCollision = 0;
        int countGoodNonPVtracks = 0;
        for (auto const& trkFromColl : trksFromColl) {
          countTracksPerCollision++;
          if (!trkFromColl.isPVContributor()) {
            countGoodNonPVtracks++;
            continue;
          }
        }

        // fill info for reconstructed collision
        runNumber = collFromMcColl.runNumber();
        bc = collFromMcColl.globalBC();
        nTrks[0] = countTracksPerCollision;
        nTrks[1] = collFromMcColl.numContrib();
        nTrks[2] = countGoodNonPVtracks;
        vtxPos[0] = collFromMcColl.posX();
        vtxPos[1] = collFromMcColl.posY();
        vtxPos[2] = collFromMcColl.posZ();
        recoMode = collFromMcColl.flags();
        occupancy = collFromMcColl.occupancyInTime();
        hadronicRate = collFromMcColl.hadronicRate();
        bcSels[0] = collFromMcColl.trs();
        bcSels[1] = collFromMcColl.trofs();
        bcSels[2] = collFromMcColl.hmpr();
        bcSels[3] = collFromMcColl.tfb();
        bcSels[4] = collFromMcColl.itsROFb();
        bcSels[5] = collFromMcColl.sbp();
        bcSels[6] = collFromMcColl.zVtxFT0vPV();
        bcSels[7] = collFromMcColl.vtxITSTPC();
        amplitudesFIT[0] = collFromMcColl.totalFT0AmplitudeA();
        amplitudesFIT[1] = collFromMcColl.totalFT0AmplitudeC();
        amplitudesFIT[2] = collFromMcColl.totalFV0AmplitudeA();
        timesFIT[0] = collFromMcColl.timeFT0A();
        timesFIT[1] = collFromMcColl.timeFT0C();
        timesFIT[2] = collFromMcColl.timeFV0A();

        // get particles associated to generated collision
        auto const& partsFromMcColl = parts.sliceBy(partPerMcCollision, mccoll.globalIndex());
        int countMothers = 0;
        for (const auto& particle : partsFromMcColl) {
          // select only tauons with checking if particle has no mother
          if (particle.has_mothers())
            continue;
          countMothers++;
          // check the generated collision does not have more than 2 tauons
          if (countMothers > 2) {
            if (verboseInfo)
              printLargeMessage("Truth collision has more than 2 no mother particles. Breaking the particle loop.");
            histos.get<TH1>(HIST("Truth/hTroubles"))->Fill(2);
            problem = true;
            break;
          }
          // fill info for each tau
          trueTauX[countMothers - 1] = particle.px();
          trueTauY[countMothers - 1] = particle.py();
          trueTauZ[countMothers - 1] = particle.pz();

          // get daughters of the tau
          const auto& daughters = particle.daughters_as<aod::UDMcParticles>();
          int countDaughters = 0;
          for (const auto& daughter : daughters) {
            // check if it is the charged particle (= no pi0 or neutrino)
            if (enumMyParticle(daughter.pdgCode()) == -1)
              continue;
            countDaughters++;
            // check there is only 1 charged daughter related to 1 tau
            if (countDaughters > 1) {
              if (verboseInfo)
                printLargeMessage("Truth collision has more than 1 charged daughters of no mother particles. Breaking the daughter loop.");
              histos.get<TH1>(HIST("Truth/hTroubles"))->Fill(3);
              problem = true;
              break;
            }
            // fill info for each daughter
            trueDaugX[countMothers - 1] = daughter.px();
            trueDaugY[countMothers - 1] = daughter.py();
            trueDaugZ[countMothers - 1] = daughter.pz();
            trueDaugPdgCode[countMothers - 1] = daughter.pdgCode();

            // get tracks associated to MC daughter (how well the daughter was reconstructed)
            auto const& tracksFromDaughter = trks.sliceBy(trackPerMcParticle, daughter.globalIndex());
            // check there is exactly 1 track per 1 particle
            if (tracksFromDaughter.size() > 1) {
              if (verboseInfo)
                printLargeMessage("Daughter has more than 1 associated track. Skipping this daughter.");
              histos.get<TH1>(HIST("Truth/hTroubles"))->Fill(4);
              problem = true;
              continue;
            }
            // grab the track and fill info for reconstructed track (should be done twice)
            const auto& trk = tracksFromDaughter.iteratorAt(0);
            px[countMothers - 1] = trk.px();
            py[countMothers - 1] = trk.py();
            pz[countMothers - 1] = trk.pz();
            sign[countMothers - 1] = trk.sign();
            dcaxy[countMothers - 1] = trk.dcaXY();
            dcaz[countMothers - 1] = trk.dcaZ();
            trkTimeRes[countMothers - 1] = trk.trackTimeRes();
            if (countMothers == 1) {
              itsClusterSizesTrk1 = trk.itsClusterSizes();
            } else {
              itsClusterSizesTrk2 = trk.itsClusterSizes();
            }
            tpcSignal[countMothers - 1] = trk.tpcSignal();
            tpcEl[countMothers - 1] = trk.tpcNSigmaEl();
            tpcMu[countMothers - 1] = trk.tpcNSigmaMu();
            tpcPi[countMothers - 1] = trk.tpcNSigmaPi();
            tpcKa[countMothers - 1] = trk.tpcNSigmaKa();
            tpcPr[countMothers - 1] = trk.tpcNSigmaPr();
            tpcIP[countMothers - 1] = trk.tpcInnerParam();
            tofSignal[countMothers - 1] = trk.tofSignal();
            tofEl[countMothers - 1] = trk.tofNSigmaEl();
            tofMu[countMothers - 1] = trk.tofNSigmaMu();
            tofPi[countMothers - 1] = trk.tofNSigmaPi();
            tofKa[countMothers - 1] = trk.tofNSigmaKa();
            tofPr[countMothers - 1] = trk.tofNSigmaPr();
            tofEP[countMothers - 1] = trk.tofExpMom();
          } // daughters
        } // particles
      } else { // get only the truth information. The reco-level info is left on default
        // get particles associated to generated collision
        auto const& partsFromMcColl = parts.sliceBy(partPerMcCollision, mccoll.globalIndex());
        int countMothers = 0;
        for (const auto& particle : partsFromMcColl) {
          // select only tauons with checking if particle has no mother
          if (particle.has_mothers())
            continue;
          countMothers++;
          // check the generated collision does not have more than 2 tauons
          if (countMothers > 2) {
            if (verboseInfo)
              printLargeMessage("Truth collision has more than 2 no mother particles. Breaking the particle loop.");
            histos.get<TH1>(HIST("Truth/hTroubles"))->Fill(12);
            problem = true;
            break;
          }
          // fill info for each tau
          trueTauX[countMothers - 1] = particle.px();
          trueTauY[countMothers - 1] = particle.py();
          trueTauZ[countMothers - 1] = particle.pz();

          // get daughters of the tau
          const auto& daughters = particle.daughters_as<aod::UDMcParticles>();
          int countDaughters = 0;
          for (const auto& daughter : daughters) {
            // select only the charged particle (= no pi0 or neutrino)
            if (enumMyParticle(daughter.pdgCode()) == -1)
              continue;
            countDaughters++;
            // check there is only 1 charged daughter related to 1 tau
            if (countDaughters > 1) {
              if (verboseInfo)
                printLargeMessage("Truth collision has more than 1 charged daughters of no mother particles. Breaking the daughter loop.");
              histos.get<TH1>(HIST("Truth/hTroubles"))->Fill(13);
              problem = true;
              break;
            }
            // fill info for each daughter
            trueDaugX[countMothers - 1] = daughter.px();
            trueDaugY[countMothers - 1] = daughter.py();
            trueDaugZ[countMothers - 1] = daughter.pz();
            trueDaugPdgCode[countMothers - 1] = daughter.pdgCode();
          } // daughters
        } // particles
      } // collisions

      // decide the channel and set the variable. Only two cahnnels suported now.
      if ((enumMyParticle(trueDaugPdgCode[0]) == P_ELECTRON) && (enumMyParticle(trueDaugPdgCode[1]) == P_ELECTRON))
        trueChannel = CH_EE;
      if ((enumMyParticle(trueDaugPdgCode[0]) == P_ELECTRON) && ((enumMyParticle(trueDaugPdgCode[1]) == P_PION) || (enumMyParticle(trueDaugPdgCode[1]) == P_MUON)))
        trueChannel = CH_EMUPI;
      if ((enumMyParticle(trueDaugPdgCode[1]) == P_ELECTRON) && ((enumMyParticle(trueDaugPdgCode[0]) == P_PION) || (enumMyParticle(trueDaugPdgCode[0]) == P_MUON)))
        trueChannel = CH_EMUPI;

      trueTauTwoTracks(runNumber, bc, nTrks[0], nTrks[1], nTrks[2], vtxPos[0], vtxPos[1], vtxPos[2],
                       recoMode, occupancy, hadronicRate, bcSels[0], bcSels[1], bcSels[2],
                       bcSels[3], bcSels[4], bcSels[5], bcSels[6], bcSels[7],
                       amplitudesFIT[0], amplitudesFIT[1], amplitudesFIT[2], -999., -999., // no ZDC info in MC
                       timesFIT[0], timesFIT[1], timesFIT[2], -999., -999.,                // no ZDC info in MC
                       px, py, pz, sign, dcaxy, dcaz, trkTimeRes,
                       itsClusterSizesTrk1, itsClusterSizesTrk2,
                       tpcSignal, tpcEl, tpcMu, tpcPi, tpcKa, tpcPr, tpcIP,
                       tofSignal, tofEl, tofMu, tofPi, tofKa, tofPr, tofEP,
                       trueChannel, trueHasRecoColl, mccoll.posX(), mccoll.posY(), mccoll.posZ(),
                       trueTauX, trueTauY, trueTauZ, trueDaugX, trueDaugY, trueDaugZ, trueDaugPdgCode, problem);
    } // mccollisions
  }
  PROCESS_SWITCH(TauEventTableProducer, processMonteCarlo, "Iterate UD tables with simulated data created by SG-Candidate-Producer.", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TauEventTableProducer>(cfgc)};
}
