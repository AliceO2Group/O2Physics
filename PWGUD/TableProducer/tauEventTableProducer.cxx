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
#include <set>
#include <utility>
#include <algorithm>
#include <vector>
#include <random>

// O2 headers
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

// O2Physics headers
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/DataModel/UDIndex.h" // for UDMcParticles2UDTracks table
#include "PWGUD/Core/SGSelector.h"

// ROOT headers
#include "TLorentzVector.h"
#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

namespace o2::aod
{
namespace tau_tree
{
// event info
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(Bc, bc, int);
DECLARE_SOA_COLUMN(TotalTracks, totalTracks, int);
DECLARE_SOA_COLUMN(NumContrib, numContrib, int);
DECLARE_SOA_COLUMN(GlobalNonPVtracks, globalNonPVtracks, int);
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(RecoMode, recoMode, int);
DECLARE_SOA_COLUMN(OccupancyInTime, occupancyInTime, int);
DECLARE_SOA_COLUMN(HadronicRate, hadronicRate, double);
DECLARE_SOA_COLUMN(Trs, trs, int);
DECLARE_SOA_COLUMN(Trofs, trofs, int);
DECLARE_SOA_COLUMN(Hmpr, hmpr, int);
DECLARE_SOA_COLUMN(Tfb, tfb, int);
DECLARE_SOA_COLUMN(ItsRofb, itsRofb, int);
DECLARE_SOA_COLUMN(Sbp, sbp, int);
DECLARE_SOA_COLUMN(ZvtxFT0vsPv, zvtxFT0vsPv, int);
DECLARE_SOA_COLUMN(VtxITSTPC, vtxITSTPC, int);
// FIT info
DECLARE_SOA_COLUMN(TotalFT0AmplitudeA, totalFT0AmplitudeA, float);
DECLARE_SOA_COLUMN(TotalFT0AmplitudeC, totalFT0AmplitudeC, float);
DECLARE_SOA_COLUMN(TotalFV0AmplitudeA, totalFV0AmplitudeA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNA, energyCommonZNA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNC, energyCommonZNC, float);
DECLARE_SOA_COLUMN(TimeFT0A, timeFT0A, float);
DECLARE_SOA_COLUMN(TimeFT0C, timeFT0C, float);
DECLARE_SOA_COLUMN(TimeFV0A, timeFV0A, float);
DECLARE_SOA_COLUMN(TimeZNA, timeZNA, float);
DECLARE_SOA_COLUMN(TimeZNC, timeZNC, float);
// tracks
DECLARE_SOA_COLUMN(TrkPx, trkPx, float[2]);
DECLARE_SOA_COLUMN(TrkPy, trkPy, float[2]);
DECLARE_SOA_COLUMN(TrkPz, trkPz, float[2]);
DECLARE_SOA_COLUMN(TrkSign, trkSign, int[2]);
DECLARE_SOA_COLUMN(TrkDCAxy, trkDCAxy, float[2]);
DECLARE_SOA_COLUMN(TrkDCAz, trkDCAz, float[2]);
DECLARE_SOA_COLUMN(TrkTimeRes, trkTimeRes, float[2]);
DECLARE_SOA_COLUMN(Trk1ITSclusterSizes, trk1ITSclusterSizes, uint32_t);
DECLARE_SOA_COLUMN(Trk2ITSclusterSizes, trk2ITSclusterSizes, uint32_t);
DECLARE_SOA_COLUMN(TrkTPCsignal, trkTPCsignal, float[2]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaEl, trkTPCnSigmaEl, float[2]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaMu, trkTPCnSigmaMu, float[2]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaPi, trkTPCnSigmaPi, float[2]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaKa, trkTPCnSigmaKa, float[2]);
DECLARE_SOA_COLUMN(TrkTPCnSigmaPr, trkTPCnSigmaPr, float[2]);
DECLARE_SOA_COLUMN(TrkTPCinnerParam, trkTPCinnerParam, float[2]);
DECLARE_SOA_COLUMN(TrkTOFsignal, trkTOFsignal, float[2]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaEl, trkTOFnSigmaEl, float[2]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaMu, trkTOFnSigmaMu, float[2]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaPi, trkTOFnSigmaPi, float[2]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaKa, trkTOFnSigmaKa, float[2]);
DECLARE_SOA_COLUMN(TrkTOFnSigmaPr, trkTOFnSigmaPr, float[2]);
DECLARE_SOA_COLUMN(TrkTOFexpMom, trkTOFexpMom, float[2]);

} // namespace tau_tree
DECLARE_SOA_TABLE(TauTwoTracks, "AOD", "TAUTWOTRACK",
                  tau_tree::RunNumber, tau_tree::Bc, tau_tree::TotalTracks, tau_tree::NumContrib, tau_tree::GlobalNonPVtracks, tau_tree::PosX, tau_tree::PosY, tau_tree::PosZ,
                  tau_tree::RecoMode, tau_tree::OccupancyInTime, tau_tree::HadronicRate,
                  tau_tree::Trs, tau_tree::Trofs, tau_tree::Hmpr, tau_tree::Tfb, tau_tree::ItsRofb, tau_tree::Sbp, tau_tree::ZvtxFT0vsPv, tau_tree::VtxITSTPC,
                  tau_tree::TotalFT0AmplitudeA, tau_tree::TotalFT0AmplitudeC, tau_tree::TotalFV0AmplitudeA, tau_tree::EnergyCommonZNA, tau_tree::EnergyCommonZNC,
                  tau_tree::TimeFT0A, tau_tree::TimeFT0C, tau_tree::TimeFV0A, tau_tree::TimeZNA, tau_tree::TimeZNC,
                  tau_tree::TrkPx, tau_tree::TrkPy, tau_tree::TrkPz, tau_tree::TrkSign, tau_tree::TrkDCAxy, tau_tree::TrkDCAz, tau_tree::TrkTimeRes,
                  tau_tree::Trk1ITSclusterSizes, tau_tree::Trk2ITSclusterSizes,
                  tau_tree::TrkTPCsignal, tau_tree::TrkTPCnSigmaEl, tau_tree::TrkTPCnSigmaMu, tau_tree::TrkTPCnSigmaPi, tau_tree::TrkTPCnSigmaKa, tau_tree::TrkTPCnSigmaPr, tau_tree::TrkTPCinnerParam,
                  tau_tree::TrkTOFsignal, tau_tree::TrkTOFnSigmaEl, tau_tree::TrkTOFnSigmaMu, tau_tree::TrkTOFnSigmaPi, tau_tree::TrkTOFnSigmaKa, tau_tree::TrkTOFnSigmaPr, tau_tree::TrkTOFexpMom);

} // namespace o2::aod

struct TauEventTableProducer {
  Produces<o2::aod::TauTwoTracks> tauTwoTracks;

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
    Configurable<float> cutTrueGapSideFV0{"cutTrueGapSideFV0", 180000, "FV0A threshold for SG selector"};
    Configurable<float> cutTrueGapSideFT0A{"cutTrueGapSideFT0A", 150., "FT0A threshold for SG selector"};
    Configurable<float> cutTrueGapSideFT0C{"cutTrueGapSideFT0C", 50., "FT0C threshold for SG selector"};
    Configurable<float> cutTrueGapSideZDC{"cutTrueGapSideZDC", 10000., "ZDC threshold for SG selector. 0 is <1n, 4.2 is <2n, 6.7 is <3n, 9.5 is <4n, 12.5 is <5n"};
    Configurable<float> cutFITtime{"cutFITtime", 40., "Maximum FIT time allowed. Default is 40ns"};
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
  using UDMcParticlesWithUDTracks = soa::Join<aod::UDMcParticles, aod::UDMcParticlesToUDTracks>;
  using UDMcCollisionsWithUDCollisions = soa::Join<aod::UDMcCollisions, aod::UDMcCollisionsToUDCollisions>;
  using UDMcCollisionsWithUDCollision = UDMcCollisionsWithUDCollisions::iterator;

  // init
  void init(InitContext&)
  {
    if (verboseInfo)
      printMediumMessage("INIT METHOD");

    mySetITShitsRule(cutGlobalTrack.cutITShitsRule);

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
  bool isGoodROFtime(C const& coll)
  {

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
        return false; // TPC chi2
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

    int gapSide = collision.gapSide();
    int trueGapSide = sgSelector.trueGap(collision, cutSample.cutTrueGapSideFV0, cutSample.cutTrueGapSideFT0A, cutSample.cutTrueGapSideFT0C, cutSample.cutTrueGapSideZDC);

    if (cutSample.useTrueGap)
      gapSide = trueGapSide;

    if (!isGoodROFtime(collision))
      return;

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

  void processMonteCarlo(UDMcCollisionsWithUDCollision const& mccollision,
                         FullMCSGUDCollisions const&,
                         FullUDTracks const&,
                         UDMcParticlesWithUDTracks const&)
  {
    LOGF(info, "mccollision idx %i", mccollision.globalIndex());
    if (mccollision.has_udcollisions()) {
      auto const& collFromMcColl = mccollision.udcollisions_as<FullMCSGUDCollisions>();
      LOGF(info, "collision size %i ", collFromMcColl.size());
    }
  }
  PROCESS_SWITCH(TauEventTableProducer, processMonteCarlo, "Iterate UD tables with simulated data created by SG-Candidate-Producer.", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TauEventTableProducer>(cfgc)};
}
