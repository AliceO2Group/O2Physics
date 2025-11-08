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
// Phi meson spin alignment task
// sourav.kundu@cern.ch

#include "PWGMM/Mult/DataModel/Index.h" // for Particles2Tracks table

#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/V0.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TF1.h"
#include "TRandom3.h"
#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace o2::aod::rctsel;
struct phispectrapbpbqa {
  double bz = 0.;

  // Enable access to the CCDB for the offset and correction constants and save them in dedicated variables.
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::framework::O2DatabasePDG> pdg;
  struct : ConfigurableGroup {
    Configurable<std::string> cfgURL{"cfgURL", "http://alice-ccdb.cern.ch", "Address of the CCDB to browse"};
    Configurable<int64_t> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "Latest acceptable timestamp of creation for the object"};
  } cfgCcdbParam;

  struct : ConfigurableGroup {
    Configurable<bool> requireRCTFlagChecker{"requireRCTFlagChecker", true, "Check event quality in run condition table"};
    Configurable<std::string> cfgEvtRCTFlagCheckerLabel{"cfgEvtRCTFlagCheckerLabel", "CBT_hadronPID", "Evt sel: RCT flag checker label"};
    Configurable<bool> cfgEvtRCTFlagCheckerLimitAcceptAsBad{"cfgEvtRCTFlagCheckerLimitAcceptAsBad", true, "Evt sel: RCT flag checker treat Limited Acceptance As Bad"};
  } rctCut;

  TH3D* hTPCCallib;
  TH3D* hTOFCallib;

  Configurable<std::string> ConfPathTPC{"ConfPathTPC", "Users/s/skundu/My/Object/PIDcallib/TPC", "Weight path TPC"};
  Configurable<std::string> ConfPathTOF{"ConfPathTOF", "Users/s/skundu/My/Object/PIDcallib/TOF", "Weight path TOF"};

  // events
  Configurable<bool> applyStrictEvSel{"applyStrictEvSel", true, "Apply strict event selection"};
  Configurable<bool> applyMCsel8{"applyMCsel8", false, "Apply sel8 in MC"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 10, "Number of event mixing"};
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  // ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {8, 0, 80}, "multiplicity percentile for bin"};

  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentrality{"cfgCutCentrality", 80.0f, "Accepted maximum Centrality"};

  // track
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "cut on Charge"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};

  Configurable<int> cfgITScluster{"cfgITScluster", 4, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 80, "Number of TPC cluster"};
  Configurable<int> cfgTPCPIDcluster{"cfgTPCPIDcluster", 80, "Number of TPC PID cluster"};
  Configurable<int> cfgTPCcrossedRows{"cfgTPCcrossedRows", 90, "Number of TPC crossed Rows"};

  Configurable<bool> cfgUpdatePID{"cfgUpdatePID", false, "Update PID callibration"};
  Configurable<bool> applyPID{"applyPID", true, "Apply PID"};
  Configurable<bool> applyPIDCluster{"applyPIDCluster", true, "Apply PID cluster"};
  Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  Configurable<bool> timeFrameMC{"timeFrameMC", false, "time frame cut in MC"};
  Configurable<bool> readOutFrameMC{"readOutFrameMC", true, "ITS read out frame cut in MC"};
  Configurable<bool> ispTdepPID{"ispTdepPID", false, "pT dependent PID"};
  Configurable<float> cfgCutTOFBeta{"cfgCutTOFBeta", 0.5, "cut TOF beta"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 2.0, "Value of the TPC Nsigma cut"};
  Configurable<bool> applyTOF{"applyTOF", true, "Apply TOF"};
  ConfigurableAxis axisOccupancy{"axisOccupancy", {VARIABLE_WIDTH, -1.0, 200.0, 500.0, 1000.0, 2000.0f, 4000.0, 10000.0f, 100000.0f}, "occupancy axis"};
  struct : ConfigurableGroup {
    ConfigurableAxis configThnAxisInvMass{"configThnAxisInvMass", {90, 0.98, 1.07}, "#it{M} (GeV/#it{c}^{2})"};
    ConfigurableAxis configThnAxisPt{"configThnAxisPt", {100, 0.0, 10.}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis configThnAxisCentrality{"configThnAxisCentrality", {8, 0., 80}, "Centrality"};
    ConfigurableAxis configThnAxisSector{"configThnAxisSector", {2, 0.0, 2.0}, "TPC sector"};

  } cnfgaxis;
  Configurable<bool> isMC{"isMC", false, "use MC"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", false, "avoid split track in MC"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = nabs(aod::cent::centFT0C) < cfgCutCentrality;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter PIDcutFilter = nabs(aod::pidtpc::tpcNSigmaKa) < nsigmaCutTPC;
  // Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::TPCMults, aod::CentFT0Cs, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCFullKa, aod::pidTOFFullKa>>;

  using CollisionMCTrueTable = aod::McCollisions;
  using TrackMCTrueTable = aod::McParticles;

  using CollisionMCRecTableCentFT0C = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::CentFT0Cs, aod::EvSels>>;
  using TrackMCRecTable = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::McTrackLabels, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCFullKa, aod::pidTOFFullKa>;
  using FilTrackMCRecTable = soa::Filtered<TrackMCRecTable>;
  Preslice<TrackMCRecTable> perCollision = aod::track::collisionId;

  SliceCache cache;
  Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  RCTFlagsChecker rctChecker;
  // Event selection cuts - Alex
  // TF1* fMultPVCutLow = nullptr;

  void init(o2::framework::InitContext&)
  {
    rctChecker.init(rctCut.cfgEvtRCTFlagCheckerLabel, rctCut.cfgEvtRCTFlagCheckerLimitAcceptAsBad);
    histos.add("hphiSE", "hphiSE", HistType::kTHnSparseF, {cnfgaxis.configThnAxisInvMass, cnfgaxis.configThnAxisPt, cnfgaxis.configThnAxisCentrality, axisOccupancy, cnfgaxis.configThnAxisSector}, true);
    histos.add("hphiME", "hphiME", HistType::kTHnSparseF, {cnfgaxis.configThnAxisInvMass, cnfgaxis.configThnAxisPt, cnfgaxis.configThnAxisCentrality, axisOccupancy, cnfgaxis.configThnAxisSector}, true);
    histos.add("hphiGen", "hphiGen", HistType::kTHnSparseF, {cnfgaxis.configThnAxisInvMass, cnfgaxis.configThnAxisPt, cnfgaxis.configThnAxisCentrality, axisOccupancy}, true);

    histos.add("hNsigmaTPC", "NsigmaKaon TPC", HistType::kTHnSparseF, {{200, -10.0f, 10.0f}, {100, 0.0, 10.0}, axisOccupancy, cnfgaxis.configThnAxisCentrality});
    histos.add("hNsigmaTOF", "NsigmaKaon TOF", HistType::kTHnSparseF, {{200, -10.0f, 10.0f}, {100, 0.0, 10.0}, axisOccupancy, cnfgaxis.configThnAxisCentrality});

    histos.add("hPhiMommentum", "hPhiMommentum", kTH3F, {{36, 0, 6.283}, {200, -10.0, 10.0}, axisOccupancy});

    histos.add("hNsigmaTPCBeforeCut", "NsigmaKaon TPC Before Cut", kTH3F, {{200, -10.0f, 10.0f}, {100, 0.0, 10.0}, axisOccupancy});
    histos.add("hNsigmaTOFBeforeCut", "NsigmaKaon TOF Before Cut", kTH3F, {{200, -10.0f, 10.0f}, {100, 0.0, 10.0}, axisOccupancy});
    histos.add("hNsigmaTPCAfterCut", "NsigmaKaon TPC After Cut", kTH3F, {{200, -10.0f, 10.0f}, {100, 0.0, 10.0}, axisOccupancy});
    histos.add("hNsigmaTOFAfterCut", "NsigmaKaon TOF After Cut", kTH3F, {{200, -10.0f, 10.0f}, {100, 0.0, 10.0}, axisOccupancy});

    histos.add("hNsigmaKaonTOFTPC", "NsigmaKaon TOFTPC distribution", kTH3F, {{200, -10.0f, 10.0f}, {200, -10.0f, 10.0f}, {100, 0.0, 10.0}});

    histos.add("hdcaxy", "DCA xy dist.", kTH2F, {{1000, -1.0f, 1.0f}, {200, -10.0, 10.0}});
    histos.add("hdcaz", "DCA z dist.", kTH2F, {{1000, -1.0f, 1.0f}, {200, -10.0, 10.0}});

    histos.add("hCentrality", "hCentrality", kTH1F, {{8, 0.0f, 80.0f}});
    histos.add("hVtxZ", "hVtxZ", kTH1F, {{8, 0.0f, 80.0f}});
    histos.add("hOccupancy", "hOccupancy", kTH2F, {axisOccupancy, cnfgaxis.configThnAxisCentrality});
    histos.add("hMC", "hMC", kTH1F, {{20, 0.0f, 20.0f}});
    histos.add("h1PhiRecsplit", "h1PhiRecsplit", kTH1F, {{100, 0.0f, 10.0f}});
    histos.add("hData", "hData", kTH1F, {{20, 0.0f, 20.0f}});
    ccdb->setURL(cfgCcdbParam.cfgURL);
    ccdbApi.init("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    if (cfgUpdatePID) {
      hTPCCallib = ccdb->getForTimeStamp<TH3D>(ConfPathTPC.value, cfgCcdbParam.nolaterthan.value);
      hTOFCallib = ccdb->getForTimeStamp<TH3D>(ConfPathTOF.value, cfgCcdbParam.nolaterthan.value);
    }
  }

  double massKa = o2::constants::physics::MassKPlus;

  int getMagneticField(uint64_t timestamp)
  {
    // Get the magnetic field
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("/GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (!(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster && candidate.tpcNClsCrossedRows() > cfgTPCcrossedRows)) {
      return false;
    }
    if (applyPIDCluster && candidate.tpcNClsPID() < cfgTPCPIDcluster) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPIDpTdependent(const T& candidate, double nsigmaTPC, double nsigmaTOF)
  {
    if (candidate.p() < 0.7 && TMath::Abs(nsigmaTPC) < nsigmaCutTPC) {
      return true;
    }
    if (candidate.p() > 0.7 && candidate.hasTOF() && TMath::Abs(nsigmaTPC) < nsigmaCutTPC) {
      if (candidate.p() > 0.7 && candidate.p() < 1.6 && candidate.beta() > cfgCutTOFBeta && nsigmaTOF > -5.0 && nsigmaTOF < 10.0) {
        return true;
      }
      if (candidate.p() >= 1.6 && candidate.p() < 2.0 && candidate.beta() > cfgCutTOFBeta && nsigmaTOF > -3.0 && nsigmaTOF < 10.0) {
        return true;
      }
      if (candidate.p() >= 2.0 && candidate.p() < 2.5 && candidate.beta() > cfgCutTOFBeta && nsigmaTOF > -3.0 && nsigmaTOF < 6.0) {
        return true;
      }
      if (candidate.p() >= 2.5 && candidate.p() < 4.0 && candidate.beta() > cfgCutTOFBeta && nsigmaTOF > -2.5 && nsigmaTOF < 4.0) {
        return true;
      }
      if (candidate.p() >= 4.0 && candidate.p() < 5.0 && candidate.beta() > cfgCutTOFBeta && nsigmaTOF > -4.0 && nsigmaTOF < 3.0) {
        return true;
      }
      if (candidate.p() >= 5.0 && candidate.p() < 6.0 && candidate.beta() > cfgCutTOFBeta && nsigmaTOF > -4.0 && nsigmaTOF < 2.5) {
        return true;
      }
      if (candidate.p() >= 6.0 && candidate.beta() > cfgCutTOFBeta && nsigmaTOF > -3.0 && nsigmaTOF < 3.0) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selectionPID(const T& candidate, double nsigmaTPC, double nsigmaTOF)
  {
    if (applyTOF) {
      if (!candidate.hasTOF() && TMath::Abs(nsigmaTPC) < nsigmaCutTPC) {
        return true;
      }
      if (candidate.p() > 0.5 && candidate.hasTOF() && TMath::Abs(nsigmaTPC) < nsigmaCutTPC) {
        if (candidate.p() > 0.5 && candidate.p() < 1.6 && candidate.beta() > cfgCutTOFBeta && nsigmaTOF > -5.0 && nsigmaTOF < 10.0) {
          return true;
        }
        if (candidate.p() >= 1.6 && candidate.p() < 2.0 && candidate.beta() > cfgCutTOFBeta && nsigmaTOF > -3.0 && nsigmaTOF < 10.0) {
          return true;
        }
        if (candidate.p() >= 2.0 && candidate.p() < 2.5 && candidate.beta() > cfgCutTOFBeta && nsigmaTOF > -3.0 && nsigmaTOF < 6.0) {
          return true;
        }
        if (candidate.p() >= 2.5 && candidate.p() < 4.0 && candidate.beta() > cfgCutTOFBeta && nsigmaTOF > -2.5 && nsigmaTOF < 4.0) {
          return true;
        }
        if (candidate.p() >= 4.0 && candidate.p() < 5.0 && candidate.beta() > cfgCutTOFBeta && nsigmaTOF > -4.0 && nsigmaTOF < 3.0) {
          return true;
        }
        if (candidate.p() >= 5.0 && candidate.p() < 6.0 && candidate.beta() > cfgCutTOFBeta && nsigmaTOF > -4.0 && nsigmaTOF < 2.5) {
          return true;
        }
        if (candidate.p() >= 6.0 && candidate.beta() > cfgCutTOFBeta && nsigmaTOF > -3.0 && nsigmaTOF < 3.0) {
          return true;
        }
      }
    } else if (TMath::Abs(nsigmaTPC) < nsigmaCutTPC) {
      return true;
    }
    return false;
  }

  // deep angle cut on pair to remove photon conversion
  template <typename T1, typename T2>
  bool selectionPair(const T1& candidate1, const T2& candidate2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = candidate1.pt();
    pt2 = candidate2.pt();
    pz1 = candidate1.pz();
    pz2 = candidate2.pz();
    p1 = candidate1.p();
    p2 = candidate2.p();
    angle = TMath::ACos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (isDeepAngle && angle < cfgDeepAngle) {
      return false;
    }
    return true;
  }

  // Keep a track only if its azimuth is NOT inside the periodic boundary window.
  // phi           : track azimuth (radians) at your chosen reference radius (e.g. Rout ~ 247 cm)
  // badHalfWidth  : half-width of the bad strip near each boundary (rad). Your example: 0.15
  // nSectors      : 18 for ALICE TPC
  bool keepTrackNoTPCBoundary(float phi,
                              float badWidth = 0.15f,
                              int nSectors = 18)
  {
    constexpr float TwoPi = 6.283185307179586f;
    float ph = std::fmod(phi, TwoPi);
    if (ph < 0.f)
      ph += TwoPi; // wrap phi into [0, 2π)

    const float secW = TwoPi / static_cast<float>(nSectors); // ≈ 0.349 rad

    float rel = std::fmod(ph, secW); // position inside sector [0, secW)

    return rel > badWidth; // keep if NOT inside boundary strip [0, badWidth]
  }

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, o2::aod::evsel::NumTracksInTimeRange>;
  ROOT::Math::PxPyPzMVector PhiMesonMother, KaonPlus, KaonMinus, fourVecDauCM;
  int currentRunNumber = -999;
  int lastRunNumber = -999;
  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const&, aod::BCsWithTimestamps const&)
  {
    histos.fill(HIST("hData"), 1);
    if (!collision.sel8() || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return;
    }
    histos.fill(HIST("hData"), 2);
    if (rctCut.requireRCTFlagChecker) {
      if (!rctChecker(collision)) {
        return;
      }
    }
    histos.fill(HIST("hData"), 3);
    if (applyStrictEvSel && (!collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll) || !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow))) {
      return;
    }
    histos.fill(HIST("hData"), 4);
    auto centrality = collision.centFT0C();
    int occupancy = collision.trackOccupancyInTimeRange();
    histos.fill(HIST("hCentrality"), centrality);
    histos.fill(HIST("hVtxZ"), collision.posZ());

    /*auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    currentRunNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    if (currentRunNumber != lastRunNumber) {
      bz = getMagneticField(bc.timestamp());
    }
    lastRunNumber = currentRunNumber;
    */

    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    int Npostrack = 0;
    histos.fill(HIST("hOccupancy"), occupancy, centrality);
    for (auto track1 : posThisColl) {
      if (!selectionTrack(track1)) {
        continue;
      }
      histos.fill(HIST("hdcaxy"), track1.dcaXY(), track1.p() / track1.sign());
      histos.fill(HIST("hdcaz"), track1.dcaZ(), track1.p() / track1.sign());
      double nSigmaTPC = track1.tpcNSigmaKa();
      double nSigmaTOF = track1.tofNSigmaKa();
      if (!track1.hasTOF()) {
        nSigmaTOF = -9999.99;
      }

      if (cfgUpdatePID) {
        // update PID
        nSigmaTPC = (nSigmaTPC - hTPCCallib->GetBinContent(hTPCCallib->FindBin(track1.p(), centrality, occupancy))) / hTPCCallib->GetBinError(hTPCCallib->FindBin(track1.p(), centrality, occupancy));
        if (track1.hasTOF()) {
          nSigmaTOF = (nSigmaTOF - hTOFCallib->GetBinContent(hTOFCallib->FindBin(track1.p(), centrality, occupancy))) / hTOFCallib->GetBinError(hTOFCallib->FindBin(track1.p(), centrality, occupancy));
        }
      }
      if (track1.p() < 0.6) {
        histos.fill(HIST("hNsigmaTPC"), nSigmaTPC, track1.p(), occupancy, centrality);
      } else if (track1.p() > 0.6 && track1.hasTOF() && std::abs(nSigmaTOF) < 2.5) {
        histos.fill(HIST("hNsigmaTPC"), nSigmaTPC, track1.p(), occupancy, centrality);
      }
      if (track1.hasTOF()) {
        histos.fill(HIST("hNsigmaTOF"), nSigmaTOF, track1.p(), occupancy, centrality);
      }
      histos.fill(HIST("hNsigmaTPCBeforeCut"), nSigmaTPC, track1.p(), occupancy);
      histos.fill(HIST("hNsigmaTOFBeforeCut"), nSigmaTOF, track1.p(), occupancy);
      if (applyPID) {
        if (ispTdepPID && !selectionPIDpTdependent(track1, nSigmaTPC, nSigmaTOF)) {
          continue;
        }
        if (!ispTdepPID && !selectionPID(track1, nSigmaTPC, nSigmaTOF)) {
          continue;
        }
      }
      histos.fill(HIST("hNsigmaTPCAfterCut"), nSigmaTPC, track1.p(), occupancy);
      histos.fill(HIST("hNsigmaTOFAfterCut"), nSigmaTOF, track1.p(), occupancy);
      if (track1.hasTOF()) {
        histos.fill(HIST("hNsigmaKaonTOFTPC"), nSigmaTOF, nSigmaTPC, track1.p());
      }
      Npostrack = Npostrack + 1;

      // 1) φ at a chosen radius (e.g., outer pad rows ~247 cm)
      histos.fill(HIST("hPhiMommentum"), track1.phi(), track1.p(), occupancy);

      for (auto track2 : negThisColl) {
        if (track1.sign() * track2.sign() > 0.0) {
          continue;
        }
        if (!selectionTrack(track2)) {
          continue;
        }
        if (Npostrack == 1) {
          histos.fill(HIST("hdcaxy"), track2.dcaXY(), track2.p() / track2.sign());
          histos.fill(HIST("hdcaz"), track2.dcaZ(), track2.p() / track2.sign());
        }
        double nSigmaTPC2 = track2.tpcNSigmaKa();
        double nSigmaTOF2 = track2.tofNSigmaKa();
        if (!track2.hasTOF()) {
          nSigmaTOF2 = -9999.9;
        }
        if (cfgUpdatePID) {
          // update PID
          nSigmaTPC2 = (nSigmaTPC2 - hTPCCallib->GetBinContent(hTPCCallib->FindBin(track2.p(), centrality, occupancy))) / hTPCCallib->GetBinError(hTPCCallib->FindBin(track2.p(), centrality, occupancy));
          if (track2.hasTOF()) {
            nSigmaTOF2 = (nSigmaTOF2 - hTOFCallib->GetBinContent(hTOFCallib->FindBin(track2.p(), centrality, occupancy))) / hTOFCallib->GetBinError(hTOFCallib->FindBin(track2.p(), centrality, occupancy));
          }
        }
        if (Npostrack == 1) {
          if (track2.p() < 0.6) {
            histos.fill(HIST("hNsigmaTPC"), nSigmaTPC2, track2.p(), occupancy, centrality);
          } else if (track2.p() > 0.6 && track2.hasTOF() && std::abs(nSigmaTOF2) < 2.5) {
            histos.fill(HIST("hNsigmaTPC"), nSigmaTPC2, track2.p(), occupancy, centrality);
          }
          if (track2.hasTOF()) {
            histos.fill(HIST("hNsigmaTOF"), nSigmaTOF2, track2.p(), occupancy, centrality);
          }
          histos.fill(HIST("hNsigmaTPCBeforeCut"), nSigmaTPC2, track2.p(), occupancy);
          histos.fill(HIST("hNsigmaTOFBeforeCut"), nSigmaTOF2, track2.p(), occupancy);
        }
        if (applyPID) {
          if (ispTdepPID && !selectionPIDpTdependent(track2, nSigmaTPC2, nSigmaTOF2)) {
            continue;
          }
          if (!ispTdepPID && !selectionPID(track2, nSigmaTPC2, nSigmaTOF2)) {
            continue;
          }
        }
        if (Npostrack == 1) {
          histos.fill(HIST("hNsigmaTPCAfterCut"), nSigmaTPC2, track2.p(), occupancy);
          histos.fill(HIST("hNsigmaTOFAfterCut"), nSigmaTOF2, track2.p(), occupancy);
          if (track2.hasTOF()) {
            histos.fill(HIST("hNsigmaKaonTOFTPC"), nSigmaTOF2, nSigmaTPC2, track2.p());
          }
        }
        if (Npostrack == 1) {
          histos.fill(HIST("hPhiMommentum"), track2.phi(), track2.p(), occupancy);
        }
        auto track1ID = track1.globalIndex();
        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        bool crossed1 = keepTrackNoTPCBoundary(track1.phi(), 0.15, 18);
        bool crossed2 = keepTrackNoTPCBoundary(track2.phi(), 0.15, 18);
        float passsector = -999.0;
        if (crossed1 && crossed2) {
          passsector = 1.5;
        } else {
          passsector = 0.5;
        }

        KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        PhiMesonMother = KaonPlus + KaonMinus;

        if (TMath::Abs(PhiMesonMother.Rapidity()) < 0.5) {
          histos.fill(HIST("hphiSE"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, occupancy, passsector);
        }
      }
    }
  }
  PROCESS_SWITCH(phispectrapbpbqa, processSameEvent, "Process Same event", true);

  void processMixedEventOpti(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    auto tracksTuple = std::make_tuple(tracks);
    BinningTypeVertexContributor binningOnPositions{{axisVertex, cnfgaxis.configThnAxisCentrality, axisOccupancy}, true};
    SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor> pair{binningOnPositions, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};
    for (auto& [collision1, tracks1, collision2, tracks2] : pair) {
      if (!collision1.sel8() || !collision1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision1.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision1.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision1.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        continue;
      }
      if (!collision2.sel8() || !collision2.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision2.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision2.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision2.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision2.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        continue;
      }
      if (rctCut.requireRCTFlagChecker) {
        if (!rctChecker(collision1) || !rctChecker(collision2)) {
          continue;
        }
      }
      if (applyStrictEvSel && (!collision1.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll) || !collision1.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow))) {
        return;
      }
      if (applyStrictEvSel && (!collision2.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll) || !collision2.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow))) {
        return;
      }
      int occupancy = collision1.trackOccupancyInTimeRange();
      auto centrality = collision1.centFT0C();
      for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (track1.sign() * track2.sign() > 0) {
          continue;
        }
        if (!selectionTrack(track1)) {
          continue;
        }
        if (!selectionTrack(track2)) {
          continue;
        }
        // PID track 1
        double nSigmaTPC = track1.tpcNSigmaKa();
        double nSigmaTOF = track1.tofNSigmaKa();
        if (!track1.hasTOF()) {
          nSigmaTOF = -9999.99;
        }
        if (cfgUpdatePID) {
          nSigmaTPC = (nSigmaTPC - hTPCCallib->GetBinContent(hTPCCallib->FindBin(track1.p(), centrality, occupancy))) / hTPCCallib->GetBinError(hTPCCallib->FindBin(track1.p(), centrality, occupancy));
          if (track1.hasTOF()) {
            nSigmaTOF = (nSigmaTOF - hTOFCallib->GetBinContent(hTOFCallib->FindBin(track1.p(), centrality, occupancy))) / hTOFCallib->GetBinError(hTOFCallib->FindBin(track1.p(), centrality, occupancy));
          }
        }
        // PID track 2
        double nSigmaTPC2 = track2.tpcNSigmaKa();
        double nSigmaTOF2 = track2.tofNSigmaKa();
        if (!track2.hasTOF()) {
          nSigmaTOF2 = -9999.99;
        }
        if (cfgUpdatePID) {
          nSigmaTPC2 = (nSigmaTPC2 - hTPCCallib->GetBinContent(hTPCCallib->FindBin(track2.p(), centrality, occupancy))) / hTPCCallib->GetBinError(hTPCCallib->FindBin(track2.p(), centrality, occupancy));
          if (track2.hasTOF()) {
            nSigmaTOF2 = (nSigmaTOF2 - hTOFCallib->GetBinContent(hTOFCallib->FindBin(track2.p(), centrality, occupancy))) / hTOFCallib->GetBinError(hTOFCallib->FindBin(track2.p(), centrality, occupancy));
          }
        }
        if (applyPID) {
          if (ispTdepPID && !selectionPIDpTdependent(track1, nSigmaTPC, nSigmaTOF)) {
            continue;
          }
          if (!ispTdepPID && !selectionPID(track1, nSigmaTPC, nSigmaTOF)) {
            continue;
          }

          if (ispTdepPID && !selectionPIDpTdependent(track2, nSigmaTPC2, nSigmaTOF2)) {
            continue;
          }
          if (!ispTdepPID && !selectionPID(track2, nSigmaTPC2, nSigmaTOF2)) {
            continue;
          }
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        bool crossed1 = keepTrackNoTPCBoundary(track1.phi(), 0.15, 18);
        bool crossed2 = keepTrackNoTPCBoundary(track2.phi(), 0.15, 18);
        float passsector = -999.0;
        if (crossed1 && crossed2) {
          passsector = 1.5;
        } else {
          passsector = 0.5;
        }
        KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
        PhiMesonMother = KaonPlus + KaonMinus;
        if (TMath::Abs(PhiMesonMother.Rapidity()) < 0.5) {
          histos.fill(HIST("hphiME"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, occupancy, passsector);
        }
      }
    }
  }
  PROCESS_SWITCH(phispectrapbpbqa, processMixedEventOpti, "Process Mixed event new", true);

  void processMC(CollisionMCTrueTable::iterator const&, CollisionMCRecTableCentFT0C const& RecCollisions, TrackMCTrueTable const& GenParticles, FilTrackMCRecTable const& RecTracks)
  {
    histos.fill(HIST("hMC"), 0);
    if (RecCollisions.size() == 0) {
      histos.fill(HIST("hMC"), 1);
      return;
    }
    if (RecCollisions.size() > 1) {
      histos.fill(HIST("hMC"), 2);
      return;
    }
    for (auto& RecCollision : RecCollisions) {
      if (applyMCsel8 && !RecCollision.sel8()) {
        histos.fill(HIST("hMC"), 3);
        continue;
      }
      if (!applyMCsel8 && !RecCollision.selection_bit(aod::evsel::kIsTriggerTVX)) {
        histos.fill(HIST("hMC"), 3);
        continue;
      }
      if (!RecCollision.selection_bit(aod::evsel::kNoSameBunchPileup) || !RecCollision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !RecCollision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
        histos.fill(HIST("hMC"), 4);
        continue;
      }
      if (TMath::Abs(RecCollision.posZ()) > cfgCutVertex) {
        histos.fill(HIST("hMC"), 6);
        continue;
      }

      if (timeFrameMC && !RecCollision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        histos.fill(HIST("hMC"), 7);
        continue;
      }
      if (readOutFrameMC && !RecCollision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
        histos.fill(HIST("hMC"), 8);
        continue;
      }
      histos.fill(HIST("hMC"), 9);
      if (rctCut.requireRCTFlagChecker) {
        if (!rctChecker(RecCollision)) {
          continue;
        }
      }
      histos.fill(HIST("hMC"), 10);
      if (applyStrictEvSel && (!RecCollision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll) || !RecCollision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow))) {
        return;
      }
      histos.fill(HIST("hMC"), 11);
      auto centrality = RecCollision.centFT0C();
      int occupancy = RecCollision.trackOccupancyInTimeRange();
      histos.fill(HIST("hOccupancy"), occupancy, centrality);

      auto oldindex = -999;
      auto Rectrackspart = RecTracks.sliceBy(perCollision, RecCollision.globalIndex());
      // loop over reconstructed particle
      int ntrack1 = 0;
      for (auto track1 : Rectrackspart) {
        if (!selectionTrack(track1)) {
          continue;
        }
        if (!track1.has_mcParticle()) {
          continue;
        }
        auto track1ID = track1.index();
        // PID track 1
        double nSigmaTPC = track1.tpcNSigmaKa();
        double nSigmaTOF = track1.tofNSigmaKa();
        if (!track1.hasTOF()) {
          nSigmaTOF = -9999.99;
        }
        if (cfgUpdatePID) {
          nSigmaTPC = (nSigmaTPC - hTPCCallib->GetBinContent(hTPCCallib->FindBin(track1.p(), centrality, occupancy))) / hTPCCallib->GetBinError(hTPCCallib->FindBin(track1.p(), centrality, occupancy));
          if (track1.hasTOF()) {
            nSigmaTOF = (nSigmaTOF - hTOFCallib->GetBinContent(hTOFCallib->FindBin(track1.p(), centrality, occupancy))) / hTOFCallib->GetBinError(hTOFCallib->FindBin(track1.p(), centrality, occupancy));
          }
        }
        histos.fill(HIST("hNsigmaTPCBeforeCut"), nSigmaTPC, track1.p(), occupancy);
        histos.fill(HIST("hNsigmaTOFBeforeCut"), nSigmaTOF, track1.p(), occupancy);

        if (applyPID) {
          if (ispTdepPID && !selectionPIDpTdependent(track1, nSigmaTPC, nSigmaTOF)) {
            continue;
          }
          if (!ispTdepPID && !selectionPID(track1, nSigmaTPC, nSigmaTOF)) {
            continue;
          }
          histos.fill(HIST("hNsigmaTPCAfterCut"), nSigmaTPC, track1.p(), occupancy);
          histos.fill(HIST("hNsigmaTOFAfterCut"), nSigmaTOF, track1.p(), occupancy);
        }
        ntrack1 = ntrack1 + 1;
        for (auto track2 : Rectrackspart) {
          auto track2ID = track2.index();
          if (track2ID <= track1ID) {
            continue;
          }
          if (!selectionTrack(track2)) {
            continue;
          }
          if (!track2.has_mcParticle()) {
            continue;
          }
          if (!selectionPair(track1, track2)) {
            continue;
          }
          if (track1.sign() * track2.sign() > 0) {
            continue;
          }

          // PID track 2
          double nSigmaTPC2 = track2.tpcNSigmaKa();
          double nSigmaTOF2 = track2.tofNSigmaKa();
          if (!track2.hasTOF()) {
            nSigmaTOF2 = -9999.99;
          }
          if (cfgUpdatePID) {
            nSigmaTPC2 = (nSigmaTPC2 - hTPCCallib->GetBinContent(hTPCCallib->FindBin(track2.p(), centrality, occupancy))) / hTPCCallib->GetBinError(hTPCCallib->FindBin(track2.p(), centrality, occupancy));
            if (track2.hasTOF()) {
              nSigmaTOF2 = (nSigmaTOF2 - hTOFCallib->GetBinContent(hTOFCallib->FindBin(track2.p(), centrality, occupancy))) / hTOFCallib->GetBinError(hTOFCallib->FindBin(track2.p(), centrality, occupancy));
            }
          }
          if (ntrack1 == 1) {
            histos.fill(HIST("hNsigmaTPCBeforeCut"), nSigmaTPC2, track2.p(), occupancy);
            histos.fill(HIST("hNsigmaTOFBeforeCut"), nSigmaTOF2, track2.p(), occupancy);
          }
          if (applyPID) {
            if (ispTdepPID && !selectionPIDpTdependent(track2, nSigmaTPC2, nSigmaTOF2)) {
              continue;
            }
            if (!ispTdepPID && !selectionPID(track2, nSigmaTPC2, nSigmaTOF2)) {
              continue;
            }
          }
          if (ntrack1 == 1) {
            histos.fill(HIST("hNsigmaTPCAfterCut"), nSigmaTPC2, track2.p(), occupancy);
            histos.fill(HIST("hNsigmaTOFAfterCut"), nSigmaTOF2, track2.p(), occupancy);
          }

          const auto mctrack1 = track1.mcParticle();
          const auto mctrack2 = track2.mcParticle();
          int track1PDG = TMath::Abs(mctrack1.pdgCode());
          int track2PDG = TMath::Abs(mctrack2.pdgCode());
          if (!mctrack1.isPhysicalPrimary()) {
            continue;
          }
          if (!mctrack2.isPhysicalPrimary()) {
            continue;
          }
          if (!(track1PDG == 321 && track2PDG == 321)) {
            continue;
          }
          for (auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
            for (auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
              if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
                continue;
              }
              if (mothertrack1 != mothertrack2) {
                continue;
              }
              if (TMath::Abs(mothertrack1.pdgCode()) != 333) {
                continue;
              }
              if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
                histos.fill(HIST("h1PhiRecsplit"), mothertrack1.pt());
                continue;
              }
              oldindex = mothertrack1.globalIndex();
              bool crossed1 = keepTrackNoTPCBoundary(track1.phi(), 0.15, 18);
              bool crossed2 = keepTrackNoTPCBoundary(track2.phi(), 0.15, 18);
              float passsector = -999.0;
              if (crossed1 && crossed2) {
                passsector = 1.5;
              } else {
                passsector = 0.5;
              }
              KaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
              KaonMinus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKa);
              PhiMesonMother = KaonPlus + KaonMinus;
              if (TMath::Abs(PhiMesonMother.Rapidity()) < 0.5) {
                histos.fill(HIST("hphiSE"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, occupancy, passsector);
              }
            }
          }
        }
      }
      // loop over generated particle
      for (auto& mcParticle : GenParticles) {
        if (TMath::Abs(mcParticle.y()) > 0.5) {
          continue;
        }
        if (mcParticle.pdgCode() != 333) {
          continue;
        }
        auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
        if (kDaughters.size() != 2) {
          continue;
        }
        auto daughtp = false;
        auto daughtm = false;
        for (auto kCurrentDaughter : kDaughters) {
          if (!kCurrentDaughter.isPhysicalPrimary()) {
            continue;
          }
          if (kCurrentDaughter.pdgCode() == +321) {
            if (kCurrentDaughter.pt() > cfgCutPT && TMath::Abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtp = true;
            }
            KaonPlus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          } else if (kCurrentDaughter.pdgCode() == -321) {
            if (kCurrentDaughter.pt() > cfgCutPT && TMath::Abs(kCurrentDaughter.eta()) < cfgCutEta) {
              daughtm = true;
            }
            KaonMinus = ROOT::Math::PxPyPzMVector(kCurrentDaughter.px(), kCurrentDaughter.py(), kCurrentDaughter.pz(), massKa);
          }
        }
        if (daughtp && daughtm) {
          PhiMesonMother = KaonPlus + KaonMinus;
          if (TMath::Abs(PhiMesonMother.Rapidity()) < 0.5) {
            histos.fill(HIST("hphiGen"), PhiMesonMother.M(), PhiMesonMother.Pt(), centrality, occupancy);
          }
        }
      }
    } // rec collision loop

  } // process MC
  PROCESS_SWITCH(phispectrapbpbqa, processMC, "Process MC", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phispectrapbpbqa>(cfgc, TaskName{"phispectrapbpbqa"})};
}
