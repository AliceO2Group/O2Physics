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

/// \file resonancesGfwFlow.cxx
/// \brief PID flow for resonances using the generic framework
/// \author Preet Bhanjan Pati <preet.bhanjan.pati@cern.ch>

#include <CCDB/BasicCCDBManager.h>
#include <cmath>
#include <vector>
#include <utility>
#include <array>
#include <string>

#include "Math/Vector4D.h"

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/StepTHn.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Multiplicity.h"
#include "PWGLF/DataModel/EPCalibrationTables.h"
#include "CommonConstants/PhysicsConstants.h"

#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/PID.h"

#include <TProfile.h>
#include <TRandom3.h>
#include <TF1.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace std;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct ResonancesGfwFlow {
  Service<ccdb::BasicCCDBManager> ccdb;
  Configurable<int64_t> noLaterThan{"noLaterThan", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgUseNch, bool, false, "Use Nch for flow observables")
  O2_DEFINE_CONFIGURABLE(cfgNbootstrap, int, 10, "Number of subsamples")
  O2_DEFINE_CONFIGURABLE(cfgTpcNsigmaCut, float, 3.0f, "TPC N-sigma cut for pions, kaons, protons")
  O2_DEFINE_CONFIGURABLE(cfgTofNsigmaCut, float, 3.0f, "TOF N-sigma cut for pions, kaons, protons")
  O2_DEFINE_CONFIGURABLE(cfgTofPtCut, float, 0.5f, "Minimum pt to use TOF N-sigma")

  O2_DEFINE_CONFIGURABLE(cfgUsePVContributor, bool, true, "Use PV contributor for tracks")
  O2_DEFINE_CONFIGURABLE(cfgFakePartCut, float, 0.1f, "Maximum difference in measured momentum and TPC inner ring momentum of particle")
  O2_DEFINE_CONFIGURABLE(cfgTpcCluster, int, 70, "Number of TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgUsePointingAngle, bool, false, "Use Pointing angle for resonances")
  O2_DEFINE_CONFIGURABLE(cfgPointingAnglePhi, float, 0.04f, "Minimum Pointing angle for Phi")
  O2_DEFINE_CONFIGURABLE(cfgPointingAngleKo, float, 0.97f, "Minimum Pointing angle for K0")
  O2_DEFINE_CONFIGURABLE(cfgPointingAngleLambda, float, 0.995f, "Minimum Pointing angle for Lambda")
  O2_DEFINE_CONFIGURABLE(cfgUseVoRadius, bool, true, "Use V0 radius for particle identification")
  O2_DEFINE_CONFIGURABLE(cfgVoRadiusMin, float, 0.5f, "Minimum V0 radius in cm")
  O2_DEFINE_CONFIGURABLE(cfgVoRadiusMax, float, 200.0f, "Maximum V0 radius in cm")
  O2_DEFINE_CONFIGURABLE(cfgUseProperLifetime, bool, true, "Use proper lifetime for particle identification")
  O2_DEFINE_CONFIGURABLE(cfgProperLtK0, float, 20.0f, "Minimum lifetime for K0 in cm")
  O2_DEFINE_CONFIGURABLE(cfgProperLtLambda, float, 30.0f, "Minimum lifetime for Lambda in cm")
  O2_DEFINE_CONFIGURABLE(cfgUseDCA, bool, true, "Use dca for daughter tracks")
  O2_DEFINE_CONFIGURABLE(cfgDCAtoPV, float, 0.06f, "Minimum DCA of daughter tracks to primary vertex in cm")
  O2_DEFINE_CONFIGURABLE(cfgDCABetDaug, int, 1, "Maximum DCA between daughter tracks")

  O2_DEFINE_CONFIGURABLE(cfgCutDCAxy, float, 2.0f, "DCAxy range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutDCAz, float, 2.0f, "DCAz range for tracks")
  O2_DEFINE_CONFIGURABLE(additionalEvsel, bool, false, "Additional event selcection")
  O2_DEFINE_CONFIGURABLE(useGlobalTrack, bool, true, "use Global track")
  O2_DEFINE_CONFIGURABLE(cfgITScluster, int, 0, "Number of ITS cluster")
  O2_DEFINE_CONFIGURABLE(cfgCutTOFBeta, float, 0.0, "cut TOF beta")
  O2_DEFINE_CONFIGURABLE(cfgCutOccupancy, int, 3000, "Occupancy cut")
  O2_DEFINE_CONFIGURABLE(ispTdepPID, bool, true, "pT dependent PID")
  O2_DEFINE_CONFIGURABLE(removefaketrack, bool, true, "Remove fake track from momentum difference")
  O2_DEFINE_CONFIGURABLE(confRapidity, float, 0.5, "Rapidity cut")

  // Defining configurable axis
  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for histograms"};
  ConfigurableAxis axisPhi{"axisPhi", {60, 0.0, constants::math::TwoPI}, "phi axis for histograms"};
  ConfigurableAxis axisEta{"axisEta", {40, -1., 1.}, "eta axis for histograms"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.40, 2.60, 2.80, 3.00, 3.50, 4.00, 5.00, 6.00, 8.00, 10.00}, "pt axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90}, "centrality axis for histograms"};
  ConfigurableAxis axisNsigmaTPC{"axisNsigmaTPC", {80, -5, 5}, "nsigmaTPC axis"};
  ConfigurableAxis axisNsigmaTOF{"axisNsigmaTOF", {80, -5, 5}, "nsigmaTOF axis"};
  ConfigurableAxis axisParticles{"axisParticles", {3, 0, 3}, "axis for different hadrons"};
  ConfigurableAxis axisPhiMass{"axisPhiMass", {50000, 0, 5}, "axis for invariant mass distibution for Phi"};
  ConfigurableAxis axisKoMass{"axisKoMass", {50000, 0, 5}, "axis for invariant mass distibution for K0"};
  ConfigurableAxis axisLambdaMass{"axisLambdaMass", {50000, 0, 5}, "axis for invariant mass distibution for Lambda"};
  ConfigurableAxis axisTPCsignal{"axisTPCsignal", {10000, 0, 1000}, "axis for TPC signal"};
  ConfigurableAxis axisTOFsignal{"axisTOFsignal", {10000, 0, 1000}, "axis for TOF signal"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz) && (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtPOIMin) && (aod::track::pt < cfgCutPtPOIMax);

  using AodCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;
  using AodTracksWithoutBayes = soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection, aod::TracksExtra, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFbeta, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;

  SliceCache cache;
  Partition<AodTracksWithoutBayes> posTracks = aod::track::signed1Pt > 0.0f;
  Partition<AodTracksWithoutBayes> negTracks = aod::track::signed1Pt < 0.0f;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  TAxis* fPtAxis;
  TRandom3* fRndm = new TRandom3(0);

  // Event selection cuts - Alex
  TF1* fMultPVCutLow = nullptr;
  TF1* fMultPVCutHigh = nullptr;
  TF1* fMultCutLow = nullptr;
  TF1* fMultCutHigh = nullptr;
  TF1* fMultMultPVCut = nullptr;

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setCreatedNotAfter(noLaterThan.value);

    histos.add("hPhi", "", {HistType::kTH1D, {axisPhi}});
    histos.add("hEta", "", {HistType::kTH1D, {axisEta}});
    histos.add("hVtxZ", "", {HistType::kTH1D, {axisVertex}});
    histos.add("hMult", "", {HistType::kTH1D, {{3000, 0.5, 3000.5}}});
    histos.add("hCent", "", {HistType::kTH1D, {{90, 0, 90}}});

    histos.add("KaplusTPC", "", {HistType::kTH2D, {{axisPt, axisTPCsignal}}});
    histos.add("KaminusTPC", "", {HistType::kTH2D, {{axisPt, axisTPCsignal}}});
    histos.add("KaplusTOF", "", {HistType::kTH2D, {{axisPt, axisTOFsignal}}});
    histos.add("KaminusTOF", "", {HistType::kTH2D, {{axisPt, axisTOFsignal}}});
    histos.add("hPhiMass_sparse", "", {HistType::kTHnSparseD, {{axisPhiMass, axisPt, axisMultiplicity}}});

    if (additionalEvsel) {
      fMultPVCutLow = new TF1("fMultPVCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutLow->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultPVCutHigh = new TF1("fMultPVCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 2.5*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x)", 0, 100);
      fMultPVCutHigh->SetParameters(2834.66, -87.0127, 0.915126, -0.00330136, 332.513, -12.3476, 0.251663, -0.00272819, 1.12242e-05);
      fMultCutLow = new TF1("fMultCutLow", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 2.5*([4]+[5]*x)", 0, 100);
      fMultCutLow->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultCutHigh = new TF1("fMultCutHigh", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 3.*([4]+[5]*x)", 0, 100);
      fMultCutHigh->SetParameters(1893.94, -53.86, 0.502913, -0.0015122, 109.625, -1.19253);
      fMultMultPVCut = new TF1("fMultMultPVCut", "[0]+[1]*x+[2]*x*x", 0, 5000);
      fMultMultPVCut->SetParameters(-0.1, 0.785, -4.7e-05);
    }
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
    angle = std::acos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (cfgUsePointingAngle && angle < cfgPointingAnglePhi) {
      return false;
    }
    return true;
  }

  template <typename TCollision>
  bool eventSelected(TCollision collision, const float& centrality)
  {
    auto multNTracksPV = collision.multNTracksPV();
    if (multNTracksPV < fMultPVCutLow->Eval(centrality))
      return 0;
    if (multNTracksPV > fMultPVCutHigh->Eval(centrality))
      return 0;

    return 1;
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (useGlobalTrack && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTpcCluster)) {
      return false;
    }
    if (!useGlobalTrack && !(candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool selectionPIDpTdependent(const T& candidate)
  {
    if (candidate.pt() < 0.5 && std::abs(candidate.tpcNSigmaKa()) < cfgTpcNsigmaCut) {
      return true;
    }
    if (candidate.pt() >= 0.5 && candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && std::abs(candidate.tpcNSigmaKa()) < cfgTpcNsigmaCut && std::abs(candidate.tofNSigmaKa()) < cfgTofNsigmaCut) {
      return true;
    }
    if (!useGlobalTrack && !candidate.hasTPC()) {
      return true;
    }
    return false;
  }
  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < cfgTpcNsigmaCut) {
      return true;
    }
    if (candidate.hasTOF() && candidate.beta() > cfgCutTOFBeta && std::abs(candidate.tpcNSigmaKa()) < cfgTpcNsigmaCut && std::abs(candidate.tofNSigmaKa()) < cfgTofNsigmaCut) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool isFakeKaon(T const& track)
  {
    const auto pglobal = track.p();
    const auto ptpc = track.tpcInnerParam();
    if (std::abs(pglobal - ptpc) > cfgFakePartCut) {
      return true;
    }
    return false;
  }

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C, o2::aod::evsel::NumTracksInTimeRange>;
  ROOT::Math::PxPyPzMVector phiMom, kaonPlus, kaonminus;
  double massKaPlus = o2::constants::physics::MassKPlus;

  void process(AodCollisions::iterator const& collision, aod::BCsWithTimestamps const&, AodTracksWithoutBayes const& tracks)
  {
    if (!collision.sel8() || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV) || !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return;
    }

    const auto cent = collision.centFT0C();
    int nTot = tracks.size();
    float vtxz = collision.posZ();
    int occupancy = collision.trackOccupancyInTimeRange();

    if (occupancy > cfgCutOccupancy) {
      return;
    }

    if (additionalEvsel && !eventSelected(collision, cent)) {
      return;
    }

    histos.fill(HIST("hVtxZ"), vtxz);
    histos.fill(HIST("hMult"), nTot);
    histos.fill(HIST("hCent"), cent);

    auto posSlicedTracks = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negSlicedTracks = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    for (auto const& track1 : posSlicedTracks) {
      // track selection
      if (!selectionTrack(track1)) {
        continue;
      }
      // PID check
      if (ispTdepPID && !selectionPIDpTdependent(track1)) {
        continue;
      }
      if (!ispTdepPID && !selectionPID(track1)) {
        continue;
      }
      histos.fill(HIST("KaplusTPC"), track1.pt(), track1.tpcNSigmaKa());
      histos.fill(HIST("KaplusTOF"), track1.pt(), track1.tofNSigmaKa());
      auto track1ID = track1.globalIndex();

      for (auto const& track2 : negSlicedTracks) {
        // track selection
        if (!selectionTrack(track2)) {
          continue;
        }
        // PID check
        if (ispTdepPID && !selectionPIDpTdependent(track2)) {
          continue;
        }
        if (!ispTdepPID && !selectionPID(track2)) {
          continue;
        }
        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        if (removefaketrack && isFakeKaon(track1)) {
          continue;
        }
        if (removefaketrack && isFakeKaon(track2)) {
          continue;
        }
        histos.fill(HIST("KaminusTPC"), track2.pt(), track2.tpcNSigmaKa());
        histos.fill(HIST("KaminusTOF"), track2.pt(), track2.tofNSigmaKa());
        kaonPlus = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKaPlus);
        kaonminus = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massKaPlus);
        phiMom = kaonPlus + kaonminus;
        if (std::abs(phiMom.Rapidity()) < confRapidity) {
          histos.fill(HIST("hPhiMass_sparse"), phiMom.M(), phiMom.Pt(), cent);
          histos.fill(HIST("hPhi"), phiMom.Phi());
          histos.fill(HIST("hEta"), phiMom.Eta());
        }
      } //  end of track 2
    } // end of track 1
  } // end of process
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<ResonancesGfwFlow>(cfgc)};
}
