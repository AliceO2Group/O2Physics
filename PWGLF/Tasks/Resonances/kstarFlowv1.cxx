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
///
/// \file kstarFlowv1.cxx
/// \brief first order flow harmonic for resonance
/// \author  Prottay Das <prottay.das@cern.ch>, Dukhishyam Mallick <dukhishyam.mallick@cern.ch>
///

#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>

#include <string>
#include <vector>
// #include <TDatabasePDG.h>
#include "PWGLF/DataModel/SPCalibrationTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/GenVector/Boost.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "TF1.h"
#include "TRandom3.h"

#include <array>
#include <cmath>
#include <cstdlib>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
using namespace o2::constants::physics;
struct KstarFlowv1 {

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  // Service<o2::framework::O2DatabasePDG> pdg;

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentrality{"cfgCutCentrality", 80.0f, "Accepted maximum Centrality"};
  // track
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "cut on Charge"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<bool> useGlobalTrack{"useGlobalTrack", true, "use Global track"};
  Configurable<float> nsigmaCutTOF{"nsigmaCutTOF", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<float> nsigmaCutTPC{"nsigmaCutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  // Configurable<bool> removefaketrak{"removefaketrack", true, "Remove fake track from momentum difference"};
  Configurable<float> confFakeKaonCut{"confFakeKaonCut", 0.1, "Cut based on track from momentum difference"};
  Configurable<bool> ispTdepPID{"ispTdepPID", true, "pT dependent PID"};
  Configurable<bool> onlyTOF{"onlyTOF", true, "onlyTOF"};
  Configurable<int> strategyPID{"strategyPID", 2, "PID strategy"};
  Configurable<float> cfgCutTOFBeta{"cfgCutTOFBeta", 0.0, "cut TOF beta"};
  Configurable<float> confMinRot{"confMinRot", 5.0 * o2::constants::math::PI / 6.0, "Minimum of rotation"};
  Configurable<float> confMaxRot{"confMaxRot", 7.0 * o2::constants::math::PI / 6.0, "Maximum of rotation"};
  Configurable<int> nBkgRotations{"nBkgRotations", 9, "Number of rotated copies (background) per each original candidate"};
  Configurable<bool> fillRotation{"fillRotation", true, "fill rotation"};
  Configurable<bool> like{"like", true, "fill rotation"};
  Configurable<int> spNbins{"spNbins", 2000, "Number of bins in sp"};
  Configurable<float> lbinsp{"lbinsp", -1.0, "lower bin value in sp histograms"};
  Configurable<float> hbinsp{"hbinsp", 1.0, "higher bin value in sp histograms"};

  ConfigurableAxis configcentAxis{"configcentAxis", {VARIABLE_WIDTH, 0.0, 10.0, 30.0, 50.0, 80.0}, "Cent V0M"};
  ConfigurableAxis configthnAxisPt{"configthnAxisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 50.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configetaAxis{"configetaAxis", {VARIABLE_WIDTH, -0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8}, "Eta"};
  ConfigurableAxis configIMAxis{"configIMAxis", {VARIABLE_WIDTH, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2}, "IM"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = nabs(aod::cent::centFT0C) < cfgCutCentrality;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter dCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::SPCalibrationTables, aod::Mults>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTOFbeta, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi>>;

  SliceCache cache;
  Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    AxisSpec spAxis = {spNbins, lbinsp, hbinsp, "Sp"};

    histos.add("hpQxtQxpvscent", "hpQxtQxpvscent", HistType::kTHnSparseF, {configcentAxis, spAxis}, true);
    histos.add("hpQytQypvscent", "hpQytQypvscent", HistType::kTHnSparseF, {configcentAxis, spAxis}, true);
    histos.add("hpQxytpvscent", "hpQxytpvscent", HistType::kTHnSparseF, {configcentAxis, spAxis}, true);
    histos.add("hpQxtQypvscent", "hpQxtQypvscent", HistType::kTHnSparseF, {configcentAxis, spAxis}, true);
    histos.add("hpQxpQytvscent", "hpQxpQytvscent", HistType::kTHnSparseF, {configcentAxis, spAxis}, true);

    histos.add("hpv1vscentpteta", "hpv1vscentpteta", HistType::kTHnSparseF, {configIMAxis, configcentAxis, configthnAxisPt, configetaAxis, spAxis}, true);
    histos.add("hpv1vscentptetarot", "hpv1vscentptetarot", HistType::kTHnSparseF, {configIMAxis, configcentAxis, configthnAxisPt, configetaAxis, spAxis}, true);
    histos.add("hpv1vscentptetalike", "hpv1vscentptetalike", HistType::kTHnSparseF, {configIMAxis, configcentAxis, configthnAxisPt, configetaAxis, spAxis}, true);
  }

  double massKa = o2::constants::physics::MassKPlus;
  double massPi = o2::constants::physics::MassPiMinus;

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (useGlobalTrack && !(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    }
    if (!useGlobalTrack && !(candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool strategySelectionPID(const T& candidate, int PID, int strategy)
  {
    if (PID == 0) {
      if (strategy == 0) {
        if (candidate.pt() < 0.5 && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.5 && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOF && candidate.beta() > 0.5) {
          return true;
        }
      } else if (strategy == 1) {
        if (candidate.pt() < 0.5 && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.5 && std::sqrt(candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa() + candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) < nsigmaCutTOF && candidate.beta() > 0.5) {
          return true;
        }
      } else if (strategy == 2) {
        if (candidate.pt() < 0.5 && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.5 && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOF && candidate.beta() > 0.5) {
          return true;
        }
        if (candidate.pt() >= 0.5 && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPC && !candidate.hasTOF()) {
          return true;
        }
      }
    }
    if (PID == 1) {
      if (strategy == 0) {
        if (candidate.pt() < 0.5 && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.5 && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOF && candidate.beta() > 0.5) {
          return true;
        }
      } else if (strategy == 1) {
        if (candidate.pt() < 0.5 && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.5 && std::sqrt(candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi() + candidate.tofNSigmaPi() * candidate.tofNSigmaPi()) < nsigmaCutTOF && candidate.beta() > 0.5) {
          return true;
        }
      } else if (strategy == 2) {
        if (candidate.pt() < 0.5 && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC) {
          return true;
        }
        if (candidate.pt() >= 0.5 && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOF && candidate.beta() > 0.5) {
          return true;
        }
        if (candidate.pt() >= 0.5 && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPC && !candidate.hasTOF()) {
          return true;
        }
      }
    }
    return false;
  }

  double getPhiInRange(double phi)
  {
    double pi = o2::constants::math::PI;
    double twoPi = 2.0 * pi;
    double result = phi;
    // Normalize to [-pi, pi]
    // double result = std::fmod(phi + pi, twoPi);
    // if (result < 0) {
    //  result += twoPi;
    // }
    // result -= pi; // Now result is in [-pi, pi]

    // Convert from [-pi, pi] to [0, pi]
    if (result < 0) {
      result += pi; // Shift negative values to positive
    }

    // If phi > 2π, subtract π instead of normalizing by 2π
    if (phi > twoPi) {
      result -= pi;
    }

    return result; // Ensures range is [0, π]
  }

  template <typename T>
  bool isFakeKaon(T const& track, int /*PID*/)
  {
    const auto pglobal = track.p();
    const auto ptpc = track.tpcInnerParam();
    if (std::abs(pglobal - ptpc) > confFakeKaonCut) {
      return true;
    }
    return false;
  }

  ROOT::Math::PxPyPzMVector KstarMother, daughter1, daughter2, kaonrot, kstarrot;

  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    if (!collision.sel8() || !collision.triggereventsp() || !collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder) || !collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return;
    }
    auto centrality = collision.centFT0C();

    auto qxZDCA = collision.qxZDCA();
    auto qxZDCC = collision.qxZDCC();
    auto qyZDCA = collision.qyZDCA();
    auto qyZDCC = collision.qyZDCC();

    auto proQxtQxp = qxZDCA * qxZDCC;
    auto proQytQyp = qyZDCA * qyZDCC;
    auto proQxytp = proQxtQxp + proQytQyp;
    auto proQxpQyt = qxZDCA * qyZDCC;
    auto proQxtQyp = qxZDCC * qyZDCA;

    histos.fill(HIST("hpQxtQxpvscent"), centrality, proQxtQxp);
    histos.fill(HIST("hpQytQypvscent"), centrality, proQytQyp);
    histos.fill(HIST("hpQxytpvscent"), centrality, proQxytp);
    histos.fill(HIST("hpQxpQytvscent"), centrality, proQxpQyt);
    histos.fill(HIST("hpQxtQypvscent"), centrality, proQxtQyp);

    for (const auto& track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      bool track1kaon = false;
      auto track1ID = track1.globalIndex();
      if (!strategySelectionPID(track1, 0, strategyPID)) {
        continue;
      }
      track1kaon = true;
      for (const auto& track2 : tracks) {
        if (!selectionTrack(track2)) {
          continue;
        }
        bool track2pion = false;
        auto track2ID = track2.globalIndex();
        if (!strategySelectionPID(track2, 1, strategyPID)) {
          continue;
        }
        track2pion = true;
        if (track2ID == track1ID) {
          continue;
        }
        if (!track1kaon || !track2pion) {
          continue;
        }
        daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa);
        daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi);
        KstarMother = daughter1 + daughter2;
        if (std::abs(KstarMother.Rapidity()) > 0.9) {
          continue;
        }

        // // constrain angle to 0 -> [0,0+2pi]
        //       auto phi = RecoDecay::constrainAngle(KstarMother.Phi(), 0,o2::constants::math::TwoPI);

        auto ux = std::cos(getPhiInRange(KstarMother.Phi()));
        auto uy = std::sin(getPhiInRange(KstarMother.Phi()));
        auto v1 = ux * (qxZDCA - qxZDCC) + uy * (qyZDCA - qyZDCC);

        // unlike sign
        if (track1.sign() * track2.sign() < 0) {
          histos.fill(HIST("hpv1vscentpteta"), KstarMother.M(), centrality, KstarMother.Pt(), KstarMother.Rapidity(), v1);
          if (fillRotation) {
            for (int nrotbkg = 0; nrotbkg < nBkgRotations; nrotbkg++) {
              auto anglestart = confMinRot;
              auto angleend = confMaxRot;
              auto anglestep = (angleend - anglestart) / (1.0 * (nBkgRotations - 1));
              auto rotangle = anglestart + nrotbkg * anglestep;
              auto rotkaonPx = track1.px() * std::cos(rotangle) - track1.py() * std::sin(rotangle);
              auto rotkaonPy = track1.px() * std::sin(rotangle) + track1.py() * std::cos(rotangle);
              kaonrot = ROOT::Math::PxPyPzMVector(rotkaonPx, rotkaonPy, track1.pz(), massKa);
              kstarrot = kaonrot + daughter2;

              if (std::abs(kstarrot.Rapidity()) > 0.9) {
                continue;
              }

              auto uxrot = std::cos(getPhiInRange(KstarMother.Phi()));
              auto uyrot = std::sin(getPhiInRange(KstarMother.Phi()));

              auto v1rot = uxrot * (qxZDCA - qxZDCC) + uyrot * (qyZDCA - qyZDCC);

              histos.fill(HIST("hpv1vscentptetarot"), kstarrot.M(), centrality, kstarrot.Pt(), kstarrot.Rapidity(), v1rot);
            }
          }
        }
        // like sign
        if (track1.sign() * track2.sign() > 0) {
          histos.fill(HIST("hpv1vscentptetalike"), kstarrot.M(), centrality, kstarrot.Pt(), kstarrot.Rapidity(), v1);
        }
      }
    }
  }
  PROCESS_SWITCH(KstarFlowv1, processSE, "Process Same event", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<KstarFlowv1>(cfgc)};
}
