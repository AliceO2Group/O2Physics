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

/// \file mlModelGenHists
/// \brief Generate momentum TH1Fs for mlAccepted mcParticles by ML model and for MC mcParticles.
///
/// \author Michał Olędzki <mioledzk@cern.ch>
/// \author Marek Mytkowski <mmytkows@cern.ch>

#include <Framework/AnalysisDataModel.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CCDB/CcdbApi.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Tools/PIDML/pidOnnxModel.h"

#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define EPS 0.00001f
#define ETA_CUT 0.8f
#define TOF_ACCEPTANCE_THRESHOLD (-999.0f)
#define TRD_ACCEPTANCE_THRESHOLD (0.0f)

struct MlModelGenHists {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  PidONNXModel pidModel; // One instance per model, e.g., one per each pid to predict
  Configurable<uint32_t> cfgDetector{"detector", kTPCTOFTRD, "What detectors to use: 0: TPC only, 1: TPC + TOF, 2: TPC + TOF + TRD"};
  Configurable<int> cfgPid{"pid", 211, "PID to predict"};
  Configurable<double> cfgNSigmaCut{"n-sigma-cut", 3.0f, "TPC and TOF PID nSigma cut"};
  Configurable<double> cfgTofPtCut{"tof-pt-cut", 0.5f, "From what pT TOF is used"};
  Configurable<double> cfgCertainty{"certainty", 0.5, "Min certainty of the model to accept given mcPart to be of given kind"};

  Configurable<std::string> cfgPathCCDB{"ccdb-path", "Users/m/mkabus/PIDML", "base path to the CCDB directory with ONNX models"};
  Configurable<std::string> cfgCCDBURL{"ccdb-url", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  Configurable<bool> cfgUseCCDB{"useCCDB", true, "Whether to autofetch ML model from CCDB. If false, local file will be used."};
  Configurable<std::string> cfgPathLocal{"local-path", "/home/mkabus/PIDML", "base path to the local directory with ONNX models"};

  Configurable<bool> cfgUseFixedTimestamp{"use-fixed-timestamp", false, "Whether to use fixed timestamp from configurable instead of timestamp calculated from the data"};
  Configurable<uint64_t> cfgTimestamp{"timestamp", 1524176895000, "Hardcoded timestamp for tests"};


  o2::ccdb::CcdbApi ccdbApi;
  int currentRunNumber = -1;

  Filter trackFilter = requireGlobalTrackInFilter() &&
    (nabs(aod::pidtofsignal::tofSignal - TOF_ACCEPTANCE_THRESHOLD) > EPS) &&
    (nabs(aod::track::trdSignal - TRD_ACCEPTANCE_THRESHOLD) > EPS);

  using BigTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal, aod::McTrackLabels,
    aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>>;
  
  bool IsNSigmaAccept(const BigTracks::iterator &track)
  {
    double nsigmaTPC, nsigmaTOF;

    switch(cfgPid) {
      case 211: // pion
        nsigmaTOF = track.tofNSigmaPi();
        nsigmaTPC = track.tpcNSigmaPi();
        break;
      case 321: // kaon
        nsigmaTOF = track.tofNSigmaKa();
        nsigmaTPC = track.tpcNSigmaKa();
        break;
      case 2212: // proton
        nsigmaTOF = track.tofNSigmaPr();
        nsigmaTPC = track.tpcNSigmaPr();
        break;
      default:
        return false;
    }

    histos.fill(HIST("hPtTOFNSigma"), track.pt(), nsigmaTOF);
    histos.fill(HIST("hPtTPCNSigma"), track.pt(), nsigmaTPC);

    if (cfgDetector == 0) {
      if (TMath::Abs(nsigmaTPC) < cfgNSigmaCut)
        return true;
      return false;
    } else {
      if (TMath::Hypot(nsigmaTOF, nsigmaTPC) < cfgNSigmaCut)
        return true;
      return false;
    }
  }

  void init(InitContext const&)
  {
    if (cfgUseCCDB) {
      ccdbApi.init(cfgCCDBURL);
    } else {
      pidModel = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, -1, cfgPid.value, static_cast<PidMLDetector>(cfgDetector.value), cfgCertainty.value);
    }

    const AxisSpec axisPt{100, 0, 3.1, "pt"};
    const AxisSpec axisBeta{100, 0, 1.0, "beta"};
    const AxisSpec axisTPCSignal{100, 0, 120.0, "dEdx"};
    const AxisSpec axisNSigma{100, -5.0, 5.0, "n-sigma"};

    // Monte Carlo
    histos.add("hPtMCPositive", "hPtMCPositive", kTH1F, {axisPt});
    histos.add("hPtMCTracked", "hPtMCTracked", kTH1F, {axisPt});

    // Machine learning PID
    histos.add("hPtMLPositive", "hPtMLPositive", kTH1F, {axisPt});
    histos.add("hPtMLTruePositive", "hPtMLTruePositive", kTH1F, {axisPt});

    // NSigma PID
    histos.add("hPtNSigmaPositive", "hPtNSigmaPositive", kTH1F, {axisPt});
    histos.add("hPtNSigmaTruePositive", "hPtNSigmaTruePositive", kTH1F, {axisPt});

    // Context detectors' data
    histos.add("hPtTOFNSigma", "hPtTOFNSigma", kTH2F, {axisPt, axisNSigma});
    histos.add("hPtTOFBeta", "hPtTOFBeta", kTH2F, {axisPt, axisBeta});
    histos.add("hPtTPCNSigma", "hPtTPCNSigma", kTH2F, {axisPt, axisNSigma});
    histos.add("hPtTPCSignal", "hPtTPCSignal", kTH2F, {axisPt, axisTPCSignal});
  }

  void process(aod::Collisions const& collisions, BigTracks const& tracks, aod::BCsWithTimestamps const&, aod::McParticles const& mcParticles)
  {
    auto bc = collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
    if (cfgUseCCDB && bc.runNumber() != currentRunNumber) {
      uint64_t timestamp = cfgUseFixedTimestamp ? cfgTimestamp.value : bc.timestamp();
      pidModel = PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value, ccdbApi, timestamp, cfgPid.value, static_cast<PidMLDetector>(cfgDetector.value), cfgCertainty.value);
    }

    for (auto& mcPart : mcParticles) {
      // eta cut is included in requireGlobalTrackInFilter() so we cut it only here
      if(mcPart.isPhysicalPrimary() && TMath::Abs(mcPart.eta()) < ETA_CUT && mcPart.pdgCode() == pidModel.mPid) {
        histos.fill(HIST("hPtMCPositive"), mcPart.pt());
      }
    }

    for (auto& track : tracks) {
      if(track.has_mcParticle()) {
        auto mcPart = track.mcParticle();
        if(mcPart.isPhysicalPrimary()) {
          bool mlAccepted = pidModel.applyModelBoolean(track);
          bool nSigmaAccepted = IsNSigmaAccept(track);
          bool tofPtCutted = cfgDetector != 0 && track.pt() <= cfgTofPtCut;

          LOGF(info, "collision id: %d track id: %d mlAccepted: %d nSigmaAccepted: %d tofPtCutted: %d p: %.3f; x: %.3f, y: %.3f, z: %.3f",
              track.collisionId(), track.index(), mlAccepted, nSigmaAccepted, tofPtCutted, track.p(), track.x(), track.y(), track.z());

          if(mcPart.pdgCode() == pidModel.mPid) {
            histos.fill(HIST("hPtMCTracked"), mcPart.pt());
          }
          
          histos.fill(HIST("hPtTOFBeta"), track.pt(), track.beta());
          histos.fill(HIST("hPtTPCSignal"), track.pt(), track.tpcSignal());

          if(mlAccepted && !tofPtCutted) {
            if(mcPart.pdgCode() == pidModel.mPid) {
              histos.fill(HIST("hPtMLTruePositive"), mcPart.pt());
            }
            histos.fill(HIST("hPtMLPositive"), mcPart.pt());
          }

          if(nSigmaAccepted && !tofPtCutted) {
            if(mcPart.pdgCode() == pidModel.mPid) {
              histos.fill(HIST("hPtNSigmaTruePositive"), mcPart.pt());
            }
            histos.fill(HIST("hPtNSigmaPositive"), mcPart.pt());
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MlModelGenHists>(cfgc)};
}
