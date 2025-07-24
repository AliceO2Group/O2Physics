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

/// \file pidMlEffAndPurProducer.cxx
/// \brief Produce pt histograms for tracks accepted by ML network and for MC mcParticles.
///
/// \author Michał Olędzki <m.oledzki@cern.ch>
/// \author Marek Mytkowski <marek.mytkowski@cern.ch>

#include "Tools/PIDML/pidOnnxModel.h"
#include "Tools/PIDML/pidUtils.h"
//
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/CcdbApi.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace pidml::pidutils;

struct PidMlEffAndPurProducer {
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<int> pdgPid{"pdgPid", 211, "PID to predict"};
  Configurable<double> nSigmaCut{"nSigmaCut", 3.0f, "TPC and TOF PID nSigma cut"};
  Configurable<MomentumLimitsMatrix> detectorMomentumLimits{"detectorMomentumLimits", MomentumLimitsMatrix(pidml_pt_cuts::defaultModelPLimits), "use {detector} when p >= y_{detector}: array of 3 doubles [y_TPC, y_TOF, y_TRD]"};
  Configurable<double> mlIdentCertaintyThreshold{"mlIdentCertaintyThreshold", 0.5, "Min certainty of the model to accept given mcPart to be of given kind"};

  Configurable<std::string> ccdbPath{"ccdbPath", "Users/m/mkabus/PIDML", "base path to the CCDB directory with ONNX models"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  Configurable<bool> useCcdb{"useCcdb", true, "Whether to autofetch ML model from CCDB. If false, local file will be used."};
  Configurable<std::string> localPath{"localPath", "/home/mkabus/PIDML", "base path to the local directory with ONNX models"};

  Configurable<bool> useFixedTimestamp{"useFixedTimestamp", false, "Whether to use fixed timestamp from configurable instead of timestamp calculated from the data"};
  Configurable<uint64_t> fixedTimestamp{"fixedTimestamp", 1524176895000, "Hardcoded timestamp for tests"};

  o2::ccdb::CcdbApi ccdbApi;
  int currentRunNumber = -1;

  Filter trackFilter = requireGlobalTrackInFilter();

  using BigTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal, aod::McTrackLabels,
                                            aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullMu,
                                            aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu>>;
  PidONNXModel<BigTracks> pidModel;

  typedef struct nSigma_t {
    double tpc, tof;
  } nSigma_t;

  nSigma_t getNSigma(const BigTracks::iterator& track)
  {
    nSigma_t nSigma;

    switch (std::abs(pdgPid)) {
      case 11: // electron
        nSigma.tof = track.tofNSigmaEl();
        nSigma.tpc = track.tpcNSigmaEl();
        break;
      case 13: // muon
        nSigma.tof = track.tofNSigmaMu();
        nSigma.tpc = track.tpcNSigmaMu();
        break;
      case 211: // pion
        nSigma.tof = track.tofNSigmaPi();
        nSigma.tpc = track.tpcNSigmaPi();
        break;
      case 321: // kaon
        nSigma.tof = track.tofNSigmaKa();
        nSigma.tpc = track.tpcNSigmaKa();
        break;
      case 2212: // proton
        nSigma.tof = track.tofNSigmaPr();
        nSigma.tpc = track.tpcNSigmaPr();
        break;
    }

    return nSigma;
  }

  bool isNSigmaAccept(const BigTracks::iterator& track, const nSigma_t& nSigma)
  {
    // FIXME: for current particles it works, but there are some particles,
    //  which can have different sign and pdgSign
    int sign = pdgPid > 0 ? 1 : -1;
    if (track.sign() != sign)
      return false;

    if (!inPLimit(track, detectorMomentumLimits.value[kTPCTOF]) || tofMissing(track)) {
      if (std::abs(nSigma.tpc) >= nSigmaCut)
        return false;
    } else {
      if (std::hypot(nSigma.tof, nSigma.tpc) >= nSigmaCut)
        return false;
    }

    return true;
  }

  void init(InitContext const&)
  {
    if (useCcdb) {
      ccdbApi.init(ccdbUrl);
    } else {
      pidModel = PidONNXModel<BigTracks>(localPath.value, ccdbPath.value, useCcdb.value, ccdbApi, -1,
                                         pdgPid.value, mlIdentCertaintyThreshold.value, &detectorMomentumLimits.value[0]);
    }

    const AxisSpec axisPt{100, 0, 5.0, "pt"};
    const AxisSpec axisP{100, 0, 5.0, "p"};
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
    histos.add("full/hPtTOFBeta", "full/hPtTOFBeta", kTH2F, {axisP, axisBeta});
    histos.add("full/hPtTPCSignal", "full/hPtTPCSignal", kTH2F, {axisP, axisTPCSignal});
    histos.add("full/hPtTOFNSigma", "full/hPtTOFNSigma", kTH2F, {axisP, axisNSigma});
    histos.add("full/hPtTPCNSigma", "full/hPtTPCNSigma", kTH2F, {axisP, axisNSigma});

    histos.add("hPtTOFNSigma", "hPtTOFNSigma", kTH2F, {axisP, axisNSigma});
    histos.add("hPtTPCNSigma", "hPtTPCNSigma", kTH2F, {axisP, axisNSigma});
  }

  void process(aod::Collisions const& collisions, BigTracks const& tracks, aod::BCsWithTimestamps const&, aod::McParticles const& mcParticles)
  {
    auto bc = collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
    if (useCcdb && bc.runNumber() != currentRunNumber) {
      uint64_t timestamp = useFixedTimestamp ? fixedTimestamp.value : bc.timestamp();
      pidModel = PidONNXModel<BigTracks>(localPath.value, ccdbPath.value, useCcdb.value, ccdbApi, timestamp,
                                         pdgPid.value, mlIdentCertaintyThreshold.value, &detectorMomentumLimits.value[0]);
    }

    for (const auto& mcPart : mcParticles) {
      // eta cut is included in requireGlobalTrackInFilter() so we cut it only here
      if (mcPart.isPhysicalPrimary() && std::abs(mcPart.eta()) < kGlobalEtaCut && mcPart.pdgCode() == pidModel.mPid) {
        histos.fill(HIST("hPtMCPositive"), mcPart.pt());
      }
    }

    for (const auto& track : tracks) {
      if (track.has_mcParticle()) {
        auto mcPart = track.mcParticle();
        if (mcPart.isPhysicalPrimary()) {
          bool mlAccepted = pidModel.applyModelBoolean(track);
          nSigma_t nSigma = getNSigma(track);
          bool nSigmaAccepted = isNSigmaAccept(track, nSigma);

          LOGF(debug, "collision id: %d track id: %d mlAccepted: %d nSigmaAccepted: %d p: %.3f; x: %.3f, y: %.3f, z: %.3f",
               track.collisionId(), track.index(), mlAccepted, nSigmaAccepted, track.p(), track.x(), track.y(), track.z());

          if (mcPart.pdgCode() == pidModel.mPid) {
            histos.fill(HIST("full/hPtTOFNSigma"), track.p(), nSigma.tof);
            histos.fill(HIST("full/hPtTPCNSigma"), track.p(), nSigma.tpc);
            histos.fill(HIST("hPtMCTracked"), track.pt());
          }

          histos.fill(HIST("full/hPtTOFBeta"), track.pt(), track.beta());
          histos.fill(HIST("full/hPtTPCSignal"), track.pt(), track.tpcSignal());

          if (mlAccepted) {
            if (mcPart.pdgCode() == pidModel.mPid) {
              histos.fill(HIST("hPtMLTruePositive"), track.pt());
            }
            histos.fill(HIST("hPtMLPositive"), track.pt());
          }

          if (nSigmaAccepted) {
            histos.fill(HIST("hPtTOFNSigma"), track.p(), nSigma.tof);
            histos.fill(HIST("hPtTPCNSigma"), track.p(), nSigma.tpc);

            if (mcPart.pdgCode() == pidModel.mPid) {
              histos.fill(HIST("hPtNSigmaTruePositive"), track.pt());
            }
            histos.fill(HIST("hPtNSigmaPositive"), track.pt());
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PidMlEffAndPurProducer>(cfgc)};
}
