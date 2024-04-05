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

/// \file pidEffAndPurProducer
/// \brief Batch PID execution task, which prepares data for efficiency and purity analysis of ML Model PID
///
/// \author Michał Olędzki <mioledzk@cern.ch>
/// \author Marek Mytkowski <mmytkows@cern.ch>

#include <Framework/AnalysisDataModel.h>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "CCDB/CcdbApi.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Tools/PIDML/pidOnnxInterface.h"

#include <string_view>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod {
  namespace effandpurpidresult {
    DECLARE_SOA_INDEX_COLUMN(Track, track);               //! Track index
    DECLARE_SOA_COLUMN(Pid, pid, int);                    //! Pid to be tested by the model
    DECLARE_SOA_COLUMN(MlCertainty, mlCertainty, float);  //! Machine learning model cartainty value for track and pid
    DECLARE_SOA_COLUMN(TpcNSigma, tpcNSigma, float);      //! tpcNSigma value for track and pid
    DECLARE_SOA_COLUMN(TofNSigma, tofNSigma, float);      //! tofNSigma value for track and pid
    DECLARE_SOA_COLUMN(IsPidMC, isPidMc, bool);           //! Is track's mcParticle marked as pid
  } // namespace effandpurpidresult

  DECLARE_SOA_TABLE(EffAndPurPidResult, "AOD", "PIDEFFANDPURRES", o2::soa::Index<>,
    effandpurpidresult::TrackId, effandpurpidresult::Pid, effandpurpidresult::MlCertainty,
    effandpurpidresult::TpcNSigma, effandpurpidresult::TofNSigma, effandpurpidresult::IsPidMC);
} // namespace o2::aod


struct PidEffAndPurProducerONNXInterface {
  Produces<o2::aod::EffAndPurPidResult> effAndPurPIDResult;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  PidONNXInterface pidInterface;

  Configurable<uint32_t> cfgDetector{"detector", kTPCTOFTRD, "What detectors to use: 0: TPC only, 1: TPC + TOF, 2: TPC + TOF + TRD"};
  Configurable<double> cfgNSigmaCut{"n-sigma-cut", 3.0f, "TPC and TOF PID nSigma cut"};
  Configurable<double> cfgTofPtCut{"tof-pt-cut", 0.5f, "From what pT TOF is used"};

  Configurable<LabeledArray<double>> cfgPTCuts{"pT_cuts", {pidml_pt_cuts::cuts[0], pidml_pt_cuts::nPids, pidml_pt_cuts::nCutVars, pidml_pt_cuts::pidLabels, pidml_pt_cuts::cutVarLabels}, "pT cuts for each output pid and each detector configuration"};
  Configurable<std::vector<int>> cfgPids{"pids", std::vector<int>{pidml_pt_cuts::pids_v}, "PIDs to predict"};
  Configurable<bool> cfgAutoMode{"autoMode", true, "Use automatic model matching: default pT cuts and min certainties"};

  Configurable<std::string> cfgPathCCDB{"ccdb-path", "Users/m/mkabus/PIDML", "base path to the CCDB directory with ONNX models"};
  Configurable<std::string> cfgCCDBURL{"ccdb-url", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  Configurable<bool> cfgUseCCDB{"useCCDB", true, "Whether to autofetch ML model from CCDB. If false, local file will be used."};
  Configurable<std::string> cfgPathLocal{"local-path", "/home/mkabus/PIDML/", "base path to the local directory with ONNX models"};

  Configurable<bool> cfgUseFixedTimestamp{"use-fixed-timestamp", false, "Whether to use fixed timestamp from configurable instead of timestamp calculated from the data"};
  Configurable<uint64_t> cfgTimestamp{"timestamp", 1524176895000, "Hardcoded timestamp for tests"};

  o2::ccdb::CcdbApi ccdbApi;
  int currentRunNumber = -1;

  static constexpr float eps = 0.0001f;
  static constexpr float etaCut = 0.8f;
  static constexpr float tofAcceptThreshold = -999.0f;
  static constexpr float trdAcceptThreshold = 0.0f;

  static constexpr int nPids = 6;
  static constexpr int pids[nPids] = { 211, 321, 2212, -211, -321, -2212 };

  static constexpr std::string_view mcTrackedHistLabels[nPids] = { "211/hPtMCTracked", "321/hPtMCTracked", "2212/hPtMCTracked", "0211/hPtMCTracked", "0321/hPtMCTracked", "02212/hPtMCTracked" };
  template <std::size_t i, typename T>
  void fillMcParticleTrackedHistos(const T& mcPart) {
    histos.fill(HIST(mcTrackedHistLabels[i]), mcPart.pt());
  }

  static constexpr std::string_view mcPositiveHistLabels[nPids] = { "211/hPtMCPositive", "321/hPtMCPositive", "2212/hPtMCPositive", "0211/hPtMCPositive", "0321/hPtMCPositive", "02212/hPtMCPositive" };
  template <std::size_t i, typename T>
  void fillMcParticlePositiveHistos(const T& mcPart) {
    histos.fill(HIST(mcPositiveHistLabels[i]), mcPart.pt());
  }

  template <typename T>
  void fillMcParticlePositiveHist(const T& mcPart) {
    switch(mcPart.pdgCode()) {
      case pids[0]:
        fillMcParticlePositiveHistos<0>(mcPart);
        break;
      case pids[1]:
        fillMcParticlePositiveHistos<1>(mcPart);
        break;
      case pids[2]:
        fillMcParticlePositiveHistos<2>(mcPart);
        break;
      case pids[3]:
        fillMcParticlePositiveHistos<3>(mcPart);
        break;
      case pids[4]:
        fillMcParticlePositiveHistos<4>(mcPart);
        break;
      case pids[5]:
        fillMcParticlePositiveHistos<5>(mcPart);
        break;
      default:
        return;
    }
  }

  template <typename T>
  void fillMcParticleTrackedHist(const T& mcPart) {
    switch(mcPart.pdgCode()) {
      case pids[0]:
        fillMcParticleTrackedHistos<0>(mcPart);
        break;
      case pids[1]:
        fillMcParticleTrackedHistos<1>(mcPart);
        break;
      case pids[2]:
        fillMcParticleTrackedHistos<2>(mcPart);
        break;
      case pids[3]:
        fillMcParticleTrackedHistos<3>(mcPart);
        break;
      case pids[4]:
        fillMcParticleTrackedHistos<4>(mcPart);
        break;
      case pids[5]:
        fillMcParticleTrackedHistos<5>(mcPart);
        break;
      default:
        return;
    }
  }

  void init(InitContext const&) {
    if(cfgUseCCDB) {
      ccdbApi.init(cfgCCDBURL);
    } else {
      pidInterface = PidONNXInterface(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value,
          ccdbApi, -1, cfgPids.value, cfgPTCuts.value, pidml_pt_cuts::certainties_v, cfgAutoMode.value);
    }

    const AxisSpec axisPt{100, 0, 3.1, "pt"};

    // Monte Carlo Model Simulation Tracking and PID truth info
    for(int i = 0; i < nPids; ++i) {
      histos.add(mcTrackedHistLabels[i].data(), mcTrackedHistLabels[i].data(), kTH1F, {axisPt});
      histos.add(mcPositiveHistLabels[i].data(), mcPositiveHistLabels[i].data(), kTH1F, {axisPt});
    }
  }

  Filter trackFilter = requireGlobalTrackInFilter() &&
    (nabs(aod::pidtofsignal::tofSignal - tofAcceptThreshold) > eps) &&
    (nabs(aod::track::trdSignal - trdAcceptThreshold) > eps);

  using BigTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal, aod::McTrackLabels,
      aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullMu, 
      aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu>>;

  typedef struct nSigma_t {
    double tpc, tof; 
  } nSigma_t;

  const nSigma_t GetNSigma(const BigTracks::iterator &track, const int &cfgPid) {
    nSigma_t nSigma;

    switch(TMath::Abs(cfgPid)) {
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
  void process(aod::Collisions const& collisions, BigTracks const& tracks, aod::BCsWithTimestamps const&, aod::McParticles const& mcParticles) {
    auto bc = collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
    if(cfgUseCCDB && bc.runNumber() != currentRunNumber) {
      uint64_t timestamp = cfgUseFixedTimestamp ? cfgTimestamp.value : bc.timestamp();
      pidInterface = PidONNXInterface(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value,
          ccdbApi, timestamp, cfgPids.value, cfgPTCuts.value, pidml_pt_cuts::certainties_v, cfgAutoMode.value);
    }

    for(auto& mcPart : mcParticles) {
      // eta cut is included in requireGlobalTrackInFilter() so we cut it only here
      if(mcPart.isPhysicalPrimary() && TMath::Abs(mcPart.eta()) < etaCut) {
        fillMcParticlePositiveHist(mcPart);
      }
    }

    for(auto& track : tracks) {
      if(track.has_mcParticle()) {
        auto mcPart = track.mcParticle();
        if(mcPart.isPhysicalPrimary()) {
          fillMcParticleTrackedHist(mcPart);

          for(auto& pid : cfgPids.value) {
            float mlCertainty = pidInterface.applyModel(track, pid);
            nSigma_t nSigma = GetNSigma(track, pid);
            bool isMCPid = mcPart.pdgCode() == pid;
            
            effAndPurPIDResult(track.index(), pid, mlCertainty, nSigma.tpc, nSigma.tof, isMCPid); 
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) {
  return WorkflowSpec{
    adaptAnalysisTask<PidEffAndPurProducerONNXInterface>(cfgc)};
}
