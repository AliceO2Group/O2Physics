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

/// \file pidMLBatchEffAndPurProducer
/// \brief Batch PID execution task. It produces derived data needed for ROOT script which
/// generates efficiency (recall) and purity (precision) analysis of ML Model PID
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

#include <string_view>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace effandpurpidresult
{
DECLARE_SOA_INDEX_COLUMN(Track, track);              //! Track index
DECLARE_SOA_COLUMN(Pid, pid, int);                   //! PDG particle ID to be tested by the model
DECLARE_SOA_COLUMN(Pt, pt, float);                   //! particle's pt
DECLARE_SOA_COLUMN(MlCertainty, mlCertainty, float); //! Machine learning model certainty value for track and pid
DECLARE_SOA_COLUMN(NSigma, nSigma, float);     //! nSigma value for track and pid
DECLARE_SOA_COLUMN(IsPidMC, isPidMc, bool);          //! Is track's mcParticle recognized as "Pid"
} // namespace effandpurpidresult

DECLARE_SOA_TABLE(EffAndPurPidResult, "AOD", "PIDEFFANDPURRES", o2::soa::Index<>,
                  effandpurpidresult::TrackId, effandpurpidresult::Pid, effandpurpidresult::Pt, effandpurpidresult::MlCertainty,
                  effandpurpidresult::NSigma, effandpurpidresult::IsPidMC);
} // namespace o2::aod

struct PidMlBatchEffAndPurProducer {
  Produces<o2::aod::EffAndPurPidResult> effAndPurPIDResult;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  int currentRunNumber = -1;
  static constexpr float eps = 1e-10;
  static constexpr float etaCut = 0.8f;
  static constexpr float nSigmaTofPtCut = 0.5f;
  static constexpr float tofMissing = -999.0f;
  static constexpr int nPids = 6;
  static constexpr int pids[nPids] = {211, 321, 2212, -211, -321, -2212};

  o2::ccdb::CcdbApi ccdbApi;
  std::vector<PidONNXModel> models;

  Configurable<std::vector<int>> cfgPids{"pids", std::vector<int>(pids, pids + nPids), "PIDs to predict"};

  Configurable<std::string> cfgPathCCDB{"ccdb-path", "Users/m/mkabus/PIDML", "base path to the CCDB directory with ONNX models"};
  Configurable<std::string> cfgCCDBURL{"ccdb-url", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  Configurable<bool> cfgUseCCDB{"use-ccdb", true, "Whether to autofetch ML model from CCDB. If false, local file will be used."};
  Configurable<std::string> cfgPathLocal{"local-path", "/home/mkabus/PIDML/", "base path to the local directory with ONNX models"};

  Configurable<bool> cfgUseFixedTimestamp{"use-fixed-timestamp", false, "Whether to use fixed timestamp from configurable instead of timestamp calculated from the data"};
  Configurable<uint64_t> cfgTimestamp{"timestamp", 1524176895000, "Hardcoded timestamp for tests"};

  static constexpr std::string_view mcTrackedHistLabels[nPids] = {
    "211/hPtMCTracked",
    "321/hPtMCTracked",
    "2212/hPtMCTracked",
    "0211/hPtMCTracked",
    "0321/hPtMCTracked",
    "02212/hPtMCTracked"};
  template <std::size_t i, typename T>
  void fillMcParticleTrackedHistos(const T& track)
  {
    histos.fill(HIST(mcTrackedHistLabels[i]), track.pt());
  }

  static constexpr std::string_view mcPositiveHistLabels[nPids] = {
    "211/hPtMCPositive",
    "321/hPtMCPositive",
    "2212/hPtMCPositive",
    "0211/hPtMCPositive",
    "0321/hPtMCPositive",
    "02212/hPtMCPositive"};
  template <std::size_t i, typename T>
  void fillMcParticlePositiveHistos(const T& mcPart)
  {
    histos.fill(HIST(mcPositiveHistLabels[i]), mcPart.pt());
  }

  template <typename T>
  void fillMcParticlePositiveHist(const T& mcPart)
  {
    switch (mcPart.pdgCode()) {
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
  void fillMcParticleTrackedHist(const T& track, int pdgCode)
  {
    switch (pdgCode) {
      case pids[0]:
        fillMcParticleTrackedHistos<0>(track);
        break;
      case pids[1]:
        fillMcParticleTrackedHistos<1>(track);
        break;
      case pids[2]:
        fillMcParticleTrackedHistos<2>(track);
        break;
      case pids[3]:
        fillMcParticleTrackedHistos<3>(track);
        break;
      case pids[4]:
        fillMcParticleTrackedHistos<4>(track);
        break;
      case pids[5]:
        fillMcParticleTrackedHistos<5>(track);
        break;
      default:
        return;
    }
  }

  void init(InitContext const&)
  {
    if (cfgUseCCDB) {
      ccdbApi.init(cfgCCDBURL);
    } else {
    }

    const AxisSpec axisPt{50, 0, 3.1, "pt"};

    // Monte Carlo Model Simulation Tracking and PID truth info
    for (int i = 0; i < nPids; ++i) {
      histos.add(mcTrackedHistLabels[i].data(), mcTrackedHistLabels[i].data(), kTH1F, {axisPt});
      histos.add(mcPositiveHistLabels[i].data(), mcPositiveHistLabels[i].data(), kTH1F, {axisPt});
    }
  }

  Filter trackFilter = requireGlobalTrackInFilter();

  using BigTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal, aod::McTrackLabels,
                                            aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullMu,
                                            aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu>>;

  typedef struct nSigma_t {
    double tpc, tof, composed;
  } nSigma_t;

  const nSigma_t getNSigma(const BigTracks::iterator& track, const int& cfgPid)
  {
    nSigma_t nSigma;

    switch (TMath::Abs(cfgPid)) {
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

    if(track.pt() < nSigmaTofPtCut || (track.tofSignal() - tofMissing) < eps) {
      nSigma.composed = TMath::Abs(nSigma.tpc);
    } else {
      nSigma.composed = TMath::Hypot(nSigma.tof, nSigma.tpc);
    }

    int sign = cfgPid > 0 ? 1 : -1;
    if(sign != track.sign()) {
      nSigma.composed = std::numeric_limits<float>::max();
    }

    return nSigma;
  }

  void process(aod::Collisions const& collisions, BigTracks const& tracks, aod::BCsWithTimestamps const&, aod::McParticles const& mcParticles)
  {
    auto bc = collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
    if (cfgUseCCDB && bc.runNumber() != currentRunNumber) {
      uint64_t timestamp = cfgUseFixedTimestamp ? cfgTimestamp.value : bc.timestamp();
      for (const int& pid : cfgPids.value)
        models.emplace_back(PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value,
                                         ccdbApi, timestamp, pid, 1.1));
    } else {
      for (int& pid : cfgPids.value)
        models.emplace_back(PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value,
                                         ccdbApi, -1, pid, 1.1));
    }

    for (auto& mcPart : mcParticles) {
      // eta cut is done in requireGlobalTrackInFilter() so we cut it only here
      if (mcPart.isPhysicalPrimary() && TMath::Abs(mcPart.eta()) < etaCut) {
        fillMcParticlePositiveHist(mcPart);
      }
    }

    for (auto& track : tracks) {
      if (track.has_mcParticle()) {
        auto mcPart = track.mcParticle();
        if (mcPart.isPhysicalPrimary()) {
          fillMcParticleTrackedHist(track, mcPart.pdgCode());

          for (int i = 0; i < cfgPids.value.size(); ++i) {
            float mlCertainty = models[i].applyModel(track);
            nSigma_t nSigma = getNSigma(track, cfgPids.value[i]);
            bool isMCPid = mcPart.pdgCode() == cfgPids.value[i];

            effAndPurPIDResult(track.index(), cfgPids.value[i], track.pt(), mlCertainty, nSigma.composed, isMCPid);
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PidMlBatchEffAndPurProducer>(cfgc)};
}
