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
/// \brief Batch PID execution task. It produces derived data needed for ROOT script, which
/// generates efficiency (recall) and purity (precision) analysis of ML Model PID
///
/// \author Michał Olędzki <m.oledzki@cern.ch>
/// \author Marek Mytkowski <marek.mytkowski@cern.ch>

#include <cstddef>
#include <string_view>
#include <algorithm>

#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/StaticFor.h"
#include "CCDB/CcdbApi.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "Tools/PIDML/pidOnnxModel.h"
#include "Tools/PIDML/pidUtils.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace pidml::pidutils;

namespace o2::aod
{
namespace effandpurpidresult
{
DECLARE_SOA_INDEX_COLUMN(Track, track);              //! Track index
DECLARE_SOA_COLUMN(Pid, pid, int32_t);               //! PDG particle ID to be tested by the model
DECLARE_SOA_COLUMN(Pt, pt, float);                   //! particle's pt
DECLARE_SOA_COLUMN(MlCertainty, mlCertainty, float); //! Machine learning model certainty value for track and pid
DECLARE_SOA_COLUMN(NSigma, nSigma, float);           //! nSigma value for track and pid
DECLARE_SOA_COLUMN(IsPidMC, isPidMc, bool);          //! Is track's mcParticle recognized as "Pid"
DECLARE_SOA_COLUMN(HasTOF, hasTof, bool);            //! Does track have TOF detector signal
DECLARE_SOA_COLUMN(HasTRD, hasTrd, bool);            //! Does track have TRD detector signal
} // namespace effandpurpidresult

DECLARE_SOA_TABLE(EffAndPurPidResult, "AOD", "PIDEFFANDPURRES", o2::soa::Index<>,
                  effandpurpidresult::TrackId, effandpurpidresult::Pid, effandpurpidresult::Pt, effandpurpidresult::MlCertainty,
                  effandpurpidresult::NSigma, effandpurpidresult::IsPidMC, effandpurpidresult::HasTOF, effandpurpidresult::HasTRD);
} // namespace o2::aod

struct PidMlBatchEffAndPurProducer {
  Produces<o2::aod::EffAndPurPidResult> effAndPurPIDResult;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr int32_t currentRunNumber = -1;
  static constexpr uint32_t kNPids = 6;
  static constexpr int32_t kPids[kNPids] = {2212, 321, 211, -211, -321, -2212};
  static constexpr std::string_view kParticleLabels[kNPids] = {"2212", "321", "211", "0211", "0321", "02212"};
  static constexpr std::string_view kParticleNames[kNPids] = {"proton", "kaon", "pion", "antipion", "antikaon", "antiproton"};

  std::array<std::shared_ptr<TH1>, kNPids> hTracked;
  std::array<std::shared_ptr<TH1>, kNPids> hMCPositive;

  o2::ccdb::CcdbApi ccdbApi;
  std::vector<PidONNXModel> models;

  Configurable<std::vector<int32_t>> cfgPids{"pids", std::vector<int32_t>(kPids, kPids + kNPids), "PIDs to predict"};
  Configurable<std::array<double, kNDetectors>> cfgDetectorsPLimits{"detectors-p-limits", std::array<double, kNDetectors>(pidml_pt_cuts::defaultModelPLimits), "\"use {detector} when p >= y_{detector}\": array of 3 doubles [y_TPC, y_TOF, y_TRD]"};
  Configurable<std::string> cfgPathCCDB{"ccdb-path", "Users/m/mkabus/PIDML", "base path to the CCDB directory with ONNX models"};
  Configurable<std::string> cfgCCDBURL{"ccdb-url", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  Configurable<bool> cfgUseCCDB{"use-ccdb", true, "Whether to autofetch ML model from CCDB. If false, local file will be used."};
  Configurable<std::string> cfgPathLocal{"local-path", "/home/mkabus/PIDML/", "base path to the local directory with ONNX models"};
  Configurable<bool> cfgUseFixedTimestamp{"use-fixed-timestamp", false, "Whether to use fixed timestamp from configurable instead of timestamp calculated from the data"};
  Configurable<uint64_t> cfgTimestamp{"timestamp", 1524176895000, "Hardcoded timestamp for tests"};

  Filter trackFilter = requireGlobalTrackInFilter();

  using BigTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal, aod::McTrackLabels,
                                            aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullMu,
                                            aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu>>;

  void initHistos()
  {
    static const AxisSpec axisPt{50, 0, 3.1, "pt"};

    static_for<0, kNPids - 1>([&](auto i) {
      if (std::find(cfgPids.value.begin(), cfgPids.value.end(), kPids[i]) != cfgPids.value.end()) {
        hTracked[i] = histos.add<TH1>(Form("%s/hPtMCTracked", kParticleLabels[i].data()), Form("Tracked %ss vs pT", kParticleNames[i].data()), kTH1F, {axisPt});
        hMCPositive[i] = histos.add<TH1>(Form("%s/hPtMCPositive", kParticleLabels[i].data()), Form("MC Positive %ss vs pT", kParticleNames[i].data()), kTH1F, {axisPt});
      }
    });
  }

  void init(InitContext const&)
  {
    if (cfgUseCCDB) {
      ccdbApi.init(cfgCCDBURL);
    }

    initHistos();
  }

  std::optional<size_t> getPartIndex(int32_t pdgCode)
  {
    std::optional<size_t> index;

    if (std::find(cfgPids.value.begin(), cfgPids.value.end(), pdgCode) != cfgPids.value.end()) {
      switch (pdgCode) {
        case 2212:
          index = 0;
          break;
        case 321:
          index = 1;
          break;
        case 211:
          index = 2;
          break;
        case -211:
          index = 3;
          break;
        case -321:
          index = 4;
          break;
        case -2212:
          index = 5;
          break;
      }
    }

    return index;
  }

  void fillTrackedHist(int32_t pdgCode, float pt)
  {
    auto ind = getPartIndex(pdgCode);
    if (ind) {
      hTracked[ind.value()]->Fill(pt);
    }
  }

  void fillMCPositiveHist(int32_t pdgCode, float pt)
  {
    auto ind = getPartIndex(pdgCode);
    if (ind) {
      hMCPositive[ind.value()]->Fill(pt);
    }
  }

  typedef struct nSigma_t {
    double tpc, tof, composed;
  } nSigma_t;

  const nSigma_t getNSigma(const BigTracks::iterator& track, const int32_t& cfgPid)
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

    if (!inPLimit(track, cfgDetectorsPLimits.value[kTPCTOF]) || tofMissing(track)) {
      nSigma.composed = TMath::Abs(nSigma.tpc);
    } else {
      nSigma.composed = TMath::Hypot(nSigma.tof, nSigma.tpc);
    }

    int32_t sign = cfgPid > 0 ? 1 : -1;
    if (sign != track.sign()) {
      nSigma.composed = std::numeric_limits<float>::max();
    }

    return nSigma;
  }

  void process(aod::Collisions const& collisions, BigTracks const& tracks, aod::BCsWithTimestamps const&, aod::McParticles const& mcParticles)
  {
    effAndPurPIDResult.reserve(mcParticles.size());

    auto bc = collisions.iteratorAt(0).bc_as<aod::BCsWithTimestamps>();
    if (cfgUseCCDB && bc.runNumber() != currentRunNumber) {
      uint64_t timestamp = cfgUseFixedTimestamp ? cfgTimestamp.value : bc.timestamp();
      for (const int32_t& pid : cfgPids.value)
        models.emplace_back(PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value,
                                         ccdbApi, timestamp, pid, 1.1, &cfgDetectorsPLimits.value[0]));
    } else {
      for (int32_t& pid : cfgPids.value)
        models.emplace_back(PidONNXModel(cfgPathLocal.value, cfgPathCCDB.value, cfgUseCCDB.value,
                                         ccdbApi, -1, pid, 1.1, &cfgDetectorsPLimits.value[0]));
    }

    for (auto& mcPart : mcParticles) {
      // eta cut is done in requireGlobalTrackInFilter() so we cut it only here
      if (mcPart.isPhysicalPrimary() && TMath::Abs(mcPart.eta()) < kGlobalEtaCut) {
        fillMCPositiveHist(mcPart.pdgCode(), mcPart.pt());
      }
    }

    for (auto& track : tracks) {
      if (track.has_mcParticle()) {
        auto mcPart = track.mcParticle();
        if (mcPart.isPhysicalPrimary()) {
          fillTrackedHist(mcPart.pdgCode(), track.pt());

          for (size_t i = 0; i < cfgPids.value.size(); ++i) {
            float mlCertainty = models[i].applyModel(track);
            nSigma_t nSigma = getNSigma(track, cfgPids.value[i]);
            bool isMCPid = mcPart.pdgCode() == cfgPids.value[i];

            effAndPurPIDResult(track.index(), cfgPids.value[i], track.pt(), mlCertainty, nSigma.composed, isMCPid, track.hasTOF(), track.hasTRD());
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
