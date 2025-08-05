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

/// \file pidMlBatchEffAndPurProducer.cxx
/// \brief Batch PID execution task. It produces derived data needed for ROOT script, which
/// generate PIDML neural network performance benchmark plots.
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
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

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
DECLARE_SOA_COLUMN(IsPidMC, isPidMC, bool);          //! Is track's mcParticle recognized as "Pid"
DECLARE_SOA_COLUMN(HasTOF, hasTOF, bool);            //! Does track have TOF detector signal
DECLARE_SOA_COLUMN(HasTRD, hasTRD, bool);            //! Does track have TRD detector signal
} // namespace effandpurpidresult

DECLARE_SOA_TABLE(EffAndPurPidResult, "AOD", "PIDEFFANDPURRES", o2::soa::Index<>,
                  effandpurpidresult::TrackId, effandpurpidresult::Pid, effandpurpidresult::Pt, effandpurpidresult::MlCertainty,
                  effandpurpidresult::NSigma, effandpurpidresult::IsPidMC, effandpurpidresult::HasTOF, effandpurpidresult::HasTRD);
} // namespace o2::aod

struct PidMlBatchEffAndPurProducer {
  Produces<o2::aod::EffAndPurPidResult> effAndPurPIDResult;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  static constexpr int32_t CurrentRunNumber = -1;
  static constexpr uint32_t KNPids = 6;
  static constexpr int32_t KPids[KNPids] = {2212, 321, 211, -211, -321, -2212};
  static constexpr std::string_view KParticleLabels[KNPids] = {"2212", "321", "211", "0211", "0321", "02212"};
  static constexpr std::string_view KPatricleNames[KNPids] = {"proton", "kaon", "pion", "antipion", "antikaon", "antiproton"};

  std::array<std::shared_ptr<TH1>, KNPids> hTracked;
  std::array<std::shared_ptr<TH1>, KNPids> hMCPositive;

  o2::ccdb::CcdbApi ccdbApi;

  Configurable<std::vector<int32_t>> pdgPids{"pdgPids", std::vector<int32_t>(KPids, KPids + KNPids), "list of PDG ids of particles to predict. Every subset of {2212, 321, 211, -211, -321, -2212} is correct value."};
  Configurable<MomentumLimitsMatrix> detectorMomentumLimits{"detectorMomentumLimits", MomentumLimitsMatrix(pidml_pt_cuts::defaultModelPLimits), "\"use {detector} when p >= y_{detector}\": array of 3 doubles [y_TPC, y_TOF, y_TRD]"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Users/m/mkabus/PIDML", "Base path to the CCDB directory with ONNX models"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "URL of the CCDB repository"};
  Configurable<bool> useCcdb{"useCcdb", true, "Whether to autofetch ML model from CCDB. If false, local file will be used."};
  Configurable<std::string> localPath{"localPath", "/home/mkabus/PIDML/", "Base path to the local directory with ONNX models"};
  Configurable<bool> useFixedTimestamp{"useFixedTimestamp", false, "Whether to use fixed timestamp from configurable instead of timestamp calculated from the data"};
  Configurable<uint64_t> fixedTimestamp{"fixedTimestamp", 1524176895000, "Hardcoded timestamp for tests"};

  Filter trackFilter = requireGlobalTrackInFilter();

  using BigTracks = soa::Filtered<soa::Join<aod::FullTracks, aod::TracksDCA, aod::pidTOFbeta, aod::TrackSelection, aod::TOFSignal, aod::McTrackLabels,
                                            aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullMu,
                                            aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu>>;
  std::vector<PidONNXModel<BigTracks>> models;

  void initHistos()
  {
    static const AxisSpec axisPt{50, 0, 3.1, "pt"};

    static_for<0, KNPids - 1>([&](auto i) {
      if (std::find(pdgPids.value.begin(), pdgPids.value.end(), KPids[i]) != pdgPids.value.end()) {
        hTracked[i] = histos.add<TH1>(Form("%s/hPtMCTracked", KParticleLabels[i].data()), Form("Tracked %ss vs pT", KPatricleNames[i].data()), kTH1F, {axisPt});
        hMCPositive[i] = histos.add<TH1>(Form("%s/hPtMCPositive", KParticleLabels[i].data()), Form("MC Positive %ss vs pT", KPatricleNames[i].data()), kTH1F, {axisPt});
      }
    });
  }

  void init(InitContext const&)
  {
    if (useCcdb) {
      ccdbApi.init(ccdbUrl);
    }

    initHistos();
  }

  std::optional<size_t> getPartIndex(int32_t pdgCode)
  {
    std::optional<size_t> index;

    if (std::find(pdgPids.value.begin(), pdgPids.value.end(), pdgCode) != pdgPids.value.end()) {
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

    switch (std::abs(cfgPid)) {
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

    if (!inPLimit(track, detectorMomentumLimits.value[kTPCTOF]) || tofMissing(track)) {
      nSigma.composed = std::abs(nSigma.tpc);
    } else {
      nSigma.composed = std::hypot(nSigma.tof, nSigma.tpc);
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
    if (useCcdb && bc.runNumber() != CurrentRunNumber) {
      uint64_t timestamp = useFixedTimestamp ? fixedTimestamp.value : bc.timestamp();
      for (const int32_t& pid : pdgPids.value)
        models.emplace_back(PidONNXModel<BigTracks>(localPath.value, ccdbPath.value, useCcdb.value,
                                                    ccdbApi, timestamp, pid, 1.1, &detectorMomentumLimits.value[0]));
    } else {
      for (const int32_t& pid : pdgPids.value)
        models.emplace_back(PidONNXModel<BigTracks>(localPath.value, ccdbPath.value, useCcdb.value,
                                                    ccdbApi, -1, pid, 1.1, &detectorMomentumLimits.value[0]));
    }

    for (const auto& mcPart : mcParticles) {
      // eta cut is done in requireGlobalTrackInFilter() so we cut it only here
      if (mcPart.isPhysicalPrimary() && std::abs(mcPart.eta()) < kGlobalEtaCut) {
        fillMCPositiveHist(mcPart.pdgCode(), mcPart.pt());
      }
    }

    for (const auto& track : tracks) {
      if (track.has_mcParticle()) {
        auto mcPart = track.mcParticle();
        if (mcPart.isPhysicalPrimary()) {
          fillTrackedHist(mcPart.pdgCode(), track.pt());

          for (size_t i = 0; i < pdgPids.value.size(); ++i) {
            float mlCertainty = models[i].applyModel(track);
            nSigma_t nSigma = getNSigma(track, pdgPids.value[i]);
            bool isMCPid = mcPart.pdgCode() == pdgPids.value[i];

            effAndPurPIDResult(track.index(), pdgPids.value[i], track.pt(), mlCertainty, nSigma.composed, isMCPid, track.hasTOF(), track.hasTRD());
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
