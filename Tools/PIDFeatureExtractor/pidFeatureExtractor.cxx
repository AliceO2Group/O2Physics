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

/// \file pidFeatureExtractor.cxx
/// \brief Produce flat, ML-ready PID feature tables from ALICE Run 3 Pb-Pb
///        AO2D data, for both MC (reconstructed + truth) and real/raw data.
///        DPG track cuts and Bayesian PID are both optional and
///        configuration-driven (off/wide-open by default); CSV export is
///        available alongside the table output for convenience.
///
/// \author Robert Forynski

#include "pidFeatureExtractor.h"
//
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/runDataProcessing.h>

#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace
{
constexpr float kNaN = std::numeric_limits<float>::quiet_NaN();

/// Detector-presence helpers, local to this project.
template <typename T>
bool tofMissing(T const& track)
{
  return !track.hasTOF();
}

template <typename T>
bool trdMissing(T const& track)
{
  return !track.hasTRD();
}

template <typename T>
float getTofMass(T const& track)
{
  return tofMissing(track) ? kNaN : track.mass();
}

/// itsClusterSizes packs 7 ITS layers into 4 bits each (cluster size per
/// layer, 0 = no hit on that layer). Number of ITS clusters is the count of
/// non-zero nibbles, not the raw column value itself.
template <typename T>
int getItsNClusters(T const& track)
{
  auto v = static_cast<uint32_t>(track.itsClusterSizes());
  int n = 0;
  for (int layer = 0; layer < 7; layer++) {
    if ((v >> (layer * 4)) & 0xF) {
      n++;
    }
  }
  return n;
}

/// One row of reconstructed (non-MC) features, computed once per track and
/// shared between the table fill and the optional CSV row - keeps the two
/// outputs from ever being able to drift apart.
struct FeatureRow {
  float p, pt, px, py, pz, eta, phi, sign;
  int trackType;
  float vz, centFT0C, dcaXY, dcaZ;
  bool hasTpc;
  float tpcSignal, tpcNSigmaPi, tpcNSigmaKa, tpcNSigmaPr, tpcNSigmaEl;
  int tpcNClsFound;
  float tpcChi2NCl;
  bool hasTof;
  float tofMass, beta, tofNSigmaPi, tofNSigmaKa, tofNSigmaPr, tofNSigmaEl;
  bool hasTrd;
  float trdSignal, trdChi2;
  int trdPattern;
  int itsClusterSizes;
  float itsChi2NCl;
  bool hasEmcal;
  float trackEtaEmcal, trackPhiEmcal;
  bool hasHmpid;
  float hmpidSignal, hmpidQMip;
  int hmpidNPhotons, hmpidClusSize;
  float hmpidMom;
  float bayesProbPi, bayesProbKa, bayesProbPr, bayesProbEl;
};

/// CSV header, kept in exactly the same column order as FeatureRow / the
/// reconstructed part of the tables above, so a diff between the ROOT table
/// and the CSV is trivial if the two are ever compared.
constexpr const char* kCsvHeader =
  "p,pt,px,py,pz,eta,phi,sign,trackType,"
  "vz,centFT0C,dcaXY,dcaZ,"
  "hasTPC,tpcSignal,tpcNSigmaPi,tpcNSigmaKa,tpcNSigmaPr,tpcNSigmaEl,tpcNClsFound,tpcChi2NCl,"
  "hasTOF,tofMass,beta,tofNSigmaPi,tofNSigmaKa,tofNSigmaPr,tofNSigmaEl,"
  "hasTRD,trdSignal,trdChi2,trdPattern,"
  "itsClusterSizes,itsChi2NCl,"
  "hasEMCal,trackEtaEmcal,trackPhiEmcal,"
  "hasHMPID,hmpidSignal,hmpidQMip,hmpidNPhotons,hmpidClusSize,hmpidMom,"
  "bayesProbPi,bayesProbKa,bayesProbPr,bayesProbEl";

void writeCsvRow(std::ofstream& out, FeatureRow const& r)
{
  out << r.p << ',' << r.pt << ',' << r.px << ',' << r.py << ',' << r.pz << ','
      << r.eta << ',' << r.phi << ',' << r.sign << ',' << r.trackType << ','
      << r.vz << ',' << r.centFT0C << ',' << r.dcaXY << ',' << r.dcaZ << ','
      << r.hasTpc << ',' << r.tpcSignal << ',' << r.tpcNSigmaPi << ',' << r.tpcNSigmaKa << ','
      << r.tpcNSigmaPr << ',' << r.tpcNSigmaEl << ',' << r.tpcNClsFound << ',' << r.tpcChi2NCl << ','
      << r.hasTof << ',' << r.tofMass << ',' << r.beta << ',' << r.tofNSigmaPi << ','
      << r.tofNSigmaKa << ',' << r.tofNSigmaPr << ',' << r.tofNSigmaEl << ','
      << r.hasTrd << ',' << r.trdSignal << ',' << r.trdChi2 << ',' << r.trdPattern << ','
      << r.itsClusterSizes << ',' << r.itsChi2NCl << ','
      << r.hasEmcal << ',' << r.trackEtaEmcal << ',' << r.trackPhiEmcal << ','
      << r.hasHmpid << ',' << r.hmpidSignal << ',' << r.hmpidQMip << ','
      << r.hmpidNPhotons << ',' << r.hmpidClusSize << ',' << r.hmpidMom << ','
      << r.bayesProbPi << ',' << r.bayesProbKa << ',' << r.bayesProbPr << ',' << r.bayesProbEl << '\n';
}
} // namespace

/// PidFeatureExtractor: flat PID feature table for ML training/inference.
///
/// Mode (MC vs. real data) is a runtime PROCESS_SWITCH choice, so one
/// executable and one pair of output tables serve both use cases.
///
/// - DPG track cuts (eta/pT/DCA/TPC-cluster/ITS-cluster) are optional and off by
///   default (wide-open ranges) - tighten them in the config if you want
///   them applied here instead of in Python post-processing.
/// - Bayesian PID combination is optional (computeBayesianPid, default
///   true) and configurable priors (bayesianPriors, default flat). Valid
///   whenever TPC is present; TOF is folded in too if also present, but is
///   not required - a TPC-only track still gets a real posterior, not NaN.
/// - CSV export is optional (exportCsv, default false) and writes the same
///   reconstructed-feature rows as the ROOT table, alongside it.
/// - No histogramming - add a companion QA task later if needed, rather
///   than folding histograms into this producer.
struct PidFeatureExtractor {
  Produces<aod::PidFeaturesData> pidFeaturesData;
  Produces<aod::PidFeaturesMc> pidFeaturesMc;

  Filter trackFilter = requireGlobalTrackInFilter();

  // --- DPG cuts: wide-open by default, i.e. effectively disabled -----------
  Configurable<float> etaMin{"etaMin", -99.f, "Minimum track eta (DPG cut; wide-open = disabled)"};
  Configurable<float> etaMax{"etaMax", 99.f, "Maximum track eta (DPG cut; wide-open = disabled)"};
  Configurable<float> ptMin{"ptMin", 0.f, "Minimum track pT, GeV/c (DPG cut; wide-open = disabled)"};
  Configurable<float> ptMax{"ptMax", 9999.f, "Maximum track pT, GeV/c (DPG cut; wide-open = disabled)"};
  Configurable<float> dcaXYMax{"dcaxyMax", 9999.f, "Maximum |DCAxy|, cm (DPG cut; wide-open = disabled)"};
  Configurable<float> dcaZMax{"dcazMax", 9999.f, "Maximum |DCAz|, cm (DPG cut; wide-open = disabled)"};
  Configurable<int> itsMinClusters{"itsMinClusters", 0, "Minimum number of ITS clusters (DPG cut; 0 = disabled)"};
  Configurable<int> tpcMinClusters{"tpcMinClusters", 0, "Minimum TPC clusters (DPG cut; 0 = disabled)"};

  // --- Bayesian PID ----------------------------------------------------------
  Configurable<bool> computeBayesianPid{"computeBayesianPid", true, "Compute Bayesian PID posteriors (else NaN)"};
  Configurable<std::vector<float>> bayesianPriors{"bayesianPriors", {1.f, 1.f, 1.f, 1.f}, "Priors [pi,ka,pr,el]; default flat"};

  // --- CSV export --------------------------------------------------------------
  Configurable<bool> exportCsv{"exportCsv", false, "Also write reconstructed features to CSV alongside the table"};
  Configurable<std::string> csvOutputPath{"csvOutputPath", "pid_features", "CSV output file base name (no extension), if exportCsv"};

  std::unique_ptr<std::ofstream> csvFile;

  using PidTracks = soa::Filtered<soa::Join<
    aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
    aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCEl,
    aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFEl,
    aod::pidTOFmass, aod::pidTOFbeta>>;

  using PidTracksMc = soa::Filtered<soa::Join<
    aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
    aod::pidTPCPi, aod::pidTPCKa, aod::pidTPCPr, aod::pidTPCEl,
    aod::pidTOFPi, aod::pidTOFKa, aod::pidTOFPr, aod::pidTOFEl,
    aod::pidTOFmass, aod::pidTOFbeta, aod::McTrackLabels>>;

  using PidCollision = soa::Join<aod::Collisions, aod::CentFT0Cs>::iterator;

  void init(InitContext const&)
  {
    if (!exportCsv.value) {
      return;
    }
    // Only one of processData/processMc is expected to be active; open the
    // CSV that matches whichever is (doprocessXxx is generated by
    // PROCESS_SWITCH). If both or neither are on, nothing here enforces
    // that - see the README note on this.
    std::string suffix = doprocessMc ? "_mc.csv" : "_data.csv";
    csvFile = std::make_unique<std::ofstream>(csvOutputPath.value + suffix);
    *csvFile << kCsvHeader << '\n';
  }

  /// DPG-style track quality cuts. Wide-open defaults mean this is a no-op
  /// unless the config tightens them.
  template <typename TTrack>
  bool passesDpgCuts(TTrack const& track) const
  {
    if (track.pt() < ptMin.value || track.pt() > ptMax.value)
      return false;
    if (track.eta() < etaMin.value || track.eta() > etaMax.value)
      return false;
    if (std::abs(track.dcaXY()) > dcaXYMax.value)
      return false;
    if (std::abs(track.dcaZ()) > dcaZMax.value)
      return false;
    if (track.tpcNClsFound() < tpcMinClusters.value)
      return false;
    if (getItsNClusters(track) < itsMinClusters.value)
      return false;
    return true;
  }

  /// Bayesian PID: requires TPC (a TPC-only track still gets a real
  /// posterior); folds in TOF too when also present. NaN in all four
  /// outputs if TPC is absent or computeBayesianPid is false.
  void computeBayesianProbs(bool hasTpc, float nsTPC[4], bool hasTof, float nsTOF[4], float out[4]) const
  {
    if (!computeBayesianPid.value || !hasTpc) {
      out[0] = out[1] = out[2] = out[3] = kNaN;
      return;
    }
    auto const& priors = bayesianPriors.value;
    float sum = 0.f;
    for (int i = 0; i < 4; i++) {
      float logL = -0.5f * nsTPC[i] * nsTPC[i];
      if (hasTof) {
        logL += -0.5f * nsTOF[i] * nsTOF[i];
      }
      out[i] = std::exp(logL) * priors[i];
      sum += out[i];
    }
    for (int i = 0; i < 4; i++) {
      out[i] = sum > 0.f ? out[i] / sum : 0.25f;
    }
  }

  /// HMPID is sparse (~0.1% of tracks matched) and linked by track global
  /// index rather than joinable 1:1, so it's looked up once per collision
  /// instead of per track.
  static std::unordered_map<int64_t, aod::HMPIDs::iterator> buildHmpidMap(aod::HMPIDs const& hmpids)
  {
    std::unordered_map<int64_t, aod::HMPIDs::iterator> map;
    for (auto h = hmpids.begin(); h != hmpids.end(); ++h) {
      map[h.trackId()] = h;
    }
    return map;
  }

  template <typename TTrack>
  struct HmpidRow {
    bool has = false;
    float signal = kNaN, qMip = kNaN, mom = kNaN;
    int nPhotons = 0, clusSize = 0;

    static HmpidRow lookup(TTrack const& track, std::unordered_map<int64_t, aod::HMPIDs::iterator> const& map)
    {
      HmpidRow row;
      if (auto it = map.find(track.globalIndex()); it != map.end()) {
        row.has = true;
        row.signal = it->second.hmpidSignal();
        row.qMip = it->second.hmpidQMip();
        row.mom = it->second.hmpidMom();
        row.nPhotons = it->second.hmpidNPhotons();
        row.clusSize = it->second.hmpidClusSize();
      }
      return row;
    }
  };

  /// Builds the shared reconstructed-feature row for one track. Identical
  /// for MC and data - the only thing that differs between the two modes is
  /// whether MC-truth columns get appended afterwards.
  template <typename TTrack>
  FeatureRow buildFeatureRow(TTrack const& track, float vz, float centFT0C,
                             std::unordered_map<int64_t, aod::HMPIDs::iterator> const& hmpidMap) const
  {
    FeatureRow r{};
    r.p = track.p();
    r.pt = track.pt();
    r.px = track.px();
    r.py = track.py();
    r.pz = track.pz();
    r.eta = track.eta();
    r.phi = track.phi();
    r.sign = static_cast<float>(track.sign());
    r.trackType = track.trackType();
    r.vz = vz;
    r.centFT0C = centFT0C;
    r.dcaXY = track.dcaXY();
    r.dcaZ = track.dcaZ();

    r.hasTpc = track.hasTPC();
    r.tpcSignal = track.tpcSignal();
    r.tpcNSigmaPi = track.tpcNSigmaPi();
    r.tpcNSigmaKa = track.tpcNSigmaKa();
    r.tpcNSigmaPr = track.tpcNSigmaPr();
    r.tpcNSigmaEl = track.tpcNSigmaEl();
    r.tpcNClsFound = track.tpcNClsFound();
    r.tpcChi2NCl = track.tpcChi2NCl();

    r.hasTof = !tofMissing(track);
    r.tofMass = getTofMass(track);
    r.beta = track.beta();
    r.tofNSigmaPi = track.tofNSigmaPi();
    r.tofNSigmaKa = track.tofNSigmaKa();
    r.tofNSigmaPr = track.tofNSigmaPr();
    r.tofNSigmaEl = track.tofNSigmaEl();

    r.hasTrd = !trdMissing(track);
    r.trdSignal = track.trdSignal();
    r.trdChi2 = track.trdChi2();
    r.trdPattern = track.trdPattern();

    r.itsClusterSizes = track.itsClusterSizes();
    r.itsChi2NCl = track.itsChi2NCl();

    r.hasEmcal = track.trackEtaEmcal() > -900.f;
    r.trackEtaEmcal = track.trackEtaEmcal();
    r.trackPhiEmcal = track.trackPhiEmcal();

    auto hm = HmpidRow<TTrack>::lookup(track, hmpidMap);
    r.hasHmpid = hm.has;
    r.hmpidSignal = hm.signal;
    r.hmpidQMip = hm.qMip;
    r.hmpidNPhotons = hm.nPhotons;
    r.hmpidClusSize = hm.clusSize;
    r.hmpidMom = hm.mom;

    float nsTPC[4] = {r.tpcNSigmaPi, r.tpcNSigmaKa, r.tpcNSigmaPr, r.tpcNSigmaEl};
    float nsTOF[4] = {r.tofNSigmaPi, r.tofNSigmaKa, r.tofNSigmaPr, r.tofNSigmaEl};
    float bayes[4];
    computeBayesianProbs(r.hasTpc, nsTPC, r.hasTof, nsTOF, bayes);
    r.bayesProbPi = bayes[0];
    r.bayesProbKa = bayes[1];
    r.bayesProbPr = bayes[2];
    r.bayesProbEl = bayes[3];

    return r;
  }

  void processData(PidCollision const& collision, PidTracks const& tracks, aod::HMPIDs const& hmpids)
  {
    auto hmpidMap = buildHmpidMap(hmpids);
    for (auto const& track : tracks) {
      if (!passesDpgCuts(track)) {
        continue;
      }
      auto r = buildFeatureRow(track, collision.posZ(), collision.centFT0C(), hmpidMap);

      pidFeaturesData(
        r.p, r.pt, r.px, r.py, r.pz, r.eta, r.phi, r.sign, r.trackType,
        r.vz, r.centFT0C, r.dcaXY, r.dcaZ,
        r.hasTpc, r.tpcSignal, r.tpcNSigmaPi, r.tpcNSigmaKa, r.tpcNSigmaPr, r.tpcNSigmaEl,
        r.tpcNClsFound, r.tpcChi2NCl,
        r.hasTof, r.tofMass, r.beta, r.tofNSigmaPi, r.tofNSigmaKa, r.tofNSigmaPr, r.tofNSigmaEl,
        r.hasTrd, r.trdSignal, r.trdChi2, r.trdPattern,
        r.itsClusterSizes, r.itsChi2NCl,
        r.hasEmcal, r.trackEtaEmcal, r.trackPhiEmcal,
        r.hasHmpid, r.hmpidSignal, r.hmpidQMip, r.hmpidNPhotons, r.hmpidClusSize, r.hmpidMom,
        r.bayesProbPi, r.bayesProbKa, r.bayesProbPr, r.bayesProbEl);

      if (exportCsv.value) {
        writeCsvRow(*csvFile, r);
      }
    }
  }
  PROCESS_SWITCH(PidFeatureExtractor, processData, "Produce PID features for real/raw data (no MC truth)", true);

  void processMc(PidCollision const& collision, PidTracksMc const& tracks, aod::McParticles const&, aod::HMPIDs const& hmpids)
  {
    auto hmpidMap = buildHmpidMap(hmpids);
    for (auto const& track : tracks) {
      if (!passesDpgCuts(track)) {
        continue;
      }
      if (!track.has_mcParticle()) {
        continue;
      }
      auto mcParticle = track.mcParticle();
      auto r = buildFeatureRow(track, collision.posZ(), collision.centFT0C(), hmpidMap);

      pidFeaturesMc(
        r.p, r.pt, r.px, r.py, r.pz, r.eta, r.phi, r.sign, r.trackType,
        r.vz, r.centFT0C, r.dcaXY, r.dcaZ,
        r.hasTpc, r.tpcSignal, r.tpcNSigmaPi, r.tpcNSigmaKa, r.tpcNSigmaPr, r.tpcNSigmaEl,
        r.tpcNClsFound, r.tpcChi2NCl,
        r.hasTof, r.tofMass, r.beta, r.tofNSigmaPi, r.tofNSigmaKa, r.tofNSigmaPr, r.tofNSigmaEl,
        r.hasTrd, r.trdSignal, r.trdChi2, r.trdPattern,
        r.itsClusterSizes, r.itsChi2NCl,
        r.hasEmcal, r.trackEtaEmcal, r.trackPhiEmcal,
        r.hasHmpid, r.hmpidSignal, r.hmpidQMip, r.hmpidNPhotons, r.hmpidClusSize, r.hmpidMom,
        r.bayesProbPi, r.bayesProbKa, r.bayesProbPr, r.bayesProbEl,
        mcParticle.pdgCode(), static_cast<uint8_t>(mcParticle.isPhysicalPrimary()));

      if (exportCsv.value) {
        writeCsvRow(*csvFile, r);
      }
    }
  }
  PROCESS_SWITCH(PidFeatureExtractor, processMc, "Produce PID features for MC (reconstructed + truth)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PidFeatureExtractor>(cfgc)};
}
