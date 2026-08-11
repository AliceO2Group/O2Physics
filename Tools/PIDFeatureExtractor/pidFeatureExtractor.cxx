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
/// \brief Produce flat, ML-ready PID feature files (ROOT TTree and/or CSV)
///        from ALICE Run 3 Pb-Pb AO2D data, for both MC (reconstructed +
///        truth) and real/raw data.
///
///        Output is written via manual TFile/TTree/ofstream rather than
///        O2's DECLARE_SOA_TABLE/Produces<> table mechanism. That's a
///        deliberate reversion: two Produces<> tables sharing a column
///        prefix in one struct triggered a reproducible framework-level
///        compile failure (ASoA.h/MetadataTrait constraint-satisfaction
///        errors, and a StructToTuple reflection failure) against this O2
///        build, independent of table description tag naming. This
///        TFile/TTree approach is the same pattern the original two-file
///        (MC/RAW) version of this task used successfully.
///
/// \author Robert Forynski

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>

#include <TFile.h>
#include <TTree.h>

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
constexpr int kNumItsLayers = 7;
constexpr int kBitsPerItsLayer = 4;
constexpr uint32_t kItsLayerMask = 0xF;
constexpr int kNumSpecies = 4; // pi, ka, pr, el
constexpr float kEmcalEtaOutOfAcceptance = -900.f;

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
  for (int layer = 0; layer < kNumItsLayers; layer++) {
    if ((v >> (layer * kBitsPerItsLayer)) & kItsLayerMask) {
      n++;
    }
  }
  return n;
}
} // namespace

/// PidFeatureExtractor: flat PID feature file (ROOT TTree and/or CSV) for
/// ML training/inference.
///
/// Mode (MC vs. real data) is a runtime PROCESS_SWITCH choice, so one
/// executable serves both use cases.
///
/// - DPG track cuts (eta/pT/DCA/TPC-cluster/ITS-cluster) are optional and off
///   by default (wide-open ranges) - tighten them in the config if you want
///   them applied here instead of in Python post-processing.
/// - Bayesian PID combination is optional (computeBayesianPid, default
///   true) and configurable priors (bayesianPriors, default flat). Valid
///   whenever TPC is present; TOF is folded in too if also present, but is
///   not required - a TPC-only track still gets a real posterior, not NaN.
/// - Output: ROOT TTree (exportROOT, default true) and/or CSV (exportCsv,
///   default false), written via a single row of member variables bound as
///   TTree branches - not an O2 AOD table.
struct PidFeatureExtractor {
  std::unique_ptr<TFile> outputFile;
  std::unique_ptr<TTree> featureTree;
  std::ofstream csvFile;

  // --- output row (bound as TTree branches; also used to build CSV rows) ---
  float p = 0, pt = 0, px = 0, py = 0, pz = 0, eta = 0, phi = 0, sign = 0;
  int trackType = 0;
  float vz = 0, centFT0C = 0, dcaXY = 0, dcaZ = 0;
  bool hasTpc = false;
  float tpcSignal = 0, tpcNSigmaPi = 0, tpcNSigmaKa = 0, tpcNSigmaPr = 0, tpcNSigmaEl = 0;
  int tpcNClsFound = 0;
  float tpcChi2NCl = 0;
  bool hasTof = false;
  float tofMass = 0, beta = 0, tofNSigmaPi = 0, tofNSigmaKa = 0, tofNSigmaPr = 0, tofNSigmaEl = 0;
  bool hasTrd = false;
  float trdSignal = 0, trdChi2 = 0;
  int trdPattern = 0;
  int itsClusterSizes = 0;
  float itsChi2NCl = 0;
  bool hasEmcal = false;
  float trackEtaEmcal = 0, trackPhiEmcal = 0;
  bool hasHmpid = false;
  float hmpidSignal = 0, hmpidQMip = 0;
  int hmpidNPhotons = 0, hmpidClusSize = 0;
  float hmpidMom = 0;
  float bayesProbPi = 0, bayesProbKa = 0, bayesProbPr = 0, bayesProbEl = 0;
  int mcPdg = 0;
  uint8_t mcIsPhysicalPrimary = 0;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Filter trackFilter = requireGlobalTrackInFilter();

  // --- output configuration --------------------------------------------------
  Configurable<std::string> outputPath{"outputPath", "pid_features", "Output file base name (no extension)"};
  Configurable<bool> exportROOT{"exportROOT", true, "Write a ROOT TTree"};
  Configurable<bool> exportCsv{"exportCsv", false, "Also write a CSV, alongside the ROOT output"};

  // --- DPG cuts: wide-open by default, i.e. effectively disabled -----------
  Configurable<float> etaMin{"etaMin", -99.f, "Minimum track eta (DPG cut; wide-open = disabled)"};
  Configurable<float> etaMax{"etaMax", 99.f, "Maximum track eta (DPG cut; wide-open = disabled)"};
  Configurable<float> ptMin{"ptMin", 0.f, "Minimum track pT, GeV/c (DPG cut; wide-open = disabled)"};
  Configurable<float> ptMax{"ptMax", 9999.f, "Maximum track pT, GeV/c (DPG cut; wide-open = disabled)"};
  Configurable<float> dcaXYMax{"dcaXYMax", 9999.f, "Maximum |DCAxy|, cm (DPG cut; wide-open = disabled)"};
  Configurable<float> dcaZMax{"dcaZMax", 9999.f, "Maximum |DCAz|, cm (DPG cut; wide-open = disabled)"};
  Configurable<int> itsMinClusters{"itsMinClusters", 0, "Minimum number of ITS clusters (DPG cut; 0 = disabled)"};
  Configurable<int> tpcMinClusters{"tpcMinClusters", 0, "Minimum TPC clusters (DPG cut; 0 = disabled)"};

  // --- Bayesian PID ----------------------------------------------------------
  Configurable<bool> computeBayesianPid{"computeBayesianPid", true, "Compute Bayesian PID posteriors (else NaN)"};
  Configurable<std::vector<float>> bayesianPriors{"bayesianPriors", std::vector<float>{1.f, 1.f, 1.f, 1.f}, "Priors [pi,ka,pr,el]; default flat"};

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
    std::string base = outputPath.value;

    if (exportROOT.value) {
      outputFile = std::make_unique<TFile>((base + ".root").c_str(), "RECREATE");
      featureTree = std::make_unique<TTree>("pid_features", "PID features");

      featureTree->Branch("p", &p);
      featureTree->Branch("pt", &pt);
      featureTree->Branch("px", &px);
      featureTree->Branch("py", &py);
      featureTree->Branch("pz", &pz);
      featureTree->Branch("eta", &eta);
      featureTree->Branch("phi", &phi);
      featureTree->Branch("sign", &sign);
      featureTree->Branch("trackType", &trackType);
      featureTree->Branch("vz", &vz);
      featureTree->Branch("centFT0C", &centFT0C);
      featureTree->Branch("dcaXY", &dcaXY);
      featureTree->Branch("dcaZ", &dcaZ);
      featureTree->Branch("hasTPC", &hasTpc);
      featureTree->Branch("tpcSignal", &tpcSignal);
      featureTree->Branch("tpcNSigmaPi", &tpcNSigmaPi);
      featureTree->Branch("tpcNSigmaKa", &tpcNSigmaKa);
      featureTree->Branch("tpcNSigmaPr", &tpcNSigmaPr);
      featureTree->Branch("tpcNSigmaEl", &tpcNSigmaEl);
      featureTree->Branch("tpcNClsFound", &tpcNClsFound);
      featureTree->Branch("tpcChi2NCl", &tpcChi2NCl);
      featureTree->Branch("hasTOF", &hasTof);
      featureTree->Branch("tofMass", &tofMass);
      featureTree->Branch("beta", &beta);
      featureTree->Branch("tofNSigmaPi", &tofNSigmaPi);
      featureTree->Branch("tofNSigmaKa", &tofNSigmaKa);
      featureTree->Branch("tofNSigmaPr", &tofNSigmaPr);
      featureTree->Branch("tofNSigmaEl", &tofNSigmaEl);
      featureTree->Branch("hasTRD", &hasTrd);
      featureTree->Branch("trdSignal", &trdSignal);
      featureTree->Branch("trdChi2", &trdChi2);
      featureTree->Branch("trdPattern", &trdPattern);
      featureTree->Branch("itsClusterSizes", &itsClusterSizes);
      featureTree->Branch("itsChi2NCl", &itsChi2NCl);
      featureTree->Branch("hasEMCal", &hasEmcal);
      featureTree->Branch("trackEtaEmcal", &trackEtaEmcal);
      featureTree->Branch("trackPhiEmcal", &trackPhiEmcal);
      featureTree->Branch("hasHMPID", &hasHmpid);
      featureTree->Branch("hmpidSignal", &hmpidSignal);
      featureTree->Branch("hmpidQMip", &hmpidQMip);
      featureTree->Branch("hmpidNPhotons", &hmpidNPhotons);
      featureTree->Branch("hmpidClusSize", &hmpidClusSize);
      featureTree->Branch("hmpidMom", &hmpidMom);
      featureTree->Branch("bayesProbPi", &bayesProbPi);
      featureTree->Branch("bayesProbKa", &bayesProbKa);
      featureTree->Branch("bayesProbPr", &bayesProbPr);
      featureTree->Branch("bayesProbEl", &bayesProbEl);
      if (doprocessMc) {
        featureTree->Branch("mcPdg", &mcPdg);
        featureTree->Branch("mcIsPhysicalPrimary", &mcIsPhysicalPrimary);
      }
    }

    if (exportCsv.value) {
      csvFile.open(base + (doprocessMc ? "_mc.csv" : "_data.csv"));
      csvFile << "p,pt,px,py,pz,eta,phi,sign,trackType,"
                 "vz,centFT0C,dcaXY,dcaZ,"
                 "hasTPC,tpcSignal,tpcNSigmaPi,tpcNSigmaKa,tpcNSigmaPr,tpcNSigmaEl,tpcNClsFound,tpcChi2NCl,"
                 "hasTOF,tofMass,beta,tofNSigmaPi,tofNSigmaKa,tofNSigmaPr,tofNSigmaEl,"
                 "hasTRD,trdSignal,trdChi2,trdPattern,"
                 "itsClusterSizes,itsChi2NCl,"
                 "hasEMCal,trackEtaEmcal,trackPhiEmcal,"
                 "hasHMPID,hmpidSignal,hmpidQMip,hmpidNPhotons,hmpidClusSize,hmpidMom,"
                 "bayesProbPi,bayesProbKa,bayesProbPr,bayesProbEl";
      if (doprocessMc) {
        csvFile << ",mcPdg,mcIsPhysicalPrimary";
      }
      csvFile << "\n";
    }

    const AxisSpec axisPt{200, 0, 10, "pT"};
    const AxisSpec axisEta{60, -1.5, 1.5, "eta"};
    const AxisSpec axisdEdx{300, 0, 300, "dE/dx"};
    const AxisSpec axisBeta{120, 0, 1.2, "beta"};
    const AxisSpec axisMass{100, -0.2, 2.0, "mass"};
    histos.add("QC/nTracks", "Tracks", kTH1F, {{10000, 0, 100000}});
    histos.add("QC/pt", "pT", kTH1F, {axisPt});
    histos.add("QC/eta", "eta", kTH1F, {axisEta});
    histos.add("QC/tpcDedxVsPt", "dE/dx vs pT", kTH2F, {axisPt, axisdEdx});
    histos.add("QC/tofBetaVsP", "beta vs p", kTH2F, {axisPt, axisBeta});
    histos.add("QC/massVsP", "mass vs p", kTH2F, {axisPt, axisMass});
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
  void computeBayesianProbs(bool hasTpcIn, const float nsTPC[4], bool hasTofIn, const float nsTOF[4], float out[4]) const
  {
    if (!computeBayesianPid.value || !hasTpcIn) {
      out[0] = out[1] = out[2] = out[3] = kNaN;
      return;
    }
    auto const& priors = bayesianPriors.value;
    float sum = 0.f;
    for (int i = 0; i < kNumSpecies; i++) {
      float logL = -0.5f * nsTPC[i] * nsTPC[i];
      if (hasTofIn) {
        logL += -0.5f * nsTOF[i] * nsTOF[i];
      }
      out[i] = std::exp(logL) * priors[i];
      sum += out[i];
    }
    for (int i = 0; i < kNumSpecies; i++) {
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

  /// Fills the member "output row" for one track. Identical for MC and
  /// data - the only thing that differs between the two modes is whether
  /// mcPdg/mcIsPhysicalPrimary get set afterwards.
  template <typename TTrack>
  void fillRow(TTrack const& track, float vzIn, float centFT0CIn,
               std::unordered_map<int64_t, aod::HMPIDs::iterator> const& hmpidMap)
  {
    p = track.p();
    pt = track.pt();
    px = track.px();
    py = track.py();
    pz = track.pz();
    eta = track.eta();
    phi = track.phi();
    sign = static_cast<float>(track.sign());
    trackType = track.trackType();
    vz = vzIn;
    centFT0C = centFT0CIn;
    dcaXY = track.dcaXY();
    dcaZ = track.dcaZ();

    hasTpc = track.hasTPC();
    tpcSignal = track.tpcSignal();
    tpcNSigmaPi = track.tpcNSigmaPi();
    tpcNSigmaKa = track.tpcNSigmaKa();
    tpcNSigmaPr = track.tpcNSigmaPr();
    tpcNSigmaEl = track.tpcNSigmaEl();
    tpcNClsFound = track.tpcNClsFound();
    tpcChi2NCl = track.tpcChi2NCl();

    hasTof = !tofMissing(track);
    tofMass = getTofMass(track);
    beta = track.beta();
    tofNSigmaPi = track.tofNSigmaPi();
    tofNSigmaKa = track.tofNSigmaKa();
    tofNSigmaPr = track.tofNSigmaPr();
    tofNSigmaEl = track.tofNSigmaEl();

    hasTrd = !trdMissing(track);
    trdSignal = track.trdSignal();
    trdChi2 = track.trdChi2();
    trdPattern = track.trdPattern();

    itsClusterSizes = track.itsClusterSizes();
    itsChi2NCl = track.itsChi2NCl();

    hasEmcal = track.trackEtaEmcal() > kEmcalEtaOutOfAcceptance;
    trackEtaEmcal = track.trackEtaEmcal();
    trackPhiEmcal = track.trackPhiEmcal();

    hasHmpid = false;
    hmpidSignal = kNaN;
    hmpidQMip = kNaN;
    hmpidNPhotons = 0;
    hmpidClusSize = 0;
    hmpidMom = kNaN;
    if (auto it = hmpidMap.find(track.globalIndex()); it != hmpidMap.end()) {
      hasHmpid = true;
      hmpidSignal = it->second.hmpidSignal();
      hmpidQMip = it->second.hmpidQMip();
      hmpidNPhotons = it->second.hmpidNPhotons();
      hmpidClusSize = it->second.hmpidClusSize();
      hmpidMom = it->second.hmpidMom();
    }

    float nsTPC[4] = {tpcNSigmaPi, tpcNSigmaKa, tpcNSigmaPr, tpcNSigmaEl};
    float nsTOF[4] = {tofNSigmaPi, tofNSigmaKa, tofNSigmaPr, tofNSigmaEl};
    float bayes[4];
    computeBayesianProbs(hasTpc, nsTPC, hasTof, nsTOF, bayes);
    bayesProbPi = bayes[0];
    bayesProbKa = bayes[1];
    bayesProbPr = bayes[2];
    bayesProbEl = bayes[3];
  }

  void fillOutputs()
  {
    if (exportROOT.value) {
      featureTree->Fill();
    }
    if (exportCsv.value) {
      csvFile << p << ',' << pt << ',' << px << ',' << py << ',' << pz << ','
              << eta << ',' << phi << ',' << sign << ',' << trackType << ','
              << vz << ',' << centFT0C << ',' << dcaXY << ',' << dcaZ << ','
              << hasTpc << ',' << tpcSignal << ',' << tpcNSigmaPi << ',' << tpcNSigmaKa << ','
              << tpcNSigmaPr << ',' << tpcNSigmaEl << ',' << tpcNClsFound << ',' << tpcChi2NCl << ','
              << hasTof << ',' << tofMass << ',' << beta << ',' << tofNSigmaPi << ','
              << tofNSigmaKa << ',' << tofNSigmaPr << ',' << tofNSigmaEl << ','
              << hasTrd << ',' << trdSignal << ',' << trdChi2 << ',' << trdPattern << ','
              << itsClusterSizes << ',' << itsChi2NCl << ','
              << hasEmcal << ',' << trackEtaEmcal << ',' << trackPhiEmcal << ','
              << hasHmpid << ',' << hmpidSignal << ',' << hmpidQMip << ','
              << hmpidNPhotons << ',' << hmpidClusSize << ',' << hmpidMom << ','
              << bayesProbPi << ',' << bayesProbKa << ',' << bayesProbPr << ',' << bayesProbEl;
      if (doprocessMc) {
        csvFile << ',' << mcPdg << ',' << static_cast<int>(mcIsPhysicalPrimary);
      }
      csvFile << '\n';
    }
  }

  void fillQcHistos()
  {
    histos.fill(HIST("QC/nTracks"), 1);
    histos.fill(HIST("QC/pt"), pt);
    histos.fill(HIST("QC/eta"), eta);
    if (hasTpc) {
      histos.fill(HIST("QC/tpcDedxVsPt"), pt, tpcSignal);
    }
    if (hasTof) {
      histos.fill(HIST("QC/tofBetaVsP"), p, beta);
      histos.fill(HIST("QC/massVsP"), p, tofMass);
    }
  }

  void processData(PidCollision const& collision, PidTracks const& tracks, aod::HMPIDs const& hmpids)
  {
    auto hmpidMap = buildHmpidMap(hmpids);
    for (auto const& track : tracks) {
      if (!passesDpgCuts(track)) {
        continue;
      }
      fillRow(track, collision.posZ(), collision.centFT0C(), hmpidMap);
      fillOutputs();
      fillQcHistos();
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
      fillRow(track, collision.posZ(), collision.centFT0C(), hmpidMap);
      mcPdg = mcParticle.pdgCode();
      mcIsPhysicalPrimary = static_cast<uint8_t>(mcParticle.isPhysicalPrimary());
      fillOutputs();
      fillQcHistos();
    }
  }
  PROCESS_SWITCH(PidFeatureExtractor, processMc, "Produce PID features for MC (reconstructed + truth)", false);

  void finalize()
  {
    if (exportROOT.value) {
      outputFile->cd();
      featureTree->Write();
      outputFile->Close();
    }
    if (exportCsv.value) {
      csvFile.close();
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PidFeatureExtractor>(cfgc)};
}
