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
/// \file onTheFlyTrackerPid.cxx
///
/// \brief This task simulates the PID signal that can be obtained from the tracker layers as ToT (Time-over-Threshold).
///        The dE/dx of every track is sampled from a per-track LUT made from the full simulation.
///        It contemplates 9 particle types: electrons, muons, pions, kaons, protons, deuterons, tritons,
///        helium-3 and alphas/helium-4.
///
/// \author Henrik Fribert TUM
/// \author Nicolò Jacazio Università del Piemonte Orientale
/// \since  May 22, 2025
///

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TRandom3.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <format>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2;
using namespace o2::framework;

static constexpr double kMinMomentumForLogBins = 0.05;
static constexpr double kMaxMomentumForLogBins = 10.0;

/// Local path of a LUT file; "ccdb:<path>" is downloaded once to /tmp/ToTLUTs/<path>/snapshot.root
std::string getLocalLutFile(const std::string& filename, o2::ccdb::BasicCCDBManager* ccdbManager)
{
  if (!filename.starts_with("ccdb:")) {
    return filename;
  }
  const std::string ccdbPath = filename.substr(5);
  const std::string outPath = "/tmp/ToTLUTs/";
  const std::string localFilename = outPath + ccdbPath + "/snapshot.root";
  if (!std::ifstream(localFilename).is_open()) {
    if (!ccdbManager) {
      LOG(fatal) << "CCDB manager not set. Please set it before loading LUT from CCDB.";
    }
    std::map<std::string, std::string> metadata;
    LOG(info) << "Retrieving " << localFilename << " from CCDB path: " << ccdbPath;
    ccdbManager->getCCDBAccessor().retrieveBlob(ccdbPath, outPath, metadata, 1);
    if (!std::ifstream(localFilename).is_open()) {
      LOG(fatal) << "Could not find downloaded CCDB file " << localFilename;
    }
  }
  return localFilename;
}

/// Per-track dE/dx LUT; for every species, |eta| bin and rigidity bin the
/// distribution of the track dE/dx in e-, as a histogram. The dE/dx of a track is sampled from it,
/// interpolated linearly in log(rigidity) between the rigidity bin centres with the same random number.
class TrackDedxLut
{
 public:
  template <size_t N>
  bool load(const std::string& filename, o2::ccdb::BasicCCDBManager* ccdbManager, const std::array<int, N>& pdgCodes)
  {
    std::unique_ptr<TFile> f(TFile::Open(getLocalLutFile(filename, ccdbManager).c_str()));
    if (!f || f->IsZombie()) {
      LOG(fatal) << "Failed to open dE/dx LUT file: " << filename;
      return false;
    }
    const auto* etaBins = f->Get<TH1>("etaBins");
    const auto* rigBins = f->Get<TH1>("rigBins");
    const auto* dedxBins = f->Get<TH1>("dedxBins");
    if (!etaBins || !rigBins || !dedxBins) {
      LOG(fatal) << "No etaBins/rigBins/dedxBins in dE/dx LUT file: " << filename;
      return false;
    }
    mEtaEdges = binEdges(etaBins);
    mRigEdges = binEdges(rigBins);
    mDedxEdges = binEdges(dedxBins);
    mLogRigCentres.clear();
    for (size_t i = 0; i + 1 < mRigEdges.size(); ++i) {
      mLogRigCentres.push_back(0.5 * (std::log(mRigEdges[i]) + std::log(mRigEdges[i + 1])));
    }
    const size_t nBins = (mEtaEdges.size() - 1) * (mRigEdges.size() - 1);
    const size_t nDedx = mDedxEdges.size() - 1;

    bool allLoaded = true;
    for (const int pdg : pdgCodes) {
      std::unique_ptr<std::vector<float>> h(f->Get<std::vector<float>>(Form("h_%d", pdg)));
      if (!h || h->size() != nBins * nDedx) {
        LOG(warning) << "No dE/dx LUT for PDG " << pdg << " in " << filename;
        allLoaded = false;
        continue;
      }
      // cumulative distributions, only for the filled bins
      auto& offset = mOffset[pdg];
      auto& cdf = mCdf[pdg];
      offset.assign(nBins, -1);
      for (size_t b = 0; b < nBins; ++b) {
        const float* hb = h->data() + b * nDedx;
        double sum = 0.;
        for (size_t k = 0; k < nDedx; ++k) {
          sum += hb[k];
        }
        if (sum <= 0.) {
          continue;
        }
        offset[b] = cdf.size();
        double cum = 0.;
        for (size_t k = 0; k < nDedx; ++k) {
          cum += hb[k];
          cdf.push_back(cum / sum);
        }
      }
    }
    LOG(info) << "Loaded dE/dx LUT " << filename << " (" << mEtaEdges.size() - 1 << " |eta| x " << mRigEdges.size() - 1
              << " rigidity x " << mDedxEdges.size() - 1 << " dE/dx bins)";
    return allLoaded;
  }

  /// Sampled dE/dx for uniform random number u, or -1 if there is no LUT entry
  float sample(int pdg, float eta, float rigidity, double u) const
  {
    const auto it = mOffset.find(pdg);
    const float absEta = std::abs(eta);
    if (it == mOffset.end() || mEtaEdges.empty() || absEta < mEtaEdges.front() || absEta >= mEtaEdges.back() ||
        rigidity < mRigEdges.front() || rigidity >= mRigEdges.back()) {
      return -1.f;
    }
    const int nRig = mRigEdges.size() - 1;
    const int iEta = std::upper_bound(mEtaEdges.begin(), mEtaEdges.end(), absEta) - mEtaEdges.begin() - 1;
    const int iRig = std::upper_bound(mRigEdges.begin(), mRigEdges.end(), rigidity) - mRigEdges.begin() - 1;
    const int offset = it->second[iEta * nRig + iRig];
    if (offset < 0) {
      return -1.f;
    }
    const float* cdf = mCdf.at(pdg).data();
    float dedx = inverseCdf(cdf + offset, u);

    // neighbouring bin on the side of the track, no interpolation if it is empty
    const double logRig = std::log(rigidity);
    const int iNext = logRig < mLogRigCentres[iRig] ? iRig - 1 : iRig + 1;
    if (iNext >= 0 && iNext < nRig && it->second[iEta * nRig + iNext] >= 0) {
      const float w = (logRig - mLogRigCentres[iRig]) / (mLogRigCentres[iNext] - mLogRigCentres[iRig]);
      dedx += w * (inverseCdf(cdf + it->second[iEta * nRig + iNext], u) - dedx);
    }
    return dedx;
  }

 private:
  /// dE/dx at which the cumulative distribution reaches u, uniform within the dE/dx bin
  float inverseCdf(const float* cdf, double u) const
  {
    const int nDedx = mDedxEdges.size() - 1;
    const int k = std::min(static_cast<int>(std::upper_bound(cdf, cdf + nDedx, u) - cdf), nDedx - 1);
    const double lo = k == 0 ? 0. : cdf[k - 1];
    const double frac = cdf[k] > lo ? (u - lo) / (cdf[k] - lo) : 0.5;
    return mDedxEdges[k] + frac * (mDedxEdges[k + 1] - mDedxEdges[k]);
  }

  static std::vector<double> binEdges(const TH1* h)
  {
    std::vector<double> edges;
    for (int i = 1; i <= h->GetNbinsX() + 1; ++i) {
      edges.push_back(h->GetXaxis()->GetBinLowEdge(i));
    }
    return edges;
  }

  std::vector<double> mEtaEdges;
  std::vector<double> mRigEdges;
  std::vector<double> mLogRigCentres;
  std::vector<double> mDedxEdges;
  std::unordered_map<int, std::vector<int>> mOffset; // |pdg| -> start of the cdf of bin etaBin * nRig + rigBin, -1 if empty
  std::unordered_map<int, std::vector<float>> mCdf;  // |pdg| -> cumulative dE/dx distributions of the filled bins
};

struct OnTheFlyTrackerPid {
  static constexpr int kNumSpecies = 9;
  static constexpr std::array<int, kNumSpecies> kPdgCodes = {11, 13, 211, 321, 2212, 1000010020, 1000010030, 1000020030, 1000020040};
  static constexpr std::array<const char*, kNumSpecies> kNames = {"Elec", "Muon", "Pion", "Kaon", "Prot", "Deut", "Trit", "He3", "Al"};
  static constexpr std::array<const char*, kNumSpecies> kPrettyNames = {"#it{e}", "#it{#mu}", "#it{#pi}", "#it{K}", "#it{p}", "#it{d}", "#it{t}", "#it{^{3}He}", "#it{^{4}He}"};

  std::array<std::shared_ptr<TH2>, kNumSpecies> h2dToTvsPperParticle;
  std::array<std::shared_ptr<TH2>, kNumSpecies> h2dToTvsPperParticleLin;
  std::array<std::shared_ptr<TH3>, kNumSpecies> h3dToTvsPEtaPerParticle;
  std::array<std::shared_ptr<TH2>, kNumSpecies> h2dTracks;
  std::array<std::shared_ptr<TH2>, kNumSpecies> h2dTracksWithDedx;

  Service<o2::framework::O2DatabasePDG> pdgDatabase{};
  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  TrackDedxLut mTrackLut;
  TRandom3 mRandom;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<std::string> trackLut{"trackLut", "", "Per-track dE/dx LUT from MakeToTLUT.C (local file or ccdb:<path>)"};
  Configurable<int> seed{"seed", 0, "Seed for the dE/dx sampling (0: different in every job)"};
  Configurable<bool> primariesOnly{"primariesOnly", false, "Only physical primaries (otherwise all tracks, including decay products)"};
  Configurable<bool> useTrueKinematics{"useTrueKinematics", false, "Fill the histograms vs the true instead of the reconstructed p/z and |eta| (the dE/dx is always sampled at the true ones)"};

  // dE/dx axis, logarithmic
  Configurable<int> numDedxBins{"numDedxBins", 300, "Number of logarithmic dE/dx bins"};
  Configurable<float> dedxMin{"dedxMin", 1e2f, "Minimum of the dE/dx axis in e-"};
  Configurable<float> dedxMax{"dedxMax", 1e5f, "Maximum of the dE/dx axis in e-"};

  ConfigurableAxis axisEta{"axisEta", {VARIABLE_WIDTH, 0.0f, 0.5f, 1.0f, 1.5f, 2.0f, 2.5f, 3.0f, 4.0f}, "|#eta| binning of the eta-differential histograms"};

  // Momentum axis binning choice
  Configurable<int> numLogBins{"numLogBins", 200, "Number of logarithmic momentum bins"};
  Configurable<bool> enableLogPBins{"enableLogPBins", true, "Enable logarithmic momentum binning for dE/dx vs p histograms"};
  Configurable<bool> enableLinearPBins{"enableLinearPBins", false, "Enable linear momentum binning for dE/dx vs p histograms"};
  Configurable<int> numLinearBins{"numLinearBins", 500, "Number of linear momentum bins"};
  Configurable<float> linearPMin{"linearPMin", 0.0f, "Minimum momentum for linear-binned histograms (GeV/c)"};
  Configurable<float> linearPMax{"linearPMax", 10.0f, "Maximum momentum for linear-binned histograms (GeV/c)"};

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setTimestamp(-1);

    mRandom.SetSeed(seed.value);
    if (!mTrackLut.load(trackLut.value, ccdb.operator->(), kPdgCodes)) {
      LOG(warning) << "dE/dx LUT incomplete, species without LUT get no tracker PID.";
    }
    bookHistograms();
  }

  static std::vector<double> logBinning(int nBins, double min, double max)
  {
    std::vector<double> edges;
    const double logMin = std::log10(min);
    const double dLog = (std::log10(max) - logMin) / nBins;
    for (int i = 0; i <= nBins; ++i) {
      edges.push_back(std::pow(10, logMin + i * dLog));
    }
    return edges;
  }

  void bookHistograms()
  {
    if (!enableLogPBins.value && !enableLinearPBins.value) {
      LOG(fatal) << "At least one of enableLogPBins or enableLinearPBins must be true.";
    }

    const std::vector<double> logBins = logBinning(numLogBins.value, kMinMomentumForLogBins, kMaxMomentumForLogBins);
    const std::string signalTitle = "d#it{E}/d#it{x} (e^{-})";
    const AxisSpec axisDedx{logBinning(numDedxBins.value, dedxMin.value, dedxMax.value), signalTitle};
    const AxisSpec axisMomentumLog{logBins, "#it{p/z} (GeV/#it{c})"};
    const AxisSpec axisMomentumLin{numLinearBins.value, static_cast<double>(linearPMin.value), static_cast<double>(linearPMax.value), "#it{p/z} (GeV/#it{c})"};
    const AxisSpec axisAbsEta{axisEta, "|#eta|"};

    if (enableLogPBins.value) {
      histos.add("hToTvsP", ("Signal vs #it{p/z}; #it{p/z} (GeV/#it{c}); " + signalTitle).c_str(), kTH2F, {axisMomentumLog, axisDedx});
      histos.add("hToTvsMom", ("Signal vs #it{p}; #it{p} (GeV/#it{c}); " + signalTitle).c_str(), kTH2F, {logBins, axisDedx});
    }
    if (enableLinearPBins.value) {
      histos.add("hToTvsPLin", ("Signal vs #it{p/z} (linear); #it{p/z} (GeV/#it{c}); " + signalTitle).c_str(), kTH2F, {axisMomentumLin, axisDedx});
      histos.add("hToTvsMomLin", ("Signal vs #it{p} (linear); #it{p} (GeV/#it{c}); " + signalTitle).c_str(), kTH2F, {axisMomentumLin, axisDedx});
    }

    for (int i = 0; i < kNumSpecies; ++i) {
      if (enableLogPBins.value) {
        const std::string name = std::format("ToTvsP/hToTvsP_{}", kNames[i]);
        const std::string title = std::format("Signal vs #it{{p/z}} for {}; #it{{p/z}} (GeV/#it{{c}}); {}", kPrettyNames[i], signalTitle);
        h2dToTvsPperParticle[i] = histos.add<TH2>(name.c_str(), title.c_str(), kTH2F, {axisMomentumLog, axisDedx});
      }
      if (enableLinearPBins.value) {
        const std::string name = std::format("ToTvsPLin/hToTvsPLin_{}", kNames[i]);
        const std::string title = std::format("Signal vs #it{{p/z}} for {} (linear); #it{{p/z}} (GeV/#it{{c}}); {}", kPrettyNames[i], signalTitle);
        h2dToTvsPperParticleLin[i] = histos.add<TH2>(name.c_str(), title.c_str(), kTH2F, {axisMomentumLin, axisDedx});
      }

      const std::string nameEta = std::format("ToTvsPEta/hToTvsPEta_{}", kNames[i]);
      const std::string titleEta = std::format("Signal vs #it{{p/z}} and |#eta| for {}; #it{{p/z}} (GeV/#it{{c}}); |#eta|; {}", kPrettyNames[i], signalTitle);
      h3dToTvsPEtaPerParticle[i] = histos.add<TH3>(nameEta.c_str(), titleEta.c_str(), kTH3F, {axisMomentumLog, axisAbsEta, axisDedx});

      const std::string nameTracks = std::format("Bookkeeping/hTracks_{}", kNames[i]);
      const std::string titleTracks = std::format("Reconstructed {}; #it{{p/z}} (GeV/#it{{c}}); |#eta|", kPrettyNames[i]);
      h2dTracks[i] = histos.add<TH2>(nameTracks.c_str(), titleTracks.c_str(), kTH2F, {axisMomentumLog, axisAbsEta});

      const std::string nameWithDedx = std::format("Bookkeeping/hTracksWithDedx_{}", kNames[i]);
      const std::string titleWithDedx = std::format("Reconstructed {} with dE/dx; #it{{p/z}} (GeV/#it{{c}}); |#eta|", kPrettyNames[i]);
      h2dTracksWithDedx[i] = histos.add<TH2>(nameWithDedx.c_str(), titleWithDedx.c_str(), kTH2F, {axisMomentumLog, axisAbsEta});
    }
  }

  void process(soa::Join<aod::Tracks, aod::McTrackLabels> const& tracks,
               aod::McParticles const&)
  {
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      const auto& mcParticle = track.mcParticle();
      if (primariesOnly.value && !mcParticle.isPhysicalPrimary()) {
        continue;
      }
      const int pdg = std::abs(mcParticle.pdgCode());
      const auto itSpecies = std::find(kPdgCodes.begin(), kPdgCodes.end(), pdg);
      if (itSpecies == kPdgCodes.end()) {
        continue;
      }
      const auto* pdgInfo = pdgDatabase->GetParticle(mcParticle.pdgCode());
      if (!pdgInfo || pdgInfo->Charge() == 0) {
        continue;
      }
      const int iSpecies = std::distance(kPdgCodes.begin(), itSpecies);
      const float absCharge = std::abs(pdgInfo->Charge()) / 3.f;
      const float rigidityTrue = mcParticle.p() / absCharge;

      // the OTF tracks carry p/z
      const float rigidity = useTrueKinematics.value ? rigidityTrue : track.p();
      const float absEta = std::abs(useTrueKinematics.value ? mcParticle.eta() : track.eta());
      const float p = rigidity * absCharge;
      h2dTracks[iSpecies]->Fill(rigidity, absEta);

      const float dedx = mTrackLut.sample(pdg, mcParticle.eta(), rigidityTrue, mRandom.Rndm());
      if (dedx < 0.f) {
        continue;
      }
      h2dTracksWithDedx[iSpecies]->Fill(rigidity, absEta);
      h3dToTvsPEtaPerParticle[iSpecies]->Fill(rigidity, absEta, dedx);
      if (enableLogPBins.value) {
        histos.fill(HIST("hToTvsMom"), p, dedx);
        histos.fill(HIST("hToTvsP"), rigidity, dedx);
        h2dToTvsPperParticle[iSpecies]->Fill(rigidity, dedx);
      }
      if (enableLinearPBins.value) {
        histos.fill(HIST("hToTvsMomLin"), p, dedx);
        histos.fill(HIST("hToTvsPLin"), rigidity, dedx);
        h2dToTvsPperParticleLin[iSpecies]->Fill(rigidity, dedx);
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyTrackerPid>(cfgc)};
}
