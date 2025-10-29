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
/// \brief This task produces the PID information that can be obtained from the tracker layers (i.e. ToT and possibly cluster size).
///        So far only ToT implemented. It currently contemplates 5 (9) particle types: electrons, muons, pions, kaons and
///        protons (as well as deuterons, tritons, helium-3 and alphas/helium-4 if added via the event generator).
///
/// \author Henrik Fribert TUM
/// \author Nicolò Jacazio Università del Piemonte Orientale
/// \since  May 22, 2025
///

#include "TableHelper.h"

#include "ALICE3/Core/DelphesO2TrackSmearer.h"
#include "ALICE3/Core/TrackUtilities.h"
#include "ALICE3/DataModel/OTFPIDTrk.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/GeomConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <CommonUtils/NameConf.h>
#include <DataFormatsCalibration/MeanVertexObject.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <ReconstructionDataFormats/HelixHelper.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>

#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TVector3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

static constexpr std::array<float, 11> kTrackerRadii = {0.5f, 1.2f, 2.5f, 3.75f, 7.0f, 12.0f, 20.0f, 30.0f, 45.0f, 60.0f, 80.0f};

// Constants for magic numbers
static constexpr int kMinEntriesForProjection = 10;
static constexpr double kMinMomentumForLogBins = 0.05;
static constexpr double kMaxMomentumForLogBins = 10.0;
static constexpr int kMinLayerForTruncation = 3;

// Constants for validation conditions
static constexpr size_t kMinValidHits = 1;
static constexpr size_t kLowValidHits = 2;
static constexpr size_t kMidLowValidHits = 3;
static constexpr size_t kMidHighValidHits = 4;
static constexpr size_t kHighValidHits1 = 5;
static constexpr size_t kHighValidHits2 = 6;
static constexpr size_t kMaxValidHits = 7;

// Constants for truncated mean calculation
static constexpr size_t kMaxValidHitsForTruncation7Plus = 4;
static constexpr size_t kMaxValidHitsForTruncation56 = 3;
static constexpr size_t kMaxValidHitsForTruncation34 = 2;
static constexpr size_t kMaxValidHitsForTruncation12 = 1;

// Constants for LUT binning
// To do: Include in LUT header or similar
static constexpr int kLUTEtaBins = 50;
static constexpr float kLUTEtaMin = -2.5f;
static constexpr float kLUTEtaMax = 2.5f;
static constexpr int kLUTPtBins = 500;
static constexpr float kLUTPtMin = 0.0f;
static constexpr float kLUTPtMax = 10.0f;

class ToTLUT
{
 public:
  explicit ToTLUT(int maxLayers, float analysisEtaMinVal = 0.0f, float analysisEtaMaxVal = 1.0f, float analysisPtMinVal = 0.0f, float analysisPtMaxVal = 10.0f)
    : mMaxLayers(maxLayers),
      mAnalysisEtaMin(analysisEtaMinVal),
      mAnalysisEtaMax(analysisEtaMaxVal),
      mAnalysisPtMin(analysisPtMinVal),
      mAnalysisPtMax(analysisPtMaxVal),
      mEtaBins(kLUTEtaBins),
      mEtaMin(kLUTEtaMin),
      mEtaMax(kLUTEtaMax),
      mPtBins(kLUTPtBins),
      mPtMin(kLUTPtMin),
      mPtMax(kLUTPtMax),
      mEtaBinWidth((kLUTEtaMax - kLUTEtaMin) / kLUTEtaBins),
      mPtBinWidth((kLUTPtMax - kLUTPtMin) / kLUTPtBins)
  {
    mPdgToIndexMap.reserve(10);
    mIndexToPdgMap.reserve(10);
  }
  ToTLUT() = delete;

  ~ToTLUT()
  {
    for (const auto& hist_ptr : mLUTHistogramFlat) {
      if (hist_ptr) {
        delete hist_ptr;
      }
    }
  }

  void setCcdbManager(o2::ccdb::BasicCCDBManager* mgr) { mCcdbManager = mgr; }

  bool load(int pdg, const std::string& filename)
  {
    if (!filename.empty() && strncmp(filename.c_str(), "ccdb:", 5) == 0) {
      std::string basePath = std::string(filename).substr(5);
      std::string path = basePath + "/PDG_" + std::to_string(pdg);
      const std::string outPath = "/tmp/ToTLUTs/";

      std::string localFilename = Form("%s/lut_tot_%d.root", outPath.c_str(), pdg);
      std::ifstream checkFile(localFilename);
      if (!checkFile.is_open()) {
        if (!mCcdbManager) {
          LOG(fatal) << "CCDB manager not set. Please set it before loading LUT from CCDB.";
        }
        std::map<std::string, std::string> metadata;
        mCcdbManager->getCCDBAccessor().retrieveBlob(path, outPath, metadata, 1);

        std::string foundFile = Form("%s/%s/snapshot.root", outPath.c_str(), path.c_str());
        std::ifstream testFile(foundFile);
        if (!testFile.is_open()) {
          LOG(error) << "Could not find downloaded CCDB file for PDG " << pdg;
          return false;
        }
        testFile.close();

        return load(pdg, foundFile);
      } else {
        checkFile.close();
        return load(pdg, localFilename);
      }
    }

    TFile* f = TFile::Open(filename.c_str());
    if (!f || f->IsZombie()) {
      LOG(error) << "Failed to open LUT file: " << filename;
      return false;
    }

    int currentPdgIdx;
    auto it = mPdgToIndexMap.find(pdg);
    if (it == mPdgToIndexMap.end()) {
      currentPdgIdx = mIndexToPdgMap.size();
      mPdgToIndexMap[pdg] = currentPdgIdx;
      mIndexToPdgMap.push_back(pdg);

      size_t totalSize = (currentPdgIdx + 1) * mMaxLayers * mEtaBins * mPtBins;
      if (mLUTHistogramFlat.size() < totalSize) {
        mLUTHistogramFlat.resize(totalSize, nullptr);
      }
    } else {
      currentPdgIdx = it->second;
    }

    bool success = true;
    for (int layer = 0; layer < mMaxLayers; ++layer) {
      for (int etaBin = 0; etaBin < mEtaBins; ++etaBin) {
        float etaMinBin = mEtaMin + etaBin * mEtaBinWidth;
        float etaMaxBin = etaMinBin + mEtaBinWidth;

        float etaCenter = (etaMinBin + etaMaxBin) / 2.0f;
        if (std::abs(etaCenter) < mAnalysisEtaMin || std::abs(etaCenter) > mAnalysisEtaMax) {
          continue;
        }

        for (int ptBin = 0; ptBin < mPtBins; ++ptBin) {
          float ptMinBin = mPtMin + ptBin * mPtBinWidth;
          float ptMaxBin = ptMinBin + mPtBinWidth;
          float ptCenter = (ptMinBin + ptMaxBin) / 2.0f;

          if (ptCenter < mAnalysisPtMin || ptCenter >= mAnalysisPtMax) {
            continue;
          }

          TString histName = Form("tot_%d_barrel%d_eta%.2f-%.2f_pt%.2f-%.2f", pdg, layer, etaMinBin, etaMaxBin, ptMinBin, ptMaxBin);

          TH1F* histFromFile = dynamic_cast<TH1F*>(f->Get(histName));
          if (histFromFile) {
            TH1F* clonedHist = static_cast<TH1F*>(histFromFile->Clone());
            clonedHist->SetDirectory(nullptr);

            size_t flatIdx = getFlatIndex(currentPdgIdx, layer, etaBin, ptBin);
            mLUTHistogramFlat[flatIdx] = clonedHist;
          } else {
            size_t flatIdx = getFlatIndex(currentPdgIdx, layer, etaBin, ptBin);
            mLUTHistogramFlat[flatIdx] = nullptr;
            success = false;
          }
        }
      }
    }

    f->Close();
    delete f;
    return success;
  }

  TH1F* getHistogramForSampling(int pdgIdx, int layer, int etaBin, int ptBin) const
  {
    if (pdgIdx < 0 || static_cast<size_t>(pdgIdx) >= getNumPdgTypes() ||
        layer < 0 || layer >= mMaxLayers ||
        etaBin < 0 || etaBin >= mEtaBins || ptBin < 0 || ptBin >= mPtBins) {
      return nullptr;
    }
    size_t flatIdx = getFlatIndex(pdgIdx, layer, etaBin, ptBin);
    return (flatIdx < mLUTHistogramFlat.size()) ? mLUTHistogramFlat[flatIdx] : nullptr;
  }

  int getPdgIndex(int pdgCode) const
  {
    auto it = mPdgToIndexMap.find(pdgCode);
    if (it != mPdgToIndexMap.end()) {
      return it->second;
    }
    return -1;
  }

  inline int getEtaBin(float eta) const
  {
    const float clampedEta = std::max(mEtaMin, std::min(eta, mEtaMax - 1e-6f));
    return std::min(static_cast<int>((clampedEta - mEtaMin) / mEtaBinWidth), mEtaBins - 1);
  }

  inline int getPtBin(float pt) const
  {
    const float clampedPt = std::max(mPtMin, std::min(pt, mPtMax - 1e-6f));
    return std::min(static_cast<int>((clampedPt - mPtMin) / mPtBinWidth), mPtBins - 1);
  }

  inline size_t getFlatIndex(int pdgIdx, int layer, int etaBin, int ptBin) const
  {
    return ((pdgIdx * mMaxLayers + layer) * mEtaBins + etaBin) * mPtBins + ptBin;
  }

  size_t getNumPdgTypes() const
  {
    return mIndexToPdgMap.size();
  }

  std::vector<TH1F*> mLUTHistogramFlat;

  std::unordered_map<int, int> mPdgToIndexMap;
  std::vector<int> mIndexToPdgMap;

  int mMaxLayers;
  float mAnalysisEtaMin;
  float mAnalysisEtaMax;
  float mAnalysisPtMin;
  float mAnalysisPtMax;
  int mEtaBins;
  float mEtaMin;
  float mEtaMax;
  int mPtBins;
  float mPtMin;
  float mPtMax;

  float mEtaBinWidth;
  float mPtBinWidth;

 private:
  o2::ccdb::BasicCCDBManager* mCcdbManager = nullptr;
};

static constexpr int kNumHypothesisParticles = 9;
std::array<std::array<std::shared_ptr<TH2>, kNumHypothesisParticles>, kNumHypothesisParticles> h2dBarrelNsigmaTrue;
std::array<std::shared_ptr<TH2>, kNumHypothesisParticles> h2dHitsPerTrackVsP;
std::array<std::shared_ptr<TH2>, kNumHypothesisParticles> h2dToTvsPperParticle;

struct OnTheFlyTrackerPid {

  float calculateNsigma(float measuredToT, float expectedToT, float resolution)
  {
    if (resolution <= 0)
      return 999.f;
    return (measuredToT - expectedToT) / resolution;
  }

  float getToTMeanFromMomentumSlice(std::shared_ptr<TH2> hist, float momentum)
  {
    if (!hist)
      return -1.f;
    int binX = hist->GetXaxis()->FindBin(momentum);
    TH1D* proj = hist->ProjectionY("temp", binX, binX);
    if (proj->GetEntries() < kMinEntriesForProjection) {
      delete proj;
      return -1.f;
    }
    float mean = proj->GetMean();
    delete proj;
    return mean;
  }

  float getToTResolutionFromMomentumSlice(std::shared_ptr<TH2> hist, float momentum)
  {
    if (!hist)
      return -1.f;
    int binX = hist->GetXaxis()->FindBin(momentum);
    TH1D* proj = hist->ProjectionY("temp", binX, binX);
    if (proj->GetEntries() < kMinEntriesForProjection) {
      delete proj;
      return -1.f;
    }
    float stddev = proj->GetStdDev();
    delete proj;
    return stddev;
  }

  float computeTrackLength(o2::track::TrackParCov track, float radius, float magneticField)
  {
    float length = -100;
    o2::math_utils::CircleXYf_t trcCircle;
    float sna, csa;
    track.getCircleParams(magneticField, trcCircle, sna, csa);

    const float centerDistance = std::hypot(trcCircle.xC, trcCircle.yC);

    if (centerDistance < trcCircle.rC + radius && centerDistance > std::fabs(trcCircle.rC - radius)) {
      length = 0.0f;
      const float ux = trcCircle.xC / centerDistance;
      const float uy = trcCircle.yC / centerDistance;
      const float vx = -uy;
      const float vy = +ux;
      const float radical = (centerDistance * centerDistance - trcCircle.rC * trcCircle.rC + radius * radius) / (2.0f * centerDistance);
      const float displace = (0.5f / centerDistance) * std::sqrt(
                                                         (-centerDistance + trcCircle.rC - radius) *
                                                         (-centerDistance - trcCircle.rC + radius) *
                                                         (-centerDistance + trcCircle.rC + radius) *
                                                         (centerDistance + trcCircle.rC + radius));

      const float point1[2] = {radical * ux + displace * vx, radical * uy + displace * vy};
      const float point2[2] = {radical * ux - displace * vx, radical * uy - displace * vy};

      std::array<float, 3> mom;
      track.getPxPyPzGlo(mom);
      const float scalarProduct1 = point1[0] * mom[0] + point1[1] * mom[1];
      const float scalarProduct2 = point2[0] * mom[0] + point2[1] * mom[1];

      std::array<float, 3> startPoint;
      track.getXYZGlo(startPoint);

      float cosAngle = -1000, modulus = -1000;

      if (scalarProduct1 > scalarProduct2) {
        modulus = std::hypot(point1[0] - trcCircle.xC, point1[1] - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
        cosAngle = (point1[0] - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (point1[1] - trcCircle.yC) * (startPoint[1] - trcCircle.yC);
      } else {
        modulus = std::hypot(point2[0] - trcCircle.xC, point2[1] - trcCircle.yC) * std::hypot(startPoint[0] - trcCircle.xC, startPoint[1] - trcCircle.yC);
        cosAngle = (point2[0] - trcCircle.xC) * (startPoint[0] - trcCircle.xC) + (point2[1] - trcCircle.yC) * (startPoint[1] - trcCircle.yC);
      }
      cosAngle /= modulus;
      length = trcCircle.rC * std::acos(cosAngle);
      length *= std::sqrt(1.0f + track.getTgl() * track.getTgl());
    }
    return length;
  }

  Produces<aod::UpgradeTrkPidSignals> tableUpgradeTrkPidSignals;
  Produces<aod::UpgradeTrkPids> tableUpgradeTrkPids;

  Service<o2::framework::O2DatabasePDG> pdg;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  std::unique_ptr<ToTLUT> mToTLUT;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<std::string> lutTotEl{"lutTotEl", "ccdb:Users/h/hfribert/ToT_LUTs", "ToT LUT for electrons"};
  Configurable<std::string> lutTotMu{"lutTotMu", "ccdb:Users/h/hfribert/ToT_LUTs", "ToT LUT for muons"};
  Configurable<std::string> lutTotPi{"lutTotPi", "ccdb:Users/h/hfribert/ToT_LUTs", "ToT LUT for pions"};
  Configurable<std::string> lutTotKa{"lutTotKa", "ccdb:Users/h/hfribert/ToT_LUTs", "ToT LUT for kaons"};
  Configurable<std::string> lutTotPr{"lutTotPr", "ccdb:Users/h/hfribert/ToT_LUTs", "ToT LUT for protons"};
  Configurable<std::string> lutTotDe{"lutTotDe", "ccdb:Users/h/hfribert/ToT_LUTs", "ToT LUT for deuteron"};
  Configurable<std::string> lutTotTr{"lutTotTr", "ccdb:Users/h/hfribert/ToT_LUTs", "ToT LUT for triton"};
  Configurable<std::string> lutTotHe{"lutTotHe", "ccdb:Users/h/hfribert/ToT_LUTs", "ToT LUT for helium-3"};
  Configurable<std::string> lutTotAl{"lutTotAl", "ccdb:Users/h/hfribert/ToT_LUTs", "ToT LUT for alphas"};

  Configurable<float> dBz{"dBz", 20, "magnetic field (kilogauss) for track propagation"};
  Configurable<int> maxBarrelLayers{"maxBarrelLayers", 11, "Maximum number of barrel layers"};
  Configurable<int> numLogBins{"numLogBins", 200, "Number of logarithmic momentum bins"};
  Configurable<float> analysisEtaMin{"analysisEtaMin", 0.0f, "Minimum |eta| for LUT loading optimization"};
  Configurable<float> analysisEtaMax{"analysisEtaMax", 1.0f, "Maximum |eta| for LUT loading optimization"};
  Configurable<float> analysisPtMin{"analysisPtMin", 0.0f, "Minimum pT (GeV/c) for LUT loading optimization"};
  Configurable<float> analysisPtMax{"analysisPtMax", 10.0f, "Maximum pT (GeV/c) for LUT loading optimization"};

  std::vector<double> mLogBins;

  std::array<int, kNumHypothesisParticles> mHypothesisPdgCodes = {
    11,         // Electron
    13,         // Muon
    211,        // Pion
    321,        // Kaon
    2212,       // Proton
    1000010020, // Deuteron
    1000010030, // Triton
    1000020030, // Helium-3
    1000020040  // Alpha
  };

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setTimestamp(-1);

    if (static_cast<size_t>(maxBarrelLayers.value) > kTrackerRadii.size()) {
      LOG(fatal) << "Configured maxBarrelLayers (" << maxBarrelLayers.value
                 << ") exceeds the size of kTrackerRadii (" << kTrackerRadii.size()
                 << "). Please adjust maxBarrelLayers.";
    }

    mToTLUT = std::make_unique<ToTLUT>(maxBarrelLayers.value, analysisEtaMin.value, analysisEtaMax.value, analysisPtMin.value, analysisPtMax.value);

    mToTLUT->setCcdbManager(ccdb.operator->());

    bool loaded = true;
    loaded &= mToTLUT->load(11, lutTotEl.value);
    loaded &= mToTLUT->load(13, lutTotMu.value);
    loaded &= mToTLUT->load(211, lutTotPi.value);
    loaded &= mToTLUT->load(321, lutTotKa.value);
    loaded &= mToTLUT->load(2212, lutTotPr.value);
    loaded &= mToTLUT->load(1000010020, lutTotDe.value);
    loaded &= mToTLUT->load(1000010030, lutTotTr.value);
    loaded &= mToTLUT->load(1000020030, lutTotHe.value);
    loaded &= mToTLUT->load(1000020040, lutTotAl.value);

    if (!loaded) {
      LOG(warning) << "Failed to load one or more ToT LUTs. PID results might be incomplete.";
    }

    // Logarithmic momentum bins
    mLogBins.clear();
    double pMin = kMinMomentumForLogBins;
    double pMax = kMaxMomentumForLogBins;
    double logMin = std::log10(pMin);
    double logMax = std::log10(pMax);
    double dLog = (logMax - logMin) / numLogBins.value;
    for (int i = 0; i <= numLogBins.value; ++i) {
      mLogBins.push_back(std::pow(10, logMin + i * dLog));
    }

    const AxisSpec axisMomentum{mLogBins, "#it{p/z} (GeV/#it{c})"};
    const AxisSpec axisToT{600, 0., 300., "ToT (#mus/10#mum)"};
    const AxisSpec axisNsigma{200, -10., 10., "N#sigma"};
    const AxisSpec axisLayer{maxBarrelLayers.value, -0.5, static_cast<double>(maxBarrelLayers.value) - 0.5, "Layer"};
    const AxisSpec axisHitsPerTrack{maxBarrelLayers.value + 1, -0.5, static_cast<double>(maxBarrelLayers.value) + 0.5, "# Hits per track"};

    histos.add("hToTvsP", "ToT vs #it{p/z}; #it{p/z} (GeV/#it{c}); ToT (#mus/10#mum)", kTH2F, {axisMomentum, axisToT});
    histos.add("hToTvsPt", "ToT vs #it{p}; #it{p} (GeV/#it{c}); ToT (#mus/10#mum)", kTH2F, {mLogBins, axisToT});
    histos.add("hHitLayers", "Number of hits on each detector layer;Layer;Counts", kTH1F, {axisLayer});
    histos.add("hHitMultiplicity", "Hit multiplicity along the track; # Hits per track;Counts", kTH1F, {axisHitsPerTrack});

    std::vector<std::pair<int, std::string>> particleInfo = {
      {11, "Elec"}, {13, "Muon"}, {211, "Pion"}, {321, "Kaon"}, {2212, "Prot"}, {1000010020, "Deut"}, {1000010030, "Trit"}, {1000020030, "He3"}, {1000020040, "Al"}};

    for (size_t iTrue = 0; iTrue < particleInfo.size(); ++iTrue) {
      std::string trueName = particleInfo[iTrue].second;
      std::string trueNamePretty = trueName; // Fallback
      if (trueName == "Elec")
        trueNamePretty = "#it{e}";
      else if (trueName == "Muon")
        trueNamePretty = "#it{#mu}";
      else if (trueName == "Pion")
        trueNamePretty = "#it{#pi}";
      else if (trueName == "Kaon")
        trueNamePretty = "#it{K}";
      else if (trueName == "Prot")
        trueNamePretty = "#it{p}";
      else if (trueName == "Deut")
        trueNamePretty = "#it{d}";
      else if (trueName == "Trit")
        trueNamePretty = "#it{t}";
      else if (trueName == "He3")
        trueNamePretty = "#it{^{3}He}";
      else if (trueName == "Al")
        trueNamePretty = "#it{^{4}He}";

      std::string hitsVsPName = "HitsPerTrack/hHitsPerTrackVsP_" + trueName;
      std::string hitsVsPTitle = "N_hits vs #it{p/z} for " + trueNamePretty + "; #it{p/z} (GeV/#it{c}); N_hits";
      h2dHitsPerTrackVsP[iTrue] = histos.add<TH2>(hitsVsPName.c_str(), hitsVsPTitle.c_str(), kTH2F, {axisMomentum, axisHitsPerTrack});

      std::string totVsPName = "ToTvsP/hToTvsP_" + trueName;
      std::string totVsPTitle = "ToT vs #it{p/z} for " + trueNamePretty + "; #it{p/z} (GeV/#it{c}); ToT (#mus/10#mum)";
      h2dToTvsPperParticle[iTrue] = histos.add<TH2>(totVsPName.c_str(), totVsPTitle.c_str(), kTH2F, {axisMomentum, axisToT});

      for (size_t iHyp = 0; iHyp < particleInfo.size(); ++iHyp) {
        std::string hypName = particleInfo[iHyp].second;
        std::string hypNamePretty = hypName; // Fallback
        if (hypName == "Elec")
          hypNamePretty = "#it{e}";
        else if (hypName == "Muon")
          hypNamePretty = "#it{#mu}";
        else if (hypName == "Pion")
          hypNamePretty = "#it{#pi}";
        else if (hypName == "Kaon")
          hypNamePretty = "#it{K}";
        else if (hypName == "Prot")
          hypNamePretty = "#it{p}";
        else if (hypName == "Deut")
          hypNamePretty = "#it{d}";
        else if (hypName == "Trit")
          hypNamePretty = "#it{t}";
        else if (hypName == "He3")
          hypNamePretty = "#it{^{3}He}";
        else if (hypName == "Al")
          hypNamePretty = "#it{^{4}He}";

        std::string histName = "NSigma/BarrelNsigmaTrue" + trueName + "Vs" + hypName + "Hypothesis";
        std::string histTitle = "Nsigma (True " + trueNamePretty + " vs Hyp " + hypNamePretty + "); #it{p/z} (GeV/#it{c}); N#sigma";
        h2dBarrelNsigmaTrue[iTrue][iHyp] = histos.add<TH2>(histName.c_str(), histTitle.c_str(), kTH2F, {axisMomentum, axisNsigma});
      }
    }
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision,
               soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels> const& tracks,
               aod::McParticles const& /*mcParticles*/,
               aod::McCollisions const& /*mcCollisions*/)
  {
    o2::dataformats::VertexBase mcPvVtx({0.0f, 0.0f, 0.0f}, {0.});

    if (collision.has_mcCollision()) {
      const auto& mcCollisionObject = collision.mcCollision();
      mcPvVtx.setX(mcCollisionObject.posX());
      mcPvVtx.setY(mcCollisionObject.posY());
      mcPvVtx.setZ(mcCollisionObject.posZ());
    }

    for (const auto& track : tracks) {
      float truncatedMeanToT = -1.0f;
      std::array<float, kNumHypothesisParticles> nSigmaValues;
      nSigmaValues.fill(999.f);

      if (!track.has_mcParticle()) {
        tableUpgradeTrkPidSignals(truncatedMeanToT);
        tableUpgradeTrkPids(nSigmaValues[0], nSigmaValues[1], nSigmaValues[2], nSigmaValues[3],
                            nSigmaValues[4], nSigmaValues[5], nSigmaValues[6], nSigmaValues[7], nSigmaValues[8]);
        continue;
      }

      const auto& mcParticle = track.mcParticle();

      const auto& pdgInfo = pdg->GetParticle(mcParticle.pdgCode());
      if (!pdgInfo) {
        tableUpgradeTrkPidSignals(truncatedMeanToT);
        tableUpgradeTrkPids(nSigmaValues[0], nSigmaValues[1], nSigmaValues[2], nSigmaValues[3],
                            nSigmaValues[4], nSigmaValues[5], nSigmaValues[6], nSigmaValues[7], nSigmaValues[8]);
        continue;
      }

      const float pt = mcParticle.pt();
      const float p = mcParticle.p();
      const float eta = mcParticle.eta();
      const int truePdgCode = std::abs(mcParticle.pdgCode());

      float rigidity = p;
      if (pdgInfo) {
        const float trueCharge = std::abs(pdgInfo->Charge()) / 3.0f;
        if (trueCharge > 0) {
          rigidity = p / trueCharge;
        }
      }

      int truePdgIdx = mToTLUT->getPdgIndex(truePdgCode);
      if (truePdgIdx == -1) {
        tableUpgradeTrkPidSignals(truncatedMeanToT);
        tableUpgradeTrkPids(nSigmaValues[0], nSigmaValues[1], nSigmaValues[2], nSigmaValues[3],
                            nSigmaValues[4], nSigmaValues[5], nSigmaValues[6], nSigmaValues[7], nSigmaValues[8]);
        continue;
      }

      const int binnedPt = mToTLUT->getPtBin(pt);
      const int binnedEta = mToTLUT->getEtaBin(std::abs(eta));

      uint16_t hitMap = 0;
      int nHitLayers = 0;
      o2::track::TrackParCov o2track = o2::upgrade::convertMCParticleToO2Track(mcParticle, pdg);

      float xPv = -100.f;
      static constexpr float kTrkXThreshold = -99.f;
      if (o2track.propagateToDCA(mcPvVtx, dBz)) {
        xPv = o2track.getX();
      }

      if (xPv > kTrkXThreshold) {
        for (int layer = 0; layer < maxBarrelLayers.value; ++layer) {
          float layerRadius = kTrackerRadii[layer];
          float trackLength = computeTrackLength(o2track, layerRadius, dBz);

          if (trackLength > 0) {
            hitMap |= (1 << layer);
            histos.fill(HIST("hHitLayers"), layer);
            nHitLayers++;
          }
        }
      }

      histos.fill(HIST("hHitMultiplicity"), nHitLayers);
      h2dHitsPerTrackVsP[truePdgIdx]->Fill(rigidity, nHitLayers);

      std::vector<float> validToTs;

      for (int layer = kMinLayerForTruncation; layer < maxBarrelLayers.value; ++layer) {
        if ((hitMap >> layer) & 0x1) {
          TH1F* totHist = mToTLUT->getHistogramForSampling(truePdgIdx, layer, binnedEta, binnedPt);

          if (totHist && totHist->GetEntries() > 1) {
            float sampledToT = totHist->GetRandom();
            validToTs.push_back(sampledToT);
          }
        }
      }

      truncatedMeanToT = -1.0f;
      const size_t nValid = validToTs.size();
      size_t nUse = 0;

      if (nValid == kMidLowValidHits || nValid == kMidHighValidHits)
        nUse = kMaxValidHitsForTruncation34;
      else if (nValid == kMinValidHits || nValid == kLowValidHits)
        nUse = kMaxValidHitsForTruncation12;
      else if (nValid == kHighValidHits1 || nValid == kHighValidHits2)
        nUse = kMaxValidHitsForTruncation56;
      else if (nValid >= kMaxValidHits)
        nUse = kMaxValidHitsForTruncation7Plus;

      if (nUse > 0 && nValid >= nUse) {
        std::sort(validToTs.begin(), validToTs.end());
        float sum = 0.0f;
        for (size_t i = 0; i < nUse; ++i) {
          sum += validToTs[i];
        }
        truncatedMeanToT = sum / static_cast<float>(nUse);

        histos.fill(HIST("hToTvsPt"), p, truncatedMeanToT);
        histos.fill(HIST("hToTvsP"), rigidity, truncatedMeanToT);
        h2dToTvsPperParticle[truePdgIdx]->Fill(rigidity, truncatedMeanToT);
      }

      nSigmaValues.fill(999.f);

      if (truncatedMeanToT > 0) {
        for (size_t iHyp = 0; iHyp < mHypothesisPdgCodes.size(); ++iHyp) {
          int hypPdgCode = mHypothesisPdgCodes[iHyp];
          int hypPdgIdx = mToTLUT->getPdgIndex(hypPdgCode);

          if (hypPdgIdx == -1) {
            nSigmaValues[iHyp] = 999.f;
            continue;
          }

          float expectedToT = getToTMeanFromMomentumSlice(h2dToTvsPperParticle[hypPdgIdx], rigidity);
          float resolution = getToTResolutionFromMomentumSlice(h2dToTvsPperParticle[hypPdgIdx], rigidity);

          if (expectedToT > 0 && resolution > 0) {
            nSigmaValues[iHyp] = calculateNsigma(truncatedMeanToT, expectedToT, resolution);
            h2dBarrelNsigmaTrue[truePdgIdx][iHyp]->Fill(rigidity, nSigmaValues[iHyp]);
          } else {
            nSigmaValues[iHyp] = 999.f;
          }
        }
      }

      tableUpgradeTrkPidSignals(truncatedMeanToT);
      tableUpgradeTrkPids(nSigmaValues[0], nSigmaValues[1], nSigmaValues[2], nSigmaValues[3],
                          nSigmaValues[4], nSigmaValues[5], nSigmaValues[6], nSigmaValues[7], nSigmaValues[8]);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyTrackerPid>(cfgc)};
}
