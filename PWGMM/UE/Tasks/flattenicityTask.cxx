// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General
// Public License v3 (GPL Version 3), copied verbatim in the file "COPYING".

/// \file flattenicityTask.cxx
/// \brief Run-3 calibration task for FT0 flattenicity
/// \author Eisha Rani
/// \since September 2026

#include "PWGLF/Utils/inelGt.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/EndOfStreamContext.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TString.h>
#include <TStyle.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

using FullTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::TrackSelectionExtension>;

using CollisionsWithMult = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::PVMults>;
using CollisionsNoMult = soa::Join<aod::Collisions, aod::EvSels>;
// Reco collisions grouped per MC collision (as in the DeepSeek/Monday version)
using ColMCRecTable = soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions, aod::EvSels>>;

struct FlattenicityTask {

  Service<o2::framework::O2DatabasePDG> pdg{};

  // ========================================================================
  // FT0 geometry
  // ========================================================================
  static constexpr int NSectorsA = 24;
  static constexpr int NSectorsC = 28;
  static constexpr int NSectorsCombined = NSectorsA + NSectorsC;
  static constexpr int ChannelsPerSector = 4;
  static constexpr int NchA = NSectorsA * ChannelsPerSector;
  static constexpr int NchC = NSectorsC * ChannelsPerSector;
  static constexpr float FT0AEtaMin = 3.5;
  static constexpr float FT0AEtaMax = 4.9;
  static constexpr float FT0CEtaMin = -3.3;
  static constexpr float FT0CEtaMax = -2.1;
  static constexpr int NPhysicalPrimaryBit = 0x4;
  static constexpr float ChargeEpsilon = 1e-6f;

  // ------------------------------------------------------------------------
  // Real per-module (x, y[, z]) positions from the O2 FT0 geometry. All values
  // verified against AliceO2 (dev):
  //  - x,y: Detectors/FIT/FT0/base/src/Geometry.cxx setAsideModules() (A) and
  //    setCsideModules() (C, ZdetC = 82 cm, compensated radius 82 + 5.5/2);
  //    ZdetA = 335.5 cm (Geometry.h).
  //  - Global placement, Detectors/FIT/FT0/simulation/src/Detector.cxx:
  //    A: translation only in x,y; z of the module centre is about 333.7 cm
  //    (frame at zdetA + 0.63, local z = -5 + 2.6). We use the nominal
  //    ZdetA = 335.5 cm; the eta difference is ~0.01-0.02 (negligible).
  //    C: module local z = zc2 - 80 (line 155, marked "fix later" by the
  //    authors), frame at z = -zdetC with rotC = axes (1,0,0),(0,1,0),
  //    (0,0,-1), i.e. a pure z -> -z reflection (x,y unchanged):
  //    global z = -82 - (zc2 - 80) = -(zc2 + 2)  ->  -84.8 ... -86.2 cm.
  //  - Detector.cxx also places both mother volumes with a +30 cm y
  //    translation (TGeoTranslation(0, 30., ...)). This is deliberately NOT
  //    applied here: applying it moves the eta coverage to A ~2.7-3.9 and
  //    C ~-2.7 to -1.3 (nominal: A 3.5-4.9, C -3.3 to -2.1), and a test with
  //    the shift applied lowered the truth-vs-reco correlation from
  //    r = 0.79 to 0.46 (Average) and 0.63 to 0.42 (Combined).
  //  Uncertainty: module-centre vs quartz-face depth is ~+-2 cm (eta shift
  //  < 0.03); the two odd/even C-module reflections give a mirror-symmetric
  //  set of positions, so flattenicity does not depend on module labelling.
  // ------------------------------------------------------------------------
  // Strict mode (ON by default, cfgRequireModuleHit = true): only count a
  // particle if its projected hit point lies on the module face (square,
  // half-width = 2 x half-channel width = 2.65 cm, from
  // calculateChannelCenter() in O2). Particles pointing into the gaps between
  // modules are dropped. Set to false to assign every particle inside the eta
  // window to its nearest module instead.
  Configurable<bool> cfgRequireModuleHit{"cfgRequireModuleHit", true, "Require projected hit inside module face"};
  static constexpr float kModuleHalfWidth = 2.65f;

  static constexpr float kZDetA = 335.5f;
  static constexpr float kModuleAx[NSectorsA] = {-12.2500f, -6.1500f, -0.0500f, 6.1500f, 12.2500f, -12.2500f, -6.1500f, -0.0500f, 6.1500f, 12.2500f, -13.5800f, -7.4800f, 7.4800f, 13.5800f, -12.2500f, -6.1500f, 0.0500f, 6.1500f, 12.2500f, -12.2500f, -6.1500f, 0.0500f, 6.1500f, 12.2500f};
  static constexpr float kModuleAy[NSectorsA] = {12.2000f, 12.2000f, 13.5300f, 12.2000f, 12.2000f, 6.1000f, 6.1000f, 7.4300f, 6.1000f, 6.1000f, 0.0000f, 0.0000f, 0.0000f, 0.0000f, -6.1000f, -6.1000f, -7.4300f, -6.1000f, -6.1000f, -12.2000f, -12.2000f, -13.5300f, -12.2000f, -12.2000f};
  static constexpr float kModuleCx[NSectorsC] = {-15.5426f, 15.5426f, -15.5069f, 15.5069f, -9.3296f, 9.3296f, -9.3296f, 9.3296f, -9.3296f, 9.3296f, -3.1011f, 3.1011f, -3.1286f, 3.1286f, 3.1010f, -3.1010f, 3.1286f, -3.1286f, 9.3295f, -9.3295f, 9.3296f, -9.3296f, 9.3296f, -9.3296f, 15.5425f, -15.5425f, 15.5068f, -15.5068f};
  static constexpr float kModuleCy[NSectorsC] = {3.2659f, -3.2659f, 9.4726f, -9.4726f, 3.2436f, -3.2436f, 9.4726f, -9.4726f, 15.6494f, -15.6494f, 9.4725f, -9.4725f, 15.6982f, -15.6982f, 9.4724f, -9.4724f, 15.6981f, -15.6981f, 3.2435f, -3.2435f, 9.4725f, -9.4725f, 15.6493f, -15.6493f, 3.2659f, -3.2659f, 9.4725f, -9.4725f};
  static constexpr float kModuleCz[NSectorsC] = {-85.2486f, -85.2486f, -84.7790f, -84.7790f, -86.1724f, -86.1724f, -85.7006f, -85.7006f, -84.7685f, -84.7685f, -86.1618f, -86.1618f, -85.2246f, -85.2246f, -86.1619f, -86.1619f, -85.2247f, -85.2247f, -86.1725f, -86.1725f, -85.7006f, -85.7006f, -84.7685f, -84.7685f, -85.2486f, -85.2486f, -84.7791f, -84.7791f};

  // Official FT0 channel ID -> geometry module index, built by inverting
  // localChannelOrder[] from AliceO2 Detectors/FIT/FT0/base/include/FT0Base/
  // Geometry.h (localChannelOrder[4*module + k] = official channel; A modules
  // 0-23, C modules 24-51). Module index = same order as kModuleA*/kModuleC*
  // above, so truth-level and reco-level sectors are the same physical
  // modules. FT0-C table is indexed by (official channel - 96), i.e. the
  // channelC() value stored in AO2D (see QA/hFT0ChannelC to confirm).
  static constexpr int kChannelToModuleA[NchA] = {12, 12, 12, 12, 8, 8, 8, 8, 7, 7, 7, 7, 6, 6, 6, 6, 11, 11, 11, 11, 15, 15, 15, 15, 16, 16, 16, 16, 17, 17, 17, 17, 13, 13, 13, 13, 9, 9, 9, 9, 4, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 5, 5, 5, 5, 10, 10, 10, 10, 14, 14, 14, 14, 19, 19, 19, 19, 20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 18, 18, 18, 18};
  static constexpr int kChannelToModuleC[NchC] = {18, 18, 18, 18, 20, 20, 20, 20, 14, 14, 14, 14, 10, 10, 10, 10, 6, 6, 6, 6, 4, 4, 4, 4, 19, 19, 19, 19, 21, 21, 21, 21, 15, 15, 15, 15, 11, 11, 11, 11, 7, 7, 7, 7, 5, 5, 5, 5, 24, 24, 24, 24, 26, 26, 26, 26, 22, 22, 22, 22, 16, 16, 16, 16, 12, 12, 12, 12, 8, 8, 8, 8, 2, 2, 2, 2, 0, 0, 0, 0, 25, 25, 25, 25, 27, 27, 27, 27, 23, 23, 23, 23, 17, 17, 17, 17, 13, 13, 13, 13, 9, 9, 9, 9, 3, 3, 3, 3, 1, 1, 1, 1};

  static constexpr int NFlatBins = 10000;
  static constexpr int NCorrTitleUpdateInterval = 500; // refresh Pearson r in the titles every N filled events

  // ========================================================================
  // Histograms
  // ========================================================================
  HistogramRegistry histos{
    "histos",
    {
      {"QA/hEvents", "Event selection;selection;Entries", {HistType::kTH1F, {{5, 0, 5}}}},
      {"QA/hNchINEL", "INEL>0 multiplicity (|eta|<0.8);N_{ch};Entries", {HistType::kTH1F, {{100, -0.5, 99.5}}}},

      // ---- MC truth: 52-sector flattenicity only ----
      {"MCTruth/hFlattenicityAvgSector_fine",
       "MC Truth - Average, 24+28 sectors;1-#rho;Entries",
       {HistType::kTH1F, {{NFlatBins, 0.0, 1.0}}}},
      {"MCTruth/hFlattenicityCombSector_fine",
       "MC Truth - Combined, 52 sectors;1-#rho;Entries",
       {HistType::kTH1F, {{NFlatBins, 0.0, 1.0}}}},

      // ---- MC reco ----
      {"MCReco/hFlattenicityFT0Average_fine",
       "MC Reco - FT0 Average of FT0A and FT0C;1-#rho;Entries",
       {HistType::kTH1F, {{NFlatBins, 0.0, 1.0}}}},
      {"MCReco/hFlattenicityFT0Combined_fine",
       "MC Reco - FT0 Combined FT0A+FT0C;1-#rho;Entries",
       {HistType::kTH1F, {{NFlatBins, 0.0, 1.0}}}},

      // ---- Data ----
      {"Data/hFlattenicityFT0Average_fine",
       "Data - FT0 Average of FT0A and FT0C;1-#rho;Entries",
       {HistType::kTH1F, {{NFlatBins, 0.0, 1.0}}}},
      {"Data/hFlattenicityFT0Combined_fine",
       "Data - FT0 Combined FT0A+FT0C;1-#rho;Entries",
       {HistType::kTH1F, {{NFlatBins, 0.0, 1.0}}}},

      // ---- Correlation matrices: x = measured (reco) 1-rho, y = true (MC particles) 1-rho ----
      {"MCReco/hTruthVsRecoCombined",
       "Combined: reco vs truth;1-#rho (reco, FT0 amplitude);1-#rho (truth, MC particles)",
       {HistType::kTH2F, {{200, 0.0, 1.0}, {200, 0.0, 1.0}}}},
      {"MCReco/hTruthVsRecoAverage",
       "Average: reco vs truth;1-#rho (reco, FT0 amplitude);1-#rho (truth, MC particles)",
       {HistType::kTH2F, {{200, 0.0, 1.0}, {200, 0.0, 1.0}}}},
    }};

  // ------------------------------------------------------------------------
  // pT axis (GeV/c) built from segments: number of bins, first edge, bin width.
  // The last edge is PtMax.
  struct PtSegment {
    int nBins;
    double start;
    double step;
  };
  static constexpr PtSegment PtSegments[] = {{20, 0.0, 0.05}, {20, 1.0, 0.1}, {14, 3.0, 0.5}, {20, 10.0, 2.0}};
  static constexpr double PtMax = 50.0;
  // pT vs flattenicity (flattenicity on x, finely binned; pT on y). The flattenicity
  // class (e.g. 0-1%) is NOT fixed here: it is extracted offline from the event
  // distribution of flattenicity (percentiles), see plotPtSpectraFlatClass.C.
  std::shared_ptr<TH2> hPtVsFlatTruthAvg;
  std::shared_ptr<TH2> hPtVsFlatTruthComb;
  std::shared_ptr<TH2> hPtVsFlatRecoAvg;
  std::shared_ptr<TH2> hPtVsFlatRecoComb;

  // ========================================================================
  // Configurables
  // ========================================================================
  Configurable<float> cfgPtMin{"cfgPtMin", 0.0, "Minimum pT (GeV/c) - require pT > 0"};
  Configurable<float> cfgEtaMax{"cfgEtaMax", 0.8, "Maximum |eta| for INEL>0 and track selection"};
  Configurable<float> cfgVzMax{"cfgVzMax", 10.0, "Maximum absolute collision vertex z (cm)"};
  Configurable<bool> cfgApplySel8{"cfgApplySel8", true, "Apply sel8 event selection"};
  Configurable<bool> cfgApplyGoodITSLayersAll{"cfgApplyGoodITSLayersAll", true, "Apply kIsGoodITSLayersAll event selection"};
  Configurable<bool> cfgRequireInelGt0{"cfgRequireInelGt0", true, "Require INEL>0"};
  Configurable<bool> cfgPileupRejection{"cfgPileupRejection", true, "Reject same-bunch pileup (kNoSameBunchPileup)"};
  Configurable<bool> cfgGoodZvertex{"cfgGoodZvertex", true, "Require good z-vertex agreement PV-vs-FT0 (kIsGoodZvtxFT0vsPV)"};

  Configurable<int> cfgITSreq{"cfgITSreq", 0, "ITS requirement (0 = Run3ITSibAny)"};

  Configurable<bool> cfgIsMCReco{"cfgIsMCReco", true, "true = fill MCReco/, false = fill Data/"};

  Configurable<int> cfgNCrossedRowsTPC{"cfgNCrossedRowsTPC", 70, "Minimum TPC crossed rows"};
  Configurable<float> cfgChi2PerClusterTPC{"cfgChi2PerClusterTPC", 4.0, "Maximum TPC chi2 per cluster"};
  Configurable<float> cfgChi2PerClusterITS{"cfgChi2PerClusterITS", 36.0, "Maximum ITS chi2 per cluster"};
  Configurable<float> cfgDCAZ{"cfgDCAZ", 0.1, "Maximum absolute DCAz (cm)"};
  Configurable<bool> cfgRequireGoldenChi2{"cfgRequireGoldenChi2", true, "Require golden chi2"};

  TrackSelection myTrkSel;

  // Counter and helper to write the Pearson r into the truth-vs-reco titles.
  // Called periodically during the job and again at the end of the stream, so the
  // title is correct even if endOfStream runs after the histograms are written.
  int nCorrFills = 0;

  void updateCorrTitles()
  {
    auto hComb = histos.get<TH2>(HIST("MCReco/hTruthVsRecoCombined"));
    if (hComb && hComb->GetEntries() > 0) {
      hComb->SetTitle(Form("Combined: reco vs truth (Pearson r = %.3f);1-#rho (reco, FT0);1-#rho (truth, MC)", hComb->GetCorrelationFactor()));
      hComb->SetOption("COLZ");
    }
    auto hAvg = histos.get<TH2>(HIST("MCReco/hTruthVsRecoAverage"));
    if (hAvg && hAvg->GetEntries() > 0) {
      hAvg->SetTitle(Form("Average: reco vs truth (Pearson r = %.3f);1-#rho (reco, FT0);1-#rho (truth, MC)", hAvg->GetCorrelationFactor()));
      hAvg->SetOption("COLZ");
    }
  }

  // ========================================================================
  // init
  // ========================================================================
  void init(InitContext&)
  {
    gStyle->SetOptStat("em");
    auto itsMatching = static_cast<TrackSelection::GlobalTrackRun3ITSMatching>(cfgITSreq.value);
    myTrkSel = getGlobalTrackSelectionRun3ITSMatch(itsMatching);

    myTrkSel.SetPtRange(cfgPtMin.value, 1e10f);
    myTrkSel.SetEtaRange(-cfgEtaMax.value, cfgEtaMax.value);
    myTrkSel.SetRequireGoldenChi2(cfgRequireGoldenChi2.value);
    myTrkSel.SetMaxChi2PerClusterTPC(cfgChi2PerClusterTPC.value);
    myTrkSel.SetMaxChi2PerClusterITS(cfgChi2PerClusterITS.value);
    myTrkSel.SetMinNCrossedRowsTPC(cfgNCrossedRowsTPC.value);
    myTrkSel.SetMaxDcaZ(cfgDCAZ.value);
    myTrkSel.SetMaxDcaXYPtDep([](float /*pt*/) { return 10000.f; });

    myTrkSel.print();

    // pT axis: 0-1 GeV/c in 0.05, 1-3 in 0.1, 3-10 in 0.5, 10-50 in 2
    std::vector<double> ptEdges;
    for (const auto& segment : PtSegments) {
      for (int i = 0; i < segment.nBins; ++i) {
        ptEdges.push_back(segment.start + segment.step * i);
      }
    }
    ptEdges.push_back(PtMax);
    const AxisSpec ptAxisSpec{ptEdges, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec flatAxisSpecTruthAvg{NFlatBins, 0.0, 1.0, "flattenicity (MC truth, average)"};
    const AxisSpec flatAxisSpecTruthComb{NFlatBins, 0.0, 1.0, "flattenicity (MC truth, combined)"};
    const AxisSpec flatAxisSpecRecoAvg{NFlatBins, 0.0, 1.0, "flattenicity (reco, FT0 average)"};
    const AxisSpec flatAxisSpecRecoComb{NFlatBins, 0.0, 1.0, "flattenicity (reco, FT0 combined)"};
    hPtVsFlatTruthAvg = histos.add<TH2>("MCTruth/hPtVsFlatTruthAvg", "MC truth: pT vs flattenicity (average);flattenicity;#it{p}_{T}",
                                        HistType::kTH2F, {flatAxisSpecTruthAvg, ptAxisSpec});
    hPtVsFlatTruthComb = histos.add<TH2>("MCTruth/hPtVsFlatTruthComb", "MC truth: pT vs flattenicity (combined);flattenicity;#it{p}_{T}",
                                         HistType::kTH2F, {flatAxisSpecTruthComb, ptAxisSpec});
    hPtVsFlatRecoAvg = histos.add<TH2>("MCReco/hPtVsFlatRecoAvg", "MC reco: pT vs flattenicity (FT0 average);flattenicity;#it{p}_{T}",
                                       HistType::kTH2F, {flatAxisSpecRecoAvg, ptAxisSpec});
    hPtVsFlatRecoComb = histos.add<TH2>("MCReco/hPtVsFlatRecoComb", "MC reco: pT vs flattenicity (FT0 combined);flattenicity;#it{p}_{T}",
                                        HistType::kTH2F, {flatAxisSpecRecoComb, ptAxisSpec});
  }

  // ========================================================================
  // Helpers
  // ========================================================================
  bool isChargedParticle(int pdgCode)
  {
    const auto* pdgParticle = pdg->GetParticle(pdgCode);
    if (!pdgParticle)
      return false;
    if (std::abs(pdgParticle->Charge()) < ChargeEpsilon)
      return false;
    return true;
  }

  template <std::size_t N>
  float computeRho(const std::array<float, N>& counts)
  {
    float total = 0.0f;
    for (std::size_t i = 0; i < N; ++i)
      total += counts[i];
    if (total <= 0.0f)
      return -1.0f;
    const float mean = total / static_cast<float>(N);
    if (mean <= 0.0f)
      return -1.0f;
    float sumSq = 0.0f;
    for (std::size_t i = 0; i < N; ++i) {
      const float diff = counts[i] - mean;
      sumSq += diff * diff;
    }
    return std::sqrt(sumSq) / (static_cast<float>(N) * mean);
  }

  template <std::size_t N>
  float computeFlattenicity(const std::array<float, N>& counts)
  {
    const float rho = computeRho(counts);
    if (rho < 0.0f)
      return -1.0f;
    return 1.0f - rho;
  }

  // channel -> physical module (see kChannelToModule* above). Returns -1 for
  // an out-of-range channel ID.
  int getFT0ASector(int channel)
  {
    if (channel < 0 || channel >= NchA)
      return -1;
    return kChannelToModuleA[channel];
  }
  int getFT0CSector(int channel)
  {
    if (channel < 0 || channel >= NchC)
      return -1;
    return kChannelToModuleC[channel];
  }

  // Project a trajectory of given (eta, phi) from the IP onto the plane
  // z = zPlane (signed - negative for FT0-C) and return the (x, y) it hits.
  static void projectToPlane(float eta, float phi, float zPlane, float vz, float& x, float& y)
  {
    // vz = z of the (MC) primary vertex: the particle flies from (0,0,vz).
    const float theta = 2.0f * std::atan(std::exp(-eta));
    const float r = (zPlane - vz) * std::tan(theta);
    x = r * std::cos(phi);
    y = r * std::sin(phi);
  }

  // Assign a primary charged particle to the FT0 module (real geometry
  // module center, nearest-neighbour) whose face it actually points at,
  // rather than an assumed uniform eta x phi grid. Returns the module
  // index (0..NSectorsA-1 for FT0-A, 0..NSectorsC-1 for FT0-C) or -1 if
  // outside both arrays' eta acceptance. isFT0A is set accordingly.
  int assignToFT0Module(float eta, float phi, float vz, bool& isFT0A)
  {
    const bool inFT0A = (eta > FT0AEtaMin && eta < FT0AEtaMax);
    const bool inFT0C = (eta > FT0CEtaMin && eta < FT0CEtaMax);
    if (!inFT0A && !inFT0C)
      return -1;
    isFT0A = inFT0A;

    float x = 0.0f, y = 0.0f;
    if (inFT0A) {
      projectToPlane(eta, phi, kZDetA, vz, x, y);
      int best = -1;
      float bestD2 = std::numeric_limits<float>::max();
      for (int m = 0; m < NSectorsA; ++m) {
        const float dx = x - kModuleAx[m];
        const float dy = y - kModuleAy[m];
        const float d2 = dx * dx + dy * dy;
        if (d2 < bestD2) {
          bestD2 = d2;
          best = m;
        }
      }
      if (cfgRequireModuleHit && best >= 0 &&
          (std::abs(x - kModuleAx[best]) > kModuleHalfWidth || std::abs(y - kModuleAy[best]) > kModuleHalfWidth))
        return -1;
      return best;
    }
    // FT0-C: project using a representative plane z (module z varies by
    // < 1.5 cm across the 84.8-86.2 cm concave support, negligible for the
    // nearest-module search below).
    projectToPlane(eta, phi, -85.6f, vz, x, y);
    int best = -1;
    float bestD2 = std::numeric_limits<float>::max();
    for (int m = 0; m < NSectorsC; ++m) {
      const float dx = x - kModuleCx[m];
      const float dy = y - kModuleCy[m];
      const float d2 = dx * dx + dy * dy;
      if (d2 < bestD2) {
        bestD2 = d2;
        best = m;
      }
    }
    if (cfgRequireModuleHit && best >= 0 &&
        (std::abs(x - kModuleCx[best]) > kModuleHalfWidth || std::abs(y - kModuleCy[best]) > kModuleHalfWidth))
      return -1;
    return best;
  }

  template <typename FT0>
  bool computeFT0Flattenicities(const FT0& ft0, float& flattenicityAverage, float& flattenicityCombined, float& sumAmplitude)
  {
    std::array<float, NSectorsA> countsA{};
    std::array<float, NSectorsC> countsC{};
    std::array<float, NSectorsCombined> countsCombined{};
    sumAmplitude = 0.0f;

    for (std::size_t i = 0; i < ft0.amplitudeA().size(); ++i) {
      const int channel = static_cast<uint8_t>(ft0.channelA()[i]);
      const int sector = getFT0ASector(channel);
      const float amplitude = ft0.amplitudeA()[i];
      if (sector >= 0 && sector < NSectorsA) {
        countsA[sector] += amplitude;
        countsCombined[sector] += amplitude;
      }
      sumAmplitude += amplitude;
    }
    for (std::size_t i = 0; i < ft0.amplitudeC().size(); ++i) {
      const int channel = static_cast<uint8_t>(ft0.channelC()[i]);
      const int sector = getFT0CSector(channel);
      const float amplitude = ft0.amplitudeC()[i];
      if (sector >= 0 && sector < NSectorsC) {
        countsC[sector] += amplitude;
        countsCombined[NSectorsA + sector] += amplitude;
      }
      sumAmplitude += amplitude;
    }

    // Both outputs default to -1 (= "not valid"); callers check >= 0 before filling.
    flattenicityAverage = -1.0f;
    flattenicityCombined = -1.0f;

    // Combined: all 52 modules pooled. Valid as soon as ANY module has signal
    // (same behaviour as MC truth), so one-sided events are kept.
    const float rhoCombined = computeRho(countsCombined);
    if (rhoCombined >= 0.0f) {
      flattenicityCombined = 1.0f - rhoCombined;
    }

    // Average: needs both FT0-A and FT0-C valid (same behaviour as MC truth).
    const float rhoA = computeRho(countsA);
    const float rhoC = computeRho(countsC);
    if (rhoA >= 0.0f && rhoC >= 0.0f) {
      flattenicityAverage = 1.0f - 0.5f * (rhoA + rhoC);
    }

    // Event is usable if at least one of the two estimators is valid.
    return (flattenicityAverage >= 0.0f || flattenicityCombined >= 0.0f);
  }

  template <typename Track>
  bool isTrackSelected(const Track& track)
  {
    if (track.pt() <= cfgPtMin)
      return false;
    if (std::abs(track.eta()) > cfgEtaMax)
      return false;

    for (int i = 0; i < static_cast<int>(TrackSelection::TrackCuts::kNCuts); ++i) {
      if (i == static_cast<int>(TrackSelection::TrackCuts::kDCAxy))
        continue;
      if (!myTrkSel.IsSelected(track, static_cast<TrackSelection::TrackCuts>(i)))
        return false;
    }

    const float dcaXyCut = 0.0105f + 0.0350f / std::pow(track.pt(), 1.1f);
    return std::abs(track.dcaXY()) <= dcaXyCut;
  }

  template <typename Tracks>
  int countSelectedTracks(const Tracks& tracks, int64_t collIdx)
  {
    int n = 0;
    for (const auto& track : tracks) {
      if (track.collisionId() != collIdx)
        continue;
      if (!isTrackSelected(track))
        continue;
      ++n;
    }
    return n;
  }

  // Fill reco-level pT (selected tracks of this collision) vs the event flattenicity.
  // Same event conditions as the 1D fine flattenicity histograms, so the percentile
  // classes taken from those histograms are consistent with these pT spectra.
  template <typename Tracks>
  void fillRecoPtVsFlat(const Tracks& tracks, int64_t collIdx, float flatAvg, float flatComb)
  {
    const bool okAvg = (flatAvg >= 0.0f && flatAvg <= 1.0f);
    const bool okComb = (flatComb >= 0.0f && flatComb <= 1.0f);
    if (!okAvg && !okComb)
      return;
    for (const auto& track : tracks) {
      if (track.collisionId() != collIdx)
        continue;
      if (!isTrackSelected(track))
        continue;
      if (okAvg)
        hPtVsFlatRecoAvg->Fill(flatAvg, track.pt());
      if (okComb)
        hPtVsFlatRecoComb->Fill(flatComb, track.pt());
    }
  }

  // ========================================================================
  // processMC - MC truth per-module plots, truth flattenicity, and
  // truth-vs-reco correlation. Selection criteria as in the DeepSeek/Monday
  // version: truth uses only |vz| < cfgVzMax (no truth INEL>0, no both-sides
  // requirement); the correlation uses the first reco collision of the MC
  // collision with the standard reco event selection. Geometry: real FT0
  // module positions (see kModule* / kChannelToModule* above).
  // ========================================================================
  void processMC(aod::McCollisions::iterator const& mcCollision,
                 ColMCRecTable const& recCols,
                 aod::McParticles const& mcParticles,
                 aod::FT0s const& ft0s)
  {
    if (std::abs(mcCollision.posZ()) > cfgVzMax)
      return;
    const float vzMC = mcCollision.posZ();

    std::array<float, NSectorsA> truthA{};
    std::array<float, NSectorsC> truthC{};
    std::array<float, NSectorsCombined> truthCombArr{};

    int nChINEL = 0;
    std::vector<float> truthPt; // primary charged, |eta| < cfgEtaMax, pT > cfgPtMin
    for (const auto& particle : mcParticles) {
      if (particle.mcCollisionId() != mcCollision.globalIndex())
        continue;
      if ((particle.flags() & NPhysicalPrimaryBit) == 0)
        continue;
      if (!isChargedParticle(particle.pdgCode()))
        continue;
      if (particle.pt() <= cfgPtMin)
        continue;
      if (std::abs(particle.eta()) < cfgEtaMax) {
        ++nChINEL;
        truthPt.push_back(particle.pt());
      }

      bool isFT0A = false;
      const int module = assignToFT0Module(particle.eta(), particle.phi(), vzMC, isFT0A);
      if (module < 0)
        continue;

      if (isFT0A) {
        truthA[module] += 1.0f;
        truthCombArr[module] += 1.0f;
      } else {
        truthC[module] += 1.0f;
        truthCombArr[NSectorsA + module] += 1.0f;
      }
    }

    histos.fill(HIST("QA/hEvents"), 1);
    histos.fill(HIST("QA/hNchINEL"), nChINEL);

    const float rhoA_sec = computeRho(truthA);
    const float rhoC_sec = computeRho(truthC);
    const float truthAvg = (rhoA_sec >= 0.0f && rhoC_sec >= 0.0f)
                             ? 1.0f - 0.5f * (rhoA_sec + rhoC_sec)
                             : -1.0f;
    const float truthCombined = computeFlattenicity(truthCombArr);

    if (truthAvg >= 0.0f && truthAvg <= 1.0f)
      histos.fill(HIST("MCTruth/hFlattenicityAvgSector_fine"), truthAvg);
    if (truthCombined >= 0.0f)
      histos.fill(HIST("MCTruth/hFlattenicityCombSector_fine"), truthCombined);

    // pT vs truth flattenicity (same event conditions as the fine histograms above)
    for (const auto& pt : truthPt) {
      if (truthAvg >= 0.0f && truthAvg <= 1.0f)
        hPtVsFlatTruthAvg->Fill(truthAvg, pt);
      if (truthCombined >= 0.0f)
        hPtVsFlatTruthComb->Fill(truthCombined, pt);
    }

    // ---- correlation matrix: measured (reco) vs true ----
    if (recCols.size() == 0)
      return;
    auto recCol = *recCols.begin();

    if (cfgApplySel8 && !recCol.sel8())
      return;
    if (cfgPileupRejection && !recCol.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return;
    if (cfgGoodZvertex && !recCol.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return;
    if (std::abs(recCol.posZ()) > cfgVzMax)
      return;
    if (cfgApplyGoodITSLayersAll &&
        !recCol.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))
      return;

    aod::FT0s::iterator ft0;
    bool foundFT0 = false;
    for (const auto& f : ft0s) {
      if (f.bcId() == recCol.bcId()) {
        ft0 = f;
        foundFT0 = true;
        break;
      }
    }
    if (!foundFT0)
      return;

    float flatAvg = -1.0f, flatComb = -1.0f, sumAmp = 0.0f;
    if (!computeFT0Flattenicities(ft0, flatAvg, flatComb, sumAmp))
      return;

    if (flatComb >= 0.0f && truthCombined >= 0.0f)
      histos.fill(HIST("MCReco/hTruthVsRecoCombined"), flatComb, truthCombined);
    if (flatAvg >= 0.0f && truthAvg >= 0.0f)
      histos.fill(HIST("MCReco/hTruthVsRecoAverage"), flatAvg, truthAvg);

    if (++nCorrFills % NCorrTitleUpdateInterval == 0)
      updateCorrTitles();
  }

  PROCESS_SWITCH(FlattenicityTask, processMC, "Process MC truth", true);

  // ========================================================================
  // processDataFull  (Hyperloop production schema)
  // ========================================================================
  void processDataFull(CollisionsWithMult::iterator const& collision,
                       aod::FT0s const& ft0s,
                       FullTracks const& tracks)
  {
    if (cfgRequireInelGt0 && !collision.isInelGt0())
      return;
    if (cfgApplySel8 && !collision.sel8())
      return;
    if (cfgPileupRejection && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return;
    if (cfgGoodZvertex && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return;
    if (std::abs(collision.posZ()) > cfgVzMax)
      return;
    if (cfgApplyGoodITSLayersAll &&
        !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))
      return;

    const int nChINEL = countSelectedTracks(tracks, collision.globalIndex());
    histos.fill(HIST("QA/hEvents"), 3);
    histos.fill(HIST("QA/hNchINEL"), nChINEL);

    aod::FT0s::iterator ft0;
    bool foundFT0 = false;
    for (const auto& f : ft0s) {
      if (f.bcId() == collision.bcId()) {
        ft0 = f;
        foundFT0 = true;
        break;
      }
    }
    if (!foundFT0)
      return;

    float flatAvg = -1.0f, flatComb = -1.0f, sumFT0M = 0.0f;
    if (!computeFT0Flattenicities(ft0, flatAvg, flatComb, sumFT0M))
      return;

    if (cfgIsMCReco.value) {
      if (flatAvg >= 0.0f && flatAvg <= 1.0f)
        histos.fill(HIST("MCReco/hFlattenicityFT0Average_fine"), flatAvg);
      if (flatComb >= 0.0f && flatComb <= 1.0f)
        histos.fill(HIST("MCReco/hFlattenicityFT0Combined_fine"), flatComb);
      fillRecoPtVsFlat(tracks, collision.globalIndex(), flatAvg, flatComb);
    } else {
      if (flatAvg >= 0.0f && flatAvg <= 1.0f)
        histos.fill(HIST("Data/hFlattenicityFT0Average_fine"), flatAvg);
      if (flatComb >= 0.0f && flatComb <= 1.0f)
        histos.fill(HIST("Data/hFlattenicityFT0Combined_fine"), flatComb);
    }
  }

  PROCESS_SWITCH(FlattenicityTask, processDataFull, "Process data / MC-reco with full selection", false);

  // ========================================================================
  // processDataNoMult  (local truth-only MC — no O2mult / O2pvmult)
  // ========================================================================
  void processDataNoMult(CollisionsNoMult::iterator const& collision,
                         aod::FT0s const& ft0s,
                         FullTracks const& tracks)
  {
    if (cfgApplySel8 && !collision.sel8())
      return;
    if (cfgPileupRejection && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      return;
    if (cfgGoodZvertex && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      return;
    if (std::abs(collision.posZ()) > cfgVzMax)
      return;
    if (cfgApplyGoodITSLayersAll &&
        !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll))
      return;
    // NOTE: collision.isInelGt0() intentionally not applied here — this
    // reduced schema (CollisionsNoMult) has no aod::Mults/aod::PVMults.

    const int nChINEL = countSelectedTracks(tracks, collision.globalIndex());
    histos.fill(HIST("QA/hEvents"), 3);
    histos.fill(HIST("QA/hNchINEL"), nChINEL);

    aod::FT0s::iterator ft0;
    bool foundFT0 = false;
    for (const auto& f : ft0s) {
      if (f.bcId() == collision.bcId()) {
        ft0 = f;
        foundFT0 = true;
        break;
      }
    }
    if (!foundFT0)
      return;

    float flatAvg = -1.0f, flatComb = -1.0f, sumFT0M = 0.0f;
    if (!computeFT0Flattenicities(ft0, flatAvg, flatComb, sumFT0M))
      return;

    if (cfgIsMCReco.value) {
      if (flatAvg >= 0.0f && flatAvg <= 1.0f)
        histos.fill(HIST("MCReco/hFlattenicityFT0Average_fine"), flatAvg);
      if (flatComb >= 0.0f && flatComb <= 1.0f)
        histos.fill(HIST("MCReco/hFlattenicityFT0Combined_fine"), flatComb);
      fillRecoPtVsFlat(tracks, collision.globalIndex(), flatAvg, flatComb);
    } else {
      if (flatAvg >= 0.0f && flatAvg <= 1.0f)
        histos.fill(HIST("Data/hFlattenicityFT0Average_fine"), flatAvg);
      if (flatComb >= 0.0f && flatComb <= 1.0f)
        histos.fill(HIST("Data/hFlattenicityFT0Combined_fine"), flatComb);
    }
  }

  PROCESS_SWITCH(FlattenicityTask, processDataNoMult, "Process local truth-only MC (no O2mult)", true);

  // ========================================================================
  // endOfStream - runs once at the end of the analysis
  // ========================================================================
  void endOfStream(EndOfStreamContext&)
  {
    updateCorrTitles();
  }
};

// ============================================================================
// Workflow
// ============================================================================
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<FlattenicityTask>(cfgc)};
}
