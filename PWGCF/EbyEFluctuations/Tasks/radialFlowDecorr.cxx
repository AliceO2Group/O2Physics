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

/// \file radialFlowDecorr.cxx
/// \brief Analysis task for event-by-event radial-flow decorrelation measurement.
/// \author Somadutta Bhatta

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <DetectorsCommonDataFormats/AlignParam.h>
#include <FT0Base/Geometry.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TRandom3.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

struct RadialFlowDecorr {

  // --- fixed constants ---------------------------------------------------------
  static constexpr int KnFt0cCell = 96;
  static constexpr int KIntM = 3; // pT-moment order used in the sums
  static constexpr int KIntK = 3; // weight-power order used in the sums

  // KNEtaMax sizes every fixed-length eta array. The *active* count nEta is a
  // runtime value: 17 for the 0.1-wide observable binning (16 bins + the index-0
  // full-range reference bin), 9 for the 0.2-wide binning (8 bins + reference).
  static constexpr int KNEtaMax = 17;

  static constexpr float KFloatEpsilon = 1e-6f;
  static constexpr float KBinOffset = 0.5f;
  static constexpr float KPhiMin = 0.f;
  static constexpr int KNbinsZvtx = 240;
  static constexpr float KZvtxMin = -12.f;
  static constexpr float KZvtxMax = 12.f;
  static constexpr float KPMin = 0.f;
  static constexpr float KPMax = 10.f;
  static constexpr int KNbinsPt = 200;
  static constexpr float KPtMin = 0.15f;
  static constexpr float KPtMax = 10.f;
  static constexpr float KEtaMin = -1.2f;
  static constexpr float KEtaMax = 1.2f;
  static constexpr int KNbinsPhi = 64;
  static constexpr int KNbinsPtRes = 50;
  static constexpr int KNbinsEtaRes = 100;
  static constexpr int KNbinsVz = 80;
  static constexpr float KVzMin = -40.f;
  static constexpr float KVzMax = 40.f;
  static constexpr int KNbinsEtaFine = 20;
  static constexpr float KEtaFineMax = 1.f;
  static constexpr float KCentMax = 90;

  // Bootstrap: KMaxBoot is the compile-time storage cap; the number actually
  // filled is nBoot = min(cfgNBootstrap, KMaxBoot).
  static constexpr int KMaxBoot = 64;

  enum ECentralityEstimator {
    kCentFT0C = 1,
    kCentFT0M = 2,
    kCentFDDM = 3,
    kCentFV0A = 4
  };
  enum SystemType {
    kPbPb = 1,
    kNeNe = 2,
    kOO = 3,
    kpp = 4
  };

  // Systematic-variation selector. Base is the main measurement; the others each
  // read their own DataMean object for the correlation step and run at 0.2.
  enum ESystType {
    kSystBase = 0,
    kSystDCA,
    kSystEff,
    kSystFlat,
    kSystNEta,
    kSystNITS,
    kSystNTPC,
    kSystPileup,
    kSystVz,
    kNSystType
  };
  // Suffix appended to the DataMean CCDB path per systematic (Base -> no suffix).
  inline static const std::vector<std::string> systSuffix = {
    "", "_systDCA", "_systEff", "_systFlat", "_systNEta",
    "_systNITS", "_systNTPC", "_systPileup", "_systVz"};

  static constexpr float KinvalidCentrality = -1.0f;

  // Observable eta binning, filled at init(): index 0 = full-range reference bin.
  std::vector<float> etaLw;
  std::vector<float> etaUp;
  int nEta = 9; // active eta-bin count (set in init())

  int nBoot = 0;       // active bootstrap samples (set in init())
  bool doBoot = false; // bootstrap active for this run (base data fluc only)

  // --- configurables -----------------------------------------------------------
  Configurable<float> cfgVtxZCut{"cfgVtxZCut", 10.f, "|z_{vtx}| acceptance (cm): collision filter + explicit event/particle vertex checks"};
  Configurable<float> cfgPtMin{"cfgPtMin", 0.2f, "min pT for observables"};
  Configurable<float> cfgPtMax{"cfgPtMax", 5.0f, "max pT for observables"};
  Configurable<float> cfgEtaCut{"cfgEtaCut", 0.8f, "|eta| cut"};
  Configurable<float> cfgCutTracKDcaMaxZ{"cfgCutTracKDcaMaxZ", 2.0f, "Maximum DcaZ"};
  Configurable<float> cfgCutTracKDcaMaxXY{"cfgCutTracKDcaMaxXY", 0.2f, "Maximum DcaXY"};

  Configurable<bool> cfgPtDepDCAxy{"cfgPtDepDCAxy", false, "Use pt-dependent DCAxy cut"};
  Configurable<float> cfgDcaXyP0{"cfgDcaXyP0", 0.0026f, "p0 for DCAxy"};
  Configurable<float> cfgDcaXyP1{"cfgDcaXyP1", 0.005f, "p1 for DCAxy"};
  Configurable<float> cfgDcaXyP2{"cfgDcaXyP2", 1.01f, "p2 for DCAxy"};

  Configurable<bool> cfgPtDepDCAz{"cfgPtDepDCAz", false, "Use pt-dependent DCAz cut"};
  Configurable<float> cfgDcaZP0{"cfgDcaZP0", 0.0026f, "p0 for DCAz"};
  Configurable<float> cfgDcaZP1{"cfgDcaZP1", 0.005f, "p1 for DCAz"};
  Configurable<float> cfgDcaZP2{"cfgDcaZP2", 1.01f, "p2 for DCAz"};

  Configurable<int> cfgITScluster{"cfgITScluster", 1, "Minimum Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 80, "Minimum Number of TPC cluster"};
  Configurable<int> cfgTPCnCrossedRows{"cfgTPCnCrossedRows", 70, "Minimum Number of TPC crossed-rows"};

  Configurable<float> cfgCutPtLower{"cfgCutPtLower", 0.2f, "Lower pT cut (track selection)"};
  Configurable<float> cfgCutPtUpper{"cfgCutPtUpper", 10.0f, "Higher pT cut (track selection)"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "absolute Eta cut"};
  Configurable<int> cfgMinTracksPerEtaBin{"cfgMinTracksPerEtaBin", 2, "Min weighted-track sum required in every narrow eta bin (systNEta; 0 = disabled)"};
  Configurable<int> cfgCentralityChoice{"cfgCentralityChoice", 1, "Which centrality estimator? 1-->FT0C, 2-->FT0M, 3-->FDDM, 4-->FV0A"};
  Configurable<bool> cfgEvSelNoSameBunchPileup{"cfgEvSelNoSameBunchPileup", true, "Pileup removal"};
  Configurable<bool> cfgUseGoodITSLayerAllCut{"cfgUseGoodITSLayerAllCut", true, "Remove time interval with dead ITS zone"};
  Configurable<bool> cfgIsGoodZvtxFT0VsPV{"cfgIsGoodZvtxFT0VsPV", true, "Good Vertexing cut"};

  Configurable<float> cfgPupnSig{"cfgPupnSig", 6.0f, "Additional Pileup Cut"};
  Configurable<bool> cfgApplySigPupCut{"cfgApplySigPupCut", 0, "nSig Pileup Cut"};
  Configurable<bool> cfgApplyLinPupCut{"cfgApplyLinPupCut", 0, "Lin Pileup Cut"};
  Configurable<float> cfgLinPupParam0{"cfgLinPupParam0", 3.0f, "(Upper) Linear Pileup Cut Const"};
  Configurable<float> cfgLinPupParam1{"cfgLinPupParam1", 3.0f, "(Upper) Linear Pileup Slope"};
  Configurable<float> cfgLinPupParam2{"cfgLinPupParam2", 3.0f, "(Lower) Linear Pileup Cut Const"};
  Configurable<float> cfgLinPupParam3{"cfgLinPupParam3", 3.0f, "(Lower) Linear Pileup Slope"};

  Configurable<int> cfgNchPbMax{"cfgNchPbMax", 5000, "Max Nch range for PbPb collisions"};
  Configurable<int> cfgNchOMax{"cfgNchOMax", 800, "Max Nch range for OO collisions"};

  Configurable<int> cfgSys{"cfgSys", 1, "Which collision system? 1-->PbPb, 2-->NeNe, 3-->OO, 4-->pp"};
  Configurable<int> cfgSystType{"cfgSystType", 0, "Systematic variation: 0=Base,1=systDCA,2=systEff,3=systFlat,4=systNEta,5=systNITS,6=systNTPC,7=systPileup,8=systVz"};
  Configurable<int> cfgNBootstrap{"cfgNBootstrap", 30, "Number of Poisson bootstrap samples (base data run only)"};
  Configurable<int> cfgBootstrapSeed{"cfgBootstrapSeed", 0, "TRandom3 seed for bootstrap (0 = machine-random per job)"};

  Configurable<bool> cfgFlat{"cfgFlat", false, "Whether to use flattening weights"};
  Configurable<bool> cfgEff{"cfgEff", false, "Whether to use Efficiency weights"};
  Configurable<bool> cfgZDC{"cfgZDC", false, "Whether to use ZDC for pileup histograms"};

  Configurable<std::string> cfgCCDBurl{"cfgCCDBurl", "https://alice-ccdb.cern.ch", "ccdb url"};
  Configurable<std::string> cfgCCDBUserPath{"cfgCCDBUserPath", "/Users/s/somadutt", "Base CCDB path"};

  ConfigurableAxis cfgAxisCent{"cfgAxisCent", {0.0, 1.0, 5.0, 10, 20, 40, 60, 80, 100}, "centrality axis (percentile)"};

  // --- axes --------------------------------------------------------------------
  AxisSpec centAxis{cfgAxisCent, "Centrality (%)"};
  AxisSpec centAxis1Per{100, 0.0, 100.0, "Centrality (%)"};
  AxisSpec nChAxis{1, 0., 1., "Nch", "Nch"};
  AxisSpec nChAxis2{1, 0., 1., "Nch", "Nch"};

  AxisSpec vzAxis{5, -12.5, 12.5, "Vz"};
  AxisSpec chgAxis{3, -1.5, 1.5};
  AxisSpec pTAxis{{0.0, 0.2, 0.4, 0.6, 0.8, 1, 3, 5, 7, 10}, "pT Axis"};
  AxisSpec phiAxis{KNbinsPhi, KPhiMin, TwoPI, "#phi"};

  // etaFlatAxis is the *fixed* granularity used only for the flattening map, so
  // that a flattening map is reusable across the 0.1/0.2 observable binnings.
  AxisSpec etaFlatAxis{{-0.8, -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8}, "#eta"};

  // etaAxis (physical) and etaBinAxis (integer index) follow the observable
  // binning and are rebuilt in init(). Placeholders here.
  AxisSpec etaAxis{9, -0.9, 0.9, "#eta"};
  AxisSpec etaBinAxis{10, -0.5, 9.5, "#eta bin Number"};

  AxisSpec gapAxis{{-1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1,
                    0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5},
                   "Gap"};
  AxisSpec sumAxis{{-1.5, -1.3, -1.1, -0.9, -0.7, -0.5, -0.3, -0.1,
                    0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5},
                   "Sum"};

  // --- process switches --------------------------------------------------------
  Configurable<bool> cfgRunGetEff{"cfgRunGetEff", false, "Run MC pass to build efficiency/fake maps"};
  Configurable<bool> cfgRunGetMCFlat{"cfgRunGetMCFlat", false, "Run MC to get flattening weights"};
  Configurable<bool> cfgRunMCMean{"cfgRunMCMean", false, "Run MC mean(pT)"};
  Configurable<bool> cfgRunMCFluc{"cfgRunMCFluc", false, "Run MC fluctuations (C2, subevent)"};

  Configurable<bool> cfgRunGetDataFlat{"cfgRunGetDataFlat", false, "Run data get flattening weights"};
  Configurable<bool> cfgRunDataMean{"cfgRunDataMean", false, "Run DATA mean(pT)"};
  Configurable<bool> cfgRunDataFluc{"cfgRunDataFluc", false, "Run DATA fluctuations (C2, subevent)"};

  Service<ccdb::BasicCCDBManager> ccdb{};
  Service<o2::framework::O2DatabasePDG> pdg{};
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  TRandom3 rng; // bootstrap Poisson weights

  // --- persistent state --------------------------------------------------------
  struct InternalState {
    TH3F* hEff = nullptr;
    TH3F* hFake = nullptr;
    THnSparseF* hFlatWeight = nullptr;

    std::vector<std::pair<float, float>> mLimitsNchCent;
    float mMinXNchCent = 0, mMaxXNchCent = 0;

    // MC mean maps (per eta bin): truth / reco / reco-eff-corrected
    TProfile2D* pmeanTruNchEtabinStep2 = nullptr;
    TProfile2D* pmeanRecoNchEtabinStep2 = nullptr;
    TProfile2D* pmeanRecoEffcorrNchEtabinStep2 = nullptr;
    TProfile2D* pmeanMultTruNchEtabinStep2 = nullptr;
    TProfile2D* pmeanMultRecoNchEtabinStep2 = nullptr;
    TProfile2D* pmeanMultRecoEffcorrNchEtabinStep2 = nullptr;

    // Data mean maps (per eta bin)
    TProfile2D* pmeanNchEtabinStep2 = nullptr;
    TProfile2D* pmeanMultNchEtabinStep2 = nullptr;

    TProfile* pmeanFT0AmultpvStep2 = nullptr;
    TProfile* pmeanFT0CmultpvStep2 = nullptr;
  } state;
  o2::ft0::Geometry ft0Det;

  // --- bootstrap replica storage (base data fluc only) -------------------------
  struct BootstrapHists {
    std::array<std::shared_ptr<TProfile2D>, KMaxBoot> meanpTCent{};
    std::array<std::shared_ptr<TProfile2D>, KMaxBoot> meanpTMult{};
    std::array<std::shared_ptr<TProfile2D>, KMaxBoot> c2Cent{};
    std::array<std::shared_ptr<TProfile2D>, KMaxBoot> c2Mult{};
    std::array<std::shared_ptr<TProfile2D>, KMaxBoot> c2SubCent{};
    std::array<std::shared_ptr<TProfile2D>, KMaxBoot> c2SubMult{};
    std::array<std::shared_ptr<TProfile2D>, KMaxBoot> covCent{};
    std::array<std::shared_ptr<TProfile2D>, KMaxBoot> covMult{};
    std::array<std::shared_ptr<TProfile2D>, KMaxBoot> covFT0ACent{};
    std::array<std::shared_ptr<TProfile2D>, KMaxBoot> covFT0AMult{};
    std::array<std::shared_ptr<TProfile2D>, KMaxBoot> covFT0CCent{};
    std::array<std::shared_ptr<TProfile2D>, KMaxBoot> covFT0CMult{};
    std::array<std::shared_ptr<TProfile3D>, KMaxBoot> c2Sub2D{};
    std::array<std::shared_ptr<TProfile3D>, KMaxBoot> gapSum2D{};
    std::array<std::shared_ptr<TProfile3D>, KMaxBoot> cov2D{};
    std::array<std::shared_ptr<TProfile3D>, KMaxBoot> covFT0A2D{};
    std::array<std::shared_ptr<TProfile3D>, KMaxBoot> covFT0C2D{};
  } bs;

  // ===========================================================================
  // Selection helpers
  // ===========================================================================
  template <typename T>
  bool isEventSelected(const T& col)
  {
    histos.fill(HIST("hEvtCount"), 0.5);
    if (!col.sel8()) {
      return false;
    }
    histos.fill(HIST("hEvtCount"), 1.5);
    if (std::abs(col.posZ()) > cfgVtxZCut) {
      return false;
    }
    histos.fill(HIST("hEvtCount"), 2.5);
    if (cfgEvSelNoSameBunchPileup && !col.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    histos.fill(HIST("hEvtCount"), 3.5);
    if (cfgUseGoodITSLayerAllCut && !col.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    histos.fill(HIST("hEvtCount"), 4.5);
    if (cfgIsGoodZvtxFT0VsPV && !col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    histos.fill(HIST("hEvtCount"), 5.5);
    return true;
  }

  bool isPassAddPileup(float multPV, int trksize, float cent)
  {
    auto checkLimits = [](float x, float y, const std::vector<std::pair<float, float>>& limits, float xM, float xMx) {
      if (limits.empty()) {
        return true;
      }
      int bin = 1 + static_cast<int>((x - xM) / (xMx - xM) * (limits.size() - 2));
      if (bin < 1 || bin >= static_cast<int>(limits.size() - 1)) {
        return false;
      }
      return (y >= limits[bin].first && y <= limits[bin].second);
    };
    if (cfgApplySigPupCut) {
      if (!checkLimits(cent, trksize, state.mLimitsNchCent, state.mMinXNchCent, state.mMaxXNchCent)) {
        return false;
      }
      histos.fill(HIST("hEvtCount"), 6.5);
    }
    if (cfgApplyLinPupCut) {
      if (trksize > (cfgLinPupParam0 + cfgLinPupParam1 * multPV)) {
        return false;
      }
      histos.fill(HIST("hEvtCount"), 7.5);
      if (trksize < (cfgLinPupParam2 + cfgLinPupParam3 * multPV)) {
        return false;
      }
      histos.fill(HIST("hEvtCount"), 8.5);
    }
    return true;
  }

  // Minimum weighted-track requirement in every narrow eta bin (systNEta). Loops
  // only over the active bins [1, nEta); the array is sized to KNEtaMax.
  template <std::size_t NetaT, std::size_t NkT>
  bool hasMinTracksInAllEtaBins(const std::array<std::array<double, NkT>, NetaT>& sw)
  {
    const int minTracks = cfgMinTracksPerEtaBin;
    if (minTracks <= 0) {
      return true;
    }
    for (int ieta = 1; ieta < nEta; ++ieta) {
      if (sw[ieta][1] < static_cast<double>(minTracks)) {
        return false;
      }
    }
    return true;
  }
  template <std::size_t NetaT>
  bool hasMinTracksInAllEtaBins(const std::array<double, NetaT>& sw)
  {
    const int minTracks = cfgMinTracksPerEtaBin;
    if (minTracks <= 0) {
      return true;
    }
    for (int ieta = 1; ieta < nEta; ++ieta) {
      if (sw[ieta] < static_cast<double>(minTracks)) {
        return false;
      }
    }
    return true;
  }

  template <typename T>
  bool isTrackSelected(const T& trk)
  {
    histos.fill(HIST("hTrkCount"), 0.5);
    if (trk.sign() == 0) {
      return false;
    }
    histos.fill(HIST("hTrkCount"), 1.5);
    if (!trk.has_collision()) {
      return false;
    }
    histos.fill(HIST("hTrkCount"), 2.5);
    if (!trk.isPVContributor()) {
      return false;
    }
    histos.fill(HIST("hTrkCount"), 3.5);
    if (!(trk.itsNCls() > cfgITScluster)) {
      return false;
    }
    histos.fill(HIST("hTrkCount"), 4.5);
    if (!(trk.tpcNClsFound() >= cfgTPCcluster)) {
      return false;
    }
    histos.fill(HIST("hTrkCount"), 5.5);
    if (!(trk.tpcNClsCrossedRows() >= cfgTPCnCrossedRows)) {
      return false;
    }
    histos.fill(HIST("hTrkCount"), 6.5);
    if (trk.pt() < cfgCutPtLower || trk.pt() > cfgCutPtUpper || std::abs(trk.eta()) > cfgCutEta) {
      return false;
    }
    histos.fill(HIST("hTrkCount"), 7.5);
    if (!trk.isGlobalTrack()) {
      return false;
    }
    histos.fill(HIST("hTrkCount"), 8.5);

    if (cfgPtDepDCAxy) {
      float maxDcaXY = cfgDcaXyP0 + cfgDcaXyP1 / std::pow(trk.pt(), cfgDcaXyP2);
      if (std::abs(trk.dcaXY()) > maxDcaXY) {
        return false;
      }
      histos.fill(HIST("hTrkCount"), 9.5);
    } else {
      if (std::abs(trk.dcaXY()) > cfgCutTracKDcaMaxXY) {
        return false;
      }
      histos.fill(HIST("hTrkCount"), 9.5);
    }
    if (cfgPtDepDCAz) {
      float maxDcaZ = cfgDcaZP0 + cfgDcaZP1 / std::pow(trk.pt(), cfgDcaZP2);
      if (std::abs(trk.dcaZ()) > maxDcaZ) {
        return false;
      }
      histos.fill(HIST("hTrkCount"), 10.5);
    } else {
      if (std::abs(trk.dcaZ()) > cfgCutTracKDcaMaxZ) {
        return false;
      }
      histos.fill(HIST("hTrkCount"), 10.5);
    }
    return true;
  }

  template <typename T>
  bool isParticleSelected(const T& particle)
  {
    auto* pd = pdg->GetParticle(particle.pdgCode());
    if (!pd) {
      return false;
    }
    if (std::abs(pd->Charge()) == 0) {
      return false;
    }
    if (particle.pt() < cfgCutPtLower || particle.pt() > cfgCutPtUpper || std::abs(particle.eta()) > cfgCutEta) {
      return false;
    }
    if (std::abs(particle.vz()) > cfgVtxZCut) {
      return false;
    }
    return true;
  }

  float getCentrality(const auto& col) const
  {
    if (cfgCentralityChoice.value == kCentFT0C) {
      return col.centFT0C();
    }
    if (cfgCentralityChoice.value == kCentFT0M) {
      return col.centFT0M();
    }
    if (cfgCentralityChoice.value == kCentFDDM) {
      return col.centFDDM();
    }
    if (cfgCentralityChoice.value == kCentFV0A) {
      return col.centFV0A();
    }
    return KinvalidCentrality;
  }

  // Inclusive efficiency/fake lookup (no species dependence).
  float getEfficiency(float mult, float pt, float eta, int effidx, bool useEff) const
  {
    if (!useEff) {
      return (effidx == 0) ? 1.0f : 0.0f;
    }
    TH3F* h = (effidx == 0) ? state.hEff : state.hFake;
    if (!h) {
      return -1;
    }
    int ibx = h->GetXaxis()->FindBin(mult);
    int iby = h->GetYaxis()->FindBin(pt);
    int ibz = h->GetZaxis()->FindBin(eta);
    float val = h->GetBinContent(ibx, iby, ibz);
    if (effidx == 0) {
      return (val > 0.f) ? val : 1.0f;
    }
    return val;
  }

  float getFlatteningWeight(float vz, float chg, float pt, float eta, float phi, bool useFlat) const
  {
    if (!useFlat) {
      return 1.0;
    }
    THnSparseF* h = state.hFlatWeight;
    if (!h) {
      return 0.0;
    }
    std::array<int, 5> bins{};
    bins[0] = h->GetAxis(0)->FindBin(vz);
    bins[1] = h->GetAxis(1)->FindBin(chg);
    bins[2] = h->GetAxis(2)->FindBin(pt);
    bins[3] = h->GetAxis(3)->FindBin(eta);
    bins[4] = h->GetAxis(4)->FindBin(phi);
    return h->GetBinContent(bins.data());
  }

  std::vector<o2::detectors::AlignParam>* offsetFT0 = nullptr;
  uint64_t mLastTimestamp = 0;
  double getEtaFT0(uint64_t globalChno, int i)
  {
    if (i > 1 || i < 0) {
      LOGF(fatal, "kFIT Index %d out of range", i);
    }
    auto chPos = ft0Det.getChannelCenter(globalChno);
    auto x = chPos.X() + (*offsetFT0)[i].getX();
    auto y = chPos.Y() + (*offsetFT0)[i].getY();
    auto z = chPos.Z() + (*offsetFT0)[i].getZ();
    if (i == 1) {
      z = -std::abs(z);
    } else if (i == 0) {
      z = std::abs(z);
    }
    auto r = std::sqrt(x * x + y * y);
    auto theta = std::atan2(r, z);
    return -std::log(std::tan(0.5 * theta));
  }

  void loadAlignParam(uint64_t timestamp)
  {
    if (timestamp == mLastTimestamp && offsetFT0 != nullptr) {
      return;
    }
    offsetFT0 = ccdb->getForTimeStamp<std::vector<o2::detectors::AlignParam>>("FT0/Calib/Align", timestamp);
    if (!offsetFT0) {
      LOGF(fatal, "Failed to load valid FT0 alignment from CCDB!");
      return;
    }
    mLastTimestamp = timestamp;
    LOGF(info, "Loaded FT0 alignment for timestamp %llu", timestamp);
  }

  // Two-particle pT correlator from the per-event power sums.
  template <int M, int K>
  std::pair<float, float> calculateMeanAndC2FromSums(const std::array<std::array<double, K>, M>& sumpmwk, const std::array<double, K>& sumwk, float referenceMeanPt) const
  {
    if (sumwk[1] == 0.) {
      return {0.f, 0.f};
    }

    double tau1 = sumwk[2] / (sumwk[1] * sumwk[1]);
    double denom2 = 1. - tau1;

    if (std::abs(denom2) < KFloatEpsilon) {
      double pmk11safe = sumpmwk[1][1] / sumwk[1];
      return {static_cast<float>(pmk11safe), 0.f};
    }

    double pmk11 = sumpmwk[1][1] / sumwk[1];
    double pmk12 = (sumwk[2] != 0.f) ? sumpmwk[1][2] / sumwk[2] : 0.f;
    double pmk22 = (sumwk[2] != 0.f) ? sumpmwk[2][2] / sumwk[2] : 0.f;

    float calculatedMeanPt = pmk11;

    double p1kBar1 = pmk11 - referenceMeanPt;
    double p2kBar2 = pmk22 - 2.0f * pmk12 * referenceMeanPt + referenceMeanPt * referenceMeanPt;

    double p1kBar1sq = p1kBar1 * p1kBar1;
    double numerator2 = p1kBar1sq - (tau1 * p2kBar2);

    float twopcorr = numerator2 / denom2;
    return {calculatedMeanPt, twopcorr};
  }

  // ===========================================================================
  // Table joins
  // ===========================================================================
  using GeneralCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults,
                                      aod::FT0sCorrected,
                                      aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFDDMs, aod::CentFV0As,
                                      aod::CentNTPVs>;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgVtxZCut;
  using AodCollisionsSel = soa::Filtered<GeneralCollisions>;

  using UnfilteredTracks = soa::Join<
    aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>;
  Filter trackFilter = aod::track::pt > KPtMin&& aod::track::pt < KPtMax&& requireGlobalTrackInFilter();
  using AodTracksSel = soa::Filtered<UnfilteredTracks>;
  using TCs = soa::Join<UnfilteredTracks, aod::McTrackLabels>;
  using FilteredTCs = soa::Filtered<TCs>;
  using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

  using MyRun3MCCollisions = soa::Join<
    aod::Collisions, aod::EvSels, aod::Mults, aod::MultsExtra,
    aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFDDMs, aod::CentFV0As,
    aod::CentNGlobals, aod::McCollisionLabels>;

  PresliceUnsorted<MyRun3MCCollisions> colPerMcCollision = aod::mccollisionlabel::mcCollisionId;

  // ===========================================================================
  // Histogram declarations
  // ===========================================================================
  void declareCommonQA()
  {
    histos.add("hVtxZ_after_sel", ";z_{vtx} (cm)", kTH1F, {{KNbinsZvtx, KZvtxMin, KZvtxMax}});
    histos.add("hVtxZ", ";z_{vtx} (cm)", kTH1F, {{KNbinsZvtx, KZvtxMin, KZvtxMax}});
    histos.add("hCentrality", ";centrality (%)", kTH1F, {{centAxis1Per}});
    histos.add("Hist2D_globalTracks_PVTracks", ";N_{global};N_{PV}", kTH2F, {{nChAxis}, {nChAxis}});
    histos.add("Hist2D_cent_nch", ";N_{PV};cent (%)", kTH2F, {{nChAxis}, {centAxis1Per}});
    histos.add("Hist2D_globalTracks_cent", "cent (%);N_{global}", kTH2F, {{centAxis1Per}, {nChAxis}});
    histos.add("Hist2D_PVTracks_cent", "cent (%);N_{PV}", kTH2F, {{centAxis1Per}, {nChAxis}});

    histos.add("hP", ";p (GeV/c)", kTH1F, {{KNbinsPt, KPMin, KPMax}});
    histos.add("hPt", ";p_{T} (GeV/c)", kTH1F, {{KNbinsPt, KPtMin, KPtMax}});
    histos.add("hEta", ";#eta", kTH1F, {{KNbinsEtaFine, KEtaMin, KEtaMax}});
    histos.add("hPhi", ";#phi", kTH1F, {{KNbinsPhi, KPhiMin, TwoPI}});

    histos.add("hEvtCount", "Number of Event;; Count", kTH1F, {{9, 0, 9}});
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(1, "all Events");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(2, "after sel8");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(3, "after VertexZ Cut");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(4, "after kNoSameBunchPileup");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(5, "after kIsGoodZvtxFT0vsPV");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(6, "after kIsGoodITSLayersAll");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(7, "after PVTracksCent Pileup Cut");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(8, "after Linear Pileup Cut (Up)");
    histos.get<TH1>(HIST("hEvtCount"))->GetXaxis()->SetBinLabel(9, "after Linear Pileup Cut (Lw)");

    histos.add("hTrkCount", "Number of Tracks;; Count", kTH1F, {{11, 0, 11}});
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(1, "all Tracks");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(2, "after sign!=0");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(3, "after has_collision");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(4, "after isPVContributor");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(5, "after itsNCls");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(6, "after tpcNClsFound");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(7, "after tpcNClsCrossedRows");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(8, "after pT,#eta selections");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(9, "after isGlobalTrack");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(10, "after dcaXY");
    histos.get<TH1>(HIST("hTrkCount"))->GetXaxis()->SetBinLabel(11, "after dcaZ");
  }

  void declareMCCommonHists()
  {
    histos.add("h3_AllPrimary", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("h3_RecoMatchedToPrimary", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("h3_AllReco", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("h3_RecoUnMatchedToPrimary_Secondary", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("h3_RecoUnMatchedToPrimary_Fake", ";N_{PV};p_{T};#eta", kTH3F, {{nChAxis2}, {KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});

    histos.add("ptResolution", ";p_{T}^{MC};(p_{T}^{reco}-p_{T}^{MC})/p_{T}^{MC}", kTH2F, {{KNbinsPtRes, KPtMin, KPtMax}, {100, -0.2, 0.2}});
    histos.add("etaResolution", ";#eta^{MC};#eta^{reco}-#eta^{MC}", kTH2F, {{KNbinsEtaRes, -KEtaFineMax, KEtaFineMax}, {100, -0.02, 0.02}});
    histos.add("etaTruthReco", ";#eta^{MC};#eta^{reco}", kTH2F, {{KNbinsEtaRes, -KEtaFineMax, KEtaFineMax}, {KNbinsEtaRes, -KEtaFineMax, KEtaFineMax}});
    histos.add("TruthTracKVz", ";Vz^{MC};Vz^{Reco}", kTH2F, {{KNbinsVz, KVzMin, KVzMax}, {KNbinsVz, KVzMin, KVzMax}});
    histos.add("vzResolution", ";Vz^{MC};(Vz^{reco}-Vz^{MC})/Vz^{MC}", kTH2F, {{KNbinsVz, KVzMin, KVzMax}, {100, -0.1, 0.1}});
  }

  void declareMCGetFlatHists()
  {
    histos.add("MCGen/hEtaPhiReco", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});
    histos.add("MCGen/hEtaPhiRecoEffWtd", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});
    histos.add("MCGen/hEtaPhiRecoWtd", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});
  }

  void declareMCMeanHists()
  {
    histos.add("Eff_cent", ";cent", kTProfile, {centAxis1Per});
    histos.add("Eff_Ntrk", ";N_{PV}", kTProfile, {nChAxis2});
    histos.add("Eff_pT", ";p_{T}", kTProfile, {{KNbinsPtRes, KPtMin, KPtMax}});
    histos.add("Eff_eta", ";#eta", kTProfile, {{KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("Fake_cent", ";cent", kTProfile, {centAxis1Per});
    histos.add("Fake_Ntrk", ";N_{PV}", kTProfile, {nChAxis2});
    histos.add("Fake_pT", ";p_{T}", kTProfile, {{KNbinsPtRes, KPtMin, KPtMax}});
    histos.add("Fake_eta", ";#eta", kTProfile, {{KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});
    histos.add("wgt_cent", ";cent", kTProfile, {centAxis1Per});
    histos.add("wgt_Ntrk", ";N_{PV}", kTProfile, {nChAxis2});
    histos.add("wgt_pT", ";p_{T}", kTProfile, {{KNbinsPtRes, KPtMin, KPtMax}});
    histos.add("wgt_eta", ";#eta", kTProfile, {{KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}});

    histos.add("pmeanFT0Amultpv", ";N_{PV};Ampl", kTProfile, {nChAxis});
    histos.add("pmeanFT0Cmultpv", ";N_{PV};Ampl", kTProfile, {nChAxis});
    histos.add("pmeanFT0A_cent", ";cent;Ampl", kTProfile, {centAxis1Per});
    histos.add("pmeanFT0C_cent", ";cent;Ampl", kTProfile, {centAxis1Per});
    histos.add<TProfile3D>("pmean_cent_id_eta_FT0", ";cent;id;#eta", kTProfile3D, {{centAxis1Per}, {200, -0.5, 199.5}, {100, -5.0, 5.0}});
    histos.add("h3_cent_id_eta_FT0", ";cent;id;#eta", kTH3F, {{centAxis1Per}, {200, -0.5, 199.5}, {100, -5.0, 5.0}});

    histos.add<TProfile>("MCGen/Prof_Cent_Nchrec", ";cent;#LT N#GT", kTProfile, {centAxis1Per});
    histos.add<TProfile>("MCGen/Prof_Mult_Nchrec", ";mult;#LT N#GT", kTProfile, {nChAxis});
    histos.add<TProfile>("MCGen/Prof_Cent_MeanpT", ";cent;#LT p_{T}#GT", kTProfile, {centAxis1Per});
    histos.add<TProfile>("MCGen/Prof_Mult_MeanpT", ";mult;#LT p_{T}#GT", kTProfile, {nChAxis});

    histos.add<TProfile2D>("pmeanTru_nch_etabin", ";mult;eta", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("pmeanReco_nch_etabin", ";mult;eta", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("pmeanRecoEffcorr_nch_etabin", ";mult;eta", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("pmeanMultTru_nch_etabin", ";mult;eta", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("pmeanMultReco_nch_etabin", ";mult;eta", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("pmeanMultRecoEffcorr_nch_etabin", ";mult;eta", kTProfile2D, {{nChAxis}, {etaBinAxis}});

    histos.add("MCGen/hEtaPhiReco", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});
    histos.add("MCGen/hEtaPhiRecoEffWtd", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});
    histos.add("MCGen/hEtaPhiRecoWtd", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});

    histos.add<TProfile3D>("Prof2D_MeanpTSub_Tru", ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {etaBinAxis}});
    histos.add<TProfile3D>("Prof2D_MeanpTSub_Reco", ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {etaBinAxis}});
    histos.add<TProfile3D>("Prof2D_MeanpTSub_RecoEffCorr", ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {etaBinAxis}});
  }

  void declareMCFlucHists()
  {
    histos.add<TProfile2D>("MCGen/Prof_Cent_NEta_Nchrec", ";cent;eta", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_Mult_NEta_Nchrec", ";mult;eta", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_Cent_NEta_MeanpT", ";cent;eta", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_Mult_NEta_MeanpT", ";mult;eta", kTProfile2D, {{nChAxis}, {etaBinAxis}});

    histos.add<TProfile2D>("MCGen/Prof_MeanpT_Cent_etabin", ";cent;eta", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_C2_Cent_etabin", ";cent;eta", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_C2Sub_Cent_etabin", ";cent;eta", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_Cov_Cent_etabin", ";cent;eta", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_CovFT0A_Cent_etabin", ";cent;eta", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_CovFT0C_Cent_etabin", ";cent;eta", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});

    histos.add<TProfile2D>("MCGen/Prof_MeanpT_Mult_etabin", ";mult;eta", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_C2_Mult_etabin", ";mult;eta", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_C2Sub_Mult_etabin", ";mult;eta", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_Cov_Mult_etabin", ";mult;eta", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_CovFT0A_Mult_etabin", ";mult;eta", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("MCGen/Prof_CovFT0C_Mult_etabin", ";mult;eta", kTProfile2D, {{nChAxis}, {etaBinAxis}});

    histos.add("MCGen/hEtaPhiReco", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});
    histos.add("MCGen/hEtaPhiRecoEffWtd", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});
    histos.add("MCGen/hEtaPhiRecoWtd", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});

    histos.add<TProfile3D>("MCGen/Prof_C2Sub2D_Cent_etaA_etaC", ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
    histos.add<TProfile3D>("MCGen/Prof_GapSum2D", ";cent;gap;sum", kTProfile3D, {{centAxis1Per}, {gapAxis}, {sumAxis}});
    histos.add<TProfile3D>("MCGen/Prof_Cov2D_Cent_etaA_etaC", ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
    histos.add<TProfile3D>("MCGen/Prof_CovFT0A2D_Cent_etaA_etaC", ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
    histos.add<TProfile3D>("MCGen/Prof_CovFT0C2D_Cent_etaA_etaC", ";cent;etaA;etaC", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
  }

  void declareDataGetFlatHists()
  {
    histos.add("hEtaPhiReco", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});
    histos.add("hEtaPhiRecoEffWtd", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});
    histos.add("hEtaPhiRecoWtd", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});
    histos.add("hnTrkPVZDC", ";N_{PV};ZDC_{A+C}", kTH2F, {{nChAxis2}, {200, 0, 3000}});
    histos.add("hNchZDC", ";N_{trk};ZDC_{A+C}", kTH2F, {{nChAxis2}, {200, 0, 30000}});
  }

  void declareDataMeanHists()
  {
    histos.add("pmeanFT0Amultpv", "N_{PV}; AmplitudeA", kTProfile, {nChAxis});
    histos.add("pmeanFT0A_cent", "cent; AmplitudeA", kTProfile, {centAxis1Per});
    histos.add("pmeanFT0Cmultpv", "N_{PV}; AmplitudeC", kTProfile, {nChAxis});
    histos.add("pmeanFT0C_cent", "cent; AmplitudeC", kTProfile, {centAxis1Per});

    histos.add<TProfile3D>("pmean_cent_id_eta_FT0", ";cent;channel id; #eta;amplitude", kTProfile3D, {{centAxis1Per}, {200, -0.5, 199.5}, {100, -5.0, 5.0}});
    histos.add("h3_cent_id_eta_FT0", ";cent;channel id; #eta", kTH3F, {{centAxis1Per}, {200, -0.5, 199.5}, {100, -5.0, 5.0}});

    histos.add<TProfile>("Prof_Cent_Nchrec", ";cent;#LT N_{PV}#GT", kTProfile, {centAxis1Per});
    histos.add<TProfile>("Prof_Mult_Nchrec", ";N_{PV};#LT N_{PV}#GT", kTProfile, {nChAxis});
    histos.add<TProfile>("Prof_Cent_MeanpT", ";cent;#LT p_{T}#GT", kTProfile, {centAxis1Per});
    histos.add<TProfile>("Prof_Mult_MeanpT", ";N_{PV};#LT p_{T}#GT", kTProfile, {nChAxis});

    histos.add<TProfile2D>("pmean_nch_etabin", ";N_{PV};#eta-bin", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("pmeanMult_nch_etabin", ";N_{PV};#eta-bin", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("pmean_cent_etabin", ";Centrality (%);#eta-bin", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
    histos.add<TProfile2D>("pmeanMult_cent_etabin", ";Centrality (%);#eta-bin", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});

    histos.add("hEtaPhiReco", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});
    histos.add("hEtaPhiRecoEffWtd", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});
    histos.add("hEtaPhiRecoWtd", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});

    histos.add<TProfile3D>("Prof2D_MeanpTSub", ";cent;#eta_{A} bin;#eta_{C} bin", kTProfile3D, {{centAxis1Per}, {etaBinAxis}, {etaBinAxis}});

    histos.add<TProfile3D>("pEffWeight_pt_eta_cent", ";p_{T} (GeV/c);#eta;cent;#LT eff#GT", kTProfile3D, {{KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}, {centAxis1Per}});
    histos.add<TProfile3D>("pFakeWeight_pt_eta_cent", ";p_{T} (GeV/c);#eta;cent;#LT fake#GT", kTProfile3D, {{KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}, {centAxis1Per}});
    histos.add<TProfile3D>("pFlatWeight_pt_eta_cent", ";p_{T} (GeV/c);#eta;cent;#LT w_{#phi}#GT", kTProfile3D, {{KNbinsPtRes, KPtMin, KPtMax}, {KNbinsEtaFine, -KEtaFineMax, KEtaFineMax}, {centAxis1Per}});
  }

  void declareDataFlucHists()
  {
    histos.add<TProfile2D>("Prof_MeanpT_Cent_etabin", ";cent;#eta-bin", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
    histos.add<TProfile2D>("Prof_MeanpT_Mult_etabin", ";N_{PV};#eta-bin", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("Prof_C2_Cent_etabin", ";cent;#eta-bin", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
    histos.add<TProfile2D>("Prof_C2_Mult_etabin", ";N_{PV};#eta-bin", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("Prof_C2Sub_Cent_etabin", ";Centrality;#eta-bin", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
    histos.add<TProfile2D>("Prof_C2Sub_Mult_etabin", ";N_{PV};#eta-bin", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("Prof_Cov_Cent_etabin", ";Centrality;#eta-bin", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
    histos.add<TProfile2D>("Prof_Cov_Mult_etabin", ";N_{PV};#eta-bin", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("Prof_CovFT0A_Cent_etabin", ";Centrality;#eta-bin", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
    histos.add<TProfile2D>("Prof_CovFT0A_Mult_etabin", ";N_{PV};#eta-bin", kTProfile2D, {{nChAxis}, {etaBinAxis}});
    histos.add<TProfile2D>("Prof_CovFT0C_Cent_etabin", ";Centrality;#eta-bin", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
    histos.add<TProfile2D>("Prof_CovFT0C_Mult_etabin", ";N_{PV};#eta-bin", kTProfile2D, {{nChAxis}, {etaBinAxis}});

    histos.add("hEtaPhiReco", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});
    histos.add("hEtaPhiRecoEffWtd", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});
    histos.add("hEtaPhiRecoWtd", ";vz;sign;pt;eta;phi", kTHnSparseF, {{vzAxis}, {chgAxis}, {pTAxis}, {etaFlatAxis}, {phiAxis}});

    histos.add<TProfile3D>("Prof_C2Sub2D_Cent_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
    histos.add<TProfile3D>("Prof_GapSum2D", ";cent;#Delta#eta (Gap);#Sigma#eta/2 (Sum)", kTProfile3D, {{centAxis1Per}, {gapAxis}, {sumAxis}});
    histos.add<TProfile3D>("Prof_Cov2D_Cent_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
    histos.add<TProfile3D>("Prof_CovFT0A2D_Cent_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
    histos.add<TProfile3D>("Prof_CovFT0C2D_Cent_etaA_etaC", ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
  }

  // 30 Poisson-bootstrap replicas of every final fluctuation observable
  // (base data run only). Filled by pointer to bypass the compile-time HIST()
  // macro, which cannot take a runtime sample index.
  void declareBootstrapHists()
  {
    for (int s = 0; s < nBoot; ++s) {
      bs.meanpTCent[s] = histos.add<TProfile2D>(Form("Bootstrap/Prof_MeanpT_Cent_etabin_sample%d", s), ";cent;#eta-bin", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
      bs.meanpTMult[s] = histos.add<TProfile2D>(Form("Bootstrap/Prof_MeanpT_Mult_etabin_sample%d", s), ";N_{PV};#eta-bin", kTProfile2D, {{nChAxis}, {etaBinAxis}});
      bs.c2Cent[s] = histos.add<TProfile2D>(Form("Bootstrap/Prof_C2_Cent_etabin_sample%d", s), ";cent;#eta-bin", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
      bs.c2Mult[s] = histos.add<TProfile2D>(Form("Bootstrap/Prof_C2_Mult_etabin_sample%d", s), ";N_{PV};#eta-bin", kTProfile2D, {{nChAxis}, {etaBinAxis}});
      bs.c2SubCent[s] = histos.add<TProfile2D>(Form("Bootstrap/Prof_C2Sub_Cent_etabin_sample%d", s), ";cent;#eta-bin", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
      bs.c2SubMult[s] = histos.add<TProfile2D>(Form("Bootstrap/Prof_C2Sub_Mult_etabin_sample%d", s), ";N_{PV};#eta-bin", kTProfile2D, {{nChAxis}, {etaBinAxis}});
      bs.covCent[s] = histos.add<TProfile2D>(Form("Bootstrap/Prof_Cov_Cent_etabin_sample%d", s), ";cent;#eta-bin", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
      bs.covMult[s] = histos.add<TProfile2D>(Form("Bootstrap/Prof_Cov_Mult_etabin_sample%d", s), ";N_{PV};#eta-bin", kTProfile2D, {{nChAxis}, {etaBinAxis}});
      bs.covFT0ACent[s] = histos.add<TProfile2D>(Form("Bootstrap/Prof_CovFT0A_Cent_etabin_sample%d", s), ";cent;#eta-bin", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
      bs.covFT0AMult[s] = histos.add<TProfile2D>(Form("Bootstrap/Prof_CovFT0A_Mult_etabin_sample%d", s), ";N_{PV};#eta-bin", kTProfile2D, {{nChAxis}, {etaBinAxis}});
      bs.covFT0CCent[s] = histos.add<TProfile2D>(Form("Bootstrap/Prof_CovFT0C_Cent_etabin_sample%d", s), ";cent;#eta-bin", kTProfile2D, {{centAxis1Per}, {etaBinAxis}});
      bs.covFT0CMult[s] = histos.add<TProfile2D>(Form("Bootstrap/Prof_CovFT0C_Mult_etabin_sample%d", s), ";N_{PV};#eta-bin", kTProfile2D, {{nChAxis}, {etaBinAxis}});
      bs.c2Sub2D[s] = histos.add<TProfile3D>(Form("Bootstrap/Prof_C2Sub2D_Cent_etaA_etaC_sample%d", s), ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
      bs.gapSum2D[s] = histos.add<TProfile3D>(Form("Bootstrap/Prof_GapSum2D_sample%d", s), ";cent;gap;sum", kTProfile3D, {{centAxis1Per}, {gapAxis}, {sumAxis}});
      bs.cov2D[s] = histos.add<TProfile3D>(Form("Bootstrap/Prof_Cov2D_Cent_etaA_etaC_sample%d", s), ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
      bs.covFT0A2D[s] = histos.add<TProfile3D>(Form("Bootstrap/Prof_CovFT0A2D_Cent_etaA_etaC_sample%d", s), ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
      bs.covFT0C2D[s] = histos.add<TProfile3D>(Form("Bootstrap/Prof_CovFT0C2D_Cent_etaA_etaC_sample%d", s), ";cent;#eta_{A};#eta_{C}", kTProfile3D, {{centAxis1Per}, {etaAxis}, {etaAxis}});
    }
  }

  // ===========================================================================
  // CCDB helpers
  // ===========================================================================
  THnSparseF* buildWeightMapFromRaw(THnSparseF* hRaw, const char* mapName)
  {
    if (!hRaw) {
      LOGF(error, "Raw eta-phi map for '%s' is null; no flattening will be applied.", mapName);
      return nullptr;
    }
    auto hWMap = dynamic_cast<THnSparseF*>(hRaw->Clone(mapName));
    hWMap->SetTitle(Form("Flattening Weight Map %s (w_{#phi} = <N_{#phi}> / N_{#phi})", mapName));
    hWMap->Reset();
    auto axV = hRaw->GetAxis(0);
    auto axChg = hRaw->GetAxis(1);
    auto axPt = hRaw->GetAxis(2);
    auto axE = hRaw->GetAxis(3);
    auto axP = hRaw->GetAxis(4);

    std::array<int, 5> bins{};
    for (int iv = 1; iv <= axV->GetNbins(); ++iv) {
      bins[0] = iv;
      for (int ichg = 1; ichg <= axChg->GetNbins(); ++ichg) {
        bins[1] = ichg;
        for (int ipt = 1; ipt <= axPt->GetNbins(); ++ipt) {
          bins[2] = ipt;
          for (int ie = 1; ie <= axE->GetNbins(); ++ie) {
            bins[3] = ie;
            double sum = 0.0;
            int nphi = axP->GetNbins();
            for (int ip = 1; ip <= nphi; ++ip) {
              bins[4] = ip;
              sum += hRaw->GetBinContent(bins.data());
            }
            const double avg = (nphi > 0 ? sum / nphi : 0.0);
            for (int ip = 1; ip <= nphi; ++ip) {
              bins[4] = ip;
              const double raw = hRaw->GetBinContent(bins.data());
              const double w = (avg > 0.0 && raw > 0.0) ? (avg / raw) : 1.0;
              hWMap->SetBinContent(bins.data(), w);
            }
          }
        }
      }
    }
    LOGF(info, "Flattening weight map '%s' built.", mapName);
    return hWMap;
  }

  template <typename TP>
  void loadProfileFromList(TList* src, const char* name, TP*& target)
  {
    if (!src) {
      return;
    }
    auto* obj = src->FindObject(name);
    if (!obj) {
      LOGF(error, "Profile %s missing in CCDB TList", name);
      return;
    }
    auto* tp = dynamic_cast<TP*>(obj);
    if (!tp) {
      LOGF(error, "%s is not the expected profile type (it is %s)", name, obj->ClassName());
      return;
    }
    target = dynamic_cast<TP*>(tp->Clone());
    target->SetDirectory(nullptr);
    LOGF(info, "Loaded %s from list", name);
  }

  // ===========================================================================
  // init
  // ===========================================================================
  void init(InitContext&)
  {
    // Nch axes by system
    if (cfgSys == kPbPb) {
      nChAxis = {cfgNchPbMax / 2, KBinOffset, cfgNchPbMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
      nChAxis2 = {cfgNchPbMax / 4, KBinOffset, cfgNchPbMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
    } else {
      nChAxis = {cfgNchOMax, KBinOffset, cfgNchOMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
      nChAxis2 = {cfgNchOMax, KBinOffset, cfgNchOMax + KBinOffset, "Nch", "PV-contributor track multiplicity"};
    }

    // ---- observable eta binning: 0.1 for the base data run, 0.2 otherwise ----
    const bool isDataRun = (cfgRunGetDataFlat || cfgRunDataMean || cfgRunDataFluc);
    const bool useFineBinning = isDataRun && (cfgSystType == kSystBase);
    {
      const float lo = -cfgCutEta;
      const float hi = cfgCutEta;
      const float width = useFineBinning ? 0.1f : 0.2f;
      const int nbins = static_cast<int>(std::lround((hi - lo) / width));
      etaLw.clear();
      etaUp.clear();
      etaLw.push_back(lo); // index 0: full-range reference bin
      etaUp.push_back(hi);
      for (int i = 0; i < nbins; ++i) {
        etaLw.push_back(lo + i * width);
        etaUp.push_back(lo + (i + 1) * width);
      }
      nEta = nbins + 1;
      if (nEta > KNEtaMax) {
        LOGF(fatal, "nEta=%d exceeds KNEtaMax=%d", nEta, KNEtaMax);
      }

      std::vector<double> obsEdges;
      obsEdges.reserve(nEta);
      obsEdges.push_back(etaLw[1]);
      for (int i = 1; i < nEta; ++i) {
        obsEdges.push_back(etaUp[i]);
      }
      etaAxis = AxisSpec{obsEdges, "#eta"};
      etaBinAxis = AxisSpec{nEta + 1, -0.5, static_cast<double>(nEta) + 0.5, "#eta bin Number"};
      LOGF(info, "Observable eta binning: %d bins of width %.2f (+ reference), nEta=%d", nbins, width, nEta);
    }

    // bootstrap active only for the base data fluctuation pass
    doBoot = cfgRunDataFluc && (cfgSystType == kSystBase) && (cfgNBootstrap > 0);
    nBoot = std::min<int>(cfgNBootstrap, KMaxBoot);
    rng.SetSeed(cfgBootstrapSeed);

    // ---- CCDB ----
    ccdb->setURL(cfgCCDBurl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    loadAlignParam(now);
    ft0Det.calculateChannelCenter();

    std::string sysDir;
    switch (cfgSys) {
      case kPbPb:
        sysDir = "PbPbTest";
        break;
      case kNeNe:
        sysDir = "NeNeTest";
        break;
      case kOO:
        sysDir = "OOTest";
        break;
      case kpp:
        sysDir = "ppTest";
        break;
      default:
        LOGF(fatal, "Invalid cfgSys value: %d", cfgSys.value);
    }
    std::string pathEff = cfgCCDBUserPath.value + "/" + sysDir + "/Job1_EffMaps";
    std::string pathMCFlat = cfgCCDBUserPath.value + "/" + sysDir + "/Job1_MCFlatMaps";
    std::string pathMCMean = cfgCCDBUserPath.value + "/" + sysDir + "/Job2_MCMean";
    std::string pathDataFlat = cfgCCDBUserPath.value + "/" + sysDir + "/Job1_DataFlatMaps";
    std::string pathDataMean = cfgCCDBUserPath.value + "/" + sysDir + "/Job2_DataMean";

    // ---- declarations ----
    declareCommonQA();
    if (cfgRunMCMean || cfgRunMCFluc || cfgRunGetEff) {
      declareMCCommonHists();
    }
    if (cfgRunGetMCFlat) {
      declareMCGetFlatHists();
      histos.addClone("MCGen/", "MCReco/");
      histos.addClone("MCGen/", "MCRecoEffCorr/");
    }
    if (cfgRunMCMean) {
      declareMCMeanHists();
      histos.addClone("MCGen/", "MCReco/");
      histos.addClone("MCGen/", "MCRecoEffCorr/");
    }
    if (cfgRunMCFluc) {
      declareMCFlucHists();
      histos.addClone("MCGen/", "MCReco/");
      histos.addClone("MCGen/", "MCRecoEffCorr/");
    }
    if (cfgRunGetDataFlat) {
      declareDataGetFlatHists();
    }
    if (cfgRunDataMean) {
      declareDataMeanHists();
    }
    if (cfgRunDataFluc) {
      declareDataFlucHists();
      if (doBoot) {
        declareBootstrapHists();
      }
    }

    // ---- efficiency/fake maps (inclusive) ----
    if (!cfgRunGetEff && cfgEff) {
      auto* lst = ccdb->getForTimeStamp<TList>(pathEff, now);
      if (!lst) {
        LOGF(fatal, "Efficiency maps required but CCDB list is null at %s!", pathEff.c_str());
        return;
      }

      auto* hNum = dynamic_cast<TH3F*>(lst->FindObject("h3_RecoMatchedToPrimary"));
      auto* hDen = dynamic_cast<TH3F*>(lst->FindObject("h3_AllPrimary"));
      if (hNum && hDen) {
        state.hEff = dynamic_cast<TH3F*>(hNum->Clone("hEff"));
        state.hEff->SetDirectory(nullptr);
        state.hEff->Divide(hDen);
      } else {
        LOGF(error, "Missing CCDB objects for efficiency (h3_RecoMatchedToPrimary / h3_AllPrimary).");
      }

      auto* hNumS = dynamic_cast<TH3F*>(lst->FindObject("h3_RecoUnMatchedToPrimary_Secondary"));
      auto* hNumF = dynamic_cast<TH3F*>(lst->FindObject("h3_RecoUnMatchedToPrimary_Fake"));
      auto* hDenF = dynamic_cast<TH3F*>(lst->FindObject("h3_AllReco"));
      if (hNumS && hNumF && hDenF) {
        state.hFake = dynamic_cast<TH3F*>(hNumS->Clone("hFake"));
        state.hFake->Add(hNumF);
        state.hFake->SetDirectory(nullptr);
        state.hFake->Divide(hDenF);
      } else {
        LOGF(error, "Missing CCDB objects for fakes.");
      }
    }

    // ---- flattening maps (inclusive) ----
    if (!cfgRunGetEff && cfgFlat) {
      if (cfgRunDataMean || cfgRunDataFluc) {
        auto* lstDataFlat = ccdb->getForTimeStamp<TList>(pathDataFlat, now);
        if (lstDataFlat) {
          std::string hName = cfgEff ? "hEtaPhiRecoWtd" : "hEtaPhiReco";
          auto* hRaw = dynamic_cast<THnSparseF*>(lstDataFlat->FindObject(hName.c_str()));
          if (hRaw) {
            state.hFlatWeight = buildWeightMapFromRaw(hRaw, "hFlatWeight");
          } else {
            LOGF(error, "Data flattening map '%s' not found.", hName.c_str());
          }
        } else {
          LOGF(error, "Could not retrieve Data Flattening TList from: %s", pathDataFlat.c_str());
        }
      }
      if (cfgRunMCMean || cfgRunMCFluc) {
        auto* lstMCFlat = ccdb->getForTimeStamp<TList>(pathMCFlat, now);
        if (lstMCFlat) {
          std::string hName = cfgEff ? "MCReco/hEtaPhiRecoWtd" : "MCReco/hEtaPhiReco";
          auto* hRaw = dynamic_cast<THnSparseF*>(lstMCFlat->FindObject(hName.c_str()));
          if (hRaw) {
            state.hFlatWeight = buildWeightMapFromRaw(hRaw, "hFlatWeight");
          } else {
            LOGF(warning, "MC flattening source '%s' not found.", hName.c_str());
          }
        } else {
          LOGF(error, "Could not retrieve MC Flattening TList from: %s", pathMCFlat.c_str());
        }
      }
    }

    // ---- sigma-pileup limits ----
    // These were produced by the (now removed) nSigma pass. They are re-sourced
    // from the flattening-map file, which also stores Hist2D_globalTracks_cent.
    // Loaded only when the sigma-pileup cut is actually requested.
    if (cfgApplySigPupCut) {
      const bool mcSide = (cfgRunGetEff || cfgRunGetMCFlat || cfgRunMCMean || cfgRunMCFluc);
      std::string limPath = mcSide ? pathMCFlat : pathDataFlat;
      auto* limList = ccdb->getForTimeStamp<TList>(limPath, now);
      if (limList) {
        auto loadLimits = [&](const char* name, std::vector<std::pair<float, float>>& limits, float& xMin, float& xMax) {
          auto* h2 = dynamic_cast<TH2*>(limList->FindObject(name));
          if (!h2) {
            return;
          }
          std::unique_ptr<TProfile> prof(h2->ProfileX("ptmp", 1, -1, "S"));
          int nBins = prof->GetNbinsX();
          xMin = prof->GetXaxis()->GetXmin();
          xMax = prof->GetXaxis()->GetXmax();
          limits.assign(nBins + 2, {-99999.f, 999999.f});
          for (int i = 1; i <= nBins; ++i) {
            float mean = prof->GetBinContent(i);
            float rms = prof->GetBinError(i);
            limits[i] = {mean - cfgPupnSig * rms, mean + cfgPupnSig * rms};
          }
        };
        loadLimits("Hist2D_globalTracks_cent", state.mLimitsNchCent, state.mMinXNchCent, state.mMaxXNchCent);
      } else {
        LOGF(warning, "sigma-pileup limits source list missing at %s; sigma-pileup cut effectively disabled.", limPath.c_str());
      }
    }

    // ---- MC mean profiles for MC fluc ----
    if (cfgRunMCFluc) {
      LOGF(info, "Loading MC Mean profiles from CCDB path: %s", pathMCMean.c_str());
      auto* lstMCMean = ccdb->getForTimeStamp<TList>(pathMCMean, now);
      if (lstMCMean) {
        loadProfileFromList(lstMCMean, "pmeanFT0Amultpv", state.pmeanFT0AmultpvStep2);
        loadProfileFromList(lstMCMean, "pmeanFT0Cmultpv", state.pmeanFT0CmultpvStep2);
        loadProfileFromList(lstMCMean, "pmeanTru_nch_etabin", state.pmeanTruNchEtabinStep2);
        loadProfileFromList(lstMCMean, "pmeanReco_nch_etabin", state.pmeanRecoNchEtabinStep2);
        loadProfileFromList(lstMCMean, "pmeanRecoEffcorr_nch_etabin", state.pmeanRecoEffcorrNchEtabinStep2);
        loadProfileFromList(lstMCMean, "pmeanMultTru_nch_etabin", state.pmeanMultTruNchEtabinStep2);
        loadProfileFromList(lstMCMean, "pmeanMultReco_nch_etabin", state.pmeanMultRecoNchEtabinStep2);
        loadProfileFromList(lstMCMean, "pmeanMultRecoEffcorr_nch_etabin", state.pmeanMultRecoEffcorrNchEtabinStep2);
      } else {
        LOGF(error, "Could not retrieve TList for MC Mean from: %s", pathMCMean.c_str());
      }
    }

    // ---- Data mean profiles for data fluc (one per systematic variation) ----
    if (cfgRunDataFluc) {
      int st = std::clamp<int>(cfgSystType, 0, kNSystType - 1);
      std::string meanPath = pathDataMean + systSuffix[st];
      LOGF(info, "Loading Data Mean profiles for systematic '%s' from: %s",
           (st == kSystBase ? "Base" : systSuffix[st].c_str()), meanPath.c_str());
      auto* lstDataMean = ccdb->getForTimeStamp<TList>(meanPath, now);
      if (lstDataMean) {
        loadProfileFromList(lstDataMean, "pmeanFT0Amultpv", state.pmeanFT0AmultpvStep2);
        loadProfileFromList(lstDataMean, "pmeanFT0Cmultpv", state.pmeanFT0CmultpvStep2);
        loadProfileFromList(lstDataMean, "pmean_nch_etabin", state.pmeanNchEtabinStep2);
        loadProfileFromList(lstDataMean, "pmeanMult_nch_etabin", state.pmeanMultNchEtabinStep2);
      } else {
        LOGF(error, "Could not retrieve TList for Data Mean from: %s", meanPath.c_str());
      }
    }
    LOGF(info, "CCDB initialization complete for RadialFlowDecorr.");
  }

  // ===========================================================================
  // MC: build efficiency / fake maps (inclusive)
  // ===========================================================================
  void processGetEffHists(MyRun3MCCollisions::iterator const& mcCollision, FilteredTCs const& mcTracks, aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("hVtxZ"), mcCollision.posZ());
    if (!mcCollision.has_mcCollision() || !isEventSelected(mcCollision)) {
      return;
    }
    float cent = getCentrality(mcCollision);
    if (cent > KCentMax) {
      return;
    }
    float multPV = mcCollision.multNTracksPV();
    float vz = mcCollision.posZ();
    if (!isPassAddPileup(multPV, mcTracks.size(), cent)) {
      return;
    }

    histos.fill(HIST("hVtxZ_after_sel"), mcCollision.posZ());
    histos.fill(HIST("hCentrality"), cent);
    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, mcTracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), mcTracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, mcTracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, multPV);

    for (const auto& particle : mcParticles) {
      if (particle.mcCollisionId() != mcCollision.mcCollisionId()) {
        continue;
      }
      if (!isParticleSelected(particle) || !particle.isPhysicalPrimary()) {
        continue;
      }
      histos.fill(HIST("h3_AllPrimary"), multPV, particle.pt(), particle.eta());
    }

    for (const auto& track : mcTracks) {
      if (track.collisionId() != mcCollision.index()) {
        continue;
      }
      if (!isTrackSelected(track)) {
        continue;
      }

      float pt = track.pt(), eta = track.eta();
      histos.fill(HIST("hPt"), pt);
      histos.fill(HIST("hEta"), eta);
      histos.fill(HIST("h3_AllReco"), multPV, pt, eta);

      if (track.has_mcParticle()) {
        auto mcP = track.mcParticle();
        if (mcP.isPhysicalPrimary()) {
          histos.fill(HIST("ptResolution"), mcP.pt(), (pt - mcP.pt()) / mcP.pt());
          histos.fill(HIST("etaResolution"), mcP.eta(), eta - mcP.eta());
          histos.fill(HIST("etaTruthReco"), mcP.eta(), eta);
          histos.fill(HIST("vzResolution"), mcP.vz(), (vz - mcP.vz()) / mcP.vz());
          histos.fill(HIST("TruthTracKVz"), mcP.vz(), vz);
          histos.fill(HIST("h3_RecoMatchedToPrimary"), multPV, mcP.pt(), mcP.eta());
        } else {
          histos.fill(HIST("h3_RecoUnMatchedToPrimary_Secondary"), multPV, pt, eta);
        }
      } else {
        histos.fill(HIST("h3_RecoUnMatchedToPrimary_Fake"), multPV, pt, eta);
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processGetEffHists, "process MC to calculate EffWeights", cfgRunGetEff);

  // ===========================================================================
  // MC: build flattening maps (inclusive)
  // ===========================================================================
  void processMCFlat(MyRun3MCCollisions::iterator const& mcCollision, FilteredTCs const& mcTracks)
  {
    histos.fill(HIST("hVtxZ"), mcCollision.posZ());
    if (!mcCollision.has_mcCollision() || !isEventSelected(mcCollision)) {
      return;
    }
    float cent = getCentrality(mcCollision);
    if (cent > KCentMax) {
      return;
    }
    float multPV = mcCollision.multNTracksPV();
    float vz = mcCollision.posZ();
    if (!isPassAddPileup(multPV, mcTracks.size(), cent)) {
      return;
    }

    histos.fill(HIST("hVtxZ_after_sel"), mcCollision.posZ());
    histos.fill(HIST("hCentrality"), cent);
    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, mcTracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), mcTracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, mcTracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, multPV);

    for (const auto& track : mcTracks) {
      if (track.collisionId() != mcCollision.index()) {
        continue;
      }
      if (!isTrackSelected(track)) {
        continue;
      }

      float pt = track.pt(), eta = track.eta(), phi = track.phi();
      auto sign = track.sign();
      histos.fill(HIST("hPt"), pt);
      histos.fill(HIST("hEta"), eta);
      histos.fill(HIST("hPhi"), phi);

      float eff = getEfficiency(multPV, pt, eta, 0, cfgEff);
      float fake = getEfficiency(multPV, pt, eta, 1, cfgEff);
      float w = (eff > KFloatEpsilon) ? (1.0f - fake) / eff : 0.0f;
      if (std::isfinite(w) && w > 0.f) {
        histos.fill(HIST("MCReco/hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, w);
        histos.fill(HIST("MCReco/hEtaPhiReco"), vz, sign, pt, eta, phi, 1.0);
        histos.fill(HIST("MCReco/hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCFlat, "process MC to calculate FlatWeights", cfgRunGetMCFlat);

  // ===========================================================================
  // MC: mean pT (truth / reco / reco-eff-corrected)
  // ===========================================================================
  void processMCMean(MyRun3MCCollisions::iterator const& mcCollision, FilteredTCs const& mcTracks, aod::FT0s const&, aod::McParticles const& mcParticles)
  {
    std::array<double, KNEtaMax> sumWiTruth{}, sumWiptiTruth{};
    std::array<double, KNEtaMax> sumWiReco{}, sumWiptiReco{};
    std::array<double, KNEtaMax> sumWiRecoEffCorr{}, sumWiptiRecoEffCorr{};

    histos.fill(HIST("hVtxZ"), mcCollision.posZ());
    if (!mcCollision.has_mcCollision() || !isEventSelected(mcCollision)) {
      return;
    }
    float cent = getCentrality(mcCollision);
    if (cent > KCentMax) {
      return;
    }
    float multPV = mcCollision.multNTracksPV();
    float vz = mcCollision.posZ();
    if (!isPassAddPileup(multPV, mcTracks.size(), cent)) {
      return;
    }

    histos.fill(HIST("hVtxZ_after_sel"), mcCollision.posZ());
    histos.fill(HIST("hCentrality"), cent);
    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, mcTracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), mcTracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, mcTracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, multPV);

    // --- truth ---
    for (const auto& particle : mcParticles) {
      if (particle.mcCollisionId() != mcCollision.mcCollisionId()) {
        continue;
      }
      if (!isParticleSelected(particle) || !particle.isPhysicalPrimary()) {
        continue;
      }
      float pt = particle.pt(), eta = particle.eta();
      if (pt <= cfgPtMin || pt > cfgPtMax) {
        continue;
      }
      for (int ieta = 0; ieta < nEta; ++ieta) {
        if (eta <= etaLw[ieta] || eta > etaUp[ieta]) {
          continue;
        }
        sumWiTruth[ieta]++;
        sumWiptiTruth[ieta] += pt;
      }
    }

    histos.fill(HIST("MCGen/Prof_Cent_Nchrec"), cent, sumWiTruth[0]);
    histos.fill(HIST("MCGen/Prof_Mult_Nchrec"), multPV, sumWiTruth[0]);
    if (sumWiTruth[0] > 1.0f) {
      histos.fill(HIST("MCGen/Prof_Cent_MeanpT"), cent, sumWiptiTruth[0] / sumWiTruth[0]);
      histos.fill(HIST("MCGen/Prof_Mult_MeanpT"), multPV, sumWiptiTruth[0] / sumWiTruth[0]);
    }

    // --- reco ---
    for (const auto& track : mcTracks) {
      if (track.collisionId() != mcCollision.index()) {
        continue;
      }
      if (!isTrackSelected(track)) {
        continue;
      }
      float pt = track.pt(), eta = track.eta(), phi = track.phi();
      if (pt <= cfgPtMin || pt > cfgPtMax) {
        continue;
      }
      auto sign = track.sign();
      histos.fill(HIST("hPt"), pt);
      histos.fill(HIST("hEta"), eta);
      histos.fill(HIST("hPhi"), phi);

      float eff = getEfficiency(multPV, pt, eta, 0, cfgEff);
      float fake = getEfficiency(multPV, pt, eta, 1, cfgEff);
      float flatW = getFlatteningWeight(vz, sign, pt, eta, phi, cfgFlat);
      float w = flatW * (1.0 - fake) / eff;
      if (!std::isfinite(w) || w <= 0.f || eff <= KFloatEpsilon) {
        continue;
      }

      for (int ieta = 0; ieta < nEta; ++ieta) {
        if (eta <= etaLw[ieta] || eta > etaUp[ieta]) {
          continue;
        }
        sumWiReco[ieta]++;
        sumWiptiReco[ieta] += pt;
        sumWiRecoEffCorr[ieta] += w;
        sumWiptiRecoEffCorr[ieta] += w * pt;
      }

      histos.fill(HIST("Eff_cent"), cent, eff);
      histos.fill(HIST("Fake_cent"), cent, fake);
      histos.fill(HIST("wgt_cent"), cent, w);
      histos.fill(HIST("Eff_Ntrk"), multPV, eff);
      histos.fill(HIST("Fake_Ntrk"), multPV, fake);
      histos.fill(HIST("wgt_Ntrk"), multPV, w);
      histos.fill(HIST("Eff_pT"), pt, eff);
      histos.fill(HIST("Fake_pT"), pt, fake);
      histos.fill(HIST("wgt_pT"), pt, w);
      histos.fill(HIST("Eff_eta"), eta, eff);
      histos.fill(HIST("Fake_eta"), eta, fake);
      histos.fill(HIST("wgt_eta"), eta, w);

      histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
      histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
      histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
    }

    if (!hasMinTracksInAllEtaBins(sumWiTruth) || !hasMinTracksInAllEtaBins(sumWiReco)) {
      return;
    }

    // subevent mean-pT maps
    for (int ietaA = 0; ietaA < nEta; ++ietaA) {
      for (int ietaC = 0; ietaC < nEta; ++ietaC) {
        float nTruAB = sumWiTruth[ietaA] + sumWiTruth[ietaC];
        float nRecoAB = sumWiReco[ietaA] + sumWiReco[ietaC];
        float nCorrAB = sumWiRecoEffCorr[ietaA] + sumWiRecoEffCorr[ietaC];

        if (nTruAB > 0) {
          histos.fill(HIST("Prof2D_MeanpTSub_Tru"), cent, ietaA, ietaC, (sumWiptiTruth[ietaA] + sumWiptiTruth[ietaC]) / nTruAB);
        }
        if (nRecoAB > 0) {
          histos.fill(HIST("Prof2D_MeanpTSub_Reco"), cent, ietaA, ietaC, (sumWiptiReco[ietaA] + sumWiptiReco[ietaC]) / nRecoAB);
        }
        if (nCorrAB > 0) {
          histos.fill(HIST("Prof2D_MeanpTSub_RecoEffCorr"), cent, ietaA, ietaC, (sumWiptiRecoEffCorr[ietaA] + sumWiptiRecoEffCorr[ietaC]) / nCorrAB);
        }
      }

      if (sumWiTruth[ietaA] > 0) {
        histos.fill(HIST("pmeanTru_nch_etabin"), multPV, ietaA, sumWiptiTruth[ietaA] / sumWiTruth[ietaA]);
        histos.fill(HIST("pmeanMultTru_nch_etabin"), multPV, ietaA, sumWiTruth[ietaA]);
      }
      if (sumWiReco[ietaA] > 0) {
        histos.fill(HIST("pmeanReco_nch_etabin"), multPV, ietaA, sumWiptiReco[ietaA] / sumWiReco[ietaA]);
        histos.fill(HIST("pmeanMultReco_nch_etabin"), multPV, ietaA, sumWiReco[ietaA]);
      }
      if (sumWiRecoEffCorr[ietaA] > 0) {
        histos.fill(HIST("pmeanRecoEffcorr_nch_etabin"), multPV, ietaA, sumWiptiRecoEffCorr[ietaA] / sumWiRecoEffCorr[ietaA]);
        histos.fill(HIST("pmeanMultRecoEffcorr_nch_etabin"), multPV, ietaA, sumWiRecoEffCorr[ietaA]);
      }
    }

    // FT0
    double amplFT0A = 0, amplFT0C = 0;
    if (mcCollision.has_foundFT0()) {
      const auto& ft0 = mcCollision.foundFT0();
      for (std::size_t iCh = 0; iCh < ft0.channelA().size(); iCh++) {
        auto chanelid = ft0.channelA()[iCh];
        float ampl = ft0.amplitudeA()[iCh];
        amplFT0A += ampl;
        auto eta = getEtaFT0(chanelid, 0);
        histos.fill(HIST("pmean_cent_id_eta_FT0"), cent, chanelid, eta, ampl);
        histos.fill(HIST("h3_cent_id_eta_FT0"), cent, chanelid, eta, ampl);
      }
      for (std::size_t iCh = 0; iCh < ft0.channelC().size(); iCh++) {
        auto chanelid = ft0.channelC()[iCh];
        auto globalId = chanelid + KnFt0cCell;
        float ampl = ft0.amplitudeC()[iCh];
        auto eta = getEtaFT0(globalId, 1);
        amplFT0C += ampl;
        histos.fill(HIST("pmean_cent_id_eta_FT0"), cent, globalId, eta, ampl);
        histos.fill(HIST("h3_cent_id_eta_FT0"), cent, globalId, eta, ampl);
      }
    }
    histos.fill(HIST("pmeanFT0Amultpv"), multPV, amplFT0A);
    histos.fill(HIST("pmeanFT0A_cent"), cent, amplFT0A);
    histos.fill(HIST("pmeanFT0Cmultpv"), multPV, amplFT0C);
    histos.fill(HIST("pmeanFT0C_cent"), cent, amplFT0C);
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCMean, "process MC to calculate mean pt", cfgRunMCMean);

  // ===========================================================================
  // MC: fluctuations (C2, subevent) at three levels
  // ===========================================================================
  void processMCFluc(MyRun3MCCollisions::iterator const& mcCollision, FilteredTCs const& mcTracks, aod::FT0s const&, aod::McParticles const& mcParticles)
  {
    if (!state.pmeanTruNchEtabinStep2 || !state.pmeanRecoNchEtabinStep2 || !state.pmeanRecoEffcorrNchEtabinStep2 ||
        !state.pmeanMultTruNchEtabinStep2 || !state.pmeanMultRecoNchEtabinStep2 || !state.pmeanMultRecoEffcorrNchEtabinStep2) {
      LOGF(warning, "MC fluc: mean pT or mult map missing");
      return;
    }

    std::array<std::array<std::array<double, KIntK>, KIntM>, KNEtaMax> sumPmwkTru{};
    std::array<std::array<double, KIntK>, KNEtaMax> sumWkTru{};
    std::array<std::array<std::array<double, KIntK>, KIntM>, KNEtaMax> sumPmwkReco{};
    std::array<std::array<double, KIntK>, KNEtaMax> sumWkReco{};
    std::array<std::array<std::array<double, KIntK>, KIntM>, KNEtaMax> sumPmwkRecoEffCor{};
    std::array<std::array<double, KIntK>, KNEtaMax> sumWkRecoEffCor{};

    std::array<double, KNEtaMax> meanTru{}, c2Tru{};
    std::array<double, KNEtaMax> meanReco{}, c2Reco{};
    std::array<double, KNEtaMax> meanRecoEffCor{}, c2RecoEffCor{};

    std::array<double, KNEtaMax> meanTruMult{}, meanRecoMult{}, meanRecoEffCorMult{};
    std::array<double, KNEtaMax> p1kBarTru{}, p1kBarReco{}, p1kBarRecoEffCor{};
    std::array<double, KNEtaMax> p1kBarTruMult{}, p1kBarRecoMult{}, p1kBarRecoEffCorMult{};

    if (!mcCollision.has_mcCollision() || !isEventSelected(mcCollision)) {
      return;
    }
    float cent = getCentrality(mcCollision);
    if (cent > KCentMax) {
      return;
    }
    float multPV = mcCollision.multNTracksPV();
    float vz = mcCollision.posZ();
    if (!isPassAddPileup(multPV, mcTracks.size(), cent)) {
      return;
    }

    histos.fill(HIST("hVtxZ_after_sel"), mcCollision.posZ());
    histos.fill(HIST("hCentrality"), cent);
    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), multPV, mcTracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), mcTracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, mcTracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, multPV);

    double p1kBarFt0A = 0.0, p1kBarFt0C = 0.0;

    // --- truth sums ---
    for (const auto& particle : mcParticles) {
      if (particle.mcCollisionId() != mcCollision.mcCollisionId()) {
        continue;
      }
      if (!isParticleSelected(particle) || !particle.isPhysicalPrimary()) {
        continue;
      }
      float pt = particle.pt();
      if (pt <= cfgPtMin || pt > cfgPtMax) {
        continue;
      }
      float eta = particle.eta();
      for (int ieta = 0; ieta < nEta; ++ieta) {
        if (eta <= etaLw[ieta] || eta > etaUp[ieta]) {
          continue;
        }
        for (int k = 0; k < KIntK; ++k) {
          for (int m = 0; m < KIntM; ++m) {
            sumPmwkTru[ieta][m][k] += std::pow(pt, m);
          }
          sumWkTru[ieta][k]++;
        }
      }
    }

    // --- reco sums ---
    for (const auto& track : mcTracks) {
      if (track.collisionId() != mcCollision.index()) {
        continue;
      }
      if (!isTrackSelected(track)) {
        continue;
      }
      float pt = track.pt();
      if (pt <= cfgPtMin || pt > cfgPtMax) {
        continue;
      }
      float eta = track.eta(), phi = track.phi();
      auto sign = track.sign();
      histos.fill(HIST("hPt"), pt);
      histos.fill(HIST("hEta"), eta);
      histos.fill(HIST("hPhi"), phi);

      float eff = getEfficiency(multPV, pt, eta, 0, cfgEff);
      float fake = getEfficiency(multPV, pt, eta, 1, cfgEff);
      float flatW = getFlatteningWeight(vz, sign, pt, eta, phi, cfgFlat);
      float w = flatW * (1.0 - fake) / eff;
      if (!std::isfinite(w) || w <= 0.f || eff <= KFloatEpsilon) {
        continue;
      }

      for (int ieta = 0; ieta < nEta; ++ieta) {
        if (eta <= etaLw[ieta] || eta > etaUp[ieta]) {
          continue;
        }
        for (int k = 0; k < KIntK; ++k) {
          for (int m = 0; m < KIntM; ++m) {
            sumPmwkReco[ieta][m][k] += std::pow(1.0, k) * std::pow(pt, m);
            sumPmwkRecoEffCor[ieta][m][k] += std::pow(w, k) * std::pow(pt, m);
          }
          sumWkReco[ieta][k] += std::pow(1.0, k);
          sumWkRecoEffCor[ieta][k] += std::pow(w, k);
        }
      }

      histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
      histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
      histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0 - fake) / eff);
    }

    if (!hasMinTracksInAllEtaBins(sumWkTru) || !hasMinTracksInAllEtaBins(sumWkReco)) {
      return;
    }

    for (int ieta = 0; ieta < nEta; ++ieta) {
      const int ibx = state.pmeanTruNchEtabinStep2->GetXaxis()->FindBin(multPV);
      const int iby = ieta + 1;

      meanTruMult[ieta] = sumWkTru[ieta][1];
      meanRecoMult[ieta] = sumWkReco[ieta][1];
      meanRecoEffCorMult[ieta] = sumWkRecoEffCor[ieta][1];

      float mmptTru = state.pmeanTruNchEtabinStep2->GetBinContent(ibx, iby);
      float mmptReco = state.pmeanRecoNchEtabinStep2->GetBinContent(ibx, iby);
      float mmptRecoEffCor = state.pmeanRecoEffcorrNchEtabinStep2->GetBinContent(ibx, iby);

      float mmMultTru = state.pmeanMultTruNchEtabinStep2->GetBinContent(ibx, iby);
      float mmMultReco = state.pmeanMultRecoNchEtabinStep2->GetBinContent(ibx, iby);
      float mmMultRecoEffCor = state.pmeanMultRecoEffcorrNchEtabinStep2->GetBinContent(ibx, iby);

      if (std::isfinite(mmptTru)) {
        std::tie(meanTru[ieta], c2Tru[ieta]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkTru[ieta], sumWkTru[ieta], mmptTru);
      }
      if (std::isfinite(mmptReco)) {
        std::tie(meanReco[ieta], c2Reco[ieta]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkReco[ieta], sumWkReco[ieta], mmptReco);
      }
      if (std::isfinite(mmptRecoEffCor)) {
        std::tie(meanRecoEffCor[ieta], c2RecoEffCor[ieta]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumPmwkRecoEffCor[ieta], sumWkRecoEffCor[ieta], mmptRecoEffCor);
      }

      if (mmptTru != 0.0f) {
        p1kBarTru[ieta] = meanTru[ieta] - mmptTru;
      }
      if (mmptReco != 0.0f) {
        p1kBarReco[ieta] = meanReco[ieta] - mmptReco;
      }
      if (mmptRecoEffCor != 0.0f) {
        p1kBarRecoEffCor[ieta] = meanRecoEffCor[ieta] - mmptRecoEffCor;
      }

      if (mmMultTru != 0.0f) {
        p1kBarTruMult[ieta] = meanTruMult[ieta] - mmMultTru;
      }
      if (mmMultReco != 0.0f) {
        p1kBarRecoMult[ieta] = meanRecoMult[ieta] - mmMultReco;
      }
      if (mmMultRecoEffCor != 0.0f) {
        p1kBarRecoEffCorMult[ieta] = meanRecoEffCorMult[ieta] - mmMultRecoEffCor;
      }
    }

    double amplFT0A = 0, amplFT0C = 0;
    if (mcCollision.has_foundFT0()) {
      const auto& ft0 = mcCollision.foundFT0();
      for (std::size_t iCh = 0; iCh < ft0.channelA().size(); iCh++) {
        amplFT0A += ft0.amplitudeA()[iCh];
      }
      for (std::size_t iCh = 0; iCh < ft0.channelC().size(); iCh++) {
        amplFT0C += ft0.amplitudeC()[iCh];
      }
    }
    p1kBarFt0A = amplFT0A - state.pmeanFT0AmultpvStep2->GetBinContent(state.pmeanFT0AmultpvStep2->GetXaxis()->FindBin(multPV));
    p1kBarFt0C = amplFT0C - state.pmeanFT0CmultpvStep2->GetBinContent(state.pmeanFT0CmultpvStep2->GetXaxis()->FindBin(multPV));

    // per-eta counts & means
    for (int ieta = 0; ieta < nEta; ++ieta) {
      histos.fill(HIST("MCGen/Prof_Cent_NEta_Nchrec"), cent, ieta, sumWkTru[ieta][1]);
      histos.fill(HIST("MCGen/Prof_Mult_NEta_Nchrec"), multPV, ieta, sumWkTru[ieta][1]);
      histos.fill(HIST("MCReco/Prof_Cent_NEta_Nchrec"), cent, ieta, sumWkReco[ieta][1]);
      histos.fill(HIST("MCReco/Prof_Mult_NEta_Nchrec"), multPV, ieta, sumWkReco[ieta][1]);
      histos.fill(HIST("MCRecoEffCorr/Prof_Cent_NEta_Nchrec"), cent, ieta, sumWkRecoEffCor[ieta][1]);
      histos.fill(HIST("MCRecoEffCorr/Prof_Mult_NEta_Nchrec"), multPV, ieta, sumWkRecoEffCor[ieta][1]);

      if (sumWkTru[ieta][1] > 1.0f) {
        histos.fill(HIST("MCGen/Prof_Cent_NEta_MeanpT"), cent, ieta, meanTru[ieta]);
        histos.fill(HIST("MCGen/Prof_Mult_NEta_MeanpT"), multPV, ieta, meanTru[ieta]);
      }
      if (sumWkReco[ieta][1] > 1.0f) {
        histos.fill(HIST("MCReco/Prof_Cent_NEta_MeanpT"), cent, ieta, meanReco[ieta]);
        histos.fill(HIST("MCReco/Prof_Mult_NEta_MeanpT"), multPV, ieta, meanReco[ieta]);
      }
      if (sumWkRecoEffCor[ieta][1] > 1.0f) {
        histos.fill(HIST("MCRecoEffCorr/Prof_Cent_NEta_MeanpT"), cent, ieta, meanRecoEffCor[ieta]);
        histos.fill(HIST("MCRecoEffCorr/Prof_Mult_NEta_MeanpT"), multPV, ieta, meanRecoEffCor[ieta]);
      }
    }

    // meanpT & C2 vs eta bin
    for (int ieta = 0; ieta < nEta; ++ieta) {
      if (std::isfinite(meanTru[ieta])) {
        histos.fill(HIST("MCGen/Prof_MeanpT_Cent_etabin"), cent, ieta, meanTru[ieta]);
        histos.fill(HIST("MCGen/Prof_MeanpT_Mult_etabin"), multPV, ieta, meanTru[ieta]);
      }
      if (std::isfinite(c2Tru[ieta])) {
        histos.fill(HIST("MCGen/Prof_C2_Cent_etabin"), cent, ieta, c2Tru[ieta]);
        histos.fill(HIST("MCGen/Prof_C2_Mult_etabin"), multPV, ieta, c2Tru[ieta]);
      }
      if (std::isfinite(meanReco[ieta])) {
        histos.fill(HIST("MCReco/Prof_MeanpT_Cent_etabin"), cent, ieta, meanReco[ieta]);
        histos.fill(HIST("MCReco/Prof_MeanpT_Mult_etabin"), multPV, ieta, meanReco[ieta]);
      }
      if (std::isfinite(c2Reco[ieta])) {
        histos.fill(HIST("MCReco/Prof_C2_Cent_etabin"), cent, ieta, c2Reco[ieta]);
        histos.fill(HIST("MCReco/Prof_C2_Mult_etabin"), multPV, ieta, c2Reco[ieta]);
      }
      if (std::isfinite(meanRecoEffCor[ieta])) {
        histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Cent_etabin"), cent, ieta, meanRecoEffCor[ieta]);
        histos.fill(HIST("MCRecoEffCorr/Prof_MeanpT_Mult_etabin"), multPV, ieta, meanRecoEffCor[ieta]);
      }
      if (std::isfinite(c2RecoEffCor[ieta])) {
        histos.fill(HIST("MCRecoEffCorr/Prof_C2_Cent_etabin"), cent, ieta, c2RecoEffCor[ieta]);
        histos.fill(HIST("MCRecoEffCorr/Prof_C2_Mult_etabin"), multPV, ieta, c2RecoEffCor[ieta]);
      }
    }

    // mirror-pair subevent (C2Sub) & covariances vs eta bin
    for (int ietaA = 1; ietaA <= (nEta - 1) / 2; ++ietaA) {
      int ietaC = nEta - ietaA;

      float c2SubTru = p1kBarTru[ietaA] * p1kBarTru[ietaC];
      float c2SubReco = p1kBarReco[ietaA] * p1kBarReco[ietaC];
      float c2SubRecoEffCor = p1kBarRecoEffCor[ietaA] * p1kBarRecoEffCor[ietaC];

      float covTru = p1kBarTruMult[ietaA] * p1kBarTru[ietaC];
      float covReco = p1kBarRecoMult[ietaA] * p1kBarReco[ietaC];
      float covRecoEffCor = p1kBarRecoEffCorMult[ietaA] * p1kBarRecoEffCor[ietaC];

      if (std::isfinite(c2SubTru)) {
        histos.fill(HIST("MCGen/Prof_C2Sub_Cent_etabin"), cent, ietaA, c2SubTru);
        histos.fill(HIST("MCGen/Prof_C2Sub_Mult_etabin"), multPV, ietaA, c2SubTru);
      }
      if (std::isfinite(c2SubReco)) {
        histos.fill(HIST("MCReco/Prof_C2Sub_Cent_etabin"), cent, ietaA, c2SubReco);
        histos.fill(HIST("MCReco/Prof_C2Sub_Mult_etabin"), multPV, ietaA, c2SubReco);
      }
      if (std::isfinite(c2SubRecoEffCor)) {
        histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Cent_etabin"), cent, ietaA, c2SubRecoEffCor);
        histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub_Mult_etabin"), multPV, ietaA, c2SubRecoEffCor);
      }
      if (std::isfinite(covTru)) {
        histos.fill(HIST("MCGen/Prof_Cov_Cent_etabin"), cent, ietaA, covTru);
        histos.fill(HIST("MCGen/Prof_Cov_Mult_etabin"), multPV, ietaA, covTru);
      }
      if (std::isfinite(covReco)) {
        histos.fill(HIST("MCReco/Prof_Cov_Cent_etabin"), cent, ietaA, covReco);
        histos.fill(HIST("MCReco/Prof_Cov_Mult_etabin"), multPV, ietaA, covReco);
      }
      if (std::isfinite(covRecoEffCor)) {
        histos.fill(HIST("MCRecoEffCorr/Prof_Cov_Cent_etabin"), cent, ietaA, covRecoEffCor);
        histos.fill(HIST("MCRecoEffCorr/Prof_Cov_Mult_etabin"), multPV, ietaA, covRecoEffCor);
      }
    }

    // FT0 covariance vs narrow eta bin (full range, indexed by the actual pT bin)
    for (int ieta = 1; ieta < nEta; ++ieta) {
      float covFT0ATru = p1kBarFt0A * p1kBarTru[ieta];
      float covFT0AReco = p1kBarFt0A * p1kBarReco[ieta];
      float covFT0ARecoEffCor = p1kBarFt0A * p1kBarRecoEffCor[ieta];
      float covFT0CTru = p1kBarFt0C * p1kBarTru[ieta];
      float covFT0CReco = p1kBarFt0C * p1kBarReco[ieta];
      float covFT0CRecoEffCor = p1kBarFt0C * p1kBarRecoEffCor[ieta];

      if (std::isfinite(covFT0ATru)) {
        histos.fill(HIST("MCGen/Prof_CovFT0A_Cent_etabin"), cent, ieta, covFT0ATru);
        histos.fill(HIST("MCGen/Prof_CovFT0A_Mult_etabin"), multPV, ieta, covFT0ATru);
      }
      if (std::isfinite(covFT0AReco)) {
        histos.fill(HIST("MCReco/Prof_CovFT0A_Cent_etabin"), cent, ieta, covFT0AReco);
        histos.fill(HIST("MCReco/Prof_CovFT0A_Mult_etabin"), multPV, ieta, covFT0AReco);
      }
      if (std::isfinite(covFT0ARecoEffCor)) {
        histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A_Cent_etabin"), cent, ieta, covFT0ARecoEffCor);
        histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A_Mult_etabin"), multPV, ieta, covFT0ARecoEffCor);
      }
      if (std::isfinite(covFT0CTru)) {
        histos.fill(HIST("MCGen/Prof_CovFT0C_Cent_etabin"), cent, ieta, covFT0CTru);
        histos.fill(HIST("MCGen/Prof_CovFT0C_Mult_etabin"), multPV, ieta, covFT0CTru);
      }
      if (std::isfinite(covFT0CReco)) {
        histos.fill(HIST("MCReco/Prof_CovFT0C_Cent_etabin"), cent, ieta, covFT0CReco);
        histos.fill(HIST("MCReco/Prof_CovFT0C_Mult_etabin"), multPV, ieta, covFT0CReco);
      }
      if (std::isfinite(covFT0CRecoEffCor)) {
        histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C_Cent_etabin"), cent, ieta, covFT0CRecoEffCor);
        histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C_Mult_etabin"), multPV, ieta, covFT0CRecoEffCor);
      }
    }

    // full 2D subevent map
    for (int ietaA = 1; ietaA < nEta; ++ietaA) {
      for (int ietaC = 1; ietaC < nEta; ++ietaC) {
        float etaValA = (etaLw[ietaA] + etaUp[ietaA]) / 2.0f;
        float etaValB = (etaLw[ietaC] + etaUp[ietaC]) / 2.0f;
        float gap = etaValA - etaValB;
        float sum = (etaValA + etaValB);

        float c2SubTru = (ietaA == ietaC) ? static_cast<float>(c2Tru[ietaA]) : p1kBarTru[ietaA] * p1kBarTru[ietaC];
        float c2SubReco = (ietaA == ietaC) ? static_cast<float>(c2Reco[ietaA]) : p1kBarReco[ietaA] * p1kBarReco[ietaC];
        float c2SubRecoEffCor = (ietaA == ietaC) ? static_cast<float>(c2RecoEffCor[ietaA]) : p1kBarRecoEffCor[ietaA] * p1kBarRecoEffCor[ietaC];

        float covTru = p1kBarTruMult[ietaA] * p1kBarTru[ietaC];
        float covReco = p1kBarRecoMult[ietaA] * p1kBarReco[ietaC];
        float covRecoEffCor = p1kBarRecoEffCorMult[ietaA] * p1kBarRecoEffCor[ietaC];

        float covFT0ATru = p1kBarFt0A * p1kBarTru[ietaC];
        float covFT0AReco = p1kBarFt0A * p1kBarReco[ietaC];
        float covFT0ARecoEffCor = p1kBarFt0A * p1kBarRecoEffCor[ietaC];

        float covFT0CTru = p1kBarFt0C * p1kBarTru[ietaA];
        float covFT0CReco = p1kBarFt0C * p1kBarReco[ietaA];
        float covFT0CRecoEffCor = p1kBarFt0C * p1kBarRecoEffCor[ietaA];

        if (std::isfinite(c2SubTru)) {
          histos.fill(HIST("MCGen/Prof_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2SubTru);
          histos.fill(HIST("MCGen/Prof_GapSum2D"), cent, gap, sum, c2SubTru);
        }
        if (std::isfinite(c2SubReco)) {
          histos.fill(HIST("MCReco/Prof_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2SubReco);
          histos.fill(HIST("MCReco/Prof_GapSum2D"), cent, gap, sum, c2SubReco);
        }
        if (std::isfinite(c2SubRecoEffCor)) {
          histos.fill(HIST("MCRecoEffCorr/Prof_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2SubRecoEffCor);
          histos.fill(HIST("MCRecoEffCorr/Prof_GapSum2D"), cent, gap, sum, c2SubRecoEffCor);
        }

        if (std::isfinite(covTru)) {
          histos.fill(HIST("MCGen/Prof_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covTru);
        }
        if (std::isfinite(covReco)) {
          histos.fill(HIST("MCReco/Prof_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covReco);
        }
        if (std::isfinite(covRecoEffCor)) {
          histos.fill(HIST("MCRecoEffCorr/Prof_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covRecoEffCor);
        }

        if (std::isfinite(covFT0ATru)) {
          histos.fill(HIST("MCGen/Prof_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0ATru);
        }
        if (std::isfinite(covFT0AReco)) {
          histos.fill(HIST("MCReco/Prof_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0AReco);
        }
        if (std::isfinite(covFT0ARecoEffCor)) {
          histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0ARecoEffCor);
        }

        if (std::isfinite(covFT0CTru)) {
          histos.fill(HIST("MCGen/Prof_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0CTru);
        }
        if (std::isfinite(covFT0CReco)) {
          histos.fill(HIST("MCReco/Prof_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0CReco);
        }
        if (std::isfinite(covFT0CRecoEffCor)) {
          histos.fill(HIST("MCRecoEffCorr/Prof_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0CRecoEffCor);
        }
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processMCFluc, "process MC to calculate pt fluc", cfgRunMCFluc);

  // ===========================================================================
  // DATA: build flattening maps (inclusive)
  // ===========================================================================
  void processGetDataFlat(AodCollisionsSel::iterator const& coll, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcsData*/, AodTracksSel const& tracks)
  {
    histos.fill(HIST("hVtxZ"), coll.posZ());
    if (!isEventSelected(coll)) {
      return;
    }
    float cent = getCentrality(coll);
    if (cent > KCentMax) {
      return;
    }
    if (!isPassAddPileup(coll.multNTracksPV(), tracks.size(), cent)) {
      return;
    }

    histos.fill(HIST("hVtxZ_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), cent);
    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), tracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), tracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, tracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, coll.multNTracksPV());

    int ntrk = 0;
    float vz = coll.posZ();

    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      float pt = track.pt();
      if (pt <= cfgPtMin || pt > cfgPtMax) {
        continue;
      }
      float eta = track.eta(), phi = track.phi();
      auto sign = track.sign();

      histos.fill(HIST("hPt"), pt);
      histos.fill(HIST("hEta"), eta);
      histos.fill(HIST("hPhi"), phi);
      if (eta > etaLw[0] && eta < etaUp[0]) {
        ntrk++;
      }

      float eff = getEfficiency(coll.multNTracksPV(), pt, eta, 0, cfgEff);
      if (eff <= KFloatEpsilon) {
        continue;
      }
      float fake = getEfficiency(coll.multNTracksPV(), pt, eta, 1, cfgEff);
      float w = (1.0f - fake) / eff;
      if (!std::isfinite(w) || w <= 0.f) {
        continue;
      }

      histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
      histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
      histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);
    }

    if (cfgZDC) {
      const auto& foundBC = coll.foundBC_as<BCsRun3>();
      if (!foundBC.has_zdc()) {
        return;
      }
      auto zdc = foundBC.zdc();
      auto zdcAmp = zdc.energyCommonZNA() + zdc.energyCommonZNC();
      histos.fill(HIST("hnTrkPVZDC"), coll.multNTracksPV(), zdcAmp);
      histos.fill(HIST("hNchZDC"), ntrk, zdcAmp);
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processGetDataFlat, "process data to calculate Flattening maps", cfgRunGetDataFlat);

  // ===========================================================================
  // DATA: mean pT
  // ===========================================================================
  void processDataMean(AodCollisionsSel::iterator const& coll, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcsData*/, aod::FT0s const&, AodTracksSel const& tracks)
  {
    std::array<double, KNEtaMax> sumWi{}, sumWipti{};

    if (!isEventSelected(coll)) {
      return;
    }
    float cent = getCentrality(coll);
    if (cent > KCentMax) {
      return;
    }
    if (!isPassAddPileup(coll.multNTracksPV(), tracks.size(), cent)) {
      return;
    }

    histos.fill(HIST("hVtxZ_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), cent);
    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), tracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), tracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, tracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, coll.multNTracksPV());

    float vz = coll.posZ();

    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      float p = track.p();
      float pt = track.pt();
      float eta = track.eta(), phi = track.phi();
      auto sign = track.sign();
      if (p < KFloatEpsilon) {
        continue;
      }
      if (pt <= cfgPtMin || pt > cfgPtMax) {
        continue;
      }

      histos.fill(HIST("hP"), p);
      histos.fill(HIST("hPt"), pt);
      histos.fill(HIST("hEta"), eta);
      histos.fill(HIST("hPhi"), phi);

      float eff = getEfficiency(coll.multNTracksPV(), pt, eta, 0, cfgEff);
      float fake = getEfficiency(coll.multNTracksPV(), pt, eta, 1, cfgEff);
      float flatWeight = getFlatteningWeight(vz, sign, pt, eta, phi, cfgFlat);

      histos.fill(HIST("pEffWeight_pt_eta_cent"), pt, eta, cent, eff);
      histos.fill(HIST("pFakeWeight_pt_eta_cent"), pt, eta, cent, fake);
      histos.fill(HIST("pFlatWeight_pt_eta_cent"), pt, eta, cent, flatWeight);

      if (eff <= KFloatEpsilon) {
        continue;
      }

      float w = flatWeight * (1.0f - fake) / eff;
      if (!std::isfinite(w) || w <= 0.f) {
        continue;
      }

      histos.fill(HIST("hEtaPhiReco"), vz, sign, pt, eta, phi);
      histos.fill(HIST("hEtaPhiRecoEffWtd"), vz, sign, pt, eta, phi, (1.0f - fake) / eff);
      histos.fill(HIST("hEtaPhiRecoWtd"), vz, sign, pt, eta, phi, w);

      for (int ieta = 0; ieta < nEta; ++ieta) {
        if (eta <= etaLw[ieta] || eta > etaUp[ieta]) {
          continue;
        }
        sumWi[ieta] += w;
        sumWipti[ieta] += w * pt;
      }
    }

    if (!hasMinTracksInAllEtaBins(sumWi)) {
      return;
    }

    if (sumWi[0] >= 1.0f) {
      histos.fill(HIST("Prof_Cent_Nchrec"), cent, sumWi[0]);
      histos.fill(HIST("Prof_Mult_Nchrec"), coll.multNTracksPV(), sumWi[0]);
      histos.fill(HIST("Prof_Cent_MeanpT"), cent, sumWipti[0] / sumWi[0]);
      histos.fill(HIST("Prof_Mult_MeanpT"), coll.multNTracksPV(), sumWipti[0] / sumWi[0]);
    }

    for (int ietaA = 0; ietaA < nEta; ++ietaA) {
      for (int ietaC = 0; ietaC < nEta; ++ietaC) {
        if ((sumWi[ietaA] < 1.0f) || (sumWi[ietaC] < 1.0f)) {
          continue;
        }
        double wCorrAB = sumWi[ietaA] + sumWi[ietaC];
        if (wCorrAB > 0) {
          float mptsub = (sumWipti[ietaA] + sumWipti[ietaC]) / wCorrAB;
          histos.fill(HIST("Prof2D_MeanpTSub"), cent, ietaA, ietaC, mptsub);
        }
      }
      if (sumWi[ietaA] >= 1.0f) {
        double mpt = sumWipti[ietaA] / sumWi[ietaA];
        if (std::isfinite(mpt)) {
          histos.fill(HIST("pmean_nch_etabin"), coll.multNTracksPV(), ietaA, mpt);
          histos.fill(HIST("pmeanMult_nch_etabin"), coll.multNTracksPV(), ietaA, sumWi[ietaA]);
          histos.fill(HIST("pmean_cent_etabin"), cent, ietaA, mpt);
          histos.fill(HIST("pmeanMult_cent_etabin"), cent, ietaA, sumWi[ietaA]);
        }
      }
    }

    double amplFT0A = 0, amplFT0C = 0;
    if (coll.has_foundFT0()) {
      const auto& ft0 = coll.foundFT0();
      for (std::size_t iCh = 0; iCh < ft0.channelA().size(); iCh++) {
        auto chanelid = ft0.channelA()[iCh];
        float ampl = ft0.amplitudeA()[iCh];
        amplFT0A += ampl;
        auto eta = getEtaFT0(chanelid, 0);
        histos.fill(HIST("pmean_cent_id_eta_FT0"), cent, chanelid, eta, ampl);
        histos.fill(HIST("h3_cent_id_eta_FT0"), cent, chanelid, eta, ampl);
      }
      for (std::size_t iCh = 0; iCh < ft0.channelC().size(); iCh++) {
        auto chanelid = ft0.channelC()[iCh];
        auto globalId = chanelid + KnFt0cCell;
        float ampl = ft0.amplitudeC()[iCh];
        amplFT0C += ampl;
        auto eta = getEtaFT0(globalId, 1);
        histos.fill(HIST("pmean_cent_id_eta_FT0"), cent, globalId, eta, ampl);
        histos.fill(HIST("h3_cent_id_eta_FT0"), cent, globalId, eta, ampl);
      }
    }
    histos.fill(HIST("pmeanFT0Amultpv"), coll.multNTracksPV(), amplFT0A);
    histos.fill(HIST("pmeanFT0A_cent"), cent, amplFT0A);
    histos.fill(HIST("pmeanFT0Cmultpv"), coll.multNTracksPV(), amplFT0C);
    histos.fill(HIST("pmeanFT0C_cent"), cent, amplFT0C);
  }
  PROCESS_SWITCH(RadialFlowDecorr, processDataMean, "process data to calculate mean pT", cfgRunDataMean);

  // ===========================================================================
  // DATA: fluctuations (C2, subevent) + Poisson bootstrap for the base run
  // ===========================================================================
  void processDataFluc(AodCollisionsSel::iterator const& coll, BCsRun3 const& /*bcs*/, aod::Zdcs const& /*zdcsData*/, aod::FT0s const&, AodTracksSel const& tracks)
  {
    if (!isEventSelected(coll)) {
      return;
    }
    float cent = getCentrality(coll);
    if (cent > KCentMax) {
      return;
    }
    if (!isPassAddPileup(coll.multNTracksPV(), tracks.size(), cent)) {
      return;
    }

    histos.fill(HIST("hVtxZ_after_sel"), coll.posZ());
    histos.fill(HIST("hCentrality"), cent);
    histos.fill(HIST("Hist2D_globalTracks_PVTracks"), coll.multNTracksPV(), tracks.size());
    histos.fill(HIST("Hist2D_cent_nch"), tracks.size(), cent);
    histos.fill(HIST("Hist2D_globalTracks_cent"), cent, tracks.size());
    histos.fill(HIST("Hist2D_PVTracks_cent"), cent, coll.multNTracksPV());

    if (!state.pmeanNchEtabinStep2 || !state.pmeanMultNchEtabinStep2) {
      LOGF(warning, "Data fluc: mean pT or mult map missing");
      return;
    }
    if ((cfgEff || cfgFlat) && (!state.hEff || !state.hFake || !state.hFlatWeight)) {
      LOGF(warning, "Data fluc: correction maps requested but not all present.");
      return;
    }

    std::array<std::array<std::array<double, KIntK>, KIntM>, KNEtaMax> sumpmwk{};
    std::array<std::array<double, KIntK>, KNEtaMax> sumwk{};
    std::array<double, KNEtaMax> mean{}, c2{}, p1kBar{};
    std::array<double, KNEtaMax> meanMult{}, p1kBarMult{};

    float vz = coll.posZ();

    for (const auto& track : tracks) {
      if (!isTrackSelected(track)) {
        continue;
      }
      float p = track.p();
      float pt = track.pt();
      float eta = track.eta(), phi = track.phi();
      auto sign = track.sign();
      if (p < KFloatEpsilon) {
        continue;
      }
      if (pt <= cfgPtMin || pt > cfgPtMax) {
        continue;
      }

      float eff = getEfficiency(coll.multNTracksPV(), pt, eta, 0, cfgEff);
      if (eff <= KFloatEpsilon) {
        continue;
      }
      float fake = getEfficiency(coll.multNTracksPV(), pt, eta, 1, cfgEff);
      float flatWeight = getFlatteningWeight(vz, sign, pt, eta, phi, cfgFlat);
      float w = flatWeight * (1.0f - fake) / eff;
      if (!std::isfinite(w) || w <= 0.f) {
        continue;
      }

      for (int ieta = 0; ieta < nEta; ++ieta) {
        if (eta <= etaLw[ieta] || eta > etaUp[ieta]) {
          continue;
        }
        for (int k = 0; k < KIntK; ++k) {
          for (int m = 0; m < KIntM; ++m) {
            sumpmwk[ieta][m][k] += std::pow(w, k) * std::pow(pt, m);
          }
          sumwk[ieta][k] += std::pow(w, k);
        }
      }
    }

    if (!hasMinTracksInAllEtaBins(sumwk)) {
      return;
    }

    double amplFT0A = 0, amplFT0C = 0;
    if (coll.has_foundFT0()) {
      const auto& ft0 = coll.foundFT0();
      for (std::size_t iCh = 0; iCh < ft0.channelA().size(); iCh++) {
        amplFT0A += ft0.amplitudeA()[iCh];
      }
      for (std::size_t iCh = 0; iCh < ft0.channelC().size(); iCh++) {
        amplFT0C += ft0.amplitudeC()[iCh];
      }
    }
    double p1kBarFt0A = amplFT0A - state.pmeanFT0AmultpvStep2->GetBinContent(state.pmeanFT0AmultpvStep2->GetXaxis()->FindBin(coll.multNTracksPV()));
    double p1kBarFt0C = amplFT0C - state.pmeanFT0CmultpvStep2->GetBinContent(state.pmeanFT0CmultpvStep2->GetXaxis()->FindBin(coll.multNTracksPV()));

    for (int ieta = 0; ieta < nEta; ++ieta) {
      const int ibx = state.pmeanNchEtabinStep2->GetXaxis()->FindBin(coll.multNTracksPV());
      const int iby = ieta + 1;

      float mmpt = state.pmeanNchEtabinStep2->GetBinContent(ibx, iby);
      float mmMult = state.pmeanMultNchEtabinStep2->GetBinContent(ibx, iby);

      mean[ieta] = sumpmwk[ieta][1][1] / sumwk[ieta][1];
      meanMult[ieta] = sumwk[ieta][1];

      if (std::isfinite(mmpt) && mmpt != 0) {
        std::tie(mean[ieta], c2[ieta]) = calculateMeanAndC2FromSums<KIntM, KIntK>(sumpmwk[ieta], sumwk[ieta], mmpt);
        p1kBar[ieta] = mean[ieta] - mmpt;
      }
      p1kBarMult[ieta] = meanMult[ieta] - mmMult;
    }

    // --- Poisson bootstrap: one weight per sample per event ---
    std::array<double, KMaxBoot> poisW{};
    if (doBoot) {
      for (int s = 0; s < nBoot; ++s) {
        poisW[s] = rng.Poisson(1.0);
      }
    }
    auto fillBS2D = [&](std::array<std::shared_ptr<TProfile2D>, KMaxBoot>& arr, double x, double y, double val) {
      if (!doBoot) {
        return;
      }
      for (int s = 0; s < nBoot; ++s) {
        arr[s]->Fill(x, y, val, poisW[s]);
      }
    };
    auto fillBS3D = [&](std::array<std::shared_ptr<TProfile3D>, KMaxBoot>& arr, double x, double y, double z, double val) {
      if (!doBoot) {
        return;
      }
      for (int s = 0; s < nBoot; ++s) {
        arr[s]->Fill(x, y, z, val, poisW[s]);
      }
    };

    // meanpT & C2 vs eta bin
    for (int ieta = 0; ieta < nEta; ++ieta) {
      if (std::isfinite(mean[ieta])) {
        histos.fill(HIST("Prof_MeanpT_Cent_etabin"), cent, ieta, mean[ieta]);
        histos.fill(HIST("Prof_MeanpT_Mult_etabin"), coll.multNTracksPV(), ieta, mean[ieta]);
        fillBS2D(bs.meanpTCent, cent, ieta, mean[ieta]);
        fillBS2D(bs.meanpTMult, coll.multNTracksPV(), ieta, mean[ieta]);
      }
      if (std::isfinite(c2[ieta])) {
        histos.fill(HIST("Prof_C2_Cent_etabin"), cent, ieta, c2[ieta]);
        histos.fill(HIST("Prof_C2_Mult_etabin"), coll.multNTracksPV(), ieta, c2[ieta]);
        fillBS2D(bs.c2Cent, cent, ieta, c2[ieta]);
        fillBS2D(bs.c2Mult, coll.multNTracksPV(), ieta, c2[ieta]);
      }
    }

    // mirror-pair subevent (C2Sub) & covariances vs eta bin
    for (int ietaA = 1; ietaA <= (nEta - 1) / 2; ++ietaA) {
      int ietaC = nEta - ietaA;
      float c2Sub = p1kBar[ietaA] * p1kBar[ietaC];
      float covAC = p1kBarMult[ietaA] * p1kBar[ietaC];
      float covCA = p1kBar[ietaA] * p1kBarMult[ietaC];

      if (std::isfinite(c2Sub)) {
        histos.fill(HIST("Prof_C2Sub_Cent_etabin"), cent, ietaA, c2Sub);
        histos.fill(HIST("Prof_C2Sub_Mult_etabin"), coll.multNTracksPV(), ietaA, c2Sub);
        fillBS2D(bs.c2SubCent, cent, ietaA, c2Sub);
        fillBS2D(bs.c2SubMult, coll.multNTracksPV(), ietaA, c2Sub);
      }
      if (std::isfinite(covAC)) {
        histos.fill(HIST("Prof_Cov_Cent_etabin"), cent, ietaA, covAC);
        histos.fill(HIST("Prof_Cov_Mult_etabin"), coll.multNTracksPV(), ietaA, covAC);
        fillBS2D(bs.covCent, cent, ietaA, covAC);
        fillBS2D(bs.covMult, coll.multNTracksPV(), ietaA, covAC);
      }
      if (std::isfinite(covCA)) {
        histos.fill(HIST("Prof_Cov_Cent_etabin"), cent, ietaC, covCA);
        histos.fill(HIST("Prof_Cov_Mult_etabin"), coll.multNTracksPV(), ietaC, covCA);
        fillBS2D(bs.covCent, cent, ietaC, covCA);
        fillBS2D(bs.covMult, coll.multNTracksPV(), ietaC, covCA);
      }
    }

    // FT0 covariance vs narrow eta bin (full range, indexed by the actual pT bin)
    for (int ieta = 1; ieta < nEta; ++ieta) {
      float covFT0Aeta = p1kBarFt0A * p1kBar[ieta];
      float covFT0Ceta = p1kBarFt0C * p1kBar[ieta];
      if (std::isfinite(covFT0Aeta)) {
        histos.fill(HIST("Prof_CovFT0A_Cent_etabin"), cent, ieta, covFT0Aeta);
        histos.fill(HIST("Prof_CovFT0A_Mult_etabin"), coll.multNTracksPV(), ieta, covFT0Aeta);
        fillBS2D(bs.covFT0ACent, cent, ieta, covFT0Aeta);
        fillBS2D(bs.covFT0AMult, coll.multNTracksPV(), ieta, covFT0Aeta);
      }
      if (std::isfinite(covFT0Ceta)) {
        histos.fill(HIST("Prof_CovFT0C_Cent_etabin"), cent, ieta, covFT0Ceta);
        histos.fill(HIST("Prof_CovFT0C_Mult_etabin"), coll.multNTracksPV(), ieta, covFT0Ceta);
        fillBS2D(bs.covFT0CCent, cent, ieta, covFT0Ceta);
        fillBS2D(bs.covFT0CMult, coll.multNTracksPV(), ieta, covFT0Ceta);
      }
    }

    // full 2D subevent map
    for (int ietaA = 1; ietaA < nEta; ++ietaA) {
      for (int ietaC = 1; ietaC < nEta; ++ietaC) {
        float etaValA = (etaLw[ietaA] + etaUp[ietaA]) / 2.0f;
        float etaValB = (etaLw[ietaC] + etaUp[ietaC]) / 2.0f;
        float gap = etaValA - etaValB;
        float sum = (etaValA + etaValB);

        float c2Sub = (ietaA == ietaC) ? static_cast<float>(c2[ietaA]) : p1kBar[ietaA] * p1kBar[ietaC];
        float cov = p1kBarMult[ietaA] * p1kBar[ietaC];
        float covFT0A = p1kBarFt0A * p1kBar[ietaC];
        float covFT0C = p1kBarFt0C * p1kBar[ietaA];

        if (std::isfinite(c2Sub)) {
          histos.fill(HIST("Prof_C2Sub2D_Cent_etaA_etaC"), cent, etaValA, etaValB, c2Sub);
          histos.fill(HIST("Prof_GapSum2D"), cent, gap, sum, c2Sub);
          fillBS3D(bs.c2Sub2D, cent, etaValA, etaValB, c2Sub);
          fillBS3D(bs.gapSum2D, cent, gap, sum, c2Sub);
        }
        if (std::isfinite(cov)) {
          histos.fill(HIST("Prof_Cov2D_Cent_etaA_etaC"), cent, etaValA, etaValB, cov);
          fillBS3D(bs.cov2D, cent, etaValA, etaValB, cov);
        }
        if (std::isfinite(covFT0A)) {
          histos.fill(HIST("Prof_CovFT0A2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0A);
          fillBS3D(bs.covFT0A2D, cent, etaValA, etaValB, covFT0A);
        }
        if (std::isfinite(covFT0C)) {
          histos.fill(HIST("Prof_CovFT0C2D_Cent_etaA_etaC"), cent, etaValA, etaValB, covFT0C);
          fillBS3D(bs.covFT0C2D, cent, etaValA, etaValB, covFT0C);
        }
      }
    }
  }
  PROCESS_SWITCH(RadialFlowDecorr, processDataFluc, "process data to calculate fluc pT", cfgRunDataFluc);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<RadialFlowDecorr>(cfgc)};
  return workflow;
}
