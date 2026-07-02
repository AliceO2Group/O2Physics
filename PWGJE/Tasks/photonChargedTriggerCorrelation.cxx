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

/// \file photonChargedTriggerCorrelation.cxx
/// \author Julius Kinner
/// \brief photon-jet angular correlation analysis
///
/// Analysis for angular correlations between jets and photons via two-particle correlations with charged high-pt triggers
/// Associated hadrons (tracks), pipm, photons (PCM), pi0 (PCM)
/// Also contains checks and monte-carlo (efficiency, purity, mc-true correlation,...)
/// End goal of studying correlations between direct photons and jets

#include "PWGJE/DataModel/PhotonChargedTriggerCorrelation.h"

#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
//
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TableHelper.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TPDGCode.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <format>
#include <functional>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

using CorrCollisions = soa::Join<aod::JetCollisions, aod::CollisionsExtraCorr>;
using CorrCollision = CorrCollisions::iterator;
using CorrMcDCollisions = soa::Join<aod::JetCollisionsMCD, aod::CollisionsExtraCorr>;
using CorrMcDCollision = CorrMcDCollisions::iterator;
using CorrMcCollisions = soa::Join<aod::JetMcCollisions, aod::McCollisionsExtraCorr>;
using CorrMcCollision = CorrMcCollisions::iterator;

using BinningZPvMult = ColumnBinningPolicy<aod::jcollision::PosZ, aod::collision_extra_corr::NGlobalTracks>;

// correlation analysis ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// =/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/

const double absEtaMaxDefault = 0.8;

struct PhotonChargedTriggerCorrelation {
  // configurables

  // general (kenobi)
  Configurable<std::string> pathCcdbCorrection{"pathCcdbCorrection", "Users/j/jkinner/correction/set_in_config", "base path to the ccdb corrections"};
  Configurable<std::string> urlCcdb{"urlCcdb", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> noLaterThanCcdbConfig{"noLaterThanCcdbConfig", -1, "latest acceptable timestamp of creation for the object (-1 for task start time)"};
  Configurable<int> splitMcEvN{"splitMcEvN", 1, "number of equal event fractions for selection in MC (max 10) | (select 1 for data)"};
  Configurable<int> splitMcEvSelect{"splitMcEvSelect", 1, "selected event fraction number in MC (starts at 1) | (select 1 for data)"};

  // analysis
  Configurable<bool> applyCTotTrigger{"applyCTotTrigger", false, "whether to apply on-the-fly total correction factor for triggers"};
  Configurable<bool> doTrigEvMixing{"doTrigEvMixing", false, "whether to use trigger events for trigger mixing"};

  Configurable<int> nTriggerSavedForMixing{"nTriggerSavedForMixing", 8192, "number of triggers that are saved for mixing with other events"};
  Configurable<int> nTriggerBinMinThreshold{"nTriggerBinMinThreshold", 128, "threshold minimum number of triggers in mixing bin required"};
  Configurable<double> iTriggerStartFraction{"iTriggerStartFraction", 0.5, "ratio of bin size for randomised iTrigger start"};
  Configurable<double> nTriggerScaleCoefficient{"nTriggerScaleCoefficient", 2, "coeffcient a in mixing-number pt power scaling"};
  Configurable<double> nTriggerScaleExponent{"nTriggerScaleExponent", 2, "exponent a in mixing-number pt power scaling"};
  Configurable<int> nMixingAt0McTrue{"nMixingAt0McTrue", 4, "number of triggers that are used for mc true mixing"};
  Configurable<int> nMixingAt0Hadron{"nMixingAt0Hadron", 16, "number of triggers that are used for hadron mixing"};
  Configurable<int> nMixingAt0Pipm{"nMixingAt0Pipm", 16, "number of triggers that are used for pipm mixing"};
  Configurable<int> nMixingAt0PhotonPCM{"nMixingAt0PhotonPCM", 256, "number of triggers that are saved for photonPCM mixing"};
  Configurable<int> nMixingAt0H0PCM{"nMixingAt0H0PCM", 1024, "number of triggers that are saved for h0PCM (pi0, eta) mixing"};

  Configurable<int> nNeighboursMixingPhotonPCMPair{"nNeighboursMixingPhotonPCMPair", 32, "number neighbours used for for photonPCM pair mixing"};
  Configurable<std::vector<double>> pi0PCMPeakMassRange{"pi0PCMPeakMassRange", {0.10, 0.15}, "photon-pair mass integration range for pi0PCM"};
  Configurable<std::vector<double>> pi0PCMSideMassRange{"pi0PCMSideMassRange", {0.16, 0.22}, "photon-pair mass integration range outside pi0PCM region"};
  Configurable<std::vector<double>> etaPCMPeakMassRange{"etaPCMPeakMassRange", {0.51, 0.56}, "photon-pair mass integration range for etaPCM"};
  Configurable<std::vector<double>> etaPCMLowSideMassRange{"etaPCMLowSideMassRange", {0.44, 0.48}, "photon-pair mass integration range below etaPCM region"};
  Configurable<std::vector<double>> etaPCMHighSideMassRange{"etaPCMHighSideMassRange", {0.58, 0.62}, "photon-pair mass integration range above etaPCM region"};

  // for histograms
  Configurable<int> nBinsZPv{"nBinsZPv", 28, "number zPv bins in histos"};
  Configurable<int> nBinsMult{"nBinsMult", 64, "number multiplicity bins in histos"};
  Configurable<int> nBinsOccupancy{"nBinsOccupancy", 2000, "number occupancy bins in histos for QA"};

  Configurable<int> nBinsPhi{"nBinsPhi", 72, "number phi bins"};
  Configurable<int> nBinsEta{"nBinsEta", 32, "number eta bins"};
  Configurable<int> nBinsDCAz{"nBinsDCAz", 100, "number DCAz bins"};
  Configurable<int> nBinsMgg{"nBinsMgg", 240, "number mass-photon-pair bins"};

  Configurable<std::vector<double>> binsPtTrig{"binsPtTrig", {0, 5, 10, 25, 50, 100}, "correlation ptTrig bins"};
  Configurable<std::vector<double>> binsPtAssoc{"binsPtAssoc",
                                                {0, 0.2, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10,
                                                 12.5, 15, 17.5, 20, 25, 35, 50, 100},
                                                "correlation ptAssoc bins"};
  Configurable<std::vector<double>> binsDPhi{"binsDPhi",
                                             {0.00 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.05 * constants::math::TwoPI - constants::math::PIHalf, 0.10 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.14 * constants::math::TwoPI - constants::math::PIHalf, 0.17 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.20 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.22 * constants::math::TwoPI - constants::math::PIHalf, 0.24 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.26 * constants::math::TwoPI - constants::math::PIHalf, 0.28 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.30 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.33 * constants::math::TwoPI - constants::math::PIHalf, 0.36 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.40 * constants::math::TwoPI - constants::math::PIHalf, 0.45 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.50 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.55 * constants::math::TwoPI - constants::math::PIHalf, 0.60 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.64 * constants::math::TwoPI - constants::math::PIHalf, 0.68 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.71 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.74 * constants::math::TwoPI - constants::math::PIHalf, 0.76 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.79 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.82 * constants::math::TwoPI - constants::math::PIHalf, 0.86 * constants::math::TwoPI - constants::math::PIHalf,
                                              0.90 * constants::math::TwoPI - constants::math::PIHalf, 0.95 * constants::math::TwoPI - constants::math::PIHalf,
                                              1.00 * constants::math::TwoPI - constants::math::PIHalf},
                                             "correlation bins DeltaPhi"};
  Configurable<std::vector<double>> binsDEta{"binsDEta",
                                             {0 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault,
                                              2 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault, 4 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault,
                                              6 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault, 8 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault,
                                              9.5 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault, 11 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault,
                                              12.5 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault, 14 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault,
                                              15.5 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault, 16.5 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault,
                                              18 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault,
                                              19.5 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault, 21 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault,
                                              22.5 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault, 24 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault,
                                              26 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault, 28 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault,
                                              30 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault, 32 / 32. * 4 * absEtaMaxDefault - 2 * absEtaMaxDefault},
                                             "correlation bins DeltaEta"};
  Configurable<std::vector<double>> binsZPvBinning{"binsZPvBinning",
                                                   {-7, -5, -3, -1, 1, 3, 5, 7},
                                                   "zPv mixing bins"};
  Configurable<std::vector<double>> binsMultBinning{"binsMultBinning",
                                                    {-0.5, 10.5, 15.5, 20.5, 27.5, 42.5},
                                                    "multiplicity mixing bins"};
  Configurable<std::vector<double>> binsZPvBinningMcTrue{"binsZPvBinningMcTrue",
                                                         {-10000, 10000},
                                                         "zPv mixing bins for mc true"};
  Configurable<std::vector<double>> binsMultBinningMcTrue{"binsMultBinningMcTrue",
                                                          {-0.5, 16.5, 24.5, 31.5, 39.5, 64.5},
                                                          "multiplicity mixing bins for mc true"};

  // configurables from other tasks

  double absEtaMax = -1;

  // further variables

  std::mt19937 randomMt{3141u};
  int randomIntInInterval(int const lowEdgeInclusive, int const upEdgeInclusive)
  {
    return lowEdgeInclusive + (randomMt() % (upEdgeInclusive - lowEdgeInclusive + 1)); // speeds up task by factor > 2; modulo bias should not be relevant here
    // return std::uniform_int_distribution<int>(lowEdgeInclusive, upEdgeInclusive)(randomMt);
  }

  // objects to hold histograms
  HistogramRegistry histos{"histogramRegistry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  // ccdb calls
  int64_t noLaterThanCcdb =
    noLaterThanCcdbConfig == -1 ? std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count() : noLaterThanCcdbConfig;
  Service<ccdb::BasicCCDBManager> ccdb{};
  // for mc
  Service<framework::O2DatabasePDG> pdg{};

  // partitions
  SliceCache cache;

  // prepare for per collision slicing
  Preslice<aod::JetTracks> perColTracks = aod::jtrack::collisionId;
  Preslice<aod::Triggers> perColTriggers = aod::corr_particle::jetCollisionId;
  Preslice<aod::Hadrons> perColHadrons = aod::corr_particle::jetCollisionId;
  Preslice<aod::Pipms> perColPipms = aod::corr_particle::jetCollisionId;
  Preslice<aod::PhotonPCMs> perColPhotonPCMs = aod::corr_particle::jetCollisionId;
  Preslice<aod::PhotonPCMPairs> perColPhotonPCMPairs = aod::corr_particle::jetCollisionId;
  Preslice<aod::JetParticles> perColMcParticles = aod::jmcparticle::mcCollisionId;
  Preslice<aod::TriggerParticles> perColTriggerParticles = aod::corr_particle::jetMcCollisionId;

  // combinations binning
  // cumbersome, but still better than having extra configurable or figuring out how to init binningZPvMult later while declaring it here
  std::function<std::vector<double>(std::vector<double> const&, double const)> prependValueToVector =
    [](std::vector<double> const& vec, double const value) {
      std::vector<double> resultVec = {value};
      resultVec.insert(resultVec.end(), vec.begin(), vec.end());
      return resultVec;
    };
  BinningZPvMult binningZPvMult{{prependValueToVector(binsZPvBinning.value, VARIABLE_WIDTH), prependValueToVector(binsMultBinning.value, VARIABLE_WIDTH)}, true};

  // declare analysis variables

  // correction histograms
  TH1D* h1PtCorrectionTrigger = nullptr;

  // mixing trigger memory
  struct MixingTrigger {
    MixingTrigger(float pt, float phi, float eta)
      : fPt(pt), fPhi(phi), fEta(eta) {}
    float fPt, fPhi, fEta;
    [[nodiscard]] float pt() const { return fPt; }
    [[nodiscard]] float phi() const { return fPhi; }
    [[nodiscard]] float eta() const { return fEta; }
  };
  // class to handle trigger info from previous collisions (beyond single dataframe)
  // organised as zPv- and mult-bin matrix of deque (pt, phi, eta) to save trigger info beyond single dataframe
  // last mult bin open to account for overflow and adjusted zVtx limits to account for rounding errors -> all events accounted for
  // (possibly replace by some advanced derived data method and O2 event mixing in future?)
  class MixingTriggerMemory
  {
   public:
    // finds bin that value belongs to (assumes ordered bin edges) (starts at 0; includes underflow (return -1) and overflow (return binEdges.size() - 1))
    // should be faster than some std binary search due to small number of bins (zPv, mult)
    static int findIntervalBin(double value, const std::vector<double>& binEdges)
    {
      if (binEdges.empty()) {
        throw std::invalid_argument("binEdges in findIntervalBin");
      }
      const int n = binEdges.size() - 1;
      if (value < binEdges[0]) {
        return -1; // underflow
      }
      for (int i_bin = 0; i_bin < n; i_bin++) {
        if (value < binEdges[i_bin + 1]) {
          return i_bin;
        }
      }
      return n; // overflow
    }

    MixingTriggerMemory(std::vector<double> const& binEdgesZPvIn, std::vector<double> const& binEdgesMultIn, size_t const nTriggerPerBinLimitIn)
      : binEdgesZPv(binEdgesZPvIn),
        nBinsZPv(binEdgesZPvIn.size() - 1),
        binEdgesMult(binEdgesMultIn),
        nBinsMult(binEdgesMultIn.size() - 1),
        nTriggerPerBinLimit(nTriggerPerBinLimitIn)
    {
      int const nEdgesMin = 2;
      if (binEdgesZPv.size() < nEdgesMin || binEdgesMult.size() < nEdgesMin) {
        throw std::invalid_argument("too few bin edges to define zPv or mult bins");
      }
      if (!std::is_sorted(binEdgesZPv.begin(), binEdgesZPv.end()) || !std::is_sorted(binEdgesMult.begin(), binEdgesMult.end())) {
        throw std::invalid_argument("bin edges must be sorted");
      }
      // prevent zPv over/underflow due to rounding errors
      binEdgesZPv.front() = binEdgesZPv.front() < 0 ? binEdgesZPv.front() * zPvRoundingErrorAdjust : binEdgesZPv.front() / zPvRoundingErrorAdjust;
      binEdgesZPv.back() = binEdgesZPv.back() > 0 ? binEdgesZPv.back() * zPvRoundingErrorAdjust : binEdgesZPv.back() / zPvRoundingErrorAdjust;
      // init correct size of zPv-mult matrix
      savedTriggersZPvMult.resize(nBinsZPv);
      for (size_t i_zPv = 0; i_zPv < nBinsZPv; i_zPv++) {
        savedTriggersZPvMult[i_zPv].resize(nBinsMult);
      }
    }

    // save trigger for mixing
    // up to nTriggerPerBinLimit stored in rolling buffer for each zPv-mult bin
    void saveTrigger(float const pt, float const phi, float const eta, double const zPv, double const mult)
    {
      int const iBinCorrZPv = getZPvBin(zPv);
      int const iBinCorrMult = getMultBin(mult);
      // fill
      std::deque<MixingTrigger>& triggerBin = savedTriggersZPvMult[iBinCorrZPv][iBinCorrMult];
      triggerBin.push_front(MixingTrigger{pt, phi, eta});
      if (triggerBin.size() > nTriggerPerBinLimit) {
        triggerBin.pop_back();
        return; // skip nTriggerBinMin update when unnecessary
      }
      // update nTriggerBinMin
      nTriggerBinMin = nTriggerPerBinLimit;
      for (size_t i_zPv = 0; i_zPv < nBinsZPv; i_zPv++) {
        for (size_t i_mult = 0; i_mult < nBinsMult; i_mult++) {
          nTriggerBinMin = std::min(nTriggerBinMin, savedTriggersZPvMult[i_zPv][i_mult].size());
        }
      }
    }

    // return deques of trigger pt, phi, eta in the given zPv/mult bin
    [[nodiscard]] std::deque<MixingTrigger> const& getTriggers(double const zPv, double const mult) const
    {
      int const iBinCorrZPv = getZPvBin(zPv);
      int const iBinCorrMult = getMultBin(mult);
      return savedTriggersZPvMult[iBinCorrZPv][iBinCorrMult];
    }

    [[nodiscard]] size_t getBinSizeMin() const
    {
      return nTriggerBinMin;
    }

   private:
    double zPvRoundingErrorAdjust = 1.0001;
    std::vector<double> binEdgesZPv;
    size_t nBinsZPv;
    std::vector<double> binEdgesMult;
    size_t nBinsMult;
    std::vector<std::vector<std::deque<MixingTrigger>>> savedTriggersZPvMult;
    size_t nTriggerPerBinLimit;
    size_t nTriggerBinMin = 0;

    [[nodiscard]] size_t getZPvBin(double const zPv) const
    {
      int const iBinInit = findIntervalBin(zPv, binEdgesZPv);
      if (iBinInit == -1 || iBinInit == static_cast<int>(nBinsZPv)) {
        throw std::runtime_error("zPv underflow or overflow in MixingTriggerMemory");
      }
      return iBinInit;
    }
    [[nodiscard]] size_t getMultBin(double const mult) const
    {
      int const iBinInit = findIntervalBin(mult, binEdgesMult);
      if (iBinInit == static_cast<int>(nBinsMult)) {
        return iBinInit - 1;
      }
      if (iBinInit == -1) {
        throw std::runtime_error("mult underflow in MixingTriggerMemory");
      }
      return iBinInit;
    }
  };

  // different buffers for mixing modes
  enum class TriggerMemoryMode { Reco,
                                 True,
                                 TrueAssocEv,
                                 TrueAssocEvRecoPtTrig };

  std::unique_ptr<MixingTriggerMemory> triggerMemoryReco;
  std::unique_ptr<MixingTriggerMemory> triggerMemoryTrue;
  std::unique_ptr<MixingTriggerMemory> triggerMemoryTrueAssocEv;
  std::unique_ptr<MixingTriggerMemory> triggerMemoryTrueAssocEvRecoPtTrig;

  // functions ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ==========================================================================================================================================================================================

  // general (kenobi) /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // get histograms from ccdb
  // save corrections from ccdb in histogram registry
  void initCcdbHistograms()
  {
    // trigger
    h1PtCorrectionTrigger = nullptr;
    if (applyCTotTrigger) {
      h1PtCorrectionTrigger = ccdb->getForTimeStamp<TH1D>(pathCcdbCorrection.value + "/trigger", noLaterThanCcdb);

      const double* correctionBinsTrigger = h1PtCorrectionTrigger->GetXaxis()->GetXbins()->GetArray();
      const AxisSpec axisPtCorrectionTrigger{std::vector<double>(correctionBinsTrigger, correctionBinsTrigger + h1PtCorrectionTrigger->GetNbinsX() + 1), "#it{p}_{T}"};
      histos.add("usedCorrection/h1_pt_correction_trigger_ccdb", "h1_pt_correction_trigger_ccdb", kTH1D, {axisPtCorrectionTrigger}, true);
      for (int iBin = 1; iBin <= h1PtCorrectionTrigger->GetNbinsX(); iBin++) {
        histos.get<TH1>(HIST("usedCorrection/h1_pt_correction_trigger_ccdb"))->SetBinContent(iBin, h1PtCorrectionTrigger->GetBinContent(iBin));
        histos.get<TH1>(HIST("usedCorrection/h1_pt_correction_trigger_ccdb"))->SetBinError(iBin, h1PtCorrectionTrigger->GetBinError(iBin));
      }
    }
  }

  // create histograms
  void initHistograms()
  {
    // define axes
    const AxisSpec axisN{1, 0., 1., "#it{N}_{something}"};
    const AxisSpec axisCategories{16, 0., 16., "categories"};

    const AxisSpec axisZPv{nBinsZPv, -7, 7, "#it{z}_{pv}"};
    const AxisSpec axisMult{nBinsMult + 1, -0.5, nBinsMult + 0.5, "multiplicity"};
    const AxisSpec axisOccupancy{nBinsOccupancy + 1, -0.5, nBinsOccupancy + 0.5, "occupancy"};

    const AxisSpec axisPhi{nBinsPhi, 0, constants::math::TwoPI, "#it{#varphi}"};
    const AxisSpec axisEta{nBinsEta, -absEtaMax, absEtaMax, "#it{#eta}"};
    const AxisSpec axisDCAz{nBinsDCAz, -5, 5, "DCA_{z}"};
    const AxisSpec axisMgg{nBinsMgg, 0, 0.8, "#it{m}_{#gamma#gamma}"};

    const AxisSpec axisPtTrig{binsPtTrig, "#it{p}_{T}^{trig}"};
    const AxisSpec axisPtAssoc{binsPtAssoc, "#it{p}_{T}^{assoc}"};
    const AxisSpec axisDPhi{binsDPhi, "#Delta#it{#varphi}"};
    const AxisSpec axisDEta{binsDEta, "#Delta#it{#eta}"};
    const AxisSpec axisZPvBinning{binsZPvBinning, "#it{z}_{pv} correlation binning"};
    const AxisSpec axisMultBinning{binsMultBinning, "multiplicity correlation binning"};
    const AxisSpec axisZPvBinningMcTrue{binsZPvBinningMcTrue, "#it{z}_{pv} correlation binning for mc true"};
    const AxisSpec axisMultBinningMcTrue{binsMultBinningMcTrue, "multiplicity correlation binning for mc true"};

    // reco info
    histos.add("reco/info/h1_nEvents", "h1_nEvents", kTH1D, {axisCategories}, true);
    histos.add("reco/info/h3_ptTrigZPvMult", "h3_ptTrigZPvMult", kTHnSparseD, {axisPtTrig, axisZPv, axisMult}, true);
    histos.add("reco/info/h2_ptTrigOccupancy", "h2_ptTrigOccupancy", kTHnSparseD, {axisPtTrig, axisOccupancy}, true);

    // reco (correlation) analysis
    histos.add("reco/corr/h3_ptPhiEta_trig", "h3_ptPhiEta_trig", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);

    auto const add_corrHists =
      [&](std::string const& name_id) {
        histos.add(std::format("reco/corr/h4_ptTrigPtAssocPhiEta_assoc_{}", name_id).data(), std::format("h4_ptTrigPtAssocPhiEta_assoc_{}", name_id).data(),
                   kTHnSparseD, {axisPtTrig, axisPtAssoc, axisPhi, axisEta}, true);
        histos.add(std::format("reco/corr/h6_corr_{}", name_id).data(), std::format("h6_corr_{}", name_id).data(),
                   kTHnSparseF, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
        histos.add(std::format("reco/corr/h6_mix_{}", name_id).data(), std::format("h6_mix_{}", name_id).data(),
                   kTHnSparseF, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
      };
    // hadron
    histos.add("reco/plain/h3_ptPhiEta_hadron", "h3_ptPhiEta_hadron", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
    add_corrHists("hadron");
    // pipm
    histos.add("reco/plain/h3_ptPhiEta_pipm", "h3_ptPhiEta_pipm", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
    add_corrHists("pipm");
    // photonPCM
    histos.add("reco/plain/h3_ptPhiEta_photonPCM", "h3_ptPhiEta_photonPCM", kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
    add_corrHists("photonPCM");
    // photonPCM pairs
    histos.add("reco/plain/h5_ptTrigPtAssocMggZPvMult_photonPCMPair", "h5_ptTrigPtAssocMggZPvMult_photonPCMPair",
               kTHnSparseD, {axisPtTrig, axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);
    histos.add("reco/corr/h5_ptTrigPtAssocMggZPvMult_assoc_photonPCMPair", "h5_ptTrigPtAssocMggZPvMult_assoc_photonPCMPair",
               kTHnSparseD, {axisPtTrig, axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);
    // pi0PCM
    add_corrHists("pi0PCMPeak");
    add_corrHists("pi0PCMSide");
    // etaPCM
    add_corrHists("etaPCMPeak");
    add_corrHists("etaPCMSide");

    // event mixing for photon pairs
    histos.add("reco/plain/h2_zPvMult_photonPCMPair_evMix", "h2_zPvMult_photonPCMPair_evMix", kTHnSparseD, {axisZPv, axisMult}, true);
    histos.add("reco/plain/h4_ptMggZPvMult_photonPCMPair_evMix", "h4_ptMggZPvMult_photonPCMPair_evMix", kTHnSparseD, {axisPtAssoc, axisMgg, axisZPvBinning, axisMultBinning}, true);

    // mc info
    histos.add("mc/info/h3_ptTrigZPvMult_true", "h3_ptTrigZPvMult_true", kTHnSparseD, {axisPtTrig, axisZPv, axisMult}, true);
    histos.add("mc/info/h3_ptTrigZPvMult_trueAssocEv", "h3_ptTrigZPvMult_trueAssocEv", kTHnSparseD, {axisPtTrig, axisZPv, axisMult}, true);
    histos.add("mc/info/h3_ptTrigZPvMult_trueAssocEvRecoPtTrig", "h3_ptTrigZPvMult_trueAssocEvRecoPtTrig", kTHnSparseD, {axisPtTrig, axisZPv, axisMult}, true);

    // reco and true collision correlations
    const std::vector<std::string> assocMcCorrNamesMcAll = {"hadron", "pipm", "photon", "pi0", "eta"};
    for (auto const& correlationType : {"true", "trueAssocEv", "trueAssocEvRecoPtTrig"}) {
      histos.add(std::format("mc/corr/h3_ptPhiEta_trig_{}", correlationType).data(), std::format("h3_ptPhiEta_trig_{}", correlationType).data(),
                 kTHnSparseD, {axisPtAssoc, axisPhi, axisEta}, true);
      for (auto const& assocName : assocMcCorrNamesMcAll) {
        histos.add(std::format("mc/corr/h6_corr_{}_{}", correlationType, assocName).data(), std::format("h6_corr_{}_{}", correlationType, assocName).data(),
                   kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinningMcTrue, axisMultBinningMcTrue}, true);
        histos.add(std::format("mc/corr/h6_mix_{}_{}", correlationType, assocName).data(), std::format("h6_mix_{}_{}", correlationType, assocName).data(),
                   kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinningMcTrue, axisMultBinningMcTrue}, true);
      }
    }

    // decay correlation extra info (just true level)
    const std::vector<std::string> assocMcCorrNamesMcDecayAddition = {"photonDecay", "photonDirect", "photonPi0", "photonEta",
                                                                      "omega", "photonOmega", "etaPrime", "photonEtaPrime", "photonOtherMother"};
    for (auto const& assocName : assocMcCorrNamesMcDecayAddition) {
      histos.add(std::format("mc/corr/h6_corr_true_{}", assocName).data(), std::format("h6_corr_true_{}", assocName).data(),
                 kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinningMcTrue, axisMultBinningMcTrue}, true);
      histos.add(std::format("mc/corr/h6_mix_true_{}", assocName).data(), std::format("h6_mix_true_{}", assocName).data(),
                 kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinningMcTrue, axisMultBinningMcTrue}, true);
    }

    // extra reco correlations with MC info
    auto const addMcRecoCorrHists =
      [&](std::string const& correlationType, std::string const& assocName) {
        histos.add(std::format("mc/corr/h6_corr_{}_{}", correlationType, assocName).data(), std::format("h6_corr_{}_{}", correlationType, assocName).data(),
                   kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
        histos.add(std::format("mc/corr/h6_mix_{}_{}", correlationType, assocName).data(), std::format("h6_mix_{}_{}", correlationType, assocName).data(),
                   kTHnSparseD, {axisDPhi, axisDEta, axisPtTrig, axisPtAssoc, axisZPvBinning, axisMultBinning}, true);
      };
    // matchable associated particles
    const std::vector<std::string> assocCorrNamesMcExtra = {"hadron", "pipm", "photonPCM"};
    for (auto const& correlationType : {"recoAssocEv", "recoPure", "recoPureTruePtAssoc"}) {
      for (auto const& assocName : assocCorrNamesMcExtra) {
        addMcRecoCorrHists(correlationType, assocName);
      }
    }
    // invariant mass associated particles
    // mc pseudo yield
    addMcRecoCorrHists("pseudoReco", "pi0PCM");
    addMcRecoCorrHists("pseudoRecoAssocEv", "pi0PCM");
    addMcRecoCorrHists("pseudoRecoPure", "pi0PCM");
    addMcRecoCorrHists("pseudoRecoPureTruePtAssoc", "pi0PCM");
    // efficiency advanced info
    addMcRecoCorrHists("trueGeoAcc", "pi0");
    addMcRecoCorrHists("trueMeasDecay", "pi0");

    // mc correction
    auto const addCorrectionHistsResolved =
      [&](std::string const& nameId) {
        histos.add(std::format("mc/correction/resol/h4_ptTrigPtAssocPhiEta_{}", nameId).data(), std::format("h4_ptTrigPtAssocPhiEta_{}", nameId).data(),
                   kTHnSparseD, {axisPtTrig, axisPtAssoc, axisPhi, axisEta}, true);
        histos.add(std::format("mc/correction/resol/h4_ptTrigPtAssocZPvMult_{}", nameId).data(), std::format("h4_ptTrigPtAssocZPvMult_{}", nameId).data(),
                   kTHnSparseD, {axisPtTrig, axisPtAssoc, axisZPv, axisMult}, true);
      };
    auto const addCorrectionHistsResolvedBasicTypes =
      [&](std::string const& nameReco, std::string const& nameTrue) {
        std::vector<std::string> const recoTypes = {"reco", "recoAssocEv", "recoPure"};
        std::vector<std::string> const trueTypes = {"true", "trueAssocEv", "trueAssocEvRecoPtTrig"};
        for (std::string const& recoType : recoTypes) {
          addCorrectionHistsResolved(std::format("{}_{}", recoType, nameReco));
        }
        for (std::string const& trueType : trueTypes) {
          addCorrectionHistsResolved(std::format("{}_{}", trueType, nameTrue));
        }
      };

    auto const addCorrectionHistPt =
      [&](std::string const& nameId) {
        histos.add(std::format("mc/correction/h2_ptTrigPtAssoc_{}", nameId).data(), std::format("h2_ptTrigPtAssoc_{}", nameId).data(),
                   kTHnSparseD, {axisPtTrig, axisPtAssoc}, true);
      };
    auto const addCorrectionHistPtCategories =
      [&](std::string const& nameId) {
        histos.add(std::format("mc/correction/h3_ptTrigPtAssocCategory_{}", nameId).data(), std::format("h3_ptTrigPtAssocCategory_{}", nameId).data(),
                   kTHnSparseD, {axisPtTrig, axisPtAssoc, axisCategories}, true);
      };
    auto const addCorrectionHistsV0 =
      [&](std::string const& nameReco, std::string const& nameTrue) {
        addCorrectionHistPt(std::format("true_{}", nameTrue));
        addCorrectionHistPt(std::format("trueAssocEv_{}", nameTrue));
        addCorrectionHistPt(std::format("trueAssocEvRecoPtTrig_{}", nameTrue));
        // mc pseudo yield
        addCorrectionHistPt(std::format("pseudoReco_{}", nameReco));
        addCorrectionHistPt(std::format("pseudoRecoAssocEv_{}", nameReco));
        addCorrectionHistPt(std::format("pseudoRecoPure_{}", nameReco));
        // efficiency advanced info
        addCorrectionHistPt(std::format("trueGeoAcc_{}", nameTrue));
        addCorrectionHistPt(std::format("trueMeasDecay_{}", nameTrue));
      };

    // matchables
    // tracks
    addCorrectionHistsResolvedBasicTypes("hadron", "hadron");
    // pipm PID
    addCorrectionHistsResolvedBasicTypes("pipm", "pipm");
    // photonPCM
    addCorrectionHistsResolvedBasicTypes("photonPCM", "photon");
    // impurities
    addCorrectionHistPtCategories("recoImpurities_photonPCM");

    // invariant mass reconstruction
    // pi0
    addCorrectionHistsV0("pi0PCM", "pi0");
    // impurities
    addCorrectionHistPtCategories("pseudoRecoImpurities_pi0PCM");
    // eta
    addCorrectionHistsV0("etaPCM", "eta");

    // tests

    // DCAz
    histos.add("reco/plain/h5_ptTrigPtAssocDCAzZPvMult_photonPCM", "h5_ptTrigPtAssocDCAzZPvMult_photonPCM",
               kTHnSparseD, {axisPtTrig, axisPtAssoc, axisDCAz, axisZPvBinning, axisMultBinning}, true);
    histos.add("mc/plain/h5_ptTrigPtAssocDCAzZPvMult_photonPCM", "h5_ptTrigPtAssocDCAzZPvMult_photonPCM",
               kTHnSparseD, {axisPtTrig, axisPtAssoc, axisDCAz, axisZPvBinning, axisMultBinning}, true);
  }

  // selections ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // mc split event selection based in the thrid decimal of mc-true posZ
  template <typename T_collision>
  bool checkSplitMcEventSelection(T_collision const& collision)
  {
    // select based on configurables
    return collision.thirdDecimalTruePosZ() % splitMcEvN == splitMcEvSelect - 1;
  }

  // total event selection (basic selections from producer and analysis level selections)
  template <typename T_collision>
  bool totalEvSel(T_collision const& collision)
  {
    bool isSelEv = true;

    if constexpr (requires { collision.selEv(); }) {
      isSelEv = collision.selEv();
    }
    return isSelEv && checkSplitMcEventSelection(collision);
  }

  // checks if mcParticle is charged
  template <typename T_mcParticle>
  bool checkChargedMc(T_mcParticle const& mcParticle)
  {
    auto const pdgParticle = pdg->GetParticle(mcParticle.pdgCode());
    return pdgParticle && pdgParticle->Charge() != 0;
  }
  // checks if fast decaying mcParticle is 'primary'
  template <typename T_mcParticle>
  bool checkDecayPrimary(T_mcParticle const& mcParticle)
  {
    // identify decaying primary
    return mcParticle.producedByGenerator();
  }
  // checks if mcParticle daughters are two photons
  template <typename T_mcParticle>
  bool checkToGG(T_mcParticle const& mcParticle)
  {
    auto const& daughters = mcParticle.template daughters_as<aod::JetParticles>();
    constexpr int NDaughtersToGG = 2;
    if (daughters.size() != NDaughtersToGG) {
      return false;
    }
    auto daughter = daughters.begin();
    auto daughterEnd = daughters.end();
    while (daughter != daughterEnd) {
      if ((*daughter).pdgCode() != PDG_t::kGamma) {
        return false;
      }
      daughter++;
    }
    return true;
  }
  // check if particle has mother in parent tree
  template <typename T_mcParticle>
  bool checkForMother(T_mcParticle const& mcParticle, int const pdgCode, bool const checkAntiParticle)
  {
    if (!mcParticle.has_mothers()) {
      return false;
    }
    auto const mothers = mcParticle.template mothers_as<aod::JetParticles>();
    auto mother = mothers.begin();
    auto motherEnd = mothers.end();
    while (mother != motherEnd) {

      // LOGF(info, "searchPdgCode: %i, current[ pdgCode: %i, status: %i, primary: %i ], mother[ pdgCode: %i, status: %i, primary: %i ]",
      //      pdgCode,
      //      mcParticle.pdgCode(), mcParticle.getGenStatusCode(), mcParticle.isPhysicalPrimary(),
      //      (*mother).pdgCode(), (*mother).getGenStatusCode(), (*mother).isPhysicalPrimary());

      if ((*mother).pdgCode() == pdgCode || (checkAntiParticle && (*mother).pdgCode() == -pdgCode)) {
        return true;
      }
      if (checkForMother(*mother, pdgCode, checkAntiParticle)) {
        return true;
      }
      mother++;
    }
    return false;
  }

  // checks if daughters are in acceptance
  template <typename T_mcParticle>
  bool checkDaughtersInAcceptance(T_mcParticle const& mcParticle, std::pair<double, double> const etaRange, std::pair<double, double> const phiRange = {0, constants::math::TwoPI})
  {
    auto const& daughterPhotons = mcParticle.template daughters_as<aod::JetParticles>();
    auto daughterPhoton = daughterPhotons.begin();
    auto daughterPhotonEnd = daughterPhotons.end();
    while (daughterPhoton != daughterPhotonEnd) {
      if ((*daughterPhoton).eta() < etaRange.first || (*daughterPhoton).eta() > etaRange.second) {
        return false;
      }
      if ((*daughterPhoton).phi() < phiRange.first || (*daughterPhoton).phi() > phiRange.second) {
        return false;
      }
      daughterPhoton++;
    }
    return true;
  }

  // checks if tracks come from photon conversion
  template <typename T_track>
  bool isConversionPhoton(T_track const& posTrack, T_track const& negTrack)
  {
    // check same mother
    auto const& posMothers = posTrack.mcParticle().template mothers_as<aod::JetParticles>();
    auto const& negMothers = negTrack.mcParticle().template mothers_as<aod::JetParticles>();
    if (posMothers.size() != 1 || negMothers.size() != 1) {
      return false;
    }
    if (posMothers.begin()->globalIndex() != negMothers.begin()->globalIndex()) {
      return false;
    }
    // check photon
    if (posMothers.begin()->pdgCode() != PDG_t::kGamma) {
      return false;
    }

    return true;
  };
  // checks if tracks come from double conversion
  template <typename T_track>
  bool isGGFromDoubleConversion(T_track const& posTrack1, T_track const& negTrack1, T_track const& posTrack2, T_track const& negTrack2, int const pdgCode)
  {
    if (!isConversionPhoton(posTrack1, negTrack1) || !isConversionPhoton(posTrack2, negTrack2)) {
      return false;
    }
    // check same mother
    auto const& mothers1 = (*(posTrack1.mcParticle().template mothers_as<aod::JetParticles>().begin())).template mothers_as<aod::JetParticles>();
    auto const& mothers2 = (*(posTrack2.mcParticle().template mothers_as<aod::JetParticles>().begin())).template mothers_as<aod::JetParticles>();
    constexpr int NMothersPhotonFromH0 = 2; // for some reason two mothers (same particle) for h0 decays (contradicts PYTHIA documentation, but whatever)
    if (mothers1.size() != NMothersPhotonFromH0 || mothers2.size() != NMothersPhotonFromH0) {
      return false;
    }
    if (mothers1.begin()->globalIndex() != mothers2.begin()->globalIndex()) {
      return false;
    }
    // check particle type
    if (mothers1.begin()->pdgCode() != pdgCode) {
      return false;
    }

    return true;
  };

  // analysis helpers /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  template <typename T_h1>
  double getH1ValueAt(T_h1 const* const h1, double const value)
  {
    return h1->GetBinContent(h1->FindFixBin(value));
  }
  // correction helpers
  enum class CorrectionParticleType { Trigger,
                                      Hadron,
                                      Pipm,
                                      PhotonPCM };

  template <CorrectionParticleType T_correctionParticleType>
  double getPtCorrection(double const value)
  {
    if constexpr (T_correctionParticleType == CorrectionParticleType::Trigger) {
      return applyCTotTrigger ? getH1ValueAt(h1PtCorrectionTrigger, value) : 1;
    } else {
      return 1;
    }
  }

  // performs 'phi1 - phi2' and pushes it into the interval [-pi/2, 3pi/2]
  inline double getDeltaPhi(double const phi1, double const phi2)
  {
    return RecoDecay::constrainAngle(phi1 - phi2, -1 * constants::math::PIHalf);
  }

  // checks that two values belong to the same category (assumes ordered bins)
  // returns -1 for negative result (also for under/overflow values) and bin number (starting at 0) otherwise
  int checkSameBin(double const value1, double const value2, std::vector<double> const& bins)
  {
    // reject underflow
    if (value1 < bins[0]) {
      return -1;
    }
    // loop over bins
    const int n = bins.size() - 1;
    for (int i_bin = 0; i_bin < n; i_bin++) {
      if (value1 < bins[i_bin + 1]) {
        if (value2 < bins[i_bin + 1] && value2 >= bins[i_bin]) {
          return i_bin;
        }
        return -1;
      }
    }
    // reject overflow
    return -1;
  }

  // check if invariant mass range
  enum class MassRange { Pi0PCMPeak,
                         Pi0PCMSide,
                         EtaPCMPeak,
                         EtaPCMSide };

  template <MassRange range>
  bool checkMassRange(double const mgg)
  {
    if constexpr (range == MassRange::Pi0PCMPeak) {
      return mgg > pi0PCMPeakMassRange.value[0] && mgg < pi0PCMPeakMassRange.value[1];
    } else if constexpr (range == MassRange::Pi0PCMSide) {
      return mgg > pi0PCMSideMassRange.value[0] && mgg < pi0PCMSideMassRange.value[1];
    } else if constexpr (range == MassRange::EtaPCMPeak) {
      return mgg > etaPCMPeakMassRange.value[0] && mgg < etaPCMPeakMassRange.value[1];
    } else if constexpr (range == MassRange::EtaPCMSide) {
      return (mgg > etaPCMLowSideMassRange.value[0] && mgg < etaPCMLowSideMassRange.value[1]) ||
             (mgg > etaPCMHighSideMassRange.value[0] && mgg < etaPCMHighSideMassRange.value[1]);
    }
    return false;
  }

  // analysis /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // generalised correlation functions
  // per collision

  // plain info
  template <typename T_collision, typename T_associatedThisEvent,
            typename T_funcPlain>
  void corrProcessPlain(T_collision const& collision, T_associatedThisEvent const& associatedThisEvent,
                        T_funcPlain const& funcPlain)
  {
    // normal spectra (per event - not per trigger)
    for (auto const& associated : associatedThisEvent) {
      funcPlain(collision, associated);
    }
  }

  // correlation
  template <typename T_collision, typename T_triggersThisEvent, typename T_associatedThisEvent,
            typename T_funcCorrelation>
  void corrProcessCorrelation(T_collision const& collision, T_triggersThisEvent const& triggersThisEvent, T_associatedThisEvent const& associatedThisEvent,
                              T_funcCorrelation const& funcCorrelation)
  {
    // correlation combinations
    for (auto const& [trigger, associated] : soa::combinations(soa::CombinationsFullIndexPolicy(triggersThisEvent, associatedThisEvent))) {
      funcCorrelation(collision, trigger, associated);
    }
  }

  // mixing-number pt scaling with power law
  size_t nMixingPtPowerScaling(double const pt, size_t const n0)
  {
    double const rawScale = 1 + nTriggerScaleCoefficient * std::pow(pt, nTriggerScaleExponent);
    return n0 * static_cast<size_t>(rawScale);
  }

  // mixing
  template <TriggerMemoryMode triggerMemoryMode,
            typename T_collision, typename T_associatedThisEvent,
            typename T_funcMixing>
  void corrProcessMixing(T_collision const& collision, T_associatedThisEvent const& associatedThisEvent,
                         T_funcMixing const& funcMixing,
                         size_t const nTriggerMixingAt0)
  {
    size_t triggerMemoryBinSizeMin = 0;

    // mixing loops (more dynamic than O2 mixing)
    std::deque<MixingTrigger> const& savedTriggers =
      [&]() -> std::deque<MixingTrigger> const& {
      if constexpr (triggerMemoryMode == TriggerMemoryMode::Reco) {
        triggerMemoryBinSizeMin = triggerMemoryReco->getBinSizeMin();
        return triggerMemoryReco->getTriggers(collision.posZ(), collision.nGlobalTracks());
      } else if constexpr (triggerMemoryMode == TriggerMemoryMode::True) {
        triggerMemoryBinSizeMin = triggerMemoryTrue->getBinSizeMin();
        return triggerMemoryTrue->getTriggers(collision.posZ(), collision.nChargedInEtaRange());
      } else if constexpr (triggerMemoryMode == TriggerMemoryMode::TrueAssocEv) {
        triggerMemoryBinSizeMin = triggerMemoryTrueAssocEv->getBinSizeMin();
        return triggerMemoryTrueAssocEv->getTriggers(collision.posZ(), collision.nChargedInEtaRange());
      } else if constexpr (triggerMemoryMode == TriggerMemoryMode::TrueAssocEvRecoPtTrig) {
        triggerMemoryBinSizeMin = triggerMemoryTrueAssocEvRecoPtTrig->getBinSizeMin();
        return triggerMemoryTrueAssocEvRecoPtTrig->getTriggers(collision.posZ(), collision.nChargedInEtaRange());
      }
    }();
    // require enough mixing statistics
    if (triggerMemoryBinSizeMin < static_cast<size_t>(nTriggerBinMinThreshold)) {
      return;
    }
    const size_t nSavedTriggers = savedTriggers.size();
    const size_t iStartTriggerParticle = randomIntInInterval(0, iTriggerStartFraction * nSavedTriggers);
    // associated loop
    for (auto const& associated : associatedThisEvent) {
      // trigger numbers
      const size_t nTriggerMixing = nMixingPtPowerScaling(associated.pt(), nTriggerMixingAt0);
      const size_t mixUpToTriggerN = std::min(nSavedTriggers, iStartTriggerParticle + nTriggerMixing);
      const float perTriggerWeight = 1. / (mixUpToTriggerN - iStartTriggerParticle); // mixUpToTriggerN <= iStartTriggerParticle caught by break statement of loop
      // trigger loop
      for (size_t i_mixingTrigger = iStartTriggerParticle; i_mixingTrigger < mixUpToTriggerN; i_mixingTrigger++) {
        MixingTrigger const& mixingTrigger = savedTriggers[i_mixingTrigger];
        funcMixing(collision, mixingTrigger, associated, perTriggerWeight);
      }
    }
  }

  // analysis tasks ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // =/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/=/

  void init(InitContext& initContext)
  {
    // analysis info
    ccdb->setURL(urlCcdb);
    // enabling object caching (otherwise each call goes to CCDB server)
    ccdb->setCaching(true);
    // ccdb->setLocalObjectValidityChecking();
    // not later than (avoids replacing objects while a train is running)
    ccdb->setCreatedNotAfter(noLaterThanCcdb);

    // init analysis variables

    // get variables from other tasks
    o2::common::core::getTaskOptionValue(initContext, "photon-charged-trigger-producer", "absEtaMax", absEtaMax, false);

    // histograms from ccdb
    initCcdbHistograms();

    // create analysis histograms
    initHistograms();

    // mixing trigger memory
    triggerMemoryReco = std::make_unique<MixingTriggerMemory>(binsZPvBinning, binsMultBinning, nTriggerSavedForMixing);
    triggerMemoryTrue = std::make_unique<MixingTriggerMemory>(binsZPvBinningMcTrue, binsMultBinningMcTrue, nTriggerSavedForMixing);
    triggerMemoryTrueAssocEv = std::make_unique<MixingTriggerMemory>(binsZPvBinningMcTrue, binsMultBinningMcTrue, nTriggerSavedForMixing);
    triggerMemoryTrueAssocEvRecoPtTrig = std::make_unique<MixingTriggerMemory>(binsZPvBinningMcTrue, binsMultBinningMcTrue, nTriggerSavedForMixing);
  }

  // reconstructed ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ==========================================================================================================================================================================================

  void processCorrHadron(CorrCollisions const& collisions, aod::Triggers const& triggers, aod::Hadrons const& hadrons)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const hadronsThisEvent = hadrons.sliceBy(perColHadrons, collision.globalIndex());

      auto const funcPlain = [this]([[maybe_unused]] auto const& collision, auto const& associated) {
        histos.fill(HIST("reco/plain/h3_ptPhiEta_hadron"),
                    associated.pt(), associated.phi(), associated.eta());
      };
      corrProcessPlain(collision, hadronsThisEvent, funcPlain);

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // exclude self correlation
        if (trigger.jetTrackId() == associated.jetTrackId()) {
          return;
        }

        histos.fill(HIST("reco/corr/h4_ptTrigPtAssocPhiEta_assoc_hadron"),
                    collision.ptMax(), associated.pt(), associated.phi(), associated.eta(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        histos.fill(HIST("reco/corr/h6_corr_hadron"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
      };
      corrProcessCorrelation(collision, triggersThisEvent, hadronsThisEvent, funcCorrelation);

      // select mixing events
      if (doTrigEvMixing && !collision.trigEv()) {
        continue;
      }

      auto const funcMixing = [this](auto const& collision, auto const& trigger, auto const& associated, auto const perTriggerWeight) {
        histos.fill(HIST("reco/corr/h6_mix_hadron"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
      };
      corrProcessMixing<TriggerMemoryMode::Reco>(collision, hadronsThisEvent, funcMixing, nMixingAt0Hadron);
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrHadron, "process standard correlation for associated hardons", false);

  void processCorrPipm(CorrCollisions const& collisions, aod::Triggers const& triggers, aod::Pipms const& pipms)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const pipmsThisEvent = pipms.sliceBy(perColPipms, collision.globalIndex());

      auto const funcPlain = [this]([[maybe_unused]] auto const& collision, auto const& associated) {
        histos.fill(HIST("reco/plain/h3_ptPhiEta_pipm"),
                    associated.pt(), associated.phi(), associated.eta());
      };
      corrProcessPlain(collision, pipmsThisEvent, funcPlain);

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // exclude self correlation
        if (trigger.jetTrackId() == associated.jetTrackId()) {
          return;
        }

        histos.fill(HIST("reco/corr/h4_ptTrigPtAssocPhiEta_assoc_pipm"),
                    collision.ptMax(), associated.pt(), associated.phi(), associated.eta(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        histos.fill(HIST("reco/corr/h6_corr_pipm"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
      };
      corrProcessCorrelation(collision, triggersThisEvent, pipmsThisEvent, funcCorrelation);

      // select mixing events
      if (doTrigEvMixing && !collision.trigEv()) {
        continue;
      }

      auto const funcMixing = [this](auto const& collision, auto const& trigger, auto const& associated, auto const perTriggerWeight) {
        histos.fill(HIST("reco/corr/h6_mix_pipm"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
      };
      corrProcessMixing<TriggerMemoryMode::Reco>(collision, pipmsThisEvent, funcMixing, nMixingAt0Pipm);
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrPipm, "process standard correlation for associated pipm", false);

  void processCorrPhotonPCM(CorrCollisions const& collisions, aod::Triggers const& triggers, aod::PhotonPCMs const& photonPCMs)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const photonPCMsThisEvent = photonPCMs.sliceBy(perColPhotonPCMs, collision.globalIndex());

      auto const funcPlain = [this]([[maybe_unused]] auto const& collision, auto const& associated) {
        histos.fill(HIST("reco/plain/h3_ptPhiEta_photonPCM"),
                    associated.pt(), associated.phi(), associated.eta());
      };
      corrProcessPlain(collision, photonPCMsThisEvent, funcPlain);

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // exclude self correlation
        if (trigger.jetTrackId() == associated.posJetTrackId() || trigger.jetTrackId() == associated.negJetTrackId()) {
          return;
        }

        histos.fill(HIST("reco/corr/h4_ptTrigPtAssocPhiEta_assoc_photonPCM"),
                    collision.ptMax(), associated.pt(), associated.phi(), associated.eta(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        histos.fill(HIST("reco/corr/h6_corr_photonPCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
      };
      corrProcessCorrelation(collision, triggersThisEvent, photonPCMsThisEvent, funcCorrelation);

      // select mixing events
      if (doTrigEvMixing && !collision.trigEv()) {
        continue;
      }

      auto const funcMixing = [this](auto const& collision, auto const& trigger, auto const& associated, auto const perTriggerWeight) {
        histos.fill(HIST("reco/corr/h6_mix_photonPCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
      };
      corrProcessMixing<TriggerMemoryMode::Reco>(collision, photonPCMsThisEvent, funcMixing, nMixingAt0PhotonPCM);
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrPhotonPCM, "process standard correlation for associated photonPCM", false);

  void processCorrPhotonPCMPair(CorrCollisions const& collisions, aod::Triggers const& triggers, aod::PhotonPCMPairs const& photonPCMPairs)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const photonPCMPairsThisEvent = photonPCMPairs.sliceBy(perColPhotonPCMPairs, collision.globalIndex());

      auto const funcPlain = [this](auto const& collision, auto const& associated) {
        histos.fill(HIST("reco/plain/h5_ptTrigPtAssocMggZPvMult_photonPCMPair"), collision.ptMax(), associated.pt(), associated.mgg(), collision.posZ(), collision.nGlobalTracks());
      };
      corrProcessPlain(collision, photonPCMPairsThisEvent, funcPlain);

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // exclude self correlation
        if (trigger.jetTrackId() == associated.posJetTrack1Id() || trigger.jetTrackId() == associated.negJetTrack1Id() ||
            trigger.jetTrackId() == associated.negJetTrack2Id() || trigger.jetTrackId() == associated.posJetTrack2Id()) {
          return;
        }

        histos.fill(HIST("reco/corr/h5_ptTrigPtAssocMggZPvMult_assoc_photonPCMPair"),
                    trigger.pt(), associated.pt(), associated.mgg(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));

        // pi0
        if (checkMassRange<MassRange::Pi0PCMPeak>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h4_ptTrigPtAssocPhiEta_assoc_pi0PCMPeak"),
                      collision.ptMax(), associated.pt(), associated.phi(), associated.eta(),
                      getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
          histos.fill(HIST("reco/corr/h6_corr_pi0PCMPeak"),
                      getDeltaPhi(trigger.phi(), associated.phi()),
                      trigger.eta() - associated.eta(),
                      trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
          return;
        }
        if (checkMassRange<MassRange::Pi0PCMSide>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h4_ptTrigPtAssocPhiEta_assoc_pi0PCMSide"),
                      collision.ptMax(), associated.pt(), associated.phi(), associated.eta(),
                      getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
          histos.fill(HIST("reco/corr/h6_corr_pi0PCMSide"),
                      getDeltaPhi(trigger.phi(), associated.phi()),
                      trigger.eta() - associated.eta(),
                      trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
          return;
        }
        // eta
        if (checkMassRange<MassRange::EtaPCMPeak>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h4_ptTrigPtAssocPhiEta_assoc_etaPCMPeak"),
                      collision.ptMax(), associated.pt(), associated.phi(), associated.eta(),
                      getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
          histos.fill(HIST("reco/corr/h6_corr_etaPCMPeak"),
                      getDeltaPhi(trigger.phi(), associated.phi()),
                      trigger.eta() - associated.eta(),
                      trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
          return;
        }
        if (checkMassRange<MassRange::EtaPCMSide>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h4_ptTrigPtAssocPhiEta_assoc_etaPCMSide"),
                      collision.ptMax(), associated.pt(), associated.phi(), associated.eta(),
                      getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
          histos.fill(HIST("reco/corr/h6_corr_etaPCMSide"),
                      getDeltaPhi(trigger.phi(), associated.phi()),
                      trigger.eta() - associated.eta(),
                      trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
          return;
        }
      };
      corrProcessCorrelation(collision, triggersThisEvent, photonPCMPairsThisEvent, funcCorrelation);

      // select mixing events
      if (doTrigEvMixing && !collision.trigEv()) {
        continue;
      }

      auto const funcMixing = [this](auto const& collision, auto const& trigger, auto const& associated, auto const perTriggerWeight) {
        // pi0
        if (checkMassRange<MassRange::Pi0PCMPeak>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h6_mix_pi0PCMPeak"),
                      getDeltaPhi(trigger.phi(), associated.phi()),
                      trigger.eta() - associated.eta(),
                      trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
          return;
        }
        if (checkMassRange<MassRange::Pi0PCMSide>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h6_mix_pi0PCMSide"),
                      getDeltaPhi(trigger.phi(), associated.phi()),
                      trigger.eta() - associated.eta(),
                      trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
          return;
        }
        // eta
        if (checkMassRange<MassRange::EtaPCMPeak>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h6_mix_etaPCMPeak"),
                      getDeltaPhi(trigger.phi(), associated.phi()),
                      trigger.eta() - associated.eta(),
                      trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
          return;
        }
        if (checkMassRange<MassRange::EtaPCMSide>(associated.mgg())) {
          histos.fill(HIST("reco/corr/h6_mix_etaPCMSide"),
                      getDeltaPhi(trigger.phi(), associated.phi()),
                      trigger.eta() - associated.eta(),
                      trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                      perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
          return;
        }
      };
      corrProcessMixing<TriggerMemoryMode::Reco>(collision, photonPCMPairsThisEvent, funcMixing, nMixingAt0H0PCM);
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrPhotonPCMPair, "process standard correlation for associated pi0PCM", false);

  void processCorrPhotonPCMPairMix(CorrCollisions const& collisions, aod::PhotonPCMs const& photonPCMs)
  {
    auto photonPCMsTuple = std::make_tuple(photonPCMs);
    SameKindPair<CorrCollisions, aod::PhotonPCMs, BinningZPvMult> pairs{binningZPvMult, nNeighboursMixingPhotonPCMPair, -1, collisions, photonPCMsTuple, &cache};

    // mixed events
    for (auto pair = pairs.begin(); pair != pairs.end(); pair++) {
      auto const& [collision1, photonPCMs1, collision2, photonPCMs2] = *pair;

      // // check that current und mixing-trigger event are from the same zPv/mult bins
      // if (checkSameBin(collision1.posZ(), collision2.posZ(), binsZPvBinning) == -1) {
      //   std::printf("ERROR: zPv bins do not match\n"); continue;
      // }
      // if (checkSameBin(collision1.nGlobalTracks(), collision2.nGlobalTracks(), binsMultBinning) == -1) {
      //   std::printf("ERROR: multiplicity bins do not match\n"); continue;
      // }

      // event selection
      if (!totalEvSel(collision1) || !totalEvSel(collision2)) {
        continue;
      }
      // event info
      histos.fill(HIST("reco/plain/h2_zPvMult_photonPCMPair_evMix"), collision1.posZ(), collision1.nGlobalTracks());
      // mixing loop
      for (auto const& [photonPCM1, photonPCM2] : soa::combinations(soa::CombinationsFullIndexPolicy(photonPCMs1, photonPCMs2))) {
        ROOT::Math::PtEtaPhiMVector const p4photonPCM1(photonPCM1.pt(), photonPCM1.eta(), photonPCM1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector const p4photonPCM2(photonPCM2.pt(), photonPCM2.eta(), photonPCM2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector const p4photonPCMPair = p4photonPCM1 + p4photonPCM2;

        histos.fill(HIST("reco/plain/h4_ptMggZPvMult_photonPCMPair_evMix"), p4photonPCMPair.pt(), p4photonPCMPair.M(), collision1.posZ(), collision1.nGlobalTracks());
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processCorrPhotonPCMPairMix, "process gamma-gamma mixing for photonPCM", false);

  void processPCMDCAz(CorrCollision const& collision, aod::PhotonPCMs const& photonPCMs, aod::V0PhotonsKF const&)
  {
    // event selection
    if (!totalEvSel(collision)) {
      return;
    }

    for (auto const& photonPCM : photonPCMs) {
      histos.fill(HIST("reco/plain/h5_ptTrigPtAssocDCAzZPvMult_photonPCM"),
                  collision.ptMax(), photonPCM.pt(), photonPCM.v0PhotonKF().dcaZtopv(), collision.posZ(), collision.nGlobalTracks());
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processPCMDCAz, "process to test DCAz distribution of photons in trigger collisions", false);

  // MC ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ==========================================================================================================================================================================================

  // correlations /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  // (sad) attempt at reducing code duplication
  enum class McCorrTruthLevel : int { True = 0,
                                      TrueAssocEv = 1,
                                      TrueAssocEvRecoPtTrig = 2,
                                      TrueMeasDecay = 3,
                                      TrueGeoAcc = 4 };
  enum class McCorrCorrelationType : int { Correlation = 0,
                                           Mixing = 1 };
  enum class McCorrAssocType : int { Hadron = 0,
                                     Pipm = 1,
                                     Photon = 2,
                                     PhotonDecay = 3,
                                     PhotonDirect = 4,
                                     Pi0 = 5,
                                     PhotonPi0 = 6,
                                     Eta = 7,
                                     PhotonEta = 8,
                                     Omega = 9,
                                     PhotonOmega = 10,
                                     EtaPrime = 11,
                                     PhotonEtaPrime = 12,
                                     PhotonOtherMother = 13 };
  static constexpr std::array<std::array<std::array<char const*, 14>, 2>, 5> McTrueCorrHistPaths = {{{{{"mc/corr/h6_corr_true_hadron", "mc/corr/h6_corr_true_pipm",
                                                                                                        "mc/corr/h6_corr_true_photon", "mc/corr/h6_corr_true_photonDecay", "mc/corr/h6_corr_true_photonDirect",
                                                                                                        "mc/corr/h6_corr_true_pi0", "mc/corr/h6_corr_true_photonPi0", "mc/corr/h6_corr_true_eta", "mc/corr/h6_corr_true_photonEta",
                                                                                                        "mc/corr/h6_corr_true_omega", "mc/corr/h6_corr_true_photonOmega", "mc/corr/h6_corr_true_etaPrime", "mc/corr/h6_corr_true_photonEtaPrime", "mc/corr/h6_corr_true_photonOtherMother"},
                                                                                                       {"mc/corr/h6_mix_true_hadron", "mc/corr/h6_mix_true_pipm",
                                                                                                        "mc/corr/h6_mix_true_photon", "mc/corr/h6_mix_true_photonDecay", "mc/corr/h6_mix_true_photonDirect",
                                                                                                        "mc/corr/h6_mix_true_pi0", "mc/corr/h6_mix_true_photonPi0", "mc/corr/h6_mix_true_eta", "mc/corr/h6_mix_true_photonEta",
                                                                                                        "mc/corr/h6_mix_true_omega", "mc/corr/h6_mix_true_photonOmega", "mc/corr/h6_mix_true_etaPrime", "mc/corr/h6_mix_true_photonEtaPrime", "mc/corr/h6_mix_true_photonOtherMother"}}},
                                                                                                     {{{"mc/corr/h6_corr_trueAssocEv_hadron", "mc/corr/h6_corr_trueAssocEv_pipm",
                                                                                                        "mc/corr/h6_corr_trueAssocEv_photon", "", "",
                                                                                                        "mc/corr/h6_corr_trueAssocEv_pi0", "", "mc/corr/h6_corr_trueAssocEv_eta", "",
                                                                                                        "", "", "", "", ""},
                                                                                                       {"mc/corr/h6_mix_trueAssocEv_hadron", "mc/corr/h6_mix_trueAssocEv_pipm",
                                                                                                        "mc/corr/h6_mix_trueAssocEv_photon", "", "",
                                                                                                        "mc/corr/h6_mix_trueAssocEv_pi0", "", "mc/corr/h6_mix_trueAssocEv_eta", "",
                                                                                                        "", "", "", "", ""}}},
                                                                                                     {{{"mc/corr/h6_corr_trueAssocEvRecoPtTrig_hadron", "mc/corr/h6_corr_trueAssocEvRecoPtTrig_pipm",
                                                                                                        "mc/corr/h6_corr_trueAssocEvRecoPtTrig_photon", "", "",
                                                                                                        "mc/corr/h6_corr_trueAssocEvRecoPtTrig_pi0", "", "mc/corr/h6_corr_trueAssocEvRecoPtTrig_eta", "",
                                                                                                        "", "", "", "", ""},
                                                                                                       {"mc/corr/h6_mix_trueAssocEvRecoPtTrig_hadron", "mc/corr/h6_mix_trueAssocEvRecoPtTrig_pipm",
                                                                                                        "mc/corr/h6_mix_trueAssocEvRecoPtTrig_photon", "", "",
                                                                                                        "mc/corr/h6_mix_trueAssocEvRecoPtTrig_pi0", "", "mc/corr/h6_mix_trueAssocEvRecoPtTrig_eta", "",
                                                                                                        "", "", "", "", ""}}},
                                                                                                     {{{"", "",
                                                                                                        "", "", "",
                                                                                                        "mc/corr/h6_corr_trueMeasDecay_pi0", "", "", "",
                                                                                                        "", "", "", "", ""},
                                                                                                       {"", "",
                                                                                                        "", "", "",
                                                                                                        "mc/corr/h6_mix_trueMeasDecay_pi0", "", "", "",
                                                                                                        "", "", "", "", ""}}},
                                                                                                     {{{"", "",
                                                                                                        "", "", "",
                                                                                                        "mc/corr/h6_corr_trueGeoAcc_pi0", "", "", "",
                                                                                                        "", "", "", "", ""},
                                                                                                       {"", "",
                                                                                                        "", "", "",
                                                                                                        "mc/corr/h6_mix_trueGeoAcc_pi0", "", "", "",
                                                                                                        "", "", "", "", ""}}}}};
  static constexpr const char* getMcTrueCorrHistPath(McCorrTruthLevel truthLevel, McCorrCorrelationType correlationType, McCorrAssocType assocType)
  {
    return McTrueCorrHistPaths[static_cast<int>(truthLevel)][static_cast<int>(correlationType)][static_cast<int>(assocType)];
  }

  // fill mc correlation histograms based on given associated mc particle
  template <McCorrTruthLevel truthLevel, McCorrCorrelationType correlationType>
  void fillMcCorrHists(auto const& mcCollision, auto const& triggerIn, auto const& associated, double const weight)
  {
    // check mc
    if constexpr (requires { triggerIn.template jetTrack_as<aod::JetTracksMCD>(); }) {
      if (!triggerIn.template jetTrack_as<aod::JetTracksMCD>().has_mcParticle()) {
        return;
      }
    }
    // trigger info separation
    auto const triggerParticle = [&]() {
      if constexpr (correlationType == McCorrCorrelationType::Correlation) {
        if constexpr (truthLevel == McCorrTruthLevel::TrueAssocEvRecoPtTrig) {
          static_assert(requires { triggerIn.template jetTrack_as<aod::JetTracksMCD>(); }, "triggerIn must be reco level for TrueAssocEvRecoPtTrig");
          return triggerIn.template jetTrack_as<aod::JetTracksMCD>().mcParticle();
        } else if constexpr (truthLevel == McCorrTruthLevel::TrueAssocEv || truthLevel == McCorrTruthLevel::True) {
          static_assert(requires { triggerIn.template jetMcParticle_as<aod::JetParticles>(); }, "triggerIn must be true level for TrueAssocEv and True");
          return triggerIn;
        }
      } else if constexpr (correlationType == McCorrCorrelationType::Mixing) {
        static_assert(std::same_as<std::remove_cvref_t<decltype(triggerIn)>, MixingTrigger>, "triggerIn must be MixingTrigger for mixing");
        return triggerIn;
      }
    }();
    auto const& triggerForPt = triggerIn;

    // exclude self correlation
    if constexpr (correlationType == McCorrCorrelationType::Correlation) {
      if constexpr (truthLevel == McCorrTruthLevel::TrueAssocEvRecoPtTrig) {
        if (triggerParticle.globalIndex() == associated.globalIndex()) {
          return;
        }
      } else if constexpr (truthLevel == McCorrTruthLevel::TrueAssocEv || truthLevel == McCorrTruthLevel::True) {
        if (triggerParticle.jetMcParticleId() == associated.globalIndex()) {
          return;
        }
      }
    }

    double const dPhi = getDeltaPhi(triggerParticle.phi(), associated.phi());
    double const dEta = triggerParticle.eta() - associated.eta();
    double const ptTrig = triggerForPt.pt();
    double const ptAssoc = associated.pt();
    double const posZ = mcCollision.posZ();
    double const mult = mcCollision.nChargedInEtaRange();

    if (std::abs(associated.eta()) > absEtaMax) {
      return;
    }

    // standard particles (marked physical primary)
    if (associated.isPhysicalPrimary()) {
      // charged primary ('hadron')
      if (checkChargedMc(associated)) {
        histos.fill(HIST(getMcTrueCorrHistPath(truthLevel, correlationType, McCorrAssocType::Hadron)),
                    dPhi, dEta, ptTrig, ptAssoc, posZ, mult, weight);
      }
      // pipm
      if (std::abs(associated.pdgCode()) == PDG_t::kPiPlus) {
        histos.fill(HIST(getMcTrueCorrHistPath(truthLevel, correlationType, McCorrAssocType::Pipm)),
                    dPhi, dEta, ptTrig, ptAssoc, posZ, mult, weight);
        return;
      }
      // photon
      if (associated.pdgCode() == PDG_t::kGamma) {
        histos.fill(HIST(getMcTrueCorrHistPath(truthLevel, correlationType, McCorrAssocType::Photon)),
                    dPhi, dEta, ptTrig, ptAssoc, posZ, mult, weight);

        // extra info for decay correlation only for total true level
        if constexpr (truthLevel != McCorrTruthLevel::True) {
          return;
        }

        // decay and direct
        int const statusDecayLow = 91;
        int const statusDecayUp = 99;
        int const statusDirect = 62;

        if (associated.getGenStatusCode() >= statusDecayLow && associated.getGenStatusCode() <= statusDecayUp) {
          histos.fill(HIST(getMcTrueCorrHistPath(truthLevel, correlationType, McCorrAssocType::PhotonDecay)),
                      dPhi, dEta, ptTrig, ptAssoc, posZ, mult,
                      weight);
          // decays from different mothers
          int const pdgMother = associated.template mothers_as<aod::JetParticles>().begin()->pdgCode();
          switch (pdgMother) {
            case PDG_t::kPi0:
              histos.fill(HIST(getMcTrueCorrHistPath(truthLevel, correlationType, McCorrAssocType::PhotonPi0)),
                          dPhi, dEta, ptTrig, ptAssoc, posZ, mult, weight);
              break;
            case constants::physics::Pdg::kEta:
              histos.fill(HIST(getMcTrueCorrHistPath(truthLevel, correlationType, McCorrAssocType::PhotonEta)),
                          dPhi, dEta, ptTrig, ptAssoc, posZ, mult, weight);
              break;
            case constants::physics::Pdg::kOmega:
              histos.fill(HIST(getMcTrueCorrHistPath(truthLevel, correlationType, McCorrAssocType::PhotonOmega)),
                          dPhi, dEta, ptTrig, ptAssoc, posZ, mult, weight);
              break;
            case constants::physics::Pdg::kEtaPrime:
              histos.fill(HIST(getMcTrueCorrHistPath(truthLevel, correlationType, McCorrAssocType::PhotonEtaPrime)),
                          dPhi, dEta, ptTrig, ptAssoc, posZ, mult, weight);
              break;
            default:
              histos.fill(HIST(getMcTrueCorrHistPath(truthLevel, correlationType, McCorrAssocType::PhotonOtherMother)),
                          dPhi, dEta, ptTrig, ptAssoc, posZ, mult, weight);
              break;
          }
        } else {
          histos.fill(HIST(getMcTrueCorrHistPath(truthLevel, correlationType, McCorrAssocType::PhotonDirect)),
                      dPhi, dEta, ptTrig, ptAssoc, posZ, mult,
                      weight);

          if (associated.getGenStatusCode() != statusDirect) {
            LOGF(info, "filled primary photon with status: %i", associated.getGenStatusCode());
          }
        }
        return;
      }
      return;
    }
    // decaying particles (not marked physical primary)
    if (!checkDecayPrimary(associated)) {
      return;
    }
    // pi0
    if (associated.pdgCode() == PDG_t::kPi0) {
      histos.fill(HIST(getMcTrueCorrHistPath(truthLevel, correlationType, McCorrAssocType::Pi0)),
                  dPhi, dEta, ptTrig, ptAssoc, posZ, mult, weight);
      // efficincy extra info
      if constexpr (truthLevel != McCorrTruthLevel::TrueAssocEvRecoPtTrig) {
        return;
      }
      // chosen decay
      if (!checkToGG(associated)) {
        return;
      }
      histos.fill(HIST(getMcTrueCorrHistPath(McCorrTruthLevel::TrueMeasDecay, correlationType, McCorrAssocType::Pi0)),
                  dPhi, dEta, ptTrig, ptAssoc, posZ, mult, weight);
      // daughters in acceptance
      if (!checkDaughtersInAcceptance(associated, {-absEtaMax, absEtaMax})) {
        return;
      }
      histos.fill(HIST(getMcTrueCorrHistPath(McCorrTruthLevel::TrueGeoAcc, correlationType, McCorrAssocType::Pi0)),
                  dPhi, dEta, ptTrig, ptAssoc, posZ, mult, weight);
      return;
    }
    // eta
    if (associated.pdgCode() == constants::physics::Pdg::kEta) {
      histos.fill(HIST(getMcTrueCorrHistPath(truthLevel, correlationType, McCorrAssocType::Eta)),
                  dPhi, dEta, ptTrig, ptAssoc, posZ, mult, weight);
      return;
    }
    // extra info for decay correlation only for total true level
    if constexpr (truthLevel != McCorrTruthLevel::True) {
      return;
    }
    // omega
    if (associated.pdgCode() == constants::physics::Pdg::kOmega) {
      histos.fill(HIST(getMcTrueCorrHistPath(truthLevel, correlationType, McCorrAssocType::Omega)),
                  dPhi, dEta, ptTrig, ptAssoc, posZ, mult, weight);
      return;
    }
    // etaPrime
    if (associated.pdgCode() == constants::physics::Pdg::kEtaPrime) {
      histos.fill(HIST(getMcTrueCorrHistPath(truthLevel, correlationType, McCorrAssocType::EtaPrime)),
                  dPhi, dEta, ptTrig, ptAssoc, posZ, mult, weight);
      return;
    }
  }

  void processMcTrueCorr(CorrMcCollisions const& mcCollisions, aod::TriggerParticles const& triggerParticles, aod::JetParticles const& mcParticles)
  {
    for (auto const& mcCollision : mcCollisions) {
      // event selection
      if (!totalEvSel(mcCollision)) {
        continue;
      }

      // group collision
      auto const triggerParticlesThisEvent = triggerParticles.sliceBy(perColTriggerParticles, mcCollision.globalIndex());
      auto const mcParticlesThisEvent = mcParticles.sliceBy(perColMcParticles, mcCollision.globalIndex());

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        fillMcCorrHists<McCorrTruthLevel::True, McCorrCorrelationType::Correlation>(collision, trigger, associated, 1);
      };
      corrProcessCorrelation(mcCollision, triggerParticlesThisEvent, mcParticlesThisEvent, funcCorrelation);

      // select mixing events
      if (doTrigEvMixing && !mcCollision.trigEv()) {
        continue;
      }

      auto const funcMixing = [this](auto const& collision, auto const& trigger, auto const& associated, auto const perTriggerWeight) {
        fillMcCorrHists<McCorrTruthLevel::True, McCorrCorrelationType::Mixing>(collision, trigger, associated, perTriggerWeight);
      };
      corrProcessMixing<TriggerMemoryMode::True>(mcCollision, mcParticlesThisEvent, funcMixing, nMixingAt0McTrue);
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcTrueCorr, "process mc-true (all collisions) correlation for multiple associated particles", false);

  void processMcTrueAssocEvCorr(CorrMcDCollisions const& collisions, CorrMcCollisions const&, aod::TriggerParticles const& triggerParticles, aod::JetParticles const& mcParticles)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      auto const& mcCollision = collision.mcCollision_as<CorrMcCollisions>();

      // group collision
      auto const triggerParticlesThisEvent = triggerParticles.sliceBy(perColTriggerParticles, collision.mcCollisionId());
      auto const mcParticlesThisEvent = mcParticles.sliceBy(perColMcParticles, collision.mcCollisionId());

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        fillMcCorrHists<McCorrTruthLevel::TrueAssocEv, McCorrCorrelationType::Correlation>(collision, trigger, associated, 1);
      };
      corrProcessCorrelation(mcCollision, triggerParticlesThisEvent, mcParticlesThisEvent, funcCorrelation);

      // select mixing events
      if (doTrigEvMixing && !mcCollision.trigEv()) {
        continue;
      }

      auto const funcMixing = [this](auto const& collision, auto const& trigger, auto const& associated, auto const perTriggerWeight) {
        fillMcCorrHists<McCorrTruthLevel::TrueAssocEv, McCorrCorrelationType::Mixing>(collision, trigger, associated, perTriggerWeight);
      };
      corrProcessMixing<TriggerMemoryMode::TrueAssocEv>(mcCollision, mcParticlesThisEvent, funcMixing, nMixingAt0McTrue);
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcTrueAssocEvCorr, "process mc-true (reco collisions) correlation for multiple associated particles", false);

  void processMcTrueAssocEvRecoPtTrigCorr(CorrMcDCollisions const& collisions, aod::Triggers const& triggers, aod::JetTracksMCD const&, CorrMcCollisions const&, aod::JetParticles const& mcParticles)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      auto const& mcCollision = collision.mcCollision_as<CorrMcCollisions>();

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const mcParticlesThisEvent = mcParticles.sliceBy(perColMcParticles, collision.mcCollisionId());

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        fillMcCorrHists<McCorrTruthLevel::TrueAssocEvRecoPtTrig, McCorrCorrelationType::Correlation>(collision, trigger, associated, 1);
      };
      corrProcessCorrelation(mcCollision, triggersThisEvent, mcParticlesThisEvent, funcCorrelation);

      // select mixing events
      if (doTrigEvMixing && !mcCollision.trigEv()) {
        continue;
      }

      auto const funcMixing = [this](auto const& collision, auto const& trigger, auto const& associated, auto const perTriggerWeight) {
        fillMcCorrHists<McCorrTruthLevel::TrueAssocEvRecoPtTrig, McCorrCorrelationType::Mixing>(collision, trigger, associated, perTriggerWeight);
      };
      corrProcessMixing<TriggerMemoryMode::TrueAssocEvRecoPtTrig>(mcCollision, mcParticlesThisEvent, funcMixing, nMixingAt0McTrue);
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcTrueAssocEvRecoPtTrigCorr, "process mc-true (reco collisions with reco ptTrig) correlation for multiple associated particles", false);

  void processMcRecoCorrHadron(CorrMcDCollisions const& collisions, aod::Triggers const& triggers, aod::Hadrons const& hadrons,
                               aod::JetTracksMCD const&, aod::JetParticles const&)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const hadronsThisEvent = hadrons.sliceBy(perColHadrons, collision.globalIndex());

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // exclude self correlation
        if (trigger.jetTrackId() == associated.jetTrackId()) {
          return;
        }

        // check mc
        if (!associated.template jetTrack_as<aod::JetTracksMCD>().has_mcParticle()) {
          return;
        }
        auto const& associatedMcParticle = associated.template jetTrack_as<aod::JetTracksMCD>().mcParticle();
        // collision association
        if (associatedMcParticle.mcCollisionId() != collision.mcCollisionId()) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_corr_recoAssocEv_hadron"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // purity
        if (!checkChargedMc(associatedMcParticle) || !associatedMcParticle.isPhysicalPrimary()) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_corr_recoPure_hadron"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // true pt
        histos.fill(HIST("mc/corr/h6_corr_recoPureTruePtAssoc_hadron"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associatedMcParticle.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
      };
      corrProcessCorrelation(collision, triggersThisEvent, hadronsThisEvent, funcCorrelation);

      // select mixing events
      if (doTrigEvMixing && !collision.trigEv()) {
        continue;
      }

      auto const funcMixing = [this](auto const& collision, auto const& trigger, auto const& associated, auto const perTriggerWeight) {
        // check mc
        if (!associated.template jetTrack_as<aod::JetTracksMCD>().has_mcParticle()) {
          return;
        }
        auto const& associatedMcParticle = associated.template jetTrack_as<aod::JetTracksMCD>().mcParticle();
        // collision association
        if (associatedMcParticle.mcCollisionId() != collision.mcCollisionId()) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_mix_recoAssocEv_hadron"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // purity
        if (!checkChargedMc(associatedMcParticle) || !associatedMcParticle.isPhysicalPrimary()) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_mix_recoPure_hadron"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // true pt
        histos.fill(HIST("mc/corr/h6_mix_recoPureTruePtAssoc_hadron"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associatedMcParticle.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
      };
      corrProcessMixing<TriggerMemoryMode::Reco>(collision, hadronsThisEvent, funcMixing, nMixingAt0Hadron);
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcRecoCorrHadron, "process correlation for associated hardons with additional mc information", false);

  void processMcRecoCorrPipm(CorrMcDCollisions const& collisions, aod::Triggers const& triggers, aod::Pipms const& pipms,
                             aod::JetTracksMCD const&, aod::JetParticles const&)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const pipmsThisEvent = pipms.sliceBy(perColPipms, collision.globalIndex());

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // exclude self correlation
        if (trigger.jetTrackId() == associated.jetTrackId()) {
          return;
        }

        // check mc
        if (!associated.template jetTrack_as<aod::JetTracksMCD>().has_mcParticle()) {
          return;
        }
        auto const& associatedMcParticle = associated.template jetTrack_as<aod::JetTracksMCD>().mcParticle();
        // collision association
        if (associatedMcParticle.mcCollisionId() != collision.mcCollisionId()) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_corr_recoAssocEv_pipm"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // purity
        if (std::abs(associatedMcParticle.pdgCode()) != PDG_t::kPiPlus || !associatedMcParticle.isPhysicalPrimary()) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_corr_recoPure_pipm"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // true pt
        histos.fill(HIST("mc/corr/h6_corr_recoPureTruePtAssoc_pipm"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associatedMcParticle.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
      };
      corrProcessCorrelation(collision, triggersThisEvent, pipmsThisEvent, funcCorrelation);

      // select mixing events
      if (doTrigEvMixing && !collision.trigEv()) {
        continue;
      }

      auto const funcMixing = [this](auto const& collision, auto const& trigger, auto const& associated, auto const perTriggerWeight) {
        // check mc
        if (!associated.template jetTrack_as<aod::JetTracksMCD>().has_mcParticle()) {
          return;
        }
        auto const& associatedMcParticle = associated.template jetTrack_as<aod::JetTracksMCD>().mcParticle();
        // collision association
        if (associatedMcParticle.mcCollisionId() != collision.mcCollisionId()) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_mix_recoAssocEv_pipm"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // purity
        if (std::abs(associatedMcParticle.pdgCode()) != PDG_t::kPiPlus || !associatedMcParticle.isPhysicalPrimary()) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_mix_recoPure_pipm"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // true pt
        histos.fill(HIST("mc/corr/h6_mix_recoPureTruePtAssoc_pipm"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associatedMcParticle.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
      };
      corrProcessMixing<TriggerMemoryMode::Reco>(collision, pipmsThisEvent, funcMixing, nMixingAt0Pipm);
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcRecoCorrPipm, "process correlation for associated pipms with additional mc information", false);

  void processMcRecoCorrPhotonPCM(CorrMcDCollisions const& collisions, aod::Triggers const& triggers, aod::PhotonPCMs const& photonPCMs,
                                  aod::JetTracksMCD const&, aod::JetParticles const&)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const photonPCMsThisEvent = photonPCMs.sliceBy(perColPhotonPCMs, collision.globalIndex());

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // exclude self correlation
        if (trigger.jetTrackId() == associated.posJetTrackId() || trigger.jetTrackId() == associated.negJetTrackId()) {
          return;
        }

        // check mc
        auto const& posTrack = associated.template posJetTrack_as<aod::JetTracksMCD>();
        auto const& negTrack = associated.template negJetTrack_as<aod::JetTracksMCD>();
        if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle()) {
          return;
        }
        // collision association
        if (posTrack.mcParticle().mcCollisionId() != collision.mcCollisionId() || negTrack.mcParticle().mcCollisionId() != collision.mcCollisionId()) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_corr_recoAssocEv_photonPCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // purity
        auto const& photons = posTrack.mcParticle().template mothers_as<aod::JetParticles>();
        if (!isConversionPhoton(posTrack, negTrack) || !photons.begin()->isPhysicalPrimary()) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_corr_recoPure_photonPCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // true pt
        histos.fill(HIST("mc/corr/h6_corr_recoPureTruePtAssoc_photonPCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), photons.begin()->pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
      };
      corrProcessCorrelation(collision, triggersThisEvent, photonPCMsThisEvent, funcCorrelation);

      // select mixing events
      if (doTrigEvMixing && !collision.trigEv()) {
        continue;
      }

      auto const funcMixing = [this](auto const& collision, auto const& trigger, auto const& associated, auto const perTriggerWeight) {
        // check mc
        auto const& posTrack = associated.template posJetTrack_as<aod::JetTracksMCD>();
        auto const& negTrack = associated.template negJetTrack_as<aod::JetTracksMCD>();
        if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle()) {
          return;
        }
        // collision association
        if (posTrack.mcParticle().mcCollisionId() != collision.mcCollisionId() || negTrack.mcParticle().mcCollisionId() != collision.mcCollisionId()) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_mix_recoAssocEv_photonPCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // purity
        auto const& photons = posTrack.mcParticle().template mothers_as<aod::JetParticles>();
        if (!isConversionPhoton(posTrack, negTrack) || !photons.begin()->isPhysicalPrimary()) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_mix_recoPure_photonPCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // true pt
        histos.fill(HIST("mc/corr/h6_mix_recoPureTruePtAssoc_photonPCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), photons.begin()->pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
      };
      corrProcessMixing<TriggerMemoryMode::Reco>(collision, photonPCMsThisEvent, funcMixing, nMixingAt0PhotonPCM);
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcRecoCorrPhotonPCM, "process correlation for associated photonPCMs with additional mc information", false);

  void processMcRecoCorrPi0PCM(CorrMcDCollisions const& collisions, aod::Triggers const& triggers, aod::PhotonPCMPairs const& photonPCMPairs,
                               aod::JetTracksMCD const&, aod::JetParticles const&)
  {
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());
      auto const photonPCMPairsThisEvent = photonPCMPairs.sliceBy(perColPhotonPCMPairs, collision.globalIndex());

      auto const funcCorrelation = [this](auto const& collision, auto const& trigger, auto const& associated) {
        // exclude self correlation
        if (trigger.jetTrackId() == associated.posJetTrack1Id() || trigger.jetTrackId() == associated.negJetTrack1Id() ||
            trigger.jetTrackId() == associated.posJetTrack2Id() || trigger.jetTrackId() == associated.negJetTrack2Id()) {
          return;
        }

        // check mc
        auto const& posTrack1 = associated.template posJetTrack1_as<aod::JetTracksMCD>();
        auto const& negTrack1 = associated.template negJetTrack1_as<aod::JetTracksMCD>();
        auto const& posTrack2 = associated.template posJetTrack2_as<aod::JetTracksMCD>();
        auto const& negTrack2 = associated.template negJetTrack2_as<aod::JetTracksMCD>();
        if (!posTrack1.has_mcParticle() || !negTrack1.has_mcParticle() || !posTrack2.has_mcParticle() || !negTrack2.has_mcParticle()) {
          return;
        }
        // pseudo yield
        if (!isGGFromDoubleConversion(posTrack1, negTrack1, posTrack2, negTrack2, PDG_t::kPi0)) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_corr_pseudoReco_pi0PCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // collision association
        if (posTrack1.mcParticle().mcCollisionId() != collision.mcCollisionId()) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_corr_pseudoRecoAssocEv_pi0PCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // purity (just secondaries)
        auto const& photons1 = posTrack1.mcParticle().template mothers_as<aod::JetParticles>();
        auto const& mothersOfPhoton = photons1.begin()->template mothers_as<aod::JetParticles>();
        if (!checkDecayPrimary(*(mothersOfPhoton.begin()))) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_corr_pseudoRecoPure_pi0PCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // true pt
        histos.fill(HIST("mc/corr/h6_corr_pseudoRecoPureTruePtAssoc_pi0PCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), mothersOfPhoton.begin()->pt(), collision.posZ(), collision.nGlobalTracks(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
      };
      corrProcessCorrelation(collision, triggersThisEvent, photonPCMPairsThisEvent, funcCorrelation);

      // select mixing events
      if (doTrigEvMixing && !collision.trigEv()) {
        continue;
      }

      auto const funcMixing = [this](auto const& collision, auto const& trigger, auto const& associated, auto const perTriggerWeight) {
        // check mc
        auto const& posTrack1 = associated.template posJetTrack1_as<aod::JetTracksMCD>();
        auto const& negTrack1 = associated.template negJetTrack1_as<aod::JetTracksMCD>();
        auto const& posTrack2 = associated.template posJetTrack2_as<aod::JetTracksMCD>();
        auto const& negTrack2 = associated.template negJetTrack2_as<aod::JetTracksMCD>();
        if (!posTrack1.has_mcParticle() || !negTrack1.has_mcParticle() || !posTrack2.has_mcParticle() || !negTrack2.has_mcParticle()) {
          return;
        }
        // pseudo yield
        if (!isGGFromDoubleConversion(posTrack1, negTrack1, posTrack2, negTrack2, PDG_t::kPi0)) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_mix_pseudoReco_pi0PCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // collision association
        if (posTrack1.mcParticle().mcCollisionId() != collision.mcCollisionId()) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_mix_pseudoRecoAssocEv_pi0PCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // purity (just secondaries)
        auto const& photons1 = posTrack1.mcParticle().template mothers_as<aod::JetParticles>();
        auto const& mothersOfPhoton = photons1.begin()->template mothers_as<aod::JetParticles>();
        if (!checkDecayPrimary(*(mothersOfPhoton.begin()))) {
          return;
        }
        histos.fill(HIST("mc/corr/h6_mix_pseudoRecoPure_pi0PCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), associated.pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
        // true pt
        histos.fill(HIST("mc/corr/h6_mix_pseudoRecoPureTruePtAssoc_pi0PCM"),
                    getDeltaPhi(trigger.phi(), associated.phi()),
                    trigger.eta() - associated.eta(),
                    trigger.pt(), mothersOfPhoton.begin()->pt(), collision.posZ(), collision.nGlobalTracks(),
                    perTriggerWeight * getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));
      };
      corrProcessMixing<TriggerMemoryMode::Reco>(collision, photonPCMPairsThisEvent, funcMixing, nMixingAt0PhotonPCM);
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcRecoCorrPi0PCM, "process correlation for associated pi0PCMs with additional mc information", false);

  // other ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  enum class McCorrectionObservable : int { PhiEta = 0,
                                            ZPvMult = 1 };
  enum class McResolvedTruthLevel : int { Reco = 0,
                                          RecoAssocEv = 1,
                                          RecoPure = 2,
                                          TrueAssocEvRecoPtTrig = 3,
                                          TrueAssocEv = 4,
                                          True = 5 };
  enum class CorrAssocType : int { Hadron = 0,
                                   Pipm = 1,
                                   PhotonPCM = 2 };
  static constexpr std::array<std::array<std::array<char const*, 3>, 6>, 2> McResolvedCorrectionHistPaths = {
    {{{{"mc/correction/resol/h4_ptTrigPtAssocPhiEta_reco_hadron",
        "mc/correction/resol/h4_ptTrigPtAssocPhiEta_reco_pipm",
        "mc/correction/resol/h4_ptTrigPtAssocPhiEta_reco_photonPCM"},
       {"mc/correction/resol/h4_ptTrigPtAssocPhiEta_recoAssocEv_hadron",
        "mc/correction/resol/h4_ptTrigPtAssocPhiEta_recoAssocEv_pipm",
        "mc/correction/resol/h4_ptTrigPtAssocPhiEta_recoAssocEv_photonPCM"},
       {"mc/correction/resol/h4_ptTrigPtAssocPhiEta_recoPure_hadron",
        "mc/correction/resol/h4_ptTrigPtAssocPhiEta_recoPure_pipm",
        "mc/correction/resol/h4_ptTrigPtAssocPhiEta_recoPure_photonPCM"},
       {"mc/correction/resol/h4_ptTrigPtAssocPhiEta_trueAssocEvRecoPtTrig_hadron",
        "mc/correction/resol/h4_ptTrigPtAssocPhiEta_trueAssocEvRecoPtTrig_pipm",
        "mc/correction/resol/h4_ptTrigPtAssocPhiEta_trueAssocEvRecoPtTrig_photon"},
       {"mc/correction/resol/h4_ptTrigPtAssocPhiEta_trueAssocEv_hadron",
        "mc/correction/resol/h4_ptTrigPtAssocPhiEta_trueAssocEv_pipm",
        "mc/correction/resol/h4_ptTrigPtAssocPhiEta_trueAssocEv_photon"},
       {"mc/correction/resol/h4_ptTrigPtAssocPhiEta_true_hadron",
        "mc/correction/resol/h4_ptTrigPtAssocPhiEta_true_pipm",
        "mc/correction/resol/h4_ptTrigPtAssocPhiEta_true_photon"}}},
     {{{"mc/correction/resol/h4_ptTrigPtAssocZPvMult_reco_hadron",
        "mc/correction/resol/h4_ptTrigPtAssocZPvMult_reco_pipm",
        "mc/correction/resol/h4_ptTrigPtAssocZPvMult_reco_photonPCM"},
       {"mc/correction/resol/h4_ptTrigPtAssocZPvMult_recoAssocEv_hadron",
        "mc/correction/resol/h4_ptTrigPtAssocZPvMult_recoAssocEv_pipm",
        "mc/correction/resol/h4_ptTrigPtAssocZPvMult_recoAssocEv_photonPCM"},
       {"mc/correction/resol/h4_ptTrigPtAssocZPvMult_recoPure_hadron",
        "mc/correction/resol/h4_ptTrigPtAssocZPvMult_recoPure_pipm",
        "mc/correction/resol/h4_ptTrigPtAssocZPvMult_recoPure_photonPCM"},
       {"mc/correction/resol/h4_ptTrigPtAssocZPvMult_trueAssocEvRecoPtTrig_hadron",
        "mc/correction/resol/h4_ptTrigPtAssocZPvMult_trueAssocEvRecoPtTrig_pipm",
        "mc/correction/resol/h4_ptTrigPtAssocZPvMult_trueAssocEvRecoPtTrig_photon"},
       {"mc/correction/resol/h4_ptTrigPtAssocZPvMult_trueAssocEv_hadron",
        "mc/correction/resol/h4_ptTrigPtAssocZPvMult_trueAssocEv_pipm",
        "mc/correction/resol/h4_ptTrigPtAssocZPvMult_trueAssocEv_photon"},
       {"mc/correction/resol/h4_ptTrigPtAssocZPvMult_true_hadron",
        "mc/correction/resol/h4_ptTrigPtAssocZPvMult_true_pipm",
        "mc/correction/resol/h4_ptTrigPtAssocZPvMult_true_photon"}}}}};
  static constexpr char const* getMcResolvedCorrectionHistPath(McCorrectionObservable observable, McResolvedTruthLevel truthLevel, CorrAssocType assocType)
  {
    return McResolvedCorrectionHistPaths[static_cast<int>(observable)][static_cast<int>(truthLevel)][static_cast<int>(assocType)];
  }
  template <McResolvedTruthLevel truthLevel, CorrAssocType assocType>
  void fillResolvedCorrectionHists(auto const& collision, auto const& associated)
  {
    if constexpr (truthLevel == McResolvedTruthLevel::Reco || truthLevel == McResolvedTruthLevel::RecoAssocEv || truthLevel == McResolvedTruthLevel::RecoPure) {
      histos.fill(HIST(getMcResolvedCorrectionHistPath(McCorrectionObservable::PhiEta, truthLevel, assocType)),
                  collision.ptMax(), associated.pt(), associated.phi(), associated.eta());
      histos.fill(HIST(getMcResolvedCorrectionHistPath(McCorrectionObservable::ZPvMult, truthLevel, assocType)),
                  collision.ptMax(), associated.pt(), collision.posZ(), collision.nGlobalTracks());
    }
    if constexpr (truthLevel == McResolvedTruthLevel::TrueAssocEvRecoPtTrig) {
      histos.fill(HIST(getMcResolvedCorrectionHistPath(McCorrectionObservable::PhiEta, truthLevel, assocType)),
                  collision.ptMax(), associated.pt(), associated.phi(), associated.eta());
      histos.fill(HIST(getMcResolvedCorrectionHistPath(McCorrectionObservable::ZPvMult, truthLevel, assocType)),
                  collision.ptMax(), associated.pt(), collision.template mcCollision_as<CorrMcCollisions>().posZ(), collision.template mcCollision_as<CorrMcCollisions>().nChargedInEtaRange());
    }
    if constexpr (truthLevel == McResolvedTruthLevel::True || truthLevel == McResolvedTruthLevel::TrueAssocEv) {
      histos.fill(HIST(getMcResolvedCorrectionHistPath(McCorrectionObservable::PhiEta, truthLevel, assocType)),
                  collision.ptMax(), associated.pt(), associated.phi(), associated.eta());
      histos.fill(HIST(getMcResolvedCorrectionHistPath(McCorrectionObservable::ZPvMult, truthLevel, assocType)),
                  collision.ptMax(), associated.pt(), collision.posZ(), collision.nChargedInEtaRange());
    }
  }

  enum class PhotonImpurity { WrongEv,
                              Misid,
                              K0Short,
                              K0Long,
                              Lambda0,
                              Other };
  static constexpr double binValuePi0Impurities(PhotonImpurity const impurity)
  {
    return static_cast<int>(impurity) + 0.5;
  }

  void processMcRecoCorrection(CorrMcDCollision const& collision, aod::JetTracksMCD const&,
                               aod::Hadrons const& hadrons, aod::Pipms const& pipms, aod::PhotonPCMs const& photonPCMs, aod::PhotonPCMPairs const& photonPCMPairs,
                               CorrMcCollisions const&, aod::JetParticles const& mcParticles)
  {
    // event selection
    if (!totalEvSel(collision)) {
      return;
    }

    auto const& mcCollision = collision.mcCollision_as<CorrMcCollisions>();

    // group collision
    auto const mcParticlesThisEvent = mcParticles.sliceBy(perColMcParticles, collision.mcCollisionId());

    // hadrons
    for (auto const& hadron : hadrons) {
      // reconstructed
      fillResolvedCorrectionHists<McResolvedTruthLevel::Reco, CorrAssocType::Hadron>(collision, hadron);
      // check mc
      if (!hadron.jetTrack_as<aod::JetTracksMCD>().has_mcParticle()) {
        continue;
      }
      auto const& hadronParticle = hadron.jetTrack_as<aod::JetTracksMCD>().mcParticle();
      // collision association
      if (hadronParticle.mcCollisionId() != collision.mcCollisionId()) {
        continue;
      }
      fillResolvedCorrectionHists<McResolvedTruthLevel::RecoAssocEv, CorrAssocType::Hadron>(collision, hadron);
      // purity
      if (!checkChargedMc(hadronParticle) || !hadronParticle.isPhysicalPrimary()) {
        continue;
      }
      fillResolvedCorrectionHists<McResolvedTruthLevel::RecoPure, CorrAssocType::Hadron>(collision, hadron);
    }

    // pipm
    for (auto const& pipm : pipms) {
      // reconstructed
      fillResolvedCorrectionHists<McResolvedTruthLevel::Reco, CorrAssocType::Pipm>(collision, pipm);
      // check mc
      if (!pipm.jetTrack_as<aod::JetTracksMCD>().has_mcParticle()) {
        continue;
      }
      auto const& pipmParticle = pipm.jetTrack_as<aod::JetTracksMCD>().mcParticle();
      // collision association
      if (pipmParticle.mcCollisionId() != collision.mcCollisionId()) {
        continue;
      }
      fillResolvedCorrectionHists<McResolvedTruthLevel::RecoAssocEv, CorrAssocType::Pipm>(collision, pipm);
      // purity
      if (std::abs(pipmParticle.pdgCode()) != PDG_t::kPiPlus || !pipmParticle.isPhysicalPrimary()) {
        continue;
      }
      fillResolvedCorrectionHists<McResolvedTruthLevel::RecoPure, CorrAssocType::Pipm>(collision, pipm);
    }

    // photonPCM
    for (auto const& photonPCM : photonPCMs) {
      // reconstructed
      fillResolvedCorrectionHists<McResolvedTruthLevel::Reco, CorrAssocType::PhotonPCM>(collision, photonPCM);
      // check mc
      auto const& posTrack = photonPCM.posJetTrack_as<aod::JetTracksMCD>();
      auto const& negTrack = photonPCM.negJetTrack_as<aod::JetTracksMCD>();
      if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle()) {
        continue;
      }
      // collision association
      if (posTrack.mcParticle().mcCollisionId() != collision.mcCollisionId() || negTrack.mcParticle().mcCollisionId() != collision.mcCollisionId()) {
        histos.fill(HIST("mc/correction/h3_ptTrigPtAssocCategory_recoImpurities_photonPCM"),
                    collision.ptMax(), photonPCM.pt(), binValuePi0Impurities(PhotonImpurity::WrongEv));
        continue;
      }
      fillResolvedCorrectionHists<McResolvedTruthLevel::RecoAssocEv, CorrAssocType::PhotonPCM>(collision, photonPCM);
      // purity
      if (!isConversionPhoton(posTrack, negTrack)) {
        histos.fill(HIST("mc/correction/h3_ptTrigPtAssocCategory_recoImpurities_photonPCM"),
                    collision.ptMax(), photonPCM.pt(), binValuePi0Impurities(PhotonImpurity::Misid));
        continue;
      }
      auto const& photons = posTrack.mcParticle().mothers_as<aod::JetParticles>();
      if (!photons.begin()->isPhysicalPrimary()) {
        if (checkForMother(*(photons.begin()), PDG_t::kK0Short, true)) {
          histos.fill(HIST("mc/correction/h3_ptTrigPtAssocCategory_recoImpurities_photonPCM"),
                      collision.ptMax(), photonPCM.pt(), binValuePi0Impurities(PhotonImpurity::K0Short));
          continue;
        }
        if (checkForMother(*(photons.begin()), PDG_t::kK0Long, true)) {
          histos.fill(HIST("mc/correction/h3_ptTrigPtAssocCategory_recoImpurities_photonPCM"),
                      collision.ptMax(), photonPCM.pt(), binValuePi0Impurities(PhotonImpurity::K0Long));
          continue;
        }
        if (checkForMother(*(photons.begin()), PDG_t::kLambda0, true)) {
          histos.fill(HIST("mc/correction/h3_ptTrigPtAssocCategory_recoImpurities_photonPCM"),
                      collision.ptMax(), photonPCM.pt(), binValuePi0Impurities(PhotonImpurity::Lambda0));
          continue;
        }
        histos.fill(HIST("mc/correction/h3_ptTrigPtAssocCategory_recoImpurities_photonPCM"),
                    collision.ptMax(), photonPCM.pt(), binValuePi0Impurities(PhotonImpurity::Other));
        continue;
      }
      fillResolvedCorrectionHists<McResolvedTruthLevel::RecoPure, CorrAssocType::PhotonPCM>(collision, photonPCM);
    }

    // h0 PCM
    for (auto const& photonPCMPair : photonPCMPairs) {
      // check mc
      auto const& posTrack1 = photonPCMPair.posJetTrack1_as<aod::JetTracksMCD>();
      auto const& negTrack1 = photonPCMPair.negJetTrack1_as<aod::JetTracksMCD>();
      auto const& posTrack2 = photonPCMPair.posJetTrack2_as<aod::JetTracksMCD>();
      auto const& negTrack2 = photonPCMPair.negJetTrack2_as<aod::JetTracksMCD>();
      if (!posTrack1.has_mcParticle() || !negTrack1.has_mcParticle() || !posTrack2.has_mcParticle() || !negTrack2.has_mcParticle()) {
        continue;
      }

      // pi0PCM
      if (isGGFromDoubleConversion(posTrack1, negTrack1, posTrack2, negTrack2, PDG_t::kPi0)) {
        // pseudo reconstructed
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_pseudoReco_pi0PCM"), collision.ptMax(), photonPCMPair.pt());
        // collision association
        if (posTrack1.mcParticle().mcCollisionId() != collision.mcCollisionId()) {
          histos.fill(HIST("mc/correction/h3_ptTrigPtAssocCategory_pseudoRecoImpurities_pi0PCM"),
                      collision.ptMax(), photonPCMPair.pt(), binValuePi0Impurities(PhotonImpurity::WrongEv));
          continue;
        }
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_pseudoRecoAssocEv_pi0PCM"), collision.ptMax(), photonPCMPair.pt());
        // purity
        // note: do not dereference 'mothers' and store element in 'const&' (leads to seg fault)
        auto const& photons1 = posTrack1.mcParticle().mothers_as<aod::JetParticles>();
        auto const& mothersOfPhoton = photons1.begin()->mothers_as<aod::JetParticles>();
        auto const& grandmothersOfPhoton = mothersOfPhoton.begin()->mothers_as<aod::JetParticles>();
        if (!checkDecayPrimary(*(mothersOfPhoton.begin()))) {
          if (std::abs(grandmothersOfPhoton.begin()->pdgCode()) == PDG_t::kK0Short) {
            histos.fill(HIST("mc/correction/h3_ptTrigPtAssocCategory_pseudoRecoImpurities_pi0PCM"),
                        collision.ptMax(), photonPCMPair.pt(), binValuePi0Impurities(PhotonImpurity::K0Short));
            continue;
          }
          if (std::abs(grandmothersOfPhoton.begin()->pdgCode()) == PDG_t::kK0Long) {
            histos.fill(HIST("mc/correction/h3_ptTrigPtAssocCategory_pseudoRecoImpurities_pi0PCM"),
                        collision.ptMax(), photonPCMPair.pt(), binValuePi0Impurities(PhotonImpurity::K0Long));
            continue;
          }
          if (std::abs(grandmothersOfPhoton.begin()->pdgCode()) == PDG_t::kLambda0) {
            histos.fill(HIST("mc/correction/h3_ptTrigPtAssocCategory_pseudoRecoImpurities_pi0PCM"),
                        collision.ptMax(), photonPCMPair.pt(), binValuePi0Impurities(PhotonImpurity::Lambda0));
            continue;
          }
          histos.fill(HIST("mc/correction/h3_ptTrigPtAssocCategory_pseudoRecoImpurities_pi0PCM"),
                      collision.ptMax(), photonPCMPair.pt(), binValuePi0Impurities(PhotonImpurity::Other));
          continue;
        }
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_pseudoRecoPure_pi0PCM"),
                    collision.ptMax(), photonPCMPair.pt());
        continue;
      }
      // etaPCM
      if (isGGFromDoubleConversion(posTrack1, negTrack1, posTrack2, negTrack2, constants::physics::Pdg::kEta)) {
        // pseudo reconstructed
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_pseudoReco_etaPCM"), collision.ptMax(), photonPCMPair.pt());
        // collision association
        if (posTrack1.mcParticle().mcCollisionId() != collision.mcCollisionId()) {
          continue;
        }
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_pseudoRecoAssocEv_etaPCM"), collision.ptMax(), photonPCMPair.pt());
        // purity
        auto const& photons1 = posTrack1.mcParticle().mothers_as<aod::JetParticles>();
        auto const& mothersOfPhoton = photons1.begin()->mothers_as<aod::JetParticles>();
        if (!checkDecayPrimary(*(mothersOfPhoton.begin()))) {
          continue;
        }
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_pseudoRecoPure_etaPCM"), collision.ptMax(), photonPCMPair.pt());
        continue;
      }
    }

    // mcParticle loop
    for (auto const& mcParticle : mcParticlesThisEvent) {
      if (std::abs(mcParticle.eta()) > absEtaMax) {
        continue;
      }

      // standard particles (marked physical primary)
      if (mcParticle.isPhysicalPrimary()) {
        // charged
        // hadrons
        if (checkChargedMc(mcParticle)) {
          fillResolvedCorrectionHists<McResolvedTruthLevel::TrueAssocEvRecoPtTrig, CorrAssocType::Hadron>(collision, mcParticle);
          fillResolvedCorrectionHists<McResolvedTruthLevel::TrueAssocEv, CorrAssocType::Hadron>(mcCollision, mcParticle);
        }
        // pipm
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus) {
          fillResolvedCorrectionHists<McResolvedTruthLevel::TrueAssocEvRecoPtTrig, CorrAssocType::Pipm>(collision, mcParticle);
          fillResolvedCorrectionHists<McResolvedTruthLevel::TrueAssocEv, CorrAssocType::Pipm>(mcCollision, mcParticle);
          continue;
        }
        // non charged
        // photons
        if (mcParticle.pdgCode() == PDG_t::kGamma) {
          fillResolvedCorrectionHists<McResolvedTruthLevel::TrueAssocEvRecoPtTrig, CorrAssocType::PhotonPCM>(collision, mcParticle);
          fillResolvedCorrectionHists<McResolvedTruthLevel::TrueAssocEv, CorrAssocType::PhotonPCM>(mcCollision, mcParticle);
          continue;
        }
        continue;
      }

      // decaying particles (not marked physical primary)
      if (!checkDecayPrimary(mcParticle)) {
        continue;
      }
      // pi0
      if (mcParticle.pdgCode() == PDG_t::kPi0) {
        // true
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_trueAssocEv_pi0"), mcCollision.ptMax(), mcParticle.pt());
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_trueAssocEvRecoPtTrig_pi0"), collision.ptMax(), mcParticle.pt());
        // chosen decay
        if (!checkToGG(mcParticle)) {
          continue;
        }
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_trueMeasDecay_pi0"), collision.ptMax(), mcParticle.pt());
        // daughters in acceptance
        if (!checkDaughtersInAcceptance(mcParticle, {-absEtaMax, absEtaMax})) {
          continue;
        }
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_trueGeoAcc_pi0"), collision.ptMax(), mcParticle.pt());
        continue;
      }
      // eta
      if (mcParticle.pdgCode() == constants::physics::Pdg::kEta) {
        // true
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_trueAssocEv_eta"), mcCollision.ptMax(), mcParticle.pt());
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_trueAssocEvRecoPtTrig_eta"), collision.ptMax(), mcParticle.pt());
        // chosen decay
        if (!checkToGG(mcParticle)) {
          continue;
        }
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_trueMeasDecay_eta"), collision.ptMax(), mcParticle.pt());
        // daughters in acceptance
        if (!checkDaughtersInAcceptance(mcParticle, {-absEtaMax, absEtaMax})) {
          continue;
        }
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_trueGeoAcc_eta"), collision.ptMax(), mcParticle.pt());
        continue;
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcRecoCorrection, "process mc reconstruction to calculate corrections", false);

  void processMcTrueCorrection(CorrMcCollision const& mcCollision, aod::JetParticles const& mcParticles)
  {
    // event selection
    if (!totalEvSel(mcCollision)) {
      return;
    }

    for (auto const& mcParticle : mcParticles) {
      if (std::abs(mcParticle.eta()) > absEtaMax) {
        continue;
      }

      // standard particles (marked physical primary)
      if (mcParticle.isPhysicalPrimary()) {
        // charged
        // hadrons
        if (checkChargedMc(mcParticle)) {
          fillResolvedCorrectionHists<McResolvedTruthLevel::True, CorrAssocType::Hadron>(mcCollision, mcParticle);
        }
        // pipm
        if (std::abs(mcParticle.pdgCode()) == PDG_t::kPiPlus) {
          fillResolvedCorrectionHists<McResolvedTruthLevel::True, CorrAssocType::Pipm>(mcCollision, mcParticle);
          continue;
        }
        // non charged
        // photons
        if (mcParticle.pdgCode() == PDG_t::kGamma) {
          fillResolvedCorrectionHists<McResolvedTruthLevel::True, CorrAssocType::PhotonPCM>(mcCollision, mcParticle);
          continue;
        }
        continue;
      }

      // decaying particles (not marked physical primary)
      if (!checkDecayPrimary(mcParticle)) {
        continue;
      }
      // pi0
      if (mcParticle.pdgCode() == PDG_t::kPi0) {
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_true_pi0"), mcCollision.ptMax(), mcParticle.pt());
        continue;
      }
      // eta
      if (mcParticle.pdgCode() == constants::physics::Pdg::kEta) {
        histos.fill(HIST("mc/correction/h2_ptTrigPtAssoc_true_eta"), mcCollision.ptMax(), mcParticle.pt());
        continue;
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcTrueCorrection, "process mc true data to calculate corrections", false);

  void processPCMDCAzTrue(CorrMcDCollision const& collision, aod::PhotonPCMs const& photonPCMs, aod::V0PhotonsKF const&, aod::JetTracksMCD const&,
                          CorrMcCollisions const&, aod::JetParticles const&)
  {
    // event selection
    if (!totalEvSel(collision)) {
      return;
    }

    for (auto const& photonPCM : photonPCMs) {
      // check mc
      auto const& posTrack = photonPCM.posJetTrack_as<aod::JetTracksMCD>();
      auto const& negTrack = photonPCM.negJetTrack_as<aod::JetTracksMCD>();
      if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle()) {
        continue;
      }
      // collision association
      if (posTrack.mcParticle().mcCollisionId() != collision.mcCollisionId() || negTrack.mcParticle().mcCollisionId() != collision.mcCollisionId()) {
        continue;
      }

      histos.fill(HIST("mc/plain/h5_ptTrigPtAssocDCAzZPvMult_photonPCM"),
                  collision.ptMax(), photonPCM.pt(), photonPCM.v0PhotonKF().dcaZtopv(), collision.posZ(), collision.nGlobalTracks());
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processPCMDCAzTrue, "process to test true DCAz distribution of photons in trigger collisions", false);

  // last /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  void processRecoLast(CorrCollisions const& collisions, aod::Triggers const& triggers)
  {
    // do at end of each data frame (after other reco correlation process functions)
    // (PROCESS_SWITCH of this process has to be declared last)

    for (auto const& collision : collisions) {
      // all events
      histos.fill(HIST("reco/info/h1_nEvents"), 0.5);

      // mc split
      if (!checkSplitMcEventSelection(collision)) {
        continue;
      }
      histos.fill(HIST("reco/info/h1_nEvents"), 1.5);

      // standard event selection
      if (!collision.selEv()) {
        continue;
      }
      histos.fill(HIST("reco/info/h1_nEvents"), 2.5);

      histos.fill(HIST("reco/info/h3_ptTrigZPvMult"), collision.ptMax(), collision.posZ(), collision.nGlobalTracks());
      histos.fill(HIST("reco/info/h2_ptTrigOccupancy"), collision.ptMax(), collision.trackOccupancyInTimeRange());

      // group collision
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());

      // trigger loop
      for (auto const& trigger : triggersThisEvent) {
        // trigger info
        histos.fill(HIST("reco/corr/h3_ptPhiEta_trig"), trigger.pt(), trigger.phi(), trigger.eta(),
                    getPtCorrection<CorrectionParticleType::Trigger>(trigger.pt()));

        // save trigger for mixing
        triggerMemoryReco->saveTrigger(trigger.pt(), trigger.phi(), trigger.eta(), collision.posZ(), collision.nGlobalTracks());
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processRecoLast, "process general info for reconstructed collisions and later correlations last", true);

  void processMcLast(CorrMcDCollisions const& collisions, aod::Triggers const& triggers, aod::JetTracksMCD const&,
                     CorrMcCollisions const& mcCollisions, aod::TriggerParticles const& triggerParticles, aod::JetParticles const&)
  {
    for (auto const& mcCollision : mcCollisions) {
      // event selection
      if (!totalEvSel(mcCollision)) {
        continue;
      }

      // group collision
      auto const triggerParticlesThisEvent = triggerParticles.sliceBy(perColTriggerParticles, mcCollision.globalIndex());

      histos.fill(HIST("mc/info/h3_ptTrigZPvMult_true"), mcCollision.ptMax(), mcCollision.posZ(), mcCollision.nChargedInEtaRange());

      // trigger loop
      for (auto const& trigger : triggerParticlesThisEvent) {
        // trigger info
        histos.fill(HIST("mc/corr/h3_ptPhiEta_trig_true"), trigger.pt(), trigger.phi(), trigger.eta());

        // save trigger for mixing
        triggerMemoryTrue->saveTrigger(trigger.pt(), trigger.phi(), trigger.eta(), mcCollision.posZ(), mcCollision.nChargedInEtaRange());
      }
    }
    for (auto const& collision : collisions) {
      // event selection
      if (!totalEvSel(collision)) {
        continue;
      }

      auto const& mcCollision = collision.mcCollision_as<CorrMcCollisions>();

      // group collision
      auto const triggerParticlesThisEvent = triggerParticles.sliceBy(perColTriggerParticles, mcCollision.globalIndex());
      auto const triggersThisEvent = triggers.sliceBy(perColTriggers, collision.globalIndex());

      histos.fill(HIST("mc/info/h3_ptTrigZPvMult_trueAssocEv"), mcCollision.ptMax(), mcCollision.posZ(), mcCollision.nChargedInEtaRange());
      histos.fill(HIST("mc/info/h3_ptTrigZPvMult_trueAssocEvRecoPtTrig"), collision.ptMax(), mcCollision.posZ(), mcCollision.nChargedInEtaRange());

      // trigger loops
      for (auto const& triggerParticle : triggerParticlesThisEvent) {
        // trigger info
        histos.fill(HIST("mc/corr/h3_ptPhiEta_trig_trueAssocEv"), triggerParticle.pt(), triggerParticle.phi(), triggerParticle.eta());

        // save trigger for mixing
        triggerMemoryTrueAssocEv->saveTrigger(triggerParticle.pt(), triggerParticle.phi(), triggerParticle.eta(), mcCollision.posZ(), mcCollision.nChargedInEtaRange());
      }
      for (auto const& trigger : triggersThisEvent) {
        if (!trigger.jetTrack_as<aod::JetTracksMCD>().has_mcParticle()) {
          continue;
        }

        auto const& triggerParticle = trigger.jetTrack_as<aod::JetTracksMCD>().mcParticle();

        // trigger info
        histos.fill(HIST("mc/corr/h3_ptPhiEta_trig_trueAssocEvRecoPtTrig"), trigger.pt(), triggerParticle.phi(), triggerParticle.eta());

        // save trigger for mixing
        triggerMemoryTrueAssocEvRecoPtTrig->saveTrigger(trigger.pt(), triggerParticle.phi(), triggerParticle.eta(), mcCollision.posZ(), mcCollision.nChargedInEtaRange());
      }
    }
  }
  PROCESS_SWITCH(PhotonChargedTriggerCorrelation, processMcLast, "process general MC info for collisions and later correlations last", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& configContext)
{
  return WorkflowSpec{adaptAnalysisTask<PhotonChargedTriggerCorrelation>(configContext)};
}
