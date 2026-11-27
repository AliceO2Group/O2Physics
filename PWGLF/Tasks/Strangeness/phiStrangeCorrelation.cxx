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
/// \file phiStrangeCorrelation.cxx
/// \brief Analysis task for phi-strangeness rapidity correlations analysis
/// \author Stefano Cannito (stefano.cannito@cern.ch)

#include "PWGLF/DataModel/LFPhiStrangeCorrelationTables.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"

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
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/StaticFor.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TPDGCode.h>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <functional>
#include <iterator>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

enum AnalysisMode {
  kOldNormalization = 0,
  kMassvsMass,
  kDeltaYvsDeltaPhi
};

enum AssociatedParticleType {
  kK0S = 0,
  kLambda,
  kAntiLambda,
  kXi,
  kOmega,
  kPion,
  kAssocPartSize
};

enum PhiMassRegion {
  kSignalRegion = 0,
  kSidebandRegion,
  kPhiMassRegions
};

template <AssociatedParticleType PartType>
constexpr int getPdgCode()
{
  if constexpr (PartType == kK0S) {
    return PDG_t::kK0Short;
  } else if constexpr (PartType == kLambda) {
    return PDG_t::kLambda0;
  } else if constexpr (PartType == kAntiLambda) {
    return PDG_t::kLambda0Bar;
  } else if constexpr (PartType == kXi) {
    return PDG_t::kXiMinus;
  } else if constexpr (PartType == kOmega) {
    return PDG_t::kOmegaMinus;
  } else if constexpr (PartType == kPion) {
    return PDG_t::kPiPlus;
  } else {
    static_assert(PartType == kK0S || PartType == kLambda || PartType == kAntiLambda ||
                    PartType == kXi || PartType == kOmega || PartType == kPion,
                  "Unsupported particle type in getPdgCode");
    return 0;
  }
}

// Non-owning: CCDB owns its objects. With caching on (the default) the manager keeps each map and deletes it when a
// newer one replaces it, so wrapping the pointer in a smart pointer deletes the map while the cache still serves it
using EffMapPtr = std::variant<const TH2*, const TH3*>;

struct BoundEfficiencyMap {
  using CoordsTuple = std::tuple<float, float, float>;

  const EffMapPtr* effMap;
  CoordsTuple coords;

  BoundEfficiencyMap(const EffMapPtr& effMap, float x, float y, float z) : effMap(&effMap), coords(x, y, z) {}
  BoundEfficiencyMap(const EffMapPtr& effMap, const CoordsTuple& coords) : effMap(&effMap), coords(coords) {}

  [[nodiscard]] float getBinEfficiency() const
  {
    return std::visit(
      [this](auto&& mapPtr) -> float {
        if (!mapPtr) {
          return 1.0f;
        }

        const auto& [x, y, z] = coords;

        // Compile-time branching on the alternative held (TH2 or TH3): generates the exact correct function call
        if constexpr (std::is_same_v<std::decay_t<decltype(mapPtr)>, const TH2*>) {
          return mapPtr->GetBinContent(mapPtr->FindFixBin(y, z)); // 2D case only
        } else {
          return mapPtr->GetBinContent(mapPtr->FindFixBin(x, y, z)); // Full 3D case
        }
      },
      *effMap);
  }

  [[nodiscard]] float interpolateEfficiency() const
  {
    return std::visit(
      [this](auto&& mapPtr) -> float {
        if (!mapPtr) {
          return 1.0f;
        }

        const auto& [x, y, z] = coords;

        if constexpr (std::is_same_v<std::decay_t<decltype(mapPtr)>, const TH2*>) {
          return mapPtr->Interpolate(y, z); // Native 2D interpolation
        } else {
          return mapPtr->Interpolate(x, y, z); // Native 3D interpolation
        }
      },
      *effMap);
  }
};

template <AssociatedParticleType PartType, typename TCands>
struct AssocInput {
  static constexpr AssociatedParticleType Type = PartType;
  TCands const* candidates;
};

template <AssociatedParticleType PartType, typename TCands>
constexpr auto makeAssocInput(TCands const& cands)
{
  return AssocInput<PartType, TCands>{&cands};
}

struct PhiStrangeCorrelation {
  HistogramRegistry histos{"phiStrangenessCorrelation", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Configurable for selection type
  Configurable<int> eventSelectionType{"eventSelectionType", 0, "Event selection type: 0 - default selection only, 1 - default + phi meson selection"};

  // Configurable for analysis mode
  Configurable<int> analysisMode{"analysisMode", kDeltaYvsDeltaPhi, "Analysis mode: 0 - old method with online normalization, 1 - new method with offline normlization, 2 - deltay vs deltaphi"};

  // Configurable for event selection
  Configurable<float> cutZVertex{"cutZVertex", 10.0f, "Accepted z-vertex range (cm)"};

  // Configurable on multiplicity bins
  Configurable<std::vector<double>> binsMult{"binsMult", {0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0}, "Multiplicity bin limits"};

  // Configurables on phi mass selection
  struct : ConfigurableGroup {
    Configurable<std::pair<float, float>> rangeMPhiSignal{"rangeMPhiSignal", {1.0f, 1.05f}, "Phi mass range for signal extraction"};
    Configurable<std::pair<float, float>> rangeMPhiSideband{"rangeMPhiSideband", {1.06f, 1.08f}, "Phi mass range for sideband extraction"};
  } phiConfigs;

  // Configurables for K0s selection
  struct : ConfigurableGroup {
    Configurable<bool> selectK0sInSigRegion{"selectK0sInSigRegion", true, "Select K0s candidates in signal region"};
    Configurable<std::pair<float, float>> rangeMK0sSignal{"rangeMK0sSignal", {0.475f, 0.515f}, "K0S mass range for signal extraction"};
  } k0sConfigs;

  // Configurables for Lambda selection
  struct : ConfigurableGroup {
    Configurable<bool> selectLambdaInSigRegion{"selectLambdaInSigRegion", true, "Select Lambda candidates in signal region"};
    Configurable<std::pair<float, float>> rangeMLambdaSignal{"rangeMLambdaSignal", {1.105f, 1.128f}, "Lambda mass range for signal extraction"};
  } lambdaConfigs;

  // Configurables for Xi selection
  struct : ConfigurableGroup {
    Configurable<bool> selectXiInSigRegion{"selectXiInSigRegion", true, "Select Xi candidates in signal region"};
    Configurable<std::pair<float, float>> rangeMXiSignal{"rangeMXiSignal", {1.305f, 1.335f}, "Xi mass range for signal extraction"};
  } xiConfigs;

  // Configurables for Omega selection
  struct : ConfigurableGroup {
    Configurable<bool> selectOmegaInSigRegion{"selectOmegaInSigRegion", true, "Select Omega candidates in signal region"};
    Configurable<std::pair<float, float>> rangeMOmegaSignal{"rangeMOmegaSignal", {1.655f, 1.685f}, "Omega mass range for signal extraction"};
  } omegaConfigs;

  // Configurables for Pions selection
  struct : ConfigurableGroup {
    Configurable<bool> selectPionInSigRegion{"selectPionInSigRegion", true, "Select Pion candidates in signal region"};
    Configurable<float> pidTPCMax{"pidTPCMax", 3.0f, "Maximum nSigma TPC"};
    Configurable<float> pidTOFMax{"pidTOFMax", 3.0f, "Maximum nSigma TOF"};
  } pionConfigs;

  // Configurables on phi pT bins
  Configurable<std::vector<double>> binspTPhi{"binspTPhi", {0.4, 0.8, 1.4, 2.0, 2.8, 4.0, 6.0, 10.0}, "pT bin limits for Phi"};
  Configurable<std::vector<double>> binspTPhiExt{"binspTPhiExt", {0.0, 0.4, 0.8, 1.4, 2.0, 2.8, 4.0, 6.0, 10.0}, "pT bin limits for Phi extended for MC Gen"};

  // Configurable on K0S pT bins
  Configurable<std::vector<double>> binspTK0S{"binspTK0S", {0.0, 0.3, 0.5, 0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 4.2, 6.4}, "pT bin limits for K0S"};
  Configurable<std::vector<double>> binspTK0SExt{"binspTK0SExt", {0.0, 0.1, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for K0S extended for MC Gen"};

  // Configurable on Lambda pT bins
  Configurable<std::vector<double>> binspTLambda{"binspTLambda", {0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 4.2, 6.4}, "pT bin limits for Lambda"};
  Configurable<std::vector<double>> binspTLambdaExt{"binspTLambdaExt", {0.0, 0.1, 0.5, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for Lambda extended for MC Gen"};

  // Configurable on Xi pT bins
  Configurable<std::vector<double>> binspTXi{"binspTXi", {0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 4.2, 6.5}, "pT bin limits for Xi"};
  Configurable<std::vector<double>> binspTXiExt{"binspTXiExt", {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for Xi extended for MC Gen"};

  // Configurable on Omega pT bins
  Configurable<std::vector<double>> binspTOmega{"binspTOmega", {0.8, 1.2, 1.6, 2.0, 2.4, 3.0, 4.2, 6.5}, "pT bin limits for Omega"};
  Configurable<std::vector<double>> binspTOmegaExt{"binspTOmegaExt", {0.0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 6.0}, "pT bin limits for Omega extended for MC Gen"};

  // Configurable on pion pT bins
  Configurable<std::vector<double>> binspTPi{"binspTPi", {0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0}, "pT bin limits for Pions"};
  Configurable<std::vector<double>> binspTPiExt{"binspTPiExt", {0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0}, "pT bin limits for Pions extended for MC Gen"};

  // The species correlated with the phi, by name (those of AssocParticleLabels): a train configuration then says which
  // species were on, which an array of switches, written as bare elements, does not
  Configurable<std::vector<std::string>> activeCorrelationTypes{"activeCorrelationTypes", {"K0S", "Lambda", "AntiLambda", "Xi", "Omega", "Pi"}, "Species correlated with the phi, any of: K0S, Lambda, AntiLambda, Xi, Omega, Pi"};

  // Configurables for delta y selection
  struct : ConfigurableGroup {
    Configurable<int> nBinsY{"nBinsY", 20, "Number of bins in y axis"};
    Configurable<int> nBinsDeltaY{"nBinsDeltaY", 20, "Number of bins in deltay axis"};
    Configurable<float> cfgYAcceptance{"cfgYAcceptance", 0.5f, "Rapidity acceptance"};
    Configurable<std::vector<float>> cfgDeltaYAcceptanceBins{"cfgDeltaYAcceptanceBins", {0.5f}, "Rapidity acceptance bins"};
  } yConfigs;

  // Configurables to apply efficiency online and how to
  struct : ConfigurableGroup {
    Configurable<bool> applyEfficiency{"applyEfficiency", false, "Use efficiency for filling histograms"};
    Configurable<bool> useEffInterpolation{"useEffInterpolation", false, "If true, interpolates efficiency map, else uses bin center"};
    Configurable<bool> applyPhiEfficiency{"applyPhiEfficiency", false, "Apply efficiency for Phi candidates"};
    Configurable<bool> perPeriodEfficiency{"perPeriodEfficiency", false, "Use per period efficiency maps"};
  } efficiencyConfigs;

  // Configurable for event mixing
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};

  // Configurables for CCDB
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository to use"};
  Configurable<std::string> ccdbEfficiencyPath{"ccdbEfficiencyPath", "Users/s/scannito/Efficiencies/h3EffMap", "Correction path to file"};

  // Configurables for minimum pt selection in McGen
  struct : ConfigurableGroup {
    Configurable<float> minPtPhi{"minPtPhi", 0.4f, "Minimum pT for Phi candidates"};
    Configurable<float> minPtK0S{"minPtK0S", 0.0f, "Minimum pT for K0S candidates"};
    Configurable<float> minPtLambda{"minPtLambda", 0.8f, "Minimum pT for Lambda candidates"};
    Configurable<float> minPtXi{"minPtXi", 0.8f, "Minimum pT for Xi candidates"};
    Configurable<float> minPtOmega{"minPtOmega", 0.8f, "Minimum pT for Omega candidates"};
    Configurable<float> minPtPion{"minPtPion", 0.2f, "Minimum pT for Pion candidates"};
    Configurable<bool> bypassPtCut{"bypassPtCut", false, "Bypass the minimum pt cut at MCGen level"};
  } minPtMcGenConfigs;

  // The generated phi of the efficiency: every channel by default, so that the phi efficiency holds the BR to K+K- as
  // every other efficiency holds that of its channel; or phi -> K+K- only, the BR then left to be applied explicitly
  Configurable<bool> phiEfficiencyKKOnly{"phiEfficiencyKKOnly", false, "Generated phi of the efficiency: false - every decay channel (the efficiency includes the BR to K+K-), true - phi -> K+K- only"};

  // Filter on the selection level: eventSelectionType 0 keeps kDefault and above, 1 keeps kDefaultWithPhi
  Filter collisionFilter = aod::lf_selection_event::selectionLevel > as<uint8_t>(eventSelectionType);

  // Defining the type of the collisions for data and MC
  using SelCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults, aod::PhiStrangeEvtSelDataLike>>;
  using SimCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::PVMults, aod::PhiStrangeEvtSelDataLike, aod::McCollisionLabels>>;
  using MCCollisions = soa::Filtered<soa::Join<aod::McCollisions, aod::McCentFT0Ms, aod::PhiStrangeEvtSelMcGen>>;

  // Slice cache and Preslices for table slicing
  SliceCache cache;

  struct : PresliceGroup {
#define DEFINE_PRESLICE(PREF, TAB_PREF, COL_ID)                                              \
  Preslice<aod::TAB_PREF##CandidatesData> PREF##DataPerCollision = aod::COL_ID::collisionId; \
  Preslice<aod::TAB_PREF##CandidatesMcReco> PREF##McRecoPerCollision = aod::COL_ID::collisionId;

    // Trigger
    DEFINE_PRESLICE(phiCand, Phimeson, lf_selection_phi_candidate)

    // Associated Strangeness
    DEFINE_PRESLICE(k0s, K0sReduced, v0)
    DEFINE_PRESLICE(lambda, LambdaReduced, v0)
    DEFINE_PRESLICE(antiLambda, AntiLambdaReduced, v0)
    DEFINE_PRESLICE(xi, XiReduced, cascade)
    DEFINE_PRESLICE(omega, OmegaReduced, cascade)
#undef DEFINE_PRESLICE

    Preslice<aod::PionTracksData> pionTrackDataPerCollision = aod::track::collisionId;
    Preslice<aod::PionTracksMcReco> pionTrackMcRecoPerCollision = aod::track::collisionId;

    Preslice<aod::McParticles> mcPartPerMcCollision = aod::mcparticle::mcCollisionId;
  } preslices;

  // Necessary service to retrieve efficiency maps from CCDB
  Service<ccdb::BasicCCDBManager> ccdb{};

  EffMapPtr effMapPhi;
  std::array<EffMapPtr, kAssocPartSize> effMapsAssoc{};
  std::array<bool, kAssocPartSize> isCorrelationTypeActive{}; // from activeCorrelationTypes, in init

  // Binning policy and axes for mixed event
  ConfigurableAxis axisVertexMixing{"axisVertexMixing", {10, -10.0f, 10.0f}, "Z vertex axis binning for mixing"};
  ConfigurableAxis axisCentralityMixing{"axisCentralityMixing", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 15.0f, 20.0f, 30.0f, 40.0f, 50.0f, 70.0f, 100.0f}, "Multiplicity percentage binning for mixing"};

  using BinningTypeVertexCent = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;
  BinningTypeVertexCent binningOnVertexAndCent{{axisVertexMixing, axisCentralityMixing}};

  static constexpr std::array<std::string_view, kPhiMassRegions> PhiMassRegionLabels{"Signal", "Sideband"};
  static constexpr std::array<std::string_view, kAssocPartSize> AssocParticleLabels{"K0S", "Lambda", "AntiLambda", "Xi", "Omega", "Pi"};

  // Light structures to store only the necessary information for the correlation analysis at MCGen level
  struct MiniParticle {
    float pt{};
    float y{};
    float phi{};
  };

  struct MiniEvent {
    float multiplicity{};
    std::vector<MiniParticle> phiParticles;
    std::array<std::vector<MiniParticle>, kAssocPartSize> assocParticles;
  };

  // Buffer for mixed event, organized as a vector of deques, one for each multiplicity bin, containing the past events with their particles of interest needed for mixing
  std::vector<std::deque<MiniEvent>> eventBuffer;

  int mLastRunNumber{-1};

  void init(InitContext&)
  {
    // From the names to the positions of AssocParticleLabels, which the loops below index
    for (const auto& name : activeCorrelationTypes.value) {
      const auto label = std::ranges::find(AssocParticleLabels, name);
      if (label == AssocParticleLabels.end()) {
        LOG(fatal) << "activeCorrelationTypes: unknown species '" << name << "', expected any of " << fmt::format("{}", fmt::join(AssocParticleLabels, ", "));
      }
      isCorrelationTypeActive[std::distance(AssocParticleLabels.begin(), label)] = true;
    }

    AxisSpec vertexZAxis = {100, -cutZVertex, cutZVertex, "vrtx_{Z} [cm]"};
    AxisSpec yAxis = {yConfigs.nBinsY, -yConfigs.cfgYAcceptance, yConfigs.cfgYAcceptance, "#it{y}"};
    AxisSpec deltayAxis = {yConfigs.nBinsDeltaY, -1.0f, 1.0f, "#Delta#it{y}"};
    AxisSpec deltaphiAxis = {72, -o2::constants::math::PIHalf, o2::constants::math::PIHalf * 3, "#Delta#varphi"};
    AxisSpec binnedmultAxis{std::vector<double>(binsMult), "centFT0M"};
    AxisSpec massPhiAxis = {200, 0.9f, 1.2f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec binnedpTPhiAxis{std::vector<double>(binspTPhi), "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTPhiAxisExt{std::vector<double>(binspTPhiExt), "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec massK0SAxis = {200, 0.45f, 0.55f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec binnedpTK0SAxis{std::vector<double>(binspTK0S), "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTK0SAxisExt{std::vector<double>(binspTK0SExt), "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec massLambdaAxis = {200, 1.1f, 1.13f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec binnedpTLambdaAxis{std::vector<double>(binspTLambda), "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTLambdaAxisExt{std::vector<double>(binspTLambdaExt), "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec massXiAxis = {200, 1.3f, 1.4f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec binnedpTXiAxis{std::vector<double>(binspTXi), "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTXiAxisExt{std::vector<double>(binspTXiExt), "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec massOmegaAxis = {200, 1.65f, 1.7f, "#it{M}_{inv} [GeV/#it{c}^{2}]"};
    AxisSpec binnedpTOmegaAxis{std::vector<double>(binspTOmega), "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTOmegaAxisExt{std::vector<double>(binspTOmegaExt), "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec nSigmaPiAxis = {100, -10.0f, 10.0f, "N#sigma #pi"};
    AxisSpec binnedpTPiAxis{std::vector<double>(binspTPi), "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec binnedpTPiAxisExt{std::vector<double>(binspTPiExt), "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pTBeyondMapAxis = {100, 0.0f, 50.0f, "#it{p}_{T} (GeV/#it{c})"}; // the binned pT axes end where these candidates begin

    std::array<AxisSpec, kAssocPartSize> massAxes = {massK0SAxis, massLambdaAxis, massLambdaAxis, massXiAxis, massOmegaAxis, nSigmaPiAxis};
    std::array<AxisSpec, kAssocPartSize> binnedpTAxes = {binnedpTK0SAxis, binnedpTLambdaAxis, binnedpTLambdaAxis, binnedpTXiAxis, binnedpTOmegaAxis, binnedpTPiAxis};
    std::array<AxisSpec, kAssocPartSize> binnedpTExtAxes = {binnedpTK0SAxisExt, binnedpTLambdaAxisExt, binnedpTLambdaAxisExt, binnedpTXiAxisExt, binnedpTOmegaAxisExt, binnedpTPiAxisExt};

    std::array<std::string, kAssocPartSize> assocTitles = {"K0Short", "Lambda", "AntiLambda", "Xi", "Omega", "Pion"};
    std::array<std::string, kAssocPartSize> dirNames = {"k0s", "lambda", "antilambda", "xi", "omega", "pi"};

    histos.add("phi/h3PhiData", "Invariant mass of Phi in Data", kTH3F, {binnedmultAxis, binnedpTPhiAxis, massPhiAxis}, true);

    // Entries dropped because an efficiency could not be computed (see computeWeight), in the coordinates of the
    // particle whose efficiency failed: SE counts phi candidates (each with all its pairs) and pairs lost to the
    // associated efficiency, ME the pairs lost to the associated efficiency (its phi are those of SE)
    histos.add("zeroEfficiency/SE/h3Phi", "Phi candidates dropped for zero Phi efficiency, SE", kTH3F, {binnedmultAxis, binnedpTPhiAxis, yAxis});

    // Entries dropped because the particle is at or beyond the last pT bin of its efficiency map (isPhiBeyondMap,
    // isAssocBeyondMap), counted as zeroEfficiency/ is
    histos.add("outOfMap/SE/h3Phi", "Phi candidates dropped beyond the Phi efficiency map, SE", kTH3F, {binnedmultAxis, pTBeyondMapAxis, yAxis});

    histos.add("phi/h4PhiMCReco", "Phi in MC Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTPhiAxis, yAxis});
    histos.add("phi/h3PhiMCGen", "Phi in MC Gen", kTH3F, {binnedmultAxis, binnedpTPhiAxis, yAxis});
    histos.add("phi/h4PhiMCGenAssocReco", "Phi in MC Gen Assoc Reco", kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTPhiAxis, yAxis});

    histos.add("phi/h3PhiMCClosureGen", "Phi in MC Gen for MC Closure Test", kTH3F, {binnedmultAxis, binnedpTPhiAxisExt, yAxis});

    for (size_t i = 0; i < kAssocPartSize; ++i) {
      std::string pName(AssocParticleLabels[i]);
      std::string pTitle = assocTitles[i];
      std::string dir = dirNames[i];

      if (i == kPion) {
        histos.add("phiPi/h6PhiPiTPCData", "Invariant mass of Phi vs nSigmaTPC of Pion in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTAxes[i], deltayAxis, massPhiAxis, massAxes[i]}, true);
        histos.add("phiPi/h6PhiPiTPCDataME", "Invariant mass of Phi vs nSigmaTPC of Pion in Data ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTAxes[i], deltayAxis, massPhiAxis, massAxes[i]}, true);
        histos.add("phiPi/h6PhiPiTOFData", "Invariant mass of Phi vs nSigmaTOF of Pion in Data", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTAxes[i], deltayAxis, massPhiAxis, massAxes[i]}, true);
        histos.add("phiPi/h6PhiPiTOFDataME", "Invariant mass of Phi vs nSigmaTOF of Pion in Data ME", kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTAxes[i], deltayAxis, massPhiAxis, massAxes[i]}, true);
      } else {
        histos.add(fmt::format("phi{}/h6Phi{}Data", pName, pName).c_str(), fmt::format("Invariant mass of Phi vs Invariant mass of {} in Data", pTitle).c_str(), kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTAxes[i], deltayAxis, massPhiAxis, massAxes[i]}, true);
        histos.add(fmt::format("phi{}/h6Phi{}DataME", pName, pName).c_str(), fmt::format("Invariant mass of Phi vs Invariant mass of {} in Data ME", pTitle).c_str(), kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTAxes[i], deltayAxis, massPhiAxis, massAxes[i]}, true);
      }

      for (const auto& label : PhiMassRegionLabels) {
        histos.add(fmt::format("phi{}/h5Phi{}Data{}", pName, pName, label).c_str(), fmt::format("Deltay vs deltaphi for Phi and {} in Data", pTitle).c_str(), kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTAxes[i], deltayAxis, deltaphiAxis}, true);
        histos.add(fmt::format("phi{}/h5Phi{}DataME{}", pName, pName, label).c_str(), fmt::format("Deltay vs deltaphi for Phi and {} in Data ME", pTitle).c_str(), kTHnSparseF, {binnedmultAxis, binnedpTPhiAxis, binnedpTAxes[i], deltayAxis, deltaphiAxis}, true);
      }

      histos.add(fmt::format("zeroEfficiency/SE/h3{}", pName).c_str(), fmt::format("Phi-{} pairs dropped for zero {} efficiency, SE", pTitle, pTitle).c_str(), kTH3F, {binnedmultAxis, binnedpTAxes[i], yAxis});
      histos.add(fmt::format("zeroEfficiency/ME/h3{}", pName).c_str(), fmt::format("Phi-{} pairs dropped for zero {} efficiency, ME", pTitle, pTitle).c_str(), kTH3F, {binnedmultAxis, binnedpTAxes[i], yAxis});
      histos.add(fmt::format("outOfMap/SE/h3{}", pName).c_str(), fmt::format("Phi-{} pairs dropped beyond the {} efficiency map, SE", pTitle, pTitle).c_str(), kTH3F, {binnedmultAxis, pTBeyondMapAxis, yAxis});
      histos.add(fmt::format("outOfMap/ME/h3{}", pName).c_str(), fmt::format("Phi-{} pairs dropped beyond the {} efficiency map, ME", pTitle, pTitle).c_str(), kTH3F, {binnedmultAxis, pTBeyondMapAxis, yAxis});

      // h4<X>MCReco is cut in the signal region (mass window, pion nSigma window) only if select<X>InSigRegion is on,
      // as the pairs are (isAssocValid). With selectPionInSigRegion off, the pions' is the efficiency without the
      // window, whose own efficiency is then taken from the producer's nSigma histograms (MC) or from the data fits
      histos.add(fmt::format("{}/h4{}MCReco", dir, pName).c_str(), fmt::format("{} in MC Reco", pTitle).c_str(), kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTAxes[i], yAxis});
      histos.add(fmt::format("{}/h3{}MCGen", dir, pName).c_str(), fmt::format("{} in MC Gen", pTitle).c_str(), kTH3F, {binnedmultAxis, binnedpTAxes[i], yAxis});
      histos.add(fmt::format("{}/h4{}MCGenAssocReco", dir, pName).c_str(), fmt::format("{} in MC Gen Assoc Reco", pTitle).c_str(), kTHnSparseF, {vertexZAxis, binnedmultAxis, binnedpTAxes[i], yAxis});

      histos.add(fmt::format("phi{}/h5Phi{}ClosureMCGen", pName, pName).c_str(), fmt::format("Deltay vs deltaphi for Phi and {} in MCGen", pTitle).c_str(), kTHnSparseF, {binnedmultAxis, binnedpTPhiAxisExt, binnedpTExtAxes[i], deltayAxis, deltaphiAxis});
      histos.add(fmt::format("phi{}/h5Phi{}ClosureMCGenME", pName, pName).c_str(), fmt::format("Deltay vs deltaphi for Phi and {} in MCGen ME", pTitle).c_str(), kTHnSparseF, {binnedmultAxis, binnedpTPhiAxisExt, binnedpTExtAxes[i], deltayAxis, deltaphiAxis});
    }

    histos.add("event/hRecoMCMultiplicityPercent", "RecoMC Multiplicity Percentile", kTH1F, {binnedmultAxis});
    histos.add("event/h2RecoMCVertexZvsMult", "RecoMC Vertex Z vs Multiplicity Percentile", kTH2F, {vertexZAxis, binnedmultAxis});
    histos.add("event/hSplitVertexZ", "Split in z-vtx", kTH1F, {{100, -5.0f, 5.0f}});
    histos.add("event/hGenMCMultiplicityPercent", "Generated MC Multiplicity Percentile", kTH1F, {binnedmultAxis});
    histos.add("event/hGenMCAssocRecoMultiplicityPercent", "Generated MC associated Multiplicity Percentile", kTH1F, {binnedmultAxis});
    histos.add("event/h2GenMCAssocRecoVertexZvsMult", "Generated MC associated reco Vertex Z vs multiplicity", kTH2F, {vertexZAxis, binnedmultAxis});

    // Load efficiency maps from CCDB
    if (efficiencyConfigs.applyEfficiency) {
      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(false);

      if (!efficiencyConfigs.perPeriodEfficiency) {
        loadEfficiencyMaps();
      }
    }

    eventBuffer.resize(binsMult->size() - 1);
  }

  void fetchSingleEfficiencyMapFromCCDB(EffMapPtr& effMap, std::string_view particleName, std::optional<int> runNumber = std::nullopt)
  {
    std::string path = fmt::format("{}{}", ccdbEfficiencyPath.value, particleName);

    // Fetch once as TH1, the common base, and dispatch on the actual type: the CCDB cache is keyed by path only and
    // hands back whatever it holds cast to the requested type, so asking the same path as TH3D and then as TH2D
    // would return a TH2 as a non-null TH3 on a cache hit
    const TH1* map = runNumber ? ccdb->getForRun<TH1>(path, *runNumber) : ccdb->get<TH1>(path);

    if (const auto* map3D = dynamic_cast<const TH3*>(map)) {
      effMap = map3D;
      LOG(debug) << "Efficiency map (TH3) for " << particleName << " loaded from CCDB";
      return;
    }

    if (const auto* map2D = dynamic_cast<const TH2*>(map)) {
      effMap = map2D;
      LOG(debug) << "Efficiency map (TH2) for " << particleName << " loaded from CCDB";
      return;
    }

    std::string runSuffix = runNumber ? fmt::format(" for run {}", *runNumber) : "";
    LOG(fatal) << "Could not load efficiency map (neither TH3 nor TH2) for " << particleName << " from CCDB" << runSuffix << "!";
  }

  void loadEfficiencyMaps(std::optional<int> runNumber = std::nullopt)
  {
    if (runNumber) {
      if (*runNumber == mLastRunNumber) {
        return;
      }
      mLastRunNumber = *runNumber;
    }

    // Load the Trigger (Phi) map if requested by analysis method
    if (efficiencyConfigs.applyPhiEfficiency) {
      fetchSingleEfficiencyMapFromCCDB(effMapPhi, "Phi", runNumber);
    }

    // Only load the associated maps that are explicitly enabled
    for (size_t i = 0; i < kAssocPartSize; ++i) {
      if (isCorrelationTypeActive[i]) {
        fetchSingleEfficiencyMapFromCCDB(effMapsAssoc[i], AssocParticleLabels[i], runNumber);
      }
    }

    if (runNumber) {
      LOG(info) << "Efficiency maps successfully updated from CCDB for run " << *runNumber;
    } else {
      LOG(info) << "Global efficiency maps successfully loaded from CCDB";
    }
  }

  // Weight of an entry: the inverse of the product of the efficiencies of the maps it reads. Only the weight goes into
  // the histograms, whose Sumw2 then gives the statistical variance, sum of w^2. The uncertainty of the maps is not
  // added per entry: it is the same for every entry reading a map bin, so it adds linearly within the bin, which a
  // per-entry sum in quadrature underestimates by sqrt(N); it is propagated offline.
  //
  // Empty when an efficiency cannot be computed: a map bin without reconstructed MC, or a candidate outside the map,
  // reads 0, and interpolation near such bins can give a non-positive or non-finite value. Such an entry is dropped and
  // counted (zeroEfficiency/), as the candidates beyond the last pT bin of the maps are (outOfMap/): a weight of 1
  // instead would mix corrected and uncorrected entries in one bin
  template <typename... BoundEffMaps>
  std::optional<float> computeWeight(const BoundEffMaps&... boundEffMaps)
  {
    if (!efficiencyConfigs.applyEfficiency) {
      return 1.0f;
    }

    float effTot = 1.0f;

    auto processMap = [&](const auto& boundMap) {
      const float eff = efficiencyConfigs.useEffInterpolation ? boundMap.interpolateEfficiency() : boundMap.getBinEfficiency();
      if (!(eff > 0.0f) || !std::isfinite(eff)) {
        return false;
      }
      effTot *= eff;
      return true;
    };

    if (!(processMap(boundEffMaps) && ...)) { // Fold expression, stops at the first efficiency that cannot be computed
      return std::nullopt;
    }

    return 1.0f / effTot;
  }

  template <AssociatedParticleType PartType>
  std::optional<float> computePhiAssocWeight(float multiplicity, const auto& phiCand, const auto& assocCand)
  {
    return computeWeight(BoundEfficiencyMap(effMapPhi, multiplicity, phiCand.pt(), phiCand.y()),
                         BoundEfficiencyMap(effMapsAssoc[PartType], multiplicity, assocCand.pt(), assocCand.y()));
  }

  float getDeltaPhi(float phiTrigger, float phiAssociated)
  {
    return RecoDecay::constrainAngle(phiTrigger - phiAssociated, -o2::constants::math::PIHalf);
  }

  int getCentBin(float multiplicity)
  {
    if (multiplicity < binsMult->front() || multiplicity >= binsMult->back()) {
      return -1;
    }

    auto it = std::upper_bound(binsMult->begin(), binsMult->end(), multiplicity);
    return std::distance(binsMult->begin(), it) - 1;
  }

  template <AssociatedParticleType PartType>
  auto getPreslice()
  {
    if constexpr (PartType == kK0S) {
      return preslices.k0sMcRecoPerCollision;
    } else if constexpr (PartType == kLambda) {
      return preslices.lambdaMcRecoPerCollision;
    } else if constexpr (PartType == kAntiLambda) {
      return preslices.antiLambdaMcRecoPerCollision;
    } else if constexpr (PartType == kXi) {
      return preslices.xiMcRecoPerCollision;
    } else if constexpr (PartType == kOmega) {
      return preslices.omegaMcRecoPerCollision;
    } else if constexpr (PartType == kPion) {
      return preslices.pionTrackMcRecoPerCollision;
    } else {
      static_assert(PartType == kK0S || PartType == kLambda || PartType == kAntiLambda ||
                      PartType == kXi || PartType == kOmega || PartType == kPion,
                    "Unsupported particle type in getPreslice");
      return preslices.k0sMcRecoPerCollision;
    }
  }

  template <AssociatedParticleType PartType>
  float getMinGenPt() const
  {
    if constexpr (PartType == kK0S) {
      return minPtMcGenConfigs.minPtK0S.value;
    } else if constexpr (PartType == kLambda || PartType == kAntiLambda) {
      return minPtMcGenConfigs.minPtLambda.value;
    } else if constexpr (PartType == kXi) {
      return minPtMcGenConfigs.minPtXi.value;
    } else if constexpr (PartType == kOmega) {
      return minPtMcGenConfigs.minPtOmega.value;
    } else if constexpr (PartType == kPion) {
      return minPtMcGenConfigs.minPtPion.value;
    } else {
      static_assert(PartType == kK0S || PartType == kLambda || PartType == kAntiLambda ||
                      PartType == kXi || PartType == kOmega || PartType == kPion,
                    "Unsupported particle type in getMinGenPt");
      return 0.f;
    }
  }

  // A candidate at or beyond the last pT bin of its efficiency map has no efficiency to read: with the efficiencies
  // applied it is dropped, and counted in outOfMap/
  bool isPhiBeyondMap(float pt)
  {
    return efficiencyConfigs.applyEfficiency && efficiencyConfigs.applyPhiEfficiency && pt >= binspTPhi.value.back();
  }

  template <AssociatedParticleType PartType>
  bool isAssocBeyondMap(float pt)
  {
    if (!efficiencyConfigs.applyEfficiency) {
      return false;
    }
    if constexpr (PartType == kK0S) {
      return pt >= binspTK0S.value.back();
    } else if constexpr (PartType == kLambda || PartType == kAntiLambda) {
      return pt >= binspTLambda.value.back();
    } else if constexpr (PartType == kXi) {
      return pt >= binspTXi.value.back();
    } else if constexpr (PartType == kOmega) {
      return pt >= binspTOmega.value.back();
    } else {
      static_assert(PartType == kPion, "Unsupported particle type in isAssocBeyondMap");
      return pt >= binspTPi.value.back();
    }
  }

  // The signal-region selection (mass or nSigma window) of an associated candidate, for the pairs and for the MC
  // reconstructed of the efficiency alike; the pT range of the efficiency maps is checked apart (isAssocBeyondMap), only
  // where an efficiency is applied
  template <AssociatedParticleType PartType, typename TAssoc>
  bool isAssocValid(TAssoc const& assoc)
  {
    if constexpr (PartType == kK0S) {
      const bool applyMassCut = (analysisMode == kDeltaYvsDeltaPhi) && k0sConfigs.selectK0sInSigRegion;
      const auto& [minMass, maxMass] = k0sConfigs.rangeMK0sSignal.value;
      return !applyMassCut || assoc.inMassRegion(minMass, maxMass);
    } else if constexpr (PartType == kLambda || PartType == kAntiLambda) {
      const bool applyMassCut = (analysisMode == kDeltaYvsDeltaPhi) && lambdaConfigs.selectLambdaInSigRegion;
      const auto& [minMass, maxMass] = lambdaConfigs.rangeMLambdaSignal.value;
      return !applyMassCut || assoc.inMassRegion(minMass, maxMass);
    } else if constexpr (PartType == kXi) {
      const bool applyMassCut = (analysisMode == kDeltaYvsDeltaPhi) && xiConfigs.selectXiInSigRegion;
      const auto& [minMass, maxMass] = xiConfigs.rangeMXiSignal.value;
      return !applyMassCut || assoc.inMassRegion(minMass, maxMass);
    } else if constexpr (PartType == kOmega) {
      const bool applyMassCut = (analysisMode == kDeltaYvsDeltaPhi) && omegaConfigs.selectOmegaInSigRegion;
      const auto& [minMass, maxMass] = omegaConfigs.rangeMOmegaSignal.value;
      return !applyMassCut || assoc.inMassRegion(minMass, maxMass);
    } else if constexpr (PartType == kPion) {
      const bool applyNSigmaCut = (analysisMode == kDeltaYvsDeltaPhi) && pionConfigs.selectPionInSigRegion;
      return !applyNSigmaCut || assoc.inNSigmaRegion(pionConfigs.pidTPCMax, pionConfigs.pidTOFMax);
    } else {
      static_assert(PartType == kK0S || PartType == kLambda || PartType == kAntiLambda ||
                      PartType == kXi || PartType == kOmega || PartType == kPion,
                    "Unsupported particle type in isAssocValid");
      return false;
    }
  }

  template <AssociatedParticleType PartType>
  bool isGenSpeciesValid(const auto& mcParticle, bool allowBypass) const
  {
    return mcParticle.isPhysicalPrimary() && (allowBypass || mcParticle.pt() >= getMinGenPt<PartType>());
  }

  template <AssociatedParticleType PartType, bool IsME, typename TPhiCand, typename TAssoc>
  void fillPhiAssocCorrelation(TPhiCand const& phiCand, TAssoc const& assoc, float multiplicity, float weight)
  {
    const std::array<std::pair<float, float>, kPhiMassRegions> phiMassRegions = {phiConfigs.rangeMPhiSignal, phiConfigs.rangeMPhiSideband};

    float deltaY = phiCand.y() - assoc.y();
    float deltaPhi = getDeltaPhi(phiCand.phi(), assoc.phi());

    if (analysisMode == kMassvsMass) {

#define FILL_MASS_CASE(PART_ENUM, PART_STR)                                                                                                         \
  if constexpr (PartType == PART_ENUM) {                                                                                                            \
    if constexpr (IsME) {                                                                                                                           \
      histos.fill(HIST("phi" PART_STR "/h6Phi" PART_STR "DataME"), multiplicity, phiCand.pt(), assoc.pt(), deltaY, phiCand.m(), assoc.m(), weight); \
    } else {                                                                                                                                        \
      histos.fill(HIST("phi" PART_STR "/h6Phi" PART_STR "Data"), multiplicity, phiCand.pt(), assoc.pt(), deltaY, phiCand.m(), assoc.m(), weight);   \
    }                                                                                                                                               \
  }

      FILL_MASS_CASE(kK0S, "K0S");
      FILL_MASS_CASE(kLambda, "Lambda");
      FILL_MASS_CASE(kAntiLambda, "AntiLambda");
      FILL_MASS_CASE(kXi, "Xi");
      FILL_MASS_CASE(kOmega, "Omega");
      if constexpr (PartType == kPion) {
        if constexpr (IsME) {
          histos.fill(HIST("phiPi/h6PhiPiTPCDataME"), multiplicity, phiCand.pt(), assoc.pt(), deltaY, phiCand.m(), assoc.nSigmaTPC(), weight);
          histos.fill(HIST("phiPi/h6PhiPiTOFDataME"), multiplicity, phiCand.pt(), assoc.pt(), deltaY, phiCand.m(), assoc.nSigmaTOF(), weight);
        } else {
          histos.fill(HIST("phiPi/h6PhiPiTPCData"), multiplicity, phiCand.pt(), assoc.pt(), deltaY, phiCand.m(), assoc.nSigmaTPC(), weight);
          histos.fill(HIST("phiPi/h6PhiPiTOFData"), multiplicity, phiCand.pt(), assoc.pt(), deltaY, phiCand.m(), assoc.nSigmaTOF(), weight);
        }
      }
#undef FILL_MASS_CASE

    } else if (analysisMode == kDeltaYvsDeltaPhi) {
      std::function<void()> fillSignal;
      std::function<void()> fillSideband;

#define FILL_DYDP_CASE(PART_ENUM, PART_STR)                                                                                                                             \
  if constexpr (PartType == PART_ENUM) {                                                                                                                                \
    if constexpr (IsME) {                                                                                                                                               \
      fillSignal = [&]() { histos.fill(HIST("phi" PART_STR "/h5Phi" PART_STR "DataMESignal"), multiplicity, phiCand.pt(), assoc.pt(), deltaY, deltaPhi, weight); };     \
      fillSideband = [&]() { histos.fill(HIST("phi" PART_STR "/h5Phi" PART_STR "DataMESideband"), multiplicity, phiCand.pt(), assoc.pt(), deltaY, deltaPhi, weight); }; \
    } else {                                                                                                                                                            \
      fillSignal = [&]() { histos.fill(HIST("phi" PART_STR "/h5Phi" PART_STR "DataSignal"), multiplicity, phiCand.pt(), assoc.pt(), deltaY, deltaPhi, weight); };       \
      fillSideband = [&]() { histos.fill(HIST("phi" PART_STR "/h5Phi" PART_STR "DataSideband"), multiplicity, phiCand.pt(), assoc.pt(), deltaY, deltaPhi, weight); };   \
    }                                                                                                                                                                   \
  }

      FILL_DYDP_CASE(kK0S, "K0S");
      FILL_DYDP_CASE(kLambda, "Lambda");
      FILL_DYDP_CASE(kAntiLambda, "AntiLambda");
      FILL_DYDP_CASE(kXi, "Xi");
      FILL_DYDP_CASE(kOmega, "Omega");
      FILL_DYDP_CASE(kPion, "Pi");
#undef FILL_DYDP_CASE

      static_for<0, kPhiMassRegions - 1>([&](auto i_idx) {
        constexpr unsigned int Idx = i_idx.value;
        const auto& [minMassPhi, maxMassPhi] = phiMassRegions[Idx];
        if (!phiCand.inMassRegion(minMassPhi, maxMassPhi)) {
          return;
        }
        if constexpr (Idx == kSignalRegion)
          fillSignal();
        else
          fillSideband();
      });
    }
  }

  template <AssociatedParticleType PartType, bool IsME, typename TPhiCand, typename TAssoc>
  void processOneAssocPair(TPhiCand const& phiCand, TAssoc const& assoc, float multiplicity)
  {
    if (isAssocBeyondMap<PartType>(assoc.pt())) {
      if constexpr (IsME) {
        histos.fill(HIST("outOfMap/ME/h3") + HIST(AssocParticleLabels[PartType]), multiplicity, assoc.pt(), assoc.y());
      } else {
        histos.fill(HIST("outOfMap/SE/h3") + HIST(AssocParticleLabels[PartType]), multiplicity, assoc.pt(), assoc.y());
      }
      return;
    }

    if (!isAssocValid<PartType>(assoc)) {
      return;
    }

    // The callers have already checked the phi efficiency, so a weight that cannot be computed here comes from the
    // associated one
    const auto weight = computePhiAssocWeight<PartType>(multiplicity, phiCand, assoc);
    if (!weight) {
      if constexpr (IsME) {
        histos.fill(HIST("zeroEfficiency/ME/h3") + HIST(AssocParticleLabels[PartType]), multiplicity, assoc.pt(), assoc.y());
      } else {
        histos.fill(HIST("zeroEfficiency/SE/h3") + HIST(AssocParticleLabels[PartType]), multiplicity, assoc.pt(), assoc.y());
      }
      return;
    }

    fillPhiAssocCorrelation<PartType, IsME>(phiCand, assoc, multiplicity, *weight);
  }

  template <AssociatedParticleType PartType, typename TPhiCand, typename TAssocCands>
  void processOneAssocSpeciesSE(TPhiCand const& phiCand, TAssocCands const& assocCands, float multiplicity)
  {
    if (!isCorrelationTypeActive[PartType]) {
      return;
    }

    for (const auto& assoc : assocCands) {
      processOneAssocPair<PartType, false>(phiCand, assoc, multiplicity);
    }
  }

  template <typename TCollision, typename TPhiCands, typename... TAssocInputs>
  void processPhiAssocSE(TCollision const& collision, TPhiCands const& phiCandidates, TAssocInputs const&... assocInputs)
  {
    float multiplicity = collision.centFT0M();

    if (efficiencyConfigs.applyEfficiency && efficiencyConfigs.perPeriodEfficiency) {
      auto bc = collision.template bc_as<aod::BCs>();
      loadEfficiencyMaps(bc.runNumber());
    }

    for (const auto& phiCand : phiCandidates) {
      if (isPhiBeyondMap(phiCand.pt())) {
        histos.fill(HIST("outOfMap/SE/h3Phi"), multiplicity, phiCand.pt(), phiCand.y());
        continue;
      }

      // Checked once per candidate, before the trigger count and the pairs, so that both see the same phi candidates
      const auto weightPhi = computeWeight(BoundEfficiencyMap(effMapPhi, multiplicity, phiCand.pt(), phiCand.y()));
      if (!weightPhi) {
        histos.fill(HIST("zeroEfficiency/SE/h3Phi"), multiplicity, phiCand.pt(), phiCand.y());
        continue;
      }

      histos.fill(HIST("phi/h3PhiData"), multiplicity, phiCand.pt(), phiCand.m(), *weightPhi);

      // Fold expression: elaborate the correlation for each associated particle type in the parameter pack
      (processOneAssocSpeciesSE<TAssocInputs::Type>(phiCand, *assocInputs.candidates, multiplicity), ...);
    }
  }

  void processPhiAssocSEDataLike(SelCollisions::iterator const& collision,
                                 aod::PhimesonCandidatesData const& phiCandidates,
                                 aod::K0sReducedCandidatesData const& k0sReduced,
                                 aod::LambdaReducedCandidatesData const& lambdaReduced,
                                 aod::AntiLambdaReducedCandidatesData const& antilambdaReduced,
                                 aod::XiReducedCandidatesData const& xiReduced,
                                 aod::OmegaReducedCandidatesData const& omegaReduced,
                                 aod::PionTracksData const& pionTracks,
                                 aod::BCs const&)
  {
    processPhiAssocSE(collision, phiCandidates,
                      makeAssocInput<kK0S>(k0sReduced),
                      makeAssocInput<kLambda>(lambdaReduced),
                      makeAssocInput<kAntiLambda>(antilambdaReduced),
                      makeAssocInput<kXi>(xiReduced),
                      makeAssocInput<kOmega>(omegaReduced),
                      makeAssocInput<kPion>(pionTracks));
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processPhiAssocSEDataLike, "Process function for Phi-Associated 2D Correlations in Data or MC w/o PDG SE", true);

  void processPhiAssocSEMCWithPDG(SimCollisions::iterator const& collision,
                                  aod::PhimesonCandidatesMcReco const& phiCandidates,
                                  aod::K0sReducedCandidatesMcReco const& k0sReduced,
                                  aod::LambdaReducedCandidatesMcReco const& lambdaReduced,
                                  aod::AntiLambdaReducedCandidatesMcReco const& antilambdaReduced,
                                  aod::XiReducedCandidatesMcReco const& xiReduced,
                                  aod::OmegaReducedCandidatesMcReco const& omegaReduced,
                                  aod::PionTracksMcReco const& pionTracks,
                                  aod::BCs const&)
  {
    processPhiAssocSE(collision, phiCandidates,
                      makeAssocInput<kK0S>(k0sReduced),
                      makeAssocInput<kLambda>(lambdaReduced),
                      makeAssocInput<kAntiLambda>(antilambdaReduced),
                      makeAssocInput<kXi>(xiReduced),
                      makeAssocInput<kOmega>(omegaReduced),
                      makeAssocInput<kPion>(pionTracks));
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processPhiAssocSEMCWithPDG, "Process function for Phi-Associated 2D Correlations in MC with PDG SE", true);

  template <AssociatedParticleType PartType, typename TCollisions, typename TPhiCands, typename TAssocCands>
  void processPhiAssocME(TCollisions const& collisions, TPhiCands const& phiCandidates, TAssocCands const& assocCands)
  {
    // As in the same event: a species left out of activeCorrelationTypes has no pairs, and no efficiency map loaded,
    // whose absence would read as an efficiency of 1
    if (!isCorrelationTypeActive[PartType]) {
      return;
    }

    auto tuplePhiAssoc = std::make_tuple(phiCandidates, assocCands);
    Pair<TCollisions, TPhiCands, TAssocCands, BinningTypeVertexCent> pairPhiAssoc{binningOnVertexAndCent, cfgNoMixedEvents, -1, collisions, tuplePhiAssoc, &cache};

    for (const auto& [c1, phiCands, c2, assocRed] : pairPhiAssoc) {
      float multiplicity = c1.centFT0M();

      if (efficiencyConfigs.applyEfficiency && efficiencyConfigs.perPeriodEfficiency) {
        auto bc = c1.template bc_as<aod::BCs>();
        loadEfficiencyMaps(bc.runNumber());
      }

      for (const auto& [phiCand, assoc] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(phiCands, assocRed))) {
        // Not counted: the phi of c1 are those of the same event, already counted there once each
        if (isPhiBeyondMap(phiCand.pt()) || !computeWeight(BoundEfficiencyMap(effMapPhi, multiplicity, phiCand.pt(), phiCand.y()))) {
          continue;
        }

        processOneAssocPair<PartType, true>(phiCand, assoc, multiplicity);
      }
    }
  }

#define DO_PROCESS_SWITCH(CLASS, FUNC_NAME, DESC) \
  PROCESS_SWITCH(CLASS, FUNC_NAME, DESC, true)

#define DEFINE_ME_PROCESSES(PART_NAME, PART_ENUM, TABLE_DATA, TABLE_MC)                                                     \
  void processPhi##PART_NAME##MEDataLike(SelCollisions const& collisions,                                                   \
                                         aod::PhimesonCandidatesData const& phiCandidates,                                  \
                                         TABLE_DATA const& assocReduced,                                                    \
                                         aod::BCs const&)                                                                   \
  {                                                                                                                         \
    processPhiAssocME<PART_ENUM>(collisions, phiCandidates, assocReduced);                                                  \
  }                                                                                                                         \
                                                                                                                            \
  DO_PROCESS_SWITCH(PhiStrangeCorrelation, processPhi##PART_NAME##MEDataLike, "Phi-" #PART_NAME " ME Data  or MC w/o PDG"); \
                                                                                                                            \
  void processPhi##PART_NAME##MEMCWithPDG(SimCollisions const& collisions,                                                  \
                                          aod::PhimesonCandidatesMcReco const& phiCandidates,                               \
                                          TABLE_MC const& assocReduced,                                                     \
                                          aod::BCs const&)                                                                  \
  {                                                                                                                         \
    processPhiAssocME<PART_ENUM>(collisions, phiCandidates, assocReduced);                                                  \
  }                                                                                                                         \
                                                                                                                            \
  DO_PROCESS_SWITCH(PhiStrangeCorrelation, processPhi##PART_NAME##MEMCWithPDG, "Phi-" #PART_NAME " ME MC w/ PDG");

  DEFINE_ME_PROCESSES(K0S, kK0S, aod::K0sReducedCandidatesData, aod::K0sReducedCandidatesMcReco)
  DEFINE_ME_PROCESSES(Lambda, kLambda, aod::LambdaReducedCandidatesData, aod::LambdaReducedCandidatesMcReco)
  DEFINE_ME_PROCESSES(AntiLambda, kAntiLambda, aod::AntiLambdaReducedCandidatesData, aod::AntiLambdaReducedCandidatesMcReco)
  DEFINE_ME_PROCESSES(Xi, kXi, aod::XiReducedCandidatesData, aod::XiReducedCandidatesMcReco)
  DEFINE_ME_PROCESSES(Omega, kOmega, aod::OmegaReducedCandidatesData, aod::OmegaReducedCandidatesMcReco)
  DEFINE_ME_PROCESSES(Pion, kPion, aod::PionTracksData, aod::PionTracksMcReco)
#undef DEFINE_ME_PROCESSES
#undef DO_PROCESS_SWITCH

  void processParticleEfficiency(MCCollisions const& mcCollisions,
                                 SimCollisions const& collisions,
                                 aod::PhimesonCandidatesMcReco const& phiCandidates,
                                 aod::K0sReducedCandidatesMcReco const& k0sReduced,
                                 aod::LambdaReducedCandidatesMcReco const& lambdaReduced,
                                 aod::AntiLambdaReducedCandidatesMcReco const& antilambdaReduced,
                                 aod::XiReducedCandidatesMcReco const& xiReduced,
                                 aod::OmegaReducedCandidatesMcReco const& omegaReduced,
                                 aod::PionTracksMcReco const& pionTracks, aod::McParticles const& mcParticles)
  {
    std::unordered_map<int, std::vector<int>> collsGrouped;
    collsGrouped.reserve(mcCollisions.size());

    for (const auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }
      const auto& mcCollision = collision.mcCollision_as<MCCollisions>();
      collsGrouped[mcCollision.globalIndex()].push_back(collision.globalIndex());
    }

    std::vector<float> zVtxs;
    zVtxs.reserve(3); // Reasonable number of associated collisions to expect at most

    for (const auto& mcCollision : mcCollisions) {
      uint16_t numberAssocColls{0};
      zVtxs.clear();

      // const auto& collIndexesThisMcColl = collsGrouped[mcCollision.globalIndex()];
      auto it = collsGrouped.find(mcCollision.globalIndex());
      if (it != collsGrouped.end()) {
        const auto& collIndexesThisMcColl = it->second;

        for (const auto& collIndex : collIndexesThisMcColl) {
          const auto& collision = collisions.rawIteratorAt(collIndex);

          histos.fill(HIST("event/hRecoMCMultiplicityPercent"), mcCollision.centFT0M());
          histos.fill(HIST("event/h2RecoMCVertexZvsMult"), collision.posZ(), mcCollision.centFT0M());

          zVtxs.push_back(collision.posZ());

          if (eventSelectionType == 0) {
            const auto phiCandidatesThisColl = phiCandidates.sliceBy(preslices.phiCandMcRecoPerCollision, collision.globalIndex());
            for (const auto& phiCand : phiCandidatesThisColl) {
              histos.fill(HIST("phi/h4PhiMCReco"), collision.posZ(), mcCollision.centFT0M(), phiCand.pt(), phiCand.y());
            }
          }

          auto fillRecoAssocSpecies = [&]<typename TAssocInputs>(TAssocInputs const& assocCands, auto histoKey) {
            if (!isCorrelationTypeActive[TAssocInputs::Type]) {
              return;
            }

            const auto assocThisColl = assocCands.candidates->sliceBy(getPreslice<TAssocInputs::Type>(), collision.globalIndex());

            for (const auto& assoc : assocThisColl) {
              if (!isAssocValid<TAssocInputs::Type>(assoc)) {
                continue;
              }

              histos.fill(histoKey, collision.posZ(), mcCollision.centFT0M(), assoc.pt(), assoc.y());
            }
          }; // NOLINT(readability/braces)

          fillRecoAssocSpecies(makeAssocInput<kK0S>(k0sReduced), HIST("k0s/h4K0SMCReco"));
          fillRecoAssocSpecies(makeAssocInput<kLambda>(lambdaReduced), HIST("lambda/h4LambdaMCReco"));
          fillRecoAssocSpecies(makeAssocInput<kAntiLambda>(antilambdaReduced), HIST("antilambda/h4AntiLambdaMCReco"));
          fillRecoAssocSpecies(makeAssocInput<kXi>(xiReduced), HIST("xi/h4XiMCReco"));
          fillRecoAssocSpecies(makeAssocInput<kOmega>(omegaReduced), HIST("omega/h4OmegaMCReco"));
          fillRecoAssocSpecies(makeAssocInput<kPion>(pionTracks), HIST("pi/h4PiMCReco"));

          numberAssocColls++;
        }
      }

      histos.fill(HIST("event/hGenMCMultiplicityPercent"), mcCollision.centFT0M());

      const bool hasAssoc = (numberAssocColls > 0);
      const float zVtxRef = hasAssoc ? zVtxs[0] : 0.0f;

      if (hasAssoc) {
        if (zVtxs.size() > 1) {
          for (size_t i = 1; i < zVtxs.size(); ++i) {
            histos.fill(HIST("event/hSplitVertexZ"), zVtxs[i] - zVtxRef);
          }
        }

        histos.fill(HIST("event/hGenMCAssocRecoMultiplicityPercent"), mcCollision.centFT0M());
        histos.fill(HIST("event/h2GenMCAssocRecoVertexZvsMult"), zVtxRef, mcCollision.centFT0M());
      }

      const auto mcParticlesThisMcColl = mcParticles.sliceBy(preslices.mcPartPerMcCollision, mcCollision.globalIndex());

      for (const auto& mcParticle : mcParticlesThisMcColl) {
        if (std::abs(mcParticle.y()) > yConfigs.cfgYAcceptance) {
          continue;
        }

        auto fillGenHistos = [&](auto h3Key, auto h4Key) {
          histos.fill(h3Key, mcCollision.centFT0M(), mcParticle.pt(), mcParticle.y());
          if (hasAssoc) {
            histos.fill(h4Key, zVtxRef, mcCollision.centFT0M(), mcParticle.pt(), mcParticle.y());
          }
        };

        auto fillGenAssocSpecies = [&]<AssociatedParticleType PartType>(auto h3Key, auto h4Key) {
          if (!isCorrelationTypeActive[PartType] || !isGenSpeciesValid<PartType>(mcParticle, false)) {
            return;
          }

          fillGenHistos(h3Key, h4Key);
        }; // NOLINT(readability/braces)

        switch (std::abs(mcParticle.pdgCode())) {
          case o2::constants::physics::Pdg::kPhi:
            if (eventSelectionType == 0 && mcParticle.pt() >= minPtMcGenConfigs.minPtPhi &&
                (!phiEfficiencyKKOnly || RecoDecay::isMatchedMCGen(mcParticles, mcParticle, o2::constants::physics::Pdg::kPhi, std::array{-PDG_t::kKPlus, +PDG_t::kKPlus}))) {
              fillGenHistos(HIST("phi/h3PhiMCGen"), HIST("phi/h4PhiMCGenAssocReco"));
            }
            break;
          case getPdgCode<kK0S>():
            fillGenAssocSpecies.template operator()<kK0S>(HIST("k0s/h3K0SMCGen"), HIST("k0s/h4K0SMCGenAssocReco"));
            break;
          case getPdgCode<kLambda>(): // the switch is on |PDG|: the sign tells Lambda from AntiLambda
            if (mcParticle.pdgCode() == getPdgCode<kLambda>()) {
              fillGenAssocSpecies.template operator()<kLambda>(HIST("lambda/h3LambdaMCGen"), HIST("lambda/h4LambdaMCGenAssocReco"));
            } else {
              fillGenAssocSpecies.template operator()<kAntiLambda>(HIST("antilambda/h3AntiLambdaMCGen"), HIST("antilambda/h4AntiLambdaMCGenAssocReco"));
            }
            break;
          case getPdgCode<kXi>():
            fillGenAssocSpecies.template operator()<kXi>(HIST("xi/h3XiMCGen"), HIST("xi/h4XiMCGenAssocReco"));
            break;
          case getPdgCode<kOmega>():
            fillGenAssocSpecies.template operator()<kOmega>(HIST("omega/h3OmegaMCGen"), HIST("omega/h4OmegaMCGenAssocReco"));
            break;
          case getPdgCode<kPion>():
            fillGenAssocSpecies.template operator()<kPion>(HIST("pi/h3PiMCGen"), HIST("pi/h4PiMCGenAssocReco"));
            break;
          default:
            break;
        }
      }
    }
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processParticleEfficiency, "Process function for Efficiency Computation for Particles of Interest", false);

  void processMCGenClosure(MCCollisions::iterator const& mcCollision, aod::McParticles const& mcParticles)
  {
    float multiplicity = mcCollision.centFT0M();

    std::vector<MiniParticle> phiParticles;
    std::array<std::vector<MiniParticle>, kAssocPartSize> assocParticles;

    // Preliminary loop to fill vectors of particles of interest for the current event, applying pt and y cuts
    for (const auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.y()) > yConfigs.cfgYAcceptance) {
        continue;
      }

      auto fillPartCollection = [&]<AssociatedParticleType PartType>() {
        if (isCorrelationTypeActive[PartType] && isGenSpeciesValid<PartType>(mcParticle, minPtMcGenConfigs.bypassPtCut.value)) {
          assocParticles[PartType].emplace_back(mcParticle.pt(), mcParticle.y(), mcParticle.phi());
        }
      }; // NOLINT(readability/braces)

      switch (std::abs(mcParticle.pdgCode())) {
        case o2::constants::physics::Pdg::kPhi:
          // Only phi -> K+K-, the trigger of data: a phi in another channel would be correlated with its own daughters
          // (the K0S of K0S K0L, the pions of pi+pi-pi0), pairs data do not have
          if ((minPtMcGenConfigs.bypassPtCut || mcParticle.pt() >= minPtMcGenConfigs.minPtPhi) &&
              RecoDecay::isMatchedMCGen(mcParticles, mcParticle, o2::constants::physics::Pdg::kPhi, std::array{-PDG_t::kKPlus, +PDG_t::kKPlus})) {
            phiParticles.emplace_back(mcParticle.pt(), mcParticle.y(), mcParticle.phi());
          }
          break;
        case getPdgCode<kK0S>():
          fillPartCollection.template operator()<kK0S>();
          break;
        case getPdgCode<kLambda>(): // the switch is on |PDG|: the sign tells Lambda from AntiLambda
          if (mcParticle.pdgCode() == getPdgCode<kLambda>()) {
            fillPartCollection.template operator()<kLambda>();
          } else {
            fillPartCollection.template operator()<kAntiLambda>();
          }
          break;
        case getPdgCode<kXi>():
          fillPartCollection.template operator()<kXi>();
          break;
        case getPdgCode<kOmega>():
          fillPartCollection.template operator()<kOmega>();
          break;
        case getPdgCode<kPion>():
          fillPartCollection.template operator()<kPion>();
          break;
        default:
          break;
      }
    }

    bool hasAnyAssoc = false;
    for (size_t i = 0; i < kAssocPartSize; ++i) {
      if (!assocParticles[i].empty()) {
        hasAnyAssoc = true;
      }
    }

    if (phiParticles.empty() && !hasAnyAssoc) {
      return;
    }

    int multBin = getCentBin(multiplicity);
    if (multBin < 0) {
      return;
    }

    // Same Event Correlations
    for (const auto& phiParticle : phiParticles) {
      histos.fill(HIST("phi/h3PhiMCClosureGen"), multiplicity, phiParticle.pt, phiParticle.y);

      static_for<0, AssocParticleLabels.size() - 1>([&](auto i_idx) {
        constexpr unsigned int Idx = i_idx.value;
        if (!isCorrelationTypeActive[Idx]) {
          return;
        }

        for (const auto& assocParticle : assocParticles[Idx]) {
          histos.fill(HIST("phi") + HIST(AssocParticleLabels[Idx]) + HIST("/h5Phi") + HIST(AssocParticleLabels[Idx]) + HIST("ClosureMCGen"),
                      multiplicity, phiParticle.pt, assocParticle.pt,
                      phiParticle.y - assocParticle.y,
                      getDeltaPhi(phiParticle.phi, assocParticle.phi));
        }
      });
    }

    // Mixed Event Correlations
    // Loop over past events in the same multiplicity bin and fill histograms with all combinations of current phi particles and past associated particles
    for (const auto& pastEvent : eventBuffer[multBin]) {
      for (const auto& phiParticle : phiParticles) {
        static_for<0, AssocParticleLabels.size() - 1>([&](auto i_idx) {
          constexpr unsigned int Idx = i_idx.value;
          if (!isCorrelationTypeActive[Idx]) {
            return;
          }

          for (const auto& assocParticle : pastEvent.assocParticles[Idx]) {
            histos.fill(HIST("phi") + HIST(AssocParticleLabels[Idx]) + HIST("/h5Phi") + HIST(AssocParticleLabels[Idx]) + HIST("ClosureMCGenME"),
                        multiplicity, phiParticle.pt, assocParticle.pt,
                        phiParticle.y - assocParticle.y,
                        getDeltaPhi(phiParticle.phi, assocParticle.phi));
          }
        });
      }
    }

    MiniEvent currentEvent;
    currentEvent.multiplicity = multiplicity;
    currentEvent.phiParticles = std::move(phiParticles);
    currentEvent.assocParticles = std::move(assocParticles);

    eventBuffer[multBin].push_front(std::move(currentEvent));
    if (eventBuffer[multBin].size() > static_cast<std::size_t>(cfgNoMixedEvents.value)) {
      eventBuffer[multBin].pop_back();
    }
  }

  PROCESS_SWITCH(PhiStrangeCorrelation, processMCGenClosure, "Process function for MC Gen Closure Test in SE and ME", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<PhiStrangeCorrelation>(cfgc)};
}
